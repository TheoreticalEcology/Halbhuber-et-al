##############################   Bi-GRU NO TRAITS   ######################################
library(readr)
library(dplyr)
library(torch)
source("~/paper/EcoToxTM/Scripts/functions.R") # for R squared calculation
source("~/paper/EcoToxTM/Scripts/BiGRU_models.R") # for model
df <- readRDS("~/paper/EcoToxTM/Data/df_final_tokenizedSMILES.rds")
df <- na.omit(df)
df <- df[,-1] # exclude SMILES codes 
df[, 22:198] <- lapply(df[, 22:198], as.numeric)

# flip the padded sequence so that the important information is at the back (then backpropagation is better)
remaining_columns_traits <- df[c(1:21)]
reversed_columns <- df[, 22:198][, ncol(df[, 22:198]):1]
df <- cbind(reversed_columns, remaining_columns_traits)


#### Make Tensor ####
# Chemical and Species Traits
data_matrix <- apply(as.matrix(df[,1:177]), 2, as.numeric) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())

# LC50
data_matrix_LC50 <- as.numeric(df[,195])
LC50_tensor <- torch_tensor(data_matrix_LC50, dtype = torch_float())

##### Training and Test datasets ####
set.seed(123)
# Indices
train_indices <- sample(1:nrow(tokens_tensor), 0.8 * nrow(tokens_tensor))
test_indices <- setdiff(1:nrow(tokens_tensor), train_indices)

# Tokens
train_tokens_tensor <- tokens_tensor[train_indices, , drop = FALSE]
test_tokens_tensor <- tokens_tensor[test_indices, , drop = FALSE]
# LC50
train_LC50_tensor <- LC50_tensor[train_indices]
test_LC50_tensor <- LC50_tensor[test_indices]


####  setup  ####
# Device
device <- torch_device("cuda:3")

set.seed(123)

# Hyperparameter grid
embedding <- c(10L, 20L)
lr <- c(0.001, 0.005)
weight_decay <- c(0.0, 0.0001)
grid <- expand.grid(embedding = embedding,
                    lr = lr,
                    weight_decay = weight_decay)

# Fold columns
fold_columns <- list(
  randomFold = df$randomFold,
  speciesFold = df$speciesFold,
  chemicalFold = df$chemicalFold
)

# Dataloader
create_dataloader <- function(data_tensor, batch_size) {
  indices <- sample(1:data_tensor$size(1))
  split(indices, ceiling(seq_along(indices) / batch_size))
}



####  Training & Evaluation ####

train_and_eval_model <- function(train_tokens_tensor, train_LC50_tensor, test_tokens_tensor, test_LC50_tensor,
                                 embedding, lr, weight_decay, num_epochs = 100, batch_size = 128) {

  model <- BiGRU_NoTRAITS(embedding = embedding)
  model$to(device)

  optimizer <- optim_adam(model$parameters, lr = lr, weight_decay = weight_decay)
  loss_fn <- nn_mse_loss()

  train_losses <- c()
  test_r2_scores <- c()

  for (epoch in 1:num_epochs) {
    model$train()
    dataloader_train <- create_dataloader(train_tokens, batch_size)

    for (batch_indices in dataloader_train) {
      x <- train_tokens[batch_indices, , drop = FALSE]
      y <- train_LC50[batch_indices]

      out <- model(x)
      loss <- loss_fn(out, y$unsqueeze(2L))

      loss$backward()
      optimizer$step()
      optimizer$zero_grad()
    }

    # Evaluation
    model$eval()
    dataloader_test <- create_dataloader(test_tokens, batch_size)

    with_no_grad({
      preds <- c()
      trues <- c()
      for (batch_indices in dataloader_test) {
        x <- test_tokens[batch_indices, , drop = FALSE]
        y <- test_LC50[batch_indices]
        out <- model(x)
        preds <- c(preds, as.numeric(out$squeeze()))
        trues <- c(trues, as.numeric(y))
      }
      test_r2_scores <- c(test_r2_scores, calculate_r_squared(trues, preds))
    })
  }

  return(list(model = model, r2 = mean(test_r2_scores)))
}

#### MAIN LOOP ####
results_grid <- list()

for (fold_type in names(fold_columns)) {
  fold_col <- fold_columns[[fold_type]]
  results_grid[[fold_type]] <- list()

  for (g in 1:nrow(grid)) {
    params <- grid[g, ]
    r2_scores <- c()

    # Phase 1: Cross-validation on folds 1-8
    for (fold_id in 1:8) {
      train_idx <- which(fold_col != fold_id & fold_col < 9)
      test_idx  <- which(fold_col == fold_id)

      train_tokens <- tokens_tensor[train_idx, , drop = FALSE]$to(device)
      train_LC50   <- log(LC50_tensor[train_idx])$to(device)
      test_tokens  <- tokens_tensor[test_idx, , drop = FALSE]$to(device)
      test_LC50    <- log(LC50_tensor[test_idx])$to(device)

      res <- train_and_eval_model(
        train_tokens, train_LC50,
        test_tokens, test_LC50,
        embedding = params$embedding,
        lr = params$lr,
        weight_decay = params$weight_decay,
        num_epochs = 100,
        batch_size = 128
      )

      r2_scores <- c(r2_scores, res$r2)
    }

    # Save mean R² for this parameter combo
    mean_r2 <- mean(r2_scores, na.rm = TRUE)
    results_grid[[fold_type]][[g]] <- list(
      hyperparams = params,
      mean_cv_r2 = mean_r2
    )
  }

  #### VALIDATION ####
  best_index <- which.max(sapply(results_grid[[fold_type]], function(x) x$mean_cv_r2))
  best_params <- results_grid[[fold_type]][[best_index]]$hyperparams

  # 2 Phase: pedict on folds 9 und 10
  train_idx_val <- which(fold_col < 9)
  val_idx <- which(fold_col >= 9)

  train_tokens_val <- tokens_tensor[train_idx_val, , drop = FALSE]$to(device)
  train_LC50_val   <- log(LC50_tensor[train_idx_val])$to(device)
  val_tokens       <- tokens_tensor[val_idx, , drop = FALSE]$to(device)
  val_LC50         <- log(LC50_tensor[val_idx])$to(device)

  val_result <- train_and_eval_model(
    train_tokens_val, train_LC50_val,
    val_tokens, val_LC50,
    embedding = best_params$embedding,
    lr = best_params$lr,
    weight_decay = best_params$weight_decay,
    num_epochs = 100,
    batch_size = 128
  )

  results_grid[[fold_type]][[best_index]]$val_r2 <- val_result$r2
  results_grid[[fold_type]][[best_index]]$final_model <- torch::torch_serialize(val_result$model)
}

# library(torch)
# 
# # Vorbereitung
# results_grid <- list()
# fold_columns <- list(
#   randomFold = df$randomFold,
#   speciesFold = df$speciesFold,
#   chemicalFold = df$chemicalFold
# )
# 
# num_epochs <- 1000L
# batch_size <- 400L
# device_used <- torch_device("cuda:3")
# 
# # Starte Schleifen
# for (fold_type in names(fold_columns)) {
#   cat("Running fold type:", fold_type, "\n")
#   fold_column <- fold_columns[[fold_type]]
#   results_grid[[fold_type]] <- list()
#   
#   for (fold_id in 1:10) {
#     cat("  Fold ID:", fold_id, "\n")
#     
#     # Index-Splits
#     train_idx <- which(!fold_column %in% c(fold_id))        # Trainingsdaten für Crossval
#     val_idx   <- which(fold_column %in% c(fold_id))          # Validierungsdaten
#     
#     train_tokens_fold <- tokens_tensor[train_idx, , drop = FALSE]$to(device_used)
#     train_LC50_fold <- LC50_tensor[train_idx]$to(device_used)$log()
#     val_tokens_fold <- tokens_tensor[val_idx, , drop = FALSE]$to(device_used)
#     val_LC50_fold <- LC50_tensor[val_idx]$to(device_used)$log()
#     
#     # Initialisierung
#     model <- BiGRU_NoTRAITS(embedding = 10L)$to(device_used)
#     optimizer <- optim_adam(model$parameters, lr = 0.01, weight_decay = 0.00001)
#     loss_fn <- nn_mse_loss()
#     
#     train_r2_scores <- c()
#     val_r2_scores <- c()
#     best_val_r2 <- -Inf
#     best_model_serialized <- NULL
#     
#     for (epoch in 1:num_epochs) {
#       model$train()
#       
#       # Batching
#       indices <- sample(1:train_tokens_fold$size(1))
#       batches <- split(indices, ceiling(seq_along(indices) / batch_size))
#       
#       for (batch in batches) {
#         x_batch <- train_tokens_fold[batch, , drop = FALSE]
#         y_batch <- train_LC50_fold[batch]
#         
#         pred <- model(x_batch)
#         loss <- loss_fn(pred, y_batch$unsqueeze(2L))
#         
#         loss$backward()
#         optimizer$step()
#         optimizer$zero_grad()
#       }
#       
#       # Evaluation
#       model$eval()
#       with_no_grad({
#         pred_val <- model(val_tokens_fold)
#         r2_val <- calculate_r_squared(as.matrix(val_LC50_fold), as.matrix(pred_val$squeeze()))
#         val_r2_scores <- c(val_r2_scores, r2_val)
#         
#         # Speichere bestes Modell
#         if (!is.na(r2_val) && r2_val > best_val_r2) {
#           best_val_r2 <- r2_val
#           best_model_serialized <- torch_serialize(model)
#         }
#       })
#       
#       cat(sprintf("    Epoch [%d/%d] - Fold: %s_%d | Val R²: %.4f\n",
#                   epoch, num_epochs, fold_type, fold_id, r2_val))
#     }
#     
#     results_grid[[fold_type]][[paste0("fold_", fold_id)]] <- list(
#       val_r2_scores = val_r2_scores,
#       final_model = best_model_serialized,
#       mean_cv_r2 = mean(val_r2_scores, na.rm = TRUE)
#     )
#   }
# }
# 


