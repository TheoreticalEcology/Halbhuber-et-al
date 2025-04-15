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

index_tensor <- function(tensor, indices) {
  # R-Indizes sind 1-basiert, also direkt verwenden
  torch_indices <- torch_tensor(indices, dtype = torch_long())
  tensor$index_select(1, torch_indices) # Umwandlung auf 0-basierte Indizes für index_select()
}






####  Training & Evaluation ####

create_dataloader <- function(data_tensor, batch_size) {
  num_samples <- data_tensor$size(1) 
  indices <- sample(1:num_samples) 
  batches <- split(indices, ceiling(seq_along(indices) / batch_size))
  return(batches) 
}

# ------------------------- #
# Hyperparameter Grid Setup
# ------------------------- #
hyper_grid <- function(lr, dropout) {
  expand.grid(
    lr = lr,
    dropout = dropout,
    train = NA,
    val = NA,
    test = NA,
    stringsAsFactors = FALSE
  )
}

folds <- c("randomFold", "chemicalFold", "speciesFold")

hyper_grid_object <- hyper_grid(lr = c(0.001, 0.005, 0.007, 0.01), dropout = c(0.1, 0.5))
hyperparameters <- list(BiGRU = list())

for (fold in folds) {
  hyperparameters$BiGRU[[fold]] <- hyper_grid_object
}

# ----------------------------------- #
# Train/Evaluate Function mit Dropout
# ----------------------------------- #
train_and_evaluate <- function(train_indices, test_indices, fold_name, lr, dropout, num_epochs = 2L, batch_size = 128) {
  loss_fn <- nn_mse_loss()
  model <- BiGRU_NoTRAITS(embedding = 10L, dropout_prob = dropout)
  model$to(device = "cuda:3")
  
  optimizer <- optim_adam(model$parameters, lr = lr, weight_decay = 0.00001)
  
  train_tokens_fold <- index_tensor(tokens_tensor, train_indices)$to(device = "cuda:3")
  train_LC50_fold <- index_tensor(LC50_tensor, train_indices)$to(device = "cuda:3")
  test_tokens_fold <- index_tensor(tokens_tensor, test_indices)$to(device = "cuda:3")
  test_LC50_fold <- index_tensor(LC50_tensor, test_indices)$to(device = "cuda:3")
  
  
  train_losses <- c()
  test_losses <- c()
  train_r2_scores <- c()
  test_r2_scores <- c()
  
  for (epoch in 1:num_epochs) {
    model$train()
    dataloader_train <- create_dataloader(train_tokens_fold, batch_size = batch_size)
    
    for (batch_indices in dataloader_train) {
      optimizer$zero_grad()
      batch_tokens_train <- train_tokens_fold[batch_indices,,drop=FALSE]
      batch_LC50_train <- train_LC50_fold[batch_indices]
      train_output <- model(batch_tokens_train)
      train_loss <- loss_fn(train_output, batch_LC50_train$unsqueeze(2L))
      train_loss$backward()
      optimizer$step()
      
      train_losses <- c(train_losses, train_loss$item())
      train_r2_scores <- c(train_r2_scores, calculate_r_squared(as.matrix(batch_LC50_train), as.matrix(train_output$squeeze())))
    }
    
    model$eval()
    dataloader_test <- create_dataloader(test_tokens_fold, batch_size = batch_size)
    with_no_grad({
      test_loss_total <- 0
      test_r2_total <- 0
      for (batch_indices in dataloader_test) {
        batch_tokens_test <- test_tokens_fold[batch_indices,,drop=FALSE]
        batch_LC50_test <- test_LC50_fold[batch_indices]
        test_output <- model(batch_tokens_test)
        test_loss_total <- test_loss_total + loss_fn(test_output, batch_LC50_test$unsqueeze(2L))$item()
        test_r2_total <- test_r2_total + calculate_r_squared(as.matrix(batch_LC50_test), as.matrix(test_output$squeeze()))
      }
    })
    test_losses <- c(test_losses, test_loss_total / length(dataloader_test))
    test_r2_scores <- c(test_r2_scores, test_r2_total / length(dataloader_test))
    
    cat(sprintf("Fold: %s | Epoch %d/%d | Train R²: %.4f | Test R²: %.4f\n", fold_name, epoch, num_epochs, tail(train_r2_scores, 1), tail(test_r2_scores, 1)))
  }
  
  list(
    train_losses = train_losses,
    test_losses = test_losses,
    train_r2_scores = train_r2_scores,
    test_r2_scores = test_r2_scores,
    model = torch::torch_serialize(model)
  )
}

# --------------------------------------------- #
# Tuning-Schleife: Folds 3–10 (Train+Val)       #
# --------------------------------------------- #
run_hyperparameter_pipeline <- function(fold_name) {
  grid <- hyperparameters$BiGRU[[fold_name]]
  best_val_r2 <- -Inf
  best_model <- NULL
  best_params <- NULL
  
  for (i in 1:nrow(grid)) {
    lr <- grid$lr[i]
    dropout <- grid$dropout[i]
    
    r2_train_list <- c()
    r2_test_list <- c()
    
    for (fold_id in 3:10) {
      fold_column <- df[[fold_name]]
      train_indices <- which(fold_column != fold_id)
      test_indices <- which(fold_column == fold_id)
      
      result <- train_and_evaluate(train_indices, test_indices, paste0(fold_name, "_", fold_id), lr = lr, dropout = dropout)
      r2_train_list <- c(r2_train_list, tail(result$train_r2_scores, 1))
      r2_test_list <- c(r2_test_list, tail(result$test_r2_scores, 1))
    }
    
    grid$train[i] <- mean(r2_train_list)
    grid$val[i] <- mean(r2_test_list)
    
    if (grid$val[i] > best_val_r2) {
      best_val_r2 <- grid$val[i]
      best_model <- result  # letzte aus Fold 10
      best_params <- list(lr = lr, dropout = dropout)
    }
  }
  
  # ------------------------------------ #
  # Finaltraining auf Fold 3–10         #
  # Evaluation auf Fold 1–2 (Holdout)   #
  # ------------------------------------ #
  train_indices <- which(!(df[[fold_name]] %in% c(1,2)))
  test_indices <- which(df[[fold_name]] %in% c(1,2))
  
  final_model <- train_and_evaluate(train_indices, test_indices = NULL, fold_name = paste0(fold_name, "_finalTrain"), lr = best_params$lr, dropout = best_params$dropout)
  
  # Holdout-Eval
  model_loaded <- torch::torch_unserialize(final_model$model)
  model_loaded$eval()
  
  test_tokens_holdout <- test_tokens_tensor[test_indices,,drop=FALSE]$to(device = "cuda:3")
  test_LC50_holdout <- (test_LC50_tensor[test_indices]$to(device = "cuda:3"))$log()
  test_output <- model_loaded(test_tokens_holdout)
  
  test_loss <- nn_mse_loss()(test_output, test_LC50_holdout$unsqueeze(2L))$item()
  test_r2 <- calculate_r_squared(as.matrix(test_LC50_holdout), as.matrix(test_output$squeeze()))
  
  # Update Grid
  idx <- which(grid$lr == best_params$lr & grid$dropout == best_params$dropout)
  grid$test[idx] <- test_r2
  
  # Save model
  model_path <- paste0("model_final_", fold_name, ".pt")
  torch::torch_save(model_loaded, model_path)
  
  hyperparameters$MLP[[trait]] <<- grid
  
  list(
    best_hyperparams = best_params,
    final_model = model_path,
    holdout_eval = list(test_loss = test_loss, test_r2 = test_r2),
    updated_hyper_grid = hyperparameters
  )
}

# ------------------------------- #
# Run Everything for a Fold Name
# ------------------------------- #
pipeline_result <- run_hyperparameter_pipeline(fold_name = "randomFold")

saveRDS(pipeline_result, file = "~/paper/EcoToxTM/Results/BiGRU_noTraits.rds")


pipeline_result$best_hyperparams      # Beste lr + dropout
pipeline_result$final_model           # Pfad zur .pt-Datei
pipeline_result$holdout_eval          # Test-Loss und Test-R²
pipeline_result$updated_hyper_grid    # Deine Hyperparametertabelle mit train/val/test-R²

### for all Fold names 
results_all <- list()
for (f in folds) {
  results_all[[f]] <- run_hyperparameter_pipeline(fold_name = f)
}