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


####  Training & Evaluation ####

create_dataloader <- function(data_tensor, batch_size) {
  num_samples <- data_tensor$size(1) 
  indices <- sample(1:num_samples) 
  batches <- split(indices, ceiling(seq_along(indices) / batch_size))
  return(batches) 
}

#batches <- create_dataloader(train_tokens_tensor,128)

train_and_evaluate <- function(train_indices, test_indices, fold_name, num_epochs = 1000L, batch_size = 128) {
  loss_fn <- nn_mse_loss()
  model <- BiGRU_NoTRAITS(embedding = 10L)
  model$to(device = "cuda:3") 
  
  train_tokens_fold = train_tokens_tensor$to(device = "cuda:3")
  train_LC50_fold = (train_LC50_tensor$to(device = "cuda:3"))$log()
  test_tokens_fold = test_tokens_tensor$to(device = "cuda:3")
  test_LC50_fold = (test_LC50_tensor$to(device = "cuda:3"))$log()
  
  optimizer <- optim_adam(model$parameters, lr = 0.01, weight_decay = 0.00001)
  
  train_losses <- c()
  test_losses <- c()
  train_r2_scores <- c()
  test_r2_scores <- c()
  
  for (epoch in 1:num_epochs) {
    dataloader_train <- create_dataloader(train_tokens_fold, batch_size = batch_size)
    
    
    model$train()
    
    for (batch_indices in dataloader_train) {
      batch_tokens_train <- train_tokens_fold[batch_indices, , drop = FALSE]
      batch_LC50_train <- train_LC50_fold[batch_indices]
      
      batch_tokens_train <- batch_tokens_train$to(device = "cuda:3")
      batch_LC50_train <- batch_LC50_train$to(device = "cuda:3")
      
      # forward
      train_output <- model(batch_tokens_train)
      train_loss <- loss_fn(train_output, batch_LC50_train$unsqueeze(2L))
      train_losses <- c(train_losses, train_loss$item())
      
      train_r2 <- calculate_r_squared(as.matrix(batch_LC50_train), as.matrix(train_output$squeeze()))
      train_r2_scores <- c(train_r2_scores, train_r2)
      
      # backward
      train_loss$backward()
      optimizer$step()
      optimizer$zero_grad()
    }
    
    model$eval()
    dataloader_test <- create_dataloader(test_tokens_fold, batch_size = batch_size)

    with_no_grad({
      test_loss_total <- 0
      test_r2_total <- 0
      num_batches <- length(dataloader_test)
      
      for (batch_indices in dataloader_test) {
        batch_tokens_test <- test_tokens_fold[batch_indices, , drop = FALSE]
        batch_LC50_test <- test_LC50_fold[batch_indices]
        
        batch_tokens_test <- batch_tokens_test$to(device = "cuda:3")
        batch_LC50_test <- batch_LC50_test$to(device = "cuda:3")
        
        # forward
        test_output <- model(batch_tokens_test)
        test_loss <- loss_fn(test_output, batch_LC50_test$unsqueeze(2L))
        test_loss_total <- test_loss_total + test_loss$item()
        
        test_r2 <- calculate_r_squared(as.matrix(batch_LC50_test), as.matrix(test_output$squeeze()))
        test_r2_total <- test_r2_total + test_r2
      }
      
      avg_test_loss <- test_loss_total / num_batches
      avg_test_r2 <- test_r2_total / num_batches
      
      test_losses <- c(test_losses, avg_test_loss)
      test_r2_scores <- c(test_r2_scores, avg_test_r2)
    })
    
    cat(sprintf("Fold: %s, Epoch [%d/%d], Train Loss: %.4f, Test R2: %.4f\n",
                fold_name, epoch, num_epochs, train_loss$item(), avg_test_r2))
    
    
  }
  
  return(list(train_losses = train_losses, test_losses = test_losses,
              train_r2_scores = train_r2_scores, test_r2_scores = test_r2_scores, model = torch::torch_serialize(model)))
}


fold_columns <- list(
  randomFold = df$randomFold,
  speciesFold = df$speciesFold,
  chemicalFold = df$chemicalFold
)

results <- list()

for (fold_name in names(fold_columns)) {
  fold_column <- fold_columns[[fold_name]]
  results[[fold_name]] <- list()
  
  for (fold_id in 1:10) {
    train_indices <- which(fold_column != fold_id)
    test_indices <- which(fold_column == fold_id) 
    
    fold_results <- train_and_evaluate(train_indices, test_indices, fold_name = paste0(fold_name, "_", fold_id), num_epochs = 1000L, batch_size = 400L)
    
    results[[fold_name]][[paste0("fold_", fold_id)]] <- fold_results
  }
}

get_r2 <- function(fold_result_list) {
  sapply(fold_result_list, function(res) {
    c(
      train_r2 = tail(res$train_r2_scores, 1),
      test_r2 = tail(res$test_r2_scores, 1)
    )
  })
}

r2 <- lapply(results, get_r2)

saveRDS(r2, file = "~/paper/EcoToxTM/Results/BiGRU_noTraits.rds")














