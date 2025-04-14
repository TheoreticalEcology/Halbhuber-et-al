##############################   Bi-GRU TRAITS   ######################################
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
data_matrix <- apply(as.matrix(df[,1:194]), 2, as.numeric) 
# Chemical Traits and Species Traits 
data_matrix_tensor <- torch_tensor(data_matrix, dtype = torch_long()) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())
traits_tensor <- torch_tensor(scale(data_matrix[,178:194]), dtype = torch_float32())

# LC50
data_matrix_LC50 <- as.numeric(df[,195]) 
LC50_tensor <- torch_tensor(data_matrix_LC50, dtype = torch_float())

#### Training and Test datasets ####
set.seed(123)
# Indices
train_indices <- sample(1:nrow(LC50_tensor), 0.8 * nrow(LC50_tensor))
test_indices <- setdiff(1:nrow(LC50_tensor), train_indices)

# Tokens
train_tokens_tensor <- tokens_tensor[train_indices, , drop = FALSE]
test_tokens_tensor <- tokens_tensor[test_indices, , drop = FALSE]
# Traits
train_traits_tensor <- traits_tensor[train_indices, , drop = FALSE]
test_traits_tensor <- traits_tensor[test_indices, , drop = FALSE]


# LC50
train_LC50_tensor <- LC50_tensor[train_indices]
test_LC50_tensor <- LC50_tensor[test_indices]


####  Training and Evaluation   ####

create_dataloader <- function(data_tensor, batch_size) {
  num_samples <- data_tensor$size(1) 
  indices <- sample(1:num_samples) 
  batches <- split(indices, ceiling(seq_along(indices) / batch_size))
  
  return(batches) 
}


train_and_evaluate <- function(train_indices, test_indices, fold_name, num_epochs = 1000L, batch_size = 128) {
  
  loss_fn <- nn_mse_loss()
  model <- BiGRU_TRAITS(embedding = 10L)
  model$to(device = "cuda:3")
  
  train_tokens_fold <- train_tokens_tensor$to(device = "cuda:3")
  train_traits_fold <- train_traits_tensor$to(device = "cuda:3")
  train_LC50_fold <- train_LC50_tensor$to(device = "cuda:3")
  
  test_tokens_fold <- test_tokens_tensor$to(device = "cuda:3")
  test_traits_fold <- test_traits_tensor$to(device = "cuda:3")
  test_LC50_fold <- test_LC50_tensor$to(device = "cuda:3")
  
  optimizer <- optim_adam(model$parameters, lr = 0.01, weight_decay = 0.00001)
  
  train_losses <- c()
  test_losses <- c()
  train_r2_scores <- c()
  test_r2_scores <- c()
  
  for (epoch in 1:num_epochs) {
    dataloader_train <- create_dataloader(train_tokens_fold, batch_size = batch_size)
    
    # Set model to training mode (enables dropout)
    model$train()
    
    for (batch_indices in dataloader_train) {
      batch_tokens_train <- train_tokens_fold[batch_indices, , drop = FALSE]
      batch_traits_train <- train_traits_fold[batch_indices, , drop = FALSE]
      batch_LC50_train <- train_LC50_fold[batch_indices]
      
      batch_tokens_train <- batch_tokens_train$to(device = "cuda:3")
      batch_traits_train <- batch_traits_train$to(device = "cuda:3")
      batch_LC50_train <- batch_LC50_train$to(device = "cuda:3")
      
      # Forward pass
      train_output <- model(batch_tokens_train, batch_traits_train)
      
      train_loss <- loss_fn(train_output, batch_LC50_train$unsqueeze(2L))
      train_losses <- c(train_losses, train_loss$item())
      
      train_r2 <- calculate_r_squared(as.matrix(batch_LC50_train), as.matrix(train_output$squeeze()))
      train_r2_scores <- c(train_r2_scores, train_r2)
      
      # Backward pass
      train_loss$backward()
      optimizer$step()
      optimizer$zero_grad()
    }
    
    # Set model to evaluation mode (disables dropout)
    model$eval()
    dataloader_test <- create_dataloader(test_tokens_fold, batch_size = batch_size)
    
    # Disable gradient calculation during evaluation
    with_no_grad({
      test_loss_total <- 0
      test_r2_total <- 0
      num_batches <- length(dataloader_test)
      
      for (batch_indices in dataloader_test) {
        batch_tokens_test <- test_tokens_fold[batch_indices, , drop = FALSE]
        batch_traits_test <- test_traits_fold[batch_indices, , drop = FALSE]
        batch_LC50_test <- test_LC50_fold[batch_indices]
        
        batch_tokens_test <- batch_tokens_test$to(device = "cuda:3")
        batch_traits_test <- batch_traits_test$to(device = "cuda:3")
        batch_LC50_test <- batch_LC50_test$to(device = "cuda:3")
        
        # Forward pass
        test_output <- model(batch_tokens_test, batch_traits_test)
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
    
    cat(sprintf("Fold: %s, Epoch [%d/%d], Train Loss: %.4f, Test Loss: %.4f\n",
                fold_name, epoch, num_epochs, train_loss$item(), avg_test_loss))
  }
  
  # save model parameters 
  final_model_weights <- lapply(1:length(model$parameters), function(l) model$parameters[[l]]$cpu() %>% as.matrix())
  saveRDS(final_model_weights, "/home/isabellehalbhuber/Toxicology/Scripts/RNN_parameters.RDS")
  
  
  return(list(train_losses = train_losses, test_losses = test_losses,
              train_r2_scores = train_r2_scores, test_r2_scores = test_r2_scores, model = model))
}



fold_columns <- list(
  randomFold = unique_df$randomFold,
  speciesFold = unique_df$speciesFold,
  chemicalFold = unique_df$chemicalFold
)

results <- list()

for (fold_name in names(fold_columns)) {
  fold_column <- fold_columns[[fold_name]]
  
  train_indices <- which(fold_column != 1)
  test_indices <- which(fold_column == 1)
  
  fold_results <- train_and_evaluate(train_indices, test_indices, fold_name, num_epochs = 1000L, batch_size =128)
  results[[fold_name]] <- fold_results
}

average_train_r2 <- sapply(results, function(res) mean(res$train_r2_scores))
average_test_r2 <- sapply(results, function(res) mean(res$test_r2_scores))



results_r_squared_rnn_with <- as.data.frame(rbind(average_train_r2, average_test_r2))

saveRDS(results_r_squared_rnn_with, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_rnn_with.rds")
