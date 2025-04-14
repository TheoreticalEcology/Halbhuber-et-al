#################################################################################
############################  NO ATTENTION  #####################################
#################################################################################

#################################################################################
######################   RNN with Species Traits   ##############################
#################################################################################

unique_df <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embedding_ready_df_withSpeciesCV.rds")
unique_df <- na.omit(unique_df)
## gepaddete sequenz umdrehen, sodass die wichtigen informationen hinten liegen (dann backpropagiation besser)
unique_df[, 2:178] <- lapply(unique_df[, 2:178], as.numeric)

remaining_columns_traits <- unique_df[c(179:ncol(unique_df))]
remaining_columns_smiles <- unique_df[1]
reversed_columns <- unique_df[, 2:178][, ncol(unique_df[, 2:178]):1]

unique_df <- cbind(remaining_columns_smiles, reversed_columns, remaining_columns_traits)


library(torch)

#################################################################################
#################################   RNN    ######################################
#################################################################################
BiRNN <- nn_module(
  initialize = function(embedding = 1L, rnn_hidden_size = 20L, rnn_layers = 3L, dropout_prob = 0.3) {
    self$embedding <- nn_embedding(num_embeddings = 425L, 
                                   embedding_dim = embedding, 
                                   sparse = FALSE, 
                                   scale_grad_by_freq = TRUE) # [batchsize, time, embedding]
    
    # RNN with Dropout
    self$rnn <- nn_gru(input_size = embedding, 
                       hidden_size = rnn_hidden_size, 
                       num_layers = rnn_layers, 
                       bidirectional = TRUE,
                       dropout = dropout_prob,  # Dropout between RNN layers
                       batch_first = TRUE)
    
    # Token NN after RNN
    self$token_full = nn_linear(in_features = 2*rnn_hidden_size, out_features = 20L)
    self$dropout_token = nn_dropout(p = dropout_prob)  # Dropout after linear layer
    
    # Species Trait NN
    self$trait_full1 = nn_linear(in_features = 17L, out_features = 20L)
    self$dropout_trait1 = nn_dropout(p = dropout_prob)  # Dropout after linear layer
    self$trait_full2 = nn_linear(in_features = 20L, out_features = 3L)
    
    # Combined output layers
    self$output_layer1 <- nn_linear(in_features = 20L+3L, out_features = 50L)
    self$dropout_output1 = nn_dropout(p = dropout_prob)  # Dropout after linear layer
    self$output_layer2 <- nn_linear(in_features = 50L, out_features = 1L)
  },
  
  forward = function(tokens, traits) {
    #### Token RNN
    x <- self$embedding(tokens)
    x <- self$rnn(x)[[1]]  # only the first element, which is the output
    output_token <- self$token_full(x[, dim(x)[2], ]) %>% nnf_selu()
    output_token <- self$dropout_token(output_token)  # Apply dropout
    
    #### Trait NN
    x <- self$trait_full1(traits) %>% nnf_selu()
    x <- self$dropout_trait1(x)  # Apply dropout
    output_traits <- self$trait_full2(x) %>% nnf_selu()  # [batch_size, 3]
    
    #### Combine outputs
    input <- torch::torch_cat(list(output_token, output_traits), dim = 2L)  # [batchsize, rnn_hidden_size + 3L]
    output <- self$output_layer1(input) %>% nnf_selu()
    output <- self$dropout_output1(output)  # Apply dropout
    output <- self$output_layer2(output)
    
    return(output)
  }
)


# BiRNN <- nn_module(
#   initialize = function(embedding = 1L, rnn_hidden_size = 20L, rnn_layers = 3L) {
#     self$embedding <- nn_embedding(num_embeddings = 425L, 
#                                    embedding_dim = embedding, 
#                                    sparse = FALSE, 
#                                    scale_grad_by_freq = TRUE) # [batchsize, time, embedding]
#     
# 
#     self$rnn <- nn_gru(input_size = embedding, 
#                        hidden_size = rnn_hidden_size, 
#                        num_layers = rnn_layers, 
#                        bidirectional =TRUE,
#                        batch_first = TRUE)
#     
#     
#     # Token NN after RNN
#     self$token_full = nn_linear(in_features = 2*rnn_hidden_size, out_features = 20L)
#     
#     
#     # Species Trait NN 
#     self$trait_full1 = nn_linear(in_features = 17L, out_features = 20L)
#     self$trait_full2 = nn_linear(in_features = 20L, out_features = 3L)
#     
#     
#     # Combined output layers
#     self$output_layer1 <- nn_linear(in_features = 20L+3L, out_features = 50L) 
#     self$output_layer2 <- nn_linear(in_features = 50L, out_features = 1L) 
#   },
#   
#   forward = function(tokens, traits) {
#     #### Token RNN
#     x <- self$embedding(tokens)
#     x <- self$rnn(x)[[1]] # only the first element, which is the output
#     output_token <- self$token_full(x[,dim(x)[2], ]) %>% nnf_selu() %>% torch_squeeze()
#     
#     #### Trait NN
#     x <- self$trait_full1(traits) %>% nnf_selu()
#     output_traits <- self$trait_full2(x) %>% nnf_selu() # [batch_size, 3]
#     
#     ##
#     #print(output_token$shape)
#     #print(output_traits$shape)
#     #### Combine outputs
#     input <- torch::torch_cat(list(output_token, output_traits), dim = 2L) # [batchsize, rnn_hidden_size + 3L]
#     output <- self$output_layer2(self$output_layer1(input) %>% nnf_selu())
#     return(output)
#   }
# )
# 
# 


#################################################################################
##################################   DNN   ######################################
#################################################################################
# 
# DNN <- nn_module(
#   initialize = function(embedding = 1L) {
#     self$embedding <- nn_embedding(425L, embedding, sparse = FALSE,scale_grad_by_freq = TRUE) # [batchsize, time, embedding]
#     self$full1 = nn_linear(in_features = embedding, out_features = 20L)
#     self$full2 = nn_linear(in_features = 20L, out_features = 1L)
#     self$trait_full1 = nn_linear(in_features = 17L, out_features = 20L)
#     self$trait_full2 = nn_linear(in_features = 20L, out_features = 3L)
#     
#     self$output_layer1 <- nn_linear(in_features = 177L+3L, out_features = 50L) 
#     self$output_layer2 <- nn_linear(in_features = 50L, out_features = 1L) 
#   },
#   
#   forward = function(tokens, traits) {
#     #### Token NN
#     x <- self$embedding(tokens) 
#     x = self$full1(x) %>% nnf_selu()
#     output_token = self$full2(x) %>% nnf_selu() %>% torch_squeeze()# [batch_size, 177]
#     
#     #### Trait NN
#     x = self$trait_full1(traits) %>% nnf_selu()
#     output_traits = self$trait_full2(x) %>% nnf_selu() # [batch_size, 8]
#     
#     #print(output_token$shape)
#     #print(output_traits$shape)
#     
#     #### Combine outputs
#     input = torch::torch_cat(list(output_token, output_traits), dim = 2L) # [batchsize, 177+8L]
#     output = self$output_layer2( self$output_layer1(input) %>% nnf_selu() )
#     return(output)
#     
#   }
# )


#################################################################################
###############################   Input Data   ##################################
#################################################################################


# Chemical and Species Traits
data_matrix <- apply(as.matrix(unique_df[,2:195]), 2, as.numeric) 
data_matrix_tensor <- torch_tensor(data_matrix, dtype = torch_long()) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())
traits_tensor <- torch_tensor(scale(data_matrix[,178:194]), dtype = torch_float32())

# LC50
data_matrix_LC50 <- as.numeric(unique_df[,196]) 
LC50_tensor <- torch_tensor(data_matrix_LC50, dtype = torch_float())

# Training and Test datasets
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

# # Chemical and Species Traits
# data_matrix <- apply(as.matrix(unique_df[,2:196]), 2, as.numeric) 
# data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
# tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())
# traits_tensor <- torch_tensor(scale(data_matrix[,178:194]), dtype = torch_float32())
# 
# # LC50
# data_matrix_LC50 <- as.numeric(data_matrix[,195]) 
# LC50_tensor <- torch_tensor(data_matrix_LC50, dtype = torch_float())
# 
# # Training and Test datasets
# set.seed(123)
# # Indices
# train_indices <- sample(1:nrow(tokens_tensor), 0.9 * nrow(tokens_tensor))
# test_indices <- setdiff(1:nrow(tokens_tensor), train_indices)
# 
# # Tokens
# train_tokens_tensor <- tokens_tensor[train_indices, , drop = FALSE]
# test_tokens_tensor <- tokens_tensor[test_indices, , drop = FALSE]
# # Traits
# train_traits_tensor <- traits_tensor[train_indices, , drop = FALSE]
# test_traits_tensor <- traits_tensor[test_indices, , drop = FALSE]
# # LC50
# train_LC50_tensor <- LC50_tensor[train_indices]
# test_LC50_tensor <- LC50_tensor[test_indices]
# 
# # Train Folds
# train_tokens_fold <- tokens_tensor[train_indices, , drop = FALSE]
# train_traits_fold <- traits_tensor[train_indices, , drop = FALSE]
# train_LC50_fold <- LC50_tensor[train_indices]
# 
# #Test Folds
# test_tokens_fold <- tokens_tensor[test_indices, , drop = FALSE]
# test_traits_fold <- traits_tensor[test_indices, , drop = FALSE]
# test_LC50_fold <- LC50_tensor[test_indices]



#################################################################################
########################   Training and Evaluation   ############################
#################################################################################

# calculate_r2 <- function(y_true, y_pred) {
#   res_ss <- sum((y_true - y_pred)^2)
#   total_ss <- sum((y_true - mean(y_true))^2)
#   r2 <- 1 - res_ss / total_ss
#   return(r2)
# }
# 
# train_and_evaluate <- function(train_indices, test_indices, fold_name, num_epochs = 1000L) {
#   
#   model <- BiRNN(embedding = 10L)
#   model$to(device = "cuda:3")
#   train_tokens_fold = train_tokens_tensor$to(device = "cuda:3")
#   train_traits_fold = train_traits_tensor$to(device = "cuda:3")
#   train_LC50_fold = (train_LC50_tensor$to(device = "cuda:3"))
#   test_tokens_fold = test_tokens_tensor$to(device = "cuda:3")
#   test_traits_fold = test_traits_tensor$to(device = "cuda:3")
#   test_LC50_fold = (test_LC50_tensor$to(device = "cuda:3"))
#   
# 
#   loss_fn <- nn_mse_loss()
#   model <- BiRNN(embedding = 10L)
#   model$to(device = "cuda:3")
#   optimizer <- optim_adam(model$parameters, lr = 0.01, weight_decay = 0.00001)
#   
# 
#   train_losses <- c()
#   test_losses <- c()
#   train_r2_scores <- c()
#   test_r2_scores <- c()
#   
#   for (epoch in 1:num_epochs) {
#     model$train()
# 
#     sample_size <- min(300L, length(train_indices))
#     ind <- sample(length(train_indices), sample_size)
#     
#     # Forward propagation
#     train_output <- model(train_tokens_fold[ind, , drop = FALSE], train_traits_fold[ind, , drop = FALSE])$relu()
#     
#     train_loss <- loss_fn(train_output, train_LC50_fold[ind]$unsqueeze(2L))
#     train_losses <- c(train_losses, train_loss$item())
# 
#     train_r2 <- calculate_r2(as.matrix(train_LC50_fold[ind]), as.matrix(train_output$squeeze()))
#     train_r2_scores <- c(train_r2_scores, train_r2)
#     
#     # Backward propagation
#     train_loss$backward()
#     optimizer$step()
#     optimizer$zero_grad()
#     
#     # Evaluation
#     model$eval()
#     with_no_grad({
#       test_output <- model(test_tokens_fold, test_traits_fold)$relu()
#       
#       test_loss <- loss_fn(test_output, test_LC50_fold$unsqueeze(2L))
#       test_losses <- c(test_losses, test_loss$item())
#       plot(as.matrix(test_output), as.matrix(test_LC50_tensor))
#       
#       test_r2 <- calculate_r2(as.matrix(test_LC50_fold), as.matrix(test_output$squeeze()))
#     })
#     test_r2_scores <- c(test_r2_scores, test_r2)
#     
#     cat(sprintf("Fold: %s, Epoch [%d/%d], Train Loss: %.4f, Test Loss:%.4f\n", 
#                 fold_name, epoch, num_epochs, train_loss$item(), test_loss$item()))
#   }
#   
#   return(list(train_losses = train_losses, test_losses = test_losses, 
#               train_r2_scores = train_r2_scores, test_r2_scores = test_r2_scores))
# }
# 
# 
# fold_columns <- list(
#   randomFold = unique_df$randomFold,
#   speciesFold = unique_df$speciesFold,
#   chemicalFold = unique_df$chemicalFold
# )
# results <- list()
# 
# for (fold_name in names(fold_columns)) {
#   fold_column <- fold_columns[[fold_name]]
#  
#   train_indices <- which(fold_column != 1)
#   test_indices <- which(fold_column == 1)
#   
#   fold_results <- train_and_evaluate(train_indices, test_indices, fold_name)
#   results[[fold_name]] <- fold_results
# }
# 
# 
# average_train_r2 <- sapply(results, function(res) mean(res$train_r2_scores))
# average_test_r2 <- sapply(results, function(res) mean(res$test_r2_scores))

library(torch)
calculate_r2 <- function(y_true, y_pred) {
  res_ss <- sum((y_true - y_pred)^2)
  total_ss <- sum((y_true - mean(y_true))^2)
  r2 <- 1 - res_ss / total_ss
  return(r2)
}

create_dataloader <- function(data_tensor, batch_size) {
  num_samples <- data_tensor$size(1) 
  indices <- sample(1:num_samples) 
  batches <- split(indices, ceiling(seq_along(indices) / batch_size))
  
  return(batches) 
}

# train_and_evaluate <- function(train_indices, test_indices, fold_name, num_epochs = 1000L, batch_size = 128) {
#   
#   loss_fn <- nn_mse_loss()
#   model <- BiRNN(embedding = 10L)
#   model$to(device = "cuda:3")
#   
#   train_tokens_fold = train_tokens_tensor$to(device = "cuda:3")
#   train_traits_fold = train_traits_tensor$to(device = "cuda:3")
#   train_LC50_fold = train_LC50_tensor$to(device = "cuda:3")
#   
#   test_tokens_fold = test_tokens_tensor$to(device = "cuda:3")
#   test_traits_fold = test_traits_tensor$to(device = "cuda:3")
#   test_LC50_fold = test_LC50_tensor$to(device = "cuda:3")
#   
#   optimizer <- optim_adam(model$parameters, lr = 0.01, weight_decay = 0.00001)
#   
#   train_losses <- c()
#   test_losses <- c()
#   train_r2_scores <- c()
#   test_r2_scores <- c()
#   
#   for (epoch in 1:num_epochs) {
#     dataloader_train <- create_dataloader(train_tokens_fold, batch_size = batch_size)
#     
#     model$train() 
#     
#     for (batch_indices in dataloader_train) {
#       batch_tokens_train <- train_tokens_fold[batch_indices, , drop = FALSE]
#       batch_traits_train <- train_traits_fold[batch_indices, , drop = FALSE]
#       batch_LC50_train <- train_LC50_fold[batch_indices]
#       
#       batch_tokens_train = batch_tokens_train$to(device = "cuda:3")
#       batch_traits_train = batch_traits_train$to(device = "cuda:3")
#       batch_LC50_train = batch_LC50_train$to(device = "cuda:3")
#       # Forward pass
#       train_output <- model(batch_tokens_train, batch_traits_train)
#       
#       train_loss <- loss_fn(train_output, batch_LC50_train$unsqueeze(2L))
#       train_losses <- c(train_losses, train_loss$item())
#       
#       train_r2 <- calculate_r2(as.matrix(batch_LC50_train), as.matrix(train_output$squeeze()))
#       train_r2_scores <- c(train_r2_scores, train_r2)
#       
#       # Backward pass
#       train_loss$backward()
#       optimizer$step()
#       optimizer$zero_grad()
#     }
#     
#     # Evaluation after every 10th epoch
#     #if (epoch %% 10 == 0) {
#     model$eval()
#     dataloader_test <- create_dataloader(test_tokens_fold, batch_size = batch_size)
#     
#     with_no_grad({
#       test_loss_total <- 0
#       test_r2_total <- 0
#       num_batches <- length(dataloader_test)
#       
#       for (batch_indices in dataloader_test) {
#         batch_tokens_test <- test_tokens_fold[batch_indices, , drop = FALSE]
#         batch_traits_test <- test_traits_fold[batch_indices, , drop = FALSE]
#         batch_LC50_test <- test_LC50_fold[batch_indices]
#         
#         batch_tokens_test = batch_tokens_test$to(device = "cuda:3")
#         batch_traits_test = batch_traits_test$to(device = "cuda:3")
#         batch_LC50_test = batch_LC50_test$to(device = "cuda:3")
#         
#         test_output <- model(batch_tokens_test, batch_traits_test)
#         test_loss <- loss_fn(test_output, batch_LC50_test$unsqueeze(2L))
#         test_loss_total <- test_loss_total + test_loss$item()
#         
#         test_r2 <- calculate_r2(as.matrix(batch_LC50_test), as.matrix(test_output$squeeze()))
#         test_r2_total <- test_r2_total + test_r2
#       }
#       
#       avg_test_loss <- test_loss_total / num_batches
#       avg_test_r2 <- test_r2_total / num_batches
#       
#       test_losses <- c(test_losses, avg_test_loss)
#       test_r2_scores <- c(test_r2_scores, avg_test_r2)
#     })
#     
#     cat(sprintf("Fold: %s, Epoch [%d/%d], Train Loss: %.4f, Test Loss: %.4f\n",
#                 fold_name, epoch, num_epochs, train_loss$item(), avg_test_loss))
#   }
#   # }
#   
#   return(list(train_losses = train_losses, test_losses = test_losses,
#               train_r2_scores = train_r2_scores, test_r2_scores = test_r2_scores, model = model))
# }
# 

train_and_evaluate <- function(train_indices, test_indices, fold_name, num_epochs = 1000L, batch_size = 128) {
  
  loss_fn <- nn_mse_loss()
  model <- BiRNN(embedding = 10L)
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
      
      train_r2 <- calculate_r2(as.matrix(batch_LC50_train), as.matrix(train_output$squeeze()))
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
        
        test_r2 <- calculate_r2(as.matrix(batch_LC50_test), as.matrix(test_output$squeeze()))
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

saveRDS(results_r_squared_rnn_with, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_rnn_with_modelSave.rds")
