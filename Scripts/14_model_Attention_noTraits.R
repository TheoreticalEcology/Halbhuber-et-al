#################################################################################
######################   RNN with Species Traits   ##############################
#################################################################################

unique_df <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embedding_ready_df_withSpeciesCV.rds")
unique_df <- na.omit(unique_df)

#embeddedTokens <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embeddedTokens.rds")
library(torch)


#################################################################################
#################################   RNN    ######################################
#################################################################################

BiRNN <- nn_module(
  initialize = function(embedding = 40L, rnn_hidden_size = 20L, rnn_layers = 3L) {
    self$embedding <- nn_embedding(
      num_embeddings = 425L, 
      embedding_dim = embedding, 
      sparse = FALSE, 
      scale_grad_by_freq = TRUE
    )
    
    self$rnn <- nn_gru(
      input_size = embedding, 
      hidden_size = rnn_hidden_size, 
      num_layers = rnn_layers, 
      bidirectional = TRUE,
      batch_first = TRUE,
      dropout = 0
    )
    
    self$self_att = torch::nn_multihead_attention(embed_dim = 40L, num_heads = 1L, dropout = 0.1, batch_first = TRUE)
    
    # Self-Attention QKV Layer
    self$qkv_size <- rnn_hidden_size * 2  
    self$q_linear <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = self$qkv_size)
    self$k_linear <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = self$qkv_size)
    self$v_linear <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = self$qkv_size)
    
    # Final linear layers for token features
    self$token_full <- nn_linear(in_features = self$qkv_size, out_features = 20L)
    
    # Combined output layers
    self$output_layer1 <- nn_linear(in_features = 20L, out_features = 50L)
    self$output_layer2 <- nn_linear(in_features = 50L, out_features = 1L)
  },
  
  forward = function(tokens) {
    #### Token RNN
    x <- self$embedding(tokens)
    rnn_out <- self$rnn(x)[[1]]  # [batch_size, seq_len, 2 * rnn_hidden_size]
    
    # QKV for Self-Attention
    q <- self$q_linear(rnn_out)  # [batch_size, seq_len, qkv_size]
    k <- self$k_linear(rnn_out)  # [batch_size, seq_len, qkv_size]
    v <- self$v_linear(rnn_out)  # [batch_size, seq_len, qkv_size]
    
#     attn_weights <- nnf_softmax(torch_bmm(q, k$transpose(2, 3)) / sqrt(self$qkv_size), dim = -1)  # [batch_size, seq_len, seq_len]
#     attn_out <- torch_bmm(attn_weights, v)  # [batch_size, seq_len, qkv_size]
#     
#     self$attn_weights <- attn_weights
#     
#     # Take the mean along the sequence length dimension to summarize the attention output
#     attn_out <- torch_mean(attn_out, dim = 2)  # [batch_size, qkv_size]
#     
#     # Token Feature NN
#     output_token <- self$token_full(attn_out) %>% nnf_relu()
#     
#     #### Combine outputs
#     output <- self$output_layer1(output_token) %>% nnf_relu()
#     output <- self$output_layer2(output)
#     
#     return(output)
#   }
# )

    #### Masking for Padding Tokens
    padding_mask <- tokens$clone()
    #print(tokens$shape)
    
    mask <- padding_mask == 1L  # Assuming padding tokens are marked with 1 [batchsize, 177]
    mask <- mask$to(dtype=torch_bool(), device=tokens$device)
    
  
    #print(rnn_out$shape)
    #print(traits_embedded$shape)
    #print(rnn_out_traits$shape)
    #print(attn_out$shape)
    
    #### Multi-Head Self-Attention
    self$mask = mask
    attention_out <- self$self_att(q, k, v, key_padding_mask = !mask )
    self$att_weights = attention_out[[2]]
    # Take the mean along the sequence length dimension to summarize the attention output
    attn_out <- torch_mean(attention_out[[1]], dim = 2)  # [batch_size, embed_dim]
    #### Token Feature NN
    output_token <- self$token_full(attn_out) %>% nnf_relu()
    #### Combine outputs
    output <- self$output_layer1(output_token) %>% nnf_relu()
    output <- self$output_layer2(output)
    
    return(output)
  }
)



#mask <- mask$unsqueeze(2)$unsqueeze(3) 
    
    

# model = BiRNN(embedding = 10L)
# #output <- model(train_tokens_fold[ind, ], train_traits_fold[ind, ])
# 
# model$parameters[[1]]
# # turn off training of embeddings
# model$parameters[[1]]$requires_grad_(FALSE)
# # load pretrained embeddings into embedding layer
# embd = data.matrix(embeddedTokens)
# #embd = matrix(1.0, 425, 10) # mit vortrainierten embeddings von GloVe ersetzen
# model$parameters[[1]]$set_data(embd)
# model$parameters[[1]]
# # train model

#################################################################################
###############################   Input Data   ##################################
#################################################################################

# Chemical and Species Traits
data_matrix <- apply(as.matrix(unique_df[,2:195]), 2, as.numeric) 
data_matrix_tensor <- torch_tensor(data_matrix, dtype = torch_long()) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())

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
# LC50
train_LC50_tensor <- LC50_tensor[train_indices]
test_LC50_tensor <- LC50_tensor[test_indices]

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
# 
# create_dataloader <- function(data_tensor, batch_size) {
#   num_samples <- length((data_tensor)[,1])
#   indices <- sample(1:num_samples)  
#   batches <- split(indices, ceiling(seq_along(indices) / batch_size))
#   
#   return(batches) 
# }
# 
# #batches <- create_dataloader(train_tokens_tensor, 10)
# 
# calculate_r2 <- function(y_true, y_pred) {
#   res_ss <- sum((y_true - y_pred)^2)
#   total_ss <- sum((y_true - mean(y_true))^2)
#   r2 <- 1 - res_ss / total_ss
#   return(r2)
# }
# 
# train_and_evaluate <- function(train_indices, test_indices, fold_name, num_epochs = 100L, batch_size = 500) {
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
#   optimizer <- optim_adam(model$parameters, lr = 0.01, weight_decay = 0.0001)
#   
#   train_losses <- c()
#   test_losses <- c()
#   train_r2_scores <- c()
#   test_r2_scores <- c()
#   
#   for (epoch in 1:num_epochs) {
#     
#     #train_dataset <- tensor_dataset(list(train_tokens_fold, train_traits_fold, train_LC50_fold))
#   
#     dataloader_train <- create_dataloader(train_tokens_fold, batch_size = batch_size)
#     model$train()
#     
#     for (batch_indices in dataloader_train) {
#       
#       batch_tokens <- train_tokens_fold[batch_indices, , drop = FALSE]
#       batch_traits <- train_traits_fold[batch_indices, , drop = FALSE]
#       batch_LC50 <- train_LC50_fold[batch_indices]
#       
#       batch_tokens = batch_tokens$to(device = "cuda:3")
#       batch_traits = batch_traits$to(device = "cuda:3")
#       batch_LC50 = batch_LC50$to(device = "cuda:3")
#       
#       
#       train_output <- model(batch_tokens, batch_traits)
#       train_loss <- loss_fn(train_output, batch_LC50$unsqueeze(2L))
#       train_losses <- c(train_losses, train_loss$item())
#       
#       train_r2 <- calculate_r2(as.matrix(batch_LC50), as.matrix(train_output$squeeze()))
#       train_r2_scores <- c(train_r2_scores, train_r2)
#       
# 
#       train_loss$backward()
#       optimizer$step()
#       optimizer$zero_grad()
#     }
#   
#     
#     
#     # Evaluation
#     model$eval()
#     with_no_grad({
#       test_output <- model(test_tokens_fold, test_traits_fold)
#       test_loss <- loss_fn(test_output, test_LC50_fold$unsqueeze(2L))
#       test_losses <- c(test_losses, test_loss$item())
#       
#       # Calculate R2 score for test data
#       test_r2 <- calculate_r2(as.matrix(test_LC50_fold), as.matrix(test_output$squeeze()))
#       test_r2_scores <- c(test_r2_scores, test_r2)
#     })
#     
#     cat(sprintf("Fold: %s, Epoch [%d/%d], Train Loss: %.4f, Test Loss: %.4f\n",
#                 fold_name, epoch, num_epochs, train_loss$item(), test_loss$item()))
#   }
#   
#   return(list(train_losses = train_losses, test_losses = test_losses,
#               train_r2_scores = train_r2_scores, test_r2_scores = test_r2_scores, model = model))
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
#   fold_results <- train_and_evaluate(train_indices, test_indices, fold_name, num_epochs = 100L, batch_size = 500L)
#   results[[fold_name]] <- fold_results
# }
# 
# 
# # random_CV_train_RNN_with <- results[["randomFold"]][[1]]$train_r2_scores  
# # ranom_CV_test_RNN_with <- results[["randomFold"]][[1]]$test_r2_scores   
# # chemicals_blocked_CV_train_RNN <-results[["chemicalFold"]][[1]]$train_r2_scores  
# # chemicals_blocked_CV_test_RNN <- results[["chemicalFold"]][[1]]$test_r2_scores   
# # species_blocked_CV_train_RNN <- results[["speciesFold"]][[1]]$train_r2_scores  
# # species_blocked_CV_test_RNN <- results[["speciesFold"]][[1]]$test_r2_scores  
# 
# average_train_r2 <- sapply(results, function(res) mean(res$train_r2_scores))
# average_test_r2 <- sapply(results, function(res) mean(res$test_r2_scores))
# 
# fold_results <- results[[fold_name]]
# 
# results_r_squared_with_speciesTraits_RNN_attention <- as.data.frame(cbind(average_train_r2, average_test_r2))
# 
# saveRDS(results_r_squared_with_speciesTraits_RNN_attention, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_with_speciesTraits_RNN_attention.rds") 
# saveRDS(fold_results, file = "/home/isabellehalbhuber/Toxicology/Results/selfAttention_fold_results.rds") 
# 
# 
# 
# 
# 
# 
# 




# library(torch)
# calculate_r2 <- function(y_true, y_pred) {
#   res_ss <- sum((y_true - y_pred)^2)
#   total_ss <- sum((y_true - mean(y_true))^2)
#   r2 <- 1 - res_ss / total_ss
#   return(r2)
# }
# 
# 
# create_dataloader <- function(data_tensor, batch_size) {
#   num_samples <- data_tensor$size(1) 
#   indices <- sample(1:num_samples) 
#   batches <- split(indices, ceiling(seq_along(indices) / batch_size))
#   
#   return(batches) 
# }
# 
# 
# train_and_evaluate <- function(train_indices, test_indices, fold_name, num_epochs = 5000L, batch_size = 256) {
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
#   optimizer <- optim_adam(model$parameters, lr = 0.01, weight_decay = 0.0001)
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
#       batch_tokens <- train_tokens_fold[batch_indices, , drop = FALSE]
#       batch_traits <- train_traits_fold[batch_indices, , drop = FALSE]
#       batch_LC50 <- train_LC50_fold[batch_indices]
#       
#       # Forward pass
#       train_output <- model(batch_tokens, batch_traits)
#       
#       train_loss <- loss_fn(train_output, batch_LC50$unsqueeze(2L))
#       train_losses <- c(train_losses, train_loss$item())
# 
#       train_r2 <- calculate_r2(as.matrix(batch_LC50), as.matrix(train_output$squeeze()))
#       train_r2_scores <- c(train_r2_scores, train_r2)
#     
#       # Backward pass
#       train_loss$backward()
#       optimizer$step()
#       optimizer$zero_grad()
#     }
#     
#     # Evaluation after every 10th epoch
#     if (epoch %% 10 == 0) {
#       model$eval()
#       with_no_grad({
#         test_output <- model(test_tokens_fold, test_traits_fold)
#         test_loss <- loss_fn(test_output, test_LC50_fold$unsqueeze(2L))
#         test_losses <- c(test_losses, test_loss$item())
# 
#         test_r2 <- calculate_r2(as.matrix(test_LC50_fold), as.matrix(test_output$squeeze()))
#         test_r2_scores <- c(test_r2_scores, test_r2)
#       })
#       
#       cat(sprintf("Fold: %s, Epoch [%d/%d], Train Loss: %.4f, Test Loss: %.4f\n",
#                   fold_name, epoch, num_epochs, train_loss$item(), test_loss$item()))
#     }
#   }
#   
#   return(list(train_losses = train_losses, test_losses = test_losses,
#               train_r2_scores = train_r2_scores, test_r2_scores = test_r2_scores, model = model))
# }
# 
# 
# fold_columns <- list(
#   randomFold = unique_df$randomFold,
#   speciesFold = unique_df$speciesFold,
#   chemicalFold = unique_df$chemicalFold
# )
# 
# results <- list()
# 
# for (fold_name in names(fold_columns)) {
#   fold_column <- fold_columns[[fold_name]]
#   
#   train_indices <- which(fold_column != 1)
#   test_indices <- which(fold_column == 1)
#   
#   fold_results <- train_and_evaluate(train_indices, test_indices, fold_name, num_epochs = 5000L, batch_size =256)
#   results[[fold_name]] <- fold_results
# }
# 
# average_train_r2 <- sapply(results, function(res) mean(res$train_r2_scores))
# average_test_r2 <- sapply(results, function(res) mean(res$test_r2_scores))
# 
# fold_results <- results[[fold_name]]
# 
# results_r_squared_with_speciesTraits_RNN_attention <- as.data.frame(cbind(average_train_r2, average_test_r2))
# 
# saveRDS(results_r_squared_with_speciesTraits_RNN_attention, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_with_speciesTraits_RNN_attention.rds") 
# saveRDS(fold_results, file = "/home/isabellehalbhuber/Toxicology/Results/selfAttention_fold_results.rds") 
# 
# 
# 
# 

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

train_and_evaluate <- function(train_indices, test_indices, fold_name, num_epochs = 1000L, batch_size = 128) {
  
  loss_fn <- nn_mse_loss()
  model <- BiRNN(embedding = 10L)
  model$to(device = "cuda:3")
  
  train_tokens_fold = train_tokens_tensor$to(device = "cuda:3")
   train_LC50_fold = train_LC50_tensor$to(device = "cuda:3")
  
  test_tokens_fold = test_tokens_tensor$to(device = "cuda:3")
  test_LC50_fold = test_LC50_tensor$to(device = "cuda:3")
  
  optimizer <- optim_adam(model$parameters, lr = 0.001, weight_decay = 0.00001)
  
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
      
      batch_tokens_train = batch_tokens_train$to(device = "cuda:3")
      batch_LC50_train = batch_LC50_train$to(device = "cuda:3")
      # Forward pass
      train_output <- model(batch_tokens_train)
      
      train_loss <- loss_fn(train_output, batch_LC50_train$unsqueeze(2L))
      train_losses <- c(train_losses, train_loss$item())
      
      train_r2 <- calculate_r2(as.matrix(batch_LC50_train), as.matrix(train_output$squeeze()))
      train_r2_scores <- c(train_r2_scores, train_r2)
      
      # Backward pass
      train_loss$backward()
      optimizer$step()
      optimizer$zero_grad()
    }
    
    # Evaluation after every 10th epoch
    #if (epoch %% 10 == 0) {
    model$eval()
    dataloader_test <- create_dataloader(test_tokens_fold, batch_size = batch_size)
    
    with_no_grad({
      test_loss_total <- 0
      test_r2_total <- 0
      num_batches <- length(dataloader_test)
      
      for (batch_indices in dataloader_test) {
        batch_tokens_test <- test_tokens_fold[batch_indices, , drop = FALSE]
        batch_LC50_test <- test_LC50_fold[batch_indices]
        
        batch_tokens_test = batch_tokens_test$to(device = "cuda:3")
        batch_LC50_test = batch_LC50_test$to(device = "cuda:3")
        
        test_output <- model(batch_tokens_test)
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
  # }
  
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

# > saveRDS(lapply(1:length(mm), function(l) mm$parameters[[l]] %>% as.matrix()), "My_model_weights.RDS")
# > loaded_weights = readRDS("My_model_weights.RDS")
# > mm = BiRNN$new()
# > mm$parameters[[1]]
# > loaded_weights[[1]]
# > mm$parameters[[1]]
# > mm$parameters[[1]]$set_data(new_data = loaded_weights[[1]])
# > mm$parameters[[1]]
# > loaded_weights[[1]] %>% head()




















results_r_squared_Attention_without <- as.data.frame(cbind(average_train_r2, average_test_r2))

saveRDS(results_r_squared_Attention_without, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_Attention_without.rds")

