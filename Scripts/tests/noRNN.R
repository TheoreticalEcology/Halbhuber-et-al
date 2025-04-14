#################################################################################
######################   RNN with Species Traits   ##############################
#################################################################################
library(readr)
library(torch)
library(dplyr)

unique_df <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embedding_ready_df_withSpeciesCV.rds")
unique_df <- na.omit(unique_df)

#################################################################################
##############################   Attention    ###################################
#################################################################################
BiRNN <- nn_module(
  initialize = function(embedding = 40L, max_len = 177L) {
    # Embedding Layer
    self$embedding <- nn_embedding(
      num_embeddings = 425L, 
      embedding_dim = embedding, 
      sparse = FALSE, 
      scale_grad_by_freq = TRUE
    )
    
    # Positional Encoding
    self$pos_encoding <- self$generate_positional_encoding(max_len, embedding)
    
    # Multi-Head Self-Attention
    self$self_att <- torch::nn_multihead_attention(embed_dim = embedding, num_heads = 1L, dropout = 0.1, batch_first = TRUE)
    
    # Trait Embedding Layers
    self$trait_full1 <- nn_linear(in_features = 2L, out_features = 100L)
    self$dropout1 <- nn_dropout(0.3)
    self$trait_full2 <- nn_linear(in_features = 100L, out_features = 100L)
    self$dropout2 <- nn_dropout(0.3)
    self$trait_full3 <- nn_linear(in_features = 100L, out_features = embedding)
    
    # Final output layers
    self$output_layer1 <- nn_linear(in_features = embedding, out_features = 50L)
    self$output_layer2 <- nn_linear(in_features = 50L, out_features = 1L)
  },
  
  generate_positional_encoding = function(max_len, embedding_dim) {
    pos_enc <- torch_zeros(max_len, embedding_dim)
    positions <- torch_arange(1, max_len)$unsqueeze(2)
    div_term <- torch_exp(torch_arange(1, embedding_dim, 2) * (-torch_log(torch_tensor(10000)) / embedding_dim))
    
    pos_enc[, seq(1, embedding_dim, 2)] <- torch_sin(positions * div_term)
    pos_enc[, seq(2, embedding_dim, 2)] <- torch_cos(positions * div_term)
    
    return(pos_enc$unsqueeze(1))  # [max_len, 1, embedding_dim]
  },
  
  forward = function(tokens, traits) {
    #### Token Embedding
    x <- self$embedding(tokens)
    
    # Add Positional Encoding to the Token Embeddings
    seq_len <- tokens$size(2)
    pos_enc <- self$pos_encoding[1:seq_len, , ]$to(dtype = torch_float32(), device = x$device)
    x <- x + pos_enc$squeeze(1)
    
    #### Trait Embedding
    traits <- traits$unsqueeze(3L)
    traits_position <- torch_tensor(matrix(scale(1:17), 17L, 1L),  # Skalieren nicht vergessen
                                    device = traits$device,
                                    dtype = traits$dtype)$`repeat`(list(tokens$shape[1], 1L, 1L))
    traits <- torch_cat(list(traits, traits_position), dim = 3L)
    
    traits_embedded <- self$trait_full1(traits) %>% nnf_relu() %>% self$dropout1()
    traits_embedded <- self$trait_full2(traits_embedded) %>% nnf_relu() %>% self$dropout2()
    traits_embedded <- self$trait_full3(traits_embedded)
    
    #### Concatenation of Token and Trait Embeddings
    token_trait_combined <- torch_cat(list(x, traits_embedded), dim = 2L)
    
    #### Self-Attention
    q <- token_trait_combined  # Query is the concatenated output of token and trait embeddings
    k <- token_trait_combined  # Key
    v <- token_trait_combined  # Value
    
    #### Masking for Padding Tokens
    padding_mask <- tokens$clone()
    mask <- padding_mask == 1L  # Assuming padding tokens are marked with 1
    mask <- torch_cat(list(mask, torch_zeros(x$shape[1], traits$shape[2], dtype = torch_bool(), device = x$device)), dim = 2)
    
    # Apply Multi-Head Attention
    attention_out <- self$self_att(q, k, v, key_padding_mask = !mask)
    attn_out <- torch_mean(attention_out[[1]], dim = 2)  # [batch_size, embedding_dim]
    
    #### Token Feature NN
    output_token <- self$output_layer1(attn_out) %>% nnf_relu()
    
    #### Combine outputs
    output <- self$output_layer2(output_token)
    
    return(output)
  }
)



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
  train_traits_fold = train_traits_tensor$to(device = "cuda:3")
  train_LC50_fold = train_LC50_tensor$to(device = "cuda:3")
  
  test_tokens_fold = test_tokens_tensor$to(device = "cuda:3")
  test_traits_fold = test_traits_tensor$to(device = "cuda:3")
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
      batch_traits_train <- train_traits_fold[batch_indices, , drop = FALSE]
      batch_LC50_train <- train_LC50_fold[batch_indices]
      
      batch_tokens_train = batch_tokens_train$to(device = "cuda:3")
      batch_traits_train = batch_traits_train$to(device = "cuda:3")
      batch_LC50_train = batch_LC50_train$to(device = "cuda:3")
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
        batch_traits_test <- test_traits_fold[batch_indices, , drop = FALSE]
        batch_LC50_test <- test_LC50_fold[batch_indices]
        
        batch_tokens_test = batch_tokens_test$to(device = "cuda:3")
        batch_traits_test = batch_traits_test$to(device = "cuda:3")
        batch_LC50_test = batch_LC50_test$to(device = "cuda:3")
        
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
  # }
  
  final_model_weights <- lapply(1:length(model$parameters), function(l) model$parameters[[l]]$cpu() %>% as.matrix())
  saveRDS(final_model_weights, "/home/isabellehalbhuber/Toxicology/Scripts/Final_TRAIT_model_weights_noRNN.RDS")
  
  return(list(train_losses = train_losses, test_losses = test_losses,
              train_r2_scores = train_r2_scores, test_r2_scores = test_r2_scores))
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



results_r_squared_Attention_with <- as.data.frame(cbind(average_train_r2, average_test_r2))

saveRDS(results_r_squared_Attention_with, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_Attention_with_noRNN.rds")

