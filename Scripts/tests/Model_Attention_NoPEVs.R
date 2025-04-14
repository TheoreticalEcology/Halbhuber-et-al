#################################################################################
######################   RNN with Species Traits   ##############################
#################################################################################
library(readr)
library(torch)
library(dplyr)

device = "cuda:3"

unique_df <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embedding_ready_df_withSpeciesCV.rds")
unique_df <- na.omit(unique_df)
#embeddedTokens <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embeddedTokens.rds")
#################################################################################
##############################   Attention    ###################################
#################################################################################
BiRNN <- nn_module(
  initialize = function(embedding = 40L, rnn_hidden_size = 20L, rnn_layers = 3L, max_len = 177L) {
    self$embedding <- nn_embedding(
      num_embeddings = 425L, 
      embedding_dim = embedding, 
      sparse = FALSE, 
      scale_grad_by_freq = TRUE
    )
    
    # self$pos_encoding <- self$generate_positional_encoding(max_len, embedding)
    
    self$rnn <- nn_gru(
      input_size = embedding,
      hidden_size = rnn_hidden_size,
      num_layers = rnn_layers,
      bidirectional = TRUE,
      batch_first = TRUE,
      dropout = 0 # war davor 0
    )
    
    self$self_att = torch::nn_multihead_attention(embed_dim = 40L, num_heads = 1L, dropout = 0.1, batch_first = TRUE) # war davor 0.1
    
    # Self-Attention QKV Layer
    self$qkv_size <- rnn_hidden_size * 2  # eigentlich mal 2
    self$q_linear <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = self$qkv_size) # eigentlich mal 2
    self$k_linear <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = self$qkv_size) # eigentlich mal 2
    self$v_linear <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = self$qkv_size) # eigentlich mal 2
    
    self$token_full <- nn_linear(in_features = rnn_hidden_size * 2, out_features = 20L) # davor mal 2
    
    self$trait_full1 <- nn_linear(in_features = 2L, out_features = 100L)
    self$dropout1 = nn_dropout(0.3)
    self$trait_full2 <- nn_linear(in_features = 100L, out_features = 100L)
    self$dropout2 = nn_dropout(0.3)
    self$trait_full3 <- nn_linear(in_features = 100L, out_features = 40L) # 40 is embedding size 
    
    # Combined output layers
    self$output_layer1 <- nn_linear(in_features = 20L, out_features = 50L)
    self$output_layer2 <- nn_linear(in_features = 50L, out_features = 1L)
  },
  # 
  # generate_positional_encoding = function(max_len, embedding_dim) {
  #   pos_enc = torch_zeros(max_len, embedding_dim)
  #   positions = torch_arange(1, max_len)$unsqueeze(2)
  #   div_term = torch_exp(torch_arange(1, embedding_dim, 2) * (-torch_log(torch_tensor(10000)) / embedding_dim))
  #   
  #   pos_enc[, seq(1, embedding_dim, 2)] = torch_sin(positions * div_term)
  #   pos_enc[, seq(2, embedding_dim, 2)] = torch_cos(positions * div_term)
  #   
  #   return(pos_enc$unsqueeze(1))  # [max_len, 1, embedding_dim]
  # },
  
  forward = function(tokens, traits) {
    #### Token Embedding
    x <- self$embedding(tokens)
    
    # Add Positional Encoding to the Token Embeddings
    # seq_len <- tokens$size(2)
    # pos_enc = self$pos_encoding[1:seq_len, , ]$to(dtype = torch_float32(), device = x$device)
    # x <- x + pos_enc$squeeze(1)
    #print(dim(x))
    #### RNN Output
    rnn_out <- self$rnn(x)[[1]]  # [batch_size, seq_len, 2 * rnn_hidden_size]
    #print(dim(rnn_out))
    
    #### Trait Embedding
    traits = traits$unsqueeze(3L)
    traits_position = torch_tensor(matrix(scale(1:10), 10L, 1L),  # skalieren nicht vergessen
                                   device=traits$device,
                                   dtype=traits$dtype)$`repeat`(list(tokens$shape[1], 1L, 1L))
    traits = torch_cat(list(traits, traits_position), dim = 3L)
    
    traits_embedded = self$trait_full1(traits) %>% nnf_relu() %>% self$dropout1()
    traits_embedded = self$trait_full2(traits_embedded) %>% nnf_relu() %>% self$dropout2()
    traits_embedded = self$trait_full3(traits_embedded)
    #print(dim(traits_embedded))
    
    
    #### Concatenation of Token and Trait Embeddings
    #rnn_out_traits = torch_cat(list(rnn_out, traits_embedded), dim = 2L)
    rnn_out_traits = torch_cat(list(rnn_out, traits_embedded), dim = 2L)
    #print(dim(rnn_out_traits))
    
    
    # QKV for Self-Attention
    q <- self$q_linear(rnn_out_traits)  # [batch_size, seq_len, qkv_size]
    k <- self$k_linear(rnn_out_traits)  # [batch_size, seq_len, qkv_size]
    v <- self$v_linear(rnn_out_traits)  # [batch_size, seq_len, qkv_size]
    
    #### Masking for Padding Tokens
    padding_mask <- tokens$clone()
    #print(tokens$shape)
    
    mask <- padding_mask == 1L  # Assuming padding tokens are marked with 1 [batchsize, 177]
    mask = torch_cat(list(mask, torch_zeros(x$shape[1], traits$shape[2], dtype=torch_bool(), device=x$device)), dim = 2)
    
    
    #print(rnn_out$shape)
    #print(traits_embedded$shape)
    #print(rnn_out_traits$shape)
    #print(attn_out$shape)
    
    #### Multi-Head Self-Attention
    self$mask = mask
    attention_out <- self$self_att(q, k, v, key_padding_mask = mask )
    #print(dim(attention_out))
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

#################################################################################
###############################   Input Data   ##################################
#################################################################################

# Chemical and Species Traits
data_matrix <- apply(as.matrix(unique_df[,2:188]), 2, as.numeric) 
#data_matrix_tensor <- torch_tensor(data_matrix, dtype = torch_long()) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
tokens_tensor <- torch_tensor(data_matrix_tokens, dtype = torch_long())
traits_tensor <- torch_tensor(scale(data_matrix[,178:187]), dtype = torch_float32())

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
  model$to(device = device)
  
  train_tokens_fold = train_tokens_tensor$to(device = device)
  train_traits_fold = train_traits_tensor$to(device = device)
  train_LC50_fold = train_LC50_tensor$to(device = device)
  
  test_tokens_fold = test_tokens_tensor$to(device = device)
  test_traits_fold = test_traits_tensor$to(device = device)
  test_LC50_fold = test_LC50_tensor$to(device = device)
  
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
      
      batch_tokens_train = batch_tokens_train$to(device = device)
      batch_traits_train = batch_traits_train$to(device = device)
      batch_LC50_train = batch_LC50_train$to(device = device)
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
        
        batch_tokens_test = batch_tokens_test$to(device = device)
        batch_traits_test = batch_traits_test$to(device = device)
        batch_LC50_test = batch_LC50_test$to(device = device)
        
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
  saveRDS(final_model_weights, "/home/isabellehalbhuber/Toxicology/Scripts/Final_TRAIT_model_weights_noPEV.RDS")
  
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

saveRDS(results_r_squared_Attention_with, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_Attention_with_noPEV.rds")

