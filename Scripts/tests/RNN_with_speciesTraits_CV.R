#################################################################################
#########################  Data Preparatation   #################################
#################################################################################
# load data
library(readr) 
library(dplyr)
Tokens <- readRDS("/home/isabellehalbhuber/Toxicology/Data/Tokens.rds")
df_fish <- readRDS("/home/isabellehalbhuber/Toxicology/Data/df_fish_cv.rds") # df_fish_rf weil da die folds dabei sind 
df_fish <- df_fish %>% select(-QSAR_READY_SMILES, - Fingerprint, -(bit1:bit881), -ALogP, -ALogp2, -AMR) # ohne chemical traits, nur Tokens als input 
token_distanceVectors <- readRDS("/home/isabellehalbhuber/Toxicology/Data/token_distanceVectors.rds")
token_distanceVectors <- as.data.frame(token_distanceVectors)
token_distanceVectors$Index <- rownames(token_distanceVectors)
colnames(token_distanceVectors)[6] <- "token"

# join tokens to initial dataset 
df = Tokens |>
  full_join(df_fish, by = c(SMILES = "SMILES"))

# 1. RNN with species traits 
# exclude species traits 
df_rnn_1 =  df
# only tokens and SMILES
df_rnn_2 <- df_rnn_1 %>% select(-CAS, -scientificNameStd)
# only token 
df_rnn_3 <- df_rnn_2 %>% select(-SMILES, -logmeanConc)


#########################################
### how many unique tokens in dataset with toxicity values without NAs 
all_tokens_noNA <- df_rnn_2[!is.na(df_rnn_2[, 179:196]), ]
all_tokens_noNA <- unlist(all_tokens_noNA[, 2:178])
all_tokens_noNA <- all_tokens_noNA[!is.na(all_tokens_noNA)]
all_unique_tokens <- unique(all_tokens_noNA)
#num_all_unique_tokens <- length(all_unique_tokens)

## lookuptable for token embeddings 
lookup_table <- all_unique_tokens
filtered_df <- token_distanceVectors[token_distanceVectors$token %in% lookup_table, ]

saveRDS(filtered_df, file = "/home/isabellehalbhuber/Toxicology/Data/embeddedTokens.rds")




## jetzt die integer values für alle unique tokens erstellen 
token_to_int <- as.data.frame(all_unique_tokens)
colnames(token_to_int)[1] <- "token"
token_to_int$int <- 1:424
colnames(token_to_int)[2] <- "int"



numbers <- as.character(1:177)
new_colnames <- paste0("position_", numbers)
colnames(df_rnn_2)[2:178] <- new_colnames

na_rows <- df_rnn_2[apply(df_rnn_2, 1, function(row) any(is.na(row))), ]
#na_count <- sum(is.na(df_rnn_2[,179]))
position_df <- df_rnn_2[!is.na(df_rnn_2$logmeanConc), ]
#na_count <- sum(is.na(position_df[,179]))
int_position_df <- position_df
dist_position_df <- position_df
#unique(unlist(position_df[, 2:178]))

######################### Convert Token to distance ########################

# final_df <- dist_position_df
# 
# for (j in 1:(ncol(dist_position_df) - 1)) {  
#   for (v in c("V1", "V2", "V3", "V4", "V5")) {  
#     final_df[[paste0("position_", j, "_", v)]] <- NA
#   }
# }
# 
# 
# for (i in 1:nrow(dist_position_df)) {
#   for (j in 2:ncol(dist_position_df[, 1:178])) {  
#     token <- dist_position_df[i, j]
#     
#    
#     if (!is.na(token) && token %in% token_distanceVectors$token) {
#       vektoren <- token_distanceVectors[token_distanceVectors$token == token, 1:5]
#       
#       
#       for (v in 1:5) {
#         final_df[i, paste0("position_", j-1, "_V", v)] <- vektoren[[v]]
#       }
#     }
#   }
# }
# 


### 1189 Variablen... für jede Position V1-V5 




# for (i in 1:nrow(dist_position_df)) {
#   for (j in 2:ncol(dist_position_df[, 1:178])) {
#     token <- dist_position_df[i, j]
# 
#     if (!is.na(token) && token %in% token_distanceVectors$token) {
#       vektoren <- token_distanceVectors[token_distanceVectors$token == token, 1:5]
#       colnames(vektoren) <- paste0(colnames(vektoren), "_", j-1)
#       dist_position_df <- cbind(dist_position_df, vektoren)
#     }
#   }
# }
# 


# 
#library(reshape2)
# tolower(dist_position_df)
#long_df <- melt(dist_position_df, id.vars = "SMILES", value.name = "token")
# # 
# # 
#df = long_df |>
#full_join(token_distanceVectors, by = c(token = "token"))
# # 
# # df <- df %>% select(-variable, - token)
# # 
# # df = df |>
# #   full_join(unique_df, by = c(SMILES = "SMILES"))

######################### Convert Token to integer ########################

for (i in 1:nrow(position_df)) {
  for (j in 2:ncol(position_df[, 1:178])) {  
    token <- position_df[i, j]
    if (!is.na(token)) {
      int_value <- token_to_int$int[token_to_int$token == token]
      int_position_df[i, j] <- ifelse(length(int_value) > 0, int_value, 0)
    } else {
      int_position_df[i, j] <- 0  
    }
  }
}

#unique(unlist(int_position_df[, 2:178]))

# df without doppelgänger 
unique_df_withSpecies <- int_position_df %>% 
  distinct()
#unique(unlist(unique_df[, 2:178]))

saveRDS(unique_df_withSpecies, file = "/home/isabellehalbhuber/Toxicology/Data/embedding_ready_df_withSpeciesCV.rds")


unique_df_withSpecies_withVectors <- final_df  %>%
  distinct()

saveRDS(unique_df_withSpecies_withVectors, file = "/home/isabellehalbhuber/Toxicology/Data/embedding_ready_df_withSpecies_withVectorsCV.rds")


#################################################################################
######################   RNN with Species Traits   ##############################
#################################################################################

unique_df <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embedding_ready_df_withSpeciesCV.rds")
unique_df <- int_position_df
unique_df <- na.omit(unique_df)
embeddedTokens <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embeddedTokens.rds")


library(torch)


#################################################################################
#################################   RNN    ######################################
#################################################################################

BiRNN <- nn_module(
  initialize = function(embedding = 1L, rnn_hidden_size = 20L, rnn_layers = 3L) {
    self$embedding <- nn_embedding(num_embeddings = 425L, 
                                   embedding_dim = embedding, 
                                   sparse = FALSE, 
                                   scale_grad_by_freq = TRUE) # [batchsize, time, embedding]
    
    # RNN layer with Tokens 
    self$rnn <- nn_gru(input_size = embedding, 
                       hidden_size = rnn_hidden_size, 
                       num_layers = rnn_layers, 
                       bidirectional =TRUE,
                       batch_first = TRUE)
    
    
    # Token NN after RNN
    self$token_full = nn_linear(in_features = 2*rnn_hidden_size, out_features = 20L)
    
    
    # Species Trait NN 
    self$trait_full1 = nn_linear(in_features = 17L, out_features = 20L)
    self$trait_full2 = nn_linear(in_features = 20L, out_features = 3L)
    
    
    # Combined output layers
    self$output_layer1 <- nn_linear(in_features = 20L+3L, out_features = 50L) 
    self$output_layer2 <- nn_linear(in_features = 50L, out_features = 1L) 
  },
  
  forward = function(tokens, traits) {
    #### Token RNN
    x <- self$embedding(tokens)
    x <- self$rnn(x)[[1]] # only the first element, which is the output
    output_token <- self$token_full(x[,dim(x)[2], ]) %>% nnf_relu() %>% torch_squeeze()
    
    #### Trait NN
    x <- self$trait_full1(traits) %>% nnf_relu()
    output_traits <- self$trait_full2(x) %>% nnf_relu() # [batch_size, 3]
    
    ##
    #print(output_token$shape)
    #print(output_traits$shape)
    #### Combine outputs
    input <- torch::torch_cat(list(output_token, output_traits), dim = 2L) # [batchsize, rnn_hidden_size + 3L]
    output <- self$output_layer2(self$output_layer1(input) %>% nnf_relu())
    return(output)
  }
)
model = BiRNN(embedding = 10L)
model$parameters[[1]]
# turn off training of embeddings
model$parameters[[1]]$requires_grad_(FALSE)
# load pretrained embeddings into embedding layer
embd = data.matrix(embeddedTokens)
#embd = matrix(1.0, 425, 10) # mit vortrainierten embeddings von GloVe ersetzen
model$parameters[[1]]$set_data(embd)
model$parameters[[1]]
# train model



#################################################################################
##################################   DNN   ######################################
#################################################################################

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
# 
# model = DNN(embedding = 10L)
# model$parameters[[1]]
# # turn off training of embeddings
# model$parameters[[1]]$requires_grad_(FALSE)
# # load pretrained embeddings into embedding layer
# embd = data.matrix(embeddedTokens)
# #embd = matrix(1.0, 425, 10) # mit vortrainierten embeddings von GloVe ersetzen
# model$parameters[[1]]$set_data(embd)
# model$parameters[[1]]
# # train model
library(torch)

# Definiere das DNN-Modul
DNN <- nn_module(
  initialize = function(embedding = 1L) {
    self$embedding <- nn_embedding(425L, embedding, sparse = FALSE, scale_grad_by_freq = TRUE) # [batchsize, time, embedding]
    self$full1 = nn_linear(in_features = embedding, out_features = 20L)
    self$full2 = nn_linear(in_features = 20L, out_features = 1L)
    self$trait_full1 = nn_linear(in_features = 17L, out_features = 20L)
    self$trait_full2 = nn_linear(in_features = 20L, out_features = 3L)
    
    self$output_layer1 <- nn_linear(in_features = 177L+3L, out_features = 50L) 
    self$output_layer2 <- nn_linear(in_features = 50L, out_features = 1L) 
  },
  
  forward = function(tokens, traits) {
    #### Token NN
    x <- self$embedding(tokens) 
    x = self$full1(x) %>% nnf_selu()
    output_token = self$full2(x) %>% nnf_selu() %>% torch_squeeze()# [batch_size, 177]
    
    #### Trait NN
    x = self$trait_full1(traits) %>% nnf_selu()
    output_traits = self$trait_full2(x) %>% nnf_selu() # [batch_size, 8]
    
    #### Combine outputs
    input = torch::torch_cat(list(output_token, output_traits), dim = 2L) # [batchsize, 177+8L]
    output = self$output_layer2(self$output_layer1(input) %>% nnf_selu())
    return(output)
  }
)


model <- DNN(embedding = 10L)
model$parameters[[1]]$requires_grad_(FALSE)
embd <- data.matrix(embeddedTokens)  
model$embedding$weight <- torch::torch_tensor(embd)
print(model$parameters[[1]])


#################################################################################
###############################   Input Data   ##################################
#################################################################################

# Chemical and Species Traits
data_matrix <- apply(as.matrix(unique_df[,2:195]), 2, as.numeric) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())
traits_tensor <- torch_tensor(scale(data_matrix[,178:194]), dtype = torch_float32())

# LC50
data_matrix_LC50 <- as.numeric(unique_df[,196]) 
LC50_tensor <- torch_tensor(data_matrix_LC50, dtype = torch_float())

# Training and Test datasets
set.seed(123)
# Indices
train_indices <- sample(1:nrow(tokens_tensor), 0.9 * nrow(tokens_tensor))
test_indices <- setdiff(1:nrow(tokens_tensor), train_indices)

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
# 
# model <- DNN(embedding = 10L)
# model$to(device = "cuda:3")
# train_tokens_fold = train_tokens_tensor$to(device = "cuda:3")
# train_traits_fold = train_traits_tensor$to(device = "cuda:3")
# train_LC50_fold = (train_LC50_tensor$to(device = "cuda:3"))
# test_tokens_fold = test_tokens_tensor$to(device = "cuda:3")
# test_traits_fold = test_traits_tensor$to(device = "cuda:3")
# test_LC50_fold = (test_LC50_tensor$to(device = "cuda:3"))


calculate_r2 <- function(y_true, y_pred) {
  res_ss <- sum((y_true - y_pred)^2)
  total_ss <- sum((y_true - mean(y_true))^2)
  r2 <- 1 - res_ss / total_ss
  return(r2)
}

train_and_evaluate <- function(train_indices, test_indices, fold_name) {

  loss_fn <- nn_mse_loss()
  model <- BiRNN(embedding = 10L)
  model$to(device = "cuda:3")
  train_tokens_fold = train_tokens_tensor$to(device = "cuda:3")
  train_traits_fold = train_traits_tensor$to(device = "cuda:3")
  train_LC50_fold = (train_LC50_tensor$to(device = "cuda:3"))
  test_tokens_fold = test_tokens_tensor$to(device = "cuda:3")
  test_traits_fold = test_traits_tensor$to(device = "cuda:3")
  test_LC50_fold = (test_LC50_tensor$to(device = "cuda:3"))
  
  
  optimizer <- optim_adam(model$parameters, lr = 0.001, weight_decay = 0.0001)
  
  num_epochs <- 100L
  train_losses <- c()
  test_losses <- c()
  train_r2_scores <- c()
  test_r2_scores <- c()
  
  for (epoch in 1:num_epochs) {
    model$train()

    sample_size <- min(3000L, length(train_indices))
    ind <- sample(length(train_indices), sample_size)
    
    # Forward propagation
    train_output <- model(train_tokens_fold[ind, , drop = FALSE], train_traits_fold[ind, , drop = FALSE])$relu()
    
    train_loss <- loss_fn(train_output, train_LC50_fold[ind]$unsqueeze(2L))
    train_losses <- c(train_losses, train_loss$item())

    train_r2 <- calculate_r2(as.matrix(train_LC50_fold[ind]), as.matrix(train_output$squeeze()))
    train_r2_scores <- c(train_r2_scores, train_r2)
    
    # Backward propagation
    train_loss$backward()
    optimizer$step()
    optimizer$zero_grad()
    
    # Evaluation
    model$eval()
    with_no_grad({
      test_output <- model(test_tokens_fold, test_traits_fold)$relu()
      
      test_loss <- loss_fn(test_output, test_LC50_fold$unsqueeze(2L))
      test_losses <- c(test_losses, test_loss$item())
      plot(as.matrix(test_output), as.matrix(test_LC50_tensor))
      
      test_r2 <- calculate_r2(as.matrix(test_LC50_fold), as.matrix(test_output$squeeze()))
    })
    test_r2_scores <- c(test_r2_scores, test_r2)
    
    cat(sprintf("Fold: %s, Epoch [%d/%d], Train Loss: %.4f, Test Loss:%.4f\n", 
                fold_name, epoch, num_epochs, train_loss$item(), test_loss$item()))
  }
  
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
  
  fold_results <- train_and_evaluate(train_indices, test_indices, fold_name)
  results[[fold_name]] <- fold_results
}



average_train_loss <- sapply(results, function(res) mean(res$train_losses))
average_test_loss <- sapply(results, function(res) mean(res$test_losses))
average_train_r2 <- sapply(results, function(res) mean(res$train_r2_scores))
average_test_r2 <- sapply(results, function(res) mean(res$test_r2_scores))

cat(sprintf("Average Train Loss per Fold: %s\n", paste(average_train_loss, collapse = ", ")))
cat(sprintf("Average Test Loss per Fold: %s\n", paste(average_test_loss, collapse = ", ")))
cat(sprintf("Average Train R^2 per Fold: %s\n", paste(average_train_r2, collapse = ", ")))
cat(sprintf("Average Test R^2 per Fold: %s\n", paste(average_test_r2, collapse = ", ")))

results_r_squared_with_speciesTraits_DNN <- as.data.frame(cbind(average_train_r2, average_test_r2))results_r_squared_with_speciesTraits_RNN <- as.data.frame(cbind(average_train_r2, average_test_r2))

saveRDS(results_r_squared_with_speciesTraits_DNN, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_with_speciesTraits_DNN.rds") 
saveRDS(results_r_squared_with_speciesTraits_RNN, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_with_speciesTraits_RNN.rds") 


