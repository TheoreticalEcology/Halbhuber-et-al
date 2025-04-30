##############################   Transformer TRAITS   ######################################
library(readr) 
library(dplyr)
library(torch)
source("~/paper/EcoToxTM/Scripts/functions.R") 
source("~/paper/EcoToxTM/Scripts/transformer-utils.R") 
df <- readRDS("~/paper/EcoToxTM/Data/df_final_tokenizedSMILES.rds")

#### Data preparation ####
df <- na.omit(df)
df <- df[,-1] # exclude SMILES codes 
df[, 22:198] <- lapply(df[, 22:198], as.numeric)

# flip the padded sequence so that the important information is at the back (then backpropagation is better)
remaining_columns_traits <- df[c(1:21)]
reversed_columns <- df[, 22:198][, ncol(df[, 22:198]):1]
df <- cbind(reversed_columns, remaining_columns_traits)
# check if there are rows that contain no tokens --> exclude... 
df_num <- df
df_num <- lapply(df_num[, 1:177], as.numeric)
rows_all_zeros <- (df_num[[177]] == 0)
which(rows_all_zeros)

df <- df[!rows_all_zeros, ]

# Chemical and Species Traits
data_matrix <- apply(as.matrix(df[,1:177]), 2, as.numeric) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())
data_matrix_traits <- apply(as.matrix(df[,178:194]), 2, as.numeric)  
traits_tensor <- torch_tensor(data_matrix_traits, torch_float32())
# LC50
data_matrix_LC50 <- as.numeric(df[,195])
LC50_tensor <- torch_tensor(data_matrix_LC50, dtype = torch_float())

##### Helper functions ####

index_tensor <- function(tensor, indices) {
  # R-Indizes sind 1-basiert, also direkt verwenden
  torch_indices <- torch_tensor(indices, dtype = torch_long())
  tensor$index_select(1, torch_indices) # 1. Dimension --> Rows!Und Umwandlung auf 0-basierte Indizes für index_select()
}

create_dataloader <- function(data_tensor, batch_size) {
  num_samples <- data_tensor$size(1) 
  indices <- sample(1:num_samples) 
  batches <- split(indices, ceiling(seq_along(indices) / batch_size))
  return(batches) 
}

#### Model configurations ####
transformer = nn_module(
  "transformer",
  initialize = function(max_token = 615L, 
                        traits = NULL,
                        dropout = 0.3,
                        emb_dim = 16L,
                        num_layers = 1L,
                        num_heads = 1L,
                        dim_feedforward = 512L) {
    
    self$emb_dim = emb_dim
    self$token_embedder = nn_embedding(max_token, emb_dim)
    
    self$traits = FALSE
    if(!is.null(traits)) {
      self$traits = TRUE
      self$traits_embedder = nn_sequential(nn_linear(2L, 50L), 
                                           nn_gelu(), 
                                           nn_linear(50L, emb_dim))
    }
    self$positional_encoding = PositionalEncoding(emb_dim, dropout, max_len = 1000L)
    encoder_layer = TransformerEncoderLayer(d_model=emb_dim, nhead=num_heads,batch_first = TRUE, dim_feedforward = dim_feedforward)
    self$transformer_encoder = TransformerEncoder(encoder_layer, num_layers=num_layers)
    self$output = nn_sequential(nn_linear(emb_dim, 100L),
                                nn_gelu(),
                                nn_dropout(dropout),
                                nn_linear(100, 1L))
  }, 
  forward = function(tokens, traits = NULL) {
    embedded_token = self$token_embedder(tokens)
    if(!is.null(traits)) {
      traits_position = torch_tensor(matrix(scale(1:traits$shape[2]), traits$shape[2], 1L),  # skalieren nicht vergessen
                                     device=traits$device,
                                     dtype=traits$dtype)$`repeat`(list(tokens$shape[1], 1L, 1L))
      traits_pos = torch_cat(list(traits$unsqueeze(3L), traits_position), dim = 3L)
      embedded_traits = self$traits_embedder(traits_pos)
      embedded_token = torch_cat(list(embedded_token, embedded_traits), dim = 2L)
    }
    embedded_token = self$positional_encoding(embedded_token)
    padding_mask <- tokens$clone()
    mask <- padding_mask == 1L  # Assuming padding tokens are marked with 1 [batchsize, 177]
    mask = torch_cat(list(mask, torch_zeros(tokens$shape[1], traits$shape[2], dtype=torch_bool(), device=tokens$device)), dim = 2)
    transformer_output = self$transformer_encoder(embedded_token, src_key_padding_mask = mask)
    return(self$output(transformer_output$mean(2L)))
  }
  
)


#### Hyperparameter Grid Setup ####

hyper_grid <- function(lr, weight_decay) {
  expand.grid(
    lr = lr,
    weight_decay = weight_decay,
    train = NA,
    test = NA,
    val_train = NA,
    val_test = NA,
    stringsAsFactors = FALSE
  )
}

fold_columns <- c("randomFold", "chemicalFold", "speciesFold")

#hyper_grid_object <- hyper_grid(lr = c(0.001, 0.005, 0.007, 0.01), weight_decay = c(0.0001, 0.00001, 0.000001))
hyper_grid_object <- hyper_grid(lr = c(0.001, 0.01), weight_decay = c(0.0001, 0.001))


hyperparameters <- list(Transformer = list())

for (fold_column in fold_columns) {
  hyperparameters$Transformer[[fold_column]] <- hyper_grid_object
}

##
device = "cuda:3"
k = 10
torch::cuda_is_available()

# save data
average_r2_summary <- data.frame(fold_column = character(), 
                                 #best_hp_mean_r2 = numeric(), 
                                 r2_test_val = numeric(),
                                 r2_train_val = numeric(),
                                 best_lr = numeric(),
                                 best_wd = numeric(),
                                 stringsAsFactors = FALSE)

for (fold_column in fold_columns) {
  cat("Processing:", fold_column, "\n")
  fold_column <- "chemicalFold"
  hp_list <- hyperparameters$Transformer[[fold_column]]
  nhyper <- nrow(hp_list)
  nepochs <- 1000
  
  R_squared_test_values <- matrix(NA, nrow = nhyper, ncol = k)
  R_squared_train_values <- matrix(NA, nrow = nhyper, ncol = k)
  ALL_R_squared_train_values <- array(NA, dim = c(nhyper, nepochs, k))
  ALL_R_squared_test_values <- array(NA, dim = c(nhyper, nepochs, k))
  R_squared_train_val_values <- matrix(NA, nrow = nhyper, ncol = k)
  R_squared_test_val_values <- matrix(NA, nrow = nhyper, ncol = k)
  ALL_R_squared_train_val_values <- array(NA, dim = c(nhyper, nepochs, k))
  ALL_R_squared_test_val_values <- array(NA, dim = c(nhyper, nepochs, k))
  
  ##### PHASE 1: Hyperparameter-Tuning (Fold 1–8) #####
  for (j in 1:nhyper) {
    lr <- hp_list$lr[j]
    weight_decay <- hp_list$weight_decay[j]
    
    for (i in 1:(k - 2)) {
      
      tokens_tensor <- tokens_tensor$to(device = device)
      traits_tensor <- traits_tensor$to(device = device)
      LC50_tensor <- LC50_tensor$to(device = device)
      
      train_indices <- torch_tensor(which(df[[fold_column]] != i & df[[fold_column]] <= (k - 2)),
                                    dtype = torch_long(), device = device)
      test_indices <- torch_tensor(which(df[[fold_column]] == i),
                                   dtype = torch_float(), device = device)
      
      train_tokens <- index_tensor(tokens_tensor, train_indices)
      test_tokens  <- index_tensor(tokens_tensor, test_indices)
      train_traits <- index_tensor(traits_tensor, train_indices)
      test_traits  <- index_tensor(traits_tensor, test_indices)
      train_LC50 <- index_tensor(LC50_tensor, train_indices)$log()
      test_LC50  <- index_tensor(LC50_tensor, test_indices)$log()
      
      #### CONFIG ####
      
      model <- transformer(traits = 17L, num_layers = 5L, num_heads = 3L, emb_dim = 36L)
      model <- model$to(device = device)
      optimizer <- optim_adam(model$parameters, lr = lr, weight_decay = weight_decay)
      loss_fn <- nn_mse_loss()
      
      r2_test_best <- -Inf
      r2_train_best <- -Inf
      
      
      for (epoch in 1:nepochs) {
        cat(sprintf("Epoch: %d | Fold column: %s | Test-Fold: %d | Learning rate: %.6f | Weight decay: %.6f\n",
                    epoch, fold_column, i, hp_list$lr[j], hp_list$weight_decay[j]))
        
        #### TRAIN ####
        model$train()
        
        dataloader_train <- create_dataloader(train_tokens, 150)
        
        for (l in seq_along(dataloader_train)) {
          batch_tokens <- (train_tokens[dataloader_train[[l]], , drop = FALSE])$to(device = device)
          batch_traits <- (train_traits[dataloader_train[[l]], , drop = FALSE])$to(device = device)
          batch_LC50 <- (train_LC50[dataloader_train[[l]]])$to(device = device)
          output <- model(batch_tokens, batch_traits)
          r2_train <- calculate_r_squared(as.matrix(batch_LC50), as.matrix(output$squeeze()))
          ALL_R_squared_train_values[j, epoch, i] <- r2_train
          if (!is.na(r2_train) && r2_train > r2_train_best) {
            r2_train_best <- r2_train
          }
          
          loss <- loss_fn(output, batch_LC50$unsqueeze(2L))
          loss$backward()
          optimizer$step()
          optimizer$zero_grad()
        }
        
        
        #### TEST ####
        model$eval()
        
        with_no_grad({
          dataloader_test <- create_dataloader(test_tokens, batch_size = 150)
          for (m in seq_along(dataloader_test)) {
            batch_tokens <- (test_tokens[dataloader_test[[m]], , drop = FALSE])$to(device = device)
            batch_traits <- (test_traits[dataloader_test[[m]], , drop = FALSE])$to(device = device)
            batch_LC50 <- (test_LC50[dataloader_test[[m]]])$to(device = device)
            output <- model(batch_tokens, batch_traits)
            r2_test <- calculate_r_squared(as.matrix(batch_LC50), as.matrix(output$squeeze()))
            ALL_R_squared_test_values[j, epoch, i] <- r2_test 
            if (!is.na(r2_test) && r2_test > r2_test_best) {
              r2_test_best <- r2_test
              
            }
          }
          cat(sprintf("Epoch %d: %.6f\n", epoch, r2_test))
        })
      }
      
      # Max R² over epochs per fold # reminder: j is the length of hp_list and i is the number of folds
      R_squared_train_values[j, i] <- r2_train_best
      R_squared_test_values[j, i] <- r2_test_best
      
    }
  }
  #ALL_R_squared_test_values[j*nepochs, i] <- r2_test 
  #ALL_R_squared_train_values[j*nepochs, i] <- r2_train
  
  # mean over folds (not Epochs!) of test R2
  mean_r2_per_hp_test <- rowMeans(R_squared_test_values, na.rm = TRUE)
  mean_r2_per_hp_train <- rowMeans(R_squared_train_values, na.rm = TRUE)
  best_index <- which.max(mean_r2_per_hp_test)
  best_lr <- hp_list$lr[best_index]
  best_wd <- hp_list$weight_decay[best_index]
  
  ##### PHASE 2: Final Training (Fold 1–8) Evaluaiton  (Fold 9 & 10) #####
  
  tokens_tensor <- tokens_tensor$to(device = device)
  traits_tensor <- traits_tensor$to(device = device)
  LC50_tensor <- LC50_tensor$to(device = device)
  
  
  val_indices <- torch_tensor(which(df[[fold_column]] <= (k - 2)),
                              dtype = torch_long(), device = device)
  val_test_indices <- torch_tensor(which(df[[fold_column]] %in% c(k - 1, k)),
                                   dtype = torch_long(), device = device)
  
  
  val_train_tokens <- index_tensor(tokens_tensor, val_indices)
  val_test_tokens  <- index_tensor(tokens_tensor, val_test_indices)
  val_train_traits <- index_tensor(traits_tensor, val_indices)
  val_test_traits  <- index_tensor(traits_tensor, val_test_indices)
  val_train_LC50 <- index_tensor(LC50_tensor, val_indices)$log()
  val_test_LC50  <- index_tensor(LC50_tensor, val_test_indices)$log()
  
  # train best model
  model_val <- transformer(traits = 17L, num_layers = 5L, num_heads = 3L, emb_dim = 36L)
  model_val <- model_val$to(device = device)
  optimizer_val <- optim_adam(model_val$parameters, lr = best_lr, weight_decay = best_wd)
  loss_fn <- nn_mse_loss()
  
  val_r2_best <- -Inf
  
  for (epoch in 1:nepochs) {
    cat(sprintf("LAST Epoch: %d | Fold column: %s | FINAL Training | Learning rate: %.6f | Weight decay: %.6f\n",
                epoch, fold_column, best_lr, best_wd))
    
    # Train
    model_val$train()
    dataloader_val_train <- create_dataloader(val_train_tokens, batch_size = 150)
    
    for (n in seq_along(dataloader_val_train)) {
      batch_tokens <- val_train_tokens[dataloader_val_train[[n]], , drop = FALSE]$to(device = device)
      batch_traits <- val_train_traits[dataloader_val_train[[n]], , drop = FALSE]$to(device = device)
      batch_LC50   <- val_train_LC50[dataloader_val_train[[n]] ]$to(device = device)
      
      output <- model_val(batch_tokens, batch_traits)
      loss <- loss_fn(output, batch_LC50$unsqueeze(2L))
      loss$backward()
      optimizer_val$step()
      optimizer_val$zero_grad()
      
      r2_train_val <- calculate_r_squared(as.matrix(batch_LC50), as.matrix(output$squeeze()))
      ALL_R_squared_train_val_values[best_index, epoch, 1] <- r2_train_val
    }
    
    # Test
    model_val$eval()
    
    with_no_grad({
      dataloader_val_test <- create_dataloader(val_test_tokens, batch_size = 150)
      for (h in seq_along(dataloader_val_test)) {
        batch_tokens <- val_test_tokens[dataloader_val_test[[h]], , drop = FALSE]$to(device = device)
        batch_traits <- val_test_traits[dataloader_val_test[[h]], , drop = FALSE]$to(device = device)
        batch_LC50   <- val_test_LC50[dataloader_val_test[[h]] ]$to(device = device)
        
        output <- model_val(batch_tokens, batch_traits)
        
        r2_test_val <- calculate_r_squared(as.matrix(batch_LC50), as.matrix(output$squeeze()))
        ALL_R_squared_test_val_values[best_index, epoch, 1] <- r2_test_val
      }
      cat(sprintf("Epoch %d: %.6f\n", epoch, r2_test_val))
    })
  }
  
  
  hyperparameters$Transformer[[fold_column]]$train <- mean_r2_per_hp_train
  hyperparameters$Transformer[[fold_column]]$test  <- mean_r2_per_hp_test
  hyperparameters$Transformer[[fold_column]]$val_train[best_index] <- r2_train_val
  hyperparameters$Transformer[[fold_column]]$val_test[best_index]  <- r2_test_val
  
}

saveRDS(hyperparameters$Transformer,file = "~/paper/EcoToxTM/Results/Transformer_Traits_SMILES.rds")
