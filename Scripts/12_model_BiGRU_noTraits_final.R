##############################   Bi-GRU NO TRAITS   ######################################
library(readr) 
library(dplyr)
library(torch)
source("~/paper/EcoToxTM/Scripts/functions.R") 
source("~/paper/EcoToxTM/Scripts/BiGRU_models.R") 
df <- readRDS("~/paper/EcoToxTM/Data/df_final_tokenizedSMILES.rds")

#### Data preparation ####
df <- na.omit(df)
df <- df[,-1] # exclude SMILES codes 
df[, 22:198] <- lapply(df[, 22:198], as.numeric)

# flip the padded sequence so that the important information is at the back (then backpropagation is better)
remaining_columns_traits <- df[c(1:21)]
reversed_columns <- df[, 22:198][, ncol(df[, 22:198]):1]
df <- cbind(reversed_columns, remaining_columns_traits)

# Chemical and Species Traits
data_matrix <- apply(as.matrix(df[,1:177]), 2, as.numeric) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())

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
hyper_grid_object <- hyper_grid(lr = c(0.001, 0.01), weight_decay = c(0.000001))


hyperparameters <- list(BiGRU = list())

for (fold_column in fold_columns) {
  hyperparameters$BiGRU[[fold_column]] <- hyper_grid_object
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
  hp_list <- hyperparameters$BiGRU[[fold_column]]
  nhyper <- nrow(hp_list)
  nepochs <- 3
  
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
      LC50_tensor <- LC50_tensor$to(device = device)
      
      train_indices <- torch_tensor(which(df[[fold_column]] != i & df[[fold_column]] <= (k - 2)),
                                    dtype = torch_long(), device = device)
      test_indices <- torch_tensor(which(df[[fold_column]] == i),
                                   dtype = torch_long(), device = device)
      
      train_tokens <- index_tensor(tokens_tensor, train_indices)
      test_tokens  <- index_tensor(tokens_tensor, test_indices)
      train_LC50 <- index_tensor(LC50_tensor, train_indices)$log()
      test_LC50  <- index_tensor(LC50_tensor, test_indices)$log()
      
      #### CONFIG ####
      
      model <- BiGRU_NoTRAITS()
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

        dataloader_train <- create_dataloader(train_tokens, 128)
        
        for (l in seq_along(dataloader_train)) {
          batch_tokens <- (train_tokens[dataloader_train[[l]], , drop = FALSE])$to(device = device)
          batch_LC50 <- (train_LC50[dataloader_train[[l]]])$to(device = device)
          output <- model(batch_tokens)
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
          dataloader_test <- create_dataloader(test_tokens, batch_size = 128)
          for (m in seq_along(dataloader_test)) {
            batch_tokens <- (test_tokens[dataloader_test[[m]], , drop = FALSE])$to(device = device)
            batch_LC50 <- (test_LC50[dataloader_test[[m]]])$to(device = device)
            output <- model(batch_tokens)
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
  LC50_tensor <- LC50_tensor$to(device = device)
  
  val_indices <- torch_tensor(which(df[[fold_column]] <= (k - 2)),
                              dtype = torch_long(), device = device)
  val_test_indices <- torch_tensor(which(df[[fold_column]] %in% c(k - 1, k)),
                                   dtype = torch_long(), device = device)
  
  val_train_tokens <- index_tensor(tokens_tensor, val_indices)
  val_test_tokens  <- index_tensor(tokens_tensor, val_test_indices)
  val_train_LC50 <- index_tensor(LC50_tensor, val_indices)$log()
  val_test_LC50  <- index_tensor(LC50_tensor, val_test_indices)$log()
  
  # train best model
  model_val <- BiGRU_NoTRAITS()
  model_val <- model_val$to(device = device)
  optimizer_val <- optim_adam(model_val$parameters, lr = best_lr, weight_decay = best_wd)
  loss_fn <- nn_mse_loss()
  
  val_r2_best <- -Inf
  
  for (epoch in 1:nepochs) {
    cat(sprintf("LAST Epoch: %d | Fold column: %s | FINAL Training | Learning rate: %.6f | Weight decay: %.6f\n",
                epoch, fold_column, best_lr, best_wd))
    
    # Train
    model_val$train()
    dataloader_val_train <- create_dataloader(val_train_tokens, batch_size = 128)
    
    for (n in seq_along(dataloader_val_train)) {
      batch_tokens <- val_train_tokens[dataloader_val_train[[n]], , drop = FALSE]$to(device = device)
      batch_LC50   <- val_train_LC50[dataloader_val_train[[n]] ]$to(device = device)
      output <- model_val(batch_tokens)
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
      dataloader_val_test <- create_dataloader(val_test_tokens, batch_size = 128)
      for (h in seq_along(dataloader_val_test)) {
        batch_tokens <- val_test_tokens[dataloader_val_test[[h]], , drop = FALSE]$to(device = device)
        batch_LC50   <- val_test_LC50[dataloader_val_test[[h]] ]$to(device = device)
        output <- model_val(batch_tokens)
        r2_test_val <- calculate_r_squared(as.matrix(batch_LC50), as.matrix(output$squeeze()))
        ALL_R_squared_test_val_values[best_index, epoch, 1] <- r2_test_val
      }
      cat(sprintf("Epoch %d: %.6f\n", epoch, r2_test_val))
    })
  }
  

  hyperparameters$BiGRU[[fold_column]]$train <- mean_r2_per_hp_train
  hyperparameters$BiGRU[[fold_column]]$test  <- mean_r2_per_hp_test
  hyperparameters$BiGRU[[fold_column]]$val_train[best_index] <- r2_train_val
  hyperparameters$BiGRU[[fold_column]]$val_test[best_index]  <- r2_test_val
   
}

#saveRDS(average_r2_summary, file = "~/paper/EcoToxTM/Results/BiGRU_noTraits_r2_summary.rds")
#saveRDS(ALL_R_squared_test_val_values, file = "~/paper/EcoToxTM/Results/ALL_R_squared_test_val_values.rds")
#saveRDS(ALL_R_squared_train_val_values, file = "~/paper/EcoToxTM/Results/ALL_R_squared_train_val_values.rds")
#saveRDS(ALL_R_squared_test_values, file = "~/paper/EcoToxTM/Results/ALL_R_squared_test_values.rds")
#saveRDS(ALL_R_squared_train_values, file = "~/paper/EcoToxTM/Results/ALL_R_squared_train_values.rds")

saveRDS(hyperparameters$BiGRU,file = "~/paper/EcoToxTM/Results/BiGRU_noTraits_SMILES.rds")


