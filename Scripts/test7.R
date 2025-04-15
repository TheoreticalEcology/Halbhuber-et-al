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

library(torch)

# Dataset Setup (Wird bereits im Setup bereitgestellt)
data_matrix <- apply(as.matrix(df[,1:177]), 2, as.numeric) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, weil Token-Indizes nicht bei 0 beginnen sollen
tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())

# LC50
data_matrix_LC50 <- as.numeric(df[,195])
LC50_tensor <- torch_tensor(data_matrix_LC50, dtype = torch_float())

##### Training und Test-Datensätze ####
set.seed(123)
# Indices
train_indices <- sample(1:nrow(tokens_tensor), 0.8 * nrow(tokens_tensor))
test_indices <- setdiff(1:nrow(tokens_tensor), train_indices)

index_tensor <- function(tensor, indices) {
  # R-Indizes sind 1-basiert, daher werden sie direkt verwendet
  torch_indices <- torch_tensor(indices, dtype = torch_long())
  return(tensor$index_select(1, torch_indices)) # Umwandlung auf 0-basierte Indizes für index_select()
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

# BiGRU Modell
BiGRU_NoTRAITS <- nn_module(
  initialize = function(embedding = 40L, rnn_hidden_size = 20L, rnn_layers = 3L, dropout_prob = 0.3) { 
    self$embedding <- nn_embedding(615L, 
                                   embedding_dim = embedding,
                                   sparse = FALSE,
                                   scale_grad_by_freq = TRUE)
    
    self$rnn <- nn_gru(input_size = embedding, 
                       hidden_size = rnn_hidden_size,
                       num_layers = rnn_layers,  
                       batch_first = TRUE,
                       dropout = ifelse(rnn_layers > 1, dropout_prob, 0), 
                       bidirectional = TRUE)
    
    self$output_layer <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = 1L)
    self$dropout_output <- nn_dropout(p = dropout_prob) 
  },
  
  forward = function(x) {
    x <- self$embedding(x)
    rnn_out <- self$rnn(x)[[1]][, dim(x)[2], ]  
    rnn_out <- self$dropout_output(rnn_out)  
    output <- self$output_layer(rnn_out)
    
    return(output)
  }
)

# Funktion zum Training mit den besten Hyperparametern
train_with_best_hyperparameters <- function(lr, dropout, train_indices, test_indices) {
  loss_fn <- nn_mse_loss()
  model <- BiGRU_NoTRAITS(embedding = 10L, dropout_prob = dropout)
  model$to(device = "cuda:3")
  
  # Train und Test Tensoren für Fold
  train_tokens_fold = index_tensor(tokens_tensor, train_indices)$to(device = "cuda:3")
  train_LC50_fold = LC50_tensor[train_indices]$log()$to(device = "cuda:3")
  test_tokens_fold = index_tensor(tokens_tensor, test_indices)$to(device = "cuda:3")
  test_LC50_fold = LC50_tensor[test_indices]$log()$to(device = "cuda:3")
  
  optimizer <- optim_adam(model$parameters, lr = lr, weight_decay = 0.00001)
  
  train_losses <- c()
  test_losses <- c()
  train_r2_scores <- c()
  test_r2_scores <- c()
  
  # Trainieren auf folds 1-8
  for (epoch in 1:1000) {
    dataloader_train <- create_dataloader(train_tokens_fold, batch_size = 128)
    model$train()
    
    for (batch_indices in dataloader_train) {
      batch_tokens_train <- train_tokens_fold[batch_indices, , drop = FALSE]
      batch_LC50_train <- train_LC50_fold[batch_indices]
      
      batch_tokens_train <- batch_tokens_train$to(device = "cuda:3")
      batch_LC50_train <- batch_LC50_train$to(device = "cuda:3")
      
      # Vorwärtsdurchgang
      train_output <- model(batch_tokens_train)
      train_loss <- loss_fn(train_output, batch_LC50_train$unsqueeze(2L))
      train_losses <- c(train_losses, train_loss$item())
      
      train_r2 <- calculate_r_squared(as.matrix(batch_LC50_train), as.matrix(train_output$squeeze()))
      train_r2_scores <- c(train_r2_scores, train_r2)
      
      # Rückwärtsdurchgang
      train_loss$backward()
      optimizer$step()
      optimizer$zero_grad()
    }
    
    model$eval()
    dataloader_test <- create_dataloader(test_tokens_fold, batch_size = 128)
    
    with_no_grad({
      test_loss_total <- 0
      test_r2_total <- 0
      num_batches <- length(dataloader_test)
      
      for (batch_indices in dataloader_test) {
        batch_tokens_test <- test_tokens_fold[batch_indices, , drop = FALSE]
        batch_LC50_test <- test_LC50_fold[batch_indices]
        
        batch_tokens_test <- batch_tokens_test$to(device = "cuda:3")
        batch_LC50_test <- batch_LC50_test$to(device = "cuda:3")
        
        # Vorwärtsdurchgang
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
    
    cat(sprintf("Epoch [%d/1000], Train Loss: %.4f, Test R²: %.4f\n", epoch, train_loss$item(), avg_test_r2))
  }
  
  return(list(train_losses = train_losses, test_losses = test_losses,
              train_r2_scores = train_r2_scores, test_r2_scores = test_r2_scores, model = torch::torch_serialize(model)))
}

# Funktion zum Hyperparameter-Tuning
hyperparameter_tuning <- function(folds, hyperparameters) {
  best_hyperparameters <- NULL
  best_test_r2 <- -Inf  # Setze als schlechtesten Anfangswert für R²
  
  # Gehe alle Hyperparameter-Kombinationen durch
  for (params in 1:nrow(hyperparameters)) {
    lr <- hyperparameters[params, "lr"]
    dropout <- hyperparameters[params, "dropout"]
    cat("Training mit lr =", lr, "und dropout =", dropout, "\n")
    
    # Trainiere und teste auf folds 1-8
    fold_r2_scores <- c()
    for (fold_id in 1:8) {
      # Trainiere und teste mit den aktuellen Hyperparametern
      fold_results <- train_with_best_hyperparameters(lr, dropout, train_indices, test_indices)
      fold_r2_scores <- c(fold_r2_scores, fold_results$test_r2_scores)
    }
    
    # Berechne den durchschnittlichen R²-Wert für diese Hyperparameter-Kombination
    avg_test_r2 <- mean(fold_r2_scores)
    cat("Durchschnittlicher Test R² für diese Kombination: ", avg_test_r2, "\n")
    
    # Speichere die besten Hyperparameter
    if (avg_test_r2 > best_test_r2) {
      best_test_r2 <- avg_test_r2
      best_hyperparameters <- list(lr = lr, dropout = dropout)
    }
  }
  
  return(best_hyperparameters)
}

# Beispielaufruf für Hyperparameter-Tuning und Training:
best_hyperparameters <- hyperparameter_tuning(folds = 1:8, hyperparameters = hyper_grid_object)

# Verwende die besten Hyperparameter für Training
best_lr <- best_hyperparameters$lr
best_dropout <- best_hyperparameters$dropout

# Training mit den besten Hyperparametern
train_results <- train_with_best_hyperparameters(best_lr, best_dropout, train_indices, test_indices)

# Evaluierung auf Fold 9 und Fold 10
fold_9_indices <- sample(1:nrow(tokens_tensor), 0.1 * nrow(tokens_tensor)) # Beispiel für Fold 9
fold_10_indices <- setdiff(1:nrow(tokens_tensor), fold_9_indices)

final_results <- evaluate_on_folds_9_and_10(train_results$model, fold_9_indices, fold_10_indices)

cat("Finale Test R² für Fold 9: ", final_results$test_r2_9, "\n")
cat("Finale Test R² für Fold 10: ", final_results$test_r2_10, "\n")

