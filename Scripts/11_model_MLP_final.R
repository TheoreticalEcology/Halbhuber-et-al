############################## MLP ###############################
# Set up four MLP models: traits & Mol. Descriptors /traits & Fingerprints/ no traits & Mol. Descriptors / no traits & Fingerprints
# Each model undergoes blocked validation: random, species and chemicals blocked  
library(readr)
library(dplyr)
library(cito)
source("~/paper/EcoToxTM/Scripts/functions.R") # for R squared calculation
allData <- readRDS("~/paper/EcoToxTM/Data/df_final&converted.rds")
# Folds from factor to numeric
allData[, 225:227] <- sapply(allData[, 225:227], function(x) as.numeric(as.character(x)))


# scale variables
variables_to_scale <- c(colnames(allData[,1:223])) # scale MDs and species traits 
scale <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}
allData  <- allData %>%
  mutate(across(all_of(variables_to_scale), scale))
# some columns show no variance in their values and therefore lead to NaN after scaling. 
allData  <- allData [, colSums(is.na(allData )) == 0]

# prepare data for the loop 
dat = list(
  df_FP <- allData[, -c(1:158)],
  df_MD <- allData[, -c(180:345)]  
)

#### Hyperparameter Grid ####
hyper_grid <- function(lr, dropout) {
  expand.grid(
    lr = lr,
    dropout = dropout,
    train = NA,
    test = NA,
    val = NA,
    stringsAsFactors = FALSE
  )
}

fold_columns <- c("randomFold", "chemicalFold", "speciesFold")
descriptors <- c("Fingerprints", "MolDescriptor")
trait_levels <- c("noTraits", "Traits")

hyper_grid_object <- hyper_grid(lr = c(0.001, 0.005, 0.007, 0.01), dropout = c(0.1, 0.5))

hyperparameters <- list(MLP = list())

for (trait in trait_levels) {
  hyperparameters$MLP[[trait]] <- list()
  
  for (desc in descriptors) {
    hyperparameters$MLP[[trait]][[desc]] <- list()
    
    # Loop over all folds
    for (fold_column in fold_columns) {
      hyperparameters$MLP[[trait]][[desc]][[fold_column]] <- hyper_grid_object
    }
  }
}


#### Model-Loop ####

for (desc in descriptors) {
  df <- if (desc == "MolDescriptor") df_MD else df_FP
  
  for (fold_column in fold_columns) {
    for (trait in trait_levels) {
      
      k <- 10
      hp_list <- hyperparameters$MLP[[trait]][[desc]][[fold_column]]
      nhyper <- nrow(hp_list)
      
      trait_cols_to_remove <- c("BEl", "VEp", "REs", "OGp", "BLs", "PFs", "PFv", "CPt", "RMl",
                                "mem_1", "mem_2", "mem_3", "mem_4", "mem_5", "mem_6", "mem_7", "body_size")
      base_cols_to_remove <- c("randomFold", "chemicalFold", "speciesFold")
      cols_to_remove <- if (trait == "Traits") base_cols_to_remove else c(base_cols_to_remove, trait_cols_to_remove)
      
      # Phase 1: Training/Test on folds 1:(k-2)
      for (j in 1:nhyper) {
        lr_val <- hp_list$lr[j]
        dropout_val <- hp_list$dropout[j]
        
        R_squared_train_values <- rep(NA, (k - 2))
        R_squared_test_values  <- rep(NA, (k - 2))
        
        for (i in 1:(k - 2)) {
          train_data_full <- df[df[[fold_column]] != i & df[[fold_column]] <= (k - 2), ]
          test_data_full  <- df[df[[fold_column]] == i, ]
          
          train_data <- train_data_full %>% select(-all_of(cols_to_remove))
          test_data  <- test_data_full %>% select(-all_of(cols_to_remove))
          
          mlp_model <- dnn(
            logmeanConc ~ .,
            hidden = c(50L, 50L),
            batchsize = 200L,
            data = train_data,
            early_stopping = 10L,
            loss = "mae",
            burnin = Inf,
            activation = "relu",
            lr = lr_val,
            dropout = dropout_val,
            device = "cuda",
            validation = 0.2,
            epochs = 500L,
            optimizer = config_optimizer("sgd", weight_decay = 0.001)
          )
          
          train_pred <- predict(mlp_model, newdata = train_data)
          test_pred  <- predict(mlp_model, newdata = test_data)
          
          R_squared_train_values[i] <- calculate_r_squared(train_data$logmeanConc, train_pred)
          R_squared_test_values[i]  <- calculate_r_squared(test_data$logmeanConc, test_pred)
        }
        
        # save
        hyperparameters$MLP[[trait]][[desc]][[fold_column]]$train[j] <- mean(R_squared_train_values, na.rm = TRUE)
        hyperparameters$MLP[[trait]][[desc]][[fold_column]]$test[j]  <- mean(R_squared_test_values, na.rm = TRUE)
      }
      
      # Phase 2: Select best hyperparameter-set (highest test-R²)
      best_index <- which.max(hyperparameters$MLP[[trait]][[desc]][[fold_column]]$test)
      best_lr <- hp_list$lr[best_index]
      best_dropout <- hp_list$dropout[best_index]
 
      
      # Train on folds 1:(k-2), evaluate on folds k-1 und k
      val_train <- df[df[[fold_column]] <= (k - 2), ]
      val_data  <- df[df[[fold_column]] %in% c(k - 1, k), ]
      
      val_train <- val_train %>% select(-all_of(cols_to_remove))
      val_data  <- val_data %>% select(-all_of(cols_to_remove))
      
      mlp_model_val <- dnn(
        logmeanConc ~ .,
        hidden = c(50L, 50L),
        batchsize = 200L,
        data = val_train,
        early_stopping = 10L,
        loss = "mae",
        burnin = Inf,
        activation = "relu",
        lr = best_lr,
        dropout = best_dropout,
        device = "cuda",
        validation = 0.2,
        epochs = 500L
      )
      
      val_pred <- predict(mlp_model_val, newdata = val_data)
      R_squared_val <- calculate_r_squared(val_data$logmeanConc, val_pred)
      
      # Save final Val-values with best hyperparameters
      hyperparameters$MLP[[trait]][[desc]][[fold_column]]$val <- rep(NA, nhyper)
      hyperparameters$MLP[[trait]][[desc]][[fold_column]]$val[best_index] <- R_squared_val
      
      #cat("Final eval for", trait, "|", desc, "|", fold_column, "| Best combination: lr =", best_lr, ", dropout =", best_dropout, ", lambda =", best_lambda, "\n")
    }
  }
}





#### Save results ####
saveRDS(hyperparameters$MLP$noTraits$MolDescriptor,file = "~/paper/EcoToxTM/Results/MLP_noTraits_MolDescriptor.rds")
saveRDS(hyperparameters$MLP$noTraits$Fingerprints,file = "~/paper/EcoToxTM/Results/MLP_noTraits_Fingerprints.rds")
saveRDS(hyperparameters$MLP$Traits$MolDescriptor,file = "~/paper/EcoToxTM/Results/MLP_Traits_MolDescriptor.rds")
saveRDS(hyperparameters$MLP$Traits$Fingerprints,file = "~/paper/EcoToxTM/Results/MLP_Traits_Fingerprints.rds")


