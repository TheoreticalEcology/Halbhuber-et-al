############################## MLP ###############################
# Set up four MLP models: traits & Mol. Descriptors /traits & Fingerprints/ no traits & Mol. Descriptors / no traits & Fingerprints
# Each model undergoes blocked validation: random, species and chemicals blocked  
library(readr)
library(dplyr)
library(cito)
source("~/paper/EcoToxTM/Scripts/functions.R") # for R squared calculation
allData <- readRDS("~/paper/EcoToxTM/Data/df_final&converted.rds")

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
hyper_grid <- function(hidden, lr, batchsize) {
  expand.grid(
    hidden = hidden,
    lr = lr,
    batchsize = batchsize,
    train = NA,
    test = NA,
    val = NA,
    best_epoch = NA,
    stringsAsFactors = FALSE
  )
}

folds <- c("randomFold", "chemicalFold", "speciesFold")
descriptors <- c("Fingerprints", "MolDescriptor")
trait_levels <- c("noTraits", "Traits")



# hyper_grid_object <- hyper_grid(hidden = c(c(64, 64),c(128, 64),c(128, 128, 64),c(100L, 150L)), lr = c(0.01, 0.001, 0.0005), batchsize = c(32, 50, 64, 128))

# hyper_grid_object <- hyper_grid(hidden = c(64, 100, 128, 150),lr = c(0.01, 0.001, 0.0005), batchsize = c(32, 50, 64, 128))

hyper_grid_object <- hyper_grid(hidden = c(64),lr = c(0.01), batchsize = c(32))

hyperparameters <- list(MLP = list())

for (trait in trait_levels) {
  hyperparameters$MLP[[trait]] <- list()
  
  for (desc in descriptors) {
    hyperparameters$MLP[[trait]][[desc]] <- list()
    
    # Loop through all folds
    for (fold in folds) {
      hyperparameters$MLP[[trait]][[desc]][[fold]] <- hyper_grid_object
    }
  }
}


#### Model-Loop ####
for (desc in descriptors) {
  df <- if (desc == "MolDescriptor") df_MD else df_FP
  
  for (fold_column in folds) {
    for (trait in trait_levels) {
      
      k <- 10
      hp_list <- hyperparameters$MLP[[trait]][[desc]][[fold_column]]
      nhyper <- nrow(hp_list)
      
      trait_cols_to_remove <- c("BEl", "VEp", "REs", "OGp", "BLs", "PFs", "PFv", "CPt", "RMl",
                                "mem_1", "mem_2", "mem_3", "mem_4", "mem_5", "mem_6", "mem_7", "body_size")
      base_cols_to_remove <- c("randomFold", "chemicalFold", "speciesFold")
      cols_to_remove <- if (trait == "Traits") base_cols_to_remove else c(base_cols_to_remove, trait_cols_to_remove)
      
      for (j in 1:nhyper) {
        hidden_val <- hp_list$hidden[j]
        lr_val <- hp_list$lr[j]
        batchsize_val <- hp_list$batchsize[j]
        
        R_squared_train_values <- rep(NA, (k - 2))
        R_squared_test_values  <- rep(NA, (k - 2))
        
        for (i in 1:(k - 2)) {
          train_data_full <- df[df[[fold_column]] != i, ]
          test_data_full  <- df[df[[fold_column]] == i, ]
          
          train_data <- train_data_full %>% select(-all_of(cols_to_remove))
          test_data  <- test_data_full %>% select(-all_of(cols_to_remove))
          
          best_epoch <- NA
          
          for (ep in 1:100) {
            mlp_model <- dnn(
              logmeanConc ~ .,
              hidden = hidden_val,
              batchsize = batchsize_val,
              data = train_data,
              epochs = 100,  
              loss = "mae",
              burnin = 5,
              activation = "relu",
              lr = lr_val,
              device = "cuda"
            )
            
            best_epoch <- ep
            
            # Prediction for test and train data
            train_predictions <- predict(mlp_model, newdata = train_data)
            R_squared_train_values[i] <- calculate_r_squared(train_data$logmeanConc, train_predictions)
            
            test_predictions <- predict(mlp_model, newdata = test_data)
            R_squared_test_values[i] <- calculate_r_squared(test_data$logmeanConc, test_predictions)
          }
        }
        
        val_train <- df[!df[[fold_column]] %in% c(k - 1, k), ] %>% select(-all_of(cols_to_remove))
        val_data  <- df[df[[fold_column]] %in% c(k - 1, k), ] %>% select(-all_of(cols_to_remove))
        
          model_epoch <- dnn(
            logmeanConc ~ .,
            hidden = hidden_val,
            batchsize = batchsize_val,
            data = val_train,
            epochs = 100,
            loss = "mae",
            burnin = 5,
            activation = "relu",
            lr = lr_val,
            device = "cuda"
          )
          
          val_pred <- predict(model_epoch, newdata = val_data)
          best_val_r2 <- calculate_r_squared(val_data$logmeanConc, val_pred)
        }
      }
      
      # Save to hyperparameter grid
      hyperparameters$MLP[[trait]][[desc]][[fold_column]]$train[j] <- mean(R_squared_train_values, na.rm = TRUE)
      hyperparameters$MLP[[trait]][[desc]][[fold_column]]$test[j]  <- mean(R_squared_test_values, na.rm = TRUE)
      hyperparameters$MLP[[trait]][[desc]][[fold_column]]$val[j]   <- best_val_r2
      hyperparameters$MLP[[trait]][[desc]][[fold_column]]$best_epoch[j] <- best_epoch  
    }
  }
}



#### Save results ####
saveRDS(hyperparameters$MLP$noTraits$MolDescriptor,file = "~/paper/EcoToxTM/Results/MLP_noTraits_MolDescriptor.rds")
saveRDS(hyperparameters$MLP$noTraits$Fingerprints,file = "~/paper/EcoToxTM/Results/MLP_noTraits_Fingerprints.rds")
saveRDS(hyperparameters$MLP$Traits$MolDescriptor,file = "~/paper/EcoToxTM/Results/MLP_Traits_MolDescriptor.rds")
saveRDS(hyperparameters$MLP$Traits$Fingerprints,file = "~/paper/EcoToxTM/Results/MLP_Traits_Fingerprints.rds")

