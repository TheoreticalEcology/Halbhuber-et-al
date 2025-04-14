############################## Random Forest ###############################
# Set up four Random Forest models: traits & Mol. Descriptors /traits & Fingerprints/ no traits & Mol. Descriptors / no traits & Fingerprints
# Each model undergoes blocked validation: random, species and chemicals blocked  
library(readr)
library(dplyr)
library(ranger)
source("~/paper/EcoToxTM/Scripts/functions.R") # for R squared calculation
allData <- readRDS("~/paper/EcoToxTM/Data/df_final&converted.rds")
# Folds from factor to numeric
allData[, 225:227] <- sapply(allData[, 225:227], function(x) as.numeric(as.character(x)))

dat = list(
  df_FP <- allData[, -c(1:206)],
  df_MD <- allData[, -c(228:393)]  
)

#### Hyperparameter Grid ####
hyper_grid <- function(mtry, max.depth) {
  expand.grid(
    mtry = mtry,
    max.depth = max.depth,
    train = NA,
    test = NA,
    val = NA,
    stringsAsFactors = FALSE
  )
}

folds <- c("randomFold", "chemicalFold", "speciesFold")
descriptors <- c("Fingerprints", "MolDescriptor")
trait_levels <- c("noTraits", "Traits")

# hyper_grid_objects
hyper_grid_noTraits_FP <- hyper_grid(c(1, 2, 55, 166), c(3, 5, 10, 15))
hyper_grid_noTraits_MD <- hyper_grid(c(1, 2, 68, 206), c(3, 5, 10, 15))

hyper_grid_Traits_FP <- hyper_grid(c(1, 2, 61, 183), c(3, 5, 10, 15))
hyper_grid_Traits_MD <- hyper_grid(c(1, 2, 74, 223), c(3, 5, 10, 15))

hyperparameters <- list(randomForest = list())

for (trait in trait_levels) {
  hyperparameters$randomForest[[trait]] <- list()
  
  for (desc in descriptors) {
    hyperparameters$randomForest[[trait]][[desc]] <- list()
    
    # Select the correct hyper_grid based on trait and descriptor
    if (trait == "noTraits" && desc == "Fingerprints") {
      hyper_grid_object <- hyper_grid_noTraits_FP
    } else if (trait == "noTraits" && desc == "MolDescriptor") {
      hyper_grid_object <- hyper_grid_noTraits_MD
    } else if (trait == "Traits" && desc == "Fingerprints") {
      hyper_grid_object <- hyper_grid_Traits_FP
    } else if (trait == "Traits" && desc == "MolDescriptor") {
      hyper_grid_object <- hyper_grid_Traits_MD
    }
    
    # Loop through all folds
    for (fold in folds) {
      hyperparameters$randomForest[[trait]][[desc]][[fold]] <- hyper_grid_object
    }
  }
}


#### Model-Loop ####

for (desc in descriptors) {
  df <- if (desc == "MolDescriptor") df_MD else df_FP
  
  for (fold_column in folds) {
    for (trait in trait_levels) {
      
      k <- 10
      hp_list <- hyperparameters$randomForest[[trait]][[desc]][[fold_column]]
      nhyper <- nrow(hp_list)
      
      trait_cols_to_remove <- c("BEl", "VEp", "REs", "OGp", "BLs", "PFs", "PFv", "CPt", "RMl",
                                "mem_1", "mem_2", "mem_3", "mem_4", "mem_5", "mem_6", "mem_7", "body_size")
      base_cols_to_remove <- c("randomFold", "chemicalFold", "speciesFold")
      cols_to_remove <- if (trait == "Traits") base_cols_to_remove else c(base_cols_to_remove, trait_cols_to_remove)
      
      # Phase 1: Training/Test on folds 1:(k-2)
      for (j in 1:nhyper) {
        mtry_val <- hp_list$mtry[j]
        max_depth_val <- hp_list$max.depth[j]
        
        R_squared_train_values <- rep(NA, (k - 2))
        R_squared_test_values  <- rep(NA, (k - 2))
        
        for (i in 1:(k - 2)) {
          train_data_full <- df[df[[fold_column]] != i & df[[fold_column]] <= (k - 2), ]
          test_data_full  <- df[df[[fold_column]] == i, ]
          
          train_data <- train_data_full %>% select(-all_of(cols_to_remove))
          test_data  <- test_data_full %>% select(-all_of(cols_to_remove))
          
          rf_model <- ranger(logmeanConc ~ ., data = train_data,
                             min.node.size = 5,
                             mtry = mtry_val,
                             max.depth = max_depth_val,
                             importance = "impurity")
          
          train_pred <- predict(rf_model, data = train_data)$predictions
          test_pred  <- predict(rf_model, data = test_data)$predictions
          
          R_squared_train_values[i] <- calculate_r_squared(train_data$logmeanConc, train_pred)
          R_squared_test_values[i]  <- calculate_r_squared(test_data$logmeanConc, test_pred)
        }
        
        # save
        hyperparameters$randomForest[[trait]][[desc]][[fold_column]]$train[j] <- mean(R_squared_train_values, na.rm = TRUE)
        hyperparameters$randomForest[[trait]][[desc]][[fold_column]]$test[j]  <- mean(R_squared_test_values, na.rm = TRUE)
      }
      
      # Phase 2: Select best hyperparameter-set (highest test-R²)
      best_index <- which.max(hyperparameters$randomForest[[trait]][[desc]][[fold_column]]$test)
      best_mtry <- hp_list$mtry[best_index]
      best_depth <- hp_list$max.depth[best_index]
      
      # Train on folds 1:(k-2), evaluate on folds k-1 und k
      val_train <- df[df[[fold_column]] <= (k - 2), ]
      val_data  <- df[df[[fold_column]] %in% c(k - 1, k), ]
      
      val_train <- val_train %>% select(-all_of(cols_to_remove))
      val_data  <- val_data %>% select(-all_of(cols_to_remove))
      
      rf_model_val <- ranger(logmeanConc ~ ., data = val_train,
                             min.node.size = 5,
                             mtry = best_mtry,
                             max.depth = best_depth,
                             importance = "impurity")
      
      val_pred <- predict(rf_model_val, data = val_data)$predictions
      R_squared_val <- calculate_r_squared(val_data$logmeanConc, val_pred)
      
      # Save final Val-values with best hyperparameters
      hyperparameters$randomForest[[trait]][[desc]][[fold_column]]$val <- rep(NA, nhyper)
      hyperparameters$randomForest[[trait]][[desc]][[fold_column]]$val[best_index] <- R_squared_val
      
      cat("Final eval for", trait, "|", desc, "|", fold_column, "| Best combination: mtry =", best_mtry, ", depth =", best_depth, "\n")
    }
  }
}


#### Save results ####
saveRDS(hyperparameters$randomForest$noTraits$MolDescriptor, file = "~/paper/EcoToxTM/Results/RF_noTraits_MolDescriptor.rds")
saveRDS(hyperparameters$randomForest$noTraits$Fingerprints, file = "~/paper/EcoToxTM/Results/RF_noTraits_Fingerprints.rds")
saveRDS(hyperparameters$randomForest$Traits$MolDescriptor,   file = "~/paper/EcoToxTM/Results/RF_Traits_MolDescriptor.rds")
saveRDS(hyperparameters$randomForest$Traits$Fingerprints,    file = "~/paper/EcoToxTM/Results/RF_Traits_Fingerprints.rds")