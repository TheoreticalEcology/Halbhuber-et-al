############################## Random Forest ###############################
# Set up four Random Forest models: traits & Mol. Descriptors /traits & Fingerprints/ no traits & Mol. Descriptors / no traits & Fingerprints
# Each model undergoes blocked validation: random, species and chemicals blocked  
library(readr)
library(dplyr)
library(ranger)

##### Data #####
df <- readRDS("Data/df_final&converted.rds")
df_FP <- df[, -c(1:206)]
df_MD <- df[, -c(228:393)]

# Function to calculate R-squared
calculate_r_squared <- function(true_values, predictions) {
  res_ss <- sum((true_values - predictions)^2)
  total_ss <- sum((true_values - mean(true_values))^2)
  R_squared <- 1 - res_ss / total_ss
  return(R_squared)
}

##### Hyperparameter #####
hyperparameters = list()
hyperparameters$randomForest = list()
hyperparameters$randomForest$noTraits$MolDescriptor = expand.grid(mtry = c(1,2,68,206), max.depth = c(3,5,10,15), train = NA, test = NA, val = NA, CV = NA)
hyperparameters$randomForest$noTraits$Fingerprints = expand.grid(mtry = c(1,2,55,166), max.depth = c(3,5,10,15), train = NA, test = NA, val = NA, CV = NA)
hyperparameters$randomForest$Traits$MolDescriptor = expand.grid(mtry = c(1,2,74,223), max.depth = c(3,5,10,15), train = NA, test = NA, val = NA, CV = NA)
hyperparameters$randomForest$Traits$Fingerprints = expand.grid(mtry = c(1,2,61,183), max.depth = c(3,5,10,15), train = NA, test = NA, val = NA, CV = NA)

##### Model #####
run_RF_model <- function(df, fold_column, descriptor_type, with_traits = FALSE) {
  
  k <- 10
  trait_key <- if (with_traits) "Traits" else "noTraits"
  hp_list <- hyperparameters$randomForest[[trait_key]][[descriptor_type]]
  nhyper <- nrow(hp_list)
  
  R_squared_train_values <- numeric(k-2)
  R_squared_test_values <- numeric(k-2)
  R_squared_val_values <- numeric(2)
  
  # Columns that have to be removed, if traits_key == "noTraits" 
  trait_cols_to_remove <- c("BEl", "VEp", "REs", "OGp", "BLs", "PFs", "PFv", "CPt", "RMl",
                            "mem_1", "mem_2", "mem_3", "mem_4", "mem_5", "mem_6", "mem_7", "body_size")
  
  for (j in 1:nhyper) {
    mtry_val <- hp_list$mtry[j]
    max_depth_val <- hp_list$max.depth[j]
    
    for (i in 1:(k-2)) {
      train_data_full <- df[df[[fold_column]] != i, ]
      test_data_full <- df[df[[fold_column]] == i, ]
      
      # Cols that need overall need to be removed
      base_cols_to_remove <- c("randomFold", "chemicalFold", "speciesFold")
      cols_to_remove <- if (with_traits) base_cols_to_remove else c(base_cols_to_remove, trait_cols_to_remove)
      
      train_data <- train_data_full %>% select(-all_of(cols_to_remove))
      test_data  <- test_data_full  %>% select(-all_of(cols_to_remove))
      
      rf_model <- ranger(logmeanConc ~ ., data = train_data,
                         min.node.size = 5,
                         mtry = mtry_val,
                         max.depth = max_depth_val,
                         importance = "impurity")
      
      train_predictions <- predict(rf_model, data = train_data)$predictions
      R_squared_train_values[i] <- 1 - sum((train_data$logmeanConc - train_predictions)^2) / 
        sum((train_data$logmeanConc - mean(train_data$logmeanConc))^2)
      
      test_predictions <- predict(rf_model, data = test_data)$predictions
      R_squared_test_values[i] <- 1 - sum((test_data$logmeanConc - test_predictions)^2) /
        sum((test_data$logmeanConc - mean(test_data$logmeanConc))^2)
    }
    
    val_data_full <- df[df[[fold_column]] %in% c(k-1, k), ]
    val_data <- val_data_full %>% select(-all_of(cols_to_remove))
    
    val_predictions <- predict(rf_model, data = val_data)$predictions
    R_squared_val_values[j] <- 1 - sum((val_data$logmeanConc - val_predictions)^2) /
      sum((val_data$logmeanConc - mean(val_data$logmeanConc))^2)
    
    hyperparameters$randomForest[[trait_key]][[descriptor_type]]$train[j] <- mean(R_squared_train_values, na.rm = TRUE)
    hyperparameters$randomForest[[trait_key]][[descriptor_type]]$test[j] <- mean(R_squared_test_values, na.rm = TRUE)
    hyperparameters$randomForest[[trait_key]][[descriptor_type]]$val[j] <- mean(R_squared_val_values, na.rm = TRUE)
    hyperparameters$randomForest[[trait_key]][[descriptor_type]]$CV[j] <- fold_column
  }
  
  return(hyperparameters)
}

#### TRAITS ####
#### Mol. Descriptor ####
rsqrt_random_rf_MD <- run_RF_model(df_MD, "randomFold", "MolDescriptor", with_traits = FALSE)
rsqrt_chemical_rf_MD <- run_RF_model(df_MD, "chemicalFold", "MolDescriptor", with_traits = FALSE)
rsqrt_species_rf_MD <- run_RF_model(df_MD, "speciesFold", "MolDescriptor", with_traits = FALSE)

combined_MD <- rbind(rsqrt_random_rf_MD$randomForest$noTraits$MolDescriptor,
                     rsqrt_chemical_rf_MD$randomForest$noTraits$MolDescriptor,
                     rsqrt_species_rf_MD$randomForest$noTraits$MolDescriptor)

saveRDS(combined_MD, file = "Results/results_randomForest_noTraits_MD.rds")

##### NO TRAITS ####
rsqrt_random_rf_MD <- run_RF_model(df_MD, "randomFold", "MolDescriptor", with_traits = TRUE)
rsqrt_chemical_rf_MD <- run_RF_model(df_MD, "chemicalFold", "MolDescriptor", with_traits = TRUE)
rsqrt_species_rf_MD <- run_RF_model(df_MD, "speciesFold", "MolDescriptor", with_traits = TRUE)

combined_MD <- rbind(rsqrt_random_rf_MD$randomForest$Traits$MolDescriptor,
                     rsqrt_chemical_rf_MD$randomForest$Traits$MolDescriptor,
                     rsqrt_species_rf_MD$randomForest$Traits$MolDescriptor)

saveRDS(combined_MD, file = "Results/results_randomForest_Traits_MD.rds")
