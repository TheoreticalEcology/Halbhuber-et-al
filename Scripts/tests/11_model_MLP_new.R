############################## MLP ###############################
# Set up four MLP models: traits & Mol. Descriptors /traits & Fingerprints/ no traits & Mol. Descriptors / no traits & Fingerprints
# Each model undergoes blocked validation: random, species and chemicals blocked  
library(readr)
library(dplyr)
library(cito)
source("Scripts/functions.R") # for R squared calculation

##### Data #####
df <- readRDS("Data/df_final&converted.rds")
df_FP <- df[, -c(1:206)]
df_MD <- df[, -c(228:393)]

# scale variables 
variables_to_scale <- c(colnames(df_MD[,1:223]))
scale <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

df_MD <- df_MD %>%
  mutate(across(all_of(variables_to_scale), scale))

### some columns show no variance in their values and therefore lead to NaN after scaling. 
df_MD <- df_MD[, colSums(is.na(df_MD)) == 0]

##### Hyperparameter #####
# hyperparameters = list()
# hyperparameters$MLP = list()
# hyperparameters$MLP$noTraits$MolDescriptor = expand.grid(hidden = c(c(64, 64), c(128, 64), c(128, 128, 64), c(100L, 150L)), lr = c(0.01, 0.001, 0.0005), batchsize = c(32, 50, 64, 128), train = NA, test = NA, val = NA, CV = NA)
# hyperparameters$MLP$noTraits$Fingerprints = expand.grid(hidden = c(c(64, 64), c(128, 64), c(128, 128, 64), c(100L, 150L)), lr = c(0.01, 0.001, 0.0005), batchsize = c(32, 50, 64, 128), train = NA, test = NA, val = NA, CV = NA)
# hyperparameters$MLP$Traits$MolDescriptor = expand.grid(hidden = c(c(64, 64), c(128, 64), c(128, 128, 64), c(100L, 150L)), lr = c(0.01, 0.001, 0.0005), batchsize = c(32, 50, 64, 128), train = NA, test = NA, val = NA, CV = NA)
# hyperparameters$MLP$Traits$Fingerprints = expand.grid(hidden = c(c(64, 64), c(128, 64), c(128, 128, 64), c(100L, 150L)), lr = c(0.01, 0.001, 0.0005), batchsize = c(32, 50, 64, 128), train = NA, test = NA, val = NA, CV = NA)

# grid, if the same hyperparameters are used 
# mlp_grid <- expand.grid(hidden = c(c(64, 64),c(128, 64),c(128, 128, 64),c(100L, 150L)), lr = c(0.01, 0.001, 0.0005), batchsize = c(32, 50, 64, 128),
#                         train = NA, test = NA, val = NA, fold_column = NA)

mlp_grid <- expand.grid(hidden = c(c(64, 64),c(128, 64)), lr = c(0.01), batchsize = c(32),
                        train = NA, test = NA, val = NA, fold_column = NA)
hyperparameters <- list()
hyperparameters$MLP <- list()
hyperparameters$MLP$noTraits <- list()
hyperparameters$MLP$Traits <- list()

hyperparameters$MLP$noTraits$MolDescriptor = mlp_grid
hyperparameters$MLP$noTraits$Fingerprints = mlp_grid
hyperparameters$MLP$Traits$MolDescriptor = mlp_grid
hyperparameters$MLP$Traits$Fingerprints = mlp_grid


######################## Model ################################

run_MLP_model <- function(df, fold_column, descriptor_type, with_traits = FALSE) {
  
  k <- 10
  trait_key <- if (with_traits) "Traits" else "noTraits"
  hp_list <- hyperparameters$MLP[[trait_key]][[descriptor_type]]
  nhyper <- nrow(hp_list)
  
  # Container for R²-values
  R_squared_train_values <- numeric(k-2)
  R_squared_test_values <- numeric(k-2)
  R_squared_val_values <- numeric(2)
  
  # Columns to remove when "noTraits" is selected
  trait_cols_to_remove <- c("BEl", "VEp", "REs", "OGp", "BLs", "PFs", "PFv", "CPt", "RMl",
                            "mem_1", "mem_2", "mem_3", "mem_4", "mem_5", "mem_6", "mem_7", "body_size")
  
  for (j in 1:nhyper) {
    hidden_val <- hp_list$hidden[j]
    lr_val <- hp_list$lr[j]
    batchsize_val <- hp_list$batchsize[j]
    
    for (i in 1:(k-2)) {
      train_data_full <- df[df[[fold_column]] != i, ]
      test_data_full <- df[df[[fold_column]] == i, ]
      
      base_cols_to_remove <- c("randomFold", "chemicalFold", "speciesFold")
      cols_to_remove <- if (with_traits) base_cols_to_remove else c(base_cols_to_remove, trait_cols_to_remove)
      
      train_data <- train_data_full %>% select(-all_of(cols_to_remove))
      test_data  <- test_data_full  %>% select(-all_of(cols_to_remove))
      
      dnn_model <- dnn(
        logmeanConc ~ .,
        hidden = hidden_val,
        batchsize = batchsize_val,
        data = train_data,
        epochs = 3,  
        loss = "mae",
        burnin = 5,
        activation = "relu",
        lr = lr_val,
        device = "cuda"
      )
      
      # Prediction for test and train data
      train_predictions <- predict(dnn_model, newdata = train_data)
      R_squared_train_values[i] <- calculate_r_squared(train_data$logmeanConc, train_predictions)
      
      test_predictions <- predict(dnn_model, newdata = test_data)
      R_squared_test_values[i] <- calculate_r_squared(test_data$logmeanConc, test_predictions)
    }
    
    # R² for Validation data 
    val_data_full <- df[df[[fold_column]] %in% c(k-1, k), ]
    val_data <- val_data_full %>% select(-all_of(cols_to_remove))
    
    val_predictions <- predict(dnn_model, newdata = val_data)
    R_squared_val_values[j] <- calculate_r_squared(val_data$logmeanConc, val_predictions)

    hyperparameters$MLP[[trait_key]][[descriptor_type]]$train[j] <- mean(R_squared_train_values, na.rm = TRUE)
    hyperparameters$MLP[[trait_key]][[descriptor_type]]$test[j] <- mean(R_squared_test_values, na.rm = TRUE)
    hyperparameters$MLP[[trait_key]][[descriptor_type]]$val[j] <- mean(R_squared_val_values, na.rm = TRUE)
    hyperparameters$MLP[[trait_key]][[descriptor_type]]$fold_column[j]  <- fold_column
  }
  
  return(hyperparameters)
}

##### Loop over all combinations ####

fold_columns     <- c("randomFold", "chemicalFold", "speciesFold")
descriptor_types <- c("MolDescriptor", "Fingerprints")

for (fold in fold_columns) {
  for (desc in descriptor_types) {
    message("Running: ", fold, " - ", desc, " - noTraits")
    hyperparameters <- run_MLP_model(df_MD, fold_column = fold, descriptor_type = desc, with_traits = FALSE)
    
    message("Running: ", fold, " - ", desc, " - Traits")
    hyperparameters <- run_MLP_model(df_MD, fold_column = fold, descriptor_type = desc, with_traits = TRUE)
  }
}


saveRDS(hyperparameters$MLP$noTraits$MolDescriptor, file = "Results/MLP_noTraits_MolDescriptor.rds")
saveRDS(hyperparameters$MLP$noTraits$Fingerprints, file = "Results/MLP_noTraits_Fingerprints.rds")
saveRDS(hyperparameters$MLP$Traits$MolDescriptor, file = "Results/MLP_Traits_MolDescriptor.rds")
saveRDS(hyperparameters$MLP$Traits$Fingerprints, file = "Results/MLP_Traits_Fingerprints.rds")

  
  