############################## Random Forest ###############################
# Set up four Random Forest models: traits & Mol. Descriptors /traits & Fingerprints/ no traits & Mol. Descriptors / no traits & Fingerprints
# Each model undergoes blocked validation: random, species and chemicals blocked  
library(readr)
library(dplyr)
library(ranger)
source("Scripts/functions.R") # for R squared calculation
df <- readRDS("Data/df_final&converted.rds")
df_FP <- df[, -c(1:206)]
df_MD <- df[, -c(228:393)]


#### Hyperparameters ####
# Grid 
hyper_grid <- function(mtry, max.depth) {
  expand.grid(mtry = mtry, max.depth = max.depth, train = NA, test = NA, val = NA, fold_column = NA, stringsAsFactors = FALSE)
}


# # Global object
# hyperparameters <- list(randomForest = list(
#   noTraits = list(
#     MolDescriptor  = hyper_grid(c(1, 2, 68, 206), c(3, 5, 10, 15)),
#     Fingerprints   = hyper_grid(c(1, 2, 55, 166), c(3, 5, 10, 15))
#   ),
#   Traits = list(
#     MolDescriptor  = hyper_grid(c(1, 2, 74, 223), c(3, 5, 10, 15)),
#     Fingerprints   = hyper_grid(c(1, 2, 61, 183), c(3, 5, 10, 15))
#   )
# ))
# Global object

hyperparameters <- list(randomForest = list(
  noTraits = list(
    MolDescriptor  = hyper_grid(c(1), c(1)),
    Fingerprints   = hyper_grid(c(1), c(1))
  ),
  Traits = list(
    MolDescriptor  = hyper_grid(c(1), c(1)),
    Fingerprints   = hyper_grid(c(1), c(1))
  )
))
#### Model ####
run_RF_model <- function(df, fold_column, descriptor_type, with_traits) {
  
  cat("\n Modelling starts with :", descriptor_type, "| Traits:", with_traits, "| Fold:", fold_column, "\n")
  
  k <- 10
  trait_key <- if (with_traits) "Traits" else "noTraits"
  hp_list <- hyperparameters$randomForest[[trait_key]][[descriptor_type]]
  
  # Debugging
  if (is.null(hp_list) || nrow(hp_list) == 0) {
    message("No hyperparameter for  ", trait_key, " / ", descriptor_type, " found.")
    return(NULL)
  }
  
  nhyper <- nrow(hp_list)
  
  trait_cols_to_remove <- c("BEl", "VEp", "REs", "OGp", "BLs", "PFs", "PFv", "CPt", "RMl",
                            "mem_1", "mem_2", "mem_3", "mem_4", "mem_5", "mem_6", "mem_7", "body_size")
  base_cols_to_remove <- c("randomFold", "chemicalFold", "speciesFold")
  cols_to_remove <- if (with_traits) base_cols_to_remove else c(base_cols_to_remove, trait_cols_to_remove)
  
  for (j in 1:nhyper) {
    mtry_val <- hp_list$mtry[j]
    max_depth_val <- hp_list$max.depth[j]
    
    R_squared_train_values <- numeric(k - 2)
    R_squared_test_values  <- numeric(k - 2)
    
    for (i in 1:(k - 2)) {
      train_data_full <- df[df[[fold_column]] != i, ]
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
    
    val_data_full <- df[df[[fold_column]] %in% c(k - 1, k), ]
    val_data <- val_data_full %>% select(-all_of(cols_to_remove))
    val_pred <- predict(rf_model, data = val_data)$predictions
    R_squared_val <- calculate_r_squared(val_data$logmeanConc, val_pred)
    
    hyperparameters$randomForest[[trait_key]][[descriptor_type]]$train[j] <- mean(R_squared_train_values, na.rm = TRUE)
    hyperparameters$randomForest[[trait_key]][[descriptor_type]]$test[j]  <- mean(R_squared_test_values, na.rm = TRUE)
    hyperparameters$randomForest[[trait_key]][[descriptor_type]]$val[j]   <- R_squared_val
    hyperparameters$randomForest[[trait_key]][[descriptor_type]]$fold_column[j] <- fold_column
    
    # Debugging 
    cat("Saved for", trait_key, " | ", descriptor_type, " | Fold: ", fold_column, "\n")
  }
  
  return(NULL)
}

#### Loop over all Combinations ####
combinations <- expand.grid(
  fold = c("randomFold", "speciesFold", "chemicalFold"),
  desc = c("MolDescriptor", "Fingerprints"),
  traits = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(combinations))) {
  combo <- combinations[i, ]
  # Debugging
  print(combo)
  df_input <- if (combo$desc == "MolDescriptor") df_MD else df_FP
  # Debugging: show df for combo
  print(head(df_input))  
  trait_label <- if (combo$traits) "Traits" else "noTraits"
  # Debugging: which traits_label?
  print(trait_label)
  # Debugging
  message("Model running for: ", combo$fold, " - ", combo$desc, " - ", trait_label)
  
  tryCatch({
    run_RF_model(
      df = df_input,
      fold_column = combo$fold,
      descriptor_type = combo$desc,
      with_traits = combo$traits
    )
  }, error = function(e) {
    message("Error ", combo$fold, " - ", combo$desc, " - ", trait_label)
    message("Details: ", e$message)
  })
  
  #Debugging
  message("completed", combo$fold, " - ", combo$desc, " - ", trait_label)
}

saveRDS(hyperparameters$randomForest$noTraits$MolDescriptor, file = "Results/RF_noTraits_MolDescriptor.rds")
saveRDS(hyperparameters$randomForest$noTraits$Fingerprints, file = "Results/RF_noTraits_Fingerprints.rds")
saveRDS(hyperparameters$randomForest$Traits$MolDescriptor,   file = "Results/RF_Traits_MolDescriptor.rds")
saveRDS(hyperparameters$randomForest$Traits$Fingerprints,    file = "Results/RF_Traits_Fingerprints.rds")

