############################## MLP ###############################
# Set up four MLP models: traits & Mol. Descriptors /traits & Fingerprints/ no traits & Mol. Descriptors / no traits & Fingerprints
# Each model undergoes blocked validation: random, species and chemicals blocked  
library(readr)
library(dplyr)
library(cito)

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
hyperparameters$dnn = list()
hyperparameters$dnn$noTraits = data.frame(lr = 1, batchsize = 50, train = NA, test = NA, val = NA, CV = "")
hyperparameters$dnn$Traits = data.frame(lr = 1, batchsize = 50, train = NA, test = NA, val = NA, CV = "")

########################      DNN         ################################
# MDs
# scale variables 
variables_to_scale <- c(colnames(df_MD[,1:223]))

scale <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

df_MD_scaled <- df_MD %>%
  mutate(across(all_of(variables_to_scale), scale))

### some columns show no variance in their values and therefore lead to NaN after scaling. 
df_MD_scaled<- df_MD_scaled[, colSums(is.na(df_MD_scaled)) == 0]

m_dnn_MD_TRAITS <- function(df, fold_column, k = 10, epochs = 2) {
  
  # Speichern des besten R²-Werts über alle Epochen für jeden Fold
  best_r_squared_train_values <- numeric(k)
  best_r_squared_test_values  <- numeric(k)
  
  for (i in 1:k) {
    # Trainings- und Testdaten
    train_data_full <- df[df[[fold_column]] != i, ]
    test_data_full  <- df[df[[fold_column]] == i, ]
    
    train_data <- train_data_full %>% select(-randomFold, -chemicalFold, -speciesFold)
    test_data  <- test_data_full  %>% select(-randomFold, -chemicalFold, -speciesFold)
    
    train_x <- as.matrix(select(train_data, -logmeanConc))
    train_y <- as.matrix(train_data$logmeanConc)
    test_x  <- as.matrix(select(test_data, -logmeanConc))
    test_y  <- as.matrix(test_data$logmeanConc)
    
    # Modell trainieren
    dnn_model <- dnn(
      logmeanConc ~ .,
      hidden = c(100L, 150L),
      batchsize = 50L,
      data = train_data,
      epochs = epochs,  # hier z.B. 300 Epochen
      loss = "mae",
      burnin = 5,
      activation = "relu",
      lr = 0.01,
      device = "cuda"
    )
    
    # Berechne R² für jede Epoche (absoluter Bestwert)
    r_squared_values_train <- numeric(epochs)
    r_squared_values_test  <- numeric(epochs)
    
    for (epoch in 1:epochs) {
      predictions_train <- predict(dnn_model, train_x)
      predictions_test  <- predict(dnn_model, test_x)
      
      r_squared_values_train[epoch] <- calculate_r_squared(train_y, predictions_train)
      r_squared_values_test[epoch]  <- calculate_r_squared(test_y, predictions_test)
    }
    
    # Besten R²-Wert für jedes Fold aus den Epochen wählen
    best_r_squared_train_values[i] <- max(r_squared_values_train)
    best_r_squared_test_values[i]  <- max(r_squared_values_test)
  }
  
  # Endgültige R²-Werte: besten Wert (max) aus allen Folds
  final_best_R_squared_train_max <- max(best_r_squared_train_values)
  final_best_R_squared_test_max  <- max(best_r_squared_test_values)
  
  # Optional: Median der besten Werte (falls du robustere Schätzung willst)
  final_best_R_squared_train_median <- median(best_r_squared_train_values)
  final_best_R_squared_test_median  <- median(best_r_squared_test_values)
  
  cat("Finaler bester R² (Train) Maximalwert: ", final_best_R_squared_train_max, "\n")
  cat("Finaler bester R² (Test) Maximalwert: ", final_best_R_squared_test_max, "\n")
  cat("Finaler bester R² (Train) Median: ", final_best_R_squared_train_median, "\n")
  cat("Finaler bester R² (Test) Median: ", final_best_R_squared_test_median, "\n")
  
  
  rsqrt_random_dnn_with_MD <- m_dnn_MD_TRAITS(df_MD_scaled, "randomFold")
  rsqrt_chemical_dnn_with_MD <- m_dnn_MD_TRAITS(df_MD_scaled, "chemicalFold")
  rsqrt_species_dnn_with_MD <- m_dnn_MD_TRAITS(df_MD_scaled, "speciesFold")
  
  results_r_squared_dnn_with <- as.data.frame(rbind(rsqrt_random_dnn_with_MD, rsqrt_chemical_dnn_with_MD, rsqrt_species_dnn_with_MD))
  
  saveRDS(results_r_squared_dnn_with_MD, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_dnn_with_MD_noRegularization.rds")
  
  #FP
  # scale species traits, FPs don't need to be scaled because they are binary
  variables_to_scale <- c(colnames(df_FP[,1:223]))
  
  scale <- function(x) {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  }
  
  df_FP_scaled <- df_FP %>%
    mutate(across(all_of(variables_to_scale), scale))
  
  ### some columns show no variance in their values and therefore lead to NaN after scaling. 
  df_FP_scaled<- df_FP_scaled[, colSums(is.na(df_FP_scaled)) == 0]
  
  m_dnn_FP_TRAITS <- function(df, fold_column, k = 10, epochs = 300) {
    
    R_squared_train_values <- numeric(k)
    R_squared_test_values <- numeric(k)
    
    for (i in 1:k) {
      train_data_full <- df[df[[fold_column]] != i, ]
      test_data_full <- df[df[[fold_column]] == i, ]
      
      train_data <- train_data_full %>% select(-randomFold, -chemicalFold, -speciesFold)
      test_data <- test_data_full %>% select(-randomFold, -chemicalFold, -speciesFold)
      
      train_x <- as.matrix(select(train_data, -logmeanConc))
      train_y <- as.matrix(train_data$logmeanConc)
      test_x <- as.matrix(select(test_data, -logmeanConc))
      test_y <- as.matrix(test_data$logmeanConc)
      
      dnn_model <- dnn(logmeanConc ~ .,
                       hidden = c(100L, 150L),
                       activation = "relu",
                       lr = 0.01,
                       #lambda = 0.001,
                       #alpha = 0.1,
                       device = "cuda")
      
      predictions_train <- predict(dnn_model, train_x)
      R_squared_train_values[i] <- calculate_r_squared(train_y, predictions_train)
      
      predictions_test <- predict(dnn_model, test_x)
      R_squared_test_values[i] <- calculate_r_squared(test_y, predictions_test)
    }
    
    
    mean_R_squared_train <- mean(R_squared_train_values, na.rm = TRUE) 
    mean_R_squared_test <- mean(R_squared_test_values, na.rm = TRUE)
    
    return(list(train_R_squared = mean_R_squared_train, test_R_squared = mean_R_squared_test))
  }
  
  
  rsqrt_random_dnn_with <- m_dnn_FP_TRAITS(df_FP_scaled, "randomFold")
  rsqrt_chemical_dnn_with <- m_dnn_FP_TRAITS(df_FP_scaled, "chemicalFold")
  rsqrt_species_dnn_with <- m_dnn_FP_TRAITS(df_FP_scaled, "speciesFold")
  
  results_r_squared_dnn_with <- as.data.frame(rbind(rsqrt_random_dnn_with, rsqrt_chemical_dnn_with, rsqrt_species_dnn_with))
  
  saveRDS(results_r_squared_dnn_with, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_dnn_with_MD_noRegularization.rds")
  
  
  ######################## DNN without ################################
  
  # df und calculate_r_squared function von oben
  
  m_dnn_MD <- function(df, fold_column, k = 10, epochs = 500) {
    
    R_squared_train_values <- numeric(k)
    R_squared_test_values <- numeric(k)
    
    for (i in 1:k) {
      
      train_data_full <- df[df[[fold_column]] != i, ]
      test_data_full <- df[df[[fold_column]] == i, ]
      
      
      train_data <- train_data_full %>% select(-randomFold, -chemicalFold, -speciesFold, -BEl, -VEp, -REs, -OGp, -BLs, -PFs, -PFv, -CPt, -RMl, -mem_1, -mem_2, -mem_3, -mem_4, -mem_5, -mem_6, -mem_7, -body_size)
      test_data <- test_data_full %>% select(-randomFold, -chemicalFold, -speciesFold, -BEl, -VEp, -REs, -OGp, -BLs, -PFs, -PFv, -CPt, -RMl, -mem_1, -mem_2, -mem_3, -mem_4, -mem_5, -mem_6, -mem_7, -body_size)
      
      
      train_x <- as.matrix(select(train_data, -logmeanConc))
      train_y <- as.matrix(train_data$logmeanConc)
      test_x <- as.matrix(select(test_data, -logmeanConc))
      test_y <- as.matrix(test_data$logmeanConc)
      
      
      dnn_model <- dnn(logmeanConc~ .,
                       hidden = c(100L, 150L),
                       batchsize = 50L,
                       data = train_data,
                       epochs = epochs,
                       loss = 'mae',
                       burnin = 5,
                       activation = "relu",
                       lr = 0.01,
                       device = "cuda")
      # #validation = 0.2,
      # early_stopping = TRUE,
      # hidden = c(20L, 20L),
      # batchsize = 300L,
      # data = train_data,
      # epochs = epochs,
      # loss = 'mae',
      # lr = 0.01,
      # burnin = 200,
      # activation = "relu",
      # #dropout= 0.2, # Wahrscheinlichkeit eines dropouts pro layer
      #lambda = 0.001,
      #alpha = 0.1)
      # #lr_scheduler = config_lr_scheduler('step',gamma = 0.7, step_size= 10),
      # device = "cuda")
      predictions_train <- predict(dnn_model, train_x)
      R_squared_train_values[i] <- calculate_r_squared(train_y, predictions_train)
      
      predictions_test <- predict(dnn_model, test_x)
      R_squared_test_values[i] <- calculate_r_squared(test_y, predictions_test)
    }
    
    
    mean_R_squared_train <- mean(R_squared_train_values, na.rm = TRUE)
    mean_R_squared_test <- mean(R_squared_test_values, na.rm = TRUE)
    
    return(list(train_R_squared = mean_R_squared_train, test_R_squared = mean_R_squared_test))
  }
  
  rsqrt_random_dnn_without <- m4(df_fish_dnn, "randomFold")
  rsqrt_chemical_dnn_without <- m4(df_fish_dnn, "chemicalFold")
  rsqrt_species_dnn_without <- m4(df_fish_dnn, "speciesFold")
  
  results_r_squared_dnn_without <- as.data.frame(rbind(rsqrt_random_dnn_without, rsqrt_chemical_dnn_without, rsqrt_species_dnn_without))
  
  saveRDS(results_r_squared_dnn_without, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_dnn_without_MACCS_10.rds")
  
  
  