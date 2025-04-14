# Set up four models: 1. RF without Species traits 2. RF with species traits 3. DNN without Species traits 2. DNN with species traits 
# each model undergoes validation: random, species and chemicals blocked  
# chemical representation: MACCS fingerprints (FP) and Molecular Descriptors (MD)

#################################################################################
##############################  Random Forest ###################################
#################################################################################
library(readr)
library(dplyr)
library(ranger)
library(cito)

df <- readRDS("Data/df_final.rds")

#################################################################################
# Data preparation for the models 
# convert fingerprints as character strings to numeric representation as a matrix

convert_to_vector <- function(fingerprint) {as.integer(unlist(strsplit(fingerprint, split = "")))}

fingerprints_list <- lapply(df$Fingerprint, convert_to_vector)

fingerprints_matrix <- do.call(rbind, fingerprints_list)

df <- cbind(df, fingerprints_matrix)


# give each fingerprint bit a categorical name (numeric not proccessable in ranger)
numbers <- as.character(1:166)
new_colnames <- paste0("bit", numbers)
colnames(df)[233:398] <- new_colnames

df <- as.data.frame(df)
df <- df %>% select (-Fingerprint, -CAS, -SMILES, -scientificNameStd, -meanConc)

# ## How many 0en are in each fold --> nearly 30 perscent
# sapply(1:10, function(fold) {
# mean(sapply(df_fish_rf[df_fish_rf$chemicalFold!=fold,-c(1:21, 902:905)], function(p) mean(p)) == 0)
# })
# table(df_fish_rf[df_fish_rf$chemicalFold!=1,-(1:21)]$bit831)
# mean(table(df_fish_rf[df_fish_rf$chemicalFold==1,-(1:21)]$bit831))
#
# ## How many 0en are in each fold --> nearly 30 perscent
# sapply(1:10, function(fold) {
#   mean(sapply(xc[df_fish_rf$chemicalFold!=fold,-c(1:21)], function(p) mean(p)) == 0)
# })
# table(df_fish_rf[df_fish_rf$chemicalFold!=1,-(1:21)]$bit831)
# mean(table(df_fish_rf[df_fish_rf$chemicalFold==1,-(1:21)]$bit831))
#
# hist(df_fish_rf[df_fish_rf$chemicalFold==1,5])
#
#
# dim(df_fish_rf[!duplicated(df_fish_rf$VEp),])


#################################################################################
################################# MODELLING #####################################
#################################################################################
## test: no chemical traits, only species traits 
#df_fish_rf <- df_fish_rf %>% select(logmeanConc, randomFold, chemicalFold, speciesFold, BEl, VEp, REs, OGp, BLs, PFs, PFv, CPt, RMl, mem_1, mem_2, mem_3, mem_4, mem_5, mem_6, mem_7, body_size)

# Function to calculate R-squared

calculate_r_squared <- function(true_values, predictions) {
  res_ss <- sum((true_values - predictions)^2)
  total_ss <- sum((true_values - mean(true_values))^2)
  R_squared <- 1 - res_ss / total_ss
  return(R_squared)
}

###################################################################################################
############################## Random Forest without species traits ###############################
###################################################################################################
# MDs 
# exclude Fingerprints
df1 <- df[, -c(228:393)]

m1 <- function(df1, fold_column, k = 10) {

  R_squared_train_values <- numeric(k)
  R_squared_test_values <- numeric(k)

  for (i in 1:k) {

    train_data_full <- df1[df1[[fold_column]] != i, ]
    test_data_full <- df1[df1[[fold_column]] == i, ]

    train_data <- train_data_full %>% select(-randomFold, -chemicalFold, -speciesFold, -BEl, -VEp, -REs, -OGp, -BLs, -PFs, -PFv, -CPt, -RMl, -mem_1, -mem_2, -mem_3, -mem_4, -mem_5, -mem_6, -mem_7, -body_size)
    test_data <- test_data_full %>% select(-randomFold, -chemicalFold, -speciesFold, -BEl, -VEp, -REs, -OGp, -BLs, -PFs, -PFv, -CPt, -RMl, -mem_1, -mem_2, -mem_3, -mem_4, -mem_5, -mem_6, -mem_7, -body_size)

    rf_model <- ranger(logmeanConc ~ .,
                       data = train_data,
                       min.node.size = 5,       
                       mtry = 206,  # Var per split 186 - 20
                       max.depth = 15,    
                       #num.trees = 1000,           
                       importance = "impurity")
    #



    train_predictions <- predict(rf_model, data = train_data)$predictions
    res_ss_train <- sum((train_data$logmeanConc - train_predictions)^2)
    total_ss_train <- sum((train_data$logmeanConc - mean(train_data$logmeanConc))^2)
    R_squared_train <- 1 - res_ss_train / total_ss_train
    R_squared_train_values[i] <- R_squared_train


    test_predictions <- predict(rf_model, data = test_data)$predictions
    res_ss_test <- sum((test_data$logmeanConc - test_predictions)^2)
    total_ss_test <- sum((test_data$logmeanConc - mean(test_data$logmeanConc))^2)
    R_squared_test <- 1 - res_ss_test / total_ss_test
    R_squared_test_values[i] <- R_squared_test
  }


  mean_R_squared_train <- mean(R_squared_train_values, na.rm = TRUE)
  mean_R_squared_test <- mean(R_squared_test_values, na.rm = TRUE)

  return(list(train_R_squared = mean_R_squared_train, test_R_squared = mean_R_squared_test,
              all_train_R_squared = R_squared_train_values, all_test_R_squared = R_squared_test_values))
}



rsqrt_random_rf_without <- m1(df1, "randomFold")
rsqrt_chemical_rf_without <- m1(df1, "chemicalFold")
rsqrt_species_rf_without <- m1(df1, "speciesFold")


results_r_squared_rf_without <- as.data.frame(rbind(rsqrt_random_rf_without, rsqrt_chemical_rf_without, rsqrt_species_rf_without))

#colnames(results_r_squared_rf_without) <- c("Train_R_squared", "Test_R_squared")

saveRDS(results_r_squared_rf_without, file = "/home/isabellehalbhuber/Toxicology/Results/results_randomForest_without_MD.rds")

# FP 
# exclude MD
df2 <- df[, -c(1:206)]

m1 <- function(df2, fold_column, k = 10) {
  
  R_squared_train_values <- numeric(k)
  R_squared_test_values <- numeric(k)
  
  for (i in 1:k) {
    
    train_data_full <- df2[df2[[fold_column]] != i, ]
    test_data_full <- df2[df2[[fold_column]] == i, ]
    
    train_data <- train_data_full %>% select(-randomFold, -chemicalFold, -speciesFold, -BEl, -VEp, -REs, -OGp, -BLs, -PFs, -PFv, -CPt, -RMl, -mem_1, -mem_2, -mem_3, -mem_4, -mem_5, -mem_6, -mem_7, -body_size)
    test_data <- test_data_full %>% select(-randomFold, -chemicalFold, -speciesFold, -BEl, -VEp, -REs, -OGp, -BLs, -PFs, -PFv, -CPt, -RMl, -mem_1, -mem_2, -mem_3, -mem_4, -mem_5, -mem_6, -mem_7, -body_size)
    
    rf_model <- ranger(logmeanConc ~ .,
                       data = train_data,
                       min.node.size = 5,       
                       mtry = 206,  # Var per split 186 - 20
                       max.depth = 15,    
                       #num.trees = 1000,           
                       importance = "impurity")
    #
    
    
    
    train_predictions <- predict(rf_model, data = train_data)$predictions
    res_ss_train <- sum((train_data$logmeanConc - train_predictions)^2)
    total_ss_train <- sum((train_data$logmeanConc - mean(train_data$logmeanConc))^2)
    R_squared_train <- 1 - res_ss_train / total_ss_train
    R_squared_train_values[i] <- R_squared_train
    
    
    test_predictions <- predict(rf_model, data = test_data)$predictions
    res_ss_test <- sum((test_data$logmeanConc - test_predictions)^2)
    total_ss_test <- sum((test_data$logmeanConc - mean(test_data$logmeanConc))^2)
    R_squared_test <- 1 - res_ss_test / total_ss_test
    R_squared_test_values[i] <- R_squared_test
  }
  
  
  mean_R_squared_train <- mean(R_squared_train_values, na.rm = TRUE)
  mean_R_squared_test <- mean(R_squared_test_values, na.rm = TRUE)
  
  return(list(train_R_squared = mean_R_squared_train, test_R_squared = mean_R_squared_test,
              all_train_R_squared = R_squared_train_values, all_test_R_squared = R_squared_test_values))
}



rsqrt_random_rf_without <- m1(df2, "randomFold")
rsqrt_chemical_rf_without <- m1(df2, "chemicalFold")
rsqrt_species_rf_without <- m1(df2, "speciesFold")


results_r_squared_rf_without <- as.data.frame(rbind(rsqrt_random_rf_without, rsqrt_chemical_rf_without, rsqrt_species_rf_without))

#colnames(results_r_squared_rf_without) <- c("Train_R_squared", "Test_R_squared")

saveRDS(results_r_squared_rf_without, file = "/home/isabellehalbhuber/Toxicology/Results/results_randomForest_without_FP.rds")

# #################################################################################
# ######################  Random Forest With Traits ###############################
# #################################################################################
# # 2. RF with species traits
# MDs 

m2 <- function(df1, fold_column, k = 10) {

  R_squared_train_values <- numeric(k)
  R_squared_test_values <- numeric(k)

  for (i in 1:k) {

    train_data_full <- df1[df1[[fold_column]] != i, ]
    test_data_full <- df1[df1[[fold_column]] == i, ]

    train_data <- train_data_full %>% select(-randomFold, -chemicalFold, -speciesFold)
    test_data <- test_data_full %>% select(-randomFold, -chemicalFold, -speciesFold)


    rf_model <- ranger(logmeanConc ~ .,
                       data = train_data,
                       min.node.size = 5,       
                       mtry = 223,  
                       max.depth = 15,            
                       #num.trees = 1000,           
                       importance = "impurity")

    train_predictions <- predict(rf_model, data = train_data)$predictions
    res_ss_train <- sum((train_data$logmeanConc - train_predictions)^2)
    total_ss_train <- sum((train_data$logmeanConc - mean(train_data$logmeanConc))^2)
    R_squared_train <- 1 - res_ss_train / total_ss_train
    R_squared_train_values[i] <- R_squared_train


    test_predictions <- predict(rf_model, data = test_data)$predictions
    res_ss_test <- sum((test_data$logmeanConc - test_predictions)^2)
    total_ss_test <- sum((test_data$logmeanConc - mean(test_data$logmeanConc))^2)
    R_squared_test <- 1 - res_ss_test / total_ss_test
    R_squared_test_values[i] <- R_squared_test
  }


  mean_R_squared_train <- mean(R_squared_train_values, na.rm = TRUE)
  mean_R_squared_test <- mean(R_squared_test_values, na.rm = TRUE)

  return(list(train_R_squared = mean_R_squared_train, test_R_squared = mean_R_squared_test,
              all_train_R_squared = R_squared_train_values, all_test_R_squared = R_squared_test_values))
}


rsqrt_random_rf_with <- m2(df1, "randomFold")
rsqrt_chemical_rf_with <- m2(df1, "chemicalFold") # chemikalien mehr out of distribution # roberts paper, strukturen, blocking... wahren prediction fehler finden... darum gehts. chemkalien vielleicht mehr out of sample... bei predicitive performance ertmal anova aufsetzen
rsqrt_species_rf_with <- m2(df1, "speciesFold")

results_r_squared_rf_with <- as.data.frame(rbind(rsqrt_random_rf_with, rsqrt_chemical_rf_with, rsqrt_species_rf_with))


saveRDS(results_r_squared_rf_with, file = "/home/isabellehalbhuber/Toxicology/Results/results_randomForest_with_MD.rds")

##### FP
m2 <- function(df2, fold_column, k = 10) {
  
  R_squared_train_values <- numeric(k)
  R_squared_test_values <- numeric(k)
  
  for (i in 1:k) {
    
    train_data_full <- df2[df2[[fold_column]] != i, ]
    test_data_full <- df2[df2[[fold_column]] == i, ]
    
    train_data <- train_data_full %>% select(-randomFold, -chemicalFold, -speciesFold)
    test_data <- test_data_full %>% select(-randomFold, -chemicalFold, -speciesFold)
    
    
    rf_model <- ranger(logmeanConc ~ .,
                       data = train_data,
                       min.node.size = 5,       
                       mtry = 176,  
                       max.depth = 15,            
                       #num.trees = 1000,           
                       importance = "impurity")
    
    train_predictions <- predict(rf_model, data = train_data)$predictions
    res_ss_train <- sum((train_data$logmeanConc - train_predictions)^2)
    total_ss_train <- sum((train_data$logmeanConc - mean(train_data$logmeanConc))^2)
    R_squared_train <- 1 - res_ss_train / total_ss_train
    R_squared_train_values[i] <- R_squared_train
    
    
    test_predictions <- predict(rf_model, data = test_data)$predictions
    res_ss_test <- sum((test_data$logmeanConc - test_predictions)^2)
    total_ss_test <- sum((test_data$logmeanConc - mean(test_data$logmeanConc))^2)
    R_squared_test <- 1 - res_ss_test / total_ss_test
    R_squared_test_values[i] <- R_squared_test
  }
  
  
  mean_R_squared_train <- mean(R_squared_train_values, na.rm = TRUE)
  mean_R_squared_test <- mean(R_squared_test_values, na.rm = TRUE)
  
  return(list(train_R_squared = mean_R_squared_train, test_R_squared = mean_R_squared_test,
              all_train_R_squared = R_squared_train_values, all_test_R_squared = R_squared_test_values))
}


rsqrt_random_rf_with <- m2(df2, "randomFold")
rsqrt_chemical_rf_with <- m2(df2, "chemicalFold") 
rsqrt_species_rf_with <- m2(df2, "speciesFold")

results_r_squared_rf_with <- as.data.frame(rbind(rsqrt_random_rf_with, rsqrt_chemical_rf_with, rsqrt_species_rf_with))


saveRDS(results_r_squared_rf_with, file = "/home/isabellehalbhuber/Toxicology/Results/results_randomForest_with_FP.rds")


##########################################################################
############################# DNN with ###################################
##########################################################################


#df_fish_dnn <- readRDS("/home/isabellehalbhuber/Toxicology/Data/df_fish_rf.rds")


### scale variables for pubchem fingerprints and AlogP
# # variables_to_scale <- c("BEl", "VEp", "REs", "OGp",
# #                         "RMl", "BLs", "PFv", "PFs", "CPt", "body_size",
# #                         "mem_1", "mem_2", "mem_3", "mem_4", "mem_5", "mem_6", "mem_7")
# 
# scale <- function(x) {
#   return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
# }
# df_fish_dnn <- df_fish_dnn %>%
#   mutate(across(all_of(variables_to_scale), scale))



variables_to_scale <- c("BEl", "VEp", "REs", "OGp",
                        "RMl", "BLs", "PFv", "PFs", "CPt", "body_size",
                        "mem_1", "mem_2", "mem_3", "mem_4", "mem_5", "mem_6", "mem_7")

scale <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

df_dnn <- df %>%
  mutate(across(all_of(variables_to_scale), scale))

df_dnn <- df_dnn[, -c(228:393)]

variables_to_scale <- c(colnames(df_dnn[,1:206]))
scale <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

df_dnn <- df_dnn %>%
  mutate(across(all_of(variables_to_scale), scale))




## model MDDs

m3 <- function(df, fold_column, k = 10, epochs = 300) {

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
                     batchsize = 50L,
                     data = train_data,
                     epochs = epochs,
                     loss = 'mae',
                     burnin = 5,
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


  mean_R_squared_train <- mean(R_squared_train_values, na.rm = TRUE) ## hier nicht den mean nehmen, sondern die liste ausgeben lassen 
  mean_R_squared_test <- mean(R_squared_test_values, na.rm = TRUE)

  return(list(train_R_squared = mean_R_squared_train, test_R_squared = mean_R_squared_test))
}


rsqrt_random_dnn_with <- m3(df_dnn[, !sapply(df_dnn, function(x) any(is.na(x)))], "randomFold", k = 1, epochs = 300L)
rsqrt_chemical_dnn_with <- m3(df_fish_dnn, "chemicalFold")
rsqrt_species_dnn_with <- m3(df_fish_dnn, "speciesFold")

results_r_squared_dnn_with <- as.data.frame(rbind(rsqrt_random_dnn_with, rsqrt_chemical_dnn_with, rsqrt_species_dnn_with))

saveRDS(results_r_squared_dnn_with, file = "/home/isabellehalbhuber/Toxicology/Results/results_r_squared_dnn_with_MD_noRegularization.rds")


#####################################################################
######################## DNN without ################################
#####################################################################
# df und calculate_r_squared function von oben

m4 <- function(df, fold_column, k = 10, epochs = 500) {

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


