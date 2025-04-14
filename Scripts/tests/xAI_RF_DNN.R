# 1. Variable Importance for RF
# 2. ALE Plots for DNN
# 3. ALE Plots for RF

###########################################################################
######################## Variable importance RF ###########################
###########################################################################
library(readr)
library(dplyr)
library(ranger)
library(cito)
library(ggplot2)
library(iml)

### exclude smiles codes
df_fish_rf <- readRDS("/home/isabellehalbhuber/Toxicology/Data/df_fish_rf.rds")
df_fish_rf <- df_fish_rf %>% select(-SMILES, -QSAR_READY_SMILES)

#############################################
### Variable importance for RF without traits 
#############################################


data = df_fish_rf %>% select(-randomFold, -chemicalFold, -speciesFold, -BEl, -VEp, -REs, -OGp, -BLs, -PFs, -PFv, -CPt, -RMl, -mem_1, -mem_2, -mem_3, -mem_4, -mem_5, -mem_6, -mem_7, -body_size)


rf_without <- ranger(logmeanConc ~ ., 
                     data = data,
                     min.node.size = 5, 
                     mtry = 884, 
                     max.depth = 15,
                     importance = "impurity")


# Variable importance
imp_values <- as.numeric(rf_without$variable.importance)
# scale values, weil im ranger paket der gini-index berechnet wird, und dieser keine prozentangabe ich sondern die Zahl, wie viel unreinheit beseitigt wird 
scaled_importance_values <- (imp_values - min(imp_values)) / (max(imp_values) - min(imp_values))
#filtered_importances_1 <- sorted_imp[sorted_imp$importance_values >= 1 & sorted_imp$importance_values <= 100, ]

# plot it 
importance_values <- data.frame(
  Variable = c(colnames(data[, -4])),
  Importance = c(scaled_importance_values) 
)

importance_values <- importance_values[order(-importance_values$Importance), ]

png("/home/isabellehalbhuber/Toxicology/Results/feature_importance_rf_without.png", width = 800, height = 600)

barplot(importance_values$Importance[1:30], 
        names.arg = importance_values$Variable[1:30],
        las = 2, 
        main = "Feature Importance",  
        xlab = "Variable",
        ylab = "Importance",
        cex.names = 0.8,  
        ylim = c(0, max(importance_values$Importance) * 1.1))
dev.off()

##########################################
### Variable importance for RF with traits 
##########################################

data = df_fish_rf %>% select(-randomFold, -chemicalFold, -speciesFold)
rf_with <- ranger(logmeanConc ~ ., 
                  data = data,
                  min.node.size = 100, 
                  mtry = 901, 
                  max.depth = 15,
                  importance = "impurity")

# Variable importance
imp_values <- as.numeric(rf_with$variable.importance)
# scale values, weil im ranger paket der gini-index berechnet wird, und dieser keine prozentangabe ich sondern die Zahl, wie viel unreinheit beseitigt wird 
scaled_importance_values <- (imp_values - min(imp_values)) / (max(imp_values) - min(imp_values))
#filtered_importances_1 <- sorted_imp[sorted_imp$importance_values >= 1 & sorted_imp$importance_values <= 100, ]

# plot it 
importance_values <- data.frame(
  Variable = c(colnames(data[, -21])),
  Importance = c(scaled_importance_values) 
)

importance_values <- importance_values[order(-importance_values$Importance), ]

png("/home/isabellehalbhuber/Toxicology/Results/feature_importance_rf_with.png", width = 800, height = 600)

barplot(importance_values$Importance[1:30], 
        names.arg = importance_values$Variable[1:30],
        las = 2, 
        main = "Feature Importance",  
        xlab = "Variable",
        ylab = "Importance",
        cex.names = 0.8,  
        ylim = c(0, max(importance_values$Importance) * 1.1))
dev.off()






###################################################
### Friednmans H-statistic for RF using imp package 
###################################################

data = df_fish_rf %>% select(-randomFold, -chemicalFold, -speciesFold)

library(iml)

predict_wrapper = function(model, newdata) predict(model, data=newdata)$predictions
predictor = Predictor$new(model = rf_with,
                          data = data[, -21],
                          y = data$logmeanConc,
                          predict.function = predict_wrapper)

predictor$task = "regression" # set task to regression

# Friedman's H-statistic
interact = Interaction$new(predictor, grid.size = 5L)
plot(interact)
results_interaction = as.data.frame(interact$results)
results_interaction <- results_interaction[order(-results_interaction$.interaction), ]

barplot(results_interaction$.interaction[1:50],
        names.arg = results_interaction$.feature[1:50],
        las = 2,
        main = "Overall interaction",
        xlab = "Variable",
        ylab = "Interaction",
        cex.names = 0.8,
        ylim = c(0, max(results_interaction$.interaction) * 1.1))


########################################################################################
################################# ALE plots for the DNN ################################
########################################################################################


df_fish_rf <- readRDS("/home/isabellehalbhuber/Toxicology/Data/df_fish_rf.rds")
df_fish_rf <- df_fish_rf %>% select(-SMILES, -QSAR_READY_SMILES)
# scale data for DNN 
variables_to_scale <- c("ALogP", "ALogp2", "AMR", "BEl", "VEp", "REs", "OGp", 
                        "RMl", "BLs", "PFv", "PFs", "CPt", "body_size", 
                        "mem_1", "mem_2", "mem_3", "mem_4", "mem_5", "mem_6", "mem_7")

scale <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

df_fish_rf <- df_fish_rf %>%
  mutate(across(all_of(variables_to_scale), scale))

data_with = df_fish_rf %>% select(-randomFold, -chemicalFold, -speciesFold)
data_without = df_fish_rf %>% select(-randomFold, -chemicalFold, -speciesFold,  -BEl, -VEp, -REs, -OGp, -BLs, -PFs, -PFv, -CPt, -RMl, -mem_1, -mem_2, -mem_3, -mem_4, -mem_5, -mem_6, -mem_7, -body_size)


# dnn_model_with <- dnn(logmeanConc ~ ., 
#                       hidden = c(20L, 20L),
#                       batchsize = 300L,
#                       data = data_with,
#                       epochs = 500,
#                       loss = 'mae',
#                       burnin = 200,
#                       activation = "relu",
#                       #dropout = 0.1,
#                       lr = 0.001,     
#                       #lambda = 0.1,
#                       #alpha = 0.01,
#                       #lr_scheduler = config_lr_scheduler('step', gamma = 0.7, step_size = 10),
#                       device = "cuda")
# 
# 
# # without species traits
# 
# dnn_model_without <- dnn(logmeanConc ~ ., 
#                          hidden = c(20L, 20L),
#                          batchsize = 300L,
#                          data = data_without,
#                          epochs = 500,
#                          loss = 'mae',
#                          burnin = 200,
#                          activation = "relu",
#                          #dropout = 0.1,
#                          lr = 0.001,     
#                          #lambda = 0.1,
#                          #alpha = 0.01,
#                          #lr_scheduler = config_lr_scheduler('step', gamma = 0.7, step_size = 10),
#                          device = "cuda")
# 
# 
# saveRDS(dnn_model_without, file = "/home/isabellehalbhuber/Toxicology/Scripts/xAI_dnn_model_without.rds") 
# saveRDS(dnn_model_with, file = "/home/isabellehalbhuber/Toxicology/Scripts/xAI_dnn_model_with.rds") 


dnn_model_without <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/xAI_dnn_model_without.rds")
dnn_model_with <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/xAI_dnn_model_with.rds")


plot_ale <- function(model_with, model_without, data_with, data_without, variables, output_pdf) {
  
  pdf(output_pdf, width = 12, height = 12)
  par(mfrow = c(3, 3))  
  
  for (variable in variables) {
    ale_with_DNN <- tryCatch({
      ALE(model_with, data = as.data.frame(data_with), plot = FALSE, variable = variable)[[1]]$data
    }, error = function(e) {
      warning(paste("Error", variable, "for", e$message))
      return(NULL)
    })
    
    ale_without_DNN <- tryCatch({
      ALE(model_without, data = as.data.frame(data_without), plot = FALSE, variable = variable)[[1]]$data
    }, error = function(e) {
      warning(paste("Error", variable, "for", e$message))
      return(NULL)
    })
    
    if (is.null(ale_with_DNN) && is.null(ale_without_DNN)) {
      next
    }
    
    if (!is.null(ale_with_DNN)) {
      ale_with_DNN$Type <- "with"
    }
    if (!is.null(ale_without_DNN)) {
      ale_without_DNN$Type <- "without"
    }
    
    if (!is.null(ale_with_DNN) && !is.null(ale_without_DNN)) {
      ale_combined <- rbind(ale_with_DNN, ale_without_DNN)
    } else if (!is.null(ale_with_DNN)) {
      ale_combined <- ale_with_DNN
    } else {
      ale_combined <- ale_without_DNN
    }
    
    types <- unique(ale_combined$Type)
    
    plot(NULL, xlim = range(ale_combined$x), ylim = range(ale_combined$y),
         xlab = variable, ylab = "ALE", 
         main = paste("ALE for", variable), type = "n")
    
    for (type in types) {
      subset_data <- subset(ale_combined, Type == type)
      lines(subset_data$x, subset_data$y, col = ifelse(type == "with", "orange", "lightblue"), lwd = 2)
    }
    # for species traits and therefore no "dnn_model_without"
    if (length(types) == 1) {
      if (types == "with") {
        legend("topright", legend = "with", col = "orange", lwd = 2, title = "Type")
      }
    } else {
      legend("topright", legend = types, col = c("orange", "lightblue"), lwd = 2, title = "Type")
    }
  }
  

  dev.off()
}


variables <- c("ALogP", "AMR", "bit435", "bit17", "ALogp2", "bit664", "PFv", "REs", "body_size")
output_pdf <- "/home/isabellehalbhuber/Toxicology/Results/DNN_ALE_plots.pdf"  # Pfad zur Ausgabedatei

plot_ale(dnn_model_with, dnn_model_without, data_with, data_without, variables, output_pdf)



########################################################################################
################################# ALE plots for the RF #################################
########################################################################################

df_fish_rf <- readRDS("/home/isabellehalbhuber/Toxicology/Data/df_fish_rf.rds")
df_fish_rf <- df_fish_rf %>% select(-SMILES, -QSAR_READY_SMILES)
data_with = df_fish_rf %>% select(-randomFold, -chemicalFold, -speciesFold)
data_without = df_fish_rf %>% select(-randomFold, -chemicalFold, -speciesFold,  -BEl, -VEp, -REs, -OGp, -BLs, -PFs, -PFv, -CPt, -RMl, -mem_1, -mem_2, -mem_3, -mem_4, -mem_5, -mem_6, -mem_7, -body_size)

# rf_without <- ranger(logmeanConc ~ ., 
#                      data = data_without,
#                      min.node.size = 100, 
#                      mtry = 884, 
#                      max.depth = 15,
#                      importance = "impurity")
# 
# rf_with <- ranger(logmeanConc ~ ., 
#                   data = data_with,
#                   min.node.size = 100, 
#                   mtry = 901, 
#                   max.depth = 15,
#                   importance = "impurity")

# save the models 
# saveRDS(rf_with, file = "/home/isabellehalbhuber/Toxicology/Scripts/xAI_rf_with.rds") 
# saveRDS(rf_without, file = "/home/isabellehalbhuber/Toxicology/Scripts/xAI_rf_without.rds") 

rf_without <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/xAI_rf_without.rds")
rf_with <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/xAI_rf_with.rds")




generate_ale_plots <- function(rf_with, rf_without, data_with, data_without, output_pdf = "RF_ALE_Plots.pdf") {

  pdf(output_pdf, width = 12, height = 12)  
  par(mfrow = c(3, 3))  

  predict_wrapper <- function(model, newdata) predict(model, data = newdata)$predictions
  
  predictor_with <- Predictor$new(model = rf_with,
                                  data = data_with[, -21],
                                  y = data_with$logmeanConc,
                                  predict.function = predict_wrapper)
  predictor_with$task <- "regression"  
  
  predictor_without <- Predictor$new(model = rf_without,
                                     data = data_without[, -4],
                                     y = data_without$logmeanConc,
                                     predict.function = predict_wrapper)
  predictor_without$task <- "regression" 
  

  plot_ale_feature <- function(feature_name, predictor_with, predictor_without, xlab, ylab, main) {
    ale_with <- FeatureEffect$new(predictor_with, feature = feature_name, method = "ale")$results
    ale_without <- FeatureEffect$new(predictor_without, feature = feature_name, method = "ale")$results
    
    ale_with$Type <- "with"
    ale_without$Type <- "without"
    
    colnames(ale_with)[3] <- "x"
    colnames(ale_with)[2] <- "y"
    colnames(ale_without)[3] <- "x"
    colnames(ale_without)[2] <- "y"
    
    ale_combined <- rbind(ale_with, ale_without)
    types <- unique(ale_combined$Type)
    
    plot(NULL, xlim = range(ale_combined$x), ylim = range(ale_combined$y),
         xlab = xlab, ylab = ylab,
         main = main, type = "n")
    
    for (type in types) {
      subset_data <- subset(ale_combined, Type == type)
      lines(subset_data$x, subset_data$y, col = ifelse(type == "with", "orange", "lightblue"), lwd = 2)
    }
    legend("topright", legend = types, col = c("orange", "lightblue"), lwd = 2, title = "Type")
  }
  

  features <- c("ALogP", "AMR", "bit435", "bit17", "ALogp2", "bit664")
  

  for (feature in features) {
    plot_ale_feature(feature, predictor_with, predictor_without, xlab = feature, ylab = "ALE", main = paste("ALE for", feature))
  }
  
# plot for species traits (no plot for "rf_without")
  plot_single_ale <- function(feature_name, predictor, xlab, ylab, main) {
    ale_result <- FeatureEffect$new(predictor, feature = feature_name, method = "ale")$results
    plot(NULL, xlim = range(ale_result[[feature_name]]), ylim = range(ale_result$.value),
         xlab = xlab, ylab = ylab,
         main = main, type = "n")
    lines(ale_result[[feature_name]], ale_result$.value, col = "orange", lwd = 2, type = "l")
  }
  

  plot_single_ale("PFv", predictor_with, "PFv", "ALE", "ALE for PFv")
  plot_single_ale("REs", predictor_with, "REs", "ALE", "ALE for REs")
  plot_single_ale("body_size", predictor_with, "body_size", "ALE", "ALE for body_size")
  

  dev.off()
}


generate_ale_plots(rf_with, rf_without, data_with, data_without, output_pdf = "/home/isabellehalbhuber/Toxicology/Results/RF_ALE_plots.pdf")



