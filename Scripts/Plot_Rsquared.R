### make plots for performance comparison (r-squared)
# without traits and with traits 
library(readr)

#### MACCS 

# RF with
RF_with_MACCS <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_rf_with_MACCS.rds")
RF_with_MACCS <- (RF_with_MACCS)[1:2]
RF_with_MACCS$Model = "RF Fingerprints"
RF_with_MACCS$Descriptor = "MACCS"
colnames(RF_with_MACCS)[1:2] <- c("train", "test")
rownames(RF_with_MACCS)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
RF_with_MACCS$CV = rownames(RF_with_MACCS)[1:3] 
rownames(RF_with_MACCS) <- NULL

# RF without 
RF_without_MACCS <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_rf_without_MACCS.rds")
RF_without_MACCS <- (RF_without_MACCS)[1:2]
RF_without_MACCS$Model = "RF Fingerprints"
RF_without_MACCS$Descriptor = "MACCS"
colnames(RF_without_MACCS)[1:2] <- c("train", "test")
rownames(RF_without_MACCS)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
RF_without_MACCS$CV = rownames(RF_without_MACCS)[1:3] 
rownames(RF_without_MACCS) <- NULL

# DNN 
DNN_with_MACCS <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_dnn_with_MACCS_noRegularization.rds")
DNN_with_MACCS$Model = "DNN Fingerprints"
DNN_with_MACCS$Descriptor = "MACCS"
colnames(DNN_with_MACCS)[1:2] <- c("train", "test")
rownames(DNN_with_MACCS)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
DNN_with_MACCS$CV = rownames(DNN_with_MACCS)[1:3] 
rownames(DNN_with_MACCS) <- NULL

DNN_without_MACCS <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_dnn_without_MACCS.rds")
DNN_without_MACCS$Model = "DNN Fingerprints"
DNN_without_MACCS$Descriptor = "MACCS"
colnames(DNN_without_MACCS)[1:2] <- c("train", "test")
rownames(DNN_without_MACCS)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
DNN_without_MACCS$CV = rownames(DNN_without_MACCS)[1:3] 
rownames(DNN_without_MACCS) <- NULL



####################
#### PubChem

# RF with
RF_with_Pubchem <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_rf_with_Pubchem.rds")
RF_with_Pubchem <- (RF_with_Pubchem)[1:2]
RF_with_Pubchem$Model = "RF"
RF_with_Pubchem$Descriptor = "PubChem"
colnames(RF_with_Pubchem)[1:2] <- c("train", "test")
rownames(RF_with_Pubchem)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
RF_with_Pubchem$CV = rownames(RF_with_Pubchem)[1:3] 
rownames(RF_with_Pubchem) <- NULL

# RF without 
RF_without_Pubchem <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_rf_without_MACCS.rds")
RF_without_Pubchem <- (RF_without_Pubchem)[1:2]
RF_without_Pubchem$Model = "RF"
RF_without_Pubchem$Descriptor = "PubChem"
colnames(RF_without_Pubchem)[1:2] <- c("train", "test")
rownames(RF_without_Pubchem)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
RF_without_Pubchem$CV = rownames(RF_without_Pubchem)[1:3] 
rownames(RF_without_Pubchem) <- NULL

# DNN
DNN_with_Pubchem <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_dnn_with_MACCS.rds")
DNN_with_Pubchem$Model = "DNN"
DNN_with_Pubchem$Descriptor = "PubChem"
colnames(DNN_with_Pubchem)[1:2] <- c("train", "test")
rownames(DNN_with_Pubchem)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
DNN_with_Pubchem$CV = rownames(DNN_with_Pubchem)[1:3] 
rownames(DNN_with_Pubchem) <- NULL

DNN_without_Pubchem <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_dnn_without_MACCS.rds")
DNN_without_Pubchem$Model = "DNN"
DNN_without_Pubchem$Descriptor = "PubChem"
colnames(DNN_without_Pubchem)[1:2] <- c("train", "test")
rownames(DNN_without_Pubchem)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
DNN_without_Pubchem$CV = rownames(DNN_without_Pubchem)[1:3] 
rownames(DNN_without_Pubchem) <- NULL

####################
#### MDs

# RF with
RF_with_MD <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_rf_with_MD_207.rds")
RF_with_MD <- (RF_with_MD)[1:2]
RF_with_MD$Model = "RF Mol.Descriptor"
RF_with_MD$Descriptor = "MolecularDescriptor"
colnames(RF_with_MD)[1:2] <- c("train", "test")
rownames(RF_with_MD)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
RF_with_MD$CV = rownames(RF_with_MD)[1:3] 
rownames(RF_with_MD) <- NULL

# RF without 
RF_without_MD <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_rf_without_207MD.rds")
RF_without_MD <- (RF_without_MD)[1:2]
RF_without_MD$Model = "RF Mol.Descriptor"
RF_without_MD$Descriptor = "MolecularDescriptor"
colnames(RF_without_MD)[1:2] <- c("train", "test")
rownames(RF_without_MD)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
RF_without_MD$CV = rownames(RF_without_MD)[1:3] 
rownames(RF_without_MD) <- NULL

# DNN
DNN_with_MD <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_dnn_with_MD_207.rds")
DNN_with_MD$Model = "DNN Mol.Descriptor"
DNN_with_MD$Descriptor = "MolecularDescriptor"
colnames(DNN_with_MD)[1:2] <- c("train", "test")
rownames(DNN_with_MD)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
DNN_with_MD$CV = rownames(DNN_with_MD)[1:3] 
rownames(DNN_with_MD) <- NULL

DNN_without_MD <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_dnn_without_MD_207.rds")
DNN_without_MD$Model = "DNN Mol.Descriptor"
DNN_without_MD$Descriptor = "MolecularDescriptor"
colnames(DNN_without_MD)[1:2] <- c("train", "test")
rownames(DNN_without_MD)[1:3] <- c("Random", "Chemical_name_blocked", "Species_name_blocked")
DNN_without_MD$CV = rownames(DNN_without_MD)[1:3] 
rownames(DNN_without_MD) <- NULL




#### SMILES 

# RNN with 
RNN_with<- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_rnn_with.rds")
RNN_with <- as.data.frame(t(RNN_with))
RNN_with$Model = "RNN SMILES"
RNN_with$Descriptor = "SMILES"
colnames(RNN_with)[1:2] <- c("train", "test")
rownames(RNN_with)[1:3] <- c("Random", "Species_name_blocked", "Chemical_name_blocked")
RNN_with$CV = rownames(RNN_with)[1:3] 
rownames(RNN_with) <- NULL

RNN_without <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_rnn_without.rds")
RNN_without <- as.data.frame(t(RNN_without))
RNN_without$Model = "RNN SMILES"
RNN_without$Descriptor = "SMILES"
colnames(RNN_without)[1:2] <- c("train", "test")
rownames(RNN_without)[1:3] <- c("Random", "Species_name_blocked", "Chemical_name_blocked")
RNN_without$CV = rownames(RNN_without)[1:3] 
rownames(RNN_without) <- NULL

# Attention
Attention_with<- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_Attention_with_test_1.rds")
Attention_with$Model = "Self-attention SMILES"
Attention_with$Descriptor = "SMILES"
colnames(Attention_with)[1:2] <- c("train", "test")
rownames(Attention_with)[1:3] <- c("Random", "Species_name_blocked", "Chemical_name_blocked")
Attention_with$CV= rownames(Attention_with)[1:3] 
rownames(Attention_with) <- NULL

Attention_without <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_Attention_without.rds")
Attention_without$Model = "Self-attention SMILES"
Attention_without$Descriptor = "SMILES"
colnames(Attention_without)[1:2] <- c("train", "test")
rownames(Attention_without)[1:3] <- c("Random", "Species_name_blocked", "Chemical_name_blocked")
Attention_without$CV = rownames(Attention_without)[1:3] 
rownames(Attention_without) <- NULL




#### combine all 


combinedRsquared_with <- as.data.frame(rbind(RF_with_MACCS, DNN_with_MACCS, RF_with_MD, DNN_with_MD, RNN_with, Attention_with))
combinedRsquared_with$test <- as.numeric(combinedRsquared_with$test)
combinedRsquared_without <- rbind(RF_without_MACCS, DNN_without_MACCS, RF_without_MD, DNN_without_MD, RNN_without, Attention_without)
combinedRsquared_without$test <- as.numeric(combinedRsquared_without$test)

combinedRsquared_with$CV <- factor(combinedRsquared_with$CV, levels = unique(combinedRsquared_with$CV))
combinedRsquared_with$Model <- factor(combinedRsquared_with$Model, levels = unique(combinedRsquared_with$Model))
combinedRsquared_with$Descriptor <- factor(combinedRsquared_with$Descriptor, levels = unique(combinedRsquared_with$Descriptor))
combinedRsquared_without$CV <- factor(combinedRsquared_without$CV, levels = unique(combinedRsquared_without$CV))
combinedRsquared_without$Model <- factor(combinedRsquared_without$Model, levels = unique(combinedRsquared_without$Model))

saveRDS(combinedRsquared_with, file = "/home/isabellehalbhuber/Toxicology/Results/combinedRsquared_with.rds")
saveRDS(combinedRsquared_without, file = "/home/isabellehalbhuber/Toxicology/Results/combinedRsquared_without.rds")

combinedRsquared_with <- readRDS("/home/isabellehalbhuber/Toxicology/Results/combinedRsquared_with.rds")
combinedRsquared_with <- as.data.frame(combinedRsquared_with)
# class(combinedRsquared_with)
# sapply(combinedRsquared_with, class)
combinedRsquared_with$train <- unlist(combinedRsquared_with$train)
write.csv(combinedRsquared_with, file = "/home/isabellehalbhuber/Toxicology/Results/combinedRsquared_with.csv")

combinedRsquared_without <- readRDS("/home/isabellehalbhuber/Toxicology/Results/combinedRsquared_without.rds")
combinedRsquared_without <- as.data.frame(combinedRsquared_without)
combinedRsquared_without$train <- unlist(combinedRsquared_without$train)
write.csv(combinedRsquared_without, file = "/home/isabellehalbhuber/Toxicology/Results/combinedRsquared_without.csv")


colors <- c("#f8dac0","#f8dac0", "#cfebce","#cfebce","lavender", "lavender","lightyellow","lightyellow", "thistle","thistle","lightblue","lightblue")
colors_legend <- c("#f8dac0", "#cfebce","lavender","lightyellow","thistle" ,"lightblue")

combinedRsquared_without_random <- combinedRsquared_without[combinedRsquared_without$CV == "Random",]
combinedRsquared_without_chemical <- combinedRsquared_without[combinedRsquared_without$CV == "Chemical_name_blocked",]
combinedRsquared_without_species <- combinedRsquared_without[combinedRsquared_without$CV == "Species_name_blocked",]



combinedRsquared_with_random <- combinedRsquared_with[combinedRsquared_with$CV == "Random",]
combinedRsquared_with_chemical <- combinedRsquared_with[combinedRsquared_with$CV == "Chemical_name_blocked",]
combinedRsquared_with_species <- combinedRsquared_with[combinedRsquared_with$CV == "Species_name_blocked",]

combinedRsquared_without_random$trait <- "without"
combinedRsquared_with_random$trait <- "with"
with_without_random <- rbind(combinedRsquared_without_random, combinedRsquared_with_random)
sorted_data_random <- with_without_random[order(with_without_random$Model, with_without_random$Descriptor), ]


combinedRsquared_without_chemical$trait <- "without"
combinedRsquared_with_chemical$trait <- "with"
with_without_chemical <- rbind(combinedRsquared_without_chemical, combinedRsquared_with_chemical)
sorted_data_chemical <- with_without_chemical[order(with_without_chemical$Model, with_without_chemical$Descriptor), ]

combinedRsquared_without_species$trait <- "without"
combinedRsquared_with_species$trait <- "with"
with_without_species <- rbind(combinedRsquared_without_species, combinedRsquared_with_species)
sorted_data_species <- with_without_species[order(with_without_species$Model, with_without_species$Descriptor), ]



#density <- rep(c(NA, 90), times = ncol(with_without_species))  # Jeder zweite Balken schraffiert
#angle <- rep(c(NA, 45), times = ncol(with_without_species))    # Schraffierung im Winkel von 45°




pdf("barplot_sorted_random_new.pdf", width = 10, height = 8)
par(lwd = 5)
p <- barplot(
  height = matrix(sorted_data_random$test, byrow = TRUE),
  beside = TRUE,
  #names.arg = rev(c("Random CV", "Chemical blocked CV", "Species blocked CV")),
  ylim = c(0, 1),
  #horiz = TRUE,                
  ylab = expression(R^2), 
  xlab = "Models",
  main = expression(R^2 ~ "for the QSAR-ML models with species traits"),
  col = colors,
  #density = density,
  #angle = angle,
  cex.main = 2.5,
  cex.axis = 1.5 ,
  #space = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,1,0.2,0.2,0.2)
  space = c(0.15,1) 
)

dev.off()
legend(
  "topright",
  legend = levels(combinedRsquared_with$Model),
  fill = colors_legend, 
  title = "Models & Chemical Traits",
  xpd = TRUE,
  horiz = FALSE,
  bty = "n",
  cex = 0.75,
)
# text(
#   x = p,                                          
#   y = matrix(combinedRsquared_without$test, ncol = 3, byrow = TRUE) + 0.01, 
#   labels = round(matrix(combinedRsquared_without$test, ncol = 3, byrow = TRUE), 2), 
#   pos = 3                              
# )


dev.off()

