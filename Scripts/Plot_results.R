#### Plot results ####
library(readr)
#### Random Forest ####
RF_noTraits_Fingerprints <- readRDS("~/paper/EcoToxTM/Results/RF_noTraits_Fingerprints.rds")
RF_noTraits_MolDescriptor <- readRDS("~/paper/EcoToxTM/Results/RF_noTraits_MolDescriptor.rds")
RF_Traits_Fingerprints <- readRDS("~/paper/EcoToxTM/Results/RF_Traits_Fingerprints.rds")
RF_Traits_MolDescriptor <- readRDS("~/paper/EcoToxTM/Results/RF_Traits_MolDescriptor.rds")

#### MLP ####
MLP_noTraits_Fingerprints <- readRDS("~/paper/EcoToxTM/Results/MLP_noTraits_Fingerprints.rds")
MLP_noTraits_MolDescriptor <- readRDS("~/paper/EcoToxTM/Results/MLP_noTraits_MolDescriptor.rds")
MLP_Traits_Fingerprints <- readRDS("~/paper/EcoToxTM/Results/MLP_Traits_Fingerprints.rds")
MLP_Traits_MolDescriptor <- readRDS("~/paper/EcoToxTM/Results/MLP_Traits_MolDescriptor.rds")

#### Bi-GRU ####
BiGRU_noTraits <- readRDS("~/paper/EcoToxTM/Results/BiGRU_noTraits.rds")
