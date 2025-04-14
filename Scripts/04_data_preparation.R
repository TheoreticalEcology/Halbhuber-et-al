################################################
############ Combine all DFs ###################
################################################
library(readr)
library(dplyr)

df_fish_traits <- readRDS("Data/fish_traits_dataframe.rds")
# filter for relevant trait information 
df_fish_traits <- df_fish_traits[,-c(1,2,3,4,5,6,17,18, 27)]
df_fish_traits <- subset(df_fish_traits, select=c(11, 1:10,12:18))
df_toxicity <- readRDS(file = "Data/LC50_dataframe.rds")
df_toxicity <- df_toxicity[,-c(2,4,5,6)]
colnames(df_toxicity)[1] <- "CAS"
colnames(df_toxicity)[2] <- "scientificNameStd"
colnames(df_toxicity)[3] <- "Concentration"



## Data for sequence-based DL model: SMILES, LC50 and species traits 
df_SMILES<- readRDS(file="Data/chemical_traits_dataframe.rds")
df_SMILES$CAS <- as.numeric(df_SMILES$CAS)  # NAs for "NoCAS" values
df_SMILES <- df_SMILES[, -c(1,2,4)]
# merge SMILES and LC50
SMILES_LC50 <- df_toxicity  |>
  full_join(df_SMILES, by = c( CAS = "CAS"))

SMILES_LC50 <- na.omit(SMILES_LC50)
nrow(unique(SMILES_LC50[,c("SMILES","scientificNameStd")])) 

#merge SMILES_LC50 and traits
SMILES_LC50_traits <- SMILES_LC50  |>
  full_join(df_fish_traits, by = c( scientificNameStd = "scientificNameStd"))
SMILES_LC50_traits <- na.omit(SMILES_LC50_traits)
nrow(unique(SMILES_LC50_traits[,c("SMILES","scientificNameStd")])) 

saveRDS(SMILES_LC50_traits,file= "Data/SMILES_LC50_traits.rds")


## Data for tabular-data based DL model: FPs, MDs, LC50 and species traits 

######################################  MDs
MDs <- readRDS("Data/molecularDescriptors_CAS.rds")
MDs <- MDs[,c(5:211)]
MDs$CAS <- as.numeric(MDs$CAS)

# merge MDs and LC50
df_MD_LC50 <- df_toxicity  |>
  full_join(MDs, by = c( CAS = "CAS"))

df_MD_LC50 <- na.omit(df_MD_LC50)

# merge df_MD_LC50 and traits
df_MD_LC50_traits <- df_MD_LC50  |>
  full_join(df_fish_traits, by = c( scientificNameStd = "scientificNameStd"))
df_MD_LC50_traits <- na.omit(df_MD_LC50_traits)
nrow(unique(df_MD_LC50_traits[,c("CAS","scientificNameStd")])) 

df_MD_LC50_traits <- df_MD_LC50_traits[,-c(2,4,5,6)]
saveRDS(SMILES_LC50_traits,file= "Data/SMILES_LC50_traits.rds")

###################################### FPs
MACCS <- readRDS("Data/CAS_SMILES_MACCS.rds")
MACCS <- MACCS[,c('CAS', 'Fingerprint')]
MACCS$CAS <- as.numeric(MACCS$CAS)

# merge MACCS and LC50
df_MACCS_LC50 <- df_toxicity  |>
  full_join(MACCS, by = c( CAS = "CAS"))

df_MACCS_LC50 <- na.omit(df_MACCS_LC50)

nrow(unique(df_MACCS_LC50[,c("CAS","scientificNameStd")]))
# merge df_MACCS_LC50 and traits 
df_MACCS_LC50_traits <- df_MACCS_LC50  |>
  full_join(df_fish_traits, by = c( scientificNameStd = "scientificNameStd"))
df_MACCS_LC50_traits <- na.omit(df_MACCS_LC50_traits)
nrow(unique(df_MACCS_LC50_traits[,c("CAS","scientificNameStd")])) 

saveRDS(df_MACCS_LC50_traits,file= "Data/df_MACCS_LC50_traits")

###################################### Full Dataset
# CHEMICAL TRAITS
MACCS_SMILES <- readRDS("Data/CAS_SMILES_MACCS.rds")
MACCS_SMILES <- MACCS_SMILES[,c('CAS', 'SMILES', 'Fingerprint')]

MDs <- readRDS("Data/molecularDescriptors_CAS.rds")
MDs <- MDs[,c(5:211)]

chemical_traits <- MACCS_SMILES |>
  full_join(MDs, by = c(CAS = "CAS"))

chemical_traits$CAS <- as.numeric(chemical_traits$CAS)

# join LC50 values and chemical traits
df_all_chems_LC50 = df_toxicity |>
  full_join(chemical_traits, by = c(CAS = "CAS"))

df_all_chems_LC50 <- na.omit(df_all_chems_LC50)

# merge df_all_chems_LC50 with traits

df_complete = df_all_chems_LC50 |>
  full_join(df_fish_traits, by = c(scientificNameStd = "scientificNameStd"))

colnames(df_complete)
df_complete <- na.omit(df_complete)

nrow(unique(df_complete[,c("CAS","scientificNameStd")])) 
length(unique(df_complete$CAS)) 

saveRDS(df_complete, file = "Data/df_complete.rds")



