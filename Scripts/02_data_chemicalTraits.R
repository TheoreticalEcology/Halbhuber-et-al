#########################
###  Chemical Trait  #### 
#########################

### SMILES Codes ###


# Prerequisite: download data from source: https://clowder.edap-cluster.com/files/6616d8d7e4b063812d70fc95?dataset=61147fefe4b0856fdc65639b&space=&folder=6616d85ce4b063812d70fc8f
library(readr)
dsstox_data <- read_delim("Data/DSSToxDump1.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

# remove the '-' from the 'CASRN' values to make the values compatible with the CAS values from EPA-Ecotox toxicity data 
library(dplyr)
dsstox_data <- dsstox_data %>% mutate(CAS = gsub("-", "", CASRN))

# filter out relevant variables: "CASRN", "IUPAC_NAME", SMILES", "QSAR_READY_SMILES" "CAS"     
dsstox_data_filtered <- dsstox_data[ , c(3, 5, 6, 10,  13)]
# save df 
saveRDS(dsstox_data_filtered, file="Data/chemical_traits_dataframe.rds")

### SMILES ###

# how many unique SMILES-Species combinations do we have? 
df_chemicals<- readRDS(file="Data/chemical_traits_dataframe.rds")
df_toxicity <- readRDS(file = "Data/LC50_dataframe.rds")
colnames(df_toxicity)[1] <- "CAS"
df_chemicals$CAS <- as.numeric(df_chemicals$CAS)  # NAs for "NoCAS" values
library('dplyr')
df_merged <- df_toxicity  |>
  full_join(df_chemicals, by = c( CAS = "CAS"))
# filter for smiles, species and LC50
df_merged_SMILES_Species_LC50 <- df_merged[,c(3, 7, 10)]
df_merged_SMILES_Species_LC50 <- na.omit(df_merged_SMILES_Species_LC50)
length(unique(df_merged_SMILES_Species_LC50$SMILES))
length(unique(df_merged_SMILES_Species_LC50$`Species Scientific Name`))
nrow(unique(df_merged_SMILES_Species_LC50[,c("SMILES","Species Scientific Name")])) ### unique SMILES chemical combinations: 24744

saveRDS(df_merged_SMILES_Species_LC50, file = "Data/SMILES_Toxicity.rds")

### Exkurs, wie viele smiles und species traits? 
fish_traits <-readRDS("Data/only_fishTraits.rds")
colnames(df_merged_SMILES_Species_LC50)[1] <- "scientificNameStd"
df_merged <- df_merged_SMILES_Species_LC50  |>
  full_join(fish_traits, by = c( scientificNameStd = "scientificNameStd"))
df_merged <- na.omit(df_merged)
nrow(unique(df_merged[,c("SMILES","scientificNameStd")])) ### unique SMILES chemical combinations: 24744
length(unique(df_merged$SMILES))
length(unique(df_merged$scientificNameStd))


### Fingerprints ###
install.packages("rcdklibs")
install.packages("rcdk")
install.packages("rJava")

# Attention: Remove and install the packages to make them work
#remove.packages("rcdk") ## habe ich gemacht, weil fehler bei der "get.desc.names" funktion aufgetreten sind 
#remove.packages("rcdklibs")
#remove.packages("rJava")
#install.packages("rcdklibs")
#install.packages("rcdk")
#install.packages("rJava")
library(rJava)
library(rcdklibs)
library(rcdk)

# To save memory and time, we only calculate fingerprints for chemicals (SMILES codes) for which we have toxicity LC50 values
SMILES_EPA_DSS <- readRDS("Data/SMILES_Toxicity.rds")
SMILES_EPA_DSS <- c(unique(SMILES_EPA_DSS[, 3]))
# calculate three different fingerprints 
mols_EPA_DSS <- parse.smiles(SMILES_EPA_DSS$SMILES, kekulise = TRUE, omit.nulls = TRUE, smiles.parser = NULL)#  1 out of 3316 SMILES were not successfully parsed, resulting in NULLs
#fps_circular <- lapply(mols_EPA_DSS, get.fingerprint, type='circular') ## type = ECFP6
#fps_pubchem <- lapply(mols_EPA_DSS, get.fingerprint, type='pubchem')
fps_maccs <- lapply(mols_EPA_DSS, get.fingerprint, type='maccs')

# save fingerprints 
# saveRDS(fps_circular, file = "Data/circular_fingerprints_EPA_DSS.rds")
# saveRDS(fps_pubchem, file = "Data/fps_pubchem_EPA_DSS.rds")
# saveRDS(fps_maccs, file = "Data/fps_maccs_EPA_DSS.rds")

### Merge Fingerprints with SMILES codes dataset ###


chemicals <- readRDS("Data/chemical_traits_dataframe.rds")

# test a single SMILES code 
#character_fps_pubchem_1 <- as.character(fps_pubchem$`OC1=C(C=C(C=C1[N+]([O-])=O)[N+]([O-])=O)[N+]([O-])=O`)
#character_fps_pubchem_2 <- as.character(fps_pubchem[[1]])


# 
#character_fps_pubchem_all <- sapply(fps_pubchem, as.character)
character_fps_maccs_all <- sapply(fps_maccs, as.character)


# make dataframe for pubchem fingerprints 
#character_fps_pubchem_all <- data.frame(SMILES = names(character_fps_pubchem_all),Fingerprint = as.character(character_fps_pubchem_all),stringsAsFactors = FALSE)
#saveRDS(character_fps_pubchem_all, file = "Data/binary_fingerprints_EPA_DSS_pubchem.rds")

# make dataframe for MACCS fingerprints 

character_fps_maccs_all <- data.frame(SMILES = names(character_fps_maccs_all),Fingerprint = as.character(character_fps_maccs_all), stringsAsFactors = FALSE)

#saveRDS(character_fps_maccs_all, file = "Data/binary_fingerprints_EPA_DSS_maccs.rds")


# join pubchem fingerprints
#fps_pubchem_with_CAS = chemicals |>
#  full_join(character_fps_pubchem_all, by = c(SMILES = "SMILES"))

# join MACCS fingerprints 
fps_maccs_with_CAS = chemicals |>
  full_join(character_fps_maccs_all, by = c(SMILES = "SMILES"))

# we have decided to model MACCS fingerprints
saveRDS(fps_maccs_with_CAS, file = "Data/CAS_SMILES_MACCS.rds")



### Molecular Descriptors ###

# define molecular descriptors
descNames <- get.desc.names("all")

# test
# # 1. FractionalCSP3Descriptor
FractionalCSP3Descriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor", descNames)])
SMILES <- rownames(FractionalCSP3Descriptor_EPA_DSS)
rownames(FractionalCSP3Descriptor_EPA_DSS) <- NULL
FractionalCSP3Descriptor<- cbind(SMILES,FractionalCSP3Descriptor_EPA_DSS)

# # 2. SmallRingDescriptor
SmallRingDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.SmallRingDescriptor", descNames)])
SMILES <- rownames(SmallRingDescriptor_EPA_DSS)
rownames(SmallRingDescriptor_EPA_DSS) <- NULL
SmallRingDescriptor <- cbind(SMILES,SmallRingDescriptor_EPA_DSS)

# # 3. FractionalPSADescriptor
FractionalPSADescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.FractionalPSADescriptor", descNames)])
SMILES <- rownames(FractionalPSADescriptor_EPA_DSS)
rownames(FractionalPSADescriptor_EPA_DSS) <- NULL
FractionalPSADescriptor <- cbind(SMILES,FractionalPSADescriptor_EPA_DSS)

# # 4. ZagrebIndexDescriptor
ZagrebIndexDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor", descNames)])
SMILES <- rownames(ZagrebIndexDescriptor_EPA_DSS)
rownames(ZagrebIndexDescriptor_EPA_DSS) <- NULL
ZagrebIndexDescriptor <- cbind(SMILES,ZagrebIndexDescriptor_EPA_DSS)

# 5. XLogPDescriptor
XLogPDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor", descNames)])
SMILES <- rownames(XLogPDescriptor_EPA_DSS)
rownames(XLogPDescriptor_EPA_DSS) <- NULL
XLogPDescriptor <- cbind(SMILES,XLogPDescriptor_EPA_DSS)

# # 6. WienerNumbersDescriptor
WienerNumbersDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor", descNames)])
SMILES <- rownames(WienerNumbersDescriptor_EPA_DSS)
rownames(WienerNumbersDescriptor_EPA_DSS) <- NULL
WienerNumbersDescriptor <- cbind(SMILES,WienerNumbersDescriptor_EPA_DSS)

# # 9. WeightDescriptor
WeightDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor", descNames)])
SMILES <- rownames(WeightDescriptor_EPA_DSS)
rownames(WeightDescriptor_EPA_DSS) <- NULL
WeightDescriptor <- cbind(SMILES,WeightDescriptor_EPA_DSS)

# # 10. VAdjMaDescriptor
VAdjMaDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor", descNames)])
SMILES <- rownames(VAdjMaDescriptor_EPA_DSS)
rownames(VAdjMaDescriptor_EPA_DSS) <- NULL
VAdjMaDescriptor <- cbind(SMILES,VAdjMaDescriptor_EPA_DSS)

# # 27. HBondDonorCountDescriptor
HBondDonorCountDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor", descNames)])
SMILES <- rownames(HBondDonorCountDescriptor_EPA_DSS)
rownames(HBondDonorCountDescriptor_EPA_DSS) <- NULL
HBondDonorCountDescriptor<- cbind(SMILES,HBondDonorCountDescriptor_EPA_DSS)

# # 28. HBondAcceptorCountDescriptor
HBondAcceptorCountDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor", descNames)])
SMILES <- rownames(HBondAcceptorCountDescriptor_EPA_DSS)
rownames(HBondAcceptorCountDescriptor_EPA_DSS) <- NULL
HBondAcceptorCountDescriptor<- cbind(SMILES,HBondAcceptorCountDescriptor_EPA_DSS)

# # 53. AminoAcidCountDescriptor
AminoAcidCountDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor", descNames)])
SMILES <- rownames(AminoAcidCountDescriptor_EPA_DSS)
rownames(AminoAcidCountDescriptor_EPA_DSS) <- NULL
AminoAcidCountDescriptor_EPA_DSS<- cbind(SMILES,AminoAcidCountDescriptor_EPA_DSS)

# # 51. AcidicGroupCountDescriptor
AcidicGroupCountDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor", descNames)])
SMILES <- rownames(AcidicGroupCountDescriptor_EPA_DSS)
rownames(AcidicGroupCountDescriptor_EPA_DSS) <- NULL
AcidicGroupCountDescriptor_EPA_DSS<- cbind(SMILES,AcidicGroupCountDescriptor_EPA_DSS)

# # 46. AtomCountDescriptor
AtomCountDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor", descNames)])
SMILES <- rownames(AtomCountDescriptor_EPA_DSS)
rownames(AtomCountDescriptor_EPA_DSS) <- NULL
AtomCountDescriptor_EPA_DSS<- cbind(SMILES,AtomCountDescriptor_EPA_DSS)

# # 44. AutocorrelationDescriptorMass
AutocorrelationDescriptorMass_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass", descNames)])
SMILES <- rownames(AutocorrelationDescriptorMass_EPA_DSS)
rownames(AutocorrelationDescriptorMass_EPA_DSS) <- NULL
AutocorrelationDescriptorMass_EPA_DSS<- cbind(SMILES,AutocorrelationDescriptorMass_EPA_DSS)


# # 43. AutocorrelationDescriptorPolarizability
AutocorrelationDescriptorPolarizability_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability", descNames)])
SMILES <- rownames(AutocorrelationDescriptorPolarizability_EPA_DSS)
rownames(AutocorrelationDescriptorPolarizability_EPA_DSS) <- NULL
AutocorrelationDescriptorPolarizability_EPA_DSS<- cbind(SMILES,AutocorrelationDescriptorPolarizability_EPA_DSS)


# # 42. BasicGroupCountDescriptor
BasicGroupCountDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor", descNames)])
SMILES <- rownames(BasicGroupCountDescriptor_EPA_DSS)
rownames(BasicGroupCountDescriptor_EPA_DSS) <- NULL
BasicGroupCountDescriptor_EPA_DSS<- cbind(SMILES,BasicGroupCountDescriptor_EPA_DSS)


# # 40. BondCountDescriptor
BondCountDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor", descNames)])
SMILES <- rownames(BondCountDescriptor_EPA_DSS)
rownames(BondCountDescriptor_EPA_DSS) <- NULL
BondCountDescriptor_EPA_DSS<- cbind(SMILES,BondCountDescriptor_EPA_DSS)

# # 38. CarbonTypesDescriptor
CarbonTypesDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor", descNames)])
SMILES <- rownames(CarbonTypesDescriptor_EPA_DSS)
rownames(CarbonTypesDescriptor_EPA_DSS) <- NULL
CarbonTypesDescriptor_EPA_DSS<- cbind(SMILES,CarbonTypesDescriptor_EPA_DSS)

# # 37. ChiChainDescriptor
ChiChainDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor", descNames)])
SMILES <- rownames(ChiChainDescriptor_EPA_DSS)
rownames(ChiChainDescriptor_EPA_DSS) <- NULL
ChiChainDescriptor_EPA_DSS<- cbind(SMILES,ChiChainDescriptor_EPA_DSS)

# # 36. ChiClusterDescriptor
ChiClusterDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor", descNames)])
SMILES <- rownames(ChiClusterDescriptor_EPA_DSS)
rownames(ChiClusterDescriptor_EPA_DSS) <- NULL
ChiClusterDescriptor_EPA_DSS<- cbind(SMILES,ChiClusterDescriptor_EPA_DSS)

# # 35. ChiPathClusterDescriptor
ChiPathClusterDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor", descNames)])
SMILES <- rownames(ChiPathClusterDescriptor_EPA_DSS)
rownames(ChiPathClusterDescriptor_EPA_DSS) <- NULL
ChiPathClusterDescriptor_EPA_DSS<- cbind(SMILES,ChiPathClusterDescriptor_EPA_DSS)

# # 34. ChiPathDescriptor
ChiPathDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor", descNames)])
SMILES <- rownames(ChiPathDescriptor_EPA_DSS)
rownames(ChiPathDescriptor_EPA_DSS) <- NULL
ChiPathDescriptor_EPA_DSS<- cbind(SMILES,ChiPathDescriptor_EPA_DSS)


# # 32. EccentricConnectivityIndexDescriptor
EccentricConnectivityIndexDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor", descNames)])
SMILES <- rownames(EccentricConnectivityIndexDescriptor_EPA_DSS)
rownames(EccentricConnectivityIndexDescriptor_EPA_DSS) <- NULL
EccentricConnectivityIndexDescriptor_EPA_DSS<- cbind(SMILES,EccentricConnectivityIndexDescriptor_EPA_DSS)


# # 31. FMFDescriptor
FMFDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor", descNames)])
SMILES <- rownames(FMFDescriptor_EPA_DSS)
rownames(FMFDescriptor_EPA_DSS) <- NULL
FMFDescriptor_EPA_DSS<- cbind(SMILES,FMFDescriptor_EPA_DSS)


# # 30. FragmentComplexityDescriptor
FragmentComplexityDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor", descNames)])
SMILES <- rownames(FragmentComplexityDescriptor_EPA_DSS)
rownames(FragmentComplexityDescriptor_EPA_DSS) <- NULL
FragmentComplexityDescriptor_EPA_DSS<- cbind(SMILES,FragmentComplexityDescriptor_EPA_DSS)

# # 25. KappaShapeIndicesDescriptor
KappaShapeIndicesDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor", descNames)])
SMILES <- rownames(KappaShapeIndicesDescriptor_EPA_DSS)
rownames(KappaShapeIndicesDescriptor_EPA_DSS) <- NULL
KappaShapeIndicesDescriptor_EPA_DSS<- cbind(SMILES,KappaShapeIndicesDescriptor_EPA_DSS)

# # 24. KierHallSmartsDescriptor
KierHallSmartsDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor", descNames)])
SMILES <- rownames(KierHallSmartsDescriptor_EPA_DSS)
rownames(KierHallSmartsDescriptor_EPA_DSS) <- NULL
KierHallSmartsDescriptor_EPA_DSS<- cbind(SMILES,KierHallSmartsDescriptor_EPA_DSS)

# # 23. LargestChainDescriptor
LargestChainDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor", descNames)])
SMILES <- rownames(LargestChainDescriptor_EPA_DSS)
rownames(LargestChainDescriptor_EPA_DSS) <- NULL
LargestChainDescriptor_EPA_DSS<- cbind(SMILES,LargestChainDescriptor_EPA_DSS)

# # 22. LargestPiSystemDescriptor
LargestPiSystemDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor", descNames)])
SMILES <- rownames(LargestPiSystemDescriptor_EPA_DSS)
rownames(LargestPiSystemDescriptor_EPA_DSS) <- NULL
LargestPiSystemDescriptor_EPA_DSS<- cbind(SMILES,LargestPiSystemDescriptor_EPA_DSS)


# # 19. MannholdLogPDescriptor
MannholdLogPDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor", descNames)])
SMILES <- rownames(MannholdLogPDescriptor_EPA_DSS)
rownames(MannholdLogPDescriptor_EPA_DSS) <- NULL
MannholdLogPDescriptor_EPA_DSS<- cbind(SMILES,MannholdLogPDescriptor_EPA_DSS)


# # 18. MDEDescriptor
MDEDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor", descNames)])
SMILES <- rownames(MDEDescriptor_EPA_DSS)
rownames(MDEDescriptor_EPA_DSS) <- NULL
MDEDescriptor_EPA_DSS<- cbind(SMILES,MDEDescriptor_EPA_DSS)

# # 16. PetitjeanNumberDescriptor
PetitjeanNumberDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor", descNames)])
SMILES <- rownames(PetitjeanNumberDescriptor_EPA_DSS)
rownames(PetitjeanNumberDescriptor_EPA_DSS) <- NULL
PetitjeanNumberDescriptor_EPA_DSS<- cbind(SMILES,PetitjeanNumberDescriptor_EPA_DSS)

# # 14. RotatableBondsCountDescriptor
RotatableBondsCountDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor", descNames)])
SMILES <- rownames(RotatableBondsCountDescriptor_EPA_DSS)
rownames(RotatableBondsCountDescriptor_EPA_DSS) <- NULL
RotatableBondsCountDescriptor_EPA_DSS<- cbind(SMILES,RotatableBondsCountDescriptor_EPA_DSS)



MD <- FractionalCSP3Descriptor |>
   inner_join(SmallRingDescriptor, by = c(SMILES = "SMILES"))

 MD <- MD |>
   inner_join(FractionalPSADescriptor, by = c(SMILES = "SMILES"))

 MD <- MD |>
   inner_join(ZagrebIndexDescriptor, by = c(SMILES = "SMILES"))

 MD <- MD |>
   inner_join(XLogPDescriptor, by = c(SMILES = "SMILES"))

 MD <- MD |>
   inner_join(WienerNumbersDescriptor, by = c(SMILES = "SMILES"))

 MD <- MD |>
   inner_join(WeightDescriptor, by = c(SMILES = "SMILES"))

 MD <- MD |>
   inner_join(VAdjMaDescriptor, by = c(SMILES = "SMILES"))


 MD <- MD |>
   inner_join(HBondDonorCountDescriptor, by = c(SMILES = "SMILES"))

 MD <- MD |>
   inner_join(HBondAcceptorCountDescriptor, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(AminoAcidCountDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(AcidicGroupCountDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(AtomCountDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(AutocorrelationDescriptorMass_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(AutocorrelationDescriptorPolarizability_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(BasicGroupCountDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(BondCountDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(ChiChainDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(ChiClusterDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(ChiPathClusterDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(ChiPathDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join( EccentricConnectivityIndexDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(FMFDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join( FragmentComplexityDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(FragmentComplexityDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join( KappaShapeIndicesDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(  KierHallSmartsDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(LargestChainDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(  LargestPiSystemDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(MannholdLogPDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join( MDEDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join(   PetitjeanNumberDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))

MD <- MD |>
  inner_join( RotatableBondsCountDescriptor_EPA_DSS, by = c(SMILES = "SMILES"))


saveRDS(MD, file = "Data/molecularDescriptors_new.rds")

# join MACCS fingerprints 
MDs_with_CAS = chemicals |>
  full_join(MD, by = c(SMILES = "SMILES"))

saveRDS(MDs_with_CAS, file = "Data/molecularDescriptors_CAS.rds")
chemical_traits <- readRDS("Data/molecularDescriptors_CAS.rds")

### The following Molecular Descriptors did not work 

# 15. PetitjeanShapeIndexDescriptor
#PetitjeanShapeIndexDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor", descNames)])
## too many NAs


# 21. LengthOverBreadthDescriptor
#LengthOverBreadthDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor", descNames)])
## too many warnings

# 29. GravitationalIndexDescriptor
#GravitationalIndexDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor", descNames)])
# too many warnings 

# 20. LongestAliphaticChainDescriptor
# LongestAliphaticChainDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor", descNames)])
# SMILES <- rownames(LongestAliphaticChainDescriptor_EPA_DSS)
# rownames(LongestAliphaticChainDescriptor_EPA_DSS) <- NULL
# LongestAliphaticChainDescriptor<- cbind(SMILES,LongestAliphaticChainDescriptor_EPA_DSS) # Müll


# 33. CPSADescriptor
#CPSADescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor", descNames)])
## warnings 
# # 48. AromaticAtomsCountDescriptor
# AromaticAtomsCountDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor", descNames)])
### nur 0en 


# # 17. MomentOfInertiaDescriptor
# MomentOfInertiaDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor", descNames)])

# # 26. HybridizationRatioDescriptor
# HybridizationRatioDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor", descNames)])
# all_descriptors_EPA_DSS <- cbind(all_descriptors_EPA_DSS, HybridizationRatioDescriptor_EPA_DSS)


# # 39. BPolDescriptor
# BPolDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor", descNames)])
# zu viele warnings 

# 41. BCUTDescriptor
# BCUTDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor", descNames)])
# all_descriptors_EPA_DSS <- cbind(all_descriptors_EPA_DSS, BCUTDescriptor_EPA_DSS)
# ## zu viele wanrnings 

# 45. AutocorrelationDescriptorCharge
#AutocorrelationDescriptorCharge_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge", descNames)])
# zu viele warnings 

# 47. AromaticBondsCountDescriptor
#AromaticBondsCountDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor", descNames)])
# nur 0en 

# 49. APolDescriptor
# APolDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor", descNames)])
## zu viele warnings 
# # 50. ALOGPDescriptor
# ALOGPDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor", descNames)])
## zu viele warnings 

# 52. TaeAminoAcidDescriptor
#TaeAminoAcidDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.protein.TaeAminoAcidDescriptor", descNames)])
## zu viele warnings


# 11. VABCDescriptor
# VABCDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor", descNames)])
# SMILES <- rownames(VABCDescriptor_EPA_DSS)    
# rownames(VABCDescriptor_EPA_DSS) <- NULL
# VABCDescriptor <- cbind(SMILES,VABCDescriptor_EPA_DSS) ### not working well 


# 12. TPSADescriptor
# TPSADescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor", descNames)])
# SMILES <- rownames(TPSADescriptor_EPA_DSS)
# rownames(TPSADescriptor_EPA_DSS) <- NULL
# TPSADescriptor_EPA_DSS_SMILES <- cbind(SMILES,TPSADescriptor_EPA_DSS) # not working well 


# 13. RuleOfFiveDescriptor
# RuleOfFiveDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor", descNames)])
# SMILES <- rownames(RuleOfFiveDescriptor_EPA_DSS)
# rownames(RuleOfFiveDescriptor_EPA_DSS) <- NULL
# RuleOfFiveDescriptor <- cbind(SMILES,RuleOfFiveDescriptor_EPA_DSS) # binär, wollen wir nicht unbedingt 

# # 7. WHIMDescriptor
# WHIMDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor", descNames)])
# SMILES <- rownames(WHIMDescriptor_EPA_DSS)
# rownames(WHIMDescriptor_EPA_DSS) <- NULL
# WHIMDescriptor_EPA_DSS_SMILES <- cbind(SMILES,WHIMDescriptor_EPA_DSS) ## didn't work well
# 
# # 8. WeightedPathDescriptor
# WeightedPathDescriptor_EPA_DSS <- eval.desc(mols_EPA_DSS, descNames[grep("org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor", descNames)])
# SMILES <- rownames(WeightedPathDescriptor_EPA_DSS)
# rownames(WeightedPathDescriptor_EPA_DSS) <- NULL
# WeightedPathDescriptor_EPA_DSS_SMILES <- cbind(SMILES,WeightedPathDescriptor_EPA_DSS) ### didn't work as well 






