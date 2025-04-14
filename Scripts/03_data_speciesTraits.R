#########################################
########## phylogenetic EV ##############
#########################################

LC50 <- readRDS("Data/LC50_all.rds")
colnames(LC50)[1] <- "CAS"
colnames(LC50)[3] <- "scientificNameStd"

install.packages("taxize")
library(taxize)
species <- unique(LC50$scientificNameStd)
### for taxsize v0.9.103
resolved_names1 <- gnr_resolve(sci = species, data_source_ids = 179) # 179 is id for Open Tree of Life Reference Taxonomy
gbfid_id <- get_gbifid(sci = resolved_names1$user_supplied_name)
classification_df <- classification(sci_id = resolved_names1$user_supplied_name ,id = gbfid_id, db = 'gbif', return_id = TRUE) # GBIF phylogeny
### click enter in Console, then NA will be filled in 
#### new (31.03.25)
resolved_names1 <- gna_verifier(names = species, data_source_ids = 179) # 179 is id for Open Tree of Life Reference Taxonomy
gbfid_id <- get_gbifid(sci = resolved_names1$submittedName)
classification_df <- classification(sci_id = resolved_names1$submittedName ,id = gbfid_id, db = 'gbif', return_id = TRUE) # GBIF phylogeny


#saveRDS(classification_df, file = "Data/classification_df.rds")
classification_df <- readRDS(file= "Data/classification_df.rds")

Phylo_tree <- class2tree(unique(classification_df), varstep = TRUE, check = TRUE)
str(Phylo_tree)
#saveRDS(Phylo_tree, file = "Data/Phylo_tree.rds")

plot(Phylo_tree)

## 
# Calculate Eigenvectors from distance matrix 
## load calculated PhyloTree 

Phylo_tree <- readRDS("Data/Phylo_tree.rds")
Phylo_matrix <- as.matrix(Phylo_tree$distmat)

# Calculation of Ähnlichkeitsmatrix 
#n <- nrow(Phylo_matrix)
#H <- diag(n) - (1/n) * matrix(1, n, n)
#S <- -0.5 * H %*% (Phylo_matrix^2) %*% H

# Calculation of Eigenvectors and Eigenvalues 
#Eigen<- eigen(S)
#eigen_vectors <- Eigen$vectors
#eigen_values <- Eigen$values

#View(eigen_vectors)

##
# Use a package to calculate Eigenvectors and Eigenvalues 

#install.packages("spatialRF")
library(spatialRF)

weights_distmatrix <- weights_from_distance_matrix(distance.matrix = Phylo_matrix)
results = mem(weights_distmatrix)
# 1# way: calculate weighted distance matrix, center weights, use eigen()
#double_centered_distmatrix <- double_center_distance_matrix (distance.matrix = weights_distmatrix)
#Eigen_matrix <- eigen(double_centered_distmatrix, symmetric = TRUE)
#eigenvectors <- (Eigen_matrix$vectors)[,1827:1877]
#eigenvalues <- (Eigen_matrix$values)
#View(eigenvectors)
#str(eigen_matrix)
# 2# way: use mem() --> da kommen unterschiedeliche Ergebnisse raus... 
#Eigenvector_map <- mem(distance.matrix = 1/Phylo_matrix)
#View(Eigenvector_map)
#str(Eigenvector_map)


######## map the two dataframes ###### 
# Sets dataframe with species names as rows 
eigenvectors_with_species<- data.frame(
  Taxon = rownames(weights_distmatrix),
  results
)

#########################################
############ other TRAITS ###############
#########################################

# Install devtools if not available
# if(!"remotes" %in% installed.packages()[,"Package"]) install.packages("remotes")

# Install traitdata package from Github
#remotes::install_github("RS-eco/traitdata", build_vignettes = T, force=T)
library(traitdata)
library(dplyr, quiet=TRUE)
library(ggplot2)
data(trait_glossary)
#vignette("data_info")
#vignette("trait_glossary")
#vignette("access-data")
require(utils)
try(data(package = "traitdata"), silent = TRUE)
df_amphibio <- as.data.frame(amphibio)
colnames(df_thermal)
df_arthropods <- as.data.frame(arthropods)
df_fish <- as.data.frame(fishmorph)
df_thermal <- as.data.frame(globTherm)

df_all_traits <- c(df_amphibio, df_arthropods, df_fish, df_thermal)
#summary(df_all_traits)

####### check dataframes 
# 1. was bedeuten die jeweiligen Variablen und was bedeuten die NAs? 
# 1.1 AmphiBIO Datset 
# https://www.nature.com/articles/sdata2017123/tables/2
# NAs sind vorsichtig zu interpretieren, weil  NA nicht bedeutet, dass bestimmtes Merkmal nicht unbedingt fehlen muss,
# aber zumindest (nach Wissen der Authoren) in der Literatur noch nie als solches beschrieben wurde. 
# Die interessantesten Variablen: Habitat eher nicht (eh nur aquatisch), Diet, Diel, Seasonality, Body_mass_g, Body_size_mm, Longevity_max_y (wie alt Organismen werden), Breeding strategy, Offspring_size_max_mm
# 1.2 Arthropods dataset 
# https://www.nature.com/articles/sdata201513/tables/3
# keine Angaben zu NAs 
# Die interessantesten Variablen: Body_size (mm), Feeding_guild, Feeding_mode, Feeding_plant_part
# 1.3 Fish dataset 
# https://onlinelibrary.wiley.com/doi/10.1111/geb.13395
# Nicht gemessene Traits werden als “NA” gecoded 
# Die interessantesten Variablen: MBI: Maximum Body length (in cm) sonst nur Fisch-spezifische, morphologische Traits
# 1.4 Thermal dataset 
# Interesting Variable: critical thermal maximum (Tmax)

####### rename columns of dataframes 
colnames(df_amphibio)[26] <- "body_size"
# impute NAs using missRanger
#install.packages("missRanger")
#library(missRanger)
#str(df_amphibio)
#df_amphibio_imputed <- missRanger(
#  data = df_amphibio, 
#  formula = .~.
#)
#df_amphibio_imputed$body_size = df_amphibio$body_size
#sum(is.na(df_amphibio$body_size))


colnames(df_arthropods)[7] <- "body_size"
colnames(df_fish)[7] <- "body_size"
df_fish$body_size <- df_fish$body_size * 10 # transform from cm to mm 


#amphibio_columns <- df_amphibio[, c(26, 39)]
#arthropods_columns <- df_arthropods[, c(7, 19)]
#fish_columns <- df_fish[, c(7, 19)]
#df_body_size <- rbind(amphibio_columns, arthropods_columns, fish_columns)



df_all_traits_1 = df_amphibio |>
  full_join(df_arthropods, by = c(scientificNameStd = "scientificNameStd", Species = "Species", Order = "Order", Family = "Family", Genus = "Genus",  body_size = "body_size"))
str(df_all_traits_1)
df_all_traits_2 = df_all_traits_1 |>
  full_join(df_fish, by = c(scientificNameStd = "scientificNameStd", Species = "Species", Order = "Order", Family = "Family", Genus = "Genus", body_size = "body_size"))
df_all_traits_3 = df_all_traits_2 |>
  full_join(df_thermal, by = c(scientificNameStd = "scientificNameStd", Species = "Species", Order = "Order", Family = "Family", Genus = "Genus"))
summary(df_all_traits_3)

#
## add phylogenetic Eigenvectors 
colnames(eigenvectors_with_species)[1] <- "scientificNameStd"

## only fish traits and Eigenvectors 
df_all_traits_4 = df_all_traits_3 |>
  full_join(eigenvectors_with_species, by = c(scientificNameStd = "scientificNameStd"))

saveRDS(df_all_traits_4, file = "Data/species_traits_dataframe.rds")


## only df_fish and eigenvectors 
df_fish_traits = df_fish |>
  full_join(eigenvectors_with_species, by = c(scientificNameStd = "scientificNameStd"))

saveRDS(df_fish_traits, file = "Data/fish_traits_dataframe.rds")
saveRDS(df_fish, file = "Data/only_fishTraits.rds")
#fish_traits <-readRDS("Data/fish_traits_dataframe.rds")


