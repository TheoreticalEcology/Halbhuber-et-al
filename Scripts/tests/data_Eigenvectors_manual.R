#########################################
########## phylogenetic EV ##############
#########################################
LC50 <- readRDS("/home/isabellehalbhuber/Toxicology/Data/LC50_all.rds")
colnames(LC50)[1] <- "CAS"
colnames(LC50)[3] <- "scientificNameStd"
# Bibliotheken laden
library(taxize)

classification_df <- readRDS("/home/isabellehalbhuber/Toxicology/Data/classification_df.rds")
# Phylogenetischen Baum erstellen
install.packages("phytools")

library(phytools) 
Phylo_tree <- class2tree(unique(classification_df), varstep = TRUE, check = TRUE)
str(Phylo_tree)

# Phylogenetischen Baum speichern (optional)
saveRDS(Phylo_tree, file = "/home/isabellehalbhuber/Toxicology/Data/Phylo_tree_test.rds")

# Phylogenetischen Baum visualisieren
plot(Phylo_tree)

# Distanzmatrix aus dem phylogenetischen Baum extrahieren
Phylo_matrix <- as.matrix(Phylo_tree$distmat)

# Manuelle Berechnung der Ähnlichkeitsmatrix mit Doppeltzentrierung
n <- nrow(Phylo_matrix)
H <- diag(n) - (1/n) * matrix(1, n, n)  # Zentrierungsmatrix
S <- -0.5 * H %*% (Phylo_matrix^2) %*% H  # Ähnlichkeitsmatrix


Eigen <- eigen(S)
eigen_vectors <- Eigen$vectors  
eigen_values <- Eigen$values   


# Eigenvektoren mit Artnamen verknüpfen
eigenvectors_with_species <- data.frame(
  Taxon = rownames(Phylo_matrix),
  eigen_vectors
)

# Ausgabe der ersten paar Zeilen zur Überprüfung
head(eigenvectors_with_species)

# Die Ergebnisse als Dataframe anzeigen
str(eigenvectors_with_species)
View(eigenvectors_with_species)


saveRDS(eigenvectors_with_species, file = "/home/isabellehalbhuber/Toxicology/Data/eigenvectors_with_species_manual.rds")
