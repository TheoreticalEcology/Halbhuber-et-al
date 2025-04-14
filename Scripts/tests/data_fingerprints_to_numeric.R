# load data
library(readr)
df_fish <- readRDS("/home/isabellehalbhuber/Toxicology/Data/df_fish.rds")
## convert columns to numeric 
convert_to_num <- function(df, columns) {
  for (col in columns) {
    df[[col]] <- as.numeric(df[[col]])
  }
  return(df)
}

columns <- c("BEl", "VEp", "REs", "OGp", "BLs", "PFs", "PFv", "CPt", "RMl") # only continous variables
df_fish <- convert_to_num(df_fish, columns)
df_fish[19] <- NULL # exclude meanConc

## convert fingerprints as character strings to numeric representation as a matrix 

convert_to_vector <- function(fingerprint) {as.integer(unlist(strsplit(fingerprint, split = "")))}

fingerprints_list <- lapply(df_fish$Fingerprint, convert_to_vector)

fingerprints_matrix <- do.call(rbind, fingerprints_list)

df_fish <- cbind(df_fish, fingerprints_matrix)


# give each fingerprint bit a categorical name (numeric not proccessable in ranger) 
numbers <- as.character(1:881)
new_colnames <- paste0("bit", numbers)
colnames(df_fish)[27:907] <- new_colnames

df_fish <- as.data.frame(df_fish)

saveRDS(df_fish, file = "/home/isabellehalbhuber/Toxicology/Data/df_fish_fingerprints.rds")