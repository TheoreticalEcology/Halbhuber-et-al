#################################################################################
#######################      SMILES preparation       ###########################
#################################################################################
library(readr)
library(dplyr)
Tokens <- readRDS("Data/Tokens.rds")
df <- readRDS("Data/df_final.rds")
df <- df[, -c(4:210)] # exclude FP and MD 
df <- df %>% select(-meanConc)

# join tokens and df 
df = df |>
   inner_join(Tokens, by = c(SMILES = "SMILES"))

df1 <- df %>% select(-CAS, -scientificNameStd)

## only token
 df2 <- df1 %>% select(-SMILES, -logmeanConc, -meanConc, -randomFold, -speciesFold, -chemicalFold)

 
##########################################
#### how many unique tokens in dataset with toxicity values without NAs
all_tokens_noNA <- unlist(df1[, 23:199])
all_tokens_noNA <- all_tokens_noNA[!is.na(all_tokens_noNA)]
all_unique_tokens <- unique(all_tokens_noNA)

### jetzt die integer values für alle unique tokens erstellen
token_to_int <- as.data.frame(all_unique_tokens)
colnames(token_to_int)[1] <- "token"
token_to_int$int <- 1:614
colnames(token_to_int)[2] <- "int"
 
write.csv(token_to_int,"Results/TokensLookUpTable.csv", row.names = FALSE)

### give each position of the token an integer value 
numbers <- as.character(1:177)
new_colnames <- paste0("position_", numbers)
colnames(df1)[23:199] <- new_colnames

position_df <- df1
int_position_df <- position_df

 for (i in 1:nrow(position_df)) {
   for (j in 23:ncol(position_df)) {
     token <- position_df[i, j]
     if (!is.na(token)) {
       int_value <- token_to_int$int[token_to_int$token == token]
       int_position_df[i, j] <- ifelse(length(int_value) > 0, int_value, 0)
     } else {
       int_position_df[i, j] <- 0
     }
   }
 }

# df without doppelgänger
#unique_df_withSpecies <- int_position_df %>% # 921 unique Tokens und 0
#distinct()

str(int_position_df) ## Folds and tokens are factors and characters respecitvely -> need to be converted to numeric 

convert_to_num <- function(df, columns) {
  for (col in columns) {
    df[[col]] <- as.numeric(df[[col]])
  }
  return(df)
}
colnames(int_position_df[, 20: 199])
columns <- c(colnames(int_position_df[, 20: 199])) 
int_position_df <- convert_to_num(int_position_df, columns)
str(int_position_df)


saveRDS(int_position_df, file = "Data/df_final_tokenizedSMILES.rds")

