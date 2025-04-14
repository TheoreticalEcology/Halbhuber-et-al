#################################################################################
############################   Tokenization   ###################################
#################################################################################

### from smiles.csv to smiles-smi and tokenization done in python 

library(readr)
tokenizedSMILES <- read.table("/home/isabellehalbhuber/Toxicology/Data/smiles_tokenized_dataframe.csv",
                              header = FALSE,
                              sep = ",",
                              stringsAsFactors =TRUE,
                              fill = TRUE,
                              check.names = TRUE
                              )
View(tokenizedSMILES)
colnames(tokenizedSMILES)[1] <- "SMILES"
colnames(tokenizedSMILES)[2] <- "Tokens"
tokenizedSMILES <- tokenizedSMILES[-1,]

tokenized_smiles <- as.character(tokenizedSMILES$Tokens)


# function to split through ", ,"
split_tokens <- function(smiles) {
  if (!is.character(smiles)) {
    stop("no character.")
  }
  tokens <- strsplit(smiles, ", ,")[[1]]
  # remove empty spaces at the beginning and end of tokens 
  tokens <- trimws(tokens)
  return(tokens)
}

# split tokenized SMILES
split_smiles <- lapply(tokenized_smiles, function(smiles) {
  tryCatch({
    split_tokens(smiles)
  }, error = function(e) {
    cat("Error in:", smiles, "\n")
    cat("Error", e$message, "\n")
    return(character())  
  })
})

# calculate max tokens
max_tokens <- max(sapply(split_smiles, function(tokens) {
  if (length(tokens) == 0) {
    return(0)
  } else {
    return(length(tokens))
  }
}), na.rm = TRUE)

# generate df with SMILES 
original_smiles <- as.character(tokenizedSMILES$SMILES)
df <- data.frame("SMILES" = original_smiles)

# add tokens 
for (i in 1:max_tokens) {
  df[[paste0("Token", i)]] <- sapply(split_smiles, function(tokens) {
    if (i <= length(tokens)) {
      return(tokens[i])
    } else {
      return(NA)
    }
  })
}

# remove additional kommas and empty spaces 
df <- data.frame(lapply(df, function(col) {
  if (is.character(col)) {
    col <- gsub(",", "", col)  
    col <- gsub(" ", "", col)  
  }
  return(col)
}), stringsAsFactors = FALSE)

View(df)
#### vectorization

smiles_vectors <- list()

for (i in 1:nrow(df)) {
  smiles_name <- df[i, 1]  
  smiles_vector <- as.vector(df[i, 2:ncol(df)])  
  names(smiles_vector) <- paste0("Token", 1:length(smiles_vector))  
  smiles_vectors[[smiles_name]] <- smiles_vector
}

saveRDS(df, file = "/home/isabellehalbhuber/Toxicology/Data/Tokens.rds")
write.csv(df, file = "/home/isabellehalbhuber/Toxicology/Data/Tokens.csv")
saveRDS(smiles_vectors, file = "/home/isabellehalbhuber/Toxicology/Data/vectorTokens.rds")


