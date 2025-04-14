# source model
library(readr)
library(torch)
library(dplyr)
## final RNN + Self-Attention model
source("/home/isabellehalbhuber/Toxicology/Scripts/BiRNNattn_model.R")

#################################################################################
###############################   Input Data   ##################################
#################################################################################

unique_df <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embedding_ready_df_withSpeciesCV.rds")

filtererd <- unique_df %>% select(SMILES, logmeanConc)
View(filtererd)
lowesttoxicity <- unique_df[unique_df$SMILES == "OCCOCCOCCO",]
highesttoxicity <- unique_df[unique_df$SMILES == "ClC1=C(Cl)C2(Cl)C3COS(=O)OCC3C1(Cl)C2(Cl)Cl",] ## 1988

filtererd_2 <- highesttoxicity %>% select(SMILES, logmeanConc)

unique_df <- na.omit(unique_df)
#unique_df[, 179:195] <- lapply(unique_df[, 179:195], sample)
#embeddedTokens <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embeddedTokens.rds")

# Chemical and Species Traits
data_matrix <- apply(as.matrix(unique_df[,2:195]), 2, as.numeric) 
data_matrix_tensor <- torch_tensor(data_matrix, dtype = torch_long()) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())
traits_tensor <- torch_tensor(scale(data_matrix[,178:194]), dtype = torch_float32())


#################################################################################
###########################   Re-Initialize Model   #############################
#################################################################################

## die gespeicherten modell parameter: gewichte und biases 
loaded_weights_1 <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/Final_TRAIT_model_weights_test_1.RDS")
loaded_weights_2 <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/Final_TRAIT_model_weights_test_2.RDS")
loaded_weights_3 <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/Final_TRAIT_model_weights_test_3.RDS")
loaded_weights_4 <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/Final_TRAIT_model_weights_test_4.RDS")
loaded_weights_5 <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/Final_TRAIT_model_weights_test_5.RDS")
loaded_weights_pre <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/Final_TRAIT_model_weights_prelayernormalization.RDS")

## modell neu initialisieren mit den gespeicherten Parametern (dann muss es nicht noch mal trainiert werden)
model <- BiRNNattn(embedding = 10L)
#model <- BiRNN$new()
#model$to(device = "cuda:3")
for (i in 1:length(loaded_weights_3)) {
  if(!stringr::str_detect(names(model$parameters[i]), "bias")) model$parameters[[i]]$set_data(new_data = loaded_weights_3[[i]])
  else model$parameters[[i]]$set_data(new_data = loaded_weights_3[[i]] %>% as.vector())
  #model$parameters[[i]]$data <- torch_tensor(loaded_weights[[i]], device = model$parameters[[i]]$device)
}

## train model, to get attention weigths matrix 
test_output_1 <- model(tokens_tensor, traits_tensor)
test_output_2 <- model(tokens_tensor, traits_tensor)
test_output_3 <- model(tokens_tensor, traits_tensor)
test_output_4 <- model(tokens_tensor, traits_tensor)
test_output_pre <- model(tokens_tensor, traits_tensor)


#as.matrix(model$attn_weights[1,,])[,-c(9:176)] %>% t() %>% fields::image.plot()


#################################################################################
############################   Attention Matrix  ################################
#################################################################################

#Attention_weights <- as_array(model$attn_weights)
Attention_weights_1 <- as_array(model$att_weights)
Attention_weights_2 <- as_array(model$att_weights)
Attention_weights_3 <- as_array(model$att_weights)
Attention_weights_4 <- as_array(model$att_weights)
Attention_weights_5 <- as_array(model$att_weights)


#test_attn <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/test_attn.RDS")

###### filter for highest and lowest prediction 
dim(Attention_weights)
lowestToxicity_attn <- Attention_weights_4[1988,,]
as.matrix(lowestToxicity_attn ) %>% t() %>% fields::image.plot()
heatmap(as.matrix(lowestToxicity_attn ))

highestToxicity_attn <- Attention_weights_4[879,,]
as.matrix(highestToxicity_attn ) %>% t() %>% fields::image.plot()
heatmap(as.matrix(highestToxicity_attn ))



## plot single chem x species interaction observations
# as.matrix(Attention_weights_5[42,,])[-c(178:194),-c(1:177)] %>% t() %>% fields::image.plot()

#calculate sum of attention weights for one observation
sum(as.vector((Attention_weights_5[42,,])[-c(178:194),-c(1:177)]))

# calculate sum of all weights per observation 
sum(as.vector((Attention_weights_1[42,,])))

# calculate sum of chem x chem weights 
sum(as.vector((Attention_weights_5[999,,])[-c(178:194),-c(178:194)]))

as.matrix(Attention_weights_5[642,,])[-c(178:194),-c(178:194)] %>% t() %>% fields::image.plot()



#### jetzt summe aller chem x chem attention weights 
# Annahme: Attention_weights_1 hat Dimensionen [N, M, M]
# Erstellt einen Vektor, der die Summen speichert
sums <- numeric(dim(Attention_weights_5)[1]) 
for (i in 1:dim(Attention_weights_5)[1]) {
  matrix_i <- Attention_weights_5[i,,]
  filtered_matrix <- matrix_i[-c(178:194),-c(1:177)]
  sums[i] <- sum(as.vector(filtered_matrix))
}
sum(sums)

tmp = Attention_weights_3
dim(tmp)
## name columns and rows of the matrix  
matrices_list <- list()
for (i in 1:dim(tmp)[1]) {
  
  matrix_i <- tmp[i, , ]
  
  tokens <- as.numeric(tokens_tensor[i, ])
  row_col_names <- c(as.character(tokens), paste0("trait_", 1:17)) 
  
  
  colnames(matrix_i) <- row_col_names
  rownames(matrix_i) <- row_col_names
  
  matrices_list[[i]] <- matrix_i
}

#all_row_names <- unique(unlist(lapply(tmp, rownames)))


## delete rownames and colnames == 1 

for (i in 1:length(matrices_list)) {
  current_matrix <- matrices_list[[i]]
  
  col_names <- colnames(current_matrix)
  row_names <- rownames(current_matrix)
  
  cols_to_delete <- which(col_names == "1")
  rows_to_delete <- which(row_names == "1")
  
  if (length(cols_to_delete) > 0) {
    current_matrix <- current_matrix[, -cols_to_delete, drop = FALSE]  
  }
  if (length(rows_to_delete) > 0) {
    current_matrix <- current_matrix[-rows_to_delete, , drop = FALSE]  
  }
  
  matrices_list[[i]] <- current_matrix
}

attention_interactions <- matrices_list

#saveRDS(attention_interactions, file = "/home/isabellehalbhuber/Toxicology/Results/attention_interactions.rds")

#Attention_scores <- readRDS( file = "/home/isabellehalbhuber/Toxicology/Results/attention_interactions.rds")

## sum all matrices up 
all_row_names <- unique(unlist(lapply(attention_interactions, rownames)))
all_col_names <- unique(unlist(lapply(attention_interactions, colnames)))


# now for all entries
# generate empty results matrix, ready to be filled up with sum of each entry
result_matrix <- matrix(0, nrow = length(all_row_names), ncol = length(all_col_names),
                        dimnames = list(all_col_names, all_row_names))

# for all matrices 
for (mat in attention_interactions) {
  # rownames() are rownames of all matrices
  for (row_name in rownames(mat)) {
    # colnames() are colnames of all matrices 
    for (col_name in colnames(mat)) {
      # search for entries with matchin row-col names combinations 
      if (row_name %in% rownames(result_matrix) && col_name %in% colnames(result_matrix)) {
        # sum up the entries
        result_matrix[row_name, col_name] <- result_matrix[row_name, col_name] + mat[row_name, col_name]
      } else {
        # if the combination does not exist, then its going to be added 
        if (!(row_name %in% rownames(result_matrix))) {
          result_matrix <- rbind(result_matrix, matrix(0, nrow = 1, ncol = ncol(result_matrix),
                                                       dimnames = list(row_name, colnames(result_matrix))))
        }
        if (!(col_name %in% colnames(result_matrix))) {
          result_matrix <- cbind(result_matrix, matrix(0, nrow = nrow(result_matrix), ncol = 1,
                                                       dimnames = list(rownames(result_matrix), col_name)))
        }
        result_matrix[row_name, col_name] <- mat[row_name, col_name]
        
      }
    }
  }
}


#saveRDS(result_matrix, file = "/home/isabellehalbhuber/Toxicology/Results/trait_interaction_matrix.rds")
#result_matrix <- readRDS(file = "/home/isabellehalbhuber/Toxicology/Results/trait_interaction_matrix.rds")

#################################################################################
########################   Normalize the sumed up ones  #########################
#################################################################################

##### hier muss ich jetzt alle colums Werte durch die Anzahl an Tokens teilen
# 1. Welche Tokens kommen wie oft in meiner Liste von Matrizen vor? 
#all_result_matrices <- readRDS(file = "/home/isabellehalbhuber/Toxicology/Results/attention_interactions.rds")

all_result_matrices <- attention_interactions


get_unique_colnames_and_frequencies <- function(all_result_matrices) {
  all_colnames <- c()
  for (matrix in all_result_matrices) {
    colnames_matrix <- attr(matrix, "dimnames")[[2]]
    all_colnames <- c(all_colnames, colnames_matrix)
  }
  
  colname_frequencies <- table(all_colnames)
  return(as.data.frame(colname_frequencies))
}
token_frequencies <- get_unique_colnames_and_frequencies(attention_interactions)
#saveRDS(token_frequencies, file = "/home/isabellehalbhuber/Toxicology/Results/token_frequencies.rds")


#print(unique_colnames_frequencies)

# 2. Teile die jeweiligen Spalten durch die Anzahl des jeweiligen Tokens 
library(dplyr)
InteractionMatrix <- result_matrix
#InteractionMatrix <- readRDS(file = "/home/isabellehalbhuber/Toxicology/Results/trait_interaction_matrix.rds")

#InteractionMatrix <- result_matrix
colnames(InteractionMatrix)[-(3:19)] <- as.numeric(colnames(InteractionMatrix)[-(3:19)])
rownames(InteractionMatrix)[-(3:19)] <- as.numeric(rownames(InteractionMatrix)[-(3:19)])
#token_frequencies <- readRDS(file = "/home/isabellehalbhuber/Toxicology/Results/token_frequencies.rds")
# filter out the traits 
#filtered_data <- token_frequencies[!grepl("^trait_([1-9]|1[0-7])$", token_frequencies$all_colnames), ]
##### divide interactionMatrix through number of each token 
freq_vector <- setNames(token_frequencies$Freq, token_frequencies$all_colnames)
matching_colnames <- intersect(colnames(InteractionMatrix), names(freq_vector))
###
InteractionMatrix[, matching_colnames] <- InteractionMatrix[, matching_colnames] / freq_vector[matching_colnames]





### Rename the tokens to SMILES codes 
TokensLookUpTable <- read.csv(file = "/home/isabellehalbhuber/Toxicology/Results/TokensLookUpTable.csv")
TokensLookUpTable$int <- as.integer(TokensLookUpTable$int +1)

## replace numeric representation of tokens with real chemical structures  for (i in 1:nrow(position_df)) {

replace_dimnames_with_lookup <- function(mat, lookup) {
  dimnames_list <- dimnames(mat)
  
  lookup_replace <- function(x) {
    replaced_values <- character(length(x))
    
    for (i in seq_along(x)) {
      
      if (grepl("^\\d+$", x[i])) {
        int_val <- as.numeric(x[i])  
        token <- lookup$token[lookup$int == int_val]
        if (length(token) > 0) {
          replaced_values[i] <- token
        } else {
          replaced_values[i] <- x[i]  
        }
      } else {
        replaced_values[i] <- x[i]  
      }
    }
    return(replaced_values)
  }
  
  new_dimnames <- lapply(dimnames_list, lookup_replace)
  
  dimnames(mat) <- new_dimnames
  
  return(mat)
}

matrix_data_new <- replace_dimnames_with_lookup(InteractionMatrix, TokensLookUpTable)

InteractionMatrix_SMILES <- matrix_data_new


saveRDS(InteractionMatrix_SMILES, file = "/home/isabellehalbhuber/Toxicology/Results/InteractionMatrix_SMILES_test_2.rds")
InteractionMatrix_SMILES <- readRDS(file = "/home/isabellehalbhuber/Toxicology/Results/InteractionMatrix_SMILES_test_3.rds")

max(InteractionMatrix_SMILES)
# # order matrix 
# 
desired_order <- paste0("trait_", 1:17) # 1:17
existing_columns <- colnames(InteractionMatrix_SMILES)
traits_to_keep <- desired_order[desired_order %in% existing_columns]
other_columns <- existing_columns[!(existing_columns %in% traits_to_keep)]
new_column_order <- c(traits_to_keep, other_columns)
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[, new_column_order, drop = FALSE]

desired_order_rows <- paste0("trait_", 1:17) #1:17
existing_rows <- rownames(InteractionMatrix_SMILES)
rows_to_keep <- desired_order_rows[desired_order_rows %in% existing_rows]
other_rows <- existing_rows[!(existing_rows %in% rows_to_keep)]
new_row_order <- c(rows_to_keep, other_rows)
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[new_row_order, , drop = FALSE]

new_names <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", "Pectoral.Fin.Size",
               "Caudal.Peduncle.Throttling", "body_size", paste0("pEV_", 1:7))
colnames(InteractionMatrix_SMILES)[1:17] <- new_names #17
rownames(InteractionMatrix_SMILES)[1:17] <- new_names #17

### wie sind die verhältnisse?
part_1 <- InteractionMatrix_SMILES[18:441,1:17]
sum(as.vector(part_2))
part_2 <- InteractionMatrix_SMILES[18:441, 18:441]



colnames(InteractionMatrix_SMILES)
# ##### chem x chem
## nur ein bestimmter Bereich für den plot

#InteractionMatrix_SMILES <- InteractionMatrix_SMILES[(73:100), (73:100)]
undesired_cols <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position",
                  "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position",
                  "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
#desired_rows <- c("OC(C)C)", "CCCN)", "C2CC2)", "CCN(C(=O)", "(=O)=O", "C(NC(=O)")
desired_rows <- rownames(InteractionMatrix_SMILES)[!(rownames(InteractionMatrix_SMILES) %in% undesired_cols)]
#desired_rows <- rownames(InteractionMatrix_SMILES)[(rownames(InteractionMatrix_SMILES) %in% desired_rows)]

desired_cols <- desired_rows

InteractionMatrix_SMILES <- InteractionMatrix_SMILES[rownames(InteractionMatrix_SMILES) %in% desired_rows, ]
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[,colnames(InteractionMatrix_SMILES) %in% desired_cols]





####### chem x spec 
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[(73:100), (1:17)]
desired_cols <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position",
                  "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position",
                  "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
desired_rows <- rownames(InteractionMatrix_SMILES)[!(rownames(InteractionMatrix_SMILES) %in% desired_cols)]
# #desired_cols <- c("Cl","C1=","CC","=C(","C=","C1)","C(","Cl)","C(Cl)(Cl)","C1=C(","=C1","O","C=C(","C=C1","OCC","CC1=","CC2","CO","P(=O)(","OC)","P(","=S)","(","CC(C)","OCC)","C(=O)O","C1","C(Cl)","N","CC=","=","CC(","S","=C", "2)" ,"C(C)","=N","[N+]([O-])=O)","C1=N","C(N")
#
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[rownames(InteractionMatrix_SMILES) %in% desired_rows, ]
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[,colnames(InteractionMatrix_SMILES) %in% desired_cols]

# ######## spec x chem
 desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position",
                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position",
                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
# 
# desired_cols <- colnames(InteractionMatrix_SMILES)[!(colnames(InteractionMatrix_SMILES) %in% desired_rows)]
# # #desired_cols <- c("Cl","C1=","CC","=C(","C=","C1)","C(","Cl)","C(Cl)(Cl)","C1=C(","=C1","O","C=C(","C=C1","OCC","CC1=","CC2","CO","P(=O)(","OC)","P(","=S)","(","CC(C)","OCC)","C(=O)O","C1","C(Cl)","N","CC=","=","CC(","S","=C", "2)" ,"C(C)","=N","[N+]([O-])=O)","C1=N","C(N")
# #
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[rownames(InteractionMatrix_SMILES) %in% desired_rows, ]
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[,colnames(InteractionMatrix_SMILES) %in% desired_rows]
# # 
# 
# #filter and save column names
filtered_rows <- list()
filtered_row_colnames <- list()

# Berechne das Threshold für die Top 10% Werte
threshold <- quantile(unlist(InteractionMatrix_SMILES), 0.999)

# Filter-Schleife über jede Zeile der Matrix
for (i in 1:nrow(InteractionMatrix_SMILES)) {
  # Filtere Werte, die über dem Threshold liegen
  filtered_row <- InteractionMatrix_SMILES[i, InteractionMatrix_SMILES[i, ] > threshold]
  
  if (length(filtered_row) > 0) {
    # Speichere gefilterte Werte und deren Spaltennamen
    filtered_rows[[rownames(InteractionMatrix_SMILES)[i]]] <- filtered_row
    filtered_row_colnames[[rownames(InteractionMatrix_SMILES)[i]]] <- colnames(InteractionMatrix_SMILES)[InteractionMatrix_SMILES[i, ] > threshold]
  }
}

# Erstelle eine gefilterte Matrix mit der gleichen Spaltenanzahl
max_columns <- ncol(InteractionMatrix_SMILES)

filtered_matrix <- do.call(rbind, lapply(names(filtered_rows), function(row_name) {
  row <- filtered_rows[[row_name]]
  col_names <- filtered_row_colnames[[row_name]]
  #col_names <- colnames(InteractionMatrix_SMILES)
  
  # Erstelle eine Zeile mit NA-Werten
  full_row <- rep(NA, max_columns)
  
  # Füge gefilterte Werte an die korrekten Positionen ein
  full_row[match(col_names, colnames(InteractionMatrix_SMILES))] <- row
  return(full_row)
}))

# Benenne Zeilen und Spalten
rownames(filtered_matrix) <- names(filtered_rows)
colnames(filtered_matrix) <- colnames(InteractionMatrix_SMILES)
### save rownmaes with highest interactions, to use them then as 'desired_rows'
rows_99 <- rownames(filtered_matrix)

InteractionMatrix_SMILES <- InteractionMatrix_SMILES[rownames(InteractionMatrix_SMILES) %in% rows_99, colnames(InteractionMatrix_SMILES) ]
#InteractionMatrix_SMILES <- InteractionMatrix_SMILES[1:17,1:17]

# color_palette <- colorRampPalette(my_colors)
# chemXspec_5_gr0 <- InteractionMatrix_SMILES
# hist(as.vector(chemXspec_5_gr0))
# chemXspec_4_gr0 <- InteractionMatrix_SMILES
# hist(as.vector(chemXspec_4))
# chemXspec_3_gr0 <- InteractionMatrix_SMILES
# hist(as.vector(chemXspec_3_gr0))
# chemXspec_2_gr0 <- InteractionMatrix_SMILES
# hist(as.vector(chemXspec_2_gr0))
# chemXspec_1_gr0 <- InteractionMatrix_SMILES
# hist(as.vector(chemXspec_1_gr0))


my_colors <- c("white", "#edf8b1","#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58", "black")
my_colors <- c("white", "#edf8b1")
# wie werden die Farben in der Matrix verteilt? 
value_range <- range(InteractionMatrix_SMILES, na.rm = TRUE) 
min_value <- value_range[1]
max_value <- value_range[2]

n_colors <- length(my_colors)
value_intervals <- seq(min_value, max_value, length.out = n_colors + 1)

# color_mapping <- data.frame(
#   Color = my_colors,
#   Min_Value = value_intervals[-(n_colors + 1)],
#   Max_Value = value_intervals[-1]
# )


color_mapping <- data.frame(
  Min_Value = value_intervals[-(n_colors + 2)],
  Color = c(my_colors, NA) 
)




layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1)) # Plot links, Skala rechts
value_centers <- (value_intervals[-1] + value_intervals[-length(value_intervals)]) / 2



pdf("chemxchem_highest10.pdf", width = 6, height = 5) 
heatmap_chemxspec <- heatmap(
  InteractionMatrix_SMILES,
  main = "Trait Interactions",
  Colv = NA,
  Rowv = NA,
  las = 1,
  col = my_colors,
  cexRow = 1,
  cexCol = 1,
  labRow = rownames(InteractionMatrix_SMILES), 
  labCol = colnames(InteractionMatrix_SMILES) 
)


dev.off()

pdf("Farbskala_test5_speciesxchem.pdf", width = 6, height = 2 ) 
image(
  x = seq(min(color_mapping$Min_Value), max(color_mapping$Min_Value), length.out = nrow(color_mapping)),
  y = 1,
  z = t(matrix(seq_along(my_colors), nrow = 1)),
  col = my_colors,
  yaxt = "n", 
  xaxt = "n", 
  xlab = "Interaction intensity",
  ylab = ""
)

axis(
  side = 1,                                
  at = color_mapping$Min_Value,           
  labels = round(color_mapping$Min_Value, 2),
  las = 1                                  
)

dev.off()
