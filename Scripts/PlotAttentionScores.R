### initialize model with saved weights 
library(torch)
## Model 
BiRNN <- nn_module(
  initialize = function(embedding = 40L, rnn_hidden_size = 20L, rnn_layers = 3L, max_len = 177L) {
    self$embedding <- nn_embedding(
      num_embeddings = 425L, 
      embedding_dim = embedding, 
      sparse = FALSE, 
      scale_grad_by_freq = TRUE
    )
    
    # self$pos_encoding <- self$generate_positional_encoding(max_len, embedding)
    
    self$rnn <- nn_gru(
      input_size = embedding,
      hidden_size = rnn_hidden_size,
      num_layers = rnn_layers,
      bidirectional = TRUE,
      batch_first = TRUE,
      dropout = 0 # war davor 0
    )
    
    self$self_att = torch::nn_multihead_attention(embed_dim = 40L, num_heads = 1L, dropout = 0.1, batch_first = TRUE) # war davor 0.1
    
    # Self-Attention QKV Layer
    self$qkv_size <- rnn_hidden_size * 2  # eigentlich mal 2
    self$q_linear <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = self$qkv_size) # eigentlich mal 2
    self$k_linear <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = self$qkv_size) # eigentlich mal 2
    self$v_linear <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = self$qkv_size) # eigentlich mal 2
    
    self$token_full <- nn_linear(in_features = rnn_hidden_size * 2, out_features = 20L) # davor mal 2
    
    self$trait_full1 <- nn_linear(in_features = 2L, out_features = 100L)
    self$dropout1 = nn_dropout(0.3)
    self$trait_full2 <- nn_linear(in_features = 100L, out_features = 100L)
    self$dropout2 = nn_dropout(0.3)
    self$trait_full3 <- nn_linear(in_features = 100L, out_features = 40L) # 40 is embedding size 
    
    # Combined output layers
    self$output_layer1 <- nn_linear(in_features = 20L, out_features = 50L)
    self$output_layer2 <- nn_linear(in_features = 50L, out_features = 1L)
  },
  # 
  # generate_positional_encoding = function(max_len, embedding_dim) {
  #   pos_enc = torch_zeros(max_len, embedding_dim)
  #   positions = torch_arange(1, max_len)$unsqueeze(2)
  #   div_term = torch_exp(torch_arange(1, embedding_dim, 2) * (-torch_log(torch_tensor(10000)) / embedding_dim))
  #   
  #   pos_enc[, seq(1, embedding_dim, 2)] = torch_sin(positions * div_term)
  #   pos_enc[, seq(2, embedding_dim, 2)] = torch_cos(positions * div_term)
  #   
  #   return(pos_enc$unsqueeze(1))  # [max_len, 1, embedding_dim]
  # },
  
  forward = function(tokens, traits) {
    #### Token Embedding
    x <- self$embedding(tokens)
    
    # Add Positional Encoding to the Token Embeddings
    # seq_len <- tokens$size(2)
    # pos_enc = self$pos_encoding[1:seq_len, , ]$to(dtype = torch_float32(), device = x$device)
    # x <- x + pos_enc$squeeze(1)
    #print(dim(x))
    #### RNN Output
    rnn_out <- self$rnn(x)[[1]]  # [batch_size, seq_len, 2 * rnn_hidden_size]
    #print(dim(rnn_out))
    
    #### Trait Embedding
    traits = traits$unsqueeze(3L)
    traits_position = torch_tensor(matrix(scale(1:17), 17L, 1L),  # skalieren nicht vergessen
                                   device=traits$device,
                                   dtype=traits$dtype)$`repeat`(list(tokens$shape[1], 1L, 1L))
    traits = torch_cat(list(traits, traits_position), dim = 3L)
    
    traits_embedded = self$trait_full1(traits) %>% nnf_relu() %>% self$dropout1()
    traits_embedded = self$trait_full2(traits_embedded) %>% nnf_relu() %>% self$dropout2()
    traits_embedded = self$trait_full3(traits_embedded)
    #print(dim(traits_embedded))
    
    
    #### Concatenation of Token and Trait Embeddings
    #rnn_out_traits = torch_cat(list(rnn_out, traits_embedded), dim = 2L)
    rnn_out_traits = torch_cat(list(rnn_out, traits_embedded), dim = 2L)
    #print(dim(rnn_out_traits))
    
    
    # QKV for Self-Attention
    q <- self$q_linear(rnn_out_traits)  # [batch_size, seq_len, qkv_size]
    k <- self$k_linear(rnn_out_traits)  # [batch_size, seq_len, qkv_size]
    v <- self$v_linear(rnn_out_traits)  # [batch_size, seq_len, qkv_size]
    
    #### Masking for Padding Tokens
    padding_mask <- tokens$clone()
    #print(tokens$shape)
    
    mask <- padding_mask == 1L  # Assuming padding tokens are marked with 1 [batchsize, 177]
    mask = torch_cat(list(mask, torch_zeros(x$shape[1], traits$shape[2], dtype=torch_bool(), device=x$device)), dim = 2)
    
    
    #print(rnn_out$shape)
    #print(traits_embedded$shape)
    #print(rnn_out_traits$shape)
    #print(attn_out$shape)
    
    #### Multi-Head Self-Attention
    self$mask = mask
    attention_out <- self$self_att(q, k, v, key_padding_mask = mask )
    #print(dim(attention_out))
    self$att_weights = attention_out[[2]]
    # Take the mean along the sequence length dimension to summarize the attention output
    attn_out <- torch_mean(attention_out[[1]], dim = 2)  # [batch_size, embed_dim]
    #### Token Feature NN
    output_token <- self$token_full(attn_out) %>% nnf_relu()
    
    #### Combine outputs
    output <- self$output_layer1(output_token) %>% nnf_relu()
    output <- self$output_layer2(output)
    
    return(output)
  }
)



#################################################################################
###############################   Input Data   ##################################
#################################################################################


unique_df <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embedding_ready_df_withSpeciesCV.rds")
unique_df <- na.omit(unique_df)
#embeddedTokens <- readRDS("/home/isabellehalbhuber/Toxicology/Data/embeddedTokens.rds")

# Chemical and Species Traits
data_matrix <- apply(as.matrix(unique_df[,2:195]), 2, as.numeric) 
data_matrix_tensor <- torch_tensor(data_matrix, dtype = torch_long()) 
data_matrix_tokens <- data_matrix[,1:177] + 1 # +1, because Token indices shall not start at 0 
tokens_tensor <- torch_tensor(data_matrix_tokens[,1:177], dtype = torch_long())
traits_tensor <- torch_tensor(scale(data_matrix[,178:194]), dtype = torch_float32())

# LC50
data_matrix_LC50 <- as.numeric(unique_df[,196]) 
LC50_tensor <- torch_tensor(data_matrix_LC50, dtype = torch_float())

# Training and Test datasets
set.seed(123)
# Indices
train_indices <- sample(1:nrow(LC50_tensor), 0.8 * nrow(LC50_tensor))
test_indices <- setdiff(1:nrow(LC50_tensor), train_indices)

# Tokens
train_tokens_tensor <- tokens_tensor[train_indices, , drop = FALSE]
test_tokens_tensor <- tokens_tensor[test_indices, , drop = FALSE]
# Traits
train_traits_tensor <- traits_tensor[train_indices, , drop = FALSE]
test_traits_tensor <- traits_tensor[test_indices, , drop = FALSE]


# LC50
train_LC50_tensor <- LC50_tensor[train_indices]
test_LC50_tensor <- LC50_tensor[test_indices]


test_tokens_fold = train_tokens_tensor#$to(device = "cuda:3")
test_traits_fold = train_traits_tensor#$to(device = "cuda:3")
test_LC50_fold = train_LC50_tensor#$to(device = "cuda:3")

#################################################################################
###########################   Re-Initialize Model   #############################
#################################################################################

## die gespeicherten modell parameter: gewichte und biases 

loaded_weights <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/Final_TRAIT_model_weights_test_2.RDS")
## modell neu initialisieren mit den gespeicherten Parametern (dann muss es nicht noch mal trainiert werden)
model <- BiRNN(embedding = 10L)
#model <- BiRNN$new()
#model$to(device = "cuda:3")
for (i in 1:length(loaded_weights)) {
  if(!stringr::str_detect(names(model$parameters[i]), "bias")) model$parameters[[i]]$set_data(new_data = loaded_weights[[i]])
  else model$parameters[[i]]$set_data(new_data = loaded_weights[[i]] %>% as.vector())
  #model$parameters[[i]]$data <- torch_tensor(loaded_weights[[i]], device = model$parameters[[i]]$device)
}

## train model, to get attention weigths matrix 
test_output <- model(tokens_tensor, traits_tensor)
#as.matrix(model$attn_weights[1,,])[,-c(9:176)] %>% t() %>% fields::image.plot()


#################################################################################
############################   Attention Matrix  ################################
#################################################################################

#Attention_weights <- as_array(model$attn_weights)
Attention_weights <- as_array(model$att_weights)
tmp = Attention_weights
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



### trials
# res <- numeric(4409)  
# 
# for (i in 1:4409) {
#   res[i] <- attention_interactions[[i]]["trait_1", "trait_1"]
# }
# 
# print(res)
# 
# sum(res)


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

#saveRDS(InteractionMatrix_SMILES, file = "/home/isabellehalbhuber/Toxicology/Results/InteractionMatrix_SMILES.rds")
#InteractionMatrix_SMILES <- readRDS(file = "/home/isabellehalbhuber/Toxicology/Results/InteractionMatrix_SMILES.rds")

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
#


# #filter and save column names
filtered_rows <- list()
filtered_row_colnames <- list()
threshold <- quantile(unlist(InteractionMatrix_SMILES), 0.99)


for (i in 1:nrow(InteractionMatrix_SMILES)) {
  filtered_row <- InteractionMatrix_SMILES[i, InteractionMatrix_SMILES[i, ] > threshold]

  if (length(filtered_row) > 0) {

    filtered_rows[[rownames(InteractionMatrix_SMILES)[i]]] <- filtered_row
    filtered_row_colnames[[rownames(InteractionMatrix_SMILES)[i]]] <- colnames(InteractionMatrix_SMILES)[InteractionMatrix_SMILES[i, ] > 0.1]
  }
}

# save filtered matrix mit ncol = ncol of initial matrix
max_columns <- ncol(InteractionMatrix_SMILES)

filtered_matrix <- do.call(rbind, lapply(names(filtered_rows), function(row_name) {
  row <- filtered_rows[[row_name]]
  col_names <- filtered_row_colnames[[row_name]]

  full_row <- rep(NA, max_columns)
  full_row[match(col_names, colnames(InteractionMatrix_SMILES))] <- row

  return(full_row)
}))

rownames(filtered_matrix) <- names(filtered_rows)
colnames(filtered_matrix) <- colnames(InteractionMatrix_SMILES)
InteractionMatrix_SMILES <- filtered_matrix

## remove rows, that show only NAs
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[!apply(InteractionMatrix_SMILES, 1, function(x) all(is.na(x))), ]


####### chem x spec 
 desired_cols <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position",
                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position",
                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))

 desired_rows <- rownames(InteractionMatrix_SMILES)[!(rownames(InteractionMatrix_SMILES) %in% desired_cols)]
# #desired_cols <- c("Cl","C1=","CC","=C(","C=","C1)","C(","Cl)","C(Cl)(Cl)","C1=C(","=C1","O","C=C(","C=C1","OCC","CC1=","CC2","CO","P(=O)(","OC)","P(","=S)","(","CC(C)","OCC)","C(=O)O","C1","C(Cl)","N","CC=","=","CC(","S","=C", "2)" ,"C(C)","=N","[N+]([O-])=O)","C1=N","C(N")
#
 InteractionMatrix_SMILES <- InteractionMatrix_SMILES[rownames(InteractionMatrix_SMILES) %in% desired_rows, ]
 InteractionMatrix_SMILES <- InteractionMatrix_SMILES[,colnames(InteractionMatrix_SMILES) %in% desired_cols]
 
# ######## spec x chem 
 # desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position",
 #                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position",
 #                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
 # 
 # desired_cols <- colnames(InteractionMatrix_SMILES)[!(colnames(InteractionMatrix_SMILES) %in% desired_rows)]
 # # #desired_cols <- c("Cl","C1=","CC","=C(","C=","C1)","C(","Cl)","C(Cl)(Cl)","C1=C(","=C1","O","C=C(","C=C1","OCC","CC1=","CC2","CO","P(=O)(","OC)","P(","=S)","(","CC(C)","OCC)","C(=O)O","C1","C(Cl)","N","CC=","=","CC(","S","=C", "2)" ,"C(C)","=N","[N+]([O-])=O)","C1=N","C(N")
 # #
 # InteractionMatrix_SMILES <- InteractionMatrix_SMILES[rownames(InteractionMatrix_SMILES) %in% desired_rows, ]
 # InteractionMatrix_SMILES <- InteractionMatrix_SMILES[,colnames(InteractionMatrix_SMILES) %in% desired_cols]
 # 
 # 



 

pdf("heatmap_chemxspec.pdf", width = 10, height = 8) 
heatmap_chemxspec <- heatmap(
  InteractionMatrix_SMILES,
  main = "Trait Interactions",
  Colv = NA,
  Rowv = NA,
  #scale = "none",
  #col = rc,
  #margins = c(10, 10),
  # RowSideColors = rc,
  # scale = "row",
  # cexRow = 0.7,
  # cexCol = 0.7,
  las = 2,
  #labCol = colnames(filtered_matrix_reversed),
  #labRow = rownames(filtered_matrix_reversed),
)
dev.off()





pdf("/home/isabellehalbhuber/Toxicology/Results/heatmap_specxchem.pdf", width = 10, height = 8) 
heatmap_specxchem <- heatmap(
  InteractionMatrix_SMILES,
  main = "Trait Interactions",
  Colv = NA,
  Rowv = NA,
  #scale = "none",
  #col = rc,
  #margins = c(10, 10),
  # RowSideColors = rc,
  # scale = "row",
  # cexRow = 0.7,
  # cexCol = 0.7,
  las = 2,
  #labCol = colnames(filtered_matrix_reversed),
  #labRow = rownames(filtered_matrix_reversed),
)
dev.off()







library(reshape)
library(ggplot2)
library(RColorBrewer)
melted_data <- melt(InteractionMatrix_SMILES)

colnames(melted_data) <- c("Token", "Trait", "Value")

## only take 20 % highest values (80 Perzentil)
threshold <- quantile(unlist(melted_data$Value), 0.99)

df_filtered <- melted_data[melted_data$Value > threshold, ]
all_values <- unlist(lapply(df_filtered$Value, function(x) x))
global_min <- min(all_values)
global_max <- max(all_values)

labels_x <- levels(df_filtered$Trait)
breaks_x <- seq_along(labels_x)
labels_y <- levels(df_filtered$Token)
labels_y <- labels_y[labels_y %in% Tokens_99]


#labels_y <- labels_y[labels_y == Token]

labels_y <- labels_y[labels_y %in% (df_filtered$Token)[1:73]]

breaks_y <- seq_along(labels_y)

pdf("/home/isabellehalbhuber/Toxicology/Results/heatmap_99.pdf", width = 15, height = 20) 

p <- ggplot(df_filtered, aes(x = as.numeric(Trait), y = Token, fill = Value)) +
  geom_tile(color = "white") +
  scale_x_continuous(breaks = breaks_x, labels = labels_x) +
  #scale_y_continuous(breaks = breaks_y, labels = labels_y) +
  scale_fill_viridis_c(option = "viridis", name = "Interaction") +
  #scale_fill_viridis_c(option = "plasma", limits = c(global_min, global_max))+
  #scale_color_brewer(palette = "YlGnBu", limits = c(global_min, global_max)) +
  theme_linedraw()+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.line = element_line(linewidth = 0.5, colour = "grey80"),
    panel.grid.major.x = element_blank()
  )+
  #scale_y_continuous(breaks = breaks_y, labels = labels_y) +
  labs(
    x = NULL,
    y = NULL)
#   )+
# scale_fill_gradient(name = "Interaction",
#                      low = "#edf8b1",
#                      high = "#2c7fb8",
#                     na.value = "grey",
#                     limits = c(global_min, global_max))
# 
print(p) 

dev.off()


plot(NULL, 
     xlim = c(0.5, length(labels_x) + 0.5), 
     ylim = c(0.5, length(labels_y) + 0.5), 
     xaxt = 'n', yaxt = 'n', 
     xlab = "Trait", ylab = "Token",
     main = "Trait Interactions")

for (i in 1:nrow(df_filtered)) {
  x <- which(labels_x == df_filtered$Trait[i])
  y <- which(labels_y == df_filtered$Token[i])
  
  color_value <- (df_filtered$Value[i] - min(df_filtered$Value)) / (max(df_filtered$Value) - min(df_filtered$Value))
  fill_color <- rgb(1 - color_value, 1 - color_value, 1)  # einfache Farbskala

  rect(x - 0.5, y - 1, x + 0.5, y + 1, col = fill_color, border = "white")
}

axis(1, at = 1:length(labels_x), labels = labels_x, las = 2)
axis(2, at = 1:length(labels_y))

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
# Ich habe meine InteractonMatrix_SMILES, sie enthält alle summierten und dann normalisierten Einträge 
# 1. mache zwei neue Matrizen: species traits (rows) x chemical traits (cols) und chemicaltraits x chemical taits 
# 2. ich teile alle 424 Tokens in Gruppen chemisch sinnvolle Gruppen ein 
# 3. ich teile die Interaktionsstärken in 8 kategorien ein
# einen großen Plot, wo zu sehen ist, welche gruppen am meisten Interaktionen aufweisen 
# viele kleine Plots, wo nochmal auf die chemische Struktur geschaut werden kann 

matrix <- InteractionMatrix_SMILES

# order matrix 

desired_order <- paste0("trait_", 1:17)
existing_columns <- colnames(matrix)
traits_to_keep <- desired_order[desired_order %in% existing_columns]
other_columns <- existing_columns[!(existing_columns %in% traits_to_keep)]
new_column_order <- c(traits_to_keep, other_columns)
matrix <- matrix[, new_column_order, drop = FALSE]

desired_order_rows <- paste0("trait_", 1:17)
existing_rows <- rownames(matrix)
rows_to_keep <- desired_order_rows[desired_order_rows %in% existing_rows]
other_rows <- existing_rows[!(existing_rows %in% rows_to_keep)]
new_row_order <- c(rows_to_keep, other_rows)
matrix <- matrix[new_row_order, , drop = FALSE]

# give traits exact names

new_names <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", "Pectoral.Fin.Size",
               "Caudal.Peduncle.Throttling", "body_size", paste0("pEV_", 1:7))
colnames(matrix)[1:17] <- new_names
rownames(matrix)[1:17] <- new_names
# 
# 
# #############################
# ##### make first matrix
# 
desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position",
                  "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position",
                  "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))

desired_cols <- colnames(matrix)[!(colnames(matrix) %in% desired_rows)]
#desired_cols <- c("Cl","C1=","CC","=C(","C=","C1)","C(","Cl)","C(Cl)(Cl)","C1=C(","=C1","O","C=C(","C=C1","OCC","CC1=","CC2","CO","P(=O)(","OC)","P(","=S)","(","CC(C)","OCC)","C(=O)O","C1","C(Cl)","N","CC=","=","CC(","S","=C", "2)" ,"C(C)","=N","[N+]([O-])=O)","C1=N","C(N")

species_matrix_1 <- matrix[rownames(matrix) %in% desired_rows, ]
species_matrix <- species_matrix_1[,colnames(species_matrix_1) %in% desired_cols]

# #############################
# ##### make second matrix
# desired_rows_1 <- colnames(matrix)[!(colnames(matrix) %in% desired_rows)]
# desired_cols_1 <- colnames(matrix)[!(colnames(matrix) %in% desired_rows)]
# chemical_matrix_1 <- matrix[rownames(matrix) %in% desired_rows_1, ]
# chemical_matrix <- chemical_matrix_1[,colnames(chemical_matrix_1) %in% desired_cols_1]
# 
# ##############################
# ##### make groups 
# 
# 
# # "Cl", "Br", "F", "(Cl)", "Cl)", "Br)", "F)", "(F)", "C(Cl)(Cl)", "C(Cl)", "C(F)", 
# #"C(F)(F)", "C(F)(F)F", "C(F)(F)F)", "C(F)(F)C(F)(F)", "FC(F)(", "FC(F)(F)", "CCl)", 
# #"(C(F)(F)F)", "FC(", "CC(F)(F)F)", "(F)F)", "(F)(F)F)", "C(F)(",
# # Liste der funktionellen Gruppen und deren zugehörige Tokens
# group_patterns <- list(
#   
#   "Halogenated compounds" = c("FC(F)(F)","F", "C(Cl)" ,"C(F)(F)F","(F)","(C(F)(F)F)","C(F)F)","C(F)","Cl","Cl)", 
#                               "C(Cl)(Cl)","Br)","I)","CC(F)(F)F)","FC(","I","Br","(Cl)","FC(F)(","C(F)(F)C(F)(F)", 
#                               "(F)F)","(F)(F)F)","C(F)(", "C(F)(F)F)"),
#   
#   "Nitro and nitroso compounds" = c("[N+](=O)", "[N+](=O)", "[N+](=O)", "[N+](=O)", "[N+]([O-])=O", "[N+]([O-])=O)", "=N", "=N1",
#                                        "=N/", "N(C=O)", "N(C)=O", "N(C)C)", "N(C)C", "=N2", "=C(N)", "N1", "N(C=O)", 
#                                        "N2", "N2C(=O)", "N=C1", "N=C(N)", "N(C(=O)", "N(", "CCN(", "CNC(=O)", "NC1=N", "(N", "C(=N"),
#   
#   "Sulfur compounds" = c("S", "S(=O)", "SC", "SC)", "SC1", "SC)", "S(=O)(=O)", "S(=O)", "S(C)", "S)", "S1", 
#                              "SS", "S(C)(=O)=O", "S(=O)(=O)N", "CS", "CS)", "C(=S)", "C(=S)N", "C(=S)S", "CS(=O)(=O)", 
#                              "CS(=O)(=O)N", "=S)", "=S", "CCS", "CCCS"),
#   
#   "Carboxylic acids and derivatives" = c("C(=O)", "C(O)", "C(=O)O", "C(=O)OC", "C(=O)N", "C(O)=O", "C(=O)NC", "C(=O)N(C)", 
#                                   "C(O)C", "C(=O)C", "C1=O", "C(=O)N1", "C(=O)N2", "C(=O)NS(=O)(=O)", "C(=O)C1", "C(=O)C2", 
#                                   "C2=O", "O=C(", "C(=O)OC)", "C(=O)N(CC)", "O=C(N", "O=C1", "O=C1N", "O=C1N", "O=C1N1", "C1=O)", 
#                                   "O=C(N)", "C(=O)OCC", "O=C1", "OC(=O)", "OCO2)", "OC)", "CC(O)=O)", "CC(O)", "CC(C)O", "CC(O)(",
#                                   "O)", "C=O)", "OC1", "CC(O", "CCCO", "CO1", "CC(=O)OC", "CCOC1", "OCC3", "OC(", "CCOC(=O)", "O2)",
#                                   "CC2=O", "C(O)(", "CC(=O)", "=O)", "C(=O)C(", "C(C(O)=O)", "COC(=O)", "(O"),
#   
#   "Amine compounds" = c("N", "CN", "CN(", "CN(C)", "CN1", "CN2", "CN3", "CN(C)C", "N(C)", "NC", "N=C", "N(C)C)", "N(C)=O", 
#                          "N(CC)", "N(CC)CC", "C(N)=O", "NC(=S)", "NC(N)", "NCC", "NC(=N)N", "NC(=N)", "NC(=O)", "NC(=O)N", 
#                          "NC(=O)C", "NC(=O)C2", "N1C(=O)", "N1CCN(", "C(N)=O)", "N2", "NC(C=O)", "C(N)", "OCCN", "N=", "C(=O)N(C)C", "C(N",
#                         "CC(C)N", "CCC(=O)N", "C1CN(","N1)", "COC(=O)N", "C(=O)N(", "CN1C(=O)"),
#   
#   "Phosphorus compounds" = c("P(=O)", "P(=O)(O)", "P(=O)(OCC)", "P(=O)(O)", "P(=O)(O)=O)", "P(O)(O)", "P(O)(=O)", 
#                              "P(O)(O)=O", "P(=O)(OCC)", "P(=O)(", "P(", "P(=O)(O"),
#   
#   "Aliphatic compounds" = c("CC", "CCC", "CCCC", "CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC", 
#                             "CCCCCCCCCCC", "CC(C)", "CC(C(", "CC(C)C", "CCC(", "CCC1", "CCC(C)", 
#   "CCC(O)", "CCCCC1", "C(CC)", "CCC(C)(C)", "CC(CO", "CCC(O", "CCCC3", "CC(C)", "C(CC)", "CC(C)(C)", 
#   "CCO", "CC(", "C(C)C", "C(C(", "CC)", "CCCC)", "C(C)", "(C)C", "C(C)C)", "C(C)(C)C)", "(C)C)", "C(C)(C)C", "(CC)", "C)", "CCC)",
#   "CC(C)C)"),
#   
#   "Unsaturated hydrogens" = c("CC(C)=","CCCC=","CC=","C=C(", "C(C)=C(", "=", "=C", "=C)", "C=CC(", "=C(C)", "C(=", "CC=C)", "=C(", "C="),
#   
#   "Aromatic and heterocyclic compounds" = c("C1=", "C2=", "C3=", "C1C", "C1=C(", "C1=C", "C2=C(", "C2=O)", "C1=O", "C1=N", "=C1", 
#                                  "=C2", "C1(C)C", "C2=N", "C=C", "C=C1", "C1(C)", "C1CC2", "C1CO", "C1CCC(", "C1CC1", 
#                                  "C1CCCCC1", "C1CCCCC1)", "C2CC2", "C2(C)C", "C2C3", "C2C(", "C(C(=O)", "C(C)(C)", 
#                                  "C1CC1)", "C2=C", "C3=C(", "C1=N", "C1O", "C2=N", "C1", "C2", "C3", "C4", "N1", "N2", "N3", "C1C", "C2=C", "C1=O", "C2CC2","CC3", "CC2", 
#                            "C1CN", "C1CC2", "C1CCCCC1", "N1C(=O)", "C1CN", "C3=", "C2N", "N=C(N", "CCC2", "CC1=", "2", "C1(", "C2(",
#                            "CC1=C(", "CC1", "C12", "C3(", "C1C2", "=C3", "C1C(", "1)", "C=C2", "C1)", "2)"),
#   
#   "Alcohols and ethers" = c("O", "OH", "OC", "CO", "COC", "COC(", "COCC", "OCC", "OCC(", "OCC)", "OCO", "OCCO", "OCCOCC", 
#                             "OCCO)", "OCCCO", "OCCCC", "COH", "C(CO)", "CCO", "O="),
#   
#   "Cyan compounds" = c("C#N", "CN", "C=N", "C=N1", "C1N", "C(N)", "CCN"),
#   
#   "Other compounds" = c("[O-]", "[Si]", "[Se]", "(=O)", "(N)", "(O)", "=O", "O=", "(O)=O", "(O)=O)", "(CO)", 
#                         "/N", "\\", "[N+]1", "C1=O)", "[Si]", "[Se]", "[O-]", "C", "C(", "(")
# )
# 
# 
# 
# col_group <- rep(NA, ncol(species_matrix))
# 
# for (group in names(group_patterns)) {
#   matching_cols <- colnames(species_matrix) %in% group_patterns[[group]]
#   col_group[matching_cols] <- group
# }
# 
# 
# annotation_col <- data.frame(Group = factor(col_group))
# rownames(annotation_col) <- colnames(species_matrix)
# 
# na_rows <- rownames(annotation_col[is.na(annotation_col$Group), ])
# 
# 
# unique_tokens <- colnames(species_matrix)
# write.table(unique_tokens, file = "/home/isabellehalbhuber/Toxicology/Results/unique_tokens.txt")
# 
# 

################################
# library(dplyr)
# library(readr)
# clusteredTokens <- read.csv(file = "/home/isabellehalbhuber/Toxicology/Data/clusteredTokens.csv")
# clusteredTokens <- clusteredTokens[-1,]
### make distinct data matrices 
# min_val <- min(species_matrix, na.rm = TRUE)
# max_val <- max(species_matrix, na.rm = TRUE)
# 
# MiscellaneousToxicity <- clusteredTokens[clusteredTokens[,3] == "Miscellaneous Toxicity",]
# desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", 
#                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", 
#                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
# 
# desired_cols <- c(MiscellaneousToxicity$Token)
# 
# species_matrix_1 <- matrix[rownames(matrix) %in% desired_rows, ]
# MiscellaneousToxicity <- species_matrix_1[,colnames(species_matrix_1) %in% desired_cols]
# 
# library(heatmaply)
# heatmaply(species_matrix,
#           dendrogram = "none",
#           # xlab = "", ylab = "",
#           # main = "",
#           # scale = "column",
#           colors = viridis(n = 10, alpha = 1, begin = 0, end = 1, option = "plasma"),
#           #colors = brewer.pal(9,"YlGnBu"),
#           margins = c(60,100,40,20),
#           grid_color = "black",
#           grid_width = 0.0000001,
#           titleX = TRUE,
#           # hide_colorbar = FALSE,
#           branches_lwd = 0.1,
#           # #label_names = c("Country", "Feature:", "Value"),
#           fontsize_row = 10, fontsize_col = 10,
#           #labCol = colnames(new_matrix),
#           #labRow = rownames(new_matrix),
#           heatmap_layers = theme(axis.line = element_blank())
# )
# 
# 
# CarcinogenicityMutagenicity <- clusteredTokens[clusteredTokens[,3] == "Carcinogenicity and Mutagenicity",]
# 
# desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", 
#                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", 
#                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
# 
# desired_cols <- c(CarcinogenicityMutagenicity$Token)
# 
# species_matrix_1 <- matrix[rownames(matrix) %in% desired_rows, ]
# species_matrix <- species_matrix_1[,colnames(species_matrix_1) %in% desired_cols]
# 
# library(heatmaply)
# heatmaply(species_matrix,
#           dendrogram = "none",
#           # xlab = "", ylab = "",
#           # main = "",
#           # scale = "column",
#           colors = viridis(n = 10, alpha = 1, begin = 0, end = 1, option = "plasma"),
#           #colors = brewer.pal(9,"YlGnBu"),
#           margins = c(60,100,40,20),
#           grid_color = "black",
#           grid_width = 0.0000001,
#           titleX = TRUE,
#           # hide_colorbar = FALSE,
#           branches_lwd = 0.1,
#           # #label_names = c("Country", "Feature:", "Value"),
#           fontsize_row = 10, fontsize_col = 10,
#           #labCol = colnames(new_matrix),
#           #labRow = rownames(new_matrix),
#           heatmap_layers = theme(axis.line = element_blank())
# )
# 
# Lipophilicity <- clusteredTokens[clusteredTokens[,3] == "Lipophilicity (Hydrophobic)",]
# 
# desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", 
#                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", 
#                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
# 
# desired_cols <- c(Lipophilicity$Token)
# 
# species_matrix_1 <- matrix[rownames(matrix) %in% desired_rows, ]
# species_matrix <- species_matrix_1[,colnames(species_matrix_1) %in% desired_cols]
# 
# library(heatmaply)
# heatmaply(species_matrix,
#           dendrogram = "none",
#           # xlab = "", ylab = "",
#           # main = "",
#           # scale = "column",
#           colors = viridis(n = 10, alpha = 1, begin = 0, end = 1, option = "plasma"),
#           #colors = brewer.pal(9,"YlGnBu"),
#           margins = c(60,100,40,20),
#           grid_color = "black",
#           grid_width = 0.0000001,
#           titleX = TRUE,
#           # hide_colorbar = FALSE,
#           branches_lwd = 0.1,
#           # #label_names = c("Country", "Feature:", "Value"),
#           fontsize_row = 10, fontsize_col = 10,
#           #labCol = colnames(new_matrix),
#           #labRow = rownames(new_matrix),
#           heatmap_layers = theme(axis.line = element_blank())
# )
# 
# 
# UnsaturatedBonds <- clusteredTokens[clusteredTokens[,3] == "Unsaturated Bonds (Reactive)",]
# desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", 
#                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", 
#                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
# 
# desired_cols <- c(UnsaturatedBonds$Token)
# 
# species_matrix_1 <- matrix[rownames(matrix) %in% desired_rows, ]
# species_matrix <- species_matrix_1[,colnames(species_matrix_1) %in% desired_cols]
# 
# library(heatmaply)
# heatmaply(species_matrix,
#           dendrogram = "none",
#           # xlab = "", ylab = "",
#           # main = "",
#           # scale = "column",
#           colors = viridis(n = 10, alpha = 1, begin = 0, end = 1, option = "plasma"),
#           #colors = brewer.pal(9,"YlGnBu"),
#           margins = c(60,100,40,20),
#           grid_color = "black",
#           grid_width = 0.0000001,
#           titleX = TRUE,
#           # hide_colorbar = FALSE,
#           branches_lwd = 0.1,
#           # #label_names = c("Country", "Feature:", "Value"),
#           fontsize_row = 10, fontsize_col = 10,
#           #labCol = colnames(new_matrix),
#           #labRow = rownames(new_matrix),
#           heatmap_layers = theme(axis.line = element_blank())
# )
# 
# Hydrophilicity <- clusteredTokens[clusteredTokens[,3] == "Hydrophilicity (Hydrophilic)",]
# 
desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position",
                  "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position",
                  "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))

desired_cols <- c(Hydrophilicity$Token)

species_matrix_1 <- matrix[rownames(matrix) %in% desired_rows, ]
species_matrix <- species_matrix_1[,colnames(species_matrix_1) %in% desired_cols]
# 
# library(heatmaply)
# heatmaply(species_matrix,
#           dendrogram = "none",
#           # xlab = "", ylab = "",
#           # main = "",
#           # scale = "column",
#           colors = viridis(n = 10, alpha = 1, begin = 0, end = 1, option = "plasma"),
#           #colors = brewer.pal(9,"YlGnBu"),
#           margins = c(60,100,40,20),
#           grid_color = "black",
#           grid_width = 0.0000001,
#           titleX = TRUE,
#           # hide_colorbar = FALSE,
#           branches_lwd = 0.1,
#           # #label_names = c("Country", "Feature:", "Value"),
#           fontsize_row = 10, fontsize_col = 10,
#           #labCol = colnames(new_matrix),
#           #labRow = rownames(new_matrix),
#           heatmap_layers = theme(axis.line = element_blank())
# )
# 
# 
# ElectrophilicGroups <- clusteredTokens[clusteredTokens[,3] == "Electrophilic Functional Groups",]
# 
# desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", 
#                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", 
#                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
# 
# desired_cols <- c(ElectrophilicGroups$Token)
# 
# species_matrix_1 <- matrix[rownames(matrix) %in% desired_rows, ]
# species_matrix <- species_matrix_1[,colnames(species_matrix_1) %in% desired_cols]
# 
# library(heatmaply)
# heatmaply(species_matrix,
#           dendrogram = "none",
#           # xlab = "", ylab = "",
#           # main = "",
#           # scale = "column",
#           colors = viridis(n = 10, alpha = 1, begin = 0, end = 1, option = "plasma"),
#           #colors = brewer.pal(9,"YlGnBu"),
#           margins = c(60,100,40,20),
#           grid_color = "black",
#           grid_width = 0.0000001,
#           titleX = TRUE,
#           # hide_colorbar = FALSE,
#           branches_lwd = 0.1,
#           # #label_names = c("Country", "Feature:", "Value"),
#           fontsize_row = 10, fontsize_col = 10,
#           #labCol = colnames(new_matrix),
#           #labRow = rownames(new_matrix),
#           heatmap_layers = theme(axis.line = element_blank())
# )
# 
# RingStructures <- clusteredTokens[clusteredTokens[,3] == "Ring Structures",]
# desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", 
#                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", 
#                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
# 
# desired_cols <- c(RingStructures$Token)
# 
# species_matrix_1 <- matrix[rownames(matrix) %in% desired_rows, ]
# species_matrix <- species_matrix_1[,colnames(species_matrix_1) %in% desired_cols]
# 
# library(heatmaply)
# heatmaply(species_matrix,
#           dendrogram = "none",
#           # xlab = "", ylab = "",
#           # main = "",
#           # scale = "column",
#           colors = viridis(n = 10, alpha = 1, begin = 0, end = 1, option = "plasma"),
#           #colors = brewer.pal(9,"YlGnBu"),
#           margins = c(60,100,40,20),
#           grid_color = "black",
#           grid_width = 0.0000001,
#           titleX = TRUE,
#           # hide_colorbar = FALSE,
#           branches_lwd = 0.1,
#           # #label_names = c("Country", "Feature:", "Value"),
#           fontsize_row = 10, fontsize_col = 10,
#           #labCol = colnames(new_matrix),
#           #labRow = rownames(new_matrix),
#           heatmap_layers = theme(axis.line = element_blank())
# )
# 
# 
# Halogens <- clusteredTokens[clusteredTokens[,3] == "Halogens",]
# 
# desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", 
#                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", 
#                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
# 
# desired_cols <- c(Halogens$Token)
# 
# species_matrix_1 <- matrix[rownames(matrix) %in% desired_rows, ]
# species_matrix <- species_matrix_1[,colnames(species_matrix_1) %in% desired_cols]
# 
# library(heatmaply)
# heatmaply(species_matrix,
#           dendrogram = "none",
#           # xlab = "", ylab = "",
#           # main = "",
#           # scale = "column",
#           colors = viridis(n = 10, alpha = 1, begin = 0, end = 1, option = "plasma"),
#           #colors = brewer.pal(9,"YlGnBu"),
#           margins = c(60,100,40,20),
#           grid_color = "black",
#           grid_width = 0.0000001,
#           titleX = TRUE,
#           # hide_colorbar = FALSE,
#           branches_lwd = 0.1,
#           # #label_names = c("Country", "Feature:", "Value"),
#           fontsize_row = 10, fontsize_col = 10,
#           #labCol = colnames(new_matrix),
#           #labRow = rownames(new_matrix),
#           heatmap_layers = theme(axis.line = element_blank())
# )
# 
# 
# 
# 
# 
# ##################################
# library(heatmaply)
# library(RColorBrewer)
# 
# 
# desired_rows <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", 
#                   "Oral.Gape.Position", "Relative.Maxillary.Length", "Body.Lateral.Shape", 
#                   "Pectoral.Fin.vertical.Position", "Pectoral.Fin.Size", 
#                   "Caudal.Peduncle.Throttling", "body_size", paste0("pEV_", 1:7))
# 
# categories <- c(unique(clusteredTokens$Toxicity.Cluster))
# 
# 
# species_matrices <- list()
# 
# for (category in categories) {
#   category_data <- clusteredTokens[clusteredTokens[,3] == category,]
#   
# 
#   if (nrow(category_data) > 0) {
# 
#     desired_cols <- category_data$Token
#     
#     species_matrix_1 <- matrix[rownames(matrix) %in% desired_rows, ]
#     species_matrix <- species_matrix_1[, colnames(species_matrix_1) %in% desired_cols]
#     
#     if (ncol(species_matrix) > 0) {
#       species_matrices[[category]] <- species_matrix
#     } 
#   } 
# }
# 
# 
# pdf("heatmaps.pdf", width = 11, height = 8.5)  
# 
# for (category in names(species_matrices)) {
#   species_matrix <- species_matrices[[category]]
#   
#   heatmap(
#     species_matrix)
# }
# 
# dev.off() 
############################################################################################
############################ Plot chemicals x traits interactions ##########################
############################################################################################
library(readr)
library(ggplot2)
library(reshape2)
library(viridis)
library(dplyr)
library(RColorBrewer)
#TokensLookUpTable <- read.csv(file = "/home/isabellehalbhuber/Toxicology/Results/TokensLookUpTable.csv")
#TokensLookUpTable$int <- as.integer(TokensLookUpTable$int +1)

#InteractionMatrix_SMILES <- readRDS(file = "/home/isabellehalbhuber/Toxicology/Results/InteractionMatrix_SMILES.rds")

# order matrix
# 
# desired_order <- paste0("trait_", 1:17)
# existing_columns <- colnames(InteractionMatrix_SMILES)
# traits_to_keep <- desired_order[desired_order %in% existing_columns]
# other_columns <- existing_columns[!(existing_columns %in% traits_to_keep)]
# new_column_order <- c(traits_to_keep, other_columns)
# InteractionMatrix_SMILES <- InteractionMatrix_SMILES[, new_column_order, drop = FALSE]
# 
# desired_order_rows <- paste0("trait_", 1:17)
# existing_rows <- rownames(InteractionMatrix_SMILES)
# rows_to_keep <- desired_order_rows[desired_order_rows %in% existing_rows]
# other_rows <- existing_rows[!(existing_rows %in% rows_to_keep)]
# new_row_order <- c(rows_to_keep, other_rows)
# InteractionMatrix_SMILES <- InteractionMatrix_SMILES[new_row_order, , drop = FALSE]
# 
# new_names <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", "Pectoral.Fin.Size",
#                "Caudal.Peduncle.Throttling", "body_size", paste0("pEV_", 1:7))
# colnames(InteractionMatrix_SMILES)[1:17] <- new_names
# rownames(InteractionMatrix_SMILES)[1:17] <- new_names
# 
# 
# desired_cols <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position",
#                   "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position",
#                   "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))
# 
# desired_rows <- rownames(InteractionMatrix_SMILES)[!(rownames(InteractionMatrix_SMILES) %in% desired_cols)]
# InteractionMatrix_SMILES <- InteractionMatrix_SMILES[rownames(InteractionMatrix_SMILES) %in% desired_rows, ]
# InteractionMatrix_SMILES <- InteractionMatrix_SMILES[,colnames(InteractionMatrix_SMILES) %in% desired_cols]
# matrix <- InteractionMatrix_SMILES
clusteredTokens <- read.csv(file = "/home/isabellehalbhuber/Toxicology/Data/clusteredTokens.csv")
clusteredTokens <- clusteredTokens[-1,-1]



## only take 20 % highest values (99 Perzentil)
melted_data <- melt(InteractionMatrix_SMILES)
colnames(melted_data) <- c("Token", "Trait", "Value")
threshold <- quantile(unlist(melted_data$Value), 0.99)
df_filtered <- melted_data[melted_data$Value > threshold, ]
df_filtered_matrix <- as.matrix(df_filtered)

#saveRDS(df_filtered_matrix, file = "/home/isabellehalbhuber/Toxicology/Results/99PercentilTokens.rds")


all_values <- unlist(lapply(df_filtered$Value, function(x) x))
global_min <- min(all_values)
global_max <- max(all_values)

clusteredTokens <- as.data.frame(df_filtered_matrix)
clusteredTokens$Toxicity.Cluster <- c("Interactions")
clusteredTokens <- as.data.frame(clusteredTokens)

## all Tokens
clusteredTokens <- as.data.frame(TokensLookUpTable$token)
clusteredTokens$Toxicity.Cluster <- c("Interactions")
clusteredTokens <- as.data.frame(clusteredTokens)


# hier weiter

desired_cols <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", 
                  "Oral.Gape.Position", "Relative.Maxillary.Length", "Body.Lateral.Shape", 
                  "Pectoral.Fin.vertical.Position", "Pectoral.Fin.Size", 
                  "Caudal.Peduncle.Throttling", "body_size")#, paste0("pEV_", 1:7))


categories <- c("Miscellaneous Toxicity", 
                "Carcinogenicity and Mutagenicity", 
                "Lipophilicity (Hydrophobic)", 
                "Unsaturated Bonds (Reactive)", 
                "Hydrophilicity (Hydrophilic)", 
                "Electrophilic Functional Groups", 
                "Ring Structures", 
                "Halogens",
                "Interactions")
#matrix <- InteractionMatrix_SMILES

matrix <- InteractionMatrix_SMILES
species_matrices <- list()
for (category in categories) {
  category_data <- clusteredTokens[clusteredTokens[,4] == category,] # change colums 2 for all, 4 
  
  
  if (nrow(category_data) > 0) {
    
    desired_rows <- category_data$Token
    
    species_matrix_1 <- matrix[colnames(matrix) %in% desired_cols, ]
    species_matrix <- species_matrix_1[rownames(species_matrix_1) %in% desired_rows, ]
    
    if (ncol(species_matrix) > 0) {
      species_matrices[[category]] <- species_matrix
    } 
  } 
}

all_values <- unlist(lapply(species_matrices, function(x) x))
global_min <- min(all_values)
global_max <- max(all_values)

#my_colors <- c("#FFFFFF", "#A0DA39", "#4AC16D", "#1FA187", "#277F8E", "#365C8D", "#46327E", "#440154")
my_colors <- c("#FFFFFF","#ffffd9", "#edf8b1","#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58")



### plot for 0.99 perzentil
pdf("heatmaps_ggplot_test_SAMPLE.pdf", width = 25, height = 20)

for (category in names(species_matrices)) {
  species_matrix <- species_matrices[[category]]

  melted_data <- melt(species_matrix)
  colnames(melted_data) <- c("Token", "Trait", "Value")

  p <- ggplot(melted_data, aes(x = Trait, y = Token, fill = Value)) +
    geom_tile(color = "lightgrey") +
    scale_fill_gradientn(colors = my_colors, limits = c(global_min, global_max), name = "Interaction") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 25 ),
      axis.text.y = element_text(size = 20),
      plot.title = element_text(hjust = 0.5, size = 30),
      panel.grid = element_blank()
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = paste("Token - Trait", category)
    )

  print(p)
}

dev.off()

### plot for all tokens
pdf("heatmaps_ggplot_all.pdf", width = 25, height = 20)

for (category in names(species_matrices)) {
  species_matrix <- species_matrices[[category]]
  
  melted_data <- melt(species_matrix)
  colnames(melted_data) <- c("Token", "Trait", "Value")
  
  p <- ggplot(melted_data, aes(x = Trait, y = Token, fill = Value)) +
    geom_tile(color = "#F0F8FF") +
    scale_fill_gradientn(colors = my_colors, limits = c(global_min, global_max), name = "Interaction") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 25 ),
      axis.text.y = element_text(size = 3),
      plot.title = element_text(hjust = 0.5, size = 30),
      panel.grid = element_blank()
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = paste("Token - Trait", category)
    )
  
  print(p)
}

dev.off()

############################################################################################
########################## Plot chemicals x chemicals interactions #########################
############################################################################################
InteractionMatrix_SMILES <- readRDS(file = "/home/isabellehalbhuber/Toxicology/Results/InteractionMatrix_SMILES.rds")

# # order matrix 
# 
desired_order <- paste0("trait_", 1:17)
existing_columns <- colnames(InteractionMatrix_SMILES)
traits_to_keep <- desired_order[desired_order %in% existing_columns]
other_columns <- existing_columns[!(existing_columns %in% traits_to_keep)]
new_column_order <- c(traits_to_keep, other_columns)
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[, new_column_order, drop = FALSE]

desired_order_rows <- paste0("trait_", 1:17)
existing_rows <- rownames(InteractionMatrix_SMILES)
rows_to_keep <- desired_order_rows[desired_order_rows %in% existing_rows]
other_rows <- existing_rows[!(existing_rows %in% rows_to_keep)]
new_row_order <- c(rows_to_keep, other_rows)
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[new_row_order, , drop = FALSE]

new_names <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position", "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position", "Pectoral.Fin.Size",
               "Caudal.Peduncle.Throttling", "body_size", paste0("pEV_", 1:7))
colnames(InteractionMatrix_SMILES)[1:17] <- new_names
rownames(InteractionMatrix_SMILES)[1:17] <- new_names


excluded_cols <- c("Body.Elongation", "Vertical.Eye.Position", "Relative.Eye.Size", "Oral.Gape.Position",
                  "Relative.Maxillary.Length", "Body.Lateral.Shape", "Pectoral.Fin.vertical.Position",
                  "Pectoral.Fin.Size", "Caudal.Peduncle.Throttling", "body_size",paste0("pEV_", 1:7))

desired_RowsCols<- rownames(InteractionMatrix_SMILES)[!(rownames(InteractionMatrix_SMILES) %in% excluded_cols)]
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[rownames(InteractionMatrix_SMILES) %in% desired_RowsCols, ]
InteractionMatrix_SMILES <- InteractionMatrix_SMILES[,colnames(InteractionMatrix_SMILES) %in% desired_RowsCols]
matrix <- InteractionMatrix_SMILES

### 99 Percentil tokens
df_filtered_matrix <- readRDS(file = "/home/isabellehalbhuber/Toxicology/Results/99PercentilTokens.rds")

clusteredTokens <- as.data.frame(df_filtered_matrix)
clusteredTokens$Toxicity.Cluster <- c("Interactions")
clusteredTokens <- as.data.frame(clusteredTokens)

## all tokens
clusteredTokens <- as.data.frame(TokensLookUpTable$token)
clusteredTokens$Toxicity.Cluster <- c("Interactions")
clusteredTokens <- as.data.frame(clusteredTokens)

#### weiter 
categories <- c("Miscellaneous Toxicity", 
                "Carcinogenicity and Mutagenicity", 
                "Lipophilicity (Hydrophobic)", 
                "Unsaturated Bonds (Reactive)", 
                "Hydrophilicity (Hydrophilic)", 
                "Electrophilic Functional Groups", 
                "Ring Structures", 
                "Halogens",
                "Interactions")

matrices <- list()
for (category in categories) {
  category_data <- clusteredTokens[clusteredTokens[,4] == category,] 
  
  
  if (nrow(category_data) > 0) {
    
    desired_rows <- category_data$Token
    
     matrix <- matrix[colnames(matrix) %in% desired_rows, ]
     matrix <- matrix[rownames(matrix) %in% desired_rows, ]
    
    if (ncol(matrix) > 0) {
      matrices[[category]] <- matrix
    } 
  } 
}

all_values <- unlist(lapply(matrices, function(x) x))
global_min <- min(all_values)
global_max <- max(all_values)

#my_colors <- c("#FFFFFF", "#A0DA39","#4AC16D","#1FA187","#277F8E","#365C8D","#46327E","#440154")
my_colors <- c("#FFFFFF","#ffffd9", "#edf8b1","#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58", "black")



### plot for all token
pdf("heatmaps_ggplot_99_chemicals.pdf", width = 25, height = 20)

for (category in names(matrices)) {
  matrix <- matrices[[category]]
  
  melted_data <- melt(matrix)
  colnames(melted_data) <- c("Token1", "Token2", "Value")
  
  p <- ggplot(melted_data, aes(x = Token2, y = Token1, fill = Value)) +
    geom_tile(color = "lightgrey") +
    scale_fill_gradientn(colors = my_colors, limits = c(global_min, global_max), name = "Interaction") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 3),
      axis.text.y = element_text(size = 3),
      plot.title = element_text(hjust = 0.5, size = 30),
      panel.grid = element_blank()
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = paste("Token - Token", category)
    )
  
  print(p)
}

dev.off()
