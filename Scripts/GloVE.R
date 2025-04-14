#################################################################################
########################  Pre-Train Embeddings   ################################
#################################################################################
#install.packages("text2vec")
library(text2vec)
Tokens <- readRDS("/home/isabellehalbhuber/Toxicology/Data/Tokens.rds")
Tokens <- Tokens[,2:178]
Tokens[is.na(Tokens)] <- 0
Tokens <- rbind(Tokens)
Tokens <- apply(Tokens, 1, paste, collapse = " ")
Tokens <- space_tokenizer(Tokens)
it <- itoken(Tokens, progressbar = FALSE)
vocab <- create_vocabulary(it)
vocab <- prune_vocabulary(vocab, term_count_min = 1L)
vectorizer <- vocab_vectorizer(vocab)

tcm_example <- create_tcm(it, vectorizer, skip_grams_window = 177L)
glove = GlobalVectors$new(rank =5, x_max = 10, learning_rate = 0.001)
wv_main = glove$fit_transform(tcm_example, n_iter = 1000, convergence_tol = 0.0001, n_threads = 8)
wv_context = glove$components
word_vectors = wv_main + t(wv_context)

saveRDS(word_vectors, file = "/home/isabellehalbhuber/Toxicology/Data/token_distanceVectors.rds") 

d = dist2(wv_main, method="cosine")  #Smaller values means closer
print(dim(d))



# Tokens <- readRDS("/home/isabellehalbhuber/Toxicology/Scripts/Tokens.rds")
# Tokens <- Tokens[,2:178]
# Tokens <- na.omit(unlist(Tokens))
# Tokens <- list(c(Tokens))
# it <- itoken(Tokens, progressbar = FALSE)
# vocab <- create_vocabulary(it)
# vocab <- prune_vocabulary(vocab)
# vectorizer <- vocab_vectorizer(vocab)
# 
# tcm_example <- create_tcm(it, vectorizer, skip_grams_window = 177L)
# glove = GlobalVectors$new(rank = 50, x_max = 100, learning_rate = 0.0001, alpha = 0.075, lambda = 0.01)
# wv_main = glove$fit_transform(tcm, n_iter = 50, convergence_tol = 0.1, n_threads = 3)

#install.packages("text2vec")
library(text2vec)
Tokens <- readRDS("/home/isabellehalbhuber/Toxicology/Data/atomwiseTokens.rds")
Tokens <- Tokens[,2:257]
Tokens[is.na(Tokens)] <- 0
Tokens <- rbind(Tokens)
Tokens <- apply(Tokens, 1, paste, collapse = " ")
Tokens <- space_tokenizer(Tokens)
it <- itoken(Tokens, progressbar = FALSE)
vocab <- create_vocabulary(it)
vocab <- prune_vocabulary(vocab, term_count_min = 2L)
vectorizer <- vocab_vectorizer(vocab)
tcm_example <- create_tcm(it, vectorizer, skip_grams_window =200L)
glove = GlobalVectors$new(rank =5, x_max = 10, learning_rate = 0.0001)
wv_main = glove$fit_transform(tcm_example, n_iter = 1000, convergence_tol = 0.0001, n_threads = 8)
wv_context = glove$components
word_vectors = wv_main + t(wv_context)

saveRDS(word_vectors, file = "/home/isabellehalbhuber/Toxicology/Data/token_distanceVectors.rds") 

d = dist2(wv_main, method="cosine")  #Smaller values means closer
print(dim(d))


