#### Bi-GRU model ####
#### NO TRAITS 
BiGRU_NoTRAITS <- nn_module(
  initialize = function(embedding = 40L, rnn_hidden_size = 20L, rnn_layers = 3L, dropout_prob = 0.3) { 
    self$embedding <- nn_embedding(615L, 
                                   embedding_dim = embedding,
                                   sparse = FALSE,
                                   scale_grad_by_freq = TRUE)
    
    self$rnn <- nn_gru(input_size = embedding, 
                       hidden_size = rnn_hidden_size,
                       num_layers = rnn_layers,  
                       batch_first = TRUE,
                       dropout = ifelse(rnn_layers > 1, dropout_prob, 0), 
                       bidirectional = TRUE)
    
    self$output_layer <- nn_linear(in_features = 2 * rnn_hidden_size, out_features = 1L)
    self$dropout_output <- nn_dropout(p = dropout_prob) 
  },
  
  forward = function(x) {
    x <- self$embedding(x)
    rnn_out <- self$rnn(x)[[1]][, dim(x)[2], ]  
    rnn_out <- self$dropout_output(rnn_out)  
    output <- self$output_layer(rnn_out)
    
    return(output)
  }
)

#### TRAITS ####

BiGRU_TRAITS <- nn_module(
  initialize = function(embedding = 1L, rnn_hidden_size = 20L, rnn_layers = 3L, dropout_prob = 0.3) {
    self$embedding <- nn_embedding(num_embeddings = 425L, 
                                   embedding_dim = embedding, 
                                   sparse = FALSE, 
                                   scale_grad_by_freq = TRUE) # [batchsize, time, embedding]
    
    # RNN with Dropout
    self$rnn <- nn_gru(input_size = embedding, 
                       hidden_size = rnn_hidden_size, 
                       num_layers = rnn_layers, 
                       bidirectional = TRUE,
                       dropout = dropout_prob,  # Dropout between RNN layers
                       batch_first = TRUE)
    
    # Token NN after RNN
    self$token_full = nn_linear(in_features = 2*rnn_hidden_size, out_features = 20L)
    self$dropout_token = nn_dropout(p = dropout_prob)  # Dropout after linear layer
    
    # Species Trait NN
    self$trait_full1 = nn_linear(in_features = 17L, out_features = 20L)
    self$dropout_trait1 = nn_dropout(p = dropout_prob)  # Dropout after linear layer
    self$trait_full2 = nn_linear(in_features = 20L, out_features = 3L)
    
    # Combined output layers
    self$output_layer1 <- nn_linear(in_features = 20L+3L, out_features = 50L)
    self$dropout_output1 = nn_dropout(p = dropout_prob)  # Dropout after linear layer
    self$output_layer2 <- nn_linear(in_features = 50L, out_features = 1L)
  },
  
  forward = function(tokens, traits) {
    #### Token RNN
    x <- self$embedding(tokens)
    x <- self$rnn(x)[[1]]  # only the first element, which is the output
    output_token <- self$token_full(x[, dim(x)[2], ]) %>% nnf_selu()
    output_token <- self$dropout_token(output_token)  # Apply dropout
    
    #### Trait NN
    x <- self$trait_full1(traits) %>% nnf_selu()
    x <- self$dropout_trait1(x)  # Apply dropout
    output_traits <- self$trait_full2(x) %>% nnf_selu()  # [batch_size, 3]
    
    #### Combine outputs
    input <- torch::torch_cat(list(output_token, output_traits), dim = 2L)  # [batchsize, rnn_hidden_size + 3L]
    output <- self$output_layer1(input) %>% nnf_selu()
    output <- self$dropout_output1(output)  # Apply dropout
    output <- self$output_layer2(output)
    
    return(output)
  }
)