BiRNNattn <- nn_module(
  initialize = function(embedding = 40L, rnn_hidden_size = 20L, rnn_layers = 3L, max_len = 177L) {
    self$embedding <- nn_embedding(
      num_embeddings = 425L, 
      embedding_dim = embedding, 
      sparse = FALSE, 
      scale_grad_by_freq = TRUE
    )
    
    self$rnn <- nn_gru(
      input_size = embedding,
      hidden_size = rnn_hidden_size,
      num_layers = rnn_layers,
      bidirectional = TRUE,
      batch_first = TRUE,
      dropout = 0 
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
  
  forward = function(tokens, traits) {
    #### Token Embedding
    x <- self$embedding(tokens)

    #### RNN Output
    rnn_out <- self$rnn(x)[[1]]  # [batch_size, seq_len, 2 * rnn_hidden_size]

    
    #### Trait Embedding
    traits = traits$unsqueeze(3L)
    traits_position = torch_tensor(matrix(scale(1:17), 17L, 1L),  # skalieren nicht vergessen
                                   device=traits$device,
                                   dtype=traits$dtype)$`repeat`(list(tokens$shape[1], 1L, 1L))
    traits = torch_cat(list(traits, traits_position), dim = 3L)
    
    traits_embedded = self$trait_full1(traits) %>% nnf_relu() %>% self$dropout1()
    traits_embedded = self$trait_full2(traits_embedded) %>% nnf_relu() %>% self$dropout2()
    traits_embedded = self$trait_full3(traits_embedded)

    
    
    #### Concatenation of Token and Trait Embeddings
    rnn_out_traits = torch_cat(list(rnn_out, traits_embedded), dim = 2L)

    # QKV for Self-Attention
    q <- self$q_linear(rnn_out_traits)  # [batch_size, seq_len, qkv_size]
    k <- self$k_linear(rnn_out_traits)  # [batch_size, seq_len, qkv_size]
    v <- self$v_linear(rnn_out_traits)  # [batch_size, seq_len, qkv_size]
    
    #### Masking for Padding Tokens
    padding_mask <- tokens$clone()
    mask <- padding_mask == 1L  # Assuming padding tokens are marked with 1 [batchsize, 177]
    mask = torch_cat(list(mask, torch_zeros(x$shape[1], traits$shape[2], dtype=torch_bool(), device=x$device)), dim = 2)
    
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
