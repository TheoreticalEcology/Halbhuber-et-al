library(readr)

# RF with
Attention_GloVe <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_Attention_with_GloVe.rds")
Attention_GloVe$Model = "GloVe"
Attention_GloVe$CV = rownames(Attention_GloVe)[1:3] 
Attention_noPEV <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_Attention_with_noPEV.rds")
Attention_noPEV$Model = "noPEV"
Attention_noPEV$CV = rownames(Attention_noPEV)[1:3] 
Attention_noRNN <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_Attention_with_noRNN.rds")
Attention_noRNN$Model = "noRNN"
Attention_noRNN$CV = rownames(Attention_noRNN)[1:3] 
Attention_sample <- readRDS("/home/isabellehalbhuber/Toxicology/Results/results_r_squared_Attention_with_sample_2.rds")
Attention_sample$Model = "sample"
Attention_sample$CV = rownames(Attention_sample)[1:3] 

Attention_trials <- rbind(Attention_GloVe, Attention_noPEV, Attention_noRNN, Attention_sample)

write_csv(Attention_trials, file = "/home/isabellehalbhuber/Toxicology/Results/Attention_trials.csv")
