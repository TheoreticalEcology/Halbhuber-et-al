from SmilesPE.pretokenizer import atomwise_tokenizer

## Basic Tokenizers

# Atom-level Tokenizer
smi = 'CC[N+](C)(C)Cc1ccccc1Br'
toks = atomwise_tokenizer(smi)
print(toks)
from SmilesPE.pretokenizer import kmer_tokenizer

#K-mer Tokenzier
smi = 'CC[N+](C)(C)Cc1ccccc1Br'
toks = kmer_tokenizer(smi, ngram=4)
print(toks)


# ####### Use the Pre-trained SmilesPE Tokenizer

# import codecs
# from SmilesPE.tokenizer import SPE_Tokenizer

# # open voc-file
# with codecs.open('/Users/isi/Desktop/SMILESpe_tokenize/SPE_ChEMBL.txt', 'r', 'utf-8') as spe_vob:
#     spe = SPE_Tokenizer(spe_vob)

# # read SMILES
# smi_file_path = '/Users/isi/Desktop/SMILESpe_tokenize/SMILES_EPA_DSS.smi'

# with codecs.open(smi_file_path, 'r', 'utf-8') as smi_file:
#     smi_content = smi_file.readlines()

# # tokenization
# toks = [spe.tokenize(smiles.strip()) for smiles in smi_content]

# # write tokenized files
# file_path = '/Users/isi/Desktop/SMILESpe_tokenize/tokenizedSMI.txt'

# with codecs.open(file_path, 'w', 'utf-8') as file:
#     for tokenized_smiles in toks:
#         file.write(' '.join(tokenized_smiles) + '\n')


######################## generate df
import codecs
import pandas as pd
from SmilesPE.tokenizer import SPE_Tokenizer

# load voc
vocab_file_path = '/Users/isi/Desktop/SMILESpe_tokenize/SPE_ChEMBL.txt'
with codecs.open(vocab_file_path, 'r', 'utf-8') as spe_vob:
    spe = SPE_Tokenizer(spe_vob)

# read SMILES
smi_file_path = '/Users/isi/Desktop/SMILESpe_tokenize/SMILES_EPA_DSS.smi'
with codecs.open(smi_file_path, 'r', 'utf-8') as smi_file:
    smi_content = smi_file.readlines()

# tokenize smiles
toks = [spe.tokenize(smiles.strip()) for smiles in smi_content]

# generate DF
df = pd.DataFrame({
    'Originaler SMILES-Code': [smiles.strip() for smiles in smi_content],
    'Tokenisierte SMILES-Codes': [','.join(tokens) for tokens in toks]
})

# save df 
df.to_csv('/Users/isi/Desktop/SMILESpe_tokenize/smiles_tokenized_dataframe.csv', index=False, encoding='utf-8')

