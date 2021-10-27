#importing the required modules
import numpy as np
import pandas as pd

#reading the numt csv into a df
numts=pd.read_csv('../results/3numt_array.csv')

#size for the regions where we would like to check the level of conservation
size=5

#writing 5 bp upstream and downstream numt sequences into a fasta file
numt_ends=numts['g_sequence'].apply(lambda sequence: sequence[:size].upper()+sequence[-size:].upper())
numt_ends.index=np.arange(0,len(numt_ends))
with open('../results/numt_g_ends.fa','w')as outfile:
    numt_ends.apply(lambda seq: outfile.write(seq+'\n'))