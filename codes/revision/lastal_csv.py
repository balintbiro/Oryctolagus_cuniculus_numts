import numpy as np
import pandas as pd

df=pd.read_csv('../data/revision/chr_scf_alignments.tab',skiprows=30,sep='\t')
df.columns=['score','chr_id','chr_start','chr_no_letters','chr_strand','chr_size','scf_id','scf_start','scf_no_letters','scf_strand','scf_size','blocks','EG2_val','E-val']

df=df.dropna()

numts=pd.read_csv('../results/3numt_array.csv')

def chr_spec_dfs(chr_id):
	global df_chunks
	alignment_subdf=df.loc[df['chr_id']==chr_id]
	numt_subdf=numts.loc[numts['g_id']==str(int(chr_id))]
	numt_ranges=np.concatenate(numt_subdf.apply(lambda row:np.arange(row['g_start'],(row['g_start']+row['g_length'])),axis=1).tolist())
	chr_fil=alignment_subdf['chr_start'].apply(lambda chr_start:chr_start in numt_ranges)
	df_chunks.append(df_chunks.append(alignment_subdf[chr_fil]))

df_chunks=[]
pd.Series(df['chr_id'].unique()).apply(chr_spec_dfs)

filtered_df=pd.concat(df_chunks)

filtered_df.to_csv('../data/revision/filtered_alignments.csv')