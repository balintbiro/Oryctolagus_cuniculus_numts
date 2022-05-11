import os
import pandas as pd
from subprocess import call

numts=pd.read_csv('../results/3numt_array.csv')

chr_fil=numts['g_id'].apply(lambda g_id:len(g_id)<3)

chr_ids=pd.Series(numts[chr_fil]['g_id'].unique())

chr_ids.apply(
	lambda chr_id: call(f'samtools faidx ../data/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa {chr_id} >> ../data/revision/chr_sequences.fasta',shell=True)
	)

scf_ids=pd.Series(numts[~chr_fil]['g_id'].unique())

scf_ids.apply(
	lambda scf_id: call(f'samtools faidx ../data/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa {scf_id} >> ../data/revision/scf_sequences.fasta',shell=True)
	)

call('cd ../data/revision/',shell=True)
call('lastdb db chr_sequences.fasta', shell=True)
call('lastal -r1 -q1 -a7 -b1 db scf_sequences.fasta > chr_scf_alignments.afa',shell=True)