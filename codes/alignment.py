import os
import pandas as pd
from subprocess import call

numts=pd.read_csv('../../results/3numt_array.csv')

fil=numts['g_id'].apply(lambda g_id:len(g_id)<3)
chr_ids=pd.Series(numts[fil]['g_id'].unique())

def alignment(chr_id):
	call(f'lastdb db {chr_id}.fa',shell=True)
	call(f'lastal -r1 -q1 -a7 -b1 db scaffolds.fasta > alignments/{chr_id}.afa',shell=True)
	call('rm db.bck db.des db.prj db.sds db.ssp db.suf db.tis',shell=True)

chr_ids.apply(alignment)