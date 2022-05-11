#importing the modules
import numpy as np
import pandas as pd
from Bio import pairwise2
from subprocess import run

#read in the filtered alignments
alignments=pd.read_csv('../data/revision/filtered_alignments.csv',index_col=0)

#function for aligning the sequences
def pairwise_alignments(row):
	try:
		global scores
		chr_id=row['chr_id']
		chr_start=row['chr_start']
		chr_end=row['chr_start']+row['chr_no_letters']
		chr_part_size=row['chr_no_letters']

		scf_id=row['scf_id']
		scf_start=row['scf_start']
		scf_end=row['scf_start']+row['scf_no_letters']
		scf_part_size=row['scf_no_letters']

		numtChr=''.join(str(run(f'samtools faidx ../data/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa {str(int(chr_id))}:{str(int(chr_start))}-{str(int(chr_end))}',shell=True,capture_output=True).stdout).split('\\n')[1:-1])
		numtScf=''.join(str(run(f'samtools faidx ../data/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa {scf_id}:{str(int(scf_start))}-{str(int(scf_end))}',shell=True,capture_output=True).stdout).split('\\n')[1:-1])

		upstreamChr=''.join(str(run(f'samtools faidx ../data/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa {str(int(chr_id))}:{str(int(chr_start-chr_part_size))}-{str(int(chr_start))}',shell=True,capture_output=True).stdout).split('\\n')[1:-1])
		upstreamScf=''.join(str(run(f'samtools faidx ../data/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa {scf_id}:{str(int(scf_start-scf_part_size))}-{str(int(scf_start))}',shell=True,capture_output=True).stdout).split('\\n')[1:-1])

		downstreamChr=''.join(str(run(f'samtools faidx ../data/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa {str(int(chr_id))}:{str(int(chr_end))}-{str(int(chr_end+chr_part_size))}',shell=True,capture_output=True).stdout).split('\\n')[1:-1])
		downstreamScf=''.join(str(run(f'samtools faidx ../data/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa {scf_id}:{str(int(scf_end))}-{str(int(scf_end+scf_part_size))}',shell=True,capture_output=True).stdout).split('\\n')[1:-1])

		upstream_score=pairwise2.align.globalms(upstreamChr,upstreamScf,1,-1,-7,-1,one_alignment_only=True)[0].score
		numt_score=pairwise2.align.globalms(numtChr,numtScf,1,-1,-7,-1,one_alignment_only=True)[0].score
		downstream_score=pairwise2.align.globalms(downstreamChr,downstreamScf,1,-1,-7,-1,one_alignment_only=True)[0].score

		scores.append([str(chr_id),str(scf_id),upstream_score,numt_score,downstream_score])
	except:
		print(row)

scores=[]
alignments.apply(pairwise_alignments,axis=1)
df=pd.DataFrame(
	data=scores,
	columns=['chr_id','scf_id','upstream_score','numt_score','downstream_score']
	)

df.to_csv('../data/revision/alignment_scores.csv')