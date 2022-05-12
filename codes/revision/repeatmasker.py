#importing the modules
import numpy as np
import pandas as pd
from Bio import pairwise2
from subprocess import run,call

#read in the filtered alignments
alignments=pd.read_csv('../data/revision/filtered_alignments.csv',index_col=0)

#function for aligning the sequences
def seq_preparation(row,outfile):
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

	outfile.write('>'+f'{str(int(chr_id))}:{str(int(chr_start))}-{str(int(chr_end))}_ChrNumt'+'\n'+numtChr+'\n')
	outfile.write('>'+f'{scf_id}:{str(int(scf_start))}-{str(int(scf_end))}_ScfNumt'+'\n'+numtScf+'\n')

	outfile.write('>'+f'{str(int(chr_id))}:{str(int(chr_start-chr_part_size))}-{str(int(chr_start))}_ChrUp'+'\n'+upstreamChr+'\n')
	outfile.write('>'f'{scf_id}:{str(int(scf_start-scf_part_size))}-{str(int(scf_start))}_ScfUp'+'\n'+upstreamScf+'\n')

	outfile.write('>'+f'{str(int(chr_id))}:{str(int(chr_end))}-{str(int(chr_end+chr_part_size))}_ChrDown'+'\n'+downstreamChr+'\n')
	outfile.write('>'+f'{scf_id}:{str(int(scf_end))}-{str(int(scf_end+scf_part_size))}_ScfDown'+'\n'+downstreamScf+'\n')

with open('../data/revision/rm_input.fa','w')as outfile:
	alignments.apply(seq_preparation,args=(outfile,),axis=1)

call('RepeatMasker -species "oryctolagus cuniculus" ../data/revision/rm_input.fa',shell=True)