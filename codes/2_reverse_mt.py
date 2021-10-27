#import the required modules
import pandas as pd

#define the genome filename
mt_filename='Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.MT.fa'

#reverse mt dna
header=[]
mt_sequence=[]
with open(f'../data/{mt_filename}')as infile, open(f'../data/r_mt.fa','w')as outfile:
    content=pd.Series(infile.readlines())
    content.apply(lambda line: header.append(line) if '>' in line else mt_sequence.append(line[:-1]))
    mt_sequence=''.join(mt_sequence)
    r_mt_sequence=mt_sequence[::-1]#actual reversal step
    outfile.write(header[0])
    outfile.write(r_mt_sequence)