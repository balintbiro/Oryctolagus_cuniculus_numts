#import the required modules
import os
import pandas as pd

#write the genome in a form that a given sequence is in one line
with open('../data/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa')as infile, open('../results/processed_g_dna.fa','w')as outfile:
    for index, line in enumerate(infile):
        if line[0] == '>':
            outfile.write('\n')
            outfile.write(line[:-1] + '\t')
        else:
            outfile.write(line[:-1])
            
#function for getting the required genomic ids (chromosomes and contigs)
def get_g_ids(line):
    global g_ids
    if len(line.rsplit())>2 and line[0]=='s':
        g_id=line.rsplit()[1]
        if g_id not in g_ids:
            g_ids.append(g_id)
                
#get the file content and use the previously defined (get_g_ids) function to get the required g ids
g_ids=[]#the global variable that is necessary for our function
with open('../results/significant_alignments.fa')as infile:
    content=pd.Series(infile.readlines())
    content.apply(get_g_ids)

#create directory for fasta files
if os.path.exists('../data/fasta_files/'):
    pass
else:
    os.mkdir('../data/fasta_files/')

#write individual fasta files for the required sequences
with open('../results/processed_g_dna.fa')as infile:
    for line in infile:
        if line != '\n':
            header = line.split('REF')[0]
            g_id=header[1:].rsplit()[0]
            sequence = line.split('REF')[1].rsplit()[0]
            if g_id in g_ids:
                filename=f'{g_id}.fa'
                with open(f'../data/fasta_files/{filename}','w')as outfile:
                    outfile.write(header+'\n')
                    outfile.write(sequence)
            