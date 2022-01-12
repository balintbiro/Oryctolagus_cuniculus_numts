import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import anderson, wilcoxon, ttest_ind

mt_sequence=''
with open('../data/Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.MT.fa')as infile:
    for line in infile:
        if '>' not in line:
            mt_sequence+=line[:-1]
        
mt_starts=pd.Series(np.round((np.linspace(0,len(mt_sequence),80)),))[:-1]
mt_ends=list(mt_starts)
mt_ends.append(len(mt_sequence))

def mt_gc_samples(start, ends, sequence):
    start=int(start)
    index=mt_starts[mt_starts==start].index[0]
    end=int(ends[index+1])
    sample_sequence=sequence[start:end].upper().replace('N','')
    gc_content=(sample_sequence.count('G')+sample_sequence.count('C'))/len(sample_sequence)
    return gc_content
    
mt_samples=mt_starts.apply(mt_gc_samples, args=(mt_ends,mt_sequence,)).tolist()

def sampling_genome(filename):
    global genome_gcs
    filepath=os.path.join('../data/fasta_files/'+filename)
    sequence=''
    with open(filepath)as infile:
        content=infile.readlines()
        sequence=content[1]
    starts=pd.Series(np.arange(0,len(sequence),3))
    gcs=[]
    sample_sequences=starts.apply(lambda start: sequence[start:(start+218)].upper().replace('N',''))
    try:
         sample_sequences.apply(lambda sample_sequence: gcs.append((sample_sequence.count('G')+sample_sequence.count('C'))/len(sample_sequence)))
         genome_gcs+=gcs
    except ZeroDivisionError:
         pass
    print(files[files==filename].index[0],'/',len(files))

genome_gcs=[]
files=pd.Series(os.listdir('../data/fasta_files'))
files.apply(sampling_genome)

np.random.seed(0)
nuclear_samples=np.random.choice(genome_gcs,len(mt_samples))

#normality testing of mt gc
mt_normality = anderson(mt_samples)
mt_stat = mt_normality[0]
mt_critical_value = mt_normality[1][2]

#normality testing of genome gc
g_normality = anderson(nuclear_samples)
g_stat = g_normality[0]
g_critical_value = g_normality[1][2]

if (mt_stat>mt_critical_value) or (g_stat>g_critical_value):
    significance=wilcoxon(mt_samples, nuclear_samples)
else:
    significance=ttest_ind(mt_samples, nuclear_samples)

#writing the output for visualisation
with open('../results/gcs_for_visualisation.txt','a')as outfile:
    outfile.write(f'genomes_gc\nsignificance: {str(significance[1])}\n')
    mt_samples_gc=pd.Series(mt_samples)
    mt_samples_gc.apply(lambda value:outfile.write(str(value)+','))
    outfile.write('\n')
    nuclear_samples_gc=pd.Series(nuclear_samples)
    nuclear_samples_gc.apply(lambda value:outfile.write(str(value)+','))
    outfile.write('\n')