#import required modules
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import anderson, wilcoxon, ttest_ind

#read numts csv file into a df
numts = pd.read_csv('../results/3numt_array.csv')

#get the gene ids
g_ids=pd.Series(numts['g_id'].unique())

#create vector for the numtless sequences
def create_dictionary(g_id):
    global numtless_sequences
    filename = f'{g_id}.fa'
    sequence = ''
    with open (f'../data/numtless_sequences/{filename}') as infile:
        content = infile.readlines()
        try:
            sequence = content[0]
        except IndexError:
            pass
    numtless_sequences[g_id] = sequence
    
#create vector for numtless sequences
numtless_sequences = pd.Series(index = g_ids, dtype=str)
g_ids.apply(create_dictionary)

#get the length of every numts of each chromosomes/scaffolds
lengths = g_ids.apply(lambda g_id : numts.loc[numts['g_id'] == g_id]['g_length'].tolist())
lengths.index = g_ids

#function for sampling the genome reagrding the corresponding numt size
def nuge_sampling(g_id):
    sequence_to_sample = numtless_sequences[g_id]
    numt_sizes = lengths[g_id]
    samples = []
    for index, numt_size in enumerate(numt_sizes):
        seed_value=numt_size
        np.random.seed(seed_value)
        start = np.random.randint(0, len(sequence_to_sample) - numt_size)
        end = start + numt_size
        sample_sequence = sequence_to_sample[start:end]
        unknown_nucleotides = sample_sequence.count('N') / len(sample_sequence)
        if unknown_nucleotides < 0.05:
            samples.append(sample_sequence)
        else:
            seed_value += 1
            np.random.seed(seed_value)
            start = np.random.randint(0, len(sequence_to_sample) - numt_size)
            end = start + numt_size
            sample_sequence = sequence_to_sample[start:end]
            samples.append(sample_sequence)
    return samples
    
#get genom samples based on numts
nuge_samples = g_ids.apply(nuge_sampling)
nuge_samples.index = g_ids

#function for calculating gc contents of sequence samples
def gc_content(sample_sequences):
    gc_contents = []
    sample_sequences=pd.Series(sample_sequences)
    sample_sequences.apply(lambda sequence: gc_contents.append((sequence.upper().count('G')+sequence.upper().count('C'))/len(sequence)))
    return gc_contents
    
#calculate gc content of sample sequences
gc_content_samples = nuge_samples.apply(gc_content)
gc_content_samples.index = g_ids

#calculate the gc content of numts
def numts_gc(g_id):
    df = numts.loc[numts['g_id'] == g_id]
    gc_content = df['g_sequence'].apply(lambda sequence : (sequence.upper().count('G') + sequence.upper().count('C')) / len(sequence))
    return gc_content.tolist()
    
#calculate the gc content of numts
numts_gc = g_ids.apply(numts_gc)
numts_gc.index = g_ids

#function for merging gc contents
def merge(gc_content):
    global sum_gc
    sum_gc += gc_content
    
#get merged gc content of numts
sum_gc = []
numts_gc.apply(merge)
numts_gc = sum_gc
sum_gc = []
gc_content_samples.apply(merge)
nuge_gc = sum_gc

#normality testing of numts gc
numts_normality = anderson(numts_gc)
numts_stat = numts_normality[0]
numts_critical_value = numts_normality[1][2]

#normality testing of genome gc
nuge_normality = anderson(nuge_gc)
nuge_stat = nuge_normality[0]
nuge_critical_value = nuge_normality[1][2]

#statistics (normality testing and significance)
#numt vs genome
if (numts_stat > numts_critical_value) or (nugee_stat > nuge_critical_value):
    nuge_significance = wilcoxon(numts_gc, nuge_gc)
else:
    nuge_significance = ttest_ind(numts_gc, nuge_gc)
nuge_significance[1]

#create a function for the statistical annotation of the graph
def statistical_annotation(data, significance, positions, height):
    x1, x2 = positions[0],positions[1]
    maximum = max([max(data[0]),max(data[1])])
    y, h, col = maximum + height + 0.03, 0.03, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c = col)
    if significance < 0.05:
        plt.text((x1+x2)*.5, y+h, '$\it{P < 0.05}$', ha='center', va='bottom', color = col, fontsize = 14)
    else:
        plt.text((x1+x2)*.5, y+h, "n.s.", ha='center', va='bottom', color = col, fontsize = 14)
        
#creating the figure
gcs = [numts_gc, nuge_gc]
sns.set_style("darkgrid")
fig, ax=plt.subplots(figsize=(3,2))
sns.violinplot(data = gcs)
ax.set_ylabel('GC content',
             fontsize = 14)
ax.set_xticklabels(['numts', 'genome'],
                  fontsize = 14)
statistical_annotation([numts_gc, nuge_gc], nuge_significance[1], [0,1], 0.1)
ax.set_ylim(0.1, 1.05)
plt.tight_layout()
plt.savefig('../results/gc_contents_numt_vs_genome.png', dpi = 600)