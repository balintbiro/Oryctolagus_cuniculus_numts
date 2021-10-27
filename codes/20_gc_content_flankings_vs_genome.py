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

#get the flanking sequences
flanking_sequences=numts['g_up_flanking']+numts['g_down_flanking']
flanking_sequences.index = numts['g_id']
flanking_sequences=flanking_sequences.dropna()

#get the gc content of flanking regions
flankings_gc = flanking_sequences.apply(lambda sequence : (sequence.upper().count('G') + sequence.upper().count('C')) / len(sequence))

#define function for sampling genome based on flanking size
def flanking_based_sampling(g_id):
    global seed_value
    sequence_to_sample = numtless_sequences[g_id]
    seed_value += 1
    np.random.seed(seed_value)
    flanking_size = 400
    start = np.random.randint(0, len(sequence_to_sample) - flanking_size)
    end = start + flanking_size
    sample_sequence = sequence_to_sample[start:end]
    unknown_nucleotides = sample_sequence.count('N') / len(sample_sequence)
    if unknown_nucleotides < 0.05:
        return sample_sequence
    else:
        seed_value += 1
        np.random.seed(seed_value)
        start = np.random.randint(0, len(sequence_to_sample) - flanking_size)
        end = start + flanking_size
        sample_sequence = sequence_to_sample[start:end]
        return sample_sequence
        
#flanking based genomic samples (flanking/genome-->flage)
seed_value=0
flage_samples = pd.Series(flankings_gc.index).apply(flanking_based_sampling)
flage_samples.index = flankings_gc.index.values

#claculate the gc content of the flanking based genomic samples
flages_gc = flage_samples.apply(lambda sequence : (sequence.upper().count('G') + sequence.upper().count('C')) / len(sequence))

#normality testing of flanking gc
flanking_normality = anderson(flankings_gc)
flanking_stat = flanking_normality[0]
flanking_critical_value = flanking_normality[1][2]

#normality testing of flanking based genom samples (flage) gc
flage_normality = anderson(flages_gc)
flage_stat = flage_normality[0]
flage_critical_value = flage_normality[1][2]

#statistics (normality testing and significance)
#flanking vs genome
if (flage_stat > flage_critical_value) or (flanking_stat > flanking_critical_value):
    flage_significance = wilcoxon(flages_gc, flankings_gc)
else:
    flage_significance = ttest_ind(flages_gc, flankings_gc)
flage_significance[1]

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
gcs = [flankings_gc, flages_gc]
sns.set_style("darkgrid")
fig, ax=plt.subplots(figsize=(3,2))
sns.violinplot(data = gcs)
ax.set_ylabel('GC content',
             fontsize = 14)
ax.set_xticklabels(['flankings', 'genome'],
                  fontsize = 14)
statistical_annotation([flankings_gc, flages_gc], flage_significance[1], [0,1], 0.1)
ax.set_ylim(0, 1.1)
plt.tight_layout()
plt.savefig('../results/gc_contents_flanking_vs_genome.png', dpi = 400)