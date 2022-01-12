#importing the required modules
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import anderson, wilcoxon, ttest_ind

#load the csv file into a pandas df
numts=pd.read_csv('../results/3numt_array.csv')

#get the g ids
g_ids = numts['g_id']

flank_size = 5000

#get the actual length of the chromosomes
chromosome_lengths = numts['g_size']
chromosome_lengths.index = g_ids

#the start positions of the upstream flankings
upstream_starts = numts['g_start'].apply(lambda start : start - flank_size)
upstream_starts[upstream_starts < 0] = 1
upstream_starts.index = g_ids

#the end positions of the upstream flankings (which are the genomic_starts of the numts)
upstream_ends = numts['g_start']
upstream_ends.index = g_ids

#start position of downstream_flanking (which is equal to the ends of the numts)
numts['g_length'].index = g_ids
downstream_starts = numts['g_start'] + numts['g_length']
downstream_starts.index = g_ids

#end position of downstream_flanking
numt_ends = numts['g_start'] + numts['g_length']
downstream_ends = numt_ends.apply(lambda end : end + flank_size)

#get the regions
indices = np.arange(0,len(upstream_starts),1)
upstream_flankings = pd.Series(index = indices, dtype=str)
downstream_flankings = pd.Series(index = indices, dtype=str)
for index, upstream_start in enumerate(upstream_starts):
    upstream_rm_item = (str(g_ids[index]) + ':' +str(upstream_start) + '-' + str(upstream_ends[index]))
    upstream_flankings[index] = upstream_rm_item
    downstream_rm_item = ''
    if chromosome_lengths[index] < downstream_ends[index]:
        downstream_rm_item = (str(g_ids[index]) + ':' +str(downstream_starts[index]) +
                          '-' + str(chromosome_lengths[index]))
    else:
        downstream_rm_item = (str(g_ids[index]) + ':' +str(downstream_starts[index]) +
                              '-' + str(downstream_ends[index]))
    downstream_flankings[index] = downstream_rm_item

#writing flanking output
with open('../results/repeatmasker_input_upstream.bed','w') as upstream_outfile, open('../results/repeatmasker_input_downstream.bed','w') as downstream_outfile:
    upstream_flankings.apply(lambda line : upstream_outfile.write(line + '\n'))
    downstream_flankings.apply(lambda line : downstream_outfile.write(line + '\n'))
    
#get the ids for the numtless sequences dictionary
numtless_ids = pd.Series(numts['g_id'].unique(), dtype=str)

#get numtless sequences
def get_numtless_sequences(g_id):
    filename = f'{g_id}.fa'
    with open(f'../data/numtless_sequences/{filename}') as infile:
        content = infile.readlines()
        try:
            sequence = content[0]
        except IndexError:
            sequence = np.nan
    return sequence
    
#create the dictionary (vector, series) for the numtless sequences
numtless_sequences = numtless_ids.apply(get_numtless_sequences)
numtless_sequences.index = numtless_ids
numtless_sequence=numtless_sequences.dropna()

#get the number of numts per each genomic region (chromosomes or scaffolds)
identifiers = pd.Series(numtless_sequences.index.values)
numt_count = identifiers.apply(lambda g_id : list(numts['g_id']).count(g_id))
numt_count.index = numtless_sequences.index.values

#get the lengths of the numtless sequences
lengths = numtless_sequences.apply(lambda sequence : len(sequence))

#function for sampling the genome corresponding to the number of numts per a given genomic region (chromosome or scaffold)
def sampling_genome(g_id):
    global sum_ranges
    number = numt_count[g_id]
    sample_ranges = []
    for i in range(0,number,1):
        np.random.seed(i + 5)
        length = lengths[g_id]
        sequence = numtless_sequences[g_id]
        range_max = length - flank_size
        if range_max < 0:
            range_max = 2
        range_start = np.random.randint(0, range_max)
        range_end = range_start + flank_size
        if range_end >= length:
            range_end == length-1
        sample_ranges.append(g_id + ':' + str(range_start) + '-' + str(range_end))
    sum_ranges += sample_ranges
    return sample_ranges
    
#get the sample ranges
sum_ranges = []
identifiers.index = identifiers
samples = identifiers.apply(sampling_genome)

#writing output file
with open('../results/repeatmasker_input_genomic_samples.bed', 'w') as outfile:
    pd.Series(sum_ranges).apply(lambda sample_range : outfile.write(sample_range + '\n'))


#################################################
#These three files                              #
#repeatmasker_input_downstream.bed,             #
#repeatmasker_input_upstream.bed                #
#and repeatmasker_input_genomic_samples.bed     #
#should be uploaded to the RepeatMasker server! #
#And the output from the server should be placed#
#into the results/ folder!                      #
#################################################

#reading the RepeatMasker output files into pandas dfs and set the genome parts (chromosomes and scaffolds) as indices
upstream_repeats=pd.read_csv('../results/upstream_repeats.tsv', sep='\t')
upstream_repeats=upstream_repeats.set_index('genoName')
downstream_repeats=pd.read_csv('../results/downstream_repeats.tsv', sep='\t')
downstream_repeats=downstream_repeats.set_index('genoName')
sample_repeats=pd.read_csv('../results/genomic_repeats.tsv', sep='\t')
sample_repeats=sample_repeats.set_index('genoName')

#get common repetitions in all three dataset (upstream, downstream, genome)
common_repeats = pd.Series(list(set(upstream_repeats['repName'].unique()) &
                             set(downstream_repeats['repName'].unique()) &
                             set(sample_repeats['repName'].unique())))
                             
#get the number of common repetitions per chromosomes and scaffolds
def get_repeatnumber(chromosome, repname):
    upstream_subdf = upstream_repeats.loc[chromosome]
    downstream_subdf = downstream_repeats.loc[chromosome]
    sample_subdf = sample_repeats.loc[chromosome]
    return [list(upstream_subdf['repName']).count(repname),
           list(downstream_subdf['repName']).count(repname),
           list(sample_subdf['repName']).count(repname)]

#empty series for the repetitions
repeatnumbers = pd.Series(index = common_repeats, dtype = str)

#get the repeatnumber of every common repeats
for common_repeat in common_repeats:
    chromosomes = pd.Series(list(set(np.unique(sample_repeats.index.values))&
                           set(np.unique(upstream_repeats.index.values))&
                           set(np.unique(downstream_repeats.index.values))))
    samples = chromosomes.apply(get_repeatnumber, args = (common_repeat,))
    upstream_reps = []
    samples.apply(lambda sample : upstream_reps.append(sample[0]))
    downstream_reps = []
    samples.apply(lambda sample : downstream_reps.append(sample[1]))
    sample_reps = []
    samples.apply(lambda sample : sample_reps.append(sample[2]))
    repeatnumbers[common_repeat] = [upstream_reps, downstream_reps, sample_reps]
    
#function for the statistical analysis of repetitions
#upstream_repetitions = nested_list[0]
#downstream_repetitions = nested_list[1]
#sample_repetitions = nested_list[2]
def statistical_analysis(nested_list, which):
    sample1 = nested_list[which[0]]
    sample2 = nested_list[which[1]]
    norm1 = anderson(sample1)
    norm2 = anderson(sample2)
    stat1 = norm1[0]
    stat2 = norm2[0]
    critical_value1 = norm1[1][2]
    critical_value2 = norm2[1][2]
    if (stat1 > critical_value1) or (stat2 > critical_value2):
        try:
            significance = wilcoxon(sample1, sample2)[1]
            return significance
        except ValueError:
            pass
    else:
        try:
            significance = ttest_ind(sample1, sample2)[1]
            return significance
        except ValueError:
            pass
            
#statistics of upstream and downstream repetitions
upstream_downstream_stat = repeatnumbers.apply(statistical_analysis, args = ([0,1],))

#statistics of upstream and sample repetitions
upstream_sample_stat = repeatnumbers.apply(statistical_analysis, args = ([0,2],))

#statistics of downstream and sample repetitions
downstream_sample_stat = repeatnumbers.apply(statistical_analysis, args = ([1,2],))

#get the repeat category for every repeats
repeat_categories = pd.Series(repeatnumbers.index.values)
repeat_classes = repeat_categories.apply(lambda category : upstream_repeats.loc[upstream_repeats['repName'] == category]['repClass'].unique()[0])
repeat_classes.index = repeatnumbers.index.values

#invert repeat categories
repeats = pd.Series(repeat_classes.index.values)
repeats.index = repeat_classes

#get individual repeats for all repeat classes
repeat_families = pd.Series(repeat_classes.unique()).apply(lambda repeat_class:repeats[repeat_class].tolist())
repeat_families.index=repeat_classes.unique()

#get the new ordered header
def get_header(repeat_list):
    global header
    header += repeat_list
    
header = []
repeat_families.apply(get_header)

#create dataframe from the significance values
df = pd.DataFrame([upstream_downstream_stat,
upstream_sample_stat,
downstream_sample_stat])
df = df[header]
df.index = (['upstream_downstream',
'upstream_sample',
'downstream_sample'])

#writing the output for visualisation
df.to_csv('../results/repeatmasker_results.csv')
with open('../results/repeat_families.txt','w')as outfile:
    outfile.write(','.join(list(repeat_families.index.values)))
    outfile.write('\n')
    repeat_families.apply(lambda repeat_list: outfile.write(','.join(repeat_list)+'\n'))