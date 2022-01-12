#importing the required modules
import os
import pandas as pd

#reading the csv into a df
numts=pd.read_csv('../results/3numt_array.csv')

#get the g_ids
g_ids=pd.Series(numts['g_id'].unique())

#function for getting the chr or scf regions without numts
def numt_coordinates(g_id):
    subdf=numts.loc[numts['g_id']==g_id]
    subdf=subdf.sort_values(by='g_start')
    g_length=subdf['g_size'].tolist()[0]
    numt_starts=[]
    subdf.apply(lambda row: numt_starts.append(row['g_start']), axis=1)
    numt_ends=[]
    subdf.apply(lambda row: numt_ends.append(row['g_start']+row['g_length']), axis=1)
    regions_to_include=[]
    if len(numt_starts)==1:
        regions_to_include.append([0,numt_starts[0]])
        regions_to_include.append([numt_ends[0],g_length])
    else:
        for index, start in enumerate(numt_starts):
            end=numt_ends[index]
            if index==0:
                regions_to_include.append([0,start])
                regions_to_include.append([end,numt_starts[index+1]])
            elif index==len(numt_starts)-1:
                regions_to_include.append([end,g_length])
            else:
                regions_to_include.append([end,numt_starts[index+1]])
    return regions_to_include
    
#get the regions without numts
numtless_regions=g_ids.apply(numt_coordinates)
numtless_regions.index=g_ids

#function for selecting the numtless sequence regions for each chr and/or scf and merge them together
def get_numtless_sequences(regions):
    global g_id_index
    try:
        filename=f'{numtless_regions.index.values[g_id_index]}.fa'
        sequence=''
        with open(f'../data/fasta_files/{filename}')as infile:
            content=infile.readlines()
            sequence=content[1]
        regions=pd.Series(regions)
        sequence_parts=[]
        regions.apply(lambda region: sequence_parts.append(sequence[region[0]:region[1]]))
        numtless_sequence=''.join(sequence_parts)
        g_id_index+=1
        return numtless_sequence
    except FileNotFoundError:
        print(numtless_regions.index.values[g_id_index])
        g_id_index+=1
    
#get the numtless sequences
g_id_index=0
numtless_sequences=numtless_regions.apply(get_numtless_sequences)
numtless_sequences.index=g_ids

#create directory for the numtless sequences
if os.path.exists('../data/numtless_sequences/'):
    pass
else:
    os.mkdir('../data/numtless_sequences/')
    
#define a function  for writing numtless output files
def write_outputs(sequence):
    global g_id_index
    if type(sequence)==str:
        filename=f'{numtless_regions.index.values[g_id_index]}.fa'
        with open(f'../data/numtless_sequences/{filename}','w')as outfile:
            outfile.write(sequence)
    g_id_index+=1
    
#writing the actual outputs (numtless sequences)
g_id_index=0
numtless_sequences.apply(write_outputs)