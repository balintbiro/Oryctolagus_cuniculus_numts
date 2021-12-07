#import the required modules
import random
import numpy as np
import pandas as pd

#load the numt related csv into a dataframe
numts=pd.read_csv('../results/3numt_array.csv', index_col=0)

#get the cromosomes
chr_mask=numts['g_id'].apply(lambda g_id: len(g_id) < 3 and g_id!='MT' and g_id!='X')
chr_numts=numts[chr_mask]

#define function for getting the chromosome specific sizes
def chr_prep(g_id):
    global informations
    subdf=chr_numts.loc[chr_numts['g_id']==g_id]
    
    chr_start=1
    chr_end=list(subdf['g_size'])[0]
    fill=969696
    species='gDNA'
    size=12
    color=252525
    informations.append([int(g_id),chr_start,chr_end,fill,species,size,color])

#define global variable
informations=[]

#append the data for the chromosomes that do not contain any numts
informations.append([5,1,37992211,969696,'gDNA',12,252525])#chr 5
informations.append([6,1,27502587,969696,'gDNA',12,252525])#chr 6

#create vector for the individual g ids
pd.Series(np.unique(chr_numts['g_id'])).apply(chr_prep)

#create dataframe and sort the df by the chromosome id and then transform it back to nested list
informations=pd.DataFrame(informations)
informations.columns=['Chr','Start','End','fill','species','size','color']
informations=informations.sort_values(by='Chr')
informations=informations.values.tolist()

#get the mitochondrial annotation (MITOS)
mt_annotation=pd.read_csv('../data/mitos_annotation.bed',
                         sep='\t',
                         header=None)
mt_annotation.columns=['Chr','Start','End','Name','Sig','Strand']

#calculate the length (in bp) of each mt part
non_annotated_length=(17245-sum(mt_annotation.apply(lambda row: row['End']-row['Start'], axis=1)))
mt_annotation.apply(lambda row: informations.append(
    ['',1,row['End']-row['Start'],969696,'mtDNA',12,252525]), axis=1)
informations.append(['',1,non_annotated_length,969696,'mtDNA',12,252525])
informations=pd.DataFrame(informations)
informations.columns=['Chr','Start','End','fill','species','size','color']

#write the output file
informations.to_csv('../results/gene_synteny_dual_karyotype.csv', index=False)

#define function for creating sets for every mt part (based on the MITOS annotation)
def prepare_mt_genes(row):
    global mt_genes
    mt_genes[row['Name']]=set(np.arange(int(row['Start']),int(row['End'])))

#global variable for mt genes
mt_genes={}

#get the mt genes
mt_annotation.apply(prepare_mt_genes, axis=1)
mt_genes=pd.Series(mt_genes)

#get the numts that are inserted into genes based on Ensemble
gene_mask=numts['ensembl_description'].apply(lambda g_id: type(g_id)==str)
gene_numts=numts[gene_mask]

#filter out the one numt that is inserted into a scaffold localized gene
gene_mask=gene_numts['g_id'].apply(lambda g_id: g_id!='GL018744')
gene_numts=gene_numts[gene_mask]

gene_numts=gene_numts.sort_values(by=['mt_start'])

#define colors in a randomized way for synteny analyses
color_part=lambda: random.randint(0,255)
color_ids=pd.Series(np.arange(1,100,1))
colors=color_ids.apply(lambda number: '%02X%02X%02X' % (color_part(),color_part(),color_part())).tolist()
       
#define function for getting the synteny blocks
def common_nucleotides(row):
    global dual_synteny
    global color_index
    numt_range=set(np.arange(row['mt_start'],(row['mt_start']+row['mt_length'])))
    g_id=row['g_id']
    g_start=row['g_start']
    common_collections=[]
    for index, mt_range in enumerate(mt_genes):
        common_elements=numt_range&mt_range
        if len(common_elements)>0:
            mt_start_index=list(mt_range).index(list(common_elements)[0])
            mt_end_index=list(mt_range).index(list(common_elements)[-1])
            common_collections.append([int(g_id),g_start,(g_start+len(common_elements)),
                                      index+1,mt_start_index,mt_end_index,colors[color_index]])#colors[color_index]])
    color_index+=1
    dual_synteny+=common_collections
    
dual_synteny=[]
color_index=0
gene_numts.apply(common_nucleotides, axis=1)
dual_synteny=pd.DataFrame(dual_synteny)
dual_synteny.columns=['Species_1','Start_1','End_1','Species_2','Start_2','End_2','fill']

#writing the synteny blocks into a csv file
dual_synteny.to_csv('../results/gene_synteny_dual_synteny.csv',index=False)