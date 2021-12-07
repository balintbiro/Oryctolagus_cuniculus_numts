#importing the required modules
import numpy as np
import pandas as pd

#read the numts csv file into a dataframe
numts=pd.read_csv('../results/3numt_array.csv', index_col=0)

#get the scaffolds
scf_mask=numts['g_id'].apply(lambda g_id: len(g_id) > 3 and g_id!='MT')
scf_numts=numts[scf_mask]

#read the mt annotation (MITOS) into a df
mt_annotation=pd.read_csv('../data/mitos_annotation.bed',
                         sep='\t',
                         header=None)
mt_annotation.columns=['Chr','Start','End','Name','Sig','Strand']

#calculate the length (in bp) of each mt part
informations=[]
non_annotated_length=(17245-sum(mt_annotation.apply(lambda row: row['End']-row['Start'], axis=1)))
mt_annotation.apply(lambda row: informations.append(['',1,row['End']-row['Start'],969696,'mtDNA',12,252525]), axis=1)

#add the mt part without annotation to the informations variable
informations.append(['',1,non_annotated_length,969696,'mtDNA',12,252525])

#define function for getting the scaffold specific sizes
def scf_prep(g_id):
    global informations
    subdf=scf_numts.loc[scf_numts['g_id']==g_id]
    
    chr_start=1
    chr_end=list(subdf['g_size'])[0]
    fill=969696
    species='Scaffolds'
    size=12
    color=252525
    informations.append(['',chr_start,chr_end,fill,species,size,color])

#get the scaffold specific informations
scf_numts=scf_numts.sort_values(by='mt_start')
pd.Series(scf_numts['g_id'].unique()).apply(scf_prep)
informations=pd.DataFrame(informations)
informations.columns=['Chr','Start','End','fill','species','size','color']

#writing the output
informations.to_csv('../results/scaffold_synteny_dual_karyotype.csv', index=False)

#define function for creating sets for every mt part (based on the MITOS annotation)
def prepare_mt_genes(row):
    global mt_genes
    mt_genes[row['Name']]=set(np.arange(int(row['Start']),int(row['End'])))
    
#global variable for mt genes
mt_genes={}

#get the mt genes
mt_annotation.apply(prepare_mt_genes, axis=1)
mt_genes=pd.Series(mt_genes)

#colorcodes for synteny
#rapidtables.com/web/color/RGB_Color.html
#https://r-charts.com/colors/
colors=['FF0000',#red
       'FFA500',#orange
       'DAA520',#golden rod
       'EEE8AA',#pale golden rod
       'BDB76B',#dark khaki
       '6B8E23',#olive drab
       '228B22',#lime
       '8FBC8F',#dark sea green
       '00FFFF',#aqua
       'AFEEEE',#pale turquoise
       '7FFFD4',#aqua marine
       'B0E0E6',#powder blue
       '5F9EA0',#cadet blue
       '6495ED',#corn flower blue
       '00BFFF',#deep sky blue
       '0000FF',#blue
       '8A2BE2',#blue violet
       'BA55D3',#medium orchid
       '800080',#purple
       'D8BFD8',#thistle
       'DDA0DD',#plum
       'FF00FF',#magenta
       'FFC0CB',#pink
       'F5F5DC',#beige
       'F5DEB3',#wheat
       'D2691E',#chocolate
       'F4A460',#sandy brown
       'BC8F8F',#rosy brown
       'B0C4DE',#light steel blue
       '000000',#black
       'A9A9A9',#dark grey
       'D3D3D3',#light grey
       'FFDEAD',#navajo white
       'FFE4E1',#misty rose
       'FF69B4',#hot pink
       'E0FFFF',#light cyan
       'FFFF00',#yellow
       'FFFACD',#lemon chiffon
       'CD853F',#peru
       'E6E6FA',#lavender
       'C0C0C0',#silver
       '1E90FF',#dodger blue
       '98FB98',#pale green
       '808000',#olive
       'CD5C5C',#indian red
       'BDB76B',#dark khaki
       '2E8B57',#sea green
       '20B2AA',#light sea green
       'FF6347',#tomato
       '191970',#midnight blue
       '4169E1'#royal blue
       ]

#define colors in a randomized way
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
            if mt_end_index>mt_start_index:
                common_collections.append([index+1,mt_start_index,mt_end_index,
                                   list(scf_numts['g_id'].unique()).index(g_id)+1,g_start,
                                   (g_start+len(common_elements)),colors[color_index%len(colors)]])
    dual_synteny+=common_collections
    color_index+=1
    
dual_synteny=[]
color_index=0
scf_numts.apply(common_nucleotides, axis=1)
dual_synteny=pd.DataFrame(dual_synteny)
dual_synteny.columns=['Species_1','Start_1','End_1','Species_2','Start_2','End_2','fill']

dual_synteny.to_csv('../results/scaffold_synteny_dual_synteny.csv',index=False)