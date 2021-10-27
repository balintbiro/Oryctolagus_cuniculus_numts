#importing the required modules
import pandas as pd
import requests, sys

#read master array into a dataframe
master_array = pd.read_csv('../results/2numt_array.csv')

#function for obtaining data from Ensembl database
def get_ensembl_data(g_id, g_start, g_length):
    end = g_start + g_length
    ext = '/overlap/region/rabbit/%s:%s-%s?feature=gene' % (g_id, g_start, end)
    server = 'https://rest.ensembl.org'
    req = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not req.ok:
        req.raise_for_status()
        sys.exit()
    if len(req.json()) == 0:
        return [float('Nan'), float('Nan')]
    elif 'gene_id' in req.json()[0] and 'description' in req.json()[0]:
        gene_id = req.json()[0]['gene_id']
        description = req.json()[0]['description']
        return [gene_id, description]
    elif 'gene_id' in req.json()[0] and 'description' not in req.json()[0]:
        gene_id = req.json()[0]['gene_id']
        description = float('Nan')
        return [gene_id, description]
    elif 'gene_id' not in req.json()[0] and 'description' in req.json()[0]:
        gene_id = float('Nan')
        description = req.json()[0]['description']
        return [gene_id, description]
        
#get the gene ids and descriptions for the specified genomic regions
gene_ids = []
descriptions = []
for index, row in master_array.iterrows():
    result = get_ensembl_data(row['g_id'], row['g_start'], row['g_length'])
    gene_id = result[0]
    description = result[1]
    gene_ids.append(gene_id)
    if description == None:
        descriptions.append(float('Nan'))
    else:
        descriptions.append(description)
        
#update the master array with the gene_id and description columns
master_array['ensembl_gene_id'] = gene_ids
master_array['ensembl_description'] = descriptions

#writing output file
master_array.to_csv('../results/3numt_array.csv', index = False)