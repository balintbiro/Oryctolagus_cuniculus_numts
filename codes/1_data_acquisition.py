#import the required modules
import os
import gzip
import shutil
import pandas as pd
import urllib.request

#define the genome filename
genome_filename='Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz'
mt_filename='Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.MT.fa.gz'

#create a vector for the genomes
genomes=pd.Series([genome_filename, mt_filename])

#function for downloading the data and decompress the files
def get_data(filename):
    #download the gunzipped files
    url=f'http://ftp.ensembl.org/pub/release-104/fasta/oryctolagus_cuniculus/dna/{filename}'
    output_dir=os.path.join('../data/')
    filepath=output_dir+filename
    urllib.request.urlretrieve(url, filepath)

#get the data for each genome
genomes.apply(get_data)

#function for decompressing gunzipped files
def decompress(filename):
    decompressed_filename=filename[:-3]
    with gzip.open(f'../data/{filename}', 'rb')as infile, open(f'../data/{decompressed_filename}', 'wb')as outfile:
        shutil.copyfileobj(infile, outfile)
    os.remove(f'../data/{filename}')
    
#decompress the files    
genomes.apply(decompress)