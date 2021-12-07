#import the required modules
import os
from subprocess import call

#change working directory to the actual 'data' directory
os.chdir('../data/')

#do the alignment with Python embedded Linux commands
call('lastdb ocdb Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa', shell=True)  #create "database files" for the alignment
call('lastal -r1 -q1 -a7 -b1 ocdb d_mt.fa > duplicated_mt_alignment.fa', shell=True)#align the duplicated mt