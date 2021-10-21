#import the required modules
import pandas as pd

#get the e values from the reversed mt alignment
e_values=[]
lines_not_in_use=[]
with open('../data/reversed_mt_alignment.fa')as infile:
    content=pd.Series(infile.readlines())
    content.apply(lambda line: e_values.append(float(line.rsplit()[3].split('=')[1])) if 'EG2' in line else lines_not_in_use.append(line))
    
#e value threshold is going to be the lowest value of the r mt alignment
e_threshold=min(e_values)

#process the duplicated alignment and filter the signifcant alignments, then write them into a separated file
with open('../data/duplicated_mt_alignment.fa')as infile, open('../results/significant_alignments.fa','w')as outfile:
    content=infile.readlines()
    for index, line in enumerate(content):
        if 'score' in line:
            e_value = float(line.rsplit()[3].split('=')[1])
            g_sequence = content[index + 1]
            mt_sequence = content[index + 2]
            if e_value < e_threshold:
                outfile.write(line)
                outfile.write(g_sequence)
                outfile.write(mt_sequence)
                outfile.write('\n')