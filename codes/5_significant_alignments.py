#import the required modules
import pandas as pd
 
#e value threshold is going to be the lowest value of the r mt alignment
e_values=[]
with open('../data/r_mt_alignment.fa')as infile:
    content=pd.Series(infile.readlines())
    mask=content.apply(lambda line: 'EG2' in line)
    content[mask].apply(lambda line: e_values.append(float(line.rsplit()[3].split('=')[1])))
e_values.sort()
e_threshold=e_values[0]

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