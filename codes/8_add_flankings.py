#import the required modules
import os
import pandas as pd

#define class for getting the flanking regions
class Flanking:
    def __init__(self, g_start, g_length, mt_start, mt_length):
        self.g_start = int(g_start)
        self.g_length = int(g_length)
        self.mt_start = int(mt_start)
        self.mt_length = int(mt_length)
        self.directory = os.path.join('../data/fasta_files/')
        
    def g_sequence(self, g_id):
        """Get the genomic sequence."""
        g_filename = '%s.fa' % g_id
        g_filepath = self.directory + g_filename
        with open (g_filepath) as g_file:
            g_seq = g_file.readlines()[1]
        return g_seq
        
    def g_up_flanking(self, g_seq):
        """Get the upstream flanking of the genomic sequence."""
        global g_up_flankings
        flank_start = self.g_start - 200
        flank_end = self.g_start
        if flank_start < 0:
            up_flanking = g_seq[:flank_end]
        else:
            up_flanking = g_seq[flank_start:flank_end]
        g_up_flankings.append(up_flanking)
        
    def g_down_flanking(self, g_seq, g_size):
        """Get the downstream flanking of the genomic sequence."""
        global g_down_flankings
        flank_start = self.g_start + self.g_length
        flank_end = flank_start + 200
        if flank_end > g_size:
            down_flanking = g_seq[flank_start:]
        else:
            down_flanking = g_seq[flank_start:flank_end]
        g_down_flankings.append(down_flanking)
        
    def mt_up_flanking(self, mt_seq, strand):
        """Get the upstream flanking of the mitochondrial sequence."""
        global mt_up_flankings
        if strand == '-':
            flank_end = (len(mt_seq) - self.mt_start) - self.mt_length
            flank_start = flank_end - 200
            if flank_start < 0:
                up_flanking = mt_seq[:flank_end] + mt_seq[flank_start:]
            else:
                up_flanking = mt_seq[flank_start:flank_end]
            mt_up_flankings.append(up_flanking)
        else:
            flank_start = self.mt_start - 200
            flank_end = self.mt_start
            if flank_start < 0:
                up_flanking = mt_seq[:flank_end] + mt_seq[flank_start:]
            else:
                up_flanking = mt_seq[flank_start:flank_end]
            mt_up_flankings.append(up_flanking)
            
    def mt_down_flanking(self, mt_seq, strand):
        """Get the downstream flanking of the mitochondrial sequence."""
        global mt_down_flankings
        if strand == '-':
            flank_start = (len(mt_seq) - self.mt_start)
            flank_end = flank_start + 200
            if flank_end > len(mt_seq):
                circular_correction = (flank_end - flank_start) + len(mt_seq[flank_start:])
                down_flanking = mt_seq[flank_start:] + mt_seq[:circular_correction]
            else:
                down_flanking = mt_seq[flank_start:flank_end]
            mt_down_flankings.append(down_flanking)
        else:
            flank_start = self.mt_start + mt_length
            flank_end = flank_start + 200
            if flank_end > len(mt_seq):
                circular_correction = (flank_end - flank_start) + len(mt_seq[flank_start:])
                down_flanking = mt_seq[flank_start:] + mt_seq[:circular_correction]
            else:
                down_flanking = mt_seq[flank_start:flank_end]
            mt_down_flankings.append(down_flanking)
            
#reading the previously made table into an array
alignment_df = pd.read_csv('../results/1numt_array.csv')

#get the double mt sequence
mt_seq = ''
with open('../data/d_mt.fa')as mt_file:
    mt_seq = mt_file.readlines()[1]
    
#define global variables
g_up_flankings = []
g_down_flankings = []
mt_up_flankings = []
mt_down_flankings = []

#get all the flankings
for index, row in alignment_df.iterrows():
    g_id = row['g_id']
    g_start = row['g_start']
    g_length = row['g_length']
    g_size = row['g_size']
    mt_start = row['mt_start']
    mt_length = row['mt_length']
    mt_strand = row['mt_strand']
    flanking = Flanking(g_start, g_length, mt_start, mt_length)
    g_sequence = flanking.g_sequence(g_id)
    g_up_flanking = flanking.g_up_flanking(g_sequence)
    g_down_flanking = flanking.g_down_flanking(g_sequence, g_size)
    mt_up_flanking = flanking.mt_up_flanking(mt_seq, mt_strand)
    mt_down_flanking = flanking.mt_down_flanking(mt_seq, mt_strand)
    
#append the dataframe with the flanking lists
alignment_df['g_up_flanking'] = g_up_flankings
alignment_df['g_down_flanking'] = g_down_flankings
alignment_df['mt_up_flanking (if strand is negative, the coordinates are corrected)'] = mt_up_flankings
alignment_df['mt_down_flanking (if strand is negative, the coordinates are corrected)'] = mt_down_flankings

alignment_df.to_csv('../results/2numt_array.csv', index = False)