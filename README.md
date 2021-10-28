# numt
 Patterns of numtogenesis in the rabbit (Oryctolagus cuniculus) genome.

This repository contains all the codes for identifying nuclear insertions with mitochondrial origin in the rabbit (Oryctolagus cuniculus) genome.
The initial folder structure is made up of three folders:
- codes
- data
- results

The alignments were performed with LASTAL.
The codes are mainly written in Python but some of the visualizations are done in R.

Used Python packages (with version number if it is available):
- os
- gzip
- shutil
- pandas (1.2.4)
- urllib.request
- subprocess
- requests (2.26.0)
- numpy (1.20.2)
- seaborn (0.11.2)
- matplotlib (3.4.3)
- scipy (1.6.2)

Used R packages (with version number if it is available):
- RIdeogram (0.2.2)
- Biostrings (2.58.0)
- TFBSTools (1.28.0)

To annotate the mitochondrion, MITOS server (http://mitos.bioinf.uni-leipzig.de/index.py) was used with the following setting(s):
- Genetic Code: 02 - Vertebrate

RepeatMasker was also run on server (http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1196065643_I1OZNebCc0toA4hkQLRwVjIcJUTg) with the following setting(s):
- clade: Mammal
- genome: Rabbit
- group: Variation and Repeats
- track: RepeatMasker
- table: rmsk
- define regions: upload the mt FASTA

Brief description of the programs:

1_data_acquisition

- download genomic dna and mitochondrial dna in .gz format
- decompress the .gz files

2_reverse_mt

- actual reversal of the mitochondrial dna

3_duplicate_mt_dna

- prepare the duplicated mt dna for the later alignment

4_alignments

- prepare the database with LASTAL from the rabbit genome
- align reverse mt dna with nuclear genome
- align double mt dna with nuclear genome

5_significant_alignments

- get the lowest e value from the reverse mt dna alignment
- use this threshold for masking the double mt dna alignment

6_writing_fasta_files

- process nuclear dna file into a format where the sequence is in one line
- get the ids of the genome parts (chromosome, scaffold) that have significant numts
- write individual FASTA files for the genome parts that have significant numts

7_get_lastal_csv

- create a dataframe from the signifcant alignments (score, eg2_value, e_value, g_id, g_start, mt_start, g_length, mt_length, g_strand, mt_strand, g_size, g_sequence, mt_sequence)
- mask the artifacts that are the results of using double mt dna for the alignment
- write dataframe into a csv

8_add_flankings

- add genomic and mitochondrial flanking regions to the previously defined csv file

9_get_ensemble_data

- get the gene id (if available) of the genomic region where a significant numt is inserted
- get the gene description (if available) of the genomic region where a significant numt is inserted

10_chromosome_synteny_data_prep

- create csv file as the input for synteny visualization
- this data is corresponding to the numts which are inserted into chromosomes
- the actual visualization is done in R

11_chromosome_synteny

- synteny visualization of numts that are inserted into chromosomes

12_scaffold_synteny_data_prep

- create csv file as the input for synteny visualization
- this data is corresponding to the numts which are inserted into scaffolds
- the actual visualization is done in R

13_scaffold_synteny

- synteny visualization of numts that are inserted into scaffolds

14_gene_synteny_data_prep

- create csv file as the input for synteny visualization
- this data is corresponding to the numts which are inserted into genes
- the actual visualization is done in R

15_gene_synteny

- synteny visualization of numts that are inserted into genes

16_seq_logo_data_prep

- creating the input file for the sequence logo
- the actual visualization is done in R

17_seq_logo

- sequence logo visualization

18_write_numtless_sequences

- get numt positions
- get the sequences of the genomic parts that contains numts
- write sequences without numts into individual FASTA files

19_gc_content_numts_vs_genome

- sample each genomic parts based on the number and length of the corresponding numts (genomic samples)
- calculate gc contents of the genomic samples
- calculate gc contents of the numts
- comapre the gc contents
- visualize

20_gc_content_flankings_vs_genome

- sample each genomic parts based on the number and length of the corresponding numts (genomic samples)
- calculate gc contents of the genomic samples
- get the flankings of the numts
- calculate flankings gc contents
- comapre the gc contents
- visualize

21_repeatmasker

- sample each genomic parts based on the number and length of the corresponding numts (genomic samples)
- get the flankings of the numts
- prepare input files for RepeatMasker server
- compare the frequency of repetitive elements
- visualize