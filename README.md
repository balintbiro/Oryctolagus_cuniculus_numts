# Oryctolagus_cuniculus_numts
 Patterns of numtogenesis in the rabbit (Oryctolagus cuniculus) genome.

This repository contains all the codes for identifying nuclear insertions with mitochondrial origin in the rabbit (Oryctolagus cuniculus) genome.

Setting up the environment:
---
```bash
mkdir codes data results
conda install -c bioconda last
```

The alignments were performed with LASTAL.
The codes are mainly written in Python (3.7.10).

Used Python packages and external programs (if the module is not built in, the version number and conda installation are provided):
---
- os
- gzip
- shutil
- pandas (1.2.4) `conda install pandas`
- urllib.request
- subprocess
- requests (2.26.0) `conda install -c anaconda requests`
- numpy (1.20.2) `conda install numpy`
- seaborn (0.11.2) `conda install seaborn`
- matplotlib (3.4.3) `conda install -c conda-forge matplotlib`
- scipy (1.6.2) `conda install -c anaconda scipy`

To annotate the mitochondrion, MITOS server (http://mitos.bioinf.uni-leipzig.de/index.py) was used with the following setting(s):
- Genetic Code: 02 - Vertebrate

RepeatMasker was also run on server (http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1196065643_I1OZNebCc0toA4hkQLRwVjIcJUTg) with the following setting(s):
- clade: Mammal
- genome: Rabbit
- group: Variation and Repeats
- track: RepeatMasker
- table: rmsk
- define regions: upload the mt BED files
- please note that:
	- the input/output filenames are the followings:
		- repeatmasker_input_downstream.bed --> downstream_repeats.tsv
		- repeatmasker_input_upstream.bead --> upstream_repeats.tsv
		- repeatmasker_input_genomic_samples.bed --> genomic_repeats.tsv
	- for the statistical analysis and visualization, the output files should be placed into the results/ folder

Brief description of the programs:
---

:one:data_acquisition

- download genomic dna and mitochondrial dna in .gz format
- decompress the .gz files

:two:reverse_mt

- actual reversal of the mitochondrial dna

:three:duplicate_mt_dna

- prepare the duplicated mt dna for the later alignment

:four:alignments

- prepare the database with LASTAL from the rabbit genome
- align reverse mt dna with nuclear genome
- align double mt dna with nuclear genome

:five:significant_alignments

- get the lowest e value from the reverse mt dna alignment
- use this threshold for masking the double mt dna alignment

:six:writing_fasta_files

- process nuclear dna file into a format where the sequence is in one line
- get the ids of the genome parts (chromosome, scaffold) that have significant numts
- write individual FASTA files for the genome parts that have significant numts

:seven:get_lastal_csv

- create a dataframe from the signifcant alignments (score, eg2_value, e_value, g_id, g_start, mt_start, g_length, mt_length, g_strand, mt_strand, g_size, g_sequence, mt_sequence)
- mask the artifacts that are the results of using double mt dna for the alignment
- write dataframe into a csv

:eight:add_flankings

- add genomic and mitochondrial flanking regions to the previously defined csv file

:nine:get_ensemble_data

- get the gene id (if available) of the genomic region where a significant numt is inserted
- get the gene description (if available) of the genomic region where a significant numt is inserted

:ten:write_numtless_sequences

- get numt positions
- get the sequences of the genomic parts that contains numts
- write sequences without numts into individual FASTA files

:eleven:gc_content_numts_vs_genome

- sample each genomic parts based on the number and length of the corresponding numts (genomic samples)
- calculate gc contents of the genomic samples
- calculate gc contents of the numts
- comapre the gc contents
- write the output for visualisation

:twelve:gc_content_flankings_vs_genome

- sample each genomic parts based on the number and length of the corresponding numts (genomic samples)
- calculate gc contents of the genomic samples
- get the flankings of the numts
- calculate flankings gc contents
- comapre the gc contents
- write output

:thirteen:repeatmasker

- sample each genomic parts based on the number and length of the corresponding numts (genomic samples)
- get the flankings of the numts
- prepare input files for RepeatMasker server
	- repeatmasker_input_downstream.bed
	- repeatmasker_input_upstream.bed
	- repeatmsker_input_genomic_samples.bed
- the RepeatMasker output files should be named properly (downstream_repeats.tsv, upstream_repeats.tsv and genomic_repeats.tsv) and placed into the results/ folder!
- compare the frequency of repetitive elements
- write output

:warning:DISCLAIMER! Codes for generating figures are not updated recently!