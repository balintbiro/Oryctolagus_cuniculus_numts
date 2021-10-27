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