#authentication
login as: birobalint
birobalint@193.225.92.115's password:
Last login: Mon May  3 12:06:15 2021 from brutus.abc.hu

#navigation
birobalint@login:~$ cd numt
birobalint@login:~/numt$ cd sequences

#content checking
birobalint@login:~/numt/sequences$ ls
double_mt_alignment.fasta
ocdb.bck
ocdb.des
ocdb.prj
ocdb.sds
ocdb.ssp
ocdb.suf
ocdb.tis
Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.MT.fa
Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa
Reversed_Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.MT.fa
reverse_mt_alignment.fasta
two_ch_MT.fasta

birobalint@login:~/numt/sequences$ lastal -r1 -q1 -a7 -b1 ocdb two_ch_MT.fasta > double_mt_alignment.fasta #align double mtDNA. Scores = match, mismatch, gap existence, gap extension
birobalint@login:~/numt/sequences$