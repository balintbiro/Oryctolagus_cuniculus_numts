#authentication
login as: birobalint
birobalint@193.225.92.115's password:
Last login: Fri Apr 30 14:06:40 2021 from brutus.abc.hu

#navigation
birobalint@login:~$ cd numt
birobalint@login:~/numt$ cd sequences

#check the content of that folder
birobalint@login:~/numt/sequences$ ls
Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.MT.fa
Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa
Reversed_Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.MT.fa
two_ch_MT.fasta

birobalint@login:~/numt/sequences$ lastdb ocdb Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa  #create "database files" for the alignment
birobalint@login:~/numt/sequences$ lastal ocdb Reversed_Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.MT.fa > reverse_mt_alignment.fasta #align the reversed mt DNA and write the result into a fasta file
birobalint@login:~/numt/sequences$
