#importing the required libraries
library(Biostrings)
library(TFBSTools)

#reading the file into a vector
numt_ends<-scan('../results/numt_g_ends.fa', what="")

numt_ends<-DNAStringSet(numt_ends)

#create consensus
consensusMatrix(numt_ends, baseOnly=TRUE)

#position frequency matrix
pfm.count<-consensusMatrix(numt_ends, baseOnly=T)[1:4,]
pfm<-PFMatrix(name='numt_ends', profileMatrix=pfm.count)

#information content matrix
icm<-toICM(pfm, pseudocounts=0)

#create sequence logo; ic.scale refers to information content scale in bits
seqLogo(icm, ic.scale=F)

