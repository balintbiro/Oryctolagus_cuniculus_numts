#installing the RIdeogram package (just for the first time)
#install.packages("RIdeogram")
#load the required module
require(RIdeogram)
#get the karyotype header
karyotype_header<-read.table('../results/scaffold_synteny_dual_karyotype.csv',
                             nrows=1,
                             header=FALSE,
                             sep=',')
#get the karyotype dataframe
karyotype_dual_comparison<-read.table('../results/scaffold_synteny_dual_karyotype.csv',
                                      header=FALSE,
                                      sep=',',
                                      skip=1)
#add header to the dataframe
colnames(karyotype_dual_comparison)<-karyotype_header

#get the synteny header
synteny_header<-read.table('../results/scaffold_synteny_dual_synteny.csv',
                           nrows=1,
                           header=FALSE,
                           sep=',')
#get the synteny dataframe
synteny_dual_comparison<-read.table('../results/scaffold_synteny_dual_synteny.csv',
                                    header=FALSE,
                                    sep=',',
                                    skip=1)
#add header to synteny dataframe
colnames(synteny_dual_comparison)<-synteny_header
ideogram(karyotype = karyotype_dual_comparison,
         synteny = synteny_dual_comparison)
#set the working directory
setwd('../results/')
convertSVG("chromosome.svg", device = "png", dpi=600)#the svg filename cannot be modified, so it has to be changed after the creation,
#otherwise it is gonna be overwritten!
