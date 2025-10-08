#setwd to alignment folder
setwd()

#install packages required
install.packages("evobiR")
install.packages("ape")
install.packages("msaR")
install.packages("phangorn")

#load the packages 
library(evobiR)
library(ape)
library(msaR)
library(phangorn)

#read in our alignment data
COI<-read.FASTA("COI alignment")
cytb<-read.FASTA("cytb_ptr_COI_alignment")
ptr<-read.FASTA("ptr&COI alignment")
length(COI)
length(unique(names(COI)))

#Look at the MSA
#DO THIS WITH EACH GENE FIRST TO MAKE SURE THE DATA IS GOOD TO USE
#BEFORE making the super matrix
msaR(COI)
msaR(cytb)
msaR(ptr)

#create a new object to not write over og file
#phyDat format
COI2<-phyDat(COI, type="DNA")
cytb2<-phyDat(cytb, type="DNA")
ptr2<-phyDat(ptr, type="DNA")

#trimming
COI3<-COI2[, colMeans(as.character(COI2)=="-") < 0.8]
cytb3<-cytb2[, colMeans(as.character(cytb2)=="-") < 0.8]
ptr3<-ptr2[, colMeans(as.character(ptr2)=="-") < 0.8]

#put back into format MSA viewer can read
#view MSA
COI4<-as.DNAbin(COI3)
msaR(COI4)
cytb4<-as.DNAbin(cytb3)
msaR(cytb4)
ptr4<-as.DNAbin(ptr3)
msaR(ptr4)

#write new trimmed data out to FASTA file
write.FASTA(COI4, "COI_trimmed.fasta")
write.FASTA(cytb4, "cytb_trimmed.fasta")
write.FASTA(ptr4, "ptr_trimmed.fasta")
