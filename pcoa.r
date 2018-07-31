setwd("/Users/whitneyware/IER_Project/IER")

library("vegan")

rm(list=ls())

project <- "SleeveGastro"
taxa <- "genus"
infile <- "SG_genus_taxa.txt"
metaFile <- 'SleeveGastro_Metadata.txt'

myT <- read.table(infile, header=TRUE,sep="\t", row.names = 1)
myT <- myT[ order(row.names(myT)), ]

meta <-read.table(metaFile, header=TRUE,sep="\t")
meta$Sample.ID <-NULL

myPCOA <- capscale(myT~1,distance="bray")

merged <- cbind(myT,meta)
myMerge <- cbind(row.names(myT), merged,myPCOA$CA$u)

#write.table(myPCOA$CA$u, sep="\t", file=paste(project,"_pcoa_",taxa,".txt",sep=""))
#write.table(myPCOA$CA$eig,file=paste(project,"_eigenValues_pcoa_",taxa,".txt", sep=""), sep="\t")	
write.table(myMerge , sep="\t", file=paste(project,"_",taxa,"_pcoa_metaMerged.txt", sep=""),row.names=FALSE)