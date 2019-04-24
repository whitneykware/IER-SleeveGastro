#input from DADA2

setwd("/Users/whitneyware/IER_Sleeve/DADA2_analysis")
myT <- read.table("dada2_merged_pVals.txt", header=TRUE,sep="\t")

plot(myT$Sleeve_DIO_FOb, myT$IER_DIO_FOb,pch=15,cex=1.2,xlab="Sleeve dio vs FoB", ylab= "ier dio vs Fob")