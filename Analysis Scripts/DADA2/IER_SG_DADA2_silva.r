vals <- read.table("/Users/whitneyware/IER_Sleeve/DADA2_analysis/IER_SG_dada_merged_pVals.txt", header = TRUE, sep = "\t")


vals <- vals[!is.na(vals$ier_pValuesFobVsDIOAdjusted), ]
vals <- vals[!is.na(vals$sle_pValuesFobVsDIOAdjusted), ]

log10IER_F_D<- log10(vals$IER_FOBvsDIO)
log10SG_F_D <- log10(vals$sleeveFOBvsDIO)

mergedPvals <- read.table("/Users/whitneyware/IER_Sleeve/dada2_analysis/dada2_merged_pVals.txt", header = TRUE, sep = "\t")
plot(mergedPvals$IERpValuesLog10VolumeWithDirection, mergedPvals$SleevePValuesLog10VolumeWithDirection,pch=15,cex=1.2,xlab="IER Log10 Volume", ylab= "Sleeve Log10 Volume")


myT <- read.table("dada2_merged_pVals.txt", header=TRUE,sep="\t")
plot(mergedPvals$IERGroupPValue, mergedPvals$IERGroupPValue,pch=15,cex=1.2,xlab="IER group", ylab= "Sleeve group")
plot(log10IER_F_D, log10SG_F_D,pch=15,cex=1.2,xlab="IER group", ylab= "Sleeve group")

library(dada2); packageVersion("dada2")

myT <- read.table("sleeve_DADA2_model.txt", sep="\t", header=TRUE)	

seqs <- myT$Sequences
seqs <- as.character(seqs)
set.seed(100) # Initialize random number generator for reproducibility
taxa <- assignTaxonomy(seqs, "silva_nr_v132_train_set.fa.gz", multithread=FALSE)

write.table(taxa, "sleeve_taxa_silva.txt", sep = "\t")


