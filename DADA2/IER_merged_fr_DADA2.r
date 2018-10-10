library(dada2); packageVersion("dada2")

#IER fwd and rev, merged runs 1 and 2
#did not work due to poor quality

#Run 1
# File parsing
pathF <- "/Users/whitneyware/IER_Sleeve/IER_Run_1_seqs/fwd"
pathR <- "/Users/whitneyware/IER_Sleeve/IER_Run_1_seqs/rev"
filtpathF <- file.path(pathF, "filtered_160_140") 
filtpathR <- file.path(pathR, "filtered_160_140") 
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(160,140), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, trimLeft = c(20,18))

# File parsing
filtpathF <- "/Users/whitneyware/IER_Sleeve/IER_Run_1_seqs/fwd/filtered_160_140" 
filtpathR <- "/Users/whitneyware/IER_Sleeve/IER_Run_1_seqs/rev/filtered_160_140"
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "-"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs), "-"), `[`, 1) 
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

errF <- learnErrors(filtFs)
errR <- learnErrors(filtRs)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

table(nchar(getSequences(seqtab)))
saveRDS(seqtab, "/Users/whitneyware/IER_Sleeve/IER_final/IER_run_1_seqtab_160_140.rds")

#run 2
# File parsing
pathF <- "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/fwd"
pathR <- "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/rev"
filtpathF <- file.path(pathF, "filtered_105_140") 
filtpathR <- file.path(pathR, "filtered_105_140") 
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(105,140), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, trimLeft = c(20,18))

# File parsing
filtpathF <- "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/fwd/filtered_105_140" 
filtpathR <- "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/rev/filtered_105_140"
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "-"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs), "-"), `[`, 1) 
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

errF <- learnErrors(filtFs)
errR <- learnErrors(filtRs)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

seqtab_2 <- makeSequenceTable(mergers)
saveRDS(seqtab_2, "/Users/whitneyware/IER_Sleeve/IER_final/IER_run_2_seqtab_105_140.rds")

seqtab_1 <- readRDS("/Users/whitneyware/IER_Sleeve/IER_final/IER_run_1_seqtab_160_140.rds")
table(nchar(getSequences(seqtab_1)))
table(nchar(getSequences(seqtab_2)))

seqtab_1<-seqtab_1 [,nchar(colnames(seqtab_1)) %in% seq(120,190)]
seqtab_2<-seqtab_2 [,nchar(colnames(seqtab_2)) %in% seq(120,190)]

merged_seqtab <- mergeSequenceTables(seqtab_1, seqtab_2)
merged.nochim <- removeBimeraDenovo(merged_seqtab, method="consensus", multithread=TRUE)

saveRDS(merged.nochim, "/Users/whitneyware/IER_Sleeve/IER_final/IER_seqtab_adjustedSeqRange.rds")

# Assign taxonomy
tax <- assignTaxonomy(merged.nochim, "/Users/whitneyware/IER_Sleeve/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(tax, "/Users/whitneyware/IER_Sleeve/IER_final/IER_tax_adjustedSeqRange.rds")

setwd('/Users/whitneyware/IER_Sleeve/IER_final')
write.table(merged.nochim, file="IER_seqtab_adjustedSeqRange.txt", sep = "\t")
write.table(tax, file="IER_tax_adjustedSeqRange.txt", sep = "\t")

