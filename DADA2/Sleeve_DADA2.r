library(dada2); packageVersion("dada2")


# File parsing
pathF <- "/Users/whitneyware/IER_Sleeve/sleeve_seqs/fwd"
pathR <- "/Users/whitneyware/IER_Sleeve/sleeve_seqs/rev"
filtpathF <- file.path(pathF, "filtered") 
filtpathR <- file.path(pathR, "filtered") 
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(210,160), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, trimLeft = c(20, 18))


# File parsing
filtpathF <- "/Users/whitneyware/IER_Sleeve/sleeve_seqs/fwd/filtered" 
filtpathR <- "/Users/whitneyware/IER_Sleeve/sleeve_seqs/rev/filtered"
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "-"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs), "-"), `[`, 1) 
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

# Learn forward error rates
errF <- learnErrors(filtFs)
# Learn reverse error rates
errR <- learnErrors(filtRs)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
seqtab <- removeBimeraDenovo(seqtab, method="consensus")

# Assign taxonomy and species
tax <- assignTaxonomy(seqtab, "/Users/whitneyware/IER_Sleeve/silva_nr_v132_train_set.fa.gz")
species <- addSpecies(tax, "/Users/whitneyware/IER_Sleeve/silva_species_assignment_v132.fa.gz")

# Write to disk
saveRDS(seqtab, "/Users/whitneyware/IER_Sleeve/sleeve_seqs/sleeve_seqtab_final.rds")
saveRDS(tax, "/Users/whitneyware/IER_Sleeve/sleeve_seqs/sleeve_tax_final.rds") 
saveRDS(species, "/Users/whitneyware/IER_Sleeve/sleeve_seqs/sleeve_species_final.rds") 

# Read in RDS if needed
seqtab <- readRDS("/Users/whitneyware/IER_Sleeve/sleeve_seqs/sleeve_seqtab_final.rds")
tax <- readRDS("/Users/whitneyware/IER_Sleeve/sleeve_seqs/sleeve_tax_final.rds")
species <- readRDS("/Users/whitneyware/IER_Sleeve/sleeve_seqs/sleeve_species_final.rds")

# Save to table
setwd('/Users/whitneyware/IER_Sleeve/sleeve_seqs')
write.table(seqtab, file="sleeve_septab.txt")
write.table(tax, file="sleeve_tax_tab.txt", sep = '\t')
write.table(species, file="sleeve_species_tab.txt", sep = '\t')

# Check species
species.print <- species # Removing sequence rownames for display only
rownames(species.print) <- NULL
head(species.print)


