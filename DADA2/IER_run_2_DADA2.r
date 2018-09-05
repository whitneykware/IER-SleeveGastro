library(dada2); packageVersion("dada2")


# File parsing
pathF <- "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/fwd"
pathR <- "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/rev"
filtpathF <- file.path(pathF, "filtered") 
filtpathR <- file.path(pathR, "filtered") 
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(120,140), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, trimLeft = c(20, 18))

# File parsing
filtpathF <- "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/fwd/filtered" 
filtpathR <- "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/rev/filtered"
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
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/IER_run_2_seqtab.rds") # CHANGE ME to where you want sequence table saved

# Merge multiple runs (if necessary)
st1 <- readRDS("/Users/whitneyware/IER_Sleeve/IER_Run_1_seqs/IER_run_1_seqtab.rds")
st2 <- readRDS("/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/IER_run_2_seqtab.rds")
st.all <- mergeSequenceTables(st1, st2)

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/Users/whitneyware/IER_Sleeve/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
species <- addSpecies(tax, "/Users/whitneyware/IER_Sleeve/silva_species_assignment_v132.fa.gz")

# Write to disk
saveRDS(seqtab, "/Users/whitneyware/IER_Sleeve/IER_combined_seqtab_final.rds") 
saveRDS(tax, "/Users/whitneyware/IER_Sleeve/IER_combined_tax_final.rds") 
saveRDS(species, "/Users/whitneyware/IER_Sleeve/IER_combined_species_final.rds")

# Check species
species.print <- species # Removing sequence rownames for display only
rownames(species.print) <- NULL
head(species.print)

