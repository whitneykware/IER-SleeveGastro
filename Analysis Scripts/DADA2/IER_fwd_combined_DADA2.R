#forward reads only for IER run 1 

library(dada2); packageVersion("dada2")

# Filename parsing
path <- "/Users/whitneyware/IER_Sleeve/IER_Run_1_seqs/fwd"
filtpath <- file.path(path, "filtered_fwdOnly_trimPrimer_110Both_2_11") 
fns <- list.files(path, pattern="fastq.gz") 
# Filtering
filterAndTrim(file.path(path,fns), file.path(filtpath,fns), 
              truncLen=110, maxEE=1, truncQ=11, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, trimLeft = 20)

# File parsing
filtpath <- "/Users/whitneyware/IER_Sleeve/IER_Run_1_seqs/fwd/filtered_fwdOnly_trimPrimer_110Both_2_11"
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) 
sample.names <- sapply(strsplit(basename(filts), "-"), `[`, 1) # Assumes filename = sample_XXX.fastq.gz
names(filts) <- sample.names
# Learn error rates
set.seed(100)
err <- learnErrors(filts, multithread=TRUE, randomize=TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err)
}
# Construct sequence table and write to disk
seqtab_run1 <- makeSequenceTable(dds)
saveRDS(seqtab_run1, "/Users/whitneyware/IER_Sleeve//IER_run_1_fwdOnly_final.rds")


#forward reads only for IER run 2

# Filename parsing
path <- "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/fwd"
filtpath <- file.path(path, "filtered_fwdOnly_trimPrimer_110Both_2_11") 
fns <- list.files(path, pattern="fastq.gz") 
# Filtering
filterAndTrim(file.path(path,fns), file.path(filtpath,fns), 
              truncLen=110, maxEE=1, truncQ=11, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, trimLeft = 20)

# File parsing
filtpath <- "/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/fwd/filtered_fwdOnly_trimPrimer_110Both_2_11"
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) 
sample.names <- sapply(strsplit(basename(filts), "-"), `[`, 1) # Assumes filename = sample_XXX.fastq.gz
names(filts) <- sample.names
# Learn error rates
set.seed(100)
err <- learnErrors(filts, multithread=TRUE, randomize=TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err)
}
# Construct sequence table and write to disk
seqtab_run2 <- makeSequenceTable(dds)
saveRDS(seqtab_run2, "/Users/whitneyware/IER_Sleeve/IER_run_2_fwdOnly_final.rds") 


# Merge multiple runs 
seqtab_run1 <- readRDS("/Users/whitneyware/IER_Sleeve/IER_run_1_fwdOnly_final.rds")
#st2 <- readRDS("/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs/IER_run_2_fwdOnly_trimPrimer_110Both.rds")
st.all <- mergeSequenceTables(seqtab_run1, seqtab_run2)

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/Users/whitneyware/IER_Sleeve/silva_nr_v132_train_set.fa.gz")

# Write to disk
saveRDS(seqtab, "/Users/whitneyware/IER_Sleeve/IER_fwdOnly_final_seqtab.rds") 
saveRDS(tax, "/Users/whitneyware/IER_Sleeve/IER_fwdOnly_final_tax.rds") 

# Save to table
setwd('/Users/whitneyware/IER_Sleeve')
write.table(seqtab, file="IER_fwdOnly_final_seqtab.txt", sep = "\t")
write.table(tax, file="IER_fwdOnly_final_tax.txt", sep = "\t")

