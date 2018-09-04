library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

meta <- read.table("/Users/whitneyware/IER_Sleeve/meta/SleeveGastro_Metadata.txt", sep = "\t", header = TRUE, row.names = 1)
meta <- as.data.frame(meta)

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), sample_data(meta), tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample


plot_richness(ps, x="Time.point..weeks.", measures=c("Shannon", "Simpson"), color="Treatment.Group")

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Treatment.Group", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

plot_bar(ps, x="Treatment.Group", fill="Family") + facet_wrap(~Time.point..weeks., scales="free_x")
plot_bar(ps.top20, x="Treatment.Group", fill="Family") + facet_wrap(~Time.point..weeks., scales="free_x")
plot_bar(ps, x="Treatment.Group", fill="Genus") + facet_wrap(~Time.point..weeks., scales="free_x")
plot_bar(ps.top20, x="Treatment.Group", fill="Genus") + facet_wrap(~Time.point..weeks., scales="free_x")


