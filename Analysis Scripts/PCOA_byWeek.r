# pcoa for SG by week
# used for pcoa plots

library("vegan")

infile <- '/Users/whitneyware/IER_Sleeve/project_info/metadata/SG_genus_taxa.txt'
metaFile <- '/Users/whitneyware/IER_Sleeve/project_info/metadata/SleeveGastro_Metadata.txt'

myT <- read.table(infile, header=TRUE,sep="\t", row.names = 1)
myT <- myT[ order(row.names(myT)), ]

meta <-read.table(metaFile, header=TRUE,sep="\t")

# PCOA 15 weeks only
wk15 <- cbind(myT,meta)
wk15 <- subset(wk15, wk15$Time.point..weeks. == "15")
wk15_meta <- wk15[, 112:118]
wk15 <- wk15[,-c(112:118)]
wk15_PCOA <- capscale(wk15~1,distance="bray")
myMerge <- cbind(row.names(wk15), wk15_meta, wk15_PCOA$CA$u)
write.table(myMerge , sep="\t", "SG_wk15_PCOA.txt",row.names=FALSE)

# PCOA 25 weeks only
wk25 <- cbind(myT,meta)
wk25 <- subset(wk25, wk25$Time.point..weeks. == "25")
wk25_meta <- wk25[, 112:118]
wk25 <- wk25[,-c(112:118)]
wk25_PCOA <- capscale(wk25~1,distance="bray")
myMerge25 <- cbind(row.names(wk25), wk25_meta, wk25_PCOA$CA$u)
write.table(myMerge25 , sep="\t", "SG_wk25_PCOA.txt",row.names=FALSE)

# pcoa for IER by week
# used for pcoa plots

infile <- '/Users/whitneyware/IER_Sleeve/project_info/metadata/IER_genus_taxa.txt'
metaFile <- '/Users/whitneyware/IER_Sleeve/project_info/metadata/IER_meta.txt'

myT <- read.table(infile, header=TRUE,sep="\t", row.names = 1)
myT <- myT[ order(row.names(myT)), ]

meta <-read.table(metaFile, header=TRUE,sep="\t")

# PCOA 15 weeks only
wk15 <- cbind(myT,meta)
wk15 <- subset(wk15, wk15$Time.point..weeks. == "15")
wk15_meta <- wk15[, 224:233]
wk15 <- wk15[,-c(224:233)]
wk15_PCOA <- capscale(wk15~1,distance="bray")
myMerge <- cbind(row.names(wk15), wk15_meta, wk15_PCOA$CA$u)
write.table(wk15 , sep="\t", "/Users/whitneyware/IER_Sleeve/project_info/metadata/Ier_wk15_PCOA_test.txt",row.names=FALSE)

# PCOA 25 weeks only
wk25 <- cbind(myT,meta)
wk25 <- subset(wk25, wk25$Time.point..weeks. == "25")
wk25_meta <- wk25[, 224:233]
wk25 <- wk25[,-c(224:233)]
wk25_PCOA <- capscale(wk25~1,distance="bray")
myMerge25 <- cbind(row.names(wk25), wk25_meta, wk25_PCOA$CA$u)
write.table(myMerge25 , sep="\t", "SG_wk25_PCOA.txt",row.names=FALSE)
