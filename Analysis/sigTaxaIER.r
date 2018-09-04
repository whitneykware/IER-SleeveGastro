setwd('/Users/whitneyware/IER_Sleeve')

rm(list=ls())
SG <- read.table('SleeveGastro_genus_Model_noMDS.txt', sep="\t", header=TRUE)	
IER <- read.table('IER_genus_Model_noMDS.txt', sep = '\t', header = TRUE, row.names = 1)
both <- read.table('merged_pVals.txt', sep = '\t', header = TRUE)
orderBoth <- both[order(both$IERGroupPValue), ]


sigFObvsDIO <- IER[(IER$pValuesFobVsDIOAdjusted < 0.05),]
sigFObvsDIO <- cbind(row.names(sigFObvsDIO), sigFObvsDIO$pValuesFobVsDIOAdjusted)
sigFObvsHCCR <- row.names(IER[(IER$pValuesFobVsHCCRAdjusted < 0.05),])
sigFObvsICCR <- row.names(IER[(IER$pValuesFobVsHCCRAdjusted < 0.05),])
sigHCCRvsICCR <- row.names(IER[(IER$pValuesHCCRvsICCRAdjusted < 0.05),])
sigFObvsICR <- row.names(IER[(IER$pValuesFobVsICRAdjusted < 0.05),])


overlappingTaxa <- Reduce(intersect, list(sigFObvsDIO, sigFObvsHCCR, sigFObvsICCR, sigFObvsICR))