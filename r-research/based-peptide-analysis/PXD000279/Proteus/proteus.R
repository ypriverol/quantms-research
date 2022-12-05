### peptide

library(proteus)

setwd('D:/dataset/R downstream analysis/0-paper/data_benchmark/0-reviewer/peptide based/PXD000279')

# out_proteus
evi <- read.csv("out_proteus.csv", row.names = NULL)
meta <- read.csv("metadata.csv")
colnames(meta)[1] <- 'experiment'
pepdat <- makePeptideTable(evi, meta)
prodat <- makeProteinTable(pepdat)
prodat.med <- normalizeData(prodat)
res <- limmaDE(prodat.med, sig.level=0.05)

r <- res[which(res$significant), c("protein", "logFC", "adj.P.Val")]
write.csv(res, "proteus-quantms.csv")


# proteinGroups
proteinGroupsFile <- "D:/dataset/R downstream analysis/0-paper/data_benchmark/0-reviewer/peptide based/PXD000279/proteinGroups.txt"

meta <- read.csv("./metadata.csv")
prot.MQ <- readProteinGroups(proteinGroupsFile, meta)
#It is possible to read these data directly into Proteus and skip peptide and protein aggregation steps.
#equalize medians (EM), quantile (Q), no normalization (NN)
prodat.EM <- normalizeData(prot.MQ)
res <- limmaDE(prodat.EM, sig.level=0.05)

r <- res[which(res$significant), c("protein", "logFC", "adj.P.Val")]
write.csv(res, "proteus-maxquant.csv")