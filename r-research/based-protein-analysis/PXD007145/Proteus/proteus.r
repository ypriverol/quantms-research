#PXD007145
library(proteus)
setwd('D:/dataset/R downstream analysis/0-paper/data_benchmark/0-reviewer/protein based/PXD007145')
proteinGroupsFile <- "D:/dataset/R downstream analysis/0-paper/data_benchmark/0-reviewer/protein based/PXD007145/proteinGroups.txt"

meta <- read.csv("./metadata.csv")
prot.MQ <- readProteinGroups(proteinGroupsFile, meta)

#It is possible to read these data directly into Proteus and skip peptide and protein aggregation steps.
#equalize medians (EM), quantile (Q), no normalization (NN)
prodat.EM <- normalizeData(prot.MQ)
prodat.Q <- normalizeData(prot.MQ, norm.fun=limma::normalizeQuantiles)
prodat.NN <- prot.MQ

# 10:1 fold
#EM
res <- limmaDE(prodat.EM, sig.level=0.05, conditions=c("10","1"))
r <- res[which(res$significant), c("protein", "logFC", "adj.P.Val")]
write.csv(res, "proteus_median_10-1.csv")
#Q
res <- limmaDE(prodat.Q, sig.level=0.05, conditions=c("10","1"))
r <- res[which(res$significant), c("protein", "logFC", "adj.P.Val")]
write.csv(res, "proteus_quantile_10-1.csv")
#NN
res <- limmaDE(prodat.NN, sig.level=0.05, conditions=c("10","1"))
r <- res[which(res$significant), c("protein", "logFC", "adj.P.Val")]
write.csv(res, "proteus_none_10-1.csv")

# 4:1 fold
#EM
res <- limmaDE(prodat.EM, sig.level=0.05, conditions=c("4","1"))
r <- res[which(res$significant), c("protein", "logFC", "adj.P.Val")]
write.csv(res, "proteus_median_4-1.csv")
#Q
res <- limmaDE(prodat.Q, sig.level=0.05, conditions=c("4","1"))
r <- res[which(res$significant), c("protein", "logFC", "adj.P.Val")]
write.csv(res, "proteus_quantile_4-1.csv")
#NN
res <- limmaDE(prodat.NN, sig.level=0.05, conditions=c("4","1"))
r <- res[which(res$significant), c("protein", "logFC", "adj.P.Val")]
write.csv(res, "proteus_none_4-1.csv")