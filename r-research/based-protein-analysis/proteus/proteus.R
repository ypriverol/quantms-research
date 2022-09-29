library(proteus)

# 把maxquant数据格式转为out_msstats
setwd('D:/dataset/R downstream analysis/0-paper/data_benchmark/PXD000279(proteome benchmark+dynamic range dataset)/maxquant output/1-proteomic')

proteinGroupsFile <- "D:/dataset/R downstream analysis/0-paper/data_benchmark/PXD000279(proteome benchmark+dynamic range dataset)/maxquant output/2-dynamic/proteinGroups.txt"
meta <- read.csv("./metadata.csv")
prot.MQ <- readProteinGroups(proteinGroupsFile, meta)

#It is possible to read these data directly into Proteus and skip peptide and protein aggregation steps.
prodat.med <- normalizeData(prot.MQ)

res <- limmaDE(prodat.med, sig.level=0.05)
r <- res[which(res$significant), c("protein", "logFC", "adj.P.Val")]
write.csv(res, "proteus.csv")