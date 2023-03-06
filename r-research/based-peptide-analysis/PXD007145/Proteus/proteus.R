### peptide

library(proteus)

setwd('D:/dataset/R downstream analysis/0-paper/data_benchmark/0-reviewer/peptide based/PXD007145')

### quantms
evi <- read.csv("out_proteus.csv", row.names = NULL)
meta <- read.csv("metadata.csv")
colnames(meta)[1] <- 'experiment'
pepdat <- makePeptideTable(evi, meta)
prodat <- makeProteinTable(pepdat)

# EM (Equalize median) 
prodat.med <- normalizeData(prodat)
write.csv(prodat.med$tab, "proteus-quantms-EM-proteinIntensity.csv")

res1 <- limmaDE(prodat.med, sig.level=0.05, conditions=c("fold4","fold1"))
res2 <- limmaDE(prodat.med, sig.level=0.05, conditions=c("fold10","fold1"))

write.csv(res1, "proteus-quantms-EM-4_1.csv")
write.csv(res2, "proteus-quantms-EM-10_1.csv")


# Q (quantile)
prodat.quant <- normalizeData(prodat, norm.fun=limma::normalizeQuantiles)
write.csv(prodat.quant$tab, "proteus-quantms-Q-proteinIntensity.csv")

res1 <- limmaDE(prodat.quant, sig.level=0.05, conditions=c("fold4","fold1"))
res2 <- limmaDE(prodat.quant, sig.level=0.05, conditions=c("fold10","fold1"))

write.csv(res1, "proteus-quantms-Q-4_1.csv")
write.csv(res2, "proteus-quantms-Q-10_1.csv")

# NN (None)
write.csv(prodat$tab, "proteus-quantms-NN-proteinIntensity.csv")

res1 <- limmaDE(prodat, sig.level=0.05, conditions=c("fold4","fold1"))
res2 <- limmaDE(prodat, sig.level=0.05, conditions=c("fold10","fold1"))

write.csv(res1, "proteus-quantms-NN-4_1.csv")
write.csv(res2, "proteus-quantms-NN-10_1.csv")


##########################################################


### maxquant
evidenceFile <- "D:/dataset/R downstream analysis/0-paper/data_benchmark/0-reviewer/peptide based/PXD007145/evidence_filter.txt"

evi <- readEvidenceFile(evidenceFile)
meta <- read.csv("./metadata-maxquant.csv")
colnames(meta)[1] <- 'experiment'
pepdat <- makePeptideTable(evi, meta)
prodat <- makeProteinTable(pepdat)

# EM (Equalize median) 
prodat.med <- normalizeData(prodat)
res1 <- limmaDE(prodat.med, sig.level=0.05, conditions=c("fold4","fold1"))
res2 <- limmaDE(prodat.med, sig.level=0.05, conditions=c("fold10","fold1"))

write.csv(res1, "proteus-maxquant-EM-4_1.csv")
write.csv(res2, "proteus-maxquant-EM-10_1.csv")

# Q (quantile)
prodat.quant <- normalizeData(prodat, norm.fun=limma::normalizeQuantiles)
res1 <- limmaDE(prodat.quant, sig.level=0.05, conditions=c("fold4","fold1"))
res2 <- limmaDE(prodat.quant, sig.level=0.05, conditions=c("fold10","fold1"))

write.csv(res1, "proteus-maxquant-Q-4_1.csv")
write.csv(res2, "proteus-maxquant-Q-10_1.csv")

# NN (None)
res1 <- limmaDE(prodat, sig.level=0.05, conditions=c("fold4","fold1"))
res2 <- limmaDE(prodat, sig.level=0.05, conditions=c("fold10","fold1"))

write.csv(res1, "proteus-maxquant-NN-4_1.csv")
write.csv(res2, "proteus-maxquant-NN-10_1.csv")
