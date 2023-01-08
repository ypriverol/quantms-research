library(tidyverse)
library(limma)
library(QFeatures)
library(msqrob2)
library(plotly)
library(gridExtra)
library(msdata)
library(openxlsx)
library(MSnbase)

setwd('D:/dataset/R downstream analysis/0-paper/data_benchmark/0-reviewer/peptide based')

### MaxQuant

## CM (center.median)
MQ_pep = read.csv("./PXD007145/peptides_filter.txt", sep = "\t")
ecols <- grep("LFQ.intensity\\.", names(MQ_pep))

pe <- readQFeatures(
  table = MQ_pep, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)

colData(pe)$condition <- as.factor(c("1", "1", "1", "1", "1", "1", "10", "10", 
                                     "10", "10", "10", "10", "4", "4", "4", "4",
                                     "4", "4"))


rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
  xlab("Peptide index (ordered by data completeness)")

pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
limma::plotDensities(assay(pe[["peptideLog"]]))

Protein_filter <- rowData(pe[["peptideLog"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
pe <- pe[Protein_filter,,]

pe <- filterFeatures(pe, ~ Reverse != "+")
pe <- filterFeatures(pe, ~ Potential.contaminant != "+")


pe <- filterFeatures(pe, ~ nNonZero >= 2)

pe <- normalize(pe,
                i = "peptideLog",
                name = "peptideNorm",
                method = "center.median"
)

limma::plotDensities(assay(pe[["peptideNorm"]]))
boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distribtutions after normalisation", ylab = "intensity"
)

limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))

# Summarization to protein level
pe <- aggregateFeatures(pe,
                        i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
                        name = "protein"
)

plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[2]])


L4 <- makeContrast("condition4=0", parameterNames = c("condition4"))
pe4 <- hypothesisTest(object = pe, i = "protein", contrast = L4)


volcano <- ggplot(
  rowData(pe4[["protein"]])$condition4,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe4[["protein"]])$condition4, file = "PXD007145-4FC-maxquant-dep-CM.csv")


L10 <- makeContrast("condition10=0", parameterNames = c("condition10"))
pe10 <- hypothesisTest(object = pe, i = "protein", contrast = L10)

volcano <- ggplot(
  rowData(pe10[["protein"]])$condition10,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe10[["protein"]])$condition10, file = "PXD007145-10FC-maxquant-dep-CM.csv")



## Q (quantile)
MQ_pep = read.csv("./PXD007145/peptides_filter.txt", sep = "\t")
ecols <- grep("LFQ.intensity\\.", names(MQ_pep))

pe <- readQFeatures(
  table = MQ_pep, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)

colData(pe)$condition <- as.factor(c("1", "1", "1", "1", "1", "1", "10", "10", 
                                     "10", "10", "10", "10", "4", "4", "4", "4",
                                     "4", "4"))


rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
  xlab("Peptide index (ordered by data completeness)")

pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
limma::plotDensities(assay(pe[["peptideLog"]]))

Protein_filter <- rowData(pe[["peptideLog"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
pe <- pe[Protein_filter,,]

pe <- filterFeatures(pe, ~ Reverse != "+")
pe <- filterFeatures(pe, ~ Potential.contaminant != "+")


pe <- filterFeatures(pe, ~ nNonZero >= 2)

pe <- normalize(pe,
                i = "peptideLog",
                name = "peptideNorm",
                method = "quantiles"
)

limma::plotDensities(assay(pe[["peptideNorm"]]))
boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distribtutions after normalisation", ylab = "intensity"
)

limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))

# Summarization to protein level
pe <- aggregateFeatures(pe,
                        i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
                        name = "protein"
)

plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[2]])


L4 <- makeContrast("condition4=0", parameterNames = c("condition4"))
pe4 <- hypothesisTest(object = pe, i = "protein", contrast = L4)


volcano <- ggplot(
  rowData(pe4[["protein"]])$condition4,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe4[["protein"]])$condition4, file = "PXD007145-4FC-maxquant-dep-Q.csv")


L10 <- makeContrast("condition10=0", parameterNames = c("condition10"))
pe10 <- hypothesisTest(object = pe, i = "protein", contrast = L10)

volcano <- ggplot(
  rowData(pe10[["protein"]])$condition10,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe10[["protein"]])$condition10, file = "PXD007145-10FC-maxquant-dep-Q.csv")




## NN (none)
MQ_pep = read.csv("./PXD007145/peptides_filter.txt", sep = "\t")
ecols <- grep("LFQ.intensity\\.", names(MQ_pep))

pe <- readQFeatures(
  table = MQ_pep, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)

colData(pe)$condition <- as.factor(c("1", "1", "1", "1", "1", "1", "10", "10", 
                                     "10", "10", "10", "10", "4", "4", "4", "4",
                                     "4", "4"))


rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
  xlab("Peptide index (ordered by data completeness)")

pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
limma::plotDensities(assay(pe[["peptideLog"]]))

Protein_filter <- rowData(pe[["peptideLog"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
pe <- pe[Protein_filter,,]

pe <- filterFeatures(pe, ~ Reverse != "+")
pe <- filterFeatures(pe, ~ Potential.contaminant != "+")


pe <- filterFeatures(pe, ~ nNonZero >= 2)

#pe <- normalize(pe,
#                i = "peptideLog",
#                name = "peptideNorm",
#                method = "center.median"
#)
pe <- addAssay(pe,pe[["peptideLog"]],"peptideNorm")

limma::plotDensities(assay(pe[["peptideNorm"]]))
boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distribtutions after normalisation", ylab = "intensity"
)

limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))

# Summarization to protein level
pe <- aggregateFeatures(pe,
                        i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
                        name = "protein"
)

plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[2]])


L4 <- makeContrast("condition4=0", parameterNames = c("condition4"))
pe4 <- hypothesisTest(object = pe, i = "protein", contrast = L4)


volcano <- ggplot(
  rowData(pe4[["protein"]])$condition4,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe4[["protein"]])$condition4, file = "PXD007145-4FC-maxquant-dep-NN.csv")


L10 <- makeContrast("condition10=0", parameterNames = c("condition10"))
pe10 <- hypothesisTest(object = pe, i = "protein", contrast = L10)

volcano <- ggplot(
  rowData(pe10[["protein"]])$condition10,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe10[["protein"]])$condition10, file = "PXD007145-10FC-maxquant-dep-NN.csv")



###################################################################



### quantms

## CM (center.median)
mzTab_pep = read.csv("./PXD007145/only_pep_filterPXD007145-Th.sdrf_openms_design_openms.mzTab", sep="\t")
#mzTab_pep = mzTab_pep[mzTab_pep$opt_global_cv_MS.1002217_decoy_peptide == 0,]
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                         "sumIntensity_3")]
mzTab_pep_sum <- mzTab_pep %>% group_by(sequence, accession) %>% summarise(sumIntensity_1 = sum(sumIntensity_1),
                                                                           sumIntensity_2 = sum(sumIntensity_2),
                                                                           sumIntensity_3 = sum(sumIntensity_3))

mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3]) > 0)&(rowSums(mzTab_pep_sum[,4]) > 0)&(rowSums(mzTab_pep_sum[,5]) > 0)),]

names(mzTab_pep_sum)[names(mzTab_pep_sum) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep_sum))

pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("1", "4", "10"))

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
  xlab("Peptide index (ordered by data completeness)")

pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
limma::plotDensities(assay(pe[["peptideLog"]]))

Protein_filter <- rowData(pe[["peptideLog"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
pe <- pe[Protein_filter,,]


pe <- filterFeatures(pe, ~ nNonZero >= 2)


pe <- normalize(pe,
                i = "peptideLog",
                name = "peptideNorm",
                method = "center.median"
)


limma::plotDensities(assay(pe[["peptideNorm"]]))
boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distribtutions after normalisation", ylab = "intensity"
)

limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))


# Summarization to protein level
pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")


plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

L4 <- makeContrast("condition4=0", parameterNames = c("condition4"))
pe4 <- hypothesisTest(object = pe, i = "protein", contrast = L4)


volcano <- ggplot(
  rowData(pe4[["protein"]])$condition4,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe4[["protein"]])$condition4, file = "PXD007145-4FC-quantms-dep-CM.csv")


L10 <- makeContrast("condition10=0", parameterNames = c("condition10"))
pe10 <- hypothesisTest(object = pe, i = "protein", contrast = L10)

volcano <- ggplot(
  rowData(pe10[["protein"]])$condition10,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe10[["protein"]])$condition10, file = "PXD007145-10FC-quantms-dep-CM.csv")




## Q (quantile)
mzTab_pep = read.csv("./PXD007145/only_pep_filterPXD007145-Th.sdrf_openms_design_openms.mzTab", sep="\t")
mzTab_pep = mzTab_pep[mzTab_pep$opt_global_cv_MS.1002217_decoy_peptide == 0,]
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                         "sumIntensity_3")]
mzTab_pep_sum <- mzTab_pep %>% group_by(sequence, accession) %>% summarise(sumIntensity_1 = sum(sumIntensity_1),
                                                                           sumIntensity_2 = sum(sumIntensity_2),
                                                                           sumIntensity_3 = sum(sumIntensity_3))

mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3:6]) > 0)|(rowSums(mzTab_pep_sum[7:10]) > 0)),]

names(mzTab_pep)[names(mzTab_pep) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep))

pe <- readQFeatures(
  table = mzTab_pep, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("1", "4", "10"))

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
  xlab("Peptide index (ordered by data completeness)")

pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
limma::plotDensities(assay(pe[["peptideLog"]]))

Protein_filter <- rowData(pe[["peptideLog"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
pe <- pe[Protein_filter,,]


pe <- filterFeatures(pe, ~ nNonZero >= 2)


pe <- normalize(pe,
                i = "peptideLog",
                name = "peptideNorm",
                method = "quantiles"
)


limma::plotDensities(assay(pe[["peptideNorm"]]))
boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distribtutions after normalisation", ylab = "intensity"
)

limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))


# Summarization to protein level
pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")


plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

L4 <- makeContrast("condition4=0", parameterNames = c("condition4"))
pe4 <- hypothesisTest(object = pe, i = "protein", contrast = L4)


volcano <- ggplot(
  rowData(pe4[["protein"]])$condition4,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe4[["protein"]])$condition4, file = "PXD007145-4FC-quantms-dep-Q.csv")


L10 <- makeContrast("condition10=0", parameterNames = c("condition10"))
pe10 <- hypothesisTest(object = pe, i = "protein", contrast = L10)

volcano <- ggplot(
  rowData(pe10[["protein"]])$condition10,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe10[["protein"]])$condition10, file = "PXD007145-10FC-quantms-dep-Q.csv")





## NN (none)
mzTab_pep = read.csv("./PXD007145/only_pep_filterPXD007145-Th.sdrf_openms_design_openms.mzTab", sep="\t")
mzTab_pep = mzTab_pep[mzTab_pep$opt_global_cv_MS.1002217_decoy_peptide == 0,]
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                         "sumIntensity_3")]
mzTab_pep_sum <- mzTab_pep %>% group_by(sequence, accession) %>% summarise(sumIntensity_1 = sum(sumIntensity_1),
                                                                           sumIntensity_2 = sum(sumIntensity_2),
                                                                           sumIntensity_3 = sum(sumIntensity_3))

mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3:6]) > 0)|(rowSums(mzTab_pep_sum[7:10]) > 0)),]

names(mzTab_pep)[names(mzTab_pep) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep))

pe <- readQFeatures(
  table = mzTab_pep, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("1", "4", "10"))

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
  xlab("Peptide index (ordered by data completeness)")

pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
limma::plotDensities(assay(pe[["peptideLog"]]))

Protein_filter <- rowData(pe[["peptideLog"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
pe <- pe[Protein_filter,,]


pe <- filterFeatures(pe, ~ nNonZero >= 2)


#pe <- normalize(pe,
#                i = "peptideLog",
#                name = "peptideNorm",
#                method = "none"
#)
pe <- addAssay(pe,pe[["peptideLog"]],"peptideNorm")


limma::plotDensities(assay(pe[["peptideNorm"]]))
boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distribtutions after normalisation", ylab = "intensity"
)

limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))


# Summarization to protein level
pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")


plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

L4 <- makeContrast("condition4=0", parameterNames = c("condition4"))
pe4 <- hypothesisTest(object = pe, i = "protein", contrast = L4)


volcano <- ggplot(
  rowData(pe4[["protein"]])$condition4,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe4[["protein"]])$condition4, file = "PXD007145-4FC-quantms-dep-NN.csv")


L10 <- makeContrast("condition10=0", parameterNames = c("condition10"))
pe10 <- hypothesisTest(object = pe, i = "protein", contrast = L10)

volcano <- ggplot(
  rowData(pe10[["protein"]])$condition10,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe10[["protein"]])$condition10, file = "PXD007145-10FC-quantms-dep-NN.csv")















## test
mzTab_pep = read.csv("./PXD007145/only_pep_filterPXD007145-Th.sdrf_openms_design_openms.mzTab", sep="\t")
#mzTab_pep = mzTab_pep[mzTab_pep$opt_global_cv_MS.1002217_decoy_peptide == 0,]
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                         "sumIntensity_3")]
mzTab_pep_sum <- mzTab_pep %>% group_by(sequence, accession) %>% summarise(sumIntensity_1 = sum(sumIntensity_1),
                                                                           sumIntensity_2 = sum(sumIntensity_2),
                                                                           sumIntensity_3 = sum(sumIntensity_3))

#mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3]) > 0)|(rowSums(mzTab_pep_sum[,4]) > 0)|(rowSums(mzTab_pep_sum[,5]) > 0)),]
mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3]) > 0)&(rowSums(mzTab_pep_sum[,4]) > 0)&(rowSums(mzTab_pep_sum[,5]) > 0)),]


names(mzTab_pep_sum)[names(mzTab_pep_sum) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep_sum))

pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("1", "4", "10"))

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
  xlab("Peptide index (ordered by data completeness)")

pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
limma::plotDensities(assay(pe[["peptideLog"]]))

Protein_filter <- rowData(pe[["peptideLog"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
pe <- pe[Protein_filter,,]


pe <- filterFeatures(pe, ~ nNonZero >= 2)


pe <- normalize(pe,
                i = "peptideLog",
                name = "peptideNorm",
                method = "center.median"
)


limma::plotDensities(assay(pe[["peptideNorm"]]))
boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distribtutions after normalisation", ylab = "intensity"
)

limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))


# Summarization to protein level
pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")


plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

L4 <- makeContrast("condition4=0", parameterNames = c("condition4"))
pe4 <- hypothesisTest(object = pe, i = "protein", contrast = L4)


volcano <- ggplot(
  rowData(pe4[["protein"]])$condition4,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe4[["protein"]])$condition4, file = "PXD007145-4FC-quantms-dep-CM.csv")


L10 <- makeContrast("condition10=0", parameterNames = c("condition10"))
pe10 <- hypothesisTest(object = pe, i = "protein", contrast = L10)

volcano <- ggplot(
  rowData(pe10[["protein"]])$condition10,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe10[["protein"]])$condition10, file = "PXD007145-10FC-quantms-dep-CM.csv")
