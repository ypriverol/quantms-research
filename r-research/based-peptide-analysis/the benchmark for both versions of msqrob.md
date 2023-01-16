# MSqRob
We use MSqRob 0.7.7 version, and use the following codes to start the shiny program
```r
library(MSqRob)
shiny::runApp(system.file("App-MSqRob", package="MSqRob"))
```
## PXD000279 dataset
### Input type: MaxQuant
<br>

**The options that are not described are default.**
- [ Input ] tab
    - Input type: MaxQuant 
    - Peptides file: peptides_filter.txt (could be downloaded from https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/peptides_filter.txt)
    - Annotation file: msqrob-maxquant_experimental_annotation.xlsx (could be downloaded from https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/msqrob-maxquant_experimental_annotation.xlsx)

- [ Preprocessing ] tab 
    - Default options

        - Group by: Proteins 

        - √ Log-transform data 

        - Base: 2 

        - √ Remove only identified by site (upload proteinGroups.txt file, could be downloaded from https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/proteinGroups.txt) 

        - Minimum number of peptides: 2 

        - Filter columns: Reverse, Potential.contaminant 

    - Changed options 

        - Normalization: center.median / quantiles / none 

- [ Quantification ] tab 

    - Fixed effects: treatment 

    - Analysis type: standard 
    - Number of contrasts: 1
    - Contrast annotation: UPS1 = 1, UPS2 = -1 


### Input type: quantms
<br>

**The options that are not described are default.**
- [ Input ] tab
    - Input type: mzTab
    - Peptides file: onlyPEP-filter-PXD000279.dynamic.sdrf_openms_design_openms.mzTab (could be downloaded from https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/onlyPEP-filter-PXD000279.dynamic.sdrf_openms_design_openms.mzTab)
    - Annotation file: msqrob-quantms_experimental_annotation.xlsx (could be downloaded from https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/msqrob-quantms_experimental_annotation.xlsx)

- [ Preprocessing ] tab 
    - Default options

        - Group by: accession

        - √ Log-transform data 

        - Base: 2 

        - Minimum number of peptides: 2 

    - Changed options 

        - Normalization: center.median / quantiles / none 

- [ Quantification ] tab 

    - Fixed effects: treatment 

    - Analysis type: standard 
    - Number of contrasts: 1
    - Contrast annotation: UPS1 = 1, UPS2 = -1 

<br>

## PXD007145 dataset
### Input type: MaxQuant
<br>

**The options that are not described are default.**
- [ Input ] tab
    - Input type: MaxQuant 
    - Peptides file: peptides_filter.txt (could be downloaded from https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145/peptides_filter.txt)
    - Annotation file: msqrob-maxquant-experimentalDesign.xlsx (could be downloaded from https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145/msqrob-maxquant-experimentalDesign.xlsx)

- [ Preprocessing ] tab 
    - Default options

        - Group by: Proteins 

        - √ Log-transform data 

        - Base: 2 

        - √ Remove only identified by site (upload proteinGroups.txt file, could be downloaded from https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145/proteinGroups.txt) 

        - Minimum number of peptides: 2 

        - Filter columns: Reverse, Potential.contaminant 

    - Changed options 

        - Normalization: center.median / quantiles / none 

- [ Quantification ] tab 

    - Fixed effects: treatment 

    - Analysis type: standard 
    - Number of contrasts: 2
    - Contrast annotation: (1)fold1 = 1, fold4 = -1; (2)fold1 = 1, fold10 = -1 


### Input type: quantms
<br>

**The options that are not described are default.**
- [ Input ] tab
    - Input type: mzTab
    - Peptides file: only_pep_filterPXD007145-Th.sdrf_openms_design_openms.mzTab (could be downloaded from https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145/only_pep_filterPXD007145-Th.sdrf_openms_design_openms.mzTab)
    - Annotation file: msqrob-quantms_experimental_annotation.xlsx (could be downloaded from https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145/msqrob-quantms_experimental_annotation.xlsx)

- [ Preprocessing ] tab 
    - Default options

        - Group by: accession

        - √ Log-transform data 

        - Base: 2 

        - Minimum number of peptides: 2 

    - Changed options 

        - Normalization: center.median / quantiles / none 

- [ Quantification ] tab 

    - Fixed effects: treatment 

    - Analysis type: standard 
    - Number of contrasts: 2
    - Contrast annotation: (1)fold1 = 1, fold4 = -1; (2)fold1 = 1, fold10 = -1 

<br>
<hr>

# msqrob2
We use msqrob2 1.5.3 version, and both datasets are analyzed using the R scripts.

## PXD000279 dataset
```r
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
MQ_pep = read.csv("./PXD000279/peptides_filter.txt", sep = "\t")
ecols <- grep("LFQ.intensity\\.", names(MQ_pep))

pe <- readQFeatures(
  table = MQ_pep, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)

colData(pe)$condition <- as.factor(c("UPS1", "UPS1", "UPS1", "UPS1", "UPS2", "UPS2", "UPS2", "UPS2"))


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
getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

L <- makeContrast("conditionUPS2=0", parameterNames = c("conditionUPS2"))
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

volcano <- ggplot(
  rowData(pe[["protein"]])$conditionUPS2,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe[["protein"]])$conditionUPS2, file = "ups-maxquant-dep-CM.csv")




## Q (quantile)
MQ_pep = read.csv("./PXD000279/peptides_filter.txt", sep = "\t")
ecols <- grep("LFQ.intensity\\.", names(MQ_pep))

pe <- readQFeatures(
  table = MQ_pep, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)

colData(pe)$condition <- as.factor(c("UPS1", "UPS1", "UPS1", "UPS1", "UPS2", "UPS2", "UPS2", "UPS2"))


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
getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

L <- makeContrast("conditionUPS2=0", parameterNames = c("conditionUPS2"))
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

volcano <- ggplot(
  rowData(pe[["protein"]])$conditionUPS2,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe[["protein"]])$conditionUPS2, file = "ups-maxquant-dep-Q.csv")




## NN (none)
MQ_pep = read.csv("./PXD000279/peptides_filter.txt", sep = "\t")
ecols <- grep("LFQ.intensity\\.", names(MQ_pep))

pe <- readQFeatures(
  table = MQ_pep, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)

colData(pe)$condition <- as.factor(c("UPS1", "UPS1", "UPS1", "UPS1", "UPS2", "UPS2", "UPS2", "UPS2"))


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
pe <- aggregateFeatures(pe,
                        i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
                        name = "protein"
)

plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)
getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

L <- makeContrast("conditionUPS2=0", parameterNames = c("conditionUPS2"))
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

volcano <- ggplot(
  rowData(pe[["protein"]])$conditionUPS2,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe[["protein"]])$conditionUPS2, file = "ups-maxquant-dep-NN.csv")





##################################################

### quantms

## CM (center.median)

#mzTab_pep = read.csv("./PXD000279/onlyPEP-filter-PXD000279.dynamic.sdrf_openms_design_openms.mzTab", sep="\t")
mzTab_pep = read.csv("./PXD000279/PXD000279-for-msqrob2.mzTab", sep="\t")

#mzTab_pep = mzTab_pep[mzTab_pep$opt_global_cv_MS.1002217_decoy_peptide == 0,]
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

#peptide_abundance_study variable[1] 就是 sumIntensity_1
mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                       "sumIntensity_3", "sumIntensity_4", "sumIntensity_5", "sumIntensity_6",
                       "sumIntensity_7", "sumIntensity_8")]
mzTab_pep_sum <- mzTab_pep %>% group_by(sequence, accession) %>% summarise(sumIntensity_1 = sum(sumIntensity_1),
                                                                           sumIntensity_2 = sum(sumIntensity_2),
                                                                           sumIntensity_3 = sum(sumIntensity_3),
                                                                           sumIntensity_4 = sum(sumIntensity_4),
                                                                           sumIntensity_5 = sum(sumIntensity_5),
                                                                           sumIntensity_6 = sum(sumIntensity_6),
                                                                           sumIntensity_7 = sum(sumIntensity_7),
                                                                           sumIntensity_8 = sum(sumIntensity_8))

mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3:6]) > 0)|(rowSums(mzTab_pep_sum[7:10]) > 0)),]

names(mzTab_pep_sum)[names(mzTab_pep_sum) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep_sum))

  pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("UPS1", "UPS1", "UPS1", "UPS1", "UPS2", "UPS2", "UPS2", "UPS2"))

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

L <- makeContrast("conditionUPS2=0", parameterNames = c("conditionUPS2"))
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

volcano <- ggplot(
  rowData(pe[["protein"]])$conditionUPS2,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe[["protein"]])$conditionUPS2, file = "ups-quantms-dep-CM.csv")




## Q (quantile)

mzTab_pep = read.csv("./PXD000279/onlyPEP-filter-PXD000279.dynamic.sdrf_openms_design_openms.mzTab", sep="\t")

#mzTab_pep = mzTab_pep[mzTab_pep$opt_global_cv_MS.1002217_decoy_peptide == 0,]
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

#peptide_abundance_study variable[1] 就是 sumIntensity_1
mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                         "sumIntensity_3", "sumIntensity_4", "sumIntensity_5", "sumIntensity_6",
                         "sumIntensity_7", "sumIntensity_8")]
mzTab_pep_sum <- mzTab_pep %>% group_by(sequence, accession) %>% summarise(sumIntensity_1 = sum(sumIntensity_1),
                                                                           sumIntensity_2 = sum(sumIntensity_2),
                                                                           sumIntensity_3 = sum(sumIntensity_3),
                                                                           sumIntensity_4 = sum(sumIntensity_4),
                                                                           sumIntensity_5 = sum(sumIntensity_5),
                                                                           sumIntensity_6 = sum(sumIntensity_6),
                                                                           sumIntensity_7 = sum(sumIntensity_7),
                                                                           sumIntensity_8 = sum(sumIntensity_8))

mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3:6]) > 0)|(rowSums(mzTab_pep_sum[7:10]) > 0)),]

names(mzTab_pep_sum)[names(mzTab_pep_sum) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep_sum))

pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("UPS1", "UPS1", "UPS1", "UPS1", "UPS2", "UPS2", "UPS2", "UPS2"))

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

L <- makeContrast("conditionUPS2=0", parameterNames = c("conditionUPS2"))
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

volcano <- ggplot(
  rowData(pe[["protein"]])$conditionUPS2,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe[["protein"]])$conditionUPS2, file = "ups-quantms-dep-Q.csv")




## NN (none)

mzTab_pep = read.csv("./PXD000279/onlyPEP-filter-PXD000279.dynamic.sdrf_openms_design_openms.mzTab", sep="\t")

#mzTab_pep = mzTab_pep[mzTab_pep$opt_global_cv_MS.1002217_decoy_peptide == 0,]
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

#peptide_abundance_study variable[1] 就是 sumIntensity_1
mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                         "sumIntensity_3", "sumIntensity_4", "sumIntensity_5", "sumIntensity_6",
                         "sumIntensity_7", "sumIntensity_8")]
mzTab_pep_sum <- mzTab_pep %>% group_by(sequence, accession) %>% summarise(sumIntensity_1 = sum(sumIntensity_1),
                                                                           sumIntensity_2 = sum(sumIntensity_2),
                                                                           sumIntensity_3 = sum(sumIntensity_3),
                                                                           sumIntensity_4 = sum(sumIntensity_4),
                                                                           sumIntensity_5 = sum(sumIntensity_5),
                                                                           sumIntensity_6 = sum(sumIntensity_6),
                                                                           sumIntensity_7 = sum(sumIntensity_7),
                                                                           sumIntensity_8 = sum(sumIntensity_8))

mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3:6]) > 0)|(rowSums(mzTab_pep_sum[7:10]) > 0)),]

names(mzTab_pep_sum)[names(mzTab_pep_sum) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep_sum))

pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("UPS1", "UPS1", "UPS1", "UPS1", "UPS2", "UPS2", "UPS2", "UPS2"))

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

L <- makeContrast("conditionUPS2=0", parameterNames = c("conditionUPS2"))
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

volcano <- ggplot(
  rowData(pe[["protein"]])$conditionUPS2,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() +
  ggtitle("Default workflow")

volcano

write.csv(rowData(pe[["protein"]])$conditionUPS2, file = "ups-quantms-dep-NN.csv")
```

<br>

## PXD007145 dataset
**Successfully run version**
```r
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
#mzTab_pep = read.csv("./PXD007145/only_pep_filterPXD007145-Th.sdrf_openms_design_openms.mzTab", sep="\t")
mzTab_pep = read.csv("./PXD007145/PXD007145-for-msqrob2-3.mzTab", sep="\t")
#mzTab_pep = mzTab_pep[mzTab_pep$opt_global_cv_MS.1002217_decoy_peptide == 0,]
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                           "sumIntensity_3", "sumIntensity_4",
                           "sumIntensity_5", "sumIntensity_6", 
                           "sumIntensity_7", "sumIntensity_8",
                           "sumIntensity_9", "sumIntensity_10",
                           "sumIntensity_11", "sumIntensity_12",
                           "sumIntensity_13", "sumIntensity_14",
                           "sumIntensity_15", "sumIntensity_16",
                           "sumIntensity_17", "sumIntensity_18")]
mzTab_pep_sum <- mzTab_pep %>% group_by(sequence, accession) %>% summarise(sumIntensity_1 = sum(sumIntensity_1),
                                                                             sumIntensity_2 = sum(sumIntensity_2),
                                                                             sumIntensity_3 = sum(sumIntensity_3),
                                                                             sumIntensity_4 = sum(sumIntensity_4),
                                                                             sumIntensity_5 = sum(sumIntensity_5),
                                                                             sumIntensity_6 = sum(sumIntensity_6),
                                                                             sumIntensity_7 = sum(sumIntensity_7),
                                                                             sumIntensity_8 = sum(sumIntensity_8),
                                                                             sumIntensity_9 = sum(sumIntensity_9),
                                                                             sumIntensity_10 = sum(sumIntensity_10),
                                                                             sumIntensity_11 = sum(sumIntensity_11),
                                                                             sumIntensity_12 = sum(sumIntensity_12),
                                                                             sumIntensity_13 = sum(sumIntensity_13),
                                                                             sumIntensity_14 = sum(sumIntensity_14),
                                                                             sumIntensity_15 = sum(sumIntensity_15),
                                                                             sumIntensity_16 = sum(sumIntensity_16),
                                                                             sumIntensity_17 = sum(sumIntensity_17),
                                                                             sumIntensity_18 = sum(sumIntensity_18))


mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3:8]) > 0)|(rowSums(mzTab_pep_sum[9:14]) > 0)|(rowSums(mzTab_pep_sum[15:20]) > 0)),]


names(mzTab_pep_sum)[names(mzTab_pep_sum) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep_sum))

pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("1", "1", "1", "1", "1", "1",
                                      "4", "4", "4", "4", "4", "4",
                                      "10", "10", "10", "10", "10", "10"))

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
#mzTab_pep = read.csv("./PXD007145/only_pep_filterPXD007145-Th.sdrf_openms_design_openms.mzTab", sep="\t")
mzTab_pep = read.csv("./PXD007145/PXD007145-for-msqrob2-3.mzTab", sep="\t")
#mzTab_pep = mzTab_pep[mzTab_pep$opt_global_cv_MS.1002217_decoy_peptide == 0,]
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                         "sumIntensity_3", "sumIntensity_4",
                         "sumIntensity_5", "sumIntensity_6", 
                         "sumIntensity_7", "sumIntensity_8",
                         "sumIntensity_9", "sumIntensity_10",
                         "sumIntensity_11", "sumIntensity_12",
                         "sumIntensity_13", "sumIntensity_14",
                         "sumIntensity_15", "sumIntensity_16",
                         "sumIntensity_17", "sumIntensity_18")]
mzTab_pep_sum <- mzTab_pep %>% group_by(sequence, accession) %>% summarise(sumIntensity_1 = sum(sumIntensity_1),
                                                                           sumIntensity_2 = sum(sumIntensity_2),
                                                                           sumIntensity_3 = sum(sumIntensity_3),
                                                                           sumIntensity_4 = sum(sumIntensity_4),
                                                                           sumIntensity_5 = sum(sumIntensity_5),
                                                                           sumIntensity_6 = sum(sumIntensity_6),
                                                                           sumIntensity_7 = sum(sumIntensity_7),
                                                                           sumIntensity_8 = sum(sumIntensity_8),
                                                                           sumIntensity_9 = sum(sumIntensity_9),
                                                                           sumIntensity_10 = sum(sumIntensity_10),
                                                                           sumIntensity_11 = sum(sumIntensity_11),
                                                                           sumIntensity_12 = sum(sumIntensity_12),
                                                                           sumIntensity_13 = sum(sumIntensity_13),
                                                                           sumIntensity_14 = sum(sumIntensity_14),
                                                                           sumIntensity_15 = sum(sumIntensity_15),
                                                                           sumIntensity_16 = sum(sumIntensity_16),
                                                                           sumIntensity_17 = sum(sumIntensity_17),
                                                                           sumIntensity_18 = sum(sumIntensity_18))


mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3:8]) > 0)|(rowSums(mzTab_pep_sum[9:14]) > 0)|(rowSums(mzTab_pep_sum[15:20]) > 0)),]


names(mzTab_pep_sum)[names(mzTab_pep_sum) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep_sum))

pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("1", "1", "1", "1", "1", "1",
                                     "4", "4", "4", "4", "4", "4",
                                     "10", "10", "10", "10", "10", "10"))

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
#mzTab_pep = read.csv("./PXD007145/only_pep_filterPXD007145-Th.sdrf_openms_design_openms.mzTab", sep="\t")
mzTab_pep = read.csv("./PXD007145/PXD007145-for-msqrob2-3.mzTab", sep="\t")
#mzTab_pep = mzTab_pep[mzTab_pep$opt_global_cv_MS.1002217_decoy_peptide == 0,]
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                         "sumIntensity_3", "sumIntensity_4",
                         "sumIntensity_5", "sumIntensity_6", 
                         "sumIntensity_7", "sumIntensity_8",
                         "sumIntensity_9", "sumIntensity_10",
                         "sumIntensity_11", "sumIntensity_12",
                         "sumIntensity_13", "sumIntensity_14",
                         "sumIntensity_15", "sumIntensity_16",
                         "sumIntensity_17", "sumIntensity_18")]
mzTab_pep_sum <- mzTab_pep %>% group_by(sequence, accession) %>% summarise(sumIntensity_1 = sum(sumIntensity_1),
                                                                           sumIntensity_2 = sum(sumIntensity_2),
                                                                           sumIntensity_3 = sum(sumIntensity_3),
                                                                           sumIntensity_4 = sum(sumIntensity_4),
                                                                           sumIntensity_5 = sum(sumIntensity_5),
                                                                           sumIntensity_6 = sum(sumIntensity_6),
                                                                           sumIntensity_7 = sum(sumIntensity_7),
                                                                           sumIntensity_8 = sum(sumIntensity_8),
                                                                           sumIntensity_9 = sum(sumIntensity_9),
                                                                           sumIntensity_10 = sum(sumIntensity_10),
                                                                           sumIntensity_11 = sum(sumIntensity_11),
                                                                           sumIntensity_12 = sum(sumIntensity_12),
                                                                           sumIntensity_13 = sum(sumIntensity_13),
                                                                           sumIntensity_14 = sum(sumIntensity_14),
                                                                           sumIntensity_15 = sum(sumIntensity_15),
                                                                           sumIntensity_16 = sum(sumIntensity_16),
                                                                           sumIntensity_17 = sum(sumIntensity_17),
                                                                           sumIntensity_18 = sum(sumIntensity_18))


mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3:8]) > 0)|(rowSums(mzTab_pep_sum[9:14]) > 0)|(rowSums(mzTab_pep_sum[15:20]) > 0)),]


names(mzTab_pep_sum)[names(mzTab_pep_sum) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep_sum))

pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("1", "1", "1", "1", "1", "1",
                                     "4", "4", "4", "4", "4", "4",
                                     "10", "10", "10", "10", "10", "10"))

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
```

**Failed to run version**
```r
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

# MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
#   xlab("Peptide index (ordered by data completeness)")

pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
# limma::plotDensities(assay(pe[["peptideLog"]]))

Protein_filter <- rowData(pe[["peptideLog"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
pe <- pe[Protein_filter,,]


pe <- filterFeatures(pe, ~ nNonZero >= 2)


pe <- normalize(pe,
                i = "peptideLog",
                name = "peptideNorm",
                method = "center.median"
)


# Summarization to protein level
pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")

# plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))



# The following code will report an error.
pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)
```

<br>
<hr>

## The reasons for not using msqrob2 in the main text

1. When using the R scripts to run PXD007145 dataset with the .mzTab input type, we encountered an error.<br>
After trying to debug the source code, we found the reason: the number of Intensity columns is equal to the number of condition columns, resulting in the reduction of the two to 0, thus causing an error for the following an if() statement. <br>
About this problem could be found in this link (https://github.com/statOmics/msqrob2/issues/49). And this link (https://github.com/statOmics/MSqRob/issues/92) also shows that MSqRob does not support .mztab input type well, and msqrob2 reuses many MSqRob codes.
2. Please refer to the Supplementary Note 2.
