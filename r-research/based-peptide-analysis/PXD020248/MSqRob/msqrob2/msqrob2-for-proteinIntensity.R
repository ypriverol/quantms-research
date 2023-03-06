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


### quantms

## div.median
mzTab_pep = read.csv("./PXD020248/onlyPEP-filter-PXD020248.mzTab", sep="\t")
mzTab_pep = mzTab_pep[!grepl("DECOY", mzTab_pep$accession),]

mzTab_pep = mzTab_pep[,c("sequence", "accession", "sumIntensity_1", "sumIntensity_2",
                         "sumIntensity_3", "sumIntensity_4",
                         "sumIntensity_5", "sumIntensity_6", 
                         "sumIntensity_7", "sumIntensity_8",
                         "sumIntensity_9", "sumIntensity_10",
                         "sumIntensity_11", "sumIntensity_12")]
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
                                                                           sumIntensity_12 = sum(sumIntensity_12))


mzTab_pep_sum <- mzTab_pep_sum[which((rowSums(mzTab_pep_sum[,3:8]) > 0)|(rowSums(mzTab_pep_sum[9:14]) > 0)),]

names(mzTab_pep_sum)[names(mzTab_pep_sum) =="accession"] <-"Proteins"

ecols <- grep("sumIntensity_", names(mzTab_pep_sum))

pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)


colData(pe)$condition <- as.factor(c("BAP", "BAP", "BAP", "BAP", "BAP", "BAP",
                                     "DMSO", "DMSO", "DMSO", "DMSO", "DMSO", "DMSO"))

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA


Protein_filter <- rowData(pe[["peptideRaw"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideRaw"]])$Proteins)
pe <- pe[Protein_filter,,]

pe <- filterFeatures(pe, ~ nNonZero >= 2)

pe <- normalize(pe,
                i = "peptideRaw",
                name = "peptideNorm",
                method = "div.median"
)

pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])


write.csv(assay(pe[["protein"]]), file = "msqobr2-protein-intensity-divM.csv")




## diff.median
pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)

colData(pe)$condition <- as.factor(c("BAP", "BAP", "BAP", "BAP", "BAP", "BAP",
                                     "DMSO", "DMSO", "DMSO", "DMSO", "DMSO", "DMSO"))

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA


Protein_filter <- rowData(pe[["peptideRaw"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideRaw"]])$Proteins)
pe <- pe[Protein_filter,,]

pe <- filterFeatures(pe, ~ nNonZero >= 2)

pe <- normalize(pe,
                i = "peptideRaw",
                name = "peptideNorm",
                method = "diff.median"
)

pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])


write.csv(assay(pe[["protein"]]), file = "msqobr2-protein-intensity-diffM.csv")



## CM
pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)
colData(pe)$condition <- as.factor(c("BAP", "BAP", "BAP", "BAP", "BAP", "BAP",
                                     "DMSO", "DMSO", "DMSO", "DMSO", "DMSO", "DMSO"))

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA


Protein_filter <- rowData(pe[["peptideRaw"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideRaw"]])$Proteins)
pe <- pe[Protein_filter,,]

pe <- filterFeatures(pe, ~ nNonZero >= 2)

pe <- normalize(pe,
                i = "peptideRaw",
                name = "peptideNorm",
                method = "center.median"
)

pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])


write.csv(assay(pe[["protein"]]), file = "msqobr2-protein-intensity-CM.csv")





## Q
pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)

colData(pe)$condition <- as.factor(c("BAP", "BAP", "BAP", "BAP", "BAP", "BAP",
                                     "DMSO", "DMSO", "DMSO", "DMSO", "DMSO", "DMSO"))

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA


Protein_filter <- rowData(pe[["peptideRaw"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideRaw"]])$Proteins)
pe <- pe[Protein_filter,,]

pe <- filterFeatures(pe, ~ nNonZero >= 2)

pe <- normalize(pe,
                i = "peptideRaw",
                name = "peptideNorm",
                method = "quantiles"
)

pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])


write.csv(assay(pe[["protein"]]), file = "msqobr2-protein-intensity-Q.csv")




## NN
pe <- readQFeatures(
  table = mzTab_pep_sum, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)

colData(pe)$condition <- as.factor(c("BAP", "BAP", "BAP", "BAP", "BAP", "BAP",
                                     "DMSO", "DMSO", "DMSO", "DMSO", "DMSO", "DMSO"))

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA


Protein_filter <- rowData(pe[["peptideRaw"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideRaw"]])$Proteins)
pe <- pe[Protein_filter,,]

pe <- filterFeatures(pe, ~ nNonZero >= 2)

pe <- addAssay(pe,pe[["peptideRaw"]],"peptideNorm")

pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")

pe <- msqrob2::msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])


write.csv(assay(pe[["protein"]]), file = "msqobr2-protein-intensity-NN.csv")