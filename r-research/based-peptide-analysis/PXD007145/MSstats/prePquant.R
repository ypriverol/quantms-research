# peptide analyasis
library(MSstats)
library(reticulate)

setwd('D:/dataset/R downstream analysis/0-paper/data_benchmark/0-reviewer/peptide based/PXD007145')

### quantms
raw.om <- read.csv('./out_msstats_filter.csv')
raw.om <- raw.om[!grepl("DECOY", raw.om$ProteinName),]
fileData <- OpenMStoMSstatsFormat(raw.om)

# equalizeMedians (EM)
DDA2009.proposed <- MSstats::dataProcess(raw = fileData,
                                         normalization = 'equalizeMedians',
                                         summaryMethod = 'TMP',
                                         censoredInt = "NA",
                                         MBimpute = TRUE)

write.csv(DDA2009.proposed$ProteinLevelData, file="MSstats_output-quantms-EM-proteinIntensity.csv")

len <- length(levels(DDA2009.proposed$FeatureLevelData$GROUP))

tmp <- t(combn(len,2))
matrix_len = length(t(combn(len,2))) / 2

ourMatrix <- matrix(c(0:0),nrow=matrix_len,ncol=len)

for(i in 1:matrix_len){
  ourMatrix[i, tmp[i]] = -1
  ourMatrix[i, tmp[i + matrix_len]] = 1
}

ourCondition <- levels(DDA2009.proposed$ProteinLevelData$GROUP)
tmp_name <- matrix(ourCondition, nr=len, nc=1)
name <- matrix(nr=matrix_len, nc=1)
for(i in 1:matrix_len){
  name[i,1] <- sprintf('%s-%s', tmp_name[tmp[i+matrix_len]], tmp_name[tmp[i]])
}
row.names(ourMatrix) <- name

colnames(ourMatrix) <- ourCondition

DDA2009.comparisons <- MSstats::groupComparison(contrast.matrix = ourMatrix,
                                                data = DDA2009.proposed)

write.csv(DDA2009.comparisons$ComparisonResult, file="MSstats_output-quantms-EM.csv")



# Quantile (Q)
DDA2009.proposed <- MSstats::dataProcess(raw = fileData,
                                         normalization = 'Quantile',
                                         summaryMethod = 'TMP',
                                         censoredInt = "NA",
                                         MBimpute = TRUE)

write.csv(DDA2009.proposed$ProteinLevelData, file="MSstats_output-quantms-Q-proteinIntensity.csv")

len <- length(levels(DDA2009.proposed$FeatureLevelData$GROUP))

tmp <- t(combn(len,2))
matrix_len = length(t(combn(len,2))) / 2

ourMatrix <- matrix(c(0:0),nrow=matrix_len,ncol=len)

for(i in 1:matrix_len){
  ourMatrix[i, tmp[i]] = -1
  ourMatrix[i, tmp[i + matrix_len]] = 1
}

ourCondition <- levels(DDA2009.proposed$ProteinLevelData$GROUP)
tmp_name <- matrix(ourCondition, nr=len, nc=1)
name <- matrix(nr=matrix_len, nc=1)
for(i in 1:matrix_len){
  name[i,1] <- sprintf('%s-%s', tmp_name[tmp[i+matrix_len]], tmp_name[tmp[i]])
}
row.names(ourMatrix) <- name

colnames(ourMatrix) <- ourCondition

DDA2009.comparisons <- MSstats::groupComparison(contrast.matrix = ourMatrix,
                                                data = DDA2009.proposed)

write.csv(DDA2009.comparisons$ComparisonResult, file="MSstats_output-quantms-Q.csv")



# None (NN)
DDA2009.proposed <- MSstats::dataProcess(raw = fileData,
                                         normalization = FALSE,
                                         summaryMethod = 'TMP',
                                         censoredInt = "NA",
                                         MBimpute = TRUE)

write.csv(DDA2009.proposed$ProteinLevelData, file="MSstats_output-quantms-NN-proteinIntensity.csv")

len <- length(levels(DDA2009.proposed$FeatureLevelData$GROUP))

tmp <- t(combn(len,2))
matrix_len = length(t(combn(len,2))) / 2

ourMatrix <- matrix(c(0:0),nrow=matrix_len,ncol=len)

for(i in 1:matrix_len){
  ourMatrix[i, tmp[i]] = -1
  ourMatrix[i, tmp[i + matrix_len]] = 1
}

ourCondition <- levels(DDA2009.proposed$ProteinLevelData$GROUP)
tmp_name <- matrix(ourCondition, nr=len, nc=1)
name <- matrix(nr=matrix_len, nc=1)
for(i in 1:matrix_len){
  name[i,1] <- sprintf('%s-%s', tmp_name[tmp[i+matrix_len]], tmp_name[tmp[i]])
}
row.names(ourMatrix) <- name

colnames(ourMatrix) <- ourCondition

DDA2009.comparisons <- MSstats::groupComparison(contrast.matrix = ourMatrix,
                                                data = DDA2009.proposed)

write.csv(DDA2009.comparisons$ComparisonResult, file="MSstats_output-quantms-NN.csv")


######################################################


### maxquant
#evidence
mq_ev = data.table::fread("evidence_filter.txt")
#proteingroups
mq_pg = data.table::fread("proteinGroups_filter.txt")
#annotation
annot = data.table::fread("2-experimentalDesign.txt")
maxq_imported = MaxQtoMSstatsFormat(mq_ev, annot, mq_pg, use_log_file = FALSE)
fileData <- maxq_imported


# equalizeMedians (EM)
DDA2009.proposed <- MSstats::dataProcess(raw = fileData,
                                         normalization = 'equalizeMedians',
                                         summaryMethod = 'TMP',
                                         censoredInt = "NA",
                                         MBimpute = TRUE)

len <- length(levels(DDA2009.proposed$FeatureLevelData$GROUP))

tmp <- t(combn(len,2))
matrix_len = length(t(combn(len,2))) / 2

ourMatrix <- matrix(c(0:0),nrow=matrix_len,ncol=len)

for(i in 1:matrix_len){
  ourMatrix[i, tmp[i]] = -1
  ourMatrix[i, tmp[i + matrix_len]] = 1
}

ourCondition <- levels(DDA2009.proposed$ProteinLevelData$GROUP)
tmp_name <- matrix(ourCondition, nr=len, nc=1)
name <- matrix(nr=matrix_len, nc=1)
for(i in 1:matrix_len){
  name[i,1] <- sprintf('%s-%s', tmp_name[tmp[i+matrix_len]], tmp_name[tmp[i]])
}
row.names(ourMatrix) <- name

colnames(ourMatrix) <- ourCondition

DDA2009.comparisons <- MSstats::groupComparison(contrast.matrix = ourMatrix,
                                                data = DDA2009.proposed)

write.csv(DDA2009.comparisons$ComparisonResult, file="MSstats_output-maxquant-EM.csv")



# Quantile (Q)
DDA2009.proposed <- MSstats::dataProcess(raw = fileData,
                                         normalization = 'Quantile',
                                         summaryMethod = 'TMP',
                                         censoredInt = "NA",
                                         MBimpute = TRUE)

len <- length(levels(DDA2009.proposed$FeatureLevelData$GROUP))

tmp <- t(combn(len,2))
matrix_len = length(t(combn(len,2))) / 2

ourMatrix <- matrix(c(0:0),nrow=matrix_len,ncol=len)

for(i in 1:matrix_len){
  ourMatrix[i, tmp[i]] = -1
  ourMatrix[i, tmp[i + matrix_len]] = 1
}

ourCondition <- levels(DDA2009.proposed$ProteinLevelData$GROUP)
tmp_name <- matrix(ourCondition, nr=len, nc=1)
name <- matrix(nr=matrix_len, nc=1)
for(i in 1:matrix_len){
  name[i,1] <- sprintf('%s-%s', tmp_name[tmp[i+matrix_len]], tmp_name[tmp[i]])
}
row.names(ourMatrix) <- name

colnames(ourMatrix) <- ourCondition

DDA2009.comparisons <- MSstats::groupComparison(contrast.matrix = ourMatrix,
                                                data = DDA2009.proposed)

write.csv(DDA2009.comparisons$ComparisonResult, file="MSstats_output-maxquant-Q.csv")



# None (NN)
DDA2009.proposed <- MSstats::dataProcess(raw = fileData,
                                         normalization = FALSE,
                                         summaryMethod = 'TMP',
                                         censoredInt = "NA",
                                         MBimpute = TRUE)

len <- length(levels(DDA2009.proposed$FeatureLevelData$GROUP))

tmp <- t(combn(len,2))
matrix_len = length(t(combn(len,2))) / 2

ourMatrix <- matrix(c(0:0),nrow=matrix_len,ncol=len)

for(i in 1:matrix_len){
  ourMatrix[i, tmp[i]] = -1
  ourMatrix[i, tmp[i + matrix_len]] = 1
}

ourCondition <- levels(DDA2009.proposed$ProteinLevelData$GROUP)
tmp_name <- matrix(ourCondition, nr=len, nc=1)
name <- matrix(nr=matrix_len, nc=1)
for(i in 1:matrix_len){
  name[i,1] <- sprintf('%s-%s', tmp_name[tmp[i+matrix_len]], tmp_name[tmp[i]])
}
row.names(ourMatrix) <- name

colnames(ourMatrix) <- ourCondition

DDA2009.comparisons <- MSstats::groupComparison(contrast.matrix = ourMatrix,
                                                data = DDA2009.proposed)

write.csv(DDA2009.comparisons$ComparisonResult, file="MSstats_output-maxquant-NN.csv")