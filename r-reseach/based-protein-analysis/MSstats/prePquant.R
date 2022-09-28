#library('MSstats', warn.conflicts = F, quietly = T, verbose = F)
library(MSstats)
library(reticulate)


# 把maxquant数据格式转为out_msstats
setwd('D:/dataset/R downstream analysis/0-paper/data_benchmark/PXD000279(proteome benchmark+dynamic range dataset)/maxquant output/1-proteomic')
#evidence
mq_ev = data.table::fread("evidence.txt")
#proteingroups
#mq_pg = data.table::fread("proteinGroups.txt")
mq_pg = data.table::fread("原-proteinGroups.txt")
#annotation
annot = data.table::fread("2-experimentalDesign.txt")
maxq_imported = MaxQtoMSstatsFormat(mq_ev, annot, mq_pg, use_log_file = FALSE)
#maxq_imported = MaxQtoMSstatsFormat(mq_ev, annot, mq_pg, use_log_file = FALSE, cartesian=TRUE)
fileData <- maxq_imported
# ---------------------

# If run dataProcess() occuring an error message, please change "summaryMethod = 'TMP'" to "summaryMethod = 'linear'"
DDA2009.proposed <- MSstats::dataProcess(raw = fileData,
                                         normalization = 'equalizeMedians',
                                         summaryMethod = 'TMP',
                                         #summaryMethod = 'linear',  
                                         censoredInt = "NA",
                                         MBimpute = TRUE)

# Automatically create the manually created matrix in MSstats, user manual p23
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

#----------End of creation-----------
#names(ourMatrix) <- ourCondition
#ourMatrix.columns <- ourCondition
#colnames(ourMatrix) <- ourCondition


ourMatrix[1] <- 1
ourMatrix[1] <- -1
row.names(ourMatrix) <- "UPS1-UPS2"
ourCondition[1] = "UPS2"
ourCondition[2] = "UPS1"
colnames(ourMatrix) <- ourCondition



DDA2009.comparisons <- MSstats::groupComparison(contrast.matrix = ourMatrix,
                                                data = DDA2009.proposed)


write.csv(DDA2009.comparisons$ComparisonResult, file="MSstats_output.csv")
