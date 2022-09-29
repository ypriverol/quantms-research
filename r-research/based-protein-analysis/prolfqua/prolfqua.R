library(prolfqua)
library(tidyverse)

setwd('D:/dataset/R downstream analysis/0-paper/data_benchmark/PXD000279(proteome benchmark+dynamic range dataset)/2-dynamicrangebenchmark')
fileData <- read.delim("proteinGroups.txt", sep = '\t')
annotation <- read.delim("prolfqua-inputAnnotation.txt", sep = '\t')   #根据startdata里的raw.file列来确定annotation的

startdata <- prolfqua::tidyMQ_ProteinGroups(fileData)
annotation <- annotation %>% mutate(raw.file = tolower(raw.file))
startdata <- dplyr::inner_join(annotation, startdata, by = "raw.file")
startdata <- dplyr::filter(startdata, nr.peptides > 1) #remove all proteins identified only by a single peptide.


#Then you need to tell prolfqua which columns in the data frame contain what information. 
#You do it using the AnalysisTableAnnotation class.
atable <- AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("proteinID")
atable$factors[["Condition_"]] = "condition"
atable$factors[["replicate"]] = "replicate"      
atable$factorDepth <- 1
#startdata$Run_ID <- as.integer(startdata$Run_ID)
atable$setWorkIntensity("mq.protein.intensity")


config <- AnalysisConfiguration$new(atable)
adata <- setup_analysis(startdata, config)



#Create the LFQData class instance and remove zeros from data (MaxQuant encodes missing values with zero).
lfqdata <- LFQData$new(adata, config)
lfqdata$remove_small_intensities()
lfqdata$factors()

#You can convert the data into wide format.
lfqdata$to_wide()$data
wr <- lfqdata$get_Writer()
wr$write_wide(".")




## Visualization of not normalized data
lfqplotter <- lfqdata$get_Plotter()
density_nn <- lfqplotter$intensity_distribution_density()

lfqdata$get_Summariser()$plot_missingness_per_group()
lfqplotter$missigness_histogram()

stats <- lfqdata$get_Stats()
stats$violin()

prolfqua::table_facade( stats$stats_quantiles()$wide, paste0("quantile of ",stats$stat ))

stats$density_median()


## Normalize protein intensities
#We normalize the data by log2 transforming and then z−scaling.
lt <- lfqdata$get_Transformer()
transformed <- lt$log2()$robscale()$lfq
transformed$config$table$is_intensity_transformed

pl <- transformed$get_Plotter()
density_norm <- pl$intensity_distribution_density()

gridExtra::grid.arrange(density_nn, density_norm)

pl$pairs_smooth()

p <- pl$heatmap_cor()
p



## Fitting a linear model
transformed$config$table$getWorkIntensity()

formula_Condition <-  strategy_lm("transformedIntensity ~ Condition_")
# specify model definition
modelName  <- "Model"
unique(transformed$data$Condition_)
# [1] "30 ug" "10 ug"

#Contrasts <- c("30 ug - 10 ug" = "Condition_30 ug - Condition_10 ug")
#Contrasts <- c("HvsL" = "Condition_30 ug - Condition_10 ug")
#Contrasts <- c("30 ug - 10 ug" = "Condition_30ug - Condition_10ug")
#Contrasts <- c("HvsL" = "Condition_30ug - Condition_10ug")

Contrasts <- c("ups1 - ups2" = "Condition_ups1 - Condition_ups2")

#Contrasts <- c("LvsH" = "Condition_10ug - Condition_30ug")

mod <- prolfqua::build_model(
  transformed$data,
  formula_Condition,
  subject_Id = transformed$config$table$hierarchyKeys() )

mod$anova_histogram("FDR.Pr..F.")


aovtable <- mod$get_anova()
head(aovtable)
dim(aovtable)
xx <- aovtable |> dplyr::filter(FDR.Pr..F. < 0.2)
dim(xx)

signif <- transformed$get_copy()
signif$data <- signif$data |> dplyr::filter(protein_Id %in% xx$protein_Id)
hmSig <- signif$get_Plotter()$heatmap()
hmSig


## Compute contrasts
contr <- prolfqua::Contrasts$new(mod, Contrasts)
v1 <- contr$get_Plotter()$volcano()

contr <- prolfqua::ContrastsModerated$new(contr)
#contr$get_contrasts_sides()
contrdf <- contr$get_contrasts()
#View(contrdf)

plotter <- contr$get_Plotter()
v2 <- plotter$volcano()
gridExtra::grid.arrange(v1$FDR,v2$FDR, ncol = 1)

plotter$ma_plotly()


#write.csv(v1$FDR$data, "prolfqua-WaldTest.csv")
write.csv(v2$FDR$data, "prolfqua-WaldTest_moderated.csv")