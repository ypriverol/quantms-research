
# Differential expression analysis using R packages

This folder contains the data for the evaluation and bernchmark of R-pakages for differential expresion analyisis. In summary, we used 3 datasets and two softwares [MaxQuant](https://www.maxquant.org/) and [quantms](https://github.com/bigbio/quantms) to evaluate of 8 different R-packages when performing differential expression analysis. 

### Datasets

The following datasets were used to perform the analysis:

- [PXD000279 - UPS spiked dataset](https://www.ebi.ac.uk/pride/archive/projects/PXD000279): Consists of two E. coli digested samples (with 4 replicates for each sample); each half of the samples are enriched with one of two “Universal Protein Standards” (UPS1 and UPS2). Both samples contained the same 48 recombinant human proteins, which were either mixed in equal amounts (UPS1) or spanned multiple orders of magnitude at a determined ratio (UPS2). Based on the experimental design, 40 of the 48 UPS proteins and none of the E. coli proteins should be detected as differentially expressed.

- [PXD007145 - large-scale mix dataset](https://www.ebi.ac.uk/pride/archive/projects/PXD007145): The large-scale mix dataset contains multiple species mixture, in which Yeast proteome was diluted into fixed ratios of 1:4:10 and added to a background of 1:1:1 human proteome to simulate the real experimental data. Six technical replicates were used for each sample to measure the coefficient of variations.

- [PXD020248 - toxicology dataset](https://www.ebi.ac.uk/pride/archive/projects/PXD020248): The toxicology dataset is a cell line hepatocytes sample (HepG2) treated with benzo[a]pyrene (BaP) using a concentration of 2 μM. The original work (31) benchmarked TMT and LFQ analytical methods using the same sample.

### Quantification softwares:

- **MaxQuant**: To evaluate each tool's parameter combinations, algorithms for data processing and protein quantification analysis, we analyzed the datasets with MaxQuant. Raw data were processed with MaxQuant (version v1.6.10.43) before the DE analysis with each tool. We used default parameters except that “the min ratio of LFQ” was set as 1 and “matching between runs” was enabled. 

- **quantms**: quantms (and its predecessor proteomicsLFQ) is a cloud-based workflow that uses OpenMS tools and DIA-NN to enable quantitative analysis of LFQ data-dependent (LFQ-DDA) and independent acquisition (LFQ-DIA) and TMT data (https://quantms.readthedocs.io/en/latest/). In the present work, we used the LFQ-DDA sub-workflow of the pipeline on the three datasets to do the peptide quantification benchmarking. The sub-workflow performs peptide identifications using Comet and MSGF+ and feature detection using proteomicsLFQ in OpenMS.
 
## R-packages evaludated

- [MSstats](https://github.com/Vitek-Lab/MSstats) is an open-source R-package for peptide and protein quantification in mass spectrometry-based proteomics. MSstats supports multiple data acquisition types: data-dependent acquisition (DDA) -both LFQ and label-based workflows-, data-independent acquisition (DIA) and targeted approaches. 

- [Proteus and limma](https://github.com/bartongroup/Proteus) is used for DE analysis of MaxQuant output data, and differential expression analysis based on the popular algorithm/package limma (37). Proteus supports two normalization methods: equalize median and quantile, and it uses a mean-variance relationship to estimate variance (limma) where data are missing. The Proteus Shiny application allows users to perform the analysis with one click if the data is provided in MaxQuant format. 

- [prolfqua](https://github.com/fgcz/prolfqua) integrates the basic steps of differential expression analysis workflow: quality control, data normalization, protein aggregation, statistical modelling, hypothesis testing, and sample size estimation. The modular design of prolfqua enables users to select the optimal differential expression analysis algorithm. 

- [ProVision](https://github.com/JamesGallant/ProVision) is an R-shiny web application to facilitate the analysis of LFQ and TMT proteomics experiments. ProVision is designed for end-users (e.g., biologists), with a set of graphical interfaces to guide the users through data processing, parameter selection, and result presentation.

- [LFQ-Analyst](https://github.com/MonashBioinformaticsPlatform/LFQ-Analyst) is an interactive, R-Shiny-based platform for quickly and easily analyzing and visualizing unlabeled proteomics data preprocessed with MaxQuant. LFQ-Analyst can process LFQ intensity, and its quality control report contains multiple visualization plots (volcano plots, heatmaps, and box plots) of differentially expressed. 

- [Eatomics](https://github.com/Millchmaedchen/Eatomics) is also an R-shiny application for the interactive exploration of quantitative proteomics data from MaxQuant, integrating quality control, differential abundance analysis, and enrichment analysis. It has a variety of interactive exploration possibilities and a unique experimental design setup module that interactively transforms a given research hypothesis into a differential abundance and enrichment analysis formula. 

- [DAPAR and ProStaR](http://www.prostar-proteomics.org/) are two tools dedicated to the discovery of differential analysis of quantitative data generated by proteomic experiments. DAPAR is an R-package which provides five processing steps (filtering, normalization, imputation, aggregation, and difference analysis), based on those functions, ProStaR provides an R-Shiny web platform for interactive exploring.

- [MSqRob](https://github.com/statOmics/MSqRob) s a free and open-source R package that can handle virtually any experimental proteomics design. MSqRob supports multiple types of inputs, including MaxQuant, moFF mzTab and Progenesis. 

- [Perseus](https://maxquant.net/perseus/) is designed for DE analysis of quantitative results in the MaxQuant ecosystem. Perseus is a desktop application which offers a wide variety of algorithms for MaxQuant data normalization, imputation, batch correction and differential expression analysis.  Users need to manually annotate the condition during data processing and can choose different types of intensity: raw intensities, LFQ intensity or IBAQ values. 

## Quantification Results

### [PXD000279](https://github.com/ypriverol/quantms-research/tree/main/r-research/based-peptide-analysis/PXD000279)

- MaxQuant output: [peptides.txt](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/peptides.txt), [proteinGroups.txt](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/proteinGroups.txt)

- quantms output: [mzTab](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/onlyPEP-filter-PXD000279.dynamic.sdrf_openms_design_openms.mzTab), [msstats input](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/out_msstats.csv)

### [PXD007145](https://github.com/ypriverol/quantms-research/tree/main/r-research/based-peptide-analysis/PXD007145)

- MaxQuant output: [peptides.txt](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145/peptides.txt), [proteinGroups.txt](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145/proteinGroups.txt)

- quantms output: [mzTab](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145/only_pep_filterPXD007145-Th.sdrf_openms_design_openms.mzTab), [msstats input](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145-Th.sdrf_openms_design_msstats_in.csv)

### [PXD020248](https://github.com/ypriverol/quantms-research/tree/main/r-research/based-peptide-analysis/PXD020248)

- MaxQuant output: [peptides.txt](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD020248/peptides.txt), [proteinGroups.txt](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD020248/proteinGroups.txt)

- quantms output: [mzTab](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD020248/onlyPEP-filter-PXD020248.mzTab), [msstats input](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/PXD020248/out_msstats.csv)

## Scripts and Protocol used by each R-package

### MSstats

The main analysis scripts include PXD000279 [peptide-level](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/MSstats/prePquant.R), PXD007145 [peptide-level](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145/MSstats/MSstats_filter.ipynb) and PXD020248 [peptide-level](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-peptide-analysis/PXD020248/MSstats/prePquant.R). The msstats_in file generated by quantms was input. Only unique peptides and proteins with at least 2 unique peptides were considered. Three normalization methods and two imputation strategies were evaluated. Tukey’s median polish summarymethod was used, which is a robust estimation method. The protein intensity can be exported.

### Proteus and limma

The main analysis scripts include PXD000279 [peptide-level](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-peptide-analysis/PXD000279/Proteus/proteus.R)/[protein-level](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-protein-analysis/PXD000279/proteus/proteus.R), PXD007145 [peptide-level](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-peptide-analysis/PXD007145/Proteus/proteus_filter.ipynb)/[protein-level](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-protein-analysis/PXD007145/Proteus/proteus.r) and PXD020248 [peptide-level](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-peptide-analysis/PXD020248/Proteus/proteus.R). The evidence.txt generated by MaxQuant or msstats_in file generated by quantms were as input. The quantms need pre-processing. Only unique peptides and proteins with at least 2 unique peptides were considered.

### prolfqua

The main analysis scripts include PXD000279 [protein-level](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-protein-analysis/PXD000279/prolfqua/prolfqua.R) and PXD007145 [protein-level](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-protein-analysis/PXD007145/prolfqua/prolfqua.r).

### ProVision

ProVision provided an online analysis [service](https://provision.shinyapps.io/provision/). We select LFQ intensity column to quantify. Default parameters were selected and then export results.

### LFQ-Analyst

LFQ-Analyst provided an online analysis [service](https://bioinformatics.erc.monash.edu/apps/LFQ-Analyst/). We uploaded proteinGroups.txt generated by MaxQuant and experimental design table. We select LFQ intensity column to quantify. Different imputation methods were tested, and then Benjamini Hochberg FDR correction was used.

### Eatomics

[Eatomics](https://github.com/Millchmaedchen/Eatomics) provided web application developed using Shiny. We performed Eatomics by shiny. We uploaded the proteinGroups.txt and experimental design table eg [PXD000279](https://github.com/daichengxin/quantms-research/blob/main/r-research/based-protein-analysis/PXD000279/eatomics/clinicaldata-1-experimentalDesign.txt). The LFQ quantification type was selected. Different imputation methods were tested as reported paper.

### DAPAR and ProStaR

The DAPAR2(1.28.0) R package was installed, and then performed shiny application. Different normalization and imputation methods
combination were tested.

### MSqRob

MSqRob 0.7.7 provided a shiny app (graphical user interface). The peptides.txt generated by MaxQuant and experimental design table were uploaded. For quantms, the peptide part of mzTab file was used as input data. Different normalization methods were evaluated and the result was export. Only unique peptides and proteins with at least 2 unique peptides were considered. In addition, we also tried to test msqrob2. But msqrob2 is incompatible with mzTab file from quantms. https://github.com/statOmics/MSqRob/issues/92 also hints that MSqRob doesn’t support. mztab input files well from other search engines, and msqrob2 reuses MSqRob codes. There are duplicate rows of peptide sequences due to different modification, charge and retention time combinations. But the peptide sequences are unique from MaxQuant peptide output, which is tested  in previous paper including MSqRob paper. So the output information is incompatible between them. Even if we can do some manual conversions (from msstats_in file to msqrob2 input) beforehand, different data conversion processes can potentially affect the results. We have tried different conversion steps, the first without any aggregation of abundance, the second with aggregation of a sequence with different charge only and the third with aggregation of both different charges and modifications. None of which have produced accepted results. In summary, MSqRob/msqrob2 may does not consider this type of data (incompatible). Also we can see from the results. The Scripts and Protocol of msqrob2 can be found in [here](https://github.com/ypriverol/quantms-research/blob/main/r-research/based-peptide-analysis/the%20benchmark%20for%20both%20versions%20of%20msqrob.md).

### Perseus

Perseus 1.6.13.0 software was installed, which is desktop application. LFQ Intensity was selected as quantification type. The contaminant and decoy were filtered. Then log2 transform was performed. Different imputation methods were tested, and then two -sample test and Benjamini Hochberg FDR correction were performed.

## Citation:

Bai M, Deng J, Dai C, Pfeuffer J, Perez-Riverol Y. LFQ-based peptide and protein intensity downstream analysis. ChemRxiv. Cambridge: Cambridge Open Engage; 2022; This content is a preprint and has not been peer-reviewed. [preprint manuscript](https://chemrxiv.org/engage/chemrxiv/article-details/6337378ffee74e5821507b75)
