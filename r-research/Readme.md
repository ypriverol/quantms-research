
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
 

## Citation: 

Bai M, Deng J, Dai C, Pfeuffer J, Perez-Riverol Y. LFQ-based peptide and protein intensity downstream analysis. ChemRxiv. Cambridge: Cambridge Open Engage; 2022; This content is a preprint and has not been peer-reviewed. [preprint manuscript](https://chemrxiv.org/engage/chemrxiv/article-details/6337378ffee74e5821507b75)
