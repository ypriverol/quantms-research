# Validation experiments for quantms

This repository contains a set of JupyterNotebooks and code with analysis and plots for the pipeline https://github.com/bigbio/quantms. All experiments were performed with quantms, and the pipeline should reproduce all analysis results from our manuscript. Here we compared the quantms results with MaxQuant/original results based on same parameters. This repo was organized with data, dataset description by types TMT, LFQ-DDA and LFQ-DIA

## Datasets for quantms experiments

All datasets are directly from quantms and original result/MaxQuant based on same parameters. Very few files are not uploaded due to exceeding github size limit (Please contact author to if they are needed).

### LFQ-DDA

- [PXD000279](https://www.ebi.ac.uk/pride/archive/projects/PXD000279): is an E. coli spiked-in label-free dataset. E. coli proteome was prepared into two conditions with 1:3 ratio (10 ug versus 30 ug) in the triplicates and added into equal amounts of human proteome background (50 ug Hela cells protein extract) and focusing into 24 fractions. The quantms results is [here](./datasets/LFQ-DDA/PXD000279/)

- [PXD014415](https://www.ebi.ac.uk/pride/archive/projects/PXD014415): is a mixture sample of Human cell lysate and yeast cell lysate prepared for label-free LC-MS/MS analysis.The quantms results is [here](./datasets/LFQ-DDA/PXD014415/)

- [PXD002099](https://www.ebi.ac.uk/pride/archive/projects/PXD002099): consists of 48 UPS1 proteins spiked into a yeast proteome digest in five different concentrations: 2, 4, 10, 25 and 50 fmol/ul. The quantms results is [here](./datasets/LFQ-DDA/PXD002099/)

- [CPTAC S06](https://cptac-data-portal.georgetown.edu/cptac/dataPublic/list/LTQ-Orbitrap%4086?currentPath=%2FPhase_I_Data%2FStudy6): The CPTAC Study 6 dataset includes 48 UPS1 proteins-spiked into a yeast proteome digest in five different concentrations: 0.25, 0.74, 2.2, 6.7 and 20 fmol/ul. The quantms results is [here](./datasets/LFQ-DDA/CPTACS06/)

- [PXD001819](https://www.ebi.ac.uk/pride/archive/projects/PXD001819): contains 48 UPS1 proteins spiked into yeast proteome digest in nine different concentrations: 0.05, 0.125, 0.25, 0.5, 2.5, 5, 12.5, 25 and 50 fmol/ul. The quantms results is [here](./datasets/LFQ-DDA/PXD001819/)

- [PASS00589](ftp://PASS00589:WF6554orn@ftp.peptideatlas.org/): A shotgun standard data set consists of 12 nonhuman proteins spiked into a constant human background (HEK-293), including 8 different sample groups with known concentrations of spike-in proteins in three master mixes. The quantms results is [here](./datasets/LFQ-DDA/PASS00589/)

- [PXD020248](https://www.ebi.ac.uk/pride/archive/projects/PXD020248): HepG2 cell lines are treated with DMSO and BaP separately, and LFQ and TMT quantfication method were used. Six biological replicates were prepared. The quantms results is [here](./datasets/LFQ-DDA/PXD020248/)

### TMT

- [PXD005486](https://www.ebi.ac.uk/pride/archive/projects/PXD005486): is E. Coli background with 12 spike-in human proteins and a bovin protein using various known concentrations. The analytical method employed tandem mass tags (TMT), proteins were spiked twice using the same concentration in different channels and only once for the two highest concentrations. A total of 12 peptide fractions were prepared. The quantms results is [here](./datasets/TMT/PXD005486/)

- [PXD013277](https://repository.jpostdb.org/entry/JPST000562): contains samples with different spiked in amounts of E.Coli protein extract (3 biological replicates) in MCF-7 background. Each sample were labeled with an isobaric TMT-tag (TMT10-plex). All search settings are the same as original paper. The quantms results is [here](./datasets/TMT/PXD013277/)

- [PXD014414](https://www.ebi.ac.uk/pride/archive/projects/PXD014414): Using multiplex quantitative tandem mass tag-based proteomics to quantify proteins in MBC, TNBC, and normal breast from 27 patients. The quantms results is [here](./datasets/TMT/PXD014414/)

- [PXD002875](https://www.ebi.ac.uk/pride/archive/projects/PXD002875): TMT9-plex analysis of S. cerevisiae grown on three carbon sources. The quantms results is [here](./datasets/TMT/PXD002875/)

- [PXD020248](https://www.ebi.ac.uk/pride/archive/projects/PXD020248): see LFQ part

### LFQ-DIA

- [PXD026600](http://proteomecentral.proteomexchange.org/dataset/PXD026600): This large dataset includes 96 DIA .raw files acquired from a complex proteomic standard composed of an E.coli protein background spiked-in with 8 different concentrations of 48 human proteins (UPS1 Sigma). These 8 samples were analyzed in triplicates on an Orbitrap mass spectrometer with 4 different DIA window schemes.

## Description of notebooks

Original results and figures which were included in the manuscript are available in the following notebook files. These can be directly viewed in Github or downloaded and re-run with the re-created analysis on your own machine.
All analysis results are directly from `datasets/` folder.

### LFQ-DDA

- [PXD000279 Benchmark](./notebooks/LFQ-DDA/PXD000279Benchmark.ipynb)

- [PXD014415 Benchmark](./notebooks/LFQ-DDA/PXD014415Benchmark.ipynb)

- [PXD002099 Benchmark](./notebooks/LFQ-DDA//PXD002099Benchmark.ipynb)
- [CPTACS06 Benchmark](./notebooks/LFQ-DDA/CPTACS06Benchmark.ipynb)

- [PXD001819 Benchmark](./notebooks/LFQ-DDA/PXD001819Benchmark.ipynb)

- [PASS00589 Benchmark](./notebooks/LFQ-DDA/PASS00589.ipynb)

- [PXD020248 Benchmark](./notebooks/LFQ-DDA/PXD020248benchmarking.ipynb)

### TMT

- [PXD005486 Benchmark](./notebooks/TMT/PXD005486Benchmark.ipynb)

- [PXD013277 Benchmark](./notebooks/TMT/PXD013277Benchmark.ipynb)

- [PXD014414 Benchmark](./notebooks/TMT/PXD014414Benchmark.ipynb)

- [PXD02875 Benchmark](./notebooks/TMT/PXD002875Benchmark.ipynb)

### LFQ-DIA

- [PXD026600 Benchmark](./notebooks/LFQ-DIA/PXD026600Benchmark.ipynb)