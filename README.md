### GSEA-based Programs for NGS data
#### This repository contains Python programs and example input file(s) from publicly available NGS datasets acquired from Gene Expression Omnibus to reproduce the Gene Set Enrichment Analysis (GSEA)-based approach for gene identification, verification, and comparison presented in the Park & Harris papers (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8342995/ and https://www.cscjournals.org/library/manuscriptinfo.php?mc=IJBB-262). 

### Gene Identification:
#### GSEA1datasetGeneID.py inputs one dataset containing clean and complete NGS data in tab-delimited .txt file format
#### GSEA2datasetsGeneID.py inputs two independent datasets containing clean and complete NGS data both in tab-delimited .txt file format

#### Both GSEAGeneID programs perform the following:
##### 1) Z-score normalization of data across all samples for each gene in a dataset
##### 2) Generation of gene signatures via T-score - returns two gene lists, one per dataset, with associated T-scores and p-values
##### 3) Generation of query gene sets - returns a .gmt file (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) that contains 4 query sets representing the most over- and under-expressed genes from each gene signature (default is 500 genes per query set)
##### 4) GSEA between gene signatures and query gene sets - returns enrichment data (normalized scores, nominal p-values, plots) and leading-edge gene identification

#### NOTE: These algorithms compute T-score per scipy.stats.ttest_ind (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html). This algorithm does not adjust T-score values using a minimum value for σ of .2 * absolute(μ), where μ=0 is adjusted to μ=1, like performed in the GSEA desktop version (https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html).

### Gene Verification/Comparison:
#### GSEAItemCompare.zip inputs your identified gene panels and one dataset containing clean and complete NGS data in tab-delimited .txt file format

#### GSEAItemCompare performs the following:
##### 1) Z-score normalization of data across all samples for each gene in the dataset
##### 2) Generation of gene signature via T-score - returns gene list with associated T-scores and p-values
##### 3) GSEA between gene signature and identified gene panels - returns enrichment data (normalized scores, nominal p-values, plots) and leading-edge gene identification
##### 4) Generation of random query gene sets - returns a .gmt file (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) that contains specified number (1000 default - reduce number if memory issues occur) of query sets representing randomly selected genes from the gene signature
##### 5) GSEA between gene signature and randomly generated gene panels - returns enrichment data (normalized scores and nominal p-values only)
##### 6) Box and whiskars plot of randomly generated normalized enrichment scores compared to achieved scores from identified gene panels

#### Required dependencies: numpy, pandas, scipy, and gseapy (https://gseapy.readthedocs.io/)

#### Limitations of these programs include:
##### 1) Manual preparation of input files prior to program use to ensure clean and complete data that consistently uses the same gene identifier (_e.g._, probes, gene symbols, etc.) across all datasets
##### 2) No GUI = manual manipulation of file names, column indicies, and query set sizes in the code as designated by the "#change ... as needed" comments
#### Please note that Python uses 0 for the index of the first column in a dataframe. Adjust your column counts accordingly.
##### 3) Manual analysis of identified leading-edge genes and preparation of gene panels file, similar to .gmt format, for verification and comparison stages
##### 4) Manual analysis to verify individual panel genes from identified leading-edge genes

### Coming soon!
#### GSEAPathSigGen.py, which inputs a dataset and uses it to define pathway signatures for use in identifying and verifying pathways of interest
