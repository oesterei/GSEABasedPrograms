### GSEA-based Programs for NGS data
#### This repository contains Python programs and example input file(s) from publicly available NGS datasets acquired from Gene Expression Omnibus to reproduce the Gene Set Enrichment Analysis (GSEA)-based approach for gene identification and verification presented in the Park & Harris papers (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8342995/ and https://www.cscjournals.org/library/manuscriptinfo.php?mc=IJBB-262). 

#### Now available!
#### GSEA1datasetGeneIDonly.py inputs one dataset containing clean and complete NGS data in tab-delimited .txt file format then performs the following:
#### 1) Z-score normalization of data across all samples for each gene in the dataset - returns two z-scored datasets
#### 2) Generation of gene signatures via T-score - returns two gene lists, one per dataset, with associated T-scores and p-values
#### 3) Generation of query gene sets - returns a .txt file similar to a .gmt (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) that contains 4 query sets representing the most over- and under-expressed genes from each gene signature (default is 500 genes per query set)
#### 4) GSEA between gene signatures and query gene sets - returns enrichment data (normalized scores, nominal p-values, plots) and leading-edge gene identification

#### NOTE: This algorithm computes T-score per scipy.stats.ttest_ind (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html). This algorithm does not adjust T-score values using a minimum value for σ of .2 * absolute(μ), where μ=0 is adjusted to μ=1, like performed in the GSEA desktop version (https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html).
#### Required dependencies: numpy, pandas, scipy, and gseapy (https://gseapy.readthedocs.io/)

#### Limitations of these programs include:
#### 1) Manual preparation of input files prior to program use to ensure clean and complete data that consistently uses the same gene identifier (_e.g._, probes, gene symbols, etc.) across all datasets
#### 2) No GUI = manual manipulation of file names, column indicies, and query set sizes in the code as designated by the "#change ... as needed" comments
#### Please note that Python uses 0 for the index of the first column in a dataframe. Adjust your column counts accordingly.
#### 3) Manual analysis of identified leading-edge genes and preparation of gene panels file, similar to .gmt format, for verification and comparison stages
#### 4) Manual analysis to verify individual panel genes from identified leading-edge genes

#### Coming soon!
#### GSEAGeneID.py, which is similar to GSEA1datasetGeneIDonly.py except it inputs two datasets and defines one gene signature from each dataset separately rather than inputing one dataset and dividing it to generate the two gene signatures.
