### GSEA-based Programs for NGS data
#### This repository contains Python programs and example input files to reproduce the Gene Set Enrichment Analysis (GSEA)-based approach for gene identification and verification presented in the Park & Harris paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8342995/). 

#### GSEAGeneID.py = Program inputs two datasets containing clean and complete data then performs the following:
#### 1) Z-score normalization of data across all samples for each gene in the dataset - returns two z-scored datasets
#### 2) Generation of gene signatures via T-score - returns two gene lists, one per dataset, with associated T-scores and p-values
#### 3) Generation of query gene sets - returns a .txt file similar to a .gmt (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) that contains 10 query sets representing a range of the five most over- and under-expressed genes from each gene signature
#### 4) GSEA between gene signatures and query gene sets - returns enrichment data (normalized scores, nominal p-values, plots) and leading-edge gene identification

#### Required dependencies: numpy, pandas, scipy, and gseapy (https://gseapy.readthedocs.io/)

#### Limitations of these programs include:
#### 1) Manual preparation of input files prior to program use to ensure clean and complete data that consistently uses the same gene identifier (_e.g._, probes, gene symbols, etc.) across all datasets
#### 2) No GUI = manual manipulation of file names and column indicies in the code as designated by the "#change ... as needed" comments
#### 3) Manual analysis of identified leading-edge genes and preparation of gene panels file, similar to .gmt format, for verification and comparison stages
