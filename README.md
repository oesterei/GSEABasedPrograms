### GSEA-based Programs for NGS data
#### This repository contains Python programs and example input files to reproduce the GSEA-based approach for gene identification, verification, and comparison presented in the Park & Harris paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8342995/). 
#### GSEAGeneID.py = Program inputs user specified datasets then prompts the user for the following:
#### 1) z-score across all samples in each dataset - returns two z-scored datasets
#### 2) T-score for user specified samples - returns two gene lists, one per dataset, with associated T-scores and p-values
#### 3) Generate query sets - produces a .txt file similar to a .gmt (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) that contains XX query sets representing the most over- and under-expressed genes from each gene signature

#### Required dependencies: numpy, pandas, scipy, and gseapy (https://gseapy.readthedocs.io/)


#### Limitations of these programs include:
#### 1) Manual preparation of input files to ensure clean, complete data that consistently uses the same gene identifier (_e.g._, probes, gene symbols, etc.) across all datasets
#### 2) Manual identification of column index for selected samples in the code 
#### 3) No GUI = manual manipulation of variable and/or file names marked by the "#change as needed" comment in the code
