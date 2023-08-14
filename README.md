# SiPSiC
TL;DR: Cell-specific pathway score calculation from single-cell RNA sequencing data.

The Single Pathway analysis in Single Cells (SiPSiC) package is used to infer biological pathway activity of the individual cells included in a single-cell RNA-seq dataset. 
The dataset and list of genes comprising the pathway of interest are provided by the user, then per-cell pathway scores are calculated for all cells by the 'getPathwayScores' function.
SiPSiC depends on R's 'Matrix' package available on CRAN and on Bioconductor's SingleCellExperiment package. For further information, go to:

[https://www.biorxiv.org/content/10.1101/2023.03.27.534310v1](https://doi.org/10.1101/2023.03.27.534310)

Install the package directly from Bioconductor by executing the following commands in an R session:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SiPSiC")

SiPSiC was developed at the Drier lab, the Lautenberg Center for Immunology and Cancer Research, IMRIC, Faculty of Medicine, Hebrew University of Jerusalem, Israel.
