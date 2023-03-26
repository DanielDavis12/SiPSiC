# SiPSiC
TL;DR: Cell-specific pathway score calculation from single-cell RNA sequencing data.

The Single Pathway analysis in Single Cells (SiPSiC) package is used to infer biological pathway activity of the individual cells included in a single-cell RNA-seq dataset. 
The dataset and list of genes comprising the pathway of interest are provided by the user, then per-cell pathway scores are calculated for all cells by the 'getPathwayScores' function.
SiPSiC depends on R's 'Matrix' package available on CRAN and on Bioconductor's SingleCellExperiment package.

To install the package use:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SiPSiC")



Developed at the Drier lab, Lautenberg center for immunology and cancer, the Hebrew University of Jerusalem
