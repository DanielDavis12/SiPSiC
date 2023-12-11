\name{getPathwayScores}
\alias{getPathwayScores}
\title{Calculate Pathway Scores for all Cells in a scRNA-Seq Data}
\description{
Calculate a pathway score for each cell included in the input scRNA-seq data, using the cell's transcription level of the pathway's genes.
}
\usage{
getPathwayScores(dataMatrix, pathwayGenes, percentForNormalization)
}
\arguments{a list of three items:
  \item{dataMatrix}{
a matrix whose rows and columns represent genes and cells, respectively, containing the scRNA-seq data; Counts-Per-Million (CPM) or Transcripts-Per-Million (TPM) units should be used for best results. The matrix should be of type sparse matrix (dgCMatrix), otherwise errors might be raised.
}
  \item{pathwayGenes}{
a character vector of the gene names of which the relevant biological pathway consists.
}
  \item{percentForNormalization}{
The percent of top cells for each gene whose median is used to normalize the gene's expression values (5 by default).
}
}
\value{
\item{pathwayScores}{an array (type double) of the pathway score of each cell in the input dataMatrix, unless none of the pathway genes is found in the data, in which case NA is returned.}
\item{index}{a numeric array of the row indices in the dataMatrix where genes of the pathway were found.}
}
\references{
\url{https://medicine.ekmd.huji.ac.il/en/research/yotamd/Pages/default.aspx}
}
\author{
Daniel Davis, Yotam Drier
}
\examples{
library(SiPSiC)
geneCountsMatrix <- matrix(rpois(16, lambda = 10), ncol = 4, nrow = 4)
geneCountsMatrix <- as(geneCountsMatrix, "dgCMatrix")
rownames(geneCountsMatrix) <- c("Gene1", "Gene2", "Gene3", "Gene4")
colnames(geneCountsMatrix) <- c("Cell1", "Cell2", "Cell3", "Cell4")
assayData <- SingleCellExperiment(assays = list(counts = geneCountsMatrix))
pathwayGenesList <- c("Gene1", "Gene2", "Gene4")
percentForNormalization <- 7
scoresAndIndices <- getPathwayScores(counts(assayData), pathwayGenesList, percentForNormalization)
pathwayScoresOfCells <- scoresAndIndices$pathwayScores
pathwayGeneIndices <- scoresAndIndices$index
}