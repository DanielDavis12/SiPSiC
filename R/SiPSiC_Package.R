# SiPSiC package - developed at the Drier's lab, Lautenberg center for immunology and cancer, the Hebrew University of Jerusalem

PERCENTAGE_FOR_COUNTS_NORMALIZATION <- 5


#' Gene counts normalization
#'
#' Get the counts of a single gene normalized by the median of the top 5 percent cells.
#' @param expressionValues An array of type double, containing the counts (in any units, e.g. CPM or TPM) of a single gene across different cells.
#' @return An array (type double) of the input counts divided by the median of the top 5 percent cells, unless it's zero; in this case, the counts are all divided by the maximum value across all cells. If all counts are zeros, they are returned untouched.
normalizeCountsForGene <- function(expressionValues)
{
  geneCountsRanked <- rank(x = expressionValues, na.last = "keep", ties.method = 'max')
  rankOfNormalizationPercentage <- ((100 - PERCENTAGE_FOR_COUNTS_NORMALIZATION) / 100) * 
    length(expressionValues)
  topCellCounts <- expressionValues[geneCountsRanked >= rankOfNormalizationPercentage]
  normalizationScore <- stats::median(topCellCounts)
  
  # In case too few cells express the gene to any extent - the normalization score is adapted
  if (normalizationScore == 0)
  {
    maxExpressionValue <- max(expressionValues)
    
    # In case non of the cells express the gene - the original expression values are returned
    if (maxExpressionValue == 0)
    {
      geneScoresForCells <- expressionValues
    }
    else
    {
      geneScoresForCells <- expressionValues / maxExpressionValue
    }
  }
  else
  {
    geneScoresForCells <- (expressionValues / normalizationScore)   
  }
  
  return (as.vector(geneScoresForCells))
}

#' Pathway scores calculation
#'
#' Calculate the scores of a given pathway for every cell in a single-cell RNA-seq data.
#' @param dataMatrix a matrix whose rows are genes and columns are cells, containing the gene counts; Counts-Per-Million (CPM) or Transcripts-Per-Million (TPM) units are recommended. The matrix should be of type sparse matrix ("dgCMatrix"), otherwise errors might be raised.
#' @param pathwayGenes a character vector of the gene names of which the relevant biological pathway consists.
#' @return a list containing two arrays: "$pathwayScores", an array (type double) of the pathway score of each cell in the input dataMatrix, and "$index", a numeric array of the indices in the dataMatrix where genes of the pathway were found.
#' @examples 
#' rowIndex <- colIndex <- matValues <- c(1:4)
#' geneCountsMatrix <- Matrix::sparseMatrix(i = rowIndex, j = colIndex, x = matValues)
#' rownames(geneCountsMatrix) <- c("Gene1", "Gene2", "Gene3", "Gene4")
#' colnames(geneCountsMatrix) <- c("Cell1", "Cell2", "Cell3", "Cell4")
#' pathwayGenesList <- c("Gene1", "Gene2", "Gene4")
#' scoresAndIndices <- getPathwayScores(geneCountsMatrix, pathwayGenesList)
#' pathwayScoresOfCells <- scoresAndIndices$pathwayScores
#' pathwayGeneIndices <- scoresAndIndices$index
#' @export
getPathwayScores <- function(dataMatrix, pathwayGenes)
{
  # Fetching only those genes that belong to the pathway
  allGenesList <- rownames(dataMatrix)
  pathwayGeneIndices <- match(pathwayGenes,allGenesList)
  pathwayGeneIndices <- pathwayGeneIndices[!is.na(pathwayGeneIndices)]
  if(length(pathwayGeneIndices) == 0) {return(NA)}
  pathwayOnlyData <- dataMatrix[pathwayGeneIndices,]
  
  totalCellsCount <- ncol(pathwayOnlyData)
  
  # Calculating the cell scores for each gene separately
  numOfPathwayGenes <- nrow(pathwayOnlyData)
  geneScoresMatrix <- matrix(0, numOfPathwayGenes, totalCellsCount)
  rownames(geneScoresMatrix) <- rownames(pathwayOnlyData)
  colnames(geneScoresMatrix) <- colnames(pathwayOnlyData)
  geneScoresLists <- lapply(X = split(pathwayOnlyData, row(pathwayOnlyData)), normalizeCountsForGene)
  geneScoresMatrix[] <- do.call("rbind", geneScoresLists)
  
  # Below - calculating the pathway score of each cell as a weighted average of all its gene scores
  
  # Using the total counts of all the genes to rank them
  totalCountsPerGene <- apply(pathwayOnlyData, 1, sum)
  genesRanking <- rank(x = totalCountsPerGene, na.last = "keep", ties.method = 'max')
  genesRanking <- genesRanking / numOfPathwayGenes
  geneScoresMatrix <- geneScoresMatrix * genesRanking
  
  # Setting the average of all weighted genes of a cell as its pathway score
  geneScoreSums <- apply(geneScoresMatrix, 2, sum, na.rm = TRUE)
  pathwayScores <- geneScoreSums / numOfPathwayGenes
  rm(geneScoresMatrix)
  gc()
  
  outputList <- list("index"= pathwayGeneIndices, "pathwayScores" = pathwayScores)
  return (outputList)
}