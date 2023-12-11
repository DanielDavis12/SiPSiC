\name{normalizeCountsForGene}
\alias{normalizeCountsForGene}
\title{
Gene counts normalization
}
\description{
Get the counts of a single gene normalized by the median of the top 5 percent cells, unless it's zero; In this case, the counts are all divided by the maximum value across all cells. If all counts are zeros, they are returned untouched.
}
\usage{
normalizeCountsForGene(expressionValues, percentForNormalization)
}
\arguments{
  \item{expressionValues}{
An array of type double, containing the counts (in any units, e.g. CPM or TPM) of a single gene across different cells.
}
  \item{percentForNormalization}{
The percent of top cells for each gene whose median is used to normalized the gene's expression values.
}
}
\value{
An array (type double) of the normalized input counts.
}
\author{
Daniel Davis, Yotam Drier
}