\name{normalizeCountsForGene}
\alias{normalizeCountsForGene}
\title{
Gene counts normalization
}
\description{
Get the counts of a single gene normalized by the median of the top 5 percent cells, unless it's zero; In this case, the counts are all divided by the maximum value across all cells. If all counts are zeros, they are returned untouched.
}
\usage{
normalizeCountsForGene(expressionValues)
}
\arguments{
  \item{expressionValues}{
An array of type double, containing the counts (in any units, e.g. CPM or TPM) of a single gene across different cells.
}
}
\value{
An array (type double) of the normalized input counts.
}
\author{
Daniel Davis, Yotam Drier
}
\alias{PERCENTAGE_FOR_COUNTS_NORMALIZATION}
\docType{data}
\format{
  The format is:
 num 5
}