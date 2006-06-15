\name{copaPerm}
\alias{copaPerm}
\title{ Measure Significance of COPA by Permutation}
\description{
  This function can be used to determine the significance of the results
  that one gets from running \code{copa} on a particular dataset, based
  on permuting the class assignments.
}
\usage{
copaPerm(object, copa, outlier.num, gene.pairs, B = 100, pval = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An \code{exprSet}, \code{ExpressionSet}, or a matrix or \code{data.frame}.}
  \item{copa}{ An object of class 'copa', produced by running
    \code{copa} on a set of microarray data.}
  \item{outlier.num}{The number of outliers to test for. See details for
    more information}
  \item{gene.pairs}{The number of gene pairs to test for. See details
    for more information}
  \item{B}{The number of permutations to perform. Defaults to 100. This
    may be too many for interactive use.}
  \item{pval}{Boolean. Output an estimated p-value and false discovery
    rate? Defaults to \code{FALSE}. This result will only be reasonable
    for large numbers of permutations (500 - 1000). See details.}
  \item{verbose}{Boolean. Print out the permutation number at each of
    100, 200, etc. Defaults to \code{TRUE}}
}
\details{
  Running \code{copa} on a set of microarray data will result in the
  output of an object of class 'copa', which is a list containing (among
  other things) an ordered vector that lists the number of mutually
  exclusive outlier samples for various gene pairs. This vector is
  ordered from smallest to largest following the assumption that the
  gene pairs with the most mutually exclusive outliers are probably more
  likely to be involved in some sort of recurrent fusion.

  One can see how many pairs of genes resulted in a given number of
  outliers by calling \code{tableCopa}. One may then
  want to determine how significant a certain number of pairs is (e.g.,
  how likely is it to get that many pairs if there is no recurrent
  fusion occuring). The most straightforward way to estimate the
  significance of a given result is to repeatedly permute the
  classlabels and see how many times one gets a result as large or
  larger than what was observed.

  Technically speaking, to get a reasonable estimate of significance and
  a false discovery rate, one would need to permute 500 - 1000
  times. However, this can take an inordinate amount of time (best left
  for an overnight run). To get a quick idea of significance, one could
  simply permute maybe 10 times (with pval = FALSE) to see how likely it
  is to get a certain number of outliers.
}
\value{
  \item{out}{A vector listing the number of gene pairs with at least as
    many outliers as 'num.outlier'.}
  \item{p.value}{A permuted p-value, only output if pval = TRUE. Note
    that the size of the p-value is determined by both the number of
    outliers >= 'num.outlier' as well as the number of permutations, so
    too  few permutations may result in
    a p-value that doesn't look very significant even if it is.}
  \item{fdr}{The expected number of gene pairs with at least as many
    outliers as 'num.outlier'. This can be converted to a \%FDR by
    dividing by the observed value.}
}
\references{
Tomlins, SA, et al. Recurrent fusion of TMPRSS2 and ETS transcription
factor genes in prostate cancer.
Science. 2005 Oct 28;310(5748):644-8. 
}
\author{ James W. MacDonald}
\keyword{ manip }% at least one, from doc/KEYWORDS

