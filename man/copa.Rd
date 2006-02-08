\name{copa}
\alias{copa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate COPA Scores from a Set of Microarrays}
\description{
 This function calculates COPA scores from a set of microarrays. Input
 can be an \code{exprSet}, \code{PLMset}, \code{MAList},
 \code{marrayNorm} object, or a matrix or \code{data.frame}.
}
\usage{
copa(object, cl, cutoff = 5, max.overlap = 0, norm.count = 0, pct = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An \code{exprSet}, \code{PLMset}, \code{MAList},
    \code{marrayNorm} object, or a matrix or \code{data.frame}.}
  \item{cl}{A vector of classlabels indicating sample status (normal =
    1, tumor = 2).}
  \item{cutoff}{The cutoff to determine 'outlier' status. See details
    for more information.}
  \item{max.overlap}{ The maximum number of samples that can be
    considered 'outliers' when comparing two genes. The default is 0,
    indicating that there can be no overlap. See details for more
    information.}
  \item{norm.count}{The number of normal samples that can be considered
    'outliers'. The default is 0, meaning that no normals may be
    outliers.}
  \item{pct}{The percentile to use for pre-filtering the data. A
    preliminary step is to compute the number of outlier samples for
    each gene. All genes with a number of outlier samples less than the
    (default 95th) percentile will be removed from further consideration.}
}
\details{
  Cancer Outlier Profile Analysis is a method that is intended to find
  pairs of genes that may be involved in recurrent gene fusion with a
  third (unknown) gene. The underlying idea here is that in certain
  cancers it may be common for the promoter region of one gene to become
  fused to certain oncogenes. For instance, Tomlins et. al. showed that
  the promoter region of TMPRSS2 fused to either ERG or ETV1 in the
  majority of prostate cancer tumors tested.

  Since this fusion should only happen with one oncogene in a given sample, we
  look for pairs of genes where some samples have much higher expression
  values, but the samples for gene 'A' are mutually exclusive from the
  samples for gene 'B'.

  The cutoff argument for this function is used to determine how high
  the centered and scaled expression value has to be in order to be
  considered an outlier. The max.overlap argument allows one to relax
  the requirement of mutual exclusivity, although in practice this is
  probably not advisable.

  Note that this function computes all row-wise comparisons, which gets
  very large very quickly. The function will throw a warning for any
  data set containing > 1000 rows and query the user to see if he/she
  really wants to proceed. The number of genes to be considered can be
  adjusted by increasing/decreasing the 'pct' argument.
}
\value{
\item{ord.prs}{A matrix with two columns containing the ordered row
  numbers from the original matrix of gene expression values.}
\item{pr.sums}{A numeric vector with the number of mutually exclusive
  outliers for each gene pair. This is the criterion for ranking the
  gene pairs; the assumption being that a pair of genes with more
  mutually exclusive outliers will be more interesting than a pair with
  relatively fewer outliers.}
\item{mat}{A matrix containing the filtered gene expression values.}
\item{cl}{The classlabel vector passed to \code{copa}}
\item{cutoff}{The cutoff used}
\item{max.overlap}{The value of max.overlap used}
\item{norm.count}{The value of norm.count used}
\item{pct}{The percentile used in the pre-filtering step}
}
\references{
Tomlins, SA, et.al. Recurrent fusion of TMPRSS2 and ETS transcription
factor genes in prostate cancer.
Science. 2005 Oct 28;310(5748):644-8. 
}
\author{ James W. MacDonald}

\examples{
library(Biobase)
data(example.exprSet)
tmp <- copa(example.exprSet)
}

\keyword{ univar }% at least one, from doc/KEYWORDS

