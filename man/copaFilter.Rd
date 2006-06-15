\name{copaFilter}
\docType{methods}
\alias{copaFilter}
\alias{do.copaFilter}
\alias{copaFilter-methods}
\alias{copaFilter, matrix-method}
\alias{copaFilter, data.frame-method}
\alias{copaFilter, exprSet-method}
\alias{copaFilter, ExpressionSet-method}
\title{ Pre-filter Genes for COPA Analysis}
\description{
  This function is used to pre-filter genes prior to doing a COPA
  analysis. The filtering is based on the nth percentile of the outlier
  samples for each gene. This function is an internal function and not
  intended to be called by the end user.
}
\usage{
\S4method{copaFilter}{matrix}(object, cl, cutoff, norm.count, pct)
\S4method{copaFilter}{data.frame}(object, cl, cutoff, norm.count, pct)
\S4method{copaFilter}{exprSet}(object, cl, cutoff, norm.count, pct)
\S4method{copaFilter}{ExpressionSet}(object, cl, cutoff, norm.count, pct)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object} {An \code{exprSet}, \code{ExpressionSet}, or a matrix or \code{data.frame}.}
  \item{cl}{ A vector of classlabels indicating sample status (normal =
    1, tumor = 2).}
  \item{cutoff}{The cutoff to determine 'outlier' status. See details
    for more information.}
  \item{norm.count}{ The number of normal samples that can be considered
    'outliers'. The default is 0, meaning that no normals may be
    outliers.}
  \item{pct}{ The percentile to use for pre-filtering the data. A
    preliminary step is to compute the number of outlier samples for
    each gene. All genes with a number of outlier samples less than the
    (default 95th) percentile will be removed from further consideration.}
}

\value{
  \item{mat}{A matrix containing the gene expression values for the
    filtered genes.}
}
\references{
  Tomlins, SA, et al. Recurrent fusion of TMPRSS2 and ETS transcription
  factor genes in prostate cancer.
  Science. 2005 Oct 28;310(5748):644-8. 
}
\author{ James W. MacDonald}


\keyword{ internal}% at least one, from doc/KEYWORDS
