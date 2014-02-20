\name{summaryCopa}
\alias{summaryCopa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create Summary Showing Top Gene Pairs }
\description{
 This function can be used to output a data.frame containing the ID and
 optionally the gene symbol for the top gene pairs, based on the number
 of outliers.
}
\usage{
summaryCopa(copa, pairnum, lib = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{copa}{An object of class 'copa', resulting from a call to the \code{copa}
    function. }
  \item{pairnum}{ The maximum number of outlier pairs to be output. A
    table can be output first using \code{tableCopa}}
  \item{lib}{For Affymetrix data that have an annotation package, this
    can be specified and the table will then also contain the gene symbol}
}

\value{
 The output from this function is a \code{data.frame} with the number of
 outliers, the manufacturer identifiers, and optionally, the gene symbol
 for the genes.
}
\references{
  Tomlins, SA, et al. Recurrent fusion of TMPRSS2 and ETS transcription
  factor genes in prostate cancer.
  Science. 2005 Oct 28;310(5748):644-8.}
\author{James W. MacDonald <jmacdon@u.washington.edu>}

\examples{
if(interactive()){
library(Biobase)
data(sample.ExpressionSet)
cl <- abs(3 - as.numeric(pData(sample.ExpressionSet)[,2]))
tmp <- copa(sample.ExpressionSet, cl)
summaryCopa(tmp, 6)
}
}
\keyword{ manip}

