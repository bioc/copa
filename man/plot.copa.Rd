\name{plot.copa}
\alias{plot.copa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Gene Pairs fom the Results of Running copa}
\description{
  This function can be used to visualize pairs of genes that may be
  involved in recurrent gene fusion in cancer.
}
\usage{
plot(copa, idx, lib = NULL,  sort = TRUE, col = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{copa}{ An object of class 'copa', resulting from a call to the \code{copa}
    function.}
  \item{idx}{ A numeric vector listing the gene pairs to plot (e.g., idx
    = 1:3 will plot the first three gene pairs).}
   \item{lib}{If the underlying data are Affymetrix expression values,
    one can specify an annotation package and the plot labels will be
    extracted from the xxxSYMBOL environment. If \code{NULL}, the
    \code{row.names} of the gene expression matrix will be used.}
  \item{sort}{Boolean. Should the data be sorted before plotting?
    Defaults to \code{TRUE}.}
  \item{col}{ A vector of color names or numbers to be used for coloring
    the different samples in the resulting barplot.}
}
\details{
  Note that this function will output all the gene pairs in the idx
  vector without pausing. This can be controlled by either setting
  par(ask = TRUE), or by redirecting the output to a file (using e.g.,
  \code{pdf}, \code{ps}, etc.).
}
\value{
  This function is called solely for outputting plots. No values are returned.
}
\references{
  Tomlins, SA, et.al. Recurrent fusion of TMPRSS2 and ETS transcription
  factor genes in prostate cancer.
  Science. 2005 Oct 28;310(5748):644-8.}
\author{ James W. MacDonald}
\examples{
if(interactive()){
library(Biobase)
data(example.exprSet)
tmp <- copa(example.exprSet)
plot(tmp, 1, pData(example.exprSet)[,2], col = c("red", "blue"))
}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS


