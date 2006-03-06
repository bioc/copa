\name{tableCopa}
\alias{tableCopa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize copa results }
\description{
  This function will output a table showing the number of gene pairs at
  each number of outliers.
}
\usage{
tableCopa(copa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{copa}{A 'copa' object, the result of a call to \code{copa} }
}

\value{
 This function simply prints a table to the screen, useful for
 summarizing the output from a call to \code{copa}.
}

\author{James W. MacDonald}
\examples{
library(Biobase)
data(sample.exprSet)
cl <- abs(3 - as.numeric(pData(sample.exprSet)[,2]))
tmp <- copa(sample.exprSet, cl)
tableCopa(tmp)
}
\keyword{ manip }% at least one, from doc/KEYWORDS

