\name{getans}
\alias{getans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Interactive Function}
\description{
  A function to query the end user. This is an internal function and not
  intended to be called directly by the end user.
}
\usage{
getans(msg, allowed = c("y", "n"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{msg}{ The query.}
  \item{allowed}{Allowed responses}
}

\value{
  The response is returned.
}

\author{ James W. MacDonald}

\keyword{ internal}

