\name{paretoFront}
\alias{paretoFront}
\alias{paretoSet}
\alias{paretoFilter}
\title{Pareto Front and pareto set accessor}
\description{
  Extract the pareto front or pareto set from an mco result object.

  Filter an mco result and extract the pareto-optimal solutions.
}
\usage{
paretoFront(x, ...)
paretoSet(x, ...)
paretoFilter(x, ...)
}
\arguments{
  \item{x}{matrix or mco result object}
  \item{\dots}{Ignored}
}
\value{
  A matrix containing the pareto front or pareto set.

  \code{paretoFilter} returns those values in \code{x} which are not
  dominated by any other solution.
}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
