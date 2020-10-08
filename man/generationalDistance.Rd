\name{generationalDistance}
\alias{generationalDistance}
\alias{generalizedSpread}
\alias{dominatedHypervolume}
\alias{epsilonIndicator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quality measures for MCO solutions}
\description{
  Functions to evaluate the quality of the estimated pareto front.
}
\usage{
generationalDistance(x, o)
generalizedSpread(x, o)
epsilonIndicator(x, o)
dominatedHypervolume(x, ref)
}
\arguments{
  \item{x}{Estimated pareto front or an object which has a paretoFront method}
  \item{o}{True pareto front or an object which has a paretoFront
    method}
  \item{ref}{Reference point (may be omitted).}
}
\details{
  Instead of the pareto front, one can also pass an object for which a
  paretoFront method exists to both methods.

  For \code{dominatedHypervolume}, if no reference point is given, the
  maximum in each dimension is used as the reference point.
}
\note{
  This code uses version 1.3 of the hypervolume code available from
  \url{http://iridia.ulb.ac.be/~manuel/hypervolume}. For a description
  of the algorithm see

  Carlos M. Fonseca, Luis Paquete, and Manuel Lopez-Ibanez.
  \emph{An improved dimension-sweep algorithm for the hypervolume indicator}.
  In IEEE Congress on Evolutionary Computation, pages 1157-1163,
  Vancouver, Canada, July 2006.
}
\value{
  The respective quality measure.
}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\references{
  Carlos M. Fonseca, Luis Paquete, and Manuel Lopez-Ibanez.
  \emph{An improved dimension-sweep algorithm for the hypervolume indicator}.
  In IEEE Congress on Evolutionary Computation, pages 1157-1163,
  Vancouver, Canada, July 2006.

  Nicola Beume, Carlos M. Fonseca, Manuel Lopez-Ibanez, Luis Paquete,
  and J. Vahrenhold. \emph{On the complexity of computing the
  hypervolume indicator}. IEEE Transactions on Evolutionary Computation,
  13(5):1075-1082, 2009.
  
  Zitzler, E., Thiele, L., Laumanns, M., Fonseca, C., and Grunert da
  Fonseca, V (2003): \emph{Performance Assessment of Multiobjective
  Optimizers: An Analysis and Review}. IEEE Transactions on Evolutionary
  Computation, 7(2), 117-132.
}
\examples{
## Estimate true front:
\dontrun{
tf <- nsga2(fonseca2, 2, 2,
            lower.bounds=c(-4, -4), upper.bounds=c(4, 4),
            popsize=1000, generations=100)
res <- nsga2(fonseca2, 2, 2,
             lower.bounds=c(-4, -4), upper.bounds=c(4, 4),
             popsize=16, generations=c(2, 4, 6, 8, 10, 20, 50))
n <- length(res)
sapply(1:n, function(i) dominatedHypervolume(res[[i]], c(1, 1)))
sapply(1:n, function(i) generationalDistance(res[[i]], tf))
sapply(1:n, function(i) generalizedSpread(res[[i]], tf))
sapply(1:n, function(i) epsilonIndicator(res[[i]], tf))
}
}
