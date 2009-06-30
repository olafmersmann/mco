\name{nsga2}
\alias{nsga2}
\title{NSGA II MOEA}
\description{
  Multicriterion optimization algorithm
}
\usage{
nsga2(fn, idim, odim, ..., constraints = NULL, cdim = 0, lower.bounds = rep(-Inf, idim), upper.bounds = rep(Inf, idim), popsize = 100, generations = 100, cprob = 0.7, cdist = 5, mprob = 0.2, mdist = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{Function}
  \item{idim}{Input dimension}
  \item{odim}{Output dimension}
  \item{\dots}{Arguments passed through to 'fn'}
  \item{constraints}{Constraint function}
  \item{cdim}{Constraint dimension}
  \item{lower.bounds}{Lower bound of input}
  \item{upper.bounds}{Upper bound of input}
  \item{popsize}{Size of population}
  \item{generations}{Number of generations to breed. If a vector, then
    the result will contain the population at each given generation.}
  \item{cprob}{Crossing probability}
  \item{cdist}{Crossing distribution index}
  \item{mprob}{Mutation probability}
  \item{mdist}{Mutation distribution index}
}
\details{
  
}
\value{
  If generation is an integer, a list describing the final population with components \code{par},
  \code{value} and \code{pareto.optimal}. If generations is a vector, a
  list is returned. The i-th element of the list contains the population
  after generations[i] generations.
}
\references{
  Deb, K., Pratap, A., and Agarwal, S..
  A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II.
  \emph{IEEE Transactions on Evolutionary Computation}, \bold{6 (8)}
  (2002), 182-197.
}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\examples{
## Binh 1 problem:
binh1 <- function(x) {
  y <- numeric(2)
  y[1] <- crossprod(x, x)
  y[2] <- crossprod(x - 5, x - 5)
  return (y)
}
r1 <- nsga2(binh1, 2, 2,
           generations=150, popsize=100,
           cprob=0.7, cdist=20,
           mprob=0.2, mdist=20,
           lower.bounds=rep(-5, 2),
           upper.bounds=rep(10, 2))
plot(r1)

## VNT problem:
vnt <- function(x) {  
  y <- numeric(3)
  xn <- crossprod(x, x)
  y[1] <- xn/2 + sin(xn);
  y[2] <- (crossprod(c(3, -2), x) + 4)^2/8 + (crossprod(c(1, -1), x) + 1)^2/27 + 15
  y[3] <- 1/(xn + 1) - 1.1*exp(-xn)
  return (y)
}

r2 <- nsga2(vnt, 2, 3,
           generations=150, popsize=100,
           lower.bounds=rep(-3, 2),
           upper.bounds=rep(3, 2))
plot(r2)
}
\keyword{optimize}
