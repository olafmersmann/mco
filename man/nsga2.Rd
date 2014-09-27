\name{nsga2}
\alias{nsga2}
\title{NSGA II MOEA}
\description{
  The NSGA-II algorithm minimizes a multidimensional function to approximate its Pareto front and Pareto set.
  It does this by successive sampling of the search space, each such sample is called a \emph{population}. 
  The number of samples taken is governed by the \code{generations} parameter, the size of the sample by the \code{popsize} parameter.
  Each population is obtained by creating so called offspring search points from the best individuals in the previous population. 
  The best individuals are calculated by non-dominated sorting breaking ties using the crowding distance.
  The total number of function evaluations used is 
    \deqn{n_eval = popsize * (generations + 1)}
  when \code{generations} is a single number and
    \deqn{n_eval = popsize * (max(generations) + 1)}
  when \code{generations} is a vector of numbers. 
}
\usage{
  nsga2(fn, idim, odim, ...,
        constraints = NULL, cdim = 0,
        lower.bounds = rep(-Inf, idim), upper.bounds = rep(Inf, idim),
        popsize = 100, generations = 100,
        cprob = 0.7, cdist = 5,
        mprob = 0.2, mdist = 10,
        vectorized=FALSE)
}
\arguments{
  \item{fn}{Function to be minimized}
  \item{idim}{Input dimension}
  \item{odim}{Output dimension}
  \item{\dots}{Arguments passed through to 'fn'}
  \item{constraints}{Constraint function}
  \item{cdim}{Constraint dimension}
  \item{lower.bounds}{Lower bound of parameters}
  \item{upper.bounds}{Upper bound of parameters}
  \item{popsize}{Size of population}
  \item{generations}{Number of generations to breed. If a vector, then
    the result will contain the population at each given generation.}
  \item{cprob}{Crossover probability}
  \item{cdist}{Crossover distribution index}
  \item{mprob}{Mutation probability}
  \item{mdist}{Mutation distribution index}
  \item{vectorized}{If \code{TRUE}, the objective and constraint
    functions must be vectorized, i.e. accept a matrix instead of a
    vector and return a matrix instead of a vector. The matrix is
    structured such that one individual parameter combination is
    contained in each row (the matrix has shape \code{popsize * idim})
    and each objective is stored in a column of the returned matrix (the
    returned matrix must have shape \code{popsize * odim})}
}
\value{
  If generation is an integer, a list describing the final population
  with components \code{par}, \code{value} and \code{pareto.optimal}. If
  generations is a vector, a list is returned. The i-th element of the
  list contains the population after generations[i] generations.  
}
\references{
  Deb, K., Pratap, A., and Agarwal, S..
  A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II.
  \emph{IEEE Transactions on Evolutionary Computation}, \bold{6 (8)}
  (2002), 182-197.
}
\seealso{\code{\link{zdt1}} for more examples and a list of
  multiobjective test functions.
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

## Example using constraints:
## minimize     f(x) = (x[1]^2, x[2]^2)
## subject to   g(x) = (sum(x) - 5) >= 0
f <- function(x) { x^2 }
g <- function(x) { sum(x) - 5 }
res <- nsga2(f, 2, 2, generations=500,
             lower.bounds=c(0, 0), upper.bounds=c(10, 10),
             constraints=g, cdim=1)
opar <-par(mfrow=c(1,2))
plot(res, xlab="y1", ylab="y2", main="Objective space")
plot(res$par, xlab="x1", ylab="x2", main="Parameter space")
par(opar)
}
\keyword{optimize}
