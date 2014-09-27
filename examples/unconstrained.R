library("mco")

## Binh 1 problem:
binh1 <- function(x) {
  y <- numeric(2)
  y[1] <- crossprod(x, x)
  y[2] <- crossprod(x - 5, x - 5)
  return (y)
}

vectorized_binh1 <- function(X) {
  apply(X, 1, binh1)
}


r1v <- nsga2(vectorized_binh1, 2, 2,
             generations=150, popsize=100,
             cprob=0.7, cdist=20,
             mprob=0.2, mdist=20,
             vectorized=TRUE)
