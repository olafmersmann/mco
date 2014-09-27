library("mco")

f <- function(x) {
  x^2
}

g <- function(x) {
  sum(x) - 5
}

h <- function(x) {
  c(sum(x) - 5, prod(x) - 10)
}

res <- nsga2(f, 2, 2, generations=500,
             lower.bounds=c(0, 0),
             upper.bounds=c(10, 10),
             constraints=g, cdim=1)

res <- nsga2(f, 2, 2, generations=500,
             lower.bounds=c(0, 0),
             upper.bounds=c(10, 10),
             constraints=h, cdim=2)
