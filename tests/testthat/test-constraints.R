context("constraint")

test_that("Simple constraints work", {
  f <- function(x) c(sum(x^2), sum((x - 2)^2))
  g <- function(x) 5 - sum(x)

  res <- nsga2(f, idim=2, odim=2, generations=500,
               lower.bounds=c(0, 0), upper.bounds=c(10, 10),
               constraints=g, cdim=1)

  expect_true(all(apply(res$par, 1, g) >= 0))
})

test_that("Complex constraints work", {
  f <- function(x) c(sum(x^2), sum((x - 2)^2))
  g <- function(x) c(1 - sqrt(sum((x-2)^2)), 4 - sum(x))

  res <- nsga2(f, idim=2, odim=2, generations=500,
               lower.bounds=c(0, 0), upper.bounds=c(10, 10),
               constraints=g, cdim=2)

  expect_true(all(apply(res$par, 1, g) >= 0))
})
