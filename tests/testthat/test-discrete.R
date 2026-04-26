context("Discrete margins")

set.seed(5)

test_that("svinecop accepts doubled-column discrete data", {
  n <- 200
  lambda <- 3
  d <- 2

  x <- matrix(rpois(n * d, lambda), n, d)
  udisc <- cbind(
    ppois(x, lambda),
    ppois(x - 1, lambda)
  )

  expect_no_error(
    svinecop(
      udisc,
      p = 1,
      var_types = c("d", "d"),
      family_set = "gaussian"
    )
  )
})
