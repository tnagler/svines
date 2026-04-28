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

test_that("svinecop_dist accepts discrete var_types", {
  bc <- bicop_dist("gaussian", parameters = 0.5)
  sv <- svinecop_dist(
    pair_copulas = list(list(bc)),
    cs_structure  = cvine_structure(1:2),
    p             = 0,
    out_vertices  = 1:2,
    in_vertices   = 1:2,
    var_types     = c("d", "d")
  )
  expect_s3_class(sv, "svinecop_dist")
})

test_that("svinecop with cs_structure and omitted var_types defaults to continuous of length cs_dim", {
  set.seed(1)
  n <- 200
  d <- 2
  u <- matrix(runif(n * d), n, d)
  cs <- cvine_structure(1:d)
  sv <- svinecop(u, p = 1, cs_structure = cs,
                 out_vertices = 1:d, in_vertices = 1:d)
  expect_length(sv$var_types, dim(cs)[1])
})

test_that("svinecop rejects discrete-shaped data when var_types is omitted", {
  set.seed(1)
  n <- 200
  lambda <- 3
  d <- 2
  x <- matrix(rpois(n * d, lambda), n, d)
  udisc <- cbind(ppois(x, lambda), ppois(x - 1, lambda))
  cs <- cvine_structure(1:d)
  expect_error(
    svinecop(udisc, p = 1, cs_structure = cs,
             out_vertices = 1:d, in_vertices = 1:d),
    regexp = "wrong number of columns"
  )
})
