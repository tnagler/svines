context("Discrete margins")

set.seed(5)

# make_poisson_margin builds a discrete svine_margin object that
# satisfies the full β margin contract:
#   $p     — F(x), the CDF
#   $q     — F^{-1}(p), the quantile function
#   $p_sub — F(x-), the CDF at the largest support point strictly
#            less than x; for integer-valued margins this is
#            F(x - 1), hence ppois(x - 1, lambda)
# Plus the four svine_margin attributes that downstream code reads:
#   type    — "discrete", the dispatch key for pmargin / qmargin /
#             pmargin_sub / logLik.svine_margin
#   model   — printable name surfaced by summary() methods
#   logLik  — model log-likelihood at fitted data; NA here because
#             the margin is hand-constructed, not fitted
#   df      — parameter count (Poisson has one parameter, lambda)
# class includes "svine_margin" so the S3 logLik method dispatches.
make_poisson_margin <- function(lambda) {
  fit <- list(
    p     = function(x) ppois(x, lambda),
    q     = function(p) qpois(p, lambda),
    p_sub = function(x) ppois(x - 1, lambda)
  )
  attr(fit, "type")   <- "discrete"
  attr(fit, "model")  <- "poisson"
  attr(fit, "logLik") <- NA_real_
  attr(fit, "df")     <- 1L
  class(fit) <- "svine_margin"
  fit
}

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

test_that("svine_dist builds end-to-end with hand-constructed Poisson margins", {
  set.seed(7)
  d <- 2

  margins <- list(make_poisson_margin(3), make_poisson_margin(5))
  bc <- bicop_dist("gaussian", parameters = 0.4)
  cop <- svinecop_dist(
    pair_copulas = list(list(bc)),
    cs_structure = cvine_structure(1:d),
    p            = 0,
    out_vertices = 1:d,
    in_vertices  = 1:d,
    var_types    = c("d", "d")
  )
  m <- svine_dist(margins = margins, copula = cop)
  expect_s3_class(m, "svine_dist")

  # to_unif on integer-valued data must produce the doubled-column
  # n x 2d layout: F(x) in columns 1..d, F(x-) in columns d+1..2d.
  n <- 50
  x <- cbind(rpois(n, 3), rpois(n, 5))
  u <- to_unif(x, m$margins)
  expect_equal(ncol(u), 2 * d)

  # logLik on a discrete margin dispatches through the discrete arm
  # of logLik.svine_margin and surfaces the df attribute (1 for
  # Poisson). logLik on the svine_dist object itself is not
  # supported pre-α (no logLik.svine_dist method exists; the
  # fitted-svine path defines logLik.svine separately).
  expect_equal(attr(logLik(margins[[1]]), "df"), 1L)
})
