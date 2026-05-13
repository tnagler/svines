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

test_that("check_margin rejects discrete margin lacking $p_sub", {
  bad <- list(
    p = function(x) ppois(x, 3),
    q = function(p) qpois(p, 3)
    # intentionally omits $p_sub
  )
  attr(bad, "type")  <- "discrete"
  attr(bad, "model") <- "poisson_bad"
  class(bad) <- "svine_margin"

  expect_error(check_margin(bad, 1), regexp = "p_sub")
})

test_that("pmargin_sub routes raw univariateML-discrete margins through pml(x - 1)", {
  set.seed(11)
  lambda <- 3
  m <- univariateML::mlpois(rpois(100, lambda))

  # Raw mlpois object — never went through select_margin, so
  # attr(., "type") is unset. The robustness expansion in the first
  # arm keys on isFALSE(attr(., "continuous")) instead, which works
  # for both auto-fit and user-constructed mlpois/mlnbinom/mlgeom.
  expect_equal(pmargin_sub(5, m), ppois(4, m[["lambda"]]))
})

test_that("pmargin_sub routes list-form discrete margins through $p_sub", {
  set.seed(11)
  lambda <- 3
  m <- make_poisson_margin(lambda)

  expect_equal(pmargin_sub(5, m), ppois(4, lambda))
})

test_that("pmargin_sub falls through to pmargin for continuous margins", {
  set.seed(11)
  m <- univariateML::mlnorm(rnorm(100))

  # Continuous margins satisfy F(x-) = F(x); the third arm delegates
  # to pmargin without shifting the input.
  expect_equal(pmargin_sub(0.5, m), pmargin(0.5, m))
})

test_that("pmargin and qmargin dispatch on inherits(model, 'univariateML')", {
  set.seed(11)

  # univariateML route delegates to pml / qml. Compare against base-R
  # pnorm / qnorm rather than pml / qml so the assertion verifies the
  # numeric output, not just that two equivalent dispatches agree.
  m_uni <- univariateML::mlnorm(rnorm(100))
  expect_equal(pmargin(0.5, m_uni), pnorm(0.5, m_uni[["mean"]], m_uni[["sd"]]))
  expect_equal(qmargin(0.5, m_uni), qnorm(0.5, m_uni[["mean"]], m_uni[["sd"]]))

  # List-form route invokes the user-supplied $p / $q callables.
  m_list <- make_poisson_margin(3)
  expect_equal(pmargin(5, m_list), ppois(5, 3))
  expect_equal(qmargin(0.5, m_list), qpois(0.5, 3))
})

test_that("select_margin tags univariateML-discrete fits with attr(., 'type') == 'discrete'", {
  set.seed(11)

  # Discrete fit: the line 89 tag derives from univariateML's
  # `continuous` flag (FALSE -> "discrete"). This is the dispatch
  # key downstream code reads (to_unif, logLik.svine_margin).
  m_d <- select_margin(rpois(200, 3), "pois", "aic")
  expect_equal(attr(m_d, "type"), "discrete")

  # Continuous fit: falls through to the input-arm `type` variable,
  # which is "univariateML" for all non-empirical family sets.
  m_c <- select_margin(rnorm(200), "norm", "aic")
  expect_equal(attr(m_c, "type"), "univariateML")
})

test_that("svine rejects factor columns with an informative error", {
  set.seed(13)
  n <- 50
  df <- data.frame(
    x = rpois(n, 3),
    y = factor(sample(letters[1:3], n, replace = TRUE))
  )
  expect_error(
    svine(df, p = 1),
    regexp = "factor"
  )
})

test_that("svine does not reject an integer matrix at the factor gate", {
  set.seed(13)
  n <- 100
  x <- matrix(rpois(n * 2, 3), n, 2)

  # Narrow gate test: we're verifying svine.R:40-43's factor check
  # doesn't fire on an integer matrix, not that the full discrete
  # pipeline works (Phase 2.7 var_types threading is not yet in).
  # If the call errors for any reason, the assertion is that the
  # error didn't come from the factor gate.
  err <- tryCatch(
    svine(x, p = 1),
    error = identity
  )
  if (inherits(err, "error")) {
    expect_false(grepl("factor", conditionMessage(err)))
  } else {
    expect_true(TRUE)
  }
})
