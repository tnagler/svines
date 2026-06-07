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
  n_disc <- sum(vapply(m$margins, function(mm) identical(attr(mm, "type"), "discrete"), logical(1)))
  expect_equal(ncol(u), d + n_disc)

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
  m_d <- select_margin(rpois(200, 3), "pois", "aic", var_type = "d")
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

test_that("select_margin restricts to the discrete panel under var_type = 'd'", {
  set.seed(17)
  m <- select_margin(rpois(200, 3), c("pois", "norm"), "aic", var_type = "d")
  expect_equal(attr(m, "model"), "Poisson")
  expect_equal(attr(m, "type"), "discrete")
})

test_that("select_margin errors when families and var_type produce an empty intersection", {
  set.seed(17)
  expect_error(
    select_margin(rpois(50, 3), "norm", "aic", var_type = "d"),
    regexp = "var_type"
  )
})

test_that("empirical margins are supported on discrete data end-to-end", {
  set.seed(17)
  n <- 200
  x <- rpois(n, 3)

  m <- select_margin(x, "empirical", "aic", var_type = "d")
  expect_equal(attr(m, "type"), "discrete")
  expect_true(is.function(m$p_sub))

  # F(x-) for an empirical discrete margin is the ecdf one integer lower,
  # on the same (0, 1)-scaled convention pmargin uses.
  F_n <- ecdf(x)
  expect_equal(pmargin_sub(5L, m), F_n(4L) * n / (n + 1))

  # full svine_dist (margins + copula) with two empirical discrete margins
  m2 <- select_margin(rpois(n, 5), "empirical", "aic", var_type = "d")
  cop <- svinecop_dist(
    pair_copulas = list(list(bicop_dist("gaussian", parameters = 0.4))),
    cs_structure = cvine_structure(1:2),
    p            = 0,
    out_vertices = 1:2,
    in_vertices  = 1:2,
    var_types    = c("d", "d")
  )
  expect_no_error(svine_dist(margins = list(m, m2), copula = cop))
})

test_that("svine accepts var_types = c('d', 'd') end-to-end on integer data", {
  set.seed(17)
  n <- 100
  x <- matrix(rpois(n * 2, 3), n, 2)
  fit <- svine(x, p = 1, var_types = c("d", "d"),
               margin_families = "pois", family_set = "gaussian")
  expect_s3_class(fit, "svine")
  expect_equal(attr(fit$margins[[1]], "type"), "discrete")
  expect_equal(attr(fit$margins[[2]], "type"), "discrete")
})

test_that("svine errors when var_types length does not match NCOL(data)", {
  set.seed(17)
  x <- matrix(rpois(100 * 2, 3), 100, 2)
  expect_error(svine(x, p = 1, var_types = c("d", "d", "d")))
})

test_that("svinecop_scores and svinecop_hessian return well-shaped finite values on discrete data (p = 0)", {
  set.seed(101)
  n <- 200
  d <- 2

  margins <- list(make_poisson_margin(3), make_poisson_margin(5))
  x <- cbind(rpois(n, 3), rpois(n, 5))
  u <- to_unif(x, margins)
  fit <- svinecop(u, p = 0, var_types = c("d", "d"), family_set = "gaussian")

  S <- svinecop_scores(u, fit)
  H <- svinecop_hessian(u, fit)
  k_cop <- fit$npars

  expect_equal(dim(S), c(n, k_cop))
  expect_equal(dim(H), c(k_cop, k_cop))
  expect_true(all(is.finite(S)))
  expect_true(all(is.finite(H)))

  # Score-at-MLE smoke probe: at the maximum-likelihood fit, the mean of
  # each score column should be near zero by construction. Tolerance is
  # loose (0.5) because n = 200 is finite and the fit is family-restricted
  # to Gaussian; this is a sanity bound, not an asymptotic-optimality test.
  expect_lt(mean(abs(colMeans(S))), 0.5)
})

test_that("svinecop_scores and svinecop_hessian return well-shaped finite values on discrete data (p = 1)", {
  set.seed(102)
  n <- 200
  d <- 2

  margins <- list(make_poisson_margin(3), make_poisson_margin(5))
  x <- cbind(rpois(n, 3), rpois(n, 5))
  u <- to_unif(x, margins)
  fit <- svinecop(u, p = 1, var_types = c("d", "d"), family_set = "gaussian")

  S <- svinecop_scores(u, fit)
  H <- svinecop_hessian(u, fit)
  k_cop <- fit$npars

  # spread_lag drops one row per Markov lag in the C++ scores path, so the
  # output row count is n - p (here 200 - 1 = 199).
  expect_equal(dim(S), c(n - 1, k_cop))
  expect_equal(dim(H), c(k_cop, k_cop))
  expect_true(all(is.finite(S)))
  expect_true(all(is.finite(H)))

  expect_lt(mean(abs(colMeans(S))), 0.5)
})

test_that("svine_scores returns well-shaped finite values on discrete data (p = 0)", {
  set.seed(201)
  n <- 200

  x <- cbind(rpois(n, 3), rpois(n, 5))
  fit <- svine(x, p = 0, var_types = c("d", "d"),
               margin_families = "pois", family_set = "gaussian")
  S <- svine_scores(x, fit)
  k_total <- sum(sapply(fit$margins, length)) + fit$copula$npars

  expect_equal(dim(S), c(n, k_total))
  expect_true(all(is.finite(S)))

  # Score-at-MLE smoke probe. Observed value is ~2.4e-4 on this seed; the
  # 5e-3 bound leaves ample headroom (~20x) while still catching the O(0.01)
  # bias a buggy mixed-block computation would produce on the copula-score
  # column.
  expect_lt(mean(abs(colMeans(S))), 5e-3)
})

test_that("svine_hessian returns well-shaped finite values on discrete data (p = 0)", {
  set.seed(202)
  n <- 200

  x <- cbind(rpois(n, 3), rpois(n, 5))
  fit <- svine(x, p = 0, var_types = c("d", "d"),
               margin_families = "pois", family_set = "gaussian")
  H <- svine_hessian(x, fit)
  k_total <- sum(sapply(fit$margins, length)) + fit$copula$npars

  # Structural-only: dim and finite. The Hessian assembled by svine_hessian
  # places a zero matrix in the bottom-left block (only the top-right mixed
  # block is computed via hessian_mxd), so the full H is asymmetric by
  # construction and isSymmetric() is not a meaningful check here.
  expect_equal(dim(H), c(k_total, k_total))
  expect_true(all(is.finite(H)))
})

test_that("svine_hessian mixed block matches external finite differences on discrete data", {
  # Seed chosen empirically: the magnitude of the shadow-column-omission
  # bias is data-dependent (varies an order of magnitude across seeds),
  # and seed 205 produces a ~1.2e-2 bias here — well above the tolerance
  # below, giving the regression test a reliable signal.
  set.seed(205)
  n <- 200

  x <- cbind(rpois(n, 3), rpois(n, 5))
  fit <- svine(x, p = 0, var_types = c("d", "d"),
               margin_families = "pois", family_set = "gaussian")
  H <- svine_hessian(x, fit)
  k_mrg <- sum(sapply(fit$margins, length))

  # Reproduce hessian_mxd's perturbation externally for (margin 1 param 1,
  # copula param 1). eps matches hessian_mxd's internal 1e-3 so the
  # comparison sees the same parameter-space neighbourhood; mismatched eps
  # introduces second-order disagreement that obscures the bug signal.
  eps <- 1e-3
  margin_1_plus  <- fit$margins[[1]]
  margin_1_plus[1]  <- fit$margins[[1]][1] + eps
  margin_1_minus <- fit$margins[[1]]
  margin_1_minus[1] <- fit$margins[[1]][1] - eps
  margins_plus  <- list(margin_1_plus,  fit$margins[[2]])
  margins_minus <- list(margin_1_minus, fit$margins[[2]])
  u_plus  <- to_unif(x, margins_plus)
  u_minus <- to_unif(x, margins_minus)
  score_plus  <- svinecop_scores(u_plus,  fit$copula)
  score_minus <- svinecop_scores(u_minus, fit$copula)

  # Column 1 of svinecop_scores corresponds to the first copula parameter
  # (canonical order: outer loop over trees, inner loop over edges; here
  # one tree with one edge under d=2, p=0).
  estimate <- mean(score_plus[, 1] - score_minus[, 1]) / (2 * eps)

  # Hessian layout from svine_hessian: rows/cols 1..k_mrg are marginal
  # params, k_mrg+1..end are copula params. H[1, k_mrg + 1] is the
  # (margin 1 param 1, copula param 1) mixed-block entry.
  hessian_entry <- H[1, k_mrg + 1]

  # Tolerance: 5e-3. After the shadow-column update is included in
  # hessian_mxd, the internal perturbation matches the external one
  # (both route through pml(x, .) for F(x) and pml(x-1, .) for F(x-)),
  # so the agreement is at floating-point precision (~1e-12). The 5e-3
  # bound leaves three orders of magnitude of headroom while still
  # failing reliably on the pre-fix O(1e-2) bias.
  expect_lt(abs(hessian_entry - estimate), 5e-3)
})

test_that("bootstrap_pobs reproduces original behaviour on continuous-only input (no var_types)", {
  set.seed(401)
  n <- 50
  d <- 2
  u <- matrix(runif(n * d), n, d)
  xi <- rexp(n)

  # Manually replicate the pre-fix continuous-only algorithm: when
  # var_types is NULL, the function must reduce to its prior behaviour.
  expected <- u
  for (j in seq_len(ncol(expected))) {
    s <- sort(expected[, j], index.return = TRUE)
    expected[, j] <- cumsum(xi[s$ix])[order(s$ix)]
    expected[, j] <- expected[, j] / (max(expected[, j]) + 1e-10)
  }

  expect_identical(bootstrap_pobs(u, xi), expected)
})

test_that("bootstrap_pobs preserves F(x-) <= F(x) pairing under discrete input", {
  set.seed(402)
  n <- 100
  lambda <- c(3, 5)
  x <- cbind(rpois(n, lambda[1]), rpois(n, lambda[2]))
  u <- cbind(
    ppois(x[, 1], lambda[1]),
    ppois(x[, 2], lambda[2]),
    ppois(x[, 1] - 1, lambda[1]),
    ppois(x[, 2] - 1, lambda[2])
  )
  xi <- rexp(n)

  ub <- bootstrap_pobs(u, xi, var_types = c("d", "d"))

  # F_tilde_j is non-decreasing as an empirical-CDF estimate, so applying
  # it to both halves of the (real, shadow) pair preserves the
  # inequality u[, j + d] <= u[, j] that holds by construction for
  # F(x-) <= F(x) on integer-supported margins.
  expect_true(all(ub[, 3] <= ub[, 1]))
  expect_true(all(ub[, 4] <= ub[, 2]))
  expect_true(all(ub >= 0 & ub <= 1))
})

test_that("bootstrap_pobs applies the same F_tilde_j to both halves of each discrete pair", {
  set.seed(403)
  n <- 50
  # Construct u so that u[i_test, 3] (shadow margin 1) exactly equals
  # u[i_other, 1] (real margin 1) at a different row. If the fix applies
  # a single bootstrapped step function F_tilde_1 to both halves, then
  # equal inputs must produce equal outputs across the halves.
  u <- matrix(0, n, 4)
  u[, 1] <- sort(runif(n, 0.01, 0.99))
  u[, 2] <- runif(n, 0.01, 0.99)
  u[, 3] <- runif(n, 0.005, u[, 1])
  u[, 4] <- runif(n, 0.005, u[, 2])

  i_other <- 10L
  i_test  <- 30L
  u[i_test, 3] <- u[i_other, 1]

  xi <- rexp(n)
  ub <- bootstrap_pobs(u, xi, var_types = c("d", "d"))

  expect_equal(ub[i_test, 3], ub[i_other, 1])
})

test_that("svine fits a mixed continuous/discrete model at p = 1", {
  set.seed(42)
  n <- 200
  x <- cbind(rnorm(n), rpois(n, 3))
  fit <- svine(x, p = 1, var_types = c("c", "d"),
               margin_families = list(univariateML::univariateML_models, "pois"),
               family_set = "gaussian")
  expect_s3_class(fit, "svine")
  expect_equal(ncol(to_unif(x, fit$margins)), 3L)
})

test_that("svine loglik and scores are numerically correct for mixed c/d model", {
  set.seed(42)
  n <- 200
  x <- cbind(rnorm(n), rpois(n, 3))
  fit <- svine(x, p = 1, var_types = c("c", "d"),
               margin_families = list(univariateML::univariateML_models, "pois"),
               family_set = "gaussian")

  # log-likelihood should be finite
  ll <- svine_loglik(x, fit)
  expect_true(is.finite(ll))

  # scores at the MLE should be near zero (mean score test)
  sc <- svine_scores(x, fit)
  expect_true(is.matrix(sc))
  expect_true(all(is.finite(sc)))
  expect_lt(max(abs(colMeans(sc))), 0.5)
})
