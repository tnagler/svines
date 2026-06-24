context("Discrete margins")

set.seed(5)

make_poisson_margin <- function(lambda) {
  fit <- list(p = function(x) ppois(x, lambda),
              q = function(p) qpois(p, lambda),
              p_sub = function(x) ppois(x - 1L, lambda))
  attr(fit, "type") <- "discrete"
  attr(fit, "model") <- "poisson"
  attr(fit, "logLik") <- NA_real_
  attr(fit, "df") <- 1L
  class(fit) <- "svine_margin"
  fit
}

n <- 200
x2 <- matrix(rpois(n * 2, 3), ncol = 2)
xm <- cbind(rnorm(n), rpois(n, 3))
bc <- bicop_dist("gaussian", parameters = 0.5)
mp <- list(make_poisson_margin(3), make_poisson_margin(5))

test_that("creating discrete models", {
  cop <- svinecop_dist(list(list(bc)), cvine_structure(1:2), p = 0,
                       out_vertices = 1:2, in_vertices = 1:2, var_types = c("d", "d"))
  expect_s3_class(cop, "svinecop_dist")

  m <- svine_dist(mp, cop)
  expect_s3_class(m, "svine_dist")
  expect_equal(ncol(to_unif(x2, m$margins)), 4)
  expect_equal(attr(logLik(mp[[1]]), "df"), 1L)

  cs <- cvine_structure(1:2)
  expect_error(svinecop(cbind(ppois(x2, 3), ppois(x2 - 1, 3)), p = 1, cs_structure = cs,
                        out_vertices = 1:2, in_vertices = 1:2), "wrong number of columns")
  sv <- svinecop(matrix(runif(n * 2), ncol = 2), p = 1, cs_structure = cs,
                 out_vertices = 1:2, in_vertices = 1:2)
  expect_length(sv$var_types, 2)
})

test_that("fitting discrete models", {
  fit <- svine(x2, p = 1, var_types = c("d", "d"), margin_families = "pois", family_set = "gaussian")
  expect_s3_class(fit, "svine")
  expect_equal(unname(dim(svine_sim(10, 1, fit))), c(10, 2))

  u <- cbind(ppois(x2, 3), ppois(x2 - 1, 3))
  expect_silent(svinecop(u, p = 1, var_types = c("d", "d"), family_set = "gaussian"))

  expect_s3_class(svine(xm, p = 1, var_types = c("c", "d"),
                        margin_families = list(univariateML::univariateML_models, "pois"),
                        family_set = "gaussian"), "svine")

  df <- data.frame(a = rpois(n, 3), b = factor(sample(letters[1:3], n, replace = TRUE)))
  expect_error(svine(df, p = 1), "factor")
})

test_that("discrete margin dispatch", {
  m_uni <- univariateML::mlpois(rpois(n, 3))
  expect_equal(pmargin_sub(5, m_uni), ppois(4, m_uni[["lambda"]]))
  expect_equal(pmargin_sub(5, mp[[1]]), ppois(4, 3))
  m_cont <- univariateML::mlnorm(rnorm(n))
  expect_equal(pmargin_sub(0.5, m_cont), pmargin(0.5, m_cont))
  expect_equal(pmargin(5, mp[[1]]), ppois(5, 3))
  expect_equal(qmargin(0.5, mp[[1]]), qpois(0.5, 3))

  m_d <- select_margin(rpois(n, 3), c("pois", "norm"), "aic", var_type = "d")
  expect_equal(attr(m_d, "type"), "discrete")
  expect_error(select_margin(rpois(n, 3), "norm", "aic", var_type = "d"), "var_type")
  expect_equal(attr(select_margin(rpois(n, 3), "empirical", "aic", var_type = "d"), "type"), "discrete")

  bad <- mp[[1]]; bad$p_sub <- NULL
  expect_error(check_margin(bad, 1), "p_sub")
})

test_that("discrete scores and hessian", {
  u <- to_unif(x2, mp)
  fit0 <- svinecop(u, p = 0, var_types = c("d", "d"), family_set = "gaussian")
  expect_equal(dim(svinecop_scores(u, fit0)), c(n, fit0$npars))
  expect_true(all(is.finite(svinecop_hessian(u, fit0))))

  fit1 <- svinecop(u, p = 1, var_types = c("d", "d"), family_set = "gaussian")
  expect_equal(dim(svinecop_scores(u, fit1)), c(n - 1, fit1$npars))

  fit <- svine(x2, p = 0, var_types = c("d", "d"), margin_families = "pois", family_set = "gaussian")
  k <- sum(sapply(fit$margins, length)) + fit$copula$npars
  H <- svine_hessian(x2, fit)
  expect_equal(dim(svine_scores(x2, fit)), c(n, k))
  expect_true(all(is.finite(H)))

  eps <- 1e-3
  k_mrg <- sum(sapply(fit$margins, length))
  a <- fit$margins[[1]]; a[1] <- a[1] + eps
  b <- fit$margins[[1]]; b[1] <- b[1] - eps
  sp <- svinecop_scores(to_unif(x2, list(a, fit$margins[[2]])), fit$copula)
  sm <- svinecop_scores(to_unif(x2, list(b, fit$margins[[2]])), fit$copula)
  expect_lt(max(abs(H[1, k_mrg + 1] - mean(sp[, 1] - sm[, 1]) / (2 * eps))), 5e-3)
})

test_that("discrete bootstrap", {
  xi <- rexp(n)
  u_c <- matrix(runif(n * 2), ncol = 2)
  expected <- u_c
  for (j in seq_len(ncol(expected))) {
    s <- sort(expected[, j], index.return = TRUE)
    expected[, j] <- cumsum(xi[s$ix])[order(s$ix)]
    expected[, j] <- expected[, j] / (max(expected[, j]) + 1e-10)
  }
  expect_identical(bootstrap_pobs(u_c, xi), expected)

  u_d <- cbind(ppois(x2[, 1], 3), ppois(x2[, 2], 3), ppois(x2[, 1] - 1, 3), ppois(x2[, 2] - 1, 3))
  ub <- bootstrap_pobs(u_d, xi, var_types = c("d", "d"))
  expect_true(all(ub[, 3] <= ub[, 1]) && all(ub[, 4] <= ub[, 2]))

  u_eq <- matrix(0, n, 4)
  u_eq[, 1] <- sort(runif(n, 0.01, 0.99))
  u_eq[, 2] <- runif(n, 0.01, 0.99)
  u_eq[, 3] <- runif(n, 0.005, u_eq[, 1])
  u_eq[, 4] <- runif(n, 0.005, u_eq[, 2])
  u_eq[30, 3] <- u_eq[10, 1]
  ub2 <- bootstrap_pobs(u_eq, xi, var_types = c("d", "d"))
  expect_equal(ub2[30, 3], ub2[10, 1])
})

test_that("mixed continuous/discrete models", {
  fit <- svine(xm, p = 1, var_types = c("c", "d"),
               margin_families = list(univariateML::univariateML_models, "pois"),
               family_set = "gaussian")
  expect_equal(ncol(to_unif(xm, fit$margins)), 3L)
  expect_true(is.finite(svine_loglik(xm, fit)))
  expect_lt(max(abs(colMeans(svine_scores(xm, fit)))), 0.5)
})
