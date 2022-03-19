context("S-vine distributions")

set.seed(5)

bc <- bicop_dist("clay", 0, 3)
mrg <- univariateML::mlnorm(rnorm(20))

test_that("creating custom models (1d)", {
  cop <- svinecop_dist(
    list(),
    rvine_structure(1),
    p = 0,
    out_vertices = 1,
    in_vertices = 1
  )
  m <- svine_dist(list(mrg), cop)
  expect_length(svine_sim(10, 1, m), 10)
  
  cop <- svinecop_dist(
    list(list(bc)),
    rvine_structure(1),
    p = 1,
    out_vertices = 1,
    in_vertices = 1
  )
  m <- svine_dist(list(mrg), cop)
  expect_length(svine_sim(10, 1, m), 10)
  
  cop <- svinecop_dist(
    lapply(1:3, function(j) replicate(4 - j, bc, simplify = FALSE)),
    rvine_structure(1),
    p = 3,
    out_vertices = 1,
    in_vertices = 1
  )
  m <- svine_dist(list(mrg), cop)
  expect_length(svine_sim(10, 1, m), 10)
})


test_that("creating custom models (multivariate)", {
  mrgs <- lapply(1:4, function(j) mrg)
  cop <- svinecop_dist(
    lapply(1:3, function(j) replicate(4 - j, bc, simplify = FALSE)),
    dvine_structure(1:4),
    p = 0,
    out_vertices = 4:1,
    in_vertices = 1:4
  )
  m <- svine_dist(mrgs, cop)
  expect_length(svine_sim(10, 1, m), 40)
  
  cop <- svinecop_dist(
    lapply(1:11, function(j) replicate(12 - j, bc, simplify = FALSE)),
    dvine_structure(1:4),
    p = 2,
    out_vertices = 1:4,
    in_vertices = 1:4
  )
  m <- svine_dist(mrgs, cop)
  expect_length(svine_sim(10, 1, m), 40)
})

test_that("fitting models (1d)", {
  x <- rnorm(100)
  fit <- svine(x, p = 0)
  expect_equal(unname(dim(svine_sim(10, 1, fit))), c(10, 1))
  expect_equal(unname(dim(svine_sim(10, 1, fit, x))), c(10, 1))
  expect_equal(unname(dim(svine_sim(10, 3, fit, x))), c(10, 1, 3))
  
  fit <- svine(x, p = 3)
  expect_equal(unname(dim(svine_sim(10, 1, fit))), c(10, 1))
  expect_equal(unname(dim(svine_sim(10, 1, fit, x))), c(10, 1))
  expect_equal(unname(dim(svine_sim(10, 3, fit, x))), c(10, 1, 3))
})


test_that("fitting models (multivariate)", {
  u <- rbicop(50, bc)
  x <- qexp(u)
  fit <- svine(x, margin_families = c("exp", "norm"), p = 0)
  AIC(fit)
  expect_equal(unname(dim(svine_sim(10, 1, fit))), c(10, 2))
  expect_equal(unname(dim(svine_sim(10, 1, fit, x))), c(10, 2))
  expect_equal(unname(dim(svine_sim(10, 3, fit, x))), c(10, 2, 3))
  expect_length(svine_loglik(x, fit), 1)
  
  fit <- svine(x, margin_families = c("exp", "norm"), p = 3)
  AIC(fit)
  expect_gt(min(svine_sim(10, 1, fit)), 0)
  expect_gt(min(svine_sim(10, 1, fit, x)), 0)
  expect_gt(min(svine_sim(10, 5, fit, x)), 0)
  
  fit2 <- svine(x,
                p = 3, 
                margin_families = c("exp", "norm"),
                cs_structure = fit$copula$cs_structure,
                out_vertices = fit$copula$out_vertices,
                in_vertices = fit$copula$in_vertices)
  expect_equal(fit, fit2, tolerance = 0.01)
  
  expect_silent(svine(x, p = 1, type = "D", margin_families = c("exp", "norm")))
  expect_silent(svine(x, p = 2, type = "M", margin_families = c("exp", "norm")))
  expect_error(svine(x, p = 1, type = "R", margin_families = c("exp", "norm")))
})

test_that("standard errors", {
  u <- rbicop(50, bc)
  x <- qexp(u)
  fit <- svine(x, margin_families = c("std", "norm", "exp"), p = 0)
  models <- svine_bootstrap_models(2, fit)
  expect_identical(summary(models[[1]])$family, summary(fit)$family)

  fit <- svine(x, margin_families = "empirical", p = 0)
  models <- svine_bootstrap_models(2, fit)
  expect_identical(summary(models[[1]])$family, summary(fit)$family)
})

test_that("print/summary", {
  u <- rbicop(50, bc)
  x <- qexp(u)
  colnames(x) <- letters[1:2]
  fit <- svine(x, margin_families = c("exp", "norm"), p = 0)
  expect_output(print(fit))
  expect_silent(sm <- summary(fit))
  expect_length(sm, 2)
  expect_s3_class(sm[[1]], "summary_df")
  expect_s3_class(sm[[2]], "summary_df")
  expect_equal(nrow(sm[[1]]), 2)
  expect_equal(nrow(sm[[2]]), 1)
})
