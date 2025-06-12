context("S-vine copulas")

set.seed(5)

bc <- bicop_dist("clay", 0, 3)

test_that("creating custom models (1d)", {
  sv <- svinecop_dist(
    list(),
    rvine_structure(1),
    p = 0,
    out_vertices = 1,
    in_vertices = 1
  )
  expect_equal(unname(dim(sv)), c(1, 0))

  sv <- svinecop_dist(
    list(list(bc)),
    rvine_structure(1),
    p = 1,
    out_vertices = 1,
    in_vertices = 1
  )
  expect_equal(unname(dim(sv)), c(2, 1))

  sv <- svinecop_dist(
    lapply(1:3, function(j) replicate(4 - j, bc, simplify = FALSE)),
    rvine_structure(1),
    p = 3,
    out_vertices = 1,
    in_vertices = 1
  )
  expect_equal(unname(dim(sv)), c(4, 3))
})


test_that("creating custom models (multivariate)", {
  sv <- svinecop_dist(
    lapply(1:3, function(j) replicate(4 - j, bc, simplify = FALSE)),
    dvine_structure(1:4),
    p = 0,
    out_vertices = 4:1,
    in_vertices = 1:4
  )
  expect_equal(unname(dim(sv)), c(4, 3))

  sv <- svinecop_dist(
    lapply(1:11, function(j) replicate(12 - j, bc, simplify = FALSE)),
    dvine_structure(1:4),
    p = 2,
    out_vertices = 1:4,
    in_vertices = 1:4
  )
  expect_equal(unname(dim(sv)), c(12, 11))

})

test_that("fitting models (1d)", {
  u <- runif(100)
  sv <- svinecop(u, p = 0)
  expect_equal(unname(dim(svinecop_sim(10, 1, sv))), c(10, 1))
  expect_equal(unname(dim(svinecop_sim(10, 1, sv, u))), c(10, 1))
  expect_equal(unname(dim(svinecop_sim(10, 3, sv, u))), c(10, 1, 3))
  
  sv <- svinecop(u, p = 3)
  expect_equal(unname(dim(svinecop_sim(10, 1, sv))), c(10, 1))
  expect_equal(unname(dim(svinecop_sim(10, 1, sv, u))), c(10, 1))
  expect_equal(unname(dim(svinecop_sim(10, 3, sv, u))), c(10, 1, 3))
})


test_that("fitting models (multivariate)", {
  set.seed(1)
  u <- rbicop(50, bc)
  sv <- svinecop(u, p = 0)
  AIC(sv)
  expect_equal(unname(dim(svinecop_sim(10, 1, sv))), c(10, 2))
  expect_equal(unname(dim(svinecop_sim(10, 1, sv, u))), c(10, 2))
  expect_equal(unname(dim(svinecop_sim(10, 3, sv, u))), c(10, 2, 3))
  
  sv <- svinecop(u, p = 3)
  AIC(sv)
  expect_equal(unname(dim(svinecop_sim(10, 1, sv))), c(10, 2))
  expect_equal(unname(dim(svinecop_sim(10, 1, sv, u))), c(10, 2))
  expect_equal(unname(dim(svinecop_sim(10, 3, sv, u))), c(10, 2, 3))
  
  sv2 <- svinecop(u, p = 3,
                  cs_structure = sv$cs_structure,
                  out_vertices = sv$out_vertices,
                  in_vertices = sv$in_vertices)
  expect_equal(sv, sv2, tolerance = 0.01)
  
  expect_silent(svinecop(u, p = 1, type = "D"))
  expect_silent(svinecop(u, p = 2, type = "M"))
  expect_error(svinecop(u, p = 1, type = "R"))
  expect_equal(svinecop_loglik(u, sv), as.numeric(logLik(sv)), 0.1)
  
  # just check whether sim works with only one sample
  expect_silent(svinecop_sim(1, 1, sv))
  expect_silent(svinecop_sim(1, 1, sv, u))
})

