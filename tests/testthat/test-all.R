bc <- bicop_dist("clay", 0, 3)

test_that("creating custom models (1d)", {
  tv <- svinecop_dist(
    list(),
    rvine_structure(1),
    p = 0,
    in_vertices = 1,
    out_vertices = 1
  )
  expect_equal(unname(dim(tv)), c(1, 0))

  tv <- svinecop_dist(
    list(list(bc)),
    rvine_structure(1),
    p = 1,
    in_vertices = 1,
    out_vertices = 1
  )
  expect_equal(unname(dim(tv)), c(2, 1))

  tv <- svinecop_dist(
    lapply(1:3, function(j) replicate(4 - j, bc, simplify = FALSE)),
    rvine_structure(1),
    p = 3,
    in_vertices = 1,
    out_vertices = 1
  )
  expect_equal(unname(dim(tv)), c(4, 3))
})


test_that("creating custom models (multivariate)", {
  tv <- svinecop_dist(
    lapply(1:3, function(j) replicate(4 - j, bc, simplify = FALSE)),
    dvine_structure(1:4),
    p = 0,
    in_vertices = 1:4,
    out_vertices = 4:1
  )
  expect_equal(unname(dim(tv)), c(4, 3))

  tv <- svinecop_dist(
    lapply(1:11, function(j) replicate(12 - j, bc, simplify = FALSE)),
    dvine_structure(1:4),
    p = 2,
    in_vertices = 1:4,
    out_vertices = 1:4
  )
  expect_equal(unname(dim(tv)), c(12, 11))

})

test_that("fitting models (1d)", {
  u <- runif(100)
  tv <- svinecop(u, p = 0)
  expect_length(svinecop_sim(10, tv), 10)
  expect_length(svinecop_sim_ahead(10, u, tv), 10)
  expect_length(svinecop_sim_conditional(10, u, tv), 10)

  tv <- svinecop(u, p = 3)
  expect_length(svinecop_sim(10, tv), 10)
  expect_length(svinecop_sim_ahead(10, u, tv), 10)
  expect_length(svinecop_sim_conditional(10, u, tv), 10)
})


test_that("fitting models (multivariate)", {
  u <- rbicop(50, bc)
  tv <- svinecop(u, p = 0)
  AIC(tv)
  expect_equal(dim(svinecop_sim(10, tv)), c(10, 2))
  expect_equal(dim(svinecop_sim_ahead(10, u, tv)), c(10, 2))
  expect_equal(dim(svinecop_sim_conditional(10, u, tv)), c(10, 2))

  tv <- svinecop(u, p = 3)
  AIC(tv)
  expect_equal(dim(svinecop_sim(10, tv)), c(10, 2))
  expect_equal(dim(svinecop_sim_ahead(10, u, tv)), c(10, 2))
  expect_equal(dim(svinecop_sim_conditional(10, u, tv)), c(10, 2))

  tv2 <- svinecop(u, p = 3,
                  cs_structure = tv$cs_structure,
                  in_vertices = tv$in_vertices,
                  out_vertices = tv$out_vertices)
  expect_equal(tv, tv2)

  expect_silent(svinecop(u, p = 1, type = "D"))
  expect_silent(svinecop(u, p = 2, type = "M"))
  expect_error(svinecop(u, p = 1, type = "R"))

})
