

#' Simulate from a S-vine model
#'
#' [svinecop_sim()] simulates a time series of length `n`;
#' [svinecop_sim_conditional()] simulates `n` observations for the next time
#' point conditional on past observations;
#' [svinecop_sim_ahead()]
#' simulates the next `n` time points conditional on the past.
#'
#' @aliases svinecop_sim svinecop_sim_conditional svinecop_sim_ahead
#'
#' @param n number of observations.
#' @param svinecop a S-vine model object.
#' @param qrng if `TRUE`, generates quasi-random numbers using the multivariate
#' Generalized Halton sequence up to dimension 300 and the Generalized Sobol
#' sequence in higher dimensions (default `qrng = FALSE`).
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#'
#' @export
#'
svinecop_sim <- function(n, svinecop, qrng = FALSE) {
  assert_that(
    is.number(n),
    inherits(svinecop, "svinecop_dist"),
    is.flag(qrng)
  )
  
  U <- svinecop_sim_cpp(svinecop, n, qrng, rvinecopulib:::get_seeds())
  if (!is.null(svinecop$names)) {
    colnames(U) <- simplify_names(svinecop)
  }
  
  U
}

simplify_names <- function(svinecop) {
  nms <- svinecop$names[seq_along(svinecop$in_vertices)]
  nms <- strsplit(nms, "-")
  sapply(nms, function(n) paste(n[-length(n)], collapse = ""))
}

#' @rdname svinecop_sim
#' @param data time series of past observations (more recent observations are
#'    at the bottom of the matrix.)
#' @export
svinecop_sim_conditional <- function(n, data, svinecop, qrng = FALSE, cores = 1) {
  assert_that(
    is.number(n),
    inherits(svinecop, "svinecop_dist"),
    is.flag(qrng),
    is.number(cores)
  )
  data <- rvinecopulib:::if_vec_to_matrix(data,
                                          dim(svinecop$cs_structure)[1] == 1)
  
  U <- svinecop_sim_conditional_cpp(
    svinecop, n, data, qrng, cores,
    rvinecopulib:::get_seeds()
  )
  if (!is.null(svinecop$names)) {
    colnames(U) <- simplify_names(svinecop)
  }
  
  U
}

#' @rdname svinecop_sim
#' @export
svinecop_sim_ahead <- function(n, data, svinecop, qrng = FALSE) {
  assert_that(
    is.number(n),
    inherits(svinecop, "svinecop_dist"),
    is.flag(qrng)
  )
  data <- rvinecopulib:::if_vec_to_matrix(data, dim(svinecop$cs_structure)[1] == 1)
  
  U <- svinecop_sim_ahead_cpp(
    svinecop, n, data, qrng,
    rvinecopulib:::get_seeds()
  )
  if (!is.null(svinecop$names)) {
    colnames(U) <- simplify_names(svinecop)
  }
  
  U
}

#' @export
svinecop_loglik <- function(u, svinecop, cores = 1) {
  assert_that(inherits(svinecop, "svinecop_dist"))
  u <- rvinecopulib:::if_vec_to_matrix(u, dim(svinecop$cs_structure)[1] == 1)
  svinecop_loglik_cpp(u, svinecop, cores)
}

# #' @export
# svinecop_scores <- function(u, svinecop, cores = 1) {
#   assert_that(inherits(svinecop, "svinecop_dist"))
#   u <- rvinecopulib:::if_vec_to_matrix(u, dim(svinecop$cs_structure)[1] == 1)
#   svinecop_scores_cpp(u, svinecop, cores)
# }

# #' @export
# svinecop_hessian <- function(u, svinecop, cores = 1) {
#   assert_that(inherits(svinecop, "svinecop_dist"))
#   u <- rvinecopulib:::if_vec_to_matrix(u, dim(svinecop$cs_structure)[1] == 1)
#   svinecop_hessian_cpp(u, svinecop, cores)
# }

# svinecop_cond_cdf <- function(u, conditioned, svinecop, cores = 1) {
#   assert_that(
#     is.count(conditioned),
#     inherits(svinecop, "svinecop_dist"),
#     conditioned <= dim(svinecop$cs_structure)[1]
#   )
#   u <- rvinecopulib:::if_vec_to_matrix(u, dim(svinecop$cs_structure)[1] == 1)
#   svinecop_cond_cdf_cpp(u, conditioned - 1, svinecop, cores)
# }