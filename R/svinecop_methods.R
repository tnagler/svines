#' Simulate from a S-vine copula model
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
#' @param model a S-vine copula model object (inheriting from [svinecop_dist]).
#' @param qrng if `TRUE`, generates quasi-random numbers using the multivariate
#' Generalized Halton sequence up to dimension 300 and the Generalized Sobol
#' sequence in higher dimensions (default `qrng = FALSE`).
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#'
#' @export
#'
svinecop_sim <- function(n, model, qrng = FALSE) {
  assert_that(
    is.number(n),
    inherits(model, "svinecop_dist"),
    is.flag(qrng)
  )
  
  U <- svinecop_sim_cpp(model, n, qrng, rvinecopulib:::get_seeds())
  if (!is.null(model$names)) {
    colnames(U) <- simplify_names(model)
  }
  
  U
}

simplify_names <- function(model) {
  nms <- model$names[seq_along(model$in_vertices)]
  nms <- strsplit(nms, "-")
  sapply(nms, function(n) paste(n[-length(n)], collapse = ""))
}

#' @rdname svinecop_sim
#' @param data time series of past observations (more recent observations are
#'    at the bottom of the matrix.)
#' @export
svinecop_sim_conditional <- function(n, data, model, qrng = FALSE, cores = 1) {
  assert_that(
    is.number(n),
    inherits(model, "svinecop_dist"),
    is.flag(qrng),
    is.number(cores)
  )
  data <- rvinecopulib:::if_vec_to_matrix(data, dim(model$cs_structure)[1] == 1)
  
  U <- svinecop_sim_conditional_cpp(model, n, data, qrng, cores,
                                    rvinecopulib:::get_seeds())
  if (!is.null(model$names)) {
    colnames(U) <- simplify_names(model)
  }
  
  U
}

#' @rdname svinecop_sim
#' @export
svinecop_sim_ahead <- function(n, data, model, qrng = FALSE) {
  assert_that(
    is.number(n),
    inherits(model, "svinecop_dist"),
    is.flag(qrng)
  )
  data <- rvinecopulib:::if_vec_to_matrix(data, dim(model$cs_structure)[1] == 1)
  
  U <- svinecop_sim_ahead_cpp(model, n, data, qrng, rvinecopulib:::get_seeds())
  if (!is.null(model$names)) {
    colnames(U) <- simplify_names(model)
  }
  
  U
}

#' Log-likelihood for S-vine models
#' 
#' @param u the data; should have approximately uniform margins..
#' @param model model inherting from class [svinecop_dist].
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @export
svinecop_loglik <- function(u, model, cores = 1) {
  assert_that(inherits(model, "svinecop_dist"))
  u <- rvinecopulib:::if_vec_to_matrix(u, dim(model$cs_structure)[1] == 1)
  svinecop_loglik_cpp(u, model, cores)
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