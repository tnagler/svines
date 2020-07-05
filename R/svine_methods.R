#' Simulate from a S-vine model
#'
#' [svine_sim()] simulates a time series of length `n`;
#' [svine_sim_conditional()] simulates `n` observations for the next time
#' point conditional on past observations;
#' [svine_sim_ahead()]
#' simulates the next `n` time points conditional on the past.
#'
#' @aliases svinecop_sim svinecop_sim_conditional svinecop_sim_ahead
#'
#' @param n number of observations.
#' @param svine a S-vine model.
#' @param qrng if `TRUE`, generates quasi-random numbers using the multivariate
#' Generalized Halton sequence up to dimension 300 and the Generalized Sobol
#' sequence in higher dimensions (default `qrng = FALSE`).
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#'
#' @export
#'
svine_sim <- function(n, svine, qrng = FALSE) {
  assert_that(
    is.number(n),
    inherits(svine, "svine_dist"),
    is.flag(qrng)
  )
  
  U <- svinecop_sim_cpp(svine$copula, n, qrng, rvinecopulib:::get_seeds())
  to_quantiles(U, svine$margins)
}


#' @rdname svine_sim
#' @param data time series of past observations (more recent observations are
#'    at the bottom of the matrix.)
#' @export
svine_sim_conditional <- function(n, data, svine, qrng = FALSE, cores = 1) {
  assert_that(
    is.number(n),
    inherits(svine, "svine_dist"),
    is.flag(qrng),
    is.number(cores)
  )
  data <- rvinecopulib:::if_vec_to_matrix(data, length(svine$margins) == 1)
  data <- to_unif(data, svine$margins)
  
  U <- svinecop_sim_conditional_cpp(
    svine$copula, n, data, qrng, cores,
    rvinecopulib:::get_seeds()
  )
  to_quantiles(U, svine$margins)
}


#' @rdname svine_sim
#' @export
svine_sim_ahead <- function(n, data, svine, qrng = FALSE) {
  assert_that(
    is.number(n),
    inherits(svine, "svine_dist"),
    is.flag(qrng)
  )
  data <- rvinecopulib:::if_vec_to_matrix(data, length(svine$margins) == 1)
  data <- to_unif(data, svine$margins)
  
  U <- svinecop_sim_ahead_cpp(
    svine$copula, n, data, qrng,
    rvinecopulib:::get_seeds()
  )
  to_quantiles(U, svine$margins)
}

#' @export
svine_loglik <- function(x, svine, cores = 1) {
  assert_that(inherits(svine, "svine_dist"))
  x <- rvinecopulib:::if_vec_to_matrix(x, length(svine$margins) == 1)
  asser_that(ncol(x) == length(svine$margins))
  
  ll_marg <- sapply(
    seq_len(ncol(x)), 
    function(j) sum(dml(x[, j], svine$margins[[j]], log = TRUE))
  )
  u <- to_unif(x, svine$margins)
  ll_cop <- svinecop_loglik_cpp(u, svine$copula, cores)
  ll_marg + ll_cop
}

# #' @export
# svine_scores <- function(x, svine, cores = 1) {
#   assert_that(inherits(svine, "svine_dist"))
#   x <- rvinecopulib:::if_vec_to_matrix(x, length(svine$margins) == 1)
#   u <- to_unif(x, svine$margins)
#   svinecop_scores_cpp(u, svine$copula, cores)
# }

# #' @export
# svine_hessian <- function(x, svine, cores = 1) {
#   assert_that(inherits(svine, "svine_dist"))
#   x <- rvinecopulib:::if_vec_to_matrix(x, length(svine$margins) == 1)
#   u <- to_unif(x, svine$margins)
#   svinecop_hessian_cpp(u, svine$copula, cores)
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

#' @export
logLik.svine <- function(object, ...) {
  ll <- -2 * object$loglik
  attr(ll, "df") <- object$npars
  ll
}

