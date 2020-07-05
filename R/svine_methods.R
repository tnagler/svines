#' Simulate from a S-vine model
#'
#' [svine_sim()] simulates a time series of length `n`;
#' [svine_sim_conditional()] simulates `n` observations for the next time
#' point conditional on past observations;
#' [svine_sim_ahead()]
#' simulates the next `n` time points conditional on the past.
#'
#' @aliases svine_sim svine_sim_conditional svine_sim_ahead
#'
#' @param n number of observations.
#' @param model a S-vine model (inheriting from [svine_dist]).
#' @param qrng if `TRUE`, generates quasi-random numbers using the multivariate
#' Generalized Halton sequence up to dimension 300 and the Generalized Sobol
#' sequence in higher dimensions (default `qrng = FALSE`).
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#'
#' @export
#'
svine_sim <- function(n, model, qrng = FALSE) {
  assert_that(
    is.number(n),
    inherits(model, "svine_dist"),
    is.flag(qrng)
  )
  
  U <- svinecop_sim_cpp(model$copula, n, qrng, rvinecopulib:::get_seeds())
  to_quantiles(U, model$margins)
}


#' @rdname svine_sim
#' @param data time series of past observations (more recent observations are
#'    at the bottom of the matrix.)
#' @export
svine_sim_conditional <- function(n, data, model, qrng = FALSE, cores = 1) {
  assert_that(
    is.number(n),
    inherits(model, "svine_dist"),
    is.flag(qrng),
    is.number(cores)
  )
  data <- rvinecopulib:::if_vec_to_matrix(data, length(model$margins) == 1)
  data <- to_unif(data, model$margins)
  
  U <- svinecop_sim_conditional_cpp(model$copula, n, data, qrng, cores,
                                    rvinecopulib:::get_seeds())
  to_quantiles(U, model$margins)
}


#' @rdname svine_sim
#' @export
svine_sim_ahead <- function(n, data, model, qrng = FALSE) {
  assert_that(
    is.number(n),
    inherits(model, "svine_dist"),
    is.flag(qrng)
  )
  data <- rvinecopulib:::if_vec_to_matrix(data, length(model$margins) == 1)
  data <- to_unif(data, model$margins)
  
  U <- svinecop_sim_ahead_cpp(model$copula, n, data, qrng, 
                              rvinecopulib:::get_seeds())
  to_quantiles(U, model$margins)
}

#' Log-likelihood for S-vine models
#' 
#' @param x the data.
#' @param model model inherting from class [svine_dist].
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @export
svine_loglik <- function(x, model, cores = 1) {
  assert_that(inherits(model, "svine_dist"))
  x <- rvinecopulib:::if_vec_to_matrix(x, length(model$margins) == 1)
  assert_that(ncol(x) == length(model$margins))
  
  ll_marg <- sapply(
    seq_len(ncol(x)), 
    function(j) sum(univariateML::dml(x[, j], model$margins[[j]], log = TRUE))
  )
  u <- to_unif(x, model$margins)
  ll_cop <- svinecop_loglik_cpp(u, model$copula, cores)
  ll_marg + ll_cop
}

# #' @export
# svine_scores <- function(x, model, cores = 1) {
#   assert_that(inherits(model, "svine_dist"))
#   x <- rvinecopulib:::if_vec_to_matrix(x, length(model$margins) == 1)
#   u <- to_unif(x, model$margins)
#   svinecop_scores_cpp(u, model$copula, cores)
# }

# #' @export
# svine_hessian <- function(x, model, cores = 1) {
#   assert_that(inherits(model, "svine_dist"))
#   x <- rvinecopulib:::if_vec_to_matrix(x, length(model$margins) == 1)
#   u <- to_unif(x, model$margins)
#   svinecop_hessian_cpp(u, model$copula, cores)
# }

# svinecop_cond_cdf <- function(u, conditioned, model, cores = 1) {
#   assert_that(
#     is.count(conditioned),
#     inherits(model, "svinecop_dist"),
#     conditioned <= dim(model$cs_structure)[1]
#   )
#   u <- rvinecopulib:::if_vec_to_matrix(u, dim(model$cs_structure)[1] == 1)
#   svinecop_cond_cdf_cpp(u, conditioned - 1, model, cores)
# }

#' @export
logLik.svine <- function(object, ...) {
  ll <- -2 * object$loglik
  attr(ll, "df") <- object$npars
  ll
}

