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
#' @param n how many steps of the time series to simulate.
#' @param model a S-vine copula model object (inheriting from [svinecop_dist]).
#' @param past (optional) matrix of past observations. If provided, time series 
#'   are simulated conditional on the past. 
#' @param rep number of replications; `rep` time series of length `n` are 
#'   generated.
#' @param qrng if `TRUE`, generates quasi-random numbers using the multivariate
#' Generalized Halton sequence up to dimension 300 and the Generalized Sobol
#' sequence in higher dimensions (default `qrng = FALSE`).
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel.
#'   
#' @return An `n`-by-`d`-by`rep` arrray, where `d` is the cross-sectional 
#'   dimension of the model. This reduces to an `n`-by-`d` matrix if `rep == 1`. 
#'
#' @export
#'
svine_sim <- function(n, model, past = NULL, rep = 1, qrng = FALSE, cores = 1) {
  assert_that(
    is.number(n),
    inherits(model, "svine_dist"),
    is.flag(qrng)
  )
  
  if (!is.null(past)) {
    d0 <- length(model$margins)
    past <- rvinecopulib:::if_vec_to_matrix(past, length(past) != d0)
    past <- to_unif(past, model$margins)
  }
  U <- svinecop_sim(n, model$copula, past, rep, qrng, cores)
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
  sum(ll_marg) + ll_cop
}


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

