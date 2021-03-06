#' Simulate from a S-vine model
#'
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
#'   done parallel over replications.
#'   
#' @return An `n`-by-`d`-by`rep` array, where `d` is the cross-sectional 
#'   dimension of the model. This reduces to an `n`-by-`d` matrix if `rep == 1`. 
#'
#' @export
#' @examples
#' # load data set
#' data(returns)  
#' returns <- returns[1:100, 1:3]
#'
#' # fit parametric S-vine model with Markov order 1
#' fit <- svine(returns, p = 1, family_set = "parametric")
#' 
#' pairs(returns)  # original data
#' pairs(svine_sim(100, rep = 1, model = fit))   # simulated data
#' 
#' # simulate the next day conditionally on the past 500 times
#' pairs(t(svine_sim(1, rep = 100, model = fit, past = returns)[1, , ]))
svine_sim <- function(n, rep, model, past = NULL, qrng = FALSE, cores = 1) {
  assert_that(
    is.number(n),
    inherits(model, "svine_dist"),
    is.flag(qrng)
  )
  if (!is.null(past)) {
    assert_that(NCOL(past) == length(model$margins))
  }

  if (!is.null(past)) {
    d0 <- length(model$margins)
    past <- rvinecopulib:::if_vec_to_matrix(past, length(past) != d0)
    past <- to_unif(past, model$margins)
  }
  U <- svinecop_sim(n, rep, model$copula, past, qrng, cores)
  to_quantiles(U, model$margins)
}

#' Log-likelihood for S-vine models
#' 
#' @param x the data.
#' @param model model inheriting from class [svine_dist].
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @export
#' @examples 
#' # load data set
#' data(returns)  
#'
#' # fit parametric S-vine model with Markov order 1
#' fit <- svine(returns[1:100, 1:3], p = 1, family_set = "parametric")
#' 
#' svine_loglik(returns[1:100, 1:3], fit)
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
  ll <- object$loglik
  attr(ll, "df") <- object$npars
  ll
}

