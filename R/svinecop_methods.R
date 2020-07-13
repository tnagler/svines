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
#'   done in parallel on `cores` batches.
#'   
#' @return An `n`-by-`d`-by`rep` arrray, where `d` is the cross-sectional 
#'   dimension of the model. This reduces to an `n`-by-`d` matrix if `rep == 1`. 
#'
#' @export
#'
svinecop_sim <- function(n, model, past = NULL, rep = 1, qrng = FALSE, cores = 1) {
  assert_that(
    is.number(n),
    inherits(model, "svinecop_dist"),
    is.flag(qrng)
  )
  d0 <- dim(model$cs_structure)[1]
  if (is.null(past)) {
    past <- matrix(NA, 0, 0)
  } else {
    past <- rvinecopulib:::if_vec_to_matrix(past, NCOL(past) == d0)
  }
  
  U <- svinecop_sim_cpp(
    model, n, rep, past, qrng, cores, rvinecopulib:::get_seeds())
  if (rep > 1)
    U <- array(U, dim = c(n, d0, rep))
  if (!is.null(model$names))
    colnames(U) <- simplify_names(model)
  
  U
}

simplify_names <- function(model) {
  nms <- model$names[seq_along(model$in_vertices)]
  nms <- strsplit(nms, "-")
  sapply(nms, function(n) paste(n[-length(n)], collapse = ""))
}

#' Log-likelihood for S-vine copula models
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

#' Log-likelihood scores for S-vine copula models
#' 
#' @param u the data; should have approximately uniform margins..
#' @param model model inherting from class [svinecop_dist].
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @export
svinecop_scores <- function(u, model, cores = 1) {
  assert_that(inherits(model, "svinecop_dist"))
  u <- rvinecopulib:::if_vec_to_matrix(u, dim(model$cs_structure)[1] == 1)
  svinecop_scores_cpp(u, model, cores)
}

#' Expected hessian for S-vine copula models
#' 
#' @param u the data; should have approximately uniform margins..
#' @param model model inherting from class [svinecop_dist].
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @export
svinecop_hessian <- function(u, model, cores = 1) {
  assert_that(inherits(model, "svinecop_dist"))
  u <- rvinecopulib:::if_vec_to_matrix(u, dim(model$cs_structure)[1] == 1)
  svinecop_hessian_cpp(u, model, cores)
}

# svinecop_cond_cdf <- function(u, conditioned, svinecop, cores = 1) {
#   assert_that(
#     is.count(conditioned),
#     inherits(svinecop, "svinecop_dist"),
#     conditioned <= dim(svinecop$cs_structure)[1]
#   )
#   u <- rvinecopulib:::if_vec_to_matrix(u, dim(svinecop$cs_structure)[1] == 1)
#   svinecop_cond_cdf_cpp(u, conditioned - 1, svinecop, cores)
# }