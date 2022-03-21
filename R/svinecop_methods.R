#' Simulate from a S-vine copula model
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
#' @return An `n`-by-`d`-by-`rep` array, where `d` is the cross-sectional 
#'   dimension of the model. This reduces to an `n`-by-`d` matrix if `rep == 1`. 
#'
#' @export
#' 
#' @examples
#' # load data set
#' data(returns)  
#' 
#' # convert to uniform margins
#' u <- pseudo_obs(returns[1:100, 1:3])
#'
#' # fit parametric S-vine copula model with Markov order 1
#' fit <- svinecop(u, p = 1, family_set = "parametric")
#' 
#' pairs(u)   # original data
#' pairs(svinecop_sim(100, rep = 1, model = fit))   # simulated data
#' 
#' # simulate the next day conditionally on the past 500 times
#' pairs(t(svinecop_sim(1, rep = 100, model = fit, past = u)[1, , ]))
svinecop_sim <- function(n, rep, model, past = NULL, qrng = FALSE, cores = 1) {
  assert_that(
    is.number(n),
    inherits(model, "svinecop_dist"),
    is.flag(qrng)
  )

  d0 <- dim(model$cs_structure)[1]
  if (is.null(past)) {
    past <- matrix(NA, 0, 0)
  } else {
    past <- if_vec_to_matrix(past, NCOL(past) == d0)
  }
  
  U <- svinecop_sim_cpp(model, n, rep, past, qrng, cores, get_seeds())
  if (rep > 1)
    U <- array(U, dim = c(n, d0, rep))
  if (!is.null(model$names))
    colnames(U) <- simplify_names(model)
  
  U
}

get_seeds <- function() {
  as.numeric(sprintf("%20.0f", stats::runif(20, 1e+06, 1e+07)))
}

simplify_names <- function(model) {
  nms <- model$names[seq_along(model$in_vertices)]
  nms <- strsplit(nms, "-")
  sapply(nms, function(n) paste(n[-length(n)], collapse = ""))
}

#' Log-likelihood for S-vine copula models
#' 
#' @param u the data; should have approximately uniform margins..
#' @param model model inheriting from class [svinecop_dist].
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @return Returns the log-likelihood of the data for the model.  
#'   
#' @export
#' @examples 
#' # load data set
#' data(returns)  
#' 
#' # convert to uniform margins
#' u <- pseudo_obs(returns[1:100, 1:3])
#'
#' # fit parametric S-vine copula model with Markov order 1
#' fit <- svinecop(u, p = 1, family_set = "parametric")
#' 
#' svinecop_loglik(u, fit)
#' svinecop_scores(u, fit)
#' svinecop_hessian(u, fit)
svinecop_loglik <- function(u, model, cores = 1) {
  assert_that(inherits(model, "svinecop_dist"))
  u <- if_vec_to_matrix(u, dim(model$cs_structure)[1] == 1)
  svinecop_loglik_cpp(u, model, cores)
}

#' Log-likelihood scores for S-vine copula models
#' 
#' @param u the data; should have approximately uniform margins..
#' @param model model inheriting from class [svinecop_dist].
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#'   
#' @returns A matrix containing the score vectors in its rows, where each 
#'   row corresponds to one observation (row in `u`). The columns correspond 
#'   to model parameters in the order: 
#'   copula parameters of first tree, copula parameters of
#'   second tree, etc. Duplicated parameters in the copula model are omitted.
#'   
#' @seealso [svinecop_hessian]
#' 
#' @export
#' @examples 
#' # load data set
#' data(returns)  
#' 
#' # convert to uniform margins
#' u <- pseudo_obs(returns[1:100, 1:3])
#'
#' # fit parametric S-vine copula model with Markov order 1
#' fit <- svinecop(u, p = 1, family_set = "parametric")
#' 
#' svinecop_loglik(u, fit)
#' svinecop_scores(u, fit)
#' svinecop_hessian(u, fit)
svinecop_scores <- function(u, model, cores = 1) {
  assert_that(inherits(model, "svinecop_dist"))
  u <- if_vec_to_matrix(u, dim(model$cs_structure)[1] == 1)
  svinecop_scores_cpp(u, model, cores)
}

#' Expected hessian for S-vine copula models
#'
#' @param u the data; should have approximately uniform margins..
#' @param model model inheriting from class [svinecop_dist].
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @return Returns the observed Hessian matrix. Rows/columns correspond to to
#'   model parameters in the order: copula parameters of first tree, copula
#'   parameters of second tree, etc. Duplicated parameters in the copula model
#'   are omitted.
#'
#' @seealso [svinecop_scores]
#'
#' @export
#' @examples
#' # load data set
#' data(returns)
#'
#' # convert to uniform margins
#' u <- pseudo_obs(returns[1:100, 1:3])
#'
#' # fit parametric S-vine copula model with Markov order 1
#' fit <- svinecop(u, p = 1, family_set = "parametric")
#'
#' svinecop_loglik(u, fit)
#' svinecop_scores(u, fit)
#' svinecop_hessian(u, fit)
svinecop_hessian <- function(u, model, cores = 1) {
  assert_that(inherits(model, "svinecop_dist"))
  u <- if_vec_to_matrix(u, dim(model$cs_structure)[1] == 1)
  svinecop_hessian_cpp(u, model, cores)
}

# svinecop_cond_cdf <- function(u, conditioned, svinecop, cores = 1) {
#   assert_that(
#     is.count(conditioned),
#     inherits(svinecop, "svinecop_dist"),
#     conditioned <= dim(svinecop$cs_structure)[1]
#   )
#   u <- if_vec_to_matrix(u, dim(svinecop$cs_structure)[1] == 1)
#   svinecop_cond_cdf_cpp(u, conditioned - 1, svinecop, cores)
# }
