#' Stationary vine distribution models
#'
#' Automated fitting or creation of custom S-vine distribution models
#'
#' @param data a matrix or data.frame of data.
#' @param p the Markov order.
#' @param margin_families either a vector of [univariateML] families to select 
#'   from (used for every margin) or a list with one entry for every variable.
#' @param ... arguments passed to `svinecop()`.
#'
#' @importFrom assertthat assert_that is.scalar is.string is.number is.flag
#' @importFrom rvinecopulib as_rvine_structure
#' @export
#' 
#' @examples
#' # load data set
#' data(returns)  
#'
#' # convert to pseudo observations with empirical cdf for marginal distributions
#' u <- pseudo_obs(returns[, 1:3]) 
#'
#' # fit parametric S-vine copula model with Markov order 1
#' fit <- svinecop(u, p = 1, family_set = "parametric")
#' fit 
#' summary(fit)
#' plot(fit)
#' contour(fit)
#' logLik(fit)

svine <- function(data, p, margin_families = "norm", selcrit = "aic", ...) {
  if (is.list(data)) {
    if (any(sapply(data, is.factor)))
      stop("discrete data not yet yupported")
  }
  data <- as.matrix(data)
  d <- ncol(data)
  if (!is.list(margin_families))
    margin_families <- lapply(seq_len(d), function(j) margin_families)
  assert_that(length(margin_families) == ncol(data))
  
  margin_crit <- ifelse(selcrit == "mbicv", "bic", selcrit)
  margins <- lapply(
    seq_len(d),
    function(j) select_margin(data[, j], margin_families[[j]], margin_crit)
  )
  names(margins) <- colnames(data)
  
  u <- sapply(
    seq_len(d),
    function(j) pml(data[, j], margins[[j]])
  )
  
  copula <- svinecop(u, p = p, selcrit = selcrit, ...)
  loglik <- sum(sapply(margins, logLik)) + logLik(copula)
  npars  <- sum(sapply(margins, length)) + copula$npars
  
  structure(
    list(
      margins = margins, 
      copula = copula,
      d = d,
      loglik = loglik,
      npars = npars
    ),
    class = c("svine", "svine_dist")
  )
}


#' Custom S-vine distribution models
#'
#' @param margins A list of length `d` containing `univariateML` objects.
#' @param copula the copula model; an object of class `svinecop_dist` with 
#'   cross-sectional dimension `d`.
#'
#' @export
svine_dist <- function(margins, copula) {
  assert_that(
    is.list(margins),
    inherits(copula, "vinecop_dist"),
    length(margins) == dim(copula$cs_structure)["dim"]
  )
  structure(
    list(
      margins = margins, 
      copula = copula,
      loglik = NA,
      npars = NA
    ),
    class = c("svine_dist")
  )
}

#' @export
print.svine_dist <- function(x, ...) {
  cat(
    length(x$margins),
    "-dimensional S-vine distribution model of order p = ", x$copula$p,
    " ('svine_dist')\n",
    sep = ""
  )
  invisible(x)
}

#' @export
summary.svine_dist <- function(object, ...) {
  list(
    margins = get_svine_dist_margin_summary(object$margins),
    copula = summary(object$copula)
  )
}

get_svine_dist_margin_summary <- function(margins) {
  nms <- names(margins)
  if (is.null(nms))
    nms <- paste0("V", seq_along(margins))
  df <- data.frame(
    margin = seq_along(nms),
    name = nms,
    model = sapply(margins, function(x) attr(x, "model")),
    parameters = NA,
    loglik = sapply(margins, function(x) attr(x, "logLik"))
  )
  for (m in seq_along(margins)) 
    df$parameters[m] <- list(coef(margins[[m]]))
  class(df) <- c("summary_df", class(df))
  df
}


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


to_quantiles <- function(u, margins) {
  x <- sapply(seq_len(ncol(u)), function(j) qml(u[, j], margins[[j]]))
  var_names <- names(margins)
  if (!is.null(var_names))
    colnames(x) <- var_names
  x
}

to_unif <- function(x, margins) {
  u <- sapply(seq_len(ncol(x)), function(j) pml(x[, j], margins[[j]]))
  var_names <- names(margins)
  if (!is.null(var_names))
    colnames(u) <- var_names
  u
}