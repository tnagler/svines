#' Stationary vine distribution models
#'
#' Automated fitting or creation of custom S-vine distribution models
#'
#' @param data a matrix or data.frame of data.
#' @param p the Markov order.
#' @param margin_families either a vector of [univariateML] families to select 
#'   from (used for every margin) or a list with one entry for every variable.
#' @param selcrit criterion for family selection, either `"loglik"`, `"aic"`,
#'   `"bic"`, `"mbicv"`.
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

svine <- function(data, p, margin_families = NA, selcrit = "aic", ...) {
  if (is.list(data)) {
    if (any(sapply(data, is.factor)))
      stop("discrete data not yet yupported")
  }
  if (any(is.na(margin_families))) {
    margin_families <- allowed_margin_families
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
  
  u <- to_unif(data, margins)
  copula <- svinecop(u, p = p, selcrit = selcrit, ...)
  loglik <- sum(sapply(margins, stats::logLik)) + stats::logLik(copula)
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
  for (m in seq_along(margins)) {
    df$parameters[m] <- list(
      stats::setNames(as.numeric(margins[[m]]), names(margins[[m]]))
    )
  }
  class(df) <- c("summary_df", class(df))
  df
}


to_quantiles <- function(u, margins) {
  x <- sapply(
    seq_len(ncol(u)),
    function(j) univariateML::qml(u[, j], margins[[j]])
  )
  var_names <- names(margins)
  if (!is.null(var_names))
    colnames(x) <- var_names
  x
}

to_unif <- function(x, margins) {
  u <- sapply(
    seq_len(ncol(x)), 
    function(j) univariateML::pml(x[, j], margins[[j]])
  )
  var_names <- names(margins)
  if (!is.null(var_names))
    colnames(u) <- var_names
  u
}