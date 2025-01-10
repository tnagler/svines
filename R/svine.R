#' Stationary vine distribution models
#'
#' Automated fitting or creation of custom S-vine distribution models
#'
#' @param data a matrix or data.frame of data.
#' @param p the Markov order.
#' @param margin_families either a vector of [univariateML::univariateML_models] to select 
#'   from (used for every margin) or a list with one entry for every variable. 
#'   Can also be `"empirical"` for empirical cdfs.
#' @param selcrit criterion for family selection, either `"loglik"`, `"aic"`,
#'   `"bic"`, `"mbicv"`.
#' @param ... arguments passed to `svinecop()`.
#' 
#' @return Returns the fitted model as an object with classes 
#'   `svine` and [svine_dist]. A list with entries 
#'   - `$margins`: list of marginal models from [univariateML::univariateML_models],
#'   - `$copula`: an object of `svinecop_dist`.
#'   
#' @seealso [svine_dist], [svine_loglik], [svine_sim], [svine_bootstrap_models]
#'
#' @importFrom assertthat assert_that is.scalar is.string is.number is.flag
#' @importFrom rvinecopulib as_rvine_structure
#' @export
#' 
#' @examples
#' # load data set
#' data(returns)  
#'
#' # fit parametric S-vine model with Markov order 1
#' fit <- svine(returns[1:100, 1:3], p = 1, family_set = "parametric")
#' fit 
#' summary(fit)
#' plot(fit$copula)
#' contour(fit$copula)
#' logLik(fit)
#' 
#' pairs(svine_sim(500, rep = 1, fit))
svine <- function(data, p, margin_families = univariateML::univariateML_models, 
                  selcrit = "aic", ...) {
  if (is.list(data)) {
    if (any(sapply(data, is.factor)))
      stop("discrete data not yet yupported")
  }
  data <- as.matrix(data)
  d <- ncol(data)
  if (!is.list(margin_families))
    margin_families <- lapply(seq_len(d), function(j) margin_families)
  assert_that(length(margin_families) == d)
  
  margin_crit <- ifelse(selcrit == "mbicv", "bic", selcrit)
  margins <- lapply(
    seq_len(d),
    function(j) select_margin(data[, j], margin_families[[j]], margin_crit)
  )
  names(margins) <- colnames(data)
  
  u <- to_unif(data, margins)
  copula <- svinecop(u, p = p, selcrit = selcrit, ...)
  if (!all(unlist(margin_families) == "empirical")) {
    loglik <- sum(sapply(margins, logLik.svine_margin)) + stats::logLik(copula)
    npars <- sum(sapply(margins, length)) + copula$npars
  } else {
    loglik <- npars <- NA
  }
  
  structure(
    list(
      margins = margins, 
      copula = copula,
      data = data,
      d = d,
      loglik = loglik,
      npars = npars
    ),
    class = c("svine", "svine_dist")
  )
}


#' Custom S-vine distribution models
#' 
#'
#' @param margins A list of length `d` containing `univariateML` objects.
#' @param copula the copula model; an object of class `svinecop_dist` with 
#'   cross-sectional dimension `d`.
#'
#' @return Returns the model as an object with class  `svine_dist`. 
#'   A list with entries 
#'   - `$margins`: list of marginal models from [univariateML::univariateML_models],
#'   - `$copula`: an object of `svinecop_dist`.
#'   
#' @seealso [svine_dist], [svine_loglik], [svine_sim], [svine_bootstrap_models]
#' 
#' @export
#' @examples 
#' ## marginal objects
#' # create dummy univariateML models
#' univ1 <- univ2 <- univariateML::mlnorm(rnorm(10))
#' 
#' # modify the parameters to N(5, 10) and N(0, 2) distributions
#' univ1[] <- c(5, 10)
#' univ2[] <- c(0, 2)
#' 
#' ## copula óbject
#' cs_struct <- cvine_structure(1:2)
#' pcs <- list(
#'   list(  # first tree
#'     bicop_dist("clayton", 0, 3), # cross sectional copula
#'     bicop_dist("gaussian", 0, -0.1)  # serial copula
#'   ),
#'   list(  # second tree
#'     bicop_dist("gaussian", 0, 0.2), bicop_dist("indep")  
#'   ),
#'   list( # third tree
#'     bicop_dist("indep")
#'   )
#' )
#' 
#' cop <- svinecop_dist(
#'   pcs, cs_struct, p = 1, out_vertices = 1:2, in_vertices = 1:2)
#'     
#' model <- svine_dist(margins = list(univ1, univ2), copula = cop)
#' summary(model)
#' 
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
    if (attr(margins[[m]], "model") != "empirical") {
      df$parameters[m] <- list(
        stats::setNames(as.numeric(margins[[m]]), names(margins[[m]]))
      )
    }
  }
  class(df) <- c("summary_df", class(df))
  df
}


to_quantiles <- function(u, margins) {
  x <- u
  for (j in seq_len(ncol(u))) {
    if (is.matrix(u)) {
      x[, j] <- qmargin(u[, j], margins[[j]])
    } else {
      x[, j, ] <- qmargin(u[, j, ], margins[[j]])
    }
  }
  
  var_names <- names(margins)
  if (!is.null(var_names))
    colnames(x) <- var_names
  
  x
}

to_unif <- function(x, margins) {
  u <- sapply(
    seq_len(ncol(x)), 
    function(j) pmargin(x[, j], margins[[j]])
  )
  var_names <- names(margins)
  if (!is.null(var_names))
    colnames(u) <- var_names
  u
}