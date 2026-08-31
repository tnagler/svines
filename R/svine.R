#' Stationary vine distribution models
#'
#' Automated fitting or creation of custom S-vine distribution models
#'
#' @param data a matrix or data.frame of data.
#' @param p the Markov order.
#' @param var_types variable types; a character vector with one entry per
#'   variable: `"c"` for continuous, `"d"` for discrete. Defaults to
#'   all-continuous when missing.
#' @param margin_families either a character vector of
#'   [univariateML::univariateML_models] to select from for every margin, or a
#'   list with one entry for every variable. For a discrete variable, the
#'   corresponding entry must contain only suitable discrete families. Use
#'   `"empirical"` for empirical CDFs.
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
#' @references Nagler, T., Krüger, D., and Min, A. (2022). Stationary vine
#'   copula models for multivariate time series. *Journal of Econometrics*,
#'   227(2), 305--324. \doi{10.1016/j.jeconom.2021.11.015}.
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
svine <- function(data, p, var_types,
                  margin_families = univariateML::univariateML_models,
                  selcrit = "aic", ...) {
  if (is.list(data)) {
    if (any(sapply(data, is.factor)))
      stop("factor columns are not supported; convert to integer counts before fitting.")
  }
  data <- as.matrix(data)
  d <- ncol(data)
  if (missing(var_types))
    var_types <- rep("c", d)
  assert_that(
    is.character(var_types),
    length(var_types) == d,
    all(var_types %in% c("c", "d"))
  )
  if (!is.list(margin_families))
    margin_families <- lapply(seq_len(d), function(j) margin_families)
  assert_that(length(margin_families) == d)

  margin_crit <- ifelse(selcrit == "mbicv", "bic", selcrit)
  margins <- lapply(
    seq_len(d),
    function(j) select_margin(data[, j], margin_families[[j]], margin_crit,
                              var_type = var_types[j], j = j)
  )
  names(margins) <- colnames(data)

  u <- to_unif(data, margins)
  copula <- svinecop(u, p = p, var_types = var_types, selcrit = selcrit, ...)
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
#' @param margins A list with one marginal distribution per variable.
#'   Each element is either a `univariateML` object or a list with
#'   callables `$p` (CDF), `$q` (quantile function), and for discrete
#'   margins `$p_sub` (left-limit CDF, i.e. F(x-)). 
#'   Discrete list-form margins must also carry
#'   `attr(., "type") = "discrete"` and `class` including
#'   `"svine_margin"`.
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
#' ## copula object
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
  for (j in seq_along(margins)) check_margin(margins[[j]], j)
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
  u <- vapply(
    seq_len(ncol(x)),
    function(j) pmargin(x[, j], margins[[j]]),
    FUN.VALUE = numeric(nrow(x))
  )
  if (!is.matrix(u))
    u <- matrix(u, nrow = nrow(x))
  var_names <- names(margins)
  if (!is.null(var_names))
    colnames(u) <- var_names

  # When any margin is discrete, the copula evaluation requires the
  # compact n x (d + n_discrete) input layout: F(x) in columns 1..d,
  # followed by one F(x-) left-limit column per discrete margin, in
  # var_types order. Continuous-only margin sets fall through with the
  # original n x d output unchanged.
  is_disc <- vapply(
    margins,
    function(m) identical(attr(m, "type"), "discrete"),
    logical(1)
  )
  if (any(is_disc)) {
    disc_idx <- which(is_disc)
    u_sub <- vapply(
      disc_idx,
      function(j) pmargin_sub(x[, j], margins[[j]]),
      FUN.VALUE = numeric(nrow(x))
    )
    if (!is.matrix(u_sub))
      u_sub <- matrix(u_sub, nrow = nrow(x))
    u <- cbind(u, u_sub)
    if (!is.null(var_names))
      colnames(u) <- c(var_names, paste0(var_names[disc_idx], "_sub"))
  }
  u
}
