find_hamilton_path <- function(u) {
  d <- NCOL(u)
  W <- abs(wdm::wdm(u, method = "kendall"))
  diag(W) <- 0
  order <- which.max(apply(W, 2, max))
  for (k in seq_len(d - 1)) {
    ind <- c(order[1], order[k])
    m <- apply(W[, ind], 2, max)
    v <- apply(W[, ind], 2, which.max)
    imax <- which.max(m)
    W[order, v[imax]] <- W[v[imax], order] <- 0
    order <- if (imax == 1) c(v[imax], order) else c(order, v[imax])
  }

  unname(order)
}

select_mvine <- function(u) {
  d0 <- ncol(u)
  css <- rvinecopulib::dvine_structure(find_hamilton_path(u))
  tau <- wdm::wdm(
    u[-nrow(u), css$order[c(1, d0)]], # t
    u[-1, css$order[c(1, d0)]], # t + 1
    method = "kendall"
  )

  if (abs(tau[1, 1]) < abs(tau[2, 2])) {
    css$order <- rev(css$order)
  }
  list(
    cs_structure = css,
    out_vertices = css$order,
    in_vertices = css$order
  )
}

select_dvine <- function(u) {
  d0 <- ncol(u)
  css <- rvinecopulib::dvine_structure(find_hamilton_path(u))
  tau <- wdm::wdm(
    u[-nrow(u), css$order[c(1, d0)]], # t
    u[-1, css$order[c(1, d0)]], # t + 1
    method = "kendall"
  )

  if (abs(tau[1, 2]) > abs(tau[2, 1])) {
    css$order <- rev(css$order)
  }
  list(
    cs_structure = css,
    out_vertices = rev(css$order),
    in_vertices = css$order
  )
}


#' @importFrom stats sd
select_margin <- function(x, families, criterion) {
  type <- if (all(families == "empirical")) "empirical" else "univariateML"
  out <- if (type == "empirical") {
    F_n <- stats::ecdf(x)
    n <- length(x)
    fit <- list(
      p = function(x) F_n(x) * n / (n + 1),
      q = function(p) stats::quantile(F_n, probs = p)
    )
    attr(fit, "model") <- "empirical"
    attr(fit, "logLik") <- NA
    fit
  } else if (all(families == "std")) {
    par <- c(mean(x), sd(x), 10)
    fn <- function(par) -sum(log(fGarch::dstd(x, par[1], par[2], par[3])))
    opt <- stats::optim(par, fn,
      lower = c(min(x), 0.01 * sd(x), 2.0001),
      upper = c(max(x), 100 * sd(x), 100),
      method = "L-BFGS-B"
    )
    fit <- univariateML::mlstd(1:2)
    attr(fit, "n") <- length(x)
    fit[] <- opt$par
    attr(fit, "logLik") <- -opt$value
    fit
  } else {
    families <- setdiff(families, "empirical")
    fit <- univariateML::model_select(x, families, criterion) |>
      suppressWarnings()
    fit
  }
  out_type <- if (inherits(out, "univariateML") && isFALSE(attr(out, "continuous"))) {
    "discrete"
  } else {
    type
  }
  structure(out, type = out_type, class = c(class(out), "svine_margin"))
}

#' @importFrom stats logLik
#' @exportS3Method
logLik.svine_margin <- function(object, ...) {
  type <- attr(object, "type")
  if (identical(type, "empirical")) {
    structure(NA, df = NA)
  } else if (identical(type, "discrete")) {
    ll <- attr(object, "logLik")
    df <- attr(object, "df")
    structure(if (is.null(ll)) NA else ll,
              df = if (is.null(df)) NA else df)
  } else {
    logLik(object)
  }
}

pmargin <- function(x, model) {
  if (inherits(model, "univariateML")) {
    univariateML::pml(x, model)
  } else {
    model$p(x)
  }
}

# pmargin_sub returns F(x-): the CDF at the largest support point
# strictly less than x — for integer-valued discrete margins this is
# F(x - 1). Companion to pmargin for the doubled-column input layout
# that discrete-margin copula evaluation expects: an n x 2d matrix
# with F(x) in columns 1..d and F(x-) in columns d+1..2d. Three arms:
# univariateML-classed discrete margins (auto-fit path, or user-
# constructed via mlpois/mlnbinom/mlgeom directly) use F(x - 1) via
# pml on shifted inputs; list-form discrete margins (β contract)
# supply $p_sub directly; everything else is continuous and falls
# through to pmargin, emitting redundant shadow columns that the
# copula collapses internally.
#' @noRd
pmargin_sub <- function(x, model) {
  if (inherits(model, "univariateML") && isFALSE(attr(model, "continuous"))) {
    univariateML::pml(x - 1, model)
  } else if (!inherits(model, "univariateML") && is.function(model$p_sub)) {
    model$p_sub(x)
  } else {
    pmargin(x, model)
  }
}

qmargin <- function(p, model) {
  if (inherits(model, "univariateML")) {
    univariateML::qml(p, model)
  } else {
    model$q(p)
  }
}

# check_margin validates a margin object against the β contract.
# univariateML objects pass on class (dispatch goes through
# univariateML::pml / qml without touching $p / $q). All other margins
# must be lists with $p and $q callables, plus $p_sub when
# attr(., "type") is "discrete". j is the column index, included in
# error messages so users can locate the offending margin.
#' @noRd
check_margin <- function(m, j) {
  if (inherits(m, "univariateML")) return(invisible(TRUE))

  if (!is.list(m)) {
    stop(sprintf("margin [%d]: expected a univariateML object or a list with $p and $q callables, got %s.", j, typeof(m)))
  }
  if (!is.function(m$p) || !is.function(m$q)) {
    stop(sprintf("margin [%d]: list-form margins must provide both $p and $q as callables.", j))
  }
  if (identical(attr(m, "type"), "discrete") && !is.function(m$p_sub)) {
    stop(sprintf("margin [%d]: discrete margins must additionally provide $p_sub (a callable returning F(x-)).", j))
  }
  invisible(TRUE)
}
