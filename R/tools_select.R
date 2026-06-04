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


# univariateML_panel returns the family names whose ml<name> entries
# in univariateML::univariateML_metadata carry a support interval over
# Z (discrete) or R (continuous), per the var_type argument. The names
# are returned without the "ml" prefix, matching univariateML_models'
# convention and the family-set inputs that model_select expects.
#' @noRd
univariateML_panel <- function(var_type) {
  meta <- univariateML::univariateML_metadata
  marker <- if (var_type == "d") "Z" else "R"
  fams <- names(meta)[vapply(
    meta,
    function(m) identical(attr(m[["support"]], "type"), marker),
    logical(1)
  )]
  sub("^ml", "", fams)
}

#' @importFrom stats sd
select_margin <- function(x, families, criterion, var_type = "c", j = NA_integer_) {
  var_type <- match.arg(var_type, c("c", "d"))
  type <- if (all(families == "empirical")) "empirical" else "univariateML"

  prefix <- if (is.na(j)) "" else sprintf("column %d: ", j)

  if (var_type == "d") {
    if (all(families == "std")) {
      stop(sprintf("%sfamilies = \"std\" is not supported for var_type = \"d\".", prefix))
    }
  }

  out <- if (type == "empirical") {
    F_n <- stats::ecdf(x)
    n <- length(x)
    fit <- list(
      p = function(x) F_n(x) * n / (n + 1),
      q = function(p) stats::quantile(F_n, probs = p),
      p_sub = if (var_type == "d") function(x) F_n(x - 1L) * n / (n + 1) else NULL
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
    families <- intersect(families, univariateML_panel(var_type))
    if (length(families) == 0L) {
      stop(sprintf("%sno univariateML families remain after filtering to var_type = \"%s\".",
                   prefix, var_type))
    }
    fit <- univariateML::model_select(x, families, criterion) |>
      suppressWarnings()
    fit
  }
  out_type <- if (inherits(out, "univariateML") && isFALSE(attr(out, "continuous"))) {
    "discrete"
  } else if (type == "empirical" && var_type == "d") {
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

# pmargin_sub returns F(x-): the left limit of the CDF at x.
# For integer-valued discrete margins, F(x-) = F(x - 1).
# Used alongside pmargin to build the doubled input matrix that
# discrete copula evaluation requires: F(x) and F(x-) for each
# variable, side by side. Three dispatch arms:
# (1) univariateML discrete margins: F(x - 1) via pml on shifted inputs.
# (2) List-form discrete margins (e.g. hand-built or empirical):
#     delegates to $p_sub, which the margin object supplies directly.
# (3) Continuous margins: returns F(x) unchanged — the shadow column
#     is redundant but the copula collapses it correctly.
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

# Validates a margin object before use in svine_dist().
# univariateML objects pass automatically. All other margins must be
# lists with callable $p and $q; discrete margins additionally need $p_sub.
# j is the column index, used in error messages.
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
