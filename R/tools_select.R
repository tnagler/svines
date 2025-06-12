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
  structure(out, type = type, class = c(class(out), "svine_margin"))
}

#' @importFrom stats logLik
#' @exportS3Method 
logLik.svine_margin <- function(object, ...) {
  if (attr(object, "type") == "empirical") {
    structure(NA, df = NA)
  } else {
    logLik(object)
  }
}

pmargin <- function(x, model) {
  if (identical(attr(model, "type"), "empirical")) {
    model$p(x)
  } else {
    univariateML::pml(x, model)
  }
}

qmargin <- function(p, model) {
  if (identical(attr(model, "type"), "empirical")) {
    model$q(p)
  } else {
    univariateML::qml(p, model)
  }
}
