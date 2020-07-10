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
    u[-1, css$order[c(1, d0)]],       # t + 1
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
    u[-1, css$order[c(1, d0)]],       # t + 1
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

get_allowed_margin_families <- function() {
  uml_funs <- univariateML:::densities
  uml_funs <- uml_funs[grepl("^ml", uml_funs)]
  substring(uml_funs, first = 3)
}

allowed_margin_families <- get_allowed_margin_families()

select_margin <- function(x, families, criterion, na.rm = FALSE) {
  mlf <- sapply(
    paste0("univariateML::ml", families),
    function(x) eval(parse(text = x))
  )
  fits <- lapply(mlf, function(f) try(f(x, na.rm = na.rm), silent = TRUE))

  ## catch out-of-bounds errors (and similar)
  error_inds <- sapply(fits, function(fit) inherits(fit, "try-error"))
  error_msgs <- sapply(fits[error_inds], as.character)
  if (all(error_inds)) {
    details <- paste0("(", names(error_msgs), ") ", error_msgs)
    stop("couldn't fit any model.\n", details)
  }

  ## select best model
  fits <- fits[!error_inds]
  crits <- switch(
    criterion,
    "loglik" = -sapply(fits, stats::logLik),
    "aic"    = sapply(fits, stats::AIC),
    "bic"    = sapply(fits, stats::BIC),
    "mbicv"  = sapply(fits, stats::BIC)
  )
  fits[[which.min(crits)]]
}
