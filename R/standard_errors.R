#' Asymptotic variance of parameter estimates of an S-vine model
#'
#' @param x the data.
#' @param model S-vine model (inheriting from [svine_dist]).
#' @param n_lags the number of lags to use when computing the Fisher information
#'   matrix.
#' @param cores the number of cores to use.
#'
#' @return A k-by-k matrix, where k is the total number of parameters in the
#'   model. Parameters are ordered as follows:
#'   marginal parameters, copula parameters of first tree, copula parameters of
#'   second tree, etc. Duplicated parameters in the copula model are omitted.
#'
#'
#' @export
#' @examples
#' # load data set
#' data(returns)
#' dat <- returns[1:100, 1:2]
#'
#' # fit parametric S-vine model with Markov order 1
#' fit <- svine(dat, p = 1, family_set = "parametric")
#'
#' # asymptotic variance matrix
#' avar <- svine_avar(dat, fit)
#'
#' # standard errors
#' sqrt(diag(avar))
svine_avar <- function(x, model, n_lags = floor(sqrt(NROW(x))), cores = 1) {
  assert_that(inherits(model, "svine_dist"))
  scores <- svine_scores(x, model, cores)
  I <- stats::cov(scores)
  if (n_lags > 0) {
    lagged_covs <- lapply(
      seq_len(n_lags) - 1,
      function(l) {
        sig <- stats::cov(
          scores[-seq_len(l + 1), ],
          scores[-((NROW(scores) - l):NROW(scores)), ]
        )
        sig + t(sig)
      }
    )
    I <- I + Reduce("+", lagged_covs)
  }
  H <- svine_hessian(x, model, cores)
  Hi <- solve(H)
  avar <- Hi %*% I %*% t(Hi) / NROW(x)
  
  # make sure it's positive definite (rounding errors!)
  eig <- eigen(avar)
  if (any(eig$values <= 0)) {
    eig$values <- abs(eig$values)
    avar <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  }
  
  avar
}

#' Score function of S-vine distribution models
#'
#' @param x the data.
#' @param model S-vine model (inheriting from [svine_dist]).
#' @param cores number of cores to use.
#'
#' @return A returns a `n`-by-`k` matrix, where `n = NROW(x)` and `k` is the
#'   total number of parameters in the
#'   model. Parameters are ordered as follows:
#'   marginal parameters, copula parameters of first tree, copula parameters of
#'   second tree, etc. Duplicated parameters in the copula model are omitted.
#'
#' @examples
#' data(returns)
#' dat <- returns[1:100, 1:2]
#'
#' # fit parametric S-vine model with Markov order 1
#' model <- svine(dat, p = 1, family_set = "parametric")
#'
#' # Implementation of asymptotic variances
#' I <- cov(svine_scores(dat, model))
#' H <- svine_hessian(dat, model)
#' Hi <- solve(H)
#' Hi %*% I %*% t(Hi) / nrow(dat)
#' @export
svine_scores <- function(x, model, cores = 1) {
  assert_that(inherits(model, "svine_dist"))
  x <- rvinecopulib:::if_vec_to_matrix(x, length(model$margins) == 1)
  u <- as.matrix(to_unif(x, model$margins))
  
  S_mrg <- scores_mrg(x, model)
  S_cop <- svinecop_scores(u, model$copula, cores = cores)
  if (model$copula$p > 0) {
    S_mrg <- S_mrg[-seq_len(model$copula$p), ]
  }
  cbind(S_mrg, S_cop)
}

#' Expected hessian of an S-vine distribution models
#'
#' @param x the data.
#' @param model S-vine model (inheriting from [svine_dist]).
#' @param cores number of cores to use.
#'
#' @return A returns a `k`-by-`k` matrix, where `k` is the
#'   total number of parameters in the
#'   model. Parameters are ordered as follows:
#'   marginal parameters, copula parameters of first tree, copula parameters of
#'   second tree, etc. Duplicated parameters in the copula model are omitted.
#'
#' @examples
#' data(returns)
#' dat <- returns[1:100, 1:2]
#'
#' # fit parametric S-vine model with Markov order 1
#' model <- svine(dat, p = 1, family_set = "parametric")
#'
#' # Implementation of asymptotic variances
#' I <- cov(svine_scores(dat, model))
#' H <- svine_hessian(dat, model)
#' Hi <- solve(H)
#' Hi %*% I %*% t(Hi) / nrow(dat)
#' @export
svine_hessian <- function(x, model, cores = 1) {
  assert_that(inherits(model, "svine_dist"))
  x <- rvinecopulib:::if_vec_to_matrix(x, length(model$margins) == 1)
  u <- as.matrix(to_unif(x, model$margins))
  
  H_mrg <- hessian_mrg(x, model)
  H_mxd <- hessian_mxd(x, model, cores = cores)
  H_cop <- svinecop_hessian(u, model$copula, cores = cores)
  rbind(
    cbind(H_mrg, H_mxd),
    cbind(matrix(0, ncol(H_mxd), nrow(H_mxd)), H_cop)
  )
}

#' Simulate models from the asymptotic distribution
#'
#' Simulates `n` \emph{iid} models, where each model
#' has parameters drawn from the asymptotic distribution.
#' @param n number of models to simulate.
#' @param model S-vine model (inheriting from [svine_dist]).
#' @param cores the number of cores used for computation of the asymptotic
#'   variance.
#' @param ... passed to [`svine_avar()`].

#' @examples
#' data(returns)
#' dat <- returns[1:100, 1:2]
#'
#' # fit parametric S-vine model with Markov order 0
#' fit <- svine(dat, p = 0, family_set = "parametric")
#'
#' new <- svine_sim_se_models(2, fit)
#' summary(new[[1]])
#' summary(new[[2]])
#' @export
svine_sim_se_models <- function(n, model, cores = 1, ...) {
  assert_that(inherits(model, "svine_dist"))
  V <- svine_avar(model$data, model = model, cores = cores, ...)
  
  # make sure it's positive definite before simualting
  eig <- eigen(V)
  V <- V + diag(diag(V)) * 0.1
  
  par <- svine_get_pars(model)
  replicate(
    n,
    svine_set_pars(model, MASS::mvrnorm(1, par, V)),
    simplify = FALSE
  )
}


sim_multipliers <- function(n, ell) {
  b <- ceiling((ell + 1) / 2)
  Z <- rnorm(n + 2 * b - 2)
  w <- kern((1:ell - b) / b)
  w <- w / sqrt(sum(w^2))
  xi <- sapply(1:n, function(i) sum(w * Z[i - 1 + 1:ell]))
  xi + 1 - mean(xi)
}

kern <- function(x) {
  K <- (1 - 6 * x^2 + 6 * abs(x)^3) * (abs(x) <= 1/2)
  K + 2 * (1 - abs(x))^3 * (abs(x) > 1 / 2) * (abs(x) <= 1)
}

bootstrap_pobs <- function(u, xi) {
  u <- as.matrix(u)
  for (j in seq_len(ncol(u))) {
    s <- sort(u[, j], index.return = TRUE)
    u[, j] <- cumsum(xi[s$ix])[order(s$ix)]
    u[, j] <- pmax(u[, j], 1e-10) / (max(u[, j]) + 1e-10)
  }
  u
}


svine_bootstrap_semipar <- function(n_models, model, ...) {
  n <- nrow(model$data)
  np <- n - model$copula$p
  ell <- ceiling(6 * n ^ (1 / 5))
  
  u <- to_unif(model$data, model$margins)
  par <- svinecop_get_pars(model$copula)
  phi <- svinecop_scores(u, model$copula)
  Hs <- svinecop_hessian_sep_cpp(u, model$copula, 1)
  Hs <- array(unlist(Hs), dim = c(ncol(phi), ncol(phi), length(Hs)))
  Hi <- solve(apply(Hs, 1:2, mean))

  model$data <- NULL
  models <- replicate(n_models, model, simplify = FALSE)
  for (b in seq_along(models)) {
    xi <- sim_multipliers(n, ell)
    u_tilde <- bootstrap_pobs(u, xi)
    
    xi <- c(xi[seq_len(np)])
    # Hi <- apply(Hs, 1:2, function(x) mean(x * xi))
    # Hi <- solve(Hi + diag(1e-10, nrow(Hi)))
    phi_tilde <- colMeans(xi * svinecop_scores(u_tilde, model$copula))
    
    par_b <- par - Hi %*% phi_tilde
    models[[b]]$copula <- svinecop_set_pars(models[[b]]$copula, par_b)
    models[[b]]$margins <- lapply(
      seq_along(model$margins),
      function(j) {
        x_tilde <- model$margins[[j]]$q(u_tilde[, j])
        select_margin(x_tilde, "empirical", "")
      }
    )
  }
  models
}

svine_bootstrap_par <- function(n_models, model, ...) {
  n <- nrow(model$data)
  np <- n - model$copula$p
  ell <- ceiling(5 * n ^ (1 / 5))
  
  par <- svine_get_pars(model)
  phi <- svine_scores(model$data, model)
  Hi <- solve(svine_hessian(model$data, model))
  
  model$data <- NULL
  models <- replicate(n_models, model, simplify = FALSE)
  for (b in seq_along(models)) {
    xi <- sim_multipliers(n, ell)
    phi_tilde <- c(xi[seq_len(np)]) * phi
    par_b <- par - Hi %*% colMeans(phi_tilde)
    models[[b]] <- svine_set_pars(models[[b]], par_b)
  }
  models
}

#' @export
svine_bootstrap_models <- function(n_models, model, ...) {
  assert_that(inherits(model, "svine_dist"))
  if (attr(model$margins[[1]], "type") == "empirical") {
    svine_bootstrap_semipar(n_models, model, ...)
  } else {
    svine_bootstrap_par(n_models, model, ...)
  }
}

svine_get_pars <- function(model) {
  assert_that(inherits(model, "svine_dist"))
  
  pars_mrg <- tryCatch(
    lapply(model$margins, as.numeric),
    error = function(e) NULL,
    silent = TRUE
  )
  pars_cop <- lapply(
    model$copula$pair_copulas,
    function(tree) {
      lapply(
        seq_len(min(length(model$margins), length(tree))),
        function(j) tree[[j]]$parameters
      )
    }
  )
  
  unname(c(unlist(pars_mrg), unlist(pars_cop)))
}

svinecop_get_pars <- function(model) {
  assert_that(inherits(model, "svinecop_dist"))
  pars_cop <- lapply(
    model$pair_copulas,
    function(tree) {
      lapply(
        seq_len(min(model$cs_structure$d, length(tree))),
        function(j) tree[[j]]$parameters
      )
    }
  )
  unname(unlist(pars_cop))
}

svinecop_set_pars <- function(model, parameters) {
  assert_that(
    inherits(model, "svinecop_dist"),
    length(parameters) == model$npars
  )
  copula <- with_parameters_cop_cpp(model, parameters)
  svinecop_dist(
    pair_copulas = copula$pair_copulas,
    p = copula$p,
    cs_structure = copula$cs_structure,
    out_vertices = copula$out_vertices,
    in_vertices = copula$in_vertices
  )
}

svine_set_pars <- function(model, parameters) {
  assert_that(
    inherits(model, "svine_dist"),
    length(parameters) == model$npars
  )
  
  margins <- with_parameters_mrg(model$margins, parameters)
  parameters <- parameters[-seq_len(sum(sapply(model$margins, length)))]
  copula <- with_parameters_cop_cpp(model$copula, parameters)
  copula <- svinecop_dist(
    pair_copulas = copula$pair_copulas,
    p = copula$p,
    cs_structure = copula$cs_structure,
    out_vertices = copula$out_vertices,
    in_vertices = copula$in_vertices
  )
  
  svine_dist(margins, copula)
}


with_parameters_mrg <- function(margins, parameters) {
  npars_mrg <- sapply(margins, length)
  ixs <- c(0, cumsum(npars_mrg))
  for (m in seq_along(margins)) {
    a <- attributes(margins[[m]])
    margins[[m]] <- parameters[(ixs[m] + 1):ixs[m + 1]]
    if (a$density %in% c("fGarch::dstd", "fGarch::dsstd")) {
      margins[[m]][3] <- max(margins[[m]][3], 2.001)
    } 
    attributes(margins[[m]]) <- a
    attr(margins[[m]], "logLik") <- NA
  }
  margins
}

scores_mrg_1 <- function(x, model) {
  npars <- length(model)
  scores <- matrix(NA, length(x), npars)
  for (p in seq_len(npars)) {
    tmp_model <- model
    
    tmp_model[p] <- model[p] - 1e-3
    f_lwr <- univariateML::dml(x, tmp_model, log = TRUE)
    
    tmp_model[p] <- model[p] + 1e-3
    f_upr <- univariateML::dml(x, tmp_model, log = TRUE)
    
    scores[, p] <- (f_upr - f_lwr) / 2e-3
  }
  scores
}

scores_mrg <- function(x, model) {
  d <- NCOL(x)
  s <- lapply(seq_len(d), function(m) scores_mrg_1(x[, m], model$margins[[m]]))
  do.call(cbind, s)
}

hessian_mrg_1 <- function(x, model, weights) {
  if (!length(weights)) weights <- rep(1, length(x))
  npars <- length(model)
  hessian <- matrix(NA, npars, npars)
  for (p in seq_len(npars)) {
    tmp_model <- model
    
    tmp_model[p] <- model[p] - 1e-3
    s_lwr <- scores_mrg_1(x, tmp_model) 
    
    tmp_model[p] <- model[p] + 1e-3
    s_upr <- scores_mrg_1(x, tmp_model)
    
    hessian[p, ] <- colMeans((s_upr - s_lwr) * weights) / 2e-3
  }
  hessian
}

bdiag <- function(blocks) {
  blocks <- lapply(blocks, as.matrix)
  rows <- sapply(blocks, NROW)
  cols <- sapply(blocks, NCOL)
  m <- matrix(0, sum(rows), sum(cols))
  
  rows <- c(0, cumsum(rows))
  cols <- c(0, cumsum(cols))
  for (i in seq_along(blocks)) {
    m[(rows[i] + 1):rows[i + 1], (cols[i] + 1):cols[i + 1]] <- blocks[[i]]
  }
  m
}

hessian_mrg <- function(x, model) {
  d <- NCOL(x)
  Hs <- lapply(
    seq_len(d), 
    function(m) hessian_mrg_1(x[, m], model$margins[[m]], model$copula$weights)
  )
  bdiag(Hs)
}

hessian_mxd <- function(x, model, cores = 1) {
  u <- to_unif(x, model$margins)
  npars_mrg <- sapply(model$margins, length)
  hessian <- matrix(NA, sum(npars_mrg), model$copula$npars)
  i_p <- 1
  w <- model$copula$weights
  if (!length(w))
    w <- 1
  for (m in seq_along(model$margins)) {
    for (p in seq_along(model$margins[[m]])) {
      tmp_model <- model$margins[[m]]
      tmp_u <- u
      
      tmp_model[p] <- model$margins[[m]][p] - 1e-3
      tmp_u[, m] <- univariateML::pml(x[, m], tmp_model)
      s_lwr <- svinecop_scores(tmp_u, model$copula, cores = cores)
      
      
      tmp_model[p] <- model$margins[[m]][p] + 1e-3
      tmp_u[, m] <- univariateML::pml(x[, m], tmp_model)
      s_upr <- svinecop_scores(tmp_u, model$copula, cores = cores)
      hessian[i_p, ] <- colMeans((s_upr - s_lwr) * w) / 2e-3
      i_p <- i_p + 1
    }
  }
  
  hessian
}