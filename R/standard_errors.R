#' Asymptotic variance of parameter estimates of an S-vine model
#'
#' @param x the data.
#' @param model S-vine model (inheriting from [svine_dist]).
#'
#' @return A k-by-k matrix, where k is the total number of parameters in the 
#'   model. Parameters are ordered as follows: 
#'   marginal parameters, copula parameters of first tree, copula parameters of 
#'   second tree, etc. Duplicated parameters in the copula model are omitted.
#' @export
#' @examples
#' # load data set
#' data(returns)  
#' dat <- returns[1:100, 1:2]
#' 
#' # fit parametric S-vine model with Markov order 1
#' fit <- svine(x, p = 1, family_set = "parametric")
#' svine_avar(u, fit)
svine_avar <- function(x, model) {
  x <- as.matrix(x)
  S_mrg <- scores_mrg(x, model)
  H_mrg <- hessian_mrg(x, model)
  H_mxd <- hessian_mxd(x, model)
  
  u <- as.matrix(to_unif(x, model$margins))
  S_cop <- svinecop_scores(u, model$copula)
  H_cop <- svinecop_hessian(u, model$copula)
  
  I <- cov(cbind(S_mrg[-seq_len(model$copula$p), ], S_cop))
  H <- rbind(
    cbind(H_mrg, H_mxd), 
    cbind(matrix(0, ncol(H_mxd), nrow(H_mxd)), H_cop)
  )
  Hi <- solve(H)
  Hi %*% I %*% t(Hi) / nrow(u)
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

hessian_mrg_1 <- function(x, model) {
  npars <- length(model)
  hessian <- matrix(NA, npars, npars)
  for (p in seq_len(npars)) {
    tmp_model <- model
    
    tmp_model[p] <- model[p] - 1e-2
    s_lwr <- scores_mrg_1(x, tmp_model)
    
    tmp_model[p] <- model[p] + 1e-2
    s_upr <- scores_mrg_1(x, tmp_model)
    
    hessian[p, ] <- colMeans(s_upr - s_lwr) / 2e-2
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
    m[(rows[i] + 1):rows[i + 1], (cols[i] + 1):cols[i +  1]] <- blocks[[i]]
  }
  m
}

hessian_mrg <- function(x, model) {
  d <- NCOL(x)
  Hs <- lapply(seq_len(d), function(m) hessian_mrg_1(x[, m], model$margins[[m]]))
  bdiag(Hs)
}

hessian_mxd <- function(x, model) {
  u <- to_unif(x, model$margins)
  npars_mrg <- sapply(model$margins, length)
  hessian <- matrix(NA, sum(npars_mrg), model$copula$npars)
  i_p <- 1
  for (m in seq_along(model$margins)) {
    for (p in seq_along(model$margins[[m]])) {
      tmp_model <- model$margins[[m]]
      tmp_u <- u
      
      tmp_model[p] <- model$margins[[m]][p] - 1e-3
      tmp_u[, m] <- univariateML::pml(x[, m], tmp_model)
      s_lwr <- svinecop_scores(tmp_u, model$copula)
      
      tmp_model[p] <- model$margins[[m]][p] + 1e-3
      tmp_u[, m] <- univariateML::pml(x[, m], tmp_model)
      s_upr <- svinecop_scores(tmp_u, model$copula)
      
      hessian[i_p, ] <- colMeans(s_upr - s_lwr) / 2e-3
      i_p <- i_p + 1
    }
  }
  
  hessian
}
