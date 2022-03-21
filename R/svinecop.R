#' Stationary vine copula models
#'
#' Automated fitting or creation of custom S-vine copula models
#'
#' @param data a matrix or data.frame (copula data should have approximately
#'   uniform margins).
#' @param p the Markov order.
#' @param var_types variable types; discrete variables not (yet) allowed.
#' @param out_vertices the out-vertex; if `NA`, the out-vertex is selected
#'   automatically if no structure is provided, and is equivalent to 1 if a
#'   structure is provided.
#' @param in_vertices the in-vertex; if `NA`, the in-vertex is selected
#'   automatically if no structure is provided, and is equivalent to 1 if a
#'   structure is provided.
#' @param type type of stationary vine; `"S"` (default) for general S-vines,
#'  `"D"` for Smith's long D-vine, `"M"` for Beare and Seo's M-vine.
#' @param family_set a character vector of families; see [rvinecopulib::bicop()]
#'   for additional options.
#' @param cs_structure the cross-sectional vine structure (see
#'   [rvinecopulib::rvine_structure()]; `cs_structure = NA` performs automatic
#'   structure selection.
#' @param par_method the estimation method for parametric models, either `"mle"`
#'   for sequential maximum likelihood, `"itau"` for inversion of Kendall's tau
#'   (only available for one-parameter families and `"t"`.
#' @param nonpar_method the estimation method for nonparametric models, either
#'   `"constant"` for the standard transformation estimator, or
#'   `"linear"`/`"quadratic"` for the local-likelihood approximations of order
#'   one/two.
#' @param mult multiplier for the smoothing parameters of nonparametric
#'   families. Values larger than 1 make the estimate more smooth, values less
#'   than 1 less smooth.
#' @param selcrit criterion for family selection, either `"loglik"`, `"aic"`,
#'   `"bic"`, `"mbic"`. For `vinecop()` there is the additional option
#'   `"mbicv"`.
#' @param weights optional vector of weights for each observation.
#' @param psi0 prior probability of a non-independence copula (only used for
#'   `selcrit = "mbic"` and `selcrit = "mbicv"`).
#' @param presel whether the family set should be thinned out according to
#'   symmetry characteristics of the data.
#' @param keep_data whether the data should be stored (necessary for using
#'   [fitted()]).
#' @param trunc_lvl currently unsupported.
#' @param tree_crit the criterion for tree selection, one of `"tau"`, `"rho"`,
#'   `"hoeffd"`, or `"mcor"` for Kendall's \eqn{\tau}, Spearman's \eqn{\rho},
#'   Hoeffding's \eqn{D}, and maximum correlation, respectively.
#' @param threshold for thresholded vine copulas; `NA` indicates that the
#'   threshold should be selected automatically by [rvinecopulib::mBICV()].
#' @param show_trace logical; whether a trace of the fitting progress should be
#'   printed.
#' @param cores number of cores to use; if more than 1, estimation of pair
#'   copulas within a tree is done in parallel.
#'   
#' @return Returns the fitted model as an object with classes 
#'   `svinecop` and `svinecop_dist`. Also inherits from `vinecop`, `vinecop_dist`
#'   such that many functions from [rvinecopulib] can be called.
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
#' u <- pseudo_obs(returns[1:100, 1:3]) 
#'
#' # fit parametric S-vine copula model with Markov order 1
#' fit <- svinecop(u, p = 1, family_set = "parametric")
#' fit 
#' summary(fit)
#' plot(fit)
#' contour(fit)
#' logLik(fit)
#' 
#' pairs(svinecop_sim(500, rep = 1, fit))
svinecop <- function(data, p, var_types = rep("c", NCOL(data)),
                     family_set = "all", cs_structure = NA,
                     out_vertices = NA, in_vertices = NA,
                     type = "S",
                     par_method = "mle", nonpar_method = "constant", mult = 1,
                     selcrit = "aic", weights = numeric(), psi0 = 0.9,
                     presel = TRUE, trunc_lvl = Inf, tree_crit = "tau",
                     threshold = 0, keep_data = FALSE, show_trace = FALSE,
                     cores = 1) {
  assert_that(
    is.character(family_set),
    inherits(cs_structure, "matrix") ||
      inherits(cs_structure, "rvine_structure") ||
      (is.scalar(cs_structure) && is.na(cs_structure)),
    is.character(var_types),
    is.string(par_method),
    is.string(nonpar_method),
    is.number(mult), mult > 0,
    is.string(selcrit),
    is.numeric(weights),
    is.number(psi0), psi0 > 0, psi0 < 1,
    is.flag(presel),
    is.scalar(trunc_lvl),
    is.string(tree_crit),
    is.scalar(threshold),
    is.flag(keep_data),
    is.number(cores), cores > 0,
    is.numeric(in_vertices) || all(is.na(in_vertices)),
    is.numeric(out_vertices) || all(is.na(out_vertices))
  )
  
  if (any(var_types != "c"))
    stop("discrete variables not yet implemented.")

  if ((type != "S") & (NCOL(data) > 1)) {
    if (!is.na(cs_structure))
      warning("fixed cross-sectional structure; type argument is ignored.")
    if (type == "M") {
      sel <- select_mvine(data)
    } else if (type == "D") {
      sel <- select_dvine(data)
    } else {
      stop("'type' must be one of 'S', 'M', 'D'.")
    }
    cs_structure <- sel$cs_structure
    out_vertices <- sel$out_vertices
    in_vertices  <- sel$in_vertices
  }

  # check if families known (w/ partial matching) and expand convenience defs
  family_set <- rvinecopulib::vinecop(
    0.5, family_set = family_set, par_method = par_method
  )$controls$family_set

  ## pre-process input
  is_structure_provided <- !(is.scalar(cs_structure) && is.na(cs_structure))
  if (is_structure_provided) {
    cs_structure <- as_rvine_structure(cs_structure)
  }

  ## fit and select copula model
  vinecop <- svinecop_select_cpp(
    data = as.matrix(data),
    p = p,
    var_types = var_types,
    out_vertices = out_vertices,
    in_vertices = in_vertices,
    is_structure_provided = is_structure_provided,
    structure = cs_structure,
    family_set = family_set,
    par_method = par_method,
    nonpar_method = nonpar_method,
    mult = mult,
    selection_criterion = selcrit,
    weights = weights,
    psi0 = psi0,
    preselect_families = presel,
    truncation_level = ifelse( # Inf cannot be passed to C++
      is.finite(trunc_lvl),
      trunc_lvl,
      .Machine$integer.max
    ),
    tree_criterion = tree_crit,
    threshold = threshold,
    select_truncation_level = is.na(trunc_lvl),
    select_threshold = is.na(threshold),
    show_trace = show_trace,
    num_threads = cores
  )

  ## make all pair-copulas bicop objects
  vinecop$pair_copulas <- lapply(
    vinecop$pair_copulas,
    function(tree) lapply(tree, as.bicop)
  )

  ## add information about the fit
  nms0 <- colnames(data)
  if (is.null(nms0))
    nms0 <- paste0("V", 1:NCOL(data))
  vinecop$names <- paste(rep(nms0, p + 1),
                         rep(seq_len(p + 1), each = NCOL(data)),
                         sep = "-")
  if (keep_data)
    vinecop$data <- data

  vinecop$controls <- list(
    family_set = family_set,
    par_method = par_method,
    nonpar_method = nonpar_method,
    mult = mult,
    selcrit = selcrit,
    weights = weights,
    presel = presel,
    trunc_lvl = trunc_lvl,
    tree_crit = tree_crit,
    threshold = threshold
  )
  vinecop$nobs <- NROW(data)

  structure(
    vinecop,
    class = c("svinecop", "svinecop_dist", "vinecop", "vinecop_dist")
  )
}


#' Custom S-vine models
#'
#' @param pair_copulas A nested list of 'bicop_dist' objects, where
#'   `pair_copulas[[t]][[e]]` corresponds to the pair-copula at edge `e` in tree
#'   `t`. Only the most-left unique pair copulas are used, others can be omitted.
#' @param cs_structure The cross-sectional structure. Either a matrix, or an
#'    `rvine_structure` object; see `rvinecopulib::rvine_structure()`
#' @param p the Markov order.
#' @param out_vertices the out-vertex; if `NA`, the out-vertex is selected
#'   automatically if no structure is provided, and is equivalent to 1 if a
#'    structure is provided.
#' @param in_vertices the in-vertex; if `NA`, the in-vertex is selected
#'   automatically if no structure is provided, and is equivalent to 1 if a
#'    structure is provided.
#' @param var_types variable types; discrete variables not (yet) allowed.
#' 
#' @return Returns the model as an object with classes 
#'   `svinecop_dist`. Also inherits from `vinecop_dist`
#'   such that many functions from [rvinecopulib] can be called.
#'   
#' @seealso [svinecop_loglik], [svinecop_sim], [svinecop_hessian], 
#'   [svinecop_scores]
#'
#'
#' @export
#' @importFrom assertthat is.count
#' @importFrom rvinecopulib is.rvine_structure
#' 
#' @examples 
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
svinecop_dist <- function(pair_copulas, cs_structure, p,
                          out_vertices, in_vertices,
                          var_types = rep("c", dim(cs_structure)[1])) {
  assert_that(
    is.list(pair_copulas),
    is.matrix(cs_structure) | is.rvine_structure(cs_structure),
    length(pair_copulas) == (p + 1) * dim(cs_structure)[1] - 1,
    is.count(p) || p == 0,
    all(in_vertices <= dim(cs_structure)[1]),
    all(out_vertices <= dim(cs_structure)[1])
  )
  if (any(var_types != "c"))
    stop("discrete variables not yet implemented.")
  
  d0 <- dim(cs_structure)[1]
  d <- d0 * (p + 1)
  if (p > 0) {
    for (t in seq_len(d - 1)) {
      # only use unique pair copulas
      pair_copulas[[t]] <- pair_copulas[[t]][seq_len(min(d0, d - t))]

      # replicate pair copulas for other lags
      pair_copulas[[t]] <- rep(pair_copulas[[t]], p + 1)[seq_len(d - t)]
    }
  }

  ## set up model in C++
  model <- list(
    pair_copulas = pair_copulas,
    cs_structure = as_rvine_structure(cs_structure),
    p = p,
    out_vertices = out_vertices,
    in_vertices = in_vertices,
    var_types = var_types
  )

  model <- svinecop_create_cpp(model)

  ## make all pair-copulas bicop objects
  model$pair_copulas <- lapply(
    model$pair_copulas,
    function(tree) lapply(tree, as.bicop)
  )

  ## fix class of structure objects
  class(model$structure) <- c("rvine_structure", class(model$structure))
  class(model$cs_structure) <- c("rvine_structure", class(model$cs_structure))

  structure(model, class = c("svinecop_dist", "vinecop_dist"))
}

#' @export
print.svinecop_dist <- function(x, ...) {
  cat(dim(x$cs_structure)[1],
    "-dimensional S-vine copula model of order p = ", x$p,
    " ('svinecop_dist')\n",
    sep = ""
  )
  invisible(x)
}

#' @importFrom rvinecopulib as_rvine_matrix get_structure par_to_ktau
#' @export
summary.svinecop_dist <- function(object, 
                                  trees = seq_len(dim(object)["trunc_lvl"]), 
                                  ...) {
  mat <- as_rvine_matrix(get_structure(object))
  d <- unname(dim(object)[1])
  cs_dim <- unname(dim(object$cs_structure)[1])
  
  trees <- intersect(trees, seq_len(dim(object)["trunc_lvl"]))
  n_pcs <- cs_dim^2 * object$p + choose(cs_dim, 2)
  mdf <- as.data.frame(matrix(NA, n_pcs, 10))
  names(mdf) <- c(
    "tree", "edge",
    "conditioned", "conditioning", "var_types",
    "family", "rotation", "parameters", "df", "tau"
  )
  k <- 1
  for (t in trees) {
    for (e in seq_len(min(d - t, cs_dim))) {
      mdf$tree[k] <- t
      mdf$edge[k] <- e
      mdf$conditioned[k] <- list(c(mat[d - e + 1, e], mat[t, e]))
      min_lag <- min(ceiling(mdf$conditioned[[k]] / cs_dim))
      if (all(mdf$conditioned[[k]] >= cs_dim)) {
        mdf$conditioned[[k]] <- mdf$conditioned[[k]] - (min_lag - 1) * cs_dim
      }
      mdf$conditioning[k] <- list(mat[rev(seq_len(t - 1)), e])
      pc <- object$pair_copulas[[t]][[e]]
      mdf$var_types[k] <- paste(pc$var_types, collapse = ",")
      mdf$family[k] <- pc$family
      mdf$rotation[k] <- pc$rotation
      mdf$parameters[k] <- list(pc$parameters)
      if (pc$family %in% "tll")
        mdf$parameters[k] <- list("[30x30 grid]")
      mdf$df[k] <- pc$npars
      mdf$tau[k] <- par_to_ktau(pc)
      k <- k + 1
    }
  }
  class(mdf) <- c("summary_df", class(mdf))
  mdf
}
