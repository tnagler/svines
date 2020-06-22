#' Stationary vine copula models
#'
#' Automated fitting or creation of custom S-vine copula models
#'
#' @param data a matrix or data.frame (copula data should have approximately
#'   uniform margins).
#' @param p the Markov order.
#' @param var_types variable types; discrete variables not (yet) allowed.
#' @param in_vertices the in-vertex; if `NA`, the in-vertex is selected
#'   automatically if no structure is provided, and is equivalent to 1 if a
#'   structure is provided.
#' @param out_vertices the out-vertex; if `NA`, the out-vertex is selected
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
#' @importFrom assertthat assert_that is.scalar is.string is.number is.flag
#' @importFrom rvinecopulib as_rvine_structure
#' @export
svinecop <- function(data, p, var_types = rep("c", NCOL(data)),
                     family_set = "all", cs_structure = NA,
                     in_vertices = NA, out_vertices = NA,
                     type = "S",
                     par_method = "mle", nonpar_method = "constant", mult = 1,
                     selcrit = "bic", weights = numeric(), psi0 = 0.9,
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
    in_vertices  <- sel$in_vertices
    out_vertices <- sel$out_vertices
  }



  # check if families known (w/ partial matching) and expand convenience defs
  family_set <- rvinecopulib:::process_family_set(family_set, par_method)

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
    in_vertices = in_vertices,
    out_vertices = out_vertices,
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
    function(tree) lapply(tree, rvinecopulib:::as.bicop)
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
#'   `t`. Only the most-left unique pair copulas are used, others can be ommitted.
#' @param cs_structure The cross-strectional structure. Either a matrix, or an
#'    `rvine_structure` object; see `rvinecopulib::rvine_structure()`
#' @param p the Markov order.
#' @param in_vertices the in-vertex; if `NA`, the in-vertex is selected
#'   automatically if no structure is provided, and is equivalent to 1 if a
#'    structure is provided.
#' @param out_vertices the out-vertex; if `NA`, the out-vertex is selected
#'   automatically if no structure is provided, and is equivalent to 1 if a
#'    structure is provided.
#'
#' @export
#' @importFrom assertthat is.count
#' @importFrom rvinecopulib is.rvine_structure
svinecop_dist <- function(pair_copulas, cs_structure, p,
                          in_vertices, out_vertices,
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
      pair_copulas[[t]] <- c(pair_copulas[[t]], rep(pair_copulas[[t]], p))

      # remove superfluous pcs
      pair_copulas[[t]] <- pair_copulas[[t]][seq_len(d - t)]
    }
  }

  ## set up model in C++
  model <- list(
    pair_copulas = pair_copulas,
    cs_structure = as_rvine_structure(cs_structure),
    p = p,
    in_vertices = in_vertices,
    out_vertices = out_vertices,
    var_types = var_types
  )

  model <- svinecop_create_cpp(model)

  ## make all pair-copulas bicop objects
  model$pair_copulas <- lapply(
    model$pair_copulas,
    function(tree) lapply(tree, rvinecopulib:::as.bicop)
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
    " ('svinecop_dist')",
    sep = ""
  )
  invisible(x)
}


#' Simulate from a S-vine model
#'
#' [svinecop_sim()] simulates a time series of length `n`;
#' [svinecop_sim_conditional()] simulates `n` observations for the next time
#' point conditional on past observations;
#' [svinecop_sim_ahead()]
#' simulates the next `n` time points conditional on the past.
#'
#' @aliases svinecop_sim svinecop_sim_conditional svinecop_sim_ahead
#'
#' @param n number of observations.
#' @param svinecop a S-vine model object.
#' @param qrng if `TRUE`, generates quasi-random numbers using the multivariate
#' Generalized Halton sequence up to dimension 300 and the Generalized Sobol
#' sequence in higher dimensions (default `qrng = FALSE`).
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#'
#' @export
#'
svinecop_sim <- function(n, svinecop, qrng = FALSE) {
  assert_that(
    is.number(n),
    inherits(svinecop, "svinecop_dist"),
    is.flag(qrng)
  )

  U <- svinecop_sim_cpp(svinecop, n, qrng, rvinecopulib:::get_seeds())
  if (!is.null(svinecop$names)) {
    colnames(U) <- simplify_names(svinecop)
  }

  U
}

simplify_names <- function(svinecop) {
  nms <- svinecop$names[seq_along(svinecop$in_vertices)]
  nms <- strsplit(nms, "-")
  sapply(nms, function(n) paste(n[-length(n)], collapse = ""))
}

#' @rdname svinecop_sim
#' @param data time series of past observations (more recent observations are
#'    at the bottom of the matrix.)
#' @export
svinecop_sim_conditional <- function(n, data, svinecop, qrng = FALSE, cores = 1) {
  assert_that(
    is.number(n),
    inherits(svinecop, "svinecop_dist"),
    is.flag(qrng),
    is.number(cores)
  )
  data <- rvinecopulib:::if_vec_to_matrix(data,
                                          dim(svinecop$cs_structure)[1] == 1)

  U <- svinecop_sim_conditional_cpp(
    svinecop, n, data, qrng, cores,
    rvinecopulib:::get_seeds()
  )
  if (!is.null(svinecop$names)) {
    colnames(U) <- simplify_names(svinecop)
  }

  U
}

#' @rdname svinecop_sim
#' @export
svinecop_sim_ahead <- function(n, data, svinecop, qrng = FALSE) {
  assert_that(
    is.number(n),
    inherits(svinecop, "svinecop_dist"),
    is.flag(qrng)
  )
  data <- rvinecopulib:::if_vec_to_matrix(data, dim(svinecop$cs_structure)[1] == 1)

  U <- svinecop_sim_ahead_cpp(
    svinecop, n, data, qrng,
    rvinecopulib:::get_seeds()
  )
  if (!is.null(svinecop$names)) {
    colnames(U) <- simplify_names(svinecop)
  }

  U
}

#' @export
svinecop_loglik <- function(u, svinecop, cores = 1) {
  assert_that(inherits(svinecop, "svinecop_dist"))
  u <- rvinecopulib:::if_vec_to_matrix(u, dim(svinecop$cs_structure)[1] == 1)
  tryCatch(svinecop_loglik_cpp(u, svinecop, cores), error = function(e) browser())
}

# svinecop_cond_cdf <- function(u, conditioned, svinecop, cores = 1) {
#   assert_that(
#     is.count(conditioned),
#     inherits(svinecop, "svinecop_dist"),
#     conditioned <= dim(svinecop$cs_structure)[1]
#   )
#   u <- rvinecopulib:::if_vec_to_matrix(u, dim(svinecop$cs_structure)[1] == 1)
#   svinecop_cond_cdf_cpp(u, conditioned - 1, svinecop, cores)
# }