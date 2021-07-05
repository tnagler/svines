as.bicop <- function(object) {
  if (!all(c("family", "rotation", "parameters", "npars") %in% 
           names(object))) {
    stop("object cannot be coerced to class 'bicop'")
  }
  structure(object, class = c("bicop", "bicop_dist"))
}

if_vec_to_matrix <- function(u, to_col = FALSE) {
  if (is.null(u)) 
    return(NULL)
  assert_that(is.numeric(u) | is.data.frame(u))
  if (NCOL(u) == 1) {
    if (to_col) {
      u <- matrix(u, length(u), 1)
    }
    else {
      u <- matrix(u, 1, length(u))
    }
  }
  if (!is.matrix(u)) {
    u <- as.matrix(u)
  }
  u
}

process_family_set <- function(family_set, par_method) 
{
  family_set <- check_and_match_family_set(family_set)
  family_set <- expand_family_set(family_set)
  if (par_method == "itau") {
    if (any(!(family_set %in% family_set_itau))) {
      warning("Only families (", paste(family_set_itau, 
                                       collapse = ", "), ") can be used with ", "'par_method = ", 
              "\"itau\"", "'; ", "reducing family set.", call. = FALSE)
      family_set <- intersect(family_set, family_set_itau)
    }
  }
  family_set
}

expand_family_set <- function(family_set) 
{
  unique(unlist(lapply(family_set, expand_family)))
}