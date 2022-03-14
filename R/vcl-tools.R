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
