#' @importFrom TSP TSP insert_dummy solve_TSP cut_tour
find_hamilton_path <- function(u) {
  d0 <- ncol(u)
  w <- 1 - abs(cor(cbind(u[-1, ], u[-nrow(u), ]), method = "kendall"))
  tsp <- TSP(w[1:d0, 1:d0])
  hamilton <- insert_dummy(tsp, label = "cut")
  sol <- solve_TSP(hamilton, method = "repetitive_nn")
  as.numeric(cut_tour(sol, "cut"))
}


select_mvine <- function(u) {
  d0 <- ncol(u)
  css <- dvine_structure(find_hamilton_path(u))
  w <- 1 - abs(cor(cbind(u[-1, ], u[-nrow(u), ]), method = "kendall"))

  if (w[css$order[d0], d0 + css$order[d0]] < w[css$order[1], d0 + css$order[1]])
    css$order <- rev(css$order)
  list(
    cs_structure = css,
    in_vertices = css$order,
    out_vertices = css$order
  )
}

select_dvine <- function(u) {
  d0 <- ncol(u)
  css <- dvine_structure(find_hamilton_path(u))
  w <- 1 - abs(cor(cbind(u[-1, ], u[-nrow(u), ]), method = "kendall"))
  if (w[css$order[1], d0 + css$order[d0]] < w[css$order[d0], d0 + css$order[1]])
    css$order <- rev(css$order)
  list(
    cs_structure = css,
    in_vertices = css$order,
    out_vertices = rev(css$order)
  )
}

