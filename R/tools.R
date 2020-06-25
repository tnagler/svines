#' @importFrom TSP TSP insert_dummy solve_TSP cut_tour
find_hamilton_path <- function(u) {
  d0 <- ncol(u)
  w <- 1 - abs(wdm::wdm(cbind(u[-1, ], u[-nrow(u), ]), method = "kendall"))
  tsp <- TSP(w[1:d0, 1:d0])
  hamilton <- insert_dummy(tsp, label = "cut")
  
  # withSeed to make it reproducible
  sol <- R.utils::withSeed(solve_TSP(hamilton, method = "repetitive_nn"), 5)
  as.numeric(cut_tour(sol, "cut"))
}


select_mvine <- function(u) {
  d0 <- ncol(u)
  css <- dvine_structure(find_hamilton_path(u))
  tau <- wdm::wdm(u[-nrow(u), css$order[c(1, d0)]],  # t
                  u[-1, css$order[c(1, d0)]],        # t + 1
                  method = "kendall")
  
  if (abs(tau[1, 1]) < abs(tau[2, 2]))
    css$order <- rev(css$order)
  list(
    cs_structure = css,
    out_vertices = css$order,
    in_vertices = css$order
  )
}

select_dvine <- function(u) {
  d0 <- ncol(u)
  css <- dvine_structure(find_hamilton_path(u))
  tau <- wdm::wdm(u[-nrow(u), css$order[c(1, d0)]],  # t
                  u[-1, css$order[c(1, d0)]],        # t + 1
                  method = "kendall")

  if (abs(tau[1, 2]) > abs(tau[2, 1]))
    css$order <- rev(css$order)
  list(
    cs_structure = css,
    out_vertices = rev(css$order),
    in_vertices = css$order
  )
}

