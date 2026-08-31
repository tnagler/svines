# Custom S-vine models

Custom S-vine models

## Usage

``` r
svinecop_dist(
  pair_copulas,
  cs_structure,
  p,
  out_vertices,
  in_vertices,
  var_types = rep("c", dim(cs_structure)[1])
)
```

## Arguments

- pair_copulas:

  A nested list of 'bicop_dist' objects, where `pair_copulas[[t]][[e]]`
  corresponds to the pair-copula at edge `e` in tree `t`. Only the
  most-left unique pair copulas are used, others can be omitted.

- cs_structure:

  The cross-sectional structure. Either a matrix, or an
  `rvine_structure` object; see
  [`rvinecopulib::rvine_structure()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.html)

- p:

  the Markov order.

- out_vertices:

  a permutation of `1, ..., d` specifying the out-vertices.

- in_vertices:

  a permutation of `1, ..., d` specifying the in-vertices.

- var_types:

  variable types; a character vector with one entry per variable: `"c"`
  for continuous, `"d"` for discrete. For discrete variables, methods on
  the returned object expect data in the compact layout: all regular CDF
  `F(x)` columns first, then one left-limit column `F(x-)` per discrete
  variable.

## Value

Returns the model as an object with classes `svinecop_dist`. Also
inherits from `vinecop_dist` such that many functions from
`rvinecopulib` can be called.

## See also

[svinecop_loglik](https://tnagler.github.io/svines/reference/svinecop_loglik.md),
[svinecop_sim](https://tnagler.github.io/svines/reference/svinecop_sim.md),
[svinecop_hessian](https://tnagler.github.io/svines/reference/svinecop_hessian.md),
[svinecop_scores](https://tnagler.github.io/svines/reference/svinecop_scores.md)

## Examples

``` r
cs_struct <- cvine_structure(1:2)
pcs <- list(
  list(  # first tree
    bicop_dist("clayton", 0, 3), # cross sectional copula
    bicop_dist("gaussian", 0, -0.1)  # serial copula
  ),
  list(  # second tree
    bicop_dist("gaussian", 0, 0.2), bicop_dist("indep")  
  ),
  list( # third tree
    bicop_dist("indep")
  )
)

cop <- svinecop_dist(
  pcs, cs_struct, p = 1, out_vertices = 1:2, in_vertices = 1:2)
```
