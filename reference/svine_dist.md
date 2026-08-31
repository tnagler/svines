# Custom S-vine distribution models

Custom S-vine distribution models

## Usage

``` r
svine_dist(margins, copula)
```

## Arguments

- margins:

  A list with one marginal distribution per variable. Each element is
  either a `univariateML` object or a list with callables `$p` (CDF),
  `$q` (quantile function), and for discrete margins `$p_sub`
  (left-limit CDF, i.e. F(x-)). Discrete list-form margins must also
  carry `attr(., "type") = "discrete"` and `class` including
  `"svine_margin"`.

- copula:

  the copula model; an object of class `svinecop_dist` with
  cross-sectional dimension `d`.

## Value

Returns the model as an object with class `svine_dist`. A list with
entries

- `$margins`: list of marginal models from
  [univariateML::univariateML_models](https://jonasmoss.github.io/univariateML/reference/univariateML_models.html),

- `$copula`: an object of `svinecop_dist`.

## See also

svine_dist,
[svine_loglik](https://tnagler.github.io/svines/reference/svine_loglik.md),
[svine_sim](https://tnagler.github.io/svines/reference/svine_sim.md),
[svine_bootstrap_models](https://tnagler.github.io/svines/reference/svine_bootstrap_models.md)

## Examples

``` r
## marginal objects
# create dummy univariateML models
univ1 <- univ2 <- univariateML::mlnorm(rnorm(10))

# modify the parameters to N(5, 10) and N(0, 2) distributions
univ1[] <- c(5, 10)
univ2[] <- c(0, 2)

## copula object
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
    
model <- svine_dist(margins = list(univ1, univ2), copula = cop)
summary(model)
#> $margins
#> # A data.frame: 2 x 5 
#>  margin name  model parameters loglik
#>       1   V1 Normal      5, 10    -16
#>       2   V2 Normal       0, 2    -16
#> 
#> $copula
#> # A data.frame: 5 x 10 
#>  tree edge conditioned conditioning var_types   family rotation parameters df
#>     1    1        3, 1                    c,c gaussian        0       -0.1  1
#>     1    2        2, 1                    c,c  clayton        0          3  1
#>     2    1        4, 1            3       c,c gaussian        0        0.2  1
#>     2    2        3, 2            1       c,c    indep        0             0
#>     3    1        4, 2         1, 3       c,c    indep        0             0
#>     tau
#>  -0.064
#>   0.600
#>   0.128
#>   0.000
#>   0.000
#> 
```
