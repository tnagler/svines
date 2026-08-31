# Score function of parametric S-vine models

Score function of parametric S-vine models

## Usage

``` r
svine_scores(x, model, cores = 1)
```

## Arguments

- x:

  the data.

- model:

  S-vine model (inheriting from
  [svine_dist](https://tnagler.github.io/svines/reference/svine_dist.md)).

- cores:

  number of cores to use.

## Value

An `n`-by-`k` matrix, where `n = NROW(x) - model$copula$p` and `k` is
the total number of model parameters. The rows correspond to the
effective time points after the first `p` observations. Parameters are
ordered as follows: marginal parameters, copula parameters of first
tree, copula parameters of second tree, etc. Duplicated parameters in
the copula model are omitted.

## Examples

``` r
data(returns)
dat <- returns[1:100, 1:2]

# fit parametric S-vine model with Markov order 1
model <- svine(dat, p = 1, family_set = "parametric")

# Implementation of asymptotic variances
I <- cov(svine_scores(dat, model))
H <- svine_hessian(dat, model)
Hi <- solve(H)
n_eff <- nrow(dat) - model$copula$p
Hi %*% I %*% t(Hi) / n_eff
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  1.663430e-06 -7.345803e-08  1.263944e-06 -1.929607e-07 -1.180159e-05
#>  [2,] -7.345803e-08  2.770406e-07 -9.746336e-08  1.595365e-07  1.513774e-06
#>  [3,]  1.263944e-06 -9.746336e-08  2.336204e-06 -2.406471e-07  2.508189e-05
#>  [4,] -1.929607e-07  1.595365e-07 -2.406471e-07  3.838785e-07 -7.006600e-06
#>  [5,] -1.180159e-05  1.513774e-06  2.508189e-05 -7.006600e-06  2.930565e-03
#>  [6,]  7.395518e-05  1.476121e-04  2.485149e-04  1.081140e-04  2.803207e-02
#>  [7,]  2.223716e-05  3.223684e-08  2.475948e-06 -4.471570e-08  5.204545e-04
#>  [8,] -5.252907e-04  2.693214e-04 -5.935837e-04 -4.655805e-06 -1.578135e-02
#>  [9,]  2.030250e-05  6.309618e-07  1.169314e-05  1.269796e-05 -1.280164e-03
#>                [,6]          [,7]          [,8]          [,9]
#>  [1,]  7.395518e-05  2.223716e-05 -5.252907e-04  2.030250e-05
#>  [2,]  1.476121e-04  3.223684e-08  2.693214e-04  6.309618e-07
#>  [3,]  2.485149e-04  2.475948e-06 -5.935837e-04  1.169314e-05
#>  [4,]  1.081140e-04 -4.471570e-08 -4.655805e-06  1.269796e-05
#>  [5,]  2.803207e-02  5.204545e-04 -1.578135e-02 -1.280164e-03
#>  [6,]  2.523648e+00  2.532996e-02 -5.715128e-02  2.020138e-02
#>  [7,]  2.532996e-02  1.294684e-02 -1.089806e-01  8.937680e-04
#>  [8,] -5.715128e-02 -1.089806e-01  6.862175e+00 -1.027885e-02
#>  [9,]  2.020138e-02  8.937680e-04 -1.027885e-02  9.113254e-03
```
