# Expected Hessian of a parametric S-vine model

Expected Hessian of a parametric S-vine model

## Usage

``` r
svine_hessian(x, model, cores = 1)
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

A `k`-by-`k` estimate of the expected Hessian, where `k` is the total
number of model parameters. The estimate averages the per-observation
Hessian contributions. Parameters are ordered as follows: marginal
parameters, copula parameters of first tree, copula parameters of second
tree, etc. Duplicated parameters in the copula model are omitted.

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
#>  [1,]  1.663430e-06 -7.345804e-08  1.263944e-06 -1.929607e-07 -1.180159e-05
#>  [2,] -7.345804e-08  2.770406e-07 -9.746338e-08  1.595365e-07  1.513773e-06
#>  [3,]  1.263944e-06 -9.746338e-08  2.336204e-06 -2.406471e-07  2.508189e-05
#>  [4,] -1.929607e-07  1.595365e-07 -2.406471e-07  3.838785e-07 -7.006600e-06
#>  [5,] -1.180159e-05  1.513773e-06  2.508189e-05 -7.006600e-06  2.930565e-03
#>  [6,]  7.395517e-05  1.476121e-04  2.485149e-04  1.081140e-04  2.803207e-02
#>  [7,]  2.223717e-05  3.223243e-08  2.475950e-06 -4.471535e-08  5.204546e-04
#>  [8,] -5.252911e-04  2.693218e-04 -5.935841e-04 -4.655831e-06 -1.578136e-02
#>  [9,]  2.030250e-05  6.309622e-07  1.169314e-05  1.269796e-05 -1.280164e-03
#>                [,6]          [,7]          [,8]          [,9]
#>  [1,]  7.395517e-05  2.223717e-05 -5.252911e-04  2.030250e-05
#>  [2,]  1.476121e-04  3.223243e-08  2.693218e-04  6.309622e-07
#>  [3,]  2.485149e-04  2.475950e-06 -5.935841e-04  1.169314e-05
#>  [4,]  1.081140e-04 -4.471535e-08 -4.655831e-06  1.269796e-05
#>  [5,]  2.803207e-02  5.204546e-04 -1.578136e-02 -1.280164e-03
#>  [6,]  2.523648e+00  2.532996e-02 -5.715132e-02  2.020138e-02
#>  [7,]  2.532996e-02  1.294684e-02 -1.089807e-01  8.937678e-04
#>  [8,] -5.715132e-02 -1.089807e-01  6.862184e+00 -1.027884e-02
#>  [9,]  2.020138e-02  8.937678e-04 -1.027884e-02  9.113254e-03
```
