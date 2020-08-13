
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R build
status](https://github.com/tnagler/svines/workflows/R-CMD-check/badge.svg)](https://github.com/tnagler/svines/actions)
[![Codecov test
coverage](https://codecov.io/gh/tnagler/svines/branch/master/graph/badge.svg)](https://codecov.io/gh/tnagler/svines?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/svines)](https://CRAN.R-project.org/package=svines)
<!-- badges: end -->

An R package that provides functionality to fit and simulate from
[stationary vine copula models for time
series](https://github.com/tnagler/tnagler.github.io/raw/master/svines.pdf).

The package is build on top of
[rvinecopulib](https://github.com/vinecopulib/rvinecopulib) and
[univariateML](https://github.com/JonasMoss/univariateML).

## Installation

Install the development version from Github.

``` r
# install.packages("devtools")
devtools::install_github("tnagler/svines")
```

## Usage

For detailed documentation and examples, see the [API
documentation](https://tnagler.github.io/svines/).

``` r
library(svines)
#> Loading required package: rvinecopulib
data(returns)  # data set of stock returns
returns <- returns[1:500, 1:3]
```

### Fitting models

``` r
fit <- svine(returns, p = 1)  # Markov order 1
#> Warning in f(x, na.rm = na.rm): The iteration limit (iterlim = 100) was reached
#> before the relative tolerance requirement (rel.tol = 0.0001220703125).

#> Warning in f(x, na.rm = na.rm): The iteration limit (iterlim = 100) was reached
#> before the relative tolerance requirement (rel.tol = 0.0001220703125).

#> Warning in f(x, na.rm = na.rm): The iteration limit (iterlim = 100) was reached
#> before the relative tolerance requirement (rel.tol = 0.0001220703125).
summary(fit)
#> $margins
#> # A data.frame: 3 x 5 
#>  margin     name          model                         parameters loglik
#>       1  Allianz Skew Student-t 0.00039, 0.01589, 5.45533, 0.91785   1382
#>       2      AXA Skew Student-t 0.00052, 0.02089, 4.35198, 0.90611   1260
#>       3 Generali      Student-t         -0.00043, 0.02090, 4.00012   1268
#> 
#> $copula
#> # A data.frame: 15 x 11 
#>  tree edge conditioned conditioning var_types family rotation      parameters
#>     1    1        4, 5                    c,c      t        0      0.86, 3.48
#>     1    2        5, 6                    c,c    bb1      180      0.41, 2.13
#>     1    3        6, 1                    c,c    joe      180             1.1
#>     1    4        1, 2                    c,c      t        0      0.86, 3.48
#>     1    5        2, 3                    c,c    bb1      180      0.41, 2.13
#>     2    1        4, 6            5       c,c  frank        0             1.1
#>     2    2        5, 1            6       c,c gumbel      180               1
#>     2    3        6, 2            1       c,c    joe       90               1
#>     2    4        1, 3            2       c,c  frank        0             1.1
#>     3    1        4, 1         6, 5       c,c      t        0 -0.0076, 9.6999
#>  df     tau loglik
#>   2  0.6624    NaN
#>   2  0.6097    NaN
#>   1  0.0506    NaN
#>   2  0.6624    NaN
#>   2  0.6097    NaN
#>   1  0.0420    NaN
#>   1  0.0194    NaN
#>   1 -0.0154    NaN
#>   1  0.0420    NaN
#>   2 -0.0048    NaN
#> # ... with 5 more rows
```

``` r
contour(fit$copula)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

### Simulation

`svine_sim()` can be used in two different ways:

#### Generate a new time series of length 500

``` r
sim <- svine_sim(n = 500, rep = 1, model = fit)
pairs(sim)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

``` r
pairs(returns)
```

<img src="man/figures/README-unnamed-chunk-5-2.png" width="100%" />

#### Sample conditionally on the past given the past

``` r
sim <- svine_sim(n = 1, rep = 100, model = fit, past = returns)
pairs(t(sim[1, , ]))
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

### Standard errors

``` r
# compute standard errors for each parameter
se <- sqrt(diag(svine_avar(returns, fit)))
```

The standard errors correspond to (in this order): parameters of first
margin, parameters of second margin, …., copula parameters in first tree
(excluding duplicates), …

The asymptotic variance is the basis of the bootstrap procedure in the
paper. To generate new models from the asymptotic distribution, use

``` r
models <- svine_sim_se_models(2, fit)
summary(models[[1]])
#> $margins
#> # A data.frame: 3 x 5 
#>  margin     name          model                          parameters loglik
#>       1  Allianz Skew Student-t     -0.0011, 0.0173, 2.9614, 0.9345     NA
#>       2      AXA Skew Student-t -0.00061, 0.02205, 2.70703, 0.94818     NA
#>       3 Generali      Student-t             -0.0023, 0.0214, 3.6512     NA
#> 
#> $copula
#> # A data.frame: 15 x 10 
#>  tree edge conditioned conditioning var_types family rotation     parameters df
#>     1    1        4, 5                    c,c      t        0     0.86, 3.55  2
#>     1    2        5, 6                    c,c    bb1      180     0.39, 2.15  2
#>     1    3        6, 1                    c,c    joe      180            1.1  1
#>     1    4        1, 2                    c,c      t        0     0.86, 3.55  2
#>     1    5        2, 3                    c,c    bb1      180     0.39, 2.15  2
#>     2    1        4, 6            5       c,c  frank        0           0.93  1
#>     2    2        5, 1            6       c,c gumbel      180              1  1
#>     2    3        6, 2            1       c,c    joe       90              1  1
#>     2    4        1, 3            2       c,c  frank        0           0.93  1
#>     3    1        4, 1         6, 5       c,c      t        0 -0.071, 11.814  2
#>      tau
#>   0.6589
#>   0.6099
#>   0.0341
#>   0.6589
#>   0.6099
#>   0.0355
#>   0.0057
#>  -0.0235
#>   0.0355
#>  -0.0452
#> # ... with 5 more rows
summary(models[[1]])
#> $margins
#> # A data.frame: 3 x 5 
#>  margin     name          model                          parameters loglik
#>       1  Allianz Skew Student-t     -0.0011, 0.0173, 2.9614, 0.9345     NA
#>       2      AXA Skew Student-t -0.00061, 0.02205, 2.70703, 0.94818     NA
#>       3 Generali      Student-t             -0.0023, 0.0214, 3.6512     NA
#> 
#> $copula
#> # A data.frame: 15 x 10 
#>  tree edge conditioned conditioning var_types family rotation     parameters df
#>     1    1        4, 5                    c,c      t        0     0.86, 3.55  2
#>     1    2        5, 6                    c,c    bb1      180     0.39, 2.15  2
#>     1    3        6, 1                    c,c    joe      180            1.1  1
#>     1    4        1, 2                    c,c      t        0     0.86, 3.55  2
#>     1    5        2, 3                    c,c    bb1      180     0.39, 2.15  2
#>     2    1        4, 6            5       c,c  frank        0           0.93  1
#>     2    2        5, 1            6       c,c gumbel      180              1  1
#>     2    3        6, 2            1       c,c    joe       90              1  1
#>     2    4        1, 3            2       c,c  frank        0           0.93  1
#>     3    1        4, 1         6, 5       c,c      t        0 -0.071, 11.814  2
#>      tau
#>   0.6589
#>   0.6099
#>   0.0341
#>   0.6589
#>   0.6099
#>   0.0355
#>   0.0057
#>  -0.0235
#>   0.0355
#>  -0.0452
#> # ... with 5 more rows
```
