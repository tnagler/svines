# Getting started with svines

``` r

library(svines)
```

## Model layers

An S-vine model combines marginal distributions with a stationary vine
copula. The package exposes these two layers separately:

- [`svine()`](https://tnagler.github.io/svines/reference/svine.md) fits
  a complete distribution model to observed data. It estimates the
  marginal distributions, transforms the observations to the unit
  hypercube, and fits the copula.
- [`svinecop()`](https://tnagler.github.io/svines/reference/svinecop.md)
  fits only the copula and therefore expects approximately uniform
  pseudo-observations.
- [`svine_dist()`](https://tnagler.github.io/svines/reference/svine_dist.md)
  and
  [`svinecop_dist()`](https://tnagler.github.io/svines/reference/svinecop_dist.md)
  construct models from specified margins, pair copulas, and an S-vine
  structure.

The argument `p` is the Markov order. An order-one model relates the
current observation to the previous observation, while larger values
include additional lags.

## Fitting a continuous model

The `returns` data contain daily log returns of 20 companies. We use two
series and restrict the candidate families to keep this example short.

``` r

data(returns)
x <- returns[1:200, 1:2]

fit <- svine(
  x,
  p = 1,
  margin_families = c("norm", "std"),
  family_set = c("gaussian", "t")
)
fit
#> 2-dimensional S-vine distribution model of order p = 1 ('svine_dist')
summary(fit)
#> $margins
#> # A data.frame: 2 x 5 
#>  margin    name     model             parameters loglik
#>       1 Allianz    Normal       0.00034, 0.01536    551
#>       2     AXA Student-t 0.0019, 0.0197, 4.0004    522
#> 
#> $copula
#> # A data.frame: 5 x 10 
#>  tree edge conditioned conditioning var_types   family rotation    parameters
#>     1    1        3, 1                    c,c        t        0 -0.084, 4.781
#>     1    2        2, 1                    c,c        t        0    0.84, 4.87
#>     2    1        4, 1            3       c,c gaussian        0        -0.064
#>     2    2        3, 2            1       c,c gaussian        0         0.074
#>     3    1        4, 2         1, 3       c,c gaussian        0        -0.072
#>  df    tau
#>   2 -0.054
#>   2  0.637
#>   1 -0.041
#>   1  0.047
#>   1 -0.046
```

The fitted object contains the marginal models in `fit$margins` and the
copula model in `fit$copula`. Standard `rvinecopulib` methods can be
applied to the copula component.

## Simulation and diagnostics

Without a conditioning history,
[`svine_sim()`](https://tnagler.github.io/svines/reference/svine_sim.md)
generates a new stationary time series. Supplying `past` instead
generates paths conditional on the observed history.

``` r

sim <- svine_sim(n = 100, rep = 1, model = fit)
dim(sim)
#> [1] 100   2

next_obs <- svine_sim(n = 1, rep = 100, model = fit, past = x)
dim(next_obs)
#>     dim     
#>   1   2 100
```

Pseudo-residuals are conditional Rosenblatt transforms. For a fitted
model of order `p`, the result has `NROW(x) - p` rows.

``` r

residuals <- svine_pseudo_residuals(x, fit)
dim(residuals)
#> [1] 199   2
```

## Discrete variables

For discrete variables, specify `var_types = "d"` and restrict
`margin_families` to suitable discrete families. The following model
uses two Poisson margins.

``` r

counts <- cbind(
  claims = rpois(250, lambda = 2),
  events = rpois(250, lambda = 4)
)

fit_discrete <- svine(
  counts,
  p = 1,
  var_types = c("d", "d"),
  margin_families = "pois",
  family_set = "gaussian"
)
fit_discrete
#> 2-dimensional S-vine distribution model of order p = 1 ('svine_dist')
svine_sim(5, rep = 1, model = fit_discrete)
#>      claims events
#> [1,]      2      5
#> [2,]      0      1
#> [3,]      5      5
#> [4,]      6      7
#> [5,]      2      3
```

[`svine()`](https://tnagler.github.io/svines/reference/svine.md)
evaluates both the CDF, `F(x)`, and its left limit, `F(x-)`, and
constructs the copula data automatically. When calling
[`svinecop()`](https://tnagler.github.io/svines/reference/svinecop.md)
directly, supply the regular CDF columns first, followed by one
left-limit column for each discrete variable.

## Copula-only models

When the marginal transformation is handled separately, fit the copula
layer directly.

``` r

u <- pseudo_obs(x)
copula_fit <- svinecop(
  u,
  p = 1,
  family_set = c("gaussian", "t")
)
copula_fit
#> 2-dimensional S-vine copula model of order p = 1 ('svinecop_dist')
```

Use
[`svinecop_loglik()`](https://tnagler.github.io/svines/reference/svinecop_loglik.md),
[`svinecop_scores()`](https://tnagler.github.io/svines/reference/svinecop_scores.md),
and
[`svinecop_hessian()`](https://tnagler.github.io/svines/reference/svinecop_hessian.md)
for copula-level inference. The corresponding `svine_*` functions
include the marginal parameters.
