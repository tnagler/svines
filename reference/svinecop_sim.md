# Simulate from a S-vine copula model

Simulate from a S-vine copula model

## Usage

``` r
svinecop_sim(n, rep, model, past = NULL, qrng = FALSE, cores = 1)
```

## Arguments

- n:

  how many steps of the time series to simulate.

- rep:

  number of replications; `rep` time series of length `n` are generated.

- model:

  a S-vine copula model object (inheriting from
  [svinecop_dist](https://tnagler.github.io/svines/reference/svinecop_dist.md)).

- past:

  (optional) matrix of past observations. If provided, time series are
  simulated conditional on the past.

- qrng:

  if `TRUE`, generates quasi-random numbers using the multivariate
  Generalized Halton sequence up to dimension 300 and the Generalized
  Sobol sequence in higher dimensions (default `qrng = FALSE`).

- cores:

  number of cores to use; if larger than one, computations are performed
  in parallel over replications.

## Value

An `n`-by-`d`-by-`rep` array, where `d` is the cross-sectional dimension
of the model. This reduces to an `n`-by-`d` matrix if `rep == 1`.

## Examples

``` r
# load data set
data(returns)  

# convert to uniform margins
u <- pseudo_obs(returns[1:100, 1:3])

# fit parametric S-vine copula model with Markov order 1
fit <- svinecop(u, p = 1, family_set = "parametric")

pairs(u)   # original data

pairs(svinecop_sim(100, rep = 1, model = fit))   # simulated data


# simulate the next day conditionally on the past 500 times
pairs(t(svinecop_sim(1, rep = 100, model = fit, past = u)[1, , ]))
```
