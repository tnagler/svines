# Pseudo-residuals of S-vine copula models

Pseudo-residuals are defined as the Rosenblatt transform of the data,
conditional on the past. Under a correctly specified model, they are
approximately *iid* uniform on \\\[0, 1\]^d\\.

## Usage

``` r
svinecop_pseudo_residuals(u, model, cores = 1)
```

## Arguments

- u:

  the data; should have approximately uniform margins.

- model:

  model inheriting from class
  [svinecop_dist](https://tnagler.github.io/svines/reference/svinecop_dist.md).

- cores:

  number of cores to use; if larger than one, computations are performed
  in parallel on `cores` batches.

## Value

An `n`-by-`d` matrix of pseudo-residuals, where `n = NROW(u) - model$p`
and `d` is the cross-sectional dimension.

## Examples

``` r
# load data set
data(returns)  

# convert to pseudo observations with empirical cdf for marginal distributions
u <- pseudo_obs(returns[1:100, 1:3]) 

# fit parametric S-vine copula model with Markov order 1
fit <- svinecop(u, p = 1, family_set = "parametric")

# compute pseudo-residuals
# (should be independent uniform across variables and time)
v <- svinecop_pseudo_residuals(u, fit)
pairs(cbind(v[-1, ], v[-nrow(v), ]))

```
