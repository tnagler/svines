# Pseudo-residuals of S-vine models

Pseudo-residuals are defined as the Rosenblatt transform of the data,
conditional on the past. Under a correctly specified model, they are
approximately *iid* uniform on \\\[0, 1\]^d\\.

## Usage

``` r
svine_pseudo_residuals(x, model, cores = 1)
```

## Arguments

- x:

  the data.

- model:

  model inheriting from class
  [svine_dist](https://tnagler.github.io/svines/reference/svine_dist.md).

- cores:

  number of cores to use; if larger than one, computations are performed
  in parallel on `cores` batches.

## Value

An `n`-by-`d` matrix of pseudo-residuals, where
`n = NROW(x) - model$copula$p` and `d` is the cross-sectional dimension.

## Examples

``` r
# load data set
data(returns)  

dat <- returns[1:100, 1:3]

# fit parametric S-vine model with Markov order 1
fit <- svine(dat, p = 1, family_set = "parametric")

# compute pseudo-residuals
# (should be independent uniform across variables and time)
v <- svine_pseudo_residuals(dat, fit)
pairs(cbind(v[-1, ], v[-nrow(v), ]))

```
