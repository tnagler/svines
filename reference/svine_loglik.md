# Log-likelihood for S-vine models

Log-likelihood for S-vine models

## Usage

``` r
svine_loglik(x, model, cores = 1)
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

The log-likelihood of the data under the model.

## Examples

``` r
# load data set
data(returns)  

# fit parametric S-vine model with Markov order 1
fit <- svine(returns[1:100, 1:3], p = 1, family_set = "parametric")

svine_loglik(returns[1:100, 1:3], fit)
#> [1] 942.4027
```
