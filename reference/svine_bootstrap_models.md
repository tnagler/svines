# Bootstrap S-vine models

Computes bootstrap replicates of a given model using the one-step block
multiplier bootstrap of Nagler et al. (2022).

## Usage

``` r
svine_bootstrap_models(n_models, model)
```

## Arguments

- n_models:

  number of bootstrap replicates.

- model:

  the initial fitted model

## Value

A list of length `n_models`, with each entry representing one
bootstrapped model as object of class
[svine](https://tnagler.github.io/svines/reference/svine.md).

## References

Nagler, T., Krüger, D., and Min, A. (2022). Stationary vine copula
models for multivariate time series. *Journal of Econometrics*, 227(2),
305–324.
[doi:10.1016/j.jeconom.2021.11.015](https://doi.org/10.1016/j.jeconom.2021.11.015)
.

## Examples

``` r
data(returns)
dat <- returns[1:100, 1:2]

# fit parametric S-vine model with Markov order 1
model <- svine(dat, p = 1, family_set = "parametric")

# compute 10 bootstrap replicates of the model
boot_models <- svine_bootstrap_models(10, model)

# compute bootstrap replicates of 90%-quantile of X_1 + X_2.
mu_boot <- sapply(
  boot_models,
  function(m) {
    xx <- rowSums(t(svine_sim(1, 10^2, m, past = dat)[1, ,]))
    quantile(xx, 0.9)
  }
) 
```
