library(svines)
library(rugarch)
library(fGarch)

# load data set
data(returns)  
dat <- returns[1:500, 1:3]

# fit parametric S-vine model with Markov order 1
model <- svine(dat, p = 1, margin_families = "norm", family_set = "gauss")


# Implementation of asymptotic variances
I <- cov(svine_scores(x, model, cores))
H <- svine_hessian(x, model, cores)
Hi <- solve(H)
Hi %*% I %*% t(Hi) / nrow(x)


V <- svine_avar(dat, fit, n_lags = 30)

data(returns)  
dat <- returns[1:100, 1:2]

# fit parametric S-vine model with Markov order 0
fit <- svine(dat, p = 0, family_set = "parametric")

new <- svine_sim_se_models(2, fit)
summary(new[[1]])
summary(new[[2]])
