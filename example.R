library(svines)
library(rugarch)
library(fGarch)

# load data set
data(returns)  
dat <- returns[1:500, 1:2]

# fit parametric S-vine model with Markov order 1
fit <- svine(dat, p = 1, family_set = "parametric")

V <- svine_avar(dat, fit)

par <- svine_get_pars(fit)

par_new <- MASS::mvrnorm(1, par, V)
summary(svine_set_pars(fit, par_new))
