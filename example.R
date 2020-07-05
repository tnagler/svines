library(svines)
library(rugarch)
library(fGarch)

# load data set
data(returns)  
dat <- returns[1:500, 1:2]

# fit parametric S-vine model with Markov order 1
fit <- svine(dat, p = 1, family_set = "parametric")

svine_avar(dat, fit)
