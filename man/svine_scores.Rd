% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standard_errors.R
\name{svine_scores}
\alias{svine_scores}
\title{Score function of parametric S-vine models}
\usage{
svine_scores(x, model, cores = 1)
}
\arguments{
\item{x}{the data.}

\item{model}{S-vine model (inheriting from \link{svine_dist}).}

\item{cores}{number of cores to use.}
}
\value{
A returns a \code{n}-by-\code{k} matrix, where \code{n = NROW(x)} and \code{k} is the
total number of parameters in the
model. Parameters are ordered as follows:
marginal parameters, copula parameters of first tree, copula parameters of
second tree, etc. Duplicated parameters in the copula model are omitted.
}
\description{
Score function of parametric S-vine models
}
\examples{
data(returns)
dat <- returns[1:100, 1:2]

# fit parametric S-vine model with Markov order 1
model <- svine(dat, p = 1, family_set = "parametric")

# Implementation of asymptotic variances
I <- cov(svine_scores(dat, model))
H <- svine_hessian(dat, model)
Hi <- solve(H)
Hi \%*\% I \%*\% t(Hi) / nrow(dat)
}
