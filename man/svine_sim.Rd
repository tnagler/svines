% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svine_methods.R
\name{svine_sim}
\alias{svine_sim}
\title{Simulate from a S-vine model}
\usage{
svine_sim(n, rep, model, past = NULL, qrng = FALSE, cores = 1)
}
\arguments{
\item{n}{how many steps of the time series to simulate.}

\item{rep}{number of replications; \code{rep} time series of length \code{n} are
generated.}

\item{model}{a S-vine copula model object (inheriting from \link{svinecop_dist}).}

\item{past}{(optional) matrix of past observations. If provided, time series
are simulated conditional on the past.}

\item{qrng}{if \code{TRUE}, generates quasi-random numbers using the multivariate
Generalized Halton sequence up to dimension 300 and the Generalized Sobol
sequence in higher dimensions (default \code{qrng = FALSE}).}

\item{cores}{number of cores to use; if larger than one, computations are
done parallel over replications.}
}
\value{
An \code{n}-by-\code{d}-by\code{rep} array, where \code{d} is the cross-sectional
dimension of the model. This reduces to an \code{n}-by-\code{d} matrix if \code{rep == 1}.
}
\description{
Simulate from a S-vine model
}
\examples{
# load data set
data(returns)  
returns <- returns[1:100, 1:3]

# fit parametric S-vine model with Markov order 1
fit <- svine(returns, p = 1, family_set = "parametric")

pairs(returns)  # original data
pairs(svine_sim(100, rep = 1, model = fit))   # simulated data

# simulate the next day conditionally on the past 500 times
pairs(t(svine_sim(1, rep = 100, model = fit, past = returns)[1, , ]))
}
