% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popsim_superfast.R
\name{multivariate_multinomial_postpred}
\alias{multivariate_multinomial_postpred}
\title{Sample from the posterior predictive of a multivariate multinomial distribution}
\usage{
multivariate_multinomial_postpred(obs, N, nsamples, priors)
}
\arguments{
\item{obs}{a \code{data.frame} representing the observed sample.}

\item{N}{the size of the predictive samples to create.}

\item{nsamples}{the number of predictive samples to create.}

\item{priors}{either a list of named numeric vectors containing the Dirichlet parameters for individual priors on
each variable or an array describing the joint Dirichlet prior. Can also be a single number, which is interpreted
as a flat joint prior with the given height. See \code{vignette("mmpp", package = "lemur.pack")} for an example of
setting priors.}
}
\description{
\loadmathjax
\code{multivariate_multinomial_postpred} samples from the posterior predictive distribution given an observation from
a multivariate multinomial distribution and a prior on the probabilities of each category.
}
\details{
I will complete this later!
}
