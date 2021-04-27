#' Sample from the posterior of a 1-d binary population (quickly)
#'
#' @param obs numeric vector representing the observed sample.
#' @param N size of the population.
#' @param prior choice of prior.
#' @param control a list of parameters for tuning the mcmc algorithm. Passed to \code{\link{lemur.control}}.
#' @param ... arguments to be used to form the \code{control} object if it is not supplied directly.
#' @importFrom coda mcmc
#' @importFrom stats runif
#' @importFrom stats rbinom
#' @importFrom stats window
#' @export
mcmc_binary_fast <- function(obs, N, prior = "beta", control = list(...), ...) {
  control <- do.call("lemur.control", control)

  # return(control)

  if(prior == "beta") {
    lprior <- logprior_binary_beta_fast_
    formals(lprior)$a <- control$a
    formals(lprior)$b <- control$b
  } else {
    stop("Unknown prior specified.")
  }

  n_1 <- sum(obs)
  n_0 <- length(obs) - n_1

  n_row <- control$nsteps + control$burnin

  stor <- matrix(0, nrow = n_row, ncol = N)
  stor[1, ] <- rbinom(N, size = 1, prob = mean(obs))

  lp0 <- logpost_binary_fast_(mean(stor[1, ]), n_0, n_1, lprior)

  for(i in 2:n_row) {
    Ycur <- stor[i-1, ]
    for(k in 1:N) {
      can <- Ycur
      can[k] <- ifelse(Ycur[k] == 1, 0, 1)
      lpcan <- logpost_binary_fast_(mean(can), n_0, n_1, lprior)
      lhastTo <- loghastings_binomial_(Ycur, can, k)
      lhastFrom <- loghastings_binomial_(can, Ycur, k)

      u <- runif(1)
      if(log(u) < ((lpcan - lp0) + (lhastFrom - lhastTo))) {
        Ycur <- can
        lp0 <- lpcan
      }
    }
    stor[i, ] <- Ycur
  }

  raw <- coda::mcmc(stor)
  samples <- window(raw, start = control$burnin + 1, thin = control$thin)
  out <- list(samples = samples, control = control)
  return(out)
}


#' Calculate the log-posterior of a 1-d binary population (quickly)
#'
#' @param pop numeric vector representing the population
#' @param obs numeric vector representing the observed sample
#' @param lprior function uses to calculate the log-prior
#' @param a fix later
#' @param b fix later
#' @param ... other arguments, passed to the lprior function
#' @export
logpost_binary_fast <- function(pop, obs, lprior, a, b, ...) {

  P <- mean(pop)
  n_1 <- sum(obs)
  n_0 <- length(obs) - sum(obs)

  logpost_binary_fast_(P, n_0, n_1, lprior, a, b, ...)
}


#' Calculate the log-posterior for a 1-d binary population (quickly)
#'
#' @param P the population proportion
#' @param n_0 the number of failures in the sample
#' @param n_1 the number of successes in the sample
#' @param lprior function uses to calculate the log-prior
#' @param ... other arguments, passed to the lprior function
logpost_binary_fast_ <- function(P, n_0, n_1, lprior, ...) {

  log_lik <- loglik_binary_fast_(P, n_0, n_1)
  log_prior <- lprior(P, ...)

  log_lik + log_prior
}


#' Calculate the log-likelihood for a 1-d binary population (quickly)
#'
#' @param P the population proportion
#' @param n_0 the number of failures in the sample
#' @param n_1 the number of successes in the sample
#' @importFrom stats dbeta
loglik_binary_fast_ <- function(P, n_0, n_1) {
  dbeta(P, n_1 + 1, n_0 + 1, log = TRUE)
}


#' Calculates the beta log-prior for a 1-d binary population (quickly)
#'
#' @param P the population proportion
#' @param a first shape parameter for beta prior
#' @param b second shape parameter for beta prior
#' @param ... other arguments (unused)
#' @importFrom stats dbeta
logprior_binary_beta_fast_ <- function(P, a, b, ...) {
  dbeta(P, a, b, log = TRUE)
}


#' Calculate the log hastings adjustment for metropolis algorithm (quickly)
#'
#' @param Ycur numeric vector for current value of population
#' @param Yprop numeric vector for proposed value of population
#' @param k position of the changing value
#' @importFrom stats dbinom
loghastings_binomial_ <- function(Ycur, Yprop, k) {
  p <- mean(Ycur)
  dbinom(Ycur[k], size = 1, prob = p, log = TRUE)
}
