#' Calculate the log-posterior of a 1-d binary population
#'
#' @param pop numeric vector representing the population
#' @param obs numeric vector representing the observed sample
#' @param lprior function uses to calculate the log-prior
#' @export
logpost_binary <- function(pop, obs, lprior) {

  N <- length(pop)
  n <- length(obs)

  Y <- sum(pop)
  y <- sum(obs)

  logpost_binary_(N, n, Y, y, lprior)
}

#' Calculate the log-posterior for a 1-d binary population
#'
#' @param N the population size
#' @param n the sample size
#' @param Y the number of successes in the population
#' @param y the number of successes in the sample
#' @param lprior function uses to calculate the log-prior
logpost_binary_ <- function(N, n, Y, y, lprior) {

  if(y > Y | (n - y) > (N - Y)) {
    return(-9e9)
  }

  log_lik <- loglik_binary_(N, n, Y, y)
  log_prior <- lprior(N, Y)

  return(log_lik + log_prior)
}

#' Calculate the log-likelihood for a 1-d binary population
#'
#' @param N the population size
#' @param n the sample size
#' @param Y the number of successes in the population
#' @param y the number of successes in the sample
#' @importFrom stats dhyper
loglik_binary_ <- function(N, n, Y, y) {
  return(dhyper(y, Y, N-Y, n, log = TRUE))
}


#' Calculates the uninformative log-prior for a 1-d binary population
#'
#' @param N the population size
#' @param Y the number of successes in the population
#' @param ... other arguments (unused)
logprior_binary_uninformative_ <- function(N, Y, ...) {
  return(-lchoose(N, Y))
}

#' Calculates the flat log-prior for a 1-d binary population
#'
#' @param ... other arguments (unused)
logprior_binary_flat_ <- function(...) {
  return(0)
}
