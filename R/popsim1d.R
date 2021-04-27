#' Sample from the posterior of a 1-d binary population
#'
#' @param obs numeric vector representing the observed sample.
#' @param N size of the population.
#' @param prior choice of prior.
#' @param control a list of parameters for tuning the mcmc algorithm. Passed to \code{\link{lemur.control}}.
#' @param ... arguments to be used to form the \code{control} object if it is not supplied directly.
#' @importFrom coda mcmc
#' @export
mcmc_binary <- function(obs, N, prior = "uninformative", control = list(...), ...) {
  control <- do.call("lemur.control", control)

  if(prior == "uninformative") {
    lprior <- logprior_binary_uninformative_
  } else if(prior == "flat") {
    lprior <- logprior_binary_flat_
  }

  n <- length(obs)
  y <- sum(obs)
  n_row <- control$nsteps + control$burnin

  stor <- matrix(0, nrow = n_row, ncol = N)

  lp0 <- logpost_binary_(N, n, 0, y)

  for(i in 2:n_row) {
    Ycur <- stor[i-1, ]
    for(k in 1:N) {
      can <- Ycur
      can[k] <- ifelse(Ycur[k] == 1, 0, 1)
      lpcan <- logpost_binary_(N, n, sum(can), y, lprior)

      u <- runif(1)
      if(log(u) < (lpcan - lp0)){
        Ycur <- can
        lp0 <- lpcan
      }
    }
    stor[i, ] <- Ycur
  }

  raw <- coda::mcmc(stor)
  out <- window(raw, start = control$burnin + 1, thin = control$thin)
  return(out)
}


#' Sample from the posterior of a 1-d multi-class population
#'
#' @param obs numeric vector representing the observed sample.
#' @param N size of the population.
#' @param prior choice of prior.
#' @param control a list of parameters for tuning the mcmc algorithm. Passed to \code{\link{lemur.control}}.
#' @param ... arguments to be used to form the \code{control} object if it is not supplied directly.
#' @importFrom coda mcmc
#' @export
mcmc_multiclass <- function(obs, N, prior = "uninformative", control = list(...), ...) {
  control <- do.call("lemur.control", control)

  if(prior == "uninformative") {
    lprior <- logprior_multiclass_uninformative_
  } else if(prior == "flat") {
    lprior <- logprior_multiclass_flat_
  }

  n <- length(obs)
  y <- as.numeric(table(obs))
  n_class <- length(y)
  n_row <- control$nsteps + control$burnin

  stor <- matrix(1, nrow = n_row, ncol = N)

  lp0 <- logpost_multiclass_(N, n, as.numeric(table(factor(stor[1, ], levels = 1:n_class))), y)

  for(i in 2:n_row) {
    Ycur <- stor[i-1, ]
    for(k in 1:N) {
      can <- Ycur
      can[k] <- sample((1:n_class)[-can[k]], 1)
      lpcan <- logpost_multiclass_(N, n, as.numeric(table(factor(can, levels = 1:n_class))), y, lprior)

      u <- runif(1)
      if(log(u) < (lpcan - lp0)){
        Ycur <- can
        lp0 <- lpcan
      }
    }
    stor[i, ] <- Ycur
  }

  raw <- coda::mcmc(stor)
  out <- window(raw, start = control$burnin + 1, thin = control$thin)
  return(out)
}

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


#' Calculate the log-posterior of a 1-d multi-class population
#'
#' @param pop numeric vector representing the population
#' @param obs numeric vector representing the observed sample
#' @param lprior function uses to calculate the log-prior
#' @export
logpost_multiclass <- function(pop, obs, lprior) {

  N <- length(pop)
  n <- length(obs)

  Y <- as.numeric(table(pop))
  y <- as.numeric(table(obs))

  logpost_multiclass_(N, n, Y, y, lprior)
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

  log_lik + log_prior
}


#' Calculate the log-posterior for a 1-d multi-class population
#'
#' @param N the population size
#' @param n the sample size
#' @param Y a vector indicating the amount of each class in the population.
#' @param y a vector indicating the amount of each class in the sample.
#' @param lprior function uses to calculate the log-prior
logpost_multiclass_ <- function(N, n, Y, y, lprior) {

  if(any(Y < y)) {
    return(-9e9)
  }

  log_lik <- loglik_multiclass_(n, Y, y)
  log_prior <- lprior(N, Y)

  log_lik + log_prior
}


#' Calculate the log-likelihood for a 1-d binary population
#'
#' @param N the population size
#' @param n the sample size
#' @param Y the number of successes in the population
#' @param y the number of successes in the sample
#' @importFrom stats dhyper
loglik_binary_ <- function(N, n, Y, y) {
  dhyper(y, Y, N-Y, n, log = TRUE)
}


#' Calculate the log-likelihood for a 1-d multi-class population
#'
#' @param n the sample size.
#' @param Y a vector indicating the amount of each class in the population.
#' @param y a vector indicating the amount of each class in the sample.
#' @importFrom extraDistr dmvhyper
loglik_multiclass_ <- function(n, Y, y) {
  dmvhyper(y, Y, n, log = TRUE)
}


#' Calculates the uninformative log-prior for a 1-d binary population
#'
#' @param N the population size
#' @param Y the number of successes in the population
#' @param ... other arguments (unused)
logprior_binary_uninformative_ <- function(N, Y, ...) {
  -lchoose(N, Y)
}


#' Calculates the uninformative log-prior for a 1-d multi-class population
#'
#' @param N the population size
#' @param Y a vector indicating the amount of each class in the population.
logprior_multiclass_uninformative_ <- function(N, Y) {
  -lmnchoose(N, Y)
}


#' Calculates the flat log-prior for a 1-d binary population
#'
#' @param ... other arguments (unused)
logprior_binary_flat_ <- function(...) {
  0
}


#' Calculates the flat log-prior for a 1-d multi-class population
#'
#' @param ... other arguments (unused)
logprior_multiclass_flat_ <- function(...) {
  0
}
