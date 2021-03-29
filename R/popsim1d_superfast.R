#' Sample from the posterior predictive of a 1-d binary population (very quickly)
#'
#' @param obs numeric vector representing the observed sample.
#' @param N size of the population.
#' @param nsamples number of samples to draw
#' @param a first hyperparameter for beta prior
#' @param b second hyperparameter for beta prior
#' @importFrom stats rbeta
#' @export
mcmc_binary_superfast <- function(obs, N, nsamples, a, b) {

  n_1 <- sum(obs)
  n_0 <- length(obs) - n_1

  P <- rbeta(nsamples, n_1 + a, n_0 + b)
  Y <- t(sapply(P, function(x) sample(c(0, 1), size = N, replace = TRUE, prob = c(1 - x, x))))

  return(Y)
}


#' Sample from the posterior predictive of a 1-d multinomial population (very quickly)
#'
#' @param obs numeric vector representing the observed sample.
#' @param N size of the population.
#' @param nsamples number of samples to draw
#' @param alpha hyperparameter for dirichlet prior; can be either a single number or a numeric vector
#' @importFrom extraDistr rdirichlet
#' @export
mcmc_multinomial_superfast <- function(obs, N, nsamples, alpha) {

  x <- unique(obs)
  n <- as.numeric(table(obs))

  if(length(alpha) == 1) {
    alpha <- rep(alpha, length(n))
  } else if(length(alpha) != length(n)) {
    stop("Length of alpha hyperparameter must be 1 or same as number of categories in observation.")
  }

  P <- rdirichlet(nsamples, n + alpha)
  Y <- t(apply(P, 1, function(p) sample(x, size = N, replace = TRUE, prob = p)))

  return(Y)
}
