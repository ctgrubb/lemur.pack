# Functional implementation of posterior predictive sampling
# The posterior (assuming you use a beta prior) is a beta distribution, and thus
# the posterior predictive is beta-binomial
beta_binomial <- function(obs, N, nsteps, a = 1, b = 1) {
  n_1 <- sum(obs)
  n_0 <- length(obs) - n_1

  P <- rbeta(nsteps, n_1 + a, n_0 + b)
  Y <- t(sapply(P, function(x) sample(c(0, 1), size = N, replace = TRUE, prob = c(1 - x, x))))
}
