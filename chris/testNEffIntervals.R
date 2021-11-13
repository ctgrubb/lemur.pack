fun1 <- function(samples = 10000, N = 1000, nEff = 1000, alpha = 0.5) {

  priorDraws <- t(apply(matrix(rdirichlet(samples, (nEff + alpha) / (N + 1) * rep(1, N + 1)), nrow = samples)[, -(N + 1), drop = FALSE], 1, cumsum))

  tns <- matrix(qtnorm(priorDraws, mean = 100000, sd = 30000, 0, Inf), nrow = samples)



  return(out)
}




