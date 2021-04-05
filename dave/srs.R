# Functional implementation of Dave's srs2 code
# Minor modification to make first sample 50/50 0/1 instead of all 0
srs2 <- function(obs, N, nsteps) {

  logpost = function(Y, y) {
    n <- length(obs)
    X <- sum(Y)
    x <- sum(obs)

    if(x > X) return(-9e9)
    if(n - x > N - X) return(-9e9)

    llike <- dhyper(x, X, N - X, n, log = T)
    lprior <- -lchoose(N, X)

    llike + lprior
  }

  n <- length(obs)
  Ysamp <- matrix(c(0, 1), ncol = N, nrow = nsteps, byrow = TRUE)
  lp0 <- logpost(Ysamp[1, ], obs)

  for(i in 2:nsteps) {
    Ycur <- Ysamp[i - 1, ]
    for(k in 1:N) {
      can <- Ycur
      can[k] <- ifelse(Ycur[k] == 1, 0, 1)
      lpcan <- logpost(can, obs)
      u <- runif(1)
      if(log(u) < (lpcan - lp0)) {
        Ycur <- can
        lp0 <- lpcan
      } else {
        Ycur <- Ycur
      }
    }
    Ysamp[i, ] <- Ycur
  }
  return(Ysamp)
}


# Functional implementation of Dave's srs3 code
# Minor modification to make first sample 50/50 0/1 instead of all 0
# Minor modification to check that the population contains at least as many 0 and 1 as the sample
srs3 <- function(obs, N, nsteps) {

  logpost <- function(Y, n0 = sum(obs == 0), n1 = sum(obs == 1)) {
    X <- sum(Y)

    if(n1 > X) return(-9e9)
    if(n0 > N - X) return(-9e9)

    p <- X / N
    logDensP <- dbeta(p, n1 + 1, n0 + 1, log = TRUE)
    logPrior <- dbeta(p, 1, 1, log = TRUE)
    logDensP + logPrior
  }

  qhast <- function(Y, Ycan, i) {
    p <- mean(Y)
    dbinom(Y[i], size = 1, prob = p, log = TRUE)
  }

  Ysamp <- matrix(c(0, 1), ncol = N, nrow = nsteps, byrow = TRUE)
  lp0 <- logpost(Ysamp[1, ])

  for(i in 2:nsteps) {
    Ycur <- Ysamp[i - 1, ]
    for(k in 1:N) {
      can <- Ycur
      can[k] <- ifelse(Ycur[k] == 1, 0, 1)
      lpcan <- logpost(can)
      lq01 <- qhast(Ycur, can, k)
      lq10 <- qhast(can, Ycur, k)
      u <- runif(1)
      if(log(u) < (lpcan - lp0 + lq10 - lq01)) {
        Ycur <- can
        lp0 <- lpcan
      }
    }
    Ysamp[i, ] <- Ycur
  }
  return(Ysamp)
}
