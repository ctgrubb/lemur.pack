# Functional implementation of Dave's srs2 code
# Minor modification to make first sample 50/50 0/1 instead of all 0
srs2 <- function(obs, N, nsteps) {

  logpost = function(Y, y) {
    N <- length(Y)
    n <- length(obs)
    X <- sum(Y)
    x <- sum(obs)

    if(x > X) return(-9e9)
    if(n - x > N - X) return(-9e9)

    llike <- dhyper(x, X, N - X, n, log = T)
    lprior <- -lchoose(N, X)

    return(llike + lprior)
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
srs3 <- function(obs, N, nsteps) {



  return(Ysamp)
}
