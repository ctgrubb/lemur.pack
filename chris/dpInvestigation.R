library(extraDistr)

iters <- 10000
n <- 1000

storage <- list(
  dpMethod = matrix(NA, nrow = iters, ncol = n),
  daveMethod = matrix(NA, nrow = iters, ncol = n)
)

for(i in 1:iters) {
  # Dirichlet process
  c <- 10000 # alpha from DP
  G_0 <- function(n) rnorm(n, 0, 1) # base function
  b <- rbeta(n, 1, c)
  p <- numeric(n)
  p[1] <- b[1]
  p[2:n] <- sapply(2:n, function(i) b[i] * prod(1 - b[1:(i-1)]))
  y <- G_0(n)
  storage$dpMethod[i, ] <- sample(y, prob = p, replace = TRUE)

  # Dave's idea
  nEff <- 10000 # n effective
  alpha <- 0.5
  priorDrawsRaw <- as.numeric(rdirichlet(1, (nEff + alpha) / (n + 1) * rep(1, n + 1)))
  priorDraws <- cumsum(priorDrawsRaw[1:n])
  storage$daveMethod[i, ] <- qnorm(priorDraws, 0, 1) # base function
}

par(mfrow = c(2, 2))

plot(density(apply(storage$dpMethod, 1, mean)), col = "red", main = "Mean")
lines(density(apply(storage$daveMethod, 1, mean)), col = "blue")

plot(density(apply(storage$dpMethod, 1, sd)), col = "red", main = "StDev")
lines(density(apply(storage$daveMethod, 1, sd)), col = "blue")

plot(density(apply(storage$dpMethod, 1, min)), col = "red", main = "Min")
lines(density(apply(storage$daveMethod, 1, min)), col = "blue")

plot(density(apply(storage$dpMethod, 1, max)), col = "red", main = "Max")
lines(density(apply(storage$daveMethod, 1, max)), col = "blue")

apply(storage$dpMethod, 1, function(x) length(unique(x)))
apply(storage$daveMethod, 1, function(x) length(unique(x)))




