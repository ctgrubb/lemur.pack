iters <- 50000

storage <- matrix(0, nrow = iters, ncol = 2)

for(i in 1:iters) {
  pop.others <- rnorm(90)
  # samp <- rnorm(10)
  samp <- qnorm(1:10/11)
  pop <- c(pop.others, samp)

  ord.samp <- sort(samp)

  # P(X_(1) <= x_(1))
  storage[i, 1] <- 1 - choose(sum(pop > ord.samp[1]), 10) / choose(100, 10)

  # P(X_(2) <= x_(2) | X_(1) <= x_(1))
  storage[i, 2] <- 1 - choose(sum(pop > ord.samp[2]), 9) * choose(sum(pop <= ord.samp[1]), 1) / choose(100, 10) / storage[i, 1]

}

hist(storage[, 1])
hist(storage[, 2])


# Testing
pop.others <- rnorm(90)
samp <- rnorm(10)
pop <- c(pop.others, samp)

ord.samp <- sort(samp)

# P(X_(1) <= x_(1))
step_1 <- 1 - choose(sum(pop > ord.samp[1]), 10) / choose(100, 10)

# P(X_(2) <= x_(2) | X_(1) <= x_(1))
step_2 <- 1 - choose(sum(pop > ord.samp[2]), 9) * choose(sum(pop <= ord.samp[1]), 1) / choose(100, 10) / step_1



# "Ideal" sample
pop.others <- rnorm(90)
samp <- qnorm(1:10/11)
pop <- c(pop.others, samp)

ord.samp <- sort(samp)

# P(X_(1) <= x_(1))
step_1 <- 1 - choose(sum(pop > ord.samp[1]), 10) / choose(100, 10)

# P(X_(2) <= x_(2) | X_(1) <= x_(1))
step_2 <- 1 - choose(sum(pop > ord.samp[2]), 9) * choose(sum(pop <= ord.samp[1]), 1) / choose(100, 10) / step_1


pop <- rnorm(100)
samp <- sample(pop, 10)

dave.func <- function(pop, samp) {
  samp <- sort(samp)

  N <- length(pop)
  n <- length(samp)

  gt.tallies <- sapply(samp, function(x) sum(pop > x))
  lte.tallies <- N - gt.tallies

  probs <- sapply(1:n, function(x) choose(gt.tallies[x], n - (x - 1)))
  for(i in 1:(n-1)) {
    probs <- probs * c(rep(1, i), rep(choose(lte.tallies[i], 1), n - i))
  }
  probs <- probs / choose(N, n)
  probs <- probs / c(1, probs[-n])
  probs <- 1 - probs

  1 - choose(gt.tallies[1], n) / choose(N, n)
  1 - choose(gt.tallies[2], n - 1) * choose(lte.tallies[1], 1) / choose(N, n)
  1 - choose(gt.tallies[3], n - 2) * choose(lte.tallies[2], 1) * choose(lte.tallies[1], 1) / choose(N, n)
  # ...
  1 - choose(gt.tallies[10], n - 9)
}

