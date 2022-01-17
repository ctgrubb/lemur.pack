library(latex2exp)

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
N <- 10
n <- 2
pop <- sort(rnorm(N))
samp <- sample(pop, n, replace = TRUE)

ord.samp <- sort(samp)

# P(X_(1) = x_(1))
nb <- sum(pop > ord.samp[1])
step_1 <- (1 - (nb / N) ** 2) - (1 - ((nb + 1) / N) ** 2)

# P(X_(2) == x_(2) | X_(1) = x_(1))
step_2 <- 1 / (nb + 1)


# Matrix heatmap
storage <- matrix(0, nrow = N, ncol = N)

for(i in 1:N) {
  for(j in i:N) {
    ord.samp <- sort(pop[c(i, j)])

    nb <- sum(pop > ord.samp[1])
    step_1 <- (1 - (nb / N) ** 2) - (1 - ((nb + 1) / N) ** 2)
    step_2 <- 1 / (nb + 1)

    storage[i, j] <- step_1 * step_2
  }
}

sum(storage)

which.max(storage)

image(x = 1:10, y = 1:10, z = t(apply(storage, 2, rev)))
image(x = 1:10, y = 1:10, z = storage, xlab = TeX("$X_{(1)}$"), ylab = TeX("$X_{(2)}$"))
heatmap(storage, Colv = NA, Rowv = NA, symm = TRUE)



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

  step_2 <- 1 - choose(sum(pop > samp[2]), 9) * choose(sum(pop <= samp[1]), 1) / choose(100, 10) / step_1

  probs <- sapply(1:n, function(x) choose(gt.tallies[x], n - (x - 1)))
  for(i in 1:(n-1)) {
    probs <- probs * c(rep(1, i), rep(choose(lte.tallies[i], 1), n - i))
  }
  probs <- probs / choose(N, n)
  probs[1] <- 1 - probs[1]
  for(i in 2:n) {
    probs[i] <- 1 - probs[i] / probs[i - 1]
  }


  test <- rep(0, 10)
  test[1] <- 1 - choose(gt.tallies[1], n) / choose(N, n)
  test[2] <- 1 - choose(gt.tallies[2], n - 1) * choose(lte.tallies[1], 1) / choose(N, n) / test[1]
  test[3] <- 1 - choose(gt.tallies[3], n - 2) * choose(lte.tallies[2], 1) * choose(lte.tallies[1], 1) / choose(N, n) /
    test[2]
  test[4] <- 1 - choose(gt.tallies[4], n - 3) * choose(lte.tallies[3], 1) * choose(lte.tallies[2], 1) *
    choose(lte.tallies[1], 1) / choose(N, n) / test[3]
  test[5] <- 1 - choose(gt.tallies[5], n - 4) * choose(lte.tallies[4], 1) * choose(lte.tallies[3], 1) *
    choose(lte.tallies[2], 1) * choose(lte.tallies[1], 1) / choose(N, n) / test[4]
  test[6] <- 1 - choose(gt.tallies[6], n - 5) * choose(lte.tallies[5], 1) * choose(lte.tallies[4], 1) *
    choose(lte.tallies[3], 1) * choose(lte.tallies[2], 1) * choose(lte.tallies[1], 1) / choose(N, n) / test[5]
  test[7] <- 1 - choose(gt.tallies[7], n - 6) * choose(lte.tallies[6], 1) * choose(lte.tallies[5], 1) *
    choose(lte.tallies[4], 1) * choose(lte.tallies[3], 1) * choose(lte.tallies[2], 1) * choose(lte.tallies[1], 1) /
    choose(N, n) / test[6]
  test[8] <- 1 - choose(gt.tallies[8], n - 7) * choose(lte.tallies[7], 1) * choose(lte.tallies[6], 1) *
    choose(lte.tallies[5], 1) * choose(lte.tallies[4], 1) * choose(lte.tallies[3], 1) * choose(lte.tallies[2], 1) *
    choose(lte.tallies[1], 1) / choose(N, n) / test[7]
  test[9] <- 1 - choose(gt.tallies[9], n - 8) * choose(lte.tallies[8], 1) * choose(lte.tallies[7], 1) *
    choose(lte.tallies[6], 1) * choose(lte.tallies[5], 1) * choose(lte.tallies[4], 1) *
    choose(lte.tallies[3], 1) * choose(lte.tallies[2], 1) * choose(lte.tallies[1], 1) / choose(N, n) / test[8]
  test[10] <- 1 - choose(gt.tallies[10], n - 9) * choose(lte.tallies[9], 1) * choose(lte.tallies[8], 1) *
    choose(lte.tallies[7], 1) * choose(lte.tallies[6], 1) * choose(lte.tallies[5], 1) * choose(lte.tallies[4], 1) *
    choose(lte.tallies[3], 1) * choose(lte.tallies[2], 1) * choose(lte.tallies[1], 1) / choose(N, n) / test[10]
}


iters <- 100000
pop <- 1:10
n <- 3
storage <- matrix(0, nrow = iters, ncol = n)

for(i in 1:iters) {
  storage[i, ] <- samp <- sample(pop, n, replace = TRUE)
}

storage2 <- t(apply(storage, 1, sort))
storage3 <- storage2[, 1] == 3
storage4 <- storage[storage3, 2:3]

table(storage4[, 1], storage4[, 2])

table(storage[storage[, 1] == 3, 2], storage[storage[, 1] == 3, 3])
