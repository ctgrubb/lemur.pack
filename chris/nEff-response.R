options(scipen = 5)

library(extraDistr)
library(dplyr)
library(tidyr)
library(ggplot2)

# x <- c(0.14, 0.21, 0.23, 0.37, 0.41, 0.43, 0.51, 0.61, 0.84, 0.94)
# x <- c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1)
# mean(x)
# N <- 100
# samples <- 1000
# nEff <- 1
# alpha <- 0.5

graph_prior <- function(x, N, samples = 1000, nEff = length(sample), alpha = 0.5) {
  priorDraws <- t(apply(matrix(rdirichlet(samples, (nEff + alpha) / (N + 1) * rep(1, N + 1)), nrow = samples)[, -(N + 1), drop = FALSE], 1, cumsum))

  hist(apply(priorDraws, 1, function(x) sum(x < 0.5)))
}

# graph_prior(x, N = 100, samples = 1000, nEff = 5, alpha = 0.1)

fun1 <- function(x, N, samples = 1000, nEff = length(sample), alpha = 0.5) {

  n <- length(x)
  nLow <- sum(x < 0.5)

  priorDraws <- t(apply(matrix(rdirichlet(samples, (nEff + alpha) / (N + 1) * rep(1, N + 1)), nrow = samples)[, -(N + 1), drop = FALSE], 1, cumsum))
  NLow <- apply(priorDraws, 1, function(x) sum(x < 0.5))
  weights <- dhyper(nLow, NLow, N - NLow, n)

  draws <- sample(1:samples, replace = TRUE, prob = weights)

  return(priorDraws[draws, ])
}

fun2 <- function(x, N, samples = 1000, nEff = length(sample), alpha = 0.5) {

  lpost <- function(nLow, pop, N, nEff, alpha) {
    NLow <- sum(pop < 0.5)
    diffs <- diff(c(0, sort(pop), 1))
    lprior <- ddirichlet(diffs, (nEff + alpha) / (N + 1) * rep(1, N + 1), log = TRUE)
    llik <- dhyper(nLow, NLow, N - NLow, n, log = TRUE)
    return(llik + lprior)
  }

  n <- length(x)
  nLow <- sum(x < 0.5)

  out <- matrix(NA, nrow = samples, ncol = N)
  out[1, ] <- runif(N)

  lp0 <- lpost(nLow, out[1, ], N, nEff, alpha)

  for(i in 2:samples) {
    current <- out[i - 1, ]
    for(k in 1:N) {
      proposal <- current
      proposal[k] <- runif(1)
      lp1 <- lpost(nLow, proposal, N, nEff, alpha)

      u <- runif(1)
      if(is.finite(lp1) && log(u) < (lp1 - lp0)) {
        current <- proposal
        lp0 <- lp1
      }
    }
    out[i, ] <- current
  }
  return(out)
}


plot(density(rowMeans(fun1(x, N = 100, nEff = 10, alpha = 0.5, samples = 1000))), col = "red",
     main = "Density Comparison of Importance Sampling vs Markov Chain Approach", xlab = "Mean of Population")
lines(density(rowMeans(fun2(x, N = 100, nEff = 10, alpha = 0.5, samples = 1000))), col = "blue")

par(mfrow = c(1, 2))
hist(rowMeans(fun1(x, N = 100, nEff = 10, alpha = 0.5, samples = 10000)))
hist(rowMeans(fun2(x, N = 100, nEff = 10, alpha = 0.5, samples = 10000)))
par(mfrow = c(1, 1))


fun1(x, 100, 1, 1, 0.5)

fun2_test <- fun2(x, N = 40, nEff = 1, samples = 10000)
hist(rowMeans(fun2_test))
mean(fun2_test)

plot(1:10000, rowMeans(fun2_test))

qqplot(rowMeans(fun1_test), rowMeans(fun2_test), xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1)

test1 <- fun1(c(rep(0, 20), rep(1, 0)), N = 20, samples = 10000, nEff = 4, alpha = 0.5)
test2 <- fun1(c(rep(0, 8), rep(1, 2)), N = 20, samples = 100000, nEff = 4, alpha = 0.5)

test3 <- fun2(c(rep(0, 20), rep(1, 0)), N = 20, samples = 1000, nEff = 4, alpha = 0.5)
test4 <- fun2(c(rep(0, 8), rep(1, 2)), N = 20, samples = 100000, nEff = 4, alpha = 0.5)

par(mfrow = c(1, 2))
plot(rowMeans(test1))
plot(rowMeans(test4))

plot(apply(test3, 1, quantile, prob = 0.1))
plot(apply(test3, 1, quantile, prob = 0.9))

par(mfrow = c(2, 2))
hist(rowMeans(test1), xlim = c(0, 1))
hist(rowMeans(test2), xlim = c(0, 1))
hist(rowMeans(test3), xlim = c(0, 1))
hist(rowMeans(test4), xlim = c(0, 1))

mean(test1); mean(test2)
mean(test3); mean(test4)


# Comparison (Importance Sampling vs MCMC) ------------------------------------------------------------------------

samp <- 100000

x1 <- c(0.11, 0.19, 0.27, 0.39, 0.42, 0.48, 0.59, 0.71, 0.84, 0.92)
x2 <- c(0.07, 0.12, 0.19, 0.26, 0.31, 0.36, 0.41, 0.49, 0.67, 0.89)

fun1_x1 <- fun1(x = x1, N = 40, samples = samp, nEff = 4, alpha = 0.5)
fun2_x1 <- fun2(x = x1, N = 40, samples = samp, nEff = 4, alpha = 0.5)
fun1_x2 <- fun1(x = x2, N = 40, samples = samp, nEff = 4, alpha = 0.5)
fun2_x2 <- fun2(x = x2, N = 40, samples = samp, nEff = 4, alpha = 0.5)

df <- bind_rows(data.frame(fun1_x1), data.frame(fun2_x1), data.frame(fun1_x2), data.frame(fun2_x2))

means <- apply(df, 1, mean)
medians <- apply(df, 1, median)
P10 <- apply(df, 1, quantile, prob = 0.1)
P25 <- apply(df, 1, quantile, prob = 0.25)
P75 <- apply(df, 1, quantile, prob = 0.75)
P90 <- apply(df, 1, quantile, prob = 0.9)

df2 <- data.frame(
  X = c(rep("X1", 2 * samp), rep("X2", 2 * samp)),
  Method = c(rep("Acceptance", samp), rep("MCMC", samp), rep("Acceptance", samp), rep("MCMC", samp)),
  Mean = means,
  Median = medians,
  P10 = P10,
  P25 = P25,
  P75 = P75,
  P90 = P90
)

ggplot(data = df2, mapping = aes(x = P75)) +
  facet_wrap(Method ~ X) +
  geom_histogram()



