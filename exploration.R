# obs <- c(1, 2, 1, 2, 3, 3, 1, 1, 1)
obs <- c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1)
table(obs) / length(obs)
N <- 100


test <- mcmc_multiclass(obs, N)

table(as.matrix(test))

hist(apply(as.matrix(test), 1, function(x) table(factor(x, levels = 1:3)))[1, ])
hist(apply(as.matrix(test), 1, function(x) table(factor(x, levels = 1:3)))[2, ])
hist(apply(as.matrix(test), 1, function(x) table(factor(x, levels = 1:3)))[3, ])



# 3/19 comparisons


obs <- rbinom(10, 1, 0.4)
mean(obs)
N <- 100

test <- mcmc_binary_fast(obs, N, nsteps = 10000)
test2 <- mcmc_binary_fast(obs, N, a = 0.5, b = 0.5, nsteps = 10000)
test3 <- mcmc_binary_superfast(obs, N, 10000, 0.5, 0.5)
test4 <- mcmc_binary_superfast(obs, N, 1000000, 1, 1)
test5 <- mcmc_binary(obs, N, nsteps = 10000)

# a = 1, b = 1
hist(apply(as.matrix(test5), 1, sum)) # slower documented function
d <- density(apply(test4, 1, sum))
d2 <- density(apply(Ysamp, 2, sum))
plot(d)
lines(d2)
# polygon(d, col = "red")
# lines(d$x, dbeta(d$x / 100, sum(obs) + 1, length(obs) - sum(obs) + 1) / (max(d$x) - min(d$x)))
lines(d$x, dbeta(d$x / 100, sum(obs) + 1, length(obs) - sum(obs) + 1) / 100)

hist(apply(as.matrix(test$samples), 1, sum)) # my documented function
hist(apply(Ysamp, 2, sum)) # dave's code
hist(apply(test4, 1, sum)) # my posterior -> predictive function

# a = 0.5, b = 0.5
hist(apply(as.matrix(test2$samples), 1, sum)) # my documented function
hist(apply(test3, 1, sum)) # my posterior -> predictive function


# 3/23 meeting
obs <- c(0, 1, 1, 0, 2, 3, 1, 1, 0, 2, 3, 1, 0, 2, 1)
table(obs)/length(obs)

Ysamp <- mcmc_multinomial_superfast(obs, 100, 1000, 0.5)



obs <- data.frame(
  A = c(sample(1:4, 16, replace = TRUE), 1:4),
  B = c(sample(1:3, 17, replace = TRUE, prob = c(0.5, 0.2, 0.3)), 1:3),
  C = c(sample(1:6, 14, replace = TRUE, prob = c(0.1, 0.2, 0.3, 0.2, 0.1, 0.1)), 1:6)
)

Ysamp <- multivariate_multinomial_postpred(obs, 1000, 100, 0.5)



# Synthetic Data
grouped2d_pop <- data.frame(
  Group = sample(LETTERS[1:3], size = 5000, replace = TRUE, prob = c(0.2, 0.45, 0.35)),
  X1 = NA,
  X2 = NA
)

grouped2d_pop$X1[grouped2d_pop$Group == "A"] <- rt(sum(grouped2d_pop$Group == "A"), 4, ncp = 350)
grouped2d_pop$X1[grouped2d_pop$Group == "B"] <- rt(sum(grouped2d_pop$Group == "B"), 6, ncp = 250)
grouped2d_pop$X1[grouped2d_pop$Group == "C"] <- rt(sum(grouped2d_pop$Group == "C"), 8, ncp = 150)

grouped2d_pop$X2[grouped2d_pop$Group == "A"] <- 0.6 * grouped2d_pop$X1[grouped2d_pop$Group == "A"] + rnorm(sum(grouped2d_pop$Group == "A"), 0, 100)
grouped2d_pop$X2[grouped2d_pop$Group == "B"] <- 0.5 * grouped2d_pop$X1[grouped2d_pop$Group == "B"] + rnorm(sum(grouped2d_pop$Group == "B"), 0, 70)
grouped2d_pop$X2[grouped2d_pop$Group == "C"] <- 0.4 * grouped2d_pop$X1[grouped2d_pop$Group == "C"] + rnorm(sum(grouped2d_pop$Group == "C"), 0, 50)

hist(grouped2d_pop$X1, probability = TRUE, breaks = 20)
hist(grouped2d_pop$X2, probability = TRUE, breaks = 20)

plot(grouped2d_pop$X1, grouped2d_pop$X2, col = as.numeric(as.factor(grouped2d_pop$Group)))





