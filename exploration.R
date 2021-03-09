obs <- c(1, 2, 1, 2, 3, 3, 1, 1, 1)
obs <- c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1)
table(obs) / length(obs)
N <- 100


test <- mcmc_multiclass(obs, N)

table(as.matrix(test))

hist(apply(as.matrix(test), 1, function(x) table(factor(x, levels = 1:3)))[1, ])
hist(apply(as.matrix(test), 1, function(x) table(factor(x, levels = 1:3)))[2, ])
hist(apply(as.matrix(test), 1, function(x) table(factor(x, levels = 1:3)))[3, ])
