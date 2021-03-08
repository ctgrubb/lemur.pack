obs <- c(1, 0, 1, 1, 0, 1, 1, 1)

results <- mcmc_binary(obs, 100)

plot(results)
