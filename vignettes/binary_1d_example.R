## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message = FALSE, warning = FALSE-----------------------------------------
library(lemur.pack)
library(dplyr)
library(tidyr)
library(ggplot2)

## -----------------------------------------------------------------------------
p <- 0.5
obs <- c(1, 0, 0, 0, 0, 0, 1, 1, 0, 1)
N <- 100
n <- length(obs)

y <- sum(obs)

## -----------------------------------------------------------------------------
df <- data.frame(Y = y:(N-n+y), LogLikelihood = NA, LogPriorFlat = NA, LogPriorUninformative = NA)

for(i in 1:nrow(df)) {
  df$LogLikelihood[i] <- lemur.pack:::loglik_binary_(N = N, n = n, Y = df$Y[i], y = y)
}

## -----------------------------------------------------------------------------
for(i in 1:nrow(df)) {
  df$LogPriorFlat[i] <- lemur.pack:::logprior_binary_flat_()
}

## -----------------------------------------------------------------------------
for(i in 1:nrow(df)) {
  df$LogPriorUninformative[i] <- lemur.pack:::logprior_binary_uninformative_(N, Y = df$Y[i])
}

## -----------------------------------------------------------------------------
df <- df %>%
  pivot_longer(contains("LogPrior"), names_to = "Prior", names_prefix = "LogPrior", 
               values_to = "LogPrior") %>%
  mutate(LogPosterior = LogLikelihood + LogPrior) %>%
  mutate(
    Likelihood = exp(LogLikelihood),
    Posterior = exp(LogPosterior)
  )
df <- as.data.frame(df)

## ----echo = FALSE, fig.width = 6----------------------------------------------
  ggplot(data = df, mapping = aes(x = Y, y = LogPrior, color = Prior)) +
  geom_line() + 
  labs(x = "Y", y = "LogPrior")

## -----------------------------------------------------------------------------
flat_samples <- mcmc_binary(obs, N, prior = "flat", nsteps = 2000)
Y_samples <- data.frame(Y = rowSums(flat_samples))

## ----echo = FALSE, fig.width = 6----------------------------------------------
ggplot(data = Y_samples, mapping = aes(x = Y)) + 
  geom_density() + 
  scale_x_continuous(limits = c(0, N)) + 
  labs(y = "Density")

## -----------------------------------------------------------------------------
uninformative_samples <- mcmc_binary(obs, N, prior = "uninformative", nsteps = 2000)
Y_samples <- data.frame(Y = rowSums(uninformative_samples))

## ----echo = FALSE, fig.width = 6----------------------------------------------
ggplot(data = Y_samples, mapping = aes(x = Y)) + 
  geom_density() + 
  scale_x_continuous(limits = c(0, N)) + 
  labs(y = "Density")

