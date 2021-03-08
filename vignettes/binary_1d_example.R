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
p <- 0.5; n <- 10; N <- 100
obs <- sample(c(0, 1), size = n, replace = TRUE, prob = c(1-p, p))
y <- sum(obs)
y <- 5

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

## ---- echo = FALSE------------------------------------------------------------
  ggplot(data = df, mapping = aes(x = Y, y = LogPrior, color = Prior)) +
  geom_line()

