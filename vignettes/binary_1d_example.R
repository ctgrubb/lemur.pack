## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(lemur.pack)
library(ggplot2)

## -----------------------------------------------------------------------------
p <- 0.5; n <- 10; N <- 100
obs <- sample(c(0, 1), size = n, replace = TRUE, prob = c(1-p, p))
y <- sum(obs)

## -----------------------------------------------------------------------------
df <- data.frame(Y = y:(N-n+y), LogLikelihood = NA, LogPrior = NA)

for(i in 1:nrow(df)) {
  df$LogLikelihood[i] <- lemur.pack:::loglik_binary_(N = N, n = n, Y = df$Y[i], y = y)
}


## -----------------------------------------------------------------------------

for(i in 1:nrow(df)) {
  df$LogPrior[i] <- lemur.pack:::logprior_binary_flat_()
}


## -----------------------------------------------------------------------------
df$LogPosterior <- df$LogLikelihood + df$LogPrior

ggplot(data = df, mapping = aes(x = Y, y = LogPosterior)) +
  geom_line()

