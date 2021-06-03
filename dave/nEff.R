# Functional implementation of Dave's unifPopDraws code
library(extraDistr)

n <- 41; nEff <- 41; nSamples = 100; alpha = 1; plot = TRUE

rDirUnif <- function(n, nEff, nSamples = 1, alpha = 0.5, plot = TRUE) {
  tmp <- rdirichlet(nSamples, (nEff + alpha) / (n + 1) * rep(1, n + 1))
  tmp2 <- apply(tmp[, -(n + 1)], 1, cumsum)
  if(!plot) {
    return(tmp2)
  } else {
    tmp3 <- apply(tmp2 < .5, 2, mean)
    par(mfrow = c(1, 2))

    plot(c(0, 1), c(0, 1), type = 'l', xlab = 'u', ylab = 'Cumulative Probability')
    matlines(tmp2, cumsum(rep(1/n, n)), lty = 1, type = 's')
    mtext(paste0('Draws of Population CDF (nEff = ', nEff, ')'))

    hist(tmp3, probability = TRUE, xlim = c(0, 1), breaks = seq(0 - .5 / (nEff + 1), 1 + .5 / (nEff + 1), length = 40),
         main = "", xlab = "P(u < 0.5)")
    x <- seq(0, 1, length = 100)
    dx <- dnorm(x, .5, sqrt(.5 * (1 - .5)) / sqrt(nEff))
    lines(x, dx, col = 'blue')
    mtext(paste0('Draws of Population Proportion < 0.5 (nEff = ', nEff, ')'))
  }
}

samples_41_5 <- rDirUnif(41, 4, 1000, 0.5, plot = TRUE)
