options(scipen = 5)

x <- c(0.14, 0.21, 0.23, 0.37, 0.41, 0.48, 0.59, 0.75, 0.84, 0.94)
N <- 100
samples <- 1000

dist_spacing <- function(x, N, samples) {
  xAug <- c(0, sort(x), 1)
  diffs <- diff(xAug)

  spacings <- sample(1:length(diffs), size = N * samples, replace = TRUE)
  lBound <- xAug[spacings]
  uBound <- xAug[spacings + 1]

  draws <- runif(N * samples, lBound, uBound)

  out <- matrix(draws, nrow = samples, ncol = N)
}

spacing_test <- dist_spacing(x, 100, 1000)
hist(spacing_test)
