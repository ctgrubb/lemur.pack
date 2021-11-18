# Idea from 11/11/2021 meeting

iters <- 10000
N <- 100
n <- 10

pop <- rnorm(N)
samp <- sort(sample(pop, n))

weight.calc <- function(synth.pop, ord.sample) {
  tabs <- table(findInterval(synth.pop, c(-Inf, samp, Inf), left.open = TRUE))
  add1 <- as.numeric(tabs)
  names(add1) <- names(tabs)
  add2 <- rep(0, n + 1)
  names(add2) <- 1:(n + 1)
  add <- list(add1, add2)

  tab <- tapply(unlist(add), names(unlist(add)), sum)
  tab <- tab[order(factor(names(tab), levels = 1:(length(ord.sample) + 1)))]

  F <- tab / N

  st <- tab * log(F)
  st[!is.finite(st)] <- 0
  return(exp(sum(st)))
  # return(tab)
}

synth.pops <- matrix(0, nrow = iters, ncol = N)
for(i in 1:iters) {
  synth.pops[i, ] <- c(samp, rnorm(N - n))
}

synth.pop <- synth.pops[1, ]
ord.sample <- samp

plot(synth.pop, rep(0, N))
abline(v = ord.sample)

weights <- apply(synth.pops, 1, weight.calc, ord.sample = samp)
index <- sample(1:iters, size = iters / 10, replace = TRUE, prob = weights)

table(index)
synth.pop <- synth.pops[X, ]

