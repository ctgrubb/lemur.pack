pop.others <- rnorm(90)
samp <- rnorm(10)
pop <- c(pop.others, samp)

ord.samp <- sort(samp)

# P(xmin <= observed min)
1 - choose(sum(pop > ord.samp[1]), 10) / choose(100, 10)
