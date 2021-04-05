library(parallel)
library(doFuture)
library(ggplot2)

source("dave/srs.R")

clust <- makeCluster(6)


# Test sample 1
obs <- c(1, 0, 1, 0, 1, 1, 1)
N <- 100
nsamples <- 100000

test_srs2 <- srs2(obs, N, nsamples)
test_srs3 <- srs3(obs, N, nsamples)
test_bb <- beta_binomial(obs, N, nsamples)

df_post <- data.frame(
  Y <- seq(0, 100, length.out = 1000),
  Density = dbeta(seq(0, 1, length.out = 1000), 1 + sum(obs), 1 + length(obs) - sum(obs)) / 100
)
df_srs2 <- data.frame(
  Method = "srs2",
  Y = rowSums(test_srs2)
)
df_srs3 <- data.frame(
  Method = "srs3",
  Y = rowSums(test_srs3)
)
df_bb <- data.frame(
  Method = "beta-binomial",
  Y = rowSums(test_bb)
)

df <- bind_rows(df_srs2, df_srs3, df_bb)

ggplot(data = df, mapping = aes(x = Y)) +
  geom_density(mapping = aes(x = Y, color = Method)) +
  geom_line(data = df_post, mapping = aes(x = Y, y = Density))
