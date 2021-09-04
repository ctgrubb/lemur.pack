library(extraDistr)

table(ex1_df1[, -c(1, 2, 4)])

hist(pop[, 1])
table(pop[, -1])

hist(samp[, 1])
table(samp[, -1])

# We basically just have to fill out the table, satisfying the margins, and select pops with higher likelihood of getting our sample

basicSampler <- function(samp, pop, iters) {

  tab <- table(samp[, -1], exclude = FALSE)
  mar <- as.numeric(table(pop$HousePriceCat, exclude = FALSE))
  tal <- as.numeric(table(samp$HousePriceCat, exclude = FALSE))
  len <- length(table(samp$IncomeCat, exclude = FALSE))

  out <- array(NA, dim = c(iters, len, length(tal)))
  lprobs <- rep(0, iters)

  for(i in 1:iters) {
    out[i, , ] <- table(samp[, -1], exclude = FALSE)
    for(k in 1:length(mar)) {
      out[i, , k] <- out[i, , k] + rmultinom(1, mar[k] - tal[k], prob = rep(1 / len, len))
      lprobs[i] <- lprobs[i] + dmvhyper(tab[, k], out[i, , k], tal[k], log = TRUE)
    }
  }

  probs <- exp(lprobs)
  samps <- sample(1:iters, size = iters, replace = TRUE, prob = probs)

  out2 <- out[samps, , ]

}

basicSampler(samp, pop, 1000)


advSampler <- function(samp, pop, iters) {

  tab <- table(samp[, -1], exclude = FALSE)
  mar <- as.numeric(table(pop$HousePriceCat, exclude = FALSE))
  tal <- as.numeric(table(samp$HousePriceCat, exclude = FALSE))
  len <- length(table(samp$IncomeCat, exclude = FALSE))

  numsI <- matrix(NA, nrow = iters, ncol = nrow(pop))
  numsHP <- matrix(NA, nrow = iters, ncol = nrow(pop))

  cats <- array(NA, dim = c(iters, len, length(tal)))

  lprobs <- rep(0, iters)

  numsI[, 1:nrow(samp)] <- matrix(samp$Income, nrow = iters, ncol = nrow(samp), byrow = TRUE)


  for(i in 1:iters) {

    numsI[i, (nrow(samp) + 1):nrow(pop)] <- rnorm(nrow(pop) - nrow(samp), mean = mean(samp$Income), sd = sd(samp$Income))
    numsHP[i, ] <- sample(pop$HousePrice)

    catsI <- case_when(
      numsI[i, ] < 40000 ~ "<40000",
      numsI[i, ] < 65000 ~ "40000-64999",
      numsI[i, ] < 110000 ~ "65000-109999",
      TRUE ~ ">=110000"
    )

    catsHP <- case_when(
      numsHP[i, ] < 80000 ~ "<80000",
      numsHP[i, ] < 140000 ~ "80000-139999",
      numsHP[i, ] < 220000 ~ "140000-219999",
      numsHP[i, ] < 360000 ~ "220000-359999",
      TRUE ~ ">=360000"
    )

    cats[i, , ] <- table(catsI, catsHP)

    for(k in 1:dim(cats)[3]) {
      lprobs[i] <- lprobs[i] + dmvhyper(tab[, k], cats[i, , k], tal[k], log = TRUE)
    }

  }

  probs <- exp(lprobs)
  samps <- sample(1:iters, size = iters, replace = TRUE, prob = probs)

  out <- list(
    cats = cats[samps, , ],
    numsI = numsI[samps, ],
    numsHP = numsHP[samps, ]
  )

  return(out)

}

basicSampler(samp, pop, 1000)
test <- advSampler(samp, pop, 10000)


table(ex1_df1[, -c(1, 2, 4)])
apply(test$cats, 2:3, sum)
table(samp[, -1])

