library(dplyr)
library(lemur.pack)
library(extraDistr)

pop <- ex1pop
samp <- ex1samp

table(pop[, -c(1, 2, 4)])
hist(pop[, 2])

table(samp[, -c(1, 2, 4)])
hist(samp[, 2])

sample <- samp[, -c(1, 4)]
population <- pop[, -c(1, 2, 3)]

iterations <- 1000

mcmcSampler <- function(sample, population, iterations) {

  tab <- table(sample[, -1], exclude = FALSE)
  mar <- as.numeric(table(population$HousePriceCat, exclude = FALSE))
  tal <- as.numeric(table(sample$HousePriceCat, exclude = FALSE))

  nCatsI <- length(table(sample$IncomeCat, exclude = FALSE))
  nCatsH <- length(table(sample$HousePriceCat, exclude = FALSE))

  weights <- rep(0, iterations)

  out <- list()

  fulldata <- list()
  incomes <- matrix(NA, nrow = iterations, ncol = nrow(population))
  hprices <- matrix(NA, nrow = iterations, ncol = nrow(population))
  cats <- array(NA, dim = c(iterations, nCatsI, nCatsH))

  for(i in 1:iterations) {
    # Get sample
    unique <- 0
    while(unique != 100) {
      init <- left_join(population, sample, by = "HousePriceCat") %>%
        group_by(Income) %>%
        sample_n(1)

      unique <- length(unique(init$HousePrice))
    }

    init2 <- anti_join(population, init, by = "HousePrice")

    init2$Income <- rnorm(nrow(init2), mean(sample$Income) + sd(sample$Income) / sd(population$HousePrice) *
                            cor(init$Income, init$HousePrice) * (init2$HousePrice - mean(population$HousePrice)),
                          sqrt((1 - cor(init$Income, init$HousePrice)) ** 2) * sd(sample$Income))
    init2 <- init2 %>%
      mutate(
        IncomeCat = case_when(
          Income < 40000 ~ "<40000",
          Income < 65000 ~ "40000-64999",
          Income < 110000 ~ "65000-109999",
          TRUE ~ ">=110000"
        )
      ) %>%
      mutate(IncomeCat = factor(IncomeCat, levels = c("<40000", "40000-64999", "65000-109999", ">=110000")))

    init <- bind_rows(init, init2)
    rm(init2)

    fulldata[[i]] <- init
    incomes[i, ] <- init$Income
    hprices[i, ] <- init$HousePrice
    cats[i, , ] <- table(init$IncomeCat, init$HousePriceCat)

    for(k in 1:nCatsH) {
      weights[i] <- weights[i] + dmvhyper(tab[, k], cats[i, , k], tal[k], log = TRUE)
    }

  }

  weights <- exp(weights)
  samps <- sample(1:iterations, replace = TRUE, prob = weights)

  out <- list(
    fulldata = fulldata[samps],
    incomes = incomes[samps, ],
    hprices = hprices[samps, ],
    cats = cats[samps, , ]
  )

}

apply(out$cats, c(2, 3), sum) / 100

table(samp[, -c(1, 2, 4)])
table(pop$HousePriceCat)

apply(out$cats, c(2, 3), min)
apply(out$cats, c(2, 3), max)

table(pop[, -c(1, 2, 4)])


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

