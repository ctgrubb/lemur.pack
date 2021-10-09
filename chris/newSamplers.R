library(tidyverse)
library(proftools)
library(lemur.pack)
library(extraDistr)

data("ex1pop")
data("ex1samp")

pop <- ex1pop %>%
  select(HousePrice, HousePriceCat)

samp <- ex1samp %>%
  select(Income, IncomeCat, HousePriceCat)

indNormalMCMC <- function(sample, population, iterations = 1000) {

  ranges <- list(
    mins = list(
      Income = c(0, 40000, 65000, 110000),
      HousePrice = c(0, 80000, 140000, 220000, 360000)
    ),
    maxs = list(
      Income = c(39999, 64999, 109999, Inf),
      HousePrice = c(79999, 139999, 219999, 359999, Inf)
    )
  )

  modTable <- function(tabs, old, new) {
    oldElem <- matrix(0, nrow = 1, ncol = length(old))
    for(i in 1:length(oldElem)) {
      oldElem[1, i] <- findInterval(old[i], vec = c(ranges$mins[[i]], Inf), all.inside = TRUE)
    }

    newElem <- matrix(0, nrow = 1, ncol = length(new))
    for(i in 1:length(newElem)) {
      newElem[1, i] <- findInterval(new[i], vec = c(ranges$mins[[i]], Inf), all.inside = TRUE)
    }

    tabs[oldElem] <- tabs[oldElem] - 1
    tabs[newElem] <- tabs[newElem] + 1

    return(tabs)
  }

  toLLik <- function(tabs) {
    out <- 0
    for(i in 1:nCatsH) {
      out <- out + suppressWarnings(dmvhyper(tab[, i], tabs[, i], tal[i], log = TRUE))
    }
    return(out)
  }

  l1 <- levels(sample$IncomeCat)
  l2 <- levels(sample$HousePriceCat)

  tab <- table(sample[, -1], exclude = FALSE)
  mar <- as.numeric(table(population$HousePriceCat, exclude = FALSE))
  tal <- as.numeric(table(sample$HousePriceCat, exclude = FALSE))

  nCatsI <- length(l1)
  nCatsH <- length(l2)

  N <- nrow(population)
  mu <- 100000
  sig <- 30000
  step <- 30000

  current <- data.frame(
    Income = rnorm(N, mu, sig),
    HousePrice = sample(population$HousePrice)
  )

  currentMat <- matrix(0, nrow = nCatsI, ncol = nCatsH)

  while(!is.finite(toLLik(currentMat))) {
    current$Income <- rnorm(N, mu, sig)

    currentMat <- current %>%
      mutate(
        IncomeCat = case_when(
          Income < 40000 ~ "<40000",
          Income < 65000 ~ "40000-64999",
          Income < 110000 ~ "65000-109999",
          TRUE ~ ">=110000"
        ),
        HousePriceCat = case_when(
          HousePrice < 80000 ~ "<80000",
          HousePrice < 140000 ~ "80000-139999",
          HousePrice < 220000 ~ "140000-219999",
          HousePrice < 360000 ~ "220000-359999",
          TRUE ~ ">=360000"
        )
      ) %>%
      transmute(
        IncomeCat = factor(IncomeCat, levels = c("<40000", "40000-64999", "65000-109999", ">=110000")),
        HousePriceCat = factor(HousePriceCat, levels = c("<80000", "80000-139999", "140000-219999", "220000-359999", ">=360000"))
      ) %>%
      table(., exclude = FALSE) %>%
      matrix(., ncol = nCatsH)
  }

  out <- list()

  for(i in 1:iterations) {
    propDf <- current
    for(k in 1:N) {
      proposal <- rtnorm(1, current[k, 1], step, 0, Inf)
      propDf[k, 1] <- proposal
      propMat <- modTable(currentMat, current[k, ], propDf[k, ])

      # compare
      alpha <- toLLik(propMat) + dnorm(proposal, mu, sig, log = TRUE) -
        toLLik(currentMat) - dnorm(current[k, 1], mu, sig, log = TRUE) +
        dtnorm(current[k, 1], proposal, step, 0, Inf, log = TRUE) - dtnorm(proposal, current[k, 1], step, 0, Inf, log = TRUE)
      if(is.finite(alpha) & log(runif(1)) < alpha) {
        # Acceptance
        current <- propDf
        currentMat <- propMat
      } else {
        # Rejection
        propDf <- current
      }
    }
    out[[i]] <- current
  }
  return(out)
}


samples_1 <- indNormalMCMC(samp, pop, iterations = 1000)

