library(lemur.pack)
library(extraDistr)
library(dplyr)

data("ex1pop", "ex1samp")

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

pop <- ex1pop %>%
  slice_sample(n = 100)

samp <- pop %>%
  slice_sample(n = 20) %>%
  select(
    Group, IncomeCat, HousePriceCat
  )

# table(pop[, c(3, 5)])
# table(samp[, -1])

# sample <- samp
# population <- pop
# iterations <- 1000
# proposalDist <- "unifMix"
# adjust <- TRUE
# a <- 0
# b <- 200000

newTry <- function(sample, population, ranges, iterations = 1000, proposalDist = "mixUnif",
                   adjust = TRUE, ...) {

  params <- list(...)
  if(!is.null(params$a)) a <- params$a
  if(!is.null(params$b)) b <- params$b
  if(!is.null(params$mu)) mu <- params$mu
  if(!is.null(params$sd)) sd <- params$sd

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

  getLLik <- function(tabs) {
    out <- 0
    for(i in 1:nCatsH) {
      out <- out + suppressWarnings(dmvhyper(tab[, i], tabs[, i], tal[i], log = TRUE))
    }
    return(out)
  }

  corrFactor <- function(currentMat, propMat) {
    if(all(currentMat == propMat)) {
      return(0)
    } else {
      ind <- as.numeric(which(currentMat != propMat, arr.ind = TRUE)[1, 2])
      return(lmnchoose(mar[ind], propMat[, ind]) - lmnchoose(mar[ind], currentMat[, ind]))
    }
  }

  l1 <- levels(sample$IncomeCat)
  l2 <- levels(sample$HousePriceCat)

  tab <- table(sample[, -1], exclude = FALSE)
  mar <- as.numeric(table(population$HousePriceCat, exclude = FALSE))
  tal <- as.numeric(table(sample$HousePriceCat, exclude = FALSE))

  nCatsI <- length(l1)
  nCatsH <- length(l2)

  N <- nrow(population)

  current <- data.frame(
    Income = rtnorm(N, 100000, 30000, 0, Inf),
    HousePrice = sample(population$HousePrice)
  )

  currentMat <- matrix(0, nrow = nCatsI, ncol = nCatsH)

  while(!is.finite(getLLik(currentMat))) {
    current$Income <- rtnorm(N, 100000, 30000, 0, Inf)

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

  incomes <- matrix(NA, nrow = iterations, ncol = N)
  membership <- array(NA, dim = c(iterations, nCatsI, nCatsH))
  llik <- rep(0, iterations)
  acceptance <- rep(FALSE, iterations * N)

  for(i in 1:iterations) {
    propDf <- current
    for(k in 1:N) {

      # Proposal switch
      if(proposalDist == "unifMix") {
        pick <- sample(1:nCatsI, 1)
        proposal <- case_when(
          pick == 1 ~ runif(1, min = a, max = 40000),
          pick == 2 ~ runif(1, min = 40000, max = 65000),
          pick == 3 ~ runif(1, min = 65000, max = 110000),
          pick == 4 ~ runif(1, min = 110000, max = b)
        )
        hastingsAdj_1 <- case_when(
          current[k, 1] < 40000 ~ log(1) - log(40000 - a),
          current[k, 1] < 65000 ~ log(1) - log(25000),
          current[k, 1] < 110000 ~ log(1) - log(45000),
          TRUE ~ log(1) - log(b - 110000)
        )
        hastingsAdj_2 <- case_when(
          proposal < 40000 ~ log(1) - log(40000 - a),
          proposal < 65000 ~ log(1) - log(25000),
          proposal < 110000 ~ log(1) - log(45000),
          TRUE ~ log(1) - log(b - 110000)
        )
        hastingsAdj <- hastingsAdj_1 - hastingsAdj_2
      } else if(proposalDist == "unif") {
        proposal <- runif(1, min = a, max = b)
        hastingsAdj <- 0
      } else if(proposalDist == "normal") {
        proposal <- rtnorm(1, mu = mu, sd = sd, a = a, b = b)
        hastingsAdj <- dtnorm(current[k, 1], proposal, mu = mu, sd = sd, a = a, b = b, log = TRUE) -
          dtnorm(proposal, current[k, 1], mu = mu, sd = sd, a = a, b = b, log = TRUE)
      }

      propDf[k, 1] <- proposal
      propMat <- modTable(currentMat, current[k, ], propDf[k, ])

      # Calculate Acceptance Probability
      if(adjust) {
        alpha <- getLLik(propMat) - getLLik(currentMat) + hastingsAdj + corrFactor(currentMat, propMat)
      } else {
        alpha <- getLLik(propMat) - getLLik(currentMat) + hastingsAdj
      }

      # Accept / Reject
      if(is.finite(alpha) & log(runif(1)) < alpha) {
        # Acceptance
        acceptance[(i - 1) * N + k] <- TRUE
        current <- propDf
        currentMat <- propMat
      } else {
        # Rejection
        propDf <- current
      }

    }

    # Storage
    incomes[i, ] <- current$Income
    membership[i, , ] <- currentMat
    llik[i] <- getLLik(currentMat)
  }


  out <- list(
    incomes = incomes,
    houseprices = current$HousePrice,
    membership = membership,
    acceptance = acceptance,
    loglikelihood = llik
  )
  return(out)
}

unif_runs <- newTry(sample = samp, population = pop, ranges = ranges, iterations = 10000, proposalDist = "unif", adjust = TRUE, a = 0, b = 200000)
unif_runs_f <- newTry(sample = samp, population = pop, ranges = ranges, iterations = 10000, proposalDist = "unif", adjust = FALSE, a = 0, b = 200000)

unifMix_runs <- newTry(sample = samp, population = pop, ranges = ranges, iterations = 10000, proposalDist = "unifMix", adjust = TRUE, a = 0, b = 200000)


c(40000, 25000, 45000, 90000) / 200000 * 2
c(40000, 25000, 45000, 90000) / 200000 * 8

apply(unif_runs$membership, c(2, 3), mean)
apply(unif_runs_f$membership, c(2, 3), mean)

apply(unifMix_runs$membership, c(2, 3), mean)

table(samp[, -1])
apply(membership, c(2, 3), mean, na.rm = TRUE)
mean(acceptance)
