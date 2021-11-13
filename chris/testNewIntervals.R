library(lemur.pack)
library(extraDistr)
library(dplyr)
library(ggplot2)

data("ex1pop")
data("ex1samp")

population <- ex1pop %>%
  select(HousePrice, HousePriceCat)

sample <- ex1samp %>%
  select(Income, IncomeCat, HousePriceCat)

iterations <- 1000

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

  out1 <- out2 <- out3 <- array(NA, dim = c(iterations, nCatsI, nCatsH))
  currentMat1 <- currentMat2 <- currentMat3 <- currentMat
  current1 <- current2 <- current3 <- current

  for(i in 1:iterations) {
    propDf1 <- current1
    propDf2 <- current2
    propDf3 <- current3
    for(k in 1:N) {

      proposal1 <- rtnorm(1, current1[k, 1], step, 0, Inf)
      propDf1[k, 1] <- proposal1
      propMat1 <- modTable(currentMat1, current1[k, ], propDf1[k, ])

      proposal2 <- rtnorm(1, current2[k, 1], step, 0, Inf)
      propDf2[k, 1] <- proposal2
      propMat2 <- modTable(currentMat2, current2[k, ], propDf2[k, ])

      proposal3 <- rtnorm(1, current3[k, 1], step, 0, Inf)
      propDf3[k, 1] <- proposal3
      propMat3 <- modTable(currentMat3, current3[k, ], propDf3[k, ])

      # calculate alpha
      alpha1 <- toLLik(propMat1) + dnorm(proposal1, mu, sig, log = TRUE) -
        toLLik(currentMat1) - dnorm(current1[k, 1], mu, sig, log = TRUE) +
        dtnorm(current1[k, 1], proposal1, step, 0, Inf, log = TRUE) -
        dtnorm(proposal1, current1[k, 1], step, 0, Inf, log = TRUE)

      alpha2 <- dnorm(proposal2, mu, sig, log = TRUE) -
        dnorm(current2[k, 1], mu, sig, log = TRUE) +
        dtnorm(current2[k, 1], proposal2, step, 0, Inf, log = TRUE) -
        dtnorm(proposal2, current2[k, 1], step, 0, Inf, log = TRUE)

      alpha3 <- toLLik(propMat3) -
        toLLik(currentMat3) +
        dtnorm(current3[k, 1], proposal3, step, 0, Inf, log = TRUE) -
        dtnorm(proposal3, current3[k, 1], step, 0, Inf, log = TRUE)

      # compare
      if(is.finite(alpha1) & log(runif(1)) < alpha1) {
        # Acceptance
        current1 <- propDf1
        currentMat1 <- propMat1
      } else {
        # Rejection
        propDf1 <- current1
      }

      if(is.finite(alpha2) & log(runif(1)) < alpha2) {
        # Acceptance
        current2 <- propDf2
        currentMat2 <- propMat2
      } else {
        # Rejection
        propDf2 <- current2
      }

      if(is.finite(alpha3) & log(runif(1)) < alpha3) {
        # Acceptance
        current3 <- propDf3
        currentMat3 <- propMat3
      } else {
        # Rejection
        propDf3 <- current3
      }

    }

    # Storage
    out1[i, , ] <- currentMat1
    out2[i, , ] <- currentMat2
    out3[i, , ] <- currentMat3

  }
  return(list(out1, out2, out3))
}

system.time(
  overnightRun <- indNormalMCMC(sample, population, iterations = 20000)
)

samps <- data.frame(Dist = "Sample", as.data.frame(table(sample[, -1])))

means1 <- apply(overnightRun[[1]][-(1:5000), , ], c(2, 3), mean)
dimnames(means1) <- dimnames(table(sample[, -1]))
means1 <- data.frame(Dist = "Combined", as.data.frame(as.table(means1)))

means2 <- apply(overnightRun[[2]][-(1:5000), , ], c(2, 3), mean)
dimnames(means2) <- dimnames(table(sample[, -1]))
means2 <- data.frame(Dist = "Prior", as.data.frame(as.table(means2)))

means3 <- apply(overnightRun[[3]][-(1:5000), , ], c(2, 3), mean)
dimnames(means3) <- dimnames(table(sample[, -1]))
means3 <- data.frame(Dist = "Likelihood", as.data.frame(as.table(means3)))

df <- bind_rows(samps, means1, means2, means3) %>%
  mutate(
    HJust = case_when(
      Dist == "Sample" ~ "right",
      Dist == "Combined" ~ "left",
      Dist == "Prior" ~ "right",
      Dist == "Likelihood" ~ "left"
    ),
    VJust = case_when(
      Dist == "Sample" ~ "bottom",
      Dist == "Combined" ~ "bottom",
      Dist == "Prior" ~ "top",
      Dist == "Likelihood" ~ "top"
    ),
    Freq = sprintf("%05.1f", round(Freq, digits = 1)),
    Dist = factor(Dist, levels = c("Sample", "Combined", "Prior", "Likelihood"))
  )

ggplot(data = df, mapping = aes(x = HousePriceCat, y = IncomeCat)) +
  geom_label(
    mapping = aes(color = Dist, label = Freq, hjust = HJust, vjust = VJust),
    label.padding = unit(0.5, "lines")
  ) +
  scale_y_discrete(limits = rev(levels(df$IncomeCat))) +
  scale_color_discrete(name = "Distribution") +
  labs(
    x = "House Price (Categorical)",
    y = "Income (Categorical)",
    title = "Mean Membership"
  )

