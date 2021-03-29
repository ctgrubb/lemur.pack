#' @title Sample from the posterior predictive of a multivariate multinomial distribution
#'
#' @description
#' \loadmathjax
#' \code{multivariate_multinomial_postpred} samples from the posterior predictive distribution given an observation from
#' a multivariate multinomial distribution and a prior on the probabilities of each category.
#'
#' @param obs a \code{data.frame} representing the observed sample.
#' @param N the size of the predictive samples to create.
#' @param nsamples the number of predictive samples to create.
#' @param priors either a list of named numeric vectors containing the Dirichlet parameters for individual priors on
#' each variable or an array describing the joint Dirichlet prior. Can also be a single number, which is interpreted
#' as a flat joint prior with the given height. See \code{vignette("mmpp", package = "lemur.pack")} for an example of
#' setting priors.
#'
#' @details
#' I will complete this later!
#'
#' @import dplyr
#' @export
multivariate_multinomial_postpred <- function(obs, N, nsamples, priors) {

  priors = 0.5
  prior_grid <- expand.grid(priors)
  colnames(prior_grid) <- paste0(colnames(prior_grid), "_Prior")
  prior_grid$Product_Prior <- apply(prior_grid, 1, prod)

  Ysamp <- array(NA, dim = c(ncol(obs), N, nsamples),
                 dimnames = list("Variable" = colnames(obs), "Individual" = 1:N, "Sample" = 1:nsamples))

  obs_factors <- as.data.frame(lapply(obs, as.factor))

  obs_df <- obs_factors %>%
    group_by(across(), .drop = FALSE) %>%
    summarize(Count = n(), .groups = "drop")

  df <- cbind(obs_df, prior_grid)
  df$Posterior <- df$Count + df$Product_Prior

  n <- nrow(df)

  P <- rdirichlet(nsamples, df$Posterior)
  Y <- t(apply(P, 1, function(p) sample(1:n, size = N, replace = TRUE, prob = p)))

  obs_df <- as.matrix(as.data.frame(lapply(obs_df, as.numeric)))

  for(i in 1:nsamples) {
      Ysamp[, , i] <- obs_df[Y[i, ], 1:ncol(obs)]
  }

  return(Ysamp)

}
