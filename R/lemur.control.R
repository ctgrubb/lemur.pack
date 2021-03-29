#' Construct a control object used by lemur.pack functions
#'
#' @param nchains number of chains to construct
#' @param nsteps number of steps to take in each chain (excluding burnin)
#' @param burnin number of steps to take during burn-in
#' @param thin thinning interval; rounded to nearest integer
#' @param ... other arguments, typically used for prior hyperparameters
#' @export
lemur.control <- function(nchains = 1, nsteps = 1000, burnin = 100, thin = 1, a = 1, b = 1) {
  if(!is.numeric(nchains) || nchains <= 0)
    stop("Number of chains must be > 0")
  if(!is.numeric(nsteps) || nsteps <= 0)
    stop("Number of steps must be > 0")
  if(!is.numeric(burnin) || burnin <= 0)
    stop("Number of chains must be > 0")
  if(!is.numeric(thin) || thin <= 0)
    stop("Number of chains must be > 0")
  if(!is.numeric(a) || a <= 0)
    stop("Beta parameter a must be > 0")
  if(!is.numeric(b) || b <= 0)
    stop("Beta parameter b must be > 0")
  list(nchains = nchains, nsteps = nsteps, burnin = burnin, thin = thin, a = a, b = b)
}
