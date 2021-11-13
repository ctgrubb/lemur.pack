#' @name Special
#' @rdname Special
#'
#' @title (More) Special Functions of Mathematics
#'
#' @description
#' \loadmathjax
#' Special mathematical functions not included in base R, specifically the multivariate form of the beta function and the
#' multinomial coefficient.
#'
#' @usage
#' mvbeta(a)
#' lmvbeta(a)
#'
#' mnchoose(n, k)
#' lmnchoose(n, k)
#'
#' @details
#' The functions \code{mvbeta} and \code{lmvbeta} return the multivariate beta function and the natural logarithm of the
#' multivariate beta function,
#' \mjdeqn{\beta(a) = \Gamma(a_1)\Gamma(a_2) ... \Gamma(a_m)/\Gamma(a_1 + a_2 + \dots + a_m).}{}
#'
#' The functions \code{mnchoose} and \code{lmnchoose} return the multinomial coefficient and the natural logarithm of the
#' multinomial coefficent,
#' \mjdeqn{{n \choose \mathbb{k}} = \frac{n!}{k_1! k_2! \dots k_m!}.}{}
#'
#' These calculations are difficult (or impossible) for large \code{n}, so it is calculated in log space, which is much
#' easier. If the logarithm is not desired, the result is exponentiated (and rounded, because
#' \mjeqn{{n \choose \mathbb{k}}}{} should be an integer), but beware this can also fail if \code{n} is large enough.
#'
#' @param a non-negative numeric vector.
#' @param n non-negative integer.
#' @param k non-negative integer vector.
#'
#' @references
#' \url{https://en.wikipedia.org/wiki/Beta_function#Multivariate_beta_function}
#' \url{https://en.wikipedia.org/wiki/Multinomial_theorem#Multinomial_coefficients}
#'
#' @seealso
#' See \code{\link[base]{Special}} for the univariate beta function and the binomial coefficient.
#'
#' @examples
#' mvbeta(c(0.5, 1, 0.5, 2))
#' lmvbeta(c(0.5, 1, 0.5, 2))
#'
#' mnchoose(12, c(3, 3, 4, 2))
#' lmnchoose(12, c(3, 3, 4, 2))
#'
#' @export
mvbeta <- function(a) {
  prod(gamma(a)) / gamma(sum(a))
}

#' @rdname Special
#' @export
lmvbeta <- function(a) {
  sum(lgamma(a)) - lgamma(sum(a))
}

#' @rdname Special
#' @export
mnchoose <- function(n, k) {
  round(exp(sum(log(1:n)) - sum(sapply(k, function(x) sum(log(1:x))))))
}

#' @rdname Special
#' @export
lmnchoose <- function(n, k) {
  k[k == 0] <- 1
  sum(log(1:n)) - sum(sapply(k, function(x) sum(log(1:x))))
}
