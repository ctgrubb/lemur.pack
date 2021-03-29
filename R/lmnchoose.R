#' @name Special
#' @rdname Special
#'
#' @title Calculate the multinomial coefficient.
#'
#' @description
#' \loadmathjax
#' \code{mnchoose} returns the multinomial coefficient.
#' \code{lmnchoose} returns the natural logarithm of the multinomial coefficient.
#'
#' @details
#' The multinomail coefficient is the
#' generalization of the binomial coefficient, and represents the number of ways it is possible to choose
#' \mjeqn{k_1}{k_1} objects of class 1, \mjeqn{k_2}{k_2} objects of class 2, up to \mjeqn{k_m}{k_m} objects of class m,
#' from a sample of size \code{n}.
#'
#' The multinomail coefficient can be expressed as
#' \mjdeqn{{n \choose \mathbb{k}} = \frac{n!}{k_1! k_2! ... k_m!}}{}
#'
#' This calculation is difficult (or impossible) for large \code{n}, so it is calculated in log space, which is much
#' easier. If the logarithm is not desired, the result is exponentiated, but beware this can also fail if \code{n} is
#' large enough.
#'
#' @param n number representing the size of the sample.
#' @param k numeric vector representing the number of times each category was chosen.
#' @export
mnchoose <- function(n, k) {
  exp(sum(log(1:n)) - sum(sapply(k, function(x) sum(log(1:x)))))
}

#' @rdname Special
#' @export
lmnchoose <- function(n, k) {
  sum(log(1:n)) - sum(sapply(k, function(x) sum(log(1:x))))
}
