#' @name LinAlgFunctions
#' @rdname LinAlgFunctions
#'
#' @title Linear Algebra Functions
#'
#' @description
#' \loadmathjax
#' A select few functions to make constructing and calculating the determinant of tridiagnonal matrices more efficient.
#'
#' @usage
#' tridiag(a, b, c)
#' tridiag(a, b)
#'
#' det.tridiag(x)
#'
#' @details
#' To be completed later.

#' @param a A vector of size n, whose entries will form a diagonal. No default value.
#' @param b A vector of size n-1, whose entries will form a sub-diagonal. No default value.
#' @param c A vector of size n-1, whose entries will form a super-diagonal. Defaults is c = b.
#'
#' @examples
#' tridiag(c(1, 2, 3, 4), c(5, 6, 7), c(8, 9, 10))
#' tridiag(c(1, 2, 3), c(4, 5), c(6, 7))
#'
#' @export
tridiag = function(a, b, c = b) {
  if(length(a) == 1 & (missing(b) | missing(c))) {
    return(matrix(a))
  } else if(length(a) == 2){
    return(matrix(c(a[1], b[1], c[1], a[2]), 2))
  } else {
    M = diag(a)
    diag(M[-1, ]) = b[1:(length(a) - 1)]
    diag(M[, -1]) = c[1:(length(a) - 1)]
    return(M)
  }
}

#' @rdname LinAlgFunctions
#' @export
det.tridiag = function(x) {
  m = dim(x)[1]
  val = numeric(m)
  val[1] = x[1, 1]
  val[2] = x[2, 2] * val[1] - x[2, 1] ** 2
  for(i in 3:m) {
    val[i] = x[i, i] * val[i - 1] - x[i, i - 1] ** 2 * val[i - 2]
  }
  return(val[m])
}
