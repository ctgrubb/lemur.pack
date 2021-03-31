#' @name Datasets
#' @rdname Datasets
#'
#' @title
#' Synthetic grouped two dimensional numeric data
#'
#' @description
#' A synthetically created two dimensional numeric dataset containing 5000 observations (in the population) and 100
#' observations (in the sample), as well as a grouping variable with three levels. The *true* correlations between the
#' two numeric variables are 0.6, 0.5, and 0.4 for group A, B, and C, respectively.
#'
#' @format
#' \code{grouped2d_pop} is a data frame with 5000 cases (rows) and 3 variables (columns) named \code{Group}, \code{X1},
#' and \code{X2}.
#'
#' @keywords datasets
"grouped2d_pop"

#' @rdname Datasets
#'
#' @format
#' \code{grouped2d_sample} is a data frame with 100 cases (rows) and 3 variables (columns) named \code{Group},
#' \code{X1}, and \code{X2}.
"grouped2d_sample"
