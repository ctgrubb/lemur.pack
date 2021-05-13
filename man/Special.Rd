\name{Special}
\alias{Special}
\alias{mvbeta}
\alias{lmvbeta}
\alias{mnchoose}
\alias{lmnchoose}
\title{(More) Special Functions of Mathematics}
\usage{
mvbeta(a)
lmvbeta(a)

mnchoose(n, k)
lmnchoose(n, k)
}
\arguments{
\item{a}{non-negative numeric vector.}

\item{n}{non-negative integer.}

\item{k}{non-negative integer vector.}
}
\description{
\loadmathjax
Special mathematical functions not included in base R, specifically the multivariate form of the beta function and the
multinomial coefficient.
}
\details{
The functions \code{mvbeta} and \code{lmvbeta} return the multivariate beta function and the natural logarithm of the
multivariate beta function,
\mjdeqn{\beta(a) = \Gamma(a_1)\Gamma(a_2) ... \Gamma(a_m)/\Gamma(a_1 + a_2 + \dots + a_m).}{}

The functions \code{mnchoose} and \code{lmnchoose} return the multinomial coefficient and the natural logarithm of the
multinomial coefficent,
\mjdeqn{{n \choose \mathbb{k}} = \frac{n!}{k_1! k_2! \dots k_m!}.}{}

These calculations are difficult (or impossible) for large \code{n}, so it is calculated in log space, which is much
easier. If the logarithm is not desired, the result is exponentiated (and rounded, because
\mjeqn{{n \choose \mathbb{k}}}{} should be an integer), but beware this can also fail if \code{n} is large enough.
}
\examples{
mvbeta(c(0.5, 1, 0.5, 2))
lmvbeta(c(0.5, 1, 0.5, 2))

mnchoose(12, c(3, 3, 4, 2))
lmnchoose(12, c(3, 3, 4, 2))

}
\references{
\url{https://en.wikipedia.org/wiki/Beta_function#Multivariate_beta_function}
\url{https://en.wikipedia.org/wiki/Multinomial_theorem#Multinomial_coefficients}
}
\seealso{
See \code{\link[base]{Special}} for the univariate beta function and the binomial coefficient.
}
