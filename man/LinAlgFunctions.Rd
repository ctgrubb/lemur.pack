\name{LinAlgFunctions}
\alias{LinAlgFunctions}
\alias{tridiag}
\alias{det.tridiag}
\title{Linear Algebra Functions}
\usage{
tridiag(a, b, c)

det.tridiag(x)
}
\arguments{
\item{a}{A vector of size n, whose entries will form a diagonal. No default value.}

\item{b}{A vector of size n-1, whose entries will form a sub-diagonal. No default value.}

\item{c}{A vector of size n-1, whose entries will form a super-diagonal. Defaults is c = b.}

\item{x}{A square matrix.}
}
\description{
\loadmathjax
A select few functions to make constructing and calculating the determinant of tridiagnonal matrices more efficient.
}
\details{
To be completed later.
}
\examples{
tridiag(c(1, 2, 3, 4), c(5, 6, 7), c(8, 9, 10))
tridiag(c(1, 2, 3), c(4, 5), c(6, 7))

}
