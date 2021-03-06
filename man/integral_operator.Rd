% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{integral_operator}
\alias{integral_operator}
\title{Integral transformation of a curve using an integral operator}
\usage{
integral_operator(operator_kernel, curve, v)
}
\arguments{
\item{operator_kernel}{Matrix with the values
of the kernel surface of the integral operator.
The dimension of the matrix is \eqn{(g x m)}, where
\eqn{g} is the number of discretization points of the
input curve and \eqn{m} is the number of discretization
points of the output curve.}

\item{curve}{Vector containing the discretized values
of a functional observation. The dimension of the
matrix is \eqn{(1 x m)}, where \eqn{m} is the
number of points observed in the curve.}

\item{v}{Numerical vector specifying the
discretization points of the curves.}
}
\value{
Returns a matrix the same size as
\code{curve} with the transformed values.
}
\description{
Compute the integral transform of the
curve \eqn{Y_i} with respect to a given integral
operator \eqn{\Psi}. The transformation is given by
\deqn{\Psi(Y_{i})(v) = \int \psi(u,v)Y_{i}(u)du}
}
\examples{
# Example 1

v <- seq(from = 0, to = 1, length.out = 20)
set.seed(10)
curve <- sin(v) + rnorm(length(v))
operator_kernel <- 0.6*(v \%*\% t(v))
hat_curve <- integral_operator(operator_kernel,curve,v)

}
