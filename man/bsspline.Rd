% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsspline.R
\name{bsspline}
\alias{bsspline}
\title{Evaluate the functions \eqn{b} and \eqn{s} at \code{x}}
\usage{
bsspline(x, bsvec, alpha, d, n.ints, natural)
}
\arguments{
\item{x}{A value or vector of values at which the functions \eqn{b}
and \eqn{s} are to be evaluated}

\item{bsvec}{The \eqn{(2q-1)}-vector
\deqn{\big(b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h) \big),}
where \eqn{q}=ceiling(\eqn{d}/0.75) and \eqn{h=d/q}.
This vector specifies the CIUUPI, for all possible values of the random
error variance and the observed response vector}

\item{alpha}{The desired minimum coverage probability is \eqn{1 - \alpha}}

\item{d}{The functions \eqn{b} and \eqn{s} are specified by cubic splines
on the interval \eqn{[-d, d]}}

\item{n.ints}{The number of equal-length intervals in \eqn{[0, d]}, where
the endpoints of these intervals specify the knots,
belonging to \eqn{[0, d]}, of the cubic spline interpolations
that specify the functions \eqn{b} and \eqn{s}. In the description
of \code{bsvec}, \code{n.ints} is also called \eqn{q}.}

\item{natural}{Equal to 1 (default) if the \eqn{b} and \eqn{s} functions are obtained by
natural cubic spline interpolation or 0 if obtained by clamped cubic spline
interpolation}
}
\value{
A data frame containing \code{x} and the corresponding values of the
functions \eqn{b} and \eqn{s}.
}
\description{
Evaluate the functions \eqn{b} and \eqn{s}, as specified by
\code{bsvec},
\code{alpha}, \code{d}, \code{n.ints} and \code{natural}, at \code{x}.
}
\examples{
x <- seq(0, 8, by = 1)
alpha <- bs.list.example$alpha
natural <- bs.list.example$natural
d <- bs.list.example$d
n.ints <- bs.list.example$n.ints
bsvec <- bs.list.example$bsvec
bs <- bsspline(x, bsvec, alpha, d, n.ints, natural)

}
