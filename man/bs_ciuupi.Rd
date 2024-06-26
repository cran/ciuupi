% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs_ciuupi.R
\name{bs_ciuupi}
\alias{bs_ciuupi}
\title{Computes the the functions \eqn{b} and \eqn{s}
that specify the CIUUPI for all possible
values of \eqn{\sigma} and the observed response vector}
\usage{
bs_ciuupi(alpha, rho, natural = 1)
}
\arguments{
\item{alpha}{The desired minimum coverage probability
is \eqn{1-\alpha}}

\item{rho}{The known correlation \eqn{\rho} between
\eqn{\widehat{\theta}} and \eqn{\widehat{\tau}}}

\item{natural}{Equal to 1 (default) if the functions \eqn{b}
and \eqn{s} are specified by natural cubic spline interpolation
or 0 if these functions are specified by clamped cubic
spline interpolation in an interval \eqn{[-d, d]}, where \eqn{d}
is computed by \code{bs_ciuupi} using a specified
function of
\code{alpha} and \code{rho}}
}
\value{
A list with the following components.

alpha, rho, natural: the inputs

d: a 'goldilocks' value of \eqn{d} that is not too large
and not too small

n.ints: number of equal-length consecutive
intervals whose union is \eqn{[0,d]},
this is the same as \eqn{q}

lambda.star: the computed value of \eqn{\lambda^*}

bsvec: the vector
\eqn{\big(b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h) \big)} that determines
the functions \eqn{b} and \eqn{s} that specify the CIUUPI for all possible
values of \eqn{\sigma} and observed response vector

comp.time: the computation time in seconds
}
\description{
Chooses the positive number \eqn{d} and the positive integer \eqn{q}, sets
\eqn{h=d/q}, and then computes the
\eqn{(2q-1)}-vector
\eqn{\big(b(h),...,b((q-1)h),
s(0),s(h)...,s((q-1)h)\big)}
that determines, via cubic spline interpolation, the functions
\eqn{b} and \eqn{s} which specify
the confidence interval for \eqn{\theta}
that utilizes the uncertain prior information (CIUUPI),
for all possible values of \eqn{\sigma} and the observed response vector.
To an excellent approximation, this confidence interval
has minimum coverage probability
\eqn{1-\alpha}.
}
\details{
Suppose that \deqn{y = X \beta + \varepsilon} where \eqn{y}
is a random \eqn{n}-vector of
responses, \eqn{X} is a known \eqn{n} by \eqn{p} matrix
with linearly
independent columns, \eqn{\beta} is an unknown parameter
\eqn{p}-vector and
\eqn{\varepsilon} is the random error with components that are iid normally distributed
with zero mean and known variance \eqn{\sigma^2}.
The parameter of interest is
\eqn{\theta = a^{\top} \beta}.
Also let \eqn{\tau = c^{\top}\beta -t}, where \eqn{a}
and \eqn{c} are specified linearly independent
vectors and \eqn{t} is a specified number.
The uncertain prior information is that \eqn{\tau = 0}.

Let \code{rho} denote the known
correlation between the \eqn{\widehat{\theta}} and \eqn{\widehat{\tau}}.
We can compute \code{rho}
from given values of \eqn{a}, \eqn{c} and \eqn{X}
using the function \code{acX_to_rho}.

The confidence interval for \eqn{\theta},
with minimum coverage probability
1\eqn{-}\code{alpha}, that utilizes the uncertain prior
information that
\eqn{\tau = } 0 belongs to a class of confidence
intervals indexed
by the functions \eqn{b} and \eqn{s}.
The function \eqn{b} is an odd continuous function and
the function \eqn{s} is an even
continuous function. In addition, \eqn{b(x)=0} and
\eqn{s(x)} is equal to the
1\eqn{-}\code{alpha}\eqn{/2}
quantile of the
standard normal distribution for all \eqn{|x| \ge d}, where
\eqn{d} is a given positive number.
Extensive numerical explorations
have been used to find a formula (in terms of
\code{alpha} and \code{rho}) for a 'goldilocks'
value of \eqn{d} that is neither too large nor too small.
Then let \eqn{q}=ceiling(\eqn{d}/0.75) and \eqn{h=d/q}.
The values of the functions \eqn{b} and \eqn{s} in
the interval \eqn{[-d,d]}
are specified by the \eqn{(2q-1)}-vector

\eqn{\big(b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h) \big)}.

The values of \eqn{b(kh)} and \eqn{s(kh)} for \eqn{k=-q,...,q} are
deduced from this vector using the assumptions made about
the functions \eqn{b} and \eqn{s}.
The values of \eqn{b(x)} and \eqn{s(x)} for any \eqn{x} in the interval
\eqn{[-d, d]}
are then found using cube spline interpolation using the
values of \eqn{b(kh)} and \eqn{s(kh)} for \eqn{k=-q,...,q}.
For \code{natural}=1 (default) this is 'natural' cubic
spline interpolation and for \code{natural}=0 this is
'clamped' cubic spline interpolation.

The vector \eqn{\big(b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h)\big)}
is found by numerical nonlinear constrained optimization
so that the confidence interval has minimum
coverage probability 1\eqn{-}\code{alpha} and utilizes
the uncertain prior information
through its desirable expected length properties.
This optimization is performed using the
\code{slsqp} function
in the \code{nloptr} package.
}
\examples{
\donttest{alpha <- 0.05
rho <- - 1 / sqrt(2)
bs.list <- bs_ciuupi(alpha, rho)}


}
