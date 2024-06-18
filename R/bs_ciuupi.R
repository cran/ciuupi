#' Computes the the functions \eqn{b} and \eqn{s}
#' that specify the CIUUPI for all possible
#' values of \eqn{\sigma} and the observed response vector
#'
#' Chooses the positive number \eqn{d} and the positive integer \eqn{q}, sets
#' \eqn{h=d/q}, and then computes the
#' \eqn{(2q-1)}-vector
#' \eqn{\big(b(h),...,b((q-1)h),
#' s(0),s(h)...,s((q-1)h)\big)}
#' that determines, via cubic spline interpolation, the functions
#' \eqn{b} and \eqn{s} which specify
#' the confidence interval for \eqn{\theta}
#' that utilizes the uncertain prior information (CIUUPI),
#' for all possible values of \eqn{\sigma} and the observed response vector.
#' To an excellent approximation, this confidence interval
#'  has minimum coverage probability
#' \eqn{1-\alpha}.
#'
#'
#' @param alpha    The desired minimum coverage probability
#' is \eqn{1-\alpha}
#'
#' @param rho      The known correlation \eqn{\rho} between
#' \eqn{\widehat{\theta}} and \eqn{\widehat{\tau}}
#'
#' @param natural  Equal to 1 (default) if the functions \eqn{b}
#' and \eqn{s} are specified by natural cubic spline interpolation
#' or 0 if these functions are specified by clamped cubic
#' spline interpolation in an interval \eqn{[-d, d]}, where \eqn{d}
#' is computed by \code{bs_ciuupi} using a specified
#' function of
#' \code{alpha} and \code{rho}
#'
#' @details
#' Suppose that \deqn{y = X \beta + \varepsilon} where \eqn{y}
#' is a random \eqn{n}-vector of
#' responses, \eqn{X} is a known \eqn{n} by \eqn{p} matrix
#' with linearly
#' independent columns, \eqn{\beta} is an unknown parameter
#'  \eqn{p}-vector and
#' \eqn{\varepsilon} is the random error with components that are iid normally distributed
#' with zero mean and known variance \eqn{\sigma^2}.
#' The parameter of interest is
#' \eqn{\theta = a^{\top} \beta}.
#' Also let \eqn{\tau = c^{\top}\beta -t}, where \eqn{a}
#' and \eqn{c} are specified linearly independent
#' vectors and \eqn{t} is a specified number.
#' The uncertain prior information is that \eqn{\tau = 0}.
#'
#'  Let \code{rho} denote the known
#' correlation between the \eqn{\widehat{\theta}} and \eqn{\widehat{\tau}}.
#' We can compute \code{rho}
#' from given values of \eqn{a}, \eqn{c} and \eqn{X}
#' using the function \code{acX_to_rho}.
#'
#' The confidence interval for \eqn{\theta},
#' with minimum coverage probability
#'  1\eqn{-}\code{alpha}, that utilizes the uncertain prior
#'  information that
#'  \eqn{\tau = } 0 belongs to a class of confidence
#'  intervals indexed
#'  by the functions \eqn{b} and \eqn{s}.
#' The function \eqn{b} is an odd continuous function and
#' the function \eqn{s} is an even
#' continuous function. In addition, \eqn{b(x)=0} and
#' \eqn{s(x)} is equal to the
#' 1\eqn{-}\code{alpha}\eqn{/2}
#' quantile of the
#' standard normal distribution for all \eqn{|x| \ge d}, where
#' \eqn{d} is a given positive number.
#' Extensive numerical explorations
#' have been used to find a formula (in terms of
#' \code{alpha} and \code{rho}) for a 'goldilocks'
#' value of \eqn{d} that is neither too large nor too small.
#' Then let \eqn{q}=ceiling(\eqn{d}/0.75) and \eqn{h=d/q}.
#' The values of the functions \eqn{b} and \eqn{s} in
#' the interval \eqn{[-d,d]}
#' are specified by the \eqn{(2q-1)}-vector
#'
#' \eqn{\big(b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h) \big)}.
#'
#' The values of \eqn{b(kh)} and \eqn{s(kh)} for \eqn{k=-q,...,q} are
#' deduced from this vector using the assumptions made about
#' the functions \eqn{b} and \eqn{s}.
#' The values of \eqn{b(x)} and \eqn{s(x)} for any \eqn{x} in the interval
#' \eqn{[-d, d]}
#' are then found using cube spline interpolation using the
#' values of \eqn{b(kh)} and \eqn{s(kh)} for \eqn{k=-q,...,q}.
#' For \code{natural}=1 (default) this is 'natural' cubic
#' spline interpolation and for \code{natural}=0 this is
#' 'clamped' cubic spline interpolation.
#'
#' The vector \eqn{\big(b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h)\big)}
#' is found by numerical nonlinear constrained optimization
#' so that the confidence interval has minimum
#' coverage probability 1\eqn{-}\code{alpha} and utilizes
#' the uncertain prior information
#' through its desirable expected length properties.
#' This optimization is performed using the
#' \code{slsqp} function
#' in the \code{nloptr} package.
#'
#' @return A list with the following components.
#'
#' alpha, rho, natural: the inputs
#'
#' d: a 'goldilocks' value of \eqn{d} that is not too large
#'               and not too small
#'
#' n.ints: number of equal-length consecutive
#'         intervals whose union is \eqn{[0,d]},
#'         this is the same as \eqn{q}
#'
#' lambda.star: the computed value of \eqn{\lambda^*}
#'
#' bsvec: the vector
#' \eqn{\big(b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h) \big)} that determines
#' the functions \eqn{b} and \eqn{s} that specify the CIUUPI for all possible
#' values of \eqn{\sigma} and observed response vector
#'
#' comp.time: the computation time in seconds
#'
#'
#' @export
#'
#' @examples
#' \donttest{alpha <- 0.05
#' rho <- - 1 / sqrt(2)
#' bs.list <- bs_ciuupi(alpha, rho)}
#'
#'
bs_ciuupi <- function(alpha, rho, natural=1){
  # This module computes the value of the vector
  # bsvec = (b(d/n.ints),...,b((n.ints-1)d/n.ints),
  #      s(0),s(d/n.ints)...,s((n.ints-1)d/n.ints))
  #      For d=6 and n.ints=6, this vector is
  #         (b(1),...,b(5),s(0),...,s(5)).
  # that specifies the CIUUPI. This vector is found
  # by numerical nonlinear constrained optimization.
  #
  # Inputs
  # alpha: the desired minimum coverage probability is 1-alpha
  # rho: correlation between the least squares estimators of
  #      theta and tau
  # natural: 1 (default) when the functions b and s are
  #          specified by natural cubic spline interpolation
  #          or 0 if these functions are specified by clamped
  #          cubic spline interpolation
  #
  # Output
  # A list with the following components.
  # alpha, rho, natural: the inputs
  # d: the value of d that is not too large
  #               and not too small
  # n.ints: number of equal-length consecutive
  #         intervals whose union is [0,d]
  # lambda.star:
  # bsvec = (b(d/n.ints),...,b((n.ints-1)d/n.ints),
  #      s(0),s(d/n.ints)...,s((n.ints-1)d/n.ints))
  # comp.time: the computation time in seconds
  #
  # Written by P. Kabaila in May 2024

  # n.nodes: the number of nodes for the Gauss Legendre quadrature
  #          used for the evaluation of the coverage probability
  #          and the objective function
  n.nodes <- 5

  # nl.info: if TRUE then the following additional information
  #          is printed to the console
  #          Number of iterations
  #          Optimal value of objective function
  #          If, on the other hand, nl.info is FALSE then this
  #          additional information is not printed to the console
  nl.info <- FALSE

  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  d <- d_goldilocks(alpha, rho)

  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  n.ints <- ceiling(d / 0.75)

  named.comp.times1 <- system.time(
    temp <-
      eta_star(rho, d, n.ints, natural, alpha, n.nodes, nl.info)
  )
  named.user.time1 <- named.comp.times1[1]
  user.time1 <- unname(named.user.time1)

  eta.star <- temp$eta.star
  start <- temp$bsvec.matrix[,3]

  # gams: vector of values of the parameter gam at
  #       which the coverage is required to be greater
  #       than or equal to 1 - alpha
  gams <- seq(0, (d+2), by = 0.05)
  named.comp.times2 <- system.time(
  bsvec.g.list <- g_fn(eta.star, rho, alpha, gams,
                  d, n.ints, n.nodes, natural, start, nl.info)
  )
  named.user.time2 <- named.comp.times2[1]
  user.time2 <- unname(named.user.time2)

  bsvec <- bsvec.g.list$bsvec

  comp.time <- user.time1 + user.time2

  lambda.star <- exp(eta.star)

  list(alpha=alpha, rho=rho, natural=natural, d=d, n.ints=n.ints,
       lambda.star=lambda.star, bsvec=bsvec, comp.time=comp.time)

}
