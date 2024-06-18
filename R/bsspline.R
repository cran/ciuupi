#' Evaluate the functions \eqn{b} and \eqn{s} at \code{x}
#'
#' Evaluate the functions \eqn{b} and \eqn{s}, as specified by
#' \code{bsvec},
#' \code{alpha}, \code{d}, \code{n.ints} and \code{natural}, at \code{x}.
#'
#' @param x A value or vector of values at which the functions \eqn{b}
#' and \eqn{s} are to be evaluated
#'
#' @param bsvec  The \eqn{(2q-1)}-vector
#'       \deqn{\big(b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h) \big),}
#'       where \eqn{q}=ceiling(\eqn{d}/0.75) and \eqn{h=d/q}.
#'       This vector specifies the CIUUPI, for all possible values of the random
#'       error variance and the observed response vector
#'
#' @param alpha The desired minimum coverage probability is \eqn{1 - \alpha}
#'
#' @param d The functions \eqn{b} and \eqn{s} are specified by cubic splines
#' on the interval \eqn{[-d, d]}
#'
#' @param n.ints The number of equal-length intervals in \eqn{[0, d]}, where
#'         the endpoints of these intervals specify the knots,
#'         belonging to \eqn{[0, d]}, of the cubic spline interpolations
#'         that specify the functions \eqn{b} and \eqn{s}. In the description
#'         of \code{bsvec}, \code{n.ints} is also called \eqn{q}.
#'
#' @param natural Equal to 1 (default) if the \eqn{b} and \eqn{s} functions are obtained by
#' natural cubic spline interpolation or 0 if obtained by clamped cubic spline
#' interpolation
#'
#' @return A data frame containing \code{x} and the corresponding values of the
#' functions \eqn{b} and \eqn{s}.
#'
#' @export
#'
#' @examples
#' x <- seq(0, 8, by = 1)
#' alpha <- bs.list.example$alpha
#' natural <- bs.list.example$natural
#' d <- bs.list.example$d
#' n.ints <- bs.list.example$n.ints
#' bsvec <- bs.list.example$bsvec
#' bs <- bsspline(x, bsvec, alpha, d, n.ints, natural)
#'
bsspline <- function(x, bsvec, alpha, d, n.ints, natural){
  # This module computes the values of the functions b
  # and s for the vector x of values, whose elements are
  # assumed to be in increasing order.
  #
  # Inputs
  # x: a vector of values, whose elements are assumed
  #    to be in increasing order
  # bsvec: the vector (b(d/n.ints),...,b((n.ints-1)d/n.ints),
  #                s(0),s(d/n.ints)...,s((n.ints-1)d/n.ints))
  #    For d=6 and n.ints=6, this vector is
  #    (b(1),...,b(5),s(0),...,s(5)).
  # alpha: the desired minimum coverage is 1 - alpha
  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  # natural: 1 when the functions b and s are
  #          specified by natural cubic spline interpolation
  #          or 0 if these functions are specified by clamped
  #          cubic spline interpolation.
  # Output
  # A data frame with x, values of b and values of s
  #
  # Written by P. Kabaila in January 2023, based
  # on earlier R code of R. Mainzer from 2017

  c.alpha <- stats::qnorm(1 - alpha/2)

  # Find b and s functions
  sspl <- spline_s(bsvec, d, n.ints, c.alpha, natural)
  bspl <- spline_b(bsvec, d, n.ints, c.alpha, natural)

  x1 <- x[which(x <= -d)]
  x2 <- x[which(x > -d & x < d)]
  x3 <- x[which(x >= d)]

  bspl.res <- c(rep(0, length(x1)), bspl(x2), rep(0, length(x3)))
  sspl.res <- c(rep(c.alpha, length(x1)), sspl(x2), rep(c.alpha, length(x3)))

  data.frame(x = x, b = bspl.res, s = sspl.res)

}
