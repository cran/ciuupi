#' Compute the scaled expected length of the CIUUPI
#'
#' Evaluate the scaled expected length of the confidence interval that
#' utilizes uncertain prior information (CIUUPI) at \code{gam}.
#' This scaled expected length is defined to be the expected length of the
#' CIUUPI divided
#' by the expected length of the standard \eqn{1 - \alpha} confidence interval.
#' The input \code{bs.list} determines the functions \eqn{b} and \eqn{s}
#' that specify the confidence interval that utilizes the uncertain
#' prior information (CIUUPI), for all possible values of \eqn{\sigma}
#' and observed response vector.
#'
#' @param gam A value of \eqn{\gamma} or vector of \eqn{\gamma} values at which
#' the coverage probability function is evaluated
#'
#' @param n.nodes The number of nodes for the Gauss Legendre quadrature used for
#' the evaluation of the scaled expected length
#'
#' @param bs.list A list that includes the following
#' components.
#'
#' \code{alpha}: \eqn{1 - \alpha} is the desired minimum coverage probability of the
#' confidence interval
#'
#' \code{rho}: The known correlation \eqn{\rho} between
#' \eqn{\widehat{\theta}} and \eqn{\widehat{\tau}}
#
#'
#' \code{natural}: 1 when the functions \eqn{b} and \eqn{s} are
#'          specified by natural cubic spline interpolation
#'          or 0 if these functions are specified by clamped
#'          cubic spline interpolation
#'
#' \code{d}: the functions \eqn{b} and \eqn{s} are specified by
#'    cubic splines on the interval \eqn{[-d, d]}
#'
#' \code{n.ints}: number of equal-length intervals in \eqn{[0, d]}, where
#'         the endpoints of these intervals specify the knots,
#'         belonging to \eqn{[0, d]}, of the cubic spline interpolations
#'         that specify the functions \eqn{b} and \eqn{s}. In the description
#'         of \code{bsvec}, \code{n.ints} is also called \eqn{q}.
#'
#' \code{bsvec}: the \eqn{(2q-1)}-vector
#'       \deqn{\big(b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h)\big),}
#'       where \eqn{q}=ceiling(\eqn{d}/0.75) and \eqn{h=d/q}.
#'
#' @return The value(s) of the scaled expected length at \code{gam}.
#'
#'
#' @export
#'
#' @examples
#' \donttest{gam <- seq(0, 10, by = 0.2)
#' n.nodes <- 10
#' sel <- selciuupi(gam, n.nodes, bs.list.example)}
#'
selciuupi <- function(gam, n.nodes, bs.list){
  # Use this program to compute the scaled expected length of the
  # CIUUPI.
  #
  # Written by R Mainzer, September 2017
  # Revised by P. Kabaila in June 2024

  # Find the nodes and weights of the Gauss Legendre quadrature
  quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")

  alpha <- bs.list$alpha
  rho <- bs.list$rho
  natural <- bs.list$natural
  d <- bs.list$d
  n.ints <- bs.list$n.ints
  bsvec <- bs.list$bsvec

  # Specify the function s
  c.alpha <- stats::qnorm(1 - alpha/2)
  s.spl <- spline_s(bsvec, d, n.ints, c.alpha, natural)

  # Compute the scaled expected length
  res <- rep(0, length(gam))
  for(i in 1:length(gam)){
    res[i] <- compute_sel_v1(gam[i], bsvec, d, n.ints, quad.info, c.alpha, s.spl)
  }

  res

}
