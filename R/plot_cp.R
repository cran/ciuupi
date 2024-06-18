#' Plot the graph of the coverage probability
#' of the CIUUPI
#'
#' The input \code{bs.list} determines the functions \eqn{b} and \eqn{s}
#' that specify the confidence interval that utilizes the uncertain
#' prior information (CIUUPI), for all possible values of \eqn{\sigma}
#' and observed response vector.
#' The coverage probability of the CIUUPI is an even function of
#' the unknown parameter
#' \eqn{\gamma = \tau \big/ \big(\sigma \, v_{\tau}^{1/2} \big)}.
#' The R function \code{plot_cp} plots the graph of the coverage probability
#' of the CIUUPI, as a function of \eqn{|\gamma|}.
#' To provide a stringent
#' assessment of this coverage probability, we use a fine equally-spaced grid
#' \code{seq(0, (d+4), by = 0.01)} of values of
#' \eqn{\gamma} and Gauss Legendre quadrature
#' using 10 nodes in the relevant integrals. By contrast,
#' for the computation
#' of the CIUUPI, implemented in \code{bs_ciuupi}, we
#' require that the coverage probability of this confidence
#' interval is greater than or equal to \eqn{1-\alpha}
#' for the equally-spaced grid \code{seq(0, (d+2), by = 0.05)} of values of
#' \eqn{\gamma} and we use Gauss Legendre quadrature
#' with 5 nodes in the relevant integrals.
#'
#' @param bs.list A list that includes the following
#' components.
#'
#' \code{alpha}: the desired minimum coverage is \eqn{1 - \alpha}.
#'
#' \code{rho}: the known correlation
#' between \eqn{\widehat{\theta}}
#' and \eqn{\widehat{\tau}}. This correlation is computed from the
#' \eqn{p}-vectors
#' \eqn{a} and \eqn{c} and the \eqn{n \times p}
#' design matrix \eqn{X} using the formula
#' \eqn{\rho=a^{\top}(X^{\top}X)^{-1}c /(v_{\theta} \, v_{\tau})^{1/2}}, where
#' \eqn{v_{\theta}
#' =a^{\top}(X^{\top}X)^{-1} a} and
#' \eqn{v_{\tau}
#' =c^{\top}(X^{\top}X)^{-1} c}.
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
#'         that specify the functions b and s. In the description
#'         of \code{bsvec}, \code{n.ints} is also called \eqn{q}.
#'
#' \code{bsvec}: the \eqn{(2q-1)}-vector
#'       \deqn{\big(b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h) \big),}
#'       where \eqn{q}=ceiling(\eqn{d}/0.75) and \eqn{h=d/q}.
#'
#'
#'
#' @return
#' A plot of the graph of the coverage probability of the
#' CIUUPI as a function of \eqn{|\gamma|}, where
#' \eqn{\gamma} denotes the unknown parameter
#' \eqn{\tau \big/ \big(\sigma \, v_{\tau}^{1/2} \big)}.
#'
#' @export
#'
#' @examples
#' plot_cp(bs.list.example)
#'
plot_cp <- function(bs.list){
  # For the functions b and s specified by the input bs.list
  # this module plots the graph of the coverage probability
  # (cp) of the confidence interval CI(b,s).
  #
  # Input
  # bs.list which is a list that includes the following
  # components.
  #
  # alpha: the desired minimum coverage is 1 - alpha
  # rho: the correlation between the least squares estimators of theta
  # and tau
  # natural: 1 when the functions b and s are
  #          specified by natural cubic spline interpolation
  #          or 0 if these functions are specified by clamped
  #          cubic spline interpolation
  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  # bsvec: the vector (b(d/n.ints),...,b((n.ints-1)d/n.ints),
  #                s(0),s(d/n.ints)...,s((n.ints-1)d/n.ints))
  #    For d=6 and n.ints=6, this vector is
  #    (b(1),...,b(5),s(0),...,s(5)).
  # comp.time: computation time in seconds
  #
  # Output
  # Graph of the coverage probability (cp)
  # of the confidence interval CI(b,s)
  #
  # Written by P. Kabaila in June 2024

  # Reset old par parameters when the current function
  # exits (either naturally or as a result of an error)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  alpha <- bs.list$alpha
  rho <- bs.list$rho
  natural <- bs.list$natural
  d <- bs.list$d
  n.ints <- bs.list$n.ints
  bsvec <- bs.list$bsvec

  # Compute the coverage probability for
  # a grid of values of gamma
  gams <- seq(0, (d+4), by = 0.01)

  n.nodes <- 10

  cp <- cpciuupi(gams, n.nodes, bs.list)

  # Plot the coverage probability
  graphics::par(col.main="blue", cex.main=0.9)
  plot(gams, cp - (1 - alpha), type = "l", lwd = 1, ylab = "",
       las = 1, xaxs = "i", main =
         paste("cov. prob. - (1 - alpha) for alpha=", alpha,
        ", rho=", signif(rho, digits=5), ", natural=", natural),
       col = "blue", xlab = expression(paste("|", gamma, "|")))
  graphics::abline(h = 0, lty = 2)
  min.cp <- min(cp)
  graphics::mtext(paste("Desired min CP=", 1 - alpha, ". Actual min CP =",
              signif(min.cp, digits=10)), cex = 0.9)

}

