#' Plot the graph of the squared scaled expected length
#' of the CIUUPI
#'
#' The input \code{bs.list} determines the functions \eqn{b} and \eqn{s}
#' that specify the confidence interval that utilizes the uncertain
#' prior information (CIUUPI), for all possible values of \eqn{\sigma}
#' and observed response vector.
#' The scaled expected length of the CIUUPI is an even function of
#' the unknown parameter
#' \eqn{\gamma = \tau \big/ \big(\sigma \, v_{\tau}^{1/2} \big)}.
#' The R function \code{plot_squared_sel} plots the graph of the squared scaled
#' expected length (i.e. squared SEL)
#' of the CIUUPI, as a function of \eqn{|\gamma|}.
#' To provide a stringent
#' assessment of this squared SEL, we use a grid
#' \code{seq(0, (d+4), by = 0.01)} of values of
#' \eqn{\gamma} and Gauss Legendre quadrature
#' with 10 nodes in the relevant integrals. By contrast,
#' for the computation
#' of the CIUUPI, implemented in \code{bs_ciuupi}, we
#' use Gauss Legendre quadrature
#' with 5 nodes in the relevant integrals.
#'
#' @param bs.list A list that specifies the CIUUPI and includes the following
#' components.
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
#' A plot of the graph of the squared scaled
#' expected length (i.e. squared SEL)
#' of the CIUUPI as a function of \eqn{|\gamma|}, where
#' \eqn{\gamma} denotes the unknown parameter
#' \eqn{\tau \big/ \big(\sigma \, v_{\tau}^{1/2} \big)}.
#'
#' @export
#'
#' @examples
#' plot_squared_sel(bs.list.example)
#'
plot_squared_sel <- function(bs.list){
  # For the functions b and s specified by the input bs.list
  # this module plots the graph of the squared SEL i.e. the
  # the square of the scaled expected length
  # of the confidence interval CI(b,s).
  #
  # Input
  # bs.list which is a list that specifies the CIUUPI and
  # includes the following components.
  #
  # alpha: the desired minimum coverage is 1 - alpha
  #
  # natural: 1 when the functions b and s are
  #          specified by natural cubic spline interpolation
  #          or 0 if these functions are specified by clamped
  #          cubic spline interpolation
  #
  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  #
  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  #
  # bsvec: the vector (b(d/n.ints),...,b((n.ints-1)d/n.ints),
  #                s(0),s(d/n.ints)...,s((n.ints-1)d/n.ints))
  #    For d=6 and n.ints=6, this vector is
  #    (b(1),...,b(5),s(0),...,s(5)).
  #
  # Output
  # A plot of the graph of the squared SEL of the CIUUPI
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

  # Compute the scaled expected length for
  # a grid of values of gamma
  gams <- seq(0, (d+4), by = 0.01)

  n.nodes <- 10
  sel <- selciuupi(gams, n.nodes, bs.list)

  # Plot the graph of the squared scaled expected length
  graphics::par(col.main="blue", cex.main=0.9)
  plot(gams, sel^2, type = "l", lwd = 1, ylab = "", las = 1, xaxs = "i",
       main = paste("squared SEL for alpha=", alpha,
                    ", rho=", signif(rho, digits=5), ", natural=", natural), col = "blue",
       xlab = expression(paste("|", gamma, "|")))

  # c.alpha is the 1-alpha/2 quantile of the N(0,1) distribution
  c.alpha <- stats::qnorm(1 - alpha/2)

  # Find the nodes and weights of the Gauss Legendre quadrature
  quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")

  temp <- sel_min_max(bsvec, d, n.ints, quad.info, c.alpha,
                      natural)
  sel.min.squared <- temp$sel.min^2
  sel.max.squared <- temp$sel.max^2

  # # Code check
  # sel.min.graph.squared <- min(sel^2)
  # sel.max.graph.squared <- max(sel^2)
  # cat(paste("sel.min.graph.squared - sel.min.squared =",
  #           sel.min.graph.squared - sel.min.squared), "\n")
  # cat(paste("sel.max.graph.squared - sel.max.squared =",
  #          sel.max.graph.squared - sel.max.squared), "\n")

  graphics::mtext(paste("If prior info true, gain=",
              signif(1 - sel.min.squared, digits=6),
              ". Max possible loss=",
              signif(sel.max.squared - 1, digits=6)), cex=0.9)
  graphics::abline(h = 1, lty = 2)

}

