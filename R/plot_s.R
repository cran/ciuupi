#' Plot the graph of the even function \eqn{s} used in the specification of the CIUUPI
#'
#' The input \code{bs.list} determines the functions \eqn{b} and \eqn{s}
#' that specify the confidence interval that utilizes the uncertain
#' prior information (CIUUPI), for all possible values of \eqn{\sigma}
#' and observed response vector. The R function \code{plot_s}
#' plots the graph of the odd function \eqn{s}.
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
#' @return
#' A plot of the graph of the even function \eqn{s} used in the
#' specification of the CIUUPI.
#'
#' @export
#'
#' @examples
#' plot_s(bs.list.example)
#'
plot_s <- function(bs.list){
    # For the function s specified by the input bs.list
    # plot the graph of the function s.
    #
    # Input
    # bs.list which is a list that includes the following
    # components.
    #
    # bsvec: the vector (b(d/n.ints),...,b((n.ints-1)d/n.ints),
    #                s(0),s(d/n.ints)...,s((n.ints-1)d/n.ints))
    #    For d=6 and n.ints=6, this vector is
    #    (b(1),...,b(5),s(0),...,s(5)).
    # d: the functions b and s are specified by
    #    cubic splines on the interval [-d, d]
    # n.ints: number of equal-length intervals in [0, d], where
    #         the endpoints of these intervals specify the knots,
    #         belonging to [0,d], of the cubic spline interpolations
    #         that specify the functions b and s
    # alpha: the desired minimum coverage is 1 - alpha
    # natural: 1 when the functions b and s are
    #          specified by natural cubic spline interpolation
    #          or 0 if these functions are specified by clamped
    #          cubic spline interpolation
    # comp.time: computation time in seconds
    #
    # Output
    # Graph of the function s
    #
    # Written by P. Kabaila in June 2024

    # Reset old par parameters when the current function
    # exits (either naturally or as a result of an error)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    bsvec <- bs.list$bsvec
    d <- bs.list$d
    n.ints <- bs.list$n.ints
    alpha <- bs.list$alpha
    natural <- bs.list$natural
    rho <- bs.list$rho

    # Compute the values of the functions b and s
    # on x.grid the grid of x-values
    x.grid <- seq(0, (d+2), by = 0.02)
    splineval <-
      bsspline(x.grid, bsvec, alpha, d, n.ints, natural)

    # Compute the values of the function s at the
    # knots
    xseq <- seq(0, d, by = d/n.ints)
    svec <- c(bsvec[n.ints:(2*n.ints-1)], stats::qnorm(1 - alpha/2))

    graphics::par(col.main="blue", cex.main=0.9)
    # Plot the graph of the function s
    plot(x.grid, splineval[, 3], type = "l", xlab = "x",
         ylab = "", las = 1, lwd = 1, xaxs = "i",
         main = paste("s function for alpha=",alpha,", rho=",
                      signif(rho, digits=5),
                      ", natural=", natural),
         col = "blue")
    graphics::points(xseq, svec, pch = 19, col = "blue", cex=0.55)
    c.alpha <- stats::qnorm(1 - alpha/2)
    graphics::abline(h=c.alpha)

}
