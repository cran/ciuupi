spline_b <- function(y, d, n.ints, c.alpha, natural){
  # This module computes the R function that specifies
  # the function b in the interval [-d,d] as a cubic spline.
  #
  # Inputs
  # y: the vector (b(d/n.ints),...,b((n.ints-1)d/n.ints),
  #                s(0),s(d/n.ints)...,s((n.ints-1)d/n.ints))
  #    For d=6 and n.ints=6, this vector is
  #    (b(1),...,b(5),s(0),...,s(5)).
  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  # c.alpha: 1 - alpha/2 quantile of the standard normal distribution
  # natural: 1 when the functions b and s are
  #          specified by natural cubic spline interpolation
  #          or 0 if these functions are specified by clamped
  #          cubic spline interpolation in the interval [-d,d].
  #
  # Output
  # R function that specifies the function b in the interval
  # [-d,d] as a cubic spline
  #
  # Written by R Mainzer, March 2017
  # Revised by P. Kabaila in January 2023
  # Changes made by P. Kabaila in January 2023
  # are highlighted in yellow.
  # What this module does remains unchanged.

  b.knots <- seq(0, d, by = d/n.ints)
  y.rev <- rev(y[1:n.ints - 1])
  b.vals <- c(0, y[1:n.ints - 1], 0)

  b.knots.all <- seq(-d, d, by = d/n.ints)
  b.vals.all <- c(0, -y.rev, b.vals)

  # If natural = 1 use natural cubic spline, otherwise use clamped cubic
  # spline
  if(natural == 1){
    b.spl <- stats::splinefun(b.knots.all, b.vals.all, method = "natural")
  } else {
    b.spl.pp <- pracma::cubicspline(b.knots.all, b.vals.all, endp2nd = TRUE)
    b.spl <- function(x) pracma::ppval(b.spl.pp, x)
  }

  b.spl

}
