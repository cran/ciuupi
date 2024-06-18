spline_s <- function(y, d, n.ints, c.alpha, natural){
  # This module computes the R function that specifies
  # the function s in the interval [-d,d] as a cubic spline.
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
  # R function that specifies the function s in the interval
  # [-d,d] as a cubic spline
  #
  # Written by R Mainzer, March 2017
  # Revised by P. Kabaila in January 2023
  # Changes made by P. Kabaila in January 2023
  # are highlighted in yellow.

  s.knots <- seq(0, d, d/n.ints)
  s.vals <- c(y[n.ints:(2 * n.ints - 1)], c.alpha)

  s.knots.all <- seq(-d, d, d/n.ints)
  s.vals.all <- c(rev(s.vals), s.vals[2:(n.ints+1)])

  if(natural == 1){
    s.spl <- stats::splinefun(s.knots.all, s.vals.all, method = "natural")
  } else {
    s.spl.pp <- pracma::cubicspline(s.knots.all, s.vals.all, endp2nd = TRUE)
    s.spl <- function(x) pracma::ppval(s.spl.pp, x)
  }

  s.spl

}
