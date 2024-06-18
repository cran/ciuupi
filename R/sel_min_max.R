sel_min_max <- function(y, d, n.ints, quad.info, c.alpha, natural)
  {
  # This module computes the minimum and the maximum of the
  # scaled expected length of the confidence interval
  # CI(b,s) for given function s.
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
  # quad.info: list of Gauss Legendre nodes and weights
  # c.alpha: 1 - alpha/2 quantile of the standard normal distribution
  # natural:
  #
  # Output
  # A list with two components: min.sel and max.sel
  #
  # Written by P. Kabaila in January 2023

  # s.spl: R function that specifies the function s in the interval
  #        [-d,d] as a cubic spline
  s.spl <- spline_s(y, d, n.ints, c.alpha, natural)

  sel.max <- stats::optimize(sel, c(0, d), maximum = TRUE,
                      s.spl = s.spl, d = d, n.ints = n.ints,
                      quad.info = quad.info, c.alpha = c.alpha)$objective
  gam <- 0
  sel.min <- sel(gam, s.spl, d, n.ints, quad.info, c.alpha)

 list(sel.min=sel.min, sel.max=sel.max)

}
