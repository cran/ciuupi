integrand_sel_v1 <- function(x, gam, y, d, n.ints, c.alpha, s.spl){
  # This function evaluates
  #   (s(x) - c_alpha) * phi(x-gamma)
  # for a vector x. This is the integrand in the expression
  # for the scaled expected length
  #
  # Inputs
  # x: vector of values at which SEL integrand is to be evaluated
  # gam: parameter
  # rho: correlation between the least squares estimators of
  #      theta and tau
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
  # s.spl: R function that specifies the function s in the interval
  #        [-d,d] as a cubic spline
  #
  # Output:
  # A vector of values of the function with the same dimension
  # as x.
  #
  # Written by P.Kabaila in June 2008.
  # Rewritten in R by R Mainzer, March 2017
  # Revised by P. Kabaila in January 2023

  tmp1 <- s.spl(x) - c.alpha
  tmp2 <- stats::dnorm(x - gam, 0, 1)
  tmp1 * tmp2

}
