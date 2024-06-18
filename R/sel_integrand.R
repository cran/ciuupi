sel_integrand <- function(x, gam, s.spl, d, n.ints, c.alpha){
  # This function evaluates
  #   (s(x) - c_alpha) * phi(x-gamma)
  # for a vector x. This is the integrand in the expression
  # for the scaled expected length of the confidence interval
  # CI(b,s).
  #
  # Inputs
  # x: vector of values at which SEL integrand is to be evaluated
  # gam: parameter
  # rho: correlation between the least squares estimators of
  #      theta and tau
  # s.spl: R function that specifies the function s in the interval
  #        [-d,d] as a cubic spline
  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  # c.alpha: 1 - alpha/2 quantile of the standard normal distribution
  #
  # Output
  # A vector of values of the function with the same dimension
  # as x.
  #
  # Written by P. Kabaila in January 2023

  tmp1 <- s.spl(x) - c.alpha
  tmp2 <- stats::dnorm(x - gam, 0, 1)

  tmp1 * tmp2

}
