start_standard_ci <- function(d, n.ints, alpha){
  # This module calculates the vector of values
  # of b and s functions evaluated at the knots
  # for the standard 1-alpha confidence interval
  # for theta.
  # These are used as starting values for the
  # numerical nonlinear constrained optimization.
  #
  # The main use of this function is to provide
  # a starting value for the optimization problem.
  #
  # Inputs
  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  # alpha: the desired minimum coverage probability is 1-alpha
  #
  # Output
  # For the standard 1-alpha confidence interval for theta,
  # the vector (b(d/n.ints),...,b((n.ints-1)d/n.ints),
  #                s(0),s(d/n.ints)...,s((n.ints-1)d/n.ints))
  #    For d=6 and n.ints=6, this vector is
  #    (b(1),...,b(5),s(0),...,s(5)).
  #
  # Written by P. Kabaila in January 2023.

  c.alpha <- stats::qnorm(1 - alpha/2)

  c(rep(0, n.ints - 1), rep(c.alpha, n.ints))

}
