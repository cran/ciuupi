constraints_slsqp_gausslegendre_v1 <- function(gams, rho, y, d, n.ints,
                                            alpha, quad.info, natural){
  # This module computes the coverage probability
  # inequality constraints.
  #
  # Inputs
  # gams: vector of values of the parameter gam at which the coverage
  #       is required to be greater than or equal to 1 - alpha
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
  # alpha: the desired minimum coverage is 1 - alpha
  # quad.info: list of Gauss Legendre nodes and weights
  # natural: 1 when the functions b and s are
  #          specified by natural cubic spline interpolation
  #          or 0 if these functions are specified by clamped
  #          cubic spline interpolation.
  #
  # Output:
  # A vector of coverage probability inequality constraints
  #
  # Written by P.Kabaila in June 2008
  # Rewritten in R by R Mainzer, March 2017
  # Revised by P. Kabaila in January 2023
  # Changes made by P. Kabaila in January 2023
  # are highlighted in yellow.

  len.gams <- length(gams)
  covs <- rep(0, len.gams)

  c.alpha <- stats::qnorm(1 - alpha/2)

  b.spl <- spline_b(y, d, n.ints, c.alpha, natural)
  s.spl <- spline_s(y, d, n.ints, c.alpha, natural)

  for(i in 1:len.gams){
    covs[i] <- compute_cov_legendre_v1(gams[i], rho, y, d, n.ints,
                                       alpha, quad.info, b.spl, s.spl)
  }

  covs - (1 - alpha)

}
