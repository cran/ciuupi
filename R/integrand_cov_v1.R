integrand_cov_v1 <- function(x, gam, rho, y, d, n.ints,
                             c.alpha, b.spl, s.spl){
  # This function evaluates
  # (k(x, gam, rho) - k_dag(x, gam, rho)) * phi(x - gam)
  # for a vector x. This is the integrand in the expression
  # for the coverage probability.
  #
  # Inputs
  # x: vector at which CP integrand is to be evaluated
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
  # b.spl: R function that specifies the function b in the interval
  #        [-d,d] as a cubic spline
  # s.spl: R function that specifies the function s in the interval
  #        [-d,d] as a cubic spline
  #
  # Output
  # A vector of values of the integrand of the CP with the same
  # dimension as the vector x.
  #
  # Written by P.Kabaila in June 2008
  # Rewritten in R by R Mainzer in March 2017
  # Revised by P. Kabaila in January 2023
  # Changes made by P. Kabaila in January 2023
  # are highlighted in yellow.

#  c.alpha <- stats::qnorm(1 - alpha/2)

  mu1 <- rho * (x - gam)
  var <- 1 - rho^2
  k.dag1 <- Psi(-c.alpha, c.alpha, mu1, var)

  term.a1 <- b.spl(x)
  term.b1 <- s.spl(x)

  lh <- term.a1 - term.b1
  uh <- term.a1 + term.b1

  k1 <- Psi(lh, uh, mu1, var)
  term1 <- stats::dnorm(x - gam, 0, 1)

  mu2 <- rho * (-x - gam)
  k.dag2 <- Psi(-c.alpha, c.alpha, mu2, var)

  term.a2 <- b.spl(-x)
  term.b2 <- s.spl(-x)

  lh2 <- term.a2 - term.b2
  uh2 <- term.a2 + term.b2

  k2 <- Psi(lh2, uh2, mu2, var)
  term2 <- stats::dnorm(x + gam, 0, 1)

  (k1 - k.dag1) * term1 + (k2 - k.dag2) * term2

}
