compute_cov_legendre_v1 <- function(gam, rho, y, d, n.ints, alpha,
                                 quad.info, b.spl, s.spl){
  # This function computes the coverage probability of CI(b, s).
  # The integral from [0, d] is expressed as the sum of n.int
  # integrals whose endpoints are the successive knots in [0,d].
  # Each of these n.int integrals is computed using Gauss Legendre
  # quadrature.
  #
  # Inputs
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
  # alpha: the desired minimum coverage is 1 - alpha
  # quad.info: list of Gauss Legendre nodes and weights
  # b.spl: R function that specifies the function b in the interval
  #        [-d,d] as a cubic spline
  # s.spl: R function that specifies the function s in the interval
  #        [-d,d] as a cubic spline
  #
  # Written by R Mainzer, May 2017
  # Revised by P. Kabaila in January 2023
  # Changes made by P. Kabaila in January 2023
  # are highlighted in yellow.

  c.alpha <- stats::qnorm(1 - alpha/2)

  # Specify where the knots of the cubic splines are located
  knots <- seq(0, d, by = d/n.ints)

  #  Set up a vector to store the n.ints integral evaluations
  int <- rep(0, n.ints)

  # Find the nodes and weights of the Gauss Legendre quadrature
  # quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
  nodes <- quad.info$nodes
  weights <- quad.info$weights

  for(i in c(1:n.ints)){
    # Specify bounds of the integral
    a <- knots[i]
    b <- knots[i+1]

    # Find the approximate integral
    adj.nodes <- ((b - a) / 2) * nodes + (a + b) / 2
    # On 6 Apr 2024 try integrand_cov_v2
    q <- integrand_cov_v1(adj.nodes, gam, rho, y, d, n.ints,
                          c.alpha, b.spl, s.spl)
    int[i] <- ((b - a) / 2) * sum(weights * q)
  }

  (1 - alpha) + sum(int)

}
