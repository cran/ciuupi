objective_v1 <- function(y, lambda, d, n.ints, quad.info, c.alpha, natural){
  # This module evaluates the objective function,
  # for given functions b and s, for the numerical
  # nonlinear constrained optimization.
  # In other words, this module computes
  #
  # int_0^d (s(x) - c_alpha)(lambda + phi(x)) dx
  #
  # Inputs
  # y: the vector (b(d/n.ints),...,b((n.ints-1)d/n.ints),
  #                s(0),s(d/n.ints)...,s((n.ints-1)d/n.ints))
  #    For d=6 and n.ints=6, this vector is
  #    (b(1),...,b(5),s(0),...,s(5)).
  # lambda: tuning constant for the objective function
  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  # quad.info: list of Gauss Legendre nodes and weights
  # c.alpha: 1 - alpha/2 quantile of the standard normal distribution
  # natural: 1 (default) when the functions b and s are
  #          specified by natural cubic spline interpolation
  #          or 0 if these functions are specified by clamped
  #          cubic spline interpolation.
  #
  # Output:
  # The objective function.
  #
  # Written by P.Kabaila in June 2008.
  # Rewritten in R by R Mainzer, March 2017
  # Revised by P.Kabaila in January 2023
  # Changes made by P. Kabaila in January 2023
  # are highlighted in yellow.

  # c.alpha <- stats::qnorm(1 - alpha/2)

  s.spl <- spline_s(y, d, n.ints, c.alpha, natural)

  # Specify where the knots of the cubic splines are located
  knots <- seq(0, d, by = d/n.ints)

  # Set up a vector to store the n.ints integral evaluations
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
    q <- integrand_obj(adj.nodes, y, lambda, d, n.ints, c.alpha, s.spl)
    int[i] <- ((b - a) / 2) * sum(weights * q)
  }

  sum(int)

}
