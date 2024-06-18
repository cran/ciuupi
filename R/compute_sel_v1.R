compute_sel_v1 <- function(gam, y, d, n.ints, quad.info, c.alpha, s.spl){
  # This function computes the value of the scaled expected
  # length of the confidence interval CI(b,s) for given
  # functions b and s.
  # In other words, this function computes
  #
  # 1 + (1/c_alpha) int_{-d}^d (s(x) - c_alpha) phi(x-gamma) dx
  #
  # Inputs
  # gam: parameter
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
  # s.spl: R function that specifies the function s in the interval
  #        [-d,d] as a cubic spline
  #
  # Output:
  # The scaled expected length for given functions b and s.
  #
  # Written by P.Kabaila in June 2008.
  # Rewritten in R by R Mainzer, March 2017
  # Revised by P. Kabaila in January 2023

  # c.alpha <- stats::qnorm(1 - alpha/2)

  # Specify where the knots are located
  knots <- seq(-d, d, by = d/n.ints)

  #  Set up a vector to store the n.ints integral evaluations
  int <- rep(0, length(knots))

  # Find the nodes and weights of the Gauss Legendre quadrature
  # quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
  nodes <- quad.info$nodes
  weights <- quad.info$weights

  for(i in 1:(length(knots) - 1)){
    # Specify bounds of the integral
    a <- knots[i]
    b <- knots[i+1]

    # Find the approximate integral
    adj.nodes <- ((b - a) / 2) * nodes + (a + b) / 2
    q <- integrand_sel_v1(adj.nodes, gam, y, d, n.ints, c.alpha, s.spl)
    int[i] <- ((b - a) / 2) * sum(weights * q)
  }

  1 + (sum(int) / c.alpha)

}
