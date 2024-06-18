sel <- function(gam, s.spl, d, n.ints, quad.info, c.alpha){
  # This function computes the value of the scaled expected
  # length of the confidence interval CI(b,s) for given
  # functions b and s.
  # In other words, this function computes
  #
  # 1 + (1/c_alpha) int_{-d}^d (s(x) - c_alpha) phi(x-gamma) dx
  #
  # Inputs
  # gam: parameter
  # s.spl: R function that specifies the function s in the interval
  #        [-d,d] as a cubic spline
  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  # quad.info: list of Gauss Legendre nodes and weights
  # c.alpha: 1 - alpha/2 quantile of the standard normal distribution
  #
  # Output:
  # The scaled expected length for given functions b and s.
  #
  # Written by P. Kabaila in January 2023

  # Specify where the knots are located
  knots <- seq(-d, d, by = d/n.ints)

  #  Set up a vector to store the integral evaluations
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
    q <- sel_integrand(adj.nodes, gam, s.spl, d, n.ints, c.alpha)
    int[i] <- ((b - a) / 2) * sum(weights * q)
  }

  1 + (sum(int) / c.alpha)

}
