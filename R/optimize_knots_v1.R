optimize_knots_v1 <- function(lambda, rho, alpha, gams,
                            d, n.ints, n.nodes, natural){
  # This module computes the value of the vector
  # y = (b(d/n.ints),...,b((n.ints-1)d/n.ints),
  #      s(0),s(d/n.ints)...,s((n.ints-1)d/n.ints))
  #      For d=6 and n.ints=6, this vector is
  #         (b(1),...,b(5),s(0),...,s(5)).
  # that specifies the CIUUPI, for a given value of lambda.
  # This vector is found by numerical nonlinear constrained
  # optimization.
  #
  # Inputs
  # lambda: tuning constant for the objective function
  # rho: correlation between the least squares estimators of
  #      theta and tau
  # alpha: the desired minimum coverage probability is 1-alpha
  # gams: vector of values of the parameter gam at which the coverage
  #       is required to be greater than or equal to 1 - alpha
  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  # n.nodes: the number of nodes for the Gauss Legendre quadrature
  #          used for the evaluation of the coverage probability
  #          and the objective function
  # natural: 1 when the functions b and s are
  #          specified by natural cubic spline interpolation
  #          or 0 if these functions are specified by clamped
  #          cubic spline interpolation.
  #
  # Output
  # A vector with the values at the knots of the b and s
  # functions
  #
  # Written by P.Kabaila in June 2008.
  # Rewritten in R by R Mainzer, March 2017
  # Revised by P. Kabaila in January 2023
  # Changes made by P. Kabaila in January 2023
  # are highlighted in yellow.

  # Compute the 1-alpha/2 quantile of the standard normal distribution
  c.alpha <- stats::qnorm(1 - alpha/2)

  # Find a starting value for the optimization function
  start <- standard_CI(d, n.ints, alpha)

  # Specify lower and upper bounds on the vector of values
  # of the b and s functions evaluated at the knots
  low <- c(rep(-100, n.ints - 1), rep(0.5, n.ints))
  up <- c(rep(100, n.ints - 1), rep(200, n.ints))

  # Find the nodes and weights of the Gauss Legendre quadrature
  quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")

  # Make the objective function a function of one argument, y
  obj_fun <- functional::Curry(objective_v1, lambda = lambda, d = d,
                               n.ints = n.ints, quad.info = quad.info,
                               c.alpha = c.alpha, natural = natural)

  # Make the constraint function a function of one argument, y
  cons_fun <- functional::Curry(constraints_slsqp_gausslegendre_v1,
                                gams = gams, rho = rho, d = d,
                                n.ints = n.ints, alpha = alpha,
                                quad.info = quad.info, natural = natural)

  # Find the values of the knots using the optimization function
  # When nl.info=TRUE, the output includes the following
  # information:
  #    Number of iterations
  #    Optimal value of objective funcion
  # When nl.info=FALSE, the output does not include this information

  res <- nloptr::slsqp(start, obj_fun, hin = cons_fun, lower = low,
               upper = up, nl.info = TRUE)
  new.par <- res$par

  # Output the vector with knot values which specifies the new
  # confidence interval
  new.par

}

