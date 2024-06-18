eta_star <- function(rho, d, n.ints, natural,
                     alpha, n.nodes, nl.info){
  # This module computes eta.star = log(lambda.star).
  # The search for eta.star defined to be a tight upper
  # bound on the likely values of eta.
  #
  # Inputs
  # rho: correlation between the least squares estimators of
  #      theta and tau
  # d: the functions b and s are specified by
  #    cubic splines on the interval [-d, d]
  # n.ints: number of equal-length intervals in [0, d], where
  #         the endpoints of these intervals specify the knots,
  #         belonging to [0,d], of the cubic spline interpolations
  #         that specify the functions b and s
  # natural: 1 when the functions b and s are
  #          specified by natural cubic spline interpolation
  #          or 0 if these functions are specified by clamped
  #          cubic spline interpolation in the interval [-d,d].
  # alpha: the desired minimum coverage is 1 - alpha
  # n.nodes: the number of nodes of the Gauss Legendre quadrature
  # nl.info: if TRUE then the following additional information
  #          is printed to the console
  #          Number of iterations
  #          Optimal value of objective function
  #          If, on the other hand, nl.info is FALSE then this
  #          additional information is not printed to the console
  #
  # Output
  # A list with the following components
  # The inputs to this module, which are the following:
  # eta.ub, rho, d, n.ints, natural, alpha, n.nodes, nl.info
  # eta.vec: a 3-vector of values of eta
  # g.vec: a 3-vector of values of g(eta) for the values of
  #        eta in eta.vec
  # bsvec.matrix: a matrix with length(bsvec) rows and 3
  #               columns, where each column is the value
  #               bsvec for the corresponding value of eta
  #               from eta.vec
  # eta.star
  # eta.ub.holds: takes the value 1 if eta.star < eta.ub.holds;
  #               otherwise ieta.ub.holds takes the value 0
  #
  # Written by P. Kabaila in January 2023

  # eta.ub: a tight upper bound on the likely values
  #              of eta.star. Even if this is not an upper bound
  #              the method will still provide a good
  #              approximation to eta.star
  eta.ub <- -2.13

  gams <- seq(0, (d+2), by = 0.05)

  eta.vec <- rep(0, 3)
  g.vec <- rep(0, 3)

  # 1st evaluation of g
  eta.vec[1] <- eta.ub
  start <- start_standard_ci(d, n.ints, alpha)
  len.bsvec <- length(start)
  bsvec.matrix <- matrix(0, nrow = len.bsvec, ncol = 3)
  g.info <- g_fn(eta.vec[1], rho, alpha, gams,
                 d, n.ints, n.nodes, natural, start, nl.info)
  g.vec[1] <- g.info$g
  bsvec.matrix[,1] <- g.info$bsvec

  if (g.vec[1] > 0){
    eta.ub.holds <- 1
    eta.vec[2] <- eta.vec[1] - 0.15
  }else{
    eta.ub.holds <- 0
    eta.vec[2] <- eta.vec[1] + 0.05
  }

  # 2nd evaluation of g
  start <- bsvec.matrix[,1]
  g.info <- g_fn(eta.vec[2], rho, alpha, gams,
                 d, n.ints, n.nodes, natural, start, nl.info)
  g.vec[2] <- g.info$g
  bsvec.matrix[,2] <- g.info$bsvec

  # Secant method
  # The linear interpolation step is a weighted
  # average when g.vec[1] and g.vec[2] have
  # opposite signs.
  wt1 <- g.vec[2] / (g.vec[2] - g.vec[1])
  wt2 <- - g.vec[1] / (g.vec[2] - g.vec[1])
  # cat("eta.vec[1]=", eta.vec[1], ",  eta.vec[2]=", eta.vec[2], "\n")
  # cat("g.vec[1]=", g.vec[1], ",  g.vec[2]=", g.vec[2], "\n")
  # cat("wt1=", wt1, ", wt2=", wt2, "\n")
  eta.lin.interp <- wt1 * eta.vec[1] + wt2 * eta.vec[2]
  # cat("eta.lin.interp=", eta.lin.interp, "\n")
  eta.vec[3] <- eta.lin.interp

  # 3rd evaluation of g
  start <- wt1 * bsvec.matrix[,1] + wt2 * bsvec.matrix[,2]
  # cat("For 3rd evaluation of g, start=", start, "\n")
  g.info <- g_fn(eta.vec[3], rho, alpha, gams,
                 d, n.ints, n.nodes, natural, start, nl.info)
  g.vec[3] <- g.info$g
  bsvec.matrix[,3] <- g.info$bsvec

  # Apply Muller's method
  eta.star <- mullers_method(eta.vec, g.vec)

  list(eta.ub=eta.ub, rho=rho, d=d, n.ints=n.ints,
       natural=natural, alpha=alpha, n.nodes=n.nodes,
       nl.info=nl.info, eta.vec=eta.vec, g.vec=g.vec,
       bsvec.matrix=bsvec.matrix, eta.star=eta.star,
       eta.ub.holds=eta.ub.holds)

}



