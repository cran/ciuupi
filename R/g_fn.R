g_fn <- function(eta, rho, alpha, gams,
                 d, n.ints, n.nodes, natural, start, nl.info){
  # This module computes bsvec for given eta.
  # This is then used to evaluate g(eta) for given eta,
  # where g(eta) = gain - loss.
  #
  # Output
  # list with elements bsvec and g(eta)
  #
  # Written by P. Kabaila in January 2023

  lambda <- exp(eta)
  list.incl.bsvec <-
    optimize_b_s_given_lambda(lambda, rho, alpha, gams,
                              d, n.ints, n.nodes, natural, start, nl.info)
  bsvec <- list.incl.bsvec$bsvec

  quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
  c.alpha <- stats::qnorm(1 - alpha/2)
  temp <- sel_min_max(bsvec, d, n.ints, quad.info, c.alpha, natural)
  sel.min <- temp$sel.min
  sel.max <- temp$sel.max
  gain <- 1 - sel.min^2
  loss <- sel.max^2 - 1
  g <- gain - loss

  list(bsvec = bsvec, g = g)

}
