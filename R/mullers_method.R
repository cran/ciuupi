mullers_method <- function(x.vec, y.vec){
  # This module computes an approximation to
  # the value of x such that y(x) = 0, using
  # quadratic interpolation through the 3 values
  # (x.vec[1], y.vec[1]), (x.vec[2], y.vec[2])
  # & (x.vec[3], y.vec[3]). Here
  # y.vec[i] = y(x.vec[i]); i=1,2,3.
  # This is Muller's method.
  # Usually (x.vec[3], y.vec[3]) has been obtained by
  # linear interpolation from
  # (x.vec[1], y.vec[1]) and (x.vec[2], y.vec[2]).
  # We use the formula for this method given on page 214 of
  # Schwarz, H.R. (1989) Numerical analysis: a comprehensive
  # introduction. John Wiley, Chichester.
  #
  # Inputs
  # x.vec: vector of length 3
  # y.vec: vector of length 3
  #        y.vec[i] = y(x.vec[i]); i=1,2,3
  #
  # Output
  # Muller's method approximation to the value of
  # x such that y(x) = 0.
  #
  # Written by P. Kabaila in January 2023

  h2 <- x.vec[1] - x.vec[3]
  d2 <- y.vec[1] - y.vec[3]

  h1 <- x.vec[2] - x.vec[3]
  d1 <- y.vec[2] - y.vec[3]

  num.A <- h1 * d2 - h2 * d1
  denom.A <- h2 * h1 * (h2 - h1)
  A <- num.A / denom.A

  num.B <- h2^2 * d1 - h1^2 * d2
  denom.B <- h2 * h1 * (h2 - h1)
  B <- num.B / denom.B

  C <- y.vec[3]

  num <- -2 * sign(B) * C
  denom <- abs(B) + sqrt(B^2 - 4 * A * C)
  h <- num / denom

  x.vec[3] + h

}
