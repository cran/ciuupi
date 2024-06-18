d_goldilocks <- function(alpha, rho.vec){
  # This module computes the "goldilocks" value
  # of d (i.e. it computes the value of d that
  # is neither too large nor too small), for a
  # given value of alpha and rho.vec a vector
  # of values of rho.
  #
  # Inputs
  # alpha: the desired minimum coverage is 1-alpha
  # rho.vec: a vector of values of rho
  #
  # Output
  # The "goldilocks" value of d for the given value
  # of alpha and rho.vec, the vector of values of rho

  if (alpha <= 0.1){
    b.alpha <- 4.747653 - 17.386735 * alpha
  }else{
    b.alpha <- 3.00898
  }

  4.1 + b.alpha * abs(rho.vec)

}
