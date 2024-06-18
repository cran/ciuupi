#' Computes the known correlation \eqn{\rho}
#' between  \eqn{\widehat{\theta}} and \eqn{\widehat{\tau}} from the
#' \eqn{p}-vectors
#' \eqn{a} and \eqn{c} and the
#' design matrix \eqn{X}
#'
#' Computes the known correlation \eqn{\rho}
#' between  \eqn{\widehat{\theta}}
#' and \eqn{\widehat{\tau}}. This correlation is computed from the
#' \eqn{p}-vectors
#' \eqn{a} and \eqn{c} and the \eqn{n \times p}
#' design matrix \eqn{X}, with linearly independent columns, using the formula
#' \eqn{\rho=a^{\top}(X^{\top} X)^{-1} c
#' /(v_{\theta} \, v_{\tau})^{1/2}}, where
#' \eqn{v_{\theta}
#' =a^{\top}(X^{\top} X)^{-1}a} and
#' \eqn{v_{\tau}
#' =c^{\top}(X^{\top} X)^{-1}c}.
#'
#' @param a The \eqn{p}-vector \eqn{a} that specifies the parameter of interest
#' \eqn{\theta =a^{\top}\beta}
#'
#' @param c The \eqn{p}-vector \eqn{c} used in the specification of the parameter
#' \eqn{\tau=c^{\top} \beta-t}. The uncertain prior
#'  information is that \eqn{\tau=0}
#'
#' @param X The \eqn{n \times p} design matrix \eqn{X}, with linearly
#' independent columns
#'
#' @return
#' The known correlation \eqn{\rho} between \eqn{\widehat{\theta}}
#' and \eqn{\widehat{\tau}}.
#'
#' @export
#'
#' @examples
#' a <- c(0, 2, 0, -2)
#' c <- c(0, 0, 0, 1)
#' x1 <- c(-1, 1, -1, 1)
#' x2 <- c(-1, -1, 1, 1)
#' X <- cbind(rep(1, 4), x1, x2, x1*x2)
#' rho <- acX_to_rho(a, c, X)
#' print(rho)
#'
acX_to_rho <- function(a, c, X){
  # This module computes rho, defined to be the
  # known correlation between the least squares
  # estimators of theta and tau, from a, c and X.
  #
  # Inputs
  # a: A p-vector used to specify the parameter of
  #    interest theta
  # c: A p-vector used to specify the parameter tau about which
  #    we have uncertain prior information that tau=0
  # X: The n by p design matrix
  #
  # Output
  # rho
  #
  # Written by P. Kabaila in June 2024,
  # using R code written by R. Mainzer

# The design matrix X is assumed to have linearly independent
# columns. This is not possible when the number of columns of X
# exceeds the number of rows of X.
  if (ncol(X) > nrow(X)){
    stop("p > n")
  }

# Use the QR decomposition of the matrix X
# to find the inverse of X'X
qrstr <- qr(X)
R <- qr.R(qrstr)
XTXinv <- chol2inv(R)

# Compute rho
rho <- (t(a) %*% XTXinv %*% c) /
  sqrt( t(a) %*% XTXinv %*% a %*% t(c) %*% XTXinv %*% c)

as.vector(rho)

}
