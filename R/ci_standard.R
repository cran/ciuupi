#' For given observed response vector \eqn{y}, compute
#' the standard \eqn{1 - \alpha} confidence interval
#'
#' If \eqn{\sigma} is provided then compute the standard \eqn{1 - \alpha}
#' confidence interval for \eqn{\theta}. If \eqn{\sigma}
#' is not provided
#' then, as long as \eqn{n-p \ge 30}, replace \eqn{\sigma} by its estimate
#' to compute an approximate \eqn{1 - \alpha} confidence interval for \eqn{\theta}.
#'
#' @param a The vector used to specify the parameter of interest
#' \eqn{\theta = a^{\top} \beta}
#'
#' @param X The known \eqn{n \times p} design matrix, with linearly
#' independent columns
#'
#' @param y The \eqn{n}-vector of observed responses
#'
#' @param alpha \eqn{1 - \alpha} is the coverage
#' probability of the standard confidence interval
#'
#' @param sig Standard deviation of the random error.
#' If a value is not specified then, provided that \eqn{n-p \ge 30},
#' \code{sig} is estimated from the data.
#'
#' @return If \eqn{\sigma} is provided then a data frame of the lower and upper
#' endpoints of the standard \eqn{1 - \alpha} confidence interval
#' for \eqn{\theta}. If \eqn{\sigma}
#' is not provided then, as long as \eqn{n-p \ge 30}, a data frame of the
#' lower and upper endpoints of
#' an approximation to this confidence interval.
#'
#' @details
#' Suppose that \deqn{y = X \beta + \varepsilon,}
#' where \eqn{y} is a random \eqn{n}-vector of responses, \eqn{X}
#' is a known \eqn{n \times p} matrix with linearly independent columns,
#' \eqn{\beta} is an unknown parameter \eqn{p}-vector, and
#' \eqn{\varepsilon \sim N(0, \, \sigma^2 \, I)}, with \eqn{\sigma^2} assumed known.
#' Suppose that the parameter of interest is \eqn{\theta = a^{\top} \beta}.
#' The R function \code{ci_standard} computes the standard \eqn{1 - \alpha}
#' confidence interval for \eqn{\theta}.
#'
#' The example below is described in Discussion 5.8 on
#' p.3426 of Kabaila and Giri (2009). This example is obtained
#' by extracting a \eqn{2 \times 2} factorial data set from the
#' \eqn{2^3} factorial data set described in Table 7.5
#' of Box et al. (1963).
#'
#' @examples
#' y <- c(87.2, 88.4, 86.7, 89.2)
#' x1 <- c(-1, 1, -1, 1)
#' x2 <- c(-1, -1, 1, 1)
#' X <- cbind(rep(1, 4), x1, x2, x1*x2)
#' a <- c(0, 2, 0, -2)
#' ci_standard(a, X, y, 0.05, sig = 0.8)
#'
#' @references
#' Box, G.E.P., Connor, L.R., Cousins, W.R., Davies, O.L., Hinsworth, F.R., Sillitto, G.P. (1963)
#' The Design and Analysis
#' of Industrial Experiments, 2nd edition, reprinted. Oliver and Boyd, London.
#'
#' Kabaila, P. and Giri, K. (2009) Confidence intervals in regression utilizing
#' prior information.  Journal of Statistical Planning and Inference, 139,
#' 3419 - 3429.
#'
#' @export

ci_standard <- function(a, X, y, alpha, sig = NULL){

  # The design matrix X is assumed to have linearly independent
  # columns. This is not possible when p, the number of columns of X,
  # exceeds n, the number of rows of X.
  p <- ncol(X)
  n <- nrow(X)
  if (p > n){
    stop("p > n")
  }

  # Use the QR decomposition of the design matrix X
  # to find the inverse of X'X
  qrstr <- qr(X)
  R <- qr.R(qrstr)
  XTXinv <- chol2inv(R)

  # Find beta hat and theta hat
  beta.hat <- XTXinv %*% t(X) %*% y
  theta.hat <- as.numeric(t(a) %*% beta.hat)

  # Find variance of theta hat on sigma squared
  v.theta <- as.numeric(t(a) %*% XTXinv %*% a)

  # If sigma is not specified, find an estimate of sigma
  # if n - p >= 30
  if (is.null(sig)){
    if (n - p >=30){
       warning(paste("Random error sd not supplied by user and,",
                     "since n-p >= 30, this sd is estimated from the data,",
                     "resulting in an approximate confidence interval"))
       sigsq.hat <- (t(y - X %*% beta.hat) %*% (y - X %*% beta.hat)) / (n - p)
       sig.hat <- as.numeric(sqrt(sigsq.hat))
       approx.ci <-
         theta.hat + c(-1, 1) * sig.hat * sqrt(v.theta) * stats::qnorm(1 - alpha/2)
       return(data.frame(lower = approx.ci[1], upper = approx.ci[2],
                  row.names = c("approximate")))
    }else{
      stop(paste("Random error sd not supplied by user and",
                 "n-p < 30, i.e. n-p is too small"))
    }
  }

  # Find the standard confidence interval
  standard.ci <-
    theta.hat + c(-1, 1) * sig * sqrt(v.theta) * stats::qnorm(1 - alpha/2)

  data.frame(lower = standard.ci[1], upper = standard.ci[2],
             row.names = c("standard"))

}
