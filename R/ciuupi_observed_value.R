#' For given observed response vector \eqn{y},
#' compute the confidence interval that utilizes the
#' uncertain prior information (CIUUPI)
#'
#' If \eqn{\sigma} is provided then, for given observed response
#' vector \eqn{y},
#' compute the confidence interval, with minimum coverage
#' probability \eqn{1-\alpha}, for the parameter
#' \eqn{\theta =a^{\top}\beta} that
#' utilizes the uncertain prior information that the parameter
#' \eqn{\tau=c^{\top} \beta-t}
#' (specified by the vector \eqn{c} and the number
#' \eqn{t}) takes the value 0. If \eqn{\sigma}
#' is not provided
#' then, as long as \eqn{n-p \ge 30}, replace \eqn{\sigma} by its estimate
#' to compute an approximation to the CIUUPI for \eqn{\theta}.
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
#' @param alpha \eqn{1 - \alpha} is the desired minimum coverage probability of the
#' confidence interval for \eqn{\theta}
#'
#' @param bs.list A list that includes the following
#' components:
#' natural, d, q and the vector
#' bsvec (b(h),...,b((q-1)h), s(0),s(h)...,s((q-1)h)), where
#' h=d/q,
#' that specifies the CIUUPI for all possible values of the random
#' error variance and the observed response vector
#'
#' @param t The number \eqn{t} used to specify the parameter \eqn{\tau=c^{\top} \beta-t}.
#' The uncertain prior information is that \eqn{\tau = 0}
#'
#' @param y The \eqn{n}-vector of observed responses
#'
#' @param sig Standard deviation of the random error.
#' If a value is not specified then, provided that \eqn{n-p \ge 30},
#' \code{sig} is estimated from the data.
#'
#' @return If \eqn{\sigma} is provided then a data frame of the lower and upper
#' endpoints of
#' the confidence interval, with minimum coverage
#' probability \eqn{1-\alpha}, for the parameter
#' \eqn{\theta} that utilizes the
#' uncertain prior information that \eqn{\tau = 0}.
#' If \eqn{\sigma}
#' is not provided then, as long as \eqn{n-p \ge 30}, a data frame of the
#' lower and upper endpoints of
#' an approximation to this confidence interval.
#'
#'
#' @details
#' Suppose that \deqn{y = X \beta + \varepsilon}
#' where \eqn{y} is a random \eqn{n}-vector of
#' responses, \eqn{X} is a known \eqn{n \times p}
#' matrix with linearly
#' independent columns, \eqn{\beta} is an unknown parameter
#'  \eqn{p}-vector and
#' \eqn{\varepsilon} has components that are iid normally distributed
#' with zero mean and known variance.
#' Suppose that
#' \eqn{\theta=}\code{a}\eqn{^{\top}} \eqn{\beta} is the
#' parameter of interest, where \code{a} is a specified
#' vector. Let
#' \eqn{\tau=}\code{c}\eqn{^{\top} \beta -}\code{t},
#' where \code{c} is a specified vector,
#' \code{t} is a specified number and
#' \code{a} and \code{c} are
#' linearly independent vectors. Also suppose that we have
#' uncertain prior information that \eqn{\tau = 0}.
#' For given observed response
#' vector \code{y} and a design matrix \code{X},
#' \code{ciuupi_observed_value} computes the
#' confidence interval, with minimum coverage probability
#' 1\eqn{-}\code{alpha}, for \eqn{\theta}
#' that utilizes the uncertain prior information that
#' \eqn{\tau = 0}.
#'
#' The example below is described in Discussion 5.8 on
#' p.3426 of Kabaila and Giri (2009). This example is obtained
#' by extracting a \eqn{2 \times 2} factorial data set from the
#' \eqn{2^3} factorial data set described in Table 7.5
#' of Box et al. (1963).
#'
#'
#' @examples
#' a <- c(0, 2, 0, -2)
#' c <- c(0, 0, 0, 1)
#' x1 <- c(-1, 1, -1, 1)
#' x2 <- c(-1, -1, 1, 1)
#' X <- cbind(rep(1, 4), x1, x2, x1*x2)
#' alpha <- 0.05
#' t <- 0
#' y <- c(87.2, 88.4, 86.7, 89.2)
#' sig <- 0.8
#' ciuupi_observed_value(a, c, X, alpha, bs.list.example, t, y, sig=sig)
#'
#' @references
#'
#' Box, G.E.P., Connor, L.R., Cousins, W.R., Davies, O.L., Hinsworth, F.R., Sillitto, G.P. (1963)
#' The Design and Analysis
#' of Industrial Experiments, 2nd edition, reprinted. Oliver and Boyd, London.
#'
#' Kabaila, P. and Giri, K. (2009) Confidence intervals in regression utilizing
#' prior information.  Journal of Statistical Planning and Inference, 139,
#' 3419 - 3429.
#'
#' @export

ciuupi_observed_value <-
  function(a, c, X, alpha, bs.list, t, y, sig = NULL){

   # The design matrix X is assumed to have linearly independent
   # columns. This is not possible when p, the number of columns of X,
   # exceeds n, the number of rows of X.
    p <- ncol(X)
    n <- nrow(X)
    if (p > n){
      stop("p > n")
    }


  # Retrieve the information needed to specify the functions
  # b and s
  bsvec <- bs.list$bsvec
  d <- bs.list$d
  n.ints <- bs.list$n.ints
  natural <- bs.list$natural

  # Use the QR decomposition of the design matrix X
  # to find the inverse of X'X
  qrstr <- qr(X)
  R <- qr.R(qrstr)
  XTXinv <- chol2inv(R)

  # Find beta hat, theta hat and tau hat
  beta.hat <- XTXinv %*% t(X) %*% y
  theta.hat <- as.numeric(t(a) %*% beta.hat)
  tau.hat <- as.numeric(t(c) %*% beta.hat - t)

  # Find (variance of theta hat) / (sigma squared)
  v.theta <- as.numeric(t(a) %*% XTXinv %*% a)
  # Find (variance of tau hat) / (sigma squared)
  v.tau <- as.numeric(t(c) %*% XTXinv %*% c)

  # If sigma is not specified, find an estimate of sigma
  # if n - p >= 30
  if (is.null(sig)){
    if (n - p >=30){
      warning(paste("Random error sd not supplied by user and,",
                    "since n-p >= 30, this sd is estimated from the data,",
                    "resulting in an approximation to the CIUUPI"))
      sigsq.hat <- (t(y - X %*% beta.hat) %*% (y - X %*% beta.hat)) / (n - p)
      sig.hat <- as.numeric(sqrt(sigsq.hat))
      # Find approximation to gamma hat
      gam.hat <- tau.hat / (sig.hat * sqrt(v.tau))
      bsfuns <-
        bsspline(gam.hat, bsvec, alpha, d, n.ints, natural)
      new.ci.approx <- theta.hat - sig.hat * sqrt(v.theta) * bsfuns[, 2] +
        c(-1, 1) * sig.hat * sqrt(v.theta) * bsfuns[, 3]
      return(data.frame(lower =  new.ci.approx[1], upper =  new.ci.approx[2],
                 row.names = c("approx. ciuupi")))
    }else{
      stop(paste("Random error sd not supplied by user and",
                 "n-p < 30, i.e. n-p is too small"))
    }
  }

  # Find gamma hat
  gam.hat <- tau.hat / (sig * sqrt(v.tau))

  bsfuns <-
    bsspline(gam.hat, bsvec, alpha, d, n.ints, natural)

  new.ci <- theta.hat - sig * sqrt(v.theta) * bsfuns[, 2] +
             c(-1, 1) * sig * sqrt(v.theta) * bsfuns[, 3]

  data.frame(lower = new.ci[1], upper = new.ci[2],
                    row.names = c("ciuupi"))

}
