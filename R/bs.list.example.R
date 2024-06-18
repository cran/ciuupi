#' The list that specifies the CIUUPI for the example
#'
#' In this example, the dataset described in Table 7.5 of Box et al. (1963) is used.
#' The design matrix X is
#' specified by the command
#' \code{X <- cbind(rep(1,4), c(-1, 1, -1, 1), c(-1, -1, 1, 1), c(1, -1, -1, 1))}.
#' A description of the parameter of interest is given in Discussion 5.8, p.3426 of
#' Kabaila and Giri (2009).
#' The parameter of interest is \eqn{\theta = a^{\top} \beta},
#' where the column vector \eqn{a} is specified by the command
#' \code{a <- c(0, 2, 0, -2)}.
#' For this example, we have uncertain prior information that
#' \eqn{\tau = c^{\top} \beta = 0}, where the column vector \eqn{c} is specified by the command
#' \code{c <- c(0, 0, 0, 1)}.
#' The known correlation \eqn{\rho} between \eqn{\widehat{\theta}} and \eqn{\widehat{\tau}} is
#' computed using the command
#' \code{rho <- acX_to_rho(a, c, X)}.
#' The desired minimum coverage probability of the CIUUPI is \eqn{1 - \alpha}, where
#' \eqn{\alpha = 0.05}, which is specified by the command
#' \code{alpha <- 0.05}.
#' The CIUUPI is determined by \eqn{\alpha} and \eqn{\rho} and is found using the
#' command
#' \code{bs.list.example <- bs_ciuupi(alpha, rho)},
#' which takes about 5 minutes to run.
#'
#'
#' @references
#' Box, G.E.P., Connor, L.R., Cousins, W.R., Davies, O.L., Hinsworth, F.R., Sillitto, G.P. (1963)
#' The Design and Analysis
#' of Industrial Experiments, 2nd edition, reprinted. Oliver and Boyd, London.
#'
#' Kabaila, P. and Giri, K. (2009) Confidence intervals in regression utilizing
#' prior information.  Journal of Statistical Planning and Inference, 139,
#' 3419 - 3429.
"bs.list.example"
