# test.R

# How long does it take to compute the Gauss
# Legendre quadrature nodes and weights for
# 5 nodes?

# system.time(
# for (i in c(1:10000)){
# temp <- statmod::gauss.quad(5, kind="legendre")
# }
# )

# What does the function spline_s do?

# y <- c(0.1, 0.2, 0.3, 0.4, 0.3, 1, 2, 1, 1, 1, 1)
# d <- 6
# n.ints <- 6
# alpha <- 0.05
# c.alpha <- qnorm(1 - alpha/2)
# natural <- 1
#
# temp <- spline_s(y, d, n.ints, c.alpha, natural)
# temp
#
# # Does the function objective_v1.R seem to work?
#
# y <- c(0.1, 0.2, 0.3, 0.4, 0.3, 1, 2, 1, 1, 1, 1)
# lambda <- 0.1
# d <- 6
# n.ints <- 6
# alpha <- 0.05
# c.alpha <- qnorm(1 - alpha/2)
# natural <- 1
#
# n.nodes <- 5
# quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
#
# temp <- objective_v1(y, lambda, d, n.ints, quad.info, c.alpha, natural)
# temp
#

# Does the function compute_cov_legendre_v1.R seem to work?

# gam <- 0
# rho <- 0.7
# y <- c(0.1, 0.2, 0.3, 0.4, 0.3, 1, 2, 1, 1, 1, 1)
# d <- 6
# n.ints <- 6
# alpha <- 0.05
#
# natural <- 1
# n.nodes <- 5
# quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
# b.spl <- spline_b(y, d, n.ints, c.alpha, natural)
# s.spl <- spline_s(y, d, n.ints, c.alpha, natural)
#
# temp <- compute_cov_legendre_v1(gam, rho, y, d, n.ints, alpha,
#                                     quad.info, b.spl, s.spl)
# temp

# Does constraints_slsqp_gausslegendre_v1.R seem to work?

# gams <- seq(0,10, length.out=100)
# y <- c(0.1, 0.2, 0.3, 0.4, 0.3, 1, 2, 1, 1, 1, 1)
# rho <- 0.7
# d <- 6
# n.ints <- 6
# alpha <- 0.05
#
# n.nodes <- 5
# quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
#
# natural <- 1
#
# test.constr <- constraints_slsqp_gausslegendre_v1(gams, rho, y, d, n.ints,
#                                     alpha, quad.info, natural)
# plot(gams, test.constr, xlab="gamma",
#      ylab="CP constraint value", las=1)

# Does optimize_knots_v1.R seem to work?

# lambda <- 0.1
# rho <- 0.7
# alpha <- 0.05
# gams <- seq(0, 10, length.out=50)
# d <- 6
# n.ints <- 6
# n.nodes <- 5
# natural <- 1
#
# system.time(test <- optimize_knots(lambda, rho, alpha, gams,
#                            d, n.ints, n.nodes, natural))
# test
#
# system.time(test.v1 <- optimize_knots_v1(lambda, rho, alpha, gams,
#                                    d, n.ints, n.nodes, natural))
# test.v1

# # Compute the vector that specifies the CIUUPI and
# # prepare graphs of the b and s functions, the coverage
# # probability and the scaled expected length
#
# lambda <- 0.1
# rho <- 0.7
# alpha <- 0.05
# gams <- seq(0, 10, length.out=50)
# d <- 6
# n.ints <- 6
# n.nodes <- 5
# natural <- 1
#
# system.time(bsvec <- optimize_knots_v1(lambda, rho, alpha, gams,
#                                    d, n.ints, n.nodes, natural))
# bsvec
#
# # Compute the functions b and s that specify the CIUUPI on a grid of values
# splineval <- bsspline(seq(0, 8, by = 0.1), bsvec, alpha)
#
# # The first 5 values of bsvec are b(1),b(2),...,b(5).
# # The last 6 values are s(0),s(1),...,s(5).
# xseq <- seq(0, 6, by = 1)
# bvec <- c(0, bsvec[1:5], 0)
# svec <- c(bsvec[6:11], qnorm(1 - alpha/2))
# # Plot the functions b and s
# plot(seq(0, 8, by = 0.1), splineval[, 2], type = "l", main = "b function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, bvec, pch = 19, col = "blue")
# abline(h=0)
# plot(seq(0, 8, by = 0.1), splineval[, 3], type = "l", main = "s function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, svec, pch = 19, col = "blue")
# abline(h=qnorm(1-alpha/2))
#
# # Compute the coverage probability and scaled expected for a grid of values of gamma
# gam <- seq(0, 10, by = 0.1)
# cp <- cpciuupi(gam, bsvec, alpha, rho=rho)
# sel <- selciuupi(gam, bsvec, alpha, rho=rho)
# # Plot the coverage probability and squared scaled expected length
# plot(gam, cp, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Coverage Probability", col = "blue",
#      xlab = expression(paste("|", gamma, "|")), ylim = c(0.9495, 0.9505))
# abline(h = 1-alpha, lty = 2)
# plot(gam, sel^2, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Squared SEL", col = "blue",
#      xlab = expression(paste("|", gamma, "|")), ylim = c(0.83, 1.17))
# abline(h = 1, lty = 2)
#
# # Try out these computations with rho close to 1
#
# lambda <- 0.1
# rho <- 0.9
# alpha <- 0.05
# gams <- seq(0, 10, length.out=50)
# d <- 6
# n.ints <- 6
# n.nodes <- 5
# natural <- 1
#
# system.time(bsvec <- optimize_knots_v1(lambda, rho, alpha, gams,
#                                        d, n.ints, n.nodes, natural))
# bsvec
#
# # Compute the functions b and s that specify the CIUUPI on a grid of values
# splineval <- bsspline(seq(0, 8, by = 0.1), bsvec, alpha)
#
# # The first 5 values of bsvec are b(1),b(2),...,b(5).
# # The last 6 values are s(0),s(1),...,s(5).
# xseq <- seq(0, 6, by = 1)
# bvec <- c(0, bsvec[1:5], 0)
# svec <- c(bsvec[6:11], qnorm(1 - alpha/2))
# # Plot the functions b and s
# plot(seq(0, 8, by = 0.1), splineval[, 2], type = "l", main = "b function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, bvec, pch = 19, col = "blue")
# abline(h=0)
# plot(seq(0, 8, by = 0.1), splineval[, 3], type = "l", main = "s function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, svec, pch = 19, col = "blue")
# abline(h=qnorm(1-alpha/2))
#
# # Compute the coverage probability and scaled expected for a grid of values of gamma
# gam <- seq(0, 10, by = 0.1)
# cp <- cpciuupi(gam, bsvec, alpha, rho=rho)
# sel <- selciuupi(gam, bsvec, alpha, rho=rho)
# # Plot the coverage probability and squared scaled expected length
# plot(gam, cp, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Coverage Probability", col = "blue",
#      xlab = expression(paste("|", gamma, "|")), ylim = c(0.9495, 0.9505))
# abline(h = 1-alpha, lty = 2)
# plot(gam, sel^2, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Squared SEL", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# abline(h = 1, lty = 2)

# # What happens when we increase n.ints to 7?
#
# lambda <- 0.1
# rho <- 0.9
# alpha <- 0.05
# gams <- seq(0, 10, length.out=50)
# d <- 6
# n.ints <- 7
# n.nodes <- 5
# natural <- 1
#
# system.time(bsvec <- optimize_knots_v1(lambda, rho, alpha, gams,
#                                        d, n.ints, n.nodes, natural))
# bsvec
#
# # Compute the functions b and s that specify the CIUUPI on a grid of values
# splineval <- bsspline_v1(seq(0, 8, by = 0.1), bsvec, alpha, d, n.ints, natural)
#
# xseq <- seq(0, d, by = d/n.ints)
# bvec <- c(0, bsvec[1:(n.ints-1)], 0)
# svec <- c(bsvec[n.ints:(2*n.ints-1)], qnorm(1 - alpha/2))
# # Plot the functions b and s
# plot(seq(0, 8, by = 0.1), splineval[, 2], type = "l", main = "b function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, bvec, pch = 19, col = "blue")
# abline(h=0)
# plot(seq(0, 8, by = 0.1), splineval[, 3], type = "l", main = "s function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, svec, pch = 19, col = "blue")
# abline(h=qnorm(1-alpha/2))
#
# # Compute the coverage probability and scaled expected for
# # a grid of values of gamma
# gam <- seq(0, 10, by = 0.1)
# n.nodes <- 10
# cp <- cpciuupi_v1(gam, bsvec, alpha, n.nodes, d, n.ints, natural,
#                         rho = rho, a, c, x)
# sel <- selciuupi_v1(gam, bsvec, alpha, n.nodes, d, n.ints, natural, rho=rho)
# # Plot the coverage probability and squared scaled expected length
# plot(gam, cp, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Coverage Probability", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# abline(h = 1-alpha, lty = 2)
# plot(gam, sel^2, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Squared SEL", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# abline(h = 1, lty = 2)
#
# # What happens when we increase n.ints to 7 and d to 7?
#
# lambda <- 0.1
# rho <- 0.9
# alpha <- 0.05
# gams <- seq(0, 10, length.out=50)
# d <- 7
# n.ints <- 7
# n.nodes <- 5
# natural <- 1
#
# system.time(bsvec <- optimize_knots_v1(lambda, rho, alpha, gams,
#                                        d, n.ints, n.nodes, natural))
# bsvec
#
# # Compute the functions b and s that specify the CIUUPI on a grid of values
# splineval <- bsspline_v1(seq(0, 8, by = 0.1), bsvec, alpha, d, n.ints, natural)
#
# xseq <- seq(0, d, by = d/n.ints)
# bvec <- c(0, bsvec[1:(n.ints-1)], 0)
# svec <- c(bsvec[n.ints:(2*n.ints-1)], qnorm(1 - alpha/2))
# # Plot the functions b and s
# plot(seq(0, 8, by = 0.1), splineval[, 2], type = "l", main = "b function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, bvec, pch = 19, col = "blue")
# abline(h=0)
# plot(seq(0, 8, by = 0.1), splineval[, 3], type = "l", main = "s function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, svec, pch = 19, col = "blue")
# abline(h=qnorm(1-alpha/2))
#
# # Compute the coverage probability and scaled expected for
# # a grid of values of gamma
# gam <- seq(0, 10, by = 0.1)
# n.nodes <- 10
# cp <- cpciuupi_v1(gam, bsvec, alpha, n.nodes, d, n.ints, natural,
#                   rho = rho, a, c, x)
# sel <- selciuupi_v1(gam, bsvec, alpha, n.nodes, d, n.ints, natural, rho=rho)
# # Plot the coverage probability and squared scaled expected length
# plot(gam, cp, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Coverage Probability", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# abline(h = 1-alpha, lty = 2)
# plot(gam, sel^2, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Squared SEL", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# abline(h = 1, lty = 2)
#
# # What happens when we set n.ints to 8 and d to 7?
#
# lambda <- 0.1
# rho <- 0.9
# alpha <- 0.05
# gams <- seq(0, 10, length.out=50)
# d <- 7
# n.ints <- 8
# n.nodes <- 5
# natural <- 1
#
# system.time(bsvec <- optimize_knots_v1(lambda, rho, alpha, gams,
#                                        d, n.ints, n.nodes, natural))
# bsvec
#
# # Compute the functions b and s that specify the CIUUPI on a grid of values
# splineval <- bsspline_v1(seq(0, 8, by = 0.1), bsvec, alpha, d, n.ints, natural)
#
# xseq <- seq(0, d, by = d/n.ints)
# bvec <- c(0, bsvec[1:(n.ints-1)], 0)
# svec <- c(bsvec[n.ints:(2*n.ints-1)], qnorm(1 - alpha/2))
# # Plot the functions b and s
# plot(seq(0, 8, by = 0.1), splineval[, 2], type = "l", main = "b function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, bvec, pch = 19, col = "blue")
# abline(h=0)
# plot(seq(0, 8, by = 0.1), splineval[, 3], type = "l", main = "s function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, svec, pch = 19, col = "blue")
# abline(h=qnorm(1-alpha/2))
#
# # Compute the coverage probability and scaled expected for
# # a grid of values of gamma
# gam <- seq(0, 10, by = 0.1)
# n.nodes <- 10
# cp <- cpciuupi_v1(gam, bsvec, alpha, n.nodes, d, n.ints, natural,
#                   rho = rho, a, c, x)
# sel <- selciuupi_v1(gam, bsvec, alpha, n.nodes, d, n.ints, natural, rho=rho)
# # Plot the coverage probability and squared scaled expected length
# plot(gam, cp, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Coverage Probability", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# abline(h = 1-alpha, lty = 2)
# plot(gam, sel^2, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Squared SEL", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# abline(h = 1, lty = 2)
#
# # What happens when we set n.ints to 9 and d to 7?
#
# lambda <- 0.1
# rho <- 0.9
# alpha <- 0.05
#
# The following was used in the original version
# of ciuupi
# gams <- seq(0, 8, by = 0.05)
#
# gams <- seq(0, 10, length.out=50)
# d <- 7
# n.ints <- 9
# n.nodes <- 5
# natural <- 1
#
# system.time(bsvec <- optimize_knots_v1(lambda, rho, alpha, gams,
#                                        d, n.ints, n.nodes, natural))
# bsvec
#
# # Compute the functions b and s that specify the CIUUPI on a grid of values
# splineval <- bsspline_v1(seq(0, 8, by = 0.1), bsvec, alpha, d, n.ints, natural)
#
# xseq <- seq(0, d, by = d/n.ints)
# bvec <- c(0, bsvec[1:(n.ints-1)], 0)
# svec <- c(bsvec[n.ints:(2*n.ints-1)], qnorm(1 - alpha/2))
# # Plot the functions b and s
# plot(seq(0, 8, by = 0.1), splineval[, 2], type = "l", main = "b function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, bvec, pch = 19, col = "blue")
# abline(h=0)
# plot(seq(0, 8, by = 0.1), splineval[, 3], type = "l", main = "s function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, svec, pch = 19, col = "blue")
# abline(h=qnorm(1-alpha/2))
#
# # Compute the coverage probability and scaled expected for
# # a grid of values of gamma
# gam <- seq(0, 10, by = 0.1)
# n.nodes <- 10
# cp <- cpciuupi_v1(gam, bsvec, alpha, n.nodes, d, n.ints, natural,
#                   rho = rho, a, c, x)
# sel <- selciuupi_v1(gam, bsvec, alpha, n.nodes, d, n.ints, natural, rho=rho)
# # Plot the coverage probability and squared scaled expected length
# plot(gam, cp, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Coverage Probability", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# abline(h = 1-alpha, lty = 2)
# plot(gam, sel^2, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Squared SEL", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# abline(h = 1, lty = 2)

# #********************************************
# # Check that plots_b_s_cp_sel.R seems to work
#
# lambda <- 0.1
# rho <- 0.98
# alpha <- 0.05
#
# # The following was used in the original version
# # of ciuupi
# # gams <- seq(0, 8, by = 0.05)
# #gams <- seq(0, 10, length.out=50)
#
# d <- 6
# n.ints <- 6
# gams <- seq(0, (d+2), by = 0.05)
# n.nodes <- 5
# natural <- 1
# nl.info <- TRUE
#
# # Specify a starting value for the numerical nonlinear constrained
# # optimization
# start <- start_standard_ci(d, n.ints, alpha)
#
# system.time(list.incl.bs <- optimize_b_s_given_lambda(lambda, rho, alpha, gams,
#                                   d, n.ints, n.nodes, natural, start, nl.info))
# list.incl.bs
# bsvec <- list.incl.bs$bsvec
#
# system.time(plots_b_s_cp_sel(bsvec, d, n.ints, alpha, natural))
#
# # Compute the functions b and s that specify the CIUUPI on a grid of values
# splineval <- bsspline_v1(seq(0, 8, by = 0.02), bsvec, alpha, d, n.ints, natural)
#
# xseq <- seq(0, d, by = d/n.ints)
# bvec <- c(0, bsvec[1:(n.ints-1)], 0)
# svec <- c(bsvec[n.ints:(2*n.ints-1)], qnorm(1 - alpha/2))
# # Plot the functions b and s
# plot(seq(0, 8, by = 0.02), splineval[, 2], type = "l", main = "b function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, bvec, pch = 19, col = "blue")
# abline(h=0)
# plot(seq(0, 8, by = 0.02), splineval[, 3], type = "l", main = "s function",
#      ylab = " ", las = 1, lwd = 2, xaxs = "i", col = "blue", xlab = "x")
# points(xseq, svec, pch = 19, col = "blue")
# abline(h=qnorm(1-alpha/2))
#
# # Compute the coverage probability and scaled expected for
# # a grid of values of gamma
# gam <- seq(0, (d+4), by = 0.01)
#
# n.nodes <- 10
#
# cp <- cpciuupi_v1(gam, bsvec, alpha, n.nodes, d, n.ints, natural,
#                   rho = rho, a, c, x)
# sel <- selciuupi_v1(gam, bsvec, alpha, n.nodes, d, n.ints, natural, rho=rho)
# # Plot the coverage probability and squared scaled expected length
# plot(gam, cp, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Coverage Probability", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# abline(h = 1-alpha, lty = 2)
# min.cp <- min(cp)
# mtext(paste("minimum CP =", min.cp))
#
# plot(gam, sel^2, type = "l", lwd = 2, ylab = "", las = 1, xaxs = "i",
#      main = "Squared SEL", col = "blue",
#      xlab = expression(paste("|", gamma, "|")))
# max.sel.squared <- max(sel)^2
# min.sel.squared <- min(sel)^2
# mtext(paste("gain=", 1 - min.sel.squared, ",  loss=", max.sel.squared -1))
# abline(h = 1, lty = 2)
#


# Does sel_min_max.R seem to work?

# quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
# c.alpha <- qnorm(1 - alpha/2)
#
#
# system.time(sel.min.max <-
#               sel_min_max(bsvec, d, n.ints, quad.info, c.alpha, natural))
#
# sel.min.max


#







# #--------------------
#
# lambda <- 0.1
# rho <- 0.98
# alpha <- 0.05
#
# d <- 6
# n.ints <- 6
# gams <- seq(0, (d+2), by = 0.05)
# n.nodes <- 5
# natural <- 1
# nl.info <- TRUE
#
# # Choose the starting value of bsvec for the numerical nonlinear
# # constrained optimization
# start <- start_standard_ci(d, n.ints, alpha)
#
# system.time(list.incl.bsvec <-
#               optimize_b_s_given_lambda(lambda, rho, alpha, gams,
#                     d, n.ints, n.nodes, natural, start, nl.info))
# bsvec <- list.incl.bsvec$bsvec
# system.time(plots_b_s_cp_sel(bsvec, d, n.ints, alpha, natural))
#
# #----------------
#
# d.new <- 7
# n.ints.new <- 9
# start.new <- start_from_bsvec(bsvec, d, n.ints, alpha, natural,
#                              d.new, n.ints.new)
# system.time(list.incl.bsvec.new <-
#               optimize_b_s_given_lambda(lambda, rho, alpha, gams,
#                     d.new, n.ints.new, n.nodes, natural, start.new, nl.info))
# bsvec.new <- list.incl.bsvec.new$bsvec
#
# system.time(plots_b_s_cp_sel(bsvec.new, d.new, n.ints.new, alpha, natural))
#
# #----------------
#
# lambda <- 0.1
# rho <- 0.98
# alpha <- 0.05
# gams <- seq(0, (d+2), by = 0.05)
# n.nodes <- 5
# natural <- 1
# nl.info <- TRUE
#
# d <- 7
# n.ints <- 9
# start <- start_standard_ci(d, n.ints, alpha)
# system.time(list.incl.bsvec <-
#               optimize_b_s_given_lambda(lambda, rho, alpha, gams,
#                     d, n.ints, n.nodes, natural, start, nl.info))
# bsvec <- list.incl.bsvec$bsvec

# n.nodes <- 10
# quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
# c.alpha <- qnorm(1 - alpha/2)
# sel_min_max(bsvec, d, n.ints, quad.info, c.alpha, natural)
# sel.min <- sel.min.max$sel.min
# sel.max <- sel.min.max$sel.max
#
# gain <-  1 - sel.min^2
# loss <- sel.max^2 - 1
#
# options(digits=10)
# gain
# loss

# system.time(plots_b_s_cp_sel(bsvec, d, n.ints, alpha, natural))
#
# lambda <- 0.15
#
# start <- bsvec
# system.time(list.incl.bsvec.new <-
#               optimize_b_s_given_lambda(lambda, rho, alpha, gams,
#                               d, n.ints, n.nodes, natural, start, nl.info))
# bsvec.new <- list.incl.bsvec.new$bsvec

# n.nodes <- 10
# quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
# c.alpha <- qnorm(1 - alpha/2)
# sel_min_max(y, d, n.ints, quad.info, c.alpha, natural)
# sel.min <- sel.min.max$sel.min
# sel.max <- sel.min.max$sel.max
#
# gain <-  1 - sel.min^2
# loss <- sel.max^2 - 1
#
# options(digits=10)
# gain
# loss

# system.time(plots_b_s_cp_sel(bsvec.new, d, n.ints, alpha, natural))

# #-------------------
# rho <- 0.7
# alpha <- 0.05
# n.nodes <- 5
# natural <- 1
# nl.info <- TRUE
#
# d <- 6
# gams <- seq(0, (d+2), by = 0.05)
# n.ints <- 6
# start <- start_standard_ci(d, n.ints, alpha)
#
# c.alpha <- qnorm(1 - alpha/2)
#
# n.nodes <- 5
# quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
#
# lambda.grid <- seq(0.05, 0.3, length.out=10)
# len.lambda.grid <- length(lambda.grid)
#
# gain.grid <- rep(0, len.lambda.grid)
# loss.grid <- rep(0, len.lambda.grid)
#
# system.time(
# for (i in c(1:len.lambda.grid)){
#   lambda <- lambda.grid[i]
#   list.incl.bsvec <-
#                 optimize_b_s_given_lambda(lambda, rho, alpha, gams,
#                             d, n.ints, n.nodes, natural, start, nl.info)
#   bsvec <- list.incl.bsvec$bsvec
#   plots_b_s_cp_sel(bsvec, d, n.ints, alpha, natural)
#   temp <- sel_min_max(bsvec, d, n.ints, quad.info, c.alpha, natural)
#   sel.min <- temp$sel.min
#   sel.max <- temp$sel.max
#
#   gain.grid[i] <-  1 - sel.min^2
#   loss.grid[i] <- sel.max^2 - 1
# }
# )
#
#
# ymin <- min(gain.grid, loss.grid)
# ymax <- max(gain.grid, loss.grid)
#
#
#
# plot(lambda.grid, gain.grid, type="l", las=1, xlab="lambda",
#      ylab="gain in black, loss in red", ylim=c(ymin, ymax))
# lines(lambda.grid, loss.grid, type="l", col="red")
# mtext(paste("rho=", rho, ",  alpha=", alpha,
#       ",  d=", d, ",  n.ints=", n.ints, ", n.nodes=", n.nodes,
#       ",  natural=", natural), cex=0.75)
#
#
# plot(lambda.grid, (gain.grid - loss.grid), type="l", las=1, xlab="lambda",
#      ylab="gain - loss")
# mtext(paste("rho=", rho, ",  alpha=", alpha,
#             ",  d=", d, ",  n.ints=", n.ints, ", n.nodes=", n.nodes,
#             ",  natural=", natural), cex=0.75)
# abline(h = 0)
#
#
# plot(log(lambda.grid), (gain.grid - loss.grid), type="l", las=1,
#      xlab="log(lambda)",
#      ylab="gain - loss")
# mtext(paste("rho=", rho, ",  alpha=", alpha,
#             ",  d=", d, ",  n.ints=", n.ints, ", n.nodes=", n.nodes,
#             ",  natural=", natural), cex=0.75)
# abline(h = 0)
#
#
# # #
# Time taken to compute the optimum functions b and s for
# given tuning constant lambda


#
# rho <- 0.7
# alpha <- 0.05
# n.nodes <- 5
# natural <- 1
# nl.info <- TRUE
#
# d <- 6
# gams <- seq(0, (d+2), by = 0.05)
# n.ints <- 6
# start <- start_standard_ci(d, n.ints, alpha)
#
# c.alpha <- qnorm(1 - alpha/2)
#
# n.nodes <- 5
# quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
#
# lambda <- 0.107
#
# system.time(list.incl.bsvec <-
#   optimize_b_s_given_lambda(lambda, rho, alpha, gams,
#                             d, n.ints, n.nodes, natural, start, nl.info))
# bsvec <- list.incl.bsvec$bsvec
# plots_b_s_cp_sel(bsvec, d, n.ints, alpha, natural)

#----

#
#
# rho <- 0.7
# alpha <- 0.05
# n.nodes <- 5
# natural <- 1
# nl.info <- TRUE
#
# d <- 6
# gams <- seq(0, (d+2), by = 0.05)
# n.ints <- 6
#
# c.alpha <- qnorm(1 - alpha/2)
#
# n.nodes <- 5
# quad.info <- statmod::gauss.quad(n.nodes, kind="legendre")
#
# eta.initial <- -2.15
#
# lambda <- exp(eta.initial)
# start <- start_standard_ci(d, n.ints, alpha)
# start
#
# system.time(list.incl.bsvec <-
#               optimize_b_s_given_lambda(lambda, rho, alpha, gams,
#                     d, n.ints, n.nodes, natural, start, nl.info))
# bsvec <- list.incl.bsvec$bsvec
#
# temp <- sel_min_max(bsvec, d, n.ints, quad.info, c.alpha, natural)
# sel.min <- temp$sel.min
# sel.max <- temp$sel.max
#
# gain <- 1 - sel.min^2
# loss <- sel.max^2 - 1
# gain
# loss
#
# eta.vec <- rep(0, 3)
# g.vec <- rep(0, 3)
# len.bsvec <- length(bsvec)
# bsvec.matrix <- matrix(0, nrow = len.bsvec, ncol = 3)
#
# eta.vec[1] <- eta.initial
# g.vec[1] <- gain - loss
#
# bsvec.matrix[,1] <- bsvec
# bsvec.matrix
#
# g.vec[1]
#
# if (g.vec[1] > 0){
#   eta.vec[2] <- eta.vec[1] - 0.15
# }else{
#   eta.vec[2] <- eta.vec[1] + 0.15
# }
#
# lambda <- exp(eta.vec[2])
# # start <- start_standard_ci(d, n.ints, alpha)
# start <- bsvec.matrix[,1]
#
# system.time(list.incl.bsvec <-
#               optimize_b_s_given_lambda(lambda, rho, alpha, gams,
#                                         d, n.ints, n.nodes, natural, start, nl.info))
# bsvec <- list.incl.bsvec$bsvec
#
# temp <- sel_min_max(bsvec, d, n.ints, quad.info, c.alpha, natural)
# sel.min <- temp$sel.min
# sel.max <- temp$sel.max
#
# gain <- 1 - sel.min^2
# loss <- sel.max^2 - 1
# gain
# loss
#
# g.vec[2] <- gain - loss
# bsvec.matrix[,2] <- bsvec
#
# plot(eta.vec[1:2], g.vec[1:2], type = "b", xlab = "eta",
#      ylab = "g(eta)")
# abline(h = 0)
#
# # Express the linear interpolation as a weighted average
# wt1 <- g.vec[2] / (g.vec[2] - g.vec[1])
# wt2 <- - g.vec[1] / (g.vec[2] - g.vec[1])
# eta.lin.interp1 <- wt1 * eta.vec[1] + wt2 * eta.vec[2]
#
# # Check that the usual formula gives the same result
# # term <- (eta.vec[2] - eta.vec[1]) / (g.vec[2] - g.vec[1])
# # eta.lin.interp <- eta.vec[2] - g.vec[2] * term
# # eta.lin.interp1 - eta.lin.interp
#
# eta.vec[3] <- eta.lin.interp1
#
# abline(v = eta.vec[3])
#
# #---------------
#
# lambda <- exp(eta.vec[3])
# # start <- start_standard_ci(d, n.ints, alpha)
# # start <- bsvec.matrix[,2]
# start <- wt1 * bsvec.matrix[,1] + wt2 * bsvec.matrix[,2]
#
# system.time(list.incl.bsvec <-
#               optimize_b_s_given_lambda(lambda, rho, alpha, gams,
#                                         d, n.ints, n.nodes, natural, start, nl.info))
# bsvec <- list.incl.bsvec$bsvec
#
# temp <- sel_min_max(bsvec, d, n.ints, quad.info, c.alpha, natural)
# sel.min <- temp$sel.min
# sel.max <- temp$sel.max
#
# gain <- 1 - sel.min^2
# loss <- sel.max^2 - 1
# gain
# loss
#
# g.vec[3] <- gain - loss
#
# eta.vec[1:3]
# g.vec[1:3]
#
# plot(sort(eta.vec[1:3]), sort(g.vec[1:3]), type = "b", xlab = "eta",
#      ylab = "g(eta)")
# abline(h = 0)
#
# plots_b_s_cp_sel(bsvec, d, n.ints, alpha, natural)
#
# eta.star <- mullers_method(eta.vec[1:3], g.vec[1:3])
#
# eta.vec
# g.vec
# eta.star
#
# test <- g_fn(eta.star, rho, alpha, gams,
#         d, n.ints, n.nodes, natural, start, nl.info)
#
# test$g
#
# test <- g_fn(eta.vec[3], rho, alpha, gams,
#              d, n.ints, n.nodes, natural, start, nl.info)
#
# test$g
#


# eta.ub <- -2.15
# rho <- 0.7
# d <- 6
# n.ints <- 6
# alpha <- 0.05
# n.nodes <- 5
# natural <- 1
# n.nodes <- 5
# nl.info <- TRUE
#
#
# system.time(test <-
#     eta_star(eta.ub, rho, d, n.ints, natural, alpha, n.nodes, nl.info))
#
# test
#
# eta <- test$eta.star
# start <- test$bsvec.matrix[,3]
# system.time(test1 <- g_fn(eta, rho, alpha, gams,
#                              d, n.ints, n.nodes, natural, start, nl.info))
# test1
#
# bsvec <- test1$bsvec
# plots_b_s_cp_sel(bsvec, d, n.ints, alpha, natural)

#****************************
# What does the graph of g(eta)
# look like for eta in [-2.3, -2]
# for small rho?

# rho <- 0.05
# alpha <- 0.02
# d <- 7
# n.ints <- 9
# gams <- seq(0, (d+2), by = 0.05)
# n.nodes <- 5
# natural <- 1
#
# eta.vec <- seq(-2.3, -2, length.out=10)
# len.eta.vec <- length(eta.vec)
#
# g.vec <- rep(0, len.eta.vec)
#
# start <- start_standard_ci(d, n.ints, alpha)
# system.time(
# for (k in c(1:len.eta.vec)){
#   eta <- eta.vec[k]
#   tmp <- g_fn(eta, rho, alpha, gams,
#                  d, n.ints, n.nodes, natural, start, nl.info)
#   g.vec[k] <- tmp$g
# #  start <- tmp$bsvec
#   cat("Finished step", k, " of", len.eta.vec, "steps", "\n")
# }
# )
#
# par(pty="s")
# plot(eta.vec, g.vec, type="l", xlab="eta", ylab="", las=1)
# mtext(paste("rho=", rho, ", alpha=", alpha, ", d=", d, ", n.ints=", n.ints,
#             ", n.nodes=", n.nodes, ", natural=", natural))
# title(main="g(eta)")
# abline(h=0)

#*************************
# For each alpha in {0.1, 0.05, 0.02, 0.01}
# compute eta.star for
# rho in {0.1, 0.2, 0.5, 0.7, 0.9, 0.98}.
# We set d=7, n.ints=9, natural=1 and
# nl.info=FALSE

# alpha <- 0.2
# rho.vec <- c(0.1, 0.2, 0.5, 0.7, 0.9, 0.98)
#
# n.nodes <- 5
# nl.info <- FALSE
# # nl.info <- TRUE
# eta.ub <- -2.15
# d <- 7
# n.ints <- 9
# natural <- 1
#
# len.rho.vec <- length(rho.vec)
# eta.star.vec <- rep(0, len.rho.vec)
#
# cat("alpha=", alpha, "\n")
# system.time(
# for (k in c(1:len.rho.vec)){
#   rho <- rho.vec[k]
#   cat("rho=", rho, "\n")
#   tmp <-
#     eta_star(eta.ub, rho, d, n.ints, natural, alpha, n.nodes, nl.info)
#   eta.star <- tmp$eta.star
#   eta.star.vec[k] <- eta.star
#   cat("Computation", k, "of", len.rho.vec, "computations completed",
#       ", eta.star=", eta.star, "\n", "\n")
# }
# )
#
# par(pty="m")
# plot(rho.vec, eta.star.vec, type="b", xlim=c(0,1),
#      xaxs="i", xlab="rho", ylab="eta.star", las=1)
# mtext(paste("alpha=", alpha, ", d=", d, ", n.ints=", n.ints, ", natural=",
#             natural))
#
# saveRDS(eta.star.vec, "eta.star.alpha20.d7.nints9", ascii = TRUE)
#
# #------------------------------------
#
# rho.vec <- c(0.1, 0.2, 0.5, 0.7, 0.9, 0.98)
#
# eta.star.alpha01 <- readRDS("eta.star.alpha01.d7.nints9")
# eta.star.alpha01
#
# eta.star.alpha02 <- readRDS("eta.star.alpha02.d7.nints9")
# eta.star.alpha02
#
# eta.star.alpha05 <- readRDS("eta.star.alpha05.d7.nints9")
# eta.star.alpha05
#
# eta.star.alpha07 <- readRDS("eta.star.alpha07.d7.nints9")
# eta.star.alpha07
#
# eta.star.alpha10 <- readRDS("eta.star.alpha10.d7.nints9")
# eta.star.alpha10
#
# eta.star.max <- max(eta.star.alpha01, eta.star.alpha02, eta.star.alpha05,
#                     eta.star.alpha07, eta.star.alpha10)
# eta.star.max
#
# eta.star.min <- min(eta.star.alpha01, eta.star.alpha02, eta.star.alpha05,
#                     eta.star.alpha07, eta.star.alpha10)
# eta.star.min
#
# plot(rho.vec, eta.star.alpha01, type="b", xlim=c(0,1),
#      ylim=c(eta.star.min, eta.star.max), las=1, xlab="rho", ylab="eta.star",
#      col="red")
#
# lines(rho.vec, eta.star.alpha02, type="b", col="purple")
# lines(rho.vec, eta.star.alpha05, type="b", col="blue")
# lines(rho.vec, eta.star.alpha07, type="b", col="green")
# lines(rho.vec, eta.star.alpha10, type="b", col="orange")
# mtext("alpha: 0.01(red), 0.02(purple), 0.05(blue), 0.07(green), 0.10(orange)")
#
# eta.star.max - eta.star.min
#
# #------------------
#
# eta.star.alpha20 <- readRDS("eta.star.alpha02.d7.nints9")
# eta.star.alpha20
#
# eta.star.max > eta.star.alpha20
# eta.star.min < eta.star.alpha20
#
# #------------------
#
# max(eta.star.alpha05)
#
# max(eta.star.alpha05) - eta.star.min
#
# eta.ub <- -2.13
#
# eta.ub - eta.star.min
# eta.star.max - eta.ub
#

# #*********************************
# # Check that the module eta.star
# # seems to be working OK. Also
# # check the timing.
#
# alpha <- 0.05
# rho <- 0.1
#
# d <- 5
# n.ints <- 6
# natural <- 1
#
# n.nodes <- 5
# nl.info <- TRUE
#
# system.time(
# temp <- eta_star(rho, d, n.ints, natural, alpha, n.nodes, nl.info)
# )
#
# temp$eta.star
#
#
# eta <- temp$eta.star
# start <- temp$bsvec.matrix[,3]
#
# system.time(test <- g_fn(eta, rho, alpha, gams,
#                              d, n.ints, n.nodes, natural, start, nl.info)
# )
#
# bsvec <- test$bsvec
#
# # plots_b_s_cp_sel(bsvec, d, n.ints, alpha, natural)
#
# #*****************************
#
# d1 <- 5
# n.ints1 <- 8
# natural1 <- 1
#
# system.time(
#   temp1 <- eta_star(rho, d1, n.ints1, natural1, alpha, n.nodes, nl.info)
# )
#
# temp1$eta.star
#
# eta1 <- temp1$eta.star
# start <- temp1$bsvec.matrix[,3]
#
# system.time(test1 <- g_fn(eta1, rho, alpha, gams,
#                           d1, n.ints1, n.nodes, natural1, start, nl.info)
# )
#
#
# bsvec1 <- test1$bsvec
#
# compare_2_bsvecs(bsvec, d, n.ints, natural,
#                              bsvec1, d1, n.ints1, natural1, alpha)
#
#

# alpha <- 0.05
#
# natural <- 1
# n.nodes <- 5
# nl.info <- FALSE
# n.ints = 15
#
# d <- 6
# n.ints <- 6
# temp <-
#   eta_star(rho, d, n.ints, natural, alpha, n.nodes, nl.info)
# eta <- temp$eta.star
# start <- temp$bsvec.matrix[,3]
# bsvec.g.list <- g_fn(eta, rho, alpha, gams,
#                      d, n.ints, n.nodes, natural, start, nl.info)
# bsvec <- bsvec.g.list$bsvec
#
# dev.new(width = 21, height = 30 , unit = cm,
# noRStudioGD = TRUE)
# # plots_b_s(bsvec, d, n.ints, alpha, natural)
# plots_cp_sel(bsvec, d, n.ints, alpha, natural)
# dev.off()

#**************************
# Choosing a good value of d -
# neither too small nor too large.
# For a given alpha, this is done
# by considering pairs (rho, d),
# where, in each case, the value
# of d exceeds the "goldylocks"
# value of d for that rho. We
# then look at where the optimized
# b and s functions flatten out
# and check that otherwsise performance
# is adequate. In each case, we choose
# n.ints = 15, which leads to a
# computation time of 12.4 minutes for
# each pair (rho, d).

# alpha <- 0.01
#
# natural <- 1
# n.nodes <- 5
# nl.info <- FALSE
#
# rho.d.matrix <- matrix(0, ncol=10, nrow=2)
#
# rho.d.matrix[,1] <- c(0.02, 5)
# rho.d.matrix[,2] <- c(0.1, 5.5)
# rho.d.matrix[,3] <- c(0.2, 5.5)
# rho.d.matrix[,4] <- c(0.4, 6.5)
# rho.d.matrix[,5] <- c(0.5, 7)
# rho.d.matrix[,6] <- c(0.6, 7)
# rho.d.matrix[,7] <- c(0.7, 8)
# rho.d.matrix[,8] <- c(0.8, 8)
# rho.d.matrix[,9] <- c(0.9, 8.5)
# rho.d.matrix[,10] <- c(0.95, 8.5)
#
# rho.d.matrix
# ncol(rho.d.matrix)
#
# n.ints <- 15
#
# system.time(
# for (k in c(1: ncol(rho.d.matrix))){
#   rho <- rho.d.matrix[1, k]
#   d <- rho.d.matrix[2, k]
#   temp <-
#     eta_star(rho, d, n.ints, natural, alpha, n.nodes, nl.info)
#   eta <- temp$eta.star
#   start <- temp$bsvec.matrix[,3]
#   bsvec.g.list <- g_fn(eta, rho, alpha, gams,
#               d, n.ints, n.nodes, natural, start, nl.info)
#   bsvec <- bsvec.g.list$bsvec
#   dev.new(width = 21, height = 30 , unit = cm,
#           noRStudioGD = TRUE)
#   plots_b_s(bsvec, d, n.ints, alpha, natural)
#   dev.new(width = 21, height = 30 , unit = cm,
#           noRStudioGD = TRUE)
#   plots_cp_sel(bsvec, d, n.ints, alpha, natural)
#   cat("Computation", k, "of", ncol(rho.d.matrix),
#       "computations completed", "\n")
# }
# )
#
# rho.alpha01.grid <- c(0.02, 0.1, 0.2, 0.4, 0.5)
# b.opt.alpha01.grid <- c(4.5, 4.58333, 4.7667, 5.85, 6.5333)
# model.alpha01 <- lm(b.opt.alpha01.grid ~ rho.alpha01.grid)
# model.alpha01
# model.alpha01$coefficients
# ymin.alpha01 <- model.alpha01$coefficients[1]
# ymin.alpha01
# ymax.alpha01 <- model.alpha01$coefficients[1] +
#   model.alpha01$coefficients[2]
# ymax.alpha01
#
# rho.alpha02.grid <- c(0.02, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8)
# b.opt.alpha02.grid <- c(4.333, 4.76667, 4.7667, 5.85, 6.06667, 6.5333,
#                         7.4667, 8)
# model.alpha02 <- lm(b.opt.alpha02.grid ~ rho.alpha02.grid)
# model.alpha02
# model.alpha02$coefficients
# ymin.alpha02 <- model.alpha02$coefficients[1]
# ymin.alpha02
# ymax.alpha02 <- model.alpha02$coefficients[1] +
#   model.alpha02$coefficients[2]
# ymax.alpha02
#
#
# rho.alpha05.grid <- c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# b.opt.alpha05.grid <- c(4.5, 4.76667, 5.4, 5.85, 6.06667, 6.716667, 7.2333, 7.5)
# model.alpha05 <- lm(b.opt.alpha05.grid ~ rho.alpha05.grid)
# model.alpha05
# model.alpha05$coefficients
# ymin.alpha05 <- model.alpha05$coefficients[1]
# ymin.alpha05
# ymax.alpha05 <- model.alpha05$coefficients[1] +
#   model.alpha05$coefficients[2]
# ymax.alpha05
#
# rho.alpha10.grid <- c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# b.opt.alpha10.grid <- c(4.5, 4.76667, 5.2, 5.6333, 5.85, 6.28333,
#                         6.5333, 7)
# model.alpha10 <- lm(b.opt.alpha10.grid ~ rho.alpha10.grid)
# model.alpha10
# model.alpha10$coefficients
# ymin.alpha10 <- model.alpha10$coefficients[1]
# ymin.alpha10
# ymax.alpha10 <- model.alpha10$coefficients[1] +
#   model.alpha10$coefficients[2]
# ymax.alpha10
#
# ymin <- min(ymin.alpha02, ymin.alpha05, ymin.alpha10)
# ymax <- max(ymax.alpha02, ymax.alpha05, ymax.alpha10)
#
# plot(rho.alpha05.grid, b.opt.alpha05.grid, las=1, xlab="rho",
#      xlim = c(0,1),
#      xaxs = "i", ylim = c(ymin, ymax),
#      yaxs = "i", ylab="optimum value of d")
# points(rho.alpha01.grid, b.opt.alpha01.grid, col="green")
# points(rho.alpha02.grid, b.opt.alpha02.grid, col="red")
# points(rho.alpha10.grid, b.opt.alpha10.grid, col="blue")
# mtext(paste("alpha: 0.01 (green), 0.02 (red), 0.05 (black), 0.1 (blue)"), cex=0.7)
# abline(lm(b.opt.alpha02.grid ~ rho.alpha02.grid), col="red")
# abline(lm(b.opt.alpha05.grid ~ rho.alpha05.grid))
# abline(lm(b.opt.alpha10.grid ~ rho.alpha10.grid), col="blue")
#
# (4.105632 + 3.971 + 4.084) / 3
# # 4.053544
#
# # Just suppose that the intercept term is 4.1 for all
# # alpha
#
# b.opt.alpha01.minus4point1.grid <- b.opt.alpha01.grid - 4.1
# model.alpha01.minus4point1 <-
#   lm(b.opt.alpha01.minus4point1.grid ~ rho.alpha01.grid - 1)
# model.alpha01.minus4point1
#
# b.opt.alpha02.minus4point1.grid <- b.opt.alpha02.grid - 4.1
# model.alpha02.minus4point1 <-
#   lm(b.opt.alpha02.minus4point1.grid ~ rho.alpha02.grid - 1)
# model.alpha02.minus4point1
#
# b.opt.alpha05.minus4point1.grid <- b.opt.alpha05.grid - 4.1
# model.alpha05.minus4point1 <-
#   lm(b.opt.alpha05.minus4point1.grid ~ rho.alpha05.grid - 1)
# model.alpha05.minus4point1
#
# b.opt.alpha10.minus4point1.grid <- b.opt.alpha10.grid - 4.1
# model.alpha10.minus4point1 <-
#   lm(b.opt.alpha10.minus4point1.grid ~ rho.alpha10.grid - 1)
# model.alpha10.minus4point1
#
# alphavec <- c(0.01, 0.02, 0.05, 0.1)
# slopevec <- c(4.575, 4.525, 3.676, 3.085)
# plot(alphavec, slopevec,
#      xlab = "alpha", ylab="slope", ylim = c(min(slopevec)-0.2,
#                                           max(slopevec)+0.2),
#      xlim = c(0,0.1), xaxs="i", las=1)
# abline(lm(slopevec ~ alphavec))
# model <- lm(slopevec ~ alphavec)
# model
# model$coefficients
#
# tmp.u <- model$coefficients[1] + 0.01 * model$coefficients[2]
# tmp.u
# abline(h=tmp.u, col="red")
# tmp.l <- model$coefficients[1] + 0.1 * model$coefficients[2]
# tmp.l
# abline(h=tmp.l, col="red")

#*************************
# Check that the "goldilocks" value
# of d is correctly given by the formulas
# that I have adopted.
#
# rho.alpha01.grid <- c(0.02, 0.1, 0.2, 0.4, 0.5)
# b.opt.alpha01.grid <- c(4.5, 4.58333, 4.7667, 5.85, 6.5333)
# model.alpha01 <- lm(b.opt.alpha01.grid ~ rho.alpha01.grid)
# model.alpha01
# model.alpha01$coefficients
# ymin.alpha01 <- model.alpha01$coefficients[1]
# ymin.alpha01
# ymax.alpha01 <- model.alpha01$coefficients[1] +
#   model.alpha01$coefficients[2]
# ymax.alpha01
#
# rho.alpha02.grid <- c(0.02, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8)
# b.opt.alpha02.grid <- c(4.333, 4.76667, 4.7667, 5.85, 6.06667, 6.5333,
#                         7.4667, 8)
# model.alpha02 <- lm(b.opt.alpha02.grid ~ rho.alpha02.grid)
# model.alpha02
# model.alpha02$coefficients
# ymin.alpha02 <- model.alpha02$coefficients[1]
# ymin.alpha02
# ymax.alpha02 <- model.alpha02$coefficients[1] +
#   model.alpha02$coefficients[2]
# ymax.alpha02
#
#
# rho.alpha05.grid <- c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# b.opt.alpha05.grid <- c(4.5, 4.76667, 5.4, 5.85, 6.06667, 6.716667, 7.2333, 7.5)
# model.alpha05 <- lm(b.opt.alpha05.grid ~ rho.alpha05.grid)
# model.alpha05
# model.alpha05$coefficients
# ymin.alpha05 <- model.alpha05$coefficients[1]
# ymin.alpha05
# ymax.alpha05 <- model.alpha05$coefficients[1] +
#   model.alpha05$coefficients[2]
# ymax.alpha05
#
# rho.alpha10.grid <- c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# b.opt.alpha10.grid <- c(4.5, 4.76667, 5.2, 5.6333, 5.85, 6.28333,
#                         6.5333, 7)
# model.alpha10 <- lm(b.opt.alpha10.grid ~ rho.alpha10.grid)
# model.alpha10
# model.alpha10$coefficients
# ymin.alpha10 <- model.alpha10$coefficients[1]
# ymin.alpha10
# ymax.alpha10 <- model.alpha10$coefficients[1] +
#   model.alpha10$coefficients[2]
# ymax.alpha10
#
# ymin <- min(ymin.alpha02, ymin.alpha05, ymin.alpha10)
# ymax <- max(ymax.alpha02, ymax.alpha05, ymax.alpha10)
#
# plot(rho.alpha05.grid, b.opt.alpha05.grid, las=1, xlab="rho",
#      xlim = c(0,1),
#      xaxs = "i", ylim = c(ymin, ymax),
#      yaxs = "i", ylab="goldilocks value of d")
# points(rho.alpha01.grid, b.opt.alpha01.grid, col="green")
# points(rho.alpha02.grid, b.opt.alpha02.grid, col="red")
# points(rho.alpha10.grid, b.opt.alpha10.grid, col="blue")
# mtext(paste("alpha: 0.01 (green), 0.02 (red)",
#             "0.05 (black), 0.1 (blue)"), cex=0.7)
#
# rho.vec <- c(0,1)
#
# alpha <- 0.01
# d.goldilocks <- d_goldilocks(alpha, rho.vec)
# d.goldilocks
# lines(rho.vec, d.goldilocks, col="green")
#
# alpha <- 0.02
# d.goldilocks <- d_goldilocks(alpha, rho.vec)
# lines(rho.vec, d.goldilocks, col="red")
#
# alpha <- 0.05
# d.goldilocks <- d_goldilocks(alpha, rho.vec)
# lines(rho.vec, d.goldilocks)
#
# alpha <- 0.1
# d.goldilocks <- d_goldilocks(alpha, rho.vec)
# lines(rho.vec, d.goldilocks, col="blue")
#
# help(bsciuupi)
#
#

#********************
# Using system.time for getting
# the computation time
#
# sleep_for_10secs <- function() { Sys.sleep(10) }
# test.time1 <- system.time(
# sleep_for_10secs()
# )
# test.time1
# named.user.time <- test.time1[1]
# user.time <- unname(named.user.time)
# user.time
#
# natural <- 1
# alpha <- 0.05
# rho <- 0.7
# system.time(
# temp <- bs_ciuupi(alpha, rho, natural)
# )
# temp
#
# bsvec <- temp$bsvec
# d <- temp$d
# n.ints <- temp$n.ints
# comp.time <- temp$comp.time
# dev.new(width = 21, height = 30 , unit = cm,
#         noRStudioGD = TRUE)
# plots_b_s(bsvec, d, n.ints, alpha, natural, comp.time)
# dev.new(width = 21, height = 30 , unit = cm,
#         noRStudioGD = TRUE)
# plots_cp_sel(bsvec, d, n.ints, alpha, natural)


#************************
# Testing of the module
# bs_ciuupi.R

# natural <- 1
#
# alpha <- 0.1
#
# rho.vec <-
# c(0.01, 0.02, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7,
#   0.8, 0.9, 0.95, 0.98, 0.99)
# len.rho.vec <- length(rho.vec)
#
# for (k in c(1: len.rho.vec)){
#   rho <- rho.vec[k]
#   cat("rho =", rho, "\n")
#   temp <- bs_ciuupi(alpha, rho, natural)
#   cat("\n")
#   bsvec <- temp$bsvec
#   d <- temp$d
#   n.ints <- temp$n.ints
#   comp.time <- temp$comp.time
#   dev.new(width = 21, height = 30 , unit = cm,
#         noRStudioGD = TRUE)
#   plots_b_s(bsvec, d, n.ints, alpha, natural, comp.time)
#   dev.new(width = 21, height = 30 , unit = cm,
#           noRStudioGD = TRUE)
#   plots_cp_sel(bsvec, d, n.ints, alpha, natural)
# }
#

#**************************
# Test plot_b.R

# natural <- 1
# alpha <- 0.05
# rho <- -1/sqrt(2)
# system.time(
# bs.list <- bs_ciuupi(alpha, rho, natural)
# )
#
# saveRDS(bs.list, "bs.list.example")
#
# plot_squared_sel(bs.list)

#*****************************
# Test ciuupi_from_data

# Specify a, c and design matrix X
# a <- c(0, 2, 0, -2)
# c <- c(0, 0, 0, 1)
# x1 <- c(-1, 1, -1, 1)
# x2 <- c(-1, -1, 1, 1)
# X <- cbind(rep(1, 4), x1, x2, x1*x2)
# #
# # Compute rho
# # The exact value of rho is -1/sqrt(2)
# # The computed value of rho is -0.7071068
# rho <- acX_to_rho(a, c, X)
# is.array(rho)
#
# # Specify alpha
# alpha <- 0.05
#
# #
# # Execution of the command
# #      bs.list.example <- bs_ciuupi(alpha, rho)
# # takes roughly 5 minutes.
# # The result of executing this command
# # has been stored in the RDA file
# #       bs.list.example
# bs.list.example
#
# # Specify t and and the observed response y
# t <- 0
# y <- c(87.2, 88.4, 86.7, 89.2)
#
# # Specify sig, the known value of the standard
# # deviation of the random error
# sig <- 0.8
#
# # Find the observed value of the CIUUPI
# ciuupi_observed_value(a, c, X, alpha, bs.list.example, t, y, sig=sig)
# #
# # For comparison, find the standard 1-alpha
# # confidence interval
# ci_standard(a, X, y, alpha, sig=sig)
#

# plot_b(bs.list.example)
# d <- bs.list.example$d
# n.ints <- bs.list.example$n.ints
# graphics::mtext(paste("d=", signif(d, digits=6), ", q=n.ints=",
#                       n.ints), cex = 0.9)
#
# plot_s(bs.list.example)
# d <- bs.list.example$d
# n.ints <- bs.list.example$n.ints
# graphics::mtext(paste("d=", signif(d, digits=6), ", q=n.ints=",
#                       n.ints), cex = 0.9)
#
# alpha <- 0.05
# rho <- - 1 / sqrt(2)
# system.time(bs.list.example.clamped <- bs_ciuupi(alpha, rho, natural = 0))

# a <- c(0, 2, 0, -2)
# c <- c(0, 0, 0, 1)
# x1 <- c(-1, 1, -1, 1)
# x2 <- c(-1, -1, 1, 1)
# X <- cbind(rep(1, 4), x1, x2, x1*x2)
# rho <- acX_to_rho(a, c, X)

# x <- seq(0, 8, by = 0.5)
# alpha <- bs.list.example$alpha
# natural <- bs.list.example$natural
# d <- bs.list.example$d
# n.ints <- bs.list.example$n.ints
# bsvec <- bs.list.example$bsvec
# bs <- bsspline(x, bsvec, alpha, d, n.ints, natural)

# y <- c(87.2, 88.4, 86.7, 89.2)
# x1 <- c(-1, 1, -1, 1)
# x2 <- c(-1, -1, 1, 1)
# X <- cbind(rep(1, 4), x1, x2, x1*x2)
# a <- c(0, 2, 0, -2)
# ci.standard <- ci_standard(a, X, y, 0.05, sig = 0.8)

# a <- c(0, 2, 0, -2)
# c <- c(0, 0, 0, 1)
# x1 <- c(-1, 1, -1, 1)
# x2 <- c(-1, -1, 1, 1)
# X <- cbind(rep(1, 4), x1, x2, x1*x2)
# alpha <- 0.05
# t <- 0
# y <- c(87.2, 88.4, 86.7, 89.2)
# sig <- 0.8
# ciuupi.observed <- ciuupi_observed_value(a, c, X, alpha, bs.list.example,
# t, y, sig=sig)

# gam <- seq(0, 10, by = 0.1)
# n.nodes <- 10
# cp <- cpciuupi(gam, n.nodes, bs.list.example)
# cp

# y <- c(87.2, 88.4, 86.7, 89.2)
# x1 <- c(-1, 1, -1, 1)
# x2 <- c(-1, -1, 1, 1)
# X <- cbind(rep(1, 4), x1, x2, x1*x2)
# a <- c(0, 2, 0, -2)
# ci.standard <-ci_standard(a, X, y, 0.05, sig = 0.8)
# print(ci.standard)

# Test revised version of acX_to_rho.R
# a <- c(0, 2, 0, -2)
# c <- c(0, 0, 0, 1)
# x1 <- c(-1, 1, -1, 1)
# x2 <- c(-1, -1, 1, 1)
# X <- cbind(rep(1, 4), x1, x2, x1*x2)
# acX_to_rho(a, c, X) + (1/sqrt(2))

# alpha <- 0.05
# qnorm(1-alpha/2)
# qt(1-alpha/2, 30)

# y <- c(87.2, 88.4, 86.7, 89.2)
# x1 <- c(-1, 1, -1, 1)
# x2 <- c(-1, -1, 1, 1)
# X <- cbind(rep(1, 4), x1, x2, x1*x2)
# a <- c(0, 2, 0, -2)
# ci_standard(a, X, y, 0.05, sig = 0.8)
#
# a <- c(0, 2, 0, -2)
# c <- c(0, 0, 0, 1)
# x1 <- c(-1, 1, -1, 1)
# x2 <- c(-1, -1, 1, 1)
# X <- cbind(rep(1, 4), x1, x2, x1*x2)
# alpha <- 0.05
# t <- 0
# y <- c(87.2, 88.4, 86.7, 89.2)
# sig <- 0.8
# ciuupi_observed_value(a, c, X, alpha, bs.list.example, t, y, sig=sig)


