## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ciuupi)

## -----------------------------------------------------------------------------
x1 <- c(-1, 1, -1, 1)
x2 <- c(-1, -1, 1, 1)
X <- cbind(rep(1, 4), x1, x2, x1*x2)

## -----------------------------------------------------------------------------
a <- c(0, 2, 0, -2)

## -----------------------------------------------------------------------------
c <- c(0, 0, 0, 1)

## -----------------------------------------------------------------------------
# Compute rho
rho <- acX_to_rho(a, c, X)
print(rho)

## ----fig.align='center', fig.width=5, fig.height=4----------------------------
plot_b(bs.list.example)

## ----fig.align='center', fig.width=5, fig.height=4----------------------------
plot_s(bs.list.example)

## ----fig.align="center", fig.width=5, fig.height=4----------------------------
plot_cp(bs.list.example)

## ----fig.align="center", fig.width=5, fig.height=4----------------------------
plot_squared_sel(bs.list.example)

## -----------------------------------------------------------------------------
y <- c(87.2, 88.4, 86.7, 89.2)
sig <- 0.8

## -----------------------------------------------------------------------------
alpha <- 0.05
t <- 0
ciuupi_observed_value(a, c, X, alpha, bs.list.example, t, y, sig = sig)

## -----------------------------------------------------------------------------
alpha <- 0.05
ci_standard(a, X, y, alpha, sig)

