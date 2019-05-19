# simulate data and parameters
sim_X <- function(n, p) matrix(rnorm(n*p), n, p)
sim_w <- function(n) abs(rnorm(n))
sim_beta <- function(p) rnorm(p)
