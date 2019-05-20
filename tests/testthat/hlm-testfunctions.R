# simulate data and parameters
sim_X <- function(n, p) matrix(rnorm(n*p), n, p)
sim_w <- function(n) abs(rnorm(n))
sim_beta <- function(p) rnorm(p)
sim_Z <- function(n, p) matrix(rnorm(n*p, sd = .1), n, p)
sim_w <- function(n) abs(rnorm(n, mean = 1, sd = .2))
sim_gamma <- function(p) rnorm(p)

#' @param y2 Vector of squared normal responses.
#' @param Z Matrix of covariates.
#' @param gamma0 Optional initial parameter vector.  If missing an OLS estimator is used.
#' @param maxit Maximum number of Fisher scoring iterations.
#' @param epsilon Error tolerance.
#' @return A list with elements:
#' \describe{
#'   \item{\code{coefficients}}{The fitted coefficient vector.}
#'   \item{\code{iter}}{The number of iterations of the Fisher scoring algorithm.}
#' }
hlm_fit <- function(y2, Z, gamma0, maxit = 25, epsilon = 1e-8) {
  # constants and "memory allocation"
  rho <- digamma(1/2) + log(2) # mean of log-chisq(1)
  ZtZ <- crossprod(Z) # Z'Z
  # deviance function
  dev_fun <- function(gam) {
    zg <- c(Z %*% gam)
    -.5 * sum(y2/exp(zg) + zg)
  }
  # initial value
  if(missing(gamma0)) {
    gamma0 <- c(solve(ZtZ, crossprod(Z, log(y2) - rho)))
  }
  gamma <- gamma0
  dev <- dev_fun(gamma)
  # subsequent steps
  for(ii in 1:maxit) {
    # score function (up to factor of .5)
    sc <- colSums((y2/exp(c(Z %*% gamma)) - 1) * Z)
    # fisher scoring step
    gamma <- gamma + c(solve(ZtZ, sc))
    # check tolerance
    dev_new <- dev_fun(gamma)
    tol <- abs(dev_new - dev)/(.1 + abs(dev_new))
    if(tol < epsilon) break else dev <- dev_new
  }
  if(ii == maxit && tol > epsilon) warning("Fisher scoring did not converge.")
  return(list(coefficients = gamma, iter = ii))
}
