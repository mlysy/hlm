# simulate data and parameters
sim_X <- function(n, p) matrix(rnorm(n*p), n, p)
sim_w <- function(n) abs(rnorm(n))
sim_beta <- function(p) rnorm(p)
sim_Z <- function(n, p) matrix(rnorm(n*p, sd = .1), n, p)
sim_w <- function(n) abs(rnorm(n, mean = 1, sd = .2))
sim_gamma <- function(p) rnorm(p)

# relative error
relerr <- function(x_new, x_old) {
  abs(x_new - x_old)/(.1 + abs(x_new))
}


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
lvlm_fitR <- function(y2, Z, gamma0, maxit = 25, epsilon = 1e-8) {
  # constants and "memory allocation"
  rho <- digamma(1/2) + log(2) # mean of log-chisq(1)
  ZtZ <- crossprod(Z) # Z'Z
  # deviance function
  loglik <- function(gam) {
    zg <- c(Z %*% gam)
    -.5 * sum(y2/exp(zg) + zg)
  }
  # initial value
  if(missing(gamma0)) {
    gamma0 <- c(solve(ZtZ, crossprod(Z, log(y2) - rho)))
  }
  gamma <- gamma0
  ll_old <- loglik(gamma)
  # subsequent steps
  for(ii in 1:maxit) {
    # score function (up to factor of .5)
    sc <- colSums((y2/exp(c(Z %*% gamma)) - 1) * Z)
    # fisher scoring step
    gamma <- gamma + c(solve(ZtZ, sc))
    # check tolerance
    ll_new <- loglik(gamma)
    tol <- relerr(ll_new, ll_old)
    if(tol < epsilon) break else ll_old <- ll_new
  }
  if(ii == maxit && tol > epsilon) warning("Fisher scoring did not converge.")
  return(list(coefficients = gamma, iter = ii))
}

hlm_loglik <- function(beta, gamma, y, X, Z) {
    mu <- c(X %*% beta)
    sig <- exp(.5 * Z %*% gamma)
    sum(dnorm(y, mean = mu, sd = sig, log = TRUE))
}

hlm_fitR <- function(y, X, Z, beta0, gamma0, maxit = 25, epsilon = 1e-8) {
  C <- length(y)/2 * log(2*pi)
  loglik <- function(beta, gamma) {
    hlm_loglik(beta = beta, gamma = gamma, y = y, X = X, Z = Z) + C
  }
  beta <- beta0
  gamma <- gamma0
  ll_old <- loglik(beta, gamma)
  for(ii in 1:maxit) {
    # update beta
    w <- exp(-c(Z %*% gamma))
    beta <- coef(lm.wfit(x = X, y = y, w = w))
    # update gamma
    y2 <- (y - c(X %*% beta))^2
    gamma <- coef(lvlm_fitR(y2 = y2, Z = Z))
    ll_new <- loglik(beta, gamma)
    tol <- relerr(ll_new, ll_old)
    ## message("ll_old = ", ll_old, ", ll_new = ", ll_new)
    ## message("error[",ii,"] = ", tol)
    if(tol < epsilon) break else ll_old <- ll_new
  }
  list(beta = setNames(beta, NULL), gamma = setNames(gamma, NULL),
       loglik = ll_new, niter = ii, error = tol)
}
