# test the fisher scoring algorithm

#' @param y Vector of normal responses
#' @param Z Matrix of covariates
hlm_fit <- function(y, Z, maxit = 25, epsilon = 1e-8) {
  # constants and "memory allocation"
  rho <- digamma(1/2) + log(2) # mean of log-chisq(1)
  ZtZ <- crossprod(Z) # Z'Z
  y2 <- y^2
  # deviance function
  dev_fun <- function(gam) {
    zg <- c(Z %*% gam)
    -.5 * sum(y2/exp(zg) + zg)
  }
  # initial value
  gamma <- c(solve(ZtZ, crossprod(Z, log(y2) - rho)))
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
  if(ii == maxit && tol > epsilon) stop("Fisher scoring did not converge.")
  return(list(coefficients = gamma, iter = ii))
}

# simulate data
sim_Z <- function(n, p) matrix(rnorm(n*p), n, p)
sim_gam <- function(p) rnorm(p)
# loglikelihood
hlm_loglik <- function(gamma, y, Z) {
  zg <- c(Z %*% gamma)
  - .5 * sum(y^2/exp(zg) + zg)
}

# random data
p <- sample(20:30, 1)
n <- sample(1000:2000,1)
Z <- sim_Z(n, p)/10
gam <- sim_gam(p)
sig <- exp(.5 * c(Z %*% gam))
y <- rnorm(n, sd = sig) # normal responses
# glm fit
gfit <- glm.fit(x = Z, y = y^2, family = Gamma("log"))
# hlm fit
hfit <- hlm_fit(y, Z)
# difference in loglikelihood
hlm_loglik(coef(gfit), y, Z) - hlm_loglik(coef(hfit), y, Z)
# difference in number of iterations
gfit$iter - hfit$iter

# check timing
nreps <- 10
system.time({
  replicate(nreps, glm.fit(x = Z, y = y^2, family = Gamma("log")))
})
system.time({
  replicate(nreps, hlm_fit(y, Z))
})
