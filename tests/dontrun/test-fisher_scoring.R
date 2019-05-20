# test the fisher scoring algorithm

require(hlm)

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

# simulate data
sim_Z <- function(n, p) matrix(rnorm(n*p), n, p)
sim_w <- function(n) abs(rnorm(n, mean = 1, sd = .1))
sim_gam <- function(p) rnorm(p)
# loglikelihood
hlm_loglik <- function(gamma, wy2, Z) {
  zg <- c(Z %*% gamma)
  - .5 * sum(wy2/exp(zg) + zg)
}

# random data
p <- sample(20:30, 1)
n <- sample(1000:2000,1)
Z <- sim_Z(n, p)/10
w <- sim_w(n)
gam <- sim_gam(p)
sig <- exp(.5 * c(Z %*% gam))/sqrt(w)
y <- rnorm(n, sd = sig) # normal responses
wy2 <- w*y^2
# glm fit
gfit <- glm.fit(x = Z, y = wy2, family = Gamma("log"))
# hlm fit
hfit <- hlm_fit(wy2, Z)
# difference in loglikelihood
hlm_loglik(coef(gfit), wy2, Z) - hlm_loglik(coef(hfit), wy2, Z)
# difference in number of iterations
gfit$iter - hfit$iter
# check c++ code...
maxit <- sample(1:10, 1)
epsilon <- 10^runif(1, -10, 0)
gamma0 <- sim_gam(p)
gamma_cpp <- hlm:::lvlm_fit(y2 = wy2, Z = Z, gamma0 = gamma0, maxit = maxit, epsilon = epsilon)
hfit_r <- hlm_fit(y2 = wy2, Z = Z, gamma0 = gamma0, maxit = maxit, epsilon = epsilon)
gamma_r <- coef(hfit_r)
range(gamma_cpp - gamma_r)

# OK
range(hlm:::lvlm_fitLS(logY2 = log(wy2), Z = Z) - qr.solve(Z, log(wy2) - (digamma(1/2) + log(2))))

# check mode
## optimCheck::optim_proj(xsol = coef(hfit),
##                        fun = function(gamma) {
##                          hlm_loglik(gamma, wy2 = wy2, Z = Z)
##                        })

# check timing
nreps <- 10
system.time({
  replicate(nreps, glm.fit(x = Z, y = wy2, family = Gamma("log")))
})
system.time({
  replicate(nreps, hlm_fit(wy2, Z))
})
