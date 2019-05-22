require(hlm)
source("hlm-testfunctions.R")


n <- sample(1000:2000, 1)
p <- sample(5:10, 1)
q <- sample(5:10, 1)
X <- sim_X(n, p)
Z <- sim_Z(n, q)
beta <- sim_beta(p)
gamma <- sim_gamma(q)
mu <- c(X %*% beta)
sig <- exp(.5 * c(Z %*% gamma))
y <- rnorm(n, mean = mu, sd = sig)
delta <- sample(c(TRUE, FALSE),
                size = n, replace = TRUE, prob = c(.9, .1))

beta0 <- rep(0, p)
gamma0 <- rep(0, q)
maxit <- 1000
epsilon <- 1e-6
system.time({
  chfit <- chlm_fitR(y = y, delta = delta, X = X, Z = Z,
                     maxit = maxit, epsilon = epsilon)
})


optimCheck::optim_proj(xsol = c(chfit$beta, chfit$gamma),
                       fun = function(theta) {
                         chlm_loglik(theta[1:p], theta[p+1:q],
                                     y = y, delta = delta, X = X, Z = Z)
                       })

## system.time({
##   chfit2 <- chlm_fit_old(y = y, delta = delta,
##                          X = X, W = Z, nreps = maxit, tol = epsilon)
## })
## optimCheck::optim_proj(xsol = c(chfit2$coef$beta, chfit2$coef$gamma),
##                        fun = function(theta) {
##                          chlm_loglik(theta[1:p], theta[p+1:q],
##                                      y = y, delta = delta, X = X, Z = Z)
##                        })
