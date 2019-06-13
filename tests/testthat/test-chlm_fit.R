## require(hlm)
source("hlm-testfunctions.R")

context("chlm_fit")

test_that("chlm_fit converges to the MLE", {
  ntest <- 5
  for(ii in 1:ntest) {
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
    maxit <- 1000
    epsilon <- 1e-6
    chfit <- chlm_fit(y = y, delta = delta, X = X, Z = Z,
                      maxit = maxit, epsilon = epsilon)
    ocheck <- optimCheck::optim_proj(xsol = c(chfit$beta, chfit$gamma),
                                     fun = function(theta) {
                                       chlm_loglik(theta[1:p], theta[p+1:q],
                                                   y = y, delta = delta, X = X, Z = Z)
                                     }, plot = FALSE, npts = 20)
    expect_lt(max_xdiff(ocheck), .01)
  }
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
