source("hlm-testfunctions.R")

context("chlm_hess")

test_that("chlm_hess and numDeriv::hessian are the same.", {
  ntest <- 20
  for(ii in 1:ntest) {
    n <- sample(10:20, 1)
    p <- sample(1:4, 1)
    q <- sample(1:4, 1)
    X <- sim_X(n,p)
    beta <- sim_beta(p)
    Z <- sim_Z(n,q)
    gamma <- sim_gamma(q)
    mu <- c(X %*% beta)
    sig <- exp(c(Z %*% gamma)/2)
    y <- rnorm(n, mu, sig)
    delta <- sample(c(TRUE, FALSE), n, replace = TRUE)
    # calculation
    # analytic
    H_an <- chlm_hess(beta, gamma, y, delta, X, Z)
    # numerical
    H_nu <- numDeriv::hessian(fun = function(theta) {
      chlm_loglik(beta = theta[1:p], gamma = theta[p+1:q], y, delta, X, Z)
    }, x = c(beta, gamma))
    expect_equal(H_an, H_nu, tol = 1e-6)
  }
})
