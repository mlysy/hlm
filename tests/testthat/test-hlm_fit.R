context("hlm_fit")

source("hlm-testfunctions.R")

test_that("hlm_fit is same in C++ and R", {
  ntest <- 20
  for(ii in 1:ntest) {
    n <- sample(100:200, 1)
    p <- sample(3:5, 1)
    q <- sample(2:4, 1)
    method <- sample(c("Fisher", "IRLS"), 1)
    X <- sim_X(n, p)
    Z <- sim_Z(n, q)
    beta <- sim_beta(p)
    gamma <- sim_gamma(q)
    mu <- c(X %*% beta)
    sig <- exp(.5 * c(Z %*% gamma))
    y <- rnorm(n, mean = mu, sd = sig)
    maxit <- sample(1:10, 1)
    epsilon <- 10^runif(1, -10, 0)
    beta0 <- sim_beta(p)
    gamma0 <- sim_gamma(q)
    suppressWarnings({
      hfit_r <- hlm_fit_R(y = y, X = X, Z = Z,
                          beta0 = beta0, gamma0 = gamma0,
                          maxit = maxit, epsilon = epsilon,
                          method = method)
    })
    hfit_cpp <- hlm_fit(y = y, X = X, Z = Z,
                        beta0 = beta0, gamma0 = gamma0,
                        maxit = maxit, epsilon = epsilon,
                        method = method)
    expect_equal(hfit_r$beta, hfit_cpp$beta)
    expect_equal(hfit_r$gamma, hfit_cpp$gamma)
  }
})

## optimCheck::optim_proj(xsol = c(hfit_r$beta, hfit_r$gamma),
##                        fun = function(theta) {
##                          hlm_loglik(theta[1:p], theta[p+1:q], y, X, Z)
##                        })
