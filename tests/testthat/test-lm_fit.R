context("lm_fit")

source("hlm-testfunctions.R")

test_that("lm_fit is same in C++ and R", {
  ntest <- 20
  for(ii in 1:ntest) {
    n <- sample(10:20, 1)
    p <- sample(3:5, 1)
    X <- sim_X(n, p)
    w <- sim_w(n)
    beta <- sim_beta(p)
    y <- rnorm(n, mean = X %*% beta, sd = 1/sqrt(w))
    beta_r <- c(solve(crossprod(X, X*w), crossprod(X, y*w)))
    beta_cpp <- hlm:::lm_fit(y = y, X = X, weights = w)
    expect_equal(beta_r, beta_cpp)
  }
})
