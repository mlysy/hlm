context("lvlm_fit")

source("hlm-testfunctions.R")

test_that("lvlm_fit is same in C++ and R", {
  ntest <- 20
  for(ii in 1:ntest) {
    n <- sample(100:200, 1)
    p <- sample(3:5, 1)
    Z <- sim_Z(n, p)
    w <- sim_w(n)
    initLS <- sample(c(TRUE, FALSE), 1)
    gamma <- sim_gamma(p)
    sigma <- exp(.5 * c(Z %*% gamma))/sqrt(w)
    y <- rnorm(n, sd = sigma)
    wy2 <- w*y^2
    maxit <- sample(1:10, 1)
    epsilon <- 10^runif(1, -10, 0)
    gamma0 <- sim_gamma(p)
    if(initLS) {
      hfit <- suppressWarnings(lvlm_fitR(y2 = wy2, Z = Z,
                                         maxit = maxit, epsilon = epsilon))
    } else {
      hfit <- suppressWarnings(lvlm_fitR(y2 = wy2, Z = Z, gamma0 = gamma0,
                                         maxit = maxit, epsilon = epsilon))
    }
    gamma_r <- coef(hfit)
    gamma_cpp <- hlm:::lvlm_fit(y2 = wy2, Z = Z, gamma0 = gamma0,
                                maxit = maxit, epsilon = epsilon,
                                initLS = initLS)
    expect_equal(gamma_r, gamma_cpp)
  }
})
