context("lvlm_fit")

source("hlm-testfunctions.R")

test_that("lvlm_fitFS is same in C++ and R", {
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
    suppressWarnings({
      if(initLS) {
        hfit <- lvlm_fitFS_R(y2 = wy2, Z = Z,
                             maxit = maxit, epsilon = epsilon)
        gamma_cpp <- lvlm_fit(y2 = wy2, Z = Z, method = "Fisher",
                              maxit = maxit, epsilon = epsilon)
      } else {
        hfit <- lvlm_fitFS_R(y2 = wy2, Z = Z, gamma0 = gamma0,
                             maxit = maxit, epsilon = epsilon)
        gamma_cpp <- lvlm_fit(y2 = wy2, Z = Z, method = "Fisher",
                              gamma0 = gamma0, maxit = maxit, epsilon = epsilon)
      }
    })
    gamma_r <- coef(hfit)
    expect_equal(gamma_r, gamma_cpp)
  }
})

test_that("lvlm Fisher Scoring and IRLS are sufficiently close", {
  ntest <- 20
  for(ii in 1:ntest) {
    n <- sample(100:200, 1)
    p <- sample(3:5, 1)
    Z <- sim_Z(n, p)
    w <- sim_w(n)
    ## initLS <- sample(c(TRUE, FALSE), 1)
    gamma <- sim_gamma(p)
    sigma <- exp(.5 * c(Z %*% gamma))/sqrt(w)
    y <- rnorm(n, sd = sigma)
    wy2 <- w*y^2
    ## maxit <- sample(1:10, 1)
    ## epsilon <- 10^runif(1, -10, 0)
    gamma0 <- sim_gamma(p)
    suppressWarnings({
      fit_fs <- lvlm_fitFS_R(y2 = wy2, Z = Z, gamma0 = gamma0,
                             maxit = 100, epsilon = 1e-10)
      fit_irls <- lvlm_fitIRLS_R(y2 = wy2, Z = Z, gamma0 = gamma0,
                                 maxit = 100, epsilon = 1e-10)
    })
    gamma_fs <- coef(fit_fs)
    gamma_irls <- coef(fit_irls)
    expect_equal(gamma_fs, gamma_irls, tolerance = 1e-4)
  }
})

test_that("lvlm_fitIRLS is same in C++ and R", {
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
    suppressWarnings({
      if(initLS) {
        hfit <- lvlm_fitIRLS_R(y2 = wy2, Z = Z,
                               maxit = maxit, epsilon = epsilon)
        gamma_cpp <- lvlm_fit(y2 = wy2, Z = Z, method = "IRLS",
                              maxit = maxit, epsilon = epsilon)
      } else {
        hfit <- lvlm_fitIRLS_R(y2 = wy2, Z = Z, gamma0 = gamma0,
                               maxit = maxit, epsilon = epsilon)
        gamma_cpp <- lvlm_fit(y2 = wy2, Z = Z, method = "IRLS",
                              gamma0 = gamma0, maxit = maxit, epsilon = epsilon)
      }
    })
    gamma_r <- coef(hfit)
    expect_equal(gamma_r, gamma_cpp)
  }
})
