context("lvlm_fitLS")

source("hlm-testfunctions.R")

test_that("lvlm_fitLS is same in C++ and R", {
  ntest <- 20
  for(ii in 1:ntest) {
    n <- sample(100:200, 1)
    p <- sample(3:5, 1)
    Z <- sim_Z(n, p)
    w <- sim_w(n)
    gamma <- sim_gamma(p)
    sigma <- exp(.5 * c(Z %*% gamma))/sqrt(w)
    y <- rnorm(n, sd = sigma)
    wy2 <- w*y^2
    gamma_r <- qr.solve(Z, log(wy2) - (digamma(1/2) + log(2)))
    gamma_cpp <- lvlm_fit(y2 = wy2, Z = Z, method = "LS")
    expect_equal(gamma_r, gamma_cpp)
  }
})
