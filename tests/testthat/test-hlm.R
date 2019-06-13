## require(hlm)
source("hlm-testfunctions.R")

context("hlm_fit")

test_that("hlm formula parsing is correct", {
  ntest <- 10
  for(ii in 1:ntest) {
    # construct dataset
    n <- sample(1000:2000, 1)
    X <- as.data.frame(sim_X(n, 3)) # continuous covariates
    # discrete covariates
    X <- cbind(X,
               data.frame(R1 = sample(letters[1:3], n, replace = TRUE),
                          R2 = sample(LETTERS[10:11], n, replace = TRUE)))
    # add missing data
    for(ii in 1:length(X)) {
      X[sample(n, 100), ii] <- NA
    }
    # response
    y <- abs(rnorm(n))
    y[sample(n, 100)] <- NA
    delta <- sample(c(TRUE, FALSE), size = n, replace = TRUE, prob = c(.9, .1))
    # full dataset
    dat <- cbind(y = y, status = delta, X)
    # create covariate basis
    cov_names <- c(names(X), # main effects
                   paste0(c("sin(", "exp("), names(X)[1:3], ")"), # nonlinear effects
                   apply(combn(names(X), 2), 2, paste0, collapse = ":")) # interactions
    # formulas for hlm model
    resp <- c("log(y)", "status")[1:sample(1:2,1)]
    cov_mean <- sample(cov_names, 3)
    cov_var <- sample(cov_names, 3)
    form <- get_form(resp, cov_mean, cov_var)
    form_mean <-  get_form(resp, cov_mean)
    form_var <- get_form(resp, cov_var)
    # automated fit
    control <- chlm_control(epsilon = runif(1, 1e-5, 2e-5),
                            maxit = sample(100:200, 1),
                            nIRLS = sample(5:10, 1),
                            splitE = sample(c(TRUE, FALSE), 1))
    hfit <- hlm(form, data = dat, control = control)
    # manual fit
    Xm <- get_mm(form_mean, dat)
    Zm <- get_mm(form_var, dat)
    ym <- log(y)
    dm <- if(length(resp) == 1) rep(TRUE, n) else delta
    ikeep <- apply(!is.na(Xm), 1, all)
    ikeep <- ikeep & apply(!is.na(Zm), 1, all)
    ikeep <- ikeep & !is.na(y)
    cfit <- chlm_fit(y = ym[ikeep],
                     delta = dm[ikeep],
                     X = Xm[ikeep,], Z = Zm[ikeep,],
                     maxit = control$maxit,
                     epsilon = control$epsilon,
                     splitE = control$splitE,
                     nIRLS = control$nIRLS)
    expect_equal(setNames(cfit$beta, NULL),
                 setNames(coefficients(hfit)$beta, NULL))
    expect_equal(setNames(cfit$gamma, NULL),
                 setNames(coefficients(hfit)$gamma, NULL))
  }
})
