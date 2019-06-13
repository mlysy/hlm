# simulate data and parameters
sim_X <- function(n, p) matrix(rnorm(n*p), n, p)
sim_beta <- function(p) rnorm(p)
sim_Z <- function(n, p) matrix(rnorm(n*p, sd = .1), n, p)
sim_w <- function(n) abs(rnorm(n, mean = 1, sd = .2))
sim_gamma <- function(p) rnorm(p)

# relative error
relerr <- function(x_new, x_old) {
  abs(x_new - x_old)/(.1 + abs(x_new))
}

# max of min of abs and rel error
max_xdiff <- function(x) {
  xdiff <- abs(diff(x))
  max(pmin(xdiff[,1], xdiff[,2]))
}


lvlm_loglik <- function(gamma, y2, Z) {
  Zg <- c(Z %*% gamma)
  -.5 * sum(y2/exp(Zg) + Zg)
}

lvlm_fitIRLS_R <- function(y2, Z, gamma0, maxit = 25, epsilon = 1e-8) {
  rho <- digamma(1/2) + log(2) # mean of log-chisq(1)
  loglik <- function(gamma) lvlm_loglik(gamma, y2, Z)
  logY2 <- log(y2)
  # initial value
  if(missing(gamma0)) {
    gamma0 <- c(solve(crossprod(Z), crossprod(Z, log(y2) - rho)))
  }
  gamma <- gamma0
  ll_old <- loglik(gamma)
  for(ii in 1:maxit) {
    mu <- logY2 - c(Z %*% gamma)
    wgt <- (exp(mu) - 1)/mu
    Zw <- Z * wgt
    gamma <- c(solve(crossprod(Zw, Z), crossprod(Zw, logY2)))
    # check tolerance
    ll_new <- loglik(gamma)
    tol <- relerr(ll_new, ll_old)
    if(tol < epsilon) break else ll_old <- ll_new
  }
  if(ii == maxit && tol > epsilon) warning("IRLS did not converge.")
  return(list(coefficients = gamma, iter = ii, loglik = ll_new))
}

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
lvlm_fitFS_R <- function(y2, Z, gamma0, maxit = 25, epsilon = 1e-8) {
  # constants and "memory allocation"
  rho <- digamma(1/2) + log(2) # mean of log-chisq(1)
  ZtZ <- crossprod(Z) # Z'Z
  # deviance function
  loglik <- function(gamma) {
    lvlm_loglik(gamma, y2, Z)
    ## zg <- c(Z %*% gam)
    ## -.5 * sum(y2/exp(zg) + zg)
  }
  # initial value
  if(missing(gamma0)) {
    gamma0 <- c(solve(ZtZ, crossprod(Z, log(y2) - rho)))
  }
  gamma <- gamma0
  ll_old <- loglik(gamma)
  # subsequent steps
  for(ii in 1:maxit) {
    # score function (up to factor of .5)
    sc <- colSums((y2/exp(c(Z %*% gamma)) - 1) * Z)
    # fisher scoring step
    gamma <- gamma + c(solve(ZtZ, sc))
    # check tolerance
    ll_new <- loglik(gamma)
    tol <- relerr(ll_new, ll_old)
    if(tol < epsilon) break else ll_old <- ll_new
  }
  if(ii == maxit && tol > epsilon) warning("Fisher scoring did not converge.")
  return(list(coefficients = gamma, iter = ii, loglik = ll_new))
}

hlm_loglik <- function(beta, gamma, y, X, Z) {
    mu <- c(X %*% beta)
    sig <- exp(.5 * Z %*% gamma)
    sum(dnorm(y, mean = mu, sd = sig, log = TRUE))
}

hlm_fit_R <- function(y, X, Z, beta0, gamma0,
                      maxit = 25, epsilon = 1e-8,
                      method = c("IRLS", "Fisher")) {
  method <- match.arg(method)
  method <- switch(method, Fisher = 0, IRLS = 1)
  C <- length(y)/2 * log(2*pi)
  loglik <- function(beta, gamma) {
    hlm_loglik(beta = beta, gamma = gamma, y = y, X = X, Z = Z) + C
  }
  beta <- beta0
  gamma <- gamma0
  ll_old <- loglik(beta, gamma)
  for(ii in 1:maxit) {
    # update beta
    w <- exp(-c(Z %*% gamma))
    beta <- coef(lm.wfit(x = X, y = y, w = w))
    # update gamma
    y2 <- (y - c(X %*% beta))^2
    if(method == 0) {
      gamma <- coef(lvlm_fitFS_R(y2 = y2, Z = Z))
    } else if(method == 1) {
      gamma <- coef(lvlm_fitIRLS_R(y2 = y2, Z = Z))
    }
    ll_new <- loglik(beta, gamma)
    tol <- relerr(ll_new, ll_old)
    ## message("ll_old = ", ll_old, ", ll_new = ", ll_new)
    ## message("error[",ii,"] = ", tol)
    if(tol < epsilon) break else ll_old <- ll_new
  }
  list(beta = setNames(beta, NULL), gamma = setNames(gamma, NULL),
       loglik = ll_new, niter = ii, error = tol)
}

## chlm_loglik <- function(beta, gamma, y, delta, X, Z) {
##   mu <- c(X %*% beta)
##   sig <- exp(.5 * c(Z %*% gamma))
##   ans <- rep(NA, length(y))
##   ans[delta] <- dnorm(y[delta], mu[delta], sig[delta], log = TRUE)
##   ans[!delta] <- pnorm(y[!delta], mu[!delta], sig[!delta],
##                        lower.tail = FALSE, log.p = TRUE)
##   sum(ans)
## }

chlm_fit_R <- function(y, delta, X, Z,
                       beta0, gamma0,
                       maxit = 25, epsilon = 1e-8, splitE = FALSE) {
  # precomputations
  ## C <- length(y)/2 * log(2*pi)
  N <- length(y)
  Xc <- X[!delta,,drop=FALSE]
  yc <- y[!delta]
  R <- y
  S <- y^2
  # helper functions
  loglik <- function(beta, gamma) {
    chlm_loglik(beta, gamma, y, delta, X, Z)
  }
  # E[Z|Z>a]
  etnorm <- function(a) {
    return(dnorm(a)/pnorm(-a))
  }
  # E[Z^2|Z>a]
  e2tnorm <- function(a) {
    return(1+a*dnorm(a)/pnorm(-a))
  }
  # E-step
  Estep <- function(beta, sigmac) {
    muc <- c(Xc %*% beta)
    zc <- (yc - muc)/sigmac
    fsigc <- sigmac * etnorm(zc)
    Rc <- fsigc + muc
    ## Sc1 <- sigmac^2 * e2tnorm(zc) + 2*muc*fsigc + muc^2
    Sc <- sigmac * (sigmac + zc*fsigc) + muc * (2*Rc - muc)
    ## list(Rc = Rc, Uc = Uc)
    R[!delta] <<- Rc
    S[!delta] <<- Sc
    invisible(NULL)
  }
  # fitting functions
  .lm_fit <- function(y, X) {
    ## coef(lm.fit(x = X, y = y))
    hlm:::lm_fit(y = y, X = X)
  }
  .wlm_fit <- function(y, X, w) {
    ## coef(lm.wfit(x = X, y = y, w = w))
    hlm:::wlm_fit(y = y, X = X, w = w)
  }
  .lvlm_fit <- function(y2, Z, gamma0) {
    ## coef(glm.fit(x = Z, y = y2, family = Gamma(link = "log")))
    if(missing(gamma0)) gamma0 <- hlm:::lvlm_fitLS(logY2 = log(y2), Z = Z)
    hlm:::lvlm_fitIRLS(y2 = y2, Z = Z, gamma0 = gamma0)
  }
  # initialize algorithm
  if(missing(beta0)) beta0 <- .lm_fit(y = y, X = X)
  if(missing(gamma0)) gamma0 <- .lvlm_fit(y2 = c(y-X%*%beta0)^2, Z = Z)
  ## if(missing(gamma0)) gamma0 <- coef(lvlm_fitR(y2 = c(y-X%*%beta0)^2, Z = Z))
  beta <- beta0
  gamma <- gamma0
  ll_old <- loglik(beta, gamma)
  for(ii in 1:maxit) {
    # E-step
    Zg <- c(Z %*% gamma)
    sigmac <- exp(.5 * Zg[!delta])
    Estep(beta, sigmac)
    # M-step: beta
    w <- exp(-Zg)
    ## beta <- coef(lm.wfit(x = X, y = R, w = w))
    beta <- .wlm_fit(y = R, X = X, w = w)
    if(splitE) {
      # recompute E-step
      Estep(beta, sigmac)
    }
    # M-step: gamma
    Xb <- c(X %*% beta)
    ## U <- S - 2*R*Xb + Xb^2
    ## gamma <- coef(lvlm_fitR(y2 = U, Z = Z))
    U <- S + Xb * (Xb - 2*R)
    gamma <- .lvlm_fit(y2 = U, Z = Z)
    ll_new <- loglik(beta, gamma)
    tol <- relerr(ll_new, ll_old)
    ## message("ll_old = ", ll_old, ", ll_new = ", ll_new)
    ## message("error[",ii,"] = ", tol)
    if(tol < epsilon) break else ll_old <- ll_new
  }
  list(beta = setNames(beta, NULL), gamma = setNames(gamma, NULL),
       loglik = ll_new, niter = ii, error = tol)
}


# by and large, the original version (for comparison purposes)
chlm_fit_old <- function(y, delta, X, W, nreps = 1e3, tol = 1e-5,
                         multi.cycle = FALSE) {
  # Helper functions
  # E[Z|Z>a]
  etnorm <- function(a) {
    return(dnorm(a)/pnorm(-a))
  }
  # E[Z^2|Z>a]
  e2tnorm <- function(a) {
    return(1+a*dnorm(a)/pnorm(-a))
  }
  # relative error
  rel.err <- function(theta1, theta2) {
    abs(theta1-theta2)/abs((theta1+theta2)/2)
  }
  # loglikelihood function
  loglik <- function(theta) {
    ## t <- exp(y)
    beta <- theta[1:px]
    gamma <- theta[px+1:pw]
    chlm_loglik(beta, gamma, y, delta, X, W)
    ## mu <- c(X%*%beta)
    ## sigma <- exp(c(W%*%gamma)/2)
    ## ll    <- ifelse(delta, dnorm(y, mean = mu, sd = sigma, log = TRUE),
    ##                 pnorm(y, mean = mu, sd = sigma, lower.tail = FALSE, log.p = TRUE))
    ## return(sum(ll))
  }
  # E-step
  Estep <- function(beta, sigmac) {
    ## mut <- c(Xc %*% betat)
    ## Wg <- c(W %*% gammat)
    ## sigmat <- exp(Wg[!delta]/2)
    ## if(ii > 5 &&
    ##    ((max(mut) > 12 || min(exp(mut + 2 * sigmat) - exp(mut - 2* sigmat)) < 365/12))) {
    ##   bad.model <- TRUE
    ##   break
    ## }
    ## zt <- (yc - mut)/sigmat
    ## sft <- sigmat * f(zt)
    ## T[!delta] <- sft + mut
    ## U[!delta] <- sigmat^2 * g(zt) + 2*mut * sft + mut^2
    muc <- c(Xc %*% beta)
    zc <- (yc - muc)/sigmac
    fsigc <- sigmac * etnorm(zc)
    Rc <- fsigc + muc
    Sc <- sigmac^2 * e2tnorm(zc) + 2*muc*fsigc + muc^2
    ## Sc <- sigmac * (sigmac + zc*fsigc) + muc * (2*tc - muc)
    ## list(Rc = Rc, Uc = Uc)
    T[!delta] <<- Rc
    U[!delta] <<- Sc
    invisible(NULL)
  }
  # fitting functions
  .lm_fit <- function(y, X) {
    ## coef(lm.fit(x = X, y = y))
    hlm:::lm_fit(y = y, X = X)
  }
  .wlm_fit <- function(y, X, w) {
    ## coef(lm.wfit(x = X, y = y, w = w))
    hlm:::wlm_fit(y = y, X = X, w = w)
  }
  .lvlm_fit <- function(y2, Z, gamma0) {
    ## coef(glm.fit(x = Z, y = y2, family = Gamma(link = "log")))
    if(missing(gamma0)) gamma0 <- hlm:::lvlm_fitLS(logY2 = log(y2), Z = Z)
    hlm:::lvlm_fit(y2 = y2, Z = Z, gamma0 = gamma0)
  }
  # obsolete?
  loglik.orig <- function(theta) {
    t <- exp(y)
    beta <- theta[1:px]
    gamma <- theta[px+1:pw]
    mu <- c(X%*%beta)
    sigma <- exp(c(W%*%gamma)/2)
    ll.orig <- ifelse(delta, dlnorm(t, meanlog = mu, sdlog = sigma, log = TRUE),
                      plnorm(t, meanlog = mu, sdlog = sigma, lower.tail = FALSE, log.p = TRUE))
    return(sum(ll.orig))
  }
  # some constants
  X <- as.matrix(X)
  W <- as.matrix(W)
  # covariate dimensions
  px <- ncol(X)
  pw <- ncol(W)
  # covariate names
  if(is.null(colnames(X))) colnames(X) <- 1:px
  if(is.null(colnames(W))) colnames(W) <- 1:pw
  XW.names <- c(paste0("X", colnames(X)), paste0("W", colnames(W)))
  # covariates for censored observations
  Xc <- as.matrix(X[!delta,])
  Wc <- as.matrix(W[!delta,])
  yc <- y[!delta]
  # T.tilde and U.tilde.  some of these values will never be updated
  T <- y
  U <- y^2
  # initialize the model parameters
  betat <- .lm_fit(y = y, X = X)
  gammat <- .lvlm_fit(y2 = (y-X%*%betat)^2, Z = W)
  ll_old <- loglik(c(betat, gammat))
  ## betat <- coef(lm.fit(x = X, y = y))
  ## gammat <- coef(glm.fit(x = W, y = (y-X%*%betat)^2, family = Gamma(link = "log")))
  # main loop
  bad.model <- FALSE
  for(ii in 1:nreps) {
    # E-step
    Wg <- c(W %*% gammat)
    sigmat <- exp(Wg[!delta]/2)
    Estep(betat, sigmat)
    ## mut <- c(Xc %*% betat)
    ## Wg <- c(W %*% gammat)
    ## sigmat <- exp(Wg[!delta]/2)
    ## if(ii > 5 &&
    ##    ((max(mut) > 12 || min(exp(mut + 2 * sigmat) - exp(mut - 2* sigmat)) < 365/12))) {
    ##   bad.model <- TRUE
    ##   break
    ## }
    ## zt <- (yc - mut)/sigmat
    ## sft <- sigmat * f(zt)
    ## T[!delta] <- sft + mut
    ## U[!delta] <- sigmat^2 * g(zt) + 2*mut * sft + mut^2
    # M-step: beta
    beta <- .wlm_fit(y = T, X = X, w = exp(-Wg))
    ## beta <- coef(lm.wfit(x = X, y = T, w = exp(-Wg)))
    if(multi.cycle) {
      # recompute E-step
      Estep(beta, sigmat)
      ## mut <- c(Xc %*% beta)
      ## zt <- (yc - mut)/sigmat
      ## sft <- sigmat * f(zt)
      ## T[!delta] <- sft + mut
      ## U[!delta] <- sigmat^2 * g(zt) + 2*mut * sft + mut^2
    }
    # M-step: gamma
    Xb <- c(X %*% beta)
    R <- U - 2*T*Xb + Xb^2
    gamma <- .lvlm_fit(y2 = R, Z = W)
    ## gamma <- coef(glm.fit(x = W, y = R, family = Gamma(link = "log")))
    # update
    betat <- beta
    gammat <- gamma
    # relative error
    ## theta.rel <- rel.err(c(beta, gamma), c(betat, gammat))
    ## if(max(theta.rel) < tol) break
    ll_new <- loglik(c(beta, gamma))
    if(relerr(ll_new, ll_old) < tol) break else ll_old <- ll_new
  }
  # return some relevant values
  names(betat) <- colnames(X)
  names(gammat) <- colnames(W)
  theta.hat <- setNames(c(betat, gammat), nm = XW.names)
  ## theta.hat <- c(betat, gammat)
  ## names(theta.hat) <- XW.names
  ll.max <- loglik(theta.hat)
  ll.max.orig <- loglik.orig(theta.hat)
  if(!bad.model) {
    grad.hat <- numDeriv::grad(loglik, theta.hat)
    names(grad.hat) <- XW.names
    var.hat <- solve(-numDeriv::hessian(func=loglik,x=theta.hat))
    colnames(var.hat) <- XW.names
    rownames(var.hat) <- XW.names
    aic <- 2*length(theta.hat) - 2*ll.max
  } else {
    grad.hat <- rep(NA, length(theta.hat))
    names(grad.hat) <- XW.names
    var.hat <- matrix(NA, length(theta.hat), length(theta.hat))
    colnames(var.hat) <- XW.names
    rownames(var.hat) <- XW.names
    aic <- Inf
  }
  out <- list(coef = list(beta = betat, gamma = gammat), vcov = var.hat, score = grad.hat,
              maxll = ll.max, maxll.orig = ll.max.orig, AIC = aic, niter = ii,
              X = X, W = W, y = y, delta = delta, dist = "hlm")
  return(out)
}
