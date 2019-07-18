#' Low-level fitting functions for the censored HLM model.
#'
#' @name chlm_fit
#' @aliases chlm_control
#'
#' @template param-y
#' @template param-delta
#' @template param-X
#' @template param-Z
#' @template param-beta0
#' @template param-gamma0
#' @template param-maxit
#' @template param-epsilon
#' @param splitE If \code{TRUE}, perform the E-step after each conditional M-step (see \strong{Details}).
#' @param nIRLS Number of IRLS steps to take before switching to Fisher scoring.  Can be 0 or greater than \code{maxit} to do only one or the other.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{beta}}{The MLE of the mean parameter vector.}
#'   \item{\code{gamma}}{The MLE of the variance parameter vector.}
#'   \item{\code{loglik}}{The value of the loglikelihood at the fitted parameter values.}
#'   \item{\code{iter}}{The number of steps taken by the algorithm.}
#'   \item{\code{error}}{The value of the loglikelihood relative error at the end of the algorithm.}
#' }
#'
#' @template details-hlm
#'
#' @details The fitting algorithm is an Expectation-Conditional-Maximization (ECM) algorithm extending the alternating weighted-LM/GLM updates of \code{beta} and \code{gamma}, proposed by Smyth (1989) for the uncensored setting.  The ECM algorithm terminates when either \code{maxit} iterations have been reached, or when
#' \preformatted{
#' |ll_curr - ll_prev| / (0.1 + |ll_curr|) < epsilon,
#' }
#' where \code{ll_curr} and \code{ll_prev} are the loglikelihood values at the current and previous iterations.
#'
#' \strong{TODO:}
#' \itemize{
#'   \item Input checking.
#'   \item \code{print}, \code{summary}, \code{vcov} methods.
#'   \item \code{residual} method.  Perhaps use expected lifetime for the censored observations?
#'   \item Separate into \code{chlm} and \code{hlm} classes?
#' }
#' @references Smyth, G.K. "Generalized Linear Models with Varying Dispersion." \emph{Journal of the Royal Statistical Society Series B} 51:1 (1989): 47-60.  \url{https://www.jstor.org/stable/2345840}.
#' @export
chlm_fit <- function(y, delta, X, Z, beta0, gamma0,
                     maxit = 100, epsilon = 1e-8,
                     splitE = FALSE, nIRLS = 5) {
  # helper functions
  # loglikelihood
  loglik <- function(beta, gamma) {
    chlm_loglik(beta, gamma, y, delta, X, Z)
  }
  # E-step
  Estep <- function(beta, sigmac) {
    muc <- c(Xc %*% beta)
    zc <- (yc - muc)/sigmac
    ## fsigc <- sigmac * etnorm(zc)
    fsigc <- sigmac * dnorm(zc)/pnorm(-zc)
    Rc <- fsigc + muc
    ## Sc1 <- sigmac^2 * e2tnorm(zc) + 2*muc*fsigc + muc^2
    Sc <- sigmac * (sigmac + zc*fsigc) + muc * (2*Rc - muc)
    ## list(Rc = Rc, Uc = Uc)
    R[!delta] <<- Rc
    S[!delta] <<- Sc
    invisible(NULL)
  }
  # initialize the algorithm
  if(missing(beta0)) {
    beta <- lm_fit(y = y, X = X)
  } else beta <- beta0
  if(missing(gamma0)) {
    gamma <- lvlm_fit(y2 = c(y - X %*% beta)^2, Z = Z, method = "LS")
  } else gamma <- gamma0
  # intermediate variables for censoring
  Xc <- X[!delta,,drop=FALSE]
  yc <- y[!delta]
  R <- y
  S <- y^2
  ll_old <- loglik(beta, gamma)
  if(all(delta)) {
    # no censoring
    out <- hlm_fit(y = y, X = X, Z = Z,
                   beta0 = beta, gamma0 = gamma,
                   maxit = maxit, epsilon = epsilon, method = "IRLS")
    ## class(out) <- "hlm"
    return(out)
  }
  for(ii in 1:maxit) {
    method <- if(ii <= nIRLS) "IRLS" else "Fisher"
    # E-step
    Zg <- c(Z %*% gamma)
    sigmac <- exp(.5 * Zg[!delta])
    Estep(beta, sigmac)
    # M-step: beta
    w <- exp(-Zg)
    beta <- lm_fit(y = R, X = X, w = w)
    if(splitE) Estep(beta, sigmac) # recompute E-step
    # M-step: gamma
    Xb <- c(X %*% beta)
    U <- S + Xb * (Xb - 2*R)
    U <- pmax(U, 1e-10) # FUDGE!!!
    gamma <- lvlm_fit(y2 = U, Z = Z, method = method)
    ll_new <- loglik(beta, gamma)
    ## tol <- rel_err(ll_new, ll_old)
    tol <- abs(ll_new - ll_old)/(.1 + abs(ll_new))
    if(tol < epsilon) break else ll_old <- ll_new
  }
  out <- list(beta = setNames(beta, NULL),
              gamma = setNames(gamma, NULL),
              loglik = ll_new, iter = ii, error = tol)
  ## class(out) <- "hlm"
  out
}

#' @rdname chlm_fit
#' @export
chlm_control <- function(epsilon = 1e-5, maxit = 100,
                         nIRLS = 5, splitE = TRUE) {
  list(epsilon = epsilon, maxit = maxit, nIRLS = nIRLS, splitE = splitE)
}


#--- helper functions ----------------------------------------------------------

## # relative error function used by stats::glm.fit
## rel_err <- function(x_new, x_old) {
##   abs(x_new - x_old)/(.1 + abs(x_new))
## }

## # truncated expectation of standard normal: E[Z|Z>a]
## etnorm <- function(a) {
##   return(dnorm(a)/pnorm(-a))
## }
## # E[Z^2|Z>a]
## e2tnorm <- function(a) {
##   return(1+a*dnorm(a)/pnorm(-a))
## }
## # fitting functions
## .lm_fit <- function(y, X) {
##   lm_fit(y = y, X = X)
## }
## .wlm_fit <- function(y, X, w) {
##   wlm_fit(y = y, X = X, w = w)
## }
## .lvlm_fit <- function(y2, Z, gamma0) {
##   if(missing(gamma0)) gamma0 <- lvlm_fitLS(logY2 = log(y2), Z = Z)
##   if(method == 0) {
##     # Fisher Scoring
##     lvlm_fitFS(y2 = y2, Z = Z, gamma0 = gamma0)
##   } else if(method == 1) {
##     # IRLS
##     lvlm_fitIRLS(y2 = y2, Z = Z, gamma0 = gamma0)
##   }
## }
