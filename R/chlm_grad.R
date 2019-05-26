#' Gradient of the HLM loglikelihood.
#'
#' @template param-beta
#' @template param-gamma
#' @template param-y
#' @template param-delta
#' @template param-X
#' @template param-Z
#'
#' @return Vector of length \code{p+q} returning the gradient of the loglikelihood at the parameter values.
#' @export
chlm_grad <- function(beta, gamma, y, delta, X, Z) {
  Zg <- c(Z %*% gamma)
  Xb <- c(X %*% beta)
  sig <- exp(.5 * Zg)
  U <- (y - Xb)/sig
  # weights
  wgt <- matrix(NA, length(y), 2)
  if(any(delta)) {
    wgt[delta,1] <- U[delta]/sig[delta]
    wgt[delta,2] <- .5 * (U[delta] * U[delta] - 1)
  }
  if(any(!delta)) {
    A <- dnorm(U[!delta])/pnorm(U[!delta], lower.tail = FALSE)
    wgt[!delta,1] <- A/sig[!delta]
    wgt[!delta,2] <- .5 * A * U[!delta]
  }
  p <- length(beta)
  q <- length(gamma)
  grad <- rep(NA, p+q)
  grad[1:p] <- colSums(wgt[,1] * X)
  grad[p+1:q] <- colSums(wgt[,2] * Z)
  grad
}

#' Hessian of the HLM loglikelihood.
#'
#' @template param-beta
#' @template param-gamma
#' @template param-y
#' @template param-delta
#' @template param-X
#' @template param-Z
#'
#' @return Matrix of size \code{(p+q) x (p+q)} containg the Hessian of the loglikelihood at the parameter values.
#' @export
chlm_hess <- function(beta, gamma, y, delta, X, Z) {
  Zg <- c(Z %*% gamma)
  Xb <- c(X %*% beta)
  sig <- exp(.5 * Zg)
  U <- (y - Xb)/sig
  # weights
  wgt <- matrix(NA, length(y), 3)
  if(any(delta)) {
    wgt[delta,1] <- 1/sig[delta]
    wgt[delta,2:3] <- U[delta]
    wgt[delta,3] <- .5 * wgt[delta,3] * U[delta]
    wgt[delta,1:2] <-  wgt[delta,1:2] * wgt[delta,1]
  }
  if(any(!delta)) {
    A <- dnorm(U[!delta])/pnorm(U[!delta], lower.tail = FALSE)
    B <- A * (A - U[!delta])
    wgt[!delta,1] <- 1/sig[!delta]
    wgt[!delta,2] <- .5 * (B * U[!delta] + A)
    wgt[!delta,3] <- .5 * U[!delta] * wgt[!delta,2]
    wgt[!delta,1:2] <-  wgt[!delta,1:2] * wgt[!delta,1]
    wgt[!delta,1] <-  B * wgt[!delta,1]
  }
  # hessian matrix
  p <- length(beta)
  q <- length(gamma)
  Hess <- matrix(NA, p+q, p+q)
  Hess[1:p,1:p] <- crossprod(X, wgt[,1] * X)
  Hess[1:p,p+1:q] <- crossprod(X, wgt[,2] * Z)
  Hess[p+1:q,1:p] <- t(Hess[1:p,p+1:q])
  Hess[p+1:q,p+1:q] <- crossprod(Z, wgt[,3] * Z)
  -Hess
}
