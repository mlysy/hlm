#' Loglikelihood of the HLM model.
#'
#' @template param-beta
#' @template param-gamma
#' @template param-y
#' @template param-delta
#' @template param-X
#' @template param-Z
#'
#' @return The loglikelihood evaluated at the parameters values (a scalar).
#' @export
chlm_loglik <- function(beta, gamma, y, delta, X, Z) {
  n <- length(y)
  if(missing(delta)) delta <- rep(TRUE, n)
  mu <- c(X %*% beta)
  sig <- exp(.5 * c(Z %*% gamma))
  ll <- rep(NA, length(y))
  if(any(delta)) {
    ll[delta] <- dnorm(y[delta], mu[delta], sig[delta], log = TRUE)
  }
  if(any(!delta)) {
    ll[!delta] <- pnorm(y[!delta], mu[!delta], sig[!delta],
                         lower.tail = FALSE, log.p = TRUE)
  }
  sum(ll)
}

