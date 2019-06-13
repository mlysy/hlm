#' Low-level fitting function for the LM model.
#'
#' @template param-y
#' @template param-X
#' @param w Optional positive weights vector of length \code{n}.
#'
#' @return MLE estimate of \code{beta} as a vector of length \code{p}.
#'
#' @details The LM model is defined as
#' \deqn{
#' \boldsymbol{y}_i \mid \boldsymbol{x}_i \stackrel{\mathrm{ind}}{\sim} \mathcal{N}(\boldsymbol{x}_i'\boldsymbol{\beta}, 1/w_i),
#' }{
#' y_i | x_i ~ind N(x_i'\beta, 1/w_i),
#' }
#' where ...
#'
#' @template warn-cpp
#' @export
lm_fit <- function(y, X, w) {
  if(missing(w)) {
    beta <- LM_Fit(y = y, X = X)
  } else {
    beta <- WLM_Fit(y = y, X = X, w = w)
  }
  beta
}

#' Low-level fitting function for the LVLM model.
#'
#' @param y2 Square of response vector of length \code{n}.
#' @template param-Z
#' @param method Which fitting algorithm to use.  See \strong{Details}.
#' @param gamma0 Initial variance parameter vector of length \code{q}.  If missing a least-squares estimate is used (see \strong{Details}).
#' @template param-maxit
#' @template param-epsilon
#'
#' @return The MLE (or least-squares estimate) of \code{gamma} as a vector of length \code{q}.
#'
#' @details the log-variance linear model (LVLM) is defined as...
#'
#' Three types of fitting algorithms are provided.  \code{method = Fisher} and \code{IRLS} are Fisher Scoring and Iteratively Reweighted Least-Squares MLE-finding algorithms, respectively.  The former is faster while the latter is more stable.  \code{method = LS} is a least-squares estimator, which is the fastest.  It is a consistent estimator but not as efficient as the MLE.
#'
#' @template warn-cpp
#' @export
lvlm_fit <- function(y2, Z, method = c("IRLS", "Fisher", "LS"), gamma0,
                     maxit = 25, epsilon = 1e-8) {
  method <- match.arg(method)
  if(method == "Fisher") {
    if(missing(gamma0)) {
      gamma <- LVLM_FitFS(y2 = y2, Z = Z,
                          maxit = maxit, epsilon = epsilon,
                          gamma0 = 0, initLS = TRUE)$gamma
    } else {
      gamma <- LVLM_FitFS(y2 = y2, Z = Z, maxit = maxit, epsilon = epsilon,
                          gamma0 = gamma0, initLS = FALSE)$gamma
    }
  } else if(method == "IRLS") {
    if(missing(gamma0)) {
      gamma <- LVLM_FitIRLS(y2 = y2, Z = Z,
                            maxit = maxit, epsilon = epsilon,
                            gamma0 = 0, initLS = TRUE)$gamma
    } else {
      gamma <- LVLM_FitIRLS(y2 = y2, Z = Z,
                            maxit = maxit, epsilon = epsilon,
                            gamma0 = gamma0, initLS = FALSE)$gamma
    }
  } else if(method == "LS") {
    gamma <- LVLM_FitLS(logY2 = log(y2), Z = Z)
  }
  gamma
}

#' Low-level fitting function for the HLM model.
#'
#' @template param-y
#' @template param-X
#' @template param-Z
#' @param method Which method to use for fitting the conditional LVLM model.  See \code{\link{lvlm_fit}}.
#' @template param-beta0
#' @template param-gamma0
#' @template param-maxit
#' @template param-epsilon
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{beta}}{The MLE of the mean parameters as a vector of length \code{p}.}
#'   \item{\code{gamma}}{The MLE of the variance parameters as a vector of length \code{q}.}
#'   \item{\code{loglik}}{The loglikelihood at the final step of the algorithm.}
#'   \item{\code{iter}}{The number of iterations of the fitting algorithm.}
#'   \item{\code{tolerance}}{The loglikelihood relative error at the last step.}
#' }
#'
#' @details The tuning parameters of the LVLM fitting methods are tuned to their default values in \code{\link{lvlm_fit}}.
#'
#' @template warn-cpp
#' @export
hlm_fit <- function(y, X, Z, method = c("IRLS", "Fisher"),
                    beta0, gamma0, maxit = 100, epsilon = 1e-5) {
  if(missing(beta0)) beta0 <- LM_Fit(y = y, X = X)
  if(missing(gamma0)) {
    gamma0 <- LVLM_FitLS(logY2 = log(c(y - X %*% beta0)^2), Z = Z)
  }
  method <- match.arg(method)
  method <- switch(method, Fisher = 0, IRLS = 1)
  HLM_Fit(y = y, X = X, Z = Z, beta0 = beta0, gamma0 = gamma0,
          maxit = maxit, epsilon = epsilon, method = method)
}
