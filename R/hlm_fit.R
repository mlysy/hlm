#' Low-level fitting function for the LM model.
#'
#' Calculates the MLE of the coefficients of the usual linear regression model (LM).
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
#' where for each subject \eqn{i}, \eqn{y_i} is the response, \eqn{\boldsymbol{x}_i \in \mathbb{R}^p}{x_i \in R^p} is the mean covariate vector, and \eqn{w_i} is an optional positive weight (defaults to 1).
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
#' Estimates the coefficients of a log-variance linear model (LVLM).  See \strong{Details}.
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
#' @details The log-variance linear model (LVLM) is defined as
#' \deqn{
#' y_i \mid \boldsymbol{z}_i \stackrel{\mathrm{ind}}{\sim} \mathcal N\big(0, \exp(\boldsymbol{z}_i'\boldsymbol{\gamma})\big),
#' }{
#' y_i | z_i ~ind N(0, exp(z_i'\gamma)),
#' }
#' where for each subject \eqn{i}, \eqn{y_i} is the response, and \eqn{\boldsymbol{z}_i \in \mathbb{R}^q}{z_i \in R^q} is the variance covariate vector.
#'
#' Three types of fitting algorithms for \eqn{\boldsymbol{\gamma}}{\gamma} are provided.  \code{method = Fisher} and \code{IRLS} are Fisher Scoring and Iteratively Reweighted Least-Squares MLE-finding algorithms, respectively.  The former is faster while the latter is more stable.  \code{method = LS} is a least-squares estimator, which is the fastest.  It is a consistent estimator but not as efficient as the MLE.
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

#' Low-level fitting function for the uncensored HLM model.
#'
#' Calculates the MLE of the coefficients of the heteroscedastic linear model (HLM).  See \strong{Details}.
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
#' @template details-hlm
#'
#' @details The low-level function \code{hlm_fit} assumes that the response vector is fully observed (uncensored).  See \code{\link{chlm_fit}} for the corresponding function with censoring, or the higher-level interface \code{\link{hlm}}.
#'
#' The tuning parameters of the LVLM fitting methods are tuned to their default values in \code{\link{lvlm_fit}}.
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
