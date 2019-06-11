#' Fit the heteroscedastic linear model.
#'
#' @param formula An object of class \code{formula} (or one that can be coerced to that class).  See \strong{Details}.
#' @param data An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the variables are taken from \code{environment(formula)}.
#' @param subset An optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights Currently ignored.
#' @param na.action A function which indicates what should happen when the data contain \code{NA}s.  Same default behavior as in \code{stats::lm}.
#' @param control List of parameters to control the fitting process.  See \strong{Details}.
#' @param y,x,model Logical values indicating which data elements should be returned in the output.  These correspond to the response, the covariance matrices, and the \code{model.frame}, respectively.  See \code{stats::lm}.
#' @param offset Currently ignored, as are \code{offset} terms in formula.
#'
#' @details The heteroscedastic linear model (HLM) is of the form
#' \deqn{
#' y_i \mid \boldsymbol{x}_i, \boldsymbol{z}_i \stackrel{\mathrm{ind}}{\sim} \mathcal N\big(\boldsymbol{x}_i'\boldsymbol{\beta}, \exp(\boldsymbol{z}_i'\boldsymbol{\gamma})\big),
#' }{
#' y_i | x_i, z_i ~ind N(x_i'\beta, \exp(z_i'\gamma)),
#' }
#' where for each subject \eqn{i}, \eqn{y_i} is the response, and \eqn{\boldsymbol{x}_i \in \mathbb{R}^p}{x_i \in R^p} and \eqn{\boldsymbol{z}_i \in \mathbb{R}^q}{z_i \in R^q} are mean and variance covariate vectors, respectively.
#'
#' The \code{formula} term is specified as
#' \preformatted{
#' y ~ x1 + x2 | z1 + z2
#' }
#' where the vertical bar separates mean and variance components.  If no bar is found, an intercept variance term of the form \code{y ~ x1 + x2 | 1} is assumed, corresponding to the usual linear model (but with a different parametrization).
#'
#' Right censoring of observations is supported by specifying the response as a two-column matrix, where the first column is the response and the second column is a censoring status indicator with `0`: censored and `1`: uncensored.
#'
#' Fitting the \code{hlm} model is done by blockwise coordinate ascent, alternatively maximing the mean parameters by weighted least-squares, and the variance parameters either via Fisher scoring or Iteratively Reweighted Least Squares.  When there is right-censoring, these maximization steps are embedded within an Expectation-Conditional-Maximization algorithm.
#'
#' @return An object of class \code{hlm} with the following elements:
#' \describe{
#'   \item{\code{coefficients}}{A list with elements \code{beta} and \code{gamma} containing the coefficient MLEs.}
#'   \item{\code{loglik}}{The loglikelihood at the MLE.}
#'   \item{\code{model}}{If requested (the default), the model frame.}
#'   \item{\code{df.residual}}{The residual degrees of freedom.}
#'   \item{\code{iter}}{The number of iterations used in the fitting algorithm.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{terms}}{The \code{terms} object used.}
#'   \item{\code{y}}{If requested, the response used.}
#'   \item{\code{x}}{If requested, a list with elements \code{X} and \code{Z} giving the mean and variance model matrices uses.}
#'   \item{\code{model}}{If requested, the \code{model.frame} used.}
#'   \item{\code{na.action}}{(Where relevant) information returned by \code{model.frame} on the special handling of \code{NA}s.}
#' }
hlm <- function(formula, data, subset, weights, na.action,
                control, y, x, model, offset) {
  forms <- parse_formula(formula)
  # call parsing copied from lm with "| <- +" substitution
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf[[match("formula", names(mf), 0L)]] <- forms$full
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  ## mt <- attr(mf, "terms")
  # extract relevant terms for fitting
  y <- model.response(mf, "numeric")
  if(is.matrix(y)) {
    if(ncol(y) == 2) {
      delta <- y[,2]
      y <- y[,1]
    } else {
      stop("Response term in `formula` must be a numeric vector or two-column matrix.")
    }
  } else {
    delta <- rep(TRUE, length(y))
  }
  w <- as.vector(model.weights(mf))
  if(!is.null(w)) {
    warning("'weights' argument currently ignored.")
  }
  offset <- model.offset(mf)
  if(!is.null(offset)) {
    warning("'offset' term(s) currently ignored.")
  }
  # mean and variance model terms
  mtX <- terms(forms$X, data = data)
  mtZ <- terms(forms$Z, data = data)
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
}

#--- helper functions ----------------------------------------------------------

# parse hlm formula
# a list with three elements:
# - original formula, with | replaced by +.
#   this is used to get the original model frame.
# - mean and variance formulas, i.e., y ~ {left or right of | removed}.
#   the response is kept so that the "." formula argument gets parsed correctly.
parse_formula <- function(formula) {
  form_full <- as.character(formula)
  # extract formulas for each covariate
  form_y <- form_full[2]
  form_cov <- form_full[3]
  if(!grepl("[|]", form_cov)) {
    form_cov <- paste0(form_cov, " | 1")
  }
  form_cov <- strsplit(form_cov, "[|]")[[1]]
  form_cov <- paste0(form_y, " ~ ", form_cov)
  # convert | to + in full formula
  form_full  <- paste0(gsub("[|]", "+", form_full)[c(2,1,3)], collapse = " ")
  sapply(c(full = form_full, X = form_cov[1], Z = form_cov[2]),
         function(form) as.formula(form, env = environment(formula)))
}
