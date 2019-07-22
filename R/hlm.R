#' Fit the heteroscedastic linear model.
#'
#' Calculates the MLE of a heteroscedastic linear model (HLM) with possibly right-censored responses.  See \strong{Details}.
#'
#' @param formula An object of class \code{formula} (or one that can be coerced to that class).  See \strong{Details}.
#' @param data An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the variables are taken from \code{environment(formula)}.
#' @param subset An optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights Currently ignored.
#' @param na.action A function which indicates what should happen when the data contain \code{NA}s.  Same default behavior as in \code{stats::lm}.
#' @param control List of parameters to control the fitting process.  See \strong{Details}.
#' @param model,x,y,qr Logical values indicating which data elements should be returned in the output.  These correspond to the response, the covariance matrices, the \code{model.frame}, and the QR decomposition of the hessian of the negative loglikelihood, respectively.  See \code{stats::lm}.
#' @param offset Currently ignored, as are \code{offset} terms in the formula.
#'
#' @template details-hlm
#'
#' @details The \code{formula} term is specified as e.g.,
#' \preformatted{
#' y ~ x1 + x2 | z1 + z2
#' }
#' where the vertical bar separates mean and variance components.  If no bar is provided, an intercept variance term of the form \code{y ~ x1 + x2 | 1} is assumed, corresponding to the usual linear model (but with a different parametrization).
#'
#' Right censoring of observations is supported by specifying the response as a two-column matrix, where the first column is the response and the second column is a censoring status indicator with \code{0}: right-censored and \code{1}: uncensored.
#'
#' Fitting the \code{hlm} model is done by blockwise coordinate ascent, alternatively maximing the mean parameters by weighted least-squares, and the variance parameters either via Fisher scoring or Iteratively Reweighted Least Squares.  When there is right-censoring, these maximization steps are embedded within an Expectation-Conditional-Maximization algorithm.
#'
#' @return An object of class \code{hlm} with the following elements:
#' \describe{
#'   \item{\code{coefficients}}{A list with elements \code{beta} and \code{gamma} containing the coefficient MLEs.}
#'   \item{\code{loglik}}{The loglikelihood at the MLE.}
#'   \item{\code{df.residual}}{The residual degrees of freedom.}
#'   \item{\code{iter}}{The number of iterations used in the fitting algorithm.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{terms}}{A list with elements \code{X} and \code{Z} giving the \code{terms} object for mean and variance models.}
#'   \item{\code{y}}{If requested, the response used.}
#'   \item{\code{x}}{If requested, a list with elements \code{X} and \code{Z} giving the mean and variance model matrices uses.}
#'   \item{\code{model}}{If requested, the \code{model.frame} used.}
#'   \item{\code{na.action}}{(Where relevant) information returned by \code{model.frame} on the special handling of \code{NA}s.}
#'   \item{\code{qr}}{If requested, the QR decomposition of the observed Fisher information matrix of size \code{(p+q) x (p+q)}.}
#' }
#'
#' @note \strong{Warning:} At present \code{hlm} cannot handle pure LM or LVLM models.  For datasets without censoring, please use the low-level functions \code{\link{lm_fit}} and \code{\link{lvlm_fit}} instead.
#'
#' @seealso Current methods for \code{hlm} objects are: \code{print}, \code{nobs}, \code{vcov}, and \code{\link{summary}}.
#'
#' @example examples/hlm.R
#'
#' @references Wang, Y., You, T., and Lysy, M. "A heteroscedastic accelerated failure time model for survival analysis" (2019): \url{https://arxiv.org/abs/1508.05137}.
#' @export
hlm <- function(formula, data, subset, weights, na.action,
                control, model = TRUE, qr = TRUE,
                x = FALSE, y = FALSE, offset) {
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
  ret_y <- y
  ret_x <- x
  resp <- model.response(mf, "numeric")
  if(is.matrix(resp)) {
    if(ncol(resp) == 2) {
      delta <- resp[,2]
      if(!all(delta %in% c(TRUE, FALSE))) {
        stop("Second column of response must be logical vector.")
      }
      delta <- as.logical(delta)
      y <- resp[,1]
    } else {
      stop("Response term in `formula` must be a numeric vector or two-column matrix.")
    }
  } else {
    y <- resp
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
  if(missing(control)) control <- chlm_control()
  # mean and variance model terms
  mtX <- if(!missing(data)) terms(forms$X, data = data) else terms(forms$X)
  mtZ <- if(!missing(data)) terms(forms$Z, data = data) else terms(forms$Z)
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  Xnames <- colnames(X)
  Znames <- colnames(Z)
  cfit <- chlm_fit(y = y, delta = delta, X = X, Z = Z,
                   maxit = control$maxit,
                   epsilon = control$epsilon,
                   splitE = control$splitE,
                   nIRLS = control$nIRLS)
  if(cfit$iter >= control$maxit || cfit$error > control$epsilon) {
    warning("Fitting algorithm did not converge.")
  }
  # format output
  beta <- setNames(cfit$beta, Xnames)
  gamma <- setNames(cfit$gamma, Znames)
  out <- list(coefficients = list(beta = beta, gamma = gamma),
              loglik = cfit$loglik, iter = cfit$iter,
              df.residual = length(y) - length(Xnames) - length(Znames),
              call = cl, terms = list(X = mtX, Z = mtZ))
  if(ret_y) out$y <- resp
  if(ret_x) out$x <- list(X = X, Z = Z)
  if(model) out$model <- mf
  out$na.action <- attr(mf, "na.action")
  if(qr) out$qr <- qr(-chlm_hess(beta = beta, gamma = gamma,
                                 y = y, delta = delta, X = X, Z = Z))
  class(out) <- "hlm"
  out
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
