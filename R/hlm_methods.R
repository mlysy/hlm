# hlm methods

#' @export
nobs.hlm <- function(object, ...) {
  object$df.residual + sum(sapply(object$coefficients, length))
}

#' @export
vcov.hlm <- function(object, ...) {
  V <- qr.solve(object$qr)
  dn <- c(paste0("X_", names(object$coefficients$beta)),
          paste0("Z_", names(object$coefficients$gamma)))
  colnames(V) <- rownames(V) <- dn
  V
}

#' @export
print.hlm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  ## Call formula
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  ## Show beta
  if (length(coefficients(x)$beta)) {
    cat("Mean Coefficients:\n")
    print.default(format(coefficients(x)$beta, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No mean coefficients\n")
  cat("\n")
  ## invisible(x)
  ## Show gamma
  if (length(coefficients(x)$gamma)) {
    cat("Variance Coefficients:\n")
    print.default(format(coefficients(x)$gamma, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No variance coefficients\n")
  cat("\n")
  invisible(x)
}

#' Summary method for \code{hlm} objects.
#'
#' @name summary.hlm
#' @param object An object of class \code{hlm}.
#' @param correlation If \code{TRUE}, calculate the correlation matrix of the estimated parameters.
#' @param symbolic.cor If \code{TRUE}, print correlations in symbolic form (see \code{stats::symnum}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{call}}{The call issuing the HLM \code{object}.}
#'   \item{\code{terms}}{The terms component from \code{object}.}
#'   \item{\code{df.residual}}{The component from \code{object}.}
#'   \item{\code{coefficients}}{List of two matrices of coefficients, standard errors, z-values, and p-values.}
#'   \item{\code{cov.unscaled}}{Estimated covariance matrix of the estimated coefficients.}
#'   \item{\code{correlation}}{If \code{correlation = TRUE}, the estimate correlations between the estimated coefficients.}
#'   \item{\code{symbolic.cor}}{If \code{correlation = TRUE}, the value of \code{symbolic.cor}.}
#'   \item{\code{aic}}{The AIC statistic.}
#'   \item{\code{iter}}{The component from \code{object}.}
#'   \item{\code{na.action}}{The component from \code{object}.}
#' }
#' @export
summary.hlm <- function(object, correlation = FALSE,
                        symbolic.cor = FALSE, ...) {
  z <- object
  if (is.null(z$terms))
    stop("Invalid 'hlm' object:  no 'terms' component.")
  ans <- z[c("call", "terms", "df.residual", "iter", "na.action")]
  # coefficients
  R <- vcov(object)
  dn <- rownames(R)
  theta <- coefficients(object)
  nbeta <- length(theta$beta)
  ngamma <- length(theta$gamma)
  est <- unlist(theta)
  se <- sqrt(diag(R))
  zval <- est/se
  coef.table <- cbind(Estimate = est, `Std. Error` = se,
                      `z value` = zval, `Pr(>|z|)` = 2 * pnorm(-abs(zval)))
  coefX <- coef.table[1:nbeta,,drop=FALSE]
  rownames(coefX) <- names(theta$beta)
  coefZ <- coef.table[nbeta+1:ngamma,,drop=FALSE]
  rownames(coefZ) <- names(theta$gamma)
  ## rownames(coef.table) <- dn
  ## ans$coefficients <- coef.table
  ans$coefficients <- list(beta = coefX, gamma = coefZ)
  ans$cov.unscaled <- R
  # aic
  ans$aic <- -2 * z$loglik + 2 * z$df.residual
  if(correlation) {
    dd <- sqrt(diag(R))
    ans$correlation <- R/outer(dd, dd)
    ans$symbolic.cor <- symbolic.cor
  }
  class(ans) <- "summary.hlm"
  ans
}

#' @rdname summary.hlm
#'
#' @param x An object of class \code{summary.hlm}.
#' @param digits The number of significant digits to use when printing.
#' @param signif.stars If \code{TRUE}, significance stars are printed for each coefficient.
#' @export
print.summary.hlm <- function(x, digits = max(3L, getOption("digits") - 3L),
                              symbolic.cor = x$symbolic.cor,
                              signif.stars = getOption("show.signif.stars"),
                              ...) {
  ## Call formula
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  ## coefficients
  cat("\nMean Coefficients:\n")
  coefs <- x$coefficients$beta
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  cat("\nVariance Coefficients:\n")
  coefs <- x$coefficients$gamma
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  if (nzchar(mess <- naprint(x$na.action)))
    cat("\n  (", mess, ")\n", sep = "")
  # aic
  cat("\nAIC: ", format(x$aic, digits = max(4L, digits + 1L)),
      "\n", "Number of algorithm iterations: ", x$iter,
      "\n", sep = "")
  # correlation
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2L), nsmall = 2L,
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}
