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

