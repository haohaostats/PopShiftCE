#' @keywords internal
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

#' @keywords internal
safe_mse <- function(est, truth) {
  d <- est - truth
  d <- d[is.finite(d)]
  if (length(d) == 0L) return(NA_real_)
  mean(d^2)
}

#' Heteroskedasticity-consistent covariance with safe fallbacks
#'
#' Attempts HC3, then HC1, then HC0. If all fail, falls back to OLS covariance.
#'
#' @param fit An `lm` object.
#' @return Covariance matrix.
#' @keywords internal
safe_vcov <- function(fit) {
  for (ty in c("HC3", "HC1", "HC0")) {
    vc <- try(sandwich::vcovHC(fit, type = ty), silent = TRUE)
    if (!inherits(vc, "try-error") && all(is.finite(vc))) return(vc)
  }
  sigma2 <- sum(stats::residuals(fit)^2) / stats::df.residual(fit)
  sigma2 * summary(fit)$cov.unscaled
}

#' @keywords internal
.kernel_weights <- function(x, x0, h, kernel = c("gaussian", "tricube")) {
  kernel <- match.arg(kernel)
  u <- (x - x0) / h
  if (kernel == "gaussian") {
    exp(-0.5 * u^2)
  } else {
    au <- abs(u)
    ifelse(au < 1, (1 - au^3)^3, 0)
  }
}

#' @keywords internal
.neff <- function(w) {
  sw <- sum(w)
  sw2 <- sum(w^2)
  if (!is.finite(sw) || !is.finite(sw2) || sw <= 0 || sw2 <= 0) return(0)
  (sw^2) / sw2
}
