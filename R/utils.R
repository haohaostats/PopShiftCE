
#' @keywords internal
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

#' @keywords internal
safe_mse <- function(est, truth) {
  d <- (est - truth)
  d <- d[is.finite(d)]
  if (length(d) == 0L) return(NA_real_)
  mean(d^2)
}

#' HC covariance with safe fallbacks
#' @keywords internal
safe_vcov <- function(fit) {
  for (ty in c("HC3","HC1","HC0")) {
    vc <- try(sandwich::vcovHC(fit, type = ty), silent = TRUE)
    if (!inherits(vc, "try-error") && all(is.finite(vc))) return(vc)
  }
  sigma2 <- sum(stats::residuals(fit)^2) / stats::df.residual(fit)
  vc_ols <- sigma2 * summary(fit)$cov.unscaled
  return(vc_ols)
}
