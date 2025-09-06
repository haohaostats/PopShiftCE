
#' Final analysis and one-sided Wald z-statistic for the design-fixed estimand
#'
#' Fits OLS with heteroskedasticity-consistent covariance (HC3 fallback to HC1/HC0)
#' to the model Y ~ T + Z + T:Z, and computes
#' \deqn{\widehat{\Delta}^\dagger_{\mathrm{marg}}=\hat\delta + \hat\eta\,\pi_Z^\dagger}
#' with variance \eqn{[1,\ \pi_Z^\dagger]\ \widehat\Sigma_{(\hat\delta,\hat\eta)}\ [1,\ \pi_Z^\dagger]^T}.
#'
#' @param df_full data.frame with columns: Y, T, Z, stage
#' @param pi_fixed design-fixed mixture prevalence in [0,1]
#' @return list with elements: Zf (Wald z), Delta_hat, se
#' @export
final_analysis <- function(df_full, pi_fixed = 0.5) {
  fit <- stats::lm(Y ~ T + Z + T:Z, data = df_full)
  vc  <- safe_vcov(fit)
  cf  <- stats::coef(fit)
  
  delta_hat <- if ("T"   %in% names(cf)) unname(cf["T"])   else 0
  eta_hat   <- if ("T:Z" %in% names(cf)) unname(cf["T:Z"]) else 0
  
  Delta_hat <- delta_hat + eta_hat * pi_fixed
  
  pick <- intersect(c("T","T:Z"), colnames(vc))
  if (length(pick) == 0) {
    var_D <- NA_real_; se_D <- NA_real_; Zf <- NA_real_
  } else {
    Sigma <- vc[pick, pick, drop = TRUE]
    if (is.null(dim(Sigma))) Sigma <- matrix(Sigma, nrow = 1, ncol = 1,
                                             dimnames = list(pick, pick))
    g_all <- c("T" = 1, "T:Z" = pi_fixed)
    g     <- matrix(unname(g_all[pick]), nrow = 1)
    var_D <- as.numeric(g %*% Sigma %*% t(g))
    if (!is.finite(var_D) || var_D < 0) var_D <- NA_real_
    se_D <- if (is.na(var_D)) NA_real_ else sqrt(var_D)
    Zf   <- if (is.na(se_D) || se_D == 0) NA_real_ else as.numeric(Delta_hat / se_D)
  }
  list(Zf = Zf, Delta_hat = Delta_hat, se = se_D)
}

#' One-sided CE-consistent LCL (final stage)
#' @param Delta_hat numeric
#' @param se numeric
#' @param c_t1 conditional critical value c(t1)
#' @return numeric
#' @export
lcl_final <- function(Delta_hat, se, c_t1) {
  if (is.na(se) || is.na(Delta_hat) || is.na(c_t1)) return(NA_real_)
  Delta_hat - c_t1 * se
}

#' One-sided CE-consistent LCL (early stop, surrogate scale)
#' @param est1 interim estimate (surrogate-projected)
#' @param se1  interim se
#' @param b_ref reference boundary
#' @return numeric
#' @export
lcl_interim <- function(est1, se1, b_ref) {
  if (is.na(se1) || is.na(est1) || is.na(b_ref)) return(NA_real_)
  est1 - b_ref * se1
}
