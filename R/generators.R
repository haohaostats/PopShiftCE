
#' Stage-1 surrogate-only interim statistic T1
#' @param n1 per-arm Stage-1 size
#' @param muX_C control-arm mean of X
#' @param muX_T treatment-arm mean of X
#' @param sigmaX sd of X
#' @param gamma0,gamma1 external calibration coefficients for Yhat = gamma0 + gamma1 * X
#' @return list(T1, est1, se1)
#' @keywords internal
gen_stage1_and_T1 <- function(n1, muX_C, muX_T, sigmaX, gamma0, gamma1) {
  Xc <- stats::rnorm(n1, muX_C, sigmaX)
  Xt <- stats::rnorm(n1, muX_T, sigmaX)
  Yhat_c <- gamma0 + gamma1 * Xc
  Yhat_t <- gamma0 + gamma1 * Xt
  mC <- mean(Yhat_c); mT <- mean(Yhat_t)
  vC <- stats::var(Yhat_c); vT <- stats::var(Yhat_t)
  se1  <- sqrt(vT/n1 + vC/n1)
  est1 <- mT - mC
  T1   <- est1 / se1
  list(T1 = T1, est1 = est1, se1 = se1)
}

#' Stage-2 new participants under partial population shift
#' @param error_type one of "normal","t","skew"
#' @return data.frame with columns stage, T, Z, Y, is_new
#' @keywords internal
gen_stage2_new <- function(n2, mu0, delta, theta, eta, sigmaY,
                           piZ, error_type = c("normal","t","skew")) {
  error_type <- match.arg(error_type)
  gen_err <- function(n) {
    if (error_type == "normal") stats::rnorm(n, 0, sigmaY)
    else if (error_type == "t")  sigmaY * (stats::rt(n, df = 5) / sqrt(5/3))
    else                         sigmaY * (stats::rexp(n, 1) - 1)
  }
  Zc2 <- stats::rbinom(n2, 1, piZ)
  Zt2 <- stats::rbinom(n2, 1, piZ)
  Yc2 <- mu0 + theta*Zc2 + gen_err(n2)
  Yt2 <- mu0 + delta + theta*Zt2 + eta*Zt2 + gen_err(n2)
  data.frame(stage = "S2",
             T = c(rep(0, n2), rep(1, n2)),
             Z = c(Zc2,        Zt2),
             Y = c(Yc2,        Yt2),
             is_new = 1)
}
