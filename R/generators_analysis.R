#' Stage-1 surrogate-only interim statistic
#'
#' Generates surrogate values `X` for control and treatment arms, constructs the
#' external projection `Yhat = gamma0 + gamma1 * X`, and returns the interim
#' two-sample z-statistic.
#'
#' @param n1 Per-arm Stage-1 sample size.
#' @param muX_C Control-arm mean of the surrogate.
#' @param muX_T Treatment-arm mean of the surrogate.
#' @param sigmaX Standard deviation of the surrogate.
#' @param gamma0,gamma1 External calibration coefficients.
#'
#' @return A list with `T1`, `est1`, and `se1`.
#' @export
gen_stage1_and_T1 <- function(n1, muX_C, muX_T, sigmaX, gamma0, gamma1) {
  Xc <- stats::rnorm(n1, muX_C, sigmaX)
  Xt <- stats::rnorm(n1, muX_T, sigmaX)
  Yhat_c <- gamma0 + gamma1 * Xc
  Yhat_t <- gamma0 + gamma1 * Xt
  mC <- mean(Yhat_c); mT <- mean(Yhat_t)
  vC <- stats::var(Yhat_c); vT <- stats::var(Yhat_t)
  se1  <- sqrt(vT / n1 + vC / n1)
  est1 <- mT - mC
  T1   <- est1 / se1
  list(T1 = T1, est1 = est1, se1 = se1)
}

#' Generate Stage-2 new participants under partial population shift
#'
#' @param n2 Per-arm Stage-2 sample size.
#' @param mu0 Baseline intercept for the primary outcome.
#' @param delta Treatment main effect in the `Z=0` stratum.
#' @param theta Main effect of `Z`.
#' @param eta Treatment-by-`Z` interaction.
#' @param sigmaY Standard deviation scale for outcome errors.
#' @param piZ Stage-2 prevalence of the post-amendment subgroup (`Z=1`).
#' @param error_type One of `"normal"`, `"t"`, or `"skew"`.
#'
#' @return A data.frame with columns `stage`, `T`, `Z`, `Y`, `is_new`.
#' @export
gen_stage2_new <- function(n2, mu0, delta, theta, eta, sigmaY,
                           piZ, error_type = c("normal", "t", "skew")) {
  error_type <- match.arg(error_type)

  gen_err <- function(n) {
    if (error_type == "normal") {
      stats::rnorm(n, 0, sigmaY)
    } else if (error_type == "t") {
      sigmaY * (stats::rt(n, df = 5) / sqrt(5 / 3))
    } else {
      sigmaY * (stats::rexp(n, 1) - 1)
    }
  }

  Zc2 <- stats::rbinom(n2, 1, piZ)
  Zt2 <- stats::rbinom(n2, 1, piZ)
  Yc2 <- mu0 + theta * Zc2 + gen_err(n2)
  Yt2 <- mu0 + delta + theta * Zt2 + eta * Zt2 + gen_err(n2)

  data.frame(
    stage = "S2",
    T = c(rep(0, n2), rep(1, n2)),
    Z = c(Zc2, Zt2),
    Y = c(Yc2, Yt2),
    is_new = 1
  )
}

#' Final analysis and one-sided Wald z-statistic for the marginal estimand
#'
#' Fits `Y ~ T + Z + T:Z` using OLS and heteroskedasticity-consistent covariance
#' (HC3 fallback to HC1/HC0/OLS). Supports two target definitions:
#'
#' * `pi_target = "fixed"`: design-fixed target prevalence `pi_fixed`
#' * `pi_target = "asObservedS2"`: uses the observed Stage-2 prevalence `pi_hat`
#'
#' For `pi_target = "asObservedS2"`, the variance includes a simple additional
#' term `eta_hat^2 * Var(pi_hat)` with `Var(pi_hat) = pi_hat(1-pi_hat)/n_S2`.
#' This function returns the final-stage test statistic and estimation output;
#' interpretation of null rejection frequencies (e.g., one-sided Type I error rate
#' or power) is handled at the simulation-summary level.
#'
#' @param df_full A data.frame with columns `Y`, `T`, `Z`, and `stage`.
#' @param pi_target Target type: `"fixed"` or `"asObservedS2"`.
#' @param pi_fixed Design-fixed prevalence used when `pi_target = "fixed"`.
#'
#' @return A list with elements `Zf`, `Delta_hat`, `se`, and `pi_hat`.
#' @export
final_analysis <- function(df_full,
                           pi_target = c("fixed", "asObservedS2"),
                           pi_fixed = NULL) {
  pi_target <- match.arg(pi_target)

  fit <- stats::lm(Y ~ T + Z + T:Z, data = df_full)
  vc  <- safe_vcov(fit)
  cf  <- stats::coef(fit)

  delta_hat <- if ("T" %in% names(cf))   unname(cf["T"])   else 0
  eta_hat   <- if ("T:Z" %in% names(cf)) unname(cf["T:Z"]) else 0

  if (pi_target == "fixed") {
    if (is.null(pi_fixed)) stop("pi_fixed must be provided when pi_target='fixed'.")
    pi_hat <- pi_fixed
    var_pi <- 0
  } else {
    zS2 <- df_full$Z[df_full$stage == "S2"]
    if (length(zS2) == 0L) {
      pi_hat <- 0
      var_pi <- 0
    } else {
      pi_hat <- mean(zS2)
      var_pi <- pi_hat * (1 - pi_hat) / length(zS2)
    }
  }

  Delta_hat <- delta_hat + eta_hat * pi_hat

  pick <- intersect(c("T", "T:Z"), colnames(vc))
  if (length(pick) == 0L) {
    var_D <- NA_real_
    se_D  <- NA_real_
    Zf    <- NA_real_
  } else {
    Sigma <- vc[pick, pick, drop = TRUE]
    if (is.null(dim(Sigma))) {
      Sigma <- matrix(Sigma, nrow = 1, ncol = 1, dimnames = list(pick, pick))
    }
    g_all <- c("T" = 1, "T:Z" = pi_hat)
    g     <- matrix(unname(g_all[pick]), nrow = 1)
    var_D_core <- as.numeric(g %*% Sigma %*% t(g))
    var_D      <- var_D_core + (eta_hat^2) * var_pi
    if (!is.finite(var_D) || var_D < 0) var_D <- NA_real_
    se_D <- if (is.na(var_D)) NA_real_ else sqrt(var_D)
    Zf   <- if (is.na(se_D) || se_D == 0) NA_real_ else as.numeric(Delta_hat / se_D)
  }

  list(Zf = Zf, Delta_hat = Delta_hat, se = se_D, pi_hat = pi_hat)
}

#' One-sided CE-consistent lower confidence limit (final stage)
#'
#' @param Delta_hat Final-stage estimated marginal effect.
#' @param se Final-stage standard error.
#' @param c_t1 Conditional critical value evaluated at the realized interim statistic.
#'
#' @return Numeric scalar.
#' @export
lcl_final <- function(Delta_hat, se, c_t1) {
  if (is.na(se) || is.na(Delta_hat) || is.na(c_t1)) return(NA_real_)
  Delta_hat - c_t1 * se
}

#' One-sided CE-consistent lower confidence limit (early stop; surrogate scale)
#'
#' @param est1 Interim surrogate-projected estimate.
#' @param se1 Interim standard error.
#' @param b_ref Reference boundary.
#'
#' @return Numeric scalar.
#' @export
lcl_interim <- function(est1, se1, b_ref) {
  if (is.na(se1) || is.na(est1) || is.na(b_ref)) return(NA_real_)
  est1 - b_ref * se1
}
