
#' Simulate a single two-stage trial using the MC-CE decision rule
#'
#' This simulates one seamless two-stage trial:
#' Stage 1 uses a surrogate-only interim (via \code{gen_stage1_and_T1}).
#' If early efficacy is declared (T1 >= b_ref), it stops and returns the
#' surrogate-scale one-sided LCL. Otherwise it proceeds to Stage 2, fits the
#' full-data primary model (via \code{final_analysis}), and applies the
#' conditional critical value \code{c(T1)} from the CE lookup.
#'
#' @param n1,n2 Stage-1/Stage-2 per-arm sample sizes.
#' @param delta Treatment main effect in the Z=0 stratum.
#' @param piZ Stage-2 shift prevalence for Z=1.
#' @param theta Main effect of Z.
#' @param eta Interaction T:Z.
#' @param sigmaY SD of the primary outcome error.
#' @param error_type One of \code{"normal"}, \code{"t"}, \code{"skew"}.
#' @param muX_C Stage-1 surrogate mean in control.
#' @param muX_T Stage-1 surrogate mean in treatment. If \code{NULL}, it is set
#'   internally to \code{muX_C + delta/gamma1} to align the interim mean with the
#'   primary effect as in the paper.
#' @param sigmaX SD of the Stage-1 surrogate.
#' @param gamma0,gamma1 External linear calibration coefficients for the surrogate.
#' @param mu0 Baseline intercept for the primary outcome.
#' @param lookup A CE lookup object returned by \code{build_ce_lookup()}.
#' @param pi_fixed Design-fixed mixture prevalence used in the estimand.
#'
#' @return A list with fields:
#' \describe{
#'   \item{reject}{Logical, whether H0 was rejected (either early or final).}
#'   \item{early_stop}{Logical, whether early efficacy was triggered.}
#'   \item{final_est}{Final-stage estimate of the marginal effect (NA if early stop).}
#'   \item{final_se}{Final-stage SE (NA if early stop).}
#'   \item{lcl_final}{Final CE-consistent one-sided lower bound (NA if early stop).}
#'   \item{lcl_interim}{Interim one-sided lower bound on the surrogate scale (NA if not early stop).}
#'   \item{coverage_final}{Indicator for final LCL covering the true effect (NA if early stop).}
#'   \item{coverage_overall}{Overall coverage (interim or final bound, depending on stop).}
#'   \item{sample_used_total}{Total sample size used across arms.}
#'   \item{sample_used_per_arm}{Per-arm sample size used.}
#'   \item{z1}{Observed interim statistic T1 (for plotting/diagnostics).}
#'   \item{zf}{Final-stage z-statistic S_final (NA if early stop).}
#'   \item{cz1}{Conditional critical value c(T1) used at final, or b_ref if early stop.}
#' }
#' @export
simulate_trial_ce <- function(n1, n2,
                              delta, piZ, theta, eta,
                              sigmaY, error_type = c("normal","t","skew"),
                              muX_C, muX_T = NULL, sigmaX,
                              gamma0, gamma1,
                              mu0,
                              lookup,
                              pi_fixed = 0.5) {
  
  error_type <- match.arg(error_type)
  
  # Auto-align interim signal if muX_T not provided:
  # muX_T = muX_C + delta/gamma1 (consistent with paper & examples)
  if (is.null(muX_T)) muX_T <- muX_C + delta / gamma1
  
  # ---- Stage 1: surrogate-only interim ----
  s1 <- gen_stage1_and_T1(n1, muX_C, muX_T, sigmaX, gamma0, gamma1)
  z1  <- s1$T1
  est1 <- s1$est1
  se1  <- s1$se1
  
  # Early stop branch (surrogate-scale bound)
  if (!is.na(z1) && z1 >= lookup$b_ref) {
    lcl1 <- lcl_interim(est1, se1, lookup$b_ref)
    # True marginal effect for design-fixed pi_fixed
    true_delta_marg <- delta + eta * pi_fixed
    covg_sur <- is.finite(lcl1) && (true_delta_marg >= lcl1)
    
    return(list(
      reject = TRUE,
      early_stop = TRUE,
      final_est = NA_real_,
      final_se  = NA_real_,
      lcl_final = NA_real_,
      lcl_interim = lcl1,
      coverage_final = NA,
      coverage_overall = covg_sur,
      sample_used_total   = 2 * n1,
      sample_used_per_arm = n1,
      z1  = z1,
      zf  = NA_real_,
      cz1 = lookup$b_ref
    ))
  }
  
  # ---- Continue to Stage 2 ----
  gen_err <- function(n) {
    if (error_type == "normal") {
      stats::rnorm(n, 0, sigmaY)
    } else if (error_type == "t") {
      sigmaY * (stats::rt(n, df = 5) / sqrt(5/3))
    } else { # "skew"
      sigmaY * (stats::rexp(n, 1) - 1)
    }
  }
  
  # Stage-1 matured Y (Z=0), both arms
  df1 <- data.frame(
    stage = "S1",
    T = c(rep(0, n1), rep(1, n1)),
    Z = 0,
    Y = c(mu0 + gen_err(n1),
          mu0 + delta + gen_err(n1)),
    is_new = 0
  )
  
  # Stage-2 new participants with partial shift
  df2 <- gen_stage2_new(n2, mu0, delta, theta, eta, sigmaY, piZ, error_type)
  
  df <- rbind(df1, df2)
  
  # Final analysis (design-fixed mixture target)
  fa <- final_analysis(df, pi_fixed = pi_fixed)
  Zf  <- fa$Zf
  est <- fa$Delta_hat
  se  <- fa$se
  
  # Conditional critical value c(T1) from lookup
  cz1 <- lookup$c_fun(z1)
  
  # Final decision & final LCL
  reject <- is.finite(Zf) && is.finite(cz1) && (Zf >= cz1)
  lcl_f  <- lcl_final(est, se, cz1)
  
  # True marginal effect for design-fixed pi_fixed
  true_delta_marg <- delta + eta * pi_fixed
  covg_final <- is.finite(lcl_f) && (true_delta_marg >= lcl_f)
  
  list(
    reject = reject,
    early_stop = FALSE,
    final_est = est,
    final_se  = se,
    lcl_final = lcl_f,
    lcl_interim = NA_real_,
    coverage_final = covg_final,
    coverage_overall = covg_final,  # overall = final here
    sample_used_total   = 2 * (n1 + n2),
    sample_used_per_arm = (n1 + n2),
    z1  = z1,
    zf  = Zf,
    cz1 = cz1
  )
}

#' Simulate many trials and summarize operating characteristics
#'
#' Returns two summaries:
#' \itemize{
#'   \item \code{summary}: numeric values (good for downstream calculations);
#'   \item \code{summary_pretty}: fixed-decimal character values (good for display).
#' }
#'
#' Default decimal places:
#' \itemize{
#'   \item Rates (rejection_rate, early_stop_rate, coverage): 3 decimals;
#'   \item Effects (bias, mse): 4 decimals;
#'   \item Sample sizes (ASN_total, ASN_per_arm): 1 decimal.
#' }
#'
#' @param R Number of replicates.
#' @inheritParams simulate_trial_ce
#' @param digits_rate Integer decimals for rate-type outputs.
#' @param digits_effect Integer decimals for effect-type outputs.
#' @param digits_asn Integer decimals for sample-size outputs.
#'
#' @return A list with
#' \describe{
#'   \item{results}{Replicate-level data.frame.}
#'   \item{summary}{Numeric one-row data.frame.}
#'   \item{summary_pretty}{Character one-row data.frame with unified decimals.}
#' }
#' @export
simulate_trials_ce <- function(R,
                               n1, n2,
                               delta, piZ, theta, eta,
                               sigmaY, error_type = c("normal","t","skew"),
                               muX_C, muX_T = NULL, sigmaX,
                               gamma0, gamma1,
                               mu0,
                               lookup,
                               pi_fixed = 0.5,
                               digits_rate = 3,
                               digits_effect = 4,
                               digits_asn = 1) {
  error_type <- match.arg(error_type)
  
  # Run R replicates
  out <- replicate(
    R,
    simulate_trial_ce(
      n1, n2, delta, piZ, theta, eta, sigmaY, error_type,
      muX_C, muX_T, sigmaX, gamma0, gamma1, mu0, lookup, pi_fixed
    ),
    simplify = FALSE
  )
  
  DF <- do.call(rbind, lapply(out, as.data.frame))
  non_early <- subset(DF, early_stop == FALSE)
  
  # Overall coverage: early-stop -> surrogate LCL; otherwise final LCL
  coverage_overall <- mean(DF$coverage_overall %in% TRUE, na.rm = TRUE)
  
  # Truth for bias/MSE on non-early-stop replicates (design-fixed π)
  truth <- delta + eta * pi_fixed
  bias_val <- {
    x <- non_early$final_est - truth
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
  }
  mse_val <- {
    x <- non_early$final_est - truth
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x^2)
  }
  
  summary <- data.frame(
    rejection_rate  = mean((DF$reject %in% TRUE) | (DF$early_stop %in% TRUE), na.rm = TRUE),
    early_stop_rate = mean(DF$early_stop %in% TRUE, na.rm = TRUE),
    coverage        = coverage_overall,
    bias            = bias_val,
    mse             = mse_val,
    ASN_total       = mean(DF$sample_used_total, na.rm = TRUE),
    ASN_per_arm     = mean(DF$sample_used_per_arm, na.rm = TRUE)
  )
  
  # Pretty (fixed-decimal) display version
  fmt <- function(x, k) sprintf(paste0("%.", k, "f"), x)
  summary_pretty <- data.frame(
    rejection_rate  = fmt(summary$rejection_rate,  digits_rate),
    early_stop_rate = fmt(summary$early_stop_rate, digits_rate),
    coverage        = fmt(summary$coverage,        digits_rate),
    bias            = fmt(summary$bias,            digits_effect),
    mse             = fmt(summary$mse,             digits_effect),
    ASN_total       = fmt(summary$ASN_total,       digits_asn),
    ASN_per_arm     = fmt(summary$ASN_per_arm,     digits_asn),
    check.names = FALSE
  )
  
  list(results = DF, summary = summary, summary_pretty = summary_pretty)
}
