#' Simulate one two-stage trial using a calibrated CE lookup
#'
#' Stage 1 uses a surrogate-only interim. If `T1 >= b_ref`, the trial stops for
#' early efficacy and returns a surrogate-scale one-sided lower confidence limit.
#' Otherwise, the trial proceeds to Stage 2, analyzes the primary endpoint using
#' [final_analysis()], and applies the conditional critical value `c(T1)`.
#'
#' @param n1,n2 Per-arm Stage-1 / Stage-2 sample sizes.
#' @param delta Treatment main effect in the `Z=0` stratum.
#' @param piZ Stage-2 prevalence of `Z=1` in the simulation scenario.
#' @param theta Main effect of `Z`.
#' @param eta Treatment-by-`Z` interaction.
#' @param sigmaY Standard deviation scale of the primary outcome error.
#' @param error_type One of `"normal"`, `"t"`, or `"skew"`.
#' @param muX_C Stage-1 control-arm surrogate mean.
#' @param muX_T Stage-1 treatment-arm surrogate mean. If `NULL`, set internally to
#'   `muX_C + delta/gamma1`.
#' @param sigmaX Standard deviation of the surrogate.
#' @param gamma0,gamma1 External surrogate projection coefficients.
#' @param mu0 Baseline intercept for the primary outcome.
#' @param lookup A `ce_lookup` object from [build_ce_lookup()].
#' @param pi_target Target definition for the final analysis: `"fixed"` or
#'   `"asObservedS2"`.
#' @param pi_fixed Design-fixed prevalence used if `pi_target = "fixed"`.
#'
#' @details For coverage under `pi_target = "asObservedS2"`:
#' - if the trial proceeds to Stage 2, coverage is assessed against the realized
#'   `pi_hat` returned by [final_analysis()];
#' - if the trial stops early (before Stage 2 is observed), coverage is assessed
#'   against the scenario prevalence `piZ`.
#'
#' @return A list containing decision indicators, estimates, lower bounds, sample
#'   usage, and diagnostic quantities (`z1`, `zf`, `cz1`).
#' @export
simulate_trial_ce <- function(n1, n2,
                              delta, piZ, theta, eta,
                              sigmaY,
                              error_type = c("normal", "t", "skew"),
                              muX_C, muX_T = NULL, sigmaX,
                              gamma0, gamma1,
                              mu0,
                              lookup,
                              pi_target = c("fixed", "asObservedS2"),
                              pi_fixed = 0.5) {
  error_type <- match.arg(error_type)
  pi_target <- match.arg(pi_target)

  if (is.null(muX_T)) muX_T <- muX_C + delta / gamma1

  # Stage 1 (surrogate-only interim)
  s1 <- gen_stage1_and_T1(n1, muX_C, muX_T, sigmaX, gamma0, gamma1)
  z1 <- s1$T1
  est1 <- s1$est1
  se1 <- s1$se1

  # Early efficacy stop
  if (!is.na(z1) && z1 >= lookup$b_ref) {
    lcl1 <- lcl_interim(est1, se1, lookup$b_ref)
    true_delta_marg <- if (pi_target == "fixed") delta + eta * pi_fixed else delta + eta * piZ
    covg_sur <- is.finite(lcl1) && (true_delta_marg >= lcl1)

    return(list(
      reject = TRUE,
      early_stop = TRUE,
      final_est = NA_real_,
      final_se = NA_real_,
      final_pi_hat = NA_real_,
      lcl_final = NA_real_,
      lcl_interim = lcl1,
      coverage_final = NA,
      coverage_overall = covg_sur,
      sample_used_total = 2 * n1,
      sample_used_per_arm = n1,
      z1 = z1,
      zf = NA_real_,
      cz1 = lookup$b_ref
    ))
  }

  # Continue to Stage 2
  gen_err <- function(n) {
    if (error_type == "normal") {
      stats::rnorm(n, 0, sigmaY)
    } else if (error_type == "t") {
      sigmaY * (stats::rt(n, df = 5) / sqrt(5 / 3))
    } else {
      sigmaY * (stats::rexp(n, 1) - 1)
    }
  }

  # Stage-1 matured Y (Z=0)
  df1 <- data.frame(
    stage = "S1",
    T = c(rep(0, n1), rep(1, n1)),
    Z = 0,
    Y = c(mu0 + gen_err(n1), mu0 + delta + gen_err(n1)),
    is_new = 0
  )

  # Stage-2 new participants
  df2 <- gen_stage2_new(n2, mu0, delta, theta, eta, sigmaY, piZ, error_type)
  df <- rbind(df1, df2)

  fa <- final_analysis(df, pi_target = pi_target, pi_fixed = pi_fixed)
  Zf <- fa$Zf
  est <- fa$Delta_hat
  se <- fa$se

  cz1 <- lookup$c_fun(z1)
  reject <- is.finite(Zf) && is.finite(cz1) && (Zf >= cz1)
  lcl_f <- lcl_final(est, se, cz1)

  true_delta_marg <- if (pi_target == "fixed") delta + eta * pi_fixed else delta + eta * fa$pi_hat
  covg_final <- is.finite(lcl_f) && (true_delta_marg >= lcl_f)

  list(
    reject = reject,
    early_stop = FALSE,
    final_est = est,
    final_se = se,
    final_pi_hat = fa$pi_hat,
    lcl_final = lcl_f,
    lcl_interim = NA_real_,
    coverage_final = covg_final,
    coverage_overall = covg_final,
    sample_used_total = 2 * (n1 + n2),
    sample_used_per_arm = (n1 + n2),
    z1 = z1,
    zf = Zf,
    cz1 = cz1
  )
}

#' Simulate many trials and summarize operating characteristics
#'
#' This is a serial (non-parallel) simulation wrapper designed for robustness
#' across different user environments. Users can control runtime by adjusting
#' `R`, `B_ref`, `B_val`, and related Monte Carlo sizes.
#'
#'
#' In mixed null/alternative settings, the primary rate reported by this function is
#' `rejection_rate` (overall rejection probability, including early stops). Under
#' null scenarios, this is interpreted as the achieved one-sided Type I error rate;
#' under alternative scenarios, it is interpreted as power.
#' @param R Number of trial replicates.
#' @inheritParams simulate_trial_ce
#' @param digits_rate Decimal places for rate-type outputs in `summary_pretty`.
#' @param digits_effect Decimal places for effect-type outputs in `summary_pretty`.
#' @param digits_asn Decimal places for ASN outputs in `summary_pretty`.
#' @param verbose Logical; if `TRUE`, prints periodic progress messages.
#'
#' @return A list with:
#' - `results`: replicate-level data frame;
#' - `summary`: numeric one-row summary (including `rejection_rate`,
#'   `early_stop_rate`, coverage, bias/MSE, and ASN);
#' - `summary_pretty`: fixed-decimal character version of `summary`.
#' @export
simulate_trials_ce <- function(R,
                               n1, n2,
                               delta, piZ, theta, eta,
                               sigmaY,
                               error_type = c("normal", "t", "skew"),
                               muX_C, muX_T = NULL, sigmaX,
                               gamma0, gamma1,
                               mu0,
                               lookup,
                               pi_target = c("fixed", "asObservedS2"),
                               pi_fixed = 0.5,
                               digits_rate = 3,
                               digits_effect = 4,
                               digits_asn = 1,
                               verbose = FALSE,
                               progress_callback = NULL) {
  error_type <- match.arg(error_type)
  pi_target <- match.arg(pi_target)
  R <- as.integer(R)
  if (!is.finite(R) || R <= 0) stop("R must be a positive integer.")

  out <- vector("list", R)
  for (rr in seq_len(R)) {
    out[[rr]] <- simulate_trial_ce(
      n1 = n1, n2 = n2,
      delta = delta, piZ = piZ, theta = theta, eta = eta,
      sigmaY = sigmaY, error_type = error_type,
      muX_C = muX_C, muX_T = muX_T, sigmaX = sigmaX,
      gamma0 = gamma0, gamma1 = gamma1,
      mu0 = mu0,
      lookup = lookup,
      pi_target = pi_target, pi_fixed = pi_fixed
    )
    if (is.function(progress_callback) &&
        ((rr %% max(1L, floor(R / 100L))) == 0L || rr == R)) {
      try(progress_callback(rr, R), silent = TRUE)
    }
    if (isTRUE(verbose) && (rr %% max(1L, floor(R / 10))) == 0L) {
      message(sprintf("simulate_trials_ce: %d/%d replicates complete", rr, R))
    }
  }

  DF <- do.call(rbind, lapply(out, as.data.frame))
  non_early <- subset(DF, early_stop == FALSE)

  coverage_final_conditional <- mean(non_early$coverage_final %in% TRUE, na.rm = TRUE)
  coverage_overall <- mean(DF$coverage_overall %in% TRUE, na.rm = TRUE)

  truth <- if (pi_target == "fixed") {
    delta + eta * pi_fixed
  } else {
    # For as-observed target, a single scalar "truth" is not fixed across final-stage
    # replicates. We use the realized target `delta + eta * final_pi_hat` replicate-wise.
    NA_real_
  }

  if (pi_target == "fixed") {
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
  } else {
    # replicate-wise target uses final_pi_hat on non-early replicates
    targ <- non_early$final_pi_hat
    x <- non_early$final_est - (delta + eta * targ)
    x <- x[is.finite(x)]
    bias_val <- if (length(x) == 0L) NA_real_ else mean(x)
    mse_val  <- if (length(x) == 0L) NA_real_ else mean(x^2)
  }

  summary <- data.frame(
    rejection_rate = mean((DF$reject %in% TRUE) | (DF$early_stop %in% TRUE), na.rm = TRUE),
    early_stop_rate = mean(DF$early_stop %in% TRUE, na.rm = TRUE),
    coverage_final_conditional = coverage_final_conditional,
    coverage_overall = coverage_overall,
    bias_final_conditional = bias_val,
    mse_final_conditional = mse_val,
    ASN_total = mean(DF$sample_used_total, na.rm = TRUE),
    ASN_per_arm = mean(DF$sample_used_per_arm, na.rm = TRUE)
  )

  fmt <- function(x, k) ifelse(is.na(x), NA_character_, sprintf(paste0("%.", k, "f"), x))
  summary_pretty <- data.frame(
    rejection_rate = fmt(summary$rejection_rate, digits_rate),
    early_stop_rate = fmt(summary$early_stop_rate, digits_rate),
    coverage_final_conditional = fmt(summary$coverage_final_conditional, digits_rate),
    coverage_overall = fmt(summary$coverage_overall, digits_rate),
    bias_final_conditional = fmt(summary$bias_final_conditional, digits_effect),
    mse_final_conditional = fmt(summary$mse_final_conditional, digits_effect),
    ASN_total = fmt(summary$ASN_total, digits_asn),
    ASN_per_arm = fmt(summary$ASN_per_arm, digits_asn),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  list(results = DF, summary = summary, summary_pretty = summary_pretty)
}
