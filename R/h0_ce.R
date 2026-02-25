#' Generate one H0 pair `(T1, S_final)` for CE calibration
#'
#' This internal helper simulates the joint Stage-1 surrogate and matured primary
#' outcomes under the null, then computes the surrogate-based interim statistic
#' `T1` and the final-stage Wald statistic `S_final` under the chosen target.
#'
#' @keywords internal
one_H0_pair <- function(n1, n2, muX_C, sigmaX, gamma0, gamma1,
                        mu0, sigmaY, theta, eta, piZ,
                        pi_target = c("fixed", "asObservedS2"),
                        pi_fixed = 0.5,
                        error_type = c("normal", "t", "skew"),
                        rho_XY = 0.3) {
  pi_target <- match.arg(pi_target)
  error_type <- match.arg(error_type)

  # Stage-2 subgroup draws (also needed for as-observed calibration target)
  Zc2 <- stats::rbinom(n2, 1, piZ)
  Zt2 <- stats::rbinom(n2, 1, piZ)
  pi_hat_S2 <- mean(c(Zc2, Zt2))

  # Null alignment for the selected marginal target
  delta0 <- if (pi_target == "fixed") {
    -eta * pi_fixed
  } else {
    -eta * pi_hat_S2
  }

  gen_joint <- function(n, muX, muY, sigmaX, sigmaY, rho) {
    Sigma <- matrix(c(sigmaX^2, rho * sigmaX * sigmaY,
                      rho * sigmaX * sigmaY, sigmaY^2), 2, 2)
    M <- MASS::mvrnorm(n, mu = c(muX, muY), Sigma = Sigma)
    list(X = M[, 1], Y = M[, 2])
  }

  # Stage-1 joint (X, Y) under H0
  S1C <- gen_joint(n1, muX = muX_C,                 muY = mu0,
                   sigmaX = sigmaX, sigmaY = sigmaY, rho = rho_XY)
  S1T <- gen_joint(n1, muX = muX_C + delta0 / gamma1, muY = mu0 + delta0,
                   sigmaX = sigmaX, sigmaY = sigmaY, rho = rho_XY)

  # Interim statistic T1 from the surrogate projection only
  Yhat_c <- gamma0 + gamma1 * S1C$X
  Yhat_t <- gamma0 + gamma1 * S1T$X
  mC <- mean(Yhat_c); mT <- mean(Yhat_t)
  vC <- stats::var(Yhat_c); vT <- stats::var(Yhat_t)
  se1 <- sqrt(vT / n1 + vC / n1)
  T1  <- (mT - mC) / se1

  # Stage-1 matured Y (Z=0)
  df1 <- data.frame(
    stage = "S1",
    T = c(rep(0, n1), rep(1, n1)),
    Z = 0,
    Y = c(S1C$Y, S1T$Y),
    is_new = 0
  )

  gen_err <- function(n) {
    if (error_type == "normal") {
      stats::rnorm(n, 0, sigmaY)
    } else if (error_type == "t") {
      sigmaY * (stats::rt(n, df = 5) / sqrt(5 / 3))
    } else {
      sigmaY * (stats::rexp(n, 1) - 1)
    }
  }

  # Stage-2 outcomes under H0 aligned to the chosen marginal target
  Yc2 <- mu0 + theta * Zc2 + gen_err(n2)
  Yt2 <- mu0 + delta0 + theta * Zt2 + eta * Zt2 + gen_err(n2)
  df2 <- data.frame(
    stage = "S2",
    T = c(rep(0, n2), rep(1, n2)),
    Z = c(Zc2, Zt2),
    Y = c(Yc2, Yt2),
    is_new = 1
  )

  fa <- final_analysis(rbind(df1, df2), pi_target = pi_target, pi_fixed = pi_fixed)
  c(z1 = T1, zf = fa$Zf)
}

#' Simulate null `(T1, S_final)` pairs for CE diagnostics or plotting
#'
#' This helper is useful for plotting CE mappings (e.g., the joint null cloud and
#' the marginal density of `T1`). It uses the same `one_H0_pair()` generator as
#' `build_ce_lookup()`.
#'
#' @param B Number of null pairs to simulate.
#' @param ... Arguments passed to [one_H0_pair()].
#' @param seed Optional integer seed for reproducibility. If `NULL`, the current
#'   RNG state is used.
#'
#' @return A data.frame with columns `z1` and `zf`.
#' @export
simulate_h0_pairs <- function(B, ..., seed = NULL) {
  B <- as.integer(B)
  if (!is.finite(B) || B <= 0) stop("B must be a positive integer.")
  if (!is.null(seed)) set.seed(seed, kind = "L'Ecuyer-CMRG")

  mat <- replicate(B, one_H0_pair(...))
  out <- as.data.frame(t(mat))
  colnames(out) <- c("z1", "zf")
  out <- out[is.finite(out$z1) & is.finite(out$zf), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Build a CE lookup table via H0 Monte Carlo calibration
#'
#' Calibrates a reference boundary `b_ref` for the rule
#' `T1 >= b_ref` or (`T1 < b_ref` and `S_final >= b_ref`) at a one-sided level
#' `alpha_one_sided`, then estimates the conditional error function `e(t1)` and
#' conditional critical value function `c(t1)`.
#'
#' The default implementation uses kernel-smoothed estimation with adaptive
#' bandwidths based on a minimum effective sample size and optional monotonicity
#' constraints (isotonic regression) for numerical stability. A hard-window
#' fallback (`ce_smoother = "binWindow"`) is also available.
#'
#' Optionally, the function performs an independent null validation step and
#' stores independent-validation summaries of the achieved one-sided null rejection probability (i.e., achieved one-sided Type I error rate) in `lookup$meta$validation`.
#'
#' @param n1,n2 Per-arm Stage-1 / Stage-2 sample sizes.
#' @param muX_C,sigmaX,gamma0,gamma1 Stage-1 surrogate and projection parameters.
#' @param mu0,sigmaY,theta,eta,piZ Primary-outcome and shift parameters.
#' @param pi_target Target definition: `"fixed"` or `"asObservedS2"`.
#' @param pi_fixed Design-fixed prevalence used if `pi_target = "fixed"`.
#' @param error_type Error distribution (`"normal"`, `"t"`, `"skew"`).
#' @param rho_XY Stage-1 surrogate-primary dependence used for H0 calibration.
#' @param alpha_one_sided Target one-sided significance level for the reference rule calibration.
#' @param B_ref Number of null pairs for calibration.
#' @param batch_size Batch size used during serial simulation to limit memory use.
#' @param z1_grid Grid of interim values where the lookup is estimated.
#' @param ce_smoother Either `"kernel"` (default) or `"binWindow"`.
#' @param kernel Kernel for the kernel smoother (`"gaussian"` or `"tricube"`).
#' @param h0 Initial bandwidth.
#' @param h_max Maximum bandwidth.
#' @param min_neff Minimum effective sample size for local smoothing.
#' @param enforce_monotone Logical; if `TRUE`, enforces monotone `e(t1)` and
#'   monotone-decreasing `c(t1)` via isotonic regression.
#' @param seed_lookup Optional seed used to make the calibration sample
#'   reproducible. Deterministic batch-specific seeds are derived from this value.
#' @param do_validation Logical; if `TRUE`, runs an independent validation sample.
#' @param B_val Number of null pairs for validation if `do_validation = TRUE`.
#' @param seed_val Optional seed for the validation sample.
#' @param verbose Logical; if `TRUE`, prints progress summaries.
#'
#' @return An object of class `"ce_lookup"` with elements `b_ref`, `z1_grid`,
#'   `e_fun`, `c_fun`, and `meta`.
#' @export
build_ce_lookup <- function(n1, n2,
                            muX_C, sigmaX, gamma0, gamma1,
                            mu0, sigmaY, theta, eta, piZ,
                            pi_target = c("fixed", "asObservedS2"),
                            pi_fixed = 0.5,
                            error_type = c("normal", "t", "skew"),
                            rho_XY = 0.3,
                            alpha_one_sided = 0.05,
                            B_ref = 150000,
                            batch_size = 5000,
                            z1_grid = seq(-3.5, 3.0, by = 0.05),
                            ce_smoother = c("kernel", "binWindow"),
                            kernel = c("gaussian", "tricube"),
                            h0 = 0.15,
                            h_max = 1.5,
                            min_neff = 800,
                            enforce_monotone = TRUE,
                            seed_lookup = 20250901,
                            do_validation = TRUE,
                            B_val = 50000,
                            seed_val = 20250911,
                            verbose = TRUE) {
  pi_target <- match.arg(pi_target)
  error_type <- match.arg(error_type)
  ce_smoother <- match.arg(ce_smoother)
  kernel <- match.arg(kernel)

  B_ref <- suppressWarnings(as.integer(B_ref))
  batch_size <- suppressWarnings(as.integer(batch_size))
  B_val <- suppressWarnings(as.integer(B_val))
  if (!is.finite(B_ref) || B_ref <= 0) stop("B_ref must be a positive integer.")
  if (!is.finite(batch_size) || batch_size <= 0) stop("batch_size must be a positive integer.")
  if (!is.finite(B_val) || B_val < 0) stop("B_val must be a nonnegative integer.")

  # 1) Generate H0 sample for lookup (serial, deterministic by batch)
  if (!is.null(seed_lookup)) set.seed(seed_lookup, kind = "L'Ecuyer-CMRG")
  n_batches <- ceiling(B_ref / batch_size)
  batch_seeds <- if (!is.null(seed_lookup)) as.integer(seed_lookup) + seq_len(n_batches) else rep(NA_integer_, n_batches)

  H0_list <- vector("list", n_batches)
  for (bi in seq_len(n_batches)) {
    if (!is.na(batch_seeds[bi])) set.seed(batch_seeds[bi], kind = "L'Ecuyer-CMRG")
    k <- if (bi < n_batches) batch_size else (B_ref - batch_size * (n_batches - 1L))
    mat <- replicate(k, one_H0_pair(
      n1 = n1, n2 = n2,
      muX_C = muX_C, sigmaX = sigmaX, gamma0 = gamma0, gamma1 = gamma1,
      mu0 = mu0, sigmaY = sigmaY, theta = theta, eta = eta, piZ = piZ,
      pi_target = pi_target, pi_fixed = pi_fixed,
      error_type = error_type, rho_XY = rho_XY
    ))
    H0_list[[bi]] <- t(mat)
    if (isTRUE(verbose)) {
      message(sprintf("CE calibration batch %d/%d complete", bi, n_batches))
    }
  }

  H0_mat <- do.call(rbind, H0_list)
  colnames(H0_mat) <- c("z1", "zf")
  H0_dt <- as.data.frame(H0_mat)
  H0_dt <- H0_dt[is.finite(H0_dt$z1) & is.finite(H0_dt$zf), , drop = FALSE]

  # 2) Calibrate reference boundary b_ref
  alpha_hat <- function(b) {
    if (nrow(H0_dt) == 0L) return(0)
    cond <- (H0_dt$z1 >= b) | ((H0_dt$z1 < b) & (H0_dt$zf >= b))
    cond[is.na(cond)] <- FALSE
    mean(cond)
  }

  b_lo <- 0.5; b_hi <- 4.5
  while (alpha_hat(b_lo) < alpha_one_sided && b_lo > -1) b_lo <- b_lo - 0.2
  while (alpha_hat(b_hi) > alpha_one_sided && b_hi < 6)  b_hi <- b_hi + 0.2
  root <- stats::uniroot(function(b) alpha_hat(b) - alpha_one_sided,
                         lower = b_lo, upper = b_hi, tol = 1e-4)
  b_ref <- root$root

  if (isTRUE(verbose)) {
    message(sprintf("Calibrated b_ref = %.4f (reference-rule null rejection probability = %.4f)", b_ref, alpha_hat(b_ref)))
  }

  # 3) Estimate e(t1) and c(t1)
  idx_sub <- H0_dt$z1 < b_ref
  z1_sub <- H0_dt$z1[idx_sub]
  zf_sub <- H0_dt$zf[idx_sub]
  y_bin  <- as.numeric(zf_sub >= b_ref)

  e_vec <- c_vec <- numeric(length(z1_grid))
  h_used <- rep(NA_real_, length(z1_grid))

  if (ce_smoother == "kernel") {
    for (ii in seq_along(z1_grid)) {
      t1 <- z1_grid[ii]
      h  <- h0

      repeat {
        w <- .kernel_weights(z1_sub, t1, h, kernel = kernel)
        ne <- .neff(w)
        if (ne >= min_neff || h >= h_max) break
        h <- h * 1.35
      }
      h_used[ii] <- h

      sw <- sum(w)
      if (!is.finite(sw) || sw <= 0) {
        e_vec[ii] <- 0.5
        c_vec[ii] <- b_ref
        next
      }

      e_hat <- sum(w * y_bin) / sw
      e_hat <- min(max(e_hat, 1e-6), 1 - 1e-6)
      e_vec[ii] <- e_hat

      mu_hat <- sum(w * zf_sub) / sw
      v_hat  <- sum(w * (zf_sub - mu_hat)^2) / sw
      v_hat  <- max(v_hat, 1e-8)
      sd_hat <- sqrt(v_hat)

      c_vec[ii] <- mu_hat + sd_hat * stats::qnorm(1 - e_hat)
    }

    if (isTRUE(enforce_monotone)) {
      e_iso <- stats::isoreg(z1_grid, e_vec)
      e_vec <- pmin(pmax(e_iso$yf, 1e-6), 1 - 1e-6)

      c_iso <- stats::isoreg(z1_grid, -c_vec)
      c_vec <- -c_iso$yf
    }
  } else {
    # Legacy hard-window fallback
    min_in_bin <- 400
    h0_bin <- 0.05

    for (ii in seq_along(z1_grid)) {
      z1g <- z1_grid[ii]
      h <- h0_bin
      idx <- which((H0_dt$z1 < b_ref) & (abs(H0_dt$z1 - z1g) <= h))
      while (length(idx) < min_in_bin && h < 0.6) {
        h <- h * 1.5
        idx <- which((H0_dt$z1 < b_ref) & (abs(H0_dt$z1 - z1g) <= h))
      }
      h_used[ii] <- h
      zf_cond <- H0_dt$zf[idx]
      zf_cond <- zf_cond[is.finite(zf_cond)]
      if (length(zf_cond) < 50L) {
        zf_pool <- H0_dt$zf[H0_dt$z1 < b_ref]
        zf_cond <- zf_pool[is.finite(zf_pool)]
      }
      e_hat <- mean(zf_cond >= b_ref)
      e_hat <- min(max(e_hat, 1e-5), 1 - 1e-5)
      e_vec[ii] <- e_hat
      c_vec[ii] <- as.numeric(stats::quantile(zf_cond, probs = 1 - e_hat, names = FALSE, type = 7))
    }
  }

  e_fun <- stats::approxfun(z1_grid, e_vec, rule = 2)
  c_fun <- stats::approxfun(z1_grid, c_vec, rule = 2)

  # 4) Optional independent validation
  val_out <- NULL
  if (isTRUE(do_validation) && B_val > 0L) {
    if (!is.null(seed_val)) set.seed(seed_val, kind = "L'Ecuyer-CMRG")
    Vmat <- replicate(B_val, one_H0_pair(
      n1 = n1, n2 = n2,
      muX_C = muX_C, sigmaX = sigmaX, gamma0 = gamma0, gamma1 = gamma1,
      mu0 = mu0, sigmaY = sigmaY, theta = theta, eta = eta, piZ = piZ,
      pi_target = pi_target, pi_fixed = pi_fixed,
      error_type = error_type, rho_XY = rho_XY
    ))
    Vdt <- as.data.frame(t(Vmat))
    colnames(Vdt) <- c("z1", "zf")
    Vdt <- Vdt[is.finite(Vdt$z1) & is.finite(Vdt$zf), , drop = FALSE]

    rej_ref <- (Vdt$z1 >= b_ref) | ((Vdt$z1 < b_ref) & (Vdt$zf >= b_ref))
    rej_ref[is.na(rej_ref)] <- FALSE
    a_ref <- mean(rej_ref)
    se_ref <- sqrt(a_ref * (1 - a_ref) / nrow(Vdt))
    ci_ref <- c(a_ref - 1.96 * se_ref, a_ref + 1.96 * se_ref)

    rej_ce <- (Vdt$z1 >= b_ref) | ((Vdt$z1 < b_ref) & (Vdt$zf >= c_fun(Vdt$z1)))
    rej_ce[is.na(rej_ce)] <- FALSE
    a_ce <- mean(rej_ce)
    se_ce <- sqrt(a_ce * (1 - a_ce) / nrow(Vdt))
    ci_ce <- c(a_ce - 1.96 * se_ce, a_ce + 1.96 * se_ce)

    val_out <- list(
      B_val = nrow(Vdt),
      seed_val = seed_val,
      achieved_alpha_ref = a_ref,
      achieved_alpha_ref_CI = ci_ref,
      achieved_alpha_CE = a_ce,
      achieved_alpha_CE_CI = ci_ce
    )

    if (isTRUE(verbose)) {
      message(sprintf(
        paste0("Validation (independent null sample): ref-rule rejection probability = %.4f [%.4f, %.4f], ",
               "CE-rule rejection probability = %.4f [%.4f, %.4f]"),
        a_ref, ci_ref[1], ci_ref[2], a_ce, ci_ce[1], ci_ce[2]
      ))
    }
  }

  structure(list(
    b_ref   = b_ref,
    z1_grid = z1_grid,
    e_fun   = e_fun,
    c_fun   = c_fun,
    meta    = list(
      alpha = alpha_one_sided,
      B_ref = B_ref,
      batch_size = batch_size,
      pi_target = pi_target,
      pi_fixed = pi_fixed,
      rho_XY = rho_XY,
      error_type = error_type,
      ce_smoother = ce_smoother,
      kernel = kernel,
      h0 = h0,
      h_max = h_max,
      min_neff = min_neff,
      enforce_monotone = enforce_monotone,
      h_used = h_used,
      seed_lookup = seed_lookup,
      validation = val_out
    )
  ), class = "ce_lookup")
}

#' Convenience wrapper for `build_ce_lookup()` using a defaults list
#'
#' @param n1,n2 Per-arm Stage-1 / Stage-2 sample sizes.
#' @param piZ Stage-2 prevalence for `Z=1`.
#' @param base A list, typically from [popshiftce_defaults()].
#' @param alpha_one_sided One-sided significance level used in CE lookup calibration.
#' @param ... Additional arguments passed to [build_ce_lookup()].
#'
#' @return A `ce_lookup` object.
#' @export
get_lookup <- function(n1, n2, piZ, base = popshiftce_defaults(),
                       alpha_one_sided = if (!is.null(base$alpha_one_sided)) base$alpha_one_sided else 0.05,
                       ...) {
  build_ce_lookup(
    n1 = n1, n2 = n2,
    muX_C = base$muX_C, sigmaX = base$sigmaX, gamma0 = base$gamma0, gamma1 = base$gamma1,
    mu0 = base$mu0, sigmaY = base$sigmaY, theta = base$theta, eta = base$eta, piZ = piZ,
    pi_target = base$pi_target, pi_fixed = base$pi_fixed,
    error_type = base$error_type,
    rho_XY = base$rho_XY,
    alpha_one_sided = alpha_one_sided,
    ...
  )
}

#' @export
print.ce_lookup <- function(x, ...) {
  cat(sprintf("CE lookup: b_ref=%.4f; one-sided level=%.3f; B_ref=%d; smoother=%s\n",
              x$b_ref,
              ifelse(is.null(x$meta$alpha), NA_real_, x$meta$alpha),
              ifelse(is.null(x$meta$B_ref), NA_integer_, x$meta$B_ref),
              ifelse(is.null(x$meta$ce_smoother), "?", x$meta$ce_smoother)))
  if (!is.null(x$meta$validation)) {
    v <- x$meta$validation
    cat(sprintf("Validation (CE-rule null rejection probability): %.4f [%0.4f, %0.4f] (B_val=%d)\n",
                v$achieved_alpha_CE, v$achieved_alpha_CE_CI[1], v$achieved_alpha_CE_CI[2], v$B_val))
  }
  invisible(x)
}
