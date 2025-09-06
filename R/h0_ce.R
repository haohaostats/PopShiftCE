
#' One H0 pair (T1, S_final) for CE calibration (design-fixed target)
#' @keywords internal
one_H0_pair <- function(n1, n2, muX_C, sigmaX, gamma0, gamma1,
                        mu0, sigmaY, theta, eta, piZ,
                        pi_fixed = 0.5,
                        error_type = c("normal","t","skew"),
                        rho_XY = 0.3) {
  
  error_type <- match.arg(error_type)
  
  # Under H0 of design-fixed target: delta0 = - eta * pi_fixed
  delta0 <- - eta * pi_fixed
  
  # Stage-2 Z draws (used later)
  Zc2 <- stats::rbinom(n2, 1, piZ)
  Zt2 <- stats::rbinom(n2, 1, piZ)
  
  # Joint (X, Y) for Stage-1 participants (Z = 0)
  gen_joint <- function(n, muX, muY, sigmaX, sigmaY, rho){
    Sigma <- matrix(c(sigmaX^2, rho*sigmaX*sigmaY,
                      rho*sigmaX*sigmaY, sigmaY^2), 2, 2)
    M <- MASS::mvrnorm(n, mu = c(muX, muY), Sigma = Sigma)
    list(X = M[,1], Y = M[,2])
  }
  S1C <- gen_joint(n1, muX = muX_C,                 muY = mu0,          sigmaX, sigmaY, rho_XY)
  S1T <- gen_joint(n1, muX = muX_C + delta0/gamma1, muY = mu0 + delta0, sigmaX, sigmaY, rho_XY)
  
  # T1 from surrogate projection
  Yhat_c <- gamma0 + gamma1 * S1C$X
  Yhat_t <- gamma0 + gamma1 * S1T$X
  mC <- mean(Yhat_c); mT <- mean(Yhat_t)
  vC <- stats::var(Yhat_c); vT <- stats::var(Yhat_t)
  se1 <- sqrt(vT/n1 + vC/n1)
  T1  <- (mT - mC) / se1
  
  # Stage-1 matured Y (Z=0) + Stage-2
  df1 <- data.frame(stage = "S1",
                    T = c(rep(0, n1),  rep(1, n1)),
                    Z = 0,
                    Y = c(S1C$Y,       S1T$Y),
                    is_new = 0)
  
  gen_err <- function(n){
    if (error_type == "normal") stats::rnorm(n, 0, sigmaY)
    else if (error_type == "t")  sigmaY * (stats::rt(n, df = 5) / sqrt(5/3))
    else                         sigmaY * (stats::rexp(n, 1) - 1)
  }
  Yc2 <- mu0 + theta*Zc2 + gen_err(n2)
  Yt2 <- mu0 + delta0 + theta*Zt2 + eta*Zt2 + gen_err(n2)
  df2 <- data.frame(stage = "S2",
                    T = c(rep(0, n2), rep(1, n2)),
                    Z = c(Zc2,       Zt2),
                    Y = c(Yc2,       Yt2),
                    is_new = 1)
  
  df  <- rbind(df1, df2)
  fa  <- final_analysis(df, pi_fixed = pi_fixed)
  c(z1 = T1, zf = fa$Zf)
}

#' Build CE lookup via H0 Monte Carlo (design-fixed target)
#'
#' Calibrates a single reference boundary b_ref such that the reference rule
#' \eqn{\{T_1 \ge b\} \cup \{T_1 < b, S_{\mathrm{final}}\ge b\}} matches a one-sided alpha.
#' Then estimates the conditional error \eqn{e(t_1)} and the conditional critical value
#' \eqn{c(t_1)} from the conditional null of \eqn{S_{\mathrm{final}}|T_1=t_1}.
#'
#' @param alpha_one_sided target one-sided alpha
#' @param B_ref Monte Carlo size for H0 pairing
#' @param z1_grid grid of t1 values for lookup
#' @param min_in_bin minimum conditional sample per grid point
#' @param h0 initial half-width for local neighborhoods in t1
#' @return an object of class 'ce_lookup' with elements b_ref, e_fun, c_fun, z1_grid, meta
#' @export
build_ce_lookup <- function(n1, n2,
                            muX_C, sigmaX, gamma0, gamma1,
                            mu0, sigmaY, theta, eta, piZ,
                            pi_fixed = 0.5,
                            error_type = c("normal","t","skew"),
                            rho_XY = 0.3,
                            alpha_one_sided = 0.05,
                            B_ref = 150000, batch_size = 5000,
                            z1_grid = seq(-3.5, 3.0, by = 0.05),
                            min_in_bin = 400, h0 = 0.05) {
  error_type <- match.arg(error_type)
  
  n_batches <- ceiling(B_ref / batch_size)
  H0_list <- lapply(seq_len(n_batches), function(bi) {
    k <- if (bi < n_batches) batch_size else (B_ref - batch_size*(n_batches-1))
    mat <- replicate(k, one_H0_pair(
      n1, n2, muX_C, sigmaX, gamma0, gamma1,
      mu0, sigmaY, theta, eta, piZ,
      pi_fixed = pi_fixed, error_type = error_type, rho_XY = rho_XY
    ))
    t(mat)
  })
  H0_mat <- do.call(rbind, H0_list)
  colnames(H0_mat) <- c("z1","zf")
  H0_dt <- as.data.frame(H0_mat)
  H0_dt <- H0_dt[is.finite(H0_dt$z1) & is.finite(H0_dt$zf), , drop = FALSE]
  
  alpha_hat <- function(b){
    if (nrow(H0_dt) == 0) return(0)
    cond <- (H0_dt$z1 >= b) | ((H0_dt$z1 < b) & (H0_dt$zf >= b))
    cond[is.na(cond)] <- FALSE
    mean(cond)
  }
  b_lo <- 0.5; b_hi <- 4.5
  while (alpha_hat(b_lo) < alpha_one_sided && b_lo > -1) b_lo <- b_lo - 0.2
  while (alpha_hat(b_hi) > alpha_one_sided && b_hi < 6)   b_hi <- b_hi + 0.2
  b_ref <- uniroot(function(b) alpha_hat(b) - alpha_one_sided,
                   lower = b_lo, upper = b_hi, tol = 1e-4)$root
  
  e_vec <- c_vec <- numeric(length(z1_grid))
  for (ii in seq_along(z1_grid)) {
    z1g <- z1_grid[ii]
    h <- h0
    idx <- which((H0_dt$z1 < b_ref) & (abs(H0_dt$z1 - z1g) <= h))
    while (length(idx) < min_in_bin && h < 0.6) {
      h <- h * 1.5
      idx <- which((H0_dt$z1 < b_ref) & (abs(H0_dt$z1 - z1g) <= h))
    }
    zf_cond <- H0_dt$zf[idx]
    zf_cond <- zf_cond[is.finite(zf_cond)]
    if (length(zf_cond) < 50) {
      zf_pool <- H0_dt$zf[H0_dt$z1 < b_ref]
      zf_cond <- zf_pool[is.finite(zf_pool)]
    }
    e_hat <- mean(zf_cond >= b_ref)
    e_hat <- min(max(e_hat, 1e-5), 1 - 1e-5)
    e_vec[ii] <- e_hat
    c_vec[ii] <- as.numeric(stats::quantile(zf_cond, probs = 1 - e_hat, names = FALSE, type = 7))
  }
  
  structure(list(
    b_ref   = b_ref,
    z1_grid = z1_grid,
    e_fun   = stats::approxfun(z1_grid, e_vec, rule = 2),
    c_fun   = stats::approxfun(z1_grid, c_vec, rule = 2),
    meta    = list(B_ref = B_ref, min_in_bin = min_in_bin, h0 = h0,
                   rho_XY = rho_XY, alpha = alpha_one_sided,
                   pi_fixed = pi_fixed, error_type = error_type)
  ), class = "ce_lookup")
}

#' @export
print.ce_lookup <- function(x, ...) {
  cat(sprintf("CE lookup: b_ref=%.4f; alpha=%.3f; B_ref=%d\n",
              x$b_ref, x$meta$alpha, x$meta$B_ref))
  invisible(x)
}
