#' Default parameter template for PopShiftCE examples
#'
#' Returns a named list of baseline parameters similar to the paper's baseline
#' setting. Users can modify selected entries and pass them into
#' `build_ce_lookup()` / `simulate_trial_ce()` / `simulate_trials_ce()`.
#'
#' @return A named list.
#' @export
popshiftce_defaults <- function() {
  list(
    n1 = 150,
    n2 = 150,
    sigmaX = 0.15,
    gamma0 = 2.0,
    gamma1 = 5.5,
    muX_C = 3.0,
    mu0 = 22,
    sigmaY = 1.0,
    theta = 1.0,
    eta = 0.0,
    error_type = "normal",
    pi_target = "fixed",
    pi_fixed = 0.5,
    rho_XY = 0.3,
    alpha_one_sided = 0.05
  )
}
