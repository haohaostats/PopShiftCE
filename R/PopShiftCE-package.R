#' PopShiftCE: Monte Carlo Conditional Error for Two-Stage Trials with Partial Population Shift
#'
#' Implements Monte Carlo-calibrated Conditional Error (CE) for two-stage seamless trials
#' with surrogate-only interim and full-data final analysis under partial population shift
#' (design-fixed mixture prevalence). Includes CE lookup calibration under the joint Stage-1 (X, Y)
#' dependence and a final OLS + heteroskedasticity-consistent (HC) covariance analysis
#' targeting the marginal estimand \eqn{\Delta_{marg}^\dagger = \delta + \eta * \pi_Z^\dagger}.
#'
#' @keywords internal
"_PACKAGE"
