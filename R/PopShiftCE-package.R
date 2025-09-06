#' PopShiftCE: Monte Carlo Conditional Error for Two-Stage Trials with Partial Population Shift
#'
#' Implements Monte Carlo-calibrated Conditional Error (CE) for two-stage seamless trials
#' with surrogate-only interim and full-data final analysis under partial population shift
#' (design-fixed mixture prevalence). Includes CE lookup calibration under the joint Stage-1 (X, Y)
#' dependence and a final OLS + heteroskedasticity-consistent (HC) covariance analysis
#' targeting the marginal estimand Δ_marg^\u2020 = δ + η * π_Z^\u2020.
#'
#' @keywords internal
"_PACKAGE"
