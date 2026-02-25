#' PopShiftCE: Monte Carlo Conditional Error for Two-Stage Trials
#'
#' `PopShiftCE` implements a Monte Carlo calibrated conditional error (CE) method
#' for two-stage seamless adaptive trials with a surrogate-only interim and a
#' primary-endpoint final analysis under partial population shift.
#'
#' Core package capabilities include:
#' - CE lookup calibration under the null via `build_ce_lookup()`;
#' - final analysis for design-fixed or as-observed Stage-2 mixture targets via
#'   `final_analysis()`;
#' - trial simulation and operating-characteristic summaries via
#'   `simulate_trial_ce()` and `simulate_trials_ce()` (rejection probability,
#'   early stopping, coverage, and ASN);
#' - diagnostic plotting for CE mappings and decision geometry.
#'
#' The package is designed as a **method implementation package**. It focuses on
#' making the proposed method easy to use and understand, rather than providing a
#' full paper-reproduction workflow for all manuscript tables and figures.
#'
#' @keywords internal
"_PACKAGE"
