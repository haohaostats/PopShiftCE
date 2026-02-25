#' Plot CE mapping diagnostics
#'
#' Creates CE calibration diagnostic plots from a null sample of `(T1, S_final)`
#' pairs and a calibrated lookup object. These plots help visualize the calibrated
#' conditional error mapping and the reference-rule geometry underlying the
#' targeted one-sided null rejection probability.
#'
#' By default the function returns a named list of plots:
#' - `pA`: joint null cloud `(T1, S_final)`;
#' - `pB`: conditional error function `e(t1)`;
#' - `pC`: conditional critical value function `c(t1)`;
#' - `pD`: marginal density of `T1` under `H0`.
#'
#' If `combine = TRUE`, a patchwork composite (2x2 panel) is returned instead.
#'
#' @param H0_dt Data frame with columns `z1` and `zf`.
#' @param lookup A `ce_lookup` object from [build_ce_lookup()].
#' @param combine Logical; if `TRUE`, combine the four panels using `patchwork`.
#' @param bins Number of bins for the 2D histogram in panel A.
#'
#' @return A named list of `ggplot` objects, or a patchwork object if
#'   `combine = TRUE`.
#' @export
plot_ce_mapping <- function(H0_dt, lookup, combine = FALSE, bins = 120) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  stopifnot(all(c("z1", "zf") %in% names(H0_dt)))
  stopifnot(is.list(lookup), !is.null(lookup$b_ref), !is.null(lookup$z1_grid),
            is.function(lookup$e_fun), is.function(lookup$c_fun))

  gg <- asNamespace("ggplot2")

  pA <- gg$ggplot(H0_dt, gg$aes(x = z1, y = zf)) +
    gg$stat_bin2d(bins = bins) +
    gg$geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    gg$geom_hline(yintercept = lookup$b_ref, linetype = 2) +
    gg$labs(x = expression(T[1]), y = expression(S[final]),
            title = "Joint Distribution under H0: (T1, Sfinal)") +
    gg$theme_minimal()

  tgrid <- lookup$z1_grid
  pB <- gg$ggplot(data.frame(t1 = tgrid, e = lookup$e_fun(tgrid)), gg$aes(t1, e)) +
    gg$geom_line() +
    gg$geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    gg$labs(x = expression(t[1]), y = expression(e(t[1])),
            title = "Conditional Error Function e(t1)") +
    gg$theme_minimal()

  pC <- gg$ggplot(data.frame(t1 = tgrid, c = lookup$c_fun(tgrid)), gg$aes(t1, c)) +
    gg$geom_line() +
    gg$geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    gg$geom_hline(yintercept = lookup$b_ref, linetype = 3) +
    gg$labs(x = expression(t[1]), y = expression(c(t[1])),
            title = "Conditional Critical Value c(t1)") +
    gg$theme_minimal()

  pD <- gg$ggplot(H0_dt, gg$aes(x = z1)) +
    gg$geom_density(fill = "gray80", alpha = 0.7) +
    gg$geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    gg$labs(x = expression(T[1]), y = "Density",
            title = "Marginal Density of T1") +
    gg$theme_minimal()

  plots <- list(pA = pA, pB = pB, pC = pC, pD = pD)

  if (!isTRUE(combine)) return(plots)

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required when combine = TRUE.")
  }

  panel <- (pA + pB) / (pC + pD) +
    patchwork::plot_annotation(tag_levels = "A")
  panel
}

#' Plot decision geometry under an alternative scenario
#'
#' Uses replicate-level results from [simulate_trials_ce()] and overlays the
#' conditional critical value curve `c(t1)`.
#'
#' @param results Data frame, typically `simulate_trials_ce(...)$results`.
#' @param lookup A `ce_lookup` object.
#'
#' @return A `ggplot` object.
#' @export
plot_decision_geometry <- function(results, lookup) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  stopifnot(all(c("z1", "zf", "early_stop", "reject") %in% names(results)))

  df_non <- subset(results, early_stop == FALSE & is.finite(z1) & is.finite(zf))
  gg <- asNamespace("ggplot2")

  gg$ggplot(df_non, gg$aes(x = z1, y = zf, color = reject)) +
    gg$geom_point(alpha = 0.35, na.rm = TRUE) +
    gg$stat_function(fun = lookup$c_fun, n = 300, color = "black") +
    gg$geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    gg$labs(x = expression(T[1]), y = expression(S[final]),
            title = "Decision Geometry under Alternative") +
    gg$theme_minimal() +
    gg$guides(color = gg$guide_legend(title = "Rejected?"))
}

#' Four-panel diagnostic panel (CE mapping + decision geometry)
#'
#' Combines CE mapping plots from a null run with the decision geometry under an
#' alternative run.
#'
#' @param h0_results Replicate-level results from a null simulation.
#' @param alt_results Replicate-level results from an alternative simulation.
#' @param lookup A `ce_lookup` object.
#'
#' @return A patchwork composite object.
#' @export
plot_diagnostic_panel <- function(h0_results, alt_results, lookup) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for this plot.")
  }
  H0_dt <- subset(
    h0_results,
    early_stop == FALSE & is.finite(z1) & is.finite(zf),
    select = c(z1, zf)
  )
  ce_plots <- plot_ce_mapping(H0_dt, lookup, combine = FALSE)
  p_decision <- plot_decision_geometry(alt_results, lookup)

  (ce_plots$pA | ce_plots$pB) / (ce_plots$pC | p_decision) +
    patchwork::plot_annotation(
      title = "PopShiftCE Method: Calibration and Decision Geometry",
      subtitle = "Panels A-C from H0 simulation; panel D from alternative simulation",
      tag_levels = "A"
    )
}
