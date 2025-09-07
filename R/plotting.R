
#' Plot CE mapping (diagnostics): joint (T1, S_final), e(t1), c(t1)
#'
#' @param H0_dt data.frame with columns z1, zf (H0 pairs)
#' @param lookup object from [build_ce_lookup()]
#' @return a named list of ggplot objects: pA, pB, pC
#' @export
plot_ce_mapping <- function(H0_dt, lookup) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install.packages('ggplot2').")
  }
  stopifnot(all(c("z1", "zf") %in% names(H0_dt)))
  stopifnot(is.list(lookup), !is.null(lookup$b_ref), !is.null(lookup$z1_grid),
            is.function(lookup$e_fun), is.function(lookup$c_fun))
  
  # (A) Joint under H0: (T1, S_final)
  pA <- ggplot2::ggplot(H0_dt, ggplot2::aes(x = z1, y = zf)) +
    ggplot2::stat_bin2d(bins = 120) +
    ggplot2::geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    ggplot2::geom_hline(yintercept = lookup$b_ref, linetype = 2) +
    ggplot2::labs(x = expression(T[1]), y = expression(S[final]),
                  title = "Joint under H0: (T1, Sfinal) with reference boundary") +
    ggplot2::theme_minimal()
  
  # (B) e(t1)
  tgrid <- lookup$z1_grid
  e_hat <- lookup$e_fun(tgrid)
  pB <- ggplot2::ggplot(data.frame(t1 = tgrid, e = e_hat),
                        ggplot2::aes(x = t1, y = e)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    ggplot2::labs(x = expression(t[1]), y = expression(e(t[1])),
                  title = "Conditional error function e(t1) under H0") +
    ggplot2::theme_minimal()
  
  # (C) c(t1)
  c_hat <- lookup$c_fun(tgrid)
  pC <- ggplot2::ggplot(data.frame(t1 = tgrid, c = c_hat),
                        ggplot2::aes(x = t1, y = c)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    ggplot2::geom_hline(yintercept = lookup$b_ref, linetype = 3) +
    ggplot2::labs(x = expression(t[1]), y = expression(c(t[1])),
                  title = "Conditional critical value c(t1) under H0") +
    ggplot2::theme_minimal()
  
  list(pA = pA, pB = pB, pC = pC)
}

#' Plot the complete 4-in-1 methodology and diagnostics panel
#'
#' This function creates a 2x2 panel plot that visualizes the core mechanics
#' of the MC-CE method. It combines the three CE calibration plots (joint H0
#' distribution, e(t1), c(t1)) with the decision geometry under an alternative.
#'
#' @param h0_results The full results data.frame from a null scenario (delta=0)
#'   simulation, as returned by `simulate_trials_ce$results`.
#' @param alt_results The full results data.frame from an alternative scenario
#'   (delta > 0), as returned by `simulate_trials_ce$results`.
#' @param lookup The CE lookup object returned by `build_ce_lookup()`.
#' @return A composite ggplot object created by `patchwork`.
#' @export
#' @examples
#' \dontrun{
#'   # Assuming `lookup`, `res0`, and `res1` have been generated:
#'   diagnostic_panel <- plot_diagnostic_panel(
#'     h0_results = res0$results,
#'     alt_results = res1$results,
#'     lookup = lookup
#'   )
#'   diagnostic_panel
#' }
plot_diagnostic_panel <- function(h0_results, alt_results, lookup) {
  # --- 1. Generate the three CE mapping plots from H0 data ---
  H0_dt <- subset(
    h0_results,
    early_stop == FALSE & is.finite(z1) & is.finite(zf),
    select = c(z1, zf)
  )
  ce_plots <- plot_ce_mapping(H0_dt, lookup)
  
  # --- 2. Generate the decision geometry plot from alternative data ---
  p_decision <- ggplot2::ggplot(alt_results, ggplot2::aes(x = z1, y = zf, color = reject)) +
    ggplot2::geom_point(alpha = 0.35, na.rm = TRUE) +
    ggplot2::stat_function(fun = lookup$c_fun, n = 300, color = "black") +
    ggplot2::geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    ggplot2::labs(x = expression(T[1]),
                  y = expression(S[final]),
                  title = "Decision Geometry under Alternative") +
    ggplot2::theme_minimal() +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Rejected?"))
  
  # --- 3. Combine all four plots using patchwork ---
  final_panel <- (ce_plots$pA | ce_plots$pB) / (ce_plots$pC | p_decision)
  
  # --- 4. Add overall annotations ---
  final_panel_with_titles <- final_panel +
    patchwork::plot_annotation(
      title = 'PopShiftCE Method: Calibration & Decision Geometry',
      subtitle = 'Panels A, B, C from H0 simulation; Panel D from Alternative',
      tag_levels = 'A'
    )
  
  return(final_panel_with_titles)
}