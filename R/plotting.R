
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
                  title = "(A) Joint under H0: (T1, Sfinal) with reference boundary") +
    ggplot2::theme_minimal()
  
  # (B) e(t1)
  tgrid <- lookup$z1_grid
  e_hat <- lookup$e_fun(tgrid)
  pB <- ggplot2::ggplot(data.frame(t1 = tgrid, e = e_hat),
                        ggplot2::aes(x = t1, y = e)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    ggplot2::labs(x = expression(t[1]), y = expression(e(t[1])),
                  title = "(B) Conditional error function e(t1) under H0") +
    ggplot2::theme_minimal()
  
  # (C) c(t1)
  c_hat <- lookup$c_fun(tgrid)
  pC <- ggplot2::ggplot(data.frame(t1 = tgrid, c = c_hat),
                        ggplot2::aes(x = t1, y = c)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = lookup$b_ref, linetype = 2) +
    ggplot2::geom_hline(yintercept = lookup$b_ref, linetype = 3) +
    ggplot2::labs(x = expression(t[1]), y = expression(c(t[1])),
                  title = "(C) Conditional critical value c(t1) under H0") +
    ggplot2::theme_minimal()
  
  list(pA = pA, pB = pB, pC = pC)
}
