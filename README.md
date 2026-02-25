# PopShiftCE

<p align="center">
  <b>Monte Carlo Conditional Error for Two-Stage Seamless Adaptive Trials under Partial Population Shift</b>
</p>

<p align="center">
  <a href="https://github.com/haohaostats/PopShiftCE"><img alt="GitHub repo" src="https://img.shields.io/badge/GitHub-PopShiftCE-black?logo=github"></a>
  <img alt="R" src="https://img.shields.io/badge/language-R-276DC3?logo=r">
  <img alt="License" src="https://img.shields.io/badge/license-MIT-green">
  <img alt="Method package" src="https://img.shields.io/badge/focus-method%20implementation-blue">
  <img alt="Monte Carlo CE" src="https://img.shields.io/badge/CE-Monte%20Carlo%20calibrated-purple">
</p>

---

**PopShiftCE** implements the method proposed in the paper:

> *Two-Stage Seamless Adaptive Trials under Partial Population Shift: Surrogate-Only Interim and Monte Carlo Conditional-Error Control*

PopShiftCE is a **method implementation package**. It is designed to help users:

- calibrate a Monte Carlo conditional-error (CE) lookup under the null,
- simulate two-stage seamless trials with a surrogate-only interim and a primary-endpoint final analysis,
- summarize operating characteristics (rejection probability, early stopping, coverage, ASN), and
- visualize CE mappings and decision geometry.

---

## ✨ What this README example is (and is not)

This README provides a **method implementation example** that walks through the main workflow of the PopShiftCE method.

It is **not**:

- a "minimum code that merely runs",
- a paper-table reproduction script, or
- a full simulation pipeline for every sensitivity analysis in the manuscript.

### ✅ Goal of this example
Show the **core method workflow** in a way that GitHub users can understand and adapt:

1. start from baseline defaults,
2. calibrate a CE lookup under the null,
3. inspect the calibrated lookup and (optionally) its independent validation,
4. simulate a null scenario and an alternative scenario,
5. interpret `rejection_rate` correctly by scenario,
6. generate diagnostic plots.

> [!NOTE]
> **Computation note:** The package is **serial by default** (no parallel backend required). This is intentional for robustness across user environments. For interactive use, reduce `B_ref`, `B_val`, and `R`; for manuscript-scale studies, increase them according to your precision target.

---

## 📦 Installation

### Install from a local source directory

```r
# install.packages("devtools")
devtools::install_local("/path/to/PopShiftCE")
```

### Install from GitHub (if hosted)

```r
# install.packages("remotes")
remotes::install_github("haohaostats/PopShiftCE")
```

---

## 🚀 Method implementation example

```r
library(PopShiftCE)

# 1) Baseline parameter template
cfg <- popshiftce_defaults()
str(cfg)

# 2) Calibrate a CE lookup under H0
#    Use demo-size Monte Carlo numbers here for runtime.

lookup <- build_ce_lookup(
  n1 = cfg$n1,
  n2 = cfg$n2,
  muX_C = cfg$muX_C,
  sigmaX = cfg$sigmaX,
  gamma0 = cfg$gamma0,
  gamma1 = cfg$gamma1,
  mu0 = cfg$mu0,
  sigmaY = cfg$sigmaY,
  theta = cfg$theta,
  eta = cfg$eta,
  piZ = 0.5,
  pi_target = cfg$pi_target,
  pi_fixed = cfg$pi_fixed,
  error_type = cfg$error_type,
  rho_XY = cfg$rho_XY,
  alpha_one_sided = cfg$alpha_one_sided,
  B_ref = 150000,     # demo size
  batch_size = 2000,
  do_validation = FALSE,
  B_val = 50000,      # demo size
  seed_lookup = 20250901,
  seed_val = 20250911,
  verbose = TRUE
)

print(lookup)

# Optional: independent null validation summary (if do_validation = TRUE)
# lookup$meta$validation
# This includes achieved null rejection probabilities on an independent sample:
# - achieved_alpha_ref (reference rule)
# - achieved_alpha_CE  (implemented CE rule)
# In the paper's terminology, these are achieved one-sided Type I error rates
# under the calibration null model.

# 3) Generate H0 pairs for CE mapping diagnostics (separate from lookup calibration)
H0_pairs <- simulate_h0_pairs(
  B = 10000,
  n1 = cfg$n1,
  n2 = cfg$n2,
  muX_C = cfg$muX_C,
  sigmaX = cfg$sigmaX,
  gamma0 = cfg$gamma0,
  gamma1 = cfg$gamma1,
  mu0 = cfg$mu0,
  sigmaY = cfg$sigmaY,
  theta = cfg$theta,
  eta = cfg$eta,
  piZ = 0.5,
  pi_target = cfg$pi_target,
  pi_fixed = cfg$pi_fixed,
  error_type = cfg$error_type,
  rho_XY = cfg$rho_XY,
  seed = 20250908
)

# Returns a named list of ggplot objects by default (pA, pB, pC, pD)
ce_plots <- plot_ce_mapping(H0_pairs, lookup)

# Or request a combined 2x2 panel (requires patchwork)
# ce_panel <- plot_ce_mapping(H0_pairs, lookup, combine = TRUE)
# ce_panel

# 4) Simulate a null scenario (delta = 0)
#    Here, rejection_rate is interpreted as the achieved one-sided Type I error rate.
res_null <- simulate_trials_ce(
  R = 10000,  # demo size
  n1 = cfg$n1,
  n2 = cfg$n2,
  delta = 0.0,
  piZ = 0.5,
  theta = cfg$theta,
  eta = cfg$eta,
  sigmaY = cfg$sigmaY,
  error_type = cfg$error_type,
  muX_C = cfg$muX_C,
  sigmaX = cfg$sigmaX,
  gamma0 = cfg$gamma0,
  gamma1 = cfg$gamma1,
  mu0 = cfg$mu0,
  lookup = lookup,
  pi_target = cfg$pi_target,
  pi_fixed = cfg$pi_fixed,
  verbose = TRUE
)

res_null$summary
res_null$summary_pretty

# 5) Simulate an alternative scenario (delta = 0.3)
#    Here, rejection_rate is interpreted as power.
res_alt <- simulate_trials_ce(
  R = 10000,  # demo size
  n1 = cfg$n1,
  n2 = cfg$n2,
  delta = 0.3,
  piZ = 0.5,
  theta = cfg$theta,
  eta = cfg$eta,
  sigmaY = cfg$sigmaY,
  error_type = cfg$error_type,
  muX_C = cfg$muX_C,
  sigmaX = cfg$sigmaX,
  gamma0 = cfg$gamma0,
  gamma1 = cfg$gamma1,
  mu0 = cfg$mu0,
  lookup = lookup,
  pi_target = cfg$pi_target,
  pi_fixed = cfg$pi_fixed,
  verbose = TRUE
)

res_alt$summary
res_alt$summary_pretty

# 6) Decision geometry plot under the alternative (requires ggplot2)
p_decision <- plot_decision_geometry(res_alt$results, lookup)
p_decision

# Optional: combined diagnostics panel (requires patchwork)
# panel <- plot_diagnostic_panel(res_null$results, res_alt$results, lookup)
# panel
```

---

## 📊 Interpreting simulation summaries (paper terminology)

`simulate_trials_ce()` returns:

- `results`: replicate-level results (useful for custom diagnostics/plots)
- `summary`: numeric summary
- `summary_pretty`: formatted character summary for display

### Key fields in `summary`

- `rejection_rate`: **overall rejection probability** (including early stops)
  - Under null scenarios: interpret as the **achieved one-sided Type I error rate**
  - Under alternative scenarios: interpret as **power**
- `early_stop_rate`: interim early-efficacy stopping probability
- `coverage_final_conditional`: final-stage coverage among trials that continue to Stage 2
- `coverage_overall`: overall one-sided coverage combining interim and final branches
- `ASN_total`, `ASN_per_arm`: average sample number (total / per arm)
- `bias_final_conditional`, `mse_final_conditional`: final-stage performance among non-early-stop trials

> [!TIP]
> This follows the paper's revised reporting convention: in mixed null/alternative settings, use **rejection probability** (or **rejection rate**) as the primary label and interpret it by scenario.

---

## 🎯 Design-fixed vs as-observed target

The following functions support both target definitions:

- `final_analysis()`
- `build_ce_lookup()`
- `simulate_trial_ce()` / `simulate_trials_ce()`

### Target options

- `pi_target = "fixed"` (design-fixed mixture target; requires `pi_fixed`)
- `pi_target = "asObservedS2"` (uses the realized Stage-2 subgroup prevalence)

For the paper's primary analysis setup, start with:

```r
cfg <- popshiftce_defaults()
cfg$pi_target   # "fixed"
cfg$pi_fixed    # 0.5
```

---

## ⏱️ Runtime guidance

- The package is **serial by default** (no parallel backend).
- For interactive exploration, use smaller Monte Carlo sizes (e.g., `B_ref = 5000-20000`, `R = 500-2000`).
- For manuscript-quality operating characteristics, use larger values consistent with your protocol and precision target.

---

## 📚 Citation

If you use **PopShiftCE** in methodological work or trial-design simulations, please cite the associated paper and the package repository.
