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

## ✨ What this package is

**PopShiftCE** implements the core **PopShiftCE method** proposed in the paper:

> *Two-Stage Seamless Adaptive Trials under Partial Population Shift: Surrogate-Only Interim and Monte Carlo Conditional-Error Control*

### ✅ What PopShiftCE is for

This package helps users:

- calibrate a Monte Carlo conditional-error (CE) lookup under the null,
- run two-stage seamless trial simulations with a **surrogate-only interim**,
- perform a **primary-endpoint final analysis** under partial population shift,
- summarize operating characteristics (rejection probability, early stopping, coverage, ASN), and
- visualize CE mappings and decision geometry.

---

## 🧭 Quick navigation

- [Installation](#-installation)
- [Before you run: parameter guide (important)](#-before-you-run-parameter-guide-important)
- [Quick start (explicit values, first-time users)](#-quick-start-explicit-values-first-time-users)
- [Interpreting outputs](#-interpreting-outputs)
- [Plotting diagnostics](#-plotting-diagnostics)
- [Design-fixed vs as-observed target](#-design-fixed-vs-as-observed-target)
- [Shorter workflow using defaults template (advanced)](#-shorter-workflow-using-defaults-template-advanced)
- [Runtime tips](#-runtime-tips)
- [Reproducibility notes](#-reproducibility-notes)
- [Citation](#-citation)

---

## 📦 Installation

### Install from a local source directory

```r
# install.packages("devtools")
devtools::install_local("/path/to/PopShiftCE")
```

### Install from GitHub

```r
# install.packages("remotes")
remotes::install_github("haohaostats/PopShiftCE")
```

---

## 📘 Before you run: parameter guide (important)

If this is your **first time** using the package, read this section first.

A common source of confusion is that functions like `build_ce_lookup()` and `simulate_trials_ce()` expose many arguments. The easiest way to understand them is to group them by **what part of the method they control**.

> [!TIP]
> For your first run, keep most parameters at paper-like baseline values and change only:
> **`n1`, `n2`, `delta`, `piZ`, `B_ref`, `B_val`, `R`**.

---

### 1) Trial size parameters (usually the first ones to change)

- **`n1`**: Stage 1 per-arm sample size (before interim)
- **`n2`**: Stage 2 additional per-arm sample size (after interim)

📌 If a trial does **not** stop early, the total sample size across two arms is:

- `2 * (n1 + n2)`

---

### 2) Stage-1 surrogate / interim parameters

These control the **surrogate-only interim** statistic \(T_1\).

- **`muX_C`**: control-arm mean of Stage-1 surrogate \(X\)
- **`sigmaX`**: standard deviation of Stage-1 surrogate \(X\)
- **`gamma0`, `gamma1`**: external projection coefficients used to construct
  \[
  \widehat{Y} = \gamma_0 + \gamma_1 X
  \]
  The interim statistic is based on this projected surrogate quantity.

---

### 3) Final primary-endpoint model parameters

These control the final-stage primary outcome \(Y\) model used in simulation / calibration.

- **`mu0`**: baseline intercept for the primary endpoint
- **`sigmaY`**: standard deviation scale of the primary-endpoint error
- **`theta`**: main effect of subgroup indicator \(Z\)
- **`eta`**: treatment-by-\(Z\) interaction (population-shift interaction)
- **`error_type`**: error distribution for the primary endpoint
  - `"normal"`
  - `"t"` (heavy-tailed)
  - `"skew"`

---

### 4) Population-shift and estimand target parameters (most commonly confused)

These parameters are related but **not the same**:

- **`piZ`**: the **scenario prevalence** of the post-amendment subgroup (`Z = 1`) in **Stage 2 simulation**
- **`pi_target`**: the target estimand definition used in the final analysis
  - `"fixed"` = design-fixed target mixture
  - `"asObservedS2"` = observed Stage-2 mixture target
- **`pi_fixed`**: the fixed target prevalence (used only when `pi_target = "fixed"`)

✅ Intuition:

- `piZ` = the prevalence in the **simulated trial scenario**
- `pi_fixed` = the prevalence used in the **target estimand** when a fixed target is chosen

---

### 5) Effect-size parameter (simulation scenario)

- **`delta`**: treatment main effect in the \(Z = 0\) stratum

Typical README examples use:

- `delta = 0.0` for a null scenario (interpret rejection rate as achieved one-sided Type I error rate)
- `delta = 0.3` for an alternative scenario (interpret rejection rate as power)

---

### 6) Monte Carlo CE calibration parameters (runtime/precision controls)

These affect **speed** and **numerical precision**.

- **`alpha_one_sided`**: one-sided significance level (e.g., `0.05`)
- **`rho_XY`**: Stage-1 surrogate–primary dependence used for **H0 CE calibration**
- **`B_ref`**: Monte Carlo size for CE lookup calibration (larger = more stable, slower)
- **`B_val`**: Monte Carlo size for independent null validation (if `do_validation = TRUE`)
- **`batch_size`**: batch size used during calibration simulation
- **`seed_lookup`, `seed_val`**: RNG seeds for reproducibility

---

## 🔤 Paper notation ↔ R arguments (quick mapping)

- \(n_1\), \(n_2\) → `n1`, `n2`
- Stage-2 prevalence \(\pi_Z\) (scenario value) → `piZ`
- Design-fixed target prevalence \(\pi_Z^\dagger\) → `pi_fixed` (when `pi_target = "fixed"`)
- One-sided level \(\alpha\) → `alpha_one_sided`
- Stage-1 surrogate–primary dependence \(\rho_{XY}\) for CE calibration → `rho_XY`
- Treatment main effect in \(Z=0\) stratum \(\delta\) → `delta`
- Treatment-by-shift interaction \(\eta\) → `eta`

---

## 🚀 Quick start (explicit values, first-time users)

This example uses **explicit numeric values** (instead of `cfg$n1`, `cfg$n2`, etc.) so that first-time users can immediately see what each parameter is.

> [!NOTE]
> The Monte Carlo sizes below are **demo values** for faster runtime. Increase them for more stable / publication-grade results.

```r
library(PopShiftCE)

# =========================================
# 1) Build a CE lookup under H0 (demo-size MC)
# =========================================
lookup <- build_ce_lookup(
  # ---- trial sizes ----
  n1 = 150,
  n2 = 150,

  # ---- stage-1 surrogate / projection ----
  muX_C = 3.0,       # control-arm mean of surrogate X
  sigmaX = 0.15,     # SD of surrogate X
  gamma0 = 2.0,      # projection intercept: Yhat = gamma0 + gamma1 * X
  gamma1 = 5.5,      # projection slope

  # ---- primary endpoint model / calibration model ----
  mu0 = 22,          # baseline intercept for Y
  sigmaY = 1.0,      # SD scale of primary-outcome error
  theta = 1.0,       # main effect of subgroup Z
  eta = 0.0,         # treatment-by-Z interaction (set to 0 in this baseline demo)
  piZ = 0.5,         # Stage-2 prevalence of post-amendment subgroup (Z=1)

  # ---- estimand target ----
  pi_target = "fixed",
  pi_fixed = 0.5,

  # ---- calibration settings ----
  error_type = "normal",
  rho_XY = 0.3,      # Stage-1 surrogate-primary dependence used in H0 calibration
  alpha_one_sided = 0.05,

  # ---- Monte Carlo sizes (demo values) ----
  B_ref = 150000,
  batch_size = 2000,
  do_validation = FALSE,
  B_val = 50000,

  # ---- reproducibility ----
  seed_lookup = 20260101,
  seed_val = 20260101,
  verbose = TRUE
)

print(lookup)

# Optional: inspect independent null validation summary
# lookup$meta$validation
# Includes achieved null rejection probabilities (reference rule and CE rule)
# on an independent H0 validation sample.

# =========================================
# 2) Simulate a null scenario (delta = 0.0)
#    -> rejection_rate = achieved one-sided Type I error rate
# =========================================
res_null <- simulate_trials_ce(
  R = 10000,

  # trial design
  n1 = 150,
  n2 = 150,

  # scenario effect + population shift
  delta = 0.0,
  piZ = 0.5,
  theta = 1.0,
  eta = 0.0,

  # primary endpoint error model
  sigmaY = 1.0,
  error_type = "normal",

  # stage-1 surrogate / projection
  muX_C = 3.0,
  sigmaX = 0.15,
  gamma0 = 2.0,
  gamma1 = 5.5,

  # primary endpoint baseline
  mu0 = 22,

  # CE lookup and target definition
  lookup = lookup,
  pi_target = "fixed",
  pi_fixed = 0.5,

  verbose = TRUE
)

# =========================================
# 3) Simulate an alternative scenario (delta = 0.3)
#    -> rejection_rate = power
# =========================================
res_alt <- simulate_trials_ce(
  R = 2000,

  # trial design
  n1 = 150,
  n2 = 150,

  # scenario effect + population shift
  delta = 0.3,
  piZ = 0.5,
  theta = 1.0,
  eta = 0.0,

  # primary endpoint error model
  sigmaY = 1.0,
  error_type = "normal",

  # stage-1 surrogate / projection
  muX_C = 3.0,
  sigmaX = 0.15,
  gamma0 = 2.0,
  gamma1 = 5.5,

  # primary endpoint baseline
  mu0 = 22,

  # CE lookup and target definition
  lookup = lookup,
  pi_target = "fixed",
  pi_fixed = 0.5,

  verbose = TRUE
)

# Inspect summaries
res_null$summary
res_null$summary_pretty

res_alt$summary
res_alt$summary_pretty
```

---

## 📊 Interpreting outputs

`simulate_trials_ce()` returns a list with:

- **`results`**: replicate-level results (for custom diagnostics and plots)
- **`summary`**: numeric one-row summary
- **`summary_pretty`**: formatted summary (character columns for display)

### Key fields in `summary`

- **`rejection_rate`**: overall rejection probability (including early stops)
  - under a **null** scenario → interpret as the **achieved one-sided Type I error rate**
  - under an **alternative** scenario → interpret as **power**
- **`early_stop_rate`**: interim early-efficacy stopping probability
- **`coverage_final_conditional`**: final-stage coverage among trials that continue to Stage 2
- **`coverage_overall`**: overall one-sided coverage combining interim and final branches
- **`ASN_total`**, **`ASN_per_arm`**: average sample number (total / per arm)
- **`bias_final_conditional`**, **`mse_final_conditional`**: final-stage performance among non-early-stop trials

> [!TIP]
> In mixed null/alternative grids, use **“rejection probability”** (or `rejection_rate`) as the primary label, and interpret it by scenario.

---

## 🖼️ Plotting diagnostics

The package provides plotting helpers for CE diagnostics and decision geometry.

### A) CE mapping diagnostics from H0 pairs

```r
# install.packages("ggplot2")
library(ggplot2)

H0_pairs <- simulate_h0_pairs(
  B = 5000,
  n1 = 150,
  n2 = 150,
  muX_C = 3.0,
  sigmaX = 0.15,
  gamma0 = 2.0,
  gamma1 = 5.5,
  mu0 = 22,
  sigmaY = 1.0,
  theta = 1.0,
  eta = 0.0,
  piZ = 0.5,
  pi_target = "fixed",
  pi_fixed = 0.5,
  error_type = "normal",
  rho_XY = 0.3,
  seed = 20260103
)

ce_plots <- plot_ce_mapping(H0_pairs, lookup, combine = FALSE)
ce_plots$pA  # joint null cloud
ce_plots$pB  # e(t1)
ce_plots$pC  # c(t1)
ce_plots$pD  # marginal density of T1
```

### B) Decision geometry under an alternative scenario

```r
# install.packages("ggplot2")
library(ggplot2)

p_decision <- plot_decision_geometry(res_alt$results, lookup)
p_decision
```

### C) Combined diagnostic panel (paper-style workflow panel)

```r
# install.packages("patchwork")
library(patchwork)

panel <- plot_diagnostic_panel(res_null$results, res_alt$results, lookup)
panel
```

---

## 🎯 Design-fixed vs as-observed target

The following functions support both target definitions:

- `final_analysis()`
- `build_ce_lookup()`
- `simulate_trial_ce()`
- `simulate_trials_ce()`

### Option 1: design-fixed target (recommended starting point)

```r
pi_target = "fixed"
pi_fixed = 0.5
```

### Option 2: as-observed Stage-2 target

```r
pi_target = "asObservedS2"
# pi_fixed is ignored in this mode
```

📌 Practical advice:

- Start with **`pi_target = "fixed"`** for your first workflow.
- Use **`"asObservedS2"`** when you specifically want the realized Stage-2 prevalence target.

---

## 🧪 Shorter workflow using defaults template (advanced)

Once you understand the parameters, you can use the built-in defaults template to reduce typing.

```r
library(PopShiftCE)

cfg <- popshiftce_defaults()

# Change only what you need
cfg$n1 <- 100
cfg$n2 <- 200
cfg$eta <- 0.3
cfg$pi_target <- "fixed"
cfg$pi_fixed <- 0.5

# Convenience wrapper around build_ce_lookup()
lookup2 <- get_lookup(
  n1 = cfg$n1,
  n2 = cfg$n2,
  piZ = 0.5,
  base = cfg,
  B_ref = 20000,
  B_val = 5000,
  do_validation = TRUE,
  seed_lookup = 20260111,
  seed_val = 20260112,
  verbose = TRUE
)

print(lookup2)
```

> [!NOTE]
> The defaults-template style is convenient for repeated simulation studies, but the **explicit-value example above is recommended for first-time users**.

---

## ⚙️ Runtime tips

- The package is **serial by default** (no parallel backend required).
- This is intentional for robustness across user environments.

### For interactive exploration (fast demo)

Use smaller Monte Carlo sizes, for example:

- `B_ref = 5000 ~ 20000`
- `B_val = 2000 ~ 5000`
- `R = 500 ~ 2000`

### For higher-precision operating characteristics

Increase Monte Carlo sizes according to your precision target and available compute.

---

## 🔁 Reproducibility notes

For reproducible runs:

- set `seed_lookup` for CE lookup calibration,
- set `seed_val` for independent H0 validation,
- set your own `set.seed(...)` before repeated simulation calls if needed.

The package uses deterministic calibration procedures once the seeds and parameters are fixed.

---

## ❓Common first-time questions

### “Why are there so many parameters?”
Because the package exposes the main components of the method (surrogate interim, primary endpoint model, target estimand, and MC-CE calibration) so users can adapt it to different trial-design scenarios.

### “Do I need to understand every parameter before running the package?”
No. Start from the **explicit-value quick start** above. Then change only:
- `n1`, `n2`
- `delta`
- `piZ`
- `B_ref`, `B_val`, `R`

---

## 📚 Citation

If you use **PopShiftCE** in methodological work or trial-design simulations, please cite:

1. the associated paper (method description), and
2. the package repository (implementation).


---

## ❤️ Acknowledgement of scope

This package is intentionally focused on making the **PopShiftCE method easy to use and adapt**. It prioritizes a clean implementation workflow over reproducing every result from the manuscript.
