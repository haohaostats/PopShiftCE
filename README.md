
# PopShiftCE: Monte Carlo Conditional Error for Adaptive Trials 🛡️

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/haohaostats/PopShiftCE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haohaostats/PopShiftCE/actions/workflows/R-CMD-check.yaml)

**An R package to design and simulate two-stage seamless adaptive trials with partial population shifts, using a Monte Carlo-calibrated Conditional Error (CE) framework.**

This package is the official implementation for the methodology described in the paper:

> Chen, H. (2025). *Two-Stage Seamless Adaptive Trials under Partial Population Shift: Surrogate-Only Interim and Monte Carlo Conditional-Error Control*.

The `PopShiftCE` package provides functions to calibrate the trial design, simulate its operating characteristics (Type I error, power, ASN), and visualize the decision geometry.

---
## 🎯 Overview

Traditional adaptive designs face challenges when the enrolled population shifts between stages. This package implements a robust solution where:

* **Stage 1 (Interim) ⏳:** A decision (stop early for efficacy or continue) is made based on a short-term surrogate endpoint.
* **Stage 2 (Final) 🏁:** A full-data analysis is performed on the primary endpoint, accounting for the population shift.
* **Type I Error Control ✅:** The conditional error principle is used to strictly control the overall one-sided Type I error via the conditional error principle. The calibration is performed via intensive Monte Carlo simulation, making the method robust and adaptable to complex final statistics (e.g., using heteroskedasticity-consistent covariance).

---
## 🛠️ Installation

You can install the development version of `PopShiftCE` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("haohaostats/PopShiftCE")
```

---
## 🚀 Example Workflow

Here is a complete workflow to design a trial, check its operating characteristics under the null and an alternative, and visualize the results.

### **Step 1: Build the CE Lookup Table 🧮**

This is the core calibration step based on `H0` simulation. It may take a few moments to run.

```r
library(PopShiftCE)

# Use a fixed seed for reproducibility
set.seed(12345)
RNGkind("L'Ecuyer-CMRG")

lookup <- build_ce_lookup(
  n1 = 150, n2 = 150,
  muX_C = 3.0, sigmaX = 0.15, gamma0 = 2.0, gamma1 = 5.5,
  mu0 = 22, sigmaY = 1.0, theta = 1.0, eta = 0.0, piZ = 0.5,
  pi_fixed = 0.5, error_type = "normal", rho_XY = 0.3,
  alpha_one_sided = 0.05,
  B_ref = 50000, # Use a smaller B_ref for quick examples
  batch_size = 5000
)

print(lookup)
```
```
CE lookup: b_ref=1.9416; alpha=0.050; B_ref=50000
```

### **Step 2: Evaluate Operating Characteristics 📊**

#### **(A) Under the Null (Type I Error Control)**

We verify that the error rate is controlled at the nominal level (e.g., 0.05).

```r
set.seed(12345)
res0 <- simulate_trials_ce(
  R = 10000,
  n1 = 150, n2 = 150,
  delta = 0.0, piZ = 0.5, theta = 1.0, eta = 0.0,
  sigmaY = 1.0, error_type = "normal",
  muX_C = 3.0, sigmaX = 0.15,
  gamma0 = 2.0, gamma1 = 5.5, mu0 = 22,
  lookup = lookup, pi_fixed = 0.5
)

print(res0$summary_pretty)
```
```
  rejection_rate early_stop_rate coverage    bias    mse ASN_total ASN_per_arm
1          0.050           0.025    0.950 -0.0001 0.0090     592.4       296.2
```

#### **(B) Under an Alternative (Power)**

We check the trial's power to detect a true effect (e.g., `delta = 0.3`).

```r
set.seed(12345)
res1 <- simulate_trials_ce(
  R = 10000,
  n1 = 150, n2 = 150,
  delta = 0.3, piZ = 0.5, theta = 1.0, eta = 0.0,
  sigmaY = 1.0, error_type = "normal",
  muX_C = 3.0, sigmaX = 0.15,
  gamma0 = 2.0, gamma1 = 5.5, mu0 = 22,
  lookup = lookup, pi_fixed = 0.5
)

print(res1$summary_pretty)
```
```
  rejection_rate early_stop_rate coverage    bias    mse ASN_total ASN_per_arm
1          0.987           0.886    0.971 -0.0015 0.0089     334.1       167.1
```

### **Step 3: Visualize the Design 📈**

The `plot_diagnostic_panel()` function provides a comprehensive "four-in-one" view of the method's mechanics.

```r
# This single function creates a complete diagnostic panel
diagnostic_panel <- plot_diagnostic_panel(
  h0_results = res0$results,
  alt_results = res1$results,
  lookup = lookup
)

# Display the panel
print(diagnostic_panel)

# Save the panel as a vector graphic for best quality
# ggsave("figures/Fig_Diagnostic_Panel.svg", diagnostic_panel, width = 12, height = 10)
```

![PopShiftCE Figure: Methodology Diagnostics Panel for the Monte Carlo Calibration.](figures/Fig_Diagnostic_Panel.svg)

The panel visualizes:
* **(A) Joint Distribution:** The relationship between the interim (`T1`) and final (`S_final`) statistics under `H0`.
* **(B) Conditional Error:** The calculated conditional error function `e(t1)`.
* **(C) Critical Value:** The final critical value `c(t1)` as a function of the interim result.
* **(D) Decision Geometry:** How the calibrated boundaries perform on data simulated under the alternative hypothesis.

---
## Functions at a Glance

| Function | Description |
| :--- | :--- |
| `build_ce_lookup()` | Calibrates the reference boundary and conditional error functions. |
| `simulate_trial_ce()` | Simulates a single trial from start to finish. |
| `simulate_trials_ce()` | A wrapper to simulate many trials and summarize operating characteristics. |
| `plot_diagnostic_panel()` | Creates the main "four-in-one" diagnostic plot panel. |
| `plot_ce_mapping()` | Creates the three diagnostic plots for the CE calibration. |

---
## License

This package is licensed under the MIT License.

## Maintainer
- hao hao ([@haohaostats](https://github.com/haohaostats))