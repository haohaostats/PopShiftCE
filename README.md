
# PopShiftCE: Monte Carlo Conditional Error for Adaptive Trials

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**An R package to design and simulate two-stage seamless adaptive trials with partial population shifts, using a Monte Carlo-calibrated Conditional Error (CE) framework.**

This package is the official implementation for the methodology described in the paper:

> Chen, H. & Grazian, C. (2025). *Two-Stage Seamless Adaptive Trials under Partial Population Shift: Surrogate-Only Interim and Monte Carlo Conditional-Error Control*. (Your Journal Here).

The `PopShiftCE` package provides functions to calibrate the trial design, simulate its operating characteristics (Type I error, power, ASN), and visualize the decision geometry.

## Overview

Traditional adaptive designs face challenges when the enrolled population shifts between stages. This package implements a robust solution where:

-   **Stage 1 (Interim):** A decision (stop early for efficacy or continue) is made based on a short-term surrogate endpoint.
-   **Stage 2 (Final):** A full-data analysis is performed on the primary endpoint, accounting for the population shift.
-   **Type I Error Control:** The conditional error principle is used to strictly control the family-wise error rate. The calibration is performed via intensive Monte Carlo simulation, making the method robust and adaptable to complex final statistics (e.g., using heteroskedasticity-consistent covariance).



## Installation

You can install the development version of `PopShiftCE` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("haohaostats/PopShiftCE")
```

## Example Workflow

Here is a complete workflow to design a trial, check its operating characteristics under the null and an alternative, and visualize the results.

```r
library(PopShiftCE)

# 1. Build the CE Lookup Table
# This is the core calibration step based on H0 simulation.
# It may take a few moments to run.
set.seed(12345)
RNGkind("L'Ecuyer-CMRG")

lookup <- build_ce_lookup(
  n1 = 150, n2 = 150,
  muX_C = 3.0, sigmaX = 0.15, gamma0 = 2.0, gamma1 = 5.5,
  mu0 = 22, sigmaY = 1.0, theta = 1.0, eta = 0.0, piZ = 0.5,
  pi_fixed = 0.5, error_type = "normal", rho_XY = 0.3,
  alpha_one_sided = 0.05,
  B_ref = 50000 # Use a smaller B_ref for quick examples
)

print(lookup)
#> CE lookup: b_ref=1.9548; alpha=0.050; B_ref=50000

# 2. Evaluate Operating Characteristics
# (A) Under the Null (Type I Error Control)
set.seed(12345)
res0 <- simulate_trials_ce(
  R = 10000,
  n1 = 150, n2 = 150, delta = 0.0, piZ = 0.5,
  # ... other parameters matching the lookup
  lookup = lookup
)
print(res0$summary)

# (B) Under an Alternative (Power)
set.seed(12345)
res1 <- simulate_trials_ce(
  R = 10000,
  n1 = 150, n2 = 150, delta = 0.3, piZ = 0.5,
  # ... other parameters matching the lookup
  lookup = lookup
)
print(res1$summary)

# 3. Visualize the Design
# (A) CE Mapping Diagnostics (using H0 results)
h0_pairs <- subset(res0$results, !early_stop, select = c(z1, zf))
plots_ce <- plot_ce_mapping(h0_pairs, lookup)
plots_ce$pC # Plot the conditional critical value c(t1)

# (B) Decision Geometry (using H1 results)
library(ggplot2)
df_alt <- res1$results
ggplot(df_alt, aes(x = z1, y = zf, color = reject)) +
  geom_point(alpha = 0.2, na.rm = TRUE) +
  stat_function(fun = lookup$c_fun, color = "black") +
  geom_vline(xintercept = lookup$b_ref, linetype = 2) +
  theme_minimal() +
  labs(title = "Decision Geometry: Reject if S_final >= c(T1)")
```

## Core Functions

-   `build_ce_lookup()`: Calibrates the reference boundary and conditional error functions.
-   `simulate_trial_ce()`: Simulates a single trial from start to finish.
-   `simulate_trials_ce()`: A wrapper to simulate many trials and summarize operating characteristics.
-   `plot_ce_mapping()`: Creates diagnostic plots for the CE calibration.

## License

This package is licensed under the MIT License.