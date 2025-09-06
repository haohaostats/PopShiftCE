
# PopShiftCE

Monte Carlo Conditional Error (CE) for two-stage seamless trials **with
partial population shift**: surrogate-only interim (projected primary)
and full-data final analysis (OLS + HC covariance), **design-fixed**
mixture prevalence $\pi_Z^\dagger$.

- 📦 Install:

  ``` r
  # install.packages("devtools")
  devtools::install_github("YOUR_GH_USERNAME/PopShiftCE")
  ```

- 📖 Paper (preprint/Manuscript): *add link when available*  

- 🐛 Issues: <https://github.com/YOUR_GH_USERNAME/PopShiftCE/issues>

## Quick start

``` r
library(PopShiftCE)
set.seed(12345); RNGkind("L'Ecuyer-CMRG")

# Build CE lookup under H0 (design-fixed target)
lookup <- build_ce_lookup(
  n1=150, n2=150,
  muX_C=3.0, sigmaX=0.15, gamma0=2.0, gamma1=5.5,
  mu0=22, sigmaY=1.0, theta=1.0, eta=0.0, piZ=0.5,
  pi_fixed=0.5, error_type="normal", rho_XY=0.3,
  alpha_one_sided=0.05, B_ref=100000, batch_size=5000
)
lookup

# One OC point (delta = 0.3)
muX_T <- 3.0 + 0.3/5.5
res <- simulate_trials_ce(
  R = 1000,
  n1=150, n2=150,
  delta=0.3, piZ=0.5, theta=1, eta=0,
  sigmaY=1, error_type="normal",
  muX_C=3.0, muX_T=muX_T, sigmaX=0.15,
  gamma0=2.0, gamma1=5.5, mu0=22,
  lookup=lookup, pi_fixed=0.5
)
res$summary
```

## Methods (brief)

- Interim statistic:
  $T_1=\dfrac{\bar{\hat Y}_T-\bar{\hat Y}_C}{\sqrt{S^2_{\hat Y,T}/n_1 + S^2_{\hat Y,C}/n_1}}$,
  $\hat Y=\gamma_0+\gamma_1 X$ (external calibration treated as fixed).

- Final model:
  $Y=\mu_0+\delta\,I(T=1)+\theta Z + \eta\,I(T=1)Z + \epsilon$ (OLS +
  HC).

- Estimand: $\Delta_{\mathrm{marg}}^\dagger=\delta+\eta\,\pi_Z^\dagger$
  (design-fixed).

- MC-CE: calibrate single $b_{\mathrm{ref}}$ to attain one-sided
  $\alpha$; estimate conditional error $e(t_1)$ and derive $c(t_1)$ from
  the conditional null of the *actual* final statistic.

## Citation

``` r
citation("PopShiftCE")
```
