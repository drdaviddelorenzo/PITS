# PITS News

## PITS 0.1.0 (2026-03-26)

### Initial release

**Parameter estimation**

* `estimate_its_params()` — all-in-one wrapper: estimates `baseline`,
  `sigma`, `rho`, and `trend_pre` jointly by maximum likelihood from a
  single GLS model with AR(1) errors (`nlme::gls` + `nlme::corAR1`), with
  an automatic OLS fallback if GLS fails to converge. The estimation
  method used is reported in the returned `method` element. This matches
  the simulation engine and avoids the upward bias in `sigma` (and the
  unreliable `rho`) that the OLS two-step produces under autocorrelation.
  Accepts a data frame or a plain numeric vector; handles date-indexed
  time columns and missing values with informative warnings.
* Individual functions: `estimate_baseline()`, `estimate_sigma()`,
  `estimate_rho()`, `estimate_trend()`. Each uses a simple OLS two-step
  and is intended for quick single-parameter checks; for power
  calculation, prefer the jointly estimated values from
  `estimate_its_params()`.

**Power simulation**

* `simulate_its_data()` — data-generating process for a single ITS site:
  Gaussian outcome with linear trend, step change, slope change, and AR(1)
  autocorrelated errors. Supports negative rho.
* `fit_its_model()` — segmented regression via `nlme::gls()` with AR(1)
  correction; returns p-value for `"level"`, `"slope"`, or `"both"` effects.
* `calculate_power()` — Monte Carlo power estimation for single-site ITS.
  Returns a `pits_power_result` object.
* `calculate_power_multi()` — power for multi-site designs via
  `nlme::lme()` with site random intercepts and per-site AR(1). Returns a
  `pits_power_result` object.
* `power_sweep()` — runs `calculate_power()` across a vector of `n_post`
  values; returns a `pits_sweep_result` data frame.
* `build_param_grid()` — constructs a factorial parameter grid for
  sensitivity analyses.
* `run_power_grid()` — applies `calculate_power()` to every row of a
  grid; appends power estimates to the grid data frame.

**Plots and diagnostics**

* `plot_power_curve()` — line plot of power vs `n_post`, with 80% target
  line and minimum adequate duration marker.
* `plot_power_heatmap()` — colour-coded power grid across two parameters,
  with cell labels and adequacy marker.
* `plot_its_example()` — simulated ITS series with fitted segmented
  regression overlaid; useful for figures in papers and protocols.
* `diagnose_params()` — 2×2 diagnostic panel: observed series + trend,
  residuals over time, Q-Q plot, and residual ACF.

**S3 methods**

* `print.pits_power_result()` — formatted console output for
  `calculate_power()` and `calculate_power_multi()` results.
* `summary.pits_power_result()` — extended output including p-value
  quantile distribution.
* `print.pits_sweep_result()` — formatted table with adequacy markers
  for `power_sweep()` output.

**Utilities and wrappers**

* `interpret_power()` — converts a numeric power estimate to
  `"Adequate (>= 80%)"`, `"Borderline (60-79%)"`, or
  `"Underpowered (< 60%)"`.
* `validate_params()` — pre-flight parameter checks with informative
  errors and warnings.
* `simulate_predata()` — generates synthetic pre-intervention data with
  known parameters; useful for testing and vignette examples.
* `export_results()` — saves `pits_power_result` or `pits_sweep_result`
  objects to timestamped CSV and plain-text summary files. The output
  directory (`dir`) must be supplied by the user; there is no default path,
  so nothing is ever written to the working directory or home filespace
  unless explicitly requested.
* `run_its_power()` — full single-site workflow with console output and
  optional sweep and file saving.
* `estimate_and_calculate()` — chains parameter estimation and power
  calculation in a single call.

**Data**

* `example_cfr_data` — 24 months of synthetic monthly case fatality rate
  data (`baseline` ≈ 14.7%, `sigma` ≈ 1.51, `rho` ≈ 0.37) for use in
  vignettes and examples.

**Vignettes**

* *PITS: Power of an ITS — CDSS/CFR worked example*: complete workflow
  from parameter estimation through design optimisation, sensitivity
  analysis, multi-site power, and export.
* *Estimating ITS parameters from pre-intervention data*: detailed guide
  to individual and all-in-one estimation functions, diagnostic
  interpretation, and handling of missing data, date columns, and short
  pre-periods.
