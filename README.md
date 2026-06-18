# PITS — Power of an Interrupted Time Series

[![R-CMD-check](https://github.com/drdaviddelorenzo/PITS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/drdaviddelorenzo/PITS/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D%204.0-brightgreen)](https://www.r-project.org/)

An R package for estimating the **statistical power of Interrupted Time Series (ITS) studies**, with a focus on healthcare applications.

Use PITS to answer:
> *"How many months of post-intervention follow-up do I need to reliably detect a meaningful effect in my ITS study?"*

---

## Installation

```r
# From GitHub (once published):
# install.packages("remotes")
# remotes::install_github("drdaviddelorenzo/PITS")

# From source:
install.packages("path/to/PITS", repos = NULL, type = "source")
```

**Dependency:** `nlme` (ships with base R). No other dependencies required.

---

## Workflow

```
Pre-intervention data (CSV or data frame)
             ↓
  estimate_its_params()
             ↓
  baseline, sigma, rho
             ↓
  + level_change  ← your clinical hypothesis
             ↓
  calculate_power() / power_sweep()
             ↓
  Choose n_post that achieves ≥ 80% power
```

---

## Quick start

```r
library(PITS)

# 1. Load your pre-intervention data (or use the built-in example)
data("example_cfr_data")

# 2. Estimate nuisance parameters
params <- estimate_its_params(example_cfr_data)
#   n_pre     = 24
#   baseline  = 15.0   (mean pre-intervention CFR, %)
#   sigma     = 2.15   (residual SD)
#   rho       = 0.22   (AR(1) autocorrelation)

# 3. Calculate power for a candidate design
result <- calculate_power(
  n_pre        = params$n_pre,
  n_post       = 30,             # planned follow-up months
  baseline     = params$baseline,
  level_change = -3,             # minimum clinically meaningful effect
  sigma        = params$sigma,
  rho          = params$rho,
  n_sim        = 1000,
  seed         = 123
)
result$power_pct      # e.g. 78.4%
result$interpretation # "Borderline (60-79%)"

# 4. Optimise with a sweep
sweep <- power_sweep(
  sweep_post   = c(12, 18, 24, 30, 36, 48),
  n_pre        = params$n_pre,
  baseline     = params$baseline,
  level_change = -3,
  sigma        = params$sigma,
  rho          = params$rho,
  n_sim        = 1000
)
plot_power_curve(sweep)

# 5. One-call shortcut
result <- estimate_and_calculate(
  data         = example_cfr_data,
  level_change = -3,
  n_post       = 30
)
```

---

## Function reference

### Parameter estimation

| Function | Description |
|---|---|
| `estimate_its_params(data)` | Estimate all nuisance parameters at once |
| `estimate_baseline(outcome)` | Baseline (intercept at *t* = 1) |
| `estimate_sigma(outcome)` | Residual standard deviation |
| `estimate_rho(outcome)` | AR(1) autocorrelation |
| `estimate_trend(outcome)` | Pre-intervention trend (slope) |

### Power simulation

| Function | Description |
|---|---|
| `calculate_power(...)` | Monte Carlo power — single site |
| `calculate_power_multi(sites, ...)` | Monte Carlo power — multiple sites |
| `power_sweep(sweep_post, ...)` | Power across a range of `n_post` values |
| `build_param_grid(...)` | Factorial parameter grid |
| `run_power_grid(grid, ...)` | Power across a full parameter grid |

### Plots and diagnostics

| Function | Description |
|---|---|
| `plot_power_curve(sweep)` | Line plot of power vs `n_post` |
| `plot_power_heatmap(grid)` | Colour grid of power over two parameters |
| `plot_its_example(...)` | Simulated ITS series with fitted model |
| `diagnose_params(data)` | 2×2 diagnostic panel for pre-intervention data |

### Utilities and wrappers

| Function | Description |
|---|---|
| `interpret_power(power)` | Qualitative label (Adequate / Borderline / Underpowered) |
| `validate_params(...)` | Check parameters before simulation |
| `simulate_predata(...)` | Generate synthetic pre-intervention data |
| `export_results(result)` | Save results to CSV and text file |
| `run_its_power(...)` | Full single-site workflow with console output |
| `estimate_and_calculate(data, ...)` | Parameter estimation + power in one call |

---

## Key parameters

| Parameter | What it is | How to set it |
|---|---|---|
| `n_pre` | Pre-intervention time points | Fixed — use your historical data |
| `n_post` | Post-intervention time points | **Design lever** — use `power_sweep()` to optimise |
| `baseline` | Mean pre-intervention outcome | From `estimate_its_params()` |
| `sigma` | Residual SD (noise) | From `estimate_its_params()` |
| `rho` | AR(1) autocorrelation | From `estimate_its_params()`; use 0.4 if unknown |
| `level_change` | Minimum meaningful step change | **Clinical hypothesis** — set based on expert judgement |
| `slope_change` | Minimum meaningful trend change | **Clinical hypothesis** — set to 0 if testing level only |
| `test` | Effect to test: `"level"`, `"slope"`, or `"both"` | Depends on intervention type |

---

## Multi-site designs

```r
sites <- list(
  list(name = "Hospital A", n_pre = 24, n_post = 30,
       baseline = 15, level_change = -3, slope_change = 0,
       sigma = 2.5, rho = 0.4),
  list(name = "Hospital B", n_pre = 24, n_post = 30,
       baseline = 18, level_change = -3, slope_change = 0,
       sigma = 3.0, rho = 0.4)
)
result <- calculate_power_multi(sites, n_sim = 1000, seed = 123)
result$power_pct
```

---

## References

- Lopez Bernal J, et al. (2017). Interrupted time series regression for the evaluation of public health interventions: a tutorial. *Int J Epidemiol* 46:348–355.
- Zhang F, et al. (2011). Simulation-based power calculation for designing interrupted time series analyses. *J Clin Epidemiol* 64:1252–1261.
- Wagner AK, et al. (2002). Segmented regression analysis of interrupted time series studies. *J Clin Pharm Ther* 27:299–309.

---

## Citing PITS

```bibtex
@software{PITS_2026,
  title  = {PITS: Power of Interrupted Time Series Studies},
  author = {de Lorenzo, David},
  year   = {2026},
  url    = {https://github.com/drdaviddelorenzo/PITS}
}
```

---

## Licence

MIT — see [LICENSE](LICENSE).
