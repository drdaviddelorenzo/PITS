# =============================================================================
#  PITS - Package-level documentation and dataset documentation
# =============================================================================


#' PITS: Power of Interrupted Time Series studies
#'
#' Tools for estimating the statistical power of Interrupted Time Series (ITS)
#' designs, with a focus on healthcare applications.
#'
#' @description
#' The package provides a complete workflow for ITS power analysis:
#'
#' **Step 1 - Estimate nuisance parameters from pre-intervention data**
#'
#' Use [estimate_its_params()] to estimate `baseline`, `sigma` (residual SD),
#' and `rho` (AR(1) autocorrelation) from a pre-intervention time series.
#' Individual functions [estimate_baseline()], [estimate_sigma()],
#' [estimate_rho()], and [estimate_trend()] are also available.
#'
#' **Step 2 - Calculate power via Monte Carlo simulation**
#'
#' Use [calculate_power()] (single site) or [calculate_power_multi()]
#' (multiple sites) to estimate the probability of detecting a specified
#' intervention effect for a given study design.
#'
#' **Step 3 - Optimise the design**
#'
#' Use [power_sweep()] to evaluate power across a range of post-intervention
#' durations, and [run_power_grid()] for full factorial sensitivity analyses.
#' Visualise results with [plot_power_curve()] and [plot_power_heatmap()].
#'
#' **Convenience wrappers**
#'
#' [run_its_power()] replicates the interactive experience of the original
#' `its_power_tool.R` script. [estimate_and_calculate()] chains parameter
#' estimation and power calculation in a single call.
#'
#' @section Key references:
#' - Lopez Bernal J, et al. (2017). Interrupted time series regression for
#'   the evaluation of public health interventions: a tutorial.
#'   *Int J Epidemiol* 46:348-355. \doi{10.1093/ije/dyw098}
#' - Zhang F, et al. (2011). Simulation-based power calculation for designing
#'   interrupted time series analyses of health policy interventions.
#'   *J Clin Epidemiol* 64:1252-1261.
#' - Wagner AK, et al. (2002). Segmented regression analysis of interrupted
#'   time series studies in medication use research.
#'   *J Clin Pharm Ther* 27:299-309.
#'
#' @docType package
#' @name PITS-package
#' @aliases PITS
"_PACKAGE"


# -----------------------------------------------------------------------------
#' Example pre-intervention case fatality rate data
#'
#' Monthly case fatality rate (CFR) observations from a hypothetical hospital
#' over a 24-month pre-intervention period. This dataset is used in the package
#' vignettes to illustrate the full PITS workflow: parameter estimation followed
#' by power calculation. It represents the motivating example from the
#' accompanying paper, in which a hospital evaluates whether a clinical decision
#' support system (CDSS) reduces case fatality rate and wishes to determine how
#' long post-intervention follow-up is needed to detect a meaningful effect.
#'
#' The pre-intervention CFR is approximately 3.6 per cent, with low residual
#' variability (sigma approximately 0.2) and moderate positive autocorrelation
#' (rho approximately 0.3 to 0.5), consistent with monthly hospital data.
#'
#' @format A data frame with 24 rows and 2 variables:
#' \describe{
#'   \item{time}{Integer. Sequential monthly time index, 1 to 24.}
#'   \item{outcome}{Numeric. Case fatality rate, expressed as a percentage
#'     (for example, 3.5 represents 3.5 per cent).}
#' }
#'
#' @source Simulated dataset generated for illustrative purposes.
#'
#' @examples
#' data("example_cfr_data")
#' head(example_cfr_data)
#' plot(example_cfr_data$time, example_cfr_data$outcome, type = "o",
#'      xlab = "Month", ylab = "CFR (per cent)",
#'      main = "Pre-intervention CFR")
#'
#' # Estimate parameters:
#' params <- estimate_its_params(example_cfr_data)
"example_cfr_data"
