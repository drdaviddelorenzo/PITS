# =============================================================================
#  PITS - Utility functions and high-level wrappers
#
#  interpret_power()       - qualitative power label
#  validate_params()       - check a parameter list before simulation
#  simulate_predata()      - generate synthetic pre-intervention data
#  export_results()        - save results to CSV / text file
#  run_its_power()         - full single-site workflow (matches script UX)
#  estimate_and_calculate()- one-call: CSV -> power estimate
# =============================================================================


# -----------------------------------------------------------------------------
#' Interpret a power estimate qualitatively
#'
#' Converts a numeric power estimate to a short descriptive label using
#' conventional thresholds.
#'
#' @param power Numeric. Power estimate in \eqn{[0, 1]}.
#'
#' @return A character string: \code{"Adequate (>= 80\%)"}, \code{"Borderline
#'   (60-79\%)"}, or \code{"Underpowered (< 60\%)"}.
#'
#' @examples
#' interpret_power(0.85)
#' interpret_power(0.70)
#' interpret_power(0.45)
#'
#' @export
interpret_power <- function(power) {
  if (!is.numeric(power) || length(power) != 1)
    stop("'power' must be a single numeric value.")
  pct <- power * 100
  if (pct >= 80) {
    "Adequate (>= 80%)"
  } else if (pct >= 60) {
    "Borderline (60-79%)"
  } else {
    "Underpowered (< 60%)"
  }
}


# -----------------------------------------------------------------------------
#' Validate ITS parameter values before simulation
#'
#' Checks that a set of ITS parameters is internally consistent and within
#' plausible ranges. Issues warnings for unusual but permitted values.
#'
#' @param n_pre Integer. Pre-intervention time points.
#' @param n_post Integer. Post-intervention time points.
#' @param baseline Numeric. Baseline outcome.
#' @param level_change Numeric. Expected level change.
#' @param sigma Numeric. Residual SD.
#' @param rho Numeric. AR(1) autocorrelation.
#' @param alpha Numeric. Significance threshold.
#' @param n_sim Integer. Number of simulations.
#'
#' @return Invisibly \code{TRUE} if all checks pass; stops with an error on
#'   critical failures.
#'
#' @examples
#' validate_params(n_pre = 24, n_post = 24, baseline = 15,
#'                 level_change = -3, sigma = 2.5, rho = 0.4)
#'
#' @export
validate_params <- function(n_pre, n_post, baseline, level_change,
                            sigma, rho, alpha = 0.05, n_sim = 1000L) {

  errors   <- character(0)
  warnings <- character(0)

  if (!is.numeric(n_pre)  || n_pre  < 2)
    errors <- c(errors, "n_pre must be a positive integer >= 2.")
  if (!is.numeric(n_post) || n_post < 1)
    errors <- c(errors, "n_post must be a positive integer >= 1.")
  if (!is.numeric(sigma)  || sigma  <= 0)
    errors <- c(errors, "sigma must be positive.")
  if (!is.numeric(rho)    || abs(rho) >= 1)
    errors <- c(errors, "rho must be in (-1, 1).")
  if (!is.numeric(alpha)  || alpha <= 0 || alpha >= 1)
    errors <- c(errors, "alpha must be in (0, 1).")
  if (!is.numeric(n_sim)  || n_sim < 10)
    errors <- c(errors, "n_sim must be >= 10.")

  if (length(errors)) stop(paste(errors, collapse = "\n"))

  if (n_pre < 12)
    warnings <- c(warnings, "n_pre < 12: parameter estimates will be unreliable.")
  if (abs(level_change) > abs(baseline))
    warnings <- c(warnings,
      "level_change exceeds baseline in magnitude - check units.")
  if (abs(level_change) / abs(sigma) < 0.5)
    warnings <- c(warnings,
      "Standardised effect size (level_change/sigma) is < 0.5 - study may be underpowered.")
  if (abs(rho) >= 0.7)
    warnings <- c(warnings,
      "rho >= 0.7 (high autocorrelation). Consider aggregating to a lower frequency.")
  if (n_sim < 500)
    warnings <- c(warnings, "n_sim < 500: power estimate will be imprecise.")

  if (length(warnings)) {
    for (w in warnings) warning(w, call. = FALSE)
  }

  invisible(TRUE)
}


# -----------------------------------------------------------------------------
#' Generate synthetic pre-intervention data
#'
#' Simulates a pre-intervention time series with known parameters. Useful for
#' testing and for vignette examples when real data are unavailable.
#'
#' @param n Integer. Number of time points to generate. Default 24.
#' @param baseline Numeric. Mean outcome at \eqn{t = 1}. Default 15.
#' @param sigma Numeric. Residual standard deviation. Default 2.5.
#' @param rho Numeric. AR(1) autocorrelation. Default 0.4.
#' @param trend Numeric. Linear trend per time unit. Default 0.
#' @param seed Integer or \code{NULL}. Random seed. Default 42.
#'
#' @return A data frame with columns \code{time} and \code{outcome}.
#'
#' @examples
#' pre <- simulate_predata(n = 24, baseline = 15, sigma = 2.5, rho = 0.4)
#' plot(pre$time, pre$outcome, type = "o")
#'
#' @export
simulate_predata <- function(n        = 24L,
                             baseline = 15,
                             sigma    = 2.5,
                             rho      = 0.4,
                             trend    = 0,
                             seed     = 42L) {

  if (!is.null(seed)) set.seed(seed)

  mu <- baseline + trend * seq_len(n)
  innovation_sd <- sigma * sqrt(max(0, 1 - rho^2))
  eps <- numeric(n)
  for (i in seq_len(n)) {
    prev   <- if (i > 1L) eps[i - 1L] else 0
    eps[i] <- rho * prev + stats::rnorm(1L, 0, innovation_sd)
  }

  data.frame(time = seq_len(n), outcome = mu + eps)
}


# -----------------------------------------------------------------------------
#' Export PITS results to CSV and plain-text summary
#'
#' Writes the output of [calculate_power()] or [power_sweep()] to timestamped
#' files in the specified directory.
#'
#' @param result Output from [calculate_power()] or [power_sweep()].
#' @param dir Character. Output directory. Created if it does not exist.
#'   Default \code{"pits_output"}.
#' @param prefix Character. File name prefix. Default \code{"pits"}.
#'
#' @return Invisibly, a named character vector of file paths written.
#'
#' @examples
#' \dontrun{
#' result <- calculate_power(n_pre = 24, n_post = 24, baseline = 15,
#'                           level_change = -3, sigma = 2.5, rho = 0.4)
#' export_results(result, dir = "my_study")
#' }
#'
#' @export
export_results <- function(result, dir = "pits_output", prefix = "pits") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

  paths <- character(0)

  # Detect whether this is a sweep result (data.frame) or a single result
  if (is.data.frame(result)) {
    csv_path <- file.path(dir, sprintf("%s_sweep_%s.csv", prefix, ts))
    utils::write.csv(result, csv_path, row.names = FALSE)
    paths["sweep_csv"] <- csv_path
    cat(sprintf("Sweep results saved: %s\n", csv_path))

  } else if (is.list(result) && "p_values" %in% names(result)) {
    # p-value CSV
    pval_df <- data.frame(
      simulation  = seq_along(result$p_values),
      p_value     = round(result$p_values, 6),
      significant = as.integer(result$p_values < result$params$alpha),
      converged   = as.integer(!is.na(result$p_values))
    )
    csv_path <- file.path(dir, sprintf("%s_pvalues_%s.csv", prefix, ts))
    utils::write.csv(pval_df, csv_path, row.names = FALSE)
    paths["pvalues_csv"] <- csv_path

    # Text summary
    txt_lines <- c(
      "PITS - Power of an ITS",
      sprintf("Run date:        %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      strrep("-", 50),
      sprintf("Power:           %.1f%%", result$power_pct),
      sprintf("Interpretation:  %s", result$interpretation),
      sprintf("N simulations:   %d", result$n_sim),
      sprintf("N converged:     %d", result$n_converged),
      strrep("-", 50),
      "Parameters:",
      sprintf("  n_pre:         %s", result$params$n_pre),
      sprintf("  n_post:        %s", result$params$n_post),
      sprintf("  baseline:      %s", result$params$baseline),
      sprintf("  level_change:  %s", result$params$level_change),
      sprintf("  slope_change:  %s", result$params$slope_change),
      sprintf("  sigma:         %s", result$params$sigma),
      sprintf("  rho:           %s", result$params$rho),
      sprintf("  test:          %s", result$params$test),
      sprintf("  alpha:         %s", result$params$alpha),
      sprintf("  seed:          %s", result$params$seed)
    )
    txt_path <- file.path(dir, sprintf("%s_summary_%s.txt", prefix, ts))
    writeLines(txt_lines, txt_path)
    paths["summary_txt"] <- txt_path

    cat(sprintf("Results saved:\n  %s\n  %s\n", csv_path, txt_path))

  } else {
    stop("'result' must be output from calculate_power() or power_sweep().")
  }

  invisible(paths)
}


# -----------------------------------------------------------------------------
#' Full single-site ITS power workflow
#'
#' Convenience wrapper that replicates the interactive experience of
#' \code{its_power_tool.R}. Runs the power simulation and optionally a design
#' sweep, prints formatted output to the console, and optionally saves results.
#'
#' @param n_pre Integer. Pre-intervention time points.
#' @param n_post Integer. Post-intervention time points (primary design).
#' @param baseline Numeric. Baseline outcome.
#' @param level_change Numeric. Expected level change (your clinical hypothesis).
#' @param slope_change Numeric. Expected slope change. Default 0.
#' @param sigma Numeric. Residual SD.
#' @param rho Numeric. AR(1) autocorrelation.
#' @param test Character. Effect to test. Default \code{"level"}.
#' @param alpha Numeric. Significance threshold. Default 0.05.
#' @param n_sim Integer. Monte Carlo replications. Default 1000.
#' @param seed Integer. Random seed. Default 123.
#' @param sweep Logical. If \code{TRUE}, also run [power_sweep()].
#' @param sweep_post Integer vector. \code{n_post} values for the sweep.
#' @param save_output Logical. If \code{TRUE}, save results via
#'   [export_results()].
#' @param output_dir Character. Directory for saved files.
#'
#' @return Invisibly, a list with elements \code{result} (from
#'   [calculate_power()]) and, if \code{sweep = TRUE}, \code{sweep}
#'   (from [power_sweep()]).
#'
#' @examples
#' \dontrun{
#' run_its_power(
#'   n_pre = 24, n_post = 24,
#'   baseline = 15, level_change = -3,
#'   sigma = 2.5, rho = 0.4,
#'   sweep = TRUE
#' )
#' }
#'
#' @seealso [calculate_power()], [power_sweep()], [estimate_its_params()]
#' @export
run_its_power <- function(n_pre,
                          n_post,
                          baseline,
                          level_change,
                          slope_change = 0,
                          sigma,
                          rho,
                          test         = c("level", "slope", "both"),
                          alpha        = 0.05,
                          n_sim        = 1000L,
                          seed         = 123L,
                          sweep        = FALSE,
                          sweep_post   = c(6L, 12L, 18L, 24L, 30L, 36L),
                          save_output  = FALSE,
                          output_dir   = "pits_output") {

  test <- match.arg(test)
  validate_params(n_pre, n_post, baseline, level_change, sigma, rho,
                  alpha, n_sim)

  cat("\nRunning ITS power simulation (n_sim =", n_sim, ")...\n\n")
  result <- calculate_power(
    n_pre        = n_pre,
    n_post       = n_post,
    baseline     = baseline,
    level_change = level_change,
    slope_change = slope_change,
    sigma        = sigma,
    rho          = rho,
    test         = test,
    alpha        = alpha,
    n_sim        = n_sim,
    seed         = seed
  )
  print(result)

  sweep_result <- NULL
  if (sweep) {
    cat("Running design optimisation sweep...\n")
    sweep_result <- power_sweep(
      sweep_post   = sweep_post,
      n_pre        = n_pre,
      baseline     = baseline,
      level_change = level_change,
      slope_change = slope_change,
      sigma        = sigma,
      rho          = rho,
      test         = test,
      alpha        = alpha,
      n_sim        = n_sim,
      seed         = seed,
      verbose      = TRUE
    )
  }

  if (save_output) {
    export_results(result, dir = output_dir)
    if (!is.null(sweep_result))
      export_results(sweep_result, dir = output_dir)
  }

  invisible(list(result = result, sweep = sweep_result))
}


# -----------------------------------------------------------------------------
#' Estimate parameters and calculate power in one step
#'
#' Convenience function that chains [estimate_its_params()] and
#' [calculate_power()]. Supply a pre-intervention dataset, a
#' \code{level_change}, and a target \code{n_post}, and receive a power
#' estimate.
#'
#' @param data Pre-intervention data: a data frame with columns \code{time}
#'   and \code{outcome}, or a numeric vector of outcome values.
#' @param level_change Numeric. Minimum clinically meaningful effect size
#'   (your clinical hypothesis). This is \strong{not} estimated from data.
#' @param n_post Integer. Planned post-intervention follow-up duration.
#' @param outcome_col Character. Name of the outcome column. Default
#'   \code{"outcome"}.
#' @param time_col Character. Name of the time column. Default \code{"time"}.
#' @param verbose Logical. Print parameter and power summaries. Default
#'   \code{TRUE}.
#' @param ... Additional arguments passed to [calculate_power()].
#'
#' @return Output from [calculate_power()].
#'
#' @examples
#' data("example_cfr_data")
#' result <- estimate_and_calculate(
#'   data         = example_cfr_data,
#'   level_change = -1.0,
#'   n_post       = 24
#' )
#' result$power_pct
#'
#' @seealso [estimate_its_params()], [calculate_power()], [run_its_power()]
#' @export
estimate_and_calculate <- function(data,
                                   level_change,
                                   n_post,
                                   outcome_col = "outcome",
                                   time_col    = "time",
                                   verbose     = TRUE,
                                   ...) {

  params <- estimate_its_params(data,
                                outcome_col = outcome_col,
                                time_col    = time_col,
                                verbose     = verbose)

  calculate_power(
    n_pre        = params$n_pre,
    n_post       = n_post,
    baseline     = params$baseline,
    level_change = level_change,
    sigma        = params$sigma,
    rho          = params$rho,
    ...
  )
}


# =============================================================================
#  Internal helpers
# =============================================================================

#' @keywords internal
.print_power_result <- function(result) {
  sep <- strrep("=", 55)
  cat(sep, "\n")
  cat(sprintf("  Power:           %.1f%%\n", result$power_pct))
  cat(sprintf("  Interpretation:  %s\n",     result$interpretation))
  cat(sprintf("  N simulations:   %d\n",     result$n_sim))
  cat(sprintf("  N converged:     %d (%.1f%%)\n",
              result$n_converged,
              100 * result$n_converged / result$n_sim))
  cat(sep, "\n\n")
}
