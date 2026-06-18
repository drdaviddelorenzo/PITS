# =============================================================================
#  PITS - S3 methods for result objects
#
#  print.pits_power_result()   - formatted console output for calculate_power()
#  summary.pits_power_result() - extended summary with p-value distribution
#  print.pits_sweep_result()   - formatted sweep table
# =============================================================================


# -----------------------------------------------------------------------------
#' @keywords internal
.make_power_result <- function(result) {
  class(result) <- c("pits_power_result", "list")
  result
}

#' @keywords internal
.make_sweep_result <- function(df) {
  class(df) <- c("pits_sweep_result", "data.frame")
  df
}


# -----------------------------------------------------------------------------
#' Print method for PITS power results
#'
#' Displays a formatted summary of the output from [calculate_power()] or
#' [calculate_power_multi()].
#'
#' @param x A \code{pits_power_result} object.
#' @param ... Ignored.
#' @return Invisibly \code{x}.
#' @export
print.pits_power_result <- function(x, ...) {
  sep <- strrep("=", 60)
  cat("\n", sep, "\n", sep = "")
  cat("  PITS - ITS Power Estimate\n")
  cat(sep, "\n\n")

  p <- x$params
  cat(sprintf("  Design:          %s\n",
              if (!is.null(x$n_sites)) sprintf("Multi-site (%d sites)", x$n_sites)
              else "Single site"))
  cat(sprintf("  Test:            %s\n",  p$test))
  cat(sprintf("  n_pre:           %s\n",  p$n_pre))
  cat(sprintf("  n_post:          %s\n",  p$n_post))
  cat(sprintf("  Baseline:        %.3f\n",p$baseline))
  cat(sprintf("  Level change:    %.3f\n",p$level_change))
  if (!is.null(p$slope_change) && p$slope_change != 0)
    cat(sprintf("  Slope change:    %.3f\n",p$slope_change))
  cat(sprintf("  Sigma:           %.3f\n",p$sigma))
  cat(sprintf("  Rho:             %.3f\n",p$rho))
  cat(sprintf("  Alpha:           %.3f\n",p$alpha))
  cat(sprintf("  N simulations:   %d\n",  x$n_sim))
  cat(sprintf("  N converged:     %d (%.1f%%)\n",
              x$n_converged, 100 * x$n_converged / x$n_sim))
  cat("\n", strrep("-", 60), "\n", sep = "")
  cat(sprintf("  POWER = %.1f%%  -  %s\n", x$power_pct, x$interpretation))
  cat(strrep("-", 60), "\n\n", sep = "")
  invisible(x)
}


# -----------------------------------------------------------------------------
#' Summary method for PITS power results
#'
#' Returns an extended summary including the p-value distribution across
#' Monte Carlo replications.
#'
#' @param object A \code{pits_power_result} object.
#' @param ... Ignored.
#' @return Invisibly, a list with \code{power}, \code{params}, and
#'   \code{pvalue_quantiles}.
#' @export
summary.pits_power_result <- function(object, ...) {
  print(object)

  p_vals <- object$p_values[!is.na(object$p_values)]
  if (length(p_vals) > 0) {
    cat("  P-value distribution (converged replications):\n")
    qs <- stats::quantile(p_vals, c(0.05, 0.25, 0.50, 0.75, 0.95))
    cat(sprintf("    5th pctile: %.4f\n",  qs[1]))
    cat(sprintf("    25th pctile: %.4f\n", qs[2]))
    cat(sprintf("    Median:      %.4f\n", qs[3]))
    cat(sprintf("    75th pctile: %.4f\n", qs[4]))
    cat(sprintf("    95th pctile: %.4f\n", qs[5]))
    cat(sprintf("    P < alpha:   %d / %d (%.1f%%)\n",
                sum(p_vals < object$params$alpha),
                length(p_vals),
                100 * mean(p_vals < object$params$alpha)))
    cat("\n")
  }

  invisible(list(
    power            = object$power,
    params           = object$params,
    pvalue_quantiles = if (length(p_vals) > 0) qs else NULL
  ))
}


# -----------------------------------------------------------------------------
#' Print method for PITS sweep results
#'
#' Displays a formatted table of power estimates across post-intervention
#' durations, as returned by [power_sweep()].
#'
#' @param x A \code{pits_sweep_result} data frame.
#' @param ... Ignored.
#' @return Invisibly \code{x}.
#' @export
print.pits_sweep_result <- function(x, ...) {
  sep <- strrep("=", 55)
  cat("\n", sep, "\n", sep = "")
  cat("  PITS - Design optimisation sweep\n")
  cat(sep, "\n")
  cat(sprintf("  %-8s  %8s  %s\n", "n_post", "Power", "Interpretation"))
  cat(strrep("-", 55), "\n", sep = "")
  for (i in seq_len(nrow(x))) {
    adequate <- x$power[i] >= 0.80
    cat(sprintf("  %-8d  %7.1f%%  %s%s\n",
                x$n_post[i],
                x$power_pct[i],
                x$interpretation[i],
                if (adequate) "  <-- adequate" else ""))
  }
  cat(sep, "\n\n", sep = "")
  invisible(x)
}
