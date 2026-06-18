# =============================================================================
#  PITS - Parameter estimation functions
#
#  These functions extract the nuisance parameters needed for power
#  calculation from pre-intervention time-series data.
#
#  Individual functions:  estimate_baseline(), estimate_sigma(),
#                         estimate_rho(), estimate_trend()
#  Convenience wrapper:   estimate_its_params()
# =============================================================================


# -----------------------------------------------------------------------------
#' Estimate baseline outcome from pre-intervention data
#'
#' Fits a linear trend to the pre-intervention series and returns the intercept
#' at \eqn{t = 1}, which represents the baseline (mean) outcome at the start
#' of the observation period.
#'
#' @param outcome Numeric vector of outcome values in the pre-intervention
#'   period, ordered chronologically. Minimum 12 observations; 24+ recommended.
#' @param time Integer vector of time indices. Defaults to
#'   \code{seq_along(outcome)}.
#'
#' @return A single numeric value: the estimated baseline outcome.
#'
#' @details
#' The baseline is the OLS intercept from the linear model
#' \eqn{outcome_t = \beta_0 + \beta_1 \cdot t + \varepsilon_t}.
#' When the pre-intervention trend is close to zero (as expected for a
#' stable outcome), this is approximately equal to the mean of \code{outcome}.
#'
#' @examples
#' data("example_cfr_data")
#' estimate_baseline(example_cfr_data$outcome)
#'
#' @seealso [estimate_its_params()] for extracting all parameters at once.
#' @export
estimate_baseline <- function(outcome, time = seq_along(outcome)) {
  outcome <- .check_outcome(outcome)
  time    <- .check_time(time, length(outcome))
  fit <- stats::lm(outcome ~ time)
  unname(stats::coef(fit)[["(Intercept)"]])
}


# -----------------------------------------------------------------------------
#' Estimate residual standard deviation from pre-intervention data
#'
#' Fits a linear trend to the pre-intervention series and returns the standard
#' deviation of the residuals, used as the noise parameter \eqn{\sigma} in
#' ITS power calculations.
#'
#' @inheritParams estimate_baseline
#'
#' @return A single positive numeric value: the estimated residual SD
#'   (\eqn{\sigma}).
#'
#' @details
#' \eqn{\sigma} is the standard deviation of the detrended pre-intervention
#' series. It captures how much the outcome varies from one time point to the
#' next after accounting for any underlying trend. Larger \eqn{\sigma} means
#' noisier data and lower power.
#'
#' As a rough guide, \eqn{\sigma} is typically 10-20\% of the baseline for
#' monthly hospital rates. If your estimate is much larger, consider whether
#' the outcome series is stable or whether aggregation to a lower frequency
#' would reduce noise.
#'
#' @examples
#' data("example_cfr_data")
#' estimate_sigma(example_cfr_data$outcome)
#'
#' @seealso [estimate_its_params()] for extracting all parameters at once.
#' @export
estimate_sigma <- function(outcome, time = seq_along(outcome)) {
  outcome <- .check_outcome(outcome)
  time    <- .check_time(time, length(outcome))
  fit <- stats::lm(outcome ~ time)
  stats::sd(stats::residuals(fit))
}


# -----------------------------------------------------------------------------
#' Estimate AR(1) autocorrelation from pre-intervention data
#'
#' Estimates the first-order autocorrelation coefficient (\eqn{\rho}) from the
#' residuals of a linear trend fit to the pre-intervention series.
#'
#' @inheritParams estimate_baseline
#'
#' @return A single numeric value in \eqn{(-1, 1)}: the estimated AR(1)
#'   autocorrelation coefficient.
#'
#' @details
#' The estimate is the Pearson correlation between consecutive residuals:
#' \deqn{\hat{\rho} = \text{cor}(\hat{\varepsilon}_t,\, \hat{\varepsilon}_{t+1})}
#'
#' Positive autocorrelation is nearly universal in routine health data and
#' reduces effective sample size, lowering power relative to a naive
#' calculation that assumes independence.
#'
#' Typical ranges by aggregation frequency:
#' \itemize{
#'   \item Daily: 0.7-0.9
#'   \item Weekly: 0.5-0.8
#'   \item Monthly: 0.3-0.5
#'   \item Quarterly: 0.1-0.4
#' }
#'
#' If \code{outcome} has fewer than 3 observations, \code{NA} is returned with
#' a warning.
#'
#' @examples
#' data("example_cfr_data")
#' estimate_rho(example_cfr_data$outcome)
#'
#' @seealso [estimate_its_params()] for extracting all parameters at once.
#' @export
estimate_rho <- function(outcome, time = seq_along(outcome)) {
  outcome <- .check_outcome(outcome)
  time    <- .check_time(time, length(outcome))
  fit   <- stats::lm(outcome ~ time)
  resid <- stats::residuals(fit)
  n     <- length(resid)
  if (n < 3) {
    warning("estimate_rho: fewer than 3 residuals - returning NA.")
    return(NA_real_)
  }
  stats::cor(resid[-n], resid[-1])
}


# -----------------------------------------------------------------------------
#' Estimate pre-intervention trend from pre-intervention data
#'
#' Fits a linear trend to the pre-intervention series and returns the slope
#' (change per time unit). Ideally this should be close to zero, indicating
#' a stable pre-intervention period.
#'
#' @inheritParams estimate_baseline
#'
#' @return A single numeric value: the estimated trend (slope) per time unit.
#'
#' @details
#' A non-trivial pre-intervention trend (e.g. \eqn{|trend| > 0.05 \times baseline})
#' may indicate that the outcome was not stable before the intervention. This
#' can bias ITS estimates and should be discussed in study planning.
#'
#' @examples
#' data("example_cfr_data")
#' estimate_trend(example_cfr_data$outcome)
#'
#' @seealso [estimate_its_params()] for extracting all parameters at once.
#' @export
estimate_trend <- function(outcome, time = seq_along(outcome)) {
  outcome <- .check_outcome(outcome)
  time    <- .check_time(time, length(outcome))
  fit <- stats::lm(outcome ~ time)
  unname(stats::coef(fit)[["time"]])
}


# -----------------------------------------------------------------------------
#' Estimate all ITS nuisance parameters from pre-intervention data
#'
#' Convenience wrapper that estimates all four nuisance parameters needed for
#' ITS power calculation from a pre-intervention data frame or numeric vector.
#' The output can be passed directly to [calculate_power()] or
#' [run_its_power()].
#'
#' @param data Either:
#'   \itemize{
#'     \item A data frame with columns named by \code{time_col} and
#'           \code{outcome_col}; or
#'     \item A numeric vector of outcome values (in which case \code{time_col}
#'           is ignored and a sequential time index is used).
#'   }
#' @param outcome_col Character. Name of the outcome column when \code{data}
#'   is a data frame. Default \code{"outcome"}.
#' @param time_col Character. Name of the time column when \code{data} is a
#'   data frame. Default \code{"time"}.
#' @param verbose Logical. If \code{TRUE} (default), prints a formatted
#'   parameter summary and guidance to the console.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{n_pre}{Integer. Number of pre-intervention observations.}
#'     \item{baseline}{Numeric. Estimated baseline (model intercept).}
#'     \item{sigma}{Numeric. Estimated residual standard deviation.}
#'     \item{rho}{Numeric. Estimated AR(1) autocorrelation.}
#'     \item{trend_pre}{Numeric. Estimated pre-intervention trend (slope).}
#'     \item{method}{Character. Estimation method actually used: \code{"GLS-ML"}
#'       (the default) or \code{"OLS"} (the fallback, used only when GLS fails
#'       to converge).}
#'   }
#'
#' @details
#' All four nuisance parameters are estimated jointly by maximum likelihood
#' from a single generalised least squares model with AR(1) errors, fitted with
#' [nlme::gls()] and [nlme::corAR1()]:
#' \deqn{y_t = \beta_0 + \beta_1 t + \varepsilon_t, \quad
#'   \varepsilon_t = \rho\,\varepsilon_{t-1} + \nu_t}
#' This is the recommended approach because \eqn{\sigma} and \eqn{\rho} are
#' mutually dependent: ordinary least squares estimates \eqn{\sigma} without
#' accounting for autocorrelation (biased upward when \eqn{\rho > 0}) and then
#' estimates \eqn{\rho} from serially correlated residuals. Fitting both from
#' the same likelihood avoids this inconsistency and yields parameters that are
#' directly compatible with the simulation engine used by [calculate_power()].
#' If GLS fails to converge - rare, and usually only for very short series
#' (fewer than about 12 observations) - the function falls back to the OLS
#' two-step estimator and records this in the returned \code{method} element.
#'
#' The individual helpers [estimate_baseline()], [estimate_sigma()],
#' [estimate_rho()] and [estimate_trend()] use the simpler OLS two-step and are
#' provided for quick, single-parameter checks; for power calculation, prefer
#' the jointly estimated values returned here.
#'
#' This function does \strong{not} estimate \code{level_change} or
#' \code{slope_change}. Those are clinical hypotheses - the minimum effects
#' you would consider meaningful to detect - and must be set based on clinical
#' judgement or published evidence, not derived from your own data.
#'
#' @examples
#' data("example_cfr_data")
#' params <- estimate_its_params(example_cfr_data, verbose = TRUE)
#' str(params)
#'
#' # Use directly in calculate_power():
#' # calculate_power(
#' #   n_pre        = params$n_pre,
#' #   n_post       = 24,
#' #   baseline     = params$baseline,
#' #   level_change = -1.0,   # your clinical hypothesis
#' #   sigma        = params$sigma,
#' #   rho          = params$rho
#' # )
#'
#' @seealso [estimate_baseline()], [estimate_sigma()], [estimate_rho()],
#'   [estimate_trend()], [calculate_power()], [diagnose_params()]
#' @export
estimate_its_params <- function(data,
                                outcome_col = "outcome",
                                time_col    = "time",
                                verbose     = TRUE) {

  # Accept a plain numeric vector
  if (is.numeric(data)) {
    outcome  <- data
    time_idx <- seq_along(outcome)
  } else {
    data <- as.data.frame(data)
    if (!outcome_col %in% names(data)) {
      stop(sprintf(
        "Column '%s' not found. Available columns: %s",
        outcome_col, paste(names(data), collapse = ", ")
      ))
    }
    outcome <- data[[outcome_col]]
    time_idx <- if (time_col %in% names(data)) {
      .coerce_time(data[[time_col]])
    } else {
      seq_along(outcome)
    }
  }

  # Remove missing values
  keep <- !is.na(outcome)
  if (any(!keep)) {
    warning(sprintf(
      "estimate_its_params: %d missing values removed.", sum(!keep)
    ))
    outcome  <- outcome[keep]
    time_idx <- time_idx[keep]
  }

  n <- length(outcome)
  if (n < 12) {
    warning(sprintf(
      paste("Only %d observations available.",
            "Estimates will be unreliable; aim for >= 24."), n
    ))
  }

  # Joint maximum-likelihood estimation via GLS with AR(1) errors. This keeps
  # sigma and rho internally consistent and aligned with the simulation engine.
  df_pre <- data.frame(outcome = outcome, time_idx = as.numeric(time_idx))

  fit_gls <- tryCatch(
    nlme::gls(
      outcome ~ time_idx,
      data        = df_pre,
      correlation = nlme::corAR1(form = ~ time_idx),
      method      = "ML"
    ),
    error = function(e) NULL
  )

  if (!is.null(fit_gls)) {
    method    <- "GLS-ML"
    baseline  <- unname(stats::coef(fit_gls)[["(Intercept)"]])
    trend_pre <- unname(stats::coef(fit_gls)[["time_idx"]])
    sigma     <- fit_gls$sigma
    rho       <- unname(stats::coef(fit_gls$modelStruct$corStruct,
                                    unconstrained = FALSE))
  } else {
    warning("estimate_its_params: GLS did not converge; ",
            "falling back to OLS two-step estimates.")
    method    <- "OLS"
    baseline  <- estimate_baseline(outcome, time_idx)
    sigma     <- estimate_sigma(outcome, time_idx)
    rho       <- estimate_rho(outcome, time_idx)
    trend_pre <- estimate_trend(outcome, time_idx)
  }

  if (verbose) {
    .print_params(n, baseline, sigma, rho, trend_pre, method)
  }

  list(
    n_pre     = n,
    baseline  = baseline,
    sigma     = sigma,
    rho       = rho,
    trend_pre = trend_pre,
    method    = method
  )
}


# =============================================================================
#  Internal helpers
# =============================================================================

#' @keywords internal
.check_outcome <- function(x) {
  if (!is.numeric(x)) stop("'outcome' must be a numeric vector.")
  if (length(x) < 2)  stop("'outcome' must have at least 2 observations.")
  x
}

#' @keywords internal
.check_time <- function(t, n) {
  if (!is.numeric(t)) stop("'time' must be a numeric vector.")
  if (length(t) != n) stop("'time' and 'outcome' must have the same length.")
  t
}

#' @keywords internal
.coerce_time <- function(x) {
  # If dates, convert to sequential integer index
  if (inherits(x, c("Date", "POSIXct", "POSIXlt"))) {
    return(as.numeric(x - min(x)) + 1L)
  }
  as.numeric(x)
}

#' @keywords internal
.print_params <- function(n, baseline, sigma, rho, trend_pre,
                          method = "GLS-ML") {
  sep <- strrep("-", 60)
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("  PITS - ITS Parameter Estimation\n")
  cat(strrep("=", 60), "\n\n")
  cat(sprintf("  Estimation method:      %s\n", method))
  cat(sprintf("  Observations (n_pre):   %d\n", n))
  cat(sprintf("  Baseline:               %.4f\n", baseline))
  cat(sprintf("  Sigma (residual SD):    %.4f  (%.1f%% of baseline)\n",
              sigma, abs(sigma / baseline) * 100))
  cat(sprintf("  Rho (AR1):              %.4f\n", rho))
  cat(sprintf("  Pre-trend (per unit):   %.4f\n", trend_pre))
  cat("\n", sep, "\n")

  if (abs(trend_pre) > 0.05 * abs(baseline)) {
    cat("  NOTE: Pre-intervention trend is non-trivial (>5% of baseline).\n")
    cat("  Check whether the outcome was stable before the intervention.\n\n")
  }

  cat("  Copy these into calculate_power():\n\n")
  cat(sprintf("    n_pre        = %d\n", n))
  cat(sprintf("    baseline     = %.4f\n", baseline))
  cat(sprintf("    sigma        = %.4f\n", sigma))
  cat(sprintf("    rho          = %.4f\n", rho))
  cat( "    level_change = ???  # YOUR clinical hypothesis\n")
  cat( "    slope_change = 0    # set if testing slope\n")
  cat("\n", strrep("=", 60), "\n\n", sep = "")
}
