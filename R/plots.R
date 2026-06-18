# =============================================================================
#  PITS - Plotting and diagnostic functions
#
#  plot_power_curve()   - power vs n_post line plot
#  plot_power_heatmap() - power across two parameters
#  plot_its_example()   - simulated ITS series with breakpoint
#  diagnose_params()    - diagnostic plots for pre-intervention data
# =============================================================================


# -----------------------------------------------------------------------------
#' Plot power as a function of post-intervention duration
#'
#' Produces a line plot of estimated power against \code{n_post}, as returned
#' by [power_sweep()]. Highlights the 80\% threshold and indicates the minimum
#' adequate duration.
#'
#' @param sweep_result A data frame as returned by [power_sweep()], with
#'   columns \code{n_post} and \code{power_pct}.
#' @param target Numeric. Power target line (0-100). Default 80.
#' @param xlab Character. x-axis label. Default
#'   \code{"Post-intervention time points (n_post)"}.
#' @param ylab Character. y-axis label. Default \code{"Estimated power (\%)"}.
#' @param main Character. Plot title. Default \code{"ITS power curve"}.
#' @param col Character. Line colour. Default \code{"steelblue"}.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @return Invisibly \code{NULL}. Called for its side-effect (plot).
#'
#' @examples
#' sweep <- power_sweep(
#'   sweep_post = c(12, 18, 24, 30, 36),
#'   n_pre = 24, baseline = 15, level_change = -3,
#'   sigma = 2.5, rho = 0.4, n_sim = 200, seed = 42
#' )
#' plot_power_curve(sweep)
#'
#' @seealso [power_sweep()], [plot_power_heatmap()]
#' @export
plot_power_curve <- function(sweep_result,
                             target = 80,
                             xlab   = "Post-intervention time points (n_post)",
                             ylab   = "Estimated power (%)",
                             main   = "ITS power curve",
                             col    = "steelblue",
                             ...) {

  if (!all(c("n_post", "power_pct") %in% names(sweep_result)))
    stop("'sweep_result' must have columns 'n_post' and 'power_pct'.")

  x <- sweep_result$n_post
  y <- sweep_result$power_pct

  graphics::plot(x, y,
       type = "o", pch = 16, lwd = 2, col = col,
       xlim = range(x), ylim = c(0, 100),
       xlab = xlab, ylab = ylab, main = main,
       las  = 1, bty = "l", ...)

  # Target line
  graphics::abline(h = target, lty = 2, col = "firebrick", lwd = 1.5)
  graphics::text(x = min(x), y = target + 3,
                 labels = sprintf("%d%% target", target),
                 col = "firebrick", adj = 0, cex = 0.85)

  # Shade adequate region
  graphics::rect(xleft   = min(x) - 1,
                 xright  = max(x) + 1,
                 ybottom = target,
                 ytop    = 100,
                 col     = grDevices::rgb(0.18, 0.55, 0.34, 0.08),
                 border  = NA)

  # Mark first adequate n_post
  adequate_idx <- which(sweep_result$power_pct >= target)
  if (length(adequate_idx)) {
    n_min <- x[min(adequate_idx)]
    graphics::abline(v = n_min, lty = 3, col = "darkgreen", lwd = 1.2)
    graphics::text(x = n_min + 0.5, y = 10,
                   labels = sprintf("n_post = %d", n_min),
                   col = "darkgreen", adj = 0, cex = 0.85)
  }

  graphics::grid(NULL, NULL, lwd = 0.8, col = grDevices::rgb(0, 0, 0, 0.1))
  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' Plot a power heatmap across two parameters
#'
#' Displays estimated power as a colour-coded grid, with one parameter on each
#' axis. Useful for visualising how power responds to combinations of effect
#' size, follow-up duration, noise, or autocorrelation.
#'
#' @param grid_result A data frame as returned by [run_power_grid()], with a
#'   \code{power} column and at least the columns specified by \code{x} and
#'   \code{y}.
#' @param x Character. Name of the column to use as the x-axis variable.
#' @param y Character. Name of the column to use as the y-axis variable.
#' @param xlab Character. x-axis label. Defaults to \code{x}.
#' @param ylab Character. y-axis label. Defaults to \code{y}.
#' @param main Character. Plot title. Default \code{"ITS power heatmap"}.
#' @param palette Character vector of colours for the gradient from 0 to 100\%
#'   power. Defaults to a white-steelblue ramp.
#'
#' @return Invisibly \code{NULL}. Called for its side-effect (plot).
#'
#' @examples
#' \dontrun{
#' grid <- build_param_grid(
#'   n_post       = c(12, 18, 24, 30, 36),
#'   level_change = c(-1, -2, -3),
#'   sigma        = 2.5,
#'   rho          = 0.4
#' )
#' results <- run_power_grid(grid, n_sim = 300)
#' plot_power_heatmap(results, x = "n_post", y = "level_change")
#' }
#'
#' @seealso [run_power_grid()], [build_param_grid()]
#' @export
plot_power_heatmap <- function(grid_result,
                               x       = "n_post",
                               y       = "level_change",
                               xlab    = x,
                               ylab    = y,
                               main    = "ITS power heatmap",
                               palette = NULL) {

  if (!all(c(x, y, "power") %in% names(grid_result)))
    stop(sprintf("'grid_result' must have columns '%s', '%s', and 'power'.", x, y))

  if (is.null(palette)) {
    palette <- grDevices::colorRampPalette(
      c("white", "#cce5f5", "#5ba4d1", "#1a5f8a", "#0a2e4f")
    )(101)
  }

  x_vals <- sort(unique(grid_result[[x]]))
  y_vals <- sort(unique(grid_result[[y]]))

  power_mat <- matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals),
                      dimnames = list(y_vals, x_vals))

  for (i in seq_len(nrow(grid_result))) {
    xi <- as.character(grid_result[[x]][i])
    yi <- as.character(grid_result[[y]][i])
    power_mat[yi, xi] <- grid_result$power[i] * 100
  }

  # Base image
  graphics::image(
    x    = seq_along(x_vals),
    y    = seq_along(y_vals),
    z    = t(power_mat),
    col  = palette,
    zlim = c(0, 100),
    xaxt = "n", yaxt = "n",
    xlab = xlab, ylab = ylab, main = main
  )

  graphics::axis(1, at = seq_along(x_vals), labels = x_vals)
  graphics::axis(2, at = seq_along(y_vals), labels = y_vals, las = 1)

  # Cell labels
  for (xi in seq_along(x_vals)) {
    for (yi in seq_along(y_vals)) {
      pwr <- power_mat[yi, xi]
      if (!is.na(pwr)) {
        text_col <- if (pwr > 60) "white" else "black"
        graphics::text(xi, yi, sprintf("%.0f%%", pwr),
                       col = text_col, cex = 0.8, font = 2)
      }
    }
  }

  # 80% adequacy border highlight
  adequate_xi <- which(apply(power_mat, 2, function(col) any(!is.na(col) & col >= 80)))
  if (length(adequate_xi)) {
    graphics::abline(v = min(adequate_xi) - 0.5, lty = 2, col = "firebrick", lwd = 1.5)
  }

  # Legend strip along x-axis label area
  graphics::mtext("  < 60% underpowered  |  60-79% borderline  |  >= 80% adequate",
                  side = 1, line = 3, cex = 0.72, col = "grey40")

  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' Plot a simulated ITS example
#'
#' Generates and plots one realisation of an ITS dataset, showing the
#' pre-intervention trend, the post-intervention trajectory, the intervention
#' breakpoint, and the fitted segmented regression lines. Useful for
#' illustrating the ITS model in papers and presentations.
#'
#' @param n_pre Integer. Pre-intervention time points.
#' @param n_post Integer. Post-intervention time points.
#' @param baseline Numeric. Baseline outcome.
#' @param level_change Numeric. Level change at the intervention.
#' @param slope_change Numeric. Slope change after the intervention. Default 0.
#' @param sigma Numeric. Residual SD.
#' @param rho Numeric. AR(1) autocorrelation.
#' @param seed Integer. Random seed for reproducibility. Default 42.
#' @param xlab Character. x-axis label.
#' @param ylab Character. y-axis label.
#' @param main Character. Plot title.
#' @param pre_col Character. Colour for pre-intervention points.
#' @param post_col Character. Colour for post-intervention points.
#'
#' @return Invisibly, the simulated data frame.
#'
#' @examples
#' plot_its_example(
#'   n_pre = 24, n_post = 24,
#'   baseline = 15, level_change = -3,
#'   sigma = 2.5, rho = 0.4
#' )
#'
#' @export
plot_its_example <- function(n_pre        = 24L,
                             n_post       = 24L,
                             baseline     = 15,
                             level_change = -3,
                             slope_change = 0,
                             sigma        = 2.5,
                             rho          = 0.4,
                             seed         = 42L,
                             xlab         = "Time",
                             ylab         = "Outcome",
                             main         = "Simulated ITS example",
                             pre_col      = "steelblue",
                             post_col     = "firebrick") {

  dat <- simulate_its_data(
    n_pre        = n_pre,
    n_post       = n_post,
    baseline     = baseline,
    level_change = level_change,
    slope_change = slope_change,
    sigma        = sigma,
    rho          = rho
  )

  # Fitted model for overlay
  fit <- tryCatch(
    nlme::gls(y ~ time + D + time_after, data = dat,
              correlation = nlme::corAR1(form = ~ time), method = "ML"),
    error = function(e) NULL
  )

  y_range <- range(dat$y, na.rm = TRUE)
  y_pad   <- diff(y_range) * 0.15

  graphics::plot(dat$time, dat$y,
       type = "n",
       xlim = c(1, n_pre + n_post),
       ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
       xlab = xlab, ylab = ylab, main = main,
       las = 1, bty = "l")

  # Shaded regions
  graphics::rect(xleft   = 0.5,
                 xright  = n_pre + 0.5,
                 ybottom = -1e9, ytop = 1e9,
                 col     = grDevices::rgb(0.18, 0.45, 0.70, 0.05),
                 border  = NA)
  graphics::rect(xleft   = n_pre + 0.5,
                 xright  = n_pre + n_post + 0.5,
                 ybottom = -1e9, ytop = 1e9,
                 col     = grDevices::rgb(0.70, 0.18, 0.18, 0.05),
                 border  = NA)

  # Intervention line
  graphics::abline(v = n_pre + 0.5, lty = 2, col = "grey40", lwd = 1.5)

  # Data points and lines
  pre_idx  <- dat$D == 0
  post_idx <- dat$D == 1
  graphics::lines(dat$time[pre_idx],  dat$y[pre_idx],  col = pre_col,  lwd = 1.2)
  graphics::lines(dat$time[post_idx], dat$y[post_idx], col = post_col, lwd = 1.2)
  graphics::points(dat$time[pre_idx],  dat$y[pre_idx],  pch = 16, col = pre_col,  cex = 0.9)
  graphics::points(dat$time[post_idx], dat$y[post_idx], pch = 16, col = post_col, cex = 0.9)

  # Fitted trend lines
  if (!is.null(fit)) {
    graphics::lines(dat$time[pre_idx],  stats::fitted(fit)[pre_idx],
                    col = "grey20", lwd = 2, lty = 1)
    graphics::lines(dat$time[post_idx], stats::fitted(fit)[post_idx],
                    col = "grey20", lwd = 2, lty = 1)
  }

  # Labels
  graphics::text(x = n_pre / 2, y = y_range[2] + y_pad * 0.8,
                 labels = "Pre-intervention", col = pre_col, cex = 0.85)
  graphics::text(x = n_pre + n_post / 2, y = y_range[2] + y_pad * 0.8,
                 labels = "Post-intervention", col = post_col, cex = 0.85)

  graphics::grid(NULL, NULL, lwd = 0.8, col = grDevices::rgb(0, 0, 0, 0.1))

  invisible(dat)
}


# -----------------------------------------------------------------------------
#' Diagnostic plots for pre-intervention data
#'
#' Produces a 2x2 panel of diagnostic plots for pre-intervention data: the
#' observed series with fitted trend, residuals over time, a Q-Q normality
#' plot, and the residual ACF. These help assess whether the ITS model
#' assumptions are plausible.
#'
#' @param data Pre-intervention data: a data frame or numeric vector. See
#'   [estimate_its_params()] for format details.
#' @param outcome_col Character. Column name for the outcome. Default
#'   \code{"outcome"}.
#' @param time_col Character. Column name for time. Default \code{"time"}.
#' @param main Character. Panel title prefix.
#'
#' @return Invisibly, a list with elements \code{params} (estimated parameters)
#'   and \code{residuals} (residual vector).
#'
#' @examples
#' data("example_cfr_data")
#' diagnose_params(example_cfr_data)
#'
#' @seealso [estimate_its_params()]
#' @export
diagnose_params <- function(data,
                            outcome_col = "outcome",
                            time_col    = "time",
                            main        = "Pre-intervention diagnostics") {

  if (is.numeric(data)) {
    outcome  <- data
    time_idx <- seq_along(outcome)
  } else {
    data <- as.data.frame(data)
    outcome  <- data[[outcome_col]]
    time_idx <- if (time_col %in% names(data)) {
      .coerce_time(data[[time_col]])
    } else {
      seq_along(outcome)
    }
  }

  keep     <- !is.na(outcome)
  outcome  <- outcome[keep]
  time_idx <- time_idx[keep]

  fit   <- stats::lm(outcome ~ time_idx)
  resid <- stats::residuals(fit)

  old_par <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(old_par), add = TRUE)

  # 1. Raw series with fitted trend
  graphics::plot(time_idx, outcome,
       type = "o", pch = 16, col = "steelblue", cex = 0.8,
       xlab = "Time", ylab = outcome_col,
       main = paste(main, "- observed + trend"), las = 1, bty = "l")
  graphics::lines(time_idx, stats::fitted(fit), col = "firebrick", lwd = 2)
  graphics::legend("topright", legend = c("Observed", "Fitted trend"),
                   col = c("steelblue", "firebrick"),
                   lty = 1, pch = c(16, NA), bty = "n", cex = 0.8)
  graphics::grid(NULL, NULL, lwd = 0.6)

  # 2. Residuals over time
  graphics::plot(time_idx, resid,
       type = "o", pch = 16, col = "steelblue", cex = 0.8,
       xlab = "Time", ylab = "Residuals",
       main = paste(main, "- residuals"), las = 1, bty = "l")
  graphics::abline(h = 0, col = "firebrick", lty = 2, lwd = 1.5)
  graphics::grid(NULL, NULL, lwd = 0.6)

  # 3. Q-Q plot
  stats::qqnorm(resid, main = paste(main, "- Q-Q plot"), las = 1,
                pch = 16, col = "steelblue", cex = 0.8)
  stats::qqline(resid, col = "firebrick", lwd = 2)

  # 4. ACF of residuals
  acf_obj <- stats::acf(resid, lag.max = min(12, length(resid) - 1),
                        plot = FALSE)
  graphics::barplot(
    height = acf_obj$acf[-1],
    names.arg = seq_len(length(acf_obj$acf) - 1),
    xlab = "Lag", ylab = "ACF",
    main = paste(main, "- residual ACF"),
    col  = ifelse(acf_obj$acf[-1] > 0, "steelblue", "firebrick"),
    border = NA, las = 1
  )
  ci <- stats::qnorm(0.975) / sqrt(length(resid))
  graphics::abline(h = c(-ci, ci), lty = 2, col = "grey40")

  params <- list(
    n_pre     = length(outcome),
    baseline  = estimate_baseline(outcome, time_idx),
    sigma     = estimate_sigma(outcome, time_idx),
    rho       = estimate_rho(outcome, time_idx),
    trend_pre = estimate_trend(outcome, time_idx)
  )

  invisible(list(params = params, residuals = resid))
}
