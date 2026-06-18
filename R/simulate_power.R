# =============================================================================
#  PITS - Power simulation functions
#
#  Core Monte Carlo engine for estimating the statistical power of ITS designs.
#
#  simulate_its_data()      - data-generating process (single site)
#  fit_its_model()          - segmented regression with AR(1) correction
#  calculate_power()        - Monte Carlo power for a single-site design
#  calculate_power_multi()  - Monte Carlo power for a multi-site design
#  power_sweep()            - power across a range of n_post values
#  build_param_grid()       - factorial parameter grid
#  run_power_grid()         - run calculate_power() across a grid
# =============================================================================


# -----------------------------------------------------------------------------
#' Simulate a single ITS dataset
#'
#' Generates one realisation of an ITS time series under the alternative
#' hypothesis (i.e. with a true intervention effect), with AR(1)
#' autocorrelated errors. Used internally by [calculate_power()].
#'
#' @param n_pre Integer. Number of pre-intervention time points.
#' @param n_post Integer. Number of post-intervention time points.
#' @param baseline Numeric. Mean outcome at \eqn{t = 1} (before any trend).
#' @param level_change Numeric. Immediate step change in outcome at the
#'   intervention point. Set to 0 if testing slope only.
#' @param slope_change Numeric. Change in trend per time unit after the
#'   intervention. Set to 0 if testing level only.
#' @param sigma Numeric. Residual standard deviation (noise).
#' @param rho Numeric. AR(1) autocorrelation coefficient. Must be in
#'   \eqn{[0, 1)}.
#' @param pre_trend Numeric. Pre-intervention trend (slope) per time unit.
#'   Default 0 (stable pre-period).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{time}{Integer time index, 1 to \eqn{n_{pre} + n_{post}}.}
#'     \item{D}{Binary intervention indicator (0 = pre, 1 = post).}
#'     \item{time_after}{Time elapsed since intervention (0 in pre-period).}
#'     \item{y}{Simulated outcome values.}
#'   }
#'
#' @details
#' The data-generating process is the standard segmented-regression ITS model:
#' \deqn{y_t = \beta_0 + \beta_1 t + \delta D_t + \gamma T_t^* + \varepsilon_t}
#' where \eqn{D_t = \mathbf{1}(t > n_{pre})},
#' \eqn{T_t^* = (t - n_{pre}) \cdot D_t},
#' \eqn{\delta} = \code{level_change}, \eqn{\gamma} = \code{slope_change},
#' and \eqn{\varepsilon_t} follows an AR(1) process with innovation SD
#' \eqn{\sigma \sqrt{1 - \rho^2}}.
#'
#' @examples
#' df <- simulate_its_data(
#'   n_pre = 24, n_post = 24,
#'   baseline = 15, level_change = -3,
#'   slope_change = 0, sigma = 2.5, rho = 0.4
#' )
#' plot(df$time, df$y, type = "o", pch = 16,
#'      xlab = "Time", ylab = "Outcome")
#' abline(v = 24.5, lty = 2, col = "red")
#'
#' @seealso [calculate_power()], [fit_its_model()]
#' @export
simulate_its_data <- function(n_pre,
                              n_post,
                              baseline,
                              level_change,
                              slope_change,
                              sigma,
                              rho,
                              pre_trend = 0) {

  .validate_sim_params(n_pre, n_post, sigma, rho)

  n_total    <- n_pre + n_post
  time       <- seq_len(n_total)
  D          <- as.integer(time > n_pre)
  time_after <- (time - n_pre) * D

  mu <- baseline +
        pre_trend    * time +
        level_change * D +
        slope_change * time_after

  innovation_sd <- sigma * sqrt(max(0, 1 - rho^2))
  eps <- numeric(n_total)
  for (i in seq_len(n_total)) {
    prev   <- if (i > 1L) eps[i - 1L] else 0
    eps[i] <- rho * prev + stats::rnorm(1L, 0, innovation_sd)
  }

  data.frame(
    time       = time,
    D          = D,
    time_after = time_after,
    y          = mu + eps
  )
}


# -----------------------------------------------------------------------------
#' Fit a segmented regression model to an ITS dataset
#'
#' Fits a Gaussian segmented-regression model with AR(1) autocorrelation
#' correction using [nlme::gls()], and returns the p-value for the
#' coefficient(s) specified by \code{test}.
#'
#' @param data A data frame as returned by [simulate_its_data()], containing
#'   columns \code{time}, \code{D}, \code{time_after}, and \code{y}.
#' @param test Character. Which effect to test:
#'   \describe{
#'     \item{\code{"level"}}{P-value for the immediate step change (\eqn{\delta}).
#'       Use for acute interventions.}
#'     \item{\code{"slope"}}{P-value for the slope change (\eqn{\gamma}).
#'       Use for gradual interventions.}
#'     \item{\code{"both"}}{Likelihood-ratio test p-value for the joint
#'       null \eqn{\delta = \gamma = 0}. Uses 2 degrees of freedom.}
#'   }
#'
#' @return A single numeric p-value, or \code{NA} if the model failed to
#'   converge.
#'
#' @details
#' The model fitted is:
#' \deqn{y_t = \beta_0 + \beta_1 t + \delta D_t + \gamma T_t^* + \varepsilon_t,
#'   \quad \varepsilon_t \sim AR(1)}
#' Estimation is by maximum likelihood (\code{method = "ML"}) to support
#' likelihood-ratio tests. For the \code{"both"} option, the full model is
#' compared against a model with only a linear trend (\eqn{\delta = \gamma = 0}).
#'
#' @examples
#' df <- simulate_its_data(
#'   n_pre = 24, n_post = 24,
#'   baseline = 15, level_change = -3,
#'   slope_change = 0, sigma = 2.5, rho = 0.4
#' )
#' fit_its_model(df, test = "level")
#'
#' @seealso [simulate_its_data()], [calculate_power()]
#' @export
fit_its_model <- function(data, test = c("level", "slope", "both")) {
  test <- match.arg(test)

  tryCatch({
    fit <- nlme::gls(
      y ~ time + D + time_after,
      data        = data,
      correlation = nlme::corAR1(form = ~ time),
      method      = "ML"
    )

    if (test == "level") {
      summary(fit)$tTable["D", "p-value"]

    } else if (test == "slope") {
      summary(fit)$tTable["time_after", "p-value"]

    } else {   # "both" - likelihood-ratio test
      fit_null <- nlme::gls(
        y ~ time,
        data        = data,
        correlation = nlme::corAR1(form = ~ time),
        method      = "ML"
      )
      stats::anova(fit_null, fit)$`p-value`[2]
    }
  }, error = function(e) NA_real_)
}


# -----------------------------------------------------------------------------
#' Estimate statistical power for a single-site ITS design
#'
#' Runs a Monte Carlo simulation to estimate the probability of detecting a
#' specified intervention effect in a single-site ITS study.
#'
#' @param n_pre Integer. Number of pre-intervention time points.
#' @param n_post Integer. Number of post-intervention time points. This is
#'   the primary design lever: use [power_sweep()] to find the minimum
#'   \code{n_post} that achieves \eqn{\ge 80\%} power.
#' @param baseline Numeric. Mean outcome in the pre-intervention period.
#' @param level_change Numeric. Expected immediate step change at the
#'   intervention point (your minimum clinically meaningful effect). Set to
#'   0 when \code{test = "slope"}.
#' @param slope_change Numeric. Expected change in trend per time unit after
#'   the intervention. Default 0. Set to 0 when \code{test = "level"}.
#' @param sigma Numeric. Residual standard deviation. Estimate from
#'   pre-intervention data using [estimate_sigma()] or [estimate_its_params()].
#' @param rho Numeric. AR(1) autocorrelation coefficient in \eqn{[0, 1)}.
#'   Use 0.4 as a conservative default for monthly data if unknown.
#' @param test Character. Effect to test: \code{"level"} (default),
#'   \code{"slope"}, or \code{"both"}. See [fit_its_model()] for details.
#' @param alpha Numeric. Significance threshold. Default 0.05.
#' @param n_sim Integer. Number of Monte Carlo replications. Use 500 for a
#'   quick check, 1000 for a reportable estimate, 2000+ for publication.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#' @param pre_trend Numeric. Pre-intervention trend per time unit. Default 0.
#'
#' @return A named list:
#'   \describe{
#'     \item{power}{Numeric. Estimated power (proportion of simulations
#'       with \eqn{p < \alpha}).}
#'     \item{power_pct}{Numeric. Power as a percentage.}
#'     \item{interpretation}{Character. Qualitative label (see
#'       [interpret_power()]).}
#'     \item{p_values}{Numeric vector. Raw p-values from all \code{n_sim}
#'       replications (\code{NA} = non-convergence).}
#'     \item{n_converged}{Integer. Number of replications that converged.}
#'     \item{n_sim}{Integer. Total replications requested.}
#'     \item{params}{Named list. Input parameters, for traceability.}
#'   }
#'
#' @examples
#' result <- calculate_power(
#'   n_pre = 24, n_post = 24,
#'   baseline = 15, level_change = -3,
#'   sigma = 2.5, rho = 0.4,
#'   n_sim = 500, seed = 42
#' )
#' print(result$power_pct)
#'
#' @seealso [power_sweep()], [calculate_power_multi()], [estimate_its_params()]
#' @export
calculate_power <- function(n_pre,
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
                            pre_trend    = 0) {

  test <- match.arg(test)
  .validate_power_params(n_pre, n_post, sigma, rho, alpha, n_sim)
  if (!is.null(seed)) set.seed(seed)

  p_values <- replicate(n_sim, {
    dat <- simulate_its_data(
      n_pre        = n_pre,
      n_post       = n_post,
      baseline     = baseline,
      level_change = level_change,
      slope_change = slope_change,
      sigma        = sigma,
      rho          = rho,
      pre_trend    = pre_trend
    )
    fit_its_model(dat, test = test)
  })

  power       <- mean(p_values < alpha, na.rm = TRUE)
  n_converged <- sum(!is.na(p_values))

  result <- list(
    power          = power,
    power_pct      = round(power * 100, 1),
    interpretation = interpret_power(power),
    p_values       = p_values,
    n_converged    = n_converged,
    n_sim          = n_sim,
    params         = list(
      n_pre        = n_pre,
      n_post       = n_post,
      baseline     = baseline,
      level_change = level_change,
      slope_change = slope_change,
      sigma        = sigma,
      rho          = rho,
      test         = test,
      alpha        = alpha,
      seed         = seed
    )
  )
  class(result) <- c("pits_power_result", "list")
  result
}


# -----------------------------------------------------------------------------
#' Estimate statistical power for a multi-site ITS design
#'
#' Runs a Monte Carlo simulation to estimate the probability of detecting a
#' common intervention effect across multiple sites (e.g. hospitals), using
#' a mixed-effects segmented regression model.
#'
#' @param sites A list of named parameter lists, one per site. Each element
#'   must contain: \code{name} (character), \code{n_pre}, \code{n_post},
#'   \code{baseline}, \code{level_change}, \code{slope_change},
#'   \code{sigma}, \code{rho}. All sites must have the same \code{n_pre}
#'   and \code{n_post}.
#' @inheritParams calculate_power
#'
#' @return Same structure as [calculate_power()], plus:
#'   \describe{
#'     \item{n_sites}{Integer. Number of sites.}
#'     \item{site_names}{Character vector. Site names.}
#'   }
#'
#' @details
#' The model is a linear mixed-effects model with site-specific random
#' intercepts and AR(1) autocorrelation within each site:
#' \deqn{y_{it} = \beta_0 + \beta_1 t + \delta D_t + \gamma T_t^* + u_i + \varepsilon_{it}}
#' where \eqn{u_i \sim N(0, \tau^2)} are site random intercepts and
#' \eqn{\varepsilon_{it}} follows AR(1) within each site.
#'
#' @examples
#' sites <- list(
#'   list(name = "Hospital A", n_pre = 24, n_post = 24,
#'        baseline = 15, level_change = -3, slope_change = 0,
#'        sigma = 2.5, rho = 0.4),
#'   list(name = "Hospital B", n_pre = 24, n_post = 24,
#'        baseline = 18, level_change = -3, slope_change = 0,
#'        sigma = 3.0, rho = 0.4)
#' )
#' result <- calculate_power_multi(sites, n_sim = 200, seed = 42)
#' result$power_pct
#'
#' @seealso [calculate_power()], [power_sweep()]
#' @export
calculate_power_multi <- function(sites,
                                  test  = c("level", "slope", "both"),
                                  alpha = 0.05,
                                  n_sim = 1000L,
                                  seed  = 123L) {

  test <- match.arg(test)
  if (!is.list(sites) || length(sites) < 2)
    stop("'sites' must be a list of at least 2 site parameter lists.")
  .validate_power_params(
    sites[[1]]$n_pre, sites[[1]]$n_post,
    sites[[1]]$sigma, sites[[1]]$rho,
    alpha, n_sim
  )
  if (!is.null(seed)) set.seed(seed)

  p_values <- replicate(n_sim, {
    multi_dat <- .simulate_multi(sites)
    tryCatch({
      fit <- nlme::lme(
        y ~ time + D + time_after,
        data        = multi_dat,
        random      = ~ 1 | site,
        correlation = nlme::corAR1(form = ~ time | site),
        method      = "ML"
      )
      coefs <- summary(fit)$tTable
      if (test == "level") {
        coefs["D", "p-value"]
      } else if (test == "slope") {
        coefs["time_after", "p-value"]
      } else {
        fit_null <- nlme::lme(
          y ~ time,
          data        = multi_dat,
          random      = ~ 1 | site,
          correlation = nlme::corAR1(form = ~ time | site),
          method      = "ML"
        )
        stats::anova(fit_null, fit)$`p-value`[2]
      }
    }, error = function(e) NA_real_)
  })

  power       <- mean(p_values < alpha, na.rm = TRUE)
  n_converged <- sum(!is.na(p_values))

  result <- list(
    power          = power,
    power_pct      = round(power * 100, 1),
    interpretation = interpret_power(power),
    p_values       = p_values,
    n_converged    = n_converged,
    n_sim          = n_sim,
    n_sites        = length(sites),
    site_names     = sapply(sites, `[[`, "name"),
    params         = list(sites = sites, test = test, alpha = alpha, seed = seed)
  )
  class(result) <- c("pits_power_result", "list")
  result
}


# -----------------------------------------------------------------------------
#' Design optimisation sweep: power across a range of n_post values
#'
#' Runs [calculate_power()] for a vector of post-intervention durations,
#' returning a data frame of power estimates. Use this to identify the minimum
#' \code{n_post} needed to achieve \eqn{\ge 80\%} power.
#'
#' @param sweep_post Integer vector. \code{n_post} values to evaluate.
#'   Default \code{c(6, 12, 18, 24, 30, 36)}.
#' @param ... All other arguments passed to [calculate_power()]. Note that
#'   \code{n_post} is overridden by \code{sweep_post} internally - do not
#'   pass it separately.
#' @param verbose Logical. If \code{TRUE} (default), prints a formatted
#'   sweep table to the console.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{n_post}{Integer. Follow-up duration evaluated.}
#'     \item{power}{Numeric. Estimated power (0-1).}
#'     \item{power_pct}{Numeric. Power as a percentage.}
#'     \item{interpretation}{Character. Qualitative label.}
#'     \item{n_converged}{Integer. Number of converged replications.}
#'   }
#'
#' @examples
#' # Small n_sim for a fast example; use 1000+ for a reportable estimate.
#' sweep <- power_sweep(
#'   sweep_post   = c(12, 24, 36),
#'   n_pre        = 24,
#'   baseline     = 15,
#'   level_change = -3,
#'   sigma        = 2.5,
#'   rho          = 0.4,
#'   n_sim        = 100,
#'   seed         = 42,
#'   verbose      = FALSE
#' )
#' plot_power_curve(sweep)
#'
#' @seealso [calculate_power()], [plot_power_curve()]
#' @export
power_sweep <- function(sweep_post = c(6L, 12L, 18L, 24L, 30L, 36L),
                        ...,
                        verbose = TRUE) {

  args <- list(...)
  if ("n_post" %in% names(args)) {
    warning("'n_post' in ... is ignored; use 'sweep_post' to specify values.")
    args$n_post <- NULL
  }

  results <- lapply(seq_along(sweep_post), function(i) {
    np      <- sweep_post[i]
    seed_i  <- if (!is.null(args$seed)) args$seed + np else NULL
    call_args <- args                   # args already has seed if supplied
    call_args$n_post <- np
    call_args$seed   <- seed_i          # override with per-sweep seed
    do.call(calculate_power, call_args)
  })

  out <- data.frame(
    n_post         = sweep_post,
    power          = sapply(results, `[[`, "power"),
    power_pct      = sapply(results, `[[`, "power_pct"),
    interpretation = sapply(results, `[[`, "interpretation"),
    n_converged    = sapply(results, `[[`, "n_converged"),
    stringsAsFactors = FALSE
  )
  class(out) <- c("pits_sweep_result", "data.frame")

  if (verbose) print(out)

  invisible(out)
}


# -----------------------------------------------------------------------------
#' Build a factorial parameter grid for power calculations
#'
#' Creates a data frame of all combinations of the supplied parameter vectors,
#' suitable for passing to [run_power_grid()]. Useful for sensitivity analyses
#' and generating figures showing how power varies across parameter space.
#'
#' @param n_post Integer vector. Post-intervention durations to evaluate.
#' @param level_change Numeric vector. Effect sizes to evaluate.
#' @param sigma Numeric vector. Residual SDs to evaluate.
#' @param rho Numeric vector. Autocorrelation values to evaluate.
#' @param n_pre Integer. Pre-intervention duration (fixed). Default 24.
#' @param baseline Numeric. Baseline outcome (fixed). Default 15.
#' @param slope_change Numeric. Slope change (fixed). Default 0.
#'
#' @return A data frame with one row per parameter combination.
#'
#' @examples
#' grid <- build_param_grid(
#'   n_post       = c(12, 18, 24, 30),
#'   level_change = c(-1, -2, -3),
#'   sigma        = c(1.5, 2.5, 3.5),
#'   rho          = c(0.2, 0.4, 0.6)
#' )
#' nrow(grid)   # 4 * 3 * 3 * 3 = 108 combinations
#'
#' @seealso [run_power_grid()]
#' @export
build_param_grid <- function(n_post,
                             level_change,
                             sigma,
                             rho,
                             n_pre        = 24L,
                             baseline     = 15,
                             slope_change = 0) {
  expand.grid(
    n_pre        = n_pre,
    n_post       = n_post,
    baseline     = baseline,
    level_change = level_change,
    slope_change = slope_change,
    sigma        = sigma,
    rho          = rho,
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
#' Run power calculations across a parameter grid
#'
#' Applies [calculate_power()] to each row of a parameter grid as produced by
#' [build_param_grid()], returning a data frame with the estimated power for
#' each combination.
#'
#' @param grid A data frame as returned by [build_param_grid()], or any data
#'   frame with columns \code{n_pre}, \code{n_post}, \code{baseline},
#'   \code{level_change}, \code{slope_change}, \code{sigma}, \code{rho}.
#' @param test Character. Effect to test: \code{"level"}, \code{"slope"}, or
#'   \code{"both"}.
#' @param alpha Numeric. Significance threshold. Default 0.05.
#' @param n_sim Integer. Monte Carlo replications per cell. Default 500.
#' @param seed Integer. Base random seed; each row adds its row index.
#' @param verbose Logical. If \code{TRUE}, prints a progress counter.
#'
#' @return The input \code{grid} data frame with columns appended:
#'   \code{power}, \code{power_pct}, \code{interpretation}, \code{n_converged}.
#'
#' @examples
#' \donttest{
#' grid <- build_param_grid(
#'   n_post       = c(12, 24, 36),
#'   level_change = c(-2, -3),
#'   sigma        = 2.5,
#'   rho          = 0.4
#' )
#' results <- run_power_grid(grid, n_sim = 100, verbose = FALSE)
#' plot_power_heatmap(results, x = "n_post", y = "level_change")
#' }
#'
#' @seealso [build_param_grid()], [plot_power_heatmap()]
#' @export
run_power_grid <- function(grid,
                           test    = c("level", "slope", "both"),
                           alpha   = 0.05,
                           n_sim   = 500L,
                           seed    = 123L,
                           verbose = TRUE) {

  test <- match.arg(test)
  n_rows <- nrow(grid)

  if (verbose) {
    cat(sprintf("Running power grid: %d combinations, %d sims each...\n",
                n_rows, n_sim))
  }

  results <- lapply(seq_len(n_rows), function(i) {
    if (verbose && i %% 10 == 0)
      cat(sprintf("  Row %d / %d\n", i, n_rows))
    r <- grid[i, ]
    calculate_power(
      n_pre        = r$n_pre,
      n_post       = r$n_post,
      baseline     = r$baseline,
      level_change = r$level_change,
      slope_change = r$slope_change,
      sigma        = r$sigma,
      rho          = r$rho,
      test         = test,
      alpha        = alpha,
      n_sim        = n_sim,
      seed         = seed + i
    )
  })

  grid$power          <- sapply(results, `[[`, "power")
  grid$power_pct      <- sapply(results, `[[`, "power_pct")
  grid$interpretation <- sapply(results, `[[`, "interpretation")
  grid$n_converged    <- sapply(results, `[[`, "n_converged")

  if (verbose) cat("Done.\n")
  grid
}


# =============================================================================
#  Internal helpers
# =============================================================================

#' @keywords internal
.simulate_multi <- function(sites) {
  do.call(rbind, lapply(seq_along(sites), function(s_idx) {
    sp  <- sites[[s_idx]]
    dat <- simulate_its_data(
      n_pre        = sp$n_pre,
      n_post       = sp$n_post,
      baseline     = sp$baseline,
      level_change = sp$level_change,
      slope_change = sp$slope_change,
      sigma        = sp$sigma,
      rho          = sp$rho
    )
    dat$site      <- s_idx
    dat$site_name <- sp$name
    dat
  }))
}

#' @keywords internal
.validate_sim_params <- function(n_pre, n_post, sigma, rho) {
  if (!is.numeric(n_pre)  || n_pre  < 2) stop("'n_pre' must be an integer >= 2.")
  if (!is.numeric(n_post) || n_post < 1) stop("'n_post' must be an integer >= 1.")
  if (!is.numeric(sigma)  || sigma  <= 0) stop("'sigma' must be positive.")
  if (!is.numeric(rho) || abs(rho) >= 1)
    stop("'rho' must be in (-1, 1).")
}

#' @keywords internal
.validate_power_params <- function(n_pre, n_post, sigma, rho, alpha, n_sim) {
  .validate_sim_params(n_pre, n_post, sigma, rho)
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be in (0, 1).")
  if (!is.numeric(n_sim) || n_sim < 10)
    stop("'n_sim' must be an integer >= 10.")
}
