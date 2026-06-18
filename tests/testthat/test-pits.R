## =============================================================================
##  PITS test suite  (testthat 3)
## =============================================================================

library(testthat)
library(PITS)


# =============================================================================
#  estimate_params.R
# =============================================================================

test_that("estimate_baseline returns a numeric scalar", {
  b <- estimate_baseline(c(14, 15, 13, 16, 14, 15))
  expect_length(b, 1)
  expect_true(is.numeric(b))
})

test_that("estimate_baseline is close to mean when trend is flat", {
  set.seed(1)
  x <- rnorm(24, mean = 15, sd = 0.5)
  expect_equal(estimate_baseline(x), mean(x), tolerance = 1)
})

test_that("estimate_sigma returns a positive scalar", {
  s <- estimate_sigma(c(14, 15, 13, 16, 14, 15, 14, 16))
  expect_length(s, 1)
  expect_gt(s, 0)
})

test_that("estimate_rho returns a value in (-1, 1)", {
  set.seed(2)
  x <- as.numeric(arima.sim(list(ar = 0.5), n = 24)) + 15
  r <- estimate_rho(x)
  expect_length(r, 1)
  expect_gt(r, -1)
  expect_lt(r, 1)
})

test_that("estimate_rho warns and returns NA for fewer than 3 obs", {
  expect_warning(
    r <- estimate_rho(c(1, 2)),
    regexp = "fewer than 3"
  )
  expect_true(is.na(r))
})

test_that("estimate_trend returns a numeric scalar", {
  tr <- estimate_trend(seq_len(24) + rnorm(24))
  expect_length(tr, 1)
  expect_true(is.numeric(tr))
})

test_that("estimate_its_params returns a correctly named list", {
  data("example_cfr_data", package = "PITS")
  params <- estimate_its_params(example_cfr_data, verbose = FALSE)
  expect_type(params, "list")
  expect_named(params,
               c("n_pre", "baseline", "sigma", "rho", "trend_pre", "method"))
  expect_equal(params$n_pre, nrow(example_cfr_data))
  expect_gt(params$sigma, 0)
  expect_gt(params$rho, -1)
  expect_lt(params$rho,  1)
  expect_true(params$method %in% c("GLS-ML", "OLS"))
})

test_that("estimate_its_params accepts a plain numeric vector", {
  set.seed(5)
  params <- estimate_its_params(rnorm(24, 15, 2), verbose = FALSE)
  expect_equal(params$n_pre, 24)
})

test_that("estimate_its_params warns when fewer than 12 obs", {
  expect_warning(
    estimate_its_params(rnorm(8, 15, 2), verbose = FALSE),
    regexp = "unreliable"
  )
})

test_that("estimate_its_params removes NAs and warns", {
  data("example_cfr_data", package = "PITS")
  d <- example_cfr_data
  d$outcome[c(3, 7)] <- NA
  expect_warning(
    p <- estimate_its_params(d, verbose = FALSE),
    regexp = "missing"
  )
  expect_equal(p$n_pre, nrow(d) - 2L)
})


# =============================================================================
#  simulate_power.R
# =============================================================================

test_that("simulate_its_data returns correct structure", {
  df <- simulate_its_data(12, 12, 15, -3, 0, 2, 0.4)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("time", "D", "time_after", "y"))
  expect_equal(nrow(df), 24L)
  expect_equal(sum(df$D == 0L), 12L)
  expect_equal(sum(df$D == 1L), 12L)
  expect_equal(sum(df$time_after[df$D == 0]), 0)
})

test_that("simulate_its_data pre-period mean is near baseline (Monte Carlo)", {
  means <- replicate(300, {
    df <- simulate_its_data(24, 24, 15, 0, 0, 0.5, 0)
    mean(df$y[df$D == 0])
  })
  expect_equal(mean(means), 15, tolerance = 0.4)
})

test_that("simulate_its_data accepts negative rho", {
  expect_no_error(simulate_its_data(12, 12, 15, -3, 0, 2, -0.5))
})

test_that("fit_its_model returns a numeric p-value or NA", {
  df <- simulate_its_data(24, 24, 15, -3, 0, 2, 0.4)
  p  <- fit_its_model(df, test = "level")
  expect_true(is.numeric(p))
  expect_true(is.na(p) || (p >= 0 && p <= 1))
})

test_that("fit_its_model accepts all three test options", {
  df <- simulate_its_data(24, 24, 15, -3, -0.1, 2, 0.4)
  for (t in c("level", "slope", "both")) {
    p <- fit_its_model(df, test = t)
    expect_true(is.na(p) || (p >= 0 && p <= 1),
                label = paste("test =", t))
  }
})

test_that("calculate_power returns a pits_power_result with correct structure", {
  r <- calculate_power(12, 12, 15, -3, 0, 2, 0.4, n_sim = 50, seed = 1)
  expect_s3_class(r, "pits_power_result")
  expect_named(r, c("power", "power_pct", "interpretation",
                    "p_values", "n_converged", "n_sim", "params"))
  expect_gte(r$power, 0)
  expect_lte(r$power, 1)
  expect_equal(r$power_pct, round(r$power * 100, 1))
  expect_length(r$p_values, 50L)
  expect_equal(r$n_sim, 50L)
})

test_that("calculate_power is reproducible with the same seed", {
  args <- list(n_pre = 12, n_post = 12, baseline = 15,
               level_change = -3, sigma = 2, rho = 0.4,
               n_sim = 30, seed = 77)
  r1 <- do.call(calculate_power, args)
  r2 <- do.call(calculate_power, args)
  expect_equal(r1$power, r2$power)
})

test_that("power increases monotonically with effect size (stochastic check)", {
  r_small <- calculate_power(12, 12, 15, -1, sigma = 2, rho = 0.4, n_sim = 300, seed = 5)
  r_large <- calculate_power(12, 12, 15, -5, sigma = 2, rho = 0.4, n_sim = 300, seed = 5)
  expect_gte(r_large$power, r_small$power)
})

test_that("power increases with longer follow-up (stochastic check)", {
  r_short <- calculate_power(24, 12, 15, -2, sigma = 2, rho = 0.4, n_sim = 300, seed = 10)
  r_long  <- calculate_power(24, 36, 15, -2, sigma = 2, rho = 0.4, n_sim = 300, seed = 10)
  expect_gte(r_long$power, r_short$power)
})

test_that("power_sweep returns a pits_sweep_result with correct columns", {
  sw <- power_sweep(
    sweep_post = c(12, 24), n_pre = 12, baseline = 15,
    level_change = -3, sigma = 2, rho = 0.4,
    n_sim = 30, seed = 1, verbose = FALSE
  )
  expect_s3_class(sw, c("pits_sweep_result", "data.frame"))
  expect_named(sw, c("n_post", "power", "power_pct",
                     "interpretation", "n_converged"))
  expect_equal(nrow(sw), 2L)
})

test_that("power_sweep power values increase with n_post", {
  sw <- power_sweep(
    sweep_post = c(12, 24, 36), n_pre = 24, baseline = 15,
    level_change = -3, sigma = 2, rho = 0.4,
    n_sim = 200, seed = 1, verbose = FALSE
  )
  expect_gte(sw$power[2], sw$power[1])
  expect_gte(sw$power[3], sw$power[2])
})

test_that("build_param_grid returns the correct row count", {
  g <- build_param_grid(
    n_post       = c(12, 24, 36),
    level_change = c(-1, -2, -3),
    sigma        = c(1, 2),
    rho          = c(0.2, 0.4)
  )
  expect_equal(nrow(g), 3L * 3L * 2L * 2L)
  expect_true(all(c("n_pre", "n_post", "baseline", "level_change",
                    "slope_change", "sigma", "rho") %in% names(g)))
})

test_that("run_power_grid appends power columns to the grid", {
  g <- build_param_grid(n_post = c(12, 24), level_change = -3,
                        sigma = 2, rho = 0.4)
  gr <- run_power_grid(g, n_sim = 50, seed = 1, verbose = FALSE)
  expect_equal(nrow(gr), 2L)
  expect_true(all(c("power", "power_pct", "interpretation", "n_converged")
                  %in% names(gr)))
  expect_true(all(gr$power >= 0 & gr$power <= 1))
})

test_that("calculate_power_multi returns pits_power_result with n_sites", {
  sites <- list(
    list(name = "A", n_pre = 12, n_post = 12,
         baseline = 15, level_change = -3, slope_change = 0,
         sigma = 2, rho = 0.4),
    list(name = "B", n_pre = 12, n_post = 12,
         baseline = 18, level_change = -3, slope_change = 0,
         sigma = 2.5, rho = 0.4)
  )
  r <- calculate_power_multi(sites, n_sim = 30, seed = 1)
  expect_s3_class(r, "pits_power_result")
  expect_equal(r$n_sites, 2L)
  expect_gte(r$power, 0)
  expect_lte(r$power, 1)
  expect_equal(length(r$site_names), 2L)
})


# =============================================================================
#  S3 methods  (methods.R)
# =============================================================================

test_that("print.pits_power_result produces output without error", {
  r <- calculate_power(12, 12, 15, -3, 0, 2, 0.4, n_sim = 30, seed = 1)
  expect_output(print(r), "POWER")
  expect_output(print(r), "%")
})

test_that("summary.pits_power_result produces extended output", {
  r <- calculate_power(12, 12, 15, -3, 0, 2, 0.4, n_sim = 50, seed = 1)
  expect_output(summary(r), "P-value distribution")
})

test_that("print.pits_sweep_result produces tabular output", {
  sw <- power_sweep(c(12, 24), n_pre = 12, baseline = 15,
                   level_change = -3, sigma = 2, rho = 0.4,
                   n_sim = 30, seed = 1, verbose = FALSE)
  expect_output(print(sw), "n_post")
  expect_output(print(sw), "%")
})


# =============================================================================
#  utils.R
# =============================================================================

test_that("interpret_power returns correct labels", {
  expect_equal(interpret_power(0.85), "Adequate (>= 80%)")
  expect_equal(interpret_power(0.80), "Adequate (>= 80%)")
  expect_equal(interpret_power(0.79), "Borderline (60-79%)")
  expect_equal(interpret_power(0.60), "Borderline (60-79%)")
  expect_equal(interpret_power(0.59), "Underpowered (< 60%)")
  expect_equal(interpret_power(0.00), "Underpowered (< 60%)")
})

test_that("interpret_power errors on non-numeric input", {
  expect_error(interpret_power("high"))
  expect_error(interpret_power(c(0.5, 0.8)))
})

test_that("validate_params passes for valid inputs", {
  expect_true(validate_params(24, 24, 15, -3, 2.5, 0.4))
})

test_that("validate_params errors on invalid sigma", {
  expect_error(validate_params(24, 24, 15, -3, -1, 0.4), "sigma")
})

test_that("validate_params errors on invalid rho", {
  expect_error(validate_params(24, 24, 15, -3, 2, 1.5), "rho")
})

test_that("validate_params errors on rho = -1.5", {
  expect_error(validate_params(24, 24, 15, -3, 2, -1.5), "rho")
})

test_that("validate_params warns for small n_sim", {
  expect_warning(validate_params(24, 24, 15, -3, 2, 0.4, n_sim = 50),
                 "imprecise")
})

test_that("simulate_predata returns a data frame with correct dimensions", {
  pre <- simulate_predata(n = 24, baseline = 15, sigma = 2, rho = 0.4)
  expect_s3_class(pre, "data.frame")
  expect_named(pre, c("time", "outcome"))
  expect_equal(nrow(pre), 24L)
  expect_equal(pre$time, seq_len(24L))
})

test_that("simulate_predata is reproducible with the same seed", {
  p1 <- simulate_predata(n = 24, seed = 7)
  p2 <- simulate_predata(n = 24, seed = 7)
  expect_identical(p1$outcome, p2$outcome)
})

test_that("estimate_and_calculate returns a valid power result", {
  data("example_cfr_data", package = "PITS")
  r <- estimate_and_calculate(
    data = example_cfr_data, level_change = -3,
    n_post = 24, n_sim = 50, seed = 1, verbose = FALSE
  )
  expect_s3_class(r, "pits_power_result")
  expect_gte(r$power, 0)
  expect_lte(r$power, 1)
})

test_that("export_results writes two files and returns paths", {
  r <- calculate_power(12, 12, 15, -3, 0, 2, 0.4, n_sim = 30, seed = 1)
  paths <- export_results(r, dir = tempdir(), prefix = "pits_test")
  expect_length(paths, 2L)
  expect_true(all(file.exists(paths)))
})

test_that("export_results handles a sweep data frame", {
  sw <- power_sweep(c(12, 24), n_pre = 12, baseline = 15,
                   level_change = -3, sigma = 2, rho = 0.4,
                   n_sim = 30, seed = 1, verbose = FALSE)
  paths <- export_results(sw, dir = tempdir(), prefix = "pits_sweep_test")
  expect_length(paths, 1L)
  expect_true(file.exists(paths))
})


# =============================================================================
#  plots.R — smoke tests (no visual inspection, just no errors)
# =============================================================================

test_that("plot_power_curve runs without error", {
  sw <- power_sweep(c(12, 24, 36), n_pre = 12, baseline = 15,
                   level_change = -3, sigma = 2, rho = 0.4,
                   n_sim = 30, seed = 1, verbose = FALSE)
  f <- tempfile(fileext = ".pdf")
  pdf(f); on.exit({ dev.off(); unlink(f) }, add = TRUE)
  expect_no_error(plot_power_curve(sw))
})

test_that("plot_power_heatmap runs without error", {
  g  <- build_param_grid(n_post = c(12, 24), level_change = c(-2, -3),
                         sigma = 2, rho = 0.4)
  gr <- run_power_grid(g, n_sim = 30, seed = 1, verbose = FALSE)
  f  <- tempfile(fileext = ".pdf")
  pdf(f); on.exit({ dev.off(); unlink(f) }, add = TRUE)
  expect_no_error(plot_power_heatmap(gr, x = "n_post", y = "level_change"))
})

test_that("plot_its_example returns a data frame invisibly", {
  f <- tempfile(fileext = ".pdf")
  pdf(f); on.exit({ dev.off(); unlink(f) }, add = TRUE)
  df <- plot_its_example(n_pre = 12, n_post = 12, baseline = 15,
                         level_change = -3, sigma = 2, rho = 0.4, seed = 1)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 24L)
})

test_that("diagnose_params returns list with params and residuals", {
  data("example_cfr_data", package = "PITS")
  f <- tempfile(fileext = ".pdf")
  pdf(f); on.exit({ dev.off(); unlink(f) }, add = TRUE)
  result <- diagnose_params(example_cfr_data)
  expect_type(result, "list")
  expect_named(result, c("params", "residuals"))
  expect_equal(result$params$n_pre, nrow(example_cfr_data))
})
