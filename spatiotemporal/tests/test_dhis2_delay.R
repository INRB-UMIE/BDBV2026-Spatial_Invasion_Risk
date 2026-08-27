# =============================================================================
# tests/test_dhis2_delay.R — DHIS2 onset->sample delay estimation
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Exercises 04c_dhis2_delay_windows.R:
#   - .make_cens_df()  interval-censoring convention (d=0->[0,0.5]; d>0->[d-.5,d+.5])
#   - .fit_all_censored() recovers a known gamma from synthetic delays
#   - estimate_dhis2_onset_sample_delay() windows + fits + summarises
#   - write_onset_sample_long() emits the pipeline long format WITHOUT a `rate`
#     name collision (Exp-rate summary vs the family's native `rate`)
# All tests use SYNTHETIC data (seeded) so they are fast and deterministic.
# =============================================================================

skip_if_missing <- function(fn) {
  if (!exists(fn, mode = "function")) testthat::skip(paste0(fn, "() not available"))
}

# ---------------------------------------------------------------------------
# 1. Interval-censoring convention
# ---------------------------------------------------------------------------
test_that(".make_cens_df applies the daily-rounding censoring convention", {
  skip_if_missing(".make_cens_df")
  cd <- .make_cens_df(c(0L, 1L, 4L))
  # d = 0 -> [0, 0.5]; d > 0 -> [d - 0.5, d + 0.5]
  expect_equal(cd$left,  c(0, 0.5, 3.5))
  expect_equal(cd$right, c(0.5, 1.5, 4.5))
  # left must never be negative and left <= right (a valid censoring interval)
  expect_true(all(cd$left >= 0))
  expect_true(all(cd$left <= cd$right))
})

# ---------------------------------------------------------------------------
# 2. Censored MLE recovers a known gamma
# ---------------------------------------------------------------------------
test_that(".fit_all_censored recovers a known gamma delay (AIC-best) ", {
  skip_if_missing(".fit_all_censored")
  set.seed(1)
  # true gamma: shape 2, rate 0.4 -> mean 5 d; round to integer days (censoring)
  x <- pmax(round(stats::rgamma(1500, shape = 2, rate = 0.4)), 0L)
  fits <- .fit_all_censored(x, "test")
  expect_true(is.list(fits) && length(fits) >= 1)
  best <- .cens_best_fam(fits)
  expect_true(best %in% c("gamma", "weibull"))          # both are close; gamma should win
  # implied mean within ~0.5 d of the truth (5 d)
  mu <- .implied_mean(best, fits[[best]]$estimate)
  expect_true(is.finite(mu) && abs(mu - 5) < 0.6)
})

# ---------------------------------------------------------------------------
# 3. estimate_dhis2_onset_sample_delay: windowing + fit + summary
# ---------------------------------------------------------------------------
test_that("estimate_dhis2_onset_sample_delay windows and fits the onset->sample delay", {
  skip_if_missing("estimate_dhis2_onset_sample_delay")
  set.seed(2)
  n  <- 1200
  ob <- get0("OUTBREAK_START", ifnotfound = as.Date("2026-04-30"))
  # onsets spread over ~8 weeks from the outbreak start; delays ~ gamma(mean 5 d)
  onset  <- ob + sample(0:55, n, replace = TRUE)
  delay  <- pmax(round(stats::rgamma(n, shape = 2, rate = 0.4)), 0L)
  sample_d <- onset + delay
  ll <- data.frame(date_of_symptom_onset   = onset,
                   date_of_sample_collection = sample_d,
                   final_mve_case_classification = "confirmed_case",
                   stringsAsFactors = FALSE)
  fit <- estimate_dhis2_onset_sample_delay(ll, analysis_date = max(sample_d))
  expect_true(is.list(fit))
  expect_true(fit$best_family %in% c("gamma", "lnorm", "weibull", "exp"))
  expect_true(fit$n_fit > 0 && fit$n_fit <= n)
  expect_true(is.finite(fit$implied_mean) && fit$implied_mean > 2 && fit$implied_mean < 9)
  # rate summary is a positive Exponential-scale 1/mean
  expect_true(is.finite(fit$rate) && fit$rate > 0)
  expect_equal(unname(fit$rate), 1 / fit$implied_mean, tolerance = 1e-6)
  # windowing drops the last 5 days: no windowed onset should exceed max(sample) - 5
  expect_true(fit$trunc_date <= max(sample_d) - 5 + 1e-9)
})

# ---------------------------------------------------------------------------
# 4. write_onset_sample_long: pipeline format, no `rate` name collision
# ---------------------------------------------------------------------------
test_that("write_onset_sample_long emits long params with prefixed native params", {
  skip_if_missing("write_onset_sample_long")
  skip_if_missing("estimate_dhis2_onset_sample_delay")
  set.seed(3)
  n <- 800
  ob <- get0("OUTBREAK_START", ifnotfound = as.Date("2026-04-30"))
  onset <- ob + sample(0:45, n, replace = TRUE)
  sample_d <- onset + pmax(round(stats::rgamma(n, shape = 2, rate = 0.4)), 0L)
  ll  <- data.frame(date_of_symptom_onset = onset, date_of_sample_collection = sample_d)
  fit <- estimate_dhis2_onset_sample_delay(ll, analysis_date = max(sample_d))
  tmp <- tempfile(fileext = ".csv")
  out <- write_onset_sample_long(fit, tmp, source_label = "unit-test")
  on.exit(unlink(tmp), add = TRUE)
  expect_true(file.exists(tmp))
  # exactly ONE row with quantity == "rate" (the Exp-rate summary); the family's own
  # `rate` parameter is stored as `param_rate`, so deframe() is unambiguous.
  expect_equal(sum(out$quantity == "rate"), 1L)
  expect_true(any(grepl("^param_", out$quantity)))
  expect_true("family" %in% out$quantity && "implied_mean_fit" %in% out$quantity)
  # gamma native params, when gamma is selected, are param_shape / param_rate
  if (identical(fit$best_family, "gamma")) {
    expect_true(all(c("param_shape", "param_rate") %in% out$quantity))
    expect_false("shape" %in% out$quantity)   # unprefixed native names must NOT leak
  }
})
