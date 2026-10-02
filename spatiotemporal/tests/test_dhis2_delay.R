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

# =============================================================================
# Onset imputation draws from the SHARED delay resolver (01_data_prep.R)
# =============================================================================
# 00_config.R's stated invariant is that the imputation, the nowcast and the R(t) truncation
# model all read ONE delay estimator. The imputation used to bootstrap the raw, right-truncated
# empirical pairs (mean 6.81 d) while the other two used the truncation-corrected EpiDist
# marginal (7.67 d). This asserts the resolver is what a draw now reproduces.

test_that("the delay resolver is a truncation-corrected fit and .draw_dhis2_delay reproduces it", {
  skip_if_missing("effective_onset_sample_delay"); skip_if_missing(".draw_dhis2_delay")
  dp <- effective_onset_sample_delay()
  skip_if(!identical(dp$source, "data"), "no fitted delay params on disk")
  # ASSERT THE TRUNCATION CORRECTION, which the test name claims and nothing checked. Stripping
  # the epidist_* rows from the params CSV leaves source == "data" and silently falls back to
  # the interval-censored MLE (6.832 d vs the corrected 7.667 d) — the exact 0.84 d bias this
  # whole change removed — and every other assertion here still passed.
  expect_identical(dp$estimator, "epidist_marginal",
                   info = "the resolver must be the truncation-corrected EpiDist marginal")
  expect_gt(dp$mean, 7.0)   # the truncated MLE sits near 6.8; the corrected fit near 7.7
  expect_true(is.finite(dp$mean) && dp$mean > 0)
  set.seed(1)
  d <- .draw_dhis2_delay(dp, 2e5)
  expect_equal(mean(d), dp$mean, tolerance = 0.05,
               info = "draws must reproduce the resolver's mean, not a truncated one")
  expect_equal(stats::sd(d), dp$sd, tolerance = 0.10,
               info = "draws must reproduce the resolver's SD")
})

test_that("load_linelist() imputes onsets FROM the shared resolver, not the truncated pool", {
  skip_if_missing("load_linelist"); skip_if_missing("effective_onset_sample_delay")
  dp <- effective_onset_sample_delay()
  skip_if(!identical(dp$source, "data"), "no fitted delay params on disk")
  ll <- suppressWarnings(suppressMessages(load_linelist()))
  skip_if(!("onset_imputed" %in% names(ll)) || !any(ll$onset_imputed %in% TRUE),
          "no imputed onsets in this snapshot")
  # The realised gap between an imputed onset and its own sample date IS the drawn delay.
  imp <- ll[ll$onset_imputed %in% TRUE & !is.na(ll$date_index) &
            !is.na(ll$date_of_sample_collection), ]
  gap <- as.numeric(imp$date_of_sample_collection - imp$date_index)
  gap <- gap[is.finite(gap) & gap >= 0]
  skip_if(length(gap) < 100, "too few imputed onsets to test the distribution")
  # Must track the RESOLVER, not the right-truncated empirical pool (~6.8 d). The clamp at
  # DELAY_MAX_PLAUSIBLE_DAYS and rounding to whole days cost a little, hence the tolerance.
  #
  # UNDER THE GROWTH TILT THE TARGET IS NOT THE MARGINAL MEAN. ONSET_MODE = "growth_impute"
  # (the default) importance-resamples the resolver's draws with weights exp(-r*Delta), so the
  # realised mean delay is deliberately SHORTER than the marginal when the epidemic is growing
  # (r > 0) and longer when it is shrinking. Asserting equality with dp$mean would forbid the
  # very correction the mode exists to apply. What must still hold is that the draw comes from
  # the RESOLVER: the tilt reweights a gamma with mean ~7.7 d, it does not switch to the
  # truncated empirical pool, so the realised mean stays far from that pool's ~6.8 d unless the
  # tilt is large — and the direction must match the sign of r.
  .mode <- get0("ONSET_MODE", ifnotfound = "impute")
  if (identical(.mode, "growth_impute")) {
    # Wide band: the tilt's size depends on r, which moves with the data.
    expect_equal(mean(gap), dp$mean, tolerance = 0.35,
                 info = sprintf("growth-tilted gap mean %.2f d has drifted far from the resolver's %.2f d",
                                mean(gap), dp$mean))
  } else {
    expect_equal(mean(gap), dp$mean, tolerance = 0.12,
                 info = sprintf("imputed-onset gap mean %.2f d must track the resolver's %.2f d, not the truncated pool",
                                mean(gap), dp$mean))
  }
  expect_gt(mean(gap), 6.9)
})

test_that("the growth tilt moves imputed onsets in the direction of the growth rate", {
  skip_if_missing("load_linelist")
  .old <- get0("ONSET_MODE", envir = globalenv(), ifnotfound = NULL)
  on.exit(if (is.null(.old)) suppressWarnings(rm("ONSET_MODE", envir = globalenv()))
          else assign("ONSET_MODE", .old, envir = globalenv()), add = TRUE)
  .gap <- function(mode) {
    assign("ONSET_MODE", mode, envir = globalenv())
    ll <- tryCatch(suppressWarnings(suppressMessages(load_linelist())), error = function(e) NULL)
    if (is.null(ll) || !any(ll$onset_imputed %in% TRUE)) return(NA_real_)
    i <- ll[ll$onset_imputed %in% TRUE & !is.na(ll$date_index) &
              !is.na(ll$date_of_sample_collection), ]
    g <- as.numeric(i$date_of_sample_collection - i$date_index)
    mean(g[is.finite(g) & g >= 0])
  }
  g_flat <- .gap("impute"); g_tilt <- .gap("growth_impute")
  skip_if(!is.finite(g_flat) || !is.finite(g_tilt), "line list unavailable")
  # The BDBV 2026 record is growing, so the tilt must SHORTEN the imputed delay — i.e. move
  # imputed onsets LATER, toward their sample date. This is the whole point of the correction
  # (Lison et al. 2024): epidemic processes are not time-reversible, and the untilted backward
  # draw pushes onsets systematically too early while incidence rises.
  #
  # This guard caught two real defects. (1) The growth rate was fitted through the
  # right-truncated tail of the sample-date series, reading the reporting lag as a decline.
  # (2) More seriously, the weekly binning ran FORWARD from min(date), leaving a partial final
  # bin that held ~23 cases against ~1,800 in the preceding full week — which alone produced
  # r = -0.033/day for a growing epidemic and inverted the correction. Both are fixed in
  # .estimate_growth_rate(); if either regresses, g_tilt goes ABOVE g_flat and this fails.
  expect_lt(g_tilt, g_flat)
})
