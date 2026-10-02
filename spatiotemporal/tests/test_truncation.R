# =============================================================================
# tests/test_truncation.R — the estimated right-truncation, and the cap that
#   nearly clipped it
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# The nowcast used to be driven by the onset->SAMPLE delay (~7.7 d). The series
# it corrects is truncated by onset->APPEARANCE-IN-EXTRACT (~11.4 d measured
# against 54 archived vintages). Wrong leg, one direction: the most recent week
# — the week every forecast is issued FROM — was inflated 2.63x where 5.92x is
# right. epinow2_truncation() (02_epi_params.R section 4b) fits the right leg,
# per regime, and every pipeline caller now passes it explicitly.
#
# These tests pin the three things that can silently undo that:
#   1. the adapter that turns a fitted EpiNow2 dist into a delay spec,
#   2. the multiplier cap, which was a fixed 5x and so would have CLIPPED the
#      5.92x correction to 5x — defeating the fix at exactly the week it was
#      meant to fix, behind a warning,
#   3. the wrong-leg fallback staying loud.
# =============================================================================

.skip_unless <- function(fn) if (!exists(fn, mode = "function")) testthat::skip(paste0(fn, "() unavailable"))

# A LogNormal with the moments of the fitted extract-regime truncation
# (meanlog 2.0280, sdlog 0.7150 -> mean 9.81 d). Built by hand rather than read
# from the cache so the test pins the arithmetic even on a machine that has
# never run the fit.
.trunc_fixture <- function(meanlog = 2.0280, sdlog = 0.7150) {
  list(family = "lnorm", params = c(meanlog = meanlog, sdlog = sdlog),
       mean = exp(meanlog + sdlog^2 / 2),
       sd   = sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2)),
       rate = 1 / exp(meanlog + sdlog^2 / 2),
       estimator = "test_fixture", source = "fitted_truncation")
}

test_that(".trunc_as_delay_spec returns a spec compute_truncation_weights() will accept", {
  .skip_unless(".trunc_as_delay_spec")
  # NULL in, NULL out — the callers rely on this to fall back rather than error.
  expect_null(.trunc_as_delay_spec(NULL, "extract"))

  d <- .trunc_fixture()
  # compute_truncation_weights() VALIDATES `rate` even for a lnorm it never reads it for;
  # a spec without a positive finite rate is rejected there, so the adapter must set one.
  expect_true(is.numeric(d$rate) && is.finite(d$rate) && d$rate > 0)
  expect_equal(d$rate, 1 / d$mean)
  # LogNormal moments, not meanlog/sdlog passed through as if they were mean/sd.
  expect_equal(d$mean, exp(2.0280 + 0.7150^2 / 2), tolerance = 1e-9)
  expect_true(d$mean > exp(2.0280))  # E[X] > median for a LogNormal
  expect_silent(invisible(suppressMessages(
    compute_truncation_weights(as.Date("2026-09-01"), as.Date("2026-09-14"), delay = d))))
})

test_that("a longer truncation gives a strictly smaller weight and a bigger correction", {
  .skip_unless("compute_truncation_weights")
  wk <- as.Date("2026-09-08"); ad <- as.Date("2026-09-14")
  short <- suppressMessages(compute_truncation_weights(wk, ad, delay = .trunc_fixture(1.70, 0.7150)))
  long  <- suppressMessages(compute_truncation_weights(wk, ad, delay = .trunc_fixture(2.30, 0.7150)))
  expect_true(is.finite(short) && is.finite(long))
  expect_true(long < short)          # longer delay -> less of the week is in yet
  expect_true(long > 0 && short <= 1)
})

test_that("the multiplier cap is derived from the delay and does NOT clip the fitted correction", {
  .skip_unless("apply_nowcast_correction")
  # Enough stable history that the cap is switched on at all (it is Inf until there
  # are stable weeks), plus the anchor week that needs the deep correction.
  # THE ANCHOR WEEK STARTS AT ad - 6, NOT AT ad. The weekly grid is anchored so the final
  # week ENDS on the analysis date (00_config.R WEEK_ANCHOR); a week starting ON ad has mean
  # lag ad - (ad + 3) = -3 and is dropped as a future week, leaving trunc_weight NA. Getting
  # this wrong is how the first draft of this test "passed" nothing.
  ad    <- as.Date("2026-09-14")
  weeks <- seq(ad - 6 - 7 * 11, ad - 6, by = "week")
  zw <- tibble::tibble(
    health_zone = rep(c("Bunia", "Nizi"), each = length(weeks)),
    week_start  = rep(weeks, 2L),
    confirmed   = rep(10L, 2L * length(weeks)),
    suspected   = rep(10L, 2L * length(weeks))
  )
  d   <- .trunc_fixture()
  out <- suppressWarnings(suppressMessages(
    apply_nowcast_correction(zw, analysis_date = ad, min_lag_days = 0, delay = d)))

  anchor <- out[out$week_start == max(weeks) & out$health_zone == "Bunia", , drop = FALSE]
  expect_equal(nrow(anchor), 1L)
  expect_true(is.finite(anchor$trunc_weight))
  realised <- as.numeric(anchor$confirmed_nc / anchor$confirmed)
  implied  <- as.numeric(1 / anchor$trunc_weight)
  # THE POINT OF THIS TEST. Under the fitted truncation the anchor week needs ~5.9x.
  # The retired fixed 5x cap would have clipped it; the delay-derived cap must not.
  expect_true(implied > 5,
              info = "fixture should need MORE than the retired 5x cap, else this proves nothing")
  expect_equal(realised, implied, tolerance = 1e-6,
               info = "the anchor-week correction was clipped by the multiplier cap")
  # The live value, pinned: this fixture carries the fitted extract-regime moments, and the
  # deployed anchor week measured 0.1688 -> 5.923x on 2026-09-21.
  expect_equal(as.numeric(anchor$trunc_weight), 0.168829, tolerance = 1e-5)
  expect_equal(as.numeric(realised), 5.923152, tolerance = 1e-5)
})

test_that("no admissible week is EVER clipped, at any depth the grid can reach", {
  .skip_unless("apply_nowcast_correction")
  # A week survives only if lag_days = ad - (week_start + 3) is both > 0 and >= min_lag_days,
  # so with min_lag_days = 0 the deepest reachable week starts at ad - 4 (lag 1), not ad - 3
  # (lag 0, suppressed as future). The live grid is week-aligned and only ever reaches ad - 6,
  # but an unaligned origin reaches deeper and that is where clipping would first appear —
  # the retired fixed 5x cap was justified by a "4.577x deepest reachable" figure computed
  # under the retired delay, which the fitted truncation invalidates.
  ad <- as.Date("2026-09-14")
  d  <- .trunc_fixture()
  for (off in c(6L, 5L, 4L)) {
    last  <- ad - off
    weeks <- seq(last - 7 * 11, last, by = "week")
    zw <- tibble::tibble(
      health_zone = rep("Bunia", length(weeks)), week_start = weeks,
      confirmed = rep(10L, length(weeks)), suspected = rep(10L, length(weeks))
    )
    out <- suppressWarnings(suppressMessages(
      apply_nowcast_correction(zw, analysis_date = ad, min_lag_days = 0, delay = d)))
    a <- out[out$week_start == max(weeks), , drop = FALSE]
    expect_true(is.finite(a$trunc_weight),
                info = sprintf("week starting ad-%d should be admissible, not NA", off))
    expect_equal(as.numeric(a$confirmed_nc / a$confirmed), as.numeric(1 / a$trunc_weight),
                 tolerance = 1e-6,
                 info = sprintf("correction clipped by the cap at depth ad-%d", off))
  }
})

test_that("the cap is loose enough for the deepest week but still catches a pathological weight", {
  .skip_unless("compute_truncation_weights")
  d <- .trunc_fixture()
  # The deepest admissible week contributes lags 0..4 over a denominator of 7.
  w_deepest <- sum(delay_cdf(d, seq(0L, 4L) + 0.5)) / 7
  cap <- max(5, ceiling(1.5 / w_deepest))
  expect_true(cap > 1 / w_deepest, info = "the cap must not bind on the deepest legitimate week")
  # ...and must still be far below what the 1e-6 weight floor in compute_truncation_weights()
  # would otherwise produce, which is the only thing the cap is now for.
  expect_true(cap < 1e6)
})

test_that("the delay-derived cap reproduces the historic 4.577x bound on the retired delay", {
  .skip_unless("compute_truncation_weights"); .skip_unless("effective_onset_sample_delay")
  dly <- tryCatch(effective_onset_sample_delay(), error = function(e) NULL)
  if (is.null(dly) || !is.finite(dly$mean)) testthat::skip("delay params unavailable")
  # 04_nowcasting.R documented "the deepest reachable non-suppressed multiplier is 4.577x"
  # under the onset->sample delay. Recomputing it from the depth formula is an independent
  # check that the formula identifies the right week: get the depth wrong by one day and
  # this lands on 3.384x (ad-5) or 6.0x+ (ad-3), not 4.577x.
  expect_equal(7 / sum(delay_cdf(dly, seq(0L, 4L) + 0.5)), 4.577, tolerance = 2e-3)
})

test_that("omitting `delay` warns that it is the wrong leg", {
  .skip_unless("compute_truncation_weights")
  if (!exists(".nowcast_state")) testthat::skip(".nowcast_state unavailable")
  # The warning is once-per-session by design (the LFO would otherwise emit one per
  # fold); clear the flag so this test sees it regardless of what ran first.
  old <- .nowcast_state$warned_default_delay
  on.exit(.nowcast_state$warned_default_delay <- old, add = TRUE)
  rm(list = intersect("warned_default_delay", ls(.nowcast_state)), envir = .nowcast_state)

  expect_warning(
    suppressMessages(compute_truncation_weights(as.Date("2026-09-01"), as.Date("2026-09-14"))),
    "wrong leg")
  # ...and exactly once, not once per call.
  expect_warning(
    suppressMessages(compute_truncation_weights(as.Date("2026-09-01"), as.Date("2026-09-14"))),
    NA)
})

test_that("the regime is part of the EXTRACT cache key, so two panels cannot alias", {
  .skip_unless(".trunc_fingerprint")
  # Only the extract regime has a panel and a cache; the as-of regime is a deterministic
  # sub-second cohort fit with neither. The fingerprint must still carry the regime so a
  # future second panel-based regime cannot be served the extract fit.
  panel <- list(
    data.frame(date = seq(as.Date("2026-07-01"), as.Date("2026-08-01"), by = "day"), confirm = 1L),
    data.frame(date = seq(as.Date("2026-07-01"), as.Date("2026-08-08"), by = "day"), confirm = 1L)
  )
  expect_false(identical(.trunc_fingerprint(panel, "extract"),
                         .trunc_fingerprint(panel, "asof")))
  expect_identical(.trunc_fingerprint(panel, "extract"), .trunc_fingerprint(panel, "extract"))
  # The estimator is keyed too: a panel alone must not decide the fit.
  expect_true(grepl("negbin", .trunc_fingerprint(panel, "extract"), fixed = TRUE))
})

test_that("the as-of cohort fit reproduces the completeness it was measured from", {
  .skip_unless(".trunc_fit_asof_cohort")
  # Synthetic line list with a KNOWN delay, right-truncated the way the real one is: a record
  # is present only once its observation date has arrived. The cohort estimator must be immune
  # to that truncation by construction — for lag L it only ever uses onsets at least L days
  # old, for which "observed within L days?" is fully determined.
  set.seed(11L)
  end <- as.Date("2026-09-07"); start <- end - 120L
  n <- 6000L
  ons <- start + sample(0:120, n, replace = TRUE)
  dly <- stats::rgamma(n, shape = 1.05, rate = 0.11)
  obs <- ons + round(dly)
  vis <- obs <= end & obs >= ons
  ll <- data.frame(confirmed = TRUE,
                   date_of_symptom_onset = as.character(ons[vis]),
                   date_of_sample_collection = as.character(obs[vis]),
                   stringsAsFactors = FALSE)
  # .trunc_fit_asof_cohort() reads the observation date through linelist_observation_date();
  # stub it to the sample column, which is what it resolves to for 100% of confirmed records.
  old_fn <- if (exists("linelist_observation_date", envir = globalenv()))
    get("linelist_observation_date", envir = globalenv()) else NULL
  on.exit({ if (!is.null(old_fn)) assign("linelist_observation_date", old_fn, envir = globalenv()) },
          add = TRUE)
  assign("linelist_observation_date",
         function(ll, end, caller = "test") as.Date(ll$date_of_sample_collection),
         envir = globalenv())

  d <- suppressMessages(suppressWarnings(.trunc_fit_asof_cohort(ll, end, start)))
  expect_false(is.null(d))
  sp <- .trunc_as_delay_spec(d, "asof")
  expect_identical(sp$family, "gamma")

  # The measured anchor-week completeness, computed the same cohort way, with NO model.
  o  <- as.Date(ll$date_of_symptom_onset)
  dl <- as.numeric(as.Date(ll$date_of_sample_collection) - o)
  measured <- mean(vapply(0:6, function(L) mean(dl[o <= end - L] <= L), numeric(1)))
  fitted_w <- suppressMessages(compute_truncation_weights(end - 6, end, delay = sp))
  # Within 5%: the gamma must track what was measured, not merely be plausible. On the live
  # line list this agreement is -1.24% (0.2902 fitted vs 0.2939 measured).
  expect_equal(as.numeric(fitted_w), measured, tolerance = 0.05)

  # ...and it must beat the naive estimate that ignores truncation, which is the failure mode
  # this estimator exists to avoid.
  naive <- mean(vapply(0:6, function(L) mean(dl <= L), numeric(1)))
  expect_true(abs(fitted_w - measured) < abs(naive - measured))
})

test_that("the as-of regime needs no archive and no Stan", {
  .skip_unless("epinow2_truncation")
  # Regression: the as-of regime was briefly routed through estimate_truncation() on
  # deterministic pseudo-vintages, whose negbin dispersion was unidentified -- 129 s/chain,
  # then failure. It must not acquire a panel dependency again.
  expect_false(exists(".trunc_panel_asof", mode = "function"))
  expect_true(exists(".trunc_fit_asof_cohort", mode = "function"))
})
