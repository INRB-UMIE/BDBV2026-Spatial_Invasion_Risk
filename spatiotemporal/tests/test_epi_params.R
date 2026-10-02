# =============================================================================
# tests/test_epi_params.R — Unit tests for 02_epi_params.R functions
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Tests cover:
#   - make_gt_pmf()               generation time PMF correctness
#   - compute_truncation_weights() Exponential-CDF weight computation
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: skip gracefully when a function is unavailable
# ---------------------------------------------------------------------------
skip_if_missing <- function(fn_name) {
  if (!existsFunction(fn_name)) {
    skip(paste("Function", fn_name, "not available — source 02_epi_params.R first"))
  }
}

existsFunction <- function(fn_name) {
  exists(fn_name, mode = "function", envir = .GlobalEnv, inherits = TRUE)
}

# =============================================================================
# 1. make_gt_pmf() — generation time PMF
# =============================================================================

test_that("GT PMF sums to exactly 1 for all three profiles", {
  skip_if_missing("make_gt_pmf")

  for (profile_name in c("short", "medium", "long")) {
    p   <- GT_PROFILES[[profile_name]]
    pmf <- make_gt_pmf(p$mean, p$sd, p$max_tau)

    expect_true(
      is.numeric(pmf),
      info = paste(profile_name, ": GT PMF must be numeric")
    )
    expect_equal(
      length(pmf), p$max_tau,
      info = paste(profile_name, ": GT PMF length must equal max_tau")
    )
    # make_gt_pmf() internally calls stopifnot(abs(sum - 1) < 1e-9) so this
    # is guaranteed; we verify it from the test side too.
    expect_equal(
      sum(pmf), 1,
      tolerance = 1e-9,
      info = paste(profile_name, ": GT PMF must sum to 1")
    )
    expect_true(
      all(pmf >= 0),
      info = paste(profile_name, ": GT PMF must be non-negative")
    )
  }
})

test_that("GT PMF first element is near-zero for slow Gamma (medium profile)", {
  skip_if_missing("make_gt_pmf")

  # medium: mean=9d, sd=4.5d → shape=4, rate=0.444
  # P(Gamma(4, 0.444) <= 1) = pgamma(1, 4, 0.444) ≈ 0.00013; should be tiny
  pmf <- make_gt_pmf(9, 4.5, 30)
  expect_true(
    pmf[1] < 0.01,
    info = "P(GT=1 day) for mean-9d Gamma should be < 1%"
  )
})

test_that("GT PMF first element is far smaller for a long than a short generation time", {
  skip_if_missing("make_gt_pmf")

  # THRESHOLD CORRECTED (was: pmf_long[1] < 1e-4, which had failed since make_gt_pmf()
  # switched to the DOUBLE-INTERVAL-CENSORED discretisation). 1e-4 was calibrated for the
  # naive single-interval form, where P(GT = 1) = F(1) - F(0) = 1.85e-05 for this Gamma.
  # Under primary censoring the recorded day-difference is floor(U + D) with U ~ Unif(0,1),
  # so tau = 1 collects the mass with D in (1 - U, 2 - U) — an order of magnitude more,
  # 2.02e-04. The code is right (primarycensored and the base-R integral fallback agree to
  # 5e-15); the constant was stale. Asserting a RATIO rather than an absolute floor also
  # makes the test scale-free, so it keeps testing the intended property — a long GT puts
  # far less mass on a 1-day interval than a short one — instead of pinning a magic number
  # that any future change to the discretisation convention would silently invalidate.
  pmf_long   <- make_gt_pmf(13, 5.5, 40)   # shape ~ 5.6
  pmf_medium <- make_gt_pmf(9,  4.5, 40)   # shape = 4

  expect_lt(pmf_long[1], pmf_medium[1] / 10)   # observed ratio ~0.038
  expect_lt(pmf_long[1], 1e-3)
  expect_gt(pmf_long[1], 0)
})

test_that("GT PMF mode is near expected value for medium profile", {
  skip_if_missing("make_gt_pmf")

  # medium: mean=9, sd=4.5 → shape=4, rate=1/4.5*9 = 0.444
  # Continuous Gamma mode = (shape-1)/rate = 3/0.444 ≈ 6.75 d
  # Discretised mode should be tau=7 (±2 allowed for discretisation rounding)
  pmf      <- make_gt_pmf(9, 4.5, 30)
  mode_tau <- which.max(pmf)
  expect_true(
    mode_tau >= 5 & mode_tau <= 10,
    info = sprintf(
      "Medium GT PMF mode should be 5-10 days; got tau=%d (pmf=%.4f)",
      mode_tau, pmf[mode_tau]
    )
  )
})

test_that("GT PMF mode shifts right as mean increases across profiles", {
  skip_if_missing("make_gt_pmf")

  mode_short  <- which.max(make_gt_pmf(5.5,  2.0, 20))
  mode_medium <- which.max(make_gt_pmf(9.0,  4.5, 30))
  mode_long   <- which.max(make_gt_pmf(13.0, 5.5, 40))

  expect_true(
    mode_short < mode_medium,
    info = "Short GT mode must be earlier than medium GT mode"
  )
  expect_true(
    mode_medium < mode_long,
    info = "Medium GT mode must be earlier than long GT mode"
  )
})

test_that("make_gt_pmf rejects invalid inputs", {
  skip_if_missing("make_gt_pmf")

  expect_error(make_gt_pmf(mean = -1, sd = 4, max_tau = 30),
               info = "Negative mean must error")
  expect_error(make_gt_pmf(mean = 9, sd = 0, max_tau = 30),
               info = "Zero sd must error")
  expect_error(make_gt_pmf(mean = 9, sd = 4, max_tau = 0),
               info = "Zero max_tau must error")
})

test_that("GT PMF is correctly normalised when truncation removes tail mass", {
  skip_if_missing("make_gt_pmf")

  # With max_tau = 5 on a mean-9d Gamma, a lot of mass is truncated.
  # The result must still sum exactly to 1 (renormalisation).
  pmf_truncated <- make_gt_pmf(9, 4.5, 5)
  expect_equal(sum(pmf_truncated), 1, tolerance = 1e-9,
               info = "Hard-truncated GT PMF must renormalise to 1")
  expect_equal(length(pmf_truncated), 5L,
               info = "Hard-truncated GT PMF must have length = max_tau")
})

# =============================================================================
# 3. compute_truncation_weights() — fitted-delay completeness weights
#
# Signature: compute_truncation_weights(week_start_dates, analysis_date, delay, rate)
#
# The weight is E_o[P(delay <= analysis_date - o)] over the seven onset days o of the
# bucket, i.e. the MEAN of the delay CDF across the week — NOT the CDF evaluated once at
# the week midpoint (which, F being concave, overstates completeness). P(Delta <= L) is
# F(L + 0.5) under the daily-rounding convention the delay is fit under. Passing `rate`
# pins a plain Exp(rate); the default resolves the fitted delay in force.
#
# Reference implementation used by the analytic tests below:
.wk_weight_ref <- function(week_start, analysis_date, rate, week_days = 7L) {
  lags <- as.numeric(analysis_date - week_start) - seq.int(0L, week_days - 1L)
  cdf  <- ifelse(lags < 0, 0, 1 - exp(-rate * (lags + 0.5)))
  mean(cdf)
}
# =============================================================================

test_that("Truncation weight equals the mean Exp CDF across the bucket's onset days", {
  skip_if_missing("compute_truncation_weights")

  analysis_date <- as.Date("2026-05-15")
  week_start_7  <- analysis_date - 11L   # onset-day lags 11..5
  week_start_14 <- analysis_date - 18L   # onset-day lags 18..12

  w <- compute_truncation_weights(
    week_start_dates = c(week_start_7, week_start_14),
    analysis_date    = analysis_date,
    rate             = 0.2286
  )

  expect_equal(as.numeric(w[1]), .wk_weight_ref(week_start_7,  analysis_date, 0.2286),
               tolerance = 1e-12,
               info = "week-averaged weight must equal mean_j (1 - exp(-rate*(lag_j+0.5)))")
  expect_equal(as.numeric(w[2]), .wk_weight_ref(week_start_14, analysis_date, 0.2286),
               tolerance = 1e-12)

  # Jensen: F is concave, so the week-AVERAGED weight must sit strictly below F evaluated
  # at the week's MEAN lag. (Note this is NOT a comparison against the retired week-midpoint
  # form: that used week_start + 3.5, i.e. a lag a full day SHORTER than the true mean lag
  # of the seven integer onset days, so the two errors partly offset and the retired value
  # can land on either side.)
  mean_lag <- as.numeric(analysis_date - week_start_7) - 3 + 0.5
  expect_lt(as.numeric(w[1]), 1 - exp(-0.2286 * mean_lag))
})

test_that("Truncation weight increases monotonically with lag (older weeks = more complete)", {
  skip_if_missing("compute_truncation_weights")

  analysis_date <- as.Date("2026-06-01")
  # Five successive weeks, oldest to most recent
  week_starts <- analysis_date - c(35L, 28L, 21L, 14L, 7L)

  w <- compute_truncation_weights(week_starts, analysis_date, rate = 0.2286)

  expect_true(
    all(diff(w) < 0),   # as lag decreases, weight decreases
    info = "Weights must be monotonically decreasing as weeks get more recent"
  )
})

test_that("Truncation weights are in (0, 1] for all past weeks", {
  skip_if_missing("compute_truncation_weights")

  analysis_date <- as.Date("2026-06-15")
  past_weeks <- seq(as.Date("2026-05-01"), analysis_date - 7L, by = "week")

  w <- compute_truncation_weights(past_weeks, analysis_date, rate = 0.2286)

  valid <- !is.na(w)
  expect_true(
    all(w[valid] > 0 & w[valid] <= 1),
    info = "All observed-week weights must be in (0, 1]"
  )
})

test_that("Future weeks receive NA truncation weight", {
  skip_if_missing("compute_truncation_weights")

  analysis_date  <- as.Date("2026-06-01")
  future_week    <- analysis_date + 7L   # clearly in the future (midpoint > analysis_date)
  past_week      <- analysis_date - 14L

  w <- suppressWarnings(
    compute_truncation_weights(c(past_week, future_week), analysis_date, rate = 0.2286)
  )

  expect_false(is.na(w[1]), info = "Past week must have a non-NA weight")
  expect_true(is.na(w[2]),  info = "Future week must receive NA weight")
})

test_that("Truncation weight formula exactly matches the week-averaged Exponential CDF", {
  skip_if_missing("compute_truncation_weights")

  rate          <- 0.2286
  analysis_date <- as.Date("2026-07-01")

  # week_start = analysis_date - 14 -> onset-day lags 14, 13, ..., 8
  week_start_a <- analysis_date - 14L
  w_a <- compute_truncation_weights(week_start_a, analysis_date, rate = rate)
  expect_equal(as.numeric(w_a), .wk_weight_ref(week_start_a, analysis_date, rate),
               tolerance = 1e-12,
               info = "Weight must equal mean_j (1 - exp(-rate*(lag_j + 0.5))), j = 0..6")

  # week_start = analysis_date - 21 -> onset-day lags 21, 20, ..., 15
  week_start_b <- analysis_date - 21L
  w_b <- compute_truncation_weights(week_start_b, analysis_date, rate = rate)
  expect_equal(as.numeric(w_b), .wk_weight_ref(week_start_b, analysis_date, rate),
               tolerance = 1e-12)
})

test_that("Truncation weights default to the fitted delay, not the fixed lab Exp rate", {
  skip_if_missing("compute_truncation_weights")
  skip_if_missing("effective_onset_sample_delay")

  analysis_date <- as.Date("2026-07-31")
  week_start    <- analysis_date - 6L        # the current (most truncated) week

  dly <- effective_onset_sample_delay()
  w_default  <- compute_truncation_weights(week_start, analysis_date)
  w_explicit <- compute_truncation_weights(week_start, analysis_date, delay = dly)

  expect_equal(as.numeric(w_default), as.numeric(w_explicit), tolerance = 1e-12,
               info = "the default must BE effective_onset_sample_delay(), not a lab constant")

  # REGRESSION GUARD for the bug this replaced: whenever the fitted delay is slower than
  # the fixed lab reference, the default weight must be strictly below the lab weight —
  # i.e. the fold correction must be LARGER, not silently pinned to the fast lab delay.
  if (isTRUE(dly$source == "data") && dly$mean > 1 / DELAY_ONSET_SAMPLE_RATE) {
    w_lab <- compute_truncation_weights(week_start, analysis_date,
                                        rate = DELAY_ONSET_SAMPLE_RATE)
    expect_lt(as.numeric(w_default), as.numeric(w_lab))
  }
})

test_that("compute_truncation_weights rejects an ambiguous delay specification", {
  skip_if_missing("compute_truncation_weights")
  skip_if_missing("effective_onset_sample_delay")
  expect_error(
    compute_truncation_weights(as.Date("2026-07-01") - 14L, as.Date("2026-07-01"),
                               delay = effective_onset_sample_delay(), rate = 0.2),
    "not both"
  )
})

test_that("Higher rate gives higher weight for same lag (faster reporting)", {
  skip_if_missing("compute_truncation_weights")

  analysis_date <- as.Date("2026-06-15")
  week_start    <- analysis_date - 14L   # lag ≈ 10.5 d

  w_slow <- compute_truncation_weights(week_start, analysis_date, rate = 0.10)
  w_fast <- compute_truncation_weights(week_start, analysis_date, rate = 0.50)

  expect_true(
    as.numeric(w_fast) > as.numeric(w_slow),
    info = "Higher reporting rate → higher truncation weight for same lag"
  )
})
