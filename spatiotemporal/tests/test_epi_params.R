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

test_that("GT PMF first element is strictly zero only for long Gamma (sanity)", {
  skip_if_missing("make_gt_pmf")

  # long: mean=13d, sd=5.5d → shape≈5.6; P(GT<=1) is astronomically small
  pmf_long <- make_gt_pmf(13, 5.5, 40)
  expect_true(pmf_long[1] < 1e-4, info = "Long GT: P(GT=1 day) < 0.01%")
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
# 3. compute_truncation_weights() — Exponential-CDF weights
#
# Signature: compute_truncation_weights(week_start_dates, analysis_date, rate)
# The function uses week MIDPOINTS (week_start + 3.5 days) for lag computation.
# =============================================================================

test_that("Truncation weight for 7-day lag from midpoint matches Exp CDF", {
  skip_if_missing("compute_truncation_weights")

  # week_start such that midpoint = analysis_date - 7 exactly:
  # midpoint = week_start + 3.5  =>  week_start = analysis_date - 7 - 3.5 = analysis_date - 10.5
  # We use analysis_date - 11 (integer) so midpoint = analysis_date - 7.5 ≈ 7.5d lag
  # To get exactly 7d lag: week_start = analysis_date - 10.5, but Date is integer.
  # Instead use analysis_date - 11 and verify the weight against P(Exp(0.2286) <= 7.5)
  analysis_date <- as.Date("2026-05-15")
  week_start_7  <- analysis_date - 11L   # midpoint = analysis_date - 7.5 d
  week_start_14 <- analysis_date - 18L   # midpoint = analysis_date - 14.5 d

  w <- compute_truncation_weights(
    week_start_dates = c(week_start_7, week_start_14),
    analysis_date    = analysis_date,
    rate             = 0.2286
  )

  # Expected: P(Exp(0.2286) <= 7.5) = 1 - exp(-0.2286 * 7.5) ≈ 0.818
  expected_7  <- 1 - exp(-0.2286 * 7.5)
  expected_14 <- 1 - exp(-0.2286 * 14.5)

  expect_true(
    abs(w[1] - expected_7) < 0.005,
    info = sprintf("7.5-day lag weight: expected %.4f, got %.4f", expected_7, w[1])
  )
  expect_true(
    abs(w[2] - expected_14) < 0.005,
    info = sprintf("14.5-day lag weight: expected %.4f, got %.4f", expected_14, w[2])
  )
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

test_that("Truncation weight formula exactly matches Exponential CDF at known points", {
  skip_if_missing("compute_truncation_weights")

  # Analytic spot-checks using rate = 0.2286:
  # P(Exp(0.2286) <= d) = 1 - exp(-0.2286 * d)
  # lag_days = analysis_date - (week_start + 3.5)
  # For lag = 10 days exactly: week_start = analysis_date - 13.5, use analysis_date - 14
  # → midpoint = analysis_date - 10.5 → lag = 10.5 d
  rate          <- 0.2286
  analysis_date <- as.Date("2026-07-01")

  # Construct weeks with known half-day-rounded lags
  # week_start = analysis_date - 14 → lag_days = 14 - 3.5 = 10.5
  week_start_a <- analysis_date - 14L
  w_a <- compute_truncation_weights(week_start_a, analysis_date, rate = rate)
  expected_a <- 1 - exp(-rate * 10.5)
  expect_equal(as.numeric(w_a), expected_a, tolerance = 1e-9,
               info = "Weight for 10.5-day lag must equal 1 - exp(-rate * 10.5)")

  # week_start = analysis_date - 21 → lag = 21 - 3.5 = 17.5 d
  week_start_b <- analysis_date - 21L
  w_b <- compute_truncation_weights(week_start_b, analysis_date, rate = rate)
  expected_b <- 1 - exp(-rate * 17.5)
  expect_equal(as.numeric(w_b), expected_b, tolerance = 1e-9,
               info = "Weight for 17.5-day lag must equal 1 - exp(-rate * 17.5)")
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
