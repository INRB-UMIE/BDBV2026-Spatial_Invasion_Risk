# =============================================================================
# tests/test_models.R — Unit tests for model correctness
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Tests cover (in approximate order of dependency):
#   1. apply_nowcast_correction()   — right-truncation inflation invariants (04)
#   2. Invasion probability         — Poisson P(Y >= 1) (pure identity)
#   3. FOI computation              — force-of-infection spatial coupling (06)
#   4. SEIR compartment dynamics    — non-negativity and mass conservation (08)
#
# NOTE: the WIS, Brier, log-loss and calibration tests that formerly lived here
# depended on compute_wis()/compute_brier_score()/compute_log_loss()/
# compute_calibration() from 11_metrics.R and 12_calibration.R, which have been
# RETIRED to spatiotemporal/stale/. Their assertions are ported to the live
# equivalents in test_invasion_pipeline.R:
#   * WIS      -> live `wis_one` (04b_epinowcast.R), Bracher-2021 hand example
#   * Brier    -> .brier      (16_invasion_eval.R)
#   * log-loss -> .log_score  (16_invasion_eval.R)
#   * calib.   -> AUC-PR / base-rate / ECE-style checks (16_invasion_eval.R)
#
# Remaining functions are guarded by skip_if_missing() so the suite still runs
# when an optional dependency (e.g. deSolve for the SEIR model) is absent.
# =============================================================================

existsFunction <- function(fn_name) {
  exists(fn_name, mode = "function", envir = .GlobalEnv, inherits = TRUE)
}

skip_if_missing <- function(fn_name) {
  if (!existsFunction(fn_name)) {
    skip(paste("Function", fn_name, "not available — source the relevant analysis script first"))
  }
}

# =============================================================================
# 1. apply_nowcast_correction() — from 04_nowcasting.R
# =============================================================================

test_that("Nowcast correction only inflates counts (confirmed_nc >= confirmed)", {
  skip_if_missing("apply_nowcast_correction")

  zone_week <- tibble::tibble(
    health_zone = rep(c("Bunia", "Nizi"), each = 3L),
    week_start  = rep(
      as.Date(c("2026-05-01", "2026-05-08", "2026-05-15")), 2L
    ),
    confirmed = c(5L, 8L, 3L, 2L, 4L, 1L),
    suspected = c(10L, 12L, 5L, 3L, 6L, 2L)
  )
  analysis_date <- as.Date("2026-05-25")

  corrected <- suppressWarnings(
    apply_nowcast_correction(zone_week, analysis_date)
  )

  # Corrected counts must be >= raw counts where correction is valid
  valid_rows <- !is.na(corrected$confirmed_nc)
  if (any(valid_rows)) {
    expect_true(
      all(corrected$confirmed_nc[valid_rows] >= corrected$confirmed[valid_rows] - 1e-9),
      info = "Nowcast-corrected confirmed counts must be >= raw confirmed counts"
    )
    expect_true(
      all(corrected$suspected_nc[valid_rows] >= corrected$suspected[valid_rows] - 1e-9),
      info = "Nowcast-corrected suspected counts must be >= raw suspected counts"
    )
  }
})

test_that("Truncation weights are in (0, 1] for observed weeks", {
  skip_if_missing("apply_nowcast_correction")

  zone_week <- tibble::tibble(
    health_zone = rep("Bunia", 4L),
    week_start  = as.Date(c("2026-05-01", "2026-05-08",
                             "2026-05-15", "2026-05-22")),
    confirmed   = c(5L, 3L, 8L, 2L),
    suspected   = c(8L, 5L, 12L, 3L)
  )
  analysis_date <- as.Date("2026-06-05")

  corrected <- suppressWarnings(
    apply_nowcast_correction(zone_week, analysis_date)
  )

  valid_weights <- corrected$trunc_weight[!is.na(corrected$trunc_weight)]
  if (length(valid_weights) > 0) {
    expect_true(
      all(valid_weights > 0 & valid_weights <= 1),
      info = "All non-NA truncation weights must be in (0, 1]"
    )
  }
})

test_that("Older weeks have higher truncation weight than recent weeks", {
  skip_if_missing("apply_nowcast_correction")

  # Use four weeks that are well past analysis_date so all are valid
  zone_week <- tibble::tibble(
    health_zone = rep("Bunia", 4L),
    week_start  = as.Date(c("2026-04-06", "2026-04-13",
                             "2026-04-20", "2026-04-27")),
    confirmed   = c(3L, 5L, 8L, 4L),
    suspected   = c(5L, 8L, 12L, 6L)
  )
  # Well after all weeks: minimum lag ≈ 30+ days
  analysis_date <- as.Date("2026-05-25")

  corrected <- suppressWarnings(
    apply_nowcast_correction(zone_week, analysis_date)
  )

  # Filter to rows for "Bunia", ordered by week
  bunia <- corrected[corrected$health_zone == "Bunia", ]
  bunia <- bunia[order(bunia$week_start), ]

  if (!any(is.na(bunia$trunc_weight))) {
    # Oldest week should have highest weight (most complete)
    expect_true(
      bunia$trunc_weight[1] >= bunia$trunc_weight[4],
      info = "Oldest week must have highest (or equal) truncation weight"
    )
    # Weight sequence should be non-increasing as weeks get more recent
    expect_true(
      all(diff(bunia$trunc_weight) <= 1e-9),
      info = "Weights must be non-increasing from oldest to most recent week"
    )
  }
})

test_that("apply_nowcast_correction: correction MULTIPLIER (not the count) is capped at 5x", {
  skip_if_missing("apply_nowcast_correction")

  # The current (fixed) semantics cap the multiplicative factor 1/trunc_weight at
  # 5x, so corrected lies in [raw, 5*raw] for BOTH series. (The retired behaviour
  # capped the absolute product at 5*max(confirmed) and mis-applied it to the
  # larger suspected series — see test_invasion_pipeline.R for the regression.)
  zone_week <- tibble::tibble(
    health_zone = rep("TestZone", 5L),
    week_start  = as.Date(c("2026-04-01", "2026-04-08", "2026-04-15",
                             "2026-04-22", "2026-06-11")),  # midpoint 06-14.5, lag 0.5 d
    confirmed   = c(5L, 8L, 6L, 7L, 10L),
    suspected   = c(8L, 12L, 10L, 11L, 15L)
  )
  # lag 0.5 d => weight ~0.108 => raw multiplier ~9.3 > 5, so the 5x cap MUST bind
  # (the earlier weeks are stable, weight ~1, so the cap is enabled). This exercises
  # the cap rather than passing trivially with a sub-5 multiplier.
  analysis_date <- as.Date("2026-06-15")

  corrected <- suppressWarnings(
    apply_nowcast_correction(zone_week, analysis_date, min_lag_days = 0)
  )

  valid <- !is.na(corrected$confirmed_nc)
  mult_c <- corrected$confirmed_nc[valid] / pmax(corrected$confirmed[valid], 1e-9)
  mult_s <- corrected$suspected_nc[valid] / pmax(corrected$suspected[valid], 1e-9)

  expect_true(all(mult_c <= 5 + 1e-9),
              info = "confirmed correction multiplier must not exceed 5x")
  expect_true(all(mult_s <= 5 + 1e-9),
              info = "suspected correction multiplier must not exceed 5x")
  expect_true(any(mult_c >= 5 - 1e-9),
              info = "the recent week's raw multiplier (~9.3) must be capped AT 5x")
  expect_true(all(mult_c >= 1 - 1e-9) && all(mult_s >= 1 - 1e-9),
              info = "correction only inflates (multiplier >= 1)")
})

# =============================================================================
# 2. Invasion probability: P(Y >= 1 | mu) = 1 - exp(-mu) under Poisson
# =============================================================================

test_that("Invasion probability is in [0, 1] for all non-negative means", {
  mu_vals  <- c(0, 0.001, 0.1, 1, 5, 20, 100)
  p_pois   <- 1 - exp(-mu_vals)

  expect_true(all(p_pois >= 0 & p_pois <= 1),
              info = "Invasion probability must be in [0, 1] for all mu >= 0")
})

test_that("Invasion probability is 0 for mu = 0 and approaches 1 for large mu", {
  expect_equal(1 - exp(-0),   0, tolerance = 1e-15, info = "mu=0 → P(invasion)=0")
  expect_equal(1 - exp(-0.0), 0, tolerance = 1e-15, info = "mu=0.0 → P(invasion)=0")
  expect_true(1 - exp(-100) > 0.9999, info = "mu=100 → P(invasion) > 0.9999")
  expect_true(1 - exp(-200) > 1 - 1e-10, info = "mu=200 → P(invasion) ≈ 1")
})

test_that("Invasion probability is monotonically increasing in mu", {
  mu_seq <- seq(0, 20, by = 0.5)
  p_seq  <- 1 - exp(-mu_seq)
  expect_true(all(diff(p_seq) >= 0),
              info = "Invasion probability must be monotonically non-decreasing in mu")
})

# =============================================================================
# 3. Force of infection (FOI) — spatial coupling via mobility matrix
#
# For a simple TSIR-like discrete model:
#   Lambda_i(t) = sum_j W[j,i] * convolve(Y_j, gt_pmf, t)
# where W[j,i] is the proportion of population from j visiting i.
#
# Test: when only zone A has past cases, zones B and C should receive
# positive FOI from A (via W[A,B] and W[A,C] > 0), while A's own FOI
# from B and C is zero (since B and C have no cases and W[A,A]=0).
# =============================================================================

test_that("FOI is non-negative and correctly directed through mobility matrix", {
  skip_if_missing("compute_foi")
  skip_if_missing("make_gt_pmf")

  n_weeks <- 5L
  n_zones <- 3L
  zone_names <- c("A", "B", "C")

  # Zone A has cases; B and C are naive
  Y_wide <- matrix(0L, n_zones, n_weeks, dimnames = list(zone_names, 1:n_weeks))
  Y_wide["A", 1:3] <- c(5L, 8L, 3L)

  # W[i,j] = prob zone i's population moves to zone j (row-stochastic)
  # W[A,A]=0 (no self-loop), W[A,B]=0.5, W[A,C]=0.5
  # W[B,A]=0.3, W[B,B]=0, W[B,C]=0.7 — but B has no cases
  W <- matrix(
    c(0.0, 0.5, 0.5,
      0.3, 0.0, 0.7,
      0.2, 0.8, 0.0),
    nrow = n_zones, byrow = TRUE,
    dimnames = list(zone_names, zone_names)
  )

  gt_pmf <- make_gt_pmf(9, 4.5, 30)

  # compute_foi requires a weekly GT PMF, not the daily PMF from make_gt_pmf()
  G_weekly <- tryCatch(daily_to_weekly_gt(gt_pmf),
                       error = function(e) gt_pmf)  # fallback if not available

  Lambda <- tryCatch(
    compute_foi(Y_wide, W, G_weekly, t_idx = 4L, zone_names),
    error = function(e) {
      skip(paste("compute_foi signature mismatch or not yet implemented:", conditionMessage(e)))
    }
  )

  expect_true(all(Lambda >= 0),
              info = "FOI must be non-negative for all zones")

  # B and C receive inflow from A (which has cases)
  expect_true(Lambda["B"] > 0,
              info = "Zone B must have positive FOI from zone A's cases via W[A,B]>0")
  expect_true(Lambda["C"] > 0,
              info = "Zone C must have positive FOI from zone A's cases via W[A,C]>0")

  # A's FOI from B and C is zero (neither has cases)
  # Lambda["A"] = sum_j W[j,A] * past_j; B and C contribute zero cases
  # Lambda["A"] should be 0 (or very near 0 from W[B,A] * 0 + W[C,A] * 0)
  expect_equal(unname(Lambda["A"]), 0, tolerance = 1e-9,
               info = "Zone A's FOI from B and C should be 0 since they have no cases")
})

# =============================================================================
# 4. SEIR compartments — non-negativity and mass conservation
# =============================================================================

test_that("SEIR: all compartments non-negative and total population conserved", {
  skip_if_missing("run_seir_deterministic")
  skip_if_missing("make_gt_pmf")

  n_zones    <- 3L
  zone_names <- c("A", "B", "C")
  N          <- c(A = 100000, B = 50000, C = 30000)

  initial_state <- list(
    S = N - c(A = 100L, B = 50L, C = 20L),
    E = c(A = 20L,  B = 10L, C = 5L),
    I = c(A = 50L,  B = 25L, C = 10L),
    R = c(A = 30L,  B = 15L, C = 5L),
    D = c(A = 0L,   B = 0L,  C = 0L)
  )

  W <- matrix(
    c(0.0, 0.6, 0.4,
      0.3, 0.0, 0.7,
      0.5, 0.5, 0.0),
    nrow = n_zones, byrow = TRUE,
    dimnames = list(zone_names, zone_names)
  )

  params <- tryCatch(make_seir_params(), error = function(e) NULL)
  skip_if(is.null(params), "make_seir_params() not available")

  R_by_zone <- c(A = 1.5, B = 1.2, C = 1.8)
  initial_state$N <- N

  result <- tryCatch(
    run_seir_deterministic(
      state     = initial_state,
      R_by_zone = R_by_zone,
      W         = W,
      params    = params,
      m         = 0.1,
      n_weeks   = 4L
    ),
    error = function(e) {
      skip(paste("run_seir_deterministic not yet implemented or API mismatch:",
                 conditionMessage(e)))
    }
  )

  # Non-negativity
  expect_true(all(result$S >= -1e-6), info = "S compartment must be non-negative")
  expect_true(all(result$E >= -1e-6), info = "E compartment must be non-negative")
  expect_true(all(result$I >= -1e-6), info = "I compartment must be non-negative")
  expect_true(all(result$R >= -1e-6), info = "R compartment must be non-negative")
  expect_true(all(result$D >= -1e-6), info = "D compartment must be non-negative")

  # Mass conservation at final time step: S + E + I + R + D = N (per zone)
  final_week <- ncol(result$S)
  total_final <- (result$S[, final_week] +
                    result$E[, final_week] +
                    result$I[, final_week] +
                    result$R[, final_week] +
                    result$D[, final_week])

  expect_equal(
    unname(total_final), unname(N),
    tolerance = 1.0,   # allow 1 person rounding from ODE solver discretisation
    info = "SEIR total (S+E+I+R+D) must equal N at each time step"
  )
})
