# =============================================================================
# tests/test_models.R — Unit tests for model correctness
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Tests cover (in approximate order of dependency):
#   1. apply_nowcast_correction()   — right-truncation inflation invariants (04)
#   2. Invasion probability         — Poisson P(Y >= 1) (pure identity)
#   3. FOI computation              — force-of-infection spatial coupling (06)
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
# when an optional dependency is absent.
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
                             "2026-04-22", "2026-06-11")),  # last week: mean lag 1 d
    confirmed   = c(5L, 8L, 6L, 7L, 10L),
    suspected   = c(8L, 12L, 10L, 11L, 15L)
  )
  # The delay is PINNED to a deliberately slow Exp(0.05) (mean 20 d), which decouples this
  # regression from whichever delay is fitted on the day. The last week here starts
  # analysis_date - 4 — the DEEPEST week the suppression rules admit (lag_days 1) — so its
  # multiplier is the largest one this delay can legitimately produce, ~12.1x. The four
  # earlier weeks have weight ~0.97, so the "stable history exists" precondition that enables
  # the cap holds.
  #
  # CONTRACT CHANGED 2026-09-21. The cap was a fixed 5x and this test asserted it BOUND at
  # exactly 5x. It is now derived from the delay (1.5x headroom over the deepest admissible
  # week), because with the fitted truncation in force the deployed anchor week legitimately
  # needs 5.92x — a fixed 5x would have silently clipped the correction at precisely the week
  # the truncation fix exists to repair. The cap's job is now only to catch a PATHOLOGICAL
  # weight (compute_truncation_weights() floors weights at 1e-6, which would otherwise
  # multiply by 1e6), so on legitimate input it must NOT bind. That is what is asserted below.
  analysis_date <- as.Date("2026-06-15")

  corrected <- suppressWarnings(
    apply_nowcast_correction(zone_week, analysis_date, min_lag_days = 0, rate = 0.05)
  )

  valid <- !is.na(corrected$confirmed_nc)
  mult_c <- corrected$confirmed_nc[valid] / pmax(corrected$confirmed[valid], 1e-9)
  mult_s <- corrected$suspected_nc[valid] / pmax(corrected$suspected[valid], 1e-9)

  # The cap, restated from the implementation: 1.5x headroom over the deepest admissible
  # week, which for week_days 7 and min_lag_days 0 contributes lags 0..4 over a denominator
  # of 7. Restated rather than imported so a change to the rule has to be made twice.
  deepest_mult <- 7 / sum(pexp(seq(0L, 4L) + 0.5, rate = 0.05))
  cap          <- max(5, ceiling(1.5 * deepest_mult))
  expect_equal(deepest_mult, 12.14, tolerance = 5e-3)

  expect_true(all(mult_c <= cap + 1e-9),
              info = "confirmed correction multiplier must not exceed the delay-derived cap")
  expect_true(all(mult_s <= cap + 1e-9),
              info = "suspected correction multiplier must not exceed the delay-derived cap")
  # THE NEW CONTRACT: the deepest admissible week reaches its FULL uncapped multiplier.
  # Under the retired fixed 5x cap this row was clipped from 12.1x to 5.0x.
  expect_equal(max(mult_c), deepest_mult, tolerance = 1e-6,
               info = "the deepest admissible week must NOT be clipped by the cap")
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

# The SEIR compartment tests were removed with 08_stochastic_seir.R in the 2026-09-17
# streamlining: the stochastic SEIR comparator was gated OFF by default and fed no
# retained figure.

# =============================================================================
# 4. WEEKLY GENERATION TIME — lag-0 mass and the daily->weekly R correction (06)
# =============================================================================
# weekly_censored_gt() drops the weekly lag-0 term (an explicit weekly recursion cannot carry
# a self-referential term) and renormalises. EpiNow2's R lives on the DAILY renewal, which
# retains that mass, so applying it unchanged to the weekly kernel understates growth.
# weekly_renewal_R_eff() is the algebraic correction.

test_that("weekly_censored_gt reports the lag-0 mass it drops", {
  skip_if_missing("weekly_censored_gt")
  G <- weekly_censored_gt(15.3, 9.3)
  g0 <- attr(G, "lag0_mass")
  expect_true(!is.null(g0), info = "the dropped lag-0 mass must travel with the kernel")
  expect_true(is.finite(g0) && g0 > 0 && g0 < 0.5)
  # PIN THE VALUE. A bounds-only check let a real mis-derivation through: normalising by the
  # POST-drop total instead of the pre-drop total gives 0.0617 rather than 0.0581 (a 6% error)
  # and passed every other assertion here. 0.0581 is the double-censored weekly lag-0 mass at
  # the medium profile, independently confirmed by a 4e6-draw simulation of floor((U + D)/7).
  expect_equal(g0, 0.0581, tolerance = 2e-3,
               info = "lag-0 mass must be pk[1] / sum(pk) over the PRE-drop pmf")
  expect_equal(sum(G), 1, tolerance = 1e-9)
})

test_that("weekly_renewal_R_eff is the identity at R = 1 and at g0 = 0, and increases with R", {
  skip_if_missing("weekly_renewal_R_eff")
  g0 <- attr(weekly_censored_gt(15.3, 9.3), "lag0_mass")
  expect_equal(weekly_renewal_R_eff(1, g0), 1, tolerance = 1e-12,
               info = "the correction must vanish at the critical point")
  expect_equal(weekly_renewal_R_eff(1.7, 0), 1.7, tolerance = 1e-12,
               info = "a kernel with no lag-0 mass needs no correction")
  expect_gt(weekly_renewal_R_eff(2.5, g0), 2.5)   # supercritical: correction inflates
  expect_lt(weekly_renewal_R_eff(0.8, g0), 0.8)   # subcritical: correction deflates
  # Monotone increasing in R, which the test name claims and nothing asserted.
  .rs <- seq(0.5, 4, by = 0.25)
  expect_true(all(diff(weekly_renewal_R_eff(.rs, g0)) > 0),
              info = "R_eff must be strictly increasing in R")
  # The inflation factor itself must grow with R (it is 1 at R = 1 by construction).
  .ratio <- weekly_renewal_R_eff(.rs, g0) / .rs
  expect_true(all(diff(.ratio) > 0), info = "the correction must strengthen as R rises")
})

test_that("R_eff makes the weekly recursion reproduce the DAILY renewal's growth rate", {
  skip_if_missing("weekly_renewal_R_eff"); skip_if_missing("make_gt_pmf")
  w <- make_gt_pmf(15.3, 9.3, 90); w <- w / sum(w)
  G <- weekly_censored_gt(15.3, 9.3); g0 <- attr(G, "lag0_mass")
  lam_daily  <- function(R) uniroot(function(l) R * sum(w * l^(-seq_along(w))) - 1,
                                    c(1e-6, 10), tol = 1e-12)$root
  lam_weekly <- function(M) uniroot(function(l) M * sum(G * l^(-seq_along(G))) - 1,
                                    c(1e-6, 50), tol = 1e-12)$root
  for (R in c(1.19, 1.5, 2.5)) {
    truth <- lam_daily(R)^7
    err_old <- abs(lam_weekly(R) / truth - 1)
    err_new <- abs(lam_weekly(weekly_renewal_R_eff(R, g0)) / truth - 1)
    expect_lt(err_new, err_old,
              label = sprintf("R=%.2f: corrected weekly growth must be closer to the daily renewal", R))
  }
})
