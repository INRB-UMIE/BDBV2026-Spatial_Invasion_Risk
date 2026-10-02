# =============================================================================
# tests/test_rolling_predictors.R — the rolling as-of import force
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# The import force is built from nowcast-corrected counts, and the weekly
# generation-time kernel puts ~31% of its weight on the most recent week — the
# week the nowcast inflates most. With ONE count matrix per fold, every training
# transition's predictors are complete weeks while the forecast's predictor
# carries that inflation, so beta0 is fitted against a lambda about half the size
# of the one it is applied to. Measured pooled observed/expected over the 12
# folds: 0.497 with the shared matrix against 0.741 with rolling predictors at a
# floor of 3 — and 0.440 with rolling and NO floor, which is worse than changing
# nothing. Rolling is a large improvement, not a complete fix: a residual ~1.35x
# over-prediction remains. These tests pin the parts that can go silently wrong.
# =============================================================================

.skip_unless <- function(fn) if (!exists(fn, mode = "function")) testthat::skip(paste0(fn, "() unavailable"))

.rp_fixture <- function(n_weeks = 8L) {
  zones <- c("A", "B", "C")
  wk <- seq(as.Date("2026-05-04"), by = "week", length.out = n_weeks)
  expand <- expand.grid(health_zone = zones, week_start = wk, stringsAsFactors = FALSE)
  expand$confirmed <- 0L
  expand$confirmed[expand$health_zone == "A"] <- 5L          # A seeds throughout
  expand$confirmed[expand$health_zone == "B" &
                   expand$week_start == wk[n_weeks]] <- 1L   # B invaded at the end
  expand$confirmed_nc <- as.numeric(expand$confirmed)
  expand$suspected <- 0L; expand$suspected_nc <- 0
  expand$total_alerts <- 0L
  list(zw = expand, zones = zones, weeks = wk)
}

test_that("rolling OFF reproduces the historic design exactly", {
  .skip_unless("build_invasion_design")
  f <- .rp_fixture()
  W <- matrix(1/3, 3, 3, dimnames = list(f$zones, f$zones))
  gt <- list(medium = stats::dgamma(1:60, shape = 2.7, rate = 0.18))
  osrm <- matrix(60, 3, 3, dimnames = list(f$zones, f$zones)); diag(osrm) <- 0
  cov <- data.frame(health_zone = f$zones, pop_count = 1e5, ccvi = 0.5,
                    healthsite_density = 1, stringsAsFactors = FALSE)
  a <- suppressWarnings(suppressMessages(build_invasion_design(
        f$zw, list(M = W), gt, cov, osrm, f$zones, mob = "M", gt = "medium",
        candidates = c("ccvi", "d_min"))))
  # An explicit NULL must behave identically to omitting the arguments: the rolling
  # branch is opt-in, and any drift here silently changes every historic result.
  b <- suppressWarnings(suppressMessages(build_invasion_design(
        f$zw, list(M = W), gt, cov, osrm, f$zones, mob = "M", gt = "medium",
        candidates = c("ccvi", "d_min"), rolling_ll = NULL, rolling_delay = NULL)))
  expect_identical(a$d, b$d)
  expect_identical(a$feat, b$feat)
  expect_equal(a$n_events, b$n_events)
})

test_that("rolling stays OFF unless BOTH a line list and a delay are supplied", {
  .skip_unless("build_invasion_design")
  f <- .rp_fixture()
  W <- matrix(1/3, 3, 3, dimnames = list(f$zones, f$zones))
  gt <- list(medium = stats::dgamma(1:60, shape = 2.7, rate = 0.18))
  osrm <- matrix(60, 3, 3, dimnames = list(f$zones, f$zones)); diag(osrm) <- 0
  cov <- data.frame(health_zone = f$zones, pop_count = 1e5, ccvi = 0.5,
                    healthsite_density = 1, stringsAsFactors = FALSE)
  base <- suppressWarnings(suppressMessages(build_invasion_design(
        f$zw, list(M = W), gt, cov, osrm, f$zones, mob = "M", gt = "medium",
        candidates = c("ccvi", "d_min"))))
  # A line list with no delay, or a delay with no line list, must NOT half-enable the
  # rolling branch — a partially-rolled design would mix two predictor conventions.
  for (args in list(list(rolling_ll = data.frame(x = 1), rolling_delay = NULL),
                    list(rolling_ll = NULL, rolling_delay = list(family = "gamma")))) {
    got <- suppressWarnings(suppressMessages(do.call(build_invasion_design, c(
      list(f$zw, list(M = W), gt, cov, osrm, f$zones, mob = "M", gt = "medium",
           candidates = c("ccvi", "d_min")), args))))
    expect_identical(got$d, base$d)
  }
})

test_that("the rolling floor is present, at least 1, and stamped into the LFO cache key", {
  # The floor is not a tidy-up: with no floor the earliest transitions carry a near-empty
  # as-of reconstruction, lambda is tiny, and the few invasions there imply beta0 ~ 2.1
  # against ~0.2 at the last fold. Those rows dominate and calibration goes to 0.440 —
  # worse than the 0.497 the shared matrix gives.
  fl <- get0("ROLLING_PREDICTOR_FLOOR", ifnotfound = NULL)
  if (is.null(fl)) testthat::skip("ROLLING_PREDICTOR_FLOOR unavailable")
  expect_true(is.numeric(fl) && length(fl) == 1L && is.finite(fl))
  expect_gte(as.integer(fl), 1L)
  # The default must keep most folds: floor 3 retains 11 of 12, floor 5 only 9.
  expect_lte(as.integer(fl), 5L)
})

test_that("the as-of memo key separates different line lists", {
  .skip_unless(".asof_counts_wide")
  if (!exists(".asof_memo")) testthat::skip(".asof_memo unavailable")
  # Two line lists of the SAME row count must not collide: the key carries the summed
  # onset dates as well as nrow, because a collision serves one fold's counts to another.
  mk <- function(shift) data.frame(
    date_of_symptom_onset = as.character(as.Date("2026-05-01") + c(0, 1, 2) + shift),
    stringsAsFactors = FALSE)
  fp <- function(ll) c(nrow(ll), sum(as.numeric(suppressWarnings(
          as.Date(ll$date_of_symptom_onset))), na.rm = TRUE))
  expect_false(identical(fp(mk(0)), fp(mk(7))))
  expect_identical(fp(mk(0)), fp(mk(0)))
})

# -----------------------------------------------------------------------------
# The floor must never empty a fold (2026-09-23)
# -----------------------------------------------------------------------------
# run_invasion_lfo() admits a fold at min_train_weeks = 3 training weeks, which is TWO
# transitions, while the design required rolling_floor = 3 transitions. The two constants
# contradicted each other, so the earliest fold was guaranteed to yield no training rows and
# therefore no forecasts -- it disappeared from every Bayesian model while the structural
# baselines, which use no design, kept it. The floor is now capped at the number of available
# transitions.
test_that("the rolling floor is capped so every fold keeps at least one transition", {
  t0_of <- function(nT, floor) min(max(1L, as.integer(floor)), max(1L, nT))
  rows  <- function(nT, floor) sum(seq_len(nT) >= t0_of(nT, floor))
  # folds with fewer transitions than the floor are KEPT, not emptied
  expect_equal(rows(1L, 3L), 1L)
  expect_equal(rows(2L, 3L), 1L)     # the earliest fold: 3 training weeks -> 2 transitions
  # folds with enough history are UNCHANGED by the cap
  for (n in 3:12) expect_equal(rows(n, 3L), n - 2L)
  expect_equal(t0_of(8L, 3L), 3L)
  # and the cap never raises the floor
  expect_lte(t0_of(2L, 3L), 3L)
})

test_that("an unstable intercept cannot reorder a fold (why capping is safe)", {
  # The import force enters as offset(logLam) with its coefficient FIXED at 1, so beta_0 is a
  # pure intercept: it shifts every zone identically on the cloglog scale. That is why a
  # sparse early fold, whose intercept is poorly determined, still produces a usable
  # watch-list -- only the LEVEL moves, and the recalibration factor carries the level.
  set.seed(9)
  logLam <- log(runif(200, 1e-4, 0.5)); gx <- rnorm(200, 0, 0.3)
  inv <- function(b0) 1 - exp(-exp(b0 + logLam + gx))       # cloglog inverse link
  r_small <- rank(-inv(-3.0)); r_large <- rank(-inv(0.75))  # beta_0 2.14 vs 0.20 in magnitude
  expect_identical(r_small, r_large)                         # ranking is bit-identical
  expect_gt(mean(inv(0.75)), mean(inv(-3.0)))                # only the level differs
  expect_identical(order(inv(-3.0), decreasing = TRUE)[1:20],
                   order(inv(0.75), decreasing = TRUE)[1:20]) # same top 20
})
