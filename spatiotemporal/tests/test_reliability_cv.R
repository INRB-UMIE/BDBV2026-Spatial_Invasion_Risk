# =============================================================================
# tests/test_reliability_cv.R — the reliability diagram's uncertainty is CROSS-VALIDATED
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Review 2.25 asked for uncertainty derived from cross-validation rather than asymptotic
# binomial intervals: a calibration diagram is about OUT-OF-SAMPLE forecasts, and the
# zone-weeks in a bin are not independent draws — they arrive in LFO folds that share a
# refit, an epidemic state and a week of reporting. reliability_curve() therefore reports a
# between-fold interval (a Wilson interval on N / design effect) beside the Wilson one, and
# reliability_fold_cells() writes the per-origin cells the figure draws behind the points.
#
# Pinned here because each of these can regress silently — the figure would still render:
#   * the `fold` argument being dropped at the call site, which reverts FigR1 to binomial
#     intervals with nothing in the panel saying so;
#   * the design effect losing its not-estimable guards, which sends a zero-variance bin's
#     effective sample size to Inf and prints a zero-width interval on a bin holding one
#     zone-week;
#   * the fold cells drifting out of agreement with the curve they decompose.
# =============================================================================

.rel_fns <- function() all(vapply(c("reliability_curve", "reliability_fold_cells"),
                                  exists, logical(1), mode = "function"))

test_that("without folds the curve keeps its original contract", {
  if (!.rel_fns()) testthat::skip("23_prospective_eval.R not sourced")
  set.seed(11)
  p <- runif(400); o <- rbinom(400, 1, p)
  r <- reliability_curve(p, o)
  expect_identical(names(r), c("bin", "bin_lo", "bin_hi", "p_mid", "mean_pred", "n", "k",
                               "obs_freq", "obs_lo", "obs_hi"))
  # The edges travel WITH the curve. With expected-balanced bins (the default) a midpoint no
  # longer implies them, so the curve must carry the very edges reliability_breaks() returns.
  b <- reliability_breaks(p, o)
  expect_identical(r$bin_lo, b[r$bin]); expect_identical(r$bin_hi, b[r$bin + 1L])
  # Equal-width bins are still available on request, and then "bin b spans 10*(b-1)% to 10*b%".
  w <- reliability_curve(p, o, n_bins = 10L, binning = "width")
  expect_identical(w$bin_lo, seq(0, 0.9, by = 0.1)[w$bin])
})

test_that("the defaults are the suite's convention, so a new caller cannot drift from it", {
  if (!.rel_fns()) testthat::skip("23_prospective_eval.R not sourced")
  # The defaults used to be 95% intervals and equal-width deciles while every published table
  # and figure passes RELIABILITY_CONF (90%) and expected-balanced RELIABILITY_N_BINS bins, so a
  # caller that relied on the defaults silently produced a differently binned curve at another
  # interval level.
  set.seed(4)
  p <- c(runif(2000, 0, 0.02), runif(150, 0.02, 1)); o <- rbinom(length(p), 1, p)
  f <- rep(1:12, length.out = length(p))
  expect_identical(reliability_breaks(p, o),
                   reliability_breaks(p, o, RELIABILITY_N_BINS, RELIABILITY_BINNING))
  expect_identical(reliability_curve(p, o, fold = f),
                   reliability_curve(p, o, fold = f, conf = RELIABILITY_CONF,
                                     breaks = reliability_breaks(p, o, RELIABILITY_N_BINS,
                                                                 RELIABILITY_BINNING)))
  expect_identical(reliability_fold_cells(p, o, f),
                   reliability_fold_cells(p, o, f, breaks = reliability_breaks(
                     p, o, RELIABILITY_N_BINS, RELIABILITY_BINNING)))
})

test_that("independent folds give a design effect near 1 and an interval near Wilson", {
  if (!.rel_fns()) testthat::skip("23_prospective_eval.R not sourced")
  # Outcomes drawn from the predicted probability, folds assigned AT RANDOM: there is no
  # between-fold structure, so the clustered variance must land on the binomial one.
  set.seed(23)
  p <- runif(6000); o <- rbinom(6000, 1, p); f <- sample(1:12, 6000, replace = TRUE)
  r <- reliability_curve(p, o, fold = f)
  expect_true(all(r$cv_basis == "between-fold"))
  expect_equal(median(r$deff), 1, tolerance = 0.25)
  expect_equal(r$obs_lo_cv, r$obs_lo, tolerance = 0.06)
  expect_equal(r$obs_hi_cv, r$obs_hi, tolerance = 0.06)
})

test_that("folds that disagree widen the interval beyond the binomial one", {
  if (!.rel_fns()) testthat::skip("23_prospective_eval.R not sourced")
  # One bin, ten folds, the SAME pooled frequency as the independent case but every event
  # concentrated in two folds — exactly the "the forecast worked in some rounds and not in
  # others" pattern a binomial interval cannot see.
  p <- rep(0.45, 1000)
  f <- rep(1:10, each = 100)
  o <- c(rep(1L, 100), rep(1L, 100), rep(0L, 800))   # folds 1-2 all events, rest none
  r <- reliability_curve(p, o, fold = f)
  expect_equal(nrow(r), 1L)
  expect_identical(r$cv_basis, "between-fold")
  expect_gt(r$deff, 10)
  expect_gt(r$obs_hi_cv - r$obs_lo_cv, 3 * (r$obs_hi - r$obs_lo))
  expect_equal(r$obs_freq, 0.2)
  expect_identical(r$n_folds, 10L)
  expect_identical(r$n_folds_with_events, 2L)
})

test_that("bins with no between-fold contrast fall back to the binomial interval and say so", {
  if (!.rel_fns()) testthat::skip("23_prospective_eval.R not sourced")
  # Bin 1: two folds, no events at all. Bin 10: a single fold. Neither can support a design
  # effect; both must report the Wilson interval unchanged rather than a zero-width one.
  p <- c(rep(0.05, 20), rep(0.95, 3))
  o <- c(rep(0L, 20), 1L, 0L, 1L)
  f <- c(rep(c("a", "b"), each = 10), rep("a", 3))
  r <- reliability_curve(p, o, fold = f, n_bins = 10L, binning = "width")
  z <- r[r$bin == 1L, ]; s <- r[r$bin == 10L, ]
  expect_identical(z$cv_basis, "binomial (no outcome contrast)")
  expect_identical(s$cv_basis, "binomial (single fold)")
  expect_true(is.na(z$deff) && is.na(s$deff))
  expect_equal(c(z$obs_lo_cv, z$obs_hi_cv), c(z$obs_lo, z$obs_hi))
  expect_equal(c(s$obs_lo_cv, s$obs_hi_cv), c(s$obs_lo, s$obs_hi))
  expect_gt(z$obs_hi_cv, 0)   # a zero-event bin still has an upper bound
})

test_that("the fold cells decompose exactly the curve they are drawn behind", {
  if (!.rel_fns()) testthat::skip("23_prospective_eval.R not sourced")
  set.seed(7)
  p <- runif(2000); o <- rbinom(2000, 1, p); f <- rep(1:13, length.out = 2000)
  r <- reliability_curve(p, o, fold = f)
  cells <- reliability_fold_cells(p, o, f)
  agg <- stats::aggregate(cbind(n, k) ~ bin, data = cells, FUN = sum)
  agg <- agg[match(r$bin, agg$bin), ]
  expect_equal(agg$n, r$n)
  expect_equal(agg$k, r$k)
  expect_equal(vapply(split(cells$fold, cells$bin), function(x) length(unique(x)), integer(1),
                      USE.NAMES = FALSE), r$n_folds)
  # Fold 10 must not sort between 1 and 2: the table's whole purpose is the per-origin view.
  b1 <- cells$fold[cells$bin == cells$bin[1]]
  expect_identical(b1, as.character(sort(as.numeric(b1))))
})

test_that("expected-balanced edges split the expected events, not the rows", {
  if (!exists("reliability_breaks", mode = "function")) testthat::skip("23_prospective_eval.R not sourced")
  # The shape this binning exists for: a huge mass of near-zero predictions carrying a few
  # expected events, and a thin tail carrying the rest. Equal-width or equal-count edges collapse
  # one side or the other; expected-balanced edges give every bin a comparable sum(p).
  set.seed(31)
  p <- c(runif(5000, 0, 0.02), runif(200, 0.02, 1))
  o <- rbinom(length(p), 1, p)
  b <- reliability_breaks(p, o, n_bins = 5L, binning = "expected")
  r <- reliability_curve(p, o, breaks = b)
  expect_equal(nrow(r), 5L)
  e <- r$mean_pred * r$n                               # expected events per bin
  expect_lt(max(e) - min(e), 0.2 * mean(e))
  expect_identical(r$bin_lo[1], 0); expect_identical(r$bin_hi[nrow(r)], 1)
  expect_true(all(diff(r$bin_lo) > 0))
  # And the contrast with the conventional binning, on the same data.
  w <- reliability_curve(p, o, n_bins = 10L, binning = "width")
  expect_gt(max(w$n) / sum(w$n), 0.9)               # one bin swallows the record
  expect_lt(max(r$n) / sum(r$n), max(w$n) / sum(w$n))
})

test_that("the edges do not depend on the outcomes", {
  if (!exists("reliability_breaks", mode = "function")) testthat::skip("23_prospective_eval.R not sourced")
  # Outcome-dependent edges (the retired observed-event binning closed each bin at its k-th
  # observed invasion) bias every bin's observed frequency upward. The edges must be a function
  # of the forecasts alone.
  set.seed(8)
  p <- c(runif(3000, 0, 0.01), runif(150, 0.01, 0.8))
  b1 <- reliability_breaks(p, rbinom(length(p), 1, p), 5L, "expected")
  b2 <- reliability_breaks(p, rbinom(length(p), 1, p), 5L, "expected")
  expect_identical(b1, b2)
  expect_identical(b1, reliability_breaks(p, NULL, 5L, "expected"))
  expect_error(reliability_breaks(p, rep(0L, length(p)), 5L, "event"))
})

test_that("tied predictions are never split across an edge", {
  if (!exists("reliability_breaks", mode = "function")) testthat::skip("23_prospective_eval.R not sourced")
  # Every prediction identical: there is no edge that separates them, so one bin is the only
  # honest answer. Splitting on the sort order would invent a difference between two
  # identical forecasts.
  b <- reliability_breaks(rep(0.3, 100), rep(c(0L, 1L), 50), n_bins = 5L, binning = "expected")
  expect_identical(b, c(0, 1))
  # Half tied, half spread: the tied block stays whole.
  p <- c(rep(0.2, 200), seq(0.3, 0.9, length.out = 40))
  o <- c(rbinom(200, 1, 0.2), rbinom(40, 1, 0.5))
  r <- reliability_curve(p, o, breaks = reliability_breaks(p, o, 4L, "expected"))
  expect_equal(sum(r$n[r$bin_lo < 0.2 & r$bin_hi > 0.2]), 200)
})

test_that("a record with no events still gets populated bins", {
  if (!exists("reliability_breaks", mode = "function")) testthat::skip("23_prospective_eval.R not sourced")
  p <- runif(500, 0, 0.3)
  b <- reliability_breaks(p, rep(0L, 500), n_bins = 5L, binning = "expected")
  expect_length(b, 6L)
  expect_true(all(diff(b) > 0))
  r <- reliability_curve(p, rep(0L, 500), breaks = b)
  expect_equal(nrow(r), 5L)
  expect_true(all(r$k == 0))
  expect_true(all(r$cv_basis == "binomial (no outcome contrast)" | is.na(r$deff)))
  # all-zero forecasts: nothing to balance, equal-count fallback, still valid edges
  expect_identical(reliability_breaks(rep(0, 50), NULL, 5L, "expected"), c(0, 1))
})

test_that("the curve and its per-origin cells are binned on the same edges", {
  if (!.rel_fns()) testthat::skip("23_prospective_eval.R not sourced")
  set.seed(5)
  p <- c(runif(2000, 0, 0.02), runif(120, 0.02, 1)); o <- rbinom(length(p), 1, p)
  f <- rep(1:13, length.out = length(p))
  b <- reliability_breaks(p, o, 5L, "expected")
  r <- reliability_curve(p, o, fold = f, breaks = b)
  cells <- reliability_fold_cells(p, o, f, breaks = b)
  agg <- stats::aggregate(cbind(n, k) ~ bin, data = cells, FUN = sum)
  expect_equal(agg$n[match(r$bin, agg$bin)], r$n)
  expect_equal(agg$k[match(r$bin, agg$bin)], r$k)
})

test_that("the binning is a named config constant, not an inline literal", {
  expect_true(exists("RELIABILITY_N_BINS") && exists("RELIABILITY_BINNING"))
  expect_true(is.numeric(RELIABILITY_N_BINS) && RELIABILITY_N_BINS >= 3L)
  expect_identical(RELIABILITY_BINNING, "expected")
  # One interval level for the whole calibration figure. Panel C used the estimator's 95%
  # default while the per-round panels carried 90% Poisson intervals, so a single figure
  # reported two levels.
  expect_true(is.numeric(RELIABILITY_CONF) && RELIABILITY_CONF > 0 && RELIABILITY_CONF < 1)
  # Both ends must READ it. The published table and the figure panel describe the same curve
  # only while they bin identically, and a literal at either end is how that drifts apart.
  for (f in c("run_all.R", "make_si_model_figures.R")) {
    src <- readLines(file.path(here::here(), "spatiotemporal", f), warn = FALSE)
    expect_true(any(grepl("RELIABILITY_N_BINS", src, fixed = TRUE)),
                info = paste(f, "must read the named bin count"))
    expect_true(any(grepl("RELIABILITY_BINNING", src, fixed = TRUE)),
                info = paste(f, "must read the named binning"))
    expect_true(any(grepl("conf = RELIABILITY_CONF", src, fixed = TRUE)),
                info = paste(f, "must read the named interval level"))
  }
})

test_that("the pipeline clusters on the fold and the calibration figure draws the result", {
  ra <- file.path(here::here(), "spatiotemporal", "run_all.R")
  rf <- file.path(here::here(), "spatiotemporal", "make_si_model_figures.R")
  if (!file.exists(ra) || !file.exists(rf)) testthat::skip("pipeline sources unavailable")
  # Matched on the CALL, not on the whole file: the argument must sit inside the
  # reliability_curve() call (which wraps onto a second line), and a whole-file match would
  # also print the whole file on failure.
  src  <- readLines(ra, warn = FALSE)
  i    <- grep("rc  <- reliability_curve\\(", src)
  expect_length(i, 1L)
  call <- paste(src[i:(i + 1L)], collapse = " ")
  expect_match(call, "fold\\s*=\\s*fvo\\$fold_id",
               info = "run_all.R must cluster the reliability curve on the LFO fold")
  expect_match(call, "breaks\\s*=\\s*\\.brk",
               info = "the curve must use the shared event-balanced edges")
  expect_true(any(grepl("reliability_fold_cells\\(", src)),
              info = "run_all.R must write the per-origin cells beside the curve")
  expect_true(any(grepl("forecast_reliability_folds_h%d.csv", src, fixed = TRUE)))
  # FigR1 is GONE: its content is panel C of FigureS7_calibration. A stray call would publish
  # a second, differently binned reliability diagram of the same model.
  rev <- readLines(file.path(here::here(), "spatiotemporal", "24_review_figures.R"), warn = FALSE)
  expect_false(any(grepl("plot_reliability_review|FigR1_reliability_h", rev)),
               info = "the review suite must no longer draw a reliability diagram")
  fig <- readLines(rf, warn = FALSE)
  expect_true(any(grepl("obs_lo_cv", fig, fixed = TRUE)),
              info = "the calibration figure must draw the cross-validated interval")
  expect_true(any(grepl("23_prospective_eval.R", fig, fixed = TRUE)),
              info = "the figure must use the pipeline's own estimator, not its own copy")
})
