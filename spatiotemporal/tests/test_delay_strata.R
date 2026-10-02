# =============================================================================
# tests/test_delay_strata.R — the onset-to-sample delay, fitted per case class
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# The onset imputation imputes onsets for CONFIRMED records, so it must draw from
# the confirmed-case delay. Until 2026-09-22 it drew from a fit pooled over every
# classification, and on this line list that pool is majority TEST-NEGATIVE:
# 6,823 `not_a_case` windowed pairs and 2,038 with NO final classification,
# against 4,981 confirmed and 145 suspected. Both groups are swabbed faster than
# confirmed cases (raw means 6.40 d and 3.56 d against 8.80 d), so the pooled EpiDist
# marginal returns 7.67 d where the CONFIRMED stratum gives 10.14 d — a 2.47 d
# (32%) gap, and imputed onsets landed that much late for the ~23% of confirmed
# records that carry one. (The interval-censored MLE puts the same contrast at
# 7.67 vs 9.11 d; quote the estimator with the number.) Because the invasion
# outcome is the FIRST onset in a zone, that shift can move a zone's invasion week.
#
# These tests pin the read/write contract, because the failure mode is silent: a
# stratum that does not resolve simply gives back the pooled delay.
# =============================================================================

.skip_unless <- function(fn) if (!exists(fn, mode = "function")) testthat::skip(paste0(fn, "() unavailable"))

# The numbers below are ILLUSTRATIVE, not the shipped fit: a fixture that happened to match
# the live parameter file could pass by reading the wrong file. They are of the right shape
# and ordering (confirmed slowest, suspected fastest), which is all these tests need.
.ds_write <- function() {
  tmp <- file.path(tempdir(), paste0("delay_params_", as.integer(runif(1, 1, 1e8)), ".csv"))
  fit <- list(params = c(shape = 0.81, rate = 0.106), best_family = "gamma",
              window = "w", n_fit = 100L, implied_mean = 7.665, rate = 0.1304)
  epi <- list(family = "gamma", mean = 7.665, sd = 8.507, mean_lo = 7.2, mean_hi = 8.1, n = 12043)
  str <- list(
    confirmed  = list(family="gamma", mean=9.12, sd=9.30, mean_lo=8.6, mean_hi=9.7, n=5006),
    not_a_case = list(family="gamma", mean=6.34, sd=6.80, mean_lo=6.0, mean_hi=6.7, n=6894),
    suspected  = list(family="gamma", mean=4.43, sd=5.10, mean_lo=3.4, mean_hi=5.6, n=147))
  invisible(suppressWarnings(write_onset_sample_long(fit, tmp, source_label = "test",
                                                     epidist = epi, strata = str)))
  tmp
}

test_that("stratum rows are ADDITIVE: the pooled read is unchanged by their presence", {
  .skip_unless("write_onset_sample_long"); .skip_unless(".load_dhis2_delay_params")
  tmp <- .ds_write(); on.exit(unlink(tmp), add = TRUE)
  p <- .load_dhis2_delay_params(tmp)
  # Every consumer that does not know about strata must resolve exactly what it did before.
  expect_equal(p$mean, 7.665, tolerance = 1e-9)
  expect_identical(p$family, "gamma")
  expect_identical(p$estimator, "epidist_marginal")
})

test_that("each stratum resolves its own fit, not the pooled one", {
  .skip_unless("write_onset_sample_long"); .skip_unless(".load_dhis2_delay_params")
  tmp <- .ds_write(); on.exit(unlink(tmp), add = TRUE)
  want <- c(confirmed = 9.12, not_a_case = 6.34, suspected = 4.43)
  for (s in names(want)) {
    p <- .load_dhis2_delay_params(tmp, stratum = s)
    expect_equal(p$mean, unname(want[[s]]), tolerance = 1e-9,
                 info = paste("stratum", s, "did not resolve its own mean"))
    expect_identical(p$estimator, paste0("epidist_marginal__", s))
    expect_identical(p$stratum, s)
  }
  # The confirmed delay must be materially LONGER than the pooled one — that difference is
  # the entire reason this exists. If it ever inverts, the fit or the strata are wrong.
  expect_gt(.load_dhis2_delay_params(tmp, stratum = "confirmed")$mean,
            .load_dhis2_delay_params(tmp)$mean)
})

test_that("a missing stratum WARNS rather than silently returning the pooled delay", {
  .skip_unless("write_onset_sample_long"); .skip_unless(".load_dhis2_delay_params")
  tmp <- .ds_write(); on.exit(unlink(tmp), add = TRUE)
  # Silence here would reinstate the majority-test-negative pooled fit under a name that
  # claims to be confirmed-only — the exact bug this change removes.
  expect_warning(.load_dhis2_delay_params(tmp, stratum = "nosuch"), "stratum")
  p <- suppressWarnings(.load_dhis2_delay_params(tmp, stratum = "nosuch"))
  expect_equal(p$mean, 7.665, tolerance = 1e-9)   # falls back, but only after warning
})

test_that("the resolver cache keys on the stratum", {
  .skip_unless("effective_onset_sample_delay"); .skip_unless("write_onset_sample_long")
  tmp <- .ds_write(); on.exit(unlink(tmp), add = TRUE)
  # The memo is keyed on the file's identity; without the stratum in the key the FIRST
  # stratum read would be served to every later caller on the same file.
  a <- suppressWarnings(effective_onset_sample_delay(tmp, stratum = "confirmed"))
  b <- suppressWarnings(effective_onset_sample_delay(tmp, stratum = "not_a_case"))
  d <- suppressWarnings(effective_onset_sample_delay(tmp))
  expect_equal(a$mean, 9.12,  tolerance = 1e-9)
  expect_equal(b$mean, 6.34,  tolerance = 1e-9)
  expect_equal(d$mean, 7.665, tolerance = 1e-9)
  # ...and reading the first one again must still give the first one.
  expect_equal(suppressWarnings(effective_onset_sample_delay(tmp, stratum = "confirmed"))$mean,
               9.12, tolerance = 1e-9)
})

test_that("the onset imputation asks for the confirmed stratum", {
  # Pinned on the source text: the call is inside load_linelist() and is not reachable
  # without a line list, but getting it wrong is silent and moves every invasion date.
  # testthat runs from tests/, so resolve the module relative to this file as well as the
  # working directory — whichever exists.
  cand <- c("01_data_prep.R", file.path("..", "01_data_prep.R"))
  f <- cand[file.exists(cand)][1]
  if (is.na(f)) testthat::skip("01_data_prep.R not reachable from the test working directory")
  src <- readLines(f, warn = FALSE)
  hit <- grep("effective_onset_sample_delay\\(stratum = \"confirmed\"\\)", src)
  expect_gte(length(hit), 1L)
})

# -----------------------------------------------------------------------------
# The panel that makes the contrast visible (04d panel E)
# -----------------------------------------------------------------------------
# Fitting the strata and never looking at them is how the pooled-delay bias
# survived in the first place. Panel E is the only view of the contrast, and it
# has two silent failure modes: the retained-figure gate can drop it, and a
# stratum whose rows are absent resolves to the POOLED fit (with a warning) and
# would be drawn three times under three names.

test_that("the by-classification panel is on the retained-figure allow-list", {
  .skip_unless("figure_is_kept")
  # .dfg_save() gates on the full path and tolerates the _YYYYMMDD stamp.
  expect_true(figure_is_kept(
    file.path("outputs", "diagnostics", "delay_fits", "dhis2_onset_sample_by_class_20260921")))
})

test_that("a stratum that falls back to POOLED is excluded from the panel, not drawn", {
  .skip_unless("write_onset_sample_long"); .skip_unless(".load_dhis2_delay_params")
  tmp <- .ds_write()
  # `probable` was never fitted: the loader warns and hands back the pooled fit with
  # $stratum = NA. That is exactly the row the panel must drop -- drawing it would put the
  # pooled curve on the figure under a stratum's name and invent a contrast that is not there.
  got <- suppressWarnings(.load_dhis2_delay_params(tmp, stratum = "probable"))
  expect_true(is.na(got$stratum))
  expect_equal(got$mean, 7.665, tolerance = 1e-6)          # the pooled mean, not a stratum's
  keep <- !is.null(got) && isTRUE(!is.na(got$stratum)) && is.finite(got$mean) && got$mean > 0
  expect_false(keep)                                        # the panel's own filter
  # A stratum that IS present survives the same filter.
  ok <- .load_dhis2_delay_params(tmp, stratum = "confirmed")
  expect_identical(ok$stratum, "confirmed")
  expect_true(!is.null(ok) && isTRUE(!is.na(ok$stratum)) && is.finite(ok$mean) && ok$mean > 0)
})

test_that("the panel's completeness curve is a CDF of the delay the pipeline resolves", {
  .skip_unless("write_onset_sample_long"); .skip_unless(".load_dhis2_delay_params")
  tmp <- .ds_write()
  sp <- .load_dhis2_delay_params(tmp, stratum = "confirmed")
  # The panel reparameterises the published (mean, sd) exactly as the loader does; it does
  # not fit. Pin that the params it plots ARE the loader's, and that the curve is a CDF.
  expect_identical(sp$family, "gamma")
  expect_equal(unname(sp$params[["shape"]]), (sp$mean / sp$sd)^2, tolerance = 1e-8)
  expect_equal(unname(sp$params[["rate"]]),   sp$mean / sp$sd^2,  tolerance = 1e-8)
  xs <- seq(0, 40, length.out = 200)
  cdf <- stats::pgamma(xs, shape = sp$params[["shape"]], rate = sp$params[["rate"]])
  expect_true(all(diff(cdf) >= -1e-12))          # monotone
  expect_true(all(cdf >= 0 & cdf <= 1))
  expect_equal(cdf[1], 0, tolerance = 1e-12)
  # Shape < 1 is why the panel plots a CDF and not a density: the density is unbounded at 0.
  expect_lt(sp$params[["shape"]], 1)
})
