# =============================================================================
# tests/test_dhis2_delay.R — DHIS2 onset->sample delay estimation
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Exercises 04c_dhis2_delay_windows.R:
#   - .make_cens_df()  interval-censoring convention (d=0->[0,0.5]; d>0->[d-.5,d+.5])
#   - .fit_all_censored() recovers a known gamma from synthetic delays
#   - estimate_dhis2_onset_sample_delay() windows + fits + summarises
#   - write_onset_sample_long() emits the pipeline long format WITHOUT a `rate`
#     name collision (Exp-rate summary vs the family's native `rate`), keeps or drops an
#     earlier EpiDist fit by the same-records rule, writes selected_param_* of the selected
#     family, and persists EpiDist convergence diagnostics
#   - build_dhis2_delay_populations() shares the estimator's window; 04d's ECDF gap
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
  # the interval-censored MLE — the bias this whole change removed — and every other assertion
  # here still passed.
  expect_identical(dp$estimator, "epidist_marginal",
                   info = "the resolver must be the truncation-corrected EpiDist marginal")
  # Relative to the SAME file's windowed MLE, not a fixed number of days: a fixed threshold
  # goes stale with every snapshot (and with the sample-date rule). Correcting the right-
  # truncation the windowed MLE still carries moves the mean UP.
  tab <- utils::read.csv(DELAY_PARAMS_PATH, colClasses = "character")
  mle <- suppressWarnings(as.numeric(tab$value[tab$quantity == "implied_mean_fit"][1]))
  skip_if(!is.finite(mle), "params file carries no MLE row")
  expect_gt(dp$mean, mle)
  expect_true(is.finite(dp$mean) && dp$mean > 0)
  set.seed(1)
  d <- .draw_dhis2_delay(dp, 2e5)
  expect_equal(mean(d), dp$mean, tolerance = 0.05,
               info = "draws must reproduce the resolver's mean, not a truncated one")
  expect_equal(stats::sd(d), dp$sd, tolerance = 0.10,
               info = "draws must reproduce the resolver's SD")
})

test_that("load_linelist() imputes missing onsets from same-class, same-week donor intervals", {
  # Since 2026-10-15 a missing onset is the sample date minus an interval drawn from the
  # recorded (sample - onset) intervals of CONFIRMED records sampled in the SAME week
  # (.impute_onset_delays(), 01_data_prep.R). On the real line list, in weeks whose own donors
  # suffice (no borrowing) and late enough that no draw can hit the outbreak floor, the mean
  # imputed gap must equal the recipient-weighted donor mean up to sampling error. Aggregates
  # only: nothing record-level is printed.
  skip_if_missing("load_linelist")
  skip_if(!identical(get0("ONSET_MODE", ifnotfound = "impute"), "impute"), "complete-case mode")
  ll <- tryCatch(suppressWarnings(suppressMessages(load_linelist())), error = function(e) NULL)
  skip_if(is.null(ll) || !any(ll$onset_imputed %in% TRUE), "line list unavailable")
  s    <- ll$date_of_sample_collection
  wk   <- as.character(lubridate::floor_date(s, "week", week_start = get0("WEEK_ANCHOR", ifnotfound = 1L)))
  conf <- ll$confirmed %in% TRUE
  don  <- conf & ll$onset_usable %in% TRUE & !(ll$onset_imputed %in% TRUE) & !is.na(s)
  rec  <- conf & ll$onset_imputed %in% TRUE & !is.na(s) & !is.na(ll$date_index)
  ivl  <- pmax(as.integer(s[don] - ll$date_index[don]), 0L)
  gap  <- as.integer(s[rec] - ll$date_index[rec])
  expect_true(all(gap >= 0L))
  expect_true(all(gap <= max(ivl)))
  dw <- split(ivl, wk[don]); rw <- split(gap, wk[rec])
  ok <- intersect(names(dw)[lengths(dw) >= get0("ONSET_DONOR_MIN", ifnotfound = 20L)], names(rw))
  ok <- ok[as.Date(ok) >= OUTBREAK_START + get0("DELAY_MAX_PLAUSIBLE_DAYS", ifnotfound = 60L)]
  skip_if(length(ok) < 3L, "too few donor-rich weeks")
  n  <- sum(lengths(rw[ok]))
  E  <- sum(vapply(ok, function(w) length(rw[[w]]) * mean(dw[[w]]), 0)) / n
  V  <- sum(vapply(ok, function(w) length(rw[[w]]) * stats::var(dw[[w]]), 0)) / n^2
  expect_lt(abs(mean(unlist(rw[ok])) - E), 4 * sqrt(V))
})

test_that("the CANONICAL delay file is refused when it names another snapshot; an explicit file is not", {
  skip_if_missing(".load_dhis2_delay_params")
  cur <- tryCatch(jsonlite::fromJSON(LINELIST_JSON)$folder, error = function(e) NA_character_)
  skip_if(is.na(cur), "no current snapshot pointer")
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(quantity = c("source", "epidist_family", "epidist_mean", "epidist_sd"),
                              value = c(paste0(cur, "_OTHER"), "gamma", "9", "7")),
                   tmp, row.names = FALSE)
  old <- get("DELAY_PARAMS_PATH", envir = globalenv())
  on.exit(assign("DELAY_PARAMS_PATH", old, envir = globalenv()), add = TRUE)
  assign("DELAY_PARAMS_PATH", tmp, envir = globalenv())        # tmp is now the canonical file
  expect_warning(r <- .load_dhis2_delay_params(tmp), "NOT used")
  expect_null(r)
  assign("DELAY_PARAMS_PATH", old, envir = globalenv())        # tmp is an explicit file again
  w <- character(0)
  withCallingHandlers(.load_dhis2_delay_params(tmp),
                      warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("NOT used", w)))
})

# =============================================================================
# write_onset_sample_long(): which EpiDist rows survive a write WITHOUT a fit, and what the
# `selected_*` rows then say. Round trips through the REAL loader, on temp files only.
# =============================================================================

.wl_fit <- function(window = "analytical_2026-10-03", n = 100L, mean = 6.667)
  list(params = c(shape = 0.8, rate = 0.8 / mean), best_family = "gamma", window = window,
       n_fit = n, implied_mean = mean, rate = 1 / mean)
.wl_epi <- function(mean = 7.5, sd = 8, family = "gamma", ...)
  c(list(family = family, mean = mean, sd = sd, mean_lo = mean - 0.3, mean_hi = mean + 0.3,
         n = 100L), list(...))
.wl_q <- function(out, q) { v <- out$value[out$quantity == q]; if (length(v)) v[1] else NA_character_ }
# Run `code` with the given globals swapped in, restoring them afterwards.
.wl_with <- function(vals, code) {
  old <- lapply(names(vals), function(n) get0(n, envir = globalenv(), ifnotfound = NULL))
  names(old) <- names(vals)
  on.exit(for (n in names(old)) assign(n, old[[n]], envir = globalenv()), add = TRUE)
  for (n in names(vals)) assign(n, vals[[n]], envir = globalenv())
  force(code)
}

test_that("another snapshot's EpiDist rows are DROPPED, not re-stamped as this snapshot's", {
  skip_if_missing("write_onset_sample_long"); skip_if_missing(".load_dhis2_delay_params")
  tmp <- tempfile(fileext = ".csv"); js <- tempfile(fileext = ".json")
  on.exit(unlink(c(tmp, js)), add = TRUE)
  writeLines('{"folder": "LINELIST_NEW"}', js)
  .wl_with(list(RUN_EPIDIST = TRUE, DELAY_PARAMS_PATH = tmp, LINELIST_JSON = js), {
    # Snapshot OLD: pooled EpiDist 7, confirmed stratum 10.
    write_onset_sample_long(.wl_fit("analytical_2026-09-26", 90L, 6.0), tmp, "LINELIST_OLD",
                            epidist = .wl_epi(7),
                            strata = list(confirmed = .wl_epi(10)))
    # Snapshot NEW, EpiDist failed (epidist = NULL): the OLD fit must not survive.
    expect_warning(out <- write_onset_sample_long(.wl_fit(), tmp, "LINELIST_NEW"), "DROPPED")
    expect_false(any(grepl("^epidist_", out$quantity)))
    expect_identical(.wl_q(out, "selected_estimator"), "interval_censored_mle")
    r <- .load_dhis2_delay_params(tmp)          # tmp is the canonical file; NEW is current
    expect_identical(r$estimator, "interval_censored_mle")
    expect_equal(r$mean, 6.667, tolerance = 1e-9)
    # No confirmed row survived: the stratum read warns and resolves the pooled MLE.
    expect_warning(rc <- .load_dhis2_delay_params(tmp, stratum = "confirmed"), "stratum")
    expect_identical(rc$estimator, "interval_censored_mle")
  })
})

test_that("a fit of the SAME records is preserved when EpiDist was requested but failed", {
  skip_if_missing("write_onset_sample_long"); skip_if_missing(".load_dhis2_delay_params")
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  .wl_with(list(RUN_EPIDIST = TRUE), {
    write_onset_sample_long(.wl_fit(), tmp, "LINELIST_NEW", epidist = .wl_epi(7.5),
                            strata = list(confirmed = .wl_epi(10)))
    expect_warning(out <- write_onset_sample_long(.wl_fit(), tmp, "LINELIST_NEW"), "PRESERVING")
    expect_identical(.wl_q(out, "selected_estimator"), "epidist_marginal_truncation_corrected")
    expect_equal(as.numeric(.wl_q(out, "selected_mean")), 7.5)
    expect_identical(anyDuplicated(out$quantity), 0L)
    expect_equal(.load_dhis2_delay_params(tmp)$mean, 7.5)
    expect_equal(.load_dhis2_delay_params(tmp, stratum = "confirmed")$mean, 10)
  })
  # Same window and n but a different MLE mean = different records (e.g. another sample-date
  # rule): dropped.
  .wl_with(list(RUN_EPIDIST = TRUE),
           expect_warning(out2 <- write_onset_sample_long(.wl_fit(mean = 5.9), tmp, "LINELIST_NEW"),
                          "DROPPED"))
  expect_false(any(grepl("^epidist_", out2$quantity)))
})

test_that("RUN_EPIDIST = FALSE falls back to the censored MLE, as 00_config.R documents", {
  skip_if_missing("write_onset_sample_long"); skip_if_missing(".load_dhis2_delay_params")
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  write_onset_sample_long(.wl_fit(), tmp, "LINELIST_NEW", epidist = .wl_epi(7.5))
  .wl_with(list(RUN_EPIDIST = FALSE),
           expect_warning(out <- write_onset_sample_long(.wl_fit(), tmp, "LINELIST_NEW"),
                          "RUN_EPIDIST is off"))
  expect_identical(.wl_q(out, "selected_estimator"), "interval_censored_mle")
  expect_identical(.load_dhis2_delay_params(tmp)$estimator, "interval_censored_mle")
})

test_that("preserved STRATUM rows alone do not make selected_* claim an EpiDist pooled fit", {
  skip_if_missing("write_onset_sample_long"); skip_if_missing(".load_dhis2_delay_params")
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  .wl_with(list(RUN_EPIDIST = TRUE), {
    write_onset_sample_long(.wl_fit(), tmp, "LINELIST_NEW", strata = list(confirmed = .wl_epi(10)))
    expect_warning(out <- write_onset_sample_long(.wl_fit(), tmp, "LINELIST_NEW"), "PRESERVING")
  })
  expect_identical(.wl_q(out, "selected_estimator"), "interval_censored_mle")
  expect_equal(as.numeric(.wl_q(out, "selected_mean")), 6.667)
  expect_identical(.wl_q(out, "selected_family"), "gamma")
  w <- character(0)
  r <- withCallingHandlers(.load_dhis2_delay_params(tmp),
                           warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_identical(r$estimator, "interval_censored_mle")
  expect_length(w, 0L)                                     # no selected_* self-check warning
  expect_equal(.load_dhis2_delay_params(tmp, stratum = "confirmed")$mean, 10)
})

test_that("selected_param_* follow the SELECTED family (lognormal EpiDist), as the loader derives", {
  skip_if_missing("write_onset_sample_long"); skip_if_missing(".load_dhis2_delay_params")
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  out <- write_onset_sample_long(.wl_fit(), tmp, "LINELIST_NEW",
                                 epidist = .wl_epi(7.123, 9.456, family = "lnorm"))
  expect_identical(.wl_q(out, "selected_family"), "lnorm")
  expect_false(any(c("selected_param_shape", "selected_param_rate") %in% out$quantity))
  r <- .load_dhis2_delay_params(tmp)
  expect_identical(r$family, "lnorm")
  expect_equal(as.numeric(.wl_q(out, "selected_param_meanlog")), unname(r$params[["meanlog"]]), tolerance = 1e-5)
  expect_equal(as.numeric(.wl_q(out, "selected_param_sdlog")),   unname(r$params[["sdlog"]]),   tolerance = 1e-5)
  # The same holds after a preserving rewrite.
  .wl_with(list(RUN_EPIDIST = TRUE),
           expect_warning(out2 <- write_onset_sample_long(.wl_fit(), tmp, "LINELIST_NEW"), "PRESERVING"))
  expect_equal(as.numeric(.wl_q(out2, "selected_param_meanlog")), unname(r$params[["meanlog"]]), tolerance = 1e-5)
  expect_false(any(c("selected_param_shape", "selected_param_rate") %in% out2$quantity))
  # A gamma pick with an unusable sd writes NO selected_param_* rather than the MLE's.
  out3 <- write_onset_sample_long(.wl_fit(), tempfile(fileext = ".csv"), "LINELIST_NEW",
                                  epidist = .wl_epi(7.5, NA_real_))
  expect_false(any(grepl("^selected_param_", out3$quantity)))
})

test_that("EpiDist convergence diagnostics are persisted, pooled and per stratum, and the loader tolerates them", {
  skip_if_missing("write_onset_sample_long"); skip_if_missing(".load_dhis2_delay_params")
  skip_if_missing(".epidist_converged")
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  ok  <- .wl_epi(7.5, max_rhat = 1.003, min_ess_bulk = 800, min_ess_tail = 700, n_divergent = 0)
  bad <- .wl_epi(10,  max_rhat = 1.05,  min_ess_bulk = 800, min_ess_tail = 700, n_divergent = 3)
  out <- write_onset_sample_long(.wl_fit(), tmp, "LINELIST_NEW", epidist = ok,
                                 strata = list(confirmed = bad, not_a_case = .wl_epi(6)))
  expect_identical(anyDuplicated(out$quantity), 0L)
  expect_equal(as.numeric(.wl_q(out, "epidist_max_rhat")), 1.003)
  expect_identical(.wl_q(out, "epidist_converged"), "TRUE")
  expect_identical(.wl_q(out, "epidist_n_divergent__confirmed"), "3")
  expect_identical(.wl_q(out, "epidist_converged__confirmed"), "FALSE")
  expect_false(any(grepl("__not_a_case$", out$quantity) & grepl("rhat|ess|diverg|converged", out$quantity)))
  expect_equal(.load_dhis2_delay_params(tmp)$mean, 7.5)
  expect_equal(.load_dhis2_delay_params(tmp, stratum = "confirmed")$mean, 10)
  expect_equal(.load_dhis2_delay_params(tmp, stratum = "not_a_case")$mean, 6)
  expect_true(is.na(.epidist_converged(1.0, NA, 700, 0)))
})

test_that("the figure builder and the estimator use the SAME window when the latest sample has no onset", {
  skip_if_missing("build_dhis2_delay_populations"); skip_if_missing("estimate_dhis2_onset_sample_delay")
  set.seed(4)
  ob <- get0("OUTBREAK_START", ifnotfound = as.Date("2026-04-30"))
  onset <- ob + sample(0:50, 400, replace = TRUE)
  samp  <- onset + pmax(round(stats::rgamma(400, shape = 2, rate = 0.4)), 0L)
  ll <- data.frame(date_of_symptom_onset = c(as.character(onset), NA),
                   date_of_sample_collection = c(as.character(samp), as.character(max(samp) + 4L)),
                   stringsAsFactors = FALSE)
  asof <- max(samp) + 4L
  pop <- build_dhis2_delay_populations(ll, analysis_date = asof)
  fit <- suppressMessages(estimate_dhis2_onset_sample_delay(ll, analysis_date = asof))
  expect_identical(pop$trunc_date, fit$trunc_date)
  expect_identical(sort(pop$delays$onset_sample$delay), sort(fit$delays_windowed))
})

test_that("the ECDF gap measures misfit, not the size of a day's step", {
  f04d <- file.path(here::here(), "spatiotemporal", "04d_delay_fit_figures.R")
  skip_if(!file.exists(f04d), "04d not found")
  skip_if_missing("build_dhis2_delay_populations")       # 04d would otherwise source 04c
  env <- new.env(parent = globalenv())
  suppressMessages(sys.source(f04d, envir = env))
  set.seed(5)
  # Recorded = round(true delay): exactly the censored-MLE convention, P(recorded <= d) = F(d+0.5).
  x <- round(stats::rgamma(20000, shape = 0.9, rate = 0.15))
  p <- c(shape = 0.9, rate = 0.15)
  expect_lt(env$.dfg_ecdf_gap(x, "gamma", p), 0.02)       # true model: sampling noise only
  # A visibly wrong model is still caught.
  expect_gt(env$.dfg_ecdf_gap(x, "gamma", c(shape = 0.9, rate = 0.30)), 0.1)
})
