# =============================================================================
# tests/test_delay_strata.R — the onset-to-sample delay, fitted per case class
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# The REPORTED onset->sample delay is the CONFIRMED-case delay (since 2026-10-15 it no longer
# drives the onset imputation, which draws same-class donor intervals). A fit pooled over every
# classification is mostly TEST-NEGATIVE and unclassified records on this line list (the
# measured composition is in estimate_onset_sample_strata()'s docstring, 04c), so the pooled
# and confirmed-stratum delays differ.
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

test_that("the reported delay resolves the CONFIRMED stratum by default", {
  # BEHAVIOURAL, not a source-text grep. This used to grep 01_data_prep.R for the literal
  # `effective_onset_sample_delay(stratum = "confirmed")`; the call became
  # `stratum = .imp_stratum` (resolved from ONSET_IMPUTE_STRATUM) and the test failed while
  # the behaviour was correct. Assert what actually matters: the configured default is the
  # confirmed stratum, and the resolver returns THAT stratum's fit rather than the pooled one.
  .skip_unless("effective_onset_sample_delay")
  expect_identical(get0("ONSET_IMPUTE_STRATUM", ifnotfound = NULL), "confirmed")

  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(
    quantity = c("epidist_family", "epidist_mean", "epidist_sd", "epidist_n",
                 "epidist_family__confirmed", "epidist_mean__confirmed",
                 "epidist_sd__confirmed", "epidist_n__confirmed",
                 "selected_family", "selected_mean"),
    value    = c("gamma", "7.665", "8.507", "13987",
                 "gamma", "10.136", "10.416", "4981",
                 "gamma", "7.665")),
    tmp, row.names = FALSE)

  pooled <- effective_onset_sample_delay(tmp, stratum = NULL)
  conf   <- effective_onset_sample_delay(tmp, stratum = "confirmed")
  expect_equal(pooled$mean, 7.665,  tolerance = 1e-9)
  expect_equal(conf$mean,  10.136, tolerance = 1e-9)
  expect_identical(conf$stratum, "confirmed")
  # The two must be DIFFERENT distributions, or the stratification is inert — which is
  # precisely the state this machinery was in before 2026-09-24.
  expect_false(isTRUE(all.equal(pooled$mean, conf$mean)))
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
  # This fixture's shape is below 1 (unbounded density at 0); the live fits sit near or below
  # 1, which is why the panel plots a CDF and not a density.
  expect_lt(sp$params[["shape"]], 1)
})

# ---------------------------------------------------------------------------
# POSTERIOR-DRAW REGISTRY KEYS
#
# 04c caches each EpiDist fit's draws under delay|model_type|family so 04d can draw the
# posterior the pipeline actually used without a second Stan run. A put overwrites silently.
# Between 2026-09-22 and 2026-10-02 every stratum fit reused the POOLED key, so the last
# stratum fitted (suspected/probable, n=153) replaced the pooled posterior and
# dhis2_delay_epidist_fits_*.png drew a 4.59 d onset->sample marginal while the pipeline
# applied 7.41 d pooled / 9.81 d confirmed. Nothing modelled was wrong -- the published
# parameters come from the fit objects, not the registry -- but the diagnostic figure, the
# fit-summary CSV and the auto-generated sentence in BAYESIAN_INVASION_REPORT.md all
# described a delay no downstream step applies.
# ---------------------------------------------------------------------------

test_that("a stratum fit never shares the pooled delay's registry key", {
  .skip_unless(".epidist_draw_key")
  expect_identical(.epidist_draw_key("onset_sample"), "onset_sample")
  expect_identical(.epidist_draw_key("onset_sample", classes = NULL, stratum = "confirmed"),
                   "onset_sample")                       # no subset -> no qualification
  expect_identical(.epidist_draw_key("onset_sample", classes = "confirmed_case",
                                     stratum = "confirmed"), "onset_sample__confirmed")
  # The qualification comes from `classes`, so forgetting the label cannot restore the
  # collision -- that omission is exactly how a future caller would reintroduce it.
  expect_identical(.epidist_draw_key("onset_sample",
                                     classes = c("suspected_case", "probable_case")),
                   "onset_sample__probable_case+suspected_case")
  keys <- c(.epidist_draw_key("onset_sample"),
            .epidist_draw_key("onset_sample", "confirmed_case", "confirmed"),
            .epidist_draw_key("onset_sample", "not_a_case", "not_a_case"),
            .epidist_draw_key("onset_sample", c("suspected_case", "probable_case"), "suspected"))
  expect_identical(anyDuplicated(keys), 0L)
})

test_that("fit_epidist_both forwards draw_key to every sub-fit", {
  .skip_unless("fit_epidist_both"); .skip_unless("fit_epidist_model")
  # fit_epidist_both() resolves fit_epidist_model() in the global environment, so the stub
  # has to be installed there and restored afterwards.
  seen <- character(0)
  .orig <- get("fit_epidist_model", envir = globalenv())
  assign("fit_epidist_model", function(df_pairs, obs_date_val, model_type, family_fn,
                                       family_name, delay_name, draw_key = delay_name) {
    seen <<- c(seen, .epidist_draws_key(draw_key, model_type, family_name))
    NULL
  }, envir = globalenv())
  on.exit(assign("fit_epidist_model", .orig, envir = globalenv()), add = TRUE)

  pairs <- tibble::tibble(onset = as.Date("2026-05-01") + 0:9,
                          sample = as.Date("2026-05-01") + 0:9 + 3L)
  invisible(fit_epidist_both(pairs, "onset_sample", as.Date("2026-07-01")))
  pooled <- seen; seen <- character(0)
  invisible(fit_epidist_both(pairs, "onset_sample", as.Date("2026-07-01"),
                             draw_key = "onset_sample__confirmed"))
  strat <- seen

  expect_length(pooled, 4L)                               # naive/marginal x lognormal/gamma
  expect_length(strat, 4L)
  expect_true(all(startsWith(pooled, "onset_sample|")))
  expect_true(all(startsWith(strat, "onset_sample__confirmed|")))
  expect_identical(length(intersect(pooled, strat)), 0L)  # the collision that broke 04d
})
