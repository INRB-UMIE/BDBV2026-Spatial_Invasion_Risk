# =============================================================================
# tests/test_onset_censoring.R — an onset never postdates its own specimen
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# A confirmed case whose recorded onset falls a day or two AFTER its own specimen
# is recording noise, not presymptomatic detection: on the 2026-09-07 snapshot
# 39.5% of the -2/-1 d records were deceased at swab against 22.4% of the
# positive-delay records, and half the affected confirmed rows are death alerts.
# Such a record is KEPT (rejecting it replaces a 1-2 d error with an imputed onset
# ~8 d too early, because an unusable onset is imputed rather than dropped) but its
# onset is CENSORED at the specimen date.
#
# Three things can silently revert this, which is why they are pinned here rather
# than left to the loader's own comments:
#   * the tolerance going back to an inline literal, so config and code can drift;
#   * the censoring being dropped, putting date_index after the specimen again and
#     re-opening the deliberate disagreement with the delay-fitting pool, which
#     requires a NON-NEGATIVE delay (01_data_prep.R `.complete`);
#   * the outbreak-floor clamp being dropped, which a specimen dated before
#     OUTBREAK_START would turn into a pre-outbreak week in the zone-week grid.
# =============================================================================

test_that("ONSET_SAMPLE_NEG_TOL_DAYS is a named config constant, not an inline literal", {
  expect_true(exists("ONSET_SAMPLE_NEG_TOL_DAYS"))
  expect_true(is.numeric(ONSET_SAMPLE_NEG_TOL_DAYS))
  expect_length(ONSET_SAMPLE_NEG_TOL_DAYS, 1L)
  expect_gte(ONSET_SAMPLE_NEG_TOL_DAYS, 0L)

  src <- readLines(file.path(here::here(), "spatiotemporal", "01_data_prep.R"), warn = FALSE)
  rule <- grep("date_of_symptom_onset <= date_of_sample_collection", src, value = TRUE)
  expect_length(rule, 1L)
  expect_match(rule, "\\.neg_tol",
               info = "the usability rule must read the named tolerance, not a literal")
})

test_that("no USABLE onset survives after its own specimen, and the floor holds", {
  if (!exists("load_linelist", mode = "function")) testthat::skip("load_linelist() unavailable")
  ll <- tryCatch(suppressMessages(load_linelist()), error = function(e) NULL)
  if (is.null(ll)) testthat::skip("line list unavailable")

  # Scoped to USABLE onsets on purpose. Censoring governs the branch where a recorded onset
  # is believed; it says nothing about the imputed branch, which is governed by the outbreak
  # floor clamp below and legitimately produces a date_index after the specimen.
  us <- ll$onset_usable %in% TRUE & !is.na(ll$date_of_sample_collection) & !is.na(ll$date_index)
  expect_true(all(ll$date_index[us] <= ll$date_of_sample_collection[us]),
              info = "a usable onset must be censored at, never left after, its specimen")

  expect_true(all(ll$date_index >= OUTBREAK_START, na.rm = TRUE),
              info = "date_index must never precede the outbreak floor (zone-week grid invariant)")
  asof <- suppressWarnings(as.Date(get0("ANALYSIS_DATE", ifnotfound = NA)))
  if (!is.na(asof)) expect_true(all(ll$date_index <= asof, na.rm = TRUE))
})

test_that("the only records dated after their specimen are imputed ones lifted by the floor", {
  # This is pre-existing, intended behaviour of `pmax(sample - delay, .ob_floor)`: a specimen
  # collected before OUTBREAK_START cannot yield an in-range onset, so the record is lifted to
  # the floor. Pinned rather than merely excluded above, so that if the set ever grows — or
  # ever contains a CONFIRMED case, which would put a fabricated date into the invasion
  # outcome — the suite says so instead of the exclusion quietly absorbing it.
  if (!exists("load_linelist", mode = "function")) testthat::skip("load_linelist() unavailable")
  ll <- tryCatch(suppressMessages(load_linelist()), error = function(e) NULL)
  if (is.null(ll)) testthat::skip("line list unavailable")

  after <- which(!is.na(ll$date_of_sample_collection) & !is.na(ll$date_index) &
                 ll$date_index > ll$date_of_sample_collection)
  if (!length(after)) succeed("no record is dated after its specimen at all")
  expect_true(all(ll$onset_imputed[after] %in% TRUE),
              info = "only the imputed branch may be dated after its specimen")
  expect_true(all(ll$date_index[after] == OUTBREAK_START),
              info = "and only by being clamped to the outbreak floor")
  expect_true(all(ll$date_of_sample_collection[after] < OUTBREAK_START),
              info = "which can only happen when the specimen predates the floor")
  expect_false(any(ll$confirmed[after] %in% TRUE),
               info = "a CONFIRMED case lifted to the floor would inject a fabricated invasion date")
})

test_that("censoring moves only the intended records, and only by the tolerance", {
  if (!exists("load_linelist", mode = "function")) testthat::skip("load_linelist() unavailable")
  ll <- tryCatch(suppressMessages(load_linelist()), error = function(e) NULL)
  if (is.null(ll)) testthat::skip("line list unavailable")
  expect_true("onset_censored" %in% names(ll))

  del <- as.numeric(ll$date_of_sample_collection - ll$date_of_symptom_onset)
  expected <- ll$onset_usable %in% TRUE & !is.na(del) & del < 0
  expect_identical(which(ll$onset_censored %in% TRUE), which(expected),
                   info = "censored set must be exactly the usable records with a negative delay")

  cz <- which(ll$onset_censored %in% TRUE)
  if (length(cz)) {
    shift <- as.numeric(ll$date_of_symptom_onset[cz] - ll$date_index[cz])
    expect_true(all(shift >= 1 & shift <= ONSET_SAMPLE_NEG_TOL_DAYS),
                info = "no censored record may move by more than the tolerance")
    expect_false(any(ll$onset_imputed[cz] %in% TRUE),
                 info = "a censored record is kept, never routed into imputation")
  }
})

test_that("censoring changes no zone's first confirmed onset", {
  # The whole forecast target is the WEEK of a zone's first confirmed onset. If this
  # ever stops holding, the tolerance has become a modelling choice rather than a
  # data-cleaning one and has to be justified as such in the manuscript.
  if (!exists("load_linelist", mode = "function")) testthat::skip("load_linelist() unavailable")
  ll <- tryCatch(suppressMessages(load_linelist()), error = function(e) NULL)
  if (is.null(ll)) testthat::skip("line list unavailable")

  idx_uncensored <- ifelse(ll$onset_usable %in% TRUE,
                           ll$date_of_symptom_onset, ll$date_index)
  idx_uncensored <- as.Date(idx_uncensored, origin = "1970-01-01")
  keep <- ll$confirmed %in% TRUE & !is.na(ll$health_zone) & !is.na(ll$date_index)
  if (!any(keep)) testthat::skip("no confirmed records")

  z  <- ll$health_zone[keep]
  fa <- tapply(ll$date_index[keep],   z, min)
  fb <- tapply(idx_uncensored[keep],  z, min)
  expect_identical(unname(fa), unname(fb),
                   info = "a zone's first confirmed onset must not depend on the censoring rule")
})
