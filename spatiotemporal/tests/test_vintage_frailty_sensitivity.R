# =============================================================================
# test_vintage_frailty_sensitivity.R — the 2026-10-14 additions:
#   * the real-snapshot ("vintage") information state (22_daily_reissue.R section 3b)
#   * the data-definition switches: CASE_ZONE_BASIS, SITREP_APPEND_SHIFT_DAYS (date_known)
#   * the gamma-frailty hazard and its reference log-likelihood (21c_frailty_models.R)
#   * the earliest-sample-date helper (00_config.R)
#   * the zone bootstrap of the model selection (16_invasion_eval.R)
# Synthetic data only: nothing here reads the line list.
# =============================================================================
library(testthat)
suppressMessages(source(file.path(here::here(), "spatiotemporal", "00_config.R")))
suppressMessages(source(file.path(here::here(), "spatiotemporal", "01_data_prep.R")))
suppressMessages(source(file.path(here::here(), "spatiotemporal", "02_epi_params.R")))
suppressMessages(source(file.path(here::here(), "spatiotemporal", "16_invasion_eval.R")))
suppressMessages(source(file.path(here::here(), "spatiotemporal", "16b_invasion_recalibration.R")))
suppressMessages(source(file.path(here::here(), "spatiotemporal", "22_daily_reissue.R")))
suppressMessages(source(file.path(here::here(), "spatiotemporal", "21c_frailty_models.R")))

# Temporarily replace global bindings (functions or settings) for one test.
with_globals <- function(vals, code) {
  g <- globalenv()
  had <- vapply(names(vals), exists, logical(1), envir = g, inherits = FALSE)
  old <- lapply(names(vals), function(n) if (had[[n]]) get(n, envir = g) else NULL)
  names(old) <- names(vals)
  on.exit(for (n in names(vals)) if (had[[n]]) assign(n, old[[n]], envir = g)
          else rm(list = n, envir = g), add = TRUE)
  for (n in names(vals)) assign(n, vals[[n]], envir = g)
  force(code)
}
snaps <- data.frame(folder = c("LINELIST_09062026", "LINELIST_12062026", "LINELIST_18062026"),
                    date = as.Date(c("2026-06-09", "2026-06-12", "2026-06-18")),
                    stringsAsFactors = FALSE)

# ---- which snapshot stands in for an issue date --------------------------------------------------
test_that("vintage_for_issue takes the newest snapshot on or before the issue date, never a later one", {
  v <- vintage_for_issue(as.Date("2026-06-19"), source = "vintage", max_gap = 2L, snaps = snaps,
                         analysis_date = as.Date("2026-10-09"))
  expect_identical(v$folder, "LINELIST_18062026"); expect_identical(v$gap, 1L)
  v <- vintage_for_issue(as.Date("2026-06-12"), source = "vintage", max_gap = 2L, snaps = snaps,
                         analysis_date = as.Date("2026-10-09"))
  expect_identical(v$folder, "LINELIST_12062026"); expect_identical(v$gap, 0L)
  # 2026-06-11: the 12 Jun snapshot is LATER and must not be used; 09 Jun is 2 days older
  v <- vintage_for_issue(as.Date("2026-06-11"), source = "vintage", max_gap = 2L, snaps = snaps,
                         analysis_date = as.Date("2026-10-09"))
  expect_identical(v$folder, "LINELIST_09062026")
})

test_that("vintage_for_issue returns NULL off vintage mode, beyond the gap, before the archive and at the analysis date", {
  ad <- as.Date("2026-10-09")
  expect_null(vintage_for_issue(as.Date("2026-06-19"), source = "recensor", snaps = snaps, analysis_date = ad))
  expect_null(vintage_for_issue(as.Date("2026-06-21"), source = "vintage", max_gap = 2L, snaps = snaps, analysis_date = ad))
  expect_null(vintage_for_issue(as.Date("2026-06-05"), source = "vintage", max_gap = 2L, snaps = snaps, analysis_date = ad))
  expect_null(vintage_for_issue(as.Date("2026-06-18"), source = "vintage", max_gap = 2L, snaps = snaps,
                                analysis_date = as.Date("2026-06-18")))
})

test_that("list_linelist_snapshots keeps only dated folders that hold the line-list CSV, oldest first", {
  d <- file.path(tempdir(), paste0("snaps_", Sys.getpid())); on.exit(unlink(d, recursive = TRUE))
  for (f in c("LINELIST_18062026", "LINELIST_09062026", "LINELIST_12062026", "LINELIST_bad", "OTHER_01012026"))
    dir.create(file.path(d, f), recursive = TRUE)
  for (f in c("LINELIST_18062026", "LINELIST_09062026", "LINELIST_bad", "OTHER_01012026"))
    writeLines("x", file.path(d, f, "dhis2_processed_linelist.csv"))
  s <- list_linelist_snapshots(d)
  expect_identical(s$folder, c("LINELIST_09062026", "LINELIST_18062026"))
  expect_equal(s$date, as.Date(c("2026-06-09", "2026-06-18")))
  expect_equal(nrow(list_linelist_snapshots(file.path(d, "nope"))), 0L)
})

test_that("load_linelist_vintage processes a snapshot AS OF its date, restores ANALYSIS_DATE and memoises", {
  calls <- 0L
  fake_load <- function(folder = NULL) {
    calls <<- calls + 1L
    data.frame(folder = folder, asof = as.Date(get("ANALYSIS_DATE", envir = globalenv())))
  }
  on.exit(rm(list = grep("^ll_LINELIST_TEST_", ls(.vintage_memo), value = TRUE), envir = .vintage_memo),
          add = TRUE)
  ad0 <- ANALYSIS_DATE
  with_globals(list(load_linelist = fake_load), {
    x <- load_linelist_vintage("LINELIST_TEST_01062026", as.Date("2026-06-01"))
    expect_equal(x$asof, as.Date("2026-06-01"))
    expect_equal(attr(x, "vintage_date"), as.Date("2026-06-01"))
    expect_equal(ANALYSIS_DATE, ad0)
    invisible(load_linelist_vintage("LINELIST_TEST_01062026", as.Date("2026-06-01")))
    expect_identical(calls, 1L)
  })
  # restored even when the load fails
  with_globals(list(load_linelist = function(folder = NULL) stop("boom")), {
    expect_error(load_linelist_vintage("LINELIST_TEST_02062026", as.Date("2026-06-02")), "boom")
    expect_equal(ANALYSIS_DATE, ad0)
  })
})

test_that("load_linelist_vintage leaves the run's onset-imputation summary globals untouched", {
  # load_linelist() publishes .PCT_ONSET_IMPUTED etc. for the RUN; a snapshot replay must not
  # overwrite them (the report quotes them).
  g <- globalenv()
  on.exit(rm(list = grep("^ll_LINELIST_TEST_", ls(.vintage_memo), value = TRUE), envir = .vintage_memo),
          add = TRUE)
  fake_load <- function(folder = NULL) {
    assign(".PCT_ONSET_IMPUTED", -1, envir = g); assign(".N_SITREP_APPENDED", -1L, envir = g)
    assign(".ONSET_IMPUTE_DESC", "vintage", envir = g)
    data.frame(folder = folder)
  }
  with_globals(list(load_linelist = fake_load, .PCT_ONSET_IMPUTED = 24.5,
                    .N_SITREP_APPENDED = 908L, .ONSET_IMPUTE_DESC = "run"), {
    invisible(load_linelist_vintage("LINELIST_TEST_03062026", as.Date("2026-06-03")))
    expect_equal(get(".PCT_ONSET_IMPUTED", envir = g), 24.5)
    expect_identical(get(".N_SITREP_APPENDED", envir = g), 908L)
    expect_identical(get(".ONSET_IMPUTE_DESC", envir = g), "run")
  })
})

test_that("asof_counts_source: snapshot + its real-time truncation, else recensor with the reason", {
  on.exit({ if (exists("snapshots", envir = .vintage_memo)) rm("snapshots", envir = .vintage_memo) },
          add = TRUE)
  assign("snapshots", snaps, envir = .vintage_memo)
  loaded <- character(0); fitted_at <- character(0)
  fake_vintage <- function(folder, vdate) { loaded <<- c(loaded, folder); data.frame(src = "vintage") }
  # A real-time truncation exists only from 2026-06-12 in this fixture.
  fake_trunc <- function(vdate) { fitted_at <<- c(fitted_at, format(vdate))
    if (as.Date(vdate) >= as.Date("2026-06-12")) list(dist = "DIST", spec = list(mean = 6)) else NULL }
  stub <- function(ll, zones, issue_date, week_spine = NULL, censor = TRUE)
    data.frame(src = ll$src[1], issue = issue_date, censor = censor)
  final_ll <- data.frame(src = "final")
  with_globals(list(reaggregate_asof = stub, load_linelist_vintage = fake_vintage,
                    vintage_truncation = fake_trunc,
                    LFO_TRAINING_SOURCE = "vintage", ANALYSIS_DATE = as.Date("2026-10-09")), {
    r <- asof_counts_source(final_ll, "z", as.Date("2026-06-19"), delay_recensor = "ASOF")
    expect_identical(r$source, "vintage"); expect_identical(r$raw$src, "vintage")
    expect_identical(r$delay, list(mean = 6)); expect_equal(r$nowcast_date, as.Date("2026-06-18"))
    expect_equal(r$raw$issue, as.Date("2026-06-19"))                     # the grid still ends at the issue date
    expect_false(r$raw$censor)                                           # a snapshot is NOT re-censored
    expect_identical(loaded, "LINELIST_18062026")
    expect_identical(fitted_at, "2026-06-18")                            # fitted on the SNAPSHOT's date
    # a snapshot whose date predates a feasible real-time fit -> recensor, with the reason
    r <- asof_counts_source(final_ll, "z", as.Date("2026-06-11"), delay_recensor = "ASOF")
    expect_identical(r$source, "recensor"); expect_identical(r$raw$src, "final")
    expect_true(r$raw$censor)                                            # the final line list IS
    expect_identical(r$delay, "ASOF"); expect_identical(r$why, "no real-time truncation")
    # before the archive -> recensor at the issue date
    r <- asof_counts_source(final_ll, "z", as.Date("2026-06-05"), delay_recensor = "ASOF")
    expect_identical(r$source, "recensor"); expect_equal(r$nowcast_date, as.Date("2026-06-05"))
    expect_identical(r$why, "no snapshot")
  })
  with_globals(list(reaggregate_asof = stub, LFO_TRAINING_SOURCE = "recensor"), {
    r <- asof_counts_source(final_ll, "z", as.Date("2026-06-19"), delay_recensor = "ASOF")
    expect_identical(r$source, "recensor"); expect_identical(r$delay, "ASOF")
    expect_identical(r$why, "recensor requested")
    expect_equal(nrow(prefetch_vintages(as.Date("2026-06-19"))), 0L)
  })
})

test_that("a snapshot is aggregated as it stands; only the final line list is re-censored", {
  # A snapshot dated 2026-06-12 holding a record whose sample date was mistyped AFTER the
  # snapshot (06-14, inside the run's analysis window). Re-censoring at the issue date deleted it,
  # although the live forecast of that day counted it.
  iss <- as.Date("2026-06-12")
  wa  <- as.integer(lubridate::wday(iss + 1L, week_start = 1L))       # the issue date ends its week
  vll <- data.frame(health_zone = "Bunia",
                    date_index = as.Date(c("2026-06-01", "2026-06-08", "2026-06-10")),
                    confirmed = TRUE, suspected = FALSE, final_mve_case_classification = "Confirme",
                    date_of_sample_collection = as.Date(c("2026-06-03", "2026-06-09", "2026-06-14")),
                    stringsAsFactors = FALSE)
  attr(vll, "vintage_date") <- iss
  with_globals(list(WEEK_ANCHOR = wa, ANALYSIS_DATE = as.Date("2026-10-09")), {
    kept  <- suppressMessages(suppressWarnings(reaggregate_asof(vll, "Bunia", iss, censor = FALSE)))
    cens  <- suppressMessages(suppressWarnings(reaggregate_asof(vll, "Bunia", iss)))
    live  <- suppressMessages(suppressWarnings(aggregate_to_zone_week(vll, "Bunia", asof = iss)))
    expect_equal(sum(kept$confirmed), 3L)
    expect_equal(sum(kept$confirmed), sum(live$confirmed))
    expect_equal(sum(cens$confirmed), 2L)                                 # what the old snapshot branch did
    # The plausibility ceiling of a replayed snapshot is its OWN date, not the run's analysis date.
    obs <- suppressMessages(suppressWarnings(linelist_observation_date(vll, iss)))
    expect_true(is.na(obs[3]))
    attr(vll, "vintage_date") <- NULL
    expect_equal(linelist_observation_date(vll, iss)[3], as.Date("2026-06-14"))
    expect_error(reaggregate_asof(vll, "Bunia", iss, censor = NA), "TRUE or FALSE")
  })
  with_globals(list(OUTBREAK_START = NULL),
               expect_error(linelist_observation_date(vll, iss), "OUTBREAK_START is not set"))
})

test_that("issue-date windows carry the model's own p when they ARE its horizon windows", {
  cut <- as.Date("2026-10-03")
  fc <- data.frame(method = "M", health_zone = rep(c("A", "B", "C"), 2), horizon = rep(1:2, each = 3),
                   mu_forecast = c(0.2, 0.5, NA, 0.35, 0.9, NA),
                   p_invasion = c(0.17, 0.36, NA, 0.27, 0.55, NA),      # posterior means != 1 - exp(-mu)
                   was_active_before = rep(c(FALSE, FALSE, TRUE), 2), stringsAsFactors = FALSE)
  an <- anchor_windows_from_analysis_date(fc, cut + 6, cut)              # a == 6: the deployed geometry
  expect_equal(an$p_invasion, fc$p_invasion)
  expect_equal(an$p_case_invasion, fc$p_invasion)
  expect_equal(an$mu_forecast, fc$mu_forecast)
  expect_equal(an$window_end, rep(cut + 6 + c(7L, 14L), each = 3))
  # mid-week, or without a model probability: the plug-in at the windowed hazard
  an5 <- suppressWarnings(anchor_windows_from_analysis_date(fc, cut + 5, cut))
  expect_equal(an5$p_invasion, 1 - exp(-an5$mu_forecast))
  an0 <- anchor_windows_from_analysis_date(fc[, names(fc) != "p_invasion"], cut + 6, cut)
  expect_equal(an0$p_invasion, 1 - exp(-an0$mu_forecast))
  expect_error(anchor_windows_from_analysis_date(fc, cut + 7, cut), "must be 0-6")
})

# ---- extract truncation panel: the date_index series of each replayed vintage --------------------
test_that(".trunc_frame_vintage tabulates CONFIRMED records by date_index inside [start, vintage]", {
  d <- file.path(tempdir(), paste0("arch_", Sys.getpid())); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  dir.create(file.path(d, "LINELIST_10062026"), recursive = TRUE)
  writeLines("x", file.path(d, "LINELIST_10062026", "dhis2_processed_linelist.csv"))
  fake <- data.frame(
    confirmed = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    date_index = as.Date(c("2026-06-01", "2026-06-01", "2026-06-09", "2026-06-05", "2026-05-20", NA)),
    date_of_symptom_onset = as.Date(NA))                     # no recorded onsets: basis must be date_index
  got <- NULL
  with_globals(list(LINELIST_DIR = d, load_linelist_vintage = function(folder, vdate) { got <<- c(folder, format(vdate)); fake }), {
    fr <- .trunc_frame_vintage(d, as.Date("2026-06-10"), as.Date("2026-06-01"))
    expect_identical(got, c("LINELIST_10062026", "2026-06-10"))
    expect_equal(fr$date, seq(as.Date("2026-06-01"), as.Date("2026-06-10"), by = "day"))
    expect_equal(fr$confirm, c(2L, rep(0L, 7), 1L, 0L))      # not-confirmed, pre-start and NA excluded
    expect_null(.trunc_frame_vintage(d, as.Date("2026-06-11"), as.Date("2026-06-01")))   # no such folder
  })
  with_globals(list(LINELIST_DIR = d, load_linelist_vintage = function(folder, vdate) stop("unreadable")), {
    expect_message(expect_null(.trunc_frame_vintage(d, as.Date("2026-06-10"), as.Date("2026-06-01"))), "dropped")
  })
  with_globals(list(LINELIST_DIR = tempdir(), load_linelist_vintage = function(folder, vdate) fake),
    expect_error(.trunc_frame_vintage(d, as.Date("2026-06-10"), as.Date("2026-06-01")), "LINELIST_DIR"))
})

# ---- data-definition switches --------------------------------------------------------------------
test_that("the published definition is the default and stamps nothing", {
  expect_identical(CASE_ZONE_BASIS, "reporting")
  expect_identical(SITREP_APPEND_SHIFT_DAYS, 0L)
  expect_identical(LFO_TRAINING_SOURCE, "vintage")
  expect_identical(ONSET_MODE, "impute")
  expect_identical(DATA_DEFINITION_TAG, "")
  expect_null(nondefault_data_definition())
  expect_true(stop_if_data_definition_sensitivity("x.R"))
  with_globals(list(DATA_DEFINITION_TAG = "residence"),
               expect_error(stop_if_data_definition_sensitivity("x.R"), "published tree"))
})

test_that("date_known, where present, is the observation date; elsewhere the usual priority applies", {
  ll <- data.frame(date_of_sample_collection = as.Date(c("2026-06-01", "2026-06-05", NA)),
                   date_of_notification = as.Date(c(NA, NA, "2026-06-07")),
                   date_known = as.Date(c(NA, "2026-06-10", NA)))
  obs <- linelist_observation_date(ll, as.Date("2026-06-30"))
  expect_equal(obs, as.Date(c("2026-06-01", "2026-06-10", "2026-06-07")))
  ll$date_known <- NULL
  expect_equal(linelist_observation_date(ll, as.Date("2026-06-30")),
               as.Date(c("2026-06-01", "2026-06-05", "2026-06-07")))
})

test_that("residence basis moves only confirmed cases with an unambiguous in-spine residence", {
  spine <- readr::read_csv(file.path(WORLDPOP_DIR, "worldpop__pop_count__static.csv"),
                           col_types = readr::cols(nom = "c", pop_count = "d"), show_col_types = FALSE)$nom
  z <- spine[1:3]
  lookup <- data.frame(health_area = c("A", "A ", "B", "B", "C"),
                       health_area_unmatched = c(FALSE, FALSE, FALSE, FALSE, TRUE),
                       health_zone = c(z[1], z[1], z[2], z[3], z[2]), stringsAsFactors = FALSE)
  ll <- data.frame(
    health_zone = c(z[2], z[2], z[3], z[2], z[3], z[2]),
    province    = c("P2", "P2", "P3", "P2", "P3", "P2"),
    final_mve_case_classification = c(CONFIRMED_STATUS, CONFIRMED_STATUS, CONFIRMED_STATUS,
                                      "suspect", CONFIRMED_STATUS, CONFIRMED_STATUS),
    permanent_residence_health_area = c("A", "B", " A", "A", NA, "C"), stringsAsFactors = FALSE)
  ll <- rbind(ll, data.frame(health_zone = z[1], province = "P1",
                             final_mve_case_classification = "suspect",
                             permanent_residence_health_area = NA))
  out <- suppressMessages(.apply_residence_zone(ll, lookup))
  # A -> z1 (unambiguous, "A " trimmed); B ambiguous; C's only pair is unmatched; suspects stay
  expect_identical(out$health_zone, c(z[1], z[2], z[1], z[2], z[3], z[2], z[1]))
  expect_identical(out$health_zone_reporting, ll$health_zone)
  expect_identical(out$province[c(1, 3)], c("P1", "P1"))                 # province follows the zone
  ll$permanent_residence_health_area <- NULL
  expect_warning(out <- .apply_residence_zone(ll, lookup), "keeps its reporting zone")
  expect_identical(out$health_zone, ll$health_zone)
})

test_that("residence basis is NA-safe: a confirmed record with no reporting zone takes its residence", {
  spine <- readr::read_csv(file.path(WORLDPOP_DIR, "worldpop__pop_count__static.csv"),
                           col_types = readr::cols(nom = "c", pop_count = "d"), show_col_types = FALSE)$nom
  z <- spine[1:2]
  lookup <- data.frame(health_area = c("A", "B"), health_area_unmatched = FALSE,
                       health_zone = z, stringsAsFactors = FALSE)
  ll <- data.frame(health_zone = c(NA, NA, z[2]), province = c(NA, "P", "P2"),
                   final_mve_case_classification = c(CONFIRMED_STATUS, CONFIRMED_STATUS, CONFIRMED_STATUS),
                   permanent_residence_health_area = c("A", NA, "A"), stringsAsFactors = FALSE)
  out <- suppressMessages(.apply_residence_zone(ll, lookup))     # used to error on the NA zone
  expect_identical(out$health_zone, c(z[1], NA, z[1]))
  expect_identical(out$health_zone_reporting, ll$health_zone)
})

# ---- frailty -------------------------------------------------------------------------------------
test_that("the frailty hazard tends to the fixed-beta hazard as theta -> 0 and never exceeds it", {
  set.seed(4); h1 <- rexp(50); h2 <- rexp(50); C <- rexp(50, 0.2)
  h0 <- frailty_cum_hazard(h1, h2, C, theta = NULL)
  expect_equal(h0$H1, h1); expect_equal(h0$H2, h1 + h2)
  ht <- frailty_cum_hazard(h1, h2, C, theta = 1e-9)
  expect_equal(ht$H1, h0$H1, tolerance = 1e-6); expect_equal(ht$H2, h0$H2, tolerance = 1e-6)
  hf <- frailty_cum_hazard(h1, h2, C, theta = 1.5)
  expect_true(all(hf$H1 <= h0$H1 + 1e-12))
  expect_true(all(hf$H2 >= hf$H1))                                        # cumulative horizons nest
  hC <- frailty_cum_hazard(h1, h2, C + 5, theta = 1.5)
  expect_true(all(hC$H1 < hf$H1))                                          # depletion: more past exposure, lower hazard
  # the 2-week window is the joint survival: S(2)/S(0) = (1 + theta (h1 + h2)/(1 + theta C))^(-1/theta)
  expect_equal(exp(-hf$H2), (1 + 1.5 * (h1 + h2) / (1 + 1.5 * C))^(-1 / 1.5))
})

test_that("the reference log-likelihood reduces to the cloglog Bernoulli and matches the closed form", {
  set.seed(5); S <- 4L; Z <- 6L; K <- 1L
  L  <- matrix(rexp(Z * S), Z, S)
  XA <- lapply(seq_len(S), function(s) matrix(rnorm(Z * K), Z, K))     # covariates vary by week
  rows <- expand.grid(zone = seq_len(Z), tr = seq_len(S))
  lam <- L[cbind(rows$zone, rows$tr)]
  X <- matrix(vapply(seq_len(nrow(rows)), function(n) XA[[rows$tr[n]]][rows$zone[n], 1], 0), ncol = 1)
  y <- rbinom(nrow(rows), 1, 0.3)
  par <- list(b0 = -1, beta = 0.4)
  H <- exp(par$b0 + drop(X %*% par$beta)) * lam
  expect_equal(frailty_loglik_R(par, y, X, lam, rows$zone, rows$tr, L, XA),
               sum(ifelse(y == 1, log(1 - exp(-H)), -H)))
  # frailty + drift: each past week enters C at ITS OWN covariates and drift level
  par$lth <- log(0.8); par$u <- c(0, 0.2, -0.1, 0.3)
  h <- exp(par$b0 + par$u[rows$tr] + drop(X %*% par$beta)) * lam
  Cz <- vapply(seq_len(nrow(rows)), function(k) {
    ss <- seq_len(rows$tr[k] - 1L)
    sum(vapply(ss, function(s) exp(par$b0 + par$u[s] + XA[[s]][rows$zone[k], 1] * par$beta) *
                                 L[rows$zone[k], s], 0)) }, 0)
  Hf <- frailty_cum_hazard(h, 0, Cz, theta = 0.8)$H1
  expect_equal(frailty_loglik_R(par, y, X, lam, rows$zone, rows$tr, L, XA),
               sum(ifelse(y == 1, log(-expm1(-Hf)), -Hf)))
})


test_that("a confirmed record takes its confirming test's lab date; other classes keep theirs", {
  es <- data.frame(alert_id = c("a", "b", "c"),
                   min_sample = as.Date(c("2026-06-01", "2026-06-20", "2026-06-03")),
                   min_lab_confirmed = as.Date(c("2026-06-02", NA, "2026-06-04")))
  attr(es, "n_snapshots") <- 3L
  df <- data.frame(alert_id = c("a", "b", "c", "d"),
                   date_of_sample_collection = as.Date(c("2026-06-05", "2026-06-10", NA, "2026-06-07")),
                   lab_analysis_date = as.Date(c("2026-06-15", "2026-06-11", "2026-06-20", "2026-06-08")),
                   final_mve_case_classification = c(CONFIRMED_STATUS, CONFIRMED_STATUS, "suspect",
                                                     CONFIRMED_STATUS), stringsAsFactors = FALSE)
  out <- suppressMessages(apply_earliest_sample_dates(df, es = es))
  # a: confirmed, earlier confirming lab date -> replaced; b: none recorded -> kept;
  # c: not confirmed now -> kept; d: not in the archive -> kept
  expect_equal(out$lab_analysis_date, as.Date(c("2026-06-02", "2026-06-11", "2026-06-20", "2026-06-08")))
  expect_equal(out$lab_analysis_date_extract, df$lab_analysis_date)
  # an archive written before the lab rule (no min_lab_confirmed) leaves the lab date alone
  out0 <- suppressMessages(apply_earliest_sample_dates(df, es = es[, c("alert_id", "min_sample")]))
  expect_equal(out0$lab_analysis_date, df$lab_analysis_date)
})

# ---- earliest sample date ------------------------------------------------------------------------
test_that("the earliest sample date replaces only later or missing dates, for Date and character columns", {
  es <- data.frame(alert_id = c("a", "b", "c"), min_sample = as.Date(c("2026-06-01", "2026-06-20", "2026-06-03")))
  attr(es, "n_snapshots") <- 3L
  df <- data.frame(alert_id = c("a", "b", "c", "d"),
                   date_of_sample_collection = as.Date(c("2026-06-05", "2026-06-10", NA, "2026-06-07")),
                   final_mve_case_classification = CONFIRMED_STATUS, stringsAsFactors = FALSE)
  out <- suppressMessages(apply_earliest_sample_dates(df, es = es))
  expect_equal(out$date_of_sample_collection, as.Date(c("2026-06-01", "2026-06-10", "2026-06-03", "2026-06-07")))
  expect_equal(out$date_of_sample_collection_extract, df$date_of_sample_collection)
  dc <- df; dc$date_of_sample_collection <- format(dc$date_of_sample_collection, "%Y-%m-%d")
  outc <- suppressMessages(apply_earliest_sample_dates(dc, es = es))
  expect_identical(outc$date_of_sample_collection, c("2026-06-01", "2026-06-10", "2026-06-03", "2026-06-07"))
})

# ---- selection bootstrap -------------------------------------------------------------------------
test_that("the selection bootstrap is reproducible, sums to one, favours the better model and spares the RNG", {
  set.seed(6)
  cuts <- as.Date("2026-06-06") + 7L * (0:7)                # post-burn-in rounds
  zones <- sprintf("Z%02d", 1:60)
  base <- expand.grid(health_zone = zones, fold_id = seq_along(cuts), horizon = 1:2, stringsAsFactors = FALSE)
  base$cutoff <- cuts[base$fold_id]
  risk <- stats::setNames(rexp(length(zones)), zones)
  base$is_new_invasion <- rbinom(nrow(base), 1, pmin(0.9, 0.04 * risk[base$health_zone] * base$horizon))
  mk <- function(m, good) transform(base, method = m,
    p_invasion = pmin(0.95, 0.04 * (if (good) risk[base$health_zone] else rexp(nrow(base))) * base$horizon))
  lfo <- rbind(mk("Bayes-M8-med", TRUE), mk("Bayes-M14-fill-geo", FALSE))
  lfo$eval_reliable <- TRUE; lfo$was_active_before <- FALSE
  set.seed(99); r0 <- runif(1)
  set.seed(99); b1 <- invasion_selection_bootstrap(lfo, B = 30L, select_on_recal = FALSE, seed = 1L); r1 <- runif(1)
  expect_identical(r0, r1)
  b2 <- invasion_selection_bootstrap(lfo, B = 30L, select_on_recal = FALSE, seed = 1L)
  expect_identical(b1, b2)
  expect_equal(sum(b1$share), 1)
  expect_identical(b1$method[1], "Bayes-M8-med")
})

test_that(".vintage_settings_key moves with every setting a snapshot replay reads", {
  k0 <- .vintage_settings_key()
  for (v in list(list(APPEND_SITREP_CONFIRMED = !isTRUE(get0("APPEND_SITREP_CONFIRMED", ifnotfound = TRUE))),
                 list(WEEK_ANCHOR = (get0("WEEK_ANCHOR", ifnotfound = 1L) %% 7L) + 1L),
                 list(ONSET_SAMPLE_NEG_TOL_DAYS = get0("ONSET_SAMPLE_NEG_TOL_DAYS", ifnotfound = 2L) + 1L),
                 list(DELAY_MAX_PLAUSIBLE_DAYS = get0("DELAY_MAX_PLAUSIBLE_DAYS", ifnotfound = 60L) + 1L),
                 list(ONSET_IMPUTE_SEED = get0("ONSET_IMPUTE_SEED", ifnotfound = 1L) + 1L),
                 list(CASE_ZONE_BASIS = if (identical(get0("CASE_ZONE_BASIS"), "residence")) "reporting"
                                        else "residence")))
    with_globals(v, expect_false(identical(.vintage_settings_key(), k0), info = names(v)))
  expect_identical(.vintage_settings_key(), k0)
})
