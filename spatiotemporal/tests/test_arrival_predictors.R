# =============================================================================
# tests/test_arrival_predictors.R — Figure 1C / S1 must only VISUALISE
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Figure 1C reports an R^2, a Spearman rho and an n for each predictor of
# arrival time. Those are published results. They used to be computed inside
# make_manuscript_figures.R from raw files it re-read itself — including zone
# coordinates scraped by regex out of an HTML page — so the printed correlation
# could not be reproduced from any table the pipeline writes, and in fact
# disagreed with every candidate table on both the zone set and the arrival
# dates. 43_spread_kinematics.R now publishes both the per-zone predictors and
# the per-predictor fit, and the figure reads them.
#
# These tests pin the contract that makes the panel auditable:
#   1. every statistic on the panel is recomputable from the per-zone table;
#   2. the drawn regression line is the published (intercept, slope);
#   3. n is per predictor, because the predictors are NOT on common support;
#   4. the arrival origin is one date for the whole outbreak;
#   5. the epicentre set enters no fit, and the *_common columns are one zone set;
#   6. the right-truncated tail (in_fit_window = FALSE) enters no fit;
#   7. the short-trip shares come from the raw-resolved profile at build_M1's snapshot.
# They read the published CSVs and skip when the pipeline has not been run.
# ARRIVAL_PREDICTORS_DIR points them at another key_outputs directory (e.g. a scratch run
# of 43) without touching the published tree.
# =============================================================================

.arr_paths <- function() {
  d <- Sys.getenv("ARRIVAL_PREDICTORS_DIR",
                  file.path(here::here(), "spatiotemporal", "outputs", "key_outputs"))
  list(dat = file.path(d, "arrival_predictors.csv"),
       fit = file.path(d, "arrival_predictor_fits.csv"))
}
.skip_unless_arr <- function() {
  p <- .arr_paths()
  if (!file.exists(p$dat) || !file.exists(p$fit))
    testthat::skip("arrival_predictor*.csv not present — run 43_spread_kinematics.R")
  lapply(p, function(f) suppressWarnings(readr::read_csv(f, show_col_types = FALSE)))
}

# THE EPICENTRE SET IS IN THE TABLE BUT IN NO FIT (43_spread_kinematics.R, Section 5b). Its
# zones are the source every predictor is measured from — distance 0, arrival day 0-2, outflow
# share undefined — and they used to enter the distance and travel-time fits while dropping out
# of the short-trip one as NA, so the predictors were compared on different supports. A table
# without the `epicentre` flag predates that fix: re-run 43_spread_kinematics.R.
.need_epi_flag <- function(dat) {
  if (!"epicentre" %in% names(dat))
    stop("arrival_predictors.csv has no `epicentre` column: it predates the ",
         "epicentre-exclusion fix. Re-run 43_spread_kinematics.R.", call. = FALSE)
  as.logical(dat$epicentre) %in% TRUE
}

# THE RIGHT-TRUNCATED TAIL IS IN THE TABLE BUT IN NO FIT (43, Section 5b): a zone first seen
# after the last fit-window day is flagged in_fit_window = FALSE, as in spread_kinematics_zones.csv.
# The fits used to include it although 43's header said truncated weeks enter no fit. A table
# without the flag predates that fix: re-run 43_spread_kinematics.R.
.need_fit_flag <- function(dat) {
  if (!"in_fit_window" %in% names(dat))
    stop("arrival_predictors.csv has no `in_fit_window` column: it predates the ",
         "right-truncation exclusion. Re-run 43_spread_kinematics.R.", call. = FALSE)
  as.logical(dat$in_fit_window) %in% TRUE
}

test_that("every published arrival-predictor statistic is reproducible from the per-zone table", {
  a <- .skip_unless_arr(); dat <- a$dat; fit <- a$fit
  expect_true(nrow(fit) > 0L)
  epi <- .need_epi_flag(dat); infit <- .need_fit_flag(dat)
  for (i in seq_len(nrow(fit))) {
    p <- fit$predictor[i]
    expect_true(p %in% names(dat), info = paste("predictor column missing:", p))
    y <- dat[[p]]; x <- dat$arrival_days
    y[epi | !infit] <- NA_real_           # neither the epicentre set nor the truncated tail is fitted
    k <- is.finite(x) & is.finite(y)
    # n is the COMPLETE-PAIR count for THIS predictor, not the table's row count
    expect_identical(as.integer(fit$n[i]), as.integer(sum(k)), info = p)
    if (sum(k) < 3L) next
    m <- stats::lm(y[k] ~ x[k])
    expect_equal(fit$r2[i],        summary(m)$r.squared,                      tolerance = 1e-5, info = p)
    expect_equal(fit$intercept[i], unname(stats::coef(m)[1L]),                tolerance = 1e-5, info = p)
    expect_equal(fit$slope[i],     unname(stats::coef(m)[2L]),                tolerance = 1e-5, info = p)
    expect_equal(fit$pearson[i],   stats::cor(x[k], y[k]),                    tolerance = 1e-5, info = p)
    expect_equal(fit$spearman[i],
                 suppressWarnings(stats::cor(x[k], y[k], method = "spearman")), tolerance = 1e-5, info = p)
  }
})

test_that("the predictors are not on common support, so n must be published per predictor", {
  a <- .skip_unless_arr(); dat <- a$dat; fit <- a$fit
  epi <- .need_epi_flag(dat); infit <- .need_fit_flag(dat)
  # The mobility predictors are -log(share) and are NA wherever no flow out of the
  # epicentre was released for the zone. That is selection ON the predictor and it is
  # correlated with arrival, so a single table-wide n would misstate every facet but
  # the geographic ones. If this ever becomes a single common n, the panel's "n ="
  # annotation has stopped carrying information and the guard should be revisited,
  # not deleted.
  expect_true(all(fit$n <= sum(!epi & infit)))
  expect_true(any(fit$n == sum(!epi & infit)))  # geography is complete on the fitted zones
  expect_gt(dplyr::n_distinct(fit$n), 1L)       # mobility is not
})

test_that("the epicentre set is excluded from every fit, and the common-support fit is one zone set", {
  a <- .skip_unless_arr(); dat <- a$dat; fit <- a$fit
  epi <- .need_epi_flag(dat); infit <- .need_fit_flag(dat)
  # The flag is the configured epicentre, resolved against the zone spine — no more, no fewer.
  if (exists("EPICENTRE_ZONES")) expect_setequal(dat$health_zone[epi], EPICENTRE_ZONES)
  expect_true(all(c("support", "n_common", "r2_common", "spearman_common") %in% names(fit)))
  expect_true(all(grepl("outside the epicentre set", fit$support, fixed = TRUE)))
  # ONE common support: the fitted zones (outside the epicentre, inside the fit window) on
  # which EVERY predictor is defined.
  x <- dat$arrival_days
  common <- is.finite(x) & !epi & infit &
    Reduce(`&`, lapply(fit$predictor, function(p) is.finite(dat[[p]])))
  expect_equal(dplyr::n_distinct(fit$n_common), 1L)
  expect_identical(as.integer(fit$n_common[1]), as.integer(sum(common)))
  for (i in seq_len(nrow(fit))) {
    y <- dat[[fit$predictor[i]]]
    if (sum(common) < 3L) next
    expect_equal(fit$r2_common[i], stats::cor(x[common], y[common])^2, tolerance = 1e-5,
                 info = fit$predictor[i])
    expect_equal(fit$spearman_common[i],
                 suppressWarnings(stats::cor(x[common], y[common], method = "spearman")),
                 tolerance = 1e-5, info = fit$predictor[i])
  }
})

test_that("the right-truncated tail is flagged by the fit-window cutoff and enters no fit", {
  a <- .skip_unless_arr(); dat <- a$dat; fit <- a$fit
  infit <- .need_fit_flag(dat)
  # The flag is a calendar cutoff on the first onset: every zone in the window arrived no later
  # than every zone outside it.
  if (any(infit) && any(!infit))
    expect_lt(max(as.Date(dat$first_onset[infit])), min(as.Date(dat$first_onset[!infit])))
  expect_true(all(grepl("right-truncated tail is excluded", fit$support, fixed = TRUE)))
  # No tail zone is inside any fit's n: n never exceeds the in-window, non-epicentre count.
  epi <- .need_epi_flag(dat)
  expect_true(all(fit$n <= sum(infit & !epi)))
})

test_that("the short-trip shares are the raw-resolved profile at build_M1's snapshot", {
  a <- .skip_unless_arr(); dat <- a$dat; fit <- a$fit
  if (!"shorttrip_snapshot" %in% names(dat))
    stop("arrival_predictors.csv has no `shorttrip_snapshot` column: it predates the ",
         "short-trip source fix. Re-run 43_spread_kinematics.R.", call. = FALSE)
  snap <- unique(as.Date(dat$shorttrip_snapshot))
  expect_length(snap, 1L)
  # build_M1's rule: the latest snapshot on or before the analysis date, else the earliest.
  if (exists("FLOWMINDER_ST_DATES") && exists("ANALYSIS_DATE") && length(FLOWMINDER_ST_DATES)) {
    el <- FLOWMINDER_ST_DATES[FLOWMINDER_ST_DATES <= ANALYSIS_DATE]
    expect_equal(snap, if (length(el)) max(el) else min(FLOWMINDER_ST_DATES))
  }
  # Shares are proportions of ONE renormalised profile.
  sh <- dat$share_shorttrip[is.finite(dat$share_shorttrip)]
  expect_true(all(sh > 0 & sh <= 1))
  # The live code reads the profile loader, never the lossy processed matrix (which silently
  # lacks real spine destinations and dropped invaded zones out of this predictor as NA).
  f43 <- file.path(here::here(), "spatiotemporal", "43_spread_kinematics.R")
  if (file.exists(f43)) {
    code <- paste(sub("#.*$", "", readLines(f43, warn = FALSE)), collapse = "\n")
    expect_true(grepl("load_short_trip_profile(", code, fixed = TRUE))
    expect_false(grepl("load_short_trip_snapshot(", code, fixed = TRUE))
  }
})

test_that("arrival days are measured from ONE origin, and that origin is the earliest onset", {
  a <- .skip_unless_arr(); dat <- a$dat
  expect_equal(dplyr::n_distinct(dat$arrival_origin), 1L)
  expect_equal(min(dat$arrival_days), 0)
  expect_true(all(dat$arrival_days >= 0))
  expect_equal(as.numeric(as.Date(dat$first_onset) - as.Date(dat$arrival_origin)),
               as.numeric(dat$arrival_days))
})

test_that("the published fit is finite and orientated predictor-on-arrival", {
  a <- .skip_unless_arr(); fit <- a$fit
  ok <- is.finite(fit$r2)
  expect_true(any(ok))
  expect_true(all(fit$r2[ok] >= 0 & fit$r2[ok] <= 1))
  expect_true(all(abs(fit$spearman[is.finite(fit$spearman)]) <= 1))
  # R^2 is the SQUARE of the Pearson correlation only under a simple linear fit of one
  # variable on the other; this pins that the published pair describes the same fit.
  expect_equal(fit$r2[ok], fit$pearson[ok]^2, tolerance = 1e-6)
})

test_that("figure scripts do not fit the arrival-predictor model themselves", {
  fig <- file.path(here::here(), "spatiotemporal", "make_manuscript_figures.R")
  if (!file.exists(fig)) testthat::skip("make_manuscript_figures.R not present")
  src <- paste(readLines(fig, warn = FALSE), collapse = "\n")
  b0 <- regexpr("build_fig1c <- function", src, fixed = TRUE)
  expect_gt(b0, 0)
  b1 <- regexpr("\nbuild_fig1 <- function", substring(src, b0), fixed = TRUE)
  body <- substring(src, b0, b0 + (if (b1 > 0) b1 else nchar(src)))
  # No model fitting and no smoother that refits: the line comes from the published
  # (intercept, slope) via geom_abline.
  expect_false(grepl("\\blm\\s*\\(", body))
  expect_false(grepl("geom_smooth", body, fixed = TRUE))
  expect_false(grepl("\\bcor\\s*\\(", body))
  expect_true(grepl("geom_abline", body, fixed = TRUE))
})

# =============================================================================
# Sitrep top-up rows must carry a province
# =============================================================================
# The sitrep cumulative file has only (nom, date, cumulative_confirmed_cases), so
# .build_sitrep_confirmed_appends() has to source province from elsewhere. It used
# to leave it NA on all 606 appended rows across 28 zones. Every consumer that
# groups by (health_zone, province) — Figure 1A's bivariate choropleth among them —
# then split a topped-up zone into TWO groups, one holding the line-list cases and
# one the appended cases: the zone was counted twice in the zone total and each
# half was binned on a fraction of its true burden.
test_that("appended sitrep rows carry a province, so a topped-up zone stays one group", {
  if (!exists("load_linelist", mode = "function")) testthat::skip("load_linelist() not available")
  ll <- tryCatch(suppressMessages(suppressWarnings(load_linelist())),
                 error = function(e) NULL)
  if (is.null(ll)) testthat::skip("line list unavailable in this environment")
  app <- ll[grepl("^SITREP-CONF", if (is.null(ll$alert_id)) "" else ll$alert_id), , drop = FALSE]
  if (!nrow(app)) testthat::skip("no sitrep top-up rows on this data frame")

  expect_true("province" %in% names(ll))
  # Every appended row whose zone also appears in the line list must carry that
  # zone's province. Zones present ONLY in the sitrep fall back to the shapefile;
  # if even that fails the builder warns, so NA is tolerated only there.
  ll_prov <- ll[!grepl("^SITREP-CONF", if (is.null(ll$alert_id)) "" else ll$alert_id) & !is.na(ll$province), ]
  known <- unique(ll_prov$health_zone)
  bad <- app[app$health_zone %in% known & is.na(app$province), ]
  expect_equal(nrow(bad), 0L)

  # The grouping invariant that actually matters downstream: one confirmed zone,
  # one (health_zone, province) group.
  cf <- ll[ll$confirmed %in% TRUE & !is.na(ll$date_index) & !is.na(ll$health_zone), ]
  grp <- unique(cf[, c("health_zone", "province")])
  expect_equal(nrow(grp), length(unique(cf$health_zone)))
})

# =============================================================================
# The three structural baselines on Figure 2 / Figure 2-cascade
# =============================================================================
# The user's specification: three baselines that are NOT the renewal process —
# (1) a gravity model, (2) the Flowminder inflow from the epicentre ranked,
# (3) a travel-time matrix, labelled as such. They differ ONLY in the connectivity
# matrix, so the comparison isolates the notion of connectivity.
#
# The set they replaced was wrong on all three counts: "Gravity-B4" has no
# destination mass and unfitted exponents; the inflow baseline ran on M8-fill,
# whose epicentre rows are 28-39% gravity fill (so it overlapped the gravity
# baseline); and "Adjacency-B7" runs the nearest-affected DISTANCE branch because
# no contiguity matrix exists anywhere in this pipeline.

test_that("no baseline scorer touches the renewal machinery", {
  st <- file.path(here::here(), "spatiotemporal")
  src <- paste(readLines(file.path(st, "05_baseline_models.R"), warn = FALSE), collapse = "\n")
  code <- paste(sub("#.*$", "", readLines(file.path(st, "05_baseline_models.R"), warn = FALSE)),
                collapse = "\n")
  # Generation-time weighting, force-of-infection, the fitted import coefficient and R(t) are
  # the renewal model's machinery. A structural baseline must use none of it.
  for (tok in c("gt_pmfs", "compute_foi", "\\.gweighted_own", "make_gt_pmf",
                "weekly_censored_gt", "bayes_rt_week_draws", "epinow", "predict_bayes"))
    expect_false(grepl(tok, code), info = paste("05_baseline_models.R references", tok))
})

test_that("the epicentre scorers are renewal-free and share one origin resolution", {
  for (fn in c("naive_epicentre_inflow_scores", "epicentre_travel_time_scores",
               ".resolve_epicentre_origins"))
    if (!exists(fn, mode = "function")) testthat::skip(paste0(fn, "() unavailable"))
  st <- file.path(here::here(), "spatiotemporal")
  code <- paste(sub("#.*$", "", readLines(file.path(st, "20_forecast_detail.R"), warn = FALSE)),
                collapse = "\n")
  # Both scorers must resolve the epicentre through the SHARED helper — if one harmonised
  # zone spellings and the other did not, the two baselines would not share an origin set and
  # would stop being comparable.
  for (fn in c("naive_epicentre_inflow_scores", "epicentre_travel_time_scores")) {
    body_txt <- paste(deparse(body(get(fn))), collapse = "\n")
    expect_true(grepl(".resolve_epicentre_origins", body_txt, fixed = TRUE),
                info = paste(fn, "does not use the shared epicentre resolver"))
    # and must read no incidence at all
    for (tok in c("gt_pmfs", "compute_foi", "confirmed_nc", "p_recal"))
      expect_false(grepl(tok, body_txt, fixed = TRUE),
                   info = paste(fn, "references", tok))
  }
})

test_that("the travel-time score is a decreasing function of travel time, epicentre zeroed", {
  if (!exists("epicentre_travel_time_scores", mode = "function"))
    testthat::skip("epicentre_travel_time_scores() unavailable")
  z <- c("E1", "A", "B", "C")
  # travel time in minutes from E1: A=10, B=60, C=Inf (unroutable)
  m <- matrix(c(0, 10, 60, Inf,
                10, 0, 20, Inf,
                60, 20, 0, Inf,
                Inf, Inf, Inf, 0), nrow = 4, byrow = TRUE, dimnames = list(z, z))
  s <- epicentre_travel_time_scores(m, "E1", z)
  expect_equal(unname(s[["E1"]]), 0)          # origin is never an at-risk target
  expect_gt(s[["A"]], s[["B"]])               # nearer = higher score
  expect_equal(unname(s[["A"]]), 1 / 11)
  expect_equal(unname(s[["B"]]), 1 / 61)
  expect_equal(unname(s[["C"]]), 0)           # unroutable = no connectivity, not NA
  expect_true(all(is.finite(s)))
})

test_that("injected static baselines are declared rank-only", {
  if (!exists("append_naive_detection_curve_model", mode = "function"))
    testthat::skip("append_naive_detection_curve_model() unavailable")
  lfo <- data.frame(method = "M", horizon = 1L, fold_id = 1L,
                    health_zone = c("A", "B", "C"), p_invasion = c(.1, .2, .3),
                    is_new_invasion = c(0L, 1L, 0L), was_active_before = FALSE,
                    prob_calibrated = TRUE, stringsAsFactors = FALSE)
  out <- append_naive_detection_curve_model(lfo, c(A = 3, B = 2, C = 1), "Baseline-x")
  inj <- out[out$method == "Baseline-x", ]
  expect_equal(nrow(inj), 3L)
  # A connectivity score is not a probability: proper scores must be suppressed for it, or
  # evaluate_invasion() would compute a log score on a mobility inflow.
  expect_true(all(inj$prob_calibrated == FALSE))
  expect_equal(inj$p_invasion[match(c("A", "B", "C"), inj$health_zone)], c(3, 2, 1))
  # idempotent
  expect_equal(nrow(append_naive_detection_curve_model(out, c(A = 9), "Baseline-x")), nrow(out))
})

test_that("an injected baseline is hung on the WIDEST method, not on row order", {
  if (!exists("append_naive_detection_curve_model", mode = "function"))
    testthat::skip("append_naive_detection_curve_model() unavailable")
  # `narrow` is listed FIRST and covers one fold; `wide` covers two. The baseline must inherit
  # `wide`'s two folds. Under the old `unique(method)[1]` rule it inherited `narrow`'s single
  # fold, so Panel 2A's structural nulls were drawn on half the cells of the model they exist
  # to beat -- and which method sorted first was an accident of row order.
  mk <- function(m, folds) data.frame(
    method = m, horizon = 1L,
    fold_id = rep(folds, each = 3L),
    health_zone = rep(c("A", "B", "C"), times = length(folds)),
    p_invasion = 0.1, is_new_invasion = 0L, was_active_before = FALSE,
    prob_calibrated = TRUE, stringsAsFactors = FALSE)
  lfo <- rbind(mk("narrow", 1L), mk("wide", 1:2))
  inj <- append_naive_detection_curve_model(lfo, c(A = 3, B = 2, C = 1), "Baseline-x")
  inj <- inj[inj$method == "Baseline-x", ]
  expect_equal(nrow(inj), 6L)
  expect_setequal(unique(inj$fold_id), 1:2)
  # Deterministic: reversing the row order must not change the skeleton.
  inj2 <- append_naive_detection_curve_model(lfo[rev(seq_len(nrow(lfo))), ],
                                             c(A = 3, B = 2, C = 1), "Baseline-x")
  inj2 <- inj2[inj2$method == "Baseline-x", ]
  expect_equal(nrow(inj2), 6L)
  expect_setequal(unique(inj2$fold_id), 1:2)
})

test_that("both Figure 2 scripts draw the comparators from ONE definition", {
  # THE INVARIANT IS THE SOURCE, NOT A LIST OF STRINGS. Figure 2 (make_manuscript_figures.R) and
  # its cascade analogue (make_manuscript_figure2_cascade.R) are the same prioritisation panel
  # over two horizons, so a reader comparing them must see one comparator set, one name and one
  # colour per comparator. Both take the MEMBERSHIP from FIGURE_BASELINE_METHODS through
  # fs_baseline_methods(), and the labels and colours from BASELINE_LBL / BASELINE_COL through
  # fs_panelA_scales() (figure_style.R). This test used to grep both scripts for three typed
  # labels — "Gravity model (fitted flows)" among them — which neither script contains any more,
  # because neither types a label: it failed on code that was correct, and would have passed on
  # code that typed the right three strings into a private palette.
  st <- file.path(here::here(), "spatiotemporal")
  scripts <- c(fig2 = file.path(st, "make_manuscript_figures.R"),
               cascade = file.path(st, "make_manuscript_figure2_cascade.R"))
  if (!all(file.exists(scripts))) testthat::skip("figure scripts not present")
  # Live code only: the explanatory comments deliberately name what the code no longer does.
  code <- lapply(scripts, function(f) paste(sub("#.*$", "", readLines(f, warn = FALSE)),
                                            collapse = "\n"))
  for (nm in names(code)) {
    expect_true(grepl("fs_baseline_methods()", code[[nm]], fixed = TRUE),
                info = paste(nm, "does not take its comparator set from fs_baseline_methods()"))
    expect_true(grepl("fs_panelA_scales(", code[[nm]], fixed = TRUE),
                info = paste(nm, "does not take its labels/colours from fs_panelA_scales()"))
    expect_false(grepl("BASELINE_(LBL|COL)(_PLAIN)?\\s*<-", code[[nm]]),
                 info = paste(nm, "defines its own BASELINE_LBL/BASELINE_COL"))
  }
  # And the shared helper returns exactly the house-style labels and colours for the published
  # set, featured series first.
  fs <- new.env()
  sys.source(file.path(st, "figure_style.R"), envir = fs)
  meth <- get0("FIGURE_BASELINE_METHODS", ifnotfound = fs$fs_baseline_methods())
  sc <- fs$fs_panelA_scales("FEATURED_KEY", "Featured", meth)
  expect_identical(names(sc$values), c("FEATURED_KEY", meth))
  expect_identical(unname(sc$values[meth]), unname(fs$BASELINE_COL[meth]))
  expect_identical(unname(sc$labels[meth]), unname(fs$BASELINE_LBL[meth]))
  expect_identical(unname(sc$values["FEATURED_KEY"]), fs$FEATURED_COL)
  # The plain-language variant changes only the words.
  scp <- fs$fs_panelA_scales("FEATURED_KEY", "Featured", meth, plain = TRUE)
  expect_identical(scp$values, sc$values)
  expect_identical(unname(scp$labels[meth]), unname(fs$BASELINE_LBL_PLAIN[meth]))
  # The retired set must be gone from the LIVE CODE of the cascade figure.
  for (dead in c("Proximity to epicentre", "Epicentre mobility inflow"))
    expect_false(grepl(dead, code[["cascade"]], fixed = TRUE),
                 info = paste("Figure2_cascade still references the retired baseline:", dead))
})
