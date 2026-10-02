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
#   4. the arrival origin is one date for the whole outbreak.
# They read the published CSVs and skip when the pipeline has not been run.
# =============================================================================

.arr_paths <- function() {
  d <- file.path(here::here(), "spatiotemporal", "outputs", "key_outputs")
  list(dat = file.path(d, "arrival_predictors.csv"),
       fit = file.path(d, "arrival_predictor_fits.csv"))
}
.skip_unless_arr <- function() {
  p <- .arr_paths()
  if (!file.exists(p$dat) || !file.exists(p$fit))
    testthat::skip("arrival_predictor*.csv not present — run 43_spread_kinematics.R")
  lapply(p, function(f) suppressWarnings(readr::read_csv(f, show_col_types = FALSE)))
}

test_that("every published arrival-predictor statistic is reproducible from the per-zone table", {
  a <- .skip_unless_arr(); dat <- a$dat; fit <- a$fit
  expect_true(nrow(fit) > 0L)
  for (i in seq_len(nrow(fit))) {
    p <- fit$predictor[i]
    expect_true(p %in% names(dat), info = paste("predictor column missing:", p))
    y <- dat[[p]]; x <- dat$arrival_days
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
  # The mobility predictors are -log(share) and are NA wherever no flow out of the
  # epicentre was released for the zone. That is selection ON the predictor and it is
  # correlated with arrival, so a single table-wide n would misstate every facet but
  # the geographic ones. If this ever becomes a single common n, the panel's "n ="
  # annotation has stopped carrying information and the guard should be revisited,
  # not deleted.
  expect_true(all(fit$n <= nrow(dat)))
  expect_true(any(fit$n == nrow(dat)))          # geography is complete
  expect_gt(dplyr::n_distinct(fit$n), 1L)       # mobility is not
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

test_that("both Figure 2 scripts name the same three baselines", {
  st <- file.path(here::here(), "spatiotemporal")
  f2 <- file.path(st, "make_manuscript_figures.R")
  fc <- file.path(st, "make_manuscript_figure2_cascade.R")
  if (!file.exists(f2) || !file.exists(fc)) testthat::skip("figure scripts not present")
  a <- paste(readLines(f2, warn = FALSE), collapse = "\n")
  b <- paste(readLines(fc, warn = FALSE), collapse = "\n")
  # The rendered LABELS must match between the two figures, or a reader comparing them sees
  # two different baseline sets.
  for (lbl in c("Gravity model \\(fitted flows\\)",
                "Flowminder cohort inflow from epicentre",
                "Road travel time from epicentre")) {
    expect_true(grepl(lbl, a), info = paste("manuscript Figure 2 is missing label:", lbl))
    expect_true(grepl(lbl, b), info = paste("Figure2_cascade is missing label:", lbl))
  }
  # And the retired set must be gone from the LIVE CODE of both. Comments are stripped first:
  # the removal notes deliberately name what they replaced, and a guard that forbade the words
  # outright would force those explanations to be deleted too.
  b_code <- paste(sub("#.*$", "", readLines(fc, warn = FALSE)), collapse = "\n")
  for (dead in c("Proximity to epicentre", "Epicentre mobility inflow"))
    expect_false(grepl(dead, b_code, fixed = TRUE),
                 info = paste("Figure2_cascade still references the retired baseline:", dead))
})
