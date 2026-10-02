# =============================================================================
# tests/test_bayes_projection.R — Bayesian h >= 2 source projection and the fill-only grid
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Covers (21_bayesian_renewal.R / 22_daily_reissue.R):
#   * .bayes_project_mu(): the per-draw projection against a hand-computed closed form,
#     R draws never touching h = 1, the logit link, and the missing-projection error;
#   * .pair_rt_draws(): input validation and restoration of the caller's RNG stream;
#   * bayes_forecast_offsets(): refuses h >= 2 without EpiNow2 R draws;
#   * bayes_fill_twin_label() and bayes_default_grid(INCLUDE_UNFILLED_MODELS = FALSE);
#   * linelist_observation_date(): priority order and implausible-date rejection.
# Run:  Rscript spatiotemporal/tests/run_tests.R
#   or: testthat::test_file("spatiotemporal/tests/test_bayes_projection.R")
# =============================================================================

suppressPackageStartupMessages(library(testthat))

.st <- if (requireNamespace("here", quietly = TRUE)) file.path(here::here(), "spatiotemporal") else "."
.need <- c(compute_foi = "06_simple_models.R", epinow2_daily_confirmed = "02_epi_params.R",
           .bayes_project_mu = "21_bayesian_renewal.R",
           linelist_observation_date = "22_daily_reissue.R")
if (!exists("OUT_DIR")) try(suppressMessages(source(file.path(.st, "00_config.R"))), silent = TRUE)
for (.fn in names(.need))
  if (!exists(.fn, mode = "function"))
    try(suppressWarnings(suppressMessages(source(file.path(.st, .need[[.fn]])))), silent = TRUE)
.skip_unless <- function(fn) if (!exists(fn, mode = "function")) skip(paste(fn, "not available"))

# ---- a small synthetic system -----------------------------------------------------------
.zones <- c("A", "B", "C")
.Y <- matrix(c(4, 6, 8,    # A: weeks 1..3
               0, 1, 2,    # B
               1, 0, 3),   # C
             nrow = 3, byrow = TRUE, dimnames = list(.zones, NULL))
.W <- matrix(c(0,   0.7, 0.3,
               0.5, 0,   0.5,
               0.2, 0.8, 0),
             nrow = 3, byrow = TRUE, dimnames = list(.zones, .zones))
.G <- c(0.6, 0.4)
.nd <- function() {
  L1 <- compute_foi(.Y, .W, .G, t_idx = 4L, .zones)
  data.frame(health_zone = rep(.zones, 2), horizon = rep(1:2, each = 3),
             .off = c(log(pmax(L1, 1e-12)), rep(NA_real_, 3)))
}
.proj <- function(rt) list(Y = .Y, W = .W, G = .G, zones = .zones, rt_draws = rt)

test_that(".bayes_project_mu reproduces the closed-form projection draw by draw", {
  .skip_unless(".bayes_project_mu"); .skip_unless("compute_foi")
  b <- 0.08; R <- c(0.5, 1, 1.5, 2); S <- length(R)
  nd <- .nd()
  mu <- .bayes_project_mu(matrix(log(b), S, nrow(nd)), nd, .proj(R), "cloglog", 2L)
  # The pairing is a seeded random permutation of R: recover it the same way.
  Rp <- .pair_rt_draws(R, S, get0("RANDOM_SEED", ifnotfound = 20260704L))
  src1 <- .G[1] * .Y[, 3] + .G[2] * .Y[, 2]
  L1   <- as.numeric(t(.W) %*% src1)
  for (s in seq_len(S)) {
    Y4   <- Rp[s] * src1 + b * L1                      # R * own + that draw's introductions
    src2 <- .G[1] * Y4 + .G[2] * .Y[, 3]
    L2   <- as.numeric(t(.W) %*% src2)
    expect_equal(mu[s, 1:3], b * L1, tolerance = 1e-12)
    expect_equal(mu[s, 4:6], b * L2, tolerance = 1e-12)
  }
  expect_setequal(Rp, R)
})

test_that("R draws never change h = 1, and raise h = 2 monotonically", {
  .skip_unless(".bayes_project_mu")
  nd <- .nd(); eta0 <- matrix(log(0.05), 5, nrow(nd))
  lo <- .bayes_project_mu(eta0, nd, .proj(rep(0.5, 5)), "cloglog", 2L)
  hi <- .bayes_project_mu(eta0, nd, .proj(rep(2.0, 5)), "cloglog", 2L)
  expect_identical(lo[, 1:3], hi[, 1:3])
  expect_true(all(hi[, 4:6] > lo[, 4:6]))
})

test_that("a non-finite kernel cell stays local at h = 2, as in compute_foi", {
  .skip_unless(".bayes_project_mu"); .skip_unless("compute_foi")
  Wna <- .W; Wna["A", "B"] <- NA
  L1 <- compute_foi(.Y, Wna, .G, t_idx = 4L, .zones)
  nd <- data.frame(health_zone = rep(.zones, 2), horizon = rep(1:2, each = 3),
                   .off = c(log(pmax(L1, 1e-12)), rep(NA_real_, 3)))
  pr <- list(Y = .Y, W = Wna, G = .G, zones = .zones, rt_draws = rep(1, 2))
  mu <- .bayes_project_mu(matrix(log(0.05), 2, nrow(nd)), nd, pr, "cloglog", 2L)
  expect_true(all(is.na(mu[, 2])) && all(is.finite(mu[, c(1, 3)])))       # h = 1: B only
  expect_true(all(is.na(mu[, 5])) && all(is.finite(mu[, c(4, 6)])))       # h = 2: still B only
})

test_that("the logit link is reconciled on the per-week hazard scale", {
  .skip_unless(".bayes_project_mu")
  nd <- .nd(); eta0 <- matrix(-4, 2, nrow(nd))
  mu <- .bayes_project_mu(eta0, nd, NULL, "logit", 1L)
  expect_equal(mu[1, 1:3], -log1p(-stats::plogis(-4 + nd$.off[1:3])), tolerance = 1e-12)
  expect_true(all(is.na(mu[, 4:6])))
})

test_that("h >= 2 without a projection or without R draws is an error", {
  .skip_unless(".bayes_project_mu"); .skip_unless("bayes_forecast_offsets")
  nd <- .nd(); eta0 <- matrix(log(0.05), 3, nrow(nd))
  expect_error(.bayes_project_mu(eta0, nd, NULL, "cloglog", 2L), "projection")
  expect_error(.bayes_project_mu(eta0, nd, .proj(numeric(0)), "cloglog", 2L), "EpiNow2")
  expect_error(bayes_forecast_offsets(NULL, NULL, NULL, NULL, NULL, .zones, "M8", "medium",
                                      character(0), c(1L, 2L)), "EpiNow2")
  expect_error(bayes_forecast_offsets(NULL, NULL, NULL, NULL, NULL, .zones, "M8", "medium",
                                      character(0), c(1L, 2L), rt_draws = c(1, NA)), "EpiNow2")
})

test_that(".pair_rt_draws validates input and restores the caller's RNG stream", {
  .skip_unless(".pair_rt_draws")
  set.seed(99); before <- .Random.seed
  x <- .pair_rt_draws(c(1.1, 1.2, 1.3, 1.4), 3L, seed = 7L)
  expect_identical(.Random.seed, before)
  expect_length(x, 3L); expect_true(!anyDuplicated(x))          # without replacement when S <= n
  expect_length(.pair_rt_draws(c(1, 2), 5L, seed = 7L), 5L)     # with replacement when S > n
  expect_error(.pair_rt_draws(c(1, -0.1), 2L, 1L))
  expect_error(.pair_rt_draws(c(1, Inf), 2L, 1L))
  expect_error(.pair_rt_draws(numeric(0), 2L, 1L))
})

test_that(".rt_prior_draws samples EPINOW2_R_PRIOR on the natural scale and restores the RNG", {
  .skip_unless(".rt_prior_draws")
  skip_if(!exists("EPINOW2_R_PRIOR"), "EPINOW2_R_PRIOR not defined")
  set.seed(5); before <- .Random.seed
  x <- .rt_prior_draws(200000L, seed = 11L)
  expect_identical(.Random.seed, before)
  expect_equal(mean(x), EPINOW2_R_PRIOR$mean, tolerance = 0.01)
  expect_equal(stats::sd(x), EPINOW2_R_PRIOR$sd, tolerance = 0.02)
  expect_identical(x[1:10], .rt_prior_draws(200000L, seed = 11L)[1:10])
  if (requireNamespace("EpiNow2", quietly = TRUE)) {   # same conversion as EpiNow2::LogNormal
    p <- unclass(EpiNow2::LogNormal(mean = EPINOW2_R_PRIOR$mean, sd = EPINOW2_R_PRIOR$sd))$parameters
    sdlog <- sqrt(log1p((EPINOW2_R_PRIOR$sd / EPINOW2_R_PRIOR$mean)^2))
    expect_equal(p$sdlog, sdlog, tolerance = 1e-8)
    expect_equal(p$meanlog, log(EPINOW2_R_PRIOR$mean) - sdlog^2 / 2, tolerance = 1e-8)
  }
})

# Temporarily set global objects (config flags, or a mocked function) for one expression.
.with_flags <- function(flags, code) {
  old <- lapply(names(flags), function(n) if (exists(n, envir = .GlobalEnv)) get(n, envir = .GlobalEnv) else NULL)
  names(old) <- names(flags)
  for (n in names(flags)) assign(n, flags[[n]], envir = .GlobalEnv)
  on.exit(for (n in names(old)) {
    if (is.null(old[[n]])) rm(list = n, envir = .GlobalEnv) else assign(n, old[[n]], envir = .GlobalEnv)
  }, add = TRUE)
  force(code)
}

.synthetic_ll <- function(n_confirmed) {
  lo <- as.Date(get0("OUTBREAK_START", ifnotfound = as.Date("2026-04-30")))
  d  <- lo + seq_len(max(n_confirmed, 1L))
  data.frame(confirmed = seq_along(d) <= n_confirmed,
             date_of_symptom_onset = d, date_of_sample_collection = d + 2L,
             lab_analysis_date = as.Date(NA), date_of_notification = as.Date(NA),
             date_index = d)
}

test_that("epinow2_daily_confirmed dates by date_index, censors as-of, and is gapless", {
  .skip_unless("epinow2_daily_confirmed"); .skip_unless("linelist_observation_date")
  lo <- as.Date(get0("OUTBREAK_START", ifnotfound = as.Date("2026-04-30")))
  ll <- data.frame(confirmed = c(TRUE, TRUE, TRUE, FALSE, TRUE),
                   date_of_symptom_onset = as.Date(rep(NA, 5)),
                   date_of_sample_collection = lo + c(20L, 12L, 30L, 5L, 9L),
                   lab_analysis_date = as.Date(NA), date_of_notification = as.Date(NA),
                   date_index = lo + c(13L, 5L, 23L, 1L, NA))
  s <- epinow2_daily_confirmed(ll, lo + 25, caller = "test")
  # onset-less cases counted on date_index, never on their later sample date; the unconfirmed
  # record and the record without a date_index are excluded; the series runs to end_date
  expect_equal(s$date[1], lo + 5); expect_equal(max(s$date), lo + 25)   # Date value, not storage type
  expect_identical(nrow(s), 21L); expect_identical(sum(s$confirm), 3L)
  expect_identical(s$confirm[match(lo + c(5L, 13L, 23L), s$date)], c(1L, 1L, 1L))
  # as of day 25 the record sampled on day 30 is not yet observable
  a <- suppressMessages(epinow2_daily_confirmed(ll, lo + 25, issue_date = lo + 25, caller = "test"))
  expect_identical(sum(a$confirm), 2L)
  expect_identical(nrow(epinow2_daily_confirmed(ll, lo + 3, caller = "test")), 0L)
  expect_error(epinow2_daily_confirmed(ll[, setdiff(names(ll), "date_index")], lo + 25), "date_index")
})

test_that(".rt_recent_posterior never carries a posterior with a different or missing specification", {
  .skip_unless(".rt_recent_posterior")
  tmp <- tempfile("carry_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  gt <- c(15.3, 9.3, 90); ws <- as.Date("2026-08-11")
  mk <- function(spec, wk, vals, tag) {
    d <- vals
    attr(d, "gt") <- c(mean = gt[1], sd = gt[2], max = gt[3])
    attr(d, "end_date") <- wk + 6; attr(d, "source") <- "epinow2"
    if (!is.null(spec)) attr(d, "rt_spec") <- spec
    saveRDS(d, file.path(tmp, sprintf("epinow2_rtweek_%s_%s.rds", format(wk, "%Y%m%d"), tag)))
  }
  mk(NULL,  ws - 7, c(1.5, 1.6), "a")   # unstamped (pre-v5 dating): must be ignored
  mk("OLD", ws - 7, c(1.4, 1.4), "b")   # different specification: must be ignored
  expect_null(.rt_recent_posterior(tmp, ws, gt, spec = "NEW"))
  mk("NEW", ws - 14, c(0.8, 0.9), "c")
  got <- .rt_recent_posterior(tmp, ws, gt, spec = "NEW")
  expect_equal(as.numeric(got), c(0.8, 0.9))
  expect_identical(attr(got, "carried_gap_weeks"), 2L)
  expect_identical(attr(got, "source"), "carried-forward")
  expect_error(.rt_recent_posterior(tmp, ws, gt), "spec")
})

test_that("bayes_rt_week_draws always returns an R: prior when no case, prior after two failed fits", {
  .skip_unless("bayes_rt_week_draws")
  skip_if(!isTRUE(get0(".HAVE_EPINOW2", ifnotfound = FALSE)), "EpiNow2 not installed")
  skip_if(!exists("effective_onset_sample_delay", mode = "function"), "delay spec not available")
  pmf <- structure(rep(1 / 90, 90), gt_mean = 15.3, gt_sd = 9.3)
  lo  <- as.Date(get0("OUTBREAK_START", ifnotfound = as.Date("2026-04-30")))
  calls <- 0L
  mock <- function(...) { calls <<- calls + 1L; stop("MOCK_FIT_FAILURE") }
  tmp <- tempfile("rtdraws_"); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  # (a) no observable confirmed case -> prior, without attempting a fit
  r0 <- .with_flags(list(.epinow2_rt_fit = mock), expect_warning(
    bayes_rt_week_draws(.synthetic_ll(0L), pmf, lo + 21, lo + 28, cache_dir = tmp,
                        n_prior_draws = 500L), "prior"))
  expect_identical(attr(r0, "source"), "prior"); expect_length(r0, 500L)
  expect_true(all(is.finite(r0) & r0 > 0)); expect_identical(calls, 0L)

  # (b) cases observable but the fit fails -> retried once, then prior; nothing cached
  r1 <- .with_flags(list(.epinow2_rt_fit = mock), expect_warning(
    bayes_rt_week_draws(.synthetic_ll(20L), pmf, lo + 14, lo + 21, cache_dir = tmp,
                        n_prior_draws = 500L), "failed twice"))
  expect_identical(attr(r1, "source"), "prior"); expect_identical(calls, 2L)
  expect_length(list.files(tmp, "\\.rds$"), 0L)
  expect_length(list.files(tmp, "\\.lock$", all.files = TRUE), 0L)
  log <- utils::read.csv(file.path(tmp, "rt_draws_source_log.csv"))
  expect_identical(log$source, c("prior", "prior"))
  # 20 confirmed onsets on days 1..20 with samples 2 days later: by the issue date (day 21) only
  # onsets up to day 19 have been sampled, so 19 cases are observable.
  expect_identical(log$n_cases, c(0L, 19L))
})

test_that("bayes_fill_twin_label maps exactly the composites that have a fill twin", {
  .skip_unless("bayes_fill_twin_label")
  lab <- c("Bayes-M8-med", "Bayes-M14-geo", "Bayes-M17-tvweek", "Bayes-M4-med",
           "Bayes-M14-dist", "Bayes-M17-dist-med", "Bayes-M14-fill-med", "Bayes-ens-mean",
           "Bayes-M10-med", "Bayes-M13-split-med", "Bayes-M14-dist-split-geo",
           "Bayes-M8-full-susp", "Bayes-M16-med")
  expect_identical(bayes_fill_twin_label(lab),
                   c("Bayes-M8-fill-med", "Bayes-M14-fill-geo", "Bayes-M17-fill-tvweek", NA,
                     # -dist models DO have fill twins (M8/M10/M13/M14/M17-dist-fill).
                     "Bayes-M14-dist-fill", "Bayes-M17-dist-fill-med",
                     NA,                    # already filled
                     NA,                    # ensemble: no single kernel
                     "Bayes-M10-fill-med",  # M10-fill exists (fill is the default)
                     NA, NA,                # -split models are filled already
                     "Bayes-M8-fill-full-susp",   # multi-token suffix preserved
                     "Bayes-M16-fill-med"))
  # M4 has no fill twin in MOBILITY_FILL_IDS, and an unknown kernel is never mapped.
  expect_identical(bayes_fill_twin_label(c("Bayes-M4-dist-geo", "Bayes-M99-med", NA)),
                   c(NA_character_, NA_character_, NA_character_))
})

.grid_flags <- list(INCLUDE_SOURCEFILL_MODELS = TRUE, MOBILITY_SOURCE_FILL = "unmeasured",
                    INCLUDE_OSRM_DIST_MODELS = TRUE, INCLUDE_GEO_COV_MODELS = TRUE,
                    INCLUDE_FULL_COV_MODELS = FALSE, INCLUDE_FLOWSTATIC_MODELS = TRUE,
                    INCLUDE_TV_BETA_MODELS = TRUE, INCLUDE_SUSPECTED_COV_MODELS = TRUE,
                    INCLUDE_COHORT_MODELS = TRUE, INCLUDE_M11_MODELS = FALSE,
                    INCLUDE_LOGIT_SENS_MODELS = FALSE, INCLUDE_M9_MODELS = FALSE,
                    INCLUDE_M15_MODELS = FALSE, INCLUDE_COHORT_SPLIT_MODELS = TRUE,
                    # PINNED so these tests exercise the unfilled-drop logic rather than
                    # whichever split families 00_config.R happens to select. Set inside
                    # list() on purpose: c(list(...), X = c("a","b")) flattens the vector
                    # into X1/X2 and the flag would never be seen.
                    SPLIT_FAMILIES = c("M8", "M10", "M13", "M14", "M16", "M17"),
                    # THE TIME-VARYING AND GENERATION-TIME ARMS RIDE ONE BASE MODEL, and by
                    # default that model is read from outputs/key_outputs/model_selection.json
                    # — a file the pipeline REWRITES on every run. Left unpinned, this fixture
                    # would compose a different grid depending on which model the last run
                    # happened to select, and these tests would pass or fail for reasons that
                    # have nothing to do with the code under test. Pin the base and the
                    # process list.
                    TV_BETA_BASE_MODEL = "Bayes-M14-fill-med",
                    INCLUDE_GT_SENSITIVITY_MODELS = TRUE,
                    TV_BETA_TYPES = c("trend", "week", "rw1", "ar1", "gp"),
                    GT_SENSITIVITY_PROFILES = c("short", "long"))
# The kernel set a DEFAULT build produces (03_mobility_matrices.R): every fill twin,
# including M10-fill and the -dist-fill family, plus the origin-split composites. An
# incomplete fixture would let unfilled kernels survive the filter only because their
# twin was missing from the test, which is the opposite of what these tests check.
.all_mobs <- c("M4", "M8", "M10", "M13", "M14", "M16", "M17",
               "M4-dist", "M8-dist", "M10-dist", "M13-dist", "M14-dist", "M17-dist",
               "M8-fill", "M10-fill", "M13-fill", "M14-fill", "M16-fill", "M17-fill",
               "M8-dist-fill", "M10-dist-fill", "M13-dist-fill", "M14-dist-fill",
               "M17-dist-fill",
               "M13-split", "M14-split", "M16-split", "M17-split",
               "M13-dist-split", "M14-dist-split", "M17-dist-split")

test_that("INCLUDE_UNFILLED_MODELS = FALSE drops exactly the kernels with a built fill twin", {
  .skip_unless("bayes_default_grid")
  mm <- stats::setNames(as.list(seq_along(.all_mobs)), .all_mobs)
  full <- .with_flags(c(.grid_flags, INCLUDE_UNFILLED_MODELS = TRUE), bayes_default_grid(mm))
  only <- .with_flags(c(.grid_flags, INCLUDE_UNFILLED_MODELS = FALSE), bayes_default_grid(mm))
  mob_full <- vapply(full, `[[`, "", "mob"); mob_only <- vapply(only, `[[`, "", "mob")
  # Every kernel whose "-fill" twin is in the fixture, travel-time AND road-distance.
  parents <- c("M8", "M10", "M13", "M14", "M16", "M17",
               "M8-dist", "M10-dist", "M13-dist", "M14-dist", "M17-dist")
  expect_true(all(parents %in% mob_full))
  expect_false(any(parents %in% mob_only))
  # everything that is not on a parent kernel survives, with identical specs
  keep <- full[!(mob_full %in% parents)]
  expect_identical(only, keep)
  lab_only <- vapply(only, `[[`, "", "label")
  # Filled twins, the split family and kernels with no twin (M4) all survive.
  # The time-varying and generation-time arms are variants of the ONE pinned base model
  # (TV_BETA_BASE_MODEL above), so their labels extend that model's label — they are no
  # longer one pair per kernel.
  expect_true(all(c("Bayes-M14-fill-med", "Bayes-M14-fill-med-tvweek",
                    "Bayes-M14-fill-med-tvrw1", "Bayes-M14-fill-gtshort", "Bayes-M4-med",
                    "Bayes-M4-dist", "Bayes-M14-dist-fill", "Bayes-M17-dist-fill-med",
                    "Bayes-M13-split-med") %in% lab_only))
  # Unfilled parents are gone, including the road-distance ones.
  expect_false(any(c("Bayes-M14-med", "Bayes-M8-susp", "Bayes-M14-tvweek",
                     "Bayes-M14-dist", "Bayes-M17-dist-med") %in% lab_only))
})

test_that("the time-varying and generation-time arms are built on ONE base model", {
  .skip_unless("bayes_default_grid")
  mm <- stats::setNames(as.list(seq_along(.all_mobs)), .all_mobs)
  g <- .with_flags(c(.grid_flags, INCLUDE_UNFILLED_MODELS = FALSE), bayes_default_grid(mm))
  lab <- vapply(g, `[[`, "", "label")
  tv  <- vapply(g, function(x) x$tv %||% "none", "")
  sen <- vapply(g, function(x) isTRUE(x$sensitivity), logical(1))

  # Exactly five time-varying specs, one per process, all on the pinned base model.
  expect_equal(sort(tv[tv != "none"]), sort(c("trend", "week", "rw1", "ar1", "gp")))
  expect_true(all(vapply(g[tv != "none"], function(x) identical(x$mob, "M14-fill"), logical(1))))
  expect_true(all(startsWith(lab[tv != "none"], "Bayes-M14-fill-med-tv")))

  # tv-trend carries the week_idx covariate (it IS the trend term) and rides the fixed-beta
  # fitter; the structural processes carry the base model's covariates unchanged.
  .sp <- function(l) g[[which(lab == l)]]
  expect_true("week_idx" %in% .sp("Bayes-M14-fill-med-tvtrend")$cov)
  for (ty in c("week", "rw1", "ar1", "gp"))
    expect_false("week_idx" %in% .sp(paste0("Bayes-M14-fill-med-tv", ty))$cov)

  # Two generation-time arms, on the base kernel, at short and long, MARKED as sensitivity so
  # best_invasion_model() cannot feature one — selecting among them would make the generation
  # time a fitted axis.
  expect_equal(sort(lab[sen]), c("Bayes-M14-fill-gtlong", "Bayes-M14-fill-gtshort"))
  expect_equal(sort(vapply(g[sen], `[[`, "", "gt")), c("long", "short"))
  expect_true(all(vapply(g[sen], function(x) identical(x$mob, "M14-fill"), logical(1))))
  # ... and nothing else is marked.
  expect_false(any(sen[tv != "none"]))

  # Every label the grid can produce must still parse back to its kernel, or the cascade and
  # the over-folds refit silently fall back to a different matrix.
  expect_false(any(is.na(vapply(lab, mobility_kernel_from_method, ""))))
})

test_that("with no built fill twin, INCLUDE_UNFILLED_MODELS = FALSE drops nothing and warns", {
  .skip_unless("bayes_default_grid")
  mm <- stats::setNames(as.list(seq_along(.all_mobs)), .all_mobs)
  mm <- mm[!grepl("-fill$", names(mm))]
  ref <- .with_flags(c(.grid_flags, INCLUDE_UNFILLED_MODELS = TRUE), bayes_default_grid(mm))
  expect_warning(
    got <- .with_flags(c(.grid_flags, INCLUDE_UNFILLED_MODELS = FALSE), bayes_default_grid(mm)),
    "no model was dropped")
  expect_identical(got, ref)
  # a parent whose twin is NOT built keeps its models even when other twins are
  mm2 <- stats::setNames(as.list(seq_along(.all_mobs)), .all_mobs)
  mm2[["M13-fill"]] <- NULL
  got2 <- .with_flags(c(.grid_flags, INCLUDE_UNFILLED_MODELS = FALSE), bayes_default_grid(mm2))
  mob2 <- vapply(got2, `[[`, "", "mob")
  expect_true("M13" %in% mob2)
  expect_false("M14" %in% mob2)
})

test_that("linelist_observation_date prefers sample date and rejects implausible dates", {
  .skip_unless("linelist_observation_date")
  lo <- as.Date(get0("OUTBREAK_START", ifnotfound = as.Date("2026-04-30")))
  ll <- data.frame(
    date_of_sample_collection = as.Date(c(NA, "1999-01-01", NA)) ,
    lab_analysis_date         = as.Date(c(NA, NA, NA)),
    date_of_notification      = as.Date(c(NA, NA, NA)),
    confirmed = c(TRUE, FALSE, FALSE))
  ll$date_of_sample_collection[1] <- lo + 10
  ll$date_of_notification[2]      <- lo + 20       # the 1999 sample date is implausible
  out <- suppressMessages(suppressWarnings(
    linelist_observation_date(ll, lo + 60, caller = "test")))
  expect_s3_class(out, "Date")
  expect_identical(out, c(lo + 10, lo + 20, as.Date(NA)))
  expect_error(linelist_observation_date(ll, "2026-06-01"), "scalar Date")
})
