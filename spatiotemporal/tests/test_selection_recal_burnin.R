# =============================================================================
# test_selection_recal_burnin.R — the 2026-10-14 changes:
#   * per-date AUC-PR (.auc_pr_by_date) and the selection view (invasion_selection_view)
#   * the ascertainment burn-in (invasion_week_index / invasion_is_burnin)
#   * the single recalibration factor drawn per posterior draw (.suite_delta_for /
#     .scale_invasion_draws)
#   * the cascade calibration function is defined when 33b is sourced
# =============================================================================
library(testthat)
suppressMessages(source(file.path(here::here(), "spatiotemporal", "00_config.R")))
suppressMessages(source(file.path(here::here(), "spatiotemporal", "16_invasion_eval.R")))
suppressMessages(source(file.path(here::here(), "spatiotemporal", "16b_invasion_recalibration.R")))
suppressMessages(source(file.path(here::here(), "spatiotemporal", "21_bayesian_renewal.R")))

# ---- per-date AUC-PR -----------------------------------------------------------------------
test_that("per-date AUC-PR equals plain AUC-PR on a single date", {
  set.seed(1); p <- runif(300); y <- rbinom(300, 1, p / 4)
  r <- .auc_pr_by_date(p, y, rep(1L, 300))
  expect_equal(r[["ap"]], .auc_pr(p, y))
  expect_equal(r[["skill"]], .auc_pr(p, y) / mean(y))
  expect_equal(r[["n_dates"]], 1)
})

test_that("per-date AUC-PR ignores the level of each date; pooled AUC-PR does not", {
  set.seed(2)
  d <- do.call(rbind, lapply(1:6, function(k) {
    p <- runif(200, 0, 0.05); data.frame(date = k, p = p, y = rbinom(200, 1, p))
  }))
  d$y[d$date == 1][1] <- 1                         # every date carries at least one event
  d2 <- d; d2$p[d2$date == 3] <- d2$p[d2$date == 3] * 9   # one date's level moves; its order does not
  expect_equal(.auc_pr_by_date(d$p, d$y, d$date), .auc_pr_by_date(d2$p, d2$y, d2$date))
  expect_false(isTRUE(all.equal(.auc_pr(d$p, d$y), .auc_pr(d2$p, d2$y))))
})

test_that("the chance level of per-date AUC-PR skill is the exact random-ordering expectation", {
  set.seed(4)
  for (cfg in list(c(n = 60, E = 1), c(n = 60, E = 3), c(n = 200, E = 5))) {
    n <- cfg[["n"]]; E <- cfg[["E"]]; y <- c(rep(1L, E), rep(0L, n - E))
    sim <- mean(replicate(20000, .auc_pr(runif(n), y)))        # random ordering, no ties
    r <- .auc_pr_by_date(runif(n), y, rep(1L, n))
    expect_equal(r[["chance"]] * E / n, sim, tolerance = 0.01,
                 info = sprintf("n=%d E=%d", n, E))
    expect_gt(r[["chance"]], 1)                                 # chance is NOT skill 1
  }
})

test_that("dates without an invasion drop out of the per-date average", {
  set.seed(3); p <- runif(400); date <- rep(1:2, each = 200)
  y <- c(rbinom(200, 1, p[1:200] / 3), rep(0L, 200))
  r <- .auc_pr_by_date(p, y, date)
  expect_equal(r[["ap"]], .auc_pr(p[1:200], y[1:200]))
  expect_equal(r[["n_dates"]], 1)
  expect_true(all(is.na(.auc_pr_by_date(p, rep(0L, 400), date)[c("ap", "skill")])))
})

# ---- burn-in index -------------------------------------------------------------------------
test_that("the week index counts from the grid week containing the outbreak start", {
  st <- as.Date("2026-04-30")                     # a Thursday
  # weeks starting Saturday (anchor 6): the week containing 30 Apr starts 25 Apr
  expect_equal(invasion_week_index(as.Date("2026-04-25"), start = st, anchor = 6L), 1L)
  expect_equal(invasion_week_index(as.Date("2026-05-30"), start = st, anchor = 6L), 6L)
  expect_equal(invasion_week_index(as.Date("2026-06-06"), start = st, anchor = 6L), 7L)
  # weeks starting Monday (anchor 1): the week containing 30 Apr starts 27 Apr
  expect_equal(invasion_week_index(as.Date("2026-06-08"), start = st, anchor = 1L), 7L)
})

test_that("burn-in is week index below ASCERTAINMENT_BURNIN_WEEKS, and NA is not burn-in", {
  w0 <- lubridate::floor_date(OUTBREAK_START, "week", week_start = WEEK_ANCHOR)
  cut <- w0 + 7L * (0:9)
  expect_equal(invasion_is_burnin(cut), seq_along(cut) < ASCERTAINMENT_BURNIN_WEEKS)
  expect_false(invasion_is_burnin(as.Date(NA)))
  expect_equal(CASCADE_CALIB_MIN_WEEKS, ASCERTAINMENT_BURNIN_WEEKS)
})

# ---- selection view --------------------------------------------------------------------------
mk_ev <- function() tibble::tribble(
  # A wins every ALL-ROUND axis; B wins every POST-BURN-IN (sel_*) axis.
  ~method, ~horizon, ~auc_pr_skill, ~mean_rank_of_truth, ~log_score, ~log_score_recal,
    ~sel_auc_pr_skill_date, ~sel_mean_rank_of_truth, ~sel_log_score, ~sel_log_score_recal, ~partial_cv,
  "A", 1L, 40, 10, 0.020, 0.019,  20, 30, 0.030, 0.029, FALSE,
  "B", 1L, 30, 20, 0.030, 0.029,  25, 20, 0.020, 0.019, FALSE,
  "A", 2L, 30, 12, 0.040, 0.039,  15, 35, 0.050, 0.049, FALSE,
  "B", 2L, 20, 22, 0.050, 0.049,  18, 25, 0.040, 0.039, FALSE)

test_that("the selection reads the per-date, post-burn-in axes by default and the old ones on request", {
  ev <- mk_ev()
  expect_equal(best_invasion_model(ev), "B")
  expect_equal(best_invasion_model(ev, basis = "pooled_all_dates"), "A")
  tb <- invasion_selection_table(ev)
  expect_equal(tb$method[which.min(tb$composite)], "B")
})

test_that("the view is idempotent and reversible", {
  ev <- mk_ev()
  v1 <- invasion_selection_view(ev, "per_date_post_burnin")
  v2 <- invasion_selection_view(v1, "per_date_post_burnin")
  expect_equal(as.data.frame(v1), as.data.frame(v2), ignore_attr = TRUE)
  back <- invasion_selection_view(v1, "pooled_all_dates")
  for (a in c("auc_pr_skill", "mean_rank_of_truth", "log_score", "log_score_recal"))
    expect_equal(back[[a]], ev[[a]])
  expect_equal(v1$auc_pr_skill, ev$sel_auc_pr_skill_date)
  expect_equal(v1$orig_auc_pr_skill, ev$auc_pr_skill)
})

test_that("a table without sel_* columns falls back to the previous basis, loudly", {
  ev <- mk_ev()[, c("method", "horizon", "auc_pr_skill", "mean_rank_of_truth", "log_score", "partial_cv")]
  expect_warning(v <- invasion_selection_view(ev, "per_date_post_burnin"), "PREVIOUS basis")
  expect_equal(attr(v, "selection_basis"), "pooled_all_dates")
  expect_false("log_score_recal" %in% names(v))     # an absent column stays absent
})

test_that("evaluate_invasion's sel_* columns are the all-round metrics on the post-burn-in rounds", {
  set.seed(7)
  w0 <- lubridate::floor_date(OUTBREAK_START, "week", week_start = WEEK_ANCHOR)
  f <- do.call(rbind, lapply(1:12, function(k) {
    p <- runif(150, 0, 0.06)
    data.frame(method = "Bayes-A", horizon = 1L, fold_id = k, cutoff = w0 + 7L * (k + 1L),
               health_zone = paste0("z", 1:150), p_invasion = p,
               is_new_invasion = rbinom(150, 1, p), was_active_before = FALSE)
  }))
  ev  <- suppressWarnings(evaluate_invasion(f, n_boot = 0L))
  sub <- f[!invasion_is_burnin(f$cutoff), ]
  evs <- suppressWarnings(evaluate_invasion(sub, n_boot = 0L))
  expect_true(ev$sel_burnin_applied)
  expect_equal(ev$sel_n_folds, length(unique(sub$fold_id)))
  expect_equal(ev$sel_mean_rank_of_truth, evs$mean_rank_of_truth)
  expect_equal(ev$sel_log_score, evs$log_score)
  expect_equal(ev$sel_auc_pr_skill_date, evs$auc_pr_skill_date)
})

# ---- the single recalibration factor in the live draws ---------------------------------------
test_that(".suite_delta_for carries the Gamma posterior of the single factor", {
  tb <- data.frame(method = "M", horizon = c(1L, 2L), delta = 0.5, gamma_shape = 30.5, gamma_rate = 61)
  d <- .suite_delta_for(tb, "M", c(1L, 2L))
  expect_equal(as.numeric(d), c(0.5, 0.5))
  expect_equal(attr(d, "gamma_shape"), 30.5); expect_equal(attr(d, "gamma_rate"), 61)
})

test_that("per-draw delta: mean equals delta, nesting holds draw by draw, caller RNG untouched", {
  set.seed(11); S <- 4000L; Z <- 5L
  c1 <- matrix(rexp(S * Z, 20), S, Z); c2 <- c1 + matrix(rexp(S * Z, 25), S, Z)
  drl <- list(h1 = list(zones = paste0("z", 1:Z), horizon = 1L, cum = c1, p = 1 - exp(-c1)),
              h2 = list(zones = paste0("z", 1:Z), horizon = 2L, cum = c2, p = 1 - exp(-c2)))
  d <- c("1" = 0.5, "2" = 0.5); attr(d, "gamma_shape") <- 30.5; attr(d, "gamma_rate") <- 61
  set.seed(99); before <- runif(1); set.seed(99)
  out <- .scale_invasion_draws(drl, d)
  expect_equal(runif(1), before)                                  # stream untouched
  ratio <- out$h1$cum[, 1] / c1[, 1]
  expect_equal(mean(ratio), 0.5, tolerance = 0.02)                 # Gamma mean = delta
  expect_gt(stats::sd(ratio), 0.05)                                # genuinely drawn
  expect_equal(out$h2$cum[, 1] / c2[, 1], ratio)                   # same delta per draw, both horizons
  expect_true(all(out$h2$p >= out$h1$p))                           # nesting preserved
  # without the Gamma attributes it stays a fixed factor
  fx <- .scale_invasion_draws(drl, c("1" = 0.5, "2" = 0.5))
  expect_equal(fx$h1$cum, c1 * 0.5)
})

# ---- cascade ------------------------------------------------------------------------------------
test_that("33b defines cascade_calibration_origins() when %||% already exists", {
  e <- new.env()
  sys.source(file.path(here::here(), "spatiotemporal", "33b_cascade_calibration.R"), envir = e,
             toplevel.env = globalenv())
  expect_true(exists("cascade_calibration_origins", envir = e, inherits = FALSE))
})

# ---- the selection basis switch (2026-10-19) ---------------------------------------------------
test_that("per_date_all_rounds selects on every round: sel_* equal the all-round per-date columns", {
  set.seed(7)
  w0 <- lubridate::floor_date(OUTBREAK_START, "week", week_start = WEEK_ANCHOR)
  f <- do.call(rbind, lapply(1:12, function(k) {
    p <- runif(150, 0, 0.06)
    data.frame(method = "Bayes-A", horizon = 1L, fold_id = k, cutoff = w0 + 7L * (k + 1L),
               health_zone = paste0("z", 1:150), p_invasion = p,
               is_new_invasion = rbinom(150, 1, p), was_active_before = FALSE)
  }))
  old <- INVASION_SELECTION_BASIS
  on.exit(assign("INVASION_SELECTION_BASIS", old, envir = globalenv()), add = TRUE)
  assign("INVASION_SELECTION_BASIS", "per_date_all_rounds", envir = globalenv())
  expect_false(invasion_selection_drops_burnin())
  expect_true(invasion_selection_per_date())
  ev <- suppressWarnings(evaluate_invasion(f, n_boot = 0L))
  expect_false(ev$sel_burnin_applied)
  expect_equal(ev$sel_n_folds, ev$n_folds)
  expect_equal(ev$sel_auc_pr_skill_date, ev$auc_pr_skill_date)
  expect_equal(ev$sel_mean_rank_of_truth, ev$mean_rank_of_truth)
  expect_equal(invasion_selection_axis_label(), "per-date AUC-PR skill, all rounds")
})

test_that("the basis helpers and axis label cover every basis, and an unknown one stops", {
  expect_true(invasion_selection_drops_burnin("per_date_post_burnin"))
  expect_false(invasion_selection_drops_burnin("pooled_all_dates"))
  expect_false(invasion_selection_per_date("pooled_all_dates"))
  expect_equal(invasion_selection_axis_label("per_date_post_burnin"),
               "per-date AUC-PR skill, post-burn-in rounds")
  expect_equal(invasion_selection_axis_label("pooled_all_dates"), "pooled AUC-PR skill, all rounds")
  expect_error(invasion_selection_axis_label("nonsense"), "unknown selection basis")
  expect_setequal(INVASION_SELECTION_BASES,
                  c("per_date_post_burnin", "per_date_all_rounds", "pooled_all_dates"))
})
