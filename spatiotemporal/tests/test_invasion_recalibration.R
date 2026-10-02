# =============================================================================
# test_invasion_recalibration.R — correctness tests for 16b_invasion_recalibration.R
# Every property the pipeline RELIES on is asserted here, not assumed.
# =============================================================================
library(testthat)
suppressMessages(source(file.path(here::here(), "spatiotemporal", "16b_invasion_recalibration.R")))

set.seed(20260704L)

# ---------------------------------------------------------------------------
context_msg <- function(x) message("\n--- ", x, " ---")

# 1. IDENTITY AT delta = 1 -----------------------------------------------------
test_that("delta = 1 is exactly the identity", {
  p <- c(0, 1e-9, 0.001, 0.05, 0.5, 0.9999, 1, NA)
  expect_equal(recalibrate_invasion_p(p, 1), as.numeric(p), tolerance = 0)
})

# 2. MONOTONICITY / RANK INVARIANCE -------------------------------------------
test_that("the transform is strictly rank-preserving for any delta > 0", {
  p <- sort(runif(500, 0, 0.9))
  for (d in c(0.01, 0.25, 0.5, 1, 2, 5)) {
    q <- recalibrate_invasion_p(p, d)
    expect_true(all(diff(q) > 0))                       # strictly increasing
    expect_equal(order(q), order(p))                    # identical ordering
    expect_equal(rank(q, ties.method = "average"),
                 rank(p, ties.method = "average"))
  }
})

test_that("ties are preserved exactly (no tie-breaking is introduced)", {
  p <- c(0.1, 0.1, 0.2, 0.2, 0.3)
  q <- recalibrate_invasion_p(p, 0.4)
  expect_equal(q[1], q[2]); expect_equal(q[3], q[4])
  expect_equal(rank(q, ties.method = "max"), rank(p, ties.method = "max"))
})

# 3. RANGE / EDGE CASES --------------------------------------------------------
test_that("output stays in [0,1] and preserves NA", {
  p <- c(0, 1, NA, 0.5)
  for (d in c(0.01, 1, 5)) {
    q <- recalibrate_invasion_p(p, d)
    expect_true(all(q >= 0 & q <= 1, na.rm = TRUE))
    expect_true(is.na(q[3]))
    expect_equal(q[1], 0)                                # p = 0 stays impossible
  }
})

test_that("non-usable delta leaves probabilities unchanged (with a warning)", {
  p <- c(0.1, 0.2)
  expect_warning(q <- recalibrate_invasion_p(p, NA_real_))
  expect_equal(q, p)
  expect_warning(q0 <- recalibrate_invasion_p(p, 0))
  expect_equal(q0, p)
})

test_that("vectorised delta applies row-wise", {
  p <- c(0.2, 0.2, 0.2)
  d <- c(0.5, 1, 2)
  q <- recalibrate_invasion_p(p, d)
  expect_equal(q[2], 0.2)
  expect_equal(q, 1 - (1 - p)^d)
  # delta multiplies the HAZARD, so a smaller delta means a SMALLER probability —
  # which is the direction that matters here, since every model over-predicts.
  expect_true(q[1] < q[2] && q[2] < q[3])
})

# 4. ML ESTIMATOR RECOVERS A KNOWN delta ---------------------------------------
test_that("fit_invasion_delta recovers the generating delta on simulated data", {
  n <- 300000L
  p_raw <- runif(n, 0, 0.05)                             # rare-event scale
  for (d_true in c(0.3, 1.0, 2.5)) {
    p_true <- 1 - (1 - p_raw)^d_true
    y <- rbinom(n, 1, p_true)
    d_hat <- fit_invasion_delta(p_raw, y)
    expect_equal(d_hat, d_true, tolerance = 0.05)
  }
})

test_that("the likelihood is maximised at the returned delta", {
  n <- 20000L
  p_raw <- runif(n, 0, 0.05); y <- rbinom(n, 1, 1 - (1 - p_raw)^0.4)
  d_hat <- fit_invasion_delta(p_raw, y)
  ll <- function(d) sum(ifelse(y == 1, log(1 - (1 - p_raw)^d), d * log(1 - p_raw)))
  expect_gt(ll(d_hat), ll(d_hat * 0.9))
  expect_gt(ll(d_hat), ll(d_hat * 1.1))
})

# 5. GUARDS --------------------------------------------------------------------
# THERE IS NO EVENT FLOOR. The fit is Jeffreys-penalised, so it exists for every pool with
# positive hazard -- including one with no invasions, where the plain MLE runs to delta -> 0
# and would zero out a whole fold's forecast.
test_that("the penalised fit is finite even with NO events, where plain ML collapses", {
  set.seed(11); p <- runif(3000, 0, 0.02); y <- rep(0L, 3000)
  d_pen <- fit_invasion_delta(p, y)
  # The unpenalised fit has no interior maximum here, so it runs to the bracket floor and
  # warns. That warning IS the diagnostic; expect it rather than letting it leak.
  d_ml  <- expect_warning(fit_invasion_delta(p, y, penalty = FALSE), "edge of the numerical")
  expect_true(is.finite(d_pen) && d_pen > 0)
  expect_equal(d_ml, RECAL_BAND[1], tolerance = 1e-6)   # pinned at the floor, not estimated
  expect_lt(d_ml, d_pen / 100)
  # "Add half an event on the hazard scale": delta = (events + 1/2) / sum(-log(1 - p)).
  expect_equal(d_pen, 0.5 / sum(-log1p(-p)), tolerance = 1e-3)
  # The penalised forecast survives; the unpenalised one is annihilated. Compare the two
  # directly rather than against absolute thresholds, which depend on the bracket floor.
  rc <- function(pp, d) -expm1(d * log1p(-pp))
  expect_gt(max(rc(p, d_pen)) / max(rc(p, d_ml)), 100)
})

test_that("the penalised fit matches its closed form in the rare-event limit", {
  # delta_hat -> (n_events + 1/2) / sum_i -log(1 - p_i); plain ML -> n_events / sum_i ...
  p <- rep(1e-4, 20000); U <- sum(-log1p(-p))
  for (nev in c(1L, 3L, 10L)) {
    y <- c(rep(1L, nev), rep(0L, 20000L - nev))
    expect_equal(fit_invasion_delta(p, y),                  (nev + 0.5) / U, tolerance = 1e-3)
    expect_equal(fit_invasion_delta(p, y, penalty = FALSE),  nev        / U, tolerance = 1e-3)
  }
})

test_that("a few events still give a usable delta (no floor to fall below)", {
  set.seed(12); p <- runif(1000, 0, 0.05); y <- rep(0L, 1000); y[1:3] <- 1L
  d <- fit_invasion_delta(p, y)
  expect_true(is.finite(d) && d > 0)
})

test_that("all-zero predictions are not estimable", {
  # The ONLY unfittable case left: no row carries positive hazard, so delta is unidentified
  # at every value.
  expect_true(is.na(fit_invasion_delta(rep(0, 100), c(rep(1L, 20), rep(0L, 80)))))
})

# 6. PREQUENTIAL CAUSALITY -----------------------------------------------------
mk_slice <- function(n_folds = 9L, n_zone = 200L, horizon = 1L, seed = 1L) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_folds), function(k) {
    p <- runif(n_zone, 0, 0.05)
    data.frame(method = "M", horizon = horizon, fold_id = k,
               cutoff = as.Date("2026-05-09") + 7L * (k - 1L),
               health_zone = paste0("z", seq_len(n_zone)),
               p_invasion = p,
               is_new_invasion = rbinom(n_zone, 1, 1 - (1 - p)^0.4),
               was_active_before = FALSE, stringsAsFactors = FALSE)
  }))
}

test_that("the prequential delta trains only on folds whose outcome week has fully ELAPSED", {
  # The rule is `C' + 7h + 6 <= C`. The +6 matters: a fold's last outcome week STARTS at
  # C' + 7h and ENDS at C' + 7h + 6, so the previous rule (`C' + 7h <= C`) admitted a fold
  # whose final outcome week had only just begun — the h=1 delta at fold k was then fitted
  # partly on the week starting on fold k's own cutoff. With 7-day fold spacing the corrected
  # rule withholds exactly one further fold at each horizon.
  s1 <- mk_slice(horizon = 1L); s2 <- mk_slice(horizon = 2L)
  q1 <- prequential_invasion_delta(s1, horizon = 1L)
  q2 <- prequential_invasion_delta(s2, horizon = 2L)
  expect_equal(q1$n_train_folds, c(0, 0:7))              # h=1: fold k trains on k-2 folds
  expect_equal(q2$n_train_folds, c(0, 0, 0, 1:6))        # h=2: one further fold withheld
})

test_that("only folds with NO closed training fold fall back to delta = 1", {
  # The fallback is now an information constraint, not a threshold: with nothing to fit there
  # is no delta, and every other fold gets one. Under the old event/origin floors this left
  # over half the folds uncorrected.
  s <- mk_slice(); q <- prequential_invasion_delta(s, horizon = 1L)
  expect_equal(q$delta_preq[q$n_train_folds == 0L],
               rep(1, sum(q$n_train_folds == 0L)))
  expect_false(any(q$delta_estimable[q$n_train_folds == 0L]))
  expect_true(all(q$delta_estimable[q$n_train_folds > 0L]))
  expect_equal(sum(q$delta_estimable), sum(q$n_train_folds > 0L))
})

test_that("a training pool restricted to shared cells does not change which folds are reported", {
  # `fit` supplies the training rows; `d` still defines the fold list. Without that split a
  # fold outside the shared support would come back with no delta at all and drag the whole
  # method's recalibrated block out of the evaluation.
  s <- mk_slice(n_folds = 6L)
  sub <- s[s$fold_id != 1L, , drop = FALSE]          # fold 1 withheld from TRAINING only
  q_all <- prequential_invasion_delta(s, horizon = 1L)
  q_sub <- prequential_invasion_delta(s, horizon = 1L, fit = sub)
  expect_equal(q_sub$fold_id, q_all$fold_id)          # same folds reported
  expect_equal(nrow(q_sub), nrow(q_all))
  # fold 3 trains on fold 1 alone, which `sub` removes -> no training rows -> delta = 1
  expect_equal(q_sub$delta_preq[q_sub$fold_id == 3L], 1)
  expect_true(q_all$delta_estimable[q_all$fold_id == 3L])
})

test_that("a missing cutoff column falls back to fold_id spacing with the same causality", {
  s <- mk_slice(); s$cutoff <- NULL
  q1 <- prequential_invasion_delta(s, horizon = 1L)
  q2 <- prequential_invasion_delta(s, horizon = 2L)
  # Same corrected rule as above, reached through the fold_id spacing fallback.
  expect_equal(q1$n_train_folds, c(0, 0:7))
  expect_equal(q2$n_train_folds, c(0, 0, 0, 1:6))
})

# 7. attach_invasion_recalibration ---------------------------------------------
test_that("attach preserves row order, row count and every original column", {
  s <- rbind(mk_slice(horizon = 1L), mk_slice(horizon = 2L, seed = 2L))
  s$tag <- seq_len(nrow(s))
  out <- suppressMessages(attach_invasion_recalibration(s))
  expect_equal(nrow(out), nrow(s))
  expect_equal(out$tag, s$tag)                            # order preserved
  expect_true(all(names(s) %in% names(out)))
  expect_equal(out$p_invasion, s$p_invasion)              # raw column untouched
})

test_that("p_recal equals the documented transform row by row", {
  s <- mk_slice()
  out <- suppressMessages(attach_invasion_recalibration(s))
  expect_equal(out$p_recal, 1 - (1 - out$p_invasion)^out$delta_preq)
})

test_that("non-scorable rows get NA and never enter the fit", {
  s <- mk_slice()
  s$was_active_before[1:50] <- TRUE
  s$p_invasion[51:60] <- NA_real_
  out <- suppressMessages(attach_invasion_recalibration(s))
  expect_true(all(is.na(out$p_recal[1:60])))
  expect_true(all(is.na(out$delta_preq[1:60])))
})

test_that("recalibration leaves every RANKING metric bit-identical", {
  s <- mk_slice()
  out <- suppressMessages(attach_invasion_recalibration(s))
  o <- out[is.finite(out$p_recal), ]
  for (k in unique(o$fold_id)) {
    g <- o[o$fold_id == k, ]
    expect_equal(rank(-g$p_recal, ties.method = "average"),
                 rank(-g$p_invasion, ties.method = "average"))
    expect_equal(rank(-g$p_recal, ties.method = "max"),
                 rank(-g$p_invasion, ties.method = "max"))
  }
})

# 8. POOLED TABLE --------------------------------------------------------------
test_that("invasion_delta_table returns one row per method x horizon with sane fields", {
  s <- rbind(mk_slice(horizon = 1L), mk_slice(horizon = 2L, seed = 2L))
  tb <- suppressMessages(invasion_delta_table(s, n_boot = 60L))
  expect_equal(nrow(tb), 2L)
  expect_true(all(is.finite(tb$delta)))
  expect_true(all(tb$delta_lo <= tb$delta & tb$delta <= tb$delta_hi, na.rm = TRUE))
  expect_true(all(tb$n_events > 0))
  expect_false(any(tb$boundary_hit))
  # the moment estimator is reported alongside but is NOT what `delta` uses
  expect_true(all(is.finite(tb$delta_moment)))
})

message("\n[test_invasion_recalibration] all assertions executed.")

# =============================================================================
# INTEGRATION — the evaluator, the selector and the prediction hook
# =============================================================================
suppressMessages(source(file.path(here::here(), "spatiotemporal", "16_invasion_eval.R")))

mk_eval_frame <- function(n_folds = 9L, n_zone = 150L, seed = 7L) {
  set.seed(seed)
  grid <- expand.grid(method = c("Bayes-A", "Bayes-B"), horizon = c(1L, 2L),
                      stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    m <- grid$method[i]; h <- grid$horizon[i]
    do.call(rbind, lapply(seq_len(n_folds), function(k) {
      p <- runif(n_zone, 0, 0.06) * ifelse(m == "Bayes-A", 1, 1.8)
      p <- pmin(p, 0.95)
      data.frame(method = m, horizon = h, fold_id = k,
                 cutoff = as.Date("2026-05-09") + 7L * (k - 1L),
                 health_zone = paste0("z", seq_len(n_zone)),
                 p_invasion = p,
                 is_new_invasion = rbinom(n_zone, 1, 1 - (1 - p / 2)^1),
                 was_active_before = FALSE, stringsAsFactors = FALSE)
    }))
  }))
}

test_that("evaluate_invasion adds recal columns without altering any raw column", {
  f <- mk_eval_frame()
  ev_before <- suppressWarnings(evaluate_invasion(f, n_boot = 0L))
  ev_after  <- suppressWarnings(evaluate_invasion(
    suppressMessages(attach_invasion_recalibration(f)), n_boot = 0L))
  raw_cols <- c("auc_pr", "auc_roc", "log_score", "brier_skill", "ece",
                "calibration_in_large", "mean_rank_of_truth", "prec_at_5",
                "prec_at_10", "prec_at_15", "n_atrisk", "n_invasions", "base_rate")
  a <- ev_before[order(ev_before$method, ev_before$horizon), raw_cols]
  b <- ev_after[order(ev_after$method, ev_after$horizon), raw_cols]
  expect_equal(as.data.frame(a), as.data.frame(b))
  expect_true(all(c("log_score_recal", "brier_skill_recal", "ece_recal",
                    "calibration_in_large_recal", "recal_delta_varies") %in% names(ev_after)))
})

test_that("recal columns are NA when 16b has not been applied", {
  ev <- suppressWarnings(evaluate_invasion(mk_eval_frame(), n_boot = 0L))
  expect_true(all(is.na(ev$log_score_recal)))
  expect_true(all(is.na(ev$brier_skill_recal)))
})

test_that("a fixed delta leaves EVERY evaluate_invasion metric identical", {
  f <- mk_eval_frame()
  ev_raw <- suppressWarnings(evaluate_invasion(f, n_boot = 0L))
  f2 <- f; f2$p_invasion <- recalibrate_invasion_p(f$p_invasion, 0.37)
  ev_fix <- suppressWarnings(evaluate_invasion(f2, n_boot = 0L))
  for (cl in c("auc_pr", "auc_roc", "mean_rank_of_truth",
               "prec_at_5", "recall_at_5", "hit_at_5",
               "prec_at_10", "prec_at_15")) {
    expect_equal(ev_raw[[cl]], ev_fix[[cl]], tolerance = 1e-12,
                 info = paste("fixed-delta invariance violated for", cl))
  }
})

test_that("selection falls back loudly when the recal column is absent or incomplete", {
  ev <- suppressWarnings(evaluate_invasion(mk_eval_frame(), n_boot = 0L))
  ev$log_score_recal <- NULL
  expect_warning(m1 <- best_invasion_model(ev, select_on_recal = TRUE), "absent")
  expect_equal(m1, best_invasion_model(ev, select_on_recal = FALSE))
  ev2 <- suppressWarnings(evaluate_invasion(
    suppressMessages(attach_invasion_recalibration(mk_eval_frame())), n_boot = 0L))
  ev2$log_score_recal[1] <- NA_real_
  expect_warning(best_invasion_model(ev2, select_on_recal = TRUE), "no recalibrated log score")
})

test_that("the pick and the audit leaderboard always use the SAME log-score axis", {
  ev <- suppressWarnings(evaluate_invasion(
    suppressMessages(attach_invasion_recalibration(mk_eval_frame())), n_boot = 0L))
  for (flag in c(FALSE, TRUE)) {
    pick <- suppressMessages(best_invasion_model(ev, select_on_recal = flag))
    tab  <- suppressMessages(invasion_selection_table(ev, select_on_recal = flag))
    expect_equal(pick, tab$method[which.min(tab$composite)])
    expect_equal(unique(tab$log_score_axis),
                 if (flag) "log_score_recal" else "log_score")
  }
})

# ---- the prediction-path delta resolver --------------------------------------
suppressMessages(source(file.path(here::here(), "spatiotemporal", "21_bayesian_renewal.R")))

test_that(".bayes_recal_delta resolves every supported delta form", {
  expect_equal(.bayes_recal_delta(NULL, 1), 1)              # no recalibration
  expect_equal(.bayes_recal_delta(numeric(0), 1), 1)
  expect_equal(.bayes_recal_delta(0.5, 2), 0.5)             # bare scalar -> all horizons
  expect_equal(.bayes_recal_delta(c("1" = 0.5, "2" = 0.4), 2), 0.4)
  expect_equal(.bayes_recal_delta(c("h1" = 0.5, "h2" = 0.4), 1), 0.5)
  expect_equal(.bayes_recal_delta(c("1" = 0.5), 3), 1)      # horizon absent -> untouched
  expect_warning(expect_equal(.bayes_recal_delta(c(0.5, 0.4), 1), 1))   # unnamed vector
  expect_warning(expect_equal(.bayes_recal_delta(c("1" = -1), 1), 1))   # unusable value
  expect_warning(expect_equal(.bayes_recal_delta(c("1" = NA_real_), 1), 1))
})

test_that(".suite_delta_for selects the right rows and declines cleanly", {
  tb <- data.frame(method = c("M1", "M1", "M2"), horizon = c(1L, 2L, 1L),
                   delta = c(0.5, 0.4, 0.7))
  expect_equal(.suite_delta_for(tb, "M1", c(1L, 2L)), c("1" = 0.5, "2" = 0.4))
  expect_null(.suite_delta_for(tb, "absent", c(1L, 2L)))
  expect_null(.suite_delta_for(NULL, "M1", 1L))
  expect_warning(d <- .suite_delta_for(tb, "M2", c(1L, 2L)))   # partial horizon coverage
  expect_equal(d, c("1" = 0.7))
})

message("\n[test_invasion_recalibration] integration assertions executed.")

# 9. RANK-ONLY ROWS TAKE THE IDENTITY ------------------------------------------
test_that("rank-only rows get p_recal = p_invasion and delta = 1, never NA", {
  # Recalibration is strictly monotone, so for a comparator that emits an ordering the
  # recalibrated and raw scales are the same ordering: delta = 1 states exactly that.
  # Leaving NA instead made fs_lfo_col() (forecast_scale.R) fall back to the RAW column for
  # the WHOLE frame -- every panel in both figure suites drawn raw while captioned
  # "recalibrated".
  s <- rbind(mk_slice(horizon = 1L), mk_slice(horizon = 1L, seed = 3L))
  s$method[seq_len(nrow(s) / 2)] <- "calibrated"
  s$method[-seq_len(nrow(s) / 2)] <- "rankonly"
  s$prob_calibrated <- s$method == "calibrated"
  out <- suppressMessages(attach_invasion_recalibration(s))
  ro <- out$method == "rankonly"
  expect_true(all(is.finite(out$p_recal[ro])))
  expect_equal(out$p_recal[ro], out$p_invasion[ro])
  expect_true(all(out$delta_preq[ro] == 1))
  expect_false(any(out$delta_estimable[ro]))
  # and the calibrated arm is untouched by the fill
  expect_true(any(out$delta_estimable[!ro]))
  # No scorable row anywhere is left without a recalibrated probability.
  sc <- is.finite(out$p_invasion) & !(out$was_active_before %in% TRUE)
  expect_equal(sum(sc & !is.finite(out$p_recal)), 0L)
})

test_that("an already-affected rank-only row still gets NA, not the identity", {
  s <- mk_slice(horizon = 1L)
  s$prob_calibrated <- FALSE
  s$was_active_before[1:5] <- TRUE
  out <- suppressMessages(attach_invasion_recalibration(s))
  expect_true(all(is.na(out$p_recal[1:5])))          # not scored => never fabricated
  expect_true(all(is.finite(out$p_recal[-(1:5)])))
})

# 10. EXPANDING WINDOW + THE DELTA-STABILITY DIAGNOSTIC -------------------------
test_that("the prequential delta trains on an EXPANDING window, never a recent slice", {
  # Pinned deliberately. A window would only earn its cost if delta genuinely moved between
  # folds, and invasion_delta_stability() reports it does not (I^2 = 0, tau^2 = 0, Cochran Q
  # non-significant for every method at both horizons). Every fold whose outcome window has
  # closed must therefore be in the training pool.
  s <- mk_slice(n_folds = 9L)
  q <- prequential_invasion_delta(s, horizon = 1L)
  # h=1, 7-day spacing: fold k trains on folds 1..k-2, i.e. n_train_folds = max(0, k-2).
  expect_equal(q$n_train_folds, pmax(0L, seq_len(9L) - 2L))
  expect_true(all(diff(q$n_train_folds) >= 0))          # never shrinks
  expect_equal(max(q$n_train_folds), 7L)                 # the last fold sees ALL prior folds
})

test_that("invasion_delta_stability fits each fold ALONE and reports heterogeneity", {
  if (!exists("invasion_delta_stability", mode = "function"))
    testthat::skip("invasion_delta_stability() unavailable")
  s <- rbind(mk_slice(horizon = 1L), mk_slice(horizon = 2L, seed = 5L))
  st <- suppressWarnings(invasion_delta_stability(s))
  expect_true(all(c("per_fold", "summary") %in% names(st)))
  expect_true(all(c("fold_id","n_atrisk","n_events","delta","se_log","lo","hi")
                  %in% names(st$per_fold)))
  expect_true(all(c("tau2","I2","Q","Q_df","Q_p","trend_slope","trend_p",
                    "extent_slope","extent_p","pooled_delta") %in% names(st$summary)))
  # independence: each fold's delta must depend ONLY on that fold's rows
  one <- st$per_fold[st$per_fold$horizon == 1L & st$per_fold$fold_id == 5L, ]
  d5  <- s[s$horizon == 1L & s$fold_id == 5L, ]
  expect_equal(one$delta[1],
               fit_invasion_delta(d5$p_invasion, d5$is_new_invasion, warn_boundary = FALSE),
               tolerance = 1e-8)
  # the extent covariate is populated (it was silently all-NA when read from
  # was_active_before, which is hard-coded FALSE in the LFO frame)
  expect_true(all(is.finite(st$per_fold$n_atrisk)))
  expect_gt(length(unique(st$per_fold$n_atrisk[st$per_fold$horizon == 1L])), 0L)
  # I^2 and tau^2 are well-formed
  fin <- st$summary[is.finite(st$summary$I2), ]
  expect_true(all(fin$I2 >= 0 & fin$I2 <= 1))
  expect_true(all(fin$tau2 >= 0))
  expect_true(all(fin$Q_p >= 0 & fin$Q_p <= 1, na.rm = TRUE))
  # data simulated with a CONSTANT delta must not look heterogeneous
  expect_lt(median(fin$I2), 0.5)
})
