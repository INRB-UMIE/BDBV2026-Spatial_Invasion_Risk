# =============================================================================
# 16_invasion_eval.R — Invasion-Focused Evaluation (pooled across folds)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# The prediction task is a RARE BINARY event: does a zone with no confirmed cases
# up to the cutoff see its first case within the next h weeks? With only a handful
# of invasion events across folds (order ~8-18, data-dependent), per-fold-then-
# average scoring is degenerate, so discrimination/probability metrics are computed
# on the POOLED at-risk rows. Ranking metrics (Precision@K, rank-of-truth) are
# computed per fold then averaged (ranking is only meaningful within a single forecast).
#
# Metrics reported per (method, horizon):
#   PRIMARY : AUC-PR (+prevalence baseline), log-score (raw), mean rank-of-truth,
#             Precision@K / Recall@K / hit-rate@K, calibration-in-the-large
#   SECONDARY: AUC-ROC, Brier skill vs base rate, ECE
# Cluster-bootstrap 90% CI (resampling by zone) on AUC-PR and log-score.
#
# Count-WIS is NOT computed here and is never a model selector — it rewards
# "predict zero everywhere" over 519 mostly-zero zones; it is kept only as a
# secondary count diagnostic (computed in the reporting/WIS helpers, not here).
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))
suppressPackageStartupMessages({ library(tidyverse) })

# ---------------------------------------------------------------------------
# Core scoring primitives (self-contained)
# ---------------------------------------------------------------------------

.auc_roc <- function(p, y) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.auc_pr <- function(p, y) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  if (sum(y) == 0) return(NA_real_)
  o <- order(p, decreasing = TRUE); ps <- p[o]; y <- y[o]
  tp <- cumsum(y == 1); fp <- cumsum(y == 0); P <- sum(y == 1)
  # Tie-aware Average Precision: collapse runs of EQUAL score to a single
  # operating point (the last index of each run), so tied examples share one
  # precision/recall pair — you cannot separate them at any real threshold.
  # Without this, a positive buried in a large tie group (e.g. the many zones with
  # p_invasion == 0 exactly, no mobility pathway) would earn an arbitrary,
  # row-order-dependent AP that can misrank models. Mirrors sklearn's
  # average_precision_score.
  keep <- c(ps[-length(ps)] != ps[-1], TRUE)
  prec <- (tp / (tp + fp))[keep]; rec <- (tp / P)[keep]
  d_rec <- diff(c(0, rec))
  sum(prec * d_rec, na.rm = TRUE)
}

.log_score <- function(p, y) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- pmin(pmax(p[ok], 1e-6), 1 - 1e-6); y <- y[ok]
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

.brier <- function(p, y) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  mean((p - y)^2)
}

.brier_skill <- function(p, y) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  base <- mean(y); bs_ref <- mean((base - y)^2)
  if (bs_ref <= 0) return(NA_real_)
  1 - mean((p - y)^2) / bs_ref
}

.ece <- function(p, y, n_bins = 10L) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  if (length(p) == 0) return(NA_real_)
  br <- cut(p, breaks = seq(0, 1, length.out = n_bins + 1L), include.lowest = TRUE)
  df <- tibble::tibble(p, y, br)
  s <- df %>% dplyr::group_by(br) %>%
    dplyr::summarise(n = dplyr::n(), pm = mean(p), ym = mean(y), .groups = "drop")
  sum(s$n / sum(s$n) * abs(s$pm - s$ym))
}

# Per-fold ranking metrics, averaged over folds.
.ranking_metrics <- function(d, k_values) {
  folds <- unique(d$fold_id)
  per_fold <- lapply(folds, function(fo) {
    g <- d[d$fold_id == fo, ]
    g <- g[is.finite(g$p_invasion), ]
    if (nrow(g) == 0 || sum(g$is_new_invasion) == 0) return(NULL)
    yv <- as.integer(g$is_new_invasion)
    # Tie-aware rank of the invaded zone(s): ties share the AVERAGE rank (1 =
    # highest risk), mirroring .auc_roc — otherwise a truth tied at p=0 with
    # hundreds of zones would get an arbitrary position-dependent rank. (This is
    # the REPORTING metric "mean rank of invaded zones"; the operational detection
    # metrics — compute_detection_curve / lead_time — deliberately use ties="max"
    # because top-K membership is only credited once ALL tied zones are monitored.)
    rk <- rank(-g$p_invasion, ties.method = "average")
    ranks_of_truth <- rk[yv == 1]
    o <- order(g$p_invasion, decreasing = TRUE); y <- yv[o]
    n_pos <- sum(y)
    res <- list(mean_rank = mean(ranks_of_truth), n_atrisk = length(y), n_pos = n_pos)
    for (K in k_values) {
      topK <- head(y, K)
      res[[paste0("prec@", K)]] <- sum(topK) / length(topK)   # realised count (fewer than K at-risk zones)
      res[[paste0("recall@", K)]] <- sum(topK) / n_pos
      res[[paste0("hit@", K)]]   <- as.integer(sum(topK) > 0)
    }
    res
  })
  per_fold <- Filter(Negate(is.null), per_fold)
  if (length(per_fold) == 0) return(NULL)
  keys <- names(per_fold[[1]])
  out <- sapply(keys, function(k) mean(vapply(per_fold, function(x) x[[k]], numeric(1)), na.rm = TRUE))
  as.list(out)
}

# ---------------------------------------------------------------------------
# Master invasion evaluation
# ---------------------------------------------------------------------------

#' Evaluate invasion forecasts on the at-risk set, pooled across folds.
#'
#' @param lfo_results tibble with method, horizon, fold_id, health_zone,
#'   p_invasion, is_new_invasion, was_active_before.
#' @param k_values top-K thresholds for the ranking metrics.
#' @param n_boot cluster-bootstrap replicates (resample zones) for CIs.
#' @return tibble: one row per method × horizon with the metric set + CIs.
evaluate_invasion <- function(lfo_results, k_values = c(5L, 10L, 15L),
                              n_boot = 400L, seed = RANDOM_SEED) {
  stopifnot(all(c("method", "horizon", "fold_id", "health_zone",
                  "p_invasion", "is_new_invasion") %in% names(lfo_results)))
  d0 <- lfo_results
  if ("was_active_before" %in% names(d0)) {
    d0 <- d0[!as.logical(d0$was_active_before) %in% TRUE, ]
  }
  d0 <- d0[is.finite(d0$p_invasion), ]
  d0$.cell <- paste(d0$fold_id, d0$health_zone, sep = "\r")   # (fold x zone) cell id

  # --- Common-support alignment so AUC-PR skill DENOMINATORS are comparable across models.
  # A model that covers fewer LFO folds (e.g. a Bayesian fit that could not be estimated on
  # an event-sparse fold) would otherwise be scored on a different, often base-rate-enriched,
  # (fold x zone) subset and mis-ranked purely through its denominator. For each horizon we
  # take the maximal per-model fold coverage; models at that coverage are "eligible" and are
  # all scored on the INTERSECTION of their cells with a single shared base rate; models below
  # it keep their native rows but are flagged partial_cv and sorted out of the comparable pack.
  .cover <- function(dh, m) length(unique(dh$fold_id[dh$method == m]))
  hsupport <- lapply(sort(unique(d0$horizon)), function(h) {
    dh <- d0[d0$horizon == h, ]; ms <- unique(dh$method)
    cov <- vapply(ms, function(m) .cover(dh, m), integer(1)); maxf <- max(cov)
    elig <- ms[cov == maxf]
    list(common = Reduce(intersect, lapply(elig, function(m) unique(dh$.cell[dh$method == m]))),
         maxf = maxf)
  })
  names(hsupport) <- as.character(sort(unique(d0$horizon)))

  set.seed(seed)
  keys <- unique(d0[, c("method", "horizon")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    meth <- keys$method[i]; h <- keys$horizon[i]
    hs <- hsupport[[as.character(h)]]
    d_all <- d0[d0$method == meth & d0$horizon == h, ]
    n_folds <- length(unique(d_all$fold_id)); partial_cv <- n_folds < hs$maxf
    # Full-coverage models are scored on the shared common support; partial models keep their
    # native rows so their numbers are transparent but they do not contaminate the leaderboard.
    d <- if (!partial_cv && length(hs$common)) d_all[d_all$.cell %in% hs$common, ] else d_all
    y <- as.integer(d$is_new_invasion); p <- d$p_invasion
    n_pos <- sum(y); n <- length(y)
    if (n == 0) return(NULL)

    base_rate <- mean(y)
    aucpr <- .auc_pr(p, y); aucroc <- .auc_roc(p, y)
    lsc <- .log_score(p, y); bss <- .brier_skill(p, y); ece <- .ece(p, y)
    cal_in_large <- mean(p, na.rm = TRUE) / max(base_rate, 1e-9)  # 1 = perfect

    rk <- .ranking_metrics(d, k_values)

    # cluster bootstrap by zone for AUC-PR + log-score
    zones <- unique(d$health_zone)
    boot_aucpr <- numeric(0); boot_ls <- numeric(0)
    if (n_pos >= 2 && length(zones) > 5) {
      for (b in seq_len(n_boot)) {
        zs <- sample(zones, length(zones), replace = TRUE)
        idx <- unlist(lapply(zs, function(z) which(d$health_zone == z)))
        yb <- y[idx]; pb <- p[idx]
        if (sum(yb) >= 1) {
          boot_aucpr <- c(boot_aucpr, .auc_pr(pb, yb))
          boot_ls    <- c(boot_ls, .log_score(pb, yb))
        }
      }
    }
    # 90% bootstrap interval (q5/q95) — the pipeline's uncertainty-band convention
    # for forecast/metric intervals (never 95%); the covariate-screen regression HRs
    # separately use standard 95% Wald CIs (see METHODS §5.5/§9.5).
    ci <- function(v) if (length(v) > 10) quantile(v, c(0.05, 0.95), na.rm = TRUE) else c(NA, NA)
    aucpr_ci <- ci(boot_aucpr); ls_ci <- ci(boot_ls)

    tibble::tibble(
      method = meth, horizon = h, n_atrisk = n, n_invasions = n_pos,
      base_rate = base_rate, n_folds = n_folds,
      coverage = n_folds / hs$maxf, partial_cv = partial_cv,
      auc_pr = aucpr, auc_pr_lo = aucpr_ci[1], auc_pr_hi = aucpr_ci[2],
      auc_pr_skill = if (is.finite(aucpr)) aucpr / max(base_rate, 1e-9) else NA_real_,
      auc_roc = aucroc,
      log_score = lsc, log_score_lo = ls_ci[1], log_score_hi = ls_ci[2],
      brier_skill = bss, ece = ece, calibration_in_large = cal_in_large,
      mean_rank_of_truth = rk$mean_rank %||% NA_real_,
      prec_at_5 = rk[["prec@5"]] %||% NA_real_,
      recall_at_5 = rk[["recall@5"]] %||% NA_real_,
      hit_at_5 = rk[["hit@5"]] %||% NA_real_,
      prec_at_10 = rk[["prec@10"]] %||% NA_real_,
      recall_at_10 = rk[["recall@10"]] %||% NA_real_,
      hit_at_10 = rk[["hit@10"]] %||% NA_real_,
      prec_at_15 = rk[["prec@15"]] %||% NA_real_,
      recall_at_15 = rk[["recall@15"]] %||% NA_real_,
      hit_at_15 = rk[["hit@15"]] %||% NA_real_
    )
  })
  out <- dplyr::bind_rows(Filter(Negate(is.null), rows))
  # Primary ranking: full-CV (comparable common-support) models first, then AUC-PR skill
  # (discrimination over the shared base rate), tie-broken by lower mean rank-of-truth.
  # partial_cv models sort last so they cannot outrank the comparable pack on a thin subset.
  out %>% dplyr::arrange(horizon, partial_cv, dplyr::desc(auc_pr_skill), mean_rank_of_truth)
}

#' Pick the best model POOLED across forecast horizons: discrimination-led, spike-guarded.
#'
#' The featured model drives the operational watch-list (rank zones, monitor the top K),
#' so the primary criterion is DISCRIMINATION — total AUC-PR skill pooled across horizons
#' (higher better). Ranking on discrimination ALONE, however, would crown a "spiky" variant
#' that nails a handful of top predictions (high AUC-PR) while scattering the rest, so a
#' quality GATE first drops any model whose worst-horizon mean rank-of-truth exceeds
#' `MRT_GATE_MULT` x the field median (i.e. ranks the truly-invaded zones much worse than a
#' typical model). Among the surviving (well-targeted) models the winner is the highest total
#' AUC-PR skill; ties are broken by better (lower) total rank-of-truth, then lower log-score
#' (a proper score penalising over-confidence). Only methods present at the maximal number of
#' pooled horizons are eligible, so a model cannot win by sitting out the harder horizon; the
#' gate falls back to the ungated set if it (or missing rank-of-truth) would leave nothing.
#'
#' @param horizons optional integer vector; if given, only these horizons are
#'   pooled (default: every horizon present in eval_tbl).
#' @param restrict optional regex; if given, only methods matching it are ranked
#'   (e.g. "^Renewal" to pick the best mobility-renewal model).
MRT_GATE_MULT <- 1.5   # spiky-model gate: drop methods with worst-horizon rank-of-truth > this x field median
best_invasion_model <- function(eval_tbl, horizons = NULL, restrict = NULL) {
  e <- eval_tbl
  if (!is.null(restrict)) {
    # STRICT: if a family filter is given and nothing matches, return NA rather
    # than silently falling back to the unrestricted set (which would return an
    # out-of-family model, e.g. a frequentist model for a Bayesian request).
    # Callers supply their own defaults.
    e <- e %>% dplyr::filter(grepl(restrict, method))
  }
  if (!is.null(horizons)) e <- e %>% dplyr::filter(horizon %in% !!horizons)
  if (nrow(e) == 0) return(NA_character_)
  # Aggregate per method across horizons; keep only methods covering the most horizons.
  agg <- e %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(n_h     = dplyr::n_distinct(horizon),
                     aps     = sum(dplyr::coalesce(auc_pr_skill, 0)),      # total AUC-PR skill (lead)
                     mrt_max = max(dplyr::coalesce(mean_rank_of_truth, Inf)),  # worst-horizon targeting
                     mrt_sum = sum(dplyr::coalesce(mean_rank_of_truth, Inf)),
                     ls_sum  = sum(dplyr::coalesce(log_score, Inf)),
                     .groups = "drop") %>%
    dplyr::filter(n_h == max(n_h))
  # spiky-model gate (relative to the field), with graceful fallback to the ungated set.
  finite_mrt <- agg$mrt_max[is.finite(agg$mrt_max)]
  gate <- if (length(finite_mrt)) MRT_GATE_MULT * stats::median(finite_mrt) else Inf
  keep <- agg %>% dplyr::filter(mrt_max <= gate)
  if (nrow(keep) == 0) keep <- agg
  keep %>% dplyr::arrange(dplyr::desc(aps), mrt_sum, ls_sum) %>% dplyr::slice(1) %>% dplyr::pull(method)
}

#' Held-out (last-k-origins) confirmatory evaluation + optimism gap (review §2.2).
#'
#' The all-fold LFO both SELECTS the featured model and reports its skill, so that
#' headline is selection-optimistic (winner's curse over many candidates on few
#' events). This reserves the last `n_holdout` weekly origins as an OUTER test set
#' touched only once: the model is selected on the earlier (inner) folds by the same
#' composite as best_invasion_model(), then scored on the untouched outer folds. It
#' returns the inner-CV skill (for comparability with the all-fold number), the outer
#' held-out skill (the honest, most-recent out-of-sample estimate), and the optimism
#' gap between them. The outer origins sit in the recent window, so this doubles as the
#' formal front end of the prospective evaluation (§3.8).
#'
#' @param lfo_results pooled LFO tibble (needs method, horizon, cutoff, plus the
#'   columns evaluate_invasion() requires).
#' @param n_holdout number of most-recent origins held out (default 2).
#' @param restrict optional method regex to select WITHIN (e.g. "^Bayes").
#' @param horizon horizon at which the skill/gap is reported (default 1).
#' @return list(selected, inner_skill, outer_skill, optimism_gap, inner_eval,
#'   outer_eval, inner_cutoffs, outer_cutoffs), or NULL if too few origins.
evaluate_invasion_heldout <- function(lfo_results, n_holdout = 2L,
                                      restrict = "^Bayes", horizon = 1L) {
  if (is.null(lfo_results) || !nrow(lfo_results) || !"cutoff" %in% names(lfo_results)) return(NULL)
  cutoffs <- sort(unique(lfo_results$cutoff))
  if (length(cutoffs) < n_holdout + 2L) {
    warning("[heldout] too few origins (", length(cutoffs), ") for a ", n_holdout,
            "-origin holdout; skipping"); return(NULL)
  }
  outer_cut <- utils::tail(cutoffs, n_holdout)
  inner_cut <- setdiff(cutoffs, outer_cut)
  inner_res <- lfo_results[lfo_results$cutoff %in% inner_cut, , drop = FALSE]
  outer_res <- lfo_results[lfo_results$cutoff %in% outer_cut, , drop = FALSE]
  inner_eval <- tryCatch(evaluate_invasion(inner_res), error = function(e) NULL)
  outer_eval <- tryCatch(evaluate_invasion(outer_res), error = function(e) NULL)
  if (is.null(inner_eval) || is.null(outer_eval)) return(NULL)
  # Select on the INNER folds only (never on the held-out outer folds).
  sel_tbl <- if (!is.null(restrict)) dplyr::filter(inner_eval, grepl(restrict, method)) else inner_eval
  sel <- best_invasion_model(sel_tbl)
  if (length(sel) != 1L || is.na(sel)) return(NULL)
  gk <- function(ev, m) { r <- ev$auc_pr_skill[ev$method == m & ev$horizon == horizon]
                          if (length(r)) r[1] else NA_real_ }
  inner_skill <- gk(inner_eval, sel); outer_skill <- gk(outer_eval, sel)
  list(selected = sel, horizon = horizon, n_holdout = n_holdout,
       inner_skill = inner_skill, outer_skill = outer_skill,
       optimism_gap = inner_skill - outer_skill,
       inner_eval = inner_eval, outer_eval = outer_eval,
       inner_cutoffs = inner_cut, outer_cutoffs = outer_cut)
}

# ---------------------------------------------------------------------------
# Clean invasion LFO-CV (same folds for all models; correct at-risk outcome)
# ---------------------------------------------------------------------------

#' Leave-future-out cross-validation focused on the invasion task.
#'
#' Every model is run on the SAME principled folds. For fold cutoff C:
#'   - training data = counts as KNOWN at the real-time forecast moment C+7,
#'     nowcast-corrected as of C+7. When `linelist` is supplied these are
#'     reconstructed leakage-free by re-aggregating the line list censored to the
#'     sample-observation date <= C+7 (so late-reported cases are excluded);
#'     otherwise the FINAL onset-bucketed counts sliced to week <= C are used,
#'     which carry a mild training-side revision leak (recent weeks include cases
#'     only reported after C+7, then further inflated by the nowcast);
#'   - at-risk set  = zones with zero confirmed cases in weeks <= C;
#'   - outcome(h)   = zone is first confirmed in weeks C+1..C+h (cumulative
#'     window; matches the models' cumulative invasion probability, so a first
#'     invasion is counted once per horizon and, being affected thereafter, is
#'     excluded from later folds — no double counting).
#' Folds whose farthest-horizon eval week is not yet ~complete (delay-truncated)
#' are dropped so the ground truth is reliable.
#'
#' @param zone_week_raw raw zone-week counts (health_zone, week_start, confirmed,
#'   suspected), outbreak period.
#' @param models named list of forecasters, each `function(zw_nc, t_idx,
#'   horizons, cutoff)` returning health_zone, horizon, mu_forecast, p_invasion,
#'   method (+ optional p_infection_invasion).
#' @param nowcast_fn function(zone_week, analysis_date) -> nowcast-corrected;
#'   default apply_nowcast_correction (must be sourced).
#' @param linelist optional line list (dat$ll). When supplied, per-fold training
#'   counts are reconstructed leakage-free via reaggregate_asof() (22_daily_reissue.R)
#'   instead of slicing the final weekly counts. Recommended for a truly real-time
#'   evaluation; NULL reproduces the previous (mildly leaky) behaviour.
#' @return pooled tibble across folds×models with p_invasion, is_new_invasion,
#'   was_active_before, fold_id, cutoff, method, horizon.
run_invasion_lfo <- function(zone_week_raw, models, horizons = LFO_HORIZONS,
                             min_train_weeks = 3L, min_atrisk_events = 1L,
                             min_eval_age_days = 10L,
                             analysis_date = ANALYSIS_DATE,
                             nowcast_fn = NULL, linelist = NULL) {
  if (is.null(nowcast_fn)) nowcast_fn <- get0("apply_nowcast_correction")
  stopifnot(!is.null(nowcast_fn))
  weeks   <- sort(unique(zone_week_raw$week_start))
  zones   <- sort(unique(zone_week_raw$health_zone))

  # Cumulative-confirmed-by-week helper
  conf_by_week <- zone_week_raw %>%
    dplyr::group_by(health_zone, week_start) %>%
    dplyr::summarise(confirmed = sum(confirmed, na.rm = TRUE), .groups = "drop")

  affected_by <- function(cut) {
    conf_by_week %>% dplyr::filter(week_start <= cut, confirmed > 0) %>%
      dplyr::pull(health_zone) %>% unique()
  }
  first_case_week <- conf_by_week %>% dplyr::filter(confirmed > 0) %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(first_wk = min(week_start), .groups = "drop")

  # Candidate cutoffs (Q1 fold augmentation): a cutoff qualifies once its
  # SHORTEST-horizon eval week is old enough to be reliable ground truth (rather
  # than requiring the FARTHEST horizon complete). Each fold then scores only the
  # horizons whose eval window is itself complete — a per-horizon completeness
  # guard. This adds the recent cutoffs where h=1 is evaluable but longer
  # horizons are not, giving materially more folds for the near-term forecast
  # that matters most, with no leakage (an incomplete horizon is simply skipped).
  hmin <- min(horizons)
  cand <- weeks[weeks >= weeks[min_train_weeks] &
                  (as.numeric(analysis_date - (weeks + hmin * 7 + 7)) >= min_eval_age_days)]
  folds <- list(); fid <- 0L
  for (ci in seq_along(cand)) {
    cut <- cand[ci]                       # index (not `for x in Date` which unclasses)
    # horizons whose eval window is reliably complete at this cutoff
    complete_h <- horizons[
      as.numeric(analysis_date - (cut + horizons * 7 + 7)) >= min_eval_age_days]
    if (length(complete_h) == 0L) next
    hmax_c <- max(complete_h)
    aff <- affected_by(cut)
    atrisk <- setdiff(zones, aff)
    # invasion events within the largest COMPLETE window among at-risk zones.
    # NOTE (deliberate): folds with < min_atrisk_events invasions are dropped so
    # per-fold ranking/discrimination is defined. Because scoring is POOLED, this
    # mildly enriches the pooled base rate (event-free at-risk weeks are removed),
    # so `calibration_in_large` is a slightly optimistic read of absolute
    # over-prediction — model *rankings* (all scored on the identical pool) are
    # unaffected. Set min_atrisk_events = 0 to score every complete fold instead.
    ev <- first_case_week %>%
      dplyr::filter(health_zone %in% atrisk,
                    first_wk > cut, first_wk <= cut + hmax_c * 7)
    if (nrow(ev) < min_atrisk_events) next
    fid <- fid + 1L
    folds[[fid]] <- list(cutoff = cut, atrisk = atrisk,
                         horizons = complete_h, n_events = nrow(ev))
  }
  if (length(folds) == 0) {
    warning("[invasion_lfo] No evaluable folds."); return(NULL)
  }
  message(sprintf("[invasion_lfo] %d folds (cutoffs %s); horizons/fold: %s; events/fold: %s",
                  length(folds),
                  paste(vapply(folds, function(f) format(f$cutoff, "%m-%d"), ""), collapse = ","),
                  paste(vapply(folds, function(f) paste(f$horizons, collapse = "+"), ""), collapse = ","),
                  paste(vapply(folds, function(f) f$n_events, 0L), collapse = ",")))

  out <- list()
  for (fi in seq_along(folds)) {
    f <- folds[[fi]]; cut <- f$cutoff
    # Leakage-free training counts: reconstruct what surveillance actually KNEW at
    # the real-time forecast moment (cut+7) by re-aggregating the line list censored
    # to that observation date (reaggregate_asof, 22_daily_reissue.R), rather than
    # slicing the FINAL onset-bucketed counts — which fold in cases only reported
    # after cut+7 and would then be further inflated by the nowcast (the training-side
    # revision leak). Falls back to the final-count slice when no line list is given
    # or the re-aggregation fails, so behaviour degrades gracefully rather than erroring.
    train_raw <- NULL
    if (!is.null(linelist) && exists("reaggregate_asof")) {
      train_raw <- tryCatch(
        reaggregate_asof(linelist, zones, cut + 7, week_spine = weeks[weeks <= cut]),
        error = function(e) {
          # Warn LOUDLY: a silent NULL here drops this fold onto the final-count (revision-leaky)
          # path even under LEAKAGE_FREE_LFO, invisibly weakening the leakage guard for that fold.
          warning(sprintf("[lfo] fold %s: reaggregate_asof failed (%s); this fold falls back to the final-count (revision-leaky) training slice.",
                          format(cut), conditionMessage(e)), call. = FALSE); NULL })
      if (!is.null(train_raw))
        train_raw <- train_raw %>% dplyr::filter(week_start <= cut)
    }
    if (is.null(train_raw))
      train_raw <- zone_week_raw %>% dplyr::filter(week_start <= cut)
    train_nc  <- tryCatch(nowcast_fn(train_raw, analysis_date = cut + 7),
                          error = function(e) {
                            warning(sprintf("[lfo] fold %s: nowcast failed (%s); using UNCORRECTED counts for this fold.",
                                            format(cut), conditionMessage(e)), call. = FALSE)
                            # Ensure the count column the models read (confirmed_nc) EXISTS: without it
                            # .count_wide() returns an all-zero matrix and the whole fold silently
                            # collapses to zero-hazard. Fall back to the raw (uncorrected) confirmed count.
                            if (!"confirmed_nc" %in% names(train_raw) && "confirmed" %in% names(train_raw))
                              train_raw$confirmed_nc <- train_raw$confirmed
                            train_raw
                          })
    t_idx <- length(unique(train_nc$week_start))
    aff_set <- f$atrisk  # NB: atrisk vector; affected = complement

    # ground-truth outcomes per COMPLETE horizon (cumulative window), at-risk only
    truth_h <- lapply(f$horizons, function(h) {
      inv_zones <- first_case_week %>%
        dplyr::filter(health_zone %in% f$atrisk,
                      first_wk > cut, first_wk <= cut + h * 7) %>%
        dplyr::pull(health_zone)
      tibble::tibble(health_zone = f$atrisk,
                     horizon = h,
                     is_new_invasion = as.integer(f$atrisk %in% inv_zones))
    }) %>% dplyr::bind_rows()

    # Fit every model on this fold. Each (fold, model) fit is independent and each
    # brm() is independently seeded, so the fits can be fanned out across a worker
    # pool with results identical to sequential. Fold 1 stays SEQUENTIAL to warm the
    # cmdstanr compile cache (all distinct Stan programs) before the parallel folds
    # reuse it; PARALLEL_JOBS = 1 keeps the whole loop sequential (historical path).
    fit_one <- function(mname) {
      fc <- tryCatch(models[[mname]](train_nc, t_idx, horizons, cut),
                     error = function(e) {
                       warning(sprintf("[invasion_lfo] %s failed on fold %d (%s): %s",
                                       mname, fi, cut, conditionMessage(e))); NULL })
      if (is.null(fc) || nrow(fc) == 0) return(NULL)
      fc$method <- mname
      # keep only at-risk zones (drop affected), attach outcome
      fc %>%
        dplyr::filter(health_zone %in% f$atrisk) %>%
        dplyr::select(dplyr::any_of(c("health_zone", "horizon", "mu_forecast",
                                      "p_invasion", "p_infection_invasion", "method"))) %>%
        dplyr::inner_join(truth_h, by = c("health_zone", "horizon")) %>%
        dplyr::mutate(fold_id = fi, cutoff = cut, was_active_before = FALSE)
    }
    .par <- (get0("PARALLEL_JOBS", ifnotfound = 1L) > 1L) && fi > 1L &&
            requireNamespace("furrr", quietly = TRUE)
    .tf <- Sys.time()
    fold_res <- if (.par)
      furrr::future_map(names(models), fit_one,
                        .options = furrr::furrr_options(seed = TRUE))
    else lapply(names(models), fit_one)
    message(sprintf("[lfo-timing] fold %d/%d: %d models %s in %.1fs",
                    fi, length(folds), length(models),
                    if (.par) "PARALLEL" else "seq",
                    as.numeric(difftime(Sys.time(), .tf, units = "secs"))))
    out <- c(out, fold_res)   # future_map preserves input order; NULLs dropped by bind_rows
  }
  dplyr::bind_rows(out)
}


message("[invasion_eval] 16_invasion_eval.R loaded — pooled invasion-focused scoring.")
