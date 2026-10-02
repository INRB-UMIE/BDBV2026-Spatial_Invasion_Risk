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
  if (!length(y)) return(NA_real_)     # mean(numeric(0)) is NaN, which reads as a real score
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

.brier <- function(p, y) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  if (!length(y)) return(NA_real_)
  mean((p - y)^2)
}

.brier_skill <- function(p, y) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  # Empty after filtering => mean() is NaN => `if (bs_ref <= 0)` was `if (NA)`, a hard error
  # inside evaluate_invasion()'s per-key lapply that would abort the ENTIRE evaluation.
  if (!length(y)) return(NA_real_)
  base <- mean(y); bs_ref <- mean((base - y)^2)
  if (!is.finite(bs_ref) || bs_ref <= 0) return(NA_real_)
  1 - mean((p - y)^2) / bs_ref
}

#' Deterministic PER-HORIZON bootstrap seed, shared by every method.
#'
#' It hashed the METHOD NAME too, which gave each method a different resample stream -- the
#' exact opposite of the pairing the call site claimed to deliver. Since every method is now
#' scored on the same (fold x zone) support its zone list is identical, so seeding on the
#' horizon alone makes the resamples IDENTICAL across methods: the intervals become genuinely
#' paired and can be differenced, while staying invariant to the composition and order of the
#' grid (which was the other, already-achieved goal). Point estimates are unaffected either
#' way; only auc_pr_lo/hi, log_score_lo/hi and n_boot_aucpr move.
.key_seed <- function(seed, horizon) {
  k <- as.character(horizon)
  # Polynomial rolling hash in DOUBLE arithmetic, reduced every step so nothing can leave the
  # signed-32-bit range set.seed() accepts. (bitwXor() returns NA above .Machine$integer.max,
  # which is why an FNV-style 32-bit unsigned hash cannot be used here.)
  h <- 0
  for (cc in utils::head(as.integer(charToRaw(k)), 256L))
    h <- (h * 31 + cc) %% 2147483647
  as.integer((as.numeric(seed) + h) %% 2147483647)
}

#' Expected calibration error over EQUAL-WIDTH probability bins.
#'
#' READ THIS BEFORE QUOTING IT AS AN INDEPENDENT CALIBRATION METRIC. At this task's base rate
#' (~0.6% at h=1) almost every forecast probability is small, so equal-width decile bins put
#' nearly all the mass in the first bin and the statistic collapses onto the in-the-large
#' miscalibration that `calibration_in_large` already reports. Measured on the featured model
#' at h=1 (2026-09-19): 5,563 of 5,757 rows (96.6%) fall in [0, 0.1], ECE = 0.0063291 against
#' |mean(p) - mean(y)| = 0.0062877 — i.e. 99.3% of the "ECE" IS calibration-in-the-large.
#' Equal-frequency deciles collapse the same way, for the same reason.
#'
#' It is kept because it is a conventional summary and is comparable ACROSS models on these
#' folds, but it carries essentially no information about refinement beyond the intercept, and
#' the reliability diagrams (FigR1_reliability_h1/_h2) are what actually show the shape.
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
    n_pos <- sum(yv)
    res <- list(mean_rank = mean(ranks_of_truth), n_atrisk = length(yv), n_pos = n_pos)
    # TIE-AWARE top-K (was head(order(p), K), which breaks ties by ROW ORDER — i.e. by
    # zone name after the joins — making prec@K/recall@K/hit@K arbitrary whenever scores
    # tie, and inconsistent with mean_rank just above). A zone counts as "in the top K"
    # only if monitoring at most K zones necessarily includes it: rank(-p, "max") <= K.
    # This is the SAME convention compute_detection_curve()/lead_time already use, and it
    # is the conservative one — a tie group straddling the K boundary is credited to none
    # of its members, because you cannot pick a subset of it without monitoring all of it.
    # Consequence (intended): if every at-risk zone ties — e.g. all p == 0 because no
    # mobility pathway exists — no zone is in the top K and all three metrics are 0, which
    # correctly scores a model that cannot discriminate at all.
    rk_max <- rank(-g$p_invasion, ties.method = "max")
    for (K in k_values) {
      inK <- rk_max <= K
      nk  <- sum(inK)
      hits <- sum(yv[inK])
      res[[paste0("prec@", K)]]   <- if (nk > 0) hits / nk else 0   # realised size (may be < K)
      res[[paste0("recall@", K)]] <- if (n_pos > 0) hits / n_pos else NA_real_
      res[[paste0("hit@", K)]]    <- as.integer(hits > 0)
    }
    res
  })
  per_fold <- Filter(Negate(is.null), per_fold)
  if (length(per_fold) == 0) return(NULL)
  keys <- names(per_fold[[1]])
  out <- sapply(keys, function(k) mean(vapply(per_fold, function(x) x[[k]], numeric(1)), na.rm = TRUE))
  as.list(out)
}

#' Shared (fold x zone) support for ONE horizon, plus that horizon's widest fold coverage.
#'
#' `dh` is the already-filtered scored rows for one horizon and must carry a `.cell` column.
#' Returns list(common = <cell ids>, maxf = <widest per-method fold count>).
#'
#' AUC-PR skill, top-K recall and the random-targeting reference are all measured against the
#' base rate of the row set they are computed on, so two methods scored on different row sets
#' cannot be compared -- the difference in their denominators is read as a difference in
#' skill. Every method is therefore scored on the intersection of the cells all of them cover.
#'
#' THE ANCHOR IS EVERY METHOD, NOT THE BEST-COVERED ONES (fixed 2026-09-22). This previously
#' took the maximal per-model fold coverage, treated only the models AT that coverage as
#' eligible, intersected over those, and left every other model on its own native rows. That
#' works when one or two fits fail a fold. It inverts when a whole FAMILY covers fewer folds:
#' the rolling-predictor floor (20_forecast_detail.R) costs the Bayesian grid its earliest
#' fold, so the three structural baselines became the maximal-coverage pack and all 20
#' Bayesian models were scored on their own, smaller, differently-based support -- 12 folds
#' against 11 at h=1, 11 against 10 at h=2. That extra origin was event-poor, so it inflated
#' the baselines' AUC-PR skill by 11.3 points at h=1 and 9.1 at h=2: Gravity-B4 ranked FIRST
#' of 23 at h=1, ahead of every Bayesian model, and thirteenth once aligned. No Bayesian
#' score moves, because the shared support is the one they were already scored on.
#'
#' Nothing is discarded: every fold stays in lfo_results.rds, and the `coverage`, `n_folds`
#' and `n_folds_scored` columns report each method's native reach. Only the comparison is put
#' on one footing.
.invasion_support <- function(dh, label = NA) {
  ms <- unique(dh$method)
  if (!nrow(dh) || !length(ms)) return(list(common = character(0), maxf = 0L))
  cells_of <- function(m) unique(dh$.cell[dh$method == m])
  n_fold <- vapply(ms, function(m) length(unique(dh$fold_id[dh$method == m])), integer(1))
  maxf   <- max(n_fold)
  common <- Reduce(intersect, lapply(ms, cells_of))
  # A pathological method (one that covers almost nothing) would otherwise silently collapse
  # the support for everyone. Refuse that rather than publish a leaderboard scored on a
  # handful of cells: fall back to the best-covered pack and say so.
  widest <- max(vapply(ms, function(m) length(cells_of(m)), integer(1)))
  if (length(common) < 0.5 * widest) {
    warning(sprintf(paste0("[eval] h=%s: the all-method common support holds %d cells against ",
                           "a widest single-method coverage of %d. One method is dragging it ",
                           "down, so the comparison falls back to the best-covered pack. ",
                           "Check per-method coverage before reading the table."),
                    label, length(common), widest), call. = FALSE, immediate. = TRUE)
    common <- Reduce(intersect, lapply(ms[n_fold == maxf], cells_of))
  }
  list(common = common, maxf = maxf)
}

#' The (fold x zone) cells every scored method covers, at one horizon.
#'
#' Wraps `.invasion_support()` so evaluate_invasion() and compute_detection_curve()
#' (20_forecast_detail.R) resolve the support through THE SAME code, fallback included: the
#' manuscript panels annotate a recall from the curve beside a recall_at_K from the evaluation
#' table. When only the first applied a restriction, the curve's random reference
#' `k / n_atrisk` differed between method families by about 0.4%, which is small but enough to
#' make Figure 2's one-row-per-k assumption false and abort it.
#'
#' @return character vector of cell ids ("<fold_id>\r<health_zone>"), or NULL.
invasion_common_cells <- function(lfo_results, horizon) {
  d <- lfo_results[lfo_results$horizon == horizon, , drop = FALSE]
  if ("was_active_before" %in% names(d))
    d <- d[!(as.logical(d$was_active_before) %in% TRUE), , drop = FALSE]
  d <- d[is.finite(d$p_invasion), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$.cell <- paste(d$fold_id, d$health_zone, sep = "\r")
  common <- .invasion_support(d, label = horizon)$common
  if (!length(common)) NULL else common
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

  # Every method is scored on ONE shared (fold x zone) support so the AUC-PR skill, top-K and
  # random-targeting denominators are the same for all of them. See `.invasion_support()`.
  .horizons <- sort(unique(d0$horizon))
  hsupport  <- lapply(.horizons, function(h) .invasion_support(d0[d0$horizon == h, ], label = h))
  names(hsupport) <- as.character(.horizons)

  keys <- unique(d0[, c("method", "horizon")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    meth <- keys$method[i]; h <- keys$horizon[i]
    # PER-HORIZON SEED — NOT one set.seed(seed) before the loop, and NOT one per method.
    # A single seed made each (method, horizon) consume a different stretch of one stream, with
    # two consequences for PUBLISHED intervals: (a) two methods' CIs came from DIFFERENT zone
    # resamples, so they could not be differenced or read as a paired comparison, and (b)
    # adding, removing or reordering ONE method shifted the CIs of every method after it in
    # `keys`. Hashing the horizon alone fixes both: every method at a horizon draws the SAME
    # zone resamples (paired), and the draw does not depend on the grid's composition or order.
    set.seed(.key_seed(seed, h))
    hs <- hsupport[[as.character(h)]]
    d_all <- d0[d0$method == meth & d0$horizon == h, ]
    n_folds <- length(unique(d_all$fold_id))
    d <- if (length(hs$common)) d_all[d_all$.cell %in% hs$common, ] else d_all
    # `partial_cv` means "not comparable with the rest of the table", which since the
    # all-method alignment is SCORED ON FEWER CELLS THAN THE SHARED SUPPORT -- no longer
    # "fitted on fewer folds". The two coincided only while the anchor was the best-covered
    # pack, because a short method then kept its own rows.
    #
    # The distinction is load-bearing, because every consumer DROPS flagged rows:
    # 17_invasion_viz.R from the figures, update_bayesian_report.R from the report, and
    # best_invasion_model() / invasion_selection_table() from the featured pick -- which
    # 30_projection_config.R also reads as the cascade kernel. The rolling-predictor floor
    # costs the whole Bayesian grid its earliest fold BY DESIGN, so the old definition flagged
    # all 20 Bayesian models and left the three baselines clean: leaderboard, figures and
    # cascade would all have fallen back to a baseline over rows identical to the ones the
    # Bayesian models were scored on. Native reach stays published as n_folds/coverage.
    #
    # So this is FALSE for every method on the all-method support (the intersection is a
    # subset of each method's cells), and TRUE only under the fallback in .invasion_support(),
    # where a short method can genuinely miss cells -- the case the flag was written for.
    partial_cv <- length(hs$common) > 0L &&
      length(unique(d$.cell)) < length(hs$common)
    y <- as.integer(d$is_new_invasion); p <- d$p_invasion
    n_pos <- sum(y); n <- length(y)
    if (n == 0) return(NULL)

    base_rate <- mean(y)
    # RANK-ONLY MODELS. A comparator may emit an ordering rather than a probability
    # (prob_calibrated = FALSE; see 05_baseline_models.R — Distance-B1 has no fitted scale
    # at all, Adjacency-B7 is 1/(1 + travel-time minutes)). Every metric that reads the
    # LEVEL of p is then meaningless for it: replacing the score by any other monotone
    # transform leaves the ranking identical and moves the log score, Brier skill, ECE and
    # calibration-in-large arbitrarily. Those columns were nonetheless computed and
    # published, and Adjacency-B7 out-scored most Bayesian models on Brier skill purely
    # because 1/(1 + minutes) happens to average near the 0.6% base rate — a coincidence of
    # units, published as evidence. They are reported NA here. The RANK metrics (AUC-PR,
    # AUC-ROC, rank-of-truth, top-K) are invariant to monotone transforms and stay.
    .calib <- if ("prob_calibrated" %in% names(d)) !any(d$prob_calibrated %in% FALSE) else TRUE
    aucpr <- .auc_pr(p, y); aucroc <- .auc_roc(p, y)
    lsc <- if (.calib) .log_score(p, y) else NA_real_
    bss <- if (.calib) .brier_skill(p, y) else NA_real_
    ece <- if (.calib) .ece(p, y) else NA_real_
    cal_in_large <- if (.calib) mean(p, na.rm = TRUE) / max(base_rate, 1e-9) else NA_real_  # 1 = perfect

    rk <- .ranking_metrics(d, k_values)

    # Cluster bootstrap by ZONE (the correlated unit) for AUC-PR + log-score.
    #   * The log score is defined for EVERY replicate, including zero-event ones
    #     (-mean(log(1-p)) there). It used to be recorded only when a replicate contained
    #     >= 1 event, which silently DROPPED exactly the low-loss replicates and pushed the
    #     reported interval toward higher loss — a real bias, not a technicality.
    #   * AUC-PR is undefined without a positive (.auc_pr returns NA), so it is recorded
    #     only from replicates that contain one. That conditioning is unavoidable; it is now
    #     REPORTED (n_boot_aucpr) instead of being invisible in a shrunken quantile.
    # Zone -> row index is precomputed once: the previous which(d$health_zone == z) inside a
    # double loop was O(n_boot x n_zones x nrow(d)) (~1e9 comparisons per method x horizon).
    zones <- unique(d$health_zone)
    idx_by_zone <- split(seq_len(nrow(d)), factor(d$health_zone, levels = zones))
    boot_aucpr <- numeric(0); boot_ls <- numeric(0)
    if (n_pos >= 2 && length(zones) > 5) {
      for (b in seq_len(n_boot)) {
        zs  <- sample.int(length(zones), length(zones), replace = TRUE)
        idx <- unlist(idx_by_zone[zs], use.names = FALSE)
        yb  <- y[idx]; pb <- p[idx]
        if (.calib) boot_ls <- c(boot_ls, .log_score(pb, yb))   # always defined when p is one
        if (sum(yb) >= 1) boot_aucpr <- c(boot_aucpr, .auc_pr(pb, yb))
      }
    }
    n_boot_aucpr <- length(boot_aucpr); n_boot_ls <- length(boot_ls)

    # ---- OPTIONAL post-hoc recalibration (16b_invasion_recalibration.R) --------
    # Scored on EXACTLY the same rows `d` as everything above, so raw and recalibrated
    # numbers are directly comparable.
    #
    # READ auc_pr_recal / auc_roc_recal WITH CARE. The recalibration transform is
    # strictly monotone, so for a FIXED delta these are identical to the raw values.
    # The PREQUENTIAL delta, however, varies BY FOLD, and AUC-PR/AUC-ROC here are
    # computed on rows POOLED across folds — so they compare probabilities issued under
    # different factors and CAN differ from the raw values. That difference is a
    # fold-mixing artifact, not a change in discrimination: every metric computed per
    # fold and then averaged (mean_rank_of_truth, prec/recall/hit@K) is exactly
    # invariant, which is what the operator actually experiences within one forecast.
    # `recal_delta_varies` flags when the comparison is not like-for-like, and model
    # SELECTION always reads the raw auc_pr_skill, never these columns.
    # Recalibration maps a PROBABILITY through 1-(1-p)^delta. Applying it to a ranking is
    # meaningless for the same reason the proper scores are, so a rank-only model reports the
    # whole recalibrated block as NA rather than a delta fitted to an arbitrary scale.
    has_recal <- .calib && "p_recal" %in% names(d)
    p_rc <- if (has_recal) d$p_recal else NULL
    if (has_recal && !identical(is.finite(p_rc), is.finite(p))) {
      warning(sprintf(paste0("[eval] %s h=%s: p_recal is missing on rows where p_invasion is ",
                             "present; recalibrated metrics would be scored on a DIFFERENT ",
                             "support and are reported as NA instead."), meth, h), call. = FALSE)
      has_recal <- FALSE; p_rc <- NULL
    }
    ls_rc  <- if (has_recal) .log_score(p_rc, y)   else NA_real_
    bss_rc <- if (has_recal) .brier_skill(p_rc, y) else NA_real_
    ece_rc <- if (has_recal) .ece(p_rc, y)         else NA_real_
    cil_rc <- if (has_recal) mean(p_rc, na.rm = TRUE) / max(base_rate, 1e-9) else NA_real_
    aucpr_rc <- if (has_recal) .auc_pr(p_rc, y)  else NA_real_
    aucroc_rc <- if (has_recal) .auc_roc(p_rc, y) else NA_real_
    delta_rc <- if (has_recal && "delta_preq" %in% names(d))
                  suppressWarnings(stats::median(d$delta_preq, na.rm = TRUE)) else NA_real_
    delta_varies <- if (has_recal && "delta_preq" %in% names(d)) {
      dv <- unique(d$delta_preq[is.finite(d$delta_preq)])
      length(dv) > 1L
    } else NA
    delta_rc_last <- if (has_recal && all(c("delta_preq", "fold_id") %in% names(d))) {
      fl <- d$fold_id[is.finite(d$delta_preq)]
      if (length(fl)) d$delta_preq[is.finite(d$delta_preq)][which.max(fl)] else NA_real_
    } else NA_real_
    # 90% bootstrap interval (q5/q95) — the pipeline's uncertainty-band convention
    # for forecast/metric intervals (never 95%); the covariate-screen regression HRs
    # separately use standard 95% Wald CIs (see METHODS §5.5/§9.5).
    ci <- function(v) if (length(v) > 10) quantile(v, c(0.05, 0.95), na.rm = TRUE) else c(NA, NA)
    aucpr_ci <- ci(boot_aucpr); ls_ci <- ci(boot_ls)

    # n_invasions counts POSITIVE ROWS; n_invasion_events counts DISTINCT invaded zones.
    # They are equal at h=1 (disjoint windows) and differ at h>=2, where the cumulative
    # windows overlap and one invasion is scored from two origins (37 events -> 66 rows on
    # the shipped run). Publishing only the row count invites reading it as an effective
    # sample size; publishing both makes the dependence visible in the table itself.
    n_events_distinct <- length(unique(d$health_zone[d$is_new_invasion == 1L]))
    tibble::tibble(
      method = meth, horizon = h, n_atrisk = n, n_invasions = n_pos,
      n_invasion_events = n_events_distinct,
      base_rate = base_rate, n_folds = n_folds,
      # n_folds is this method's NATIVE reach; n_folds_scored is what the metrics in this row
      # were actually computed on, after restriction to the shared support. They differ
      # whenever a method covers folds that some other method does not, and reporting only
      # the first would suggest a baseline's 11-fold number came from 11 folds when it did
      # not. Both are published so the table cannot be misread either way.
      n_folds_scored = length(unique(d$fold_id)),
      coverage = n_folds / hs$maxf, partial_cv = partial_cv,
      # FALSE = this method emits a ranking, not a probability, so every calibration-dependent
      # column in this row is NA by construction and must not be read as a missing value.
      prob_calibrated = .calib,
      auc_pr = aucpr, auc_pr_lo = aucpr_ci[1], auc_pr_hi = aucpr_ci[2],
      auc_pr_skill = if (is.finite(aucpr)) aucpr / max(base_rate, 1e-9) else NA_real_,
      auc_roc = aucroc,
      log_score = lsc, log_score_lo = ls_ci[1], log_score_hi = ls_ci[2],
      # Retained bootstrap replicates. n_boot_aucpr < n_boot means AUC-PR replicates were
      # dropped for containing no event, so its interval is conditional on >= 1 event; the
      # log-score interval uses every replicate and is unconditional.
      n_boot_aucpr = n_boot_aucpr, n_boot_log_score = n_boot_ls,
      brier_skill = bss, ece = ece, calibration_in_large = cal_in_large,
      # Post-hoc recalibrated counterparts (NA when 16b was not applied). Purely
      # additive: no column above is altered by recalibration.
      log_score_recal = ls_rc, brier_skill_recal = bss_rc, ece_recal = ece_rc,
      calibration_in_large_recal = cil_rc,
      auc_pr_recal = aucpr_rc, auc_roc_recal = aucroc_rc,
      # Skill on the recalibrated scale, against the SAME base rate: the denominator is a
      # property of the outcome, not of the forecast, so it does not change with the
      # probability scale. Added so a figure drawn on the recalibrated probabilities can
      # quote a matching skill rather than borrowing the raw column.
      auc_pr_skill_recal = if (is.finite(aucpr_rc)) aucpr_rc / max(base_rate, 1e-9) else NA_real_,
      delta_recal_median = delta_rc, delta_recal_last = delta_rc_last,
      recal_delta_varies = delta_varies,
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
  # partial_cv models (scored on fewer cells than the shared support -- only reachable via
  # the fallback) sort last so they cannot outrank the comparable pack on a thin subset.
  out %>% dplyr::arrange(horizon, partial_cv, dplyr::desc(auc_pr_skill), mean_rank_of_truth)
}

#' Pick the best model POOLED across forecast horizons: composite-led, spike-guarded.
#'
#' The featured model drives the operational watch-list AND is the kernel the 13-week
#' cascade adopts (30_projection_config.R), so selection must reward all three qualities
#' that matter downstream, not discrimination alone. The criterion is the CV COMPOSITE:
#' within each horizon, models are ranked by (i) AUC-PR skill (higher better),
#' (ii) mean rank of the truly-invaded zones (lower better) and (iii) log score (lower
#' better; a proper score that penalises over-confidence); the three ranks are summed
#' within horizon and then summed across horizons, and the LOWEST total wins.
#'
#' WHY NOT AUC-PR SKILL ALONE (the rule this replaces): average precision is dominated by
#' the very top of the ranking, where near-front zones are ordered near-identically by every
#' mobility kernel — so it cannot separate a peaked kernel from a dispersed one (in the
#' 2026-08 run `M14-med` and `M17-med` differed by 0.5% in pooled skill). It is also blind to
#' calibration. Those two blind spots are precisely what the cascade is sensitive to: it
#' compounds the kernel row at every weekly step, so a dispersed kernel multiplies the
#' branching factor of the invasion tree, and a mis-calibrated hazard rescales delta. Adding
#' rank-of-truth and log score back into the OBJECTIVE (they were demoted to tie-breakers,
#' which continuous metrics never reach) restores the separation.
#'
#' The spiky-model GATE is retained ahead of the composite: any model whose worst-horizon
#' mean rank-of-truth exceeds `MRT_GATE_MULT` x the field median is dropped outright. Only
#' methods present at the maximal number of pooled horizons are eligible, so a model cannot
#' win by sitting out the harder horizon, and partial-CV methods (scored on an incomparable
#' support) are excluded. Both gates fall back to the ungated set if they would leave nothing.
#'
#' Ranks are computed over the SURVIVING eligible set, so the composite does not shift with
#' the presence of models that the gates already removed.
#'
#' @param horizons optional integer vector; if given, only these horizons are
#'   pooled (default: every horizon present in eval_tbl).
#' @param restrict optional regex; if given, only methods matching it are ranked
#'   (e.g. "^Renewal" to pick the best mobility-renewal model).
#' @param rule "composite" (default) or "aucpr" to recover the previous
#'   discrimination-only behaviour. Override globally with INVASION_SELECTION_RULE.
MRT_GATE_MULT <- 1.5   # spiky-model gate: drop methods with worst-horizon rank-of-truth > this x field median
# Selection objective. "composite" = summed within-horizon ranks of (AUC-PR skill,
# mean rank-of-truth, log score), pooled across horizons; "aucpr" = total AUC-PR skill only.
INVASION_SELECTION_RULE <- local({
  r <- tolower(Sys.getenv("INVASION_SELECTION_RULE", unset = "composite"))
  if (!r %in% c("composite", "aucpr")) {
    warning(sprintf("[best_model] unknown INVASION_SELECTION_RULE '%s'; using 'composite'.", r),
            call. = FALSE)
    r <- "composite"
  }
  r
})
#' Which log-score column the selection composite should use.
#'
#' Raw log score is dominated by calibration-in-the-large — a one-parameter defect that
#' 16b_invasion_recalibration.R removes — so with recalibration on, the log-score AXIS
#' measures refinement (does the model separate events from non-events) rather than how
#' badly its intercept was scaled. The other two composite axes are rank metrics and are
#' invariant to recalibration, so this is the ONLY channel through which recalibration
#' can move the featured pick.
#'
#' Fails loudly, never silently: an absent recalibrated column, or one missing on a row that
#' SHOULD carry it, reverts to the raw score with a warning. Rank-only rows are exempt from
#' that test -- their NA is structural. Shared by best_invasion_model() and
#' invasion_selection_table() so the audit table can never disagree with the pick.
.invasion_ls_col <- function(e, select_on_recal = INVASION_SELECT_ON_RECAL,
                             who = "best_model") {
  if (!isTRUE(select_on_recal)) return("log_score")
  if (!"log_score_recal" %in% names(e)) {
    warning(sprintf(paste0("[%s] INVASION_SELECT_ON_RECAL is TRUE but `log_score_recal` is ",
                           "absent; selecting on the RAW log score."), who), call. = FALSE)
    return("log_score")
  }
  # A RANK-ONLY ROW IS NA BY CONSTRUCTION, NOT BY FAILURE. A comparator that emits an ordering
  # rather than a probability (prob_calibrated = FALSE) has nothing to recalibrate, so 16b
  # never fits it a delta and this function's own caller NAs its whole calibration block. The
  # finiteness test used to run over those rows too, so a single surviving rank-only
  # comparator -- Distance-B1 clears the spiky-model gate on this run -- silently reverted the
  # ENTIRE composite to the raw log score. INVASION_SELECT_ON_RECAL was then TRUE and
  # inoperative: the full-table pick ran on the raw axis while the Bayes-only pick ran on the
  # recalibrated one, and on the shipped frame they chose DIFFERENT models
  # (Bayes-M14-fill-med vs Bayes-M14-fill-geo) -- which, through CASCADE_KERNEL, chooses the
  # 13-week projection's mobility kernel. Judge the axis on the rows that can carry it.
  .cal <- if ("prob_calibrated" %in% names(e)) !(e$prob_calibrated %in% FALSE) else rep(TRUE, nrow(e))
  if (!any(.cal)) {
    warning(sprintf(paste0("[%s] every candidate is rank-only, so no recalibrated log score ",
                           "exists; selecting on the RAW log score."), who), call. = FALSE)
    return("log_score")
  }
  if (any(!is.finite(e$log_score_recal[.cal]))) {
    warning(sprintf(paste0("[%s] INVASION_SELECT_ON_RECAL is TRUE but %d of %d CALIBRATED rows ",
                           "have no recalibrated log score; selecting on the RAW log score."),
                    who, sum(!is.finite(e$log_score_recal[.cal])), sum(.cal)), call. = FALSE)
    return("log_score")
  }
  "log_score_recal"
}

#' @param exclude regex of method labels that are SCORED but are not CANDIDATES, built in
#'   00_config.R from SENSITIVITY_ARM_SUFFIXES. Two families are covered, for the same
#'   reason: both are REFITS of a model already in the grid rather than competing
#'   hypotheses, so letting one win would make the candidate set depend on the previous
#'   winner. The generation-time arms (`-gtshort` / `-gtlong`) would additionally make the
#'   GT a selection axis, which composing the grid at a single anchor is designed to
#'   prevent. The time-varying arms (`-tv*`) were eligible until 2026-09-23; the held-out
#'   optimism check, which selects from this same pool, picked one on two consecutive runs
#'   and reported a gap of +5.1 skill points (+22.3 on a common base rate) on the second.
#'   All of them must still be cross-validated — that IS the sensitivity — so they cannot
#'   simply be left out of the suite. Applied AFTER `restrict`. Set NULL to make every
#'   scored model a candidate.
best_invasion_model <- function(eval_tbl, horizons = NULL, restrict = NULL,
                                rule = INVASION_SELECTION_RULE,
                                select_on_recal = INVASION_SELECT_ON_RECAL,
                                exclude = get0("INVASION_SELECTION_EXCLUDE",
                                               ifnotfound = "-(gtshort|gtlong|tv[a-z0-9]+)$"),
                                require_calibrated = TRUE) {
  e <- eval_tbl
  if (!is.null(exclude) && length(exclude) == 1L && !is.na(exclude) && nzchar(exclude)) {
    .drop <- grepl(exclude, e$method)
    if (any(.drop)) {
      message(sprintf("[best_model] %d sensitivity-arm method(s) scored but not eligible: %s",
                      dplyr::n_distinct(e$method[.drop]),
                      paste(sort(unique(e$method[.drop])), collapse = ", ")))
      e <- e[!.drop, , drop = FALSE]
    }
  }
  # CANDIDACY: a method that does not produce PROBABILITIES cannot be featured. best_method
  # drives the headline figures, DISPLAY_MODELS and plot_predobs_over_folds(), all of which
  # read p_invasion AS a probability -- a rank-only comparator (the structural baselines,
  # whose score is a connectivity mass or an inverse travel time) would have its score
  # plotted and reported as P(first case). Until the structural baselines entered the
  # evaluation table this was academic, because the rank-only comparators in it were weak;
  # a baseline that matches the model on AUC-PR makes it live. Being ineligible to be
  # FEATURED is not the same as being unscored: they are ranked in the table like everyone
  # else. The gate never empties the frame.
  if (isTRUE(require_calibrated) && "prob_calibrated" %in% names(e)) {
    .rank_only <- !(e$prob_calibrated %in% TRUE)
    if (any(.rank_only) && any(!.rank_only)) {
      message(sprintf("[best_model] %d rank-only comparator(s) scored but not eligible to be featured: %s",
                      dplyr::n_distinct(e$method[.rank_only]),
                      paste(sort(unique(e$method[.rank_only])), collapse = ", ")))
      e <- e[!.rank_only, , drop = FALSE]
    }
  }
  if (!is.null(restrict)) {
    # STRICT: if a family filter is given and nothing matches, return NA rather
    # than silently falling back to the unrestricted set (which would return an
    # out-of-family model, e.g. a frequentist model for a Bayesian request).
    # Callers supply their own defaults.
    e <- e %>% dplyr::filter(grepl(restrict, method))
  }
  if (!is.null(horizons)) e <- e %>% dplyr::filter(horizon %in% !!horizons)
  if (nrow(e) == 0) return(NA_character_)
  # WHICH log score feeds the composite. Raw log score is dominated by calibration-in-
  # the-large — a one-parameter defect that 16b removes — so with recalibration on, the
  # log-score AXIS measures refinement (does the model separate events from non-events)
  # rather than how badly its intercept was scaled. The other two composite axes are
  # rank metrics and are invariant to recalibration, so this is the ONLY channel through
  # which recalibration can move the featured pick. Falls back loudly, never silently.
  # The log-score column is chosen from the FINAL eligible set, below — NOT here. Resolving it
  # at this point (after only the restrict/horizon filters) was wrong: the partial-CV gate and
  # the spiky-MRT gate BOTH remove further methods afterwards, and .invasion_ls_col() reverts
  # the whole composite to the raw log score if ANY row it is shown is non-finite. So a single
  # method that was about to be gated out could flip the scoring axis for every method that
  # survived — the exact "spurious fallback" the old comment here claimed was impossible.
  # Both candidate sums are carried through the aggregation and the choice is made once the
  # eligible set is final.
  agg <- e %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(n_h     = dplyr::n_distinct(horizon),
                     aps     = sum(dplyr::coalesce(auc_pr_skill, 0)),      # total AUC-PR skill (lead)
                     mrt_max = max(dplyr::coalesce(mean_rank_of_truth, Inf)),  # worst-horizon targeting
                     mrt_sum = sum(dplyr::coalesce(mean_rank_of_truth, Inf)),
                     ls_sum_raw   = sum(dplyr::coalesce(log_score, Inf)),
                     ls_sum_recal = if ("log_score_recal" %in% names(e))
                       sum(dplyr::coalesce(log_score_recal, Inf)) else NA_real_,
                     # A method is partial if it is partial at ANY horizon.
                     any_partial = if ("partial_cv" %in% names(e))
                       any(partial_cv %in% TRUE) else FALSE,
                     .groups = "drop") %>%
    dplyr::filter(n_h == max(n_h))
  # COMPARABILITY gate. evaluate_invasion() scores every method on ONE shared (fold x zone)
  # support and sets `partial_cv` when a method still came up short of it, precisely so a
  # model is not compared on a different, usually base-rate-ENRICHED subset -- it would be
  # "mis-ranked purely through its denominator". It sorts those models last in its own output.
  # But this function -- the one that actually PICKS the featured model -- filtered only on
  # horizon coverage (n_h) and ignored the flag, so a model scored on a handful of event-rich
  # cells could take the crown on an incomparable denominator. Prefer comparable methods; fall
  # back to the flagged set only if nothing is comparable, and say so rather than silently.
  # NOTE this gate is about the SCORED support, not native fold reach: a family that covers
  # fewer folds by design (the rolling-predictor floor) is still scored on the shared cells
  # and is not excluded here. See the `partial_cv` note in evaluate_invasion().
  if (any(!agg$any_partial)) {
    if (any(agg$any_partial))
      message(sprintf("[best_model] excluding %d partial-CV method(s) from selection: %s",
                      sum(agg$any_partial), paste(agg$method[agg$any_partial], collapse = ", ")))
    agg <- agg[!agg$any_partial, , drop = FALSE]
  } else if (nrow(agg) > 0) {
    warning("[best_model] EVERY candidate is partial-CV; selecting among incomparable ",
            "supports (skill denominators differ across methods).", call. = FALSE)
  }
  # spiky-model gate (relative to the field), with graceful fallback to the ungated set.
  finite_mrt <- agg$mrt_max[is.finite(agg$mrt_max)]
  gate <- if (length(finite_mrt)) MRT_GATE_MULT * stats::median(finite_mrt) else Inf
  keep <- agg %>% dplyr::filter(mrt_max <= gate)
  if (nrow(keep) == 0) keep <- agg
  # NOW the eligible set is final: choose the log-score axis on exactly the rows that will be
  # ranked, so a gated-out method can no longer force the fallback for the survivors.
  ls_col <- .invasion_ls_col(e[e$method %in% keep$method, , drop = FALSE], select_on_recal)
  keep$ls_sum <- if (identical(ls_col, "log_score_recal")) keep$ls_sum_recal else keep$ls_sum_raw
  if (identical(rule, "aucpr")) {
    return(keep %>% dplyr::arrange(dplyr::desc(aps), mrt_sum, ls_sum) %>%
             dplyr::slice(1) %>% dplyr::pull(method))
  }
  # CV composite: rank WITHIN horizon over the surviving eligible set, sum the three
  # ranks, then sum across horizons. NA metrics rank last on their own axis (never
  # silently best) via the +/-Inf coalesce.
  comp <- e %>%
    dplyr::filter(method %in% keep$method) %>%
    dplyr::group_by(horizon) %>%
    dplyr::mutate(
      r_aps = rank(-dplyr::coalesce(auc_pr_skill, -Inf),      ties.method = "min"),
      r_mrt = rank( dplyr::coalesce(mean_rank_of_truth, Inf), ties.method = "min"),
      r_ls  = rank( dplyr::coalesce(.data[[ls_col]], Inf),    ties.method = "min")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(composite = sum(r_aps + r_mrt + r_ls), .groups = "drop")
  keep %>%
    dplyr::left_join(comp, by = "method") %>%
    # composite first; the old discrimination lead survives only as a tie-break
    dplyr::arrange(composite, dplyr::desc(aps), mrt_sum, ls_sum) %>%
    dplyr::slice(1) %>% dplyr::pull(method)
}

#' The full CV-composite leaderboard (same computation as best_invasion_model()).
#'
#' Exposed so the selection is auditable: `write_model_selection.R` and the reports can
#' show WHY a model won and by how much, rather than only naming the winner. Returns one
#' row per eligible method with the per-horizon rank triplets and the pooled composite.
#'
#' Applies the SAME eligibility filters as best_invasion_model() — max horizon coverage,
#' no partial-CV, and the spiky-model gate — and ranks over the surviving pool, so the
#' composite printed here is the one the selection actually used. (Ranking over the full
#' field instead would shift ranks by however many gated models sat above each survivor
#' on each axis, which is not guaranteed to preserve the composite ordering.)
invasion_selection_table <- function(eval_tbl, horizons = NULL, restrict = NULL,
                                    select_on_recal = INVASION_SELECT_ON_RECAL,
                                    exclude = get0("INVASION_SELECTION_EXCLUDE",
                                                   ifnotfound = "-(gtshort|gtlong|tv[a-z0-9]+)$")) {
  e <- eval_tbl
  # Same non-candidate gate as best_invasion_model(), for the same reason: an audit table
  # that ranks a model the selection could never pick describes a different competition
  # from the one that was held.
  if (!is.null(exclude) && length(exclude) == 1L && !is.na(exclude) && nzchar(exclude))
    e <- e[!grepl(exclude, e$method), , drop = FALSE]
  if (!is.null(restrict)) e <- e %>% dplyr::filter(grepl(restrict, method))
  if (!is.null(horizons)) e <- e %>% dplyr::filter(horizon %in% !!horizons)
  if (nrow(e) == 0) return(NULL)
  if ("partial_cv" %in% names(e)) {
    full <- e %>% dplyr::group_by(method) %>%
      dplyr::summarise(any_partial = any(partial_cv %in% TRUE), .groups = "drop")
    if (any(!full$any_partial))
      e <- e %>% dplyr::filter(method %in% full$method[!full$any_partial])
  }
  nh <- e %>% dplyr::group_by(method) %>%
    dplyr::summarise(n_h = dplyr::n_distinct(horizon), .groups = "drop")
  e <- e %>% dplyr::filter(method %in% nh$method[nh$n_h == max(nh$n_h)])
  # spiky-model gate, identical to best_invasion_model()
  mrt <- e %>% dplyr::group_by(method) %>%
    dplyr::summarise(mrt_max = max(dplyr::coalesce(mean_rank_of_truth, Inf)), .groups = "drop")
  fin <- mrt$mrt_max[is.finite(mrt$mrt_max)]
  gate <- if (length(fin)) MRT_GATE_MULT * stats::median(fin) else Inf
  surv <- mrt$method[mrt$mrt_max <= gate]
  if (length(surv)) e <- e %>% dplyr::filter(method %in% surv)
  # Resolved AFTER EVERY eligibility filter — the partial-CV gate, the horizon-coverage gate
  # and the spiky-MRT gate — and with the SAME helper best_invasion_model() uses, so the audit
  # table and the pick rank on the same axis. Resolving it before those gates (as this did)
  # let one method that was about to be removed flip the axis for all the survivors, and
  # because best_invasion_model() had the identical ordering bug the two flipped in lockstep
  # and the disagreement was invisible.
  ls_col <- .invasion_ls_col(e, select_on_recal, who = "selection_table")
  e %>%
    dplyr::group_by(horizon) %>%
    dplyr::mutate(
      r_aps = rank(-dplyr::coalesce(auc_pr_skill, -Inf),      ties.method = "min"),
      r_mrt = rank( dplyr::coalesce(mean_rank_of_truth, Inf), ties.method = "min"),
      r_ls  = rank( dplyr::coalesce(.data[[ls_col]], Inf),    ties.method = "min"),
      r_sum = r_aps + r_mrt + r_ls) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(log_score_used = .data[[ls_col]], log_score_axis = ls_col) %>%
    dplyr::select(method, horizon, auc_pr_skill, mean_rank_of_truth, log_score,
                  log_score_used, log_score_axis,
                  r_aps, r_mrt, r_ls, r_sum) %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(composite = sum(r_sum)) %>%
    dplyr::ungroup() %>%
    # SAME TIE-BREAK AS best_invasion_model(), so the audit table's top row IS the model the
    # pipeline features. Ordering by `method` alphabetically broke ties differently: on the
    # shipped run Bayes-M14-fill-geo and Bayes-M14-fill-med tie at composite 39, the table
    # listed fill-geo first and the pipeline featured fill-med — while this function's own
    # docstring promises "the audit table can never disagree with the pick".
    # The per-method tie-break keys are recomputed here exactly as best_invasion_model()
    # builds them (summed over horizons), then dropped again so the returned columns are
    # unchanged.
    dplyr::group_by(method) %>%
    dplyr::mutate(.aps_sum = sum(dplyr::coalesce(auc_pr_skill, 0)),
                  .mrt_sum = sum(dplyr::coalesce(mean_rank_of_truth, Inf)),
                  .ls_sum  = sum(dplyr::coalesce(log_score_used, Inf))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(composite, dplyr::desc(.aps_sum), .mrt_sum, .ls_sum, method, horizon) %>%
    dplyr::select(-.aps_sum, -.mrt_sum, -.ls_sum)
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
evaluate_invasion_heldout <- function(lfo_results, n_holdout = 3L,
                                      restrict = "^Bayes", horizon = 1L) {
  if (is.null(lfo_results) || !nrow(lfo_results) || !"cutoff" %in% names(lfo_results)) return(NULL)
  cutoffs <- sort(unique(lfo_results$cutoff))
  if (length(cutoffs) < n_holdout + 2L) {
    warning("[heldout] too few origins (", length(cutoffs), ") for a ", n_holdout,
            "-origin holdout; skipping"); return(NULL)
  }
  # THE HELD-OUT SET MUST BE EVALUABLE. This used to be the last `n_holdout` origins, full
  # stop. Since the cross-validation runs close to the analysis date, the final rounds can
  # legitimately contain NO invasion at all -- that was the case when the admission gate was
  # -1 and the grid reached the weeks beginning 2026-08-18 and 2026-08-25; under the current
  # 13-day gate every admitted round carries at least one -- and every discrimination metric is
  # undefined without a positive — so the whole optimism check would silently collapse to
  # NA. Take the last `n_holdout` origins that actually carry an event at the scored
  # horizon. Event-free rounds later than those belong to neither set: they carry no
  # information about discrimination, and putting them in the inner set would not change
  # the selection either.
  .hz <- lfo_results[lfo_results$horizon %in% horizon, , drop = FALSE]
  if ("was_active_before" %in% names(.hz))
    .hz <- .hz[!(as.logical(.hz$was_active_before) %in% TRUE), , drop = FALSE]
  .ev_by_cut <- tapply(as.integer(.hz$is_new_invasion), as.character(.hz$cutoff),
                       function(v) sum(v, na.rm = TRUE) > 0)
  .evaluable <- cutoffs[as.character(cutoffs) %in% names(.ev_by_cut)[which(.ev_by_cut)]]
  if (length(.evaluable) < n_holdout) {
    warning(sprintf("[heldout] only %d origin(s) at h=%s carry an invasion; a %d-origin holdout is not possible.",
                    length(.evaluable), paste(horizon, collapse = "/"), n_holdout),
            call. = FALSE)
    return(NULL)
  }
  outer_cut <- utils::tail(.evaluable, n_holdout)
  if (!identical(as.character(outer_cut), as.character(utils::tail(cutoffs, n_holdout))))
    message(sprintf("[heldout] held-out origins are the last %d EVALUABLE rounds (%s); %d later round(s) carry no invasion at h=%s and are used by neither set.",
                    n_holdout, paste(format(outer_cut), collapse = ", "),
                    sum(cutoffs > max(outer_cut)), paste(horizon, collapse = "/")))
  # EMBARGO between the selection set and the held-out set. Splitting on cutoff alone does NOT
  # make the outer set untouched: an inner fold's outcome window extends h weeks past its own
  # cutoff, so the last inner folds' windows reach into the outer origins' windows and the two
  # sets share EVENTS. Measured on the shipped folds, 2 of the outer set's 4 h=1 events were
  # also positives in the inner set — while this function reports the result as "touched only
  # once". Purge the inner folds whose outcome window can reach the first outer origin
  # (standard purged/embargoed CV). `hmax` is taken from the data, not assumed.
  .hmax <- suppressWarnings(max(as.integer(lfo_results$horizon), na.rm = TRUE))
  if (!is.finite(.hmax) || .hmax < 1L) .hmax <- 1L
  .embargo_days <- 7L * .hmax
  # Only rounds BEFORE the held-out ones can be inner: a round later than the held-out
  # origins is in the future relative to them, and selecting on it would invert the
  # causality this check exists to respect.
  inner_cut <- setdiff(cutoffs[cutoffs < min(outer_cut)], outer_cut)
  .first_outer <- min(outer_cut)
  .kept <- inner_cut[inner_cut + .embargo_days <= .first_outer]
  if (length(.kept) >= 2L) {
    .n_purged <- length(inner_cut) - length(.kept)
    if (.n_purged > 0L)
      message(sprintf("[heldout] embargo: purged %d inner fold(s) whose %d-week outcome window reached the held-out origins; %d inner folds remain.",
                      .n_purged, .hmax, length(.kept)))
    inner_cut <- .kept
  } else {
    warning(sprintf("[heldout] the %d-week embargo would leave only %d inner fold(s); it was NOT applied, so the inner and held-out sets SHARE events and the optimism gap is understated.",
                    .hmax, length(.kept)), call. = FALSE)
  }
  inner_res <- lfo_results[lfo_results$cutoff %in% inner_cut, , drop = FALSE]
  outer_res <- lfo_results[lfo_results$cutoff %in% outer_cut, , drop = FALSE]
  inner_eval <- tryCatch(evaluate_invasion(inner_res), error = function(e) NULL)
  outer_eval <- tryCatch(evaluate_invasion(outer_res), error = function(e) NULL)
  if (is.null(inner_eval) || is.null(outer_eval)) return(NULL)
  # Select on the INNER folds only (never on the held-out outer folds).
  sel_tbl <- if (!is.null(restrict)) dplyr::filter(inner_eval, grepl(restrict, method)) else inner_eval
  sel <- best_invasion_model(sel_tbl)
  if (length(sel) != 1L || is.na(sel)) return(NULL)
  gk <- function(ev, m, col = "auc_pr_skill") {
    r <- ev[[col]][ev$method == m & ev$horizon == horizon]
    if (length(r)) r[1] else NA_real_
  }
  inner_skill <- gk(inner_eval, sel); outer_skill <- gk(outer_eval, sel)
  # BASE-RATE DRIFT. inner_eval and outer_eval are two independent evaluate_invasion() calls,
  # so each computes its OWN base_rate, and auc_pr_skill = auc_pr / base_rate. The published
  # `optimism_gap` is therefore a difference of ratios with DIFFERENT denominators: it mixes
  # genuine selection optimism with the fact that the outer window is simply a harder or easier
  # problem. Measured 2026-09-19 on the shipped folds (Bayes-M14-fill-med, h=1): inner base rate
  # 0.006837 (AUC-PR 0.3221, skill 47.12), outer 0.004301 (AUC-PR 0.1485, skill 34.52) — a 37%
  # denominator shift. The reported gap is 12.60; on a COMMON denominator it is 25.40, i.e. the
  # headline number understates the optimism roughly twofold.
  # Both are now returned, with the inputs, so the reader can see which is which.
  inner_base <- gk(inner_eval, sel, "base_rate"); outer_base <- gk(outer_eval, sel, "base_rate")
  inner_ap   <- gk(inner_eval, sel, "auc_pr");    outer_ap   <- gk(outer_eval, sel, "auc_pr")
  list(selected = sel, horizon = horizon, n_holdout = n_holdout,
       inner_skill = inner_skill, outer_skill = outer_skill,
       # As published: each set scored against its own base rate.
       optimism_gap = inner_skill - outer_skill,
       # On the INNER base rate throughout, so the denominator cannot move the answer.
       optimism_gap_common_base = if (is.finite(inner_base) && inner_base > 0)
         (inner_ap - outer_ap) / inner_base else NA_real_,
       inner_base_rate = inner_base, outer_base_rate = outer_base,
       inner_auc_pr = inner_ap, outer_auc_pr = outer_ap,
       base_rate_ratio = if (is.finite(inner_base) && inner_base > 0)
         outer_base / inner_base else NA_real_,
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
#'   - outcome(h)   = zone is first confirmed in weeks C+1..C+h (cumulative window,
#'     matching the models' cumulative invasion probability).
#'
#'     OVERLAPPING WINDOWS AT h >= 2 — read this before quoting an h=2 number. With weekly
#'     cutoffs the h=1 windows are disjoint, so each invasion is scored exactly once: on the
#'     shipped run, 37 positives from 37 distinct zones. The h=2 windows OVERLAP by a week,
#'     so a zone invaded in week W is still at-risk at both C = W-14 and C = W-7 and is a
#'     positive in both: 66 positives from only 37 distinct invasions, 29 zones counted
#'     twice (e.g. Adja in the 2026-06-23 and 2026-06-30 folds). This docstring previously
#'     asserted the opposite ("no double counting"), which is false.
#'
#'     Each ROW is still a legitimate (forecast, outcome) pair — "will this zone be invaded
#'     within 2 weeks, as of C?" asked at two different origins — and every metric is
#'     internally consistent on those rows, including the h=2 base rate, which is the rate
#'     among the rows actually scored. What is NOT true at h=2 is independence: 66 rows carry
#'     37 events' worth of information. Interval estimates are protected (the bootstrap
#'     resamples ZONES, the correlated unit); any reading that treats n_invasions at h=2 as
#'     an effective sample size is not. `n_invasion_events` is reported alongside
#'     `n_invasions` so the two cannot be confused.
#' Folds whose farthest-horizon eval week is not yet ~complete (delay-truncated)
#' are dropped so the ground truth is reliable.
#'
#' @param zone_week_raw raw zone-week counts (health_zone, week_start, confirmed,
#'   suspected), outbreak period.
#' @param models named list of forecasters, each `function(zw_nc, t_idx,
#'   horizons, cutoff)` returning health_zone, horizon, mu_forecast, p_invasion,
#'   method.
#' @param nowcast_fn function(zone_week, analysis_date) -> nowcast-corrected;
#'   default apply_nowcast_correction (must be sourced).
#' @param linelist optional line list (dat$ll). When supplied, per-fold training
#'   counts are reconstructed leakage-free via reaggregate_asof() (22_daily_reissue.R)
#'   instead of slicing the final weekly counts. Recommended for a truly real-time
#'   evaluation; NULL reproduces the previous (mildly leaky) behaviour.
#' @return pooled tibble across folds×models with p_invasion, is_new_invasion,
#'   was_active_before, fold_id, cutoff, method, horizon.
#' @param delay the truncation spec used by the per-fold nowcast — in the pipeline the fitted
#'   AS-OF truncation (epinow2_truncation(regime = "asof"), passed by run_all.R). One spec,
#'   fitted on the FULL snapshot, shared by every fold and every model; see the "SHARED
#'   NUISANCE PARAMETER" note in the body. Default NULL falls through to whatever `nowcast_fn`
#'   resolves, which is the onset->SAMPLE delay — the wrong leg for a truncation correction —
#'   and warns.
run_invasion_lfo <- function(zone_week_raw, models, horizons = LFO_HORIZONS,
                             min_train_weeks = 3L, min_atrisk_events = 1L,
                             # ADMISSION THRESHOLD, in days: a fold is admitted at horizon h
                             # when analysis_date >= (cut + 7h + 7) + min_eval_age_days. Truth
                             # comes from first_case_week on the RAW (un-nowcast) counts
                             # bucketed by ONSET, and this run's fitted delays put onset->swab
                             # at p90 = 20 d plus swab->lab p90 = 2 d, so a young outcome week
                             # is RIGHT-TRUNCATED: a zone invaded late in it can be scored
                             # is_new_invasion = 0 because its index case is not yet confirmed.
                             #
                             # This used to be a hard-coded 25 d, which stopped the CV five
                             # weeks short of the deployed forecast. The value now comes from
                             # 00_config.R (LFO_MIN_EVAL_AGE_DAYS, default -1 = admit every
                             # cutoff whose outcome window closes by the analysis date), and
                             # the truncation is handled where it belongs — reported rather
                             # than avoided. Every returned row carries `eval_age_days` (the
                             # age of ITS OWN outcome window, per horizon) and
                             # `eval_reliable` (age >= LFO_RELIABLE_EVAL_AGE_DAYS), so any
                             # consumer can restrict to settled folds; run_all.R publishes
                             # exactly that restriction as a sensitivity table.
                             #
                             # Truncation biases the newest folds' skill DOWNWARD and their
                             # calibration-in-the-large UPWARD. It cannot bias the ranking
                             # between models: all are scored on identical rows and truth.
                             min_eval_age_days = get0("LFO_MIN_EVAL_AGE_DAYS", ifnotfound = 25L),
                             # Outcome-window age at which a fold's invasion count is treated
                             # as settled. Purely a LABEL: it never admits or excludes a fold.
                             reliable_eval_age_days = get0("LFO_RELIABLE_EVAL_AGE_DAYS",
                                                           ifnotfound = 25L),
                             analysis_date = ANALYSIS_DATE,
                             nowcast_fn = NULL, linelist = NULL, delay = NULL) {
  if (is.null(nowcast_fn)) nowcast_fn <- get0("apply_nowcast_correction")
  stopifnot(!is.null(nowcast_fn))
  # ---- SHARED NUISANCE PARAMETER: the onset->sample delay ---------------------------
  # Every fold's nowcast uses ONE delay, fitted on the whole snapshot. That is deliberate,
  # and it is NOT the same thing as the training-count leak that reaggregate_asof() closes:
  #   * A training-count leak lets a MODEL see its own future, inflating that model's skill
  #     and biasing the RANKING. That is what this evaluation exists to prevent.
  #   * The delay is applied identically to every model when `train_nc` is built, so it
  #     cannot move any model relative to another. It sits in the same class as the GT
  #     profile, the mobility kernels, the population spine and the ascertainment grid —
  #     all fixed across folds, all informed by the full record or by external sources.
  # Refitting it per fold was measured and is WORSE. A censored MLE on data available at
  # the cutoff carries the very right-truncation bias the delay model exists to remove, and
  # the bias SHRINKS as folds accrue data (fitted mean 3.41 d at the 2026-05-16 cutoff,
  # 4.71 d at 06-06, 6.48 d at 07-25, against a truncation-corrected 8.00 d). That turns a
  # constant shared nuisance into one CORRELATED WITH FOLD INDEX: early folds would be
  # under-corrected by ~34% and late folds by ~10%, so apparent skill would drift with
  # calendar time for a purely artefactual reason. The only per-fold estimator that would
  # not do this is a truncation-corrected EpiDist fit at each cutoff — a Stan fit per fold,
  # which this evaluation already declines to do for epinowcast, for the same cost reason.
  # CONSEQUENCE, stated plainly: absolute skill here is a RETROSPECTIVE estimate (a genuine
  # real-time system would not have known this delay), while model rankings are unaffected.
  # To run the sensitivity, pass `delay =` a spec (see effective_onset_sample_delay()).
  .delay_used <- delay
  if (is.null(.delay_used) && exists("effective_onset_sample_delay", mode = "function"))
    .delay_used <- tryCatch(effective_onset_sample_delay(), error = function(e) NULL)
  if (!is.null(.delay_used))
    message(sprintf("[invasion_lfo] Shared truncation for every fold: %s (mean %.2f d, %s)",
                    .delay_used$family, .delay_used$mean, .delay_used$estimator))
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

  # Candidate cutoffs. WHICH HORIZON GOVERNS ADMISSION is set by LFO_COMMON_FOLD_GRID
  # (00_config.R, which carries the full rationale and the reason no single age threshold can
  # align the two horizons' grids):
  #   TRUE  (default) — the FARTHEST horizon must be complete, so every admitted cutoff is
  #                     scored at EVERY horizon and the horizons share one set of origins.
  #   FALSE           — the Q1 fold augmentation: the SHORTEST horizon governs, adding the
  #                     recent cutoffs where h=1 is evaluable and h=2 is not.
  # The per-horizon `complete_h` guard below is kept under BOTH settings: it is what makes the
  # augmentation leak-free, and under the common grid it is simply always satisfied.
  .h_admit <- if (isTRUE(get0("LFO_COMMON_FOLD_GRID", ifnotfound = TRUE)))
                max(horizons) else min(horizons)
  cand <- weeks[weeks >= weeks[min_train_weeks] &
                  (as.numeric(analysis_date - (weeks + .h_admit * 7 + 7)) >= min_eval_age_days)]
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
    # DROPPING EVENT-FREE FOLDS IS SELECTION ON THE OUTCOME. With min_atrisk_events = 1 a
    # round in which nothing happened is deleted, which (a) enriches the pooled base rate,
    # making `calibration_in_large` read as less over-predicting than it is, and (b) puts a
    # hole in every over-rounds panel at exactly the quiet rounds. The pipeline now passes
    # min_atrisk_events = 0 (run_all.R) so every complete round is scored, including the
    # quiet ones — on the 2026-09-07 frame no zone is first confirmed in the weeks beginning
    # 2026-08-25 or 2026-09-01, so those rounds exist ONLY under 0.
    #
    # WHAT AN EVENT-FREE FOLD CONTRIBUTES. Rows, and therefore an honest denominator for the
    # pooled base rate, Brier, log score and calibration-in-the-large. It contributes NOTHING
    # to any per-fold discrimination metric: AUC-PR and AUC-ROC are undefined without a
    # positive, and the per-fold scorers return NA for it (19_spacetime_eval.R), which is why
    # skill-over-time correctly shows a gap rather than a spurious zero at such a round.
    # min_atrisk_events = 1 restores the previous behaviour exactly.
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
  # Say OUT LOUD how many folds carry a not-yet-settled outcome window, and which they are.
  # The admission threshold is a config constant now, so the one place a reader is guaranteed
  # to see the consequence is the run log.
  local({
    .prov <- vapply(folds, function(f)
      any(as.numeric(analysis_date - (f$cutoff + f$horizons * 7 + 7)) < reliable_eval_age_days),
      logical(1))
    if (any(.prov))
      message(sprintf(paste0("[invasion_lfo] %d of %d folds have an outcome window younger than ",
                             "%d d and are FLAGGED PROVISIONAL (eval_reliable = FALSE): %s. ",
                             "Their invasion counts may still rise as cases are confirmed."),
                      sum(.prov), length(folds), reliable_eval_age_days,
                      paste(vapply(folds[.prov], function(f) format(f$cutoff, "%m-%d"), ""),
                            collapse = ",")))
  })

  out <- list()
  for (fi in seq_along(folds)) {
    f <- folds[[fi]]; cut <- f$cutoff
    # Leakage-free training counts: reconstruct what surveillance actually KNEW at
    # the real-time forecast moment (cut+6) by re-aggregating the line list censored
    # to that observation date (reaggregate_asof, 22_daily_reissue.R), rather than
    # slicing the FINAL onset-bucketed counts — which fold in cases only reported
    # after cut+6 and would then be further inflated by the nowcast (the training-side
    # revision leak). Falls back to the final-count slice when no line list is given
    # or the re-aggregation fails, so behaviour degrades gracefully rather than erroring.
    train_raw <- NULL
    if (!is.null(linelist) && exists("reaggregate_asof")) {
      train_raw <- tryCatch(
        reaggregate_asof(linelist, zones, cut + 6, week_spine = weeks[weeks <= cut]),
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
    # Pass the shared delay explicitly when one was supplied, so a sensitivity run
    # (`delay =`) actually reaches the fold nowcast instead of being silently ignored.
    # FOLD ORIGIN = cut + 6, the LAST DAY OF THE LAST TRAINING WEEK — mirroring deployment,
    # where the grid is anchored so the final week ENDS on ANALYSIS_DATE and the as-of date IS
    # that day. The origin used to be cut + 7, which handed each fold ONE EXTRA DAY of
    # reporting for its most recent week: with the fitted onset->sample delay the deployed
    # current week is multiplied by 2.631x while the folds saw 2.185x, a 20.4% gap. Since the
    # recalibration delta is FITTED on the folds and APPLIED to the deployed forecast, delta
    # could not absorb a regime it never saw, and every published invasion probability carried
    # the residual. The horizon targets are unchanged: they are defined relative to `cut`.
    train_nc  <- tryCatch(if (is.null(delay)) nowcast_fn(train_raw, analysis_date = cut + 6)
                          else nowcast_fn(train_raw, analysis_date = cut + 6, delay = delay),
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
        # prob_calibrated travels with the forecast: a model that declares FALSE is emitting
        # a RANKING, not a probability, and its calibration-dependent scores are suppressed in
        # evaluate_invasion(). Absent = TRUE (every Bayesian model has a fitted hazard scale).
        dplyr::select(dplyr::any_of(c("health_zone", "horizon", "mu_forecast",
                                      "p_invasion", "method", "prob_calibrated"))) %>%
        dplyr::inner_join(truth_h, by = c("health_zone", "horizon")) %>%
        dplyr::mutate(fold_id = fi, cutoff = cut, was_active_before = FALSE,
                      # HOW SETTLED IS THIS ROW'S GROUND TRUTH? Per-HORIZON, because the h=2
                      # outcome window closes a week later than the h=1 window at the same
                      # cutoff and is correspondingly less complete. Attached to every row so
                      # no consumer has to re-derive it from `cutoff` and guess the
                      # convention. `eval_reliable` is a label, never a filter: the rows are
                      # scored either way (see min_eval_age_days).
                      eval_age_days = as.numeric(analysis_date - (cut + horizon * 7 + 7)),
                      eval_reliable = eval_age_days >= reliable_eval_age_days)
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
