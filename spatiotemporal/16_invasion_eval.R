# =============================================================================
# 16_invasion_eval.R — Invasion-Focused Evaluation (pooled across folds)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# The prediction task is a RARE BINARY event: does a zone with no confirmed cases
# up to the cutoff see its first case within the next h weeks? Each forecast date has only a
# few invasions, so probability metrics (log score, Brier, calibration) are computed on the
# POOLED at-risk rows. AUC-PR is reported both POOLED and PER DATE (average precision within
# each forecast list, weighted by that date's invasions; skill = sum_d E_d AP_d / sum_d E_d pi_d,
# pi_d the date's prevalence -- see .auc_pr_by_date()):
# the per-date version is the SELECTION axis, because a pooled ranking also rewards differences
# in probability level between dates. Ranking metrics (Precision@K, rank-of-truth) are computed
# per fold then averaged (ranking is only meaningful within a single forecast).
#
# Metrics reported per (method, horizon):
#   PRIMARY : AUC-PR (+prevalence baseline), log-score (raw), mean rank-of-truth,
#             Precision@K / Recall@K / hit-rate@K, calibration-in-the-large
#   SECONDARY: AUC-ROC, Brier skill vs base rate, ECE
# Cluster-bootstrap 90% CI (resampling by zone) on AUC-PR, log score, AUC-PR skill, the per-date
# AUC-PR skill and the selection-axis (post-burn-in) per-date skill.
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

#' PER-DATE AUC-PR: average precision computed WITHIN each forecast date, then combined.
#'
#' Pooled AUC-PR ranks every date's rows in ONE list, so it rewards a model whose overall
#' probability LEVEL happens to be high on the dates that had invasions -- a property of the
#' level, not of how well the zones are ordered within the forecast an operator acts on. On the
#' 2026-10-09 frame the time-varying-beta arms rank zones IDENTICALLY to their fixed-beta base
#' within every date, yet scored 4-6 points higher pooled; and the static Flowminder-inflow rule
#' (one fixed list) was first of 30 pooled but twentieth per date.
#'
#' Definition, with E_f invasions, b_f base rate and AP_f the tie-aware AP among date f's at-risk
#' zones (dates with no invasion carry no ranking information and drop out):
#'     ap    = sum_f E_f * AP_f / sum_f E_f             (event-weighted mean AP)
#'     skill = sum_f E_f * AP_f / sum_f E_f * b_f       (AP relative to the base rate)
#' A skill of 1 is AP equal to the base rate. It is NOT what a random ordering scores: with few
#' invasions per date a random ordering's expected AP is well above the base rate. That chance
#' level is returned too, from the exact expectation of AP under a uniformly random ordering of
#' n zones with E invasions (Bestgen 2015):
#'     E[AP_f] = (1/n) * [ H_n + (E - 1)/(n - 1) * (n - H_n) ],   H_n the n-th harmonic number,
#'     chance  = sum_f E_f * E[AP_f] / sum_f E_f * b_f.
#' Use `chance` (not 1) as the reference in any "times better than chance" statement; the
#' selection compares methods on the same dates, so it is unaffected.
#'
#' @param p,y predicted probabilities and 0/1 outcomes.
#' @param date the forecast-date key per row (fold_id or cutoff).
#' @return c(ap, skill, n_dates, chance) -- n_dates counts the dates that carry an invasion.
.auc_pr_by_date <- function(p, y, date) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y) & !is.na(date)
  p <- p[ok]; y <- y[ok]; date <- date[ok]
  if (!length(y) || sum(y) == 0) return(c(ap = NA_real_, skill = NA_real_, n_dates = 0, chance = NA_real_))
  g <- split(seq_along(y), date)
  s <- vapply(g, function(i) {
    E <- sum(y[i]); if (E == 0) return(c(0, 0, 0, 0))
    n <- length(i); Hn <- sum(1 / seq_len(n))
    eap <- (Hn + if (n > 1L) (E - 1) / (n - 1) * (n - Hn) else 0) / n
    c(E * .auc_pr(p[i], y[i]), E, E * mean(y[i]), E * eap)
  }, numeric(4))
  if (sum(s[2, ]) == 0) return(c(ap = NA_real_, skill = NA_real_, n_dates = 0, chance = NA_real_))
  c(ap = sum(s[1, ]) / sum(s[2, ]), skill = sum(s[1, ]) / sum(s[3, ]), n_dates = sum(s[2, ] > 0),
    chance = sum(s[4, ]) / sum(s[3, ]))
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
#' way; only the interval columns move (auc_pr_lo/hi, log_score_lo/hi, auc_pr_skill_lo/hi,
#' auc_pr_skill_date_lo/hi, sel_auc_pr_skill_date_lo/hi) and their counts (n_boot_aucpr,
#' n_boot_sel).
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
#' the calibration figure (FigureS7_calibration_h1/_h2) is what actually shows the shape.
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
#' THE ANCHOR IS EVERY METHOD, NOT THE BEST-COVERED ONES. An earlier version took the
#' maximal per-model fold coverage, treated only the models AT that coverage as eligible,
#' intersected over those, and left every other model on its own native rows. That works
#' when one or two fits fail a fold; it inverts when a whole FAMILY covers fewer folds,
#' because the short family is then scored on its own smaller, differently-based support and
#' the best-covered pack gains an extra — typically event-poor, therefore easy — origin.
#' Anchoring on every method removes the possibility. Whether any family is currently short
#' is visible in the `n_folds` / `n_folds_scored` / `coverage` columns.
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
#' @param k_values top-K budgets for the ranking metrics. Defaults to TOPK_VALUES
#'   (00_config.R); the emitted prec_at_<K>/recall_at_<K>/hit_at_<K> columns are generated
#'   from this vector, so the constant and the published table cannot disagree.
#' @param n_boot cluster-bootstrap replicates (resample zones) for CIs.
#' @return tibble: one row per method × horizon with the metric set + CIs.
evaluate_invasion <- function(lfo_results,
                              k_values = get0("TOPK_VALUES", ifnotfound = c(5L, 15L, 20L)),
                              n_boot = 400L, seed = RANDOM_SEED) {
  stopifnot(all(c("method", "horizon", "fold_id", "health_zone",
                  "p_invasion", "is_new_invasion") %in% names(lfo_results)))
  # The bootstrap seeds per horizon (set_fixed_seed below); the caller's RNG stream is restored,
  # as invasion_selection_bootstrap() and the other seeded helpers here do.
  .had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  .old <- if (.had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({ if (.had) assign(".Random.seed", .old, envir = globalenv())
            else suppressWarnings(rm(".Random.seed", envir = globalenv())) }, add = TRUE)
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
    set_fixed_seed(.key_seed(seed, h))
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
    # 30_projection_config.R also reads as the cascade kernel. Under the old definition a
    # family that covers fewer folds BY DESIGN (the rolling-predictor floor can do this) was
    # flagged wholesale while the baselines stayed clean, so leaderboard, figures and cascade
    # could all fall back to a baseline over rows identical to the ones the Bayesian models
    # were scored on. Native reach stays published as n_folds/coverage.
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
    # PER-DATE AUC-PR, scored within each forecast date (see .auc_pr_by_date()). The pooled
    # `auc_pr` / `auc_pr_skill` columns are kept unchanged beside it for comparison.
    apd <- .auc_pr_by_date(p, y, d$fold_id)

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
    # SORTED, so every method enumerates the zones in the same order and the same seed draws the
    # same zones for every method (paired intervals). unique() followed row order, and the
    # ensembles' rows come back grouped differently, so they were resampled differently.
    zones <- sort(unique(d$health_zone))
    idx_by_zone <- split(seq_len(nrow(d)), factor(d$health_zone, levels = zones))
    # Post-burn-in rows (the selection axes, below), so the selection-axis interval comes from
    # the same zone resamples.
    # The SELECTION rounds: post-burn-in under the default basis, every round under
    # "per_date_all_rounds" (invasion_selection_drops_burnin()). The sel_* columns always mean
    # "the rounds the selection reads", so every consumer follows the switch.
    .burnin_ok <- "cutoff" %in% names(d) && exists("invasion_is_burnin", mode = "function") &&
                  invasion_selection_drops_burnin()
    .bi <- if (.burnin_ok) invasion_is_burnin(d$cutoff) else rep(FALSE, nrow(d))
    boot_aucpr <- numeric(0); boot_ls <- numeric(0); boot_apd <- numeric(0)
    boot_skill <- numeric(0); boot_sel_apd <- numeric(0)
    if (n_pos >= 2 && length(zones) > 5) {
      for (b in seq_len(n_boot)) {
        zs  <- sample.int(length(zones), length(zones), replace = TRUE)
        idx <- unlist(idx_by_zone[zs], use.names = FALSE)
        yb  <- y[idx]; pb <- p[idx]
        if (.calib) boot_ls <- c(boot_ls, .log_score(pb, yb))   # always defined when p is one
        if (sum(yb) >= 1) {
          boot_aucpr <- c(boot_aucpr, .auc_pr(pb, yb))
          # The SKILL of the replicate: its AP over ITS OWN base rate. (AP's interval divided
          # by the full-sample base rate is not an interval for the skill: the base rate moves
          # with the resample too.)
          boot_skill <- c(boot_skill, .auc_pr(pb, yb) / mean(yb))
          # Same zone resample, so the per-date and pooled intervals are paired.
          boot_apd <- c(boot_apd, .auc_pr_by_date(pb, yb, d$fold_id[idx])[["skill"]])
          ib <- idx[!.bi[idx]]
          if (length(ib) && sum(y[ib]) >= 1)
            boot_sel_apd <- c(boot_sel_apd, .auc_pr_by_date(p[ib], y[ib], d$fold_id[ib])[["skill"]])
        }
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
    aucpr_ci <- ci(boot_aucpr); ls_ci <- ci(boot_ls); apd_ci <- ci(boot_apd)
    skill_ci <- ci(boot_skill); sel_apd_ci <- ci(boot_sel_apd)

    # ---- SELECTION AXES: the selection's rounds -- post-burn-in under the default basis
    # (ASCERTAINMENT_BURNIN_WEEKS, 00_config.R), every round under "per_date_all_rounds" ----
    # The featured-model choice is made on the rounds after the discovery phase. The all-round
    # columns above stay as published; these `sel_*` columns are what best_invasion_model()
    # reads (via invasion_selection_view()). A frame without a `cutoff` column cannot be split,
    # so every round is used and `sel_burnin_applied` says so.
    ds  <- d[!.bi, , drop = FALSE]                     # .bi: computed above the bootstrap
    ys  <- as.integer(ds$is_new_invasion); ps <- ds$p_invasion
    sel_apd <- .auc_pr_by_date(ps, ys, ds$fold_id)
    sel_rk  <- if (nrow(ds) && sum(ys) > 0) .ranking_metrics(ds, k_values) else NULL
    sel_ls    <- if (.calib && nrow(ds)) .log_score(ps, ys) else NA_real_
    sel_ls_rc <- if (has_recal && nrow(ds)) .log_score(ds$p_recal, ys) else NA_real_

    # n_invasions counts POSITIVE ROWS; n_invasion_events counts DISTINCT invaded zones.
    # They are equal at h=1 (disjoint windows) and differ at h>=2, where the cumulative
    # windows overlap and one invasion is scored from two origins (n_invasion_events against
    # n_invasions at h=2 shows how many). Publishing only the row count invites reading it as an effective
    # sample size; publishing both makes the dependence visible in the table itself.
    n_events_distinct <- length(unique(d$health_zone[d$is_new_invasion == 1L]))
    .row <- tibble::tibble(
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
      # 90% zone-bootstrap interval OF THE SKILL (each replicate's AP over its own base rate).
      auc_pr_skill_lo = skill_ci[1], auc_pr_skill_hi = skill_ci[2],
      # Per-date AUC-PR (the selection axis), all rounds; interval from the SAME zone resamples.
      auc_pr_date = apd[["ap"]], auc_pr_skill_date = apd[["skill"]],
      # The per-date skill a uniformly random ordering would score on these dates (not 1).
      auc_pr_skill_date_chance = apd[["chance"]],
      auc_pr_skill_date_lo = apd_ci[1], auc_pr_skill_date_hi = apd_ci[2],
      n_dates_with_invasion = as.integer(apd[["n_dates"]]),
      auc_roc = aucroc,
      log_score = lsc, log_score_lo = ls_ci[1], log_score_hi = ls_ci[2],
      # Retained bootstrap replicates. n_boot_aucpr < n_boot means AUC-PR replicates were
      # dropped for containing no event, so its interval is conditional on >= 1 event; the
      # log-score interval uses every replicate and is unconditional.
      n_boot_aucpr = n_boot_aucpr, n_boot_log_score = n_boot_ls,
      # Replicates behind sel_auc_pr_skill_date_lo/hi: only those with at least one
      # post-burn-in event are kept.
      n_boot_sel = length(boot_sel_apd),
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
      # The three selection axes on POST-BURN-IN rounds (see the block above).
      sel_burnin_applied = .burnin_ok,
      sel_n_folds = length(unique(ds$fold_id)), sel_n_invasions = sum(ys),
      # Of those rounds, how many carry a FITTED prequential delta. The others had no settled,
      # post-burn-in outcome to train on by their issue date, so their recalibrated
      # probability IS the raw one (delta = 1): the recalibrated axis is raw on those rounds.
      sel_n_folds_delta_fitted = if ("delta_estimable" %in% names(ds))
        length(unique(ds$fold_id[ds$delta_estimable %in% TRUE])) else NA_integer_,
      sel_base_rate = if (length(ys)) mean(ys) else NA_real_,
      sel_auc_pr_skill_date = sel_apd[["skill"]],
      # its 90% zone-bootstrap interval (the same resamples as the all-round columns)
      sel_auc_pr_skill_date_lo = sel_apd_ci[1], sel_auc_pr_skill_date_hi = sel_apd_ci[2],
      sel_auc_pr_skill_date_chance = sel_apd[["chance"]],
      # Its numerator and denominator: the event-weighted per-date AP and the event-weighted
      # base rate (skill = ap / prev), so a consumer can rescale it to another set's base rate.
      sel_auc_pr_date = sel_apd[["ap"]],
      sel_prev_date = if (is.finite(sel_apd[["skill"]]) && sel_apd[["skill"]] > 0)
        sel_apd[["ap"]] / sel_apd[["skill"]] else NA_real_,
      sel_mean_rank_of_truth = sel_rk$mean_rank %||% NA_real_,
      sel_log_score = sel_ls, sel_log_score_recal = sel_ls_rc
    )
    # TOP-K COLUMNS GENERATED FROM `k_values`, not hard-coded. Three K's used to be written
    # out by hand, so changing TOPK_VALUES left the table carrying columns for the OLD
    # budgets (filled with NA) and silently missing the new ones — and every consumer that
    # asked for `recall_at_10` by name fell back to `recall_at_5` without saying so.
    for (K in k_values) {
      .row[[paste0("prec_at_",   K)]] <- rk[[paste0("prec@",   K)]] %||% NA_real_
      .row[[paste0("recall_at_", K)]] <- rk[[paste0("recall@", K)]] %||% NA_real_
      .row[[paste0("hit_at_",    K)]] <- rk[[paste0("hit@",    K)]] %||% NA_real_
    }
    .row
  })
  out <- dplyr::bind_rows(Filter(Negate(is.null), rows))
  # Primary ranking: full-CV (comparable common-support) models first, then AUC-PR skill
  # (discrimination over the shared base rate), tie-broken by lower mean rank-of-truth.
  # partial_cv models (scored on fewer cells than the shared support -- only reachable via
  # the fallback) sort last so they cannot outrank the comparable pack on a thin subset.
  out %>% dplyr::arrange(horizon, partial_cv, dplyr::desc(auc_pr_skill), mean_rank_of_truth)
}

MRT_GATE_MULT <- 1.5   # spiky-model gate: drop methods with worst-horizon rank-of-truth > this x field median
# SELECTION BASIS (2026-10-14; all-rounds option 2026-10-19). Which version of the three
# composite axes the featured-model choice reads (environment variable INVASION_SELECTION_BASIS):
#   "per_date_post_burnin" (DEFAULT) — AUC-PR scored WITHIN each forecast date
#       (sel_auc_pr_skill_date), mean rank-of-truth and log score, all on the rounds after the
#       ascertainment burn-in (ASCERTAINMENT_BURNIN_WEEKS, 00_config.R). A DELIBERATE CHOICE:
#       in the discovery phase where first cases were confirmed partly reflected where
#       surveillance looked first, so those rounds are kept out of the choice of model.
#       PERFORMANCE IS STILL REPORTED OVER ALL ROUNDS (the figures and tables quote the
#       all-round per-date columns; the selection rounds' values are shown beside them).
#   "per_date_all_rounds" — the same per-date axes over EVERY round (no burn-in).
#   "pooled_all_dates" — the PREVIOUS rule: pooled AUC-PR, every round. Retained so the pick
#       under the old rule is reported beside the new one (write_model_selection.R).
# An unrecognised value is an error, not a fallback: it would select under a rule nobody chose.
INVASION_SELECTION_BASES <- c("per_date_post_burnin", "per_date_all_rounds", "pooled_all_dates")
INVASION_SELECTION_BASIS <- local({
  b <- tolower(trimws(Sys.getenv("INVASION_SELECTION_BASIS", unset = "per_date_post_burnin")))
  if (!b %in% INVASION_SELECTION_BASES)
    stop(sprintf("[best_model] INVASION_SELECTION_BASIS must be one of %s; got '%s'.",
                 paste(INVASION_SELECTION_BASES, collapse = ", "), b), call. = FALSE)
  b
})
#' TRUE when the selection's rounds exclude the ascertainment burn-in (the default basis).
invasion_selection_drops_burnin <- function(basis = get0("INVASION_SELECTION_BASIS",
                                                         ifnotfound = "per_date_post_burnin"))
  identical(basis, "per_date_post_burnin")
#' TRUE when the selection reads the per-date axes (the sel_* columns), on whichever rounds.
invasion_selection_per_date <- function(basis = get0("INVASION_SELECTION_BASIS",
                                                     ifnotfound = "per_date_post_burnin"))
  basis %in% c("per_date_post_burnin", "per_date_all_rounds")
#' The selection's AUC-PR axis in words, for tables and figures that report ON it (the held-out
#' optimism check, FigR4): which estimator, on which rounds. Everything else reports all rounds.
invasion_selection_axis_label <- function(basis = get0("INVASION_SELECTION_BASIS",
                                                       ifnotfound = "per_date_post_burnin"))
  switch(basis,
         per_date_post_burnin = "per-date AUC-PR skill, post-burn-in rounds",
         per_date_all_rounds  = "per-date AUC-PR skill, all rounds",
         pooled_all_dates     = "pooled AUC-PR skill, all rounds",
         stop("[best_model] unknown selection basis '", basis, "'.", call. = FALSE))

#' The evaluation table as the SELECTION sees it.
#'
#' Maps the standard axis columns -- auc_pr_skill, mean_rank_of_truth, log_score,
#' log_score_recal -- onto the chosen basis, so best_invasion_model(),
#' invasion_selection_table() and write_model_selection.R all rank on the SAME numbers through
#' ONE function. The all-round originals are preserved in `orig_*` columns (written once, never
#' overwritten), which makes the mapping idempotent and reversible: viewing an already-viewed
#' table under either basis gives the same result as viewing the original.
#'
#' A table without the `sel_*` columns (written before 2026-10-14) cannot be viewed under the
#' new basis; it is returned under the old one with a warning rather than silently.
#' @return the table with the axis columns replaced and attr "selection_basis" set.
invasion_selection_view <- function(eval_tbl, basis = INVASION_SELECTION_BASIS) {
  e <- eval_tbl
  if (is.null(e) || !nrow(e)) return(e)
  axes <- c(auc_pr_skill = "sel_auc_pr_skill_date", mean_rank_of_truth = "sel_mean_rank_of_truth",
            log_score = "sel_log_score", log_score_recal = "sel_log_score_recal")
  # Preserve each all-round original ONCE (and only if it exists: a column that was never there
  # must stay absent, so the downstream "recalibrated column absent" fallback still fires).
  for (a in names(axes)) {
    o <- paste0("orig_", a)
    if (!o %in% names(e) && a %in% names(e)) e[[o]] <- e[[a]]
  }
  if (invasion_selection_per_date(basis) && !all(axes %in% names(e))) {
    warning(sprintf(paste0("[best_model] the evaluation table has no %s column(s) (written before ",
                           "2026-10-14?); selecting on the PREVIOUS basis (pooled AUC-PR, all ",
                           "rounds)."), paste(setdiff(axes, names(e)), collapse = ", ")),
            call. = FALSE)
    basis <- "pooled_all_dates"
  }
  for (a in names(axes)) {
    src <- if (invasion_selection_per_date(basis)) axes[[a]] else paste0("orig_", a)
    e[[a]] <- if (src %in% names(e)) e[[src]] else NULL
  }
  attr(e, "selection_basis") <- basis
  e
}
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

#' Pick the best model POOLED across forecast horizons: composite-led, spike-guarded.
#'
#' THE AXES ARE READ THROUGH invasion_selection_view() (2026-10-14). By default the three axes
#' are computed on the rounds AFTER the ascertainment burn-in, and the AUC-PR axis is scored
#' WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled into one list.
#' `basis = "pooled_all_dates"` reproduces the previous rule for comparison.
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

#' @param eval_tbl an evaluate_invasion() table.
#' @param horizons optional integer vector; if given, only these horizons are pooled
#'   (default: every horizon present in eval_tbl).
#' @param restrict optional regex; if given, only methods matching it are ranked
#'   (e.g. "^Bayes" to pick the best Bayesian model). STRICT: if nothing matches, NA is
#'   returned rather than falling back to the unrestricted field.
#' @param rule "composite" (default) or "aucpr" to recover the previous
#'   discrimination-only behaviour. Override globally with INVASION_SELECTION_RULE.
#' @param select_on_recal put the composite's log-score axis on the recalibrated score
#'   (INVASION_SELECT_ON_RECAL, 00_config.R).
#' @param require_calibrated exclude rank-only comparators (prob_calibrated = FALSE) from
#'   candidacy. They remain scored and ranked in the table; they simply cannot be featured,
#'   because every headline product reads p_invasion AS a probability.
#' @param exclude regex of method labels that are SCORED but are not CANDIDATES here. The
#'   default, INVASION_SELECTION_EXCLUDE (00_config.R), is STAGE 1 of the featured-model
#'   selection: it removes the sensitivity arms (SENSITIVITY_ARM_SUFFIXES: the generation-time
#'   arms, never candidates, which would otherwise make the GT a selection axis) AND the
#'   time-varying forms (TV_FORM_SUFFIXES), which are candidates only in STAGE 2, against the
#'   fixed-beta model they generalise (select_featured_two_stage()). Applied AFTER `restrict`.
#'   Set NULL to make every scored model a candidate.
best_invasion_model <- function(eval_tbl, horizons = NULL, restrict = NULL,
                                rule = INVASION_SELECTION_RULE,
                                select_on_recal = INVASION_SELECT_ON_RECAL,
                                exclude = get0("INVASION_SELECTION_EXCLUDE",
                                               ifnotfound = "-(gtshort|gtlong|tv[a-z0-9]+)$"),
                                require_calibrated = TRUE,
                                basis = INVASION_SELECTION_BASIS) {
  # Every axis below is read through the selection view: per-date AUC-PR and post-burn-in
  # rounds by default (see invasion_selection_view()).
  e <- invasion_selection_view(eval_tbl, basis)
  if (!is.null(exclude) && length(exclude) == 1L && !is.na(exclude) && nzchar(exclude)) {
    .drop <- grepl(exclude, e$method)
    if (any(.drop)) {
      message(sprintf(paste0("[best_model] %d method(s) scored but not candidates here (sensitivity ",
                             "arms, never candidates; time-varying forms, stage-2 candidates only): %s"),
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

#' How firm is the featured-model choice? A zone-cluster bootstrap of the SELECTION.
#'
#' Each replicate resamples zones (the correlated unit) with replacement, recomputes the three
#' selection axes ON THE SELECTION BASIS (`basis`, INVASION_SELECTION_BASIS): per-date AUC-PR
#' skill on the post-burn-in rounds, or pooled AUC-PR skill on all rounds -- with the mean rank
#' of the invaded zones and the log score (recalibrated when `select_on_recal`) on the same
#' rounds -- on the shared (fold x zone) support, and runs the TWO-STAGE rule
#' (select_featured_two_stage()) through best_invasion_model() itself, so the gates, the
#' composite and the tie-breaks are the selector's own, not a copy. Stage 2 can only consider
#' forms that were cross-validated, i.e. those of THIS run's stage-1 pick: a replicate whose
#' stage 1 picks another specification features that specification (fixed beta). The prequential
#' recalibration factors are held at their fitted values (they are not refitted per replicate).
#'
#' @param lfo_results the LFO frame AFTER attach_invasion_recalibration().
#' @param B replicates.
#' @return tibble(method, kernel, share, stage1_share) -- the share of replicates each model is
#'   FEATURED, and the share it wins stage 1 -- with attr "kernel_share" (featured share by
#'   mobility kernel) and attr "B".
invasion_selection_bootstrap <- function(lfo_results, B = get0("INVASION_SELECTION_BOOT", ifnotfound = 200L),
                                         restrict = "^Bayes",
                                         select_on_recal = INVASION_SELECT_ON_RECAL,
                                         seed = RANDOM_SEED,
                                         basis = INVASION_SELECTION_BASIS) {
  per_date <- invasion_selection_per_date(basis)
  d <- lfo_results
  if ("was_active_before" %in% names(d)) d <- d[!(as.logical(d$was_active_before) %in% TRUE), , drop = FALSE]
  d <- d[is.finite(d$p_invasion), , drop = FALSE]
  cells <- lapply(sort(unique(d$horizon)), function(h) invasion_common_cells(lfo_results, h))
  names(cells) <- as.character(sort(unique(d$horizon)))
  d$.cell <- paste(d$fold_id, d$health_zone, sep = "\r")
  keep <- rep(TRUE, nrow(d))                  # a horizon without a shared support keeps its rows
  for (h in names(cells)) if (length(cells[[h]])) {
    i <- d$horizon == as.integer(h)
    keep[i] <- d$.cell[i] %in% cells[[h]]
  }
  d <- d[keep, , drop = FALSE]
  # Sensitivity arms never compete; the time-varying forms stay in (stage 2).
  arms <- get0("SENSITIVITY_ARM_REGEX", ifnotfound = NULL)
  d <- d[grepl(restrict, d$method) & !grepl("-ens-", d$method), , drop = FALSE]
  if (!is.null(arms) && nzchar(arms)) d <- d[!grepl(arms, d$method), , drop = FALSE]
  if ("prob_calibrated" %in% names(d)) d <- d[!(d$prob_calibrated %in% FALSE), , drop = FALSE]
  if (per_date && invasion_selection_drops_burnin(basis) && "cutoff" %in% names(d) &&
      exists("invasion_is_burnin", mode = "function"))
    d <- d[!invasion_is_burnin(d$cutoff), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  use_recal <- isTRUE(select_on_recal) && "p_recal" %in% names(d) && all(is.finite(d$p_recal))
  axes <- function(dd) {
    dd %>% dplyr::group_by(method, horizon) %>% dplyr::group_modify(function(g, k) {
      y <- as.integer(g$is_new_invasion)
      rot <- vapply(split(seq_len(nrow(g)), g$fold_id), function(i) {
        if (!sum(y[i])) return(NA_real_)
        mean(rank(-g$p_invasion[i], ties.method = "average")[y[i] == 1L]) }, numeric(1))
      pp <- pmin(pmax(g$p_invasion, 1e-6), 1 - 1e-6)
      pr <- if (use_recal) pmin(pmax(g$p_recal, 1e-6), 1 - 1e-6) else pp
      ap  <- if (per_date) .auc_pr_by_date(g$p_invasion, y, g$fold_id)[["skill"]]
             else .auc_pr(g$p_invasion, y) / max(mean(y), 1e-9)
      out <- tibble::tibble(auc_pr_skill = ap,
                            mean_rank_of_truth = mean(rot, na.rm = TRUE),
                            log_score = -mean(y * log(pp) + (1 - y) * log(1 - pp)),
                            log_score_recal = -mean(y * log(pr) + (1 - y) * log(1 - pr)),
                            partial_cv = FALSE, prob_calibrated = TRUE)
      # The per-date basis reads the sel_* columns (invasion_selection_view()).
      if (per_date) out <- dplyr::rename(out, sel_auc_pr_skill_date = auc_pr_skill,
                                         sel_mean_rank_of_truth = mean_rank_of_truth,
                                         sel_log_score = log_score,
                                         sel_log_score_recal = log_score_recal)
      out
    }) %>% dplyr::ungroup()
  }
  pick <- function(dd) tryCatch(suppressMessages(suppressWarnings({
    tw <- select_featured_two_stage(axes(dd), restrict = restrict, select_on_recal = use_recal,
                                    basis = basis)
    c(tw$stage1, tw$featured) })), error = function(e) c(NA_character_, NA_character_))
  # Reproducible; the caller's RNG stream is restored.
  .had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  .old <- if (.had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({ if (.had) assign(".Random.seed", .old, envir = globalenv())
            else suppressWarnings(rm(".Random.seed", envir = globalenv())) }, add = TRUE)
  set_fixed_seed(as.integer(seed) + 3301L)
  # Sorted (C locale), so the resamples -- and the published shares -- do not depend on the
  # frame's row order (evaluate_invasion() sorts for the same reason).
  zones <- sort(unique(d$health_zone), method = "radix"); idx <- split(seq_len(nrow(d)), d$health_zone)
  wins2 <- vapply(seq_len(B), function(b) {
    zs <- sample(zones, length(zones), replace = TRUE)
    rows <- unlist(idx[zs], use.names = FALSE)
    dd <- d[rows, , drop = FALSE]
    # a zone drawn twice must count as two zones inside a fold's ranking, not one
    dd$health_zone <- paste0(dd$health_zone, "#", rep(seq_along(zs), vapply(idx[zs], length, 1L)))
    pick(dd)
  }, character(2))
  ok <- !is.na(wins2[2, ])
  wins <- wins2[2, ok]; s1 <- wins2[1, ok]
  kern <- vapply(wins, function(m) as.character(mobility_kernel_from_method(m) %||% NA_character_), "")
  meths <- sort(unique(c(wins, s1)), method = "radix")
  out <- tibble::tibble(method = meths,
                        share = as.numeric(vapply(meths, function(m) sum(wins == m), 0)) / length(wins),
                        stage1_share = as.numeric(vapply(meths, function(m) sum(s1 == m), 0)) / length(wins)) %>%
    dplyr::mutate(kernel = vapply(method, function(m) as.character(mobility_kernel_from_method(m) %||% NA_character_), ""),
                  .after = method) %>%
    dplyr::arrange(dplyr::desc(share), dplyr::desc(stage1_share))
  attr(out, "kernel_share") <- sort(table(kern) / length(kern), decreasing = TRUE)
  attr(out, "B") <- length(wins)
  out
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
                                                   ifnotfound = "-(gtshort|gtlong|tv[a-z0-9]+)$"),
                                    basis = INVASION_SELECTION_BASIS) {
  e <- invasion_selection_view(eval_tbl, basis)   # same view as best_invasion_model()
  # Same non-candidate gate as best_invasion_model(), for the same reason: an audit table
  # that ranks a model the selection could never pick describes a different competition
  # from the one that was held.
  if (!is.null(exclude) && length(exclude) == 1L && !is.na(exclude) && nzchar(exclude))
    e <- e[!grepl(exclude, e$method), , drop = FALSE]
  # Rank-only comparators are scored but not eligible, exactly as in best_invasion_model()
  # (require_calibrated); without this they entered the spiky-gate median and the within-horizon
  # ranks, so the audit composite could disagree with the pick.
  if ("prob_calibrated" %in% names(e)) {
    .rank_only <- !(e$prob_calibrated %in% TRUE)
    if (any(.rank_only) && any(!.rank_only)) e <- e[!.rank_only, , drop = FALSE]
  }
  if (!is.null(restrict)) e <- e %>% dplyr::filter(grepl(restrict, method))
  if (!is.null(horizons)) e <- e %>% dplyr::filter(horizon %in% !!horizons)
  if (nrow(e) == 0) return(NULL)
  # Horizon coverage FIRST, then the partial-CV gate: the order best_invasion_model() applies
  # them in (its `agg` filters n_h == max before any_partial). In the other order a method
  # partial at one horizon could set max(n_h) and the two would rank different fields.
  nh <- e %>% dplyr::group_by(method) %>%
    dplyr::summarise(n_h = dplyr::n_distinct(horizon), .groups = "drop")
  e <- e %>% dplyr::filter(method %in% nh$method[nh$n_h == max(nh$n_h)])
  if ("partial_cv" %in% names(e)) {
    full <- e %>% dplyr::group_by(method) %>%
      dplyr::summarise(any_partial = any(partial_cv %in% TRUE), .groups = "drop")
    if (any(!full$any_partial))
      e <- e %>% dplyr::filter(method %in% full$method[!full$any_partial])
  }
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

#' The FEATURED model by the two-stage rule (00_config.R, TIME-VARYING FORMS).
#'
#' STAGE 1 picks the specification (kernel, covariates) among the fixed-beta candidates:
#' best_invasion_model() with INVASION_SELECTION_EXCLUDE, which removes the sensitivity arms AND
#' the time-varying forms. STAGE 2 puts that pick and its cross-validated SELECTABLE forms
#' ("<stage1>-tvrw1"; TV_FORM_SELECTABLE_SUFFIXES) through the SAME composite, gates and
#' tie-breaks, on their own; the winner is featured. The secondary forms (the frailty ones,
#' TV_FORM_SECONDARY_SUFFIXES) are scored but never enter the pool. Without a cross-validated
#' selectable form the stage-1 pick is featured.
#'
#' @param eval_tbl evaluate_invasion() output.
#' @param restrict candidate family regex (ensembles are always excluded).
#' @return list(stage1, forms (the cross-validated selectable forms of stage1), secondary (its
#'   cross-validated secondary forms, never candidates), pool, featured, table (the stage-2
#'   invasion_selection_table(), or NULL)).
select_featured_two_stage <- function(eval_tbl, restrict = "^Bayes",
                                      select_on_recal = INVASION_SELECT_ON_RECAL,
                                      basis = INVASION_SELECTION_BASIS) {
  e <- eval_tbl[grepl(restrict, eval_tbl$method) & !grepl("-ens-", eval_tbl$method), , drop = FALSE]
  s1 <- best_invasion_model(e, select_on_recal = select_on_recal, basis = basis,
                            exclude = INVASION_SELECTION_EXCLUDE)
  out <- list(stage1 = s1, forms = character(0), secondary = character(0), pool = s1,
              featured = s1, table = NULL)
  if (length(s1) != 1L || is.na(s1)) return(out)
  .cv_forms <- function(sfx) intersect(paste0(s1, "-", sfx), unique(e$method))
  out$forms <- .cv_forms(get0("TV_FORM_SELECTABLE_SUFFIXES", ifnotfound = character(0)))
  out$secondary <- .cv_forms(get0("TV_FORM_SECONDARY_SUFFIXES", ifnotfound = character(0)))
  if (!length(out$forms)) return(out)
  out$pool <- c(s1, out$forms)
  ep <- e[e$method %in% out$pool, , drop = FALSE]
  out$featured <- best_invasion_model(ep, select_on_recal = select_on_recal, basis = basis,
                                      exclude = SENSITIVITY_ARM_REGEX)
  out$table <- invasion_selection_table(ep, select_on_recal = select_on_recal,
                                        exclude = SENSITIVITY_ARM_REGEX, basis = basis)
  out
}

#' The time-varying forms against fixed beta, on the stage-1 specification: calibration and
#' discrimination per form, and the paired log-score gain of each form over fixed beta with
#' its uncertainty. (This was 51_tv_comparison.R's pre-registered GUARD, which decided between
#' the forms; since 2026-10-15 the selectable forms compete for the featured model on the
#' selection composite like every other candidate, and the guard is REPORTED alongside --
#' `guard_passed` records whether a form would have passed it. Since 2026-10-19 the frailty
#' forms are secondary: compared here, never candidates; see TV_FORMS, 00_config.R.)
#'
#' Every metric and the guard are computed on the post-burn-in rounds (ASCERTAINMENT_BURNIN_WEEKS),
#' both horizons, on the at-risk cells every arm covers, from the back-test the selection used
#' (`by_date` also carries the burn-in rounds, flagged `burnin`, for the figures):
#'   metrics  per form x horizon: log score (raw and recalibrated) x1000; expected/observed raw,
#'            recalibrated, and recalibrated in the early / late half of the rounds; per-date
#'            AUC-PR skill; mean rank of the invaded zones.
#'   guard    per form and probability column (p_recal, p_invasion): the horizon-mean paired
#'            log-score gain over fixed beta (> 0 = the form is better) x1000, with a zone
#'            bootstrap (`B` replicates) and an origin-level Diebold-Mariano bound with the
#'            Harvey-Leybourne-Newbold correction (t on n_origins - 1 df, autocovariance lag
#'            = longest horizon - 1, the dependence the nested outcome windows create; lag 0 if
#'            that long-run variance is not positive, as forecast::dm.test). Bounds at the
#'            Bonferroni level 5%/k (k = number of forms) and at 5% (the 90% interval).
#'            guard_passed = both Bonferroni bounds > 0 on the recalibrated forecasts.
#'   by_date  per form x horizon x round: observed, expected raw and recalibrated, delta_preq.
#' The zone bootstrap is reproducible from `seed` and leaves the caller's RNG stream untouched.
#' @param lfo_results the back-test AFTER attach_invasion_recalibration(), holding the stage-1
#'   model and its forms.
#' @param stage1 the stage-1 label; its forms are "<stage1>-<TV_FORM_SUFFIXES>".
#' @return list(metrics, guard, by_date, n_rounds, k, bonferroni_level, B, labels), or NULL when fewer
#'   than one form or no post-burn-in round is available.
invasion_tv_form_comparison <- function(lfo_results, stage1, B = 4000L,
                                        seed = get0("RANDOM_SEED", ifnotfound = 20260704L)) {
  sfx <- get0("TV_FORM_SUFFIXES", ifnotfound = character(0))
  forms_tbl <- get0("TV_FORMS", ifnotfound = list())
  labs <- c(fixed = stage1)
  for (x in sfx) labs[[forms_tbl[[x]]$form]] <- paste0(stage1, "-", x)
  labs <- labs[labs %in% unique(lfo_results$method)]
  if (!"fixed" %in% names(labs) || length(labs) < 2L) return(NULL)
  forms <- setdiff(names(labs), "fixed")
  if (!all(c("p_recal", "delta_preq") %in% names(lfo_results)))
    stop("[tv-forms] the comparison needs the recalibrated back-test (attach_invasion_recalibration()).",
         call. = FALSE)
  lsv <- function(p, y) { p <- pmin(pmax(p, 1e-6), 1 - 1e-6); -(y * log(p) + (1 - y) * log(1 - p)) }
  .active <- if ("was_active_before" %in% names(lfo_results))
    as.logical(lfo_results$was_active_before) %in% TRUE else rep(FALSE, nrow(lfo_results))
  d0 <- lfo_results[lfo_results$method %in% labs & !.active &
                      is.finite(lfo_results$p_invasion) & is.finite(lfo_results$p_recal), , drop = FALSE]
  d0$method <- names(labs)[match(d0$method, labs)]
  d0$burnin <- invasion_is_burnin(d0$cutoff)
  key <- d0 %>% dplyr::count(cutoff, health_zone, horizon) %>% dplyr::filter(n == length(labs))
  d0 <- d0 %>% dplyr::semi_join(key, by = c("cutoff", "health_zone", "horizon"))
  P <- d0[!d0$burnin, , drop = FALSE]
  if (!nrow(P)) return(NULL)
  mid <- stats::median(unique(as.Date(P$cutoff)))
  eo <- function(p, y) if (sum(y) > 0) sum(p) / sum(y) else NA_real_

  by_date <- d0 %>% dplyr::group_by(form = method, horizon, cutoff, burnin) %>%
    dplyr::summarise(observed = sum(is_new_invasion), expected_raw = sum(p_invasion),
                     expected_recal = sum(p_recal), delta = dplyr::first(delta_preq),
                     # FALSE: no factor was fitted for this round yet, so its "recalibrated"
                     # probabilities are the raw ones (label them "(uncorrected)").
                     delta_estimable = if ("delta_estimable" %in% names(d0))
                       any(delta_estimable %in% TRUE) else NA,
                     n_atrisk = dplyr::n(), .groups = "drop")
  metrics <- P %>% dplyr::group_by(form = method, horizon) %>% dplyr::group_modify(function(g, k) {
    ap <- .auc_pr_by_date(g$p_invasion, g$is_new_invasion, g$fold_id)
    rot <- vapply(split(seq_len(nrow(g)), g$fold_id), function(i) {
      y <- g$is_new_invasion[i]
      if (!sum(y)) return(NA_real_)
      mean(rank(-g$p_invasion[i], ties.method = "average")[y == 1]) }, numeric(1))
    early <- as.Date(g$cutoff) <= mid
    tibble::tibble(
      n_rounds = dplyr::n_distinct(g$cutoff), n_invasions = sum(g$is_new_invasion),
      log_score_raw_x1000 = 1000 * mean(lsv(g$p_invasion, g$is_new_invasion)),
      log_score_recal_x1000 = 1000 * mean(lsv(g$p_recal, g$is_new_invasion)),
      eo_raw = eo(g$p_invasion, g$is_new_invasion), eo_recal = eo(g$p_recal, g$is_new_invasion),
      eo_recal_early = eo(g$p_recal[early], g$is_new_invasion[early]),
      eo_recal_late  = eo(g$p_recal[!early], g$is_new_invasion[!early]),
      auc_pr_skill_per_date = ap[["skill"]], mean_rank_of_truth = mean(rot, na.rm = TRUE))
  }) %>% dplyr::ungroup()

  k <- length(forms); qlo <- 0.05 / k
  .dm_bound <- function(w) {
    hs <- sort(unique(w$horizon)); cs <- sort(unique(w$cutoff)); n <- length(cs)
    # One origin carries no between-origin variance: no bound (NA), rather than NaN on 0 df.
    if (n < 2L) return(list(d = rep(mean(tapply(w$g, w$horizon, mean)), n), se = NA_real_,
                            lag = NA_integer_, lower = function(a) NA_real_))
    Nh <- vapply(hs, function(h) sum(w$horizon == h), numeric(1))
    # One contribution per origin, d_c = (n / n_horizons) * sum_h S_ch / N_h, so mean(d) IS
    # the gain statistic (checked below).
    dd <- vapply(cs, function(c) sum(vapply(seq_along(hs), function(j)
            sum(w$g[w$cutoff == c & w$horizon == hs[j]]) / Nh[j], numeric(1))), numeric(1)) *
          n / length(hs)
    m <- mean(dd)
    lrv <- function(L) { gam <- vapply(0:L, function(l) sum((dd[(1 + l):n] - m) * (dd[1:(n - l)] - m)) / n,
                                       numeric(1)); gam[1] + 2 * sum(gam[-1]) }
    L <- min(max(hs) - 1L, n - 1L); V <- lrv(L)
    if (!(V > 0)) { L <- 0L; V <- lrv(0L) }
    hh <- L + 1L; cc <- sqrt((n + 1 - 2 * hh + hh * (hh - 1) / n) / n)
    se <- sqrt(V / n) / cc
    list(d = dd, se = se, lag = L, lower = function(a) m - stats::qt(1 - a, df = n - 1L) * se)
  }
  .had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  .old <- if (.had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({ if (.had) assign(".Random.seed", .old, envir = globalenv())
            else suppressWarnings(rm(".Random.seed", envir = globalenv())) }, add = TRUE)
  guard_one <- function(form, col) {
    w <- P %>% dplyr::filter(method %in% c(form, "fixed")) %>%
      dplyr::select(method, cutoff, health_zone, horizon, p = dplyr::all_of(col), y = is_new_invasion) %>%
      tidyr::pivot_wider(names_from = method, values_from = p)
    w$g <- lsv(w$fixed, w$y) - lsv(w[[form]], w$y)
    stat <- function(rows) { s <- w[rows, ]; mean(tapply(s$g, s$horizon, mean)) }
    pt <- stat(seq_len(nrow(w)))
    set_fixed_seed(as.integer(seed) + 5101L)
    byz <- split(seq_len(nrow(w)), w$health_zone)
    bz <- replicate(B, stat(unlist(byz[sample(length(byz), replace = TRUE)], use.names = FALSE)))
    dm <- .dm_bound(w)
    if (!isTRUE(all.equal(mean(dm$d), pt)))
      stop("[tv-forms] per-origin contributions do not reproduce the gain statistic.", call. = FALSE)
    tibble::tibble(form = form, forecasts = col, gain_x1000 = 1000 * pt,
                   zone_lo_bonferroni = 1000 * unname(stats::quantile(bz, qlo)),
                   origin_lo_bonferroni = 1000 * dm$lower(qlo),
                   zone_lo90 = 1000 * unname(stats::quantile(bz, 0.05)),
                   zone_hi90 = 1000 * unname(stats::quantile(bz, 0.95)),
                   origin_lo90 = 1000 * dm$lower(0.05),
                   origin_hi90 = 1000 * (2 * pt - dm$lower(0.05)),
                   origin_se_x1000 = 1000 * dm$se, n_origins = length(dm$d), dm_lag = dm$lag)
  }
  guard <- dplyr::bind_rows(lapply(forms, guard_one, col = "p_recal"),
                            lapply(forms, guard_one, col = "p_invasion")) %>%
    dplyr::mutate(guard_passed = forecasts == "p_recal" &
                    (pmin(zone_lo_bonferroni, origin_lo_bonferroni) > 0) %in% TRUE)
  list(metrics = metrics, guard = guard, by_date = by_date,
       n_rounds = dplyr::n_distinct(P$cutoff), k = k, bonferroni_level = qlo, B = as.integer(B),
       labels = labs)
}

#' The WHOLE FIELD's inner-vs-outer skill at one held-out split.
#'
#' WHY THIS EXISTS. `optimism_gap` conflates two effects that have nothing to do with each
#' other. (a) THE WINDOW: the outer origins are late ones, and on this outbreak they are
#' sparser and lower-base-rate than the inner ones (h=1 base rate falls 0.0080 at the first
#' round to 0.0043 at the last), so EVERY model scores worse there — the field's median drop
#' over the 18 eligible candidates is 20.2 skill points at h=1 on the 2026-09-07 frame, and
#' the drop is NEGATIVE for two models, i.e. they improve. (b) SELECTION: the model picked on
#' the inner folds may lose MORE than the others because it was picked there. Only (b) is
#' optimism. Scoring every candidate on the SAME two sets separates them: the field median is
#' the window effect, and the selected model's EXCESS over it is what selection cost.
#'
#' Measured on the shipped folds at k=3, the excess is +2.0 skill points at h=1 and +6.3 at
#' h=2 against a published conflated gap of 27.8 — i.e. the headline number was almost
#' entirely the window. The rank correlation is the third piece: at rho = 0.81 (h=1) the
#' field's ORDERING survives out of sample even though every level falls, which is the
#' property selection actually relies on.
#'
#' USE THE COMMON-BASE COLUMNS FOR THE DECOMPOSITION. `drop` is a difference of two skills
#' scored against their OWN base rates, so the base-rate shift enters it; at h=2 that pushes
#' the field median to -1.1 and makes the whole table read as if nothing degrades.
#' `drop_common` puts both sets over the inner base rate, so the difference is a difference
#' of numerators. Both are returned; the figure and the summary use the common-base one.
#'
#' @param inner_eval,outer_eval evaluate_invasion() tables for the two sets.
#' @param restrict method regex defining the candidate family (as passed to the caller).
#' @param exclude regex of scored-but-not-stage-1 methods (INVASION_SELECTION_EXCLUDE: the
#'   sensitivity arms and the time-varying forms).
#' @param featured the DEPLOYED pick, flagged separately: the inner-fold selection need not
#'   return it, and publishing only the inner pick's row would report optimism for a model
#'   nothing ships.
#' @return data.frame, one row per (method, horizon), or NULL.
.heldout_field <- function(inner_eval, outer_eval, restrict = "^Bayes",
                           exclude = get0("INVASION_SELECTION_EXCLUDE", ifnotfound = NULL),
                           selected = NA_character_, featured = NA_character_) {
  if (is.null(inner_eval) || is.null(outer_eval)) return(NULL)
  .ex_ok <- !is.null(exclude) && length(exclude) == 1L && !is.na(exclude) && nzchar(exclude)
  .arms <- get0("SENSITIVITY_ARM_REGEX", ifnotfound = NULL)
  keep <- function(e) {
    if (!is.null(restrict)) e <- e[grepl(restrict, e$method), , drop = FALSE]
    # The FIELD is the stage-1 candidate set. A time-varying form is excluded from it, but the
    # selected / featured model keeps its row even when it is a form (stage 2 can pick one), so
    # its own drop is still reported; `in_field` marks the rows the field statistics use.
    if (.ex_ok) {
      own <- e$method %in% c(selected, featured) &
             !(if (!is.null(.arms) && nzchar(.arms)) grepl(.arms, e$method) else FALSE)
      e <- e[!grepl(exclude, e$method) | own, , drop = FALSE]
    }
    # Ensembles are not candidates for the featured single model (run_all.R builds the pick
    # from .bayes_singles), so they must not sit in the field either — the field IS the
    # candidate set, and a comparator the pipeline cannot select would shift the median.
    e <- e[!grepl("-ens-", e$method), , drop = FALSE]
    # A rank-only comparator has no calibrated scale; it is still rankable, but it is also
    # ineligible to be featured, so it is not part of the field for this purpose.
    if ("prob_calibrated" %in% names(e)) e <- e[!(e$prob_calibrated %in% FALSE), , drop = FALSE]
    e
  }
  ei <- keep(inner_eval); eo <- keep(outer_eval)
  if (!nrow(ei) || !nrow(eo)) return(NULL)
  # ON THE SELECTION'S AXIS. Winner's curse can only be measured on the axis that was used to
  # select: per-date AUC-PR skill on the post-burn-in rounds (the sel_* columns). The pooled
  # all-round skill is a different quantity (it rewards a probability LEVEL; the inner set also
  # carries the burn-in rounds), so measured on it the selection's own optimism was invisible.
  .sel_cols <- c("sel_auc_pr_date", "sel_auc_pr_skill_date", "sel_prev_date", "sel_n_invasions")
  if (!all(.sel_cols %in% names(ei)) || !all(.sel_cols %in% names(eo)))
    stop("[heldout-field] the evaluations carry no selection-axis columns (",
         paste(setdiff(.sel_cols, intersect(names(ei), names(eo))), collapse = ", "), ").",
         call. = FALSE)
  .to_sel <- function(e) { e$auc_pr <- e$sel_auc_pr_date; e$auc_pr_skill <- e$sel_auc_pr_skill_date
                           e$base_rate <- e$sel_prev_date; e$n_invasions <- e$sel_n_invasions; e }
  ei <- .to_sel(ei); eo <- .to_sel(eo)
  .in_field <- function(m) if (.ex_ok) !grepl(exclude, m) else rep(TRUE, length(m))
  m <- merge(
    ei[, c("method", "horizon", "auc_pr", "auc_pr_skill", "base_rate", "n_invasions")],
    eo[, c("method", "horizon", "auc_pr", "auc_pr_skill", "base_rate", "n_invasions")],
    by = c("method", "horizon"), suffixes = c("_inner", "_outer"))
  if (!nrow(m)) return(NULL)
  m <- do.call(rbind, lapply(split(m, m$horizon), function(d) {
    # The inner base rate is one number per horizon (every method is scored on the same
    # shared support), so taking [1] is safe; assert it rather than assume it.
    bi <- unique(signif(d$base_rate_inner, 10))
    if (length(bi) != 1L)
      warning(sprintf("[heldout-field] h=%s: %d distinct inner base rates in the field; the common-base drop is computed against the first.",
                      d$horizon[1], length(bi)), call. = FALSE)
    d$drop        <- d$auc_pr_skill_inner - d$auc_pr_skill_outer
    d$drop_common <- (d$auc_pr_inner - d$auc_pr_outer) / d$base_rate_inner[1]
    # Ranks AMONG THE FIELD (ties "min"): each row is ranked against the stage-1 candidates
    # other than itself, so a selected / featured form's row neither enters nor shifts the
    # candidates' ranks. NA skill -> NA rank.
    f <- .in_field(d$method)
    .rk <- function(x) vapply(seq_along(x), function(i)
      if (is.na(x[i])) NA_integer_ else 1L + sum(f[-i] & x[-i] > x[i], na.rm = TRUE), integer(1))
    d$inner_rank  <- .rk(d$auc_pr_skill_inner)
    d$outer_rank  <- .rk(d$auc_pr_skill_outer)
    d[order(d$inner_rank), , drop = FALSE]
  }))
  m$is_selected <- m$method == selected
  m$is_featured <- m$method == featured
  m$in_field <- .in_field(m$method)
  rownames(m) <- NULL
  m
}

#' Collapse .heldout_field() to the three numbers that belong in a caption.
#' @return data.frame, one row per horizon.
.heldout_field_summary <- function(field, selected = NA_character_) {
  if (is.null(field) || !nrow(field)) return(NULL)
  do.call(rbind, lapply(split(field, field$horizon), function(d) {
    s <- d[d$method == selected, , drop = FALSE]
    # The field statistics are over the stage-1 candidates only (a selected time-varying form
    # carries its own row but is not part of the field: see .heldout_field()).
    if ("in_field" %in% names(d)) d <- d[d$in_field %in% TRUE, , drop = FALSE]
    # Spearman needs a spread of at least three candidates to mean anything.
    rho <- if (nrow(d) >= 3L)
      suppressWarnings(stats::cor(d$auc_pr_skill_inner, d$auc_pr_skill_outer,
                                  method = "spearman")) else NA_real_
    data.frame(
      horizon = d$horizon[1], n_candidates = nrow(d), selected = selected,
      # What the selected model's ranks are "out of": the field, plus the selected model itself
      # when it is not a field member (a time-varying form picked in stage 2).
      selected_rank_of = nrow(d) + as.integer(nrow(s) > 0L && "in_field" %in% names(s) &&
                                                !isTRUE(s$in_field[1])),
      n_events_inner = d$n_invasions_inner[1], n_events_outer = d$n_invasions_outer[1],
      inner_base_rate = d$base_rate_inner[1], outer_base_rate = d$base_rate_outer[1],
      # THE WINDOW EFFECT: what every candidate loses on the outer origins.
      field_median_drop = stats::median(d$drop, na.rm = TRUE),
      field_median_drop_common = stats::median(d$drop_common, na.rm = TRUE),
      field_min_drop_common = min(d$drop_common, na.rm = TRUE),
      field_max_drop_common = max(d$drop_common, na.rm = TRUE),
      # THE SELECTION EFFECT: what the picked model loses OVER AND ABOVE the field.
      selected_drop = if (nrow(s)) s$drop[1] else NA_real_,
      selected_drop_common = if (nrow(s)) s$drop_common[1] else NA_real_,
      selected_excess = if (nrow(s)) s$drop[1] - stats::median(d$drop, na.rm = TRUE) else NA_real_,
      selected_excess_common = if (nrow(s))
        s$drop_common[1] - stats::median(d$drop_common, na.rm = TRUE) else NA_real_,
      selected_rank_inner = if (nrow(s)) s$inner_rank[1] else NA_integer_,
      selected_rank_outer = if (nrow(s)) s$outer_rank[1] else NA_integer_,
      # DOES THE ORDERING SURVIVE? The level falling for everyone is survivable; the ranking
      # inverting is not, because the ranking is what selection consumes.
      rank_correlation = rho,
      stringsAsFactors = FALSE)
  }))
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
#' @param n_holdout number of most-recent origins held out (default 3).
#' @param restrict optional method regex to select WITHIN (e.g. "^Bayes").
#' @param horizon horizon at which the skill/gap is reported (default 1).
#' @return list(selected, inner_skill, outer_skill, optimism_gap, inner_eval,
#'   outer_eval, inner_cutoffs, outer_cutoffs), or NULL if too few origins.
evaluate_invasion_heldout <- function(lfo_results, n_holdout = 3L,
                                      restrict = "^Bayes", horizon = 1L,
                                      featured = NA_character_) {
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
  # The featured pick is a SINGLE model (run_all.R selects from the non-ensemble Bayesian models)
  # by the TWO-STAGE rule, so the inner pick is made the same way, from the same field. Stage 2
  # can only see the forms that were cross-validated (those of the full-data stage-1 pick).
  sel_tbl <- sel_tbl[!grepl("-ens-", sel_tbl$method), , drop = FALSE]
  sel <- select_featured_two_stage(sel_tbl, restrict = restrict %||% ".")$featured
  if (length(sel) != 1L || is.na(sel)) return(NULL)
  gk <- function(ev, m, col = "auc_pr_skill") {
    r <- ev[[col]][ev$method == m & ev$horizon == horizon]
    if (length(r)) r[1] else NA_real_
  }
  # THE SELECTION'S AXIS (per-date AUC-PR skill, post-burn-in rounds) for every published
  # number below; the pooled all-round skill is kept, labelled *_pooled.
  inner_skill <- gk(inner_eval, sel, "sel_auc_pr_skill_date")
  outer_skill <- gk(outer_eval, sel, "sel_auc_pr_skill_date")
  inner_skill_pooled <- gk(inner_eval, sel); outer_skill_pooled <- gk(outer_eval, sel)
  # BASE-RATE DRIFT. inner_eval and outer_eval are two independent evaluate_invasion() calls,
  # so each computes its OWN base_rate, and auc_pr_skill = auc_pr / base_rate. The published
  # `optimism_gap` is therefore a difference of ratios with DIFFERENT denominators: it mixes
  # genuine selection optimism with the fact that the outer window is simply a harder or easier
  # problem. Measured 2026-09-19 on the shipped folds (Bayes-M14-fill-med, h=1): inner base rate
  # 0.006837 (AUC-PR 0.3221, skill 47.12), outer 0.004301 (AUC-PR 0.1485, skill 34.52) — a 37%
  # denominator shift. The reported gap is 12.60; on a COMMON denominator it is 25.40, i.e. the
  # headline number understates the optimism roughly twofold.
  # Both are now returned, with the inputs, so the reader can see which is which.
  # Score EVERY eligible candidate on the same two sets, so the window effect the whole
  # field shares can be subtracted from the selected model's drop.
  .fld <- tryCatch(.heldout_field(inner_eval, outer_eval, restrict = restrict,
                                  selected = sel, featured = featured),
                   error = function(e) {
                     warning("[heldout] field decomposition failed: ", conditionMessage(e),
                             call. = FALSE); NULL })
  .fld_sum <- tryCatch(.heldout_field_summary(.fld, selected = sel),
                       error = function(e) NULL)
  if (!is.null(.fld_sum)) {
    .r <- .fld_sum[.fld_sum$horizon == horizon, , drop = FALSE]
    if (nrow(.r)) message(sprintf(
      paste0("[heldout] h=%s: field median drop %.1f (common base), selected %s drops %.1f ",
             "=> EXCESS attributable to selection %.1f; rank %d -> %d of %d; rho = %.2f"),
      horizon, .r$field_median_drop_common[1], sel, .r$selected_drop_common[1],
      .r$selected_excess_common[1], .r$selected_rank_inner[1], .r$selected_rank_outer[1],
      .r$selected_rank_of[1], .r$rank_correlation[1]))
  }
  # The selection axis's numerator and denominator (event-weighted per-date AP and base rate).
  inner_base <- gk(inner_eval, sel, "sel_prev_date"); outer_base <- gk(outer_eval, sel, "sel_prev_date")
  inner_ap   <- gk(inner_eval, sel, "sel_auc_pr_date"); outer_ap <- gk(outer_eval, sel, "sel_auc_pr_date")
  list(selected = sel, horizon = horizon, n_holdout = n_holdout,
       inner_skill = inner_skill, outer_skill = outer_skill,
       inner_skill_pooled = inner_skill_pooled, outer_skill_pooled = outer_skill_pooled,
       optimism_gap_pooled = inner_skill_pooled - outer_skill_pooled,
       # Each set scored against its own base rate.
       optimism_gap = inner_skill - outer_skill,
       # On the INNER base rate throughout, so the denominator cannot move the answer.
       optimism_gap_common_base = if (is.finite(inner_base) && inner_base > 0)
         (inner_ap - outer_ap) / inner_base else NA_real_,
       inner_base_rate = inner_base, outer_base_rate = outer_base,
       inner_auc_pr = inner_ap, outer_auc_pr = outer_ap,
       base_rate_ratio = if (is.finite(inner_base) && inner_base > 0)
         outer_base / inner_base else NA_real_,
       inner_eval = inner_eval, outer_eval = outer_eval,
       inner_cutoffs = inner_cut, outer_cutoffs = outer_cut,
       # THE FIELD, at every horizon. `optimism_gap` above is one model's inner-minus-outer
       # difference and cannot tell a window effect from a selection effect; these two
       # separate them. See .heldout_field(). Both are NULL-safe: a caller that only wants
       # the scalar gap is unaffected.
       field = .fld, field_summary = .fld_sum)
}

# ---------------------------------------------------------------------------
# Clean invasion LFO-CV (same folds for all models; correct at-risk outcome)
# ---------------------------------------------------------------------------

#' Candidate fold cutoffs, and the horizons complete at each.
#'
#' THE SINGLE DEFINITION of the LFO fold grid. run_invasion_lfo() builds its folds from this,
#' and write_run_info.R calls it to record the PLANNED grid — previously a hand-written
#' re-implementation in that file, which had already drifted in idiom and would drift in
#' substance the moment the admission rule changed.
#'
#' @param weeks sorted unique week_start dates of the outbreak grid.
#' @param horizons forecast horizons (weeks).
#' @param analysis_date the as-of date.
#' @param min_train_weeks minimum training weeks before the first fold.
#' @param min_eval_age_days how old an outcome window must be to be admitted.
#' @return list(cutoff = Date, horizons = integer vector) per admitted cutoff, oldest first.
#'   The caller applies any `min_atrisk_events` rule, which needs the outcome data.
invasion_fold_candidates <- function(weeks, horizons, analysis_date,
                                     min_train_weeks = get0("LFO_MIN_TRAINING_WEEKS",
                                                            ifnotfound = 3L),
                                     min_eval_age_days = get0("LFO_MIN_EVAL_AGE_DAYS",
                                                              ifnotfound = 13L)) {
  weeks <- sort(unique(as.Date(weeks)))
  if (length(weeks) < min_train_weeks) return(list())
  # The FARTHEST horizon governs admission under the common fold grid (the default), so every
  # admitted round is scored at every horizon; the SHORTEST governs under the Q1 augmentation.
  h_admit <- if (isTRUE(get0("LFO_COMMON_FOLD_GRID", ifnotfound = TRUE)))
               max(horizons) else min(horizons)
  cand <- weeks[weeks >= weeks[min_train_weeks] &
                  (as.numeric(analysis_date - (weeks + h_admit * 7 + 7)) >= min_eval_age_days)]
  out <- list()
  for (ci in seq_along(cand)) {
    cut <- cand[ci]                     # index, NOT `for (x in <Date>)`, which unclasses
    complete_h <- horizons[
      as.numeric(analysis_date - (cut + horizons * 7 + 7)) >= min_eval_age_days]
    if (!length(complete_h)) next
    out[[length(out) + 1L]] <- list(cutoff = cut, horizons = complete_h)
  }
  out
}

#' Leave-future-out cross-validation focused on the invasion task.
#'
#' Every model is run on the SAME principled folds. Weeks are labelled by their first day
#' (WEEK_ANCHOR grid). For fold cutoff C -- the first day of the last training week:
#'   - origin       = C + 6 (LFO_ORIGIN_OFFSET_DAYS), the last day of that week: the forecast
#'                    is issued on this day, as the deployed forecast is on ANALYSIS_DATE;
#'   - training data = the line list AS KNOWN ON THE ORIGIN (asof_counts_source(),
#'                    22_daily_reissue.R). Under LFO_TRAINING_SOURCE = "vintage" (default) the
#'                    archived snapshot of the origin is re-aggregated and nowcast at its own
#'                    date with the extract truncation fitted on that snapshot's history; a round
#'                    with no snapshot (or none supporting a real-time truncation fit) re-censors
#'                    the final line list to the origin and is nowcast with `delay`. Each row
#'                    records which (`training_source`, `training_fallback`, `vintage_date`);
#'   - at-risk set  = zones with no confirmed case in weeks <= C (onset-dated, final line list);
#'   - outcome(h)   = zone is first confirmed in weeks C+7 .. C+7h, i.e. days (C+6, C+6+7h]
#'                    (cumulative window, matching the models' cumulative invasion probability).
#'
#'     OVERLAPPING WINDOWS AT h >= 2 — read this before quoting an h=2 number. With weekly
#'     cutoffs the h=1 windows are disjoint, so each invasion is scored exactly once. The h=2
#'     windows OVERLAP by a week, so a zone invaded in week W is still at-risk at both
#'     C = W-14 and C = W-7 and is a positive in both: at h=2 the number of positive rows
#'     exceeds the number of distinct invasions.
#'
#'     Each ROW is still a legitimate (forecast, outcome) pair — "will this zone be invaded
#'     within 2 weeks, as of C?" asked at two different origins — and every metric is
#'     internally consistent on those rows, including the h=2 base rate, which is the rate
#'     among the rows actually scored. What is NOT true at h=2 is independence: the positive
#'     rows carry fewer distinct events' worth of information. Interval estimates are
#'     protected (the bootstrap
#'     resamples ZONES, the correlated unit); any reading that treats n_invasions at h=2 as
#'     an effective sample size is not. `n_invasion_events` is reported alongside
#'     `n_invasions` so the two cannot be confused.
#' A fold is admitted at horizon h once its outcome window has been closed for
#' `min_eval_age_days`; rows younger than `reliable_eval_age_days` are labelled
#' eval_reliable = FALSE (scored, but excluded from the deployment delta).
#'
#' @param zone_week_raw raw zone-week counts (health_zone, week_start, confirmed,
#'   suspected), outbreak period.
#' @param models named list of forecasters, each `function(zw_nc, t_idx,
#'   horizons, cutoff)` returning health_zone, horizon, mu_forecast, p_invasion,
#'   method.
#' @param nowcast_fn function(zone_week, analysis_date) -> nowcast-corrected;
#'   default apply_nowcast_correction (must be sourced).
#' @param linelist the final line list (dat$ll). Required for training data as known at the
#'   origin; NULL slices the FINAL weekly counts instead (revision-leaky; rows are labelled
#'   training_source = "final_slice").
#' @return pooled tibble across folds x models with p_invasion, is_new_invasion,
#'   was_active_before, fold_id, cutoff, method, horizon, eval_age_days, eval_reliable,
#'   training_source, training_fallback, vintage_date.
#' @param delay the truncation spec for the RE-CENSORED rounds only -- in the pipeline the
#'   fitted as-of truncation (epinow2_truncation(regime = "asof"), passed by run_all.R). A
#'   snapshot round uses its own real-time extract truncation instead. NULL falls through to
#'   whatever `nowcast_fn` resolves.
run_invasion_lfo <- function(zone_week_raw, models, horizons = LFO_HORIZONS,
                             min_train_weeks = get0("LFO_MIN_TRAINING_WEEKS",
                                                    ifnotfound = 3L),
                             min_atrisk_events = 1L,
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
                             # 00_config.R (LFO_MIN_EVAL_AGE_DAYS), and the truncation is
                             # handled where it belongs — reported rather than avoided. Every returned row carries `eval_age_days` (the
                             # age of ITS OWN outcome window, per horizon) and
                             # `eval_reliable` (age >= LFO_RELIABLE_EVAL_AGE_DAYS), so any
                             # consumer can restrict to settled folds; run_all.R publishes
                             # exactly that restriction as a sensitivity table.
                             #
                             # Truncation biases the newest folds' skill DOWNWARD and their
                             # calibration-in-the-large UPWARD. It cannot bias the ranking
                             # between models: all are scored on identical rows and truth.
                             min_eval_age_days = get0("LFO_MIN_EVAL_AGE_DAYS", ifnotfound = 13L),
                             # Outcome-window age at which a fold's invasion count is treated
                             # as settled. Purely a LABEL: it never admits or excludes a fold.
                             reliable_eval_age_days = get0("LFO_RELIABLE_EVAL_AGE_DAYS",
                                                           ifnotfound = 25L),
                             analysis_date = ANALYSIS_DATE,
                             nowcast_fn = NULL, linelist = NULL, delay = NULL) {
  if (is.null(nowcast_fn)) nowcast_fn <- get0("apply_nowcast_correction")
  stopifnot(!is.null(nowcast_fn))
  # ---- TRUNCATION USED BY EACH ROUND'S NOWCAST ------------------------------------------
  # A snapshot round is nowcast with the extract truncation fitted on that snapshot's own
  # history (vintage_truncation(), 22_daily_reissue.R), i.e. what was estimable on its origin.
  # A re-censored round (no snapshot, or none supporting a real-time fit -- the early rounds)
  # uses `delay`, the as-of truncation fitted on the final line list: that is the one input to
  # those rounds a real-time system would not have had, and the rows say which rounds they are.
  if (!is.null(delay))
    message(sprintf("[invasion_lfo] Truncation for re-censored rounds: %s (mean %.2f d, %s)",
                    delay$family, delay$mean, delay$estimator))
  else
    message("[invasion_lfo] No truncation passed: re-censored rounds use nowcast_fn's own default.")
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

  # ONE definition of the fold grid — see invasion_fold_candidates() above, which
  # write_run_info.R also calls so run_info.json can never report a planned grid the
  # cross-validation did not use.
  #
  # Candidate cutoffs. WHICH HORIZON GOVERNS ADMISSION is set by LFO_COMMON_FOLD_GRID
  # (00_config.R, which carries the full rationale and the reason no single age threshold can
  # align the two horizons' grids):
  #   TRUE  (default) — the FARTHEST horizon must be complete, so every admitted cutoff is
  #                     scored at EVERY horizon and the horizons share one set of origins.
  #   FALSE           — the Q1 fold augmentation: the SHORTEST horizon governs, adding the
  #                     recent cutoffs where h=1 is evaluable and h=2 is not.
  # The per-horizon `complete_h` guard below is kept under BOTH settings: it is what makes the
  # augmentation leak-free, and under the common grid it is simply always satisfied.
  .cand <- invasion_fold_candidates(weeks, horizons, analysis_date,
                                    min_train_weeks = min_train_weeks,
                                    min_eval_age_days = min_eval_age_days)
  folds <- list(); fid <- 0L
  for (ci in seq_along(.cand)) {
    cut        <- .cand[[ci]]$cutoff
    complete_h <- .cand[[ci]]$horizons
    hmax_c     <- max(complete_h)
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

  # SNAPSHOT MODE (the default): load every snapshot the rounds will read -- each round's origin
  # and every rolling-predictor transition (week ends up to the last cutoff) -- and fit each
  # snapshot's real-time truncation, in THIS process, so the forked workers inherit them and
  # the Stan fits run once. Then say which rounds read a snapshot and which re-censor, and why.
  if (identical(get0("LFO_TRAINING_SOURCE", ifnotfound = "vintage"), "vintage") &&
      !is.null(linelist) && exists("prefetch_vintages", mode = "function")) {
    .last_cut <- max(vapply(folds, function(f) as.numeric(f$cutoff), 0))
    .pv <- prefetch_vintages(weeks[as.numeric(weeks) <= .last_cut] + 6L)
    .orig <- vapply(folds, function(f) as.numeric(f$cutoff + 6L), 0)
    .row  <- match(.orig, as.numeric(.pv$issue_date))
    .snap <- !is.na(.row)
    .read <- .snap & is.finite(.pv$truncation_mean[ifelse(.snap, .row, 1L)])
    .lab  <- function(i) paste(vapply(folds[i], function(f) format(f$cutoff + 6L, "%m-%d"), ""),
                               collapse = ",")
    message(sprintf(paste0("[invasion_lfo] training data AS KNOWN AT EACH ORIGIN: %d of %d rounds ",
                           "read the archived snapshot with its real-time truncation (origins %s); ",
                           "%d re-censor the final line list because no snapshot exists (%s) and ",
                           "%d because the archive did not yet support a real-time truncation ",
                           "fit (%s). %d distinct snapshots loaded."),
                    sum(.read), length(folds), .lab(.read), sum(!.snap), .lab(!.snap),
                    sum(.snap & !.read), .lab(.snap & !.read), length(unique(.pv$folder))))
  }

  out <- list()
  for (fi in seq_along(folds)) {
    f <- folds[[fi]]; cut <- f$cutoff
    # Training counts AS KNOWN ON THE ORIGIN (cut + 6): asof_counts_source() (22_daily_reissue.R)
    # reads the archived snapshot of that date, nowcast with the truncation estimated on it, or
    # -- where no usable snapshot exists -- re-censors the final line list to the origin. The
    # final-count slice is used only when no line list is given or the reconstruction fails
    # (warned, and labelled training_source = "final_slice").
    train_raw <- NULL
    .fold_src <- list(delay = delay, nowcast_date = cut + 6, source = "final_slice",
                      vintage_date = as.Date(NA), why = "no line list")
    if (!is.null(linelist) && exists("reaggregate_asof")) {
      .src <- tryCatch(
        if (exists("asof_counts_source", mode = "function"))
          asof_counts_source(linelist, zones, cut + 6, week_spine = weeks[weeks <= cut],
                             delay_recensor = delay)
        else list(raw = reaggregate_asof(linelist, zones, cut + 6,
                                         week_spine = weeks[weeks <= cut]),
                  delay = delay, nowcast_date = cut + 6, source = "recensor",
                  vintage_date = as.Date(NA), why = "22_daily_reissue.R not loaded"),
        error = function(e) {
          # Warn LOUDLY: a silent NULL here drops this fold onto the final-count (revision-leaky)
          # path even under LEAKAGE_FREE_LFO, invisibly weakening the leakage guard for that fold.
          warning(sprintf("[lfo] fold %s: reaggregate_asof failed (%s); this fold falls back to the final-count (revision-leaky) training slice.",
                          format(cut), conditionMessage(e)), call. = FALSE); NULL })
      if (!is.null(.src)) {
        train_raw <- .src$raw %>% dplyr::filter(week_start <= cut)
        .fold_src <- .src[c("delay", "nowcast_date", "source", "vintage_date", "why")]
      }
    }
    if (is.null(train_raw))
      train_raw <- zone_week_raw %>% dplyr::filter(week_start <= cut)
    # FOLD ORIGIN = cut + 6, the LAST DAY OF THE LAST TRAINING WEEK -- mirroring deployment,
    # where the grid is anchored so the final week ENDS on ANALYSIS_DATE and the as-of date IS
    # that day. The horizon targets are defined relative to `cut`. Each round's truncation
    # (.fold_src$delay) is passed explicitly; a snapshot is nowcast at ITS OWN date, and a
    # snapshot older than the origin (by at most LFO_VINTAGE_MAX_GAP_DAYS) leaves the origin
    # week's last days unobserved, which the per-day completeness weights count as such.
    train_nc  <- tryCatch(if (is.null(.fold_src$delay))
                            nowcast_fn(train_raw, analysis_date = .fold_src$nowcast_date)
                          else nowcast_fn(train_raw, analysis_date = .fold_src$nowcast_date,
                                          delay = .fold_src$delay),
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
                      eval_reliable = eval_age_days >= reliable_eval_age_days,
                      # WHERE THIS ROUND'S TRAINING COUNTS CAME FROM: "vintage" (the archived
                      # snapshot of `vintage_date`, nowcast with the truncation estimated on
                      # that date), "recensor" (the final line list censored to the origin;
                      # `training_fallback` says why) or "final_slice" (no line list).
                      training_source = .fold_src$source,
                      training_fallback = .fold_src$why %||% NA_character_,
                      vintage_date = .fold_src$vintage_date)
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
