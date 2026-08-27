# =============================================================================
# 19_spacetime_eval.R — Spatiotemporal Evaluation of Invasion Forecasts (Q5)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# The pooled metrics in 16_invasion_eval.R answer "how well does each model rank
# at-risk zones overall?". This module answers the SPATIOTEMPORAL questions that
# drive iterative model improvement:
#
#   * skill-over-time  — does discrimination improve as the outbreak matures?
#                        (per-fold AUC-PR skill / hit-rate@K vs forecast date)
#   * lead-time        — how many weeks AHEAD of its first case is each newly
#                        affected zone flagged into the top-K? (anticipation)
#   * spatial error    — WHERE do the misses and chronic false alarms sit?
#                        (per-zone pooled predicted-vs-observed, joined to region)
#
# Depends on the scoring primitives (.auc_pr, .auc_roc) defined in
# 16_invasion_eval.R, which run_all.R sources first; sourced here too for
# standalone use.
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))
if (!exists(".auc_pr"))
  source(file.path(here::here(), "spatiotemporal", "16_invasion_eval.R"))
suppressPackageStartupMessages({ library(tidyverse) })

#' First confirmed week per zone (for lead-time), from raw zone-week counts.
first_case_from_zone_week <- function(zone_week_raw) {
  zone_week_raw %>%
    dplyr::group_by(health_zone, week_start) %>%
    dplyr::summarise(confirmed = sum(confirmed, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(confirmed > 0) %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(first_wk = min(week_start), .groups = "drop")
}

# ---------------------------------------------------------------------------
# 1. Skill over time — per-fold discrimination/targeting, one row per fold
# ---------------------------------------------------------------------------

#' Per-(method, horizon, fold) skill: AUC-PR + skill, hit/precision@K, rank-of-truth.
#'
#' A fold contributes only if it has >=1 invasion event (AUC-PR undefined
#' otherwise); such folds are still reported with NA discrimination so the
#' event-free periods are visible.
#'
#' @return tibble(method, horizon, fold_id, cutoff, n_atrisk, n_events,
#'   base_rate, auc_pr, auc_pr_skill, hit_at_k, prec_at_k, mean_rank_of_truth).
spatiotemporal_skill <- function(lfo_results, k = 5L) {
  stopifnot(all(c("method", "horizon", "fold_id", "cutoff", "health_zone",
                  "p_invasion", "is_new_invasion") %in% names(lfo_results)))
  d <- lfo_results %>%
    dplyr::filter(is.finite(p_invasion))
  if ("was_active_before" %in% names(d))
    d <- d %>% dplyr::filter(!as.logical(was_active_before) %in% TRUE)

  d %>%
    dplyr::group_by(method, horizon, fold_id, cutoff) %>%
    dplyr::group_modify(function(g, ...) {
      y <- as.integer(g$is_new_invasion); p <- g$p_invasion
      n <- length(y); n_pos <- sum(y, na.rm = TRUE)   # na.rm: never NA the guard
      o <- order(p, decreasing = TRUE); yo <- y[o]
      topk <- head(yo, k)
      ap <- .auc_pr(p, y)
      rk <- rank(-p, ties.method = "max")             # pessimistic on zero-ties (1 = highest)
      tibble::tibble(
        n_atrisk = n, n_events = n_pos, base_rate = n_pos / max(n, 1),
        auc_pr = ap, auc_roc = .auc_roc(p, y),
        auc_pr_skill = if (n_pos > 0) ap / (n_pos / n) else NA_real_,
        hit_at_k  = if (n_pos > 0) as.integer(sum(topk) > 0) else NA_integer_,
        prec_at_k = if (n_pos > 0) sum(topk) / k else NA_real_,
        recall_at_k = if (n_pos > 0) sum(topk) / n_pos else NA_real_,
        mean_rank_of_truth = if (n_pos > 0) mean(rk[y == 1]) else NA_real_)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cutoff = as.Date(cutoff)) %>%
    dplyr::arrange(method, horizon, cutoff)
}

# ---------------------------------------------------------------------------
# 2. Lead-time — weeks of anticipation before a zone's first case
# ---------------------------------------------------------------------------

#' For each newly affected zone, how many weeks ahead it entered the top-K.
#'
#' Using the h=1 forecast at every fold, rank the at-risk zones; for a zone whose
#' first case falls in week W, a fold at cutoff C < W gives lead (W-C)/7 weeks.
#' The zone's anticipation is the LARGEST lead at which it was in the top-K — i.e.
#' the earliest the model would have raised the alarm. best_rank is its single
#' best (lowest) rank across those pre-invasion folds.
#'
#' @param first_case_week tibble(health_zone, first_wk) from
#'   first_case_from_zone_week().
#' @return tibble(health_zone, first_wk, n_prefolds, best_rank,
#'   anticipation_weeks, ever_topk), one row per zone flagged at least once.
lead_time_analysis <- function(lfo_results, first_case_week, k = 5L, horizon = 1L,
                               method = NULL) {
  d <- lfo_results %>% dplyr::filter(horizon == !!horizon, is.finite(p_invasion))
  # Ranks are within-fold across the at-risk set of ONE model; pooling methods
  # would rank over n_methods×n_zones rows and corrupt every rank. Always reduce
  # to a single method (the caller's, or the first present with a warning).
  if (is.null(method)) {
    ms <- unique(d$method)
    if (length(ms) > 1)
      warning("[lead_time] method=NULL with ", length(ms),
              " methods present; using '", ms[1], "'.")
    method <- ms[1]
  }
  d <- d %>% dplyr::filter(method == !!method)
  if (nrow(d) == 0) return(NULL)
  d <- d %>%
    dplyr::group_by(fold_id) %>%
    # Pessimistic ranking for ties: a zone buried in a large group tied at
    # p_invasion == 0 (no mobility pathway) must NOT be credited as "detected in
    # top-K". min_rank gives every tied zone the group minimum (top) rank; use
    # ties.method="max" so tied zones take the WORST position in their run.
    dplyr::mutate(rank = rank(-p_invasion, ties.method = "max")) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(first_case_week, by = "health_zone") %>%
    dplyr::mutate(cutoff = as.Date(cutoff),
                  lead_weeks = as.numeric(first_wk - cutoff) / 7,
                  detected = rank <= k) %>%
    dplyr::filter(lead_weeks > 0)                # folds strictly before first case
  if (nrow(d) == 0) return(NULL)
  d %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(
      first_wk = dplyr::first(first_wk),
      n_prefolds = dplyr::n(),
      best_rank = min(rank),
      anticipation_weeks = { det <- lead_weeks[detected]
        if (length(det)) max(det) else 0 },
      ever_topk = any(detected),
      .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(anticipation_weeks), best_rank)
}

# ---------------------------------------------------------------------------
# 3. Spatial error — per-zone pooled predicted-vs-observed
# ---------------------------------------------------------------------------

#' Per-zone pooled behaviour across folds: how often at-risk, how often invaded,
#' mean/max predicted probability, best rank; flags chronic false alarms (high
#' mean risk, never invaded) and misses (invaded but low mean risk). Joined to
#' province for spatial mapping.
#'
#' @return tibble(health_zone, province, n_atrisk_folds, n_invasions, mean_p,
#'   max_p, best_rank, flag) sorted by mean_p.
zone_spatial_error <- function(lfo_results, province_map = NULL, horizon = 1L,
                               method = NULL, k = 5L) {
  d <- lfo_results %>% dplyr::filter(horizon == !!horizon, is.finite(p_invasion))
  # Single-method only: per-zone counts (n_atrisk_folds, n_invasions) and
  # within-fold ranks are meaningless if pooled across methods (each event/zone
  # would be multiplied by the number of methods).
  if (is.null(method)) {
    ms <- unique(d$method)
    if (length(ms) > 1)
      warning("[zone_spatial_error] method=NULL with ", length(ms),
              " methods present; using '", ms[1], "'.")
    method <- ms[1]
  }
  d <- d %>% dplyr::filter(method == !!method)
  if (nrow(d) == 0) return(NULL)
  d <- d %>%
    dplyr::group_by(fold_id) %>%
    # Pessimistic ranking for ties: a zone buried in a large group tied at
    # p_invasion == 0 (no mobility pathway) must NOT be credited as "detected in
    # top-K". min_rank gives every tied zone the group minimum (top) rank; use
    # ties.method="max" so tied zones take the WORST position in their run.
    dplyr::mutate(rank = rank(-p_invasion, ties.method = "max")) %>%
    dplyr::ungroup()
  z <- d %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(
      n_atrisk_folds = dplyr::n(),
      n_invasions = sum(is_new_invasion, na.rm = TRUE),
      mean_p = mean(p_invasion, na.rm = TRUE),
      max_p  = max(p_invasion, na.rm = TRUE),
      best_rank = min(rank, na.rm = TRUE),
      .groups = "drop")
  # thresholds relative to the pooled distribution
  hi <- stats::quantile(z$mean_p, 0.9, na.rm = TRUE)
  z <- z %>% dplyr::mutate(flag = dplyr::case_when(
    n_invasions > 0 & best_rank <= k        ~ "hit",
    n_invasions > 0 & best_rank >  k        ~ "miss",
    n_invasions == 0 & mean_p >= hi         ~ "false_alarm",
    TRUE                                     ~ "correct_negative"))
  if (!is.null(province_map))
    # distinct on the join key so a duplicated `nom` in province_map cannot inflate
    # (duplicate) a zone's single spatial-error row.
    z <- z %>% dplyr::left_join(dplyr::distinct(province_map, nom, .keep_all = TRUE),
                                by = c("health_zone" = "nom"))
  else z$province <- NA_character_
  z %>% dplyr::arrange(dplyr::desc(mean_p))
}

message("[spacetime_eval] 19_spacetime_eval.R loaded — skill-over-time, lead-time, spatial error (Q5).")
