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
#' @return tibble(method, horizon, fold_id, cutoff, n_atrisk, n_events, base_rate,
#'   auc_pr, auc_roc, auc_pr_skill, auc_pr_lift, hit_at_k, prec_at_k, recall_at_k,
#'   mean_rank_of_truth). Every column is a within-fold ranking quantity, so it is the same on
#'   the raw and the recalibrated probability scale (one monotone transform per fold).
#' @param common_support restrict to the (fold x zone) cells EVERY method covers
#'   (invasion_common_cells(), 16_invasion_eval.R), so per-fold skill is comparable ACROSS
#'   methods. TRUE by default and it must stay that way for anything published: this function
#'   feeds skill_over_time.csv and plot_skill_over_time(), which draws the featured model plus
#'   the five best other methods BY THIS MEAN (sensitivity arms and non-featured time-varying
#'   forms excluded). Scored on native rows, a method that covers an origin or a cell
#'   the others lack -- the structural baselines use no design, so they keep a round in which
#'   a Bayesian fit produced no forecast -- averages its per-fold skill over a different set
#'   of origins, and the early, easier origins can move it up the ranking on that alone. (This
#'   once put Gravity-B4 into the published panel on an origin the Bayesian models had lost.)
#'   It is the same defect fixed in 16_invasion_eval.R, which did not reach here because this
#'   module scores independently.
spatiotemporal_skill <- function(lfo_results, k = 5L, common_support = TRUE) {
  stopifnot(all(c("method", "horizon", "fold_id", "cutoff", "health_zone",
                  "p_invasion", "is_new_invasion") %in% names(lfo_results)))
  d <- lfo_results %>%
    dplyr::filter(is.finite(p_invasion))
  if ("was_active_before" %in% names(d))
    d <- d %>% dplyr::filter(!as.logical(was_active_before) %in% TRUE)
  if (isTRUE(common_support)) {
    if (!exists("invasion_common_cells", mode = "function"))
      stop("[skill] common_support = TRUE but invasion_common_cells() is not loaded ",
           "(16_invasion_eval.R). Skipping the restriction silently would publish a ",
           "per-fold skill panel whose methods are scored on different origins.",
           call. = FALSE)
    # Resolved per horizon from the SAME frame that is scored, then applied as one filter.
    .keep <- unlist(lapply(sort(unique(d$horizon)), function(h) {
      cl <- invasion_common_cells(lfo_results, h)
      if (!length(cl)) character(0) else paste(h, cl, sep = "\r")
    }), use.names = FALSE)
    if (length(.keep))
      d <- d[paste(d$horizon, d$fold_id, d$health_zone, sep = "\r") %in% .keep, , drop = FALSE]
  }

  d %>%
    dplyr::group_by(method, horizon, fold_id, cutoff) %>%
    dplyr::group_modify(function(g, ...) {
      y <- as.integer(g$is_new_invasion); p <- g$p_invasion
      n <- length(y); n_pos <- sum(y, na.rm = TRUE)   # na.rm: never NA the guard
      ap <- .auc_pr(p, y)
      # TWO tie conventions, because the two things they serve are different:
      #   rk_max ("max") is the OPERATIONAL rank — a zone counts as inside the top K only when
      #     monitoring K zones necessarily includes it. Top-K membership below uses this, as do
      #     .ranking_metrics() and lead_time_analysis().
      #   rk_avg ("average") is the REPORTING rank, and it is what
      #     .ranking_metrics()'s mean_rank_of_truth uses (16_invasion_eval.R). A truth tied at
      #     p = 0 with hundreds of other zones should sit at the middle of its tie group, not
      #     at its far end.
      # These used to be one vector: this per-fold series published a "mean rank of invaded
      # zones" computed under "max" while the pooled table published the same-named quantity
      # under "average", so the two disagreed for exactly the models whose scores tie.
      rk     <- rank(-p, ties.method = "max")         # pessimistic on zero-ties (1 = highest)
      rk_avg <- rank(-p, ties.method = "average")     # reporting rank, matches 16's
      # Top-K membership from the TIE-AWARE rank already computed above, not from
      # head(order(p), k). order() breaks ties by row position, so the previous form made
      # hit/prec/recall@K depend on zone ordering while mean_rank_of_truth beside it was
      # tie-aware — the two disagreed on the same fold. `rk <= k` credits a zone only if
      # monitoring K zones necessarily includes it, matching .ranking_metrics()
      # (16_invasion_eval.R) and lead_time_analysis() below.
      inK  <- rk <= k
      nk   <- sum(inK)
      hits <- sum(y[inK])
      tibble::tibble(
        n_atrisk = n, n_events = n_pos, base_rate = n_pos / max(n, 1),
        auc_pr = ap, auc_roc = .auc_roc(p, y),
        # PREVALENCE-NORMALISED to [0,1], not the raw lift ap / base_rate.
        #
        # The lift's attainable MAXIMUM is n / n_pos, which varies across folds with the number
        # of at-risk zones and, above all, with how many invasions happened to land in the fold
        # (one event among n zones allows a lift of n). Plotted over time it is therefore
        # dominated by the fold's event count, not by how well the model ranked, and it can
        # REVERSE the true ordering: a fold with fewer events can draw a higher lift than a fold
        # with more events and a strictly higher auc_pr. That is the retained
        # skill_over_time_auc_pr_skill figure, whose question is "does discrimination improve as
        # the outbreak matures?".
        #
        # (ap - base) / (1 - base) is the standard skill form: 0 = no better than prevalence,
        # 1 = perfect, comparable across folds, and monotone in auc_pr at a given base rate.
        # NOTE: the POOLED auc_pr_skill in 16_invasion_eval.R is left as a
        # lift, which is defensible there because every method is scored on one shared pool with
        # one base rate; only the per-fold series needed this.
        auc_pr_skill = if (n_pos > 0 && n_pos < n) (ap - n_pos / n) / (1 - n_pos / n) else NA_real_,
        # The raw lift is retained as its own column so nothing that used it is silently changed.
        auc_pr_lift  = if (n_pos > 0) ap / (n_pos / n) else NA_real_,
        hit_at_k  = if (n_pos > 0) as.integer(hits > 0) else NA_integer_,
        # Denominator is the REALISED top-K size, not k: a fold with fewer than k at-risk
        # zones (or a tie group straddling the boundary) has nk < k, and dividing by k
        # would understate precision for a reason that has nothing to do with the model.
        prec_at_k = if (n_pos > 0 && nk > 0) hits / nk else if (n_pos > 0) 0 else NA_real_,
        recall_at_k = if (n_pos > 0) hits / n_pos else NA_real_,
        mean_rank_of_truth = if (n_pos > 0) mean(rk_avg[y == 1]) else NA_real_)
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
#' Using the `horizon` forecast (default h = 1) at every fold, rank the at-risk zones; for a zone whose
#' first case falls in week W, a fold at cutoff C < W gives lead (W-C)/7 weeks.
#' The zone's anticipation is the LARGEST lead at which it was in the top-K — i.e.
#' the earliest the model would have raised the alarm. best_rank is its single
#' best (lowest) rank across those pre-invasion folds.
#'
#' @param first_case_week tibble(health_zone, first_wk) from
#'   first_case_from_zone_week().
#' @return tibble(health_zone, first_wk, n_prefolds, best_rank,
#'   anticipation_weeks, ever_topk): one row per zone of `first_case_week` scored at
#'   `horizon` in at least one fold before its first case (whether or not it was ever in the top K;
#'   ever_topk says which, and anticipation_weeks is 0 for a zone never flagged).
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
#' mean/max predicted probability, best within-fold rank. Flags (k = top-K size):
#'   hit              invaded (an is_new_invasion row at this horizon) and inside the top K
#'                    at one or more of its at-risk folds (best_rank <= k);
#'   miss             invaded but never inside the top K at any at-risk fold;
#'   false_alarm      never invaded, with mean_p in the top decile of zones (>= the 90th
#'                    percentile of mean_p across this method's zones);
#'   correct_negative everything else.
#' Joined to province for spatial mapping.
#'
#' PROBABILITY SCALE. mean_p / max_p (and therefore the false-alarm threshold) are on the
#' PRIMARY, recalibrated scale -- the prequential p_recal of 16b, which is the raw probability
#' on rounds without a fitted delta -- whenever the frame carries a usable p_recal for this
#' method; otherwise on the raw scale, with a warning. Unlike the within-fold ranks, a mean
#' over folds is NOT invariant to the per-fold recalibration factor, so the scale matters;
#' `prob_scale` on every row says which one was used.
#'
#' @return tibble(health_zone, province, n_atrisk_folds, n_invasions, mean_p,
#'   max_p, best_rank, flag, prob_scale) sorted by mean_p.
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
  # The scale of mean_p / max_p (see the docstring): the recalibrated column when every scored
  # row of this method carries it, else the raw one -- said out loud, never silently mixed.
  .rec <- "p_recal" %in% names(d) && all(is.finite(d$p_recal))
  if (!.rec)
    warning(sprintf(paste0("[zone_spatial_error] %s: %s, so mean_p / max_p and the false-alarm ",
                           "threshold are on the RAW probability scale."), method,
                    if ("p_recal" %in% names(d)) "p_recal is missing on some scored rows"
                    else "the frame carries no p_recal"), call. = FALSE)
  d$.p <- if (.rec) d$p_recal else d$p_invasion
  d <- d %>%
    dplyr::group_by(fold_id) %>%
    # Pessimistic ranking for ties: a zone buried in a large group tied at
    # p_invasion == 0 (no mobility pathway) must NOT be credited as "detected in
    # top-K". min_rank gives every tied zone the group minimum (top) rank; use
    # ties.method="max" so tied zones take the WORST position in their run.
    # (Within a fold the recalibration is one monotone transform, so this rank is the same on
    # either scale.)
    dplyr::mutate(rank = rank(-p_invasion, ties.method = "max")) %>%
    dplyr::ungroup()
  z <- d %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(
      n_atrisk_folds = dplyr::n(),
      n_invasions = sum(is_new_invasion, na.rm = TRUE),
      mean_p = mean(.p, na.rm = TRUE),
      max_p  = max(.p, na.rm = TRUE),
      best_rank = min(rank, na.rm = TRUE),
      .groups = "drop") %>%
    # A round with no fitted delta (delta_estimable FALSE) carries its RAW probability in p_recal,
    # so the stamp says how many of the rounds averaged here were left uncorrected.
    dplyr::mutate(prob_scale = if (!.rec) "raw" else {
      .nf <- dplyr::n_distinct(d$fold_id)
      .nu <- if ("delta_estimable" %in% names(d))
               dplyr::n_distinct(d$fold_id[!(d$delta_estimable %in% TRUE)]) else NA_integer_
      if (is.na(.nu) || .nu == 0L) "recalibrated (prequential)"
      else sprintf("recalibrated (prequential); %d of %d rounds uncorrected", .nu, .nf)
    })
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
