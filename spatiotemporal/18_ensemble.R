# =============================================================================
# 18_ensemble.R — Invasion Forecast Ensembles (Q3)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Combine several member forecasters into a consensus invasion forecast. For a
# rare binary event with only a handful of training signals, no single mobility/
# GT/observation formulation is reliably best, and ensembles of epidemic models
# systematically match or beat their best member out-of-sample (Reich et al.
# 2019 PNAS 116:3146-3154; Cramer et al. 2022 PNAS 119:e2113561119; Ray et al.
# 2023 Int J Forecast 39:1366-1383). We provide the two canonical combination rules used
# by the US/Euro forecast hubs:
#
#   * mean   — linear opinion pool: p_ens = mean_m p_m  (and Vincent averaging of
#              the count quantiles: q_ens(tau) = mean_m q_m(tau)).
#   * median — median-probability / median-of-quantiles: robust to a single
#              badly-behaved member, the hubs' default operational combiner.
#
# The member set is PRE-SPECIFIED (the mobility-informed renewal family), so the
# ensemble involves no selection on the test folds — it is evaluated in exactly
# the same leakage-free LFO as every other model.
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))
suppressPackageStartupMessages({ library(tidyverse) })

#' Combine member forecasts into one ensemble method.
#'
#' Works on any long forecast tibble whose rows are one (method × forecast
#' target). Forecast-target identity is the set of key columns present among
#' {fold_id, cutoff, training_cutoff, health_zone, horizon}; value columns
#' (probabilities, mu, count quantiles) are combined across members by `combine`;
#' outcome/mask attributes (is_new_invasion, was_active_before) are carried
#' through (they are constant across members for a given target).
#'
#' @param fc_long  long forecast tibble with a `method` column.
#' @param members  character vector of member method names to combine.
#' @param combine  "mean" (linear pool / Vincent) or "median".
#' @param label    ensemble method label (default "Ensemble-<combine>").
#' @param min_members minimum members that must be present at a target for it to
#'   receive an ensemble forecast (default 2 — never emit a 1-member "ensemble").
#' @return tibble of ensemble rows (method = label), or NULL if nothing to combine.
ensemble_forecasts <- function(fc_long, members, combine = c("mean", "median"),
                               label = NULL, min_members = 2L) {
  combine <- match.arg(combine)
  if (is.null(label)) label <- paste0("Ensemble-", combine)
  present <- intersect(members, unique(fc_long$method))
  if (length(present) < min_members) {
    warning(sprintf("[ensemble] only %d of %d members present (need >= %d); skipping '%s'.",
                    length(present), length(members), min_members, label))
    return(NULL)
  }
  d <- fc_long[fc_long$method %in% present, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)

  key_cols   <- intersect(c("fold_id", "cutoff", "training_cutoff",
                            "health_zone", "horizon"), names(d))
  val_cols   <- intersect(c("mu_forecast", "mu_wk0", "p_invasion", "p_case_invasion",
                            "q05", "q20", "q25", "q75", "q80", "q95"), names(d))
  # issue_date / window_days are constant within a (fold_id, cutoff, horizon) target
  # (one issue date per daily fold), so carry them through unchanged like the outcome
  # columns — otherwise the daily-backtest ensemble rows would lose them. "lead_days" is
  # accepted too: it is the former name of window_days (renamed because it carried the window
  # LENGTH, not a lead time), so an ensemble built over an older forecast frame still carries
  # the column instead of silently dropping it.
  # eval_age_days / eval_reliable are properties of the FOLD's outcome window, constant
  # within (fold_id, cutoff, horizon) exactly like is_new_invasion, so they carry through
  # unchanged. Without them the ensemble rows would be the only rows in the frame with no
  # reliability label, and every consumer that restricts to settled rounds (the deployment
  # recalibration, the truncation sensitivity) would silently drop the ensembles entirely.
  carry_cols <- intersect(c("is_new_invasion", "was_active_before",
                            "issue_date", "window_days", "lead_days",
                            "eval_age_days", "eval_reliable"), names(d))
  agg <- if (combine == "mean") function(x) mean(x, na.rm = TRUE) else
                                function(x) stats::median(x, na.rm = TRUE)

  # Count only members that actually contribute a USABLE value at each target
  # (non-NA primary probability), not merely a present row — so a target with
  # e.g. one non-NA member and two NA members is NOT emitted as a "2-member"
  # ensemble. `min_members` then genuinely guards against 1-member combinations.
  guard_col <- intersect(c("p_invasion", "p_case_invasion",
                           "mu_forecast"),
                         val_cols)[1]
  d$.valid <- if (is.na(guard_col)) 1L else as.integer(!is.na(d[[guard_col]]))

  # DROP INVALID ROWS BEFORE AGGREGATING, do not merely count them. `.valid` gated the member
  # COUNT but not which rows entered each column's aggregate, so an ensemble row could mix
  # member sets across columns: with members at p_invasion = 0.1/0.2/NA and
  # mu_forecast = 0.1/0.2/0.3, the emitted row had .n_members = 2, p_invasion = 0.15 (2-member
  # mean) and mu_forecast = 0.2 (3-member mean) — and 1 - exp(-0.2) = 0.181 != 0.15, so the two
  # columns of one published row described different ensembles.
  .n_drop <- sum(d$.valid == 0L)
  if (.n_drop > 0L) {
    message(sprintf("[ensemble] %s: %d member row(s) with no usable %s excluded from every aggregate.",
                    label, .n_drop, guard_col))
    d <- d[d$.valid == 1L, , drop = FALSE]
    if (!nrow(d)) return(NULL)
  }

  # WHY EACH COLUMN IS POOLED INDEPENDENTLY, and why `1 - exp(-mu_forecast) != p_invasion`
  # on an ensemble row. That is NOT an inconsistency to be "corrected": it is the same Jensen
  # gap every MEMBER row already carries. In 21_bayesian_renewal.R,
  #     mu_forecast = colMeans(cum)          (posterior mean of the cumulative hazard)
  #     p_invasion  = colMeans(1 - exp(-cum)) (posterior mean of the probability)
  # and for a right-skewed hazard posterior these are far apart — measured on the shipped LFO,
  # ~5% of rows (mostly h=2) have mu_forecast > 20 alongside p_invasion < 0.05, because a few
  # draws with enormous hazard dominate E[cum] while p is bounded at 1.
  #
  # Linear pooling of EACH summary is exactly right for a mixture of member posteriors: by
  # linearity of expectation the mean of the members' E[cum] IS the mixture's E[cum], and the
  # mean of their E[p] IS the mixture's E[p] (verified numerically to full precision). Deriving
  # one column from the other would replace a correct posterior mean with a Jensen-biased
  # transform of the other one.
  g <- d %>% dplyr::group_by(dplyr::across(dplyr::all_of(key_cols)))
  out <- g %>%
    dplyr::summarise(
      # Count DISTINCT contributing member methods, not valid rows: a member that
      # emits a duplicate (health_zone, horizon) row must not be counted (or later
      # weighted) twice, so `min_members` stays a true member-count guard.
      .n_members = dplyr::n_distinct(method[.valid == 1L]),
      dplyr::across(dplyr::all_of(val_cols), agg),
      # first NON-NA, not first: a member emitting is_new_invasion = NA at a target where others
      # carry the outcome would otherwise make the ensemble row unscorable.
      dplyr::across(dplyr::all_of(carry_cols),
                    ~ { .v <- .x[!is.na(.x)]; if (length(.v)) .v[1] else .x[1] }),
      .groups = "drop") %>%
    # MEMBERSHIP MUST NOT VARY BY CELL. `min_members` lets a cell through on a SUBSET of the
    # members, so an ensemble could be a 3-model average at one (fold x zone) and a 2-model
    # average at the next -- a different estimator per row, scored as one method, and the
    # cells where a member dropped out are exactly the hard ones. Verified harmless on the
    # shipped run (all 10,056 targets carry 3 of 3) precisely because nothing had yet caused a
    # member to drop out; the rule permitted it silently. Report any cell that is short, so a
    # ragged ensemble is visible instead of being averaged away.
    { .short <- dplyr::filter(., .n_members >= min_members &
                                 .n_members < max(.n_members, na.rm = TRUE))
      if (nrow(.short))
        warning(sprintf(paste0("[ensemble] %s: %d of %d target cells combine only %s of %d ",
                               "members, so this ensemble is not one estimator across the ",
                               "scored set. Check why members are missing before publishing."),
                        label, nrow(.short), nrow(.),
                        paste(sort(unique(.short$.n_members)), collapse = "/"),
                        max(.$.n_members, na.rm = TRUE)), call. = FALSE)
      . } %>%
    dplyr::filter(.n_members >= min_members) %>%
    dplyr::select(-.n_members) %>%
    dplyr::mutate(method = label,
                  mobility_id = "ensemble", gt_profile = combine)
  if (nrow(out) == 0) return(NULL)
  out
}

#' Append mean and median ensembles (over `members`) to a long forecast tibble.
#'
#' @param fc_long long forecast tibble.
#' @param members member method names to combine.
#' @param combines which combiners to add (default both).
#' @param prefix ensemble label prefix (default "Ensemble").
#' @return fc_long with the ensemble rows row-bound on (only member/value/key
#'   columns shared with the ensemble are guaranteed; extra columns are filled NA
#'   by bind_rows).
append_ensembles <- function(fc_long, members, combines = c("mean", "median"),
                             prefix = "Ensemble") {
  ens <- lapply(combines, function(cb)
    ensemble_forecasts(fc_long, members, combine = cb,
                       label = paste0(prefix, "-", cb)))
  ens <- Filter(Negate(is.null), ens)
  if (length(ens) == 0) return(fc_long)
  dplyr::bind_rows(fc_long, dplyr::bind_rows(ens))
}

message("[ensemble] 18_ensemble.R loaded — mean / median invasion ensembles (Q3).")
