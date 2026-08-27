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
                            "p_infection_invasion",
                            "q05", "q20", "q25", "q75", "q80", "q95"), names(d))
  # issue_date / lead_days are constant within a (fold_id, cutoff, horizon) target
  # (one issue date per daily fold), so carry them through unchanged like the
  # outcome columns — otherwise the daily-backtest ensemble rows would lose them.
  carry_cols <- intersect(c("is_new_invasion", "was_active_before",
                            "issue_date", "lead_days"), names(d))
  agg <- if (combine == "mean") function(x) mean(x, na.rm = TRUE) else
                                function(x) stats::median(x, na.rm = TRUE)

  # Count only members that actually contribute a USABLE value at each target
  # (non-NA primary probability), not merely a present row — so a target with
  # e.g. one non-NA member and two NA members is NOT emitted as a "2-member"
  # ensemble. `min_members` then genuinely guards against 1-member combinations.
  guard_col <- intersect(c("p_invasion", "p_case_invasion",
                           "p_infection_invasion", "mu_forecast"),
                         val_cols)[1]
  d$.valid <- if (is.na(guard_col)) 1L else as.integer(!is.na(d[[guard_col]]))

  g <- d %>% dplyr::group_by(dplyr::across(dplyr::all_of(key_cols)))
  out <- g %>%
    dplyr::summarise(
      # Count DISTINCT contributing member methods, not valid rows: a member that
      # emits a duplicate (health_zone, horizon) row must not be counted (or later
      # weighted) twice, so `min_members` stays a true member-count guard.
      .n_members = dplyr::n_distinct(method[.valid == 1L]),
      dplyr::across(dplyr::all_of(val_cols), agg),
      dplyr::across(dplyr::all_of(carry_cols), ~ dplyr::first(.x)),
      .groups = "drop") %>%
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
