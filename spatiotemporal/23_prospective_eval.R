# =============================================================================
# 23_prospective_eval.R — Prospective real-time evaluation & forecast-vs-observed
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Implements the review's evaluation asks (2026-08-06 review):
#   §3.4/§3.8  prospective_invasion_check()      — zones invaded AFTER the last CV
#              cutoff: were they highly ranked BEFORE their first case? (systematic,
#              not anecdotal; subsumes the Rethy/Kisangani/Isiro claim)
#   §3.5       reliability_curve()               — calibration curve (binned predicted
#              probability vs observed invasion frequency, Wilson intervals)
#   §3.5       aggregate_count_calibration()     — per-week predicted expected new
#              invasions (sum of probabilities) vs observed new invasions
#   §3.6       performance_over_time()           — per-fold skill NORMALISED to each
#              fold's own base rate + event counts + a trend test (fold-comparability)
#
# All functions are pure (data.frame in -> data.frame out) and carry no pipeline
# side effects, so they are unit-testable in isolation and reusable for the CV
# folds, the last-two-weeks holdout, and the prospective window alike.
# =============================================================================

#' Wilson score interval for a binomial proportion (k successes of n).
#' Vectorised; returns a data.frame(p, lo, hi). NA-safe for n == 0.
#' @param k successes; n trials; conf confidence level (default 0.95).
wilson_ci <- function(k, n, conf = 0.95) {
  z  <- stats::qnorm(1 - (1 - conf) / 2)
  p  <- ifelse(n > 0, k / n, NA_real_)
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  hw  <- (z * sqrt(pmax(0, p * (1 - p) / n + z^2 / (4 * n^2)))) / den
  data.frame(p = p, lo = ifelse(n > 0, pmax(0, ctr - hw), NA_real_),
             hi = ifelse(n > 0, pmin(1, ctr + hw), NA_real_))
}

#' Reliability / calibration curve for a rare binary forecast (review §3.5, view 1).
#'
#' Bins predicted probabilities and, per bin, reports the mean predicted probability
#' and the observed invasion frequency with a Wilson interval. A perfectly calibrated
#' forecast lies on the 45-degree line; systematic over-prediction (as the manuscript
#' reports, calibration-in-the-large 1.64x/1.94x) shows as points below it.
#'
#' @param pred numeric predicted probabilities in [0,1] (at-risk zone-forecasts).
#' @param obs  0/1 realised invasion outcomes, same length as pred.
#' @param n_bins number of equal-width probability bins (default 10).
#' @param conf Wilson confidence level.
#' @return data.frame(bin, p_mid, mean_pred, n, k, obs_freq, obs_lo, obs_hi).
reliability_curve <- function(pred, obs, n_bins = 10L, conf = 0.95) {
  ok <- is.finite(pred) & !is.na(obs)
  pred <- pred[ok]; obs <- as.integer(obs[ok])
  if (!length(pred)) return(data.frame())
  brks <- seq(0, 1, length.out = n_bins + 1L)
  bin  <- cut(pmin(pmax(pred, 0), 1), breaks = brks, include.lowest = TRUE, labels = FALSE)
  out <- lapply(seq_len(n_bins), function(b) {
    idx <- which(bin == b); n <- length(idx)
    if (!n) return(NULL)
    k <- sum(obs[idx]); ci <- wilson_ci(k, n, conf)
    data.frame(bin = b, p_mid = (brks[b] + brks[b + 1]) / 2,
               mean_pred = mean(pred[idx]), n = n, k = k,
               obs_freq = ci$p, obs_lo = ci$lo, obs_hi = ci$hi)
  })
  do.call(rbind, out)
}

#' Aggregate-count calibration over time (review §3.5, view 2).
#'
#' Per time bucket (week/origin), the predicted EXPECTED number of new invasions is
#' the sum of at-risk zones' invasion probabilities; the observed number is the count
#' of realised first cases. Systematic divergence and its drift over time expose the
#' over-prediction the fixed-beta model incurs as susceptibles depletes (links §3.6/§4.2).
#'
#' @param df data.frame with columns: time (week/origin key), p_invasion, invaded (0/1).
#' @return data.frame(time, expected, observed, n_atrisk, ratio) sorted by time.
aggregate_count_calibration <- function(df) {
  stopifnot(all(c("time", "p_invasion", "invaded") %in% names(df)))
  d <- df[is.finite(df$p_invasion) & !is.na(df$invaded), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  g <- split(d, d$time)
  out <- lapply(names(g), function(tk) {
    x <- g[[tk]]; ex <- sum(x$p_invasion); ob <- sum(x$invaded)
    data.frame(time = tk, expected = ex, observed = ob, n_atrisk = nrow(x),
               ratio = ifelse(ob > 0, ex / ob, NA_real_))
  })
  res <- do.call(rbind, out)
  res[order(res$time), , drop = FALSE]
}

#' Prospective check: were zones invaded AFTER the last CV cutoff ranked highly
#' BEFORE their first case? (review §3.4/§3.8 — replaces the anecdotal claim.)
#'
#' For each zone whose eventually-observed first-onset week falls strictly after the
#' last CV cutoff, look up its predicted-risk rank and probability from the forecast
#' issued before its onset, and summarise: how many such zones, their median
#' pre-invasion rank, and the share that were inside a top-K watch-list.
#'
#' @param forecast_df data.frame with health_zone and a ranking key. Provide either
#'   `rank_med` (posterior-median rank, preferred) or `p_invasion` (a rank is then
#'   derived, 1 = highest). One row per at-risk zone (the forecast issued before the
#'   invasions being checked).
#' @param invaded_after character vector of zones invaded after the last cutoff.
#' @param top_k watch-list size for the "captured" share (default 15).
#' @return list(summary = one-row data.frame, per_zone = data.frame). `summary` has
#'   n_zones, median_rank, mean_rank, n_in_topk, share_in_topk, n_atrisk.
prospective_invasion_check <- function(forecast_df, invaded_after, top_k = 15L) {
  stopifnot("health_zone" %in% names(forecast_df))
  fd <- forecast_df[!is.na(forecast_df$health_zone), , drop = FALSE]
  # Derive a rank (1 = highest risk) from rank_med if present, else from p_invasion.
  if ("rank_med" %in% names(fd) && any(is.finite(fd$rank_med))) {
    fd$.rank <- fd$rank_med
  } else if ("p_invasion" %in% names(fd)) {
    atrisk <- is.finite(fd$p_invasion)
    fd$.rank <- NA_real_
    fd$.rank[atrisk] <- rank(-fd$p_invasion[atrisk], ties.method = "average")
  } else stop("forecast_df needs rank_med or p_invasion")
  n_atrisk <- sum(is.finite(fd$.rank))
  hit <- fd[fd$health_zone %in% invaded_after & is.finite(fd$.rank), , drop = FALSE]
  per_zone <- data.frame(
    health_zone = hit$health_zone,
    pre_invasion_rank = hit$.rank,
    p_invasion = if ("p_invasion" %in% names(hit)) hit$p_invasion else NA_real_,
    in_topk = hit$.rank <= top_k)
  per_zone <- per_zone[order(per_zone$pre_invasion_rank), , drop = FALSE]
  summary <- data.frame(
    n_zones = nrow(per_zone),
    median_rank = if (nrow(per_zone)) stats::median(per_zone$pre_invasion_rank) else NA_real_,
    mean_rank = if (nrow(per_zone)) mean(per_zone$pre_invasion_rank) else NA_real_,
    n_in_topk = sum(per_zone$in_topk),
    share_in_topk = if (nrow(per_zone)) mean(per_zone$in_topk) else NA_real_,
    n_atrisk = n_atrisk, top_k = top_k)
  list(summary = summary, per_zone = per_zone)
}

#' Performance over time with fold-comparability controls (review §3.6).
#'
#' Raw skill is NOT comparable across folds (base rate, at-risk count, event count and
#' task difficulty all drift). This (a) carries each fold's event count and base rate,
#' (b) flags folds below a minimum event count as unreliable, and (c) fits a simple
#' weighted trend of the metric on fold index (weighted by event count) to test whether
#' there is a real slope rather than eyeballing a jagged line. Skill is expressed over
#' each fold's OWN base rate so a decline reflects the model, not the task.
#'
#' @param fold_tbl data.frame with columns: fold (ordered index or origin), metric
#'   (e.g. auc_pr_skill), n_events, base_rate (optional).
#' @param min_events folds with fewer events are flagged unreliable (default 3).
#' @return list(per_fold = annotated data.frame, trend = list(slope, p_value, note)).
performance_over_time <- function(fold_tbl, min_events = 3L) {
  stopifnot(all(c("fold", "metric", "n_events") %in% names(fold_tbl)))
  d <- fold_tbl[order(fold_tbl$fold), , drop = FALSE]
  d$fold_idx  <- seq_len(nrow(d))
  d$reliable  <- is.finite(d$n_events) & d$n_events >= min_events
  fit_d <- d[d$reliable & is.finite(d$metric), , drop = FALSE]
  trend <- list(slope = NA_real_, p_value = NA_real_, note = "insufficient reliable folds")
  if (nrow(fit_d) >= 3L) {
    wt <- pmax(fit_d$n_events, 1)
    fit <- tryCatch(stats::lm(metric ~ fold_idx, data = fit_d, weights = wt),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      sm <- summary(fit)$coefficients
      trend <- list(slope = unname(sm["fold_idx", 1]), p_value = unname(sm["fold_idx", 4]),
                    note = if (sm["fold_idx", 4] < 0.05)
                      (if (sm["fold_idx", 1] < 0) "significant decline over time" else "significant improvement over time")
                    else "no significant trend (flat within noise)")
    }
  }
  list(per_fold = d, trend = trend)
}
