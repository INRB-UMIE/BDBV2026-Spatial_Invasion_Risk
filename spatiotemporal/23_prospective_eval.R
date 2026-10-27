# =============================================================================
# 23_prospective_eval.R — Prospective real-time evaluation & forecast-vs-observed
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Implements the review's evaluation asks (2026-08-06 review):
#   §3.4/§3.8  prospective_invasion_check()      — zones first invaded AFTER every scored
#              back-test outcome window: were they highly ranked by the forecast issued at
#              the last back-test origin, before their first case? (systematic, not
#              anecdotal; subsumes the Rethy/Kisangani/Isiro claim)
#   §3.5       reliability_curve()               — calibration curve (binned predicted
#              probability vs observed invasion frequency, with BOTH a within-bin
#              binomial (Wilson) interval and a cross-validated between-fold one)
#   §3.5       reliability_fold_cells()          — the same bins per LFO fold, i.e. the
#              across-origin spread the CV interval summarises (review 2.25)
#   §3.5       aggregate_count_calibration()     — per-week predicted expected new
#              invasions (sum of probabilities) vs observed new invasions
# (Per-fold skill over time, normalised to each fold's base rate, is spatiotemporal_skill()
# in 19_spacetime_eval.R. performance_over_time(), which had no caller, was removed.)
#
# All functions are pure (data.frame in -> data.frame out) and carry no pipeline
# side effects, so they are unit-testable in isolation and reusable for the CV
# folds, the last-two-weeks holdout, and the prospective window alike.
# =============================================================================

# wilson_ci() is defined in 00_config.R (one definition for the suite).

#' Bin edges for a reliability curve (review §3.5, view 1).
#'
#' WHY THIS IS A CHOICE AND NOT A DETAIL. Invasion probabilities here span five orders of
#' magnitude and are overwhelmingly near zero, while the events are not: equal-WIDTH deciles
#' put nearly every zone-forecast in the first bin and leave the other nine to share the few
#' that remain, often from a single forecast origin. Equal-COUNT deciles are no better for the
#' opposite reason — most of them hold no invasion at all, so most points are "at most this"
#' bounds and the curve rests on the few bins that do.
#'
#' "expected" edges instead split the EXPECTED events (the sum of the predicted probabilities)
#' as evenly as the ties allow: every bin then carries a comparable share of the forecast's
#' event mass, so each has an observed frequency worth plotting and draws on most of the
#' forecast origins, and the bins spread across the probability range the model actually uses.
#' The cost is unequal bin widths, which is why the edges travel with the curve as
#' bin_lo/bin_hi rather than being inferred from a midpoint.
#'
#' THE EDGES MUST NOT DEPEND ON THE OUTCOMES. The previous "event" binning closed each bin at
#' its k-th OBSERVED invasion, i.e. inverse sampling, which biases every bin's observed
#' frequency upward: simulated from perfectly calibrated forecasts with this record's shape, the
#' per-bin observed/predicted ratios came out at 1.05-1.10 instead of 1. Balancing the expected
#' events uses the forecasts alone, so a calibrated forecast plots on the diagonal.
#'
#' Equal predictions are never split across an edge: an edge is pushed to the end of its run of
#' tied probabilities, since two identical forecasts landing in different bins would be an
#' artefact of the sort order.
#'
#' @param pred predicted probabilities; obs 0/1 outcomes (not used by either binning; kept so
#'   callers need not change, and validated for length).
#' @param n_bins target number of bins (default RELIABILITY_N_BINS, 00_config.R). For
#'   "expected" this is a TARGET: ties can yield fewer, and the returned vector says how many
#'   there really are.
#' @param binning "expected" (the default, the suite's convention, RELIABILITY_BINNING) for
#'   edges that balance the expected events; "width" for equal-width bins over [0, 1] (the
#'   conventional diagram).
#' @return numeric vector of bin edges, strictly increasing, spanning [0, 1].
reliability_breaks <- function(pred, obs = NULL,
                               n_bins = get0("RELIABILITY_N_BINS", ifnotfound = 5L),
                               binning = c("expected", "width")) {
  binning <- match.arg(binning)
  n_bins <- max(1L, as.integer(n_bins))
  if (identical(binning, "width"))
    return(seq(0, 1, length.out = n_bins + 1L))
  if (!is.null(obs)) stopifnot(length(obs) == length(pred))
  p <- pmin(pmax(pred[is.finite(pred)], 0), 1)
  tot <- sum(p)
  # Nothing to balance (no forecasts, all-zero forecasts, or one bin asked for): equal-count
  # edges keep the bins populated.
  cuts <- if (!length(p) || n_bins < 2L) numeric(0)
  else if (!(tot > 0))
    stats::quantile(p, probs = seq_len(n_bins - 1L) / n_bins, names = FALSE, type = 1)
  else {
    ps <- sort(p); cum <- cumsum(ps); n <- length(ps)
    vapply(seq_len(n_bins - 1L), function(k) {
      i <- which(cum >= k * tot / n_bins)[1]
      if (is.na(i)) return(NA_real_)
      # Push past the whole run of tied probabilities, then split halfway to the next
      # distinct value so the edge itself is not one of the observed probabilities.
      i <- max(which(ps <= ps[i]))
      if (i >= n) NA_real_ else (ps[i] + ps[i + 1L]) / 2
    }, numeric(1))
  }
  cuts <- sort(unique(cuts[is.finite(cuts) & cuts > 0 & cuts < 1]))
  c(0, cuts, 1)
}

#' Reliability / calibration curve for a rare binary forecast (review §3.5, view 1).
#'
#' Bins predicted probabilities and, per bin, reports the mean predicted probability
#' and the observed invasion frequency with TWO intervals: the within-bin binomial
#' (Wilson) interval, and — when `fold` is supplied — a CROSS-VALIDATED interval that
#' also carries the disagreement between forecast origins.
#' A perfectly calibrated forecast lies on the 45-degree line; systematic
#' over-prediction shows as points below it.
#' (The magnitude is NOT quoted here: calibration-in-the-large is published per model and
#' horizon in invasion_evaluation.csv and moves every run. The figures that show it read it
#' from there; a docstring that hard-codes a pair of values, and attributes them to "the
#' manuscript", goes stale the first time the record grows.)
#'
#' WHY A SECOND INTERVAL (review 2.25). The Wilson interval treats the zone-weeks in a bin
#' as independent draws, so it reports only the sampling error of one binomial proportion.
#' The zone-weeks are not independent: they arrive in LFO folds, and everything a fold
#' shares — the refitted parameters, the epidemic state, that week's reporting — moves the
#' whole fold's outcomes together. What a reader wants from a calibration diagram is the
#' uncertainty of an OUT-OF-SAMPLE forecast, i.e. how far the observed frequency in a bin
#' would move at another forecast origin. That is the between-fold (cross-validated)
#' uncertainty, and it can be larger or smaller than the binomial one; both are reported so
#' the figure can show how much of a bin's departure from the diagonal is within-bin
#' sampling noise and how much is disagreement between rounds.
#'
#' HOW. Clustering on fold, with k_f events of n_f zone-weeks in fold f and the pooled
#' frequency p = sum(k_f)/sum(n_f), the Taylor-linearised cluster-robust variance of that
#' ratio is v_clu = F/((F-1) N^2) * sum((k_f - p n_f)^2) over the F folds that contribute to
#' the bin. Its ratio to the binomial variance v_bin = p(1-p)/N is the design effect
#' (Rao-Scott); the cross-validated interval is then a Wilson interval on the effective
#' sample size N / deff (Korn & Graubard). No bootstrap, no RNG, no tuning constant: the
#' interval is a deterministic function of the fold-level counts.
#'
#' The design effect is NOT ESTIMABLE when a bin holds no events, holds only events, or
#' draws on a single fold — with no between-fold contrast there is nothing to estimate the
#' fold heterogeneity from. Those bins fall back to the binomial interval and say so in
#' `cv_basis`, rather than being dropped or silently reported as if cross-validated.
#'
#' @param pred numeric predicted probabilities in [0,1] (at-risk zone-forecasts).
#' @param obs  0/1 realised invasion outcomes, same length as pred.
#' @param n_bins target number of bins when `breaks` is not given (default
#'   RELIABILITY_N_BINS, 00_config.R); see reliability_breaks().
#' @param conf confidence level of BOTH intervals (default RELIABILITY_CONF, 00_config.R: the
#'   suite's 90%).
#' @param fold optional cluster key, same length as pred: the LFO fold / forecast origin
#'   each zone-forecast was issued at. When absent the CV columns are omitted entirely
#'   (so a caller that never had folds is not handed columns that duplicate the Wilson ones).
#' @param binning,breaks how the bins are formed — see reliability_breaks(). Pass `breaks`
#'   directly to bin a curve and its per-fold cells on IDENTICAL edges; expected-balanced edges
#'   depend on the row set, so recomputing them per caller from a filtered frame would
#'   silently produce two different binnings of the same forecast.
#' @return data.frame(bin, bin_lo, bin_hi, p_mid, mean_pred, n, k, obs_freq, obs_lo, obs_hi)
#'   and, with `fold`, additionally n_folds, n_folds_with_events, deff, n_eff, se_binom,
#'   se_cv, obs_lo_cv, obs_hi_cv, cv_basis.
reliability_curve <- function(pred, obs, n_bins = get0("RELIABILITY_N_BINS", ifnotfound = 5L),
                              conf = get0("RELIABILITY_CONF", ifnotfound = 0.90), fold = NULL,
                              binning = c("expected", "width"), breaks = NULL) {
  ok <- is.finite(pred) & !is.na(obs)
  if (!is.null(fold)) {
    stopifnot(length(fold) == length(pred))
    ok <- ok & !is.na(fold)
  }
  pred <- pred[ok]; obs <- as.integer(obs[ok])
  fold <- if (is.null(fold)) NULL else as.character(fold[ok])
  if (!length(pred)) return(data.frame())
  brks <- if (is.null(breaks)) reliability_breaks(pred, obs, n_bins, match.arg(binning))
          else sort(unique(breaks))
  bin  <- cut(pmin(pmax(pred, 0), 1), breaks = brks, include.lowest = TRUE, labels = FALSE)
  out <- lapply(seq_len(length(brks) - 1L), function(b) {
    idx <- which(bin == b); n <- length(idx)
    if (!n) return(NULL)
    k <- sum(obs[idx]); ci <- wilson_ci(k, n, conf)
    row <- data.frame(bin = b, bin_lo = brks[b], bin_hi = brks[b + 1],
                      p_mid = (brks[b] + brks[b + 1]) / 2,
                      mean_pred = mean(pred[idx]), n = n, k = k,
                      obs_freq = ci$p, obs_lo = ci$lo, obs_hi = ci$hi)
    if (is.null(fold)) return(row)
    cbind(row, .cv_bin_interval(obs[idx], fold[idx], conf))
  })
  do.call(rbind, out)
}

#' Cross-validated (between-fold) interval for ONE reliability bin. Internal helper for
#' reliability_curve(); see its "WHY A SECOND INTERVAL" note for the estimator.
#' @param obs 0/1 outcomes of the zone-weeks in the bin; fold their cluster key.
#' @return one-row data.frame of the CV columns.
.cv_bin_interval <- function(obs, fold, conf = get0("RELIABILITY_CONF", ifnotfound = 0.90)) {
  n_f <- tapply(obs, fold, length); k_f <- tapply(obs, fold, sum)
  F <- length(n_f); N <- sum(n_f); p <- sum(k_f) / N
  v_bin <- p * (1 - p) / N
  v_clu <- if (F >= 2L) F / ((F - 1) * N^2) * sum((k_f - p * n_f)^2) else NA_real_
  # Estimable only with a between-fold contrast to measure: >= 2 folds, a non-degenerate
  # outcome (0 < p < 1, else v_bin is 0 and the ratio is 0/0) and a positive clustered
  # variance (v_clu == 0 means every fold's count sat exactly on p*n_f, which would send
  # n_eff to Inf and print a zero-width interval).
  ok <- F >= 2L && v_bin > 0 && is.finite(v_clu) && v_clu > 0
  deff <- if (ok) v_clu / v_bin else NA_real_
  basis <- if (ok) "between-fold" else if (F < 2L) "binomial (single fold)" else
           if (v_bin <= 0) "binomial (no outcome contrast)" else "binomial (zero fold variance)"
  n_eff <- if (ok) N / deff else N
  ci <- wilson_ci(p * n_eff, n_eff, conf)
  data.frame(n_folds = F, n_folds_with_events = sum(k_f > 0),
             deff = deff, n_eff = n_eff,
             se_binom = sqrt(v_bin), se_cv = if (ok) sqrt(v_clu) else NA_real_,
             obs_lo_cv = ci$lo, obs_hi_cv = ci$hi, cv_basis = basis,
             stringsAsFactors = FALSE)
}

#' Per-fold cells of a reliability curve (review §3.5, view 1, the CV spread itself).
#'
#' The same bins as reliability_curve(), broken down by LFO fold: one row per
#' (bin, fold) with its own trials, events, observed frequency and mean predicted
#' probability. This is what the reliability figure draws behind the pooled points, so a
#' reader sees the between-origin disagreement directly rather than only its summary in
#' the interval; it is also the audit trail for `deff`.
#'
#' @param pred,obs,fold as reliability_curve(); n_bins, binning, breaks likewise. PASS THE
#'   CURVE'S OWN `breaks` whenever the bins are expected-balanced: re-deriving them here from the
#'   same data would happen to agree, but from a filtered or reordered frame it would not, and
#'   the cells would be drawn against bins they do not belong to.
#' @return data.frame(bin, fold, n, k, obs_freq, mean_pred), fold-major within bin.
reliability_fold_cells <- function(pred, obs, fold,
                                   n_bins = get0("RELIABILITY_N_BINS", ifnotfound = 5L),
                                   binning = c("expected", "width"), breaks = NULL) {
  ok <- is.finite(pred) & !is.na(obs) & !is.na(fold)
  pred <- pred[ok]; obs <- as.integer(obs[ok]); fold <- as.character(fold[ok])
  if (!length(pred)) return(data.frame())
  brks <- if (is.null(breaks)) reliability_breaks(pred, obs, n_bins, match.arg(binning))
          else sort(unique(breaks))
  bin  <- cut(pmin(pmax(pred, 0), 1), breaks = brks, include.lowest = TRUE, labels = FALSE)
  key  <- interaction(bin, fold, drop = TRUE, sep = "\r")
  out <- lapply(split(seq_along(pred), key), function(idx) {
    parts <- strsplit(as.character(key[idx[1]]), "\r", fixed = TRUE)[[1]]
    data.frame(bin = as.integer(parts[1]), fold = parts[2], n = length(idx),
               k = sum(obs[idx]), obs_freq = mean(obs[idx]),
               mean_pred = mean(pred[idx]), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  # Order the fold key NATIVELY, for the same reason aggregate_count_calibration() does:
  # the keys are characters (run_all.R passes as.character(fold_id)), so with ten or more
  # folds a plain sort would file "10","11" between "1" and "2" in a table whose whole point
  # is the per-origin breakdown. Numeric-looking keys sort numerically; anything else (ISO week
  # starts, which sort correctly as text) falls back to the character order.
  .num <- suppressWarnings(as.numeric(out$fold))
  fkey <- if (!anyNA(.num)) .num else out$fold
  out <- out[order(out$bin, fkey), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Aggregate-count calibration over time (review §3.5, view 2).
#'
#' Per time bucket (week/origin), the predicted EXPECTED number of new invasions is
#' the sum of at-risk zones' invasion probabilities; the observed number is the count
#' of realised first cases. Systematic divergence and its drift over time expose the
#' over-prediction a fixed-beta hazard incurs over a growing record (links §3.6/§4.2). NOTE:
#' susceptible depletion is NOT modelled anywhere in this pipeline — the drift is a property
#' of holding beta fixed while the at-risk set and the import force change, not of depletion.
#'
#' @param df data.frame with columns: time (week/origin key), p_invasion, invaded (0/1).
#' @return data.frame(time, expected, observed, n_atrisk, ratio) sorted by time.
aggregate_count_calibration <- function(df) {
  stopifnot(all(c("time", "p_invasion", "invaded") %in% names(df)))
  d <- df[is.finite(df$p_invasion) & !is.na(df$invaded), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  # ORDER THE TIME KEY NATIVELY. This used to split() and then order(res$time); split()
  # names are always CHARACTER, so the sort was lexicographic — with >= 10 keys that puts
  # "10","11","12" between "1" and "2". run_all.R passes as.character(fold_id) and the
  # back-test has had ten or more folds, so forecast_count_calibration_h1.csv was written
  # mis-ordered: a table whose ENTIRE purpose is to show expected-vs-observed DRIFT over time,
  # read in the wrong time order. Numeric-looking keys now sort numerically; anything else
  # (e.g. ISO "YYYY-MM-DD" week starts, which sort correctly as text) falls back to sort().
  keys <- unique(as.character(d$time))
  .num <- suppressWarnings(as.numeric(keys))
  keys <- if (!anyNA(.num)) keys[order(.num)] else sort(keys)
  tchr <- as.character(d$time)
  out <- lapply(keys, function(tk) {
    x <- d[tchr == tk, , drop = FALSE]; ex <- sum(x$p_invasion); ob <- sum(x$invaded)
    data.frame(time = tk, expected = ex, observed = ob, n_atrisk = nrow(x),
               ratio = ifelse(ob > 0, ex / ob, NA_real_))
  })
  do.call(rbind, out)
}

#' Prospective check: were zones first invaded AFTER the back-test ranked highly BEFORE
#' their first case? (review §3.4/§3.8 — replaces the anecdotal claim.)
#'
#' OUTSIDE THE CROSS-VALIDATION, or it is not prospective. The caller (run_all.R) passes as
#' `invaded_after` only the zones whose first confirmed case falls in a week starting AFTER
#' last_cut + 7 * max(LFO_HORIZONS), with last_cut the last back-test cutoff: i.e. after the
#' last round's longest outcome window (days (last_cut + 6, last_cut + 6 + 7 * max(LFO_HORIZONS)])
#' had closed, and with no scored invasion event anywhere in the back-test frame. A zone first
#' invaded inside a scored window is a back-test event the featured-model selection already
#' used, so counting it here would re-report the cross-validation.
#'
#' For those zones, look up the predicted-risk rank and probability in `forecast_df` (in the
#' pipeline: the featured model's 1-week forecast issued at the last back-test origin,
#' last_cut + 6, over the zones at risk then), and summarise: how many such zones, their
#' median pre-invasion rank, and the share inside a top-K watch-list. The function itself does
#' not choose the zones or the forecast; it ranks what it is given.
#'
#' @param forecast_df data.frame with health_zone and a ranking key. Provide either
#'   `rank_med` (posterior-median rank, preferred) or `p_invasion` (a rank is then
#'   derived, 1 = highest). One row per at-risk zone (the forecast issued before the
#'   invasions being checked). `p_invasion` is the model's RAW probability (the LFO frame's
#'   column of that name); an optional `p_recal` carries the recalibrated one. Within one
#'   forecast round the recalibration is a single monotone transform, so the rank is the same
#'   on either scale.
#' @param invaded_after character vector of the zones to check (see above).
#' @param top_k watch-list size for the "captured" share (default 15).
#' @return list(summary = one-row data.frame, per_zone = data.frame). `summary` has
#'   n_zones, median_rank, mean_rank, n_in_topk, share_in_topk, n_atrisk, top_k. `per_zone`
#'   has health_zone, pre_invasion_rank, p_invasion_raw, p_recal (only when `forecast_df`
#'   supplies it; the PRIMARY scale) and in_topk.
prospective_invasion_check <- function(forecast_df, invaded_after, top_k = 15L) {
  stopifnot("health_zone" %in% names(forecast_df))
  fd <- forecast_df[!is.na(forecast_df$health_zone), , drop = FALSE]
  # Derive a rank (1 = highest risk) from rank_med if present, else from p_invasion.
  if ("rank_med" %in% names(fd) && any(is.finite(fd$rank_med))) {
    fd$.rank <- fd$rank_med
  } else if ("p_invasion" %in% names(fd)) {
    # ties.method = "max", not "average": `in_topk` below is an OPERATIONAL claim — the
    # zone would have been on a top-K watch-list — and a zone is only genuinely inside a
    # top-K list if monitoring K zones necessarily includes it. With "average" a zone
    # buried in a tie group (e.g. the many zones at p_invasion == 0 with no mobility
    # pathway) could be credited to the watch-list it could not actually have made.
    # Matches lead_time_analysis() (19) and .ranking_metrics() (16); this feeds the
    # headline "n/N zones invaded after the back-test were in the top-15" statement.
    atrisk <- is.finite(fd$p_invasion)
    fd$.rank <- NA_real_
    fd$.rank[atrisk] <- rank(-fd$p_invasion[atrisk], ties.method = "max")
  } else stop("forecast_df needs rank_med or p_invasion")
  n_atrisk <- sum(is.finite(fd$.rank))
  hit <- fd[fd$health_zone %in% invaded_after & is.finite(fd$.rank), , drop = FALSE]
  # Both probability scales, each NAMED: an unlabelled `p_invasion` column read as the primary
  # (recalibrated) scale when it was the raw model output.
  per_zone <- data.frame(
    health_zone = hit$health_zone,
    pre_invasion_rank = hit$.rank,
    p_invasion_raw = if ("p_invasion" %in% names(hit)) hit$p_invasion else rep(NA_real_, nrow(hit)),
    stringsAsFactors = FALSE)
  if ("p_recal" %in% names(hit)) per_zone$p_recal <- hit$p_recal
  per_zone$in_topk <- hit$.rank <= top_k
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
