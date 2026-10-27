# =============================================================================
# 16b_invasion_recalibration.R — post-hoc calibration of invasion probabilities
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# WHY. The invasion suite DISCRIMINATES well but is badly calibrated in the large.
# Measured on the shipped evaluation table (2026-09-19, invasion_evaluation.csv):
# Bayesian AUC-PR runs 0.181-0.298 at h=1 and 0.193-0.379 at h=2, against base rates of
# 0.64% and 1.25%; every Bayesian model over-predicts (calibration_in_large 1.91-3.95 at
# h=1 and 2.27-6.32 at h=2), and the frequentist Distance baseline is ~45-47x. The
# resulting Brier skill is NEGATIVE for 84 of 114 model x horizon pairs (74%) — i.e. the
# probabilities score worse than a constant base-rate forecast even though the
# ranking is highly informative. That is the textbook signature of a well-ordered,
# badly-scaled forecast, and it is fixable with ONE parameter.
#
# WHAT. A single hazard-scale multiplier delta PER MODEL (2026-10-14; it was per model x horizon):
#
#     mu~_i = -log(1 - p_i)                       (cumulative hazard implied by p_i)
#     p_i(delta) = 1 - exp(-delta * mu~_i) = 1 - (1 - p_i)^delta
#
# THE RULE, in full (2026-10-14):
#   * fitted on the model's 1-WEEK forecasts and outcomes only;
#   * only from rounds AFTER the ascertainment burn-in (invasion_is_burnin(), 00_config.R);
#   * only from rounds whose 1-week outcome window had been closed for at least
#     LFO_RELIABLE_EVAL_AGE_DAYS (25) days on the date of the forecast being corrected -- the
#     pipeline's definition of a SETTLED outcome, already used for the deployed factor;
#   * Jeffreys-penalised ML (below), ~ (invasions + 1/2) / sum of forecast hazards;
#   * the SAME delta is applied to the 1-week and 2-week cumulative hazards, so
#     P(within 2 weeks) >= P(within 1 week) holds in every posterior draw. Separate per-horizon
#     factors broke that: 7,382 recalibrated back-test rows (22 models, post-burn-in) had
#     P(2 wk) < P(1 wk) on the 2026-10-09 frame;
#   * in the live forecast delta is drawn per posterior draw from Gamma(invasions + 1/2,
#     rate = (invasions + 1/2) / delta), whose mean is the fitted delta, so the published
#     intervals carry its uncertainty (21_bayesian_renewal.R, .scale_invasion_draws()).
# The back-test ("prequential") delta replays exactly this rule at each historical forecast
# date, so it never uses an outcome that was not yet settled at that date.
#
# Three properties make this the right family here, and each is relied on downstream:
#   (1) IDENTITY AT delta = 1. p_i(1) == p_i exactly, so "recalibrated" and "raw"
#       differ ONLY by delta — no functional-change confound sneaks into the
#       comparison. (Recalibrating `mu_forecast`, the posterior-mean hazard, would
#       NOT have this property: p_invasion is a posterior MEAN probability and
#       1 - exp(-mean(cum)) != mean(1 - exp(-cum)); the two differ by up to 0.045
#       in this run's LFO frame.)
#   (2) STRICT MONOTONICITY in p_i for delta > 0. For a FIXED delta every ranking
#       product — AUC-PR, AUC-ROC, precision/recall/hit@K, mean rank-of-truth, top-K
#       watchlists, lead time, detection curves — is EXACTLY invariant. Recalibration
#       cannot trade discrimination for calibration; it only rescales.
#       IMPORTANT CAVEAT, and it is easy to get wrong: the PREQUENTIAL delta varies
#       BY FOLD, so the guarantee holds WITHIN a fold (i.e. within one issued
#       forecast), not across a pool of folds recalibrated by different factors.
#       Metrics computed per fold and then averaged — every ranking metric in
#       .ranking_metrics(), spatiotemporal_skill(), lead_time_analysis() and
#       compute_detection_curve() — are therefore exactly invariant. Metrics
#       computed on POOLED rows across folds — AUC-PR and AUC-ROC in
#       evaluate_invasion() — are NOT, because pooling compares probabilities issued
#       under different deltas. Selection is unaffected: its AUC-PR axis reads the
#       RAW column. (The 16c comparison audit was removed in the 2026-09-17 streamlining.)
#   (3) HAZARD SCALE, not probability scale. The additive cumulative-hazard identity
#       p = 1 - exp(-sum mu) is preserved, so the transform composes correctly with
#       the multi-week accumulation the cascade layer performs. Isotonic or Platt
#       (logit-scale) recalibration would break it.
#
# ESTIMATOR: JEFFREYS-PENALISED MAXIMUM LIKELIHOOD, WITH NO TUNING CONSTANTS.
#   l(delta)  = sum_{y=1} log(1 - (1-p)^delta) + sum_{y=0} delta * log(1 - p)
#   l*(delta) = l(delta) + 0.5 * log I(delta)                    (Firth 1993; Jeffreys prior)
# l is strictly concave in delta, but its maximum is INTERIOR only when the pool contains at
# least one event; with none it runs to delta = 0 and would zero out a whole fold's forecast
# (measured: delta = 1e-6, maximum recalibrated probability 2e-8). The penalty is -Inf at both
# ends, so l* always has a finite interior maximum. It is tuning-free -- no strength to pick --
# and in the rare-event limit reduces to delta = (events + 1/2) / sum_i -log(1 - p_i), the
# classical "add half an event" rate ratio. See fit_invasion_delta() for the derivation, the
# numerical checks, and the closed-form verification.
#
# THE FLOORS ARE GONE (2026-09-22). This module used to require >= 10 training events AND >= 3
# training origins before fitting delta at all, and to DISCARD any fit that landed on the
# search band, falling back to delta = 1 in each case. That left over half the folds
# uncorrected -- 6 of 11 at h=1 and 5 of 10 at h=2 -- and the fallback, not the estimator, was
# the residual miscalibration: at h=2 the folds carrying a fitted delta sat at obs/exp 1.04
# while those left at delta = 1 sat at 0.70. "Assume perfectly calibrated" is the worst
# available default for a model we have measured to over-predict.
#
# Removing them is safe because the choice was never doing work. Swept prequentially on the
# shared support over the 20 Bayesian models, min_folds in {1,2,3} x min_events in {0,1,5,10}
# moves the mean LOG SCORE by 0.9% at h=1 and 0.6% at h=2, against a 3.4%/3.3% gain from
# recalibrating at all. Per-fold |log(obs/exp)| cannot separate them either: at ~3.3 events per
# fold its value under PERFECT calibration is 0.451 (90% range 0.279-0.635) and every setting
# lands inside that. A rolling training window was also tested and is no better than the
# expanding one (h=2: 0.04255 expanding vs 0.04260-0.04272 for windows of 2-4 folds).
#
# This is deliberately NOT the moment estimator 1/calibration_in_large used by
# cascade_fit_delta(): on this run's LFO frame the likelihood value is systematically smaller
# and scores better on BOTH log score and Brier skill, because a rare-event mean is dominated
# by the many near-zero cells.
#
# TWO DELTAS, DELIBERATELY. Fitting delta on the same rows it is scored on is
# optimistically biased, and the bias is LARGEST for the worst-calibrated models —
# so in-sample recalibration would systematically reward models that needed the most
# fixing. Hence:
#   * delta_preq  — PREQUENTIAL, fitted only on post-burn-in folds whose 1-week outcome was
#                   SETTLED (>= LFO_RELIABLE_EVAL_AGE_DAYS) by the target fold's issue date,
#                   on the shared support. Used for every reported score and for model
#                   selection.
#   * delta_pooled— fitted on every settled, post-burn-in fold. Used ONLY to recalibrate the
#                   live forecast, where those folds are genuinely the whole usable past.
#
# CAUSALITY OF THE PREQUENTIAL SPLIT (2026-10-14). A fold with cutoff C' has a 1-week outcome
# window whose last day is C' + 13. A forecast at cutoff C is issued on C + 6. A training fold is
# usable only if (C + 6) - (C' + 14) >= LFO_RELIABLE_EVAL_AGE_DAYS, i.e. its outcome had been
# SETTLED for 25 days by the issue date; with weekly folds that is C' <= C - 35. The previous
# rule (C' + 13 <= C) admitted folds whose outcome week had ended only 7 days before the issue
# date and scored them with their FINAL outcomes -- invasions that had mostly not yet been
# confirmed at the time, i.e. look-ahead. The cost of the stricter rule was measured on the
# 2026-10-09 frame (22 Bayesian models, post-burn-in): mean log score x1000 22.31 (h=1) / 38.14
# (h=2) under the 25-day rule against 21.98 / 37.28 under a 13-day rule, with both better than
# no recalibration in 43 and 42 of 44 model x horizon rows. The 25-day rule is used because it
# is the same settledness the deployed factor has always required: the prequential delta then
# replays deployment exactly.
#
# CLUSTERING. The point estimate treats zone-weeks as independent (a quasi-likelihood
# M-estimator; the estimating equation is unbiased regardless of correlation), but the
# interval must not: CIs come from a bootstrap that resamples ZONES, the correlated
# unit, matching evaluate_invasion()'s convention.
#
# SCOPE. This module adds columns and files and changes no EXISTING metric, but it is not
# cosmetic: INVASION_SELECT_ON_RECAL and INVASION_RECALIBRATE_DEPLOY are BOTH TRUE by default
# (00_config.R), so the recalibrated log score IS the model-selection axis and the live
# forecast IS recalibrated before publication. This block previously said both flags were
# "FALSE by default", from which a methods reader would conclude that neither the featured
# pick nor the published probabilities are affected by anything in this file. Both are.
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))
suppressPackageStartupMessages({ library(dplyr) })


# ---------------------------------------------------------------------------
# Tunables (overridable from 00_config.R / the environment)
# ---------------------------------------------------------------------------
# Standalone-load fallbacks, mirroring 00_config.R. THERE IS NO EVENT OR ORIGIN FLOOR: the
# delta fit is Jeffreys-penalised and therefore exists for every training pool, so the only
# bar is that at least `min_folds` (= 1) training fold's outcome window has CLOSED -- an
# information constraint, not a threshold. See the 00_config.R block for the sweep showing
# that the floors this module used to carry changed the log score by under 1%.
RECAL_MIN_FOLDS  <- get0("RECAL_MIN_FOLDS",  ifnotfound = 1L)
# Numerical search bracket, NOT a prior. The penalised objective is bounded and unimodal, so
# the optimum is interior and this never binds; a hit means pathological input and is reported.
# It spans both sides of 1 because delta > 1 means the model UNDER-predicts, which this
# outbreak's first LFO fold genuinely did (delta ~1.4 during explosive early growth).
RECAL_BAND       <- get0("RECAL_BAND",       ifnotfound = c(1e-4, 1e4))
RECAL_N_BOOT     <- get0("RECAL_N_BOOT",     ifnotfound = 400L)

# ---------------------------------------------------------------------------
# Core transform
# ---------------------------------------------------------------------------

#' Apply a hazard-scale recalibration factor to invasion probabilities.
#'
#' p -> 1 - (1 - p)^delta. Strictly increasing in p for delta > 0, and the identity
#' at delta = 1. `delta` may be a scalar (one factor for every element) or the same
#' length as `p` (one factor per row — the prequential case, where delta varies by
#' fold). NA in -> NA out. Non-finite or non-positive delta entries leave their
#' probabilities UNCHANGED and raise a warning, so a failed fit degrades to current
#' behaviour rather than silently zeroing a forecast.
#' @param p numeric probabilities in [0, 1].
#' @param delta positive numeric, length 1 or length(p).
#' @return numeric probabilities in [0, 1], same length as p.
recalibrate_invasion_p <- function(p, delta) {
  p <- as.numeric(p); delta <- as.numeric(delta)
  if (length(delta) == 1L) delta <- rep(delta, length(p))
  if (length(delta) != length(p))
    stop("[recal] `delta` must be length 1 or length(p).", call. = FALSE)
  bad <- !is.finite(delta) | delta <= 0
  if (any(bad, na.rm = TRUE)) {
    warning(sprintf("[recal] %d non-usable delta value(s); those probabilities are left unchanged.",
                    sum(bad, na.rm = TRUE)), call. = FALSE)
    delta[bad] <- 1
  }
  # The 1 - 1e-12 clamp keeps log1p(-p) finite at p = 1, but it must not cost us the
  # exact identity at the endpoints: (1 - 1)^delta = 0 for every delta > 0, so p = 1
  # maps to 1 exactly, and p = 0 maps to 0 exactly. Set both explicitly rather than
  # letting the clamp leak a 1e-12 error into what is meant to be an exact transform.
  lq  <- delta * log1p(-pmin(pmax(p, 0), 1 - 1e-12))   # = delta * log(1 - p)
  out <- -expm1(lq)                                    # = 1 - (1 - p)^delta
  out <- pmin(pmax(out, 0), 1)
  out[!is.na(p) & p >= 1] <- 1
  out[!is.na(p) & p <= 0] <- 0
  out[is.na(p)] <- NA_real_
  out
}

#' Jeffreys-penalised maximum-likelihood hazard-scale recalibration factor.
#'
#' theta = log(delta); u_i = -log(1 - p_i) is the raw cumulative hazard; lambda_i = delta*u_i.
#'
#'   l(theta)  = sum_{y=1} log(1 - e^-lambda_i)  -  sum_{y=0} lambda_i
#'   I(theta)  = sum_i lambda_i^2 e^-lambda_i / (1 - e^-lambda_i)       (expected information)
#'   l*(theta) = l(theta) + 0.5 * log I(theta)                          (Firth 1993; Jeffreys)
#'
#' WHY PENALISED, AND WHY THIS PENALTY. The plain MLE does not exist in the interior when the
#' training pool contains no invasions: l is then strictly decreasing and delta runs to 0, so
#' the next fold's probabilities collapse to ~0. Measured on a 3,000-row zero-event pool the
#' unpenalised fit returns delta = 1e-6 and a maximum recalibrated probability of 2e-8 -- a
#' whole fold of forecasts silently zeroed. The previous code avoided this by DISCARDING any
#' estimate that landed on the search band and reverting to delta = 1, which threw away the
#' estimate rather than fixing it.
#'
#' Firth's penalty is the Jeffreys prior's log-density, so it is tuning-free: no strength to
#' choose. It is -Inf at both delta -> 0 and delta -> Inf, so the penalised optimum is always
#' interior and finite. In the rare-event limit that holds here it has a closed form,
#'
#'     delta_hat = (events + 1/2) / sum_i u_i,
#'
#' i.e. exactly "add half an invasion" on the hazard scale -- the classical Jeffreys-corrected
#' rate ratio. (The plain MLE is events / sum_i u_i, which is 0 when events = 0.) In that limit
#' the correction multiplies delta by (events + 1/2) / events, i.e. a relative change of
#' 1 / (2 * events), so it is negligible once the pool holds a few dozen invasions. One
#' estimator therefore serves both the prequential and the pooled fits, with no gates and no
#' special cases.
#'
#' UNIMODALITY. l(theta) is concave: the binomial cloglog log-likelihood is concave in the
#' linear predictor, and theta enters every linear predictor additively (log lambda_i =
#' log u_i + theta). The penalty 0.5 * log I(theta) tends to 0.5 * theta + const in the
#' rare-event limit, so l* is concave there too, and optimize() on the bracket is safe.
#'
#' @param p numeric predicted probabilities (at-risk rows only).
#' @param y 0/1 realised outcomes, same length.
#' @param band length-2 positive NUMERICAL BRACKET for delta (not a prior).
#' @param penalty TRUE for the Jeffreys-penalised fit (the default and what the pipeline uses);
#'   FALSE gives the plain MLE, retained only so the sensitivity can be reproduced.
#' @param warn_boundary emit a warning when the optimum lands on the bracket edge. With the
#'   penalty on, that cannot happen for well-formed data and means something is badly wrong.
#' @param details TRUE returns list(delta, at_band, converged) instead of the bare delta.
#' @return delta (numeric scalar), or NA_real_ if there is nothing to fit at all.
fit_invasion_delta <- function(p, y, band = RECAL_BAND, penalty = TRUE,
                               warn_boundary = TRUE, details = FALSE) {
  stopifnot(length(p) == length(y), length(band) == 2L, band[1] > 0, band[2] > band[1])
  ok <- is.finite(p) & !is.na(y)
  p  <- as.numeric(p)[ok]; y <- as.integer(y)[ok]
  # EVERY early return must honour `details`. They used to return a bare NA_real_ regardless,
  # so prequential_invasion_delta() -- which always calls with details = TRUE and then reads
  # $delta -- raised "$ operator is invalid for atomic vectors" on any slice with no positive
  # hazard. That error escapes attach_invasion_recalibration(), is swallowed by the tryCatch in
  # run_all.R, and the run continues WITHOUT log_score_recal -- at which point .invasion_ls_col()
  # reverts the selection composite to the RAW log score and the featured model (and with it
  # CASCADE_KERNEL) can change, with only a warning in the log.
  .bail <- function() if (isTRUE(details))
    list(delta = NA_real_, at_band = NA_character_, converged = FALSE) else NA_real_
  # The ONLY unfittable case left: no row carries positive hazard, so delta is unidentified at
  # every value. An empty pool and a zero-EVENT pool are both fittable under the penalty.
  u <- -log1p(-pmin(pmax(p, 0), 1 - 1e-12))          # cumulative hazard, >= 0
  if (!length(u) || !any(u > 0)) return(.bail())
  nll <- function(log_d) {
    lam <- exp(log_d) * u
    l1  <- log(-expm1(-lam))                          # log(1 - e^-lam), stable for small lam
    l1[!is.finite(l1)] <- log(.Machine$double.xmin)
    ll  <- sum(l1[y == 1L]) - sum(lam[y == 0L])
    if (isTRUE(penalty)) {
      # lam^2 e^-lam / (1 - e^-lam), which -> lam as lam -> 0. The branch avoids 0/0 there.
      g <- ifelse(lam < 1e-8, lam, lam^2 * exp(-lam) / -expm1(-lam))
      ll <- ll + 0.5 * log(max(sum(g), .Machine$double.xmin))
    }
    -ll
  }
  opt <- tryCatch(stats::optimize(nll, log(band), tol = 1e-10), error = function(e) NULL)
  if (is.null(opt)) return(.bail())
  d <- exp(opt$minimum)
  rel <- (log(d) - log(band[1])) / (log(band[2]) - log(band[1]))
  at_band <- if (rel < 1e-4) "lower" else if (rel > 1 - 1e-4) "upper" else NA_character_
  # REPORTED, NEVER DISCARDED. The old code treated a bracket hit as "not estimable" and fell
  # back to delta = 1 -- i.e. it answered "the data say the forecast is wildly mis-scaled" with
  # "assume it is perfectly scaled". Under the penalty a bracket hit is not reachable for
  # well-formed data, so if it fires the input is pathological and must be looked at.
  if (warn_boundary && !is.na(at_band))
    warning(sprintf(paste0("[recal] delta reached the %s edge of the numerical bracket ",
                           "(%.4g in [%.4g, %.4g]). Under the Jeffreys penalty the optimum ",
                           "should be interior, so this indicates pathological input -- ",
                           "inspect before using this fit."), at_band, d, band[1], band[2]),
            call. = FALSE)
  if (isTRUE(details)) list(delta = d, at_band = at_band, converged = is.na(at_band)) else d
}

#' Zone-clustered bootstrap interval for a recalibration factor.
#'
#' Resamples ZONES with replacement (the correlated unit; the same convention
#' evaluate_invasion() uses for its AUC-PR / log-score intervals) and refits delta on
#' each replicate. Replicates with too few events to identify delta are dropped and
#' COUNTED, so a thin interval cannot masquerade as a well-determined one.
#' @return list(lo, hi, n_boot_used, se_log) — 90% interval, the pipeline's band convention.
bootstrap_invasion_delta <- function(p, y, zone, band = RECAL_BAND,
                                     n_boot = RECAL_N_BOOT,
                                     seed = get0("RANDOM_SEED", ifnotfound = 20260704L)) {
  ok <- is.finite(p) & !is.na(y) & !is.na(zone)
  p <- as.numeric(p)[ok]; y <- as.integer(y)[ok]; zone <- as.character(zone)[ok]
  if (!length(p)) return(list(lo = NA_real_, hi = NA_real_, n_boot_used = 0L, se_log = NA_real_))
  zs  <- unique(zone)
  idx <- split(seq_along(zone), factor(zone, levels = zs))
  # Restore the caller's RNG stream exactly. If there was NO stream before, remove the
  # one we created rather than leaving a seeded generator behind — otherwise calling
  # this function would silently change every subsequent random draw in the session.
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set_fixed_seed(seed)
  on.exit({
    if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  bs <- vapply(seq_len(n_boot), function(b) {
    i <- unlist(idx[sample.int(length(zs), length(zs), replace = TRUE)], use.names = FALSE)
    as.numeric(fit_invasion_delta(p[i], y[i], band = band,
                                  warn_boundary = FALSE))
  }, numeric(1))
  bs <- bs[is.finite(bs)]
  if (length(bs) < 10L)
    return(list(lo = NA_real_, hi = NA_real_, n_boot_used = length(bs), se_log = NA_real_,
                n_at_band_lo = NA_integer_, n_at_band_hi = NA_integer_))
  # CENSOR THE ENDPOINTS when too many replicates are pinned on the search band. Replicates that
  # did not converge sit exactly on band[1] or band[2], and a raw quantile then REPORTS THE BAND
  # EDGE AS A CONFIDENCE LIMIT. On the current run Distance-B1 h=2 had 34/400 replicates (8.5%)
  # at the floor, so the published delta_lo was 0.01000000 — the floor itself, not a quantile;
  # Adjacency-B7 h=2 had 114/400 (28.5%) at the ceiling and published delta_hi = 4.9999998.
  # When the tail mass at a band exceeds the quantile level, that endpoint is not identified.
  .tol <- 1e-6
  n_lo <- sum(bs <= band[1] * (1 + .tol)); n_hi <- sum(bs >= band[2] * (1 - .tol))
  q <- stats::quantile(bs, c(0.05, 0.95), names = FALSE)
  lo <- if (n_lo / length(bs) > 0.05) NA_real_ else q[1]
  hi <- if (n_hi / length(bs) > 0.05) NA_real_ else q[2]
  if (is.na(lo) || is.na(hi))
    warning(sprintf("[recal] bootstrap interval censored: %d/%d replicates at the lower band, %d/%d at the upper; the affected endpoint is not identified.",
                    n_lo, length(bs), n_hi, length(bs)), call. = FALSE)
  # se_log over the UNPINNED replicates only — the sd of a clamped set understates the spread.
  .free <- bs[bs > band[1] * (1 + .tol) & bs < band[2] * (1 - .tol)]
  list(lo = lo, hi = hi, n_boot_used = length(bs),
       se_log = if (length(.free) >= 10L) stats::sd(log(.free)) else NA_real_,
       n_at_band_lo = n_lo, n_at_band_hi = n_hi)
}

# ---------------------------------------------------------------------------
# Fold bookkeeping
# ---------------------------------------------------------------------------

#' At-risk, finitely-predicted rows — byte-for-byte the SAME filter
#' evaluate_invasion() applies, so every row it scores carries a p_recal.
#'
#' Deliberately does NOT also drop rows with a missing outcome. evaluate_invasion()
#' keeps those rows, so dropping them here would leave p_recal = NA on a row that IS
#' scored, which trips evaluate_invasion()'s support check and silently disables the
#' recalibrated metrics for the whole slice. Missing outcomes are excluded from the
#' FIT instead, inside fit_invasion_delta() (`ok <- is.finite(p) & !is.na(y)`), which
#' is where the exclusion belongs.
.recal_scorable <- function(d) {
  if ("was_active_before" %in% names(d))
    d <- d[!(as.logical(d$was_active_before) %in% TRUE), , drop = FALSE]
  # RANK-ONLY METHODS ARE NOT RECALIBRABLE. Recalibration maps a PROBABILITY through
  # 1-(1-p)^delta; a comparator that emits an ordering (prob_calibrated = FALSE — Distance-B1
  # has no fitted scale at all, Adjacency-B7 is 1/(1 + travel-time minutes)) has no probability
  # to move, and any delta fitted to it is a property of the arbitrary monotone transform, not
  # of the forecast. evaluate_invasion() already NAs the recalibrated block for these, but this
  # function feeds invasion_recalibration.csv/.json too, which published deltas for them
  # (Distance-B1 0.016, Adjacency-B7 2.48/4.88 on the shipped frame) as if they meant something.
  # Absent column = calibrated, which is right for every Bayesian model.
  if ("prob_calibrated" %in% names(d))
    d <- d[!(d$prob_calibrated %in% FALSE), , drop = FALSE]
  d[is.finite(d$p_invasion), , drop = FALSE]
}

#' Map fold_id -> cutoff date for one (method, horizon) slice.
#' Falls back to a synthetic weekly calendar keyed on fold_id when `cutoff` is absent,
#' which preserves the ORDER and the 7-day spacing the causality rule assumes.
.recal_fold_cutoffs <- function(d) {
  fid <- sort(unique(d$fold_id))
  if ("cutoff" %in% names(d)) {
    cut <- as.Date(vapply(fid, function(k) {
      v <- unique(as.Date(d$cutoff[d$fold_id == k])); v <- v[!is.na(v)]
      if (length(v)) as.character(min(v)) else NA_character_
    }, character(1)))
    if (!anyNA(cut)) return(stats::setNames(cut, as.character(fid)))
    warning("[recal] some folds carry no usable cutoff date; falling back to fold_id spacing.",
            call. = FALSE)
  }
  stats::setNames(as.Date("2020-01-06") + 7L * (seq_along(fid) - 1L), as.character(fid))
}

#' Prequential recalibration factors: the deployment rule replayed at each forecast date.
#'
#' For each fold k with cutoff C (issued on C + 6), delta is fitted on the `fit` rows of every
#' fold j that is (a) NOT in the ascertainment burn-in and (b) SETTLED at the issue date:
#' (C + 6) - (C_j + 7h + 7) >= `settled_days`, with h the horizon of the TRAINING rows (1 in
#' the pipeline). Folds with nothing usable get delta = 1 (no recalibration), which mirrors what
#' a live system holds before delta is estimable. The same per-fold delta is then applied to
#' every horizon of that fold by the caller.
#'
#' @param d the rows that define WHICH folds need a delta and their cutoffs (any horizon).
#' @param horizon horizon of the TRAINING rows in `fit` (drives the settledness lag).
#' @param fit the training rows (default `d`); the pipeline passes the 1-week rows.
#' @param settled_days outcome age (days after the window closes) at which it counts as known.
#' @param burnin_fn function(cutoff) -> TRUE for burn-in rounds (excluded from training).
#'   Ignored when `d` has no `cutoff` column (the fold_id spacing fallback has no real dates).
#' @return data.frame(fold_id, delta_preq, delta_estimable, delta_at_band, n_train_events,
#'   n_train_folds, n_train_hazard)
prequential_invasion_delta <- function(d, horizon = 1L, band = RECAL_BAND,
                                       min_folds = RECAL_MIN_FOLDS, fit = NULL,
                                       settled_days = get0("LFO_RELIABLE_EVAL_AGE_DAYS",
                                                           ifnotfound = 25L),
                                       burnin_fn = get0("invasion_is_burnin", mode = "function",
                                                        ifnotfound = NULL)) {
  # `d` defines WHICH FOLDS need a delta and their cutoffs; `fit` (default: `d`) supplies the
  # rows the delta is TRAINED on. They differ when the caller restricts training to the 1-week
  # rows on the shared (fold x zone) support: every fold must still come back with a delta --
  # including one outside that support -- while no method trains on cells its competitors were
  # denied.
  if (is.null(fit)) fit <- d
  cuts <- .recal_fold_cutoffs(d)
  fid  <- as.integer(names(cuts))
  # Burn-in needs REAL dates: the fold_id-spacing fallback calendar is synthetic.
  .real_dates <- "cutoff" %in% names(d) && !anyNA(suppressWarnings(as.Date(d$cutoff)))
  in_burnin <- if (.real_dates && is.function(burnin_fn)) burnin_fn(cuts) else rep(FALSE, length(cuts))
  out <- lapply(seq_along(fid), function(i) {
    k <- fid[i]; issue <- cuts[i] + 6L
    # SETTLED, not merely closed: the outcome window ended at least `settled_days` before the
    # issue date, so its invasions had been confirmed when this forecast was made.
    settled <- as.numeric(issue - (cuts + 7L * as.integer(horizon) + 7L)) >= settled_days
    usable  <- fid[settled & !in_burnin]
    tr <- fit[fit$fold_id %in% usable, , drop = FALSE]
    n_ev <- sum(tr$is_new_invasion == 1L, na.rm = TRUE)
    # FITTED AS SOON AS THERE IS ANYTHING TO FIT. The penalised estimate exists for any pool
    # with positive hazard, including one with no invasions, so the only bar left is that at
    # least `min_folds` (= 1) training fold has CLOSED. The old event/origin floors left over
    # half the folds at delta = 1 -- 6 of 11 at h=1 and 5 of 10 at h=2 -- and that fallback,
    # not the estimator, was the residual miscalibration: at h=2 the folds carrying a fitted
    # delta sat at obs/exp 1.04 while those left at delta = 1 sat at 0.70.
    .dfit <- if (length(usable) >= min_folds && nrow(tr))
               fit_invasion_delta(tr$p_invasion, tr$is_new_invasion, band = band,
                                  warn_boundary = FALSE, details = TRUE)
             else list(delta = NA_real_, at_band = NA_character_, converged = FALSE)
    dl <- .dfit$delta
    # A bracket hit is REPORTED, not discarded (see fit_invasion_delta): reverting to delta = 1
    # answers "the data say this forecast is wildly mis-scaled" with "assume it is perfectly
    # scaled". Under the penalty it is unreachable for well-formed data.
    if (!is.na(.dfit$at_band))
      warning(sprintf(paste0("[recal] fold %s: prequential delta reached the %s edge of the ",
                             "numerical bracket (%.4g); the fit is USED but the input is ",
                             "pathological -- inspect it."),
                      as.character(k), .dfit$at_band, as.numeric(dl)), call. = FALSE)
    .usable_dl <- is.finite(dl) && dl > 0
    data.frame(fold_id = k,
               delta_preq = if (.usable_dl) as.numeric(dl) else 1,
               delta_estimable = .usable_dl,
               delta_at_band = .dfit$at_band,
               n_train_events = n_ev,
               n_train_folds = length(usable),
               n_train_hazard = if (nrow(tr)) sum(-log1p(-pmin(pmax(tr$p_invasion, 0), 1 - 1e-12)))
                                else 0)
  })
  dplyr::bind_rows(out)
}

# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------

#' Give rank-only SCORED rows the identity recalibration (delta = 1, p_recal = p_invasion).
#'
#' A comparator that emits an ordering rather than a probability (prob_calibrated = FALSE) is
#' excluded from the delta FIT by .recal_scorable(), and that is right -- a delta fitted to an
#' arbitrary monotone scale means nothing. But leaving its p_recal at NA is not: recalibration
#' is a strictly monotone transform, so for a rank-only comparator the recalibrated and raw
#' scales are the SAME ORDERING, and delta = 1 states exactly that.
#' append_naive_detection_curve_model() (20_forecast_detail.R) already used this convention for
#' the injected baselines; attach did not, and the mismatch cost the figures dearly: 22,112
#' scored rows with NA p_recal made fs_lfo_col() (forecast_scale.R) fall back to the RAW column
#' for the WHOLE frame, so every panel in both suites was drawn raw while captioned
#' "recalibrated". evaluate_invasion() still reports their calibration block as NA -- it keys on
#' prob_calibrated, not on p_recal -- so nothing gains a proper score it should not have.
#'
#' Rows that are NOT scored (already-affected zones, non-finite p) keep NA: nothing is
#' fabricated for a row no metric reads.
.recal_identity_fill <- function(d) {
  if (!"prob_calibrated" %in% names(d)) return(d)
  ro  <- d$prob_calibrated %in% FALSE
  act <- if ("was_active_before" %in% names(d)) as.logical(d$was_active_before) %in% TRUE
         else rep(FALSE, nrow(d))
  fill <- ro & is.finite(d$p_invasion) & !act & is.na(d$p_recal)
  if (any(fill)) {
    d$delta_preq[fill]      <- 1
    d$delta_estimable[fill] <- FALSE
    d$p_recal[fill]         <- d$p_invasion[fill]
  }
  d
}

#' Attach prequential recalibration to an LFO result frame.
#'
#' Adds, per (method, horizon, fold):
#'   delta_preq       the causally-fitted recalibration factor (1 where not estimable)
#'   delta_estimable  whether a delta was actually fitted for that fold
#'   n_train_events   events in the causally-available training pool
#'   n_train_folds    origins in that pool
#'   p_recal          the recalibrated invasion probability
#' Rows excluded from scoring (already-affected zones, non-finite p) keep p_recal = NA
#' and delta_preq = NA so they can never be mistaken for scored rows.
#'
#' Idempotent: calling it on an already-recalibrated frame refits from the raw
#' probabilities and replaces the columns, rather than accumulating .x/.y duplicates.
#'
#' @param lfo_results the LFO prediction frame.
#' @return the same frame with the five columns added (row order preserved).
attach_invasion_recalibration <- function(lfo_results, band = RECAL_BAND,
                                          min_folds = RECAL_MIN_FOLDS,
                                          burnin_fn = get0("invasion_is_burnin", mode = "function",
                                                           ifnotfound = NULL)) {
  stopifnot(all(c("method", "horizon", "fold_id", "p_invasion", "is_new_invasion")
                %in% names(lfo_results)))
  d <- lfo_results
  # Idempotency: a second call must refit from p_invasion, not join alongside a stale
  # set of columns (which dplyr would silently rename to delta_preq.x / delta_preq.y).
  .added <- c("delta_preq", "delta_estimable", "n_train_events", "n_train_folds", "p_recal")
  if (any(.added %in% names(d))) d <- d[, setdiff(names(d), .added), drop = FALSE]
  d$.row <- seq_len(nrow(d))
  sc <- .recal_scorable(d)
  if (!nrow(sc)) {
    warning("[recal] no rows are eligible for a fitted delta; every scorable row that is ",
            "rank-only takes the identity (delta = 1) and the rest are returned ",
            "unrecalibrated.", call. = FALSE)
    d$delta_preq <- NA_real_; d$delta_estimable <- NA; d$p_recal <- NA_real_
    d$.row <- NULL
    # The early return used to stop here, so when EVERY method is rank-only the whole frame
    # came back with p_recal = NA -- the very state that demotes both figure suites to raw.
    return(.recal_identity_fill(d))
  }
  # THE DELTA FIT USES THE SHARED (fold x zone) SUPPORT, NOT EACH METHOD'S NATIVE ROWS.
  # `delta_preq` exists ONLY to score -- deployment uses the separate pooled delta from
  # invasion_delta_table() -- and scoring happens on the support every method covers
  # (16_invasion_eval.R). Fitting it on native rows gave whichever method reaches furthest
  # back a nuisance parameter trained on evidence its competitors were never allowed: the
  # rolling-predictor floor costs the Bayesian grid its earliest origin, so Gravity-B4 alone
  # kept that fold and began recalibrating three folds earlier, with a delta differing by up
  # to 0.63 at h=1 and 0.74 at h=2 from the same fit on shared cells. Every recalibrated
  # column in invasion_evaluation.csv was non-comparable between it and the Bayesian grid.
  .cells <- lapply(sort(unique(sc$horizon)), function(h)
    if (exists("invasion_common_cells", mode = "function"))
      tryCatch(invasion_common_cells(lfo_results, h), error = function(e) NULL) else NULL)
  names(.cells) <- as.character(sort(unique(sc$horizon)))
  if (!exists("invasion_common_cells", mode = "function"))
    warning("[recal] invasion_common_cells() is not loaded (16_invasion_eval.R), so the ",
            "prequential delta is fitted on each method's NATIVE rows. Methods reaching ",
            "further back then train their delta on evidence the others never saw.",
            call. = FALSE)
  # ONE DELTA PER METHOD (2026-10-14): fitted on the method's 1-WEEK rows (shared cells) and
  # applied to EVERY horizon of the same fold, so P(<= 2 wk) >= P(<= 1 wk) holds row by row.
  methods <- unique(sc$method)
  parts <- lapply(methods, function(m) {
    # %in% rather than == : a NA in method would make == yield NA and pull in phantom rows.
    s  <- sc[sc$method %in% m, , drop = FALSE]                 # every horizon: the rows to adjust
    s1 <- s[s$horizon %in% 1L, , drop = FALSE]                 # the 1-week rows: the training pool
    # Fit on the shared cells; ATTACH to all of this method's scorable rows, so a row outside
    # the shared support still carries the delta of the fold it belongs to rather than NA.
    cl <- .cells[["1"]]
    s_fit <- if (length(cl))
      s1[paste(s1$fold_id, s1$health_zone, sep = "\r") %in% cl, , drop = FALSE] else s1
    if (!nrow(s_fit)) s_fit <- s1
    dl <- if (nrow(s1)) {
      prequential_invasion_delta(s, horizon = 1L, band = band, min_folds = min_folds, fit = s_fit,
                                 burnin_fn = burnin_fn)
    } else {
      warning(sprintf("[recal] %s has no 1-week rows to fit delta on; it is left unrecalibrated.", m),
              call. = FALSE)
      data.frame(fold_id = sort(unique(s$fold_id)), delta_preq = 1, delta_estimable = FALSE,
                 n_train_events = 0L, n_train_folds = 0L)
    }
    s <- dplyr::left_join(s, dl, by = "fold_id")
    # delta_preq is a per-ROW column (it varies by fold), so the transform is a
    # single vectorised call — no split/reassemble, hence no row-order hazard.
    s$p_recal <- recalibrate_invasion_p(s$p_invasion, s$delta_preq)
    s[, c(".row", "delta_preq", "delta_estimable", "n_train_events", "n_train_folds", "p_recal")]
  })
  keys <- unique(sc[, c("method", "horizon")])
  add <- dplyr::bind_rows(parts)
  d <- dplyr::left_join(d, add, by = ".row")
  d$.row <- NULL
  # RANK-ONLY ROWS TAKE THE IDENTITY, NOT NA. A comparator that emits an ordering rather than
  # a probability (prob_calibrated = FALSE) is excluded from the FIT by .recal_scorable(), and
  # that is right -- a delta fitted to an arbitrary monotone scale means nothing. But leaving
  # its p_recal at NA is not: recalibration is a strictly monotone transform, so for a
  # rank-only comparator the recalibrated and raw scales are the same ordering, and delta = 1
  # states exactly that. append_naive_detection_curve_model() (20_forecast_detail.R) already
  # uses this convention for the injected baselines; attach did not, and the mismatch cost the
  # figures dearly: 22,112 scored rows with NA p_recal made fs_lfo_col() fall back to the RAW
  # column for the WHOLE frame, so every panel in both suites was drawn raw while captioned
  # "recalibrated". evaluate_invasion() still reports their calibration block as NA -- it keys
  # on prob_calibrated, not on p_recal -- so nothing gains a proper score it should not have.
  d <- .recal_identity_fill(d)
  n_rec <- sum(is.finite(d$p_recal) & d$delta_estimable %in% TRUE)
  message(sprintf("[recal] prequential recalibration attached: %d/%d scorable rows carry a fitted delta (%d method x horizon slices)",
                  n_rec, nrow(sc), nrow(keys)))
  d
}

#' Pooled (deployment) recalibration factors, with zone-clustered 90% intervals.
#'
#' ONE factor per method, fitted on its 1-week rows (shared support) from settled, post-burn-in folds and
#' reported on one row per horizon with the SAME value (fit_horizon = 1), so every consumer that
#' reads a (method, horizon) table -- the live suite, the cascade's cascade_fit_delta() -- keeps
#' working and receives the single factor. Appropriate ONLY for recalibrating a live forecast.
#' Never use these for scoring: see the module header on in-sample optimism.
#'
#' @return tibble(method, horizon, fit_horizon, delta, delta_lo, delta_hi, se_log_delta,
#'   n_rows, n_events, sum_hazard, n_folds, n_boot_used, cal_in_large, delta_moment,
#'   delta_preq_final, gamma_shape, gamma_rate, boundary_hit)
invasion_delta_table <- function(lfo_results, band = RECAL_BAND,
                                 min_folds = RECAL_MIN_FOLDS,
                                 n_boot = RECAL_N_BOOT,
                                 # FIT THE DEPLOYED FACTOR ON SETTLED OUTCOMES ONLY.
                                 # delta is essentially (events + 1/2) / sum(hazard): it is
                                 # driven by the EVENT COUNT. Since 2026-09-22 the
                                 # cross-validation runs to the last round whose outcome window
                                 # closes by the analysis date, so its newest rounds are
                                 # right-truncated — they contribute their full hazard mass to
                                 # the denominator while some of their invasions are not yet
                                 # laboratory-confirmed and so are missing from the numerator.
                                 # Including them would bias delta DOWNWARD and deflate every
                                 # deployed invasion probability, for a purely artefactual
                                 # reason. Rounds are labelled by run_invasion_lfo() via
                                 # `eval_reliable`; when that column is absent (an older LFO
                                 # frame) every row is used, reproducing the previous behaviour.
                                 # Set FALSE to fit on every scored round.
                                 reliable_only = TRUE,
                                 burnin_fn = get0("invasion_is_burnin", mode = "function",
                                                  ifnotfound = NULL)) {
  sc <- .recal_scorable(lfo_results)
  if (!nrow(sc)) return(tibble::tibble())
  # `sc_all` keeps EVERY scored round. The restriction below applies ONLY to the pooled
  # deployment factor; the prequential trace is computed from sc_all, because each of its
  # per-fold values is a genuine "what a live system held at that moment" estimate and
  # `delta_preq_final` must remain the value in force at the LAST round, not the last
  # settled one.
  sc_all <- sc
  # THE TRAINING POOL (2026-10-14): 1-WEEK rows, post-burn-in rounds, settled outcomes. The same
  # factor is then reported for -- and applied to -- every horizon of the method.
  sc <- sc[sc$horizon %in% 1L, , drop = FALSE]
  if ("cutoff" %in% names(sc) && is.function(burnin_fn)) {
    .bi <- burnin_fn(sc$cutoff)
    if (any(.bi))
      message(sprintf("[recal] deployment delta excludes %d burn-in round(s) (week index < %d).",
                      dplyr::n_distinct(sc$fold_id[.bi]),
                      as.integer(get0("ASCERTAINMENT_BURNIN_WEEKS", ifnotfound = 7L))))
    sc <- sc[!.bi, , drop = FALSE]
  }
  if (isTRUE(reliable_only) && "eval_reliable" %in% names(sc)) {
    n_all <- nrow(sc)
    sc_r  <- sc[sc$eval_reliable %in% TRUE, , drop = FALSE]
    # Never let the restriction empty the pool: a frame in which NO round is settled yet
    # would otherwise silently return an empty delta table and leave the live forecast
    # unrecalibrated without saying so.
    if (nrow(sc_r) && sum(sc_r$is_new_invasion, na.rm = TRUE) > 0) {
      if (nrow(sc_r) < n_all)
        message(sprintf(paste0("[recal] deployment delta fitted on SETTLED rounds only: ",
                               "%d of %d scorable 1-week rows (%d of %d rounds)."),
                        nrow(sc_r), n_all, dplyr::n_distinct(sc_r$fold_id),
                        dplyr::n_distinct(sc$fold_id)))
      sc <- sc_r
    } else {
      warning("[recal] no settled round carries an invasion; the deployment delta is fitted ",
              "on ALL scored post-burn-in rounds, including right-truncated ones.", call. = FALSE)
    }
  }
  methods <- unique(sc_all$method)
  .cells1 <- if (exists("invasion_common_cells", mode = "function"))
    tryCatch(invasion_common_cells(lfo_results, 1L), error = function(e) NULL) else NULL
  rows <- lapply(methods, function(m) {
    # %in% rather than == : a NA in method would make == yield NA and pull in phantom rows.
    s <- sc[sc$method %in% m, , drop = FALSE]
    # ON THE SHARED (fold x zone) SUPPORT, like delta_preq and the stability diagnostic, so the
    # deployed factor is fitted on the same cells every method is scored on. A method with no
    # rows on that support keeps its native rows (and says so).
    if (!is.null(.cells1) && nrow(s)) {
      .on <- paste(s$fold_id, s$health_zone, sep = "\r") %in% .cells1
      if (any(.on)) s <- s[.on, , drop = FALSE]
      else warning(sprintf("[recal] %s has no 1-week rows on the shared support; its deployment delta uses its native rows.", m),
                   call. = FALSE)
    }
    hs <- sort(unique(sc_all$horizon[sc_all$method %in% m]))
    if (!nrow(s)) {
      warning(sprintf("[recal] %s has no usable 1-week rows; no deployment delta.", m), call. = FALSE)
      return(NULL)
    }
    y <- as.integer(s$is_new_invasion); p <- s$p_invasion
    d  <- fit_invasion_delta(p, y, band = band)
    ci <- if (is.finite(d) && "health_zone" %in% names(s))
            bootstrap_invasion_delta(p, y, s$health_zone, band = band,
                                     n_boot = n_boot)
          else list(lo = NA_real_, hi = NA_real_, n_boot_used = 0L, se_log = NA_real_)
    base <- mean(y)
    cil  <- if (base > 0) mean(p) / base else NA_real_
    # delta_preq_final = the delta a live system would have been carrying at the LAST fold,
    # i.e. the prequential value in force at the end of the record, under the same rule.
    # Trained on the SHARED 1-week zone-rounds, exactly as attach_invasion_recalibration() trains
    # delta_preq, so this is the last value of THAT trace (it used to be fitted on the method's
    # native rows, which can extend beyond the shared support).
    s_pq <- sc_all[sc_all$method %in% m, , drop = FALSE]
    .fit1 <- s_pq[s_pq$horizon %in% 1L, , drop = FALSE]
    if (!is.null(.cells1))
      .fit1 <- .fit1[paste(.fit1$fold_id, .fit1$health_zone, sep = "\r") %in% .cells1, , drop = FALSE]
    pq <- prequential_invasion_delta(s_pq, horizon = 1L, band = band, min_folds = min_folds,
                                     fit = .fit1, burnin_fn = burnin_fn)
    # GAMMA PARAMETERS FOR THE LIVE DRAWS. Under the Jeffreys prior on log(delta) the
    # rare-event posterior of delta is Gamma(events + 1/2, sum of hazards); the rate is set so
    # the mean equals the fitted delta exactly (the Bernoulli fit and its Poisson limit
    # (E + 1/2) / sum(hazard) differ by about 1-3% on these pools).
    .shape <- sum(y) + 0.5
    tibble::tibble(
      method = m, horizon = hs, fit_horizon = 1L,
      delta = d, delta_lo = ci$lo, delta_hi = ci$hi, se_log_delta = ci$se_log,
      n_rows = nrow(s), n_events = sum(y),
      sum_hazard = sum(-log1p(-pmin(pmax(p, 0), 1 - 1e-12))),
      n_folds = dplyr::n_distinct(s$fold_id),
      n_boot_used = ci$n_boot_used,
      cal_in_large = cil,
      delta_moment = if (is.finite(cil) && cil > 0) 1 / cil else NA_real_,
      delta_preq_final = pq$delta_preq[which.max(pq$fold_id)],
      gamma_shape = .shape,
      gamma_rate  = if (is.finite(d) && d > 0) .shape / d else NA_real_,
      boundary_hit = is.finite(d) &&
        (abs(log(d) - log(band[1])) < 1e-4 || abs(log(d) - log(band[2])) < 1e-4))
  })
  dplyr::bind_rows(rows) %>% dplyr::arrange(horizon, method)
}

#' Is the calibration offset STABLE over time? (the delta-stability diagnostic)
#'
#' Fits delta INDEPENDENTLY on each fold -- not cumulatively -- because a cumulative series
#' cannot distinguish genuine drift from an expanding window converging, and that distinction
#' is the whole question. Then asks, per (method, horizon), whether the per-fold values differ
#' by more than their own sampling error.
#'
#' THE POOL IS THE DEPLOYED ONE. The deployment delta (invasion_delta_table) is trained on
#' post-burn-in rounds whose outcomes are settled (`eval_reliable`), so the folds tested here
#' are exactly those: a burn-in round is excluded from training by design, and an unsettled
#' round is missing confirmations, so either would show up as "drift" the deployed factor never
#' sees. Rows are further restricted to the shared (fold x zone) support so that every method is
#' diagnosed on the same cells. The deployed factor is fitted on the 1-WEEK rows only; the
#' h = 1 summary therefore speaks to the deployed pool, while the h = 2 summary asks the
#' separate question of whether the 2-week offset is itself stable.
#'
#' WHY IT MATTERS. delta is an intercept shift on the cloglog scale and beta_0 is the
#' intercept, so a per-fold delta IS a time-varying beta_0 fitted post hoc. This diagnostic is
#' therefore the direct evidence on two questions at once:
#'   * should the deployment delta use all history, or a recent window?  A window only earns
#'     its cost if delta genuinely moves.
#'   * is a time-varying beta_0 worth fitting?  This is the quantity such a model would be
#'     fitted to, free of the kernel and the covariates.
#'
#' The heterogeneity test is DerSimonian-Laird random effects on log delta:
#'   Q   = sum w_k (y_k - ybar)^2,  w_k = 1/se_k^2      (Cochran; chi-sq on n_folds - 1 df)
#'   tau2= max(0, (Q - (k-1)) / (sum w - sum w^2 / sum w))   genuine between-fold variance
#'   I2  = max(0, (Q - (k-1)) / Q)                            share NOT from sampling noise
#' se(log delta) comes from the EXPECTED (Fisher) information of theta = log delta,
#' I(theta) = sum lambda^2 e^-lambda / (1 - e^-lambda) with lambda = delta * u, evaluated at the
#' fitted value (it does not depend on the outcomes, so it is not the observed information).
#'
#' COVARIATES. Two quantities that move with epidemic phase are regressed on: the fold index (a
#' trend test) and the fold's at-risk count (epidemic extent -- the zone set shrinks as zones
#' are invaded; a property of the fold design, not of the predictions or the outcome, so the
#' regression is not circular). The nowcast multiplier applied to each round is NOT carried in
#' the LFO frame and is not tested; the round's `training_source` (archived snapshot vs
#' re-censored final line list) is reported per fold, so a shift between the two regimes is
#' visible in the per-fold table.
#'
#' @param burnin_fn function(cutoff) -> TRUE for burn-in rounds (excluded, as in training).
#' @return list(per_fold = tibble, summary = tibble). Written to CSV by run_all.R and READ by
#'   the figure; no plotting code recomputes any of it.
invasion_delta_stability <- function(lfo_results, band = RECAL_BAND,
                                     burnin_fn = get0("invasion_is_burnin", mode = "function",
                                                      ifnotfound = NULL)) {
  sc <- .recal_scorable(lfo_results)
  if (!nrow(sc)) return(list(per_fold = tibble::tibble(), summary = tibble::tibble()))
  # The deployed training pool: post-burn-in rounds with settled outcomes (see invasion_delta_table).
  if ("cutoff" %in% names(sc) && is.function(burnin_fn))
    sc <- sc[!burnin_fn(sc$cutoff), , drop = FALSE]
  if ("eval_reliable" %in% names(sc))
    sc <- sc[sc$eval_reliable %in% TRUE, , drop = FALSE]
  if (!nrow(sc)) return(list(per_fold = tibble::tibble(), summary = tibble::tibble()))
  # Restrict to the shared support, for the same reason the prequential fit does: a method
  # reaching further back would otherwise be diagnosed on cells its competitors never saw.
  if (exists("invasion_common_cells", mode = "function")) {
    keep <- unlist(lapply(sort(unique(sc$horizon)), function(h) {
      cl <- tryCatch(invasion_common_cells(lfo_results, h), error = function(e) NULL)
      if (!length(cl)) character(0) else paste(h, cl, sep = "\r")
    }), use.names = FALSE)
    if (length(keep))
      sc <- sc[paste(sc$horizon, sc$fold_id, sc$health_zone, sep = "\r") %in% keep, , drop = FALSE]
  }
  if (!nrow(sc)) return(list(per_fold = tibble::tibble(), summary = tibble::tibble()))

  # EPIDEMIC EXTENT. `was_active_before` is hard-coded FALSE in the LFO frame -- already-
  # affected zones are dropped upstream rather than flagged -- so counting it gives nothing.
  # The AT-RISK count per fold is the usable proxy: the zone set shrinks monotonically as
  # zones are invaded, so it is a monotone transform of extent, and it is a property of the
  # fold DESIGN rather than of the model's predictions or the outcome, so regressing delta on
  # it is not circular.

  .z90 <- stats::qnorm(0.95)
  .one <- function(d) {                      # delta + se(log delta) on ONE fold's rows
    y <- as.integer(d$is_new_invasion)
    u <- -log1p(-pmin(pmax(d$p_invasion, 0), 1 - 1e-12))
    dd <- fit_invasion_delta(d$p_invasion, y, band = band, warn_boundary = FALSE)
    if (!is.finite(dd)) return(c(NA_real_, NA_real_, sum(y, na.rm = TRUE), sum(u)))
    lam <- dd * u
    I <- sum(ifelse(lam < 1e-8, lam, lam^2 * exp(-lam) / -expm1(-lam)))
    c(dd, 1 / sqrt(max(I, 1e-12)), sum(y, na.rm = TRUE), sum(u))
  }
  keys <- unique(sc[, c("method", "horizon")])
  pf <- dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(i) {
    m <- keys$method[i]; h <- keys$horizon[i]
    s <- sc[sc$method %in% m & sc$horizon %in% h, , drop = FALSE]
    fid <- sort(unique(s$fold_id))
    v <- t(vapply(fid, function(k) .one(s[s$fold_id == k, , drop = FALSE]), numeric(4)))
    tibble::tibble(method = m, horizon = h, fold_id = fid,
                   cutoff = as.Date(vapply(fid, function(k) {
                     cc <- unique(as.Date(s$cutoff[s$fold_id == k])); cc <- cc[!is.na(cc)]
                     if (length(cc)) as.character(min(cc)) else NA_character_ }, character(1))),
                   training_source = if ("training_source" %in% names(s))
                     vapply(fid, function(k) paste(sort(unique(as.character(
                       s$training_source[s$fold_id == k]))), collapse = "+"), character(1))
                     else NA_character_,
                   n_atrisk = as.integer(vapply(fid, function(k) sum(s$fold_id == k), integer(1))),
                   n_events = as.integer(v[, 3]), expected = v[, 4],
                   delta = v[, 1], se_log = v[, 2],
                   # 90% Wald interval on log delta: the pipeline's interval convention.
                   lo = v[, 1] * exp(-.z90 * v[, 2]), hi = v[, 1] * exp(.z90 * v[, 2]))
  }))


  sm <- dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(i) {
    m <- keys$method[i]; h <- keys$horizon[i]
    z <- pf[pf$method == m & pf$horizon == h & is.finite(pf$delta) &
            is.finite(pf$se_log) & pf$se_log > 0, , drop = FALSE]
    k <- nrow(z)
    if (k < 2L) return(tibble::tibble(method = m, horizon = h, n_folds = k))
    y <- log(z$delta); v <- z$se_log^2; w <- 1 / v
    mu0 <- sum(w * y) / sum(w); Q <- sum(w * (y - mu0)^2)
    Cc  <- sum(w) - sum(w^2) / sum(w)
    tau2 <- max(0, (Q - (k - 1)) / Cc)
    w2 <- 1 / (v + tau2); mu <- sum(w2 * y) / sum(w2); se_mu <- sqrt(1 / sum(w2))
    .slope <- function(x) {
      if (length(unique(x[is.finite(x)])) < 2L) return(c(NA_real_, NA_real_))
      f <- try(stats::lm(y ~ x, weights = w), silent = TRUE)
      if (inherits(f, "try-error")) return(c(NA_real_, NA_real_))
      cf <- summary(f)$coefficients
      if (nrow(cf) < 2L) c(NA_real_, NA_real_) else c(cf[2, 1], cf[2, 4])
    }
    st <- .slope(seq_len(k)); sa <- .slope(as.numeric(z$n_atrisk))
    tibble::tibble(
      method = m, horizon = h, n_folds = k,
      pooled_delta = exp(mu), pooled_lo = exp(mu - .z90 * se_mu),
      pooled_hi = exp(mu + .z90 * se_mu),
      tau2 = tau2, I2 = max(0, (Q - (k - 1)) / Q),
      Q = Q, Q_df = k - 1L, Q_p = stats::pchisq(Q, k - 1L, lower.tail = FALSE),
      trend_slope = st[1], trend_p = st[2],
      extent_slope = sa[1], extent_p = sa[2],
      median_se_log = stats::median(z$se_log))
  }))
  list(per_fold = pf, summary = sm)
}

message("[16b] invasion recalibration loaded (Jeffreys-penalised ML, min_folds=",
        RECAL_MIN_FOLDS,
        ", band=[", RECAL_BAND[1], ", ", RECAL_BAND[2], "])")
