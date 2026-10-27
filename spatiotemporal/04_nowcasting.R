# =============================================================================
# 04_nowcasting.R — Right-Truncation Correction (Nowcasting)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: Correct the zone × week observation matrix for right-truncation
#   arising because recent onsets have not yet reached the line list.
#   The correction is re-applied at each LFO-CV fold, at that fold's own nowcast date, to
#   counts holding only what was known on the fold's origin. The TRUNCATION a fold is
#   corrected with is not always real-time: a re-censored round uses the shared "asof" fit,
#   estimated on the final line list (see the RESIDUAL note under assumption 1).
#
# ASSUMPTIONS THIS CORRECTION MAKES — stated, because none of them is verified here.
#
#   1. THE DELAY MUST MATCH THE LEG THAT ACTUALLY TRUNCATES THE SERIES, and which leg that is
#      depends on how the series was built. Since 2026-09-21 every pipeline caller passes an
#      explicitly FITTED truncation (epinow2_truncation(), 02_epi_params.R section 4b) instead
#      of relying on this file's default, and there are two regimes because there are two
#      constructions:
#        - "extract" — the DEPLOYED series (run_all.R, run_cascade.R) is the current extract.
#          A case appears in it only after sample -> laboratory result -> data entry, so the
#          truncating leg is onset -> APPEARANCE-IN-EXTRACT, deconvolved from the vintage
#          archive (EpiNow2::estimate_truncation(), a lognormal).
#          A back-test round that reads an ARCHIVED SNAPSHOT (LFO_TRAINING_SOURCE = "vintage",
#          the default) is in this regime too, with the extract truncation fitted on the
#          snapshots up to its own date (vintage_truncation(), 22_daily_reissue.R).
#        - "asof" — a back-test round or cascade calibration origin that reads no snapshot (none
#          within LFO_VINTAGE_MAX_GAP_DAYS of its origin, or too few yet for a real-time
#          truncation fit) is rebuilt from the final line list by reaggregate_asof(), which
#          censors on linelist_observation_date() — the SAMPLE-collection date for all but a
#          handful of confirmed records (those fall through to the lab or notification date).
#          The truncating leg there is onset -> OBSERVATION: shorter, because it does not wait
#          for the laboratory result and data entry. Measured by onset cohort rather than
#          deconvolved, since both dates are recorded, and fitted as a gamma
#          (.trunc_fit_asof_cohort()).
#      The weight each regime gives the newest week is reported every time the correction
#      runs (the minimum of the weight range compute_truncation_weights() messages), and
#      52_vintage_level_check.R compares both reconstructions' nowcasts with the counts those
#      weeks eventually reached.
#      HISTORY. Until 2026-09-21 both paths used the onset->SAMPLE delay. For the DEPLOYED series
#      that is the wrong leg, and the error is one-directional: measured then against 54
#      archived vintages (5,607 records matched to their first-appearance vintage),
#      onset->appearance had mean 11.43 d against the modelled 7.67 d, so the most recent week —
#      the one every forecast is issued FROM — was under-corrected, and the shortfall
#      propagated into every h=1 invasion probability. The first fitted extract truncation
#      (2026-09-21) over-corrected in the other direction: until 2026-10-14 its panel dated cases
#      by RECORDED onset while the corrected series is dated by date_index, and it multiplied
#      the newest live week about twice too much (02_epi_params.R, SERIES BASIS).
#      THE TWO REGIMES ARE DELIBERATELY DIFFERENT DISTRIBUTIONS. The invariant that matters is
#      not "same multiplier both sides" but "each series corrected for its own truncation, so
#      both estimate the same quantity — true onsets per week". Forcing one distribution on
#      both would necessarily bias one of them.
#      RESIDUAL, FOR THE RE-CENSORED ROUNDS ONLY: the as-of reconstruction is OPTIMISTIC, twice.
#      (a) At a real origin T you would not yet have known about a case sampled at T-1 but not
#      yet confirmed and entered; reaggregate_asof() assumes instant lab-and-data-entry and
#      shows it anyway. (b) Its truncation is ONE shared fit on the final line list
#      (epinow2_truncation(), "asof" regime), so it uses rows that postdate the origin —
#      look-ahead in a nuisance parameter, accepted there and stated here. Rounds that read an
#      archived snapshot (by default, every round with a snapshot within
#      LFO_VINTAGE_MAX_GAP_DAYS AND a real-time truncation fit; the earliest rounds lack one
#      or the other and re-censor) have neither gap: they see exactly what the extract held on
#      that date, corrected with the truncation estimable on that date.
#      This file's own default is STILL the onset->sample delay (it is the only delay available
#      without an archive), so leaving `delay` NULL now WARNS once per session.
#
#   2. THE FITTED TRUNCATIONS ARE CONFIRMED-ONLY, like the series they correct: the extract
#      panel tabulates confirmed records (.trunc_frame_vintage()) and the as-of cohort keeps
#      confirmed records (.trunc_fit_asof_cohort()). The wrong-leg FALLBACK
#      (effective_onset_sample_delay() with no stratum) is the POOLED fit over every case
#      classification, which is shorter than the confirmed stratum (.trunc_fit_asof_cohort()
#      documents the contrast), so a fallback correction under-corrects for a second reason.
#
#   3. THE SUSPECTED SERIES IS CORRECTED WITH THE SAME CDF. A suspected case's appearance is
#      governed by alert notification, not by specimen collection, so `suspected_nc` carries a
#      delay that was not fitted for it. It feeds only the suspected-covariate models, which
#      are OFF by default (INCLUDE_SUSPECTED_COV_MODELS).
#
#   4. NO WEEK IS "FULLY OBSERVED" UNLESS THE SPEC CARRIES ITS SUPPORT BOUND. A fitted family's
#      CDF never reaches 1, so every historical week keeps a multiplier 1/F(lag + 0.5) above 1
#      (largest for the heavier-tailed lognormal) and every confirmed_nc is non-integer.
#      EpiNow2 fits the truncation CONDITIONED on its maximum delay (its pmf is renormalised on
#      [0, max]); when the spec records that bound as `trunc_max`, .nowcast_cdf() renormalises
#      the same way, so F reaches 1 at the bound and older weeks are left uncorrected. A spec
#      without it (the onset->sample fallback, a `rate` override) is evaluated
#      unconditionally, and its residual 1 - F is an inflation applied across the whole series.
#
#   5. THE CORRECTION IS SCORED, BUT NOT INSIDE THE FORECAST. 52_vintage_level_check.R (run by
#      run_all.R by default) compares each reconstruction's nowcast of the origin week, and of
#      the last three weeks, with the count those weeks eventually reached in the final line
#      list; nothing feeds that back into the correction. compare_epinowcast_models() exists
#      but has no caller. What the LFO itself fixes is the GEOMETRY: since 2026-09-19 the fold
#      origin is cut + 6, the last day of the last training week, which is the same geometry as
#      the live grid, so a fold's most recent week occupies the same position relative to its
#      origin as the deployed current week does. (A snapshot round is nowcast at its snapshot's
#      date, 0-2 days before the origin.) The MULTIPLIERS differ between regimes, and should —
#      see assumption 1.
#
# Key statistical context:
#   F is the CDF of whatever delay the CALLER passes — in the pipeline, the fitted truncation
#   for that caller's regime (assumption 1). It is evaluated on its OWN family by delay_cdf()
#   (00_config.R), through .nowcast_cdf() below: lnorm for the extract truncation, gamma for the
#   as-of truncation and for the retired onset->sample fit. When no delay is passed this falls
#   back to effective_onset_sample_delay() (00_config.R) — the pooled EpiDist MARGINAL fit from
#   04c_dhis2_delay_windows.R (a gamma, truncation- and double-interval-censoring corrected;
#   its parameters are in the params file 04c writes), then the windowed interval-censored
#   MLE, then the fixed lab reference Exp(DELAY_ONSET_SAMPLE_RATE). That fallback is the WRONG
#   LEG for a truncation correction and warns.
#
#   For a week w, the truncation weight is the expected fraction of that week's onsets
#   already observed by the analysis date — averaged over the week's seven onset days,
#   NOT evaluated once at the week midpoint:
#
#     weight_w = (1/7) * sum_{j=0}^{6} P(delay <= analysis_date - (w + j))
#              = (1/7) * sum_{j=0}^{6} F(analysis_date - (w + j) + 0.5)
#
#   (the + 0.5 is the daily-rounding convention: P(delay <= L) = F(L + 0.5) for an integer
#   lag L. It is the convention 04c fits the onset->sample delay under and the as-of cohort fit
#   matches its gamma under. The extract truncation is EpiNow2's primary-censored
#   discretisation of its lognormal, which F(L + 0.5) approximates rather than reproduces.)
#   The nowcast-corrected count is:
#
#     confirmed_nc = confirmed / weight_w
#
#   This is a simple inverse-probability-of-reporting correction.  It is
#   analytic and fast, but ignores delay uncertainty (04b_epinowcast.R is the probabilistic
#   alternative, run only as a sensitivity arm and diagnostic; see apply_nowcast_correction()).
#
# Dependencies: 00_config.R (sourced first), tidyverse, lubridate
#
# Usage (typical; pass the fitted truncation for the series' regime):
#   source("spatiotemporal/04_nowcasting.R")
#   zw_corrected <- apply_nowcast_correction(dat$zone_week, ANALYSIS_DATE,
#                                            delay = .trunc_as_delay_spec(
#                                              epinow2_truncation(ll = dat$ll), "extract"))
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# =============================================================================
# SECTION 1: compute_truncation_weights()
# =============================================================================
#
# Computes a named vector of right-truncation (completeness) weights, one per week.
#
# DEFINITION. For a week w, the weight is the expected fraction of that week's onsets
# already sampled — and therefore already visible in the line list — by `analysis_date`:
#
#     weight_w = E_o[ P(Delta <= analysis_date - o) ]      o = onset date, uniform over w
#              = (1/7) * sum_{j=0}^{6} F(analysis_date - (w + j) + 0.5)
#
# Two corrections vs the previous implementation:
#
#   1. THE DELAY (the substantive one). It hard-defaulted to the FIXED lab reference
#      Exp(DELAY_ONSET_SAMPLE_RATE) (rate 0.228/d, mean 4.39 d) even though the header
#      claimed the data-estimated rate was used when ONSET_SAMPLE_DELAY_SOURCE = "data".
#      Nothing ever substituted it, so every fold was corrected with a mean delay of
#      4.39 d against the 7.665 d fit then in force — about 43% too short, i.e. reporting
#      assumed far more complete than it is, so the correction was too small. It now defaults
#      to effective_onset_sample_delay() (00_config.R) — the EpiDist marginal fit
#      (truncation- and double-interval-censoring corrected) when 04c wrote one, else the
#      windowed interval-censored MLE, else the lab reference — and evaluates that fit's
#      ACTUAL family (delay_cdf(); a gamma for the EpiDist fit) rather than moment-matching it
#      to an Exponential. SUPERSEDED 2026-09-21: that fixed the delay's magnitude but not its
#      LEG. The deployed series is truncated by onset->appearance, not onset->sample, so the
#      callers now pass a fitted truncation and this resolver is only the fallback (assumption
#      1). The onset->sample fit remains the right object for the missing-onset IMPUTATION,
#      which really is imputing to sample; imputation and nowcast are therefore no longer the
#      same distribution, and correctly so.
#      Validated when that fix was made, against the untruncated empirical delay (onsets
#      >= 35 d before the as-of date, n=3821, mean 7.21 d): mean |CDF error| over lags 0-30 d
#      was 0.018 for the fitted delay vs 0.087 for the lab reference.
#
#   2. THE WITHIN-WEEK AVERAGE (second-order). It evaluated F once at week_start + 3.5,
#      i.e. F(E[lag]) rather than the required E[F(lag)]. Two distinct errors were folded
#      into that single evaluation and they point OPPOSITE ways, so the retired weight was
#      not biased in a predictable direction:
#        (a) Jensen: F is strongly concave for a shape<1 delay, so F(E[lag]) > E[F(lag)]
#            — overstates completeness.
#        (b) week_start + 3.5 is the midpoint of a CONTINUOUS 7-day interval, but the mean
#            of the seven integer onset DATES is week_start + 3, and the daily-rounding
#            convention adds a further +0.5 to the lag. The retired lag was therefore a
#            full day SHORT — understates completeness.
#      Holding the delay fixed, (b) dominated (a) on the snapshot this was fixed on
#      (2026-09-19), so the convention change partly offset fix (1); the deployed current week
#      then carried weight 0.3801 under the retired onset->sample delay, a 2.631x correction.
#      The multiplier cap is derived from the delay rather than fixed at 5x (section 2) so that
#      a larger fitted correction remains a correction and not a clip.
#   We average F over the week's seven onset days, with P(Delta <= L) = F(L + 0.5) for integer
#   L (the interval-censoring convention 04c fits under, d -> [d-0.5, d+0.5]; see the header
#   for how it relates to each fitted truncation). F is .nowcast_cdf(): delay_cdf() on the
#   spec's own family, renormalised on [0, trunc_max] when the spec carries that bound
#   (assumption 4). Onset days after `analysis_date` contribute 0 (not yet observable) but stay
#   in the denominator, since they belong to the week whose total is being estimated.
#
# Arguments:
#   week_start_dates — Date vector; the start of each 7-day bucket (WEEK_ANCHOR-anchored)
#   analysis_date    — Date scalar; reference date for the correction
#                      (typically ANALYSIS_DATE from 00_config.R)
#   delay            — delay spec: either a fitted truncation (.trunc_as_delay_spec(), the
#                      pipeline path) or an effective_onset_sample_delay() spec. Leaving it
#                      NULL falls back to the latter, which is the WRONG LEG for a truncation
#                      correction, and warns once per session.
#   rate             — optional Exponential-rate OVERRIDE (day^-1). Mutually exclusive with
#                      `delay`; provided so a caller (or a unit test) can pin a plain
#                      Exp(rate) delay. Passing both is an error rather than a silent
#                      precedence rule — silent precedence is how bug (1) survived.
#   week_days        — bucket length in days (7; exposed only so the averaging is not
#                      silently wrong if the weekly grid ever changes).
#
# Returns:
#   Named numeric vector; names are character representations of week_start_dates.
#   Values are in (0, 1] for observed weeks, NA for weeks that start after the analysis date
#   (nothing in them is observable yet). apply_nowcast_correction() additionally sets the
#   corrected counts (not the weight) to NA for any week whose mean lag is below its
#   min_lag_days.
#
# Warnings:
#   - Emits a warning for any week where the mean lag < 3 days (weight is very uncertain
#     and the correction can be extremely large).

#' Lag from a week's MEAN DAY to the analysis date, in days.
#'
#' ONE definition of the week-geometry convention, shared by compute_truncation_weights()
#' (the "very uncertain" warning) and apply_nowcast_correction() (the min_lag_days
#' suppression). They used to disagree: (week_days - 1) / 2 in one, a hard-coded + 3 in the
#' other — the same number for a 7-day week, and silently different for any other.
#'
#' NOT how the truncation WEIGHT is computed. The weight averages F() over the week's seven
#' onset days — E[F(lag)], not F(E[lag]) — because F is strongly concave for a shape < 1
#' delay. See the derivation in section 1.
.week_lag_days <- function(week_start_dates, analysis_date, week_days = 7L) {
  as.numeric(analysis_date - (week_start_dates + (week_days - 1) / 2))
}

#' CDF of the delay the nowcast corrects with: delay_cdf() (00_config.R) on the spec's own
#' family, RENORMALISED on [0, trunc_max] when the spec carries that bound.
#'
#' EpiNow2 fits a truncation CONDITIONED on its maximum delay (the pmf is renormalised on
#' [0, max]), so the fitted parameters describe P(delay <= L | delay <= max), and evaluating
#' the unconditional F(L) instead understates every week's completeness by the factor F(max).
#' A spec that records the bound as `trunc_max` (element or attribute) is therefore divided by
#' its CDF at the support's upper edge, trunc_max + 0.5 under this file's daily-rounding
#' convention (P(delay <= L) = F(L + 0.5)), and reaches 1 there. A spec without it (the
#' onset->sample fallback, a `rate` override, a hand-built fixture) is evaluated exactly as
#' delay_cdf() evaluates it.
#' The CDF is taken from the spec with `trunc_max` REMOVED, so the result is the same whether
#' or not delay_cdf() itself ever learns to honour the bound: never renormalised twice.
#' Used for the weights AND for the multiplier cap, so the two cannot describe different
#' distributions.
.nowcast_cdf <- function(delay, x) {
  tm <- suppressWarnings(as.numeric(delay$trunc_max %||% attr(delay, "trunc_max")))
  if (length(tm) != 1L || !is.finite(tm) || tm <= 0) return(delay_cdf(delay, x))
  bare <- delay
  bare$trunc_max <- NULL
  attr(bare, "trunc_max") <- NULL
  z <- delay_cdf(bare, tm + 0.5)
  if (length(z) != 1L || !is.finite(z) || z <= 0)
    stop("[nowcast] the delay's CDF at its trunc_max (", tm, " d) is not positive; ",
         "cannot renormalise.", call. = FALSE)
  pmin(delay_cdf(bare, x) / z, 1)
}

# Session-scoped flags for warnings that must fire once, not once per call. parent =
# emptyenv() so a typo'd field cannot silently resolve up the search path to a global.
.nowcast_state <- new.env(parent = emptyenv())

compute_truncation_weights <- function(
    week_start_dates,
    analysis_date,
    delay     = NULL,
    rate      = NULL,
    week_days = 7L
) {

  # ---- input validation ------------------------------------------------------
  if (!inherits(week_start_dates, "Date")) {
    week_start_dates <- tryCatch(
      as.Date(week_start_dates),
      error = function(e) stop("[compute_truncation_weights] Cannot coerce week_start_dates to Date: ",
                               e$message, call. = FALSE)
    )
  }
  if (length(week_start_dates) == 0) {
    stop("[compute_truncation_weights] week_start_dates is empty.", call. = FALSE)
  }
  if (!inherits(analysis_date, "Date") || length(analysis_date) != 1) {
    stop("[compute_truncation_weights] analysis_date must be a scalar Date.", call. = FALSE)
  }
  if (!is.numeric(week_days) || length(week_days) != 1 || week_days < 1) {
    stop("[compute_truncation_weights] week_days must be a single positive integer.", call. = FALSE)
  }

  # ---- resolve the delay distribution ---------------------------------------
  if (!is.null(delay) && !is.null(rate)) {
    stop("[compute_truncation_weights] Pass `delay` OR `rate`, not both.", call. = FALSE)
  }
  if (!is.null(rate)) {
    if (!is.numeric(rate) || length(rate) != 1 || !is.finite(rate) || rate <= 0) {
      stop("[compute_truncation_weights] rate must be a single positive number.", call. = FALSE)
    }
    delay <- list(family = "exp", params = c(rate = rate), rate = rate,
                  mean = 1 / rate, sd = 1 / rate,
                  estimator = "explicit_rate_override", source = "override")
  } else if (is.null(delay)) {
    # WRONG-LEG FALLBACK. effective_onset_sample_delay() is onset->SAMPLE; what truncates
    # either series this pipeline corrects is onset->appearance (deployed) or onset->observation
    # (as-of folds). It is kept as the default only because it is the one delay obtainable
    # without the vintage archive, and because tests and ad-hoc calls need some default.
    # It WARNS because silence is how the previous version of this bug survived: the header
    # claimed a substitution that no caller ever made, and every fold was corrected with a
    # delay about 43% too short for four months with nothing in the log to say so. Once per session
    # rather than per call — the LFO would otherwise emit one per fold and bury it.
    if (is.null(.nowcast_state$warned_default_delay)) {
      .nowcast_state$warned_default_delay <- TRUE
      warning("[compute_truncation_weights] no `delay` passed; falling back to the ",
              "onset->SAMPLE delay. That is the wrong leg for a right-truncation correction ",
              "and UNDER-corrects the most recent weeks. Pipeline callers should pass the ",
              "fitted truncation for their regime (epinow2_truncation(); see assumption 1 in ",
              "04_nowcasting.R). This warning is emitted once per session.", call. = FALSE)
    }
    delay <- effective_onset_sample_delay()
  }
  if (!is.list(delay) || is.null(delay$family) || !is.numeric(delay$rate) ||
      !is.finite(delay$rate) || delay$rate <= 0) {
    stop("[compute_truncation_weights] `delay` must be a spec with a positive `rate` ",
         "(see effective_onset_sample_delay()).", call. = FALSE)
  }

  # ---- mean lag per week (used for the future-week guard and the warning) ---
  # Mean over the bucket's DAYS (j = 0..week_days-1), i.e. week_start + (week_days-1)/2.
  # For a 7-day bucket that is week_start + 3 — the mean of seven integer onset dates —
  # not the + 3.5 of a continuous-time midpoint.
  lag_days <- .week_lag_days(week_start_dates, analysis_date, week_days)

  # ---- flag weeks with NOTHING observable yet -------------------------------
  # OBSERVABILITY IS A PROPERTY OF THE WEEK'S FIRST DAY, not of its mean day. A week is in
  # the future only when none of its days has occurred by the analysis date; a week that has
  # begun is partially observable and its weight — which averages F() over the days that
  # HAVE occurred, with the rest contributing 0 — is well defined and small.
  #
  # The old rule (`mean lag <= 0`, i.e. analysis_date <= week_start + 3) returned NA for any
  # week whose midpoint had not passed, so a week with four observable days was dropped
  # rather than corrected. It is latent under the WEEK_ANCHOR grid, where the newest week
  # starts analysis_date - 6 and no week can straddle the boundary — but it is wrong, and it
  # disagreed with apply_nowcast_correction()'s own `+ 3` by half a day.
  is_future <- week_start_dates > analysis_date

  # ---- flag very uncertain corrections (lag < 3) ----------------------------
  very_uncertain <- !is_future & lag_days < 3
  if (any(very_uncertain)) {
    n_uncertain <- sum(very_uncertain)
    warning(
      "[compute_truncation_weights] ", n_uncertain,
      " week(s) have lag < 3 days: truncation correction is very uncertain ",
      "(weight may be far below 0.5). Proceed with caution.",
      call. = FALSE
    )
  }

  # ---- weight = mean over the bucket's onset days of P(delay <= lag) --------
  # P(Delta <= L) = F(L + 0.5) under the daily-rounding convention the delay was fit
  # under; onset days beyond analysis_date (L < 0) contribute 0 but remain in the
  # denominator (they belong to the week whose TOTAL is being estimated).
  weights  <- rep(NA_real_, length(week_start_dates))
  observed <- !is_future
  if (any(observed)) {
    day_off <- seq.int(0L, as.integer(week_days) - 1L)
    # rows = weeks, cols = the bucket's onset days; lag_mat[i, j] = analysis_date - (w_i + j)
    lag_mat <- outer(as.numeric(analysis_date - week_start_dates[observed]), day_off, "-")
    lag_vec <- as.numeric(lag_mat)                     # column-major, dims dropped
    cdf_vec <- .nowcast_cdf(delay, lag_vec + 0.5)      # pgamma()/plnorm() drop dims too
    cdf_vec[lag_vec < 0] <- 0                          # onset day not yet observable
    weights[observed] <- rowMeans(matrix(cdf_vec, nrow = nrow(lag_mat)))
  }

  # Weights must be in (0, 1]; clamp tiny numerical underflow at a floor of
  # 1e-6 to prevent division by zero (should never occur for positive lags).
  weights[observed & weights < 1e-6] <- 1e-6
  # Clamp to ceiling of 1 (numerical guard only)
  weights[observed & weights > 1]    <- 1.0

  names(weights) <- as.character(week_start_dates)
  # Carry the RESOLVED delay (and the bucket length it was averaged over) out with the
  # weights. apply_nowcast_correction() derives its multiplier cap from the delay, and
  # re-resolving it there would mean two copies of the `rate`/`delay`/NULL precedence rules
  # drifting apart — which is precisely how the original wrong-delay bug survived: a header
  # describing a substitution the code never made. One resolution, one object, no second rule.
  attr(weights, "delay")     <- delay
  attr(weights, "week_days") <- week_days

  if (any(observed)) {
    .tm <- suppressWarnings(as.numeric(delay$trunc_max %||% attr(delay, "trunc_max")))
    message(
      "[compute_truncation_weights] Weights for ",
      sum(observed), " observed weeks; range [",
      round(min(weights[observed], na.rm = TRUE), 4), ", ",
      round(max(weights[observed], na.rm = TRUE), 4), "]",
      "; delay = ", delay$family,
      sprintf(" (mean %.2f d, %s%s)", delay$mean, delay$estimator,
              if (length(.tm) == 1L && is.finite(.tm) && .tm > 0)
                sprintf("; renormalised on [0, %g] d", .tm) else "")
    )
  }
  if (any(is_future)) {
    message(
      "[compute_truncation_weights] ", sum(is_future),
      " future week(s) set to NA (starting after the analysis date)"
    )
  }

  weights
}


# =============================================================================
# SECTION 2: apply_nowcast_correction()
# =============================================================================
#
# Applies truncation weights to the zone × week observation matrix.
#
# Arguments:
#   zone_week_mat  — tibble from aggregate_to_zone_week() with columns:
#                    health_zone, week_start, confirmed, suspected
#   analysis_date  — Date scalar (default: ANALYSIS_DATE from config)
#   min_lag_days   — numeric; weeks with mean lag < min_lag_days are set to NA
#                    rather than corrected (default: 0 — so the most recent
#                    partial week is nowcast-CORRECTED, not suppressed; see the
#                    ROLE note in the function body for why 0 is load-bearing)
#   delay / rate   — passed through to compute_truncation_weights(). The pipeline always
#                    passes the fitted truncation for its regime; leaving both NULL falls
#                    back to the onset->sample delay and warns (see assumption 1)
#
# Returns: augmented tibble with additional columns:
#   trunc_weight         — the fitted-delay completeness weight for the week
#   confirmed_nc         — nowcast-corrected confirmed count
#   suspected_nc         — nowcast-corrected suspected count
#   correction_uncertainty — logical; TRUE if weight < 0.5 (very uncertain)
#
# The multiplicative correction FACTOR (1 / trunc_weight) is capped at 1.5x the deepest
# multiplier the delay can legitimately produce (at least 5x), so the cap catches only a
# pathological weight; it is enabled only once at least one "stable" week (weight >= 0.9)
# exists. NB: the cap is on the factor, not on the absolute count — capping the product at 5×
# the max confirmed count (an earlier behaviour) could push suspected_nc below its raw value.
#
# Assertion checks (stop on failure):
#   (a) All valid corrected counts >= their raw counts
#   (b) All valid weights are in (0, 1]

apply_nowcast_correction <- function(
    zone_week_mat,
    analysis_date = ANALYSIS_DATE,
    min_lag_days  = 0,
    delay         = NULL,
    rate          = NULL
) {
  # ROLE (revised 2026-09-17): this is now the nowcast for BOTH the deployed training data
  # and every leave-future-out fold. It used to be the fold-only estimator, with epinowcast
  # (04b) on the deployed side — but the recalibration delta is fitted on the folds and
  # applied to the deployed forecast, so two estimators meant delta was estimated in one
  # regime and applied in another. Standardising here is what makes that transfer valid.
  # epinowcast is retained as the sensitivity arm (bayes_nowcast_sensitivity).
  #
  # It applies inverse-probability-of-reporting weights derived from the truncation the CALLER
  # passes: the fitted "extract" truncation for the deployed series and for a snapshot round
  # (fitted on that snapshot's date), the fitted "asof" truncation for a re-censored round or
  # calibration origin (epinow2_truncation(), 02_epi_params.R section 4b; see assumption 1 in
  # the file header for why those must differ). Leaving `delay`
  # NULL falls back to the onset->sample delay, which is the wrong leg here, and warns.
  #   HISTORY (do not re-introduce): this argument used to default to the FIXED lab reference
  #   `rate = DELAY_ONSET_SAMPLE_RATE` (0.228/d, mean 4.39 d) while the comment claimed the
  #   data-estimated rate was substituted when ONSET_SAMPLE_DELAY_SOURCE="data". Nothing ever
  #   substituted it and no caller passed `rate`, so every fold was corrected with a delay
  #   about 43% shorter than the one then fitted. Holding the lag convention fixed, that is a strictly
  #   ONE-DIRECTIONAL error: a too-short delay overstates completeness and under-corrects.
  #   Since 2026-09-19 the LFO fold origin is cut + 6, the last day of the last training week,
  #   which is the SAME geometry as the deployed grid, so a fold's most recent week is exactly as
  #   recent as the deployed current week. (The multiplier each receives depends on its
  #   regime's truncation: a snapshot round and the live forecast use the extract truncation,
  #   a re-censored round the as-of one.) Both `delay` and
  #   `rate` are now NULL by default and resolve to the fitted delay; pass `rate` only to
  #   deliberately pin a plain Exp(rate).
  # min_lag_days=0 means the most recent (partial) onset week is nowcast-CORRECTED (upweighted
  # by 1/P(observed), with the MULTIPLIER capped at a delay-derived bound) rather than set to NA
  # — suppressing it zeroed the single week every model forecasts FROM. Two kinds of week get
  # confirmed_nc = NA: a week starting after the nowcast date (nothing observable; trunc_weight
  # NA too), and a week whose MEAN day is after it (lag_days < min_lag_days, i.e. starting 0-2
  # days before the nowcast date at the default 0; trunc_weight is still reported there). The
  # second is counted in the "suppressed" message below. No grid the pipeline builds reaches
  # either: its last week starts 4-6 days before the nowcast date.

  # ---- input validation ------------------------------------------------------
  if (!is.data.frame(zone_week_mat)) {
    stop("[apply_nowcast_correction] zone_week_mat must be a data frame.", call. = FALSE)
  }
  required_cols <- c("health_zone", "week_start", "confirmed", "suspected")
  missing <- setdiff(required_cols, names(zone_week_mat))
  if (length(missing) > 0) {
    stop(
      "[apply_nowcast_correction] Missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!inherits(zone_week_mat$week_start, "Date")) {
    stop("[apply_nowcast_correction] week_start must be a Date column.", call. = FALSE)
  }
  if (!inherits(analysis_date, "Date") || length(analysis_date) != 1) {
    stop("[apply_nowcast_correction] analysis_date must be a scalar Date.", call. = FALSE)
  }
  if (!is.numeric(min_lag_days) || length(min_lag_days) != 1 || min_lag_days < 0) {
    stop("[apply_nowcast_correction] min_lag_days must be a non-negative scalar.", call. = FALSE)
  }

  # ---- idempotency: drop any pre-existing nowcast output columns so a second
  #      pass (e.g. per-fold in LFO-CV) recomputes cleanly rather than colliding
  #      with the left_join below (which would produce trunc_weight.x/.y) -------
  nc_output_cols <- c("trunc_weight", "confirmed_nc", "suspected_nc",
                      "correction_uncertainty", "nowcast_cv")
  pre_existing <- intersect(nc_output_cols, names(zone_week_mat))
  if (length(pre_existing) > 0) {
    zone_week_mat <- zone_week_mat %>% dplyr::select(-dplyr::all_of(pre_existing))
  }

  # ---- compute truncation weights for all unique weeks ----------------------
  all_weeks <- sort(unique(zone_week_mat$week_start))
  weights_vec <- compute_truncation_weights(
    week_start_dates = all_weeks,
    analysis_date    = analysis_date,
    delay            = delay,
    rate             = rate
  )

  # ---- build a lookup tibble -------------------------------------------------
  weight_df <- tibble::tibble(
    week_start   = all_weeks,
    trunc_weight = weights_vec
  )

  # ---- join weights into zone_week_mat ----------------------------------------
  out <- zone_week_mat %>%
    dplyr::left_join(weight_df, by = "week_start")

  # ---- compute lag for each row (used for min_lag_days check) ----------------
  # MEAN lag over the bucket's seven integer onset days (week_start + 3), from the same
  # .week_lag_days() compute_truncation_weights() uses for its "very uncertain" warning (this
  # was a hard-coded `+ 3`). NOT its future-week guard, which is on the week's FIRST day
  # (week_start > analysis_date): a week starting 0-2 days before the nowcast date has a weight
  # but a negative mean lag, and is suppressed here even at the default min_lag_days = 0.
  lag_col <- .week_lag_days(out$week_start, analysis_date)

  # ---- determine validity mask:
  #  - weight must be non-NA (i.e. not a future week)
  #  - lag must be >= min_lag_days
  # The count is every row the lag gate removes from a week that HAS a weight, whatever the sign
  # of its lag: with the default min_lag_days = 0 every such row has a negative lag, and an
  # extra `lag_col > 0` condition (as this used to have) made the count identically zero.
  valid <- !is.na(out$trunc_weight) & lag_col >= min_lag_days
  n_suppressed <- sum(!is.na(out$trunc_weight) & lag_col < min_lag_days)
  if (n_suppressed > 0) {
    message(
      "[apply_nowcast_correction] ", n_suppressed,
      " zone-week rows suppressed (mean-day lag < ", min_lag_days, " days) → confirmed_nc = NA"
    )
  }

  # ---- cap on the correction MULTIPLIER (derived from the delay; see below) --
  # Cap the multiplicative correction factor 1/trunc_weight, NOT the corrected
  # count. The old behaviour capped the *product* at 5x the max CONFIRMED count and
  # applied that same absolute cap to the (typically much larger) SUSPECTED series,
  # which could push suspected_nc BELOW its raw value and trip the "corrected >= raw"
  # assertion below, halting the whole pipeline. Capping the factor (always >= 1 for
  # trunc_weight in (0,1]) guarantees corrected >= raw for both series.
  # Cap the multiplier whenever fully-observed (stable, trunc_weight >= 0.9) weeks EXIST —
  # the decision hinges on having a reliable baseline history, NOT on whether those weeks happen
  # to carry nonzero counts (an all-zero-but-present stable history — e.g. very early outbreak —
  # still warrants the cap). Basing it on max(confirmed) > 0 wrongly disabled the cap for a
  # present-but-zero stable history, and max() on an empty set also warned. Disabled only when no
  # stable week exists yet, where any cap threshold would be arbitrary.
  .n_stable <- out %>%
    dplyr::filter(!is.na(trunc_weight) & trunc_weight >= 0.9) %>%
    nrow()

  # THE CAP IS DERIVED FROM THE DELAY, not fixed at 5x.
  #
  # WHAT IT IS FOR. A week's multiplier is 1/trunc_weight, and compute_truncation_weights()
  # floors trunc_weight at 1e-6, so a pathological weight would multiply a count by up to 1e6.
  # The cap catches that. It is NOT a calibration bound: a large multiplier on a recent week is
  # the EXPECTED consequence of the delay, not an anomaly.
  #
  # WHY IT CHANGED (2026-09-21). It was a flat 5x, justified by a measured "deepest reachable
  # non-suppressed multiplier is 4.577x, so the cap does not bind today". That figure was
  # computed under the onset->SAMPLE delay. The truncation estimated from the vintage archive
  # (onset->appearance; 02_epi_params.R section 4b) is longer, and the anchor week's multiplier
  # moves with every refit — the first fit (2026-09-21) asked for 5.92x — so a flat 5x could
  # silently clip precisely the correction the truncation exists to apply, at precisely the
  # week every forecast is issued from, leaving only the warning below as evidence.
  #
  # WHAT IT IS NOW. 1.5x headroom over the deepest multiplier the delay can legitimately
  # produce, floored at the historic 5x so a very short delay cannot make the guard absurdly
  # tight. Since every admissible week's weight is at least the deepest admissible week's, the
  # cap CANNOT bind on correct input — which is the intent: it fires only on a pathological
  # weight, never on a large-but-correct correction.
  #
  # WHICH WEEK IS "DEEPEST". A week must survive BOTH suppressions:
  #     compute_truncation_weights():  week_start <= analysis_date   (a later start is future)
  #     apply_nowcast_correction():    lag_days >= min_lag_days
  # where lag_days = analysis_date - (week_start + (week_days - 1)/2). The second binds first, so
  # the deepest admissible week has lag_days = min_lag_days and reaches day-lag
  # min_lag_days + (week_days - 1)/2 (capped at the bucket length); its later days have negative
  # lag and contribute 0 to its weight.
  #   With week_days 7 and min_lag_days 0 that is the week starting analysis_date - 3: it
  #   contributes F(0.5)..F(3.5) over a denominator of 7. (Until 2026-10-15 this said
  #   analysis_date - 4, from an older rule under which a lag-0 week was suppressed as future;
  #   the cap was then sized one day too shallow. It never bound, because every grid the
  #   pipeline builds starts its last week 4-6 days before the nowcast date: the week ends ON
  #   the nowcast date for the live and re-censored grids, and 0-2 days after it for a
  #   snapshot round, which is nowcast at its snapshot's date.)
  # Deriving the bound from the most recent FULL week instead would under-bound it: a
  # week-aligned analysis_date only ever reaches analysis_date - 6, but an unaligned origin
  # reaches deeper, and clipping would appear there.
  #
  # `delay` in THIS frame may still be NULL (the caller passed `rate`, or neither), so the cap
  # reads the delay compute_truncation_weights() actually resolved and used, carried out on the
  # weights as an attribute. One resolution, one object.
  .wd  <- attr(weights_vec, "week_days"); if (is.null(.wd) || !length(.wd)) .wd <- 7L
  .dly <- attr(weights_vec, "delay")
  .kmax <- max(0L, min(as.integer(.wd) - 1L,
                       as.integer(floor(min_lag_days + (.wd - 1) / 2))))
  .implied_max_mult <- tryCatch({
    if (is.null(.dly)) NA_real_ else {
      w <- sum(.nowcast_cdf(.dly, seq(0L, .kmax) + 0.5), na.rm = TRUE) / .wd
      if (is.finite(w) && w > 0) 1 / w else NA_real_
    }
  }, error = function(e) NA_real_)
  factor_cap <- if (.n_stable > 0) {
    if (is.finite(.implied_max_mult)) max(5, ceiling(1.5 * .implied_max_mult)) else 5
  } else Inf
  message(sprintf("[apply_nowcast_correction] Correction multiplier cap: %s%s",
                  if (is.finite(factor_cap)) paste0(factor_cap, "x") else "disabled (no stable weeks yet)",
                  if (is.finite(.implied_max_mult))
                    sprintf(" (delay implies at most %.2fx for the deepest admissible week)",
                            .implied_max_mult)
                  else ""))

  # REPORT WHEN THE CAP BINDS. `trunc_weight` is the UNCAPPED weight, so wherever the cap
  # binds the applied multiplier is `factor_cap` while trunc_weight still implies 1/w — and
  # trunc_weight is what run_all.R prints as the nowcast summary and what
  # beta_weighting = "completeness" consumes. The cap is now derived from the delay with 1.5x
  # headroom (above), so on a well-fitted truncation it should not bind at all; if it does,
  # that is a signal the weight is pathological rather than merely small, and it is reported
  # here. The epinowcast backend already reports its own cap hits.
  # Only weeks the lag gate let through (`valid`): a suppressed week gets no correction at all,
  # so a cap cannot have bound on it.
  .binds <- out[valid, , drop = FALSE] %>%
    dplyr::filter(!is.na(trunc_weight), trunc_weight > 0, 1 / trunc_weight > factor_cap) %>%
    dplyr::distinct(week_start)
  if (nrow(.binds))
    warning(sprintf(paste0("[apply_nowcast_correction] the %gx multiplier cap BINDS on %d week(s) ",
                           "(%s). trunc_weight there implies a larger correction than was applied, ",
                           "so it understates the completeness actually achieved."),
                    factor_cap, nrow(.binds),
                    paste(format(sort(.binds$week_start)), collapse = ", ")),
            call. = FALSE, immediate. = TRUE)

  # ---- apply correction -------------------------------------------------------
  out <- out %>%
    dplyr::mutate(
      .cfac = pmin(1 / trunc_weight, factor_cap),
      confirmed_nc = dplyr::case_when(
        !valid                         ~ NA_real_,
        trunc_weight <= 0              ~ NA_real_,
        TRUE                           ~ confirmed * .cfac
      ),
      suspected_nc = dplyr::case_when(
        !valid                         ~ NA_real_,
        trunc_weight <= 0              ~ NA_real_,
        TRUE                           ~ suspected * .cfac
      ),
      correction_uncertainty = !is.na(trunc_weight) & trunc_weight < 0.5,
      # Schema parity with the epinowcast backend, which carries a real CV. The deterministic
      # correction is a point calculation with no uncertainty model behind it, so NA is the
      # honest value — consumers must treat NA as "unknown", never as zero.
      nowcast_cv = NA_real_
    ) %>%
    dplyr::select(-.cfac)

  # Schema parity with the epinowcast backend, which tags itself the same way, so a consumer
  # can always tell which estimator produced the counts it is holding.
  attr(out, "nowcast_method") <- "deterministic"

  # ---- assertions ------------------------------------------------------------
  # (a) valid corrected counts must be >= raw counts
  bad_confirmed <- with(out, !is.na(confirmed_nc) & confirmed_nc < confirmed - 1e-9)
  if (any(bad_confirmed, na.rm = TRUE)) {
    n_bad <- sum(bad_confirmed, na.rm = TRUE)
    stop(
      "[apply_nowcast_correction] Assertion failed: ",
      n_bad, " rows have confirmed_nc < confirmed. ",
      "Check weight computation.",
      call. = FALSE
    )
  }

  bad_suspected <- with(out, !is.na(suspected_nc) & suspected_nc < suspected - 1e-9)
  if (any(bad_suspected, na.rm = TRUE)) {
    n_bad <- sum(bad_suspected, na.rm = TRUE)
    stop(
      "[apply_nowcast_correction] Assertion failed: ",
      n_bad, " rows have suspected_nc < suspected. ",
      "Check weight computation.",
      call. = FALSE
    )
  }

  # (b) valid weights are in (0, 1]
  w_valid <- out$trunc_weight[!is.na(out$trunc_weight)]
  if (any(w_valid <= 0 | w_valid > 1 + 1e-9)) {
    n_bad <- sum(w_valid <= 0 | w_valid > 1 + 1e-9)
    stop(
      "[apply_nowcast_correction] Assertion failed: ",
      n_bad, " weight(s) outside (0, 1]. ",
      "This indicates a bug in compute_truncation_weights().",
      call. = FALSE
    )
  }

  message(
    "[apply_nowcast_correction] Done. ",
    sum(!is.na(out$confirmed_nc)), " zone-week rows have valid confirmed_nc; ",
    sum(out$correction_uncertainty, na.rm = TRUE), " flagged as uncertain (weight < 0.5)"
  )

  out
}



# REMOVED 2026-09-17: compute_effective_observations(), plot_nowcast_qa() and the
# interactive self-test block. None had a caller outside this file, and the self-test
# only exercised the two removed helpers. compute_truncation_weights() and
# apply_nowcast_correction() are the live API; both are covered by tests/.
