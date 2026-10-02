# =============================================================================
# 04_nowcasting.R — Right-Truncation Correction (Nowcasting)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: Correct the zone × week observation matrix for right-truncation
#   arising because recent onset dates have not yet been sampled/reported.
#   The correction is re-applied at each LFO-CV fold so that no future
#   information leaks into the training window.
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
#          truncating leg is onset -> APPEARANCE-IN-EXTRACT. Fitted against the vintage archive.
#        - "asof" — each LFO fold (run_invasion_lfo) and each cascade calibration origin
#          (33b_cascade_calibration.R) is rebuilt by reaggregate_asof(), which censors on
#          linelist_observation_date() — the SAMPLE-collection date for 100% of confirmed
#          records. The truncating leg there is onset -> OBSERVATION: shorter, because it does
#          not wait for the laboratory result and data entry. Measured by onset cohort rather
#          than deconvolved, since both dates are recorded (see .trunc_fit_asof_cohort()).
#          Anchor week 3.45x here against 5.92x for the extract regime.
#      Until 2026-09-21 both paths used the onset->SAMPLE delay. For the DEPLOYED series that
#      is the wrong leg, and the error is one-directional: measured against 54 archived
#      vintages (5,607 records matched to their first-appearance vintage) onset->appearance has
#      mean 11.43 d against the modelled 7.67 d, so the most recent week was inflated 2.63x
#      where the fitted truncation gives 5.92x. That week is the one every forecast is issued
#      FROM, so the shortfall propagated into every h=1 invasion probability.
#      THE TWO REGIMES ARE DELIBERATELY DIFFERENT DISTRIBUTIONS. The invariant that matters is
#      not "same multiplier both sides" but "each series corrected for its own truncation, so
#      both estimate the same quantity — true onsets per week". Forcing one distribution on
#      both would necessarily bias one of them.
#      RESIDUAL, AND IT IS ACCEPTED: the as-of reconstruction is OPTIMISTIC. At a real origin T
#      you would not yet have known about a case sampled at T-1 but not yet confirmed and
#      entered; reaggregate_asof() assumes instant lab-and-data-entry and shows it anyway. So
#      the folds see a cleaner world than deployment does, and the recalibration delta fitted
#      on them carries that optimism into the deployed forecast. Correcting each side for its
#      own truncation is the right thing to do GIVEN the reconstruction; it does not repair the
#      reconstruction. Removing that gap needs vintage-censored folds, which is out of scope by
#      explicit decision — optimistic LFO skill is acceptable, a wrong deployed R(t) was not.
#      This file's own default is STILL the onset->sample delay (it is the only delay available
#      without an archive), so leaving `delay` NULL now WARNS once per session.
#
#   2. THE DELAY IS FITTED ON ALL COMPLETE PAIRS, not only confirmed ones (04c's
#      `confirmed_only = FALSE`, ~14k pairs against ~6.7k confirmed cases). This assumes
#      suspected and confirmed cases are sampled on the same schedule. Untested.
#
#   3. THE SUSPECTED SERIES IS CORRECTED WITH THE SAME CDF. A suspected case's appearance is
#      governed by alert notification, not by specimen collection, so `suspected_nc` carries a
#      delay that was not fitted for it. It feeds only the suspected-covariate models, which
#      are OFF by default (INCLUDE_SUSPECTED_COV_MODELS).
#
#   4. NO WEEK IS EVER "FULLY OBSERVED". The gamma CDF never reaches 1, so every historical
#      week keeps a multiplier just above 1 (x1.011 at 41 d, x1.001 at ~90 d) and every
#      confirmed_nc is non-integer. The residual is bounded by 1 - F(60) = 0.001, since 04c
#      discards delays beyond DELAY_MAX_PLAUSIBLE_DAYS = 60 when fitting; the fitted
#      distribution is not renormalised on [0, 60], so that 0.1% is an inflation applied
#      uniformly across the whole series.
#
#   5. THE CORRECTION IS NEVER SCORED. Nothing in the pipeline evaluates the nowcast against
#      later-observed truth: compare_epinowcast_models() exists but has no caller, and
#      Renewal-M8-raw (which compared raw against corrected counts through the LFO) went with
#      the frequentist arm. What IS checked is the GEOMETRY: since 2026-09-19 the LFO fold
#      origin is cut + 6, the last day of the last training week, which is the same geometry as
#      the live grid, so a fold's most recent week occupies the same position relative to its
#      origin as the deployed current week does. The MULTIPLIERS differ (and should — see
#      assumption 1), because the two regimes are truncated by different legs.
#
# Key statistical context:
#   F is the CDF of whatever delay the CALLER passes — in the pipeline, the fitted truncation
#   for that caller's regime (assumption 1). It is evaluated on its OWN family: lnorm for the
#   fitted truncations, gamma for the retired onset->sample fit. When no delay is passed this
#   falls back to effective_onset_sample_delay() (00_config.R) — the EpiDist MARGINAL fit from
#   04c_dhis2_delay_windows.R (Gamma, mean 7.665 d, SD 8.504 d, truncation- and
#   double-interval-censoring corrected), then the windowed interval-censored MLE, then the
#   fixed lab reference Exp(DELAY_ONSET_SAMPLE_RATE). That fallback is the WRONG LEG for a
#   truncation correction and warns.
#
#   For a week w, the truncation weight is the expected fraction of that week's onsets
#   already sampled by the analysis date — averaged over the week's seven onset days,
#   NOT evaluated once at the week midpoint:
#
#     weight_w = (1/7) * sum_{j=0}^{6} P(delay <= analysis_date - (w + j))
#              = (1/7) * sum_{j=0}^{6} F(analysis_date - (w + j) + 0.5)
#
#   (the + 0.5 is the daily-rounding convention the delay was FIT under in 04c).
#   The nowcast-corrected count is:
#
#     confirmed_nc = confirmed / weight_w
#
#   This is a simple inverse-probability-of-reporting correction.  It is
#   analytic and fast, but ignores delay uncertainty (the full probabilistic
#   treatment is 04b_epinowcast.R, used for the primary current-week nowcast).
#
# Dependencies: 00_config.R (sourced first), tidyverse, lubridate
#
# Usage (typical):
#   source("spatiotemporal/04_nowcasting.R")
#   zw_corrected <- apply_nowcast_correction(dat$zone_week, ANALYSIS_DATE)
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
#      4.39 d against a fitted 7.665 d — 43% too short, i.e. reporting assumed far more
#      complete than it is, so the correction was too small. It now defaults to
#      effective_onset_sample_delay() (00_config.R) — the EpiDist marginal fit
#      (truncation- and double-interval-censoring corrected) when 04c wrote one, else the
#      windowed interval-censored MLE, else the lab reference — and evaluates that fit's
#      ACTUAL family (delay_cdf(); gamma here, shape 0.812, rate 0.106) rather than moment-matching it
#      to an Exponential. SUPERSEDED 2026-09-21: that fixed the delay's magnitude but not its
#      LEG. The deployed series is truncated by onset->appearance, not onset->sample, so the
#      callers now pass a fitted truncation and this resolver is only the fallback (assumption
#      1). The onset->sample fit remains the right object for the missing-onset IMPUTATION,
#      which really is imputing to sample; imputation and nowcast are therefore no longer the
#      same distribution, and correctly so.
#      Validated against the untruncated empirical delay (onsets >= 35 d before the as-of
#      date, n=3821, mean 7.21 d): mean |CDF error| over lags 0-30 d is 0.018 for the
#      fitted delay vs 0.087 for the lab reference.
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
#      Holding the delay fixed, (b) dominates (a) on this snapshot: for the current week the
#      fitted delay differs materially between the retired midpoint convention and the correct
#      within-week average, and the convention change partly offsets fix (1). Measured under the
#      retired onset->sample delay (2026-09-19) the deployed current week carried weight 0.3801,
#      a 2.631x correction; under the fitted extract truncation (2026-09-21) the same week
#      carries weight 0.1688, a 5.92x correction. The cap is derived from the delay rather than
#      fixed at 5x precisely so that this remains a correction and not a clip. We average F over the week's seven onset days, with
#      P(Delta <= L) = F(L + 0.5) for integer L (the interval-censoring convention 04c fits
#      under, d -> [d-0.5, d+0.5]). Onset days after `analysis_date` contribute 0 (not yet
#      observable) but stay in the denominator, since they belong to the week whose total is
#      being estimated.
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
#   Values are in (0, 1] for observed weeks, NA for weeks whose mean lag has not yet
#   passed (mean lag <= 0) — the same weeks the previous midpoint rule flagged as future.
#
# Warnings:
#   - Emits a warning for any week where the mean lag < 3 days (weight is very uncertain
#     and the correction can be extremely large).

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
    # delay 44% too short for four months with nothing in the log to say so. Once per session
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
  lag_days <- as.numeric(analysis_date - (week_start_dates + (week_days - 1) / 2))

  # ---- flag future weeks (mean lag <= 0) ------------------------------------
  # These weeks are (essentially) in the future relative to analysis_date; the weight is
  # undefined (NA). For integer dates this selects exactly the same weeks as the previous
  # midpoint rule (analysis_date <= week_start + 3), so the guard is unchanged.
  is_future <- lag_days <= 0

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
    cdf_vec <- delay_cdf(delay, lag_vec + 0.5)         # pgamma()/plnorm() drop dims too
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
    message(
      "[compute_truncation_weights] Weights for ",
      sum(observed), " observed weeks; range [",
      round(min(weights[observed], na.rm = TRUE), 4), ", ",
      round(max(weights[observed], na.rm = TRUE), 4), "]",
      "; delay = ", delay$family,
      sprintf(" (mean %.2f d, %s)", delay$mean, delay$estimator)
    )
  }
  if (any(is_future)) {
    message(
      "[compute_truncation_weights] ", sum(is_future),
      " future week(s) set to NA (lag <= 0)"
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
# The multiplicative correction FACTOR (1 / trunc_weight) is capped at 5× (so a
# corrected count stays in [raw, 5·raw]); the cap is enabled only once at least one
# "stable" week (weight >= 0.9) exists. NB: the cap is on the factor, not on the
# absolute count — capping the product at 5× the max confirmed count (the earlier
# behaviour) could push suspected_nc below its raw value and trip the assertions.
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
  # passes: the fitted "extract" truncation for the deployed series, the fitted "asof"
  # truncation for every fold and calibration origin (epinow2_truncation(), 02_epi_params.R
  # section 4b; see assumption 1 in the file header for why those must differ). Leaving `delay`
  # NULL falls back to the onset->sample delay, which is the wrong leg here, and warns.
  #   HISTORY (do not re-introduce): this argument used to default to the FIXED lab reference
  #   `rate = DELAY_ONSET_SAMPLE_RATE` (0.228/d, mean 4.39 d) while the comment claimed the
  #   data-estimated rate was substituted when ONSET_SAMPLE_DELAY_SOURCE="data". Nothing ever
  #   substituted it and no caller passed `rate`, so every fold was corrected with a delay
  #   44% shorter than the fitted one. Holding the lag convention fixed, that is a strictly
  #   ONE-DIRECTIONAL error: a too-short delay overstates completeness and under-corrects.
  #   Since 2026-09-19 the LFO fold origin is cut + 6, the last day of the last training week,
  #   which is the SAME geometry as the deployed grid — so a fold's most recent week and the
  #   deployed current week carry the identical 2.631x correction, and the recalibration delta
  #   is fitted in the regime it is applied in. Both `delay` and
  #   `rate` are now NULL by default and resolve to the fitted delay; pass `rate` only to
  #   deliberately pin a plain Exp(rate).
  # min_lag_days=0 means the most recent (partial) onset week is nowcast-CORRECTED (upweighted
  # by 1/P(observed), with the MULTIPLIER capped at a delay-derived bound) rather than to NA — suppressing
  # it zeroed the single week every model forecasts FROM. Only genuinely future weeks (lag<0,
  # trunc_weight NA) are dropped.

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
  # MEAN lag over the bucket's seven integer onset days (week_start + 3) — the same
  # definition compute_truncation_weights() uses for its future-week guard, so the two
  # cannot disagree about which weeks are observable. (With the default min_lag_days = 0
  # this gate never bites; it selects the same weeks as the previous + 3.5 form anyway.)
  lag_col <- as.numeric(analysis_date - (out$week_start + 3))

  # ---- determine validity mask:
  #  - weight must be non-NA (i.e. not a future week)
  #  - lag must be >= min_lag_days
  valid <- !is.na(out$trunc_weight) & lag_col >= min_lag_days
  n_suppressed <- sum(!is.na(out$trunc_weight) & lag_col < min_lag_days & lag_col > 0)
  if (n_suppressed > 0) {
    message(
      "[apply_nowcast_correction] ", n_suppressed,
      " zone-week rows suppressed (lag < ", min_lag_days, " days) → confirmed_nc = NA"
    )
  }

  # ---- cap: correction MULTIPLIER capped at 5x (so corrected in [raw, 5*raw]) --
  # Cap the multiplicative correction factor 1/trunc_weight, NOT the corrected
  # count. The old behaviour capped the *product* at 5x the max CONFIRMED count and
  # applied that same absolute cap to the (typically much larger) SUSPECTED series,
  # which could push suspected_nc BELOW its raw value and trip the "corrected >= raw"
  # assertion below, halting the whole pipeline. Capping the factor (always >= 1 for
  # trunc_weight in (0,1]) guarantees corrected >= raw for both series.
  # Cap the multiplier at 5x whenever fully-observed (stable, trunc_weight >= 0.9) weeks EXIST —
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
  # computed under the onset->SAMPLE delay. With the truncation estimated from the vintage
  # archive (onset->appearance; 02_epi_params.R section 4b) the deployed anchor week
  # legitimately needs 5.92x — so a flat 5x would have silently clipped precisely the
  # correction the truncation fix exists to apply, at precisely the week every forecast is
  # issued from, leaving only the warning below as evidence.
  #
  # WHAT IT IS NOW. 1.5x headroom over the deepest multiplier the delay can legitimately
  # produce, floored at the historic 5x so a very short delay cannot make the guard absurdly
  # tight. Since every admissible week's weight is at least the deepest admissible week's, the
  # cap CANNOT bind on correct input — which is the intent: it fires only on a pathological
  # weight, never on a large-but-correct correction.
  #
  # WHICH WEEK IS "DEEPEST". A week must survive BOTH suppressions:
  #     compute_truncation_weights():  lag_days > 0          (strict — lag 0 counts as future)
  #     apply_nowcast_correction():    lag_days >= min_lag_days
  # where lag_days = analysis_date - (week_start + (week_days - 1)/2). For integer dates and an
  # odd bucket length lag_days is an integer, so the smallest admissible lag is
  # max(min_lag_days, 1), and that week reaches day-lag max(min_lag_days, 1) + (week_days-1)/2,
  # capped at the bucket length.
  #   For the live grid (week_days 7, min_lag_days 0) that is the week starting
  #   analysis_date - 4: it contributes F(0.5)..F(4.5) over a denominator of 7. NOT
  #   analysis_date - 3 — that week has lag_days 0 and is suppressed as future.
  # Deriving the bound from the most recent FULL week instead would under-bound it: a
  # week-aligned analysis_date only ever reaches analysis_date - 6, but an unaligned origin
  # reaches deeper, and clipping would have reappeared there. Recomputing the historic 4.577x
  # from this rule under the retired delay gives 4.578x, which is what confirms the rule picks
  # out the right week; tests/test_truncation.R pins that.
  #
  # `delay` in THIS frame may still be NULL (the caller passed `rate`, or neither), so the cap
  # reads the delay compute_truncation_weights() actually resolved and used, carried out on the
  # weights as an attribute. One resolution, one object.
  .wd  <- attr(weights_vec, "week_days"); if (is.null(.wd) || !length(.wd)) .wd <- 7L
  .dly <- attr(weights_vec, "delay")
  .kmax <- max(0L, min(as.integer(.wd) - 1L,
                       as.integer(floor(max(min_lag_days, 1) + (.wd - 1) / 2))))
  .implied_max_mult <- tryCatch({
    if (is.null(.dly)) NA_real_ else {
      w <- sum(delay_cdf(.dly, seq(0L, .kmax) + 0.5), na.rm = TRUE) / .wd
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
  .binds <- out %>%
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
