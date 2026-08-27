# =============================================================================
# 04_nowcasting.R — Right-Truncation Correction (Nowcasting)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: Correct the zone × week observation matrix for right-truncation
#   arising because recent onset dates have not yet been sampled/reported.
#   The correction is re-applied at each LFO-CV fold so that no future
#   information leaks into the training window.
#
# Key statistical context:
#   The onset-to-sample delay follows an Exponential distribution with
#   rate = DELAY_ONSET_SAMPLE_RATE (≈ 0.228 day^-1, mean ≈ 4.39 d, fit to
#   BDBV-2026 Ituri lab linelist; reloaded from the CFR CSV at runtime, the source
#   of truth). For a week whose midpoint is
#   `lag` days before the analysis date, the probability that an onset from
#   that week has already been sampled (and thus appears in the linelist) is:
#
#     P(delay <= lag) = 1 - exp(-rate * lag)     [Exponential CDF]
#
#   We call this the truncation weight. The nowcast-corrected count is:
#
#     confirmed_nc = confirmed / truncation_weight
#
#   This is a simple inverse-probability-of-reporting correction.  It is
#   analytic and fast, but ignores delay uncertainty (a Bayesian treatment
#   is implemented in the bayesian/ sub-suite).
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
# Computes a named vector of Exponential-CDF truncation weights for each week.
#
# Arguments:
#   week_start_dates — Date vector; the start of each 7-day bucket (WEEK_ANCHOR-anchored)
#   analysis_date    — Date scalar; reference date for the correction
#                      (typically ANALYSIS_DATE from 00_config.R)
#   rate             — Exponential decay rate (day^-1); default is
#                      DELAY_ONSET_SAMPLE_RATE from 00_config.R
#
# Returns:
#   Named numeric vector; names are character representations of week_start_dates.
#   Values are in (0, 1] for observed weeks, NA for weeks whose midpoint has not
#   yet passed (lag <= 0).
#
# Warnings:
#   - Emits a warning for any week where lag < 3 days (weight is very uncertain
#     and the correction can be extremely large).

compute_truncation_weights <- function(
    week_start_dates,
    analysis_date,
    rate = DELAY_ONSET_SAMPLE_RATE
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
  if (!is.numeric(rate) || length(rate) != 1 || rate <= 0) {
    stop("[compute_truncation_weights] rate must be a single positive number.", call. = FALSE)
  }

  # ---- compute lag from week midpoint to analysis_date ----------------------
  # Week midpoint = bucket start + 3.5 days (centre of the 7-day bucket, anchor-agnostic)
  week_midpoints <- week_start_dates + 3.5
  lag_days <- as.numeric(analysis_date - week_midpoints)

  # ---- flag future weeks (lag <= 0) -----------------------------------------
  # These weeks are entirely in the future relative to analysis_date; no
  # observations are yet possible so the weight is undefined (NA).
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

  # ---- compute Exponential CDF weight  P(delay <= lag) ----------------------
  weights <- rep(NA_real_, length(week_start_dates))
  observed <- !is_future
  weights[observed] <- 1 - exp(-rate * lag_days[observed])

  # Weights must be in (0, 1]; clamp tiny numerical underflow at a floor of
  # 1e-6 to prevent division by zero (should never occur for positive lags).
  weights[observed & weights < 1e-6] <- 1e-6
  # Clamp to ceiling of 1 (numerical guard only)
  weights[observed & weights > 1]    <- 1.0

  names(weights) <- as.character(week_start_dates)

  if (any(observed)) {
    message(
      "[compute_truncation_weights] Weights for ",
      sum(observed), " observed weeks; range [",
      round(min(weights[observed], na.rm = TRUE), 4), ", ",
      round(max(weights[observed], na.rm = TRUE), 4), "]"
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
#   min_lag_days   — numeric; weeks with lag < min_lag_days are set to NA
#                    rather than corrected (default: 0 — so the most recent
#                    partial week is nowcast-CORRECTED, not suppressed; see the
#                    ROLE note in the function body for why 0 is load-bearing)
#
# Returns: augmented tibble with additional columns:
#   trunc_weight         — the Exponential-CDF weight for the week
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
    rate          = DELAY_ONSET_SAMPLE_RATE
) {
  # ROLE: this is the FAST, deterministic per-fold nowcast for leave-future-out CV. The
  # PRIMARY current-week nowcast is the full probabilistic epinowcast (04b_epinowcast.R,
  # nowcast_zone_week_epinowcast); it is not refit inside the ~11 LFO folds because a Stan
  # fit per fold is prohibitively expensive and fails on event-sparse early folds. Instead
  # each fold applies inverse-probability-of-reporting weights from the SAME onset->sample
  # reporting delay epinowcast models (rate = DELAY_ONSET_SAMPLE_RATE, or the data-estimated
  # rate when ONSET_SAMPLE_DELAY_SOURCE="data"), so the two are delay-consistent by construction.
  # min_lag_days=0 means the most recent (partial) onset week is nowcast-CORRECTED (upweighted
  # by 1/P(observed), with the MULTIPLIER capped at 5x) rather than suppressed to NA — suppressing
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
                      "correction_uncertainty")
  pre_existing <- intersect(nc_output_cols, names(zone_week_mat))
  if (length(pre_existing) > 0) {
    zone_week_mat <- zone_week_mat %>% dplyr::select(-dplyr::all_of(pre_existing))
  }

  # ---- compute truncation weights for all unique weeks ----------------------
  all_weeks <- sort(unique(zone_week_mat$week_start))
  weights_vec <- compute_truncation_weights(
    week_start_dates = all_weeks,
    analysis_date    = analysis_date,
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
  lag_col <- as.numeric(analysis_date - (out$week_start + 3.5))

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

  factor_cap <- if (.n_stable > 0) 5 else Inf
  message("[apply_nowcast_correction] Correction multiplier cap: ",
          if (is.finite(factor_cap)) paste0(factor_cap, "x") else "disabled (no stable weeks yet)")

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
      correction_uncertainty = !is.na(trunc_weight) & trunc_weight < 0.5
    ) %>%
    dplyr::select(-.cfac)

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


# =============================================================================
# SECTION 3: compute_effective_observations()
# =============================================================================
#
# For each health zone, determine the first and last weeks that had at least
# one confirmed case, and compute basic activity statistics. Useful for
# reporting which zones are in the "active" epidemic front.
#
# Arguments:
#   zone_week_mat — tibble from apply_nowcast_correction() or
#                   aggregate_to_zone_week(), containing:
#                   health_zone, week_start, confirmed
#
# Returns a tibble with one row per health zone:
#   health_zone, first_confirmed_week, last_confirmed_week,
#   n_active_weeks, total_confirmed, ever_active

compute_effective_observations <- function(zone_week_mat) {

  if (!is.data.frame(zone_week_mat)) {
    stop("[compute_effective_observations] Input must be a data frame.", call. = FALSE)
  }
  missing <- setdiff(c("health_zone", "week_start", "confirmed"), names(zone_week_mat))
  if (length(missing) > 0) {
    stop(
      "[compute_effective_observations] Missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  active_summary <- zone_week_mat %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(
      total_confirmed      = sum(confirmed, na.rm = TRUE),
      n_active_weeks       = sum(confirmed > 0, na.rm = TRUE),
      first_confirmed_week = dplyr::if_else(
        any(confirmed > 0, na.rm = TRUE),
        min(week_start[confirmed > 0], na.rm = TRUE),
        as.Date(NA)
      ),
      last_confirmed_week  = dplyr::if_else(
        any(confirmed > 0, na.rm = TRUE),
        max(week_start[confirmed > 0], na.rm = TRUE),
        as.Date(NA)
      ),
      ever_active = any(confirmed > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(total_confirmed), health_zone)

  n_ever_active <- sum(active_summary$ever_active)
  message(
    "[compute_effective_observations] ",
    n_ever_active, " / ", nrow(active_summary),
    " zones have ever had >=1 confirmed case"
  )

  active_summary
}


# =============================================================================
# SECTION 4: QA plot — raw vs corrected counts for top zones
# =============================================================================
#
# Generates a console-only QA visualisation comparing raw and nowcast-corrected
# confirmed counts for the top 5 most affected zones. Not saved to disk.
# Call manually after apply_nowcast_correction() to inspect corrections visually.
#
# Arguments:
#   zw_corrected — output of apply_nowcast_correction()
#   top_n        — integer; number of top zones to plot (default 5)

plot_nowcast_qa <- function(zw_corrected, top_n = 5L) {

  if (!is.data.frame(zw_corrected)) {
    warning("[plot_nowcast_qa] Input must be a data frame. Skipping plot.", call. = FALSE)
    return(invisible(NULL))
  }
  required <- c("health_zone", "week_start", "confirmed", "confirmed_nc", "trunc_weight")
  missing  <- setdiff(required, names(zw_corrected))
  if (length(missing) > 0) {
    warning(
      "[plot_nowcast_qa] Missing columns for QA plot: ",
      paste(missing, collapse = ", "),
      ". Skipping.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("[plot_nowcast_qa] ggplot2 not available — skipping QA plot.")
    return(invisible(NULL))
  }

  top_zones <- zw_corrected %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(tot = sum(confirmed, na.rm = TRUE), .groups = "drop") %>%
    dplyr::slice_max(tot, n = top_n, with_ties = FALSE) %>%
    dplyr::pull(health_zone)

  if (length(top_zones) == 0) {
    message("[plot_nowcast_qa] No zones with confirmed > 0; skipping QA plot.")
    return(invisible(NULL))
  }

  plot_df <- zw_corrected %>%
    dplyr::filter(health_zone %in% top_zones) %>%
    dplyr::select(health_zone, week_start, confirmed, confirmed_nc, trunc_weight) %>%
    tidyr::pivot_longer(
      cols      = c(confirmed, confirmed_nc),
      names_to  = "series",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      series = dplyr::case_match(series,
                                 "confirmed"    ~ "Raw",
                                 "confirmed_nc" ~ "Nowcast-corrected"),
      health_zone = factor(health_zone, levels = top_zones)
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = week_start, y = count, colour = series, linetype = series)
  ) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.5, na.rm = TRUE) +
    ggplot2::facet_wrap(~ health_zone, scales = "free_y", ncol = 1L) +
    ggplot2::scale_colour_manual(
      values = c("Raw" = "#5b86b3", "Nowcast-corrected" = "#b23b2e")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Raw" = "solid", "Nowcast-corrected" = "dashed")
    ) +
    ggplot2::labs(
      title    = paste0("QA: Raw vs Nowcast-Corrected Confirmed Cases (Top ", top_n, " Zones)"),
      subtitle = paste0("Exponential delay, rate = ", DELAY_ONSET_SAMPLE_RATE,
                        " day⁻¹; analysis_date = ", ANALYSIS_DATE),
      x        = "Week start (7-day bucket ending on the analysis weekday)",
      y        = "Confirmed cases",
      colour   = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_classic(base_size = VIZ_THEME_BASE) +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  message("[plot_nowcast_qa] QA plot rendered for zones: ",
          paste(top_zones, collapse = ", "))

  invisible(p)
}


# =============================================================================
# SECTION 6: Demonstration / self-test (runs when script is sourced directly)
# =============================================================================
#
# When this script is sourced standalone (i.e., `dat` from 01_data_prep.R is
# already in the environment), automatically runs the correction and QA plot.
# If `dat` is not available, emits an informative message and exits cleanly.

if (interactive() && exists("dat") && !is.null(dat$zone_week)) {

  message("\n[04_nowcasting.R] dat$zone_week found — running demonstration correction.")

  zw_corrected <- apply_nowcast_correction(
    zone_week_mat = dat$zone_week,
    analysis_date = ANALYSIS_DATE,
    min_lag_days  = 2L
  )

  eff_obs <- compute_effective_observations(zw_corrected)

  message("\n[04_nowcasting.R] Top 10 zones by total confirmed (effective observations):")
  print(utils::head(eff_obs, 10))

  # QA plot (console only, not saved)
  plot_nowcast_qa(zw_corrected, top_n = 5L)

  message("\n[04_nowcasting.R] Objects created:")
  message("  zw_corrected  — zone × week matrix with nowcast correction")
  message("  eff_obs       — zone-level effective observation summary")

} else {
  message(
    "\n[04_nowcasting.R] Loaded. `dat` not found in environment.",
    "\nRun 01_data_prep.R first, or call functions directly:",
    "\n  compute_truncation_weights(week_start_dates, analysis_date)",
    "\n  apply_nowcast_correction(zone_week_mat, analysis_date)",
    "\n  compute_effective_observations(zone_week_mat)",
    "\n  apply_nowcast_correction(zone_week_mat, analysis_date = fold_date)  [per-fold LFO-CV]"
  )
}
