# =============================================================================
# 22_daily_reissue.R — Daily re-issue of invasion forecasts (fixed weekly targets)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# WHY THIS FILE EXISTS
# --------------------
# The invasion suite forecasts a FIXED WEEKLY target grid: for a training cutoff
# week C (the current, possibly-partial ISO week), horizon h scores the ISO week
# starting C + 7h. Those weekly targets never move. What DOES move, when the
# pipeline is re-run each day on a fresh extract, is the information state:
#   * epinowcast is refit against the new ANALYSIS_DATE, so the partial current
#     week's expected counts sharpen day by day (essential — a frozen nowcast
#     would make a daily re-run reprint identical numbers);
#   * mid-week invasions retire zones from the at-risk set;
#   * the mobility import force Lambda_i is recomputed on the sharpened counts.
# So the SAME weekly target is re-forecast every day with better information, and
# the invasion probabilities legitimately move between issues.
#
# WHAT THIS MODULE PROVIDES (the modelling machinery is elsewhere — 16_invasion_eval.R and
# 21_bayesian_renewal.R):
#   * anchor_windows_from_analysis_date() / issue_daily_reissue(): re-anchor the weekly
#     forecasts to rolling windows measured FROM the issue date — P(first invasion within
#     7 / 14 days of the analysis date) — and report them. This is the operational product.
#   * reaggregate_asof(): rebuild zone-week counts from the line list censored to a given
#     observation date. This is the leakage-free training reconstruction that the weekly LFO
#     (16), the fold traces (20, 21) and the cascade calibration (33b) all depend on.
#   * linelist_observation_date(): the single definition of when a case became observable.
#
# REMOVED 2026-09-19: the daily BACKTEST (run_invasion_lfo_daily) and the cumulative series
# writer (add_forecast_dating / persist_daily_reissue). The backtest scored the frequentist
# models, which no longer exist, and issued at offsets where the origin week is still "future"
# to the truncation model, so it did not mirror the live operating point it claimed to. The
# weekly LFO now does: its fold origin is cut + 6, the same geometry as the live grid, so it
# exercises the identical nowcast regime.
#
# Window convention: the target window is (forecast_date, forecast_date + window_days], i.e. it
# opens the DAY AFTER the issue date and closes on `window_end`. `window_days` is its LENGTH
# (7 or 14), not a lead time — the lead to the window's start is always exactly 1 day. The field
# was called `lead_days` and documented as "days to the START of the target window" while
# carrying the window LENGTH, so the name, the docstring and the value disagreed three ways.
# =============================================================================


# ---------------------------------------------------------------------------
# 1b. Re-anchor weekly forecasts to the ANALYSIS DATE (rolling 7/14-day windows)
# ---------------------------------------------------------------------------

#' Convert weekly (cutoff-anchored) invasion forecasts into cumulative windows
#' measured from the analysis date: P(first invasion within 7 days) and within
#' 14 days of `forecast_date`.
#'
#' The renewal engine runs weekly; we interpolate its weekly expected introductions
#' to a UNIFORM daily hazard (rate = weekly mu / 7) and integrate over the day
#' windows. With offset a = forecast_date - training_cutoff (0..6, i.e. how far into
#' the current ISO week the analysis date sits), the window (d, d+7] covers the
#' (6-a) remaining days of the current week W plus the first (a+1) days of week W+1;
#' (d, d+14] additionally covers all of W+1 and the first (a+1) days of W+2. So:
#'   mu_7  = mu_wk0·(6-a)/7 + incr1·(a+1)/7
#'   mu_14 = mu_wk0·(6-a)/7 + incr1 + incr2·(a+1)/7
#' where mu_wk0 = current-week (W) expected introductions (from the workhorse; for
#' comparators lacking it, continuity fallback mu_wk0 = incr1), incr1 = week-W+1
#' introductions (= the cumulative h=1 mu), incr2 = week-W+2 introductions
#' (cumulative h=2 − h=1). p = 1 − exp(−mu_window) (Poisson invasion probability,
#' the ascertainment-agnostic primary score). Already-affected zones stay NA.
#'
#' Approximation: the uniform-within-week hazard treats introductions as evenly
#' spread across each ISO week; the partial current week's remaining-days hazard is
#' its whole-week rate scaled by the remaining fraction (unconditioned on the
#' partial non-observation). Documented, and adequate at these horizons.
#'
#' @param fc long weekly forecast tibble: method, health_zone, horizon (1,2),
#'   mu_forecast (CUMULATIVE expected introductions), was_active_before, optionally
#'   mu_wk0 / mobility_id / gt_profile.
#' @param forecast_date issue date (= ANALYSIS_DATE live, or the issue date in the
#'   backtest).
#' @param training_cutoff WEEK_ANCHOR-anchored start of the current week (the model's cutoff).
#' @return long tibble with horizon ∈ {1,2} now meaning "within 7d / within 14d of
#'   the analysis date", window_days, forecast_date, window_end, lead_days,
#'   mu_forecast (windowed), p_invasion/p_case_invasion, was_active_before.
anchor_windows_from_analysis_date <- function(fc, forecast_date, training_cutoff) {
  need <- c("method", "health_zone", "horizon", "mu_forecast")
  stopifnot(is.data.frame(fc), all(need %in% names(fc)))
  # WINDOW GEOMETRY. `a` is how far the issue date sits into its own week. The weekly grid is
  # anchored so the current week ENDS on the analysis date (WEEK_ANCHOR, 00_config.R), so on
  # the deployed path a == 6: none of the current week is still ahead of the issue date, and
  # the 7-day window is exactly the next weekly bucket. The general-a interpolation below is
  # retained so a caller issuing mid-week still gets a sensible answer, but note it can only
  # be exact when the forecast carries `mu_wk0` (the hazard for the REMAINDER of the issue
  # week). No model in the suite supplies that column — the frequentist workhorse was its only
  # producer and it has been removed — so for a < 6 the `wk0` fallback below makes the 7-day
  # window equal to the next full weekly bucket rather than a true a-shifted window. That is
  # exact at a == 6 (the deployed case) and an approximation otherwise; it is flagged rather
  # than hidden, and the warning fires if a caller ever relies on it.
  # OUT-OF-RANGE IS AN ERROR, NOT A CLAMP. `a` is the issue date's position within its own
  # weekly bucket, so it is 0..6 by construction WHEN the grid is anchored to the analysis date.
  # If the grid has gone stale (forecast_date more than 6 days past the last training week, or
  # before it), min()/max() used to pin `a` silently at 6 or 0 — and the function then reported
  # its window as "the next 7 days" when part of that window was already in the PAST. A daily
  # operational product must not mis-date its own window without saying so.
  .a_raw <- as.numeric(forecast_date - training_cutoff)
  if (!is.finite(.a_raw) || .a_raw < 0 || .a_raw > 6)
    stop(sprintf(paste0("[daily_reissue] the issue date (%s) is %.0f day(s) from the last ",
                        "training week's start (%s); it must be 0-6. The weekly grid is stale ",
                        "or the issue date is out of range — rebuild the grid for this date ",
                        "rather than forecasting a window that is partly in the past."),
                 format(forecast_date), .a_raw, format(training_cutoff)), call. = FALSE)
  a <- .a_raw
  if (a < 6)
    warning(sprintf(paste0("[daily_reissue] issue date sits %d day(s) into its week; no model ",
                           "supplies mu_wk0, so the %d-day windows are approximated by the next ",
                           "weekly bucket. Exact only when the issue date ends its week."),
                    as.integer(a), 7L), call. = FALSE)
  fr_wk0 <- (6 - a) / 7      # fraction of the current week still ahead of the analysis date
  fr_hd  <- (a + 1) / 7      # fraction of the far boundary week inside the window

  wide <- fc %>%
    dplyr::filter(.data$horizon %in% c(1L, 2L)) %>%
    dplyr::select(.data$method, .data$health_zone, .data$horizon, .data$mu_forecast) %>%
    # Guard against duplicate (method, health_zone, horizon) keys: without this a caller that
    # passes overlapping frames would make pivot_wider build list-columns and the downstream
    # mu arithmetic would throw "non-numeric argument to binary operator".
    dplyr::distinct(.data$method, .data$health_zone, .data$horizon, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = "horizon", values_from = "mu_forecast",
                       names_prefix = "mu_h")
  if (!"mu_h1" %in% names(wide)) wide$mu_h1 <- NA_real_
  if (!"mu_h2" %in% names(wide)) wide$mu_h2 <- NA_real_

  meta <- fc %>% dplyr::filter(.data$horizon == 1L) %>%
    dplyr::select(dplyr::any_of(c("method", "health_zone", "was_active_before",
                                  "mu_wk0", "mobility_id", "gt_profile"))) %>%
    dplyr::distinct(.data$method, .data$health_zone, .keep_all = TRUE)

  z <- wide %>% dplyr::left_join(meta, by = c("method", "health_zone"))
  if (!"mu_wk0" %in% names(z)) z$mu_wk0 <- NA_real_
  if (!"was_active_before" %in% names(z)) z$was_active_before <- FALSE
  z <- z %>% dplyr::mutate(
    incr1 = .data$mu_h1,
    incr2 = pmax(.data$mu_h2 - .data$mu_h1, 0),
    wk0   = dplyr::coalesce(.data$mu_wk0, .data$incr1),   # comparators: continuity fallback
    mu_7  = .data$wk0 * fr_wk0 + .data$incr1 * fr_hd,
    mu_14 = .data$wk0 * fr_wk0 + .data$incr1 + .data$incr2 * fr_hd)

  build_win <- function(muv, win) {
    muv <- ifelse(is.nan(muv), NA_real_, muv)       # affected/degenerate -> NA, not NaN
    p   <- ifelse(is.na(muv), NA_real_, 1 - exp(-muv))
    tibble::tibble(
      method            = z$method,
      health_zone       = z$health_zone,
      horizon           = as.integer(win %/% 7L),   # 1 = within 7d, 2 = within 14d
      window_days       = as.integer(win),
      forecast_date     = forecast_date,
      window_end        = forecast_date + win,
      window_days       = as.integer(win),   # LENGTH of the window, not a lead time (see header)
      mu_forecast       = muv,
      p_invasion        = p,
      p_case_invasion   = p,
      was_active_before = z$was_active_before,
      mobility_id       = if ("mobility_id" %in% names(z)) z$mobility_id else NA_character_,
      gt_profile        = if ("gt_profile" %in% names(z)) z$gt_profile else NA_character_)
  }
  dplyr::bind_rows(build_win(z$mu_7, 7L), build_win(z$mu_14, 14L))
}

# ---------------------------------------------------------------------------
# 2. Persistence: accumulate daily issues into a series (idempotent per day)
# ---------------------------------------------------------------------------

#' Convenience wrapper for the live run: date-stamp and persist in one call.
#'
#' @param fc_all_current the masked current-week forecast bind from run_all.R.
#' @param forecast_date  Date (= ANALYSIS_DATE).
#' @param training_cutoff Date (WEEK_ANCHOR-anchored start of the current week).
#' @return the dated forecast tibble (also persisted as a side effect).
issue_daily_reissue <- function(fc_all_current, forecast_date, training_cutoff,
                                out_dir = get0("OUT_FORECASTS")) {
  # Re-anchor the weekly forecasts to rolling windows measured FROM the analysis
  # date: P(first invasion within 7 / 14 days of forecast_date). This is the
  # operational product persisted to the daily series.
  fc_anch <- anchor_windows_from_analysis_date(fc_all_current, forecast_date,
                                               training_cutoff)
  # NO LONGER PERSISTED to a cumulative daily series. persist_daily_reissue() upserted into an
  # ever-growing daily_invasion_series file that (a) no script in the repo ever read, and
  # (b) accumulated rows across PROBABILITY SCALES, code versions and model grids with no
  # column recording which: the shipped copy mixed nine raw-scale July issues with one
  # recalibrated September issue, 22 methods from retired modules, and four dating columns
  # that were NA on every row. A product nothing reads and that silently mixes scales is worse
  # than no product. The dated forecast is returned to the caller as before.
  message(sprintf("[daily_reissue] Issue %s (cutoff %s) — invasion probability within:",
                  forecast_date, training_cutoff))
  for (win in c(7L, 14L))
    message(sprintf("    %d days of the analysis date (by %s)",
                    win, forecast_date + win))
  invisible(fc_anch)
}

# ---------------------------------------------------------------------------
# 3. Leakage-free re-aggregation of the linelist as of an issue date
# ---------------------------------------------------------------------------

#' Re-aggregate zone-week counts as they were KNOWN on a given issue date.
#'
#' A retrospective daily backtest must not reuse the final weekly counts: those
#' bucket every case that ever occurred in a week, including cases only reported
#' later. To reconstruct the real-time state on issue date `d` we keep only
#' records whose laboratory/sample OBSERVATION date is <= d, then bucket them by
#' ONSET week (`date_index`) exactly as `aggregate_to_zone_week()` does for the
#' live pipeline. Cases with onset in week W but not yet observed on d are absent
#' — which is precisely the right-truncation the nowcast is meant to correct.
#'
#' Observation date = the date the case becomes knowable to surveillance. We
#' prioritise **date_of_sample_collection** because the deterministic nowcast used
#' downstream models the onset->SAMPLE delay: a case is "observed" when its sample
#' is collected, so censoring on sample collection is delay-consistent with the
#' correction applied to these same counts. Falls back to lab_analysis_date then
#' date_of_notification when the sample date is missing. (Gating instead on the
#' later lab-confirmation date would impose a different, longer delay than the
#' nowcast corrects for, so it is deliberately NOT used as the primary gate.)
#'
#' @param ll linelist (output of load_linelist): must carry date_index and at
#'   least one of lab_analysis_date / date_of_sample_collection /
#'   date_of_notification, plus confirmed / suspected classification columns.
#' @param zones canonical zone spine.
#' @param issue_date Date; keep only records observed on or before this date.
#' @param week_spine optional vector of WEEK_ANCHOR-anchored week starts that MUST be present in
#'   the output (zero-filled if unobserved by `issue_date`). This guarantees the
#'   training frame spans exactly through the cutoff week even when nothing has
#'   been observed in it yet — without it, a Monday-morning issue could drop the
#'   cutoff week entirely and misalign every horizon against its target.
#' @return zone-week tibble (same schema as aggregate_to_zone_week output),
#'   censored to information available on `issue_date`.
reaggregate_asof <- function(ll, zones, issue_date, week_spine = NULL) {
  # THE ISSUE DATE MUST END ITS WEEK. aggregate_to_zone_week() pins the terminal week to
  # floor_date(issue_date, WEEK_ANCHOR), which is the week CONTAINING the issue date — and that
  # week can end up to 6 days AFTER it. The caller's `week_start <= cut` filter then keeps that
  # week as though it were complete, so a calendar-partial week would enter the model as a full
  # one. Every current caller passes cut + 6 (congruent by construction), which is why this has
  # never bitten; there was no guard to keep it that way.
  .wa <- get0("WEEK_ANCHOR", ifnotfound = NULL)
  if (!is.null(.wa)) {
    .ends_week <- identical(
      as.Date(lubridate::floor_date(as.Date(issue_date), "week", week_start = .wa)) + 6L,
      as.Date(issue_date))
    if (!.ends_week)
      warning(sprintf(paste0("[reaggregate_asof] issue_date %s does not END its week under ",
                             "WEEK_ANCHOR=%s, so the terminal bucket is calendar-partial and ",
                             "will be treated as a complete week. Pass the last day of the last ",
                             "training week (cut + 6)."),
                      format(as.Date(issue_date)), .wa), call. = FALSE)
  }
  obs_date <- linelist_observation_date(ll, issue_date, caller = "reaggregate_asof")
  ll_asof <- ll[!is.na(obs_date) & obs_date <= issue_date, , drop = FALSE]
  # aggregate_to_zone_week rebuilds week_start from date_index (onset) and returns
  # the zero-filled zone × week grid over the OBSERVED weeks, so downstream code
  # sees the identical schema it gets from the live path.
  # asof = issue_date, NOT the global ANALYSIS_DATE: this is a FOLD re-aggregation, so the
  # week grid must end at this fold's issue week, not at the run's analysis week.
  zw <- aggregate_to_zone_week(ll_asof, zones, asof = issue_date)

  # Force any required weeks (the cutoff week and every week before it) to exist,
  # even if no case has been observed in them yet on this issue date.
  if (!is.null(week_spine)) {
    miss <- setdiff(as.Date(week_spine), unique(zw$week_start))
    if (length(miss) > 0L) {
      pad <- tidyr::expand_grid(health_zone = zones, week_start = as.Date(miss)) %>%
        dplyr::mutate(confirmed = 0L, suspected = 0L, total_alerts = 0L,
                      tests_analyzed = 0, positivity = NA_real_)
      zw <- dplyr::bind_rows(zw, pad) %>%
        dplyr::arrange(.data$health_zone, .data$week_start)
    }
  }
  zw
}

#' Observation date of every line-list record: the date it became knowable to surveillance.
#'
#' The as-of censoring rule of reaggregate_asof(), factored out so the Bayesian suite's
#' as-of R(t) series (bayes_rt_week_draws(), 21_bayesian_renewal.R) censors on EXACTLY the
#' same rule as the weekly training counts it projects. Priority: sample collection, then
#' lab analysis, then notification (see reaggregate_asof() for why).
#'
#' @param ll linelist (output of load_linelist).
#' @param issue_date scalar Date; the plausibility ceiling when ANALYSIS_DATE is unset.
#' @param caller label prefixed to the diagnostics.
#' @return Date vector parallel to the rows of `ll`; NA where no plausible date exists.
linelist_observation_date <- function(ll, issue_date, caller = "reaggregate_asof") {
  stopifnot(is.data.frame(ll))
  if (!inherits(issue_date, "Date") || length(issue_date) != 1L)
    stop(sprintf("[%s] issue_date must be a scalar Date.", caller), call. = FALSE)
  date_cols <- intersect(c("date_of_sample_collection", "lab_analysis_date",
                           "date_of_notification"), names(ll))
  if (length(date_cols) == 0L)
    stop(sprintf("[%s] linelist has no observation-date column ", caller),
         "(lab_analysis_date / date_of_sample_collection / date_of_notification).",
         call. = FALSE)
  # Coalesce available observation dates in priority order, REJECTING implausible values
  # per column rather than trusting whichever column happened to be non-missing.
  # `date_of_notification` is the last-resort source and is NOT cleaned upstream (only the
  # onset field is): it carries data-entry errors spanning 1999 to 2027. Left unguarded
  # those corrupt the as-of censoring in both directions — a 1999 notification makes a
  # record visible in EVERY fold, including folds long before it could have been known
  # (injecting information the leakage reconstruction exists to withhold), and a 2027 one
  # makes it visible in NONE, silently deleting it from every training frame. An
  # observation date must lie between the outbreak start and the snapshot date; anything
  # else falls through to the next source column, and a record with no plausible date at
  # all is excluded and COUNTED rather than dated by a typo.
  # The floor is OUTBREAK_START minus a 30-day buffer, NOT OUTBREAK_START itself: the
  # earliest real sample in this line list (2026-04-28) precedes the nominal outbreak start
  # (2026-04-30), so a hard floor at OUTBREAK_START would reject genuine early samples and
  # push those records onto a later fallback column. The gate's job is to catch dates that
  # are IMPOSSIBLE (wrong year: 1999, 2023, 2024, 2025, 2027), not to adjudicate a date a
  # few days before the nominal start — a borderline-but-possible alert is left alone.
  .lo <- suppressWarnings(as.Date(get0("OUTBREAK_START", ifnotfound = as.Date("2026-04-30")))) - 30L
  .hi <- suppressWarnings(as.Date(get0("ANALYSIS_DATE",  ifnotfound = issue_date)))
  if (length(.hi) != 1L || is.na(.hi)) .hi <- issue_date
  obs_date <- as.Date(rep(NA_real_, nrow(ll)), origin = "1970-01-01")
  for (cc in date_cols) {
    v <- ll[[cc]]
    if (!inherits(v, "Date")) v <- suppressWarnings(as.Date(v))
    v[!is.na(v) & (v < .lo | v > .hi)] <- as.Date(NA)   # implausible -> try the next column
    obs_date <- dplyr::coalesce(obs_date, v)
  }
  .undated <- is.na(obs_date)
  if (any(.undated)) {
    .nc <- if ("confirmed" %in% names(ll)) sum(.undated & ll$confirmed %in% TRUE) else 0L
    msg <- sprintf(paste0("[%s] %d record(s) carry no plausible observation ",
                          "date in [%s, %s] and are excluded from the as-of reconstruction ",
                          "(%d of them CONFIRMED)."),
                   caller, sum(.undated), format(.lo), format(.hi), .nc)
    # Undated CONFIRMED cases would distort the invasion outcome itself, so those are loud.
    if (.nc > 0L) warning(msg, call. = FALSE) else message(msg)
  }
  obs_date
}

# ---------------------------------------------------------------------------
# 4. Daily-issue backtest (mirrors the LIVE operating point)
# ---------------------------------------------------------------------------


message("[daily_reissue] 22_daily_reissue.R loaded — daily re-issue + intra-week backtest.")
