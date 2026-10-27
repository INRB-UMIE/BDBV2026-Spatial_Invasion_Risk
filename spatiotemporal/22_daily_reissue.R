# =============================================================================
# 22_daily_reissue.R — As-of reconstruction, archived-snapshot dispatch, and issue-date windows
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# WHAT THIS MODULE PROVIDES (the modelling machinery is elsewhere — 16_invasion_eval.R and
# 21_bayesian_renewal.R):
#   * reaggregate_asof(): rebuild zone-week counts from a line list censored to a given
#     observation date. This is the "what was known then" reconstruction for a back-test round
#     that cannot read an archived snapshot; the weekly LFO (16), the fold traces (20, 21) and
#     the cascade calibration (33b) reach it through asof_counts_source().
#   * linelist_observation_date(): the single definition of when a case became observable.
#   * vintage_for_issue(), vintage_truncation(), prefetch_vintages(), asof_counts_source()
#     (section 3b): the archived-snapshot information state, and the ONE dispatcher that decides,
#     per round, between a snapshot (extract-regime truncation fitted on its date) and a
#     re-censored final line list (shared "asof" truncation).
#   * anchor_windows_from_analysis_date() / issue_daily_reissue(): re-express the weekly
#     forecasts as windows measured FROM the issue date — P(first invasion within 7 / 14 days).
#     run_all.R calls issue_daily_reissue() once per run, which reports the window dates and
#     returns the windows; nothing is persisted and no other module reads them.
#
# THE WEEKLY TARGETS. Weeks are WEEK_ANCHOR weeks (00_config.R), labelled by their first day and
# anchored on the analysis weekday, so the last training week ENDS on the analysis date. For a
# cutoff C (the first day of the last training week) horizon h scores the CUMULATIVE window of
# weeks C+7 .. C+7h, i.e. days (C+6, C+6+7h]. A run on a later extract anchors its own grid on
# its own analysis date and re-forecasts from there: the deterministic nowcast (04, with the
# extract truncation) is recomputed at the new date, invasions observed since retire zones from
# the at-risk set, and the mobility import force Lambda_i is recomputed on the new counts.
#
# REMOVED 2026-09-19: the daily BACKTEST (run_invasion_lfo_daily) and the cumulative series
# writer (add_forecast_dating / persist_daily_reissue). The backtest scored the frequentist
# models, which no longer exist, and issued at offsets where the origin week is still "future"
# to the truncation model, so it did not mirror the live operating point it claimed to. The
# weekly LFO now does: its fold origin is cut + 6, the same geometry as the live grid. Its
# nowcast REGIME matches the live one (extract truncation) on the rounds that read an archived
# snapshot with a real-time truncation fit; the rounds re-censored from the final line list use
# the as-of truncation instead (section 3b; 04_nowcasting.R, assumption 1).
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
#' its own WEEK_ANCHOR week the analysis date sits), the window (d, d+7] covers the
#' (6-a) remaining days of the current week W plus the first (a+1) days of week W+1;
#' (d, d+14] additionally covers all of W+1 and the first (a+1) days of W+2. So:
#'   mu_7  = mu_wk0·(6-a)/7 + incr1·(a+1)/7
#'   mu_14 = mu_wk0·(6-a)/7 + incr1 + incr2·(a+1)/7
#' where mu_wk0 = the expected introductions over the REMAINDER of week W (no model in the
#' suite supplies it; the continuity fallback mu_wk0 = incr1 is used), incr1 = week-W+1
#' introductions (= the cumulative h=1 mu), incr2 = week-W+2 introductions
#' (cumulative h=2 − h=1). Already-affected zones stay NA.
#'
#' THE PROBABILITY. At a == 6 (the deployed geometry: the issue date ends its week) the 7- and
#' 14-day windows ARE the model's h = 1 and h = 2 outcome windows, so when `fc` carries
#' `p_invasion` it is passed through unchanged (and `p_case_invasion` likewise, else it mirrors
#' p_invasion): the model's own posterior-mean probability, on whatever scale the caller passed
#' (recalibrated, for run_all.R's featured forecast). It is NOT recomputed as
#' 1 − exp(−mu_forecast): mu_forecast is the posterior-mean hazard, and 1 − exp(−E[mu]) exceeds
#' E[1 − exp(−mu)] (Jensen). For a < 6, or when `fc` has no p_invasion, there is no model
#' probability for the shifted window, and p is that plug-in, 1 − exp(−mu_window).
#'
#' Approximation (a < 6 only): the uniform-within-week hazard treats introductions as evenly
#' spread across each weekly bucket; the partial current week's remaining-days hazard is
#' its whole-week rate scaled by the remaining fraction (unconditioned on the
#' partial non-observation).
#'
#' @param fc long weekly forecast tibble: method, health_zone, horizon (1,2),
#'   mu_forecast (CUMULATIVE expected introductions), was_active_before, optionally
#'   p_invasion / p_case_invasion / mu_wk0 / mobility_id / gt_profile.
#' @param forecast_date issue date (= ANALYSIS_DATE on the live path).
#' @param training_cutoff WEEK_ANCHOR-anchored start of the current week (the model's cutoff).
#' @return long tibble with horizon ∈ {1,2} now meaning "within 7d / within 14d of
#'   the analysis date", forecast_date, window_end, window_days, mu_forecast (windowed),
#'   p_invasion/p_case_invasion (see THE PROBABILITY), was_active_before, mobility_id,
#'   gt_profile.
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
  # window equal to the next full weekly bucket rather than a true a-shifted window. The
  # geometry is exact at a == 6 (the deployed case) and an approximation otherwise; it is
  # flagged rather than hidden, and the warning fires if a caller ever relies on it.
  # OUT-OF-RANGE IS AN ERROR, NOT A CLAMP. `a` is the issue date's position within its own
  # weekly bucket, so it is 0..6 by construction WHEN the grid is anchored to the analysis date.
  # If the grid has gone stale (forecast_date more than 6 days past the last training week, or
  # before it), min()/max() used to pin `a` silently at 6 or 0 — and the function then reported
  # its window as "the next 7 days" when part of that window was already in the PAST. A
  # forecast must not mis-date its own window without saying so.
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
                           "supplies mu_wk0, so the windows are approximations (the 7-day window ",
                           "is the next weekly bucket) and their probabilities are the plug-in ",
                           "1 - exp(-mu). Exact only when the issue date ends its week."),
                    as.integer(a)), call. = FALSE)
  fr_wk0 <- (6 - a) / 7      # fraction of the current week still ahead of the analysis date
  fr_hd  <- (a + 1) / 7      # fraction of the far boundary week inside the window

  wide <- fc %>%
    dplyr::filter(.data$horizon %in% c(1L, 2L)) %>%
    dplyr::select("method", "health_zone", "horizon", "mu_forecast") %>%
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

  # The model's own probabilities, per horizon, matched to z's rows by (method, health_zone):
  # used only where the window IS the model's horizon window (a == 6; see THE PROBABILITY).
  .p_model <- if (a == 6 && "p_invasion" %in% names(fc)) {
    .pc <- if ("p_case_invasion" %in% names(fc)) "p_case_invasion" else "p_invasion"
    fc %>%
      dplyr::filter(.data$horizon %in% c(1L, 2L)) %>%
      dplyr::distinct(.data$method, .data$health_zone, .data$horizon, .keep_all = TRUE) %>%
      dplyr::transmute(method = .data$method, health_zone = .data$health_zone,
                       horizon = as.integer(.data$horizon),
                       .p = .data$p_invasion, .p_case = .data[[.pc]])
  } else NULL

  build_win <- function(muv, win) {
    h   <- as.integer(win %/% 7L)                   # 1 = within 7d, 2 = within 14d
    muv <- ifelse(is.nan(muv), NA_real_, muv)       # affected/degenerate -> NA, not NaN
    if (is.null(.p_model)) {
      p <- ifelse(is.na(muv), NA_real_, 1 - exp(-muv)); p_case <- p
    } else {
      pm <- .p_model[.p_model$horizon == h, , drop = FALSE]
      i  <- match(paste(z$method, z$health_zone, sep = "\r"),
                  paste(pm$method, pm$health_zone, sep = "\r"))
      p <- pm$.p[i]; p_case <- pm$.p_case[i]
    }
    tibble::tibble(
      method            = z$method,
      health_zone       = z$health_zone,
      horizon           = h,
      forecast_date     = forecast_date,
      window_end        = forecast_date + win,
      window_days       = as.integer(win),   # LENGTH of the window, not a lead time (see header)
      mu_forecast       = muv,
      p_invasion        = p,
      p_case_invasion   = p_case,
      was_active_before = z$was_active_before,
      mobility_id       = if ("mobility_id" %in% names(z)) z$mobility_id else NA_character_,
      gt_profile        = if ("gt_profile" %in% names(z)) z$gt_profile else NA_character_)
  }
  dplyr::bind_rows(build_win(z$mu_7, 7L), build_win(z$mu_14, 14L))
}

# ---------------------------------------------------------------------------
# 2. issue_daily_reissue(): the live run's call (nothing is persisted)
# ---------------------------------------------------------------------------

#' Re-anchor the live run's weekly forecasts to windows measured from the issue date, report
#' the window dates, and return the windows.
#'
#' @param fc_all_current the current-week forecast bind from run_all.R (the featured Bayesian
#'   forecast, plus the stacked ensemble when it is a distinct object).
#' @param forecast_date  Date (= ANALYSIS_DATE).
#' @param training_cutoff Date (WEEK_ANCHOR-anchored start of the current week).
#' @return the windowed forecast tibble (anchor_windows_from_analysis_date()), invisibly.
#'   Nothing is written to disk.
issue_daily_reissue <- function(fc_all_current, forecast_date, training_cutoff) {
  fc_anch <- anchor_windows_from_analysis_date(fc_all_current, forecast_date,
                                               training_cutoff)
  # NO LONGER PERSISTED to a cumulative daily series. persist_daily_reissue() upserted into an
  # ever-growing daily_invasion_series file that (a) no script in the repo ever read, and
  # (b) accumulated rows across PROBABILITY SCALES, code versions and model grids with no
  # column recording which: the shipped copy mixed nine raw-scale July issues with one
  # recalibrated September issue, 22 methods from retired modules, and four dating columns
  # that were NA on every row. A product nothing reads and that silently mixes scales is worse
  # than no product. The windowed forecast is returned to the caller instead.
  .n_m <- length(unique(fc_anch$method))
  .n_z <- length(unique(fc_anch$health_zone[!is.na(fc_anch$p_invasion)]))
  message(sprintf(paste0("[daily_reissue] Issue %s (cutoff %s): invasion-probability windows for ",
                         "%d method(s), %d zone(s) at risk; returned to the caller, not persisted."),
                  forecast_date, training_cutoff, .n_m, .n_z))
  for (win in c(7L, 14L))
    message(sprintf("    within %d days of the analysis date: (%s, %s]",
                    win, forecast_date, forecast_date + win))
  invisible(fc_anch)
}

# ---------------------------------------------------------------------------
# 3. Leakage-free re-aggregation of the linelist as of an issue date
# ---------------------------------------------------------------------------

#' Re-aggregate zone-week counts as they were KNOWN on a given issue date.
#'
#' A back-test round must not reuse the final weekly counts: those bucket every case that
#' ever occurred in a week, including cases only reported later. To reconstruct the real-time
#' state on issue date `d` from the FINAL line list we keep only records whose OBSERVATION
#' date (linelist_observation_date()) is <= d, then bucket them by ONSET week (`date_index`)
#' exactly as `aggregate_to_zone_week()` does for the live pipeline. Cases with onset in week
#' W but not yet observed on d are absent — which is precisely the right-truncation the
#' nowcast is meant to correct.
#'
#' Observation date = the date the case becomes knowable to surveillance. We prioritise
#' **date_of_sample_collection** because a re-censored round is nowcast with the "asof"
#' truncation, which is measured on exactly this observation date (.trunc_fit_asof_cohort(),
#' 02_epi_params.R), so the censoring and the correction describe the same lag. Falls back to
#' lab_analysis_date then date_of_notification when the sample date is missing or
#' implausible, and a sitrep top-up row's `date_known` takes precedence over all three (see
#' linelist_observation_date()). (Gating instead on the later lab-confirmation date would
#' impose a different, longer delay than the as-of truncation measures, so it is deliberately
#' NOT used as the primary gate.)
#'
#' @param ll linelist (output of load_linelist): must carry date_index and at
#'   least one of lab_analysis_date / date_of_sample_collection /
#'   date_of_notification, plus confirmed / suspected classification columns.
#' @param zones canonical zone spine.
#' @param issue_date Date; the as-of date. It bounds the week grid, and with `censor = TRUE`
#'   only records observed on or before it are kept.
#' @param week_spine optional vector of WEEK_ANCHOR-anchored week starts that MUST be present in
#'   the output (zero-filled if unobserved by `issue_date`). This guarantees the
#'   training frame spans exactly through the cutoff week even when nothing has
#'   been observed in it yet — without it, an issue date before the first observed case of
#'   the cutoff week would drop that week entirely and misalign every horizon against its
#'   target.
#' @param censor TRUE (the default) applies the observation-date censoring above. FALSE skips
#'   it, for a line list that already IS the state on the issue date — an archived snapshot
#'   (asof_counts_source()). Every record a snapshot holds was in the extract on its date, so
#'   re-censoring it could only delete records whose recorded dates postdate the snapshot
#'   (data-entry errors), which the live forecast of that day would have counted. The week
#'   grid, the end-of-week guard and the spine padding are the same either way.
#' @return zone-week tibble (same schema as aggregate_to_zone_week output),
#'   censored to information available on `issue_date`.
reaggregate_asof <- function(ll, zones, issue_date, week_spine = NULL, censor = TRUE) {
  if (!is.logical(censor) || length(censor) != 1L || is.na(censor))
    stop("[reaggregate_asof] `censor` must be TRUE or FALSE.", call. = FALSE)
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
  ll_asof <- if (censor) {
    obs_date <- linelist_observation_date(ll, issue_date, caller = "reaggregate_asof")
    ll[!is.na(obs_date) & obs_date <= issue_date, , drop = FALSE]
  } else ll
  # aggregate_to_zone_week rebuilds week_start from date_index (onset) and returns
  # the zero-filled zone × week grid over the OBSERVED weeks, so downstream code
  # sees the identical schema it gets from the live path.
  # asof = issue_date, NOT the global ANALYSIS_DATE: this is a FOLD re-aggregation, so the
  # week grid must end at this fold's issue week, not at the run's analysis week. (For a
  # snapshot this also drops any week starting after the issue week, as the live pipeline
  # drops weeks after its own as-of week.)
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
#' same rule as the weekly training counts it projects. Priority: `date_known` where the line
#' list carries it (sitrep top-up rows), then sample collection, then lab analysis, then
#' notification (see reaggregate_asof() for why), each subject to the plausibility gate below.
#'
#' @param ll linelist (output of load_linelist, or of load_linelist_vintage(), whose
#'   `vintage_date` attribute caps the plausibility ceiling at the snapshot's own date).
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
  # observation date must lie between the outbreak start and the snapshot date — the run's
  # ANALYSIS_DATE for the live line list, the snapshot's own date for a replayed one (a
  # snapshot cannot hold a genuine observation made after it was taken); anything else falls
  # through to the next source column, and a record with no plausible date at all is
  # excluded and COUNTED rather than dated by a typo.
  # The floor is OUTBREAK_START minus a 30-day buffer, NOT OUTBREAK_START itself: the
  # earliest real sample in this line list (2026-04-28) precedes the nominal outbreak start
  # (2026-04-30), so a hard floor at OUTBREAK_START would reject genuine early samples and
  # push those records onto a later fallback column. The gate's job is to catch dates that
  # are IMPOSSIBLE (wrong year: 1999, 2023, 2024, 2025, 2027), not to adjudicate a date a
  # few days before the nominal start — a borderline-but-possible alert is left alone.
  # No date is invented when OUTBREAK_START is missing: the floor is undefined, so stop.
  .start <- suppressWarnings(as.Date(get0("OUTBREAK_START", ifnotfound = NA)))
  if (length(.start) != 1L || is.na(.start))
    stop(sprintf(paste0("[%s] OUTBREAK_START is not set (00_config.R); it defines the floor of ",
                        "the observation-date plausibility gate."), caller), call. = FALSE)
  .lo <- .start - 30L
  .hi <- suppressWarnings(as.Date(get0("ANALYSIS_DATE",  ifnotfound = issue_date)))
  if (length(.hi) != 1L || is.na(.hi)) .hi <- issue_date
  .vd <- suppressWarnings(as.Date(attr(ll, "vintage_date")))
  if (length(.vd) == 1L && !is.na(.vd)) .hi <- min(.hi, .vd)
  obs_date <- as.Date(rep(NA_real_, nrow(ll)), origin = "1970-01-01")
  # `date_known`, where present, is the date a record is KNOWN to have entered surveillance and
  # takes precedence: the sitrep top-up rows carry it under SITREP_APPEND_SHIFT_DAYS > 0, which
  # back-dates their specimen but must not make them visible before the report that created them.
  if ("date_known" %in% names(ll)) date_cols <- c("date_known", date_cols)
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
# 3b. "What was known on the issue date" from the ARCHIVED SNAPSHOTS (LFO_TRAINING_SOURCE)
# ---------------------------------------------------------------------------
# A back-test round is forecast from the data as they stood on its issue date. The archive holds
# the line list exactly as it stood on most days since 2026-06-09, so for those dates that state
# is READ, not rebuilt: the snapshot, processed by load_linelist() as of its own date
# (load_linelist_vintage(), 01_data_prep.R). It is nowcast as the live forecast of that day
# would have been: with the extract-regime truncation (onset -> appearance) ESTIMATED ON THAT
# DATE from the snapshots then available (vintage_truncation()), by the deployed estimator,
# window and thinning.
# A date with no snapshot (before the archive), or whose archive does not yet support the
# truncation fit (fewer than TRUNC_MIN_SNAPSHOTS weekly snapshots in the window), falls back to
# re-censoring the FINAL line list to the issue date (reaggregate_asof()), nowcast with the
# shared "asof" truncation. That reconstruction keeps later corrections (onset and zone edits,
# de-confirmations) and the cases sampled before the issue date but entered after it, so it is
# the second-best reconstruction and every row says which one it used.
# asof_counts_source() is the ONE dispatcher; the LFO training counts (16), the rolling
# predictors (20) and the h >= 2 R(t) (21) all go through it or through vintage_for_issue() /
# vintage_truncation(), so they cannot disagree about which data a round saw.

# The snapshot memo (.vintage_memo), list_linelist_snapshots() and load_linelist_vintage() live
# in 01_data_prep.R next to load_linelist(): the extract-regime truncation fit (02_epi_params.R)
# replays every archive vintage through them, and it runs before this file is sourced.

.vintage_snapshots <- function() {
  s <- .vintage_memo[["snapshots"]]
  if (is.null(s)) { s <- list_linelist_snapshots(); assign("snapshots", s, envir = .vintage_memo) }
  s
}

#' The snapshot that stands in for what was known on `issue_date`: the newest one dated ON OR
#' BEFORE it, provided it is at most `max_gap` days older. Never a later one -- that would hand
#' the round data it did not have.
#' @return list(folder, date, gap), or NULL when the run re-censors (LFO_TRAINING_SOURCE =
#'   "recensor"), when no snapshot qualifies, or when `issue_date` is on or after the run's
#'   analysis date (the caller's own line list IS that snapshot).
vintage_for_issue <- function(issue_date,
                              source = get0("LFO_TRAINING_SOURCE", ifnotfound = "vintage"),
                              max_gap = get0("LFO_VINTAGE_MAX_GAP_DAYS", ifnotfound = 2L),
                              snaps = NULL,
                              analysis_date = get0("ANALYSIS_DATE", ifnotfound = NA)) {
  if (!identical(source, "vintage")) return(NULL)
  issue_date <- as.Date(issue_date)
  stopifnot(length(issue_date) == 1L, !is.na(issue_date))
  ad <- suppressWarnings(as.Date(analysis_date))
  if (length(ad) == 1L && !is.na(ad) && issue_date >= ad) return(NULL)
  if (is.null(snaps)) snaps <- .vintage_snapshots()
  ok <- snaps$date <= issue_date
  if (!any(ok)) return(NULL)
  i <- max(which(ok))
  gap <- as.integer(issue_date - snaps$date[i])
  if (gap > as.integer(max_gap)) return(NULL)
  list(folder = snaps$folder[i], date = snaps$date[i], gap = gap)
}

#' The extract-regime truncation AS IT COULD HAVE BEEN ESTIMATED on `vdate`: the deployed
#' estimator (epinow2_truncation(), extract regime) on the snapshots dated up to `vdate`.
#' Memoised per date in this session (the fit itself is also cached on disk by its panel).
#' @return list(dist = the EpiNow2 dist_spec (for R(t)), spec = the nowcast's delay spec), or
#'   NULL when the archive does not yet support a fit or the fit is refused.
vintage_truncation <- function(vdate) {
  key <- paste0("trunc_", format(as.Date(vdate)), "_", .vintage_settings_key())
  if (exists(key, envir = .vintage_memo, inherits = FALSE)) return(get(key, envir = .vintage_memo))
  d <- tryCatch(epinow2_truncation(issue_date = as.Date(vdate), regime = "extract"),
                error = function(e) {
                  message(sprintf("[vintage] truncation at %s could not be fitted: %s",
                                  format(as.Date(vdate)), conditionMessage(e))); NULL })
  sp <- if (is.null(d)) NULL else .trunc_as_delay_spec(d, "extract")
  out <- if (is.null(sp)) NULL else list(dist = d, spec = sp)
  assign(key, out, envir = .vintage_memo)
  out
}

#' Load, in THIS process, every snapshot that `issue_dates` resolve to, and fit each snapshot's
#' real-time truncation. Call it before a parallel fan-out so forked workers inherit both and
#' the Stan fits run once, sequentially. No-op under LFO_TRAINING_SOURCE = "recensor".
#' @return data.frame(issue_date, folder, vintage_date, gap, truncation_mean) for the dates that
#'   resolved to a snapshot (truncation_mean NA where no real-time fit exists).
prefetch_vintages <- function(issue_dates) {
  issue_dates <- sort(unique(as.Date(issue_dates)))
  res <- lapply(issue_dates, function(d) {
    v <- vintage_for_issue(d)
    if (is.null(v)) return(NULL)
    load_linelist_vintage(v$folder, v$date)
    tr <- vintage_truncation(v$date)
    data.frame(issue_date = d, folder = v$folder, vintage_date = v$date, gap = v$gap,
               truncation_mean = if (is.null(tr)) NA_real_ else tr$spec$mean,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, res)
  if (is.null(out))
    out <- data.frame(issue_date = as.Date(character(0)), folder = character(0),
                      vintage_date = as.Date(character(0)), gap = integer(0),
                      truncation_mean = numeric(0))
  out
}

#' Raw zone-week counts as known on `issue_date`, plus the truncation and nowcast date that
#' match them. The ONE dispatcher for fold-side as-of reconstructions.
#'
#' @param ll the final line list (used whenever the round re-censors).
#' @param delay_recensor the shared "asof" truncation spec, for re-censored rounds.
#' @return list(raw, delay, nowcast_date, source, vintage_date, why). `raw` is
#'   reaggregate_asof()'s schema: the snapshot aggregated as it stands (no re-censoring), or
#'   the final line list censored to `issue_date`; nowcast it with
#'   apply_nowcast_correction(raw, analysis_date = nowcast_date, delay = delay). `source` is "vintage" or "recensor"; `why` says why a round
#'   re-censored ("no snapshot", "no real-time truncation", or "recensor requested").
asof_counts_source <- function(ll, zones, issue_date, week_spine = NULL, delay_recensor = NULL) {
  issue_date <- as.Date(issue_date)
  v  <- vintage_for_issue(issue_date)
  tr <- if (!is.null(v)) vintage_truncation(v$date)
  if (!is.null(v) && !is.null(tr)) {
    vll <- load_linelist_vintage(v$folder, v$date)
    # NOT re-censored: the snapshot IS what was known on its date (reaggregate_asof(), `censor`).
    raw <- reaggregate_asof(vll, zones, issue_date, week_spine = week_spine, censor = FALSE)
    return(list(raw = raw, delay = tr$spec, nowcast_date = v$date, source = "vintage",
                vintage_date = v$date, why = NA_character_))
  }
  raw <- reaggregate_asof(ll, zones, issue_date, week_spine = week_spine)
  list(raw = raw, delay = delay_recensor, nowcast_date = issue_date, source = "recensor",
       vintage_date = as.Date(NA),
       why = if (!identical(get0("LFO_TRAINING_SOURCE", ifnotfound = "vintage"), "vintage"))
               "recensor requested"
             else if (is.null(v)) "no snapshot" else "no real-time truncation")
}

message(paste0("[daily_reissue] 22_daily_reissue.R loaded — as-of reconstruction, snapshot ",
               "dispatch and issue-date windows."))
