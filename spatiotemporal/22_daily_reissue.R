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
# This module does two things and NOTHING else (the modelling machinery is
# unchanged — it lives in 15_workhorse.R / 16_invasion_eval.R):
#   1. LIVE:     stamp each daily run's forecasts with the issue date, the fixed
#                target ISO week, and the lead time IN DAYS, then persist them so
#                successive days ACCUMULATE into a series (rather than one file
#                being overwritten). This is the operational product.
#   2. BACKTEST: `run_invasion_lfo_daily()` validates the daily operating point by
#                re-issuing on intra-week dates. Crucially it (a) mirrors the LIVE
#                convention — the partial issue-week is IN the training set and the
#                issue happens mid-week — which the weekly LFO in 16_invasion_eval.R
#                does NOT do (it issues at the clean week boundary), and (b) is
#                leakage-free: it does not reuse the final weekly counts (those
#                already contain late-reported cases) but RE-AGGREGATES the linelist
#                censored to each issue date's laboratory/sample observation date.
#
# Lead-time convention: `lead_days` is the number of days from the issue date to
# the START of the target ISO week. On a Monday issue h=1 is ~7 days out; by
# Friday the same fixed target is ~3 days out. Always report lead time in days so
# "1-week-ahead" is never mistaken for "7 days from today".
# =============================================================================

# ---------------------------------------------------------------------------
# 1. Dating: attach issue date, fixed target week, and lead time (in days)
# ---------------------------------------------------------------------------

#' Stamp invasion forecasts with issue date, fixed target ISO week, and lead days.
#'
#' The target grid is anchored on the training cutoff (current week's WEEK_ANCHOR start),
#' NOT on the issue date, so re-issuing on any weekday keeps the target fixed.
#'
#' @param fc tibble with a `horizon` column (integer weeks ahead), e.g. the
#'   current-week forecast bind for all models.
#' @param forecast_date Date; the issue date (= ANALYSIS_DATE for the live run).
#' @param training_cutoff Date; WEEK_ANCHOR-anchored START of the current (cutoff) week
#'   (the week ends on ANALYSIS_DATE).
#' @return `fc` with added columns forecast_date, target_week_start,
#'   target_week_end, lead_days (issue -> target week start), lead_days_end
#'   (issue -> target week end).
add_forecast_dating <- function(fc, forecast_date, training_cutoff) {
  stopifnot(is.data.frame(fc), "horizon" %in% names(fc))
  if (!inherits(forecast_date, "Date") || length(forecast_date) != 1L)
    stop("[daily_reissue] forecast_date must be a scalar Date.", call. = FALSE)
  if (!inherits(training_cutoff, "Date") || length(training_cutoff) != 1L)
    stop("[daily_reissue] training_cutoff must be a scalar Date.", call. = FALSE)
  # Guard the operating convention: the issue date must fall within the cutoff
  # ISO week [C, C+6]. Outside that window the "current partial week" framing and
  # the lead-day arithmetic below no longer hold.
  d_off <- as.numeric(forecast_date - training_cutoff)
  if (is.na(d_off) || d_off < 0 || d_off > 6)
    warning(sprintf(paste0("[daily_reissue] forecast_date (%s) is not within the ",
                           "cutoff ISO week starting %s (offset %s d); lead-time ",
                           "and target dating assume a same-week mid-week issue."),
                    forecast_date, training_cutoff, d_off), call. = FALSE)
  fc %>%
    dplyr::mutate(
      forecast_date     = forecast_date,
      training_cutoff   = training_cutoff,
      target_week_start = training_cutoff + 7L * .data$horizon,
      target_week_end   = training_cutoff + 7L * .data$horizon + 6L,
      lead_days         = as.numeric(training_cutoff + 7L * .data$horizon - forecast_date),
      lead_days_end     = as.numeric(training_cutoff + 7L * .data$horizon + 6L - forecast_date)
    )
}

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
  a <- min(max(as.numeric(forecast_date - training_cutoff), 0), 6)
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
      lead_days         = as.integer(win),
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

# Compact set of columns kept in the cumulative daily series. Keeping this narrow
# stops the series file from ballooning while retaining everything the daily
# monitoring / convergence views need.
DAILY_SERIES_COLS <- c(
  "forecast_date", "horizon", "window_days", "window_end", "lead_days",
  "method", "health_zone", "p_invasion", "p_case_invasion", "mu_forecast",
  "was_active_before", "mobility_id", "gt_profile"
)

#' Persist one day's invasion re-issue: dated snapshot + cumulative series.
#'
#' Writes a dated snapshot (never overwritten across days) and upserts the rows
#' into a cumulative series keyed by (forecast_date, method, health_zone,
#' horizon). Re-running the SAME `forecast_date` replaces that day's rows only —
#' so a same-day re-run is idempotent and a new day appends.
#'
#' @param fc_dated forecasts already passed through `add_forecast_dating()`.
#' @param forecast_date Date; the issue date (used for the snapshot filename and
#'   the idempotent upsert key).
#' @param out_dir directory for forecast outputs (default OUT_FORECASTS).
#' @return (invisibly) the updated cumulative series tibble.
persist_daily_reissue <- function(fc_dated, forecast_date,
                                  out_dir = get0("OUT_FORECASTS")) {
  stopifnot(is.data.frame(fc_dated))
  if (is.null(out_dir))
    stop("[daily_reissue] out_dir is NULL (OUT_FORECASTS unset).", call. = FALSE)
  if (!"forecast_date" %in% names(fc_dated))
    stop("[daily_reissue] fc_dated has no forecast_date column — run add_forecast_dating() first.",
         call. = FALSE)

  daily_dir <- file.path(out_dir, "daily")
  dir.create(daily_dir, showWarnings = FALSE, recursive = TRUE)

  # (a) immutable dated snapshot (full columns)
  snap_path <- file.path(daily_dir,
                         sprintf("invasion_forecasts_%s.rds", forecast_date))
  saveRDS(fc_dated, snap_path)

  # (b) cumulative series (compact columns), idempotent upsert on forecast_date
  keep <- intersect(DAILY_SERIES_COLS, names(fc_dated))
  new_rows <- fc_dated %>% dplyr::select(dplyr::all_of(keep))

  series_rds <- file.path(out_dir, "daily_invasion_series.rds")
  prior <- if (file.exists(series_rds)) tryCatch(readRDS(series_rds),
             error = function(e) NULL) else NULL
  if (!is.null(prior) && "forecast_date" %in% names(prior)) {
    # drop any existing rows for this issue date so a same-day re-run overwrites
    # only its own day rather than duplicating it
    prior <- prior %>% dplyr::filter(.data$forecast_date != !!forecast_date)
    series <- dplyr::bind_rows(prior, new_rows)
  } else {
    series <- new_rows
  }
  series <- series %>% dplyr::arrange(.data$forecast_date, .data$method,
                                      .data$health_zone, .data$horizon)
  saveRDS(series, series_rds)
  readr::write_csv(series, file.path(out_dir, "daily_invasion_series.csv"))

  message(sprintf("[daily_reissue] Issue %s persisted: %d rows this day; series now %d issue-dates, %d rows.",
                  forecast_date, nrow(new_rows),
                  dplyr::n_distinct(series$forecast_date), nrow(series)))
  invisible(series)
}

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
  persist_daily_reissue(fc_anch, forecast_date, out_dir)
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
  stopifnot(is.data.frame(ll))
  if (!inherits(issue_date, "Date") || length(issue_date) != 1L)
    stop("[reaggregate_asof] issue_date must be a scalar Date.", call. = FALSE)
  date_cols <- intersect(c("date_of_sample_collection", "lab_analysis_date",
                           "date_of_notification"), names(ll))
  if (length(date_cols) == 0L)
    stop("[reaggregate_asof] linelist has no observation-date column ",
         "(lab_analysis_date / date_of_sample_collection / date_of_notification).",
         call. = FALSE)
  # coalesce available observation dates in priority order
  obs_date <- as.Date(rep(NA_real_, nrow(ll)), origin = "1970-01-01")
  for (cc in date_cols) {
    v <- ll[[cc]]
    if (!inherits(v, "Date")) v <- suppressWarnings(as.Date(v))
    obs_date <- dplyr::coalesce(obs_date, v)
  }
  ll_asof <- ll[!is.na(obs_date) & obs_date <= issue_date, , drop = FALSE]
  # aggregate_to_zone_week rebuilds week_start from date_index (onset) and returns
  # the zero-filled zone × week grid over the OBSERVED weeks, so downstream code
  # sees the identical schema it gets from the live path.
  zw <- aggregate_to_zone_week(ll_asof, zones)

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

# ---------------------------------------------------------------------------
# 4. Daily-issue backtest (mirrors the LIVE operating point)
# ---------------------------------------------------------------------------

#' Analysis-date-anchored backtest: re-issue on intra-week dates, score rolling
#' windows measured FROM the issue date (within 7 / 14 days), leakage-free.
#'
#' Mirrors the LIVE operating point exactly: every issue date `d` sits INSIDE its
#' week W = floor to the WEEK_ANCHOR start of d; the partial week W is in the training set
#' (nowcast-corrected as of d); the forecast is re-anchored to cumulative windows
#' (d, d+7] and (d, d+14] via `anchor_windows_from_analysis_date`. The OUTCOME is
#' onset-defined at daily resolution — a zone "invades within N days" iff its
#' first-ever confirmed ONSET (date_index) falls in (d, d+N] — matching the
#' pipeline's onset-based target, and evaluable against actual onset dates.
#'
#' Leakage control: TRAINING counts are re-aggregated from the line list censored
#' to the issue date's sample-observation date (`reaggregate_asof`), never the final
#' weekly counts. The OUTCOME uses the FINAL onset dates (that is correct — the
#' target is what truly happened), scored only for windows old enough that
#' `analysis_date - (d + N) >= min_eval_age_days`.
#'
#' @param zone_week_raw final zone-week counts (defines the week spine only).
#' @param ll line list: health_zone, date_index (onset), confirmed, plus the
#'   sample/lab/notification dates reaggregate_asof needs.
#' @param models named list of forecasters `function(zw_nc, t_idx, horizons, cutoff)`.
#' @param zones canonical zone spine.
#' @param windows_days rolling windows in days from the issue date (default 7, 14).
#' @param issue_offsets integer days after each week's start (WEEK_ANCHOR) to issue on
#'   (default c(0,3,6) ~ Mon/Thu/Sun; use 0:6 for every day at higher cost).
#' @param nowcast_fn deterministic per-fold nowcast (epinowcast is not refit per
#'   issue-fold — intractable; the live forecast refits it once per day).
#' @return pooled tibble across (issue_date × model × window) with p_invasion,
#'   is_new_invasion, was_active_before, fold_id, issue_date, cutoff, window_days,
#'   lead_days (= window_days), method, horizon (1 = 7d, 2 = 14d).
#'
#' NB on scoring: a fixed (zone, window, cutoff-week) event recurs once per issue
#' offset, so pooling directly over `evaluate_invasion()` double-counts events and
#' the zone-cluster bootstrap treats correlated re-issues as independent. Evaluate
#' per issue offset, or cluster the bootstrap by `cutoff`/`fold_id`.
run_invasion_lfo_daily <- function(zone_week_raw, ll, models, zones,
                                   windows_days = c(7L, 14L),
                                   min_train_weeks = 3L, min_atrisk_events = 1L,
                                   min_eval_age_days = 10L,
                                   issue_offsets = c(0L, 3L, 6L),
                                   ensemble_members = NULL,
                                   analysis_date = ANALYSIS_DATE,
                                   nowcast_fn = NULL) {
  if (is.null(nowcast_fn)) nowcast_fn <- get0("apply_nowcast_correction")
  stopifnot(!is.null(nowcast_fn), exists("anchor_windows_from_analysis_date"),
            exists("reaggregate_asof"))
  if (!all(c("health_zone", "date_index", "confirmed") %in% names(ll)))
    stop("[invasion_lfo_daily] linelist needs health_zone, date_index, confirmed.", call. = FALSE)
  weeks <- sort(unique(zone_week_raw$week_start))
  if (length(weeks) < min_train_weeks + 1L) {
    warning("[invasion_lfo_daily] Too few weeks for any fold."); return(NULL)
  }
  windows_days <- sort(as.integer(windows_days))

  # First confirmed ONSET date per zone (the ground-truth outcome clock).
  first_onset <- ll %>%
    dplyr::filter(.data$confirmed > 0, !is.na(.data$date_index),
                  .data$health_zone %in% zones) %>%
    dplyr::group_by(.data$health_zone) %>%
    dplyr::summarise(onset1 = min(.data$date_index), .groups = "drop")
  onset_of <- stats::setNames(first_onset$onset1, first_onset$health_zone)  # Date by zone
  # zones invaded (by onset) on/before a date; and newly invaded within (d, d+N].
  invaded_by <- function(d) names(onset_of)[!is.na(onset_of) & onset_of <= d]
  invaded_in <- function(d, N) names(onset_of)[!is.na(onset_of) & onset_of > d &
                                                 onset_of <= d + N]

  # Model weekly horizons needed so the anchor transform has mu_wk0 (week W), the
  # h=1 increment (week W+1) and the h=2 increment (week W+2) for windows up to 14 d.
  model_horizons <- c(1L, 2L)

  cand_weeks <- weeks[weeks >= weeks[min_train_weeks]]
  out <- list(); fid <- 0L
  for (wi in seq_along(cand_weeks)) {
    W <- cand_weeks[wi]            # index, NOT `for (W in cand_weeks)` — Date unclassing hazard
    for (off in issue_offsets) {
      d <- W + off
      if (d > analysis_date) next  # never issue in the future of the extract
      # windows whose end is old enough to be reliable ground truth as of the extract
      complete_win <- windows_days[
        as.numeric(analysis_date - (d + windows_days)) >= min_eval_age_days]
      if (length(complete_win) == 0L) next
      wmax_c <- max(complete_win)

      atrisk <- setdiff(zones, invaded_by(d))          # at-risk as of the issue date (onset-based)
      inv_max <- intersect(invaded_in(d, wmax_c), atrisk)
      if (length(inv_max) < min_atrisk_events) next    # need an event to define ranking

      fid <- fid + 1L
      # leakage-free training counts as known on issue date d, nowcast as of d.
      # week_spine forces every week up to and INCLUDING W present (zero-filled),
      # so the model's last training week is exactly W (weekly horizon alignment).
      spine <- weeks[weeks <= W]
      train_raw <- tryCatch(reaggregate_asof(ll, zones, d, week_spine = spine),
                            error = function(e) {
                              warning(sprintf("[invasion_lfo_daily] reaggregate_asof failed (issue %s): %s",
                                              d, conditionMessage(e))); NULL })
      if (is.null(train_raw)) next
      train_raw <- train_raw %>% dplyr::filter(week_start <= W)
      train_nc  <- tryCatch(nowcast_fn(train_raw, analysis_date = d),
                            error = function(e) train_raw)
      # Deterministic nowcast sets confirmed_nc = NA for a week whose midpoint has not
      # passed at d; the live forecast always keeps the current week as origin, so fall
      # back to censored raw (observed by d, leakage-free, uninflated) where it declined.
      if (!"confirmed_nc" %in% names(train_nc)) train_nc$confirmed_nc <- train_nc$confirmed
      if (!"suspected_nc" %in% names(train_nc)) train_nc$suspected_nc <- train_nc$suspected
      train_nc <- train_nc %>%
        dplyr::mutate(confirmed_nc = dplyr::coalesce(.data$confirmed_nc, .data$confirmed),
                      suspected_nc = dplyr::coalesce(.data$suspected_nc, .data$suspected))
      if (max(train_nc$week_start) != W) {
        warning(sprintf("[invasion_lfo_daily] last training week %s != cutoff %s (issue %s); skipping.",
                        max(train_nc$week_start), W, d)); next
      }
      t_idx <- length(unique(train_nc$week_start))
      if (t_idx < min_train_weeks) next

      # onset-window outcome, for at-risk zones only, per complete window
      truth <- lapply(complete_win, function(N) {
        inv_N <- intersect(invaded_in(d, N), atrisk)
        tibble::tibble(health_zone = atrisk, window_days = as.integer(N),
                       is_new_invasion = as.integer(atrisk %in% inv_N))
      }) %>% dplyr::bind_rows()

      # Collect every model's WEEKLY member forecast for this issue-fold.
      weekly <- lapply(names(models), function(mname) {
        fc <- tryCatch(models[[mname]](train_nc, t_idx, model_horizons, W),
                       error = function(e) {
                         warning(sprintf("[invasion_lfo_daily] %s failed (issue %s): %s",
                                         mname, d, conditionMessage(e))); NULL })
        if (is.null(fc) || nrow(fc) == 0) return(NULL)
        fc$method <- mname; fc
      })
      weekly <- Filter(Negate(is.null), weekly)
      if (length(weekly) == 0L) next
      weekly_all <- dplyr::bind_rows(weekly)
      # Build the ensembles on the WEEKLY member forecasts, BEFORE anchoring, so the
      # backtest scores the SAME combination the live product issues. Anchoring is
      # nonlinear (p = 1 - exp(-mu_window)), so ensembling AFTER anchoring would score
      # mean_m(1 - exp(-mu_m)) instead of the live 1 - exp(-mean_m mu_m) (Jensen gap).
      if (!is.null(ensemble_members) && exists("append_ensembles"))
        weekly_all <- append_ensembles(weekly_all, ensemble_members)
      # Re-anchor every method (members + ensembles) to (d, d+7] / (d, d+14].
      anch <- tryCatch(anchor_windows_from_analysis_date(weekly_all, d, W),
                       error = function(e) NULL)
      if (is.null(anch) || nrow(anch) == 0) next
      j <- anch %>%
        dplyr::filter(.data$health_zone %in% atrisk,
                      .data$window_days %in% complete_win,
                      !is.na(.data$p_invasion)) %>%   # drop model-masked zones
        dplyr::select(dplyr::any_of(c("method", "health_zone", "horizon", "window_days",
                                      "mu_forecast", "p_invasion"))) %>%
        dplyr::inner_join(truth, by = c("health_zone", "window_days")) %>%
        dplyr::mutate(fold_id = fid, cutoff = W, issue_date = d,
                      lead_days = .data$window_days, was_active_before = FALSE)
      out[[length(out) + 1L]] <- j
    }
  }
  if (length(out) == 0L) { warning("[invasion_lfo_daily] No evaluable issue-folds."); return(NULL) }
  res <- dplyr::bind_rows(out)
  message(sprintf("[invasion_lfo_daily] %d issue-folds across %d cutoff weeks; %d pooled rows; %d 7d events.",
                  dplyr::n_distinct(res$fold_id), dplyr::n_distinct(res$cutoff),
                  nrow(res), sum(res$is_new_invasion[res$window_days == min(windows_days)], na.rm = TRUE)))
  res
}

message("[daily_reissue] 22_daily_reissue.R loaded — daily re-issue + intra-week backtest.")
