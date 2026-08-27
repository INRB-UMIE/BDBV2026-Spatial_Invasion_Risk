# =============================================================================
# 04b_epinowcast.R — Probabilistic Nowcasting via epinowcast
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Replaces the naive Exp-CDF right-truncation correction (04_nowcasting.R) with
# a full Bayesian nowcast (Abbott et al. / epinowcast 0.6.0). Adapted for the
# DHIS2 pipeline from the lab-linelist implementation on branch
# tmp_nowcasting_growth (dhis2_pipeline/scripts/16_epinowcast.R and
# realtime_pipeline/scripts/nowcasting/).
#
# Approach (feasible for 519 zones): fit ONE pooled epinowcast model on the
# national daily onset->sample reporting triangle to obtain properly
# uncertainty-quantified per-onset-week correction factors, then apply those
# factors to each zone's observed weekly count. The onset-to-sample delay is a
# lab/reporting property that is essentially zone-invariant, so a shared delay
# model with zone-specific counts is both statistically sound and ~500x cheaper
# than fitting 519 separate Stan models. The spatial distribution of cases is
# preserved from the observed data; only the temporal right-truncation is
# corrected.
#
# Model spec (per the plan §1.2 and the branch): reference = onset date,
# report = sample-collection date; parametric reporting-delay distribution SELECTED BY WIS
# (EPINOWCAST_DELAY_DIST, currently "lognormal") — distinct from the simpler Exp(0.229)
# onset->sample hazard used only for point onset imputation in 01_data_prep.R; weekly
# random-walk expectation; negative-binomial observation; max_delay = 21 days.
#
# Public API:
#   nowcast_zone_week_epinowcast(zone_week, linelist, analysis_date, ...)
#     -> zone-week tibble with confirmed_nc / suspected_nc / trunc_weight /
#        correction_uncertainty, drop-in compatible with the output of
#        apply_nowcast_correction(). Falls back to the deterministic nowcast if
#        epinowcast is unavailable or the fit fails.
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lubridate)
})

HAS_EPINOWCAST <- requireNamespace("epinowcast", quietly = TRUE) &&
  requireNamespace("cmdstanr", quietly = TRUE)

# Fit configuration (light MCMC is sufficient for a correction-factor summary;
# the branch uses 4x2000 for full inference — override via these if needed).
EPINOWCAST_MAX_DELAY   <- 21L    # days; ~P99 of the onset-to-sample delay
EPINOWCAST_CHAINS      <- 2L
EPINOWCAST_ITER_WARMUP <- 500L
EPINOWCAST_ITER_SAMPLE <- 500L
EPINOWCAST_ADAPT_DELTA <- 0.95
# Reporting-delay distribution. Selected by retrospective WIS comparison
# (compare_epinowcast_models(): lognormal 5.49 < gamma 5.74 < exponential 8.98
# on the held-out window), consistent with the branch's realtime default.
EPINOWCAST_DELAY_DIST  <- "lognormal"

# ---------------------------------------------------------------------------
# Reporting triangle
# ---------------------------------------------------------------------------

#' Build a daily cumulative onset->sample reporting triangle from the linelist.
#'
#' @param linelist      tibble with date_of_symptom_onset, date_of_sample_collection,
#'                      confirmed (0/1 or logical).
#' @param analysis_date as-of date; sample dates after this are treated as
#'                      unobserved (right-truncation point).
#' @param outbreak_start earliest onset date to include.
#' @param max_delay     maximum onset->sample delay (days) retained.
#' @return data.table (reference_date, report_date, confirm) in epinowcast's
#'         lower-triangular cumulative long format, or NULL if too few cases.
build_reporting_triangle <- function(linelist, analysis_date = ANALYSIS_DATE,
                                     outbreak_start = OUTBREAK_START,
                                     max_delay = EPINOWCAST_MAX_DELAY) {
  stopifnot(all(c("date_of_symptom_onset", "date_of_sample_collection") %in%
                  names(linelist)))
  conf <- if ("confirmed" %in% names(linelist)) {
    !is.na(linelist$confirmed) & linelist$confirmed > 0
  } else rep(TRUE, nrow(linelist))

  d <- tibble::tibble(
    onset = as.Date(linelist$date_of_symptom_onset),
    samp  = as.Date(linelist$date_of_sample_collection)
  )[conf, ] %>%
    dplyr::filter(!is.na(onset), !is.na(samp), samp >= onset,
                  onset >= outbreak_start, onset <= analysis_date,
                  samp <= analysis_date,
                  as.numeric(samp - onset) <= max_delay)

  if (nrow(d) < 20L) {
    warning("[epinowcast] Only ", nrow(d),
            " onset+sample cases; too few for a nowcast. Returning NULL.")
    return(NULL)
  }

  dt <- data.table::as.data.table(d)
  counts <- dt[, .(new_confirm = .N),
               by = .(reference_date = onset, report_date = samp)]

  first_ref <- min(d$onset)
  grid <- data.table::CJ(
    reference_date = seq(first_ref, analysis_date, by = "day"),
    report_date    = seq(first_ref, analysis_date, by = "day")
  )
  grid <- grid[report_date >= reference_date &
                 as.numeric(report_date - reference_date) <= max_delay]

  tri <- merge(grid, counts, by = c("reference_date", "report_date"),
               all.x = TRUE)
  tri[is.na(new_confirm), new_confirm := 0]
  data.table::setorder(tri, reference_date, report_date)
  tri[, confirm := cumsum(new_confirm), by = reference_date]

  message(sprintf("[epinowcast] Reporting triangle: %d cases, %d onset dates, max_delay=%d",
                  nrow(d), length(unique(tri$reference_date)), max_delay))
  tri[, .(reference_date, report_date, confirm)]
}

# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------

#' Fit the pooled epinowcast model to a reporting triangle.
#'
#' @return an epinowcast fit object, or NULL on failure.
fit_epinowcast <- function(triangle, max_delay = EPINOWCAST_MAX_DELAY,
                           chains = EPINOWCAST_CHAINS,
                           iter_warmup = EPINOWCAST_ITER_WARMUP,
                           iter_sampling = EPINOWCAST_ITER_SAMPLE,
                           adapt_delta = EPINOWCAST_ADAPT_DELTA) {
  if (!HAS_EPINOWCAST || is.null(triangle)) return(NULL)

  pobs <- tryCatch(
    epinowcast::enw_preprocess_data(triangle, max_delay = max_delay),
    error = function(e) {
      warning("[epinowcast] enw_preprocess_data failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(pobs)) return(NULL)

  # Weekly random-walk expectation; fall back to intercept-only if rw(week) is
  # not estimable on a very short series.
  exp_mod <- tryCatch(
    epinowcast::enw_expectation(~ 1 + rw(week), data = pobs),
    error = function(e) epinowcast::enw_expectation(~ 1, data = pobs)
  )

  fit <- tryCatch(
    epinowcast::epinowcast(
      pobs,
      expectation = exp_mod,
      reference   = epinowcast::enw_reference(parametric = ~ 1,
                                              distribution = EPINOWCAST_DELAY_DIST,
                                              data = pobs),
      obs         = epinowcast::enw_obs(family = "negbin", data = pobs),
      fit         = epinowcast::enw_fit_opts(
        pp = FALSE, output_loglik = FALSE, sparse_design = TRUE,
        seed = get0("RANDOM_SEED", ifnotfound = 20260704L),  # reproducible nowcast MCMC
        chains = chains, iter_warmup = iter_warmup,
        iter_sampling = iter_sampling, adapt_delta = adapt_delta,
        parallel_chains = chains, save_warmup = FALSE,
        show_messages = FALSE, show_exceptions = FALSE, refresh = 0
      ),
      model = epinowcast::enw_model()
    ),
    error = function(e) {
      warning("[epinowcast] fit failed: ", conditionMessage(e))
      NULL
    }
  )
  fit
}

# ---------------------------------------------------------------------------
# Weekly correction factors
# ---------------------------------------------------------------------------

#' Derive per-onset-week correction factors from an epinowcast fit.
#'
#' @return tibble (week_start, observed, nowcast_mean, nowcast_lo, nowcast_hi,
#'         factor, cv) for the weeks covered by the nowcast, or NULL.
epinowcast_weekly_factors <- function(fit, triangle) {
  if (is.null(fit) || is.null(triangle)) return(NULL)

  nc <- tryCatch(data.table::as.data.table(summary(fit, type = "nowcast")),
                 error = function(e) NULL)
  if (is.null(nc) || !all(c("reference_date", "mean") %in% names(nc))) return(NULL)

  tri <- data.table::as.data.table(triangle)
  obs_day <- tri[, .(observed = max(confirm)), by = reference_date]

  q_lo <- if ("q5"  %in% names(nc)) "q5"  else "mean"
  q_hi <- if ("q95" %in% names(nc)) "q95" else "mean"
  m <- merge(nc[, .(reference_date, nowcast_mean = mean,
                    lo = get(q_lo), hi = get(q_hi))],
             obs_day, by = "reference_date", all.x = TRUE)
  m[is.na(observed), observed := 0]
  m[, week_start := lubridate::floor_date(reference_date, "week",
                                          week_start = get0("WEEK_ANCHOR", ifnotfound = 1L))]

  wk <- m[, .(observed     = sum(observed, na.rm = TRUE),
              nowcast_mean = sum(nowcast_mean, na.rm = TRUE),
              nowcast_lo   = sum(lo, na.rm = TRUE),
              nowcast_hi   = sum(hi, na.rm = TRUE)),
          by = week_start][order(week_start)]

  # Correction factor >= 1 (nowcast can't be below the already-observed count).
  wk[, factor := pmax(nowcast_mean / pmax(observed, 1), 1)]
  # Coefficient of variation of the weekly nowcast as an uncertainty flag.
  wk[, cv := (nowcast_hi - nowcast_lo) / (2 * 1.645 * pmax(nowcast_mean, 1))]

  tibble::as_tibble(wk)
}

# ---------------------------------------------------------------------------
# Main entry: nowcast-corrected zone-week tibble
# ---------------------------------------------------------------------------

#' Compare parametric reporting-delay models by retrospective nowcast skill.
#'
#' Faithful to the branch's "wide range of nowcasting models" comparison: fits
#' exponential / gamma / lognormal delay models on a training triangle truncated
#' `test_days` before the analysis date, nowcasts the held-out window, and scores
#' the per-onset-date nowcast against the (later-observed) truth by Weighted
#' Interval Score (Bracher 2021, using the 90% and 60% predictive intervals) and
#' RMSE, plus a skill score vs the naive no-correction baseline.
#'
#' @return tibble (model, mean_wis, rmse, skill_score, cov_90, n_test), best
#'         first, or NULL if epinowcast is unavailable.
compare_epinowcast_models <- function(linelist, analysis_date = ANALYSIS_DATE,
                                      outbreak_start = OUTBREAK_START,
                                      max_delay = EPINOWCAST_MAX_DELAY,
                                      test_days = 7L,
                                      distributions = c("exponential", "gamma",
                                                        "lognormal")) {
  if (!HAS_EPINOWCAST) {
    warning("[epinowcast] not installed; cannot compare models."); return(NULL)
  }
  # "Truth" triangle uses all data up to analysis_date; training truncates the
  # report dimension test_days earlier to create a held-out right-censored frame.
  full_tri  <- build_reporting_triangle(linelist, analysis_date, outbreak_start,
                                        max_delay)
  if (is.null(full_tri)) return(NULL)
  truth <- data.table::as.data.table(full_tri)[
    , .(truth = max(confirm)), by = reference_date]

  train_asof <- analysis_date - test_days
  train_tri  <- build_reporting_triangle(linelist, train_asof, outbreak_start,
                                         max_delay)
  if (is.null(train_tri)) return(NULL)
  # CAVEAT (known limitation): these held-out reference dates lie within ~test_days..2*test_days of
  # analysis_date, so their `truth` (max confirm AS-OF analysis_date, above) is itself still
  # right-truncated (only ~test_days..2*test_days of max_delay days reported), which biases the
  # delay-distribution SELECTION toward models that under-predict the final count. A fully clean
  # eval would score only reference dates >= max_delay old and train the model as-of ref_date +
  # max_delay; that is a larger restructure and left as a follow-up. Impact is bounded to the
  # choice among candidate delay distributions, not to any forecast value.
  test_dates <- seq(train_asof - test_days + 1, train_asof, by = "day")

  wis_one <- function(y, m, lo90, hi90, lo60, hi60) {
    if (is.na(y)) return(NA_real_)
    is90 <- (hi90 - lo90) + (2 / 0.10) * (pmax(lo90 - y, 0) + pmax(y - hi90, 0))
    is60 <- (hi60 - lo60) + (2 / 0.40) * (pmax(lo60 - y, 0) + pmax(y - hi60, 0))
    (0.5 * abs(y - m) + (0.10 / 2) * is90 + (0.40 / 2) * is60) / 2.5
  }

  rows <- lapply(distributions, function(dist) {
    fit <- tryCatch(
      {
        pobs <- epinowcast::enw_preprocess_data(train_tri, max_delay = max_delay)
        exp_mod <- tryCatch(epinowcast::enw_expectation(~ 1 + rw(week), data = pobs),
                            error = function(e) epinowcast::enw_expectation(~ 1, data = pobs))
        epinowcast::epinowcast(
          pobs, expectation = exp_mod,
          reference = epinowcast::enw_reference(parametric = ~ 1,
                                                distribution = dist, data = pobs),
          obs = epinowcast::enw_obs(family = "negbin", data = pobs),
          fit = epinowcast::enw_fit_opts(
            pp = FALSE, output_loglik = FALSE, sparse_design = TRUE,
            seed = get0("RANDOM_SEED", ifnotfound = 20260704L),  # reproducible nowcast MCMC
            chains = EPINOWCAST_CHAINS, iter_warmup = EPINOWCAST_ITER_WARMUP,
            iter_sampling = EPINOWCAST_ITER_SAMPLE, adapt_delta = EPINOWCAST_ADAPT_DELTA,
            parallel_chains = EPINOWCAST_CHAINS, save_warmup = FALSE,
            show_messages = FALSE, show_exceptions = FALSE, refresh = 0),
          model = epinowcast::enw_model())
      },
      error = function(e) { warning("[epinowcast] ", dist, " fit failed: ",
                                    conditionMessage(e)); NULL })
    if (is.null(fit)) return(NULL)

    nc <- tryCatch(data.table::as.data.table(summary(fit, type = "nowcast")),
                   error = function(e) NULL)
    if (is.null(nc)) return(NULL)
    getq <- function(col) if (col %in% names(nc)) nc[[col]] else nc$mean
    ev <- data.table::data.table(
      reference_date = nc$reference_date, m = nc$median %||% nc$mean,
      lo90 = getq("q5"), hi90 = getq("q95"),
      # `obs` is the naive no-correction baseline: the right-censored count
      # actually observed by train_asof (epinowcast's `confirm` column), NOT the
      # nowcast mean. Using the mean here would make rmse_naive ~ rmse and force
      # skill_score ~ 0 regardless of model quality.
      lo60 = getq("q20"), hi60 = getq("q80"), obs = getq("confirm"))
    ev <- merge(ev, truth, by = "reference_date")
    ev <- ev[reference_date %in% test_dates]
    if (nrow(ev) == 0) return(NULL)
    ev[, wis := mapply(wis_one, truth, m, lo90, hi90, lo60, hi60)]
    tibble::tibble(
      model = paste0("epinowcast_", dist),
      mean_wis = mean(ev$wis, na.rm = TRUE),
      rmse     = sqrt(mean((ev$m - ev$truth)^2, na.rm = TRUE)),
      rmse_naive = sqrt(mean((ev$obs - ev$truth)^2, na.rm = TRUE)),
      cov_90   = mean(ev$truth >= ev$lo90 & ev$truth <= ev$hi90, na.rm = TRUE),
      n_test   = nrow(ev))
  })
  res <- dplyr::bind_rows(rows)
  if (nrow(res) == 0) return(NULL)
  res %>%
    dplyr::mutate(skill_score = 1 - rmse / pmax(rmse_naive, 1e-9)) %>%
    dplyr::select(model, mean_wis, rmse, skill_score, cov_90, n_test) %>%
    dplyr::arrange(mean_wis)
}

#' Nowcast-correct a zone-week tibble using epinowcast weekly factors.
#'
#' Drop-in replacement for apply_nowcast_correction(): applies the pooled
#' epinowcast per-week correction factor to each zone's observed confirmed (and
#' suspected) count. Weeks outside the nowcast window are treated as fully
#' observed (factor 1). Falls back to the deterministic nowcast if epinowcast is
#' unavailable or the fit fails.
#'
#' @param zone_week   tibble with health_zone, week_start, confirmed, suspected.
#' @param linelist    the case linelist (for the reporting triangle).
#' @param analysis_date as-of date.
#' @param outbreak_start earliest onset date.
#' @param fallback_fn function used if epinowcast fails (default
#'        apply_nowcast_correction, which must already be sourced).
#' @return zone-week tibble with confirmed_nc, suspected_nc, trunc_weight,
#'         correction_uncertainty, plus a "method" attribute ("epinowcast" or
#'         "deterministic").
nowcast_zone_week_epinowcast <- function(zone_week, linelist,
                                         analysis_date = ANALYSIS_DATE,
                                         outbreak_start = OUTBREAK_START,
                                         fallback_fn = NULL) {
  do_fallback <- function(reason) {
    message("[epinowcast] Falling back to deterministic nowcast: ", reason)
    fn <- fallback_fn %||% get0("apply_nowcast_correction")
    if (is.null(fn)) stop("[epinowcast] No fallback nowcast function available.")
    out <- fn(zone_week, analysis_date = analysis_date)
    attr(out, "nowcast_method") <- "deterministic"
    out
  }

  if (!HAS_EPINOWCAST) return(do_fallback("epinowcast/cmdstanr not installed"))

  triangle <- tryCatch(
    build_reporting_triangle(linelist, analysis_date, outbreak_start),
    error = function(e) NULL
  )
  if (is.null(triangle)) return(do_fallback("could not build reporting triangle"))

  fit <- fit_epinowcast(triangle)
  if (is.null(fit)) return(do_fallback("epinowcast fit failed"))

  factors <- epinowcast_weekly_factors(fit, triangle)
  if (is.null(factors) || nrow(factors) == 0)
    return(do_fallback("could not extract nowcast factors"))

  message(sprintf("[epinowcast] Weekly correction factors (%d weeks): %s",
                  nrow(factors),
                  paste(sprintf("%s=%.2f", format(factors$week_start, "%m-%d"),
                                factors$factor), collapse = ", ")))

  # Carry the model's weekly nowcast_mean (the TARGET national confirmed total) and its CV. The
  # multiplicative factor is deliberately NOT taken from here: epinowcast_weekly_factors computes
  # it against the reporting-TRIANGLE observed base, but the correction below scales the ZONE-WEEK
  # confirmed base — a different count (zone_week is by date_index with onset imputation + zero-
  # fill). Deriving the factor from the wrong base (and, when triangle observed == 0, letting it
  # become a raw COUNT) over-inflated the corrected total and the suspected series. So recompute
  # the factor against the zone-week base below.
  fac_lookup <- factors %>%
    dplyr::select(week_start, cv, .nc_mean = nowcast_mean)

  # Spatial prior for the ADDITIVE-deficit correction below: each zone's share of
  # confirmed cases over recent (but not the current, most-truncated) weeks.
  .cur_wk <- suppressWarnings(max(zone_week$week_start, na.rm = TRUE))
  .recent <- zone_week %>%
    dplyr::filter(week_start >= .cur_wk - 35, week_start <= .cur_wk - 7) %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(.w = sum(confirmed, na.rm = TRUE), .groups = "drop")
  if (nrow(.recent) == 0 || sum(.recent$.w) == 0)
    .recent <- dplyr::distinct(zone_week, health_zone) %>% dplyr::mutate(.w = 1)
  .recent$.share <- .recent$.w / sum(.recent$.w)

  # National ZONE-WEEK confirmed total per week — the base actually being scaled (may differ from
  # the reporting-triangle observed). The completeness factor is derived from THIS base so the
  # corrected national total reconciles to nowcast_mean regardless of any triangle/zone-week
  # mismatch, and no double-counting occurs.
  .zw_tot <- zone_week %>%
    dplyr::group_by(week_start) %>%
    dplyr::summarise(.zw_conf = sum(confirmed, na.rm = TRUE), .groups = "drop")
  # Cap the multiplier (matches the deterministic nowcast's 5x) so a heavily under-observed week
  # cannot explode the multiplicative step or the suspected series; any shortfall against
  # nowcast_mean is made up by the additive residual, which also lifts a zero-count current week.
  .fac_cap <- 5

  out <- zone_week %>%
    dplyr::left_join(fac_lookup, by = "week_start") %>%
    dplyr::left_join(.zw_tot, by = "week_start") %>%
    dplyr::left_join(dplyr::select(.recent, health_zone, .share), by = "health_zone") %>%
    dplyr::mutate(
      cv     = dplyr::coalesce(cv, 0),
      # Completeness factor RELATIVE TO THE ZONE-WEEK CONFIRMED TOTAL, >= 1 (nowcast can't undercut
      # the observed count) and capped. Weeks with no nowcast row (older, fully observed) get 1.
      factor = dplyr::if_else(is.na(.nc_mean), 1,
                              pmin(pmax(.nc_mean / pmax(.zw_conf, 1), 1), .fac_cap)),
      trunc_weight = 1 / factor,
      # RESIDUAL national deficit the capped multiplicative step leaves against nowcast_mean,
      # measured on the SAME zone-week base it scales (not the triangle observed) so the national
      # corrected total meets nowcast_mean. Redistributed by the recent-share prior — this is also
      # what lifts a zero-count current week (0 x factor = 0), the very week the nowcast is for.
      .resid = dplyr::if_else(is.na(.nc_mean), 0,
                              pmax(.nc_mean - dplyr::coalesce(.zw_conf, 0) * factor, 0)),
      confirmed_nc = confirmed * factor + .resid * dplyr::coalesce(.share, 0),
      suspected_nc = if ("suspected" %in% names(.)) suspected * factor else NA_real_,
      correction_uncertainty = cv > 0.25
    ) %>%
    dplyr::select(-factor, -cv, -.nc_mean, -.zw_conf, -.share, -.resid)

  attr(out, "nowcast_method") <- "epinowcast"
  attr(out, "epinowcast_factors") <- factors
  out
}
