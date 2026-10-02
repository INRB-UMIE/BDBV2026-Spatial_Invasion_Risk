# =============================================================================
# 04b_epinowcast.R — Probabilistic Nowcasting via epinowcast
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# ROLE (corrected 2026-09-19). This module is the SENSITIVITY arm, not the deployed nowcast.
# It used to say it "replaces the naive Exp-CDF right-truncation correction (04_nowcasting.R)",
# which has not been true since 2026-09-17: run_all.R applies apply_nowcast_correction() — the
# deterministic delay-CDF correction — on BOTH the deployed and the fold paths (that symmetry is
# the point: delta must be fitted and applied in the same regime), and epinowcast runs only
# behind RUN_EPINOWCAST_DIAGNOSTIC (default FALSE) plus the bayes_nowcast_sensitivity arm.
# It replaces nothing; it is the probabilistic cross-check on the deployed correction.

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
# ---------------------------------------------------------------------------
# Maximum onset->sample delay modelled by the nowcast — DERIVED, not hard-coded
# ---------------------------------------------------------------------------
# The bound D does two things inside epinowcast, and both matter:
#   (1) the delay distribution is modelled on [0, D] and RENORMALISED there, so the
#       expected eventual count for a reference date is observed * F(D) / F(lag).
#       A short D therefore under-states every incomplete week, systematically and in
#       one direction only;
#   (2) reference dates older than D are declared COMPLETE and receive no correction.

EPINOWCAST_MAX_DELAY_Q   <- 0.99   # quantile of the delay in force
EPINOWCAST_MAX_DELAY_MIN <- 21L    # floor: never model a shorter support than before
EPINOWCAST_MAX_DELAY_MAX <- 45L    # ceiling: matches the R(t) truncation clamp, and stays
                                   # well inside DELAY_MAX_PLAUSIBLE_DAYS (60 d)

#' The modelled delay support, in days.
#'
#' Precedence: an EPINOWCAST_MAX_DELAY environment override (for a deliberate
#' sensitivity run) > the fitted delay's `q` quantile, clamped to [lo, hi] > `lo`.
#' Never returns a non-finite or non-positive value.
epinowcast_max_delay <- function(delay = NULL, q = EPINOWCAST_MAX_DELAY_Q,
                                 lo = EPINOWCAST_MAX_DELAY_MIN,
                                 hi = EPINOWCAST_MAX_DELAY_MAX, quiet = FALSE) {
  env <- trimws(Sys.getenv("EPINOWCAST_MAX_DELAY", ""))
  if (nzchar(env)) {
    v <- suppressWarnings(as.integer(env))
    if (!is.na(v) && v >= 1L) {
      if (!quiet) message(sprintf("[epinowcast] max_delay = %d d (EPINOWCAST_MAX_DELAY override)", v))
      return(v)
    }
    warning(sprintf("[epinowcast] unreadable EPINOWCAST_MAX_DELAY='%s'; deriving it instead.", env),
            call. = FALSE)
  }
  if (is.null(delay))
    delay <- tryCatch(effective_onset_sample_delay(), error = function(e) NULL)
  qd <- if (is.null(delay) || !exists("delay_quantile", mode = "function")) NA_real_
        else tryCatch(delay_quantile(delay, q), error = function(e) NA_real_)
  if (!is.finite(qd)) {
    if (!quiet)
      message(sprintf("[epinowcast] max_delay = %d d (delay unavailable; floor used)", as.integer(lo)))
    return(as.integer(lo))
  }
  D <- as.integer(min(hi, max(lo, ceiling(qd))))
  if (!quiet)
    message(sprintf(paste0("[epinowcast] max_delay = %d d (P%.0f of the %s delay in force = %.1f d, ",
                           "clamped to [%d, %d]; covers %.1f%% of the delay mass)"),
                    D, 100 * q, if (is.null(delay$family)) "fitted" else delay$family, qd,
                    as.integer(lo), as.integer(hi),
                    100 * tryCatch(delay_cdf(delay, D), error = function(e) NA_real_)))
  D
}

EPINOWCAST_CHAINS      <- 4L
EPINOWCAST_ITER_WARMUP <- 2500L
EPINOWCAST_ITER_SAMPLE <- 2500L
EPINOWCAST_ADAPT_DELTA <- 0.965
# Reporting-delay distribution. Selected by a retrospective WIS comparison run ONCE, offline:
# compare_epinowcast_models() gave lognormal 5.49 < gamma 5.74 < exponential 8.98 on the
# held-out window, consistent with the branch's realtime default. NOT re-derived per run —
# that function is never called by the pipeline (see the note at EPINOWCAST_MAX_DELAY_Q), so
# these three numbers are a recorded historical result, not a claim about the current data.
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
                                     max_delay = epinowcast_max_delay()) {
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
fit_epinowcast <- function(triangle, max_delay = epinowcast_max_delay(),
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
  #
  # observed == 0 gives NA, NOT a ratio against the pmax(., 1) floor. With that floor the
  # "factor" for an unobserved week is the raw nowcast COUNT (e.g. 40), which downstream is
  # applied as a MULTIPLIER — the raw-count unit error this module has hit before. NA routes
  # the consumer to its `is.na(.nc_factor) -> 1` branch, i.e. no correction, which is the
  # honest answer for a week the triangle cannot see at all.
  wk[, factor := ifelse(observed > 0, pmax(nowcast_mean / observed, 1), NA_real_)]
  # Coefficient of variation of the weekly nowcast as an uncertainty flag.
  # NOTE ON `cv`: nowcast_lo/hi are SUMS OF DAILY QUANTILES, and a quantile of a sum is not the
  # sum of quantiles — summing seven daily 90% intervals overstates the weekly interval by
  # roughly sqrt(7) for near-independent days (seven days at mean 10, q5 5, q95 15 give
  # cv = 0.304 here against a true weekly cv of about 0.115). This is therefore an UPPER BOUND
  # on the weekly coefficient of variation, and it is used only to set the
  # correction_uncertainty flag, never to widen a predictive interval. Computing it properly
  # needs summary(fit, type = "nowcast_samples") and the 5th/95th percentile of the weekly sum
  # PER DRAW; left as an upper bound deliberately, since a conservative flag is the safe error.
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
                                      max_delay = epinowcast_max_delay(),
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

  # Resolve the modelled delay support ONCE and pass it to both stages. Letting each call
  # evaluate its own default would let the triangle and the model that consumes it disagree
  # if the fitted delay were refreshed between them — and a triangle built on one support
  # and pre-processed on another is silently mis-specified, not an error.
  .maxd <- epinowcast_max_delay()

  triangle <- tryCatch(
    build_reporting_triangle(linelist, analysis_date, outbreak_start, max_delay = .maxd),
    error = function(e) NULL
  )
  if (is.null(triangle)) return(do_fallback("could not build reporting triangle"))

  fit <- fit_epinowcast(triangle, max_delay = .maxd)
  if (is.null(fit)) return(do_fallback("epinowcast fit failed"))

  factors <- epinowcast_weekly_factors(fit, triangle)
  if (is.null(factors) || nrow(factors) == 0)
    return(do_fallback("could not extract nowcast factors"))

  message(sprintf("[epinowcast] Weekly correction factors (%d weeks): %s",
                  nrow(factors),
                  paste(sprintf("%s=%.2f", format(factors$week_start, "%m-%d"),
                                factors$factor), collapse = ", ")))

  # Carry the model's weekly COMPLETENESS RATIO and its CV.
  #
  # 2026-09-17 CORRECTION. This block used to carry nowcast_mean (an absolute count on the
  # reporting-TRIANGLE scale) and divide it by the ZONE-WEEK total to derive the factor. Those
  # two counts are different populations: the triangle holds only confirmed cases with BOTH an
  # onset and a sample date within max_delay, while zone_week additionally holds onset-imputed
  # and sitrep-appended cases. On the 2026-09-07 snapshot the triangle covered 4,963 of 6,669
  # confirmed cases (74%), so nowcast_mean / zone_week_total < 1 on almost every week and the
  # pmax(., 1) floor CLIPPED THE CORRECTION AWAY ENTIRELY: five of the six incomplete weeks got
  # a factor of exactly 1.000, and the raw vs epinowcast beta0 posteriors came out bit-identical
  # (outputs/forecasts/bayes_nowcast_sensitivity.csv). The pipeline was reporting a nowcast it
  # was not performing.
  #
  # A ratio of two incommensurable counts is not a completeness. The transferable quantity is
  # the DIMENSIONLESS ratio nowcast_mean / observed, both measured on the triangle's own base —
  # which epinowcast_weekly_factors() already computes (`factor`) and this block used to discard.
  # Applying that ratio to the zone-week base is scale-free and reconciles by construction.
  #
  # CAVEAT, deliberately not "corrected" away: the zone-week base is itself already partially
  # nowcast, because the onset imputation draws missing onsets backward from sample dates and so
  # populates recent onset-weeks that the triangle cannot see. The zone-week series is therefore
  # somewhat MORE complete than the triangle, and this ratio modestly over-corrects the most
  # recent week. Quantifying that is what bayes_nowcast_sensitivity is for.
  fac_lookup <- factors %>%
    dplyr::select(week_start, cv, .nc_factor = factor)

  # Spatial prior for the ADDITIVE-deficit correction below: each zone's share of
  # confirmed cases over recent (but not the current, most-truncated) weeks.
  .cur_wk <- suppressWarnings(max(zone_week$week_start, na.rm = TRUE))
  .recent <- zone_week %>%
    dplyr::filter(week_start >= .cur_wk - 35, week_start <= .cur_wk - 7) %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(.w = sum(confirmed, na.rm = TRUE), .groups = "drop")
  # FALLBACK (was: uniform .w = 1 over EVERY zone). A uniform prior spreads the additive
  # residual across all 519 zones — including zones that have NEVER had a confirmed case —
  # giving them confirmed_nc > 0 where confirmed == 0. forecast_workhorse() then reads that
  # as "already affected" and silently removes them from the at-risk set, corrupting both
  # the invasion outcome and the beta fit; it also contradicts this function's own
  # documented guarantee that it "never fabricates cases into never-recently-active zones".
  # Degrade instead to the ALL-TIME confirmed share, which is still zero for never-affected
  # zones; if even that is empty there is no epidemic to redistribute, so the additive step
  # is disabled (share 0 everywhere) rather than invented.
  if (nrow(.recent) == 0 || sum(.recent$.w) == 0) {
    .recent <- zone_week %>%
      dplyr::group_by(health_zone) %>%
      dplyr::summarise(.w = sum(confirmed, na.rm = TRUE), .groups = "drop")
    if (nrow(.recent) == 0 || sum(.recent$.w) == 0) {
      .recent <- dplyr::distinct(zone_week, health_zone) %>% dplyr::mutate(.w = 0)
      message("[epinowcast] No confirmed cases anywhere in the grid — additive residual ",
              "redistribution disabled (no zone receives a fabricated case).")
    } else {
      message("[epinowcast] No cases in the recent window; residual redistributed by the ",
              "ALL-TIME confirmed share (never-affected zones still receive nothing).")
    }
  }
  .recent$.share <- if (sum(.recent$.w) > 0) .recent$.w / sum(.recent$.w) else 0

  # National ZONE-WEEK confirmed total per week — the base actually being scaled (may differ from
  # the reporting-triangle observed). The completeness factor is derived from THIS base so the
  # corrected national total reconciles to nowcast_mean regardless of any triangle/zone-week
  # mismatch, and no double-counting occurs.
  .zw_tot <- zone_week %>%
    dplyr::group_by(week_start) %>%
    dplyr::summarise(.zw_conf = sum(confirmed, na.rm = TRUE), .groups = "drop")
  # Cap the multiplier (matches the deterministic nowcast's 5x) so a heavily under-observed week
  # cannot explode the multiplicative step or the suspected series; any shortfall against
  # nowcast_mean is made up by the additive residual below.
  # NOTE: the residual does NOT lift a zero-count week. The 2026-09-17 rewrite made it
  # proportional to the clipped EXCESS (.zw_conf * pmin(pmax(factor,1) - .fac_cap, .fac_cap)),
  # which is identically 0 whenever factor <= 5 — so a week with confirmed = 0 stays 0, as it
  # must: with no observed cases there is no spatial pattern to distribute a correction over,
  # and fabricating one would invent cases in zones with no evidence. (The previous wording
  # promised the opposite and described code that no longer exists.)
  .fac_cap <- 5

  out <- zone_week %>%
    dplyr::left_join(fac_lookup, by = "week_start") %>%
    dplyr::left_join(.zw_tot, by = "week_start") %>%
    dplyr::left_join(dplyr::select(.recent, health_zone, .share), by = "health_zone") %>%
    dplyr::mutate(
      cv     = dplyr::coalesce(cv, 0),
      # DIMENSIONLESS completeness factor from the triangle's own base (>= 1: a nowcast cannot
      # undercut what has already been reported), capped. Weeks with no nowcast row (older, fully
      # observed) get 1. Being scale-free, it applies to the zone-week base without reconciliation.
      factor = dplyr::if_else(is.na(.nc_factor), 1,
                              pmin(pmax(.nc_factor, 1), .fac_cap)),
      # trunc_weight is the reciprocal of the CAPPED multiplier, i.e. the completeness of the
      # MULTIPLICATIVE step only. If the cap ever binds (raw factor > 5) the realised national
      # correction exceeds 5x through the additive residual, so trunc_weight would read 0.2
      # and UNDERSTATE the correction actually applied — run_all.R prints mean(trunc_weight)
      # as the nowcast summary and beta_weighting = "completeness" consumes it. The schema is
      # deliberately NOT widened for this (zone_week flows into many consumers); instead the
      # cap-binding weeks are reported explicitly below, so the discrepancy cannot be silent.
      trunc_weight = 1 / factor,
      # The multiplicative step now reconciles by construction, so the additive residual carries
      # ONLY the mass the 5x cap clips off. It is redistributed by the recent-confirmed share
      # rather than proportionally, because a week that needs >5x correction is one where the
      # observed spatial pattern is least informative.
      # Residual = the mass the 5x cap clips. Bounded by a second cap so a pathological factor
      # cannot fabricate an unbounded case count through the additive channel.
      .resid = dplyr::if_else(is.na(.nc_factor), 0,
                              dplyr::coalesce(.zw_conf, 0) *
                                pmin(pmax(pmax(.nc_factor, 1) - .fac_cap, 0), .fac_cap)),
      confirmed_nc = confirmed * factor + .resid * dplyr::coalesce(.share, 0),
      suspected_nc = if ("suspected" %in% names(.)) suspected * factor else NA_real_,
      correction_uncertainty = cv > 0.25
    ) %>%
    # `cv` is RETAINED (renamed) rather than dropped: it is the coefficient of variation of
    # the weekly NATIONAL nowcast total implied by epinowcast's 90% interval, and it is the
    # only magnitude of nowcast uncertainty this pipeline produces. Discarding it forced every
    # downstream model to treat a corrected count as if it were observed data.
    dplyr::select(-factor, -.nc_factor, -.zw_conf, -.share, -.resid) %>%
    dplyr::rename(nowcast_cv = cv)

  # Report any week where the 5x cap bound, naming the raw factor: on those weeks
  # trunc_weight is not the completeness actually achieved (see the note at its definition).
  .capped <- factors[is.finite(factors$factor) & factors$factor > .fac_cap, , drop = FALSE]
  if (nrow(.capped))
    warning(sprintf(paste0("[epinowcast] the %gx multiplier cap bound on %d week(s) (raw factor ",
                           "up to %.2f: %s); trunc_weight reports %.2f there, understating the ",
                           "realised correction, which the additive residual completes."),
                    .fac_cap, nrow(.capped), max(.capped$factor, na.rm = TRUE),
                    paste(format(.capped$week_start), collapse = ", "), 1 / .fac_cap),
            call. = FALSE)

  attr(out, "nowcast_method") <- "epinowcast"
  attr(out, "epinowcast_factors") <- factors
  out
}
