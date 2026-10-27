# ============================================================================
# 04c_dhis2_delay_windows.R — DHIS2 line-list delay analysis (spatiotemporal)
#
# DHIS2 analog of the retired Ituri LAB-linelist delay analysis. It estimates the DHIS2
# reporting delays — primarily the
# onset->sample delay (REPORTED; since 2026-10-15 it no longer drives the missing-onset imputation)
# — with the SAME rigorous machinery, adapted to the DHIS2 line list and wired to
# the spatiotemporal 00_config.R (no dependency on the CFR pipeline's config_r.R).
#
# WHY a DHIS2-specific fit: the fixed lab-linelist delay (Exp rate 0.228, mean
# 4.39 d) is the wrong shape for the DHIS2 onset->sample delay (the fitted values are written
# to dhis2_onset_sample_delay_params.csv every run). The CONFIRMED-case stratum of this fit is
# the REPORTED onset->sample delay (ONSET_IMPUTE_STRATUM, 00_config.R); the pooled fit is the
# fallback of the nowcast and R(t) truncation. Neither drives the missing-onset imputation,
# which draws same-class donor intervals (01_data_prep.R).
#
# INPUT: the line list must have been passed through apply_earliest_sample_dates()
# (00_config.R), as run_all.R step 1a, the MAIN block below and 04d's standalone block do, so
# the sample date is the earliest diagnostic sample and a confirmed record's lab date is that
# of its confirming test.
#
# METHODS (mirroring the lab template):
#   * Windowing: onset in [ANALYSIS_START, TRUNC_DATE = max(sample_date) - TEST_DAYS], where
#     max(sample_date) is taken over complete onset/sample pairs on or before the as-of date
#     and TEST_DAYS = DELAY_TRUNC_BUFFER_DAYS (00_config.R). The final days before the
#     extract are heavily right-truncated (recent onsets have not yet been sampled), so delay
#     fitting stops that many days before the line-list end. Both dates are DERIVED from the
#     data (no hard-coded cutoffs). ANALYSIS_START is the Monday on or before OUTBREAK_START
#     (floor_date(OUTBREAK_START, "week", week_start = 1)).
#   * Interval-censored MLE  — fitdistcens(); d=0 -> [0,0.5]; d>0 -> [d-0.5,d+0.5]
#        (corrects daily-rounding censoring); AIC-selected across gamma/lnorm/
#        weibull/exp. Always run; the selected estimator when EpiDist is off, or when its fit
#        fails and no EpiDist fit of the same records is already on file.
#   * Bayesian EpiDist       — the DEFAULT (RUN_EPIDIST, 00_config.R, default TRUE; skipped
#        when the epidist/epiparameter/brms packages are absent). The MARGINAL model
#        corrects BOTH double interval censoring AND right-truncation, and is preferred
#        over the MLE wherever it was fitted.
#
# OUTPUTS (written next to the source line list, folder-versioned):
#   * dhis2_delay_params_censored.csv     — full per-family censored MLE table
#       (window, delay_type, family, n, aic, bic, implied_mean_d, params, best).
#   * dhis2_onset_sample_delay_params.csv — long (delay,quantity,value) params for
#       the onset->sample delay (MLE, pooled EpiDist marginal, per-classification strata,
#       their convergence diagnostics and the `selected_*` rows), mirroring the CFR
#       onset_to_sample CSV, read by effective_onset_sample_delay() (00_config.R).
#   * <OUT_DIAGNOSTICS>/dhis2_delay_fits_<TRUNC_DATE as YYYYMMDD>.pdf — top-2 censored-MLE
#       fits per delay; not on FIGURE_KEEP, so normally not saved.
#
# REUSABLE: estimate_dhis2_onset_sample_delay(ll, analysis_date) returns the fit
# list (best family, params, rate, windowed delay vector, censored table) so
# run_all.R step 1a (or an interactive session) can obtain the DHIS2 delay without
# re-running the whole script.
#
# Run: Rscript spatiotemporal/04c_dhis2_delay_windows.R
#   env RUN_EPIDIST=FALSE -> skip the Bayesian truncation-corrected model (MLE only).
# ============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(fitdistrplus)
})
# Figure packages are optional (fitting does not need them); load if present.
.HAVE_FIG <- all(vapply(c("patchwork", "scales"), requireNamespace, logical(1),
                        quietly = TRUE))
if (.HAVE_FIG) suppressPackageStartupMessages({ library(patchwork); library(scales) })
# EpiDist (Bayesian, truncation-correcting) is the DEFAULT delay estimator: 00_config.R
# (sourced above) sets RUN_EPIDIST from env RUN_EPIDIST (default TRUE). It is gated here on
# package availability, so on a box without epidist RUN_EPIDIST is FALSE and the
# interval-censored MLE is the selected delay.
.HAVE_EPIDIST <- requireNamespace("epidist", quietly = TRUE) &&
                 requireNamespace("epiparameter", quietly = TRUE) &&
                 requireNamespace("brms", quietly = TRUE)
RUN_EPIDIST <- isTRUE(RUN_EPIDIST) && .HAVE_EPIDIST
if (RUN_EPIDIST) suppressPackageStartupMessages({
  library(epidist); library(epiparameter); library(brms)
})


MAX_DELAY   <- get0("DELAY_MAX_PLAUSIBLE_DAYS", ifnotfound = 60L)  # shared plausibility ceiling (00_config.R): drop longer delays as typos/outliers before fitting
TEST_DAYS   <- get0("DELAY_TRUNC_BUFFER_DAYS", ifnotfound = 5L)  # right-truncation buffer (00_config.R): stop fitting this many days before the last sample
# Bayesian sampler settings (only used when RUN_EPIDIST=TRUE)
STAN_CHAINS <- 2L; STAN_ITER <- 1000L; STAN_WARMUP <- 500L
STAN_CORES  <- min(STAN_CHAINS, 4L)

# ============================================================================
# HELPER FUNCTIONS
# The interval-censored MLE core is REPLICATED VERBATIM from the validated retired lab
# template, so the DHIS2 fits use the identical,
# audited machinery — same censoring convention, multi-start, and AIC selection.
# ============================================================================

# Dates are parsed with .parse_date() (00_config.R): the ONE per-element parser (mixed ISO /
# d-m-y / m-d-y) that 01_data_prep.R also uses, so the delay is fitted on exactly the dates the
# line list is loaded with. (Base as.Date(tryFormats=) picks ONE format from the first non-NA
# element and NA-s every value in a different format.)

# fit_mle_naive() and naive_implied_mean() REMOVED: the naive uncensored MLE was fitted in the
# MAIN block and never used (the QA figure plots the censored fits only). The censored-MLE path
# reports `implied_mean_d` from the fit object directly (.fit_all_censored), and the EpiDist
# path reports the posterior mean.

# ── Interval-censored MLE ─────────────────────────────────────────────────────
.make_cens_df <- function(x) {
  stopifnot(is.numeric(x), all(x >= 0, na.rm = TRUE))
  data.frame(left  = ifelse(x == 0L, 0, as.numeric(x) - 0.5),
             right = as.numeric(x) + 0.5)
}

.get_cens_starts <- function(x, fam) {
  m <- mean(x, na.rm = TRUE); v <- var(x, na.rm = TRUE)
  if (!is.finite(m) || m <= 0) m <- 1
  if (!is.finite(v) || v <= 0) v <- m^2
  switch(fam,
    gamma = {
      sh <- m^2 / v; ra <- m / v
      if (!is.finite(sh) || sh <= 0) sh <- 1
      if (!is.finite(ra) || ra <= 0) ra <- 1 / m
      list(shape = sh, rate = ra)
    },
    lnorm = {
      ml <- log(m^2 / sqrt(v + m^2)); sl <- sqrt(log(1 + v / m^2))
      if (!is.finite(ml) || !is.finite(sl) || sl <= 0) { ml <- log(m); sl <- 0.5 }
      list(meanlog = ml, sdlog = sl)
    },
    weibull = {
      cv <- sqrt(v) / m
      k  <- if (is.finite(cv) && cv > 0 && cv < 10) cv^(-1.086) else 1
      sc <- tryCatch(m / gamma(1 + 1 / k), error = function(e) m)
      if (!is.finite(k) || k <= 0) k <- 1
      if (!is.finite(sc) || sc <= 0) sc <- m
      list(shape = k, scale = sc)
    },
    exp = list(rate = 1 / m)
  )
}

.fit_one_censored <- function(x, fam, delay_name = "") {
  if (length(x) < 5L) return(NULL)
  cdf <- .make_cens_df(x)
  m <- mean(x)
  # FAMILY-APPROPRIATE multi-start: gamma-shaped fallbacks (shape/rate) are invalid
  # parameter names for weibull (shape/scale) and lnorm (meanlog/sdlog). Give each
  # family its own alternative starts so a failed primary can still recover.
  sl <- switch(fam,
    exp     = list(list(rate = 1 / m)),
    gamma   = list(.get_cens_starts(x, fam), list(shape = 1, rate = 1 / m), list(shape = 2, rate = 2 / m)),
    lnorm   = list(.get_cens_starts(x, fam), list(meanlog = log(m), sdlog = 0.5),
                   list(meanlog = log(max(stats::median(x), 0.5)), sdlog = 1)),
    weibull = list(.get_cens_starts(x, fam), list(shape = 1, scale = m), list(shape = 1.5, scale = m)),
    list(.get_cens_starts(x, fam)))
  for (s in sl) {
    fit <- tryCatch(suppressWarnings(fitdistcens(cdf, fam, start = s)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      p <- fit$estimate
      # lnorm's meanlog is legitimately negative when median delay < 1 day, so test
      # only the strictly-positive parameters (sdlog); all others must be > 0.
      pos <- if (fam == "lnorm") p[["sdlog"]] else p
      if (all(is.finite(p)) && is.finite(fit$aic) && all(pos > 0)) return(fit)
    }
  }
  message(sprintf("  [%s/%s] all censored starts failed", delay_name, fam))
  NULL
}

.fit_all_censored <- function(x, delay_name) {
  x <- x[!is.na(x) & x >= 0L]
  if (length(x) < 5L) {
    cat(sprintf("  [%s] n=%d<5, skipping censored\n", delay_name, length(x)))
    return(NULL)
  }
  cat(sprintf("  %-24s n=%-5d n_zero=%-4d (%.1f%%)\n",
              delay_name, length(x), sum(x == 0L), sum(x == 0L) / length(x) * 100))
  out <- list()
  for (fam in c("gamma", "lnorm", "weibull", "exp")) {
    fit <- .fit_one_censored(x, fam, delay_name)
    if (!is.null(fit)) out[[fam]] <- fit
  }
  out
}

.implied_mean <- function(fam, params) {
  tryCatch(switch(fam,
    gamma   = unname(params["shape"] / params["rate"]),
    lnorm   = unname(exp(params["meanlog"] + 0.5 * params["sdlog"]^2)),
    weibull = unname(params["scale"] * gamma(1 + 1 / params["shape"])),
    exp     = unname(1 / params["rate"]),
    NA_real_), error = function(e) NA_real_)
}

.density_from_fit <- function(fam, params, x_grid) {
  tryCatch(switch(fam,
    gamma   = dgamma(x_grid,   shape   = params["shape"],   rate  = params["rate"]),
    lnorm   = dlnorm(x_grid,   meanlog = params["meanlog"], sdlog = params["sdlog"]),
    weibull = dweibull(x_grid, shape   = params["shape"],   scale = params["scale"]),
    exp     = dexp(x_grid,     rate    = params["rate"]),
    rep(NA_real_, length(x_grid))),
  error = function(e) rep(NA_real_, length(x_grid)))
}

.build_cens_tbl <- function(fits, delay_type, window_name, n_obs) {
  if (is.null(fits) || length(fits) == 0L) return(NULL)
  purrr::map_dfr(names(fits), function(fam) {
    fit <- fits[[fam]]; p <- fit$estimate
    tibble(window = window_name, delay_type = delay_type,
           family = toupper(fam), n = n_obs,
           aic = round(fit$aic, 4), bic = round(fit$bic, 4),
           implied_mean_d = round(.implied_mean(fam, p), 3),
           shape   = if ("shape"   %in% names(p)) round(unname(p["shape"]), 4)   else NA_real_,
           rate    = if ("rate"    %in% names(p)) round(unname(p["rate"]), 4)    else NA_real_,
           scale   = if ("scale"   %in% names(p)) round(unname(p["scale"]), 4)   else NA_real_,
           meanlog = if ("meanlog" %in% names(p)) round(unname(p["meanlog"]), 4) else NA_real_,
           sdlog   = if ("sdlog"   %in% names(p)) round(unname(p["sdlog"]), 4)   else NA_real_)
  }) %>%
    mutate(delta_aic = round(aic - min(aic), 4), best = delta_aic == 0) %>%
    arrange(aic)
}

.cens_best_fam <- function(f)
  if (!is.null(f) && length(f) > 0) names(sort(sapply(f, `[[`, "aic")))[1] else NA_character_

# ── Windowed delay populations — SINGLE definition ────────────────────────────
# The MAIN block below and 04d_delay_fit_figures.R both need the same four windowed
# delay populations. They used to be built by an inline `mk_delay()`/`delay_specs`
# pair inside the MAIN block, which meant the figure module could only get them by
# duplicating that logic — exactly the drift hazard this file already warns about
# for .parse_date(). Factored out here so the fits and the figures are guaranteed to
# describe the SAME records.
# `lbl` IS ASCII ON PURPOSE: it is written to the `label` column of
# dhis2_delay_selected.csv, so it must stay byte-identical to the previous inline
# definition. 04d_delay_fit_figures.R prettifies the arrow for display only.
# THE TWO LAB-REFERENCED DELAYS need both dates of a pair to describe the SAME test. The
# extract's sample and lab dates belong to a record's LATEST test; apply_earliest_sample_dates()
# (00_config.R), which every caller applies first, moves the sample date to the earliest
# diagnostic sample and, for a currently-confirmed record, the lab date to its confirming
# test's. A record whose sample date was moved while its lab date was NOT (every re-tested
# record that is not currently confirmed, and a confirmed one whose confirming test is the
# extract's) pairs an earlier test's sample with a later test's result, so
# build_dhis2_delay_populations() leaves it out of the two lab-referenced delays (counted in its
# message) rather than measure a delay across two tests. The other two delays keep it.
dhis2_delay_specs <- function() list(
  onset_sample        = list(from = "date_of_symptom_onset",     to = "date_of_sample_collection", lbl = "Onset -> sample"),
  onset_notification  = list(from = "date_of_symptom_onset",     to = "date_of_notification",      lbl = "Onset -> notification"),
  onset_lab_analysis  = list(from = "date_of_symptom_onset",     to = "lab_analysis_date",         lbl = "Onset -> lab analysis"),
  sample_lab_analysis = list(from = "date_of_sample_collection", to = "lab_analysis_date",         lbl = "Sample -> lab analysis"))

#' Build the windowed, plausibility-filtered delay populations for every DHIS2 delay.
#'
#' @param ll DHIS2 line list after apply_earliest_sample_dates() (character or already-parsed
#'   date columns).
#' @param analysis_date as-of date, or NULL. NULL applies no as-of bound (truncation
#'   reference = max observed date): the MAIN block's DHIS2_DELAY_NO_ASOF=1 path. When a
#'   date is supplied the same as-of rule as estimate_dhis2_onset_sample_delay() applies
#'   — secondary dates after the as-of date are dropped before the window is derived, so
#'   a back-dated re-run cannot leak future data into the figures.
#' @return list(specs, delays = named list(onset, sample, delay), analysis_start,
#'   trunc_date, obs_date, window_name, max_sample, n_records).
#'   `delays[[nm]]$delay` for "onset_sample" is IDENTICAL to
#'   estimate_dhis2_onset_sample_delay()$delays_windowed for the same arguments.
build_dhis2_delay_populations <- function(ll, analysis_date = NULL,
                                          outbreak_start = OUTBREAK_START,
                                          test_days = TEST_DAYS,
                                          max_delay = MAX_DELAY) {
  df <- as.data.frame(ll, stringsAsFactors = FALSE)
  date_cols <- c("date_of_symptom_onset", "date_of_sample_collection",
                 "date_of_notification", "reporting_date", "lab_analysis_date")
  for (dc in intersect(date_cols, names(df))) df[[dc]] <- .parse_date(df[[dc]])
  if (!"date_of_sample_collection" %in% names(df)) return(NULL)

  .asof    <- if (is.null(analysis_date)) NULL else .parse_date(analysis_date)
  has_asof <- !is.null(.asof) && length(.asof) == 1L && !is.na(.asof)
  # As-of bound on the SAMPLE date first: it defines the truncation window shared by all
  # four delays, so it must be applied before max(sample) is taken (otherwise a future-dated
  # sample typo would push the window past the as-of date).
  samp <- df$date_of_sample_collection
  if (has_asof) samp[!is.na(samp) & samp > .asof] <- NA
  # max(sample) over COMPLETE onset/sample pairs, exactly as estimate_dhis2_onset_sample_delay()
  # takes it. Over every record, a latest sample whose record has no onset would move this
  # window away from the fit's, and the figures would describe different records.
  if ("date_of_symptom_onset" %in% names(df)) samp[is.na(df$date_of_symptom_onset)] <- NA
  max_sample <- suppressWarnings(max(samp, na.rm = TRUE))
  if (!is.finite(as.numeric(max_sample))) return(NULL)
  trunc_date     <- max_sample - test_days
  analysis_start <- lubridate::floor_date(outbreak_start, unit = "week", week_start = 1L)
  # OBS_DATE (right-truncation reference for EpiDist) deliberately excludes
  # date_of_notification / reporting_date, which carry far-future data-entry typos.
  obs_date <- if (has_asof) .asof else
    suppressWarnings(max(c(df$date_of_sample_collection, df$date_of_symptom_onset), na.rm = TRUE))

  # CROSS-TEST PAIRS (see dhis2_delay_specs()): the sample date was moved by
  # apply_earliest_sample_dates() and the lab date was left as extracted, so the lab date is a
  # later test's. Identified from the `_extract` columns that function keeps; without them
  # (SAMPLE_DATE_SOURCE = "extract", or no archive) no date was moved and no pair is excluded.
  xtest <- rep(FALSE, nrow(df))
  if ("date_of_sample_collection_extract" %in% names(df)) {
    s_ex  <- .parse_date(df$date_of_sample_collection_extract)
    moved <- !is.na(df$date_of_sample_collection) &
             (is.na(s_ex) | df$date_of_sample_collection != s_ex)
    l_kept <- if ("lab_analysis_date_extract" %in% names(df) && "lab_analysis_date" %in% names(df)) {
      l_ex <- .parse_date(df$lab_analysis_date_extract)
      (is.na(df$lab_analysis_date) & is.na(l_ex)) |
        (!is.na(df$lab_analysis_date) & !is.na(l_ex) & df$lab_analysis_date == l_ex)
    } else rep(TRUE, nrow(df))
    xtest <- moved & l_kept
  }
  specs <- dhis2_delay_specs()
  delays <- list()
  n_xtest <- integer(0)
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    if (!all(c(sp$from, sp$to) %in% names(df))) next
    a <- df[[sp$from]]; b <- df[[sp$to]]
    if (identical(sp$to, "lab_analysis_date")) {
      n_xtest[nm] <- sum(xtest & !is.na(b))
      b[xtest] <- NA
    }
    if (has_asof) b[!is.na(b) & b > .asof] <- NA   # as-of: secondary event not yet observable
    d <- as.integer(b - a)
    keep <- !is.na(a) & !is.na(b) & a >= analysis_start & a <= trunc_date &
            is.finite(d) & d >= 0L & d <= max_delay
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) next
    delays[[nm]] <- list(onset = a[keep], sample = b[keep], delay = d[keep])
  }
  if (length(n_xtest) && any(n_xtest > 0L))
    message(sprintf(paste0("[delays] lab-referenced delays exclude %d record(s) whose sample date ",
                           "is an earlier test's than their lab date (before windowing)."),
                    max(n_xtest)))
  list(specs = specs, delays = delays, analysis_start = analysis_start,
       trunc_date = trunc_date, obs_date = obs_date, max_sample = max_sample, n_cross_test = n_xtest,
       window_name = sprintf("analytical_%s", format(trunc_date, "%Y-%m-%d")),
       n_records = nrow(df), as_of = if (has_asof) .asof else NA)
}

.plot_censored_fits <- function(x, fits, panel_title, xlim_max = NULL) {
  x_c <- x[!is.na(x) & x >= 0L]
  if (length(x_c) == 0) return(ggplot() + labs(title = panel_title) + theme_void())
  if (is.null(xlim_max)) xlim_max <- min(max(x_c) + 2L, MAX_DELAY)
  xg  <- seq(0.01, xlim_max, by = 0.1)
  p   <- ggplot(data.frame(delay = x_c), aes(x = delay)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                   fill = "grey85", colour = "grey60", alpha = 0.75, linewidth = 0.3) +
    # coord_cartesian, NOT xlim(): xlim() is a SCALE limit, which censors the [-0.5, 0.5) zero
    # bin to NA and DROPS the same-day bar entirely (the share is in this panel's subtitle),
    # with only a "Removed rows containing missing values" warning. 04d uses coord_cartesian
    # for exactly this reason.
    coord_cartesian(xlim = c(-0.5, xlim_max)) +
    labs(title = panel_title,
         subtitle = sprintf("n=%d  (zeros=%d, %.1f%%)", length(x_c), sum(x_c == 0L),
                            sum(x_c == 0L) / length(x_c) * 100),
         x = "Delay (days)", y = "Density") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8, colour = "grey45"))
  if (is.null(fits) || length(fits) == 0L) return(p)
  aic_ord <- sort(sapply(fits, function(f) f$aic))
  top2    <- names(aic_ord)[seq_len(min(2L, length(aic_ord)))]
  cols2   <- c("#d35400", "#1a6faf")
  leg_lbl <- function(fam) sprintf("%s (AIC=%.1f, mean=%.2fd)", toupper(fam),
                                   fits[[fam]]$aic, .implied_mean(fam, fits[[fam]]$estimate))
  for (fam in top2) {
    dens <- .density_from_fit(fam, fits[[fam]]$estimate, xg)
    p <- p + geom_line(data = data.frame(x = xg, y = dens, label = leg_lbl(fam)),
                       aes(x = x, y = y, colour = label), linewidth = 1.1, na.rm = TRUE)
  }
  p + scale_colour_manual(values = setNames(cols2[seq_along(top2)], sapply(top2, leg_lbl)),
                          name = NULL) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
}

# ── C. Bayesian EpiDist (optional; corrects censoring AND right-truncation) ────
#
# POSTERIOR-DRAW REGISTRY. fit_epidist_model() summarises each fit to a one-row tibble
# (posterior median + 95% CrI of the mean) and then discards the draws, which is all the
# pipeline needs but leaves nothing to DRAW a fitted density — let alone a credible band —
# with. Refitting for the figure would double the Stan cost and, worse, could land on a
# different posterior from the one the pipeline actually used. So every fit deposits its
# thinned (mean, sd) draws here as a side effect, and 04d_delay_fit_figures.R renders the
# EXACT fits run_all.R already performed in step 1a. Storing draws (not the brms object)
# keeps the memory cost at a few hundred KB.
#
# KEYING CONTRACT. The key is delay|model_type|family and a put OVERWRITES silently. Every
# fit of a SUBSET of a delay's records therefore needs a key of its own: pass `draw_key =`
# through fit_epidist_both()/fit_epidist_model() (the convention is "<delay>__<stratum>").
# Reusing the plain delay key for a subset replaces the population posterior with the
# subset's, and because 04d is the only reader, the corruption surfaces only as a figure
# that disagrees with the published parameters.
.EPIDIST_DRAWS <- new.env(parent = emptyenv())
.epidist_draws_key <- function(delay, model_type, family)
  paste(delay, model_type, family, sep = "|")

#' Registry key for a delay fitted on a SUBSET of its records.
#'
#' Returns the bare delay name when `classes` is NULL (the whole population), and a
#' `<delay>__<stratum>` key otherwise. Derived from `classes` rather than from `stratum`, so
#' a caller that filters the records but forgets to name the stratum still gets a key of its
#' own instead of silently overwriting the population posterior.
.epidist_draw_key <- function(delay, classes = NULL, stratum = NULL) {
  if (is.null(classes) || !length(classes)) return(delay)
  lbl <- if (length(stratum) == 1L && !is.na(stratum) && nzchar(stratum)) stratum
         else paste(sort(unique(as.character(classes))), collapse = "+")
  paste0(delay, "__", lbl)
}

# `mean_v`/`sd_v` hold ONE value per posterior draw (fit_epidist_model() passes the draws of the
# first data row; the model is intercept-only, so every row carries the same draws).
# Deterministic thinning — seq(), never sample(): a side-channel that consumed the caller's
# random-number stream would make any later seeded draw in the session depend on whether a
# fit had been registered.
.epidist_draws_put <- function(delay, model_type, family, mean_v, sd_v, n,
                               frac_fully_observable, max_keep = 2000L) {
  ok <- is.finite(mean_v) & is.finite(sd_v) & mean_v > 0 & sd_v > 0
  mean_v <- mean_v[ok]; sd_v <- sd_v[ok]
  if (!length(mean_v)) return(invisible(NULL))
  if (length(mean_v) > max_keep) {
    idx    <- unique(round(seq(1, length(mean_v), length.out = max_keep)))
    mean_v <- mean_v[idx]; sd_v <- sd_v[idx]
  }
  assign(.epidist_draws_key(delay, model_type, family),
         list(delay = delay, model_type = model_type, family = family,
              mean = as.numeric(mean_v), sd = as.numeric(sd_v),
              n = n, frac_fully_observable = frac_fully_observable),
         envir = .EPIDIST_DRAWS)
  invisible(NULL)
}

#' Retrieve captured EpiDist posterior draws.
#' @param delay optional delay name to filter on (e.g. "onset_sample").
#' @return unnamed list of list(delay, model_type, family, mean, sd, n, frac_fully_observable).
epidist_draws_get <- function(delay = NULL) {
  out <- mget(ls(.EPIDIST_DRAWS), envir = .EPIDIST_DRAWS)
  if (!is.null(delay)) out <- out[vapply(out, function(x) identical(x$delay, delay), logical(1))]
  unname(out)
}
# epidist_draws_clear() REMOVED: nothing called it. The draw registry is a per-session
# side-channel from the fits to 04d's density/CDF bands; it is rebuilt every run and there is
# no point in a session where clearing it is correct.

#' Convergence verdict for one EpiDist fit (Vehtari et al. 2021 criteria: R-hat <= 1.01, bulk
#' and tail ESS >= 100 per chain, no divergent transitions). NA when a diagnostic is missing.
.epidist_converged <- function(max_rhat, min_ess_bulk, min_ess_tail, n_divergent,
                               chains = STAN_CHAINS) {
  v <- suppressWarnings(as.numeric(c(max_rhat, min_ess_bulk, min_ess_tail, n_divergent)))
  if (length(v) != 4L || !all(is.finite(v))) return(NA)
  v[1] <= 1.01 && v[2] >= 100 * chains && v[3] >= 100 * chains && v[4] == 0
}

# Mirrors the lab template's epidist path. The MARGINAL model right-truncates at
# obs_date; the NAIVE model does not. Returns posterior mean/sd (median + 95% CrI).
#' @param delay_name the delay this fit DESCRIBES; goes in the returned tibble's `delay`
#'   column, which every downstream writer keys on.
#' @param draw_key the key this fit is STORED under in the posterior-draw registry. It
#'   defaults to `delay_name`, which is right for a fit of the whole population. A fit of a
#'   SUBSET of that population (a case-classification stratum) must pass a distinct key:
#'   see the registry's keying contract above.
fit_epidist_model <- function(df_pairs, obs_date_val, model_type, family_fn,
                              family_name, delay_name, draw_key = delay_name) {
  tryCatch({
    df_input <- df_pairs %>%
      dplyr::transmute(pdate_lwr = .data$onset, sdate_lwr = .data$sample) %>%
      dplyr::filter(!is.na(pdate_lwr), !is.na(sdate_lwr),
                    as.numeric(sdate_lwr - pdate_lwr) >= 0) %>%
      dplyr::mutate(pdate_upr = pdate_lwr + 1L, sdate_upr = sdate_lwr + 1L,
                    obs_date  = pmax(obs_date_val, sdate_lwr + 1L))
    # naive lognormal/gamma cannot take exact zero delays -> keep strictly positive
    if (model_type == "naive" && family_name %in% c("lognormal", "gamma"))
      df_input <- df_input %>% dplyr::filter(as.numeric(sdate_lwr - pdate_lwr) > 0)
    if (nrow(df_input) < 5) return(NULL)
    # RIGHT-TRUNCATION EXPOSURE. The old metric was
    #   frac_c <- sum(sdate_lwr <= obs_date_val) / nrow(df_input)
    # which is 1.0 BY CONSTRUCTION: df_pairs is already filtered to samples on or before the
    # as-of date and obs_date_val IS that date. So the `< 0.30` guard could never fire and the
    # log always announced "100% complete" for a model whose entire purpose is correcting
    # right-truncation — the single most misleading line this file could print.
    #
    # Measure the real thing instead: the share of pairs whose PRIMARY event is far enough
    # behind the observation date that the full delay support (MAX_DELAY) could have been
    # observed. Everything else is subject to truncation. On the 2026-10-09 snapshot (earliest
    # sample dates, all windowed pairs) it is 0.55 / 0.53 / 0.55 / 0.51 for onset->sample,
    # onset->notification, onset->lab and sample->lab, and 0.53-0.61 for the three strata.
    frac_obs <- mean(as.numeric(obs_date_val - df_input$pdate_lwr) >= MAX_DELAY, na.rm = TRUE)
    if (!is.finite(frac_obs)) frac_obs <- 0
    # Threshold DELIBERATELY set far below that range so no fit that runs
    # today starts being skipped: this is a floor against pathological input (a window with
    # essentially no fully-observable pair), not a quality bar. The MARGINAL model is designed
    # to correct truncation, so a high truncated share is a reason to prefer it, not to skip it.
    if (frac_obs < 0.05) {
      cat(sprintf("  [%s %s %s] skipped — only %.1f%% of pairs are fully observable (>= %d d before obs_date)\n",
                  delay_name, model_type, family_name, frac_obs * 100, MAX_DELAY))
      return(NULL)
    }
    cat(sprintf("  Fitting EpiDist %s %s %s (n=%d, %.0f%% of pairs fully observable; the rest are right-truncated and corrected by the marginal model)...\n",
                delay_name, model_type, family_name, nrow(df_input), frac_obs * 100))
    ll_obj <- df_input %>%
      epidist::as_epidist_linelist_data(
        pdate_lwr = "pdate_lwr", pdate_upr = "pdate_upr",
        sdate_lwr = "sdate_lwr", sdate_upr = "sdate_upr", obs_date = "obs_date")
    model_obj <- switch(model_type,
      naive    = epidist::as_epidist_naive_model(ll_obj),
      marginal = ll_obj %>% epidist::as_epidist_aggregate_data() %>%
                   epidist::as_epidist_marginal_model())
    # SEEDED: without a fixed seed the posterior mean/sd wander between runs (observed:
    # onset->sample marginal mean 7.85 d vs 8.00 d on two identical runs), which would make
    # the reported delay (.load_dhis2_delay_params()) non-reproducible. 04b_epinowcast.R seeds
    # its fit the same way.
    fit <- epidist::epidist(model_obj, family = family_fn, chains = STAN_CHAINS,
                            iter = STAN_ITER, warmup = STAN_WARMUP, cores = STAN_CORES,
                            seed = get0("RANDOM_SEED", ifnotfound = 20260704L),
                            refresh = 0, silent = 2,
                            control = list(adapt_delta = 0.95, max_treedepth = 12))
    # CONVERGENCE, CHECKED AND REPORTED (.epidist_converged()). A failing fit is still returned
    # -- the posterior is what it is -- but it is flagged loudly, and its diagnostics travel with
    # it: the returned tibble carries them, write_onset_sample_long() persists them as
    # epidist_* rows of the delay params file, and 04d joins them to dhis2_delay_fit_summary.csv.
    .cv <- tryCatch({
      sm <- posterior::summarise_draws(posterior::as_draws_df(fit), "rhat", "ess_bulk", "ess_tail")
      sm <- sm[!grepl("^lp__$|^lprior$", sm$variable), , drop = FALSE]
      np <- brms::nuts_params(fit)
      c(max_rhat = max(sm$rhat, na.rm = TRUE), min_ess_bulk = min(sm$ess_bulk, na.rm = TRUE),
        min_ess_tail = min(sm$ess_tail, na.rm = TRUE),
        n_divergent = sum(np$Value[np$Parameter == "divergent__"]))
    }, error = function(e) c(max_rhat = NA_real_, min_ess_bulk = NA_real_,
                             min_ess_tail = NA_real_, n_divergent = NA_real_))
    .ess_min <- 100 * STAN_CHAINS
    if (!isTRUE(.epidist_converged(.cv[["max_rhat"]], .cv[["min_ess_bulk"]],
                                   .cv[["min_ess_tail"]], .cv[["n_divergent"]])))
      warning(sprintf(paste0("[epidist] %s %s %s: convergence NOT established (max R-hat %.3f, ",
                             "min bulk/tail ESS %.0f/%.0f against %d, %s divergent transition(s))."),
                      delay_name, model_type, family_name, .cv[["max_rhat"]],
                      .cv[["min_ess_bulk"]], .cv[["min_ess_tail"]], .ess_min,
                      format(.cv[["n_divergent"]])), call. = FALSE)
    samps <- epidist::predict_delay_parameters(fit) %>% epidist::add_mean_sd()
    # ONE ROW PER POSTERIOR DRAW. predict_delay_parameters() returns draws x data rows (draw
    # varying fastest). The model is intercept-only (epidist()'s default mu ~ 1), so every data
    # row carries the same draws: keep the first row's. Summarising or thinning the stacked
    # vector instead weights draws unevenly, and an evenly spaced pick over draws x rows can
    # land on only a few distinct draws when the row count shares a factor with the spacing.
    if ("index" %in% names(samps)) samps <- samps[samps$index == min(samps$index), , drop = FALSE]
    # Side effect only — deposit the draws for 04d's density/CDF bands. Never alters the
    # returned tibble, so every existing consumer of this function is untouched.
    .epidist_draws_put(draw_key, model_type, family_name,
                       samps$mean, samps$sd, n = nrow(df_input), frac_fully_observable = frac_obs)
    tibble(delay = delay_name, model_type = model_type, family = family_name,
           n = nrow(df_input), frac_fully_observable = round(frac_obs, 3),
           mean_post_med = round(stats::median(samps$mean, na.rm = TRUE), 3),
           mean_post_lo  = round(stats::quantile(samps$mean, 0.025, na.rm = TRUE), 3),
           mean_post_hi  = round(stats::quantile(samps$mean, 0.975, na.rm = TRUE), 3),
           sd_post_med   = round(stats::median(samps$sd, na.rm = TRUE), 3),
           max_rhat = round(.cv[["max_rhat"]], 4), min_ess_bulk = round(.cv[["min_ess_bulk"]]),
           min_ess_tail = round(.cv[["min_ess_tail"]]), n_divergent = .cv[["n_divergent"]])
  }, error = function(e) {
    # A warning(), not only a cat(). This handler swallowed an "object 'frac_c' not found"
    # for every single EpiDist fit — raised AFTER Stan had finished, so four fits were burned
    # per run and the caller quietly fell back to whatever epidist_* rows were already on
    # disk. Nothing in the run summary said so. A cat() into a long log is not a failure
    # signal for the estimator that sets the reported onset->sample delay and the fallbacks of
    # effective_onset_sample_delay().
    cat(sprintf("  [%s %s %s] FAILED: %s\n", delay_name, model_type, family_name, e$message))
    warning(sprintf("[epidist] %s %s %s FAILED: %s", delay_name, model_type, family_name,
                    conditionMessage(e)), call. = FALSE)
    NULL
  })
}

fit_epidist_both <- function(df_pairs, delay_name, obs_date_val, draw_key = delay_name) {
  fams <- list(lognormal = quote(brms::lognormal()), gamma = quote(stats::Gamma(link = "log")))
  purrr::map_dfr(names(fams), function(fam) {
    fam_fn <- eval(fams[[fam]])
    purrr::map_dfr(c("naive", "marginal"), function(mt)
      fit_epidist_model(df_pairs, obs_date_val, mt, fam_fn, fam, delay_name,
                        draw_key = draw_key))
  })
}

#' Fit the onset->sample delay SEPARATELY by case classification.
#'
#' WHY. The REPORTED onset->sample delay is that of CONFIRMED cases (ONSET_IMPUTE_STRATUM,
#' 00_config.R; since 2026-10-15 it no longer drives the onset imputation, which draws
#' same-class donor intervals in 01_data_prep.R). A fit pooled over every classification is
#' mostly TEST-NEGATIVE and partly UNADJUDICATED records: on the 2026-10-09 snapshot (earliest
#' sample dates) the 19,968 windowed onset->sample pairs are 10,262 `not_a_case`, 2,793 with no
#' final classification, 6,736 confirmed and 177 suspected/probable, with raw means of 6.00 d,
#' 3.53 d and 6.21 d for the first three groups; the interval-censored MLE gives 5.72 d pooled
#' against 6.21 d confirmed-only on that frame. Name the estimator with any number quoted: the
#' MLE and EpiDist figures for the same contrast differ, and have been confused before.
#'
#' Each stratum is fitted INDEPENDENTLY, with its own double-interval censoring and its own
#' right-truncation correction. That is not a nicety: a stratum sampled more slowly is also
#' more right-truncated at any as-of date, and a pooled correction would distort the very
#' contrast this function exists to measure.
#'
#' `suspected` pools `suspected_case` with `probable_case` (25 windowed probable_case pairs on
#' the frame above): too few for a usable interval on their own.
#'
#' @return named list of epidist results (the shape estimate_dhis2_onset_sample_epidist()
#'   returns), one per stratum, with failures dropped. Empty list if epidist is unavailable.
estimate_onset_sample_strata <- function(ll, analysis_date = ANALYSIS_DATE,
                                         outbreak_start = OUTBREAK_START,
                                         test_days = TEST_DAYS,
                                         max_delay = MAX_DELAY) {
  if (!isTRUE(get0("RUN_EPIDIST", ifnotfound = FALSE))) return(list())
  if (!"final_mve_case_classification" %in% names(ll)) {
    warning("[04c_dhis2] no final_mve_case_classification column; delay strata not fitted.",
            call. = FALSE)
    return(list())
  }
  conf <- get0("CONFIRMED_STATUS", ifnotfound = "confirmed_case")
  spec <- list(confirmed = conf,
               not_a_case = "not_a_case",
               suspected  = c("suspected_case", "probable_case"))
  out <- list()
  for (nm in names(spec)) {
    r <- tryCatch(estimate_dhis2_onset_sample_epidist(
           ll, analysis_date = analysis_date, outbreak_start = outbreak_start,
           test_days = test_days, max_delay = max_delay, classes = spec[[nm]],
           stratum = nm),
         error = function(e) { message("[04c_dhis2] stratum '", nm, "' failed: ",
                                       conditionMessage(e)); NULL })
    if (!is.null(r) && is.finite(r$mean) && r$mean > 0) {
      out[[nm]] <- r
      message(sprintf("[04c_dhis2] delay stratum %-11s %s mean %.2f d (sd %.2f, n=%s)",
                      nm, r$family, r$mean, r$sd, format(r$n)))
    } else {
      message(sprintf("[04c_dhis2] delay stratum %-11s NOT fitted", nm))
    }
  }
  out
}

# ============================================================================
# CORE: estimate the DHIS2 onset->sample delay (windowed, interval-censored)
# ============================================================================

#' Estimate the DHIS2 onset->sample delay by windowed interval-censored MLE.
#'
#' @param ll   DHIS2 line list with `date_of_symptom_onset` and
#'   `date_of_sample_collection` (Date, or coercible via .parse_date).
#' @param analysis_date snapshot / extraction date (for window derivation).
#' @param outbreak_start window start (default OUTBREAK_START).
#' @param test_days right-truncation buffer (default TEST_DAYS = 5).
#' @param max_delay drop delays > this many days (typos/outliers).
#' @param confirmed_only if TRUE, restrict to confirmed cases (default FALSE:
#'   all complete pairs give the largest sample; note the delay is NOT invariant to
#'   classification -- confirmed cases are swabbed later than test-negatives -- which is why
#'   estimate_onset_sample_strata() also fits each classification separately).
#' @return list(best_family, params, implied_mean, rate, n_fit, window, trunc_date,
#'   delays_windowed [vector for empirical bootstrap], cens_table, cens_fits) or
#'   NULL if too few pairs.
estimate_dhis2_onset_sample_delay <- function(ll, analysis_date = ANALYSIS_DATE,
                                              outbreak_start = OUTBREAK_START,
                                              test_days = TEST_DAYS,
                                              max_delay = MAX_DELAY,
                                              confirmed_only = FALSE) {
  on <- .parse_date(ll[["date_of_symptom_onset"]])
  sa <- .parse_date(ll[["date_of_sample_collection"]])
  keep <- !is.na(on) & !is.na(sa)
  if (confirmed_only && "final_mve_case_classification" %in% names(ll))
    keep <- keep & (ll[["final_mve_case_classification"]] %in% CONFIRMED_STATUS)
  # As-of consistency: a sample collected AFTER the analysis (as-of) date is not observable
  # at that moment, so drop it before fitting. This makes `analysis_date` actually bound the
  # window (it was previously ignored — trunc_date came purely from max(sample)), preventing
  # future-data leakage on a back-dated re-run and neutralising a future-dated sample-date
  # typo. No-op when the line list is already snapshotted to the as-of date (all sa <= as-of).
  .asof <- .parse_date(analysis_date)
  if (length(.asof) == 1L && !is.na(.asof)) keep <- keep & (sa <= .asof)
  on <- on[keep]; sa <- sa[keep]
  if (length(on) < 5L) return(NULL)
  # Dynamic window: onset in [outbreak_start, TRUNC_DATE = max(sample) - test_days], where
  # max(sample) is now bounded at the as-of date by the filter above.
  trunc_date <- max(sa, na.rm = TRUE) - test_days
  ob_floor   <- lubridate::floor_date(outbreak_start, unit = "week", week_start = 1L)
  in_win <- on >= ob_floor & on <= trunc_date
  delay  <- as.integer(sa - on)
  ok     <- in_win & is.finite(delay) & delay >= 0L & delay <= max_delay
  d      <- delay[ok]
  if (length(d) < 5L) return(NULL)
  fits <- .fit_all_censored(d, "onset->sample")
  best <- .cens_best_fam(fits)
  if (is.na(best)) return(NULL)
  p    <- fits[[best]]$estimate
  window_name <- sprintf("analytical_%s", format(trunc_date, "%Y-%m-%d"))
  list(
    best_family    = best,
    params         = p,
    implied_mean   = .implied_mean(best, p),
    # Exponential-rate SUMMARY (1/mean) for pipeline consumers that key on a rate,
    # regardless of the AIC-best family — so every consumer's Exp fallback and the
    # delay CDF stay well-defined. The full best-family params are also returned.
    rate           = 1 / .implied_mean(best, p),
    n_fit          = length(d),
    window         = window_name,
    trunc_date     = trunc_date,
    delays_windowed = d,               # windowed, plausibility-filtered delay vector
    cens_table     = .build_cens_tbl(fits, "onset_sample", window_name, length(d)),
    cens_fits      = fits)
}

#' Bayesian EpiDist MARGINAL onset->sample delay for the PIPELINE path (truncation +
#' double-interval-censoring corrected). estimate_dhis2_onset_sample_delay() above returns
#' the interval-censored MLE only; this companion fits the EpiDist marginal model on the
#' SAME windowed complete pairs and returns list(family, mean, sd, mean_lo, mean_hi, n,
#' max_rhat, min_ess_bulk, min_ess_tail, n_divergent) — the shape write_onset_sample_long()
#' writes (`epidist=`, and each element of `strata=`) and .load_dhis2_delay_params() PREFERS —
#' so run_all.R writes the truncation-corrected delay to the params file, not only the
#' standalone Rscript path. Mirrors the MAIN block's §1.2 extraction (prefer the gamma
#' marginal, else the first marginal row; map lognormal->lnorm).
#' Returns NULL when RUN_EPIDIST is off, the `epidist` package is unavailable (RUN_EPIDIST
#' already folds in .HAVE_EPIDIST), there are too few pairs, or the fit fails; the writer then
#' keeps an EpiDist fit of the SAME records already on disk, else selects the censored MLE,
#' so the default remains non-fatal.
#' @param classes optional character vector of `final_mve_case_classification` values to
#'   restrict the fit to. NULL (the default) pools every classification, which is what the
#'   pipeline did until 2026-09-22. Pass CONFIRMED_STATUS to fit the delay of the population
#'   whose delay is reported (confirmed cases) -- see estimate_onset_sample_strata() for why
#'   that matters.
#' @param stratum optional display name for the `classes` subset (e.g. "confirmed"), used
#'   ONLY to key this fit's posterior draws. Ignored when `classes` is NULL. When `classes`
#'   is given but `stratum` is not, the key falls back to the sorted class values, so a
#'   caller cannot accidentally reuse the pooled key by forgetting the label.
estimate_dhis2_onset_sample_epidist <- function(ll, analysis_date = ANALYSIS_DATE,
                                                outbreak_start = OUTBREAK_START,
                                                test_days = TEST_DAYS,
                                                max_delay = MAX_DELAY,
                                                classes = NULL,
                                                stratum = NULL) {
  if (!isTRUE(get0("RUN_EPIDIST", ifnotfound = FALSE)) ||
      !exists("fit_epidist_both", mode = "function")) return(NULL)
  on <- .parse_date(ll[["date_of_symptom_onset"]])
  sa <- .parse_date(ll[["date_of_sample_collection"]])
  keep <- !is.na(on) & !is.na(sa)
  # CLASS FILTER FIRST, before the window and truncation logic, so a stratum's window is
  # derived from its OWN sample dates. Deriving it from the pooled maximum would hand a small
  # stratum a truncation reference it never reaches.
  if (!is.null(classes)) {
    if (!"final_mve_case_classification" %in% names(ll)) {
      warning("[04c_dhis2] classes= requested but final_mve_case_classification is absent; ",
              "the fit would silently pool every classification. Returning NULL.", call. = FALSE)
      return(NULL)
    }
    keep <- keep & (ll[["final_mve_case_classification"]] %in% classes)
  }
  # Right-truncation reference = the as-of / observation date (bounds what is observable);
  # also drop any sample dated after it, exactly like estimate_dhis2_onset_sample_delay().
  .asof    <- .parse_date(analysis_date)
  has_asof <- length(.asof) == 1L && !is.na(.asof)
  if (has_asof) keep <- keep & (sa <= .asof)
  on <- on[keep]; sa <- sa[keep]
  if (length(on) < 5L) return(NULL)
  obs_date   <- if (has_asof) .asof else suppressWarnings(max(sa, na.rm = TRUE))
  # SAME window as the censored-MLE fit: onset in [outbreak floor, max(sample) - test_days].
  trunc_date <- suppressWarnings(max(sa, na.rm = TRUE)) - test_days
  ob_floor   <- lubridate::floor_date(outbreak_start, unit = "week", week_start = 1L)
  delay <- as.integer(sa - on)
  ok    <- on >= ob_floor & on <= trunc_date & is.finite(delay) & delay >= 0L & delay <= max_delay
  if (sum(ok, na.rm = TRUE) < 5L) return(NULL)
  tryCatch({
    # REGISTRY KEY. The returned tibble keeps delay = "onset_sample" (every writer keys on
    # it), but a stratum fit must be STORED under its own key. Until 2026-10-02 it was not:
    # every stratum deposited its draws under the pooled "onset_sample|<model>|<family>" key,
    # so the last one fitted (suspected/probable, the smallest stratum) overwrote the pooled
    # posterior, and 04d_delay_fit_figures.R -- the registry's only reader -- drew that
    # stratum as if it were the pooled fit written to the params file. Derived from `classes`,
    # not from `stratum`, so the qualification cannot be lost by omitting the label.
    .draw_key <- .epidist_draw_key("onset_sample", classes, stratum)
    tbl <- fit_epidist_both(tibble::tibble(onset = on[ok], sample = sa[ok]),
                            "onset_sample", obs_date, draw_key = .draw_key)
    if (is.null(tbl) || !nrow(tbl)) return(NULL)
    mm <- tbl[tbl$model_type == "marginal" & is.finite(tbl$mean_post_med), , drop = FALSE]
    if (!nrow(mm)) return(NULL)
    gi   <- match("gamma", mm$family)
    pick <- if (!is.na(gi)) mm[gi, , drop = FALSE] else mm[1, , drop = FALSE]
    efam <- if (identical(as.character(pick$family[1]), "lognormal")) "lnorm"
            else as.character(pick$family[1])
    # mean_lo/mean_hi are the 95% posterior interval of the MEAN, carried so a consumer can
    # treat the delay as uncertain rather than known exactly (04d's panel E draws them).
    # max_rhat/min_ess_*/n_divergent are the picked fit's convergence diagnostics, persisted by
    # write_onset_sample_long() so a non-converged fit is visible in the params file itself.
    out <- list(family = efam, mean = pick$mean_post_med[1], sd = pick$sd_post_med[1],
                mean_lo = pick$mean_post_lo[1], mean_hi = pick$mean_post_hi[1],
                n = pick$n[1],
                max_rhat = pick$max_rhat[1], min_ess_bulk = pick$min_ess_bulk[1],
                min_ess_tail = pick$min_ess_tail[1], n_divergent = pick$n_divergent[1])
    # Full naive+marginal x lognormal+gamma table carried along as an attribute: run_all.R
    # hands it to 04d, which joins its convergence diagnostics to dhis2_delay_fit_summary.csv
    # for ALL four sub-fits (the curves themselves come from the posterior-draw registry).
    # An attribute keeps write_onset_sample_long(epidist=)'s expected list shape.
    attr(out, "epidist_table") <- tbl
    out
  }, error = function(e) {
    message("[04c_dhis2] EpiDist marginal fit failed (non-fatal): ", conditionMessage(e)); NULL
  })
}

#' Native parameters of an EpiDist fit from its (mean, sd), by the SAME conversion
#' .load_dhis2_delay_params() (00_config.R) applies: method of moments for gamma, log-moments
#' for lnorm. numeric(0) when the family has no conversion or the moments are unusable, so a
#' caller never substitutes another family's parameters.
.epidist_native_params <- function(family, m, s) {
  m <- suppressWarnings(as.numeric(m)); s <- suppressWarnings(as.numeric(s))
  if (length(family) != 1L || is.na(family) || length(m) != 1L || length(s) != 1L ||
      !is.finite(m) || !is.finite(s) || m <= 0 || s <= 0) return(numeric(0))
  switch(as.character(family),
    gamma = c(shape = (m / s)^2, rate = m / s^2),
    lnorm = { s2 <- log(1 + (s / m)^2); c(meanlog = log(m) - s2 / 2, sdlog = sqrt(s2)) },
    numeric(0))
}

#' Convergence rows of one EpiDist fit for the params file: epidist_max_rhat, _min_ess_bulk,
#' _min_ess_tail, _n_divergent and _converged (.epidist_converged()), each with suffix `sfx`
#' ("" for the pooled fit, "__<stratum>" for a stratum). NULL when the fit carries no
#' diagnostics. The loader reads none of these rows; they are the audit record.
.epidist_convergence_rows <- function(e, sfx = "") {
  q <- c("max_rhat", "min_ess_bulk", "min_ess_tail", "n_divergent")
  if (!is.list(e) || !all(q %in% names(e))) return(NULL)
  v <- vapply(q, function(k) {
    x <- suppressWarnings(as.numeric(e[[k]])); if (length(x) == 1L) x else NA_real_ }, numeric(1))
  tibble::tibble(delay = "onset_to_sample",
                 quantity = paste0("epidist_", c(q, "converged"), sfx),
                 value = c(as.character(round(v, 4)),
                           as.character(.epidist_converged(v[[1]], v[[2]], v[[3]], v[[4]]))))
}

# Write the SELECTED onset->sample family to the long (delay,quantity,value) CSV
# the pipeline consumes (mirrors data/cfr_reference/onset_to_sample_delay_params.csv).
#' @param strata optional named list of per-classification epidist fits, as returned by
#'   estimate_onset_sample_strata(). Written as `<quantity>__<stratum>` rows, e.g.
#'   `epidist_mean__confirmed`. The separator is a DOUBLE underscore so a stratum name can
#'   never be confused with an existing suffix (`epidist_mean_lo` is a pooled quantity, not
#'   the "lo" stratum). The pooled rows are untouched, so every existing reader keeps working;
#'   .load_dhis2_delay_params(stratum=) is what reads these.
write_onset_sample_long <- function(fit, path, source_label, epidist = NULL, strata = NULL) {
  p <- fit$params; fam <- fit$best_family
  rows <- tibble::tribble(
    ~delay,            ~quantity,          ~value,
    "onset_to_sample", "family",           fam,
    "onset_to_sample", "source",           source_label,
    "onset_to_sample", "window",           fit$window,
    "onset_to_sample", "estimator",        "interval_censored_mle",
    "onset_to_sample", "n_fit",            as.character(fit$n_fit),
    "onset_to_sample", "implied_mean_fit", as.character(round(fit$implied_mean, 3)),
    "onset_to_sample", "rate",             as.character(round(fit$rate, 4)))
  # Truncation-corrected EpiDist MARGINAL params (review §1.2): when the Bayesian
  # marginal model has been fit (RUN_EPIDIST=TRUE), append its (family, mean, sd) so
  # .load_dhis2_delay_params() can PREFER the right-truncation + double-interval-censored
  # estimate over the windowed interval-censored MLE (which only mitigates truncation by
  # dropping the recent tail), plus its convergence diagnostics. No-op when `epidist` is NULL.
  epi_rows <- NULL
  if (!is.null(epidist) && !is.null(epidist$mean) && is.finite(epidist$mean) && epidist$mean > 0) {
    epi_rows <- dplyr::bind_rows(tibble::tribble(
      ~delay,            ~quantity,         ~value,
      "onset_to_sample", "epidist_family",  as.character(epidist$family),
      "onset_to_sample", "epidist_mean",    as.character(round(epidist$mean, 3)),
      "onset_to_sample", "epidist_sd",      as.character(round(epidist$sd, 3)),
      "onset_to_sample", "epidist_mean_lo", as.character(if (is.null(epidist$mean_lo)) NA else round(epidist$mean_lo, 3)),
      "onset_to_sample", "epidist_mean_hi", as.character(if (is.null(epidist$mean_hi)) NA else round(epidist$mean_hi, 3)),
      "onset_to_sample", "epidist_n",       as.character(if (is.null(epidist$n)) NA else epidist$n),
      "onset_to_sample", "epidist_estimator", "epidist_marginal_truncation_corrected"),
      .epidist_convergence_rows(epidist))
  }
  # Append the best family's native parameters so a consumer can reconstruct the
  # full fitted distribution (not only the Exp-rate summary). PREFIX with `param_`
  # so a family whose native parameter is itself called `rate` (gamma, exp) does
  # NOT collide with the Exp-rate summary row above (which would make deframe()
  # ambiguous). e.g. gamma -> param_shape, param_rate.
  par_rows <- purrr::imap_dfr(as.list(round(p, 5)), function(v, nm)
    tibble(delay = "onset_to_sample", quantity = paste0("param_", nm),
           value = as.character(unname(v))))
  # SELECTED-* rows: state unambiguously which numbers the pipeline actually uses.
  # Without them this file is internally contradictory. The `family`/`estimator`/`rate`/
  # `param_*`/`implied_mean_fit` rows above describe the interval-censored MLE, but
  # .load_dhis2_delay_params() (00_config.R) PREFERS the epidist_* rows and re-derives the
  # family's native parameters from the marginal mean and sd -- so without these rows the file
  # asserted two different parameterisations and an `estimator` that is false whenever
  # EpiDist ran. These rows mirror the loader's own preference order AND its moment
  # conversion (.epidist_native_params(), on the rounded mean and sd the loader reads), so
  # `selected_*` is always what the loader resolves for the pooled fit.
  .use_epi <- !is.null(epi_rows)
  .epi_supplied <- .use_epi
  if (.use_epi) {
    sel_fam  <- as.character(epidist$family)
    sel_mean <- round(epidist$mean, 3)
    sel_sd   <- suppressWarnings(as.numeric(epidist$sd))
    sel_sd   <- if (length(sel_sd) == 1L && is.finite(sel_sd)) round(sel_sd, 3) else NA_real_
    sel_par  <- .epidist_native_params(sel_fam, sel_mean, sel_sd)
  } else {
    sel_fam <- fam; sel_mean <- fit$implied_mean; sel_sd <- NA_real_; sel_par <- p
  }
  # NO EpiDist FIT SUPPLIED (RUN_EPIDIST off, or the fit failed), but the file on disk may
  # carry an earlier one. It is KEPT only when (a) EpiDist was requested (RUN_EPIDIST, which
  # in this file also requires the packages) and (b) it is a fit of the SAME records: the same
  # snapshot and window, and the same deterministic interval-censored MLE (n_fit and
  # implied_mean_fit). Then a failed refit does not downgrade the delay to the MLE. Otherwise
  # every epidist_* row is DROPPED and the MLE is selected: rows kept from another snapshot (or
  # another sample-date rule) would be written next to THIS snapshot's `source` row and pass
  # the loader's "fit must be this snapshot's" guard, and RUN_EPIDIST=FALSE would never fall
  # back to the censored MLE as 00_config.R documents.
  if (!.use_epi && file.exists(path)) {
    .prev <- tryCatch(readr::read_csv(path, col_types = readr::cols(.default = "c"),
                                      show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(.prev) && all(c("delay", "quantity", "value") %in% names(.prev))) {
      .pv <- function(q) { v <- .prev$value[.prev$quantity == q]; if (length(v)) v[1] else NA_character_ }
      .keep <- .prev[grepl("^epidist_", .prev$quantity), c("delay", "quantity", "value"), drop = FALSE]
      # A stratum this call REWRITES (strata_rows below) must not also be preserved: the file
      # would then hold two rows per quantity, which the loader refuses.
      .fresh <- names(Filter(function(e) !is.null(e) && is.finite(e$mean) && e$mean > 0,
                             if (is.null(strata)) list() else strata))
      if (length(.fresh))
        .keep <- .keep[!grepl(paste0("__(", paste(.fresh, collapse = "|"), ")$"), .keep$quantity), ,
                       drop = FALSE]
      if (nrow(.keep)) {
        .same <- identical(.pv("source"), as.character(source_label)) &&
                 identical(.pv("window"), as.character(fit$window)) &&
                 identical(.pv("n_fit"), as.character(fit$n_fit)) &&
                 identical(.pv("implied_mean_fit"), as.character(round(fit$implied_mean, 3)))
        if (isTRUE(get0("RUN_EPIDIST", ifnotfound = FALSE)) && .same) {
          warning(sprintf(paste0("[04c_dhis2] no EpiDist fit was supplied, but %s already carries ",
                                 "one of the SAME records (snapshot %s, window %s, n = %s). ",
                                 "PRESERVING its epidist_* rows rather than downgrading the delay ",
                                 "to the interval-censored MLE."),
                          basename(path), source_label, fit$window, format(fit$n_fit)),
                  call. = FALSE, immediate. = TRUE)
          epi_rows <- .keep
          .kv  <- function(q) { v <- .keep$value[.keep$quantity == q]; if (length(v)) v[1] else NA_character_ }
          .k_m <- suppressWarnings(as.numeric(.kv("epidist_mean")))
          .k_f <- .kv("epidist_family")
          # Only a POOLED EpiDist fit changes the selection; a file that kept stratum rows alone
          # still resolves the MLE for the pooled read, and selected_* must say so.
          if (!is.na(.k_f) && is.finite(.k_m) && .k_m > 0) {
            .use_epi <- TRUE
            sel_fam  <- .k_f
            sel_mean <- .k_m
            sel_sd   <- suppressWarnings(as.numeric(.kv("epidist_sd")))
            sel_par  <- .epidist_native_params(sel_fam, sel_mean, sel_sd)
          }
        } else {
          .why <- if (!isTRUE(get0("RUN_EPIDIST", ifnotfound = FALSE)))
                    "RUN_EPIDIST is off (or the epidist packages are absent)"
                  else sprintf("it describes other records (snapshot %s, window %s, n = %s, MLE mean %s)",
                               .pv("source"), .pv("window"), .pv("n_fit"), .pv("implied_mean_fit"))
          warning(sprintf(paste0("[04c_dhis2] no EpiDist fit was supplied and the one in %s is ",
                                 "DROPPED because %s; the interval-censored MLE is the selected ",
                                 "delay."), basename(path), .why),
                  call. = FALSE, immediate. = TRUE)
        }
      }
    }
  }
  sel_rows <- tibble::tibble(
    delay = "onset_to_sample",
    # if(): paste0() of a prefix and a zero-length vector returns the bare prefix, not
    # character(0), which would leave a value-less "selected_param_" row.
    quantity = c("selected_estimator", "selected_family", "selected_mean", "selected_sd",
                 "selected_rate",
                 if (length(sel_par)) paste0("selected_param_", names(sel_par))),
    value = c(if (.use_epi) "epidist_marginal_truncation_corrected" else "interval_censored_mle",
              sel_fam,
              as.character(round(sel_mean, 3)),
              as.character(if (is.finite(sel_sd)) round(sel_sd, 3) else NA),
              as.character(round(1 / sel_mean, 5)),
              as.character(round(unname(sel_par), 5))))
  # A STRATUM WHOSE REFIT FAILED while the pooled fit succeeded (the branch above runs only when
  # the POOLED fit is missing): its earlier rows are kept under the same rule -- EpiDist
  # requested and a fit of the SAME records -- so a failed confirmed-stratum refit does not
  # downgrade the reported delay to the pooled fit; otherwise they are dropped.
  if (.epi_supplied && file.exists(path)) {
    .prev <- tryCatch(readr::read_csv(path, col_types = readr::cols(.default = "c"),
                                      show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(.prev) && all(c("delay", "quantity", "value") %in% names(.prev))) {
      .pv <- function(q) { v <- .prev$value[.prev$quantity == q]; if (length(v)) v[1] else NA_character_ }
      .fresh <- names(Filter(function(e) !is.null(e) && is.finite(e$mean) && e$mean > 0,
                             if (is.null(strata)) list() else strata))
      .ks <- .prev[grepl("^epidist_.*__[^_].*$", .prev$quantity), c("delay", "quantity", "value"),
                   drop = FALSE]
      if (length(.fresh))
        .ks <- .ks[!grepl(paste0("__(", paste(.fresh, collapse = "|"), ")$"), .ks$quantity), , drop = FALSE]
      if (nrow(.ks)) {
        .same <- identical(.pv("source"), as.character(source_label)) &&
                 identical(.pv("window"), as.character(fit$window)) &&
                 identical(.pv("n_fit"), as.character(fit$n_fit)) &&
                 identical(.pv("implied_mean_fit"), as.character(round(fit$implied_mean, 3)))
        .st <- unique(sub("^.*__", "", .ks$quantity))
        if (isTRUE(get0("RUN_EPIDIST", ifnotfound = FALSE)) && .same) {
          warning(sprintf(paste0("[04c_dhis2] stratum fit(s) %s were not refitted; %s already carries ",
                                 "fits of the SAME records, which are PRESERVED."),
                          paste(.st, collapse = ", "), basename(path)), call. = FALSE, immediate. = TRUE)
          epi_rows <- dplyr::bind_rows(epi_rows, .ks)
        } else
          warning(sprintf(paste0("[04c_dhis2] stratum fit(s) %s were not refitted and the ones in %s ",
                                 "describe other records (or EpiDist is off); they are DROPPED."),
                          paste(.st, collapse = ", "), basename(path)), call. = FALSE, immediate. = TRUE)
      }
    }
  }
  # PER-CLASSIFICATION rows. Additive only: nothing above is altered, so a reader that does
  # not know about strata resolves exactly what it resolved before. A stratum not refitted by
  # this call is kept or dropped under the same-records rule: with the pooled rows when the
  # pooled fit is also missing (the first branch), on its own otherwise (the branch above).
  strata_rows <- NULL
  if (!is.null(strata) && length(strata)) {
    strata_rows <- dplyr::bind_rows(lapply(names(strata), function(nm) {
      e <- strata[[nm]]
      if (is.null(e) || !is.finite(e$mean) || e$mean <= 0) return(NULL)
      .n <- function(x) as.character(if (is.null(x) || !is.finite(x)) NA else round(x, 3))
      dplyr::bind_rows(tibble::tibble(
        delay = "onset_to_sample",
        quantity = paste0(c("epidist_family", "epidist_mean", "epidist_sd",
                            "epidist_mean_lo", "epidist_mean_hi", "epidist_n"), "__", nm),
        value = c(as.character(e$family), .n(e$mean), .n(e$sd), .n(e$mean_lo), .n(e$mean_hi),
                  as.character(if (is.null(e$n)) NA else e$n))),
        .epidist_convergence_rows(e, sfx = paste0("__", nm)))
    }))
  }
  out <- dplyr::bind_rows(rows, par_rows, epi_rows, sel_rows, strata_rows)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(out, path)
  out
}

# ============================================================================
# MAIN — run only when executed as a script (sourcing this file for its functions --
# run_all.R step 1a, 04d_delay_fit_figures.R, the tests -- must NOT trigger the fit).
# `sys.nframe() == 0L` is TRUE only for `Rscript this.R` / `R CMD BATCH`, FALSE when source()d.
# ============================================================================
if (sys.nframe() == 0L || isTRUE(getOption("dhis2_delay_run_main"))) {
  # This block writes the canonical delay files directly; a data-definition sensitivity must not.
  stop_if_data_definition_sensitivity("04c_dhis2_delay_windows.R")
  cat(strrep("=", 78), "\n")
  cat("04c_dhis2_delay_windows.R — DHIS2 line-list delay analysis\n")
  cat(sprintf("Interval-censored MLE%s | %s\n",
              if (RUN_EPIDIST) " + EpiDist (Bayesian)" else "", Sys.time()))
  cat(strrep("=", 78), "\n\n")

  # ── Load the DHIS2 line list (raw complete pairs; NO onset imputation) ───────
  if (!file.exists(LINELIST_JSON))
    stop("[04c_dhis2] latest.json not found: ", LINELIST_JSON)
  .meta <- jsonlite::fromJSON(LINELIST_JSON)
  ll_csv <- file.path(LINELIST_DIR, .meta$folder, "dhis2_processed_linelist.csv")
  if (!file.exists(ll_csv)) stop("[04c_dhis2] line list not found: ", ll_csv)
  cat("Loading:", ll_csv, "\n")
  df_all <- readr::read_csv(ll_csv, col_types = readr::cols(.default = "c"),
                            show_col_types = FALSE, na = c("", "NA", "N/A"))
  # Same diagnostic sample date -- and, for confirmed records, the confirming test's lab date --
  # as the pipeline (SAMPLE_DATE_SOURCE, 00_config.R). Unconditional: 00_config.R is sourced at
  # the top of this file, so the function always exists, and skipping it would fit the delays
  # on the extract's repeat-test dates.
  df_all <- apply_earliest_sample_dates(df_all, label = "04c standalone")
  cat(sprintf("  Loaded: %d records\n", nrow(df_all)))

  date_cols <- c("date_of_symptom_onset", "date_of_sample_collection",
                 "date_of_notification", "reporting_date", "lab_analysis_date")
  for (dc in intersect(date_cols, names(df_all))) df_all[[dc]] <- .parse_date(df_all[[dc]])

  # AS-OF BOUND, by default, even standalone. This used to pass analysis_date = NULL, which
  # deliberately disables the bound: max_sample, TRUNC_DATE and OBS_DATE were then taken over
  # ALL rows, including any future-dated onset or sample typo (this file's own comment at the
  # builder concedes such typos exist). The os_fit call below then passed OBS_DATE straight back
  # into estimate_dhis2_onset_sample_delay(), which NEUTRALISES that function's own as-of guard
  # (`keep <- keep & (sa <= .asof)` becomes a no-op), and the resulting parameters are written
  # to data/cfr_reference/dhis2_onset_sample_delay_params.csv — the exact file
  # effective_onset_sample_delay() reads (the reported onset->sample delay, and the nowcast and
  # R(t) truncation fallbacks). So one bad date in a future snapshot could silently widen the
  # fitting window, weaken the truncation correction, and become the pipeline's delay.
  # DHIS2_DELAY_NO_ASOF=1 restores the old unbounded behaviour for a deliberate
  # whole-history characterisation; it is never the path that feeds the pipeline.
  .asof_main <- if (identical(tolower(trimws(Sys.getenv("DHIS2_DELAY_NO_ASOF", ""))), "1")) NULL
                else get0("ANALYSIS_DATE", ifnotfound = NULL)
  if (is.null(.asof_main))
    cat("  [as-of] DISABLED (DHIS2_DELAY_NO_ASOF=1 or no ANALYSIS_DATE): window uses ALL dates.\n")
  else
    cat(sprintf("  [as-of] bound = %s (matches the pipeline's window)\n", format(.asof_main)))
  .pops <- build_dhis2_delay_populations(df_all, analysis_date = .asof_main)
  if (is.null(.pops)) stop("[04c_dhis2] no usable sample dates — cannot derive the delay window.")
  OBS_DATE       <- .pops$obs_date
  max_sample     <- .pops$max_sample
  TRUNC_DATE     <- .pops$trunc_date
  ANALYSIS_START <- .pops$analysis_start
  WINDOW_NAME    <- .pops$window_name
  delay_specs    <- .pops$specs

  cat(sprintf("  ANALYSIS_START (Monday floor)   : %s\n", ANALYSIS_START))
  cat(sprintf("  max(sample_date)                : %s\n", max_sample))
  cat(sprintf("  TRUNC_DATE     (max - %dd)        : %s\n", TEST_DAYS, TRUNC_DATE))
  cat(sprintf("  OBS_DATE (as-of, else max date)  : %s\n\n", OBS_DATE))

  cat("-- INTERVAL-CENSORED MLE (d=0->[0,0.5]; d>0->[d-0.5,d+0.5]) --\n")
  cens_tables <- list(); cens_fits_all <- list(); panels <- list()
  # Retain each delay's windowed date PAIRS (not just the integer delays): the EpiDist
  # block below refits every one of them, and it needs the primary/secondary dates to
  # build the censoring intervals and the right-truncation observation date.
  delay_dat <- list()
  for (nm in names(delay_specs)) {
    sp <- delay_specs[[nm]]
    dat <- .pops$delays[[nm]]
    if (is.null(dat) || length(dat$delay) < 5L) {
      cat(sprintf("  %-24s unavailable / n<5 — skipped\n", sp$lbl)); next
    }
    delay_dat[[nm]]     <- dat
    fits <- .fit_all_censored(dat$delay, sp$lbl)
    cens_tables[[nm]]   <- .build_cens_tbl(fits, nm, WINDOW_NAME, length(dat$delay))
    cens_fits_all[[nm]] <- fits
    if (.HAVE_FIG) panels[[nm]] <- .plot_censored_fits(dat$delay, fits, sp$lbl)
  }
  cens_tbl <- dplyr::bind_rows(cens_tables)

  cat("\n-- Censored MLE results (AIC-ranked per delay) --\n")
  for (dt in unique(cens_tbl$delay_type)) {
    sub <- dplyr::filter(cens_tbl, delay_type == dt)
    cat(sprintf("\n  %s (n=%d):\n", dt, sub$n[1]))
    print(dplyr::select(sub, family, aic, delta_aic, implied_mean_d,
                        shape, rate, scale, meanlog, sdlog), n = Inf)
  }

  # ── EpiDist for EVERY delay (Bayesian; censoring + right-truncation corrected) ─
  # EpiDist is the project's designated estimator for delay parameters, so it is applied
  # to ALL the delays this script characterises — not just onset->sample. It used to be
  # fit for onset->sample alone, leaving onset->notification, onset->lab and sample->lab
  # reported from the interval-censored MLE only, i.e. corrected for daily rounding but
  # NOT for right-truncation. Only onset->sample is written to the delay params file (the
  # reported delay, and the nowcast and R(t) truncation fallbacks); the other three are
  # descriptive, but they are quoted in methods text and QA, so they get the same estimator and
  # the same corrections.
  epidist_tbl <- NULL
  if (RUN_EPIDIST && length(delay_dat)) {
    cat("\n-- EpiDist (Bayesian: naive + marginal x lognormal + gamma), all delays --\n")
    epidist_tbl <- purrr::map_dfr(names(delay_dat), function(nm) {
      d <- delay_dat[[nm]]
      cat(sprintf("\n  [%s] n=%d windowed pairs\n", delay_specs[[nm]]$lbl, length(d$delay)))
      # OBS_DATE = the as-of date (max(sample, onset) when DHIS2_DELAY_NO_ASOF=1) is the
      # truncation reference. It deliberately EXCLUDES date_of_notification / reporting_date,
      # which carry far-future data-entry typos (e.g. 2027) that would inflate the observation
      # date and switch the truncation correction off.
      fit_epidist_both(tibble::tibble(onset = d$onset, sample = d$sample), nm, OBS_DATE)
    })
    if (!is.null(epidist_tbl) && nrow(epidist_tbl)) print(epidist_tbl, n = Inf)
  }

  # ── Assemble the onset->sample fit + write outputs ───────────────────────────
  # Pass the AS-OF date, not OBS_DATE. OBS_DATE is the truncation reference derived from the
  # data; feeding it back as `analysis_date` made the estimator's own `sa <= analysis_date`
  # filter a tautology. When the as-of bound is deliberately disabled, OBS_DATE is the
  # correct (and only) reference available.
  os_fit <- estimate_dhis2_onset_sample_delay(
    df_all, analysis_date = if (is.null(.asof_main)) OBS_DATE else .asof_main)
  if (is.null(os_fit)) stop("[04c_dhis2] onset->sample fit failed (too few complete pairs).")

  out_dir <- file.path(LINELIST_DIR, .meta$folder)
  readr::write_csv(cens_tbl, file.path(out_dir, "dhis2_delay_params_censored.csv"))
  cat(sprintf("\nSaved dhis2_delay_params_censored.csv (%d rows) -> %s\n",
              nrow(cens_tbl), out_dir))
  # Truncation-corrected EpiDist MARGINAL params (best family; review §1.2): the POOLED
  # epidist_* rows read by .load_dhis2_delay_params(). Prefer the gamma marginal
  # (clean mean/sd -> shape/rate); fall back to the first marginal row. NULL when EpiDist
  # was not run or no onset->sample marginal returned a finite mean; write_onset_sample_long()
  # then keeps an EpiDist fit of the same records already on disk, else selects the MLE.
  epi_marg <- NULL
  if (!is.null(epidist_tbl) && nrow(epidist_tbl)) {
    # MUST filter to the onset->sample delay: epidist_tbl now holds all four delays, and
    # without this an onset->notification or sample->lab marginal row could be written into
    # dhis2_onset_sample_delay_params.csv and silently become the pooled onset->sample delay.
    mm <- dplyr::filter(epidist_tbl, .data$delay == "onset_sample",
                        .data$model_type == "marginal", is.finite(.data$mean_post_med))
    if (nrow(mm)) {
      pick <- mm[match("gamma", mm$family), , drop = FALSE]
      if (!nrow(pick) || is.na(pick$family[1])) pick <- mm[1, , drop = FALSE]
      efam <- if (identical(as.character(pick$family[1]), "lognormal")) "lnorm" else as.character(pick$family[1])
      epi_marg <- list(family = efam, mean = pick$mean_post_med[1],
                       sd = pick$sd_post_med[1], n = pick$n[1],
                       mean_lo = pick$mean_post_lo[1], mean_hi = pick$mean_post_hi[1],
                       max_rhat = pick$max_rhat[1], min_ess_bulk = pick$min_ess_bulk[1],
                       min_ess_tail = pick$min_ess_tail[1], n_divergent = pick$n_divergent[1])
      cat(sprintf("    [§1.2] EpiDist marginal (truncation-corrected): family=%s mean=%.2f d sd=%.2f d -> the pooled epidist_* rows\n",
                  epi_marg$family, epi_marg$mean, epi_marg$sd))
    }
  }
  # The per-classification fits, exactly as run_all.R step 1a writes them. Without them this
  # standalone run rewrote the pipeline's params file with every `__<stratum>` row deleted, and
  # the next load fell back to the pooled delay with only a warning.
  # The strata take the SAME truncation reference as the pooled fit above: under
  # DHIS2_DELAY_NO_ASOF=1 that is OBS_DATE, not each stratum's own latest sample.
  os_strata <- tryCatch(estimate_onset_sample_strata(df_all,
                          analysis_date = if (is.null(.asof_main)) OBS_DATE else .asof_main),
                        error = function(e) {
                          warning("[04c_dhis2] delay strata not fitted (", conditionMessage(e), ").",
                                  call. = FALSE); list() })
  long_path <- file.path(out_dir, "dhis2_onset_sample_delay_params.csv")
  write_onset_sample_long(os_fit, long_path, source_label = .meta$folder, epidist = epi_marg,
                          strata = os_strata)
  cat(sprintf("Saved dhis2_onset_sample_delay_params.csv -> %s\n", long_path))
  # Stable "latest" copy so consumers need not resolve the folder version.
  stable_dir <- file.path(DATA_DIR, "cfr_reference")
  write_onset_sample_long(os_fit, file.path(stable_dir, "dhis2_onset_sample_delay_params.csv"),
                          source_label = .meta$folder, epidist = epi_marg, strata = os_strata)
  # SELECTED, as the loader resolves it from the file just written: the pooled selection, and
  # the REPORTED delay (stratum ONSET_IMPUTE_STRATUM, 00_config.R). Printing the MLE here would
  # quote a delay the loader does not resolve whenever EpiDist was fitted.
  .lab_mean <- 1 / get0("DELAY_ONSET_SAMPLE_RATE", ifnotfound = 0.228)
  .sel_rd <- suppressWarnings(.load_dhis2_delay_params(long_path))
  .rep_st <- get0("ONSET_IMPUTE_STRATUM", ifnotfound = NULL)
  .rep_rd <- if (is.null(.rep_st)) .sel_rd else
    suppressWarnings(.load_dhis2_delay_params(long_path, stratum = .rep_st))
  if (!is.null(.sel_rd))
    cat(sprintf("\n*** onset->sample POOLED selection: %s (%s) | mean=%.2f d | n=%s ***\n",
                .sel_rd$estimator, .sel_rd$family, .sel_rd$mean, format(.sel_rd$n_fit)))
  if (!is.null(.rep_rd))
    cat(sprintf("*** REPORTED delay (%s): %s (%s) | mean=%.2f d | lab-linelist reference mean %.2f d ***\n",
                if (is.null(.rep_st)) "pooled"
                else if (is.na(.rep_rd$stratum)) sprintf("stratum '%s' NOT fitted -> pooled", .rep_st)
                else .rep_st,
                .rep_rd$estimator, .rep_rd$family, .rep_rd$mean, .lab_mean))
  if (!is.null(epidist_tbl) && nrow(epidist_tbl))
    readr::write_csv(epidist_tbl, file.path(out_dir, "dhis2_delay_epidist.csv"))

  # ── SELECTED estimator per delay (EpiDist marginal preferred) ────────────────
  # One table saying, for every delay this script characterises, which estimate is the
  # headline one and why — so methods text and QA quote the SAME number the pipeline uses
  # instead of re-deriving it from whichever CSV was opened first. EpiDist marginal
  # (censoring + right-truncation corrected) wins wherever it returned a finite posterior mean
  # (its convergence diagnostics are in dhis2_delay_epidist.csv); the AIC-best
  # interval-censored MLE is the documented fallback. Every delay here is POOLED over case
  # classifications; the reported onset->sample delay is its confirmed stratum. Families use
  # the params file's lower-case R names (gamma, lnorm, weibull, exp) for both estimators.
  sel <- purrr::map_dfr(names(delay_dat), function(nm) {
    cb <- dplyr::filter(cens_tbl, .data$delay_type == nm, .data$best)
    em <- if (!is.null(epidist_tbl) && nrow(epidist_tbl))
      dplyr::filter(epidist_tbl, .data$delay == nm, .data$model_type == "marginal",
                    is.finite(.data$mean_post_med)) else NULL
    if (!is.null(em) && nrow(em)) {
      pick <- em[match("gamma", em$family), , drop = FALSE]
      if (!nrow(pick) || is.na(pick$family[1])) pick <- em[1, , drop = FALSE]
      tibble::tibble(delay = nm, label = delay_specs[[nm]]$lbl, window = WINDOW_NAME,
                     estimator = "epidist_marginal_truncation_corrected",
                     family = if (identical(as.character(pick$family[1]), "lognormal")) "lnorm"
                              else as.character(pick$family[1]),
                     n = pick$n[1],
                     mean_d = pick$mean_post_med[1], sd_d = pick$sd_post_med[1],
                     mean_lo = pick$mean_post_lo[1], mean_hi = pick$mean_post_hi[1])
    } else {
      tibble::tibble(delay = nm, label = delay_specs[[nm]]$lbl, window = WINDOW_NAME,
                     estimator = "interval_censored_mle",
                     family = if (nrow(cb)) tolower(as.character(cb$family[1])) else NA_character_,
                     n = if (nrow(cb)) cb$n[1] else NA_integer_,
                     mean_d = if (nrow(cb)) cb$implied_mean_d[1] else NA_real_,
                     sd_d = NA_real_, mean_lo = NA_real_, mean_hi = NA_real_)
    }
  })
  if (nrow(sel)) {
    cat("\n-- SELECTED delay estimates (EpiDist marginal preferred) --\n")
    print(sel, n = Inf)
    readr::write_csv(sel, file.path(out_dir, "dhis2_delay_selected.csv"))
    readr::write_csv(sel, file.path(stable_dir, "dhis2_delay_selected.csv"))
    cat(sprintf("Saved dhis2_delay_selected.csv -> %s (and %s)\n", out_dir, stable_dir))
  }

  # ── QA figure: top-2 interval-censored MLE fits per delay ────────────────────
  if (.HAVE_FIG && length(panels)) {
    fig <- patchwork::wrap_plots(panels, ncol = 2) +
      patchwork::plot_annotation(
        title = sprintf("DHIS2 reporting delays — interval-censored MLE (window %s to %s)",
                        ANALYSIS_START, TRUNC_DATE),
        subtitle = sprintf("Top-2 AIC families per delay, pooled over case classifications (the reported onset->sample delay is the confirmed stratum) | %s",
                           .meta$folder))
    fig_path <- file.path(OUT_DIAGNOSTICS, sprintf("dhis2_delay_fits_%s.pdf",
                                                   format(TRUNC_DATE, "%Y%m%d")))
    # Retained-figure gate: dhis2_delay_fits is NOT on the published allow-list
    # (04d's dhis2_delay_epidist_fits is the retained delay figure).
    .fk <- get0("figure_is_kept", ifnotfound = NULL)
    if (is.function(.fk) && !.fk(fig_path)) {
      cat(sprintf("QA figure skipped (not on FIGURE_KEEP) -> %s\n", basename(fig_path)))
    } else {
      ggplot2::ggsave(fig_path, fig, width = 11, height = 8, device = figure_pdf_device())
      cat(sprintf("Saved QA figure -> %s\n", fig_path))
    }
  }

  cat(sprintf("\n%s\n04c_dhis2_delay_windows.R complete.\n%s\n",
              strrep("=", 78), strrep("=", 78)))
}
