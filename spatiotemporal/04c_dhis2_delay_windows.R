# ============================================================================
# 04c_dhis2_delay_windows.R — DHIS2 line-list delay analysis (spatiotemporal)
#
# DHIS2 analog of 04c_delay_four_windows.R (the Ituri LAB-linelist delay
# analysis). It estimates the DHIS2 reporting delays — primarily the
# onset->sample delay that drives the missing-onset imputation in 01_data_prep.R
# — with the SAME rigorous machinery, adapted to the DHIS2 line list and wired to
# the spatiotemporal 00_config.R (no dependency on the CFR pipeline's config_r.R).
#
# WHY a DHIS2-specific fit: the fixed lab-linelist delay (Exp rate 0.228, mean
# 4.39 d) understates the DHIS2 onset->sample delay, which is materially longer
# (empirically ~5.8 d). Imputing DHIS2 onsets with the lab delay biases them late.
#
# METHODS (mirroring the lab template):
#   * Windowing: onset in [ANALYSIS_START, TRUNC_DATE = max(sample_date) - 5 d].
#     The final 5 days before the extract are heavily right-truncated (recent
#     onsets have not yet been sampled), so delay fitting stops 5 days before the
#     line-list end — matching the training cutoff convention in 15/16. Both dates
#     are DERIVED from the data (no hard-coded cutoffs). ANALYSIS_START defaults to
#     OUTBREAK_START.
#   * A. Naive uncensored MLE   — fitdist(); zeros excluded; lnorm uses x+0.5.
#   * B. Interval-censored MLE  — fitdistcens(); d=0 -> [0,0.5]; d>0 -> [d-0.5,d+0.5]
#        (corrects daily-rounding censoring); AIC-selected across gamma/lnorm/
#        weibull/exp. This is the PRIMARY, always-on estimator.
#   * C. Bayesian EpiDist       — optional (RUN_EPIDIST=TRUE); the MARGINAL model
#        corrects BOTH interval censoring AND right-truncation.
#
# OUTPUTS (written next to the source line list, folder-versioned):
#   * dhis2_delay_params_censored.csv     — full per-family censored MLE table
#       (window, delay_type, family, n, aic, bic, implied_mean_d, params, best).
#   * dhis2_onset_sample_delay_params.csv — long (delay,quantity,value) params for
#       the SELECTED onset->sample family, mirroring the CFR onset_to_sample CSV,
#       consumed by 01_data_prep.R for the onset imputation.
#   * <OUT_DIAGNOSTICS>/dhis2_delay_fits_<window>.pdf — naive-vs-censored panels.
#
# REUSABLE: estimate_dhis2_onset_sample_delay(ll, analysis_date) returns the fit
# list (best family, params, rate, windowed delay vector, censored table) so
# 01_data_prep.R (or an interactive session) can obtain the DHIS2 delay without
# re-running the whole script.
#
# Run: Rscript spatiotemporal/04c_dhis2_delay_windows.R
#   env RUN_EPIDIST=TRUE  -> also fit the Bayesian truncation-corrected model.
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
# EpiDist (Bayesian, truncation-correcting) is the DEFAULT delay estimator (00_config.R
# sets RUN_EPIDIST, default TRUE), gated here on package availability. A standalone run
# that has NOT sourced 00_config.R falls back to env RUN_EPIDIST (default FALSE), so this
# script keeps working in isolation; an explicit env var overrides the config value.
.HAVE_EPIDIST <- requireNamespace("epidist", quietly = TRUE) &&
                 requireNamespace("epiparameter", quietly = TRUE) &&
                 requireNamespace("brms", quietly = TRUE)
.RUN_EPIDIST_DEFAULT <- isTRUE(get0("RUN_EPIDIST", ifnotfound = FALSE))
RUN_EPIDIST <- (tolower(trimws(Sys.getenv("RUN_EPIDIST",
                 if (.RUN_EPIDIST_DEFAULT) "true" else "false"))) %in%
                c("true", "t", "1", "yes", "y")) && .HAVE_EPIDIST
if (RUN_EPIDIST) suppressPackageStartupMessages({
  library(epidist); library(epiparameter); library(brms)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

MAX_DELAY   <- get0("DELAY_MAX_PLAUSIBLE_DAYS", ifnotfound = 60L)  # shared plausibility ceiling (00_config.R): drop longer delays as typos/outliers before fitting
TEST_DAYS   <- 5L     # right-truncation buffer: stop fitting 5 d before the extract
# Bayesian sampler settings (only used when RUN_EPIDIST=TRUE)
STAN_CHAINS <- 2L; STAN_ITER <- 1000L; STAN_WARMUP <- 500L
STAN_CORES  <- min(STAN_CHAINS, 4L)

# ============================================================================
# HELPER FUNCTIONS
# The interval-censored MLE core is REPLICATED VERBATIM from the validated
# 04c_delay_four_windows.R (lab template) so the DHIS2 fits use the identical,
# audited machinery — same censoring convention, multi-start, and AIC selection.
# ============================================================================

# Robust per-element date parser (mixed ISO / d-m-y / m-d-y), mirroring
# 01_data_prep.R::.parse_date — base as.Date(tryFormats=) picks ONE format from
# the first non-NA element and NA-s every value in a different format.
.parse_date <- function(x) {
  if (inherits(x, "Date"))   return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  raw <- trimws(as.character(x)); raw[raw == ""] <- NA_character_
  fmts <- c("%Y-%m-%d", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S", "%d/%m/%Y", "%m/%d/%Y")
  out <- as.Date(rep(NA_real_, length(raw)), origin = "1970-01-01")
  for (fmt in fmts) {
    todo <- is.na(out) & !is.na(raw)
    if (!any(todo)) break
    out[todo] <- suppressWarnings(as.Date(raw[todo], format = fmt))
  }
  out
}

# ── A. Naive uncensored MLE ───────────────────────────────────────────────────
fit_mle_naive <- function(x, delay_name) {
  x     <- x[!is.na(x) & x >= 0]
  x_pos <- x[x > 0]
  if (length(x_pos) < 5) {
    cat(sprintf("  [%s] n_pos=%d < 5, skipping naive MLE\n", delay_name, length(x_pos)))
    return(NULL)
  }
  out <- list()
  for (fam in c("gamma", "lnorm", "weibull", "exp")) {
    data_fit <- if (fam == "lnorm") x_pos + 0.5 else x_pos
    tryCatch({
      fit <- fitdist(data_fit, fam, method = "mle")
      out[[fam]] <- list(family = fam, params = fit$estimate,
                         aic = fit$aic, bic = fit$bic, fit = fit)
    }, error = function(e)
      cat(sprintf("  [%s] %s naive fit failed: %s\n", delay_name, fam, e$message)))
  }
  out
}

naive_implied_mean <- function(fam, params) {
  tryCatch(switch(fam,
    gamma   = unname(params["shape"] / params["rate"]),
    lnorm   = unname(exp(params["meanlog"] + 0.5 * params["sdlog"]^2) - 0.5),
    weibull = unname(params["scale"] * gamma(1 + 1 / params["shape"])),
    exp     = unname(1 / params["rate"]),
    NA_real_), error = function(e) NA_real_)
}

# ── B. Interval-censored MLE ──────────────────────────────────────────────────
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

.plot_censored_fits <- function(x, fits, panel_title, xlim_max = NULL) {
  x_c <- x[!is.na(x) & x >= 0L]
  if (length(x_c) == 0) return(ggplot() + labs(title = panel_title) + theme_void())
  if (is.null(xlim_max)) xlim_max <- min(max(x_c) + 2L, MAX_DELAY)
  xg  <- seq(0.01, xlim_max, by = 0.1)
  p   <- ggplot(data.frame(delay = x_c), aes(x = delay)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                   fill = "grey85", colour = "grey60", alpha = 0.75, linewidth = 0.3) +
    xlim(0, xlim_max) +
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
# Mirrors the lab template's epidist path. The MARGINAL model right-truncates at
# obs_date; the NAIVE model does not. Returns posterior mean/sd (median + 95% CrI).
fit_epidist_model <- function(df_pairs, obs_date_val, model_type, family_fn,
                              family_name, delay_name) {
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
    frac_c <- sum(df_input$sdate_lwr <= obs_date_val, na.rm = TRUE) / nrow(df_input)
    if (frac_c < 0.30) {
      cat(sprintf("  [%s %s %s] skipped — %.0f%% right-truncated\n",
                  delay_name, model_type, family_name, (1 - frac_c) * 100))
      return(NULL)
    }
    cat(sprintf("  Fitting EpiDist %s %s %s (n=%d, %.0f%% complete)...\n",
                delay_name, model_type, family_name, nrow(df_input), frac_c * 100))
    ll_obj <- df_input %>%
      epidist::as_epidist_linelist_data(
        pdate_lwr = "pdate_lwr", pdate_upr = "pdate_upr",
        sdate_lwr = "sdate_lwr", sdate_upr = "sdate_upr", obs_date = "obs_date")
    model_obj <- switch(model_type,
      naive    = epidist::as_epidist_naive_model(ll_obj),
      marginal = ll_obj %>% epidist::as_epidist_aggregate_data() %>%
                   epidist::as_epidist_marginal_model())
    fit <- epidist::epidist(model_obj, family = family_fn, chains = STAN_CHAINS,
                            iter = STAN_ITER, warmup = STAN_WARMUP, cores = STAN_CORES,
                            refresh = 0, silent = 2,
                            control = list(adapt_delta = 0.95, max_treedepth = 12))
    samps <- epidist::predict_delay_parameters(fit) %>% epidist::add_mean_sd()
    tibble(delay = delay_name, model_type = model_type, family = family_name,
           n = nrow(df_input), frac_complete = round(frac_c, 3),
           mean_post_med = round(stats::median(samps$mean, na.rm = TRUE), 3),
           mean_post_lo  = round(stats::quantile(samps$mean, 0.025, na.rm = TRUE), 3),
           mean_post_hi  = round(stats::quantile(samps$mean, 0.975, na.rm = TRUE), 3),
           sd_post_med   = round(stats::median(samps$sd, na.rm = TRUE), 3))
  }, error = function(e) {
    cat(sprintf("  [%s %s %s] FAILED: %s\n", delay_name, model_type, family_name, e$message))
    NULL
  })
}

fit_epidist_both <- function(df_pairs, delay_name, obs_date_val) {
  fams <- list(lognormal = quote(brms::lognormal()), gamma = quote(stats::Gamma(link = "log")))
  purrr::map_dfr(names(fams), function(fam) {
    fam_fn <- eval(fams[[fam]])
    purrr::map_dfr(c("naive", "marginal"), function(mt)
      fit_epidist_model(df_pairs, obs_date_val, mt, fam_fn, fam, delay_name))
  })
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
#'   the reporting delay is a lab/reporting property, ~invariant to classification,
#'   and all complete pairs give the largest, most stable sample).
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
    # regardless of the AIC-best family — so the imputation's Exp fallback and the
    # nowcast CDF stay well-defined. The full best-family params are also returned.
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
#' SAME windowed complete pairs and returns list(family, mean, sd, n) — exactly the shape
#' write_onset_sample_long(epidist=) writes and .load_dhis2_delay_params() PREFERS — so
#' run_all.R routes the truncation-corrected delay into onset imputation + the nowcast rate
#' rather than only into the standalone Rscript path. Mirrors the MAIN block's §1.2
#' extraction (prefer the gamma marginal, else the first marginal row; map lognormal->lnorm).
#' Returns NULL when RUN_EPIDIST is off, the `epidist` package is unavailable (RUN_EPIDIST
#' already folds in .HAVE_EPIDIST), there are too few pairs, or the fit fails — so the caller
#' falls back to the censored-MLE cleanly and the default remains non-fatal.
estimate_dhis2_onset_sample_epidist <- function(ll, analysis_date = ANALYSIS_DATE,
                                                outbreak_start = OUTBREAK_START,
                                                test_days = TEST_DAYS,
                                                max_delay = MAX_DELAY) {
  if (!isTRUE(get0("RUN_EPIDIST", ifnotfound = FALSE)) ||
      !exists("fit_epidist_both", mode = "function")) return(NULL)
  on <- .parse_date(ll[["date_of_symptom_onset"]])
  sa <- .parse_date(ll[["date_of_sample_collection"]])
  keep <- !is.na(on) & !is.na(sa)
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
    tbl <- fit_epidist_both(tibble::tibble(onset = on[ok], sample = sa[ok]),
                            "onset_sample", obs_date)
    if (is.null(tbl) || !nrow(tbl)) return(NULL)
    mm <- tbl[tbl$model_type == "marginal" & is.finite(tbl$mean_post_med), , drop = FALSE]
    if (!nrow(mm)) return(NULL)
    gi   <- match("gamma", mm$family)
    pick <- if (!is.na(gi)) mm[gi, , drop = FALSE] else mm[1, , drop = FALSE]
    efam <- if (identical(as.character(pick$family[1]), "lognormal")) "lnorm"
            else as.character(pick$family[1])
    list(family = efam, mean = pick$mean_post_med[1], sd = pick$sd_post_med[1], n = pick$n[1])
  }, error = function(e) {
    message("[04c_dhis2] EpiDist marginal fit failed (non-fatal): ", conditionMessage(e)); NULL
  })
}

# Write the SELECTED onset->sample family to the long (delay,quantity,value) CSV
# the pipeline consumes (mirrors data/cfr_reference/onset_to_sample_delay_params.csv).
write_onset_sample_long <- function(fit, path, source_label, epidist = NULL) {
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
  # dropping the recent tail). No-op when `epidist` is NULL (default run).
  epi_rows <- NULL
  if (!is.null(epidist) && !is.null(epidist$mean) && is.finite(epidist$mean) && epidist$mean > 0) {
    epi_rows <- tibble::tribble(
      ~delay,            ~quantity,         ~value,
      "onset_to_sample", "epidist_family",  as.character(epidist$family),
      "onset_to_sample", "epidist_mean",    as.character(round(epidist$mean, 3)),
      "onset_to_sample", "epidist_sd",      as.character(round(epidist$sd, 3)),
      "onset_to_sample", "epidist_n",       as.character(if (is.null(epidist$n)) NA else epidist$n),
      "onset_to_sample", "epidist_estimator", "epidist_marginal_truncation_corrected")
  }
  # Append the best family's native parameters so a consumer can reconstruct the
  # full fitted distribution (not only the Exp-rate summary). PREFIX with `param_`
  # so a family whose native parameter is itself called `rate` (gamma, exp) does
  # NOT collide with the Exp-rate summary row above (which would make deframe()
  # ambiguous). e.g. gamma -> param_shape, param_rate.
  par_rows <- purrr::imap_dfr(as.list(round(p, 5)), function(v, nm)
    tibble(delay = "onset_to_sample", quantity = paste0("param_", nm),
           value = as.character(unname(v))))
  out <- dplyr::bind_rows(rows, par_rows, epi_rows)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(out, path)
  out
}

# ============================================================================
# MAIN — run only when executed as a script (sourcing this file for its functions
# from 01_data_prep.R must NOT trigger the fit). `sys.nframe() == 0L` is TRUE only
# for `Rscript this.R` / `R CMD BATCH`, FALSE when source()d.
# ============================================================================
if (sys.nframe() == 0L || isTRUE(getOption("dhis2_delay_run_main"))) {
  cat(strrep("=", 78), "\n")
  cat("04c_dhis2_delay_windows.R — DHIS2 line-list delay analysis\n")
  cat(sprintf("Naive + interval-censored MLE%s | %s\n",
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
  cat(sprintf("  Loaded: %d records\n", nrow(df_all)))

  date_cols <- c("date_of_symptom_onset", "date_of_sample_collection",
                 "date_of_notification", "reporting_date", "lab_analysis_date")
  for (dc in intersect(date_cols, names(df_all))) df_all[[dc]] <- .parse_date(df_all[[dc]])

  OBS_DATE   <- suppressWarnings(max(c(df_all$date_of_sample_collection,
                                       df_all$date_of_symptom_onset), na.rm = TRUE))
  max_sample <- suppressWarnings(max(df_all$date_of_sample_collection, na.rm = TRUE))
  TRUNC_DATE <- max_sample - TEST_DAYS
  ANALYSIS_START <- lubridate::floor_date(OUTBREAK_START, unit = "week", week_start = 1L)
  WINDOW_NAME <- sprintf("analytical_%s", format(TRUNC_DATE, "%Y-%m-%d"))

  cat(sprintf("  ANALYSIS_START (outbreak floor) : %s\n", ANALYSIS_START))
  cat(sprintf("  max(sample_date)                : %s\n", max_sample))
  cat(sprintf("  TRUNC_DATE     (max - %dd)        : %s\n", TEST_DAYS, TRUNC_DATE))
  cat(sprintf("  OBS_DATE       (max all dates)   : %s\n\n", OBS_DATE))

  # ── Build windowed delay populations (onset in [start, TRUNC]) ───────────────
  in_win <- function(onset) !is.na(onset) & onset >= ANALYSIS_START & onset <= TRUNC_DATE
  mk_delay <- function(from, to) {
    if (!all(c(from, to) %in% names(df_all))) return(NULL)
    d <- as.integer(df_all[[to]] - df_all[[from]])
    keep <- in_win(df_all[[from]]) & !is.na(df_all[[from]]) & !is.na(df_all[[to]]) &
            is.finite(d) & d >= 0L & d <= MAX_DELAY
    list(onset = df_all[[from]][keep], sample = df_all[[to]][keep], delay = d[keep])
  }
  delay_specs <- list(
    onset_sample          = list(from = "date_of_symptom_onset",   to = "date_of_sample_collection", lbl = "Onset -> sample"),
    onset_notification    = list(from = "date_of_symptom_onset",   to = "date_of_notification",       lbl = "Onset -> notification"),
    onset_lab_analysis    = list(from = "date_of_symptom_onset",   to = "lab_analysis_date",          lbl = "Onset -> lab analysis"),
    sample_lab_analysis   = list(from = "date_of_sample_collection", to = "lab_analysis_date",        lbl = "Sample -> lab analysis"))

  cat("-- INTERVAL-CENSORED MLE (d=0->[0,0.5]; d>0->[d-0.5,d+0.5]) --\n")
  cens_tables <- list(); cens_fits_all <- list(); panels <- list()
  for (nm in names(delay_specs)) {
    sp <- delay_specs[[nm]]
    dat <- mk_delay(sp$from, sp$to)
    if (is.null(dat) || length(dat$delay) < 5L) {
      cat(sprintf("  %-24s unavailable / n<5 — skipped\n", sp$lbl)); next
    }
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

  # ── Naive MLE (for the naive-vs-censored comparison) ─────────────────────────
  os_dat <- mk_delay("date_of_symptom_onset", "date_of_sample_collection")
  naive_os <- fit_mle_naive(os_dat$delay, "onset->sample")

  # ── EpiDist (optional; truncation-corrected onset->sample) ───────────────────
  epidist_tbl <- NULL
  if (RUN_EPIDIST && !is.null(os_dat)) {
    cat("\n-- EpiDist (Bayesian: naive + marginal x lognormal + gamma) --\n")
    epidist_tbl <- fit_epidist_both(
      tibble::tibble(onset = os_dat$onset, sample = os_dat$sample),
      "onset_sample", OBS_DATE)
    if (!is.null(epidist_tbl) && nrow(epidist_tbl)) print(epidist_tbl, n = Inf)
  }

  # ── Assemble the onset->sample fit + write outputs ───────────────────────────
  os_fit <- estimate_dhis2_onset_sample_delay(df_all, analysis_date = OBS_DATE)
  if (is.null(os_fit)) stop("[04c_dhis2] onset->sample fit failed (too few complete pairs).")

  best_os <- dplyr::filter(cens_tbl, delay_type == "onset_sample", best)
  cat(sprintf("\n*** onset->sample SELECTED: %s | mean=%.2f d | Exp-rate summary=%.4f /d | n=%d ***\n",
              os_fit$best_family, os_fit$implied_mean, os_fit$rate, os_fit$n_fit))
  cat(sprintf("    (lab-linelist reference was Exp rate %.3f /d, mean %.2f d — DHIS2 delay is %s)\n",
              get0("DELAY_ONSET_SAMPLE_RATE", ifnotfound = 0.228),
              1 / get0("DELAY_ONSET_SAMPLE_RATE", ifnotfound = 0.228),
              if (os_fit$implied_mean > 1 / get0("DELAY_ONSET_SAMPLE_RATE", ifnotfound = 0.228))
                "LONGER" else "shorter"))

  out_dir <- file.path(LINELIST_DIR, .meta$folder)
  readr::write_csv(cens_tbl, file.path(out_dir, "dhis2_delay_params_censored.csv"))
  cat(sprintf("\nSaved dhis2_delay_params_censored.csv (%d rows) -> %s\n",
              nrow(cens_tbl), out_dir))
  # Truncation-corrected EpiDist MARGINAL params (best family; review §1.2) to route into
  # onset imputation / nowcast via .load_dhis2_delay_params(). Prefer the gamma marginal
  # (clean mean/sd -> shape/rate); fall back to the first marginal row. NULL when EpiDist
  # was not run (RUN_EPIDIST=FALSE), so the default run is byte-identical to before.
  epi_marg <- NULL
  if (!is.null(epidist_tbl) && nrow(epidist_tbl)) {
    mm <- dplyr::filter(epidist_tbl, .data$model_type == "marginal", is.finite(.data$mean_post_med))
    if (nrow(mm)) {
      pick <- mm[match("gamma", mm$family), , drop = FALSE]
      if (!nrow(pick) || is.na(pick$family[1])) pick <- mm[1, , drop = FALSE]
      efam <- if (identical(as.character(pick$family[1]), "lognormal")) "lnorm" else as.character(pick$family[1])
      epi_marg <- list(family = efam, mean = pick$mean_post_med[1],
                       sd = pick$sd_post_med[1], n = pick$n[1])
      cat(sprintf("    [§1.2] EpiDist marginal (truncation-corrected): family=%s mean=%.2f d sd=%.2f d -> routed into imputation/nowcast\n",
                  epi_marg$family, epi_marg$mean, epi_marg$sd))
    }
  }
  long_path <- file.path(out_dir, "dhis2_onset_sample_delay_params.csv")
  write_onset_sample_long(os_fit, long_path, source_label = .meta$folder, epidist = epi_marg)
  cat(sprintf("Saved dhis2_onset_sample_delay_params.csv -> %s\n", long_path))
  # Stable "latest" copy so consumers need not resolve the folder version.
  stable_dir <- file.path(DATA_DIR, "cfr_reference")
  write_onset_sample_long(os_fit, file.path(stable_dir, "dhis2_onset_sample_delay_params.csv"),
                          source_label = .meta$folder, epidist = epi_marg)
  if (!is.null(epidist_tbl) && nrow(epidist_tbl))
    readr::write_csv(epidist_tbl, file.path(out_dir, "dhis2_delay_epidist.csv"))

  # ── QA figure: naive vs censored per delay ───────────────────────────────────
  if (.HAVE_FIG && length(panels)) {
    fig <- patchwork::wrap_plots(panels, ncol = 2) +
      patchwork::plot_annotation(
        title = sprintf("DHIS2 reporting delays — interval-censored MLE (window %s to %s)",
                        ANALYSIS_START, TRUNC_DATE),
        subtitle = sprintf("Top-2 AIC families per delay | onset->sample drives onset imputation (01_data_prep.R) | %s",
                           .meta$folder))
    fig_path <- file.path(OUT_DIAGNOSTICS, sprintf("dhis2_delay_fits_%s.pdf",
                                                   format(TRUNC_DATE, "%Y%m%d")))
    ggplot2::ggsave(fig_path, fig, width = 11, height = 8, device = "pdf")
    cat(sprintf("Saved QA figure -> %s\n", fig_path))
  }

  cat(sprintf("\n%s\n04c_dhis2_delay_windows.R complete.\n%s\n",
              strrep("=", 78), strrep("=", 78)))
}
