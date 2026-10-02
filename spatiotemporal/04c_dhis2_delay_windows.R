# ============================================================================
# 04c_dhis2_delay_windows.R — DHIS2 line-list delay analysis (spatiotemporal)
#
# DHIS2 analog of the retired Ituri LAB-linelist delay analysis (the
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

# ONE canonical definition, identical in every module that declares it, so that source()
# ORDER CANNOT CHANGE SEMANTICS. This file's version used to be
#   function(a, b) if (!is.null(a) && !is.na(a)[1]) a else b
# and was declared UNCONDITIONALLY, so sourcing 06 installed it over every other module's and
# it governed the whole cascade run. It differs in two ways that silently corrupt results:
# it falls back whenever the FIRST ELEMENT of a vector is NA (blanking an entire otherwise-good
# vector, e.g. run_cascade.R's R_zone_conjugate column or 33b's pop_vec), and it ERRORS on a
# zero-length left-hand side ("missing value where TRUE/FALSE needed") instead of falling back.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

MAX_DELAY   <- get0("DELAY_MAX_PLAUSIBLE_DAYS", ifnotfound = 60L)  # shared plausibility ceiling (00_config.R): drop longer delays as typos/outliers before fitting
TEST_DAYS   <- 5L     # right-truncation buffer: stop fitting 5 d before the extract
# Bayesian sampler settings (only used when RUN_EPIDIST=TRUE)
STAN_CHAINS <- 2L; STAN_ITER <- 1000L; STAN_WARMUP <- 500L
STAN_CORES  <- min(STAN_CHAINS, 4L)

# ============================================================================
# HELPER FUNCTIONS
# The interval-censored MLE core is REPLICATED VERBATIM from the validated
# the retired lab template, so the DHIS2 fits use the identical,
# audited machinery — same censoring convention, multi-start, and AIC selection.
# ============================================================================

# Robust per-element date parser (mixed ISO / d-m-y / m-d-y), mirroring
# 01_data_prep.R::.parse_date — base as.Date(tryFormats=) picks ONE format from
# the first non-NA element and NA-s every value in a different format.
# DUPLICATE OF 01_data_prep.R's .parse_date — keep the two bodies IN STEP. run_all.R sources
# 04c first and calls estimate_dhis2_onset_sample_delay() BEFORE sourcing 01, so the delay is
# fitted with THIS parser while the line list is later loaded with 01's; 01 is sourced last, so
# its definition then wins for the rest of the session. If the two ever diverge, the delay would
# be estimated on differently-parsed dates from the ones the imputation applies it to, and the
# discrepancy would be invisible. Verified identical in behaviour on mixed ISO / d-m-y /
# datetime / blank / unparseable input.
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

# naive_implied_mean() REMOVED: nothing called it. The censored-MLE path reports
# `implied_mean_d` from the fit object directly (.fit_all_censored), and the EpiDist path
# reports the posterior mean, so a second family->mean converter was one more place for the
# two to disagree.

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
dhis2_delay_specs <- function() list(
  onset_sample        = list(from = "date_of_symptom_onset",     to = "date_of_sample_collection", lbl = "Onset -> sample"),
  onset_notification  = list(from = "date_of_symptom_onset",     to = "date_of_notification",      lbl = "Onset -> notification"),
  onset_lab_analysis  = list(from = "date_of_symptom_onset",     to = "lab_analysis_date",         lbl = "Onset -> lab analysis"),
  sample_lab_analysis = list(from = "date_of_sample_collection", to = "lab_analysis_date",         lbl = "Sample -> lab analysis"))

#' Build the windowed, plausibility-filtered delay populations for every DHIS2 delay.
#'
#' @param ll raw DHIS2 line list (character or already-parsed date columns).
#' @param analysis_date as-of date, or NULL. NULL reproduces the standalone MAIN-block
#'   window exactly (no as-of bound; truncation reference = max observed date). When a
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
  max_sample <- suppressWarnings(max(samp, na.rm = TRUE))
  if (!is.finite(as.numeric(max_sample))) return(NULL)
  trunc_date     <- max_sample - test_days
  analysis_start <- lubridate::floor_date(outbreak_start, unit = "week", week_start = 1L)
  # OBS_DATE (right-truncation reference for EpiDist) deliberately excludes
  # date_of_notification / reporting_date, which carry far-future data-entry typos.
  obs_date <- if (has_asof) .asof else
    suppressWarnings(max(c(df$date_of_sample_collection, df$date_of_symptom_onset), na.rm = TRUE))

  specs <- dhis2_delay_specs()
  delays <- list()
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    if (!all(c(sp$from, sp$to) %in% names(df))) next
    a <- df[[sp$from]]; b <- df[[sp$to]]
    if (has_asof) b[!is.na(b) & b > .asof] <- NA   # as-of: secondary event not yet observable
    d <- as.integer(b - a)
    keep <- !is.na(a) & !is.na(b) & a >= analysis_start & a <= trunc_date &
            is.finite(d) & d >= 0L & d <= max_delay
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) next
    delays[[nm]] <- list(onset = a[keep], sample = b[keep], delay = d[keep])
  }
  list(specs = specs, delays = delays, analysis_start = analysis_start,
       trunc_date = trunc_date, obs_date = obs_date, max_sample = max_sample,
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
    # bin to NA and DROPS the same-day bar entirely (12.3% of onset->sample records on the
    # current snapshot), with only a "Removed rows containing missing values" warning. 04d uses
    # coord_cartesian for exactly this reason.
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
.EPIDIST_DRAWS <- new.env(parent = emptyenv())
.epidist_draws_key <- function(delay, model_type, family)
  paste(delay, model_type, family, sep = "|")

# Deterministic thinning — seq(), never sample(): 01_data_prep.R's onset imputation draws
# from the global RNG stream, so consuming random numbers here would silently shift every
# imputed onset downstream (see the load_linelist() reseeding contract).
.epidist_draws_put <- function(delay, model_type, family, mean_v, sd_v, n,
                               frac_complete, max_keep = 2000L) {
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
              n = n, frac_complete = frac_complete),
         envir = .EPIDIST_DRAWS)
  invisible(NULL)
}

#' Retrieve captured EpiDist posterior draws.
#' @param delay optional delay name to filter on (e.g. "onset_sample").
#' @return unnamed list of list(delay, model_type, family, mean, sd, n, frac_complete).
epidist_draws_get <- function(delay = NULL) {
  out <- mget(ls(.EPIDIST_DRAWS), envir = .EPIDIST_DRAWS)
  if (!is.null(delay)) out <- out[vapply(out, function(x) identical(x$delay, delay), logical(1))]
  unname(out)
}
# epidist_draws_clear() REMOVED: nothing called it. The draw registry is a per-session
# side-channel from the fits to 04d's density/CDF bands; it is rebuilt every run and there is
# no point in a session where clearing it is correct.

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
    # RIGHT-TRUNCATION EXPOSURE. The old metric was
    #   frac_c <- sum(sdate_lwr <= obs_date_val) / nrow(df_input)
    # which is 1.0 BY CONSTRUCTION: df_pairs is already filtered to samples on or before the
    # as-of date and obs_date_val IS that date. So the `< 0.30` guard could never fire and the
    # log always announced "100% complete" for a model whose entire purpose is correcting
    # right-truncation — the single most misleading line this file could print.
    #
    # Measure the real thing instead: the share of pairs whose PRIMARY event is far enough
    # behind the observation date that the full delay support (MAX_DELAY) could have been
    # observed. Everything else is subject to truncation. Measured 2026-09-19 this is
    # 0.41 / 0.41 / 0.41 / 0.34 for the four DHIS2 delays.
    frac_obs <- mean(as.numeric(obs_date_val - df_input$pdate_lwr) >= MAX_DELAY, na.rm = TRUE)
    if (!is.finite(frac_obs)) frac_obs <- 0
    # Threshold DELIBERATELY set far below the observed range (0.34-0.41) so no fit that runs
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
    # onset->sample marginal mean 7.85 d vs 8.00 d on two identical runs), and because
    # .load_dhis2_delay_params() routes this fit into the onset imputation, the nowcast and
    # the R(t) truncation, that MCMC noise propagates into every downstream number and makes
    # a run non-reproducible. 04b_epinowcast.R already seeds its fit the same way.
    fit <- epidist::epidist(model_obj, family = family_fn, chains = STAN_CHAINS,
                            iter = STAN_ITER, warmup = STAN_WARMUP, cores = STAN_CORES,
                            seed = get0("RANDOM_SEED", ifnotfound = 20260704L),
                            refresh = 0, silent = 2,
                            control = list(adapt_delta = 0.95, max_treedepth = 12))
    samps <- epidist::predict_delay_parameters(fit) %>% epidist::add_mean_sd()
    # Side effect only — deposit the draws for 04d's density/CDF bands. Never alters the
    # returned tibble, so every existing consumer of this function is untouched.
    .epidist_draws_put(delay_name, model_type, family_name,
                       samps$mean, samps$sd, n = nrow(df_input), frac_complete = frac_obs)
    tibble(delay = delay_name, model_type = model_type, family = family_name,
           n = nrow(df_input), frac_complete = round(frac_obs, 3),
           mean_post_med = round(stats::median(samps$mean, na.rm = TRUE), 3),
           mean_post_lo  = round(stats::quantile(samps$mean, 0.025, na.rm = TRUE), 3),
           mean_post_hi  = round(stats::quantile(samps$mean, 0.975, na.rm = TRUE), 3),
           sd_post_med   = round(stats::median(samps$sd, na.rm = TRUE), 3))
  }, error = function(e) {
    # A warning(), not only a cat(). This handler swallowed an "object 'frac_c' not found"
    # for every single EpiDist fit — raised AFTER Stan had finished, so four fits were burned
    # per run and the caller quietly fell back to whatever epidist_* rows were already on
    # disk. Nothing in the run summary said so. A cat() into a long log is not a failure
    # signal for the estimator that sets the delay behind the onset imputation, the nowcast
    # weight, EpiNow2's truncation and epinowcast's max_delay.
    cat(sprintf("  [%s %s %s] FAILED: %s\n", delay_name, model_type, family_name, e$message))
    warning(sprintf("[epidist] %s %s %s FAILED: %s", delay_name, model_type, family_name,
                    conditionMessage(e)), call. = FALSE)
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

#' Fit the onset->sample delay SEPARATELY by case classification.
#'
#' WHY. The onset imputation (01_data_prep.R) imputes onsets for CONFIRMED records, so the
#' delay it draws from should be the confirmed-case delay. Until 2026-09-22 it drew from a fit
#' pooled over every classification, and on this line list that pool is majority
#' TEST-NEGATIVE and partly UNADJUDICATED: of the 13,987 windowed onset->sample pairs, 6,823
#' are finally classified `not_a_case` and 2,038 carry no final classification at all, against
#' 4,981 confirmed and 145 suspected/probable. Both groups are swabbed faster than confirmed
#' cases (raw means 6.40 d and 3.56 d against 8.80 d), so the pooled EpiDist marginal
#' returns 7.67 d where the CONFIRMED stratum gives 10.14 d -- a 2.47 d (32%) gap, and imputed
#' onsets landed that much late for the ~23% of confirmed records that carry one. Name the
#' estimator with the number: the interval-censored MLE puts the same contrast at 7.67 vs
#' 9.11 d, and the two pairs have been confused before.
#'
#' Each stratum is fitted INDEPENDENTLY, with its own double-interval censoring and its own
#' right-truncation correction. That is not a nicety: confirmed cases are sampled more slowly,
#' so at any as-of date they are also more right-truncated, and a pooled correction would
#' distort the very contrast this function exists to measure.
#'
#' `suspected` pools `suspected_case` with `probable_case` (n = 4 on this snapshot): fitting
#' four records alone would return an unusable interval.
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
           test_days = test_days, max_delay = max_delay, classes = spec[[nm]]),
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
#' @param classes optional character vector of `final_mve_case_classification` values to
#'   restrict the fit to. NULL (the default) pools every classification, which is what the
#'   pipeline did until 2026-09-22. Pass CONFIRMED_STATUS to fit the delay of the population
#'   the onset imputation actually imputes for -- see estimate_onset_sample_strata() for why
#'   that matters.
estimate_dhis2_onset_sample_epidist <- function(ll, analysis_date = ANALYSIS_DATE,
                                                outbreak_start = OUTBREAK_START,
                                                test_days = TEST_DAYS,
                                                max_delay = MAX_DELAY,
                                                classes = NULL) {
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
    tbl <- fit_epidist_both(tibble::tibble(onset = on[ok], sample = sa[ok]),
                            "onset_sample", obs_date)
    if (is.null(tbl) || !nrow(tbl)) return(NULL)
    mm <- tbl[tbl$model_type == "marginal" & is.finite(tbl$mean_post_med), , drop = FALSE]
    if (!nrow(mm)) return(NULL)
    gi   <- match("gamma", mm$family)
    pick <- if (!is.na(gi)) mm[gi, , drop = FALSE] else mm[1, , drop = FALSE]
    efam <- if (identical(as.character(pick$family[1]), "lognormal")) "lnorm"
            else as.character(pick$family[1])
    # mean_lo/mean_hi are the 95% posterior interval of the MEAN. They are carried so the
    # R(t) right-truncation model can be given an UNCERTAIN delay instead of treating a fitted
    # nuisance parameter as known exactly (EpiNow2's own guidance, and the deployed R window
    # lies entirely inside the truncation-corrected region).
    out <- list(family = efam, mean = pick$mean_post_med[1], sd = pick$sd_post_med[1],
                mean_lo = pick$mean_post_lo[1], mean_hi = pick$mean_post_hi[1],
                n = pick$n[1])
    # Full naive+marginal x lognormal+gamma table carried along as an attribute: run_all.R
    # hands it to 04d so the figure can show ALL four sub-fits (and the naive-vs-marginal
    # truncation correction) rather than only the single row routed into the imputation.
    # An attribute keeps write_onset_sample_long(epidist=)'s expected shape ($family/$mean/
    # $sd/$n) exactly as it was.
    attr(out, "epidist_table") <- tbl
    out
  }, error = function(e) {
    message("[04c_dhis2] EpiDist marginal fit failed (non-fatal): ", conditionMessage(e)); NULL
  })
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
  # dropping the recent tail). No-op when `epidist` is NULL (default run).
  epi_rows <- NULL
  if (!is.null(epidist) && !is.null(epidist$mean) && is.finite(epidist$mean) && epidist$mean > 0) {
    epi_rows <- tibble::tribble(
      ~delay,            ~quantity,         ~value,
      "onset_to_sample", "epidist_family",  as.character(epidist$family),
      "onset_to_sample", "epidist_mean",    as.character(round(epidist$mean, 3)),
      "onset_to_sample", "epidist_sd",      as.character(round(epidist$sd, 3)),
      "onset_to_sample", "epidist_mean_lo", as.character(if (is.null(epidist$mean_lo)) NA else round(epidist$mean_lo, 3)),
      "onset_to_sample", "epidist_mean_hi", as.character(if (is.null(epidist$mean_hi)) NA else round(epidist$mean_hi, 3)),
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
  # SELECTED-* rows: state unambiguously which numbers the pipeline actually uses.
  # Without them this file is internally contradictory. The `family`/`estimator`/`rate`/
  # `param_*`/`implied_mean_fit` rows above describe the interval-censored MLE, but
  # .load_dhis2_delay_params() (00_config.R) PREFERS the epidist_* rows and re-derives the
  # gamma parameters by method of moments — so the file simultaneously asserted two different
  # gamma parameterisations (6.832 d, rate 0.1464, shape 0.833 vs 7.667 d, rate 0.1059,
  # shape 0.812) and an `estimator` that is false whenever EpiDist ran. A human, a methods
  # document or a future script reading the obvious keys got the numbers the model does NOT use.
  # These rows mirror the loader's own preference order, so `selected_*` is always the truth.
  .use_epi <- !is.null(epi_rows) && nrow(epi_rows) > 0
  sel_mean <- if (.use_epi) epidist$mean else fit$implied_mean
  sel_sd   <- if (.use_epi) epidist$sd   else NA_real_
  sel_fam  <- if (.use_epi) as.character(epidist$family) else fam
  sel_par  <- if (.use_epi && identical(sel_fam, "gamma") &&
                  is.finite(sel_sd) && sel_sd > 0)
                c(shape = (sel_mean / sel_sd)^2, rate = sel_mean / sel_sd^2)
              else p
  sel_rows <- tibble::tibble(
    delay = "onset_to_sample",
    quantity = c("selected_estimator", "selected_family", "selected_mean", "selected_sd",
                 "selected_rate",
                 paste0("selected_param_", names(sel_par))),
    value = c(if (.use_epi) "epidist_marginal_truncation_corrected" else "interval_censored_mle",
              sel_fam,
              as.character(round(sel_mean, 3)),
              as.character(if (is.finite(sel_sd)) round(sel_sd, 3) else NA),
              as.character(round(1 / sel_mean, 5)),
              as.character(round(unname(sel_par), 5))))
  # NON-DESTRUCTIVE when EpiDist did not run. Writing this file with `epidist = NULL` used to
  # DROP any epidist_* rows already on disk, silently downgrading the pipeline's delay from the
  # truncation-corrected marginal fit to the interval-censored MLE (7.67 d -> 6.83 d, a 12%
  # shift that propagates into the onset imputation, the nowcast completeness and the R(t)
  # truncation model). That is a destructive side effect of merely re-running this script
  # without RUN_EPIDIST, and it happened. If the caller supplies no EpiDist fit but the target
  # already carries one, KEEP the existing rows and say so, rather than quietly regressing.
  if (is.null(epi_rows) && file.exists(path)) {
    .prev <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(.prev) && all(c("delay", "quantity", "value") %in% names(.prev))) {
      .keep <- .prev[grepl("^epidist_", .prev$quantity), , drop = FALSE]
      if (nrow(.keep)) {
        warning(sprintf(paste0("[04c_dhis2] no EpiDist fit was supplied, but %s already carries ",
                               "one. PRESERVING the existing epidist_* rows rather than ",
                               "downgrading the pipeline to the interval-censored MLE. Re-run ",
                               "with RUN_EPIDIST=TRUE to refresh them."), basename(path)),
                call. = FALSE, immediate. = TRUE)
        epi_rows <- .keep
        .use_epi <- TRUE
        sel_mean <- suppressWarnings(as.numeric(.keep$value[.keep$quantity == "epidist_mean"]))[1]
        sel_sd   <- suppressWarnings(as.numeric(.keep$value[.keep$quantity == "epidist_sd"]))[1]
        sel_fam  <- as.character(.keep$value[.keep$quantity == "epidist_family"])[1]
        if (identical(sel_fam, "gamma") && is.finite(sel_mean) && is.finite(sel_sd) && sel_sd > 0)
          sel_par <- c(shape = (sel_mean / sel_sd)^2, rate = sel_mean / sel_sd^2)
        sel_rows <- tibble::tibble(
          delay = "onset_to_sample",
          quantity = c("selected_estimator", "selected_family", "selected_mean", "selected_sd",
                       "selected_rate", paste0("selected_param_", names(sel_par))),
          value = c("epidist_marginal_truncation_corrected", sel_fam,
                    as.character(round(sel_mean, 3)), as.character(round(sel_sd, 3)),
                    as.character(round(1 / sel_mean, 5)),
                    as.character(round(unname(sel_par), 5))))
      }
    }
  }
  # PER-CLASSIFICATION rows. Additive only: nothing above is altered, so a reader that does
  # not know about strata resolves exactly what it resolved before. Note these names match
  # the `^epidist_` pattern used by the preserve-on-rerun branch above, which is deliberate --
  # re-running without RUN_EPIDIST keeps the strata too rather than dropping them.
  strata_rows <- NULL
  if (!is.null(strata) && length(strata)) {
    strata_rows <- dplyr::bind_rows(lapply(names(strata), function(nm) {
      e <- strata[[nm]]
      if (is.null(e) || !is.finite(e$mean) || e$mean <= 0) return(NULL)
      .n <- function(x) as.character(if (is.null(x) || !is.finite(x)) NA else round(x, 3))
      tibble::tibble(
        delay = "onset_to_sample",
        quantity = paste0(c("epidist_family", "epidist_mean", "epidist_sd",
                            "epidist_mean_lo", "epidist_mean_hi", "epidist_n"), "__", nm),
        value = c(as.character(e$family), .n(e$mean), .n(e$sd), .n(e$mean_lo), .n(e$mean_hi),
                  as.character(if (is.null(e$n)) NA else e$n)))
    }))
  }
  out <- dplyr::bind_rows(rows, par_rows, epi_rows, sel_rows, strata_rows)
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

  # AS-OF BOUND, by default, even standalone. This used to pass analysis_date = NULL, which
  # deliberately disables the bound: max_sample, TRUNC_DATE and OBS_DATE were then taken over
  # ALL rows, including any future-dated onset or sample typo (this file's own comment at the
  # builder concedes such typos exist). Line 784 below then passes OBS_DATE straight back into
  # estimate_dhis2_onset_sample_delay(), which NEUTRALISES that function's own as-of guard
  # (`keep <- keep & (sa <= .asof)` becomes a no-op), and the resulting parameters are written
  # to data/cfr_reference/dhis2_onset_sample_delay_params.csv — the exact file
  # effective_onset_sample_delay() reads for the onset imputation, the nowcast and the R(t)
  # truncation model. So one bad date in a future snapshot could silently widen the fitting
  # window, weaken the truncation correction, and become the pipeline's delay.
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

  cat(sprintf("  ANALYSIS_START (outbreak floor) : %s\n", ANALYSIS_START))
  cat(sprintf("  max(sample_date)                : %s\n", max_sample))
  cat(sprintf("  TRUNC_DATE     (max - %dd)        : %s\n", TEST_DAYS, TRUNC_DATE))
  cat(sprintf("  OBS_DATE       (max all dates)   : %s\n\n", OBS_DATE))

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

  # ── Naive MLE (for the naive-vs-censored comparison) ─────────────────────────
  os_dat   <- .pops$delays[["onset_sample"]]
  naive_os <- if (!is.null(os_dat)) fit_mle_naive(os_dat$delay, "onset->sample") else NULL

  # ── EpiDist for EVERY delay (Bayesian; censoring + right-truncation corrected) ─
  # EpiDist is the project's designated estimator for delay parameters, so it is applied
  # to ALL the delays this script characterises — not just onset->sample. It used to be
  # fit for onset->sample alone, leaving onset->notification, onset->lab and sample->lab
  # reported from the interval-censored MLE only, i.e. corrected for daily rounding but
  # NOT for right-truncation. Only onset->sample feeds the models (imputation / nowcast /
  # R(t) truncation); the other three are descriptive, but they are quoted in methods text
  # and QA, so they get the same estimator and the same corrections.
  epidist_tbl <- NULL
  if (RUN_EPIDIST && length(delay_dat)) {
    cat("\n-- EpiDist (Bayesian: naive + marginal x lognormal + gamma), all delays --\n")
    epidist_tbl <- purrr::map_dfr(names(delay_dat), function(nm) {
      d <- delay_dat[[nm]]
      cat(sprintf("\n  [%s] n=%d windowed pairs\n", delay_specs[[nm]]$lbl, length(d$delay)))
      # OBS_DATE = max(sample, onset) is the extraction-date proxy. It deliberately EXCLUDES
      # date_of_notification / reporting_date, which carry far-future data-entry typos (e.g.
      # 2027) that would inflate the observation date and switch the truncation correction off.
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
    # MUST filter to the onset->sample delay: epidist_tbl now holds all four delays, and
    # without this an onset->notification or sample->lab marginal row could be written into
    # dhis2_onset_sample_delay_params.csv and silently become the imputation/nowcast delay.
    mm <- dplyr::filter(epidist_tbl, .data$delay == "onset_sample",
                        .data$model_type == "marginal", is.finite(.data$mean_post_med))
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

  # ── SELECTED estimator per delay (EpiDist marginal preferred) ────────────────
  # One table saying, for every delay this script characterises, which estimate is the
  # headline one and why — so methods text and QA quote the SAME number the pipeline uses
  # instead of re-deriving it from whichever CSV was opened first. EpiDist marginal
  # (censoring + right-truncation corrected) wins wherever it converged; the AIC-best
  # interval-censored MLE is the documented fallback.
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
                     family = as.character(pick$family[1]), n = pick$n[1],
                     mean_d = pick$mean_post_med[1], sd_d = pick$sd_post_med[1],
                     mean_lo = pick$mean_post_lo[1], mean_hi = pick$mean_post_hi[1])
    } else {
      tibble::tibble(delay = nm, label = delay_specs[[nm]]$lbl, window = WINDOW_NAME,
                     estimator = "interval_censored_mle",
                     family = if (nrow(cb)) as.character(cb$family[1]) else NA_character_,
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
    # Retained-figure gate: dhis2_delay_fits is NOT on the published allow-list
    # (04d's dhis2_delay_epidist_fits is the retained delay figure).
    .fk <- get0("figure_is_kept", ifnotfound = NULL)
    if (is.function(.fk) && !.fk(fig_path)) {
      cat(sprintf("QA figure skipped (not on FIGURE_KEEP) -> %s\n", basename(fig_path)))
    } else {
      ggplot2::ggsave(fig_path, fig, width = 11, height = 8, device = "pdf")
      cat(sprintf("Saved QA figure -> %s\n", fig_path))
    }
  }

  cat(sprintf("\n%s\n04c_dhis2_delay_windows.R complete.\n%s\n",
              strrep("=", 78), strrep("=", 78)))
}
