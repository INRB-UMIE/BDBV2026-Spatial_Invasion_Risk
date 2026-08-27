# =============================================================================
# 02_epi_params.R — Epidemiological Parameter Estimation
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: Estimate and organise all epidemiological parameters required by
#   the spatiotemporal model suite.  Outputs:
#     - Generation time PMFs for each profile in GT_PROFILES
#     - National R(t) estimates via EpiNow2 (not EpiEstim), with 60% & 90% bands
#     - Serial interval from contact-tracing pairs (if available)
#     - Onset-to-sample and onset-to-death delay distributions
#     - Delay-adjusted cCFR point estimate and CI
#
# Sources / citations (documented inline):
#   - Generation time: GT_PROFILES in 00_config.R
#   - EpiNow2 >= 1.5.0 API (dist_spec Gamma()/LogNormal()); Abbott et al. 2020 Wellcome Open Res
#   - Radiation model: Simini et al. 2012 Nature 484:96-100
#   - Delay params: pipelines/cfr_analyses/tables/
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(fitdistrplus)   # for discretised gamma MLE on SI
})
# EpiNow2 (Abbott et al. 2020; NOT EpiEstim) — estimate_rt_epinow2() calls several of its
# constructors BARE (generation_time_opts/trunc_opts/obs_opts/rt_opts/stan_opts), so it must be
# attached when present. Guard the attach so a missing EpiNow2 does not abort the whole script at
# source() time (R(t) is an optional diagnostic, not on the invasion-forecast path); the routine
# returns NULL early when it is unavailable.
.HAVE_EPINOW2 <- requireNamespace("EpiNow2", quietly = TRUE)
if (.HAVE_EPINOW2) suppressPackageStartupMessages(library(EpiNow2))

# ---------------------------------------------------------------------------
# 1. GENERATION TIME PMFs
# ---------------------------------------------------------------------------

#' Compute a discretised Gamma generation-time PMF (DOUBLE-INTERVAL CENSORED).
#'
#' Two clearly-separated steps span the pipeline; this function does STEP 1
#' (continuous Gamma -> DAILY pmf). STEP 2 (daily -> weekly aggregation) is
#' `daily_to_weekly_gt()` in 06_simple_models.R.
#'
#' STEP 1 — discretisation with interval censoring. Onset DATES are recorded only
#' to the calendar day, so the generation interval between two onset days is
#' interval-censored. We use the standard primary-censoring construction (Cori et
#' al. 2013, Am J Epidemiol 178:1505; Park/Abbott et al. `primarycensored`): with
#' the infector's unknown within-day onset time U ~ Uniform(0,1) and the
#' continuous interval D ~ Gamma, the recorded integer day-difference is
#' tau = floor(U + D), so
#'     p_day(tau) = P(floor(U + D) = tau)
#'                = \int_tau^{tau+1} F(v) dv  -  \int_{tau-1}^{tau} F(w) dw ,
#' i.e. the ONCE-INTEGRATED CDF difference — NOT the naive F(tau) - F(tau-1),
#' which treats the primary onset as observed at an exact time and biases the
#' short-lag weights. Computed exactly by
#' `primarycensored::dprimarycensored(pwindow = swindow = 1)` (empirically equal
#' to the integral above to < 1e-9); a base-R `integrate()` fallback is used when
#' `primarycensored` is not installed, so this carries no hard dependency.
#'
#' tau = 0 convention: the recorded day-difference CAN be 0, but the FOI
#' convolution (`compute_foi`) and `daily_to_weekly_gt()` are indexed on lag >= 1
#' day (generation interval biologically >= 1 day). We therefore DROP the tau = 0
#' bin and renormalise over tau = 1..max_tau — exactly as the > max_tau tail mass
#' is redistributed — which preserves the g[1..max_tau] contract every caller
#' relies on (unchanged return length and indexing).
#'
#' @param mean    Mean of the Gamma distribution (days).
#' @param sd      Standard deviation of the Gamma distribution (days).
#' @param max_tau Maximum generation time (days); PMF truncated + renormalised.
#' @param method  "censored" (default; double-interval / primary-censored) or
#'                "naive" (legacy single-interval F(tau)-F(tau-1); retained ONLY
#'                for the regression test and the before/after comparison in §1.5).
#' @return Numeric vector of length max_tau on tau = 1..max_tau, summing to 1.
#'         Natural Gamma params are shape = (mean/sd)^2, rate = mean/sd^2.
make_gt_pmf <- function(mean, sd, max_tau, method = c("censored", "naive")) {
  method <- match.arg(method)
  stopifnot(mean > 0, sd > 0, max_tau >= 1)

  # Method-of-moments Gamma parameters (natural params, reported per Charniga et al.)
  shape <- (mean / sd)^2
  rate  <- mean  / sd^2

  if (method == "naive") {
    # Legacy single-interval discretisation: g[tau] = P(tau-1 < GT <= tau).
    taus <- seq_len(max_tau)
    pmf  <- pgamma(taus, shape = shape, rate = rate) -
            pgamma(taus - 1L, shape = shape, rate = rate)
  } else {
    # Double-interval (primary-censored) DAILY pmf on tau = 0..max_tau, then drop tau = 0.
    pday <- .gt_daily_censored_pmf(0:max_tau, shape = shape, rate = rate)
    pmf  <- pday[-1L]                       # keep tau = 1..max_tau
  }

  # Normalise over the retained support (redistributes tau=0 and > max_tau mass).
  pmf <- pmf / sum(pmf)

  stopifnot(length(pmf) == max_tau, abs(sum(pmf) - 1) < 1e-9, all(pmf >= 0))
  pmf
}

#' Daily double-interval-censored (primary-censored) pmf of a Gamma delay.
#'
#' Returns P(floor(U + D) = tau) for tau in `taus0` (typically 0,1,...,max_tau),
#' with U ~ Uniform(0,1) and D ~ Gamma(shape, rate). Prefers the vetted
#' `primarycensored` implementation; falls back to the exact base-R integral
#'   p(tau) = \int_tau^{tau+1} F - \int_{max(tau-1,0)}^{tau} F   (F = Gamma CDF),
#' which was verified to equal `primarycensored` to < 1e-9.
#'
#' @param taus0  Non-negative integer support to evaluate.
#' @param shape,rate  Gamma natural parameters.
#' @return Numeric vector, same length as `taus0` (NOT renormalised here).
.gt_daily_censored_pmf <- function(taus0, shape, rate) {
  if (requireNamespace("primarycensored", quietly = TRUE)) {
    d <- tryCatch(
      primarycensored::dprimarycensored(
        taus0, pdist = stats::pgamma, shape = shape, rate = rate,
        pwindow = 1, swindow = 1, D = Inf),
      error = function(e) NULL)
    if (!is.null(d) && all(is.finite(d))) return(as.numeric(d))
  }
  # Base-R fallback (no external dependency).
  Fc   <- function(q) stats::pgamma(q, shape = shape, rate = rate)
  intF <- function(a, b) stats::integrate(Fc, a, b, subdivisions = 1000L,
                                          rel.tol = 1e-10)$value
  vapply(taus0, function(tau) {
    lo <- max(tau - 1, 0)                   # F(<0) = 0, so the lower integral vanishes at tau = 0
    intF(tau, tau + 1) - intF(lo, tau)
  }, numeric(1))
}

#' Compute all GT PMFs from GT_PROFILES.
#'
#' @return Named list of numeric PMF vectors, one per profile.
compute_all_gt_pmfs <- function() {
  pmfs <- lapply(names(GT_PROFILES), function(nm) {
    p  <- GT_PROFILES[[nm]]
    pmf <- make_gt_pmf(mean = p$mean, sd = p$sd, max_tau = p$max_tau)
    # Report NATURAL Gamma parameters (shape/rate) alongside mean/sd (Charniga et al.;
    # §5.6). Discretisation is double-interval censored (§1.5, make_gt_pmf).
    .shape <- (p$mean / p$sd)^2; .rate <- p$mean / p$sd^2
    message(sprintf(
      "[gt_pmf] %s (%s): mean=%.1f d, sd=%.1f d, shape=%.3f, rate=%.4f, max_tau=%d d, len=%d, sum=%.6f (double-interval censored)",
      nm, p$label, p$mean, p$sd, .shape, .rate, p$max_tau, length(pmf), sum(pmf)
    ))
    pmf
  })
  names(pmfs) <- names(GT_PROFILES)
  pmfs
}

# ---------------------------------------------------------------------------
# 1b. GENERATION-TIME PRIOR GRID  (review §2.1 — marginalise, don't select)
# ---------------------------------------------------------------------------

#' Monte-Carlo grid over the generation-time prior.
#'
#' Discretises a GT prior (GT_PRIOR by default) into a small grid of
#' (gt_mean, gt_sd) points with renormalised (truncated) Gaussian prior weights,
#' for MARGINALISING the invasion posterior over generation-time uncertainty
#' instead of selecting a single scenario. The grid spans +/- 2 prior SDs on each
#' axis (~95% of the prior mass), clipped to the prior's truncation bounds, with an
#' ODD number of points per axis so the prior mean is always a grid point.
#'
#' @param prior list like GT_PRIOR (mean_mu, mean_sd, sd_mu, sd_sd, mean_bounds,
#'   sd_bounds, n_grid_mean, n_grid_sd, max_tau).
#' @return data.frame(gt_mean, gt_sd, weight) sorted by descending weight;
#'   weight sums to 1.
gt_prior_grid <- function(prior = GT_PRIOR) {
  .span <- function(mu, s, n, bounds) {
    if (n <= 1L) return(mu)
    g <- seq(mu - 2 * s, mu + 2 * s, length.out = n)
    pmin(pmax(g, bounds[1]), bounds[2])
  }
  ms   <- unique(.span(prior$mean_mu, prior$mean_sd, prior$n_grid_mean, prior$mean_bounds))
  ss   <- unique(.span(prior$sd_mu,   prior$sd_sd,   prior$n_grid_sd,   prior$sd_bounds))
  grid <- expand.grid(gt_mean = ms, gt_sd = ss, KEEP.OUT.ATTRS = FALSE)
  # Renormalised (truncated-)Gaussian prior weight; truncation constant cancels on renorm.
  w    <- stats::dnorm(grid$gt_mean, prior$mean_mu, prior$mean_sd) *
          stats::dnorm(grid$gt_sd,   prior$sd_mu,   prior$sd_sd)
  grid$weight <- w / sum(w)
  grid[order(-grid$weight), , drop = FALSE]
}

#' Build the daily GT PMFs and prior weights for each grid point of a GT prior.
#'
#' Returns everything the marginalised Bayesian prediction needs, reusing the
#' existing `gt_pmfs[[key]]` interface: a named list of double-interval-censored
#' daily PMFs (one per grid point) and the matching renormalised prior weights.
#'
#' @param prior list like GT_PRIOR.
#' @return list(gt_pmfs = named list of daily PMF vectors, weights = named numeric
#'   summing to 1, grid = the gt_prior_grid() data.frame with an added `key`).
make_gt_prior_pmfs <- function(prior = GT_PRIOR) {
  g      <- gt_prior_grid(prior)
  g$key  <- sprintf("gtgrid%02d", seq_len(nrow(g)))
  pmfs   <- stats::setNames(
    lapply(seq_len(nrow(g)), function(i)
      make_gt_pmf(mean = g$gt_mean[i], sd = g$gt_sd[i], max_tau = prior$max_tau)),
    g$key)
  list(gt_pmfs = pmfs, weights = stats::setNames(g$weight, g$key), grid = g)
}

# ---------------------------------------------------------------------------
# 2. R(t) ESTIMATION — NATIONAL (EpiNow2)
# ---------------------------------------------------------------------------

#' Estimate national R(t) using EpiNow2 >= 1.4.0.
#'
#' Uses daily confirmed case counts dated by onset (coalesced to sample when onset
#' is missing). A RIGHT-TRUNCATION model (trunc_opts, from the onset->sample reporting
#' delay) is applied so the most recent, still-incomplete weeks are nowcast rather than
#' read as a genuine decline — without it R(t) drops spuriously below 1 at the tail.
#' Ascertainment is modelled as a latent scaling factor (obs_opts scale) rather than fixed.
#'
#' Caches results to outputs/diagnostics/epinow2_rt_{gt_profile_name}_{analysis_date}.rds.
#' Returns NULL (with a warning) if EpiNow2 fails or if there is insufficient
#' data.
#'
#' @param cases_df            Tibble with columns `date` (Date) and `confirm` (integer).
#'                            Must cover at least 14 days with at least 10 total cases.
#' @param gt_profile_name     Character string naming a key in GT_PROFILES.
#' @param generation_time_params List with elements `mean`, `sd` (days) — used to
#'                            parameterise EpiNow2's Gamma generation-time object.
#' @param analysis_date       Date used in the cache filename.
#' @return Tibble with columns: date, R_mean, R_lo_60, R_hi_60, R_lo_90, R_hi_90
#'         (the 60% [q20/q80] and 90% [q5/q95] credible bands — the pipeline
#'         convention, never 50%/95%). Returns NULL on failure.
estimate_rt_epinow2 <- function(cases_df,
                                gt_profile_name,
                                generation_time_params,
                                analysis_date = ANALYSIS_DATE) {
  stopifnot(
    is.data.frame(cases_df),
    all(c("date", "confirm") %in% colnames(cases_df)),
    inherits(cases_df$date, "Date"),
    is.character(gt_profile_name),
    gt_profile_name %in% names(GT_PROFILES),
    is.list(generation_time_params),
    all(c("mean", "sd") %in% names(generation_time_params))
  )
  if (!isTRUE(get0(".HAVE_EPINOW2", ifnotfound = FALSE))) {
    warning("[epinow2] Package 'EpiNow2' not installed; skipping national R(t) estimation.")
    return(NULL)
  }

  # EpiNow2 requires a strictly contiguous daily series. Case counts built via
  # count(date) omit zero-case days, leaving gaps that break EpiNow2's internal
  # date bookkeeping ("N items to be assigned to group of size N-1"). Pad the
  # series to every calendar day in range, filling absent days with zero.
  cases_df <- cases_df %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(confirm = sum(confirm, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(
      date = seq(min(date), max(date), by = "day"),
      fill = list(confirm = 0)
    ) %>%
    dplyr::arrange(date)

  # Minimum data checks
  n_days   <- nrow(cases_df)
  tot_cases <- sum(cases_df$confirm, na.rm = TRUE)
  if (n_days < 14 || tot_cases < 10) {
    warning(sprintf(
      "[epinow2] Insufficient data for R(t) estimation: %d days, %d total cases. Returning NULL.",
      n_days, tot_cases
    ))
    return(NULL)
  }

  # Cache path. The key MUST fingerprint the R(t) MODEL CONFIGURATION, not just (gt, date):
  # otherwise a cached rds computed under an OLD estimate_rt_epinow2 (e.g. before the
  # right-truncation model was added) is silently served for the same (gt, analysis_date),
  # masking the code change (this exact stale-cache bug served a spurious below-1 R(t) once).
  # RT_CACHE_VERSION is bumped whenever this function's model spec changes; the GT parameters
  # are folded into the name so different runs never collide on one key.
  .rt_ver <- get0("RT_CACHE_VERSION", ifnotfound = 2L)
  # Fingerprint EVERY config-level input to the R(t) model spec below, not just (gt, date):
  # the right-truncation delay rate (DELAY_ONSET_SAMPLE_RATE, used ~line 180) and the
  # ascertainment scale prior (ASCERTAINMENT_NOMINAL, ~line 199) both enter the fit, so a
  # change to either without a RT_CACHE_VERSION bump would otherwise silently serve a stale
  # cache — the very failure this key guards against.
  .asc     <- get0("ASCERTAINMENT_NOMINAL",   ifnotfound = 0.45)
  .os_rate <- get0("DELAY_ONSET_SAMPLE_RATE", ifnotfound = 0.228)
  cache_file <- file.path(
    OUT_DIAGNOSTICS,
    sprintf("epinow2_rt_%s_%s_v%d_gt%.0f-%.0f_a%.0f_d%.0f.rds", gt_profile_name,
            format(analysis_date, "%Y%m%d"), .rt_ver,
            generation_time_params$mean * 10, generation_time_params$sd * 10,
            .asc * 100, .os_rate * 1000)
  )
  if (file.exists(cache_file)) {
    message(sprintf("[epinow2] Loading cached R(t) from %s", basename(cache_file)))
    return(readRDS(cache_file))
  }

  gt_mean   <- generation_time_params$mean
  gt_sd     <- generation_time_params$sd
  gt_max    <- GT_PROFILES[[gt_profile_name]]$max_tau

  message(sprintf(
    "[epinow2] Running EpiNow2 for GT profile '%s' (mean=%.1f d, sd=%.1f d) on %d days of data",
    gt_profile_name, gt_mean, gt_sd, n_days
  ))

  result <- tryCatch({
    # EpiNow2 >= 1.5.0 API: use EpiNow2::Gamma() inside generation_time_opts()
    gt_obj <- generation_time_opts(
      EpiNow2::Gamma(mean = gt_mean, sd = gt_sd, max = gt_max)
    )

    # Right-truncation model: a case with a recent ONSET has usually not been sampled/
    # reported yet, so the last ~1-2 weeks of the onset series are under-counted. Feeding
    # that raw incomplete tail to EpiNow2 makes it infer a decline and R(t) falls below 1
    # spuriously. Model the onset->sample reporting delay (Exp(rate) fit => Gamma shape 1,
    # mean = 1/rate) as the truncation so EpiNow2 nowcasts the recent counts.
    .os_rate  <- get0("DELAY_ONSET_SAMPLE_RATE", ifnotfound = 0.228)
    .os_mean  <- 1 / max(.os_rate, 1e-6)
    trunc_obj <- tryCatch(
      trunc_opts(dist = EpiNow2::Gamma(mean = .os_mean, sd = .os_mean, max = 21L)),
      error = function(e) { warning("[epinow2] truncation model unavailable: ", conditionMessage(e)); trunc_opts() })

    fit <- EpiNow2::epinow(
      data = cases_df,
      generation_time = gt_obj,
      delays = delay_opts(),          # onset-dated input: no infection->onset delay needed
      truncation = trunc_obj,         # correct the right-truncated recent onset tail
      rt = rt_opts(
        # EpiNow2 >= 1.5 requires a dist_spec (not a list). LogNormal() takes the
        # prior mean/sd on the natural R scale; BDBV R0 ~1.5-2.5.
        prior = LogNormal(mean = 2.0, sd = 0.5)
      ),
      obs = obs_opts(
        # dist_spec required (EpiNow2 >= 1.5); ascertainment scale prior centred on the
        # project-wide nominal (ASCERTAINMENT_NOMINAL, currently 0.45), not a stray 0.5.
        scale = Normal(mean = get0("ASCERTAINMENT_NOMINAL", ifnotfound = 0.45), sd = 0.1)
      ),
      # Summarise at the pipeline's band convention: 60% (q20/q80) and 90% (q5/q95)
      # credible intervals, so R(t) reports lower_60/upper_60 and lower_90/upper_90
      # (NOT EpiNow2's default 50% band). Extracted as R_lo_60/R_hi_60 below.
      CrIs = c(0.6, 0.9),
      stan = stan_opts(
        backend = "cmdstanr",   # rstan backend fails to compile on this platform
        seed    = get0("RANDOM_SEED", ifnotfound = 20260704L),  # reproducible R(t) MCMC
        cores   = 4L,
        chains  = 4L,
        iter_sampling = 1000L,
        iter_warmup   = 500L
      ),
      verbose = FALSE
    )

    # Extract R(t) summary table. EpiNow2 1.9.0 made fit$estimates defunct;
    # summarised parameters now come from summary(fit, type = "parameters").
    # CRITICAL: summary() RE-summarises the posterior at ITS OWN default CrIs (0.2/0.5/0.9),
    # NOT the bands set on the fit — so it must be told CrIs = c(0.6, 0.9) here, otherwise it
    # emits lower_20/50/90 (no lower_60/upper_60) and the select() below errors, silently
    # NULL-ing R(t) on every run (after a full MCMC) and masquerading as "EpiNow2 failed".
    rt_raw <- tibble::as_tibble(summary(fit, type = "parameters", CrIs = c(0.6, 0.9)))
    if (is.null(rt_raw) || nrow(rt_raw) == 0) stop("EpiNow2 returned no summarised estimates.")

    rt_tbl <- rt_raw %>%
      dplyr::filter(variable == "R") %>%
      dplyr::select(date,
                    R_mean   = mean,
                    R_lo_60  = lower_60,   # 60% CI lower (q20) — pipeline band, not 50%
                    R_hi_60  = upper_60,   # 60% CI upper (q80)
                    R_lo_90  = lower_90,   # 90% CI lower (q5)
                    R_hi_90  = upper_90) %>%  # 90% CI upper (q95)
      dplyr::mutate(gt_profile = gt_profile_name)

    message(sprintf("[epinow2] R(t) estimated for %d dates.", nrow(rt_tbl)))
    saveRDS(rt_tbl, cache_file)
    rt_tbl

  }, error = function(e) {
    warning(sprintf("[epinow2] EpiNow2 failed for profile '%s': %s. Returning NULL.",
                    gt_profile_name, conditionMessage(e)))
    NULL
  })

  result
}

# ---------------------------------------------------------------------------
# 4. SERIAL INTERVAL FROM CONTACT-TRACING PAIRS
# ---------------------------------------------------------------------------

#' Estimate the serial interval distribution from contact-tracing transmission pairs.
#'
#' Matches contact-tracing records to the linelist on alert_id to recover both
#' source and secondary case onset dates.  Fits a CONTINUOUS Gamma to the empirical
#' SI values via maximum likelihood (fitdistrplus::fitdist, method = "mle"); the
#' fitted continuous Gamma is discretised downstream (make_gt_pmf) for use in the
#' renewal convolution. (Not an interval-censored/discrete-likelihood fit.)
#'
#' Returns NULL if fewer than 10 valid pairs are found.  This is expected early
#' in the outbreak when the contact-tracing database is sparse.
#'
#' @param contacts_df  Data frame with at least: `source_case_alert_id`, and
#'                     optionally `contact_id` (used as secondary ID if it
#'                     appears in linelist).
#' @param linelist_df  Data frame with at least: `alert_id`, `date_of_symptom_onset`.
#' @return List with elements: si_pmf (numeric vector, tau=0..30), shape, rate,
#'         mean, sd, n_pairs.  Or NULL if insufficient data.
estimate_si_from_contacts <- function(contacts_df, linelist_df) {
  stopifnot(
    is.data.frame(contacts_df),
    is.data.frame(linelist_df),
    "source_case_alert_id" %in% colnames(contacts_df),
    "alert_id"             %in% colnames(linelist_df),
    "date_of_symptom_onset" %in% colnames(linelist_df)
  )

  message("[serial_interval] Matching contact-tracing pairs to linelist...")

  # Prepare linelist lookup: alert_id → onset date
  ll_lookup <- linelist_df %>%
    dplyr::select(alert_id, onset_date = date_of_symptom_onset) %>%
    dplyr::mutate(onset_date = as.Date(onset_date)) %>%
    dplyr::filter(!is.na(onset_date)) %>%
    dplyr::distinct(alert_id, .keep_all = TRUE)

  # Source onset dates
  pairs <- contacts_df %>%
    dplyr::select(source_case_alert_id, contact_id) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      ll_lookup %>% dplyr::rename(source_onset = onset_date,
                                  source_alert  = alert_id),
      by = c("source_case_alert_id" = "source_alert")
    ) %>%
    dplyr::left_join(
      ll_lookup %>% dplyr::rename(contact_onset = onset_date,
                                  contact_alert  = alert_id),
      by = c("contact_id" = "contact_alert")
    ) %>%
    dplyr::filter(!is.na(source_onset), !is.na(contact_onset)) %>%
    dplyr::mutate(si_days = as.integer(contact_onset - source_onset)) %>%
    # Biological plausibility filter: SI in [0, 30] days
    dplyr::filter(si_days >= 0L, si_days <= 30L)

  n_pairs <- nrow(pairs)
  message(sprintf("[serial_interval] %d valid transmission pairs with SI in [0,30] d.", n_pairs))

  if (n_pairs < 10L) {
    warning(sprintf(
      "[serial_interval] Only %d pairs found (need >=10). Returning NULL.", n_pairs
    ))
    return(NULL)
  }

  # Fit a Gamma distribution to the empirical SI using MLE.
  # fitdist expects positive continuous values; SI=0 is given a small offset.
  si_vals <- pmax(pairs$si_days, 0.5)   # shift 0s to 0.5 for Gamma fit

  fit_gamma <- tryCatch(
    fitdistrplus::fitdist(si_vals, "gamma", method = "mle"),
    error = function(e) {
      warning("[serial_interval] Gamma MLE failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit_gamma)) return(NULL)

  shape <- unname(fit_gamma$estimate["shape"])
  rate  <- unname(fit_gamma$estimate["rate"])
  si_mean <- shape / rate
  si_sd   <- sqrt(shape) / rate

  message(sprintf(
    "[serial_interval] Fitted Gamma SI: shape=%.3f, rate=%.3f, mean=%.2f d, sd=%.2f d",
    shape, rate, si_mean, si_sd
  ))

  # Discretise the PMF over tau = 0, 1, ..., 30
  taus   <- 0:30
  si_pmf <- pgamma(taus + 0.5, shape = shape, rate = rate) -
            pgamma(pmax(taus - 0.5, 0), shape = shape, rate = rate)
  si_pmf[1] <- pgamma(0.5, shape = shape, rate = rate)  # bin [0, 0.5)
  si_pmf     <- si_pmf / sum(si_pmf)

  list(
    si_pmf  = si_pmf,
    shape   = shape,
    rate    = rate,
    mean    = si_mean,
    sd      = si_sd,
    n_pairs = n_pairs,
    pairs   = pairs
  )
}

# ---------------------------------------------------------------------------
# 6. MAIN — Compute and summarise all parameters
# ---------------------------------------------------------------------------

if (interactive()) {
message("\n=== 02_epi_params.R: Computing all epidemiological parameters ===\n")

# --- 6a. Generation time PMFs ---
message("--- Generation time PMFs ---")
GT_PMFS <- compute_all_gt_pmfs()

for (nm in names(GT_PMFS)) {
  p   <- GT_PROFILES[[nm]]
  pmf <- GT_PMFS[[nm]]
  # Verify moments approximately match requested parameters
  tau_seq      <- seq_len(length(pmf))
  pmf_mean     <- sum(tau_seq * pmf)
  pmf_sd       <- sqrt(sum(tau_seq^2 * pmf) - pmf_mean^2)
  message(sprintf(
    "  [check] %s PMF: requested mean=%.1f d, sd=%.1f d | PMF mean=%.2f d, sd=%.2f d",
    nm, p$mean, p$sd, pmf_mean, pmf_sd
  ))
}

# --- 6d. Serial interval (optional — requires contact-tracing data) ---
message("\n--- Serial interval (contact-tracing pairs) ---")

# Locate the most recent contact-tracing processed file
ct_dirs <- list.dirs(
  file.path(DATA_DIR, "processed", "contact_tracing_processed"),
  recursive = FALSE
)
ct_dirs <- ct_dirs[grepl("CONTACT-TRACING_", basename(ct_dirs))]

if (length(ct_dirs) == 0L) {
  message("[serial_interval] No contact-tracing directories found; skipping SI estimation.")
  SI_ESTIMATE <- NULL
} else {
  # Use the most recent directory (sort by name which encodes date)
  ct_latest <- sort(ct_dirs, decreasing = TRUE)[1L]
  ct_file   <- file.path(ct_latest, "contact_tracing_processed.csv")

  # Locate linelist (use Ituri lab linelist as primary)
  ll_file   <- file.path(DATA_DIR, "processed", "lab_linelist_processed",
                          "ituri_processed_linelist_latest.csv")

  if (!file.exists(ct_file)) {
    message(sprintf("[serial_interval] Contact-tracing file not found: %s", ct_file))
    SI_ESTIMATE <- NULL
  } else if (!file.exists(ll_file)) {
    message(sprintf("[serial_interval] Linelist file not found: %s", ll_file))
    SI_ESTIMATE <- NULL
  } else {
    message(sprintf("[serial_interval] Loading from %s", ct_latest))

    contacts_df <- readr::read_csv(ct_file, col_types = readr::cols(.default = "c"),
                                   show_col_types = FALSE)
    linelist_df <- readr::read_csv(ll_file, col_types = readr::cols(.default = "c"),
                                   show_col_types = FALSE)

    # Attempt to harmonise column naming for alert_id in linelist
    # Lab linelist uses kinshasa_lab_id or bunia_lab_id; contact tracing uses source_case_alert_id
    # Check whether linelist has an alert_id column; if not, derive from bunia_lab_id
    if (!"alert_id" %in% colnames(linelist_df)) {
      if ("bunia_lab_id" %in% colnames(linelist_df)) {
        linelist_df <- linelist_df %>%
          dplyr::mutate(alert_id = bunia_lab_id)
        message("[serial_interval] Using bunia_lab_id as alert_id proxy.")
      } else if ("kinshasa_lab_id" %in% colnames(linelist_df)) {
        linelist_df <- linelist_df %>%
          dplyr::mutate(alert_id = kinshasa_lab_id)
        message("[serial_interval] Using kinshasa_lab_id as alert_id proxy.")
      } else {
        message("[serial_interval] Linelist has no alert_id column; skipping SI estimation.")
        contacts_df <- NULL
      }
    }

    if (!is.null(contacts_df) && "source_case_alert_id" %in% colnames(contacts_df)) {
      SI_ESTIMATE <- tryCatch(
        estimate_si_from_contacts(contacts_df, linelist_df),
        error = function(e) {
          warning("[serial_interval] Error estimating SI: ", conditionMessage(e))
          NULL
        }
      )
    } else {
      message("[serial_interval] contact_tracing_processed.csv missing source_case_alert_id; skipping.")
      SI_ESTIMATE <- NULL
    }
  }
}

if (!is.null(SI_ESTIMATE)) {
  message(sprintf(
    "[serial_interval] SI estimate: mean=%.2f d, sd=%.2f d (n=%d pairs)",
    SI_ESTIMATE$mean, SI_ESTIMATE$sd, SI_ESTIMATE$n_pairs
  ))
} else {
  message("[serial_interval] No SI estimate available from contact tracing.")
}

# --- 6e. Summary printout ---
message("\n=== Parameter summary ===")
message(sprintf("GT profiles computed: %s", paste(names(GT_PMFS), collapse = ", ")))
message("=== 02_epi_params.R complete ===\n")
}
