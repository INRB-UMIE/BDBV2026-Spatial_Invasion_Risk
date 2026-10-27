# =============================================================================
# 02_epi_params.R — Epidemiological Parameter Estimation
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: the epidemiological-parameter layer. What this file actually produces:
#   - generation-time PMFs for each profile in GT_PROFILES, and the GT-prior grid
#     (make_gt_pmf / compute_all_gt_pmfs / gt_prior_grid / make_gt_prior_pmfs)
#   - national R(t) via EpiNow2, with 60% and 90% bands (estimate_rt_epinow2)
#   - the RIGHT-TRUNCATION of the R(t) case series, estimated per regime
#     (epinow2_truncation and the .trunc_* helpers)
#   - the daily confirmed series every R(t) consumer shares (epinow2_daily_confirmed)
#
# It does NOT produce an onset-to-death delay or a cCFR: this suite forecasts INVASION and
# projects confirmed cases, and models no deaths (see 00_config.R, where CFR_POINT_ESTIMATE
# was removed). The onset-to-SAMPLE delay is fitted in 04c_dhis2_delay_windows.R and
# resolved through effective_onset_sample_delay() (00_config.R).
#
# Sources / citations (documented inline):
#   - Generation time: GT_PROFILES in 00_config.R
#   - EpiNow2 >= 1.5.0 API (dist_spec Gamma()/LogNormal()); Abbott et al. 2020 Wellcome Open Res
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})
# EpiNow2 (Abbott et al. 2020) — estimate_rt_epinow2() calls several of its
# constructors BARE (generation_time_opts/trunc_opts/obs_opts/rt_opts/stan_opts), so it must be
# attached when present. The attach is guarded only so that a missing EpiNow2 does not abort the
# whole script at source() time; EpiNow2 is NOT optional for the pipeline. Both truncation
# regimes need it (the extract fit runs estimate_truncation(); the as-of fit returns an
# EpiNow2::Gamma), and those truncations drive the live and back-test nowcasts; the two-week
# projection's per-draw R (bayes_rt_week_draws(), 21_bayesian_renewal.R) is an EpiNow2 fit.
# Without it each routine here that needs it returns NULL, and its caller refuses or warns.
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
#'                for the regression test in tests/test_invasion_pipeline.R).
#' @return Numeric vector of length max_tau on tau = 1..max_tau, summing to 1.
#'         Natural Gamma params are shape = (mean/sd)^2, rate = mean/sd^2.
make_gt_pmf <- function(mean, sd, max_tau, method = c("censored", "naive")) {
  method <- match.arg(method)
  stopifnot(mean > 0, sd > 0, max_tau >= 1)

  # Method-of-moments Gamma parameters (natural params, reported as recommended by Charniga et
  # al. 2024, PLOS Comput Biol, "Best practices for estimating and reporting epidemiological
  # delay distributions of infectious diseases")
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
  # Carry the GENERATING parameters with the pmf. daily_to_weekly_gt() needs them to build the
  # WEEKLY kernel by censoring the Gamma at a 7-day window, which is a different distribution
  # from summing this daily pmf into 7-day blocks. Attached last, once the vector is final:
  # SUBSETTING drops non-`names` attributes (pday[-1L] above would have), while arithmetic such
  # as the normalisation keeps them. A caller that subsets the pmf (g[1:k]) loses them too and
  # gets the binned weekly kernel. `method` travels too, so a pmf whose daily
  # discretisation was explicitly requested as naive is never silently handed a censored weekly
  # kernel — the two would then disagree about which distribution the model is using.
  attr(pmf, "gt_mean")   <- mean
  attr(pmf, "gt_sd")     <- sd
  attr(pmf, "gt_method") <- method
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
    # Report NATURAL Gamma parameters (shape/rate) alongside mean/sd (Charniga et al. 2024;
    # see make_gt_pmf). Discretisation is double-interval censored (make_gt_pmf).
    .shape <- (p$mean / p$sd)^2; .rate <- p$mean / p$sd^2
    # REALISED moments of the PMF the models actually use. These need not equal the labelled
    # (mean, sd): make_gt_pmf() drops the tau = 0 cell and truncates at max_tau, then
    # renormalises over the retained support, and chopping the far right tail removes a
    # disproportionate share of the VARIANCE.
    #
    # AT THE DEPLOYED max_tau = 90 d THE GAP IS NEGLIGIBLE. Measured 2026-09-19 with the
    # shipped GT_PROFILES: short 12.003/6.511 against 12.0/6.5 (mean +0.03%, SD +0.17%),
    # medium 15.307/9.302 against 15.3/9.3 (+0.05%, +0.02%), long 18.000/10.493 against
    # 18.0/10.5 (+0.00%, -0.07%). (This comment previously described max_tau = 45 — "the
    # realised SD runs 4-8% below", "medium -> 14.98/8.63, SD -7.2%", "raise 45 to 80 d" —
    # which has not been the configured support since the profiles moved to 90.)
    #
    # The 2% warning below is therefore DORMANT BY DESIGN, not dead: it is the guard that
    # detects a future profile whose tail outruns its support, so the gap can never sit
    # silently between the cited parameter and the fitted model. It is retained for that.
    .tau <- seq_along(pmf)
    .m   <- sum(.tau * pmf)
    .s   <- sqrt(sum((.tau - .m)^2 * pmf))
    message(sprintf(
      "[gt_pmf] %s (%s): labelled mean=%.1f d, sd=%.1f d, shape=%.3f, rate=%.4f, max_tau=%d d, len=%d, sum=%.6f (double-interval censored) | REALISED mean=%.2f d (%+.1f%%), sd=%.2f d (%+.1f%%)",
      nm, p$label, p$mean, p$sd, .shape, .rate, p$max_tau, length(pmf), sum(pmf),
      .m, 100 * (.m / p$mean - 1), .s, 100 * (.s / p$sd - 1)
    ))
    if (abs(.s / p$sd - 1) > 0.02 || abs(.m / p$mean - 1) > 0.02)
      warning(sprintf(
        "[gt_pmf] '%s': the discretised/truncated PMF realises mean %.2f d, sd %.2f d against a labelled %.1f/%.1f — max_tau = %d d is too short for this tail. Raise GT_PROFILES$%s$max_tau (90 d realises all three shipped profiles to within 0.2%%) or report the realised moments.",
        nm, .m, .s, p$mean, p$sd, p$max_tau, nm), call. = FALSE)
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

#' Estimate national R(t) using EpiNow2 (>= 1.5.0 dist_spec API).
#'
#' Uses the daily confirmed case counts passed in; in this pipeline they are built by
#' epinow2_daily_confirmed(), i.e. dated by date_index exactly as the weekly counts. A
#' RIGHT-TRUNCATION model (trunc_opts) is applied so the most recent, still-incomplete weeks
#' are nowcast rather than read as a genuine decline — without it R(t) drops spuriously below
#' 1 at the tail. The truncation is FITTED per regime by epinow2_truncation() (section 4b),
#' not derived from the onset->sample delay.
#' There is NO ascertainment scale: a constant-in-time ascertainment is confounded with the
#' seeded infections and cannot move R(t) (Lison et al. 2024, PLOS Comput Biol
#' 10.1371/journal.pcbi.1012021), so estimating one buys nothing and invites the reader to
#' think a reporting fraction was inferred. See obs_opts() in the body.
#'
#' Caches results to OUT_DIAGNOSTICS as epinow2_rt_<profile>_<analysis date>_v<RT_CACHE_VERSION>_...rds,
#' the name keyed on the GT parameters, the fitted truncation, the R prior, the Stan seed and the
#' case series, with RT_CACHE_VERSION standing for the code (see the key below). Returns NULL
#' (with a warning) if EpiNow2 fails or if there is insufficient data.
#'
#' @param cases_df            Tibble with columns `date` (Date) and `confirm` (integer).
#'                            Must cover at least 14 days with at least 10 total cases.
#' @param gt_profile_name     Character string naming a key in GT_PROFILES.
#' @param generation_time_params List with elements `mean`, `sd` (days) and optionally
#'                            `max_tau` (days; GT_PROFILES[[gt_profile_name]]$max_tau when
#'                            absent) — used to parameterise EpiNow2's Gamma generation-time
#'                            object.
#' @param analysis_date       Date: with `issue_date`, it selects the truncation regime and the
#'                            truncation panel's end date (epinow2_truncation()); it also names
#'                            the cache file.
#' @param issue_date          NULL for the deployed fit; otherwise a back-test origin
#'                            (epinow2_truncation()).
#' @param ll                  the line list, needed by the "asof" truncation regime.
#' @return Tibble with columns: date, type (EpiNow2's estimate / estimate based on partial
#'         data / forecast), R_mean, R_lo_60, R_hi_60, R_lo_90, R_hi_90 (the 60% [q20/q80] and
#'         90% [q5/q95] credible bands — the pipeline convention, never 50%/95%), gt_profile.
#'         Returns NULL on failure.
estimate_rt_epinow2 <- function(cases_df,
                                gt_profile_name,
                                generation_time_params,
                                analysis_date = ANALYSIS_DATE,
                                issue_date = NULL,
                                ll = NULL) {
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

  # NO DATED ROW, NO SERIES. The padding below takes seq(min(date), max(date)), which errors on
  # an empty (or all-NA) date column instead of reaching the insufficient-data return, so that
  # case is caught here and handled exactly as documented. A row without a date cannot be placed
  # on a daily series at all; it is dropped, and said so.
  .undated <- is.na(cases_df$date)
  if (any(.undated))
    warning(sprintf("[epinow2] %d row(s) of cases_df have no date and are dropped.", sum(.undated)),
            call. = FALSE)
  cases_df <- cases_df[!.undated, , drop = FALSE]
  if (!nrow(cases_df)) {
    warning("[epinow2] Insufficient data for R(t) estimation: no dated case rows. Returning NULL.")
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
  .rt_ver <- get0("RT_CACHE_VERSION", ifnotfound = 8L)
  # Fingerprint EVERY config-level input to the R(t) model spec below, not just (gt, date):
  # the fitted right-truncation (resolved and hashed just below), the R prior, the GT support
  # and the Stan seed all enter the fit, so a change to any of them without a RT_CACHE_VERSION
  # bump would otherwise silently serve a stale cache — the failure this key guards against.
  # ASCERTAINMENT IS NOT PART OF THE SPEC (see obs_opts below: a constant scale cannot move
  # R and is no longer passed), so it must not key the cache. This literal is a FROZEN
  # FILENAME COMPONENT, not a model parameter: it keeps the existing cache file names
  # addressable so fits computed under an identical spec are still found. Do not read it as
  # an ascertainment rate — this pipeline estimates none.
  .cache_name_salt <- 45
  # RESOLVE THE TRUNCATION FIRST — it is part of the model spec, so the cache key depends on
  # it. Resolving before the cache check costs nothing after the first run: .trunc_fit() caches
  # its own Stan fit to disk, so this is a file read on every subsequent call.
  .trunc_dist <- epinow2_truncation(issue_date = issue_date, ll = ll,
                                    analysis_date = analysis_date)

  # THE TRUNCATION IS NO LONGER A CONFIG CONSTANT, so it cannot be fingerprinted by reading
  # one. It is fitted from data (the vintage archive, or the line list for a fold), so the key
  # carries a hash of the FITTED DISTRIBUTION ITSELF. A version bump alone would not do: two
  # runs on the same code and the same analysis date can legitimately see different archives.
  #
  # HASH THE OBJECT, NOT ITS str(). The first version of this key hashed
  # capture.output(str(.trunc_dist)), which is wrong twice over:
  #   * str() prints numbers at 3 SIGNIFICANT DIGITS (the fitted meanlog 2.0280 prints as
  #     "2.03"), so two refits differing by up to ~0.005 in meanlog — ~0.05 d of mean lag, and
  #     any change at all in the parameter SDs, which print as "0.022" — produced the IDENTICAL
  #     key and a stale R(t) fit would be served for a different truncation.
  #   * str()'s formatting follows the session's strOptions()/digits, so the key was not even
  #     reproducible across environments: a user with different options refits everything.
  # rlang::hash() serialises the object, so it separates meanlog at 1e-9, carries the `max`
  # attribute, and distinguishes the uncertain dist from its fix_parameters() collapse — which
  # matters here because trunc_opts() is given the UNCERTAIN one and the posterior depends on
  # it. Verified stable across separate R sessions and under options(digits = 3).
  .trunc_k <- if (is.null(.trunc_dist)) "none" else substr(rlang::hash(.trunc_dist), 1L, 12L)
  # Two further model-spec inputs the comment above CLAIMED to fingerprint but did not:
  #   * EPINOW2_R_PRIOR  -> rt_opts(prior = ) : changing the R prior changes the posterior.
  #   * GT max_tau       -> Gamma(max = )     : changing the GT support changes the renewal
  #                                             convolution. This gap was covered only by luck
  #                                             (max_tau 45->90 happened to land in the same
  #                                             commit as an unrelated cache-version bump).
  # .num1() rather than %||% (00_config.R): besides a NULL it also rejects a non-scalar or
  # non-finite value, which %||% would pass straight into sprintf() and the fit.
  .rpri    <- get0("EPINOW2_R_PRIOR", ifnotfound = list(mean = 2.0, sd = 0.5))
  .num1    <- function(x, default) {
    if (is.null(x) || length(x) != 1L || !is.finite(suppressWarnings(as.numeric(x)))) default
    else as.numeric(x)
  }
  # ONE GT support for the key AND the fit (gt_max below): the caller's max_tau, else the
  # profile's. The key used to read the caller's (0 when absent) while the fit read the
  # profile's, so the two could describe different supports.
  .gtmax   <- .num1(generation_time_params$max_tau, GT_PROFILES[[gt_profile_name]]$max_tau)
  # THE STAN SEED. .epinow2_rt_fit() seeds the sampler with RANDOM_SEED, so two seeds give two
  # posteriors; without it in the key a run varying RANDOM_SEED alone was served the other
  # seed's posterior. (bayes_rt_week_draws() and the truncation fingerprint key it already.)
  .seed    <- as.integer(get0("RANDOM_SEED", ifnotfound = 20260704L))
  # FINGERPRINT THE CASE SERIES ITSELF. Every other term here describes the MODEL SPEC; none
  # described the DATA being fitted, even though `cases_df` is what the model is estimated on.
  # cases_df comes from epinow2_daily_confirmed(dat$ll, ...), which depends on the onset
  # imputation (each record's draw is a hash of (alert_id, ONSET_IMPUTE_SEED) inside
  # load_linelist(), with no RNG stream), on ONSET_MODE, on APPEND_SITREP_CONFIRMED, and on the
  # line-list snapshot. An imputation replicate (a different ONSET_IMPUTE_SEED) re-dates the
  # imputed onsets; before this key existed a replicate found the previous replicate's rds and
  # returned its posterior verbatim, giving R(t) exactly zero between-replicate variation.
  # RT_CACHE_VERSION is a manual discipline for code changes; it cannot see the data.
  .data_k <- local({
    cc <- cases_df[order(cases_df$date), , drop = FALSE]
    v  <- suppressWarnings(as.numeric(cc$confirm))
    v[!is.finite(v)] <- -1
    # Cheap order-sensitive digest of (n, first/last date, every count) in double arithmetic.
    h <- length(v) %% 2147483647
    h <- (h * 31 + as.numeric(min(cc$date, na.rm = TRUE))) %% 2147483647
    h <- (h * 31 + as.numeric(max(cc$date, na.rm = TRUE))) %% 2147483647
    for (x in v) h <- (h * 31 + x) %% 2147483647
    h
  })
  cache_file <- file.path(
    OUT_DIAGNOSTICS,
    # "_ini": the fit's MCMC initial values are seeded (below); a file written before that was
    # not reproducible from the seed in its name and must not be served.
    sprintf("epinow2_rt_%s_%s_v%d_gt%.0f-%.0f-m%.0f_a%.0f_t%s_r%.0f-%.0f_s%d_y%09.0f_ini.rds",
            gt_profile_name,
            format(analysis_date, "%Y%m%d"), .rt_ver,
            generation_time_params$mean * 10, generation_time_params$sd * 10,
            .gtmax,
            .cache_name_salt * 100, .trunc_k,
            .num1(.rpri$mean, 2.0) * 100, .num1(.rpri$sd, 0.5) * 100,
            .seed, .data_k)
  )
  if (file.exists(cache_file)) {
    message(sprintf("[epinow2] Loading cached R(t) from %s", basename(cache_file)))
    return(readRDS(cache_file))
  }

  gt_mean   <- generation_time_params$mean
  gt_sd     <- generation_time_params$sd
  gt_max    <- .gtmax

  message(sprintf(
    "[epinow2] Running EpiNow2 for GT profile '%s' (mean=%.1f d, sd=%.1f d) on %d days of data",
    gt_profile_name, gt_mean, gt_sd, n_days
  ))

  # EpiNow2 draws its MCMC initial values from R's global RNG, so the fit is seeded with the seed
  # the cache name records, and the caller's stream is restored afterwards (a cache miss and a
  # cache hit then leave every later random number identical), as in bayes_rt_week_draws().
  .old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(.old_seed)) assign(".Random.seed", .old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set_fixed_seed(.seed)
  result <- tryCatch({
    fit <- .epinow2_rt_fit(cases_df, gt_mean, gt_sd, gt_max, truncation = .trunc_dist)

    # Extract R(t) summary table. EpiNow2 1.9.0 made fit$estimates defunct;
    # summarised parameters now come from summary(fit, type = "parameters").
    # CRITICAL: summary() RE-summarises the posterior at ITS OWN default CrIs (0.2/0.5/0.9),
    # NOT the bands set on the fit — so it must be told CrIs = c(0.6, 0.9) here, otherwise it
    # emits lower_20/50/90 (no lower_60/upper_60) and the select() below errors, silently
    # NULL-ing R(t) on every run (after a full MCMC) and masquerading as "EpiNow2 failed".
    rt_raw <- tibble::as_tibble(summary(fit, type = "parameters", CrIs = c(0.6, 0.9)))
    if (is.null(rt_raw) || nrow(rt_raw) == 0) stop("EpiNow2 returned no summarised estimates.")

    # CARRY `type`. EpiNow2::epinow() is called without a `forecast =` argument, so
    # forecast_opts(horizon = 7) applies and summary(type = "parameters") returns SEVEN DAYS
    # OF PROJECTION past the data, tagged in the `type` column as "forecast" (vs "estimate" /
    # "estimate based on partial data"). Dropping that column and keeping every row is how
    # run_all.R came to report tail(R_mean, 1) — a 7-day-ahead projection, with its wider
    # forecast band — as the "Primary R(t) estimate".
    # Consumers filter on `type` (or on date <= analysis_date) to get estimates only.
    rt_tbl <- rt_raw %>%
      dplyr::filter(variable == "R") %>%
      dplyr::select(date,
                    type,                  # estimate | estimate based on partial data | forecast
                    R_mean   = mean,
                    R_lo_60  = lower_60,   # 60% CI lower (q20) — pipeline band, not 50%
                    R_hi_60  = upper_60,   # 60% CI upper (q80)
                    R_lo_90  = lower_90,   # 90% CI lower (q5)
                    R_hi_90  = upper_90) %>%  # 90% CI upper (q95)
      dplyr::mutate(gt_profile = gt_profile_name)

    message(sprintf("[epinow2] R(t): %d dates (%d estimate, %d forecast).", nrow(rt_tbl),
                    sum(rt_tbl$type != "forecast", na.rm = TRUE),
                    sum(rt_tbl$type == "forecast", na.rm = TRUE)))
    saveRDS(rt_tbl, cache_file)
    rt_tbl

  }, error = function(e) {
    # The message names the CONDITION, not a presumed culprit: the most common cause is now a
    # refused truncation, in which case EpiNow2 never ran at all and "EpiNow2 failed" sends the
    # reader to the wrong place.
    warning(sprintf("[epinow2] R(t) could not be estimated for profile '%s': %s. Returning NULL.",
                    gt_profile_name, conditionMessage(e)))
    NULL
  })

  result
}

#' Daily national CONFIRMED cases for EpiNow2, dated EXACTLY as the weekly case counts.
#'
#' Dating is `date_index` (load_linelist(), 01_data_prep.R): the symptom onset where it is
#' usable, otherwise an onset imputed as sample date minus a sample-to-onset interval drawn from
#' same-week donors (01_data_prep.R), clamped to the outbreak floor. This is the dating the
#' fitted right-truncation describes (02_epi_params.R, SERIES BASIS), which treats every case as
#' onset-dated. Dating onset-less cases at their SAMPLE date instead (the former
#' convention) presents already-complete cases to that model as incomplete, so the recent tail
#' was inflated a second time and R biased upward (7 Sep 2026 snapshot, weekly R for 1-7 Sep:
#' 1.50 with sample-date dating vs 0.82 with this dating). Records without a date_index (the
#' complete_case onset mode) are excluded, exactly as in the weekly counts.
#'
#' Shared by the R(t) diagnostic (run_all.R) and the Bayesian suite's R draws
#' (bayes_rt_week_draws(), 21_bayesian_renewal.R).
#'
#' @param ll line list (load_linelist() output; needs `confirmed` and `date_index`).
#' @param end_date last day of the series (inclusive); trailing days without cases are zeros.
#' @param issue_date optional as-of date: only records observable by it are kept, using the
#'   observation-date rule of the weekly training counts (linelist_observation_date(),
#'   22_daily_reissue.R).
#' @param caller label for diagnostics.
#' @return gapless daily data.frame(date, confirm) from the first dated case to `end_date`;
#'   zero rows when no confirmed case qualifies.
epinow2_daily_confirmed <- function(ll, end_date, issue_date = NULL, caller = "epinow2") {
  stopifnot(is.data.frame(ll))
  miss <- setdiff(c("confirmed", "date_index"), names(ll))
  if (length(miss))
    stop(sprintf("[%s] the line list lacks %s, so cases cannot be dated as the weekly counts are.",
                 caller, paste(miss, collapse = ", ")), call. = FALSE)
  end_date <- as.Date(end_date)
  stopifnot(length(end_date) == 1L, !is.na(end_date))
  keep <- ll$confirmed %in% TRUE
  if (!is.null(issue_date)) {
    issue_date <- as.Date(issue_date)
    obs  <- linelist_observation_date(ll, issue_date, caller = caller)
    keep <- keep & !is.na(obs) & obs <= issue_date
  }
  d <- as.Date(ll$date_index[keep])
  d <- d[!is.na(d) & d <= end_date]
  if (!length(d))
    return(data.frame(date = as.Date(character(0)), confirm = integer(0)))
  days <- seq(min(d), end_date, by = "day")
  data.frame(date = days,
             confirm = tabulate(as.integer(d - min(d)) + 1L, nbins = length(days)))
}


# ---------------------------------------------------------------------------
# 4b. RIGHT-TRUNCATION FOR THE R(t) FIT  (estimated, not assumed)
# ---------------------------------------------------------------------------
# The series EpiNow2 is fitted to is right-truncated: recent onsets have occurred but have not
# yet reached the extract. EpiNow2 corrects for that if it is told the distribution of the lag.
# It used to be told the onset->SAMPLE delay, which is the wrong leg — see RT_CACHE_VERSION v8
# in 00_config.R for the measurements. These functions estimate the right one from data.
#
# TWO REGIMES (see 00_config.R): "extract" for the deployed fit (and a snapshot round's), "asof"
# for a re-censored back-test round's. The extract regime builds a list of (date, confirm)
# frames, oldest first, as EpiNow2::estimate_truncation() requires, each frame the daily
# confirmed series AS IT STOOD at that vintage. The as-of regime builds no frames: its lag is
# recorded on every record, so it is measured per record (.trunc_fit_asof_cohort()).
#
# SERIES BASIS. Both consumers correct the DATE_INDEX series: the recorded onset where usable,
# otherwise an onset imputed from the sample date (aggregate_to_zone_week() for the nowcast,
# epinow2_daily_confirmed() for R(t)). The truncation must therefore describe how complete THAT
# series is at each lag.
#   extract  date_index, by replaying every vintage through load_linelist() as of its own date
#            (load_linelist_vintage(), 01_data_prep.R; .trunc_frame_vintage() below). Until
#            2026-10-14 this panel used RECORDED onsets, and that is not a small approximation
#            here: when a snapshot is taken only about a third of the newest week's confirmed
#            records carry a recorded onset (about 70% eventually do); the rest are in the
#            date_index series already, at an imputed onset, but were absent from the
#            recorded-onset panel until their onset was entered. So the newest week of a
#            snapshot was far more complete on the date_index basis than the recorded-onset fit
#            assumed (about 0.19 on the 2026-10-09 frame), and the live nowcast multiplied that
#            week about twice too much; the live R(t) was corrected by the same distribution.
#            52_vintage_level_check.R compares, every run, the newest-week completeness the
#            extract fit assumes with the one observed (diagnostics/vintage_level_check.csv).
#            The objection that once ruled the replay out -- imputed onsets would be re-drawn
#            with a later delay fit -- no longer arises: each vintage imputes from its OWN donors
#            (records sampled in the same week whose onset that vintage records), exactly as the
#            live pipeline imputes on its date, and a record's draw depends only on its id and
#            the seed. Imputed onsets that move between vintages (an onset later entered, or a
#            changed donor pool) are re-dating noise the negative-binomial observation model
#            absorbs like any other.
#   asof     date_index of the records whose recorded onset is USABLE (.trunc_fit_asof_cohort()):
#            the recorded onset exactly as the corrected series dates it, i.e. censored at the
#            specimen date when it was entered a day or two after it, and never a pre-outbreak
#            typo. The as-of lag is that date -> observation date, both on every record.
#            Records with an IMPUTED date_index are left out, because their lag there would be
#            their own imputation draw. .trunc_basis_guard() measures what leaving them out
#            costs on every run; 52_vintage_level_check.R compares the as-of nowcast's assumed
#            origin-week completeness with the observed one (diagnostics/vintage_level_check.csv).

#' Archive vintage dates, thinned and windowed.
#' @return Date vector, ascending, or length 0.
.trunc_vintages <- function(dir = get0("TRUNC_ARCHIVE_DIR", ifnotfound = NA_character_),
                            end, window_days = get0("TRUNC_WINDOW_DAYS", ifnotfound = 56L),
                            thin_days = get0("TRUNC_THIN_DAYS", ifnotfound = 7L)) {
  if (is.na(dir) || !dir.exists(dir)) return(as.Date(character(0)))
  f <- list.files(dir, pattern = "^LINELIST_[0-9]{8}$")
  if (!length(f)) return(as.Date(character(0)))
  d <- as.Date(sub("^LINELIST_", "", f), format = "%d%m%Y")
  d <- sort(unique(d[!is.na(d) & d <= end & d > end - window_days]))
  if (!length(d)) return(d)
  # Thin FORWARD FROM THE OLDEST, then always keep the newest: the newest vintage is the
  # reference the earlier ones are compared against, so it must never be thinned away.
  # INDEXED, not `for (x in d)`: iterating a Date vector UNCLASSES it, so `x` arrives as a
  # bare numeric and the date arithmetic below fails ("can only subtract from Date objects").
  # The same trap is flagged at 16_invasion_eval.R's fold loop.
  keep_i <- 1L
  if (length(d) > 1L) for (i in 2:length(d)) {
    if (as.numeric(d[i] - d[keep_i[length(keep_i)]]) >= thin_days) keep_i <- c(keep_i, i)
  }
  sort(unique(c(d[keep_i], d[length(d)])))
}

#' Daily confirmed series from ONE archive vintage, dated by date_index exactly as the live
#' series is.
#'
#' The vintage is replayed through the full loader as of its own date (load_linelist_vintage(),
#' 01_data_prep.R), so its date_index is what the live pipeline would have built on that day:
#' the recorded onset where usable, otherwise an onset imputed from the sample date, plus the
#' sitrep top-up rows. Confirmed records only, as in the counts the nowcast corrects
#' (aggregate_to_zone_week()) and the series the R(t) fit reads (epinow2_daily_confirmed()).
#' See SERIES BASIS above for why this is not the recorded-onset series.
.trunc_frame_vintage <- function(dir, vintage, start) {
  folder <- sprintf("LINELIST_%s", format(vintage, "%d%m%Y"))
  if (!file.exists(file.path(dir, folder, "dhis2_processed_linelist.csv"))) return(NULL)
  if (!exists("load_linelist_vintage", mode = "function"))
    stop("[trunc] load_linelist_vintage() is not loaded (01_data_prep.R): the extract panel ",
         "cannot be built on the date_index basis.", call. = FALSE)
  # load_linelist() reads LINELIST_DIR, so a different archive directory here would silently
  # build the panel from another archive than the one listed.
  if (!identical(normalizePath(dir), normalizePath(get0("LINELIST_DIR", ifnotfound = ""))))
    stop("[trunc] TRUNC_ARCHIVE_DIR must be LINELIST_DIR: vintages are replayed through ",
         "load_linelist(), which reads LINELIST_DIR.", call. = FALSE)
  # A vintage the loader cannot read is DROPPED and said so, as an unreadable file always was;
  # .trunc_validate_panel() then refuses the panel if too few vintages remain.
  vll <- tryCatch(load_linelist_vintage(folder, vintage), error = function(e) {
    message(sprintf("[trunc] extract: vintage %s dropped, it could not be loaded: %s",
                    folder, conditionMessage(e))); NULL })
  if (is.null(vll)) return(NULL)
  di <- suppressWarnings(as.Date(vll$date_index))
  keep <- (vll$confirmed %in% TRUE) & !is.na(di) & di >= start & di <= vintage
  .trunc_tabulate(di[keep], start, vintage)
}

#' Tabulate onset dates onto a GAPLESS daily spine. estimate_truncation() requires a contiguous
#' series; a count(date) with absent zero days breaks its internal date bookkeeping.
.trunc_tabulate <- function(dates, start, end) {
  start <- as.Date(start); end <- as.Date(end)
  # GUARD BEFORE seq(), not after. seq.Date() raises "wrong sign in 'by' argument" when
  # end < start, so the !length(days) check below could never be reached — the function
  # errored instead of returning NULL, and a single out-of-range vintage would have taken
  # down the whole panel build rather than being dropped from it.
  if (!is.finite(start) || !is.finite(end) || end < start) return(NULL)
  days <- seq(start, end, by = "day")
  if (!length(days)) return(NULL)
  idx <- as.integer(as.Date(dates) - as.Date(start)) + 1L
  idx <- idx[is.finite(idx) & idx >= 1L & idx <= length(days)]
  data.frame(date = days,
             confirm = as.integer(tabulate(idx, nbins = length(days))))
}

#' The onset date the AS-OF regime measures its lag from: date_index, for the records whose
#' recorded onset is USABLE (load_linelist(), 01_data_prep.R), NA for every other record.
#'
#' That is the recorded onset exactly as the corrected series dates it: censored at the specimen
#' date when it was entered a day or two after it, and never a pre-outbreak typo (those are
#' unusable, so the series dates them by imputation). The raw date_of_symptom_onset is NOT that
#' series: measuring from it drops the post-specimen onsets that the series holds at lag 0 and
#' admits the onsets the series replaces, so the fitted completeness described a different
#' series from the one corrected. Records with an IMPUTED date_index stay out: their lag would be
#' their own imputation draw (.trunc_basis_guard() measures what leaving them out costs).
#' Because only recorded onsets enter, the result does not depend on ONSET_IMPUTE_SEED.
#' @return Date vector, one per row of `ll`, or NULL when `ll` lacks the loader's columns.
.trunc_asof_onset <- function(ll) {
  miss <- setdiff(c("date_index", "onset_usable"), names(ll))
  if (length(miss)) {
    message(sprintf("[trunc] asof: the line list lacks %s (load_linelist()); cannot fit.",
                    paste(miss, collapse = ", ")))
    return(NULL)
  }
  ons <- suppressWarnings(as.Date(ll$date_index))
  ons[!(ll$onset_usable %in% TRUE)] <- as.Date(NA)
  ons
}

#' The as-of completeness curve, measured by ONSET COHORT, then fitted parametrically.
#'
#' WHY NOT estimate_truncation(). That deconvolves a lag from SNAPSHOTS, which is the only
#' option when the lag is recorded nowhere — true for the extract regime, where "appearance in
#' the extract" exists only as the difference between two archive vintages. Here the lag is
#' onset -> linelist_observation_date() and BOTH are on every record (the observation date is
#' the sample date for all but a handful of confirmed records: 9,046 of 9,054 on the 2026-10-09
#' frame), so nothing needs deconvolving. The onset is the date_index of the records with a
#' usable recorded onset (.trunc_asof_onset()).
#' An earlier version built deterministic pseudo-vintages and handed them to
#' estimate_truncation(): a record with observation <= v appears in EVERY later vintage and
#' never changes, so the panel was perfectly nested (0 of 721 shared day-cells non-monotone,
#' against 39 of 819 for the real archive), the negative-binomial dispersion had nothing to
#' estimate, and sampling took 129 s/chain before failing on CSV retrieval. No sampler setting
#' fixes an unidentified parameter.
#'
#' WHY NOT A RIGHT-TRUNCATION-CORRECTED MLE EITHER. That is the textbook answer and it is
#' measurably biased on this data. Each record is truncated at T = end - onset, so records with
#' small T carry weight 1/F(T); the parametric fit buys a heavy tail to pay for them, and the
#' tail steals mass from the short lags the anchor week is made of. Validated against onsets
#' old enough to be fully observed (>= TRUNC_MAX_DAYS before the as-of date), the anchor-week
#' completeness came out 13-25% LOW, and the bias fell monotonically to zero as the truncation
#' floor rose — i.e. it vanished exactly when the correction stopped doing anything. (The
#' likelihood itself is right: on simulated data with a known truncated LogNormal it recovered
#' mean 9.78 d against a truth of 9.92 d, where the naive estimate gave 8.55 d.)
#'
#' WHAT THIS DOES INSTEAD. The completeness curve is estimated NONPARAMETRICALLY with the
#' right-truncation-corrected reverse-hazard (Lynden-Bell) estimator,
#'     h(k) = #{delay = k, end - onset >= k} / #{delay <= k, end - onset >= k},
#'     F_hat(L) = prod_{k = L+1..TRUNC_MAX_DAYS} (1 - h(k)),
#' which needs no parametric family, and a gamma is then fitted to that curve by least squares
#' weighted by cohort size. (Until 2026-10-15 the curve was the plain cohort share
#' P(delay <= L | onset <= end - L), which is NOT free of truncation: records whose delay
#' exceeds end - onset are absent from the data, so it over-states completeness, by 1.5% at the
#' newest week on the 2026-10-09 frame and by 5-18% on synthetic data growing at 0.01-0.03/d.)
#'
#' THE CURVE IS FITTED AS THE MODEL DESCRIBES IT. F_hat lives on 0..TRUNC_MAX_DAYS with
#' F_hat(TRUNC_MAX_DAYS) = 1, and the returned Gamma is a distribution on that support,
#' renormalised there: its `max` is chosen so EpiNow2 discretises it on 0..TRUNC_MAX_DAYS
#' (.trunc_support_max()), and the nowcast divides by F(TRUNC_MAX_DAYS + 0.5)
#' (.trunc_as_delay_spec()). So F(L + 0.5) / F(TRUNC_MAX_DAYS + 0.5) is what is matched to
#' F_hat(L): F(L + 0.5) = P(delay <= L) is the daily-rounding convention
#' compute_truncation_weights() evaluates under, and the renormalisation is the one both
#' consumers apply, so the fitted spec reproduces the estimated completeness where it is used.
#'
#' GAMMA, NOT LOGNORMAL. Fitted to the same curve on the 2026-10-09 frame, gamma tracks it to a
#' maximum absolute error of 0.053 across the fitted lags against 0.090 for lognormal, and its
#' newest-week completeness is within 0.1% of the estimated one against 2.7% for lognormal.
#'
#' NO ROLLING WINDOW, unlike the extract panel. Windowing would starve the long lags (a 56-day
#' window leaves 11 days of onsets at lag 45). The cost, stated: the fit is NOT dominated by
#' recent cases. On the 2026-10-09 frame 72% of the lag-6 risk set (4,900 of 6,790 records) is
#' old enough to be fully observed, so the curve is close to an outbreak-wide average of the
#' onset -> sample delay rather than the delay at any one date. That is what this fit is for:
#' ONE shared nuisance estimate for every re-censored round (epinow2_truncation()), and those are
#' the early rounds the archive does not cover.
#'
#' CONFIRMED ONLY, deliberately: the folds correct the confirmed series, and the onset->sample
#' delay differs by classification (suspected cases are sampled on an alert-driven schedule), so
#' effective_onset_sample_delay() pooled over every classification is shorter than its
#' CONFIRMED stratum. Using the pooled delay is why the fold nowcast once under-corrected even
#' though its LEG was already right.
#'
#' @return an EpiNow2 dist_spec (Gamma), or NULL.
.trunc_fit_asof_cohort <- function(ll, end, start) {
  if (!exists("linelist_observation_date", mode = "function")) return(NULL)
  ons <- .trunc_asof_onset(ll)
  if (is.null(ons)) return(NULL)
  obs <- tryCatch(linelist_observation_date(ll, as.Date(end), caller = "trunc-asof"),
                  error = function(e) NULL)
  if (is.null(obs)) return(NULL)
  end <- as.Date(end); start <- as.Date(start)
  keep <- (ll$confirmed %in% TRUE) & !is.na(obs) & !is.na(ons) & ons >= start &
          as.Date(obs) <= end
  o  <- ons[keep]
  dl <- as.numeric(as.Date(obs[keep]) - o)
  # A NEGATIVE lag counts at lag 0, as the corrected series counts it. On this basis the
  # usual source (an onset entered a day or two after its own specimen) is already censored at
  # the specimen, so a remaining dl < 0 is a record observed before its date_index (an
  # observation date from a fallback column, or a specimen taken before the outbreak floor that
  # date_index is clamped to). reaggregate_asof() makes such a record visible from its
  # date_index onward, i.e. complete at lag 0, and that is where it belongs in the curve.
  dl <- pmax(dl, 0)
  mx <- as.integer(get0("TRUNC_MAX_DAYS", ifnotfound = 45L))
  # The model's support is 0..mx (the pmf is renormalised over it), so a delay beyond mx is
  # outside the distribution being estimated.
  ok <- is.finite(dl) & dl <= mx; o <- o[ok]; dl <- dl[ok]
  lags <- 0:mx
  # RIGHT TRUNCATION, CORRECTED EXACTLY. A record is in the as-of data only if its delay is at
  # most its truncation time Tt = end - onset, so recent onsets are represented only by their
  # short delays and a plain cohort share P(delay <= L | onset <= end - L) over-states
  # completeness (it misses the long delays still unobserved). The reverse hazard
  #     h(k) = P(delay = k | delay <= k)  =  #{dl = k, Tt >= k} / #{dl <= k, Tt >= k}
  # is estimable without bias from the truncated data (Lynden-Bell 1971; Lagakos et al. 1988),
  # and F(L) = prod_{k = L+1..mx} (1 - h(k)) with F(mx) = 1. (2026-10-09 frame, newest week:
  # cohort share 0.371, this estimator 0.366, onsets old enough to be fully observed 0.367; on
  # synthetic data the cohort share was biased by 5-18% at growth 0.01-0.03/d.)
  Tt  <- as.numeric(end - o)
  n_l <- vapply(lags, function(L) sum(Tt >= L), integer(1))
  h_k <- vapply(lags, function(k) {
    r <- dl <= k & Tt >= k
    if (!any(r)) 0 else sum(r & dl == k) / sum(r)
  }, numeric(1))
  f_l <- numeric(length(lags)); f_l[length(lags)] <- 1
  for (i in rev(seq_len(length(lags) - 1L))) f_l[i] <- f_l[i + 1L] * (1 - h_k[i + 1L])
  n_min <- get0("TRUNC_MIN_FRAME_CASES", ifnotfound = 200L)
  usable <- is.finite(f_l) & n_l >= n_min
  if (sum(usable) < 10L) {
    message(sprintf("[trunc] asof: only %d lag(s) have >= %d records; cannot fit.",
                    sum(usable), n_min))
    return(NULL)
  }
  # Weighted least squares on the CDF. Not a likelihood: the cohorts overlap, so their
  # contributions are not independent and a product of binomials would overstate precision.
  # Weighting by cohort size is what keeps the long, thin lags from dominating the short ones.
  # The model CDF is renormalised on 0..mx, as both consumers evaluate it (see the docstring).
  w <- as.numeric(n_l[usable]); fu <- f_l[usable]; lu <- lags[usable]
  cdf_on_support <- function(sh, rt)
    pmin(pgamma(lu + 0.5, shape = sh, rate = rt) / pgamma(mx + 0.5, shape = sh, rate = rt), 1)
  sse <- function(par) {
    sh <- exp(par[1]); rt <- exp(par[2])
    if (!is.finite(sh) || !is.finite(rt) || sh <= 0 || rt <= 0) return(1e10)
    cdf <- cdf_on_support(sh, rt)
    if (any(!is.finite(cdf))) return(1e10)
    sum(w * (cdf - fu)^2)
  }
  # Start from the method-of-moments gamma of the observed delays: feasible, and close enough
  # that Nelder-Mead does not wander into the flat region at large shape.
  .m <- mean(dl) + 0.5; .v <- max(stats::var(dl), 1e-6)
  fit <- tryCatch(stats::optim(c(log(max(.m^2 / .v, 1e-3)), log(max(.m / .v, 1e-4))),
                               sse, method = "Nelder-Mead",
                               control = list(maxit = 5000L, reltol = 1e-12)),
                  error = function(e) NULL)
  if (is.null(fit) || !identical(fit$convergence, 0L) || !is.finite(fit$value)) {
    warning("[trunc] asof: the cohort completeness fit did not converge; REFUSED.", call. = FALSE)
    return(NULL)
  }
  sh <- exp(fit$par[1]); rt <- exp(fit$par[2])
  err <- max(abs(cdf_on_support(sh, rt) - fu))
  message(sprintf(paste0("[trunc] asof: truncation-corrected completeness on %d confirmed records ",
                         "(lag-0 cohort %d, lag-%d cohort %d); gamma fit max |F - F_hat| = %.4f"),
                  length(dl), n_l[1], mx, n_l[mx + 1L], err))
  # A curve the family cannot follow is a wrong family, not a noisy fit — refuse rather than
  # correct every week with it. 0.08 lies between the realised gamma error and the lognormal
  # error the family comparison rejected (0.053 and 0.090 on the 2026-10-09 frame).
  if (!is.finite(err) || err > 0.08) {
    warning(sprintf("[trunc] asof: gamma cannot follow the measured completeness curve (max |F - F_hat| = %.3f); REFUSED.",
                    err), call. = FALSE)
    return(NULL)
  }
  # EPINOW2 MUST SEE THE SUPPORT THE CURVE WAS ESTIMATED ON, 0..mx. A fixed-parameter dist is
  # discretised one day shorter than its `max` (.trunc_support_max()), so the `max` that gives
  # 0..mx is found by measuring the support rather than assumed; a dist that cannot be given
  # it is refused.
  for (m in c(mx, mx + 1L)) {
    g <- tryCatch(EpiNow2::Gamma(shape = sh, rate = rt, max = m), error = function(e) NULL)
    if (!is.null(g) && identical(.trunc_support_max(g), as.numeric(mx))) return(g)
  }
  warning(sprintf("[trunc] asof: no Gamma(max = ) gives EpiNow2 the support 0..%d d; REFUSED.", mx),
          call. = FALSE)
  NULL
}

#' Extract-basis panel: one frame per archive vintage.
.trunc_panel_extract <- function(end, start,
                                 dir = get0("TRUNC_ARCHIVE_DIR", ifnotfound = NA_character_)) {
  vs <- .trunc_vintages(dir, end = end)
  if (length(vs) < 2L) return(NULL)
  lapply(vs, function(v) .trunc_frame_vintage(dir, v, start))
}

#' Reject a panel that cannot support a truncation fit. Refuses rather than approximates.
.trunc_validate_panel <- function(panel, regime) {
  bad <- function(why) { message(sprintf("[trunc] %s panel rejected: %s", regime, why)); NULL }
  if (is.null(panel)) return(bad("could not be built"))
  panel <- Filter(function(z) !is.null(z) && nrow(z) > 0L, panel)
  n_min <- get0("TRUNC_MIN_SNAPSHOTS", ifnotfound = 6L)
  if (length(panel) < n_min)
    return(bad(sprintf("%d usable vintage(s), need %d", length(panel), n_min)))
  # Ascending by last date, and strictly increasing in length — estimate_truncation() assumes
  # each snapshot extends the previous one.
  ends <- as.Date(vapply(panel, function(z) as.character(max(z$date)), character(1)))
  o <- order(ends); panel <- panel[o]; ends <- ends[o]
  if (anyDuplicated(ends)) return(bad("two vintages share an end date"))
  tot <- vapply(panel, function(z) sum(z$confirm), numeric(1))
  n_cases <- get0("TRUNC_MIN_FRAME_CASES", ifnotfound = 200L)
  if (tot[length(tot)] < n_cases)
    return(bad(sprintf("newest vintage carries %.0f cases, need %d", tot[length(tot)], n_cases)))
  starts <- as.Date(vapply(panel, function(z) as.character(min(z$date)), character(1)))
  if (length(unique(starts)) != 1L) return(bad("vintages do not share a start date"))
  # NOT CHECKED, deliberately: estimate_truncation() treats each snapshot as a prefix of the
  # final one, but individual day-cells can FALL between vintages through reclassification or
  # an onset correction (27 such cells on the live extract panel). They are absorbed by the
  # negative-binomial observation model as noise. Rejecting the panel over them would discard
  # a usable fit; asserting monotonicity would be false. Recorded so neither is attempted.
  panel
}

#' Measure, per run, what the AS-OF fit's record set costs (the extract regime is fitted on the
#' whole date_index series of each vintage; see SERIES BASIS).
#' The as-of fit measures its lag only on the records with a usable recorded onset
#' (.trunc_asof_onset()), while the series it corrects also holds every record with an IMPUTED
#' date_index. Both are dated by the same date_index; the guard compares the share observed
#' within each lag (k = 0..TRUNC_MAX_DAYS, as observed by `end`, plain shares) for the fit's
#' records against all confirmed records. No re-imputation, no RNG — so it measures the
#' record-set difference and nothing else.
.trunc_basis_guard <- function(ll, end, start, regime, tol = 0.05) {
  ons <- .trunc_asof_onset(ll)
  if (is.null(ons)) return(invisible(NA_real_))
  di  <- suppressWarnings(as.Date(ll$date_index))
  obs <- tryCatch(linelist_observation_date(ll, as.Date(end), caller = "trunc-guard"),
                  error = function(e) NULL)
  if (is.null(obs)) return(invisible(NA_real_))
  # THE TWO CURVES ARE TAKEN OVER DIFFERENT RECORD SETS, and that is the entire point. On one
  # common subset the two dates are EQUAL wherever a usable recorded onset exists, so a
  # comparison there could never fire (the first version of this guard did exactly that): the
  # records that make the difference are the ones whose date_index is imputed.
  #
  # RESTRICT TO THE RECORDS THIS REGIME ACTUALLY SEES (observable by `end`), so the number
  # describes the fit that was just made rather than the whole line list. A negative lag counts
  # at lag 0 in both curves, as it does in the fit.
  base  <- (ll$confirmed %in% TRUE) & !is.na(obs) & as.Date(obs) <= as.Date(end)
  k_fit <- base & !is.na(ons) & ons >= start
  k_all <- base & !is.na(di)  & di  >= start
  if (sum(k_fit) < 100L || sum(k_all) < 100L) return(invisible(NA_real_))
  lag_fit <- pmax(as.numeric(as.Date(obs[k_fit]) - ons[k_fit]), 0)
  lag_all <- pmax(as.numeric(as.Date(obs[k_all]) - di[k_all]), 0)
  ks <- 0:get0("TRUNC_MAX_DAYS", ifnotfound = 45L)
  z_fit <- vapply(ks, function(k) mean(lag_fit <= k), numeric(1))
  z_all <- vapply(ks, function(k) mean(lag_all <= k), numeric(1))
  d <- max(abs(z_fit - z_all), na.rm = TRUE)
  # DIRECTION MATTERS MORE THAN MAGNITUDE. Where z_all >= z_fit at every lag, the full series is
  # observed sooner after its date_index than the fit's records are, so the fitted truncation is
  # slightly too LONG and both consumers over-correct the recent weeks; the reverse means they
  # under-correct. The message says which, so the sign is reported rather than silently absorbed.
  .dir <- if (all(z_all >= z_fit - 1e-12)) " (one-signed: fitted truncation is too LONG, so both consumers over-correct)"
          else if (all(z_fit >= z_all - 1e-12)) " (one-signed: fitted truncation is too SHORT, so both consumers under-correct)"
          else ""
  if (is.finite(d) && d > tol)
    warning(sprintf(paste0("[trunc] %s: the usable-onset records the truncation is fitted on ",
                           "differ from the whole date_index series (imputed onsets included) by ",
                           "%.3f at worst (tolerance %.2f)%s. The fitted truncation may not ",
                           "describe the corrected series."), regime, d, tol, .dir), call. = FALSE)
  else
    message(sprintf("[trunc] %s: usable-onset records vs the whole date_index series agree to %.4f (tol %.2f)%s",
                    regime, d, tol, .dir))
  invisible(d)
}

#' Fingerprint a panel so a cached truncation fit cannot be served for a different one.
.trunc_fingerprint <- function(panel, regime) {
  # THE FULL DAILY SERIES of every frame, not a summary of it. The key used to be each frame's
  # first date, last date and total count, which re-dating cases leaves unchanged (a change to the
  # onset imputation moves cases between days, not in or out), so a fit made on one series was
  # served for another.
  parts <- rlang::hash(lapply(panel, function(z) list(format(z$date), as.integer(z$confirm))))
  # THE ESTIMATOR IS PART OF THE KEY, not just the panel. TRUNC_SERIES_VERSION covers the
  # panel DEFINITION, but changing the prior, the observation family, the sampler settings or
  # the seed changes the fit on an unchanged panel — and nothing else here would notice, so the
  # old fit would be served. Pinned literally so editing the estimator below without editing
  # this line is not possible without noticing.
  est <- sprintf("lnorm(N(0,1),N(1,1))|negbin|chains2|s2000|w500|seed%d",
                 get0("RANDOM_SEED", ifnotfound = 20260704L))
  paste(c(regime,
          sprintf("v%d", get0("TRUNC_SERIES_VERSION", ifnotfound = 1L)),
          sprintf("m%d", get0("TRUNC_MAX_DAYS", ifnotfound = 45L)),
          est, parts), collapse = "|")
}

#' Fit (and cache) the truncation distribution for one regime.
#' @return an EpiNow2 dist_spec, or NULL.
.trunc_fit <- function(panel, regime) {
  fp <- .trunc_fingerprint(panel, regime)
  key <- substr(rlang::hash(fp), 1L, 16L)
  cache <- file.path(OUT_DIAGNOSTICS, sprintf("epinow2_truncation_%s_%s.rds", regime, key))
  if (file.exists(cache)) {
    d <- tryCatch(readRDS(cache), error = function(e) NULL)
    if (!is.null(d)) {
      message(sprintf("[trunc] %s: cached fit %s", regime, basename(cache)))
      # THE GATES ARE RE-APPLIED TO A CACHED FIT. They are not part of the key (a gate change
      # does not change the fit, so it must not force a Stan refit), so without this a fit
      # accepted under looser TRUNC_MEAN_RANGE / TRUNC_Z3_MAX would still be served after they
      # were tightened. A refused cached fit is refused exactly as a fresh one would be.
      return(if (.trunc_sanity(d, regime)) d else NULL)
    }
  }
  mx <- get0("TRUNC_MAX_DAYS", ifnotfound = 45L)
  message(sprintf("[trunc] %s: fitting truncation on %d vintages (max %d d) — Stan, one-off",
                  regime, length(panel), mx))
  # estimate_truncation() draws its initial values from the GLOBAL RNG (create_initial_conditions
  # -> truncnorm/runif/rnorm). Seed it and restore the caller's stream, exactly as the R(t) fit
  # does, so a cache miss and a cache hit leave every downstream random number identical.
  .old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(.old)) assign(".Random.seed", .old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set_fixed_seed(get0("RANDOM_SEED", ifnotfound = 20260704L))
  est <- tryCatch(
    EpiNow2::estimate_truncation(
      panel,
      truncation = trunc_opts(EpiNow2::LogNormal(meanlog = EpiNow2::Normal(0, 1),
                                                 sdlog   = EpiNow2::Normal(1, 1),
                                                 max     = mx)),
      # A day-of-week term is NOT available here (estimate_truncation()'s obs_opts takes only
      # family and dispersion), so the fitted lag averages over the weekly pull rhythm.
      obs = obs_opts(family = "negbin"),
      stan = stan_opts(backend = "cmdstanr",
                       seed = get0("RANDOM_SEED", ifnotfound = 20260704L),
                       chains = 2L, cores = 2L, samples = 2000L, warmup = 500L),
      CrIs = c(0.5, 0.9), verbose = FALSE),
    error = function(e) { message("[trunc] ", regime, ": fit failed — ", conditionMessage(e)); NULL })
  if (is.null(est)) return(NULL)
  # NOT est$dist. `$.estimate_truncation` lifecycle::deprecate_stop()s on "dist" in
  # EpiNow2 >= 1.9.0 and names the replacement: get_parameters(x)[["truncation"]]. The old
  # accessor would abort the fit AFTER Stan had run.
  d <- tryCatch(EpiNow2::get_parameters(est)[["truncation"]], error = function(e) NULL)
  if (is.null(d)) {
    warning(sprintf("[trunc] %s: fit produced no truncation parameters; REFUSED.", regime),
            call. = FALSE)
    return(NULL)
  }
  if (!.trunc_sanity(d, regime)) return(NULL)
  # WRITE VIA tempfile() + file.rename(), as the rt-draws cache does. run_invasion_lfo() fans
  # models out with furrr::future_map(), so several workers can miss this cache at the same
  # moment; a bare saveRDS() can then be read half-written by another worker. run_all.R
  # pre-resolves both regimes before the parallel section, so in practice the file already
  # exists — this is the guard for every other entry point.
  tryCatch({
    tmp <- tempfile(tmpdir = dirname(cache), fileext = ".rds")
    saveRDS(d, tmp)
    if (!file.rename(tmp, cache)) { unlink(tmp); stop("rename failed") }
  }, error = function(e)
    warning("[trunc] could not cache the fit: ", conditionMessage(e), call. = FALSE))
  d
}

#' Refuse an implausible truncation rather than correct with it.
#' A wrong truncation moves EVERY published R(t) in one direction, so the gates are on the two
#' quantities that matter: the implied mean lag, and how complete it claims the last few days are.
.trunc_sanity <- function(d, regime) {
  # THE PMF CHAIN IS fix_parameters -> discretise -> get_pmf, in that order.
  #   * discretise() REFUSES a distribution with uncertain parameters ("Cannot discretise a
  #     distribution with uncertain parameters"), and estimate_truncation() returns exactly
  #     that, so fix_parameters() (which collapses each uncertain parameter to its mean) must
  #     come first. That is right for a SANITY CHECK; the fit itself keeps the uncertainty.
  #   * as.numeric() on a discretised dist_spec does NOT give the pmf — it returns length 1
  #     and NA. get_pmf() is the accessor.
  pmf <- tryCatch(EpiNow2::get_pmf(EpiNow2::discretise(EpiNow2::fix_parameters(d))),
                  error = function(e) NULL)
  if (is.null(pmf) || !is.numeric(pmf) || !length(pmf) || any(!is.finite(pmf)) || sum(pmf) <= 0) {
    warning(sprintf("[trunc] %s: fitted distribution has no usable pmf; REFUSED.", regime),
            call. = FALSE); return(FALSE)
  }
  pmf <- pmf / sum(pmf)
  lags <- seq_along(pmf) - 1L
  mean_lag <- sum(lags * pmf)
  z3 <- sum(pmf[lags <= 3L])
  rng <- get0("TRUNC_MEAN_RANGE", ifnotfound = c(4, 25))
  z3max <- get0("TRUNC_Z3_MAX", ifnotfound = 0.55)
  ok <- is.finite(mean_lag) && mean_lag >= rng[1] && mean_lag <= rng[2] && z3 <= z3max
  msg <- sprintf("implied mean %.2f d (allowed %.0f-%.0f), completeness at lag 3 %.3f (max %.2f)",
                 mean_lag, rng[1], rng[2], z3, z3max)
  if (!ok) warning(sprintf("[trunc] %s: REFUSED — %s.", regime, msg), call. = FALSE)
  else message(sprintf("[trunc] %s: accepted — %s", regime, msg))
  ok
}

#' The longest delay (days) EpiNow2 gives mass to when `dist` is the truncation passed to
#' trunc_opts(): the support on which it renormalises the truncation pmf.
#'
#' Mirrors EpiNow2's create_stan_delays() (verified on the installed 1.8.0.9000, and pinned by
#' tests/test_epi_params.R against that function). Both paths start from
#' discretise(fix_parameters(dist)), a pmf of length n:
#'   * UNCERTAIN parameters (the extract fit from estimate_truncation()) are passed to Stan
#'     parametrically with delay_max = n, and Stan evaluates the primary-censored pmf on 0..n,
#'     normalised by its CDF at n + 1;
#'   * FIXED parameters (the as-of Gamma, the legacy delay) are passed as that discretised pmf
#'     itself, on 0..n-1, normalised at n.
#' So the same `max` gives a support one day shorter on the fixed path; it is measured here,
#' never assumed.
#' @return numeric scalar, or NA when the distribution cannot be discretised.
.trunc_support_max <- function(dist) {
  n <- tryCatch(length(EpiNow2::get_pmf(EpiNow2::discretise(EpiNow2::fix_parameters(dist)))),
                error = function(e) NA_integer_)
  if (length(n) != 1L || !is.finite(n) || n < 1L) return(NA_real_)
  fixed <- !is.null(tryCatch(EpiNow2::discretise(dist), error = function(e) NULL))
  as.numeric(if (fixed) n - 1L else n)
}

#' Convert a fitted truncation dist_spec into the pipeline's own delay-spec shape.
#'
#' The NOWCAST (04_nowcasting.R) evaluates a delay through delay_cdf() (00_config.R), which
#' takes a list(family, params, rate, mean, sd) — not an EpiNow2 dist_spec. Both consumers must
#' describe the SAME truncation, so this converts rather than re-deriving.
#'
#' THE SUPPORT TRAVELS WITH THE SPEC (`trunc_max`). EpiNow2 renormalises the truncation pmf on
#' 0..trunc_max days (.trunc_support_max()), so the fitted parameters describe the delay
#' CONDITIONED on being at most trunc_max days; the extract fit is estimated under that
#' normalisation, and the as-of curve is fitted to it. Evaluated without it, the parametric CDF
#' leaves out the mass beyond trunc_max and reads every recent week as less complete than
#' EpiNow2 does. A consumer evaluating completeness at an integer lag L under the daily-rounding
#' convention must therefore use F(L + 0.5) / F(trunc_max + 0.5), capped at 1 (04_nowcasting.R).
#' `mean` and `sd` are the moments of that same conditioned delay (continuous delay below
#' trunc_max + 0.5), so the spec's summary describes the distribution both consumers use.
#'
#' The uncertain parameters are collapsed to their means (fix_parameters): the nowcast applies a
#' deterministic weight per week, so it has nowhere to put the uncertainty. That is a real
#' difference between the two consumers, stated here rather than hidden — the R(t) fit keeps the
#' posterior, the nowcast uses the point estimate.
#'
#' @return list in delay_cdf() shape plus `trunc_max`, or NULL if the distribution cannot be
#'   represented.
.trunc_as_delay_spec <- function(dist, regime = "unknown") {
  if (is.null(dist)) return(NULL)
  fx <- tryCatch(EpiNow2::fix_parameters(dist), error = function(e) NULL)
  if (is.null(fx)) return(NULL)
  pr <- tryCatch(EpiNow2::get_parameters(fx), error = function(e) NULL)
  if (is.null(pr) || !length(pr)) return(NULL)
  nm <- names(pr)
  tmax <- .trunc_support_max(dist)
  if (!is.finite(tmax)) {
    warning(sprintf("[trunc] %s: the truncation's support could not be determined; REFUSED.",
                    regime), call. = FALSE)
    return(NULL)
  }
  b <- tmax + 0.5                     # continuous bound of "rounded delay <= trunc_max"
  # Only the families delay_cdf() can evaluate: `lnorm` is the extract regime (estimate_truncation()
  # is fitted with a LogNormal prior family), `gamma` the as-of regime (.trunc_fit_asof_cohort()),
  # and the legacy delay takes its own family. Anything else returns NULL rather than silently
  # degrading to pexp().
  # Moments of the delay conditioned on D < b, in closed form:
  #   lognormal  E[D^j; D < b] = exp(j ml + j^2 sl^2 / 2) * pnorm((log b - ml - j sl^2) / sl)
  #   gamma      E[D^j; D < b] = (sh)_j / rt^j * pgamma(b, sh + j, rt)   ((sh)_j rising factorial)
  if (all(c("meanlog", "sdlog") %in% nm)) {
    ml <- as.numeric(pr$meanlog); sl <- as.numeric(pr$sdlog)
    if (!is.finite(ml) || !is.finite(sl) || sl <= 0) return(NULL)
    Fb <- stats::plnorm(b, ml, sl)
    m1 <- exp(ml + sl^2 / 2) * stats::pnorm((log(b) - ml - sl^2) / sl) / Fb
    m2 <- exp(2 * ml + 2 * sl^2) * stats::pnorm((log(b) - ml - 2 * sl^2) / sl) / Fb
    out <- list(family = "lnorm", params = c(meanlog = ml, sdlog = sl))
  } else if (all(c("shape", "rate") %in% nm)) {
    sh <- as.numeric(pr$shape); rt <- as.numeric(pr$rate)
    if (!is.finite(sh) || !is.finite(rt) || sh <= 0 || rt <= 0) return(NULL)
    Fb <- stats::pgamma(b, sh, rt)
    m1 <- sh / rt * stats::pgamma(b, sh + 1, rt) / Fb
    m2 <- sh * (sh + 1) / rt^2 * stats::pgamma(b, sh + 2, rt) / Fb
    out <- list(family = "gamma", params = c(shape = sh, rate = rt))
  } else return(NULL)
  if (!is.finite(Fb) || Fb <= 0 || !is.finite(m1) || m1 <= 0) return(NULL)
  out$mean <- m1
  out$sd   <- sqrt(max(m2 - m1^2, 0))
  out$trunc_max <- tmax
  # delay_cdf() falls back to pexp(q, dp$rate) when the named parameters are unusable, and
  # compute_truncation_weights() VALIDATES that `rate` is finite and positive, so it must be
  # present and sane even though the lnorm/gamma branch never reads it. 1/mean is the
  # moment-matched exponential rate — the correct fallback if it is ever reached.
  out$rate <- 1 / out$mean
  # THE LABEL NAMES THE ESTIMATOR THAT PRODUCED `dist`, as epinow2_truncation() dispatches it:
  # the legacy delay under TRUNC_SOURCE = "legacy" (either regime), else the Lynden-Bell cohort
  # gamma for "asof" and estimate_truncation() for "extract".
  out$estimator <-
    if (identical(get0("TRUNC_SOURCE", ifnotfound = "fitted"), "legacy"))
      sprintf("legacy_onset_sample_%s", regime)
    else if (identical(regime, "asof")) "lynden_bell_cohort_gamma_asof"
    else sprintf("estimate_truncation_%s", regime)
  out$source <- "fitted_truncation"
  out
}

#' The RETIRED onset->sample delay, expressed as an EpiNow2 dist_spec.
#'
#' TRUNC_SOURCE = "legacy" exists so the truncation change can be A/B'd against what the
#' pipeline used to do. It must therefore return a USABLE distribution: returning NULL would
#' hit .epinow2_rt_fit()'s refusal to run without one, so "legacy" would abort the run and the
#' rollback advice in that very error message would be circular. (It was, until 2026-09-21.)
#'
#' @return an EpiNow2 dist_spec, or NULL if the retired delay cannot be represented.
.trunc_legacy_dist <- function() {
  d <- tryCatch(effective_onset_sample_delay(), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mx <- get0("TRUNC_MAX_DAYS", ifnotfound = 45L)
  p  <- d$params
  nm <- names(p)
  tryCatch({
    if (identical(d$family, "gamma") && all(c("shape", "rate") %in% nm))
      EpiNow2::Gamma(shape = as.numeric(p[["shape"]]), rate = as.numeric(p[["rate"]]), max = mx)
    else if (identical(d$family, "lnorm") && all(c("meanlog", "sdlog") %in% nm))
      EpiNow2::LogNormal(meanlog = as.numeric(p[["meanlog"]]),
                         sdlog = as.numeric(p[["sdlog"]]), max = mx)
    else if (is.numeric(d$rate) && is.finite(d$rate) && d$rate > 0)
      # EpiNow2 exports no Exponential(); Gamma(shape = 1) IS the Exponential.
      EpiNow2::Gamma(shape = 1, rate = as.numeric(d$rate), max = mx)
    else NULL
  }, error = function(e) NULL)
}

#' THE ENTRY POINT. Truncation for one R(t) fit.
#'
#' @param issue_date NULL (or >= ANALYSIS_DATE) for the DEPLOYED fit, which reads the live
#'   extract and is truncated by onset -> appearance. A fold's own origin otherwise: that
#'   series is rebuilt by reaggregate_asof() and is truncated by onset -> sample.
#' @param ll the line list; required for the "asof" regime and for the basis guard.
#' @return an EpiNow2 dist_spec for trunc_opts(), or NULL (caller must decide).
#' @param regime NULL (default) derives the regime from the dates, which is what every R(t)
#'   call site wants. Pass "asof" explicitly to estimate the FOLD-basis truncation at the
#'   latest available date — the LFO shares ONE delay across folds (see run_invasion_lfo), so
#'   it needs the as-of estimate evaluated on all the data, which the date rule alone cannot
#'   express (issue_date = ANALYSIS_DATE would select the extract regime).
epinow2_truncation <- function(issue_date = NULL, ll = NULL,
                               analysis_date = get0("ANALYSIS_DATE", ifnotfound = Sys.Date()),
                               regime = NULL) {
  src <- get0("TRUNC_SOURCE", ifnotfound = "fitted")
  if (identical(src, "legacy")) {
    warning(paste0("[trunc] TRUNC_SOURCE='legacy': using the RETIRED onset->sample delay. That ",
                   "delay is not the lag that truncates either corrected series (onset -> ",
                   "appearance in the extract, or onset -> observation on the date_index basis), ",
                   "so the recent weeks and R(t) at the present are mis-corrected. For ",
                   "comparison only."), call. = FALSE)
    lg <- .trunc_legacy_dist()
    if (is.null(lg))
      warning("[trunc] TRUNC_SOURCE='legacy' but the retired delay could not be expressed as a ",
              "dist_spec; R(t) will refuse to fit.", call. = FALSE)
    return(lg)
  }
  start <- suppressWarnings(as.Date(get0("OUTBREAK_START", ifnotfound = NA)))
  if (is.na(start)) return(NULL)
  # REGIME: an explicit `regime =` when the caller knows it, otherwise DERIVED FROM THE DATES.
  # Never from an attribute on the data — the zero-padding pipeline in estimate_rt_epinow2()
  # drops non-standard attributes, so an attribute-borne regime would vanish silently.
  if (!is.null(regime)) {
    regime <- match.arg(regime, c("extract", "asof"))
    deployed <- identical(regime, "extract")
  } else {
    deployed <- is.null(issue_date) || as.Date(issue_date) >= as.Date(analysis_date)
    regime <- if (deployed) "extract" else "asof"
  }
  # END DATE — and the ASOF REGIME IS SHARED ACROSS FOLDS, deliberately.
  #
  # The extract panel ends at the issue date (or the analysis date when none is given).
  #
  # The AS-OF fit (a cohort estimator, no panel) always ends at the ANALYSIS DATE, whatever
  # issue_date says. The truncation
  # is a SHARED NUISANCE PARAMETER of the LFO — the same choice run_invasion_lfo() already
  # makes for its nowcast delay, and for the same reason: a nuisance parameter refitted per
  # fold is correlated with fold index, so apparent skill drifts with calendar time for a
  # purely artefactual reason. Ending every fold's as-of fit at the analysis date makes the fitted
  # truncation identical across folds by construction.
  #
  # It also removes a failure that per-fold fitting caused: THE EARLY FOLDS HAD NO TRUNCATION
  # AT ALL. A fold at issue 2026-05-18 sees 18 days of history, which was too little to
  # estimate from; the fit was refused, .epinow2_rt_fit() then refuses to run, BOTH retries
  # fail identically (the truncation is resolved once, outside the retry),
  # .rt_recent_posterior() cannot rescue it because rt_spec carries trunc="none", and the fold
  # silently fell back to the R0-scale LogNormal(2.0, 0.5) PRIOR. Measured on this snapshot
  # that was 3 of the 12 candidate cutoffs (issues 2026-05-18, -05-25, -06-01), and `delta` is
  # fitted on those folds and applied to the live forecast.
  #
  # THE COST, STATED: the as-of estimate then uses line-list rows that postdate a fold's origin.
  # That is look-ahead in the NUISANCE parameter — it is not a leak of the outcome being
  # scored, and it is the LFO's established policy for exactly this class of quantity, but it
  # is one more reason fold skill is optimistic relative to deployment. See the "RESIDUAL, FOR
  # THE RE-CENSORED ROUNDS ONLY" note in 04_nowcasting.R.
  end <- if (deployed && !is.null(issue_date)) as.Date(issue_date) else as.Date(analysis_date)

  # TWO REGIMES, TWO ESTIMATORS — because the lag is observable in one and not the other.
  #   extract : the appearance lag exists only as the difference between archive vintages,
  #             so it is DECONVOLVED from snapshots (estimate_truncation(), Stan, cached).
  #   asof    : the lag is onset (date_index of a usable recorded onset) ->
  #             linelist_observation_date(), both on every record, so the completeness curve is
  #             MEASURED by onset cohort and a gamma fitted to it (no Stan, no cache — it is
  #             deterministic and sub-second).
  if (deployed) {
    # No basis guard here: the extract panel IS the date_index series (SERIES BASIS above).
    panel <- .trunc_validate_panel(.trunc_panel_extract(end, start), regime)
    if (is.null(panel)) return(NULL)
    return(.trunc_fit(panel, regime))
  }
  if (is.null(ll)) {
    warning("[trunc] the 'asof' regime needs a line list (ll=); none was supplied.",
            call. = FALSE)
    return(NULL)
  }
  .trunc_basis_guard(ll, end, start, regime)
  d <- .trunc_fit_asof_cohort(ll, end, start)
  if (is.null(d)) return(NULL)
  if (!.trunc_sanity(d, regime)) return(NULL)
  d
}

#' Fit the national EpiNow2 renewal model behind every R(t) in the suite.
#'
#' Gamma generation time; right-truncation ESTIMATED per regime and supplied by the caller
#' (see `truncation` below — it is NOT derived from the onset->sample delay); LogNormal
#' R prior; NO ascertainment scale (a constant scale is confounded with seeded infections
#' and cannot move R — see obs_opts in the body). Factored out of
#' estimate_rt_epinow2() so the diagnostic R(t) panel and the Bayesian invasion suite's
#' posterior R draws (bayes_rt_week_draws(), 21_bayesian_renewal.R) fit the IDENTICAL
#' model and cannot drift apart. Errors propagate to the caller.
#'
#' @param truncation REQUIRED EpiNow2 dist_spec for the right-truncation, from
#'   epinow2_truncation(). NULL is refused rather than defaulted: see the body for why an
#'   un-truncated fit, or one on the retired onset->sample delay, would move every published
#'   R(t) with nothing in the output to show it.
#'
#' @param cases_df daily onset-dated confirmed counts (date, confirm), contiguous.
#' @param gt_mean,gt_sd generation-time mean and SD (days).
#' @param gt_max generation-time support (days).
#' @return the epinow() fit object.
.epinow2_rt_fit <- function(cases_df, gt_mean, gt_sd, gt_max, truncation) {
  # EpiNow2 >= 1.5.0 API: use EpiNow2::Gamma() inside generation_time_opts()
  gt_obj <- generation_time_opts(
    EpiNow2::Gamma(mean = gt_mean, sd = gt_sd, max = gt_max)
  )

  # ---- RIGHT TRUNCATION -------------------------------------------------------
  # The distribution is ESTIMATED per regime (section 4b above: from the archived vintages for
  # the extract regime, from the line list's onset cohorts for the as-of regime), not derived
  # from the onset->sample delay. That delay describes the wrong leg for the deployed series,
  # which is truncated by onset -> APPEARANCE IN THE EXTRACT, including an unrecorded data-entry
  # lag (sample->appearance measured at median 2.0 d, mean 2.6 d). While that delay was used
  # (until 2026-09-21) EpiNow2 treated the recent weeks as near-complete, read the reporting dip
  # as a decline, and returned an R(t) biased DOWN at the present — the value RT_WINDOW_WEEKS = 1
  # hands to the 2-week forecast and the cascade.
  #
  # `truncation` is supplied by the CALLER, which knows which series it is fitting:
  #   estimate_rt_epinow2()  -> deployed fit         -> extract regime (onset -> appearance)
  #   bayes_rt_week_draws()  -> snapshot round       -> extract regime fitted on that date
  #                          -> re-censored round    -> asof regime    (onset -> sample)
  # See epinow2_truncation() for why those differ and why using one for both is wrong.
  #
  # NULL is NOT silently tolerated. Falling back to an un-truncated fit would overstate the
  # recent decline even further, and falling back to the retired delay would reinstate the very
  # bias this replaced — either way every published R(t) would move with nothing in the output
  # to show it. Refuse, and let the caller's tryCatch report a failed R(t) honestly.
  if (is.null(truncation)) {
    stop("[epinow2] no truncation distribution was supplied. R(t) is NOT estimated without one: ",
         "an uncorrected fit reads the reporting lag as a decline. If the truncation fit failed, ",
         "fix it (or set TRUNC_SOURCE='legacy' deliberately, which warns) rather than running ",
         "without a truncation model.", call. = FALSE)
  }
  trunc_obj <- trunc_opts(dist = truncation)

  EpiNow2::epinow(
    data = cases_df,
    generation_time = gt_obj,
    delays = delay_opts(),          # onset-dated input: no infection->onset delay needed
    truncation = trunc_obj,         # correct the right-truncated recent onset tail
    rt = rt_opts(
      # EpiNow2 >= 1.5 requires a dist_spec (not a list). LogNormal() takes the
      # prior mean/sd on the natural R scale; EPINOW2_R_PRIOR (00_config.R) documents the choice.
      prior = LogNormal(mean = EPINOW2_R_PRIOR$mean, sd = EPINOW2_R_PRIOR$sd)
    ),
    obs = obs_opts(
      # NO ASCERTAINMENT SCALE. This used to pass scale = Normal(0.45, 0.1). With
      # rt_opts(pop = Fixed(0)) — EpiNow2's default, which this call keeps — a CONSTANT
      # scaling of the latent trajectory is exactly confounded with the seeded-infections
      # parameter, so the renewal posterior for R is invariant to it: it could not move
      # R(t) by construction, while being documented as a modelled feature and fingerprinted
      # into the R(t) cache key as though it were a model input. Lison et al. (2024, PLOS
      # Comput Biol 10.1371/journal.pcbi.1012021) make the same point directly: constant-in-
      # time ascertainment cannot bias R(t); only TIME-VARYING ascertainment can, and this
      # pipeline neither estimates nor validates a time-varying one. Default scale = Fixed(1).
      #
      # week_effect is stated EXPLICITLY rather than inherited. It is EpiNow2's default
      # (TRUE), and it is kept: recalled onset dates carry genuine weekday heaping, and R is
      # consumed downstream as a WEEKLY mean, which is exactly the quantity a within-week
      # effect should be removed from before averaging. Caveat, recorded rather than hidden:
      # ~24% of onsets are imputed (sample date minus a drawn delay), so the weekday effect
      # is estimated on a mixture of recalled and synthetic dates.
      week_effect = TRUE
    ),
    # Summarise at the pipeline's band convention: 60% (q20/q80) and 90% (q5/q95)
    # credible intervals, so R(t) reports lower_60/upper_60 and lower_90/upper_90
    # (NOT EpiNow2's default 50% band). Extracted as R_lo_60/R_hi_60 below.
    CrIs = c(0.6, 0.9),
    # SAMPLES/WARMUP, not iter_sampling/iter_warmup. EpiNow2::stan_opts() takes the TOTAL
    # post-warmup draws as `samples` and per-chain warmup as `warmup`; it DISCARDS
    # iter_sampling (warning at every fit: "Number of samples must be specified using the
    # `samples` and `warmup` arguments") and appends iter_warmup as a DUPLICATE name that
    # create_stan_args() then dedupes in favour of the package default. Verified against the
    # installed EpiNow2 1.8.0.9000: the old call ran at iter_warmup = 250, iter_sampling = 500
    # — half the intended warmup and half the draws — so every published R(t) rested on 2000
    # draws from a GP renewal model at adapt_delta 0.9, with the shortfall recorded nowhere.
    # samples = 4000 over 4 chains is 1000 post-warmup draws per chain, the documented intent.
    stan = stan_opts(
      backend = "cmdstanr",   # rstan backend fails to compile on this platform
      seed    = get0("RANDOM_SEED", ifnotfound = 20260704L),  # reproducible R(t) MCMC
      cores   = 4L,
      chains  = 4L,
      samples = 4000L,        # TOTAL post-warmup draws across chains (1000 per chain)
      warmup  = 1000L          # per-chain warmup
    ),
    verbose = FALSE
  )
}

# ---------------------------------------------------------------------------
# This file is a LIBRARY: sourcing it defines the estimators and computes nothing.
# ---------------------------------------------------------------------------
# An `if (interactive())` block here used to recompute every GT PMF and print its realised
# moments. It could not run under Rscript (interactive() is FALSE), so it was dead in every
# pipeline run, and compute_all_gt_pmfs() already reports and warns on those moments itself.
# Call compute_all_gt_pmfs() to get them.
