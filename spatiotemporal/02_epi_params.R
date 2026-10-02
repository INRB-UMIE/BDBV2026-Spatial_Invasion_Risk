# =============================================================================
# 02_epi_params.R — Epidemiological Parameter Estimation
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: Estimate and organise all epidemiological parameters required by
#   the spatiotemporal model suite.  Outputs:
#     - Generation time PMFs for each profile in GT_PROFILES
#     - National R(t) estimates via EpiNow2 (not EpiEstim), with 60% & 90% bands
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
# EpiNow2 (Abbott et al. 2020) — estimate_rt_epinow2() calls several of its
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
  # Carry the GENERATING parameters with the pmf. daily_to_weekly_gt() needs them to build the
  # WEEKLY kernel by censoring the Gamma at a 7-day window, which is a different distribution
  # from summing this daily pmf into 7-day blocks. Attached AFTER normalisation because R drops
  # non-`names` attributes through arithmetic. `method` travels too, so a pmf whose daily
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
    # Report NATURAL Gamma parameters (shape/rate) alongside mean/sd (Charniga et al.;
    # §5.6). Discretisation is double-interval censored (§1.5, make_gt_pmf).
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

#' Estimate national R(t) using EpiNow2 >= 1.4.0.
#'
#' Uses the daily confirmed case counts passed in; in this pipeline they are built by
#' epinow2_daily_confirmed(), i.e. dated by date_index exactly as the weekly counts. A RIGHT-TRUNCATION model (trunc_opts, from the onset->sample reporting
#' delay) is applied so the most recent, still-incomplete weeks are nowcast rather than
#' read as a genuine decline — without it R(t) drops spuriously below 1 at the tail.
#' There is NO ascertainment scale: a constant-in-time ascertainment is confounded with the
#' seeded infections and cannot move R(t) (Lison et al. 2024, PLOS Comput Biol
#' 10.1371/journal.pcbi.1012021), so estimating one buys nothing and invites the reader to
#' think a reporting fraction was inferred. See obs_opts() in the body.
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
  # the right-truncation delay (mean AND sd, resolved just below), the R prior and the GT
  # support all enter the fit, so a change to any of them without a RT_CACHE_VERSION bump
  # would otherwise silently serve a stale cache — the failure this key guards against.
  # Ascertainment is NOT part of the spec (see obs_opts below: a constant scale cannot move
  # R and is no longer passed), so it must not key the cache. The literal keeps the existing
  # file-naming scheme so fits cached under an identical spec stay addressable.
  .asc     <- 45
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
  # NOTE: do NOT use %||% here. It is not defined by 00_config.R, and 02_epi_params.R runs
  # (run_all.R step 2) BEFORE 03_mobility_matrices.R, which is where the suite's copy is
  # defined. Base R only gained %||% in 4.4, so relying on it would make this fail on 4.3.
  .rpri    <- get0("EPINOW2_R_PRIOR", ifnotfound = list(mean = 2.0, sd = 0.5))
  .num1    <- function(x, default) {
    if (is.null(x) || length(x) != 1L || !is.finite(suppressWarnings(as.numeric(x)))) default
    else as.numeric(x)
  }
  .gtmax   <- .num1(generation_time_params$max_tau, 0)
  # FINGERPRINT THE CASE SERIES ITSELF. Every other term here describes the MODEL SPEC; none
  # described the DATA being fitted, even though `cases_df` is what the model is estimated on.
  # cases_df comes from epinow2_daily_confirmed(dat$ll, ...), which depends on the onset
  # imputation (a stochastic per-record draw reseeded from RANDOM_SEED inside load_linelist()),
  # on ONSET_MODE, on APPEND_SITREP_CONFIRMED, and on the line-list snapshot. So varying
  # RANDOM_SEED for a replicate — the documented way to propagate imputation uncertainty —
  # found the previous replicate's rds and returned its posterior verbatim, giving R(t) exactly
  # zero between-replicate variation while rt_national.pdf claimed to describe the new draw.
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
    sprintf("epinow2_rt_%s_%s_v%d_gt%.0f-%.0f-m%.0f_a%.0f_t%s_r%.0f-%.0f_y%09.0f.rds",
            gt_profile_name,
            format(analysis_date, "%Y%m%d"), .rt_ver,
            generation_time_params$mean * 10, generation_time_params$sd * 10,
            .gtmax,
            .asc * 100, .trunc_k,
            .num1(.rpri$mean, 2.0) * 100, .num1(.rpri$sd, 0.5) * 100,
            .data_k)
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
    # forecast band — as the "Primary R(t) estimate". 44_reff_epinow2_check.R:118-125
    # documents and guards against exactly this; the diagnostic path did not.
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
#' usable, otherwise an onset imputed as sample date minus a drawn onset->sample delay, clamped
#' to the outbreak-week floor. This is the only dating consistent with EpiNow2's right-truncation
#' model (.epinow2_rt_fit()), which treats every case as onset-dated and still subject to the
#' onset->sample delay. Dating onset-less cases at their SAMPLE date instead (the former
#' convention) presents already-complete cases to that model as incomplete, so the recent tail
#' was inflated a second time and R biased upward (7 Sep 2026 snapshot, weekly R for 1-7 Sep:
#' 1.50 with sample-date dating vs 0.82 with this dating). Records without a date_index (the
#' complete_case onset mode) are excluded, exactly as in the weekly counts.
#'
#' Shared by the R(t) diagnostic (run_all.R), the cascade's EpiNow2 check
#' (44_reff_epinow2_check.R) and the Bayesian suite's R draws (bayes_rt_week_draws()).
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
# TWO REGIMES (see 00_config.R): "extract" for the deployed fit, "asof" for an LFO fold's fit.
# Both build a list of (date, confirm) frames, oldest first, as EpiNow2::estimate_truncation()
# requires, and each frame is the daily confirmed series AS IT STOOD at that vintage.
#
# SERIES BASIS: RECORDED onsets only, on both regimes, even though EpiNow2 is fitted to the full
# date_index series (recorded onset where usable, imputed otherwise). Reconstructing a vintage's
# date_index is not possible faithfully — the imputed delay is drawn positionally for all rows in
# one seeded vectorised call, row order does not replay across vintages, and only one archive
# folder carries its own delay parameters, so a June vintage would be re-imputed with a September
# fit. The recorded-onset basis has none of those dependencies. It is an approximation, and
# .trunc_basis_guard() below measures the approximation error on every run rather than assuming it.

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

#' Daily confirmed-by-recorded-onset series from ONE archive vintage.
#' Every column is read as character and only two are REQUIRED (the classification and the
#' recorded onset). That is deliberate: the archive spans 13 header signatures (26-36 columns)
#' and six early files lack lab_analysis_date, so naming a fixed column set would drop vintages.
#' Do not "restore" a column list here.
.trunc_frame_vintage <- function(dir, vintage, start) {
  f <- file.path(dir, sprintf("LINELIST_%s", format(vintage, "%d%m%Y")),
                 "dhis2_processed_linelist.csv")
  if (!file.exists(f)) return(NULL)
  x <- tryCatch(suppressWarnings(readr::read_csv(f, col_types = readr::cols(.default = "c"),
                                                 show_col_types = FALSE)),
                error = function(e) NULL)
  if (is.null(x)) return(NULL)
  need <- c("final_mve_case_classification", "date_of_symptom_onset")
  if (!all(need %in% names(x))) return(NULL)
  ons <- suppressWarnings(as.Date(x$date_of_symptom_onset))
  keep <- (x$final_mve_case_classification %in% CONFIRMED_STATUS) &
          !is.na(ons) & ons >= start & ons <= vintage
  .trunc_tabulate(ons[keep], start, vintage)
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

#' The as-of completeness curve, measured by ONSET COHORT, then fitted parametrically.
#'
#' WHY NOT estimate_truncation(). That deconvolves a lag from SNAPSHOTS, which is the only
#' option when the lag is recorded nowhere — true for the extract regime, where "appearance in
#' the extract" exists only as the difference between two archive vintages. Here the lag is
#' onset -> linelist_observation_date() and BOTH are columns on every record (the observation
#' date IS the sample date for 100% of confirmed records), so nothing needs deconvolving.
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
#' WHAT THIS DOES INSTEAD. For a lag L, every record whose onset is at least L days old has a
#' FULLY DETERMINED answer to "was it observed within L days?" — no truncation, no censoring,
#' no extrapolation. So estimate the completeness curve cohort-wise,
#'     F_hat(L) = P(delay <= L | onset <= end - L),
#' each lag using the largest unbiased sample available to it, and fit a parametric family to
#' that curve by least squares weighted by cohort size. Matching F(L + 0.5) to
#' P(delay <= L) is the same daily-rounding convention compute_truncation_weights() evaluates
#' under, so the fitted spec reproduces the measured completeness where it is used.
#'
#' GAMMA, NOT LOGNORMAL: fitted to the same curve, gamma tracks it to a maximum absolute error
#' of 0.040 across lags 0-45 against 0.083 for lognormal, and its anchor-week weight is within
#' 1.3% of the measured one against 9% for lognormal.
#'
#' NO ROLLING WINDOW, unlike the extract panel. Windowing would starve the long lags (a 56-day
#' window leaves 11 days of onsets at lag 45). It is not needed for recency either: the
#' epidemic grew, so the short-lag cohorts that determine the anchor week are already dominated
#' by recent cases -- the lag-6 cohort holds 5,005 records against the 2,882 that are old
#' enough to be fully observed.
#'
#' CONFIRMED ONLY, deliberately: the folds correct the confirmed series.
#' effective_onset_sample_delay() pooled over every classification (~14k pairs against ~5.0k
#' confirmed) and suspected cases are sampled on an alert-driven schedule, so it returns
#' 7.67 d where the CONFIRMED stratum gives 10.14 d. That is why the fold nowcast
#' under-corrected even though its LEG was already right. (Both are the EpiDist marginal; the
#' interval-censored MLE puts the same contrast at 7.67 vs 9.11 d.)
#'
#' @return an EpiNow2 dist_spec (Gamma), or NULL.
.trunc_fit_asof_cohort <- function(ll, end, start) {
  if (!exists("linelist_observation_date", mode = "function")) return(NULL)
  obs <- tryCatch(linelist_observation_date(ll, as.Date(end), caller = "trunc-asof"),
                  error = function(e) NULL)
  if (is.null(obs)) return(NULL)
  ons <- suppressWarnings(as.Date(ll$date_of_symptom_onset))
  end <- as.Date(end)
  keep <- (ll$confirmed %in% TRUE) & !is.na(obs) & !is.na(ons) & as.Date(obs) <= end
  o  <- ons[keep]
  dl <- as.numeric(as.Date(obs[keep]) - o)
  # dl < 0 is a data error (observation before onset), not a short delay. Dropped rather than
  # clamped to 0, which would pile spurious mass on the shortest lag — exactly where the
  # anchor-week weight is most sensitive.
  ok <- is.finite(dl) & dl >= 0; o <- o[ok]; dl <- dl[ok]
  mx <- as.integer(get0("TRUNC_MAX_DAYS", ifnotfound = 45L))
  lags <- 0:mx
  n_l <- vapply(lags, function(L) sum(o <= end - L), integer(1))
  f_l <- vapply(lags, function(L) {
    c_ <- o <= end - L
    if (!any(c_)) NA_real_ else mean(dl[c_] <= L)
  }, numeric(1))
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
  w <- as.numeric(n_l[usable]); fu <- f_l[usable]; lu <- lags[usable]
  sse <- function(par) {
    sh <- exp(par[1]); rt <- exp(par[2])
    if (!is.finite(sh) || !is.finite(rt) || sh <= 0 || rt <= 0) return(1e10)
    cdf <- pgamma(lu + 0.5, shape = sh, rate = rt)
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
  err <- max(abs(pgamma(lu + 0.5, shape = sh, rate = rt) - fu))
  message(sprintf(paste0("[trunc] asof: cohort completeness on %d confirmed records ",
                         "(lag-0 cohort %d, lag-%d cohort %d); gamma fit max |F - F_hat| = %.4f"),
                  length(dl), n_l[1], mx, n_l[mx + 1L], err))
  # A curve the family cannot follow is a wrong family, not a noisy fit — refuse rather than
  # correct every week with it. 0.08 is twice the realised gamma error and below the lognormal
  # error the family comparison rejected.
  if (!is.finite(err) || err > 0.08) {
    warning(sprintf("[trunc] asof: gamma cannot follow the measured completeness curve (max |F - F_hat| = %.3f); REFUSED.",
                    err), call. = FALSE)
    return(NULL)
  }
  tryCatch(EpiNow2::Gamma(shape = sh, rate = rt, max = mx), error = function(e) NULL)
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

#' Measure the recorded-onset vs date_index approximation, per run.
#' Compares the completeness curve of the panel's own basis against one dated by the LIVE
#' date_index but censored by MEMBERSHIP ONLY (which records were confirmed-and-visible at each
#' vintage). No re-imputation, no RNG — so it measures the basis difference and nothing else.
.trunc_basis_guard <- function(ll, end, start, regime, tol = 0.05) {
  ons <- suppressWarnings(as.Date(ll$date_of_symptom_onset))
  di  <- suppressWarnings(as.Date(ll$date_index))
  obs <- tryCatch(linelist_observation_date(ll, as.Date(end), caller = "trunc-guard"),
                  error = function(e) NULL)
  if (is.null(obs)) return(invisible(NA_real_))
  # THE TWO BASES ARE TAKEN OVER DIFFERENT RECORD SETS, and that is the entire point.
  # The first version required !is.na(ons) AND !is.na(di) and computed both curves on that one
  # subset. But date_index EQUALS date_of_symptom_onset wherever a usable recorded onset
  # exists, so that subset is precisely where the two bases agree by definition: it excluded
  # the ~24% of confirmed cases whose date_index is IMPUTED, which are the only records that
  # create a basis difference at all. The guard reported 0.0082 and could never have fired.
  # Compared honestly it is 0.0355 (5,106 recorded-onset records vs 6,669 date_index ones).
  #
  # RESTRICT TO THE RECORDS THIS REGIME ACTUALLY SEES (observable by `end`), so the number
  # describes the panel that was just built rather than the whole line list.
  base  <- (ll$confirmed %in% TRUE) & !is.na(obs) & as.Date(obs) <= as.Date(end)
  k_rec <- base & !is.na(ons) & ons >= start
  k_idx <- base & !is.na(di)  & di  >= start
  if (sum(k_rec) < 100L || sum(k_idx) < 100L) return(invisible(NA_real_))
  lag_rec <- as.numeric(as.Date(obs[k_rec]) - ons[k_rec])
  lag_idx <- as.numeric(as.Date(obs[k_idx]) - di[k_idx])
  ks <- 0:get0("TRUNC_MAX_DAYS", ifnotfound = 45L)
  z_rec <- vapply(ks, function(k) mean(lag_rec >= 0 & lag_rec <= k), numeric(1))
  z_idx <- vapply(ks, function(k) mean(lag_idx >= 0 & lag_idx <= k), numeric(1))
  d <- max(abs(z_rec - z_idx), na.rm = TRUE)
  # DIRECTION MATTERS MORE THAN MAGNITUDE. The difference is one-signed on this snapshot
  # (z_idx >= z_rec at every lag): the imputed-onset records are observed sooner after their
  # date_index than the recorded-onset records are after their onset, so the series EpiNow2 is
  # fitted to is LESS truncated than the panel the truncation is fitted on. A truncation
  # fitted on the recorded-onset basis is therefore slightly too long, and both consumers
  # over-correct the recent weeks rather than under-correcting them. Reported, not silently
  # absorbed: it cannot be removed without a date_index that replays per vintage, which the
  # archive does not support (see the SERIES BASIS note at the top of this section).
  .dir <- if (all(z_idx >= z_rec - 1e-12)) " (one-signed: fitted truncation is too LONG, so both consumers over-correct)"
          else if (all(z_rec >= z_idx - 1e-12)) " (one-signed: fitted truncation is too SHORT, so both consumers under-correct)"
          else ""
  if (is.finite(d) && d > tol)
    warning(sprintf(paste0("[trunc] %s: the recorded-onset basis used for the truncation fit ",
                           "differs from the date_index basis EpiNow2 is fitted to by %.3f at ",
                           "worst (tolerance %.2f)%s. The fitted truncation may not describe the ",
                           "fitted series."), regime, d, tol, .dir), call. = FALSE)
  else
    message(sprintf("[trunc] %s: recorded-onset vs date_index basis agree to %.4f (tol %.2f)%s",
                    regime, d, tol, .dir))
  invisible(d)
}

#' Fingerprint a panel so a cached truncation fit cannot be served for a different one.
.trunc_fingerprint <- function(panel, regime) {
  parts <- vapply(panel, function(z)
    sprintf("%s:%s:%.0f", format(min(z$date)), format(max(z$date)), sum(z$confirm)), character(1))
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
      return(d)
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
  set.seed(get0("RANDOM_SEED", ifnotfound = 20260704L))
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

#' Convert a fitted truncation dist_spec into the pipeline's own delay-spec shape.
#'
#' The NOWCAST (04_nowcasting.R) evaluates a delay through delay_cdf() (00_config.R), which
#' takes a list(family, params, rate, mean, sd) — not an EpiNow2 dist_spec. Both consumers must
#' describe the SAME truncation, so this converts rather than re-deriving.
#'
#' The uncertain parameters are collapsed to their means (fix_parameters): the nowcast applies a
#' deterministic weight per week, so it has nowhere to put the uncertainty. That is a real
#' difference between the two consumers, stated here rather than hidden — the R(t) fit keeps the
#' posterior, the nowcast uses the point estimate.
#'
#' @return list in delay_cdf() shape, or NULL if the distribution cannot be represented.
.trunc_as_delay_spec <- function(dist, regime = "unknown") {
  if (is.null(dist)) return(NULL)
  fx <- tryCatch(EpiNow2::fix_parameters(dist), error = function(e) NULL)
  if (is.null(fx)) return(NULL)
  pr <- tryCatch(EpiNow2::get_parameters(fx), error = function(e) NULL)
  if (is.null(pr) || !length(pr)) return(NULL)
  nm <- names(pr)
  # Only the families delay_cdf() can evaluate. estimate_truncation() is fitted with a
  # LogNormal here, so `lnorm` is the live path; gamma is accepted in case the prior family
  # is ever changed. Anything else returns NULL rather than silently degrading to pexp().
  if (all(c("meanlog", "sdlog") %in% nm)) {
    ml <- as.numeric(pr$meanlog); sl <- as.numeric(pr$sdlog)
    if (!is.finite(ml) || !is.finite(sl) || sl <= 0) return(NULL)
    mean <- exp(ml + sl^2 / 2)
    sd   <- sqrt((exp(sl^2) - 1) * exp(2 * ml + sl^2))
    out <- list(family = "lnorm", params = c(meanlog = ml, sdlog = sl),
                mean = mean, sd = sd)
  } else if (all(c("shape", "rate") %in% nm)) {
    sh <- as.numeric(pr$shape); rt <- as.numeric(pr$rate)
    if (!is.finite(sh) || !is.finite(rt) || sh <= 0 || rt <= 0) return(NULL)
    out <- list(family = "gamma", params = c(shape = sh, rate = rt),
                mean = sh / rt, sd = sqrt(sh) / rt)
  } else return(NULL)
  # delay_cdf() falls back to pexp(q, dp$rate) when the named parameters are unusable, and
  # compute_truncation_weights() VALIDATES that `rate` is finite and positive, so it must be
  # present and sane even though the lnorm/gamma branch never reads it. 1/mean is the
  # moment-matched exponential rate — the correct fallback if it is ever reached.
  out$rate <- 1 / out$mean
  out$estimator <- sprintf("estimate_truncation_%s", regime)
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
                   "delay describes the wrong reporting leg (onset->sample ~7.7 d against a ",
                   "measured onset->appearance ~11.4 d) and biases R(t) DOWN at the present. ",
                   "For comparison only."), call. = FALSE)
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
  # PANEL END DATE — and the ASOF REGIME IS SHARED ACROSS FOLDS, deliberately.
  #
  # The extract panel ends at the issue date (or the analysis date when none is given).
  #
  # The AS-OF panel always ends at the ANALYSIS DATE, whatever issue_date says. The truncation
  # is a SHARED NUISANCE PARAMETER of the LFO — the same choice run_invasion_lfo() already
  # makes for its nowcast delay, and for the same reason: a nuisance parameter refitted per
  # fold is correlated with fold index, so apparent skill drifts with calendar time for a
  # purely artefactual reason. Ending every fold's panel at the analysis date makes the fitted
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
  # is one more reason fold skill is optimistic relative to deployment. See the "RESIDUAL, AND
  # IT IS ACCEPTED" note in 04_nowcasting.R.
  end <- if (deployed && !is.null(issue_date)) as.Date(issue_date) else as.Date(analysis_date)

  # TWO REGIMES, TWO ESTIMATORS — because the lag is observable in one and not the other.
  #   extract : the appearance lag exists only as the difference between archive vintages,
  #             so it is DECONVOLVED from snapshots (estimate_truncation(), Stan, cached).
  #   asof    : the lag is onset -> linelist_observation_date() and both are columns, so the
  #             completeness curve is MEASURED by onset cohort and a gamma fitted to it (no
  #             Stan, no cache — it is deterministic and sub-second).
  if (deployed) {
    panel <- .trunc_validate_panel(.trunc_panel_extract(end, start), regime)
    if (is.null(panel)) return(NULL)
    if (!is.null(ll)) .trunc_basis_guard(ll, end, start, regime)
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
  # The distribution is ESTIMATED from line-list vintages (section 4b above), not derived from
  # the onset->sample delay. That delay describes the wrong leg: this series is truncated by
  # onset -> APPEARANCE IN THE EXTRACT, which adds an unrecorded data-entry lag (measured
  # sample->appearance median 2.0 d, mean 2.6 d). Understating it made EpiNow2 treat the recent
  # weeks as near-complete, read the reporting dip as a decline, and return an R(t) biased DOWN
  # at the present — the value RT_WINDOW_WEEKS = 1 hands to the 2-week forecast and the cascade.
  #
  # `truncation` is supplied by the CALLER, which knows which series it is fitting:
  #   estimate_rt_epinow2()  -> deployed fit  -> extract regime  (onset -> appearance)
  #   bayes_rt_week_draws()  -> one LFO fold  -> asof regime     (onset -> sample)
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

# --- 6e. Summary printout ---
message("\n=== Parameter summary ===")
message(sprintf("GT profiles computed: %s", paste(names(GT_PMFS), collapse = ", ")))
message("=== 02_epi_params.R complete ===\n")
}
