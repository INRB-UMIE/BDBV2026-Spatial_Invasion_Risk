# =============================================================================
# 07_hhh4_model.R — hhh4 Endemic-Epidemic Spatiotemporal Models (C1a–C1d)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Model (AS FITTED): mu_{i,t} = nu_it + lambda_it * Y_{i,t-1} + phi_it * sum_j w_{ji} * Y_{j,t-1}
# NB: surveillance::hhh4 here uses its DEFAULT single-week (lag-1) autoregression, NOT the
# generation-time-weighted history Y*_{i,t} = sum_tau g_w[tau] Y_{i,t-tau}. compute_gt_weighted_lag()
# below implements Y* but is NOT currently wired into the fit, and `gt_pmf_daily` is accepted for
# API symmetry only (it does not enter the model). To make hhh4 GT-sensitive, feed Y* via funct_lag.
#
# Variants:
#   C1a — power-law distance weights
#   C1b — Flowminder mobility weights
#   C1c — time-varying phi (epidemic phase interaction)
#   C1d — covariate-modulated lambda
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  # Every surveillance call here is namespaced (surveillance::sts/hhh4/weightedSumNE), so the
  # package only needs to be INSTALLED, not attached. Guard instead of library() so a missing
  # 'surveillance' does not crash the whole pipeline at source() time — 07 is sourced
  # unconditionally, even in a Bayesian-only run that never fits hhh4. The fit/simulate wrappers
  # are tryCatch-guarded and degrade to NULL when the package (or a fit) is unavailable.
  if (!requireNamespace("surveillance", quietly = TRUE))
    warning("[hhh4] Package 'surveillance' not installed; hhh4 (C1a-C1d) comparators will be skipped.")
})

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)[1]) a else b

# ---------------------------------------------------------------------------
# Helper: construct STS (surveillance time series) object
# ---------------------------------------------------------------------------

#' Build an STS object from a zone-week matrix
#'
#' @param zone_week_nc tibble with health_zone, week_start, confirmed_nc
#' @param zones_all    canonical zone name vector (defines rows)
#' @param outbreak_start Date of first week to include
#' @param n_ahead      number of future (placeholder) weeks to append after the
#'                     observed period. These carry zero observations and
#'                     carried-forward population, so that simulate.hhh4() can
#'                     forecast onto valid in-range time indices.
#' @return surveillance::sts object
make_sts <- function(zone_week_nc, zones_all, pop_vec = NULL,
                     outbreak_start = OUTBREAK_START, n_ahead = 0L) {
  weeks <- sort(unique(zone_week_nc$week_start))
  weeks <- weeks[weeks >= outbreak_start]
  n_zones <- length(zones_all)
  n_obs   <- length(weeks)

  # Append n_ahead future week dates (weekly cadence) as forecast placeholders
  if (n_ahead > 0L) {
    future_weeks <- max(weeks) + 7 * seq_len(n_ahead)
    weeks <- c(weeks, future_weeks)
  }
  n_weeks <- length(weeks)

  # Build observed count matrix (zones × weeks), integer-rounded.
  # Future placeholder columns remain zero (overwritten by simulation).
  Y_mat <- matrix(0L, nrow = n_zones, ncol = n_weeks,
                  dimnames = list(zones_all, as.character(weeks)))
  for (k in seq_len(n_obs)) {
    sub <- zone_week_nc |>
      filter(week_start == weeks[k], health_zone %in% zones_all)
    vals <- ifelse(is.na(sub$confirmed_nc), 0, sub$confirmed_nc)
    Y_mat[sub$health_zone, k] <- pmax(round(vals), 0L)
  }

  # Build population matrix (n_weeks × n_zones) for endemic offset.
  # Population is treated as static, so future rows carry the same values.
  pop_mat <- if (!is.null(pop_vec)) {
    pop_aligned <- pop_vec[zones_all]
    pop_aligned[is.na(pop_aligned)] <- median(pop_vec, na.rm = TRUE)
    pop_aligned <- pmax(pop_aligned, 1)
    matrix(rep(pop_aligned, each = n_weeks), nrow = n_weeks, ncol = n_zones,
           dimnames = list(NULL, zones_all))
  } else NULL

  # surveillance::sts expects:
  #   observed: n_weeks × n_zones matrix (time in rows, units in columns)
  sts_obj <- surveillance::sts(
    observed   = t(Y_mat),       # transpose: time × zones
    start      = c(as.integer(format(min(weeks), "%G")),   # ISO year (pairs with %V)
                   as.integer(format(min(weeks), "%V"))),  # ISO week
    frequency  = 52,
    population = pop_mat
  )
  sts_obj
}

# ---------------------------------------------------------------------------
# Helper: GT-weighted lag matrix (Y*)
# ---------------------------------------------------------------------------

#' Compute Y*_{i,t} = sum_tau g_w[tau] * Y_{i,t-tau} for all zones and times
#'
#' @note Currently UNUSED — the fitted hhh4 models use surveillance's default lag-1
#'   autoregression (see file header). Retained for a future GT-weighted variant.
#' @param Y_mat    n_weeks × n_zones observed matrix (as from observed(sts))
#' @param G_weekly weekly GT PMF vector
#' @return n_weeks × n_zones matrix of Y* values
compute_gt_weighted_lag <- function(Y_mat, G_weekly) {
  n_weeks <- nrow(Y_mat)
  n_zones <- ncol(Y_mat)
  K <- length(G_weekly)
  Y_star <- matrix(0, n_weeks, n_zones, dimnames = dimnames(Y_mat))
  for (t in seq_len(n_weeks)) {
    for (k in seq_len(K)) {
      t_past <- t - k
      if (t_past < 1) break
      Y_star[t, ] <- Y_star[t, ] + G_weekly[k] * Y_mat[t_past, ]
    }
  }
  Y_star
}

# ---------------------------------------------------------------------------
# Fit hhh4 variants
# ---------------------------------------------------------------------------

#' Fit a hhh4 endemic-epidemic model
#'
#' @param sts_obj    sts object from make_sts()
#' @param W_ne       n_zones × n_zones neighbourhood weight matrix indexed
#'                   [source j, destination i]: W_ne[j,i] = fraction of source j's
#'                   outflow that reaches i. surveillance forms the ne term as
#'                   `observed %*% weights` (weightedSumNE), i.e. unit i's pressure is
#'                   sum_j Y_j * weights[j,i], so this is W_outflow ITSELF, NOT its
#'                   transpose (passing t(W_outflow) would reverse the flow direction).
#' @param pop_vec    named population vector (for endemic offset)
#' @param covariates tibble with nom, log_pop, log_hs, ccvi, positivity (optional)
#' @param variant    "C1a" | "C1b" | "C1c" | "C1d"
#' @return list(fit, aic, bic, loglik, coef_tbl, variant)
fit_hhh4 <- function(sts_obj, W_ne, pop_vec, covariates = NULL,
                     variant = "C1b", n_fit = NULL) {
  n_zones <- ncol(observed(sts_obj))
  n_weeks <- nrow(observed(sts_obj))
  zone_names <- colnames(observed(sts_obj))

  # Number of weeks to fit on (training period). Any trailing weeks beyond n_fit
  # are forecast placeholders and must be excluded from estimation.
  if (is.null(n_fit)) n_fit <- n_weeks
  n_fit <- min(n_fit, n_weeks)

  # Population offset aligned to zone order
  pop_aligned <- pop_vec[zone_names]
  pop_aligned[is.na(pop_aligned)] <- median(pop_vec, na.rm = TRUE)
  pop_aligned <- pmax(pop_aligned, 1)

  # W_ne is expected already row-normalised (normalised by run_all_hhh4 before passing)
  W_ne_norm <- W_ne[zone_names, zone_names]

  control <- switch(variant,
    C1a = {
      # Power-law distance weights (estimated via W.powerlaw)
      # Uses surveillance's built-in power-law neighbour weight function
      # Requires a neighbourhood matrix with integer distances
      # Approximation: derive from OSRM distances via binned distance classes
      # Fallback to W_ne if not available
      list(
        ar  = list(f = ~ 1),
        ne  = list(f = ~ 1, weights = W_ne_norm, normalize = FALSE),
        end = list(f = ~ 1, offset = population(sts_obj)),
        family = "NegBin1",
        data = list(pop = pop_aligned)
      )
    },
    C1b = {
      # Flowminder-based weights
      list(
        ar  = list(f = ~ 1),
        ne  = list(f = ~ 1, weights = W_ne_norm, normalize = FALSE),
        end = list(f = ~ 1, offset = population(sts_obj)),
        family = "NegBin1"
      )
    },
    C1c = {
      # Time-varying phi: endemic baseline + week-index trend
      list(
        ar  = list(f = ~ 1),
        ne  = list(f = ~ 1 + t,  # linear week trend on log(phi)
                   weights = W_ne_norm, normalize = FALSE),
        end = list(f = ~ 1, offset = population(sts_obj)),
        family = "NegBin1",
        data = list(t = seq_len(n_weeks) / n_weeks)  # scaled [0,1]
      )
    },
    C1d = {
      # Covariate-modulated lambda via endemic + zone-level covariates
      if (!is.null(covariates)) {
        cov_aligned <- covariates[match(zone_names, covariates$nom), ]
        log_hs  <- coalesce(cov_aligned$log_healthsite_count, 0)
        ccvi_v  <- coalesce(cov_aligned$ccvi, 0)
        # hhh4 requires covariates as nTime x nUnit matrices (or length-nTime
        # vectors), NOT length-nUnit vectors — a per-zone vector throws
        # "Covariate '...' is not suitably specified" and the tryCatch below then
        # silently drops the whole C1d variant. These covariates are static per
        # zone, so tile each zone's value across all weeks (constant down columns).
        log_hs_m <- matrix(rep(log_hs, each = n_weeks), nrow = n_weeks, ncol = n_zones)
        ccvi_m   <- matrix(rep(ccvi_v, each = n_weeks), nrow = n_weeks, ncol = n_zones)
        list(
          ar  = list(f = ~ 1 + log_hs + ccvi),
          ne  = list(f = ~ 1, weights = W_ne_norm, normalize = FALSE),
          end = list(f = ~ 1, offset = population(sts_obj)),
          family = "NegBin1",
          data = list(log_hs = log_hs_m, ccvi = ccvi_m)
        )
      } else {
        message("[hhh4] C1d requires covariates; falling back to C1b")
        list(
          ar  = list(f = ~ 1),
          ne  = list(f = ~ 1, weights = W_ne_norm, normalize = FALSE),
          end = list(f = ~ 1, offset = population(sts_obj)),
          family = "NegBin1"
        )
      }
    },
    stop("Unknown hhh4 variant: ", variant)
  )

  # Restrict estimation to the training period (weeks 2..n_fit); hhh4 always
  # drops t=1 because the autoregressive term needs a lag. Trailing forecast
  # placeholder weeks (n_fit+1 .. n_weeks) are excluded from the likelihood.
  control$subset <- 2:n_fit

  message(sprintf("[hhh4] Fitting variant %s (%d zones, fit weeks 2:%d of %d)...",
                  variant, n_zones, n_fit, n_weeks))

  fit <- tryCatch(
    surveillance::hhh4(sts_obj, control = control),
    error = function(e) {
      warning(sprintf("[hhh4] %s fit failed: %s", variant, e$message))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  # Coefficient table with standard errors; se=TRUE reads fit$cov, which is
  # absent/non-numeric when the optimizer did not converge — fall back to point
  # estimates rather than erroring out of the whole run.
  converged <- isTRUE(fit$convergence)
  coef_tbl <- tryCatch(
    as_tibble(t(coef(fit, se = TRUE)), .name_repair = "unique"),
    error = function(e) tryCatch(
      as_tibble(t(coef(fit, se = FALSE)), .name_repair = "unique"),
      error = function(e2) tibble::tibble()
    )
  )
  if (!converged) {
    warning(sprintf("[hhh4] %s did not converge; forecasts may be unreliable.",
                    variant))
  }

  list(
    fit       = fit,
    converged = converged,
    aic       = suppressWarnings(AIC(fit)),
    bic       = suppressWarnings(BIC(fit)),
    loglik    = suppressWarnings(logLik(fit)),
    coef_tbl  = coef_tbl,
    variant   = variant
  )
}

# ---------------------------------------------------------------------------
# Forecast from a fitted hhh4 object
# ---------------------------------------------------------------------------

#' Generate stochastic forecasts from a fitted hhh4 model
#'
#' @param fit_obj  list from fit_hhh4()
#' @param sts_obj  the sts object used for fitting (extended with n_ahead
#'                 forecast-placeholder weeks by make_sts)
#' @param n_ahead  number of weeks to forecast ahead (1 or 2)
#' @param n_fit    number of training weeks (the observed period); the future
#'                 indices simulated are (n_fit+1)..(n_fit+n_ahead)
#' @param n_sim    number of stochastic simulations
#' @param seed     random seed
#' @return tibble: health_zone, horizon, mu_forecast, p_invasion, q05..q95
forecast_hhh4 <- function(fit_obj, sts_obj, n_ahead = 2, n_fit = NULL,
                          n_sim = N_SIMULATIONS, seed = RANDOM_SEED) {
  if (is.null(fit_obj)) return(NULL)
  fit       <- fit_obj$fit
  zone_names <- colnames(observed(sts_obj))
  n_zones   <- length(zone_names)

  set.seed(seed)

  # The sts object spans n_fit training weeks + n_ahead placeholder weeks.
  # Simulate the placeholder indices, seeding from the last training week.
  n_total <- nrow(observed(sts_obj))
  if (is.null(n_fit)) n_fit <- n_total - n_ahead
  n_train <- n_fit
  # `simulate` is the stats generic; surveillance registers simulate.hhh4 as an
  # S3 method but does NOT export a `simulate` symbol, so surveillance::simulate
  # errors. Fetch the registered method directly for robust dispatch.
  simulate_hhh4 <- getS3method("simulate", "hhh4")
  sims <- tryCatch({
    simulate_hhh4(fit,
                  nsim    = n_sim,
                  seed    = seed,
                  y.start = observed(sts_obj)[n_train, , drop = FALSE],
                  subset  = (n_train + 1):(n_train + n_ahead))
  }, error = function(e) {
    warning("[hhh4] simulate() failed: ", e$message)
    NULL
  })

  if (is.null(sims)) return(NULL)

  # Parse simulation output
  # If sims is a list of sts objects:
  if (is.list(sims) && inherits(sims[[1]], "sts")) {
    # Extract final n_ahead weeks from each simulation
    sim_arr <- array(NA_real_, dim = c(n_ahead, n_zones, n_sim))
    for (s in seq_len(n_sim)) {
      obs_s <- observed(sims[[s]])  # n_ahead × n_zones
      if (nrow(obs_s) >= n_ahead) {
        sim_arr[, , s] <- obs_s[seq_len(n_ahead), , drop = FALSE]
      }
    }
  } else if (is.array(sims) && length(dim(sims)) == 3) {
    # Array form [time, unit, sim]
    sim_arr <- sims
  } else {
    warning("[hhh4] Unexpected simulate() output format")
    return(NULL)
  }

  results <- vector("list", n_ahead)
  for (h in seq_len(n_ahead)) {
    # CUMULATIVE counts over weeks 1..h — the LFO outcome is "first confirmed case in
    # weeks (cutoff, cutoff+h]" (cumulative, counted once), and the renewal (mu_cum) and
    # SEIR (summed arrivals) models score against it with a CUMULATIVE probability. Using
    # week-h-only (marginal) counts understated hhh4's invasion probability at h>=2 (a
    # zone with a case in week 1 but none in week 2 is an invasion by h=2). Summing over
    # the first array dim also guards the n_sim==1 vector-drop.
    Y_cum <- apply(sim_arr[seq_len(h), , , drop = FALSE], c(2, 3), sum)  # n_zones × n_sim
    if (!is.matrix(Y_cum)) Y_cum <- matrix(Y_cum, nrow = n_zones)

    mu_h      <- rowMeans(Y_cum, na.rm = TRUE)
    p_inv     <- rowMeans(Y_cum > 0, na.rm = TRUE)

    results[[h]] <- tibble(
      health_zone = zone_names,
      horizon     = h,
      mu_forecast = mu_h,
      p_invasion  = p_inv,
      q05  = apply(Y_cum, 1, quantile, 0.05,  na.rm = TRUE),
      q20  = apply(Y_cum, 1, quantile, 0.20,  na.rm = TRUE),
      q25  = apply(Y_cum, 1, quantile, 0.25,  na.rm = TRUE),
      q75  = apply(Y_cum, 1, quantile, 0.75,  na.rm = TRUE),
      q80  = apply(Y_cum, 1, quantile, 0.80,  na.rm = TRUE),
      q95  = apply(Y_cum, 1, quantile, 0.95,  na.rm = TRUE),
      method = paste0("C1_", fit_obj$variant)
    )
  }
  bind_rows(results)
}

# ---------------------------------------------------------------------------
# Run all hhh4 variants
# ---------------------------------------------------------------------------

#' Run hhh4 C1a–C1d with a specified mobility matrix
#'
#' @param zone_week_nc  tibble with health_zone, week_start, confirmed_nc
#' @param W_outflow     519×519 outflow mobility matrix (row = origin, col = dest)
#' @param gt_pmf_daily  daily GT PMF
#' @param pop_vec       named population vector
#' @param covariates    tibble with zone-level covariates
#' @param zones_all     canonical zone name vector
#' @param t_idx         week index for training cutoff
#' @param horizons      forecast horizons (default c(1L, 2L))
#' @param n_sim         number of stochastic draws
#' @param mobility_id   label for the mobility matrix (for output method column)
#' @param gt_profile    label for the GT profile (for output)
#' @param cache         logical; cache fitted models to disk
#' @return combined tibble of hhh4 forecasts across all variants
run_all_hhh4 <- function(zone_week_nc, W_outflow, gt_pmf_daily,
                         pop_vec, covariates = NULL, zones_all,
                         t_idx, horizons = LFO_HORIZONS,
                         n_sim = N_SIMULATIONS, seed = RANDOM_SEED,
                         mobility_id = MOBILITY_PRIMARY,
                         gt_profile = GT_PRIMARY,
                         training_cutoff = NULL,
                         cache = TRUE) {

  # Zones present in both the linelist and the mobility matrix
  zones_fit <- zones_all[zones_all %in% rownames(W_outflow)]
  if (length(zones_fit) < 5) stop("[hhh4] Too few zones overlap mobility matrix")

  # Subset training data to t_idx weeks
  all_weeks <- sort(unique(zone_week_nc$week_start))
  train_weeks <- all_weeks[seq_len(min(t_idx, length(all_weeks)))]
  zone_week_train <- zone_week_nc |> filter(week_start %in% train_weeks)

  # Build STS object extended with max(horizons) forecast-placeholder weeks, so
  # simulate.hhh4() can forecast onto valid in-range indices. n_fit is the count
  # of training weeks actually present at/after OUTBREAK_START (make_sts filters).
  n_ahead <- max(horizons)
  sts_obj <- make_sts(zone_week_train, zones_fit, pop_vec = pop_vec,
                      n_ahead = n_ahead)
  n_fit <- nrow(observed(sts_obj)) - n_ahead

  # hhh4 needs at least 3 in-outbreak training weeks (t=1 is dropped for the AR
  # lag, and >=2 remaining points to estimate). Pre-outbreak fold cutoffs leave
  # too few weeks after make_sts's OUTBREAK_START filter → skip gracefully.
  if (n_fit < 3L) {
    message(sprintf(
      "[hhh4] Only %d in-outbreak training week(s) (need >=3); skipping hhh4 for this cutoff.",
      max(n_fit, 0L)))
    return(tibble())
  }

  # Neighbourhood weights for hhh4's epidemic component. surveillance::weightedSumNE
  # forms unit i's neighbour pressure as sum_j weights[j,i] * Y_{j,t-1}, i.e. the
  # weight is indexed [source j, destination i]. The flow from source j INTO i is
  # exactly W_outflow[j,i] (W_outflow is [origin, destination], row-stochastic), so
  # hhh4 must receive W_outflow ITSELF — NOT its transpose (which would weight j's
  # cases by how much i travels to j, reversing the direction). Row-normalise so
  # each source's outward influence sums to 1 (a no-op guard when already stochastic).
  W_ne <- W_outflow[zones_fit, zones_fit]
  row_sums <- rowSums(W_ne); row_sums[row_sums == 0] <- 1
  W_ne_norm <- W_ne / row_sums

  variants <- c("C1b", "C1c", "C1d")  # C1a approximated by C1b for now

  results <- vector("list", length(variants))
  for (v_idx in seq_along(variants)) {
    v <- variants[v_idx]

    # Check cache
    cache_path <- file.path(OUT_FORECASTS,
      sprintf("hhh4_fit_%s_%s_t%02d.rds", v, mobility_id, t_idx))
    if (cache && file.exists(cache_path)) {
      message("[hhh4] Loading cached fit: ", basename(cache_path))
      fit_obj <- readRDS(cache_path)
    } else {
      fit_obj <- fit_hhh4(sts_obj, W_ne_norm, pop_vec, covariates,
                          variant = v, n_fit = n_fit)
      if (cache && !is.null(fit_obj)) saveRDS(fit_obj, cache_path)
    }

    if (is.null(fit_obj)) next

    fc <- forecast_hhh4(fit_obj, sts_obj, n_ahead = n_ahead, n_fit = n_fit,
                        n_sim = n_sim, seed = seed)
    if (!is.null(fc)) {
      results[[v_idx]] <- fc |>
        filter(horizon %in% horizons) |>
        mutate(
          method          = paste0("C1_", v, "_", mobility_id),
          mobility_id     = mobility_id,
          gt_profile      = gt_profile,
          training_cutoff = training_cutoff %||% max(train_weeks)
        )
    }
  }

  bind_rows(Filter(Negate(is.null), results))
}
