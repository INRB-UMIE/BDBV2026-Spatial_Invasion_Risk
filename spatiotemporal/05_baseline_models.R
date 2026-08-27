# =============================================================================
# 05_baseline_models.R — Baseline Predictive Models (B1–B6)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Models: B1 distance-weighted average, B2 K-NN average, B3 spatial diffusion,
#         B4 gravity-decay invasion, B5 population null, B6 persistent-zero
#
# All return a tibble: health_zone, horizon, mu_forecast, p_invasion, method
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---------------------------------------------------------------------------
# Helper: build wide case matrix (zones × weeks)
# ---------------------------------------------------------------------------

#' Convert zone-week tibble to wide numeric matrix
#' @param zone_week tibble with health_zone, week_start, confirmed_nc columns
#' @param zones_all character vector of all canonical zone names (rows)
#' @return numeric matrix (n_zones × n_weeks), rownames = zones, colnames = week indices
zone_week_to_wide <- function(zone_week, zones_all) {
  weeks <- sort(unique(zone_week$week_start))
  n_zones <- length(zones_all)
  n_weeks <- length(weeks)

  mat <- matrix(0, nrow = n_zones, ncol = n_weeks,
                dimnames = list(zones_all, as.character(weeks)))

  for (k in seq_along(weeks)) {
    w <- weeks[k]
    sub <- zone_week[zone_week$week_start == w, c("health_zone", "confirmed_nc")]
    sub <- sub[sub$health_zone %in% zones_all, ]
    mat[sub$health_zone, k] <- ifelse(is.na(sub$confirmed_nc), 0, sub$confirmed_nc)
  }
  mat
}

# ---------------------------------------------------------------------------
# B1: Distance-weighted neighbour average (power-law decay)
# ---------------------------------------------------------------------------

#' B1: Forecast via inverse-distance-power-law weighted mean of neighbour cases
#'
#' @param Y_wide  wide case matrix (n_zones × n_weeks)
#' @param osrm_mat travel-time matrix (n_zones × n_zones) in minutes, rownames=colnames=zones
#' @param t_idx   integer week index (1-based) to use as current observation
#' @param horizons integer vector of forecast horizons in weeks (1, 2)
#' @param zones_all character vector of all zone names
#' @param alpha   power-law exponent (default 1.0)
#' @return tibble: health_zone, horizon, mu_forecast, p_invasion, method
forecast_B1 <- function(Y_wide, osrm_mat, t_idx, horizons, zones_all,
                        alpha = 1.0) {
  stopifnot(is.matrix(Y_wide), is.matrix(osrm_mat))
  stopifnot(all(rownames(Y_wide) %in% zones_all),
            all(rownames(osrm_mat) %in% zones_all))

  # Align matrices to same zone order
  ord  <- zones_all[zones_all %in% rownames(osrm_mat) & zones_all %in% rownames(Y_wide)]
  D    <- osrm_mat[ord, ord]
  Y_t  <- Y_wide[ord, min(t_idx, ncol(Y_wide))]

  results <- vector("list", length(horizons))
  Y_curr  <- Y_t

  for (h_idx in seq_along(horizons)) {
    h <- horizons[h_idx]
    # Carry forward the previous mu after the first horizon (index-based guard)
    if (h_idx > 1) {
      Y_curr <- mu_prev
    }

    # Distance matrix: NA diagonal treated as Inf → zero weight (no self-loop)
    D_h <- D
    diag(D_h) <- NA  # exclude self

    # Weight w_ji = d(j→i)^(-alpha); treat 0 or NA distances as near-zero weight
    W_decay <- (D_h + 1)^(-alpha)   # +1 avoids infinite weights at d=0 for non-diagonal
    W_decay[is.na(W_decay)] <- 0

    # Row-normalise: mu_i = sum_j w_ji * Y_j / sum_j w_ji
    row_norm <- rowSums(W_decay)
    row_norm[row_norm == 0] <- 1  # avoid /0 for isolated zones
    W_norm <- W_decay / row_norm

    mu <- as.numeric(W_norm %*% Y_curr)
    mu[mu < 0] <- 0
    mu_prev <- mu   # marginal week-h intensity seeds the next week's projection
    # The LFO outcome is cumulative ("first case in weeks (cut, cut+h*7]"), so report the
    # CUMULATIVE hazard 1-exp(-sum_k mu_k) (as the workhorse/hhh4/SEIR do); the marginal
    # week-h probability understates invasion at h>=2. Count quantiles then match too.
    mu_cum <- if (h_idx == 1) mu else mu_cum + mu

    results[[h_idx]] <- tibble(
      health_zone  = ord,
      horizon      = h,
      mu_forecast  = mu_cum,
      p_invasion   = 1 - exp(-mu_cum),
      method       = paste0("B1_alpha", alpha)
    )
  }

  bind_rows(results)
}

# ---------------------------------------------------------------------------
# B2: K-nearest neighbour average
# ---------------------------------------------------------------------------

#' B2: Forecast via K-nearest-neighbour (by OSRM travel time) average
#'
#' @param K number of nearest neighbours (excluding self)
forecast_B2 <- function(Y_wide, osrm_mat, t_idx, horizons, zones_all, K = 3) {
  stopifnot(K >= 1)

  ord  <- zones_all[zones_all %in% rownames(osrm_mat) & zones_all %in% rownames(Y_wide)]
  D    <- osrm_mat[ord, ord]
  diag(D) <- NA  # exclude self

  results <- vector("list", length(horizons))
  Y_curr  <- Y_wide[ord, min(t_idx, ncol(Y_wide))]

  for (h_idx in seq_along(horizons)) {
    h <- horizons[h_idx]
    if (h_idx > 1) Y_curr <- mu_prev

    mu <- numeric(length(ord))
    for (i in seq_along(ord)) {
      dists     <- D[i, ]
      dists[is.na(dists)] <- Inf
      nn_idx    <- order(dists)[seq_len(min(K, sum(is.finite(dists))))]
      mu[i]     <- if (length(nn_idx) > 0) mean(Y_curr[nn_idx], na.rm = TRUE) else 0
    }
    mu[mu < 0] <- 0
    mu_prev    <- mu
    mu_cum     <- if (h_idx == 1) mu else mu_cum + mu   # cumulative hazard (matches LFO outcome)

    results[[h_idx]] <- tibble(
      health_zone = ord, horizon = h,
      mu_forecast = mu_cum, p_invasion = 1 - exp(-mu_cum),
      method = paste0("B2_K", K)
    )
  }
  bind_rows(results)
}

# ---------------------------------------------------------------------------
# B3: Spatial diffusion (discrete graph Laplacian)
# ---------------------------------------------------------------------------

#' B3: Spatial diffusion based on discrete graph Laplacian
#'
#' @param D_coef  diffusion coefficient in [0, 0.5]; if NULL, estimated from training data
#' @param threshold_min OSRM threshold (minutes) to define adjacency
forecast_B3 <- function(Y_wide, osrm_mat, t_idx, horizons, zones_all,
                        D_coef = NULL, threshold_min = 120) {
  ord  <- zones_all[zones_all %in% rownames(osrm_mat) & zones_all %in% rownames(Y_wide)]
  D    <- osrm_mat[ord, ord]
  n    <- length(ord)

  # Adjacency by travel-time threshold (binary)
  A <- (D < threshold_min & !is.na(D))
  diag(A) <- FALSE
  deg  <- rowSums(A)
  # Degree matrix (avoid div/0 for isolated zones)
  deg_safe <- pmax(deg, 1)
  # Row-normalised adjacency (each row sums to 1 for connected zones, 0 for isolates)
  W_adj <- sweep(A * 1.0, 1, deg_safe, "/")
  # Graph Laplacian: L = I - W_adj
  L <- diag(n) - W_adj

  # Estimate D_coef from last 2 training weeks via LOO RMSE minimisation
  if (is.null(D_coef) && t_idx >= 3) {
    obj <- function(d) {
      err_sq <- 0
      for (tau in (t_idx - 2):(t_idx - 1)) {
        Y_tau <- Y_wide[ord, tau]
        Y_next <- Y_wide[ord, tau + 1]
        mu_hat <- Y_tau - d * (L %*% Y_tau)
        mu_hat <- pmax(mu_hat, 0)
        err_sq <- err_sq + sum((Y_next - mu_hat)^2, na.rm = TRUE)
      }
      err_sq
    }
    opt    <- optimise(obj, interval = c(0, 0.45))
    D_coef <- opt$minimum
    message(sprintf("[B3] Estimated diffusion coefficient: %.4f", D_coef))
  } else if (is.null(D_coef)) {
    D_coef <- 0.05  # conservative default if insufficient history
  }

  results <- vector("list", length(horizons))
  Y_curr  <- Y_wide[ord, min(t_idx, ncol(Y_wide))]

  for (h_idx in seq_along(horizons)) {
    h <- horizons[h_idx]
    if (h_idx > 1) Y_curr <- mu_prev

    mu <- as.numeric(Y_curr - D_coef * (L %*% Y_curr))
    mu[mu < 0] <- 0
    mu_prev <- mu
    mu_cum  <- if (h_idx == 1) mu else mu_cum + mu   # cumulative hazard (matches LFO outcome)

    results[[h_idx]] <- tibble(
      health_zone = ord, horizon = h,
      mu_forecast = mu_cum, p_invasion = 1 - exp(-mu_cum),
      method = "B3"
    )
  }
  bind_rows(results)
}

# ---------------------------------------------------------------------------
# B4: Gravity-decay invasion probability
# ---------------------------------------------------------------------------

#' B4: Gravity-weighted invasion hazard (exponential arrival model)
#'
#' @param scaling_const multiplicative scaling; if NULL, estimated from training
forecast_B4 <- function(Y_wide, pop_vec, osrm_mat, t_idx, horizons, zones_all,
                        beta = 0.5, gamma_d = 1.0, scaling_const = NULL) {
  ord  <- zones_all[zones_all %in% rownames(osrm_mat) &
                    zones_all %in% rownames(Y_wide) &
                    zones_all %in% names(pop_vec)]
  D    <- osrm_mat[ord, ord]
  diag(D) <- NA
  N    <- pop_vec[ord]

  # Estimate scaling constant from historical training data if possible
  if (is.null(scaling_const) && t_idx >= 2) {
    # Build gravity scores for each past week and regress invasion probability
    # Using simple calibration: scale so that P(invasion) at mean lambda ≈ empirical rate
    lambdas <- numeric(0)
    invasions <- logical(0)
    for (tau in seq_len(t_idx - 1)) {
      Y_tau <- Y_wide[ord, tau]
      already_active <- rowSums(Y_wide[ord, seq_len(tau), drop = FALSE] > 0) > 0
      for (i in seq_along(ord)) {
        # Calibrate on the AT-RISK set only (matching the workhorse and the invasion
        # task). Including already-active zones as forced-negatives keeps their large
        # gravity score with invaded = FALSE, dragging the MLE scaling constant down
        # and biasing B4's invasion probability low.
        if (already_active[i]) next
        d_ij <- D[i, ]; d_ij[is.na(d_ij)] <- Inf
        lam_i <- sum(N^beta * (d_ij + 1)^(-gamma_d) * Y_tau, na.rm = TRUE)
        lambdas   <- c(lambdas, lam_i)
        invasions <- c(invasions,
                       Y_wide[ord[i], min(tau + 1, ncol(Y_wide))] > 0)
      }
    }
    # MLE: p = 1 - exp(-k * lambda) → find k by grid search
    k_grid <- 10^seq(-6, 0, length.out = 60)
    ll_grid <- vapply(k_grid, function(k) {
      p <- pmin(1 - exp(-k * lambdas), 1 - 1e-9)
      p <- pmax(p, 1e-9)
      sum(invasions * log(p) + (1 - invasions) * log(1 - p), na.rm = TRUE)
    }, numeric(1))
    scaling_const <- k_grid[which.max(ll_grid)]
    message(sprintf("[B4] Estimated gravity scaling constant: %.2e", scaling_const))
  } else if (is.null(scaling_const)) {
    scaling_const <- 1e-5  # conservative default
  }

  results <- vector("list", length(horizons))
  # Source incidence is held at the current observed week for EVERY horizon; the per-week gravity
  # hazard (computed once, since Y_curr is constant) is accumulated over the window below. The
  # count-scale sibling baselines (B1/B2/B3) feed a ROW-NORMALISED projection forward to model
  # growth, but B4's raw (unnormalised) gravity kernel cannot: the previous code fed the SCALED
  # hazard mu = lambda*scaling_const back in as "incidence", so week-2's force was ~scaling_const^2
  # (≈0) and the 2-week forecast collapsed to the 1-week value. A constant weekly hazard is the
  # scale-consistent baseline behaviour and keeps p(2wk) > p(1wk) as required.
  Y_curr <- Y_wide[ord, min(t_idx, ncol(Y_wide))]
  lambda <- numeric(length(ord))
  for (i in seq_along(ord)) {
    d_ij <- D[i, ]; d_ij[is.na(d_ij)] <- Inf
    lambda[i] <- sum(N^beta * (d_ij + 1)^(-gamma_d) * Y_curr, na.rm = TRUE)
  }
  mu_week <- pmax(lambda * scaling_const, 0)   # per-week gravity invasion hazard

  for (h_idx in seq_along(horizons)) {
    h <- horizons[h_idx]
    mu_cum <- mu_week * h   # cumulative hazard over h weeks (matches the LFO cumulative outcome)
    results[[h_idx]] <- tibble(
      health_zone = ord, horizon = h,
      mu_forecast = mu_cum, p_invasion = pmin(1 - exp(-mu_cum), 1),
      method = "B4"
    )
  }
  bind_rows(results)
}

# ---------------------------------------------------------------------------
# B5: Population-proportional null
# ---------------------------------------------------------------------------

#' B5: Distribute total national case count proportionally to population
forecast_B5 <- function(Y_wide, pop_vec, t_idx, horizons, zones_all) {
  ord <- zones_all[zones_all %in% rownames(Y_wide) & zones_all %in% names(pop_vec)]
  N   <- pop_vec[ord]; N_total <- sum(N)

  results <- vector("list", length(horizons))
  Y_curr  <- Y_wide[ord, min(t_idx, ncol(Y_wide))]

  for (h_idx in seq_along(horizons)) {
    h <- horizons[h_idx]
    if (h_idx > 1) Y_curr <- mu_prev

    Y_total <- sum(Y_curr, na.rm = TRUE)
    mu <- (N / N_total) * Y_total
    mu_prev <- mu
    mu_cum  <- if (h_idx == 1) mu else mu_cum + mu   # cumulative hazard (matches LFO outcome)

    results[[h_idx]] <- tibble(
      health_zone = ord, horizon = h,
      mu_forecast = mu_cum, p_invasion = 1 - exp(-mu_cum),
      method = "B5"
    )
  }
  bind_rows(results)
}

# ---------------------------------------------------------------------------
# B6: Persistent-zero for uninvaded zones
# ---------------------------------------------------------------------------

#' B6: Zero forecast for zones never yet invaded; persist last count for active zones
#'
#' @param currently_active character vector of zones with >= 1 confirmed case to date
forecast_B6 <- function(Y_wide, t_idx, horizons, zones_all, currently_active) {
  ord <- zones_all[zones_all %in% rownames(Y_wide)]

  results <- vector("list", length(horizons))
  Y_curr  <- Y_wide[ord, min(t_idx, ncol(Y_wide))]

  for (h_idx in seq_along(horizons)) {
    h <- horizons[h_idx]
    mu <- ifelse(ord %in% currently_active, Y_curr, 0)
    mu[mu < 0] <- 0

    results[[h_idx]] <- tibble(
      health_zone = ord, horizon = h,
      mu_forecast = mu,
      p_invasion  = ifelse(ord %in% currently_active, 1 - exp(-mu), 0),
      method      = "B6"
    )
  }
  bind_rows(results)
}

# ---------------------------------------------------------------------------
# Master baseline runner
# ---------------------------------------------------------------------------
# B7: Nearest-affected / adjacency spatial-spread baseline  (review §3.1)
# ---------------------------------------------------------------------------

#' B7: Rank at-risk zones purely by proximity to the already-affected set.
#'
#' This is the "does the virus just go next door?" null the reviewer asked for
#' (2026-08-06 review, §3.1). It carries NO mobility-network information and NO
#' case-magnitude information — only how close each candidate zone is to the
#' current front — so BEATING it (not merely beating random chance) is the real
#' evidence that the mobility kernel adds signal beyond spatial contiguity.
#'
#' Two variants (method):
#'   * "distance" (default): score decreasing in the distance to the nearest
#'     already-affected zone, d_min_i = min_{j in affected} D[i, j]. We report a
#'     bounded, monotone score s_i = 1 / (1 + d_min_i); all rank-based metrics
#'     (AUC-PR, AUC-ROC, mean-rank-of-truth, detection curve) depend only on the
#'     ORDER, which is exact regardless of the monotone transform. `osrm_mat`
#'     (travel time) or a great-circle matrix both work.
#'   * "adjacency": score = number of already-affected first-order neighbours,
#'     using a supplied 0/1 contiguity matrix `adj` (shared-border), with
#'     1/(1+d_min) as the continuous tie-breaker. Requires `adj`; falls back to
#'     "distance" with a warning when `adj` is NULL.
#'
#' Already-affected zones are the origin, never an at-risk target, so they are
#' scored 0 (masked out of the invasion denominator downstream, as elsewhere).
#' The ranking is static across horizons (a pure spatial null), so the same score
#' is emitted for every horizon — matching how `naive_epicentre_inflow_scores`
#' and the detection-curve machinery treat a rank-only comparator.
#'
#' @param Y_wide   wide case matrix (n_zones × n_weeks); affected = any case up to t_idx.
#' @param dist_mat zone×zone distance matrix (OSRM travel time or great-circle),
#'                 rownames = colnames = zones. Diagonal ignored.
#' @param t_idx    integer week index (1-based) of the latest training observation.
#' @param horizons integer vector of horizons (weeks).
#' @param zones_all canonical zone order.
#' @param method   "distance" (default) or "adjacency".
#' @param adj      optional 0/1 contiguity matrix (zones × zones) for "adjacency".
#' @return tibble: health_zone, horizon, mu_forecast, p_invasion, method.
forecast_B7_adjacency <- function(Y_wide, dist_mat, t_idx, horizons, zones_all,
                                  method = c("distance", "adjacency"), adj = NULL) {
  method <- match.arg(method)
  stopifnot(is.matrix(Y_wide), is.matrix(dist_mat))

  # Affected = zones with any case in weeks 1..t_idx (same rule as the LFO at-risk set).
  tt        <- min(t_idx, ncol(Y_wide))
  Y_hist    <- Y_wide[, seq_len(tt), drop = FALSE]
  affected  <- rownames(Y_wide)[rowSums(Y_hist > 0, na.rm = TRUE) > 0]

  # Align to the zones present in BOTH the distance matrix and the canonical spine.
  ord      <- zones_all[zones_all %in% rownames(dist_mat)]
  D        <- dist_mat[ord, ord, drop = FALSE]
  aff_here <- intersect(affected, ord)

  score <- stats::setNames(rep(0, length(ord)), ord)
  if (length(aff_here) >= 1L) {
    # Distance to the nearest already-affected zone (exclude self by construction:
    # affected zones get score 0 below regardless).
    Dsub   <- D[, aff_here, drop = FALSE]
    diag_i <- match(aff_here, ord)                    # affected columns that are also rows
    for (cc in seq_along(aff_here)) {                 # blank self-distance so it is never the min
      ri <- match(aff_here[cc], ord); if (!is.na(ri)) Dsub[ri, cc] <- NA_real_
    }
    d_min  <- apply(Dsub, 1, function(z) { z <- z[is.finite(z)]; if (length(z)) min(z) else NA_real_ })
    s_dist <- 1 / (1 + d_min)
    s_dist[!is.finite(s_dist)] <- 0

    if (method == "adjacency" && !is.null(adj)) {
      A          <- adj[ord, ord, drop = FALSE]
      n_aff_nbr  <- as.numeric(A[, aff_here, drop = FALSE] %*% rep(1, length(aff_here)))
      # rank primarily by #affected neighbours, break ties by 1/(1+d_min):
      score[]    <- n_aff_nbr + s_dist / (max(s_dist, na.rm = TRUE) + 1)  # tie-breaker < 1
    } else {
      if (method == "adjacency")
        warning("[B7] adjacency requested but no contiguity matrix supplied; using distance.")
      score[] <- s_dist
    }
  }
  score[intersect(ord, affected)] <- 0                # affected zones are not at-risk targets

  meth_lab <- paste0("B7_", method)
  purrr::map_dfr(horizons, function(h) tibble::tibble(
    health_zone = ord,
    horizon     = h,
    mu_forecast = as.numeric(score),   # monotone ranking key (NOT a calibrated intensity)
    p_invasion  = as.numeric(score),   # rank-based metrics use the order only
    method      = meth_lab
  ))
}

# ---------------------------------------------------------------------------

#' Run all baseline models B1–B6 for a given training cutoff
#'
#' @param zone_week_nc tibble (health_zone, week_start, confirmed_nc, confirmed)
#' @param osrm_mat     travel-time matrix
#' @param pop_vec      named population vector
#' @param zones_all    canonical zone name vector
#' @param t_idx        week index of latest training observation
#' @param horizons     integer vector of horizons (typically c(1L, 2L))
#' @param training_cutoff Date of the final training week (WEEK_ANCHOR-anchored start)
#' @return combined tibble of all baseline forecasts
run_baselines <- function(zone_week_nc, osrm_mat, pop_vec, zones_all,
                          t_idx, horizons = LFO_HORIZONS,
                          training_cutoff = NULL) {
  Y_wide <- zone_week_to_wide(zone_week_nc, zones_all)

  # Zones with any confirmed case up to and including t_idx
  currently_active <- zone_week_nc |>
    filter(week_start <= sort(unique(zone_week_nc$week_start))[t_idx],
           !is.na(confirmed) & confirmed > 0) |>
    pull(health_zone) |> unique()

  message("[baselines] Running B1–B6 for ", length(currently_active),
          " active zones, t_idx=", t_idx)

  purrr::map_dfr(
    list(
      forecast_B1(Y_wide, osrm_mat, t_idx, horizons, zones_all, alpha = 1.0),
      forecast_B1(Y_wide, osrm_mat, t_idx, horizons, zones_all, alpha = 0.5),
      forecast_B2(Y_wide, osrm_mat, t_idx, horizons, zones_all, K = 3),
      forecast_B2(Y_wide, osrm_mat, t_idx, horizons, zones_all, K = 5),
      forecast_B3(Y_wide, osrm_mat, t_idx, horizons, zones_all),
      forecast_B4(Y_wide, pop_vec, osrm_mat, t_idx, horizons, zones_all),
      forecast_B5(Y_wide, pop_vec, t_idx, horizons, zones_all),
      forecast_B6(Y_wide, t_idx, horizons, zones_all, currently_active)
    ),
    ~ .x
  ) |>
    # Baselines are point forecasts; attach a Poisson predictive distribution
    # (the natural count-observation model) so they are WIS-scorable alongside
    # the stochastic models. Without these the baselines return NA WIS.
    mutate(
      q05 = qpois(0.05, pmax(mu_forecast, 1e-6)),
      q20 = qpois(0.20, pmax(mu_forecast, 1e-6)),
      q25 = qpois(0.25, pmax(mu_forecast, 1e-6)),
      q75 = qpois(0.75, pmax(mu_forecast, 1e-6)),
      q80 = qpois(0.80, pmax(mu_forecast, 1e-6)),
      q95 = qpois(0.95, pmax(mu_forecast, 1e-6)),
      training_cutoff = if (!is.null(training_cutoff)) training_cutoff else as.Date(NA),
      mobility_id     = "none",
      gt_profile      = "none"
    )
}
