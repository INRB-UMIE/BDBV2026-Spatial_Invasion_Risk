# =============================================================================
# 05_baseline_models.R - Structural Baselines (B1, B4, B7)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
#   Models: B1 inverse-distance weighted average, B4 gravity-decay invasion hazard,
#           B7 nearest-affected / adjacency spatial-spread null.
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

  # WEEK BY WEEK, accumulating the hazard. The outcome this is scored against is cumulative
  # ("first case in weeks (cut, cut + h*7]"), so the reported intensity must be
  # mu_cum(h) = sum_{k=1..h} mu_k.
  #
  # It used to be `mu * h` — h copies of the h-th STEP. That is only the cumulative sum when
  # mu is constant across the window, which is true of B4 (whose source incidence is held
  # fixed) but NOT of B1: B1 deliberately feeds the previous step's intensity forward
  # (Y_curr <- mu_prev) to model growth, so every step differs. At h = 2 the baseline
  # reported 2*mu_2 instead of mu_1 + mu_2. Reconstructed from the saved LFO, that
  # UNDERSTATED the baseline's own ranking (AUC-PR 0.080 -> 0.124, skill 6.51 -> 10.13,
  # AUC-ROC 0.902 -> 0.931) while OVERSTATING its hazard (median ratio 1.19) — so the
  # featured model's published h=2 margin over Distance-B1 was inflated by ~55% in AUC-PR
  # skill, in the direction that flatters the model this paper is about.
  #
  # Looping over WEEKS (not over the horizon index) also makes a non-consecutive horizon set
  # correct: with horizons c(1, 2, 4), week 4's entry accumulates four steps, not three.
  h_max   <- max(horizons)
  Y_curr  <- Y_t
  mu_cum  <- numeric(length(ord))
  by_week <- vector("list", h_max)

  # Distance matrix: NA diagonal treated as Inf -> zero weight (no self-loop). The decay
  # weights do not depend on the week, so they are built once outside the loop.
  D_h <- D
  diag(D_h) <- NA  # exclude self
  # Weight w_ji = d(j->i)^(-alpha); treat 0 or NA distances as near-zero weight
  W_decay <- (D_h + 1)^(-alpha)   # +1 avoids infinite weights at d=0 for non-diagonal
  W_decay[is.na(W_decay)] <- 0
  # Row-normalise: mu_i = sum_j w_ji * Y_j / sum_j w_ji
  row_norm <- rowSums(W_decay)
  row_norm[row_norm == 0] <- 1  # avoid /0 for isolated zones
  W_norm <- W_decay / row_norm

  for (wk in seq_len(h_max)) {
    mu <- as.numeric(W_norm %*% Y_curr)
    mu[mu < 0] <- 0
    mu_cum <- mu_cum + mu
    Y_curr <- mu            # marginal week intensity seeds the next week's projection
    by_week[[wk]] <- mu_cum
  }

  results <- vector("list", length(horizons))
  for (h_idx in seq_along(horizons)) {
    h  <- horizons[h_idx]
    mc <- by_week[[h]]
    results[[h_idx]] <- tibble(
      health_zone  = ord,
      horizon      = h,
      mu_forecast  = mc,
      p_invasion   = pmin(1 - exp(-mc), 1),
      method       = paste0("B1_alpha", alpha),
      # DECLARED RANK-ONLY. mu is a row-normalised weighted MEAN OF NEIGHBOUR CASE COUNTS
      # fed through 1 - exp(-mu) as though a case count were a weekly invasion rate. Unlike
      # B4, nothing here is fitted — there is no scale constant at all — so the probability
      # level is arbitrary and drifts upward simply as the epidemic grows (mean p rose
      # 0.22 -> 0.43 across folds). The published calibration_in_large of ~45 and
      # brier_skill of -17 measure that missing constant, not the baseline's information.
      # The ORDERING is meaningful and is what this comparator exists to provide.
      prob_calibrated = FALSE
    )
  }

  bind_rows(results)
}

# ---------------------------------------------------------------------------
# B4: source-mass x distance-decay import pressure
# ---------------------------------------------------------------------------

#' B4: SOURCE-MASS-WEIGHTED DISTANCE-DECAY import pressure — NOT a gravity model.
#'
#' lambda_i = sum_j  N_j^beta * (d_ij + 1)^(-gamma_d) * Y_j,  then mu = lambda * scaling_const.
#'
#' WHAT IT IS NOT. A gravity kernel is N_i^alpha N_j^beta / d_ij^gamma: it carries a
#' DESTINATION mass. This has none — `N` is indexed by j, the source, only. A city of a
#' million and a zone of twenty thousand at equal travel time from the same source therefore
#' receive IDENTICAL hazard, which is the one prediction a gravity model exists to make. The
#' exponents (beta = 0.5, gamma_d = 1.0) are hard-coded defaults and are never estimated —
#' only the scalar `scaling_const` is fitted (by MLE, below), and a scalar cannot change the
#' RANKING. So the ordering this comparator supplies is entirely unparameterised.
#'
#' Retaining it is deliberate: an unparameterised structural yardstick is exactly what a
#' baseline should be. What it must not do is borrow the name of a model it is not. This
#' pipeline HAS a fitted gravity kernel — M4, a censored gravity GLM with estimated exponents
#' (03_mobility_matrices.R) — and publishing this alongside it as "Gravity-B4" invited the
#' reader to believe the comparator was that model, or at least its species.
#'
#' NOTE ON THE NAME. The method STRING stays "Gravity-B4": it is a key in the cached LFO
#' results, the evaluation CSVs and the model-selection artifacts, and changing it would
#' silently orphan them. What changed is every label a reader sees, and this docstring.
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
    .best <- which.max(ll_grid)
    scaling_const <- k_grid[.best]
    # BOUNDARY CHECK: the grid spans 1e-6..1, so an optimum sitting ON either end is not
    # an interior maximum — the likelihood is still climbing and the "MLE" is really the
    # grid edge, silently clipped. Say so rather than reporting a clipped value as a fit.
    if (.best == 1L || .best == length(k_grid))
      warning(sprintf(paste0("[B4] The gravity scaling constant hit the %s edge of the search ",
                             "grid (k = %.2e, grid 1e-6..1): this is NOT an interior maximum, ",
                             "so the likelihood is still improving outside the grid and the ",
                             "value is a clipped bound, not an MLE. Widen k_grid."),
                     if (.best == 1L) "LOWER" else "UPPER", scaling_const), call. = FALSE)
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
      method = "B4",
      # B4 DOES carry a fitted scale (scaling_const, MLE'd above), so its p_invasion is on a
      # probability scale and the proper scores are meaningful for it.
      prob_calibrated = TRUE
    )
  }
  bind_rows(results)
}

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
    method      = meth_lab,
    # DECLARED RANK-ONLY. score = 1/(1 + travel-time MINUTES) is an arbitrary monotone
    # transform of distance, not a probability: 1/(1 + d^2) would order the zones
    # identically and change every calibration-dependent number at will. The docstring
    # said so, but evaluate_invasion() still ran .log_score/.brier_skill/.ece/
    # calibration_in_large on whatever sat in p_invasion, and 16b fitted a recalibration
    # delta to it. This flag makes the declaration machine-readable so those metrics are
    # reported NA instead of as if they measured information. The score is also identical
    # at h=1 and h=2 while the base rate doubles, which alone guarantees miscalibration.
    prob_calibrated = FALSE
  ))
}



