# =============================================================================
# 15_workhorse.R — Mobility-Informed Renewal Workhorse for Spatial Invasion
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# The primary deliverable: for each health zone with NO confirmed cases up to the
# forecast date, the probability it sees its first case(s) in the next 1-2 weeks
# (spatial INVASION), plus a set of relative-risk scores. Zones already affected
# carry NO invasion probability (NA everywhere).
#
# Model (renewal / force-of-infection, verified core in 06_simple_models.R):
#   Import force to at-risk zone i:  Lambda_i = sum_j W[j,i] * sum_k g(k) Y_nc[j,t-k]
#   Expected introductions:          mu_i     = beta * Lambda_i
#   Invasion probability:            p_i      = 1 - exp(-mu_i)     (Poisson arrival)
#                                    p_i      = 1 - (theta/(theta+mu_i))^theta  (NegBin)
#   Affected zones:                  p_i      = NA
#
# The import coefficient beta is LEARNED from realised invasion events via a
# complementary-log-log regression pooled over training weeks:
#   cloglog(P(invade)) = log(beta) + log(Lambda)   [offset = log Lambda]
# so p is tied directly to the mobility-driven import force and calibrated to the
# actual (rare) invasion base rate — fixing the previous 1-exp(-R*Lambda)
# overconfidence (R is a self-sustaining transmission number, NOT an import scale).
#
# Risk scores (at-risk zones only): absolute p_case, ascertainment-adjusted
# p_infection (as a band across the rho grid), relative risk within each province
# of interest (Ituri, Nord-Kivu, Haut-Uele) and relative risk nationwide.
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))
suppressPackageStartupMessages({ library(tidyverse) })

# ---------------------------------------------------------------------------
# Province / Ituri identification
# ---------------------------------------------------------------------------

#' Map canonical health-zone name -> province from the shapefile.
#'
#' @return tibble(nom, province); Ituri zones are those with province == "Ituri".
load_province_map <- function() {
  shp_path <- SHAPEFILE_PATH
  if (!file.exists(shp_path)) {
    warning("[workhorse] shapefile not found; province map unavailable.")
    return(NULL)
  }
  shp <- suppressWarnings(sf::st_read(shp_path, quiet = TRUE))
  sf::st_geometry(shp) <- NULL
  if (!all(c("Nom", "PROVINCE") %in% names(shp))) {
    warning("[workhorse] shapefile missing Nom/PROVINCE columns.")
    return(NULL)
  }
  pm <- shp %>%
    dplyr::transmute(nom = as.character(Nom), province = as.character(PROVINCE),
                     zscode = if ("ZSCode" %in% names(shp)) as.character(ZSCode) else NA) %>%
    dplyr::filter(!is.na(nom))
  # A few health-zone NAMES recur in different provinces (e.g. Bili, Lubunga). This nom->province
  # lookup must return one province per name, so it keeps the first occurrence — but WARN rather
  # than silently pick, since the choice is arbitrary. (Maps disambiguate via a (name, province)
  # join, so only free-text province labels for these zones are affected — all far from Ituri.)
  dup <- pm %>% dplyr::distinct(nom, province) %>% dplyr::count(nom) %>% dplyr::filter(n > 1)
  if (nrow(dup))
    warning(sprintf("[workhorse] %d health-zone name(s) span multiple provinces (%s); province label uses the first shapefile occurrence.",
                    nrow(dup), paste(dup$nom, collapse = ", ")))
  pm %>% dplyr::distinct(nom, .keep_all = TRUE) %>% dplyr::select(nom, province)
}

# ---------------------------------------------------------------------------
# At-risk mask
# ---------------------------------------------------------------------------

#' Zones AFFECTED (>=1 confirmed case) in any week up to and including `cutoff`.
#'
#' @param zone_week tibble with health_zone, week_start, confirmed.
#' @param cutoff    Date; the forecast-generation (training cutoff) week.
#' @return character vector of affected zone names.
affected_zones <- function(zone_week, cutoff) {
  zone_week %>%
    dplyr::filter(week_start <= cutoff, !is.na(confirmed), confirmed > 0) %>%
    dplyr::pull(health_zone) %>% unique()
}

# ---------------------------------------------------------------------------
# Learn the import coefficient beta from realised invasions
# ---------------------------------------------------------------------------

#' Fit the import coefficient beta by cloglog regression of invasion events on
#' the mobility import force, pooled over all training weeks.
#'
#' For each training week t' (with a valid next week) and each zone at-risk at t'
#' (no cases in weeks <= t'), we form (Lambda_i(t'), invaded_i(t'+1)) and fit
#'   cloglog(P(invade)) = intercept + offset(log Lambda),  beta = exp(intercept).
#' Optionally augments with standardised covariates (penalised via a ridge-like
#' data augmentation) — off by default given the small number of events.
#'
#' @return list(beta, coef_cov = named vector or NULL, n_events, n_obs).
fit_import_beta <- function(Y_wide, W, G_weekly, zones_all,
                            covariate_mat = NULL, lambda_ridge = 1,
                            week_weights = NULL) {
  n_weeks <- ncol(Y_wide)
  rows <- list()
  for (t in seq_len(n_weeks - 1L)) {
    Lambda_t <- compute_foi(Y_wide, W, G_weekly, t_idx = t + 1L, zones_all)
    # affected through week t (cumulative)
    aff <- rowSums(Y_wide[, seq_len(t), drop = FALSE] > 0) > 0
    at_risk <- !aff & Lambda_t > 0        # need positive import force to be estimable
    if (!any(at_risk)) next
    invaded <- Y_wide[, t + 1L] > 0
    # Q1 completeness weight: reliability of the outcome week t+1 (default 1).
    w_t <- if (is.null(week_weights)) 1 else
      week_weights[[min(t + 1L, length(week_weights))]]
    rows[[length(rows) + 1L]] <- tibble::tibble(
      zone    = zones_all[at_risk],
      logLam  = log(Lambda_t[at_risk]),
      invaded = as.integer(invaded[at_risk]),
      w       = as.numeric(w_t)
    )
  }
  d <- dplyr::bind_rows(rows)
  n_events <- sum(d$invaded)
  if (nrow(d) < 20L || n_events < 2L) {
    # Fall back to a completeness-weighted moment estimate:
    # beta = sum(w·events) / sum(w·Lambda over at-risk).
    beta <- if (nrow(d) > 0)
      max(sum(d$w * d$invaded) / sum(d$w * exp(d$logLam)), 1e-4) else 0.01
    return(list(beta = beta, coef_cov = NULL, n_events = n_events, n_obs = nrow(d)))
  }

  fit <- tryCatch(
    suppressWarnings(              # benign non-integer-weights note (see fit_import_model)
      glm(invaded ~ 1, family = binomial(link = "cloglog"),
          offset = d$logLam, data = d, weights = d$w)),
    error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(coef(fit)[1])) {
    beta <- max(sum(d$w * d$invaded) / sum(d$w * exp(d$logLam)), 1e-4)
    return(list(beta = beta, coef_cov = NULL, n_events = n_events, n_obs = nrow(d)))
  }
  beta <- exp(unname(coef(fit)[1]))
  list(beta = beta, coef_cov = NULL, n_events = n_events, n_obs = nrow(d))
}

# ---------------------------------------------------------------------------
# Covariate-augmented import model (Q2 covariates + Q6 reporting-rate)
# ---------------------------------------------------------------------------
# The invasion hazard becomes  p_i = 1 - exp(-beta_i * Lambda_i)  with a
# ZONE-SPECIFIC scale  beta_i = exp(intercept + sum_m gamma_m * x_{m,i}), i.e. we
# keep the mobility import force Lambda as a fixed OFFSET and let penalised
# covariates modulate the per-zone susceptibility. When covariate_spec is NULL
# this collapses exactly to the intercept-only fit_import_beta (beta_i = beta).
# Coefficients are Firth-penalised (brglm2) so a handful of events cannot make
# them explode (the failure mode of the retired free-coefficient S3 model).

#' Build a zones x weeks matrix from an arbitrary count column.
.count_wide <- function(zone_week, zones_all, col) {
  weeks <- sort(unique(zone_week$week_start))
  mat <- matrix(0, length(zones_all), length(weeks),
                dimnames = list(zones_all, as.character(weeks)))
  sub <- zone_week[zone_week$health_zone %in% zones_all, ]
  if (!col %in% names(sub)) return(mat)
  for (k in seq_along(weeks)) {
    s <- sub[sub$week_start == weeks[k], c("health_zone", col)]
    v <- s[[col]]; v[is.na(v)] <- 0
    mat[s$health_zone, k] <- v
  }
  mat
}

#' zones x weeks matrix of the suspected-but-not-confirmed leading-indicator counts.
#' Prefers the nowcast-corrected `suspected_nc` (consistent with the confirmed_nc primary
#' signal); falls back to the RAW `suspected` count when the corrected column is absent —
#' e.g. the raw-nowcast SENSITIVITY variant (.zw_raw in run_all), which by construction
#' carries confirmed_nc but no other _nc columns. Without this fallback a susp-covariate
#' model refit on that variant would silently see an all-zero (hence dropped) covariate.
.susp_wide <- function(zone_week, zones_all, prefer = "suspected_nc") {
  col <- if (prefer %in% names(zone_week)) prefer
         else if ("suspected" %in% names(zone_week)) "suspected"
         else prefer
  .count_wide(zone_week, zones_all, col)
}

#' Named list of standardised STATIC covariate vectors aligned to zones_all.
.static_features <- function(static_cov, zones_all) {
  out <- list()
  if (is.null(static_cov)) return(out)
  g <- function(nm, tf = identity) {
    if (!nm %in% names(static_cov)) return(NULL)
    v <- static_cov[[nm]][match(zones_all, static_cov$nom)]
    v <- tf(as.numeric(v)); v[!is.finite(v)] <- NA_real_
    v[is.na(v)] <- stats::median(v, na.rm = TRUE); v
  }
  out$log_pop            <- g("log_pop")
  if (is.null(out$log_pop)) out$log_pop <- g("pop_count", function(x) log(pmax(x, 1)))
  out$healthsite_density <- g("healthsite_density")
  out$ccvi               <- g("ccvi")
  if (is.null(out$ccvi)) out$ccvi <- g("socioeconomic_deprivation")
  out$positivity         <- g("positivity")
  out[!vapply(out, is.null, logical(1))]
}

#' Relative reporting-rate vector (Q6 reporting-rate structure).
#'
#' Turns a raw completeness proxy (e.g. health-site density: more infrastructure
#' -> more complete confirmed-case ascertainment) into a per-zone relative
#' reporting rate r_j, normalised to GEOMETRIC mean 1 (so it re-weights sources
#' relative to one another without rescaling the overall import force — a global
#' scale is absorbed by the learned intercept beta and is NOT identifiable) and
#' bounded to [0.25, 4] to prevent any single under-reporting source from
#' dominating. Source incidence is later inflated by 1/r_j, i.e. low-reporting
#' zones exert MORE importation pressure than their observed counts suggest.
.reporting_rate_vec <- function(report_rate, zones_all) {
  if (is.null(report_rate)) return(NULL)
  v <- as.numeric(report_rate[match(zones_all, names(report_rate))])
  v[!is.finite(v) | v <= 0] <- NA_real_
  gm <- exp(mean(log(v[is.finite(v) & v > 0])))          # geometric mean of valid
  if (!is.finite(gm) || gm <= 0) return(NULL)
  v[is.na(v)] <- gm                                       # impute missing at neutral
  r <- v / gm
  r <- pmin(pmax(r, 0.25), 4)                             # bound extreme adjustments
  stats::setNames(r, zones_all)                           # named for safe alignment
}

#' log(1 + OSRM travel time to the nearest zone active through week `t`).
#' Infinite (no active zone / unroutable) -> a large finite sentinel. `t` is
#' capped at ncol(Y_wide) so callers can pass the OBSERVED incidence matrix with a
#' future week index (the observed frontier is fixed over the forecast window).
.dmin_vec <- function(Y_wide, osrm, zones_all, t) {
  if (is.null(osrm)) return(rep(0, length(zones_all)))
  t <- min(max(t, 1L), ncol(Y_wide))
  active <- zones_all[rowSums(Y_wide[, seq_len(t), drop = FALSE] > 0) > 0]
  active <- intersect(active, colnames(osrm))
  if (length(active) == 0) return(rep(log1p(1e4), length(zones_all)))
  sub <- osrm[intersect(zones_all, rownames(osrm)), active, drop = FALSE]
  dmin <- apply(sub, 1, function(x) { x <- x[is.finite(x) & x > 0]
    if (length(x)) min(x) else 1e4 })
  v <- setNames(rep(1e4, length(zones_all)), zones_all)
  v[names(dmin)] <- dmin
  log1p(as.numeric(v[zones_all]))
}

#' Feature matrix (zones x features) for predicting invasion in week `t_for`,
#' using only information available through week t_for-1 (no leakage).
#'
#' @param Y_active incidence matrix used for the d_min frontier feature. At fit
#'   time this equals Y_wide (per-week historical frontier); at forecast time the
#'   caller passes the OBSERVED training matrix so d_min is the distance to the
#'   real frontier — NOT to zones with a small fractional PROJECTED incidence,
#'   which would otherwise all count as "active" from horizon 2 onward.
.feature_matrix <- function(t_for, Y_wide, A_wide, W, G, static, osrm,
                            zones_all, cov_spec, Y_active = Y_wide, S_wide = NULL) {
  f <- list()
  if ("alert_import" %in% cov_spec && !is.null(A_wide))
    f$alert_import <- log1p(compute_foi(A_wide, W, G, t_idx = t_for, zones_all))
  if ("alert_local" %in% cov_spec && !is.null(A_wide))
    f$alert_local  <- log1p(.gweighted_own(A_wide, G, t_for))
  # Suspected-but-not-confirmed leading indicators (S_wide = suspected-case counts, a
  # zones x weeks matrix built exactly like the alert matrix): susp_import is the
  # mobility-weighted import of OTHER zones' preceding suspected cases (same FOI kernel
  # as the confirmed import force), susp_local the own g-weighted preceding suspected
  # count. Both read only weeks < t_for (compute_foi / .gweighted_own index columns
  # t_for-1 and earlier), so they are past-only and causal in LFO — mirroring alerts.
  if ("susp_import" %in% cov_spec && !is.null(S_wide))
    f$susp_import <- log1p(compute_foi(S_wide, W, G, t_idx = t_for, zones_all))
  if ("susp_local" %in% cov_spec && !is.null(S_wide))
    f$susp_local  <- log1p(.gweighted_own(S_wide, G, t_for))
  for (nm in c("log_pop", "healthsite_density", "ccvi", "positivity"))
    if (nm %in% cov_spec && !is.null(static[[nm]])) f[[nm]] <- static[[nm]]
  if ("d_min" %in% cov_spec) f$d_min <- .dmin_vec(Y_active, osrm, zones_all, t_for - 1L)
  if (length(f) == 0) return(NULL)
  m <- do.call(cbind, f); colnames(m) <- names(f); m
}

#' Fit the (optionally covariate-augmented) import model.
#'
#' @return list(intercept, coef [named, on standardised scale], center, scale,
#'   cov_spec, n_events, n_obs). With cov_spec=NULL, coef is empty and
#'   exp(intercept) is the scalar beta of fit_import_beta.
fit_import_model <- function(Y_wide, W, G_weekly, zones_all,
                             A_wide = NULL, static = NULL, osrm = NULL,
                             cov_spec = NULL, week_weights = NULL, S_wide = NULL) {
  if (is.null(cov_spec) || length(cov_spec) == 0) {
    b <- fit_import_beta(Y_wide, W, G_weekly, zones_all, week_weights = week_weights)
    return(list(intercept = log(max(b$beta, 1e-8)), coef = numeric(0),
                center = numeric(0), scale = numeric(0), cov_spec = NULL,
                n_events = b$n_events, n_obs = b$n_obs))
  }
  n_weeks <- ncol(Y_wide); rows <- list()
  for (t in seq_len(n_weeks - 1L)) {
    Lambda_t <- compute_foi(Y_wide, W, G_weekly, t_idx = t + 1L, zones_all)
    aff <- rowSums(Y_wide[, seq_len(t), drop = FALSE] > 0) > 0
    at_risk <- !aff & Lambda_t > 0
    if (!any(at_risk)) next
    X <- .feature_matrix(t + 1L, Y_wide, A_wide, W, G_weekly, static, osrm,
                         zones_all, cov_spec, S_wide = S_wide)
    # Q1 completeness weighting: each transition t->t+1 is weighted by the
    # reporting completeness of its OUTCOME week t+1 (recent, right-truncated
    # weeks carry a less certain invasion outcome and so count for less).
    w_t <- if (is.null(week_weights)) 1 else
      week_weights[[min(t + 1L, length(week_weights))]]
    df <- data.frame(invaded = as.integer(Y_wide[at_risk, t + 1L] > 0),
                     logLam  = log(Lambda_t[at_risk]),
                     .w = as.numeric(w_t))
    if (!is.null(X)) df <- cbind(df, as.data.frame(X[at_risk, , drop = FALSE]))
    rows[[length(rows) + 1L]] <- df
  }
  d <- dplyr::bind_rows(rows)
  feat_cols <- setdiff(names(d), c("invaded", "logLam", ".w"))
  n_events <- sum(d$invaded)
  # standardise features on the training set
  center <- vapply(feat_cols, function(c) mean(d[[c]], na.rm = TRUE), numeric(1))
  scale  <- vapply(feat_cols, function(c) { s <- stats::sd(d[[c]], na.rm = TRUE)
    if (!is.finite(s) || s == 0) 1 else s }, numeric(1))
  for (c in feat_cols) d[[c]] <- (d[[c]] - center[[c]]) / scale[[c]]

  # too few events -> intercept-only moment estimate (ignore covariates)
  if (nrow(d) < 25L || n_events < 3L) {
    b <- fit_import_beta(Y_wide, W, G_weekly, zones_all, week_weights = week_weights)
    return(list(intercept = log(max(b$beta, 1e-8)), coef = numeric(0),
                center = numeric(0), scale = numeric(0), cov_spec = NULL,
                n_events = n_events, n_obs = nrow(d)))
  }
  fml <- stats::as.formula(paste("invaded ~", paste(feat_cols, collapse = " + ")))
  fit <- tryCatch(
    # suppressWarnings: fractional completeness weights trip the benign
    # "non-integer #successes" note; the fit is a valid weighted likelihood and
    # the coefficient sanity guard below independently rejects any bad fit.
    suppressWarnings(
      stats::glm(fml, family = binomial("cloglog"), offset = d$logLam, data = d,
                 weights = d$.w,
                 method = if (requireNamespace("brglm2", quietly = TRUE))
                   brglm2::brglmFit else "glm.fit")),
    error = function(e) NULL)
  cf <- if (is.null(fit)) NULL else coef(fit)
  # Reject a diverged fit: with standardised features and Firth penalisation a
  # sane coefficient is O(1) and the base-rate cloglog intercept is ~ -3; anything
  # huge (in ANY coefficient, intercept included, since exp(intercept)=beta0 feeds
  # the dispersion) means separation / collinearity with the import offset (e.g.
  # alerts ~ confirmed). Fall back to intercept-only.
  bad <- is.null(cf) || any(!is.finite(cf)) || any(abs(cf) > 15)
  if (bad) {
    b <- fit_import_beta(Y_wide, W, G_weekly, zones_all, week_weights = week_weights)
    return(list(intercept = log(max(b$beta, 1e-8)), coef = numeric(0),
                center = numeric(0), scale = numeric(0), cov_spec = NULL,
                n_events = n_events, n_obs = nrow(d)))
  }
  # Coefficient standard errors (for parameter visualisation / CIs). Firth-
  # penalised vcov where available; NA if the (co)variance is unavailable.
  se_all <- tryCatch(sqrt(diag(stats::vcov(fit))), error = function(e) NULL)
  list(intercept = unname(cf[1]), coef = cf[feat_cols],
       intercept_se = if (!is.null(se_all)) unname(se_all[1]) else NA_real_,
       se = if (!is.null(se_all)) se_all[feat_cols] else
            stats::setNames(rep(NA_real_, length(feat_cols)), feat_cols),
       center = center, scale = scale, cov_spec = feat_cols,
       n_events = n_events, n_obs = nrow(d))
}

#' Per-zone import scale beta_i = exp(intercept + standardised-features %*% coef).
.beta_vector <- function(model, t_for, Y_wide, A_wide, W, G, static, osrm,
                         zones_all, Y_active = Y_wide, S_wide = NULL) {
  if (length(model$coef) == 0)
    return(rep(exp(model$intercept), length(zones_all)))
  X <- .feature_matrix(t_for, Y_wide, A_wide, W, G, static, osrm, zones_all,
                       model$cov_spec, Y_active = Y_active, S_wide = S_wide)
  eta <- rep(model$intercept, length(zones_all))
  for (c in model$cov_spec) {
    xc <- (X[, c] - model$center[[c]]) / model$scale[[c]]
    xc[!is.finite(xc)] <- 0
    eta <- eta + model$coef[[c]] * xc
  }
  pmin(exp(eta), 1e6)
}

# ---------------------------------------------------------------------------
# National renewal R (for forward projection of the import force at h>=2)
# ---------------------------------------------------------------------------

estimate_R_local <- function(Y_wide, G_weekly, t_idx) {
  y <- colSums(Y_wide, na.rm = TRUE)
  nT <- min(t_idx, length(y))
  if (nT < 4L) return(1.5)
  # Use the last STABLE week (nT-1); the final onset week is still truncation-low
  # and would spuriously depress R. Average the two most recent stable transitions.
  Rs <- c()
  for (ref in c(nT - 1L, nT - 2L)) {
    if (ref < 2L) next
    Kk  <- min(length(G_weekly), ref - 1L)
    lam <- sum(G_weekly[seq_len(Kk)] * y[ref - seq_len(Kk)])
    if (lam > 0) Rs <- c(Rs, y[ref] / lam)
  }
  R <- if (length(Rs) > 0) mean(Rs) else 1.5
  min(max(R, 0.5), 5)
}

# ---------------------------------------------------------------------------
# The workhorse forecast
# ---------------------------------------------------------------------------

#' Forecast spatial invasion with the mobility-informed renewal workhorse.
#'
#' @param zone_week_nc nowcast-corrected zone-week tibble (health_zone,
#'   week_start, confirmed, confirmed_nc).
#' @param W            outflow mobility matrix (row-stochastic).
#' @param gt_pmf_daily daily GT PMF.
#' @param zones_all    canonical zone vector.
#' @param t_idx        week index of the training cutoff.
#' @param horizons     forecast horizons (weeks).
#' @param obs          "poisson" or "negbin".
#' @param rho          ascertainment (for the infection-scale score).
#' @param mobility_id, gt_profile labels.
#' @param training_cutoff Date of the cutoff week.
#' @return tibble: health_zone, horizon, mu_case, p_case_invasion,
#'   p_infection_invasion, was_active_before, method (+ q05..q95 count quantiles).
forecast_workhorse <- function(zone_week_nc, W, gt_pmf_daily, zones_all,
                               t_idx, horizons = LFO_HORIZONS,
                               obs = "poisson", rho = ASCERTAINMENT_NOMINAL,
                               mobility_id = MOBILITY_PRIMARY,
                               gt_profile = GT_PRIMARY, training_cutoff = NULL,
                               method_label = "Renewal",
                               covariate_spec = NULL,          # Q2 covariate options
                               nowcast_mode = "corrected",     # Q6: corrected | raw
                               alert_col = "total_alerts",
                               static_cov = NULL, osrm = NULL,
                               report_rate = NULL,             # Q6: reporting-rate
                               beta_weighting = "none",        # Q1: none|completeness
                               susp_col = "suspected_nc") {    # suspected-case covariate source
  weeks   <- sort(unique(zone_week_nc$week_start))
  cutoff  <- weeks[min(t_idx, length(weeks))]
  # Q6 nowcast/delay option: use the nowcast-corrected counts (default) or the
  # raw right-censored counts, so the value of nowcasting can be evaluated.
  count_col <- if (identical(nowcast_mode, "raw")) "confirmed" else "confirmed_nc"
  Y_wide  <- .count_wide(zone_week_nc, zones_all, count_col)
  G       <- daily_to_weekly_gt(gt_pmf_daily)
  tmax    <- min(t_idx, ncol(Y_wide))

  # Q6 reporting-rate structure (OPTIONAL): inflate the importation pressure from
  # under-reporting SOURCE zones. Since Lambda_i = sum_j W[j,i]·(G⊛Y_j), dividing
  # source row j of W by its relative reporting rate r_j is exactly equivalent to
  # using true incidence Y_j/r_j in the import force, while leaving the observed
  # Y_wide (used for the at-risk mask and the invasion OUTCOME) untouched. When
  # report_rate is NULL this is a no-op (W_use == W).
  r_rep <- .reporting_rate_vec(report_rate, zones_all)
  W_use <- W
  if (!is.null(r_rep)) {
    rr <- r_rep[rownames(W)]; rr[is.na(rr)] <- 1   # align to W's actual row order
    W_use <- W / rr                                # divides source row j by r_j
  }

  aff_vec <- rowSums(Y_wide[, seq_len(tmax), drop = FALSE] > 0) > 0
  names(aff_vec) <- zones_all

  # Alert (leading-indicator) matrix + static covariate features, only if needed.
  needs_alert  <- any(c("alert_import", "alert_local") %in% covariate_spec)
  A_train <- if (needs_alert)
    .count_wide(zone_week_nc, zones_all, alert_col)[, seq_len(tmax), drop = FALSE] else NULL
  # Forecast-week alert matrix: alerts are not projected, so pad with zero columns
  # for the future weeks (compute_foi indexes columns up to t_for-1).
  A_fore  <- if (!is.null(A_train))
    cbind(A_train, matrix(0, nrow(A_train), max(horizons),
                          dimnames = list(rownames(A_train), NULL))) else NULL
  # Suspected-but-not-confirmed leading-indicator matrix (built like the alert matrix,
  # only when a susp_* covariate is requested): the nowcast-corrected suspected series
  # (susp_col; falls back to zeros if absent). Not projected forward, so pad the forecast
  # weeks with zeros exactly as the alert matrix is, keeping compute_foi in-bounds.
  needs_susp  <- any(c("susp_import", "susp_local") %in% covariate_spec)
  S_train <- if (needs_susp)
    .susp_wide(zone_week_nc, zones_all, prefer = susp_col)[, seq_len(tmax), drop = FALSE] else NULL
  S_fore  <- if (!is.null(S_train))
    cbind(S_train, matrix(0, nrow(S_train), max(horizons),
                          dimnames = list(rownames(S_train), NULL))) else NULL
  static  <- if (any(c("log_pop", "healthsite_density", "ccvi", "positivity") %in%
                     covariate_spec)) .static_features(static_cov, zones_all) else NULL

  # Q1 completeness weighting (OPTIONAL): weight each training week by its
  # national reporting completeness (mean nowcast truncation weight), so recent
  # right-truncated weeks — whose invasion outcomes are least certain — count for
  # less when the import scale beta is calibrated. Default "none" => all weights 1.
  week_w <- NULL
  if (identical(beta_weighting, "completeness") &&
      "trunc_weight" %in% names(zone_week_nc)) {
    cw <- zone_week_nc %>%
      dplyr::group_by(week_start) %>%
      dplyr::summarise(w = mean(trunc_weight, na.rm = TRUE), .groups = "drop")
    week_w <- cw$w[match(weeks, cw$week_start)][seq_len(tmax)]
    week_w[!is.finite(week_w)] <- 1
    week_w <- pmin(pmax(week_w, 0.05), 1)   # keep strictly positive, cap at 1
  }

  Yw_train <- Y_wide[, seq_len(tmax), drop = FALSE]
  model   <- fit_import_model(Yw_train, W_use, G, zones_all, A_wide = A_train,
                              static = static, osrm = osrm,
                              cov_spec = covariate_spec, week_weights = week_w,
                              S_wide = S_train)
  beta0   <- exp(model$intercept)   # mean import scale (for dispersion + fallback)
  # Pass the TRAINING matrix (no columns beyond the cutoff), so the no-leakage
  # guarantee is structural rather than merely true of the current index arithmetic.
  R_local <- estimate_R_local(Yw_train, G, tmax)

  # NegBin dispersion from TEMPORAL residuals of the import model on training wks
  theta <- NA_real_
  if (obs == "negbin")
    theta <- .workhorse_temporal_dispersion(Yw_train, W_use, G, zones_all, beta0, tmax)

  # Current-week (cutoff week W) expected introductions per zone: the import force
  # INTO week W from observed incidence, times the calibrated beta. This is the
  # week-W analogue of the horizon-week mu below, and it lets a forecast be anchored
  # to the analysis DATE (mid-week) rather than the week boundary — a uniform daily
  # introduction hazard of mu_wk0/7 covers the days of the current week that still
  # lie ahead of the analysis date (see anchor_windows_from_analysis_date,
  # 22_daily_reissue.R). Masked to NA for already-affected zones like every other mu.
  Lambda_wk0 <- compute_foi(Yw_train, W_use, G, t_idx = tmax, zones_all)
  beta_wk0   <- .beta_vector(model, tmax, Yw_train, A_fore, W_use, G, static, osrm,
                             zones_all, Y_active = Yw_train, S_wide = S_fore)
  mu_wk0_vec <- as.numeric(beta_wk0 * Lambda_wk0)
  mu_wk0_vec[aff_vec] <- NA_real_

  # Forward-project incidence to build the import force at each horizon.
  Yw <- Yw_train
  mu_cum <- setNames(numeric(length(zones_all)), zones_all)  # cumulative introductions
  results <- vector("list", length(horizons))
  hmax <- max(horizons)
  mu_by_h <- list()
  for (h in seq_len(hmax)) {
    t_for <- ncol(Yw) + 1L
    Lambda_h <- compute_foi(Yw, W_use, G, t_idx = t_for, zones_all)
    # Zone-specific import scale from the (optionally covariate-augmented) model.
    # d_min uses the OBSERVED frontier (Yw_train), fixed over the forecast window,
    # not the projected incidence in the growing Yw (which would spuriously mark
    # fractional-risk frontier zones as already "active" from horizon 2 on).
    beta_h   <- .beta_vector(model, t_for, Yw, A_fore, W_use, G, static, osrm,
                             zones_all, Y_active = Yw_train, S_wide = S_fore)
    mu_h     <- beta_h * Lambda_h               # expected introductions this week
    mu_cum   <- mu_cum + mu_h
    mu_by_h[[h]] <- mu_cum
    # Project incidence forward for the NEXT horizon: local renewal + new imports.
    own_foi <- .gweighted_own(Yw, G, t_for)
    Y_next  <- R_local * own_foi + mu_h
    Yw <- cbind(Yw, Y_next)
  }

  for (hi in seq_along(horizons)) {
    h <- horizons[hi]
    mu_case <- mu_by_h[[h]]                      # cumulative expected introductions
    # PRIMARY forecast p_case is ascertainment-AGNOSTIC: it is P(first confirmed-case
    # onset), estimated from observed invasions, and does NOT divide by rho.
    p_case  <- .invasion_prob(mu_case, obs, theta)
    # Secondary infection-scale score: rho is UNCERTAIN, so report a BAND across the
    # ascertainment grid rather than one fixed value (Task 8). Lower rho => more
    # unobserved infections => higher implied infection risk, so rho=max(grid) gives
    # the lower band and rho=min(grid) the upper band; p_infection at nominal rho.
    rho_lo <- max(ASCERTAINMENT_GRID); rho_hi <- min(ASCERTAINMENT_GRID)
    p_inf     <- .invasion_prob(mu_case / max(rho,    1e-3), obs, theta)
    p_inf_lo  <- .invasion_prob(mu_case / max(rho_lo, 1e-3), obs, theta)
    p_inf_hi  <- .invasion_prob(mu_case / max(rho_hi, 1e-3), obs, theta)
    # Mask affected zones -> NA (no invasion probability for already-affected zones)
    p_case[aff_vec]  <- NA_real_
    p_inf[aff_vec]   <- NA_real_
    p_inf_lo[aff_vec] <- NA_real_
    p_inf_hi[aff_vec] <- NA_real_
    mu_case[aff_vec] <- NA_real_
    # Count predictive quantiles, matched to the observation model (secondary
    # count task); NegBin tail when obs="negbin" so they agree with p_case.
    # mu_case is NA for affected zones, so mu_q (and every quantile) stays NA there
    # — the "no invasion output for already-affected zones" invariant holds fully.
    mu_q <- pmax(mu_case, 1e-6)
    results[[hi]] <- tibble::tibble(
      health_zone          = zones_all,
      horizon              = h,
      mu_forecast          = mu_case,
      mu_wk0               = mu_wk0_vec,          # current-week introductions (analysis-date anchoring)
      p_invasion           = p_case,             # kept for back-compat = p_case
      p_case_invasion      = p_case,
      p_infection_invasion = p_inf,
      p_infection_lo       = p_inf_lo,   # ascertainment band (rho = max grid)
      p_infection_hi       = p_inf_hi,   # ascertainment band (rho = min grid)
      was_active_before    = unname(aff_vec),
      q05 = .qcount(0.05, mu_q, obs, theta), q20 = .qcount(0.20, mu_q, obs, theta),
      q25 = .qcount(0.25, mu_q, obs, theta), q75 = .qcount(0.75, mu_q, obs, theta),
      q80 = .qcount(0.80, mu_q, obs, theta), q95 = .qcount(0.95, mu_q, obs, theta),
      method               = method_label,
      mobility_id          = mobility_id,
      gt_profile           = gt_profile,
      training_cutoff      = training_cutoff %||% cutoff
    )
  }
  out <- dplyr::bind_rows(results)
  attr(out, "beta") <- beta0
  attr(out, "R_local") <- R_local
  attr(out, "n_events") <- model$n_events
  attr(out, "coef_cov") <- model$coef
  out
}

.gweighted_own <- function(Yw, G, t_for) {
  n <- nrow(Yw); own <- numeric(n)
  for (k in seq_along(G)) {
    tp <- t_for - k
    if (tp < 1) break
    own <- own + G[k] * Yw[, tp]
  }
  own
}

.invasion_prob <- function(mu, obs, theta) {
  mu <- pmax(mu, 0)
  if (obs == "negbin" && is.finite(theta) && theta > 0) {
    1 - (theta / (theta + mu))^theta
  } else {
    1 - exp(-mu)
  }
}

# Count predictive quantile consistent with the observation model, so the
# reported q05..q95 imply the SAME P(Y>=1) as p_case_invasion (NegBin tail when
# obs="negbin", Poisson otherwise) rather than always using the Poisson tail.
.qcount <- function(pr, mu, obs, theta) {
  if (obs == "negbin" && is.finite(theta) && theta > 0)
    stats::qnbinom(pr, size = theta, mu = mu)
  else
    stats::qpois(pr, mu)
}

# Dispersion from temporal residuals of the import model (mu vs realised new
# cases) across training weeks, on at-risk zones with positive import force.
.workhorse_temporal_dispersion <- function(Y_wide, W, G, zones_all, beta, t_max) {
  res <- list()
  for (t in seq_len(t_max - 1L)) {
    Lam <- compute_foi(Y_wide, W, G, t_idx = t + 1L, zones_all)
    aff <- rowSums(Y_wide[, seq_len(t), drop = FALSE] > 0) > 0
    keep <- !aff & Lam > 0
    if (!any(keep)) next
    res[[length(res) + 1L]] <- tibble::tibble(mu = beta * Lam[keep],
                                              y = Y_wide[keep, t + 1L])
  }
  d <- dplyr::bind_rows(res)
  if (nrow(d) < 10L) return(10)
  # Method-of-moments NB size. With heterogeneous means, E[(y-mu)^2] = E[mu] +
  # E[mu^2]/theta, so the excess residual variance over the mean estimates
  # E[mu^2]/theta and theta = E[mu^2] / excess. Using mean(mu)^2 here (instead of
  # E[mu^2] = mean(mu^2)) would understate theta whenever mu varies across zones.
  # MoM needs the mean SQUARED residual E[(y-mu)^2], not the central sample variance
  # var(y-mu) (which subtracts the mean residual and understates the second moment when
  # the forecast is biased); they coincide only when mean(y-mu)=0.
  v <- mean((d$y - d$mu)^2, na.rm = TRUE); mbar <- mean(d$mu, na.rm = TRUE)
  extra <- v - mbar
  if (!is.finite(extra) || extra <= 0) return(50)  # ~Poisson
  max(mean(d$mu^2, na.rm = TRUE) / extra, 0.1)
}

# ---------------------------------------------------------------------------
# Risk scores (at-risk zones only)
# ---------------------------------------------------------------------------

#' Compute the four relative/absolute risk scores for a workhorse forecast.
#'
#' @param fc         forecast tibble from forecast_workhorse() (one method).
#' @param province_map tibble(nom, province) from load_province_map().
#' @return fc augmented with province, rr_ituri, rr_ituri_rank, rr_nat,
#'   rr_nat_rank (all NA for affected zones).
#' Column-safe province suffix: "Nord-Kivu" -> "nordkivu", "Haut-Uele" -> "hautuele".
.prov_suffix <- function(p) tolower(gsub("[^A-Za-z]", "", p))

compute_risk_scores <- function(fc, province_map,
                                provinces = get0("PROVINCES_OF_INTEREST",
                                                 ifnotfound = c("Ituri", "Nord-Kivu", "Haut-Uele"))) {
  if (!is.null(province_map)) {
    fc <- fc %>% dplyr::left_join(province_map, by = c("health_zone" = "nom"))
  } else {
    fc$province <- NA_character_
  }
  grp <- intersect(c("method", "horizon", "training_cutoff"), names(fc))
  fc <- fc %>%
    dplyr::mutate(.mu = ifelse(was_active_before, NA_real_, mu_forecast)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
    dplyr::mutate(
      # nationwide relative risk (multiplier vs mean at-risk zone) + share + rank
      rr_nat        = .mu / mean(.mu, na.rm = TRUE),
      rr_nat_share  = .mu / sum(.mu, na.rm = TRUE),
      rr_nat_rank   = dplyr::min_rank(dplyr::desc(.mu))
    ) %>%
    dplyr::ungroup()
  # Within-province relative risk for each province of interest (Ituri, Nord-Kivu,
  # Haut-Uele, ...): each zone's expected introductions relative to the mean/sum of
  # the still-at-risk zones IN THAT PROVINCE, with an in-province rank. A zone can
  # top its province yet be modest nationally, and vice-versa. rr_ituri is retained
  # (backward-compatible) as the Ituri entry of this generalised set.
  for (prov in provinces) {
    s <- .prov_suffix(prov)
    fc <- fc %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grp)), .inp = (province == prov)) %>%
      dplyr::mutate(
        !!paste0("rr_", s)           := ifelse(.inp, .mu / mean(.mu, na.rm = TRUE), NA_real_),
        !!paste0("rr_", s, "_share") := ifelse(.inp, .mu / sum(.mu, na.rm = TRUE), NA_real_),
        !!paste0("rr_", s, "_rank")  := ifelse(.inp, dplyr::min_rank(dplyr::desc(.mu)), NA_integer_)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.inp)
  }
  fc %>% dplyr::select(-.mu)
}

message("[workhorse] 15_workhorse.R loaded — mobility-informed renewal invasion model.")
