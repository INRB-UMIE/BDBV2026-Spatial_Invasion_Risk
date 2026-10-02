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
# Risk scores (at-risk zones only): absolute p_case on the CONFIRMED-case scale,
# relative risk within each province of interest (Ituri, Nord-Kivu, Haut-Uele) and
# relative risk nationwide. There is no ascertainment adjustment: see 00_config.R.
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
  # TIME TREND on the import->invasion conversion (the `tv-trend` models). The week index
  # itself is the covariate, so beta_t = exp(beta0 + gamma_w * z(week)) — a log-linear trend
  # on the hazard. Deliberately implemented as an ORDINARY COVARIATE rather than as a new
  # model class: build_invasion_design() then standardises it, .bayes_invasion_draws()
  # re-standardises it on the forecast rows with the SAME center/scale, and the LFO closure,
  # the stacking alignment and the hazard-ratio tables all work unchanged. It is constant
  # across zones within a week (a pure time effect) and, being evaluated at t_for, it takes
  # the FORECAST week's value on forecast rows — i.e. the trend extrapolates, which is the
  # whole point. Note the offset log(Lambda) has a fixed coefficient of 1, so this term is
  # not collinear with it: it absorbs systematic drift in how Lambda converts to invasions.
  if ("week_idx" %in% cov_spec) f$week_idx <- rep(as.numeric(t_for), length(zones_all))
  if (length(f) == 0) return(NULL)
  m <- do.call(cbind, f); colnames(m) <- names(f); m
}



# ---------------------------------------------------------------------------
# National renewal R (for forward projection of the import force at h>=2)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# The workhorse forecast
# ---------------------------------------------------------------------------


.gweighted_own <- function(Yw, G, t_for) {
  n <- nrow(Yw); own <- numeric(n)
  for (k in seq_along(G)) {
    tp <- t_for - k
    if (tp < 1) break
    own <- own + G[k] * Yw[, tp]
  }
  own
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
