# =============================================================================
# 43_spread_kinematics.R — DESCRIPTIVE SPATIOTEMPORAL KINEMATICS of the BDBV 2026
# outbreak: how fast, in which direction, and how far from the epicentre the
# epidemic has moved.
#
# Everything here is DESCRIPTIVE (no model, no forecast): a small set of classical
# spatial-statistics summaries computed on the onset-dated confirmed line list.
#
#   (1) CENTRE OF GRAVITY (mean centre)
#         - of WEEKLY CASES  : case-weighted mean centre of that week's confirmed cases
#         - of NEW INVASIONS : unweighted mean centre of the zones invaded that week
#       with, per week, the STANDARD DISTANCE (weighted RMS distance to the mean
#       centre = spatial dispersion) and the standard deviational ELLIPSE
#       (semi-major/semi-minor axis, orientation) that carries the anisotropy.
#
#   (2) SPEED OF SPREAD — complementary estimators, all in km/week unless stated.
#       They deliberately do NOT agree: each answers a different question, and the
#       gap between them IS the result (long-range seeding, not a travelling wave).
#         v_radial : OLS slope of zone distance-from-epicentre on invasion week
#                    (mean radial expansion rate of the invaded set; also fitted on
#                    great-circle km and on travel-time hours)
#         v_front  : the LEADING EDGE R_max(t) = furthest zone invaded by week t
#                    (whole-window OLS slope, for reference only), and the ENVELOPE
#                    R_p90(t) = 90th-percentile reach, split into the slope while it was
#                    ESTABLISHED and the slope AFTER it first reached 90% of its final
#                    value; the saturation verdict is derived from those two slopes
#         v_reach  : OLS slope of the MEDIAN distance of invaded zones on week, and of
#                    the case-weighted mean distance (bulk, not envelope, expansion)
#         v_cog    : step speed of the weekly-case centre of gravity
#                    (median of |CoG_t - CoG_{t-1}| / Dt), plus the NET speed
#                    |CoG_T - CoG_0| / T and the straightness of the CoG track
#         v_disp   : OLS slope of the standard distance on week (dispersion growth);
#                    the companion diffusion coefficient D = slope(SD^2 ~ t)/4 km^2/wk
#         accrual  : OLS slope of log(cumulative invaded zones) on week -> doubling time
#
#   (3) REACH FROM THE EPICENTRE — per week, over the zones invaded SO FAR:
#       max / 90th-percentile / median ROAD distance from the index zone, the
#       case-weighted mean road distance of that week's cases, the number of zones and
#       provinces reached, and the share of weekly cases arising more than
#       REACH_FAR_KM by road from the index zone.
#       PLUS, and not optional: how far THAT WEEK's new invasions actually landed
#       (new_max_km / new_med_km). The cumulative maxima above are running maxima, so
#       once one zone sets a far record every later introduction is invisible in them
#       however distant — on the 2026-10-09 frame a 767 km seeding in week 13 hides
#       completely behind a 794 km one from week 3, in a different direction.
#       Reporting only the running max would say "the extent saturated" and let a
#       reader conclude that long-range spread had stopped, when it demonstrably had not.
#
# Distance is measured three ways for every zone — great-circle (km), OSRM road
# distance (km) and OSRM travel time (h) — and the three are compared on how well
# they order the observed arrival times (Spearman rho / OLS R^2). Road distance is
# the headline metric; great-circle is reported alongside.
#
# CARE TAKEN (the things that would otherwise make these numbers wrong):
#   * RIGHT TRUNCATION. Onset-dated recent weeks are incomplete. The number of
#     trailing weeks to drop is chosen from the DATA: the empirical onset->
#     confirmation delay ECDF, estimated on the sub-sample with >= DELAY_CAP_DAYS of
#     follow-up (so the ECDF itself is not truncated), gives an expected reporting
#     completeness per onset week; trailing weeks below COMPLETE_MIN are excluded
#     from EVERY fit — the kinematics AND the arrival-predictor fits of Section 5b —
#     and drawn hollow in Figure_spread_kinematics_compact. Sensitivity to that choice
#     (1/2/3 weeks) is reported in the diagnostics CSV.
#   * ONSET IMPUTATION. A share of confirmed records have no usable onset date and are
#     dated by load_linelist() from their sample date minus a donor interval (00_config.R
#     ONSET_MODE). The share is NOT typed here: load_linelist() publishes it as
#     .PCT_ONSET_IMPUTED and this script prints it at run time (a "~15%" written here was
#     ~24.5% on the frame it described). Each record's draw is a hash of (alert_id,
#     ONSET_IMPUTE_SEED), with no RNG stream, so a local set.seed() is a no-op. The primary
#     run uses the project's own ONSET_IMPUTE_SEED (default RANDOM_SEED), so its invasion
#     dates are exactly the ones the rest of the suite scores; N_IMPUTE replicates then VARY
#     that onset seed alone and re-run load_linelist() and compute_kinematics() on each
#     completion, so every summary-CSV statistic carries an imputation interval (imp_*
#     columns) except the centre-of-gravity BEARING, an angle, for which linear percentiles
#     are not meaningful; the envelope verdict's stability across completions is stated on
#     its row, and the arrival-predictor fits carry their own interval columns. The
#     headline VALUE is that primary completion, not a summary over completions (see the
#     imp_* columns of the summary CSV).
#   * PROJECTION. Mean centres are computed in a Lambert azimuthal EQUAL-AREA plane
#     centred on the outbreak; every reported DISTANCE is geodesic (WGS84) or the
#     OSRM road network — never a planar approximation.
#   * ASYMMETRIC / MISSING ROAD ROUTES. The OSRM matrices are near- but not exactly
#     symmetric (one-way and ferry routing); they are symmetrised as the mean of the
#     available directions. Two zones (Bokoro, Idjwi) are unroutable and stay NA.
#   * ZERO-EVENT WEEKS. The new-invasion centre of gravity is undefined in weeks
#     with no new invasion; step speeds divide by the ACTUAL week gap, not by 1.
#
# CAVEATS stated on the products themselves:
#   - Cases are located at zone CENTROIDS, so every dispersion measure is
#     BETWEEN-zone only and understates the true spatial spread.
#   - R_max(t) is a running maximum: monotone by construction and serially
#     dependent, so its OLS slope is a descriptive summary, not an estimate with a
#     valid confidence interval. The inferential speed is v_radial (zone bootstrap).
#     A FLAT R_max means no introduction beat the standing record — it does NOT mean
#     introductions stopped, and must never be read that way. Always read it against
#     new_max_km, which is the per-week series that shows what actually arrived.
#   - Zones whose confirmed count exists only through the INSP sitrep top-up
#     (01_data_prep.R .build_sitrep_confirmed_appends) are dated from the sitrep
#     REPORTING date minus a delay draw, so their invasion week is late-biased;
#     they are flagged in the per-zone CSV.
#   - The new-invasion centre of gravity is a mean over however many zones were invaded
#     that week — often just one or two. Its weekly "speed" is then the gap between two
#     individual zones rather than a centre-of-mass velocity, which is precisely why it
#     is one to two orders of magnitude larger than the case centre's. new_zones is in
#     the weekly CSV so any week's n can be checked.
#   - The earliest weeks rest on very few cases (week 0 carries single digits), so the
#     first centre of gravity — and therefore the NET centre-of-gravity speed, which is
#     anchored on it — is the noisiest quantity here. The onset-draw interval on it is
#     correspondingly wide; read the MEDIAN step speed as the stable summary.
#   - The completeness ECDF is fitted on cases with >= DELAY_CAP_DAYS of follow-up, i.e.
#     EARLY-epidemic cases. If notification delays have since shortened it over-drops
#     recent weeks (and vice versa); the truncation sensitivity table bounds that.
#   - Nord-Kivu was invaded within the first epidemic weeks, far from the index zone (on the
#     2026-10-09 frame four Nord-Kivu zones 190-280 km by road away in week 1, Goma at 550 km
#     in week 2; invasion_week in the zones CSV), so the invaded set is NOT a single
#     expanding wave; the radial statistics summarise a multi-focus process and should be
#     read as such.
#
# Run:  Rscript spatiotemporal/43_spread_kinematics.R
# Env:  SPREAD_SEED, SPREAD_TRUNC_WEEKS, SPREAD_IMPUTE_REPS, SPREAD_BOOT
# Out:  <OUT_DIR>/key_outputs/spread_kinematics_{weekly,zones,summary,diagnostics}.csv
#       <OUT_DIR>/key_outputs/arrival_predictors.csv, arrival_predictor_fits.csv
#       (no figure: make_spread_kinematics_compact.R and make_manuscript_figures.R draw
#       these tables; see SECTION 6)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(sf); library(geosphere); library(here)
})

ST_DIR <- Sys.getenv("SPREAD_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
source(file.path(ST_DIR, "00_config.R"))
source(file.path(ST_DIR, "01_data_prep.R"))
# Mobility loaders + name harmonisation, for the arrival-predictor table in Section 5b.
# 03 defines functions only (its sole top-level statement re-sources 00_config.R).
source(file.path(ST_DIR, "03_mobility_matrices.R"))

.envint <- function(nm, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(nm, unset = NA_character_)))
  if (is.na(v)) default else v
}

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
# Seeds the zone bootstrap only. The onset draw is NOT seeded from here — it comes
# from the project-wide ONSET_IMPUTE_SEED that load_linelist() reads (see Section 4).
SEED           <- .envint("SPREAD_SEED", 20260812L)
SEED_WEEKS     <- 3L        # weeks used to identify the index (epicentre) zone
DELAY_CAP_DAYS <- 30L       # follow-up needed for an untruncated delay ECDF (and its cap)
COMPLETE_MIN   <- 0.90      # weeks below this expected completeness are treated as truncated
TRUNC_MIN      <- 1L        # never trust the final week, whatever the delay model says
TRUNC_MAX      <- 3L        # never drop more than this (guards a pathological delay fit)
TRUNC_FORCE    <- .envint("SPREAD_TRUNC_WEEKS", NA_integer_)   # manual override
N_BOOT         <- .envint("SPREAD_BOOT", 2000L)       # zone bootstrap replicates
N_IMPUTE       <- .envint("SPREAD_IMPUTE_REPS", 100L) # onset-imputation replicates
REACH_FAR_KM   <- 100       # "far from the epicentre" threshold for the case-share statistic
TRUNC_SENS     <- 1:3       # trailing-week drops compared in the sensitivity table

KEY_DIR <- file.path(OUT_DIR, "key_outputs")
if (!dir.exists(KEY_DIR)) dir.create(KEY_DIR, recursive = TRUE, showWarnings = FALSE)



# =============================================================================
# SECTION 1 — geography: centroids, provinces, distance metrics
# =============================================================================

#' Symmetrise an OSRM matrix: mean of the two available directions; a pair routable
#' in only one direction keeps that direction; an unroutable pair stays NA.
.symmetrise <- function(M) {
  Mt <- t(M)
  out <- (M + Mt) / 2
  i <- is.na(out) & !is.na(M);  out[i] <- M[i]
  j <- is.na(out) & !is.na(Mt); out[j] <- Mt[j]
  out
}

.load_osrm_named <- function(kind) {
  f <- file.path(OSRM_DIR, sprintf("osrm__%s__static.matrix.csv", kind))
  stopifnot(file.exists(f))
  df <- readr::read_csv(f, col_types = readr::cols(.default = "d", nom = "c"), show_col_types = FALSE)
  M <- as.matrix(df[, setdiff(colnames(df), "nom")]); rownames(M) <- df$nom
  stopifnot(nrow(M) == ncol(M), identical(rownames(M), colnames(M)))
  .symmetrise(M)
}

#' Lambert azimuthal equal-area PROJ string centred on (lon0, lat0).
.laea <- function(lon0, lat0)
  sprintf("+proj=laea +lat_0=%.6f +lon_0=%.6f +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
          lat0, lon0)

#' Build the geography once: WGS84 + equal-area centroids, provinces, road/time matrices.
#' `proj_centre` (lon, lat) fixes the equal-area plane used for the mean centres; when
#' NULL the plane is centred on the middle of the national bounding box.
build_geo <- function(proj_centre = NULL) {
  shp <- sf::st_read(SHAPEFILE_PATH, quiet = TRUE)
  stopifnot(all(c("Nom", "PROVINCE") %in% names(shp)))
  # The shapefile Nom is a 1:1 match to the WorldPop modelling spine, so a name join
  # is unambiguous here — no province-composite key is needed (asserted below).
  stopifnot(anyDuplicated(shp$Nom) == 0L)

  # Pass 1: rough WGS84 centroids, only to locate the projection plane.
  c0 <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(shp))))
  if (is.null(proj_centre)) proj_centre <- c(mean(range(c0[, 1])), mean(range(c0[, 2])))

  # Pass 2: equal-area centroids on a plane centred where the outbreak actually is,
  # then transformed back to WGS84 so every distance can be computed geodesically.
  crs_laea <- .laea(proj_centre[1], proj_centre[2])
  ctr_m  <- sf::st_centroid(sf::st_geometry(sf::st_transform(shp, crs_laea)))
  xy_m   <- sf::st_coordinates(ctr_m)
  xy_ll  <- sf::st_coordinates(sf::st_transform(ctr_m, 4326))
  stopifnot(nrow(xy_m) == nrow(shp), nrow(xy_ll) == nrow(shp))

  ctr <- tibble::tibble(health_zone = as.character(shp$Nom),
                        province = as.character(shp$PROVINCE),
                        x = xy_m[, 1], y = xy_m[, 2],
                        lon = xy_ll[, 1], lat = xy_ll[, 2])

  road <- .load_osrm_named("road_distance")            # km
  time <- .load_osrm_named("travel_time") / 60         # min -> h
  stopifnot(setequal(rownames(road), ctr$health_zone),
            identical(dimnames(road), dimnames(time)))

  list(shp = shp, ctr = ctr, road = road, time = time,
       crs_laea = crs_laea, proj_centre = proj_centre)
}

# =============================================================================
# SECTION 2 — right-truncation: expected reporting completeness per onset week
# =============================================================================

#' Expected share of an onset week's confirmed cases already observed as of `asof`,
#' from the empirical onset -> confirmation delay ECDF.
#'
#' The ECDF is estimated ONLY on records whose onset is at least `cap` days before
#' `asof` — that sub-sample has complete follow-up out to `cap` days, so the ECDF is
#' not itself right-truncated (the naive all-records ECDF would over-state
#' completeness). Delays are CAPPED at `cap` days (pmin), not dropped: a longer delay
#' still counts in the denominator, as not yet arrived at any follow-up below `cap`, and
#' is treated as arrived at `cap`, so this is a slight UPPER bound on completeness only
#' at `cap` itself — conservative about how few weeks it flags, which is why TRUNC_MIN
#' also applies. (Dropping the longer delays instead would condition the ECDF on
#' D <= cap and over-state completeness at every follow-up.) The share of delays beyond
#' `cap` is returned as attr "share_over_cap" and printed by the caller, so how rare they
#' are is measured on the frame rather than asserted here. "Known" is the lab date, which
#' load_linelist() has already restored to the CONFIRMING test's (apply_earliest_sample_dates(),
#' 00_config.R), so a later repeat test does not lengthen the delay. Sitrep top-up rows
#' carry no lab/report/notification date and drop out of the ECDF automatically.
week_completeness <- function(ll, weeks, asof, cap = DELAY_CAP_DAYS) {
  cf <- ll %>% dplyr::filter(confirmed %in% TRUE, !is.na(date_index))
  known <- dplyr::coalesce(as.Date(cf$lab_analysis_date),
                           as.Date(cf$reporting_date),
                           as.Date(cf$date_of_notification))
  del <- as.numeric(known - cf$date_index)
  ok  <- is.finite(del) & del >= 0 & cf$date_index <= (asof - cap)
  if (sum(ok) < 50L) {
    warning("[kinematics] only ", sum(ok), " complete-follow-up delay pairs — ",
            "completeness not estimated; falling back to TRUNC_MIN.", call. = FALSE)
    return(rep(NA_real_, length(weeks)))
  }
  Fhat <- stats::ecdf(pmin(del[ok], cap))
  out <- vapply(weeks, function(w) {
    avail <- as.numeric(asof - (w + 0:6))              # follow-up days per onset day
    mean(ifelse(avail < 0, 0, Fhat(pmin(avail, cap)))) # onset days after `asof` contribute 0
  }, numeric(1))
  attr(out, "share_over_cap") <- mean(del[ok] > cap)
  attr(out, "n_pairs")        <- sum(ok)
  out
}

#' Number of trailing weeks to exclude, from the completeness profile.
choose_trunc <- function(comp, force = NA_integer_) {
  if (!is.na(force)) return(max(force, 0L))
  if (all(is.na(comp))) return(TRUNC_MIN)
  k <- 0L; i <- length(comp)
  while (i >= 1L && !is.na(comp[i]) && comp[i] < COMPLETE_MIN) { k <- k + 1L; i <- i - 1L }
  min(max(k, TRUNC_MIN), TRUNC_MAX)
}

# =============================================================================
# SECTION 3 — kinematics
# =============================================================================

#' Weighted mean centre, standard distance and standard deviational ellipse of a
#' point set with weights w, in the equal-area plane (metres in, km out).
#' Returns NA-filled when there is no positive weight.
.mean_centre <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  if (!any(ok))
    return(list(x = NA_real_, y = NA_real_, sd_km = NA_real_, major_km = NA_real_,
                minor_km = NA_real_, orient_deg = NA_real_, n = 0L))
  x <- x[ok]; y <- y[ok]; w <- w[ok]; sw <- sum(w)
  mx <- sum(w * x) / sw; my <- sum(w * y) / sw
  sxx <- sum(w * (x - mx)^2) / sw
  syy <- sum(w * (y - my)^2) / sw
  sxy <- sum(w * (x - mx) * (y - my)) / sw
  ev  <- eigen(matrix(c(sxx, sxy, sxy, syy), 2, 2), symmetric = TRUE)
  lam <- pmax(ev$values, 0)                       # guard tiny negative eigenvalues
  v1  <- ev$vectors[, 1]
  # Orientation as a compass bearing of the major axis, folded to [0, 180).
  orient <- (90 - atan2(v1[2], v1[1]) * 180 / pi) %% 180
  list(x = mx, y = my,
       sd_km      = sqrt(sxx + syy) / 1000,
       major_km   = sqrt(lam[1]) / 1000,
       minor_km   = sqrt(lam[2]) / 1000,
       orient_deg = orient, n = length(x))
}

#' Geodesic distance (km) between two lon/lat points; NA-safe.
.geo_km <- function(p1, p2) {
  if (any(!is.finite(c(p1, p2)))) return(NA_real_)
  geosphere::distGeo(p1, p2) / 1000
}

#' Compass bearing (deg, [0,360)) between two lon/lat points; NA-safe.
.geo_brg <- function(p1, p2) {
  if (any(!is.finite(c(p1, p2)))) return(NA_real_)
  b <- geosphere::bearing(p1, p2)
  if (!is.finite(b)) NA_real_ else (b + 360) %% 360
}

#' OLS slope / intercept / R^2 of y on x; NA when under-determined.
.ols <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 2L)
    return(list(slope = NA_real_, intercept = NA_real_, r2 = NA_real_, n = sum(ok)))
  # A constant response (e.g. a saturated running maximum) has zero slope by
  # inspection; lm() would return the same but warn about a "perfect fit", and R^2
  # is genuinely undefined when the total sum of squares is zero.
  if (stats::var(y[ok]) == 0)
    return(list(slope = 0, intercept = y[ok][1], r2 = NA_real_, n = sum(ok)))
  fit <- stats::lm(y[ok] ~ x[ok])
  list(slope = unname(stats::coef(fit)[2]), intercept = unname(stats::coef(fit)[1]),
       r2 = summary(fit)$r.squared, n = sum(ok))
}

#' Theil-Sen slope: the median of all pairwise slopes. Reported beside every OLS
#' radial rate because the invaded set is small (tens of zones) and heavy-tailed —
#' a handful of long-range seedings carry enormous leverage, and OLS alone cannot
#' distinguish "the epidemic expanded at v" from "one distant zone was infected".
#' Agreement between the two is the evidence that the rate is a property of the
#' process; divergence says the OLS slope is an artefact of a few points.
.theil_sen <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 3L) return(NA_real_)
  i <- utils::combn(n, 2L)
  dx <- x[i[2, ]] - x[i[1, ]]; dy <- y[i[2, ]] - y[i[1, ]]
  s <- dy[dx != 0] / dx[dx != 0]
  if (!length(s)) NA_real_ else stats::median(s)
}

#' NA-safe quantile: NA (not an error) when nothing is finite.
.q <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, p))
}
.mx <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else max(x) }

#' Displacement / speed / bearing between CONSECUTIVE DEFINED centres, dividing by
#' the actual week gap (the new-invasion centre is undefined in weeks with no new
#' invasion, so gaps are common there).
.steps <- function(lon, lat, wk) {
  n <- length(lon)
  disp <- rep(NA_real_, n); spd <- disp; brg <- disp
  ok <- which(is.finite(lon) & is.finite(lat))
  if (length(ok) >= 2L) for (j in seq_along(ok)[-1]) {
    a <- ok[j - 1L]; b <- ok[j]
    dt <- wk[b] - wk[a]
    disp[b] <- .geo_km(c(lon[a], lat[a]), c(lon[b], lat[b]))
    spd[b]  <- if (dt > 0) disp[b] / dt else NA_real_
    brg[b]  <- .geo_brg(c(lon[a], lat[a]), c(lon[b], lat[b]))
  }
  list(disp = disp, speed = spd, bearing = brg)
}

#' Full kinematic summary for one zone-week case table.
#'
#' Both windows are anchored to CALENDAR DATES, not to positions in this table's week
#' vector, so two tables whose week grids start on different weeks still fit the same span
#' of real time. On the current line list the grid's first week does not move between onset
#' completions — load_linelist() floors date_index at OUTBREAK_START, so no draw can open an
#' earlier week, and aggregate_to_zone_week() ends the grid on the analysis week — so the
#' anchoring is a guard on that invariant, not a correction the replicates currently need.
#'
#' @param zw          zone-week tibble (health_zone, week_start, confirmed) covering
#'                    ALL zones, zero-filled — output of aggregate_to_zone_week().
#' @param geo         build_geo() product.
#' @param fit_cutoff  last week_start (a Date) included in the fits; later weeks are
#'                    the right-truncated tail.
#' @param seed_cutoff last week_start of the window used to pick the index zone.
#' @param epi_zone    optional fixed index zone (else derived from the seed window).
compute_kinematics <- function(zw, geo, fit_cutoff, seed_cutoff, epi_zone = NULL) {

  weeks <- sort(unique(zw$week_start))
  stopifnot(length(weeks) >= 4L, inherits(fit_cutoff, "Date"), inherits(seed_cutoff, "Date"))
  wk_idx <- function(d) as.integer(round(as.numeric(d - weeks[1]) / 7))
  T_last  <- wk_idx(weeks[length(weeks)])
  stopifnot(fit_cutoff >= weeks[1], seed_cutoff >= weeks[1])
  fit_max <- max(wk_idx(weeks[weeks <= fit_cutoff]))   # last week index used in fits
  n_trunc <- T_last - fit_max
  stopifnot(fit_max >= 2L, n_trunc >= 0L)

  zw <- zw %>%
    dplyr::mutate(week = wk_idx(week_start)) %>%
    dplyr::left_join(dplyr::select(geo$ctr, health_zone, province, x, y, lon, lat),
                     by = "health_zone")
  stopifnot(!anyNA(zw$x))                         # every modelled zone must have geometry

  # ---- index zone (epicentre) ---------------------------------------------
  # Most confirmed cases in the seed window; ties broken by the earliest first
  # confirmed week, then alphabetically, so the choice is fully deterministic.
  first_wk <- zw %>% dplyr::filter(confirmed > 0) %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(inv_week = min(week), .groups = "drop")
  seed <- zw %>% dplyr::filter(week_start <= seed_cutoff) %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(seed_cases = sum(confirmed), .groups = "drop") %>%
    dplyr::filter(seed_cases > 0) %>%
    dplyr::left_join(first_wk, by = "health_zone") %>%
    dplyr::arrange(dplyr::desc(seed_cases), inv_week, health_zone)
  stopifnot(nrow(seed) > 0)
  if (is.null(epi_zone)) epi_zone <- seed$health_zone[1]
  stopifnot(epi_zone %in% geo$ctr$health_zone)
  epi <- geo$ctr[match(epi_zone, geo$ctr$health_zone), ]
  epi_ll <- c(epi$lon, epi$lat)
  stopifnot(all(is.finite(epi_ll)))

  # ---- per-zone distances from the index zone ------------------------------
  zc <- geo$ctr %>%
    dplyr::mutate(
      d_gc    = as.numeric(geosphere::distGeo(epi_ll, cbind(lon, lat))) / 1000,
      d_road  = as.numeric(geo$road[epi_zone, health_zone]),
      t_road  = as.numeric(geo$time[epi_zone, health_zone]),
      bearing = (as.numeric(geosphere::bearing(epi_ll, cbind(lon, lat))) + 360) %% 360)
  zc$bearing[zc$health_zone == epi_zone] <- NA_real_   # bearing to self is undefined
  zx <- setNames(zc$x, zc$health_zone); zy <- setNames(zc$y, zc$health_zone)

  # ---- invaded zones -------------------------------------------------------
  zones_inv <- zw %>% dplyr::filter(confirmed > 0) %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(inv_week = min(week), cases_total = sum(confirmed), .groups = "drop") %>%
    dplyr::left_join(dplyr::select(zc, health_zone, province, lon, lat,
                                   d_gc, d_road, t_road, bearing), by = "health_zone") %>%
    dplyr::mutate(inv_week_start = weeks[inv_week + 1L],
                  in_fit_window  = inv_week <= fit_max) %>%
    dplyr::arrange(inv_week, dplyr::desc(cases_total))
  fitz <- zones_inv %>% dplyr::filter(in_fit_window)

  # ---- weekly centre of gravity + reach ------------------------------------
  weekly <- dplyr::bind_rows(lapply(seq_along(weeks), function(k) {
    w <- weeks[k]; ti <- wk_idx(w)
    dw  <- zw %>% dplyr::filter(week_start == w)
    nz  <- zones_inv %>% dplyr::filter(inv_week == ti)   # newly invaded this week
    cum <- zones_inv %>% dplyr::filter(inv_week <= ti)   # invaded so far
    cg  <- .mean_centre(dw$x, dw$y, dw$confirmed)
    ig  <- .mean_centre(zx[nz$health_zone],  zy[nz$health_zone],  rep(1, nrow(nz)))
    cd  <- .mean_centre(zx[cum$health_zone], zy[cum$health_zone], rep(1, nrow(cum)))
    cs  <- dw$confirmed
    drd <- zc$d_road[match(dw$health_zone, zc$health_zone)]
    tot <- sum(cs)
    okr <- is.finite(drd) & cs > 0     # case-weighted mean uses routable zones only
    tibble::tibble(
      week = ti, week_start = w, cases = tot,
      new_zones = nrow(nz), cum_zones = nrow(cum),
      cum_provinces = dplyr::n_distinct(cum$province),
      cog_x = cg$x, cog_y = cg$y, cog_sd_km = cg$sd_km,
      cog_major_km = cg$major_km, cog_minor_km = cg$minor_km,
      cog_orient_deg = cg$orient_deg,
      inv_cog_x = ig$x, inv_cog_y = ig$y,
      invaded_sd_km = cd$sd_km,
      reach_max_km    = .mx(cum$d_road),
      reach_p90_km    = .q(cum$d_road, 0.90),
      reach_med_km    = .q(cum$d_road, 0.50),
      reach_max_gc_km = .mx(cum$d_gc),
      # How far THIS WEEK's new invasions actually landed. The running maxima above
      # are cumulative, so once one zone sets a far record every later introduction is
      # invisible in them however distant it is — on the 2026-10-09 frame a 767 km seeding in
      # week 13 hides completely behind a 794 km one from week 3, in another direction.
      # These two columns are the un-hidden version, and are what shows that long-range
      # seeding continued long after the envelope stopped growing.
      new_max_km      = .mx(nz$d_road),
      new_med_km      = .q(nz$d_road, 0.50),
      case_mean_km    = if (any(okr)) sum(cs[okr] * drd[okr]) / sum(cs[okr]) else NA_real_,
      # Road distance, and denominated over the SAME routable cases as case_mean_km,
      # so the whole reach block is one metric with one denominator.
      case_far_share  = if (any(okr)) sum(cs[okr][drd[okr] > REACH_FAR_KM]) / sum(cs[okr])
                        else NA_real_,
      truncated = ti > fit_max)
  }))

  # Back-transform the mean centres to WGS84 so every displacement is geodesic.
  .to_ll <- function(x, y) {
    out <- matrix(NA_real_, length(x), 2)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) {
      p <- sf::st_sfc(lapply(which(ok), function(i) sf::st_point(c(x[i], y[i]))),
                      crs = geo$crs_laea)
      out[ok, ] <- sf::st_coordinates(sf::st_transform(p, 4326))[, 1:2, drop = FALSE]
    }
    out
  }
  ll_c <- .to_ll(weekly$cog_x, weekly$cog_y)
  ll_i <- .to_ll(weekly$inv_cog_x, weekly$inv_cog_y)
  weekly$cog_lon <- ll_c[, 1]; weekly$cog_lat <- ll_c[, 2]
  weekly$inv_cog_lon <- ll_i[, 1]; weekly$inv_cog_lat <- ll_i[, 2]

  sc <- .steps(weekly$cog_lon, weekly$cog_lat, weekly$week)
  si <- .steps(weekly$inv_cog_lon, weekly$inv_cog_lat, weekly$week)
  weekly$cog_step_km <- sc$disp
  weekly$cog_speed_km_wk <- sc$speed; weekly$cog_bearing_deg <- sc$bearing
  weekly$inv_cog_speed_km_wk <- si$speed; weekly$inv_cog_bearing_deg <- si$bearing
  weekly$cog_dist_from_epi_km <- vapply(seq_len(nrow(weekly)), function(i)
    .geo_km(epi_ll, c(weekly$cog_lon[i], weekly$cog_lat[i])), numeric(1))

  # ---- speed estimators (fit window only) ----------------------------------
  fw <- weekly %>% dplyr::filter(!truncated)

  v_radial_road <- .ols(fitz$inv_week, fitz$d_road)
  v_radial_gc   <- .ols(fitz$inv_week, fitz$d_gc)
  v_radial_time <- .ols(fitz$inv_week, fitz$t_road)   # travel-time analogue (h/week)
  v_front       <- .ols(fw$week, fw$reach_max_km)
  v_front_p90   <- .ols(fw$week, fw$reach_p90_km)
  v_reach_med   <- .ols(fw$week, fw$reach_med_km)
  v_case_mean   <- .ols(fw$week, fw$case_mean_km)
  v_disp        <- .ols(fw$week, fw$cog_sd_km)
  v_disp2       <- .ols(fw$week, fw$cog_sd_km^2)          # SD^2 ~ 4 D t  =>  D = slope/4

  # SATURATION of the envelope. If the epidemic reaches its full extent early by
  # long-range seeding and then stops expanding, a single OLS slope over the whole window
  # is NOT a wavefront speed. Split it: the slope while the envelope was being established,
  # and the slope after it first reached 90% of its final value.
  #
  # THE ENVELOPE IS THE 90th-PERCENTILE REACH R_p90(t), NOT THE RUNNING MAXIMUM R_max(t).
  # R_max is set by ONE zone. In the run published 2026-10-13 (index zone then Mongbwalu) a single
  # seeding 1,606 km out at week 18 moved it from 948 km to 1,606 km, so the "week R_max reached 90% of its final
  # value" was simply that zone's week, the "speed to 90%" was that zone's distance divided
  # by 18, and the slope "after 90%" was zero BY CONSTRUCTION, because a running maximum is
  # flat once its final value is set. R_p90 moves only when a tenth of the invaded set moves,
  # and it is not monotone, so its post-threshold slope is a measurement, not an identity.
  #
  # Both phases are OLS slopes of the SAME series (weeks 0..wk90 and wk90..end of the fit
  # window), so the two are like for like; dividing R(wk90) by wk90 would instead assume
  # the envelope started at zero, which it does not (week 0 is the seed zones' own reach).
  # A slope needs >= 3 fit-window weeks on a side (.ols()); with fewer the verdict is
  # UNDETERMINED rather than defaulting to either answer. The verdict rule itself — the
  # post-threshold slope at most a tenth of the establishing slope — is the one used before.
  p90_fin <- fw$reach_p90_km[nrow(fw)]
  i90 <- which(is.finite(fw$reach_p90_km) & fw$reach_p90_km >= 0.9 * p90_fin)
  wk90 <- if (is.finite(p90_fin) && length(i90)) fw$week[i90[1]] else NA_real_
  early <- fw %>% dplyr::filter(is.finite(wk90), week <= wk90)
  late  <- fw %>% dplyr::filter(is.finite(wk90), week >= wk90)
  v_env_early <- .ols(early$week, early$reach_p90_km)
  v_env_late  <- .ols(late$week,  late$reach_p90_km)
  envelope_verdict <- if (!is.finite(v_env_early$slope) || !is.finite(v_env_late$slope)) {
    sprintf(paste0("UNDETERMINED: %d fit-window week(s) up to the 90%% threshold and %d from it; ",
                   "a slope needs at least 3 on each side"), v_env_early$n, v_env_late$n)
  } else if (abs(v_env_late$slope) <= 0.1 * abs(v_env_early$slope)) {   # <=: a flat-from-the-start envelope (0 vs 0) is saturated
    sprintf("SATURATED: %+.1f km/wk from week %d against %+.1f km/wk before it; no travelling wavefront",
            v_env_late$slope, as.integer(wk90), v_env_early$slope)
  } else {
    sprintf("NOT SATURATED (%s): %+.1f km/wk from week %d against %+.1f km/wk before it",
            if (v_env_late$slope > 0) "still expanding" else "contracting",
            v_env_late$slope, as.integer(wk90), v_env_early$slope)
  }
  # A saturated envelope does NOT mean long-range introductions stopped — only that
  # the bulk of the invaded set stopped moving outward. Quantify what kept arriving after
  # the threshold week.
  after <- fw %>% dplyr::filter(is.finite(wk90), week > wk90, new_zones > 0)
  newmax_med_after <- stats::median(after$new_max_km, na.rm = TRUE)
  newmax_max_after <- .mx(after$new_max_km)
  n_weeks_after    <- nrow(after)

  cg_ok <- fw %>% dplyr::filter(is.finite(cog_lon))
  n_ok  <- nrow(cg_ok)
  net_km <- if (n_ok >= 2L)
    .geo_km(c(cg_ok$cog_lon[1], cg_ok$cog_lat[1]),
            c(cg_ok$cog_lon[n_ok], cg_ok$cog_lat[n_ok])) else NA_real_
  net_wk <- if (n_ok >= 2L) cg_ok$week[n_ok] - cg_ok$week[1] else NA_real_
  path_km <- sum(fw$cog_step_km, na.rm = TRUE)
  cog_straight <- if (is.finite(net_km) && is.finite(path_km) && path_km > 0)
    net_km / path_km else NA_real_
  cog_bearing <- if (n_ok >= 2L)
    .geo_brg(c(cg_ok$cog_lon[1], cg_ok$cog_lat[1]),
             c(cg_ok$cog_lon[n_ok], cg_ok$cog_lat[n_ok])) else NA_real_

  # Invaded-zone accrual: log-linear growth of the cumulative invaded count.
  acc_df <- fw %>% dplyr::filter(cum_zones > 0)
  acc <- .ols(acc_df$week, log(acc_df$cum_zones))
  dbl <- if (is.finite(acc$slope) && acc$slope > 0) log(2) / acc$slope else NA_real_

  # Which geography orders the arrival times best?
  .sp <- function(v) suppressWarnings(stats::cor(fitz$inv_week, v, method = "spearman",
                                                 use = "complete.obs"))
  metric_fit <- tibble::tibble(
    metric   = c("great-circle (km)", "road distance (km)", "travel time (h)"),
    spearman = c(.sp(fitz$d_gc), .sp(fitz$d_road), .sp(fitz$t_road)),
    r2       = c(v_radial_gc$r2, v_radial_road$r2, v_radial_time$r2),
    slope_per_week = c(v_radial_gc$slope, v_radial_road$slope, v_radial_time$slope))

  list(
    weeks = weeks, week_max = T_last, fit_max = fit_max, n_trunc = n_trunc,
    epi_zone = epi_zone, epi_ll = epi_ll, epi_province = epi$province,
    zones = zones_inv, zone_dist = zc, weekly = weekly, metric_fit = metric_fit,
    seed_rank = seed,
    stats = list(
      v_radial_road = v_radial_road$slope, v_radial_road_r2 = v_radial_road$r2,
      v_radial_gc   = v_radial_gc$slope,   v_radial_gc_r2   = v_radial_gc$r2,
      v_radial_time = v_radial_time$slope, v_radial_time_r2 = v_radial_time$r2,
      v_radial_road_ts = .theil_sen(fitz$inv_week, fitz$d_road),
      v_radial_gc_ts   = .theil_sen(fitz$inv_week, fitz$d_gc),
      v_front_road  = v_front$slope,       v_front_p90_road = v_front_p90$slope,
      v_env_p90_early = v_env_early$slope, v_env_p90_late = v_env_late$slope,
      n_weeks_env_early = v_env_early$n,   n_weeks_env_late = v_env_late$n,
      week_p90_reach_90pct = wk90,         envelope_verdict = envelope_verdict,
      newmax_med_after = newmax_med_after, newmax_max_after = newmax_max_after,
      n_weeks_after = n_weeks_after,
      v_reach_med = v_reach_med$slope, v_case_mean = v_case_mean$slope,
      diffusion_r2 = v_disp2$r2,
      v_cog_med     = stats::median(fw$cog_speed_km_wk, na.rm = TRUE),
      v_cog_iqr_lo  = .q(fw$cog_speed_km_wk, 0.25),
      v_cog_iqr_hi  = .q(fw$cog_speed_km_wk, 0.75),
      v_invcog_med  = stats::median(fw$inv_cog_speed_km_wk, na.rm = TRUE),
      v_cog_net     = if (is.finite(net_km) && is.finite(net_wk) && net_wk > 0)
                        net_km / net_wk else NA_real_,
      cog_net_disp_km = net_km, cog_path_km = path_km,
      cog_straightness = cog_straight, cog_net_bearing_deg = cog_bearing,
      v_disp_km_wk  = v_disp$slope, diffusion_km2_wk = v_disp2$slope / 4,
      accrual_rate_wk = acc$slope, accrual_doubling_wk = dbl,
      n_zones_fit = nrow(fitz), n_zones_all = nrow(zones_inv),
      reach_max_km = fw$reach_max_km[nrow(fw)], reach_med_km = fw$reach_med_km[nrow(fw)],
      cases_fit = sum(fw$cases)))
}

# =============================================================================
# SECTION 4 — run: primary (seeded) analysis
# =============================================================================
# THIS FILE IS A SCRIPT, NOT A LIBRARY. Everything below EXECUTES and overwrites the
# spread-kinematics and arrival-predictor tables in outputs/key_outputs/. run_all.R launches
# it as its own subprocess for that reason. source()-ing it to get at build_geo() or
# compute_kinematics() used to run the whole analysis silently — which is how a tooling sweep
# that loaded "every module" rewrote six published files with no indication it had.
#
# Guarded so only `Rscript 43_spread_kinematics.R` runs it. Sourcing now defines the functions
# and stops, with one message saying so. Set SPREAD_FORCE_RUN=1 to run it from a source() on
# purpose.
if (!is_script_run("43_spread_kinematics.R") &&
    !identical(Sys.getenv("SPREAD_FORCE_RUN"), "1")) {
  message("[kinematics] 43_spread_kinematics.R sourced, not run: functions are defined, ",
          "nothing is written. Use `Rscript 43_spread_kinematics.R` (or SPREAD_FORCE_RUN=1) ",
          "to run the analysis.")
} else {

message("\n", strrep("=", 74))
message("[kinematics] 43_spread_kinematics.R — descriptive spatiotemporal summaries")
message(strrep("=", 74))

# load_linelist() derives each record's onset draw from a hash of (alert_id, GLOBAL
# `ONSET_IMPUTE_SEED`) (01_data_prep.R; 00_config.R defaults it to RANDOM_SEED) and uses no
# RNG stream. Two consequences:
#   1. the primary run needs no set.seed() here and is byte-identical to the line
#      list every other module in the suite sees — the invasion dates below are the
#      same events the forecasting pipeline scores;
#   2. imputation replicates must vary ONSET_IMPUTE_SEED ITSELF — the onset draw's only seed;
#      RANDOM_SEED no longer reaches it. A set.seed() around the call is a no-op and would
#      yield N identical "replicates", i.e. a zero-width uncertainty interval that looks
#      like precision.
# RANDOM_SEED_BASE still seeds this script's own bootstrap (Section 5); ONSET_SEED_BASE is the
# primary completion's onset seed, and replicate r is ONSET_SEED_BASE + r.
RANDOM_SEED_BASE <- as.integer(get0("RANDOM_SEED", ifnotfound = 20260704L))
ONSET_SEED_BASE  <- as.integer(get0("ONSET_IMPUTE_SEED", ifnotfound = RANDOM_SEED_BASE))
ll_primary <- load_linelist()
# The imputed share, as load_linelist() measured it on THIS frame (never a typed literal).
message(sprintf("[kinematics] %.1f%% of confirmed records have an imputed onset (%s)",
                get0(".PCT_ONSET_IMPUTED", ifnotfound = NA_real_),
                get0(".ONSET_IMPUTE_MODE", ifnotfound = "mode not published")))
pop        <- load_population()
zw_primary <- aggregate_to_zone_week(ll_primary, names(pop))
weeks_all  <- sort(unique(zw_primary$week_start))

# The equal-area plane is centred on the case-weighted centre of the seed window, so
# planar distortion is smallest exactly where the mean centres live.
geo0 <- build_geo()
SEED_CUTOFF <- weeks_all[min(SEED_WEEKS, length(weeks_all))]   # calendar-anchored
seed_ctr <- zw_primary %>%
  dplyr::filter(week_start <= SEED_CUTOFF, confirmed > 0) %>%
  dplyr::left_join(dplyr::select(geo0$ctr, health_zone, lon, lat), by = "health_zone") %>%
  dplyr::summarise(lon = stats::weighted.mean(lon, confirmed),
                   lat = stats::weighted.mean(lat, confirmed))
geo <- build_geo(proj_centre = c(seed_ctr$lon, seed_ctr$lat))
message(sprintf("[kinematics] equal-area plane centred at %.3fE %.3fN (seed-window case centre)",
                seed_ctr$lon, seed_ctr$lat))

# ---- right-truncation: how many trailing weeks to drop ---------------------
comp    <- week_completeness(ll_primary, weeks_all, ANALYSIS_DATE)
message(sprintf("[kinematics] delay ECDF: %s complete-follow-up pairs; %.1f%% of their delays exceed the %d-day cap (counted as arrived at the cap)",
                format(attr(comp, "n_pairs")), 100 * as.numeric(attr(comp, "share_over_cap")),
                DELAY_CAP_DAYS))
comp    <- as.numeric(comp)          # attributes reported; keep a plain vector downstream
n_trunc <- choose_trunc(comp, TRUNC_FORCE)
stopifnot(n_trunc < length(weeks_all) - 2L)
# Calendar anchor for the fit window, so every replicate and every sensitivity run
# means the same span of real time no matter how its own week grid comes out.
FIT_CUTOFF <- weeks_all[length(weeks_all) - n_trunc]
message(sprintf("[kinematics] estimated onset-week completeness: %s",
                paste(sprintf("%s=%.2f", format(weeks_all, "%m-%d"), comp), collapse = "  ")))
message(sprintf("[kinematics] fits use weeks up to %s; the last %d week(s) are dropped (right truncation)",
                FIT_CUTOFF, n_trunc))

K <- compute_kinematics(zw_primary, geo, fit_cutoff = FIT_CUTOFF, seed_cutoff = SEED_CUTOFF)
S <- K$stats
message(sprintf("[kinematics] index zone: %s (%s); %d invaded zones (%d inside the fit window)",
                K$epi_zone, K$epi_province, S$n_zones_all, S$n_zones_fit))

# ---- zone bootstrap CI on the radial expansion rate ------------------------
# Resamples INVADED ZONES with replacement (one row per zone; no clustering needed)
# and refits v_radial. The leading-edge slope has no honest bootstrap (a resampled
# running maximum is biased low by construction) and is reported point-only.
# 90% (5th-95th percentile), the suite's default level and the level of every other lo/hi pair
# in the summary CSV (week bootstrap, imputation percentiles) EXCEPT the case centre-of-gravity
# speed row, whose lo/hi are its interquartile range (its basis column says so).
boot_ci <- function(z, col, B = N_BOOT) {
  dd <- z[[col]]; tt <- z$inv_week
  ok <- is.finite(dd) & is.finite(tt); dd <- dd[ok]; tt <- tt[ok]; n <- length(dd)
  if (n < 5L) return(c(NA_real_, NA_real_))
  s <- vapply(seq_len(B), function(b) {
    i <- sample.int(n, n, replace = TRUE)
    if (length(unique(tt[i])) < 2L) return(NA_real_)
    unname(stats::coef(stats::lm(dd[i] ~ tt[i]))[2])
  }, numeric(1))
  c(.q(s, 0.05), .q(s, 0.95))
}
set_fixed_seed(SEED + 1L)
fitz_primary <- dplyr::filter(K$zones, in_fit_window)
ci_road <- boot_ci(fitz_primary, "d_road")
ci_gc   <- boot_ci(fitz_primary, "d_gc")

# ---- sensitivity to the index-zone choice ----------------------------------
# The top seed-window zones are separated by a few cases (counts in the diagnostics CSV), so the
# argmax rule is fragile by construction. Re-running the
# whole analysis from each of the top ORIGIN_SENS_N candidates shows whether that
# fragility actually moves any number.
ORIGIN_SENS_N <- 3L
orig_sens <- dplyr::bind_rows(lapply(utils::head(K$seed_rank$health_zone, ORIGIN_SENS_N),
  function(z) {
    Kz <- compute_kinematics(zw_primary, geo, fit_cutoff = FIT_CUTOFF,
                             seed_cutoff = SEED_CUTOFF, epi_zone = z)
    tibble::tibble(origin = z, primary = z == K$epi_zone,
                   seed_cases = K$seed_rank$seed_cases[match(z, K$seed_rank$health_zone)],
                   v_radial_road = Kz$stats$v_radial_road,
                   v_radial_road_r2 = Kz$stats$v_radial_road_r2,
                   v_reach_med = Kz$stats$v_reach_med,
                   reach_max_km = Kz$stats$reach_max_km,
                   reach_med_km = Kz$stats$reach_med_km)
  }))

# ---- sensitivity to the truncation choice ----------------------------------
sens <- dplyr::bind_rows(lapply(TRUNC_SENS, function(nt) {
  if (nt > length(weeks_all) - 3L) return(NULL)
  Kn <- compute_kinematics(zw_primary, geo,
                           fit_cutoff = weeks_all[length(weeks_all) - nt],
                           seed_cutoff = SEED_CUTOFF, epi_zone = K$epi_zone)
  tibble::tibble(trunc_weeks = nt, primary = nt == n_trunc,
                 v_radial_road = Kn$stats$v_radial_road,
                 v_front_road = Kn$stats$v_front_road,
                 v_cog_med = Kn$stats$v_cog_med,
                 v_disp_km_wk = Kn$stats$v_disp_km_wk,
                 n_zones_fit = Kn$stats$n_zones_fit)
}))

# Mean new invasions per week over weeks 1..fit_max of a compute_kinematics() weekly table:
# week 0 is the seeding condition (those zones already had confirmed cases when the series
# opens) and the truncated tail is incompletely observed. One definition for the primary run
# (the summary row, Section 5) and for every onset completion.
.nz_window <- function(weekly) {
  rows <- weekly$week >= 1L & !weekly$truncated & is.finite(weekly$new_zones)
  vals <- as.numeric(weekly$new_zones[rows])
  list(vals = vals, mean = if (length(vals)) mean(vals) else NA_real_,
       fit_max = if (any(rows)) max(weekly$week[rows]) else NA_integer_)
}
.verdict_class <- function(v) sub(":.*$", "", v)   # "SATURATED" / "NOT SATURATED (...)" / "UNDETERMINED"

# ---- onset-imputation replicates -------------------------------------------
# load_linelist() and compute_kinematics() (including the index-zone choice) are re-run on
# each onset completion. Every scalar statistic compute_kinematics() returns is kept, so each
# summary-CSV row it feeds carries an imputation interval (imp_map below); the zone bootstrap
# and the sensitivity tables are primary-completion analyses and are not re-run.
imp <- NULL
if (N_IMPUTE > 0L) {
  message(sprintf("[kinematics] %d onset-imputation replicate(s)...", N_IMPUTE))
  imp <- dplyr::bind_rows(lapply(seq_len(N_IMPUTE), function(r) {
    # Vary the seed load_linelist() actually reads (ONSET_IMPUTE_SEED, NOT RANDOM_SEED, which
    # no longer reaches the onset draw: varying it would give N identical replicates), and put
    # it back afterwards so the rest of the session still sees the project default.
    .old <- get0("ONSET_IMPUTE_SEED", envir = globalenv(), ifnotfound = NULL)
    assign("ONSET_IMPUTE_SEED", ONSET_SEED_BASE + r, envir = globalenv())
    on.exit(if (is.null(.old)) suppressWarnings(rm("ONSET_IMPUTE_SEED", envir = globalenv()))
            else assign("ONSET_IMPUTE_SEED", .old, envir = globalenv()), add = TRUE)
    out <- try({
      llr <- suppressMessages(suppressWarnings(load_linelist()))
      zwr <- suppressMessages(suppressWarnings(aggregate_to_zone_week(llr, names(pop))))
      Kr  <- compute_kinematics(zwr, geo, fit_cutoff = FIT_CUTOFF,
                                seed_cutoff = SEED_CUTOFF)
      dplyr::bind_cols(
        tibble::tibble(rep = r, epi_zone = Kr$epi_zone, n_weeks = length(Kr$weeks),
                       envelope_verdict = Kr$stats$envelope_verdict,
                       nz_mean = .nz_window(Kr$weekly)$mean,
                       # per-zone first confirmed onset on THIS draw, carried as a list
                       # column so Section 5b can put an imputation interval on the
                       # arrival-vs-predictor correlations too (the predictors are fixed
                       # geography and mobility; only the arrival dates are imputed).
                       arrivals = list(llr %>%
                         dplyr::filter(confirmed %in% TRUE, !is.na(date_index), !is.na(health_zone)) %>%
                         dplyr::group_by(health_zone) %>%
                         dplyr::summarise(first_onset = min(date_index), .groups = "drop"))),
        # every scalar numeric statistic of this completion, under its compute_kinematics() name
        tibble::as_tibble(Filter(function(v) is.numeric(v) && length(v) == 1L, Kr$stats)))
    }, silent = TRUE)
    if (inherits(out, "try-error")) NULL else out
  }))
  if (nrow(imp) < N_IMPUTE)
    warning("[kinematics] ", N_IMPUTE - nrow(imp), " imputation replicate(s) failed.",
            call. = FALSE)
}
# 5th / 50th / 95th percentile over the COMPLETIONS, with the primary completion counted as
# replicate 0 (it is the ONSET_SEED_BASE draw; the replicates are BASE + 1..N). The published
# point VALUE stays the primary completion, because that is the completion every other module
# scores; the median is published beside it, so a value lying outside its own replicate band
# (the primary is one draw, not a summary) is visible as such rather than looking like an error.
NZ_PRIMARY <- .nz_window(K$weekly)
imp_q <- function(col) {
  if (is.null(imp) || !nrow(imp)) return(c(NA_real_, NA_real_, NA_real_))
  stopifnot(col %in% names(imp))
  v <- c(if (identical(col, "nz_mean")) NZ_PRIMARY$mean else S[[col]], imp[[col]])
  c(.q(v, 0.05), .q(v, 0.50), .q(v, 0.95))
}

# =============================================================================
# SECTION 5 — tabular outputs
# =============================================================================

# Zones whose confirmed cases exist ONLY through the sitrep top-up: their invasion
# week is dated from a sitrep REPORT date minus a delay draw, so it is late-biased.
sit_counts <- ll_primary %>% dplyr::filter(confirmed %in% TRUE) %>%
  dplyr::group_by(health_zone) %>%
  dplyr::summarise(n_all = dplyr::n(),
                   n_sit = sum(grepl("^SITREP-CONF", alert_id)), .groups = "drop")
sitrep_only_zones <- sit_counts$health_zone[sit_counts$n_sit == sit_counts$n_all]

# Leverage check on those zones. Their invasion week is LATE-biased (sitrep report date
# minus a delay draw rather than a real onset), and a late-dated DISTANT zone steepens
# the radial slope — so dropping them can only be checked, never assumed harmless.
fitz_nosit  <- fitz_primary %>% dplyr::filter(!health_zone %in% sitrep_only_zones)
nosit_road  <- .ols(fitz_nosit$inv_week, fitz_nosit$d_road)
nosit_n     <- nrow(fitz_primary) - nrow(fitz_nosit)

# Distances and bearing here are from the single INDEX zone, and the column names say so:
# arrival_predictors.csv (Section 5b) publishes same-unit columns measured from the nearest
# of the three EPICENTRE_ZONES, so the bare names would mean two different things in two
# files of the same run.
zones_out <- K$zones %>%
  dplyr::transmute(health_zone, province,
                   invasion_week = inv_week, invasion_week_start = inv_week_start,
                   cases_total,
                   dist_greatcircle_km_from_index_zone = round(d_gc, 2),
                   dist_road_km_from_index_zone = round(d_road, 2),
                   travel_time_h_from_index_zone = round(t_road, 3),
                   bearing_from_index_zone_deg = round(bearing, 1),
                   lon = round(lon, 5), lat = round(lat, 5),
                   in_fit_window,
                   sitrep_only = health_zone %in% sitrep_only_zones)
readr::write_csv(zones_out, file.path(KEY_DIR, "spread_kinematics_zones.csv"))

# -----------------------------------------------------------------------------
# SECTION 5b — arrival time vs candidate predictors
#   publishes: arrival_predictors.csv, arrival_predictor_fits.csv
#   consumed by: make_manuscript_figures.R (Figure 1C and Figure S1),
#                45_ascertainment_figures.R (arrival_predictors.csv)
# -----------------------------------------------------------------------------
# Figure 1C asks a descriptive question — which geography orders the arrivals? —
# and answers it with a per-predictor OLS fit and a Spearman rho. Those statistics
# are computed HERE, not in the figure. The panel previously rebuilt its own
# coordinates (by regex over a raw HTML file), its own epicentre, its own road and
# great-circle distances and its own mobility shares, and so published correlations
# that no table in this pipeline could reproduce and that disagreed with
# spread_kinematics_zones.csv on both the zone set and the arrival dates.
#
# EPICENTRE. The kinematics above are radial from the single INDEX zone (the zone
# with most confirmed cases in the seed window). This table instead measures access
# to the epicentre REGION — the three zones seeded together, EPICENTRE_ZONES (M+B+R)
# — because that is the exposure the manuscript panel claims. Distance and travel
# time are the MINIMUM over the three (a zone is reached via its nearest seed);
# mobility is the flow leaving the region as a whole. The index-zone distances
# remain available in spread_kinematics_zones.csv (the *_from_index_zone columns), so the
# two are comparable rather than conflated.
#
# ARRIVAL. Days from ARRIVAL_ORIGIN, the earliest confirmed date_index in the
# reconciled line list — derived from the data, not typed. date_index is
# load_linelist()'s, i.e. the same onset draw the forecasting pipeline scores, including
# the sitrep top-up. load_linelist() floors date_index at OUTBREAK_START, so when an
# imputed onset is clamped there the origin is that floor rather than an observed onset;
# the run message says whether it is.
#
# FIT WINDOW. A zone whose first onset falls in the right-truncated tail (after the last
# day of the FIT_CUTOFF week) is in the table, flagged `in_fit_window = FALSE` exactly as
# in spread_kinematics_zones.csv, and in no fit: its week is incompletely observed, so the
# zones that arrived then and are already confirmed are the quickly reported subset.
#
# MOBILITY PREDICTORS. -log(share of flow from the epicentre region into the zone).
# That is the one-step Brockmann-Helbing effective distance up to the additive
# constant 1 (d_eff = 1 - log p) and NOT the shortest-path D_eff minimised over all
# routes, so the columns are named for what they are. A zone with no flow released
# from the epicentre has an undefined value and is carried as NA rather than floored
# at some arbitrary small share: flooring would place exactly the zones with least
# measured connectivity at a fabricated finite predictor value, at the top of the
# range being fitted. They are therefore DROPPED, the per-predictor n differs
# between columns, and every consumer must report n alongside R^2 and rho.
#
# THE EPICENTRE ZONES ARE IN THIS TABLE BUT IN NO FIT. They are flagged `epicentre`
# and excluded from every predictor's fit alike (see arrival_predictor_fits): they are
# the origin the predictors are measured from, not arrivals. A consumer that draws the
# points must draw the same set it reports the fit for, i.e. filter
# `!epicentre & in_fit_window`.
FIT_LAST_DAY <- FIT_CUTOFF + 6L   # last day of the last fit-window week (weeks run week_start..+6)
arrival_predictors <- local({
  zones_all <- names(pop)
  aliases   <- tryCatch(load_aliases(), error = function(e) { 
    warning("[kinematics] alias table unavailable; zone names used as-is.", call. = FALSE)
    character(0)
  })
  .can <- function(x) if (length(aliases)) suppressWarnings(harmonise_names(x, aliases, zones_all)) else x

  # collapse a named vector onto canonical names, summing zones that merge
  .collapse <- function(v) {
    v <- v[!is.na(names(v))]
    tapply(v, names(v), function(z) if (all(is.na(z))) NA_real_ else sum(z, na.rm = TRUE))
  }

  # ---- arrival time ---------------------------------------------------------
  arr <- ll_primary %>%
    dplyr::filter(confirmed %in% TRUE, !is.na(date_index), !is.na(health_zone)) %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(first_onset = min(date_index), cases_total = dplyr::n(), .groups = "drop")
  stopifnot(nrow(arr) > 0L)
  arrival_origin <- min(arr$first_onset)

  # ---- epicentre region -----------------------------------------------------
  epi_set <- intersect(.can(EPICENTRE_ZONES), zones_all)
  if (length(epi_set) != length(EPICENTRE_ZONES))
    warning(sprintf("[kinematics] %d of %d epicentre zone(s) resolved against the zone spine (%s).",
                    length(epi_set), length(EPICENTRE_ZONES), paste(epi_set, collapse = ", ")),
            call. = FALSE)
  stopifnot(length(epi_set) > 0L)

  # ---- distance and travel time: MIN over the epicentre set -----------------
  .min_over_epi <- function(M) {
    ez <- intersect(epi_set, rownames(M)); stopifnot(length(ez) > 0L)
    sub <- M[ez, , drop = FALSE]
    apply(sub, 2L, function(z) { z <- z[is.finite(z)]; if (length(z)) min(z) else NA_real_ })
  }
  road_km_v <- .min_over_epi(geo$road)   # km
  road_h_v  <- .min_over_epi(geo$time)   # hours (build_geo already divides by 60)

  ctr    <- geo$ctr
  epi_ll <- as.matrix(ctr[match(epi_set, ctr$health_zone), c("lon", "lat"), drop = FALSE])
  stopifnot(all(is.finite(epi_ll)))
  gc_mat <- matrix(vapply(seq_len(nrow(epi_ll)),
                          function(i) as.numeric(geosphere::distGeo(epi_ll[i, ], cbind(ctr$lon, ctr$lat))) / 1000,
                          numeric(nrow(ctr))),
                   nrow = nrow(ctr))
  gc_km_v <- setNames(apply(gc_mat, 1L, min, na.rm = TRUE), ctr$health_zone)

  # ---- relocation OD share out of the epicentre region (month: FLOWMINDER_OD_MONTH) ----
  # NA in the OD matrix means "no count released for this pair this month", not zero
  # (see load_flowminder_od). A destination for which EVERY epicentre row is NA has
  # an unobserved share and is carried as NA; where some rows are released, the
  # observed ones are summed. Shares are normalised over the observed destinations.
  M_od <- load_flowminder_od()
  rownames(M_od) <- .can(rownames(M_od)); colnames(M_od) <- .can(colnames(M_od))
  ei <- intersect(epi_set, rownames(M_od))
  if (!length(ei)) {
    warning("[kinematics] no epicentre zone present in the Flowminder OD matrix; ",
            "the relocation predictor is unavailable.", call. = FALSE)
    od_share_v <- setNames(rep(NA_real_, nrow(ctr)), ctr$health_zone)
  } else {
    sub  <- M_od[ei, , drop = FALSE]
    n_ob <- colSums(!is.na(sub))
    raw  <- colSums(sub, na.rm = TRUE)
    raw[n_ob == 0L] <- NA_real_
    raw <- .collapse(raw)
    od_share_v <- raw / sum(raw, na.rm = TRUE)
  }

  # ---- short-trip snapshot (pooled M+B+R cohort) -----------------------------
  # SNAPSHOT: build_M1's rule — the LATEST on or before the analysis date, else (a back-dated
  # run older than every snapshot) the EARLIEST, which leaks least; never the newest overall.
  # SOURCE: load_short_trip_profile(), i.e. the RAW annex table resolved against the zone
  # spine with its province column — the source build_M1 and M2 use. The processed matrix
  # (load_short_trip_snapshot()) was canonicalised upstream against the 437-zone OD header and
  # silently lacks real spine destinations (see load_short_trip_profile and
  # .warn_if_processed_source_stale), which then dropped out of this predictor as NA. The
  # profile's proportions are of the pooled epicentre cohort, so there is no epicentre
  # subsetting to do here; they are renormalised because the released percentages need not
  # sum to exactly 100. The loader itself falls back to the processed matrix, with a warning,
  # only if the raw table is absent or unusable; the source used is printed below.
  .st_elig <- which(FLOWMINDER_ST_DATES <= ANALYSIS_DATE)
  st_tag <- if (length(.st_elig)) FLOWMINDER_ST_TAGS[.st_elig[which.max(FLOWMINDER_ST_DATES[.st_elig])]]
            else FLOWMINDER_ST_TAGS[which.min(FLOWMINDER_ST_DATES)]
  st_date <- FLOWMINDER_ST_DATES[match(st_tag, FLOWMINDER_ST_TAGS)]
  st_prof <- load_short_trip_profile(st_tag, zones_all, if (length(aliases)) aliases else NULL)
  st_v    <- st_prof$props / sum(st_prof$props, na.rm = TRUE)

  .lk <- function(v, z) { out <- as.numeric(v[z]); out[is.na(match(z, names(v)))] <- NA_real_; out }
  .neglog <- function(p) ifelse(is.finite(p) & p > 0, -log(p), NA_real_)

  out <- arr %>%
    dplyr::mutate(
      province            = ctr$province[match(health_zone, ctr$health_zone)],
      # The epicentre set is the SOURCE of every predictor below, not an arrival: its
      # distance and travel time are 0 by definition and its own outflow share is undefined
      # (the short-trip cohort IS these zones). Flagged so the fits can exclude it from EVERY
      # predictor alike — see arrival_predictor_fits.
      epicentre           = health_zone %in% epi_set,
      # first onset on or before the last day of the FIT_CUTOFF week: the right-truncated tail
      # enters no fit (the same rule, and the same flag, as spread_kinematics_zones.csv)
      in_fit_window       = first_onset <= FIT_LAST_DAY,
      arrival_origin      = arrival_origin,
      arrival_days        = as.numeric(first_onset - arrival_origin),
      shorttrip_snapshot  = st_date,
      dist_greatcircle_km = .lk(gc_km_v,  health_zone),
      dist_road_km        = .lk(road_km_v, health_zone),
      travel_time_h       = .lk(road_h_v,  health_zone),
      share_shorttrip     = .lk(st_v,      health_zone),
      share_relocation    = .lk(od_share_v, health_zone),
      population          = .lk(pop,       health_zone)) %>%
    dplyr::mutate(
      neglog_share_shorttrip  = .neglog(share_shorttrip),
      neglog_share_relocation = .neglog(share_relocation),
      log10_population        = ifelse(is.finite(population) & population > 0, log10(population), NA_real_)) %>%
    dplyr::arrange(arrival_days, health_zone)

  attr(out, "epi_set")          <- epi_set
  attr(out, "shorttrip_tag")    <- st_tag
  attr(out, "shorttrip_source") <- st_prof$source
  attr(out, "arrival_origin")   <- arrival_origin
  out
})

# The fit-window flag must be the kinematics' own: same line list, same weeks, so a zone is in
# the arrival fits exactly when compute_kinematics() counted it inside the fit window.
.fw_chk <- match(K$zones$health_zone, arrival_predictors$health_zone)
stopifnot(!anyNA(.fw_chk),
          identical(arrival_predictors$in_fit_window[.fw_chk], K$zones$in_fit_window))

ARRIVAL_ORIGIN <- attr(arrival_predictors, "arrival_origin")
# Whether the origin is the date_index floor, and how many of the confirmed onsets dated on it
# are imputed (an imputed onset that would fall earlier is clamped there) — counts only.
.on_floor <- ll_primary$confirmed %in% TRUE & !is.na(ll_primary$date_index) &
  ll_primary$date_index == ARRIVAL_ORIGIN
message(sprintf("[kinematics] arrival-predictor table: %d zone(s), %d in the fit window (first onset <= %s), origin %s%s, epicentre = %s, short-trip snapshot %s (%s source)",
                nrow(arrival_predictors), sum(arrival_predictors$in_fit_window),
                format(FIT_LAST_DAY, "%Y-%m-%d"), format(ARRIVAL_ORIGIN, "%Y-%m-%d"),
                if (ARRIVAL_ORIGIN == as.Date(OUTBREAK_START))
                  sprintf(" (= OUTBREAK_START, the date_index floor; %d of the %d confirmed onset(s) on it imputed)",
                          sum(.on_floor & ll_primary$onset_imputed %in% TRUE), sum(.on_floor))
                else "",
                paste(attr(arrival_predictors, "epi_set"), collapse = "+"),
                attr(arrival_predictors, "shorttrip_tag"),
                attr(arrival_predictors, "shorttrip_source")))

# Written at FULL WORKING PRECISION, not rounded for reading. The fits below are the
# panel's published statistics and a reader must be able to recompute them from this
# table exactly; rounding the predictors to 2-4 decimals for legibility moved the
# recomputed Spearman rho by ~5e-4, which is enough to disagree with the printed
# 2-decimal annotation and to make the panel unverifiable. signif(, 12) is lossless at
# these magnitudes and still keeps the file readable.
readr::write_csv(
  arrival_predictors %>%
    dplyr::transmute(health_zone, province, epicentre, in_fit_window,
                     first_onset, arrival_origin, arrival_days,
                     dplyr::across(c(dist_greatcircle_km, dist_road_km, travel_time_h,
                                     share_shorttrip, share_relocation,
                                     neglog_share_shorttrip, neglog_share_relocation,
                                     population, log10_population),
                                   ~ signif(.x, 12)),
                     cases_total,
                     # the snapshot share_shorttrip was read from (one date for the table), so a
                     # figure can name its month without re-deriving the selection rule
                     shorttrip_snapshot),
  file.path(KEY_DIR, "arrival_predictors.csv"))

# Per-predictor fit of PREDICTOR on ARRIVAL DAYS. The orientation matters and is the
# panel's: arrival time on x, predictor on y, so `slope` is in predictor-units per day
# and (intercept, slope) draw exactly the line the panel shows — the figure plots this
# line from these two numbers rather than refitting. R^2 is symmetric in x and y, rho
# is rank-based, so neither depends on the orientation; the slope does.
arrival_predictor_fits <- local({
  # Month NAME by lookup (month.name), not format(, "%B"), which follows the session locale.
  .month_of <- function(d) if (length(d) == 1L && !is.na(d))
    paste0(month.name[as.integer(format(as.Date(d), "%m"))], " ") else ""
  preds <- c(neglog_share_shorttrip  = sprintf("-log(share of the epicentre region's %sshort-trip outflow)",
                                               .month_of(arrival_predictors$shorttrip_snapshot[1])),
             neglog_share_relocation = sprintf("-log(share of the epicentre region's %srelocation outflow)",
                                               .month_of(get0("FLOWMINDER_OD_MONTH", ifnotfound = as.Date(NA)))),
             dist_greatcircle_km     = "great-circle distance from the nearest epicentre zone (km)",
             dist_road_km            = "road distance from the nearest epicentre zone (km)",
             travel_time_h           = "road travel time from the nearest epicentre zone (h)",
             log10_population        = "population (log10 people)")
  # one predictor, one arrival vector -> the panel's statistics
  .fit_one <- function(x, y) {
    k <- is.finite(x) & is.finite(y)
    n <- sum(k)
    if (n < 3L || length(unique(x[k])) < 2L || length(unique(y[k])) < 2L)
      return(list(n = n, r2 = NA_real_, intercept = NA_real_, slope = NA_real_,
                  pearson = NA_real_, spearman = NA_real_))
    fit <- stats::lm(y[k] ~ x[k])
    list(n = n,
         r2        = summary(fit)$r.squared,
         intercept = unname(stats::coef(fit)[1L]),
         slope     = unname(stats::coef(fit)[2L]),
         pearson   = suppressWarnings(stats::cor(x[k], y[k])),
         spearman  = suppressWarnings(stats::cor(x[k], y[k], method = "spearman")))
  }

  # Imputation uncertainty. The predictors are fixed geography and mobility; only the
  # arrival dates are imputed, so each replicate re-dates the SAME zones and refits. The zone
  # set itself does not change (every confirmed record has a date_index, so the onset seed
  # moves arrival dates, never which zones are invaded), but a completion's dates can move a
  # zone's first onset across FIT_LAST_DAY. Each completion is therefore held to its OWN fit
  # window — as compute_kinematics() holds each replicate to the same calendar cutoff — and
  # its n may differ from the primary n; it is not asserted equal.
  rep_arr <- if (!is.null(imp) && nrow(imp) && "arrivals" %in% names(imp)) imp$arrivals else list()
  .rep_x <- lapply(rep_arr, function(a) {
    if (is.null(a) || !nrow(a)) return(NULL)
    fo <- a$first_onset[match(arrival_predictors$health_zone, a$health_zone)]
    xr <- as.numeric(fo - min(a$first_onset))
    xr[!(fo <= FIT_LAST_DAY)] <- NA_real_          # this completion's right-truncated tail
    xr
  })
  .rep_x <- Filter(Negate(is.null), .rep_x)

  # THE EPICENTRE SET IS EXCLUDED FROM EVERY FIT. Its zones are the source the predictors are
  # measured FROM: distance and travel time are 0 by definition and their arrival days are the
  # earliest of all (arrival_predictors.csv lists them), while their outflow share is
  # undefined. They used to sit in the distance, travel-time and relocation fits as three
  # perfectly concordant points at the origin of both axes, and to drop out of the short-trip
  # fit as NA, so the predictors were compared on DIFFERENT supports and the geographic ones
  # carried an advantage that measured nothing about arrival order. Every predictor is now fitted
  # on the invaded zones OUTSIDE the epicentre set and INSIDE the fit window (`support` says so
  # on every row), and n is published per fit, because the mobility predictors still drop the
  # zones for which no flow out of the epicentre was released.
  #
  # THE COMMON-SUPPORT COLUMNS (*_common) close the remaining gap: every predictor refitted on
  # the ONE zone set on which ALL of them are defined, so R^2 and rho can be compared across
  # predictors like for like (n_common is the same on every row by construction).
  epi   <- arrival_predictors$epicentre %in% TRUE
  infit <- arrival_predictors$in_fit_window %in% TRUE
  x_all <- arrival_predictors$arrival_days
  x_all[!infit] <- NA_real_                          # the right-truncated tail enters no fit
  .y_of <- function(p) { y <- arrival_predictors[[p]]; y[epi] <- NA_real_; y }
  common <- is.finite(x_all) & !epi &
    Reduce(`&`, lapply(names(preds), function(p) is.finite(arrival_predictors[[p]])))
  .support <- sprintf(paste0("invaded zones outside the epicentre set (%s), first onset on or before %s ",
                             "(the right-truncated tail is excluded, as from every kinematics fit), ",
                             "with this predictor defined"),
                      paste(attr(arrival_predictors, "epi_set"), collapse = "+"),
                      format(FIT_LAST_DAY, "%Y-%m-%d"))
  dplyr::bind_rows(lapply(names(preds), function(p) {
    y  <- .y_of(p)
    f  <- .fit_one(x_all, y)
    fc <- .fit_one(x_all[common], y[common])
    k  <- is.finite(x_all) & is.finite(y)
    # Imputation interval over the COMPLETIONS, the primary counted as replicate 0 — the same
    # convention as the summary CSV's imp_* columns, so a published value can no longer sit
    # outside its own interval merely because the primary draw was left out of it.
    rr <- cbind(c(f$r2, f$spearman),
                if (length(.rep_x)) vapply(.rep_x, function(xr) {
                  g <- .fit_one(xr, y); c(g$r2, g$spearman) }, numeric(2))
                else matrix(numeric(0), nrow = 2))
    .qq <- function(v, q) if (length(v) && any(is.finite(v))) unname(stats::quantile(v[is.finite(v)], q, names = FALSE)) else NA_real_
    tibble::tibble(
      predictor = p, label = unname(preds[p]), support = .support, n = f$n,
      r2 = f$r2, intercept = f$intercept, slope = f$slope,
      pearson = f$pearson, spearman = f$spearman,
      n_common = fc$n, r2_common = fc$r2, spearman_common = fc$spearman,
      n_imputation_completions = ncol(rr),
      r2_lo       = .qq(rr[1, ], 0.05), r2_med       = .qq(rr[1, ], 0.50),
      r2_hi       = .qq(rr[1, ], 0.95),
      spearman_lo = .qq(rr[2, ], 0.05), spearman_med = .qq(rr[2, ], 0.50),
      spearman_hi = .qq(rr[2, ], 0.95),
      # Beyond the epicentre exclusion the per-predictor fits are still NOT on common support,
      # and that dropping is selection ON the predictor, which is correlated with the outcome:
      # a zone is dropped precisely when no flow out of the epicentre was released for it.
      # These two columns make that selection visible (both over the non-epicentre zones
      # inside the fit window), so a reader comparing R^2 across predictors can see whether
      # the denominators describe the same outbreak; the *_common columns remove it.
      mean_arrival_days_used    = if (any(k)) mean(x_all[k]) else NA_real_,
      mean_arrival_days_dropped = { d <- is.finite(x_all) & !epi & !is.finite(y)
                                    if (any(d)) mean(x_all[d]) else NA_real_ })
  }))
})
readr::write_csv(
  arrival_predictor_fits %>%
    dplyr::mutate(dplyr::across(c(r2, intercept, slope, pearson, spearman,
                                  r2_common, spearman_common,
                                  r2_lo, r2_med, r2_hi, spearman_lo, spearman_med, spearman_hi,
                                  mean_arrival_days_used, mean_arrival_days_dropped),
                                ~ signif(.x, 12))),
  file.path(KEY_DIR, "arrival_predictor_fits.csv"))
message(sprintf("[kinematics] arrival-predictor fits exclude the epicentre set (%s) and the %d zone(s) first seen after %s (right-truncated tail); common support n = %d",
                paste(attr(arrival_predictors, "epi_set"), collapse = "+"),
                sum(!arrival_predictors$in_fit_window), format(FIT_LAST_DAY, "%Y-%m-%d"),
                arrival_predictor_fits$n_common[1]))
for (i in seq_len(nrow(arrival_predictor_fits))) with(arrival_predictor_fits[i, ],
  message(sprintf("[kinematics]   %-24s n=%2d  R2=%.3f  rho=%+.3f | common n=%2d  R2=%.3f  rho=%+.3f",
                  predictor, n, r2, spearman, n_common, r2_common, spearman_common)))

weekly_out <- K$weekly %>%
  dplyr::mutate(completeness_est = comp[match(week_start, weeks_all)]) %>%
  dplyr::transmute(week, week_start, cases, new_zones, cum_zones, cum_provinces,
                   cog_lon = round(cog_lon, 5), cog_lat = round(cog_lat, 5),
                   cog_sd_km = round(cog_sd_km, 2),
                   cog_major_km = round(cog_major_km, 2),
                   cog_minor_km = round(cog_minor_km, 2),
                   cog_orient_deg = round(cog_orient_deg, 1),
                   cog_step_km = round(cog_step_km, 2),
                   cog_speed_km_wk = round(cog_speed_km_wk, 2),
                   cog_bearing_deg = round(cog_bearing_deg, 1),
                   # from the single index zone (cf. arrival_predictors.csv's M+B+R distances)
                   cog_dist_from_index_zone_km = round(cog_dist_from_epi_km, 2),
                   inv_cog_lon = round(inv_cog_lon, 5), inv_cog_lat = round(inv_cog_lat, 5),
                   inv_cog_speed_km_wk = round(inv_cog_speed_km_wk, 2),
                   inv_cog_bearing_deg = round(inv_cog_bearing_deg, 1),
                   invaded_sd_km = round(invaded_sd_km, 2),
                   reach_max_km = round(reach_max_km, 1),
                   reach_p90_km = round(reach_p90_km, 1),
                   reach_med_km = round(reach_med_km, 1),
                   reach_max_gc_km = round(reach_max_gc_km, 1),
                   new_max_km = round(new_max_km, 1),
                   new_med_km = round(new_med_km, 1),
                   case_mean_km = round(case_mean_km, 1),
                   case_far_share = round(case_far_share, 4),
                   completeness_est = round(completeness_est, 3), truncated)
readr::write_csv(weekly_out, file.path(KEY_DIR, "spread_kinematics_weekly.csv"))

# Mean new invasions per week over weeks 1..fit_max (.nz_window()), with a 90% bootstrap
# interval that resamples WEEKS as if exchangeable. They are not: new_zones trends over the
# window and is serially correlated, so this interval describes the week-to-week spread
# around the window mean, not the sampling uncertainty of a stationary invasion rate, and is
# labelled "week-bootstrap interval" (not a CI) wherever it is shown. Published for Panel C of
# Figure_spread_kinematics_compact, which must not compute its own reference level.
.nz_vals <- NZ_PRIMARY$vals
.nz_fit_max <- NZ_PRIMARY$fit_max
.nz_mean <- NZ_PRIMARY$mean
.nz_ci <- if (length(.nz_vals) >= 3L) {
  .old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set_fixed_seed(RANDOM_SEED_BASE)
  .bs <- vapply(seq_len(2000L),
                function(i) mean(.nz_vals[sample.int(length(.nz_vals), length(.nz_vals), TRUE)]),
                numeric(1))
  if (!is.null(.old_seed)) assign(".Random.seed", .old_seed, envir = globalenv())
  unname(stats::quantile(.bs, c(0.05, 0.95), names = FALSE))
} else c(NA_real_, NA_real_)

#' The geography (named R2 vector) with the largest finite R2; exact ties are named together.
.best_geo <- function(r2) {
  r2 <- r2[is.finite(r2)]
  if (!length(r2)) return("not estimable")
  paste(names(r2)[r2 == max(r2)], collapse = " = ")
}

# How stable the envelope verdict is across the onset completions (primary included, as in imp_q()).
VERDICT_STABILITY <- if (is.null(imp) || !nrow(imp)) "" else
  sprintf(" (the same verdict class in %d of %d onset completions, primary included)",
          1L + sum(.verdict_class(imp$envelope_verdict) == .verdict_class(S$envelope_verdict)),
          nrow(imp) + 1L)

summary_out <- tibble::tribble(
  ~statistic, ~value, ~lo, ~hi, ~unit, ~basis,
  "Radial expansion rate (road)",         S$v_radial_road, ci_road[1], ci_road[2], "km/week",
    "OLS of road distance from the index zone on invasion week; lo/hi = 90% zone bootstrap",
  "Radial expansion rate (great-circle)", S$v_radial_gc, ci_gc[1], ci_gc[2], "km/week",
    "OLS of great-circle distance from the index zone on invasion week; lo/hi = 90% zone bootstrap",
  "Radial expansion rate (road, Theil-Sen)", S$v_radial_road_ts, NA_real_, NA_real_, "km/week",
    "median of pairwise slopes — outlier-resistant check on the OLS road rate above",
  "Radial expansion rate (great-circle, Theil-Sen)", S$v_radial_gc_ts, NA_real_, NA_real_, "km/week",
    "median of pairwise slopes — outlier-resistant check on the OLS great-circle rate",
  "Radial expansion rate (travel time)", S$v_radial_time, NA_real_, NA_real_, "hours/week",
    # WHICH geography fits best is DERIVED from the three R2, never typed: the sentence used to
    # assert travel time won whatever the numbers beside it said.
    sprintf("OLS of OSRM travel time on invasion week; R2 = %.2f vs %.2f road, %.2f great-circle (best-fitting geography of the three: %s)",
            S$v_radial_time_r2, S$v_radial_road_r2, S$v_radial_gc_r2,
            .best_geo(c(`travel time` = S$v_radial_time_r2, `road distance` = S$v_radial_road_r2,
                        `great-circle` = S$v_radial_gc_r2))),
  "Median-reach growth rate (road)",      S$v_reach_med, NA_real_, NA_real_, "km/week",
    "OLS of the MEDIAN distance of invaded zones on week (bulk expansion of the invaded set)",
  "Case-weighted mean distance trend",    S$v_case_mean, NA_real_, NA_real_, "km/week",
    "OLS of the case-weighted mean distance from the index zone on week",
  # THE ENVELOPE IS THE p90 REACH and the verdict is the DATA's (compute_kinematics(): the
  # running maximum let one distant zone set the threshold week, the "to 90%" speed and a
  # post-threshold slope of zero by construction, and the word SATURATED was typed here).
  "Envelope (p90 reach): week reaching 90% of final value", S$week_p90_reach_90pct, NA_real_, NA_real_, "week",
    "first fit-window week in which the 90th-percentile road reach of the invaded set, R_p90(t), attains 0.9x its fit-window final value (one distant zone cannot set it, unlike the running maximum)",
  "Envelope (p90 reach) slope, up to that week", S$v_env_p90_early, NA_real_, NA_real_, "km/week",
    sprintf("OLS of R_p90 on week over the %d fit-window week(s) up to and including the threshold week: how fast the envelope was established",
            S$n_weeks_env_early),
  "Envelope (p90 reach) slope, from that week", S$v_env_p90_late, NA_real_, NA_real_, "km/week",
    sprintf("OLS of R_p90 on week over the %d fit-window week(s) from the threshold week on. Verdict (post-threshold slope at most a tenth of the establishing slope = saturated): %s%s. A saturated envelope does NOT mean long-range seeding stopped: see the two rows below",
            S$n_weeks_env_late, S$envelope_verdict, VERDICT_STABILITY),
  "New invasions after the envelope threshold, median furthest", S$newmax_med_after, NA_real_, NA_real_, "km",
    sprintf("median, over the %d fit-window weeks after the p90 threshold week with a new invasion, of the FURTHEST zone invaded that week — the long-range seeding an envelope does not show",
            S$n_weeks_after),
  "New invasions after the envelope threshold, furthest", S$newmax_max_after, NA_real_, NA_real_, "km",
    "single furthest zone invaded in any fit-window week after the p90 threshold week",
  "Leading-edge slope, whole window (road)", S$v_front_road, NA_real_, NA_real_, "km/week",
    "OLS of R_max(t) over the whole fit window; a running max, so this is a window AVERAGE, not a wavefront speed — read it with the two rows above",
  "Leading-edge slope, p90 (road)",       S$v_front_p90_road, NA_real_, NA_real_, "km/week",
    "OLS of the 90th-percentile reach on week (outlier-robust envelope)",
  "Case centre-of-gravity speed",         S$v_cog_med, S$v_cog_iqr_lo, S$v_cog_iqr_hi, "km/week",
    "median (IQR) of weekly geodesic centre-of-gravity displacement",
  "New-invasion centre-of-gravity speed", S$v_invcog_med, NA_real_, NA_real_, "km/week",
    "median weekly displacement of the mean centre of newly invaded zones; in weeks with only 1-2 new zones (see new_zones in the weekly CSV) this is a jump between individual zones, not a centre-of-mass velocity",
  "Case centre-of-gravity net speed",     S$v_cog_net, NA_real_, NA_real_, "km/week",
    "net first-to-last displacement over the fit window divided by its length",
  "Case centre-of-gravity net displacement", S$cog_net_disp_km, NA_real_, NA_real_, "km",
    "geodesic first-to-last centre of gravity, fit window",
  "Case centre-of-gravity path length",   S$cog_path_km, NA_real_, NA_real_, "km",
    "summed weekly steps of the centre of gravity, fit window",
  "Case centre-of-gravity bearing",       S$cog_net_bearing_deg, NA_real_, NA_real_, "deg (compass)",
    "net direction of centre-of-gravity travel",
  "Case centre-of-gravity straightness",  S$cog_straightness, NA_real_, NA_real_, "0-1",
    "net displacement / path length (1 = perfectly directional)",
  "Dispersion growth rate",               S$v_disp_km_wk, NA_real_, NA_real_, "km/week",
    "OLS of the weekly case standard distance on week",
  "Effective diffusivity",                S$diffusion_km2_wk, NA_real_, NA_real_, "km^2/week",
    # The R2 is STATED, not interpreted: the sentence here used to call the row "a null result"
    # whatever the R2 was (it printed R2 = 0.87 beside it).
    sprintf("slope(SD^2 ~ t)/4, the isotropic-diffusion reading of the dispersion growth; SD^2 ~ t fit R2 = %.2f (D is meaningful only if SD^2 grows linearly in t)",
            S$diffusion_r2),
  "Invaded-zone accrual rate",            S$accrual_rate_wk, NA_real_, NA_real_, "log-zones/week",
    "OLS of log(cumulative invaded zones) on week",
  # PUBLISHED because a figure was computing it. make_spread_kinematics_compact.R drew Panel
  # C's reference line from its own mean of new_zones over the same weeks, and its comment
  # asserted the interval "is in spread_kinematics_summary.csv" — there was no such row. The
  # window is weeks 1..fit_max: week 0 is the seeding condition (those zones already had
  # confirmed cases when the series opens), and the truncated tail is incompletely observed.
  "Mean new invasions per week",          .nz_mean, .nz_ci[1], .nz_ci[2], "zones/week",
    sprintf("mean of new_zones over weeks 1-%d (week 0 = seeding, truncated tail excluded); lo/hi = 90%% week-bootstrap interval, weeks resampled as exchangeable (the series trends and is serially correlated, so this is the week-to-week spread around the mean, not a CI for a stationary rate)",
            .nz_fit_max),
  "Invaded-zone doubling time",           S$accrual_doubling_wk, NA_real_, NA_real_, "weeks",
    "log(2) / accrual rate",
  "Reach at the end of the fit window",   S$reach_max_km, NA_real_, NA_real_, "km",
    "furthest invaded zone by road from the index zone",
  "Median reach at the end of the fit window", S$reach_med_km, NA_real_, NA_real_, "km",
    "median road distance of invaded zones from the index zone",
  "Zones invaded (fit window)",           S$n_zones_fit, NA_real_, NA_real_, "zones",
    "zones with >=1 confirmed case by the end of the fit window",
  "Zones invaded (all weeks)",            S$n_zones_all, NA_real_, NA_real_, "zones",
    "includes the right-truncated tail")

# Imputation-uncertainty columns: 5th / 50th / 95th percentile across the onset COMPLETIONS
# (the primary counted as replicate 0; see imp_q()). `value` is the PRIMARY completion — the
# one every other module scores — not the median, and `value_basis` says so on every row;
# imp_med is published beside it so the two can be compared directly. EVERY row has an entry
# except the centre-of-gravity bearing: an angle, whose linear percentiles are meaningless
# near the 0/360 wrap. The stopifnot keeps a new summary row from silently going without one.
imp_map <- c("Radial expansion rate (road)" = "v_radial_road",
             "Radial expansion rate (great-circle)" = "v_radial_gc",
             "Radial expansion rate (road, Theil-Sen)" = "v_radial_road_ts",
             "Radial expansion rate (great-circle, Theil-Sen)" = "v_radial_gc_ts",
             "Radial expansion rate (travel time)" = "v_radial_time",
             "Median-reach growth rate (road)" = "v_reach_med",
             "Case-weighted mean distance trend" = "v_case_mean",
             "Envelope (p90 reach): week reaching 90% of final value" = "week_p90_reach_90pct",
             "Envelope (p90 reach) slope, up to that week" = "v_env_p90_early",
             "Envelope (p90 reach) slope, from that week" = "v_env_p90_late",
             "New invasions after the envelope threshold, median furthest" = "newmax_med_after",
             "New invasions after the envelope threshold, furthest" = "newmax_max_after",
             "Leading-edge slope, whole window (road)" = "v_front_road",
             "Leading-edge slope, p90 (road)" = "v_front_p90_road",
             "Case centre-of-gravity speed" = "v_cog_med",
             "New-invasion centre-of-gravity speed" = "v_invcog_med",
             "Case centre-of-gravity net speed" = "v_cog_net",
             "Case centre-of-gravity net displacement" = "cog_net_disp_km",
             "Case centre-of-gravity path length" = "cog_path_km",
             "Case centre-of-gravity straightness" = "cog_straightness",
             "Dispersion growth rate" = "v_disp_km_wk",
             "Effective diffusivity" = "diffusion_km2_wk",
             "Invaded-zone accrual rate" = "accrual_rate_wk",
             "Mean new invasions per week" = "nz_mean",
             "Invaded-zone doubling time" = "accrual_doubling_wk",
             "Reach at the end of the fit window" = "reach_max_km",
             "Median reach at the end of the fit window" = "reach_med_km",
             "Zones invaded (fit window)" = "n_zones_fit",
             "Zones invaded (all weeks)" = "n_zones_all")
stopifnot(setequal(c(names(imp_map), "Case centre-of-gravity bearing"), summary_out$statistic))
iq <- t(vapply(summary_out$statistic, function(s)
  if (s %in% names(imp_map)) imp_q(unname(imp_map[s])) else c(NA_real_, NA_real_, NA_real_),
  numeric(3)))
.n_compl <- if (is.null(imp) || !nrow(imp)) 1L else nrow(imp) + 1L
summary_out <- summary_out %>%
  dplyr::mutate(value = round(value, 4), lo = round(lo, 4), hi = round(hi, 4),
                imp_lo = round(iq[, 1], 4), imp_med = round(iq[, 2], 4),
                imp_hi = round(iq[, 3], 4),
                value_basis = sprintf("primary onset completion (ONSET_IMPUTE_SEED %d)%s",
                                      ONSET_SEED_BASE,
                                      ifelse(statistic %in% names(imp_map),
                                             sprintf("; imp_lo/imp_med/imp_hi = 5th/50th/95th percentile over %d completions (primary + %d replicates)",
                                                     .n_compl, .n_compl - 1L), ""))) %>%
  dplyr::relocate(imp_lo, imp_med, imp_hi, .after = hi)
readr::write_csv(summary_out, file.path(KEY_DIR, "spread_kinematics_summary.csv"))

diag_out <- dplyr::bind_rows(
  K$metric_fit %>%
    dplyr::transmute(diagnostic = "arrival-time ordering by distance metric",
                     key = metric, spearman = round(spearman, 4), r2 = round(r2, 4),
                     value = round(slope_per_week, 4), note = "slope = units per week"),
  sens %>%
    dplyr::transmute(diagnostic = "truncation sensitivity",
                     key = sprintf("drop %d trailing week(s)%s", trunc_weeks,
                                   ifelse(primary, " [primary]", "")),
                     spearman = NA_real_, r2 = NA_real_,
                     value = round(v_radial_road, 3),
                     note = sprintf("v_radial=%.1f v_front_rmax_ols=%.1f v_cog=%.1f v_disp=%.1f n_zones=%d",
                                    v_radial_road, v_front_road, v_cog_med, v_disp_km_wk,
                                    n_zones_fit)),
  orig_sens %>%
    dplyr::transmute(diagnostic = "index-zone sensitivity",
                     key = sprintf("%s (%d seed cases)%s", origin, seed_cases,
                                   ifelse(primary, " [primary]", "")),
                     spearman = NA_real_, r2 = round(v_radial_road_r2, 4),
                     value = round(v_radial_road, 3),
                     note = sprintf("v_radial=%.1f v_reach_med=%.1f max=%.0f km med=%.0f km",
                                    v_radial_road, v_reach_med, reach_max_km, reach_med_km)),
  tibble::tibble(diagnostic = "sitrep-only zone leverage",
                 key = sprintf("drop %d sitrep-only zone(s)%s", nosit_n,
                               if (nosit_n) paste0(": ", paste(
                                 intersect(fitz_primary$health_zone, sitrep_only_zones),
                                 collapse = ", ")) else ""),
                 spearman = NA_real_, r2 = round(nosit_road$r2, 4),
                 value = round(nosit_road$slope, 3),
                 note = sprintf("v_radial %.1f -> %.1f km/wk on n=%d (late-biased invasion dates)",
                                S$v_radial_road, nosit_road$slope, nosit_road$n)),
  tibble::tibble(diagnostic = "onset-week completeness (delay ECDF)",
                 key = as.character(weeks_all), spearman = NA_real_, r2 = NA_real_,
                 value = round(comp, 4),
                 note = ifelse(seq_along(weeks_all) > length(weeks_all) - n_trunc,
                               "dropped from fits", "used in fits")))
readr::write_csv(diag_out, file.path(KEY_DIR, "spread_kinematics_diagnostics.csv"))

# ---- console summary --------------------------------------------------------
message("\n", strrep("-", 74))
message(sprintf("  index zone (epicentre)      %s (%s)", K$epi_zone, K$epi_province))
message(sprintf("  fit window                  weeks 0-%d of 0-%d (last %d dropped: truncation)",
                K$fit_max, K$week_max, K$n_trunc))
message(sprintf("  radial expansion  (road)    %.1f km/wk  [90%% zone bootstrap %.1f, %.1f]",
                S$v_radial_road, ci_road[1], ci_road[2]))
message(sprintf("  radial expansion  (gc)      %.1f km/wk  [90%% zone bootstrap %.1f, %.1f]",
                S$v_radial_gc, ci_gc[1], ci_gc[2]))
message(sprintf("  robustness of that rate     Theil-Sen %.1f (road) / %.1f (gc) km/wk; %.1f km/wk dropping %d sitrep-only zone(s)",
                S$v_radial_road_ts, S$v_radial_gc_ts, nosit_road$slope, nosit_n))
message(sprintf("  median-reach growth (road)  %.1f km/wk", S$v_reach_med))
# The verdict is the one compute_kinematics() derived from the data (and wrote to the CSV);
# it is not re-decided, and not typed, here.
message(sprintf("  envelope (p90 reach)        reached 90%% of its final value by week %.0f (%+.1f km/wk to then over %d wk, %+.1f km/wk after over %d wk) -> %s",
                S$week_p90_reach_90pct, S$v_env_p90_early, S$n_weeks_env_early,
                S$v_env_p90_late, S$n_weeks_env_late, S$envelope_verdict))
message(sprintf("  new invasions after that    in the %d later fit-window week(s) with a new invasion, the furthest new zone was a median %.0f km out (max %.0f km)",
                S$n_weeks_after, S$newmax_med_after, S$newmax_max_after))
message(sprintf("  case CoG speed              %.1f km/wk (IQR %.1f-%.1f); net %.1f km/wk, bearing %.0f deg, straightness %.2f",
                S$v_cog_med, S$v_cog_iqr_lo, S$v_cog_iqr_hi, S$v_cog_net,
                S$cog_net_bearing_deg, S$cog_straightness))
message(sprintf("  new-invasion CoG speed      %.1f km/wk", S$v_invcog_med))
message(sprintf("  dispersion growth           %.1f km/wk  (effective diffusivity %.0f km2/wk)",
                S$v_disp_km_wk, S$diffusion_km2_wk))
message(sprintf("  invaded-zone doubling       %.1f weeks (%d zones by week %d)",
                S$accrual_doubling_wk, S$n_zones_fit, K$fit_max))
message(sprintf("  reach from index zone       max %.0f km, median %.0f km (road)",
                S$reach_max_km, S$reach_med_km))
message("  arrival-time ordering by metric:")
for (i in seq_len(nrow(K$metric_fit)))
  message(sprintf("     %-20s spearman %+.3f   R2 %.3f   slope %.2f/wk",
                  K$metric_fit$metric[i], K$metric_fit$spearman[i],
                  K$metric_fit$r2[i], K$metric_fit$slope_per_week[i]))
message("  index-zone sensitivity (top seed-window candidates):")
for (i in seq_len(nrow(orig_sens)))
  message(sprintf("     %-12s (%2d seed cases)  v_radial %.1f km/wk  max reach %.0f km%s",
                  orig_sens$origin[i], orig_sens$seed_cases[i], orig_sens$v_radial_road[i],
                  orig_sens$reach_max_km[i], ifelse(orig_sens$primary[i], "  [primary]", "")))
if (!is.null(imp) && nrow(imp)) {
  tb <- sort(table(imp$epi_zone), decreasing = TRUE)
  message(sprintf("  onset-imputation replicates %d; index zone stable in %d%% (modal %s)",
                  nrow(imp), round(100 * tb[1] / sum(tb)), names(tb)[1]))
  # median [5th, 95th] over the completions, primary included as replicate 0 (imp_q())
  message(sprintf("     v_radial %.1f [%.1f, %.1f]   v_cog %.1f [%.1f, %.1f]   max reach %.0f [%.0f, %.0f] km",
                  imp_q("v_radial_road")[2], imp_q("v_radial_road")[1], imp_q("v_radial_road")[3],
                  imp_q("v_cog_med")[2], imp_q("v_cog_med")[1], imp_q("v_cog_med")[3],
                  imp_q("reach_max_km")[2], imp_q("reach_max_km")[1], imp_q("reach_max_km")[3]))
  message(sprintf("     envelope verdict %s%s", .verdict_class(S$envelope_verdict), VERDICT_STABILITY))
  # A guard, not an expected event: load_linelist() floors date_index at OUTBREAK_START and the
  # grid ends on the analysis week, so on the current line list every completion has the same
  # grid. If that ever changes, the calendar-anchored fit window keeps the runs comparable.
  if (dplyr::n_distinct(imp$n_weeks) > 1L)
    message(sprintf("     note: onset draws produced %d different week-grid lengths (%s); the fit window is calendar-anchored, so the runs stay comparable",
                    dplyr::n_distinct(imp$n_weeks),
                    paste(sort(unique(imp$n_weeks)), collapse = "/")))
}
message(strrep("-", 74), "\n")

# =============================================================================
# SECTION 6 — figure (REMOVED)
# =============================================================================
# This module used to assemble a five-panel `Figure_spread_kinematics`. That stem is not
# on FIGURE_KEEP, so the whole figure — five panel builders and an inset — was computed on
# every run and then discarded by the gate. The RETAINED kinematics figure is
# `Figure_spread_kinematics_compact`, which make_spread_kinematics_compact.R builds from
# the CSVs written above; Figures 1C and S1 likewise read arrival_predictor{s,_fits}.csv.
# This module is therefore a pure ESTIMATION step and writes no figure of its own.
message("[kinematics] done — CSVs in ", KEY_DIR)

}   # end of the script-run guard opened at SECTION 4
