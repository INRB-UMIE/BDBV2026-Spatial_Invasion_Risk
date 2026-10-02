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
#         v_front  : the LEADING EDGE R_max(t) = furthest zone invaded by week t,
#                    split into the rate at which the extent was ESTABLISHED and the
#                    slope AFTER 90% of it was reached (~0 = saturated, no wavefront),
#                    plus the whole-window OLS slope for reference only
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
#       however distant — a fresh 767 km seeding into a new province hides completely
#       behind an earlier 794 km one in a different direction. Reporting only the
#       running max would say "the extent saturated" and let a reader conclude that
#       long-range spread had stopped, when it demonstrably had not.
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
#     from EVERY fit and drawn hollow in the figure. Sensitivity to that choice
#     (1/2/3 weeks) is reported in the diagnostics CSV.
#   * ONSET IMPUTATION. ~15% of confirmed records have no onset date and are dated
#     by a stochastic backward draw (00_config.R ONSET_MODE). The primary run uses
#     the project's own RANDOM_SEED, so its invasion dates are exactly the ones the
#     rest of the suite scores; N_IMPUTE replicates then VARY that seed (a local
#     set.seed() would be ignored — load_linelist() reseeds itself) and re-run the
#     whole pipeline, so every headline statistic carries an imputation interval.
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
#   - The outbreak seeded Nord-Kivu within the first epidemic week, so the invaded
#     set is NOT a single expanding wave; the radial statistics summarise a
#     multi-focus process and should be read as such (Panel C makes this visible).
#
# Run:  Rscript spatiotemporal/43_spread_kinematics.R
# Env:  SPREAD_SEED, SPREAD_TRUNC_WEEKS, SPREAD_IMPUTE_REPS, SPREAD_BOOT
# Out:  outputs/key_outputs/Figure_spread_kinematics.{pdf,png}
#       outputs/key_outputs/spread_kinematics_{weekly,zones,summary,diagnostics}.csv
# Panel titles/subtitles/captions are intentionally OMITTED (house style): figure
# identity is carried by the A-E tags and axis titles; the caption lives externally.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(sf); library(scales); library(geosphere)
  library(ggrepel); library(here)
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
# from the project-wide RANDOM_SEED that load_linelist() reads (see Section 4).
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
# Geodesic (straight-line) rings drawn on the map. NOTE they are NOT road distance —
# every numeric reach statistic is road distance — so the innermost ring is labelled
# "straight-line" to keep the two metrics from being read as one.
RINGS_KM       <- c(200, 400, 600)
TRUNC_SENS     <- 1:3       # trailing-week drops compared in the sensitivity table

KEY_DIR <- file.path(OUT_DIR, "key_outputs")
FIG_DIR <- file.path(KEY_DIR, "figures")
for (d in c(KEY_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Style (verbatim from 40_cascade_next_dominoes.R / make_publication_figures.R)
# ---------------------------------------------------------------------------
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
# Province identity, VERBATIM from 37/40_* so a province keeps the same hue across every
# figure in the suite. It deliberately re-uses the hues the series palette below also
# draws on: Panel C is the only panel with a province scale and it carries its own
# legend, and cross-figure province identity matters more than intra-figure hue purity.
# Do not "fix" this by recolouring provinces — that would desynchronise Figures 2 and 4.
PROV_COL <- c("Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73",
              "Sud-Kivu"="#CC79A7","Tshopo"="#E69F00","Bas-Uele"="#56B4E9")
SPARE_COL <- c("#882255","#332288","#661100","#999933","#44AA99","#AA4499","#DDCC77","#117733")
# ONE semantic colour scheme across every panel. The panels share a vocabulary, so the
# colours must too: a reader who learns "orange = cases, green = new invasions" from the
# map must not then meet orange standing for "furthest invaded zone" in a line panel —
# they would read that panel exactly backwards. Every series below is coloured by WHAT
# IT MEASURES, never by which panel it happens to sit in.
CASE_COL <- "#D55E00"   # anything CASE-weighted (that week's confirmed cases)
INV_COL  <- "#009E73"   # anything about the zones NEWLY invaded that week
EXT_COL  <- "#0072B2"   # the CUMULATIVE invaded set — its extent / envelope
EXT_COL2 <- "#56B4E9"   # lighter shade, for a second cumulative-set series in one panel
EPI_COL  <- "#7A0177"   # epicentre marker
base_family <- "sans"

theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    axis.title = element_text(size = base - 0.4, colour = MUTED),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4), angle = 90),
    axis.text = element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
    strip.text = element_text(size = base - 0.6, colour = INK, face = "bold"),
    legend.position = "top", legend.justification = "left",
    legend.title = element_text(size = base - 1.2, colour = MUTED),
    legend.text = element_text(size = base - 1.4, colour = INK),
    legend.key.height = unit(9, "pt"), legend.key.width = unit(15, "pt"),
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(3, "pt"),
    plot.tag = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin = margin(6, 8, 6, 6))
}
theme_map <- function(base = 8.6) {
  theme_void(base_size = base, base_family = base_family) %+replace% theme(
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    legend.position = "right", legend.justification = "centre",
    legend.title = element_text(size = base - 1.4, colour = MUTED),
    legend.text = element_text(size = base - 1.8, colour = INK),
    legend.key.height = unit(14, "pt"), legend.key.width = unit(9, "pt"),
    plot.tag = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin = margin(2, 2, 2, 2))
}
save_dual <- function(p, name, w, h, dirs = c(FIG_DIR, KEY_DIR)) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  # Gate on each destination path (these helpers write into several directories), so
  # FIGURE_DROP's "<dir>/<stem>" entries and the raw/ rule can actually match.
  if (is.function(.fk) && !any(vapply(dirs, function(.d) .fk(file.path(.d, name)),
                                      logical(1)))) return(invisible(p))
  for (d in dirs) {
    ggsave(file.path(d, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
    ggsave(file.path(d, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  }
  message(sprintf("  saved %-34s %.1f x %.1f in", name, w, h)); invisible(p)
}

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
#' completeness). Delays are capped at `cap` days; longer delays (rare: the 95th
#' percentile is ~24 d) are therefore treated as if they had arrived, so this is a
#' slight UPPER bound on completeness — conservative about how few weeks it flags,
#' which is why TRUNC_MIN also applies. Sitrep top-up rows carry no lab/report date
#' and drop out of the ECDF automatically.
week_completeness <- function(ll, weeks, asof, cap = DELAY_CAP_DAYS) {
  cf <- ll %>% dplyr::filter(confirmed %in% TRUE, !is.na(date_index))
  known <- dplyr::coalesce(as.Date(cf$lab_analysis_date),
                           as.Date(cf$reporting_date),
                           as.Date(cf$date_of_notification))
  del <- as.numeric(known - cf$date_index)
  ok  <- is.finite(del) & del >= 0 & del <= cap & cf$date_index <= (asof - cap)
  if (sum(ok) < 50L) {
    warning("[kinematics] only ", sum(ok), " complete-follow-up delay pairs — ",
            "completeness not estimated; falling back to TRUNC_MIN.", call. = FALSE)
    return(rep(NA_real_, length(weeks)))
  }
  Fhat <- stats::ecdf(del[ok])
  vapply(weeks, function(w) {
    avail <- as.numeric(asof - (w + 0:6))              # follow-up days per onset day
    mean(ifelse(avail < 0, 0, Fhat(pmin(avail, cap)))) # onset days after `asof` contribute 0
  }, numeric(1))
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
#' vector: an onset-imputation replicate can draw an onset far enough back to add an
#' earlier week to the grid, which would shift every positional index and make the
#' replicates incomparable with each other and with the primary run.
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
      # invisible in them however distant it is — a brand-new 767 km seeding into a new
      # province hides completely behind an earlier 794 km one in another direction.
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

  # SATURATION of the leading edge. R_max(t) is a running maximum: if the epidemic
  # reaches its full extent early by long-range seeding and then stops expanding, a
  # single OLS slope over the whole window is NOT a wavefront speed. Split it:
  # how quickly the envelope reached 90% of its final value, and whether it still
  # grows afterwards.
  rmax_fin <- fw$reach_max_km[nrow(fw)]
  i90 <- which(is.finite(fw$reach_max_km) & fw$reach_max_km >= 0.9 * rmax_fin)
  wk90 <- if (length(i90)) fw$week[i90[1]] else NA_real_
  v_front_early <- if (is.finite(wk90) && wk90 > 0) fw$reach_max_km[i90[1]] / wk90 else NA_real_
  late <- fw %>% dplyr::filter(is.finite(wk90), week >= wk90)
  v_front_late <- .ols(late$week, late$reach_max_km)
  # A saturated envelope does NOT mean long-range introductions stopped — only that
  # none exceeded the standing record. Quantify what kept arriving after saturation.
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
    front_intercept = v_front$intercept, radial_intercept = v_radial_road$intercept,
    disp_intercept = v_disp$intercept, reach_med_intercept = v_reach_med$intercept,
    stats = list(
      v_radial_road = v_radial_road$slope, v_radial_road_r2 = v_radial_road$r2,
      v_radial_gc   = v_radial_gc$slope,   v_radial_gc_r2   = v_radial_gc$r2,
      v_radial_time = v_radial_time$slope, v_radial_time_r2 = v_radial_time$r2,
      v_radial_road_ts = .theil_sen(fitz$inv_week, fitz$d_road),
      v_radial_gc_ts   = .theil_sen(fitz$inv_week, fitz$d_gc),
      v_front_road  = v_front$slope,       v_front_p90_road = v_front_p90$slope,
      v_front_early = v_front_early,       v_front_late = v_front_late$slope,
      week_reach_90pct = wk90,
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

# load_linelist() seeds its own onset draw from the GLOBAL `RANDOM_SEED`
# (01_data_prep.R) and restores the caller's RNG stream on exit. Two consequences:
#   1. the primary run needs no set.seed() here and is byte-identical to the line
#      list every other module in the suite sees — the invasion dates below are the
#      same events the forecasting pipeline scores;
#   2. imputation replicates must vary RANDOM_SEED ITSELF. A set.seed() around the
#      call is silently ignored and would yield N identical "replicates", i.e. a
#      zero-width uncertainty interval that looks like precision.
RANDOM_SEED_BASE <- as.integer(get0("RANDOM_SEED", ifnotfound = 20260704L))
ll_primary <- load_linelist()
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
boot_ci <- function(z, col, B = N_BOOT) {
  dd <- z[[col]]; tt <- z$inv_week
  ok <- is.finite(dd) & is.finite(tt); dd <- dd[ok]; tt <- tt[ok]; n <- length(dd)
  if (n < 5L) return(c(NA_real_, NA_real_))
  s <- vapply(seq_len(B), function(b) {
    i <- sample.int(n, n, replace = TRUE)
    if (length(unique(tt[i])) < 2L) return(NA_real_)
    unname(stats::coef(stats::lm(dd[i] ~ tt[i]))[2])
  }, numeric(1))
  c(.q(s, 0.025), .q(s, 0.975))
}
set.seed(SEED + 1L)
fitz_primary <- dplyr::filter(K$zones, in_fit_window)
ci_road <- boot_ci(fitz_primary, "d_road")
ci_gc   <- boot_ci(fitz_primary, "d_gc")

# ---- sensitivity to the index-zone choice ----------------------------------
# The top seed-window zones are separated by only a case or two AND all sit in the
# same urban focus, so the argmax rule is fragile by construction. Re-running the
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

# ---- onset-imputation replicates -------------------------------------------
# The whole pipeline (including the index-zone choice) is re-run on independent
# onset draws, so every headline number carries an imputation-uncertainty interval.
imp <- NULL
if (N_IMPUTE > 0L) {
  message(sprintf("[kinematics] %d onset-imputation replicate(s)...", N_IMPUTE))
  imp <- dplyr::bind_rows(lapply(seq_len(N_IMPUTE), function(r) {
    # Vary the seed load_linelist() actually reads, and put it back afterwards so the
    # rest of the session still sees the project default.
    .old <- get0("RANDOM_SEED", envir = globalenv(), ifnotfound = NULL)
    assign("RANDOM_SEED", RANDOM_SEED_BASE + r, envir = globalenv())
    on.exit(if (is.null(.old)) suppressWarnings(rm("RANDOM_SEED", envir = globalenv()))
            else assign("RANDOM_SEED", .old, envir = globalenv()), add = TRUE)
    out <- try({
      llr <- suppressMessages(suppressWarnings(load_linelist()))
      zwr <- suppressMessages(suppressWarnings(aggregate_to_zone_week(llr, names(pop))))
      Kr  <- compute_kinematics(zwr, geo, fit_cutoff = FIT_CUTOFF,
                                seed_cutoff = SEED_CUTOFF)
      tibble::tibble(rep = r, epi_zone = Kr$epi_zone, n_weeks = length(Kr$weeks),
                     # per-zone first confirmed onset on THIS draw, carried as a list
                     # column so Section 5b can put an imputation interval on the
                     # arrival-vs-predictor correlations too (the predictors are fixed
                     # geography and mobility; only the arrival dates are imputed).
                     arrivals = list(llr %>%
                       dplyr::filter(confirmed %in% TRUE, !is.na(date_index), !is.na(health_zone)) %>%
                       dplyr::group_by(health_zone) %>%
                       dplyr::summarise(first_onset = min(date_index), .groups = "drop")),
                     v_radial_road = Kr$stats$v_radial_road,
                     v_radial_gc = Kr$stats$v_radial_gc,
                     v_reach_med = Kr$stats$v_reach_med,
                     v_cog_med = Kr$stats$v_cog_med,
                     v_cog_net = Kr$stats$v_cog_net,
                     v_disp_km_wk = Kr$stats$v_disp_km_wk,
                     accrual_doubling_wk = Kr$stats$accrual_doubling_wk,
                     reach_max_km = Kr$stats$reach_max_km,
                     n_zones_fit = Kr$stats$n_zones_fit)
    }, silent = TRUE)
    if (inherits(out, "try-error")) NULL else out
  }))
  if (nrow(imp) < N_IMPUTE)
    warning("[kinematics] ", N_IMPUTE - nrow(imp), " imputation replicate(s) failed.",
            call. = FALSE)
}
imp_q <- function(col) {
  if (is.null(imp) || !nrow(imp)) return(c(NA_real_, NA_real_))
  c(.q(imp[[col]], 0.05), .q(imp[[col]], 0.95))
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

zones_out <- K$zones %>%
  dplyr::transmute(health_zone, province,
                   invasion_week = inv_week, invasion_week_start = inv_week_start,
                   cases_total,
                   dist_greatcircle_km = round(d_gc, 2),
                   dist_road_km = round(d_road, 2),
                   travel_time_h = round(t_road, 3),
                   bearing_from_epicentre_deg = round(bearing, 1),
                   lon = round(lon, 5), lat = round(lat, 5),
                   in_fit_window,
                   sitrep_only = health_zone %in% sitrep_only_zones)
readr::write_csv(zones_out, file.path(KEY_DIR, "spread_kinematics_zones.csv"))

# -----------------------------------------------------------------------------
# SECTION 5b — arrival time vs candidate predictors
#   publishes: arrival_predictors.csv, arrival_predictor_fits.csv
#   consumed by: make_manuscript_figures.R (Figure 1C and Figure S1)
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
# remain available in spread_kinematics_zones.csv, so the two are comparable rather
# than conflated.
#
# ARRIVAL. Days from ARRIVAL_ORIGIN, the earliest confirmed date_index in the
# reconciled line list — an observed event, not a chosen calendar date. date_index
# is load_linelist()'s, i.e. the same onset draw the forecasting pipeline scores,
# including the sitrep top-up.
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

  # ---- March relocation OD share out of the epicentre region ----------------
  # NA in the OD matrix means "no count released for this pair this month", not zero
  # (see load_flowminder_od). A destination for which EVERY epicentre row is NA has
  # an unobserved share and is carried as NA; where some rows are released, the
  # observed ones are summed. Shares are normalised over the observed destinations.
  M_od <- load_flowminder_od()
  rownames(M_od) <- .can(rownames(M_od)); colnames(M_od) <- .can(colnames(M_od))
  ei <- intersect(epi_set, rownames(M_od))
  if (!length(ei)) {
    warning("[kinematics] no epicentre zone present in the Flowminder OD matrix; ",
            "the March relocation predictor is unavailable.", call. = FALSE)
    od_share_v <- setNames(rep(NA_real_, nrow(ctr)), ctr$health_zone)
  } else {
    sub  <- M_od[ei, , drop = FALSE]
    n_ob <- colSums(!is.na(sub))
    raw  <- colSums(sub, na.rm = TRUE)
    raw[n_ob == 0L] <- NA_real_
    raw <- .collapse(raw)
    od_share_v <- raw / sum(raw, na.rm = TRUE)
  }

  # ---- latest short-trip snapshot (pooled M+B+R cohort) ---------------------
  # load_short_trip_snapshot() already returns destination PROPORTIONS for the pooled
  # epicentre cohort, so there is no epicentre subsetting to do here; it is
  # renormalised because the released percentages need not sum to exactly 100.
  st_tag <- FLOWMINDER_ST_TAGS[length(FLOWMINDER_ST_TAGS)]
  st_v   <- load_short_trip_snapshot(st_tag)
  names(st_v) <- .can(names(st_v))
  st_v <- .collapse(st_v)
  st_v <- st_v / sum(st_v, na.rm = TRUE)

  .lk <- function(v, z) { out <- as.numeric(v[z]); out[is.na(match(z, names(v)))] <- NA_real_; out }
  .neglog <- function(p) ifelse(is.finite(p) & p > 0, -log(p), NA_real_)

  out <- arr %>%
    dplyr::mutate(
      province            = ctr$province[match(health_zone, ctr$health_zone)],
      arrival_origin      = arrival_origin,
      arrival_days        = as.numeric(first_onset - arrival_origin),
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

  attr(out, "epi_set")        <- epi_set
  attr(out, "shorttrip_tag")  <- st_tag
  attr(out, "arrival_origin") <- arrival_origin
  out
})

ARRIVAL_ORIGIN <- attr(arrival_predictors, "arrival_origin")
message(sprintf("[kinematics] arrival-predictor table: %d zone(s), origin %s, epicentre = %s, short-trip snapshot %s",
                nrow(arrival_predictors), format(ARRIVAL_ORIGIN, "%Y-%m-%d"),
                paste(attr(arrival_predictors, "epi_set"), collapse = "+"),
                attr(arrival_predictors, "shorttrip_tag")))

# Written at FULL WORKING PRECISION, not rounded for reading. The fits below are the
# panel's published statistics and a reader must be able to recompute them from this
# table exactly; rounding the predictors to 2-4 decimals for legibility moved the
# recomputed Spearman rho by ~5e-4, which is enough to disagree with the printed
# 2-decimal annotation and to make the panel unverifiable. signif(, 12) is lossless at
# these magnitudes and still keeps the file readable.
readr::write_csv(
  arrival_predictors %>%
    dplyr::transmute(health_zone, province, first_onset, arrival_origin, arrival_days,
                     dplyr::across(c(dist_greatcircle_km, dist_road_km, travel_time_h,
                                     share_shorttrip, share_relocation,
                                     neglog_share_shorttrip, neglog_share_relocation,
                                     population, log10_population),
                                   ~ signif(.x, 12)),
                     cases_total),
  file.path(KEY_DIR, "arrival_predictors.csv"))

# Per-predictor fit of PREDICTOR on ARRIVAL DAYS. The orientation matters and is the
# panel's: arrival time on x, predictor on y, so `slope` is in predictor-units per day
# and (intercept, slope) draw exactly the line the panel shows — the figure plots this
# line from these two numbers rather than refitting. R^2 is symmetric in x and y, rho
# is rank-based, so neither depends on the orientation; the slope does.
arrival_predictor_fits <- local({
  preds <- c(neglog_share_shorttrip  = "-log(share of the epicentre region's short-trip outflow)",
             neglog_share_relocation = "-log(share of the epicentre region's March relocation outflow)",
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
  # arrival dates are imputed, so each replicate re-dates the SAME zones and refits.
  # A replicate can invade a zone the primary draw does not (its first usable onset
  # can land inside the window), so the replicate n is allowed to differ and is not
  # asserted equal to the primary n.
  rep_arr <- if (!is.null(imp) && nrow(imp) && "arrivals" %in% names(imp)) imp$arrivals else list()
  .rep_x <- lapply(rep_arr, function(a) {
    if (is.null(a) || !nrow(a)) return(NULL)
    as.numeric(a$first_onset[match(arrival_predictors$health_zone, a$health_zone)] -
                 min(a$first_onset))
  })
  .rep_x <- Filter(Negate(is.null), .rep_x)

  x_all <- arrival_predictors$arrival_days
  dplyr::bind_rows(lapply(names(preds), function(p) {
    y  <- arrival_predictors[[p]]
    f  <- .fit_one(x_all, y)
    k  <- is.finite(x_all) & is.finite(y)
    rr <- if (length(.rep_x)) vapply(.rep_x, function(xr) {
            g <- .fit_one(xr, y); c(g$r2, g$spearman) }, numeric(2)) else matrix(numeric(0), nrow = 2)
    .qq <- function(v, q) if (length(v) && any(is.finite(v))) unname(stats::quantile(v[is.finite(v)], q, names = FALSE)) else NA_real_
    tibble::tibble(
      predictor = p, label = unname(preds[p]), n = f$n,
      r2 = f$r2, intercept = f$intercept, slope = f$slope,
      pearson = f$pearson, spearman = f$spearman,
      n_imputation_replicates = ncol(rr),
      r2_lo       = .qq(rr[1, ], 0.05), r2_hi       = .qq(rr[1, ], 0.95),
      spearman_lo = .qq(rr[2, ], 0.05), spearman_hi = .qq(rr[2, ], 0.95),
      # The facets are NOT on common support, and the dropping is selection ON the
      # predictor, which is correlated with the outcome: a zone is dropped precisely
      # when no flow out of the epicentre was released for it. These two columns make
      # that selection visible, so a reader comparing R^2 across predictors can see
      # whether the denominators describe the same outbreak.
      mean_arrival_days_used    = if (any(k)) mean(x_all[k]) else NA_real_,
      mean_arrival_days_dropped = if (any(!k)) mean(x_all[!k], na.rm = TRUE) else NA_real_)
  }))
})
readr::write_csv(
  arrival_predictor_fits %>%
    dplyr::mutate(dplyr::across(c(r2, intercept, slope, pearson, spearman,
                                  r2_lo, r2_hi, spearman_lo, spearman_hi,
                                  mean_arrival_days_used, mean_arrival_days_dropped),
                                ~ signif(.x, 12))),
  file.path(KEY_DIR, "arrival_predictor_fits.csv"))
for (i in seq_len(nrow(arrival_predictor_fits))) with(arrival_predictor_fits[i, ],
  message(sprintf("[kinematics]   %-24s n=%2d  R2=%.3f  rho=%+.3f", predictor, n, r2, spearman)))

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
                   cog_dist_from_epi_km = round(cog_dist_from_epi_km, 2),
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

# Mean new invasions per week over the fit window, with a 90% bootstrap interval over WEEKS
# (the independent unit here: each week contributes one count). Published for Panel C of
# Figure_spread_kinematics_compact, which must not compute its own reference level.
.nz_rows <- K$weekly$week >= 1L & !K$weekly$truncated & is.finite(K$weekly$new_zones)
.nz_vals <- as.numeric(K$weekly$new_zones[.nz_rows])
.nz_fit_max <- if (any(.nz_rows)) max(K$weekly$week[.nz_rows]) else NA_integer_
.nz_mean <- if (length(.nz_vals)) mean(.nz_vals) else NA_real_
.nz_ci <- if (length(.nz_vals) >= 3L) {
  .old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(RANDOM_SEED_BASE)
  .bs <- vapply(seq_len(2000L),
                function(i) mean(.nz_vals[sample.int(length(.nz_vals), length(.nz_vals), TRUE)]),
                numeric(1))
  if (!is.null(.old_seed)) assign(".Random.seed", .old_seed, envir = globalenv())
  unname(stats::quantile(.bs, c(0.05, 0.95), names = FALSE))
} else c(NA_real_, NA_real_)

summary_out <- tibble::tribble(
  ~statistic, ~value, ~lo, ~hi, ~unit, ~basis,
  "Radial expansion rate (road)",         S$v_radial_road, ci_road[1], ci_road[2], "km/week",
    "OLS of road distance from the index zone on invasion week; 95% zone bootstrap",
  "Radial expansion rate (great-circle)", S$v_radial_gc, ci_gc[1], ci_gc[2], "km/week",
    "OLS of great-circle distance on invasion week; 95% zone bootstrap",
  "Radial expansion rate (road, Theil-Sen)", S$v_radial_road_ts, NA_real_, NA_real_, "km/week",
    "median of pairwise slopes — outlier-resistant check on the OLS road rate above",
  "Radial expansion rate (great-circle, Theil-Sen)", S$v_radial_gc_ts, NA_real_, NA_real_, "km/week",
    "median of pairwise slopes — outlier-resistant check on the OLS great-circle rate",
  "Radial expansion rate (travel time)", S$v_radial_time, NA_real_, NA_real_, "hours/week",
    sprintf("OLS of OSRM travel time on invasion week; the best-fitting geography of the three (R2 = %.2f vs %.2f road, %.2f great-circle)",
            S$v_radial_time_r2, S$v_radial_road_r2, S$v_radial_gc_r2),
  "Median-reach growth rate (road)",      S$v_reach_med, NA_real_, NA_real_, "km/week",
    "OLS of the MEDIAN distance of invaded zones on week (bulk expansion of the invaded set)",
  "Case-weighted mean distance trend",    S$v_case_mean, NA_real_, NA_real_, "km/week",
    "OLS of the case-weighted mean distance from the index zone on week",
  "Envelope: week reaching 90% of final extent", S$week_reach_90pct, NA_real_, NA_real_, "week",
    "first week in which R_max(t) attains 0.9x its fit-window final value",
  "Envelope expansion, to 90% of extent", S$v_front_early, NA_real_, NA_real_, "km/week",
    "R_max at that week divided by that week: how fast the extent was established",
  "Envelope expansion, after 90% of extent", S$v_front_late, NA_real_, NA_real_, "km/week",
    "OLS of R_max on week thereafter; ~0 means the extent has SATURATED (no wavefront). It does NOT mean long-range seeding stopped — see the two rows below, which are the statistics that show it did not",
  "New invasions after saturation, median furthest", S$newmax_med_after, NA_real_, NA_real_, "km",
    sprintf("median, over the %d post-saturation fit-window weeks with a new invasion, of the FURTHEST zone invaded that week — the recurring long-range seeding the running-maximum envelope hides",
            S$n_weeks_after),
  "New invasions after saturation, furthest", S$newmax_max_after, NA_real_, NA_real_, "km",
    "single furthest zone invaded in any post-saturation fit-window week",
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
    sprintf("slope(SD^2 ~ t)/4 — the isotropic-diffusion reading of the dispersion growth. The SD^2 ~ t fit has R2 = %.2f: a low R2 means the case cloud is not spreading out at all, so this is a null result, NOT a fitted diffusion constant to be quoted on its own.",
            S$diffusion_r2),
  "Invaded-zone accrual rate",            S$accrual_rate_wk, NA_real_, NA_real_, "log-zones/week",
    "OLS of log(cumulative invaded zones) on week",
  # PUBLISHED because a figure was computing it. make_spread_kinematics_compact.R drew Panel
  # C's reference line from its own mean of new_zones over the same weeks, and its comment
  # asserted the interval "is in spread_kinematics_summary.csv" — there was no such row. The
  # window is weeks 1..fit_max: week 0 is the seeding condition (those zones already had
  # confirmed cases when the series opens), and the truncated tail is incompletely observed.
  "Mean new invasions per week",          .nz_mean, .nz_ci[1], .nz_ci[2], "zones/week",
    sprintf("mean of new_zones over weeks 1-%d (week 0 = seeding, truncated tail excluded); 90%% week-bootstrap",
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

# Imputation-uncertainty columns (5th-95th percentile across independent onset draws).
imp_map <- c("Radial expansion rate (road)" = "v_radial_road",
             "Radial expansion rate (great-circle)" = "v_radial_gc",
             "Median-reach growth rate (road)" = "v_reach_med",
             "Case centre-of-gravity speed" = "v_cog_med",
             "Case centre-of-gravity net speed" = "v_cog_net",
             "Dispersion growth rate" = "v_disp_km_wk",
             "Invaded-zone doubling time" = "accrual_doubling_wk",
             "Reach at the end of the fit window" = "reach_max_km",
             "Zones invaded (fit window)" = "n_zones_fit")
iq <- t(vapply(summary_out$statistic, function(s)
  if (s %in% names(imp_map)) imp_q(unname(imp_map[s])) else c(NA_real_, NA_real_),
  numeric(2)))
summary_out <- summary_out %>%
  dplyr::mutate(value = round(value, 4), lo = round(lo, 4), hi = round(hi, 4),
                imp_lo = round(iq[, 1], 4), imp_hi = round(iq[, 2], 4)) %>%
  dplyr::relocate(imp_lo, imp_hi, .after = hi)
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
                     note = sprintf("v_radial=%.1f v_envelope_ols=%.1f v_cog=%.1f v_disp=%.1f n_zones=%d",
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
message(sprintf("  radial expansion  (road)    %.1f km/wk  [95%% zone bootstrap %.1f, %.1f]",
                S$v_radial_road, ci_road[1], ci_road[2]))
message(sprintf("  radial expansion  (gc)      %.1f km/wk  [95%% zone bootstrap %.1f, %.1f]",
                S$v_radial_gc, ci_gc[1], ci_gc[2]))
message(sprintf("  robustness of that rate     Theil-Sen %.1f (road) / %.1f (gc) km/wk; %.1f km/wk dropping %d sitrep-only zone(s)",
                S$v_radial_road_ts, S$v_radial_gc_ts, nosit_road$slope, nosit_n))
message(sprintf("  median-reach growth (road)  %.1f km/wk", S$v_reach_med))
message(sprintf("  envelope                    reached 90%% of its final extent by week %.0f (%.0f km/wk to then, %+.1f km/wk after) -> %s",
                S$week_reach_90pct, S$v_front_early, S$v_front_late,
                if (isTRUE(abs(S$v_front_late) < 0.1 * abs(S$v_front_early)))
                  "SATURATED: no travelling wavefront" else "still expanding"))
message(sprintf("  ...but seeding CONTINUED     in the %d post-saturation week(s) with a new invasion, the furthest new zone was a median %.0f km out (max %.0f km) — invisible in the running-max envelope",
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
  message(sprintf("     v_radial %.1f [%.1f, %.1f]   v_cog %.1f [%.1f, %.1f]   max reach %.0f [%.0f, %.0f] km",
                  stats::median(imp$v_radial_road), imp_q("v_radial_road")[1], imp_q("v_radial_road")[2],
                  stats::median(imp$v_cog_med), imp_q("v_cog_med")[1], imp_q("v_cog_med")[2],
                  stats::median(imp$reach_max_km), imp_q("reach_max_km")[1], imp_q("reach_max_km")[2]))
  if (dplyr::n_distinct(imp$n_weeks) > 1L)
    message(sprintf("     note: onset draws produced %d different week-grid lengths (%s); the fit window is calendar-anchored, so the runs stay comparable",
                    dplyr::n_distinct(imp$n_weeks),
                    paste(sort(unique(imp$n_weeks)), collapse = "/")))
}
message(strrep("-", 74), "\n")

# =============================================================================
# SECTION 6 — figure
# =============================================================================
W  <- K$weekly
WF <- W %>% dplyr::filter(!truncated)
Z  <- K$zones

# Province palette: identity hues for the outbreak provinces, disjoint spares else.
prov_palette <- function(provs) {
  provs <- sort(unique(as.character(provs[!is.na(provs)])))
  provs <- c(intersect(names(PROV_COL), provs), setdiff(provs, names(PROV_COL)))
  spare <- setdiff(SPARE_COL, unname(PROV_COL)); si <- 1L
  cols <- setNames(character(length(provs)), provs)
  for (p in provs) {
    if (p %in% names(PROV_COL)) cols[p] <- unname(PROV_COL[p])
    else { cols[p] <- spare[((si - 1L) %% length(spare)) + 1L]; si <- si + 1L }
  }
  stopifnot(!anyNA(cols), all(nzchar(cols)))
  cols
}
PCOL <- prov_palette(Z$province)
Z$province <- factor(Z$province, levels = names(PCOL))

# ggplot transforms the rect bounds with the y scale, so the default +/-Inf bounds
# become NaN/Inf under a log scale (and warn). Log-scaled panels therefore pass
# finite bounds that comfortably bracket their data instead.
trunc_band <- function(ymin = -Inf, ymax = Inf) if (K$n_trunc > 0)
  # xmax = Inf, NOT week_max + 0.5: a finite bound is data and would push the x scale
  # out to it, leaving a white strip between the band and the panel edge. Infinite
  # bounds are excluded from the scale range and drawn flush against the panel.
  annotate("rect", xmin = K$fit_max + 0.5, xmax = Inf,
           ymin = ymin, ymax = ymax, fill = "grey55", alpha = 0.10) else NULL
WK_SCALE <- scale_x_continuous(breaks = scales::breaks_width(2),
                               expand = expansion(mult = c(0.02, 0.03)))
SHAPE_TRUNC <- scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none")
#' "value [lo, hi] unit" — but with the bracket dropped when the interval is
#' unavailable (e.g. SPREAD_IMPUTE_REPS=0), so no panel can ever print "[NA, NA]".
fmt_ci <- function(v, ci, unit, digits = 1) {
  f <- paste0("%+.", digits, "f")
  if (any(!is.finite(ci))) sprintf(paste0(f, " %s"), v, unit)
  else sprintf(paste0(f, " [%.", digits, "f, %.", digits, "f] %s"), v, ci[1], ci[2], unit)
}

#' In-panel annotation anchored to a corner known to be empty in that panel.
note <- function(txt, where = c("topleft", "bottomright")) {
  where <- match.arg(where)
  if (where == "topleft")
    annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.7,
             size = 2.5, colour = MUTED, label = txt)
  else
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.9,
             size = 2.5, colour = MUTED, label = txt)
}

# ---- Panel A — map ----------------------------------------------------------
# The invaded set spans ~800 km but the case centre of gravity only ever moves a few
# tens of km, so at national scale its track collapses to a blob. INSET_HALF_DEG
# controls a zoom on the origin focus that is drawn into the empty corner of the map
# as a grob, keeping Panel A a single ggplot (so patchwork still tags it once).
INSET_HALF_DEG <- 0.85

build_core_inset <- function(shp) {
  ex <- c(K$epi_ll[1] - INSET_HALF_DEG, K$epi_ll[1] + INSET_HALF_DEG,
          K$epi_ll[2] - INSET_HALF_DEG, K$epi_ll[2] + INSET_HALF_DEG)
  cg <- W %>% dplyr::filter(is.finite(cog_lon),
                            cog_lon > ex[1], cog_lon < ex[2],
                            cog_lat > ex[3], cog_lat < ex[4])
  arw <- grid::arrow(length = unit(3, "pt"), type = "closed")
  ggplot() +
    geom_sf(data = shp, aes(fill = invasion_week), colour = "white", linewidth = 0.12) +
    geom_path(data = cg, aes(cog_lon, cog_lat), colour = CASE_COL, linewidth = 0.6,
              arrow = arw, lineend = "round") +
    geom_point(data = cg, aes(cog_lon, cog_lat), shape = 21, size = 1.7,
               fill = CASE_COL, colour = "white", stroke = 0.4) +
    ggrepel::geom_text_repel(data = cg[c(1, nrow(cg)), ],
                             aes(cog_lon, cog_lat, label = paste0("wk ", week)),
                             size = 2.1, colour = INK, seed = 3L, min.segment.length = 0,
                             segment.colour = FAINT, segment.size = 0.22,
                             bg.colour = "white", bg.r = 0.14, box.padding = 0.4) +
    geom_point(data = data.frame(lon = K$epi_ll[1], lat = K$epi_ll[2]),
               aes(lon, lat), shape = 21, size = 2.6, fill = EPI_COL,
               colour = "white", stroke = 0.5) +
    scale_fill_viridis_c(option = "mako", begin = 0.08, end = 0.92,
                         na.value = "grey92", limits = c(0, max(Z$inv_week)),
                         guide = "none") +
    coord_sf(xlim = ex[1:2], ylim = ex[3:4], expand = FALSE) +
    theme_void() +
    theme(panel.background = element_rect(fill = "white", colour = "grey45",
                                          linewidth = 0.4),
          plot.margin = margin(0, 0, 0, 0))
}

build_map <- function() {
  shp <- geo$shp
  shp$invasion_week <- Z$inv_week[match(shp$Nom, Z$health_zone)]
  # Crop to the invaded footprint plus a margin, so the epicentre rings and the
  # centre-of-gravity tracks all sit inside the frame.
  bb  <- sf::st_bbox(shp[!is.na(shp$invasion_week), ])
  pad <- 0.3
  xlim <- c(bb[["xmin"]] - pad, bb[["xmax"]] + pad)
  ylim <- c(bb[["ymin"]] - pad, bb[["ymax"]] + pad)
  box <- sf::st_as_sfc(sf::st_bbox(c(xmin = xlim[1], ymin = ylim[1],
                                     xmax = xlim[2], ymax = ylim[2]),
                                   crs = sf::st_crs(shp)))
  keep <- shp[as.vector(sf::st_intersects(shp, box, sparse = FALSE)), ]
  prov <- keep %>% dplyr::group_by(PROVINCE) %>% dplyr::summarise(.groups = "drop")

  rings <- dplyr::bind_rows(lapply(RINGS_KM, function(km) {
    p <- geosphere::destPoint(K$epi_ll, b = seq(0, 360, by = 2), d = km * 1000)
    tibble::tibble(lon = p[, 1], lat = p[, 2], km = km)
  }))
  # Ring labels due WEST of the epicentre: the southward radials run through the
  # labelled southern zones (Goma, Miti Murhesa) and the south-west corner holds the
  # inset, so west is the only quarter that is empty at every radius. Clipped to the
  # visible extent and set on a translucent white plate so they stay readable over any fill.
  rlab <- dplyr::bind_rows(lapply(RINGS_KM, function(km) {
    p <- geosphere::destPoint(K$epi_ll, b = 260, d = km * 1000)
    tibble::tibble(lon = p[1, 1], lat = p[1, 2], km = km)
  })) %>%
    dplyr::filter(lat > ylim[1], lat < ylim[2], lon > xlim[1], lon < xlim[2]) %>%
    # The qualifier goes on the SMALLEST ring: its label is the one guaranteed to be
    # inside the frame, so the straight-line/road distinction is never lost to clipping.
    dplyr::mutate(lab = ifelse(km == min(km), paste0(km, " km straight-line"),
                               paste0(km, " km")))

  cg <- W %>% dplyr::filter(is.finite(cog_lon))
  ig <- W %>% dplyr::filter(is.finite(inv_cog_lon))
  # Label the index zone plus the FURTHEST invaded zone in each province: one label
  # per direction of travel, rather than four labels on the same distant cluster.
  labz <- dplyr::bind_rows(
    Z %>% dplyr::filter(health_zone == K$epi_zone),
    Z %>% dplyr::filter(health_zone != K$epi_zone, is.finite(d_road)) %>%
      dplyr::group_by(province) %>%
      dplyr::slice_max(d_road, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()) %>%
    dplyr::distinct(health_zone, .keep_all = TRUE)
  arw <- grid::arrow(length = unit(3.4, "pt"), type = "closed")
  # Legend for the two tracks, drawn as a manual colour scale on invisible points. The
  # shared "Centre of gravity" wording lives in the legend TITLE, not in each label:
  # spelled out per key it overflows the map cell and crowds Panel D's axis title.
  trk <- tibble::tibble(lon = NA_real_, lat = NA_real_,
                        track = factor(c("Weekly cases", "New invasions"),
                                       levels = c("Weekly cases", "New invasions")))

  # Zoom on the origin focus, dropped into the empty south-west corner of the frame,
  # with a leader box marking the area it magnifies.
  inset <- ggplot2::ggplotGrob(build_core_inset(keep))
  iw <- 0.30 * diff(xlim); ih <- iw                       # square, in degrees
  ibox <- c(xlim[1] + 0.02 * diff(xlim), ylim[1] + 0.02 * diff(ylim))
  zoombox <- data.frame(
    xmin = K$epi_ll[1] - INSET_HALF_DEG, xmax = K$epi_ll[1] + INSET_HALF_DEG,
    ymin = K$epi_ll[2] - INSET_HALF_DEG, ymax = K$epi_ll[2] + INSET_HALF_DEG)

  ggplot() +
    geom_sf(data = keep, aes(fill = invasion_week), colour = "white", linewidth = 0.09) +
    geom_sf(data = prov, fill = NA, colour = "grey45", linewidth = 0.28) +
    geom_path(data = rings, aes(lon, lat, group = km),
              colour = "grey40", linewidth = 0.24, linetype = "22") +
    geom_label(data = rlab, aes(lon, lat, label = lab),
               size = 2.05, colour = MUTED, fill = "white", alpha = 0.7,
               linewidth = 0, label.padding = unit(1.2, "pt")) +
    # The new-invasion centre jumps across the frame week to week — that erratic
    # track IS the finding (no coherent wavefront), so it is drawn, but muted so it
    # does not overpower the case centre.
    geom_path(data = ig, aes(inv_cog_lon, inv_cog_lat), colour = INV_COL, linewidth = 0.32,
              linetype = "31", alpha = 0.65, lineend = "round") +
    geom_point(data = ig, aes(inv_cog_lon, inv_cog_lat), shape = 23, size = 1.4,
               fill = "white", colour = INV_COL, stroke = 0.45) +
    geom_path(data = cg, aes(cog_lon, cog_lat), colour = CASE_COL, linewidth = 0.55,
              arrow = arw, lineend = "round") +
    geom_point(data = cg, aes(cog_lon, cog_lat), shape = 21, size = 1.5,
               fill = CASE_COL, colour = "white", stroke = 0.35) +
    geom_point(data = trk, aes(lon, lat, colour = track), size = 1.4, na.rm = TRUE) +
    ggrepel::geom_text_repel(data = labz, aes(lon, lat, label = health_zone),
                             size = 2.35, colour = INK, seed = 1L, min.segment.length = 0,
                             segment.colour = FAINT, segment.size = 0.25,
                             bg.colour = "white", bg.r = 0.13,
                             box.padding = 0.34, point.padding = 0.2, max.overlaps = 30) +
    geom_point(data = data.frame(lon = K$epi_ll[1], lat = K$epi_ll[2]),
               aes(lon, lat), shape = 21, size = 2.8, fill = EPI_COL,
               colour = "white", stroke = 0.55) +
    geom_rect(data = zoombox, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, colour = "grey45", linewidth = 0.35) +
    annotation_custom(inset, xmin = ibox[1], xmax = ibox[1] + iw,
                      ymin = ibox[2], ymax = ibox[2] + ih) +
    scale_fill_viridis_c(option = "mako", begin = 0.08, end = 0.92, na.value = "grey92",
                         name = "Week of\nfirst case", breaks = scales::breaks_width(3),
                         limits = c(0, max(Z$inv_week)),
                         guide = guide_colourbar(order = 1)) +
    scale_colour_manual(values = setNames(c(CASE_COL, INV_COL), levels(trk$track)),
                        name = "Centre of\ngravity", drop = FALSE,
                        guide = guide_legend(order = 2, ncol = 1,
                                             override.aes = list(size = 2, shape = 16))) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_map(9) +
    theme(legend.box = "vertical", legend.spacing.y = unit(4, "pt"))
}

# ---- Panel B — reach from the epicentre ------------------------------------
build_reach <- function() {
  # Labels kept SHORT: four keys in a two-column legend inside a 3.5-inch panel, and a
  # long key silently overflows the cell and is clipped. "so far" vs "that week" is what
  # carries the cumulative/per-week distinction that this panel exists to make.
  lv <- c("Furthest zone so far", "Furthest new zone that week",
          "Median zone so far", "Weekly cases, mean")
  # Blues = the cumulative invaded SET (dark = its outer envelope, light = its middle);
  # green = the zones newly invaded THAT WEEK, as in Panels A and D; orange = the CASES.
  RCOL <- setNames(c(EXT_COL, INV_COL, EXT_COL2, CASE_COL), lv)
  d <- W %>%
    dplyr::transmute(week, truncated,
                     `Furthest zone so far`        = reach_max_km,
                     `Furthest new zone that week` = new_max_km,
                     `Median zone so far`          = reach_med_km,
                     `Weekly cases, mean`          = case_mean_km) %>%
    tidyr::pivot_longer(-c(week, truncated), names_to = "series", values_to = "km") %>%
    dplyr::filter(is.finite(km)) %>%
    dplyr::mutate(series = factor(series, levels = lv))
  # The fitted line is on the MEDIAN series. R_max is a running maximum that
  # saturates once the long-range seeding events have happened, so fitting a line
  # to it would advertise a travelling wavefront that the data do not show; the
  # saturation is annotated instead.
  fitline <- tibble::tibble(week = WF$week,
                            km = K$reach_med_intercept + S$v_reach_med * WF$week)
  ggplot(d, aes(week, km, colour = series)) +
    trunc_band() +
    # Two different dashed marks would otherwise be indistinguishable, so they are
    # separated on BOTH linetype and colour: a DOTTED GREY vertical rule is a time
    # marker (it belongs to no series), while the fit line takes the colour of the
    # series it is fitted to, so it can only be read against that series.
    geom_vline(xintercept = S$week_reach_90pct, colour = FAINT,
               linewidth = 0.35, linetype = "12") +
    geom_line(data = fitline, aes(week, km), inherit.aes = FALSE,
              colour = EXT_COL2, alpha = 0.75, linewidth = 0.4, linetype = "22") +
    geom_line(aes(group = series), linewidth = 0.5) +
    geom_point(aes(shape = truncated), size = 1.3, stroke = 0.45, fill = "white") +
    SHAPE_TRUNC +
    # Two rows: on one row the third key overflows the cell and crowds Panel C's legend.
    scale_colour_manual(values = RCOL, name = NULL,
                        guide = guide_legend(nrow = 2, byrow = TRUE)) +
    WK_SCALE + expand_limits(y = 0) +
    # States BOTH halves of the finding, because the flat envelope on its own reads as
    # "the spread stopped", which the green series shows it did not. The per-week series
    # sweeps the full height of the panel every other week, so there is NO empty region
    # to sit in: it goes on a translucent white plate (as the map's ring labels do) and
    # stays legible wherever a series happens to pass under it.
    annotate("label", x = 0.15, y = 0.96 * max(d$km, na.rm = TRUE),
             hjust = 0, vjust = 1, size = 2.45, colour = MUTED,
             # linewidth, NOT label.size: the latter is deprecated in ggplot2 4.x and is
             # ignored, leaving a border box on the plate that no other panel has.
             fill = scales::alpha("white", 0.78), linewidth = 0,
             label.padding = unit(1.5, "pt"), lineheight = 1.05,
             label = sprintf(paste0("extent fixed by week %.0f and flat since,\n",
                                    "but new invasions kept landing a\n",
                                    "median %.0f km out (max %.0f km).\n",
                                    "Median zone so far %s"),
                             S$week_reach_90pct, S$newmax_med_after, S$newmax_max_after,
                             fmt_ci(S$v_reach_med, imp_q("v_reach_med"), "km/wk"))) +
    labs(x = "Epidemic week", y = "Road distance from index zone (km)") +
    theme_pub()
}

# ---- Panel C — arrival week vs distance -------------------------------------
build_arrival <- function() {
  zf <- Z %>% dplyr::filter(is.finite(d_road))
  fitd <- tibble::tibble(inv_week = seq(0, K$fit_max, by = 0.1)) %>%
    dplyr::mutate(d_road = K$radial_intercept + S$v_radial_road * inv_week)
  ggplot(zf, aes(inv_week, d_road)) +
    trunc_band() +
    geom_line(data = fitd, colour = "grey35", linewidth = 0.5) +
    # Invasion weeks are integers, so identically-timed zones would overplot; a small
    # deterministic horizontal jitter separates them without moving anything a full week.
    geom_point(aes(fill = province, alpha = in_fit_window),
               shape = 21, size = 1.9, colour = "white", stroke = 0.3,
               position = position_jitter(width = 0.16, height = 0, seed = 7L)) +
    scale_fill_manual(values = PCOL, name = NULL, drop = FALSE,
                      guide = guide_legend(nrow = 2, byrow = TRUE,
                                           override.aes = list(size = 2.2, alpha = 1))) +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.95), guide = "none") +
    WK_SCALE + expand_limits(y = 0) +
    # Both R² values, because which GEOGRAPHY orders the arrivals is itself a result:
    # travel time beats road km, which beats straight-line (full table in the
    # diagnostics CSV). The bracket is the zone bootstrap, not the onset draws.
    note(sprintf("%s, 95%% CI\nR² %.2f road · %.2f travel time",
                 fmt_ci(S$v_radial_road, ci_road, "km/week", digits = 0),
                 S$v_radial_road_r2, S$v_radial_time_r2), "bottomright") +
    labs(x = "Week of first confirmed case", y = "Road distance from index zone (km)") +
    theme_pub()
}

# ---- Panel D — centre-of-gravity step speed ---------------------------------
build_cogspeed <- function() {
  lv <- c("Weekly cases", "Newly invaded zones")
  CCOL <- setNames(c(CASE_COL, INV_COL), lv)
  d <- W %>%
    dplyr::transmute(week, truncated,
                     `Weekly cases` = cog_speed_km_wk,
                     `Newly invaded zones` = inv_cog_speed_km_wk) %>%
    tidyr::pivot_longer(-c(week, truncated), names_to = "series", values_to = "kmwk") %>%
    dplyr::filter(is.finite(kmwk)) %>%
    dplyr::mutate(series = factor(series, levels = lv))
  # Log y: the invasion centre moves 1-2 orders of magnitude faster than the case
  # centre, so a linear axis flattens the case series onto the baseline.
  # The band must run the FULL panel height, as it does in every other panel. Two traps:
  #   * +/-Inf bounds (what the other panels use) are illegal here — log10(-Inf) is NaN.
  #   * finite bounds do NOT get clipped: annotate("rect") is data, so an overshooting
  #     rect EXPANDS the y scale and squashes the series into a sliver.
  # So the window is pinned with coord_cartesian (which really does clip) and the rect
  # is then free to overshoot it.
  rng <- range(d$kmwk)
  ylo <- rng[1] / 1.15; yhi <- rng[2] * 1.15
  ggplot(d, aes(week, kmwk, colour = series)) +
    trunc_band(ymin = ylo / 10, ymax = yhi * 10) +
    geom_hline(yintercept = S$v_cog_med, colour = CASE_COL, linewidth = 0.35, linetype = "22") +
    geom_hline(yintercept = S$v_invcog_med, colour = INV_COL, linewidth = 0.35, linetype = "22") +
    geom_line(aes(group = series), linewidth = 0.5) +
    geom_point(aes(shape = truncated), size = 1.3, stroke = 0.45, fill = "white") +
    SHAPE_TRUNC +
    scale_colour_manual(values = CCOL, name = NULL) +
    WK_SCALE +
    scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300), labels = c("1","3","10","30","100","300")) +
    annotation_logticks(sides = "l", colour = FAINT, linewidth = 0.2,
                        short = unit(1, "pt"), mid = unit(1.6, "pt"), long = unit(2.4, "pt")) +
    coord_cartesian(ylim = c(ylo, yhi)) +
    note(sprintf("medians: cases %.0f, invasions %.0f km/week",
                 S$v_cog_med, S$v_invcog_med)) +
    labs(x = "Epidemic week", y = "Centre-of-gravity speed (km/week, log scale)") +
    theme_pub()
}

# ---- Panel E — dispersion ---------------------------------------------------
build_disp <- function() {
  lv <- c("Weekly cases (case-weighted)", "Invaded zones (cumulative)")
  # NOT grey for the invaded set: grey is this figure's annotation ink (fit lines,
  # truncation band, time markers), so a grey data series reads as an annotation.
  DCOL <- setNames(c(CASE_COL, EXT_COL), lv)
  d <- W %>%
    dplyr::transmute(week, truncated,
                     `Weekly cases (case-weighted)` = cog_sd_km,
                     `Invaded zones (cumulative)`   = invaded_sd_km) %>%
    tidyr::pivot_longer(-c(week, truncated), names_to = "series", values_to = "km") %>%
    dplyr::filter(is.finite(km)) %>%
    dplyr::mutate(series = factor(series, levels = lv))
  fitline <- tibble::tibble(week = WF$week,
                            km = K$disp_intercept + S$v_disp_km_wk * WF$week)
  ggplot(d, aes(week, km, colour = series)) +
    trunc_band() +
    # Fitted on the CASE-weighted series only, so it carries that series' colour.
    geom_line(data = fitline, aes(week, km), inherit.aes = FALSE,
              colour = CASE_COL, alpha = 0.75, linewidth = 0.4, linetype = "22") +
    geom_line(aes(group = series), linewidth = 0.5) +
    geom_point(aes(shape = truncated), size = 1.3, stroke = 0.45, fill = "white") +
    SHAPE_TRUNC +
    scale_colour_manual(values = DCOL, name = NULL,
                        guide = guide_legend(nrow = 2, byrow = TRUE)) +
    WK_SCALE + expand_limits(y = 0) +
    # The interval spans zero: the case cloud is not measurably spreading out, and
    # the fitted line must not be read as growth. Same onset-draw basis as Panel B.
    note(paste0(fmt_ci(S$v_disp_km_wk, imp_q("v_disp_km_wk"), "km/week"), ", onset draws")) +
    labs(x = "Epidemic week", y = "Standard distance (km)") +
    theme_pub()
}

# Flat wrap_plots (NOT nested patchworks) so tag_levels tags every panel A-E.
# The map cell is sized ~square to match the study area's own aspect ratio; a taller
# cell would letterbox the map and waste half the panel.
design <- "
AAAABBCC
AAAABBCC
AAAADDEE
AAAADDEE
"
fig <- patchwork::wrap_plots(A = build_map(), B = build_reach(), C = build_arrival(),
                             D = build_cogspeed(), E = build_disp(), design = design) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold", colour = INK))

save_dual(fig, "Figure_spread_kinematics", w = 14.0, h = 7.0)
message("[kinematics] done — CSVs in ", KEY_DIR)

}   # end of the script-run guard opened at SECTION 4
