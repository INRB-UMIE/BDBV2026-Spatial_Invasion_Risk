# =============================================================================
# 03_mobility_matrices.R — Mobility Matrix Construction
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: Build the mobility kernels that parameterise between-zone transmission
#   in the spatiotemporal metapopulation model: M1, M2a, M2b, M3, M4, M4b, M5, M6a,
#   M6b, M7, M8, M10 always; M9/M15 behind their flags; the Flowminder-cohort
#   composites M13/M14/M16 and the consensus M17; the source-cell fill variants
#   (*-fill); and the OSRM road-distance twins (*-dist). M11 (inward FOI) is built
#   on demand in run_all.R.
#
#   Zone names in the short-trip and cohort SOURCES are resolved here, against the
#   519-zone spine, rather than trusting the upstream processed matrices — see
#   load_short_trip_profile() for why.
#
# Matrices are 519×519 (national health zones from WorldPop), named by
# canonical zone name, and row-stochastic (each row sums to ≤ 1; zero rows
# are permitted for zones with no observed outflow).
#
# Data sources:
#   M1 / M2a / M2b : Flowminder short-trip cohort (epicentre pooled)
#   M3             : Flowminder national RELOCATION matrix (est_flows_2026_03: monthly home-location
#                    changes Feb->Mar 2026, NOT trips)
#   M4             : Negative-binomial gravity model (calibrated on the M3 relocations)
#   M5             : Radiation model (Simini et al. 2012 Nature 484:96-100)
#   M6a / M6b      : OSRM travel-time decay (exponential / power-law)
#   M7             : IDP-augmented hybrid (M3 relocations + IDP flows)
#   M8             : Composite — short trips for epicentre, gravity elsewhere (RECOMMENDED)
#   M13/M14/M16    : Flowminder cohort presence rows over gravity / radiation / relocation OD
#                    (M16's base is M3 with its coverage gaps filled: cover_relocation_od())
#   M17            : Consensus of M3 (coverage-filled) + M4 + M5 with empirical source rows overlaid
#   *-fill         : as above, but cells the source could not observe come from the base
#                    (the default form; MOBILITY_SOURCE_FILL = "unmeasured")
#   *-split        : cohort composites whose pooled cohort profile is disaggregated per origin
#                    (split_cohort_rows()), then filled
#
# Usage:
#   source("spatiotemporal/03_mobility_matrices.R")
#   # or call build_all_mobility_matrices() directly for downstream scripts.
#
# Outputs:
#   outputs/mobility/mobility_M*.rds          — individual matrix files
#   outputs/mobility/mobility_summary.csv     — sparsity and top-destination stats
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# ONE canonical definition, identical in every module that declares it, so that source()
# ORDER CANNOT CHANGE SEMANTICS. This file's version used to be
#   function(a, b) if (!is.null(a) && !is.na(a)[1]) a else b
# and was declared UNCONDITIONALLY, so sourcing 06 installed it over every other module's and
# it governed the whole cascade run. It differs in two ways that silently corrupt results:
# it falls back whenever the FIRST ELEMENT of a vector is NA (blanking an entire otherwise-good
# vector, e.g. run_cascade.R's R_zone_conjugate column or 33b's pop_vec), and it ERRORS on a
# zero-length left-hand side ("missing value where TRUE/FALSE needed") instead of falling back.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

#' Load zone name aliases for harmonisation.
#'
#' @return Named character vector: names = observed_name, values = canonical_nom.
load_aliases <- function() {
  stopifnot(file.exists(ALIASES_PATH))
  df <- readr::read_csv(ALIASES_PATH, col_types = readr::cols(.default = "c"),
                        show_col_types = FALSE)
  stopifnot(all(c("observed_name", "canonical_nom") %in% colnames(df)))
  al <- setNames(df$canonical_nom, df$observed_name)
  message(sprintf("[aliases] Loaded %d zone name alias entries.", length(al)))
  al
}

#' Harmonise a vector of zone names against the canonical list.
#'
#' Applies the alias table (load_aliases()); a name already canonical, or with no alias
#' entry, is returned unchanged. Names that are still outside `canonical` afterwards are
#' returned AS-IS with a warning (they are dropped later by the name-indexed embedders).
#' NOTE: the alias table is province-blind, so it maps the bare duplicate noms "Lubunga"
#' and "Bili" onto one of their two zones. Sources that carry a province column must use
#' .resolve_source_labels() instead.
#'
#' @param names_vec  Character vector of (possibly non-canonical) zone names.
#' @param aliases    Named vector from load_aliases().
#' @param canonical  Character vector of canonical zone names (zones_all).
#' @return Character vector same length as names_vec with canonical names.
harmonise_names <- function(names_vec, aliases, canonical) {
  harmonised <- dplyr::if_else(names_vec %in% names(aliases),
                                aliases[names_vec],
                                names_vec)
  unresolved <- harmonised[!harmonised %in% canonical]
  if (length(unresolved) > 0L) {
    warning(sprintf(
      "[harmonise] %d zone name(s) could not be matched to canonical list: %s",
      length(unresolved),
      paste(unique(unresolved), collapse = ", ")
    ))
  }
  harmonised
}

#' Make a matrix row-stochastic.
#'
#' Divides each row by its sum.  Rows with sum == 0 are left as zero
#' (representing no known outflow from that zone).  The diagonal is
#' zeroed before normalisation to prevent self-loops.
#'
#' @param M  Numeric matrix (zones × zones).
#' @return   Row-stochastic matrix with the same dimensions and names as M.
make_row_stochastic <- function(M) {
  stopifnot(is.matrix(M), nrow(M) == ncol(M))

  # Zero diagonal first
  diag(M) <- 0

  rs  <- rowSums(M, na.rm = TRUE)
  # Avoid division by zero for zero-outflow rows
  rs[rs == 0] <- NA_real_
  M_norm <- M / rs
  M_norm[is.na(M_norm)] <- 0   # rows with no outflow remain zero

  # STRIP inherited attributes. R's arithmetic copies custom attributes from its
  # first operand, so a kernel built on top of another (M1 over the M3 fallback,
  # every composite over its base) silently carried that one's payload — M1 shipped
  # M3's full "raw" 519x519 flow matrix, and the composites shipped M4's fitted
  # "model_summary". Only dim/dimnames survive here; every builder attaches its own
  # metadata AFTER normalising.
  .dd <- dim(M_norm); .dn <- dimnames(M_norm)
  attributes(M_norm) <- NULL
  dim(M_norm) <- .dd; dimnames(M_norm) <- .dn

  M_norm
}

#' Assert that a mobility matrix satisfies required invariants.
#'
#' Stops with a descriptive error if any invariant is violated.
#'
#' @param W         Row-stochastic mobility matrix.
#' @param zones_all Expected zone name vector.
#' @param label     Human-readable label for error messages.
assert_mobility_matrix <- function(W, zones_all, label = "W") {
  n <- length(zones_all)

  if (!all(dim(W) == c(n, n)))
    stop(sprintf("[assert] %s: expected dim c(%d,%d), got c(%d,%d)",
                 label, n, n, nrow(W), ncol(W)))

  if (!identical(rownames(W), zones_all))
    stop(sprintf("[assert] %s: rownames do not match zones_all", label))

  if (!identical(colnames(W), zones_all))
    stop(sprintf("[assert] %s: colnames do not match zones_all", label))

  rs <- rowSums(W)
  if (!all(rs >= -1e-9 & rs <= 1 + 1e-9))
    stop(sprintf("[assert] %s: row sums out of [0,1]: min=%.4f, max=%.4f",
                 label, min(rs), max(rs)))

  if (!all(W >= -1e-9))
    stop(sprintf("[assert] %s: negative weights present (min=%.4f)", label, min(W)))

  if (!all(abs(diag(W)) < 1e-9))
    stop(sprintf("[assert] %s: non-zero diagonal (max abs diag=%.4g)",
                 label, max(abs(diag(W)))))

  invisible(TRUE)
}

#' Initialise a named zero matrix for all zones.
#'
#' @param zones_all  Character vector of canonical zone names.
#' @return Named zero matrix of dim c(n, n).
make_zero_matrix <- function(zones_all) {
  n <- length(zones_all)
  M <- matrix(0, nrow = n, ncol = n,
               dimnames = list(zones_all, zones_all))
  M
}

# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

#' Load WorldPop population data.
#'
#' @return Named numeric vector: names = canonical zone name, values = pop_count.
load_worldpop <- function() {
  f <- file.path(WORLDPOP_DIR, "worldpop__pop_count__static.csv")
  stopifnot(file.exists(f))
  df <- readr::read_csv(f, col_types = readr::cols(.default = "c"),
                        show_col_types = FALSE)
  # File columns: (index), nom, pop_count
  stopifnot(all(c("nom", "pop_count") %in% colnames(df)))
  pop <- setNames(as.numeric(df$pop_count), df$nom)
  message(sprintf("[worldpop] Loaded %d zones (total pop: %.0f)", length(pop), sum(pop, na.rm=TRUE)))
  pop
}

# ---------------------------------------------------------------------------
# OSRM coverage gaps: zones the road network cannot reach at all
# ---------------------------------------------------------------------------
# Two of the 519 zones -- Bokoro (Mai-Ndombe) and Idjwi (Sud-Kivu) -- have NO finite OSRM
# entry in either direction (518 of 518 off-diagonal cells NA each); every other zone has at
# most 2 NAs, which are its own cells to those two. That is not a data defect: Idjwi is an
# island in Lake Kivu and Bokoro is river-accessed, so there genuinely is no road route.
#
# The defect is the INFERENCE the kernels drew from it. build_M4 zeroes `na_dist_mask` and
# build_M5/build_M6 skip non-finite costs. In the ORIGIN direction that is right (an
# unroutable origin must not radiate uniformly -- see build_M6's note). In the DESTINATION
# direction it gives the zone an all-zero COLUMN: structurally zero import hazard, so the
# zone can never be invaded under that kernel, whatever happens around it. Measured on the
# 2026-09-22 build, Bokoro and Idjwi were the only zero-inflow zones of M4, M5, M6a and the
# whole M8/M10 family. The Flowminder cohort tables RECORD subscriber presence in both
# (people reach them by boat and by river), which is why the cohort families M13/M14/M16/M17
# never had the problem -- so the mobility data contradict the zero directly.
#
# The repair imputes the missing costs from the great-circle separation of the zone
# centroids, through the road network's OWN cost-vs-separation relation: a log-log fit on
# the routable pairs, so an imputed pair lands where a routable pair at the same separation
# would. Deliberately NO water-crossing penalty is applied -- any multiplier would be a free
# parameter with nothing behind it, and this keeps the imputation a statement about distance
# rather than a guess about ferries. Set OSRM_GAP_FILL = FALSE to restore the previous
# (zero-inflow) behaviour exactly.

.osrm_centroid_cache <- new.env(parent = emptyenv())

#' Zone centroids (lon/lat, WGS84) from the health-zone shapefile.
#'
#' The shapefile's `Nom` IS the canonical 519-zone spine -- verified identical as a SET to the
#' WorldPop `nom` column, with no duplicated names -- so no alias resolution is needed here.
#' Uses st_point_on_surface with s2 disabled (as 17_invasion_viz.R does): a few polygons have
#' invalid rings, and a true centroid can fall outside a concave zone.
#'
#' @param path Shapefile path (default SHAPEFILE_PATH).
#' @return data.frame(nom, lon, lat), one row per zone, or NULL when sf or the file is absent.
load_zone_centroids <- function(path = get0("SHAPEFILE_PATH", ifnotfound = NA_character_)) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    warning("[centroids] shapefile not found; cannot compute zone centroids.", call. = FALSE)
    return(NULL)
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    warning("[centroids] package 'sf' is not installed; cannot compute zone centroids.",
            call. = FALSE)
    return(NULL)
  }
  key <- paste0("centroids:", path)
  if (!is.null(.osrm_centroid_cache[[key]])) return(.osrm_centroid_cache[[key]])
  shp <- try(sf::st_read(path, quiet = TRUE), silent = TRUE)
  if (inherits(shp, "try-error") || !"Nom" %in% names(shp)) {
    warning("[centroids] could not read the shapefile, or it has no 'Nom' column.", call. = FALSE)
    return(NULL)
  }
  .s2 <- suppressMessages(sf::sf_use_s2())
  on.exit(suppressMessages(try(sf::sf_use_s2(.s2), silent = TRUE)), add = TRUE)
  suppressMessages(try(sf::sf_use_s2(FALSE), silent = TRUE))
  g  <- suppressWarnings(sf::st_point_on_surface(sf::st_make_valid(sf::st_geometry(shp))))
  xy <- sf::st_coordinates(sf::st_transform(g, 4326))
  out <- data.frame(nom = as.character(shp$Nom), lon = xy[, 1], lat = xy[, 2],
                    stringsAsFactors = FALSE)
  out <- out[!is.na(out$nom) & is.finite(out$lon) & is.finite(out$lat), , drop = FALSE]
  .osrm_centroid_cache[[key]] <- out
  message(sprintf("[centroids] %d zone centroids read from %s.", nrow(out), basename(path)))
  out
}

#' Great-circle distance in kilometres (haversine, R = 6371.0088 km). Vectorised.
.haversine_km <- function(lon1, lat1, lon2, lat2) {
  p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
       cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(a)))
}

#' Impute OSRM cells that carry no road route, from great-circle separation.
#'
#' Fits log(cost) ~ log(great-circle km) over the routable off-diagonal pairs and predicts
#' the unroutable ones, so the imputed cost is what a routable pair at the same separation
#' costs on this network. The diagonal is never touched, and a cell whose centroids are
#' unavailable stays NA (the callers' existing non-finite guards then apply as before).
#'
#' @param M    square named OSRM cost matrix (minutes or km).
#' @param kind "travel_time" or "road_distance"; used only in messages.
#' @return M with the imputable NA cells filled.
.fill_osrm_gaps <- function(M, kind = "travel_time") {
  na_mask <- is.na(M)
  diag(na_mask) <- FALSE
  n_gap <- sum(na_mask)
  if (n_gap == 0L) return(M)
  ctr <- load_zone_centroids()
  if (is.null(ctr)) {
    warning(sprintf("[osrm-gapfill] %s: %d unroutable cell(s) left as NA (no centroids).",
                    kind, n_gap), call. = FALSE)
    return(M)
  }
  zn  <- rownames(M)
  idx <- match(zn, ctr$nom)
  if (anyNA(idx))
    warning(sprintf("[osrm-gapfill] %s: %d OSRM zone(s) absent from the shapefile (%s%s); their cells stay NA.",
                    kind, sum(is.na(idx)),
                    paste(utils::head(zn[is.na(idx)], 5), collapse = ", "),
                    if (sum(is.na(idx)) > 5) ", ..." else ""), call. = FALSE)
  lon <- ctr$lon[idx]; lat <- ctr$lat[idx]
  ii  <- seq_along(zn)
  G   <- outer(ii, ii, function(a, b) .haversine_km(lon[a], lat[a], lon[b], lat[b]))
  dimnames(G) <- dimnames(M)
  use <- is.finite(M) & M > 0 & is.finite(G) & G > 0
  diag(use) <- FALSE
  if (sum(use) < 100L) {
    warning(sprintf("[osrm-gapfill] %s: only %d routable pair(s) to calibrate on; cells left as NA.",
                    kind, sum(use)), call. = FALSE)
    return(M)
  }
  fit <- stats::lm(log(as.numeric(M[use])) ~ log(as.numeric(G[use])))
  b   <- stats::coef(fit)
  need <- na_mask & is.finite(G) & G > 0
  M[need] <- exp(b[[1]] + b[[2]] * log(G[need]))
  left <- sum(na_mask & !need)
  message(sprintf(paste0("[osrm-gapfill] %s: imputed %d of %d unroutable cell(s) from great-circle ",
                         "separation (log-log fit on %d routable pairs, R2 = %.3f, cost ~ %.2f * km^%.2f)%s."),
                  kind, sum(need), n_gap, sum(use), summary(fit)$r.squared,
                  exp(b[[1]]), b[[2]],
                  if (left) sprintf("; %d cell(s) had no centroid and stay NA", left) else ""))
  # Name the zones the repair rescued: a zone whose whole column was NA had structurally
  # zero import hazard, which is the failure this exists to remove.
  rescued <- zn[colSums(na_mask) == length(zn) - 1L]
  if (length(rescued))
    message(sprintf("[osrm-gapfill] %s: zones with NO road route at all, now reachable: %s.",
                    kind, paste(rescued, collapse = ", ")))
  M
}

#' Load an OSRM zone-to-zone cost matrix.
#'
#' The file has: first column = `nom` (origin), remaining columns = destination zone names.
#' The diagonal is 0 and some entries are NA (unroutable pairs).
#'
#' @param kind "travel_time" (minutes; default — the deterrence used by every base
#'   mobility kernel) or "road_distance" (kilometres — used by the `-dist` mobility
#'   variants). Both are the SAME square zone x zone form.
#' @return named square numeric matrix (rows/cols = canonical zone names).
load_osrm <- function(kind = c("travel_time", "road_distance")) {
  kind <- match.arg(kind)
  f <- file.path(OSRM_DIR, sprintf("osrm__%s__static.matrix.csv", kind))
  stopifnot(file.exists(f))

  df <- readr::read_csv(f, col_types = readr::cols(.default = "d",
                                                    nom      = "c"),
                         show_col_types = FALSE)
  # First column is 'nom'
  zones <- df$nom
  M     <- as.matrix(df[, -which(colnames(df) == "nom")])
  rownames(M) <- zones
  # Ensure square
  stopifnot(nrow(M) == ncol(M))
  message(sprintf("[osrm] Loaded %dx%d %s matrix (%d NA entries).",
                  nrow(M), ncol(M), gsub("_", "-", kind), sum(is.na(M))))
  # Repair the road-network coverage gaps before ANY kernel sees the matrix, so every
  # consumer (M4/M5/M6, the -dist twins, d_min, the Distance-B1 baseline) gets the same
  # repaired costs. See the block above .fill_osrm_gaps() for why this is a fix and not a
  # fudge. OSRM_GAP_FILL = FALSE restores the previous behaviour exactly.
  if (isTRUE(get0("OSRM_GAP_FILL", ifnotfound = TRUE))) M <- .fill_osrm_gaps(M, kind)
  M
}

#' Load a Flowminder origin-destination matrix (outflow perspective).
#'
#' WHAT THE VALUES ARE. Both files hold Flowminder's national ESTIMATED RELOCATIONS — monthly
#' changes of home location, NOT trips. The provider's variable list defines each column as
#' "Estimated relocations YYYY_MM-1 to YYYY_MM"; the default March file matches the HDX column
#' est_flows_2026_03 on 99.9% of its cells (compared by value), and the April file is
#' est_flows_2026_04.
#'
#' The file is selected by FLOWMINDER_OD_FILE (00_config.R), which defaults to the March 2026
#' file (historically called the "provincial PDF extract") — 437 zones, every cell a number,
#' suppressed counts written as 0 and therefore indistinguishable from measured zeros. The
#' April 2026 national HDX export ("flowminder__outflow_202604__static.matrix.csv") is a strict
#' superset at 467 zones. It leaves EMPTY both the pairs Flowminder marks "redacted (count <15)"
#' (9,416) and the pairs listed with a BLANK value that month (47,945), and writes the pairs absent
#' from the long table as 0. NA there therefore means "no count released this month", of which only
#' a minority are documented redactions; the meaning of a blank is not documented. build_M3 carries
#' NA as a censoring mask, and the gravity fit's three-state likelihood treats every NA as a
#' redacted 1-14 count — correct for the redactions, an assumption for the blanks. The April file
#' is a sensitivity arm only.
#'
#' @return Raw (unnormalised) square numeric matrix. NA = "no count released for this pair
#'   this month": the documented redactions (count < 15) AND the undocumented blanks.
load_flowminder_od <- function(file = get0("FLOWMINDER_OD_FILE",
                                           ifnotfound = "flowminder__outflow__static.matrix.csv")) {
  f <- file.path(FLOWMINDER_DIR, file)
  stopifnot(file.exists(f))

  df <- readr::read_csv(f, col_types = readr::cols(.default = "d",
                                                    nom      = "c"),
                         show_col_types = FALSE)
  zones <- df$nom
  M     <- as.matrix(df[, -which(colnames(df) == "nom")])
  rownames(M) <- zones
  message(sprintf("[flowminder_od] Loaded %dx%d OD matrix from %s (%d redacted cell(s)).",
                  nrow(M), ncol(M), basename(f), sum(is.na(M))))
  M
}

#' Load a Flowminder short-trip snapshot.
#'
#' Returns a named numeric vector of proportions (converted from %) for
#' each destination zone.  All three origin rows (Bunia/Mongbalu/Rwampara)
#' are identical (pooled cohort); only the first data row is used.
#'
#' @param tag  Date tag string, e.g. "20260524".
#' @return Named numeric vector: names = destination zone, values = proportions [0,1].
load_short_trip_snapshot <- function(tag) {
  f <- file.path(
    FLOWMINDER_ST_DIR,
    sprintf("flowminder_short_trips__outflow_%s__static.matrix.csv", tag)
  )
  stopifnot(file.exists(f))

  df <- readr::read_csv(f, col_types = readr::cols(.default = "d",
                                                    nom      = "c"),
                         show_col_types = FALSE)
  stopifnot("nom" %in% colnames(df))

  # All rows are identical — use only the first
  dest_cols <- setdiff(colnames(df), "nom")
  proportions <- as.numeric(df[1L, dest_cols]) / 100   # % → proportion
  names(proportions) <- dest_cols

  message(sprintf(
    "[short_trips] Tag %s: %d destination zones, first row origin='%s', total prop=%.3f",
    tag, length(proportions), df$nom[1L], sum(proportions, na.rm = TRUE)
  ))
  proportions
}

#' Load a Flowminder COHORT subscriber-day presence snapshot.
#'
#' Distinct from load_short_trip_snapshot(): these are the cohort presence-day
#' matrices (avg subscriber-days per cohort member in each destination zone), NOT
#' the short-trip outflow proportions. Values are DAYS (not %), so they are NOT
#' divided by 100. All origin rows are identical (pooled cohort aggregate); only
#' the first data row is read. The processed matrix is canonicalised on the 519
#' WorldPop spine (see data/flowminder_short_trips/process.py), so names align to
#' zones_all directly. The latest tag for the requested window is used.
#'
#' @param cohort One of the cohort ids (e.g. "ituri", "nk", "tshopo").
#' @param window "followup" (during-outbreak; default) or "prior" (look-back).
#' @param analysis_date As-of date; the latest snapshot tag on or before it is used.
#' @return Named numeric vector: names = destination zone, values = avg presence days.
load_cohort_snapshot <- function(cohort, window = "followup",
                                 analysis_date = get0("ANALYSIS_DATE", ifnotfound = NA)) {
  pat <- sprintf(
    "^flowminder_short_trips__%s_subscriber_days_%s_[0-9]{8}__static\\.matrix\\.csv$",
    cohort, window)
  fs <- list.files(FLOWMINDER_ST_DIR, pattern = pat)
  if (!length(fs))
    stop(sprintf("[cohort] No snapshot for cohort='%s' window='%s' in %s",
                 cohort, window, FLOWMINDER_ST_DIR))
  tags <- sub(sprintf(
    "^flowminder_short_trips__%s_subscriber_days_%s_([0-9]{8})__static\\.matrix\\.csv$",
    cohort, window), "\\1", fs)
  # AS-OF SELECTION, matching .read_cohort_raw()/build_M1/build_M2. This used to be
  # which.max(as.integer(tags)) — always the NEWEST snapshot regardless of the analysis date.
  # This function is the FALLBACK the cohort kernels take whenever the raw subscriber-days
  # table is absent or unusable, so on that path every cohort kernel (M13/M14/M16/M17 and
  # their -fill/-split/-dist twins) used post-analysis-date mobility in any back-dated or
  # rolling-origin (LFO) evaluation, leaking future information into the fold.
  .tagd <- suppressWarnings(as.Date(tags, format = "%Y%m%d"))
  .ad   <- suppressWarnings(as.Date(analysis_date))
  .ok   <- if (length(.ad) == 1L && !is.na(.ad) && any(!is.na(.tagd)))
             !is.na(.tagd) & .tagd <= .ad else rep(TRUE, length(fs))
  if (!any(.ok)) {
    # Back-dated run that predates every snapshot: take the EARLIEST, which leaks least.
    .pick <- which.min(ifelse(is.na(.tagd), Inf, as.numeric(.tagd)))
    warning(sprintf("[cohort] %s/%s: no snapshot at or before %s; using the earliest (%s).",
                    cohort, window, format(.ad), fs[.pick]), call. = FALSE)
  } else {
    .cand <- which(.ok)
    # which.max() returns integer(0) when every candidate date is NA (an 8-digit tag that is
    # not a real date), which would silently index fs with nothing. Fall back to lexical order.
    .pick <- if (any(!is.na(.tagd[.cand]))) .cand[which.max(as.numeric(.tagd[.cand]))]
             else .cand[which.max(as.integer(tags[.cand]))]
  }
  f <- file.path(FLOWMINDER_ST_DIR, fs[.pick])
  stopifnot(length(f) == 1L, file.exists(f))

  df <- readr::read_csv(f, col_types = readr::cols(.default = "d", nom = "c"),
                        show_col_types = FALSE)
  stopifnot("nom" %in% colnames(df))
  dest_cols <- setdiff(colnames(df), "nom")
  # All origin rows are the same cohort aggregate — use only the first data row.
  vals <- as.numeric(df[1L, dest_cols])
  names(vals) <- dest_cols
  message(sprintf(
    "[cohort] %s/%s (tag %s): %d destinations, first origin='%s', total presence-days=%.3f",
    # tags[.pick], NOT max(tags): report the snapshot actually loaded, which under as-of
    # selection on a back-dated run is not the newest one on disk.
    cohort, window, tags[.pick], length(vals), df$nom[1L], sum(vals, na.rm = TRUE)
  ))
  vals
}

#' Load IDP weekly flow matrix and aggregate to a static OD matrix.
#'
#' @param aliases  Named alias vector (from load_aliases()).
#' @param zones_all  Canonical zone vector.
#' @return Numeric matrix of summed IDP flows, dim = (n_idp_zones × n_idp_zones).
#'         Rownames and colnames are canonical zone names (post-harmonisation).
load_idp_static <- function(aliases, zones_all) {
  f <- file.path(IDP_DIR, "idp__individuals__weekly.matrix.csv")
  stopifnot(file.exists(f))

  df <- readr::read_csv(f, col_types = readr::cols(date = "c", nom = "c",
                                                    .default = "d"),
                         show_col_types = FALSE)

  # Aggregate over all dates to get a static matrix
  dest_cols <- setdiff(colnames(df), c("date", "nom"))
  static_df <- df %>%
    dplyr::select(-date) %>%
    dplyr::group_by(nom) %>%
    dplyr::summarise(dplyr::across(everything(), \(x) sum(x, na.rm = TRUE)), .groups = "drop")

  origins <- harmonise_names(static_df$nom, aliases, zones_all)
  dests   <- harmonise_names(dest_cols,     aliases, zones_all)

  # Build matrix
  M_raw <- as.matrix(static_df[, dest_cols])
  rownames(M_raw) <- origins
  colnames(M_raw) <- dests

  # AGGREGATE duplicate canonical names, exactly as build_M3() and .embed_od_to_519() do.
  # Harmonisation can map two raw labels onto one canonical zone; without this the block
  # assignment in build_M7 (`W_idp[oz, dz] <- idp_raw[oz, dz]`) keeps only the FIRST match and
  # silently discards the other's entire flow — the reason the other two loaders were hardened.
  # Inert on the current 44x44 file (no duplicates), but data/ is synced from an upstream repo.
  if (anyDuplicated(rownames(M_raw))) {
    message(sprintf("[idp] collapsing %d duplicate origin name(s) by summation.",
                    sum(duplicated(rownames(M_raw)))))
    M_raw <- rowsum(M_raw, group = rownames(M_raw), reorder = FALSE)
  }
  if (anyDuplicated(colnames(M_raw))) {
    message(sprintf("[idp] collapsing %d duplicate destination name(s) by summation.",
                    sum(duplicated(colnames(M_raw)))))
    M_raw <- t(rowsum(t(M_raw), group = colnames(M_raw), reorder = FALSE))
  }

  # The observation window, so a consumer can put these CUMULATIVE counts on a rate basis.
  # Summed over every week in the file they are not comparable with a one-month flow count
  # (see build_M7, which divides by n_months).
  .d <- suppressWarnings(as.Date(df$date))
  .d <- .d[is.finite(.d)]
  n_weeks  <- length(unique(df$date))
  n_months <- if (length(.d) > 1L)
    max(1, as.numeric(difftime(max(.d), min(.d), units = "days")) / 30.44) else NA_real_
  attr(M_raw, "n_weeks")  <- n_weeks
  attr(M_raw, "n_months") <- n_months
  attr(M_raw, "date_range") <- if (length(.d)) range(.d) else as.Date(c(NA, NA))
  message(sprintf("[idp] Loaded IDP static matrix: %d origins × %d destinations (%d weeks%s).",
                  nrow(M_raw), ncol(M_raw), n_weeks,
                  if (length(.d)) sprintf(", %s to %s", min(.d), max(.d)) else ""))
  M_raw
}

# ---------------------------------------------------------------------------
# Source-label resolution for the RAW short-trip / cohort tables
# ---------------------------------------------------------------------------
# WHY THE PIPELINE READS THE RAW TABLES RATHER THAN THE PROCESSED MATRICES.
# The processed short-trip and cohort matrices are canonicalised UPSTREAM (in
# data/flowminder_short_trips/process.py) against the HEADER OF THE 437-ZONE
# Flowminder OD matrix — not against the 519-zone WorldPop spine this pipeline
# runs on. Every destination whose name is outside those 437 is dropped, which
# costs the Ituri cohort 6.97% of its follow-up presence mass and the short-trip
# annex 7.15% of its D+31 mass — in both cases almost entirely ONE zone, Kilo,
# the annex's rank-4 destination, itself invaded on 2026-05-12 and reachable
# with weight exactly zero from Bunia in M1/M8/M13/M14. data/ is synced from the
# public data repo (see run_spatiotemporal.sh PUBLIC_PATHS) and is NOT versioned
# here, so patching process.py alone would be undone by the next sync. The
# pipeline therefore resolves the RAW source tables itself, against zones_all.
# The processed matrices remain the fallback when a raw table is absent.

#' Accent/punctuation-insensitive comparison key (province matching).
.norm_key <- function(x) {
  x <- tolower(as.character(x))
  x <- chartr("àáâãäåèéêëìíîïòóôõöùúûüýÿçñ",
              "aaaaaaeeeeiiiiooooouuuuyycn", x)
  gsub("[^a-z0-9]", "", x)
}

# Source-label spellings that are NOT canonical zone names and are not covered by
# data/aliases.csv (which is itself synced and unversioned here). Each target was
# checked against the health-zone shapefile's PROVINCE field:
#   Massa, Kisantu -> Kongo-Central   (Flowminder province "Kongo")
#   Lukonga        -> Kasai-Central   (Flowminder province "Kasai")
#   Kalambayi Kabanga -> Lomami       (truncated source label)
#   Makiso Kisangani  -> Tshopo       (truncated source label)
# None of these targets appears as a separate row in the source tables, so no
# destination is double-counted by adding them.
SHORT_TRIP_LABEL_ALIASES <- c(
  "Central Massa"   = "Massa",
  "Central Lukonga" = "Lukonga",
  "Central Kisantu" = "Kisantu",
  "Kalambayi Kaban" = "Kalambayi Kabanga",
  "Makiso Kisangan" = "Makiso Kisangani"
)

#' Resolve raw source zone labels to canonical zone names, province-aware.
#'
#' Order: exact canonical match; then the parenthetical province qualifier that
#' the canonical names themselves carry (the only way to separate the two
#' same-named zone pairs in the spine, Bili and Lubunga); then the committed
#' alias table above; then data/aliases.csv. A colliding bare name whose province
#' does not identify exactly one candidate is left UNRESOLVED rather than folded
#' into an arbitrary one — the upstream script folds a Kasai-Central "Lubunga"
#' into "Lubunga (Tshopo)", which is how a Kisangani->Lubunga weight of 4.9e-07
#' reached M13/M14/M16/M17.
#'
#' @param labels    character vector of raw zone labels.
#' @param provinces character vector of the same length (or NULL).
#' @param zones_all canonical zone vector.
#' @param aliases   named vector from load_aliases() (observed -> canonical), optional.
#' @return character vector of canonical names, NA where unresolved.
.resolve_source_labels <- function(labels, provinces, zones_all, aliases = NULL) {
  labels <- trimws(as.character(labels))
  provinces <- if (is.null(provinces)) rep(NA_character_, length(labels))
               else trimws(as.character(provinces))
  stopifnot(length(provinces) == length(labels))

  zone_keys <- .norm_key(zones_all)
  if (anyDuplicated(zone_keys))
    zone_keys[zone_keys %in% zone_keys[duplicated(zone_keys)]] <- NA_character_
  paren    <- grepl(" \\([^()]*\\)$", zones_all)
  zp       <- zones_all[paren]
  base_of  <- sub(" \\([^()]*\\)$", "", zp)
  prov_of  <- .norm_key(sub("^.*\\(([^()]*)\\)$", "\\1", zp))

  out <- rep(NA_character_, length(labels))
  for (i in seq_along(labels)) {
    lab <- labels[i]
    if (is.na(lab) || !nzchar(lab)) next
    if (lab %in% zones_all) { out[i] <- lab; next }
    k <- which(base_of == lab)
    if (length(k)) {
      pk  <- .norm_key(provinces[i])
      hit <- if (is.na(pk) || !nzchar(pk)) integer(0)
             else k[!is.na(prov_of[k]) & prov_of[k] == pk]
      if (length(hit) == 1L) out[i] <- zp[hit]
      next                     # ambiguous / unmatched province -> leave unresolved
    }
    if (lab %in% names(SHORT_TRIP_LABEL_ALIASES)) {
      a <- unname(SHORT_TRIP_LABEL_ALIASES[[lab]])
      if (a %in% zones_all) { out[i] <- a; next }
    }
    if (!is.null(aliases) && lab %in% names(aliases)) {
      a <- unname(aliases[[lab]])
      if (a %in% zones_all) { out[i] <- a; next }
    }
    # Last resort: a case/punctuation-insensitive match against the spine, used only
    # when it is UNAMBIGUOUS. The 519 canonical names have 519 distinct normalised
    # keys, so this can never merge two zones; it recovers labels that differ from
    # the canonical spelling only in case or punctuation ("NGandajika" for
    # "Ngandajika"), which the exact-match path would otherwise drop.
    nk <- .norm_key(lab)
    hit2 <- which(zone_keys == nk)
    if (length(hit2) == 1L) out[i] <- zones_all[hit2]
  }
  out
}

#' Directory holding the RAW short-trip / cohort tables.
.st_raw_dir <- function() file.path(dirname(FLOWMINDER_ST_DIR), "raw")

#' Sum values onto canonical zone names (duplicates ADD, matching the model layer).
.sum_by_zone <- function(values, canon) {
  keep <- !is.na(canon) & !is.na(values)
  if (!any(keep)) return(stats::setNames(numeric(0), character(0)))
  agg <- tapply(as.numeric(values[keep]), canon[keep], sum)
  stats::setNames(as.numeric(agg), names(agg))
}

#' Warn when the PROCESSED short-trip/cohort matrices are missing destinations that
#' the raw tables do resolve.
#'
#' data/ is synced from the public data repo, so an upstream process.py that still
#' canonicalises against the 437-zone OD header will quietly reinstate the drops
#' whenever the data are refreshed. The model layer is insulated (it reads the raw
#' tables), but the processed matrices ship to other consumers, so a divergence is
#' worth reporting. Diagnostic only: never alters a kernel.
.warn_if_processed_source_stale <- function(zones_all, aliases, tag) {
  cmp <- function(raw_names, proc_names, what) {
    proc_can <- harmonise_names(proc_names, aliases, zones_all)
    missing  <- setdiff(raw_names, proc_can)
    if (length(missing))
      warning(sprintf(paste0("[source-check] the PROCESSED %s matrix is missing %d destination(s) ",
                             "that the raw table resolves (%s%s). The model layer uses the raw ",
                             "table and is unaffected; data/flowminder_short_trips/process.py is ",
                             "stale upstream."),
                      what, length(missing), paste(utils::head(missing, 5), collapse = ", "),
                      if (length(missing) > 5) ", ..." else ""), call. = FALSE)
  }
  # The loaders warn about their own unresolved labels during the real build; this
  # diagnostic re-reads the same tables, so silence the duplicate warnings here.
  tryCatch({
    pr <- suppressWarnings(load_short_trip_profile(tag, zones_all, aliases))
    if (identical(pr$source, "raw"))
      cmp(pr$measured, names(suppressMessages(load_short_trip_snapshot(tag))),
          sprintf("short-trip (%s)", tag))
  }, error = function(e) invisible(NULL))
  for (ch in names(get0("COHORT_SOURCES", ifnotfound = list()))) {
    tryCatch({
      pr <- suppressWarnings(load_cohort_profile(ch, get0("COHORT_WINDOW", ifnotfound = "followup"),
                                                 zones_all, aliases))
      if (identical(pr$source, "raw"))
        cmp(pr$measured,
            names(suppressMessages(load_cohort_snapshot(ch, get0("COHORT_WINDOW", ifnotfound = "followup")))),
            sprintf("cohort %s", ch))
    }, error = function(e) invisible(NULL))
  }
  invisible(NULL)
}

#' Short-trip annex profile for one snapshot tag, canonicalised on zones_all.
#'
#' @return list(props = named numeric PROPORTIONS (percent/100) over canonical
#'   destinations, measured = those destination names, source = "raw"/"processed",
#'   dropped = unresolved raw labels).
load_short_trip_profile <- function(tag, zones_all, aliases = NULL) {
  f <- file.path(.st_raw_dir(), "short_trips_destination_rankings.csv")
  if (file.exists(f)) {
    df <- tryCatch(readr::read_csv(f, col_types = readr::cols(.default = "c"),
                                   show_col_types = FALSE),
                   error = function(e) NULL)
    if (!is.null(df) && all(c("province", "health_zone") %in% names(df))) {
      # Map the YYYYMMDD tag to its value column through the date_* columns, so a
      # renamed or newly added snapshot column needs no hard-coded lookup.
      vcol <- NA_character_
      for (dc in grep("^date_", names(df), value = TRUE)) {
        dv <- suppressWarnings(as.Date(df[[dc]][!is.na(df[[dc]])][1]))
        if (length(dv) && !is.na(dv) && format(dv, "%Y%m%d") == tag) {
          vcol <- sub("^date_", "", dc); break
        }
      }
      if (!is.na(vcol) && vcol %in% names(df)) {
        val <- suppressWarnings(as.numeric(df[[vcol]]))
        can <- .resolve_source_labels(df$health_zone, df$province, zones_all, aliases)
        dropped <- unique(trimws(as.character(df$health_zone))[is.na(can)])
        if (length(dropped))
          warning(sprintf("[short_trips] %s: %d raw label(s) unresolved and dropped: %s",
                          tag, length(dropped), paste(utils::head(dropped, 10), collapse = ", ")),
                  call. = FALSE)
        props <- .sum_by_zone(val / 100, can)
        message(sprintf("[short_trips] %s (raw): %d destinations, total prop=%.3f%s",
                        tag, length(props), sum(props),
                        if (length(dropped)) sprintf(", %d dropped", length(dropped)) else ""))
        return(list(props = props, measured = names(props), source = "raw",
                    dropped = dropped))
      }
    }
    warning(sprintf("[short_trips] %s: raw annex table unusable; falling back to the processed matrix (upstream name drops persist).",
                    tag), call. = FALSE)
  } else {
    warning("[short_trips] raw annex table absent; falling back to the processed matrices (upstream name drops persist).",
            call. = FALSE)
  }
  # ---- fallback: the processed matrix (already canonicalised upstream) ----
  p   <- load_short_trip_snapshot(tag)
  can <- harmonise_names(names(p), aliases, zones_all)
  can[!can %in% zones_all] <- NA_character_
  props <- .sum_by_zone(p, can)
  list(props = props, measured = names(props), source = "processed",
       dropped = names(p)[is.na(can)])
}

#' Cohort presence profile for one cohort/window, canonicalised on zones_all.
#'
#' @return list(values = named numeric PRESENCE DAYS over canonical destinations,
#'   measured = those names, origins = the cohort's own origin zones (which the
#'   provider states are NOT captured: "movements within and between those zones
#'   are not captured here"), nodata = zones the provider excluded for
#'   insufficient data, source, dropped).
load_cohort_profile <- function(cohort, window = "followup", zones_all,
                                aliases = NULL,
                                analysis_date = get0("ANALYSIS_DATE", ifnotfound = NA)) {
  fs <- list.files(.st_raw_dir(),
                   pattern = sprintf("^drc-bvd_%s-cohort_subscriber-days-.*external\\.csv$",
                                     cohort),
                   full.names = TRUE)
  if (length(fs)) {
    # AS-OF SELECTION, matching build_M1/build_M2's short-trip snapshot logic. This used to be
    # sort(fs)[length(fs)] — always the NEWEST release — so every cohort kernel (M13/M14/M16/M17
    # and their -fill/-split/-dist twins, i.e. the default family and CASCADE_KERNEL's fallback)
    # used post-analysis-date mobility in any back-dated or rolling-origin evaluation. It also
    # assumed lexical order equals chronological order of the embedded release date.
    .rel <- as.Date(sub(".*subscriber-days-([0-9]{4}[-_][0-9]{2}[-_][0-9]{2}).*", "\\1",
                        basename(fs)) |> gsub(pattern = "_", replacement = "-"),
                    format = "%Y-%m-%d")
    .ad <- suppressWarnings(as.Date(analysis_date))
    .ok <- if (length(.ad) == 1L && !is.na(.ad) && any(!is.na(.rel))) !is.na(.rel) & .rel <= .ad
           else rep(TRUE, length(fs))
    if (!any(.ok)) {
      # Back-dated run that predates every release: take the EARLIEST, which leaks least.
      .pick <- which.min(ifelse(is.na(.rel), Inf, as.numeric(.rel)))
      warning(sprintf("[cohort] %s: no release at or before %s; using the earliest (%s).",
                      cohort, format(.ad), basename(fs[.pick])), call. = FALSE)
    } else {
      .cand <- which(.ok)
      .pick <- if (any(!is.na(.rel[.cand]))) .cand[which.max(as.numeric(.rel[.cand]))]
               else .cand[length(.cand)]
    }
    f  <- fs[.pick]
    df <- tryCatch(readr::read_csv(f, col_types = readr::cols(.default = "c"),
                                   show_col_types = FALSE),
                   error = function(e) NULL)
    vcol <- if (identical(window, "prior")) "Avg days (prior)" else "Avg days (follow-up)"
    if (!is.null(df) && all(c("Province", "Health Zone", vcol) %in% names(df))) {
      flag <- if ("Flag" %in% names(df)) trimws(as.character(df$Flag)) else rep("", nrow(df))
      flag[is.na(flag)] <- ""
      val  <- suppressWarnings(as.numeric(df[[vcol]]))
      can  <- .resolve_source_labels(df[["Health Zone"]], df$Province, zones_all, aliases)
      dropped <- unique(trimws(as.character(df[["Health Zone"]]))[is.na(can)])
      if (length(dropped))
        warning(sprintf("[cohort] %s/%s: %d raw label(s) unresolved and dropped: %s",
                        cohort, window, length(dropped),
                        paste(utils::head(dropped, 10), collapse = ", ")), call. = FALSE)
      is_origin <- flag == "Origin"
      is_meas   <- !is_origin & !is.na(val) & !is.na(can)
      values  <- .sum_by_zone(ifelse(is_meas, val, NA_real_), can)
      origins <- unique(stats::na.omit(can[is_origin]))
      nodata  <- setdiff(unique(stats::na.omit(can[!is_origin & is.na(val)])), names(values))
      message(sprintf("[cohort] %s/%s (raw): %d measured destinations (total %.3f days), %d origin zone(s), %d no-data zone(s)%s",
                      cohort, window, length(values), sum(values), length(origins),
                      length(nodata),
                      if (length(dropped)) sprintf(", %d dropped", length(dropped)) else ""))
      return(list(values = values, measured = names(values), origins = origins,
                  nodata = nodata, source = "raw", dropped = dropped))
    }
    warning(sprintf("[cohort] %s/%s: raw table unusable; falling back to the processed matrix (upstream name drops persist).",
                    cohort, window), call. = FALSE)
  } else {
    warning(sprintf("[cohort] %s: raw table absent; falling back to the processed matrix (upstream name drops persist).",
                    cohort), call. = FALSE)
  }
  # ---- fallback: the processed matrix ----
  v   <- load_cohort_snapshot(cohort, window)
  can <- harmonise_names(names(v), aliases, zones_all)
  can[!can %in% zones_all] <- NA_character_
  values <- .sum_by_zone(v, can)
  list(values = values, measured = names(values), origins = character(0),
       nodata = character(0), source = "processed", dropped = names(v)[is.na(can)])
}

# ---------------------------------------------------------------------------
# M1: Flowminder short trips — static (D+31, tag 20260524)
# ---------------------------------------------------------------------------

#' Build M1: epicentre short-trip matrix using the most recent (D+31) snapshot.
#'
#' For epicentre zone rows (Bunia, Mongbalu, Rwampara): fill from the pooled
#' short-trip proportions.  For all other origin zones: use fallback_M3 row
#' (if provided) or leave as zero.
#'
#' @param zones_all       Canonical zone vector (length 519).
#' @param epicentre_zones Character vector of epicentre health zone names.
#' @param aliases         Alias lookup from load_aliases().
#' @param fallback_M3     Optional row-stochastic 519×519 matrix; provides
#'                        non-epicentre rows when M3 is available.
#' @return Row-stochastic 519×519 matrix.
build_M1 <- function(zones_all, epicentre_zones, aliases, fallback_M3 = NULL,
                     analysis_date = get0("ANALYSIS_DATE", ifnotfound = Sys.Date())) {
  # Use the LATEST short-trip snapshot on or before the analysis date. Snapshots are discovered
  # dynamically in config (FLOWMINDER_ST_TAGS/DATES), so a newer snapshot is picked up
  # automatically as it comes online — no hard-coded tag.
  .elig <- which(FLOWMINDER_ST_DATES <= analysis_date)
  tag <- if (length(.elig)) {
    FLOWMINDER_ST_TAGS[.elig[which.max(FLOWMINDER_ST_DATES[.elig])]]
  } else {
    # Every snapshot post-dates the analysis date (a back-dated run). Use the EARLIEST,
    # which leaks least; never which.max over all dates, which would pick the most
    # future-dated snapshot of all.
    FLOWMINDER_ST_TAGS[which.min(FLOWMINDER_ST_DATES)]
  }
  message(sprintf("[M1] Building Flowminder short-trip matrix (latest snapshot <= %s: %s)...",
                  format(analysis_date), tag))
  # Resolve the RAW annex table against zones_all (see load_short_trip_profile:
  # the processed matrix is canonicalised upstream against a 437-zone list and
  # silently loses real spine zones, Kilo among them).
  prof        <- load_short_trip_profile(tag, zones_all, aliases)
  st_props    <- prof$props
  dest_names_can <- names(st_props)
  message(sprintf("[M1] Snapshot %s (%s source): %d canonical destinations.",
                  tag, prof$source, length(dest_names_can)))

  # Initialise with fallback or zeros
  W <- if (!is.null(fallback_M3)) {
    stopifnot(identical(dim(fallback_M3), c(length(zones_all), length(zones_all))))
    fallback_M3
  } else {
    make_zero_matrix(zones_all)
  }

  # Epicentre rows: overwrite with short-trip proportions
  filled_epi <- character(0)
  for (ez in epicentre_zones) {
    ez_can <- harmonise_names(ez, aliases, zones_all)
    if (!ez_can %in% zones_all) {
      warning(sprintf("[M1] Epicentre zone '%s' not found in zones_all; skipping.", ez))
      next
    }
    # Zero out this row first, then fill destinations
    W[ez_can, ] <- 0
    # Names are already canonical and de-duplicated (summed) by the profile loader.
    W[ez_can, dest_names_can] <- st_props[dest_names_can]
    # Ensure no self-loop
    W[ez_can, ez_can] <- 0
    filled_epi <- c(filled_epi, ez_can)
    message(sprintf("[M1] Epicentre '%s': total outflow prop = %.4f", ez_can,
                    sum(W[ez_can, ], na.rm = TRUE)))
  }

  # Row-normalise (normalises epicentre rows that may not sum to exactly 1,
  # and normalises fallback rows that are already stochastic — idempotent)
  W <- make_row_stochastic(W)

  assert_mobility_matrix(W, zones_all, "M1")
  # Carry the source's OBSERVED support and its own origin set, so a composite can
  # tell "the source measured zero here" from "the source could not look here"
  # (compose_epicentre(fill=)). Attached after make_row_stochastic(), which
  # returns a fresh matrix and would drop them.
  if (length(filled_epi)) {
    attr(W, "measured")       <- stats::setNames(rep(list(prof$measured), length(filled_epi)),
                                                 filled_epi)
    attr(W, "source_origins") <- stats::setNames(rep(list(filled_epi), length(filled_epi)),
                                                 filled_epi)
  }
  message(sprintf("[M1] Built: sparsity=%.1f%%",
                  100 * mean(W == 0)))
  W
}

# ---------------------------------------------------------------------------
# M2a / M2b: Flowminder short trips — time-evolving
# ---------------------------------------------------------------------------

#' Build M2a or M2b: time-evolving short-trip matrix.
#'
#' M2b (default): use only the most recent snapshot with date <= analysis_date.
#' M2a:           average all snapshots with date <= analysis_date (equal weights).
#'
#' @param zones_all       Canonical zone vector.
#' @param epicentre_zones Epicentre zone names.
#' @param aliases         Alias lookup.
#' @param analysis_date   Reference date (Date object).
#' @param variant         "M2a" (average) or "M2b" (most recent).
#' @param fallback_M3     Optional fallback matrix for non-epicentre rows.
#' @return Row-stochastic 519×519 matrix.
build_M2 <- function(zones_all, epicentre_zones, aliases,
                     analysis_date = ANALYSIS_DATE,
                     variant       = "M2b",
                     fallback_M3   = NULL) {
  stopifnot(variant %in% c("M2a", "M2b"))
  message(sprintf("[%s] Building time-evolving short-trip matrix (date <= %s)...",
                  variant, analysis_date))

  # Select eligible snapshot tags
  eligible_mask <- FLOWMINDER_ST_DATES <= analysis_date
  if (!any(eligible_mask)) {
    warning(sprintf("[%s] No snapshots available on or before %s; using fallback rows only.",
                    variant, analysis_date))
    # Preserve the fallback (national) structure when no snapshot is eligible, rather
    # than discarding it for an all-zero matrix (consistent with build_M1's fallback).
    base <- if (!is.null(fallback_M3)) fallback_M3 else make_zero_matrix(zones_all)
    return(make_row_stochastic(base))
  }

  eligible_tags  <- FLOWMINDER_ST_TAGS[eligible_mask]
  eligible_dates <- FLOWMINDER_ST_DATES[eligible_mask]

  if (variant == "M2b") {
    # Most recent snapshot
    use_tags <- eligible_tags[which.max(eligible_dates)]
    message(sprintf("[M2b] Using most recent snapshot: %s", use_tags))
  } else {
    # All eligible snapshots (average)
    use_tags <- eligible_tags
    message(sprintf("[M2a] Averaging %d snapshot(s): %s",
                    length(use_tags), paste(use_tags, collapse = ", ")))
  }

  # Load and average proportions across selected snapshots. Each profile is
  # resolved from the RAW annex table against zones_all (load_short_trip_profile).
  all_profs <- lapply(use_tags, load_short_trip_profile,
                      zones_all = zones_all, aliases = aliases)
  all_props <- lapply(all_profs, `[[`, "props")
  # Align on common destination names. NOTE: intersect() DROPS a destination absent from
  # any one snapshot rather than treating it as zero there. A no-op while every annex
  # snapshot carries the same destination list (the raw annex table ranks the same 142
  # zones at all five dates), but a future snapshot with a different destination set
  # would silently lose zones, so the loss is reported.
  dest_sets <- lapply(all_props, names)
  all_dests <- Reduce(intersect, dest_sets)
  .lost <- setdiff(Reduce(union, dest_sets), all_dests)
  if (length(.lost))
    warning(sprintf("[M2] %d destination(s) are absent from at least one snapshot and are DROPPED from the average (not zero-filled): %s",
                    length(.lost), paste(utils::head(.lost, 10), collapse = ", ")), call. = FALSE)
  # EQUAL-WEIGHT the snapshots. The annex proportions do NOT sum to the same total across
  # dates (0.544, 0.746, 0.888, 1.008, 1.059 over the five current snapshots, resolved against
  # zones_all), so a plain rowMeans of the RAW proportions blends the destination profiles with
  # weights proportional to each snapshot's total — the 24 May snapshot would carry roughly twice
  # the weight of the 30 Apr one. Those totals GROW because the annex reports CUMULATIVE reach
  # (the share of the cohort that had visited a zone by each date, non-exclusive across
  # destinations, so a total above 1 is expected); the weighting is an artefact of the reporting
  # window, not a modelling choice, and it is invisible in the result because the row is
  # row-normalised immediately afterwards (so the totals are discarded anyway and only the
  # blended SHAPE survives). Normalising each snapshot to a profile first makes "average
  # the snapshots", as the docstring promises, mean exactly that.
  .prof <- lapply(all_props, function(p) { v <- p[all_dests]; s <- sum(v, na.rm = TRUE)
                                           if (is.finite(s) && s > 0) v / s else v })
  avg_props <- rowMeans(do.call(cbind, .prof))

  # Build matrix in same way as M1 but with averaged proportions. Destination
  # names are already canonical and de-duplicated by load_short_trip_profile().
  dest_names_can <- all_dests

  W <- if (!is.null(fallback_M3)) fallback_M3 else make_zero_matrix(zones_all)

  filled_epi <- character(0)
  for (ez in epicentre_zones) {
    ez_can <- harmonise_names(ez, aliases, zones_all)
    if (!ez_can %in% zones_all) next
    W[ez_can, ] <- 0
    W[ez_can, dest_names_can] <- avg_props[dest_names_can]
    W[ez_can, ez_can] <- 0
    filled_epi <- c(filled_epi, ez_can)
  }

  W <- make_row_stochastic(W)
  assert_mobility_matrix(W, zones_all, variant)
  if (length(filled_epi)) {
    attr(W, "measured")       <- stats::setNames(rep(list(dest_names_can), length(filled_epi)),
                                                 filled_epi)
    attr(W, "source_origins") <- stats::setNames(rep(list(filled_epi), length(filled_epi)),
                                                 filled_epi)
  }
  message(sprintf("[%s] Built: sparsity=%.1f%%", variant, 100 * mean(W == 0)))
  W
}

# ---------------------------------------------------------------------------
# M3: Flowminder full OD — national, row-normalised
# ---------------------------------------------------------------------------

#' Build M3: national Flowminder RELOCATION matrix embedded in the 519-zone space.
#'
#' The values are monthly home-location changes (load_flowminder_od()), not trips. Zones in the
#' OD matrix are matched by canonical name. Zones not present in the OD matrix get zero rows and
#' zero columns — a statement of MISSING DATA, not of no movement — so M3 is never used as a base
#' kernel directly: cover_relocation_od() fills those gaps first. The origins and destinations the
#' table covers are attached as attributes "od_origins" and "od_dests".
#'
#' @param zones_all  Canonical zone vector.
#' @param aliases    Alias lookup.
#' @return Row-stochastic 519×519 matrix.  Also invisibly returns the
#'         raw (unnormalised) matrix as attribute "raw".
build_M3 <- function(zones_all, aliases) {
  message("[M3] Building Flowminder full OD matrix...")

  M_raw <- load_flowminder_od()   # 437x437 (March) or 467x467 (April national)

  # Redacted cells (NA) are counts SUPPRESSED below the privacy threshold, not
  # zero flow. Keep them as an explicit mask and zero the flows, so nothing
  # downstream has to reason about NA, and the gravity fit can still tell a
  # suppressed cell from a measured zero.
  cens_raw <- is.na(M_raw)
  M_raw[cens_raw] <- 0

  # Harmonise names
  orig_names_can <- harmonise_names(rownames(M_raw), aliases, zones_all)
  dest_names_can <- harmonise_names(colnames(M_raw), aliases, zones_all)

  rownames(M_raw) <- orig_names_can
  colnames(M_raw) <- dest_names_can
  dimnames(cens_raw) <- dimnames(M_raw)

  # Aggregate (SUM) duplicate canonical rows/cols rather than dropping all but the
  # first occurrence: harmonise_names() can map several OD zones onto one canonical
  # name, and keeping only the first would silently discard the others' entire flow.
  # rowsum() groups by name and sums; names stay unique and name-indexable below.
  n_dup_rows <- sum(duplicated(rownames(M_raw)))
  n_dup_cols <- sum(duplicated(colnames(M_raw)))
  if (n_dup_rows > 0L) {
    message(sprintf("[M3] Aggregating %d duplicate origin row(s) by summation.", n_dup_rows))
    cens_raw <- rowsum(cens_raw * 1, group = rownames(M_raw), reorder = FALSE)
    M_raw    <- rowsum(M_raw,        group = rownames(M_raw), reorder = FALSE)
  }
  if (n_dup_cols > 0L) {
    message(sprintf("[M3] Aggregating %d duplicate destination column(s) by summation.", n_dup_cols))
    cens_raw <- t(rowsum(t(cens_raw * 1), group = colnames(M_raw), reorder = FALSE))
    M_raw    <- t(rowsum(t(M_raw),        group = colnames(M_raw), reorder = FALSE))
  }
  # After aggregation a cell counts as censored only if it carries no observed
  # flow at all. A group mixing a redacted cell with an observed one is treated as an
  # observation: strictly its total is a LOWER BOUND (observed + something in 1..14), so
  # this slightly understates it. No such group exists with either shipped file (the March
  # file has no NA cells at all), so the approximation is currently inert.
  cens_raw <- (cens_raw > 0) & (M_raw == 0)

  # Coverage report
  od_origins <- rownames(M_raw)
  od_dests   <- colnames(M_raw)
  in_zones_o <- od_origins %in% zones_all
  in_zones_d <- od_dests   %in% zones_all
  message(sprintf(
    "[M3] OD coverage: %d/%d origins and %d/%d destinations in zones_all.",
    sum(in_zones_o), length(od_origins),
    sum(in_zones_d), length(od_dests)
  ))

  missing_zones <- zones_all[!zones_all %in% od_origins]
  if (length(missing_zones) > 0L) {
    message(sprintf("[M3] %d zones not in OD matrix (will have zero rows): e.g. %s",
                    length(missing_zones),
                    paste(head(missing_zones, 5L), collapse = ", ")))
  }

  # Embed into 519×519
  W_full <- make_zero_matrix(zones_all)
  C_full <- matrix(FALSE, length(zones_all), length(zones_all),
                   dimnames = list(zones_all, zones_all))

  shared_orig <- od_origins[in_zones_o]
  shared_dest <- od_dests[in_zones_d]

  W_full[shared_orig, shared_dest] <- M_raw[shared_orig, shared_dest]
  C_full[shared_orig, shared_dest] <- cens_raw[shared_orig, shared_dest]
  diag(W_full) <- 0
  diag(C_full) <- FALSE

  # Store raw for gravity model calibration
  M_raw_519 <- W_full   # raw (unnormalised) 519×519

  W_norm <- make_row_stochastic(W_full)
  assert_mobility_matrix(W_norm, zones_all, "M3")
  message(sprintf("[M3] Built: sparsity=%.1f%%", 100 * mean(W_norm == 0)))

  attr(W_norm, "raw") <- M_raw_519
  attr(W_norm, "censor_mask") <- C_full
  attr(W_norm, "od_origins")  <- shared_orig    # origins the relocation table covers
  attr(W_norm, "od_dests")    <- shared_dest    # destinations the relocation table covers
  if (any(C_full))
    message(sprintf("[M3] %d cell(s) carry an explicit redaction flag (suppressed count).",
                    sum(C_full)))
  W_norm
}

# ---------------------------------------------------------------------------
# M15: Flowminder SYMMETRISED OD kernel (S = O + t(O))
# ---------------------------------------------------------------------------

#' Embed a named (unnormalised) OD flow table into the 519-zone canonical space.
#'
#' Harmonises origin/destination names to the canonical spine, SUMS any duplicate
#' canonical rows/columns (so several observed zones mapping to one canonical name
#' keep their combined flow, never dropped), embeds into a zero 519x519 matrix by
#' name, and zeroes the diagonal. Returns the RAW (unnormalised) 519x519 matrix.
#' This reproduces build_M3()'s embedding exactly and is shared by build_M15().
#'
#' @param M_raw    Named numeric OD matrix (rows = origin, cols = destination).
#' @param zones_all Canonical zone vector (length n).
#' @param aliases  Alias lookup from load_aliases().
#' @return Raw 519x519 numeric matrix (dimnames = zones_all).
.embed_od_to_519 <- function(M_raw, zones_all, aliases, label = "OD") {
  # ZERO THE CENSORED CELLS FIRST. build_M3() does this explicitly (`cens_raw <- is.na(M_raw);
  # M_raw[cens_raw] <- 0`) and this helper's docstring claims to "reproduce build_M3()'s
  # embedding exactly" — but it carried NA straight through. build_M15() then computes
  # S <- Oe + t(Oe), which propagates each NA to BOTH directions of the pair, and
  # make_row_stochastic() finally turns the result into a hard zero. On the April OD file that
  # silently deleted 1,454 cells carrying a POSITIVE measured flow (41,568 relocations) purely
  # because the reciprocal cell was blank — the "assert zero where the source could not look"
  # failure this repo has fixed once already. Zeroing here also makes the two rowsum() calls
  # below safe, since they run with na.rm = FALSE and would otherwise null a whole
  # duplicate-name group from one NA.
  n_cens <- sum(is.na(M_raw))
  if (n_cens > 0L) {
    M_raw[is.na(M_raw)] <- 0
    message(sprintf("[mobility] %s: %d censored (NA) OD cell(s) treated as unobserved-zero before embedding.",
                    label, n_cens))
  }
  rownames(M_raw) <- harmonise_names(rownames(M_raw), aliases, zones_all)
  colnames(M_raw) <- harmonise_names(colnames(M_raw), aliases, zones_all)
  if (any(duplicated(rownames(M_raw))))
    M_raw <- rowsum(M_raw, group = rownames(M_raw), reorder = FALSE)
  if (any(duplicated(colnames(M_raw))))
    M_raw <- t(rowsum(t(M_raw), group = colnames(M_raw), reorder = FALSE))
  W_full <- make_zero_matrix(zones_all)
  shared_o <- intersect(rownames(M_raw), zones_all)
  shared_d <- intersect(colnames(M_raw), zones_all)
  # Name-indexed block assignment (equivalent to build_M3's per-cell loop, vectorised).
  W_full[shared_o, shared_d] <- M_raw[shared_o, shared_d]
  diag(W_full) <- 0
  W_full
}

#' Build M15: Flowminder SYMMETRISED OD kernel.
#'
#' NOT an inflow-informed kernel. No independent inflow table exists at any stage:
#' the two processed exports are byte-identical (same md5), the 437-zone raw pair
#' is byte-identical too, and the older 101-zone raw pair matches neither identity
#' nor transposition (723 of 858 non-trivial cells equal in position against 395
#' under transposition, and different totals), so it cannot be read as a measured
#' inflow either. What is well defined from one directed table is the SYMMETRISED
#' total two-way flow:
#'
#'   S[j,i] = O[j,i] + O[i,j]   (relocations OUT of j to i + relocations IN to j from i)
#'
#' computed in the 519 canonical space from the outflow/inflow tables of the SAME release
#' (the inflow filename is derived from FLOWMINDER_OD_FILE), embedded by name so pairs are
#' matched by zone, never by row position — the raw tables' row and column orders differ.
#' The combination is orientation-aware (see the block below): `Oe + t(Oe)` when the inflow
#' file duplicates the outflow file, `Oe + Ie` when it is a genuine inflow table. Values are
#' estimated RELOCATIONS, not trips. This is distinct from M3 (directed outflow, normalised as
#' is) and materially less sparse: reciprocal edges present in only one direction are
#' filled. If the inflow file is absent, t(Oe) is used as the inflow counterpart
#' (identical here, and the correct symmetrisation regardless).
#'
#' @param zones_all Canonical zone vector.
#' @param aliases   Alias lookup.
#' @return Row-stochastic 519x519 matrix; attribute "raw" is the unnormalised S.
build_M15 <- function(zones_all, aliases) {
  message("[M15] Building Flowminder symmetrised OD kernel S = O + t(O)...")

  Oe <- .embed_od_to_519(load_flowminder_od(), zones_all, aliases, "M15/outflow")

  # The inflow counterpart must come from the SAME release as the outflow table: the
  # filename is derived from FLOWMINDER_OD_FILE rather than hard-coded, otherwise an April
  # outflow would be combined with a March inflow (two different months, and in the wrong
  # orientation) whenever the sensitivity arm is selected.
  f_out <- get0("FLOWMINDER_OD_FILE", ifnotfound = "flowminder__outflow__static.matrix.csv")
  f_in  <- file.path(FLOWMINDER_DIR, sub("outflow", "inflow", basename(f_out), fixed = TRUE))
  Ie <- if (grepl("outflow", basename(f_out), fixed = TRUE) && file.exists(f_in)) {
    df <- readr::read_csv(f_in, col_types = readr::cols(.default = "d", nom = "c"),
                          show_col_types = FALSE)
    Min <- as.matrix(df[, -which(colnames(df) == "nom")]); rownames(Min) <- df$nom
    message(sprintf("[M15] Loaded %dx%d inflow OD table from %s.",
                    nrow(Min), ncol(Min), basename(f_in)))
    .embed_od_to_519(Min, zones_all, aliases, "M15/inflow")
  } else {
    message("[M15] No matching inflow file; the outflow transpose supplies the inflow counterpart.")
    t(Oe)
  }

  # ---- Total two-way exchange S[i,j] = flow(i->j) + flow(j->i) --------------------------
  # ORIENTATION-AWARE (this block previously hard-coded `S <- Oe + t(Ie)`, which is only
  # right when the inflow table happens to share the OUTFLOW orientation). The data
  # README documents the inflow table as "arrivals to ROW zone FROM columns", i.e.
  # Ie[i,j] = flow(j->i), for which the correct combination is Oe + Ie, NOT Oe + t(Ie).
  # Two errors were cancelling: the March PDF-derived inflow file is a BYTE-DUPLICATE of
  # the outflow file (verified: identical MD5), not its transpose, so `Oe + t(Ie)` reduced
  # to the symmetrisation Oe + t(Oe) — the intended object, reached by the wrong route.
  # Either error alone silently collapses M15 onto the plain outflow kernel M3:
  #   * current code + a CORRECT transpose inflow (as in the April HDX pair, where
  #     inflow == t(outflow) exactly) gives Oe + Oe = 2*Oe  -> row-normalises to M3;
  #   * a blind orientation "fix" + the current duplicate file gives Oe + Oe -> also M3.
  # So the relationship is DETECTED rather than assumed, and a duplicate is called out.
  .same <- isTRUE(all.equal(unname(Ie), unname(Oe), tolerance = 1e-9))
  .tsp  <- isTRUE(all.equal(unname(Ie), unname(t(Oe)), tolerance = 1e-9))
  S <- if (.same) {
    warning("[M15] The inflow OD table is IDENTICAL to the outflow table (a duplicate, not ",
            "a transpose), so it carries no independent inflow information. M15 is built as ",
            "the SYMMETRISED OUTFLOW kernel O + t(O); it is NOT a two-way exchange kernel ",
            "informed by measured arrivals. Fix data/flowminder/ upstream to make this real.",
            call. = FALSE)
    Oe + t(Oe)
  } else if (.tsp) {
    message("[M15] Inflow table is the exact transpose of the outflow table -> ",
            "S = Oe + Ie is the symmetrised two-way exchange.")
    Oe + Ie
  } else {
    message("[M15] Inflow table is independent of the outflow table -> combining per the ",
            "documented orientation (arrivals to row zone from columns): S = Oe + Ie.")
    Oe + Ie
  }
  diag(S) <- 0
  # A two-way exchange kernel is symmetric BEFORE row-normalisation, whichever branch ran.
  # If this ever fails the inflow table's orientation is not what the README documents.
  if (!isTRUE(all.equal(unname(S), unname(t(S)), tolerance = 1e-8)))
    warning("[M15] The combined exchange matrix is NOT symmetric — the inflow table's ",
            "orientation does not match the documented 'arrivals to row zone' convention.",
            call. = FALSE)
  # Guard the failure mode both errors lead to: M15 silently becoming a copy of M3.
  if (isTRUE(all.equal(unname(make_row_stochastic(S)),
                       unname(make_row_stochastic(Oe)), tolerance = 1e-9)))
    warning("[M15] The combined kernel row-normalises to the SAME matrix as the plain ",
            "outflow kernel M3 — M15 adds no information and must not be reported as a ",
            "distinct mobility hypothesis.", call. = FALSE)

  W_norm <- make_row_stochastic(S)
  assert_mobility_matrix(W_norm, zones_all, "M15")
  message(sprintf("[M15] Built: sparsity=%.1f%% (M3 directed-outflow sparsity reported separately).",
                  100 * mean(W_norm == 0)))
  attr(W_norm, "raw") <- S
  W_norm
}

# ---------------------------------------------------------------------------
# M4: Gravity model — NB-GLM calibrated on M3
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Left-censored gravity fit (F10)
# ---------------------------------------------------------------------------
# Flowminder suppresses OD counts below a privacy threshold and the processed matrix encodes
# suppression as 0 — verified in the data: 7,100 positive off-diagonal cells with a minimum
# of EXACTLY 15 and ZERO cells anywhere in 1..14. So the zeros among the COVERED pairs (the
# training set is restricted to the table's own destinations) are
# LEFT-CENSORED observations of "somewhere in 0..14", not observed zeros, and fitting them
# with glm.nb()/glm(poisson) is a mis-specified likelihood. Because censoring hits SMALL
# flows, and small flows are disproportionately LONG-RANGE pairs, the uncensored fit's
# deterrence coefficient is biased steeper; row-normalisation cancels a constant scale error
# but NOT that shape error, and M4 feeds M8 = MOBILITY_PRIMARY.
#
# The correct likelihood keeps the observed (y >= k) rows as NB densities and replaces each
# censored row's contribution with the interval probability P(Y <= k-1):
#     logL = sum_{y >= k} log f_NB(y | mu, size) + sum_{censored} log F_NB(k-1 | mu, size)
# with log mu = X beta (the SAME gravity mean structure), maximised by optim() from the
# uncensored fit's coefficients. Falls back to the uncensored fit, loudly, if it fails.
#' @param cens_vec optional logical vector marking the CENSORED training rows. When
#'   supplied (the April national export flags suppressed cells explicitly), the
#'   censored set is exactly those rows and every other zero is an OBSERVED zero.
#'   When NULL, censoring is implicit and every recorded zero is treated as
#'   censored — the only option for the March file, where the two are written the
#'   same way.
#' @param cens_lo lower end of the censoring interval: a censored cell is known to
#'   lie in [cens_lo, k_cens - 1]. 0 for the implicit case (a recorded zero could
#'   be a true zero), 1 when measured zeros are recorded separately.
.fit_gravity_censored <- function(train_df, gravity_formula, k_cens, start_fit, label,
                                  cens_vec = NULL, cens_lo = 0L) {
  tr  <- stats::delete.response(stats::terms(gravity_formula))
  X   <- stats::model.matrix(tr, data = train_df)
  y   <- train_df$flow
  cen <- if (is.null(cens_vec)) y < k_cens else as.logical(cens_vec)
  cen[is.na(cen)] <- FALSE
  unc <- !cen
  p   <- ncol(X)
  # Start from the uncensored fit: same design, so the coefficients are directly comparable.
  b0 <- stats::coef(start_fit)
  b0 <- b0[colnames(X)]; b0[!is.finite(b0)] <- 0
  size0 <- if (!is.null(start_fit$theta) && is.finite(start_fit$theta) && start_fit$theta > 0)
    start_fit$theta else 1
  # log(size) is BOUNDED, not free. Left unbounded, optim probes sizes like 1e111 while
  # exploring; those are finite (so a simple is.finite() guard passes) but drive
  # pnbinom(log.p = TRUE) -> pbeta(a = 1e111, b = 15, x = 1) into an underflow-to--Inf
  # branch, emitting ~180,000 warnings per fit. The likelihood value is still correctly
  # rejected below, so the ESTIMATE was never affected — but that warning flood would bury
  # any real diagnostic in the run log. size in [1e-4, 1e4] spans "wildly overdispersed" to
  # "numerically Poisson" for counts of this magnitude, so the bound cannot exclude a
  # plausible optimum; L-BFGS-B enforces it directly.
  lo <- c(rep(-Inf, p), log(1e-4)); hi <- c(rep(Inf, p), log(1e4))
  # PARAMETER SCALING. The design columns differ by orders of magnitude: under the
  # exp-deterrence variant dist_ij is in MINUTES, so its coefficient is ~1e-3 while the
  # intercept is ~-12. L-BFGS-B uses numeric gradients on the raw parameter scale, and on
  # that geometry the line search fails (convergence 52) with the coefficients still at
  # their starting values. Optimise in a standardised design instead — each non-intercept
  # column centred and scaled to unit sd — and back-transform afterwards. The likelihood
  # is identical; only the search coordinates change.
  is_int  <- apply(X, 2L, function(z) all(z == z[1L]))     # intercept column(s)
  ctr <- ifelse(is_int, 0, colMeans(X))
  scl <- apply(X, 2L, stats::sd); scl[is_int | !is.finite(scl) | scl <= 0] <- 1
  Xs <- X
  for (k in seq_len(p)) if (!is_int[k]) Xs[, k] <- (X[, k] - ctr[k]) / scl[k]
  # b (raw) -> bs (scaled): bs_k = b_k * scl_k for slopes, and the intercept absorbs
  # the centring, bs_int = b_int + sum_k b_k * ctr_k.
  to_scaled <- function(b) {
    bs <- b * scl
    bs[is_int] <- b[is_int] + sum((b * ctr)[!is_int])
    bs
  }
  from_scaled <- function(bs) {
    b <- bs / scl
    b[is_int] <- bs[is_int] - sum((bs * ctr / scl)[!is_int])
    b
  }
  nll <- function(par) {
    b <- par[seq_len(p)]; size <- exp(par[p + 1L])
    if (!is.finite(size) || size <= 0) return(1e12)
    mu <- exp(pmin(pmax(as.numeric(Xs %*% b), -700), 700))
    if (any(!is.finite(mu))) return(1e12)
    l1 <- if (any(unc)) sum(stats::dnbinom(y[unc], mu = mu[unc], size = size, log = TRUE)) else 0
    # Censored rows contribute log P(cens_lo <= Y <= k-1). With cens_lo = 0 that is
    # log F(k-1) exactly as before; with cens_lo = 1 it is log[F(k-1) - F(0)],
    # evaluated as a stable log-difference (log1p form) and floored at a finite
    # value so a numerically impossible interval cannot return NaN to optim().
    l2 <- if (any(cen)) {
      lp_hi <- stats::pnbinom(k_cens - 1L, mu = mu[cen], size = size, log.p = TRUE)
      if (cens_lo <= 0L) sum(lp_hi) else {
        lp_lo <- stats::pnbinom(cens_lo - 1L, mu = mu[cen], size = size, log.p = TRUE)
        d <- lp_hi + log1p(-exp(pmin(lp_lo - lp_hi, 0)))
        sum(pmax(ifelse(is.finite(d), d, -745), -745))
      }
    } else 0
    v <- -(l1 + l2)
    if (!is.finite(v)) 1e12 else v
  }
  t0 <- Sys.time()
  start <- c(to_scaled(b0), log(min(max(size0, 1e-4), 1e4)))
  .run <- function(par, method) {
    tryCatch(suppressWarnings(
      if (identical(method, "L-BFGS-B"))
        stats::optim(par, nll, method = "L-BFGS-B", lower = lo, upper = hi,
                     control = list(maxit = 500L, factr = 1e7))
      else
        stats::optim(par, nll, method = "Nelder-Mead",
                     control = list(maxit = 5000L, reltol = 1e-10))),
      error = function(e) NULL)
  }
  op <- .run(start, "L-BFGS-B")
  # A non-zero convergence code means optim STOPPED WITHOUT CONVERGING (e.g. 52 = line
  # search failure), typically leaving most coefficients at their starting values. That
  # must never be accepted silently: restart from the best point with a derivative-free
  # method and polish, and refuse the fit if it still has not converged.
  if (!is.null(op) && !identical(op$convergence, 0L)) {
    message(sprintf("[%s] Censored fit did not converge (code %d); restarting.",
                    label, op$convergence))
    op2 <- .run(op$par, "Nelder-Mead")
    if (!is.null(op2)) {
      op3 <- .run(pmin(pmax(op2$par, lo), hi), "L-BFGS-B")
      if (!is.null(op3) && is.finite(op3$value) &&
          (is.null(op) || !is.finite(op$value) || op3$value <= op$value)) op <- op3
    }
  }
  if (is.null(op) || !is.finite(op$value) || op$value >= 1e12 ||
      !identical(op$convergence, 0L)) {
    warning(sprintf(paste0("[%s] Censored gravity fit failed to converge (%s); keeping the ",
                           "(biased) uncensored fit."),
                    label,
                    if (is.null(op)) "optim error" else sprintf("code %d", op$convergence)),
            call. = FALSE)
    return(NULL)
  }
  cf <- from_scaled(op$par[seq_len(p)]); names(cf) <- colnames(X)
  message(sprintf("[%s] Censored NB gravity fit: converged=%s, size=%.4f, logL=%.1f (%.0fs)",
                  label, op$convergence == 0L, exp(op$par[p + 1L]), -op$value,
                  as.numeric(Sys.time() - t0, units = "secs")))
  structure(list(coefficients = cf, size = exp(op$par[p + 1L]), terms_rhs = tr,
                 convergence = op$convergence, loglik = -op$value, k_cens = k_cens,
                 cens_lo = cens_lo, explicit_mask = !is.null(cens_vec),
                 n_censored = sum(cen), n_observed = sum(unc)),
            class = "censored_gravity")
}

#' predict() for the censored gravity fit, so the caller's
#' predict(fit, newdata, type = "response") path is unchanged.
predict.censored_gravity <- function(object, newdata, type = "response", ...) {
  X   <- stats::model.matrix(object$terms_rhs, data = newdata)
  eta <- as.numeric(X %*% object$coefficients[colnames(X)])
  if (identical(type, "link")) eta else exp(pmin(pmax(eta, -700), 700))
}

#' Build M4: gravity model calibrated on the Flowminder RELOCATION flows (M3; not trips).
#'
#' Fits an unconstrained gravity GLM (negative-binomial, falling back to Poisson
#' if the NB theta iteration diverges) with log(pop_i), log(pop_j) and
#' log(dist_ij + 1).  Predicts for all 519×519 pairs and row-normalises.  A
#' per-origin intercept would cancel under row-normalisation, so origin size is
#' captured by the single log(pop_i) covariate rather than fixed effects.
#'
#' When the flows are censored (see the diagnostic below), the NB/Poisson fit is used only
#' as a starting point and the kernel comes from a REFIT under the censored likelihood
#' (.fit_gravity_censored); set M4_CENSORED_FIT = FALSE to keep the uncensored fit.
#'
#' @param M3_raw   Raw (unnormalised) 519×519 flow matrix (from attr(M3, "raw")).
#' @param pop_vec  Named population vector (length 519).
#' @param osrm_mat OSRM cost matrix (zones_all × zones_all; minutes, or km for -dist).
#' @param zones_all Canonical zone vector.
#' @param censor_mask Optional logical 519×519 mask of SUPPRESSED cells (April export).
#' @param od_dests Optional character vector of destinations the source OD table covers
#'   (attr(M3, "od_dests")). Training is restricted to them; prediction is not.
#' @return Row-stochastic 519×519 matrix.  Attribute "model_summary" holds the fit used:
#'         the NB/Poisson GLM summary, or the "censored_gravity" object after a refit.
build_M4 <- function(M3_raw, pop_vec, osrm_mat, zones_all,
                     deterrence = c("power", "exp"), label = "M4",
                     censor_mask = NULL, od_dests = NULL) {
  deterrence <- match.arg(deterrence)
  message(sprintf("[%s] Building gravity model (calibrated on M3, %s deterrence)...",
                  label, deterrence))

  n <- length(zones_all)

  # Align pop and OSRM to zones_all
  pop_aligned  <- pop_vec[zones_all]
  pop_aligned[is.na(pop_aligned)] <- 0

  # OSRM: subset and align to zones_all (zones in OSRM matrix may be a subset)
  osrm_zones <- rownames(osrm_mat)
  shared_osrm <- intersect(zones_all, osrm_zones)
  message(sprintf("[%s] OSRM covers %d / %d zones.", label, length(shared_osrm), n))

  # Build a full zones_all × zones_all distance matrix
  dist_mat <- matrix(NA_real_, nrow = n, ncol = n,
                     dimnames = list(zones_all, zones_all))
  dist_mat[shared_osrm, shared_osrm] <- osrm_mat[shared_osrm, shared_osrm]
  diag(dist_mat) <- 0

  # Max observed distance (for imputing NAs later)
  max_dist <- max(dist_mat, na.rm = TRUE)

  # ---- Training data construction ----
  # One row per (i,j) pair where:
  #   - both zones have pop > 0
  #   - OSRM distance > 0 and not NA
  #   - M3_raw[i,j] is the response (including zeros)
  # Restrict to zones that appear as origins in M3 (rows with any non-zero outflow)
  has_outflow <- rowSums(M3_raw, na.rm = TRUE) > 0
  train_zones <- zones_all[has_outflow & pop_aligned > 0]
  message(sprintf("[%s] Training data: %d origin zones with outflow and pop > 0.",
                  label, length(train_zones)))

  # DESTINATION COVERAGE. Zones outside the source OD table have an all-zero column in
  # M3_raw because the source never reported them, not because no one travels there. Left
  # in the training set they enter as measured zeros (explicit-mask arm) or as censored
  # 0..k-1 counts (implicit arm), which biases the mass and deterrence coefficients. Train
  # only on destinations the table could report; prediction still covers every zone.
  dest_pool <- zones_all[pop_aligned > 0]
  if (!is.null(od_dests)) {
    covered <- intersect(dest_pool, od_dests)
    if (length(covered) < 2L) {
      warning(sprintf("[%s] od_dests covers %d usable destination(s); training on all zones.",
                      label, length(covered)), call. = FALSE)
    } else {
      message(sprintf("[%s] Restricting training destinations to the %d zones the OD table covers (of %d).",
                      label, length(covered), length(dest_pool)))
      dest_pool <- covered
    }
  }

  train_rows <- lapply(train_zones, function(oz) {
    for_dests <- dest_pool[dest_pool != oz]
    dists     <- dist_mat[oz, for_dests]
    valid     <- !is.na(dists) & dists > 0
    for_dests <- for_dests[valid]
    dists     <- dists[valid]
    if (length(for_dests) == 0L) return(NULL)
    tibble::tibble(
      origin  = oz,
      dest    = for_dests,
      flow    = as.numeric(M3_raw[oz, for_dests]),
      # TRUE where the source suppressed the count rather than measuring zero.
      cens    = if (is.null(censor_mask)) rep(FALSE, length(for_dests))
                else as.logical(censor_mask[oz, for_dests]),
      pop_i   = pop_aligned[oz],
      pop_j   = pop_aligned[for_dests],
      dist_ij = dists
    )
  })
  train_df <- dplyr::bind_rows(train_rows)
  train_df$flow <- pmax(round(train_df$flow), 0L)   # ensure non-negative integer
  train_df$cens[is.na(train_df$cens)] <- FALSE

  message(sprintf("[%s] Training set: %d OD pairs (%.1f%% zero flows).",
                  label, nrow(train_df), 100 * mean(train_df$flow == 0)))

  # ---- LEFT-CENSORING DIAGNOSTIC (do not silence) ----------------------------------
  # Flowminder suppresses OD counts below a privacy threshold and the processed matrix
  # encodes suppression as 0. Detect it from the data itself: if there is a positive
  # minimum k > 1 and NO observed value in 1..k-1, the zeros are not observations of
  # "no travel" but LEFT-CENSORED observations of "somewhere in 0..k-1". Fitting
  # glm.nb()/glm(poisson) to them is then a mis-specified likelihood — every true flow
  # of 1..k-1 enters the fit as an exact 0. Because censoring hits SMALL flows, and
  # small flows are disproportionately LONG-RANGE pairs, the deterrence coefficient is
  # biased steeper and the mass coefficients with it. Row-normalisation cancels a constant
  # scale error but NOT this shape error, and M4 feeds M8 (MOBILITY_PRIMARY), so it reaches
  # the featured model. The statistically correct treatment is a censored likelihood
  # (flow = 0 contributes P(Y <= k-1) rather than P(Y = 0)) — a deliberate modelling change,
  # so it is surfaced here rather than made silently.
  .pos <- train_df$flow[train_df$flow > 0]
  .k_cens <- NA_integer_
  .cens_explicit <- FALSE
  if (length(.pos)) {
    .k <- min(.pos)
    # EXPLICIT censoring: the source marks suppressed cells (April national export),
    # so a recorded 0 is a MEASURED zero and only the flagged cells are intervals,
    # known to lie in [1, k-1]. Needs k > 1 for that interval to be non-empty.
    .cens_explicit <- any(train_df$cens) && .k > 1
    if (.cens_explicit) {
      .k_cens <- as.integer(.k)
      message(sprintf(paste0("[%s] OD training data carry EXPLICIT redaction flags: %d of %d ",
                             "pairs (%.1f%%) are suppressed counts in 1..%d; the remaining ",
                             "%.1f%% of zeros are measured zeros."),
                      label, sum(train_df$cens), nrow(train_df),
                      100 * mean(train_df$cens), .k - 1L,
                      100 * mean(train_df$flow == 0 & !train_df$cens)))
    } else if (.k > 1 && !any(train_df$flow > 0 & train_df$flow < .k)) {
      .k_cens <- as.integer(.k)
      message(sprintf(paste0("[%s] OD training data are LEFT-CENSORED: smallest observed ",
                             "positive flow is %d, no value occurs in 1..%d, so the %.1f%% of ",
                             "cells recorded as 0 are suppressed counts in 0..%d."),
                      label, .k, .k - 1L, 100 * mean(train_df$flow == 0), .k - 1L))
    }
  }

  # ---- Fit gravity GLM ----
  # Unconstrained gravity model: flow ~ pop_i^b1 * pop_j^b2 * (dist+1)^b3.
  # The origin-size term (pop_i) enters as a single covariate rather than 519
  # origin fixed effects: a per-origin intercept is an origin-constant
  # multiplier that cancels exactly under row-normalisation, so fixed effects
  # add ~500 dense design columns and numerical instability for no effect on
  # the final row-stochastic matrix.  Negative-binomial is attempted first
  # (overdispersion), but its theta iteration can diverge on these heavy-tailed,
  # ~97%-zero flows; we fall back to Poisson, whose mean structure is identical
  # in form and fully sufficient for a normalised mobility kernel.
  # Deterrence function: "power" → flow ∝ (dist+1)^b3 (log-distance term, the
  # classic gravity form); "exp" → flow ∝ exp(b3·dist) (linear-distance term, an
  # exponential deterrence kernel; Truscott & Ferguson 2012 PLoS Comput Biol
  # 8:e1002699; Balcan et al. 2009 PNAS 106:21484). Both share the identical
  # pop_i/pop_j mass terms and differ only in how travel time attenuates flow.
  gravity_formula <- if (deterrence == "power")
    flow ~ log(pop_i) + log(pop_j) + log(dist_ij + 1)
  else
    flow ~ log(pop_i) + log(pop_j) + dist_ij
  message(sprintf("[%s] Fitting gravity GLM (NB with Poisson fallback)...", label))
  nb_fit <- tryCatch(
    {
      m <- MASS::glm.nb(gravity_formula, data = train_df,
                        control = glm.control(maxit = 100))
      if (!isTRUE(m$converged)) stop("did not converge")
      message(sprintf("[%s] NB-GLM converged | theta=%.3f", label, m$theta))
      m
    },
    error = function(e) {
      message(sprintf("[%s] glm.nb unstable (%s); falling back to Poisson GLM.",
                      label, conditionMessage(e)))
      pm <- glm(gravity_formula, family = poisson(), data = train_df,
                control = glm.control(maxit = 100))
      message(sprintf("[%s] Poisson GLM converged: %s", label, pm$converged))
      pm
    }
  )

  # ---- Refit under the CENSORED likelihood when the data are censored (F10) ----
  # Uses the uncensored fit above purely as a starting point. Set M4_CENSORED_FIT = FALSE
  # to reproduce the previous (mis-specified, zeros-as-observed) kernel.
  if (!is.na(.k_cens) && isTRUE(get0("M4_CENSORED_FIT", ifnotfound = TRUE))) {
    .cf <- .fit_gravity_censored(train_df, gravity_formula, .k_cens, nb_fit, label,
                                 cens_vec = if (.cens_explicit) train_df$cens else NULL,
                                 cens_lo  = if (.cens_explicit) 1L else 0L)
    if (!is.null(.cf)) {
      .old <- stats::coef(nb_fit); .new <- .cf$coefficients
      .shared <- intersect(names(.old), names(.new))
      message(sprintf("[%s] Deterrence/mass coefficients, uncensored -> censored: %s",
                      label, paste(sprintf("%s %.3f->%.3f", .shared, .old[.shared], .new[.shared]),
                                   collapse = "; ")))
      nb_fit <- .cf
    }
  } else if (!is.na(.k_cens)) {
    warning(sprintf("[%s] M4_CENSORED_FIT = FALSE: the gravity model is fitted with censored zeros treated as observed zeros, which biases the deterrence coefficient steeper.",
                    label), call. = FALSE)
  }

  # ---- Predict for all 519×519 pairs ----
  message(sprintf("[%s] Predicting for all zone pairs...", label))
  W_pred <- make_zero_matrix(zones_all)

  for (oz in zones_all) {
    if (pop_aligned[oz] == 0) next   # no outflow expected from empty zones

    dests    <- zones_all[zones_all != oz & pop_aligned > 0]
    dists    <- dist_mat[oz, dests]

    # Track which destinations are unroutable before imputing NAs
    na_dist_mask <- is.na(dists)

    # Impute NA distances as max + 1 (very distant → near-zero gravity)
    dists[na_dist_mask] <- max_dist + 1

    if (length(dests) == 0L) next

    pred_df <- tibble::tibble(
      origin  = oz,
      dest    = dests,
      pop_i   = pop_aligned[oz],
      pop_j   = pop_aligned[dests],
      dist_ij = dists
    )

    # Predicted expected flow (origin-size term cancels under row-normalisation)
    fitted_vals <- tryCatch(
      suppressWarnings(
        predict(nb_fit, newdata = pred_df, type = "response")
      ),
      error = function(e) {
        # FAIL LOUD AND INERT, never quiet and wrong. This used to return the grand mean of
        # positive flows — a CONSTANT over every destination — so a failed origin came out of
        # make_row_stochastic() radiating UNIFORMLY to all ~517 populated zones. That is exactly
        # the pathology build_M6() documents and guards against ("an OSRM-missing infected origin
        # would radiate invasion force everywhere"), it passed every assertion because a uniform
        # row is a valid distribution, and it printed nothing. M4 feeds M8 = MOBILITY_PRIMARY,
        # so a silent uniform row there propagates into the featured kernel. A ZERO row is the
        # honest answer and is the case make_row_stochastic() already handles, exactly as M4/M5/M6
        # do for unroutable origins.
        warning(sprintf("[mobility] %s: gravity prediction failed for origin '%s' (%s); row left at ZERO rather than uniform.",
                        label, oz, conditionMessage(e)), call. = FALSE)
        rep(0, length(dests))
      }
    )
    # Sanitise predictions. A +Inf response under the exp-deterrence variant (M4b) is a
    # genuine HIGH-flow pair whose linear predictor overflowed (large populations, short
    # distance) — NOT a zero-flow pair — so cap it to the largest finite flow rather than
    # zeroing it (which would invert the kernel and collapse a dominant edge to nothing).
    # NA/negative predictions become zero flow.
    fitted_vals[is.na(fitted_vals) | fitted_vals < 0] <- 0
    .fin_flow <- fitted_vals[is.finite(fitted_vals)]
    fitted_vals[is.infinite(fitted_vals) & fitted_vals > 0] <-
      if (length(.fin_flow)) max(.fin_flow) else 0
    fitted_vals[!is.finite(fitted_vals)] <- 0
    fitted_vals[na_dist_mask] <- 0   # zero out truly unroutable pairs

    for (i in seq_along(dests)) {
      W_pred[oz, dests[i]] <- fitted_vals[i]
    }
  }
  diag(W_pred) <- 0

  W_norm <- make_row_stochastic(W_pred)
  assert_mobility_matrix(W_norm, zones_all, label)
  message(sprintf("[%s] Built: sparsity=%.1f%%", label, 100 * mean(W_norm == 0)))

  # summary() has no method for the censored fit, so carry the fit itself in that case.
  attr(W_norm, "model_summary") <- if (inherits(nb_fit, "censored_gravity")) nb_fit
                                   else tryCatch(summary(nb_fit), error = function(e) NULL)
  # The FITTED COEFFICIENTS, in a form that needs no summary()-parsing. build_M4c() borrows
  # log(pop_j) from here: the cohort tables identify the deterrence but not the mass exponent,
  # so M4 and M4c differ in exactly one parameter by construction.
  attr(W_norm, "coefficients") <- tryCatch(stats::coef(nb_fit), error = function(e) NULL)
  W_norm
}

# ---------------------------------------------------------------------------
# M4c: gravity with the deterrence calibrated on Flowminder COHORT PRESENCE
# ---------------------------------------------------------------------------
# WHY THIS KERNEL EXISTS. M4's deterrence is fitted to M3, Flowminder's national RELOCATION
# table -- monthly changes of home location. People who move house move far, so that source is
# rich in LONG-range signal and gives a shallow decay (log(d+1) coefficient -1.44 on the
# shipped fit). The cohort tables measure something else: average presence-DAYS per subscriber.
# Time is overwhelmingly spent at or near home, so those data are rich in SHORT-range signal
# and give a much steeper decay (-2.69). Neither is wrong; each is the best available view of
# its own slice of the mobility spectrum, and neither can see the other's. Averaging them is
# the honest response to not knowing which slice governs importation.
#
# WHY A FITTED KERNEL RATHER THAN A KERNEL MEMBER. The cohort tables cover 10 origin rows, so
# they cannot enter a kernel average AS ROWS: a member with 10 non-zero rows of 519 drops out
# of the other 509 under .consensus_base(), reproducing exactly the coverage-driven weight
# drift that the raw-M3 member used to cause (and worse -- 10 rows rather than 413). And on
# the 10 rows where it IS non-zero, compose_epicentre() overwrites the average anyway. Cohort
# data generalise as a PARAMETER, not as rows: fit the deterrence here, and all 519 origin
# rows carry cohort information at equal weight.
#
# WHAT IS FITTED AND WHAT IS BORROWED. Only the deterrence is fitted. The destination-mass
# exponent is taken from M4 (fitted on ~400 origins nationally), because three pooled cohort
# profiles cannot identify it -- per-cohort estimates are -0.03 (Ituri), +1.83 (NK), +0.15
# (Tshopo), and the two regressors are NOT collinear (within-cohort corr 0.03), so that spread
# is genuine between-cohort heterogeneity with n = 3, not an estimation artefact. Fixing it
# costs 0.016 on the deterrence, a third of a standard error.
#
# HONEST UNCERTAINTY. The regression SE (~0.055) is a WITHIN-cohort standard error over ~1,350
# destinations that are not independent -- they are three profiles. The BETWEEN-cohort spread,
# about -2.1 (NK) to -3.1 (Tshopo), is the honest measure; fit_cohort_gravity() reports it and
# it is what should be quoted. The follow-up and look-back windows agree (-2.69 vs -2.66), so
# the steep decay is NOT an artefact of outbreak response or movement restriction.

#' A zones_all x zones_all cost matrix from an OSRM matrix that may cover only a subset.
#' Unrepresented pairs stay NA; the diagonal is 0. Mirrors build_M4/build_M5's alignment.
.align_dist <- function(osrm_mat, zones_all) {
  n  <- length(zones_all)
  D  <- matrix(NA_real_, n, n, dimnames = list(zones_all, zones_all))
  sh <- intersect(intersect(zones_all, rownames(osrm_mat)), colnames(osrm_mat))
  if (length(sh)) D[sh, sh] <- osrm_mat[sh, sh]
  diag(D) <- 0
  D
}

#' Fit the gravity deterrence to Flowminder cohort presence.
#'
#' One observation per (cohort, measured destination): the cohort's row-normalised presence
#' share, its destination population, and the POPULATION-WEIGHTED travel time from the
#' cohort's origin zones (the cohort is pooled over its origins, so no single origin distance
#' exists). Fitted by quasi-Poisson with a log link and cohort fixed effects -- quasi-likelihood
#' because the response is a continuous share with exact zeros (the measured zeros are real
#' observations and belong in the fit), and fixed effects because each cohort's overall level
#' is a nuisance that cancels under row-normalisation anyway.
#'
#' @param M_cohort       build_M_cohort() output (carries "measured").
#' @param cohort_sources COHORT_SOURCES.
#' @param osrm_mat       cost matrix (travel time, or road km for the -dist twin).
#' @param b_pop          destination-mass exponent to HOLD FIXED (M4's); NULL fits it freely.
#' @return list(b_pop, b_dist, per_cohort, n, n_cohorts, mass_fixed, pseudo_r2), or NULL.
fit_cohort_gravity <- function(M_cohort, cohort_sources, pop_vec, osrm_mat, zones_all, aliases,
                               b_pop = NULL, label = "M4c") {
  measured <- attr(M_cohort, "measured")
  if (is.null(measured)) {
    warning(sprintf("[%s] the cohort kernel carries no 'measured' attribute; cannot fit.",
                    label), call. = FALSE)
    return(NULL)
  }
  D   <- .align_dist(osrm_mat, zones_all)
  pop <- pop_vec[zones_all]; names(pop) <- zones_all; pop[is.na(pop)] <- 0
  rows <- list()
  for (ch in names(cohort_sources)) {
    O <- intersect(harmonise_names(cohort_sources[[ch]], aliases, zones_all), zones_all)
    O <- intersect(O, names(measured))
    if (!length(O)) next
    A <- setdiff(intersect(measured[[O[1]]], zones_all), O)
    if (length(A) < 25L) next
    w <- as.numeric(pop[O])
    if (!all(is.finite(w)) || sum(w) <= 0) w <- rep(1, length(O))
    w <- w / sum(w)
    rows[[ch]] <- data.frame(cohort = ch, dest = A,
                             share = as.numeric(M_cohort[O[1], A]),
                             pop_j = as.numeric(pop[A]),
                             tt    = as.numeric(w %*% D[O, A, drop = FALSE]),
                             stringsAsFactors = FALSE)
  }
  df <- dplyr::bind_rows(rows)
  if (!nrow(df)) {
    warning(sprintf("[%s] no cohort produced usable training rows.", label), call. = FALSE)
    return(NULL)
  }
  n_raw <- nrow(df)
  df <- df[is.finite(df$share) & is.finite(df$tt) & df$tt > 0 &
           is.finite(df$pop_j) & df$pop_j > 0, , drop = FALSE]
  if (nrow(df) < 50L) {
    warning(sprintf("[%s] only %d usable training row(s) of %d; kernel not built.",
                    label, nrow(df), n_raw), call. = FALSE)
    return(NULL)
  }
  df$.off <- if (is.null(b_pop)) 0 else as.numeric(b_pop)[1] * log(df$pop_j)
  n_ch  <- length(unique(df$cohort))
  base_terms <- c(if (is.null(b_pop)) "log(pop_j)", "log(tt + 1)")
  fml <- stats::as.formula(paste("share ~",
           paste(c(base_terms, if (n_ch > 1L) "factor(cohort)"), collapse = " + "),
           "+ offset(.off)"))
  fit <- tryCatch(stats::glm(fml, family = stats::quasipoisson(), data = df),
                  error = function(e) {
                    warning(sprintf("[%s] the cohort deterrence fit failed: %s",
                                    label, conditionMessage(e)), call. = FALSE); NULL })
  if (is.null(fit)) return(NULL)
  cf     <- stats::coef(fit)
  b_dist <- unname(cf[["log(tt + 1)"]])
  b_pop_used <- if (is.null(b_pop)) unname(cf[["log(pop_j)"]]) else as.numeric(b_pop)[1]
  if (!is.finite(b_dist) || b_dist >= 0 || !is.finite(b_pop_used)) {
    warning(sprintf("[%s] implausible fit (deterrence %.3f, mass %.3f); kernel not built.",
                    label, b_dist, b_pop_used), call. = FALSE)
    return(NULL)
  }
  # PER-COHORT refits. The between-cohort spread -- not the regression SE, which is a
  # within-cohort quantity over non-independent destinations -- is the honest uncertainty.
  per_fml <- stats::as.formula(paste("share ~", paste(base_terms, collapse = " + "),
                                     "+ offset(.off)"))
  chs <- sort(unique(df$cohort))
  per <- vapply(chs, function(ch) {
    f1 <- tryCatch(stats::glm(per_fml, family = stats::quasipoisson(),
                              data = df[df$cohort == ch, , drop = FALSE]),
                   error = function(e) NULL)
    if (is.null(f1)) NA_real_ else unname(stats::coef(f1)[["log(tt + 1)"]])
  }, numeric(1))
  names(per) <- chs
  message(sprintf(paste0("[%s] Cohort deterrence fitted: log(d+1) = %.3f (within-cohort SE %.3f), ",
                         "mass %.3f (%s); n = %d over %d cohort(s), pseudo-R2 %.3f."),
                  label, b_dist,
                  tryCatch(summary(fit)$coefficients["log(tt + 1)", "Std. Error"],
                           error = function(e) NA_real_),
                  b_pop_used, if (is.null(b_pop)) "fitted" else "fixed at M4's",
                  nrow(df), n_ch, 1 - fit$deviance / fit$null.deviance))
  message(sprintf("[%s] Between-cohort spread of the deterrence (the honest uncertainty): %s.",
                  label, paste(sprintf("%s %.2f", names(per), per), collapse = ", ")))
  list(b_pop = b_pop_used, b_dist = b_dist, per_cohort = per, n = nrow(df),
       n_cohorts = n_ch, mass_fixed = !is.null(b_pop),
       pseudo_r2 = 1 - fit$deviance / fit$null.deviance)
}

#' Build M4c: gravity with M4's destination-mass exponent and the cohort-fitted deterrence.
#'
#' w[o, j] proportional to pop_j^b_pop * (d_oj + 1)^b_dist, row-normalised. The origin-size
#' term and the intercept are origin-constant and cancel exactly under row-normalisation, so
#' they are not carried. Unroutable pairs get zero weight, exactly as build_M4 does (after the
#' OSRM gap fill there are none).
#'
#' @return Row-stochastic 519x519 matrix carrying its coefficients as an attribute.
build_M4c <- function(pop_vec, osrm_mat, zones_all, b_pop, b_dist, label = "M4c") {
  stopifnot(is.finite(b_pop), is.finite(b_dist), b_dist < 0)
  message(sprintf("[%s] Building cohort-calibrated gravity (mass %.3f, deterrence %.3f)...",
                  label, b_pop, b_dist))
  n   <- length(zones_all)
  pop <- as.numeric(pop_vec[zones_all]); pop[is.na(pop)] <- 0
  D   <- .align_dist(osrm_mat, zones_all)
  W   <- make_zero_matrix(zones_all)
  lmass <- b_pop * log(pmax(pop, 1))
  for (i in seq_len(n)) {
    if (pop[i] <= 0) next
    d  <- as.numeric(D[i, ])
    ok <- which(seq_len(n) != i & pop > 0 & is.finite(d) & d > 0)
    if (!length(ok)) next
    eta <- lmass[ok] + b_dist * log(d[ok] + 1)
    # Subtract the row max before exponentiating: a row-constant shift, removed exactly by
    # row-normalisation, that keeps exp() away from underflow on the long-range tail.
    W[i, ok] <- exp(eta - max(eta))
  }
  diag(W) <- 0
  W <- make_row_stochastic(W)
  assert_mobility_matrix(W, zones_all, label)
  message(sprintf("[%s] Built: sparsity=%.1f%%", label, 100 * mean(W == 0)))
  attr(W, "coefficients") <- c(`log(pop_j)` = b_pop, `log(dist_ij + 1)` = b_dist)
  W
}

# ---------------------------------------------------------------------------
# M5: Radiation model (Simini et al. 2012, Nature 484:96-100)
# ---------------------------------------------------------------------------

#' Build M5: radiation model of inter-zone mobility.
#'
#' For each (i,j) pair, s_ij is the sum of populations of zones k strictly
#' closer to i than j (by OSRM time), excluding N_i and N_j.  The radiation
#' flux T_ij ∝ N_i * N_j / ((N_i + s_ij) * (N_i + N_j + s_ij)).
#'
#' Reference: Simini F, González MC, Maritan A, Barabási AL. (2012).
#'   A universal model for mobility and migration patterns. Nature, 484(7392):96-100.
#'
#' @param pop_vec   Named population vector.
#' @param osrm_mat  OSRM travel-time matrix (named).
#' @param zones_all Canonical zone vector.
#' @return Row-stochastic 519×519 matrix.
build_M5 <- function(pop_vec, osrm_mat, zones_all) {
  message("[M5] Building radiation model (Simini et al. 2012)...")

  n <- length(zones_all)

  # Align population
  pop <- pop_vec[zones_all]
  pop[is.na(pop)] <- 0

  # Build full zones_all × zones_all distance matrix
  osrm_zones  <- rownames(osrm_mat)
  shared_osrm <- intersect(zones_all, osrm_zones)
  dist_mat    <- matrix(Inf, nrow = n, ncol = n, dimnames = list(zones_all, zones_all))
  dist_mat[shared_osrm, shared_osrm] <- osrm_mat[shared_osrm, shared_osrm]
  diag(dist_mat) <- 0
  # NA → Inf (unroutable = infinitely far)
  dist_mat[is.na(dist_mat)] <- Inf

  # Radiation model flux matrix (unnormalised)
  message("[M5] Computing radiation fluxes (this may take a moment for 519 zones)...")
  T_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(zones_all, zones_all))

  for (i in seq_len(n)) {
    Ni <- pop[i]
    if (Ni == 0) next
    d_i <- dist_mat[i, ]   # distances from zone i to all others

    for (j in seq_len(n)) {
      if (i == j) next
      Nj <- pop[j]
      if (Nj == 0) next
      d_ij <- d_i[j]
      # OSRM-unreachable destinations (d_ij = Inf) must get ZERO mass. Otherwise the
      # intervening-opportunity mask d_i < Inf admits every finite zone, so s_ij ~ total
      # population and the pair still receives a small positive flux (then amplified by
      # row-normalisation). Gravity M4 already zeroes these; the radiation kernel must too.
      if (!is.finite(d_ij)) next

      # s_ij: sum of pops of zones k with d(i,k) < d(i,j), excluding i and j
      s_ij <- sum(pop[d_i < d_ij & seq_len(n) != i & seq_len(n) != j])

      # Radiation average-flux formula (Simini et al. 2012, Eq. 2; T_i via row-normalisation)
      T_mat[i, j] <- (Ni * Nj) / ((Ni + s_ij) * (Ni + Nj + s_ij))
    }
  }

  diag(T_mat) <- 0
  W_norm <- make_row_stochastic(T_mat)
  assert_mobility_matrix(W_norm, zones_all, "M5")
  message(sprintf("[M5] Built: sparsity=%.1f%%", 100 * mean(W_norm == 0)))
  W_norm
}

# ---------------------------------------------------------------------------
# M6a / M6b: OSRM travel-time decay
# ---------------------------------------------------------------------------

#' Build M6a or M6b: distance-decay mobility matrix from OSRM travel times.
#'
#' M6a (exponential): w_ij = exp(-d_ij / kappa); kappa = 120 min by default.
#' M6b (power law):   w_ij = d_ij^(-gamma); gamma = 1 by default.
#' NA OSRM values are treated as infinite distance (weight ≈ 0).
#'
#' @param osrm_mat  Named OSRM travel-time matrix (minutes).
#' @param zones_all Canonical zone vector.
#' @param variant   "exp" for M6a or "power" for M6b.
#' @param kappa     Decay length scale (minutes); used only for "exp".
#' @param gamma     Power-law exponent; used only for "power".
#' @return Row-stochastic 519×519 matrix.
build_M6 <- function(osrm_mat, zones_all,
                     variant = "exp",
                     kappa   = 120,
                     gamma   = 1) {
  stopifnot(variant %in% c("exp", "power"))
  label <- if (variant == "exp") "M6a" else "M6b"
  message(sprintf("[%s] Building travel-time decay matrix (variant=%s, kappa=%g, gamma=%g)...",
                  label, variant, kappa, gamma))

  n <- length(zones_all)

  # Build full zones_all × zones_all distance matrix; NA (unroutable / OSRM-missing)
  # pairs are zeroed after weighting (see na_mask below), not imputed as "very far".
  osrm_zones  <- rownames(osrm_mat)
  shared_osrm <- intersect(zones_all, osrm_zones)
  dist_mat    <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(zones_all, zones_all))
  dist_mat[shared_osrm, shared_osrm] <- osrm_mat[shared_osrm, shared_osrm]
  diag(dist_mat) <- 0

  # Unroutable pairs (NA distance) get ZERO weight rather than an imputed "very
  # far" distance. Imputing max_dist+1 and then row-normalising turns a zone that
  # is ENTIRELY absent from OSRM (whole row NA → all off-diagonals equal) into a
  # spurious UNIFORM outflow over all 519 zones (the tiny constant cancels under
  # normalisation), so an OSRM-missing infected origin would radiate invasion
  # force everywhere. Zeroing instead leaves such a zone with no outflow (an
  # all-zero row make_row_stochastic keeps at zero), matching M4/M5.
  na_mask <- is.na(dist_mat)

  # Compute weights
  if (variant == "exp") {
    # Exponential decay: w_ij = exp(-d_ij / kappa)
    W_raw <- exp(-dist_mat / kappa)
  } else {
    # Power-law decay: w_ij = d_ij^(-gamma). Only the SELF pair is a true 0-minute distance;
    # distinct centroids are always > 0 min, so mask just the diagonal (masking every == 0
    # would wrongly zero a genuine 0-minute off-diagonal pair, which the power law implies is
    # the MAX-weight pair). Any residual non-finite (an Inf from a hypothetical 0 off-diagonal)
    # is set to 0 defensively; the diagonal is zeroed below regardless.
    dist_nozero <- dist_mat
    diag(dist_nozero) <- NA_real_
    W_raw <- dist_nozero^(-gamma)
    # CAP the +Inf at the largest finite weight rather than zeroing it. A 0-minute off-diagonal
    # pair is, as the comment above says, the MAX-weight pair under a power law — sending it to
    # 0 inverts the kernel for exactly that pair and collapses what should be a dominant edge.
    # (build_M4 already caps its +Inf this way.) Unreachable on the current OSRM matrices, which
    # have no exact off-diagonal zeros; correctness here is about the next matrix, not this one.
    .fin <- W_raw[is.finite(W_raw)]
    if (any(is.infinite(W_raw)) && length(.fin) && max(.fin) > 0) {
      warning(sprintf("[%s] %d zero-distance off-diagonal pair(s) capped at the largest finite weight.",
                      label, sum(is.infinite(W_raw))), call. = FALSE)
      W_raw[is.infinite(W_raw)] <- max(.fin)
    }
    W_raw[!is.finite(W_raw)] <- 0   # NA (diagonal, unroutable) contributes no outflow
  }

  W_raw[na_mask] <- 0     # unroutable / OSRM-missing pairs contribute no outflow
  diag(W_raw) <- 0

  W_norm <- make_row_stochastic(W_raw)
  assert_mobility_matrix(W_norm, zones_all, label)
  message(sprintf("[%s] Built: sparsity=%.1f%%", label, 100 * mean(W_norm == 0)))
  W_norm
}

# ---------------------------------------------------------------------------
# M7: IDP-augmented hybrid
# ---------------------------------------------------------------------------

#' Build M7: relocation flows augmented by IDP displacement flows.
#'
#' Augmented flow: F_aug[i,j] = M3_raw[i,j] + theta * IDP_per_month[i,j], row-normalised.
#'
#' UNITS. M3_raw is ONE MONTH of Flowminder estimated relocations. The IOM IDP table is a
#' CUMULATIVE sum over every week it covers (~10 years), so adding it directly made the IDP
#' term dominate by roughly its number of months — the displacement history outweighed the
#' current month's relocations (e.g. Goma's top destination flipped). It is therefore divided
#' by the window length in months (attr(idp, "n_months")) before weighting, putting both
#' sources on a per-month basis so theta is a dimensionless relative weight (1 = an IDP
#' displacement counts the same as a relocation).
#'
#' @param M3_raw    Raw (unnormalised) 519×519 M3 flow matrix (one month of relocations).
#' @param pop_vec   Named population vector (unused; kept for a uniform builder signature).
#' @param osrm_mat  OSRM travel-time matrix (unused; kept for a uniform builder signature).
#' @param zones_all Canonical zone vector.
#' @param aliases   Alias lookup.
#' @param theta     IDP weight relative to relocations, per month (default 1.0).
#' @return Row-stochastic 519×519 matrix.
build_M7 <- function(M3_raw, pop_vec, osrm_mat, zones_all, aliases, theta = 1.0) {
  message(sprintf("[M7] Building IDP-augmented hybrid (theta=%.2f)...", theta))

  idp_raw <- load_idp_static(aliases, zones_all)   # origins × destinations

  # Build IDP 519×519 matrix
  W_idp <- make_zero_matrix(zones_all)

  idp_origins <- rownames(idp_raw)
  idp_dests   <- colnames(idp_raw)
  shared_orig <- intersect(idp_origins, zones_all)
  shared_dest <- intersect(idp_dests,   zones_all)

  for (oz in shared_orig) {
    for (dz in shared_dest) {
      W_idp[oz, dz] <- idp_raw[oz, dz]
    }
  }
  diag(W_idp) <- 0

  # Put the CUMULATIVE IDP counts on the same MONTHLY basis as the relocation flows.
  n_months <- attr(idp_raw, "n_months")
  if (!is.numeric(n_months) || !is.finite(n_months) || n_months <= 0) {
    warning("[M7] IDP window length unknown; treating the cumulative IDP total as one month.",
            call. = FALSE)
    n_months <- 1
  }
  W_idp <- W_idp / n_months
  message(sprintf("[M7] IDP matrix embedded: %d origin zones with IDP flows (%.1f-month window -> per-month rate).",
                  sum(rowSums(W_idp) > 0), n_months))

  # Augment
  F_aug <- M3_raw + theta * W_idp
  diag(F_aug) <- 0
  F_aug[F_aug < 0] <- 0   # safety

  W_norm <- make_row_stochastic(F_aug)
  assert_mobility_matrix(W_norm, zones_all, "M7")
  message(sprintf("[M7] Built: sparsity=%.1f%%", 100 * mean(W_norm == 0)))
  W_norm
}

# ---------------------------------------------------------------------------
# M8: Composite 
# ---------------------------------------------------------------------------

#' Build M8: composite matrix — short trips for epicentre, gravity elsewhere.
#'
#' For epicentre-origin rows: use M1 (Flowminder short trips).
#' For all other origins: use M4 (calibrated gravity model).
#' This is the recommended primary mobility matrix for the BDBV 2026 analysis.
#'
#' @param M1            Row-stochastic 519×519 M1 matrix.
#' @param M4            Row-stochastic 519×519 M4 matrix.
#' @param epicentre_zones Character vector of epicentre zone names.
#' @param zones_all     Canonical zone vector.
#' @return Row-stochastic 519×519 matrix.
build_M8 <- function(M1, M4, epicentre_zones, zones_all, fill = "none",
                     label = "M8") {
  message(sprintf("[%s] Building composite matrix (epicentre=M1, elsewhere=M4%s)...",
                  label, if (identical(fill, "none")) "" else sprintf(", fill='%s'", fill)))
  # Delegates to compose_epicentre so every composite shares one code path. With
  # fill = "none" this is exactly the previous construction: replace the epicentre
  # rows, zero the diagonal, renormalise, assert.
  compose_epicentre(M4, M1, epicentre_zones, zones_all, label, fill = fill)
}

# ---------------------------------------------------------------------------
# Composite helper — epicentre rows from one kernel, elsewhere from another
# ---------------------------------------------------------------------------

#' Replace the epicentre-origin rows of a base kernel with those of another.
#'
#' Generalises the M8 construction: `base` supplies the outflow structure for all
#' non-epicentre origins; `epi_mat` supplies it for the epicentre origins (whose
#' short-trip behaviour during an outbreak is better captured by Flowminder D+31
#' snapshots than by a static gravity/radiation kernel). Both inputs must already
#' be row-stochastic; the result is re-normalised for safety.
#' @param fill  what to do with destinations the SOURCE could not observe:
#'   "none"       (default, and the behaviour of every classic kernel) — the source
#'                row replaces the base row outright, so an unobservable cell
#'                becomes a hard zero;
#'   "origins"    — only the source's OWN origin zones are refilled from the base.
#'                Both sources exclude them by construction: the annex lists no
#'                column for Bunia/Mongbwalu/Rwampara, and the cohort release states
#'                that "movements within and between those zones are not captured".
#'                Asserting zero there removes the dominant local pathway: the
#'                weight Beni->Butembo, Makiso Kisangani->Mangobo and Makiso
#'                Kisangani->Lubunga are all exactly 0 in M13/M14 while the
#'                radiation base gives 0.051, 0.262 and 0.062;
#'   "unmeasured" — additionally refills every zone the source never measured
#'                (the cohort's "insufficient data" exclusions, and any zone below
#'                a ranked table's reporting cut).
#' With q = the base row's own mass on the refilled set U,
#'   W[o, A] = (1 - q) * p[A]   and   W[o, U] = base[o, U],
#' where p is the source profile renormalised over its observed support A. The row
#' sums to (1 - q) + q = 1 exactly. q is capped at `q_max` and the filled block is
#' RESCALED by the capped ratio (capping alone would leave the row short of 1).
compose_epicentre <- function(base, epi_mat, epicentre_zones, zones_all, label,
                              fill = "none",
                              measured       = attr(epi_mat, "measured"),
                              source_origins = attr(epi_mat, "source_origins"),
                              q_max = get0("MOBILITY_FILL_QMAX", ifnotfound = 0.9)) {
  stopifnot(identical(dim(base), c(length(zones_all), length(zones_all))),
            identical(dim(epi_mat), dim(base)))
  fill <- match.arg(fill, c("none", "origins", "unmeasured"))
  stopifnot(is.numeric(q_max), length(q_max) == 1L, is.finite(q_max),
            q_max >= 0, q_max <= 1)
  W <- base
  n_filled <- 0L
  ez_use <- character(0)
  for (ez in epicentre_zones) {
    if (!ez %in% zones_all) {
      warning(sprintf("[%s] Epicentre zone '%s' not in zones_all; skipping.", label, ez))
      next
    }
    # ez_use is appended only where the row is ACTUALLY replaced with source data — see the
    # three `next` branches below, each of which keeps the base row. Appending here meant
    # attr(W, "source_rows") over-claimed rows that n_filled correctly excluded, and
    # 27_mobility_comparison.R panel D scores TV distance over exactly source_rows, so each
    # over-claimed row injected a spurious zero into the violin.
    src <- epi_mat[ez, ]
    src[!is.finite(src)] <- 0          # a non-finite source cell must not spread
    if (identical(fill, "none")) {
      # An all-zero source row would REPLACE a good base row with "no outflow at all" — the
      # assert-zero-where-the-source-could-not-look failure. The fill path guards this ten lines
      # below (sA <= 0 keeps the base row); the unfilled path did not, and split_cohort_rows()
      # deliberately leaves a cohort's rows at zero when it has no measured mass.
      if (!is.finite(sum(src, na.rm = TRUE)) || sum(src, na.rm = TRUE) <= 0) {
        warning(sprintf("[%s] '%s': the source row carries no mass; keeping the base row.",
                        label, ez), call. = FALSE)
        next
      }
      W[ez, ] <- src; ez_use <- c(ez_use, ez); next
    }

    U <- if (identical(fill, "origins")) {
      so <- if (!is.null(source_origins)) source_origins[[ez]] else NULL
      if (is.null(so)) character(0) else setdiff(intersect(so, zones_all), ez)
    } else {
      me <- if (!is.null(measured)) measured[[ez]] else NULL
      if (is.null(me)) character(0)
      else setdiff(zones_all, c(intersect(me, zones_all), ez))
    }
    if (!length(U)) {
      if (is.null(measured[[ez]]) && is.null(source_origins[[ez]]))
        warning(sprintf("[%s] '%s': no observed-support metadata on the source kernel; fill='%s' had no effect.",
                        label, ez, fill), call. = FALSE)
      W[ez, ] <- src; ez_use <- c(ez_use, ez); next
    }
    A  <- setdiff(zones_all, c(U, ez))
    sA <- sum(src[A], na.rm = TRUE)
    if (!is.finite(sA) || sA <= 0) {
      warning(sprintf("[%s] '%s': the source has no mass outside the refilled set; keeping the base row.",
                      label, ez), call. = FALSE)
      next                       # W[ez, ] is already base[ez, ]
    }
    # NORMALISE by the base row's own total. q is used as a PROPORTION (the complement 1 - q_use
    # is handed to the measured set A), but base[ez, ] is not guaranteed to sum to 1: the M17
    # consensus is documented as summing to <= 1 ("convex combo -> row sums <= 1"). With a base
    # row summing to 0.5 and 0.3 of that on U, the true proportion is 0.6 and the raw sum says
    # 0.3 — and the final renormalisation preserves the wrong split rather than repairing it.
    .base_tot <- sum(base[ez, ], na.rm = TRUE)
    q <- if (is.finite(.base_tot) && .base_tot > 0) sum(base[ez, U], na.rm = TRUE) / .base_tot else 0
    if (!is.finite(q) || q < 0) q <- 0
    q_use <- min(q, q_max)
    if (q > q_max)
      warning(sprintf("[%s] '%s': base mass on the refilled set is %.3f, capped at %.2f.",
                      label, ez, q, q_max), call. = FALSE)
    scl <- if (q > 0) q_use / q else 0
    row <- stats::setNames(numeric(length(zones_all)), zones_all)
    row[A] <- (1 - q_use) * src[A] / sA
    # base[ez, U] / .base_tot is the normalised fill shape; scl then rescales it to q_use.
    # .base_tot == 0 is newly reachable (build_M4's fallback now emits genuinely all-zero rows),
    # and 0/0 would put NaN here. make_row_stochastic() sanitises it downstream, but relying on
    # that is fragile — be explicit.
    row[U] <- if (.base_tot > 0) (base[ez, U] / .base_tot) * scl else 0
    W[ez, ] <- row
    # RECORD THE ROW. ez_use is what attr(W, "source_rows") reports and what
    # 27_mobility_comparison.R panel D scores over; omitting it here left source_rows EMPTY for
    # every -fill kernel (including MOBILITY_PRIMARY = M8-fill) while n_filled counted correctly.
    ez_use <- c(ez_use, ez)
    n_filled <- n_filled + 1L
  }
  diag(W) <- 0
  W <- make_row_stochastic(W)
  assert_mobility_matrix(W, zones_all, label)
  message(sprintf("[%s] Built: sparsity=%.1f%%%s", label, 100 * mean(W == 0),
                  if (!identical(fill, "none"))
                    sprintf(" (source-cell fill='%s' applied to %d row(s))", fill, n_filled)
                  else ""))
  # The origin rows whose profile came from the empirical source rather than the base
  # kernel. Every other row equals the base by construction, so any comparison of a
  # composite against its base is only meaningful on these rows (27_mobility_comparison.R).
  attr(W, "source_rows") <- ez_use
  W
}

# ---------------------------------------------------------------------------
# M_cohort: Flowminder cohort presence kernel (source rows for M13/M14 composites)
# ---------------------------------------------------------------------------

#' Build the pooled Flowminder-cohort presence kernel.
#'
#' For each cohort in `cohort_sources` (Ituri, Nord-Kivu, Tshopo), the cohort's
#' subscriber-day presence vector (row-normalised) becomes the outflow row W[o, ]
#' for every origin zone `o` of that cohort. This is the mobility analogue of M1's
#' short-trip epicentre rows, but sourced from the cohort presence data and spanning
#' all three outbreak hubs. Non-origin rows are left at zero — the composites
#' (compose_epicentre) supply those from a gravity / radiation base kernel.
#'
#' The origins of a cohort share one destination profile (the cohort is pooled), so
#' their rows are identical — the same limitation as M1. Presence is a stock (time),
#' not trips: row-normalising sends most mass to home-adjacent zones (documented in
#' data/flowminder_short_trips/COHORT_INGESTION_PLAN.md).
#'
#' @param zones_all      Canonical zone vector (length n).
#' @param cohort_sources Named list: cohort id -> character vector of origin zone names.
#' @param aliases        Alias lookup from load_aliases().
#' @param window         "followup" (default) or "prior".
#' @return Row-stochastic n x n matrix; attribute "cohort_origins" lists the filled
#'         origin zones (canonical). Rows for zones with no cohort data stay zero.
build_M_cohort <- function(zones_all, cohort_sources, aliases, window = "followup",
                           analysis_date = get0("ANALYSIS_DATE", ifnotfound = NA)) {
  message(sprintf("[M_cohort] Building pooled cohort presence kernel (window=%s)...", window))
  W <- make_zero_matrix(zones_all)
  filled <- character(0)
  measured_by_origin <- list()
  origins_by_origin  <- list()

  for (cohort in names(cohort_sources)) {
    origins_cfg <- harmonise_names(cohort_sources[[cohort]], aliases, zones_all)
    prof        <- load_cohort_profile(cohort, window, zones_all, aliases,
                                       analysis_date = analysis_date)
    vals        <- prof$values
    dest_can    <- names(vals)
    # The cohort's OWN origin zones, as flagged by the provider ("movements within
    # and between those zones are not captured here"). The configured origin list
    # is the fallback when the raw table is unavailable, and the union guards
    # against the two disagreeing.
    origins_src <- union(intersect(prof$origins, zones_all),
                         intersect(origins_cfg, zones_all))

    for (ez in origins_cfg) {
      if (!ez %in% zones_all) {
        warning(sprintf("[M_cohort] Origin '%s' (cohort %s) not in zones_all; skipping.",
                        ez, cohort)); next
      }
      W[ez, ] <- 0
      # Names are already canonical and de-duplicated (summed) by the profile loader.
      W[ez, dest_can] <- vals[dest_can]
      W[ez, ez] <- 0   # no self-loop
      filled <- c(filled, ez)
      measured_by_origin[[ez]] <- dest_can
      origins_by_origin[[ez]]  <- origins_src
    }
  }

  W <- make_row_stochastic(W)
  assert_mobility_matrix(W, zones_all, "M_cohort")
  filled <- unique(filled)
  if (!length(filled))
    stop("[M_cohort] No cohort origins resolved to zones_all — check COHORT_SOURCES.")
  attr(W, "cohort_origins") <- filled
  # Observed support + the provider-excluded origin set, per origin row; see
  # compose_epicentre(fill=). Attached after make_row_stochastic().
  attr(W, "measured")       <- measured_by_origin[filled]
  attr(W, "source_origins") <- origins_by_origin[filled]
  message(sprintf("[M_cohort] Filled %d cohort-origin rows: %s",
                  length(filled), paste(filled, collapse = ", ")))
  W
}

# ---------------------------------------------------------------------------
# Base-kernel helpers: consensus mean, relocation-OD coverage fill, origin split
# ---------------------------------------------------------------------------

#' Weighted element-wise mean of base kernels (the M17 consensus base).
#' Shared by build_M17() and the M17 origin-split so both use the identical base.
.consensus_base <- function(bases, weights = NULL) {
  # weights[seq_along(bases)] takes the FIRST k weights positionally. When build_M17() drops a
  # malformed member, the survivors would then inherit weights belonging to other kernels. Inert
  # while every call site passes equal weights, wrong the first time anyone does not.
  if (is.null(weights)) weights <- rep(1, length(bases))
  if (length(weights) != length(bases))
    stop(sprintf("[consensus] %d weight(s) for %d base kernel(s); the caller must subset weights alongside bases.",
                 length(weights), length(bases)), call. = FALSE)
  w <- weights / sum(weights)
  Reduce(`+`, Map(function(M, wi) wi * M, bases, w))
}

#' The relocation OD (M3) with its COVERAGE GAPS filled, for use as a base kernel.
#'
#' M3 writes 0 wherever the Flowminder table carries no information: an origin absent from the
#' table (or whose every cell is suppressed) gets an all-zero row, and a destination absent from
#' the table gets a zero column. Used as a base, that asserts "no movement" where the data are
#' silent — under the former M16 the invaded zones Mangala and Vuhovi received no incoming weight
#' at all, so their invasion was impossible. This fills exactly those gaps and nothing else, from
#' `fill_base` (radiation, M5, at the call sites):
#'   * an origin with no observed outflow takes fill_base's row;
#'   * an observed origin keeps its measured shape over the destinations the table covers and takes
#'     fill_base's values over the destinations it does not, blended exactly as
#'     compose_epicentre()'s "unmeasured" fill: W[o, A] = (1 - q) p[A], W[o, U] = fill_base[o, U],
#'     with q = fill_base's mass on U, capped at MOBILITY_FILL_QMAX.
#' Suppressed cells INSIDE the table stay 0: they are known to be small (< 15 relocations).
#' Internal only: a base for M16/M17, not a grid kernel.
#'
#' @param M3 build_M3() output (carries "raw" and "od_dests").
#' @param fill_base row-stochastic 519x519 kernel supplying the missing rows and cells.
#' @return row-stochastic 519x519 matrix.
cover_relocation_od <- function(M3, fill_base, zones_all, label = "M3-cov") {
  od_dests <- attr(M3, "od_dests"); raw <- attr(M3, "raw")
  if (is.null(od_dests) || is.null(raw))
    stop(sprintf("[%s] M3 lacks its od_dests / raw attributes (build_M3).", label), call. = FALSE)
  stopifnot(identical(dim(fill_base), dim(M3)))
  observed <- zones_all[rowSums(raw[zones_all, zones_all, drop = FALSE]) > 0]
  src <- M3
  attr(src, "measured") <- stats::setNames(rep(list(od_dests), length(observed)), observed)
  message(sprintf(paste0("[%s] Filling the relocation OD's coverage gaps: %d origin(s) without ",
                         "observed outflow take the base row; %d destination(s) outside the table ",
                         "are filled on the %d observed rows."),
                  label, length(zones_all) - length(observed),
                  length(setdiff(zones_all, od_dests)), length(observed)))
  # Keep the fill-cap warning audible: q (the base mass landing on the unobserved set) runs
  # close to MOBILITY_FILL_QMAX on the widest rows, and a silenced cap would be invisible.
  # Only the per-row "no observed-support metadata" notice is redundant here (every row in
  # `observed` carries it by construction), so warnings are filtered, not suppressed.
  withCallingHandlers(
    compose_epicentre(fill_base, src, observed, zones_all, label, fill = "unmeasured"),
    warning = function(w) {
      if (grepl("no observed-support metadata", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
}

#' Disaggregate each POOLED cohort presence profile into one row per origin.
#'
#' A Flowminder cohort is defined by presence in ANY of its origin zones, so the release gives one
#' destination profile per cohort, and build_M_cohort() copies it to every origin (Beni, Butembo and
#' Katwa receive identical rows). This builds per-origin rows consistent with that measurement AND
#' with each origin's own geography. Over the cohort's measured destinations A (excluding every
#' origin of the cohort):
#'   r_o(j) = B[o, j] c_j * t / sum_k B[o, k] c_k       (each origin keeps the base kernel's shape)
#'   sum_o w_o r_o(j) = p(j)                            (the pooled profile is reproduced exactly)
#' where B is the base kernel, p the pooled profile normalised over A, t the mass fitted (below),
#' and w_o the origin's share of the cohort, taken as its POPULATION share — the only per-origin
#' weight available, since the release does not report members by origin. The column scalings c_j
#' are found by iterative proportional fitting (each step multiplies c_j by p(j) / current mix(j)),
#' which converges when every destination with p(j) > 0 is reachable from some origin under B.
#' Destinations reachable from NO origin (all B[o, j] = 0) carry no information to split, so they
#' keep the pooled share p(j) in every row and the fit distributes the remaining mass t. A cohort
#' whose split cannot be fitted (an origin with no population or no base mass on A, or no
#' convergence) keeps the POOLED profile for all its origins, with a warning.
#'
#' The result is a SOURCE kernel carrying M_cohort's measured / source_origins attributes, to be
#' passed to compose_epicentre(), which then fills each origin's unmeasured destinations (including
#' the other origins of its cohort) from the base. After that fill the population-weighted mean of
#' the final rows equals p only up to each origin's filled share: the fit reproduces the pooled
#' SHAPE over the measured destinations, which is what the release measures.
#'
#' @return 519x519 source kernel (rows only for the cohort origins).
split_cohort_rows <- function(M_cohort, base, cohort_sources, pop_vec, zones_all, aliases,
                              tol = 1e-10, max_iter = 10000L, label = "cohort-split",
                              origins = NULL) {
  measured    <- attr(M_cohort, "measured")
  # Works for ANY pooled source, not only the cohort tables: the Flowminder short-trip annex
  # is pooled over the three epicentre zones in exactly the same way (M1 gives Bunia,
  # Mongbwalu and Rwampara a bit-identical destination profile), and build_M1 attaches
  # "measured"/"source_origins" but no "cohort_origins". Fall back to the measured rows.
  origins_all <- origins %||% attr(M_cohort, "cohort_origins") %||% names(measured)
  if (is.null(measured) || is.null(origins_all))
    stop(sprintf("[%s] source kernel lacks its measured / origin attributes (build_M_cohort / build_M1).",
                 label), call. = FALSE)
  stopifnot(identical(dim(base), dim(M_cohort)))
  S <- make_zero_matrix(zones_all)
  for (cohort in names(cohort_sources)) {
    O <- intersect(harmonise_names(cohort_sources[[cohort]], aliases, zones_all), origins_all)
    if (!length(O)) next
    A <- setdiff(intersect(measured[[O[1]]], zones_all), O)
    p <- M_cohort[O[1], A]
    p[!is.finite(p)] <- 0
    if (!length(A) || sum(p) <= 0) {
      warning(sprintf("[%s] cohort '%s' has no measured mass; its rows stay zero.", label, cohort),
              call. = FALSE)
      next
    }
    p <- p / sum(p)
    keep_pooled <- function(why) {
      warning(sprintf("[%s] cohort '%s': %s; its origins keep the POOLED profile.",
                      label, cohort, why), call. = FALSE)
      for (o in O) S[o, A] <<- p
    }
    if (length(O) == 1L) { S[O, A] <- p; next }
    pw <- pop_vec[O]
    if (any(!is.finite(pw)) || any(pw <= 0)) { keep_pooled("an origin has no population"); next }
    w <- as.numeric(pw / sum(pw))
    B <- base[O, A, drop = FALSE]
    B[!is.finite(B)] <- 0
    act   <- A[p > 0 & colSums(B) > 0]
    fixed <- A[p > 0 & colSums(B) <= 0]
    t_fit <- 1 - sum(p[fixed])
    R <- matrix(0, length(O), length(A), dimnames = list(O, A))
    if (length(fixed)) R[, fixed] <- matrix(p[fixed], length(O), length(fixed), byrow = TRUE)
    it <- 0L; err <- 0
    if (length(act)) {
      Ba <- B[, act, drop = FALSE]
      if (any(rowSums(Ba) <= 0)) {
        keep_pooled("an origin has no base mass on the measured destinations"); next
      }
      target <- p[act]
      cj <- rep(1, length(act))
      repeat {
        Ra  <- sweep(Ba, 2, cj, `*`)
        Ra  <- t_fit * Ra / rowSums(Ra)
        mix <- colSums(w * Ra)                 # w recycles down the rows: row o times w_o
        err <- max(abs(mix - target))
        if (err < tol || it >= max_iter) break
        cj <- cj * target / mix
        it <- it + 1L
      }
      if (!(err < tol)) {
        keep_pooled(sprintf("the split did not converge (max error %.2e after %d iterations)",
                            err, it))
        next
      }
      R[, act] <- Ra
    }
    S[O, A] <- R
    dep <- max(vapply(seq_along(O), function(k) 0.5 * sum(abs(R[k, ] - p)), numeric(1)))
    message(sprintf(paste0("[%s] cohort '%s': %d origins split over %d destinations in %d ",
                           "iteration(s) (max |mix - pooled| %.1e; largest origin departure from ",
                           "the pooled profile, TVD %.3f)."),
                    label, cohort, length(O), length(A), it, err, dep))
  }
  attr(S, "measured")       <- measured
  attr(S, "source_origins") <- attr(M_cohort, "source_origins")
  attr(S, "cohort_origins") <- origins_all
  S
}

#' Origin-split of the pooled SHORT-TRIP annex profile (the M1 epicentre rows).
#'
#' The annex reports ONE ranked destination list for a cohort pooled over Bunia, Mongbwalu and
#' Rwampara, so M1 hands all three origins the identical profile — the same pooling artefact the
#' cohort split exists to correct, on the rows that drive the epicentre's onward spread. This
#' disaggregates them against the composite's own base kernel (population-weighted, reproducing
#' the pooled profile on average), exactly as split_cohort_rows does for the cohort tables.
#'
#' @param M1 build_M1() output (carries "measured" and "source_origins").
#' @param base the composite's background kernel (M4 gravity for M8, M5 radiation for M10).
#' @param epi_origins the epicentre origin zones (canonical).
#' @return Source kernel with split rows, ready for compose_epicentre().
split_shorttrip_rows <- function(M1, base, epi_origins, pop_vec, zones_all, aliases,
                                 label = "shorttrip-split") {
  split_cohort_rows(M1, base, list(shorttrip = epi_origins), pop_vec, zones_all, aliases,
                    label = label, origins = epi_origins)
}

# ---------------------------------------------------------------------------
# M9: Multi-kernel mobility ensemble (consensus)
# ---------------------------------------------------------------------------

#' Build M9: an ensemble ("consensus") mobility kernel.
#'
#' The three structural hypotheses for non-epicentre flow — calibrated gravity
#' (M4), parameter-free radiation (M5) and travel-time exponential decay (M6a) —
#' each capture the mobility field imperfectly and are individually uncertain.
#' Averaging row-stochastic kernels element-wise yields a consensus kernel that
#' is itself row-stochastic (a convex combination of stochastic rows) and is
#' typically better calibrated than any single member — the mobility analogue of
#' the multi-model forecast ensembles that outperform their members (Cramer et
#' al. 2022 PNAS 119:e2113561119; Reich et al. 2019 PNAS 116:3146). Epicentre
#' rows use Flowminder short trips (M1), exactly as in M8.
#'
#' @param M1,M4,M5,M6a Row-stochastic 519×519 kernels.
#' @param weights   Named/positional weights for c(M4, M5, M6a); default equal.
#' @return Row-stochastic 519×519 matrix.
build_M9 <- function(M1, M4, M5, M6a, epicentre_zones, zones_all,
                     weights = c(1, 1, 1)) {
  message("[M9] Building multi-kernel mobility ensemble (mean of M4/M5/M6a; epicentre=M1)...")
  w <- weights / sum(weights)
  base <- w[1] * M4 + w[2] * M5 + w[3] * M6a     # convex combo → row-stochastic
  compose_epicentre(base, M1, epicentre_zones, zones_all, "M9")
}

# ---------------------------------------------------------------------------
# M10: Radiation composite (epicentre short trips + radiation elsewhere)
# ---------------------------------------------------------------------------

#' Build M10: composite with radiation (M5) elsewhere instead of gravity (M4).
#'
#' A parameter-free alternative to M8: the radiation model (Simini et al. 2012)
#' needs no fitted distance exponent, so M10 tests whether the M8 result is
#' robust to the choice of the non-epicentre kernel. Epicentre rows use M1.
#'
#' @return Row-stochastic 519×519 matrix.
build_M10 <- function(M1, M5, epicentre_zones, zones_all, fill = "none",
                      label = "M10") {
  message(sprintf("[%s] Building radiation composite (epicentre=M1, elsewhere=M5%s)...",
                  label, if (identical(fill, "none")) "" else sprintf(", fill='%s'", fill)))
  compose_epicentre(M5, M1, epicentre_zones, zones_all, label, fill = fill)
}

# ---------------------------------------------------------------------------
# M17: Grand all-kernel consensus ensemble
# ---------------------------------------------------------------------------

#' Build M17: the grand all-kernel consensus ensemble.
#'
#' Averages the distinct structural base hypotheses for non-source flow element-wise — at the
#' call sites gravity calibrated on national relocations (M4), gravity calibrated on
#' outbreak-period cohort presence (M4c) and parameter-free radiation (M5). The relocation OD
#' itself (M3), the travel-time decay (M6a) and the symmetrised static kernel (M15) are
#' excluded; see the call site for why M3 is not a member (coverage-driven per-row weights and
#' censored zeros). The mean of row-stochastic kernels is itself a convex combination
#' (row-stochastic up to zero-outflow rows), extending the M9 three-kernel ensemble
#' to the full set of mobility data sources; multi-model consensus kernels are
#' typically better calibrated than any single member (Cramer et al. 2022 PNAS
#' 119:e2113561119; Reich et al. 2019 PNAS 116:3146). The empirical SOURCE rows are
#' then overlaid from the best-available data — the pooled Flowminder cohort presence
#' (Ituri/NK/Tshopo) when available, else the short-trip epicentre kernel (M1) — via
#' compose_epicentre(), which re-normalises the whole matrix so a consensus row that
#' drew mass from a differing number of members becomes a proper distribution.
#' Composite kernels (M8/M9/M10/M13/M13c/M14) are deliberately EXCLUDED from the average to
#' avoid double-counting their constituent bases.
#'
#' @param bases      Named list of >=2 row-stochastic 519x519 base kernels to average.
#' @param epi_mat    Row-stochastic 519x519 kernel supplying the empirical source rows.
#' @param epi_origins Character vector of source-zone origins to overlay from epi_mat.
#' @param zones_all  Canonical zone vector.
#' @param label      Matrix id for messages/assertions ("M17" or "M17-dist").
#' @param weights    Optional positive weights (default equal) over `bases`.
#' @return Row-stochastic 519x519 matrix.
build_M17 <- function(bases, epi_mat, epi_origins, zones_all,
                      label = "M17", weights = NULL, fill = "none") {
  n <- length(zones_all)
  keep <- vapply(bases, function(M)
    is.matrix(M) && all(dim(M) == c(n, n)) && identical(rownames(M), zones_all),
    logical(1))
  if (any(!keep))
    # Never drop a member silently: the consensus is defined by WHICH kernels it averages,
    # and a dropped member changes the kernel without changing its name.
    warning(sprintf("[%s] ignoring %d base kernel(s) that are not %d x %d matrices on zones_all: %s",
                    label, sum(!keep), n, n,
                    paste(names(bases)[!keep] %||% which(!keep), collapse = ", ")),
            call. = FALSE)
  # Subset the WEIGHTS alongside the bases, before anything is dropped — otherwise the
  # survivors silently inherit weights that belonged to the dropped members.
  if (!is.null(weights)) {
    if (length(weights) != length(keep))
      stop(sprintf("[%s] %d weight(s) supplied for %d base kernel(s).",
                   label, length(weights), length(keep)), call. = FALSE)
    weights <- weights[keep]
  }
  bases <- bases[keep]
  if (length(bases) < 2L)
    stop(sprintf("[%s] need >= 2 valid base kernels to form a consensus; got %d.",
                 label, length(bases)))
  if (is.null(weights)) weights <- rep(1, length(bases))
  w <- weights / sum(weights)
  message(sprintf("[%s] Averaging %d base kernels (%s)...",
                  label, length(bases), paste(names(bases), collapse = "+")))
  base <- .consensus_base(bases, w)   # convex combo -> row sums <= 1
  compose_epicentre(base, epi_mat, epi_origins, zones_all, label, fill = fill)
}

# ---------------------------------------------------------------------------
# Summary statistics helper
# ---------------------------------------------------------------------------

#' Compute summary statistics for a mobility matrix.
#'
#' @param W         Row-stochastic matrix.
#' @param matrix_id Character label (e.g., "M1").
#' @param zones_all Canonical zone vector.
#' @param epicentre_zones Epicentre zone names.
#' @return Single-row tibble with summary columns.
summarise_matrix <- function(W, matrix_id, zones_all, epicentre_zones) {
  rs <- rowSums(W)
  # ZERO-INFLOW ZONES. The OSRM-derived builders zero unroutable PAIRS, which is right in the
  # origin direction (an unroutable origin must not radiate everywhere) but has an unexamined
  # consequence in the destination direction: a zone absent from the routing network gets an
  # all-zero COLUMN, so its import hazard is structurally zero forever and it can never be
  # invaded. Measured 2026-09-21 on the saved kernels: Bokoro and Idjwi are zero-inflow under
  # the whole M4/M8/M10 family (10 of the 24 on-disk grid kernels), while the Flowminder
  # relocation table says 21 zones send Bokoro relocations. The cohort/consensus families
  # (M13/M14/M16/M17) have NO zero-inflow zones, and neither does MOBILITY_PRIMARY = M14-fill —
  # the retired comment here named MOBILITY_PRIMARY as the example, which stopped being true
  # when the primary moved off M8-fill. Count and name them so the exclusion is visible.
  .cs <- colSums(W)
  .zero_in <- names(.cs)[.cs <= 0]
  if (length(.zero_in))
    warning(sprintf("[%s] %d zone(s) have ZERO inflow and can never be invaded under this kernel: %s",
                    matrix_id, length(.zero_in),
                    paste(utils::head(.zero_in, 8), collapse = ", ")), call. = FALSE)

  # Top-5 destinations per source origin. The columns are NAMED FROM THE ORIGINS
  # THEMSELVES: the previous version hard-coded top5_Mongbalu, the pre-2026-07
  # spelling, while the rows are keyed on the canonical Mongbwalu, so that column
  # was NA for every kernel ever written.
  keep_origins <- unique(epicentre_zones[epicentre_zones %in% zones_all])
  top5 <- vapply(keep_origins, function(ez) {
    row_w <- W[ez, ]
    # Filter to POSITIVE weights first. Without this an all-zero row returns the first five
    # zone names in spine order as its "top destinations" — and top_destinations_long() does
    # filter, so the two reports contradicted each other for the same kernel.
    row_w <- row_w[is.finite(row_w) & row_w > 0]
    if (!length(row_w)) return(NA_character_)
    paste(names(sort(row_w, decreasing = TRUE))[seq_len(min(5L, length(row_w)))],
          collapse = "; ")
  }, character(1))

  out <- tibble::tibble(
    matrix_id    = matrix_id,
    n_zones      = nrow(W),
    sparsity_pct = 100 * mean(W == 0),
    n_zero_inflow = length(.zero_in),
    mean_weight  = if (any(W > 0)) mean(W[W > 0]) else NA_real_,
    n_zero_rows  = sum(rs < 1e-9)
  )
  if (length(keep_origins)) {
    extra <- as.list(unname(top5))
    names(extra) <- paste0("top5_", keep_origins)
    out <- dplyr::bind_cols(out, tibble::as_tibble(extra))
  }
  out
}

#' Long-format top-N destinations for every SOURCE origin of a kernel.
#'
#' One row per (matrix_id, origin, rank), so cohort origins are reported too
#' rather than only the three epicentre zones of the wide summary.
top_destinations_long <- function(W, matrix_id, origins, zones_all, n_top = 5L) {
  origins <- unique(origins[origins %in% zones_all])
  if (!length(origins)) return(NULL)
  dplyr::bind_rows(lapply(origins, function(ez) {
    row_w <- sort(W[ez, ], decreasing = TRUE)
    row_w <- row_w[row_w > 0]          # a zero-weight "top destination" is noise
    if (!length(row_w)) return(NULL)
    k     <- min(n_top, length(row_w))
    tibble::tibble(matrix_id = matrix_id, origin = ez, rank = seq_len(k),
                   destination = names(row_w)[seq_len(k)],
                   weight = as.numeric(row_w[seq_len(k)]))
  }))
}

# ---------------------------------------------------------------------------
# Master function
# ---------------------------------------------------------------------------

#' Build all mobility matrix variants (M1, M2a, M2b, M3, M4, M4b, M5, M6a, M6b, M7,
#' M8, M9, M10; Flowminder-cohort composites M13/M14; combined-Flowminder-static M15,
#' cohort+directed-OD M16, all-kernel consensus M17; and the road-distance -dist twins).
#'
#' Saves each matrix to outputs/mobility/mobility_{ID}.rds and writes a
#' summary CSV to outputs/mobility/mobility_summary.csv.
#'
#' @param zones_all     Canonical zone vector (length 519).
#' @param pop_vec       Named population vector.
#' @param osrm_mat      OSRM travel-time matrix (named, minutes).
#' @param analysis_date Reference date for time-evolving matrices.
#' @return Named list of row-stochastic matrices, including M8 (composite gravity,
#'   primary), M9 (short-trip kernel ensemble), M10 (radiation composite), M13/M14
#'   (cohort composites), M15 (combined-Flowminder static), M16 (cohort + relocation OD),
#'   M17 (all-kernel consensus ensemble), and their -dist road-distance twins.
build_all_mobility_matrices <- function(zones_all,
                                        pop_vec,
                                        osrm_mat,
                                        analysis_date = ANALYSIS_DATE,
                                        osrm_dist_mat = NULL) {
  # osrm_dist_mat (optional): the OSRM ROAD-DISTANCE (km) matrix. When supplied, the
  # gravity/radiation/decay/composite kernels are ALSO built keyed on road distance rather
  # than travel time, exposed as M4-dist / M8-dist / M9-dist / M10-dist (a mobility-deterrence
  # sensitivity axis; the epicentre short-trip rows are distance-agnostic and reused).
  message("\n=== Building all mobility matrices ===\n")

  aliases       <- load_aliases()
  epicentre_can <- harmonise_names(EPICENTRE_ZONES, aliases, zones_all)
  # Diagnostic: does the (synced, unversioned) processed data still drop zones the
  # raw tables resolve? See .warn_if_processed_source_stale().
  .st_tag_now <- local({
    e <- which(FLOWMINDER_ST_DATES <= analysis_date)
    if (!length(e)) e <- seq_along(FLOWMINDER_ST_DATES)
    FLOWMINDER_ST_TAGS[e[which.max(FLOWMINDER_ST_DATES[e])]]
  })
  suppressMessages(.warn_if_processed_source_stale(zones_all, aliases, .st_tag_now))

  out_list    <- list()
  summary_rows <- list()
  top_rows     <- list()
  # Every origin whose row can come from an empirical SOURCE kernel, for the
  # long-format top-destination report (the wide summary keeps the epicentre trio).
  report_origins <- unique(c(epicentre_can,
    unlist(lapply(get0("COHORT_SOURCES", ifnotfound = list()), function(z)
      harmonise_names(z, aliases, zones_all)), use.names = FALSE)))
  report_origins <- report_origins[report_origins %in% zones_all]

  # Helper: save and summarise
  save_matrix <- function(W, id) {
    rds_path <- file.path(OUT_MOBILITY, sprintf("mobility_%s.rds", id))
    saveRDS(W, rds_path)
    message(sprintf("[save] %s → %s", id, basename(rds_path)))
    out_list[[id]]     <<- W
    summary_rows[[id]] <<- summarise_matrix(W, id, zones_all, epicentre_can)
    top_rows[[id]]     <<- top_destinations_long(W, id, report_origins, zones_all)
  }

  # ---- M3 first (needed for M1/M2 fallback, M4, M7) ----
  M3 <- build_M3(zones_all, aliases)
  M3_raw  <- attr(M3, "raw")
  M3_cens <- attr(M3, "censor_mask")   # NULL/all-FALSE unless the source flags redactions
  save_matrix(M3, "M3")

  # ---- M5 radiation (built early: it supplies the M3 coverage fill used as a base below) ----
  M5 <- build_M5(pop_vec, osrm_mat, zones_all)
  save_matrix(M5, "M5")

  # ---- M3 with its coverage gaps filled (base for M1/M2 fallback rows and M16/M17; not saved) ----
  # M3 asserts zero flow for every zone the relocation table does not cover; as a BASE that made
  # those zones' invasion impossible (Mangala, Vuhovi under the former M16). Filled from radiation.
  M3c <- cover_relocation_od(M3, M5, zones_all, "M3-cov")

  # ---- M1 ----
  # Non-epicentre rows come from M3c, not raw M3: an origin the relocation table never covered
  # must not be handed to the model as a zone that travels nowhere (the same reasoning that
  # makes the source-cell fill the default).
  M1 <- build_M1(zones_all, epicentre_can, aliases, fallback_M3 = M3c,
                 analysis_date = analysis_date)
  save_matrix(M1, "M1")

  # ---- M2a / M2b (LEGACY, off by default) ----
  # Neither is read by any model spec, figure, baseline or cascade path — verified 2026-09-21,
  # every other reference is the MOBILITY_IDS registry. M2b is additionally the SAME KERNEL as
  # M1 (max |M1 - M2b| = 3.9e-16): build_M1 takes "the latest short-trip snapshot on or before
  # the analysis date" and build_M2(variant = "M2b") takes "most recent snapshot" — the 4e-16
  # gap is build_M2's per-snapshot renormalisation. M2a averages five snapshots and still lands
  # within 0.011 of M1, because the annex ranks the same 142 zones at all five dates.
  if (isTRUE(get0("INCLUDE_LEGACY_KERNELS", ifnotfound = FALSE))) {
    M2a <- build_M2(zones_all, epicentre_can, aliases,
                     analysis_date = analysis_date, variant = "M2a",
                     fallback_M3 = M3c)
    save_matrix(M2a, "M2a")

    M2b <- build_M2(zones_all, epicentre_can, aliases,
                     analysis_date = analysis_date, variant = "M2b",
                     fallback_M3 = M3c)
    save_matrix(M2b, "M2b")
  }

  # ---- M4 gravity (power deterrence) ----
  # od_dests keeps zones the source never covered out of the FIT (they would enter as
  # measured or censored zeros); prediction still spans every zone.
  M3_dests <- attr(M3, "od_dests")
  M4 <- build_M4(M3_raw, pop_vec, osrm_mat, zones_all, censor_mask = M3_cens,
                 od_dests = M3_dests)
  save_matrix(M4, "M4")

  # ---- M4b gravity, exponential deterrence (LEGACY, off by default: no consumer) ----
  if (isTRUE(get0("INCLUDE_LEGACY_KERNELS", ifnotfound = FALSE))) {
    M4b <- build_M4(M3_raw, pop_vec, osrm_mat, zones_all,
                    deterrence = "exp", label = "M4b", censor_mask = M3_cens,
                    od_dests = M3_dests)
    save_matrix(M4b, "M4b")
  }

  # ---- M6a exponential decay ----
  M6a <- build_M6(osrm_mat, zones_all, variant = "exp", kappa = 120)
  save_matrix(M6a, "M6a")

  # ---- M6b power-law decay / M7 IDP-augmented (LEGACY, off by default: no consumer) ----
  # M6a is NOT gated: build_M9() reads it, and M9 has its own flag.
  if (isTRUE(get0("INCLUDE_LEGACY_KERNELS", ifnotfound = FALSE))) {
    M6b <- build_M6(osrm_mat, zones_all, variant = "power", gamma = 1)
    save_matrix(M6b, "M6b")

    M7 <- build_M7(M3_raw, pop_vec, osrm_mat, zones_all, aliases, theta = 1.0)
    save_matrix(M7, "M7")
  }

  # ---- M8 composite (RECOMMENDED) ----
  M8 <- build_M8(M1, M4, epicentre_can, zones_all)
  save_matrix(M8, "M8")

  # ---- M9 multi-kernel ensemble (epicentre=M1; elsewhere=mean of M4/M5/M6a) ----
  # OFF by default (INCLUDE_M9_MODELS): consistently beaten by the simpler M8/M10 composites, so
  # neither built nor fit unless explicitly re-enabled.
  want_m9 <- isTRUE(get0("INCLUDE_M9_MODELS", ifnotfound = FALSE))
  if (want_m9) {
    M9 <- build_M9(M1, M4, M5, M6a, epicentre_can, zones_all)
    save_matrix(M9, "M9")
  }

  # ---- M10 radiation composite (epicentre=M1; elsewhere=M5) ----
  M10 <- build_M10(M1, M5, epicentre_can, zones_all)
  save_matrix(M10, "M10")

  # ---- M15: Flowminder SYMMETRISED OD kernel (S = O + t(O) in canonical space) ----
  # One directed OD table symmetrised; there is no independent inflow table (build_M15).
  # OFF by default (INCLUDE_M15_MODELS): among the weakest/spiky standalone kernels, and no
  # longer a dependency of any default kernel (M16 now uses the directed OD M3; M17 drops it).
  # want_flowstat still governs the M16/M17 family below.
  want_flowstat <- isTRUE(get0("INCLUDE_FLOWSTATIC_MODELS", ifnotfound = TRUE))
  want_m15 <- isTRUE(get0("INCLUDE_M15_MODELS", ifnotfound = FALSE))
  M15 <- NULL
  if (want_m15) {
    M15 <- tryCatch(build_M15(zones_all, aliases),
                    error = function(e) { warning("[M15] build failed: ", conditionMessage(e)); NULL })
    if (!is.null(M15)) save_matrix(M15, "M15")
  }

  # Source-cell fill level (MOBILITY_SOURCE_FILL; "unmeasured" by default) and the origin-split
  # switch, set once here: the split kernels below use the same fill as the -fill twins.
  want_fill  <- isTRUE(get0("INCLUDE_SOURCEFILL_MODELS", ifnotfound = TRUE))
  fill_level <- get0("MOBILITY_SOURCE_FILL", ifnotfound = "unmeasured")
  split_fill <- if (want_fill) fill_level else "none"
  want_split <- isTRUE(get0("INCLUDE_COHORT_SPLIT_MODELS", ifnotfound = TRUE))

  # ---- M_cohort + M13/M14: Flowminder-cohort composites (travel-time base) ----------------
  # M13 (cohort + gravity) and M14 (cohort + radiation) fill the Ituri/NK/Tshopo cohort-origin
  # rows from cohort subscriber-day PRESENCE (build_M_cohort) and take the M4 gravity / M5
  # radiation travel-time kernel elsewhere — the cohort-data analogue of M8/M10 (which use the
  # short-trip M1 for the Ituri epicentre only). Additive: M1/M8/M10 are untouched. Gated by
  # INCLUDE_COHORT_MODELS (default TRUE) and the presence of COHORT_SOURCES; get0() keeps a
  # config-less unit-test build working.
  cohort_sources <- get0("COHORT_SOURCES", ifnotfound = NULL)
  want_cohort <- isTRUE(get0("INCLUDE_COHORT_MODELS", ifnotfound = TRUE)) &&
                 !is.null(cohort_sources) && length(cohort_sources) > 0L
  M_cohort <- NULL; cohort_origins <- character(0)
  # M4c (gravity with the deterrence calibrated on the cohort tables) and its fit. Declared
  # here so they exist when the cohort family is off -- the M17 consensus below filters NULLs.
  M4c <- NULL; cohort_grav <- NULL
  if (want_cohort) {
    M_cohort <- tryCatch(
      # analysis_date MUST be threaded: without it build_M_cohort() falls back to the GLOBAL
      # ANALYSIS_DATE, so a back-dated build_all_mobility_matrices(analysis_date = <earlier>)
      # still picked the newest cohort release — the exact leak the as-of selection fixes.
      build_M_cohort(zones_all, cohort_sources, aliases,
                     window = get0("COHORT_WINDOW", ifnotfound = "followup"),
                     analysis_date = analysis_date),
      error = function(e) { warning("[M_cohort] build failed: ", conditionMessage(e)); NULL })
    if (!is.null(M_cohort)) {
      cohort_origins <- attr(M_cohort, "cohort_origins")
      # PERSIST THE COHORT KERNEL ITSELF. It was only ever consumed to build M13/M14/M16 and
      # never saved, so the richest observed mobility this study holds — Flowminder cohort
      # subscriber-day PRESENCE, measured DURING the outbreak, per origin rather than pooled —
      # could not be read back by anything. Manuscript Figure 2's Flowminder-inflow baseline
      # uses it directly: measured against the pooled short-trip annex (M1) it covers 305
      # destinations from the epicentre instead of 142, while ranking the same zones at the top
      # (Spearman 0.79). It is NOT added to MOBILITY_IDS: it is a base kernel and a baseline
      # source, not a candidate in the model grid.
      save_matrix(M_cohort, "M_cohort")
      # ---- M4c: gravity, deterrence calibrated on COHORT PRESENCE ------------------------
      # Fitted here because it needs M_cohort. Used twice below: as a member of the M17
      # consensus (so cohort information reaches all 519 base rows, which it cannot do as
      # rows -- see the M4c block above) and as the base kernel of the M13c composite.
      # The destination-mass exponent is HELD AT M4's; if M4 exposed no coefficients, fall
      # back to fitting it, with a warning, because three pooled cohorts identify it poorly.
      if (isTRUE(get0("INCLUDE_COHORT_GRAVITY_MODELS", ifnotfound = TRUE))) {
        .m4_cf <- attr(M4, "coefficients")
        .b_pop <- if (!is.null(.m4_cf) && "log(pop_j)" %in% names(.m4_cf))
                    unname(.m4_cf[["log(pop_j)"]]) else NULL
        if (is.null(.b_pop))
          warning(paste0("[M4c] M4 exposed no log(pop_j) coefficient; fitting the mass exponent ",
                         "from three pooled cohorts, which cannot identify it well."),
                  call. = FALSE)
        cohort_grav <- tryCatch(
          fit_cohort_gravity(M_cohort, cohort_sources, pop_vec, osrm_mat, zones_all, aliases,
                             b_pop = .b_pop, label = "M4c"),
          error = function(e) { warning("[M4c] fit failed: ", conditionMessage(e)); NULL })
        if (!is.null(cohort_grav)) {
          M4c <- build_M4c(pop_vec, osrm_mat, zones_all, cohort_grav$b_pop, cohort_grav$b_dist)
          save_matrix(M4c, "M4c")
        }
      }
      save_matrix(compose_epicentre(M4, M_cohort, cohort_origins, zones_all, "M13"), "M13")
      save_matrix(compose_epicentre(M5, M_cohort, cohort_origins, zones_all, "M14"), "M14")
      # ---- M13c: cohort source rows over the cohort-calibrated gravity base ---------------
      # The M13/M14 analogue with M4c as the background kernel, so BOTH the source rows and
      # the deterrence elsewhere come from outbreak-period mobility rather than relocations.
      if (!is.null(M4c))
        save_matrix(compose_epicentre(M4c, M_cohort, cohort_origins, zones_all, "M13c"), "M13c")
      # ---- M16: cohort + Flowminder RELOCATION OD composite (cohort where available, M3 elsewhere) --
      # The analogue of M13/M14 pairing the cohort presence rows with the empirical directed
      # Flowminder relocation flows — M3 with its coverage gaps filled (M3c), since the raw M3
      # asserts zero flow out of and into every zone the relocation table does not cover.
      if (want_flowstat)
        save_matrix(compose_epicentre(M3c, M_cohort, cohort_origins, zones_all, "M16"), "M16")
      # ---- ORIGIN-SPLIT cohort composites (M13/M14/M16-split) ------------------------------------
      # The pooled cohort profile disaggregated per origin against each composite's own base
      # (split_cohort_rows()), then filled like the -fill twins. INCLUDE_COHORT_SPLIT_MODELS.
      if (want_split) {
        .split <- function(base, lab) split_cohort_rows(M_cohort, base, cohort_sources, pop_vec,
                                                        zones_all, aliases, label = lab)
        save_matrix(compose_epicentre(M4, .split(M4, "M13-split"), cohort_origins, zones_all,
                                      "M13-split", fill = split_fill), "M13-split")
        save_matrix(compose_epicentre(M5, .split(M5, "M14-split"), cohort_origins, zones_all,
                                      "M14-split", fill = split_fill), "M14-split")
        if (!is.null(M4c))
          save_matrix(compose_epicentre(M4c, .split(M4c, "M13c-split"), cohort_origins, zones_all,
                                        "M13c-split", fill = split_fill), "M13c-split")
        if (want_flowstat)
          save_matrix(compose_epicentre(M3c, .split(M3c, "M16-split"), cohort_origins, zones_all,
                                        "M16-split", fill = split_fill), "M16-split")
      }
    }
  }

  # ---- M17: grand all-kernel consensus ensemble (travel-time base) ----
  # Mean of the distinct structural bases {M4, M4c, M5}, with the empirical source rows overlaid
  # from the cohort kernel (when available) else the short-trip epicentre M1. The travel-time
  # decay kernel (M6a) and the symmetrised static kernel (M15) are deliberately EXCLUDED, and so
  # are all composites (M8/M9/M10/M13/M13c/M14) -- see below. Gated by INCLUDE_FLOWSTATIC.
  if (want_flowstat) {
    .m17_epi_mat    <- if (!is.null(M_cohort)) M_cohort else M1
    .m17_epi_origin <- if (!is.null(M_cohort)) cohort_origins else epicentre_can
    # M3 IS NOT A MEMBER (changed 2026-09-22). It used to be, and it made the consensus two
    # different things at once:
    #
    #  (a) COVERAGE-DRIVEN WEIGHTS. M3's row is all-zero wherever the relocation table has no
    #      entry for that origin, so the member simply dropped out and the survivors shared the
    #      row. Measured on the previous build: 402 of the 509 base rows were the nominal
    #      three-way mean, 105 were M4+M5 at 1/2 each, 1 was a single member and 1 was empty.
    #      Radiation's weight therefore rose from 1/3 to 1/2 depending on whether Flowminder
    #      happened to cover that origin -- an artefact of table coverage, not of mobility.
    #  (b) CENSORED ZEROS PROPAGATED AS REAL ZEROS. M3 carries only 7,100 of 268,842 off-diagonal
    #      pairs (2.6%), and a suppressed count (<15) is written as a hard zero. Because the pool
    #      is per-row, every pair M3 did not report -- 97.4% of the pairs where M4 and M5 both
    #      give positive weight -- took a factor of exactly 2/3 against a pool of the members
    #      that could actually observe it, with the mass redirected onto M3's 2.6% support. For
    #      an INVASION kernel that is the worst place to take it from: the signal is the small
    #      long-range weight into not-yet-invaded zones, which is exactly the set M3 misses.
    #      It also contradicted M4_CENSORED_FIT, which exists because treating suppressed counts
    #      as observed zeros biases the deterrence.
    #
    # The members are now {M4, M4c, M5}: gravity calibrated on national RELOCATIONS (shallow
    # decay), gravity calibrated on outbreak-period COHORT PRESENCE (steep decay), and
    # parameter-free radiation. Three structural hypotheses, two independent empirical
    # calibration sources, no composites and no censored member -- so every row is an
    # equal-weight mean of the same three kernels, with no coverage-driven drift.
    #
    # WHY NOT SIMPLY ADD THE COMPOSITES (M8/M10/M13/M14) AS MEMBERS. They differ from their
    # bases ONLY on the source rows, and compose_epicentre() overwrites those rows afterwards.
    # mean{M4, M5, M8, M10, M13, M14} equals (M4 + M5)/2 on every base row to 8e-16 -- verified
    # -- so three of those six members contribute nothing but a misleading name. That is what
    # the "composites are deliberately EXCLUDED" rule in build_M17()'s docstring protects
    # against, and it is why cohort data enter here through M4c's DETERRENCE instead.
    #
    # M3 has not left the pipeline: it remains M4's calibration target and the base of M16.
    .b17 <- Filter(Negate(is.null), list(M4 = M4, M4c = M4c, M5 = M5))
    M17 <- tryCatch(
      build_M17(.b17, .m17_epi_mat, .m17_epi_origin, zones_all, "M17"),
      error = function(e) { warning("[M17] build failed: ", conditionMessage(e)); NULL })
    if (!is.null(M17)) save_matrix(M17, "M17")
    if (want_split && !is.null(M_cohort)) {
      M17s <- tryCatch(
        build_M17(.b17, split_cohort_rows(M_cohort, .consensus_base(.b17), cohort_sources, pop_vec,
                                          zones_all, aliases, label = "M17-split"),
                  cohort_origins, zones_all, "M17-split", fill = split_fill),
        error = function(e) { warning("[M17-split] build failed: ", conditionMessage(e)); NULL })
      if (!is.null(M17s)) save_matrix(M17s, "M17-split")
    }
  }

  # ---- SOURCE-CELL FILL variants (M8/M13/M14/M16/M17-fill) --------------------------------
  # Same composites, but the destinations the SOURCE COULD NOT OBSERVE are taken
  # from the base kernel instead of being asserted as zero (compose_epicentre's
  # `fill`; see its docstring for the construction). These are the DEFAULT kernels of the
  # model suite (INCLUDE_UNFILLED_MODELS = FALSE drops the unfilled parents from the grid);
  # the unfilled kernels above are still BUILT, with fill = "none", as a sensitivity arm.
  if (want_fill && !identical(fill_level, "none")) {
    message(sprintf("[mobility] Building source-cell fill variants (fill='%s')...", fill_level))
    save_matrix(build_M8(M1, M4, epicentre_can, zones_all, fill = fill_level,
                         label = "M8-fill"), "M8-fill")
    save_matrix(build_M10(M1, M5, epicentre_can, zones_all, fill = fill_level,
                          label = "M10-fill"), "M10-fill")
    # ---- ORIGIN-SPLIT short-trip composites (M8/M10-split) ---------------------------------
    # The annex pools its cohort over the three epicentre zones, so M1 gives Bunia, Mongbwalu
    # and Rwampara a bit-identical destination profile (verified: the measured-cell shape is
    # equal to 1e-16 in M8-fill/M10-fill; only the fill mass differs, 0.28/0.39/0.38). These
    # split that pooled profile per origin against each composite's own base, the short-trip
    # analogue of M13/M14-split.
    if (want_split) {
      save_matrix(compose_epicentre(M4, split_shorttrip_rows(M1, M4, epicentre_can, pop_vec,
                                                             zones_all, aliases, "M8-split"),
                                    epicentre_can, zones_all, "M8-split",
                                    fill = split_fill), "M8-split")
      save_matrix(compose_epicentre(M5, split_shorttrip_rows(M1, M5, epicentre_can, pop_vec,
                                                             zones_all, aliases, "M10-split"),
                                    epicentre_can, zones_all, "M10-split",
                                    fill = split_fill), "M10-split")
    }
    if (!is.null(M_cohort)) {
      save_matrix(compose_epicentre(M4, M_cohort, cohort_origins, zones_all, "M13-fill",
                                    fill = fill_level), "M13-fill")
      save_matrix(compose_epicentre(M5, M_cohort, cohort_origins, zones_all, "M14-fill",
                                    fill = fill_level), "M14-fill")
      if (!is.null(M4c))
        save_matrix(compose_epicentre(M4c, M_cohort, cohort_origins, zones_all, "M13c-fill",
                                      fill = fill_level), "M13c-fill")
      if (want_flowstat)
        save_matrix(compose_epicentre(M3c, M_cohort, cohort_origins, zones_all, "M16-fill",
                                      fill = fill_level), "M16-fill")
    }
    if (want_flowstat) {
      .f_epi_mat    <- if (!is.null(M_cohort)) M_cohort else M1
      .f_epi_origin <- if (!is.null(M_cohort)) cohort_origins else epicentre_can
      M17f <- tryCatch(
        build_M17(Filter(Negate(is.null), list(M4 = M4, M4c = M4c, M5 = M5)),  # see the M17 note above
                  .f_epi_mat, .f_epi_origin, zones_all, "M17-fill", fill = fill_level),
        error = function(e) { warning("[M17-fill] build failed: ", conditionMessage(e)); NULL })
      if (!is.null(M17f)) save_matrix(M17f, "M17-fill")
    }
  }

  # ---- OSRM ROAD-DISTANCE kernels (generic M4/M8/M9/M10-dist AND cohort M13/M14-dist) -----
  # Two independent families keyed on OSRM road DISTANCE (km) instead of travel TIME (minutes):
  #  * generic -dist (M4/M8/M9/M10-dist): the short-trip-composite deterrence axis, gated by
  #    INCLUDE_OSRM_DIST_MODELS (get0 default TRUE, so an arm without that flag keeps building it).
  #  * cohort -dist (M13/M14-dist): the geographic-distance analogues of M13/M14, built whenever
  #    the cohort kernel exists and a road-distance matrix is supplied — INDEPENDENT of the generic
  #    toggle (so "flowminder + {gravity,radiation} + geographic distance" is available by default).
  # The km gravity/radiation bases (M4_d/M5_d) are built ONCE and shared by both. The gravity GLM
  # RE-FITS on the km covariate (not a rescale); radiation re-orders opportunities by km.
  build_generic_dist  <- !is.null(osrm_dist_mat) &&
                         isTRUE(get0("INCLUDE_OSRM_DIST_MODELS", ifnotfound = TRUE))
  build_cohort_dist   <- !is.null(osrm_dist_mat) && !is.null(M_cohort)
  build_flowstat_dist <- !is.null(osrm_dist_mat) && want_flowstat
  if (build_generic_dist || build_cohort_dist || build_flowstat_dist) {
    message(sprintf("[mobility] Building OSRM road-distance kernels (generic-dist=%s, cohort-dist=%s, flowstat-dist=%s)...",
                    build_generic_dist, build_cohort_dist, build_flowstat_dist))
    M4_d <- build_M4(M3_raw, pop_vec, osrm_dist_mat, zones_all, label = "M4-dist",
                     censor_mask = M3_cens, od_dests = M3_dests)
    M5_d <- build_M5(pop_vec, osrm_dist_mat, zones_all)
    # Road-km twin of the cohort-calibrated gravity. Built inside the cohort -dist block below
    # (it needs M_cohort); declared here so M17-dist can filter it out when it is absent.
    M4c_d <- NULL
    # M6a_d (km travel-time-decay) is only needed by the generic M9-dist ensemble (M17-dist no
    # longer uses it), so build it only when M9 is enabled AND the generic -dist family is wanted.
    M6a_d <- NULL
    if (build_generic_dist && want_m9) {
      # M6a's exponential deterrence length was calibrated as 120 TRAVEL-TIME MINUTES; reusing the
      # numeral 120 as KILOMETRES for M6a_d would impose a physically different scale. Convert the
      # 120-minute length to km via the road network's own median speed (km/min) over shared finite
      # off-diagonal pairs, so M6a_d decays over the same typical trips as M6a. Fall back to 60 km.
      kappa_km <- tryCatch({
        st <- intersect(rownames(osrm_dist_mat), rownames(osrm_mat))
        st <- intersect(st, intersect(colnames(osrm_dist_mat), colnames(osrm_mat)))
        dd <- osrm_dist_mat[st, st]; tt <- osrm_mat[st, st]
        ok <- is.finite(dd) & is.finite(tt) & tt > 0; diag(ok) <- FALSE
        sp <- stats::median(dd[ok] / tt[ok], na.rm = TRUE)   # km per minute
        if (is.finite(sp) && sp > 0) 120 * sp else 60
      }, error = function(e) 60)
      message(sprintf("[mobility] M6a-dist decay length = %.0f km (120 min at the network median %.2f km/min)",
                      kappa_km, kappa_km / 120))
      M6a_d <- build_M6(osrm_dist_mat, zones_all, variant = "exp", kappa = kappa_km)
    }
    if (build_generic_dist) {
      save_matrix(M4_d, "M4-dist")
      save_matrix(build_M8(M1, M4_d, epicentre_can, zones_all), "M8-dist")
      # M9-dist only when M9 is enabled (needs the km decay kernel M6a_d, built under the same guard).
      if (want_m9)
        save_matrix(build_M9(M1, M4_d, M5_d, M6a_d, epicentre_can, zones_all), "M9-dist")
      save_matrix(build_M10(M1, M5_d, epicentre_can, zones_all), "M10-dist")
      if (want_fill && !identical(fill_level, "none")) {
        save_matrix(build_M8(M1, M4_d, epicentre_can, zones_all, fill = fill_level,
                             label = "M8-dist-fill"), "M8-dist-fill")
        save_matrix(build_M10(M1, M5_d, epicentre_can, zones_all, fill = fill_level,
                              label = "M10-dist-fill"), "M10-dist-fill")
      }
      if (want_split) {
        # Road-km twins of the short-trip origin-split composites (split against the km base).
        save_matrix(compose_epicentre(M4_d, split_shorttrip_rows(M1, M4_d, epicentre_can, pop_vec,
                                        zones_all, aliases, "M8-dist-split"),
                                      epicentre_can, zones_all, "M8-dist-split",
                                      fill = split_fill), "M8-dist-split")
        save_matrix(compose_epicentre(M5_d, split_shorttrip_rows(M1, M5_d, epicentre_can, pop_vec,
                                        zones_all, aliases, "M10-dist-split"),
                                      epicentre_can, zones_all, "M10-dist-split",
                                      fill = split_fill), "M10-dist-split")
      }
    }
    if (build_cohort_dist) {
      # M4c-dist: the cohort deterrence RE-FITTED on the km covariate (not a rescale of the
      # travel-time fit), exactly as M4-dist re-fits M4's. Mass exponent from M4-dist.
      if (isTRUE(get0("INCLUDE_COHORT_GRAVITY_MODELS", ifnotfound = TRUE))) {
        .m4d_cf  <- attr(M4_d, "coefficients")
        .b_pop_d <- if (!is.null(.m4d_cf) && "log(pop_j)" %in% names(.m4d_cf))
                      unname(.m4d_cf[["log(pop_j)"]]) else NULL
        .cg_d <- tryCatch(
          fit_cohort_gravity(M_cohort, cohort_sources, pop_vec, osrm_dist_mat, zones_all,
                             aliases, b_pop = .b_pop_d, label = "M4c-dist"),
          error = function(e) { warning("[M4c-dist] fit failed: ", conditionMessage(e)); NULL })
        if (!is.null(.cg_d)) {
          M4c_d <- build_M4c(pop_vec, osrm_dist_mat, zones_all, .cg_d$b_pop, .cg_d$b_dist,
                             "M4c-dist")
          save_matrix(M4c_d, "M4c-dist")
        }
      }
      save_matrix(compose_epicentre(M4_d, M_cohort, cohort_origins, zones_all, "M13-dist"), "M13-dist")
      save_matrix(compose_epicentre(M5_d, M_cohort, cohort_origins, zones_all, "M14-dist"), "M14-dist")
      if (!is.null(M4c_d))
        save_matrix(compose_epicentre(M4c_d, M_cohort, cohort_origins, zones_all, "M13c-dist"),
                    "M13c-dist")
      if (want_fill && !identical(fill_level, "none")) {
        save_matrix(compose_epicentre(M4_d, M_cohort, cohort_origins, zones_all, "M13-dist-fill",
                                      fill = fill_level), "M13-dist-fill")
        save_matrix(compose_epicentre(M5_d, M_cohort, cohort_origins, zones_all, "M14-dist-fill",
                                      fill = fill_level), "M14-dist-fill")
        if (!is.null(M4c_d))
          save_matrix(compose_epicentre(M4c_d, M_cohort, cohort_origins, zones_all,
                                        "M13c-dist-fill", fill = fill_level), "M13c-dist-fill")
      }
      if (want_split) {
        save_matrix(compose_epicentre(M4_d, split_cohort_rows(M_cohort, M4_d, cohort_sources,
                                        pop_vec, zones_all, aliases, label = "M13-dist-split"),
                                      cohort_origins, zones_all, "M13-dist-split",
                                      fill = split_fill), "M13-dist-split")
        save_matrix(compose_epicentre(M5_d, split_cohort_rows(M_cohort, M5_d, cohort_sources,
                                        pop_vec, zones_all, aliases, label = "M14-dist-split"),
                                      cohort_origins, zones_all, "M14-dist-split",
                                      fill = split_fill), "M14-dist-split")
        if (!is.null(M4c_d))
          save_matrix(compose_epicentre(M4c_d, split_cohort_rows(M_cohort, M4c_d, cohort_sources,
                                          pop_vec, zones_all, aliases, label = "M13c-dist-split"),
                                        cohort_origins, zones_all, "M13c-dist-split",
                                        fill = split_fill), "M13c-dist-split")
      }
    }
    if (build_flowstat_dist) {
      # M17-dist: all-kernel consensus on ROAD-KM over the bases {M4-dist, M4c-dist, M5-dist}. Gravity and
      # radiation are re-keyed on km (M4_d/M5_d); the empirical Flowminder relocation OD (M3) carries
      # no distance axis, so it is reused, with its coverage gaps filled from the km radiation kernel.
      # The travel-time decay (M6a) and symmetrised static (M15) are excluded, matching the
      # travel-time M17. Empirical source rows from the cohort kernel else M1.
      # M3 is NOT a member here either, for the reasons given at the travel-time M17 above.
      .b17d <- Filter(Negate(is.null),
                      list(`M4-dist` = M4_d, `M4c-dist` = M4c_d, `M5-dist` = M5_d))
      .m17d_epi_mat    <- if (!is.null(M_cohort)) M_cohort else M1
      .m17d_epi_origin <- if (!is.null(M_cohort)) cohort_origins else epicentre_can
      save_matrix(build_M17(.b17d, .m17d_epi_mat, .m17d_epi_origin, zones_all, "M17-dist"),
                  "M17-dist")
      if (want_fill && !identical(fill_level, "none"))
        save_matrix(build_M17(.b17d, .m17d_epi_mat, .m17d_epi_origin, zones_all, "M17-dist-fill",
                              fill = fill_level), "M17-dist-fill")
      if (want_split && !is.null(M_cohort))
        save_matrix(build_M17(.b17d, split_cohort_rows(M_cohort, .consensus_base(.b17d),
                                cohort_sources, pop_vec, zones_all, aliases,
                                label = "M17-dist-split"),
                              cohort_origins, zones_all, "M17-dist-split", fill = split_fill),
                    "M17-dist-split")
    }
  }

  # ---- Summary CSV ----
  summary_df <- dplyr::bind_rows(summary_rows)
  summary_path <- file.path(OUT_MOBILITY, "mobility_summary.csv")
  readr::write_csv(summary_df, summary_path)
  message(sprintf("[summary] Written to %s", summary_path))

  # ---- Long-format top destinations (all source origins, not just the epicentre) ----
  top_df <- dplyr::bind_rows(top_rows)
  if (nrow(top_df))
    readr::write_csv(top_df, file.path(OUT_MOBILITY, "mobility_top_destinations.csv"))

  # ---- MANIFEST: exactly the kernels THIS run built -------------------------
  # 27_mobility_comparison.R reads every mobility_*.rds in the directory, so a file
  # left behind by an earlier run (a kernel since switched off, or one from another
  # branch) silently enters the similarity matrix as if it were current. The
  # manifest is the run's own record of what it wrote.
  # CONFIG FINGERPRINT. matrix_id + file + built_at cannot distinguish a kernel built under a
  # different FLOWMINDER_OD_FILE / MOBILITY_SOURCE_FILL / COHORT_WINDOW / MOBILITY_FILL_QMAX
  # from a current one, and 30_projection_config.R adopts CASCADE_KERNEL on file EXISTENCE
  # alone — so a stale mobility_*.rds from another branch is indistinguishable from this run's.
  # Recording the governing values lets a consumer check rather than assume.
  .cfg <- list(
    # The ARGUMENT, not the global. analysis_date is the axis that now drives the as-of
    # selection of the short-trip snapshot and the cohort release, so recording the global made
    # a back-dated build and a real-time build hash IDENTICAL — defeating the one case the
    # fingerprint was added for.
    analysis_date        = format(analysis_date),
    flowminder_od_file   = as.character(get0("FLOWMINDER_OD_FILE",   ifnotfound = NA)),
    mobility_source_fill = as.character(get0("MOBILITY_SOURCE_FILL", ifnotfound = NA)),
    mobility_fill_qmax   = as.character(get0("MOBILITY_FILL_QMAX",   ifnotfound = NA)),
    cohort_window        = as.character(get0("COHORT_WINDOW",        ifnotfound = NA)),
    mobility_home_frac   = as.character(get0("MOBILITY_HOME_FRACTION", ifnotfound = NA)))
  .cfg_str <- paste(names(.cfg), unlist(.cfg), sep = "=", collapse = "; ")
  manifest <- tibble::tibble(
    matrix_id = names(out_list),
    file      = sprintf("mobility_%s.rds", names(out_list)),
    built_at  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    config    = .cfg_str,
    config_hash = if (requireNamespace("digest", quietly = TRUE))
      digest::digest(.cfg, algo = "xxhash64") else substr(
        paste0(utils::capture.output(str(.cfg)), collapse = ""), 1, 16))
  readr::write_csv(manifest, file.path(OUT_MOBILITY, "mobility_manifest.csv"))
  message(sprintf("[manifest] %d kernel(s) recorded in mobility_manifest.csv", nrow(manifest)))

  message("\n=== All mobility matrices complete ===")
  message(sprintf("Matrices built: %s", paste(names(out_list), collapse = ", ")))
  message(sprintf("Primary (M8) zero rows: %d / %d",
                  sum(rowSums(M8) < 1e-9), length(zones_all)))

  out_list
}

# ---------------------------------------------------------------------------
# MAIN — Execute only when run standalone (guarded to avoid a costly
# double-build when sourced from run_all.R, which calls the builder itself).
# ---------------------------------------------------------------------------

if (interactive() && !exists("MOBILITY_MATRICES")) {
  message("\n=== 03_mobility_matrices.R: Loading base data and building all matrices ===\n")

  # Load canonical zone list from WorldPop (the population spine defines the zone set — do NOT
  # hard-assert a fixed count, so a changed admin structure / new zones do not break the build).
  pop_raw  <- load_worldpop()
  zones_all <- names(pop_raw)
  stopifnot(length(zones_all) > 0L)
  if (length(zones_all) != 519L)
    message(sprintf("[mobility] %d zones in the WorldPop spine (reference build had 519).", length(zones_all)))

  # Load OSRM matrix
  osrm_mat <- load_osrm()

  # Build all matrices
  MOBILITY_MATRICES <- build_all_mobility_matrices(
    zones_all     = zones_all,
    pop_vec       = pop_raw,
    osrm_mat      = osrm_mat,
    analysis_date = ANALYSIS_DATE
  )

  message("\n=== 03_mobility_matrices.R complete ===\n")
}
