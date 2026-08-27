# =============================================================================
# 03_mobility_matrices.R — Mobility Matrix Construction
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: Build all 13 mobility matrix variants (M1, M2a, M2b, M3, M4, M4b, M5,
#   M6a, M6b, M7, M8, M9, M10; M11 inward-FOI is built on demand in run_all) that
#   parameterise between-zone transmission in the spatiotemporal metapopulation
#   model.
#
# Matrices are 519×519 (national health zones from WorldPop), named by
# canonical zone name, and row-stochastic (each row sums to ≤ 1; zero rows
# are permitted for zones with no observed outflow).
#
# Data sources:
#   M1 / M2a / M2b : Flowminder short-trip cohort (epicentre pooled)
#   M3             : Flowminder full origin-destination matrix
#   M4             : Negative-binomial gravity model (calibrated on M3)
#   M5             : Radiation model (Simini et al. 2012 Nature 484:96-100)
#   M6a / M6b      : OSRM travel-time decay (exponential / power-law)
#   M7             : IDP-augmented hybrid (M3 + IDP flows)
#   M8             : Composite — short trips for epicentre, gravity elsewhere (RECOMMENDED)
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

`%||%` <- function(a, b) if (!is.null(a)) a else b

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
#' First applies known aliases, then attempts direct matching.  Unresolved
#' names are returned as-is with a warning.
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

#' Load OSRM travel-time matrix.
#'
#' The file has: first column = `nom` (origin), remaining columns = destination
#' zone names.  Values are travel times in minutes; 0 on the diagonal; some
#' entries may be NA (unroutable pairs).
#'
#' @return Square numeric matrix with rownames = colnames = zone names.
#' Load an OSRM zone-to-zone cost matrix.
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
  M
}

#' Load Flowminder full OD matrix (outflow perspective).
#'
#' @return Raw (unnormalised) square numeric matrix, 438×438.
load_flowminder_od <- function() {
  f <- file.path(FLOWMINDER_DIR, "flowminder__outflow__static.matrix.csv")
  stopifnot(file.exists(f))

  df <- readr::read_csv(f, col_types = readr::cols(.default = "d",
                                                    nom      = "c"),
                         show_col_types = FALSE)
  zones <- df$nom
  M     <- as.matrix(df[, -which(colnames(df) == "nom")])
  rownames(M) <- zones
  message(sprintf("[flowminder_od] Loaded %dx%d OD matrix.", nrow(M), ncol(M)))
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
#' @return Named numeric vector: names = destination zone, values = avg presence days.
load_cohort_snapshot <- function(cohort, window = "followup") {
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
  f <- file.path(FLOWMINDER_ST_DIR, fs[which.max(as.integer(tags))])

  df <- readr::read_csv(f, col_types = readr::cols(.default = "d", nom = "c"),
                        show_col_types = FALSE)
  stopifnot("nom" %in% colnames(df))
  dest_cols <- setdiff(colnames(df), "nom")
  # All origin rows are the same cohort aggregate — use only the first data row.
  vals <- as.numeric(df[1L, dest_cols])
  names(vals) <- dest_cols
  message(sprintf(
    "[cohort] %s/%s (tag %s): %d destinations, first origin='%s', total presence-days=%.3f",
    cohort, window, max(tags), length(vals), df$nom[1L], sum(vals, na.rm = TRUE)
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

  message(sprintf("[idp] Loaded IDP static matrix: %d origins × %d destinations.",
                  nrow(M_raw), ncol(M_raw)))
  M_raw
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
  if (!length(.elig)) .elig <- seq_along(FLOWMINDER_ST_DATES)   # all in the future -> earliest available
  tag <- FLOWMINDER_ST_TAGS[.elig[which.max(FLOWMINDER_ST_DATES[.elig])]]
  message(sprintf("[M1] Building Flowminder short-trip matrix (latest snapshot <= %s: %s)...",
                  format(analysis_date), tag))
  st_props    <- load_short_trip_snapshot(tag)

  # Harmonise destination names in the short-trip data
  dest_names_raw <- names(st_props)
  dest_names_can <- harmonise_names(dest_names_raw, aliases, zones_all)

  # Report dropped destinations
  in_zones   <- dest_names_can %in% zones_all
  n_dropped  <- sum(!in_zones)
  if (n_dropped > 0L) {
    message(sprintf("[M1] Dropping %d short-trip destinations not in zones_all: %s",
                    n_dropped,
                    paste(dest_names_raw[!in_zones], collapse = ", ")))
  }

  # Initialise with fallback or zeros
  W <- if (!is.null(fallback_M3)) {
    stopifnot(identical(dim(fallback_M3), c(length(zones_all), length(zones_all))))
    fallback_M3
  } else {
    make_zero_matrix(zones_all)
  }

  # Epicentre rows: overwrite with short-trip proportions
  for (ez in epicentre_zones) {
    ez_can <- harmonise_names(ez, aliases, zones_all)
    if (!ez_can %in% zones_all) {
      warning(sprintf("[M1] Epicentre zone '%s' not found in zones_all; skipping.", ez))
      next
    }
    # Zero out this row first, then fill destinations
    W[ez_can, ] <- 0
    for (i in seq_along(dest_names_can)) {
      dest_c <- dest_names_can[i]
      if (dest_c %in% zones_all) {
        # ACCUMULATE (not overwrite): harmonise_names() can map several observed
        # destinations onto one canonical zone; overwriting would drop the earlier
        # share and understate that zone's outflow. Row was zeroed above, so += is safe.
        W[ez_can, dest_c] <- W[ez_can, dest_c] + st_props[i]
      }
    }
    # Ensure no self-loop
    W[ez_can, ez_can] <- 0
    message(sprintf("[M1] Epicentre '%s': total outflow prop = %.4f", ez_can,
                    sum(W[ez_can, ], na.rm = TRUE)))
  }

  # Row-normalise (normalises epicentre rows that may not sum to exactly 1,
  # and normalises fallback rows that are already stochastic — idempotent)
  W <- make_row_stochastic(W)

  assert_mobility_matrix(W, zones_all, "M1")
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

  # Load and average proportions across selected snapshots
  all_props <- lapply(use_tags, load_short_trip_snapshot)
  # Align on common destination names
  all_dests <- Reduce(intersect, lapply(all_props, names))
  avg_props  <- rowMeans(do.call(cbind, lapply(all_props, \(p) p[all_dests])))

  # Build matrix in same way as M1 but with averaged proportions
  dest_names_can <- harmonise_names(all_dests, aliases, zones_all)

  W <- if (!is.null(fallback_M3)) fallback_M3 else make_zero_matrix(zones_all)

  for (ez in epicentre_zones) {
    ez_can <- harmonise_names(ez, aliases, zones_all)
    if (!ez_can %in% zones_all) next
    W[ez_can, ] <- 0
    for (i in seq_along(dest_names_can)) {
      dest_c <- dest_names_can[i]
      # ACCUMULATE (not overwrite): several observed destinations can harmonise to one
      # canonical zone; overwriting would drop the earlier share. Row zeroed above.
      if (dest_c %in% zones_all) W[ez_can, dest_c] <- W[ez_can, dest_c] + avg_props[i]
    }
    W[ez_can, ez_can] <- 0
  }

  W <- make_row_stochastic(W)
  assert_mobility_matrix(W, zones_all, variant)
  message(sprintf("[%s] Built: sparsity=%.1f%%", variant, 100 * mean(W == 0)))
  W
}

# ---------------------------------------------------------------------------
# M3: Flowminder full OD — national, row-normalised
# ---------------------------------------------------------------------------

#' Build M3: national Flowminder OD matrix embedded in the 519-zone space.
#'
#' Zones in the 438-zone OD matrix are matched by canonical name.  Zones
#' not present in the OD matrix get zero rows.
#'
#' @param zones_all  Canonical zone vector.
#' @param aliases    Alias lookup.
#' @return Row-stochastic 519×519 matrix.  Also invisibly returns the
#'         raw (unnormalised) matrix as attribute "raw".
build_M3 <- function(zones_all, aliases) {
  message("[M3] Building Flowminder full OD matrix...")

  M_raw <- load_flowminder_od()   # 438×438

  # Harmonise names
  orig_names_can <- harmonise_names(rownames(M_raw), aliases, zones_all)
  dest_names_can <- harmonise_names(colnames(M_raw), aliases, zones_all)

  rownames(M_raw) <- orig_names_can
  colnames(M_raw) <- dest_names_can

  # Aggregate (SUM) duplicate canonical rows/cols rather than dropping all but the
  # first occurrence: harmonise_names() can map several OD zones onto one canonical
  # name, and keeping only the first would silently discard the others' entire flow.
  # rowsum() groups by name and sums; names stay unique and name-indexable below.
  n_dup_rows <- sum(duplicated(rownames(M_raw)))
  n_dup_cols <- sum(duplicated(colnames(M_raw)))
  if (n_dup_rows > 0L) {
    message(sprintf("[M3] Aggregating %d duplicate origin row(s) by summation.", n_dup_rows))
    M_raw <- rowsum(M_raw, group = rownames(M_raw), reorder = FALSE)
  }
  if (n_dup_cols > 0L) {
    message(sprintf("[M3] Aggregating %d duplicate destination column(s) by summation.", n_dup_cols))
    M_raw <- t(rowsum(t(M_raw), group = colnames(M_raw), reorder = FALSE))
  }

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

  shared_orig <- od_origins[in_zones_o]
  shared_dest <- od_dests[in_zones_d]

  for (oz in shared_orig) {
    for (dz in shared_dest) {
      W_full[oz, dz] <- M_raw[oz, dz]
    }
  }
  diag(W_full) <- 0

  # Store raw for gravity model calibration
  M_raw_519 <- W_full   # raw (unnormalised) 519×519

  W_norm <- make_row_stochastic(W_full)
  assert_mobility_matrix(W_norm, zones_all, "M3")
  message(sprintf("[M3] Built: sparsity=%.1f%%", 100 * mean(W_norm == 0)))

  attr(W_norm, "raw") <- M_raw_519
  W_norm
}

# ---------------------------------------------------------------------------
# M15: Flowminder combined inflow+outflow static kernel (symmetrised total flow)
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

#' Build M15: Flowminder STATIC kernel from the combined inflow + outflow OD data.
#'
#' Flowminder's processed inflow and outflow exports are the SAME directed
#' origin->destination table (verified byte-identical), so the information-adding
#' combination of "inflow + outflow" is the SYMMETRISED total two-way flow:
#'
#'   S[j,i] = O[j,i] + O[i,j]   (trips OUT of j to i  +  trips IN to j from i)
#'
#' computed as `Oe + t(Ie)` in the 519 canonical space, where Oe/Ie are the
#' outflow/inflow OD tables embedded by name (so t() is well-defined and pairs are
#' matched by zone, never by row position — the raw tables' row and column orders
#' differ). This is genuinely distinct from M3 (directed outflow, row-normalised as
#' is) and materially less sparse: reciprocal edges present in only one direction are
#' filled. If the inflow file is absent, t(Oe) is used as the inflow counterpart
#' (identical here, and the correct symmetrisation regardless).
#'
#' @param zones_all Canonical zone vector.
#' @param aliases   Alias lookup.
#' @return Row-stochastic 519x519 matrix; attribute "raw" is the unnormalised S.
build_M15 <- function(zones_all, aliases) {
  message("[M15] Building Flowminder combined inflow+outflow static kernel (symmetrised total flow)...")

  Oe <- .embed_od_to_519(load_flowminder_od(), zones_all, aliases, "M15/outflow")

  f_in <- file.path(FLOWMINDER_DIR, "flowminder__inflow__static.matrix.csv")
  Ie <- if (file.exists(f_in)) {
    df <- readr::read_csv(f_in, col_types = readr::cols(.default = "d", nom = "c"),
                          show_col_types = FALSE)
    Min <- as.matrix(df[, -which(colnames(df) == "nom")]); rownames(Min) <- df$nom
    message(sprintf("[M15] Loaded %dx%d inflow OD table.", nrow(Min), ncol(Min)))
    .embed_od_to_519(Min, zones_all, aliases, "M15/inflow")
  } else {
    message("[M15] Inflow file absent; using the outflow transpose as the inflow counterpart.")
    Oe
  }

  # Total two-way exchange between each pair: outflow (j->i) + inflow into j from i (= t(inflow OD)).
  S <- Oe + t(Ie)
  diag(S) <- 0

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

#' Build M4: gravity model calibrated on Flowminder OD flows.
#'
#' Fits an unconstrained gravity GLM (negative-binomial, falling back to Poisson
#' if the NB theta iteration diverges) with log(pop_i), log(pop_j) and
#' log(dist_ij + 1).  Predicts for all 519×519 pairs and row-normalises.  A
#' per-origin intercept would cancel under row-normalisation, so origin size is
#' captured by the single log(pop_i) covariate rather than fixed effects.
#'
#' @param M3_raw   Raw (unnormalised) 519×519 flow matrix (from attr(M3, "raw")).
#' @param pop_vec  Named population vector (length 519).
#' @param osrm_mat OSRM travel-time matrix (zones_all × zones_all, minutes).
#' @param zones_all Canonical zone vector.
#' @return Row-stochastic 519×519 matrix.  Attribute "model_summary" contains
#'         the NB-GLM summary.
build_M4 <- function(M3_raw, pop_vec, osrm_mat, zones_all,
                     deterrence = c("power", "exp"), label = "M4") {
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

  train_rows <- lapply(train_zones, function(oz) {
    for_dests <- zones_all[zones_all != oz & pop_aligned > 0]
    dists     <- dist_mat[oz, for_dests]
    valid     <- !is.na(dists) & dists > 0
    for_dests <- for_dests[valid]
    dists     <- dists[valid]
    if (length(for_dests) == 0L) return(NULL)
    tibble::tibble(
      origin  = oz,
      dest    = for_dests,
      flow    = as.numeric(M3_raw[oz, for_dests]),
      pop_i   = pop_aligned[oz],
      pop_j   = pop_aligned[for_dests],
      dist_ij = dists
    )
  })
  train_df <- dplyr::bind_rows(train_rows)
  train_df$flow <- pmax(round(train_df$flow), 0L)   # ensure non-negative integer

  message(sprintf("[%s] Training set: %d OD pairs (%.1f%% zero flows).",
                  label, nrow(train_df), 100 * mean(train_df$flow == 0)))

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
        # Fallback: grand mean of positive flows
        rep(mean(train_df$flow[train_df$flow > 0]), length(dests))
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

  attr(W_norm, "model_summary") <- summary(nb_fit)
  W_norm
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
    W_raw[!is.finite(W_raw)] <- 0
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

#' Build M7: M3 flows augmented by IDP displacement flows.
#'
#' Augmented flow: F_aug[i,j] = M3_raw[i,j] + theta * IDP[i,j].
#' Row-normalised to give a probability matrix.
#'
#' @param M3_raw    Raw (unnormalised) 519×519 M3 flow matrix.
#' @param pop_vec   Named population vector.
#' @param osrm_mat  OSRM travel-time matrix.
#' @param zones_all Canonical zone vector.
#' @param aliases   Alias lookup.
#' @param theta     IDP weight scalar (default 1.0).
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

  message(sprintf("[M7] IDP matrix embedded: %d origin zones with IDP flows.",
                  sum(rowSums(W_idp) > 0)))

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
# M8: Composite (RECOMMENDED)
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
build_M8 <- function(M1, M4, epicentre_zones, zones_all) {
  message("[M8] Building composite matrix (epicentre=M1, elsewhere=M4)...")

  stopifnot(identical(dim(M1), c(length(zones_all), length(zones_all))))
  stopifnot(identical(dim(M4), c(length(zones_all), length(zones_all))))

  W <- M4
  for (ez in epicentre_zones) {
    if (!ez %in% zones_all) {
      warning(sprintf("[M8] Epicentre zone '%s' not in zones_all; skipping.", ez))
      next
    }
    W[ez, ] <- M1[ez, ]
  }

  diag(W) <- 0
  # Re-normalise to be safe (M1 and M4 rows should already be stochastic)
  W <- make_row_stochastic(W)

  assert_mobility_matrix(W, zones_all, "M8")
  message(sprintf("[M8] Built: sparsity=%.1f%%", 100 * mean(W == 0)))
  W
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
compose_epicentre <- function(base, epi_mat, epicentre_zones, zones_all, label) {
  stopifnot(identical(dim(base), c(length(zones_all), length(zones_all))),
            identical(dim(epi_mat), dim(base)))
  W <- base
  for (ez in epicentre_zones) {
    if (!ez %in% zones_all) {
      warning(sprintf("[%s] Epicentre zone '%s' not in zones_all; skipping.", label, ez))
      next
    }
    W[ez, ] <- epi_mat[ez, ]
  }
  diag(W) <- 0
  W <- make_row_stochastic(W)
  assert_mobility_matrix(W, zones_all, label)
  message(sprintf("[%s] Built: sparsity=%.1f%%", label, 100 * mean(W == 0)))
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
build_M_cohort <- function(zones_all, cohort_sources, aliases, window = "followup") {
  message(sprintf("[M_cohort] Building pooled cohort presence kernel (window=%s)...", window))
  W <- make_zero_matrix(zones_all)
  filled <- character(0)

  for (cohort in names(cohort_sources)) {
    origins_can <- harmonise_names(cohort_sources[[cohort]], aliases, zones_all)
    vals        <- load_cohort_snapshot(cohort, window)
    dest_can    <- harmonise_names(names(vals), aliases, zones_all)
    keep        <- dest_can %in% zones_all
    if (any(!keep))
      message(sprintf("[M_cohort] %s: dropping %d destination(s) not in zones_all: %s",
                      cohort, sum(!keep),
                      paste(utils::head(names(vals)[!keep], 5L), collapse = ", ")))

    for (ez in origins_can) {
      if (!ez %in% zones_all) {
        warning(sprintf("[M_cohort] Origin '%s' (cohort %s) not in zones_all; skipping.",
                        ez, cohort)); next
      }
      W[ez, ] <- 0
      for (i in seq_along(dest_can)) {
        # ACCUMULATE (not overwrite): harmonise_names() can map several observed
        # destinations onto one canonical zone; overwriting would drop the earlier
        # share. Row zeroed above, so += is safe. (Mirrors build_M1.)
        if (keep[i]) W[ez, dest_can[i]] <- W[ez, dest_can[i]] + vals[i]
      }
      W[ez, ez] <- 0   # no self-loop
      filled <- c(filled, ez)
    }
  }

  W <- make_row_stochastic(W)
  assert_mobility_matrix(W, zones_all, "M_cohort")
  filled <- unique(filled)
  if (!length(filled))
    stop("[M_cohort] No cohort origins resolved to zones_all — check COHORT_SOURCES.")
  attr(W, "cohort_origins") <- filled
  message(sprintf("[M_cohort] Filled %d cohort-origin rows: %s",
                  length(filled), paste(filled, collapse = ", ")))
  W
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
build_M10 <- function(M1, M5, epicentre_zones, zones_all) {
  message("[M10] Building radiation composite (epicentre=M1, elsewhere=M5)...")
  compose_epicentre(M5, M1, epicentre_zones, zones_all, "M10")
}

# ---------------------------------------------------------------------------
# M17: Grand all-kernel consensus ensemble
# ---------------------------------------------------------------------------

#' Build M17: the grand all-kernel consensus ensemble.
#'
#' Averages EVERY distinct structural / empirical base hypothesis for non-source
#' flow element-wise — Flowminder directed OD (M3), calibrated gravity (M4),
#' radiation (M5), travel-time decay (M6a) and the combined Flowminder static
#' kernel (M15). The mean of row-stochastic kernels is itself a convex combination
#' (row-stochastic up to zero-outflow rows), extending the M9 three-kernel ensemble
#' to the full set of mobility data sources; multi-model consensus kernels are
#' typically better calibrated than any single member (Cramer et al. 2022 PNAS
#' 119:e2113561119; Reich et al. 2019 PNAS 116:3146). The empirical SOURCE rows are
#' then overlaid from the best-available data — the pooled Flowminder cohort presence
#' (Ituri/NK/Tshopo) when available, else the short-trip epicentre kernel (M1) — via
#' compose_epicentre(), which re-normalises the whole matrix so a consensus row that
#' drew mass from a differing number of members becomes a proper distribution.
#' Composite kernels (M8/M9/M10/M13/M14) are deliberately EXCLUDED from the average to
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
                      label = "M17", weights = NULL) {
  n <- length(zones_all)
  keep <- vapply(bases, function(M)
    is.matrix(M) && all(dim(M) == c(n, n)) && identical(rownames(M), zones_all),
    logical(1))
  bases <- bases[keep]
  if (length(bases) < 2L)
    stop(sprintf("[%s] need >= 2 valid base kernels to form a consensus; got %d.",
                 label, length(bases)))
  if (is.null(weights)) weights <- rep(1, length(bases))
  w <- weights[seq_along(bases)]; w <- w / sum(w)
  message(sprintf("[%s] Averaging %d base kernels (%s)...",
                  label, length(bases), paste(names(bases), collapse = "+")))
  base <- Reduce(`+`, Map(function(M, wi) wi * M, bases, w))   # convex combo -> row sums <= 1
  compose_epicentre(base, epi_mat, epi_origins, zones_all, label)
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

  # Top-5 destinations for each epicentre zone
  top5_list <- lapply(epicentre_zones, function(ez) {
    if (!ez %in% zones_all) return(NA_character_)
    row_w <- W[ez, ]
    top5  <- names(sort(row_w, decreasing = TRUE))[1:5]
    paste(top5, collapse = "; ")
  })
  names(top5_list) <- epicentre_zones

  tibble::tibble(
    matrix_id    = matrix_id,
    n_zones      = nrow(W),
    sparsity_pct = 100 * mean(W == 0),
    mean_weight  = mean(W[W > 0]),
    n_zero_rows  = sum(rs < 1e-9),
    top5_Bunia   = top5_list[["Bunia"]] %||% NA_character_,
    top5_Mongbalu = top5_list[["Mongbalu"]] %||% NA_character_,
    top5_Rwampara = top5_list[["Rwampara"]] %||% NA_character_
  )
}

# ---------------------------------------------------------------------------
# Master function
# ---------------------------------------------------------------------------

#' Build all mobility matrix variants (M1, M2a, M2b, M3, M4, M4b, M5, M6a, M6b, M7,
#' M8, M9, M10; Flowminder-cohort composites M13/M14; combined-Flowminder-static M15,
#' cohort+static M16, all-kernel consensus M17; and the road-distance -dist twins).
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
#'   (cohort composites), M15 (combined-Flowminder static), M16 (cohort + static),
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

  out_list    <- list()
  summary_rows <- list()

  # Helper: save and summarise
  save_matrix <- function(W, id) {
    rds_path <- file.path(OUT_MOBILITY, sprintf("mobility_%s.rds", id))
    saveRDS(W, rds_path)
    message(sprintf("[save] %s → %s", id, basename(rds_path)))
    out_list[[id]]     <<- W
    summary_rows[[id]] <<- summarise_matrix(W, id, zones_all, epicentre_can)
  }

  # ---- M3 first (needed for M1/M2 fallback, M4, M7) ----
  M3 <- build_M3(zones_all, aliases)
  M3_raw <- attr(M3, "raw")
  save_matrix(M3, "M3")

  # ---- M1 ----
  M1 <- build_M1(zones_all, epicentre_can, aliases, fallback_M3 = M3)
  save_matrix(M1, "M1")

  # ---- M2a ----
  M2a <- build_M2(zones_all, epicentre_can, aliases,
                   analysis_date = analysis_date, variant = "M2a",
                   fallback_M3 = M3)
  save_matrix(M2a, "M2a")

  # ---- M2b ----
  M2b <- build_M2(zones_all, epicentre_can, aliases,
                   analysis_date = analysis_date, variant = "M2b",
                   fallback_M3 = M3)
  save_matrix(M2b, "M2b")

  # ---- M4 gravity (power deterrence) ----
  M4 <- build_M4(M3_raw, pop_vec, osrm_mat, zones_all)
  save_matrix(M4, "M4")

  # ---- M4b gravity (exponential deterrence) ----
  M4b <- build_M4(M3_raw, pop_vec, osrm_mat, zones_all,
                  deterrence = "exp", label = "M4b")
  save_matrix(M4b, "M4b")

  # ---- M5 radiation ----
  M5 <- build_M5(pop_vec, osrm_mat, zones_all)
  save_matrix(M5, "M5")

  # ---- M6a exponential decay ----
  M6a <- build_M6(osrm_mat, zones_all, variant = "exp", kappa = 120)
  save_matrix(M6a, "M6a")

  # ---- M6b power-law decay ----
  M6b <- build_M6(osrm_mat, zones_all, variant = "power", gamma = 1)
  save_matrix(M6b, "M6b")

  # ---- M7 IDP-augmented ----
  M7 <- build_M7(M3_raw, pop_vec, osrm_mat, zones_all, aliases, theta = 1.0)
  save_matrix(M7, "M7")

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

  # ---- M15: Flowminder combined inflow+outflow static kernel (symmetrised total flow) ----
  # A base kernel using BOTH Flowminder OD perspectives (S = O + t(O) in canonical space).
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
  if (want_cohort) {
    M_cohort <- tryCatch(
      build_M_cohort(zones_all, cohort_sources, aliases,
                     window = get0("COHORT_WINDOW", ifnotfound = "followup")),
      error = function(e) { warning("[M_cohort] build failed: ", conditionMessage(e)); NULL })
    if (!is.null(M_cohort)) {
      cohort_origins <- attr(M_cohort, "cohort_origins")
      save_matrix(compose_epicentre(M4, M_cohort, cohort_origins, zones_all, "M13"), "M13")
      save_matrix(compose_epicentre(M5, M_cohort, cohort_origins, zones_all, "M14"), "M14")
      # ---- M16: cohort + directed Flowminder OD composite (cohort where available, M3 elsewhere) ----
      # The direct analogue of M13/M14 but pairing the cohort presence rows with the raw empirical
      # DIRECTED Flowminder OD flows (M3) instead of a gravity/radiation model — no symmetrisation.
      if (want_flowstat)
        save_matrix(compose_epicentre(M3, M_cohort, cohort_origins, zones_all, "M16"), "M16")
    }
  }

  # ---- M17: grand all-kernel consensus ensemble (travel-time base) ----
  # Mean of the distinct structural/empirical bases {M3, M4, M5}, with the empirical source rows
  # overlaid from the cohort kernel (when available) else the short-trip epicentre M1. The
  # travel-time-decay kernel (M6a) and symmetrised static (M15) are deliberately EXCLUDED. This is
  # the "combine all mobility matrices" ensemble option; gated by INCLUDE_FLOWSTATIC.
  if (want_flowstat) {
    .m17_epi_mat    <- if (!is.null(M_cohort)) M_cohort else M1
    .m17_epi_origin <- if (!is.null(M_cohort)) cohort_origins else epicentre_can
    M17 <- tryCatch(
      build_M17(list(M3 = M3, M4 = M4, M5 = M5),
                .m17_epi_mat, .m17_epi_origin, zones_all, "M17"),
      error = function(e) { warning("[M17] build failed: ", conditionMessage(e)); NULL })
    if (!is.null(M17)) save_matrix(M17, "M17")
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
    M4_d <- build_M4(M3_raw, pop_vec, osrm_dist_mat, zones_all, label = "M4-dist")
    M5_d <- build_M5(pop_vec, osrm_dist_mat, zones_all)
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
    }
    if (build_cohort_dist) {
      save_matrix(compose_epicentre(M4_d, M_cohort, cohort_origins, zones_all, "M13-dist"), "M13-dist")
      save_matrix(compose_epicentre(M5_d, M_cohort, cohort_origins, zones_all, "M14-dist"), "M14-dist")
    }
    if (build_flowstat_dist) {
      # M17-dist: all-kernel consensus on ROAD-KM over the bases {M3, M4-dist, M5-dist}. Gravity and
      # radiation are re-keyed on km (M4_d/M5_d); the empirical directed Flowminder OD (M3) carries no
      # distance axis, so it is reused as-is. The travel-time decay (M6a) and symmetrised static (M15)
      # are excluded, matching the travel-time M17. Empirical source rows from the cohort kernel else M1.
      .m17d_epi_mat    <- if (!is.null(M_cohort)) M_cohort else M1
      .m17d_epi_origin <- if (!is.null(M_cohort)) cohort_origins else epicentre_can
      save_matrix(build_M17(list(M3 = M3, `M4-dist` = M4_d, `M5-dist` = M5_d),
                            .m17d_epi_mat, .m17d_epi_origin, zones_all, "M17-dist"), "M17-dist")
    }
  }

  # ---- Summary CSV ----
  summary_df <- dplyr::bind_rows(summary_rows)
  summary_path <- file.path(OUT_MOBILITY, "mobility_summary.csv")
  readr::write_csv(summary_df, summary_path)
  message(sprintf("[summary] Written to %s", summary_path))

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
