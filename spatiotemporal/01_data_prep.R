# =============================================================================
# 01_data_prep.R — Spatiotemporal Data Preparation
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: Load, clean, and harmonise all data sources required for the
#   spatiotemporal modelling suite (scripts 05–16). Returns a single named
#   list `dat` containing:
#     dat$ll          — cleaned DHIS2 linelist tibble
#     dat$zone_week   — zone × week observation matrix (7-day buckets ending on ANALYSIS_DATE; zero-filled)
#     dat$pop         — named numeric vector: zone → WorldPop population count
#     dat$covariates  — static covariate tibble joined on `nom`
#     dat$sitrep      — INSP sitrep weekly zone tibble (optional)
#
# Run order: source 00_config.R first (done below), then source this file.
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(readr)
  library(jsonlite)
})

# =============================================================================
# SECTION 1: Helper utilities
# =============================================================================

# Parse a logical column that may arrive as character "TRUE"/"FALSE", logical,
# or integer 0/1. Returns a logical vector.
.parse_logical_col <- function(x) {
  if (is.logical(x)) return(x)
  raw <- trimws(as.character(x))
  case_when(
    raw %in% c("TRUE",  "true",  "1", "yes") ~ TRUE,
    raw %in% c("FALSE", "false", "0", "no")  ~ FALSE,
    TRUE                                       ~ NA
  )
}

# Safe date parser: handles character, Date, POSIXt.
# Returns a Date vector; unparseable values become NA with a warning-free path.
.parse_date <- function(x) {
  if (inherits(x, "Date"))   return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  raw <- trimws(as.character(x))
  raw[raw == ""] <- NA_character_
  # Parse PER ELEMENT, not per vector: base as.Date(tryFormats=) picks ONE format
  # from the first non-NA element and applies it to the whole column, silently
  # NA-ing every value in a different format. Instead try each format in priority
  # order (ISO first — unambiguous — then day/month slash formats) and coalesce,
  # so a mixed-format column keeps all parseable values.
  fmts <- c("%Y-%m-%d", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S", "%d/%m/%Y", "%m/%d/%Y")
  out <- as.Date(rep(NA_real_, length(raw)), origin = "1970-01-01")
  for (fmt in fmts) {
    todo <- is.na(out) & !is.na(raw)
    if (!any(todo)) break
    out[todo] <- suppressWarnings(as.Date(raw[todo], format = fmt))
  }
  out
}

# Apply zone name aliases from data/aliases.csv to a character vector of
# zone names. `aliases` is the aliases tibble (observed_name, canonical_nom).
.apply_aliases <- function(zone_vec, aliases) {
  # A duplicated observed_name would resolve silently to whichever row came first — and
  # differently at the two call sites, which pass differently-filtered alias tables.
  if (!is.null(aliases) && nrow(aliases) && anyDuplicated(aliases$observed_name))
    warning(sprintf("[aliases] %d duplicated observed_name(s) (e.g. %s); the FIRST mapping wins and may differ between call sites.",
                    sum(duplicated(aliases$observed_name)),
                    paste(utils::head(unique(aliases$observed_name[duplicated(aliases$observed_name)]), 3),
                          collapse = ", ")), call. = FALSE)
  lookup <- setNames(aliases$canonical_nom, aliases$observed_name)
  out <- zone_vec
  mask <- zone_vec %in% names(lookup)
  out[mask] <- lookup[zone_vec[mask]]
  out
}

# Check that `required` column names all exist in `df`. Stops with an
# informative message if any are absent.
.check_cols <- function(df, required, context) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      "[", context, "] Missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# NOTE: .load_dhis2_delay_params() and .draw_dhis2_delay() MOVED to 00_config.R (sourced
# above), alongside delay_cdf() and effective_onset_sample_delay(). They lived here, but
# 04_nowcasting.R and 02_epi_params.R need the same delay and could not reach a definition
# private to the data-prep module — so they silently kept using the fixed lab Exp reference
# while this file used the fitted DHIS2 delay. One shared definition removes that whole class
# of drift: the imputation draw (.draw_dhis2_delay), the nowcast weights (delay_cdf) and the
# R(t) truncation model now all read ONE resolver, effective_onset_sample_delay().


# Build sitrep-confirmed TOP-UP line-list rows that reconcile the DHIS2 line list
# to the INSP sitrep's cumulative confirmed counts, WITHOUT double-counting.
#
# WHY: invasion status and case counts in this suite are defined purely from the
# line list's confirmed cases, but the line list lags the official INSP sitrep for
# some zones (e.g. a zone the sitrep confirms while the line list still holds only
# suspects). This appends the SHORTFALL — for each zone, max(0, sitrep_cumulative -
# linelist_confirmed) confirmed rows — so every zone reaches AT LEAST its sitrep
# confirmed count. Zones where the line list already meets/exceeds the sitrep get
# nothing (no cases are ever removed; the sitrep is a floor, never a ceiling).
#
# SOURCE: the CUMULATIVE confirmed file is authoritative for magnitude — the daily
# `new_confirmed_cases` extraction grossly undercounts (e.g. Bunia 66 vs 507), so it
# is NOT used here. Per canonical zone we merge spelling variants by MAX-per-date,
# take the running maximum over dates (absorbing the sitrep's occasional cumulative
# dips/revisions), and diff it into dated positive increments — one confirmed case
# per unit increment, dated at the sitrep report date on which it appeared.
#
# PLACEMENT & DATING: the shortfall is filled with the MOST RECENT sitrep confirmed
# dates (reporting lag ⇒ the line list already holds the earlier cases). Each appended
# row's `date_of_sample_collection` is set to that sitrep report date and its onset is
# LEFT BLANK, so the caller's onset-imputation block imputes onset = sample − a delay
# DRAWN from the fitted onset→sample distribution — exactly as it treats the ~15% of
# real records lacking an onset (same convention as build_conditional_linelist.R).
#
# IDEMPOTENT vs the conditional pre-append: build_conditional_linelist.R already
# injects some sitrep zones into the conditional CSV; those raise `n_ll` so their
# shortfall is 0 here — they are not re-appended. Toggle with APPEND_SITREP_CONFIRMED.
#
# Args:  ll          — current line list (outbreak-filtered, unmatched-excluded,
#                       health_zone alias-corrected), with final_mve_case_classification.
#        aliases_all — full aliases tibble (observed_name, canonical_nom) or NULL;
#                       used to canonicalise sitrep `nom` to the modelling spine.
# Returns a tibble (0 rows if nothing to append) whose columns are a subset of `ll`'s,
#   safe to dplyr::bind_rows() onto `ll`.
#' The INSP sitrep's CUMULATIVE CONFIRMED series, canonicalised onto the modelling spine.
#'
#' ONE reader, used by the line-list reconciliation below AND by the surveillance-stream
#' figures, because "what the sitrep says" has to mean the same thing in both. The steps are
#' not cosmetic and a second implementation would not reproduce them:
#'   * the CUMULATIVE file, never the daily `new_confirmed_cases` one, which grossly
#'     undercounts (Bunia 66 against 507 on the 2026-07 snapshot);
#'   * zone names are the FRENCH `nom` and are canonicalised through aliases.csv — the spine
#'     writes "Nia Nia" where the sitrep writes "Nia-Nia", and an un-aliased join drops it;
#'   * as-of filtered to ANALYSIS_DATE, so nothing downstream can see past the snapshot;
#'   * per (zone, date) the MAX across spelling variants, then a RUNNING MAX over dates — the
#'     published cumulative series is not monotone (it is re-stated), and a bare diff of a
#'     non-monotone cumulative column yields negative "new cases".
#'
#' @param aliases_all the alias table (may be NULL: names are then left as published).
#' @return tibble(nom, date, cum, inc) — cum is the running-max cumulative count and inc its
#'   dated positive increment — or an empty tibble when the file is absent or unusable.
sitrep_cumulative_confirmed <- function(aliases_all = NULL) {
  cum_path <- file.path(SITREP_DIR, "insp_sitrep__cumulative_confirmed_cases__daily.csv")
  if (!file.exists(cum_path)) {
    warning("[sitrep] cumulative-confirmed file not found: ", cum_path, call. = FALSE)
    return(tibble::tibble())
  }
  cum <- tryCatch(
    readr::read_csv(cum_path,
      col_types = readr::cols(nom = "c", date = "c", cumulative_confirmed_cases = "c"),
      show_col_types = FALSE),
    error = function(e) { warning("[sitrep] read error: ", e$message, call. = FALSE); NULL })
  if (is.null(cum) || !all(c("nom", "date", "cumulative_confirmed_cases") %in% names(cum)))
    return(tibble::tibble())

  .asof <- suppressWarnings(as.Date(get0("ANALYSIS_DATE", ifnotfound = NA)))

  cum <- cum %>%
    dplyr::mutate(
      nom  = trimws(nom),
      date = .parse_date(date),
      cum  = suppressWarnings(as.numeric(cumulative_confirmed_cases))
    ) %>%
    dplyr::filter(!is.na(nom), nom != "", nom != "NA",
                  !is.na(date), !is.na(cum), is.finite(cum),
                  date >= OUTBREAK_START,
                  is.na(.asof) | date <= .asof)                      # as-of: no future leakage
  if (nrow(cum) == 0) return(tibble::tibble())

  if (!is.null(aliases_all) && nrow(aliases_all) > 0)               # canonicalise sitrep names
    cum$nom <- .apply_aliases(cum$nom, aliases_all)

  # Per canonical zone × date: MAX cumulative across spelling variants; then running max
  # over dates; then diff into dated positive integer increments.
  cum %>%
    dplyr::group_by(nom, date) %>%
    dplyr::summarise(cum = max(cum), .groups = "drop") %>%
    dplyr::arrange(nom, date) %>%
    dplyr::group_by(nom) %>%
    dplyr::mutate(cum = cummax(cum),
                  inc = as.integer(round(cum - dplyr::lag(cum, default = 0)))) %>%
    dplyr::ungroup()
}

.build_sitrep_confirmed_appends <- function(ll, aliases_all) {

  if (!isTRUE(get0("APPEND_SITREP_CONFIRMED", ifnotfound = TRUE))) return(tibble::tibble())

  # Canonical modelling spine (WorldPop 519 zones): a sitrep zone absent from it cannot
  # enter the zone-week grid, so we skip (and warn about) it rather than fabricate an
  # unmatched row that aggregate_to_zone_week would silently drop.
  spine <- tryCatch(
    readr::read_csv(file.path(WORLDPOP_DIR, "worldpop__pop_count__static.csv"),
      col_types = readr::cols(nom = "c", pop_count = "d"), show_col_types = FALSE)$nom,
    error = function(e) NULL)

  incs <- sitrep_cumulative_confirmed(aliases_all)
  if (!nrow(incs)) return(tibble::tibble())
  incs <- incs %>% dplyr::filter(inc > 0L) %>% dplyr::select(nom, date, inc)
  if (nrow(incs) == 0) return(tibble::tibble())

  sit_tot <- incs %>% dplyr::group_by(nom) %>%
    dplyr::summarise(n_sit = sum(inc), .groups = "drop")

  ll_tot <- ll %>%
    dplyr::filter((final_mve_case_classification == CONFIRMED_STATUS) %in% TRUE) %>%
    dplyr::count(health_zone, name = "n_ll")

  recon <- sit_tot %>%
    dplyr::left_join(ll_tot, by = c("nom" = "health_zone")) %>%
    dplyr::mutate(n_ll = as.integer(dplyr::coalesce(n_ll, 0L)),
                  shortfall = pmax(0L, n_sit - n_ll))

  if (!is.null(spine)) {
    off <- setdiff(recon$nom[recon$shortfall > 0L], spine)
    if (length(off) > 0)
      warning("[load_linelist] ", length(off),
              " sitrep zone(s) with confirmed cases are not in the WorldPop spine — appends skipped: ",
              paste(off, collapse = ", "), call. = FALSE)
    recon <- recon %>% dplyr::filter(nom %in% spine)
  }
  recon <- recon %>% dplyr::filter(shortfall > 0L)
  if (nrow(recon) == 0) return(tibble::tibble())

  # For each zone take the LATEST `shortfall` dated cases (one date per unit increment).
  pick_dates <- function(z, short) {
    ev <- incs %>% dplyr::filter(nom == z) %>% dplyr::arrange(date)
    utils::tail(rep(ev$date, ev$inc), short)         # ascending; keep the most recent `short`
  }
  # PROVINCE for the appended rows. The sitrep cumulative file carries only
  # (nom, date, cumulative_confirmed_cases), so province has to come from elsewhere.
  # Leaving it NA is not harmless: every consumer that groups by (health_zone,
  # province) — Figure 1A's bivariate choropleth among them — then splits a topped-up
  # zone into TWO groups, one carrying the line-list cases and one the appended cases,
  # so the zone is counted twice and both halves are binned on a fraction of its burden.
  # Preference order: the zone's own line-list rows (always right when it has any),
  # then the shapefile attribute table for the zones that exist only in the sitrep.
  .prov_from_ll <- ll %>%
    dplyr::filter(!is.na(province), !is.na(health_zone)) %>%
    dplyr::count(health_zone, province, sort = TRUE) %>%
    dplyr::distinct(health_zone, .keep_all = TRUE)
  .prov_lookup <- stats::setNames(.prov_from_ll$province, .prov_from_ll$health_zone)
  .need_shp <- setdiff(recon$nom, names(.prov_lookup))
  if (length(.need_shp)) {
    .shp_prov <- tryCatch({
      stopifnot(requireNamespace("sf", quietly = TRUE), file.exists(SHAPEFILE_PATH))
      a <- sf::st_drop_geometry(sf::st_read(SHAPEFILE_PATH, quiet = TRUE))
      stopifnot(all(c("Nom", "PROVINCE") %in% names(a)))
      nm <- as.character(a$Nom)
      if (!is.null(aliases_all) && nrow(aliases_all) > 0) nm <- .apply_aliases(nm, aliases_all)
      stats::setNames(as.character(a$PROVINCE), nm)
    }, error = function(e) character(0))
    .hit <- intersect(.need_shp, names(.shp_prov))
    if (length(.hit)) .prov_lookup[.hit] <- unname(.shp_prov[.hit])
    .still <- setdiff(.need_shp, .hit)
    if (length(.still))
      warning("[load_linelist] province unresolved for ", length(.still),
              " sitrep-only zone(s); their appended rows carry province = NA and will not ",
              "group with any line-list rows: ", paste(.still, collapse = ", "), call. = FALSE)
  }

  mk <- function(z, short) {
    d <- pick_dates(z, short)
    tibble::tibble(
      country                       = "Democratic Republic of the Congo",
      province                      = unname(.prov_lookup[z]),
      health_zone                   = z,
      health_zone_unmatched         = FALSE,
      health_area_unmatched         = FALSE,
      date_of_sample_collection     = as.Date(d),
      date_of_symptom_onset         = as.Date(NA),   # imputed downstream from the delay dist
      final_mve_case_classification = CONFIRMED_STATUS,
      genexpert_result              = "positive",
      # samples_received / samples_analyzed are deliberately NA, NOT 1. These rows are a
      # COUNT RECONCILIATION against the sitrep cumulative, not laboratory records: no test
      # was observed for them. Setting them to 1 made them the ONLY rows in the whole line
      # list with a non-NA samples_analyzed (every real DHIS2 record is NA), which turned
      # the derived `positivity` covariate into "confirmed cases per synthetic row" —
      # values up to 74, non-NA for exactly the zones that received a top-up. See the
      # positivity guard in aggregate_to_zone_week() below.
      samples_received              = NA_real_,
      samples_analyzed              = NA_real_,
      alert_id = sprintf("SITREP-CONF-%s-%02d", gsub("[^A-Za-z0-9]", "", z), seq_along(d))
    )
  }
  appended <- purrr::map2_dfr(recon$nom, recon$shortfall, mk)

  message(sprintf(
    "[load_linelist] Sitrep reconciliation: appended %d confirmed case(s) across %d zone(s) (top-up to sitrep cumulative): %s",
    nrow(appended), nrow(recon),
    paste(sprintf("%s+%d", recon$nom, recon$shortfall), collapse = ", ")))
  appended
}


# =============================================================================
# SECTION 2: load_linelist()
# =============================================================================
#
# Reads the DHIS2 processed linelist identified by LINELIST_JSON (latest.json),
# filters to the current BDBV 2026 outbreak (>= OUTBREAK_START), handles
# zone-name aliasing, and excludes records with unmatched health zones.
#
# Returns a tibble with:
#   - All original columns (dates coerced to Date class)
#   - `date_index`    : date_of_symptom_onset if available, else date_of_sample_collection
#   - `confirmed`     : logical — final_mve_case_classification == CONFIRMED_STATUS
#   - `suspected`     : logical — final_mve_case_classification == SUSPECTED_STATUS
#   - `health_zone`   : alias-corrected zone name

load_linelist <- function() {

  # ---- resolve path via latest.json ----------------------------------------
  if (!file.exists(LINELIST_JSON)) {
    stop("[load_linelist] latest.json not found: ", LINELIST_JSON, call. = FALSE)
  }
  meta <- tryCatch(
    jsonlite::fromJSON(LINELIST_JSON, simplifyVector = TRUE),
    error = function(e) stop("[load_linelist] Cannot parse latest.json: ", e$message, call. = FALSE)
  )
  if (!"folder" %in% names(meta)) {
    stop("[load_linelist] latest.json has no 'folder' key.", call. = FALSE)
  }
  ll_csv <- file.path(LINELIST_DIR, meta$folder, "dhis2_processed_linelist.csv")
  if (!file.exists(ll_csv)) {
    stop("[load_linelist] Linelist CSV not found: ", ll_csv, call. = FALSE)
  }

  message("[load_linelist] Reading: ", ll_csv)
  ll <- tryCatch(
    readr::read_csv(
      ll_csv,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE,
      na = c("", "NA", "N/A")
    ),
    error = function(e) stop("[load_linelist] Read error: ", e$message, call. = FALSE)
  )
  message("[load_linelist] Rows loaded: ", nrow(ll))

  # ---- check minimum required columns ---------------------------------------
  required_cols <- c(
    "health_zone", "health_zone_unmatched",
    "date_of_sample_collection", "date_of_symptom_onset",
    "final_mve_case_classification"
  )
  .check_cols(ll, required_cols, "load_linelist")

  # ---- coerce types ----------------------------------------------------------
  # Logical columns
  ll <- ll %>%
    dplyr::mutate(
      health_zone_unmatched  = .parse_logical_col(health_zone_unmatched),
      health_area_unmatched  =
        if ("health_area_unmatched" %in% names(.))
          .parse_logical_col(health_area_unmatched)
        else
          NA,
      permanent_residence_health_area_unmatched =
        if ("permanent_residence_health_area_unmatched" %in% names(.))
          .parse_logical_col(permanent_residence_health_area_unmatched)
        else
          NA
    )

  # Date columns: parse all date-like columns to Date
  date_cols <- c(
    "reporting_date", "date_of_sample_collection", "date_of_symptom_onset",
    "investigation_datetime", "date_of_hospitalisation_start",
    "date_of_hospitalisation_end", "date_of_hospitalisation_end_2",
    "date_of_notification", "lab_analysis_date", "contact_end_date"
  )
  for (dc in intersect(date_cols, names(ll))) {
    ll[[dc]] <- .parse_date(ll[[dc]])
  }

  # Numeric columns
  if ("samples_analyzed" %in% names(ll)) {
    ll$samples_analyzed <- suppressWarnings(as.numeric(ll$samples_analyzed))
  }
  if ("samples_received" %in% names(ll)) {
    ll$samples_received <- suppressWarnings(as.numeric(ll$samples_received))
  }

  # ---- filter to outbreak start date ----------------------------------------
  # Keep a record if EITHER its symptom onset OR its sample-collection date is on/after
  # OUTBREAK_START. The forecast target is onset-defined, so requiring a NON-NA sample
  # date silently dropped onset-only confirmed cases (e.g. epicentre cases with a blank
  # sample date). Records with both dates pre-outbreak, or neither present, are dropped.
  n_before <- nrow(ll)
  ll <- ll %>%
    dplyr::filter(
      (!is.na(date_of_symptom_onset)     & date_of_symptom_onset     >= OUTBREAK_START) |
      (!is.na(date_of_sample_collection) & date_of_sample_collection >= OUTBREAK_START))
  message(
    "[load_linelist] After OUTBREAK_START filter (",
    OUTBREAK_START, "): ", nrow(ll), " rows (dropped ", n_before - nrow(ll), ")"
  )

  # ---- exclude unmatched health zones ----------------------------------------
  n_unmatched <- sum(ll$health_zone_unmatched %in% TRUE, na.rm = TRUE)
  if (n_unmatched > 0) {
    warning(
      "[load_linelist] Excluding ", n_unmatched,
      " records where health_zone_unmatched == TRUE",
      call. = FALSE
    )
    ll <- ll %>% dplyr::filter(!(health_zone_unmatched %in% TRUE))
  }

  # ---- apply zone name aliases -----------------------------------------------
  if (!file.exists(ALIASES_PATH)) {
    warning("[load_linelist] aliases.csv not found — skipping alias correction.", call. = FALSE)
  } else {
    aliases <- tryCatch(
      readr::read_csv(ALIASES_PATH, col_types = readr::cols(.default = readr::col_character()),
                      show_col_types = FALSE),
      error = function(e) {
        warning("[load_linelist] Cannot read aliases.csv: ", e$message, call. = FALSE)
        NULL
      }
    )
    if (!is.null(aliases)) {
      .check_cols(aliases, c("observed_name", "canonical_nom"), "load_linelist/aliases")
      # Only apply aliases relevant to the linelist source. `source_dataset` is optional
      # (not required by .check_cols), so guard its presence — a future aliases.csv lacking
      # the column would otherwise crash load_linelist() here rather than degrade.
      # Apply linelist/epi/dhis2 aliases AND the shapefile-migration aliases (old canonical zone
      # spellings -> the current shapefile/spine names, e.g. Nyakunde -> Nyankunde). The line list
      # still carries the OLD spellings, so without the migration entries those already-affected
      # zones fail to match the spine, are dropped from the zone-week grid, and then reappear as
      # spurious high-risk "at-risk" zones. Match "shapefile*" loosely so a future dated migration
      # (shapefile_migration_YYYY-MM) is picked up automatically.
      ll_aliases <- if ("source_dataset" %in% names(aliases)) {
        aliases %>%
          # "insp_sitrep" is INCLUDED (2026-09-17). The sitrep reconciliation
          # (.build_sitrep_confirmed_appends) canonicalises its side with the FULL alias table,
          # so mappings tagged insp_sitrep (Lubunga -> Lubunga (Tshopo), Rumba -> Rimba,
          # Tchomai -> Tchomia, Manguripa -> Manguredjipa) were applied to one side only. If a
          # DHIS2 export ever used one of those spellings, the join would report n_ll = 0 for the
          # canonical zone and append the ENTIRE sitrep cumulative on top of the real cases — a
          # straight double count, breaking the "line list UNION sitrep floor" invariant in the
          # one direction it must never break. These mappings are no-ops on the current line list.
          dplyr::filter(source_dataset %in% c("linelist", "epi", "dhis2", "insp_sitrep") |
                        grepl("shapefile", source_dataset, ignore.case = TRUE) |
                        is.na(source_dataset))
      } else aliases
      if (nrow(ll_aliases) > 0) {
        old <- ll$health_zone
        ll$health_zone <- .apply_aliases(ll$health_zone, ll_aliases)
        n_aliased <- sum(old != ll$health_zone, na.rm = TRUE)
        if (n_aliased > 0) {
          message("[load_linelist] Zone aliases applied: ", n_aliased, " records renamed.")
        }
      }
    }
  }

  # ---- append sitrep-confirmed top-up rows (reconcile line list to sitrep) ----
  # Injected AFTER the outbreak filter, unmatched exclusion, and aliasing (so the
  # appends carry canonical, in-spine, in-outbreak, matched zones) and BEFORE the
  # onset-imputation block below (so their blank onset is imputed from the fitted
  # onset->sample delay, like any real onset-less record). The appends carry no
  # complete onset+sample pair, so they do NOT perturb the delay estimation; and
  # binding them at the END preserves the positional RNG draws of the existing rows.
  # `aliases` is assigned only inside the alias block above (when aliases.csv exists); guard
  # with exists() so a missing aliases.csv degrades to NULL (the helper then skips sitrep-name
  # canonicalisation) instead of raising "object 'aliases' not found".
  .sitrep_appends <- .build_sitrep_confirmed_appends(
    ll, if (exists("aliases", inherits = FALSE)) aliases else NULL)
  if (nrow(.sitrep_appends) > 0) ll <- dplyr::bind_rows(ll, .sitrep_appends)

  # ---- construct date_index and classification flags -------------------------
  # The forecast target is onset-dated (P of first symptom ONSET). Records missing an onset
  # date are imputed to onset = sample_date - delay, where the delay is DRAWN from the fitted
  # onset->sample distribution (see below) rather than fixed at its mean — so imputed onsets
  # contribute a true onset week spread by the delay's variance. Records with neither date are
  # dropped downstream (date_index = NA).
  # Onset->sample reporting delay used for the single (stochastic, per-record draw) imputation of
  # missing onsets. Default: the fixed lab-linelist Exp(rate) from config. Set
  # ONSET_SAMPLE_DELAY_SOURCE = "data" (or
  # IMPUTE_DELAY_FROM_DATA = TRUE) to instead ESTIMATE the mean delay from the CURRENT
  # linelist's complete onset+sample pairs (MLE of an Exponential = 1/mean-delay) — the right
  # choice when processing DHIS2 data whose reporting delay differs from the lab linelist.
  .delay_src <- get0("ONSET_SAMPLE_DELAY_SOURCE",
                     ifnotfound = if (isTRUE(get0("IMPUTE_DELAY_FROM_DATA", ifnotfound = FALSE))) "data" else "lab")
  .fixed_rate <- get0("DELAY_ONSET_SAMPLE_RATE", ifnotfound = NA_real_)
  # OUTBREAK_START itself, NOT its WEEK_ANCHOR-floored week. WEEK_ANCHOR is derived from the
  # weekday ANALYSIS_DATE falls on (00_config.R), so flooring here made the onset-plausibility
  # floor slide with the day of the week the pipeline was run: on the same frozen snapshot a
  # Monday run floored to 2026-04-28 and a Thursday run to 2026-04-24, flipping onsets in that
  # band between "usable" and "imputed". The grid re-anchoring is deliberate and documented; this
  # side effect on record-level usability was not. A plausibility floor is a property of the
  # outbreak, not of the run.
  .ob_floor <- OUTBREAK_START
  # WINDOWED complete onset+sample pairs for the empirical delay bootstrap: onset in
  # [outbreak-week floor, max(sample) - TEST_DAYS]. Windowing (mirroring the delay estimator
  # 04c_dhis2_delay_windows.R) drops the right-truncated final days — where recent onsets that
  # WILL be sampled at long delays are not yet observed, biasing the raw complete-pair delay
  # SHORT — and any pre-outbreak onset typos, so the bootstrap sample is the clean,
  # representative onset->sample delay. If the window is too thin (<30 pairs) it falls back to
  # ALL complete pairs, so a short line list is never left with no delay data.
  .test_days <- get0("DELAY_TRUNC_BUFFER_DAYS", ifnotfound = 5L)
  # Shared onset->sample plausibility ceiling (00_config.R) — the SAME bound the delay estimator
  # fits under (04c MAX_DELAY), so the bootstrap pairs, the imputed-delay clamp and the fit all
  # share one support instead of the previous 90-vs-60 split.
  .max_plausible <- get0("DELAY_MAX_PLAUSIBLE_DAYS", ifnotfound = 60L)
  # Clerical tolerance for an onset recorded AFTER its own specimen (00_config.R). Named
  # rather than written inline at its single use site: its value and its meaning have to be
  # read together, and an inline "+ 2L" invited being read as a biological claim.
  .neg_tol <- get0("ONSET_SAMPLE_NEG_TOL_DAYS", ifnotfound = 2L)
  .max_samp  <- suppressWarnings(max(ll$date_of_sample_collection, na.rm = TRUE))
  # As-of consistency (mirror 04c): never let a sample collected AFTER the analysis date define
  # the window edge — a no-op when the line list is snapshotted to the as-of date, but on a
  # back-dated re-run (or a future-dated sample-date typo) it stops future data leaking in.
  .asof_date <- suppressWarnings(as.Date(get0("ANALYSIS_DATE", ifnotfound = NA)))
  if (length(.asof_date) == 1L && !is.na(.asof_date) && is.finite(.max_samp))
    .max_samp <- min(.max_samp, .asof_date)
  .trunc_d   <- if (is.finite(.max_samp)) .max_samp - .test_days else as.Date(NA)
  .dd_raw    <- as.numeric(ll$date_of_sample_collection - ll$date_of_symptom_onset)
  .complete  <- !is.na(ll$date_of_symptom_onset) & !is.na(ll$date_of_sample_collection) &
                is.finite(.dd_raw) & .dd_raw >= 0 & .dd_raw <= .max_plausible &   # plausible complete pairs
                (is.na(.asof_date) | ll$date_of_sample_collection <= .asof_date)  # as-of observable only
  .in_win    <- .complete & ll$date_of_symptom_onset >= .ob_floor &
                (is.na(.trunc_d) | ll$date_of_symptom_onset <= .trunc_d)
  .windowed  <- sum(.in_win) >= 30L
  .dd        <- if (.windowed) .dd_raw[.in_win] else .dd_raw[.complete]     # windowed if enough, else all
  # (The companion `.dd_onset` vector was REMOVED. It existed for a Lynden-Bell right-truncation
  # reweighting of the empirical pool, and its comment still described those weights as being
  # applied "below" — but that machinery was deleted when the imputation moved to drawing from
  # the EpiDist resolver, which corrects truncation in the fit itself. The variable was computed
  # and never read, and it pointed an auditor at code that no longer exists.)
  # Rigorous DHIS2 delay params (windowed interval-censored MLE) from
  # 04c_dhis2_delay_windows.R, when it has been run and the source is "data": its AIC-best
  # family + Exponential-rate summary supersede the crude 1/mean for the reported rate and the
  # parametric-fallback draw (a gamma/weibull/lnorm fit rather than an Exponential).
  .dhis2_delay <- if (identical(.delay_src, "data")) .load_dhis2_delay_params() else NULL
  # The REPORTED rate follows the same resolver the imputation draws from, so the logged rate and
  # the realised draw cannot describe different distributions. 1/mean(.dd) is the right-truncated
  # empirical rate and is used only when no fitted delay exists at all — the same condition under
  # which the draw itself falls back, and it warns then.
  .est_rate <- local({
    # THE CONFIRMED-ONLY DELAY. This block imputes onsets for CONFIRMED records, so it must
    # draw from the confirmed-case delay. Until 2026-09-22 it drew from a fit pooled over
    # every classification, and on this line list that pool is majority TEST-NEGATIVE --
    # 6,823 not_a_case windowed pairs and 2,038 with NO final classification, against 4,981
    # confirmed and 145 suspected. Both of those groups are swabbed faster (raw means 6.40 d
    # and 3.56 d against 8.80 d), so the pooled EpiDist marginal returns 7.67 d
    # where the CONFIRMED stratum gives 10.14 d -- a 2.47 d (32%) gap, and imputed onsets
    # landed that much LATE for the ~23% of confirmed records that carry one. (Both numbers
    # are the truncation-corrected EpiDist marginal, which is what the resolver returns; the
    # windowed interval-censored MLE puts the same contrast at 7.67 vs 9.11 d, so quote the
    # estimator with the number.) Because the invasion outcome is the FIRST
    # onset in a zone, a shift of that size can move a zone's invasion week -- this touches
    # the outcome, not merely a covariate.
    # A missing stratum warns inside the resolver and falls back to the pooled fit.
    r <- tryCatch(effective_onset_sample_delay(stratum = "confirmed"), error = function(e) NULL)
    if (!is.null(r) && identical(r$source, "data") && is.finite(r$rate) && r$rate > 0) r$rate
    else if (!is.null(.dhis2_delay) && is.finite(.dhis2_delay$rate)) .dhis2_delay$rate
    else if (length(.dd) >= 30 && mean(.dd) > 0) 1 / mean(.dd) else NA_real_
  })
  .rate_used <- if (identical(.delay_src, "data") && is.finite(.est_rate)) .est_rate else .fixed_rate
  # --- Onset-handling MODE (review §1.1/§1.3) --------------------------------------------------
  # ONSET_MODE selects how confirmed cases lacking a usable onset are handled (see 00_config.R):
  # impute (default) / growth_impute (growth-tilted backward draw) / complete_case (drop) /
  # sample_verbatim (date at sample). Unset -> derived from IMPUTE_ONSET_FROM_SAMPLE.
  .onset_mode <- get0("ONSET_MODE", ifnotfound = NA_character_)
  if (length(.onset_mode) != 1L || is.na(.onset_mode))
    .onset_mode <- if (isTRUE(get0("IMPUTE_ONSET_FROM_SAMPLE", ifnotfound = TRUE))) "impute" else "sample_verbatim"
  if (!.onset_mode %in% c("impute", "growth_impute", "complete_case", "sample_verbatim")) {
    warning("[load_linelist] unknown ONSET_MODE '", .onset_mode, "'; using 'impute'"); .onset_mode <- "impute" }
  # .imp_active is resolved HERE, after .onset_mode, and from .onset_mode alone.
  #   (a) 00_config.R states ONSET_MODE "supersedes IMPUTE_ONSET_FROM_SAMPLE". It did not: this
  #       gated on IMPUTE_ONSET_FROM_SAMPLE as well, so ONSET_MODE="impute" with
  #       IMPUTE_ONSET_FROM_SAMPLE=FALSE silently degraded to .imp_mode="none" (onset = sample
  #       verbatim) — the documented override was inoperative. The legacy flag still selects the
  #       DEFAULT mode above when ONSET_MODE is unset, which is the compatibility that matters.
  #   (b) The finite-rate precondition was a PARAMETRIC-branch requirement applied to all
  #       branches. The empirical bootstrap needs no rate at all (13,984 windowed pairs are available on
  #       the current snapshot), so a missing rate needlessly disabled the better estimator.
  #       The rate is now required only where it is actually used.
  # complete_case is included so .imp_mode resolves to "complete_case" rather than falling
  # through to "none" (which emitted a false "no delay available" warning and made the
  # complete_case switch arm unreachable). It needs no delay, hence the || TRUE arm.
  .imp_active <- (.onset_mode == "complete_case") ||
                 (.onset_mode %in% c("impute", "growth_impute") &&
                  (length(.dd) >= 30L || (is.finite(.rate_used) && .rate_used > 0)))
  # National per-day epidemic growth rate r: log-linear slope of recent weekly confirmed counts
  # (by sample date). Used ONLY to growth-tilt the backward delay draw under "growth_impute":
  # weighting Delta by exp(-r*Delta) down-weights long delays while the epidemic grows, so imputed
  # onsets are not pushed systematically too early (epidemic processes are not time-reversible).
  .estimate_growth_rate <- function(dates, asof,
                                    window_weeks = get0("ONSET_GROWTH_WINDOW_WEEKS", ifnotfound = 8L),
                                    trunc_buffer = get0("ONSET_GROWTH_TRUNC_BUFFER_DAYS", ifnotfound = 14L)) {
    d <- dates[!is.na(dates) & (is.na(asof) | dates <= asof)]
    if (length(d) < 20L) return(NA_real_)

    # END OF THE FITTING WINDOW. These are SAMPLE-collection dates and the most recent are
    # right-truncated: a case sampled three days ago may not be in the extract yet. Fitting
    # through that tail reads the reporting lag as an epidemiological decline.
    end <- if (!is.na(asof)) asof else max(d)
    if (is.finite(trunc_buffer) && trunc_buffer > 0) end <- end - trunc_buffer
    d <- d[d <= end]
    if (length(d) < 20L) return(NA_real_)

    # BIN BACKWARD FROM `end`, SO EVERY BIN IS A FULL 7 DAYS. Binning forward from min(d) —
    # floor((d - min(d)) / 7) — leaves the FINAL bin partial: it holds only the days between
    # the last 7-day boundary and the cutoff. That bin is undercounted by construction, and a
    # log-linear fit reads it as a collapse. Measured on the 2026-09-07 frame the final forward
    # bin held 23 cases against ~1,800 in the preceding full week, and the fitted rate came out
    # at -0.033/day — a halving time of three weeks — for an epidemic that is in fact GROWING
    # at about +0.006/day. The artefact survived every truncation buffer, because shifting the
    # cutoff just moves where the partial bin falls.
    #
    # THE SIGN IS WHAT MATTERS. This rate tilts the backward onset draw by exp(-r*Delta). With
    # the true r > 0 the tilt shortens imputed delays (onsets later); a spuriously negative r
    # LENGTHENS them, actively worsening the early-shift bias the tilt exists to remove. An
    # estimator that can flip the sign of the correction is worse than no correction.
    back <- as.integer(floor(as.numeric(end - d) / 7))   # 0 = the most recent COMPLETE week
    # Keep only bins wholly inside the observed span (the oldest bin can otherwise be partial
    # for the opposite reason) and inside the requested window.
    oldest_ok <- as.integer(floor(as.numeric(end - min(d)) / 7))
    keep_bins <- back < window_weeks & back < oldest_ok
    back <- back[keep_bins]
    if (!length(back)) return(NA_real_)
    tab <- table(back)
    bk  <- as.integer(names(tab)); cnt <- as.numeric(tab)
    if (length(bk) < 3L || sum(cnt) < 10) return(NA_real_)
    fit <- tryCatch(stats::lm(log(cnt + 0.5) ~ bk), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    # `bk` counts WEEKS INTO THE PAST, so a growing epidemic has a NEGATIVE slope against it.
    # Forward per-day growth rate is therefore -slope / 7.
    r <- -unname(coef(fit)[2]) / 7
    if (!is.finite(r)) return(NA_real_)

    # PLAUSIBILITY BOUND. The tilt weight is exp(-r*Delta) over Delta up to
    # DELAY_MAX_PLAUSIBLE_DAYS; a wild |r| makes it collapse onto one end of the delay support
    # and the importance resample degenerates to a handful of distinct values. |r| <= 0.1/day
    # is a weekly growth factor of ~2 — far outside anything this outbreak shows.
    .rmax <- get0("ONSET_GROWTH_RATE_MAX", ifnotfound = 0.1)
    if (abs(r) > .rmax) {
      warning(sprintf(paste0("[load_linelist] estimated growth rate %+0.4f/day exceeds the ",
                             "plausibility bound %.2f/day; the growth tilt is DISABLED for this ",
                             "run and imputation falls back to the untilted draw."), r, .rmax),
              call. = FALSE)
      return(NA_real_)
    }
    r
  }

  .r_growth <- if (identical(.onset_mode, "growth_impute"))
    .estimate_growth_rate(ll$date_of_sample_collection, .asof_date) else NA_real_
  .dd_wts <- if (identical(.onset_mode, "growth_impute") && is.finite(.r_growth) && length(.dd))
    exp(-.r_growth * .dd) else NULL
  if (identical(.onset_mode, "growth_impute"))
    message(sprintf("[load_linelist] growth_impute: national growth rate r=%.4f/day; backward delay draw tilted by exp(-r*Delta)",
                    ifelse(is.finite(.r_growth), .r_growth, NA_real_)))
  # STOCHASTIC single imputation: draw EACH missing onset's delay from the onset->sample delay
  # distribution and subtract it from the sample date, rather than subtracting the single MEAN
  # (which piled every imputed case on one week — an artificial spike). We draw from the
  # EMPIRICAL delay distribution (bootstrap the WINDOWED complete onset+sample pairs) when >=30
  # exist: the most faithful "correct delay distribution" and, unlike a mean-matched Exponential
  # (whose mode is 0), it reproduces the true delay SHAPE. Falls back to the parametric fit
  # otherwise — the rigorous DHIS2 best-family fit (04c_dhis2_delay_windows.R) if present, else
  # Exp(rate_used). Delays clamped to [0, DELAY_MAX_PLAUSIBLE_DAYS] d (60, not 90 — the
  # ceiling is the shared constant, and it is the same bound the delay fit treats as an
  # outlier). Seeded for reproducibility.
  # CAVEATS (documented, not hidden): a single stochastic imputation does not fully propagate
  # imputation uncertainty (full multiple imputation would loop the pipeline and pool), so
  # imputation-dependent intervals are conditionally slightly narrow. Windowing to onset <=
  # max(sample) - TEST_DAYS mitigates (but does not fully remove) the right-truncation of the
  # empirical delay; the EpiDist marginal model in 04c_dhis2_delay_windows.R corrects it fully.
  # Reproducible draw WITHOUT clobbering the caller's global RNG stream: snapshot .Random.seed,
  # seed locally, draw, then restore. (Draws are positional — draw i -> row i — so reproducibility
  # also assumes a stable linelist row order, which the JSON read + dplyr preserve.)
  .rng_state <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(.rng_state)) assign(".Random.seed", .rng_state, envir = .GlobalEnv), add = TRUE)
  set.seed(get0("RANDOM_SEED", ifnotfound = 20260704L))
  # ONSET IMPUTATION DRAWS FROM THE SHARED DELAY RESOLVER (2026-09-17, revised).
  #
  # 00_config.R states that four consumers read ONE delay resolver so the estimator "can never
  # drift apart between consumers again", and names this imputation as consumer 1. It was not:
  # whenever >=30 complete pairs existed the draw bootstrapped the RAW windowed pairs, whose
  # mean is 6.81 d, while the nowcast and the EpiNow2 right-truncation model were simultaneously
  # using effective_onset_sample_delay() at 7.67 d. ~3,600 imputed onsets sat ~0.85 d too late.
  #
  # The resolver IS the truncation-corrected estimator: on this snapshot it is the EpiDist
  # MARGINAL fit, gamma(shape 0.812, rate 0.106), mean 7.667 d, SD 8.507 d, corrected for BOTH
  # right truncation and double interval censoring (04c_dhis2_delay_windows.R). Drawing from it
  # removes the bias exactly rather than approximately, and makes the single-resolver invariant
  # true instead of aspirational.
  #
  # An earlier revision kept the empirical bootstrap and applied Lynden-Bell style weights
  # 1/G(T-d). That is the right form for right-truncated data, but G must be estimated from the
  # observed onsets, which are themselves truncated — so the plug-in under-corrects (it reached
  # 7.48 d against 7.67 d, ~81% of the gap). The exact fix is the joint NPMLE; the resolver
  # already IS a truncation-corrected fit, so the bootstrap is no longer the better estimator.
  # The empirical bootstrap is retained only as the fallback when no fitted delay exists.
  .resolver <- tryCatch(effective_onset_sample_delay(), error = function(e) NULL)
  .resolver_ok <- !is.null(.resolver) && identical(.resolver$source, "data") &&
                  is.finite(.resolver$mean) && .resolver$mean > 0
  # ONSET_SAMPLE_DELAY_SOURCE = "lab" is a DELIBERATE choice, documented in 00_config.R as
  # "deliberately imputing DHIS2 onsets with the (faster) lab delay". It did not do that: the
  # resolver returns source = "lab", .resolver_ok was FALSE, and the cascade below fell through
  # to the EMPIRICAL bootstrap of the DHIS2 pairs (mean ~6.8 d) — so a sensitivity run under
  # "lab" measured the DHIS2 delay, not the lab delay (mean 4.39 d) it asked for, and the
  # warning below told the user there was "no fitted delay on disk", which was untrue.
  # Honour the explicit request; an ACCIDENTAL lab fallback (no params file) still is not
  # treated as a resolver, because that case must warn rather than silently substitute.
  .lab_requested <- !identical(get0("ONSET_SAMPLE_DELAY_SOURCE", ifnotfound = "data"), "data")
  .resolver_lab_ok <- .lab_requested && !is.null(.resolver) &&
                      identical(.resolver$source, "lab") &&
                      is.finite(.resolver$mean) && .resolver$mean > 0

  .imp_mode <- if (.onset_mode == "sample_verbatim" || !.imp_active) "none"
               else if (.onset_mode == "complete_case") "complete_case"
               else if ((.resolver_ok || .resolver_lab_ok) && .onset_mode == "growth_impute")
                 "resolver_growth"
               else if (.resolver_ok || .resolver_lab_ok) "resolver"
               else if (.onset_mode == "growth_impute" && length(.dd) >= 30L) "growth"
               else if (length(.dd) >= 30L) "empirical"
               else "parametric"

  # Parametric fallback draw (no resolver AND <30 complete pairs): Exp(rate_used).
  .param_draw <- function(n) stats::rexp(n, rate = .rate_used)

  # growth_impute tilts the BACKWARD draw by exp(-r*Delta): while incidence grows, a case
  # observed now is likelier to have a short delay, so long delays must be down-weighted.
  # With a parametric resolver there is no pool to reweight, so draw a large pool FROM the
  # resolver and importance-resample it with those weights — the tilt is applied to the
  # corrected distribution rather than to the truncated empirical one.
  .resolver_draw <- function(n, tilt = FALSE) {
    if (!tilt || !is.finite(.r_growth)) return(.draw_dhis2_delay(.resolver, n))
    pool <- .draw_dhis2_delay(.resolver, max(20000L, 10L * n))
    pool <- pool[is.finite(pool) & pool >= 0]
    if (!length(pool)) return(.draw_dhis2_delay(.resolver, n))
    w <- exp(-.r_growth * pool)
    if (!any(is.finite(w)) || sum(w[is.finite(w)]) <= 0) return(pool[sample.int(length(pool), n, TRUE)])
    w[!is.finite(w)] <- 0
    pool[sample.int(length(pool), n, replace = TRUE, prob = w)]
  }

  ll$.imp_delay <- switch(.imp_mode,
    resolver        = as.integer(pmin(pmax(round(.resolver_draw(nrow(ll))), 0L), .max_plausible)),
    resolver_growth = as.integer(pmin(pmax(round(.resolver_draw(nrow(ll), tilt = TRUE)), 0L), .max_plausible)),
    # FALLBACKS ONLY (no fitted delay on disk). These remain right-truncated; the message below
    # says so rather than letting a degraded estimator pass as the corrected one.
    empirical  = as.integer(pmin(pmax(round(sample(.dd, nrow(ll), replace = TRUE)), 0L), .max_plausible)),
    growth     = as.integer(pmin(pmax(round(sample(.dd, nrow(ll), replace = TRUE,
                                              prob = .dd_wts)), 0L), .max_plausible)),
    parametric = as.integer(pmin(pmax(round(.param_draw(nrow(ll))), 0L), .max_plausible)),
    # complete_case: onset-less confirmed records are DROPPED in the mutate below (delay unused).
    complete_case = rep(0L, nrow(ll)),
    none       = rep(0L, nrow(ll)))   # no delay info at all: degenerate (onset=sample); warned below
  if (.imp_mode %in% c("empirical", "growth"))
    warning(sprintf(paste0("[load_linelist] no usable fitted onset->sample delay was resolved%s; imputation ",
                           "fell back to the RIGHT-TRUNCATED empirical pool (mean %.2f d). It is NOT consistent ",
                           "with the nowcast or the R(t) truncation model. Run 04c_dhis2_delay_windows.R."),
                   # Do not assert "no fit on disk" without checking: the fit may exist and the
                   # resolver may simply have been overridden, which is a different problem with
                   # a different fix, and the old wording sent the user to the wrong one.
                   if (file.exists(get0("DELAY_PARAMS_PATH", ifnotfound = "")))
                     " (a params file EXISTS on disk - check ONSET_SAMPLE_DELAY_SOURCE)" else
                     " (no params file on disk)",
                   mean(.dd)),
            call. = FALSE)
  if (identical(.imp_mode, "none"))
    warning("[load_linelist] Onset imputation requested but no onset->sample delay is available ",
            "(no fixed rate and <30 complete pairs); imputed onsets fall back to the sample date.")
  ll <- ll %>%
    dplyr::mutate(
      # A recorded onset is USABLE only if it is epidemiologically plausible: not before
      # the outbreak week (data-entry YEAR TYPOS put onsets years early — onset 2020/2023
      # with a 2026 sample), not more than 2 DAYS after its own sample, and not implausibly
      # long before it (> the shared DELAY_MAX_PLAUSIBLE_DAYS ceiling — the same bound the
      # delay fit treats as an outlier). An unusable onset is treated as MISSING and imputed
      # from the sample date, so a typo can neither leak a spurious pre-outbreak week into the
      # zone-week grid (the CRITICAL corruption) nor drop the real case.
      #
      # ONSET RECORDED AFTER ITS OWN SPECIMEN: KEPT, BUT CENSORED AT THE SPECIMEN DATE.
      # `.neg_tol` (ONSET_SAMPLE_NEG_TOL_DAYS) is how far past the specimen an onset may fall
      # and still be treated as recording noise rather than a lost field. Three decisions are
      # bundled here and each has a different reason.
      #
      # WHY THE RECORD IS NOT REJECTED. An unusable onset is not dropped, it is imputed from
      # the specimen date at a mean delay of ~7-8 d. Rejecting a -1 d record therefore moves
      # it about 8 d EARLIER — a larger error, and in the wrong direction, than the 1-2 d of
      # noise it removes. This was the previous rationale for the tolerance and it still holds.
      #
      # WHY THE ONSET IS NOT CARRIED FORWARD EITHER (this is the change). These are not
      # presymptomatic detections of traced contacts, which is the only reading under which a
      # post-specimen onset is real. On the 2026-09-07 snapshot 39.5% of the -2/-1 d records
      # were deceased at swab against 22.4% of the positive-delay records, and half of the
      # affected CONFIRMED rows are death alerts — a person cannot develop symptoms after
      # being swabbed post mortem. The bin also matches the delay-0 bin on that split (41.4%),
      # i.e. it is the same-day population displaced by a day or two of transcription noise.
      # Carrying the onset forward left `date_index` after the specimen date, which is why
      # these rows had to be held out of the delay-fitting pools (.complete requires a
      # non-negative delay) and the two rules disagreed by design.
      #
      # WHY THE WINDOW IS SMALL. The empirical delay distribution decays steeply out of zero
      # (1756 at 0 d, 92 at -1, 22 at -2, 7 at -3) and then runs flat and sparse to -177.
      # Inside the window censoring costs a median of 1 day; outside it the onset field is
      # not recoverable by censoring and imputation is the honest default.
      #
      # MEASURED WHEN THIS WAS INTRODUCED: no zone's FIRST confirmed onset is contributed by
      # a negative-delay record (the nearest is 4 days after its zone's first case), so no
      # invasion label moves under any choice of window.
      onset_usable = !is.na(date_of_symptom_onset) &
                     date_of_symptom_onset >= .ob_floor &
                     (is.na(date_of_sample_collection) |
                      (date_of_symptom_onset <= date_of_sample_collection + .neg_tol &
                       as.numeric(date_of_sample_collection - date_of_symptom_onset) <= .max_plausible)),
      onset_imputed = !onset_usable & !is.na(date_of_sample_collection),
      onset_censored = onset_usable & !is.na(date_of_sample_collection) &
                       date_of_symptom_onset > date_of_sample_collection,
      date_index = dplyr::if_else(
        onset_usable,
        # Censored at the specimen date where the two disagree in the impossible direction.
        # coalesce() covers a usable onset with NO specimen (pmin would return NA there), and
        # pmax holds the invariant the zone-week grid depends on: date_index never precedes
        # the outbreak floor, which a specimen dated before the floor could otherwise break.
        pmax(dplyr::coalesce(pmin(date_of_symptom_onset, date_of_sample_collection),
                             date_of_symptom_onset), .ob_floor),
        # For onset-less records: complete_case DROPS them (NA date_index -> dropped downstream,
        # review §1.3); otherwise the imputed onset = sample - a delay DRAWN from the (optionally
        # growth-tilted) fitted distribution, clamped to not precede the outbreak week
        # (sample_verbatim uses delay 0, i.e. onset = sample).
        if (identical(.onset_mode, "complete_case")) as.Date(NA)
        else pmax(date_of_sample_collection - .imp_delay, .ob_floor)
      ),
      confirmed = (final_mve_case_classification == CONFIRMED_STATUS) %in% TRUE,
      suspected = (final_mve_case_classification == SUSPECTED_STATUS) %in% TRUE
    )
  # AS-OF UPPER BOUND. Everything else in this loader guards the as-of date (the sitrep helper,
  # the delay window, the complete-pair mask), but date_index did not: it was bounded below by
  # OUTBREAK_START and above by nothing. A record with onset AND sample both mistyped into the
  # future passes onset_usable and lands in a future week. Censor rather than drop, so the record
  # is still counted if it is merely future-DATED but otherwise valid; NA date_index is dropped
  # downstream exactly as it is for records with no usable date at all.
  .asof_idx <- suppressWarnings(as.Date(get0("ANALYSIS_DATE", ifnotfound = NA)))
  if (length(.asof_idx) == 1L && !is.na(.asof_idx)) {
    .n_future <- sum(!is.na(ll$date_index) & ll$date_index > .asof_idx, na.rm = TRUE)
    if (.n_future > 0L) {
      warning(sprintf("[load_linelist] %d record(s) have date_index AFTER the analysis date (%s); their date_index is set to NA (dropped downstream). Check for future-dated onset/sample entries.",
                      .n_future, format(.asof_idx)), call. = FALSE)
      ll$date_index[!is.na(ll$date_index) & ll$date_index > .asof_idx] <- as.Date(NA)
    }
  }
  .imp_draws <- ll$.imp_delay[ll$onset_imputed %in% TRUE]
  ll <- ll %>% dplyr::select(-dplyr::any_of(".imp_delay"))
  if (identical(.onset_mode, "complete_case")) {
    .n_cc <- sum(is.na(ll$date_index) & !is.na(ll$date_of_sample_collection) &
                 (ll$confirmed %in% TRUE), na.rm = TRUE)
    message(sprintf("[load_linelist] complete-case (§1.3): dropped %d confirmed records lacking a usable onset (non-imputed analysis)", .n_cc))
  }
  # Name the estimator actually in force. This said "interval-censored" unconditionally, but
  # .dhis2_delay may be the EpiDist MARGINAL fit (estimator == "epidist_marginal"), which is
  # what the current params CSV carries — so the log and any methods text copied from it
  # mislabelled the estimator.
  # No %||% here. 01_data_prep.R runs at run_all.R step 1, BEFORE 03_mobility_matrices.R
  # defines the suite's copy, and 00_config.R deliberately does not rely on base R's (4.4+).
  # A local scalar-or-default is also stricter than %||%, which passes NA and length-0 through
  # — and switch() on NA_character_ is an error, not a fallthrough.
  .chr1 <- function(x, default) {
    if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(as.character(x))) default
    else as.character(x)
  }
  .param_desc <- if (!is.null(.dhis2_delay)) {
    .est <- .chr1(.dhis2_delay$estimator, "censored_mle")
    sprintf("the DHIS2 %s %s fit (rate %.3f/d, window %s)",
            switch(.est,
                   epidist_marginal = "EpiDist marginal (truncation-corrected)",
                   censored_mle     = "interval-censored MLE",
                   .est),
            .chr1(.dhis2_delay$family, "unknown-family"), .rate_used,
            .chr1(.dhis2_delay$window, "n/a"))
  } else sprintf("Exp(rate %.3f/d, source %s)", .rate_used, .delay_src)
  # Publish the realised imputed share so downstream prose derives it instead of hard-coding a
  # literal (17_invasion_viz.R's Weaknesses list said "~15%" where the truth is ~24%).
  .pct_imp <- 100 * sum(ll$onset_imputed %in% TRUE & ll$confirmed %in% TRUE, na.rm = TRUE) /
              max(sum(ll$confirmed %in% TRUE, na.rm = TRUE), 1L)
  assign(".PCT_ONSET_IMPUTED", .pct_imp, envir = .GlobalEnv)
  # SPLIT BY PROVENANCE. The overall share mixes two different things and would mislead if
  # quoted as a property of DHIS2 reporting: the sitrep-reconciliation rows appended by
  # .build_sitrep_confirmed_appends() carry date_of_symptom_onset = NA BY CONSTRUCTION (they
  # are a count reconciliation against the sitrep cumulative, not case records), so they are
  # 100% "imputed" by definition. Publishing both numbers lets prose state the genuine
  # line-list missingness separately from the share of the modelled series that is imputed.
  .is_sitrep <- grepl("^SITREP-CONF-", as.character(ll$alert_id))
  .conf <- ll$confirmed %in% TRUE
  .imp  <- ll$onset_imputed %in% TRUE
  assign(".PCT_ONSET_IMPUTED_LINELIST",
         100 * sum(.conf & .imp & !.is_sitrep, na.rm = TRUE) /
           max(sum(.conf & !.is_sitrep, na.rm = TRUE), 1L), envir = .GlobalEnv)
  assign(".N_SITREP_APPENDED", sum(.conf & .is_sitrep, na.rm = TRUE), envir = .GlobalEnv)
  # Publish the MECHANISM too, so report prose describes what the code actually did rather
  # than a hand-written sentence that silently outlives the implementation it describes.
  assign(".ONSET_IMPUTE_MODE", .imp_mode, envir = .GlobalEnv)
  assign(".ONSET_IMPUTE_DESC",
         switch(.imp_mode,
                resolver        = "the shared truncation- and interval-censoring-corrected DHIS2 onset-to-sample delay fit",
                resolver_growth = "the shared truncation- and interval-censoring-corrected DHIS2 onset-to-sample delay fit, growth-tilted",
                empirical       = "a fallback empirical bootstrap of the right-truncated pool of complete onset-sample pairs",
                growth          = "a growth-tilted fallback empirical bootstrap of the right-truncated pool of complete onset-sample pairs",
                parametric      = "a parametric exponential fallback",
                complete_case   = "no imputation (complete-case: onset-less records dropped)",
                none            = "no imputation (the sample date is used verbatim)",
                sprintf("an unrecognised imputation mode (%s)", .imp_mode)),
         envir = .GlobalEnv)
  # Report the censoring alongside the imputation, so a run's log states how many recorded
  # onsets were moved and by how much rather than leaving it to be rediscovered.
  if (any(ll$onset_censored %in% TRUE)) {
    .cz <- which(ll$onset_censored %in% TRUE)
    .shift <- as.numeric(ll$date_of_symptom_onset[.cz] - ll$date_index[.cz])
    message(sprintf(paste0("[load_linelist] Onset CENSORED at the specimen date for %d record(s) ",
                           "(%d confirmed) recorded up to %d d after their own specimen ",
                           "(ONSET_SAMPLE_NEG_TOL_DAYS = %d); shift %d-%d d earlier, median %g."),
                    length(.cz), sum(ll$confirmed[.cz] %in% TRUE),
                    max(.shift), .neg_tol, min(.shift), max(.shift), stats::median(.shift)))
  }
  message("[load_linelist] Onset imputed for ", sum(ll$onset_imputed, na.rm = TRUE),
          " records via ",
          # Every reachable .imp_mode needs an arm. "growth" and "complete_case" are both
          # reachable (see the .imp_mode assignment above) and had none, so switch() returned
          # NULL and the log read "...Onset imputed for N records via ; drawn delay mean...".
          switch(.imp_mode,
                 resolver   = sprintf("a draw from the SHARED delay resolver: %s", describe_delay(.resolver)),
                 resolver_growth = sprintf("a growth-tilted (r=%.4f/d) draw from the SHARED delay resolver: %s",
                                           .r_growth, describe_delay(.resolver)),
                 empirical  = sprintf("a FALLBACK draw from the right-truncated EMPIRICAL pool (%d %s pairs, source %s)",
                                      length(.dd), if (.windowed) "windowed" else "all-complete", .delay_src),
                 growth     = sprintf("a growth-tilted FALLBACK draw from the right-truncated EMPIRICAL pool (%d %s pairs, r=%.4f/d)",
                                      length(.dd), if (.windowed) "windowed" else "all-complete", .r_growth),
                 parametric = sprintf("a parametric draw from %s", .param_desc),
                 complete_case = "n/a (complete-case: onset-less records dropped)",
                 none       = "the sample date (no delay available)",
                 sprintf("an unrecognised imputation mode (%s)", .imp_mode)),
          if (length(.imp_draws))
            sprintf("; drawn delay mean %.1f d, range %d-%d d", mean(.imp_draws),
                    as.integer(min(.imp_draws)), as.integer(max(.imp_draws)))
          else "",
          "; onset-dated target.")

  message(
    "[load_linelist] Confirmed: ", sum(ll$confirmed, na.rm = TRUE),
    "  Suspected: ", sum(ll$suspected, na.rm = TRUE),
    "  Zones represented: ", length(unique(ll$health_zone[!is.na(ll$health_zone)]))
  )

  ll
}


# =============================================================================
# SECTION 3: aggregate_to_zone_week()
# =============================================================================
#
# Aggregates the linelist to a zone × week matrix (7-day buckets, WEEK_ANCHOR-anchored) and zero-fills all
# 519 zone-week cells that have no observations, so downstream models see
# a complete rectangular data frame.
#
# Arguments:
#   ll    — output of load_linelist()
#   zones — character vector of canonical zone names (e.g. from load_population())
#
# Returns a tibble sorted by health_zone, week_start with columns:
#   health_zone, week_start (bucket-start Date; buckets end on ANALYSIS_DATE), confirmed, suspected,
#   total_alerts, tests_analyzed, positivity

#' @param asof as-of date bounding the week grid. Defaults to the global ANALYSIS_DATE for the
#'   deployed call. reaggregate_asof() (22_daily_reissue.R) MUST pass its own issue date: this
#'   function is called there with a line list censored to a FOLD cutoff, and taking the bound
#'   from the global would zero-pad each fold's grid all the way to the analysis week. Every
#'   current caller re-filters afterwards, so no number moves — but the function must be a
#'   function of its arguments, not of a global.
aggregate_to_zone_week <- function(ll, zones,
                                   asof = get0("ANALYSIS_DATE", ifnotfound = NA)) {

  # ---- input checks ----------------------------------------------------------
  stopifnot(is.data.frame(ll))
  .check_cols(ll, c("health_zone", "date_index", "confirmed", "suspected",
                    "final_mve_case_classification"), "aggregate_to_zone_week")
  if (length(zones) == 0) {
    stop("[aggregate_to_zone_week] `zones` vector is empty.", call. = FALSE)
  }

  # ---- compute week starts (anchored so weeks END on the analysis date) ------
  ll <- ll %>%
    dplyr::mutate(
      week_start = lubridate::floor_date(date_index, unit = "week",
                                         week_start = get0("WEEK_ANCHOR", ifnotfound = 1L))
    ) %>%
    dplyr::filter(!is.na(week_start))

  # Ensure samples_analyzed column exists (optional in some DHIS2 exports)
  if (!"samples_analyzed" %in% names(ll)) ll$samples_analyzed <- NA_real_

  # ---- aggregate observed zone-week cells ------------------------------------
  obs <- ll %>%
    dplyr::group_by(health_zone, week_start) %>%
    dplyr::summarise(
      confirmed     = sum(confirmed,  na.rm = TRUE),
      suspected     = sum(suspected,  na.rm = TRUE),
      total_alerts  = dplyr::n(),
      tests_analyzed = sum(suppressWarnings(as.numeric(samples_analyzed)),
                           na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # A PROPORTION or nothing. `confirmed > tests_analyzed` is not a high positivity rate,
      # it is evidence that the denominator is not the tests behind this numerator — which is
      # exactly what happened when the synthetic sitrep rows carried samples_analyzed = 1.
      # DHIS2 carries no test counts for real records, so on the current data this is NA
      # everywhere, which is the honest answer.
      positivity = dplyr::if_else(
        tests_analyzed > 0 & confirmed <= tests_analyzed,
        confirmed / tests_analyzed,
        NA_real_
      )
    )

  # ---- build complete zone × week grid and zero-fill ------------------------
  .obs_weeks <- sort(unique(obs$week_start))

  if (length(.obs_weeks) == 0) {
    stop("[aggregate_to_zone_week] No valid week_start dates derived from date_index.", call. = FALSE)
  }

  # Complete the weekly sequence (min..max in 7-day steps) rather than keeping only weeks that
  # happen to carry observations: a fully case-free INTERIOR calendar week must NOT be dropped,
  # because compute_foi's generation-time convolution and the LFO horizon->date alignment index
  # weeks by POSITION on a contiguous 7-day grid — an interior gap would misweight the GT lag and
  # desync the truth windows. Observed week_starts are already 7-day-anchored, so this only inserts
  # (zero-filled, in tidyr::complete below) any interior gaps; it is a no-op when weeks are dense.
  # TERMINAL WEEK ANCHORED TO THE ANALYSIS DATE, not to whatever the data happen to contain.
  # Two silent failure modes otherwise, both one bad row away (verified clean on the current
  # snapshot, where max(date_index) == ANALYSIS_DATE):
  #   (a) OVER-EXTENSION. load_linelist() bounds date_index below (>= OUTBREAK_START) but not
  #       above, and onset_usable only requires onset <= sample + 2. A single record with onset
  #       AND sample both mistyped into the future is therefore "usable", and the grid runs past
  #       ANALYSIS_DATE — adding thousands of structurally-zero cells and, worse, moving the
  #       grid's terminal week off the as-of week the LFO cutoff and the nowcast both assume.
  #   (b) MISSING TERMINAL WEEK. If no record anywhere has a date_index inside the final 7-day
  #       window, that week vanishes and every POSITIONAL week index shifts by one — which
  #       matters because compute_foi() indexes weeks by position, as the note above says.
  .asof <- suppressWarnings(as.Date(asof))
  .last_wk <- if (length(.asof) == 1L && !is.na(.asof))
    lubridate::floor_date(.asof, "week", week_start = get0("WEEK_ANCHOR", ifnotfound = 1L))
    else max(.obs_weeks)
  if (max(.obs_weeks) > .last_wk) {
    warning(sprintf("[aggregate_to_zone_week] %d observed week(s) start AFTER the as-of week (%s) and are dropped: check for future-dated onset/sample records.",
                    sum(.obs_weeks > .last_wk), format(.last_wk)), call. = FALSE)
    obs <- dplyr::filter(obs, week_start <= .last_wk)
    .obs_weeks <- .obs_weeks[.obs_weeks <= .last_wk]
    if (!length(.obs_weeks))
      stop("[aggregate_to_zone_week] No observed weeks at or before the as-of week.", call. = FALSE)
  }
  all_weeks <- seq(min(.obs_weeks), max(.last_wk, max(.obs_weeks)), by = 7L)

  message(
    "[aggregate_to_zone_week] Observed weeks: ",
    min(.obs_weeks), " to ", max(.obs_weeks),
    " (", length(.obs_weeks), " observed",
    if (length(all_weeks) > length(.obs_weeks))
      paste0("; +", length(all_weeks) - length(.obs_weeks), " zero-filled interior week(s) for contiguity")
    else "",
    ", ", length(all_weeks), " total)"
  )
  message("[aggregate_to_zone_week] Completing grid for ", length(zones), " zones × ",
          length(all_weeks), " weeks")

  # Warn about zone names in the linelist not in the canonical spine
  extra_zones <- setdiff(unique(obs$health_zone), zones)
  if (length(extra_zones) > 0) {
    warning(
      "[aggregate_to_zone_week] ", length(extra_zones),
      " zone(s) in linelist not in population spine — dropped: ",
      paste(extra_zones, collapse = ", "), call. = FALSE
    )
  }

  zone_week <- tidyr::complete(
    obs,
    health_zone = zones,
    week_start  = all_weeks,
    fill        = list(
      confirmed      = 0L,
      suspected      = 0L,
      total_alerts   = 0L,
      tests_analyzed = 0,
      positivity     = NA_real_
    )
  ) %>%
    dplyr::filter(health_zone %in% zones) %>%   # drop any obs zones outside canonical set
    dplyr::arrange(health_zone, week_start)

  n_total   <- nrow(zone_week)
  n_nonzero <- sum(zone_week$total_alerts > 0)
  message(
    "[aggregate_to_zone_week] Grid: ", n_total, " zone-week cells; ",
    n_nonzero, " have >=1 alert (", round(100 * n_nonzero / n_total, 1), "% active)"
  )

  zone_week
}


# =============================================================================
# SECTION 4: load_population()
# =============================================================================
#
# Reads WorldPop processed CSV and returns a named numeric vector:
#   pop_vec[zone_name] = population count
#
# The WorldPop file has 519 health zones for all of DRC.

load_population <- function() {

  pop_path <- file.path(WORLDPOP_DIR, "worldpop__pop_count__static.csv")

  if (!file.exists(pop_path)) {
    stop("[load_population] WorldPop file not found: ", pop_path, call. = FALSE)
  }

  pop <- tryCatch(
    readr::read_csv(
      pop_path,
      col_types = readr::cols(
        nom       = readr::col_character(),
        pop_count = readr::col_double()
      ),
      show_col_types = FALSE
    ),
    error = function(e) stop("[load_population] Read error: ", e$message, call. = FALSE)
  )

  .check_cols(pop, c("nom", "pop_count"), "load_population")

  pop <- pop %>%
    dplyr::filter(!is.na(nom), !is.na(pop_count))

  if (nrow(pop) == 0) {
    stop("[load_population] WorldPop CSV is empty after filtering NA rows.", call. = FALSE)
  }

  pop_vec <- setNames(pop$pop_count, pop$nom)

  message(
    "[load_population] ", length(pop_vec), " zones loaded; ",
    "total DRC population: ", format(round(sum(pop_vec)), big.mark = ",")
  )

  pop_vec
}


# =============================================================================
# SECTION 5: load_static_covariates()
# =============================================================================
#
# Reads the following static covariate files and joins them on `nom`:
#   - CCVI socioeconomic deprivation
#   - GDP per capita
#   - GRID3 healthsites count and density
#   - PCR testing capacity
#   - Cross-border passenger volumes
#
# Returns a tibble with `nom` as the join key, one row per health zone.
# Missing data for a given zone is allowed (NA in the relevant column).

load_static_covariates <- function() {

  # Helper: read a two-column (index + nom + value) static covariate CSV.
  # Returns a tibble with nom + one value column, skipping the row-index column.
  read_static <- function(path, value_col, context) {
    if (!file.exists(path)) {
      warning("[load_static_covariates] File not found, skipping: ", path, call. = FALSE)
      return(NULL)
    }
    tryCatch(
      {
        df <- readr::read_csv(
          path,
          col_types = readr::cols(.default = readr::col_character()),
          show_col_types = FALSE
        )
        .check_cols(df, c("nom", value_col), context)
        df %>%
          dplyr::select(nom, dplyr::all_of(value_col)) %>%
          dplyr::filter(!is.na(nom), nom != "") %>%
          dplyr::mutate(across(dplyr::all_of(value_col),
                               ~ suppressWarnings(as.numeric(.x))))
      },
      error = function(e) {
        warning("[load_static_covariates] Error reading ", path, ": ", e$message, call. = FALSE)
        NULL
      }
    )
  }

  # ---- Read each covariate file --------------------------------------------
  ccvi <- read_static(
    file.path(CCVI_DIR, "ccvi__socioeconomic_deprivation__static.csv"),
    "socioeconomic_deprivation",
    "CCVI"
  )

  gdp <- read_static(
    file.path(GDP_DIR, "gdp_pc__gdp_pc__static.csv"),
    "gdp_pc",
    "GDP"
  )

  healthsites_count <- read_static(
    file.path(HEALTHSITES_DIR, "grid3_healthsites__healthsite_count__static.csv"),
    "healthsite_count",
    "healthsites_count"
  )

  healthsites_density <- read_static(
    file.path(HEALTHSITES_DIR, "grid3_healthsites__healthsite_density__static.csv"),
    "healthsite_density",
    "healthsites_density"
  )

  testing <- read_static(
    file.path(TESTING_DIR, "testing_capacity__pcr_tests__static.csv"),
    "pcr_tests",
    "testing_capacity"
  )

  cross_border_path <- file.path(
    DATA_DIR, "cross-border-movements", "processed",
    "cross_border__poe_passengers__static.csv"
  )
  cross_border <- if (file.exists(cross_border_path)) {
    tryCatch(
      {
        cb <- readr::read_csv(
          cross_border_path,
          col_types = readr::cols(.default = readr::col_character()),
          show_col_types = FALSE
        )
        if ("nom" %in% names(cb) && "mean_daily_passengers" %in% names(cb)) {
          cb %>%
            dplyr::select(nom, mean_daily_passengers) %>%
            dplyr::mutate(mean_daily_passengers =
                            suppressWarnings(as.numeric(mean_daily_passengers)))
        } else if ("nom" %in% names(cb)) {
          # Pick first numeric-ish column as the passenger volume
          num_cols <- names(cb)[names(cb) != "nom"]
          if (length(num_cols) > 0) {
            cb %>%
              dplyr::select(nom, dplyr::all_of(num_cols[1])) %>%
              dplyr::rename(mean_daily_passengers = dplyr::all_of(num_cols[1])) %>%
              dplyr::mutate(mean_daily_passengers =
                              suppressWarnings(as.numeric(mean_daily_passengers)))
          } else NULL
        } else NULL
      },
      error = function(e) {
        warning("[load_static_covariates] Cross-border read error: ", e$message, call. = FALSE)
        NULL
      }
    )
  } else {
    warning("[load_static_covariates] Cross-border file not found: ", cross_border_path, call. = FALSE)
    NULL
  }

  # ---- Build full zone spine from WorldPop (519 zones) ----------------------
  pop_path <- file.path(WORLDPOP_DIR, "worldpop__pop_count__static.csv")
  if (!file.exists(pop_path)) {
    stop("[load_static_covariates] WorldPop file needed for zone spine: ", pop_path, call. = FALSE)
  }
  spine <- readr::read_csv(
    pop_path,
    col_types = readr::cols(nom = readr::col_character(), pop_count = readr::col_double()),
    show_col_types = FALSE
  ) %>%
    dplyr::select(nom, pop_count) %>%
    dplyr::filter(!is.na(nom))

  # ---- Sequential left-joins onto the zone spine ----------------------------
  cov <- spine
  join_one <- function(base, df, label) {
    if (is.null(df)) {
      message("[load_static_covariates] Skipped (NULL): ", label)
      return(base)
    }
    n_matched <- sum(df$nom %in% base$nom)
    # A duplicated `nom` in any covariate CSV would silently ROW-MULTIPLY the 519-zone spine,
    # and n_matched cannot detect it. Collapse duplicates first and say so, then assert the
    # spine is unchanged — a covariate join must never alter the number of zones.
    if (anyDuplicated(df$nom)) {
      warning(sprintf("[load_static_covariates] %s: %d duplicate zone name(s); keeping the first row of each.",
                      label, sum(duplicated(df$nom))), call. = FALSE)
      df <- dplyr::distinct(df, nom, .keep_all = TRUE)
    }
    message("[load_static_covariates] Joining ", label,
            " — ", n_matched, "/", nrow(df), " zones matched")
    out <- dplyr::left_join(base, df, by = "nom")
    if (nrow(out) != nrow(base))
      stop(sprintf("[load_static_covariates] %s changed the zone spine (%d -> %d rows).",
                   label, nrow(base), nrow(out)), call. = FALSE)
    out
  }

  cov <- join_one(cov, ccvi,              "CCVI socioeconomic_deprivation")
  cov <- join_one(cov, gdp,               "GDP gdp_pc")
  cov <- join_one(cov, healthsites_count, "healthsite_count")
  cov <- join_one(cov, healthsites_density, "healthsite_density")
  cov <- join_one(cov, testing,           "pcr_tests")
  cov <- join_one(cov, cross_border,      "cross_border mean_daily_passengers")

  # ---- Derived, model-facing covariate columns ------------------------------
  # Downstream models (hhh4 C1d, INLA C4, ZINB C5) reference these canonical
  # names. Deriving them here keeps a single source of truth and avoids each
  # model silently falling back to constant covariates (rank-deficiency risk).
  # Guard each derived column: read_static()/join_one() deliberately warn-and-skip
  # a missing covariate file (so its source column is never created); referencing
  # it unconditionally here would crash the whole data prep the moment a layer is
  # absent — exactly the situation the defensive loading anticipates.
  cov <- cov %>%
    dplyr::mutate(
      log_pop              = if ("pop_count" %in% names(cov)) log(pmax(pop_count, 1)) else NA_real_,
      ccvi                 = if ("socioeconomic_deprivation" %in% names(cov)) socioeconomic_deprivation else NA_real_,
      log_healthsite_count = if ("healthsite_count" %in% names(cov)) log1p(dplyr::coalesce(healthsite_count, 0)) else NA_real_
    )

  message("[load_static_covariates] Covariate tibble: ",
          nrow(cov), " zones × ", ncol(cov) - 1, " covariates (excl. nom)")

  cov
}


# =============================================================================
# SECTION 6: load_sitrep()
# =============================================================================
#
# Reads the INSP sitrep daily new confirmed cases file, aggregates to weekly
# totals (7-day buckets, WEEK_ANCHOR-anchored, ending on ANALYSIS_DATE), and returns a tibble with:
#   nom, week_start, sitrep_confirmed_weekly
#
# This is OPTIONAL: if the file does not exist, returns NULL with a warning.

load_sitrep <- function() {

  sitrep_path <- file.path(SITREP_DIR, "insp_sitrep__new_confirmed_cases__daily.csv")

  if (!file.exists(sitrep_path)) {
    warning(
      "[load_sitrep] INSP sitrep file not found — skipping: ", sitrep_path,
      call. = FALSE
    )
    return(NULL)
  }

  sitrep <- tryCatch(
    readr::read_csv(
      sitrep_path,
      col_types = readr::cols(
        nom                = readr::col_character(),
        date               = readr::col_character(),
        new_confirmed_cases = readr::col_double()
      ),
      show_col_types = FALSE
    ),
    error = function(e) {
      warning("[load_sitrep] Read error: ", e$message, call. = FALSE)
      return(NULL)
    }
  )

  if (is.null(sitrep)) return(NULL)

  .check_cols(sitrep, c("nom", "date", "new_confirmed_cases"), "load_sitrep")

  # HARMONISE the zone names. Every other zone-keyed loader canonicalises through aliases.csv;
  # this one only RENAMED the column, which is not the same thing — so the returned tibble
  # carried raw sitrep spellings ("Nia-Nia" against the spine's "Nia Nia") that can never join
  # the 519-zone spine, even though the alias for exactly that case exists. Also apply the as-of
  # filter that .build_sitrep_confirmed_appends() applies on the same source.
  .sit_aliases <- if (file.exists(ALIASES_PATH))
    tryCatch(readr::read_csv(ALIASES_PATH,
                             col_types = readr::cols(.default = readr::col_character()),
                             show_col_types = FALSE),
             error = function(e) {
               warning("[load_sitrep] Cannot read aliases.csv: ", e$message, call. = FALSE); NULL })
    else NULL
  .sit_asof <- suppressWarnings(as.Date(get0("ANALYSIS_DATE", ifnotfound = NA)))

  sitrep <- sitrep %>%
    dplyr::mutate(
      date       = .parse_date(date),
      nom        = if (!is.null(.sit_aliases)) .apply_aliases(nom, .sit_aliases) else nom,
      week_start = lubridate::floor_date(date, unit = "week",
                                         week_start = get0("WEEK_ANCHOR", ifnotfound = 1L))
    ) %>%
    dplyr::filter(!is.na(date), !is.na(nom), nom != "",
                  is.na(.sit_asof) | date <= .sit_asof) %>%
    dplyr::group_by(nom, week_start) %>%
    dplyr::summarise(
      sitrep_confirmed_weekly = sum(new_confirmed_cases, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(nom, week_start) %>%
    dplyr::rename(health_zone = nom)  # align key name with zone_week

  # min()/max() on an empty tibble return -Inf/Inf with a warning rather than failing cleanly.
  if (nrow(sitrep) == 0) {
    message("[load_sitrep] No sitrep rows at or before the analysis date.")
    return(NULL)
  }

  message(
    "[load_sitrep] Sitrep: ", length(unique(sitrep$health_zone)), " zones, ",
    nrow(sitrep), " zone-week rows, weeks ",
    min(sitrep$week_start), " to ", max(sitrep$week_start)
  )

  sitrep
}


# =============================================================================
# SECTION 8: prep_all_data() — master loader
# =============================================================================
#
# Calls all loaders in order, returns a named list for use by all downstream
# spatiotemporal scripts. Zone list for zero-filling is derived from WorldPop
# (519 DRC health zones), ensuring the observation matrix is always complete.

prep_all_data <- function() {

  message("\n", strrep("=", 70))
  message("[prep_all_data] Loading all spatiotemporal data sources")
  message(strrep("=", 70))

  # 1. Population (provides the canonical 519-zone list)
  pop <- load_population()
  zones <- names(pop)

  # 2. Linelist
  ll <- load_linelist()

  # 3. Zone × week aggregation (zero-filled to all 519 zones)
  zone_week <- aggregate_to_zone_week(ll, zones)

  # 4. Static covariates
  covariates <- load_static_covariates()

  # 4b. Per-zone test positivity (a model covariate; derived from the zone-week
  #     aggregation since it is not a static input file). Zones with no testing data
  #     fall back to the nominal ascertainment rate downstream.
  #     POOLED ratio (sum confirmed / sum tests), not the unweighted mean of the weekly
  #     ratios: a week with one test must not carry the same weight as a week with 500.
  #     Only weeks with a usable denominator contribute.
  zone_positivity <- zone_week %>%
    dplyr::filter(is.finite(positivity), tests_analyzed > 0) %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(.conf = sum(confirmed, na.rm = TRUE),
                     .tests = sum(tests_analyzed, na.rm = TRUE), .groups = "drop") %>%
    dplyr::transmute(health_zone,
                     positivity = dplyr::if_else(.tests > 0, .conf / .tests, NA_real_))
  covariates <- covariates %>%
    dplyr::left_join(zone_positivity, by = c("nom" = "health_zone"))

  # 5. INSP sitrep (optional)
  sitrep <- load_sitrep()

  # ---- Summary statistics ---------------------------------------------------
  message("\n", strrep("-", 70))
  message("[prep_all_data] SUMMARY STATISTICS")
  message(strrep("-", 70))

  confirmed_ll  <- sum(ll$confirmed, na.rm = TRUE)
  zones_affected <- ll %>%
    dplyr::filter(confirmed) %>%
    dplyr::pull(health_zone) %>%
    unique() %>%
    length()
  date_range <- range(ll$date_index, na.rm = TRUE)

  message(
    "  Total confirmed cases (linelist): ", confirmed_ll
  )
  message(
    "  Zones with >=1 confirmed case:    ", zones_affected,
    " / ", length(zones), " total zones"
  )
  message(
    "  Date index range:                 ",
    date_range[1], " to ", date_range[2]
  )
  message(
    "  Zone-week matrix:                 ",
    nrow(zone_week), " rows (",
    length(unique(zone_week$health_zone)), " zones × ",
    length(unique(zone_week$week_start)), " weeks)"
  )
  if (!is.null(sitrep)) {
    message(
      "  Sitrep weekly rows:               ", nrow(sitrep)
    )
  }
  message(strrep("=", 70), "\n")

  list(
    ll         = ll,
    zone_week  = zone_week,
    pop        = pop,
    zones_all  = zones,          # canonical 519-zone name vector (= names(pop))
    covariates = covariates,
    sitrep     = sitrep
  )
}

# =============================================================================
# SECTION 9: Execute
# =============================================================================

if (interactive() && !exists("dat")) {
  dat <- prep_all_data()
  message("[01_data_prep.R] Complete. Object `dat` available with components: ",
          paste(names(dat), collapse = ", "))
}
