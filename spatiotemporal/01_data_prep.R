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

# The canonical modelling spine (WorldPop zone names), or NULL when the file is unreadable.
.zone_spine_names <- function() {
  tryCatch(
    readr::read_csv(file.path(WORLDPOP_DIR, "worldpop__pop_count__static.csv"),
      col_types = readr::cols(nom = "c", pop_count = "d"), show_col_types = FALSE)$nom,
    error = function(e) NULL)
}

# Case-, accent- and punctuation-insensitive key for comparing province names
# ("Kasai Central" and "Kasaï-Central" compare equal).
.province_key <- function(x) {
  # stringi transliterates in any locale; chartr() with an accented `old` string fails under
  # LANG=C ("'old' is longer than 'new'": the literal is read as bytes there).
  x <- tolower(stringi::stri_trans_general(trimws(as.character(x)), "Latin-ASCII"))
  gsub("[^a-z0-9]", "", x)
}

# Resolve, by the record's own PROVINCE, a bare zone name that the spine carries only in its
# "Name (Province)" form. The spine separates its same-named zones that way ("Bili (Bas-Uele)"
# and "Bili (Nord-Ubangi)", likewise Lubunga), while the line list writes the bare name and
# records the province in its own column; an exact-name match therefore never finds them, and
# every such record fell out of the zone-week grid and the alert counts. A name the spine
# carries as written, a bare name whose province matches no candidate (or is missing), and
# every other name are returned unchanged. Vectorised over records.
.resolve_zone_by_province <- function(zone, province, spine) {
  out <- zone
  if (is.null(spine) || is.null(province) || !length(zone)) return(out)
  qual <- spine[grepl(" \\([^()]*\\)$", spine)]
  if (!length(qual)) return(out)
  key <- paste(sub(" \\([^()]*\\)$", "", qual),
               .province_key(sub("^.*\\(([^()]*)\\)$", "\\1", qual)), sep = "\r")
  qual <- qual[!key %in% key[duplicated(key)]]           # a (name, province) pair must be unique
  key  <- key[!key %in% key[duplicated(key)]]
  cand <- !is.na(zone) & !(zone %in% spine)
  hit  <- qual[match(paste(zone[cand], .province_key(province[cand]), sep = "\r"), key)]
  out[cand] <- ifelse(is.na(hit), zone[cand], hit)
  out
}

# The per-record key of the onset-imputation draw (.uniform_from_ids(), below): the record's
# alert_id, which the draw is designed to depend on. A record WITHOUT an alert_id would hash to
# the same uniform as every other id-less record, so in any donor week all of them would take
# the same interval; each is keyed instead on a hash of its own fields AS EXTRACTED (call this
# before any column is coerced or replaced, so the key does not move with SAMPLE_DATE_SOURCE
# or CASE_ZONE_BASIS). The key, like the id, depends neither on row order nor on other records.
.onset_draw_keys <- function(df) {
  key <- as.character(df$alert_id)
  noid <- which(is.na(key))
  if (length(noid)) {
    flds <- df[noid, setdiff(names(df), ".onset_draw_key"), drop = FALSE]
    key[noid] <- paste0("NOID-", vapply(do.call(paste, c(unname(as.list(flds)), sep = "\x1f")),
                                        rlang::hash, character(1), USE.NAMES = FALSE))
  }
  key
}

# NOTE: .load_dhis2_delay_params() and .draw_dhis2_delay() MOVED to 00_config.R (sourced
# above), alongside delay_cdf() and effective_onset_sample_delay(). They lived here, but
# 04_nowcasting.R and 02_epi_params.R need the same delay and could not reach a definition
# private to the data-prep module — so they silently kept using the fixed lab Exp reference
# while this file used the fitted DHIS2 delay. One shared definition removes that whole class
# of drift: the nowcast weights' fallback (delay_cdf) and the R(t) truncation's legacy mode read
# ONE resolver, effective_onset_sample_delay(). (The missing-onset imputation no longer uses a
# fitted delay: it draws from donor intervals, .impute_onset_delays() below.)


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
# SOURCE: sitrep_cumulative_confirmed() above — the CUMULATIVE confirmed file, spelling
# variants merged by MAX-per-date, a running maximum over dates, diffed into dated positive
# increments: one confirmed case per unit increment, dated at the sitrep report date on
# which it appeared.
#
# PLACEMENT & DATING: the shortfall is filled with the MOST RECENT sitrep confirmed
# dates (reporting lag ⇒ the line list already holds the earlier cases). Each appended
# row's `date_of_sample_collection` is set to that sitrep report date (SITREP_APPEND_SHIFT_DAYS
# earlier under that sensitivity) and its onset is LEFT BLANK, so load_linelist() imputes
# onset = sample − an interval drawn from the recorded intervals of same-class records sampled
# in the same week (.impute_onset_delays()) — exactly as it treats real records lacking a
# usable onset (same convention as build_conditional_linelist.R).
#
# IDEMPOTENT vs the conditional pre-append: build_conditional_linelist.R already
# injects some sitrep zones into the conditional CSV; those raise `n_ll` so their
# shortfall is 0 here — they are not re-appended. Toggle with APPEND_SITREP_CONFIRMED.
#
# Args:  ll          — current line list (outbreak-filtered, unmatched-excluded, health_zone
#                       canonicalised by province and aliases), with final_mve_case_classification.
#        aliases_all — full aliases tibble (observed_name, canonical_nom) or NULL;
#                       used to canonicalise sitrep `nom` to the modelling spine.
# Returns a tibble (0 rows if nothing to append), safe to dplyr::bind_rows() onto `ll`. A
#   column it carries that `ll` lacks (e.g. genexpert_result, or date_known, added only when
#   SITREP_APPEND_SHIFT_DAYS > 0) is NA on the line-list rows after the bind.
.build_sitrep_confirmed_appends <- function(ll, aliases_all) {

  if (!isTRUE(get0("APPEND_SITREP_CONFIRMED", ifnotfound = TRUE))) return(tibble::tibble())

  # Canonical modelling spine (WorldPop 519 zones): a sitrep zone absent from it cannot
  # enter the zone-week grid, so we skip (and warn about) it rather than fabricate an
  # unmatched row that aggregate_to_zone_week would silently drop.
  spine <- .zone_spine_names()

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

  # SITREP_APPEND_SHIFT_DAYS (00_config.R): the report date follows the specimen, so the
  # sensitivity dates the specimen `.shift` days earlier while the row still becomes KNOWN only on
  # the report date (date_known, read first by linelist_observation_date()). The column is added
  # only when the shift is on, so the default line list is unchanged.
  .shift <- as.integer(get0("SITREP_APPEND_SHIFT_DAYS", ifnotfound = 0L))
  mk <- function(z, short) {
    d <- pick_dates(z, short)
    out <- tibble::tibble(
      country                       = "Democratic Republic of the Congo",
      province                      = unname(.prov_lookup[z]),
      health_zone                   = z,
      health_zone_unmatched         = FALSE,
      health_area_unmatched         = FALSE,
      date_of_sample_collection     = as.Date(d) - .shift,
      date_of_symptom_onset         = as.Date(NA),   # imputed downstream from donor intervals
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
    if (.shift > 0L) out$date_known <- as.Date(d)
    out
  }
  appended <- purrr::map2_dfr(recon$nom, recon$shortfall, mk)

  if (.shift > 0L)
    message(sprintf(paste0("[load_linelist] SITREP_APPEND_SHIFT_DAYS = %d: the %d appended row(s) are ",
                           "dated %d day(s) before their sitrep report date (known on the report date)."),
                    .shift, nrow(appended), .shift))
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
# filters to the current BDBV 2026 outbreak (>= OUTBREAK_START), resolves zone names
# (by province for same-named zones, then aliases), and excludes records with unmatched
# health zones.
#
# Returns a tibble with:
#   - All original columns (dates coerced to Date class)
#   - `date_index`    : the usable recorded onset (censored at the sample date), else the
#                       sample date minus an interval drawn from same-class, same-week donors,
#                       floored at OUTBREAK_START; NA when neither exists, when no donor exists,
#                       when it falls after ANALYSIS_DATE, or (ONSET_MODE = "complete_case")
#                       when the onset is not usable
#   - `onset_usable`, `onset_imputed`, `onset_censored` : logical flags for the above
#   - `confirmed`     : logical — final_mve_case_classification == CONFIRMED_STATUS
#   - `suspected`     : logical — final_mve_case_classification == SUSPECTED_STATUS
#   - `health_zone`   : canonical zone name (province-resolved, alias-corrected)

# .read_linelist_records(): the record-level half of load_linelist() — read, type coercion,
# unmatched-zone exclusion and zone-name resolution, with NO date filter. Split out (2026-10-03)
# so that load_unplaced_alerts() sees the SAME zone matching and aliasing as the case data. The
# outbreak filter used to run before the exclusion and the aliasing; all three are row-wise
# and order-preserving, so moving it after them leaves load_linelist()'s output identical
# (verified with identical() on the LINELIST_23092026 snapshot) and only shifts the counts
# its log lines report. It also carries `.onset_draw_key` (.onset_draw_keys()), which
# load_linelist() consumes and drops.
.read_linelist_records <- function(folder = NULL) {

  # ---- resolve path: an explicit snapshot folder, else latest.json -----------
  # `folder` lets the back-test read an ARCHIVED snapshot (LFO_TRAINING_SOURCE = "vintage",
  # load_linelist_vintage(), below); the default is the live snapshot.
  if (is.null(folder)) {
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
    folder <- meta$folder
  }
  ll_csv <- file.path(LINELIST_DIR, folder, "dhis2_processed_linelist.csv")
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
  # alert_id is required: the onset-imputation draw, the earliest-sample lookup and the sitrep
  # top-up's identification all key on it.
  required_cols <- c(
    "alert_id", "health_zone", "health_zone_unmatched",
    "date_of_sample_collection", "date_of_symptom_onset",
    "final_mve_case_classification"
  )
  .check_cols(ll, required_cols, "load_linelist")
  # Before any column is coerced or replaced (see .onset_draw_keys()).
  ll$.onset_draw_key <- .onset_draw_keys(ll)

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

  # DIAGNOSTIC SAMPLE DATE (SAMPLE_DATE_SOURCE, 00_config.R): the earliest sample date recorded
  # for each alert across the snapshot archive, not the extract's (possibly repeat-test) date.
  ll <- apply_earliest_sample_dates(ll, label = "load_linelist")

  # Numeric columns
  if ("samples_analyzed" %in% names(ll)) {
    ll$samples_analyzed <- suppressWarnings(as.numeric(ll$samples_analyzed))
  }
  if ("samples_received" %in% names(ll)) {
    ll$samples_received <- suppressWarnings(as.numeric(ll$samples_received))
  }

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

  # ---- same-named zones: resolve the bare name by province --------------------
  # BEFORE the aliases, because an alias is province-blind (aliases.csv maps a bare
  # "Lubunga" to "Lubunga (Tshopo)"), while the record's own province tells the two apart.
  .spine <- .zone_spine_names()
  if ("province" %in% names(ll) && !is.null(.spine)) {
    old <- ll$health_zone
    ll$health_zone <- .resolve_zone_by_province(ll$health_zone, ll$province, .spine)
    ch <- !is.na(old) & old != ll$health_zone
    if (any(ch)) {
      tb <- table(paste0(old[ch], " -> ", ll$health_zone[ch]))
      message("[load_linelist] Same-named zones resolved by province: ", sum(ch), " records (",
              paste(sprintf("%s: %d", names(tb), as.integer(tb)), collapse = "; "), ").")
    }
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
  # A bare same-named zone name that neither its province nor an alias resolved is outside the
  # spine: every zone count would drop its records, at most with a generic off-spine warning.
  if (!is.null(.spine)) {
    .bare <- setdiff(unique(sub(" \\([^()]*\\)$", "", .spine[grepl(" \\([^()]*\\)$", .spine)])), .spine)
    .left <- ll$health_zone %in% .bare
    if (any(.left))
      warning(sprintf(paste0("[load_linelist] %d record(s) carry the bare name of a same-named zone (%s) ",
                             "that neither their province nor an alias resolves to one of its spine entries; ",
                             "they fall outside the zone spine and are dropped from every zone count."),
                      sum(.left), paste(unique(ll$health_zone[.left]), collapse = ", ")), call. = FALSE)
  }

  # The alias table travels with the records: load_linelist() needs it again for the sitrep
  # reconciliation below. NULL when aliases.csv is absent, exactly as before the split.
  attr(ll, "aliases") <- if (exists("aliases", inherits = FALSE)) aliases else NULL
  ll
}

#' CASE_ZONE_BASIS = "residence" (00_config.R): count each CONFIRMED record in the zone of its
#' recorded residence health area, where that resolves.
#'
#' The lookup is the MRT construction (mrt_affected_zones/R/02_cases.R::residence_lookup): every
#' record's REPORTING (health_area, health_zone) pair with a matched health area and an in-spine
#' zone, keeping only health areas that map to exactly one zone. A residence that is missing,
#' unrecognised or ambiguous leaves the record in its reporting zone; so do alerts that are not
#' confirmed (they describe where the system detects, not where cases live) and the sitrep top-up
#' rows (they carry no residence, and the sitrep counts them by reporting zone). The reporting
#' zone is kept in `health_zone_reporting`. A snapshot without the residence column (the archive
#' before 2026-07-17) keeps every reporting zone rather than borrowing residences from a later one.
#' So in this arm a back-test round trained on such a snapshot (LFO_TRAINING_SOURCE = "vintage")
#' sees REPORTING-zone counts, while its outcomes are scored on the final, residence-basis line list.
#' A confirmed record with no reporting zone but a resolved residence is placed in that residence.
#' @param lookup_records every record of the extract (before the outbreak filter).
.apply_residence_zone <- function(ll, lookup_records) {
  ll$health_zone_reporting <- ll$health_zone
  need <- c("health_area", "health_area_unmatched", "health_zone")
  if (!"permanent_residence_health_area" %in% names(ll) || !all(need %in% names(lookup_records))) {
    warning("[load_linelist] CASE_ZONE_BASIS = 'residence' but this extract has no residence (or ",
            "health-area) column; every case keeps its reporting zone.", call. = FALSE)
    return(ll)
  }
  spine <- .zone_spine_names()
  if (is.null(spine)) {
    warning("[load_linelist] CASE_ZONE_BASIS = 'residence' but the zone spine is unreadable; every ",
            "case keeps its reporting zone.", call. = FALSE)
    return(ll)
  }
  lk <- lookup_records[!is.na(lookup_records$health_area) &
                         lookup_records$health_area_unmatched %in% FALSE &
                         lookup_records$health_zone %in% spine, c("health_area", "health_zone")]
  lk$health_area <- trimws(lk$health_area)
  lk <- unique(lk)
  n_z <- table(lk$health_area)
  lk <- lk[lk$health_area %in% names(n_z)[n_z == 1L], , drop = FALSE]
  res  <- lk$health_zone[match(trimws(ll$permanent_residence_health_area), lk$health_area)]
  conf <- (ll$final_mve_case_classification == CONFIRMED_STATUS) %in% TRUE
  # NA-safe: `res != NA` is NA, and one NA in `move` breaks the any() and the assignment below.
  move <- conf & !is.na(res) & (is.na(ll$health_zone) | res != ll$health_zone)
  before <- unique(ll$health_zone[conf & !is.na(ll$health_zone)])
  # Province follows the zone (the commonest province among the zone's records), so nothing that
  # groups by (health_zone, province) splits a zone in two.
  if ("province" %in% names(ll) && any(move)) {
    pz <- stats::aggregate(n ~ health_zone + province,
                           data = transform(ll[!is.na(ll$province), c("health_zone", "province")], n = 1L),
                           FUN = sum)
    pz <- pz[order(-pz$n), ]; pz <- pz[!duplicated(pz$health_zone), ]
    newp <- pz$province[match(res[move], pz$health_zone)]
    ll$province[move] <- ifelse(is.na(newp), ll$province[move], newp)
  }
  ll$health_zone[move] <- res[move]
  after <- unique(ll$health_zone[conf & !is.na(ll$health_zone)])
  message(sprintf(paste0("[load_linelist] CASE_ZONE_BASIS = 'residence': %d of %d confirmed record(s) ",
                         "resolve a residence zone (lookup: %d unambiguous health areas, %d ambiguous ",
                         "dropped); %d live outside their reporting zone (or have none) and are MOVED. Zones with a ",
                         "confirmed case: %d gained, %d lost."),
                  sum(conf & !is.na(res)), sum(conf), nrow(lk), sum(n_z > 1L), sum(move),
                  length(setdiff(after, before)), length(setdiff(before, after))))
  ll
}

#' Uniform(0, 1) per record, a deterministic function of (id, seed). The imputation draw of a
#' record therefore does not depend on row order or on which other records are present, and
#' changing the seed (ONSET_IMPUTE_SEED, default RANDOM_SEED) gives an independent completion
#' (the multiple-imputation sweeps). Refuses a missing id: every NA would hash to ONE value, so
#' all id-less records would share a single draw (load_linelist() keys them by .onset_draw_keys()).
.uniform_from_ids <- function(ids, seed) {
  if (anyNA(ids))
    stop("[load_linelist] .uniform_from_ids(): ", sum(is.na(ids)), " missing id(s); each record ",
         "needs its own key (.onset_draw_keys()).", call. = FALSE)
  h <- vapply(as.character(ids), function(x) rlang::hash(c(x, as.character(seed))),
              character(1), USE.NAMES = FALSE)
  # The first 7 hex digits of the hash: an integer in [0, 16^7), always within strtoi()'s range.
  (strtoi(substr(h, 1L, 7L), base = 16L) + 0.5) / 16^7
}

#' Missing-onset imputation: the interval to subtract from a record's sample date.
#'
#' What is needed for a record SAMPLED on date s is the distribution of (sample - onset) AMONG
#' RECORDS SAMPLED AT s. That is the backward interval, and it is not the forward onset->sample
#' delay: while incidence grows, recently sampled cases have disproportionately short intervals.
#' It is observed directly. The donors are the records of the SAME CLASS (`group`: confirmed or
#' not) sampled in the SAME WEEK (the WEEK_ANCHOR grid) whose onset is recorded and usable; the
#' interval is drawn from theirs. Conditioning on the sample date, a donor's interval is fully
#' observed, so no delay family, growth-rate model or truncation correction enters.
#'
#' ASSUMPTION, stated: within a sample week and class, whether a record's onset was recorded does
#' not depend on its interval (missing at random given sample week and class).
#'
#' A sample week with fewer than `min_donors` donors borrows the adjacent weeks on both sides,
#' one week at a time, until it has enough (on the 2026-10-09 frame: the first two outbreak weeks
#' in both classes, and for the non-confirmed also two sample weeks dated before the outbreak).
#'
#' The draw for record i is the empirical quantile of its donor intervals at u_i =
#' .uniform_from_ids(id_i, seed): equivalent to picking a donor uniformly at random, but
#' reproducible record by record.
#'
#' @param anchor sample dates; @param recorded the usable recorded onsets (NA otherwise);
#' @param usable logical, recorded onset usable; @param group donor class of each record;
#' @param need logical, the records to impute; @param ids record keys, one per record and
#'   non-missing where `need` (alert_id, or .onset_draw_keys());
#' @param asof records sampled after it are not donors.
#' @return integer vector of intervals (days); NA where none is needed or no donor exists.
.impute_onset_delays <- function(anchor, recorded, usable, group, need, ids, seed, asof = NA,
                                 week_start = get0("WEEK_ANCHOR", ifnotfound = 1L),
                                 min_donors = get0("ONSET_DONOR_MIN", ifnotfound = 20L)) {
  anchor <- as.Date(anchor); recorded <- as.Date(recorded)
  if (length(ids) != length(anchor))
    stop("[load_linelist] .impute_onset_delays(): `ids` must have one entry per record.", call. = FALSE)
  out <- rep(NA_integer_, length(anchor))
  todo <- which(need %in% TRUE & !is.na(anchor))
  if (!length(todo)) return(out)
  wk <- as.numeric(lubridate::floor_date(anchor, "week", week_start = week_start))
  donor <- usable %in% TRUE & !is.na(anchor) & !is.na(recorded) &
           (is.na(asof) | anchor <= as.Date(asof))
  # >= 0 by construction (the recorded onset is censored at the sample date); pmax guards a
  # sample dated before the outbreak floor, where the floored onset would follow it.
  ivl <- pmax(as.integer(anchor - recorded), 0L)
  u <- .uniform_from_ids(ids[todo], seed)
  for (g in unique(group[todo])) {
    dg <- which(donor & group == g)
    if (!length(dg)) next
    tg <- todo[group[todo] == g]
    for (w in unique(wk[tg])) {
      k <- 0L
      repeat {
        sel <- abs(wk[dg] - w) <= 7 * k
        if (sum(sel) >= min_donors || all(sel)) break
        k <- k + 1L
      }
      pool <- sort(ivl[dg[sel]])
      ii <- tg[wk[tg] == w]
      out[ii] <- pool[floor(u[match(ii, todo)] * length(pool)) + 1L]
    }
  }
  out
}

load_linelist <- function(folder = NULL) {
  # `folder` = an archived snapshot; every as-of rule below reads the global ANALYSIS_DATE, which
  # load_linelist_vintage() (this file) sets to the snapshot's date for the duration.
  ll <- .read_linelist_records(folder)
  aliases <- attr(ll, "aliases"); attr(ll, "aliases") <- NULL
  # ---- alert_id: the key of the onset draw and of the earliest-sample lookup ----------
  # A missing id is NOT refused and NOT dropped: the archived snapshots of 2026-06-18 to
  # 2026-09-17 carry id-less records, confirmed cases among them, and dropping them would drop
  # cases. They are keyed for the draw by .onset_draw_keys() instead.
  .n_noid <- sum(is.na(ll$alert_id))
  if (.n_noid > 0L)
    warning(sprintf(paste0("[load_linelist] %d record(s) have no alert_id: they keep the extract's sample ",
                           "date (the earliest-sample lookup is by alert_id), and their onset-imputation ",
                           "draw is keyed on a hash of the record's own fields."), .n_noid), call. = FALSE)
  .dup <- !is.na(ll$alert_id) & ll$alert_id %in% ll$alert_id[duplicated(ll$alert_id)]
  if (any(.dup))
    warning(sprintf(paste0("[load_linelist] %d alert_id value(s) are shared by %d records: records sharing ",
                           "an id share one earliest-sample lookup and, if imputed, one onset-imputation ",
                           "uniform."), length(unique(ll$alert_id[.dup])), sum(.dup)), call. = FALSE)
  .residence_basis <- identical(get0("CASE_ZONE_BASIS", ifnotfound = "reporting"), "residence")
  .lookup_records <- if (.residence_basis)
    ll[, intersect(c("health_area", "health_area_unmatched", "health_zone"), names(ll)), drop = FALSE]

  # ---- filter to outbreak start date ----------------------------------------
  # Keep a record if EITHER its symptom onset OR its sample-collection date is on/after
  # OUTBREAK_START. The forecast target is onset-defined, so requiring a NON-NA sample
  # date silently dropped onset-only confirmed cases (e.g. epicentre cases with a blank
  # sample date). Records with both dates pre-outbreak, or neither present, are dropped.
  n_before <- nrow(ll)
  # .record_in_outbreak() (Section 3b) is the ONE statement of this rule, so the records it
  # rejects here and the ones load_unplaced_alerts() returns are complements by construction.
  ll <- ll[.record_in_outbreak(ll$date_of_symptom_onset, ll$date_of_sample_collection), , drop = FALSE]
  message(
    "[load_linelist] After OUTBREAK_START filter (",
    OUTBREAK_START, "): ", nrow(ll), " rows (dropped ", n_before - nrow(ll), ")"
  )

  # ---- append sitrep-confirmed top-up rows (reconcile line list to sitrep) ----
  # Injected AFTER the outbreak filter, unmatched exclusion, and aliasing (so the
  # appends carry canonical, in-spine, in-outbreak, matched zones) and BEFORE the
  # onset-imputation block below (so their blank onset is imputed from their report date, which
  # stands in for a sample date, like any real onset-less record). They carry no recorded onset,
  # so they are never donors; and each record's imputation draw depends only on its own key and
  # the seed, so appending them moves no other record's draw.
  # `aliases` is the table .read_linelist_records() attached: NULL when aliases.csv is absent or
  # unreadable, and the helper then skips sitrep-name canonicalisation.
  .sitrep_appends <- .build_sitrep_confirmed_appends(ll, aliases)
  if (nrow(.sitrep_appends) > 0) ll <- dplyr::bind_rows(ll, .sitrep_appends)
  # The appended rows' draw key is their (unique, generated) alert_id.
  ll$.onset_draw_key <- dplyr::coalesce(ll$.onset_draw_key, as.character(ll$alert_id))
  # RESIDENCE BASIS (sensitivity) AFTER the sitrep top-up: the sitrep counts by reporting zone, so
  # reconciling against residence-basis totals would top up the zones that lose residents with
  # cases the line list already holds elsewhere.
  if (.residence_basis) ll <- .apply_residence_zone(ll, .lookup_records)

  # ---- construct date_index and classification flags -------------------------
  # The forecast target is onset-dated (P of first symptom ONSET). date_index is the recorded
  # onset where it is usable; a record with a sample date but no usable onset gets an onset
  # IMPUTED as sample date minus an interval drawn from DONORS (.impute_onset_delays() below:
  # records of the same class sampled in the same week whose onset IS recorded). Records with
  # neither date get date_index = NA and are dropped downstream.
  # OUTBREAK_START itself, NOT its WEEK_ANCHOR-floored week. WEEK_ANCHOR is derived from the
  # weekday ANALYSIS_DATE falls on (00_config.R), so flooring here made the onset-plausibility
  # floor slide with the day of the week the pipeline was run: on the same frozen snapshot a
  # Monday run floored to 2026-04-28 and a Thursday run to 2026-04-24, flipping onsets in that
  # band between "usable" and "imputed". A plausibility floor is a property of the outbreak, not
  # of the run.
  .ob_floor <- OUTBREAK_START
  # Shared onset->sample plausibility ceiling (00_config.R), the same bound the delay fit uses.
  .max_plausible <- get0("DELAY_MAX_PLAUSIBLE_DAYS", ifnotfound = 60L)
  # Clerical tolerance for an onset recorded AFTER its own specimen (00_config.R).
  .neg_tol <- get0("ONSET_SAMPLE_NEG_TOL_DAYS", ifnotfound = 2L)
  .asof_date <- suppressWarnings(as.Date(get0("ANALYSIS_DATE", ifnotfound = NA)))
  .onset_mode <- get0("ONSET_MODE", ifnotfound = "impute")
  if (!.onset_mode %in% c("impute", "complete_case"))
    stop("[load_linelist] ONSET_MODE must be 'impute' or 'complete_case', got '", .onset_mode, "'.",
         call. = FALSE)
  ll <- ll %>%
    dplyr::mutate(
      # A recorded onset is USABLE only if it is epidemiologically plausible: not before the
      # outbreak start (data-entry YEAR TYPOS put onsets years early), not more than
      # ONSET_SAMPLE_NEG_TOL_DAYS after its own sample, and not more than
      # DELAY_MAX_PLAUSIBLE_DAYS before it. An unusable onset is treated as MISSING and imputed
      # from the sample date, so a typo can neither put a case in a spurious pre-outbreak week
      # nor drop it.
      #
      # ONSET RECORDED AFTER ITS OWN SPECIMEN (within the tolerance): KEPT, CENSORED AT THE
      # SPECIMEN DATE. Such a record is the same-day population displaced by a day or two of
      # transcription noise (on the 2026-09-07 snapshot 39.5% of the -2/-1 d records were
      # deceased at swab, against 22.4% of positive-delay records, so a post-swab symptom onset
      # is not a real reading for them). Rejecting it would impute it ~6 d earlier -- a larger
      # error than the 1-2 d censoring moves it.
      onset_usable = !is.na(date_of_symptom_onset) &
                     date_of_symptom_onset >= .ob_floor &
                     (is.na(date_of_sample_collection) |
                      (date_of_symptom_onset <= date_of_sample_collection + .neg_tol &
                       as.numeric(date_of_sample_collection - date_of_symptom_onset) <= .max_plausible)),
      onset_imputed = !onset_usable & !is.na(date_of_sample_collection),
      onset_censored = onset_usable & !is.na(date_of_sample_collection) &
                       date_of_symptom_onset > date_of_sample_collection,
      # The recorded onset, censored at the specimen date; coalesce() covers a usable onset with
      # no specimen, and pmax holds the invariant the zone-week grid depends on (date_index never
      # precedes the outbreak floor).
      .rec_onset = dplyr::if_else(onset_usable,
        pmax(dplyr::coalesce(pmin(date_of_symptom_onset, date_of_sample_collection),
                             date_of_symptom_onset), .ob_floor),
        as.Date(NA)),
      confirmed = (final_mve_case_classification == CONFIRMED_STATUS) %in% TRUE,
      suspected = (final_mve_case_classification == SUSPECTED_STATUS) %in% TRUE
    )
  .imp_delay <- if (identical(.onset_mode, "impute"))
    .impute_onset_delays(anchor = ll$date_of_sample_collection, recorded = ll$.rec_onset,
                         usable = ll$onset_usable, group = ll$confirmed,
                         need = ll$onset_imputed, ids = ll$.onset_draw_key,
                         seed = get0("ONSET_IMPUTE_SEED",
                                     ifnotfound = get0("RANDOM_SEED", ifnotfound = 20260704L)),
                         asof = .asof_date)
    else rep(NA_integer_, nrow(ll))
  .no_donor <- identical(.onset_mode, "impute") & ll$onset_imputed %in% TRUE & is.na(.imp_delay)
  if (any(.no_donor))
    warning(sprintf(paste0("[load_linelist] %d record(s) needing an imputed onset have no donor of their ",
                           "class (no record of that class with a recorded onset); their date_index is NA ",
                           "and they are dropped downstream."), sum(.no_donor)), call. = FALSE)
  ll$date_index <- dplyr::if_else(
    ll$onset_usable, ll$.rec_onset,
    # complete_case: onset-less records are DROPPED (NA date_index).
    pmax(as.Date(ll$date_of_sample_collection) - .imp_delay, .ob_floor))
  ll$.imp_delay_used <- .imp_delay
  ll$.rec_onset <- NULL
  ll$.onset_draw_key <- NULL
  # AS-OF UPPER BOUND. The rest of this loader guards the as-of date (the sitrep top-up, the
  # donor as-of cut, the earliest-sample lookup), and date_index must too: bounded below by
  # OUTBREAK_START, a record with onset AND sample both mistyped into the future still passes
  # onset_usable and would land in a future week. Such a record is DROPPED, not censored onto
  # the as-of week (a date it cannot have had by then would invent a case in the terminal week):
  # its date_index is set to NA, which is dropped downstream exactly as it is for records with
  # no usable date at all.
  .asof_idx <- suppressWarnings(as.Date(get0("ANALYSIS_DATE", ifnotfound = NA)))
  if (length(.asof_idx) == 1L && !is.na(.asof_idx)) {
    .n_future <- sum(!is.na(ll$date_index) & ll$date_index > .asof_idx, na.rm = TRUE)
    if (.n_future > 0L) {
      warning(sprintf("[load_linelist] %d record(s) have date_index AFTER the analysis date (%s); their date_index is set to NA (dropped downstream). Check for future-dated onset/sample entries.",
                      .n_future, format(.asof_idx)), call. = FALSE)
      ll$date_index[!is.na(ll$date_index) & ll$date_index > .asof_idx] <- as.Date(NA)
    }
  }
  .got_draw  <- ll$onset_imputed %in% TRUE & !is.na(ll$.imp_delay_used)
  .imp_draws <- ll$.imp_delay_used[.got_draw]
  ll$.imp_delay_used <- NULL
  if (identical(.onset_mode, "complete_case")) {
    .n_cc <- sum(is.na(ll$date_index) & !is.na(ll$date_of_sample_collection) &
                 (ll$confirmed %in% TRUE), na.rm = TRUE)
    message(sprintf("[load_linelist] complete-case: dropped %d confirmed records lacking a usable onset (non-imputed analysis)", .n_cc))
  }
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
  # Publish the MECHANISM too, so report prose describes what the code actually did.
  assign(".ONSET_IMPUTE_MODE", .onset_mode, envir = .GlobalEnv)
  assign(".ONSET_IMPUTE_DESC",
         if (identical(.onset_mode, "impute"))
           paste0("the sample date minus a sample-to-onset interval drawn from records of the same ",
                  "class sampled in the same week whose onset is recorded")
         else "no imputation (complete-case: onset-less records dropped)",
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
  # Counts the records that RECEIVED an interval (the no-donor ones are warned about above), split
  # so the confirmed count and its sitrep top-up part are not read as DHIS2 records. Under
  # complete_case nothing is imputed and the complete-case message above says what was dropped.
  if (identical(.onset_mode, "impute"))
    message(sprintf(paste0("[load_linelist] Onset imputed for %d record(s) (%d confirmed, %d of them ",
                           "sitrep top-up rows) via %s%s; onset-dated target."),
                    sum(.got_draw), sum(.got_draw & .conf), sum(.got_draw & .conf & .is_sitrep),
                    get(".ONSET_IMPUTE_DESC", envir = .GlobalEnv),
                    if (length(.imp_draws))
                      sprintf("; drawn interval mean %.1f d, range %d-%d d", mean(.imp_draws),
                              as.integer(min(.imp_draws)), as.integer(max(.imp_draws)))
                    else ""))

  message(
    "[load_linelist] Confirmed: ", sum(ll$confirmed, na.rm = TRUE),
    "  Suspected: ", sum(ll$suspected, na.rm = TRUE),
    "  Zones represented: ", length(unique(ll$health_zone[!is.na(ll$health_zone)]))
  )

  ll
}


# -----------------------------------------------------------------------------
# Archived snapshots: the line list AS IT STOOD on an earlier date
# -----------------------------------------------------------------------------
# Two consumers replay archive vintages through the full loader: the extract-regime truncation
# fit (02_epi_params.R, .trunc_frame_vintage()), which needs each vintage's date_index series,
# and the back-test's real-snapshot information state (LFO_TRAINING_SOURCE = "vintage",
# 22_daily_reissue.R). Defined here, next to load_linelist(), because the truncation fit runs
# before 22_daily_reissue.R is sourced (run_all.R) and in entry points that never source it.

# Session memo: the snapshot listing and every loaded snapshot. parent = emptyenv() so a typo'd
# key cannot resolve up the search path. Filled in the parent process BEFORE the LFO forks its
# workers (prefetch_vintages()), so the workers inherit it instead of each re-reading the CSVs.
.vintage_memo <- new.env(parent = emptyenv())

#' The archived line-list snapshots, oldest first.
#' @return data.frame(folder, date); zero rows when the archive is absent.
list_linelist_snapshots <- function(dir = get0("LINELIST_DIR", ifnotfound = NULL)) {
  empty <- data.frame(folder = character(0), date = as.Date(character(0)),
                      stringsAsFactors = FALSE)
  if (is.null(dir) || !dir.exists(dir)) return(empty)
  f <- list.files(dir, pattern = "^LINELIST_[0-9]{8}$")
  d <- as.Date(sub("^LINELIST_", "", f), format = "%d%m%Y")
  ok <- !is.na(d) & file.exists(file.path(dir, f, "dhis2_processed_linelist.csv"))
  out <- data.frame(folder = f[ok], date = d[ok], stringsAsFactors = FALSE)
  out[order(out$date), , drop = FALSE]
}

#' The settings a replayed snapshot depends on beyond its folder and date: the seed of the onset
#' draw, the onset mode, the onset-plausibility bounds, the donor rule and the week grid it bins
#' donors on (WEEK_ANCHOR), the sample-date rule, and the data-definition switches (case-zone
#' basis, sitrep top-up on/off and its date shift). Part of every session-memo key built from
#' snapshots (here, vintage_truncation() in 22_daily_reissue.R and .asof_counts_wide() in
#' 20_forecast_detail.R).
.vintage_settings_key <- function() {
  substr(rlang::hash(list(
    get0("RANDOM_SEED", ifnotfound = NA), get0("ONSET_IMPUTE_SEED", ifnotfound = NA),
    get0("ONSET_MODE", ifnotfound = NULL),
    get0("ONSET_SAMPLE_NEG_TOL_DAYS", ifnotfound = NULL), get0("DELAY_MAX_PLAUSIBLE_DAYS", ifnotfound = NULL),
    get0("ONSET_DONOR_MIN", ifnotfound = NULL), get0("WEEK_ANCHOR", ifnotfound = NULL),
    get0("SAMPLE_DATE_SOURCE", ifnotfound = NULL),
    get0("CASE_ZONE_BASIS", ifnotfound = NULL), get0("APPEND_SITREP_CONFIRMED", ifnotfound = NULL),
    get0("SITREP_APPEND_SHIFT_DAYS", ifnotfound = NULL))),
    1L, 12L)
}

#' Load one archived snapshot through the FULL load_linelist() pipeline, as of its own date.
#'
#' Every as-of rule inside load_linelist() -- the sitrep top-up, the donor as-of cut of the onset
#' imputation, the date_index ceiling and the earliest-sample-date lookup -- reads the global
#' ANALYSIS_DATE. It is set to the snapshot's date for the duration and restored on exit, so the
#' snapshot is processed as a run on that day would have processed it, with today's code, today's
#' reference files (aliases, zone spine, the sitrep series cut at the date) and today's shared
#' nuisance fits -- and ONE further difference: WEEK_ANCHOR, derived once from the live
#' ANALYSIS_DATE (00_config.R), is not re-derived, so the donor weeks are binned on the live run's
#' weekday grid, the grid the back-test aggregates the snapshot's counts on. load_linelist() uses
#' no RNG stream.
#'
#' Memoised on the folder, the date and .vintage_settings_key(), so a snapshot loaded under one
#' setting is never served under another.
load_linelist_vintage <- function(folder, vdate) {
  key <- paste0("ll_", folder, "_", format(as.Date(vdate)), "_", .vintage_settings_key())
  hit <- .vintage_memo[[key]]
  if (!is.null(hit)) return(hit)
  # load_linelist() publishes the RUN's onset-imputation summary as globals (read by
  # write_run_info.R and the report). A replayed snapshot must not overwrite them, so they are
  # saved and restored with ANALYSIS_DATE.
  .globals <- c("ANALYSIS_DATE", ".PCT_ONSET_IMPUTED", ".PCT_ONSET_IMPUTED_LINELIST",
                ".N_SITREP_APPENDED", ".ONSET_IMPUTE_MODE", ".ONSET_IMPUTE_DESC")
  .saved <- lapply(stats::setNames(.globals, .globals), function(g)
    if (exists(g, envir = .GlobalEnv, inherits = FALSE)) list(get(g, envir = .GlobalEnv)) else NULL)
  on.exit(for (g in .globals) {
    if (is.null(.saved[[g]])) { if (exists(g, envir = .GlobalEnv, inherits = FALSE))
                                  rm(list = g, envir = .GlobalEnv) }
    else assign(g, .saved[[g]][[1]], envir = .GlobalEnv)
  }, add = TRUE)
  assign("ANALYSIS_DATE", as.Date(vdate), envir = .GlobalEnv)
  warns <- character(0)
  ll <- withCallingHandlers(suppressMessages(load_linelist(folder)),
          warning = function(w) { warns <<- c(warns, conditionMessage(w))
                                  invokeRestart("muffleWarning") })
  if (length(warns))
    message(sprintf("[vintage] %s: %d warning(s) while loading, first: %s", folder,
                    length(warns), warns[1]))
  attr(ll, "vintage_date")   <- as.Date(vdate)
  attr(ll, "vintage_folder") <- folder
  assign(key, ll, envir = .vintage_memo)
  ll
}


# =============================================================================
# SECTION 3: aggregate_to_zone_week()
# =============================================================================
#
# Aggregates the linelist to a zone × week matrix (7-day buckets, WEEK_ANCHOR-anchored) and zero-fills
# every cell of (all `zones`) × (the contiguous week grid up to the as-of week) that has no
# observations, so downstream models see a complete rectangular data frame.
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
  # TERMINAL WEEK ANCHORED TO THE AS-OF DATE, not to whatever the data happen to contain.
  # Two silent failure modes otherwise, both one bad row away:
  #   (a) OVER-EXTENSION. load_linelist() sets a date_index after ANALYSIS_DATE to NA, but this
  #       function takes any line list (reaggregate_asof() passes one censored to a fold cutoff,
  #       with its own `asof`), and onset_usable only requires onset <= sample + 2: one record
  #       with onset AND sample both after the as-of date would run the grid past it — adding
  #       structurally-zero cells and, worse, moving the grid's terminal week off the as-of
  #       week the LFO cutoff and the nowcast both assume.
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
  .n_tail <- sum(all_weeks > max(.obs_weeks))                 # padded up to the as-of week
  .n_gap  <- length(all_weeks) - length(.obs_weeks) - .n_tail  # interior gaps

  message(
    "[aggregate_to_zone_week] Observed weeks: ",
    min(.obs_weeks), " to ", max(.obs_weeks),
    " (", length(.obs_weeks), " observed",
    if (.n_gap > 0L) paste0("; +", .n_gap, " zero-filled interior week(s) for contiguity") else "",
    if (.n_tail > 0L) paste0("; +", .n_tail, " zero-filled week(s) after the last observed one, up to the as-of week")
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
# SECTION 3b: .record_in_outbreak(), load_unplaced_alerts(), zone_alert_counts()
# =============================================================================
#
# Per-zone ALERT and LABORATORY-OUTCOME counts. ONE definition, shared by the ascertainment
# appendix (45_ascertainment_figures.R -> zone_ascertainment.csv) and the manuscript
# watch-list figure (make_manuscript_figures.R -> Figure3_rank_surveillance, panel D), so
# the two cannot publish different counts for the same zone.
#
#   alerts_total     alert records in the window (records, not de-duplicated people)
#   n_confirmed      final classification == CONFIRMED_STATUS
#   n_not_a_case     final classification == NOT_A_CASE_STATUS (laboratory-discarded)
#   alerts_resolved  n_confirmed + n_not_a_case: alerts with a laboratory OUTCOME. Suspected,
#                    probable and never-classified alerts are in neither term, including
#                    alerts that were sampled but have no recorded outcome.
#
# WHICH RECORDS ARE IN THE WINDOW. An alert is a line-list record, and it is in [from, to]
# when EITHER of its two kinds of date puts it there:
#   * its CASE date — date_index, the onset (recorded, censored at the specimen, or imputed
#     from it) that every case count uses; only records load_linelist() keeps have one;
#   * its NOTIFICATION date — when the alert was raised.
# Either, because each field carries errors the other does not. On LINELIST_23092026 (as of its
# own date, earliest-sample rule; measured 2026-10-17) 66 of the records load_linelist() keeps,
# ten of them confirmed cases, have a notification date outside [OUTBREAK_START, the extract
# date] — 40 of them 2026 dates before the outbreak start or after the extract date, 26 in
# another year (2006, 2020-2025, 2027) — and a notification-only window would drop them. And
# load_linelist() keeps only records with an onset or specimen date on or after OUTBREAK_START,
# because a CASE needs one: on the same snapshot it leaves out 2,828 records — 2,725 with
# neither date and 103 whose only dates precede the outbreak (several impossible: years 0202,
# 1997, 2006) — of which 2,819 were notified inside the window. All 2,828 are unclassified.
# They are still alerts raised during the outbreak; a case-date-only count missed them and so
# overstated the share resolved. They are read by load_unplaced_alerts() and placed by
# notification.
#
# The sitrep top-up rows (alert_id "SITREP-CONF-*", .build_sitrep_confirmed_appends) are
# EXCLUDED: they reconcile a confirmed COUNT to the official cumulative and carry no alert,
# no test and no laboratory date, so counting them would invent laboratory activity in
# exactly the zones that received a top-up. Records with no zone are dropped.
#
# A zone with no alert in the window has NO ROW. Callers left-join and coalesce to 0,
# because a zone absent here is a zone with zero alerts, not a zone with missing data.
#
# `unplaced`, `from` and `to` are required, not read from globals: 45 and the manuscript
# figure must be seen at the call site to count the same records on the same window.

#' TRUE for a record load_linelist() keeps as part of the outbreak: an onset OR a specimen
#' date on or after OUTBREAK_START. NA-safe; vectorised over the two date vectors.
.record_in_outbreak <- function(onset, sample) {
  (!is.na(onset)  & onset  >= OUTBREAK_START) |
  (!is.na(sample) & sample >= OUTBREAK_START)
}

#' The alert records load_linelist() leaves out — no onset or specimen date on or after
#' OUTBREAK_START — after the same unmatched-zone exclusion and aliasing as the case data.
load_unplaced_alerts <- function() {
  r <- .read_linelist_records()
  attr(r, "aliases") <- NULL
  r$.onset_draw_key <- NULL
  r[!.record_in_outbreak(r$date_of_symptom_onset, r$date_of_sample_collection), , drop = FALSE]
}

#' @param ll       load_linelist() output.
#' @param unplaced load_unplaced_alerts() output.
#' @param from,to  the window, inclusive at both ends.
#' @return one row per zone with an alert; attributes `n_line_list` and `n_unplaced` count the
#'   records of each source inside the window.
zone_alert_counts <- function(ll, unplaced, from, to) {
  cls <- "final_mve_case_classification"
  .check_cols(ll, c("health_zone", "date_index", "date_of_notification", "alert_id", cls),
              "zone_alert_counts")
  .check_cols(unplaced, c("health_zone", "date_of_notification", "date_of_symptom_onset",
                          "date_of_sample_collection", "alert_id", cls), "zone_alert_counts")
  from <- suppressWarnings(as.Date(from)); to <- suppressWarnings(as.Date(to))
  if (length(from) != 1L || is.na(from) || length(to) != 1L || is.na(to) || to < from)
    stop("[zone_alert_counts] `from` and `to` must be single dates with from <= to.",
         call. = FALSE)
  # The two inputs must be the two halves the definition names, or records are double-counted:
  # a record load_linelist() keeps must not also arrive as "unplaced".
  if (any(.record_in_outbreak(unplaced$date_of_symptom_onset, unplaced$date_of_sample_collection)))
    stop("[zone_alert_counts] `unplaced` holds records load_linelist() keeps; pass ",
         "load_unplaced_alerts().", call. = FALSE)
  .in <- function(d) !is.na(d) & d >= from & d <= to
  not_sitrep <- function(id) is.na(id) | !grepl("^SITREP-CONF-", id)
  a <- ll %>%
    dplyr::filter(not_sitrep(alert_id), !is.na(health_zone),
                  .in(date_index) | .in(date_of_notification)) %>%
    dplyr::select(health_zone, dplyr::all_of(cls))
  u <- unplaced %>%
    dplyr::filter(not_sitrep(alert_id), !is.na(health_zone), .in(date_of_notification)) %>%
    dplyr::select(health_zone, dplyr::all_of(cls))
  out <- dplyr::bind_rows(a, u) %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(
      alerts_total = dplyr::n(),
      n_confirmed  = sum(.data[[cls]] %in% CONFIRMED_STATUS),
      n_not_a_case = sum(.data[[cls]] %in% NOT_A_CASE_STATUS),
      .groups = "drop") %>%
    dplyr::mutate(alerts_resolved = n_confirmed + n_not_a_case)
  attr(out, "n_line_list") <- nrow(a)
  attr(out, "n_unplaced")  <- nrow(u)
  out
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
# Reads the following static covariate files and left-joins them on `nom` onto the WorldPop
# zone spine:
#   - CCVI socioeconomic deprivation
#   - GRID3 healthsite count and density
# and derives the model-facing columns log_pop and ccvi.
#
# Returns a tibble with `nom` as the join key, one row per health zone.
# Missing data for a given zone is allowed (NA in the relevant column).

load_static_covariates <- function() {

  # Helper: read a static covariate CSV that has a `nom` column and the value column.
  # Returns a tibble with nom + the value column (numeric); any other column (e.g. a row index)
  # is dropped.
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

  # GDP per capita, PCR testing capacity and cross-border passenger volumes were loaded and
  # joined here until 2026-09-27 and consumed by NOTHING — no covariate set, no vulnerability
  # pillar, no figure. Each also warned when its file was absent, so a clean run printed
  # warnings about data it did not use. Removed; the loaders are one read_static() call away
  # if a model ever wants them.

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
  cov <- join_one(cov, healthsites_count, "healthsite_count")
  cov <- join_one(cov, healthsites_density, "healthsite_density")

  # ---- Derived, model-facing covariate columns ------------------------------
  # The canonical names the covariate sets use (BAYES_GEO_COVARIATES /
  # BAYES_FULL_COVARIATES, 00_config.R) and the vulnerability index reads
  # (compute_vulnerability_index, 20_forecast_detail.R). Deriving them here keeps one source
  # of truth and stops a model silently falling back to a constant covariate.
  # Each is guarded: read_static()/join_one() deliberately warn-and-skip a missing covariate
  # file, so the source column may not exist.
  # `log_healthsite_count` was derived here and read by nothing; removed 2026-09-27.
  cov <- cov %>%
    dplyr::mutate(
      log_pop = if ("pop_count" %in% names(cov)) log(pmax(pop_count, 1)) else NA_real_,
      ccvi    = if ("socioeconomic_deprivation" %in% names(cov))
                  socioeconomic_deprivation else NA_real_
    )

  message("[load_static_covariates] Covariate tibble: ",
          nrow(cov), " zones × ", ncol(cov) - 1, " covariates (excl. nom)")

  cov
}


# =============================================================================
# SECTION 6: prep_all_data() — master loader
# =============================================================================
# load_sitrep() lived here and was removed on 2026-09-27. It was called only by
# prep_all_data(), and its result (dat$sitrep) was assigned in run_all.R and then read by
# nothing. It also read insp_sitrep__new_confirmed_cases__daily.csv — the DAILY file that
# .build_sitrep_confirmed_appends() documents as grossly undercounting (Bunia 66 against
# 507) — so reviving it as-is would be wrong anyway. The sitrep readers that ARE used are
# sitrep_cumulative_confirmed() above (the reconciliation, and FigureS11).
#
# Calls all loaders in order and returns a named list for every downstream script. The zone
# list for zero-filling comes from WorldPop (519 DRC health zones), so the observation matrix
# is always complete.

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

  # A per-zone `positivity` covariate was derived here and joined onto `covariates`. DHIS2
  # carries no test counts for real records, so aggregate_to_zone_week()'s own guard makes
  # `positivity` NA on every row — the join added an all-NA column. It was also not in any
  # covariate set: 20_forecast_detail.R removed it from the default on 2026-09-17 as
  # circular and a full-data static leak. Removed 2026-09-27.

  # ---- Summary statistics ---------------------------------------------------
  message("\n", strrep("-", 70))
  message("[prep_all_data] SUMMARY STATISTICS")
  message(strrep("-", 70))

  confirmed_ll  <- sum(ll$confirmed, na.rm = TRUE)
  .sitrep_rows  <- grepl("^SITREP-CONF-", as.character(ll$alert_id))
  confirmed_sit <- sum(ll$confirmed & .sitrep_rows, na.rm = TRUE)
  zones_affected <- ll %>%
    dplyr::filter(confirmed) %>%
    dplyr::pull(health_zone) %>%
    unique() %>%
    length()
  date_range <- range(ll$date_index, na.rm = TRUE)

  message(
    "  Total confirmed cases:            ", confirmed_ll,
    " (line list ", confirmed_ll - confirmed_sit, " + sitrep top-up ", confirmed_sit, ")"
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
  message(strrep("=", 70), "\n")

  # `zones_all` (a duplicate of names(pop)) and `sitrep` (never read) were dropped from this
  # list on 2026-09-27. Consumers take the canonical zone vector from names(dat$pop).
  list(
    ll         = ll,
    zone_week  = zone_week,
    pop        = pop,
    covariates = covariates
  )
}

# This file is a LIBRARY: sourcing it defines loaders and nothing else. It used to call
# prep_all_data() on source() in an interactive session, which is the same "source() does
# work" hazard the executable scripts are guarded against. Call prep_all_data() yourself.
