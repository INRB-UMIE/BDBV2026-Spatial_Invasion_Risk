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
#     dat$contacts    — cleaned contact-tracing tibble
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

# Load the DHIS2 onset->sample delay params written by 04c_dhis2_delay_windows.R
# (windowed interval-censored MLE). Returns list(family, rate, mean, n_fit, window,
# params = named numeric of the family's native parameters) or NULL if unavailable.
# The `rate` is the Exponential-rate SUMMARY (1/mean) of the AIC-best family; the
# native params (param_shape, param_rate, ...) let the parametric fallback draw from
# the actual best-fit family (e.g. gamma) rather than an Exponential.
.load_dhis2_delay_params <- function() {
  f <- file.path(get0("DATA_DIR", ifnotfound = file.path(here::here(), "data")),
                 "cfr_reference", "dhis2_onset_sample_delay_params.csv")
  if (!file.exists(f)) return(NULL)
  tab <- tryCatch(readr::read_csv(f, col_types = readr::cols(.default = "c"),
                                  show_col_types = FALSE), error = function(e) NULL)
  if (is.null(tab) || !all(c("quantity", "value") %in% names(tab))) return(NULL)
  g   <- function(q) { v <- tab$value[tab$quantity == q]; if (length(v)) v[1] else NA_character_ }
  num <- function(q) suppressWarnings(as.numeric(g(q)))
  # PREFER the truncation-corrected EpiDist MARGINAL estimate when 04c wrote it
  # (RUN_EPIDIST=TRUE): it corrects right-truncation AND double-interval censoring, whereas
  # the windowed interval-censored MLE only mitigates truncation by dropping the recent tail
  # (review §1.2). Derive the family's native params from the marginal mean/SD (method of
  # moments for gamma; log-moments for lnorm) so the parametric fallback draw and the reported
  # rate both use the corrected distribution. Fall back to the MLE family/rate otherwise.
  epi_fam <- g("epidist_family"); epi_mean <- num("epidist_mean"); epi_sd <- num("epidist_sd")
  if (!is.na(epi_fam) && is.finite(epi_mean) && epi_mean > 0) {
    ep <- numeric(0)
    if (identical(epi_fam, "gamma") && is.finite(epi_sd) && epi_sd > 0)
      ep <- c(shape = (epi_mean / epi_sd)^2, rate = epi_mean / epi_sd^2)
    else if (identical(epi_fam, "lnorm") && is.finite(epi_sd) && epi_sd > 0) {
      .s2 <- log(1 + (epi_sd / epi_mean)^2); ep <- c(meanlog = log(epi_mean) - .s2 / 2, sdlog = sqrt(.s2)) }
    return(list(family = epi_fam, rate = 1 / epi_mean, mean = epi_mean, sd = epi_sd,
                n_fit = num("epidist_n"), window = "epidist_marginal (truncation-corrected)",
                params = ep[is.finite(ep)], estimator = "epidist_marginal"))
  }
  fam <- g("family"); rate <- num("rate")
  if (is.na(fam) || !is.finite(rate) || rate <= 0) return(NULL)
  par_q  <- tab$quantity[grepl("^param_", tab$quantity)]
  params <- suppressWarnings(setNames(as.numeric(tab$value[match(par_q, tab$quantity)]),
                                      sub("^param_", "", par_q)))
  list(family = fam, rate = rate, mean = num("implied_mean_fit"),
       n_fit = num("n_fit"), window = g("window"), params = params[is.finite(params)],
       estimator = "interval_censored_mle")
}

# Draw n onset->sample delays from a loaded DHIS2 best-family fit (used only for the
# parametric FALLBACK when <30 windowed complete pairs exist; the empirical bootstrap
# is preferred otherwise). Falls back to Exp(rate) if the native params are absent.
.draw_dhis2_delay <- function(dp, n) {
  p  <- dp$params
  ok <- function(nm) all(nm %in% names(p)) && all(is.finite(p[nm]))
  switch(as.character(dp$family),
    gamma   = if (ok(c("shape", "rate")))    stats::rgamma(n,   shape = p[["shape"]],   rate  = p[["rate"]])    else stats::rexp(n, dp$rate),
    weibull = if (ok(c("shape", "scale")))   stats::rweibull(n, shape = p[["shape"]],   scale = p[["scale"]])   else stats::rexp(n, dp$rate),
    lnorm   = if (ok(c("meanlog", "sdlog"))) stats::rlnorm(n,   meanlog = p[["meanlog"]], sdlog = p[["sdlog"]]) else stats::rexp(n, dp$rate),
    stats::rexp(n, dp$rate))
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
.build_sitrep_confirmed_appends <- function(ll, aliases_all) {

  if (!isTRUE(get0("APPEND_SITREP_CONFIRMED", ifnotfound = TRUE))) return(tibble::tibble())

  cum_path <- file.path(SITREP_DIR, "insp_sitrep__cumulative_confirmed_cases__daily.csv")
  if (!file.exists(cum_path)) {
    warning("[load_linelist] sitrep cumulative-confirmed file not found — no sitrep reconciliation: ",
            cum_path, call. = FALSE)
    return(tibble::tibble())
  }
  cum <- tryCatch(
    readr::read_csv(cum_path,
      col_types = readr::cols(nom = "c", date = "c", cumulative_confirmed_cases = "c"),
      show_col_types = FALSE),
    error = function(e) { warning("[load_linelist] sitrep read error: ", e$message, call. = FALSE); NULL })
  if (is.null(cum) || !all(c("nom", "date", "cumulative_confirmed_cases") %in% names(cum)))
    return(tibble::tibble())

  # Canonical modelling spine (WorldPop 519 zones): a sitrep zone absent from it cannot
  # enter the zone-week grid, so we skip (and warn about) it rather than fabricate an
  # unmatched row that aggregate_to_zone_week would silently drop.
  spine <- tryCatch(
    readr::read_csv(file.path(WORLDPOP_DIR, "worldpop__pop_count__static.csv"),
      col_types = readr::cols(nom = "c", pop_count = "d"), show_col_types = FALSE)$nom,
    error = function(e) NULL)

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
  # over dates; then diff into dated positive integer increments (one row per case).
  incs <- cum %>%
    dplyr::group_by(nom, date) %>%
    dplyr::summarise(cum = max(cum), .groups = "drop") %>%
    dplyr::arrange(nom, date) %>%
    dplyr::group_by(nom) %>%
    dplyr::mutate(inc = as.integer(round(cummax(cum) - dplyr::lag(cummax(cum), default = 0)))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(inc > 0L) %>%
    dplyr::select(nom, date, inc)
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
  mk <- function(z, short) {
    d <- pick_dates(z, short)
    tibble::tibble(
      country                       = "Democratic Republic of the Congo",
      health_zone                   = z,
      health_zone_unmatched         = FALSE,
      health_area_unmatched         = FALSE,
      date_of_sample_collection     = as.Date(d),
      date_of_symptom_onset         = as.Date(NA),   # imputed downstream from the delay dist
      final_mve_case_classification = CONFIRMED_STATUS,
      genexpert_result              = "positive",
      samples_received              = 1,
      samples_analyzed              = 1,
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
          dplyr::filter(source_dataset %in% c("linelist", "epi", "dhis2") |
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
  .ob_floor <- lubridate::floor_date(OUTBREAK_START, unit = "week",
                                     week_start = get0("WEEK_ANCHOR", ifnotfound = 1L))
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
  # Rigorous DHIS2 delay params (windowed interval-censored MLE) from
  # 04c_dhis2_delay_windows.R, when it has been run and the source is "data": its AIC-best
  # family + Exponential-rate summary supersede the crude 1/mean for the reported rate and the
  # parametric-fallback draw (a gamma/weibull/lnorm fit rather than an Exponential).
  .dhis2_delay <- if (identical(.delay_src, "data")) .load_dhis2_delay_params() else NULL
  .est_rate <- if (!is.null(.dhis2_delay)) .dhis2_delay$rate
               else if (length(.dd) >= 30 && mean(.dd) > 0) 1 / mean(.dd) else NA_real_
  .rate_used <- if (identical(.delay_src, "data") && is.finite(.est_rate)) .est_rate else .fixed_rate
  .imp_active <- isTRUE(get0("IMPUTE_ONSET_FROM_SAMPLE", ifnotfound = TRUE)) &&
                 is.finite(.rate_used) && .rate_used > 0
  # --- Onset-handling MODE (review §1.1/§1.3) --------------------------------------------------
  # ONSET_MODE selects how confirmed cases lacking a usable onset are handled (see 00_config.R):
  # impute (default) / growth_impute (growth-tilted backward draw) / complete_case (drop) /
  # sample_verbatim (date at sample). Unset -> derived from IMPUTE_ONSET_FROM_SAMPLE.
  .onset_mode <- get0("ONSET_MODE", ifnotfound = NA_character_)
  if (length(.onset_mode) != 1L || is.na(.onset_mode))
    .onset_mode <- if (isTRUE(get0("IMPUTE_ONSET_FROM_SAMPLE", ifnotfound = TRUE))) "impute" else "sample_verbatim"
  if (!.onset_mode %in% c("impute", "growth_impute", "complete_case", "sample_verbatim")) {
    warning("[load_linelist] unknown ONSET_MODE '", .onset_mode, "'; using 'impute'"); .onset_mode <- "impute" }
  # National per-day epidemic growth rate r: log-linear slope of recent weekly confirmed counts
  # (by sample date). Used ONLY to growth-tilt the backward delay draw under "growth_impute":
  # weighting Delta by exp(-r*Delta) down-weights long delays while the epidemic grows, so imputed
  # onsets are not pushed systematically too early (epidemic processes are not time-reversible).
  .estimate_growth_rate <- function(dates, asof, window_weeks = get0("ONSET_GROWTH_WINDOW_WEEKS", ifnotfound = 8L)) {
    d <- dates[!is.na(dates) & (is.na(asof) | dates <= asof)]
    if (length(d) < 20L) return(NA_real_)
    wk  <- as.integer(floor(as.numeric(d - min(d)) / 7))
    tab <- table(wk); wks <- as.integer(names(tab)); cnt <- as.numeric(tab)
    keep <- wks >= (max(wks) - window_weeks); wks <- wks[keep]; cnt <- cnt[keep]
    if (length(wks) < 3L || sum(cnt) < 10) return(NA_real_)
    fit <- tryCatch(stats::lm(log(cnt + 0.5) ~ wks), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    unname(coef(fit)[2]) / 7   # per-week slope -> per-day growth rate
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
  # Exp(rate_used). Delays clamped to [0, 90] d. Seeded for reproducibility.
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
  .imp_mode <- if (.onset_mode == "sample_verbatim" || !.imp_active) "none"
               else if (.onset_mode == "complete_case") "complete_case"
               else if (.onset_mode == "growth_impute" && length(.dd) >= 30L) "growth"
               else if (length(.dd) >= 30L) "empirical"
               else "parametric"
  # Parametric fallback draw (only hit when <30 complete pairs): prefer the rigorous DHIS2
  # best-family fit (e.g. gamma), else Exp(rate_used).
  .param_draw <- if (!is.null(.dhis2_delay)) function(n) .draw_dhis2_delay(.dhis2_delay, n)
                 else function(n) stats::rexp(n, rate = .rate_used)
  ll$.imp_delay <- switch(.imp_mode,
    empirical  = as.integer(pmin(pmax(round(sample(.dd, nrow(ll), replace = TRUE)), 0L), .max_plausible)),
    # growth_impute: same empirical bootstrap but with growth-tilt weights exp(-r*Delta) (§1.1).
    growth     = as.integer(pmin(pmax(round(sample(.dd, nrow(ll), replace = TRUE, prob = .dd_wts)), 0L), .max_plausible)),
    parametric = as.integer(pmin(pmax(round(.param_draw(nrow(ll))), 0L), .max_plausible)),
    # complete_case: onset-less confirmed records are DROPPED in the mutate below (delay unused).
    complete_case = rep(0L, nrow(ll)),
    none       = rep(0L, nrow(ll)))   # no delay info at all: degenerate (onset=sample); warned below
  if (identical(.imp_mode, "none"))
    warning("[load_linelist] Onset imputation requested but no onset->sample delay is available ",
            "(no fixed rate and <30 complete pairs); imputed onsets fall back to the sample date.")
  ll <- ll %>%
    dplyr::mutate(
      # A recorded onset is USABLE only if it is epidemiologically plausible: not before
      # the outbreak week (data-entry YEAR TYPOS put onsets years early — onset 2020/2023
      # with a 2026 sample), not after its own sample, and not implausibly long before it
      # (> the shared DELAY_MAX_PLAUSIBLE_DAYS ceiling — the same bound the delay fit treats
      # as an outlier). An unusable onset is treated as MISSING and imputed from the sample
      # date, so a typo can neither leak a spurious pre-outbreak week into the zone-week grid
      # (the CRITICAL corruption) nor drop the real case.
      onset_usable = !is.na(date_of_symptom_onset) &
                     date_of_symptom_onset >= .ob_floor &
                     (is.na(date_of_sample_collection) |
                      (date_of_symptom_onset <= date_of_sample_collection + 2L &
                       as.numeric(date_of_sample_collection - date_of_symptom_onset) <= .max_plausible)),
      onset_imputed = !onset_usable & !is.na(date_of_sample_collection),
      date_index = dplyr::if_else(
        onset_usable,
        date_of_symptom_onset,
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
  .imp_draws <- ll$.imp_delay[ll$onset_imputed %in% TRUE]
  ll <- ll %>% dplyr::select(-dplyr::any_of(".imp_delay"))
  if (identical(.onset_mode, "complete_case")) {
    .n_cc <- sum(is.na(ll$date_index) & !is.na(ll$date_of_sample_collection) &
                 (ll$confirmed %in% TRUE), na.rm = TRUE)
    message(sprintf("[load_linelist] complete-case (§1.3): dropped %d confirmed records lacking a usable onset (non-imputed analysis)", .n_cc))
  }
  .param_desc <- if (!is.null(.dhis2_delay))
    sprintf("the DHIS2 interval-censored %s fit (rate %.3f/d, window %s)",
            .dhis2_delay$family, .rate_used, .dhis2_delay$window %||% "n/a")
    else sprintf("Exp(rate %.3f/d, source %s)", .rate_used, .delay_src)
  message("[load_linelist] Onset imputed for ", sum(ll$onset_imputed, na.rm = TRUE),
          " records via ",
          switch(.imp_mode,
                 empirical  = sprintf("a draw from the EMPIRICAL onset->sample delay (%d %s pairs, source %s)",
                                      length(.dd), if (.windowed) "windowed" else "all-complete", .delay_src),
                 parametric = sprintf("a parametric draw from %s", .param_desc),
                 none       = "the sample date (no delay available)"),
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

aggregate_to_zone_week <- function(ll, zones) {

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
      positivity = dplyr::if_else(
        tests_analyzed > 0,
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
  all_weeks <- seq(min(.obs_weeks), max(.obs_weeks), by = 7L)

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
    message("[load_static_covariates] Joining ", label,
            " — ", n_matched, "/", nrow(df), " zones matched")
    dplyr::left_join(base, df, by = "nom")
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

  sitrep <- sitrep %>%
    dplyr::mutate(
      date       = .parse_date(date),
      week_start = lubridate::floor_date(date, unit = "week",
                                         week_start = get0("WEEK_ANCHOR", ifnotfound = 1L))
    ) %>%
    dplyr::filter(!is.na(date), !is.na(nom), nom != "") %>%
    dplyr::group_by(nom, week_start) %>%
    dplyr::summarise(
      sitrep_confirmed_weekly = sum(new_confirmed_cases, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(nom, week_start) %>%
    dplyr::rename(health_zone = nom)  # align key name with zone_week

  message(
    "[load_sitrep] Sitrep: ", length(unique(sitrep$health_zone)), " zones, ",
    nrow(sitrep), " zone-week rows, weeks ",
    min(sitrep$week_start), " to ", max(sitrep$week_start)
  )

  sitrep
}


# =============================================================================
# SECTION 7: load_contact_tracing()
# =============================================================================
#
# Reads the processed contact-tracing CSV identified by the contact-tracing
# latest.json pointer, applies logical coercions to unmatched flags, and
# returns a cleaned tibble. Downstream scripts (02_epi_params.R) extract
# transmission pairs; here we simply load, type-cast, and report.

load_contact_tracing <- function() {

  ct_json <- file.path(
    PROC_DIR, "contact_tracing_processed", "latest.json"
  )

  if (!file.exists(ct_json)) {
    warning(
      "[load_contact_tracing] contact tracing latest.json not found: ", ct_json,
      call. = FALSE
    )
    return(NULL)
  }

  ct_meta <- tryCatch(
    jsonlite::fromJSON(ct_json, simplifyVector = TRUE),
    error = function(e) {
      warning("[load_contact_tracing] Cannot parse latest.json: ", e$message, call. = FALSE)
      return(NULL)
    }
  )

  if (is.null(ct_meta) || !"folder" %in% names(ct_meta)) {
    warning("[load_contact_tracing] latest.json has no 'folder' key.", call. = FALSE)
    return(NULL)
  }

  ct_csv <- file.path(
    PROC_DIR, "contact_tracing_processed",
    ct_meta$folder, "contact_tracing_processed.csv"
  )

  if (!file.exists(ct_csv)) {
    warning("[load_contact_tracing] Contact tracing CSV not found: ", ct_csv, call. = FALSE)
    return(NULL)
  }

  message("[load_contact_tracing] Reading: ", ct_csv)

  ct <- tryCatch(
    readr::read_csv(
      ct_csv,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE,
      na = c("", "NA", "N/A")
    ),
    error = function(e) {
      warning("[load_contact_tracing] Read error: ", e$message, call. = FALSE)
      return(NULL)
    }
  )

  if (is.null(ct)) return(NULL)

  required_ct_cols <- c(
    "health_zone", "health_zone_unmatched",
    "source_case_alert_id", "contact_id",
    "registration_date", "followup_date", "followup_result"
  )
  .check_cols(ct, required_ct_cols, "load_contact_tracing")

  # ---- coerce types ----------------------------------------------------------
  ct <- ct %>%
    dplyr::mutate(
      health_zone_unmatched = .parse_logical_col(health_zone_unmatched),
      health_area_unmatched =
        if ("health_area_unmatched" %in% names(.))
          .parse_logical_col(health_area_unmatched) else NA,
      is_last_followup_day  =
        if ("is_last_followup_day" %in% names(.))
          .parse_logical_col(is_last_followup_day) else NA,
      registration_date     = .parse_date(registration_date),
      last_contact_with_source_date =
        if ("last_contact_with_source_date" %in% names(.))
          .parse_date(last_contact_with_source_date) else as.Date(NA),
      followup_date         = .parse_date(followup_date),
      contact_age           = suppressWarnings(as.numeric(contact_age))
    )

  n_unmatched <- sum(ct$health_zone_unmatched %in% TRUE, na.rm = TRUE)
  if (n_unmatched > 0) {
    message(
      "[load_contact_tracing] Note: ", n_unmatched,
      " contacts have health_zone_unmatched == TRUE (retained; flag for downstream use)"
    )
  }

  message(
    "[load_contact_tracing] ", nrow(ct), " contact-follow-up records; ",
    length(unique(ct$source_case_alert_id)), " source cases; ",
    length(unique(ct$health_zone)), " zones"
  )

  ct
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

  # 4b. Per-zone mean test positivity (a model covariate; derived from the
  #     zone-week aggregation since it is not a static input file). Zones with
  #     no testing data fall back to the nominal ascertainment rate downstream.
  zone_positivity <- zone_week %>%
    dplyr::group_by(health_zone) %>%
    dplyr::summarise(positivity = mean(positivity, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(positivity = ifelse(is.nan(positivity), NA_real_, positivity))
  covariates <- covariates %>%
    dplyr::left_join(zone_positivity, by = c("nom" = "health_zone"))

  # 5. INSP sitrep (optional)
  sitrep <- load_sitrep()

  # 6. Contact tracing
  contacts <- load_contact_tracing()

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
  message(
    "  Contact tracing records:          ",
    if (!is.null(contacts)) nrow(contacts) else "not loaded"
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
    sitrep     = sitrep,
    contacts   = contacts
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
