# =============================================================================
# 00_config.R — Analysis Configuration
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: single source of truth for all paths, constants, and analysis
#   parameters referenced across scripts 01–16 and the test suite.
# Usage: source("spatiotemporal/00_config.R") at the top of every script.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT           <- here::here()                    # repo root
DATA_DIR       <- file.path(ROOT, "data")
PROC_DIR       <- file.path(DATA_DIR, "processed")
ST_DIR         <- file.path(ROOT, "spatiotemporal")
OUT_DIR        <- file.path(ST_DIR, "outputs")
MODELS_DIR     <- file.path(ST_DIR, "models")

# Data sub-paths (all relative to DATA_DIR)
LINELIST_DIR   <- file.path(PROC_DIR, "dhis2_linelist_processed")
LINELIST_JSON  <- file.path(LINELIST_DIR, "latest.json")

FLOWMINDER_DIR      <- file.path(DATA_DIR, "flowminder", "processed")
FLOWMINDER_ST_DIR   <- file.path(DATA_DIR, "flowminder_short_trips", "processed")
OSRM_DIR            <- file.path(DATA_DIR, "osrm", "processed")
IDP_DIR             <- file.path(DATA_DIR, "IDP", "processed")
WORLDPOP_DIR        <- file.path(DATA_DIR, "worldpop", "processed")
CCVI_DIR            <- file.path(DATA_DIR, "ccvi", "processed")
GDP_DIR             <- file.path(DATA_DIR, "gdp_pc", "processed")
HEALTHSITES_DIR     <- file.path(DATA_DIR, "grid3_healthsites", "processed")
TESTING_DIR         <- file.path(DATA_DIR, "testing_capacity", "processed")
SHAPEFILE_PATH      <- file.path(DATA_DIR, "shapefiles", "DRC_Health_zones.shp")
ALIASES_PATH        <- file.path(DATA_DIR, "aliases.csv")
SITREP_DIR          <- file.path(DATA_DIR, "insp_sitrep", "processed")

# Output sub-directories
OUT_FORECASTS    <- file.path(OUT_DIR, "forecasts")
OUT_MOBILITY     <- file.path(OUT_DIR, "mobility")
OUT_DIAGNOSTICS  <- file.path(OUT_DIR, "diagnostics")
OUT_MAPS         <- file.path(OUT_DIR, "maps")
OUT_CALIBRATION  <- file.path(OUT_DIR, "calibration")
OUT_REPORTS      <- file.path(OUT_DIR, "reports")

# ---------------------------------------------------------------------------
# Linelist constants
# ---------------------------------------------------------------------------
# Primary outcome: confirmed BDBV case classification value
CONFIRMED_STATUS   <- "confirmed_case"
SUSPECTED_STATUS   <- "suspected_case"
PROBABLE_STATUS    <- "probable_case"
NOT_A_CASE_STATUS  <- "not_a_case"

# Current outbreak start: the earliest SAMPLE-COLLECTION date among CONFIRMED cases — one
# record, Rwampara (Ituri), swabbed this day and confirmed on 31 May. Three things this is
# NOT, each of which has been mistaken for it:
#   * not the earliest sample date in the line list, which is 2026-02-23 on an alert that
#     was never classified (2026-04-28 on the earliest `not_a_case`);
#   * not the earliest recorded confirmed ONSET, which is 2026-04-20. Three confirmed
#     records carry a recorded onset before this floor; `onset_usable` rejects them and
#     their onsets are re-imputed from the specimen date (see 01_data_prep.R);
#   * not an estimate of when transmission began. `date_index` is CLAMPED to this floor, so
#     the first week of the onset-dated epidemic curve is set by this constant rather than
#     discovered from the data — on the 2026-09-07 snapshot eleven records sit on the floor,
#     eight of them by clamping. Describe it as the first confirmed case DETECTED.
OUTBREAK_START     <- as.Date("2026-04-30")

# How far an onset may be recorded AFTER its own specimen and still be believed, in days.
# Within this window the onset is CENSORED at the specimen date; beyond it the onset field
# is treated as unrecoverable and imputed (01_data_prep.R, the `date_index` block, which
# carries the evidence for the window).
ONSET_SAMPLE_NEG_TOL_DAYS <- 2L

# Forecast TARGET is P(first symptom ONSET of a confirmed case in the next 1-2
# weeks) — an onset-dated event, since transmission tracks onset, not reporting.
# ~24% of confirmed records lack an onset date (~16% genuine DHIS2 missingness plus the
# sitrep-reconciliation rows, which have none by construction); when TRUE we impute their
# onset as (sample_date - Delta), where Delta is DRAWN PER RECORD from the onset->sample
# delay — the fitted, truncation-corrected delay by default, with the empirical bootstrap
# and a parametric Exp(rate) as fallbacks (see ONSET_SAMPLE_DELAY_SOURCE below for the full
# precedence) — rather than fixed at its mean, so imputed onsets
# contribute a true ONSET week spread by the delay's shape rather than piling on one
# week. Set FALSE to fall back to sample date verbatim.
IMPUTE_ONSET_FROM_SAMPLE <- TRUE

# Onset-handling MODE for confirmed cases lacking a usable onset date (review §1.1/§1.3).
# Supersedes IMPUTE_ONSET_FROM_SAMPLE when set (which is kept for backward compatibility:
# unset ONSET_MODE derives "impute" when IMPUTE_ONSET_FROM_SAMPLE=TRUE else "sample_verbatim").
#   "growth_impute" : DEFAULT. Backward draw TILTED by the epidemic growth rate r.
#   "impute"        : untilted stochastic backward draw onset = sample - Delta.
#   "complete_case" : DROP onset-less confirmed records (a robustness arm; pair with the
#                     nowcast). NOTE: this changes the invasion event set.
#   "sample_verbatim": date the record at its sample date (old IMPUTE_ONSET_FROM_SAMPLE=FALSE path).
#
# WHY THE TILT IS THE DEFAULT, NOT A SENSITIVITY ARM.
# Naive backward imputation draws Delta from the MARGINAL onset->sample delay and sets
# onset = sample - Delta. That is only correct if the epidemic is time-reversible, which it is
# not: conditional on a case being sampled at time s, the delay's density is
#     p(Delta | sampled at s)  proportional to  f(Delta) * I(s - Delta)
# where I is incidence. Under I(t) ~ exp(r*t) this is f(Delta) * exp(-r*Delta) up to a
# constant — the marginal f(Delta) is the special case r = 0. Using the marginal therefore
# mis-dates onsets SYSTEMATICALLY, in the direction of the growth rate's sign: while incidence
# is rising it pushes them too early (over-smoothing the recent rise), and while it is falling
# it pushes them too late. Because ~24% of confirmed cases here carry an imputed onset, and
# because the invasion outcome is defined by the FIRST onset in a zone, that bias lands
# directly on the quantity this pipeline forecasts.
#
# This is the criticism in Lison et al. 2024, PLOS Comput Biol 10.1371/journal.pcbi.1012021
# (which covers imputation and nowcasting together), and the tilt above is exactly the
# growth-adjusted imputation it calls for. The implementation draws a large pool from the
# EpiDist-corrected delay and importance-resamples it with weights exp(-r*Delta), so the tilt
# is applied to the truncation-corrected distribution rather than to a raw empirical one
# (01_data_prep.R, .resolver_draw(tilt = TRUE)).
#
# LIMITS, stated rather than hidden:
#   * r is estimated NATIONALLY (log-linear slope of recent weekly counts by sample date, see
#     ONSET_GROWTH_WINDOW_WEEKS) and applied to every zone. In a spatially heterogeneous
#     outbreak a zone growing against the national trend is still mis-tilted, just far less
#     than it was under r = 0.
#   * Imputation and nowcasting remain TWO STEPS, not one joint generative model. The same
#     paper's other criticism — that stepwise pipelines are biased relative to an integrated
#     model — therefore still applies in part. epinowcast (04b_epinowcast.R) is the integrated
#     alternative and runs here as a sensitivity arm, not the primary.
#   * A single stochastic imputation does not propagate imputation uncertainty; intervals
#     conditional on the draw are slightly narrow. RANDOM_SEED varies it for replicates.
ONSET_MODE <- get0("ONSET_MODE", ifnotfound = "growth_impute")
ONSET_GROWTH_WINDOW_WEEKS <- 8L   # recent-weeks window for the log-linear growth-rate estimate
# Days of the most recent SAMPLE-date series to DROP before fitting that rate. Those days are
# right-truncated (sampled but not yet in the extract), and fitting through them reads the
# reporting lag as an epidemiological decline — which flips the sign of the imputation tilt.
# See .estimate_growth_rate() in 01_data_prep.R for the measured effect.
ONSET_GROWTH_TRUNC_BUFFER_DAYS <- 14L
ONSET_GROWTH_RATE_MAX <- 0.1      # |r| above this disables the tilt (implausible; ~2x/week)

# --- Sitrep reconciliation of the confirmed-case set --------------------------------------------
# When TRUE, load_linelist() reconciles the DHIS2 line list to the INSP sitrep's OFFICIAL
# cumulative confirmed counts before modelling (.build_sitrep_confirmed_appends). The sitrep is
# treated as a FLOOR, never a ceiling: for each canonical zone we append the SHORTFALL =
# max(0, sitrep_cumulative_confirmed - linelist_confirmed) as extra confirmed rows, so every zone
# reaches AT LEAST its officially confirmed count; zones where the line list already meets/exceeds
# the sitrep get nothing (no case is ever removed). This corrects zones the sitrep confirms while
# the DHIS2 line list still holds only suspects (e.g. Oicha, Isiro, Pawa), which would otherwise be
# scored as never-invaded — a ground-truth error for the retrospective model-selection CV run in
# THIS baseline folder (the truly-invaded zones must be scoreable events for a valid comparison).
# Assumptions & mechanics are documented in full at 01_data_prep.R (.build_sitrep_confirmed_appends)
# and in the analysis report. Distinct from build_conditional_linelist.R's speculative prospective
# appends: this is a correction to the OFFICIAL record and is applied in BOTH the baseline and
# conditional pipelines; it is idempotent w.r.t. those appends (they raise the line-list count, so
# the shortfall is 0). Set FALSE to model the raw DHIS2 line list verbatim (e.g. a sitrep-free
# sensitivity run, or to reproduce a model selection made on the unreconciled line list).
APPEND_SITREP_CONFIRMED <- TRUE

# Analysis reference ("as-of") date for right-truncation nowcasting and R(t).
# REAL-TIME: this DERIVES automatically from the processed data so it always tracks
# the latest pull — it reads `processed_at` from the DHIS2 latest.json (the snapshot
# date), then falls back to today's date, then to a fixed date. Set ANALYSIS_DATE in
# the environment before sourcing to override (e.g. for a back-dated re-run).
ANALYSIS_DATE <- local({
  # Explicit override for a back-dated re-run: an ANALYSIS_DATE environment variable wins.
  env <- Sys.getenv("ANALYSIS_DATE", "")
  if (nzchar(env)) { e <- tryCatch(as.Date(env), error = function(e) as.Date(NA))
                     if (!is.na(e)) return(e) }
  d <- NA_character_
  j <- tryCatch(jsonlite::fromJSON(LINELIST_JSON), error = function(e) NULL)
  if (!is.null(j) && !is.null(j$processed_at)) d <- substr(as.character(j$processed_at)[1], 1, 10)
  # as.Date() ERRORS (not warns) on a non-empty unparseable string, so tryCatch (not
  # suppressWarnings) is required for the fallback to actually engage; [1] guards a
  # vector-valued processed_at (JSON array) that would break the scalar is.na() below.
  d <- tryCatch(as.Date(d[1]), error = function(e) as.Date(NA))
  if (is.na(d)) d <- tryCatch(Sys.Date(), error = function(e) as.Date(NA))
  # No frozen fallback date (which would silently mis-date a real-time run): fail loudly so the
  # operator sets ANALYSIS_DATE explicitly. Reaching here needs both no processed_at AND a broken
  # Sys.Date(), i.e. a misconfigured environment.
  if (is.na(d)) stop("[config] Could not derive ANALYSIS_DATE from latest.json 'processed_at' and Sys.Date() is unavailable. Set the ANALYSIS_DATE environment variable (YYYY-MM-DD).")
  d
})

# Weekly grid anchor. The weekly zone-week grid is built as 7-day windows that END
# on the analysis date's weekday, so the FINAL (current) week ends exactly on
# ANALYSIS_DATE and the training window closes on the as-of date rather than on the
# preceding ISO Monday. Weeks stay a uniform 7 days (so the renewal / R(t) / GT-PMF
# machinery is untouched); only the anchor day shifts. Expressed as lubridate's
# `week_start` (1 = Monday … 7 = Sunday): the day AFTER the analysis weekday.
# NOTE: because the anchor tracks the analysis weekday, the grid re-anchors when the
# pipeline is run on a different weekday — a deliberate consequence of ending on the
# as-of date rather than on fixed ISO weeks. Every floor_date(..., "week") in the
# pipeline reads this constant (via get0 with a Monday fallback for config-less unit tests).
WEEK_ANCHOR <- (as.integer(lubridate::wday(ANALYSIS_DATE, week_start = 1L)) %% 7L) + 1L
message(sprintf("[config] Weekly grid anchored to end on the analysis weekday (week_start = %d; weeks end on %s)",
                WEEK_ANCHOR, format(ANALYSIS_DATE, "%A")))

# ---------------------------------------------------------------------------
# Epidemiological parameters (peer-reviewed sources documented in 02_epi_params.R)
# ---------------------------------------------------------------------------

# Onset-to-sample delay: Exponential, fit to BDBV-2026 Ituri LAB data — the FALLBACK
# reference (used when no DHIS2-specific fit is available). Method: interval-censored MLE, n=545.
# The DHIS2 line list reports MORE SLOWLY than the lab (onset->sample mean ~5.9 d, rate ~0.17
# vs 4.39 d / 0.228 here) — 04c_dhis2_delay_windows.R fits it rigorously (windowed interval-
# censored MLE), which the onset imputation prefers when ONSET_SAMPLE_DELAY_SOURCE = "data" (below).
# READ FROM THE FILE, not hard-coded. data/cfr_reference/onset_to_sample_delay_params.csv is
# the lab fit's own record (family=exp, rate, n_fit, window); the literal 0.228 here was a
# mirror of its `rate` row and would silently diverge the moment that file was refitted.
# Falls back to the literal ONLY if the file is unreadable, and says so.
LAB_DELAY_PARAMS_PATH <- file.path(DATA_DIR, "cfr_reference", "onset_to_sample_delay_params.csv")
DELAY_ONSET_SAMPLE_RATE <- local({
  .fallback <- 0.228
  r <- tryCatch({
    if (!file.exists(LAB_DELAY_PARAMS_PATH)) stop("not on disk")
    d <- utils::read.csv(LAB_DELAY_PARAMS_PATH, stringsAsFactors = FALSE)
    v <- suppressWarnings(as.numeric(d$value[d$quantity == "rate"]))[1]
    if (!is.finite(v) || v <= 0) stop("no usable `rate` row")
    v
  }, error = function(e) {
    warning(sprintf(paste0("[config] lab onset->sample reference unreadable (%s): %s. Using the ",
                           "built-in %.3f/d. This value is only ever a LAST-RESORT fallback — the ",
                           "pipeline imputes onsets from the EpiDist fit (see ",
                           "effective_onset_sample_delay below)."),
                    basename(LAB_DELAY_PARAMS_PATH), conditionMessage(e), .fallback), call. = FALSE)
    .fallback
  })
  as.numeric(r)
})

# Right-truncation buffer (days) for delay fitting/imputation windowing: onset->sample delays
# are fit on onset <= max(sample_date) - this many days, dropping the incompletely-observed
# recent tail (recent onsets that WILL be sampled at long delays are not yet present, biasing
# the raw delay short). Matches TEST_DAYS in 04c_dhis2_delay_windows.R and the train cutoff.
DELAY_TRUNC_BUFFER_DAYS <- 5L

# Single plausibility ceiling (days) for the onset->sample delay: a recorded/observed delay
# longer than this is treated as an implausible data-entry error (an onset->sample gap this
# large is not epidemiologically credible for BDBV) and is excluded UNIFORMLY — from the
# interval-censored MLE fit (04c_dhis2_delay_windows.R, MAX_DELAY), from the empirical
# bootstrap pairs and the imputed-delay clamp, and from the recorded-onset usability test
# (01_data_prep.R) — so the fitted, bootstrapped, and trusted-onset delay supports all agree.
DELAY_MAX_PLAUSIBLE_DAYS <- 60L

# Source of the onset->sample delay for the missing-onset imputation.
#
# PRECEDENCE (01_data_prep.R, .imp_mode). This block previously said the imputation "ALWAYS
# draws from an EMPIRICAL bootstrap ... regardless of this setting". That was inverted by the
# 2026-09-17 rewrite and is no longer true: the SHARED RESOLVER comes FIRST, and the empirical
# bootstrap is now only a fallback that never fires while a fitted delay is on disk. The order
# actually applied is:
#   1. effective_onset_sample_delay()  — the truncation- and censoring-corrected EpiDist
#      marginal fit (or the windowed censored MLE if EpiDist did not run). THE DEFAULT PATH.
#   2. empirical bootstrap of the CURRENT line list's WINDOWED complete pairs (onset in
#      [outbreak floor, max(sample) - DELAY_TRUNC_BUFFER_DAYS]), when >=30 pairs exist and
#      no fitted delay could be resolved. Warns.
#   3. a parametric Exponential at the reported rate, when fewer than 30 pairs exist.
# This constant selects WHICH delay the resolver returns, and hence what step 1 draws from:
#   "data" (default here — this IS the DHIS2 pipeline) — the DHIS2 fit from
#           04c_dhis2_delay_windows.R (dhis2_onset_sample_delay_params.csv), EpiDist marginal
#           preferred over the interval-censored MLE.
#   "lab"  — the fixed lab-linelist Exp(rate) above (DELAY_ONSET_SAMPLE_RATE). Use only when
#           deliberately imputing DHIS2 onsets with the (faster) lab delay. This request is
#           now HONOURED by the imputation itself, not only by the reported rate: before
#           2026-09-19 it fell through to step 2 and drew the DHIS2 empirical delay (~6.8 d)
#           while labelling it "source lab", i.e. it measured the opposite of what it asked for.
ONSET_SAMPLE_DELAY_SOURCE <- "data"

# ---------------------------------------------------------------------------
# Delay estimation: Bayesian EpiDist (truncation- AND double-interval-censoring
# corrected) as the DEFAULT delay estimator.  (04c_dhis2_delay_windows.R)
# ---------------------------------------------------------------------------
# The onset->sample delay drives (a) the missing-onset imputation (~24% of confirmed records:
# ~16% genuine DHIS2 missingness plus the sitrep-reconciliation rows, which have no onset by
# construction — see .PCT_ONSET_IMPUTED / .PCT_ONSET_IMPUTED_LINELIST in 01_data_prep.R) and
# (b) the fast delay-CDF nowcast rate used inside LFO-CV. Two estimators exist: the
# always-on interval-censored MLE (fitdistcens; corrects daily rounding), and the
# Bayesian EpiDist MARGINAL model, which additionally corrects right-truncation (recent
# onsets not yet sampled) and the double-interval censoring of onset & sample dates.
# EpiDist is the gold standard, so it is now the DEFAULT: when TRUE (and the `epidist`
# package is installed) 04c fits it and 01_data_prep.R PREFERS its truncation-corrected
# (family, mean, sd) over the censored-MLE fit. It is a ONE-TIME fit at data-prep time
# (NOT refit per LFO fold or per model), so it does not multiply with the model grid.
# Override with env RUN_EPIDIST=FALSE to fall back to the (faster) censored-MLE delay.
# The 04c scripts additionally gate this on `epidist` package availability, so a box
# without it silently falls back rather than erroring.
RUN_EPIDIST <- {
  .e <- tolower(trimws(Sys.getenv("RUN_EPIDIST", "true")))
  .e %in% c("true", "t", "1", "yes", "y")
}

# ---------------------------------------------------------------------------
# THE EFFECTIVE onset->sample DELAY — single source of truth for every consumer
# ---------------------------------------------------------------------------
# Four places need "the onset->sample delay currently in force":
#   1. 01_data_prep.R    — the parametric FALLBACK draw for the missing-onset imputation.
#   2. 04_nowcasting.R   — the per-LFO-fold right-truncation weights.
#   3. 02_epi_params.R   — EpiNow2's right-truncation model for national R(t).
#   4. reporting/figures — the delay quoted in captions and methods text.
# They used to resolve it independently, and (2)-(4) silently kept the FIXED lab
# Exp(DELAY_ONSET_SAMPLE_RATE) reference (mean 4.39 d) while (1) used the fitted DHIS2
# delay — so the nowcast and R(t) truncation were corrected with a delay ~40% too fast
# and systematically UNDER-corrected the recent weeks. These helpers make all four read
# ONE resolver, so the estimator can never drift apart between consumers again.
#
# Resolution order (ONSET_SAMPLE_DELAY_SOURCE = "data", the default):
#   (a) EpiDist MARGINAL fit  — truncation- AND double-interval-censoring corrected;
#       PREFERRED whenever 04c_dhis2_delay_windows.R wrote epidist_* rows (RUN_EPIDIST).
#   (b) Windowed interval-censored MLE (AIC-best family) — corrects daily rounding only.
#   (c) The fixed lab Exp(DELAY_ONSET_SAMPLE_RATE) reference — last resort.
# With ONSET_SAMPLE_DELAY_SOURCE = "lab", (c) is used unconditionally.
#
# Evaluated LAZILY (a function, never a constant): run_all.R refreshes the params CSV in
# Step 1a, AFTER this file is sourced, so a constant computed here would pin the previous
# run's delay. Results are memoised on the file's (mtime, size) so the hot LFO path does
# not re-read the CSV hundreds of times but a mid-run refresh is still picked up.

DELAY_PARAMS_PATH <- file.path(DATA_DIR, "cfr_reference",
                               "dhis2_onset_sample_delay_params.csv")

.DELAY_CACHE <- new.env(parent = emptyenv())

# NULL-safe "is a usable finite scalar" test (base R's %||% is 4.4+, and this file must
# source on the older R the CI runner may carry).
.fin <- function(x) length(x) == 1L && is.numeric(x) && is.finite(x)

# Load the DHIS2 onset->sample delay params written by 04c_dhis2_delay_windows.R.
# Returns list(family, rate, mean, sd, n_fit, window, params = named numeric of the
# family's native parameters, estimator) or NULL if unavailable.
# `rate` is the Exponential-rate SUMMARY (1/mean); `params` carries the native parameters
# so the draw and the CDF use the ACTUAL best-fit family (e.g. gamma), not an Exponential.
# (Moved here from 01_data_prep.R so 02/04/04b/figures share one definition. Uses
# utils::read.csv so 00_config.R stays dependency-free at source time.)
#' @param stratum optional case-classification stratum ("confirmed", "not_a_case",
#'   "suspected") written by 04c as `<quantity>__<stratum>` rows. NULL resolves the pooled
#'   fit, which is what every consumer did before 2026-09-22.
.load_dhis2_delay_params <- function(path = DELAY_PARAMS_PATH, stratum = NULL) {
  if (!file.exists(path)) return(NULL)
  tab <- tryCatch(utils::read.csv(path, colClasses = "character",
                                  stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(tab) || !all(c("quantity", "value") %in% names(tab))) return(NULL)
  # A requested stratum that is ABSENT must not silently resolve the pooled fit: the pooled
  # fit is majority test-negative on this line list, and one-seventh unadjudicated (6,823
  # not_a_case windowed pairs and 2,038 with no final classification, against 4,981
  # confirmed), so falling back would reinstate the very bias the stratum exists to remove.
  # Warn once, loudly, and then fall back so the run stays non-fatal.
  .has_stratum <- !is.null(stratum) &&
                  any(tab$quantity == paste0("epidist_mean__", stratum))
  if (!is.null(stratum) && !.has_stratum)
    warning(sprintf(paste0("[delay] stratum '%s' was requested but %s carries no ",
                           "epidist_mean__%s row; falling back to the POOLED delay, which ",
                           "mixes every case classification. Re-run 04c (step 1a) with ",
                           "RUN_EPIDIST=TRUE to write the strata."),
                    stratum, basename(path), stratum), call. = FALSE)
  .sfx <- if (.has_stratum) paste0("__", stratum) else ""
  g   <- function(q) { v <- tab$value[tab$quantity == paste0(q, .sfx)]
                       if (length(v)) v[1] else NA_character_ }
  num <- function(q) suppressWarnings(as.numeric(g(q)))
  # PREFER the truncation-corrected EpiDist MARGINAL estimate when 04c wrote it
  # (RUN_EPIDIST=TRUE): it corrects right-truncation AND double-interval censoring, whereas
  # the windowed interval-censored MLE only mitigates truncation by dropping the recent tail.
  # Derive the family's native params from the marginal mean/SD (method of moments for gamma;
  # log-moments for lnorm) so the draw, the CDF and the reported rate all use the corrected
  # distribution. Fall back to the MLE family/rate otherwise.
  epi_fam <- g("epidist_family"); epi_mean <- num("epidist_mean"); epi_sd <- num("epidist_sd")
  if (!is.na(epi_fam) && is.finite(epi_mean) && epi_mean > 0) {
    ep <- numeric(0)
    if (identical(epi_fam, "gamma") && is.finite(epi_sd) && epi_sd > 0)
      ep <- c(shape = (epi_mean / epi_sd)^2, rate = epi_mean / epi_sd^2)
    else if (identical(epi_fam, "lnorm") && is.finite(epi_sd) && epi_sd > 0) {
      .s2 <- log(1 + (epi_sd / epi_mean)^2); ep <- c(meanlog = log(epi_mean) - .s2 / 2, sdlog = sqrt(.s2)) }
    # SELF-CHECK against the file's own `selected_*` rows (written by write_onset_sample_long).
    # Those rows exist precisely so the CSV states which numbers are live; if what this loader
    # derives ever diverges from what the writer declared, the two halves of the contract have
    # drifted and every downstream number is suspect — so say so loudly rather than proceed.
    # The selected_* rows describe the POOLED selection, so they are only a valid check when
    # the pooled fit is what was resolved. Comparing them against a stratum would fire on
    # every stratified read. `.sfx` is "" exactly when no stratum was resolved.
    .sel_mean <- if (nzchar(.sfx)) NA_real_ else suppressWarnings(as.numeric(
                   tab$value[tab$quantity == "selected_mean"][1]))
    .sel_fam  <- if (nzchar(.sfx)) NA_character_ else
                   tab$value[tab$quantity == "selected_family"][1]
    if (is.finite(.sel_mean) && abs(.sel_mean - epi_mean) > 5e-3)
      warning(sprintf(paste0("[delay] %s declares selected_mean = %.3f d but this loader ",
                             "resolves %.3f d; the delay parameter file is inconsistent."),
                      basename(path), .sel_mean, epi_mean), call. = FALSE)
    if (!is.na(.sel_fam) && !identical(.sel_fam, epi_fam))
      warning(sprintf(paste0("[delay] %s declares selected_family = '%s' but this loader ",
                             "resolves '%s'."), basename(path), .sel_fam, epi_fam), call. = FALSE)
    # Posterior interval of the MEAN, when the fit recorded one: consumers that can use an
    # uncertain delay (the EpiNow2 right-truncation model) should not have to treat a fitted
    # nuisance parameter as known exactly.
    .m_lo <- num("epidist_mean_lo"); .m_hi <- num("epidist_mean_hi")
    return(list(family = epi_fam, rate = 1 / epi_mean, mean = epi_mean, sd = epi_sd,
                mean_lo = .m_lo, mean_hi = .m_hi,
                n_fit = num("epidist_n"),
                window = if (nzchar(.sfx))
                           sprintf("epidist_marginal, %s only (truncation-corrected)", stratum)
                         else "epidist_marginal (truncation-corrected)",
                stratum = if (nzchar(.sfx)) stratum else NA_character_,
                params = ep[is.finite(ep)],
                estimator = if (nzchar(.sfx)) paste0("epidist_marginal__", stratum)
                            else "epidist_marginal"))
  }
  fam <- g("family"); rate <- num("rate")
  if (is.na(fam) || !is.finite(rate) || rate <= 0) return(NULL)
  par_q  <- tab$quantity[grepl("^param_", tab$quantity)]
  params <- suppressWarnings(setNames(as.numeric(tab$value[match(par_q, tab$quantity)]),
                                      sub("^param_", "", par_q)))
  params <- params[is.finite(params)]
  out <- list(family = fam, rate = rate, mean = num("implied_mean_fit"),
              n_fit = num("n_fit"), window = g("window"), params = params,
              estimator = "interval_censored_mle")
  out$sd <- .delay_sd(out)   # the MLE rows carry no SD; derive it from the native params
  out
}

# Analytic SD of a delay spec from its native parameters (NA when they are absent, in which
# case consumers fall back to the Exponential summary where SD = mean = 1/rate).
.delay_sd <- function(dp) {
  p  <- dp$params
  ok <- function(nm) length(p) && all(nm %in% names(p)) && all(is.finite(p[nm]))
  switch(as.character(dp$family),
    gamma   = if (ok(c("shape", "rate")))    sqrt(p[["shape"]]) / p[["rate"]]                      else NA_real_,
    weibull = if (ok(c("shape", "scale")))   p[["scale"]] * sqrt(gamma(1 + 2 / p[["shape"]]) -
                                                                 gamma(1 + 1 / p[["shape"]])^2)    else NA_real_,
    lnorm   = if (ok(c("meanlog", "sdlog"))) sqrt(expm1(p[["sdlog"]]^2)) *
                                             exp(p[["meanlog"]] + p[["sdlog"]]^2 / 2)              else NA_real_,
    exp     = if (is.finite(dp$rate) && dp$rate > 0) 1 / dp$rate                                   else NA_real_,
    NA_real_)
}

#' The onset->sample delay actually in force, as a normalised spec.
#'
#' @return list(family, params, rate, mean, sd, estimator, source, n_fit, window).
#'   Never NULL — degrades to the fixed lab Exp reference so every caller always has a
#'   usable delay. `source` is "data" (fitted) or "lab" (fixed reference).
#' @param stratum optional case-classification stratum. The ONSET IMPUTATION passes
#'   "confirmed", because it imputes onsets for confirmed records and the pooled fit is
#'   majority test-negative. Every other consumer leaves it NULL.
effective_onset_sample_delay <- function(path = DELAY_PARAMS_PATH, stratum = NULL) {
  # THE LAB FALLBACK IS A DIFFERENT DISTRIBUTION, NOT A ROUNDING. It is an EXPONENTIAL with
  # mean ~4.39 d where the fitted DHIS2 delay is a GAMMA with mean 7.665 d and shape 0.81 —
  # 43% too fast and the wrong shape. That delay drives the onset imputation, the nowcast
  # multiplier, EpiNow2's truncation model and epinowcast's max_delay, so substituting it
  # silently would move every one of them in the same direction at once. An EXPLICIT request
  # (ONSET_SAMPLE_DELAY_SOURCE='lab') is a documented sensitivity arm and passes quietly; an
  # ACCIDENTAL fallback — params file missing, unreadable, or with no usable row — warns.
  .lab <- function(why, accidental = TRUE) {
    r <- get0("DELAY_ONSET_SAMPLE_RATE", ifnotfound = 0.228)
    if (accidental)
      warning(sprintf(paste0("[delay] falling back to the fixed LAB reference (exp, mean %.2f d) ",
                             "because %s. The fitted DHIS2 delay is a gamma with mean ~7.7 d, so ",
                             "this is ~43%% too fast AND the wrong shape; onset imputation, the ",
                             "nowcast weight, EpiNow2 truncation and epinowcast max_delay are all ",
                             "affected. Run 04c_dhis2_delay_windows.R / run_all.R step 1a."),
                     1 / r, why), call. = FALSE)
    list(family = "exp", params = c(rate = r), rate = r, mean = 1 / r, sd = 1 / r,
         estimator = "lab_reference_fixed", source = "lab", n_fit = NA_real_, window = why)
  }
  if (!identical(get0("ONSET_SAMPLE_DELAY_SOURCE", ifnotfound = "data"), "data"))
    return(.lab("ONSET_SAMPLE_DELAY_SOURCE='lab'", accidental = FALSE))   # deliberate arm

  # Memoise on the file's identity so the hot LFO path does not re-read the CSV, while a
  # mid-run refresh (run_all.R Step 1a rewrites it) still invalidates the entry.
  fi  <- file.info(path)
  # THE STRATUM IS PART OF THE KEY. Without it the first resolved stratum would be served to
  # every later caller on the same file -- the pooled delay to the imputation, or worse the
  # confirmed delay to a consumer that asked for the pooled one.
  key <- paste(path, fi$mtime, fi$size,
               if (is.null(stratum)) "" else as.character(stratum), sep = "|")
  if (!is.null(.DELAY_CACHE$key) && identical(.DELAY_CACHE$key, key))
    return(.DELAY_CACHE$value)

  dp <- .load_dhis2_delay_params(path, stratum = stratum)
  if (is.null(dp) || !.fin(dp$rate) || dp$rate <= 0)
    return(.lab("no fitted delay params on disk"))
  if (!.fin(dp$mean) || dp$mean <= 0) dp$mean <- 1 / dp$rate
  if (!.fin(dp$sd))                   dp$sd   <- .delay_sd(dp)
  if (!.fin(dp$sd))                   dp$sd   <- dp$mean   # Exponential summary fallback
  dp$source <- "data"
  .DELAY_CACHE$key <- key; .DELAY_CACHE$value <- dp
  dp
}

#' CDF of a delay spec, evaluated on the spec's ACTUAL family (not an Exponential
#' moment-match). Falls back to Exp(rate) when the native parameters are absent, so it
#' always agrees with .draw_dhis2_delay()'s fallback for the same spec.
delay_cdf <- function(dp, q) {
  p  <- dp$params
  ok <- function(nm) length(p) && all(nm %in% names(p)) && all(is.finite(p[nm]))
  switch(as.character(dp$family),
    gamma   = if (ok(c("shape", "rate")))    stats::pgamma(q,   shape = p[["shape"]],     rate  = p[["rate"]])  else stats::pexp(q, dp$rate),
    weibull = if (ok(c("shape", "scale")))   stats::pweibull(q, shape = p[["shape"]],     scale = p[["scale"]]) else stats::pexp(q, dp$rate),
    lnorm   = if (ok(c("meanlog", "sdlog"))) stats::plnorm(q,   meanlog = p[["meanlog"]], sdlog = p[["sdlog"]]) else stats::pexp(q, dp$rate),
    stats::pexp(q, dp$rate))
}

#' Quantile function of a delay spec, on the spec's ACTUAL family. Mirrors delay_cdf()
#' family-for-family and fallback-for-fallback, so an upper bound derived from a quantile
#' and the CDF used to correct for it can never come from different distributions.
#'
#' Used to SIZE the modelled delay support rather than hard-coding it: 02_epi_params.R
#' bounds EpiNow2's right-truncation model at the delay's 99th percentile, and
#' 04b_epinowcast.R bounds the reporting triangle the same way. A hard-coded bound goes
#' stale silently the moment the fitted delay changes — which is exactly what happened to
#' the nowcast's 21 days, the 99th percentile of a lab reference that is no longer used
#' and the 92nd of the delay now in force.
delay_quantile <- function(dp, p = 0.99) {
  par <- dp$params
  ok  <- function(nm) length(par) && all(nm %in% names(par)) && all(is.finite(par[nm]))
  q <- switch(as.character(dp$family),
    gamma   = if (ok(c("shape", "rate")))    stats::qgamma(p,   shape = par[["shape"]],     rate  = par[["rate"]])  else stats::qexp(p, dp$rate),
    weibull = if (ok(c("shape", "scale")))   stats::qweibull(p, shape = par[["shape"]],     scale = par[["scale"]]) else stats::qexp(p, dp$rate),
    lnorm   = if (ok(c("meanlog", "sdlog"))) stats::qlnorm(p,   meanlog = par[["meanlog"]], sdlog = par[["sdlog"]]) else stats::qexp(p, dp$rate),
    stats::qexp(p, dp$rate))
  if (length(q) != 1L || !is.finite(q) || q <= 0) NA_real_ else as.numeric(q)
}

#' Draw n delays from a spec. Mirrors delay_cdf() family-for-family and fallback-for-
#' fallback, so a drawn sample and the CDF used to correct for it can never disagree.
.draw_dhis2_delay <- function(dp, n) {
  p  <- dp$params
  ok <- function(nm) length(p) && all(nm %in% names(p)) && all(is.finite(p[nm]))
  switch(as.character(dp$family),
    gamma   = if (ok(c("shape", "rate")))    stats::rgamma(n,   shape = p[["shape"]],     rate  = p[["rate"]])  else stats::rexp(n, dp$rate),
    weibull = if (ok(c("shape", "scale")))   stats::rweibull(n, shape = p[["shape"]],     scale = p[["scale"]]) else stats::rexp(n, dp$rate),
    lnorm   = if (ok(c("meanlog", "sdlog"))) stats::rlnorm(n,   meanlog = p[["meanlog"]], sdlog = p[["sdlog"]]) else stats::rexp(n, dp$rate),
    stats::rexp(n, dp$rate))
}

#' One-line human-readable description of a delay spec, for logs and figure captions.
describe_delay <- function(dp = effective_onset_sample_delay()) {
  sprintf("%s (mean %.2f d, SD %.2f d; %s%s)", dp$family, dp$mean, dp$sd, dp$estimator,
          if (.fin(dp$n_fit)) sprintf(", n=%d", as.integer(dp$n_fit)) else "")
}

# Extra Bayesian PROFILE analyses (refit the featured model ~10x each): a loo-predictive
# posterior over the generation-time mean, and a nowcast-input sensitivity of beta0. Set FALSE
# to skip them (they add several minutes to a run); the core suite/forecasts are unaffected.
BAYES_PROFILE_ANALYSES <- TRUE

# ---------------------------------------------------------------------------
# R(t) AVERAGING WINDOW — ONE constant, shared by EVERY arm
# ---------------------------------------------------------------------------
# How many weeks of the EpiNow2 posterior, ending at the forecast origin, are averaged
# WITHIN each posterior draw to give "the" reproduction number a projection transmits on.
#
# WHY THIS IS ONE CONSTANT. Until 2026-09-19 the two arms silently used different windows:
# the 1-2 week invasion forecast averaged 1 week (bayes_rt_week_draws()'s default) while the
# 13-week cascade averaged 3 (CASCADE_R_WINDOW_WEEKS). On the shipped snapshot those are
# R = 0.83 and R = 1.00 — a 20% difference, and either side of 1, so the two arms of the same
# manuscript disagreed on whether transmission was growing. Worse, 36_report.R asserted they
# shared "ONE national number, from the same model the short-term arm uses". Same model,
# different estimand, different number.
#
# WHY 1 WEEK. A projection needs R AT THE FORECAST ORIGIN. A 3-week mean estimates something
# else — the average over the preceding three weeks — which lags whenever R is trending, and
# the national R has been falling. EpiNow2's posterior is already smoothed by its own GP, so
# a further multi-week mean is a second smoothing; Lison et al. 2024 (PLOS Comput Biol
# 10.1371/journal.pcbi.1012021) warn specifically that smoothing near the present biases R
# toward 1. The cascade's 3-week choice was adopted to match the window of the CONJUGATE
# zone-level estimator, which has since been retired to a diagnostic — so its justification
# no longer exists.
#
# CONSEQUENCE, stated plainly: the 13-week cascade now transmits on R ~= 0.83 rather than
# ~= 1.00. That is a materially less alarming projection, and it is the honest one: it is the
# reproduction number estimated at the origin the projection starts from.
RT_WINDOW_WEEKS <- 1L

# BUMP whenever estimate_rt_epinow2()'s model spec changes (delays/truncation/priors/GT), so the
# on-disk R(t) cache (keyed by gt + analysis_date + this version + GT params) is invalidated
# rather than silently serving a stale rds computed under the old spec. v2 = right-truncation model;
# v3 = summarise at CrIs c(0.6, 0.9) so R(t) reports the 60% (q20/q80) and 90% (q5/q95) bands;
# v4 = right-truncation delay taken from effective_onset_sample_delay() (EpiDist marginal:
#      mean AND sd, so no longer pinned to Exponential shape 1) with max = its 99th percentile,
#      replacing the fixed lab Exp(0.228) / max 21 d spec.
# v5 = the EpiNow2 case series is dated by date_index, exactly as the weekly counts
#      (epinow2_daily_confirmed(), 02_epi_params.R). Onset-less cases were previously dated at
#      their SAMPLE date, which the right-truncation model then inflated a second time, biasing the
#      recent R upward (7 Sep 2026: weekly R 1.50 -> 0.82). Every cached R(t) and Bayesian R draw
#      is invalidated, and no pre-v5 posterior can be carried forward (rt_spec, 21_bayesian_renewal.R).
# 6 (2026-09-17): estimate_rt_epinow2() now carries EpiNow2's `type` column so the 7-day
# forecast tail can be separated from the estimates. The cached column set changed, so a
# v5 cache would be served without it.
# 8 (2026-09-21): the right-truncation distribution is now ESTIMATED from the line-list
#      snapshot archive with EpiNow2::estimate_truncation(), replacing the onset->sample delay.
#      That delay described the wrong leg: the R(t) series is truncated by onset -> APPEARANCE
#      IN THE EXTRACT, which includes an unrecorded data-entry lag. Measured over 54 vintages by
#      matching 5,607 confirmed records to the vintage each first appears in:
#        onset -> sample      median  7.0 d   (what was being modelled)
#        onset -> appearance  median  8.0 d, mean 11.4 d, q90 21 d   (the truth)
#        sample -> appearance median  2.0 d, mean  2.6 d             (the missing leg)
#      The understated delay made EpiNow2 treat the recent weeks as more complete than they are
#      and read the reporting dip as a decline, biasing R(t) DOWN at the present — which is
#      exactly the value RT_WINDOW_WEEKS = 1 hands to the 2-week forecast and the cascade. On
#      the last training week the applied correction was 2.63x where the fitted truncation
#      gives 5.92x (weight 0.1688; pinned in tests/test_truncation.R).
#      Every cached R(t), every rt_draws posterior and the LFO stamp are invalidated.
RT_CACHE_VERSION <- 8L

# ---------------------------------------------------------------------------
# Right-truncation of the R(t) case series (02_epi_params.R section 4b)
# ---------------------------------------------------------------------------
# TWO REGIMES, because two different series are fitted and they are truncated by DIFFERENT
# processes. This is not a stylistic choice; giving both the same distribution mis-corrects one
# of them, and `delta` (fitted on the folds, applied to the deployed forecast) then carries the
# difference into every published probability.
#
#   "extract"  the DEPLOYED fit, on the live line list. A case is present only once it has been
#              entered, so this series is truncated by onset -> APPEARANCE. Estimated from the
#              vintage archive, which is the only place that lag is observable: no column in the
#              line list records it (reporting_date is an ALERT date and is typically EARLIER
#              than the sample; lab_analysis_date is only ~0.8 d after it).
#
#   "asof"     the LFO FOLDS' fit, on history reconstructed by reaggregate_asof(), which
#              censors on linelist_observation_date(). For confirmed cases that date IS the
#              sample date (100% of them), so this series is truncated by onset -> SAMPLE, a
#              shorter process than the extract's.
#              A DIFFERENT ESTIMATOR, because this lag is not hidden: both dates are columns.
#              The completeness curve is MEASURED by onset cohort — for lag L, every record
#              whose onset is at least L days old has a fully determined answer to "observed
#              within L days?", so there is nothing to correct for — and a gamma is fitted to
#              that curve. No Stan, no archive, deterministic, sub-second. See
#              .trunc_fit_asof_cohort() for why snapshot deconvolution and a
#              right-truncation-corrected MLE were both tried and both rejected on measurement.
#              ONE FIT, SHARED BY EVERY FOLD: it is evaluated at ANALYSIS_DATE, never at a
#              fold's own origin, so the truncation cannot become a nuisance parameter
#              correlated with fold index -- which is what run_invasion_lfo()'s SHARED NUISANCE
#              PARAMETER note exists to prevent. The cost is that it sees rows postdating a
#              fold's origin: look-ahead in a nuisance parameter, not in the scored outcome,
#              and one more reason fold skill is optimistic relative to deployment.
#
# NOTE, deliberately not "fixed": the folds consequently see cases the deployed system did not
# yet hold (~2.6 d of them), so LFO skill is optimistic relative to deployment. Correcting that
# means re-censoring the folds on true appearance date, which changes every model's training
# counts and is out of scope here. What matters for the forecast is that each R(t) fit is
# unbiased FOR THE SERIES IT IS GIVEN, so the two remain comparable and `delta` keeps doing its
# own job rather than silently absorbing a truncation mismatch.
TRUNC_ARCHIVE_DIR <- LINELIST_DIR   # the LINELIST_DDMMYYYY vintages
# Longest reporting lag the model represents. The PMF is renormalised over 0..max, so too small
# a max asserts completeness too early. Measured q90 of onset->appearance is 21 d; 45 leaves
# ample tail.
TRUNC_MAX_DAYS <- 45L
# Rolling lookback for vintage selection. Reporting is NOT stationary here — pooled completeness
# at lag 14 rose 0.752 (June) -> 0.827 (July) -> 0.873 (August) — so a long window fits an
# average dominated by the slower past and under-corrects the present.
# EXTRACT REGIME ONLY (the as-of cohort estimator uses no window — windowing would starve the
# long lags, and the short lags that set the anchor week are already dominated by recent cases
# because the epidemic grew).
TRUNC_WINDOW_DAYS <- 56L
# Thin vintages to ~weekly. The archive is near-daily (54 vintages, mean gap 1.7 d); adjacent
# pulls overlap almost completely yet enter the likelihood as independent observations.
TRUNC_THIN_DAYS <- 7L
# Bump whenever the PANEL DEFINITION changes (dating basis, the confirmed filter, the
# censoring rule, the spine). It keys the truncation cache, which is otherwise fingerprinted
# only by the vintage dates and totals — those can be identical across two different panel
# definitions, so a code change alone would silently serve the old fit.
TRUNC_SERIES_VERSION <- 1L
TRUNC_MIN_SNAPSHOTS <- 6L     # EXTRACT only; below this the fit is REFUSED, never approximated
# EXTRACT: the newest frame must carry at least this many cases.
# ASOF: a lag is used in the cohort fit only if its cohort holds at least this many records.
TRUNC_MIN_FRAME_CASES <- 200L
# "fitted" = estimate from snapshots (default). "legacy" = the retired onset->sample delta-method
# Gamma, kept ONLY as a one-line rollback for comparison; it warns loudly every fit.
TRUNC_SOURCE <- "fitted"
# Plausibility gates on the fitted distribution. Realised values today: implied mean 9-12 d;
# completeness at lag 3 is 0.14 (extract) and 0.30 (as-of). The retired spec claimed 0.41.
TRUNC_MEAN_RANGE <- c(4, 25)  # days; implied mean outside this => refuse the fit
TRUNC_Z3_MAX     <- 0.55      # completeness at lag 3 above this => refuse (too optimistic)
# PRIOR for R in the EpiNow2 renewal model (.epinow2_rt_fit(), 02_epi_params.R; EpiNow2's
# rt_opts() applies it to the initial reproduction number). LogNormal on the NATURAL scale, the
# parameterisation of EpiNow2::LogNormal(mean, sd): mean 2.0, SD 0.5 -> meanlog 0.663, sdlog 0.246,
# 90% interval ~[1.29, 2.91], matching the BDBV R0 ~1.5-2.5 range this prior was specified for.
# It is ALSO the Bayesian invasion suite's R when no EpiNow2 posterior can be obtained (no
# observable case, or a fit that fails twice; bayes_rt_week_draws() in 21_bayesian_renewal.R):
# with no information the model's R is its prior. The Bayesian R draws are keyed on it; the R(t)
# diagnostic cache is not, so changing it also requires bumping RT_CACHE_VERSION.
EPINOW2_R_PRIOR <- list(mean = 2.0, sd = 0.5)
# How far back bayes_rt_week_draws() may carry a previous week's EpiNow2 posterior forward when a
# fit fails twice (.rt_recent_posterior(), 21_bayesian_renewal.R). The prior above is a prior on
# the INITIAL reproduction number and is correct in that role; using it as a RECENT week's R would
# centre that week near 2.0 where the data run ~1.0-1.2, so the last successful posterior for the
# SAME generation time is carried forward first and the prior is reached only when no fit has ever
# succeeded. Carrying assumes R has not moved since the borrowed week, which is why it is bounded:
# 8 weeks is roughly two generation intervals at the medium GT (15.3 d), beyond which "R has not
# moved" stops being a defensible claim. Every carried value is warned and written to
# rt_draws_source_log.csv with its origin week and gap, and is never cached — so the next call
# re-attempts the real fit. NOT part of the R-draw cache key (carried values are never cached), so
# changing it does not require an RT_CACHE_VERSION bump.
RT_CARRY_MAX_WEEKS <- 8L

# Onset-to-death delay: Gamma, fit to BDBV-2026 Ituri lab data (n=77)

# BDBV generation-time profiles (discretised Gamma PMFs; built in 02_epi_params.R).
#
# NO Bundibugyo-ebolavirus-specific generation time or serial interval has ever
# been published. Towner et al. 2008 (PLoS Pathog 4(11):e1000212, doi:10.1371/journal.ppat.1000212)
# characterised the virus; MacNeil et al. 2010 and Wamala et al. 2010 (Emerg
# Infect Dis 16(12) & 16(7)) reported clinical features and an incubation period
# of ~5.7-7.4 d for the 2007 Uganda outbreak — but no generation/serial interval.
# We therefore proxy the BDBV generation time with the well-characterised Zaire
# ebolavirus serial interval, standard practice for EVD renewal models (serial
# interval approximate to generation time for filoviruses, given similar latent
# and infectious profiles). All three profiles trace to peer-reviewed sources:
#
#  * Central anchor — serial interval mean 15.3 d (SD 9.3): WHO Ebola Response
#    Team 2014, N Engl J Med 371:1481-1495, doi:10.1056/NEJMoa1411100
#    (https://www.nejm.org/doi/full/10.1056/NEJMoa1411100). Concordant with the
#    pooled random-effects serial interval of 15.4 d [95% CI 13.2-17.5] from the
#    systematic review of Nash et al. 2024, Lancet Infect Dis 24(10):e647-e657,
#    doi:10.1016/S1473-3099(24)00374-8
#    (https://www.thelancet.com/article/S1473-3099(24)00374-8/abstract). No
#    DRC-specific serial interval has been independently estimated: the analysis of
#    the 2014 DRC (Boende) outbreak by Maganga et al. 2014 (N Engl J Med
#    371:2083-2091, doi:10.1056/NEJMoa1411099,
#    https://www.nejm.org/doi/full/10.1056/NEJMoa1411099) itself ASSUMED the
#    West-Africa serial interval of 15.3 +/- 9.3 d, so it stands as a DRC precedent
#    for this anchor rather than an independent estimate.
#  * Short / long profiles are LOW / HIGH sensitivity scenarios bracketing the
#    range of published EVD serial-interval estimates compiled by Van Kerkhove
#    et al. 2015, Sci Data 2:150019, doi:10.1038/sdata.2015.19
#    (https://www.nature.com/articles/sdata201519).
# max_tau is the DAILY PMF SUPPORT — a computational bound, not a scientific parameter.
# It was 35 / 45 / 50 d, which was too short for these right-skewed Gammas: make_gt_pmf()
# truncates there and renormalises over the retained support, and chopping the far tail
# removes a disproportionate share of the VARIANCE. The profiles therefore did not realise
# the literature moments they are labelled with and cited for:
#     short  12.0 / 6.5  -> realised 11.87 / 6.23   (sd -4.1%)
#     medium 15.3 / 9.3  -> realised 14.98 / 8.63   (sd -7.2%)
#     long   18.0 / 10.5 -> realised 17.58 / 9.69   (sd -7.7%)
# i.e. the fitted model used a generation interval materially TIGHTER than WHO 2014's, so
# the FOI kernel was more concentrated on recent weeks than the cited parameters imply.
# 90 d realises all three to within 0.2% on both mean and SD (short +0.0/+0.2%, medium
# +0.0/0.0%, long +0.0/-0.1%) and covers >99.9% of each distribution. Cost: 13 weekly FOI
# lags instead of 5-8 — negligible, because compute_foi()/estimate_R_local() truncate lags
# that reach before the start of the grid anyway (those terms are genuine zeros, correctly
# NOT renormalised), so the extra lags contribute nothing until the series is long enough
# to support them. compute_all_gt_pmfs() reports the realised moments and warns past 2%,
# so any future drift between the label and the realised distribution stays visible.
GT_PROFILES <- list(
  short = list(
    label = "GT-Short (12.0 d)",
    mean  = 12.0,  # low-end sensitivity: just below the pooled 95% CI lower bound (13.2 d, Nash 2024)
    sd    = 6.5,   # CV ~= 0.54
    max_tau = 90
  ),
  medium = list(
    label = "GT-Medium (15.3 d)",
    mean  = 15.3,  # Zaire ebolavirus serial interval, WHO Ebola Response Team 2014 NEJM (proxy; SI ~= GT)
    sd    = 9.3,   # WHO Ebola Response Team 2014 NEJM (SD of the serial interval); CV ~= 0.61
    max_tau = 90
  ),
  long = list(
    label = "GT-Long (18.0 d)",
    mean  = 18.0,  # high-end sensitivity: just above the pooled 95% CI upper bound (17.5 d, Nash 2024)
    sd    = 10.5,  # CV ~= 0.58
    max_tau = 90
  )
)
GT_PRIMARY <- "medium"   # default for the CASCADE / R(t) / SEIR paths (GT_PROFILES sensitivity anchor)

# ---------------------------------------------------------------------------
# Generation-time PRIOR for the Bayesian invasion renewal model  (review §2.1)
# ---------------------------------------------------------------------------
# The reviewer (2026-08-06) flagged that SELECTING the generation time across a
# short/medium/long scenario grid by cross-validation is "akin to fitting" a
# quantity the data cannot identify. Per the response plan we therefore STOP
# selecting the GT and instead place a PRIOR over the GT distribution's
# parameters, anchored on the literature, and MARGINALISE the invasion posterior
# over it (Monte-Carlo grid; see gt_prior_grid() / make_gt_pmf() in 02_epi_params.R
# and the marginalisation in 21_bayesian_renewal.R). GT uncertainty is thus
# represented ONCE, not selected, and the resulting predictive intervals widen to
# reflect our genuine ignorance of the BDBV generation time.
#
# Prior (author decision 2026-08-07 — "span the 12-18 d envelope"):
#   mu_GT ~ Normal(15.3, 1.82)  truncated to [10, 21]  -> 90% approx [12.3, 18.3] d
#           (central anchor: Zaire-ebolavirus SI, WHO Ebola Response Team 2014;
#            width chosen so the 90% mass spans what the old short/long profiles
#            bracketed, 12-18 d, and is concordant with Nash et al. 2024 pooled
#            95% CI 13.2-17.5 d).
#   sd_GT ~ Normal(9.3, 1.50)   truncated to [4, 14]    -> spans the short/long SDs (6.5-10.5).
# n_grid Gauss-Hermite-style grid points (odd, so the anchor is included) are used
# for the Monte-Carlo marginalisation; each grid point's prior weight is the
# (truncated) bivariate-normal density, renormalised over the grid.
# GENERATION TIME IS TREATED AS KNOWN (2026-09-22). GT_MARGINALISE_FEATURED governs the
# 1-2 week featured forecast, CASCADE_GT_MARGINALISE (30_projection_config.R) the 13-week
# projection; both are FALSE, so each analysis uses its single GT anchor and neither
# integrates over GT_PRIOR.
#
# Why: the prior below is an assumption, not an estimate from this outbreak. Its widths
# (mean_sd 1.82, sd_sd 1.50) were chosen, so marginalising over them yields intervals whose
# width reflects that choice while presenting itself as propagated uncertainty. We prefer a
# limitation we can state — the reported intervals contain NO generation-time uncertainty, so
# they are too narrow by an amount we have not quantified — over a widening we cannot defend.
#
# GT_PRIOR is retained because the marginalisation machinery is retained as a sensitivity arm;
# it is also what GT_PROFILES' max_tau is sized against. Setting either flag TRUE runs it.
GT_MARGINALISE_FEATURED <- FALSE

GT_PRIOR <- list(
  mean_mu     = 15.3,  mean_sd = 1.82,   # mu_GT ~ Normal(mean_mu, mean_sd)
  sd_mu       = 9.3,   sd_sd   = 1.50,   # sd_GT ~ Normal(sd_mu,  sd_sd)
  mean_bounds = c(10, 21),               # truncation support for the GT mean (days)
  sd_bounds   = c(4, 14),                # truncation support for the GT sd  (days)
  # Raised 45 -> 90 in step with GT_PROFILES (see the note there): the grid marginalises
  # over GT means up to 21 d and SDs up to 14 d, whose tails are LONGER than the medium
  # anchor's, so a 45 d support truncated the outer grid points hardest — biasing exactly
  # the long-GT arm the marginalisation exists to represent.
  max_tau     = 90L,                     # daily PMF support (days); matches GT_PROFILES
  n_grid_mean = 5L,                      # grid points over the GT mean  (odd -> includes anchor)
  n_grid_sd   = 3L                       # grid points over the GT sd     (odd -> includes anchor)
)
# Alternative GT priors for the SENSITIVITY analysis (shorter- / longer-centred),
# reported per §2.1 step 3 (NOT selectable models — just a couple of re-runs).
GT_PRIOR_ALTS <- list(
  shorter = modifyList(GT_PRIOR, list(mean_mu = 13.0)),
  longer  = modifyList(GT_PRIOR, list(mean_mu = 17.5))
)

# ---------------------------------------------------------------------------
# Kernel-diverse Bayesian ENSEMBLE members  (review §2.4)
# ---------------------------------------------------------------------------
# The ensemble combines a small, pre-specified, STRUCTURALLY-DIVERSE set of
# mobility kernels (the genuine, irreducible uncertainty about how movement routes
# import pressure) as a proper mixture predictive (bayes_ensemble_mixture). GT is
# handled WITHIN each member by GT_PRIOR (marginalised), so it is NOT an ensemble
# axis — the reviewer's point that "some models don't need to be ensembled with a
# prior".
# NOT READ BY ANY CODE PATH. The ensemble members are chosen in run_all.R
# (BAYES_ENSEMBLE_MEMBERS), which takes the fill twin of each member that has one; this
# vector named M15, which is not even built by default. Kept only as a record of the
# original pre-specified set. Do not add a consumer without reconciling the two.

# Fraction of residents' time spent in their home zone, for the manuscript-
# motivated inward / meeting-location force-of-infection variant (Mills 2026).
# A documented modelling assumption (commuting fractions 18-40% in that work =>
# home fractions 0.60-0.82); the central 0.70 is used, sensitivity-analysable.
MOBILITY_HOME_FRACTION <- 0.70

# INCUBATION_MEAN_DAYS / INCUBATION_RATE were REMOVED (2026-09-20): they parameterised the
# E compartment of the stochastic spatial SEIR comparator (08_stochastic_seir.R), which was
# deleted in the 2026-09-17 streamlining. Verified before removal: zero references anywhere
# outside this file. This pipeline has no latent compartment — the invasion hazard is a
# discrete-time renewal/survival model on confirmed onsets, and the generation time
# (GT_PROFILES) is the only transmission-timing parameter it uses.

# ASCERTAINMENT WAS REMOVED FROM THE PIPELINE (2026-09-19). There is no ASCERTAINMENT_NOMINAL
# and no ASCERTAINMENT_GRID, and nothing divides by an assumed reporting fraction.
#
# Every quantity this suite estimates, evaluates and publishes is on the CONFIRMED-CASE scale:
# the forecast target is P(>=1 confirmed case), the models are fitted to confirmed onsets and
# scored against confirmed invasions, and the cascade projects confirmed cases. That scale is
# self-consistent and needs no ascertainment parameter.
#
# What was removed and why:
#   * EpiNow2 obs_opts(scale = Normal(0.45, 0.1)). With rt_opts(pop = Fixed(0)) a CONSTANT
#     scaling of the latent trajectory is exactly confounded with the seeded-infections
#     parameter, so it could not move R(t) by construction. Lison et al. 2024 (PLOS Comput
#     Biol 10.1371/journal.pcbi.1012021) state the same result: constant-in-time ascertainment
#     cannot bias R(t); only TIME-VARYING ascertainment can, and this pipeline neither
#     estimates nor validates one.
#   * The p_infection_* columns and figure series, which were p_case / 0.45 — a uniform
#     rescale by a guessed constant that added no information and invited reading the guess
#     as an estimate.
#   * The cascade within-zone depletion path, whose only implementation converted confirmed
#     cases to infections by dividing by the same constant (and which was permanently off).

# Provinces for which within-province relative-risk scores and maps are produced
# (in addition to the nationwide scores). Strings must match the shapefile PROVINCE
# field exactly. Ituri is the outbreak epicentre; Nord-Kivu and Haut-Uele are the
# highest-exposure neighbours.
PROVINCES_OF_INTEREST <- c("Ituri", "Nord-Kivu", "Haut-Uele")

# CFR_POINT_ESTIMATE was REMOVED (2026-09-20). Its only consumer was the stochastic spatial
# SEIR comparator (08_stochastic_seir.R, deleted 2026-09-17), where it partitioned the
# infectious outflow into death vs recovery. Verified before removal: zero references anywhere
# outside this file. This suite forecasts INVASION (first confirmed case) and projects
# confirmed cases; it models no deaths, so it needs no CFR. The delay-adjusted cCFR remains
# available from the CFR pipeline's own artifact (data/cfr_reference/cfr_summary.csv) for
# anyone who needs it, which is where it belongs — a number this pipeline never reads should
# not be mirrored into its config.

# ---------------------------------------------------------------------------
# OSRM road-network coverage gaps
# ---------------------------------------------------------------------------
# Bokoro (Mai-Ndombe) and Idjwi (Sud-Kivu) have NO finite OSRM entry in either direction --
# Idjwi is an island in Lake Kivu and Bokoro is river-accessed, so there is genuinely no road
# route. OSRM is right; the INFERENCE the kernels drew was not. A zone with an all-NA column
# ends up with an all-zero column in M4/M5/M6 and every composite built on them, i.e. zero
# import hazard forever: it can never be invaded whatever happens around it. On the 2026-09-22
# build those two were the only zero-inflow zones of M4, M5, M6a and the whole M8/M10 family.
# The Flowminder cohort tables record subscriber presence in both, so the data contradict the
# zero directly.
#
# With this ON (the default), load_osrm() imputes the unroutable cells from the great-circle
# separation of the zone centroids through the network's own log-log cost-vs-separation
# relation, so an imputed pair lands where a routable pair at the same separation does. No
# water-crossing penalty is applied -- it would be a free parameter with nothing behind it.
# FALSE restores the previous (zero-inflow) behaviour exactly. See .fill_osrm_gaps().
OSRM_GAP_FILL <- local({
  env <- toupper(trimws(Sys.getenv("OSRM_GAP_FILL", "")))
  if (nzchar(env)) env %in% c("TRUE", "T", "1", "YES", "Y") else TRUE
})

# ---------------------------------------------------------------------------
# Mobility matrix identifiers
# ---------------------------------------------------------------------------
# SOURCE-CELL FILL variants: the same composites, but destinations the empirical
# SOURCE could not observe (its own origin zones; optionally every zone it never
# measured) are taken from the base kernel instead of being asserted as zero.
# Built by 03_mobility_matrices.R under INCLUDE_SOURCEFILL_MODELS.
MOBILITY_FILL_IDS <- c("M8-fill", "M10-fill", "M13-fill", "M13c-fill", "M14-fill",
                       "M16-fill", "M17-fill",
                       "M8-dist-fill", "M10-dist-fill", "M13-dist-fill", "M13c-dist-fill",
                       "M14-dist-fill", "M17-dist-fill")
# ORIGIN-SPLIT cohort composites (split_cohort_rows()), built under INCLUDE_COHORT_SPLIT_MODELS.
# Includes the SHORT-TRIP splits (M8/M10-split): the annex pools its cohort over the three
# epicentre zones, so M1 hands Bunia, Mongbwalu and Rwampara a bit-identical profile — the same
# pooling artefact the cohort split corrects, on the rows that drive the epicentre's spread.
MOBILITY_SPLIT_IDS <- c("M8-split", "M10-split", "M13-split", "M13c-split", "M14-split",
                        "M16-split", "M17-split",
                        "M8-dist-split", "M10-dist-split",
                        "M13-dist-split", "M13c-dist-split", "M14-dist-split", "M17-dist-split")
MOBILITY_IDS <- c("M1", "M2a", "M2b", "M3", "M4", "M4b", "M4c", "M5", "M6a", "M6b",
                  "M7", "M8", "M9", "M10",   # M11 (inward FOI) built on demand in run_all
                  "M13", "M13c", "M14",      # Flowminder-cohort composites (gravity / cohort-gravity / radiation)
                  "M15", "M16", "M17",       # symmetrised Flowminder OD; cohort+relocation OD; all-kernel ensemble
                  MOBILITY_FILL_IDS, MOBILITY_SPLIT_IDS)
# OSRM ROAD-DISTANCE (km) deterrence variants of the mobility kernels, built on demand when the
# osrm road-distance matrix is present (build_all_mobility_matrices' osrm_dist_mat arg; M11-dist
# in run_all). Same construction as their base kernels but keyed on km instead of travel time.
# M13-dist / M14-dist are the geographic-distance analogues of the cohort composites; M17-dist is
# the road-km analogue of the grand all-kernel consensus ensemble (M15/M16 are empirical Flowminder
# flows and are distance-agnostic, so they have no -dist twin).
MOBILITY_DIST_IDS <- c("M4-dist", "M4c-dist", "M8-dist", "M9-dist", "M10-dist", "M11-dist",
                       "M13-dist", "M13c-dist", "M14-dist", "M17-dist")
# The primary kernel: W_primary, which supplies the M11 inward-FOI kernel, the mobility-flow
# figure, and the fallback kernel for the Bayesian over-folds trace.
#
# COHORT-BASED (M14 = cohort + radiation), not the short-trip composite M8. Two reasons:
#   * COVERAGE. M8 fills only the THREE epicentre origin rows from observed mobility (the
#     pooled short-trip annex); every other origin is modelled gravity. M14 fills the
#     Ituri / Nord-Kivu / Tshopo COHORT origin rows from Flowminder cohort subscriber-day
#     presence — a far wider set of origins carrying real movement data, which is what the
#     import force is built out of.
#   * IT IS THE SELECTED KERNEL. The featured model is chosen every run by the leave-future-out
#     CV composite (16_invasion_eval.R), and it selects an M14-fill kernel. Having W_primary
#     be a DIFFERENT kernel from the one the evaluation picks meant the M11 inward-FOI variant
#     and the mobility figure described a kernel the model itself had not chosen.
#
# FILLED form: source-cell fill is the default for every composite (INCLUDE_SOURCEFILL_MODELS,
# MOBILITY_SOURCE_FILL); the unfilled kernel asserts zero flow for destinations the source
# could not observe — the assumption this pipeline rejects everywhere else. Consumers fall
# back to the unfilled parent ("M14") when the filled kernel was not built.
#
# NOTE: this is deliberately NOT the kernel behind the manuscript figures' structural
# baselines. Those are pinned to specific matrices on purpose — M4 for gravity, M1 for the
# Flowminder inflow, OSRM for travel time (see run_all.R) — so they cannot silently follow a
# change here and stop being the three independent comparators they are labelled as.
MOBILITY_PRIMARY <- "M14-fill"

# ---------------------------------------------------------------------------
# Which Flowminder origin-destination export the national kernel (M3) reads
# ---------------------------------------------------------------------------
# WHAT BOTH FILES ARE: Flowminder's national ESTIMATED RELOCATIONS — monthly changes of home
# location, NOT trips (provider definition: "Estimated relocations YYYY_MM-1 to YYYY_MM").
# DEFAULT: the March 2026 file (437 zones; historically called the "provincial PDF extract", in
# fact the HDX column est_flows_2026_03, which it matches on 99.9% of cells by value). Every cell is
# a number there, so a suppressed count and a measured zero are indistinguishable and the gravity
# fit has to treat EVERY zero as left-censored in 0..14.
# ALTERNATIVE: "flowminder__outflow_202604__static.matrix.csv", the April 2026 (est_flows_2026_04)
# national HDX export. It is a strict superset (467 zones; spine zones ABSENT from the matrix fall
# from 82 to 52) and its names are already canonical. Its 57,361 EMPTY cells are 9,416 pairs that
# Flowminder marks "redacted (count <15)" plus 47,945 pairs listed with a BLANK value that month
# (what a blank means is not documented); its 153,499 zeros are, apart from the diagonal, pairs
# absent from the long table; 7,251 cells are positive integers with a minimum of exactly 15.
# build_M4's three-state likelihood treats every EMPTY cell as a redacted 1-14 count, which is an
# assumption for the blanks. Switching changes M3 and therefore EVERY kernel derived from it:
# the M1/M2 fallback rows, the gravity fits M4/M4b/M4-dist (and so the composites M8/M9/M10/M13/
# M14 and their -dist/-fill/-split twins), M7, M15, M16 and M17. Treat it as a sensitivity arm,
# adopt only on LFO-CV evidence, and check the other readers of the March file
# (38_urban_scenarios.R, 40_cascade_next_dominoes.R, 56_covariates_and_burden.R,
# make_manuscript_figures.R) so each follows the switch or is pinned on purpose.
FLOWMINDER_OD_FILE <- local({
  v <- trimws(Sys.getenv("FLOWMINDER_OD_FILE", ""))
  if (nzchar(v)) v else "flowminder__outflow__static.matrix.csv"
})

# ---------------------------------------------------------------------------
# Flowminder cohort mobility composites (M13/M14 + -dist)
# ---------------------------------------------------------------------------
# Four NEW composite kernels fill the Ituri / Nord-Kivu / Tshopo cohort-origin rows from
# Flowminder cohort subscriber-day PRESENCE (data/flowminder_short_trips, "followup" window)
# and take a gravity (M13) or radiation (M14) base elsewhere, keyed on travel time (M13/M14)
# or road-km geographic distance (M13-dist / M14-dist). ADDITIVE — M1/M8/M10 are untouched.
# ON by default; export INCLUDE_COHORT_MODELS=FALSE to skip. Full design + provenance:
# data/flowminder_short_trips/COHORT_INGESTION_PLAN.md.
INCLUDE_COHORT_MODELS <- local({
  env <- toupper(trimws(Sys.getenv("INCLUDE_COHORT_MODELS", "")))
  if (nzchar(env)) env %in% c("TRUE", "T", "1", "YES", "Y") else TRUE
})
# Cohort origin zones (canonical 519-spine names). Each cohort's presence vector becomes the
# outflow row for ALL its origin zones (pooled cohort → identical rows, as with M1).
COHORT_SOURCES <- list(
  ituri  = c("Bunia", "Mongbwalu", "Rwampara", "Nyankunde"),
  nk     = c("Beni", "Butembo", "Katwa"),
  tshopo = c("Lubunga (Tshopo)", "Makiso Kisangani", "Mangobo")
)
# Window driving the cohort composites: "followup" (during-outbreak) is primary; "prior"
# (look-back) is available as a sensitivity by re-sourcing with COHORT_WINDOW = "prior".
COHORT_WINDOW <- "followup"

# ---------------------------------------------------------------------------
# Cohort-CALIBRATED gravity (M4c / M13c + the -dist twins)
# ---------------------------------------------------------------------------
# M4c is the gravity form with its DETERRENCE fitted to Flowminder cohort presence instead of
# to the M3 relocation table, and its destination-mass exponent borrowed from M4. M13c is the
# cohort composite on that base (the M13/M14 analogue). Both ON by default; export
# INCLUDE_COHORT_GRAVITY_MODELS=FALSE to skip, which also drops M4c from the M17 consensus.
#
# WHY. The two Flowminder products see different parts of the distance range. Relocations are
# permanent home changes -- people who move house move far -- so M4's decay is shallow
# (log(d+1) = -1.44). Cohort presence is time spent, which is concentrated near home, so its
# decay is steep (-2.69, within-cohort SE 0.055). Neither is wrong; averaging them in M17 is
# the honest response to not knowing which slice governs importation.
#
# WHY A FITTED KERNEL RATHER THAN A KERNEL MEMBER. The cohort tables fill 10 origin rows of
# 519, so as a consensus MEMBER they would drop out of the other 509 and reproduce the
# coverage-driven weight drift that removed M3 from the consensus. Cohort data generalise as a
# PARAMETER, not as rows. See 03_mobility_matrices.R's M4c block for the full argument, for
# why only the deterrence is fitted, and for the honest (between-cohort) uncertainty.
INCLUDE_COHORT_GRAVITY_MODELS <- local({
  env <- toupper(trimws(Sys.getenv("INCLUDE_COHORT_GRAVITY_MODELS", "")))
  if (nzchar(env)) env %in% c("TRUE", "T", "1", "YES", "Y") else TRUE
})

# ---------------------------------------------------------------------------
# Flowminder symmetrised OD kernel + composites (M15/M16/M17)
# ---------------------------------------------------------------------------
# THREE mobility kernels drawn from the full Flowminder origin-destination data:
#   * M15  Flowminder SYMMETRISED OD kernel, S = O + t(O). NOT inflow-informed: no
#          independent inflow table exists at any stage (the two processed exports are
#          byte-identical, the 437-zone raw pair likewise, and the older 101-zone raw pair
#          is neither identity nor transpose), so the only well-defined object is the
#          symmetrisation of the ONE directed table. Distinct from M3 (directed outflow)
#          because it fills reciprocal edges, hence materially less sparse. OFF by default.
#   * M16  Flowminder COHORT + RELOCATION OD composite. Cohort subscriber-day presence rows
#          (Ituri/NK/Tshopo origins) where available, and the directed Flowminder RELOCATION
#          flows (M3, monthly home-location changes) elsewhere, with M3's coverage gaps filled
#          from radiation (cover_relocation_od()) — the analogue of M13 (cohort + gravity) and
#          M14 (cohort + radiation), pairing the cohort data with empirical flows rather than a
#          model. It pairs with M3, not M15: M15 is not built by default.
#   * M17  ALL-KERNEL CONSENSUS ensemble. Element-wise mean of the distinct structural /
#          empirical BASE hypotheses actually built — {Flowminder relocation OD (M3, coverage
#          gaps filled), gravity (M4), radiation (M5)} — a convex combination, so row-stochastic, with the cohort /
#          epicentre source rows overlaid. The travel-time decay (M6a) and the symmetrised
#          static (M15) are deliberately excluded. M17-dist is the road-km analogue (gravity
#          and radiation re-keyed on km; the empirical M3 carries no distance axis and is reused).
# M15/M16/M17 are ADDITIVE — every existing kernel/model is untouched. The MODEL variants that use
# them are gated by INCLUDE_FLOWSTATIC_MODELS (below); M16 additionally needs the cohort kernel
# (INCLUDE_COHORT_MODELS), and M17-dist needs the OSRM road-distance matrix. "M12" stays reserved
# for the effective-distance kernel on the parallel branch, so these ids start at M15.

# Human-readable mobility-kernel labels for figures/tables (id -> description).
MOBILITY_LABELS <- c(
  M1  = "Short-trip epicentre",
  M2a = "Short-trip epicentre (avg)",
  M2b = "Short-trip epicentre (latest)",
  M3  = "Flowminder relocation OD (national; not trips)",
  M4  = "Gravity fitted to relocations (travel-time)",
  M4b = "Gravity (exp deterrence)",
  M4c = "Gravity fitted to cohort presence (travel-time)",
  M5  = "Radiation (travel-time)",
  M6a = "Travel-time decay (exp)",
  M6b = "Travel-time decay (power)",
  M7  = "IDP-augmented relocation OD",
  M8  = "Short-trip + gravity (travel-time)",
  M9  = "Short-trip + kernel ensemble",
  M10 = "Short-trip + radiation (travel-time)",
  M11 = "Inward/meeting-location FOI",
  M13 = "Flowminder cohort + gravity (travel-time)",
  M13c = "Flowminder cohort + cohort-calibrated gravity (travel-time)",
  M14 = "Flowminder cohort + radiation (travel-time)",
  M15 = "Flowminder symmetrised OD (O+Ot)",
  M16 = "Flowminder cohort + relocation OD",
  M17 = "Consensus of gravity, cohort-gravity and radiation (travel-time)",
  `M4-dist`  = "Gravity fitted to relocations (road-km)",
  `M4c-dist` = "Gravity fitted to cohort presence (road-km)",
  `M8-dist`  = "Short-trip + gravity (road-km)",
  `M9-dist`  = "Short-trip + kernel ensemble (road-km)",
  `M10-dist` = "Short-trip + radiation (road-km)",
  `M11-dist` = "Inward/meeting-location FOI (road-km)",
  `M13-dist` = "Flowminder cohort + gravity (road-km)",
  `M13c-dist` = "Flowminder cohort + cohort-calibrated gravity (road-km)",
  `M14-dist` = "Flowminder cohort + radiation (road-km)",
  `M17-dist` = "Consensus of gravity, cohort-gravity and radiation (road-km)",
  # Source-cell fill variants: same composites, but destinations the SOURCE could
  # not observe are taken from the base kernel instead of being asserted as zero.
  `M8-fill`  = "Short-trip + gravity, source-cell fill",
  `M13-fill` = "Flowminder cohort + gravity, source-cell fill",
  `M13c-fill` = "Flowminder cohort + cohort-calibrated gravity, source-cell fill",
  `M14-fill` = "Flowminder cohort + radiation, source-cell fill",
  `M10-fill` = "Short-trip + radiation, source-cell fill",
  `M16-fill` = "Flowminder cohort + relocation OD, source-cell fill",
  `M17-fill` = "Gravity + cohort-gravity + radiation consensus, source-cell fill",
  `M8-dist-fill`  = "Short-trip + gravity (road-km), source-cell fill",
  `M10-dist-fill` = "Short-trip + radiation (road-km), source-cell fill",
  `M13-dist-fill` = "Flowminder cohort + gravity (road-km), source-cell fill",
  `M13c-dist-fill` = "Flowminder cohort + cohort-calibrated gravity (road-km), source-cell fill",
  `M14-dist-fill` = "Flowminder cohort + radiation (road-km), source-cell fill",
  `M17-dist-fill` = "Gravity + cohort-gravity + radiation consensus (road-km), source-cell fill",
  # Origin-split cohort composites: the pooled cohort profile disaggregated per origin, then filled.
  `M8-split`  = "Short-trip (origin-split) + gravity, source-cell fill",
  `M10-split` = "Short-trip (origin-split) + radiation, source-cell fill",
  `M8-dist-split`  = "Short-trip (origin-split) + gravity (road-km), source-cell fill",
  `M10-dist-split` = "Short-trip (origin-split) + radiation (road-km), source-cell fill",
  `M13-split` = "Flowminder cohort (origin-split) + gravity, source-cell fill",
  `M13c-split` = "Flowminder cohort (origin-split) + cohort-calibrated gravity, source-cell fill",
  `M14-split` = "Flowminder cohort (origin-split) + radiation, source-cell fill",
  `M16-split` = "Flowminder cohort (origin-split) + relocation OD, source-cell fill",
  `M17-split` = "All-kernel consensus ensemble, cohort origin-split, source-cell fill",
  `M13-dist-split` = "Flowminder cohort (origin-split) + gravity (road-km), source-cell fill",
  `M13c-dist-split` = "Flowminder cohort (origin-split) + cohort-calibrated gravity (road-km), source-cell fill",
  `M14-dist-split` = "Flowminder cohort (origin-split) + radiation (road-km), source-cell fill",
  `M17-dist-split` = "All-kernel consensus ensemble (road-km), cohort origin-split, source-cell fill"
)
# Annotate a model/method label (e.g. "Renewal-M13-med", "Bayes-M14-dist") with the readable
# kernel description, matching the LONGEST mobility id first so "M13-dist" wins over "M13"/"M1".
# Returns the input unchanged when no id is present (e.g. "hhh4", "Gravity-B4").
mobility_pretty_label <- function(x) {
  ord <- names(MOBILITY_LABELS)[order(nchar(names(MOBILITY_LABELS)), decreasing = TRUE)]
  vapply(x, function(s) {
    for (id in ord) {
      if (grepl(paste0("(^|[^0-9A-Za-z])", id, "([^0-9A-Za-z]|$)"), s))
        return(sprintf("%s [%s]", s, MOBILITY_LABELS[[id]]))
    }
    s
  }, character(1), USE.NAMES = FALSE)
}

# ---------------------------------------------------------------------------
# READER-FACING MODEL NAMES
# ---------------------------------------------------------------------------
# Every figure that names a model on an axis or in a legend needs a human-readable name, and
# the only way one can be trusted is if it is DERIVED from the label rather than looked up.
# A hand-written lookup (the former MODEL_LABELS in make_publication_figures.R) was missing
# 12 of the 57 cross-validated models, so raw codes like "Bayes-M8-dist-split-geo" were
# printed on the published figure whose entire purpose is to avoid them. Composing the name
# from the SAME tokens bayes_default_grid() composes the label from makes that impossible:
# a new kernel or a new variant is named automatically, and anything genuinely unknown is
# returned unchanged AND reported, rather than silently mislabelled.
#
# A name is built as   <kernel><qualifiers><covariates><beta process><assumptions>
#   Bayes-M14-fill-geo          -> "Cohort + radiation (source-filled) + covariates"
#   Bayes-M14-fill-geo-tvrw1    -> "Cohort + radiation (source-filled) + covariates, beta random walk"
#   Bayes-M17-dist-fill-med     -> "All-kernel consensus (road-km, source-filled)"
#   Bayes-M14-fill-gtshort      -> "Cohort + radiation (source-filled), short generation time"
#   Gravity-B4                  -> "Gravity baseline"
#
# SHORT NAMES, not descriptions. These go on axes where the row pitch is ~8pt: the kernel
# names below are the compact form of MOBILITY_LABELS (which stays the long-form glossary
# used in reports and captions).
MODEL_KERNEL_NAMES <- c(
  M1  = "Short-trip epicentre",       M2a  = "Short-trip epicentre (avg)",
  M2b = "Short-trip epicentre (latest)",
  M3  = "Relocation OD",              M4   = "Gravity",
  M4b = "Gravity (exp deterrence)",   M4c  = "Cohort-calibrated gravity",
  M5  = "Radiation",                  M6a  = "Travel-time decay (exp)",
  M6b = "Travel-time decay (power)",  M7   = "IDP-augmented relocation OD",
  M8  = "Short-trip + gravity",       M9   = "Multi-kernel ensemble",
  M10 = "Short-trip + radiation",     M11  = "Inward (meeting-location) FOI",
  M13 = "Cohort + gravity",           M13c = "Cohort + cohort-gravity",
  M14 = "Cohort + radiation",         M15  = "Symmetrised OD",
  M16 = "Cohort + relocation OD",     M17  = "All-kernel consensus")

# Non-kernel label tokens -> the phrase each contributes, in the order they are appended.
# `med` contributes nothing: it is the grid's generation-time ANCHOR, so naming it would put
# a qualifier on every model that distinguishes none of them.
MODEL_TOKEN_NAMES <- list(
  med     = "",
  geo     = "+ covariates",
  full    = "+ covariates (extended)",
  susp    = "+ suspected-case indicators",
  logit   = "logit link",
  inward  = "inward FOI",
  short   = "short generation time",
  long    = "long generation time",
  gtshort = "short generation time",
  gtlong  = "long generation time",
  tvtrend = "beta log-linear trend",
  tvweek  = "beta weekly random effect",
  tvrw1   = "beta random walk",
  tvar1   = "beta AR(1)",
  tvgp    = "beta Gaussian process")

# Whole-label names for methods that are not kernel+token compositions.
MODEL_FIXED_NAMES <- c(
  `Gravity-B4`        = "Gravity baseline",
  `Distance-B1`       = "Travel-time baseline",
  `Adjacency-B7`      = "Adjacency baseline",
  `Bayes-ens-mean`    = "Bayesian ensemble (mean)",
  `Bayes-ens-median`  = "Bayesian ensemble (median)",
  `Bayes-stacked`     = "Bayesian ensemble (LOO-stacked)")

#' Reader-facing name for a model/method label. Vectorised; NA in, NA out.
#'
#' @param x method labels, e.g. "Bayes-M14-fill-geo-tvrw1".
#' @param warn_unknown warn once, naming the labels that could not be decoded, so a gap shows
#'   up in the run log rather than only in the printed PDF. TRUE everywhere a figure is built.
#' @return character of the same length; an undecodable label is returned UNCHANGED.
model_pretty_label <- function(x, warn_unknown = TRUE) {
  x <- as.character(x)
  out <- vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_character_)
    if (s %in% names(MODEL_FIXED_NAMES)) return(unname(MODEL_FIXED_NAMES[[s]]))
    kern <- mobility_kernel_from_method(s)
    if (is.na(kern)) return(s)
    # Split the kernel id into its family and its qualifiers. The family is the leading
    # M-number (M13c included); -dist / -fill / -split are qualifiers shared by many families.
    fam <- sub("-(dist|fill|split)(-(dist|fill|split))*$", "", kern)
    if (!fam %in% names(MODEL_KERNEL_NAMES)) return(s)
    quals <- character(0)
    if (grepl("(^|-)dist(-|$)",  kern)) quals <- c(quals, "road-km")
    if (grepl("(^|-)split(-|$)", kern)) quals <- c(quals, "origin-split")
    # A -split kernel IS source-filled — the origin-split composites are built on the filled
    # source rows (see MOBILITY_LABELS and bayes_fill_twin_label(), which returns NA for a
    # -split label precisely because it is already filled). Naming only the split would tell
    # the reader that two kernels differ in their fill treatment when they do not.
    if (grepl("(^|-)(fill|split)(-|$)", kern)) quals <- c(quals, "source-filled")
    nm <- unname(MODEL_KERNEL_NAMES[[fam]])
    if (length(quals)) nm <- sprintf("%s (%s)", nm, paste(quals, collapse = ", "))
    # Trailing tokens, in the order they appear on the label, so the name reads in the same
    # order as the code it decodes.
    toks <- strsplit(sub(sprintf("^(Bayes|Renewal)-%s", kern), "", s), "-", fixed = TRUE)[[1]]
    toks <- toks[nzchar(toks)]
    if (!all(toks %in% names(MODEL_TOKEN_NAMES))) return(s)
    # "+ covariates" joins with a space (it continues the noun phrase); everything else is a
    # separate qualifying clause and joins with a comma.
    for (tk in toks) {
      ph <- MODEL_TOKEN_NAMES[[tk]]
      if (!nzchar(ph)) next
      nm <- if (startsWith(ph, "+")) paste(nm, ph) else paste0(nm, ", ", ph)
    }
    nm
  }, character(1), USE.NAMES = FALSE)
  if (isTRUE(warn_unknown)) {
    bad <- unique(x[!is.na(x) & nzchar(x) & out == x & grepl("^(Bayes|Renewal)-", x)])
    if (length(bad))
      warning(sprintf(paste0("[labels] model_pretty_label() could not decode %d model label(s); ",
                             "the raw code will be printed: %s"),
                      length(bad), paste(bad, collapse = ", ")), call. = FALSE)
  }
  out
}

# Mobility-kernel id underlying a model label ("Bayes-M13-dist-fill-geo" -> "M13-dist-fill").
# ONE canonical parser, used by the cascade (.cascade_kernel_from_method), the featured-model
# refits in run_all.R and the report generator: a local regex in any of those silently dropped
# the -fill / -split tokens and refit or projected the WRONG kernel (always wrong once the
# source-cell fill became the default). Strips the trailing model-variant suffixes only.
# Returns NA for non-kernel methods (ensembles, baselines) or an unrecognised id.
mobility_kernel_from_method <- function(method) {
  if (!is.character(method) || length(method) != 1L || is.na(method)) return(NA_character_)
  if (!grepl("^(Bayes|Renewal)-", method)) return(NA_character_)
  # Every non-kernel token a model label may carry, stripped from the END (repeatedly), so
  # what remains is the mobility id. MUST be kept in step with bayes_default_grid(): a token
  # missing here leaves it glued to the kernel, `k %in% known` fails, and the function
  # returns NA — which silently sends the cascade, the over-folds refit and the report to
  # the fallback kernel instead of the featured model's own.
  #   med/short/long  generation-time anchor      gtshort/gtlong  GT sensitivity arms
  #   geo/full/susp   covariate set               logit/inward    link / kernel direction
  #   tvtrend tvweek tvrw1 tvar1 tvgp             time-varying beta_t process
  # The sensitivity-arm tokens come from SENSITIVITY_ARM_SUFFIXES so this parser and the
  # selection gate cannot drift apart; the rest are covariate-set, link and anchor tokens.
  sfx <- sprintf("(med|geo|short|long|full|susp|logit|inward|%s)",
                 paste(get0("SENSITIVITY_ARM_SUFFIXES",
                            ifnotfound = c("gtshort", "gtlong", "tvtrend", "tvweek",
                                           "tvrw1", "tvar1", "tvgp")), collapse = "|"))
  k <- sub("^(Bayes|Renewal)-", "", method)
  k <- sub(sprintf("(-%s)+$", sfx), "", k)
  known <- unique(c(MOBILITY_IDS, MOBILITY_DIST_IDS, "M11", "M11-dist"))
  if (!k %in% known) return(NA_character_)
  k
}

# ---------------------------------------------------------------------------
# Model-suite composition toggles (comprehensive modelling suite)
# ---------------------------------------------------------------------------
# Each flag switches an OPTIONAL family of model variants on or off WITHOUT touching
# the CORE suite, which is always fitted: the base mobility-informed renewal kernels
# (M4/M8/M10) at the medium generation-time anchor (intercept-only import hazard —
# the featured-model class). (Generation time is NO LONGER a grid axis — it is
# marginalised over GT_PRIOR in the featured forecast; see 21_bayesian_renewal.R.)
# The defaults below configure the suite as: gravity M4, the composites M8/M10, the
# Flowminder-cohort composites M13/M14 and the OD/consensus kernels M16/M17, each in its
# SOURCE-CELL-FILL form (the unfilled parents are dropped by INCLUDE_UNFILLED_MODELS =
# FALSE), plus the cohort origin-SPLIT composites, over covariates {none, geo}, on BOTH
# distance measures (OSRM travel time and road km; INCLUDE_OSRM_DIST_MODELS is ON).
# M9, M15, M11, the OSRM-dist family, the FULL covariate set, the suspected-case and logit
# variants and the time-varying-beta families are OFF. The shipped 2026-09-07 run
# cross-validated 52 Bayesian models (verified against attr(lfo_cv_results.rds,
# "lfo_stamp")$bayes_grid, which has length 52; three docstrings said 54 and were wrong).
# Turning OSRM-dist and M11 off on 2026-09-21 takes the default grid to 26. The current count
# is printed by the [config] line at the end of this file. Set a flag TRUE/FALSE here, or export an environment
# variable of the same name, to change it. bayes_default_grid() and run_all.R read these
# via get0()/the flag itself, so a context that sources neither still sees the defaults.
#   INCLUDE_OSRM_DIST_MODELS     road-distance (-dist) deterrence twins — the MASTER switch for
#                                EVERY -dist variant in the Bayesian grid (generic M4/M8/M9/M10-dist
#                                AND cohort M13/M14-dist AND consensus M17-dist); ON by default
#   INCLUDE_GEO_COV_MODELS       the reduced "geo" covariate set (log_pop, CCVI, d_min); ON by default
#   INCLUDE_FULL_COV_MODELS      the FULL-exogenous covariate set (geo + healthsite_density); OFF by
#                                default (the extra covariate roughly halves the covariate sweep)
#   INCLUDE_M11_MODELS           the inward / meeting-location FOI kernel (M11) and its variants;
#                                ON by default since 2026-09-17. M11 encodes a DIFFERENT
#                                transmission assumption from every other kernel, so it can be
#                                selected as the featured model and, through model_selection.json,
#                                set CASCADE_KERNEL — it must not be described as unfitted.
#   INCLUDE_LOGIT_SENS_MODELS    the (slow-to-fit) logit-link observation-process sensitivity; OFF
#   INCLUDE_SUSPECTED_COV_MODELS the suspected-but-not-confirmed leading-indicator covariate
#                                models (own preceding-week + mobility-weighted import); OFF by default
#   INCLUDE_FLOWSTATIC_MODELS    the Flowminder OD family — the M15 symmetrised static kernel, the
#                                M16 cohort+directed-OD composite, and the M17/M17-dist grand
#                                all-kernel consensus ensemble; ON by default
#   INCLUDE_SOURCEFILL_MODELS    the source-cell fill kernels (M8/M10/M13/M14/M16/M17-fill and the
#                                -dist twins of M8/M10/M13/M14/M17), which take the cells the empirical
#                                source could not observe from the base kernel instead of asserting
#                                zero; ON by default, and the DEFAULT form of every composite
#                                (MOBILITY_SOURCE_FILL = "unmeasured")
#   INCLUDE_UNFILLED_MODELS      the UNFILLED originals of those fill twins; OFF by default (a
#                                sensitivity only). When OFF, every Bayesian model on a kernel whose fill
#                                twin is in the grid and built is dropped (note at the flag)
#   INCLUDE_COHORT_SPLIT_MODELS  the origin-split cohort composites (M13/M14/M16/M17-split and the
#                                -dist twins of M13/M14/M17); ON by default (note at the flag)
# Structural-baseline toggles (BASELINE_MODELS in run_all.R — always-on yardsticks). The cheap
# deterministic baselines (Gravity-B4, Distance-B1, Adjacency-B7)
# always run; the two EXPENSIVE per-fold mechanistic ones are gated so they can be skipped for speed:
.model_suite_flag <- function(name, default) {
  # Precedence: an explicit environment variable wins; otherwise the documented default below.
  # We deliberately do NOT fall back to a pre-existing global of the same name: in a persistent
  # R/RStudio session that had sourced an OLDER 00_config.R, get0(name) would return that stale
  # value and silently override a CHANGED default (this is exactly what left INCLUDE_OSRM_DIST_MODELS
  # FALSE — no -dist models — after the default was flipped to TRUE). Honour the default instead,
  # so re-sourcing this file reliably applies the current defaults even in a long-lived session.
  env <- toupper(trimws(Sys.getenv(name, "")))
  if (nzchar(env)) return(env %in% c("TRUE", "T", "1", "YES", "Y"))
  isTRUE(default)
}
# OSRM road-distance (km) deterrence kernels. OFF by default since 2026-09-21.
#
# They are a SENSITIVITY axis (does deterrence measured in km rather than travel-time minutes
# change the ranking?), not a scientific question the manuscript answers, and travel time is
# the better-supported predictor in this outbreak's own data — it is the strongest arrival
# predictor in Figure 1C (Spearman +0.637, ahead of short-trip mobility at +0.602).
#
# They also carry an internal inconsistency: `d_min` in the covariate set is computed from
# `osrm_mat`, i.e. travel-time MINUTES (15_workhorse.R .dmin_vec), on every kernel. So a
# `-dist` model pairs a km-deterrence kernel with a minutes-based frontier covariate.
#
# Measured on the shipped grid: 24 of the 52 specs carry `-dist`, so turning them off leaves
# 28 (x0.538, NOT a halving — M16 has no `-dist` twin because M15/M16 are empirical Flowminder
# flows and are distance-agnostic). The retired comment here claimed "roughly doubles".
# Combined with INCLUDE_M11_MODELS = FALSE the default grid is 26 Bayesian specs, from 52.
#
# THIS FLAG GATES THE GRID AND THE *GENERIC* -dist BUILDS, NOT THE COHORT ONES.
# 03_mobility_matrices.R:2743 gates `build_generic_dist` on this flag (so M4-dist, M5-dist and
# the M8/M10-dist composites stop being built), but `build_cohort_dist` (:2745) is gated only
# on the cohort table being present. M13-dist / M14-dist / M17-dist are therefore STILL built.
# That is deliberate and must stay: CASCADE_SENS_KERNELS (30_projection_config.R:592) names
# "M13-dist-fill", so gating the cohort builds here would silently drop a cascade sensitivity
# arm. The cost is a handful of kernels built but not scored.
INCLUDE_OSRM_DIST_MODELS     <- .model_suite_flag("INCLUDE_OSRM_DIST_MODELS",     FALSE)
INCLUDE_GEO_COV_MODELS       <- .model_suite_flag("INCLUDE_GEO_COV_MODELS",       TRUE)
# THE COVARIATE SETS THEMSELVES, defined once here because they had drifted: the LFO grid
# built them locally in bayes_default_grid() while run_all.R hard-coded a copy for the
# published beta0/covariate trace (compute_bayes_params_over_time). Two definitions of one
# thing is how the trace ends up describing different models from the ones that were scored.
# See bayes_default_grid() for the evidence behind dropping `log_pop` (2026-09-21).
BAYES_GEO_COVARIATES  <- c("ccvi", "d_min")
BAYES_FULL_COVARIATES <- c("ccvi", "d_min", "healthsite_density")
INCLUDE_FULL_COV_MODELS      <- .model_suite_flag("INCLUDE_FULL_COV_MODELS",      FALSE)
# M11 (inward / meeting-location FOI), ON by default since 2026-09-17. It is the only kernel
# family testing a different TRANSMISSION assumption (two-sided, frequency-dependent mixing at
# shared meeting locations) rather than a different distance decay, and until the same date it
# carried a specification error — the per-susceptible contact rate was used as if it were a count
# of introductions, with no receiving-population factor — so it had never been fairly evaluated.
# TURNED OFF 2026-09-21, for a reason unrelated to its specification: M11 and M11-dist are the
# ONLY grid kernels never written to disk (they are built in memory at run_all.R:436/:446 from
# W_primary). .cascade_selected_kernel() (30_projection_config.R:82-84) only adopts a
# CV-selected kernel whose .rds exists, so an M11 win would be SILENTLY discarded and the
# cascade would fall back to CASCADE_KERNEL_FALLBACK — the featured model and the 13-week
# projection would then disagree without saying so. M11 is also a deterministic TRANSFORMATION
# of the primary kernel, not an independent data source, so excluding it costs no data.
# To restore it, set this TRUE *and* persist the kernels at run_all.R:436/:446 first.
INCLUDE_M11_MODELS           <- .model_suite_flag("INCLUDE_M11_MODELS",           FALSE)
# LEGACY KERNELS (M2a, M2b, M4b, M6b, M7): built and saved every run, read by NOTHING.
# Verified 2026-09-21 — every reference to them is their own build/save call plus the
# MOBILITY_IDS registry below; no model spec, figure, baseline or cascade path consumes them.
# M2b is additionally INDISTINGUISHABLE from M1 (max |M1 - M2b| = 3.9e-16): both take the
# latest short-trip snapshot, so they are the same kernel built twice. OFF by default; the
# builders are retained and still unit-tested.
# ROLLING AS-OF PREDICTORS for the invasion design (20_forecast_detail.R).
# The import force is built from nowcast-corrected counts and the generation-time kernel puts
# 31% of its weight on the most recent week -- the week the nowcast inflates most. With one
# count matrix per fold, every TRAINING transition's predictors are complete weeks while the
# FORECAST's predictor carries that inflation, so beta0 is fitted against a lambda about 2x
# smaller than the one it is applied to. Measured pooled observed/expected over the 12 folds:
#     shared matrix (historic)      0.497
#     raw counts everywhere         0.930
#     rolling, no floor             0.440   <- WORSE than changing nothing
#     rolling, floor t >= 3         0.741   <- the realised value, see below
#     rolling, floor t >= 5         0.921
# CORRECTION (2026-09-22): the floor-3 row first read 1.071. That came from a maximum-
# likelihood harness in which four early folds hit complete separation, returned beta0 = 0 and
# predicted ZERO events while still carrying 9 observed ones -- inflating the pooled ratio.
# The pipeline fits with brms, whose prior regularises those folds, and the realised value on
# the 2026-09-22 run is 0.741 (39.2 expected against 29 observed over 11 folds, no fold
# degenerate). Rolling predictors therefore move calibration from 0.497 to 0.741: a large
# improvement, NOT a complete fix. A residual ~1.35x over-prediction remains and delta still
# carries it (median 0.551 across methods; 0.719 for the featured model).
ROLLING_PREDICTORS <- .model_suite_flag("ROLLING_PREDICTORS", TRUE)
# THE FLOOR IS LOAD-BEARING. At small t the as-of reconstruction is sparse, lambda is tiny, and
# the few invasions that did occur imply a huge beta0 (2.14 at the first fold against 0.20 at
# the last); those rows then dominate the fit. 3 keeps 11 of 12 folds; 5 keeps 9. Both sit
# within one Poisson standard error of 1.0 on 29 and 23 events, so 3 is chosen for retaining
# more folds, not because it scores better. NEVER set this below 1.
ROLLING_PREDICTOR_FLOOR <- 3L
INCLUDE_LEGACY_KERNELS       <- .model_suite_flag("INCLUDE_LEGACY_KERNELS",       FALSE)
# WHICH kernel families get an origin-SPLIT twin in the model grid. The split variants
# disaggregate a pooled mobility profile across the origins that share it — the 3 epicentre
# zones for the short-trip annex, the 10 cohort origins for the cohort tables — so they exist
# to repair a real defect: without them those origins carry bit-identical rows.
#
# They do not improve forecasts. Measured on the shipped 12 folds, paired by fold
# (split - fill, log score; negative means split is worse):
#        family     h = 1                  h = 2
#        M8         -0.00002 +/- 0.00022   -0.00103 +/- 0.00130
#        M10        -0.00041 +/- 0.00018   -0.00062 +/- 0.00106
#        M13        +0.00005 +/- 0.00017   -0.00051 +/- 0.00080
#        M14        -0.00043 +/- 0.00018   -0.00091 +/- 0.00076
#        M16        -0.00040 +/- 0.00016   -0.00070 +/- 0.00067
#        M17        -0.00033 +/- 0.00010   -0.00117 +/- 0.00029
# Negative in 11 of 12 cells, and on AUC the split wins only for M8 (+0.004 at h=1, +0.005 at
# h=2) and loses or ties elsewhere. The differences are 2-4% of the log score, and the paired
# design has large power precisely because the twins are near-identical: only 3-10 of 519
# origin rows differ at all. So this is not "split is bad" — it is that per-origin
# disaggregation does not materially change invasion forecasting here, which is a result.
#
# We keep ONE family per pooling mechanism so that result is still supported: M8 for the
# short-trip profile pooled over 3 epicentre zones, M14 for the cohort profile pooled over 10
# origins (and the primary kernel family). The other four only duplicate the same comparison.
# Set to all six to restore the full sweep; set to character(0) to drop splits entirely.
SPLIT_FAMILIES <- c("M8", "M14")
INCLUDE_LOGIT_SENS_MODELS    <- .model_suite_flag("INCLUDE_LOGIT_SENS_MODELS",    FALSE)
INCLUDE_SUSPECTED_COV_MODELS <- .model_suite_flag("INCLUDE_SUSPECTED_COV_MODELS", FALSE)
INCLUDE_FLOWSTATIC_MODELS    <- .model_suite_flag("INCLUDE_FLOWSTATIC_MODELS",    TRUE)
# SOURCE-CELL FILL family, ON by default and the DEFAULT form of every composite:
# M8/M10/M13/M14/M16/M17-fill and M8/M10/M13/M14/M17-dist-fill.
# Both empirical sources exclude their OWN origin zones by construction — the annex
# lists no column for Bunia/Mongbwalu/Rwampara, and the cohort release states that
# "movements within and between those zones are not captured here" — so the classic
# composites assert a hard ZERO on exactly the cells the provider could not measure.
# That removes the dominant local pathway: the weights Beni->Butembo, Makiso
# Kisangani->Mangobo and Makiso Kisangani->Lubunga are all exactly 0 in M13/M14 while
# the radiation base gives 0.051, 0.262 and 0.062, and Butembo, Lubunga and Kilo were
# all invaded. The fill variants take those cells from the base kernel instead. They
# are ADDITIVE: the classic kernels are untouched and LFO-CV decides between them,
# because the refilled mass q reaches 0.71-0.89 on several cohort rows and defaulting
# it would quietly let the base model displace the data.
# CROSS-VALIDATED (LINELIST_07092026, 14 folds): the fill does NOT win in aggregate —
# 1 of 10 pairings on AUC-PR skill, 5 of 10 on rank-of-truth — and the featured model is
# unchanged. Do NOT read that as "the fill does not help": the folds open on 2026-05-12,
# so Kilo (invaded 05-12) and Butembo (05-05) are already affected at every cutoff and
# are scored in NO fold, and on the one target event inside the window (Lubunga, 07-14)
# the fill lifts the rank from 42 to 32 at h=1 and 46 to 32 at h=2 of 472 at-risk zones.
# DECISION (2026-09-17): the fill is the DEFAULT. Asserting zero travel between neighbouring
# outbreak hubs, and to every destination a source did not measure, is not a defensible modelling
# assumption whatever its cross-validated score; the unfilled kernels remain available only as a
# sensitivity (INCLUDE_UNFILLED_MODELS). See METHODS.md 4.0b.
INCLUDE_SOURCEFILL_MODELS    <- .model_suite_flag("INCLUDE_SOURCEFILL_MODELS",    TRUE)
# UNFILLED composites, OFF by default (a sensitivity only). TRUE also fits the unfilled originals of
# the fill twins above — kernels that assert ZERO travel between a source's own origin zones and to
# every destination it did not measure. With FALSE, bayes_default_grid() drops every Bayesian model
# on a kernel whose "-fill" twin is in the grid and built — all covariate, suspected-case, link and
# time-varying variants. Kernels with no fill twin (M4, M4-dist, M9, M9-dist, M11, M15, and the
# origin-split family) are untouched, and nothing is dropped when no twin is available (e.g.
# INCLUDE_SOURCEFILL_MODELS = FALSE). The suspected-case (M8-susp, M8-full-susp) and logit-link
# (M8-*-logit) variants have no fill twin, so they leave the grid. The Bayesian LFO ensemble uses
# the fill twin of each member that has one (Bayes-M8-fill-med, Bayes-M10-fill-med). This changes
# only which models are FITTED: the unfilled matrices are still built (03_mobility_matrices.R does
# not read this flag), and M11 is built from M8.
INCLUDE_UNFILLED_MODELS      <- .model_suite_flag("INCLUDE_UNFILLED_MODELS",      FALSE)
# ORIGIN-SPLIT cohort composites (M13/M13c/M14/M16/M17-split and the -dist twins).
# OFF BY DEFAULT since 2026-09-22 (was ON). Setting this TRUE restores the previous behaviour
# exactly: the matrices are built again and SPLIT_FAMILIES (above) governs which reach the grid.
#
# WHAT THEY DO. A Flowminder cohort is defined by presence in ANY of its origin zones, so the
# release has one profile per cohort and the other cohort kernels copy it to every origin (Beni,
# Butembo and Katwa get identical rows). The split kernels give each origin its own row: the base
# kernel's shape for that origin, rescaled by iterative proportional fitting so the
# population-weighted mix of the cohort's rows reproduces the measured profile exactly
# (split_cohort_rows(), 03_mobility_matrices.R), then filled like the -fill kernels.
#
# WHY THEY ARE OFF. Two reasons, and the second is the one that decided it.
#  * They do not improve forecasts. Paired by fold over 12 folds, split minus fill on the log
#    score is negative in 11 of 12 cells (the table at SPLIT_FAMILIES above).
#  * THE SPLIT DISAGGREGATES USING THE BASE KERNEL. r_o(j) is proportional to B[o,j]*c_j, so
#    each origin's within-support shape is the gravity/radiation shape, merely rescaled to
#    reproduce the pooled marginal. That makes the empirical source rows LESS empirical, which
#    cuts against the reason for having source rows at all.
# NOTE the contrast with the source-cell fill, which is ON despite a similarly weak
# cross-validated case: the fill repairs an assertion of IMPOSSIBILITY (zero travel between
# neighbouring outbreak hubs), whereas identical rows for Beni/Butembo/Katwa is a smoothing
# error that forecloses nothing. The "regardless of CV score" argument applies to one and not
# the other. Measured row TVD against the -fill twin, for the record: Mongbwalu 0.466,
# Beni 0.223, Bunia 0.091, every other origin <= 0.07.
INCLUDE_COHORT_SPLIT_MODELS  <- .model_suite_flag("INCLUDE_COHORT_SPLIT_MODELS",  FALSE)
# Which unobservable cells the fill covers (default "unmeasured"):
#   "unmeasured" — every destination the source did not measure: its own origin zones, the
#                  cohort's "No data" exclusions, and every zone below the annex's ranked cut.
#   "origins"    — only the source's own origin zones (the provider-documented gap); unmeasured
#                  destinations stay at zero.
#   "none"       — disable the fill entirely (equivalent to INCLUDE_SOURCEFILL_MODELS=FALSE).
MOBILITY_SOURCE_FILL <- local({
  v <- tolower(trimws(Sys.getenv("MOBILITY_SOURCE_FILL", "")))
  if (v %in% c("none", "origins", "unmeasured")) return(v)
  if (nzchar(v))
    warning(sprintf("[config] MOBILITY_SOURCE_FILL='%s' is not one of none/origins/unmeasured; using 'unmeasured'.", v),
            call. = FALSE)
  "unmeasured"
})
# Ceiling on the refilled mass q. The filled block is RESCALED by the capped ratio, so
# the row still sums to 1; the cap only stops a base kernel that puts almost all of an
# origin's mass on the unobservable cells from erasing the measured profile.
MOBILITY_FILL_QMAX <- 0.9
# M9 (short-trip + M4/M5/M6a ensemble) and M15 (symmetrised OD kernel O + t(O)) are OFF by
# default: M9 was consistently beaten by the simpler composites and M15 was among the weakest
# (spiky) kernels, so neither is fit unless explicitly re-enabled. Turning these on also builds
# the respective mobility matrices; when off they are neither built nor added to the grid.
INCLUDE_M9_MODELS            <- .model_suite_flag("INCLUDE_M9_MODELS",             FALSE)
INCLUDE_M15_MODELS           <- .model_suite_flag("INCLUDE_M15_MODELS",            FALSE)
# ---------------------------------------------------------------------------
# TIME-VARYING import coefficient beta_t  (ON by default since 2026-09-22)
# ---------------------------------------------------------------------------
# The fixed-beta model asserts one import->invasion conversion rate for the whole epidemic.
# That is the assumption most at risk as the outbreak expands beyond the Ituri epicentre, so
# the suite now carries strict generalisations of the SELECTED model in which beta varies by
# week, and cross-validation prices them against it:
#   tvar1    Ornstein-Uhlenbeck / continuous-time AR(1)    — decays toward beta0; nests both
#            the iid-week process (lengthscale -> 0) and the random walk (-> infinity)
#   tvrw1    first-order random walk                       — PERSISTS the current level
# Three further processes (tvtrend, tvweek, tvgp) were fitted until 2026-09-23 and are
# retired; TV_BETA_TYPES below says why. The code for all five remains, so re-adding one is
# a config change.
# See the TIME-VARYING block in 21_bayesian_renewal.R for the model algebra, the forecast
# behaviour of each and why they are fitted the way they are.
#
# SCOPE. They ride ONE base model — the one the last cross-validation featured — not the
# whole kernel grid. Five extra models rather than five per kernel: the question is whether
# beta has drifted, not whether it has drifted differently under each mobility assumption.
INCLUDE_TV_BETA_MODELS       <- .model_suite_flag("INCLUDE_TV_BETA_MODELS",       TRUE)
# Which processes to fit. Dropping one here removes it from the grid and from every figure
# that reads the grid; nothing else needs changing.
#
# CUT FROM FIVE TO TWO on 2026-09-23. Five arms is a menu, and a menu invites exactly the
# selection this family is supposed to stay out of. Two answer the question.
#   * "trend" is DROPPED as mis-specified, not merely flexible. It rides the ordinary
#     covariate machinery as `week_idx`, which build_invasion_design() standardises by the
#     DESIGN's own mean and SD — in cross-validation, the fold's training window. Those
#     windows run from 4 to 16 weeks here, so gamma_w means "per training-window SD" and
#     carries a different weekly drift in every fold, and the per-week extrapolation step
#     shrinks as the window lengthens. It is not the same model across the folds it is
#     scored on, and it is the one arm that clearly loses: -3.1 skill at one week and -6.8
#     at two against the model it refits, with the worst log score of the six.
#   * "week" and "rw1" are DROPPED as redundant: ar1 nests both exactly (l -> 0 is week,
#     l -> inf is rw1), so ar1 spans the whole carry-over axis with one fitted parameter.
#   * "gp" is DROPPED: over ar1 it adds only smoothness, and with 4 to 16 weekly points the
#     lengthscale is weakly identified, so the prior does most of that work.
# RETAINED: ar1 (the general process, which estimates how much of the current level carries
# forward) and rw1 (its l -> inf endpoint, kept explicitly because a persisting level is the
# maximally adverse case for a fixed-coefficient forecast and therefore the strongest test of
# the assumption at the forecast week).
TV_BETA_TYPES <- c("ar1", "rw1")
# Base model for the time-varying variants. NULL = read the featured Bayesian model from
# outputs/key_outputs/model_selection.json (the pipeline's own record). Set an explicit
# label here to pin them to a specific model instead.
TV_BETA_BASE_MODEL <- NULL

# ---------------------------------------------------------------------------
# Generation-time SENSITIVITY arms
# ---------------------------------------------------------------------------
# The grid is composed at ONE generation time so that GT is not a selection axis. These arms
# refit the SELECTED model at the short and long GT profiles so the assumption's effect on
# (a) cross-validated skill and (b) today's invasion probabilities is measured rather than
# asserted. They are scored but are NOT selection candidates — see INVASION_SELECTION_EXCLUDE.
INCLUDE_GT_SENSITIVITY_MODELS <- .model_suite_flag("INCLUDE_GT_SENSITIVITY_MODELS", TRUE)
GT_SENSITIVITY_PROFILES       <- c("short", "long")
# Methods that are cross-validated but may never be FEATURED. A regex on the method label,
# built from ONE vocabulary (SENSITIVITY_ARM_SUFFIXES) so that adding an arm family to the
# grid cannot leave it silently eligible.
#
# BOTH families belong here, and the time-varying arms did not until 2026-09-23.
#   * The GT arms would make the generation time a fitted axis by the back door, which is
#     exactly what composing the grid at a single anchor prevents.
#   * The time-varying arms are REFITS of a model already in the grid, not competing
#     hypotheses, so selecting one makes the candidate set depend on the previous winner.
#     They were argued to be "genuine competing models"; the evidence says otherwise. Against
#     the model each refits, their one-week skill runs -3.1 to +1.2 against a bootstrap
#     interval of [34.6, 71.7], and the sign of that difference reversed between two
#     consecutive runs. Meanwhile the held-out optimism check, which selects from this same
#     pool, picked a tv arm on both runs and reported an optimism gap of +5.1 (+22.3 on a
#     common base rate) on the second. That gap is what admitting them costs.
#
# Sensitivity arms are still fully cross-validated and reported — being ineligible for
# selection is not the same as being unscored.
SENSITIVITY_ARM_SUFFIXES      <- c("gtshort", "gtlong", "tvtrend", "tvweek", "tvrw1",
                                   "tvar1", "tvgp")
INVASION_SELECTION_EXCLUDE    <- sprintf("-(%s)$", paste(SENSITIVITY_ARM_SUFFIXES,
                                                         collapse = "|"))
message(sprintf(paste0("[config] Optional model families — OSRM-dist:%s geo-cov:%s full-cov:%s M9:%s ",
                       "M11:%s M15:%s legacy-kernels:%s logit:%s susp:%s flowstatic:%s ",
                       "src-fill:%s(%s) cohort-split:%s cohort-gravity:%s osrm-gapfill:%s tv-beta:%s"),
                INCLUDE_OSRM_DIST_MODELS, INCLUDE_GEO_COV_MODELS, INCLUDE_FULL_COV_MODELS,
                INCLUDE_M9_MODELS, INCLUDE_M11_MODELS, INCLUDE_M15_MODELS,
                INCLUDE_LEGACY_KERNELS,
                INCLUDE_LOGIT_SENS_MODELS, INCLUDE_SUSPECTED_COV_MODELS,
                INCLUDE_FLOWSTATIC_MODELS, INCLUDE_SOURCEFILL_MODELS, MOBILITY_SOURCE_FILL,
                INCLUDE_COHORT_SPLIT_MODELS, INCLUDE_COHORT_GRAVITY_MODELS, OSRM_GAP_FILL,
                INCLUDE_TV_BETA_MODELS))

# Short-trip Flowminder snapshots — DISCOVERED DYNAMICALLY from the processed directory so new
# snapshots come online automatically (no hard-coded date list). Each file is named
# flowminder_short_trips__outflow_<YYYYMMDD>__static.matrix.csv; the 8-digit tag is the date.
local({
  fs   <- list.files(FLOWMINDER_ST_DIR,
                     pattern = "^flowminder_short_trips__outflow_[0-9]{8}__static\\.matrix\\.csv$")
  tags <- sub("^.*outflow_([0-9]{8})__static\\.matrix\\.csv$", "\\1", fs)
  tags <- unique(tags[nchar(tags) == 8L])
  dates <- as.Date(tags, format = "%Y%m%d")
  ok <- !is.na(dates); tags <- tags[ok]; dates <- dates[ok]; ord <- order(dates)
  if (length(tags)) {
    FLOWMINDER_ST_TAGS  <<- tags[ord]
    FLOWMINDER_ST_DATES <<- dates[ord]
    message(sprintf("[config] Flowminder short-trip snapshots discovered: %d (%s ... %s)",
                    length(tags), tags[ord][1], tags[ord][length(tags)]))
  } else {
    # Fallback if the directory is unavailable at config time — keeps a working default.
    FLOWMINDER_ST_DATES <<- as.Date(c("2026-04-30","2026-05-07","2026-05-14","2026-05-21","2026-05-24"))
    FLOWMINDER_ST_TAGS  <<- c("20260430","20260507","20260514","20260521","20260524")
    warning("[config] No Flowminder short-trip snapshots found in ", FLOWMINDER_ST_DIR,
            "; using built-in defaults.")
  }
})

# Epicentre zones used as origins in the short-trip cohort (pooled). CANONICAL
# spellings (the 2026-07 shapefile renamed Mongbalu -> Mongbwalu; data/aliases.csv
# still resolves the old spelling, but consumers that intersect these names against
# the zone spine WITHOUT harmonising — naive_epicentre_inflow_scores(), the mobility
# summary's per-origin columns — silently lost the zone.)
EPICENTRE_ZONES <- c("Bunia", "Mongbwalu", "Rwampara")

# ---------------------------------------------------------------------------
# Script guard
# ---------------------------------------------------------------------------
#' TRUE only when the current file is being run AS A SCRIPT (`Rscript foo.R`), FALSE when it
#' is `source()`d from another file or from a session.
#'
#' Several files in this suite are executable SCRIPTS, not function libraries: sourcing them
#' runs a full analysis and OVERWRITES published outputs. run_all.R launches each of them as
#' its own subprocess (`system2(rscript, ...)`) for exactly that reason. Anything that merely
#' wants their FUNCTIONS — a test, an ad-hoc check, a tooling sweep that loads "every module"
#' — must not trigger the run, and `source()` gives no warning that it has.
#'
#' NOTE ON THE IMPLEMENTATION. The familiar `sys.nframe() == 0L` idiom works only when written
#' INLINE at a file's top level (as 04d_delay_fit_figures.R does). Inside a helper function it
#' is always >= 1, because the helper's own frame counts — so a helper must not use it. This
#' asks R directly which file it was launched on.
#'
#' @param file basename of the calling script, e.g. "43_spread_kinematics.R".
is_script_run <- function(file) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) != 1L) return(FALSE)
  identical(basename(f), basename(file))
}

# ---------------------------------------------------------------------------
# Evaluation parameters
# ---------------------------------------------------------------------------
LFO_MIN_TRAINING_WEEKS  <- 4      # minimum weeks before first LFO-CV fold
LFO_HORIZONS            <- c(1L, 2L)    # 1-week ahead and 2-week ahead
N_SIMULATIONS           <- 2000L  # stochastic model simulations per forecast
RANDOM_SEED             <- 20260704L

# --- HOW FAR FORWARD THE CROSS-VALIDATION RUNS -----------------------------------------
# run_invasion_lfo() admits a fold cutoff C at horizon h once the outcome week has aged
#     age(C, h) = ANALYSIS_DATE - (C + 7h + 7)   [days]
# past LFO_MIN_EVAL_AGE_DAYS. The guard exists because invasion truth is the first
# confirmed case bucketed by ONSET, and a zone first infected in the outcome week is only
# counted once its index case has been swabbed and laboratory-confirmed. This run's own
# fitted delays put onset->swab at p90 = 20 d plus swab->lab at p90 = 2 d.
#
# TWO CONSTANTS, TWO DIFFERENT JOBS, and they are deliberately not equal:
#
#   LFO_MIN_EVAL_AGE_DAYS      ADMISSION. How old an outcome week must be to be SCORED at
#                              all. -1 admits every cutoff whose outcome window closes on
#                              or before ANALYSIS_DATE — on the 2026-09-07 frame, cutoffs
#                              through 2026-08-25 at h=1 and 2026-08-18 at h=2, both
#                              ending on 2026-09-07 itself.
#   LFO_RELIABLE_EVAL_AGE_DAYS RELIABILITY. How old an outcome week must be before its
#                              invasion count is treated as settled. Folds younger than
#                              this are scored but FLAGGED provisional, and the pipeline
#                              additionally publishes the restricted-fold metrics beside
#                              the full ones (invasion_evaluation_reliable_folds.csv).
#
# WHY ADMISSION WAS LOOSENED FROM 25 TO -1 (2026-09-22). At 25 d the newest scored round
# was 2026-07-28 while the deployed forecast was issued at 2026-09-07: six weeks of
# observed epidemic — the six weeks in which the outbreak spread furthest beyond the Ituri
# epicentre — were never cross-validated, and every over-rounds panel (rank evolution,
# front approach, beta over folds, skill over time) carried a five-week hole between the
# last CV round and the live run.
#
# WHAT IT COSTS, STATED PLAINLY. The newest folds are RIGHT-TRUNCATED: a zone invaded late
# in an outcome window that has had only days to report can be scored is_new_invasion = 0
# when it was in fact invaded. That biases the newest folds' apparent skill DOWNWARD
# (missed positives) and biases calibration-in-the-large UPWARD (predicted invasions that
# are real but not yet observed read as over-prediction). It cannot bias the RANKING
# between models: every model is scored on the identical rows with the identical truth.
# The reliable-fold table is the sensitivity that quantifies the size of the effect.
# SET TO 13 (2026-09-23). -1 admitted outcome windows that closed on the analysis date
# itself, and those are materially incomplete: measured on the featured model, folds younger
# than 25 days carry observed/expected 0.37 (h=1) and 0.42 (h=2) against 0.74 and 0.69 for
# the settled folds -- i.e. the newest rounds look over-predicted mostly because their
# invasions have not been reported yet. 13 days covers the bulk of the onset->confirmation
# delay (p90 onset->swab 20 d is longer, which is why the RELIABLE flag below stays at 25 and
# the restricted-fold sensitivity is still published).
#
# WHY 13 AND NOT 14. At 14 the h=2 round with cutoff 2026-08-04 (forecast origin 10 Aug) was
# excluded by ONE DAY: its window closes 24 Aug, so its age on the 2026-09-07 frame is 13.
# 13 admits it and both horizons then run to the 10 Aug origin. The extra day carries no
# epidemiological content; the boundary simply landed on a round.
LFO_MIN_EVAL_AGE_DAYS      <- 13L
LFO_RELIABLE_EVAL_AGE_DAYS <- 25L

# COMMON FOLD GRID (16_invasion_eval.R, and mirrored by plan_folds() in write_run_info.R).
#   TRUE  (default) — the FARTHEST horizon must be complete for a cutoff to be admitted, so
#                     every admitted round is scored at EVERY horizon and h=1/h=2 are pooled
#                     over the IDENTICAL set of origins.
#   FALSE           — the "Q1 fold augmentation": the SHORTEST horizon governs admission and
#                     each fold scores only the horizons whose own window is complete. That
#                     adds the recent cutoffs where h=1 is evaluable and h=2 is not (leak-free
#                     -- an incomplete horizon is simply skipped), at the cost of the two
#                     horizons being pooled over DIFFERENT origins.
#
# WHY NO SINGLE AGE THRESHOLD CAN ALIGN THE GRIDS. The h=1 window at cutoff C and the h=2
# window at cutoff C-7 both CLOSE on the same day, C + 13, so they always carry the IDENTICAL
# age -- no value of LFO_MIN_EVAL_AGE_DAYS admits one and excludes the other. On the
# 2026-09-07 frame at 13 days the per-horizon rule gives h=1 fourteen rounds (to origin
# 17 Aug) against h=2's thirteen (to origin 10 Aug); the common grid gives both thirteen,
# to 10 Aug.
#
# WHAT IT COSTS, AND WHAT IT BUYS. The near-term horizon gives up its most recent round --
# the augmentation's whole purpose. Bought back: every pooled metric, the recalibration delta
# and every over-rounds figure describe the same origins at both horizons, so an h=1-vs-h=2
# difference cannot be an artefact of one horizon carrying an extra, most-truncated round.
# On THIS frame it costs h=1 nothing: at 13 days the common grid reproduces exactly the 13
# rounds h=1 already had at 14, so only h=2 changes (12 rounds -> 13).
LFO_COMMON_FOLD_GRID <- TRUE

# --- FOLD CUTOFF vs FORECAST ORIGIN ------------------------------------------------------
# `cutoff` is the START of the last TRAINING week, not the date the forecast was made. The
# forecast is issued once that week has CLOSED, so the real-time origin is cutoff + 6 -- the
# value run_invasion_lfo() actually passes as the as-of date (16_invasion_eval.R, the
# "FOLD ORIGIN = cut + 6" block) and that 33b_cascade_calibration.R passes as issue_date.
#
# EVERY DATED ROUND LABEL MUST USE THIS. Labelling a round by its cutoff prints it six days
# before the forecast it names: on the 2026-09-07 frame the last 2-week round reads "28 Jul"
# when the forecast was issued on 3 Aug, and Figure 3 then draws a 34-day gap to the live
# run where the true gap is 28 days -- because the live point is anchored on ANALYSIS_DATE,
# a week END, while the CV points were anchored on week STARTS.
#
# DISPLAY ONLY. Window arithmetic stays on `cutoff`: the horizon targets are defined against
# it (first_wk in (cutoff, cutoff + 7h]), as is the admission age (cutoff + 7h + 7).
LFO_ORIGIN_OFFSET_DAYS <- 6L
lfo_origin <- function(cutoff) as.Date(cutoff) + LFO_ORIGIN_OFFSET_DAYS

# ---------------------------------------------------------------------------
# Post-hoc recalibration of invasion probabilities (16b_invasion_recalibration.R)
# ---------------------------------------------------------------------------
# The invasion suite discriminates well but over-predicts in the large (every Bayesian
# model has calibration_in_large > 1; most have NEGATIVE Brier skill). A single
# hazard-scale multiplier delta per (method, horizon), p -> 1 - (1-p)^delta, fixes the
# level without touching ANY ranking product (the transform is strictly monotone).
#
# Three switches, deliberately separate, because they carry very different risk:
#   INVASION_RECALIBRATE        purely ADDITIVE — fits delta prequentially and reports
#                               recalibrated metrics alongside the raw ones. Changes no
#                               existing number. Safe to leave on. DEFAULT: TRUE.
#   INVASION_SELECT_ON_RECAL    lets the selection composite's log-score axis use the
#                               recalibrated score. This CAN change the featured model —
#                               and the featured model sets CASCADE_KERNEL. DEFAULT: TRUE.
#   INVASION_RECALIBRATE_DEPLOY applies the pooled delta to the LIVE forecast, changing
#                               the probabilities that reach the watch-list and maps.
#                               Rankings are unchanged; absolute probabilities are not.
#                               DEFAULT: TRUE.
#
# WHY BOTH ARE ON (changed 2026-09-11, previously both FALSE). The two flags answer the
# SAME question and must not disagree: if the level error is corrected before the forecast
# is published, selection must stop penalising models for it; if it is not corrected,
# selection must keep doing so. Splitting them is the only incoherent setting, and it was
# the setting we had.
#   * Raw log score on a ~1% base rate is dominated by calibration-in-the-large, a
#     ONE-PARAMETER defect this module removes. Leaving it on the selection axis makes the
#     composite's third axis largely a second reading of a level error we then correct.
#     The recalibrated axis measures refinement instead. The other two axes are rank
#     metrics and are invariant to the transform, so this is the only channel through
#     which recalibration can move the pick (see .invasion_ls_col, 16_invasion_eval.R).
#   * The switch was verified to be a NO-OP on the 2026-09-07 evaluation before being
#     flipped: Bayes-M14-med wins under BOTH axes by 5 composite points, over Bayesian
#     singles, singles + ensembles, and the full field; the runner-up (Bayes-M14-geo)
#     shares the M14 kernel, so CASCADE_KERNEL is M14 either way. Flipping a rule while
#     it changes nothing is deliberate — the alternative is discovering the flip in a run
#     where it also moves the 13-week projections. write_model_selection.json now records
#     the pick under BOTH axes every run so a future divergence is visible immediately.
#   * Deployment: the featured model over-predicts about two-fold (calibration_in_large
#     1.92 at h=1, 2.00 at h=2); delta = 0.464 [0.34, 0.63] at h=1, no band boundary hit,
#     and the pooled value matches the final prequential value to three decimals, i.e. the
#     level error is persistent and well identified rather than an artifact of one fold.
#     Recalibrated, 44 of 56 Bayesian slices have POSITIVE Brier skill (raw: 24 of 56).
# TWO THINGS THIS DOES NOT DO, stated so they are not assumed:
#   * It does not make the featured model calibrated — a single multiplier leaves
#     calibration_in_large near 1.3, i.e. roughly 30% residual over-prediction.
#   * It does not touch the 13-week cascade, which fits its OWN factor in
#     cascade_fit_delta() (30_projection_config.R) from the RAW calibration column. The
#     insulation is deliberate (see the note at the Bayesian suite in run_all.R); the
#     correction is applied once per layer and never compounded.
# delta is phase-dependent (the first LFO fold wanted ~1.4 during explosive growth and the
# pooled value is now ~0.4), so a resurgence would make a pooled factor UNDER-predict.
# Monitor the gap between `delta` and `delta_preq_final` in invasion_recalibration.csv.
.recal_flag <- function(nm, default) {
  env <- Sys.getenv(nm, unset = NA_character_)
  if (!is.na(env) && nzchar(env)) {
    v <- tolower(env)
    if (v %in% c("1", "true", "yes", "on"))  return(TRUE)
    if (v %in% c("0", "false", "no", "off")) return(FALSE)
    warning(sprintf("[config] unreadable %s='%s'; using default %s.", nm, env, default),
            call. = FALSE)
  }
  default
}
INVASION_RECALIBRATE        <- .recal_flag("INVASION_RECALIBRATE",        TRUE)
INVASION_SELECT_ON_RECAL    <- .recal_flag("INVASION_SELECT_ON_RECAL",    TRUE)
INVASION_RECALIBRATE_DEPLOY <- .recal_flag("INVASION_RECALIBRATE_DEPLOY", TRUE)
# RECALIBRATION HAS NO TUNING CONSTANTS (2026-09-22). delta is fitted by JEFFREYS-PENALISED
# maximum likelihood (16b_invasion_recalibration.R), which is finite for every training pool
# -- including one with no invasions, where the plain MLE collapses to delta -> 0 and would
# zero out a whole fold's forecast. Because the estimate always exists, the two floors this
# block used to carry are gone:
#   RECAL_MIN_EVENTS <- 10L   # "enough events to trust delta"
#   RECAL_MIN_FOLDS  <- 3L    # "enough origins to trust delta"
# Both were chosen, and neither was doing work. Measured prequentially on the shared support
# over the 20 Bayesian models, sweeping min_folds in {1,2,3} x min_events in {0,1,5,10}, the
# mean LOG SCORE moves by 0.9% at h=1 (0.02335-0.02355) and 0.6% at h=2 (0.04243-0.04267),
# against a 3.4%/3.3% gain from recalibrating at all (raw 0.02416 / 0.04399). Per-fold
# |log(obs/exp)| cannot separate them either: at ~3.3 events per fold its value under PERFECT
# calibration is 0.451 (90% range 0.279-0.635), and every configuration lands inside that.
# The data cannot distinguish the settings, so we take the one with no free constants rather
# than present a tuned choice as if it were justified.
#
# What remains is an INFORMATION constraint, not a threshold: delta is fitted as soon as at
# least one training fold's outcome window has closed, and before that there is nothing to
# fit and the forecast is left uncorrected (delta = 1).
RECAL_MIN_FOLDS  <- 1L
# Numerical search bracket for delta, NOT a prior. The penalised objective is bounded and
# unimodal (verified on all 338 real training pools), so the optimum is interior and this
# never binds; hitting it means something is badly wrong and is reported as such. It spans
# both sides of 1 because a model can genuinely UNDER-predict (this outbreak's first LFO fold
# did, at delta ~1.4 during explosive growth).
RECAL_BAND       <- c(1e-4, 1e4)
RECAL_N_BOOT     <- 400L   # zone-clustered bootstrap replicates for the delta interval

# Parallelism for independent Bayesian fits. Fans out the per-fold LFO models,
# current-forecast suite, GT profile, nowcast sensitivity, and parameter-over-time
# refits across a `future` multicore pool (LFO fold 1 stays sequential to warm the
# cmdstanr compile cache). Each fit is independently seeded, so results are unchanged
# — only wall-clock differs. Each fit uses 2 chains/cores, so N jobs ≈ 2N cores.
# DEFAULT is now auto-parallel: ~half the physical cores, clamped to [1, 6], leaving
# headroom for the OS and each fit's 2 chains. Override with the PARALLEL_JOBS
# environment variable (e.g. PARALLEL_JOBS=1 forces the sequential historical path;
# raise it on a big box with ample RAM). NOTE: `future::multicore` forks, so it is a
# no-op (sequential) under RStudio or on Windows — run via `Rscript` in a terminal to
# get the fan-out on macOS/Linux.
PARALLEL_JOBS <- local({
  env <- suppressWarnings(as.integer(Sys.getenv("PARALLEL_JOBS", "")))
  if (!is.na(env) && env >= 1L) return(env)
  # Honour a value preset as an R variable before this file is sourced.
  v <- get0("PARALLEL_JOBS", ifnotfound = NA_integer_)
  if (is.numeric(v) && length(v) == 1L && !is.na(v) && v >= 1L) return(as.integer(v))
  nc <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
  if (is.na(nc) || nc < 1L) nc <- 2L
  max(1L, min(6L, nc %/% 2L))
})

# WIS quantile levels (Bracher et al. 2021, PLOS Comp Biol). The two central-interval
# levels the pipeline reports and scores: alpha=0.10 -> q5/q95 (90% CI) and alpha=0.40 ->
# q20/q80 (60% CI). These are exactly the quantiles the forecasters emit (q05/q20/q80/q95),
# and they honour the pipeline band convention (90%/60%; never 95%/50%). The WIS also uses
# the median (weight w_0 = 1/2) per Bracher et al.
WIS_ALPHA_LEVELS <- c(0.10, 0.40)
WIS_QUANTILE_LOWER <- WIS_ALPHA_LEVELS / 2       # c(0.05, 0.20)  -> q5,  q20
WIS_QUANTILE_UPPER <- 1 - WIS_ALPHA_LEVELS / 2   # c(0.95, 0.80)  -> q95, q80

# Top-K for precision / accuracy metrics
TOPK_VALUES <- c(3L, 5L, 10L)

# ---------------------------------------------------------------------------
# Visualisation defaults
# ---------------------------------------------------------------------------
VIZ_THEME_BASE <- 12          # ggplot base font size
VIZ_WIDTH_WIDE  <- 14         # inches
VIZ_WIDTH_NARROW <- 8
VIZ_HEIGHT_STD  <- 7
VIRIDIS_OPTION  <- "plasma"   # default colour scale for risk maps

# ---------------------------------------------------------------------------
# RETAINED-FIGURE ALLOW-LIST (2026-09-17 streamlining)
# ---------------------------------------------------------------------------
# The suite historically wrote ~250 PDF/PNG figures per run, of which a small
# curated set is actually published. FIGURE_KEEP is the definitive list of
# figure BASENAMES (no directory, no extension) that a run is allowed to write.
#
# HOW IT IS ENFORCED: every save helper in the suite (.save, .fd_save,
# .fd_save_png, .rv_save, .cascade_ggsave, and each save_dual) calls
# figure_is_kept() and returns without writing when the answer is FALSE. The
# allow-list is a BACKSTOP, not the primary mechanism — the plot CALL SITES for
# dropped figures are removed too, so the work is not done and then discarded.
# The backstop exists because several retained plotters are shared library
# functions called in loops over horizons, provinces and models, where the same
# function legitimately produces both kept and dropped files.
#
# HORIZON SUFFIXES: an entry without a trailing _h<N> matches the bare name AND
# any _h1/_h2/_h4/_h8/_h13 variant, so "bayes_priority_scatter" keeps both
# horizons. Give an explicit _h<N> suffix to keep one horizon only.
#
# TO PUBLISH A NEW FIGURE: add its basename here AND make sure its call site is
# live. Setting FIGURE_KEEP_ALL=1 in the environment disables the gate entirely
# (restores the historical behaviour for debugging).
FIGURE_KEEP <- c(
  # --- delay fits (04d) -----------------------------------------------------
  "dhis2_delay_epidist_fits",
  # The onset->sample delay is fitted per case classification (04c) because the onset
  # imputation imputes CONFIRMED onsets and the pooled fit is majority test-negative. This
  # panel is the only place that contrast is visible; without it on the allow-list, .dfg_save()
  # would build it and then silently drop it.
  "dhis2_onset_sample_by_class",
  # --- Bayesian diagnostics (17 / 20, via run_all.R) ------------------------
  # Is the calibration offset constant over time? The evidence behind pooling ALL folds for
  # the deployment delta rather than using a recent window, and behind NOT fitting a
  # time-varying beta_0 -- the per-fold delta is exactly the quantity such a model would fit.
  "delta_stability",
  # Figure 2B shows three of sixteen rounds so the main figure stays legible; this shows
  # EVERY round, including the early ones that carry no recalibration factor yet.
  "FigureS5_topk_all_rounds",
  "bayes_beta_over_folds", "bayes_params_over_time", "rt_national",
  "skill_over_time_auc_pr_skill",
  "bayes_discrimination_summary", "bayes_lfo_forecast_vs_outcome",
  "bayes_invasion_rank_map_national", "bayes_priority_scatter",
  "bayes_prob_vuln_choropleth_national",
  # --- review panels (24_review_figures.R) ----------------------------------
  "FigR1_reliability", "FigR3_prospective_ranks", "FigR4_heldout_optimism",
  # --- manuscript figures (make_manuscript_figures.R) -----------------------
  "Figure1", "Figure2", "Figure3", "Figure4",
  "FigureS1_arrival_vs_predictors", "FigureS2_operational_prob_priority",
  "FigureS3_epicentre_case_share", "FigureS4_front_rank_evolution",
  "Figure_spread_kinematics_compact",
  # --- supplementary model figures (make_si_model_figures.R) ----------------
  "FigureS5_roc_pr_curves", "FigureS6_mcmc_diagnostics",
  "FigureS7_calibration_over_time", "FigureS8_generation_time_sensitivity",
  "FigureS9_time_varying_beta",
  # --- supplementary data / mobility figures (make_si_data_figures.R) -------
  "FigureS10_invasion_by_source", "FigureS11_linelist_vs_sitrep",
  "FigureS12_mobility_kernels",
  # --- publication figures (make_publication_figures.R / make_topk15_ever.R)--
  "Figure2_labelled", "F_topk15_ever",
  # --- cascade (35 / 37 / 39 / 42 / make_manuscript_figure2_cascade.R) ------
  "cascade_corridors", "Figure_kindiv_sweep",
  "Figure2_cascade", "Figure4_cascade",
  # Intermediate: 37_cascade_figure4.R builds the h13 panel and copies it to the canonical
  # Figure4_cascade name. h4/h8 are excluded via FIGURE_DROP below.
  "Figure4_cascade_h13",
  # --- urban scenarios (38) — Kinshasa only ---------------------------------
  "Figure_urban_Kinshasa"
)

# PATH-QUALIFIED EXCLUSIONS, checked before the allow-list. The allow-list matches on the
# figure STEM, which cannot by itself distinguish two different figures that share a name in
# different directories — and three do:
#   * Figure2          — RETAINED as manuscript_figures/Figure2; the figures/ one is a
#                        second, differently-built figure that the brief does not ask for.
#   * Figure4_cascade  — RETAINED as figures/Figure4_cascade. 37_cascade_figure4.R builds it
#                        per horizon and copies h13 to that canonical name, so h13 is kept as
#                        a necessary INTERMEDIATE while h4/h8 are not wanted.
# Each entry is "<parent directory>/<stem>".
#   * Figure2 — RETAINED only as manuscript_figures/Figure2; figures/Figure2 is a second,
#                        differently-built figure that the brief does not ask for.
#
# NOTE on Figure1 / Figure3: BOTH directories are retained deliverables (the brief asks for
# manuscript_figures/Figure1 and figures/Figure1, and likewise Figure3), even though the two are
# genuinely different figures — figures/Figure3 has panels A-E, manuscript_figures/Figure3 has
# A-B. They are therefore NOT dropped. The real hazard is CITATION, not publication: a reference
# to "Figure 3C" resolves only in the publication variant. update_bayesian_report.R names the
# directory explicitly for exactly that reason.
FIGURE_DROP <- c("figures/Figure2",
                 "figures/Figure4_cascade_h4",
                 "figures/Figure4_cascade_h8")

# TRUE when `path` (a filename, basename or full path) is on the allow-list.
# Non-figure artifacts (csv/json/rds/md/log) are ALWAYS allowed: the gate governs
# figures only, per the streamlining brief.
figure_is_kept <- function(path) {
  if (identical(tolower(trimws(Sys.getenv("FIGURE_KEEP_ALL", ""))), "1")) return(TRUE)
  b <- basename(as.character(path)[1])
  ext <- tolower(tools::file_ext(b))
  # Three cases. A DATA extension is never gated. A FIGURE extension is gated on
  # its stem. NO extension means the caller passed a bare stem (every save_dual()
  # in the suite does: it appends .pdf/.png itself), which is gated as a stem —
  # treating that as "not a figure" would silently disable the gate entirely.
  if (ext %in% c("csv", "json", "rds", "rda", "md", "log", "txt", "tsv", "xlsx", "html"))
    return(TRUE)
  stem <- if (nzchar(ext)) sub("\\.[^.]+$", "", b) else b
  keep <- get0("FIGURE_KEEP", ifnotfound = character(0))
  # A raw/ sibling is a SCALE TWIN of a primary figure, not a deliverable in its own right
  # (RUN_RAW_FIGURE_TWIN, run_all.R). Never write one through the gate.
  parts <- strsplit(as.character(path)[1], .Platform$file.sep, fixed = TRUE)[[1]]
  if (length(parts) > 1L && "raw" %in% parts[-length(parts)]) return(FALSE)
  parent <- if (length(parts) > 1L) parts[length(parts) - 1L] else ""
  if (nzchar(parent) &&
      paste(parent, stem, sep = "/") %in% get0("FIGURE_DROP", ifnotfound = character(0)))
    return(FALSE)
  if (stem %in% keep) return(TRUE)
  # Allow a horizon-suffixed variant of an unsuffixed allow-list entry, and a
  # date-stamped variant (the delay-fit figures carry a _YYYYMMDD snapshot tag).
  base_stem <- sub("_h[0-9]+$", "", stem)
  base_stem <- sub("_[0-9]{8}$", "", base_stem)
  base_stem %in% keep
}

# ---------------------------------------------------------------------------
# Package requirements (checked in run_all.R)
# ---------------------------------------------------------------------------
REQUIRED_PACKAGES <- c(
  # Core
  "tidyverse", "data.table", "lubridate", "here",
  # Spatial
  "sf", "spdep", "spatstat.geom",
  # Models
  # surveillance (hhh4), pscl (ZINB) and deSolve (SEIR ODE solver) were REMOVED with the
  # comparators that used them; none has a single call left in the tree. Declaring a hard
  # dependency the code never loads makes a fresh environment install packages it does not
  # need, and makes a missing-package failure look like a real one.
  "EpiNow2",         # R(t) estimation (NOT EpiEstim)
  "MASS",            # glm.nb for the gravity calibration (M4)
  # Visualisation
  "ggplot2", "patchwork", "viridis", "ggridges", "cowplot",
  "scales", "RColorBrewer",
  # Tests
  "testthat",
  # Utilities
  "jsonlite", "readxl",
  # Parallelism (optional; only needed when PARALLEL_JOBS > 1)
  "future", "furrr"
)

# ---------------------------------------------------------------------------
# Utility: create output dirs on first source
# ---------------------------------------------------------------------------
invisible(lapply(
  c(OUT_FORECASTS, OUT_MOBILITY, OUT_DIAGNOSTICS, OUT_MAPS,
    OUT_CALIBRATION, OUT_REPORTS, MODELS_DIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

message("[config] 00_config.R loaded — analysis date: ", ANALYSIS_DATE)
