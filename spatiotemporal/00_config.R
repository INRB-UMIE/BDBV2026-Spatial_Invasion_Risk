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

# Current outbreak start: earliest sample collection date
OUTBREAK_START     <- as.Date("2026-04-30")

# Forecast TARGET is P(first symptom ONSET of a confirmed case in the next 1-2
# weeks) — an onset-dated event, since transmission tracks onset, not reporting.
# ~15% of records lack an onset date; when TRUE we impute their onset as
# (sample_date - Delta), where Delta is DRAWN PER RECORD from the onset->sample delay
# distribution (empirical bootstrap of the observed complete pairs when >=30 exist,
# else a parametric Exp(rate) draw) rather than fixed at its mean — so imputed onsets
# contribute a true ONSET week spread by the delay's shape rather than piling on one
# week. Set FALSE to fall back to sample date verbatim.
IMPUTE_ONSET_FROM_SAMPLE <- TRUE

# Onset-handling MODE for confirmed cases lacking a usable onset date (review §1.1/§1.3).
# Supersedes IMPUTE_ONSET_FROM_SAMPLE when set (which is kept for backward compatibility:
# unset ONSET_MODE derives "impute" when IMPUTE_ONSET_FROM_SAMPLE=TRUE else "sample_verbatim").
#   "impute"        : stochastic backward draw onset = sample - Delta (default; unchanged behaviour).
#   "growth_impute" : backward draw TILTED by the epidemic growth rate r (Delta weighted by
#                     exp(-r*Delta)), so onsets are not systematically pushed too early during
#                     growth — the reviewer's non-reversibility point (a growth-adjusted imputation).
#   "complete_case" : DROP onset-less confirmed records (the reviewer's recommended robustness
#                     analysis; pair with the nowcast). NOTE: this changes the invasion event set.
#   "sample_verbatim": date the record at its sample date (old IMPUTE_ONSET_FROM_SAMPLE=FALSE path).
# The three non-default modes are for the imputation-robustness comparison (rank stability across
# modes), not the primary run.
ONSET_MODE <- get0("ONSET_MODE", ifnotfound = "impute")
ONSET_GROWTH_WINDOW_WEEKS <- 8L   # recent-weeks window for the log-linear growth-rate estimate

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
DELAY_ONSET_SAMPLE_FAMILY <- "exp"
DELAY_ONSET_SAMPLE_RATE   <- 0.228     # day^-1 → mean = 4.39 d (lab CSV implied_mean_fit)

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
# NOTE (01_data_prep.R): whenever >=30 complete onset+sample pairs exist, the imputation
# ALWAYS draws from an EMPIRICAL bootstrap of the CURRENT line list's WINDOWED pairs (onset in
# [outbreak floor, max(sample) - DELAY_TRUNC_BUFFER_DAYS]) regardless of this setting — the
# most faithful delay. This constant chooses the PARAMETRIC FALLBACK (used only when <30 pairs)
# AND which rate is reported:
#   "data" (default here — this IS the DHIS2 pipeline) — the rigorous DHIS2 interval-censored
#           best-family fit from 04c_dhis2_delay_windows.R (dhis2_onset_sample_delay_params.csv)
#           when present, else the windowed Exponential MLE (1/mean) of the current pairs.
#   "lab"  — the fixed lab-linelist Exp(rate) above (DELAY_ONSET_SAMPLE_RATE). Use only when
#           deliberately imputing DHIS2 onsets with the (faster) lab delay.
ONSET_SAMPLE_DELAY_SOURCE <- "data"

# ---------------------------------------------------------------------------
# Delay estimation: Bayesian EpiDist (truncation- AND double-interval-censoring
# corrected) as the DEFAULT delay estimator.  (04c_dhis2_delay_windows.R)
# ---------------------------------------------------------------------------
# The onset->sample delay drives (a) the missing-onset imputation (~15% of cases) and
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

# Extra Bayesian PROFILE analyses (refit the featured model ~10x each): a loo-predictive
# posterior over the generation-time mean, and a nowcast-input sensitivity of beta0. Set FALSE
# to skip them (they add several minutes to a run); the core suite/forecasts are unaffected.
BAYES_PROFILE_ANALYSES <- TRUE

# BUMP whenever estimate_rt_epinow2()'s model spec changes (delays/truncation/priors/GT), so the
# on-disk R(t) cache (keyed by gt + analysis_date + this version + GT params) is invalidated
# rather than silently serving a stale rds computed under the old spec. v2 = right-truncation model;
# v3 = summarise at CrIs c(0.6, 0.9) so R(t) reports the 60% (q20/q80) and 90% (q5/q95) bands.
RT_CACHE_VERSION <- 3L

# Onset-to-death delay: Gamma, fit to BDBV-2026 Ituri lab data (n=77)
DELAY_ONSET_DEATH_FAMILY  <- "gamma"
DELAY_ONSET_DEATH_SHAPE   <- 1.8914
DELAY_ONSET_DEATH_RATE    <- 0.2698    # → mean = shape/rate ≈ 7.01 d

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
GT_PROFILES <- list(
  short = list(
    label = "GT-Short (12.0 d)",
    mean  = 12.0,  # low-end sensitivity: just below the pooled 95% CI lower bound (13.2 d, Nash 2024)
    sd    = 6.5,   # CV ~= 0.54
    max_tau = 35
  ),
  medium = list(
    label = "GT-Medium (15.3 d)",
    mean  = 15.3,  # Zaire ebolavirus serial interval, WHO Ebola Response Team 2014 NEJM (proxy; SI ~= GT)
    sd    = 9.3,   # WHO Ebola Response Team 2014 NEJM (SD of the serial interval); CV ~= 0.61
    max_tau = 45
  ),
  long = list(
    label = "GT-Long (18.0 d)",
    mean  = 18.0,  # high-end sensitivity: just above the pooled 95% CI upper bound (17.5 d, Nash 2024)
    sd    = 10.5,  # CV ~= 0.58
    max_tau = 50
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
GT_PRIOR <- list(
  mean_mu     = 15.3,  mean_sd = 1.82,   # mu_GT ~ Normal(mean_mu, mean_sd)
  sd_mu       = 9.3,   sd_sd   = 1.50,   # sd_GT ~ Normal(sd_mu,  sd_sd)
  mean_bounds = c(10, 21),               # truncation support for the GT mean (days)
  sd_bounds   = c(4, 14),                # truncation support for the GT sd  (days)
  max_tau     = 45L,                     # daily PMF support (days); matches medium profile
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
# prior". Members (all medium anchor, no covariates): M4 gravity, M10 radiation
# composite, M8 short-trip+gravity, M15 combined-static, M17 all-kernel consensus.
BAYES_ENSEMBLE_KERNELS <- c("M4", "M10", "M8", "M15", "M17")

# Fraction of residents' time spent in their home zone, for the manuscript-
# motivated inward / meeting-location force-of-infection variant (Mills 2026).
# A documented modelling assumption (commuting fractions 18-40% in that work =>
# home fractions 0.60-0.82); the central 0.70 is used, sensitivity-analysable.
MOBILITY_HOME_FRACTION <- 0.70

# BDBV incubation period (for SEIR model compartments)
# Source: Wamala et al. 2010 / MacNeil et al. 2010 Emerg Infect Dis — BDBV 2007
# Uganda incubation ~5.7 d (survivors) to 7.4 d (fatal). 5.3 d is a modelling
# choice near the lower survivor estimate; Nash et al. 2024 pool EVD incubation at
# 8.5 d, so this is a conservative (short) latent period for the SEIR comparator.
INCUBATION_MEAN_DAYS <- 5.3
INCUBATION_RATE      <- 1 / INCUBATION_MEAN_DAYS   # σ

# Nominal ascertainment (overall test positivity) and sensitivity range.
# Source: positivity_aggregated.csv (Ituri, DHIS2). Ascertainment is treated as
# UNCERTAIN, not a fixed constant: the PRIMARY invasion forecast p_case is defined
# on confirmed-case onset and is ascertainment-agnostic (it does not divide by rho);
# the secondary infection-scale score p_infection is reported as a BAND across the
# grid below (rho=0.60 -> lower bound, 0.30 -> upper bound), never a single point.
# Nominal 0.45 is the centre of the grid and matches the observed Ituri rolling test
# positivity (~0.45-0.49 in positivity_aggregated.csv); the grid brackets it +/-0.15.
ASCERTAINMENT_NOMINAL <- 0.45
ASCERTAINMENT_GRID    <- c(0.30, 0.45, 0.60)

# Provinces for which within-province relative-risk scores and maps are produced
# (in addition to the nationwide scores). Strings must match the shapefile PROVINCE
# field exactly. Ituri is the outbreak epicentre; Nord-Kivu and Haut-Uele are the
# highest-exposure neighbours.
PROVINCES_OF_INTEREST <- c("Ituri", "Nord-Kivu", "Haut-Uele")

# Delay-adjusted cCFR — in the SEIR comparator it PARTITIONS the (fixed ~8 d) infectious
# outflow into death (cCFR/mean_illness) vs recovery ((1-cCFR)/mean_illness); it does NOT
# set the infectious period itself.
# Source: cfr_summary.csv; dhis2_all, confirmed, delay-adjusted
CFR_POINT_ESTIMATE <- 0.264

# ---------------------------------------------------------------------------
# Mobility matrix identifiers
# ---------------------------------------------------------------------------
MOBILITY_IDS <- c("M1", "M2a", "M2b", "M3", "M4", "M4b", "M5", "M6a", "M6b",
                  "M7", "M8", "M9", "M10",   # M11 (inward FOI) built on demand in run_all
                  "M13", "M14",              # Flowminder-cohort composites (gravity / radiation, travel-time)
                  "M15", "M16", "M17")       # Flowminder combined inflow+outflow static; cohort+static; all-kernel ensemble
# OSRM ROAD-DISTANCE (km) deterrence variants of the mobility kernels, built on demand when the
# osrm road-distance matrix is present (build_all_mobility_matrices' osrm_dist_mat arg; M11-dist
# in run_all). Same construction as their base kernels but keyed on km instead of travel time.
# M13-dist / M14-dist are the geographic-distance analogues of the cohort composites; M17-dist is
# the road-km analogue of the grand all-kernel consensus ensemble (M15/M16 are empirical Flowminder
# flows and are distance-agnostic, so they have no -dist twin).
MOBILITY_DIST_IDS <- c("M4-dist", "M8-dist", "M9-dist", "M10-dist", "M11-dist",
                       "M13-dist", "M14-dist", "M17-dist")
MOBILITY_PRIMARY <- "M8"    # recommended composite

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
# Flowminder combined inflow+outflow static kernel + composites (M15/M16/M17)
# ---------------------------------------------------------------------------
# THREE new mobility kernels drawn from the full Flowminder origin-destination data:
#   * M15  Flowminder STATIC kernel built from BOTH the inflow AND outflow OD tables. The
#          processed inflow/outflow exports are the SAME directed origin->destination table
#          (Flowminder ships one OD matrix under both labels), so the information-adding
#          combination is the SYMMETRISED total two-way flow S = O + t(O): a zone's connection
#          to another is trips OUT to it PLUS trips IN from it. This is genuinely distinct from
#          M3 (directed outflow only) and fills reciprocal edges, so it is materially LESS sparse.
#   * M16  Flowminder COHORT + STATIC composite. Cohort subscriber-day presence rows (Ituri/NK/
#          Tshopo origins) where available, and the M15 static flows elsewhere — the direct
#          analogue of M13 (cohort + gravity) and M14 (cohort + radiation), but pairing the cohort
#          data with the empirical Flowminder static flows instead of a gravity/radiation model.
#          "Static used only where the cohort data are unavailable", as with cohort+gravity/radiation.
#   * M17  ALL-KERNEL CONSENSUS ensemble. Element-wise mean of every distinct STRUCTURAL/empirical
#          base hypothesis {Flowminder directed OD (M3), gravity (M4), radiation (M5), travel-time
#          decay (M6a), Flowminder static (M15)} — a convex combination, so row-stochastic — with
#          the cohort/epicentre source rows overlaid (best-available empirical data). M17-dist is
#          the road-km analogue (gravity/radiation/decay re-keyed on km; the empirical M3/M15 reused).
# M15/M16/M17 are ADDITIVE — every existing kernel/model is untouched. The MODEL variants that use
# them are gated by INCLUDE_FLOWSTATIC_MODELS (below); M16 additionally needs the cohort kernel
# (INCLUDE_COHORT_MODELS), and M17-dist needs the OSRM road-distance matrix. "M12" stays reserved
# for the effective-distance kernel on the parallel branch, so these ids start at M15.

# Human-readable mobility-kernel labels for figures/tables (id -> description).
MOBILITY_LABELS <- c(
  M1  = "Short-trip epicentre",
  M2a = "Short-trip epicentre (avg)",
  M2b = "Short-trip epicentre (latest)",
  M3  = "Flowminder OD (national)",
  M4  = "Gravity (travel-time)",
  M4b = "Gravity (exp deterrence)",
  M5  = "Radiation (travel-time)",
  M6a = "Travel-time decay (exp)",
  M6b = "Travel-time decay (power)",
  M7  = "IDP-augmented OD",
  M8  = "Short-trip + gravity (travel-time)",
  M9  = "Short-trip + kernel ensemble",
  M10 = "Short-trip + radiation (travel-time)",
  M11 = "Inward/meeting-location FOI",
  M13 = "Flowminder cohort + gravity (travel-time)",
  M14 = "Flowminder cohort + radiation (travel-time)",
  M15 = "Flowminder static (inflow+outflow)",
  M16 = "Flowminder cohort + directed OD",
  M17 = "All-kernel consensus ensemble (travel-time)",
  `M4-dist`  = "Gravity (road-km)",
  `M8-dist`  = "Short-trip + gravity (road-km)",
  `M9-dist`  = "Short-trip + kernel ensemble (road-km)",
  `M10-dist` = "Short-trip + radiation (road-km)",
  `M11-dist` = "Inward/meeting-location FOI (road-km)",
  `M13-dist` = "Flowminder cohort + gravity (road-km)",
  `M14-dist` = "Flowminder cohort + radiation (road-km)",
  `M17-dist` = "All-kernel consensus ensemble (road-km)"
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
# Model-suite composition toggles (comprehensive modelling suite)
# ---------------------------------------------------------------------------
# Each flag switches an OPTIONAL family of model variants on or off WITHOUT touching
# the CORE suite, which is always fitted: the base mobility-informed renewal kernels
# (M4/M8/M9/M10) at the medium generation-time anchor (intercept-only import hazard —
# the featured-model class). (Generation time is NO LONGER a grid axis — it is
# marginalised over GT_PRIOR in the featured forecast; see 21_bayesian_renewal.R.)
# The defaults below configure a RUNTIME-TRIMMED suite: the nine mobility kernels
# M4/M8/M9/M10 (base), M13/M14 (Flowminder cohort) and M15/M16/M17 (Flowminder
# static/consensus), each swept over covariates {none, geo} on a SINGLE distance
# measure (OSRM travel time). The road-distance (-dist) twins and the FULL covariate
# models are OFF by default (explore them on the best-fitting model instead), which
# keeps the selection grid small. Set a flag TRUE/FALSE here, or export an environment
# variable of the same name, to change it. bayes_default_grid() and run_all.R read these
# via get0()/the flag itself, so a context that sources neither still sees the defaults.
#   INCLUDE_OSRM_DIST_MODELS     road-distance (-dist) deterrence twins — the MASTER switch for
#                                EVERY -dist variant in the Bayesian grid (generic M4/M8/M9/M10-dist
#                                AND cohort M13/M14-dist AND consensus M17-dist); OFF by default
#   INCLUDE_GEO_COV_MODELS       the reduced "geo" covariate set (log_pop, CCVI, d_min); ON by default
#   INCLUDE_FULL_COV_MODELS      the FULL-exogenous covariate set (geo + healthsite_density); OFF by
#                                default (the extra covariate roughly halves the covariate sweep)
#   INCLUDE_M11_MODELS           the inward / meeting-location FOI kernel (M11) and its variants; OFF
#   INCLUDE_LOGIT_SENS_MODELS    the (slow-to-fit) logit-link observation-process sensitivity; OFF
#   INCLUDE_SUSPECTED_COV_MODELS the suspected-but-not-confirmed leading-indicator covariate
#                                models (own preceding-week + mobility-weighted import); OFF by default
#   INCLUDE_FLOWSTATIC_MODELS    the Flowminder combined inflow+outflow static family — the M15 static
#                                kernel, the M16 cohort+static composite, and the M17/M17-dist grand
#                                all-kernel consensus ensemble; ON by default
# Structural-baseline toggles (BASELINE_MODELS in run_all.R — always-on yardsticks, decoupled from
# RUN_FREQUENTIST_MODELS). The cheap deterministic baselines (Gravity-B4, Distance-B1, Adjacency-B7)
# always run; the two EXPENSIVE per-fold mechanistic ones are gated so they can be skipped for speed:
#   INCLUDE_SEIR_BASELINE        the stochastic spatial SEIR metapopulation comparator (C6); OFF by default
#   INCLUDE_HHH4_BASELINE        the endemic-epidemic hhh4 comparator; OFF by default
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
# OSRM road-distance (km) deterrence kernels: ON by default so the suite considers BOTH the
# travel-time and the geographic-distance deterrence axis for every kernel (roughly doubles the
# LFO-CV grid). Set INCLUDE_OSRM_DIST_MODELS=FALSE to fall back to travel-time only.
INCLUDE_OSRM_DIST_MODELS     <- .model_suite_flag("INCLUDE_OSRM_DIST_MODELS",     TRUE)
INCLUDE_GEO_COV_MODELS       <- .model_suite_flag("INCLUDE_GEO_COV_MODELS",       TRUE)
INCLUDE_FULL_COV_MODELS      <- .model_suite_flag("INCLUDE_FULL_COV_MODELS",      FALSE)
INCLUDE_M11_MODELS           <- .model_suite_flag("INCLUDE_M11_MODELS",           FALSE)
INCLUDE_LOGIT_SENS_MODELS    <- .model_suite_flag("INCLUDE_LOGIT_SENS_MODELS",    FALSE)
INCLUDE_SUSPECTED_COV_MODELS <- .model_suite_flag("INCLUDE_SUSPECTED_COV_MODELS", FALSE)
INCLUDE_FLOWSTATIC_MODELS    <- .model_suite_flag("INCLUDE_FLOWSTATIC_MODELS",    TRUE)
# M9 (short-trip + M4/M5/M6a ensemble) and M15 (symmetrised inflow+outflow static) are OFF by
# default: M9 was consistently beaten by the simpler composites and M15 was among the weakest
# (spiky) kernels, so neither is fit unless explicitly re-enabled. Turning these on also builds
# the respective mobility matrices; when off they are neither built nor added to the grid.
INCLUDE_M9_MODELS            <- .model_suite_flag("INCLUDE_M9_MODELS",             FALSE)
INCLUDE_M15_MODELS           <- .model_suite_flag("INCLUDE_M15_MODELS",            FALSE)
INCLUDE_SEIR_BASELINE        <- .model_suite_flag("INCLUDE_SEIR_BASELINE",        FALSE)
INCLUDE_HHH4_BASELINE        <- .model_suite_flag("INCLUDE_HHH4_BASELINE",        FALSE)
message(sprintf("[config] Optional model families — OSRM-dist:%s geo-cov:%s full-cov:%s M9:%s M11:%s M15:%s logit:%s susp:%s flowstatic:%s | baselines SEIR:%s hhh4:%s",
                INCLUDE_OSRM_DIST_MODELS, INCLUDE_GEO_COV_MODELS, INCLUDE_FULL_COV_MODELS,
                INCLUDE_M9_MODELS, INCLUDE_M11_MODELS, INCLUDE_M15_MODELS,
                INCLUDE_LOGIT_SENS_MODELS, INCLUDE_SUSPECTED_COV_MODELS,
                INCLUDE_FLOWSTATIC_MODELS, INCLUDE_SEIR_BASELINE, INCLUDE_HHH4_BASELINE))

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

# Epicentre zones used as origins in the short-trip cohort (pooled)
EPICENTRE_ZONES <- c("Bunia", "Mongbalu", "Rwampara")

# ---------------------------------------------------------------------------
# Evaluation parameters
# ---------------------------------------------------------------------------
LFO_MIN_TRAINING_WEEKS  <- 4      # minimum weeks before first LFO-CV fold
LFO_HORIZONS            <- c(1L, 2L)    # 1-week ahead and 2-week ahead
N_SIMULATIONS           <- 2000L  # stochastic model simulations per forecast
RANDOM_SEED             <- 20260704L

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
# Package requirements (checked in run_all.R)
# ---------------------------------------------------------------------------
REQUIRED_PACKAGES <- c(
  # Core
  "tidyverse", "data.table", "lubridate", "here",
  # Spatial
  "sf", "spdep", "spatstat.geom",
  # Models
  "surveillance",    # hhh4
  "EpiNow2",         # R(t) estimation (NOT EpiEstim)
  "MASS",            # glm.nb for gravity calibration
  "pscl",            # ZINB (zeroinfl)
  "deSolve",         # ODE solver for SEIR
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
