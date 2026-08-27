# =============================================================================
# run_all.R — Master Orchestration Script
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Usage: Rscript spatiotemporal/run_all.R
#   or:  source("spatiotemporal/run_all.R")  within an interactive session
#
# Order of operations:
#   00  Config and package checks
#   01  Data preparation
#   02  Epidemiological parameters (EpiNow2 R(t), GT PMFs, delays)
#   03  Mobility matrices (M1–M8)
#   04  Nowcast correction
#   05–09 Model fitting and forecasting
#   10  LFO-CV evaluation
#   11  Metric computation
#   12  Calibration
#   13  Visualisations
#   14  Methodology report
# =============================================================================

# ---------------------------------------------------------------------------
# 0. Config
# ---------------------------------------------------------------------------
source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(here)
})

# Check required packages
missing_pkgs <- REQUIRED_PACKAGES[!vapply(REQUIRED_PACKAGES,
  requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  warning("Missing packages (some models may be unavailable): ",
          paste(missing_pkgs, collapse = ", "))
}

# Parallel worker pool for independent Bayesian fits (opt-in via PARALLEL_JOBS;
# see 00_config.R). multicore forks on Linux, so the large shared objects (mobility
# matrices, covariates) are shared copy-on-write; a no-op (sequential) when
# PARALLEL_JOBS = 1 or `future` is unavailable.
if (isTRUE(get0("PARALLEL_JOBS", ifnotfound = 1L) > 1L) &&
    requireNamespace("future", quietly = TRUE)) {
  future::plan(future::multicore, workers = PARALLEL_JOBS)
  # furrr ships each LFO fold's fitted-model list to the workers as a captured
  # "global". That list is ~2 GiB, far over future's 500 MiB default ceiling —
  # which SILENTLY aborted a prior LFO (run_invasion_lfo tryCatch -> NULL ->
  # every downstream product left stale). Raise the ceiling generously. This is a
  # guard THRESHOLD, not an allocation, and multicore FORKS, so the globals are
  # shared copy-on-write: the higher limit costs no additional memory.
  options(future.globals.maxSize = 12 * 1024^3)  # 12 GiB export ceiling (guard only)
  message(sprintf("[parallel] future multicore plan: %d workers (each Bayesian fit uses 2 chains); future.globals.maxSize=12 GiB", PARALLEL_JOBS))
} else {
  message("[parallel] PARALLEL_JOBS = 1 (or future unavailable) — sequential fitting.")
}

# Lightweight phase timing (diagnostic): minutes elapsed per major phase.
.PH_T0 <- Sys.time(); .PH_LAST <- .PH_T0
.phase <- function(lbl) {
  now <- Sys.time()
  message(sprintf("[timing] %-42s %6.1f min  (cum %5.1f)", lbl,
                  as.numeric(difftime(now, .PH_LAST, units = "mins")),
                  as.numeric(difftime(now, .PH_T0, units = "mins"))))
  .PH_LAST <<- now
}

# ---------------------------------------------------------------------------
# 1. Data preparation
# ---------------------------------------------------------------------------
# 1a. Refresh the DHIS2-specific onset->sample delay (windowed interval-censored MLE)
#     BEFORE data prep, so the onset imputation (01_data_prep.R) draws its parametric
#     fallback / reported rate from the current line list's own delay rather than the
#     faster lab reference. Sourcing defines the estimator functions only (the script's
#     main block is Rscript-gated); the fit is fast (~5 s) and non-fatal — if it fails the
#     imputation falls back to the windowed empirical bootstrap + Exp(1/mean). Set
#     RUN_EPIDIST=TRUE to additionally run the Bayesian truncation-corrected fit standalone.
message("\n=== Step 1a: DHIS2 onset->sample delay refresh ===")
source(file.path(ST_DIR, "04c_dhis2_delay_windows.R"))   # estimate_dhis2_onset_sample_delay(), estimate_dhis2_onset_sample_epidist(), write_onset_sample_long()
tryCatch({
  .dmeta   <- jsonlite::fromJSON(LINELIST_JSON)
  .dll_csv <- file.path(LINELIST_DIR, .dmeta$folder, "dhis2_processed_linelist.csv")
  .dll_raw <- readr::read_csv(.dll_csv, col_types = readr::cols(.default = "c"),
                              show_col_types = FALSE, na = c("", "NA", "N/A"))
  .osfit   <- estimate_dhis2_onset_sample_delay(.dll_raw, analysis_date = ANALYSIS_DATE)
  # Bayesian EpiDist MARGINAL delay (truncation + double-interval-censoring corrected) — the
  # DEFAULT estimator (RUN_EPIDIST, 00_config.R). Fit ONCE here on the same windowed pairs;
  # write_onset_sample_long() appends its (family, mean, sd) as epidist_* rows so
  # .load_dhis2_delay_params() (01_data_prep.R) PREFERS it for the onset imputation AND the
  # nowcast rate. NULL (no epidist_* rows; censored-MLE used) when RUN_EPIDIST=FALSE, the
  # package is absent, or the fit fails — the run stays non-fatal either way. This is a
  # ONE-TIME data-prep fit (a handful of Stan fits), NOT repeated per LFO fold or per model.
  .osepi   <- estimate_dhis2_onset_sample_epidist(.dll_raw, analysis_date = ANALYSIS_DATE)
  if (!is.null(.osfit)) {
    write_onset_sample_long(.osfit,
      file.path(DATA_DIR, "cfr_reference", "dhis2_onset_sample_delay_params.csv"),
      source_label = .dmeta$folder, epidist = .osepi)
    write_onset_sample_long(.osfit,
      file.path(LINELIST_DIR, .dmeta$folder, "dhis2_onset_sample_delay_params.csv"),
      source_label = .dmeta$folder, epidist = .osepi)
    message(sprintf("[run_all] DHIS2 onset->sample delay: %s, mean %.2f d, rate %.3f/d (n=%d, window %s)%s",
                    .osfit$best_family, .osfit$implied_mean, .osfit$rate, .osfit$n_fit, .osfit$window,
                    if (!is.null(.osepi))
                      sprintf(" | EpiDist marginal PREFERRED: %s mean %.2f d sd %.2f d",
                              .osepi$family, .osepi$mean, .osepi$sd)
                    else " | EpiDist marginal: not used (censored-MLE)"))
  } else message("[run_all] DHIS2 delay fit returned NULL (too few pairs); imputation uses the windowed bootstrap.")
}, error = function(e) message("[run_all] DHIS2 delay refresh skipped (non-fatal): ", conditionMessage(e)))

message("\n=== Step 1: Data preparation ===")
source(file.path(ST_DIR, "01_data_prep.R"))
dat <- prep_all_data()

zone_week_nc <- dat$zone_week   # all zone-week counts (will be nowcast-corrected in step 4)
zones_all    <- names(dat$pop)  # canonical zone name vector (519 zones, from population spine)
pop_vec      <- dat$pop
covariates   <- dat$covariates
sitrep       <- dat$sitrep
contacts     <- dat$contacts

# ---------------------------------------------------------------------------
# 2. Epidemiological parameters
# ---------------------------------------------------------------------------
.phase("Step 1   data prep")
message("\n=== Step 2: Epidemiological parameters ===")
source(file.path(ST_DIR, "02_epi_params.R"))

# GT PMFs for all profiles
gt_pmfs <- compute_all_gt_pmfs()           # named list: short/medium/long → PMF vector

# EpiNow2 R(t) estimation (national-level, per GT profile). DIAGNOSTIC OUTPUT ONLY —
# it does NOT feed any forecast: the Bayesian invasion model rolls source-zone incidence
# forward with its own internally-estimated LOCAL renewal R (estimate_R_local() in
# 21_bayesian_renewal.R), and the (frequentist) SEIR uses a fixed R. It is therefore gated
# behind RUN_RT_ESTIMATION (default FALSE): EpiNow2 is a separate MCMC fit per GT profile
# and one of the slowest steps, so it is skipped unless the R(t) diagnostic figure is
# wanted. Export RUN_RT_ESTIMATION=TRUE to produce it. When skipped, Rt_primary/Rt_scalar
# stay NULL — plot_rt() no-ops on NULL and the run bundle simply carries Rt = NULL.
RUN_RT_ESTIMATION <- {
  .env <- tolower(trimws(Sys.getenv("RUN_RT_ESTIMATION", "")))
  if (nzchar(.env)) .env %in% c("true", "t", "1", "yes", "y")
  else isTRUE(get0("RUN_RT_ESTIMATION", ifnotfound = FALSE))
}
Rt_national_list <- list()
Rt_primary <- NULL
Rt_scalar  <- NULL
Rt_zone    <- NULL   # zone-level R(t): optional, never populated in the default path
if (isTRUE(RUN_RT_ESTIMATION)) {
  # Count CONFIRMED cases only (dat$ll retains all classifications, so an unfiltered count
  # would be "all alerts", inflating R(t)); date on onset, coalescing to sample only when
  # onset is missing. NB: the recent tail is right-truncated (onsets not yet reported), so
  # the last ~2 weeks of R(t) are biased low — read as indicative, not operational.
  daily_cases <- dat$ll |>
    filter(confirmed, !is.na(date_of_symptom_onset) | !is.na(date_of_sample_collection)) |>
    mutate(date = coalesce(date_of_symptom_onset, date_of_sample_collection)) |>
    count(date, name = "confirm") |>
    arrange(date) |>
    filter(date >= OUTBREAK_START, date <= ANALYSIS_DATE)
  for (gt_name in names(gt_pmfs)) {
    message(sprintf("[Rt] Estimating national R(t) with GT profile: %s", gt_name))
    gt_p <- GT_PROFILES[[gt_name]]
    Rt_national_list[[gt_name]] <- tryCatch(
      estimate_rt_epinow2(daily_cases, gt_name, gt_p, ANALYSIS_DATE),
      error = function(e) { warning("EpiNow2 failed for ", gt_name, ": ", e$message); NULL }
    )
  }
  # Use primary GT profile R(t)
  Rt_primary <- Rt_national_list[[GT_PRIMARY]]
  if (!is.null(Rt_primary)) {
    Rt_scalar <- list(
      mean = tail(Rt_primary$R_mean, 1),
      sd   = (tail(Rt_primary$R_hi_90, 1) - tail(Rt_primary$R_lo_90, 1)) / (2 * 1.645)
    )
    message(sprintf("[Rt] Primary R(t) estimate: %.2f (90%% CI: %.2f–%.2f)",
                    Rt_scalar$mean, tail(Rt_primary$R_lo_90, 1), tail(Rt_primary$R_hi_90, 1)))
  } else {
    Rt_scalar <- list(mean = 1.5, sd = 0.3)  # fallback prior
    warning("[Rt] EpiNow2 failed; using prior R(t) = 1.5")
  }
} else {
  message("[Rt] RUN_RT_ESTIMATION = FALSE — skipping EpiNow2 R(t) (diagnostic only; feeds no forecast).")
}

# Serial interval from contact tracing (optional)
si_result <- if (!is.null(contacts) && nrow(contacts) > 0 && !is.null(dat$ll)) {
  tryCatch(estimate_si_from_contacts(contacts, dat$ll),
           error = function(e) NULL)
} else NULL

# ---------------------------------------------------------------------------
# 3. Mobility matrices
# ---------------------------------------------------------------------------
.phase("Step 2   epi params + EpiNow2 R(t)")
message("\n=== Step 3: Mobility matrices ===")
source(file.path(ST_DIR, "03_mobility_matrices.R"))

osrm_mat <- load_osrm()                          # travel time (minutes) — base kernels
# OSRM road-distance (km) — drives the -dist mobility variants (M4/M8/M9/M10/M11-dist).
# Optional: if the matrix is absent the -dist variants are simply not built and the
# Bayesian grid filters them out (bayes_default_grid keeps only available kernels).
osrm_dist_mat <- tryCatch(load_osrm("road_distance"), error = function(e) {
  message("[mobility] OSRM road-distance matrix unavailable (", conditionMessage(e),
          "); road-distance (-dist) variants will be skipped."); NULL })
# Reuse the canonical population spine (dat$pop == load_population(); deterministic) so the
# base kernels here and the M11 / M11-dist inward matrices built later (which use pop_vec)
# share ONE population vector — no reliance on two independent loads agreeing.
pop_for_mob <- pop_vec

# Subset to zones_all for efficiency
shared_zones <- intersect(zones_all, rownames(osrm_mat))
message(sprintf("[mobility] %d zones with OSRM coverage (of %d total)",
                length(shared_zones), length(zones_all)))

# OSRM road-distance (-dist) kernels are OPTIONAL (INCLUDE_OSRM_DIST_MODELS, default FALSE).
# The road-distance matrix is supplied when EITHER the generic -dist family OR the cohort
# geographic-distance composites (M13/M14-dist, INCLUDE_COHORT_MODELS, default TRUE) are wanted;
# build_all_mobility_matrices gates the GENERIC -dist variants separately (on INCLUDE_OSRM_DIST_MODELS)
# so passing the matrix for the cohort composites does not switch the generic -dist family on.
mobility_matrices <- build_all_mobility_matrices(
  zones_all   = zones_all,
  pop_vec     = pop_for_mob,
  osrm_mat    = osrm_mat,
  analysis_date = ANALYSIS_DATE,
  osrm_dist_mat = if (isTRUE(INCLUDE_OSRM_DIST_MODELS) || isTRUE(INCLUDE_COHORT_MODELS) ||
                      isTRUE(INCLUDE_FLOWSTATIC_MODELS))
    osrm_dist_mat else NULL
)
message("[mobility] Built matrices: ", paste(names(mobility_matrices), collapse=", "))

# Primary mobility matrix for main analysis
W_primary <- mobility_matrices[[MOBILITY_PRIMARY]]

# ---------------------------------------------------------------------------
# 4. Nowcast correction
# ---------------------------------------------------------------------------
.phase("Step 3   mobility matrices")
message("\n=== Step 4: Nowcast correction ===")
source(file.path(ST_DIR, "04_nowcasting.R"))
source(file.path(ST_DIR, "04b_epinowcast.R"))

# Primary current-week nowcast: full probabilistic epinowcast (Bayesian delay +
# NB observation + weekly random walk), falling back to the deterministic
# Exp-CDF correction if epinowcast/cmdstanr is unavailable or the fit fails.
# NOTE: LFO-CV (step 10) intentionally keeps the fast deterministic per-fold
# nowcast — refitting epinowcast at each of the ~11 folds would be prohibitively
# slow, and for retrospective scoring the deterministic correction is adequate.
zone_week_nc <- nowcast_zone_week_epinowcast(
  zone_week      = zone_week_nc,
  linelist       = dat$ll,
  analysis_date  = ANALYSIS_DATE,
  outbreak_start = OUTBREAK_START
)
message(sprintf("[nowcast] Method: %s; mean weight: %.3f",
                attr(zone_week_nc, "nowcast_method") %||% "deterministic",
                mean(zone_week_nc$trunc_weight, na.rm = TRUE)))

# Persist the epinowcast weekly correction factors for reporting, if produced.
enw_factors <- attr(zone_week_nc, "epinowcast_factors")
if (!is.null(enw_factors)) {
  readr::write_csv(enw_factors,
                   file.path(OUT_DIAGNOSTICS, "epinowcast_weekly_factors.csv"))
}

# ---------------------------------------------------------------------------
# 5–9. Model fitting and current-week forecasts
# ---------------------------------------------------------------------------
source(file.path(ST_DIR, "05_baseline_models.R"))   # forecast_B1/B4, zone_week_to_wide
source(file.path(ST_DIR, "06_simple_models.R"))     # compute_foi, daily_to_weekly_gt
source(file.path(ST_DIR, "07_hhh4_model.R"))        # hhh4 endemic-epidemic
source(file.path(ST_DIR, "08_stochastic_seir.R"))   # stochastic SEIR (C6)
source(file.path(ST_DIR, "15_workhorse.R"))         # mobility-informed renewal workhorse
source(file.path(ST_DIR, "16_invasion_eval.R"))     # invasion LFO + evaluation
source(file.path(ST_DIR, "18_ensemble.R"))          # Q3 mean/median invasion ensembles
source(file.path(ST_DIR, "19_spacetime_eval.R"))    # Q5 spatiotemporal evaluation
source(file.path(ST_DIR, "20_forecast_detail.R"))   # spec/params/priority/uncertainty viz + build_invasion_design
source(file.path(ST_DIR, "21_bayesian_renewal.R"))  # Bayesian (brms) renewal invasion suite
source(file.path(ST_DIR, "22_daily_reissue.R"))     # daily re-issue dating/persistence + intra-week backtest
source(file.path(ST_DIR, "23_prospective_eval.R"))  # prospective + forecast-vs-observed eval (review §3.4/3.5/3.6/3.8)
source(file.path(ST_DIR, "24_review_figures.R"))    # publication figures for the review analyses (house style)

province_map <- load_province_map()

# Task 3: manuscript-motivated inward / meeting-location contact matrix (Mills
# 2026). Registered as an extra mobility matrix M11 so the existing renewal
# machinery uses a frequency-dependent, two-sided-mobility force of infection
# (a less-naive beta) with no other code change (see build_inward_contact_matrix).
# OPTIONAL family (INCLUDE_M11_MODELS, default FALSE): only built when the toggle is on,
# so the M11 renewal/Bayesian variants are absent by default (they gate on this matrix).
if (isTRUE(INCLUDE_M11_MODELS)) {
  mobility_matrices$M11 <- tryCatch(
    build_inward_contact_matrix(mobility_matrices[["M8"]], pop_vec, zones_all),
    error = function(e) { warning("inward contact matrix (M11): ", e$message); NULL })
  if (!is.null(mobility_matrices$M11))
    message("[mobility] Built M11 inward/meeting-location effective-contact matrix")
  # M11-dist: the same inward/meeting-location FOI but built on the ROAD-DISTANCE composite
  # (M8-dist), so its presence matrix inherits the km-deterrence routing. Only when the
  # road-distance family is on AND M8-dist was built.
  if (!is.null(mobility_matrices[["M8-dist"]])) {
    mobility_matrices[["M11-dist"]] <- tryCatch(
      build_inward_contact_matrix(mobility_matrices[["M8-dist"]], pop_vec, zones_all),
      error = function(e) { warning("inward contact matrix (M11-dist): ", e$message); NULL })
    if (!is.null(mobility_matrices[["M11-dist"]]))
      message("[mobility] Built M11-dist (inward FOI on the road-distance composite)")
  }
} else {
  message("[mobility] INCLUDE_M11_MODELS = FALSE — skipping the M11 inward-contact kernel(s).")
}

# Mobility-kernel comparison figure (kernel similarity, origin coverage, source-zone
# outflow profiles, and composite divergence) — reads the saved mobility_*.rds, so it
# is a cheap post-build diagnostic. tryCatch so a plotting hiccup never aborts the run.
tryCatch({
  source(file.path(ST_DIR, "27_mobility_comparison.R"))
  make_mobility_comparison_figures()
}, error = function(e) warning("[mobility] comparison figure failed: ", conditionMessage(e)))

all_weeks <- sort(unique(zone_week_nc$week_start))
t_current <- length(all_weeks)
# training_cutoff is the WEEK ANCHOR (start of the current week) that every
# weekly routine keys off (target_week_start = training_cutoff + 7*horizon, the
# forecast filename, affected_zones, the LFO fold guard). Under the analysis-date-
# anchored grid the current week ENDS on ANALYSIS_DATE, so the last calendar day
# trained on is training_cutoff + 6 (== ANALYSIS_DATE when the final week carries
# data). Report that end as the human-facing "training cutoff".
training_cutoff     <- max(all_weeks)
training_window_end <- training_cutoff + 6L
# Real-time sanity: ANALYSIS_DATE (from latest.json) should be at/after the last data
# week. If a stale/wrong processed_at lands it earlier, recent weeks would be silently
# treated as "future" and dropped by the nowcast — warn loudly rather than fail quietly.
if (ANALYSIS_DATE < training_cutoff)
  warning(sprintf("[run_all] ANALYSIS_DATE (%s) precedes the last data week (%s) — recent weeks may be dropped as future. Check latest.json 'processed_at'.",
                  ANALYSIS_DATE, training_cutoff), call. = FALSE)
gt_pmf_primary <- gt_pmfs[[GT_PRIMARY]]
message(sprintf("\n=== Fitting invasion models for cutoff: %s (t=%d) ===",
                training_cutoff, t_current))

# ── Model registry — ONLY models that work; uniform signature ──────────────
# fn(zone_week_nc, t_idx, horizons, cutoff) -> forecast tibble (health_zone,
# horizon, mu_forecast, p_invasion, method [, p_infection_invasion]).
# Primary: the mobility-informed renewal model (variants over GT / mobility /
# observation, PLUS the additional model options requested — covariate-augmented
# betas (Q2), extra mobility kernels M4b/M9/M10 (Q4), raw vs nowcast-corrected and
# reporting-rate structures (Q6), completeness-weighted calibration (Q1), and the
# NEW suspected-but-not-confirmed leading-indicator covariate models). Every
# addition is an EXTRA variant; the base formulations are untouched. OPTIONAL families
# are gated by the 00_config.R toggles: the reduced "geo" covariate model
# (INCLUDE_GEO_COV_MODELS) and the M11 inward-FOI variants (INCLUDE_M11_MODELS) are OFF
# by default; the suspected-covariate models (INCLUDE_SUSPECTED_COV_MODELS) are ON. The
# FULL-exogenous covariate model is always kept. Comparators: hhh4, gravity (B4) /
# distance (B1) baselines, stochastic SEIR (C6). Broken models (S2/S3/S4, KNN B2, null
# B3/B5/B6, INLA C4, ZINB C5) remain removed. The full grid is cross-validated; only the
# best per family is featured in figures.
# Relative reporting rate for the Q6 reporting-rate structure (health-site density
# as an ascertainment proxy: more infrastructure -> more complete confirmation).
report_rate_vec <- if ("healthsite_density" %in% names(covariates))
  setNames(covariates$healthsite_density, covariates$nom) else NULL
.wh <- function(mob, gt, obs, lbl, cov = NULL, nowcast = "corrected",
                report = NULL, betawt = "none") {
  function(zw, ti, hz, cut) {
    forecast_workhorse(zw, mobility_matrices[[mob]], gt_pmfs[[gt]], zones_all,
                       t_idx = ti, horizons = hz, obs = obs, mobility_id = mob,
                       gt_profile = gt, training_cutoff = cut, method_label = lbl,
                       covariate_spec = cov, nowcast_mode = nowcast,
                       static_cov = covariates, osrm = osrm_mat,
                       report_rate = report, beta_weighting = betawt)
  }
}
INVASION_MODELS <- list(
  # --- base mobility-informed renewal family (unchanged) --------------------
  `Renewal-M8-med`   = .wh("M8", "medium", "poisson", "Renewal-M8-med"),
  `Renewal-M8-NB`    = .wh("M8", "medium", "negbin",  "Renewal-M8-NB"),
  `Renewal-M4-med`   = .wh("M4", "medium", "poisson", "Renewal-M4-med"),
  `Renewal-M8-short` = .wh("M8", "short",  "poisson", "Renewal-M8-short"),
  `Renewal-M8-long`  = .wh("M8", "long",   "poisson", "Renewal-M8-long"),
  # --- Q4 additional mobility kernels ---------------------------------------
  `Renewal-M9-med`   = if (!is.null(mobility_matrices[["M9"]]))
    .wh("M9",  "medium", "poisson", "Renewal-M9-med") else NULL,             # ensemble kernel; OFF by default (INCLUDE_M9_MODELS)
  `Renewal-M10-med`  = .wh("M10", "medium", "poisson", "Renewal-M10-med"),  # radiation composite
  `Renewal-M4b-med`  = .wh("M4b", "medium", "poisson", "Renewal-M4b-med"),  # exp-deterrence gravity
  # --- Flowminder cohort composites (M13/M14 + geographic-distance -dist) ----
  # Cohort presence source rows (Ituri/NK/Tshopo) + gravity/radiation base, on travel-time
  # (M13/M14) or road-km geographic distance (M13/M14-dist). Guarded on the built matrix so a
  # cohort-off run (INCLUDE_COHORT_MODELS=FALSE) or missing road-distance simply omits them.
  `Renewal-M13-med`      = if (!is.null(mobility_matrices[["M13"]]))
    .wh("M13", "medium", "poisson", "Renewal-M13-med") else NULL,            # cohort + gravity (travel-time)
  `Renewal-M14-med`      = if (!is.null(mobility_matrices[["M14"]]))
    .wh("M14", "medium", "poisson", "Renewal-M14-med") else NULL,            # cohort + radiation (travel-time)
  `Renewal-M13-dist-med` = if (!is.null(mobility_matrices[["M13-dist"]]))
    .wh("M13-dist", "medium", "poisson", "Renewal-M13-dist-med") else NULL,  # cohort + gravity (road-km)
  `Renewal-M14-dist-med` = if (!is.null(mobility_matrices[["M14-dist"]]))
    .wh("M14-dist", "medium", "poisson", "Renewal-M14-dist-med") else NULL,  # cohort + radiation (road-km)
  # --- Flowminder combined inflow+outflow static family (M15/M16/M17 + M17-dist) -------------
  # M15 = combined-Flowminder static; M16 = cohort + static (the featured new composite); M17 /
  # M17-dist = grand all-kernel consensus ensemble. Guarded on the built matrix AND
  # INCLUDE_FLOWSTATIC_MODELS; M16 carries the medium/short/long GT sweep, the others medium GT
  # (the Bayesian grid carries the exhaustive covariate cross for this family).
  `Renewal-M15-med`      = if (!is.null(mobility_matrices[["M15"]]) && isTRUE(INCLUDE_FLOWSTATIC_MODELS))
    .wh("M15", "medium", "poisson", "Renewal-M15-med") else NULL,            # Flowminder static (inflow+outflow)
  `Renewal-M16-med`      = if (!is.null(mobility_matrices[["M16"]]) && isTRUE(INCLUDE_FLOWSTATIC_MODELS))
    .wh("M16", "medium", "poisson", "Renewal-M16-med") else NULL,            # cohort + static
  `Renewal-M16-short`    = if (!is.null(mobility_matrices[["M16"]]) && isTRUE(INCLUDE_FLOWSTATIC_MODELS))
    .wh("M16", "short",  "poisson", "Renewal-M16-short") else NULL,          # cohort + static (short GT)
  `Renewal-M16-long`     = if (!is.null(mobility_matrices[["M16"]]) && isTRUE(INCLUDE_FLOWSTATIC_MODELS))
    .wh("M16", "long",   "poisson", "Renewal-M16-long") else NULL,           # cohort + static (long GT)
  `Renewal-M17-med`      = if (!is.null(mobility_matrices[["M17"]]) && isTRUE(INCLUDE_FLOWSTATIC_MODELS))
    .wh("M17", "medium", "poisson", "Renewal-M17-med") else NULL,            # all-kernel consensus (travel-time)
  `Renewal-M17-dist-med` = if (!is.null(mobility_matrices[["M17-dist"]]) && isTRUE(INCLUDE_FLOWSTATIC_MODELS))
    .wh("M17-dist", "medium", "poisson", "Renewal-M17-dist-med") else NULL,  # all-kernel consensus (road-km)
  # --- Task 3: inward / meeting-location FOI (frequency-dependent, two-sided) ---
  # OPTIONAL (INCLUDE_M11_MODELS): only present when the M11 matrix was built above.
  `Renewal-M11-med`  = if (!is.null(mobility_matrices$M11))
    .wh("M11", "medium", "poisson", "Renewal-M11-med") else NULL,
  # M11 with the reduced "geo" covariate set — needs BOTH the M11 matrix and the geo family.
  `Renewal-M11-cov`  = if (!is.null(mobility_matrices$M11) && isTRUE(INCLUDE_GEO_COV_MODELS))
    .wh("M11", "medium", "poisson", "Renewal-M11-cov",
        cov = c("log_pop", "ccvi", "d_min")) else NULL,
  # --- Q2 covariate-augmented betas -----------------------------------------
  # Reduced "geo" covariate model — OPTIONAL (INCLUDE_GEO_COV_MODELS); the FULL set
  # (Renewal-M8-cov, below) is always kept.
  `Renewal-M8-geo`   = if (isTRUE(INCLUDE_GEO_COV_MODELS))
    .wh("M8", "medium", "poisson", "Renewal-M8-geo",
        cov = c("log_pop", "ccvi", "d_min")) else NULL,
  `Renewal-M8-alert` = .wh("M8", "medium", "poisson", "Renewal-M8-alert",
                           cov = c("alert_import", "alert_local")),
  # --- Suspected-but-not-confirmed leading-indicator covariates (default ON) --
  # susp_import = mobility-weighted import of OTHER zones' preceding suspected (not yet
  # confirmed) cases; susp_local = own preceding-week suspected cases. Both use the
  # nowcast-corrected suspected series and are past-only (causal in LFO), mirroring the
  # alert covariates but on the line-list suspected classification. A susp-only model and
  # a wider FULL-exogenous + susp model.
  `Renewal-M8-susp`     = if (isTRUE(INCLUDE_SUSPECTED_COV_MODELS))
    .wh("M8", "medium", "poisson", "Renewal-M8-susp",
        cov = c("susp_import", "susp_local")) else NULL,
  `Renewal-M8-cov-susp` = if (isTRUE(INCLUDE_SUSPECTED_COV_MODELS))
    .wh("M8", "medium", "poisson", "Renewal-M8-cov-susp",
        cov = c("log_pop", "ccvi", "healthsite_density", "d_min",
                "susp_import", "susp_local")) else NULL,
  # Headline-eligible covariate model uses the EXOGENOUS set only (matches the Bayesian
  # "full" grid and METHODS §5.4). `positivity` (= confirmed/tests) is REMOVED: it is
  # circular (outcome in the numerator) AND built as a full-data static per-zone mean, so
  # in an LFO fold at cutoff C it leaks confirmations from weeks > C. The alert_* signals
  # are time-varying and past-only (causal in LFO), so they stay as their own exploratory
  # model (Renewal-M8-alert); only the circular static covariate is dropped here.
  `Renewal-M8-cov`   = .wh("M8", "medium", "poisson", "Renewal-M8-cov",
                           cov = c("log_pop", "ccvi", "healthsite_density", "d_min")),
  # --- Q6 nowcast / reporting-rate structures -------------------------------
  `Renewal-M8-raw`   = .wh("M8", "medium", "poisson", "Renewal-M8-raw",
                           nowcast = "raw"),
  `Renewal-M8-report`= .wh("M8", "medium", "poisson", "Renewal-M8-report",
                           report = report_rate_vec),
  # --- Q1 completeness-weighted calibration ---------------------------------
  `Renewal-M8-cwt`   = .wh("M8", "medium", "poisson", "Renewal-M8-cwt",
                           betawt = "completeness")
  # NOTE (review §0.1): the comparators (hhh4, Gravity-B4, Distance-B1, SEIR-C6) were
  # MOVED OUT of the frequentist-gated INVASION_MODELS into the always-on BASELINE_MODELS
  # (defined below) so the model is scored against structural baselines even in a
  # Bayesian-only run. INVASION_MODELS now holds ONLY the frequentist renewal family, which
  # is retired by default (RUN_FREQUENTIST_MODELS = FALSE).
)

# Toggle: run all FREQUENTIST models (the renewal grid, comparators, and their
# ensemble)? Default FALSE — a BAYESIAN-ONLY run: no frequentist current forecast, no
# frequentist LFO members, no frequentist ensemble, and the frequentist maps / parameter
# / freq-vs-bayes products are skipped. The Bayesian suite (§ below) supplies the featured
# products; the Bayesian models are unaffected (cross-validated on the same folds regardless).
# Set TRUE here — or export RUN_FREQUENTIST_MODELS=TRUE — to also run the frequentist suite.
RUN_FREQUENTIST_MODELS <- {
  .env <- tolower(trimws(Sys.getenv("RUN_FREQUENTIST_MODELS", "")))
  if (nzchar(.env)) .env %in% c("true", "t", "1", "yes", "y")
  else isTRUE(get0("RUN_FREQUENTIST_MODELS", ifnotfound = FALSE))
}
if (!isTRUE(RUN_FREQUENTIST_MODELS)) {
  message("[models] RUN_FREQUENTIST_MODELS = FALSE — skipping ALL frequentist models (Bayesian-only run).")
  INVASION_MODELS <- list()
}

# Q3 ensemble members: a PRE-SPECIFIED, diverse set of individually strong
# mobility-informed renewal formulations (composite-gravity, plain gravity,
# multi-kernel-ensemble and radiation mobilities; Poisson and negative-binomial
# observation). Pre-specification means the ensemble does no selection on the
# test folds and is scored in exactly the same LFO as its members.
# Drop any variants whose dependencies were unavailable (e.g. M11 build failed).
INVASION_MODELS <- Filter(Negate(is.null), INVASION_MODELS)
# `.have_freq` gates every frequentist-specific PRODUCT below (current forecast, risk
# maps, parameter screen, freq-vs-bayes comparison); the Bayesian branch is independent.
.have_freq <- length(INVASION_MODELS) > 0L

# ── Always-on BASELINE comparators (review §0.1 / §3.1 / §3.2) ───────────────
# Decoupled from RUN_FREQUENTIST_MODELS so the featured model is ALWAYS scored against
# structural baselines — the reviewer's central evaluation ask — even in the Bayesian-only
# run: gravity-flow (B4), inverse-distance (B1), the NEAREST-AFFECTED / adjacency
# spatial-spread null (B7, new — "does the virus just go next door?"), and the SEIR / hhh4
# mechanistic comparators. These are yardsticks to beat, NOT the removed frequentist renewal
# family. They are scored on the identical LFO folds as the Bayesian models (below).
BASELINE_MODELS <- Filter(Negate(is.null), list(
  `Gravity-B4` = function(zw, ti, hz, cut) {
    Yw <- zone_week_to_wide(zw, zones_all)
    forecast_B4(Yw, pop_vec, osrm_mat, ti, hz, zones_all) %>% dplyr::mutate(method = "Gravity-B4")
  },
  `Distance-B1` = function(zw, ti, hz, cut) {
    Yw <- zone_week_to_wide(zw, zones_all)
    forecast_B1(Yw, osrm_mat, ti, hz, zones_all, alpha = 1) %>% dplyr::mutate(method = "Distance-B1")
  },
  `Adjacency-B7` = function(zw, ti, hz, cut) {
    Yw <- zone_week_to_wide(zw, zones_all)
    forecast_B7_adjacency(Yw, osrm_mat, ti, hz, zones_all, method = "distance") %>%
      dplyr::mutate(method = "Adjacency-B7")
  },
  # hhh4 and SEIR are the EXPENSIVE per-fold mechanistic baselines — gated OFF by default
  # (INCLUDE_HHH4_BASELINE / INCLUDE_SEIR_BASELINE, 00_config.R) to save runtime. NULL entries
  # are dropped by Filter(Negate(is.null)); the cheap deterministic baselines above always run.
  `hhh4` = if (isTRUE(get0("INCLUDE_HHH4_BASELINE", ifnotfound = FALSE))) function(zw, ti, hz, cut) {
    r <- tryCatch(run_all_hhh4(zone_week_nc = zw, W_outflow = W_primary,
                    gt_pmf_daily = gt_pmf_primary, pop_vec = pop_vec,
                    covariates = covariates, zones_all = zones_all, t_idx = ti,
                    horizons = hz, n_sim = 200L, training_cutoff = cut, cache = FALSE),
                  error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)
    variants <- unique(r$method); pick <- variants[grepl("C1c", variants)][1]
    if (is.na(pick)) pick <- variants[1]
    r %>% dplyr::filter(method == pick) %>% dplyr::mutate(method = "hhh4")
  } else NULL,
  `SEIR-C6` = if (isTRUE(get0("INCLUDE_SEIR_BASELINE", ifnotfound = FALSE))) function(zw, ti, hz, cut) {
    tryCatch(forecast_seir(zone_week_nc = zw, W = W_primary, pop_vec = pop_vec,
                    Rt_national = 1.5, zones_all = zones_all, t_idx = ti,
                    horizons = hz, n_sim = 200L, training_cutoff = cut) %>%
             dplyr::mutate(method = "SEIR-C6"), error = function(e) NULL)
  } else NULL
))
message(sprintf("[baselines] %d always-on baseline comparators: %s",
                length(BASELINE_MODELS), paste(names(BASELINE_MODELS), collapse = ", ")))

ENSEMBLE_MEMBERS <- intersect(
  c("Renewal-M8-med", "Renewal-M8-NB", "Renewal-M4-med", "Renewal-M9-med",
    "Renewal-M10-med", "Renewal-M11-med"),
  names(INVASION_MODELS))

# ── Current-week forecasts (all models) ─────────────────────────────────────
affected_now <- affected_zones(zone_week_nc, training_cutoff)
message(sprintf("[invasion] %d/%d zones already affected; %d at-risk",
                length(affected_now), length(zones_all),
                length(zones_all) - length(affected_now)))
fc_all_current <- purrr::imap_dfr(INVASION_MODELS, function(fn, nm) {
  fc <- tryCatch(fn(zone_week_nc, t_current, LFO_HORIZONS, training_cutoff),
                 error = function(e) { warning(nm, " failed: ", e$message); NULL })
  if (is.null(fc)) return(NULL)
  fc$method <- nm
  fc
})
# Enforce the invasion framing on EVERY model: affected zones carry no invasion probability,
# anywhere. GUARDED: with no frequentist models (Bayesian-only run) imap_dfr returns an empty
# 0-column tibble, so this mutate — which references health_zone — must be skipped rather than
# error with "object 'health_zone' not found". The Bayesian forecasts already mask affected
# zones inside predict_bayes_invasion().
if (nrow(fc_all_current) > 0 && "health_zone" %in% names(fc_all_current)) {
  fc_all_current <- fc_all_current %>%
    dplyr::mutate(
      was_active_before = health_zone %in% affected_now,
      p_invasion   = ifelse(was_active_before, NA_real_, p_invasion),
      mu_forecast  = ifelse(was_active_before, NA_real_, mu_forecast),
      p_infection_invasion = if ("p_infection_invasion" %in% names(.))
        ifelse(was_active_before, NA_real_, p_infection_invasion) else NA_real_
    )
}
# Q3: add mean / median ensembles of the pre-specified renewal members. Affected
# zones have NA for every member, so the combined value is NaN there — re-mask to
# NA to preserve the invasion framing (no probability for already-affected zones).
# GUARDED: with < 2 members present (e.g. a Bayesian-only run, where ENSEMBLE_MEMBERS
# is empty) there is no frequentist ensemble to build, so leave fc_all_current as-is
# rather than invoking the ensemble machinery on nothing.
if (length(ENSEMBLE_MEMBERS) >= 2) {
  fc_all_current <- append_ensembles(fc_all_current, ENSEMBLE_MEMBERS) %>%
    dplyr::mutate(dplyr::across(
      dplyr::any_of(c("p_invasion", "p_case_invasion", "p_infection_invasion",
                      "mu_forecast", "mu_wk0")),
      ~ ifelse(is.nan(.x), NA_real_, .x)))
}
message(sprintf("[models] %d forecast rows from %d models (incl. ensembles)",
                nrow(fc_all_current),
                if (nrow(fc_all_current) > 0 && "method" %in% names(fc_all_current))
                  dplyr::n_distinct(fc_all_current$method) else 0L))
saveRDS(fc_all_current,
        file.path(OUT_FORECASTS, paste0("current_forecasts_", training_cutoff, ".rds")))

# ── Daily re-issue is DEFERRED to after the Bayesian suite (§ below) so it anchors the
# FREQUENTIST current forecast AND the featured Bayesian model + stacked ensemble in a single
# upsert (the persistence keys on forecast_date, so one combined call — not two — is correct).

# ---------------------------------------------------------------------------
# 10. Invasion LFO-CV — same principled folds for ALL models
# ---------------------------------------------------------------------------
.phase("Step 4-9 nowcast (epinowcast) + current fc")
message("\n=== Step 10: Invasion LFO-CV ===")
# Use the week floor of OUTBREAK_START (a Thursday) under the analysis-date-anchored
# grid, so the first outbreak week is kept rather than silently dropped — it yields an
# extra early evaluable fold.
zone_week_outbreak <- dat$zone_week %>%
  dplyr::filter(week_start >= lubridate::floor_date(OUTBREAK_START, "week",
                                                    week_start = get0("WEEK_ANCHOR", ifnotfound = 1L)))

# Task 6: evaluate a curated set of BAYESIAN renewal models in the SAME folds as
# the frequentist models (each refit on the fold's training data, cmdstanr reusing
# the compiled Stan binary), so the over-time evaluation includes them. Kept small
# (2 specs, lighter sampling) to bound cross-validation runtime; toggle with
# RUN_BAYES_LFO.
RUN_BAYES_LFO <- TRUE
# EVERY Bayesian model in the current-forecast grid is also cross-validated on the
# SAME leave-future-out folds as the frequentist ones (single source of truth:
# bayes_default_grid). So the ranking / over-time / detection plots compare exactly the
# Bayesian grid that bayes_default_grid() composes from the CORE (mobility {M4,M8,M9,M10}
# x GT {short,medium,long}, plus the FULL-exogenous covariate and the suspected-covariate
# models) and whichever OPTIONAL families are toggled on in 00_config.R (M11, OSRM -dist
# kernels, the reduced "geo" covariate set, and the logit-link sensitivity) — and the
# FEATURED Bayesian model is chosen BY CV SKILL (best_bayes_method), not hardcoded. Each model is
# refit per fold by MCMC, so LFO uses lighter sampling (iter=600) to keep the full-
# grid cross-validation tractable; set RUN_BAYES_LFO <- FALSE to skip it entirely.
BAYES_LFO_SPECS <- bayes_default_grid(mobility_matrices)
bayes_lfo_models <- list()
if (isTRUE(RUN_BAYES_LFO) && requireNamespace("brms", quietly = TRUE)) {
  for (sp in BAYES_LFO_SPECS) {
    bayes_lfo_models[[sp$label]] <- make_bayes_lfo_model(sp$mob, sp$gt, sp$cov,
      mobility_matrices, gt_pmfs, covariates, osrm_mat, zones_all,
      iter = 600L, link = sp$link %||% "cloglog")
  }
}
message(sprintf("[bayes] cross-validating %d Bayesian models in LFO", length(bayes_lfo_models)))
# Leakage-free per-fold training reconstruction: when TRUE the weekly LFO rebuilds
# each fold's training counts from the line list censored to the forecast moment
# (cut+7), excluding cases only reported later — rather than slicing the final
# onset-bucketed counts (which carry a mild training-side revision leak). This
# CHANGES the reported LFO metrics relative to the previous behaviour; set to FALSE
# to reproduce the older (mildly leaky) numbers.
LEAKAGE_FREE_LFO <- TRUE
# REUSE_CACHED_LFO (default FALSE): load the previous run's lfo_cv_results.rds instead of
# re-running the hours-long cross-validation. Valid ONLY when nothing that affects the LFO has
# changed since that cache was written — used for the GT-marginalisation deployed-forecast
# re-run, which changes ONLY the current-forecast step (predict_bayes_gt_marginal), not the LFO.
# The cached RDS already has the ensembles appended, so the append/save block below is skipped.
.lfo_cache <- file.path(OUT_FORECASTS, "lfo_cv_results.rds")
.reuse_lfo <- (tolower(trimws(Sys.getenv("REUSE_CACHED_LFO", ""))) %in% c("true", "t", "1", "yes", "y") ||
               isTRUE(get0("REUSE_CACHED_LFO", ifnotfound = FALSE))) && file.exists(.lfo_cache)
lfo_results <- if (isTRUE(.reuse_lfo)) {
  message("[LFO] REUSE_CACHED_LFO=TRUE — loading cached lfo_cv_results.rds (cross-validation NOT re-run)")
  tryCatch(readRDS(.lfo_cache), error = function(e) { warning("cached LFO load failed: ", e$message); NULL })
} else tryCatch(
  # BASELINE_MODELS are always scored (review §3.1/§3.2), alongside the Bayesian models and
  # (only if RUN_FREQUENTIST_MODELS=TRUE) the frequentist renewal family.
  run_invasion_lfo(zone_week_outbreak, c(BASELINE_MODELS, INVASION_MODELS, bayes_lfo_models),
                   horizons = LFO_HORIZONS, analysis_date = ANALYSIS_DATE,
                   nowcast_fn = apply_nowcast_correction,
                   linelist = if (isTRUE(LEAKAGE_FREE_LFO)) dat$ll else NULL),
  error = function(e) { warning("Invasion LFO failed: ", e$message); NULL }
)
if (isTRUE(.reuse_lfo) && !is.null(lfo_results))
  message(sprintf("[LFO] cached: %d rows; %d methods (ensembles already appended)",
                  nrow(lfo_results), dplyr::n_distinct(lfo_results$method)))
if (!is.null(lfo_results) && !isTRUE(.reuse_lfo)) {
  # Q3: build the mean / median ensembles inside the SAME folds, so they are
  # evaluated identically to their members (each fold's at-risk outcome is
  # carried through the combination).
  lfo_results <- append_ensembles(lfo_results, ENSEMBLE_MEMBERS)
  # Bayesian ENSEMBLE in the SAME folds — the pre-specified, leakage-free analogue of the
  # frequentist ensemble: mean/median over the diverse mobility-kernel Bayesian members
  # (M4/M8/M9/M10/M11 at medium GT, cloglog). Labelled "Bayes-ens-{mean,median}" so it is
  # picked up by the Bayesian-restricted discrimination / detection / top-K figures and can
  # be compared against the individual Bayesian models. (The loo-STACKED ensemble remains the
  # separate full-data featured product — "Bayes-stack".)
  BAYES_ENSEMBLE_MEMBERS <- intersect(
    c("Bayes-M4-med", "Bayes-M8-med", "Bayes-M9-med", "Bayes-M10-med", "Bayes-M11-inward"),
    unique(lfo_results$method))
  if (length(BAYES_ENSEMBLE_MEMBERS) >= 2L)
    lfo_results <- append_ensembles(lfo_results, BAYES_ENSEMBLE_MEMBERS, prefix = "Bayes-ens")
  saveRDS(lfo_results, file.path(OUT_FORECASTS, "lfo_cv_results.rds"))
  message(sprintf("[LFO] %d pooled at-risk rows; %d invasion events (h1); methods: %d",
                  nrow(lfo_results),
                  sum(lfo_results$is_new_invasion[lfo_results$horizon == 1], na.rm = TRUE),
                  dplyr::n_distinct(lfo_results$method)))
}

# ── Fail loudly rather than silently emit stale products ─────────────────────
# run_invasion_lfo() returning NULL means the leave-future-out cross-validation
# ERRORED (the tryCatch above turns any failure into NULL). Everything downstream
# — the featured-model selection (best_bayes_method), risk scores, the prospective
# and held-out checks, and the figures — keys on lfo_results. If it is NULL the run
# would otherwise "complete" while leaving the PREVIOUS run's files untouched, so
# the reported featured model and the on-disk data silently disagree (exactly the
# failure a prior run hit via future.globals.maxSize). Halt instead.
if (is.null(lfo_results)) {
  stop("[LFO] run_invasion_lfo() returned NULL: the leave-future-out cross-validation FAILED ",
       "(see the 'Invasion LFO failed' warning above). Refusing to continue — every downstream ",
       "product (model selection, risk scores, prospective/held-out checks, figures) would be ",
       "left as STALE files from the previous run while the labels claim the new featured model. ",
       "Fix the LFO failure and re-run. Common cause: future.globals.maxSize too small for the ",
       "fitted-model list (raised to 12 GiB near the top of this script).")
}

# ── Daily-issue backtest: validate the LIVE operating point (partial issue-week
# in training, mid-week issue) that the weekly LFO above does NOT exercise. Uses
# the deterministic per-fold nowcast (epinowcast refit per issue-fold is
# intractable; the live forecast still refits it once per day). Bounded to a few
# issue offsets per week to keep runtime sane; set RUN_DAILY_LFO <- FALSE to skip.
# Skipped in a Bayesian-only run: it backtests the FREQUENTIST INVASION_MODELS (the daily
# re-issue anchors the frequentist current forecast, which is absent when frequentist is off).
RUN_DAILY_LFO <- TRUE
daily_lfo_results <- NULL
if (isTRUE(RUN_DAILY_LFO) && .have_freq) {
  # Ensembles are built INSIDE run_invasion_lfo_daily on the weekly member forecasts
  # BEFORE the analysis-date anchoring (so the scored ensemble matches the live one);
  # do NOT append_ensembles again on the already-anchored output.
  daily_lfo_results <- tryCatch(
    run_invasion_lfo_daily(zone_week_outbreak, dat$ll, INVASION_MODELS, zones_all,
                           windows_days = c(7L, 14L), issue_offsets = c(0L, 3L, 6L),
                           ensemble_members = ENSEMBLE_MEMBERS,
                           analysis_date = ANALYSIS_DATE,
                           nowcast_fn = apply_nowcast_correction),
    error = function(e) { warning("Daily-issue LFO failed: ", e$message); NULL })
  if (!is.null(daily_lfo_results)) {
    saveRDS(daily_lfo_results, file.path(OUT_FORECASTS, "lfo_daily_results.rds"))
    message(sprintf("[LFO-daily] %d pooled rows; %d issue-folds; methods: %d",
                    nrow(daily_lfo_results),
                    dplyr::n_distinct(daily_lfo_results$fold_id),
                    dplyr::n_distinct(daily_lfo_results$method)))
  }
}

# ---------------------------------------------------------------------------
# 11. Invasion evaluation (pooled; AUC-PR, ranking, log-score, calibration)
# ---------------------------------------------------------------------------
.phase("Step 10  LFO-CV")
message("\n=== Step 11: Invasion evaluation ===")
eval_tbl <- if (!is.null(lfo_results)) tryCatch(evaluate_invasion(lfo_results),
              error = function(e) { warning("eval failed: ", e$message); NULL }) else NULL
if (!is.null(eval_tbl)) {
  readr::write_csv(eval_tbl, file.path(OUT_DIAGNOSTICS, "invasion_evaluation.csv"))
  message("[eval] Models by AUC-PR skill (h=1):")
  print(as.data.frame(eval_tbl %>% dplyr::filter(horizon == 1) %>%
    dplyr::transmute(method, n_inv = n_invasions, auc_pr = round(auc_pr, 3),
      aucpr_skill = round(auc_pr_skill, 1), rank_of_truth = round(mean_rank_of_truth, 1),
      log_score = round(log_score, 4), cal = round(calibration_in_large, 1))))
}
# metrics_summary alias for downstream figures/report
metrics_summary <- eval_tbl

# ── Held-out (last-two-origins) confirmatory evaluation + optimism gap (§2.2) ──
# The all-fold eval_tbl above SELECTS the featured model AND reports its skill, so it is
# selection-optimistic. Re-evaluate honestly: select on the earlier folds, score once on the
# two most-recent origins (never seen by selection). Report both and the gap.
heldout_eval <- if (!is.null(lfo_results))
  tryCatch(evaluate_invasion_heldout(lfo_results, n_holdout = 2L, restrict = "^Bayes", horizon = 1L),
           error = function(e) { warning("held-out eval failed: ", e$message); NULL }) else NULL
if (!is.null(heldout_eval)) {
  readr::write_csv(
    tibble::tibble(selected = heldout_eval$selected, horizon = heldout_eval$horizon,
                   inner_cv_skill = heldout_eval$inner_skill,
                   heldout_skill  = heldout_eval$outer_skill,
                   optimism_gap   = heldout_eval$optimism_gap,
                   n_holdout_origins = heldout_eval$n_holdout),
    file.path(OUT_DIAGNOSTICS, "invasion_heldout_optimism.csv"))
  message(sprintf("[held-out §2.2] %s: inner-CV AUC-PR skill %.1fx vs held-out (last 2 origins) %.1fx (optimism gap %.1fx)",
                  heldout_eval$selected, heldout_eval$inner_skill %||% NA_real_,
                  heldout_eval$outer_skill %||% NA_real_, heldout_eval$optimism_gap %||% NA_real_))
}

# NOTE: the forecast-vs-observed + prospective block (§3.4/§3.5/§3.6/§3.8) is placed
# AFTER best_bayes_method is defined (below the headline-model selection), since it keys
# on the featured Bayesian method.

# Headline model (best overall) + best renewal model, both by a
# calibration-aware criterion (discrimination + targeting + log-score) POOLED across
# both forecast horizons, so we do not feature an over-confident variant or one that
# is strong at 1-week but weak at 2-week (or vice versa).
# Degenerate fallback (only if the whole evaluation is missing/empty): the first RENEWAL model
# in the grid — derived from INVASION_MODELS, not a hard-coded name, so it always tracks the grid.
# NA_character_ (a length-1 NA, NOT names(list())[1] which is NULL) when there is no
# frequentist model, so the `if (is.na(primary_method))` guards below never see length-0.
.renewal_default <- { .rn <- grep("^Renewal", names(INVASION_MODELS), value = TRUE)
                      if (length(.rn)) .rn[1] else NA_character_ }
best_method <- if (!is.null(eval_tbl)) best_invasion_model(eval_tbl) else .renewal_default
primary_method <- if (!is.null(eval_tbl))
  best_invasion_model(eval_tbl, restrict = "^Renewal") else .renewal_default
if (is.na(primary_method)) primary_method <- .renewal_default
if (is.na(best_method))    best_method    <- primary_method
# Featured BAYESIAN model: the best cross-validated SINGLE Bayesian model by the calibration-aware
# CV composite (summed within-horizon ranks of AUC-PR skill + mean rank-of-truth + log-score, POOLED
# across both horizons) — the same leave-future-out forecast-skill criterion as the frequentist
# featured model. This is the PRIMARY (and default) selector; the loo predictive-stacking weight is
# NOT used to pick the featured single model (it defines the loo-stacked ENSEMBLE, bayes_ensemble_*).
# Not hardcoded. Ensembles (Bayes-ens-*) are excluded so the featured single model always has a
# current-forecast row.
best_bayes_method <- if (!is.null(eval_tbl) && any(grepl("^Bayes", eval_tbl$method))) {
  .bayes_singles <- eval_tbl %>% dplyr::filter(grepl("^Bayes", method), !grepl("-ens-", method))
  if (nrow(.bayes_singles) > 0) best_invasion_model(.bayes_singles) else NA_character_
} else NA_character_
if (length(best_bayes_method) != 1 || is.na(best_bayes_method)) {
  .cand <- if (!is.null(lfo_results)) grep("^Bayes", unique(lfo_results$method), value = TRUE) else character(0)
  .cand <- .cand[!grepl("-ens-", .cand)]
  best_bayes_method <- if (length(.cand)) .cand[1] else NULL
}
# Bayesian-only / degenerate fallback: if no frequentist model produced a featured pick
# (RUN_FREQUENTIST_MODELS = FALSE, or the whole frequentist grid failed), the headline model
# is the best Bayesian one — best_method drives the skill/lead/detection/predobs figures, so it
# must be a real cross-validated method rather than NA.
if ((length(best_method) != 1 || is.na(best_method)) && !is.null(best_bayes_method))
  best_method <- best_bayes_method
message(sprintf("[eval] Headline model: %s | best renewal model: %s",
                best_method,
                if (length(primary_method) != 1 || is.na(primary_method)) "none (Bayesian-only)" else primary_method))

# ── Forecast-vs-observed suite + prospective check (review §3.4/§3.5/§3.6/§3.8) ──
# Uses the featured Bayesian model's LFO forecasts (predicted vs realised). Guarded and
# additive; produces CSVs that feed the reliability / count-calibration / prospective figures.
# Placed here (not earlier) because it keys on best_bayes_method, defined just above.
if (!is.null(lfo_results) && !is.null(best_bayes_method) &&
    best_bayes_method %in% lfo_results$method) tryCatch({
  fvo <- lfo_results %>% dplyr::filter(method == best_bayes_method, horizon == 1L,
                                       !was_active_before, is.finite(p_invasion))
  if (nrow(fvo)) {
    # (§3.5 view 1) reliability curve; (view 2) aggregate-count calibration over folds.
    rc <- reliability_curve(fvo$p_invasion, as.integer(fvo$is_new_invasion), n_bins = 10L)
    if (nrow(rc)) readr::write_csv(rc, file.path(OUT_DIAGNOSTICS, "forecast_reliability_h1.csv"))
    acc <- aggregate_count_calibration(data.frame(time = as.character(fvo$fold_id),
             p_invasion = fvo$p_invasion, invaded = as.integer(fvo$is_new_invasion)))
    if (nrow(acc)) readr::write_csv(acc, file.path(OUT_DIAGNOSTICS, "forecast_count_calibration_h1.csv"))
  }
  # (§3.4/§3.8) prospective check: zones invaded AFTER the last CV origin — were they ranked
  # highly by the forecast issued at that origin (before their first case)?
  last_cut <- max(lfo_results$cutoff, na.rm = TRUE)
  fc_last  <- lfo_results %>% dplyr::filter(method == best_bayes_method, horizon == 1L,
                                            cutoff == last_cut, !was_active_before)
  aff_at_cut    <- affected_zones(zone_week_nc, last_cut)
  invaded_after <- setdiff(affected_now, aff_at_cut)
  if (nrow(fc_last) && length(invaded_after)) {
    pc <- prospective_invasion_check(
      dplyr::transmute(fc_last, health_zone, p_invasion), invaded_after, top_k = 15L)
    readr::write_csv(pc$summary,  file.path(OUT_DIAGNOSTICS, "prospective_invasion_summary.csv"))
    readr::write_csv(pc$per_zone, file.path(OUT_DIAGNOSTICS, "prospective_invasion_per_zone.csv"))
    message(sprintf("[prospective §3.4] %d zones invaded after the last CV origin; median pre-invasion rank %.0f; %d/%d in top-15",
                    pc$summary$n_zones, pc$summary$median_rank %||% NA_real_,
                    pc$summary$n_in_topk, pc$summary$n_zones))
  }
}, error = function(e) warning("forecast-vs-observed / prospective suite: ", e$message))

# Render the review-response figures (reliability, count-calibration, prospective ranks,
# held-out optimism) in the key_outputs house style, from the CSVs just written (§3.4/3.5/2.2).
tryCatch({
  .cal <- if (!is.null(eval_tbl) && !is.null(best_bayes_method))
    eval_tbl$calibration_in_large[eval_tbl$method == best_bayes_method & eval_tbl$horizon == 1L][1] else NA_real_
  make_review_figures(OUT_DIAGNOSTICS,
                      file.path(OUT_DIR, "key_outputs", "manuscript_figures", "panels"),
                      cal_in_large = if (is.finite(.cal)) .cal else NULL, top_k = 15L)
}, error = function(e) warning("review figures: ", e$message))

# Curated display set for the (otherwise cluttered) multi-model figures: the
# best-fitting model + best renewal model, both ensembles, the best few renewal variants,
# and the comparators. The FULL grid remains in invasion_evaluation.csv.
eval_display <- eval_tbl
if (!is.null(eval_tbl)) {
  top_renew <- eval_tbl %>% dplyr::filter(horizon == 1, grepl("^Renewal", method)) %>%
    dplyr::arrange(dplyr::desc(auc_pr_skill)) %>% dplyr::pull(method) %>% head(4)
  # Bayesian-only run (no frequentist models): fold the top Bayesian models + ensemble into the
  # display so the MAIN discrimination/summary figures are not degenerate. Left empty when
  # frequentist models ran, so the default frequentist display is unchanged.
  top_bayes <- if (!.have_freq) {
    .tb <- eval_tbl %>% dplyr::filter(horizon == 1, grepl("^Bayes", method)) %>%
      dplyr::arrange(dplyr::desc(auc_pr_skill)) %>% dplyr::pull(method)
    unique(c(head(.tb, 6L), "Bayes-ens-mean", "Bayes-ens-median"))
  } else character(0)
  DISPLAY_MODELS <- unique(c(best_method, primary_method, "Ensemble-mean",
    "Ensemble-median", top_renew, top_bayes, "hhh4", "Gravity-B4", "Distance-B1", "SEIR-C6"))
  eval_display <- eval_tbl %>% dplyr::filter(method %in% DISPLAY_MODELS)
}

# ---------------------------------------------------------------------------
# 11b. Spatiotemporal evaluation (Q5): skill-over-time, lead-time, spatial error
# ---------------------------------------------------------------------------
st_skill <- NULL; lead_tbl <- NULL; zone_err <- NULL
if (!is.null(lfo_results)) {
  message("\n=== Step 11b: Spatiotemporal evaluation (Q5) ===")
  first_case_wk <- first_case_from_zone_week(zone_week_outbreak)
  st_skill <- tryCatch(spatiotemporal_skill(lfo_results, k = 5L),
                       error = function(e) { warning("skill-over-time: ", e$message); NULL })
  lead_tbl <- tryCatch(lead_time_analysis(lfo_results, first_case_wk, k = 5L,
                         horizon = 1L, method = best_method),
                       error = function(e) { warning("lead-time: ", e$message); NULL })
  zone_err <- tryCatch(zone_spatial_error(lfo_results, province_map, horizon = 1L,
                         method = best_method),
                       error = function(e) { warning("spatial-error: ", e$message); NULL })
  if (!is.null(st_skill))
    readr::write_csv(st_skill, file.path(OUT_DIAGNOSTICS, "skill_over_time.csv"))
  if (!is.null(lead_tbl))
    readr::write_csv(lead_tbl, file.path(OUT_DIAGNOSTICS, "lead_time.csv"))
  if (!is.null(zone_err))
    readr::write_csv(zone_err, file.path(OUT_DIAGNOSTICS, "zone_spatial_error.csv"))
  if (!is.null(lead_tbl))
    message(sprintf("[Q5] %d newly affected zones analysed; %d flagged in top-5 pre-invasion",
                    nrow(lead_tbl), sum(lead_tbl$ever_topk, na.rm = TRUE)))
}

# ---------------------------------------------------------------------------
# 12. Risk scores on the best renewal model (at-risk zones only)
# ---------------------------------------------------------------------------
message("\n=== Step 12: Risk scores ===")
# 20_forecast_detail.R + 21_bayesian_renewal.R already sourced (main source block)
# Frequentist featured-model risk scores. In a Bayesian-only run (no frequentist models)
# risk_scores is NULL, so every `if (!is.null(risk_scores))` frequentist product below is
# skipped and the Bayesian suite supplies the featured maps/tables instead.
primary_fc  <- if (.have_freq) fc_all_current %>% dplyr::filter(method == primary_method) else fc_all_current
risk_scores <- if (!.have_freq) NULL else
  tryCatch(compute_risk_scores(primary_fc, province_map),
           error = function(e) { warning("[run_all] compute_risk_scores failed (risk tables/maps degrade to NA): ",
                                         conditionMessage(e), call. = FALSE); primary_fc })
# Q2 + Q3: attach 0-1 relative-risk indices (1 = highest) and the vulnerability-
# and-capacity-adjusted preparedness-priority score (invasion risk x vulnerability).
# vuln_index is computed regardless of the frequentist toggle — the Bayesian suite uses it too.
vuln_index <- tryCatch(compute_vulnerability_index(covariates, zones_all, osrm_mat = osrm_mat),
                       error = function(e) { warning("vulnerability index: ", e$message); NULL })
if (!is.null(vuln_index) && !is.null(risk_scores))
  risk_scores <- tryCatch(add_risk_indices(risk_scores, vuln_index),
                          error = function(e) { warning("risk indices: ", e$message); risk_scores })
# Attach the ENSEMBLE member-spread (min/max p across members) as p_lo/p_hi so the
# frequentist featured model gets the same probability+uncertainty map (#2) as the
# Bayesian (whose p_lo/p_hi are the posterior 90% CrI).
if (!is.null(risk_scores)) {
  .unc <- tryCatch(ensemble_member_uncertainty(fc_all_current, ENSEMBLE_MEMBERS, "p_invasion"),
                   error = function(e) NULL)
  if (!is.null(.unc) && all(c("health_zone", "horizon", "p_lo", "p_hi") %in% names(.unc)) &&
      !any(c("p_lo", "p_hi") %in% names(risk_scores)))
    risk_scores <- dplyr::left_join(risk_scores,
      dplyr::select(.unc, health_zone, horizon, p_lo, p_hi), by = c("health_zone", "horizon"))
}
saveRDS(risk_scores, file.path(OUT_FORECASTS, "risk_scores_current.rds"))

# Full per-zone invasion probabilities + risk scores for ALL health zones (every
# zone x both horizons; already-affected zones carry NA by construction), exported
# to CSV alongside the RDS so the complete table is usable outside R.
.write_risk_csv <- function(rs, path) {
  if (is.null(rs) || !nrow(rs)) return(invisible(NULL))
  # p_lo/p_hi are the 5% and 95% bounds on p_case_invasion (the 90% CrI for the Bayesian scores;
  # the ensemble min-max spread for the frequentist/ensemble scores) — exported as p_case_lo /
  # p_case_hi. p_infection_invasion is intentionally NOT exported.
  keep <- c("health_zone", "province", "horizon", "was_active_before",
            "p_case_invasion", "p_lo", "p_hi", "p_median",
            # Posterior RANK credible interval (review §5.1). NOTE: for the featured
            # intercept-only model the ranking is DETERMINISTic given the fixed import
            # forces (beta rescales all zones together, preserving order), so these are a
            # point [r,r]; they carry genuine width only for covariate-modulated /
            # GT-marginalised / ensemble variants, where the posterior ranking varies.
            "rank_med", "rank_lo", "rank_hi",
            "mu_forecast", "rr_nat", "rr_nat_rank", "rr01_nat",
            "rr_ituri", "rr_ituri_rank", "rr_nordkivu", "rr_nordkivu_rank",
            "rr_hautuele", "rr_hautuele_rank", "V", "priority", "priority_rank",
            "surveillance_gap", "healthcare_gap", "access_gap", "social_vulnerability",
            "healthcare_travel_min", "method")
  out <- rs %>% dplyr::select(dplyr::any_of(keep)) %>%
    dplyr::rename(dplyr::any_of(c(p_case_lo = "p_lo", p_case_hi = "p_hi"))) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 4)))
  if (all(c("horizon", "p_case_invasion") %in% names(out)))
    out <- out %>% dplyr::arrange(horizon, dplyr::desc(dplyr::coalesce(p_case_invasion, -1)))
  readr::write_csv(out, path)
  message(sprintf("[export] %d zone-rows (%d zones x horizons) -> %s",
                  nrow(out), dplyr::n_distinct(out$health_zone), basename(path)))
  invisible(out)
}
tryCatch(.write_risk_csv(risk_scores, file.path(OUT_REPORTS, "risk_scores_all_zones.csv")),
         error = function(e) warning("risk_scores CSV export: ", e$message))

# Harmonised per-zone cumulative confirmed cases (line-list ∪ sitrep, whatever
# zone_week_nc carries) at the training cutoff — the dashboard reads this to
# colour active zones and place case markers. Computed from the SAME
# zone_week_nc / training_cutoff / `confirmed` column as affected_zones(), so the
# invariant holds by construction: cumulative_confirmed_cases > 0 iff the zone is
# in affected_now (was_active_before).
harmonised_confirmed <- zone_week_nc %>%
  dplyr::filter(week_start <= training_cutoff) %>%
  dplyr::group_by(health_zone) %>%
  dplyr::summarise(cumulative_confirmed_cases = sum(confirmed, na.rm = TRUE),
                   .groups = "drop")
readr::write_csv(harmonised_confirmed,
                 file.path(OUT_REPORTS, "harmonised_confirmed_cases.csv"))

if (!is.null(risk_scores) && "p_case_invasion" %in% names(risk_scores)) {
  rs1 <- risk_scores %>% dplyr::filter(horizon == 1, !was_active_before)
  has_pri <- all(c("rr01_nat", "priority", "priority_rank") %in% names(rs1))
  # Relative-risk columns exist only when compute_risk_scores succeeds (its error
  # fallback returns the bare forecast with p_case_invasion but no rr_* columns);
  # guard them so that fallback does not abort the whole run here.
  has_rr <- all(c("rr_ituri", "rr_nat", "rr_nat_rank") %in% names(rs1))
  ituri_top <- rs1 %>% dplyr::filter(province == "Ituri") %>%
    dplyr::arrange(dplyr::desc(p_case_invasion)) %>%
    dplyr::transmute(health_zone, province,
      p_case = round(p_case_invasion, 3), p_infection = round(p_infection_invasion, 3),
      rr01 = if (has_pri) round(rr01_nat, 3) else NA_real_,
      rr_ituri = if (has_rr) round(rr_ituri, 1) else NA_real_,
      rr_nat = if (has_rr) round(rr_nat, 1) else NA_real_,
      vulnerability = if (has_pri) round(V, 2) else NA_real_,
      priority = if (has_pri) round(priority, 3) else NA_real_) %>% head(15)
  nat_top <- rs1 %>% dplyr::arrange(dplyr::desc(p_case_invasion)) %>%
    dplyr::transmute(health_zone, province,
      p_case = round(p_case_invasion, 3), p_infection = round(p_infection_invasion, 3),
      rr01 = if (has_pri) round(rr01_nat, 3) else NA_real_,
      rr_nat = if (has_rr) round(rr_nat, 1) else NA_real_,
      rr_nat_rank = if (has_rr) rr_nat_rank else NA_integer_) %>% head(15)
  readr::write_csv(ituri_top, file.path(OUT_REPORTS, "risk_table_ituri.csv"))
  readr::write_csv(nat_top,   file.path(OUT_REPORTS, "risk_table_national.csv"))
  # Q3: standalone preparedness-priority table (top zones by the composite).
  if (has_pri) {
    priority_top <- rs1 %>% dplyr::arrange(dplyr::desc(priority)) %>%
      dplyr::transmute(health_zone, province, priority = round(priority, 3),
        rr01_invasion = round(rr01_nat, 3), vulnerability = round(V, 2),
        surveillance_gap = round(surveillance_gap, 2),
        healthcare_gap = round(healthcare_gap, 2),
        social_vulnerability = round(social_vulnerability, 2),
        priority_rank) %>% head(20)
    readr::write_csv(priority_top, file.path(OUT_REPORTS, "priority_table.csv"))
    message("[priority] Top-5 vulnerability-adjusted preparedness-priority zones:")
    print(as.data.frame(head(priority_top, 5)))
  }
  message("[risk] Top-5 at-risk Ituri zones (next-week first-case probability):")
  print(as.data.frame(head(ituri_top, 5)))
}
calibrated_forecasts <- NULL; calibration_gain <- NULL; fc_ensemble <- NULL

# ---------------------------------------------------------------------------
# 12b. Bayesian renewal suite (Task 6): posterior parameters + posterior
#      invasion probabilities WITH credible intervals, across mobility/covariate
#      assumptions, combined by loo predictive stacking. EXTENDS the frequentist
#      suite; the current-week fit is separate from fc_all_current.
# ---------------------------------------------------------------------------
bayes_suite <- NULL; bayes_stack <- NULL; bayes_risk_scores <- NULL; bayes_stack_risk_scores <- NULL
bayes_best_preds <- NULL   # featured single Bayesian model's current forecast (set inside the block)
if (requireNamespace("brms", quietly = TRUE)) {
  message("\n=== Step 12b: Bayesian renewal suite (brms) ===")
  bayes_suite <- tryCatch(
    fit_bayes_suite(zone_week_nc, mobility_matrices, gt_pmfs, covariates, osrm_mat,
                    zones_all, affected_now, LFO_HORIZONS, iter = 2000L),
    error = function(e) { warning("bayes suite: ", e$message); NULL })
  if (!is.null(bayes_suite)) {
    bayes_stack <- tryCatch(bayes_stacked_prediction(bayes_suite), error = function(e) NULL)
    saveRDS(bayes_suite$preds,  file.path(OUT_FORECASTS, "bayes_current_predictions.rds"))
    saveRDS(bayes_suite$params, file.path(OUT_FORECASTS, "bayes_parameters.rds"))
    if (!is.null(bayes_suite$weights))
      saveRDS(bayes_suite$weights, file.path(OUT_FORECASTS, "bayes_stacking_weights.rds"))
    if (!is.null(bayes_stack)) saveRDS(bayes_stack, file.path(OUT_FORECASTS, "bayes_stacked_current.rds"))
    readr::write_csv(bayes_suite$params, file.path(OUT_REPORTS, "bayes_parameters.csv"))
    # Convergence diagnostics per fitted model (review §4.4): divergent transitions,
    # bulk/tail ESS, max Rhat, post-warmup draws. Previously only Rhat was surfaced.
    bayes_diag <- tryCatch(
      dplyr::bind_rows(lapply(names(bayes_suite$fits), function(nm)
        bayes_fit_diagnostics(bayes_suite$fits[[nm]], nm))),
      error = function(e) NULL)
    if (!is.null(bayes_diag) && nrow(bayes_diag)) {
      readr::write_csv(bayes_diag, file.path(OUT_REPORTS, "bayes_convergence_diagnostics.csv"))
      message(sprintf("[bayes] convergence (§4.4): max Rhat=%.3f, min bulk-ESS=%.0f, total divergences=%d across %d models",
                      suppressWarnings(max(bayes_diag$rhat_max, na.rm = TRUE)),
                      suppressWarnings(min(bayes_diag$ess_bulk_min, na.rm = TRUE)),
                      suppressWarnings(sum(bayes_diag$n_divergent, na.rm = TRUE)), nrow(bayes_diag)))
    }
    # FEATURED single Bayesian model = the best cross-validated single model by the CV composite
    # (AUC-PR skill + mean rank-of-truth + log-score), already selected above from eval_tbl. The loo
    # predictive-stacking WEIGHTS are still computed and saved — they define the loo-stacked ENSEMBLE
    # (bayes_ensemble_*) and the bayes_stacking_weights figure — but they NO LONGER pick the featured
    # single model. Guard: only if the CV pick somehow lacks a current-forecast row do we fall back to
    # the highest-weight model that does have one (then the stacked ensemble downstream).
    if ((is.null(best_bayes_method) || !best_bayes_method %in% bayes_suite$preds$method) &&
        !is.null(bayes_suite$weights) && length(bayes_suite$weights)) {
      .valid_w <- bayes_suite$weights[names(bayes_suite$weights) %in% bayes_suite$preds$method]
      if (length(.valid_w)) best_bayes_method <- names(.valid_w)[which.max(.valid_w)]
    }
    message(sprintf("[bayes] featured single model = %s (best CV composite: AUC-PR skill + mean rank + log-score)",
                    best_bayes_method %||% "n/a"))
    # Bayesian risk products come from this featured model; the loo-stacked posterior is
    # retained separately as the Bayesian ENSEMBLE analogue (bayes_stack_risk_scores).
    bayes_best_preds <- NULL
    if (!is.null(best_bayes_method) && best_bayes_method %in% bayes_suite$preds$method) {
      # DEPLOYED-FORECAST GT MARGINALISATION (review §2.1). The mobility KERNEL is selected
      # by cross-validation at the central literature GT (so GT is NOT selected — the
      # reviewer's point); here we propagate generation-time uncertainty into the DEPLOYED
      # featured forecast by marginalising over GT_PRIOR (refit at each GT grid point, draws
      # pooled). This widens the per-zone credible intervals and gives the ranking a genuine
      # (non-degenerate) rank credible interval, without re-opening the CV skill. Gated by
      # GT_MARGINALISE_FEATURED (default TRUE); falls back to the fixed-GT prediction on any error.
      .fm_i <- which(vapply(bayes_suite$grid, function(g) identical(g$label, best_bayes_method), logical(1)))
      .fm   <- if (length(.fm_i)) bayes_suite$grid[[.fm_i[1]]] else NULL
      if (!is.null(.fm) && isTRUE(get0("GT_MARGINALISE_FEATURED", ifnotfound = TRUE))) {
        bayes_best_preds <- tryCatch({
          mp <- predict_bayes_gt_marginal(zone_week_nc, mobility_matrices, covariates, osrm_mat,
                  zones_all, mob = .fm$mob, cov_spec = .fm$cov, horizons = LFO_HORIZONS,
                  affected_zones = affected_now, gt_prior = GT_PRIOR, link = .fm$link %||% "cloglog")
          if (!is.null(mp)) {
            mp$method <- best_bayes_method
            message(sprintf("[gt-marginal] deployed featured forecast marginalised over GT_PRIOR (%d grid points) for %s",
                            length(make_gt_prior_pmfs(GT_PRIOR)$gt_pmfs), best_bayes_method))
          }
          mp
        }, error = function(e) { warning("[gt-marginal] featured marginalisation failed; using fixed-GT: ", e$message); NULL })
      }
      if (is.null(bayes_best_preds))
        bayes_best_preds <- bayes_suite$preds %>% dplyr::filter(method == best_bayes_method)
    } else if (!is.null(bayes_stack)) bayes_best_preds <- bayes_stack   # fallback
    if (!is.null(bayes_best_preds)) {
      bayes_risk_scores <- tryCatch({
        rs <- compute_risk_scores(bayes_best_preds, province_map)
        if (!is.null(vuln_index)) rs <- add_risk_indices(rs, vuln_index)
        rs
      }, error = function(e) { warning("bayes risk scores: ", e$message); NULL })
      if (!is.null(bayes_risk_scores)) {
        saveRDS(bayes_risk_scores, file.path(OUT_FORECASTS, "bayes_risk_scores_current.rds"))
        # Full per-zone Bayesian invasion probabilities + risk scores (all zones) to CSV.
        tryCatch(.write_risk_csv(bayes_risk_scores, file.path(OUT_REPORTS, "bayes_risk_scores_all_zones.csv")),
                 error = function(e) warning("bayes risk CSV: ", e$message))
        # Pairwise DIRECTED importation pressure / force of infection between zones,
        # decomposed from the featured model's renewal equation at its posterior-median
        # import coefficient. sum over origins of foi == the featured model's per-week
        # destination hazard, so this is the network behind the saved mu_forecast.
        tryCatch({
          .fm_i <- which(vapply(bayes_suite$grid,
                                function(g) identical(g$label, best_bayes_method), logical(1)))
          .fm <- if (length(.fm_i)) bayes_suite$grid[[.fm_i[1]]] else NULL
          .pr <- bayes_suite$params
          .beta_med <- if (!is.null(.pr)) {
            v <- .pr$hr[.pr$model == best_bayes_method & .pr$is_intercept]; if (length(v)) v[1] else NA_real_
          } else NA_real_
          if (!is.null(.fm) && is.finite(.beta_med)) {
            if (length(.fm$cov))
              warning(sprintf("[bayes] pairwise FOI: featured model %s carries covariates; pairwise decomposition uses the intercept-only posterior-median beta (per-zone covariate modulation omitted).",
                              best_bayes_method))
            pw <- bayes_pairwise_import_force(
              zone_week_nc, mobility_matrices, gt_pmfs, zones_all,
              mob = .fm$mob, gt = .fm$gt, beta_med = .beta_med, horizons = LFO_HORIZONS,
              beta_proj = 0.05, affected_zones = affected_now, province_map = province_map)
            readr::write_csv(pw, file.path(OUT_REPORTS, "bayes_pairwise_import_force.csv"))
            message(sprintf("[bayes] wrote bayes_pairwise_import_force.csv: %d directed pairs from %s (beta_med=%.3f)",
                            nrow(pw), best_bayes_method, .beta_med))
          }
        }, error = function(e) warning("bayes pairwise import-force CSV: ", e$message))
      }
    }
    # Bayesian ENSEMBLE (loo-stacked) risk scores — the analogue of the frequentist
    # ensemble, kept distinct from the best-single-model maps above.
    bayes_stack_risk_scores <- NULL
    if (!is.null(bayes_stack)) {
      bayes_stack_risk_scores <- tryCatch({
        rs <- compute_risk_scores(bayes_stack, province_map)
        if (!is.null(vuln_index)) rs <- add_risk_indices(rs, vuln_index)
        rs
      }, error = function(e) NULL)
      if (!is.null(bayes_stack_risk_scores))
        tryCatch(.write_risk_csv(bayes_stack_risk_scores, file.path(OUT_REPORTS, "bayes_ensemble_risk_scores_all_zones.csv")),
                 error = function(e) warning("bayes ensemble risk CSV: ", e$message))
    }
    message(sprintf("[bayes] fitted %d models; stacking weights: %s",
                    length(bayes_suite$fits),
                    if (!is.null(bayes_suite$weights))
                      paste(sprintf("%s=%.2f", names(bayes_suite$weights), bayes_suite$weights),
                            collapse = ", ") else "n/a"))
  }
}

# ---------------------------------------------------------------------------
# 12c. Daily re-issue — anchor ALL current forecasts to rolling day-windows FROM the
# analysis date (P(first case within 7 / 14 days of ANALYSIS_DATE)) and upsert into the
# accumulating daily series. ONE combined call (frequentist current forecast + featured
# Bayesian model + Bayesian stacked ensemble): the persistence keys on forecast_date, so a
# single call is correct where two would clobber. Frequentist rows carry a proper current-week
# rate (mu_wk0); the Bayesian rows have none, so the anchoring uses the documented continuity
# fallback (current-week rate = next-week rate) — the same treatment as the frequentist
# comparators that lack mu_wk0. In a Bayesian-only run this is how the daily operational product
# is produced; if no forecast at all is available (no frequentist models AND no Bayesian suite)
# it is simply skipped.
# ---------------------------------------------------------------------------
.reissue_fc <- dplyr::bind_rows(
  if (nrow(fc_all_current) > 0)   fc_all_current   else NULL,
  if (!is.null(bayes_best_preds)) bayes_best_preds else NULL,
  # Only add the stacked ensemble when it is a DISTINCT object from the featured single
  # model. best_bayes_method can fall back to the stack (line ~781), making bayes_best_preds
  # the SAME object as bayes_stack; binding both would duplicate every (method, health_zone,
  # horizon) key, which makes anchor_windows_from_analysis_date()'s pivot_wider emit list-
  # columns and abort — silently dropping the ENTIRE daily re-issue via its non-fatal tryCatch.
  if (!is.null(bayes_stack) && !identical(bayes_best_preds, bayes_stack)) bayes_stack else NULL)
# Defensive net: guarantee unique (method, health_zone, horizon) keys so the anchor's
# pivot_wider can never collapse duplicates into list-columns regardless of upstream overlap.
if (nrow(.reissue_fc) > 0 && all(c("method", "health_zone", "horizon") %in% names(.reissue_fc)))
  .reissue_fc <- dplyr::distinct(.reissue_fc, method, health_zone, horizon, .keep_all = TRUE)
if (nrow(.reissue_fc) > 0)
  invisible(tryCatch(
    issue_daily_reissue(.reissue_fc, forecast_date = ANALYSIS_DATE,
                        training_cutoff = training_cutoff),
    error = function(e) warning("[run_all] daily re-issue failed (non-fatal): ",
                                conditionMessage(e))))

# ---------------------------------------------------------------------------
# 13. Visualisations
# ---------------------------------------------------------------------------
.phase("Step 11-12 eval + risk + Bayesian suite")
message("\n=== Step 13: Visualisations (legible invasion suite) ===")
source(file.path(ST_DIR, "17_invasion_viz.R"))
# 20_forecast_detail.R already sourced in step 12 (risk indices + priority + viz)

# Persist the in-memory forecast objects so the detail visualisations (and any
# ad-hoc re-plotting) can run without repeating the expensive forecast/LFO step.
saveRDS(fc_all_current, file.path(OUT_FORECASTS, "fc_all_current.rds"))
if (!is.null(lfo_results)) saveRDS(lfo_results, file.path(OUT_FORECASTS, "lfo_results.rds"))

# Remove figures/tables from the superseded pre-rebuild pipeline so the outputs
# folder shows ONLY the current invasion products (no maps of removed models).
.stale <- c(
  file.path(OUT_DIAGNOSTICS, c("brier_comparison.pdf", "calibration_curves.pdf",
    "case_count_tile.pdf", "evaluation_heatmap_h1.pdf", "evaluation_heatmap_h2.pdf",
    "pr_curves.pdf", "roc_curves.pdf", "wis_by_model_1w.pdf",
    "metrics_summary.csv", "evaluation_heatmap.pdf")),
  list.files(OUT_MAPS, "^risk_map_", full.names = TRUE),
  file.path(OUT_REPORTS, c("methodology_report.md", "summary_dashboard.pdf",
    "table_best_per_metric.csv", "table_model_comparison.csv",
    "table_outbreak_summary.csv")))
suppressWarnings(file.remove(.stale[file.exists(.stale)]))

shapefile <- if (file.exists(SHAPEFILE_PATH))
  tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE), error = function(e) NULL) else NULL
vwrap <- function(expr) tryCatch(expr, error = function(e) warning(conditionMessage(e)))

# F1 — legible epidemic curve (windowed; top zones + Other)
vwrap(plot_epidemic_curve_legible(zone_week_nc, province_map, save = TRUE))
# F9 — epinowcast nowcast fan
vwrap(plot_nowcast_fan2(file.path(OUT_DIAGNOSTICS, "epinowcast_weekly_factors.csv")))
# National R(t) from the renewal estimate (with GT-profile sensitivity)
vwrap(plot_rt(Rt_primary, rt_all = Rt_national_list, save = TRUE))

# Task 2 — explicit fit/prediction date window annotated on every current-forecast
# figure (fit window = training data; prediction window = weeks being forecast).
# Use the OUTBREAK training start, not min(all_weeks): the zone-week grid carries
# pre-outbreak historical weeks (all structurally zero for this outbreak), so the
# epidemiologically meaningful fit window begins at the outbreak start.
.fit_start <- suppressWarnings(min(zone_week_outbreak$week_start, na.rm = TRUE))
if (!is.finite(.fit_start)) .fit_start <- OUTBREAK_START
window_txt <- .window_caption(.fit_start, training_cutoff, LFO_HORIZONS)

# ONE map/decision-SUITE renderer, called once per FEATURED model so the frequentist
# and Bayesian paradigms get an IDENTICAL, model-labelled set of figures (Tasks 1-4):
# invasion maps (national + Ituri zoom, per horizon), per-province invasion maps
# (Ituri/Nord-Kivu/Haut-Uele), prob x vulnerability choropleths (national + provinces),
# ranked at-risk bars, and the preparedness-priority scatter/bars/map. No arbitrary
# model picks; every figure's title names the exact model that produced it, and the
# Bayesian set is file-prefixed "bayes_".
render_map_suite <- function(rs, model_label, file_prefix = "") {
  if (is.null(rs)) return(invisible(NULL))
  pf <- function(b) if (nzchar(file_prefix)) paste0(file_prefix, b) else b
  # EVERY map/decision product is produced at BOTH horizons (1- and 2-week ahead).
  for (h in LFO_HORIZONS) {
    vwrap(plot_invasion_risk_map(rs, horizon = h, method_label = model_label,
            shapefile = shapefile, window_txt = window_txt, save = TRUE,
            file = pf("invasion_risk_map")))
    # Invasion-probability RANKING map (1 = highest risk) — robust operational targeting view,
    # complementing the absolute-probability map above. Produced per model at both horizons, as
    # the Ituri+national pair AND a standalone national-only panel.
    vwrap(plot_invasion_rank_map(rs, horizon = h, method_label = model_label,
            shapefile = shapefile, window_txt = window_txt, save = TRUE,
            file = pf("invasion_rank_map"), extent = "both"))
    vwrap(plot_invasion_rank_map(rs, horizon = h, method_label = model_label,
            shapefile = shapefile, window_txt = window_txt, save = TRUE,
            file = pf("invasion_rank_map"), extent = "national"))
    # #2 — probability WITH uncertainty (produced when p_lo/p_hi are present: the
    # Bayesian posterior CrI, or the frequentist ensemble spread attached below).
    # Both the Ituri zoom AND the whole-DRC national extent (#4, side-by-side prob+uncertainty).
    vwrap(plot_invasion_uncertainty_map(rs, shapefile = shapefile, horizon = h,
            model_label = model_label, window_txt = window_txt, save = TRUE,
            file = pf("invasion_uncertainty_map"), extent = "ituri"))
    vwrap(plot_invasion_uncertainty_map(rs, shapefile = shapefile, horizon = h,
            model_label = model_label, window_txt = window_txt, save = TRUE,
            file = pf("invasion_uncertainty_map_national"), extent = "national"))
    vwrap(plot_province_risk_maps(rs, shapefile = shapefile, horizon = h,
            method_label = model_label, window_txt = window_txt, save = TRUE,
            file = pf("invasion_risk_map")))
    vwrap(plot_risk_scores_bars(rs, horizon = h, save = TRUE,
            file = pf("risk_scores_bars"), model_label = model_label))
    if ("V" %in% names(rs)) {
      vwrap(plot_prob_vuln_choropleth(rs, shapefile = shapefile, horizon = h,
              province_zoom = NULL, window_txt = window_txt, save = TRUE,
              file = pf(sprintf("prob_vuln_choropleth_national_h%d", h)), model_label = model_label))
      for (prov in PROVINCES_OF_INTEREST)
        vwrap(plot_prob_vuln_choropleth(rs, shapefile = shapefile, horizon = h,
                province_zoom = prov, window_txt = window_txt, save = TRUE,
                file = pf(sprintf("prob_vuln_choropleth_%s_h%d", .prov_suffix(prov), h)),
                model_label = model_label))
    }
    if ("priority" %in% names(rs)) {
      vwrap(plot_priority_scatter(rs, horizon = h, window_txt = window_txt, save = TRUE,
              file = pf("priority_scatter"), model_label = model_label))
      vwrap(plot_priority_bars(rs, horizon = h, save = TRUE,
              file = pf("priority_bars"), model_label = model_label))
      vwrap(plot_priority_map(rs, shapefile = shapefile, horizon = h, window_txt = window_txt,
              save = TRUE, file = pf("priority_map"), model_label = model_label))
    }
  }
  invisible(TRUE)
}
# Best FREQUENTIST model (primary_method — skill-selected).
if (!is.null(risk_scores))
  render_map_suite(risk_scores, sprintf("Frequentist: %s", primary_method), "")

# Q1 — key parameter estimates: the mobility import coefficient (dominant driver)
# plus a covariate association screen (the full candidate range, each estimable
# marginally and annotated by what survives adjustment for the mobility offset).
if (.have_freq && !is.na(primary_method) && grepl("^Renewal", primary_method)) {
  .best_mob <- sub("^Renewal-(M[0-9]+[a-z]?).*$", "\\1", primary_method)
  if (!grepl("^M", .best_mob)) .best_mob <- "M8"
  .best_gt  <- if (grepl("-short", primary_method)) "short" else
               if (grepl("-long", primary_method)) "long" else "medium"
  inv_design <- tryCatch(build_invasion_design(zone_week_nc, mobility_matrices, gt_pmfs,
                           covariates, osrm_mat, zones_all, mob = .best_mob, gt = .best_gt),
                         error = function(e) { warning("invasion design: ", e$message); NULL })
  if (!is.null(inv_design)) {
    cov_assoc <- tryCatch(covariate_associations(inv_design),
                          error = function(e) { warning("covariate assoc: ", e$message); NULL })
    saveRDS(cov_assoc, file.path(OUT_FORECASTS, "covariate_associations.rds"))
    vwrap(plot_model_parameters(cov_assoc, primary_method, beta0 = inv_design$beta0,
                                n_events = inv_design$n_events, save = TRUE))
  }
}

# (The frequentist preparedness-priority products — scatter, component bars, priority
# map — are rendered by render_map_suite() above alongside the invasion maps and
# choropleths, so both paradigms produce the identical labelled set.)

# Evaluation — curated display set (best per family + ensembles + comparators).
if (!is.null(eval_display)) {
  for (h in LFO_HORIZONS) vwrap(plot_discrimination_summary(eval_display, horizon = h, save = TRUE))
}
# Bayesian-only discrimination summary (ALL cross-validated Bayesian models, from the full
# eval_tbl) so the Bayesian grid can be compared amongst itself.
if (!is.null(eval_tbl) && any(grepl("^Bayes", eval_tbl$method))) {
  for (h in LFO_HORIZONS)
    vwrap(plot_discrimination_summary(eval_tbl, horizon = h, restrict = "^Bayes",
            label = "Bayesian", file = "bayes_discrimination_summary", save = TRUE))
}
# #8 — models vs baselines/simple models: discrimination + top-10 detection over the SAME
# folds, every model coloured by family, naive-baseline band + random watch-list marked. Uses
# the FULL eval_tbl (all models, not just the curated display set) so the gap is complete.
if (!is.null(eval_tbl)) {
  for (h in LFO_HORIZONS) vwrap(plot_model_vs_baseline(eval_tbl, horizon = h, save = TRUE))
}
# Task 4 — the per-model spatial / space-time diagnostics are rendered for BOTH the
# best FREQUENTIST (primary_method) and best BAYESIAN (best_bayes_method) model, each
# labelled and file-tagged (bayes_ prefix), recomputing the per-method spatial error.
# Frequentist diagnostic model only when frequentist models ran; the Bayesian entry is
# appended so a Bayesian-only run still produces its single-model spatial diagnostics.
.diag_models <- if (.have_freq && !is.na(primary_method))
  list(list(m = primary_method, lab = sprintf("Frequentist: %s", primary_method), pfx = "")) else list()
if (!is.null(best_bayes_method) && best_bayes_method %in% lfo_results$method)
  .diag_models[[length(.diag_models) + 1L]] <- list(m = best_bayes_method,
    lab = sprintf("Bayesian: %s", best_bayes_method), pfx = "bayes_")
if (!is.null(lfo_results)) {
  # NAIVE comparator for the prioritisation (figure-3, panel B) curves: rank at-risk zones
  # purely by mobility inflow FROM THE EPICENTRE, ignoring case data. Injected as an extra
  # "method" on the SAME folds/outcomes so its detection curve overlays the model's — the
  # model must beat this structural, incidence-free baseline. Scoped to figure 3 (a local
  # augmented copy), so it does not enter the ranking / skill-over-time / balance products.
  .naive_lbl    <- "Naive-epicentre-inflow"
  .naive_scores <- tryCatch(naive_epicentre_inflow_scores(W_primary, EPICENTRE_ZONES,
                              pop_vec, zones_all), error = function(e) {
                              warning("[naive] epicentre-inflow scores failed: ", conditionMessage(e)); NULL })
  lfo_fig3 <- if (!is.null(.naive_scores))
    append_naive_detection_curve_model(lfo_results, .naive_scores, .naive_lbl) else lfo_results
  .naive_ok <- !is.null(.naive_scores) && .naive_lbl %in% lfo_fig3$method
  # All per-model evaluation figures produced at BOTH horizons (1- and 2-week ahead).
  for (hh in LFO_HORIZONS) {
    vwrap(plot_reliability(lfo_results, horizon = hh, save = TRUE))
    # Bayesian calibration/reliability + combined figure-3, featured Bayesian model. Panel A
    # (reliability) uses methods[1] = the featured Bayesian model; panel B (prioritisation)
    # overlays the naive epicentre-inflow baseline for a like-for-like comparison.
    if (!is.null(best_bayes_method) && best_bayes_method %in% lfo_results$method) {
      vwrap(plot_reliability(lfo_results, horizon = hh, methods = best_bayes_method,
              file = "bayes_reliability", save = TRUE))
      vwrap(plot_paper_figure3(lfo_fig3,
              if (.naive_ok) c(best_bayes_method, .naive_lbl) else best_bayes_method,
              horizon = hh, save = TRUE, file = "bayes_model_performance_figure3"))
    }
    for (dm in .diag_models) {
      vwrap(plot_lfo_forecast_vs_outcome(lfo_results, method = dm$m, horizon = hh, save = TRUE,
              file = paste0(dm$pfx, "lfo_forecast_vs_outcome")))
      vwrap(plot_spacetime_risk(lfo_results, method = dm$m, horizon = hh, save = TRUE,
              file = paste0(dm$pfx, "spacetime_risk")))
      ze <- tryCatch(zone_spatial_error(lfo_results, province_map, horizon = hh, method = dm$m),
                     error = function(e) NULL)
      if (!is.null(ze)) vwrap(plot_spatial_error_map(ze, shapefile = shapefile, save = TRUE,
              file = paste0(dm$pfx, "spatial_error_map"), model_label = dm$lab, horizon = hh))
    }
  }
}

# Q5 — spatiotemporal evaluation figures: skill-over-time (all models) + lead-time.
if (!is.null(st_skill)) {
  vwrap(plot_skill_over_time(st_skill, metric = "auc_pr_skill", save = TRUE))
  vwrap(plot_skill_over_time(st_skill, metric = "hit_at_k", save = TRUE))
}
if (!is.null(lead_tbl)) vwrap(plot_lead_time(lead_tbl, save = TRUE))

# Q1 — the ENSEMBLE visualised: a side-by-side best-model-vs-ensemble panel with a
# spatial member-disagreement (uncertainty) map. (A standalone ensemble call to
# plot_invasion_risk_map would clobber the best-model map — it writes a fixed,
# method-agnostic filename — so the ensemble is shown via the dedicated panel.)
ens_method <- if (!is.null(eval_tbl) &&
                  "Ensemble-mean" %in% eval_tbl$method) "Ensemble-mean" else NULL
if (!is.null(ens_method)) {
  vwrap(plot_forecast_map_panel(fc_all_current, members = ENSEMBLE_MEMBERS,
          primary_method = primary_method, ensemble_method = ens_method,
          shapefile = shapefile, horizon = 1, save = TRUE))
}

# Q4 — forecasts for the best-fitting models WITH uncertainty (member spread):
#   current top-zone probabilities with min–max bands, and the space-time
#   trajectories (predicted probability over folds) for zones that invaded.
# (These are the FREQUENTIST-ensemble member-spread products; skipped in a Bayesian-only
# run — the Bayesian posterior 90% CrI maps play the equivalent role there.)
if (.have_freq) {
  vwrap(plot_forecast_uncertainty(fc_all_current, members = ENSEMBLE_MEMBERS,
          primary_method = primary_method, province_map = province_map,
          horizon = 1, window_txt = window_txt, save = TRUE))
  if (!is.null(lfo_results))
    vwrap(plot_spacetime_forecast(lfo_results, members = ENSEMBLE_MEMBERS,
            province_map = province_map, horizon = 1, save = TRUE))
}

# Task 6 — Bayesian suite figures: posterior parameter forest (with CrI, across
# mobility/covariate assumptions), posterior invasion probabilities with 90% CrI
# (featured model + stacked ensemble), and stacking weights.
if (!is.null(bayes_suite)) {
  vwrap(plot_bayes_parameters(bayes_suite$params, weights = bayes_suite$weights,
                              window_txt = window_txt, save = TRUE))
  # Highlight the SKILL-SELECTED Bayesian model. If best_bayes_method is somehow invalid, fall
  # back DATA-DRIVENLY to the highest loo-stacking-weight model that has predictions (then the
  # first available) — never a hard-coded model name — so the figure always tracks the best fit.
  .feat_bayes <- best_bayes_method
  if (is.null(.feat_bayes) || !.feat_bayes %in% bayes_suite$preds$method) {
    .wv <- bayes_suite$weights
    .cand <- if (!is.null(.wv) && length(.wv))
               intersect(names(sort(.wv, decreasing = TRUE)), bayes_suite$preds$method) else character(0)
    .feat_bayes <- if (length(.cand)) .cand[1] else unique(bayes_suite$preds$method)[1]
  }

  # Task 2 — FULL posterior distributions of the Bayesian parameters (densities, not just
  # median + CrI), extracted from the fitted models.
  vwrap(plot_bayes_posterior_densities(bayes_posterior_draws(bayes_suite$fits),
          window_txt = window_txt, save = TRUE))
  # Tasks 2 & 3 — REFIT-based analyses for the featured model (gated; each adds ~10 brms fits):
  #   (2) a loo-predictive POSTERIOR over the generation-time mean (the GT is otherwise a fixed
  #       assumption), and (3) a SENSITIVITY of beta0 to the two-stage nowcast INPUT (raw vs
  #       epinowcast vs fast delay-CDF training counts) — a pragmatic check on whether feeding
  #       nowcast-corrected counts as a fixed second stage drives the inference.
  if (isTRUE(get0("BAYES_PROFILE_ANALYSES", ifnotfound = TRUE))) {
    .feat_spec <- Filter(function(g) g$label == .feat_bayes, bayes_default_grid(mobility_matrices))
    .feat_spec <- if (length(.feat_spec)) .feat_spec[[1]]
                  else list(mob = "M8", gt = "medium", cov = character(0), link = "cloglog")
    .gtp <- tryCatch(bayes_gt_posterior(zone_week_nc, mobility_matrices, gt_pmfs, covariates,
              osrm_mat, zones_all, mob = .feat_spec$mob, cov = .feat_spec$cov,
              link = .feat_spec$link %||% "cloglog"), error = function(e) { warning("bayes GT posterior: ", conditionMessage(e)); NULL })
    vwrap(plot_bayes_gt_posterior(.gtp, model_label = .feat_bayes, window_txt = window_txt))
    if (!is.null(.gtp)) readr::write_csv(.gtp, file.path(OUT_FORECASTS, "bayes_gt_posterior.csv"))
    .zw_raw <- dat$zone_week; .zw_raw$confirmed_nc <- .zw_raw$confirmed
    .zw_variants <- list(raw = .zw_raw, epinowcast = zone_week_nc)
    .zw_fast <- tryCatch({ z <- apply_nowcast_correction(dat$zone_week, analysis_date = ANALYSIS_DATE)
                           if (!"confirmed_nc" %in% names(z)) z$confirmed_nc <- z$confirmed; z },
                         error = function(e) NULL)
    if (!is.null(.zw_fast)) .zw_variants[["fast"]] <- .zw_fast
    .sens <- tryCatch(bayes_nowcast_sensitivity(.zw_variants, mobility_matrices, gt_pmfs, covariates,
               osrm_mat, zones_all, mob = .feat_spec$mob, gt = .feat_spec$gt,
               cov = .feat_spec$cov, link = .feat_spec$link %||% "cloglog"),
               error = function(e) { warning("bayes nowcast sensitivity: ", conditionMessage(e)); NULL })
    vwrap(plot_bayes_nowcast_sensitivity(.sens, model_label = .feat_bayes, window_txt = window_txt))
    if (!is.null(.sens)) readr::write_csv(.sens, file.path(OUT_FORECASTS, "bayes_nowcast_sensitivity.csv"))
  }
  for (hh in LFO_HORIZONS) {
    vwrap(plot_bayes_invasion_uncertainty(bayes_suite$preds, province_map = province_map,
            horizon = hh, model = .feat_bayes, window_txt = window_txt, save = TRUE,
            file = "bayes_invasion_uncertainty"))
    if (!is.null(bayes_stack))
      vwrap(plot_bayes_invasion_uncertainty(bayes_stack, province_map = province_map,
              horizon = hh, model = NULL, window_txt = window_txt, save = TRUE,
              file = "bayes_stacked_invasion_uncertainty"))
  }
  vwrap(plot_bayes_stacking(bayes_suite$weights, save = TRUE))

  # Tasks 1-4 — the FULL Bayesian map/decision suite, an exact analogue of the
  # frequentist set, from the SKILL-SELECTED best Bayesian model (best_bayes_method),
  # file-prefixed "bayes_" and labelled with the model. Plus the loo-STACKED Bayesian
  # ENSEMBLE as the analogue of the frequentist ensemble ("bayes_ensemble_").
  # Only render the best-single-model suite when a Bayesian model was actually
  # skill-selected; otherwise bayes_risk_scores == the stack and would duplicate the
  # ensemble suite below with an unparseable label.
  .map_jobs <- list()
  if (!is.null(bayes_risk_scores) && !is.null(best_bayes_method))
    .map_jobs[[length(.map_jobs) + 1L]] <- list(
      rs = bayes_risk_scores,
      label = sprintf("Bayesian: %s", .bayes_model_label(best_bayes_method)),
      prefix = "bayes_")
  if (!is.null(bayes_stack_risk_scores))
    .map_jobs[[length(.map_jobs) + 1L]] <- list(
      rs = bayes_stack_risk_scores, label = "Bayesian loo-stacked ensemble",
      prefix = "bayes_ensemble_")
  if (length(.map_jobs)) {
    .map_one <- function(j) {
      tryCatch(render_map_suite(j$rs, j$label, j$prefix),
               error = function(e) warning(conditionMessage(e)))
      TRUE
    }
    .map_par <- get0("PARALLEL_JOBS", ifnotfound = 1L) > 1L &&
                length(.map_jobs) > 1L && requireNamespace("furrr", quietly = TRUE)
    .tv <- Sys.time()
    if (.map_par)
      invisible(furrr::future_map(.map_jobs, .map_one,
        .options = furrr::furrr_options(seed = TRUE)))
    else invisible(lapply(.map_jobs, .map_one))
    message(sprintf("[viz-timing] %d static Bayesian map suites %s in %.1fs",
                    length(.map_jobs), if (.map_par) "PARALLEL" else "seq",
                    as.numeric(difftime(Sys.time(), .tv, units = "secs"))))
  }
}

# Time-evolution GIFs are intentionally not generated. They duplicated the static
# decision products while dominating the visualisation runtime. Remove artifacts
# from earlier runs so the output directory cannot retain stale animations.
.anim_dir <- file.path(OUT_MAPS, "animations")
if (dir.exists(.anim_dir)) {
  .old_anim <- list.files(.anim_dir, pattern = "\\.(gif|pdf)$",
                          full.names = TRUE, ignore.case = TRUE)
  if (length(.old_anim)) unlink(.old_anim)
}

# Task 0b — Frequentist vs Bayesian comparison (separation + robustness): model
# rankings on the same folds, and per-zone agreement of the two paradigms.
# (Only meaningful when BOTH paradigms ran; skipped in a Bayesian-only run.)
if (.have_freq && !is.null(risk_scores)) {
  for (hh in LFO_HORIZONS) {
    vwrap(plot_freq_bayes_ranking(eval_tbl, horizon = hh, save = TRUE))
    if (!is.null(bayes_stack))
      vwrap(plot_freq_bayes_agreement(risk_scores, bayes_stack, province_map = province_map,
              horizon = hh, save = TRUE))
  }
}

# Task 4 — evaluation over time: parameter estimates over folds, and per-fold
# predicted-vs-observed invasion for the featured models.
if (!is.null(lfo_results)) {
  .fold_cutoffs <- sort(unique(as.Date(lfo_results$cutoff)))
  pot <- tryCatch(compute_params_over_time(zone_week_outbreak, .fold_cutoffs,
           mobility_matrices, gt_pmfs, covariates, osrm_mat, zones_all,
           nowcast_fn = apply_nowcast_correction),
           error = function(e) { warning("params-over-time: ", e$message); NULL })
  if (!is.null(pot)) vwrap(plot_params_over_time(pot, save = TRUE))
  for (hh in LFO_HORIZONS) {
    vwrap(plot_predobs_over_folds(lfo_results, method = best_method, horizon = hh, save = TRUE))
    if (!is.null(best_bayes_method) && best_bayes_method %in% lfo_results$method)
      vwrap(plot_predobs_over_folds(lfo_results, method = best_bayes_method, horizon = hh, save = TRUE,
              file = "bayes_predicted_vs_observed_over_folds"))
  }

  # Task 6 — intuitive detection-vs-budget curve + balanced skill metrics for the
  # imbalanced invasion task (sensitivity at a fixed alert budget; balanced
  # accuracy / F1 / MCC at the Youden-optimal threshold).
  .skill_methods <- intersect(unique(c(best_method, primary_method, best_bayes_method)),
                              unique(lfo_results$method))
  # The BAYESIAN models to compare among themselves in the top-K precision view (request 2):
  # the strongest cross-validated Bayesian models by AUC-PR skill, the featured single model,
  # and the Bayesian ensemble — so the top-K precision figure contrasts DIFFERENT Bayesian
  # models rather than a single featured one against the frequentist pack.
  .bayes_skill_methods <- if (!is.null(eval_tbl) && any(grepl("^Bayes", eval_tbl$method))) {
    .bt <- eval_tbl %>% dplyr::filter(horizon == 1L, grepl("^Bayes", method)) %>%
      dplyr::arrange(dplyr::desc(auc_pr_skill)) %>% dplyr::pull(method)
    intersect(unique(c(head(.bt, 6L), best_bayes_method,
                       "Bayes-ens-mean", "Bayes-ens-median")),
              unique(lfo_results$method))
  } else character(0)
  for (hh in LFO_HORIZONS) {
    vwrap(plot_detection_curve(lfo_results, .skill_methods, horizon = hh, save = TRUE))
    # #3 — precision: % of the top-K highest-risk zones that were actually invaded.
    vwrap(plot_topk_precision(lfo_results, .skill_methods, horizon = hh, save = TRUE))
    # Bayesian-only top-K precision: compare the leading Bayesian models + ensemble (request 2).
    if (length(.bayes_skill_methods) >= 2L)
      vwrap(plot_topk_precision(lfo_results, .bayes_skill_methods, horizon = hh, save = TRUE,
              file = "bayes_topk_precision", title_suffix = "Bayesian models"))
  }
  # Combined "Figure 3" briefing panel (calibration+AUC | prioritisation-vs-random),
  # after Kraemer & Cauchemez 2017 — the headline "how good is the model" figure. Panel B
  # also overlays the naive epicentre-mobility-inflow baseline (built above; lfo_fig3), so
  # the model's prioritisation is compared against pure epicentre connectivity.
  .fig3_lfo     <- if (exists("lfo_fig3")) lfo_fig3 else lfo_results
  .fig3_methods <- if (exists(".naive_ok") && isTRUE(.naive_ok))
    c(.skill_methods, .naive_lbl) else .skill_methods
  for (hh in LFO_HORIZONS)
    vwrap(plot_paper_figure3(.fig3_lfo, .fig3_methods, horizon = hh, save = TRUE))
  balance_tbl <- tryCatch(
    purrr::map_dfr(LFO_HORIZONS, function(hh)
      purrr::map_dfr(unique(lfo_results$method),
                     function(m) invasion_balance_metrics(lfo_results, m, hh))),
    error = function(e) NULL)
  if (!is.null(balance_tbl) && nrow(balance_tbl))
    readr::write_csv(balance_tbl, file.path(OUT_DIAGNOSTICS, "invasion_balanced_skill.csv"))

  # Task 4 — reporting rate across space, and the import coefficient beta0 +
  # reporting completeness across folds.
  vwrap(plot_reporting_rate_map(covariates, shapefile = shapefile, zones_all = zones_all, save = TRUE))
  bof <- tryCatch(compute_beta_over_folds(zone_week_outbreak, .fold_cutoffs,
           mobility_matrices, gt_pmfs, covariates, osrm_mat, zones_all,
           nowcast_fn = apply_nowcast_correction),
           error = function(e) { warning("beta-over-folds: ", e$message); NULL })
  if (!is.null(bof)) vwrap(plot_beta_over_folds(bof, save = TRUE))

  # Bayesian analogues of params-over-time + beta-over-folds: refit the featured
  # Bayesian model's mobility kernel (+ geo covariates, so covariate-HR traces exist)
  # at each fold cutoff and plot the POSTERIOR covariate HRs and import coefficient
  # beta0 over folds, with proper 90% CrI (bayes_params_over_time / bayes_beta_over_folds).
  if (!is.null(best_bayes_method) && requireNamespace("brms", quietly = TRUE)) {
    # Keep the "-dist" suffix so a road-distance featured kernel (e.g. Bayes-M4-dist-geo)
    # is refit as M4-dist, not the travel-time M4 — the over-time trace must be faithful
    # to the featured model's actual mobility kernel.
    .bmob <- sub("^Bayes-(M[0-9]+[a-z]?(-dist)?).*$", "\\1", best_bayes_method)
    if (!grepl("^M", .bmob)) .bmob <- "M8"
    bpot <- tryCatch(compute_bayes_params_over_time(zone_week_outbreak, .fold_cutoffs,
              mobility_matrices, gt_pmfs, covariates, osrm_mat, zones_all,
              mob = .bmob, gt = "medium", cov_spec = c("log_pop", "ccvi", "d_min"),
              iter = 600L, nowcast_fn = apply_nowcast_correction),
              error = function(e) { warning("bayes params-over-time: ", e$message); NULL })
    if (!is.null(bpot)) {
      if (!is.null(bpot$params)) vwrap(plot_params_over_time(bpot$params, save = TRUE,
              file = "bayes_params_over_time", model_label = sprintf("Bayesian: %s + geo", .bmob)))
      if (!is.null(bpot$beta)) vwrap(plot_beta_over_folds(bpot$beta, save = TRUE,
              file = "bayes_beta_over_folds",
              model_label = sprintf("Bayesian: %s", best_bayes_method), ci_label = "90% CrI"))
    }
  }
}

# ---------------------------------------------------------------------------
# 14. Invasion report
# ---------------------------------------------------------------------------
.phase("Step 13  visualisations")
message("\n=== Step 14: Invasion report ===")
vwrap(write_invasion_report(
  eval_tbl = eval_tbl,
  # Bayesian-only run: the featured risk table + primary method fall back to the Bayesian
  # products so the report's built-in NULL fallbacks engage (rather than a stray NA).
  risk_scores = if (!is.null(risk_scores)) risk_scores else get0("bayes_risk_scores"),
  lfo_results = lfo_results,
  zone_week = zone_week_nc, training_cutoff = training_cutoff,
  best_method = best_method,
  primary_method = if (length(primary_method) == 1 && !is.na(primary_method)) primary_method else NULL,
  n_models = if (!is.null(lfo_results)) dplyr::n_distinct(lfo_results$method)
             else dplyr::n_distinct(fc_all_current$method),
  best_bayes_method = best_bayes_method,
  bayes_weights = if (!is.null(bayes_suite)) bayes_suite$weights else NULL))

# Q2/Q3 — document HOW models are selected and WHAT goes into the best model
# (structure, mobility kernel, generation time, observation process, covariates,
# calibration, nowcast). Writes model_specification.md and appends to the report.
vwrap(write_model_details_report(best_method = best_method,
        primary_method = if (length(primary_method) == 1 && !is.na(primary_method)) primary_method else NULL,
        eval_tbl = eval_tbl,
        bayes_params = if (!is.null(bayes_suite)) bayes_suite$params else NULL,
        best_bayes_method = best_bayes_method,
        freq_beta0 = if (exists("inv_design")) inv_design$beta0 else NA_real_,
        freq_cov = if (exists("cov_assoc")) cov_assoc else NULL))

# ---------------------------------------------------------------------------
# Key outputs — gather the headline Bayesian deliverables (both horizons) into a single
# key_outputs/ folder for quick sharing, once everything above has been generated. Each entry
# is the h1 file; its h2 analogue (…_h1 -> …_h2) is copied too when present.
# ---------------------------------------------------------------------------
local({
  key_dir <- file.path(OUT_DIR, "key_outputs")
  dir.create(key_dir, showWarnings = FALSE, recursive = TRUE)
  h1_files <- c(
    "maps/bayes_invasion_rank_map_national_h1.pdf",
    "maps/bayes_invasion_uncertainty_map_national_h1.pdf",
    "maps/bayes_prob_vuln_choropleth_national_h1.pdf",
    "diagnostics/bayes_discrimination_summary_h1.pdf",
    "diagnostics/bayes_lfo_forecast_vs_outcome_h1.pdf",
    "diagnostics/bayes_predicted_vs_observed_over_folds_h1.pdf",
    "diagnostics/topk_precision_h1.pdf",
    "diagnostics/bayes_topk_precision_h1.pdf",
    "reports/bayes_invasion_uncertainty_h1.pdf",
    "reports/bayes_model_performance_figure3_h1.pdf",   # "model_performance_figures" == figure3
    "reports/priority_scatter_h1.pdf",                  # preparedness-priority scatter (frequentist featured)
    "reports/bayes_priority_scatter_h1.pdf")            # preparedness-priority scatter (Bayesian featured)
  no_horizon <- c("reports/bayes_risk_scores_all_zones.csv",
                  "reports/harmonised_confirmed_cases.csv")
  wanted <- unique(c(as.vector(rbind(h1_files, sub("_h1\\.", "_h2.", h1_files))), no_horizon))
  copied <- 0L; missing <- character(0)
  for (rel in wanted) {
    src <- file.path(OUT_DIR, rel)
    if (file.exists(src)) { file.copy(src, file.path(key_dir, basename(src)), overwrite = TRUE); copied <- copied + 1L }
    else missing <- c(missing, basename(rel))
  }
  message(sprintf("[key_outputs] copied %d files -> %s", copied, key_dir))
  if (length(missing)) message("[key_outputs] not yet generated (skipped): ", paste(missing, collapse = ", "))
})

# ---------------------------------------------------------------------------
# 15. Publication + manuscript figures, 3-month cascade, and Bayesian report refresh
# ---------------------------------------------------------------------------
# These are STANDALONE tools invoked as isolated Rscript subprocesses rather than
# source()d because (a) update_bayesian_report.R calls quit() on a no-op / --check,
# which would abort THIS session if sourced, and (b) each re-sources 00_config.R and
# would clobber run_all's in-memory globals. Each call is guarded so a failure can never
# fail the pipeline; the working directory is pinned to ROOT so make_publication_figures.R's
# normalizePath(".") resolves to the repo root.
#   - make_publication_figures.R / make_manuscript_figures.R / make_topk15_ever.R read ONLY
#     the saved outputs written above (no model re-fit), so they run after step 14 + the
#     key_outputs gather. make_topk15_ever.R reuses make_publication_figures.R's house style,
#     so it runs after it.
#   - run_cascade.R is the SELF-CONTAINED 3-month spatial-invasion cascade driver (modules
#     30–40): it re-sources data prep + the modelling layer and runs its own MC simulation,
#     so it is HEAVY (full run ~hours at CASCADE_N_MC=1000). Set CASCADE_SMOKE=1 for the fast
#     path, or SKIP_CASCADE=1 to omit it from this pipeline.
message("\n=== Step 15: Publication + manuscript figures, cascade, report refresh ===")
local({
  rscript <- file.path(R.home("bin"),
                       if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  old_wd <- setwd(ROOT); on.exit(setwd(old_wd), add = TRUE)
  run_tool <- function(rel, label, args = character()) {
    tryCatch({
      message("\n--- ", label, " (", rel, ") ---")
      st <- system2(rscript, c(shQuote(file.path(ROOT, rel)), args),
                    stdout = "", stderr = "")
      if (!identical(as.integer(st), 0L))
        message(sprintf("[%s] Rscript exited with status %s (non-fatal).", label, st))
    }, error = function(e)
      message(sprintf("[%s] skipped (non-fatal): %s", label, conditionMessage(e))))
  }
  run_tool("spatiotemporal/make_publication_figures.R", "publication figures")
  run_tool("spatiotemporal/make_manuscript_figures.R",  "manuscript figures")
  run_tool("spatiotemporal/make_topk15_ever.R",         "top-15 three-way outcome figure")
  run_tool("spatiotemporal/update_bayesian_report.R",   "Bayesian report refresh")
  # 3-month spatial-invasion cascade (self-contained; heavy). Opt out with SKIP_CASCADE=1.
  if (!(tolower(trimws(Sys.getenv("SKIP_CASCADE", ""))) %in% c("true", "t", "1", "yes", "y")))
    run_tool("spatiotemporal/run_cascade.R", "3-month invasion cascade")
  else
    message("[3-month invasion cascade] SKIP_CASCADE set — skipping cascade driver.")
})

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
message("\n", paste(rep("=", 60), collapse=""))
message("BDBV 2026 Spatiotemporal Analysis — COMPLETE")
message(paste(rep("=", 60), collapse=""))
message(sprintf("Analysis date:      %s", ANALYSIS_DATE))
message(sprintf("Training cutoff:    %s (week %d: %s to %s)",
                training_window_end, t_current, training_cutoff, training_window_end))
message(sprintf("Zones included:     %d", length(zones_all)))
message(sprintf("Models run:         %d frequentist (current forecast) + %d cross-validated (LFO)%s",
                if (nrow(fc_all_current) > 0 && "method" %in% names(fc_all_current))
                  dplyr::n_distinct(fc_all_current$method) else 0L,
                if (!is.null(lfo_results) && "method" %in% names(lfo_results))
                  dplyr::n_distinct(lfo_results$method) else 0L,
                if (!.have_freq) "  [Bayesian-only run]" else ""))
message(sprintf("Forecast rows:      %d", nrow(fc_all_current)))
message(sprintf("LFO-CV rows:        %s",
                if (!is.null(lfo_results)) nrow(lfo_results) else "skipped"))
message(sprintf("Outputs in:         %s", OUT_DIR))
message(paste(rep("=", 60), collapse=""))

# ---- Run-metadata info file (outputs/key_outputs/run_info.{json,md}) --------
# Save key facts about this run — line-list cutoff, analysis date, health zones
# invaded, CV fold count, forecast windows — from THIS run's in-memory values.
# Guarded so a failure here can never fail the pipeline.
tryCatch(source(file.path(ROOT, "write_run_info.R")),
         error = function(e) message("[run_info] skipped: ", conditionMessage(e)))

# ---- Model-selection provenance (outputs/key_outputs/model_selection.{json,md}) ----
# Detailed record of the FEATURED Bayesian model (+ best renewal + headline) chosen
# by the leave-future-out CV composite, with the full scored leaderboard and a
# cross-check that this file's independent scoring reproduces the pipeline's picks.
# Guarded so a failure here can never fail the pipeline.
tryCatch(source(file.path(ROOT, "write_model_selection.R")),
         error = function(e) message("[model_selection] skipped: ", conditionMessage(e)))

invisible(list(
  fc_all        = fc_all_current,
  fc_calibrated = calibrated_forecasts,
  fc_ensemble   = fc_ensemble,
  lfo_results   = lfo_results,
  metrics       = metrics_summary,
  calibration   = calibration_gain,
  mobility      = mobility_matrices,
  gt_pmfs       = gt_pmfs,
  Rt            = Rt_scalar,
  dat           = dat
))
