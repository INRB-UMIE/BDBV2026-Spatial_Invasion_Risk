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
#   12  Risk scores (no separate calibration step; see 16b recalibration)
#   13  Visualisations
#   14  Methodology report
#
# The list above stops at 14 because it predates most of what this script now schedules.
# What actually runs, after the steps above:
#   16b Prequential recalibration (delta per method x horizon); the pooled delta is DEPLOYED
#       into the live forecast, the per-fold one is used for scoring
#   18  Mean / median / loo-stacked ensembles
#   19  Space-time evaluation      20  Forecast detail, risk tables, published CSVs
#   21  Bayesian renewal suite (the featured model)      22  Daily re-issue
#   24  Review-response figures (FigR1/R3/R4)            27  Mobility comparison
#   Then write_run_info.R and write_model_selection.R (run metadata + selection provenance),
#   and .run_downstream(), which launches the figure suites (make_publication_figures.R,
#   make_manuscript_figures.R, make_si_model_figures.R, make_si_data_figures.R,
#   make_topk15_ever.R, make_spread_kinematics_compact.R),
#   update_bayesian_report.R, and the cascade driver run_cascade.R (modules 30-43).
#   Finally the key_outputs manifest and the retained-figure archive sweep.
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

# ONE environment-flag parser for the whole file. Switches here previously disagreed:
# RUN_EPINOWCAST_DIAGNOSTIC accepted only the literal "1" (so =TRUE silently did nothing), while
# MAKE_DELAY_FIGURES disabled only on "0" and .flag() in .run_downstream accepted 1/true/t/yes/y.
# Accept the same spellings everywhere, and treat an unset/blank value as the default.
#
# allow_global: when the env var is unset, ALSO honour a variable of the same name already in
# the session. RUN_RT_ESTIMATION documents exactly that workflow ("Set TRUE here - or
# export ..."), so it is deliberate there and NOT for the rest.
# It used to be implicit, via two hand-rolled copies of this parser that each assigned back to
# a GLOBAL of the same name: on a re-source in the same session (the usage this file's header
# documents) with the env var now unset, they picked up the PREVIOUS run's value instead of the
# documented default, and nothing said so. Making it a named argument keeps the workflow and
# removes the surprise — everything else is stateless across a re-source.
.env_flag <- function(nm, default = FALSE, allow_global = FALSE) {
  v <- tolower(trimws(Sys.getenv(nm, "")))
  if (nzchar(v)) {
    .yes <- c("1", "true", "t", "yes", "y")
    .no  <- c("0", "false", "f", "no", "n")
    if (v %in% .yes) return(TRUE)
    if (v %in% .no)  return(FALSE)
    # An UNRECOGNISED value is a typo, not a request to disable. Matching only the truthy
    # spellings meant `on`, `Y E S`, or a fat-fingered `ture` silently turned a default-TRUE
    # stage OFF — losing a whole suite with no message anywhere. Say so and keep the default.
    warning(sprintf(paste0("[flags] %s=%s is not a recognised boolean (use one of %s / %s); ",
                           "falling back to the default (%s)."),
                    nm, shQuote(Sys.getenv(nm, "")), paste(.yes, collapse = "/"),
                    paste(.no, collapse = "/"), isTRUE(default)),
            call. = FALSE, immediate. = TRUE)
    return(isTRUE(default))
  }
  if (isTRUE(allow_global)) return(isTRUE(get0(nm, ifnotfound = default)))
  isTRUE(default)
}

# Lightweight phase timing (diagnostic): minutes elapsed per major phase.
# A re-source() in the SAME session (the usage this file's header documents) leaves the previous
# run's objects behind, and several blocks below gate on exists(). In a Bayesian-only re-run after
# a frequentist run that would hand the report the PRIOR run's frequentist beta0 / covariate
# table. Clear the ones that are gated by exists() rather than by is.null().
suppressWarnings(rm(list = intersect(c("lfo_fig3", ".bl_ok", "lfo_results"),
                                     ls(envir = .GlobalEnv)), envir = .GlobalEnv))

.PH_T0 <- Sys.time(); .PH_LAST <- .PH_T0
.phase <- function(lbl) {
  now <- Sys.time()
  message(sprintf("[timing] %-42s %6.1f min  (cum %5.1f)", lbl,
                  as.numeric(difftime(now, .PH_LAST, units = "mins")),
                  as.numeric(difftime(now, .PH_T0, units = "mins"))))
  if (exists(".mark", mode = "function")) .mark(sprintf("PHASE     %s", lbl))
  .PH_LAST <<- now
}

# ---------------------------------------------------------------------------
# RUN MARKER — written IMMEDIATELY, so "did run_all actually start, and how far
# did it get?" is answerable from disk alone.
# ---------------------------------------------------------------------------
# run_info.json is written near the END of the modelling block, and the downstream
# logs are written only by run_one(). A run that dies before either leaves NO trace
# at all, and the outputs on disk are then indistinguishable from a run that was
# never launched -- which is exactly how a failed run was mistaken for a silent one.
# This marker is created at startup and stamped at each phase, so the last line
# always names the furthest point reached.
.RUN_MARKER <- tryCatch({
  .d <- file.path(OUT_DIR, "logs"); dir.create(.d, recursive = TRUE, showWarnings = FALSE)
  .f <- file.path(.d, "run_all_progress.log")
  cat(sprintf("%s  START     run_all.R  pid=%s  ANALYSIS_DATE=%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Sys.getpid(),
              format(get0("ANALYSIS_DATE", ifnotfound = NA))),
      file = .f, append = FALSE)
  .f
}, error = function(e) NULL)
.mark <- function(lbl) {
  if (is.null(.RUN_MARKER)) return(invisible(NULL))
  tryCatch(cat(sprintf("%s  %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), lbl),
               file = .RUN_MARKER, append = TRUE), error = function(e) NULL)
  invisible(NULL)
}
# A run that dies leaves no DONE line, so the absence of one IS the signal. There is no
# finalizer: reg.finalizer() on the global environment does not fire reliably on a crash or
# a kill, and a marker that lies about having run is worse than no marker.

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
  # PER-CLASSIFICATION delays (2026-09-22). The onset imputation imputes onsets for CONFIRMED
  # records, so it must draw from the confirmed-case delay; the pooled fit above is majority
  # test-negative and one-seventh unadjudicated on this line list (6,823 not_a_case windowed
  # pairs and 2,038 with no final classification, against 4,981 confirmed) and runs 2.47 d
  # (32%) short. Each stratum is fitted independently, with its own interval censoring and
  # its own truncation correction, because confirmed cases are sampled more slowly and are
  # therefore also more truncated -- a pooled correction would distort the contrast.
  # The pooled rows are still written unchanged, so every other consumer is unaffected.
  .osstr   <- tryCatch(estimate_onset_sample_strata(.dll_raw, analysis_date = ANALYSIS_DATE),
                       error = function(e) {
                         warning("[run_all] delay strata not fitted (", conditionMessage(e),
                                 "); the onset imputation will warn and use the POOLED delay.",
                                 call. = FALSE); list() })
  if (!is.null(.osfit)) {
    write_onset_sample_long(.osfit,
      file.path(DATA_DIR, "cfr_reference", "dhis2_onset_sample_delay_params.csv"),
      source_label = .dmeta$folder, epidist = .osepi, strata = .osstr)
    write_onset_sample_long(.osfit,
      file.path(LINELIST_DIR, .dmeta$folder, "dhis2_onset_sample_delay_params.csv"),
      source_label = .dmeta$folder, epidist = .osepi, strata = .osstr)
    if (length(.osstr))
      message(sprintf("[run_all] delay strata written: %s",
                      paste(sprintf("%s %.2f d (n=%s)", names(.osstr),
                                    vapply(.osstr, function(e) e$mean, 0),
                                    vapply(.osstr, function(e) as.character(e$n), "")),
                            collapse = "; ")))
    message(sprintf("[run_all] DHIS2 onset->sample delay: %s, mean %.2f d, rate %.3f/d (n=%d, window %s)%s",
                    .osfit$best_family, .osfit$implied_mean, .osfit$rate, .osfit$n_fit, .osfit$window,
                    if (!is.null(.osepi))
                      sprintf(" | EpiDist marginal PREFERRED: %s mean %.2f d sd %.2f d",
                              .osepi$family, .osepi$mean, .osepi$sd)
                    else " | EpiDist marginal: not used (censored-MLE)"))
  } else message("[run_all] DHIS2 delay fit returned NULL (too few pairs); imputation uses the windowed bootstrap.")
}, error = function(e) message("[run_all] DHIS2 delay refresh skipped (non-fatal): ", conditionMessage(e)))

# 1b. VISUAL CHECK of the delay fits just made. The onset->sample distribution feeds the
#     onset imputation, the nowcast rate and the R(t) truncation, so a bad fit propagates
#     into every downstream number — and until now nothing in a run ever plotted it (04c's
#     QA figure sits inside its Rscript-gated MAIN block, which source()ing skips). This
#     renders the EpiDist naive + marginal densities and the interval-censored MLE against
#     the empirical delay histograms, on exactly the records the fits used.
#     Free by default: the onset->sample EpiDist posterior is reused from step 1a via 04c's
#     draw registry (no Stan refit); only the cheap censored MLE is re-run for the other
#     three delays. DELAY_FIG_EPIDIST_ALL=TRUE additionally fits EpiDist for those three
#     (12 more Stan fits, minutes). MAKE_DELAY_FIGURES=0 skips the step entirely.
#     Non-fatal throughout — a missing package or a failed panel never stops the run.
if (.env_flag("MAKE_DELAY_FIGURES", TRUE)) {
  message("\n=== Step 1b: Delay-fit diagnostic figures ===")
  tryCatch({
    source(file.path(ST_DIR, "04d_delay_fit_figures.R"))   # make_delay_fit_figures()
    if (exists(".dll_raw", inherits = TRUE))
      make_delay_fit_figures(
        .dll_raw, analysis_date = ANALYSIS_DATE,
        os_fit      = if (exists(".osfit", inherits = TRUE)) .osfit else NULL,
        epidist_tbl = if (exists(".osepi", inherits = TRUE)) attr(.osepi, "epidist_table") else NULL)
    else message("[run_all] delay figures skipped: step 1a did not load the line list.")
  }, error = function(e)
    message("[run_all] delay-fit figures skipped (non-fatal): ", conditionMessage(e)))
}
.phase("Step 1b  delay-fit figures")

message("\n=== Step 1: Data preparation ===")
source(file.path(ST_DIR, "01_data_prep.R"))
dat <- prep_all_data()

zone_week_nc <- dat$zone_week   # all zone-week counts (will be nowcast-corrected in step 4)
zones_all    <- names(dat$pop)  # canonical zone name vector (519 zones, from population spine)
pop_vec      <- dat$pop
covariates   <- dat$covariates
sitrep       <- dat$sitrep

# ---------------------------------------------------------------------------
# 2. Epidemiological parameters
# ---------------------------------------------------------------------------
.phase("Step 1   data prep")
message("\n=== Step 2: Epidemiological parameters ===")
source(file.path(ST_DIR, "02_epi_params.R"))

# GT PMFs for all profiles
gt_pmfs <- compute_all_gt_pmfs()           # named list: short/medium/long → PMF vector

# EpiNow2 R(t) estimation (national-level, per GT profile). DIAGNOSTIC OUTPUT ONLY — this
# panel does NOT feed any forecast. The Bayesian invasion suite fits its OWN EpiNow2 R draws
# for the h >= 2 source projection (bayes_rt_week_draws() in 21_bayesian_renewal.R: the same
# model via .epinow2_rt_fit(), but with the invasion model's generation-time parameters and
# as-of censoring).
#
# DEFAULT TRUE since the 2026-09-17 streamlining: rt_national.pdf is a RETAINED
# deliverable, and it is the sole product of this block — with the flag off,
# Rt_primary stayed NULL and plot_rt() silently no-opped, so the figure could never
# be produced by a default run. It is a separate MCMC fit per GT profile (three fits)
# and one of the slower steps; export RUN_RT_ESTIMATION=FALSE to skip it when
# iterating on downstream figures.
RUN_RT_ESTIMATION <- .env_flag("RUN_RT_ESTIMATION", TRUE, allow_global = TRUE)
Rt_national_list <- list()
Rt_primary <- NULL
Rt_scalar  <- NULL
Rt_zone    <- NULL   # zone-level R(t): optional, never populated in the default path
if (isTRUE(RUN_RT_ESTIMATION)) {
  # CONFIRMED cases only (dat$ll retains all classifications, so an unfiltered count would be
  # "all alerts", inflating R(t)), dated EXACTLY as the weekly case counts (date_index: onset,
  # or an onset imputed from the sample date) by epinow2_daily_confirmed() — the same series
  # the Bayesian suite's R draws use. Dating onset-less cases at their sample date (the former
  # convention) let the right-truncation model inflate already-complete cases and biased the
  # recent R upward; see epinow2_daily_confirmed() in 02_epi_params.R.
  daily_cases <- epinow2_daily_confirmed(dat$ll, ANALYSIS_DATE, caller = "Rt")
  for (gt_name in names(gt_pmfs)) {
    message(sprintf("[Rt] Estimating national R(t) with GT profile: %s", gt_name))
    gt_p <- GT_PROFILES[[gt_name]]
    Rt_national_list[[gt_name]] <- tryCatch(
      # ll is passed for the truncation BASIS GUARD only (it compares the recorded-onset
      # basis used for the truncation fit against the date_index basis this series uses).
      # issue_date stays NULL: this is the deployed fit, so the extract regime applies.
      estimate_rt_epinow2(daily_cases, gt_name, gt_p, ANALYSIS_DATE, ll = dat$ll),
      error = function(e) { warning("EpiNow2 failed for ", gt_name, ": ", e$message); NULL }
    )
  }
  # Use primary GT profile R(t)
  Rt_primary <- Rt_national_list[[GT_PRIMARY]]
  if (!is.null(Rt_primary)) {
    # ESTIMATES ONLY. EpiNow2 projects 7 days past the data by default, so the last row of
    # this table is a forecast, not a measurement. Reporting it as "the current R(t)" — and
    # deriving Rt_scalar$sd from its (wider) forecast band — presented a 7-day-ahead
    # projection as an independent estimate. Filter on `type` when it is present (it is,
    # from RT_CACHE_VERSION 6), else fall back to the as-of date.
    .rt_est <- if ("type" %in% names(Rt_primary))
      Rt_primary[!(Rt_primary$type %in% "forecast"), , drop = FALSE] else
      Rt_primary[as.Date(Rt_primary$date) <= ANALYSIS_DATE, , drop = FALSE]
    if (nrow(.rt_est) == 0) .rt_est <- Rt_primary   # degenerate; better than an empty tail()
    Rt_scalar <- list(
      mean = tail(.rt_est$R_mean, 1),
      sd   = (tail(.rt_est$R_hi_90, 1) - tail(.rt_est$R_lo_90, 1)) / (2 * 1.645)
    )
    message(sprintf("[Rt] Primary R(t) estimate (as of %s, forecast tail excluded): %.2f (90%% CI: %.2f-%.2f)",
                    format(tail(.rt_est$date, 1)), Rt_scalar$mean,
                    tail(.rt_est$R_lo_90, 1), tail(.rt_est$R_hi_90, 1)))
  } else {
    Rt_scalar <- list(mean = 1.5, sd = 0.3)  # fallback prior
    warning("[Rt] EpiNow2 failed; using prior R(t) = 1.5")
  }
} else {
  message("[Rt] RUN_RT_ESTIMATION = FALSE — skipping EpiNow2 R(t); rt_national.pdf will NOT be produced.")
}

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

# OSRM road-distance (-dist) kernels are OPTIONAL (INCLUDE_OSRM_DIST_MODELS, default TRUE).
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

# Primary mobility matrix for main analysis. MOBILITY_PRIMARY is the FILLED composite
# (00_config.R, now the COHORT composite M14-fill).
#
# THE FALLBACK CHAIN MUST END ON A KERNEL THAT ALWAYS EXISTS. Stripping "-fill" alone is not
# enough any more: M14 and M14-fill are BOTH built only inside `if (want_cohort)` and only
# when build_M_cohort() returned non-NULL (03_mobility_matrices.R), so switching
# INCLUDE_COHORT_MODELS off, emptying COHORT_SOURCES, or any cohort build failure would leave
# W_primary NULL and kill the run at the stopifnot below. Under the previous M8-fill primary
# the parent M8 was built unconditionally, so this hard-fail path is new. M8-fill/M8 are the
# last resort because M8 is built on every path.
.primary_chain <- unique(c(MOBILITY_PRIMARY, sub("-fill$", "", MOBILITY_PRIMARY),
                           "M8-fill", "M8"))
MOBILITY_PRIMARY_BUILT <- .primary_chain[
  vapply(.primary_chain, function(k) !is.null(mobility_matrices[[k]]), logical(1))][1]
if (is.na(MOBILITY_PRIMARY_BUILT))
  stop("[mobility] no primary kernel could be resolved; none of ",
       paste(.primary_chain, collapse = ", "), " was built. The mobility layer is broken — ",
       "check build_all_mobility_matrices()'s warnings above.", call. = FALSE)
if (!identical(MOBILITY_PRIMARY_BUILT, MOBILITY_PRIMARY))
  warning(sprintf(paste0("[mobility] %s was not built; the primary kernel falls back to %s. ",
                         "W_primary drives the M11 inward-FOI kernel and the mobility figure, ",
                         "so this is a substantive change of kernel, not a cosmetic one."),
                  MOBILITY_PRIMARY, MOBILITY_PRIMARY_BUILT), call. = FALSE)
W_primary <- mobility_matrices[[MOBILITY_PRIMARY_BUILT]]
stopifnot(!is.null(W_primary))

# ---------------------------------------------------------------------------
# 4. Nowcast correction
# ---------------------------------------------------------------------------
.phase("Step 3   mobility matrices")
message("\n=== Step 4: Nowcast correction ===")
source(file.path(ST_DIR, "04_nowcasting.R"))
source(file.path(ST_DIR, "04b_epinowcast.R"))

# TRAIN/DEPLOY CONSISTENCY (2026-09-17). The deployed training data and the LFO folds
# must be nowcast by the SAME estimator, because the recalibration delta (16b,
# INVASION_RECALIBRATE_DEPLOY) is estimated on the folds and applied to the deployed
# forecast. It previously was not: the folds used apply_nowcast_correction() (effective
# h=1 inflation x1.323) while the deployed path used epinowcast, which — through the
# base-mismatch bug since fixed in 04b — delivered x1.027. Delta was therefore fitted in
# one regime and applied in another, biasing published current-week invasion probabilities
# LOW by roughly 1.29x on the hazard scale. The cascade inherited the same mismatch.
#
# apply_nowcast_correction() is now the estimator on BOTH sides. It is the right one to
# standardise on: it scales exactly the quantity it models (onset-dated confirmed counts
# observed by the as-of date), which is the same quantity in a fold and at deployment, so
# one estimator applies identically to both. epinowcast targets a different population
# (the reporting triangle: only cases carrying BOTH an onset and a sample date, 74% of
# confirmed on the 2026-09-07 snapshot) and needs a reconciliation step to reach the
# zone-week base — which is where the bug lived.
#
# epinowcast is RETAINED as the sensitivity analysis (bayes_nowcast_sensitivity, step 13),
# which is its proper home and which — now that 04b is fixed — is finally informative
# rather than bit-identical to the raw arm.
# THE NOWCAST USES THE SAME ESTIMATED TRUNCATION AS R(t), not the onset->sample delay.
# This is the h = 1 channel: R(t) enters only at h >= 2 (.bayes_project_mu), so the 1-week
# invasion hazard is driven entirely by these corrected counts. They were corrected for
# onset->SAMPLE (~7.7 d) while the series is truncated by onset->APPEARANCE (~11.4 d measured),
# so the most recent week was inflated 2.63x where the fitted truncation gives 5.92x —
# under-correcting the very week the live forecast is issued from. Resolved once here (the fit is cached) and passed
# explicitly, so the deployed nowcast and the deployed R(t) provably use one distribution.
# NULL falls through to apply_nowcast_correction()'s own resolver, which warns.
.trunc_deployed <- tryCatch(
  .trunc_as_delay_spec(epinow2_truncation(issue_date = NULL, ll = dat$ll), "extract"),
  error = function(e) { warning("[nowcast] deployed truncation unavailable (", conditionMessage(e),
                                "); falling back to the onset->sample delay, which UNDER-corrects ",
                                "the recent weeks.", call. = FALSE); NULL })
if (!is.null(.trunc_deployed))
  message(sprintf("[nowcast] deployed truncation: %s (mean %.2f d, sd %.2f d, %s)",
                  .trunc_deployed$family, .trunc_deployed$mean, .trunc_deployed$sd,
                  .trunc_deployed$estimator))
zone_week_nc <- apply_nowcast_correction(zone_week_nc, analysis_date = ANALYSIS_DATE,
                                         delay = .trunc_deployed)
message(sprintf("[nowcast] Method: %s (train/deploy-consistent); mean weight: %.3f",
                attr(zone_week_nc, "nowcast_method") %||% "deterministic",
                mean(zone_week_nc$trunc_weight, na.rm = TRUE)))

# Weekly epinowcast factors remain a reported DIAGNOSTIC (they are what the sensitivity
# arm quantifies), but they no longer drive the training data. Fit only when asked:
# RUN_EPINOWCAST_DIAGNOSTIC=1. Non-fatal — a missing cmdstanr must never stop a run.
if (.env_flag("RUN_EPINOWCAST_DIAGNOSTIC", FALSE)) {
  tryCatch({
    .enw <- nowcast_zone_week_epinowcast(
      zone_week = dat$zone_week, linelist = dat$ll,
      analysis_date = ANALYSIS_DATE, outbreak_start = OUTBREAK_START)
    enw_factors <- attr(.enw, "epinowcast_factors")
    if (!is.null(enw_factors))
      readr::write_csv(enw_factors, file.path(OUT_DIAGNOSTICS, "epinowcast_weekly_factors.csv"))
  }, error = function(e)
    message("[nowcast] epinowcast diagnostic skipped (non-fatal): ", conditionMessage(e)))
}

# ---------------------------------------------------------------------------
# 5–9. Model fitting and current-week forecasts
# ---------------------------------------------------------------------------
source(file.path(ST_DIR, "05_baseline_models.R"))   # forecast_B1/B4, zone_week_to_wide
source(file.path(ST_DIR, "06_simple_models.R"))     # compute_foi, daily_to_weekly_gt
source(file.path(ST_DIR, "15_workhorse.R"))         # mobility-informed renewal workhorse
source(file.path(ST_DIR, "16_invasion_eval.R"))     # invasion LFO + evaluation
source(file.path(ST_DIR, "16b_invasion_recalibration.R"))  # post-hoc hazard-scale delta
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
# OPTIONAL family (INCLUDE_M11_MODELS, default FALSE since 2026-09-21): NOT built unless the
# toggle is turned on, so the M11 renewal/Bayesian variants are NOT part of the default grid
# (they gate on this matrix). It was turned off because M11 and M11-dist are the only grid
# kernels never written to disk — they are built in memory right here, so an M11 win would be
# silently discarded by .cascade_selected_kernel(), which only adopts a kernel whose .rds
# exists. See 00_config.R:1128-1140 for the full reasoning and what restoring it requires.
if (isTRUE(INCLUDE_M11_MODELS)) {
  # Built on the PRIMARY (filled) composite: the presence matrix P must not inherit the
  # unfilled kernel's asserted-zero destinations.
  mobility_matrices$M11 <- tryCatch(
    build_inward_contact_matrix(W_primary, pop_vec, zones_all),
    error = function(e) { warning("inward contact matrix (M11): ", e$message); NULL })
  if (!is.null(mobility_matrices$M11))
    message("[mobility] Built M11 inward/meeting-location effective-contact matrix")
  # M11-dist: the same inward/meeting-location FOI but built on the ROAD-DISTANCE composite,
  # so its presence matrix inherits the km-deterrence routing. Filled form where available,
  # for the same reason as M11. Only when that matrix was built.
  .m8d <- if (!is.null(mobility_matrices[["M8-dist-fill"]])) "M8-dist-fill" else "M8-dist"
  if (!is.null(mobility_matrices[[.m8d]])) {
    mobility_matrices[["M11-dist"]] <- tryCatch(
      build_inward_contact_matrix(mobility_matrices[[.m8d]], pop_vec, zones_all),
      error = function(e) { warning("inward contact matrix (M11-dist): ", e$message); NULL })
    if (!is.null(mobility_matrices[["M11-dist"]]))
      message("[mobility] Built M11-dist (inward FOI on the road-distance composite)")
  }
} else {
  message("[mobility] INCLUDE_M11_MODELS = FALSE — skipping the M11 inward-contact kernel(s).")
}

# Mobility-kernel comparison figure (kernel similarity, origin coverage, source-zone
# outflow profiles, and composite divergence) — a cheap post-build diagnostic.
# It is handed THIS RUN's kernels: reading the output directory instead would mix in
# mobility_*.rds files left by earlier runs (a kernel since switched off, or one from
# another branch) and show them as current. tryCatch so a plotting hiccup never aborts.
MOBILITY_MATRICES_CURRENT <- mobility_matrices
tryCatch({
  source(file.path(ST_DIR, "27_mobility_comparison.R"))
  make_mobility_comparison_figures(kernels = MOBILITY_MATRICES_CURRENT)
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
# horizon, mu_forecast, p_invasion, method).
# Primary: the mobility-informed renewal model (variants over GT / mobility /
# observation, PLUS the additional model options requested — covariate-augmented
# betas (Q2), extra mobility kernels M4b/M9/M10 (Q4), raw vs nowcast-corrected and
# reporting-rate structures (Q6), completeness-weighted calibration (Q1), and the
# NEW suspected-but-not-confirmed leading-indicator covariate models). Every
# addition is an EXTRA variant; the base formulations are untouched. OPTIONAL families
# are gated by the 00_config.R toggles. Defaults, verified against 00_config.R: the reduced
# "geo" covariate model (INCLUDE_GEO_COV_MODELS) is ON, the M11 inward-FOI variants
# (INCLUDE_M11_MODELS) are ON, the suspected-covariate models (INCLUDE_SUSPECTED_COV_MODELS)
# are OFF, and the FULL-exogenous covariate model (INCLUDE_FULL_COV_MODELS) is OFF. This
# paragraph previously stated all four the wrong way round. Comparators: gravity (B4) / distance (B1)
# / adjacency (B7) baselines. Broken models (S2/S3/S4, KNN B2, null
# B3/B5/B6, INLA C4, ZINB C5) remain removed. The full grid is cross-validated; only the
# best per family is featured in figures.
# THE FREQUENTIST RENEWAL FAMILY WAS REMOVED (2026-09-19).
# It was 28 model definitions (Renewal-M*) driven by forecast_workhorse(), behind a
# RUN_FREQUENTIST_MODELS flag that defaulted FALSE. Both are gone. Verified on the shipped
# evaluation table:
# ZERO Renewal-* models are scored — the published field is 52 Bayesian models plus the three
# always-on structural baselines (Gravity-B4, Distance-B1, Adjacency-B7), which are defined
# below and are unaffected. The workhorse and its private helpers (fit_import_beta,
# fit_import_model, estimate_R_local, .reporting_rate_vec, .beta_vector, .invasion_prob,
# .qcount, .workhorse_temporal_dispersion) are gone from 15_workhorse.R with them; the
# helpers the BAYESIAN design builder shares with it (.count_wide, .feature_matrix,
# .static_features, .susp_wide, .dmin_vec, .gweighted_own, compute_risk_scores) remain.
#
# The reporting-rate structure went too: `report_rate_vec` (health-site density as a
# completeness proxy) existed solely to feed Renewal-M8-report, and ascertainment has been
# removed from the pipeline entirely.
# The frequentist renewal arm is GONE, not merely switched off. Nothing defines
# INVASION_MODELS any more, so the current-week frequentist forecast, its ensemble, its
# risk tables and its map/parameter suites are all removed below rather than left behind
# guards that can never open. The LFO cache stamp still asks for `INVASION_MODELS` through
# get0(..., ifnotfound = character(0)) — the value it already had — so existing caches stay
# valid. BASELINE_MODELS (gravity-B4, distance-B1, adjacency-B7) are unaffected: they are
# structural comparators, not frequentist renewal models, and are always scored.

# ── Always-on BASELINE comparators (review §0.1 / §3.1 / §3.2) ───────────────
# The featured model is ALWAYS scored against structural baselines — the reviewer's central
# evaluation ask: source-mass x distance-decay import pressure (B4 — NOT a gravity model; it has no
# destination-mass term, see forecast_B4's docstring), inverse-distance (B1) and the
# NEAREST-AFFECTED / adjacency
# spatial-spread null (B7 — "does the virus just go next door?"). These are yardsticks to
# beat, NOT the removed frequentist renewal family. They are scored on the identical LFO
# folds as the Bayesian models (below).
#
# The hhh4 and stochastic-SEIR mechanistic comparators were removed with 07_hhh4_model.R /
# 08_stochastic_seir.R (2026-09-17 streamlining): both were gated OFF by default and neither
# fed any retained figure. Gravity-B4 and Adjacency-B7 must stay — make_manuscript_figures.R
# hard-codes them as BASELINE_METHODS for manuscript Figure 2.
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
  }
))
message(sprintf("[baselines] %d always-on baseline comparators: %s",
                length(BASELINE_MODELS), paste(names(BASELINE_MODELS), collapse = ", ")))

# ── Current-week forecasts (all models) ─────────────────────────────────────
affected_now <- affected_zones(zone_week_nc, training_cutoff)
message(sprintf("[invasion] %d/%d zones already affected; %d at-risk",
                length(affected_now), length(zones_all),
                length(zones_all) - length(affected_now)))
# The BAYESIAN suite (section 12) produces the current-week forecast; there is no
# frequentist current-week forecast to build, no member ensemble to pool over it and no
# fc_all_current.rds to write.

# ---------------------------------------------------------------------------
# 10. Invasion LFO-CV — same principled folds for ALL models
# ---------------------------------------------------------------------------
.phase("Step 4-9 nowcast (deterministic) + current fc")
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
# Bayesian grid that bayes_default_grid() composes from the CORE (mobility {M4, M8, M10} at the
# single medium GT anchor — the generation time is treated as KNOWN, neither swept across the
# grid nor marginalised over GT_PRIOR (00_config.R) — plus the
# cohort M13/M14, the OD/consensus M16/M17 and their source-cell-fill and origin-split forms)
# and whichever OPTIONAL families are toggled on in 00_config.R (M9, M11, M15, the FULL-exogenous
# and suspected-case covariate sets, the logit-link sensitivity, the time-varying-beta families,
# and — ON by default — the OSRM road-distance kernels and the reduced "geo" covariate set) — and the
# FEATURED Bayesian model is chosen BY CV SKILL (best_bayes_method), not hardcoded. Each model is
# refit per fold by MCMC, so LFO uses lighter sampling (iter=600) to keep the full-
# grid cross-validation tractable; set RUN_BAYES_LFO <- FALSE to skip it entirely.
# The AS-OF truncation for the per-fold nowcasts, resolved ONCE before the LFO.
# regime = "asof" is forced rather than derived from a date, because deriving it here would
# pick the EXTRACT regime: issue_date is absent, and absent means deployed. The fit itself is
# shared across folds by construction (epinow2_truncation() always evaluates the as-of regime
# at ANALYSIS_DATE), matching run_invasion_lfo()'s SHARED NUISANCE PARAMETER policy.
# Resolving it HERE, before the parallel LFO, also means every worker sees one finished
# estimate rather than recomputing it.
# Must sit OUTSIDE the lfo_results tryCatch(): an assignment there is not an argument.
# It is also stamped into .lfo_stamp below, so REUSE_CACHED_LFO cannot serve an LFO built
# under a different truncation.
.trunc_asof <- tryCatch(
  .trunc_as_delay_spec(epinow2_truncation(regime = "asof", ll = dat$ll), "asof"),
  error = function(e) { warning("[nowcast] as-of truncation unavailable (", conditionMessage(e),
                                "); folds fall back to the onset->sample delay.",
                                call. = FALSE); NULL })
if (!is.null(.trunc_asof))
  message(sprintf("[nowcast] fold truncation: %s (mean %.2f d, sd %.2f d, %s)",
                  .trunc_asof$family, .trunc_asof$mean, .trunc_asof$sd, .trunc_asof$estimator))

BAYES_LFO_SPECS <- bayes_default_grid(mobility_matrices)
bayes_lfo_models <- list()
if (isTRUE(RUN_BAYES_LFO) && requireNamespace("brms", quietly = TRUE)) {
  for (sp in BAYES_LFO_SPECS) {
    bayes_lfo_models[[sp$label]] <- make_bayes_lfo_model(sp$mob, sp$gt, sp$cov,
      mobility_matrices, gt_pmfs, covariates, osrm_mat, zones_all,
      iter = 600L, link = sp$link %||% "cloglog", tv = sp$tv %||% "none", linelist = dat$ll,
      # ROLLING AS-OF PREDICTORS: each training transition's import force is computed from the
      # counts known AT that transition, not from one matrix per fold. Without this, beta0 is
      # fitted against a lambda ~2x smaller than the forecast's and the recalibration delta
      # absorbs the difference (see ROLLING_PREDICTORS in 00_config.R).
      trunc_delay = .trunc_asof)
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
# The cache is stamped with the Bayesian h >= 2 source projection and the Bayesian grid it was
# produced under. A cache from a different projection (e.g. the retired crude-R / fixed-constant
# one) or grid (e.g. a different INCLUDE_UNFILLED_MODELS) is NOT reused: its selection and
# recalibration factors would not describe the models this run forecasts with.
# RT_CACHE_VERSION is part of the stamp because the h >= 2 projection's R comes from EpiNow2: any
# change to how that R is estimated (e.g. v5, the date_index dating) changes the backtest too.
# The stamp must cover EVERYTHING that changes the LFO, and it did not: the DATA SNAPSHOT and
# ANALYSIS_DATE were absent, so a cache from an older data cut was accepted wholesale — and the
# recalibration table it produces is DEPLOYED into the live forecast below, i.e. a delta fitted
# on one data cut would be applied to another. Also absent were LEAKAGE_FREE_LFO, LFO_HORIZONS,
# and the baseline/invasion model sets (so a change to the scored field was invisible).
.lfo_stamp <- list(bayes_projection = BAYES_H2_PROJECTION,
                   rt_cache_version = RT_CACHE_VERSION,
                   analysis_date    = format(ANALYSIS_DATE),
                   linelist_snapshot = tryCatch({
                     .lj <- get0("LINELIST_JSON",
                                 ifnotfound = file.path(get0("LINELIST_DIR", ifnotfound = ""),
                                                        "latest.json"))
                     .m <- jsonlite::fromJSON(.lj, simplifyVector = TRUE)
                     paste(.m$folder, .m$processed_at, sep = "|")
                   }, error = function(e) NA_character_),
                   leakage_free     = isTRUE(get0("LEAKAGE_FREE_LFO", ifnotfound = TRUE)),
                   lfo_horizons     = as.integer(LFO_HORIZONS),
                   # THE FOLD WINDOW. These decide WHICH rounds exist and therefore which
                   # outcomes the recalibration delta is fitted on and which rows every
                   # metric is pooled over. A cache built under the old 25-day admission
                   # threshold holds five fewer rounds and must not be served for this run.
                   lfo_min_eval_age  = as.integer(get0("LFO_MIN_EVAL_AGE_DAYS", ifnotfound = 25L)),
                   lfo_min_events    = 0L,
                   # STAMP WHAT THE BASELINES ACTUALLY COMPUTE, not just their wrappers.
                   # as.character(BASELINE_MODELS) deparses the three closures DEFINED IN THIS
                   # FILE, so a change to forecast_B1/B4/B7 in 05_baseline_models.R — the B1
                   # cumulative-hazard rewrite, the prob_calibrated flags — left the stamp
                   # byte-identical and REUSE_CACHED_LFO=1 would have accepted a pre-fix cache.
                   baseline_models  = sort(as.character(get0("BASELINE_MODELS",
                                                             ifnotfound = character(0))),
                                           method = "radix"),
                   baseline_bodies  = vapply(c("forecast_B1", "forecast_B4",
                                               "forecast_B7_adjacency"),
                                             function(f) if (exists(f, mode = "function"))
                                               paste(deparse(body(get(f))), collapse = "") else "",
                                             character(1)),
                   # M11 is built FROM W_primary, so the primary kernel changes what
                   # `Bayes-M11-*` denotes without changing its label. Without this, a cached
                   # LFO from a different MOBILITY_PRIMARY is accepted and the ensemble and
                   # selection silently inherit a different matrix under the same name.
                   mobility_primary = as.character(get0("MOBILITY_PRIMARY",
                                                        ifnotfound = NA_character_)),
                   invasion_models  = sort(as.character(get0("INVASION_MODELS",
                                                             ifnotfound = character(0))),
                                           method = "radix"),
                   # ROLLING AS-OF PREDICTORS and their floor. Changing either changes the
                   # import force of every training row, so a cached LFO built under the
                   # other setting describes different models entirely. The floor is stamped
                   # as well as the flag because it decides how many folds survive (3 -> 11
                   # of 12, 5 -> 9) and therefore which outcomes the delta is fitted on.
                   rolling_predictors = isTRUE(get0("ROLLING_PREDICTORS", ifnotfound = TRUE)),
                   rolling_floor      = as.integer(get0("ROLLING_PREDICTOR_FLOOR",
                                                        ifnotfound = 3L)),
                   # THE FOLD TRUNCATION. The R(t) FILE cache hashes it, but this RESULT
                   # cache did not: changing TRUNC_WINDOW_DAYS / TRUNC_THIN_DAYS /
                   # TRUNC_MAX_DAYS / TRUNC_MIN_SNAPSHOTS / TRUNC_SERIES_VERSION — or simply
                   # gaining line-list rows that move the as-of estimate — without also bumping
                   # RT_CACHE_VERSION left this stamp byte-identical, so REUSE_CACHED_LFO=1
                   # would serve an LFO produced under a different fold nowcast AND a different
                   # fold truncation. Hashing the FITTED spec rather than the config is tighter:
                   # a config edit that does not move the fit does not invalidate the cache.
                   trunc_asof = if (is.null(.trunc_asof)) "none" else
                     rlang::hash(.trunc_asof[c("family", "params", "mean", "sd", "estimator")]),
                   bayes_grid = sort(vapply(BAYES_LFO_SPECS, function(s) s$label, character(1)),
                                     method = "radix"))   # locale-independent order
if (isTRUE(.reuse_lfo)) {
  .cached_stamp <- tryCatch(attr(readRDS(.lfo_cache), "lfo_stamp"), error = function(e) NULL)
  if (!identical(.cached_stamp, .lfo_stamp)) {
    warning("[LFO] REUSE_CACHED_LFO was requested, but the cached cross-validation was produced ",
            "under a different data snapshot, analysis date, horizon set, model grid or source ",
            "projection; re-running it.", call. = FALSE)
    .reuse_lfo <- FALSE
  }
}
lfo_results <- if (isTRUE(.reuse_lfo)) {
  message("[LFO] REUSE_CACHED_LFO=TRUE — loading cached lfo_cv_results.rds (cross-validation NOT re-run)")
  tryCatch(readRDS(.lfo_cache), error = function(e) { warning("cached LFO load failed: ", e$message); NULL })
} else tryCatch(
  # BASELINE_MODELS are always scored (review §3.1/§3.2), alongside the Bayesian models.
  # delay = the AS-OF truncation. Each fold's training counts come from reaggregate_asof(),
  # which censors on linelist_observation_date(), so those series are truncated by
  # onset->SAMPLE — a shorter process than the deployed extract's onset->appearance. Handing
  # the folds the deployed truncation would over-correct them and break `delta`, which is
  # fitted here and applied to the live forecast. run_invasion_lfo() threads this to every
  # per-fold nowcast_fn() call. One shared delay across folds, exactly as before.
  run_invasion_lfo(zone_week_outbreak, c(BASELINE_MODELS, bayes_lfo_models),
                   delay = .trunc_asof,
                   horizons = LFO_HORIZONS, analysis_date = ANALYSIS_DATE,
                   # SCORE EVERY COMPLETE ROUND, INCLUDING THE QUIET ONES. Dropping rounds
                   # with no invasion is selection on the outcome: it enriches the pooled
                   # base rate (flattering calibration-in-the-large) and punches a hole in
                   # every over-rounds panel at exactly the rounds where nothing happened.
                   # Per-fold discrimination is undefined at an event-free round and
                   # returns NA, which is the correct gap. On the 2026-09-07 frame no
                   # ADMITTED round is event-free (the thinnest carry one invasion at h=1),
                   # so 0 and 1 enumerate the same grid here; the setting is what keeps a
                   # future quiet round in rather than selecting it away.
                   min_atrisk_events = 0L,
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
  # Bayesian ENSEMBLE in the SAME folds — the pre-specified, leakage-free analogue of the
  # frequentist ensemble: mean/median over the diverse mobility-kernel Bayesian members
  # (M4/M8/M9/M10/M11 at medium GT, cloglog). Labelled "Bayes-ens-{mean,median}" so it is
  # picked up by the Bayesian-restricted discrimination / detection / top-K figures and can
  # be compared against the individual Bayesian models. (The loo-STACKED ensemble remains the
  # separate full-data featured product — "Bayes-stack".)
  # "Bayes-M11-inward" REMOVED 2026-09-21 with INCLUDE_M11_MODELS = FALSE. It was already
  # protected by the intersect() below, so its absence would not have errored — it would have
  # silently shrunk the ensemble from 4 members to 3. Listing a model that cannot exist is
  # exactly the kind of quiet drift this file guards against elsewhere.
  #
  # NOTE, pre-existing and NOT changed here: this list spans M4/M8/M9/M10 but contains no
  # member of the COHORT composite family (M13/M14/M16/M17) — even though MOBILITY_PRIMARY is
  # M14-fill. With M9 and M11 both off, the realised ensemble is three members: Bayes-M4-med
  # plus the fill twins of M8 and M10. Whether the primary family should be represented is a
  # modelling decision, flagged rather than taken.
  BAYES_ENSEMBLE_MEMBERS <- c("Bayes-M4-med", "Bayes-M8-med", "Bayes-M9-med", "Bayes-M10-med")
  # INCLUDE_UNFILLED_MODELS = FALSE does not fit an unfilled composite whose fill twin was built,
  # so take that member's fill twin: the ensemble keeps the same kernel families instead of
  # shrinking. Only a twin that was actually cross-validated is swapped in — when the fill
  # kernels were not built the grid keeps the unfilled model, and so does the ensemble.
  if (!isTRUE(INCLUDE_UNFILLED_MODELS)) {
    .twin <- bayes_fill_twin_label(BAYES_ENSEMBLE_MEMBERS)
    .swap <- !is.na(.twin) & .twin %in% unique(lfo_results$method)
    BAYES_ENSEMBLE_MEMBERS[.swap] <- .twin[.swap]
  }
  BAYES_ENSEMBLE_MEMBERS <- intersect(BAYES_ENSEMBLE_MEMBERS, unique(lfo_results$method))
  if (length(BAYES_ENSEMBLE_MEMBERS) >= 2L)
    lfo_results <- append_ensembles(lfo_results, BAYES_ENSEMBLE_MEMBERS, prefix = "Bayes-ens")
  attr(lfo_results, "lfo_stamp") <- .lfo_stamp
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

# The DAILY-ISSUE BACKTEST WAS REMOVED (2026-09-19).
# run_invasion_lfo_daily() backtested the FREQUENTIST INVASION_MODELS, which no longer exist,
# so the block was gated on `.have_freq` and never ran; its output lfo_daily_results.rds was
# read by nothing. It was also not the thing it claimed to be: it issued at offsets
# c(0, 3, 6) days into the week, and on the first two the origin week is still "future" to the
# truncation model, so 2 of its 3 issue-folds trained on an UNCORRECTED origin week while the
# live system always corrects it — the opposite of "mirrors the LIVE operating point exactly".
# The weekly LFO now exercises exactly the deployed operating point: its fold origin is
# cut + 6, the last day of the last training week, which is the same geometry as the live
# grid (last week ends on ANALYSIS_DATE) and therefore the same nowcast regime (2.631x).
daily_lfo_results <- NULL

# ---------------------------------------------------------------------------
# 11. Invasion evaluation (pooled; AUC-PR, ranking, log-score, calibration)
# ---------------------------------------------------------------------------
.phase("Step 10  LFO-CV")
message("\n=== Step 11: Invasion evaluation ===")

# ── Post-hoc recalibration (16b) ─────────────────────────────────────────────
# PURELY ADDITIVE under the default flags: attaches a prequentially-fitted, strictly
# rank-preserving delta per (method, horizon, fold) plus the recalibrated probability,
# so evaluate_invasion() can report recalibrated proper scores ALONGSIDE the raw ones.
# No existing metric changes. Selection moves under INVASION_SELECT_ON_RECAL and the live
# forecast under INVASION_RECALIBRATE_DEPLOY — BOTH DEFAULT TRUE since 2026-09-11 (they
# answer the same question and must not disagree; see the rationale in 00_config.R).
recal_tbl <- NULL
if (isTRUE(INVASION_RECALIBRATE)) {
  lfo_results <- tryCatch(attach_invasion_recalibration(lfo_results),
    error = function(e) { warning("[recal] attach failed: ", e$message,
                                  " — continuing with raw probabilities only."); lfo_results })
  recal_tbl <- tryCatch(invasion_delta_table(lfo_results),
    error = function(e) { warning("[recal] delta table failed: ", e$message); NULL })
  if (!is.null(recal_tbl) && nrow(recal_tbl)) {
    readr::write_csv(recal_tbl, file.path(OUT_DIAGNOSTICS, "invasion_recalibration.csv"))
    # Deployment deltas, keyed method -> horizon, for the live forecast path.
    jsonlite::write_json(
      list(generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
           estimator = paste("Jeffreys-penalised maximum likelihood (Firth), p -> 1 - (1-p)^delta",
                             "on the LFO at-risk pool; no event or origin floors"),
           note = paste("POOLED over all folds — for recalibrating the LIVE forecast only.",
                        "Scoring and selection use the prequential delta_preq attached to",
                        "lfo_results, never these."),
           band = RECAL_BAND, min_folds = RECAL_MIN_FOLDS,
           delta = recal_tbl %>% dplyr::transmute(method, horizon, delta, delta_lo, delta_hi,
                                                  n_events, boundary_hit)),
      file.path(OUT_DIAGNOSTICS, "invasion_recalibration.json"),
      auto_unbox = TRUE, pretty = TRUE, digits = 6)
    # Built from the horizons actually present, so a change to LFO_HORIZONS cannot
    # leave this line reporting NA for a horizon that was never scored.
    .dmed <- vapply(sort(unique(recal_tbl$horizon)), function(hh)
      stats::median(recal_tbl$delta[recal_tbl$horizon == hh], na.rm = TRUE), numeric(1))
    message(sprintf("[recal] pooled delta fitted for %d method x horizon slices (median by horizon: %s)",
                    nrow(recal_tbl),
                    paste(sprintf("h%d=%.3f", sort(unique(recal_tbl$horizon)), .dmed),
                          collapse = ", ")))
  }
  # DELTA-STABILITY DIAGNOSTIC. delta is an intercept shift on the cloglog scale and beta_0 is
  # the intercept, so the per-fold delta IS a time-varying beta_0 fitted post hoc. Fitting it
  # on each fold ALONE (not cumulatively, which cannot separate drift from an expanding window
  # converging) and testing for heterogeneity therefore answers two design questions at once:
  # whether the deployment delta should use all history or a recent window, and whether a
  # time-varying beta_0 is supported. Computed HERE and written to disk; the figure reads it
  # and computes nothing.
  .dstab <- tryCatch(invasion_delta_stability(lfo_results),
                     error = function(e) { warning("[recal] delta-stability diagnostic failed: ",
                                                   conditionMessage(e), call. = FALSE); NULL })
  if (!is.null(.dstab) && nrow(.dstab$summary)) {
    readr::write_csv(.dstab$per_fold,
                     file.path(OUT_DIAGNOSTICS, "invasion_delta_stability.csv"))
    readr::write_csv(.dstab$summary,
                     file.path(OUT_DIAGNOSTICS, "invasion_delta_stability_summary.csv"))
    .hs <- sort(unique(.dstab$summary$horizon))
    jsonlite::write_json(
      list(generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
           question = paste("Does the calibration offset vary between folds? If not, the",
                            "deployment delta should pool ALL history (no window) and a",
                            "time-varying beta_0 is not supported by these data."),
           method = paste("delta fitted independently per fold; DerSimonian-Laird random",
                          "effects on log delta; se from the observed information."),
           note_not_tested = paste("The anchor-week nowcast multiplier is CONSTANT across",
                                   "folds by construction (same lag from each cutoff, shared",
                                   "as-of delay fit), so it cannot be regressed on."),
           by_horizon = lapply(.hs, function(hh) {
             z <- .dstab$summary[.dstab$summary$horizon == hh, ]
             list(horizon = hh, n_methods = nrow(z),
                  median_I2 = stats::median(z$I2, na.rm = TRUE),
                  median_tau2 = stats::median(z$tau2, na.rm = TRUE),
                  n_heterogeneous_p05 = sum(z$Q_p < 0.05, na.rm = TRUE),
                  n_trend_p05 = sum(z$trend_p < 0.05, na.rm = TRUE),
                  n_extent_p05 = sum(z$extent_p < 0.05, na.rm = TRUE),
                  median_se_log_delta = stats::median(z$median_se_log, na.rm = TRUE)) })),
      file.path(OUT_DIAGNOSTICS, "invasion_delta_stability.json"),
      auto_unbox = TRUE, pretty = TRUE, digits = 6)
    for (hh in .hs) {
      z <- .dstab$summary[.dstab$summary$horizon == hh, ]
      message(sprintf(paste0("[recal] delta stability h=%s: median I2 %.1f%%, median tau2 %.4f; ",
                             "heterogeneous (Q p<0.05) in %d/%d methods, trend in %d/%d. %s"),
                      hh, 100 * stats::median(z$I2, na.rm = TRUE),
                      stats::median(z$tau2, na.rm = TRUE),
                      sum(z$Q_p < 0.05, na.rm = TRUE), nrow(z),
                      sum(z$trend_p < 0.05, na.rm = TRUE), nrow(z),
                      if (stats::median(z$I2, na.rm = TRUE) < 0.25)
                        "=> pool ALL folds; a window and a time-varying beta_0 are unsupported."
                      else "=> delta MOVES; revisit the pooled deployment factor."))
    }
  }
}

# ---- STRUCTURAL BASELINES -> lfo_results, BEFORE the evaluation table is built --------
if (!is.null(lfo_results)) {
  # NAIVE comparator for the prioritisation (figure-3, panel B) curves: rank at-risk zones
  # purely by mobility inflow FROM THE EPICENTRE, ignoring case data. Injected as an extra
  # "method" on the SAME folds/outcomes so its detection curve overlays the model's — the
  # model must beat this structural, incidence-free baseline.
  #
  # THEY ARE SCORED, NOT JUST DRAWN (2026-09-23). These rows used to be appended to a LOCAL
  # copy (lfo_fig3) built after invasion_evaluation.csv had already been written, so the
  # baselines the manuscript displays carried no AUC-PR, no AUC-ROC and no rank-of-truth --
  # the leaderboard showed a different, superseded baseline family (B4/B1/B7) from the
  # figures, and no published table could reproduce the figures' comparator. They are now
  # appended to lfo_results BEFORE evaluate_invasion(), so one baseline vocabulary reaches
  # the leaderboard, the discrimination panels and the detection curves alike.
  #
  # Two things make that safe. (1) append_naive_detection_curve_model() zero-fills
  # unobserved destinations and reports each baseline's coverage, so a rank-only row can no
  # longer shrink the shared support for every model. (2) best_invasion_model() gates
  # candidacy on prob_calibrated, which these rows set FALSE, so a baseline is ranked but can
  # never become the featured model whose score is plotted as a probability.
  # THREE STRUCTURAL BASELINES, differing ONLY in the connectivity matrix. All three rank the
  # at-risk zones by how strongly they connect to the EPICENTRE, so the comparison isolates
  # the notion of connectivity and nothing else. None reads incidence; none touches the
  # renewal machinery (no generation time, no import coefficient beta, no R(t), no Lambda).
  #
  #   (1) Gravity      — M4, the FITTED gravity kernel. Estimated on the Flowminder relocation
  #                      OD under a censored likelihood:
  #                        log(flow) = b0 + 0.518*log(pop_i) + 0.586*log(pop_j)
  #                                       - 1.438*log(dist_ij + 1)
  #                      Both masses present, all three exponents estimated.
  #   (2) Flowminder   — M_cohort, Flowminder COHORT subscriber-day presence from the
  #                      epicentre. Observed mobility, no model, measured during the outbreak.
  #   (3) Travel time  — the OSRM road travel-time matrix, 1/(1 + minutes to the nearest
  #                      epicentre zone). Pure geography.
  #
  # WHY NOT THE PREVIOUS SET. "Gravity-B4" is not a gravity model: its population mass is
  # indexed by SOURCE only (no destination mass) and its exponents are hard-coded at 0.5/1.0
  # and never fitted, so its ranking is unparameterised. "Adjacency-B7" runs the
  # method = "distance" branch (nearest-AFFECTED travel time) — there is no contiguity matrix
  # anywhere in this pipeline, so its adjacency branch is unreachable. And the inflow baseline
  # used MOBILITY_PRIMARY = M8-fill, whose epicentre rows are Flowminder ONLY over the
  # destinations Flowminder measured: measured on the shipped kernels, 28-39% of each
  # epicentre row's mass is gravity fill (Bunia 0.2795, Mongbwalu 0.3933, Rwampara 0.3807).
  # That made the "Flowminder" and "gravity" baselines partly the SAME model. M1's epicentre
  # rows are pure Flowminder; M8's are identical to M1's (max|diff| 4e-16) but M8-FILL's are not.
  #
  # B4/B1/B7 remain SCORED in invasion_evaluation.csv as comparators; they are simply no
  # longer the baselines the manuscript figures display.
  .bl_specs <- list(
    list(lbl = "Baseline-gravity",
         sc  = function() naive_epicentre_inflow_scores(mobility_matrices[["M4"]],
                            EPICENTRE_ZONES, pop_vec, zones_all)),
    list(lbl = "Baseline-flowminder-inflow",
         # THE COHORT KERNEL, not the short-trip annex. Flowminder cohort subscriber-day
         # presence is the best mobility data this study has: measured DURING the outbreak
         # (COHORT_WINDOW = "followup"), and far less censored — on the shipped kernels it
         # reaches 305 destinations from the epicentre against M1's 142, with the same zones at
         # the top (Spearman 0.79). Same signal, twice the reach.
         #
         # IT IS POOLED, NOT PER-ORIGIN. build_M_cohort() assigns one destination vector to
         # every origin in a cohort, so Bunia / Mongbwalu / Rwampara carry BIT-IDENTICAL rows
         # (verified: max pairwise difference 0), exactly as M1 does. The population weighting
         # in naive_epicentre_inflow_scores() is therefore a pure scalar here and does not
         # change the ranking. split_cohort_rows() is the per-origin object; this baseline
         # deliberately does not use it, because the split is itself a modelling step and a
         # structural null should carry as little modelling as possible.
         # Falls back to M1 (with a warning) only if the cohort build was unavailable.
         sc  = function() {
           .W <- mobility_matrices[["M_cohort"]]
           if (is.null(.W)) {
             warning("[baseline] M_cohort unavailable; the Flowminder-inflow baseline falls ",
                     "back to the pooled short-trip kernel M1, which reaches far fewer ",
                     "destinations. Check INCLUDE_COHORT_MODELS / COHORT_SOURCES.",
                     call. = FALSE)
             .W <- mobility_matrices[["M1"]]
           }
           naive_epicentre_inflow_scores(.W, EPICENTRE_ZONES, pop_vec, zones_all)
         }),
    list(lbl = "Baseline-travel-time",
         sc  = function() epicentre_travel_time_scores(osrm_mat, EPICENTRE_ZONES, zones_all))
  )
  .bl_ok <- character(0)
  for (.b in .bl_specs) {
    .sc <- tryCatch(.b$sc(), error = function(e) {
      warning(sprintf("[baseline] %s scores failed: %s", .b$lbl, conditionMessage(e)),
              call. = FALSE); NULL })
    if (is.null(.sc) || !any(is.finite(.sc) & .sc > 0)) {
      warning(sprintf(paste0("[baseline] %s is unavailable or all-zero and will NOT appear on ",
                             "the figures. Manuscript Figure 2 expects all three structural ",
                             "baselines."), .b$lbl), call. = FALSE)
      next
    }
    lfo_results <- append_naive_detection_curve_model(lfo_results, .sc, .b$lbl)
    if (.b$lbl %in% lfo_results$method) .bl_ok <- c(.bl_ok, .b$lbl)
  }
  message(sprintf("[baseline] structural baselines scored and on the figures: %s",
                  if (length(.bl_ok)) paste(.bl_ok, collapse = ", ") else "NONE"))
}

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

# ── TRUNCATION SENSITIVITY: the same evaluation on SETTLED folds only ─────────────────
# The cross-validation now runs to the last round whose outcome window closes by the analysis
# date (LFO_MIN_EVAL_AGE_DAYS, 00_config.R). The newest rounds are therefore right-truncated:
# a zone invaded late in an outcome week that has had only days to report can be scored
# is_new_invasion = 0. This re-scores every model on the rounds whose outcome window is at
# least LFO_RELIABLE_EVAL_AGE_DAYS old — the rounds the previous, more conservative fold
# window admitted — so the size of the effect on every published metric is visible rather
# than argued about. Nothing downstream reads this table: it is a sensitivity, not the result.
local({
  if (is.null(lfo_results) || !"eval_reliable" %in% names(lfo_results)) return(invisible(NULL))
  rel <- lfo_results[lfo_results$eval_reliable %in% TRUE, , drop = FALSE]
  n_drop <- dplyr::n_distinct(lfo_results$fold_id) - dplyr::n_distinct(rel$fold_id)
  if (!nrow(rel) || n_drop == 0L) {
    message("[eval] every fold's outcome window is settled; no truncation sensitivity needed.")
    return(invisible(NULL))
  }
  et <- tryCatch(evaluate_invasion(rel),
                 error = function(e) { warning("reliable-fold eval failed: ", e$message); NULL })
  if (is.null(et)) return(invisible(NULL))
  readr::write_csv(et, file.path(OUT_DIAGNOSTICS, "invasion_evaluation_reliable_folds.csv"))
  message(sprintf(paste0("[eval] truncation sensitivity written: %d of %d rounds have a settled ",
                         "(>= %d d) outcome window -> invasion_evaluation_reliable_folds.csv"),
                  dplyr::n_distinct(rel$fold_id), dplyr::n_distinct(lfo_results$fold_id),
                  as.integer(get0("LFO_RELIABLE_EVAL_AGE_DAYS", ifnotfound = 25L))))
})

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
    # The base rates and the common-denominator gap are exported alongside the headline gap.
    # inner_cv_skill and heldout_skill are each scored against their OWN base rate, so their
    # difference mixes selection optimism with a shift in how hard the window is; on the
    # shipped folds that shift is 37% and it roughly halves the reported gap.
    tibble::tibble(selected = heldout_eval$selected, horizon = heldout_eval$horizon,
                   inner_cv_skill = heldout_eval$inner_skill,
                   heldout_skill  = heldout_eval$outer_skill,
                   optimism_gap   = heldout_eval$optimism_gap,
                   optimism_gap_common_base = heldout_eval$optimism_gap_common_base,
                   inner_base_rate = heldout_eval$inner_base_rate,
                   heldout_base_rate = heldout_eval$outer_base_rate,
                   base_rate_ratio = heldout_eval$base_rate_ratio,
                   inner_auc_pr = heldout_eval$inner_auc_pr,
                   heldout_auc_pr = heldout_eval$outer_auc_pr,
                   n_holdout_origins = heldout_eval$n_holdout),
    file.path(OUT_DIAGNOSTICS, "invasion_heldout_optimism.csv"))
  message(sprintf(paste0("[held-out §2.2] %s: inner-CV AUC-PR skill %.1fx vs held-out (last 2 origins) ",
                         "%.1fx (optimism gap %.1fx; %.1fx on a common base rate - the held-out base ",
                         "rate is %.2fx the inner one)"),
                  heldout_eval$selected, heldout_eval$inner_skill %||% NA_real_,
                  heldout_eval$outer_skill %||% NA_real_, heldout_eval$optimism_gap %||% NA_real_,
                  heldout_eval$optimism_gap_common_base %||% NA_real_,
                  heldout_eval$base_rate_ratio %||% NA_real_))
}

# NOTE: the forecast-vs-observed + prospective block (§3.4/§3.5/§3.6/§3.8) is placed
# AFTER best_bayes_method is defined (below the headline-model selection), since it keys
# on the featured Bayesian method.

# Headline model (best overall) + best renewal model, both by a
# calibration-aware criterion (discrimination + targeting + log-score) POOLED across
# both forecast horizons, so we do not feature an over-confident variant or one that
# is strong at 1-week but weak at 2-week (or vice versa).
# There is no frequentist model to fall back to, so a missing/empty evaluation leaves
# best_method NA and the Bayesian fallback immediately below supplies it.
best_method <- if (!is.null(eval_tbl)) best_invasion_model(eval_tbl) else NA_character_
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
# Degenerate fallback: if the evaluation table is missing or empty, best_invasion_model()
# returns NA. best_method drives the skill/lead/detection/predobs figures, so it must be a real
# cross-validated method — fall back to the featured Bayesian pick.
if ((length(best_method) != 1 || is.na(best_method)) && !is.null(best_bayes_method))
  best_method <- best_bayes_method
message(sprintf("[eval] Headline model: %s", best_method))

# ── Forecast-vs-observed suite + prospective check (review §3.4/§3.5/§3.6/§3.8) ──
# Uses the featured Bayesian model's LFO forecasts (predicted vs realised). Guarded and
# additive; produces CSVs that feed the reliability / count-calibration / prospective figures.
# Placed here (not earlier) because it keys on best_bayes_method, defined just above.
if (!is.null(lfo_results) && !is.null(best_bayes_method) &&
    best_bayes_method %in% lfo_results$method) tryCatch({
  .rel_dir <- function(raw) if (raw) { d <- file.path(OUT_DIAGNOSTICS, "raw")
                                       if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
                                       d } else OUT_DIAGNOSTICS
  # BOTH horizons: FigR1_reliability is a retained deliverable at h=1 AND h=2, and the
  # horizon was previously hard-coded to 1 here and again in 24_review_figures.R, so the
  # h=2 panel could never be produced. The basenames now carry the horizon.
  for (.hh in LFO_HORIZONS) {
  fvo <- lfo_results %>% dplyr::filter(method == best_bayes_method, horizon == .hh,
                                       !was_active_before, is.finite(p_invasion))
  if (nrow(fvo)) {
    # (§3.5 view 1) reliability curve; (view 2) aggregate-count calibration over folds.
    # Written on BOTH probability scales, same basenames, the raw one under diagnostics/raw/.
    # The PRIMARY files carry the PREQUENTIAL recalibration — the leakage-free retrospective
    # scale, fitted for each fold only on folds whose outcome had already closed — because
    # that is the forecast the deployed system issues. The raw twin is what motivates the
    # correction, and 24_review_figures.R overlays the two so the calibration panels show
    # the defect and the fix together rather than one without the other.
    .pcol <- if ("p_recal" %in% names(fvo) && all(is.finite(fvo$p_recal))) "p_recal" else {
      warning("[eval §3.5] no usable recalibrated probability for the featured model; the ",
              "calibration files are written on the RAW scale.", call. = FALSE); "p_invasion" }
    for (.v in list(list(col = .pcol, raw = FALSE), list(col = "p_invasion", raw = TRUE))) {
      .p  <- fvo[[.v$col]]
      rc  <- reliability_curve(.p, as.integer(fvo$is_new_invasion), n_bins = 10L)
      if (nrow(rc)) readr::write_csv(rc, file.path(.rel_dir(.v$raw),
                       sprintf("forecast_reliability_h%d.csv", .hh)))
      acc <- aggregate_count_calibration(data.frame(time = as.character(fvo$fold_id),
               p_invasion = .p, invaded = as.integer(fvo$is_new_invasion)))
      if (nrow(acc)) readr::write_csv(acc, file.path(.rel_dir(.v$raw),
                        sprintf("forecast_count_calibration_h%d.csv", .hh)))
    }
  }
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
  # Calibration-in-the-large of the PRIMARY (recalibrated) curve, so the annotation and the
  # points it annotates are on one scale; the raw value travels alongside for the overlay.
  # PER HORIZON, named by horizon. This was a scalar computed at h == 1 and handed to both
  # panels, so FigR1_reliability_h2 would have printed h=1's calibration-in-the-large against
  # h=2 data (1.32/2.03 where the truth is 1.55/2.36). make_review_figures()'s .cal_for()
  # already reads a horizon-named vector.
  .cil <- function(col) {
    if (is.null(eval_tbl) || is.null(best_bayes_method) || !col %in% names(eval_tbl))
      return(NULL)
    stats::setNames(vapply(LFO_HORIZONS, function(h) {
      v <- eval_tbl[[col]][eval_tbl$method == best_bayes_method & eval_tbl$horizon == h]
      if (length(v)) as.numeric(v[1]) else NA_real_
    }, numeric(1)), as.character(LFO_HORIZONS))
  }
  .cal_raw <- .cil("calibration_in_large")
  # Fall back to the raw value PER HORIZON where the recalibrated one is unavailable.
  .cal     <- local({
    v <- .cil("calibration_in_large_recal")
    if (is.null(v)) .cal_raw
    else { bad <- !is.finite(v); if (any(bad) && !is.null(.cal_raw)) v[bad] <- .cal_raw[bad]; v }
  })
  make_review_figures(OUT_DIAGNOSTICS,
                      file.path(OUT_DIR, "key_outputs", "manuscript_figures", "panels"),
                      cal_in_large = .cal, cal_in_large_raw = .cal_raw,
                      top_k = 15L, horizons = LFO_HORIZONS)
}, error = function(e) warning("review figures: ", e$message))

# Curated display set for the (otherwise cluttered) multi-model figures: the featured
# model, the best few Bayesian variants, both Bayesian ensembles and the three structural
# comparators. The FULL grid remains in invasion_evaluation.csv.
eval_display <- eval_tbl
if (!is.null(eval_tbl)) {
  top_bayes <- {
    .tb <- eval_tbl %>% dplyr::filter(horizon == 1, grepl("^Bayes", method)) %>%
      dplyr::arrange(dplyr::desc(auc_pr_skill)) %>% dplyr::pull(method)
    unique(c(head(.tb, 6L), "Bayes-ens-mean", "Bayes-ens-median"))
  }
  DISPLAY_MODELS <- unique(c(best_method, top_bayes,
                             "Gravity-B4", "Distance-B1", "Adjacency-B7"))
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
# STALE FREQUENTIST PRODUCTS MUST NOT SURVIVE. The frequentist arm is gone, so nothing
# refreshes these files; left alone they sit in outputs/reports/ looking current — six weeks
# stale on the audited tree, with no date stamp to warn a reader. Delete them: this file's own
# stated principle is to fail loudly rather than silently emit stale products. The Bayesian
# suite supplies the featured maps and tables (bayes_risk_scores_all_zones.csv et al.).
local({
  .stale_freq <- c(file.path(OUT_REPORTS, c("risk_scores_all_zones.csv", "risk_table_national.csv",
                                            "risk_table_ituri.csv", "priority_table.csv")),
                   file.path(OUT_FORECASTS, "covariate_associations.rds"))
  .stale_freq <- .stale_freq[file.exists(.stale_freq)]
  if (length(.stale_freq)) {
    message(sprintf("[cleanup] removing %d stale frequentist product(s) (that arm is retired): %s",
                    length(.stale_freq), paste(basename(.stale_freq), collapse = ", ")))
    suppressWarnings(file.remove(.stale_freq))
  }
  # Same reasoning for risk_scores_current.rds: it was written by the frequentist branch and
  # was read unguarded downstream, so a leftover copy fails with an opaque pipe error instead
  # of saying the arm did not run.
  suppressWarnings(file.remove(file.path(OUT_FORECASTS, "risk_scores_current.rds")))
})
# Q2 + Q3: 0-1 relative-risk indices (1 = highest) and the vulnerability-and-capacity-adjusted
# preparedness-priority score (invasion risk x vulnerability). The BAYESIAN suite consumes this.
vuln_index <- tryCatch(compute_vulnerability_index(covariates, zones_all, osrm_mat = osrm_mat),
                       error = function(e) { warning("vulnerability index: ", e$message); NULL })

# Full per-zone invasion probabilities + risk scores for ALL health zones (every
# zone x both horizons; already-affected zones carry NA by construction), exported
# to CSV alongside the RDS so the complete table is usable outside R.
# `prob_scale` stamps which probability scale the file is on. Two are in circulation and
# they are NOT interchangeable: "recalibrated-deployed" carries the pooled post-hoc factor
# (16b) that the live forecast is published on, "raw" is the uncorrected model output kept
# as a twin. Ranks, mu-based relative risks and the priority index are (near-)invariant;
# the probability columns are not. A file without the stamp predates this convention.
# Sibling location for the uncorrected twin of any primary artifact: <dir>/raw/<basename>.
# Basenames are preserved deliberately, so a raw file is the same file on the other scale
# and nothing downstream has to learn a second naming convention. Directory is created on
# demand; NULL in, NULL out.
.raw_path <- function(path) {
  if (is.null(path) || !length(path) || is.na(path)) return(NULL)
  d <- file.path(dirname(path), "raw")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, basename(path))
}

.write_risk_csv <- function(rs, path, prob_scale = "recalibrated-deployed") {
  if (is.null(rs) || !nrow(rs)) return(invisible(NULL))
  rs$prob_scale <- prob_scale
  # p_lo/p_hi are the 5% and 95% bounds on p_case_invasion (the 90% CrI for the Bayesian scores;
  # the ensemble min-max spread for the frequentist/ensemble scores) — exported as p_case_lo /
  # p_case_hi. p_infection_invasion is intentionally NOT exported.
  keep <- c("health_zone", "province", "horizon", "was_active_before", "prob_scale",
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
  # SIGNIFICANT digits, not 4 decimal places. round(x, 4) destroyed every probability below
  # 5e-5: 34 of the 458 at-risk zones were published with p_case_invasion and mu_forecast
  # EXACTLY 0, and 4 with rr01_nat exactly 0. A renewal hazard gives p = 1 - exp(-mu) > 0
  # strictly, so every one of those zeros was an artefact — and in this file NA already means
  # "already affected", so 0 read as "not at risk" while the same row carried a positive
  # rr_nat. The figure suites and the conditional-arm ratio analysis read this CSV, so the
  # artefact propagated into Figures 1-3, the top-15 panel and a 0/0 ratio. signif() keeps
  # small probabilities resolved while still trimming float noise from the large columns.
  out <- rs %>% dplyr::select(dplyr::any_of(keep)) %>%
    dplyr::rename(dplyr::any_of(c(p_case_lo = "p_lo", p_case_hi = "p_hi"))) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 6)))
  if (all(c("horizon", "p_case_invasion") %in% names(out)))
    out <- out %>% dplyr::arrange(horizon, dplyr::desc(dplyr::coalesce(p_case_invasion, -1)))
  readr::write_csv(out, path)
  message(sprintf("[export] %d zone-rows (%d zones x horizons) -> %s",
                  nrow(out), dplyr::n_distinct(out$health_zone), basename(path)))
  invisible(out)
}

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

# The per-province and preparedness-priority risk tables were built from the frequentist
# featured model's risk_scores (and still referenced the retired p_infection_invasion column).
# The Bayesian suite writes the equivalents: bayes_risk_scores_all_zones.csv and the priority
# products in section 12.
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
  # Deployment recalibration (16b), ON by default since 2026-09-11. Each suite model's
  # LIVE forecast is rescaled by its own pooled delta inside the posterior draws, so
  # rankings, top-K watch-lists and rank credible intervals are unchanged while absolute
  # probabilities move onto the observed frequency scale.
  # NOTE the cascade is deliberately NOT affected: 33_cascade_eval.R calls
  # predict_bayes_invasion() without a delta and cascade_fit_delta() derives its own
  # factor from the RAW calibration_in_large column of invasion_evaluation.csv (which is
  # computed from the LFO frame and is itself unaffected by this flag), so the correction
  # is applied once in each layer and never compounded across them.
  .deploy_delta <- NULL
  if (isTRUE(INVASION_RECALIBRATE_DEPLOY)) {
    if (is.null(recal_tbl) || !nrow(recal_tbl)) {
      warning("[recal] INVASION_RECALIBRATE_DEPLOY is TRUE but no delta table is available ",
              "(is INVASION_RECALIBRATE on, and did the LFO succeed?); the live Bayesian ",
              "forecast is NOT recalibrated.", call. = FALSE)
    } else {
      .deploy_delta <- recal_tbl[, c("method", "horizon", "delta")]
      message(sprintf("[recal] DEPLOYING pooled delta to the live Bayesian suite (%d method x horizon rows).",
                      nrow(.deploy_delta)))
    }
  }
  bayes_suite <- tryCatch(
    fit_bayes_suite(zone_week_nc, mobility_matrices, gt_pmfs, covariates, osrm_mat,
                    zones_all, affected_now, LFO_HORIZONS, iter = 2000L,
                    delta_tbl = .deploy_delta, linelist = dat$ll,
                    analysis_date = ANALYSIS_DATE,
                    # Same rolling predictors the LFO uses, and the same as-of truncation:
                    # the rolling reconstruction is built by reaggregate_asof(), which censors
                    # on the observation date, so the as-of fit is the one that describes it.
                    # See the note in fit_bayes_suite() for why this must not diverge.
                    trunc_delay = .trunc_asof),
    error = function(e) { warning("bayes suite: ", e$message); NULL })
  if (!is.null(bayes_suite)) {
    bayes_stack <- tryCatch(bayes_stacked_prediction(bayes_suite), error = function(e) NULL)
    # Raw twin of the stacked ensemble: the same weights and the same combination rule
    # applied to the uncorrected member summaries. Built from `preds_raw`, not by
    # un-transforming the stacked result — the stack averages posterior MEANS, and the
    # recalibration is nonlinear in p, so un-transforming a stacked mean is not the stack
    # of the un-transformed means.
    bayes_stack_raw <- if (!is.null(bayes_suite$preds_raw) && nrow(bayes_suite$preds_raw))
      tryCatch(bayes_stacked_prediction(utils::modifyList(bayes_suite,
                 list(preds = bayes_suite$preds_raw))), error = function(e) NULL) else NULL
    saveRDS(bayes_suite$preds,  file.path(OUT_FORECASTS, "bayes_current_predictions.rds"))
    saveRDS(bayes_suite$params, file.path(OUT_FORECASTS, "bayes_parameters.rds"))
    if (!is.null(bayes_suite$weights))
      saveRDS(bayes_suite$weights, file.path(OUT_FORECASTS, "bayes_stacking_weights.rds"))
    if (!is.null(bayes_stack)) saveRDS(bayes_stack, file.path(OUT_FORECASTS, "bayes_stacked_current.rds"))
    readr::write_csv(bayes_suite$params, file.path(OUT_REPORTS, "bayes_parameters.csv"))
    # POSTERIOR TRAJECTORY OF beta_t FOR EVERY SUITE MODEL (fixed-beta models included, as the
    # flat reference the time-varying variants are read against). This is the numeric content
    # behind the time-varying-beta figure, and the only artefact that records what the fitted
    # import->invasion conversion rate did over the epidemic — the fits themselves are not
    # persisted.
    if (!is.null(bayes_suite$beta_traj) && nrow(bayes_suite$beta_traj)) {
      readr::write_csv(bayes_suite$beta_traj %>%
          dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 12))),
        file.path(OUT_DIR, "key_outputs", "bayes_beta_trajectory.csv"))
      message(sprintf("[bayes] beta_t trajectory written for %d model(s) over %d week(s)",
                      dplyr::n_distinct(bayes_suite$beta_traj$model),
                      dplyr::n_distinct(bayes_suite$beta_traj$week)))
    }
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
    # Every primary product is on the DEPLOYED (recalibrated) scale; `*_raw` carries the
    # uncorrected twin so the parallel raw artifact set can be written from the same fits.
    bayes_best_preds <- NULL; bayes_best_preds_raw <- NULL
    if (!is.null(best_bayes_method) && best_bayes_method %in% bayes_suite$preds$method) {
      # DEPLOYED-FORECAST GT MARGINALISATION (review §2.1). The mobility KERNEL is selected
      # by cross-validation at the central literature GT (so GT is NOT selected — the
      # reviewer's point); here we propagate generation-time uncertainty into the DEPLOYED
      # featured forecast by marginalising over GT_PRIOR (refit at each GT grid point, draws
      # pooled). This widens the per-zone credible intervals and gives the ranking a genuine
      # (non-degenerate) rank credible interval, without re-opening the CV skill. Gated by
      # GT_MARGINALISE_FEATURED, which is **FALSE by default since 2026-09-22**: the generation
      # time is treated as a KNOWN distribution and this branch does not run. It is kept as a
      # sensitivity arm — see the GT_PRIOR block in 00_config.R for why an assumed prior's
      # width is not uncertainty we are willing to report as if it were measured. When the
      # flag is on, it falls back to the fixed-GT prediction on any error.
      .fm_i <- which(vapply(bayes_suite$grid, function(g) identical(g$label, best_bayes_method), logical(1)))
      .fm   <- if (length(.fm_i)) bayes_suite$grid[[.fm_i[1]]] else NULL
      if (!is.null(.fm) && isTRUE(get0("GT_MARGINALISE_FEATURED", ifnotfound = FALSE))) {
        .mp <- tryCatch({
          # The featured model's OWN deployment delta must travel with it. This branch
          # replaces the suite's prediction for the featured model, and the suite applied
          # the delta inside fit_bayes_suite(); without passing it here the featured
          # forecast would be the only deployed product left on the raw scale. scales=TRUE
          # returns both from ONE set of GT-grid fits (the expensive part), so the raw twin
          # costs a second summarisation rather than a second marginalisation.
          mp <- predict_bayes_gt_marginal(zone_week_nc, mobility_matrices, covariates, osrm_mat,
                  zones_all, mob = .fm$mob, cov_spec = .fm$cov, horizons = LFO_HORIZONS,
                  affected_zones = affected_now, gt_prior = GT_PRIOR, link = .fm$link %||% "cloglog",
                  tv = .fm$tv %||% "none", scales = TRUE,
                  linelist = dat$ll, analysis_date = ANALYSIS_DATE,
                  delta = .suite_delta_for(.deploy_delta, best_bayes_method, LFO_HORIZONS))
          if (!is.null(mp)) {
            mp$calibrated$method <- best_bayes_method
            mp$raw$method        <- best_bayes_method
            message(sprintf("[gt-marginal] deployed featured forecast marginalised over GT_PRIOR (%d grid points) for %s",
                            length(make_gt_prior_pmfs(GT_PRIOR)$gt_pmfs), best_bayes_method))
          }
          mp
        }, error = function(e) { warning("[gt-marginal] featured marginalisation failed; using fixed-GT: ", e$message); NULL })
        if (!is.null(.mp)) {
          bayes_best_preds     <- .mp$calibrated
          bayes_best_preds_raw <- .mp$raw
        }
      }
      if (is.null(bayes_best_preds)) {
        bayes_best_preds     <- bayes_suite$preds %>% dplyr::filter(method == best_bayes_method)
        bayes_best_preds_raw <- if (!is.null(bayes_suite$preds_raw) && nrow(bayes_suite$preds_raw))
          bayes_suite$preds_raw %>% dplyr::filter(method == best_bayes_method) else NULL
      }
    } else if (!is.null(bayes_stack)) {
      bayes_best_preds     <- bayes_stack                      # fallback
      bayes_best_preds_raw <- bayes_stack_raw
    }
    if (!is.null(bayes_best_preds)) {
      bayes_risk_scores <- tryCatch({
        rs <- compute_risk_scores(bayes_best_preds, province_map)
        if (!is.null(vuln_index)) rs <- add_risk_indices(rs, vuln_index)
        rs
      }, error = function(e) { warning("bayes risk scores: ", e$message); NULL })
      if (!is.null(bayes_risk_scores)) {
        saveRDS(bayes_risk_scores, file.path(OUT_FORECASTS, "bayes_risk_scores_current.rds"))
        # Full per-zone Bayesian invasion probabilities + risk scores (all zones) to CSV,
        # on the DEPLOYED (recalibrated) scale.
        .brs_path <- file.path(OUT_REPORTS, "bayes_risk_scores_all_zones.csv")
        tryCatch(.write_risk_csv(bayes_risk_scores, .brs_path),
                 error = function(e) warning("bayes risk CSV: ", e$message))
        # Uncorrected twin, same basename, under reports/raw/. Built through the SAME
        # risk-score machinery so the two files differ only in the probability scale and in
        # nothing else (province join, vulnerability index, relative risks, priority).
        if (!is.null(bayes_best_preds_raw)) tryCatch({
          .rs_raw <- compute_risk_scores(bayes_best_preds_raw, province_map)
          if (!is.null(vuln_index)) .rs_raw <- add_risk_indices(.rs_raw, vuln_index)
          saveRDS(.rs_raw, file.path(OUT_FORECASTS, "bayes_risk_scores_current_raw.rds"))
          .write_risk_csv(.rs_raw, .raw_path(.brs_path), prob_scale = "raw")
        }, error = function(e) warning("bayes risk CSV (raw twin): ", e$message))
        # Pairwise DIRECTED importation pressure / force of infection between zones,
        # decomposed from the featured model's renewal equation at its posterior-median
        # import coefficient: sum over origins of foi == that model's per-week destination
        # hazard. NOTE the scale: this uses the RAW posterior beta, whereas the deployed
        # forecast additionally carries the 16b recalibration factor, so `foi` is the
        # model-scale network behind mu_forecast rather than its deployed-scale twin. The
        # share_of_dest column, which is what this table is read for, is invariant to that.
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
              rt_draws = bayes_rt_week_draws(dat$ll, gt_pmfs[[.fm$gt]],
                                             week_start = max(zone_week_nc$week_start),
                                             issue_date = ANALYSIS_DATE),
              affected_zones = affected_now, province_map = province_map)
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
      if (!is.null(bayes_stack_risk_scores)) {
        .bes_path <- file.path(OUT_REPORTS, "bayes_ensemble_risk_scores_all_zones.csv")
        tryCatch(.write_risk_csv(bayes_stack_risk_scores, .bes_path),
                 error = function(e) warning("bayes ensemble risk CSV: ", e$message))
        if (!is.null(bayes_stack_raw)) tryCatch({
          .es_raw <- compute_risk_scores(bayes_stack_raw, province_map)
          if (!is.null(vuln_index)) .es_raw <- add_risk_indices(.es_raw, vuln_index)
          .write_risk_csv(.es_raw, .raw_path(.bes_path), prob_scale = "raw")
        }, error = function(e) warning("bayes ensemble risk CSV (raw twin): ", e$message))
      }
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

# Persist the in-memory LFO object so the detail visualisations (and any ad-hoc
# re-plotting) can run without repeating the expensive LFO step. fc_all_current.rds held
# the frequentist current-week forecast; remove any copy an older run left behind.
suppressWarnings(file.remove(file.path(OUT_FORECASTS, "fc_all_current.rds")))
if (!is.null(lfo_results)) saveRDS(lfo_results, file.path(OUT_FORECASTS, "lfo_results.rds"))

# Remove figures/tables from the superseded pre-rebuild pipeline so the outputs
# folder shows ONLY the current invasion products (no maps of removed models).
.stale <- c(
  file.path(OUT_DIAGNOSTICS, c("brier_comparison.pdf", "calibration_curves.pdf",
    "case_count_tile.pdf", "evaluation_heatmap_h1.pdf", "evaluation_heatmap_h2.pdf",
    "pr_curves.pdf", "roc_curves.pdf", "wis_by_model_1w.pdf",
    "metrics_summary.csv", "evaluation_heatmap.pdf")),
  list.files(OUT_MAPS, "^risk_map_", full.names = TRUE),
  # Written only under RUN_EPINOWCAST_DIAGNOSTIC since the deployed nowcast changed. A stale copy
  # sitting undated beside this run's products asserts that these counts were epinowcast-corrected.
  if (!.env_flag("RUN_EPINOWCAST_DIAGNOSTIC", FALSE))
    file.path(OUT_DIAGNOSTICS, "epinowcast_weekly_factors.csv") else character(0),
  file.path(OUT_REPORTS, c("methodology_report.md", "summary_dashboard.pdf",
    "table_best_per_metric.csv", "table_model_comparison.csv",
    "table_outbreak_summary.csv")))
suppressWarnings(file.remove(.stale[file.exists(.stale)]))

shapefile <- if (file.exists(SHAPEFILE_PATH))
  tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE), error = function(e) NULL) else NULL
vwrap <- function(expr) tryCatch(expr, error = function(e) warning(conditionMessage(e)))

# F1 — legible epidemic curve (windowed; top zones + Other)
vwrap(plot_epidemic_curve_legible(zone_week_nc, province_map, save = TRUE))
# F9 (REMOVED) — the epinowcast nowcast fan. The deployed nowcast is apply_nowcast_correction()
# (step 4), so epinowcast_weekly_factors.csv is written only under RUN_EPINOWCAST_DIAGNOSTIC and
# nowcast_fan.pdf is not on FIGURE_KEEP; the call rendered a figure titled "epinowcast
# right-truncation correction" from whatever copy an earlier run left behind, and then the gate
# discarded it. Draw it only when this run actually produced the factors.
if (.env_flag("RUN_EPINOWCAST_DIAGNOSTIC", FALSE))
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
    # show_title = FALSE: bayes_invasion_rank_map_national is a published panel whose
    # caption lives in the manuscript text (2026-09-17 streamlining brief).
    vwrap(plot_invasion_rank_map(rs, horizon = h, method_label = model_label,
            shapefile = shapefile, window_txt = window_txt, save = TRUE,
            file = pf("invasion_rank_map"), extent = "national", show_title = FALSE))
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
      # show_title = FALSE: published panel, caption in the manuscript text.
      vwrap(plot_prob_vuln_choropleth(rs, shapefile = shapefile, horizon = h,
              province_zoom = NULL, window_txt = window_txt, save = TRUE,
              file = pf(sprintf("prob_vuln_choropleth_national_h%d", h)),
              model_label = model_label, show_title = FALSE))
      for (prov in PROVINCES_OF_INTEREST)
        vwrap(plot_prob_vuln_choropleth(rs, shapefile = shapefile, horizon = h,
                province_zoom = prov, window_txt = window_txt, save = TRUE,
                file = pf(sprintf("prob_vuln_choropleth_%s_h%d", .prov_suffix(prov), h)),
                model_label = model_label))
    }
    if ("priority" %in% names(rs)) {
      # show_title = FALSE: published panel, caption in the manuscript text.
      vwrap(plot_priority_scatter(rs, horizon = h, window_txt = window_txt, save = TRUE,
              file = pf("priority_scatter"), model_label = model_label, show_title = FALSE))
      vwrap(plot_priority_bars(rs, horizon = h, save = TRUE,
              file = pf("priority_bars"), model_label = model_label))
      vwrap(plot_priority_map(rs, shapefile = shapefile, horizon = h, window_txt = window_txt,
              save = TRUE, file = pf("priority_map"), model_label = model_label))
    }
  }
  invisible(TRUE)
}
# render_map_suite() above is used by the BAYESIAN suite (section 12), which renders the
# invasion maps, choropleths and the preparedness-priority set for the featured model.
# The frequentist featured-model map suite and its covariate-association parameter screen
# (build_invasion_design / plot_model_parameters, and covariate_associations.rds) went with
# the renewal arm they described.

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
# Task 4 — the per-model spatial / space-time diagnostics for the featured BAYESIAN model,
# file-tagged with the bayes_ prefix, recomputing the per-method spatial error.
.diag_models <- list()
if (!is.null(best_bayes_method) && best_bayes_method %in% lfo_results$method)
  .diag_models[[length(.diag_models) + 1L]] <- list(m = best_bayes_method,
    lab = sprintf("Bayesian: %s", best_bayes_method), pfx = "bayes_")
if (!is.null(lfo_results)) {
  # The three structural baselines were appended to lfo_results before evaluate_invasion()
  # (see the STRUCTURAL BASELINES block above the eval_tbl call). lfo_fig3 is kept as the
  # name the figure code below expects; it is now simply lfo_results, baselines included.
  lfo_fig3 <- lfo_results

  # PUBLISH THE DETECTION CURVES. The manuscript and publication figure suites each drew this
  # curve by hand — re-deriving the naive epicentre-inflow score from the raw mobility kernel
  # (without the alias harmonisation naive_epicentre_inflow_scores() applies), and pooling
  # recall over folds where evaluate_invasion() averages it per fold. Three hand copies, three
  # chances to drift, and a "share of true invasions caught" on a main-text panel that no
  # published table could reproduce. compute_detection_curve() is now the single estimator and
  # its output is a published table; the figures read it. The baseline rows now DO reach
  # lfo_results (and therefore lfo_results.rds and invasion_evaluation.csv) -- see the
  # STRUCTURAL BASELINES block above the eval_tbl call -- so the curve is no longer the only
  # published trace of them. It remains the single estimator: recall is averaged per fold here
  # and in evaluate_invasion(), so a panel that annotates a recall from this curve beside a
  # recall_at_K from that table is quoting one estimand, not two hand copies of it.
  # NOTE lfo_cv_results.rds (the LFO cache, written before the append) stays model-only, which
  # is why REUSE_CACHED_LFO re-appends the baselines rather than serving them from cache.
  tryCatch({
    .dc_methods <- unique(c(best_bayes_method, .bl_ok,
                            intersect(c("Gravity-B4", "Distance-B1", "Adjacency-B7",
                                        "Bayes-ens-mean", "Bayes-ens-median"),
                                      unique(lfo_fig3$method))))
    .dc_methods <- .dc_methods[!is.na(.dc_methods)]
    # SUPPORT FROM lfo_results, which is now the SAME object lfo_fig3 points at (the baselines
    # are appended upstream of the evaluation table). The distinction this comment used to draw
    # has collapsed, and deliberately: append_naive_detection_curve_model() zero-fills a
    # baseline's unobserved destinations, so a rank-only row no longer carries NA and can no
    # longer drop cells for every method. Read each baseline's coverage line in the run log
    # beside its rank metrics -- a large zero block is a property of the data source. Resolving
    # the support here keeps the curve on exactly the rows evaluate_invasion() scored, so the
    # panel's recall and the table's recall_at_K remain one estimand. Once per horizon.
    .dc <- dplyr::bind_rows(lapply(LFO_HORIZONS, function(hh) {
      .cells <- tryCatch(invasion_common_cells(lfo_results, hh), error = function(e) NULL)
      if (!length(.cells))
        warning(sprintf(paste0("[detection] no shared support resolved from the scored table ",
                               "at h=%s: the curves there are UNRESTRICTED and their k / ",
                               "n_atrisk reference will not match invasion_evaluation.csv."), hh),
                call. = FALSE)
      # common_support = FALSE deliberately: it is the only other way to resolve a support,
      # and it would resolve it from lfo_fig3 -- the wrong table -- whenever .cells is NULL.
      dplyr::bind_rows(lapply(.dc_methods, function(m)
        tryCatch(compute_detection_curve(lfo_fig3, m, hh, ks = 1:25,
                                         common_support = FALSE, support_cells = .cells),
                 error = function(e) NULL)))
    }))
    if (!is.null(.dc) && nrow(.dc)) {
      .dc <- .dc %>% dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 12)))
      readr::write_csv(.dc, file.path(OUT_DIR, "key_outputs", "detection_curves.csv"))
      message(sprintf("[detection] published detection_curves.csv: %d method(s) x %d horizon(s) x %d budget(s)",
                      dplyr::n_distinct(.dc$method), dplyr::n_distinct(.dc$horizon),
                      dplyr::n_distinct(.dc$k)))
    }
  }, error = function(e) warning("[detection] detection_curves.csv not written: ",
                                 conditionMessage(e), call. = FALSE))
  # All per-model evaluation figures produced at BOTH horizons (1- and 2-week ahead).
  for (hh in LFO_HORIZONS) {
    vwrap(plot_reliability(lfo_results, horizon = hh, save = TRUE))
    # Delta stability: is the calibration offset constant over time? Reads the tables written
    # in the recalibration block above and computes nothing, so the panel and the published
    # numbers cannot disagree.
    vwrap(plot_delta_stability(
      file.path(OUT_DIAGNOSTICS, "invasion_delta_stability.csv"),
      file.path(OUT_DIAGNOSTICS, "invasion_delta_stability_summary.csv"),
      horizon = hh, save = TRUE))
    # Bayesian calibration/reliability + combined figure-3, featured Bayesian model. Panel A
    # (reliability) uses methods[1] = the featured Bayesian model; panel B (prioritisation)
    # overlays the naive epicentre-inflow baseline for a like-for-like comparison.
    if (!is.null(best_bayes_method) && best_bayes_method %in% lfo_results$method) {
      vwrap(plot_reliability(lfo_results, horizon = hh, methods = best_bayes_method,
              file = "bayes_reliability", save = TRUE))
      # Support from lfo_results -- the same object as lfo_fig3 now; passed explicitly so this
      # panel is drawn on exactly the rows evaluate_invasion() scored. See the detection-curve
      # block above for why the baselines no longer shrink it.
      vwrap(plot_paper_figure3(lfo_fig3, c(best_bayes_method, .bl_ok),
              horizon = hh, save = TRUE, file = "bayes_model_performance_figure3",
              support_cells = tryCatch(invasion_common_cells(lfo_results, hh),
                                       error = function(e) NULL)))
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

# The frequentist member-spread products (best-model-vs-ensemble panel, top-zone min-max
# uncertainty bands, space-time member trajectories) went with the renewal ensemble they
# visualised. The Bayesian posterior 90% CrI maps are the equivalent and are rendered by the
# Bayesian suite above.

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
    # THE ARMS MUST BE DIFFERENT ESTIMATORS. Since the deployed nowcast became
    # apply_nowcast_correction() (step 4), `zone_week_nc` IS the deterministic arm — labelling
    # it "epinowcast" made two of the three arms bit-identical and the sensitivity vacuous
    # (the previous run's bayes_nowcast_sensitivity.csv already showed the raw and "epinowcast"
    # beta0 draws agreeing to the last bit). Fit epinowcast for real here; this is the one place
    # that exercises the 04b completeness-ratio fix, and it is exactly what the arm is for.
    .zw_raw <- dat$zone_week; .zw_raw$confirmed_nc <- .zw_raw$confirmed
    .zw_variants <- list(raw = .zw_raw, deterministic = zone_week_nc)
    .zw_enw <- tryCatch(
      nowcast_zone_week_epinowcast(zone_week = dat$zone_week, linelist = dat$ll,
                                   analysis_date = ANALYSIS_DATE, outbreak_start = OUTBREAK_START),
      error = function(e) { message("[bayes-sens] epinowcast arm unavailable (non-fatal): ",
                                    conditionMessage(e)); NULL })
    if (!is.null(.zw_enw) &&
        identical(attr(.zw_enw, "nowcast_method"), "epinowcast"))
      .zw_variants[["epinowcast"]] <- .zw_enw
    else
      message("[bayes-sens] epinowcast arm omitted: the fit fell back to the deterministic ",
              "correction, which is already the `deterministic` arm.")
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

# Task 4 — evaluation over time: parameter estimates over folds, and per-fold
# predicted-vs-observed invasion for the featured models.
if (!is.null(lfo_results)) {
  # THE beta_0 TRACE RUNS PAST THE LAST SCORED ROUND, on purpose. This refits the featured
  # model at each cutoff and reads its posterior import coefficient; it needs TRAINING data,
  # not an outcome, so it is not bound by the cross-validation's fold window — which stops at
  # the last round whose outcome window closes by the analysis date. Tracing only the scored
  # rounds left the published beta_0 series ending a week before the data do, with no reason
  # a reader could see. Every week from the first scored round to the last COMPLETE week of
  # the grid is traced instead.
  #
  # "Complete" = the week whose last day is on or before ANALYSIS_DATE. A partial final week
  # would be fitted on a few days of reporting and its beta_0 would be an artefact of the
  # truncation, not a level.
  .fold_cutoffs <- local({
    cv <- sort(unique(as.Date(lfo_results$cutoff)))
    wk <- sort(unique(as.Date(zone_week_outbreak$week_start)))
    wk <- wk[wk >= min(cv) & (wk + 6L) <= as.Date(ANALYSIS_DATE)]
    out <- sort(unique(c(cv, wk)))
    if (length(out) > length(cv))
      message(sprintf("[over-time] beta_0 traced at %d cutoffs (%d scored rounds + %d later week(s) up to %s)",
                      length(out), length(cv), length(out) - length(cv),
                      format(max(out))))
    out
  })
  # The FREQUENTIST params-over-time trace (compute_params_over_time: a Firth cloglog GLM
  # refitted at every fold cutoff) is gone. Its figure, "params_over_time", is not in
  # FIGURE_KEEP, so the refits ran every time and the output was then gated out of the
  # published tree — pure cost. The RETAINED trace is the Bayesian one below
  # (bayes_params_over_time), which refits the featured model's own kernel.
  for (hh in LFO_HORIZONS) {
    vwrap(plot_predobs_over_folds(lfo_results, method = best_method, horizon = hh, save = TRUE))
    if (!is.null(best_bayes_method) && best_bayes_method %in% lfo_results$method)
      vwrap(plot_predobs_over_folds(lfo_results, method = best_bayes_method, horizon = hh, save = TRUE,
              file = "bayes_predicted_vs_observed_over_folds"))
  }

  # Task 6 — intuitive detection-vs-budget curve + balanced skill metrics for the
  # imbalanced invasion task (sensitivity at a fixed alert budget; balanced
  # accuracy / F1 / MCC at the Youden-optimal threshold).
  .skill_methods <- intersect(unique(c(best_method, best_bayes_method)),
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
  # also overlays the THREE structural baselines (built above; lfo_fig3) — gravity,
  # Flowminder epicentre inflow, and road travel time — so the model's prioritisation is
  # compared against three different notions of pure epicentre connectivity.
  .fig3_lfo     <- if (exists("lfo_fig3")) lfo_fig3 else lfo_results
  .fig3_methods <- unique(c(.skill_methods, if (exists(".bl_ok")) .bl_ok else character(0)))
  for (hh in LFO_HORIZONS)
    # Support from lfo_results, NOT .fig3_lfo -- see the bayes_model_performance_figure3 call
    # above. .fig3_lfo is lfo_fig3 whenever the appended baselines exist, and deriving the
    # support from that copy would put this panel's k / n_atrisk reference on a different
    # denominator from the published detection curve.
    vwrap(plot_paper_figure3(.fig3_lfo, .fig3_methods, horizon = hh, save = TRUE,
            support_cells = tryCatch(invasion_common_cells(lfo_results, hh),
                                     error = function(e) NULL)))
  balance_tbl <- tryCatch(
    purrr::map_dfr(LFO_HORIZONS, function(hh)
      purrr::map_dfr(unique(lfo_results$method),
                     function(m) invasion_balance_metrics(lfo_results, m, hh))),
    error = function(e) NULL)
  if (!is.null(balance_tbl) && nrow(balance_tbl))
    readr::write_csv(balance_tbl, file.path(OUT_DIAGNOSTICS, "invasion_balanced_skill.csv"))

  # REMOVED: plot_reporting_rate_map(). It mapped "the per-zone RELATIVE reporting-rate proxy
  # THE MODEL USES to up-weight under-ascertained source zones". No model uses it — the
  # reporting-rate structure (report_rate_vec, Renewal-M8-report) went with the frequentist arm
  # and ascertainment is not modelled at all. A published figure asserting a mechanism the
  # pipeline does not have is worse than no figure; "reporting_rate" is not in FIGURE_KEEP
  # either, so it was gated out of the tree while still being rendered every run.
  #
  # REMOVED: compute_beta_over_folds() + its figure, for the same reason as the params trace
  # above — "beta_and_completeness_over_folds" is not retained, so the per-fold GLM refits
  # were discarded. The Bayesian beta trace below IS retained.

  # Bayesian analogues of params-over-time + beta-over-folds: refit the featured
  # Bayesian model's mobility kernel (+ geo covariates, so covariate-HR traces exist)
  # at each fold cutoff and plot the POSTERIOR covariate HRs and import coefficient
  # beta0 over folds, with proper 90% CrI (bayes_params_over_time / bayes_beta_over_folds).
  if (!is.null(best_bayes_method) && requireNamespace("brms", quietly = TRUE)) {
    # Keep every kernel token (-dist, -fill, -split) so a featured kernel such as
    # Bayes-M13-dist-fill-geo is refit as M13-dist-fill, not the travel-time M13 — the
    # over-time trace must be faithful to the featured model's actual mobility kernel.
    # The local regex used here dropped -fill/-split, which is ALWAYS wrong under the
    # source-cell-fill default; mobility_kernel_from_method() (00_config.R) is the shared
    # parser used by the cascade too.
    .bmob <- mobility_kernel_from_method(best_bayes_method)
    if (is.na(.bmob) || !.bmob %in% names(mobility_matrices)) {
      .fallback <- if (MOBILITY_PRIMARY %in% names(mobility_matrices)) MOBILITY_PRIMARY
                   else names(mobility_matrices)[1]
      if (!is.na(.bmob))
        warning(sprintf("params-over-time: kernel %s not built; tracing %s instead.",
                        .bmob, .fallback), call. = FALSE)
      .bmob <- .fallback
    }
    bpot <- tryCatch(compute_bayes_params_over_time(zone_week_outbreak, .fold_cutoffs,
              mobility_matrices, gt_pmfs, covariates, osrm_mat, zones_all,
              # ONE definition, from 00_config.R. This was a hard-coded copy of the retired
              # 3-covariate geo set, so after the 2026-09-21 change the published beta0 trace
              # would have described a covariate set no scored model used.
              mob = .bmob, gt = "medium", cov_spec = get0("BAYES_GEO_COVARIATES",
                                            ifnotfound = c("ccvi", "d_min")),
              iter = 600L, nowcast_fn = apply_nowcast_correction, linelist = dat$ll,
              delay = .trunc_asof),
              error = function(e) { warning("bayes params-over-time: ", e$message); NULL })
    if (!is.null(bpot)) {
      # PUBLISH THE NUMBERS BEHIND THE RETAINED PANELS. bayes_params_over_time.pdf and
      # bayes_beta_over_folds.pdf are key_outputs deliverables carrying posterior hazard
      # ratios, an import coefficient and 90% CrIs — and no CSV existed behind either, so a
      # reader could not check a single value on them. The figures are drawn from exactly
      # these tables.
      tryCatch({
        .kd <- file.path(OUT_DIR, "key_outputs")
        if (!is.null(bpot$params) && nrow(bpot$params))
          readr::write_csv(bpot$params %>%
              dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 12))),
            file.path(.kd, "bayes_params_over_time.csv"))
        if (!is.null(bpot$beta) && nrow(bpot$beta))
          readr::write_csv(bpot$beta %>%
              dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 12))),
            file.path(.kd, "bayes_beta_over_folds.csv"))
        message("[over-time] published bayes_params_over_time.csv / bayes_beta_over_folds.csv")
      }, error = function(e)
        warning("[over-time] over-folds CSVs not written: ", conditionMessage(e), call. = FALSE))
      # show_title = FALSE: published panel, caption in the manuscript text.
      if (!is.null(bpot$params)) vwrap(plot_params_over_time(bpot$params, save = TRUE,
              file = "bayes_params_over_time", model_label = sprintf("Bayesian: %s + geo", .bmob),
              show_title = FALSE))
      # Title-free, larger-font, PDF + 600-dpi PNG: this panel is used as a
      # standalone manuscript figure, so its caption lives in the text.
      if (!is.null(bpot$beta)) vwrap(plot_beta_over_folds(bpot$beta, save = TRUE,
              file = "bayes_beta_over_folds",
              model_label = sprintf("Bayesian: %s", best_bayes_method), ci_label = "90% CrI",
              show_title = FALSE, base_size = 16, png = TRUE))
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
  risk_scores = get0("bayes_risk_scores"),
  lfo_results = lfo_results,
  zone_week = zone_week_nc, training_cutoff = training_cutoff,
  best_method = best_method,
  primary_method = NULL,
  n_models = if (!is.null(lfo_results)) dplyr::n_distinct(lfo_results$method) else 0L,
  best_bayes_method = best_bayes_method,
  bayes_weights = if (!is.null(bayes_suite)) bayes_suite$weights else NULL))

# Q2/Q3 — document HOW models are selected and WHAT goes into the best model
# (structure, mobility kernel, generation time, observation process, covariates,
# calibration, nowcast). Writes model_specification.md and appends to the report.
vwrap(write_model_details_report(best_method = best_method,
        primary_method = NULL,
        eval_tbl = eval_tbl,
        bayes_params = if (!is.null(bayes_suite)) bayes_suite$params else NULL,
        best_bayes_method = best_bayes_method,
        freq_beta0 = NA_real_, freq_cov = NULL))

# ---------------------------------------------------------------------------
# Key outputs — gather the headline Bayesian deliverables (both horizons) into a single
# key_outputs/ folder for quick sharing, once everything above has been generated. Each entry
# is the h1 file; its h2 analogue (…_h1 -> …_h2) is copied too when present.
# ---------------------------------------------------------------------------
local({
  key_dir <- file.path(OUT_DIR, "key_outputs")
  dir.create(key_dir, showWarnings = FALSE, recursive = TRUE)
  # THE MANIFEST MUST AGREE WITH FIGURE_KEEP (00_config.R). The 2026-09-17 streamlining gated
  # every save helper on figure_is_kept(), but this list was not updated, so it asked for SEVEN
  # per-horizon stems that can no longer be written (bayes_invasion_uncertainty_map_national,
  # bayes_predicted_vs_observed_over_folds, topk_precision, bayes_topk_precision,
  # bayes_invasion_uncertainty, bayes_model_performance_figure3, priority_scatter) — 14 names
  # that produced a permanent "not yet generated (skipped)" line every run — while FOUR figures
  # that ARE retained and ARE produced never reached key_outputs at all, including rt_national.pdf,
  # the very figure RUN_RT_ESTIMATION was flipped to TRUE to obtain.
  #
  # The `wanted` set is now FILTERED through figure_is_kept() below, so the two can never drift
  # apart again: adding a stem here without adding it to FIGURE_KEEP is a no-op, not a phantom.
  h1_files <- c(
    "maps/bayes_invasion_rank_map_national_h1.pdf",
    "maps/bayes_prob_vuln_choropleth_national_h1.pdf",
    "diagnostics/bayes_discrimination_summary_h1.pdf",
    "diagnostics/bayes_lfo_forecast_vs_outcome_h1.pdf",
    "reports/bayes_priority_scatter_h1.pdf")            # preparedness-priority scatter (Bayesian featured)
  no_horizon <- c("reports/bayes_risk_scores_all_zones.csv",
                  "reports/harmonised_confirmed_cases.csv",
                  "diagnostics/invasion_recalibration.csv",
                  "diagnostics/invasion_recalibration.json",
                  # Retained figures that this block previously never gathered.
                  "diagnostics/rt_national.pdf",
                  "diagnostics/skill_over_time_auc_pr_skill.pdf",
                  "diagnostics/bayes_beta_over_folds.pdf",
                  "diagnostics/bayes_params_over_time.pdf")
  wanted <- unique(c(as.vector(rbind(h1_files, sub("_h1\\.", "_h2.", h1_files))), no_horizon))
  # Drop anything the figure gate would refuse, so the manifest cannot ask for a file the
  # pipeline is no longer allowed to write. Non-figures (csv/json) always pass the gate.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk)) {
    .drop <- wanted[!vapply(wanted, .fk, logical(1))]
    if (length(.drop))
      message(sprintf("[key_outputs] %d manifest entr%s excluded by FIGURE_KEEP: %s",
                      length(.drop), if (length(.drop) == 1L) "y is" else "ies are",
                      paste(basename(.drop), collapse = ", ")))
    wanted <- wanted[vapply(wanted, .fk, logical(1))]
  }
  # FRESHNESS CHECK. This loop used to copy any `wanted` file that merely EXISTED. Every
  # producer upstream is wrapped in vwrap()/tryCatch, so a silently failed stage leaves the
  # PREVIOUS run's file in place and the manifest then promoted it into key_outputs/ — which
  # ci/collect_outputs.sh publishes wholesale. The result is last week's artifact shipped under
  # this week's bundle with a fresh copy mtime and nothing anywhere saying so. rt_national.pdf
  # is the clearest case: run_all.R prints "RUN_RT_ESTIMATION = FALSE ... will NOT be produced"
  # and then the old copy was published regardless.
  #
  # A file counts as fresh if it was written at or after this run started (.PH_T0, set at the
  # top of this script). The 2 s slack absorbs filesystem timestamp granularity only.
  .t_run <- get0(".PH_T0", ifnotfound = NULL)
  copied <- 0L; missing <- character(0); stale <- character(0)
  for (rel in wanted) {
    src <- file.path(OUT_DIR, rel)
    dst <- file.path(key_dir, basename(rel))
    if (!file.exists(src)) { missing <- c(missing, basename(rel)); next }
    .fresh <- is.null(.t_run) ||
              isTRUE(as.numeric(difftime(file.mtime(src), .t_run, units = "secs")) >= -2)
    if (.fresh) {
      file.copy(src, dst, overwrite = TRUE); copied <- copied + 1L
    } else {
      # Not produced by this run. Do NOT publish it, and clear any copy an earlier run left
      # in key_outputs/ so the bundle cannot present it as current.
      stale <- c(stale, basename(rel))
      if (file.exists(dst)) suppressWarnings(file.remove(dst))
    }
  }
  message(sprintf("[key_outputs] copied %d files -> %s", copied, key_dir))
  if (length(missing)) message("[key_outputs] not yet generated (skipped): ", paste(missing, collapse = ", "))
  if (length(stale))
    warning(sprintf(paste0("[key_outputs] %d manifest entr%s NOT produced by this run and were ",
                           "excluded (any stale copy in key_outputs/ was removed): %s"),
                    length(stale), if (length(stale) == 1L) "y was" else "ies were",
                    paste(stale, collapse = ", ")), call. = FALSE, immediate. = TRUE)
  # Raw twins keep the SAME basenames one level down, so a file under key_outputs/raw/ is
  # the same artifact on the uncorrected scale. Only the tables are gathered here; the raw
  # figure twins are written straight into their own raw/ directories by the second suite
  # pass, which is why they are not listed above.
  raw_dir <- file.path(key_dir, "raw")
  # EVERY horizon, not just h=1. The reliability / count-calibration writers above loop over
  # LFO_HORIZONS, so the h=2 raw twins are produced but were never gathered here — the primary
  # (recalibrated) list two blocks up already derives its h=2 entries, so the two lists
  # disagreed and key_outputs/raw/ silently held only half the scale comparison.
  raw_src <- c("reports/raw/bayes_risk_scores_all_zones.csv",
               "reports/raw/bayes_ensemble_risk_scores_all_zones.csv",
               sprintf("diagnostics/raw/forecast_reliability_h%d.csv", LFO_HORIZONS),
               sprintf("diagnostics/raw/forecast_count_calibration_h%d.csv", LFO_HORIZONS))
  n_raw <- 0L
  for (rel in raw_src) {
    src <- file.path(OUT_DIR, rel)
    if (file.exists(src)) {
      if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
      file.copy(src, file.path(raw_dir, basename(src)), overwrite = TRUE); n_raw <- n_raw + 1L
    }
  }
  message(sprintf("[key_outputs] copied %d raw twin(s) -> %s", n_raw, raw_dir))
})

# ---------------------------------------------------------------------------
# 15. Publication + manuscript figures, 3-month cascade, and Bayesian report refresh
# ---------------------------------------------------------------------------
# DELIBERATELY EMPTY. These standalone tools are launched ONCE, by .run_downstream()
# at the very END of this script. A second launcher used to sit here and ran the SAME
# seven scripts — make_publication_figures.R, make_manuscript_figures.R,
# make_si_model_figures.R, make_si_data_figures.R,
# make_topk15_ever.R, update_bayesian_report.R and run_cascade.R — so every pipeline
# run executed them twice, the cascade included (hours at CASCADE_N_MC = 1000).
# Worse than the wasted time: this position is BEFORE write_model_selection.R, so the
# first launch read the PREVIOUS run's model_selection.json — and 30_projection_config.R
# takes CASCADE_KERNEL from exactly that file, so the first cascade could simulate with a
# superseded mobility kernel. The surviving launcher runs after the selection artifact is
# written, adds the cascade figure stages, logs each child to outputs/logs/, and pins the
# child working directory to ROOT. See "DOWNSTREAM SUITES" at the foot of this file.

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
message(sprintf("Models run:         %d cross-validated (LFO)",
                if (!is.null(lfo_results) && "method" %in% names(lfo_results))
                  dplyr::n_distinct(lfo_results$method) else 0L))
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

# ---------------------------------------------------------------------------
# DOWNSTREAM SUITES — cascade analyses + publication/manuscript figures
# ---------------------------------------------------------------------------
# These modules used to be orphaned: nothing in run_all.R referenced 37, 40-43, the
# four make_*.R figure suites, or update_bayesian_report.R, and the cascade driver
# run_cascade.R had to be invoked by hand. They are now scheduled here.
#
# WHY CHILD PROCESSES rather than source(): every one of these is written as a
# STANDALONE Rscript. They source 00_config.R themselves and define their own
# top-level globals — `theme_pub`, `%||%`, `ST_DIR`, `OUTBREAK_START`, `IMP_DELAY`,
# `INK`/`MUTED`/`OKABE`, `HERE` — several of which collide with names this orchestrator
# is still holding when they would run. source()ing them into this session would
# silently rebind the orchestrator's own state (and each other's) in file order. A
# child Rscript per script gives exact isolation, and a crash in a figure cannot take
# the pipeline's outputs down with it.
#
# ORDERING is a genuine data dependency, not a preference:
#   A  read THIS run's outputs (bayes_risk_scores_all_zones.csv, lfo_results.rds, ...)
#   B  run_cascade.R -> modules 30-36, 38, 39; needs A's mobility_*.rds on disk
#   C  read the CASCADE outputs B just wrote (cascade_reach_scores_all_zones.csv, ...)
#   D  the k_indiv sweep (41) and its figure (42); 41 re-simulates the whole cascade
#      across a grid of overdispersion values, so it is OFF unless asked for.
# Each stage is independently switchable, each script is non-fatal, and every failure
# is collected and reported in one place at the end rather than scrolling past.
.run_downstream <- function() {
  # The SAME parser as the rest of the file (see .env_flag above), not a fourth private copy.
  .flag <- function(nm, default) .env_flag(nm, default)
  if (!.flag("RUN_DOWNSTREAM_SUITES", TRUE)) {
    message("\n[downstream] RUN_DOWNSTREAM_SUITES=0 — skipping cascade + figure suites.")
    return(invisible(NULL))
  }
  log_dir <- file.path(OUT_DIR, "logs"); dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  # Children inherit this process's working directory. Both figure suites now anchor on
  # here::here() and no longer care, but the child scripts are also run by hand from the
  # repo root, so ROOT is kept as the one documented working directory. Restored on exit.
  .old_wd <- setwd(ROOT); on.exit(setwd(.old_wd), add = TRUE)
  rscript <- file.path(R.home("bin"),
                       if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  results <- list()

  run_one <- function(script, label = script, dir = ST_DIR, extra_env = character()) {
    path <- file.path(dir, script)
    if (!file.exists(path)) {
      message(sprintf("[downstream] %-38s MISSING (%s)", label, path))
      results[[label]] <<- "missing"; return(invisible(NULL))
    }
    logf <- file.path(log_dir, sprintf("downstream_%s%s.log", sub("\\.R$", "", script),
                      if (any(grepl("^FORECAST_SCALE=raw$", extra_env))) "_raw" else ""))
    t0 <- Sys.time()
    # Child inherits this run's ANALYSIS_DATE explicitly: each script otherwise re-derives
    # it from latest.json, so a snapshot rotating mid-run would silently date the figures
    # differently from the forecasts they illustrate.
    # --no-save --no-restore, deliberately NOT --vanilla: --vanilla also implies
    # --no-init-file, which skips .Rprofile and therefore skips renv activation on any
    # branch or CI image that uses it. The child would then start with the wrong library paths and fail to find
    # packages the parent can see. These two flags give a clean, non-interactive session
    # while leaving .Rprofile/.Renviron — and thus library resolution — intact.
    st <- tryCatch(
      system2(rscript, c("--no-save", "--no-restore", shQuote(path)),
              stdout = logf, stderr = logf,
              env = c(sprintf("ANALYSIS_DATE=%s", format(ANALYSIS_DATE)),
                      sprintf("CASCADE_ST_DIR=%s", ST_DIR), extra_env)),
      error = function(e) { message("[downstream] ", label, " could not start: ",
                                    conditionMessage(e)); 127L })
    dt <- as.numeric(Sys.time() - t0, units = "secs")
    ok <- identical(as.integer(st), 0L)
    message(sprintf("[downstream] %-38s %-4s (%5.0fs)  log: %s",
                    label, if (ok) "OK" else sprintf("FAIL:%s", st), dt, basename(logf)))
    results[[label]] <<- if (ok) "ok" else paste0("fail:", st)
    invisible(NULL)
  }

  # ---- A. Figure/report suites that read THIS run's outputs -----------------
  # The figure suites are built TWICE: once on the primary (recalibrated) scale and once
  # on the raw scale into a raw/ sibling with the same basenames. `FS_FIGURE_SUITES` is the
  # set that takes the switch; 43_spread_kinematics.R is descriptive (observed arrivals, no
  # forecast probabilities) and update_bayesian_report.R is a report, so neither is rebuilt.
  FS_FIGURE_SUITES <- c("make_publication_figures.R", "make_manuscript_figures.R",
                        "make_topk15_ever.R")
  if (.flag("RUN_FIGURE_SUITES", TRUE)) {
    message("\n=== Downstream A: publication + manuscript figure suites (recalibrated) ===")
    # 43 FIRST. It publishes arrival_predictors.csv and arrival_predictor_fits.csv, which
    # Figure 1C and Figure S1 are drawn from — those panels no longer compute their own
    # R^2 / rho / n. Run after the figure suites, as this used to, the first run of a fresh
    # checkout fails outright and every later run draws Panel C from the PREVIOUS run's
    # arrival table: a main-text panel silently one run stale. 43 is descriptive (line list
    # + geography + mobility only) and reads nothing the figure suites write, so it is safe
    # at the head of this block.
    run_one("43_spread_kinematics.R")
    for (f in FS_FIGURE_SUITES) run_one(f)
    # Figure_spread_kinematics_compact is a RETAINED deliverable but was orphaned:
    # nothing in the pipeline ever invoked make_spread_kinematics_compact.R, so the
    # figure could only ever appear if someone ran it by hand. It reads 43's
    # spread_kinematics_{weekly,summary}.csv, so it MUST follow 43 in this order.
    run_one("make_spread_kinematics_compact.R")
    # Supplementary figure suites. They read only saved artefacts of THIS run — the
    # cross-validation, the convergence diagnostics, the beta_t trajectories, the mobility
    # matrices and both surveillance streams — so they follow the main suites and, like them,
    # must not run before the run that writes those files.
    run_one("make_si_model_figures.R")
    run_one("make_si_data_figures.R")
    run_one("update_bayesian_report.R")
    # DEFAULT FALSE since 2026-09-17. The raw twin re-ran all three figure suites to write a
    # second copy of every figure under raw/ siblings. None of those twins is on the retained
    # deliverable list, and figure_is_kept() now rejects any path with a raw/ component — so
    # leaving this on would spend a full second pass building figures that are then discarded.
    # Set RUN_RAW_FIGURE_TWIN=1 (with FIGURE_KEEP_ALL=1) to rebuild the raw-scale comparison.
    if (.flag("RUN_RAW_FIGURE_TWIN", FALSE)) {
      message("\n=== Downstream A2: the same suites on the RAW scale (raw/ siblings) ===")
      for (f in FS_FIGURE_SUITES)
        run_one(f, label = paste0(f, " [raw]"), extra_env = "FORECAST_SCALE=raw")
    }
  }

  # ---- B. The 3-month invasion cascade (modules 30-36, 38, 39) -------------
  # CASCADE_SMOKE=1 gives the fast path (M=40, heavy stages off) for CI.
  # SKIP_CASCADE is the opt-out documented in the header and in README; it must work
  # here, since this is now the only place the cascade is launched.
  .skip_cascade <- tolower(trimws(Sys.getenv("SKIP_CASCADE", ""))) %in% c("1", "true", "t", "yes", "y")
  if (.flag("RUN_CASCADE_SUITE", TRUE) && !.skip_cascade) {
    message("\n=== Downstream B: 3-month invasion cascade (run_cascade.R) ===")
    run_one("run_cascade.R")
  } else {
    message("\n[downstream] cascade skipped (",
            if (.skip_cascade) "SKIP_CASCADE set" else "RUN_CASCADE_SUITE=0", ").")
  }

  # ---- C. Figures that read the CASCADE outputs written by B ---------------
  if (.flag("RUN_CASCADE_FIGURES", TRUE)) {
    message("\n=== Downstream C: cascade figure suites ===")
    run_one("37_cascade_figure4.R")
    run_one("40_cascade_next_dominoes.R")
    run_one("make_manuscript_figure2_cascade.R")
  }

  # ---- D. k_indiv overdispersion sweep (expensive) -------------------------
  # 41 re-simulates the cascade across the whole k grid. OFF by default; 42 skips
  # cleanly when 41's outputs are absent, so requesting the figure alone is safe.
  if (.flag("RUN_KINDIV_SWEEP", FALSE)) {
    message("\n=== Downstream D: k_indiv overdispersion sweep ===")
    run_one("41_kindiv_sweep.R")
  } else {
    message("\n[downstream] RUN_KINDIV_SWEEP=0 (default) — skipping the k_indiv sweep (41).")
  }
  if (.flag("RUN_CASCADE_FIGURES", TRUE)) run_one("42_kindiv_sweep_figure.R")

  # ---- Scale manifest -------------------------------------------------------
  # Every primary artifact that has a raw twin, with a flag for whether the two actually
  # differ. Most do not: ranks are exactly invariant under the recalibration, so capture
  # curves, top-K panels and rank maps come out byte-comparable, and only the panels and
  # columns that carry probability LEVELS move. Recording that explicitly stops a reader
  # inferring a difference from the mere existence of two folders.
  tryCatch({
    key_dir <- file.path(OUT_DIR, "key_outputs")
    roots <- c(file.path(key_dir, "manuscript_figures"), file.path(key_dir, "manuscript_figures", "panels"),
               file.path(key_dir, "figures"), file.path(key_dir, "figures", "panels"), key_dir)
    rows <- list()
    for (d in roots) {
      rd <- file.path(d, "raw")
      if (!dir.exists(rd)) next
      for (f in list.files(rd, full.names = FALSE, recursive = FALSE)) {
        praw <- file.path(rd, f); pprim <- file.path(d, f)
        if (dir.exists(praw)) next
        # PDF writers embed a /CreationDate, so two byte-different PDFs can have identical
        # content: the comparison is only meaningful for formats without a timestamp.
        # Figures are written as a PDF+PNG pair, so the PNG twin carries the verdict and
        # the PDF row is left NA rather than asserting a difference that may not exist.
        .ext <- tolower(tools::file_ext(f))
        # A missing primary is NOT COMPARABLE (NA), not "differs": file.exists() && ... yields
        # FALSE, which the summary then counted under "differing".
        same <- if (.ext == "pdf") NA
                else if (!file.exists(pprim)) NA
                else identical(tools::md5sum(pprim)[[1]], tools::md5sum(praw)[[1]])
        # fixed=TRUE: OUT_DIR is a filesystem path, not a regex, and a "." in it would
        # otherwise match any character.
        .rel <- function(x) sub("^/", "", sub(OUT_DIR, "", x, fixed = TRUE))
        rows[[length(rows) + 1L]] <- data.frame(
          artifact = f, primary = .rel(pprim), raw_twin = .rel(praw),
          primary_exists = file.exists(pprim), identical_bytes = same,
          stringsAsFactors = FALSE)
      }
    }
    man <- if (length(rows)) dplyr::bind_rows(rows) else
      data.frame(artifact = character(0), primary = character(0), raw_twin = character(0),
                 primary_exists = logical(0), identical_bytes = logical(0))
    man$primary_scale <- "recalibrated (prequential for cross-validated panels; pooled/deployed for the current forecast)"
    man$raw_scale     <- "raw (uncorrected model output)"
    readr::write_csv(man, file.path(key_dir, "forecast_scale_manifest.csv"))
    message(sprintf("[scale] manifest: %d twin(s) recorded — %d identical, %d differing, %d not comparable (PDF) -> %s",
                    nrow(man), sum(man$identical_bytes %in% TRUE),
                    sum(man$identical_bytes %in% FALSE), sum(is.na(man$identical_bytes)),
                    file.path(basename(key_dir), "forecast_scale_manifest.csv")))
  }, error = function(e) message("[scale] manifest skipped: ", conditionMessage(e)))

  bad <- names(results)[!vapply(results, identical, logical(1), "ok")]
  message("\n[downstream] ", length(results) - length(bad), "/", length(results),
          " script(s) succeeded",
          if (length(bad)) paste0("; NOT OK: ", paste(bad, collapse = ", "),
                                  " (see ", log_dir, ")") else "", ".")
  invisible(results)
}
# Non-fatal by construction: the modelling outputs are already written and committed
# above, so no downstream figure can invalidate them.
downstream_results <- tryCatch(.run_downstream(),
  error = function(e) { message("[downstream] suite aborted: ", conditionMessage(e)); NULL })

# ---- Archive figures the retained-figure gate no longer produces ------------
# figure_is_kept() stops a suppressed figure from being WRITTEN, but it cannot remove one that
# an earlier run (before the allow-list) already left on disk. Those files keep current-looking
# mtimes in the live output tree and are indistinguishable from this run's deliverables to
# anyone browsing the folder. Move them — never delete — into a dated archive outside the
# published tree, so the allow-list is what the output directories actually contain.
tryCatch({
  .fig_dirs <- unique(c(file.path(OUT_DIR, "key_outputs", "figures"),
                        file.path(OUT_DIR, "key_outputs", "manuscript_figures"),
                        file.path(OUT_DIR, "cascade", "figures"),
                        OUT_MAPS, OUT_DIAGNOSTICS, file.path(OUT_DIAGNOSTICS, "delay_fits")))
  .fig_dirs <- .fig_dirs[dir.exists(.fig_dirs)]
  .cands <- unlist(lapply(.fig_dirs, function(d)
    list.files(d, pattern = "\\.(pdf|png)$", full.names = TRUE, recursive = TRUE)), use.names = FALSE)
  # Never touch anything already inside an archive, and ask the SAME gate the savers ask.
  .cands <- .cands[!grepl("_not_retained", .cands, fixed = TRUE)]
  .drop  <- .cands[!vapply(.cands, function(f) isTRUE(figure_is_kept(f)), logical(1))]
  if (length(.drop)) {
    .arch <- file.path(OUT_DIR, "_archive_not_retained",
                       format(Sys.Date(), "%Y%m%d"))
    .moved <- 0L
    for (f in .drop) {
      # substring(), not sub() with an interpolated pattern: OUT_DIR is a FILESYSTEM PATH, not
      # a regex, and a "." / "+" / "(" anywhere in the checkout path would otherwise match the
      # wrong characters and mis-strip the prefix, landing the archive at a nested absolute
      # path. (The scale-manifest block above makes the same point and uses fixed = TRUE.)
      .pfx <- paste0(OUT_DIR, .Platform$file.sep)
      rel  <- if (startsWith(f, .pfx)) substring(f, nchar(.pfx) + 1L) else basename(f)
      dst <- file.path(.arch, rel)
      dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
      if (isTRUE(suppressWarnings(file.rename(f, dst)))) .moved <- .moved + 1L
      else if (isTRUE(file.copy(f, dst, overwrite = TRUE))) {
        suppressWarnings(file.remove(f)); .moved <- .moved + 1L
      }
    }
    message(sprintf("[gate] archived %d/%d non-retained figure file(s) out of the published tree -> %s",
                    .moved, length(.drop), .arch))
  } else {
    message("[gate] no non-retained figure files left in the published tree.")
  }
}, error = function(e) message("[gate] archive sweep skipped: ", conditionMessage(e)))

invisible(list(
  fc_calibrated = calibrated_forecasts,
  fc_ensemble   = fc_ensemble,
  lfo_results   = lfo_results,
  metrics       = metrics_summary,
  calibration   = calibration_gain,
  mobility      = mobility_matrices,
  gt_pmfs       = gt_pmfs,
  Rt            = Rt_scalar,
  dat           = dat,
  downstream    = downstream_results
))

if (exists(".mark", mode = "function"))
  .mark(sprintf("DONE      run_all.R completed (%.1f min total)",
                as.numeric(difftime(Sys.time(), .PH_T0, units = "mins"))))
