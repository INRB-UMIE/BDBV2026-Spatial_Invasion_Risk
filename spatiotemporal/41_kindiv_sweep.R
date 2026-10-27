# =============================================================================
# 41_kindiv_sweep.R — COMPREHENSIVE individual-offspring overdispersion sweep
#
# Purpose: replace the 2-point k_indiv sensitivity (0.20, 0.40) with a dense grid
# spanning heavy superspreading (k=0.05) to the near-Poisson limit (k=4), densely
# sampling the Lloyd-Smith (2005) EBOV-plausible band [0.2, 0.4]. For each k we
# re-run the PRIMARY scenario cascade and record ranking stability (Spearman,
# Kendall, top-15 Jaccard) vs the production k (CASCADE_K_INDIV) plus absolute reach/establishment.
#
# FIX vs 33_cascade_eval.R::cascade_sensitivity: that routine compared k-perturbed
# runs (executed at the DEFAULT psi=CASCADE_PSI=1.0) against a baseline computed at
# the FITTED psi (=4.0) — conflating the k change with a psi change. Here every run
# (baseline AND alternatives) uses the identical (delta, psi) of the production calibration
# (psi is now fixed at 0), so the effect attributed to k is purely k.
#
# Usage:  Rscript spatiotemporal/41_kindiv_sweep.R
# Env:    CASCADE_SWEEP_M   MC iterations per k (default 2000, matches headline run)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(here) })
ST_DIR <- Sys.getenv("CASCADE_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
t0 <- Sys.time(); tick <- function(m) message(sprintf("[sweep t=%.0fs] %s",
                                     as.numeric(Sys.time() - t0, units = "secs"), m))

source(file.path(ST_DIR, "00_config.R"))
# ONE module list, from 00_config.R. Nine files used to hand-maintain their own and all nine
# differed — see load_pipeline_modules().
load_pipeline_modules("cascade",
                      extra = c("17_invasion_viz.R", "33_cascade_eval.R",
                                "33b_cascade_calibration.R", "34_conditional_queries.R",
                                "35_cascade_viz.R", "36_report.R"),
                      dir = ST_DIR)

SMOKE   <- identical(Sys.getenv("CASCADE_SWEEP_SMOKE"), "1")
M       <- as.integer(Sys.getenv("CASCADE_SWEEP_M", unset = if (SMOKE) "60" else "2000"))
FIT_ITER<- as.integer(Sys.getenv("CASCADE_SWEEP_ITER", unset = if (SMOKE) "400" else "2000"))
# The production sweep grid (30_projection_config.R), the one run_cascade.R's sweep uses.
K_GRID  <- if (SMOKE) c(0.20, 0.30, 0.40) else CASCADE_K_INDIV_SWEEP
K_BASE  <- CASCADE_K_INDIV                     # reference for ranking-stability metrics (the
                                               # production k; 42 skips a sweep built on another)
kern    <- CASCADE_KERNEL                      # whatever the LFO-CV selection picked
DIAG    <- file.path(OUT_CASCADE, if (SMOKE) "diagnostics_smoke" else "diagnostics")
if (!dir.exists(DIAG)) dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)
tick(sprintf("sourced; M=%d, kernel=%s, grid=[%s]", M, kern, paste(K_GRID, collapse=",")))

# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. run_all.R / run_cascade.R launch it as its own subprocess. Guarded so a source()
# — a test, an ad-hoc check, a sweep that loads "every module" — defines what it needs and
# stops, instead of silently rewriting files. KINDIV_FORCE_RUN=1 forces a run from a source().
if (!is_script_run("41_kindiv_sweep.R") && !identical(Sys.getenv("KINDIV_FORCE_RUN"), "1")) {
  message("[kindiv] sourced, not run: nothing is written. Use `Rscript 41_kindiv_sweep.R` to run it.")
} else {

# ---- reconstruct production setup (identical to run_cascade.R) --------------
# The default MUST be a layer carrying the line list (`ll`): the cascade's national R is the
# shared EpiNow2 posterior fitted from it, so against a layer without one cascade_reff() stops
# rather than quietly putting this sweep on a different R from the production run it claims to
# perturb. (The old requirement here was `nowcast_cv`, for the nowcast-uncertainty bootstrap;
# that bootstrap is retired — EpiNow2 models the incomplete recent tail itself.)
#
# ONLY THE CANONICAL LAYER, never CASCADE_LAYER_CACHE. run_cascade.R writes that optional build
# cache only when it REBUILDS the layer, before it stamps the run id, and never re-writes it when
# it loads it — so a layer read from there carries no cascade_run_id, the calibration-JSON check
# below could never match it, and this sweep silently recalibrated instead of using production's
# delta. The canonical file is re-written by every run_cascade.R run, after the run id is set,
# with the layer that run actually used.
CACHE <- file.path(OUT_CASCADE, "layer_cache_v2.rds")
if (!file.exists(CACHE))
  stop("[kindiv] no canonical layer cache at ", CACHE, ". Run run_cascade.R first (it writes ",
       "the layer it used there).", call. = FALSE)
layer <- readRDS(CACHE); tick("layer loaded from cache")
# PROVENANCE: the layer must be for THIS analysis date. run_cascade.R stamps the layer it used
# (cascade_key, cascade_run_id); one from another frame has a different week grid and would make
# this a sweep around a different model from the production run it claims to perturb.
.lk <- layer$cascade_key
if (is.null(.lk) || !identical(unname(.lk[["analysis_date"]]), format(as.Date(ANALYSIS_DATE))))
  stop(sprintf(paste0("[kindiv] the layer cache %s is for analysis date %s, not this run's %s. ",
                      "Re-run run_cascade.R (it rewrites the cache) before the sweep."),
               CACHE, if (is.null(.lk)) "UNKNOWN (unstamped cache)" else .lk[["analysis_date"]],
               format(as.Date(ANALYSIS_DATE))), call. = FALSE)
zones_all <- layer$zones_all; gt_pmfs <- layer$gt_pmfs
province_map <- load_province_map()
zone_province <- setNames(province_map$province, province_map$nom)
vuln <- compute_vulnerability_index(layer$covariates, zones_all, layer$osrm_mat)

stopifnot(kern %in% names(layer$mobility_matrices))
# Fit the primary kernel plus the sensitivity kernels (for the corrected kernel axis
# of cascade_sensitivity); skip the extra fits in smoke mode.
kernels <- unique(c(kern, if (!SMOKE) intersect(CASCADE_SENS_KERNELS, names(layer$mobility_matrices))))
designs <- list(); fits <- list()
# The validated model's rolling as-of design, exactly as run_cascade.R builds production's.
.trunc_asof <- cascade_asof_truncation(layer$ll)
for (kn in kernels) {
  designs[[kn]] <- cascade_invasion_design(layer$zone_week_nc, layer, kernel = kn, gt = CASCADE_GT,
                                           linelist = layer$ll, trunc_delay = .trunc_asof)
  # The production hazard's form (CASCADE_FORM), through the cascade's one fitting entry point.
  fits[[kn]]    <- cascade_fit_hazard(designs[[kn]], layer$zone_week_nc, layer, kernel = kn,
                                      gt = CASCADE_GT, linelist = layer$ll,
                                      trunc_delay = .trunc_asof, iter = FIT_ITER)
  tick(sprintf("fitted hazard: %s (%s)", kn, CASCADE_FORM))
}
if (is.null(fits[[kern]])) stop("[kindiv] the primary kernel's hazard did not fit.", call. = FALSE)
fits <- Filter(Negate(is.null), fits)
design <- designs[[kern]]; fit <- fits[[kern]]
# The sweep is a SENSITIVITY AROUND PRODUCTION, so it must use production's calibration.
# Prefer the parameters run_cascade.R just fitted and recorded; only re-fit if that record
# is missing. (It previously called the retired naive-pace cascade_fit_psi() here, so the
# sweep could be centred on a different psi from the run it claims to perturb.)
.calib_json <- file.path(OUT_CASCADE, "diagnostics", "cascade_calibration.json")
.cal <- if (file.exists(.calib_json))
  tryCatch(jsonlite::fromJSON(.calib_json), error = function(e) NULL) else NULL
# ONLY the calibration of the run that produced THIS layer: same cascade_run_id and analysis date.
# A JSON left by an earlier run (or a run whose calibration stage did not complete) is not this
# production run's calibration, and centring the sweep on it is the stale-file bug.
.cal_ok <- !is.null(.cal) && is.finite(.cal$delta %||% NA_real_) && is.finite(.cal$psi %||% NA_real_) &&
           identical(.cal$cascade_run_id %||% NA_character_, layer$cascade_run_id %||% "") &&
           identical(.cal$analysis_date %||% NA_character_, format(as.Date(ANALYSIS_DATE))) &&
           # and for the SAME hazard: form AND kernel (a calibration of fixed beta is not one of a
           # form, nor one of another kernel when the two launches' overrides differ)
           identical(.cal$form %||% "fixed", CASCADE_FORM) &&
           identical(.cal$kernel %||% NA_character_, kern)
if (.cal_ok) {
  delta <- .cal$delta; psi <- .cal$psi
  tick(sprintf("calibration read from %s (delta=%.3f, psi=%.2f, %s, run %s)", basename(.calib_json),
               delta, psi, .cal$calibration_status %||% "status unrecorded", .cal$cascade_run_id))
} else {
  message("[sweep] no cascade_calibration.json from the run that wrote this layer — recalibrating ",
          "(this costs one brms fit per origin).")
  .cc   <- cascade_calibrate(layer, zone_province = zone_province, kernel = kern, gt = CASCADE_GT)
  delta <- .cc$delta; psi <- .cc$psi
  if (is.null(.cc$delta_oos))
    warning("[sweep] the recalibration fitted NO multi-week delta; the sweep runs at 16b's ",
            "short-horizon factor, which is not production's calibrated level.", call. = FALSE)
}
# R_eff EXACTLY AS PRODUCTION BUILDS IT (run_cascade.R). The projection transmits on the national
# EpiNow2 R, which no delta enters; delta reaches only the conjugate R DIAGNOSTIC's import term
# (31_source_dynamics.R), where production passes 16b's short-horizon factor (delta_imp), not the
# K-week calibrated delta used for seeding. The same factor here keeps that diagnostic equal to
# production's; the sweep's premise is that k_indiv is the only thing that changes.
delta_imp <- tryCatch(cascade_fit_delta(cov_model = cascade_hazard_label(kern)),
                      error = function(e) NA_real_)
if (!is.finite(delta_imp) || delta_imp <= 0) {                # production's own fallback
  warning("[sweep] no recalibration factor available for the import term; using the RAW hazard ",
          "scale, as run_cascade.R does.", call. = FALSE)
  delta_imp <- 1
}
reff  <- cascade_reff(layer$zone_week_nc, zones_all, layer, fit, design,
                      delta = as.numeric(delta_imp), zone_province = zone_province,
                      gt = CASCADE_GT, kernel = kern, issue_date = ANALYSIS_DATE)
prep  <- cascade_prepare(layer, fit, design, reff, kernel = kern, gt = CASCADE_GT)

gate  <- cascade_consistency_gate(prep, fit, design, delta, layer$zone_week_nc,
           layer$mobility_matrices, gt_pmfs, layer$covariates, layer$osrm_mat,
           zones_all, n_mc = 3000L)
tick(sprintf("production setup: delta=%.3f  psi=%.2f  gate X2 p=%.3f (%d zones; rank rho=%.3f, n=%d)",
             as.numeric(delta), psi, gate$x2_p, gate$n_tested, gate$spearman_signal, gate$n_signal))

# ---- run the sweep via the shared canonical routine (33_cascade_eval.R) -----
# Every k — baseline AND alternatives — uses the fitted (delta, psi); the effect
# attributed to k is purely k (no psi conflation). See cascade_kindiv_sweep().
scP <- CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]
prov_vec <- setNames(province_map$province, province_map$nom)
res <- cascade_kindiv_sweep(prep, scP, delta = delta, psi = psi, pop_vec = layer$pop,
                            province = prov_vec, k_grid = K_GRID, k_base = K_BASE, n_mc = M)
sweep  <- res$summary
byzone <- res$byzone
print(as.data.frame(sweep), digits = 4)

# ---- write outputs ----------------------------------------------------------
readr::write_csv(sweep,  file.path(DIAG, "cascade_kindiv_sweep.csv"))
readr::write_csv(byzone, file.path(DIAG, "cascade_kindiv_sweep_byzone.csv"))
meta <- list(cascade_run_id = layer$cascade_run_id %||% NA_character_,
             analysis_date = format(as.Date(ANALYSIS_DATE)),
             calibrated = if (.cal_ok) isTRUE(.cal$calibrated) else !is.null(.cc$delta_oos),
             rolling_predictors = isTRUE(design$rolling_predictors),
             kernel = kern, form = CASCADE_FORM, M = M, k_grid = K_GRID, k_base = K_BASE,
             delta = as.numeric(delta), cal_in_large = attr(delta, "cal_in_large") %||% NA,
             psi = psi, gate_spearman = gate$spearman_signal, gate_n_signal = gate$n_signal,
             gate_x2_p = gate$x2_p, gate_pass = isTRUE(gate$pass),
             gate_level_ratio = gate$level_ratio, scenario = CASCADE_SCENARIO_PRIMARY,
             seed = CASCADE_SEED, generated = as.character(Sys.time()))
jsonlite::write_json(meta, file.path(DIAG, "cascade_kindiv_sweep_meta.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 6)

# ---- corrected kernel & psi sensitivity axes (matched-psi baseline) ---------
# Recompute cascade_sensitivity with the DEPLOYED psi so the kernel/psi axes are perturbed
# against a baseline at the SAME psi (fixes the psi-conflation bug). psi is no longer fitted —
# it is fixed at 0 — so the baseline is the no-saturation model and the psi axis measures what
# turning saturation ON would do, rather than perturbing around a fitted value. The k_indiv
# axis is owned by the comprehensive sweep above and no longer duplicated.
if (!SMOKE && length(fits) > 1L) {
  base_reach <- cascade_reach_table(
    simulate_cascade(prep, scP, n_mc = M, delta = delta, psi = psi, pop_vec = layer$pop))
  sens <- tryCatch(cascade_sensitivity(layer, designs, fits, reff, delta,
                     base_reach = base_reach, psi = psi),
                   error = function(e) { message("[sweep] sensitivity: ", conditionMessage(e)); NULL })
  if (!is.null(sens)) {
    readr::write_csv(sens, file.path(DIAG, "cascade_sensitivity.csv"))
    tick(sprintf("corrected sensitivity axes: %s",
                 paste(sprintf("%s=%.3f", sens$axis, sens$spearman), collapse = "; ")))
  }
}
tick(sprintf("DONE — wrote cascade_kindiv_sweep.csv (%d rows), byzone (%d rows), meta.json",
             nrow(sweep), nrow(byzone)))

}   # end of the script-run guard
