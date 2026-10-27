# =============================================================================
# run_cascade.R — 3-MONTH SPATIAL INVASION CASCADE: master driver
# BDBV 2026 DRC · runs modules 30-36 plus 38 (urban scenarios) and 39 (the cascade
#                 evaluation figure); see PLAN_3MONTH_INVASION.md
#
# Usage:  Rscript spatiotemporal/run_cascade.R
# Env flags:
#   CASCADE_N_MC          MC iterations (default 2000; smoke 40)
#   CASCADE_SMOKE=1       fast path (M=40, heavy stages off, calibration off)
#   CASCADE_LAYER_CACHE   optional rds path to cache/reuse the reconstructed layer (reused only
#                         when its key matches this run). That copy is written only when the
#                         layer is REBUILT, and before the run id is stamped, so it is a build
#                         cache, not a provenance-stamped product. The layer actually used is
#                         ALWAYS written, stamped with this run's id, to the canonical
#                         cascade/layer_cache_v2.rds — the ONLY layer file 41_kindiv_sweep.R and
#                         make_manuscript_figure2_cascade.R read.
#   CASCADE_CALIBRATE=0   run DELIBERATELY uncalibrated (delta = 16b's short-horizon factor).
#                         Without it a calibration that fails or fits nothing STOPS the run.
#   CASCADE_RUN_BACKTEST / CASCADE_RUN_SBC / CASCADE_RUN_SENS  (default 1; smoke 0)
#
# PROVENANCE. Every run gets a CASCADE_RUN_ID, written to diagnostics/cascade_run_stamp.json at
# the start (status "running") and at the end (status "complete"), and stamped into the
# validation CSV/JSON products and the canonical layer cache. Each stage DELETES its own previous
# outputs when it starts, so a stage that fails or is switched off leaves nothing behind for a
# later script to mistake for this run's. What the two downstream readers check:
#   * make_manuscript_figure2_cascade.R stops unless cascade_run_stamp.json has status "complete"
#     and this analysis date, and cascade_backtest_detail.csv and the canonical layer cache carry
#     that stamp's cascade_run_id.
#   * 41_kindiv_sweep.R does NOT read the stamp. It stops if the canonical layer's
#     cascade_key[["analysis_date"]] is missing or is not this analysis date, and it takes
#     delta/psi from cascade_calibration.json only when that JSON's cascade_run_id equals the
#     layer's and its analysis_date is this one — otherwise it recalibrates (one brms fit per
#     origin).
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(here) })
ST_DIR <- Sys.getenv("CASCADE_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
t0 <- Sys.time(); tick <- function(m) message(sprintf("[cascade t=%.0fs] %s",
                                     as.numeric(Sys.time() - t0, units = "secs"), m))

source(file.path(ST_DIR, "00_config.R"))
# ONE module list, from 00_config.R (PIPELINE_MODULES, through "cascade"), plus the cascade's own
# suite. This file used to hand-keep its list; it drifted when the time-varying forms' engine
# (21c_frailty_models.R) joined the model layer, so a featured form's hazard could not be fitted
# here. quiet = FALSE: 30_projection_config.R's warnings (form/kernel fallbacks) must be seen.
load_pipeline_modules("cascade",
                      extra = c("17_invasion_viz.R", "33_cascade_eval.R",
                                "33b_cascade_calibration.R", "34_conditional_queries.R",
                                "35_cascade_viz.R", "36_report.R", "38_urban_scenarios.R",
                                "39_cascade_eval_figure.R"),
                      dir = ST_DIR, quiet = FALSE)
SMOKE <- identical(Sys.getenv("CASCADE_SMOKE"), "1")

# ---- run identity, stamp, and stale-output hygiene -------------------------
CASCADE_RUN_ID <- sprintf("%s-%d", format(Sys.time(), "%Y%m%dT%H%M%S"), Sys.getpid())
.cd <- file.path(OUT_CASCADE, "diagnostics")
dir.create(.cd, recursive = TRUE, showWarnings = FALSE)
.stamp_path <- file.path(.cd, "cascade_run_stamp.json")
.write_stamp <- function(status) {
  jsonlite::write_json(list(cascade_run_id = CASCADE_RUN_ID,
                            analysis_date = format(as.Date(ANALYSIS_DATE)),
                            kernel = CASCADE_KERNEL, form = CASCADE_FORM, status = status,
                            written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
                       .stamp_path, auto_unbox = TRUE, pretty = TRUE)
}
.write_stamp("running")
# Provenance columns for every validation table this run writes. A reader compares them with
# cascade_run_stamp.json, so a file left by an earlier run can never pass as this run's.
.stamp_df <- function(d) {
  if (is.null(d) || !is.data.frame(d)) return(d)
  d$cascade_run_id <- CASCADE_RUN_ID
  d$analysis_date  <- format(as.Date(ANALYSIS_DATE))
  d
}
# Delete a stage's previous outputs at the moment the stage starts. A stage that then fails, or
# is switched off, leaves NOTHING behind — rather than last run's files, which later stages and
# figure scripts used to read as current (cascade_backtest_detail.csv from a run three weeks
# older than the reach tables beside it).
.reset_outputs <- function(paths, what) {
  paths <- paths[file.exists(paths)]
  if (length(paths)) {
    file.remove(paths)
    message(sprintf("[cascade] %s: removed %d previous output(s) before running: %s",
                    what, length(paths), paste(basename(paths), collapse = ", ")))
  }
  invisible(paths)
}
RUN_BT   <- !SMOKE && !identical(Sys.getenv("CASCADE_RUN_BACKTEST"), "0")
RUN_SBC  <- !SMOKE && !identical(Sys.getenv("CASCADE_RUN_SBC"), "0")
RUN_SENS <- !SMOKE && !identical(Sys.getenv("CASCADE_RUN_SENS"), "0")
RUN_URBAN <- !identical(Sys.getenv("CASCADE_RUN_URBAN"), "0")
tick(sprintf("sourced (M=%d, backtest=%s sbc=%s sens=%s)", CASCADE_N_MC, RUN_BT, RUN_SBC, RUN_SENS))
# A FRAILTY FORM IS NEVER FEATURED (TV_FORMS `selectable`, 00_config.R; since 2026-10-19). One that
# arrives from model_selection.json, not from the CASCADE_FORM override, comes from a selection run
# before that rule, so the main projection would be built on a model the selection can no longer
# pick. The override stays available for a deliberate sensitivity run (30 warns).
if (!nzchar(Sys.getenv("CASCADE_FORM")) && !identical(CASCADE_FORM, "fixed") &&
    !isTRUE(Filter(function(x) identical(x$form, CASCADE_FORM), TV_FORMS)[[1]]$selectable))
  stop(sprintf(paste0("[run_cascade] model_selection.json features the '%s' form, which is a secondary ",
                      "analysis and is never featured: that selection predates the rule. Re-run ",
                      "run_all.R first (or set CASCADE_FORM explicitly for a sensitivity run)."),
               CASCADE_FORM), call. = FALSE)

# ---- reconstruct modelling layer (cache option) ----------------------------
# THE CACHE IS KEYED. A layer is a function of the analysis date, the line-list snapshot, the
# onset-imputation seed/stratum and the data definition; one built for any other combination has
# a different week grid and different counts. It used to be loaded unconditionally whenever the
# file existed, so a cache from an earlier frame silently became this run's modelling layer.
.layer_key <- c(analysis_date = format(as.Date(ANALYSIS_DATE)),
                linelist_snapshot = tryCatch(rlang::hash(readLines(LINELIST_JSON, warn = FALSE)),
                                             error = function(e) NA_character_),
                random_seed = as.character(get0("RANDOM_SEED", ifnotfound = NA)),
                # The onset-imputation draw has its own seed (defaults to RANDOM_SEED, 00_config.R);
                # a layer built under another onset draw has different counts.
                onset_impute_seed = as.character(get0("ONSET_IMPUTE_SEED",
                                                      ifnotfound = get0("RANDOM_SEED", ifnotfound = NA))),
                onset_impute_stratum = as.character(get0("ONSET_IMPUTE_STRATUM",
                                                         ifnotfound = NULL) %||% "pooled"),
                data_definition_tag = as.character(get0("DATA_DEFINITION_TAG", ifnotfound = "")),
                # The CODE that builds the layer (data prep, delays and truncation, nowcast): a
                # change there (e.g. the truncation's support renormalisation) changes the layer
                # without changing any of the keys above.
                code = paste(unname(tools::md5sum(file.path(ST_DIR,
                  c("00_config.R", "01_data_prep.R", "02_epi_params.R", "03_mobility_matrices.R",
                    "04_nowcasting.R")))),
                  collapse = "|"),
                # The MOBILITY MATRICES the layer carries are read from OUT_MOBILITY: a kernel
                # rebuilt since the cache was written changes the layer under the same code.
                mobility = local({
                  .ks <- sort(unique(c(CASCADE_KERNEL, CASCADE_KERNEL_FALLBACK, CASCADE_SENS_KERNELS)),
                              method = "radix")
                  paste(unname(tools::md5sum(file.path(OUT_MOBILITY, sprintf("mobility_%s.rds", .ks)))),
                        collapse = "|")
                }))
CACHE <- Sys.getenv("CASCADE_LAYER_CACHE")
layer <- NULL
if (nzchar(CACHE) && file.exists(CACHE)) {
  layer <- readRDS(CACHE)
  if (!identical(layer$cascade_key, .layer_key)) {
    message(sprintf(paste0("[cascade] layer cache %s was built for a different frame (key %s vs this ",
                           "run's %s); REBUILDING it rather than reuse a stale layer."),
                    CACHE, paste(layer$cascade_key %||% "none", collapse = "|"),
                    paste(.layer_key, collapse = "|")))
    layer <- NULL
  } else tick("layer loaded from cache (key matches this run)")
}
if (is.null(layer)) {
  dat <- prep_all_data()
  gt_pmfs <- compute_all_gt_pmfs(); osrm_mat <- load_osrm()
  zones_all <- names(dat$pop)
  # TRAIN/DEPLOY CONSISTENCY (2026-09-17), matching run_all.R step 4. The cascade's own
  # calibration (33b_cascade_calibration.R) nowcasts each as-of origin with
  # apply_nowcast_correction(); building the deployed layer with epinowcast instead left the
  # cascade with exactly the estimator mismatch the rest of the suite just closed — delta
  # fitted under one correction and applied under another.
  # ...and, since 2026-09-21, with the SAME estimated truncation run_all.R's deployed nowcast
  # uses ("extract" regime: dat$zone_week is the current extract, truncated by onset ->
  # appearance-in-extract, not onset -> sample). Leaving this on the old onset->sample delay
  # would have re-opened the estimator mismatch from the other side — 33b's origins corrected
  # with a fitted truncation, this layer with a delay ~2 d shorter.
  .trunc_layer <- tryCatch(
    .trunc_as_delay_spec(epinow2_truncation(issue_date = NULL, ll = dat$ll), "extract"),
    error = function(e) { warning("[cascade] deployed truncation unavailable (", conditionMessage(e),
                                  "); falling back to the onset->sample delay, which UNDER-corrects ",
                                  "the recent weeks.", call. = FALSE); NULL })
  if (!is.null(.trunc_layer))
    message(sprintf("[cascade] deployed truncation: %s (mean %.2f d, sd %.2f d, %s)",
                    .trunc_layer$family, .trunc_layer$mean, .trunc_layer$sd, .trunc_layer$estimator))
  zw <- apply_nowcast_correction(dat$zone_week, analysis_date = ANALYSIS_DATE,
                                 delay = .trunc_layer)
  mob_dir <- file.path(OUT_MOBILITY)
  need <- unique(c(CASCADE_KERNEL, CASCADE_KERNEL_FALLBACK, CASCADE_SENS_KERNELS))
  # Load what is on disk. A MISSING SENSITIVITY kernel must not abort the run (the kernel
  # set moves with the mobility build, and the sensitivity arm is optional) — it is dropped
  # with a warning and the filter below sees a shorter list. A missing FEATURED kernel is
  # fatal: there is nothing to project with.
  .have <- file.exists(file.path(mob_dir, sprintf("mobility_%s.rds", need)))
  if (!all(.have))
    warning(sprintf("[cascade] mobility kernel(s) not built, skipping: %s (run run_all.R to rebuild)",
                    paste(need[!.have], collapse = ", ")), call. = FALSE)
  if (!file.exists(file.path(mob_dir, sprintf("mobility_%s.rds", CASCADE_KERNEL))))
    stop(sprintf("[cascade] featured kernel mobility_%s.rds is missing from %s; rerun run_all.R.",
                 CASCADE_KERNEL, mob_dir), call. = FALSE)
  need <- need[.have]
  mobility_matrices <- setNames(lapply(need, function(k)
    readRDS(file.path(mob_dir, sprintf("mobility_%s.rds", k)))), need)
  # `ll` is carried so cascade_backtest() can reconstruct each origin's TRAINING counts
  # as-of that origin (reaggregate_asof + per-fold nowcast), instead of slicing the final
  # revised counts — which would fold in cases only reported after the origin and make the
  # backtest optimistic. An older cached layer without `ll` still works: cascade_backtest()
  # warns and falls back to the leaky slice rather than failing.
  layer <- list(zone_week_nc = zw, zones_all = zones_all, covariates = dat$covariates,
                gt_pmfs = gt_pmfs, osrm_mat = osrm_mat, mobility_matrices = mobility_matrices,
                pop = dat$pop, ll = dat$ll)
  layer$cascade_key <- .layer_key
  # The optional build cache. Written HERE, before cascade_run_id is set below, and never
  # re-written when it is loaded, so it carries no run id; nothing downstream reads it (they read
  # the canonical copy written below).
  if (nzchar(CACHE)) saveRDS(layer, CACHE)
  tick("layer reconstructed")
}
# THE LAYER THIS RUN USED IS ALWAYS WRITTEN to the canonical path the later stages read
# (41_kindiv_sweep.R, make_manuscript_figure2_cascade.R), stamped with this run's id and week
# grid. It was written only when CASCADE_LAYER_CACHE was set — never on a default run_all.R pass —
# so those scripts rebuilt per-origin baselines from whatever layer a hand-run had left there,
# on a different week grid from the backtest they were drawn beside.
layer$cascade_run_id <- CASCADE_RUN_ID
layer$week_grid <- sort(unique(as.Date(layer$zone_week_nc$week_start)))
.canon_cache <- file.path(OUT_CASCADE, "layer_cache_v2.rds")
tryCatch(saveRDS(layer, .canon_cache),
         error = function(e) stop("[cascade] could not write the layer cache ", .canon_cache, ": ",
                                  conditionMessage(e), call. = FALSE))
zones_all <- layer$zones_all; gt_pmfs <- layer$gt_pmfs
province_map <- load_province_map()
zone_province <- setNames(province_map$province, province_map$nom)
vuln <- compute_vulnerability_index(layer$covariates, zones_all, layer$osrm_mat)

# ---- fit covariate-bearing hazards (featured + sensitivity kernels) --------
# ON THE VALIDATED MODEL'S DESIGN: rolling as-of predictors with the shared as-of truncation,
# exactly as run_all.R builds the cross-validated and deployed fits (trunc_delay = .trunc_asof).
# This used to fit on the shared count matrix, so the hazard's intercept — and therefore the
# meaning of every delta carried into it, 16b's short-horizon factor included — was on a
# different scale from the model it claims to reuse (cascade_invasion_design(), 33b).
.trunc_asof <- cascade_asof_truncation(layer$ll)
kernels <- unique(c(CASCADE_KERNEL, if (RUN_SENS) CASCADE_SENS_KERNELS))
kernels <- kernels[kernels %in% names(layer$mobility_matrices)]
iter_fit <- if (SMOKE) 800L else 2000L
designs <- list(); fits <- list()
# THE FEATURED MODEL'S FORM (CASCADE_FORM, 30_projection_config.R): fixed beta, or the
# time-varying form the two-stage selection featured, refitted with the cascade's covariates on
# every kernel through the one entry point the calibration origins use (cascade_fit_hazard()).
for (kern in kernels) {
  designs[[kern]] <- cascade_invasion_design(layer$zone_week_nc, layer, kernel = kern,
                                             gt = CASCADE_GT, linelist = layer$ll,
                                             trunc_delay = .trunc_asof)
  fits[[kern]] <- cascade_fit_hazard(designs[[kern]], layer$zone_week_nc, layer, kernel = kern,
                                     gt = CASCADE_GT, linelist = layer$ll,
                                     trunc_delay = .trunc_asof, iter = iter_fit)
}
if (is.null(fits[[CASCADE_KERNEL]]))
  stop(sprintf("[run_cascade] the %s hazard on the featured kernel %s did not fit.", CASCADE_FORM,
               CASCADE_KERNEL), call. = FALSE)
fits <- Filter(Negate(is.null), fits)
ROLLING_ON <- isTRUE(designs[[CASCADE_KERNEL]]$rolling_predictors)
tick(sprintf("fitted %d kernel hazard(s): %s (form %s; rolling as-of predictors: %s)", length(fits),
             paste(names(fits), collapse=","), CASCADE_FORM, if (ROLLING_ON) "ON" else "OFF"))

# The import term in the CONJUGATE estimator's denominator needs a hazard scale, so a delta is
# resolved first. It comes from 16b's prequential recalibration of the short-horizon model — not
# from the cascade — so there is no circularity: the cascade's own calibration consumes this
# term, it does not produce the delta that built it.
#
# Its reach is now much smaller than it was. The import term feeds estimate_zone_reff(), which
# is a DIAGNOSTIC: the reproduction number the projection transmits on comes from EpiNow2 and
# never passes through delta_imp. The note further down records any disagreement with the
# K-week delta fitted on the calibration windows rather than warning, because the two are no
# longer required to agree.
delta_imp <- tryCatch(cascade_fit_delta(cov_model = cascade_hazard_label(CASCADE_KERNEL)),
                      error = function(e) NA_real_)
if (!is.finite(delta_imp) || delta_imp <= 0) {
  warning("[run_cascade] no recalibration factor available for the import term; using the RAW ",
          "hazard scale, which over-states importation and so under-states local R.",
          call. = FALSE)
  delta_imp <- 1
}
reff  <- cascade_reff(layer$zone_week_nc, zones_all, layer,
                      fits[[CASCADE_KERNEL]], designs[[CASCADE_KERNEL]],
                      delta = delta_imp, zone_province = zone_province,
                      gt = CASCADE_GT, kernel = CASCADE_KERNEL,
                      # The production anchor is fitted as of the analysis date: this is the
                      # one call that is SUPPOSED to see the whole record.
                      linelist = layer$ll, issue_date = ANALYSIS_DATE)
prep  <- cascade_prepare(layer, fits[[CASCADE_KERNEL]], designs[[CASCADE_KERNEL]], reff,
                         kernel = CASCADE_KERNEL, gt = CASCADE_GT)

# ---- CALIBRATION (33b) -----------------------------------------------------
# ONE fitted parameter. delta is fitted against the realised new-invasion counts over K-week
# windows at the calibration origins, so the cascade's expected count matches the count that
# actually occurred. Each origin's hazard and R are trained on counts as known at that origin (or,
# where that reconstruction fails, the final-count slice — recorded per origin as
# training_source), but delta itself is FITTED to those windows' outcomes: the level calibration
# is IN-SAMPLE, and re-scoring the same windows cannot test it. psi (frontier saturation)
# is no longer fitted: it was not identifiable on this outbreak and is fixed at 0
# (30_projection_config.R). The h=1 refit and 16b's short-horizon factor are both retained as
# comparisons, so the gap between a one-week and a K-week level is reported rather than hidden.
#
# A FAILED CALIBRATION IS FATAL. It used to be caught by a tryCatch that only warned, so when
# cascade_calibration_origins() went missing (a stray `if (!exists("%||%"))` in 33b swallowed
# its definition) the run carried on with 16b's ONE-week factor, zero calibration uncertainty,
# and a report that read "fitted OUT OF SAMPLE against realised NA-week invasion counts". The
# same now holds for a calibration that returns but fits nothing (no settled K-week origin).
# CASCADE_CALIBRATE=0 runs deliberately uncalibrated — e.g. when no held-out origin is settled
# yet — and every product then says so (`calibrated = FALSE`); CASCADE_SMOKE=1 implies it.
.reset_outputs(file.path(.cd, c("cascade_calibration.json", "cascade_delta_fit.json",
                                "cascade_delta_fit_curve.csv", "cascade_psi_fit.json",
                                "cascade_psi_fit_curve.csv", "cascade_calibration_origin_reff.csv",
                                "cascade_calibration_per_origin.csv",
                                "cascade_calibration_pooled.csv",
                                "cascade_calibration_reliability.csv",
                                "cascade_psi_pace_diagnostic_curve.csv")), "calibration")
CALIB_K  <- as.integer(Sys.getenv("CASCADE_CALIB_K", unset = "6"))
CALIB_DISABLED_BY <- if (SMOKE) "CASCADE_SMOKE=1" else
                     if (identical(Sys.getenv("CASCADE_CALIBRATE"), "0")) "CASCADE_CALIBRATE=0" else
                     NA_character_
DO_CALIB <- is.na(CALIB_DISABLED_BY)
.calib_fail <- function(why) {
  .write_stamp(paste("failed: calibration -", why))
  stop(sprintf(paste0("[run_cascade] calibration FAILED (%s). The 13-week level would be ",
                      "uncalibrated at the horizon it is read at, so nothing is published. Fix ",
                      "the cause, or set CASCADE_CALIBRATE=0 to run DELIBERATELY uncalibrated ",
                      "(every product is then labelled calibrated = FALSE)."), why), call. = FALSE)
}
calib <- NULL
if (DO_CALIB) {
  calib <- tryCatch(cascade_calibrate(layer, zone_province = zone_province, K = CALIB_K,
             kernel = CASCADE_KERNEL, gt = CASCADE_GT, seed = CASCADE_SEED,
             # 600, not 400: n_mc_report DEFAULTS TO n_mc_search, so the delta search and the
             # residual check use the same Monte-Carlo size and their modelled
             # counts are comparable. They are NOT identical — the residual check runs at an
             # offset seed on purpose (33b, "INDEPENDENT SEED"), so the two differ by Monte-
             # Carlo error. Matching the SIZE is what keeps that error small; at 400 vs 600
             # the report printed 72.96 and 71.7 for one quantity.
             n_mc_search = as.integer(Sys.getenv("CASCADE_CALIB_NMC", unset = "600"))),
           error = function(e) .calib_fail(conditionMessage(e)))
  if (is.null(calib$delta_oos) || !is.finite(calib$delta_oos$delta %||% NA_real_))
    .calib_fail(sprintf("no delta could be fitted on settled %d-week held-out outcomes: %s",
                        CALIB_K, calib$note %||% "no note"))
} else {
  # DELIBERATELY uncalibrated. The short-horizon factor is carried WITH its attributes, so its
  # zone-clustered bootstrap SE still reaches delta_sd_log and the projection's intervals: this
  # path used to keep only `estimator`, which silently set delta_sd_log to 0 and delta_prior to NA.
  .d0 <- cascade_fit_delta(cov_model = cascade_hazard_label(CASCADE_KERNEL))
  calib <- list(delta = .d0, psi = CASCADE_PSI, psi_fit = NULL, delta_oos = NULL,
                delta_prior = as.numeric(.d0), report = NULL, origins = list(),
                delta_info = cascade_short_horizon_delta_info(.d0), rolling_predictors = NA,
                note = sprintf(paste0("calibration DISABLED deliberately (%s): delta is 16b's ",
                                      "SHORT-HORIZON (1-week) factor, so the 13-week level is ",
                                      "uncalibrated at the horizon it is read at; psi at its ",
                                      "configured default"), CALIB_DISABLED_BY))
  warning("[run_cascade] ", calib$note, call. = FALSE, immediate. = TRUE)
}
CALIBRATED <- !is.null(calib$delta_oos) && is.finite(calib$delta_oos$delta %||% NA_real_)
CALIB_STATUS <- if (!CALIBRATED) {
  sprintf("not calibrated (disabled: %s)", CALIB_DISABLED_BY)
} else if (isTRUE(calib$delta_oos$boundary_hit)) {
  sprintf("calibrated, CENSORED at the %s bound", calib$delta_oos$boundary_side)
} else if (!isTRUE(calib$delta_oos$converged)) "calibrated, NOT converged" else "calibrated"
delta <- calib$delta
psi   <- calib$psi
# delta_imp (16b's short-horizon factor) built the import term in the CONJUGATE estimator's
# denominator; `delta` is now fitted to the realised K-week counts at the calibration origins
# (in-sample for that level) and drives the simulated hazard. They are no longer required to agree, and the old >25% warning would fire on
# a legitimate difference.
#
# WHY THE COUPLING NO LONGER MATTERS, which is the substantive change: the import term feeds
# estimate_zone_reff(), and that estimator is now a DIAGNOSTIC. The reproduction number the
# projection transmits on comes from EpiNow2 and never passes through delta_imp at all, so a
# mismatch can no longer put R and the deployed hazard on different scales. It is reported as a
# note, because a large gap still says something about the model — that the level implied at one
# week and the level implied over K weeks disagree — and that belongs in the record.
if (is.finite(delta) && delta > 0 && is.finite(delta_imp) && delta_imp > 0 &&
    abs(log(delta / delta_imp)) > log(1.25))
  message(sprintf(paste0("[run_cascade] the %d-week delta fitted on the calibration windows ",
                         "(%.3f) differs from the short-horizon factor the import term was built ",
                         "on (%.3f) by %.0f%%. ",
                         "This no longer affects R (EpiNow2 supplies it); it records that the ",
                         "one-week and %d-week levels disagree."),
                  CALIB_K, delta, delta_imp, 100 * abs(delta / delta_imp - 1), CALIB_K))
# psi_fit must be a local too: the k_indiv sweep's provenance sidecar below refers to a
# bare `psi_fit`, which was a top-level object before calibration moved into
# cascade_calibrate(). Without this the sweep block dies with "object 'psi_fit' not
# found" AFTER the cascade has run, losing every product written past that point.
psi_fit <- calib$psi_fit
# Calibration uncertainty enters the Monte Carlo: delta is drawn per parameter group on
# the log scale using the zone-clustered bootstrap SE that came with the estimate.
delta_sd_log <- { v <- calib$delta_info$se_log %||% NA_real_
                  if (is.finite(v) && v > 0) v else 0 }
tick(sprintf("reff R_nat=%.2f | %s | delta_sd_log=%.3f", reff$R_nat, calib$note, delta_sd_log))

# The reproduction numbers the projection actually ran on. These were previously not
# written anywhere: R_eff drives every zone's onward transmission and therefore the whole
# 13-week extrapolation, and a manuscript cannot report or defend numbers that exist only
# inside a function call. The single national R is recorded with its posterior; the per-zone
# file keeps the case/importation bookkeeping and the conjugate estimator's DIAGNOSTICS (which
# do not drive the projection).
.cd <- file.path(OUT_CASCADE, "diagnostics")
dir.create(.cd, recursive = TRUE, showWarnings = FALSE)
# TWO files, because there are now two different kinds of quantity and writing them in one
# table would misrepresent the model. The reproduction number is a SINGLE NATIONAL estimate:
# repeating it down 519 rows of a per-zone file, as this did, reads as 519 per-zone estimates
# that happen to agree, which is the opposite of what the model says. The per-zone file keeps
# only what genuinely varies by zone — the case and importation bookkeeping, and the conjugate
# estimator's per-zone diagnostics, all explicitly labelled as diagnostics.
readr::write_csv(tibble::tibble(
    health_zone = zones_all, province = reff$prov,
    cases_in_window = as.numeric(reff$ncase),
    local_denominator = as.numeric(reff$loc), imported_cases = as.numeric(reff$imp),
    # Conjugate-estimator diagnostics. These do NOT drive the projection — the cascade
    # transmits on the single national EpiNow2 R in cascade_reff_national.csv — and are kept
    # so the per-zone picture behind that national number stays inspectable.
    diag_R_zone_conjugate = as.numeric(reff$R_zone_conjugate %||% rep(NA_real_, length(zones_all))),
    diag_a_post = as.numeric(reff$a_post), diag_b_post = as.numeric(reff$b_post),
    diag_prior_dominated = reff$prior_dominated, diag_implausible = reff$implausible,
    was_active_before = prep$affected0),
  file.path(.cd, "cascade_reff_by_zone.csv"))

# The number the 13-week projection actually runs on, with its posterior and its provenance.
readr::write_csv(tibble::tibble(
    quantity = c("R_national_epinow2_window_mean", "R_national_conjugate_diagnostic"),
    estimate = c(reff$R_nat, reff$R_nat_conjugate),
    lo_90 = c(unname(reff$rt_quantiles[1]), NA_real_),
    q25   = c(unname(reff$rt_quantiles[2]), NA_real_),
    median= c(unname(reff$rt_quantiles[3]), NA_real_),
    q75   = c(unname(reff$rt_quantiles[4]), NA_real_),
    hi_90 = c(unname(reff$rt_quantiles[5]), NA_real_),
    sd_log = c(as.numeric(reff$sd_zone[1]), NA_real_),
    n_draws = c(reff$rt_n_draws, NA_integer_),
    n_clamped = c(reff$rt_n_clamped, NA_integer_),
    source = c(reff$rt_source, "conjugate renewal (diagnostic only)"),
    window_start = as.Date(c(reff$rt_window_start, NA)),
    window_end   = as.Date(c(reff$rt_window_end, NA)),
    window_weeks = c(as.integer(reff$rt_window_weeks), NA_integer_),
    issue_date   = as.Date(c(reff$rt_issue_date, reff$rt_issue_date)),
    gt_profile   = CASCADE_GT),
  file.path(.cd, "cascade_reff_national.csv"))
# The draws themselves: the projection resamples these directly, so a reader can reproduce the
# R input to the Monte Carlo rather than re-deriving it from a mean and an SD.
readr::write_csv(tibble::tibble(draw = seq_along(reff$R_draws), R = as.numeric(reff$R_draws)),
                 file.path(.cd, "cascade_reff_national_draws.csv"))

# GENERATION-TIME MARGINALISATION — OFF by default since 2026-09-22 (CASCADE_GT_MARGINALISE),
# so this block writes nothing on a default run and the projection uses the single CASCADE_GT
# anchor. It is retained for the sensitivity arm. When the
# projection marginalises over GT_PRIOR it runs on FIFTEEN reproduction numbers, one per grid
# point, each paired with its own weekly kernel — and until this block existed none of that
# reached an artifact, so a manuscript could cite only the single-GT anchor while the run used
# something else. Each row is one grid point: its generation time, its prior weight, and the R
# posterior fitted AT that generation time.
# Removed first: on a default (non-marginalised) run a previous sensitivity run's files would
# otherwise sit beside this run's single-GT R as though they described it.
.reset_outputs(file.path(.cd, c("cascade_reff_by_gt.csv", "cascade_reff_gt_marginal.json")),
               "GT marginalisation")
if (!is.null(reff$R_draws_by_gt) && length(reff$R_draws_by_gt)) {
  .gq <- function(v, p) unname(stats::quantile(v, p, names = FALSE))
  readr::write_csv(tibble::tibble(
      gt_key   = reff$gt_keys,
      gt_mean  = reff$gt_grid$gt_mean,
      gt_sd    = reff$gt_grid$gt_sd,
      prior_weight = as.numeric(reff$gt_weights),
      R_mean   = vapply(reff$R_draws_by_gt, mean, numeric(1)),
      R_lo_90  = vapply(reff$R_draws_by_gt, .gq, numeric(1), p = 0.05),
      R_median = vapply(reff$R_draws_by_gt, .gq, numeric(1), p = 0.50),
      R_hi_90  = vapply(reff$R_draws_by_gt, .gq, numeric(1), p = 0.95),
      n_draws  = vapply(reff$R_draws_by_gt, length, integer(1))),
    file.path(.cd, "cascade_reff_by_gt.csv"))
  jsonlite::write_json(list(
      gt_marginalised = TRUE,
      n_grid_points = length(reff$gt_keys),
      R_prior_weighted = reff$R_nat_gt_marginal,
      R_single_gt_anchor = reff$R_nat,
      pct_apart = 100 * abs(reff$R_nat_gt_marginal / max(reff$R_nat, 1e-12) - 1),
      R_min_across_grid = min(vapply(reff$R_draws_by_gt, mean, numeric(1))),
      R_max_across_grid = max(vapply(reff$R_draws_by_gt, mean, numeric(1))),
      gt_prior = get0("GT_PRIOR", ifnotfound = NULL)),
    file.path(.cd, "cascade_reff_gt_marginal.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = 6, null = "null")
  tick(sprintf(paste0("GT marginalisation written: %d grid points, prior-weighted R = %.3f ",
                      "(single-GT anchor %.3f, %.1f%% apart), R spans %.3f-%.3f"),
               length(reff$gt_keys), reff$R_nat_gt_marginal, reff$R_nat,
               100 * abs(reff$R_nat_gt_marginal / max(reff$R_nat, 1e-12) - 1),
               min(vapply(reff$R_draws_by_gt, mean, numeric(1))),
               max(vapply(reff$R_draws_by_gt, mean, numeric(1)))))
}
tick(sprintf(paste0("R written: national EpiNow2 R = %.3f [%.3f, %.3f] over %s..%s (%s); ",
                    "conjugate diagnostic %.3f (%.0f%% %s); %d zones with cases, importation ",
                    "%.1f%% of the denominator"),
             reff$R_nat, reff$rt_quantiles[1], reff$rt_quantiles[5],
             format(as.Date(reff$rt_window_start)), format(as.Date(reff$rt_window_end)),
             reff$rt_source, reff$R_nat_conjugate,
             100 * abs(reff$R_nat / max(reff$R_nat_conjugate, 1e-12) - 1),
             if (reff$R_nat >= reff$R_nat_conjugate) "higher" else "lower",
             reff$n_with_cases, 100 * reff$import_share))

# Persist the whole calibration, including the delta search path (so a censored or
# non-converged fit is auditable rather than reduced to one scalar) and the held-out
# residual check (whose COUNT is in-sample for delta; see cascade_calibration_report()).
.cd <- file.path(OUT_CASCADE, "diagnostics")
if (!is.null(calib$delta_oos)) {
  readr::write_csv(.stamp_df(calib$delta_oos$curve), file.path(.cd, "cascade_delta_fit_curve.csv"))
  jsonlite::write_json(c(list(cascade_run_id = CASCADE_RUN_ID,
                              analysis_date = format(as.Date(ANALYSIS_DATE))),
                         calib$delta_oos[setdiff(names(calib$delta_oos), "curve")]),
                       file.path(.cd, "cascade_delta_fit.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = 6)
  if (!isTRUE(calib$delta_oos$converged))
    warning("[run_cascade] the K-week delta search did not converge — the 13-week layer is NOT ",
            "fit for publication until this is resolved (see cascade_delta_fit.json).",
            call. = FALSE)
}
# psi is no longer fitted (30_projection_config.R). This block stays so that a run which
# deliberately re-enables the retired search still writes its path rather than dropping it.
if (!is.null(calib$psi_fit)) {
  readr::write_csv(calib$psi_fit$curve, file.path(.cd, "cascade_psi_fit_curve.csv"))
  jsonlite::write_json(calib$psi_fit[setdiff(names(calib$psi_fit), "curve")],
                       file.path(.cd, "cascade_psi_fit.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = 6)
}
if (!is.null(calib$origin_diagnostics) && NROW(calib$origin_diagnostics))
  readr::write_csv(.stamp_df(calib$origin_diagnostics),
                   file.path(.cd, "cascade_calibration_origin_reff.csv"))
if (!is.null(calib$report)) {
  readr::write_csv(.stamp_df(calib$report$per_origin), file.path(.cd, "cascade_calibration_per_origin.csv"))
  readr::write_csv(.stamp_df(calib$report$pooled),     file.path(.cd, "cascade_calibration_pooled.csv"))
  if (!is.null(calib$report$reliability))
    readr::write_csv(.stamp_df(calib$report$reliability), file.path(.cd, "cascade_calibration_reliability.csv"))
  tick(paste("calibration residual check (count in-sample):", calib$report$note))
}
jsonlite::write_json(list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    cascade_run_id = CASCADE_RUN_ID, analysis_date = format(as.Date(ANALYSIS_DATE)),
    kernel = CASCADE_KERNEL, form = CASCADE_FORM, calib_K = CALIB_K,
    # THE FLAG EVERY CONSUMER BRANCHES ON. FALSE only when calibration was disabled
    # deliberately (CASCADE_CALIBRATE=0 / smoke); a FAILED calibration stops the run instead.
    calibrated = CALIBRATED, calibration_status = CALIB_STATUS,
    calibration_disabled_by = CALIB_DISABLED_BY,
    # Whether the hazards (production and every calibration origin) were refit on the validated
    # model's rolling as-of design (cascade_invasion_design()).
    rolling_predictors = ROLLING_ON,
    rolling_predictors_origins = calib$rolling_predictors %||% NA,
    settled_label_age_days = get0("LFO_RELIABLE_EVAL_AGE_DAYS", ifnotfound = NA_integer_),
    delta = as.numeric(delta), delta_info = calib$delta_info, delta_sd_log = delta_sd_log,
    delta_short_horizon = calib$delta_prior,
    delta_estimator = if (is.null(calib$delta_oos)) "short-horizon factor (no K-week fit)"
                      else sprintf("%d-week count match at held-out origins", CALIB_K),
    delta_oos_converged = if (is.null(calib$delta_oos)) NA else isTRUE(calib$delta_oos$converged),
    delta_oos_boundary_hit = if (is.null(calib$delta_oos)) NA
                             else isTRUE(calib$delta_oos$boundary_hit),
    delta_oos_modelled_cum = if (is.null(calib$delta_oos)) NA_real_
                             else calib$delta_oos$modelled_cum_at_delta,
    delta_oos_observed_cum = if (is.null(calib$delta_oos)) NA_real_
                             else calib$delta_oos$target_observed_cum,
    # psi is FIXED, not fitted (30_projection_config.R). Recorded explicitly so a reader of this
    # file can tell "held at 0 by design" from "fitted and happened to return 0".
    psi = psi, psi_fitted = !is.null(calib$psi_fit),
    psi_converged = if (is.null(calib$psi_fit)) NA else isTRUE(calib$psi_fit$converged),
    psi_boundary_hit = if (is.null(calib$psi_fit)) NA else isTRUE(calib$psi_fit$boundary_hit),
    # The pooled residual-check row (cascade_calibration_report()$pooled). RENAMED from
    # `verification`: the pooled count is the quantity delta was fitted to, so its residual is
    # in-sample by construction. Recorded under a name that says what it is, with explicit flags,
    # so nothing downstream can read it as validation. Its AUC-PR / skill columns are the part
    # delta was not fitted to.
    count_residual_in_sample = if (is.null(calib$report)) NULL else as.list(calib$report$pooled),
    count_residual_ok        = if (is.null(calib$report)) NA else isTRUE(calib$report$count_residual_ok),
    count_level_in_sample    = if (is.null(calib$report)) NA else isTRUE(calib$report$count_level_in_sample),
    note = calib$note),
  file.path(.cd, "cascade_calibration.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6,
  null = "null", na = "null")

# Naive-pace comparison, retired as an estimator (see 33_cascade_eval.R) and off by default.
if (identical(Sys.getenv("CASCADE_RUN_PACE_DIAG"), "1")) {
  pace <- tryCatch(cascade_fit_psi(prep, as.numeric(delta), layer$zone_week_nc, zones_all,
                                   pop_vec = layer$pop),
                   error = function(e) NULL)
  if (!is.null(pace)) {
    tick(paste("naive-pace diagnostic:", pace$note))
    readr::write_csv(pace$curve, file.path(.cd, "cascade_psi_pace_diagnostic_curve.csv"))
  }
}

# ---- consistency gate ------------------------------------------------------
gate <- cascade_consistency_gate(prep, fits[[CASCADE_KERNEL]], designs[[CASCADE_KERNEL]], delta,
          layer$zone_week_nc, layer$mobility_matrices, gt_pmfs, layer$covariates,
          layer$osrm_mat, zones_all, n_mc = if (SMOKE) 1000L else 3000L)
tick(paste("consistency gate:", gate$note))

# ---- run all scenarios -----------------------------------------------------
sims <- list(); reaches <- list(); rs_by_scenario <- list()
for (sk in names(CASCADE_SCENARIOS)) {
  # delta_sd_log applies to the PUBLISHED scenario products only. The comparative
  # analyses below (gateway knockout, conditional seeding, the k_indiv sweep) are
  # DIFFERENCES against a common baseline, so adding calibration noise there would
  # widen the contrast without moving it; they stay on the point estimate.
  sims[[sk]] <- simulate_cascade(prep, CASCADE_SCENARIOS[[sk]], n_mc = CASCADE_N_MC,
                                 delta = delta, psi = psi, delta_sd_log = delta_sd_log,
                                 pop_vec = layer$pop)
  reaches[[sk]] <- cascade_reach_table(sims[[sk]])
  rs_by_scenario[[sk]] <- cascade_enrich_reach(reaches[[sk]], province_map, vuln,
                            method_label = paste0("Cascade-", CASCADE_KERNEL, "-", sk))
}
scP_key <- CASCADE_SCENARIO_PRIMARY
rs_primary <- rs_by_scenario[[scP_key]]; base_sim <- sims[[scP_key]]
tick(sprintf("simulated %d scenarios x M=%d", length(sims), CASCADE_N_MC))

# ---- aggregates for the report ---------------------------------------------
# Stage-failure handler: message() ALONE is not enough for a validation stage. A message does
# not survive into warnings(), leaves no artifact, and 36_report.R simply omits the section —
# so a crashed backtest or sensitivity sweep was indistinguishable from one that passed.
# Raise a real warning as well, immediately, then return NULL so the caller degrades as before.
.stage_failed <- function(what) function(e) {
  warning(sprintf("[cascade] %s FAILED: %s -- its report section will be marked NOT AVAILABLE.",
                  what, conditionMessage(e)), call. = FALSE, immediate. = TRUE)
  NULL
}

prov_of <- function(z) zone_province[match(z, names(zone_province))]
affected_provinces0 <- unique(prov_of(zones_all[base_sim$affected0]))
agg <- lapply(names(sims), function(sk) {
  sim <- sims[[sk]]; cum <- apply(sim$new_by_week, 2, cumsum)
  # Horizons from CASCADE_REPORT_HORIZONS, not the literals 4/8/13: hard-coded row indices
  # abort with a subscript error (or silently report the wrong weeks) if the reporting
  # horizons or CASCADE_HORIZON_WEEKS change — and only AFTER every scenario has been simulated.
  .h <- sort(unique(as.integer(CASCADE_REPORT_HORIZONS)))
  .h <- .h[.h >= 1L & .h <= nrow(cum)]
  .hmax <- max(.h)
  # new provinces reached by the last reporting horizon (provinces with no affected zone at t0)
  # DROP unmapped zones rather than counting them as a province. prov_of() returns NA for any
  # zone absent from province_map; NA survives setdiff() unless it is already in
  # affected_provinces0, so a single unmapped zone being reached added a phantom "newly reached
  # province" in every iteration where it was reached. cascade_flux_network() guards the same
  # case (34_conditional_queries.R:133); this call did not.
  .provs <- function(z) { pv <- prov_of(z); unique(pv[!is.na(pv)]) }
  newprov <- vapply(seq_len(sim$n_mc), function(m) {
    zr <- zones_all[!is.na(sim$tau[, m]) & sim$tau[, m] <= .hmax]
    length(setdiff(.provs(zr), affected_provinces0)) }, numeric(1))
  # ONE consistent estimator throughout: the MEDIAN of the cumulative-new-zone
  # distribution at each horizon (+ 90% predictive band), so 4/8/13w are comparable
  # (means and medians differ for the right-skewed count; do not mix them).
  q <- function(h, p) stats::quantile(cum[h, ], p, names = FALSE)
  list(label = CASCADE_SCENARIOS[[sk]]$label,
       new4 = stats::median(cum[.h[1], ]),
       new8 = stats::median(cum[.h[min(2L, length(.h))], ]),
       new13_med = stats::median(cum[.hmax, ]),
       new13_lo = q(.hmax, .05), new13_hi = q(.hmax, .95),
       # MEDIAN, not mean: this sits in a row of medians, and the block's own comment says
       # "ONE consistent estimator throughout ... do not mix them". mean(newprov) mixed one in.
       new_prov13 = stats::median(newprov))
})
names(agg) <- names(sims)

# ---- outputs: CSVs, maps, figures ------------------------------------------
TAB <- file.path(OUT_CASCADE, "tables")
cascade_write_reach_csv(rs_primary, file.path(TAB, "cascade_reach_scores_all_zones.csv"))
for (sk in names(rs_by_scenario))
  cascade_write_reach_csv(rs_by_scenario[[sk]], file.path(TAB, sprintf("cascade_reach_%s.csv", sk)))
for (H in CASCADE_REPORT_HORIZONS) {
  cascade_map_reach(rs_primary, H, scenario_label = CASCADE_SCENARIOS[[scP_key]]$label)
  cascade_map_uncertainty(rs_primary, H)
}
cascade_map_timetoinvasion(rs_primary)
# PUBLISH THE BAND, not just the picture. cascade_fig_fanchart() returns the per-week median
# and 90% interval its ribbon draws; that result was assigned here and never used, so the
# credible band on a published figure was unreadable as data.
.reset_outputs(file.path(OUT_CASCADE, "tables", "cascade_newzones_fanchart.csv"), "fan-chart band")
fan <- cascade_fig_fanchart(sims)
tryCatch({
  if (!is.null(fan) && nrow(fan)) {
    .fanp <- file.path(OUT_CASCADE, "tables", "cascade_newzones_fanchart.csv")
    dir.create(dirname(.fanp), recursive = TRUE, showWarnings = FALSE)
    # The primary scenario flagged in the table, so a reader (37_cascade_figure4.R) can pick the
    # headline total without re-deriving the scenario config.
    fan$is_primary <- fan$scenario == CASCADE_SCENARIOS[[scP_key]]$label
    readr::write_csv(dplyr::mutate(fan, dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 12))),
                     .fanp)
    message(sprintf("[cascade] wrote the fan-chart band (%d row(s)) -> %s",
                    nrow(fan), basename(.fanp)))
  }
}, error = function(e)
  warning("[cascade] fan-chart band not written: ", conditionMessage(e), call. = FALSE))
cascade_map_scenarios(rs_by_scenario, horizon = max(CASCADE_REPORT_HORIZONS))
cascade_fig_establishment(rs_primary)
cascade_render_reuse_maps(rs_primary)
tick("maps + figures written")

# ---- onward propagation: gateway, corridors, tree, conditional -------------
Ync <- .count_wide(layer$zone_week_nc, zones_all, "confirmed_nc")
recent_inc <- rowSums(Ync[, tail(seq_len(ncol(Ync)), 4), drop = FALSE])
cand_hubs <- names(sort(recent_inc[base_sim$affected0], decreasing = TRUE))
cand_hubs <- head(cand_hubs[recent_inc[cand_hubs] > 0], CASCADE_KNOCKOUT_TOPN)
gateway <- cascade_knockout(prep, CASCADE_SCENARIOS[[scP_key]], base_sim, delta, psi = psi,
                            n_mc = min(CASCADE_N_MC, if (SMOKE) 40L else 600L),
                            candidates = cand_hubs, pop_vec = layer$pop)
readr::write_csv(gateway, file.path(TAB, "cascade_gateway_knockout.csv"))
cascade_fig_gateway(gateway)
flux_net <- cascade_flux_network(base_sim, zone_province)
readr::write_csv(flux_net$all_edges, file.path(TAB, "cascade_flux_edges.csv"))
cascade_fig_corridors(flux_net)
readr::write_csv(cascade_transmission_tree(base_sim), file.path(TAB, "cascade_transmission_tree.csv"))
tick("gateway + corridors + tree written")

# ---- designed conditional experiment: force-seed a hub NOT YET INVADED --------------
# THE HUBS MUST BE ZONES THAT CAN STILL BE INVADED. They were taken from the gateway list —
# head(gateway$health_zone, 2), the currently INFECTED sources — and simulate_cascade() forces
# only zones that are not already infected, so the "conditional" runs were the baseline under
# another name: identical across hubs, published as each hub's downstream reach.
#
# RULE (stated in the report): CASCADE_HUBS_FOR_COND when set; otherwise the two zones NOT
# invaded at the projection origin with the largest expected onward seedings in the baseline
# projection (cascade_gateway_screen()'s out-flux, the flux screen of PLAN §3.7) — the
# not-yet-invaded zones the projection expects to act as hubs. Each is force-seeded at week 1
# and its downstream effect is the PAIRED attributable increase over a same-seed, same-M
# baseline, ranked among zones whose paired 90% interval excludes zero — so the list is the
# hub's own catchment, not the zones that were already likely without it.
.reset_outputs(list.files(TAB, pattern = "^cascade_conditional_.*\\.csv$", full.names = TRUE),
               "conditional hubs")
conditional <- list()
.Hc  <- max(CASCADE_REPORT_HORIZONS)
.scr <- cascade_gateway_screen(base_sim)
.cond_explicit <- !is.null(get0("CASCADE_HUBS_FOR_COND", ifnotfound = NULL))
cond_hubs <- if (.cond_explicit) as.character(CASCADE_HUBS_FOR_COND) else
  head(.scr$health_zone[!.scr$was_active_before & .scr$out_flux > 0], 2L)
if (!length(cond_hubs)) {
  warning("[cascade] conditional experiment NOT RUN: no zone outside the current front seeded any ",
          "onward invasion in the baseline projection, so no hub qualifies under the rule.",
          call. = FALSE, immediate. = TRUE)
} else {
  .n_cond <- min(CASCADE_N_MC, if (SMOKE) 40L else 800L)
  scP_c   <- CASCADE_SCENARIOS[[scP_key]]
  # The paired baseline: same seed, same M, same (default) delta_sd_log as the conditional runs.
  .base_cond <- simulate_cascade(prep, scP_c, n_mc = .n_cond, delta = delta, psi = psi,
                                 seed = CASCADE_SEED, pop_vec = layer$pop)
  for (hub in cond_hubs) {
    cs <- cascade_conditional_seed(prep, scP_c, hub, delta, psi = psi, n_mc = .n_cond,
                                   seed = CASCADE_SEED, pop_vec = layer$pop)
    rs_c <- cascade_enrich_reach(cs$reach, province_map, vuln, method_label = paste0("cond-", hub))
    cascade_map_conditional(rs_c, hub)
    ctr <- cascade_paired_contrast(.base_cond, cs$sim, horizon = .Hc, exclude = hub)
    pz  <- ctr$per_zone[ctr$per_zone$eligible & ctr$per_zone$elevated %in% TRUE, , drop = FALSE]
    pz  <- pz[order(-pz$delta), , drop = FALSE]
    conditional[[hub]] <- list(
      top = head(pz$health_zone, 6), n_elevated = nrow(pz),
      expected_added = ctr$total, expected_added_lo = ctr$lo, expected_added_hi = ctr$hi,
      paired = isTRUE(ctr$paired), horizon = .Hc, n_mc = .n_cond,
      baseline_reach = rs_primary$p_case_invasion[rs_primary$health_zone == hub &
                                                    rs_primary$horizon == .Hc][1],
      out_flux = .scr$out_flux[match(hub, .scr$health_zone)],
      rule = if (.cond_explicit) "CASCADE_HUBS_FOR_COND" else "auto: top-2 not-yet-invaded by baseline onward seedings")
    # The conditional reach table, plus the paired attributable increase at the headline horizon
    # (joined on horizon, so the 4/8-week rows do not carry a 13-week contrast).
    .att <- tibble::tibble(health_zone = ctr$per_zone$health_zone, horizon = .Hc,
                           attributable_delta = ctr$per_zone$delta,
                           attributable_lo = ctr$per_zone$delta_lo,
                           attributable_hi = ctr$per_zone$delta_hi,
                           attributable_elevated = ctr$per_zone$elevated)
    readr::write_csv(dplyr::left_join(rs_c, .att, by = c("health_zone", "horizon")),
                     file.path(TAB, sprintf("cascade_conditional_%s.csv",
                                            gsub("[^A-Za-z0-9]+","_",hub))))
  }
}
tick("conditional hubs written")

# ---- urban-hub invasion scenarios (timing-aware) ---------------------------
# Previous tables removed whether or not the stage runs: 40_cascade_next_dominoes.R draws
# manuscript Figure 4 from them, and a disabled or failed stage must not leave last run's urban
# scenarios beside this run's reach table.
.reset_outputs(list.files(TAB, pattern = "^urban_(scenario_summary|impact_.*|cases_.*|hub_selection)\\.csv$",
                          full.names = TRUE), "urban scenarios")
urban <- NULL
if (RUN_URBAN) {
  urban <- tryCatch({
    hubs_u <- resolve_urban_hubs(URBAN_HUBS, zones_all)
    urban_write_hub_selection(TAB, hubs_u)   # seed-zone provenance (audit trail)
    sw_u   <- if (SMOKE) c(1L, 4L) else URBAN_SEED_WEEKS
    # dedicated baseline at URBAN_N_MC (matched to the urban conditional runs, so the
    # per-zone attributable difference for far zones is pure low-M noise, not a mix)
    # The baseline is kept as a SIM, not just a reach table: the contrast is paired
    # iteration by iteration, which needs both first-passage matrices. Same seed and same
    # M as every conditional run — that is the pairing.
    sim_base_u <- simulate_cascade(prep, CASCADE_SCENARIOS[[scP_key]],
                    n_mc = URBAN_N_MC, delta = delta, psi = psi, seed = CASCADE_SEED,
                    pop_vec = layer$pop)
    reach_base_u <- cascade_reach_table(sim_base_u, horizons = CASCADE_REPORT_HORIZONS)
    grid_u <- urban_scenario_grid(prep, CASCADE_SCENARIOS[[scP_key]], hubs_u, sw_u,
                delta, psi, n_mc = URBAN_N_MC, reach_base = reach_base_u,
                province_map = province_map, sim_base = sim_base_u, seed = CASCADE_SEED,
                pop_vec = layer$pop)
    urban_write_tables(grid_u, TAB)
    for (city in names(grid_u$detail)) {
      im <- grid_u$detail[[city]]$impact
      urban_map_attributable(im, city, grid_u$min_week)
      urban_fig_top_zones(im, city)
      urban_figure4(im, city, hubs_u[[city]])          # publication two-panel (map + dumbbell)
    }
    urban_fig_summary_bar(grid_u$summary, grid_u$min_week)
    urban_fig_seed_r(grid_u$summary, grid_u$min_week)
    urban_fig_cases(grid_u$summary, grid_u$min_week)
    urban_fig_timing(grid_u$summary)
    urban_map_facet(setNames(lapply(names(grid_u$detail),
                    function(c) grid_u$detail[[c]]$impact$per_zone), names(grid_u$detail)))
    grid_u
  }, error = .stage_failed("urban scenarios"))
  # Report the FAILURE as a failure. `length(if (!is.null(urban)) urban$detail else 0L)`
  # evaluates length(0L) = 1 on the error path, so a block that wrote nothing logged
  # "urban scenarios written (1 hubs x 1 seed-weeks)" — a success line with fabricated counts.
  tick(if (is.null(urban))
         "urban scenarios FAILED - nothing written (see the error above)"
       else sprintf("urban scenarios written (%d hubs x %d seed-weeks)",
                    length(urban$detail), length(urban$seed_weeks)))
}

# ---- validation: sensitivity, backtest, SBC --------------------------------
# NOTE: pass the DEPLOYED psi so the k_indiv and kernel axes are perturbed against a
# baseline at the SAME psi (that is what base_reach was built at). Omitting it ran those
# axes at a different psi from the baseline, conflating the k/kernel change with a psi
# change; the psi axis itself still sweeps psi internally. psi is now fixed at 0, so the
# baseline is the no-saturation model and the psi axis reports what enabling saturation
# would do — it is a sensitivity, not a perturbation around a fitted value.
.reset_outputs(file.path(OUT_CASCADE, "diagnostics",
                         c("cascade_sensitivity.csv", "cascade_kindiv_sweep.csv",
                           "cascade_kindiv_sweep_byzone.csv", "cascade_kindiv_sweep_meta.json")),
               "sensitivity + k_indiv sweep")
sens <- if (RUN_SENS) tryCatch(cascade_sensitivity(layer, designs, fits, reff, delta,
                        base_reach = reaches[[scP_key]], psi = psi),
                        error = .stage_failed("sensitivity")) else NULL
if (!is.null(sens)) readr::write_csv(.stamp_df(sens), file.path(OUT_CASCADE, "diagnostics", "cascade_sensitivity.csv"))

# ---- comprehensive individual-offspring overdispersion (k_indiv) sweep ------
# Dense k grid at the fitted (delta, psi); reuses the production `prep` (no refit).
ksweep <- if (RUN_SENS) tryCatch(
  cascade_kindiv_sweep(prep, CASCADE_SCENARIOS[[scP_key]], delta = delta, psi = psi,
                       pop_vec = layer$pop, province = zone_province,
                       k_grid = CASCADE_K_INDIV_SWEEP, n_mc = CASCADE_N_MC),
  error = .stage_failed("k_indiv sweep")) else NULL
if (!is.null(ksweep)) {
  readr::write_csv(ksweep$summary, file.path(OUT_CASCADE, "diagnostics", "cascade_kindiv_sweep.csv"))
  readr::write_csv(ksweep$byzone,  file.path(OUT_CASCADE, "diagnostics", "cascade_kindiv_sweep_byzone.csv"))
  # provenance sidecar (kept in step with the CSVs so kernel/delta/psi never go stale)
  ksweep_meta <- list(cascade_run_id = CASCADE_RUN_ID, analysis_date = format(as.Date(ANALYSIS_DATE)),
                      calibrated = CALIBRATED,
                      kernel = CASCADE_KERNEL, form = CASCADE_FORM, M = CASCADE_N_MC,
                      k_grid = CASCADE_K_INDIV_SWEEP, k_base = CASCADE_K_INDIV,
                      delta = as.numeric(delta), cal_in_large = attr(delta, "cal_in_large") %||% NA,
                      psi = psi,
                      psi_converged = if (is.null(psi_fit)) NA else isTRUE(psi_fit$converged),
                      psi_boundary_hit = if (is.null(psi_fit)) NA else isTRUE(psi_fit$boundary_hit),
                      psi_horizon_pace_ratio = if (is.null(psi_fit)) NA else psi_fit$horizon_pace_ratio,
                      gate_spearman = gate$spearman_signal, gate_n_signal = gate$n_signal,
                      gate_x2 = gate$x2, gate_x2_df = gate$x2_df, gate_x2_p = gate$x2_p,
                      # ~delta by construction: a plumbing check, not calibration (33).
                      gate_level_ratio = gate$level_ratio, gate_pass = isTRUE(gate$pass),
                      scenario = scP_key, seed = CASCADE_SEED)
  jsonlite::write_json(ksweep_meta, file.path(OUT_CASCADE, "diagnostics", "cascade_kindiv_sweep_meta.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = 6)
  tick(sprintf("k_indiv sweep: %d k-values (Spearman %.3f-%.3f vs k=%.2f)",
               nrow(ksweep$summary), min(ksweep$summary$spearman), max(ksweep$summary$spearman),
               CASCADE_K_INDIV))
}
# Reuses the calibration's held-out origins when they exist, so the per-origin hazards
# are fitted ONCE rather than twice, and the backtest scores exactly the origins the
# calibration was fitted on.
.reset_outputs(c(file.path(OUT_CASCADE, "diagnostics",
                           c("cascade_backtest.csv", "cascade_backtest_detail.csv",
                             "cascade_eval_figure_inputs.rds")),
                 file.path(get0("CASCADE_OUT_DIR", ifnotfound = OUT_DIR), "key_outputs",
                           "cascade_backtest_per_origin.csv")), "backtest")
backtest <- if (RUN_BT) tryCatch(cascade_backtest(layer, delta = delta, psi = psi, K = CALIB_K,
                        zone_province = zone_province, linelist = layer$ll,
                        origins = if (length(calib$origins)) calib$origins else NULL),
                        error = .stage_failed("backtest")) else NULL
# WHETHER THE BACKTEST'S LEVEL IS IN-SAMPLE. When calibration ran, the backtest scores exactly the
# origins and K delta was fitted on (calib$origins), so its count ratio / calibration-in-the-
# large is ~1 by construction; only the ranking metrics are out-of-sample. Recorded per row so no
# consumer has to infer it.
BT_LEVEL_IN_SAMPLE <- CALIBRATED && length(calib$origins) > 0L
if (!is.null(backtest)) {
  .bt_out <- .stamp_df(backtest); .bt_out$level_in_sample <- BT_LEVEL_IN_SAMPLE
  readr::write_csv(.bt_out, file.path(OUT_CASCADE, "diagnostics", "cascade_backtest.csv"))
  # per-zone detail (drives the Figure2_cascade capture-curve / mean-rank / top-K panels)
  bt_detail <- attr(backtest, "detail")
  if (!is.null(bt_detail)) {
    bt_detail <- .stamp_df(bt_detail); bt_detail$level_in_sample <- BT_LEVEL_IN_SAMPLE
    readr::write_csv(bt_detail, file.path(OUT_CASCADE, "diagnostics", "cascade_backtest_detail.csv"))
  }
}
sbc <- if (RUN_SBC) tryCatch(cascade_sbc(designs[[CASCADE_KERNEL]], n_sbc = 20L),
                        error = .stage_failed("sbc")) else NULL
tick("validation done")

# ---- report + key_outputs bundle -------------------------------------------
ctx <- list(reach_primary = rs_primary, n_mc = CASCADE_N_MC, delta = delta, psi = psi,
            calib = calib, delta_sd_log = delta_sd_log,
            # Provenance and the flags the report branches on (36_report.R never prints a
            # calibration claim that these do not support).
            cascade_run_id = CASCADE_RUN_ID, calibrated = CALIBRATED,
            # The run's hazard: a context re-rendered later describes THIS run, not the session's.
            form = CASCADE_FORM, kernel = CASCADE_KERNEL,
            calibration_status = CALIB_STATUS, calibration_disabled_by = CALIB_DISABLED_BY,
            rolling_predictors = ROLLING_ON, backtest_level_in_sample = BT_LEVEL_IN_SAMPLE,
            # Whether the CROSS-VALIDATED model rolls (ROLLING_PREDICTORS): 36_report.R needs it to
            # say whether a non-rolling cascade hazard is on that model's scale.
            rolling_predictors_requested = isTRUE(get0("ROLLING_PREDICTORS", ifnotfound = TRUE)),
            gate = gate, agg = agg, gateway = gateway, conditional = conditional,
            sens = sens, ksweep = ksweep, backtest = backtest, sbc = sbc, urban = urban,
            reff = reff, delta_imp = delta_imp)
report_path <- file.path(ST_DIR, "SPATIAL_INVASION_3MONTH_REPORT.md")
# Persist the report CONTEXT, not just the rendered report. Every number in the report is
# derived from `ctx`, so saving it means a wording or formatting fix can be re-rendered in
# seconds instead of repeating a full run to regenerate figures nobody changed.
tryCatch(saveRDS(ctx, file.path(OUT_CASCADE, "cascade_report_ctx.rds")),
         error = function(e)
           warning("[cascade] could not persist the report context: ", conditionMessage(e),
                   call. = FALSE))
cascade_write_report(report_path, ctx)

# comprehensive k_indiv sweep figure (house aesthetic) — reads the CSVs just written
if (!is.null(ksweep))
  tryCatch(source(file.path(ST_DIR, "42_kindiv_sweep_figure.R")),
           error = function(e) message("[cascade] k_indiv figure skipped: ", conditionMessage(e)))

# evaluation figure (Figure 2 style): leave-future-out backtest of the cascade —
# discrimination (A), prioritisation capture curve (B), ranking accuracy (C),
# per-round forecasts vs realised outcomes (D). Needs ctx$backtest + its detail attr.
tryCatch(cascade_fig_evaluation(ctx),
         error = function(e) message("[cascade] evaluation figure skipped: ", conditionMessage(e)))

# key_outputs bundle (mirror the short-horizon key_outputs style).
# CASCADE_OUT_DIR is OUT_DIR unless CASCADE_OUT_ROOT redirected it (30_projection_config.R):
# the bundle must follow the tables it bundles, or an imputation arm would copy its own
# sandboxed tables over the published key_outputs.
KO <- file.path(get0("CASCADE_OUT_DIR", ifnotfound = OUT_DIR), "key_outputs")
if (!dir.exists(KO)) dir.create(KO, recursive = TRUE)
figdir <- file.path(OUT_CASCADE, "figures")
ko_files <- c(file.path(TAB, "cascade_reach_scores_all_zones.csv"),
              file.path(TAB, "cascade_gateway_knockout.csv"),
              file.path(TAB, "urban_scenario_summary.csv"),
              # The national new-zone TOTAL with its interval: mean, median and the 5th-95th
              # percentiles of the SIMULATED total per week and scenario. Bundled beside the
              # per-zone reach scores because a total's interval must come from here — summing
              # the per-zone p_case_lo/p_case_hi is not an interval for the total.
              file.path(TAB, "cascade_newzones_fanchart.csv"),
              list.files(figdir, pattern = "^cascade_(reach_map_national_h13|time_to_invasion|newzones_fanchart|scenario_reach_maps_h13|gateway_knockout|corridors)", full.names = TRUE),
              list.files(figdir, pattern = "^(urban_impact_summary|urban_timing_sensitivity|urban_case_burden|urban_seed_r_sweep)", full.names = TRUE))
for (f in ko_files) {
  .dst <- file.path(KO, basename(f))
  if (file.exists(f)) file.copy(f, .dst, overwrite = TRUE)
  # A bundle entry whose source this run did not produce is STALE — remove the old copy rather
  # than leave last run's table in key_outputs beside this run's.
  else if (file.exists(.dst)) {
    file.remove(.dst)
    message("[cascade] key_outputs: removed stale ", basename(.dst), " (not produced by this run)")
  }
}

# THE CASCADE FIGURE SCRIPTS WRITE INTO THE CANONICAL TREE, whatever CASCADE_OUT_ROOT says:
# 37 targets key_outputs/figures/ and 40 targets key_outputs/manuscript_figures/Figure4 —
# and 40 DELETES the existing Figure4 when a panel fails. An onset-imputation arm must
# therefore not run them: arm 7's Figure 4 is not a deliverable, and twenty arms in
# succession would leave the published figure describing one arbitrary completion. The MI
# figure is drawn once, from the POOLED tables, by 48/49. Default is on, so an ordinary
# cascade run is unchanged.
.run_cascade_figs <- !identical(trimws(Sys.getenv("CASCADE_RUN_FIGURES", "1")), "0")
if (.run_cascade_figs) {
  # publication-style Figure 4 (cascade): relative-risk map + top-20 ranked reach,
  # written to key_outputs/figures/Figure4_cascade_h{4,8,13}.{pdf,png} (+ data CSVs).
  FIG4C_OK <- tryCatch({ source(file.path(ST_DIR, "37_cascade_figure4.R")); TRUE },
           error = function(e) { message("[cascade] Figure4_cascade skipped: ", conditionMessage(e))
                                 FALSE })

  # "next dominoes" narrative figure (main-text Figure 4): anchor frontier map + city ranking
  # + per-city watch-lists, written to key_outputs/manuscript_figures/Figure4.{pdf,png} (+ CSV).
  # LAUNCHED AS ITS OWN Rscript, not source()d. 40 is a guarded script (is_script_run()): sourced
  # from here it only defined its helpers and wrote NOTHING, while the final log line below
  # announced "Figure4 (next dominoes)" as produced. A child process is how run_all.R runs it
  # too, and it keeps 40's top-level globals out of this session. Its exit status is reported.
  .rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  .st40 <- tryCatch(
    system2(.rscript, c("--no-save", "--no-restore",
                        shQuote(file.path(ST_DIR, "40_cascade_next_dominoes.R"))),
            env = c(sprintf("ANALYSIS_DATE=%s", format(as.Date(ANALYSIS_DATE))),
                    sprintf("CASCADE_ST_DIR=%s", ST_DIR),
                    # This run's id: the stamp still reads "running" here, and 40 accepts that
                    # only for the run that launched it (cascade_live_run_id(), 35_cascade_viz.R).
                    # 37, source()d above, sees CASCADE_RUN_ID as a global instead.
                    sprintf("CASCADE_RUN_ID=%s", CASCADE_RUN_ID))),
    error = function(e) { message("[cascade] Figure4 (next dominoes) could not start: ",
                                  conditionMessage(e)); 127L })
  FIG4_OK <- identical(as.integer(.st40), 0L)
  if (!FIG4_OK)
    warning(sprintf("[cascade] Figure4 (next dominoes) FAILED (exit %s) — see its output above.",
                    format(.st40)), call. = FALSE, immediate. = TRUE)
} else {
  FIG4C_OK <- FIG4_OK <- NA
  message("[cascade] CASCADE_RUN_FIGURES=0 — Figure4_cascade and Figure4 NOT drawn ",
          "(tables only; this is the onset-imputation arm path).")
}
.write_stamp("complete")
.fig_word <- function(ok) if (is.na(ok)) "not drawn" else if (isTRUE(ok)) "written" else "FAILED"
tick(sprintf("DONE (%s; run %s) — report + %d key_outputs file(s); Figure4_cascade %s; Figure4 (next dominoes) %s",
             CALIB_STATUS, CASCADE_RUN_ID, sum(file.exists(file.path(KO, basename(ko_files)))),
             .fig_word(FIG4C_OK), .fig_word(FIG4_OK)))
