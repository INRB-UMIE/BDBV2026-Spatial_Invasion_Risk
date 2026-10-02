# =============================================================================
# run_cascade.R — 3-MONTH SPATIAL INVASION CASCADE: master driver
# BDBV 2026 DRC · runs modules 30–36 of PLAN_3MONTH_INVASION.md
#
# Usage:  Rscript spatiotemporal/run_cascade.R
# Env flags:
#   CASCADE_N_MC          MC iterations (default 1000; smoke 40)
#   CASCADE_SMOKE=1       fast path (M=40, heavy stages off)
#   CASCADE_LAYER_CACHE   optional rds path to cache/reuse the reconstructed layer
#   CASCADE_RUN_BACKTEST / CASCADE_RUN_SBC / CASCADE_RUN_SENS  (default 1; smoke 0)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(here) })
ST_DIR <- Sys.getenv("CASCADE_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
t0 <- Sys.time(); tick <- function(m) message(sprintf("[cascade t=%.0fs] %s",
                                     as.numeric(Sys.time() - t0, units = "secs"), m))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

source(file.path(ST_DIR, "00_config.R"))
for (f in c("01_data_prep.R","02_epi_params.R","03_mobility_matrices.R","04_nowcasting.R",
            "04b_epinowcast.R","05_baseline_models.R","06_simple_models.R",
            "15_workhorse.R","16_invasion_eval.R","18_ensemble.R",
            "19_spacetime_eval.R","20_forecast_detail.R","21_bayesian_renewal.R",
            "17_invasion_viz.R",
            # 16b supplies fit_invasion_delta()/bootstrap_invasion_delta()/RECAL_BAND, which the
            # cascade calibration reuses (the hazard-scale transform is the same family).
            "16b_invasion_recalibration.R","22_daily_reissue.R",
            "30_projection_config.R","31_source_dynamics.R","32_cascade_simulator.R",
            "33_cascade_eval.R","33b_cascade_calibration.R",
            "34_conditional_queries.R","35_cascade_viz.R","36_report.R",
            "38_urban_scenarios.R","39_cascade_eval_figure.R"))
  source(file.path(ST_DIR, f))
SMOKE <- identical(Sys.getenv("CASCADE_SMOKE"), "1")
RUN_BT   <- !SMOKE && !identical(Sys.getenv("CASCADE_RUN_BACKTEST"), "0")
RUN_SBC  <- !SMOKE && !identical(Sys.getenv("CASCADE_RUN_SBC"), "0")
RUN_SENS <- !SMOKE && !identical(Sys.getenv("CASCADE_RUN_SENS"), "0")
RUN_URBAN <- !identical(Sys.getenv("CASCADE_RUN_URBAN"), "0")
tick(sprintf("sourced (M=%d, backtest=%s sbc=%s sens=%s)", CASCADE_N_MC, RUN_BT, RUN_SBC, RUN_SENS))

# ---- reconstruct modelling layer (cache option) ----------------------------
CACHE <- Sys.getenv("CASCADE_LAYER_CACHE")
if (nzchar(CACHE) && file.exists(CACHE)) {
  layer <- readRDS(CACHE); tick("layer loaded from cache")
} else {
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
  if (nzchar(CACHE)) saveRDS(layer, CACHE)
  tick("layer reconstructed")
}
zones_all <- layer$zones_all; gt_pmfs <- layer$gt_pmfs
province_map <- load_province_map()
zone_province <- setNames(province_map$province, province_map$nom)
vuln <- compute_vulnerability_index(layer$covariates, zones_all, layer$osrm_mat)

# ---- fit covariate-bearing hazards (featured + sensitivity kernels) --------
kernels <- unique(c(CASCADE_KERNEL, if (RUN_SENS) CASCADE_SENS_KERNELS))
kernels <- kernels[kernels %in% names(layer$mobility_matrices)]
iter_fit <- if (SMOKE) 800L else 2000L
designs <- list(); fits <- list()
for (kern in kernels) {
  designs[[kern]] <- build_invasion_design(layer$zone_week_nc, layer$mobility_matrices, gt_pmfs,
                       layer$covariates, layer$osrm_mat, zones_all, mob = kern, gt = CASCADE_GT)
  fits[[kern]] <- fit_bayes_renewal(designs[[kern]], cov_spec = CASCADE_COV_SPEC,
                                    iter = iter_fit, chains = 2L)
}
tick(sprintf("fitted %d kernel hazard(s): %s", length(fits), paste(names(fits), collapse=",")))

# The import term in the CONJUGATE estimator's denominator needs a hazard scale, so a delta is
# resolved first. It comes from 16b's prequential recalibration of the short-horizon model — not
# from the cascade — so there is no circularity: the cascade's own calibration consumes this
# term, it does not produce the delta that built it.
#
# Its reach is now much smaller than it was. The import term feeds estimate_zone_reff(), which
# is a DIAGNOSTIC: the reproduction number the projection transmits on comes from EpiNow2 and
# never passes through delta_imp. The note further down records any disagreement with the
# out-of-sample delta rather than warning, because the two are no longer required to agree.
delta_imp <- tryCatch(cascade_fit_delta(cov_model = paste0("Bayes-", CASCADE_KERNEL, "-geo")),
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
# ONE fitted parameter. delta is fitted OUT OF SAMPLE so that the cascade's expected
# new-invasion count over held-out K-week windows matches the count that actually occurred,
# and is then verified on the same origins. psi (frontier saturation) is no longer fitted:
# it was not identifiable on this outbreak and is fixed at 0 (30_projection_config.R).
# The h=1 refit and 16b's short-horizon factor are both retained as comparisons, so the gap
# between a one-week and a K-week level is reported rather than hidden. Set
# CASCADE_CALIBRATE=0 to skip calibration entirely (delta then falls back to the
# short-horizon factor), e.g. when no held-out origin is reporting-complete yet.
CALIB_K  <- as.integer(Sys.getenv("CASCADE_CALIB_K", unset = "6"))
DO_CALIB <- !SMOKE && !identical(Sys.getenv("CASCADE_CALIBRATE"), "0")
calib <- if (DO_CALIB)
  tryCatch(cascade_calibrate(layer, zone_province = zone_province, K = CALIB_K,
             kernel = CASCADE_KERNEL, gt = CASCADE_GT, seed = CASCADE_SEED,
             # 600, not 400: n_mc_report DEFAULTS TO n_mc_search, so the delta search and the
             # out-of-sample verification use the same Monte-Carlo size and their modelled
             # counts are comparable. They are NOT identical — the verification runs at an
             # offset seed on purpose (33b, "INDEPENDENT SEED"), so the two differ by Monte-
             # Carlo error. Matching the SIZE is what keeps that error small; at 400 vs 600
             # the report printed 72.96 and 71.7 for one quantity.
             n_mc_search = as.integer(Sys.getenv("CASCADE_CALIB_NMC", unset = "600"))),
           error = function(e) { warning("[run_cascade] calibration failed (", conditionMessage(e),
                                 "); falling back to delta-only.", call. = FALSE); NULL }) else NULL
if (is.null(calib))
  calib <- list(delta = cascade_fit_delta(cov_model = paste0("Bayes-", CASCADE_KERNEL, "-geo")),
                psi = CASCADE_PSI, psi_fit = NULL, delta_oos = NULL, delta_prior = NA_real_,
                report = NULL, origins = list(),
                delta_info = list(estimator = "not calibrated this run"),
                note = paste("calibration skipped or failed; delta is 16b's SHORT-HORIZON factor",
                             "and the 13-week level is therefore uncalibrated at the horizon it",
                             "is read at; psi at its configured default"))
delta <- calib$delta
psi   <- calib$psi
# delta_imp (16b's short-horizon factor) built the import term in the CONJUGATE estimator's
# denominator; `delta` is now fitted out of sample against realised K-week counts and drives the
# simulated hazard. They are no longer required to agree, and the old >25% warning would fire on
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
  message(sprintf(paste0("[run_cascade] the out-of-sample %d-week delta (%.3f) differs from the ",
                         "short-horizon factor the import term was built on (%.3f) by %.0f%%. ",
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
# inside a function call. Recorded per zone with the posterior, the shrinkage diagnostics
# and the local/imported split, so a reader can see which zones are data-driven and which
# are carried by the province pool.
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
# non-converged fit is auditable rather than reduced to one scalar) and the
# out-of-sample verification.
.cd <- file.path(OUT_CASCADE, "diagnostics")
if (!is.null(calib$delta_oos)) {
  readr::write_csv(calib$delta_oos$curve, file.path(.cd, "cascade_delta_fit_curve.csv"))
  jsonlite::write_json(calib$delta_oos[setdiff(names(calib$delta_oos), "curve")],
                       file.path(.cd, "cascade_delta_fit.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = 6)
  if (!isTRUE(calib$delta_oos$converged))
    warning("[run_cascade] delta did not converge out of sample — the 13-week layer is NOT fit ",
            "for publication until this is resolved (see cascade_delta_fit.json).", call. = FALSE)
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
  readr::write_csv(calib$origin_diagnostics,
                   file.path(.cd, "cascade_calibration_origin_reff.csv"))
if (!is.null(calib$report)) {
  readr::write_csv(calib$report$per_origin, file.path(.cd, "cascade_calibration_per_origin.csv"))
  readr::write_csv(calib$report$pooled,     file.path(.cd, "cascade_calibration_pooled.csv"))
  if (!is.null(calib$report$reliability))
    readr::write_csv(calib$report$reliability, file.path(.cd, "cascade_calibration_reliability.csv"))
  tick(paste("calibration verification:", calib$report$note))
}
jsonlite::write_json(list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    kernel = CASCADE_KERNEL, calib_K = CALIB_K,
    delta = as.numeric(delta), delta_info = calib$delta_info, delta_sd_log = delta_sd_log,
    delta_short_horizon = calib$delta_prior,
    delta_estimator = if (is.null(calib$delta_oos)) "short-horizon factor (no OOS fit)"
                      else sprintf("out-of-sample %d-week count match", CALIB_K),
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
    verification = if (is.null(calib$report)) NULL else as.list(calib$report$pooled),
    # NOT "verification": the pooled count is the quantity delta was fitted to, so its residual
    # is in-sample by construction (see cascade_calibration_report). Recorded under a name that
    # says what it is, with an explicit flag, so nothing downstream can read it as validation.
    count_residual_ok        = if (is.null(calib$report)) NA else isTRUE(calib$report$count_residual_ok),
    count_level_in_sample    = if (is.null(calib$report)) NA else isTRUE(calib$report$count_level_in_sample),
    note = calib$note),
  file.path(.cd, "cascade_calibration.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6, null = "null")

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
  # new provinces reached by 13w (provinces with no affected zone at t0)
  # DROP unmapped zones rather than counting them as a province. prov_of() returns NA for any
  # zone absent from province_map; NA survives setdiff() unless it is already in
  # affected_provinces0, so a single unmapped zone being reached added a phantom "newly reached
  # province" in every iteration where it was reached. cascade_flux_network() guards the same
  # case (34_conditional_queries.R:133); this call did not.
  .provs <- function(z) { pv <- prov_of(z); unique(pv[!is.na(pv)]) }
  newprov <- vapply(seq_len(sim$n_mc), function(m) {
    zr <- zones_all[!is.na(sim$tau[, m]) & sim$tau[, m] <= 13L]
    length(setdiff(.provs(zr), affected_provinces0)) }, numeric(1))
  # ONE consistent estimator throughout: the MEDIAN of the cumulative-new-zone
  # distribution at each horizon (+ 90% predictive band), so 4/8/13w are comparable
  # (means and medians differ for the right-skewed count; do not mix them).
  q <- function(h, p) stats::quantile(cum[h, ], p, names = FALSE)
  # Horizons from CASCADE_REPORT_HORIZONS, not the literals 4/8/13: hard-coded row indices
  # abort with a subscript error (or silently report the wrong weeks) if the reporting
  # horizons or CASCADE_HORIZON_WEEKS change — and only AFTER every scenario has been simulated.
  .h <- sort(unique(as.integer(CASCADE_REPORT_HORIZONS)))
  .h <- .h[.h >= 1L & .h <= nrow(cum)]
  .hmax <- max(.h)
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
fan <- cascade_fig_fanchart(sims)
tryCatch({
  if (!is.null(fan) && nrow(fan)) {
    .fanp <- file.path(OUT_CASCADE, "tables", "cascade_newzones_fanchart.csv")
    dir.create(dirname(.fanp), recursive = TRUE, showWarnings = FALSE)
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

conditional <- list()
for (hub in head(gateway$health_zone, 2L)) {
  cs <- cascade_conditional_seed(prep, CASCADE_SCENARIOS[[scP_key]], hub, delta, psi = psi,
                                 n_mc = min(CASCADE_N_MC, if (SMOKE) 40L else 800L),
                                 pop_vec = layer$pop)
  rs_c <- cascade_enrich_reach(cs$reach, province_map, vuln, method_label = paste0("cond-", hub))
  cascade_map_conditional(rs_c, hub)
  r13 <- cs$reach[cs$reach$horizon == 13L & !cs$reach$was_active_before, ]
  conditional[[hub]] <- list(top = head(r13$health_zone[order(-r13$p_invasion)], 6))
  readr::write_csv(rs_c, file.path(TAB, sprintf("cascade_conditional_%s.csv",
                                    gsub("[^A-Za-z0-9]+","_",hub))))
}
tick("conditional hubs written")

# ---- urban-hub invasion scenarios (timing-aware) ---------------------------
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
sens <- if (RUN_SENS) tryCatch(cascade_sensitivity(layer, designs, fits, reff, delta,
                        base_reach = reaches[[scP_key]], psi = psi),
                        error = .stage_failed("sensitivity")) else NULL
if (!is.null(sens)) readr::write_csv(sens, file.path(OUT_CASCADE, "diagnostics", "cascade_sensitivity.csv"))

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
  ksweep_meta <- list(kernel = CASCADE_KERNEL, M = CASCADE_N_MC,
                      k_grid = CASCADE_K_INDIV_SWEEP, k_base = CASCADE_K_INDIV,
                      delta = as.numeric(delta), cal_in_large = attr(delta, "cal_in_large") %||% NA,
                      psi = psi,
                      psi_converged = if (is.null(psi_fit)) NA else isTRUE(psi_fit$converged),
                      psi_boundary_hit = if (is.null(psi_fit)) NA else isTRUE(psi_fit$boundary_hit),
                      psi_horizon_pace_ratio = if (is.null(psi_fit)) NA else psi_fit$horizon_pace_ratio,
                      gate_spearman = gate$spearman_signal, gate_n_signal = gate$n_signal,
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
backtest <- if (RUN_BT) tryCatch(cascade_backtest(layer, delta = delta, psi = psi, K = CALIB_K,
                        zone_province = zone_province, linelist = layer$ll,
                        origins = if (length(calib$origins)) calib$origins else NULL),
                        error = .stage_failed("backtest")) else NULL
if (!is.null(backtest)) {
  readr::write_csv(backtest, file.path(OUT_CASCADE, "diagnostics", "cascade_backtest.csv"))
  # per-zone detail (drives the Figure2_cascade capture-curve / mean-rank / top-K panels)
  bt_detail <- attr(backtest, "detail")
  if (!is.null(bt_detail)) readr::write_csv(bt_detail,
    file.path(OUT_CASCADE, "diagnostics", "cascade_backtest_detail.csv"))
}
sbc <- if (RUN_SBC) tryCatch(cascade_sbc(designs[[CASCADE_KERNEL]], n_sbc = 20L),
                        error = .stage_failed("sbc")) else NULL
tick("validation done")

# ---- report + key_outputs bundle -------------------------------------------
ctx <- list(reach_primary = rs_primary, n_mc = CASCADE_N_MC, delta = delta, psi = psi,
            calib = calib, delta_sd_log = delta_sd_log,
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

# key_outputs bundle (mirror the short-horizon key_outputs style)
KO <- file.path(OUT_DIR, "key_outputs")
if (!dir.exists(KO)) dir.create(KO, recursive = TRUE)
figdir <- file.path(OUT_CASCADE, "figures")
ko_files <- c(file.path(TAB, "cascade_reach_scores_all_zones.csv"),
              file.path(TAB, "cascade_gateway_knockout.csv"),
              file.path(TAB, "urban_scenario_summary.csv"),
              list.files(figdir, pattern = "^cascade_(reach_map_national_h13|time_to_invasion|newzones_fanchart|scenario_reach_maps_h13|gateway_knockout|corridors)", full.names = TRUE),
              list.files(figdir, pattern = "^(urban_impact_summary|urban_timing_sensitivity|urban_case_burden|urban_seed_r_sweep)", full.names = TRUE))
for (f in ko_files) if (file.exists(f)) file.copy(f, file.path(KO, basename(f)), overwrite = TRUE)

# publication-style Figure 4 (cascade): relative-risk map + top-20 ranked reach,
# written to key_outputs/figures/Figure4_cascade_h{4,8,13}.{pdf,png} (+ data CSVs).
tryCatch(source(file.path(ST_DIR, "37_cascade_figure4.R")),
         error = function(e) message("[cascade] Figure4_cascade skipped: ", conditionMessage(e)))

# "next dominoes" narrative figure (main-text Figure 4): anchor frontier map + city ranking
# + per-city watch-lists, written to key_outputs/manuscript_figures/Figure4.{pdf,png} (+ CSV).
tryCatch(source(file.path(ST_DIR, "40_cascade_next_dominoes.R")),
         error = function(e) message("[cascade] Figure4 (next dominoes) skipped: ", conditionMessage(e)))
tick(sprintf("DONE — report + %d key_outputs files + Figure4_cascade + Figure4 (next dominoes)", length(ko_files)))
