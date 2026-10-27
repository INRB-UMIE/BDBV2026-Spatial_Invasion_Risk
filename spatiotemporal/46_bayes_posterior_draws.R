# =============================================================================
# 46_bayes_posterior_draws.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
#
# Rebuild outputs/key_outputs/bayes_posterior_draws.csv — the FULL posterior draws of the
# Bayesian invasion suite's parameters — WITHOUT a full run_all.R run.
#
# run_all.R writes this artefact from the fits it already has, and that is the normal way it
# comes to exist (see the block above `bayes_posterior_draws.csv` there). This script exists
# for the other case: the draws are missing or stale but the rest of the run's outputs are
# current, and refitting the whole suite inside a pipeline run to recover one CSV is not
# proportionate. It refits ONLY the models asked for (default: the featured one), each on its
# own full design, as fit_bayes_suite() does; a time-varying form is refitted by the Stan engine.
#
# IT REFUSES TO WRITE A FILE THAT DISAGREES WITH THE PUBLISHED TABLE.
# Every model it refits is compared against its row in outputs/reports/bayes_parameters.csv,
# and the CSV is written only if all of them reproduce it within REPRO_TOL. That gate is not
# ceremony — it is the whole reason this script is safe to have. Reproducing the suite's
# posterior means rebuilding its TRAINING DESIGN exactly, and that design is quietly fragile:
#
#   * The rolling as-of predictors are gated on epinow2_truncation(regime = "asof"), which
#     runs .trunc_fit_asof_cohort(); that function's FIRST line returns NULL — silently, with
#     no warning and no message — when linelist_observation_date() is not in scope. It lives
#     in 22_daily_reissue.R, which has nothing else to do with the Bayesian suite. Omit that
#     one source() and every model is refit with ROLLING PREDICTORS OFF.
#   * That failure is close to invisible. bayes_forecast_offsets() never takes the rolling
#     arguments, so the FORECAST import force stays bit-identical to the published one
#     (checked zone by zone on the run where this was found: ratio 1.0000 across every at-risk
#     zone). Only the training design moves — on that run 7,595 rows / 39 events became
#     8,609 / 53 — and with it every coefficient: the then-featured model's ccvi posterior
#     median went 2.39 -> 1.67, its intercept 0.103 -> 0.137. Both look perfectly reasonable
#     on their own.
#
# So the check below is what stands between "these are the model's posteriors" and "these are
# some model's posteriors". Keep it.
#
# Usage:   Rscript 46_bayes_posterior_draws.R
#          BAYES_DRAWS_MODELS="Bayes-M14-fill-geo,Bayes-M8-fill-geo" Rscript 46_bayes_posterior_draws.R
#          BAYES_DRAWS_MODELS=all Rscript 46_bayes_posterior_draws.R     # the whole suite
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr) })

ST <- file.path(here::here(), "spatiotemporal")
# THE SOURCE LIST IS LOAD-BEARING, NOT INCIDENTAL. 22_daily_reissue.R in particular is
# required for the as-of truncation to resolve at all; see the header. Keep it in step with
# the modules run_all.R has loaded by the time it calls fit_bayes_suite().
# ONE module list, from 00_config.R. Nine files used to hand-maintain their own and all nine
# differed — see load_pipeline_modules().
source(file.path(ST, "00_config.R"))
stop_if_data_definition_sensitivity("46_bayes_posterior_draws.R")   # hard-coded published paths
load_pipeline_modules("model", dir = ST)


# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("46_bayes_posterior_draws.R")) {
  message("[bayes-draws] sourced, not run: nothing is written. Use `Rscript 46_bayes_posterior_draws.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

OUT         <- file.path(ST, "outputs")
DRAWS_PATH  <- file.path(OUT, "key_outputs", "bayes_posterior_draws.csv")
PARAMS_PATH <- file.path(OUT, "reports", "bayes_parameters.csv")
# Relative tolerance, applied to the posterior median AND both 90% CrI endpoints.
#
# THE EXPECTED DEVIATION IS EXACTLY ZERO. fit_bayes_renewal() seeds from RANDOM_SEED, so a
# refit on the same design is the same chain draw for draw; rebuilding the featured model this
# way reproduces its published median and both endpoints to 0.000%. The tolerance is not there
# to absorb sampling noise — it is a margin for a different cmdstan or toolchain version, which
# can perturb the last digits without anything being wrong. A rebuilt-design error is nothing
# like that size: the rolling-predictor failure in the header moves ccvi by 30%.
#
# So read a PASS with a non-zero deviation as a result in itself, not as a rounding artefact,
# and go and find out what moved. The pass message prints the number for that reason.
REPRO_TOL <- 0.05

if (!file.exists(PARAMS_PATH))
  stop("[draws] ", basename(PARAMS_PATH), " is missing. This script validates every refit ",
       "against the published parameter table and will not write without it; run the ",
       "pipeline instead.", call. = FALSE)
published <- read_csv(PARAMS_PATH, show_col_types = FALSE)

# WHICH MODELS. The featured one by default — it is the only one any published figure draws
# the posterior of — read from the pipeline's own record rather than named here, through the
# suite's ONE reader, read_featured_model() (figure_style.R; sourced into a private environment
# so its palette and theme objects stay out of this session). It errors, never guesses, when the
# record is missing or names no model.
read_featured_model <- local({
  .fs <- new.env()
  sys.source(file.path(ST, "figure_style.R"), envir = .fs)
  .fs$read_featured_model
})
FEATURED <- read_featured_model(file.path(OUT, "key_outputs", "model_selection.json"))
.want <- trimws(Sys.getenv("BAYES_DRAWS_MODELS", ""))
WANT <- if (!nzchar(.want)) {
  FEATURED
} else if (identical(tolower(.want), "all")) {
  NA_character_                                   # resolved once the grid is known, below
} else {
  trimws(strsplit(.want, ",")[[1]])
}

# ---------------------------------------------------------------------------
# 1. The suite's inputs, rebuilt from saved artefacts
# ---------------------------------------------------------------------------
message("\n=== 46: rebuilding the Bayesian suite's fitting context ===")
dat        <- prep_all_data()
zones_all  <- names(dat$pop)
covariates <- dat$covariates
gt_pmfs    <- compute_all_gt_pmfs()
osrm_mat   <- load_osrm()

# Nowcast — the DEPLOYED (extract-regime) truncation, exactly as run_all.R step 4.
.trunc_deployed <- tryCatch(
  .trunc_as_delay_spec(epinow2_truncation(issue_date = NULL, ll = dat$ll), "extract"),
  error = function(e) { warning("[draws] deployed truncation unavailable (",
                                conditionMessage(e), ")", call. = FALSE); NULL })
zone_week_nc <- apply_nowcast_correction(dat$zone_week, analysis_date = ANALYSIS_DATE,
                                         delay = .trunc_deployed)

# As-of truncation — what the ROLLING predictors ride on. A NULL here silently turns them off
# and changes every coefficient (see the header), so it is checked rather than tolerated.
.trunc_asof <- tryCatch(
  .trunc_as_delay_spec(epinow2_truncation(regime = "asof", ll = dat$ll), "asof"),
  error = function(e) NULL)
.roll_on <- isTRUE(get0("ROLLING_PREDICTORS", ifnotfound = TRUE)) && !is.null(.trunc_asof)
if (isTRUE(get0("ROLLING_PREDICTORS", ifnotfound = TRUE)) && is.null(.trunc_asof))
  stop("[draws] ROLLING_PREDICTORS is on but the as-of truncation did not resolve, so every ",
       "model would be refit WITHOUT the rolling predictors and would not reproduce the ",
       "published posterior. The usual cause is linelist_observation_date() not being in ",
       "scope (22_daily_reissue.R). Refusing to fit.", call. = FALSE)
message(sprintf("[draws] rolling as-of predictors: %s", if (.roll_on) "ON" else "off"))

# Mobility: the CACHED matrices this run already built. Rebuilding them takes minutes and is
# deterministic given the same inputs, so reading them back is the same kernel, faster.
#
# ONLY THE KERNELS THE RUN RECORDED BUILDING. outputs/mobility/ also holds files left by earlier
# runs and other branches (kernels since switched off: -split, -dist, M2a/M2b, M6b, M7, ...),
# and every file read here is a candidate grid member: bayes_default_grid() keeps a spec when its
# kernel is PRESENT, so a stale file whose switch is on would add a model to the grid (and a name
# to BAYES_DRAWS_MODELS=all) that the run never fitted.
# mobility_manifest.csv (03_mobility_matrices.R) is the run's own record of what it wrote, with a
# config fingerprint. Only the manifest's files are read, the manifest must describe ONE build,
# and that build must be for this analysis date and this session's governing settings — the
# same fields 03 fingerprints, compared as the human-readable `config` string it writes (a
# field added to 03's fingerprint makes this refuse until it is added here, never pass).
.mob_dir <- file.path(OUT, "mobility")
.man_f <- file.path(.mob_dir, "mobility_manifest.csv")
if (!file.exists(.man_f))
  stop("[draws] ", .man_f, " is missing, so the current kernels cannot be told from stale ",
       "files in ", .mob_dir, "; run the pipeline (03_mobility_matrices.R) first.", call. = FALSE)
.man <- read_csv(.man_f, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
stopifnot(all(c("matrix_id", "file", "config", "config_hash") %in% names(.man)), nrow(.man) > 0L)
if (dplyr::n_distinct(.man$config_hash) != 1L)
  stop("[draws] mobility_manifest.csv records more than one build configuration.", call. = FALSE)
# The fields are mirrored here one for one, in 03's order. The input-identity values (content
# hashes, the cohort origin lists) come from 03's own mobility_input_fingerprint() on this
# session's inputs — the same population vector and OSRM matrices run_all.R hands the build, the
# road-distance one by run_all.R's rule (load_osrm_dist_for_mobility()) — and the string is
# formatted by 03's mobility_fingerprint_string(), so both sides compute every value with ONE code.
.inp_now <- mobility_input_fingerprint(zones_all, dat$pop, ANALYSIS_DATE, osrm_mat,
                                       load_osrm_dist_for_mobility())
.cfg_now <- list(
  analysis_date        = format(ANALYSIS_DATE),
  flowminder_od_file   = as.character(get0("FLOWMINDER_OD_FILE",   ifnotfound = NA)),
  mobility_source_fill = as.character(get0("MOBILITY_SOURCE_FILL", ifnotfound = NA)),
  mobility_fill_qmax   = as.character(get0("MOBILITY_FILL_QMAX",   ifnotfound = NA)),
  cohort_window        = as.character(get0("COHORT_WINDOW",        ifnotfound = NA)),
  osrm_gap_fill        = as.character(isTRUE(get0("OSRM_GAP_FILL",                 ifnotfound = TRUE))),
  m4_censored_fit      = as.character(isTRUE(get0("M4_CENSORED_FIT",               ifnotfound = TRUE))),
  cohort_models        = as.character(isTRUE(get0("INCLUDE_COHORT_MODELS",         ifnotfound = TRUE))),
  cohort_gravity       = as.character(isTRUE(get0("INCLUDE_COHORT_GRAVITY_MODELS", ifnotfound = TRUE))),
  sourcefill_models    = as.character(isTRUE(get0("INCLUDE_SOURCEFILL_MODELS",     ifnotfound = TRUE))),
  cohort_sources       = .inp_now$cohort_sources,
  epicentre_zones      = .inp_now$epicentre_zones,
  population           = .inp_now$population,
  osrm_travel_time     = .inp_now$osrm_travel_time,
  osrm_road_distance   = .inp_now$osrm_road_distance,
  flowminder_od        = .inp_now$flowminder_od,
  short_trip           = .inp_now$short_trip,
  cohort_profiles      = .inp_now$cohort_profiles)
.cfg_now_str <- mobility_fingerprint_string(.cfg_now)
if (!identical(unique(.man$config), .cfg_now_str))
  stop("[draws] the cached kernels were built under\n  ", unique(.man$config),
       "\nbut this session's settings are\n  ", .cfg_now_str,
       "\nso they are not the kernels the published suite was fitted on. Rebuild them ",
       "(03_mobility_matrices.R) first.", call. = FALSE)
.mf <- .man$file
.absent <- .mf[!file.exists(file.path(.mob_dir, .mf))]
if (length(.absent))
  stop("[draws] the manifest lists kernel file(s) that are missing: ",
       paste(.absent, collapse = ", "), call. = FALSE)
mobility_matrices <- setNames(lapply(file.path(.mob_dir, .mf), readRDS), .man$matrix_id)
.stale <- setdiff(list.files(.mob_dir, "^mobility_.*\\.rds$"), .mf)
message(sprintf("[draws] %d cached mobility matrices from the manifest (%d stale file(s) in %s ignored)",
                length(mobility_matrices), length(.stale), basename(.mob_dir)))

# ---------------------------------------------------------------------------
# 2. The deployed suite's grid, and the models to refit
# ---------------------------------------------------------------------------
# THE GRID IS run_all.R's DEPLOYED SUITE GRID: the cross-validated grid (gt_arms = FALSE), the
# time-varying forms of the stage-1 pick that were cross-validated (model_selection.json,
# featured$bayesian$stage2_pool), and the generation-time arms of THIS run's featured model
# (bayes_gt_sensitivity_specs()). fit_bayes_suite() fits EVERY model on its OWN full design
# (only the stacking LOO uses the rows common to the stacked models), so a model is reproduced
# from its own design alone and only the models asked for are built.
# NOT MIRRORED: run_all.R drops the arms when their cross-validation fails to cover the shared
# support; if that happened in the published run, bayes_parameters.csv has no -gt rows and the
# reproduction check below fails rather than passing on a different model.
.sel <- jsonlite::fromJSON(file.path(OUT, "key_outputs", "model_selection.json"), simplifyVector = TRUE)
.lfo_grid <- bayes_default_grid(mobility_matrices, gt_arms = FALSE)
.pool <- as.character(unlist(.sel$featured$bayesian$stage2_pool))
.s1 <- bayes_spec_for_label(.sel$featured$bayesian$stage1_method %||% NA_character_, .lfo_grid)
.form_specs <- Filter(function(g) g$label %in% .pool, bayes_tv_forms_of(.s1))
.feat_spec <- bayes_spec_for_label(FEATURED, c(.lfo_grid, .form_specs))
.gt_arm_specs <- if (isTRUE(get0("INCLUDE_GT_SENSITIVITY_MODELS", ifnotfound = TRUE)) &&
                     !is.null(.feat_spec)) bayes_gt_sensitivity_specs(.feat_spec) else list()
grid <- Filter(function(g) !is.null(mobility_matrices[[g$mob]]), c(.lfo_grid, .form_specs, .gt_arm_specs))
.labels <- vapply(grid, `[[`, "", "label")
if (all(is.na(WANT))) WANT <- .labels
WANT <- unique(WANT)          # a repeated name would fit twice and duplicate its draws
.missing <- setdiff(WANT, .labels)
if (length(.missing))
  stop("[draws] not in the current grid: ", paste(.missing, collapse = ", "), call. = FALSE)
if (any(vapply(grid[match(WANT, .labels)], .spec_is_stan, logical(1))) && !.roll_on)
  stop("[draws] a time-varying form is requested but the rolling as-of design is off.", call. = FALSE)

# ---------------------------------------------------------------------------
# 3. Refit, validate, write
# ---------------------------------------------------------------------------
draws <- list(); checks <- list()
for (nm in WANT) {
  g <- grid[[match(nm, .labels)]]
  des <- build_invasion_design(
           zone_week_nc, mobility_matrices, gt_pmfs, covariates, osrm_mat, zones_all,
           mob = g$mob, gt = g$gt,
           candidates    = unique(c(g$cov, "log_pop", "ccvi", "d_min")),
           rolling_ll    = if (.roll_on) dat$ll else NULL,
           rolling_delay = if (.roll_on) .trunc_asof else NULL)
  message(sprintf("[draws] fitting %s — %d rows, %d events, covariates: %s%s",
                  nm, des$n_obs, des$n_events,
                  if (length(g$cov)) paste(g$cov, collapse = " + ") else "none",
                  if (.spec_is_stan(g)) paste0("; Stan form ", g$form) else ""))
  .tv <- g$tv %||% "none"
  fit <- if (.spec_is_stan(g)) {
           ex <- invasion_exposure(zone_week_nc, dat$ll, zones_all, mobility_matrices, gt_pmfs,
                                   covariates, osrm_mat, g$mob, g$gt, intersect(g$cov, des$feat),
                                   des$center, des$scale, .trunc_asof)
           stan_invasion_fit(des, ex, zones_all, g$cov, g$use_rw, g$use_frailty, iter = 2000L,
                             label = nm)
         } else if (.tv_is_structural(.tv))
           fit_bayes_renewal_tv(des, cov_spec = g$cov, iter = 2000L,
                                link = g$link %||% "cloglog", tv = .tv)
         else
           fit_bayes_renewal(des, cov_spec = g$cov, iter = 2000L,
                             link = g$link %||% "cloglog")
  if (is.null(fit)) stop("[draws] ", nm, " did not fit.", call. = FALSE)
  # EVERY DRAW (see the matching note in run_all.R): the published 90% interval is computed
  # from all of them, so thinning here would put this artefact's interval up to ~2.5% away
  # from the table's for the same parameter — and would blunt the gate below, which is sharp
  # enough with full draws to detect a rebuilt design as an EXACT zero deviation or not at all.
  d <- bayes_posterior_draws(setNames(list(fit), nm), max_draws = .Machine$integer.max)
  if (is.null(d) || !nrow(d))
    stop("[draws] ", nm, " fitted but yielded no posterior draws.", call. = FALSE)
  draws[[nm]] <- d
  # A missing prior stamp is not fatal — the figure simply draws that posterior alone — but it
  # must not pass unremarked, because the usual cause is the centring check in
  # bayes_posterior_draws() having refused the intercept prior.
  .nop <- d %>% dplyr::distinct(term, prior) %>% dplyr::filter(is.na(prior))
  if (nrow(.nop))
    message("[draws] ", nm, ": no prior stamped for ", paste(.nop$term, collapse = ", "),
            " — any prior-vs-posterior figure will draw those posteriors alone.")

  # VALIDATE THE WHOLE MARGINAL, NOT JUST ITS CENTRE. The artefact's purpose is to carry the
  # SHAPE of each posterior, so a gate on the median alone would wave through a refit that
  # lands the right centre with the wrong spread — precisely the error a density figure would
  # then publish. Median and both 90% CrI endpoints are checked against the published table.
  .pub_m <- published %>% dplyr::filter(model == nm)
  .refit <- d %>%
    dplyr::group_by(model, term) %>%
    dplyr::summarise(med = stats::median(hr),
                     lo  = stats::quantile(hr, 0.05, names = FALSE),
                     hi  = stats::quantile(hr, 0.95, names = FALSE), .groups = "drop")
  # Every fitted coefficient must HAVE a published counterpart. inner_join() alone would drop
  # an unmatched term silently, leaving it written to the artefact but never validated. (The
  # converse is expected and fine: the table also carries the brms sd_ terms and an RW1 form's
  # beta_last_training_week level, which are not in the draws. A form's frailty_theta /
  # rw1_sigma_week ARE in both and are checked like the coefficients.)
  .unmatched <- setdiff(.refit$term, .pub_m$term)
  if (length(.unmatched))
    stop("[draws] ", nm, ": refitted parameter(s) ", paste(.unmatched, collapse = ", "),
         " have no row in ", basename(PARAMS_PATH), ", so they cannot be validated. ",
         "NOTHING WAS WRITTEN.", call. = FALSE)
  checks[[nm]] <- .refit %>%
    dplyr::inner_join(.pub_m %>% dplyr::select(model, term, p_med = hr, p_lo = lo, p_hi = hi),
                      by = c("model", "term")) %>%
    dplyr::mutate(rel_med = abs(med - p_med) / p_med,
                  rel_lo  = abs(lo  - p_lo)  / p_lo,
                  rel_hi  = abs(hi  - p_hi)  / p_hi,
                  rel     = pmax(rel_med, rel_lo, rel_hi))
}

chk <- dplyr::bind_rows(checks)
if (!nrow(chk))
  stop("[draws] none of the refitted models appears in ", basename(PARAMS_PATH),
       "; there is nothing to validate against, so nothing is written.", call. = FALSE)
cat("\n--- refit vs published: posterior median and 90% CrI endpoints ---\n")
print(as.data.frame(chk), digits = 5)
.bad <- chk %>% dplyr::filter(!is.finite(rel) | rel > REPRO_TOL) %>% dplyr::arrange(dplyr::desc(rel))
if (nrow(.bad)) {
  .w1 <- .bad[1, ]
  # Name the quantity that actually failed, not just the parameter: a median that matches with
  # an interval that does not is a different diagnosis from both being off.
  .which <- c("median", "5% CrI", "95% CrI")[
    which.max(c(.w1$rel_med, .w1$rel_lo, .w1$rel_hi))]
  .got <- c(.w1$med, .w1$lo, .w1$hi)[which.max(c(.w1$rel_med, .w1$rel_lo, .w1$rel_hi))]
  .exp <- c(.w1$p_med, .w1$p_lo, .w1$p_hi)[which.max(c(.w1$rel_med, .w1$rel_lo, .w1$rel_hi))]
  stop(sprintf(paste0("[draws] the refit does NOT reproduce %s: %d of %d parameters differ by ",
                      "more than %.0f%% (worst: %s %s, %s %.4g against the published %.4g). ",
                      "The training design has been rebuilt differently from the one the suite ",
                      "fitted — see the header for the usual cause. NOTHING WAS WRITTEN."),
               basename(PARAMS_PATH), nrow(.bad), nrow(chk), 100 * REPRO_TOL,
               .w1$model, .w1$term, .which, .got, .exp), call. = FALSE)
}
message(sprintf("[draws] reproduction check passed: %d parameter(s), max deviation %.2f%%",
                nrow(chk), 100 * max(chk$rel)))

out <- dplyr::bind_rows(draws) %>% dplyr::mutate(hr = signif(hr, 8))
# MERGE, not overwrite: the run's file carries every suite model, and a refit of the featured
# one must not silently drop the others. Only the refitted models' rows are replaced; a file
# written before a column existed is refused rather than mixed.
if (file.exists(DRAWS_PATH)) {
  .old <- readr::read_csv(DRAWS_PATH, show_col_types = FALSE)
  if (!setequal(names(.old), names(out)))
    stop("[draws] the existing ", basename(DRAWS_PATH), " has columns (",
         paste(names(.old), collapse = ", "), ") that differ from the refit's (",
         paste(names(out), collapse = ", "), "); delete it or re-run the pipeline. NOTHING WAS WRITTEN.",
         call. = FALSE)
  .kept <- .old[!(.old$model %in% unique(out$model)), names(out), drop = FALSE]
  if (nrow(.kept))
    message(sprintf("[draws] keeping %d other model(s) from the existing file: %s",
                    dplyr::n_distinct(.kept$model), paste(unique(.kept$model), collapse = ", ")))
  out <- dplyr::bind_rows(.kept, out)
}
readr::write_csv(out, DRAWS_PATH)
message(sprintf("[draws] wrote %s — %d model(s) x %d parameter(s), %d draws each",
                DRAWS_PATH, dplyr::n_distinct(out$model),
                dplyr::n_distinct(paste(out$model, out$term)), max(out$.draw)))

}   # end of the script-run guard

