# =============================================================================
# 33b_cascade_calibration.R — CALIBRATION of the 13-week invasion cascade
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# WHAT THIS FILE DOES — the procedure cascade_calibrate() runs, in order
# ----------------------------------------------------------------------
# The cascade has ONE fitted level parameter, delta, a multiplier on the per-week seeding hazard
# (mu -> delta * mu, i.e. p -> 1 - (1 - p)^delta). Frontier saturation psi
# (sat = exp(-psi * max(f_inv - f_inv0, 0))) is NOT fitted: it is held at CASCADE_PSI (0, so sat
# is identically 1; 30_projection_config.R) because it was not identifiable on this outbreak (see
# cascade_fit_delta_oos()). cascade_fit_psi_oos() is kept for sensitivity sweeps only.
#
#   1. Short-horizon factor (`delta_prior`). cascade_fit_delta() (30_projection_config.R) reads
#      the prequential maximum-likelihood h = 1 hazard-scale factor that
#      16b_invasion_recalibration.R fits for the cascade's hazard, cascade_hazard_label() =
#      Bayes-<kernel>-geo[-<form suffix>] (the featured form's, falling back to its fixed-beta
#      parent), with its zone-clustered bootstrap interval and se_log, and falls back
#      with a warning at each step when that row is missing. The projection does NOT run on it
#      unless step 4 fails. It is used as: the hazard scale of the import term in each origin's
#      conjugate R diagnostic (`delta_imp`); the reference of the week-1 transfer check (step 3);
#      and the source of `se_log`, the spread delta is drawn with (step 5).
#
#   2. Origins (cascade_calibration_origins()). Cutoffs `cutoffs_from_end` weeks before the last
#      week of the record, skipping any inside the ascertainment burn-in (CASCADE_CALIB_MIN_WEEKS),
#      any with no settled outcome window, and any with fewer than 3 training events. At each,
#      the training counts, the hazard (rolling as-of design) and the national EpiNow2 R are
#      rebuilt as of that origin (LEAKAGE, below); the realised outcomes at h = 1 and h = K are
#      read from the FINAL record, each horizon only once its window is SETTLED
#      (LFO_RELIABLE_EVAL_AGE_DAYS). "Held-out" means held out from each origin's TRAINING and
#      nothing more: step 4 fits delta to these same outcomes.
#
#   3. Week-1 transfer check (cascade_delta_h1_refit()) — a DIAGNOSTIC that sets nothing. One week
#      is simulated from each origin at delta = 1 (psi is irrelevant at h = 1: no zone has been
#      seeded yet, so sat = 1) and the factor that calibrates those week-1 probabilities to the
#      realised week-1 invasions is fitted with 16b's own estimator and zone-clustered bootstrap
#      (not attempted below 10 week-1 events). It is recorded against BOTH references:
#      `h1_refit_ratio_short` (to the step-1 factor; warned beyond `transfer_tol`) and
#      `h1_refit_ratio_deployed` (to the factor step 4 deploys). The two differ whenever the
#      one-week and K-week levels do, so neither may be quoted as the other.
#
#   4. delta, THE ONE FITTED PARAMETER (cascade_fit_delta_oos()). Geometric bisection over
#      RECAL_BAND for the delta at which the cascade's pooled expected number of new invasions over
#      the K-week window (K = 6 by default), across the origins with a settled K-week outcome,
#      equals the number observed; psi held at CASCADE_PSI. Pooled deviance and AUC-PR skill are
#      recorded along the search path; the returned bracket is numerical width, not uncertainty.
#      This LEVEL CALIBRATION IS IN-SAMPLE: delta is fitted to exactly the outcomes it is later
#      re-scored against. If no delta can be fitted the step-1 factor is returned unchanged, which
#      run_cascade.R treats as a FAILED calibration and stops on.
#
#   5. Uncertainty. The count-match solver produces no standard error, so `delta_info$se_log`
#      stays the step-1 factor's zone-clustered bootstrap SE (an approximation, stated where it is
#      set). run_cascade.R passes it as delta_sd_log, and simulate_cascade() draws delta per
#      posterior-parameter group on it for the published scenario products.
#
#   6. Residual check, NOT validation (cascade_calibration_report()). The origins are re-simulated
#      at the fitted delta, at an OFFSET seed. The pooled COUNT ratio is ~1 by construction (step
#      4); only the rank/discrimination columns (AUC-PR skill, rank of truth) are something delta
#      was not fitted to. cascade_backtest() in run_cascade.R reuses these origins and this K, so
#      its count ratio is in-sample too (`level_in_sample`). Nothing in this file is an
#      out-of-sample check of the cascade's LEVEL.
#
# HISTORY — why it is built this way. Until 2026-09-11 delta was 1 / calibration-in-the-large of
# the 1-week Bayesian model by the MOMENT estimator, clamped to [0.2, 1.0] (the worse of 16b's two
# estimators on this frame, and a clamp that could not correct an UNDER-running cascade), and psi
# was fitted so the week-8 cumulative count matched "mean weekly new invasions over the last 6
# weeks x 8" — a constant-pace extrapolation, biased low by reporting truncation. In the
# 2026-09-07 run that psi fit was CENSORED at 0 with converged = FALSE and the layer was published
# anyway. At the time delta was not fitted to the backtest's K-week counts, and cascade_backtest()
# reported, at K = 6 over three origins:
#
#     origin        predicted new   observed new   ratio   AUC-PR skill
#     2026-06-23        26.7            14         1.90        17.6x
#     2026-06-30        23.8            14         1.70        17.5x
#     2026-07-07        24.3            13         1.87        17.1x
#
# i.e. good ranking and a level ~1.8x too high. A two-stage design followed (delta at h = 1 from
# 16b, psi fitted to the realised K-week count); psi then proved unidentifiable (censored at a
# bound, or far from its own root, while deviance and AUC-PR kept improving to the search
# ceiling), so it was fixed at 0 and delta alone is now fitted at K (step 4). That buys a
# multi-week level at the price of the week-1 tie to the short-horizon model, which step 3 reports.
#
# LEAKAGE. Every origin is reconstructed as-of its own cutoff: training counts are what was
# knowable at cutoff + 6, through the LFO's own dispatcher (asof_counts_source(): the archived
# snapshot of that date when one exists, else the line list re-censored to that date), nowcast
# with the matching truncation; the hazard is refit on those counts only, with the validated
# model's rolling as-of design (cascade_invasion_design()); R_eff is refit at that issue date;
# and the OUTCOME is read from the final record. If that as-of reconstruction fails (or no line
# list is passed) the origin falls back to the FINAL-COUNT SLICE, which carries cases reported
# after the origin; this is warned and recorded per origin as `training_source = "final_slice"`.
# Outcomes are used only once SETTLED (LFO_RELIABLE_EVAL_AGE_DAYS), because delta is fitted to
# them. The evaluation set is the intersection of "at risk as known at the cutoff" and "at risk
# in the final record", so a zone already invaded but not yet reported is neither scored as a
# false alarm nor counted as an event.
#
# Sourced after 30/31/32/33 and after 16b_invasion_recalibration.R (fit_invasion_delta,
# bootstrap_invasion_delta, RECAL_BAND).
# =============================================================================

# (`%||%` is defined once, in 00_config.R. A stray `if (!exists("%||%"))` left here by the
# 2026-09-28 consolidation swallowed the next expression -- the definition of
# cascade_calibration_origins() -- whenever %||% already existed, i.e. always.)

# ---------------------------------------------------------------------------
# The hazard design, built the way the VALIDATED short-horizon model builds its own
# ---------------------------------------------------------------------------

#' The as-of ("asof" regime) truncation the validated model's rolling predictors and per-fold
#' nowcasts use (run_all.R's `.trunc_asof`): onset -> OBSERVATION, the delay that truncates a
#' re-censored origin. Disk-cached by epinow2_truncation(), so resolving it again is cheap.
#' @return a delay spec, or NULL (with a warning) when it cannot be resolved.
cascade_asof_truncation <- function(linelist) {
  if (is.null(linelist) || !exists("epinow2_truncation", mode = "function") ||
      !exists(".trunc_as_delay_spec", mode = "function")) return(NULL)
  tryCatch(.trunc_as_delay_spec(epinow2_truncation(regime = "asof", ll = linelist), "asof"),
           error = function(e) {
             warning(sprintf(paste0("[cascade] as-of truncation unavailable (%s): the cascade's ",
                                    "re-censored origins fall back to the onset->sample delay and ",
                                    "its hazard designs lose their rolling as-of predictors."),
                             conditionMessage(e)), call. = FALSE)
             NULL })
}

#' Build the cascade's invasion-hazard design EXACTLY as the validated model builds its own.
#'
#' WHY THIS EXISTS. The validated short-horizon model (make_bayes_lfo_model() and the deployed
#' fit_bayes_suite(), 21_bayesian_renewal.R) builds every design with ROLLING AS-OF PREDICTORS:
#' each training transition's import force comes from the counts that were known at that
#' transition, nowcast with the regime's truncation, instead of from one fully-revised matrix.
#' The cascade used to refit its hazard on the shared matrix, which puts beta0 on a different
#' scale (calibrated against complete-week lambdas while the forecast's lambda carries the
#' nowcast inflation — the ~2x artefact build_invasion_design() documents). Any delta carried
#' between the two models — 16b's short-horizon factor, the import term's delta_imp — was then
#' a factor for a different intercept. Same gate as the validated model:
#' ROLLING_PREDICTORS on, a line list, and a truncation.
#'
#' @param trunc_delay the as-of truncation (cascade_asof_truncation()).
#' @return build_invasion_design()'s list, plus `rolling_predictors` (logical).
cascade_invasion_design <- function(zone_week_nc, layer, kernel = CASCADE_KERNEL,
                                    gt = CASCADE_GT, linelist = layer$ll, trunc_delay = NULL) {
  .want <- isTRUE(get0("ROLLING_PREDICTORS", ifnotfound = TRUE))
  .roll_on <- .want && !is.null(linelist) && !is.null(trunc_delay)
  if (.want && !.roll_on)
    warning(sprintf(paste0("[cascade] ROLLING_PREDICTORS is on but %s is missing, so the %s hazard ",
                           "is refit on the shared count matrix — NOT the design the validated ",
                           "short-horizon model uses. Its intercept, and any delta transferred to ",
                           "it, are on a different scale."),
                    if (is.null(linelist)) "the line list" else "the as-of truncation", kernel),
            call. = FALSE)
  des <- build_invasion_design(zone_week_nc, layer$mobility_matrices, layer$gt_pmfs,
           layer$covariates, layer$osrm_mat, layer$zones_all, mob = kernel, gt = gt,
           rolling_ll    = if (.roll_on) linelist else NULL,
           rolling_delay = if (.roll_on) trunc_delay else NULL)
  des$rolling_predictors <- .roll_on
  des
}

#' Fit the cascade's seeding hazard on a cascade design: fixed beta (brms), or the featured
#' model's time-varying FORM (CASCADE_FORM, 30_projection_config.R; the Stan engine of
#' 21c_frailty_models.R) with the same covariates. The ONE entry point for every cascade fit --
#' production, the calibration origins, the kernel sensitivity and the overdispersion sweep -- so
#' the hazard that is calibrated is the hazard that is deployed.
#' @param design cascade_invasion_design() output for `zone_week_nc`.
#' @param trunc_delay the as-of truncation the design's rolling predictors used.
#' @return a brmsfit or a stan_invasion_fit, or NULL when the fit fails (a form's failure is
#'   warned, as a brms failure is logged by the brms fitter).
cascade_fit_hazard <- function(design, zone_week_nc, layer, kernel = CASCADE_KERNEL,
                               gt = CASCADE_GT, linelist = layer$ll, trunc_delay = NULL,
                               iter = 2000L, cov_spec = CASCADE_COV_SPEC,
                               form = CASCADE_FORM) {
  if (identical(form, "fixed"))
    return(fit_bayes_renewal(design, cov_spec = cov_spec, iter = iter, chains = 2L))
  sw <- Filter(function(x) identical(x$form, form), get0("TV_FORMS", ifnotfound = list()))
  if (length(sw) != 1L) stop(sprintf("[cascade] unknown hazard form '%s'.", form), call. = FALSE)
  # The forms are defined on the rolling as-of design (their exposure IS the as-of import force).
  if (!isTRUE(design$rolling_predictors) || is.null(linelist) || is.null(trunc_delay))
    stop(sprintf(paste0("[cascade] the featured form '%s' needs the rolling as-of design (a line ",
                        "list and the as-of truncation); refusing to fit it on another design."),
                 form), call. = FALSE)
  # Never below the forms' validated back-test effort (TV_FORM_ITER): the walk adds a parameter
  # per week, and the calibration origins pass the lighter brms back-test effort.
  iter <- max(as.integer(iter), as.integer(get0("TV_FORM_ITER", ifnotfound = 1000L)))
  tryCatch({
    feat <- intersect(cov_spec, design$feat)
    ex <- invasion_exposure(zone_week_nc, linelist, layer$zones_all, layer$mobility_matrices,
                            layer$gt_pmfs, layer$covariates, layer$osrm_mat, kernel, gt, feat,
                            design$center, design$scale, trunc_delay)
    fit <- stan_invasion_fit(design, ex, layer$zones_all, cov_spec, sw[[1]]$use_rw, sw[[1]]$use_frailty,
                             iter = iter, label = sprintf("cascade %s %s", kernel, form))
    # Convergence is checked, not assumed: these fits are not otherwise inspected.
    if (!is.null(fit) && ((fit$diag$rhat_max > 1.05) %in% TRUE || (fit$diag$n_divergent > 0) %in% TRUE))
      warning(sprintf(paste0("[cascade] the %s hazard fit on %s did not converge cleanly: R-hat max ",
                             "%.3f, %d divergent transition(s)."), form, kernel,
                      fit$diag$rhat_max, as.integer(fit$diag$n_divergent)), call. = FALSE)
    fit
  }, error = function(e) {
    warning(sprintf("[cascade] the %s hazard fit on %s failed: %s", form, kernel,
                    conditionMessage(e)), call. = FALSE)
    NULL })
}

#' delta_info for the SHORT-HORIZON factor (16b), with every attribute it carries.
#'
#' Shared by cascade_calibrate() and run_cascade.R's deliberately-uncalibrated path, so the two
#' cannot drift. The fallback there used to keep only `estimator`, which silently set the
#' projection's delta_sd_log to 0 although the factor carries a zone-clustered bootstrap SE.
cascade_short_horizon_delta_info <- function(delta) {
  list(delta = as.numeric(delta), estimator = attr(delta, "estimator") %||% NA_character_,
       lo = attr(delta, "lo") %||% NA_real_, hi = attr(delta, "hi") %||% NA_real_,
       se_log = attr(delta, "se_log") %||% NA_real_,
       cal_in_large = attr(delta, "cal_in_large") %||% NA_real_,
       boundary_hit = attr(delta, "boundary_hit") %||% NA,
       source = attr(delta, "source") %||% NA_character_)
}

# ---------------------------------------------------------------------------
# Held-out origin construction (shared by calibration AND the backtest)
# ---------------------------------------------------------------------------

#' Build leakage-free held-out cascade origins.
#'
#' One entry per usable origin, each carrying everything needed to simulate the
#' cascade from that origin and score it against what actually happened:
#' an as-of `prep`, the evaluation zone set, and the realised outcomes at each
#' horizon in `truth_horizons`.
#'
#' @param layer the cascade layer (zone_week_nc, zones_all, mobility_matrices, ...).
#' @param zone_province named province vector (for the province-pooled R_eff).
#' @param cutoffs_from_end origins expressed as "this many weeks before the last".
#' @param truth_horizons horizons (weeks) at which realised outcomes are recorded.
#' @param min_eval_age_days a truth window is recorded only once it has been closed for
#'   this many days. Defaults to LFO_RELIABLE_EVAL_AGE_DAYS (00_config.R), the age at which
#'   the suite treats an invasion count as SETTLED and the one its own recalibration factor is
#'   trained on: delta is FITTED to these counts, so they must be settled, not merely admitted
#'   (it was 10 d, under the p90 onset->confirmation delay, so the most recent windows were
#'   still filling and pulled delta down). Settled is not "complete": a late confirmation can
#'   still arrive. Applied PER HORIZON: an origin whose 1-week window is settled but whose
#'   K-week window is not still contributes its 1-week outcome (`complete_horizons` records
#'   which were kept).
#' @param linelist the line list, for the as-of reconstruction. Without it the
#'   origins train on the FINAL revised counts (a revision leak) and the function
#'   says so loudly.
#' @return list of origins; each list(cutoff_date, Ccol, prep, eval_zones, eval_idx,
#'   truth = named list of 0/1 vectors over eval_zones, n_atrisk, n_events).
cascade_calibration_origins <- function(layer, zone_province = NULL,
                                        cutoffs_from_end = c(12L, 11L, 10L, 9L, 8L),
                                        truth_horizons = c(1L, 6L),
                                        min_eval_age_days = get0("LFO_RELIABLE_EVAL_AGE_DAYS",
                                                                 ifnotfound = 25L),
                                        min_weeks_from_start = get0("CASCADE_CALIB_MIN_WEEKS",
                                                                   ifnotfound = 7L),
                                        # The deployed hazard scale the import term in each
                                        # origin's R denominator is built on. The default of 1
                                        # is the identity (raw scale) for a direct call;
                                        # cascade_calibrate() always passes the fitted value,
                                        # which is what keeps calibration and production on
                                        # one scale.
                                        delta_imp = 1,
                                        kernel = CASCADE_KERNEL, gt = CASCADE_GT,
                                        iter_fit = 600L, linelist = layer$ll,
                                        analysis_date = get0("ANALYSIS_DATE", ifnotfound = NA)) {
  zones_all <- layer$zones_all
  stopifnot(length(truth_horizons) >= 1L, all(truth_horizons >= 1L))
  Hmax <- max(as.integer(truth_horizons))

  # Outcome clock: first week with a confirmed case in the FINAL record. The
  # nowcast is multiplicative, so confirmed_nc > 0 has exactly the same support as
  # confirmed > 0 and the first-case week is identical either way.
  Yfull    <- .count_wide(layer$zone_week_nc, zones_all, "confirmed_nc")
  first_wk <- apply(Yfull, 1, function(x) { w <- which(x > 0); if (length(w)) min(w) else NA_integer_ })
  nT       <- ncol(Yfull)
  weeks_sorted <- sort(unique(layer$zone_week_nc$week_start))
  .asof <- suppressWarnings(as.Date(analysis_date))

  if (is.null(linelist))
    warning("[cascade_calib] no `linelist`: each origin trains on the FINAL revised ",
            "zone-week counts, which contain cases reported only AFTER that origin. That is ",
            "the training-side revision leak reaggregate_asof() exists to close, and the ",
            "calibration fitted on it will be optimistic. Pass linelist = layer$ll.",
            call. = FALSE)

  # TRUNCATION REGIME: "asof" (2026-09-21). Each origin below is rebuilt with
  # reaggregate_asof(), which censors on linelist_observation_date() — so what truncates an
  # origin's most recent week is onset -> OBSERVATION, not onset -> appearance-in-extract.
  # That is a DIFFERENT distribution from the one run_all.R's deployed nowcast uses, and
  # epinow2_truncation() fits it separately; passing the deployed ("extract") fit here would
  # over-correct every calibration origin and bias `delta` low. Resolved ONCE outside the
  # loop: the fit is disk-cached, but re-resolving per origin would still re-read and re-hash
  # it. NULL falls through to apply_nowcast_correction()'s own resolver, which warns.
  # It is ALSO the rolling-predictor delay for every origin's hazard design
  # (cascade_invasion_design()), exactly as run_all.R passes its .trunc_asof to the LFO.
  .trunc_calib <- cascade_asof_truncation(linelist)
  if (!is.null(.trunc_calib))
    message(sprintf("[cascade_calib] as-of truncation: %s (mean %.2f d, sd %.2f d, %s)",
                    .trunc_calib$family, .trunc_calib$mean, .trunc_calib$sd, .trunc_calib$estimator))

  origins <- list(); n_burnin <- 0L
  for (ce in sort(unique(as.integer(cutoffs_from_end)), decreasing = TRUE)) {
    Ccol <- nT - ce
    if (Ccol < 4L) next
    cutoff_date <- weeks_sorted[Ccol]
    # ASCERTAINMENT BURN-IN. In the first weeks after detection the apparent growth is
    # dominated by the response arriving and a backlog appearing at once, which no renewal
    # estimator can separate from transmission without an explicit reporting model. At a
    # week-4 origin on this outbreak every zone with cases returned R at its ceiling and the
    # national estimate likewise, so an origin inside the ramp calibrates the cascade against
    # a reporting artefact. Excluded here rather than censored in the estimator, because the
    # estimate is not wrong about the DATA — the data are not yet about transmission.
    if (Ccol < as.integer(min_weeks_from_start)) {
      n_burnin <- n_burnin + 1L
      message(sprintf(paste0("[cascade_calib] origin %s skipped: week %d is inside the %d-week ",
                             "ascertainment burn-in."), format(cutoff_date), Ccol,
                      as.integer(min_weeks_from_start)))
      next
    }
    # PER-HORIZON completeness, not one gate on the longest window. A truth window is
    # usable only once it has both elapsed within the record AND been closed long enough
    # to be SETTLED (min_eval_age_days). Gating on the LONGEST horizon alone threw away origins
    # whose 1-week outcome was perfectly usable, which is what starved the delta transfer
    # check of events (5 across two origins, below the 10-event floor). This mirrors the
    # per-horizon guard run_invasion_lfo() applies to its folds.
    hs_ok <- Filter(function(H) {
      if (Ccol + H > nT) return(FALSE)
      if (is.na(.asof)) return(TRUE)
      as.numeric(.asof - (cutoff_date + H * 7 + 7)) >= min_eval_age_days
    }, as.integer(truth_horizons))
    if (!length(hs_ok)) {
      message(sprintf(paste0("[cascade_calib] origin %s skipped: no truth window is both within ",
                             "the record and settled (min_eval_age_days)."), format(cutoff_date)))
      next
    }

    # ---- as-of training frame (leakage-free) --------------------------------
    # ORIGIN = cutoff_date + 6, the LAST DAY OF THE LAST TRAINING WEEK — the same convention
    # as the LFO (16_invasion_eval.R), the fold diagnostics (20_forecast_detail.R), the R(t)
    # issue date (21_bayesian_renewal.R) and deployment itself, where the weekly grid is
    # anchored so the final week ENDS on ANALYSIS_DATE and the as-of date IS that day
    # (00_config.R WEEK_ANCHOR).
    #
    # It was cutoff_date + 7, which handed every calibration origin ONE EXTRA DAY of reporting
    # for its most recent week. Under the retired onset->sample delay that week was multiplied
    # by 2.631x at +6 but only 2.185x at +7 — a 20.4% weaker import force; the fitted as-of
    # truncation is longer, so the gap between the two conventions is WIDER, not narrower. `delta` is FITTED on these
    # origins and then APPLIED (run_cascade.R) to a deployed layer built at 2.631x, so delta
    # was absorbing a regime the deployment never presents: the published 13-week invasion
    # projection over-projected. The identical +7 -> +6 correction was made in the LFO on
    # 2026-09-19 and never propagated here.
    #
    # THE SAME DISPATCHER THE LFO USES (asof_counts_source(), 22_daily_reissue.R), so an origin's
    # training counts are what was knowable at it in exactly the sense run_invasion_lfo()
    # means: the archived snapshot of that date, nowcast with the truncation estimable on that
    # date, when the archive supports it (LFO_TRAINING_SOURCE = "vintage"); otherwise the final
    # line list re-censored to cutoff + 6 and nowcast with the shared as-of truncation. This
    # block used to re-censor unconditionally, so on a vintage run the cascade calibrated on a
    # different information state from the folds 16b's short-horizon factor was fitted on.
    zw_c <- NULL; .src_c <- list(source = "final_slice", vintage_date = as.Date(NA),
                                 why = "no line list")
    if (!is.null(linelist) && exists("reaggregate_asof") && exists("apply_nowcast_correction")) {
      zw_c <- tryCatch({
        .wk_c <- weeks_sorted[weeks_sorted <= cutoff_date]
        src <- if (exists("asof_counts_source", mode = "function"))
                 asof_counts_source(linelist, zones_all, cutoff_date + 6, week_spine = .wk_c,
                                    delay_recensor = .trunc_calib)
               else list(raw = reaggregate_asof(linelist, zones_all, cutoff_date + 6,
                                                 week_spine = .wk_c),
                         delay = .trunc_calib, nowcast_date = cutoff_date + 6,
                         source = "recensor", vintage_date = as.Date(NA),
                         why = "22_daily_reissue.R not loaded")
        raw <- src$raw[src$raw$week_start <= cutoff_date, , drop = FALSE]
        .src_c <- src[c("source", "vintage_date", "why")]
        apply_nowcast_correction(raw, analysis_date = src$nowcast_date, delay = src$delay)
      }, error = function(e) {
        warning(sprintf(paste0("[cascade_calib] origin %s: as-of reconstruction failed (%s); this ",
                               "origin falls back to the revision-leaky final-count slice."),
                        format(cutoff_date), conditionMessage(e)), call. = FALSE)
        NULL })
    }
    if (is.null(zw_c)) {
      zw_c <- layer$zone_week_nc[layer$zone_week_nc$week_start <= cutoff_date, , drop = FALSE]
      .src_c <- list(source = "final_slice", vintage_date = as.Date(NA),
                     why = if (is.null(linelist)) "no line list" else "as-of reconstruction failed")
    }

    # ---- hazard + source dynamics refit on that frame ONLY ------------------
    # With the validated model's design (rolling as-of predictors), so the origin's intercept
    # is on the scale delta is fitted for and deployed on (cascade_invasion_design()).
    dz <- tryCatch(cascade_invasion_design(zw_c, layer, kernel = kernel, gt = gt,
                                           linelist = linelist, trunc_delay = .trunc_calib),
             error = function(e) NULL)
    if (is.null(dz) || dz$n_events < 3L) {
      message(sprintf("[cascade_calib] origin %s skipped: design unavailable or < 3 training events.",
                      format(cutoff_date)))
      next
    }
    fz <- cascade_fit_hazard(dz, zw_c, layer, kernel = kernel, gt = gt, linelist = linelist,
                             trunc_delay = .trunc_calib, iter = iter_fit)
    if (is.null(fz)) {
      message(sprintf("[cascade_calib] origin %s skipped: hazard fit failed.", format(cutoff_date)))
      next
    }
    # Built through the SAME entry point as the production run, on the same delta, so the
    # fitted (delta, psi) describe the simulator that is actually deployed.
    # issue_date is THIS origin's own as-of moment, the same cutoff + 6 the count
    # reconstruction above uses. It is what stops the national R anchor from being fitted on
    # cases reported after the origin: without it the held-out evaluation would score a
    # simulator that had already seen its own future, and the calibration would be optimistic
    # in a way no output reports. cascade_reff() has no default for it precisely so that this
    # line cannot be omitted by accident.
    rz <- cascade_reff(zw_c, zones_all, layer, fz, dz, delta = delta_imp,
                       zone_province = zone_province, gt = gt, kernel = kernel,
                       linelist = linelist, issue_date = cutoff_date + 6)
    layer_c <- layer; layer_c$zone_week_nc <- zw_c
    prep_c  <- cascade_prepare(layer_c, fz, dz, rz, kernel = kernel, gt = gt)

    # ---- evaluation set and realised outcomes ------------------------------
    # At-risk under BOTH definitions: what the model knew at the cutoff (prep's
    # affected0, from the as-of counts) AND the final record. A zone already
    # invaded at the cutoff but not yet reported is therefore neither scored as a
    # false alarm nor counted as an event — scoring it either way would attribute a
    # reporting lag to the model. (run_invasion_lfo() applies the same intersection
    # by filtering each model's rows to the fold's final-record at-risk set.)
    atrisk_final <- !(rowSums(Yfull[, seq_len(Ccol), drop = FALSE] > 0) > 0)
    eval_mask    <- atrisk_final & !prep_c$affected0
    eval_idx     <- which(eval_mask)
    if (!length(eval_idx)) next
    fw <- first_wk; fw[is.na(fw)] <- .Machine$integer.max
    truth <- lapply(hs_ok, function(H)
      as.integer(fw[eval_idx] > Ccol & fw[eval_idx] <= Ccol + H))
    names(truth) <- as.character(hs_ok)

    origins[[length(origins) + 1L]] <- list(
      # Carried on the origin so every simulation driven from it — the delta h=1 refit, the
      # K-week delta search, the residual report — runs the SAME model as production. Indexed by
      # zone name inside simulate_cascade(), so the full national vector is correct for every
      # origin.
      # NOT "with within-zone depletion engaged": pop_vec only ARMS depletion. Whether it
      # engaged within-zone depletion, which was REMOVED (2026-09-19) and was never enabled by
      # any caller in this suite, so s_frac is NULL and pop_vec is inert in every simulation
      # here. Calibration and production are still consistent (both have it off) — which is
      # what this comment is really asserting — but the plumbing is dormant, not armed.
      pop_vec = layer$pop,
      # R_nat at this origin, kept so the burn-in rule is auditable rather than asserted:
      # measured on the 14-week frame, the national estimate ran 10.6 / 8.5 / 4.3 / 2.9 at
      # weeks 4 / 5 / 6 / 7 and only reached the published EVD range (1.5-2.5) from week 8.
      R_nat = rz$R_nat, rt_source = rz$rt_source %||% NA_character_, n_with_cases = rz$n_with_cases,
      n_prior_dominated = rz$n_prior_dominated, import_share = rz$import_share,
      cutoff_date = cutoff_date, Ccol = Ccol, prep = prep_c,
      # Where this origin's training counts came from (asof_counts_source()), and whether its
      # hazard was refit on the validated model's rolling design — auditable per origin.
      training_source = .src_c$source, vintage_date = as.Date(.src_c$vintage_date),
      rolling_predictors = isTRUE(dz$rolling_predictors),
      eval_zones = zones_all[eval_idx], eval_idx = eval_idx,
      truth = truth, complete_horizons = hs_ok, n_atrisk = length(eval_idx),
      n_events = vapply(truth, sum, integer(1)),
      n_asof_only = sum(!prep_c$affected0 & !atrisk_final))
  }

  if (!length(origins)) {
    warning("[cascade_calib] no usable held-out origins; the cascade cannot be calibrated ",
            "on this snapshot.", call. = FALSE)
    return(list())
  }
  .ev <- function(o, H) { k <- as.character(H)
                          if (k %in% names(o$truth)) as.integer(sum(o$truth[[k]])) else NA_integer_ }
  message(sprintf("[cascade_calib] %d held-out origin(s): %s | at-risk %s | events@H=%d %s | events@H=1 %s",
                  length(origins),
                  paste(vapply(origins, function(o) format(o$cutoff_date), ""), collapse = ", "),
                  paste(vapply(origins, function(o) o$n_atrisk, 0L), collapse = "/"),
                  Hmax,
                  paste(vapply(origins, .ev, 0L, H = Hmax), collapse = "/"),
                  paste(vapply(origins, .ev, 0L, H = 1L), collapse = "/")))
  if (length(origins))
    message(sprintf("[cascade_calib] R_nat by origin: %s (published EVD R sits near 1.5-2.5)",
                    paste(vapply(origins, function(o)
                      sprintf("%s %.2f", format(o$cutoff_date), o$R_nat), ""), collapse = ", ")))
  # A SPECIFIC diagnostic when the burn-in and the requested origin grid do not overlap.
  # Falling through to the generic "no held-out origins" would read as "the data are too
  # short", when in fact the operator asked for origins that all sit inside the
  # ascertainment ramp and the grid is what needs moving.
  if (!length(origins) && n_burnin > 0L)
    warning(sprintf(paste0("[cascade_calib] ALL %d requested origin(s) fall inside the %d-week ",
                           "ascertainment burn-in, so nothing is calibrated. The record has %d ",
                           "weeks and cutoffs_from_end = c(%s) puts every cutoff at week <= %d. ",
                           "Move the grid later (cutoffs_from_end no greater than %d), or lower ",
                           "CASCADE_CALIB_MIN_WEEKS if you have reason to trust the early ",
                           "origins — on the 2026 frame the national R was 10.6 at week 4 and ",
                           "only became plausible from week 8."),
                   n_burnin, as.integer(min_weeks_from_start), nT,
                   paste(sort(unique(as.integer(cutoffs_from_end))), collapse = ","),
                   nT - min(as.integer(cutoffs_from_end)),
                   nT - as.integer(min_weeks_from_start)), call. = FALSE)
  origins
}

# ---------------------------------------------------------------------------
# Simulate the held-out origins and pool (p, y)
# ---------------------------------------------------------------------------

#' Pooled (predicted probability, realised outcome) pairs over held-out origins.
#'
#' @return data.frame(origin, health_zone, p, y) restricted to each origin's
#'   evaluation set; `p` is the reach probability by week `H`.
#' SEEDING NOTE. Every origin is simulated from the SAME `seed`. That is deliberate: it
#' holds the posterior-parameter and R_eff draws fixed across successive evaluations of
#' (delta, psi), which is what makes the psi search curve smooth enough to bracket and
#' bisect. It does correlate the process noise BETWEEN origins within one evaluation, so
#' the pooled count is unbiased but its Monte-Carlo variance is not the independent-origin
#' variance; do not read the spread across origins as independent replication.
.cascade_pool_po <- function(origins, H, delta, psi, n_mc, seed, pop_vec = NULL) {
  scP <- CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]
  parts <- lapply(origins, function(o) {
    # pop_vec is threaded for symmetry with production; within-zone depletion is modelled in
    # neither (removed 2026-09-19), so it is inert, and psi is fixed at 0 (not fitted).
    sim <- simulate_cascade(o$prep, scP, n_mc = n_mc, horizon = as.integer(H),
                            delta = delta, psi = psi, seed = seed,
                            pop_vec = pop_vec %||% o$pop_vec)
    rt  <- cascade_reach_table(sim, horizons = as.integer(H))
    p   <- rt$p_invasion[match(o$eval_zones, rt$health_zone)]
    data.frame(origin = format(o$cutoff_date), health_zone = o$eval_zones,
               p = p, y = o$truth[[as.character(as.integer(H))]],
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, parts)
  d[is.finite(d$p) & !is.na(d$y), , drop = FALSE]
}

#' Pooled Bernoulli deviance (lower is better) of a set of (p, y) pairs.
.cascade_deviance <- function(p, y) {
  p <- pmin(pmax(p, 1e-9), 1 - 1e-9)
  -2 * sum(y * log(p) + (1 - y) * log(1 - p))
}

# ---------------------------------------------------------------------------
# Step 3 (header) — week-1 transfer check: refit delta from the cascade's OWN week-1 hazard
# ---------------------------------------------------------------------------

#' Maximum-likelihood hazard-scale factor for the cascade at h = 1.
#'
#' Simulates one week from each held-out origin at delta = 1 and fits the factor
#' that calibrates those week-1 probabilities to the realised week-1 invasions,
#' with the SAME estimator, band and zone-clustered bootstrap 16b uses for the
#' short-horizon model (mu -> delta*mu is exactly p -> 1-(1-p)^delta). psi is
#' irrelevant here and passed as 0: at the first simulated week no zone has been
#' seeded yet, so f_inv == f_inv0 and sat == 1 for every psi.
#'
#' @return list(delta, lo, hi, se_log, n_rows, n_events, boundary_hit, pooled)
cascade_delta_h1_refit <- function(origins, n_mc = 600L, seed = CASCADE_SEED,
                                   band = get0("RECAL_BAND", ifnotfound = c(1e-4, 1e4)),
                                   # LOCAL reporting floor, not a recalibration constant:
                                   # 16b's fit is penalised and needs no event floor. This one
                                   # only decides whether a cascade-transfer refit carries
                                   # enough evidence to be worth reporting at all.
                                   min_events = 10L) {
  if (!length(origins)) return(NULL)
  if (!exists("fit_invasion_delta"))
    stop("[cascade_calib] fit_invasion_delta() not found — source 16b_invasion_recalibration.R.")
  origins <- Filter(function(o) "1" %in% names(o$truth), origins)
  if (!length(origins)) {
    warning("[cascade_calib] no held-out origin carries a complete 1-week outcome; the delta ",
            "transfer check is skipped.", call. = FALSE)
    return(NULL)
  }
  d <- .cascade_pool_po(origins, H = 1L, delta = 1, psi = 0, n_mc = n_mc, seed = seed)
  n_ev <- sum(d$y == 1L)
  if (!nrow(d) || n_ev < min_events) {
    warning(sprintf(paste0("[cascade_calib] week-1 refit not attempted: %d event(s) across the ",
                           "held-out origins, below the %d-event floor."), n_ev, min_events),
            call. = FALSE)
    return(list(delta = NA_real_, lo = NA_real_, hi = NA_real_, se_log = NA_real_,
                n_rows = nrow(d), n_events = n_ev, boundary_hit = NA, pooled = d))
  }
  dh <- fit_invasion_delta(d$p, d$y, band = band)
  ci <- if (is.finite(dh) && exists("bootstrap_invasion_delta"))
          bootstrap_invasion_delta(d$p, d$y, d$health_zone, band = band)
        else list(lo = NA_real_, hi = NA_real_, se_log = NA_real_, n_boot_used = 0L)
  list(delta = dh, lo = ci$lo, hi = ci$hi, se_log = ci$se_log,
       n_rows = nrow(d), n_events = n_ev,
       boundary_hit = is.finite(dh) &&
         (abs(log(dh) - log(band[1])) < 1e-4 || abs(log(dh) - log(band[2])) < 1e-4),
       pooled = d)
}

# ---------------------------------------------------------------------------
# RETIRED (sensitivity sweeps only) — psi against REALISED multi-week outcomes
# ---------------------------------------------------------------------------

#' Fit the frontier-saturation parameter to the realised K-week invasion counts at the
#' held-out origins (IN-SAMPLE for those counts, exactly as for delta; retired, kept for sweeps).
#'
#' The modelled quantity is the cascade's expected number of new invasions among
#' the held-out evaluation zones over the K-week window; the target is the number
#' that ACTUALLY occurred in that window. sat = exp(-psi * max(f_inv - f_inv0, 0))
#' is non-increasing in psi, so the modelled count is monotone non-increasing and
#' the root is unique — the same property the previous solver relied on, now with
#' an observed rather than an extrapolated target.
#'
#' STRUCTURAL FLOOR, and it matters for reading a censored upper bound. The damping
#' acts only on the INCREMENT of the inflow-weighted invaded fraction over its value
#' at the origin (max(f_inv - f_inv0, 0)), so zones seeded directly by the front that
#' already existed are not damped at all, however large psi is. psi therefore cannot
#' drive the modelled count below that floor, and it cannot correct a level error
#' that is already present in week 1 — only delta can. A fit censored at the UPPER
#' bound is evidence about delta or the hazard, not about saturation.
#'
#' Monte-Carlo noise: every evaluation uses the SAME seed, so the posterior
#' parameter draws and the seeding R_eff draws are identical across psi values and
#' the residual noise is process-only. It is not eliminated (the random streams
#' desynchronise once the seeded sets differ), so `tol` is a relative tolerance on
#' the count, not an exact root.
#'
#' @return list(psi, converged, boundary_hit, boundary_side, curve, target, ...)
cascade_fit_psi_oos <- function(origins, delta, K = 6L,
                                psi_lo = 0, psi_hi = CASCADE_PSI_MAX,
                                n_coarse = 6L, n_mc = 400L, seed = CASCADE_SEED,
                                tol = CASCADE_PSI_TOL, max_refine = 12L,
                                psi_tol = get0("CASCADE_PSI_BRACKET_TOL",
                                               ifnotfound = 0.02)) {
  if (!length(origins)) return(NULL)
  K <- as.integer(K)
  key <- as.character(K)
  # Per-horizon completeness means a recent origin can carry its 1-week outcome but not
  # its K-week one. SELECT the origins that do rather than aborting: refusing to fit
  # because one origin is too recent would throw away the whole calibration.
  has_K  <- vapply(origins, function(o) key %in% names(o$truth), logical(1))
  if (!any(has_K)) {
    warning(sprintf(paste0("[cascade_calib] no held-out origin carries a settled ",
                           "%d-week outcome; psi cannot be fitted and is left at its default."), K),
            call. = FALSE)
    return(NULL)
  }
  if (!all(has_K))
    message(sprintf("[cascade_calib] psi fitted on %d of %d origin(s) — %d lack a complete %d-week outcome.",
                    sum(has_K), length(origins), sum(!has_K), K))
  origins <- origins[has_K]

  target <- sum(vapply(origins, function(o) sum(o$truth[[key]]), numeric(1)))
  if (target <= 0) {
    warning("[cascade_calib] zero realised invasions across the held-out windows; psi is not ",
            "identifiable and is left at its configured default.", call. = FALSE)
    return(NULL)
  }
  trace <- list()
  f_cum <- function(ps) {
    d <- .cascade_pool_po(origins, H = K, delta = delta, psi = ps, n_mc = n_mc, seed = seed)
    v <- sum(d$p)
    trace[[length(trace) + 1L]] <<- data.frame(
      psi = ps, modelled_cum = v, observed_cum = target,
      deviance = .cascade_deviance(d$p, d$y),
      auc_pr_skill = { ap <- .casc_auc_pr(d$p, d$y); br <- mean(d$y)
                       if (is.finite(ap) && br > 0) ap / br else NA_real_ })
    v
  }

  ladder <- if (n_coarse > 1L)
    exp(seq(log(0.25), log(max(psi_hi, 0.5)), length.out = n_coarse - 1L)) else numeric(0)
  grid  <- sort(unique(c(psi_lo, ladder, psi_hi)))
  curve <- vapply(grid, f_cum, numeric(1))
  gap   <- curve - target                      # > 0 = cascade seeds more than observed
  atol  <- tol * target

  boundary_hit <- FALSE; boundary_side <- NA_character_; converged <- TRUE; psi_hat <- NA_real_
  psi_bracket_lo <- NA_real_; psi_bracket_hi <- NA_real_
  if (all(gap <= 0)) {
    psi_hat <- grid[1L]
    if (abs(gap[1L]) > atol) {
      converged <- FALSE; boundary_hit <- TRUE; boundary_side <- "lower"
      warning(sprintf(paste0("[cascade_fit_psi_oos] at psi=%.2f (no saturation) the cascade already ",
                             "projects only %.1f new zones over the held-out %d-week window(s) vs ",
                             "%.0f observed (%.0f%% under). Saturation cannot correct an UNDER-running ",
                             "front: the miss is in delta, the hazard or R_eff."),
                     psi_hat, curve[1L], K, target,
                     100 * (1 - curve[1L] / target)), call. = FALSE)
    }
  } else if (all(gap > 0)) {
    psi_hat <- grid[length(grid)]
    boundary_hit <- TRUE; boundary_side <- "upper"; converged <- FALSE
    warning(sprintf(paste0("[cascade_fit_psi_oos] saturation did not converge: at psi=%.1f (search ",
                           "maximum) the cascade still projects %.1f new zones over the held-out ",
                           "%d-week window(s) vs %.0f observed (%.0f%% over). psi is CENSORED and the ",
                           "projection will over-run; inspect the kernel's dispersion and delta."),
                   psi_hat, curve[length(curve)], K, target,
                   100 * (curve[length(curve)] / target - 1)), call. = FALSE)
  } else {
    i_hi <- max(which(gap > 0))                # last psi still above target
    lo <- grid[i_hi]; hi <- grid[i_hi + 1L]
    g_lo <- gap[i_hi]; g_hi <- gap[i_hi + 1L]
    # STOPPING RULE — on the BRACKET, not on the first acceptable count. Stopping as soon as
    # |gap| <= atol returns the first psi that lands inside the band, which is not the estimator's
    # value. On the 2026-09 frame that rule returned the very FIRST bisection midpoint (7.2046,
    # gap -3.03, inside the 3.65 band) while the zero crossing sat near 5.5 — psi reported 30%
    # above its own root, to six decimals, and then used as a point value in every projection.
    # Refining until the bracket is narrow makes psi a property of the data rather than of where
    # the search happened to stop. Bisection uses only the SIGN of the gap, so Monte-Carlo noise
    # near the root cannot stop the bracket from closing; it only blurs which end is picked, and
    # the bracket is returned so that residual width is visible rather than implied.
    #
    # MEASURED sensitivity to the DRAW SIZE, recorded so the bracket is never mistaken for psi's
    # statistical uncertainty: with this same stopping rule and seed on the 2026-09-07 frame,
    # n_mc = 400 returned psi = 5.12 and n_mc = 600 returned psi = 4.43. The numerical bracket was
    # ~1.4% of psi in both cases, so the Monte-Carlo component dominates it by roughly an order of
    # magnitude. Raising n_mc narrows that component; the report states the distinction explicitly
    # rather than quoting the bracket as though it were psi's precision.
    for (it in seq_len(max_refine)) {
      if ((hi - lo) <= psi_tol * max(abs(hi), 1e-9)) break
      mid   <- 0.5 * (lo + hi)
      g_mid <- f_cum(mid) - target
      if (g_mid > 0) { lo <- mid; g_lo <- g_mid } else { hi <- mid; g_hi <- g_mid }
    }
    psi_hat   <- if (abs(g_lo) <= abs(g_hi)) lo else hi
    converged <- min(abs(g_lo), abs(g_hi)) <= atol
    psi_bracket_lo <- lo; psi_bracket_hi <- hi
    if (!converged)
      warning(sprintf(paste0("[cascade_fit_psi_oos] bisection exhausted %d refinements without ",
                             "reaching the %.0f%% tolerance; residual %.1f zones on a target of %.0f."),
                      max_refine, 100 * tol, min(abs(g_lo), abs(g_hi)), target), call. = FALSE)
  }

  curve_df <- do.call(rbind, trace)
  curve_df <- curve_df[order(curve_df$psi), , drop = FALSE]
  curve_df <- curve_df[!duplicated(curve_df$psi), , drop = FALSE]
  curve_df$gap <- curve_df$modelled_cum - target
  # The count match is the ESTIMATOR (monotone, unique root); the deviance-optimal
  # psi on the same path is reported beside it as an independent check. A large
  # disagreement means the level and the per-zone allocation disagree, which the
  # count match alone cannot see.
  psi_dev <- curve_df$psi[which.min(curve_df$deviance)]
  # AUC-PR optimum alongside the deviance optimum. Two proper diagnostics agreeing on a psi
  # far from the count match is a much stronger signal than either alone.
  psi_auc <- if (any(is.finite(curve_df$auc_pr_skill)))
    curve_df$psi[which.max(curve_df$auc_pr_skill)] else NA_real_
  # BOUNDARY DETECTION FOR THE PROBABILISTIC CRITERIA. `boundary_hit` above tracks only the
  # COUNT-match root; it is FALSE whenever that root lies inside the interval, even when the
  # deviance and AUC-PR are still improving at the search ceiling. On the 2026-09 frame both
  # were monotone to psi = 40 (deviance 355 -> 334, AUC-PR skill 14.8 -> 15.7) while the count
  # match sat at 11.25 — so the ceiling was doing real work and nothing said so. A criterion
  # pinned at the edge is not an optimum: it means the search interval, not the data, chose it.
  .at_edge <- function(x) is.finite(x) && abs(x - psi_hi) < 1e-9
  dev_at_edge <- .at_edge(psi_dev); auc_at_edge <- .at_edge(psi_auc)
  if (dev_at_edge || auc_at_edge)
    warning(sprintf(paste0("[cascade_fit_psi_oos] a PROBABILISTIC criterion is pinned at the ",
                           "search ceiling psi_hi=%.1f (%s%s%s). The count match returned %.2f. ",
                           "This is the documented diagnostic that frontier saturation is ",
                           "absorbing misspecification it cannot represent — most likely the ",
                           "kernel's dispersion or delta — and NOT evidence that psi should be ",
                           "%.1f. Widen psi_hi to test whether the criterion keeps improving, ",
                           "and do not adopt the edge value as a fit."),
                   psi_hi,
                   if (dev_at_edge) "deviance" else "",
                   if (dev_at_edge && auc_at_edge) " and " else "",
                   if (auc_at_edge) "AUC-PR" else "",
                   psi_hat, psi_hi), call. = FALSE)
  .near <- function(x) which.min(abs(curve_df$psi - x))
  dev_at_psi <- curve_df$deviance[.near(psi_hat)]
  dev_min    <- min(curve_df$deviance, na.rm = TRUE)
  # Evaluate the chosen psi ONCE. Every f_cum() call is a full multi-origin simulation,
  # and calling it separately for the returned scalar and for the note would both double
  # the cost and report two slightly different numbers (the random streams desynchronise
  # once the seeded sets differ). Reuse the value already on the search path when one is
  # close enough, so the common case costs nothing at all.
  modelled_at_psi <- if (abs(curve_df$psi[.near(psi_hat)] - psi_hat) < 1e-9)
    curve_df$modelled_cum[.near(psi_hat)] else f_cum(psi_hat)

  list(psi = psi_hat, converged = converged, boundary_hit = boundary_hit,
       boundary_side = boundary_side, psi_lo = psi_lo, psi_hi = psi_hi, tol = tol,
       psi_tol = psi_tol, psi_bracket_lo = psi_bracket_lo, psi_bracket_hi = psi_bracket_hi,
       K = K, target_observed_cum = target, modelled_cum_at_psi = modelled_at_psi,
       psi_min_deviance = psi_dev, deviance_at_psi = dev_at_psi, deviance_min = dev_min,
       psi_max_auc_pr = psi_auc,
       deviance_boundary_hit = dev_at_edge, auc_pr_boundary_hit = auc_at_edge,
       criteria_disagree = isTRUE(is.finite(psi_dev) && abs(psi_dev - psi_hat) >
                                    0.1 * max(abs(psi_hat), 1e-9)),
       n_mc = n_mc, n_origins = length(origins),
       curve = curve_df,
       note = sprintf(paste0("psi=%.2f (%s%s) fitted to the held-out-origin counts: modelled %.1f vs %.0f ",
                             "observed new zones over %d held-out %d-week window(s); ",
                             "deviance-optimal psi=%.2f (deviance %.0f there vs %.0f at the count match)"),
                      psi_hat, if (converged) "converged" else "NOT converged",
                      if (boundary_hit) sprintf(", CENSORED at the %s bound", boundary_side) else "",
                      modelled_at_psi, target, length(origins), K, psi_dev, dev_min, dev_at_psi))
}

# ---------------------------------------------------------------------------
# Step 4 (header) — delta against REALISED K-week outcomes (psi retired)
# ---------------------------------------------------------------------------

#' Fit the hazard scale factor to the realised K-week invasion counts at the held-out origins.
#'
#' WHAT "held-out" AND "_oos" MEAN HERE, and what they do not. Each origin's hazard and R are
#' trained only on data as known at that origin (cascade_calibration_origins(); an origin whose
#' as-of reconstruction failed falls back to the final-count slice and is recorded as such in
#' `training_source`), so its K-week window is held out from THAT training. delta itself is
#' FITTED to those windows' outcomes: the level it produces is calibrated IN-SAMPLE on them, and
#' re-scoring the pooled count on the same origins (cascade_calibration_report(), the backtest)
#' cannot test it. The "_oos" in this function's name, and in `delta_oos`, is historical.
#'
#' WHY THIS REPLACED THE (delta at h=1, psi at K) PAIR. psi was not identifiable on this
#' outbreak: across runs it came back censored at a bound, or at a value 30% from its own root,
#' and two probabilistic criteria (deviance, AUC-PR) kept improving to the search ceiling while
#' the count match sat far below it — the documented signature of a phenomenological term
#' absorbing misspecification it cannot represent. Carrying an unidentifiable parameter into a
#' published projection is worse than not having it, so psi is fixed at 0 (no frontier
#' saturation) and ONE parameter is fitted, at the horizon the projection is actually read at.
#'
#' WHAT IS GIVEN UP, stated plainly. With delta fitted at h=1 the cascade's first week was
#' exactly the recalibrated short-horizon hazard. Fitting delta at K instead buys a multi-week
#' LEVEL matched to the K-week count at the cost of that tie. The h=1 consistency gate does NOT
#' test the level either way: its level ratio is ~delta for ANY delta (the cascade and the raw
#' short-horizon model share one fit), so it is a plumbing check, not calibration
#' (cascade_consistency_gate()). The real one-week level check is cascade_delta_h1_refit(),
#' which refits delta from the cascade's own week-1 probabilities against REALISED week-1
#' invasions at held-out origins. A large gap between that and the K-week delta is a statement
#' about the model, and belongs in the report.
#'
#' MONOTONICITY RUNS THE OPPOSITE WAY FROM psi, and this is the one thing a reader of both
#' solvers must not conflate: the per-week hazard is mu = delta * exp(eta) * sat, so the modelled
#' count is non-INcreasing in psi (cascade_fit_psi_oos, and stated correctly at its own
#' docstring) and strictly INCREASING in delta. This line said "non-DEcreasing in psi's
#' absence", which inverts the very direction it exists to keep straight. The bracket is
#' therefore [last delta still BELOW target, first delta above it], and bisection moves the LOWER
#' end up when the midpoint still under-predicts — the mirror of cascade_fit_psi_oos().
#'
#' Bisection is GEOMETRIC, not arithmetic: delta is a multiplicative scale factor, its band
#' (RECAL_BAND) spans four orders of magnitude either side of 1 and its uncertainty is reported as
#' se_log, so equal steps in log
#' delta are the natural ones. An arithmetic midpoint would spend nearly every refinement in the
#' upper half of the band.
#'
#' STOPPING RULE — on the BRACKET, not on the first acceptable count. This mirrors
#' cascade_fit_psi_oos() deliberately; keep the two in step. Stopping as soon as |gap| <= atol
#' returns the first value that lands inside the band, which is not the estimator's value: on the
#' 2026-09 frame that rule returned the very first bisection midpoint while the true root sat 30%
#' away. Refining until the bracket is narrow makes the estimate a property of the data rather
#' than of where the search happened to stop.
#'
#' @return list(delta, converged, boundary_hit, boundary_side, curve, target, ...)
cascade_fit_delta_oos <- function(origins, K = 6L, psi = 0,
                                  delta_lo = NULL, delta_hi = NULL,
                                  n_coarse = 6L, n_mc = 600L, seed = CASCADE_SEED,
                                  tol = get0("CASCADE_PSI_TOL", ifnotfound = 0.05),
                                  max_refine = 12L,
                                  delta_tol = get0("CASCADE_DELTA_BRACKET_TOL",
                                                   ifnotfound = 0.01)) {
  if (!length(origins)) return(NULL)
  K <- as.integer(K); key <- as.character(K)
  band <- get0("RECAL_BAND", ifnotfound = c(1e-4, 1e4))   # the 00_config.R value
  if (is.null(delta_lo)) delta_lo <- band[1]
  if (is.null(delta_hi)) delta_hi <- band[2]
  stopifnot(delta_lo > 0, delta_hi > delta_lo)

  # Per-horizon completeness: a recent origin can carry its 1-week outcome but not its K-week
  # one. Select the origins that do rather than aborting on the most recent.
  has_K <- vapply(origins, function(o) key %in% names(o$truth), logical(1))
  if (!any(has_K)) {
    warning(sprintf(paste0("[cascade_calib] no held-out origin carries a settled ",
                           "%d-week outcome; delta cannot be fitted on the K-week windows and ",
                           "cascade_calibrate() falls back to the short-horizon factor."), K),
            call. = FALSE)
    return(NULL)
  }
  if (!all(has_K))
    message(sprintf("[cascade_calib] delta fitted on %d of %d origin(s) — %d lack a complete %d-week outcome.",
                    sum(has_K), length(origins), sum(!has_K), K))
  origins <- origins[has_K]

  target <- sum(vapply(origins, function(o) sum(o$truth[[key]]), numeric(1)))
  if (target <= 0) {
    warning("[cascade_calib] zero realised invasions across the held-out windows; delta is not ",
            "identifiable from the K-week counts and cascade_calibrate() falls back to the ",
            "short-horizon factor.", call. = FALSE)
    return(NULL)
  }

  trace <- list()
  f_cum <- function(dl) {
    d <- .cascade_pool_po(origins, H = K, delta = dl, psi = psi, n_mc = n_mc, seed = seed)
    v <- sum(d$p)
    trace[[length(trace) + 1L]] <<- data.frame(
      delta = dl, modelled_cum = v, observed_cum = target,
      deviance = .cascade_deviance(d$p, d$y),
      auc_pr_skill = { ap <- .casc_auc_pr(d$p, d$y); br <- mean(d$y)
                       if (is.finite(ap) && br > 0) ap / br else NA_real_ })
    v
  }

  grid  <- exp(seq(log(delta_lo), log(delta_hi), length.out = max(n_coarse, 2L)))
  curve <- vapply(grid, f_cum, numeric(1))
  gap   <- curve - target                 # INCREASING in delta
  atol  <- tol * target

  boundary_hit <- FALSE; boundary_side <- NA_character_; converged <- TRUE
  delta_hat <- NA_real_; br_lo <- NA_real_; br_hi <- NA_real_
  if (all(gap >= 0)) {
    # Even the smallest admissible delta over-predicts: the root is below the band.
    delta_hat <- grid[1L]
    # A solution sitting ON the band floor is a CENSORED value whatever the residual gap.
    # Flagging it only when |gap| > atol meant a floor solution that happened to land within
    # tolerance was returned with converged = TRUE and boundary_hit = FALSE — exactly the
    # "censored boundary solution returned with no convergence flag" this module's header says
    # was fixed. `converged` is left alone: the search did terminate; the value is just pinned.
    boundary_hit <- TRUE; boundary_side <- "lower"
    if (abs(gap[1L]) > atol) {
      converged <- FALSE; boundary_hit <- TRUE; boundary_side <- "lower"
      warning(sprintf(paste0("[cascade_fit_delta_oos] at delta=%.3f (the band's floor) the cascade ",
                             "still projects %.1f new zones over the held-out %d-week window(s) vs ",
                             "%.0f observed (%.0f%% over). delta is CENSORED at the lower bound; the ",
                             "excess is in the hazard, the kernel or R, not in the scale factor."),
                     delta_hat, curve[1L], K, target,
                     100 * (curve[1L] / target - 1)), call. = FALSE)
    }
  } else if (all(gap < 0)) {
    delta_hat <- grid[length(grid)]
    boundary_hit <- TRUE; boundary_side <- "upper"; converged <- FALSE
    warning(sprintf(paste0("[cascade_fit_delta_oos] at delta=%.3f (the band's ceiling) the cascade ",
                           "projects only %.1f new zones over the held-out %d-week window(s) vs %.0f ",
                           "observed (%.0f%% under). delta is CENSORED and the projection will ",
                           "under-run; inspect the hazard and the mobility kernel."),
                   delta_hat, curve[length(curve)], K, target,
                   100 * (1 - curve[length(curve)] / target)), call. = FALSE)
  } else {
    i_lo <- max(which(gap < 0))           # last delta still BELOW the observed count
    lo <- grid[i_lo]; hi <- grid[i_lo + 1L]
    g_lo <- gap[i_lo]; g_hi <- gap[i_lo + 1L]
    for (it in seq_len(max_refine)) {
      if (log(hi / lo) <= delta_tol) break
      mid   <- sqrt(lo * hi)              # geometric: delta is a scale factor
      g_mid <- f_cum(mid) - target
      # Bisection uses only the SIGN of the gap, so Monte-Carlo noise near the root cannot stop
      # the bracket closing; it blurs only which end is finally picked, and the bracket is
      # returned so that residual width stays visible rather than implied.
      if (g_mid < 0) { lo <- mid; g_lo <- g_mid } else { hi <- mid; g_hi <- g_mid }
    }
    delta_hat <- if (abs(g_lo) <= abs(g_hi)) lo else hi
    converged <- min(abs(g_lo), abs(g_hi)) <= atol
    br_lo <- lo; br_hi <- hi
    if (!converged)
      warning(sprintf(paste0("[cascade_fit_delta_oos] bisection exhausted %d refinements without ",
                             "reaching the %.0f%% tolerance; residual %.1f zones on a target of %.0f."),
                      max_refine, 100 * tol, min(abs(g_lo), abs(g_hi)), target), call. = FALSE)
  }

  curve_df <- do.call(rbind, trace)
  curve_df <- curve_df[order(curve_df$delta), , drop = FALSE]
  curve_df <- curve_df[!duplicated(curve_df$delta), , drop = FALSE]
  curve_df$gap <- curve_df$modelled_cum - target
  # The count match is the ESTIMATOR (monotone, unique root); the deviance- and AUC-PR-optimal
  # delta on the same path are reported beside it as independent checks. Agreement on a value far
  # from the count match would mean the level and the per-zone allocation disagree, which the
  # count match alone cannot see.
  d_dev <- curve_df$delta[which.min(curve_df$deviance)]
  d_auc <- if (any(is.finite(curve_df$auc_pr_skill)))
    curve_df$delta[which.max(curve_df$auc_pr_skill)] else NA_real_
  .at_edge <- function(x) is.finite(x) &&
    (abs(log(x / delta_lo)) < 1e-9 || abs(log(x / delta_hi)) < 1e-9)
  dev_edge <- .at_edge(d_dev); auc_edge <- .at_edge(d_auc)
  if (dev_edge || auc_edge)
    warning(sprintf(paste0("[cascade_fit_delta_oos] a PROBABILISTIC criterion is pinned at a band ",
                           "edge (%s%s%s) while the count match returned %.3f. A criterion at the ",
                           "edge means the band, not the data, chose it — widen the band to test ",
                           "whether it keeps improving; do not adopt the edge value as a fit."),
                   if (dev_edge) "deviance" else "",
                   if (dev_edge && auc_edge) " and " else "",
                   if (auc_edge) "AUC-PR" else "", delta_hat), call. = FALSE)
  .near <- function(x) which.min(abs(log(curve_df$delta / x)))
  dev_at <- curve_df$deviance[.near(delta_hat)]
  dev_min <- min(curve_df$deviance, na.rm = TRUE)
  # Evaluate the chosen delta ONCE: every f_cum() call is a full multi-origin simulation, and a
  # separate evaluation for the note would both double the cost and report a second, slightly
  # different number for one quantity.
  modelled_at <- if (abs(log(curve_df$delta[.near(delta_hat)] / delta_hat)) < 1e-9)
    curve_df$modelled_cum[.near(delta_hat)] else f_cum(delta_hat)

  list(delta = delta_hat, converged = converged, boundary_hit = boundary_hit,
       boundary_side = boundary_side, delta_lo = delta_lo, delta_hi = delta_hi,
       tol = tol, delta_tol = delta_tol, bracket_lo = br_lo, bracket_hi = br_hi,
       K = K, psi_held_at = psi, target_observed_cum = target,
       modelled_cum_at_delta = modelled_at,
       delta_min_deviance = d_dev, deviance_at_delta = dev_at, deviance_min = dev_min,
       delta_max_auc_pr = d_auc,
       deviance_boundary_hit = dev_edge, auc_pr_boundary_hit = auc_edge,
       criteria_disagree = isTRUE(is.finite(d_dev) &&
                                    abs(log(d_dev / max(delta_hat, 1e-12))) > 0.1),
       n_mc = n_mc, n_origins = length(origins), curve = curve_df,
       note = sprintf(paste0("delta=%.3f (%s%s) from HELD-OUT-ORIGIN counts at psi=%g: modelled %.1f ",
                             "vs %.0f observed new zones over %d held-out %d-week window(s); ",
                             "deviance-optimal delta=%.3f (deviance %.0f there vs %.0f at the count match)"),
                     delta_hat, if (converged) "converged" else "NOT converged",
                     if (boundary_hit) sprintf(", CENSORED at the %s bound", boundary_side) else "",
                     psi, modelled_at, target, length(origins), K, d_dev, dev_min, dev_at))
}

# ---------------------------------------------------------------------------
# Step 6 (header) — residual check under the FITTED delta (count in-sample)
# ---------------------------------------------------------------------------

#' Residual check of the fitted cascade on the calibration origins. The pooled COUNT is
#' IN-SAMPLE for delta (delta was fitted to it); only the rank/discrimination columns are not
#' what delta was fitted to. See the note above `resid_ok` below.
#'
#' @return list(per_origin, pooled, reliability, pass, note)
cascade_calibration_report <- function(origins, delta, psi, K = 6L,
                                       n_mc = 600L, seed = CASCADE_SEED,
                                       count_tol = 0.25, n_bins = 5L) {
  if (!length(origins)) return(NULL)
  K <- as.integer(K)
  origins <- Filter(function(o) as.character(K) %in% names(o$truth), origins)
  if (!length(origins)) return(NULL)
  d <- .cascade_pool_po(origins, H = K, delta = delta, psi = psi, n_mc = n_mc, seed = seed)
  if (!nrow(d)) return(NULL)

  per_origin <- do.call(rbind, lapply(split(d, d$origin), function(g) {
    br <- mean(g$y); ap <- .casc_auc_pr(g$p, g$y)
    data.frame(origin = g$origin[1], K = K, n_atrisk = nrow(g), n_events = sum(g$y),
               predicted_new = sum(g$p), observed_new = sum(g$y),
               count_ratio = sum(g$p) / max(sum(g$y), 1e-9),
               auc_pr = ap, auc_pr_skill = if (is.finite(ap) && br > 0) ap / br else NA_real_,
               mean_rank_of_truth = { rk <- rank(-g$p, ties.method = "average")
                                      if (any(g$y == 1L)) mean(rk[g$y == 1L]) else NA_real_ },
               stringsAsFactors = FALSE)
  }))
  rownames(per_origin) <- NULL

  br <- mean(d$y); ap <- .casc_auc_pr(d$p, d$y)
  pooled <- data.frame(
    K = K, n_origins = length(unique(d$origin)), n_atrisk = nrow(d), n_events = sum(d$y),
    predicted_new = sum(d$p), observed_new = sum(d$y),
    count_ratio = sum(d$p) / max(sum(d$y), 1e-9),
    auc_pr = ap, auc_pr_skill = if (is.finite(ap) && br > 0) ap / br else NA_real_,
    deviance = .cascade_deviance(d$p, d$y),
    delta = as.numeric(delta), psi = psi)

  # Reliability on equal-count bins: with a ~3% base rate, equal-WIDTH bins put
  # almost every zone in the first bin and say nothing.
  q  <- unique(stats::quantile(d$p, seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE))
  rel <- if (length(q) >= 3L) {
    b <- cut(d$p, breaks = q, include.lowest = TRUE, labels = FALSE)
    do.call(rbind, lapply(sort(unique(b[!is.na(b)])), function(k) {
      idx <- which(b == k)
      data.frame(bin = k, n = length(idx), mean_pred = mean(d$p[idx]),
                 obs_freq = mean(d$y[idx]), k_events = sum(d$y[idx]))
    }))
  } else NULL

  # THE LEVEL MATCH IS IN-SAMPLE FOR DELTA, and must be labelled as such.
  #
  # cascade_fit_delta_oos() chooses delta by bisecting the pooled predicted count at exactly
  # these origins and this K onto the observed count. Re-evaluating that same quantity here is
  # not a validation — with the same seed and n_mc it is the SAME deterministic computation, and
  # the shipped artifacts showed it: cascade_delta_fit.json's modelled_cum_at_delta = 73.1167 and
  # cascade_calibration_pooled.csv's predicted_new = 73.1167, identical to 4 decimal places, with
  # "verification_pass": true. The search guarantees the report; the check could not fail.
  #
  # What IS informative, and is retained: (a) the count residual as a NUMBER, computed at an
  # INDEPENDENT Monte-Carlo seed so it at least carries MC error rather than being bit-identical;
  # (b) the rank/discrimination columns (AUC-PR skill, rank-of-truth), which delta was NOT
  # fitted to. Stated precisely: at one week delta is a monotone transform of every zone's
  # hazard and cannot reorder them; over K weeks it can, but only second-order, through the
  # cascade (a zone reached in two hops gains more from a larger delta than one reached in one).
  # The search chose delta to match the COUNT, never to improve a ranking, so these columns are
  # the out-of-sample part of the check. The field is named for what it is.
  resid_ok <- is.finite(pooled$count_ratio) && abs(pooled$count_ratio - 1) <= count_tol
  if (!resid_ok)
    warning(sprintf(paste0("[cascade_calib] IN-SAMPLE level residual at K=%d is outside tolerance: ",
                           "predicted %.1f vs observed %.0f (ratio %.2f, tol +/-%.0f%%). Since delta ",
                           "was fitted to match exactly this quantity, a residual this large means ",
                           "the SEARCH did not converge, not that the model was validated."),
                   K, pooled$predicted_new, pooled$observed_new, pooled$count_ratio,
                   100 * count_tol), call. = FALSE)

  list(per_origin = per_origin, pooled = pooled, reliability = rel, pooled_po = d,
       # `pass` retained for backward compatibility with consumers, but it is a residual check.
       pass = resid_ok, count_level_in_sample = TRUE,
       count_residual_ok = resid_ok, count_tol = count_tol,
       note = sprintf(paste0("K=%d: predicted %.1f vs observed %.0f new zones (ratio %.2f, tol +/-%.0f%%) %s ",
                             "[IN-SAMPLE for delta: delta was fitted to match this count]; ",
                             "AUC-PR skill %.1fx [not fitted to: delta was chosen on the count]"),
                      K, pooled$predicted_new, pooled$observed_new, pooled$count_ratio,
                      100 * count_tol, if (resid_ok) "OK" else "OFF", pooled$auc_pr_skill))
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

#' Calibrate the cascade end to end and return the parameters to simulate with.
#'
#' @param K held-out window length (weeks) used to FIT delta and to verify.
#' @param n_mc_* Monte-Carlo sizes for the h=1 delta refit and the K-week delta
#'   search. The verification size DEFAULTS TO the search size and should stay tied to it:
#'   both pool the same origins at the same delta and psi, so an equal Monte-Carlo size keeps
#'   the search's `modelled_cum_at_delta` and the verification's `predicted_new` comparable.
#'   They are NOT one number: the verification deliberately runs at an OFFSET SEED
#'   (seed + 99991) so the reported residual is an independent Monte-Carlo realisation and MC
#'   error is exposed rather than hidden — see the "INDEPENDENT SEED" note at that call. Expect
#'   the two to differ by Monte-Carlo error, and read `count_residual_ok` with that in mind.
#'   The search evaluates the cascade ~13-20 times, so its size dominates the cost.
#' @return list(delta, delta_info, delta_oos, delta_prior, psi, psi_fit, report,
#'   origins, note)
cascade_calibrate <- function(layer, zone_province = NULL, K = 6L,
                              cutoffs_from_end = c(12L, 11L, 10L, 9L, 8L),
                              # n_mc_report FOLLOWS n_mc_search rather than carrying its own
                              # literal. Matching the two DEFAULTS was not enough: run_cascade.R
                              # sets the search size from an env var, so the call site silently
                              # restored the mismatch and the change was inert. Tying them means
                              # both move together whatever any caller overrides. While they
                              # differed, the report printed two numbers for one quantity —
                              # 72.96 from the search and 71.7 from the verification.
                              #
                              # RENAMED from n_mc_psi: it now sizes the DELTA search, psi having
                              # been retired. A parameter whose name says psi while it drives
                              # delta is exactly the kind of stale label that misleads a later
                              # reader. run_cascade.R was the only caller passing it.
                              n_mc_delta = 600L, n_mc_search = 600L, n_mc_report = n_mc_search,
                              kernel = CASCADE_KERNEL, gt = CASCADE_GT,
                              iter_fit = 600L, seed = CASCADE_SEED,
                              transfer_tol = 0.35) {
  K <- as.integer(K)
  # delta comes from the prequential maximum-likelihood factor 16b fits for the
  # cascade's own covariate-bearing model; cascade_fit_delta() carries the
  # fallbacks and the diagnostics.
  delta <- cascade_fit_delta(cov_model = cascade_hazard_label(kernel))
  delta_info <- cascade_short_horizon_delta_info(delta)

  origins <- cascade_calibration_origins(
    layer, zone_province = zone_province, cutoffs_from_end = cutoffs_from_end,
    truth_horizons = unique(c(1L, K)), kernel = kernel, gt = gt,
    iter_fit = iter_fit, linelist = layer$ll, delta_imp = as.numeric(delta))
  # No origins means NO calibration: delta_oos is NULL, which run_cascade.R treats as a failed
  # calibration (fatal unless CASCADE_CALIBRATE=0 was set deliberately). delta_prior is carried
  # so the record is complete either way.
  if (!length(origins))
    return(list(delta = delta, delta_info = delta_info, psi = CASCADE_PSI,
                psi_fit = NULL, delta_oos = NULL, delta_prior = as.numeric(delta),
                report = NULL, origins = list(), h1_refit = NULL,
                note = "no held-out origins; delta from 16b, psi at its configured default"))

  # Verify the delta transfer: refit it from the cascade's OWN week-1 hazard
  # against realised week-1 invasions. A large disagreement means the two
  # implementations do not share a level and delta is not transferable.
  #
  # TWO RATIOS, EACH NAMED FOR ITS REFERENCE. At this point `delta` is still 16b's SHORT-HORIZON
  # factor; the K-week fit below replaces it. A single `h1_refit_ratio` used to be stored here,
  # i.e. h1 / short-horizon factor, and 36_report.R printed it as the ratio "to the K-week factor"
  # and "to the deployed delta" — a different number whenever the two factors differ. The ratio
  # to the short-horizon factor is stored now (`h1_refit_ratio_short`); the ratio to the DEPLOYED
  # factor (`h1_refit_ratio_deployed`) is stored only once that factor is final, below. The
  # ambiguous `h1_refit_ratio` is no longer written.
  h1 <- cascade_delta_h1_refit(origins, n_mc = n_mc_delta, seed = seed)
  if (!is.null(h1) && is.finite(h1$delta) && is.finite(as.numeric(delta))) {
    ratio_short <- h1$delta / as.numeric(delta)
    if (abs(ratio_short - 1) > transfer_tol)
      warning(sprintf(paste0("[cascade_calib] delta TRANSFER check: 16b's short-horizon factor is ",
                             "%.3f but refitting from the cascade's own week-1 hazard on the held-out ",
                             "origins gives %.3f (ratio %.2f, tolerance +/-%.0f%%). The two ",
                             "implementations do not share a level at h=1. Neither value is the one ",
                             "projected with: delta is next fitted on the %d-week counts and that ",
                             "factor is deployed (the short-horizon factor is carried only if that fit ",
                             "fails). The discrepancy is real and is recorded."),
                     as.numeric(delta), h1$delta, ratio_short, 100 * transfer_tol, K),
              call. = FALSE)
    delta_info$h1_refit_delta       <- h1$delta
    delta_info$h1_refit_ratio_short <- ratio_short
    delta_info$h1_refit_lo    <- h1$lo
    delta_info$h1_refit_hi    <- h1$hi
    delta_info$h1_refit_n_events <- h1$n_events
  }

  # ONE fitted parameter. psi is held at its configured value (0 — no frontier saturation; see
  # 30_projection_config.R) rather than fitted, because it was not identifiable on this outbreak;
  # delta is fitted against the realised K-week counts instead of being carried over from h=1.
  # `delta_prior` (16b's short-horizon factor) is kept: it built the import term in the origins
  # above, and it is the natural comparison for what the multi-week fit returns.
  psi <- as.numeric(get0("CASCADE_PSI", ifnotfound = 0))
  delta_prior <- delta
  delta_oos <- cascade_fit_delta_oos(origins, K = K, psi = psi,
                                     n_mc = n_mc_search, seed = seed)
  if (!is.null(delta_oos) && is.finite(delta_oos$delta)) {
    delta <- delta_oos$delta
    delta_info$delta_short_horizon <- as.numeric(delta_prior)
    delta_info$delta_oos           <- delta_oos$delta
    # delta / lo / hi MUST follow the estimator label. They held the 16b SHORT-HORIZON values
    # while $estimator was overwritten to claim the K-week fit, so
    # cascade_calibration.json reported delta_info$delta = 0.430 [0.311, 0.601] for an estimator
    # that in fact returned 0.389. Rename the prior's interval rather than reuse it: an interval
    # from a different estimator on a different horizon does not describe this point estimate.
    delta_info$delta               <- as.numeric(delta_oos$delta)
    delta_info$delta_prior_lo      <- delta_info$lo
    delta_info$delta_prior_hi      <- delta_info$hi
    delta_info$lo                  <- NA_real_
    delta_info$hi                  <- NA_real_
    delta_info$delta_oos_ratio     <- delta_oos$delta / max(as.numeric(delta_prior), 1e-12)
    delta_info$estimator           <- sprintf("%d-week count match at held-out origins", K)
    # PROVENANCE follows the estimate: `source` and `cal_in_large` described the short-horizon
    # prior's row (16b), so beside the K-week delta they read as its provenance. Kept, renamed.
    delta_info$delta_prior_source       <- delta_info$source
    delta_info$delta_prior_cal_in_large <- delta_info$cal_in_large
    delta_info$source                   <- "cascade_fit_delta_oos (33b)"
    delta_info$cal_in_large             <- NA_real_
    # The bracket is NUMERICAL width, not statistical uncertainty: the Monte-Carlo component of
    # the search dominates it by roughly an order of magnitude (measured for the psi solver on
    # the same frames), so it must never be quoted as delta's precision.
    #
    # `delta_info$se_log` is therefore LEFT as it was — the zone-clustered bootstrap SE that
    # 16b's prequential short-horizon estimator reported for delta_prior, NOT an SE for the
    # K-week fit, which the count-match solver does not produce. simulate_cascade() draws delta
    # per parameter group on that se_log, so the projection still carries an estimated interval
    # rather than a fixed value. APPROXIMATION, stated: this assumes the K-week factor's
    # RELATIVE precision resembles the short-horizon factor's. Both are hazard-scale factors on
    # the same data, so that is reasonable, but it is an assumption and not a measurement.
    delta_info$bracket_lo <- delta_oos$bracket_lo
    delta_info$bracket_hi <- delta_oos$bracket_hi
    delta_info$boundary_hit <- delta_oos$boundary_hit
    if (isTRUE(delta_oos$boundary_hit))
      warning(sprintf(paste0("[cascade_calib] delta is CENSORED at its %s bound (%.3f): the band, ",
                             "not the %d-week calibration windows, set it, so the 13-week level is ",
                             "NOT supported by the calibration and the layer is not fit for ",
                             "publication until this is resolved."),
                      delta_oos$boundary_side, delta, K), call. = FALSE)
  } else {
    warning(sprintf(paste0("[cascade_calib] delta could not be fitted on the %d-week calibration ",
                           "windows; the short-horizon factor %.3f is returned unchanged and the ",
                           "13-week LEVEL is therefore uncalibrated at the horizon it is read at."),
                    K, as.numeric(delta_prior)), call. = FALSE)
  }
  # The week-1 refit against the DEPLOYED factor, computed only now that `delta` is final: the
  # K-week factor when that fit succeeded, else the short-horizon factor carried in its place
  # (identical to h1_refit_ratio_short in that case). See the note at the week-1 refit above.
  if (is.finite(delta_info[["h1_refit_delta"]] %||% NA_real_) && is.finite(as.numeric(delta)) &&
      as.numeric(delta) > 0)
    delta_info$h1_refit_ratio_deployed <- delta_info[["h1_refit_delta"]] / as.numeric(delta)
  psi_fit <- NULL   # retired: psi is no longer fitted (cascade_fit_psi_oos remains for sweeps)

  # INDEPENDENT SEED. With the search's own seed this call reproduces the bisection's objective
  # bit-for-bit (verified in the shipped artifacts). Offsetting it makes the reported residual an
  # independent Monte-Carlo realisation, which at least exposes MC error; it does not make the
  # level check out-of-sample, and the report says so.
  report <- cascade_calibration_report(origins, delta = as.numeric(delta), psi = psi,
                                       K = K, n_mc = n_mc_report,
                                       seed = as.integer(seed) + 99991L)

  # Per-origin source-dynamics diagnostics, written out so that the burn-in rule and the R
  # values each origin was calibrated on are auditable artifacts rather than log lines.
  origin_diag <- tryCatch(dplyr::bind_rows(lapply(origins, function(o) tibble::tibble(
      cutoff_date = o$cutoff_date, week_index = o$Ccol, R_nat = o$R_nat,
      rt_source = o$rt_source %||% NA_character_,
      n_with_cases = o$n_with_cases, n_prior_dominated = o$n_prior_dominated,
      import_share = o$import_share, n_atrisk = o$n_atrisk,
      horizons = paste(o$complete_horizons, collapse = "/"),
      events = paste(o$n_events, collapse = "/"),
      training_source = o$training_source %||% NA_character_,
      vintage_date = as.Date(o$vintage_date %||% NA),
      rolling_predictors = isTRUE(o$rolling_predictors)))), error = function(e) NULL)

  list(delta = delta, delta_info = delta_info, psi = psi, psi_fit = psi_fit,
       delta_oos = delta_oos, delta_prior = as.numeric(delta_prior),
       report = report, origins = origins, h1_refit = h1,
       origin_diagnostics = origin_diag,
       # TRUE only when every origin's hazard was refit on the validated model's rolling design.
       rolling_predictors = length(origins) > 0L &&
         all(vapply(origins, function(o) isTRUE(o$rolling_predictors), logical(1))),
       min_weeks_from_start = get0("CASCADE_CALIB_MIN_WEEKS", ifnotfound = 7L),
       note = sprintf("delta=%.3f (%s%s) | short-horizon delta=%.3f%s | psi=%g (FIXED, not fitted) | %s",
                      as.numeric(delta),
                      if (is.null(delta_oos)) "short-horizon factor, K-week fit unavailable"
                      else sprintf("held-out %d-week count match, %s", K,
                                   if (isTRUE(delta_oos$converged)) "converged" else "NOT converged"),
                      if (!is.null(delta_oos) && isTRUE(delta_oos$boundary_hit))
                        sprintf(", CENSORED at the %s bound", delta_oos$boundary_side) else "",
                      as.numeric(delta_prior),
                      if (!is.null(h1) && is.finite(h1$delta))
                        sprintf(", week-1 refit %.3f", h1$delta) else "",
                      psi,
                      if (is.null(report)) "no residual check" else report$note))
}

message("[33b] cascade calibration loaded — ONE fitted parameter: delta, fitted to the realised ",
        "K-week counts at as-of origins (an in-sample level fit; psi retired, fixed at 0).")
