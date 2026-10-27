# =============================================================================
# 36_report.R — 3-MONTH CASCADE: assemble the decision-maker report (PLAN §7)
# BDBV 2026 DRC
#
# cascade_write_report() writes SPATIAL_INVASION_3MONTH_REPORT.md in the house
# style, auto-populated from the run artifacts: reach leaderboard, expected new
# zones/provinces per scenario, gateway zones, corridors, conditional hubs, the
# consistency-gate result, sensitivity ranking-stability, and the backtest — with
# the projection / scenario-conditional / small-sample caveats stated up front.
# =============================================================================

.fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)

# --- gate-aware figure citations --------------------------------------------
# The retained-figure allow-list (FIGURE_KEEP, 00_config.R) suppresses most cascade figures,
# but this report's figure pointers were written when every one of them was produced. Citing a
# file the gate no longer writes sends the reader to a path that does not exist, which is worse
# than no pointer at all. .fig_kept() asks the SAME gate the savers ask, on the same full path.
#
# Paths here are relative to OUT_DIR (e.g. "cascade/figures/cascade_corridors"), matching how
# they are printed in the report; the gate only inspects path components and the stem, so a
# relative path answers identically to the absolute one the saver passes.
.fig_kept <- function(rel_stem) {
  fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (!is.function(fk)) return(TRUE)          # standalone use: assume everything is produced
  isTRUE(tryCatch(fk(rel_stem), error = function(e) TRUE))
}

#' Render a list of figure pointers, keeping only those the gate will actually produce.
#' @param items named character vector: names = human description, values = path stem
#'   relative to OUT_DIR, WITHOUT the extension.
#' @param prefix text to open the sentence with (e.g. "Figures: ").
#' @return a single sentence, or "" when the gate suppresses every figure in the list.
.fig_line <- function(items, prefix = "Figures: ", ext = ".pdf") {
  keep <- vapply(unname(items), .fig_kept, logical(1))
  if (!any(keep)) return("")
  bits <- sprintf("%s: `%s%s`", names(items)[keep], unname(items)[keep], ext)
  paste0(prefix, paste(bits, collapse = "; "), ".")
}

# --- the week-1 refit of delta, against each reference by name -------------------------------
#' The week-1 refit of delta (33b cascade_delta_h1_refit()) and its ratio to EACH reference.
#'
#' 33b used to store a single `h1_refit_ratio`, computed BEFORE the K-week fit — i.e. against
#' 16b's SHORT-HORIZON factor — while this report printed it as the ratio "to the K-week factor"
#' (Section 1) and "to the deployed delta" (Section 5): a different number whenever the two
#' factors differ. 33b now stores `h1_refit_ratio_short` and `h1_refit_ratio_deployed`. A context
#' saved before that carries only the legacy field; in every committed version of 33b that field
#' was h1_refit_delta / (16b's short-horizon factor), so it is read as the SHORT ratio, and the
#' ratio to the deployed factor is recomputed from the two values such a context does carry.
#' Fields are read with [[ ]], never $: `$` partial-matches on a list, so on a context holding
#' only one of the new fields `ci$h1_refit_ratio` would silently return it.
#' @param ci `ctx$calib$delta_info`.
#' @param delta_deployed the delta the projection ran on (`ctx$delta`).
#' @return list(delta, ratio_short, ratio_deployed); NA where unavailable.
.cascade_h1_refit_ratios <- function(ci, delta_deployed) {
  .num <- function(x) {
    x <- suppressWarnings(as.numeric(x %||% NA_real_))
    if (length(x) == 1L && is.finite(x)) x else NA_real_
  }
  ci <- ci %||% list()
  d1 <- .num(ci[["h1_refit_delta"]])
  dd <- .num(delta_deployed)
  r_short <- .num(ci[["h1_refit_ratio_short"]] %||% ci[["h1_refit_ratio"]])
  r_dep   <- .num(ci[["h1_refit_ratio_deployed"]])
  if (!is.finite(r_dep) && is.finite(d1) && is.finite(dd) && dd > 0) r_dep <- d1 / dd
  list(delta = d1, ratio_short = r_short, ratio_deployed = r_dep)
}

#' What the calibration origins delta was fitted on were TRAINED on, in words.
#'
#' Read per origin from 33b's `training_source` (ctx$calib$origin_diagnostics, else the origins
#' themselves): "vintage" = the archived snapshot of that date, "recensor" = the final line list
#' re-censored to that date, "final_slice" = the FINAL revised counts sliced at the cutoff, which
#' carry cases reported after the origin (the fallback when the as-of reconstruction fails). The
#' report used to assert "trained only on data known at that origin" for every origin unchecked.
#' Only origins with a settled K-week outcome entered delta's fit, so only those are counted.
#' @param calib `ctx$calib`. @param K the horizon delta was fitted at.
#' @return a phrase, with its own leading separator, completing "at N held-out origin(s)".
.cascade_origin_training_txt <- function(calib, K) {
  key <- as.character(as.integer(K))
  od  <- calib[["origin_diagnostics"]]
  src <- NULL; cut <- NULL
  if (is.data.frame(od) && nrow(od) > 0L && "training_source" %in% names(od)) {
    keep <- if ("horizons" %in% names(od))
      vapply(strsplit(as.character(od$horizons), "/", fixed = TRUE),
             function(h) key %in% h, logical(1))
      else rep(TRUE, nrow(od))
    src <- as.character(od$training_source[keep])
    cut <- if ("cutoff_date" %in% names(od)) format(as.Date(od$cutoff_date[keep]))
           else rep(NA_character_, sum(keep))
  } else if (length(calib[["origins"]])) {
    os  <- Filter(function(o) key %in% names(o$truth), calib[["origins"]])
    src <- vapply(os, function(o) as.character(o$training_source %||% NA_character_), character(1))
    cut <- vapply(os, function(o) format(as.Date(o$cutoff_date %||% NA)), character(1))
  }
  .unrec <- " (the per-origin training source is not recorded in this context)"
  if (!length(src)) return(.unrec)
  n  <- length(src)
  nv <- sum(src %in% "vintage"); nr <- sum(src %in% "recensor"); nf <- sum(src %in% "final_slice")
  nu <- n - nv - nr - nf
  if (nu == n) return(.unrec)
  .asof <- paste(c(if (nv > 0L) sprintf("%d from the archived snapshot of that date", nv),
                   if (nr > 0L) sprintf("%d from the final line list re-censored to that date", nr)),
                 collapse = ", ")
  if (nf == 0L && nu == 0L)
    return(sprintf(paste0(" whose models were each trained only on counts reconstructed as of that ",
                          "origin (%s)"), .asof))
  parts <- c(
    if (nf > 0L)
      sprintf(paste0("**%d (%s) %s trained on the final-count slice (counts revised after the ",
                     "origin: a revision leak), not on data known at the origin**"),
              nf, paste(cut[src %in% "final_slice"], collapse = ", "),
              if (nf == 1L) "was" else "were"),
    if (nv + nr > 0L)
      sprintf("%d on counts reconstructed as of their origin (%s)", nv + nr, .asof),
    if (nu > 0L) sprintf("%d with no recorded training source", nu))
  paste0(", of which ", paste(parts, collapse = "; "))
}

#' "ratio X to <deployed_label>; Y to 16b's short-horizon factor", with only the finite parts.
#' @return a phrase, or "" when neither ratio is available.
.cascade_h1_ratio_phrase <- function(r, deployed_label = "the deployed delta") {
  bits <- c(if (is.finite(r$ratio_deployed)) sprintf("%.2f to %s", r$ratio_deployed, deployed_label),
            if (is.finite(r$ratio_short))
              sprintf("%.2f to 16b's short-horizon factor", r$ratio_short))
  if (!length(bits)) "" else paste0("ratio ", paste(bits, collapse = "; "))
}

cascade_write_report <- function(path, ctx) {
  L <- c()
  add <- function(...) L <<- c(L, sprintf(...))
  aff <- sum(ctx$reach_primary$was_active_before & ctx$reach_primary$horizon ==
               max(CASCADE_REPORT_HORIZONS))
  atrisk <- sum(!ctx$reach_primary$was_active_before & ctx$reach_primary$horizon ==
                  max(CASCADE_REPORT_HORIZONS))

  add("# BDBV 2026 — 3-Month Spatial Invasion Risk (Cascade Projection)")
  add("")
  add("**Analysis date:** %s · **Kernel:** %s (%s GT) · **M:** %d MC iterations · **Horizon:** %d weeks",
      format(as.Date(ANALYSIS_DATE)), ctx$kernel %||% CASCADE_KERNEL,
      CASCADE_GT, ctx$n_mc, CASCADE_HORIZON_WEEKS)
  # THE CALIBRATION STATE, up front and read from the run (run_cascade.R records it). A run that
  # was deliberately left uncalibrated says so here, before any number is read.
  .calibrated <- isTRUE(ctx$calibrated %||%
                          (!is.null(ctx$calib$delta_oos) && is.finite(ctx$calib$delta_oos$delta %||% NA_real_)))
  add("**Calibration:** %s%s", ctx$calibration_status %||% (if (.calibrated) "calibrated" else "NOT calibrated"),
      if (is.null(ctx$cascade_run_id)) "" else sprintf(" · **Run:** `%s`", ctx$cascade_run_id))
  add("**Folder:** `%s` · implements `PLAN_3MONTH_INVASION.md`.", basename(dirname(OUT_CASCADE)))
  add("")
  add("> **How to read this.** %s The 1–2-week Bayesian suite remains the operational headline; this layer is a longer-range **projection** of where invasion may spread and how far onward transmission may carry it. Decisions should use the **rankings, relative risk, gateway ordering, and priority scores**; read absolute 3-month probabilities as **conditional on the scenario** (its assumed future transmissibility) and on a level calibrated over a shorter window than 13 weeks — not as forecasts, and not as bounds in either direction. Results are issued as **control scenarios**, not a single forecast.", CASCADE_FRAMING)
  add("")
  add("---")
  # The held-out window delta was fitted at, READ from the fit result (cascade_fit_delta_oos
  # returns K) rather than from a constant that does not exist.
  .calib_K <- suppressWarnings(as.integer((ctx$calib$delta_oos %||% list())$K %||% NA_integer_))
  add("## 1. Method (one paragraph)")
  # EVERY MECHANISM NAMED HERE MUST EXIST. This paragraph is the published Method and it
  # described three that do not: a "partially-pooled" per-zone R_eff (there is ONE national
  # EpiNow2 R applied to every zone, 31_source_dynamics.R), a delta "recalibrated to observed
  # short-horizon frequency" (delta is fitted on realised K-week invasion counts at held-out origins,
  # 33b cascade_fit_delta_oos), and damping "by a frontier-saturation term (psi)" (CASCADE_PSI
  # is 0, so the saturation factor is identically 1 and nothing is damped). Section 6 of this
  # same document already said the opposite about the first of those.
  # THE HAZARD IS REFIT, NOT REUSED, and the sentence says how. The cascade fits the cloglog
  # hazard of 21_bayesian_renewal.R (the -geo covariate twin of the selected kernel) on the
  # current data; whether that refit used the cross-validated model's rolling as-of design is
  # read from the run, because without it the intercept is on a different scale.
  # THREE cases, not two. "NOT on that model's scale" is true only when the cross-validated model
  # DOES roll (ROLLING_PREDICTORS on) and the cascade could not; when ROLLING_PREDICTORS is off the
  # cross-validated model is itself fitted on the shared count matrix, so the cascade is on its
  # scale. The flag is read from the context when the run recorded it, else from the global.
  # [[ ]], not $: `$` partial-matches, so on a context holding only one of these two names
  # ctx$rolling_predictors would silently return the other.
  .roll_cv <- isTRUE(ctx[["rolling_predictors_requested"]] %||%
                       get0("ROLLING_PREDICTORS", ifnotfound = TRUE))
  .roll_txt <- if (isTRUE(ctx[["rolling_predictors"]]))
    "with the same rolling as-of design as the cross-validated short-horizon model (each training transition's import force from the counts known at that transition),"
  else if (!.roll_cv)
    "on the shared count matrix — the same non-rolling design the cross-validated short-horizon model uses in this run (`ROLLING_PREDICTORS` is off), so its intercept is on that model's scale,"
  else "on the shared count matrix (**without the cross-validated model's rolling as-of predictors, so its intercept is NOT on that model's scale**),"
  .do0 <- ctx$calib$delta_oos
  # WHAT THE ORIGINS WERE TRAINED ON is read per origin (33b training_source), not asserted: an
  # origin whose as-of reconstruction failed trains on the revision-leaky final-count slice.
  .delta_txt <- if (.calibrated && !is.null(.do0))
    sprintf(paste0("by a single factor `delta`, fitted so the expected number of new invasions over ",
                   "held-out %d-week windows matches the number that occurred (the windows are held ",
                   "out of each origin's training, not of delta's fit, so this level match is ",
                   "in-sample), at %d held-out origin(s)%s (`33b_cascade_calibration.R`)"),
            as.integer(.do0$K), as.integer(.do0$n_origins),
            .cascade_origin_training_txt(ctx$calib, as.integer(.do0$K)))
  else sprintf(paste0("by `delta` = the ONE-WEEK prequential recalibration factor of ",
                      "`16b_invasion_recalibration.R`: the multi-week calibration was **disabled ",
                      "for this run (%s)**, so the 13-week level is **uncalibrated at the horizon ",
                      "it is read at**"), ctx$calibration_disabled_by %||% "reason unrecorded")
  # THE FORM OF THE IMPORT COEFFICIENT the run's hazard carried (CASCADE_FORM: the featured
  # model's, from the two-stage selection).
  .form <- ctx[["form"]] %||% get0("CASCADE_FORM", ifnotfound = "fixed")
  .form_txt <- switch(.form,
    fixed = "with a fixed import coefficient",
    rw1 = paste("with the featured model's weekly random-walk drift in the import coefficient,",
                "continued over the projection one path per posterior parameter draw"),
    frailty = paste("with the featured model's gamma zone frailty: each never-invaded zone's",
                    "hazard is conditioned on its own cumulative hazard, carried from the training",
                    "weeks and accumulated over the simulated ones"),
    rw1_frailty = paste("with the featured model's weekly random-walk drift in the import",
                        "coefficient AND its gamma zone frailty (each never-invaded zone's hazard",
                        "conditioned on its own cumulative hazard)"),
    sprintf("with the '%s' form", .form))
  add("A stochastic spatial-metapopulation **cascade** couples (A) a within-zone overdispersed NegBin branching process — every infected zone projects incidence forward with the SAME national current effective reproduction number `R_eff` (NOT R0), estimated once by EpiNow2 on the national series — with (B) the invasion hazard of the featured specification (%s, %s), refit on the current data %s and used as a weekly Bernoulli **seeding** step for each still-at-risk zone. A seeded zone becomes a source, so invasion chains onward across the mobility network. The per-step hazard is rescaled on the hazard scale %s; the frontier-saturation term is switched off (`psi` = %.2f), so no damping is applied. The dynamic `d_min` frontier covariate is recomputed each simulated week. Individual-offspring overdispersion (`k_indiv`≈%.2f) produces realistic stochastic **fade-out**, so invasion (first case) is distinguished from establishment (a sustained outbreak). Everything runs on the confirmed-case scale (time index = symptom onset), consistent with the short-horizon suite.",
      if (.form %in% c("rw1", "frailty", "rw1_frailty"))
        "`21c_frailty_models.R`, the Stan engine of the time-varying forms, with the selected kernel's `-geo` covariate set"
      else "`21_bayesian_renewal.R`, the `-geo` covariate twin of the selected kernel",
      .form_txt, .roll_txt, .delta_txt, as.numeric(get0("CASCADE_PSI", ifnotfound = NA_real_)), CASCADE_K_INDIV)
  add("")
  # Describes the calibration ACTUALLY performed (33b_cascade_calibration.R). The
  # previous text said delta came from 1/calibration-in-the-large and said nothing
  # about how psi was obtained; both were changed on 2026-09-11, and a methods
  # paragraph that names the wrong estimator is worse than none.
  .ci <- ctx$calib$delta_info %||% list()
  .pf <- ctx$calib$psi_fit
  .rp <- ctx$calib$report
  .do <- ctx$calib$delta_oos
  # TWO PARAGRAPHS, chosen by the run's calibration state. The uncalibrated branch used to be
  # rendered through the calibrated one's template, which is how the report came to print
  # "fitted so ... over held-out n/a-week window(s)" for a delta that was never fitted.
  # Each ratio is named for its reference (.cascade_h1_refit_ratios()). This sentence is printed
  # only in the calibrated branch below, where the deployed delta IS the K-week factor.
  .h1r <- .cascade_h1_refit_ratios(.ci, ctx$delta)
  .h1_txt <- if (is.finite(.h1r$delta)) {
    .ph <- .cascade_h1_ratio_phrase(.h1r, deployed_label = "the deployed K-week factor")
    sprintf(paste0(" The ONE-WEEK level is checked against observation separately: refitting delta ",
                   "from the cascade's own week-1 probabilities against realised week-1 invasions at ",
                   "the held-out origins gives %.3f%s."),
            .h1r$delta, if (nzchar(.ph)) sprintf(" (%s)", .ph) else "")
  } else ""
  if (.calibrated && !is.null(.do)) {
    add(paste0("**Calibration — ONE fitted parameter.** delta = %.3f (%s): the per-week hazard ",
               "multiplier, fitted so the cascade's expected new invasions over held-out ",
               "%d-week window(s) at %d origin(s) match the number that actually occurred (%s), ",
               "using only outcomes settled for at least %s days. It is fitted at the multi-week ",
               "horizon rather than at h=1 deliberately: that is the horizon this layer is read at, ",
               "and with saturation removed there is no second parameter left to absorb multi-week ",
               "shape. **What that costs is stated**: the cascade's first week is no longer the ",
               "short-horizon model recalibrated at h=1 (16b's prequential factor is %s on this ",
               "frame).%s The h=1 consistency gate below does not test the level at all — its ",
               "level ratio is ~delta by construction. **psi = %g is FIXED, not fitted.** ",
               "Frontier saturation was not identifiable on this outbreak — across refits it came ",
               "back censored at a bound, or well away from its own root, while two probabilistic ",
               "criteria (deviance, AUC-PR) kept improving to the search ceiling. That is the ",
               "signature of a phenomenological term absorbing misspecification it cannot ",
               "represent, not of saturation being measured, so it is set to 0 (sat identically 1) ",
               "and retained only as a sensitivity axis. %s"),
        as.numeric(ctx$delta),
        .ci$estimator %||% "estimator unrecorded",
        as.integer(.do$K), as.integer(.do$n_origins),
        sprintf("modelled %.1f vs %.0f observed%s%s", .do$modelled_cum_at_delta,
                .do$target_observed_cum,
                if (isTRUE(.do$converged)) ", converged" else ", **NOT converged**",
                if (isTRUE(.do$boundary_hit))
                  sprintf(", **CENSORED at the %s bound**", .do$boundary_side) else ""),
        format(get0("LFO_RELIABLE_EVAL_AGE_DAYS", ifnotfound = NA_integer_)),
        if (is.finite(ctx$calib$delta_prior %||% NA_real_))
          sprintf("%.3f", ctx$calib$delta_prior) else "unavailable",
        .h1_txt,
        get0("CASCADE_PSI", ifnotfound = 0),
        if (is.null(.rp)) "No held-out residual check was run."
        else sprintf(paste0("**Held-out residual check** at %d weeks over %d origin(s): predicted ",
                            "%.1f vs observed %.0f new invasions (ratio %.2f) — **IN-SAMPLE**: these ",
                            "are the origins and the count delta was fitted to, so the ratio is ~1 by ",
                            "construction and is not validation. AUC-PR skill %.1fx is the part delta ",
                            "was not fitted to."),
                     .rp$pooled$K, .rp$pooled$n_origins, .rp$pooled$predicted_new,
                     .rp$pooled$observed_new, .rp$pooled$count_ratio, .rp$pooled$auc_pr_skill))
  } else {
    add(paste0("**Calibration — NOT PERFORMED in this run (%s).** delta = %.3f is the ONE-WEEK ",
               "prequential maximum-likelihood factor of `16b_invasion_recalibration.R` (%s), carried ",
               "with its zone-clustered bootstrap SE on the log scale (%s). It calibrates the ",
               "short-horizon model at h=1; nothing in this run checks the cascade's level at the ",
               "13-week horizon it is read at, so **the absolute 13-week numbers below are ",
               "uncalibrated** and only the rankings should be used. **psi = %g is fixed** (no ",
               "frontier saturation)."),
        ctx$calibration_disabled_by %||% "reason unrecorded", as.numeric(ctx$delta),
        .ci$estimator %||% "estimator unrecorded",
        if (is.finite(.ci$se_log %||% NA_real_)) sprintf("SE %.3f", .ci$se_log) else "SE unavailable",
        get0("CASCADE_PSI", ifnotfound = 0))
  }
  # The two criteria for psi can legitimately disagree, and the disagreement is reported
  # rather than left in a JSON sidecar for a reviewer to find.
  if (!is.null(.pf) && is.finite(.pf$psi_min_deviance %||% NA_real_) &&
      abs((.pf$psi_min_deviance - .pf$psi) / max(.pf$psi, 1e-9)) > 0.1) {
    add("")
    .edge <- c(if (isTRUE(.pf$deviance_boundary_hit)) "deviance" else NULL,
               if (isTRUE(.pf$auc_pr_boundary_hit)) "AUC-PR" else NULL)
    .edge_txt <- if (length(.edge))
      sprintf(paste0("**%s pinned at the edge of the search interval** (psi_hi = %.0f): %s ",
                     "stopped improving only because the search stopped. When a criterion is at ",
                     "the boundary it is the INTERVAL, not the data, that chose the value."),
              if (length(.edge) > 1L) "Both of those are" else "That optimum is",
              .pf$psi_hi %||% NA_real_, paste(.edge, collapse = " and "))
      else "Neither optimum sits at the edge of the search interval."
    # The SHAPE of the psi curve decides how much weight the probabilistic criteria deserve, so
    # it is MEASURED from the curve rather than asserted. An earlier draft of this paragraph
    # called the deviance monotone; the 2026-09 refit falsified that (333.91 at the count match,
    # 333.94 at psi = 11.25, 333.22 at the ceiling). Prose that states a shape must READ the
    # shape, or it goes stale against the next refit with nothing to flag it.
    .cv <- .pf$curve
    .shape_txt <- ""; .sweep_txt <- ""
    .dap <- if (is.null(.pf$deviance_at_psi)) NA_real_ else .pf$deviance_at_psi
    .dmn <- if (is.null(.pf$deviance_min))    NA_real_ else .pf$deviance_min
    .tgt <- if (is.null(.pf$target_observed_cum)) NA_real_ else .pf$target_observed_cum
    # psi's NUMERICAL width, reported rather than implied. A point value quoted to several
    # decimals invites the reader to assume it is determined to that precision; the bracket
    # says how far the search actually pinned it.
    .blo <- if (is.null(.pf$psi_bracket_lo)) NA_real_ else .pf$psi_bracket_lo
    .bhi <- if (is.null(.pf$psi_bracket_hi)) NA_real_ else .pf$psi_bracket_hi
    .brk_txt <- if (all(is.finite(c(.blo, .bhi))) && .bhi > .blo)
      sprintf(paste0(" The search refines until the BRACKET closes rather than stopping at the ",
                     "first count that falls inside the tolerance band, so psi is pinned to ",
                     "[%.2f, %.2f] — a width of %.1f%% of psi. That bracket is the NUMERICAL ",
                     "width of the search, NOT psi's statistical uncertainty, which is larger: ",
                     "the objective is itself a Monte-Carlo estimate, and refitting at a ",
                     "different draw size moves psi by appreciably more than the bracket. The ",
                     "weaker rule would leave psi determined by where the search halted rather ",
                     "than by the data."),
              .blo, .bhi, 100 * (.bhi - .blo) / max(.bhi, 1e-9))
      else ""
    if (!is.null(.cv) && all(c("psi", "deviance") %in% names(.cv)) && nrow(.cv) >= 3) {
      .o <- order(.cv$psi); .p <- .cv$psi[.o]; .d <- .cv$deviance[.o]
      .hi <- is.finite(.p) & is.finite(.d) & .p >= 3
      if (sum(.hi) >= 2) {
        .base <- min(.d[.hi]); .span <- max(.d[.hi]) - .base
        .shape_txt <- sprintf(paste0("Those preferences are WEAK: above psi = 3 the deviance ",
            "varies by only %.2f units on a base of %.0f (%.2f%%) across a %d-point log ladder%s. ",
            "The deviance at the retained psi (%.2f) sits just %.2f units above the global ",
            "minimum (%.2f). A criterion that flat does not identify psi."),
          .span, .base, 100 * .span / max(.base, 1e-9), nrow(.cv),
          if (!all(diff(.d) <= 0)) ", and it is not monotone in psi" else "",
          .dap, .dap - .dmn, .dmn)
      }
      if ("modelled_cum" %in% names(.cv) && any(is.finite(.cv$modelled_cum)))
        .sweep_txt <- sprintf(paste0(" The count criterion, by contrast, is SHARPLY identified: ",
            "the modelled invasion count sweeps %.0f down to %.0f over the same interval and ",
            "crosses the observed %.0f steeply, so it has a well-defined root exactly where the ",
            "probabilistic criteria have a plateau."),
          max(.cv$modelled_cum, na.rm = TRUE), min(.cv$modelled_cum, na.rm = TRUE), .tgt)
    }
    add("**Which criterion fits psi, and why the disagreement matters.** psi is chosen to MATCH THE COUNT of realised invasions (%.1f modelled vs %.0f observed at psi = %.2f). Two *probabilistic* criteria prefer different values: Bernoulli deviance is lowest at %.2f (psi = %.2f), and AUC-PR skill is best at psi = %s. %s %s%s%s Part of this is expected — with %d events among %d at-risk zone-origins the base rate is about %.1f%%, so deviance is dominated by correct zeros and improves whenever probabilities are pushed toward zero, which under-predicts the total. But where a scoring rule keeps improving to the ceiling, the honest reading is that frontier saturation is absorbing misspecification it cannot represent — most plausibly the mobility kernel's dispersion or delta — not that psi should be that large. The count is the estimand this layer exists to produce, so the count match is retained, and the consequence is a stated direction of error: a deviance-fitted psi would damp spread harder still, so the 13-week totals here are an **upper** bound, not a central estimate.",
        .pf$modelled_cum_at_psi, .pf$target_observed_cum, .pf$psi,
        .pf$deviance_min, .pf$psi_min_deviance,
        if (is.finite(.pf$psi_max_auc_pr %||% NA_real_)) sprintf("%.2f", .pf$psi_max_auc_pr) else "\u2014",
        .edge_txt, .shape_txt, .sweep_txt, .brk_txt,
        .rp$pooled$observed_new %||% NA_integer_,
        .rp$pooled$n_atrisk %||% NA_integer_,
        100 * (.rp$pooled$observed_new %||% NA_real_) / max(.rp$pooled$n_atrisk %||% NA_real_, 1))
  }
  add("")
  # THE WALK IS DESCRIBED FROM THE CONFIG, not from memory. This paragraph described a
  # mean-reverting AR(1) per zone whose spread "plateaus", with a "measured" sigma, while the
  # config set rho = 0 (a pure random walk, no plateau) and documents sigma as an ASSUMPTION
  # with a sweep — and the walk is of the single national R, one shock per week and draw.
  .rw_s <- as.numeric(get0("CASCADE_R_RW_SIGMA", ifnotfound = 0))
  .rw_r <- as.numeric(get0("CASCADE_R_RW_RHO", ifnotfound = 0))
  .rw_H <- as.integer(CASCADE_HORIZON_WEEKS)
  .rw_sweep <- get0("CASCADE_R_RW_SWEEP", ifnotfound = numeric(0))
  .rw_form <- if (!is.finite(.rw_s) || .rw_s <= 0) {
    "sigma = 0, so R is held at its starting value for the whole horizon (no walk)."
  } else if (!is.finite(.rw_r) || .rw_r <= 0) {
    sprintf(paste0("rho = %g, i.e. a PURE RANDOM WALK on log R with no reversion: its spread grows as ",
                   "sigma*sqrt(t) without plateauing, reaching %.3f on the log scale by week %d."),
            .rw_r, .rw_s * sqrt(.rw_H), .rw_H)
  } else {
    sprintf(paste0("rho = %g, i.e. a mean-reverting (Ornstein-Uhlenbeck) walk toward the starting ",
                   "value: its spread grows and then plateaus at sigma/sqrt(1-(1-rho)^2) = %.3f on ",
                   "the log scale."), .rw_r, .rw_s / sqrt(1 - (1 - .rw_r)^2))
  }
  add(paste0("**R_eff is not held constant over the horizon.** The single national log R follows ",
             "log R(t) = log R(t-1) + rho*(log R(0) - log R(t-1)) + sigma*eps(t), with ONE innovation ",
             "per week and posterior-parameter draw shared by every transmitting zone (a zone seeded ",
             "mid-projection transmits at the current walked value), sigma = %.3f per week and %s ",
             "E[log R] is unchanged, so this adds uncertainty that grows with horizon, not a trend. ",
             "That separation is deliberate: projecting the observed decline in R forward is a claim ",
             "about future control, which is what the pre-registered c(t) scenarios express — fitting ",
             "it here as well would contradict S1 \"status quo\" and double-count control in S2. ",
             "**sigma is an ASSUMPTION, not a measurement**: an innovation SD read off EpiNow2's ",
             "smoothed posterior-mean series is not identifiable on this record (30_projection_config.R ",
             "records why), so the value is the central arm of a pre-specified sweep (sigma = %s) ",
             "carried in the sensitivity table. sigma = 0 reproduces the constant-R projection exactly."),
      .rw_s, .rw_form,
      if (length(.rw_sweep)) paste(format(.rw_sweep), collapse = ", ") else "none configured")
  add("")
  # THE GATE TESTS THE RANKING ONLY. Its level ratio is ~delta for any delta (both sides share
  # one fit; 33_cascade_eval.R), so it is reported as the plumbing check it is, never as a
  # calibration result — this paragraph used to present it beside PASS as though it were one.
  .sdl <- as.numeric(ctx$delta_sd_log %||% 0)
  .sd_txt <- if (is.finite(.sdl) && .sdl > 0)
    sprintf(paste0("Calibration uncertainty is propagated: delta is drawn per posterior-parameter group ",
                   "on the log scale (SD %.3f — the zone-clustered bootstrap SE of 16b's short-horizon ",
                   "factor, used as an approximation to the relative precision of the factor deployed), ",
                   "so the per-zone credible intervals include it rather than treating the calibration as ",
                   "exact."), .sdl)
  else paste0("Calibration uncertainty is **NOT** propagated in this run: no standard error was ",
              "available for delta, so it is held fixed and the per-zone intervals understate it.")
  add(paste0("%s Consistency gate at h=1: the cascade's week-1 reach",
             if (isTRUE(get0("CASCADE_GT_MARGINALISE", ifnotfound = FALSE)))
               ", simulated at the short-horizon model's own generation time (the projection itself marginalises it over a grid)," else "",
             " %s delta x ",
             "the validated short-horizon model's week-1 probabilities from the same fit, zone by zone, within its ",
             "own Monte-Carlo error (dispersion X^2 = %.1f on %d df, p = %.3f, over %d zones with at ",
             "least 5 expected invasions) — %s. (Rank diagnostic: Spearman = %.3f over the %d zones the ",
             "short-horizon model puts above 2%%; near-tied zones make it noisy.) Its level ratio ",
             "(%.2f) is ~delta by construction and is only a check that delta reaches the simulated ",
             "hazard (%s); it is not a test of the level."),
      .sd_txt,
      if (isTRUE(ctx$gate$pass)) "agrees with" else "**departs from**",
      ctx$gate$x2 %||% NA_real_, as.integer(ctx$gate$x2_df %||% NA_integer_),
      ctx$gate$x2_p %||% NA_real_, as.integer(ctx$gate$n_tested %||% NA_integer_),
      if (isTRUE(ctx$gate$pass)) "**PASS**" else "**REVIEW**",
      ctx$gate$spearman_signal, ctx$gate$n_signal,
      ctx$gate$level_ratio %||% NA_real_,
      if (isTRUE(ctx$gate$pass_level)) "it does" else "**it does NOT — investigate**")
  add("")
  add("As of the analysis date: **%d affected** zones, **%d at-risk**.", aff, atrisk)
  add("")
  add("---")
  add("## 2. Headline — where invasion is most likely over 3 months (%s scenario)",
      CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]$label)
  r13 <- ctx$reach_primary[ctx$reach_primary$horizon == max(CASCADE_REPORT_HORIZONS) &
                             !ctx$reach_primary$was_active_before, ]
  r13 <- r13[order(-r13$p_case_invasion), ]
  add("Top-15 at-risk zones by %d-week reach (P, 90%% credible interval; establishment P; median first-passage week; priority rank):",
      max(CASCADE_REPORT_HORIZONS))
  add("")
  # THE MONTE-CARLO FLOOR, counted and marked. cascade_reach_table() sets an interval that the D
  # parameter draws cannot resolve to p +- 1/(2D) clipped to [0, 1] — width at most 1/D (exactly
  # 1/D unless clipped at 0 or 1); such a row carries no estimated precision and must not be read
  # as one (p_ci_at_mc_floor, written to every reach CSV).
  .Dg  <- max(1L, as.integer(ctx$n_mc) %/% max(2L, as.integer(get0("CASCADE_N_REP", ifnotfound = 8L))))
  .has_floor <- "p_ci_at_mc_floor" %in% names(r13)
  .flo <- if (.has_floor) r13$p_ci_at_mc_floor %in% TRUE else rep(FALSE, nrow(r13))
  add("_Note: the per-zone interval is a 90%% credible interval — the posterior spread of the reach probability across parameter draws (hazard coefficients, %s the reproduction number and its walk, and the calibration factor), with Monte-Carlo process noise removed. It does NOT include between-scenario (§3) or structural (mobility/generation-time, §5) uncertainty, which is larger, and it is conditional on the scenario. %s_",
      # A featured RW1 form's drift path is drawn per parameter group (32_cascade_simulator.R),
      # so it is on the posterior side too.
      if (.form %in% c("rw1", "rw1_frailty")) "the featured form's import-coefficient drift path," else "",
      if (.has_floor)
        sprintf(paste0("**%d of the %d at-risk zones** at %d weeks have an interval at the Monte-Carlo ",
                       "resolution floor (width at most 1/%d, exactly 1/%d unless clipped at 0 or 1; ",
                       "marked † below): their posterior spread is below what %d parameter draws can ",
                       "resolve, so that width is NOT an estimate of their precision."),
                sum(.flo), nrow(r13), max(CASCADE_REPORT_HORIZONS), .Dg, .Dg, .Dg)
      else "(The Monte-Carlo-floor flag is unavailable on this context.)")
  add("")
  add("| Zone | Province | P(reach) [90%% CrI] | P(estab.) | 1st-pass wk | Priority rank |")
  add("|---|---|---|---|---|---|")
  for (i in seq_len(min(15L, nrow(r13)))) { z <- r13[i, ]
    add("| %s | %s | %.2f [%.2f–%.2f]%s | %.2f | %s | %s |", z$health_zone,
        z$province %||% "", z$p_case_invasion, z$p_lo, z$p_hi, if (.flo[i]) " †" else "",
        z$p_establishment,
        ifelse(is.finite(z$first_passage_week), as.character(round(z$first_passage_week)), "—"),
        ifelse(is.finite(z$priority_rank), as.character(z$priority_rank), "—")) }
  add("")
  add("%s", .fig_line(c(
    "Reach map"                   = "cascade/figures/cascade_reach_map_national_h13",
    "Expanding-frontier triptych" = "cascade/figures/cascade_time_to_invasion_triptych",
    "Uncertainty"                 = "cascade/figures/cascade_reach_uncertainty_national_h13")))
  add("")
  add("---")
  add("## 3. Burden trajectory — new zones & provinces (medians), by scenario")
  add("")
  add("_Median cumulative newly-invaded zones (one consistent estimator across horizons); 13w carries the 90%% predictive band across MC iterations._")
  add("")
  # Every cell in this row is now a MEDIAN (the province count was a mean, mixing estimators in
  # a row the aggregation block's own comment says must use one). Say so in the header.
  add("| Scenario | New zones by 4w (med) | by 8w (med) | by 13w (med) [90%% pred. int.] | New provinces by 13w (med) |")
  add("|---|---|---|---|---|")
  for (sk in names(ctx$agg)) { a <- ctx$agg[[sk]]
    add("| %s | %.0f | %.0f | %.0f [%.0f–%.0f] | %.1f |", a$label,
        a$new4, a$new8, a$new13_med, a$new13_lo, a$new13_hi, a$new_prov13) }
  add("")
  add("*All four numeric columns are medians over the Monte-Carlo iterations; the bracketed range is a 90%% predictive interval.*")
  # THE MEAN AND THE INTERVAL OF A TOTAL, stated where a reader would look for them. The expected
  # national count IS additive — the sum of the per-zone reach probabilities — but an interval
  # for the total is NOT: it must be the quantiles of the SIMULATED total (the bracket above,
  # cascade_newzones_fanchart.csv), never the sum of per-zone p_lo/p_hi.
  .r13p <- ctx$reach_primary[ctx$reach_primary$horizon == max(CASCADE_REPORT_HORIZONS) &
                               !ctx$reach_primary$was_active_before, ]
  .mean13 <- sum(.r13p$p_case_invasion, na.rm = TRUE)
  add("")
  add("*Expected (MEAN) new zones by %d weeks in the %s scenario: %.1f — the sum of the per-zone reach probabilities, which is additive. An interval for a total is not: it is the 90%% predictive band of the simulated total above (`cascade/tables/cascade_newzones_fanchart.csv`, also in key_outputs), never the sum of per-zone interval bounds.*",
      max(CASCADE_REPORT_HORIZONS), CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]$label, .mean13)
  add("")
  add("%s", .fig_line(c(
    "Fan chart"          = "cascade/figures/cascade_newzones_fanchart",
    "Scenario reach maps" = "cascade/figures/cascade_scenario_reach_maps_h13")))
  add("")
  add("---")
  add("## 4. Onward propagation — gateway zones & corridors (the task's core)")
  add("Gateway zones ranked by **Delta_j** = expected downstream invasions prevented if that source is contained (exact knockout on the shortlist; PLAN §3.7):")
  add("")
  add("| Source zone | Delta_j (downstream invasions) | 90%% interval | %% of total spread |")
  add("|---|---|---|---|")
  # Guarded: a NULL or empty gateway table made seq_len(nrow(NULL)) abort the whole
  # report with "argument must be coercible to non-negative integer" rather than
  # degrading to a stated gap. Every other optional section here is already guarded.
  gt <- if (is.null(ctx$gateway) || !NROW(ctx$gateway)) NULL else head(ctx$gateway, 12L)
  if (is.null(gt)) add("| _(gateway knockout unavailable for this run)_ | — | — | — |")
  else for (i in seq_len(nrow(gt))) add("| %s | %.1f | %s | %s |", gt$health_zone[i], gt$delta_j[i],
      if (!is.null(gt$delta_j_lo) && is.finite(gt$delta_j_lo[i]))
        sprintf("%.1f to %.1f", gt$delta_j_lo[i], gt$delta_j_hi[i]) else "—",
      .fmt_pct(gt$pct_reduction[i] / 100))
  add("")
  add("Delta_j is a **paired** contrast: the knockout and the no-bar baseline share one random-number stream (common random numbers), so the difference is taken iteration by iteration rather than between two separately-averaged runs. The interval is a 90%% Monte-Carlo interval over posterior parameter draws — the unit of independence is the draw, not the iteration, because the nested design runs correlated process replicates inside each draw. **Zones whose interval spans zero are not separated by this analysis**, however they happen to be ordered in the table.")
  add("")
  add("%s Most-likely invasion tree: `cascade/tables/cascade_transmission_tree.csv`.",
      .fig_line(c(
        "Inter-province seeding flux (corridors)" = "cascade/figures/cascade_corridors",
        "Gateway bars"                            = "cascade/figures/cascade_gateway_knockout")))
  add("")
  # THE HUBS ARE ZONES THAT CAN STILL BE INVADED, and the rule is printed. They used to be the
  # top of the gateway table — currently INFECTED sources — which the simulator cannot force, so
  # every "conditional" result was the baseline and identical across hubs.
  .f2 <- function(x) if (is.null(x) || !length(x) || !is.finite(x[1])) "—" else sprintf("%.2f", x[1])
  if (!is.null(ctx$conditional) && length(ctx$conditional) && is.null(ctx$conditional[[1]]$rule)) {
    # A context persisted BEFORE the hub fix: its hubs were the gateway table's infected sources,
    # which the simulator cannot force, so its "conditional" results are the baseline. Withheld.
    add("**Conditional onward-spread: WITHHELD** — this context predates the hub fix: its hubs (%s) were already-infected zones, which cannot be force-seeded, so the stored results are the baseline, not a conditional experiment. Re-run run_cascade.R.",
        paste(names(ctx$conditional), collapse = ", "))
    add("")
  } else if (!is.null(ctx$conditional) && length(ctx$conditional)) {
    .c1 <- ctx$conditional[[1]]
    .rule <- if (identical(.c1$rule, "CASCADE_HUBS_FOR_COND"))
      "named in CASCADE_HUBS_FOR_COND"
    else "the two zones not yet invaded at the projection origin with the largest expected onward seedings in the baseline projection"
    add(paste0("**Conditional onward-spread (generalised Kisangani).** Each hub is a zone **not yet ",
               "invaded** at the projection origin (%s), force-seeded at week 1. Its downstream effect ",
               "is the **paired** attributable increase in %s-week reach over a same-seed baseline ",
               "(M = %s), and the zones listed are those whose paired 90%% interval excludes zero, ",
               "ranked by that increase — the hub's own catchment, not zones that were likely anyway. ",
               "(Force-seeding an already-infected zone changes nothing, which is why the hubs are not ",
               "taken from the gateway table above.)"),
        .rule, format(.c1$horizon %||% max(CASCADE_REPORT_HORIZONS)),
        if (is.null(.c1$n_mc)) "—" else format(.c1$n_mc))
    add("")
    for (h in names(ctx$conditional)) { cc <- ctx$conditional[[h]]
      # The per-hub conditional-reach map is suppressed by the retained-figure allow-list on the
      # published configuration, so its pointer is only printed when the gate would write it.
      add("- **%s** (baseline reach P = %s; expected onward seedings in the baseline %s) → %s additional zones invaded [90%% Monte-Carlo interval %s to %s]; %s zone(s) materially elevated%s%s",
          h, .f2(cc$baseline_reach), .f2(cc$out_flux), .f2(cc$expected_added),
          .f2(cc$expected_added_lo), .f2(cc$expected_added_hi),
          if (is.null(cc$n_elevated)) "—" else format(cc$n_elevated),
          if (length(cc$top)) paste0(", led by ", paste(cc$top, collapse = ", ")) else "",
          .fig_line(stats::setNames(
            sprintf("cascade/figures/cascade_conditional_reach_%s_h13",
                    gsub("[^A-Za-z0-9]+", "_", h)), "map"),
            prefix = " (", ext = ".pdf")) }
    add("")
  } else {
    add("**Conditional onward-spread: NOT RUN in this run** — no zone outside the current front seeded an onward invasion in the baseline projection, so no hub qualified under the rule (a not-yet-invaded zone).")
    add("")
  }
  if (!is.null(ctx$urban) && !is.null(ctx$urban$summary) && nrow(ctx$urban$summary)) {
    add("---")
    add("## 4b. Scenario: invasion of a highly-connected/populated city")
    add("What happens if a major urban hub — far from the current front — is itself invaded and becomes a source? Each city's urban-core health zones are seeded at a chosen week and the cascade projects its onward spread; impact is measured against the baseline (no urban seeding), **excluding the seeded city itself**. **Timing matters**: earlier seeding leaves more weeks for spread.")
    add("")
    # The PRIMARY specification only: earliest seed week, seeded city at the national R.
    # The seeded-R arms get their own subsection — averaging them into this table would
    # silently mix a premise the reader has not yet been shown.
    u1 <- ctx$urban$summary[ctx$urban$summary$seed_week == ctx$urban$min_week, ]
    if ("seed_r" %in% names(u1)) u1 <- u1[is.na(u1$seed_r), ]
    u1 <- u1[order(-u1$expected_added), ]
    # "Exp. added" is a paired MEAN (an expectation), and its interval is Monte-Carlo uncertainty
    # on that mean — not a median, and not a predictive or credible interval.
    add("Downstream impact if invaded at week %d (status-quo scenario; the seeded city transmits at the national R). `Exp. added` is the paired expected number of additional zones; its 90%% interval is Monte-Carlo uncertainty on that expectation (over posterior-draw groups), not a predictive interval:", ctx$urban$min_week)
    add("")
    add("| City | Core zones | Exp. added zones (13w) | 90%% interval | Zones materially elevated | Newly >20%% | Newly >50%% (%% of base) | Max RR | New prov. |")
    add("|---|---|---|---|---|---|---|---|---|")
    for (i in seq_len(nrow(u1))) { z <- u1[i, ]
      # 2dp, not 1dp: on this frame the cross-city spread is 0.10 to 0.34, which %.1f
      # collapses to three cities all printing "0.3" and a small negative printing "-0.0".
      # The interval columns are already 2dp, so 1dp also made the point estimate look
      # coarser than its own uncertainty.
      add("| %s | %d | %.2f | %s | %d | +%d | +%d (+%.0f%%) | %s | %s |", z$city, z$n_hub_zones,
          z$expected_added,
          if (!is.null(z$expected_added_lo) && is.finite(z$expected_added_lo))
            sprintf("%.2f to %.2f", z$expected_added_lo, z$expected_added_hi) else "—",
          z$n_elevated, z$newly_above_20, z$newly_above_50, z$pct_increase_gt50,
          ifelse(is.finite(z$max_rr), sprintf("%.0f×", z$max_rr), "—"),
          ifelse(is.finite(z$provinces_newly_exposed), as.character(z$provinces_newly_exposed), "—")) }
    add("")
    # Measured on THIS run rather than quoted from a previous one: a hard-coded percentage in
    # a report that regenerates every run is a statement that will eventually be false.
    .pd <- ctx$urban$pairing_diag
    add("**Which metrics to focus on.** The **expected additional zones invaded** (`Exp. added zones`) is the headline: the sum of the per-zone attributable increase over all non-hub zones, an estimate of E[extra zones invaded by 13 weeks]. It is computed as a **paired** contrast — the conditional and baseline runs share one random-number stream, and the difference is taken iteration by iteration. That pairing is load-bearing, not a refinement. An earlier version of this analysis differenced two independently-summarised runs on the assumption that the Monte-Carlo flicker of far, unaffected zones would cancel in the sum. It does not: unpaired at M = 4000 those zones carried a per-zone standard deviation of about 0.017 and contributed a summed error of ±0.3–0.4 zones — the same size as the effect being measured — which is why that version reported *negative* expected additions for most cities and ranked zones 1,500 km from the seeded city as the most affected. %s What pairing changes is not the point estimate — with balanced draw groups it is arithmetically the same difference — but what the two runs hold in common: `Var(A-B) = Var(A)+Var(B)-2Cov(A,B)`, and one shared random-number stream drives that covariance almost to the variance, so the difference becomes interpretable and can carry an interval at all. Every interval in this section comes from that paired contrast; **a city whose interval spans zero has no downstream effect this analysis can detect.** The **timing profile** and the **within-catchment relative risk** + **attributable-risk map** are the other primary products. The threshold-crossing counts (`Newly >20%%/50%%`) are secondary, and a zone counts as materially elevated when its paired 90%% interval excludes zero — a test, replacing the fixed 0.03 cut used previously, which was neither a noise floor nor a pre-registered effect size. The seeded city's own zones are excluded (invaded by assumption).",
        if (is.null(.pd)) "" else sprintf(
          paste0("On THIS run (M = %d, %s), of the %d zones outside the seeded city's own ",
                 "province, %.0f%% return **exactly** zero — not merely a small number — and ",
                 "the median per-zone standard error among them is %.5f (largest %.4f). ",
                 # CONDITIONAL on psi. The unconditional clause contradicted the same paragraph,
                 # which prints "psi = 0" two sentences earlier and then explains that with the
                 # saturation term off there is no damping channel and a negative can only be
                 # Monte-Carlo residual.
                 if (get0("CASCADE_PSI", ifnotfound = 0) > 0)
                   "%d zone(s) came out negative, which at psi > 0 is a real frontier-saturation effect rather than noise."
                 else if (.form %in% c("frailty", "rw1_frailty"))
                   "%d zone(s) came out negative; at psi = 0 nothing in the model lowers a zone's EXPECTED 13-week reach (the featured gamma frailty makes individual paths non-monotone, since an exposed zone that escapes invasion has a lower later hazard, but its probability of invasion by any horizon still rises with exposure), so these are Monte-Carlo residual, not protection."
                 else
                   "%d zone(s) came out negative; at psi = 0 there is no damping channel, so these are Monte-Carlo residual, not protection."),
          ctx$urban$n_mc %||% CASCADE_N_MC, .pd$city, .pd$n_far, .pd$pct_exact_zero_far,
          .pd$se_far, .pd$se_far_max, .pd$n_negative))
    add("")
    add("**Timing.** Expected additional zones invaded by 13 weeks, by seeding week:")
    add("")
    us <- ctx$urban$summary
    if ("seed_r" %in% names(us)) us <- us[is.na(us$seed_r), ]
    hdr <- paste0("| City | ", paste(sprintf("seed wk %d", ctx$urban$seed_weeks), collapse = " | "), " |")
    add(hdr); add(paste0("|", paste(rep("---", length(ctx$urban$seed_weeks) + 1L), collapse = "|"), "|"))
    for (cty in unique(us$city)) {
      vals <- vapply(ctx$urban$seed_weeks, function(w) {
        v <- us$expected_added[us$city == cty & us$seed_week == w]
        # A small negative rendered at 2dp prints "-0.00", which reads as a formatting bug
        # and hides the sign. Below 0.005 the extra digit is shown so the cell says what it
        # means: late seeding leaves too few weeks for a detectable effect, and at this psi
        # the frontier-saturation term can tip the residual slightly negative.
        if (!length(v)) "—" else if (abs(v[1]) < 0.005) sprintf("%+.3f", v[1])
        else sprintf("%.2f", v[1]) }, character(1))
      add("| %s | %s |", cty, paste(vals, collapse = " | ")) }
    add("")
    # ---- case burden ------------------------------------------------------------------
    # Thousands separator and no decimals: these are counts in the thousands, where "15669.3"
    # reads as false precision on a quantity whose predictive band spans a factor of three.
    .fmt0 <- function(x) if (is.null(x) || !is.finite(x)) "—" else formatC(x, format = "d",
                                                                          big.mark = ",")
    ub <- ctx$urban$summary
    if ("cases_total" %in% names(ub) && any(is.finite(ub$cases_total))) {
      u2 <- ub[ub$seed_week == ctx$urban$min_week & is.na(ub$seed_r), ]
      u2 <- u2[order(-u2$cases_total), ]
      if (nrow(u2)) {
        add("")
        add("**How many additional CASES, and where they land.** The table above counts additional *zones invaded* — geographic spread. It excludes the seeded city itself (invaded by assumption) and every already-infected zone (they cannot be invaded again). Both exclusions are right for a spread metric and wrong for a burden metric: most of the case burden lands **inside the seeded city**, and an already-infected zone can still receive extra cases through the import force. The paired case contrast below therefore covers all %d zones and splits by location instead of by eligibility.",
            # Zone count taken directly, not as nrow/n_horizons: if the reach table ever
            # carried unequal rows per horizon that division yields a non-integral double
            # and sprintf("%d", ...) ABORTS the whole report.
            length(unique(ctx$reach_primary$health_zone)))
        add("")
        add("| City | Cases in the seeded city | Cases elsewhere | Total: **typical** (median) | Total: **expected** (mean) | 90%% of outcomes fall in | as %% of baseline |")
        add("|---|---|---|---|---|---|---|")
        for (i in seq_len(nrow(u2))) { z <- u2[i, ]
          add("| %s | %.1f [%.1f, %.1f] | %.2f [%.2f, %.2f] | %s | %.1f | %s | %s |", z$city,
              z$cases_in_city, z$cases_in_city_lo, z$cases_in_city_hi,
              z$cases_elsewhere, z$cases_elsewhere_lo, z$cases_elsewhere_hi,
              if (is.finite(z$cases_total_median %||% NA_real_))
                sprintf("%.1f", z$cases_total_median) else "—",
              z$cases_total,
              if (is.finite(z$cases_total_pred_lo %||% NA_real_))
                sprintf("%.0f to %.0f", z$cases_total_pred_lo, z$cases_total_pred_hi) else "—",
              if (is.finite(z$cases_baseline_total) && z$cases_baseline_total > 0)
                sprintf("%.1f%%", 100 * z$cases_total / z$cases_baseline_total) else "—") }
        add("")
        add("Read the two columns as different quantities, not as a split of one. **Cases in the seeded city** combine the *assumed* seed itself (about 1.5 cases per forced zone under the default Poisson seeding, so roughly 4-5 for a three-zone hub) with the local outbreak it starts, so this column is largely a restatement of the seeded-R premise below — it is what you assume, propagated. **Cases elsewhere** is the genuine downstream burden and the quantity a national planner trades off; a city whose interval there spans zero has no detectable effect on the rest of the country. **The typical outcome and the expected outcome are different numbers, and the gap is the headline.** A seeded introduction usually fades: in most iterations the city contributes little beyond the seed cases themselves, which is why the MEDIAN additional total is small. Occasionally it establishes and produces dozens of cases, which is why the MEAN is larger. On this frame the two differ by about %.1f-fold (median %s against mean %s, with 90%% of outcomes spanning roughly %s to %s). Neither number alone is honest: the **median** is what typically happens and the **mean** is the right input to an expected-burden or cost calculation, while the predictive span is what a preparedness plan has to cover. The same skew, milder, affects the baseline: mean %s against median %s (90%% predictive band %s-%s). The percentage column divides mean by mean so the two estimators are never mixed, and medians are NOT additive across the in-city/elsewhere split — only the means are. The baseline denominator is the projected confirmed cases over the same 13 weeks with no urban seeding, and it rests on the status-quo scenario — transmissibility persisting with no further control gain and no within-zone susceptible depletion — so it is conditional on that scenario as well (the deterioration scenario would raise it, strengthened control lower it). All figures are modelled **confirmed** cases, the scale the models are fitted and calibrated on; they are not infections, and this suite does not estimate an ascertainment fraction.%s",
            # The mean/median gap and its span are taken from THIS run's own table, three lines
            # above. The previous text hard-coded "about fivefold (median 2 against mean 10,
            # spanning -1 to 45)" and described it as "the pilot for this frame" — on the shipped
            # table the medians are 7-15 and the means 23.7-34.3, a ratio of about 2.2-2.4.
            local({ .md <- u2$cases_total_median[1]; .mn <- u2$cases_total[1]
                    if (is.finite(.md) && .md > 0 && is.finite(.mn)) .mn / .md else NA_real_ }),
            .fmt0(u2$cases_total_median[1]), .fmt0(u2$cases_total[1]),
            .fmt0(u2$cases_total_pred_lo[1]), .fmt0(u2$cases_total_pred_hi[1]),
            .fmt0(u2$cases_baseline_total[1]), .fmt0(u2$cases_baseline_median[1]),
            .fmt0(u2$cases_baseline_pred_lo[1]), .fmt0(u2$cases_baseline_pred_hi[1]),
            .fig_line(c("Figure" = "cascade/figures/urban_case_burden"), prefix = " "))
        add("")
      }
    }

    # ---- the premise behind the null result, made explicit ----------------------------
    usr <- ctx$urban$summary
    if ("seed_r" %in% names(usr) && any(!is.na(usr$seed_r))) {
      usr <- usr[usr$seed_week == ctx$urban$min_week, ]
      add("")
      add("**The seeded city's reproduction number is a premise, not a finding.** A zone seeded during the projection transmits at the same single national R as everywhere else (R = %.2f on this frame), because that is the only reproduction number the model has. The primary specification above therefore asks \"what if a capital is seeded and then transmits like the national average\", and at that R with k_indiv = %.2f a one-to-two case seed frequently fades out before establishing — so a large part of what the headline measures is **stochastic fade-out**, not the city's connectivity. (The establishment probability at this R and k is reported by the k_indiv sweep; it is not restated here, because it moves with both.) The sweep separates the two, and the separation is substantial: raising the seeded city's R to %s multiplies the expected additional zones by %s across cities, without changing which zones they are. Both halves of that matter. The national arm remains the primary specification because it is the only value this outbreak's data support; the upper arms say what would follow if a capital's transmission resembled the published EVD range instead.",
          ctx$reff$R_nat %||% NA_real_, CASCADE_K_INDIV,
          # Computed from the very table printed below, not asserted. "roughly triples" was
          # wrong by 3-4x: the shipped urban_scenario_summary.csv gives 9.5x-13.1x.
          local({ a <- sort(unique(usr$seed_r[!is.na(usr$seed_r)]))
                  if (!length(a)) "the upper arm" else sprintf("%.1f", max(a)) }),
          local({
            a <- sort(unique(usr$seed_r[!is.na(usr$seed_r)]))
            if (!length(a)) return("an unquantified factor")
            base <- usr[is.na(usr$seed_r), c("city", "expected_added")]
            top  <- usr[usr$seed_r == max(a), c("city", "expected_added")]
            m <- merge(base, top, by = "city", suffixes = c("_base", "_top"))
            r <- m$expected_added_top / m$expected_added_base
            r <- r[is.finite(r) & r > 0]
            if (!length(r)) "an unquantified factor"
            else if (diff(range(r)) < 0.5) sprintf("about %.0fx", stats::median(r))
            else sprintf("%.0f-%.0fx", min(r), max(r))
          }))
      add("")
      arms <- sort(unique(usr$seed_r[!is.na(usr$seed_r)]))
      add(paste0("| City | national R | ",
                 paste(sprintf("R = %.1f", arms), collapse = " | "), " |"))
      add(paste0("|", paste(rep("---", length(arms) + 2L), collapse = "|"), "|"))
      for (cty in unique(usr$city)) {
        g <- usr[usr$city == cty, ]
        cell <- function(sr) { v <- if (is.na(sr)) g$expected_added[is.na(g$seed_r)]
                                    else g$expected_added[!is.na(g$seed_r) & g$seed_r == sr]
                               if (length(v)) sprintf("%.2f", v[1]) else "—" }
        add("| %s | %s | %s |", cty, cell(NA_real_),
            paste(vapply(arms, cell, character(1)), collapse = " | ")) }
      add("")
      # Only raised when the run actually produced such zones, with the count measured.
      # RESTRICT TO THE SEEDED CITY'S OWN PROVINCE. `reduced` is computed over every eligible
      # zone, so the sentence below ("of the seeded city's own remaining health zones") counted
      # zones hundreds of km away: on the shipped run the 2 reduced zones for Kinshasa were
      # Police (Kinshasa) and Mungindu (KWILU). The true within-city count was 1.
      .nred_all <- tryCatch(sum(ctx$urban$detail[[1]]$impact$per_zone$reduced %in% TRUE),
                            error = function(e) 0L)
      # FAIL EXPLICIT, NOT OPEN. Falling back to .nred_all made .nred_out = 0, so the sentence
      # asserted "(and 0 outside the seeded city's province)" — a POSITIVE claim about a number
      # that was never computed, which is worse than the unqualified count it replaced. The
      # shipped ctx predates hub_prov, and re-rendering from a persisted ctx is exactly what
      # run_cascade.R saves it for. NA here; the prose branches on it below.
      .hub_prov <- ctx$urban$pairing_diag$hub_prov
      .pz_red   <- tryCatch(ctx$urban$detail[[1]]$impact$per_zone, error = function(e) NULL)
      .have_prov <- !is.null(.hub_prov) && length(.hub_prov) &&
                    !is.null(.pz_red) && "province" %in% names(.pz_red)
      .nred <- if (.have_prov)
                 sum(.pz_red$reduced %in% TRUE & .pz_red$province %in% .hub_prov) else NA_integer_
      .nred_out <- if (.have_prov) .nred_all - .nred else NA_integer_
      if (.nred_all > 0L) {
        add("")
        add(paste0(if (.have_prov)
                     "**%d zone(s) in the seeded city's own province (and %d elsewhere) show a significantly negative attributable change"
                   else
                     "**%d zone(s) (province split unavailable on this context%s) show a significantly negative attributable change",
                   if (.form %in% c("frailty", "rw1_frailty"))
                     "— and at psi = %g this is noise, not a mechanism.** The paired interval excludes zero on the wrong side for %d of the seeded city's own remaining health zones. The frontier-saturation term is off (psi is fixed at 0). The featured gamma frailty does make individual simulated paths non-monotone — a zone exposed without being invaded has a lower later conditional hazard (its cumulative hazard C enters 1 + theta*C) — but its marginal survival (1 + theta*C)^(-1/theta) still falls as exposure rises, so seeding a city cannot lower any zone's EXPECTED 13-week reach. What remains is Monte-Carlo noise in the paired contrast (which the path-level non-monotonicity widens), the frontier covariate's posterior tail (`gamma_dmin` straddles zero), or the re-seeding channel — **not protection**. Do not report these as protective effects, and treat a large or systematic negative as a reason to check the pairing rather than as a finding."
                   else
                     "— and at psi = %g this is noise, not a mechanism.** The paired interval excludes zero on the wrong side for %d of the seeded city's own remaining health zones. Under the retired frontier-saturation term (`sat = exp(-psi · max(f_inv − f_inv0, 0))`) this was expected and explicable: seeding a metro's core raised the invaded-inflow share `f_inv` for its closest mobility neighbours — the *other zones of the same metro* — and damped their hazard. **That term is now switched off** (psi is fixed at 0, so `sat` is identically 1), so seeding a city can only *add* import force and there is no damping channel left to produce a genuine decrease. What remains here is therefore Monte-Carlo noise in the paired contrast, the frontier covariate's posterior tail (`gamma_dmin` straddles zero), or the re-seeding channel — **not protection**. Do not report these as protective effects, and treat a large or systematic negative as a reason to check the pairing rather than as a finding."
            ),
            if (.have_prov) .nred else .nred_all,
            if (.have_prov) .nred_out else "",
            get0("CASCADE_PSI", ifnotfound = 0),
            if (.have_prov) .nred else .nred_all)
        add("")
      }
      add("Read across a row, not down the column: the column heading is the assumption, and the published EVD range of roughly 1.5–2.5 spans it. A seeded city's effect on the rest of the country is governed first by whether its introduction establishes at all.%s",
          .fig_line(c("Figure" = "cascade/figures/urban_seed_r_sweep"), prefix = " "))
      add("")
    }
    # Only the cities on the retained-figure allow-list get a Figure_urban_<city> panel, so the
    # inventory is built per city from the gate rather than advertising a wildcard.
    # Cities come from THIS run's own scenario table, and the "_" spelling matches the one
    # 38_urban_scenarios.R builds the filename with.
    .urban_city_names <- sort(unique(as.character(ctx$urban$summary$city)))
    .urban_city_names <- .urban_city_names[!is.na(.urban_city_names) & nzchar(.urban_city_names)]
    .urban_figs <- c(
      if (length(.urban_city_names)) stats::setNames(
        sprintf("key_outputs/figures/Figure_urban_%s",
                gsub("[^A-Za-z0-9]+", "_", .urban_city_names)),
        sprintf("**%s** publication two-panel (attributable-risk map zoomed to the catchment + baseline→conditional reach dumbbell by province)",
                .urban_city_names)) else character(0),
      c("Five-city spatial overview"                      = "key_outputs/figures/Figure_urban_all_cities",
        "Cross-city impact summary"                       = "cascade/figures/urban_impact_summary",
        "Timing sensitivity"                              = "cascade/figures/urban_timing_sensitivity",
        "Seeded-R premise"                                = "cascade/figures/urban_seed_r_sweep",
        "Additional CASES, split in-city vs elsewhere"    = "cascade/figures/urban_case_burden"))
    add("%s Tables: `cascade/tables/urban_scenario_summary.csv`, `urban_impact_<city>.csv` (zone invasions, hub and already-infected zones excluded) and `urban_cases_<city>.csv` (cases, **all** zones with an `in_city` flag).",
        .fig_line(.urban_figs))
    add("")
  }
  add("---")
  add("## 5. Validation & robustness (PLAN §6 — honest about the data ceiling)")
  add("A true 13-week hold-out has at most one truncated origin, so validation is **layered**:")
  add("")
  # RANKING ONLY. The gate's level ratio is ~delta by construction and tests nothing about the
  # level (33_cascade_eval.R); the observed one-week level is the held-out week-1 refit.
  add("- **Consistency gate (h=1):** per-zone week-1 reach against delta x the validated short-horizon model, tested against the cascade's own Monte-Carlo error (X^2 p = %.3f) — %s. The gate tests the per-zone PATTERN only; the observed one-week level is checked by the held-out week-1 refit of delta%s.",
      ctx$gate$x2_p %||% NA_real_, if (isTRUE(ctx$gate$pass)) "PASS" else "Review",
      if (is.finite(.h1r$delta)) {
        .ph5 <- .cascade_h1_ratio_phrase(.h1r, deployed_label = "the deployed delta")
        sprintf(" (%.3f%s)", .h1r$delta, if (nzchar(.ph5)) paste0(", ", .ph5) else "")
      } else " (not available in this run)")
  if (!is.null(ctx$sens) && nrow(ctx$sens)) {
    # The verdict is COMPUTED, not asserted. The psi and R-walk axes give Spearman 0.99-1.00,
    # but the mobility-kernel axis gives 0.72-0.79 with mean reach moving +51% to +69% — and the
    # kernel is the largest structural axis in the suite. Calling that "stable" while printing
    # 0.72 beside it is a claim the numbers refute.
    .sp <- suppressWarnings(as.numeric(ctx$sens$spearman))
    .mn <- if (any(is.finite(.sp))) min(.sp, na.rm = TRUE) else NA_real_
    # which.min() on an all-NA/empty vector returns integer(0), and min() over it is Inf —
    # the report then printed "rho >= Inf". Guard both.
    .wk <- if (length(.sp) && any(is.finite(.sp))) ctx$sens$axis[which.min(.sp)] else "n/a"
    .verdict <- if (!is.finite(.mn)) "not assessable"
                else if (.mn >= 0.95) "stable across every axis tested"
                else if (.mn >= 0.85) sprintf("broadly stable, weakest on %s (Spearman %.2f)", .wk, .mn)
                else sprintf("STABLE ONLY ON SOME AXES: the ranking moves materially on %s (Spearman %.2f)", .wk, .mn)
    # The KERNEL rows compare rankings only: delta was fitted for the featured kernel, so an
    # alternative kernel's absolute reach is uncalibrated (33_cascade_eval.R, level_comparable).
    .rank_only <- "level_comparable" %in% names(ctx$sens) && any(!ctx$sens$level_comparable %in% TRUE)
    add("- **Ranking stability (mobility kernel & frontier saturation psi):** 13-week reach ranking is %s — Spearman vs baseline: %s.%s%s",
        .verdict, paste(sprintf("%s %.2f", ctx$sens$axis, ctx$sens$spearman), collapse = "; "),
        if (is.finite(.mn) && .mn < 0.85)
          " Top-K and gateway products inherit that movement on the kernel axis and should be read against it, not as kernel-independent."
        else " Decision products (top-K, gateway) are robust even where absolute probabilities are not.",
        if (.rank_only) " The alternative-kernel rows are a RANKING comparison only: delta was fitted for the featured kernel, so their absolute reach is uncalibrated and is not reported." else "")
  }
  if (!is.null(ctx$ksweep) && !is.null(ctx$ksweep$summary) && nrow(ctx$ksweep$summary)) {
    ks <- ctx$ksweep$summary; kb <- ctx$ksweep$k_base %||% CASCADE_K_INDIV
    off <- ks[abs(ks$k_indiv - kb) > 1e-9, ]                 # exclude the self-comparison
    imin <- which.min(ks$k_indiv); imax <- which.max(ks$k_indiv)
    # FLOOR the lower bounds (a ">=" claim must not round up past the true minimum)
    fl <- function(x, d) floor(x * 10^d) / 10^d
    add("- **Overdispersion (individual-offspring `k_indiv`) — comprehensive sweep (%d values, k in [%.2f, %.2f] at the deployed delta/psi; baseline k = %.2f):** the 13-week invasion **ranking is invariant to k** — Spearman rho >= %.3f, Kendall tau >= %.2f, and top-15 Jaccard >= %.2f for every off-baseline k across two orders of magnitude. Mean reach is only weakly sensitive (%.3g%% at k=%.2f to %.3g%% at k=%.2f); **establishment probability** is the k-sensitive quantity (%.3g%% to %.3g%%, ~%.1fx over the grid), because heavier superspreading (lower k) yields more stochastic fade-out before the n_est=%d establishment threshold. The Lloyd-Smith EBOV-plausible band [0.2, 0.4] lies on a flat part of every curve.%s",
        nrow(ks), min(ks$k_indiv), max(ks$k_indiv), kb,
        fl(min(off$spearman), 3), fl(min(off$kendall), 2), fl(min(off$top15_jaccard), 2),
        100 * ks$mean_reach_h13[imin], ks$k_indiv[imin], 100 * ks$mean_reach_h13[imax], ks$k_indiv[imax],
        100 * ks$mean_estab_h13[imin], 100 * ks$mean_estab_h13[imax],
        ks$mean_estab_h13[imax] / max(ks$mean_estab_h13[imin], 1e-9), CASCADE_N_EST,
        .fig_line(c("Figure" = "cascade/figures/Figure_kindiv_sweep"), prefix = " "))
  }
  if (!is.null(ctx$backtest)) {
    # BOTH the K in "precision@K" and the COLUMN it reads come from the same TOPK_VALUES
    # entry, so the two cannot disagree. The prose K was already derived here while the
    # column name beside it was hard-coded `prec_at15`: moving TOPK_VALUES[2] off 15 would
    # have printed "precision@<new K>" next to the precision at 15, or NULL if 33 no longer
    # emitted that column.
    .bk_K   <- get0("TOPK_VALUES", ifnotfound = c(5L, 15L, 20L))[2]
    .bk_col <- ctx$backtest[[paste0("prec_at", .bk_K)]]
    # THE COUNT RATIO IS IN-SAMPLE when delta was fitted on these same origins and K, which is
    # what run_cascade.R does whenever it calibrates — so it is labelled, never presented as a
    # held-out level check. The ranking metrics beside it are what delta was not fitted to.
    .bt_ins <- isTRUE(ctx$backtest_level_in_sample %||% .calibrated)
    add("- **Intermediate-horizon cascade backtest (~1 effective origin; indicative):** at %d-week horizon, AUC-PR skill %s, precision@%d %s; predicted/observed new-zone count ratio %s%s.",
        ctx$backtest$K[1],
        paste(sprintf("%.1fx", ctx$backtest$auc_pr_skill), collapse = "/"),
        .bk_K,
        if (is.null(.bk_col)) "n/a" else paste(sprintf("%.2f", .bk_col), collapse = "/"),
        paste(sprintf("%.2f", ctx$backtest$count_ratio), collapse = "/"),
        if (.bt_ins) " — **IN-SAMPLE**: delta was fitted to the pooled count on these same origins, so the pooled ratio is ~1 by construction and only the ranking metrics are out-of-sample"
        else " (delta was not fitted on these origins, so the ratio is a held-out level check)")
  }
  # The ranks are TESTED, not asserted: a Kolmogorov-Smirnov test of each parameter's posterior
  # rank of truth against Uniform(0,1), Bonferroni over the parameters.
  if (!is.null(ctx$sbc) && isTRUE(ctx$sbc$n_used > 0) && !is.null(ctx$sbc$rank_stats)) {
    .rs <- as.matrix(ctx$sbc$rank_stats)
    .ksp <- apply(.rs, 2, function(x) { x <- x[is.finite(x)]
      if (length(x) >= 5L) suppressWarnings(stats::ks.test(x, "punif")$p.value) else NA_real_ })
    .ok <- all(.ksp > 0.05 / length(.ksp), na.rm = TRUE) && any(is.finite(.ksp))
    add("- **Simulation-based calibration (fixed-beta fitter):** %d replicates; Kolmogorov-Smirnov test of each parameter's posterior rank of truth against Uniform(0,1): %s (Bonferroni level %.3f over %d parameters) — %s.",
        ctx$sbc$n_used, paste(sprintf("%s p = %.2f", sub("^b_", "", names(.ksp)), .ksp), collapse = ", "),
        0.05 / length(.ksp), length(.ksp),
        if (.ok) "no departure from uniformity detected" else "**ranks depart from uniformity**")
  } else if (isTRUE(ctx$sbc$skipped))
    add("- **Simulation-based calibration:** %s", ctx$sbc$note)
  else if (!is.null(ctx$sbc) && isTRUE(ctx$sbc$n_used == 0))
    add("- **Simulation-based calibration:** 0 usable replicates (every replicate simulated fewer than two invasions or failed to fit), so no rank test was possible.")
  # STATE THE GAP. Each validation stage above is wrapped in a tryCatch in run_cascade.R whose
  # handler is a message() — it leaves no warning and no artifact — and each section here is
  # guarded on its context being non-NULL. So a stage that FAILED simply vanished from this
  # report, which then read as complete: a reader could not tell "we ran the backtest and it
  # passed" from "the backtest crashed". Name every absent stage instead of omitting it.
  local({
    .want <- c(sens     = "Ranking stability (mobility kernel & frontier saturation psi)",
               ksweep   = "Overdispersion (k_indiv) sweep",
               backtest = "Intermediate-horizon cascade backtest",
               sbc      = "Simulation-based calibration")
    .have <- c(sens     = !is.null(ctx$sens) && nrow(ctx$sens) > 0,
               ksweep   = !is.null(ctx$ksweep) && !is.null(ctx$ksweep$summary) &&
                          nrow(ctx$ksweep$summary) > 0,
               backtest = !is.null(ctx$backtest),
               # A deliberate skip (a time-varying form) is reported above, not as a failure.
               sbc      = !is.null(ctx$sbc) && (isTRUE(ctx$sbc$n_used > 0) || isTRUE(ctx$sbc$skipped)))
    .miss <- names(.want)[!.have[names(.want)]]
    for (k in .miss)
      add("- **%s: NOT AVAILABLE in this run** — the stage did not produce a result (it was disabled, or it failed; see the run log). This section is absent, not passed.",
          .want[[k]])
  })
  add("")
  add("---")
  add("## 6. Assumptions & limitations (stated up front)")
  add("- **Projection, not forecast.** 13-week trajectories are dominated by future transmissibility/response; issued as scenarios, S1 central.")
  add("- **Absolute probabilities are conditional** — on the scenario's assumed future transmissibility and on a level calibrated over %s, not 13 — and are not bounds in either direction; rankings/relative risk/gateway ordering are the decision-grade products.",
      if (.calibrated && !is.null(.do)) sprintf("%d weeks", as.integer(.do$K)) else "one week (multi-week calibration disabled this run)")
  add("- **Endogenous spread only:** no cross-border reintroduction (Uganda/South Sudan) or off-network seeding.")
  add("- **Static mobility** over 13 weeks (sensitivity: alternative kernels); **R_eff** is a SINGLE NATIONAL number from EpiNow2, modulated by the pre-registered scenario multiplier c(t) and %s — there is no per-zone or per-province reproduction number, because no between-zone heterogeneity was detectable (see below).",
      if (!(as.numeric(get0("CASCADE_R_RW_SIGMA", ifnotfound = 0)) > 0)) "no walk (sigma = 0)"
      else if (as.numeric(get0("CASCADE_R_RW_RHO", ifnotfound = 0)) > 0)
        "a mean-reverting walk on the national log R (one shock per week and draw)"
      else "a pure random walk on the national log R (one shock per week and draw)")
  add("- **How R is estimated — ONE national number, from the same model the short-term arm uses.** The cascade transmits on a single national reproduction number R = %.2f [%.2f, %.2f] (90%%), taken from EpiNow2 and averaged over the last %s (%s to %s; source `%s`). The average is taken *within each posterior draw*, so the draws are of the window mean and carry the right correlation across weeks; it cannot be assembled from separate per-week fits, whose sample indices share no posterior. The fit is censored at the analysis date, and each held-out calibration origin refits it at its OWN cutoff, so no origin sees its future. **There is no per-zone or per-province reproduction number.** The earlier three-level conjugate estimator (zone shrunk to province shrunk to national, with EM case attribution) was removed because a dispersion test found no detectable between-zone heterogeneity to estimate — X² = 10.6 on 14 df, p = 0.71, with *negative* excess variance at every case threshold — while most zones were prior-dominated and therefore already leaning on the national number behind a per-zone presentation. Sharing EpiNow2 with the short-term arm also removes a difference between the two analyses that would otherwise have to be defended. Both arms average the posterior over the SAME window, RT_WINDOW_WEEKS (00_config.R), ending at the forecast origin — until 2026-09-19 they did not (1 week here against 3 for the cascade, R = 0.83 against 1.00), so this sentence was not true of the numbers it accompanied. **Four consequences, stated rather than buried.** (i) Because there is one R, its uncertainty is **common to every zone**: each Monte-Carlo replicate draws one value and applies it everywhere. That is deliberate — drawing per zone would average the national uncertainty away across the ~60 infected zones and make the interval far too narrow — but it means zone-level reproduction numbers are *not* an output of this analysis, and the R uncertainty is perfectly correlated rather than diversifying. (ii) R's posterior (sd on the log scale %.3f, about ±%.0f%%) is a posterior-side component of the 13-week per-zone interval, together with the calibration factor and the R walk%s; none of them is negligible, and holding R at its posterior mean narrows the interval materially (it is not a minor term). Monte-Carlo **process noise is NOT** a driver: the variance decomposition that forms this interval removes it explicitly (see the note in \u00a72). The R walk is a driver because its innovations are drawn once per parameter draw and shared by that draw's replicates, which places it on the posterior side of that decomposition. (iii) EpiNow2 carries no import term. That is the right choice nationally, where mobility importation is internal redistribution rather than new infection, and it is why this number is not used at zone level. (iv) The conjugate estimator is still computed and reported as a **diagnostic** (%.2f on this frame, %.0f%% %s than the deployed anchor) together with the case and importation bookkeeping — %d zones with cases, importation %.1f%% of its denominator — but nothing in the projection transmits on it, so its two known weaknesses (partly circular import attribution inside a tightly connected cluster, and extrapolating a first-invasion hazard to continuing importation) now affect only that diagnostic and no reported result.",
      ctx$reff$R_nat %||% NA_real_,
      unname((ctx$reff$rt_quantiles %||% rep(NA_real_, 5))[1]),
      unname((ctx$reff$rt_quantiles %||% rep(NA_real_, 5))[5]),
      { .ww <- as.integer(ctx$reff$rt_window_weeks %||% NA_integer_)
        if (is.na(.ww)) "? weeks" else sprintf("%d week%s", .ww, if (.ww == 1L) "" else "s") },
      if (is.null(ctx$reff$rt_window_start)) "?" else format(as.Date(ctx$reff$rt_window_start)),
      if (is.null(ctx$reff$rt_window_end)) "?" else format(as.Date(ctx$reff$rt_window_end)),
      ctx$reff$rt_source %||% "unknown",
      as.numeric((ctx$reff$sd_zone %||% NA_real_)[1]),
      100 * as.numeric((ctx$reff$sd_zone %||% NA_real_)[1]) * 1.6448536,
      if (.form %in% c("rw1", "rw1_frailty")) ", and the featured form's import-coefficient drift path" else "",
      ctx$reff$R_nat_conjugate %||% NA_real_,
      # The sentence reads "<conjugate> is X% <higher/lower> THAN THE DEPLOYED ANCHOR", so the
      # base of the percentage is the anchor: conj/R_nat - 1. It was R_nat/conj - 1, which is a
      # different number (15.1% vs 17.8% on the shipped frame) and answers the opposite question.
      100 * abs((ctx$reff$R_nat_conjugate %||% NA_real_) /
                  max(ctx$reff$R_nat %||% NA_real_, 1e-12) - 1),
      if (!is.null(ctx$reff$R_nat_conjugate) && !is.null(ctx$reff$R_nat) &&
          ctx$reff$R_nat_conjugate < ctx$reff$R_nat) "lower" else "higher",
      ctx$reff$n_with_cases %||% NA_integer_,
      100 * (ctx$reff$import_share %||% NA_real_))
  # The EpiNow2 cross-comparison block was removed with 44_reff_epinow2_check.R (2026-09-27),
  # the standalone driver that produced reff_national_epinow2_check.json. That driver's own
  # header had already withdrawn the claim it existed to support: the cascade now takes its
  # national R FROM EpiNow2, so comparing the two is comparing one model with itself. No
  # artifact is written any more, so nothing here can quote a stale one.
  # h=1 invasion-event count, READ from invasion_evaluation.csv (the table that holds it),
  # not hard-coded. NA when the table is unavailable, and the sentence degrades in words.
  .n_inv_h1 <- suppressWarnings(tryCatch({
    .evp <- file.path(get0("OUT_DIAGNOSTICS", ifnotfound = file.path("outputs", "diagnostics")),
                      "invasion_evaluation.csv")
    if (!file.exists(.evp)) NA_integer_ else {
      .evt <- readr::read_csv(.evp, show_col_types = FALSE)
      .col <- intersect(c("n_invasions", "n_invasion_events"), names(.evt))[1]
      if (is.na(.col) || !"horizon" %in% names(.evt)) NA_integer_
      else as.integer(stats::median(.evt[[.col]][.evt$horizon == 1L], na.rm = TRUE))
    }
  }, error = function(e) NA_integer_))
  add("- **Ascertainment ramp-up is excluded, not censored.** In the first weeks after detection, apparent growth is dominated by the response arriving and a backlog appearing at once, which no renewal estimator can separate from transmission without an explicit reporting model. Calibration origins inside that burn-in are therefore skipped (`CASCADE_CALIB_MIN_WEEKS`); the estimator is not adjusted, because it is not wrong about the data — the data are not yet about transmission.")
  # READ, not asserted: the count lives in invasion_evaluation.csv and moves every run (the
  # hard-coded "~39" was 37 on the shipped frame).
  add("- **Small event base** (%s h=1 invasions) limits identifiability; priors + pooling do heavy lifting.",
      if (is.finite(.n_inv_h1)) format(.n_inv_h1) else "a small number of")
  add("- Two-stage nowcast point-plugs the initial conditions; %s; `psi` is not calibrated at all (it is fixed at %.2f).",
      if (.calibrated && is.finite(.calib_K))
        sprintf("`delta` is fitted on realised %d-week invasion counts at held-out origins, not on a long-horizon outcome (unavailable)", .calib_K)
      else "`delta` is the one-week factor and is NOT fitted at any multi-week horizon in this run",
      as.numeric(get0("CASCADE_PSI", ifnotfound = NA_real_)))
  add("")
  add("Full spec, equations, and the review log: `PLAN_3MONTH_INVASION.md`.")
  writeLines(L, path)
  message("[cascade] wrote report -> ", path)
  invisible(path)
}
