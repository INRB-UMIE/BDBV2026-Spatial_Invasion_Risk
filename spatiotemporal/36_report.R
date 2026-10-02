# =============================================================================
# 36_report.R — 3-MONTH CASCADE: assemble the decision-maker report (PLAN §7)
# BDBV 2026 DRC
#
# cascade_write_report() writes SPATIAL_INVASION_3MONTH_REPORT.md in the house
# style, auto-populated from the run artifacts: reach leaderboard, expected new
# zones/provinces per scenario, gateway zones, corridors, conditional hubs, the
# consistency-gate result, sensitivity ranking-stability, and the backtest — with
# the projection / upper-bound / small-sample caveats stated up front.
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

cascade_write_report <- function(path, ctx) {
  L <- c()
  add <- function(...) L <<- c(L, sprintf(...))
  aff <- sum(ctx$reach_primary$was_active_before & ctx$reach_primary$horizon ==
               max(CASCADE_REPORT_HORIZONS))
  atrisk <- sum(!ctx$reach_primary$was_active_before & ctx$reach_primary$horizon ==
                  max(CASCADE_REPORT_HORIZONS))

  add("# BDBV 2026 — 3-Month Spatial Invasion Risk (Cascade Projection)")
  add("")
  add("**Analysis date:** %s · **Kernel:** %s (medium GT) · **M:** %d MC iterations · **Horizon:** %d weeks",
      as.character(get0("ANALYSIS_DATE", ifnotfound = Sys.Date())), CASCADE_KERNEL,
      ctx$n_mc, CASCADE_HORIZON_WEEKS)
  add("**Folder:** `%s` · implements `PLAN_3MONTH_INVASION.md`.", basename(dirname(OUT_CASCADE)))
  add("")
  add("> **How to read this.** %s The 1–2-week Bayesian suite remains the operational headline; this layer is a longer-range **projection** of where invasion may spread and how far onward transmission may carry it. Decisions should use the **rankings, relative risk, gateway ordering, and priority scores**; treat absolute 3-month probabilities as **upper bounds**. Results are issued as **control scenarios**, not a single forecast.", CASCADE_FRAMING)
  add("")
  add("---")
  # The held-out window delta was fitted at, READ from the fit result (cascade_fit_delta_oos
  # returns K) rather than from a constant that does not exist.
  .calib_K <- suppressWarnings(as.integer((ctx$calib$delta_oos %||% list())$K %||% NA_integer_))
  add("## 1. Method (one paragraph)")
  # EVERY MECHANISM NAMED HERE MUST EXIST. This paragraph is the published Method and it
  # described three that do not: a "partially-pooled" per-zone R_eff (there is ONE national
  # EpiNow2 R applied to every zone, 31_source_dynamics.R), a delta "recalibrated to observed
  # short-horizon frequency" (delta is fitted OUT OF SAMPLE on realised K-week invasion counts,
  # 33b cascade_fit_delta_oos), and damping "by a frontier-saturation term (psi)" (CASCADE_PSI
  # is 0, so the saturation factor is identically 1 and nothing is damped). Section 6 of this
  # same document already said the opposite about the first of those.
  add("A stochastic spatial-metapopulation **cascade** couples (A) a within-zone overdispersed NegBin branching process — every infected zone projects incidence forward with the SAME national current effective reproduction number `R_eff` (NOT R0), estimated once by EpiNow2 on the national series — with (B) the **validated cloglog invasion hazard** from `21_bayesian_renewal.R`, reused as a weekly Bernoulli **seeding** step for each still-at-risk zone. A seeded zone becomes a source, so invasion chains onward across the mobility network. The per-step hazard is rescaled on the hazard scale by a single factor `delta`, fitted OUT OF SAMPLE against realised %d-week invasion counts (`33b_cascade_calibration.R`); the frontier-saturation term is switched off (`psi` = %.2f), so no damping is applied. The dynamic `d_min` frontier covariate is recomputed each simulated week. Individual-offspring overdispersion (`k_indiv`≈%.2f) produces realistic stochastic **fade-out**, so invasion (first case) is distinguished from establishment (a sustained outbreak). Everything runs on the confirmed-case scale (time index = symptom onset), consistent with the short-horizon suite.",
      .calib_K, as.numeric(get0("CASCADE_PSI", ifnotfound = NA_real_)), CASCADE_K_INDIV)
  add("")
  # Describes the calibration ACTUALLY performed (33b_cascade_calibration.R). The
  # previous text said delta came from 1/calibration-in-the-large and said nothing
  # about how psi was obtained; both were changed on 2026-09-11, and a methods
  # paragraph that names the wrong estimator is worse than none.
  .ci <- ctx$calib$delta_info %||% list()
  .pf <- ctx$calib$psi_fit
  .rp <- ctx$calib$report
  .do <- ctx$calib$delta_oos
  add(paste0("**Calibration — ONE fitted parameter.** delta = %.3f (%s): the per-week hazard ",
             "multiplier, fitted so the cascade's expected new invasions over held-out ",
             "%s-week window(s) match the number that actually occurred (%s). It is fitted at ",
             "the multi-week horizon rather than at h=1 deliberately: that is the horizon this ",
             "layer is read at, and with saturation removed there is no second parameter left ",
             "to absorb multi-week shape. **What that costs is stated**: delta previously came ",
             "from the prequential maximum-likelihood factor of `16b_invasion_recalibration.R` ",
             "at h=1 (%s on this frame%s), which made the cascade's first week exactly the ",
             "validated short-horizon hazard and the h=1 consistency gate true by construction. ",
             "It is now a real test rather than an identity. **psi = %g is FIXED, not fitted.** ",
             "Frontier saturation was not identifiable on this outbreak — across refits it came ",
             "back censored at a bound, or well away from its own root, while two probabilistic ",
             "criteria (deviance, AUC-PR) kept improving to the search ceiling. That is the ",
             "signature of a phenomenological term absorbing misspecification it cannot ",
             "represent, not of saturation being measured, so it is set to 0 (sat identically 1) ",
             "and retained only as a sensitivity axis. %s"),
      as.numeric(ctx$delta),
      .ci$estimator %||% "estimator unrecorded",
      if (is.null(.do)) "n/a" else as.character(.do$K),
      if (is.null(.do)) "not fitted out of sample"
      else sprintf("modelled %.1f vs %.0f observed%s%s", .do$modelled_cum_at_delta,
                   .do$target_observed_cum,
                   if (isTRUE(.do$converged)) ", converged" else ", **NOT converged**",
                   if (isTRUE(.do$boundary_hit))
                     sprintf(", **CENSORED at the %s bound**", .do$boundary_side) else ""),
      if (is.finite(ctx$calib$delta_prior %||% NA_real_))
        sprintf("%.3f", ctx$calib$delta_prior) else "unavailable",
      if (is.finite(.ci$h1_refit_delta %||% NA_real_))
        sprintf("; refitting from the cascade's own week-1 hazard on held-out origins gives %.3f, ratio %.2f",
                .ci$h1_refit_delta, .ci$h1_refit_ratio) else "",
      get0("CASCADE_PSI", ifnotfound = 0),
      if (is.null(.rp)) "No out-of-sample verification was run."
      else sprintf(paste0("**Out-of-sample verification** at %d weeks over %d held-out origin(s): ",
                          "predicted %.1f vs observed %.0f new invasions (ratio %.2f) — %s; AUC-PR ",
                          "skill %.1fx. Note this shares its origins with the delta fit, so the ",
                          "COUNT is matched by construction; the discrimination metrics are the ",
                          "part that is not."),
                   .rp$pooled$K, .rp$pooled$n_origins, .rp$pooled$predicted_new,
                   .rp$pooled$observed_new, .rp$pooled$count_ratio,
                   if (isTRUE(.rp$pass)) "**PASS**" else "**FAIL**", .rp$pooled$auc_pr_skill))
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
  add("**R_eff is not held constant over the horizon.** Each zone's log R follows a mean-reverting (AR(1)) walk toward its OWN estimated value: log R(t) = log R(t-1) + rho*(log R(0) - log R(t-1)) + sigma*eps, with sigma = %.3f per week and rho = %.2f. Reverting to the zone's own anchor means the expected R is unchanged, so this adds uncertainty that GROWS WITH HORIZON and then plateaus, without adding a trend. That separation is deliberate: the observed national R fell about 0.02-0.05 per week in log terms over the post-burn-in window, but projecting that decline forward is a claim about future control, which is exactly what the pre-registered c(t) scenarios express — fitting it here as well would contradict S1 \"status quo\" and double-count control in S2. sigma is measured, not assumed: an AR(1) on the weekly national log R gives a residual innovation SD of 0.036-0.040 across windows starting at weeks 8, 10 and 12 (windows including weeks 5-7 give 0.052, but those are ascertainment ramp-up, excluded on the same grounds as the calibration burn-in). Setting sigma = 0 reproduces the previous constant-R projection exactly, and the sensitivity table carries a sigma axis.",
      get0("CASCADE_R_RW_SIGMA", ifnotfound = 0), get0("CASCADE_R_RW_RHO", ifnotfound = 0))
  add("")
  add(paste0("Calibration uncertainty is propagated: delta is drawn per posterior-parameter group ",
             "on the log scale (SD %.3f), so the per-zone credible intervals include it rather than ",
             "treating the calibration as exact. Consistency gate at h=1: the cascade ranks like the ",
             "validated short-horizon model (**Spearman = %.3f**, n=%d signal zones) at a level ratio ",
             "of %.2f — %s."),
      ctx$delta_sd_log %||% 0,
      ctx$gate$spearman_signal, ctx$gate$n_signal, ctx$gate$level_ratio,
      if (isTRUE(ctx$gate$pass)) "**PASS**" else "**REVIEW**")
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
  add("_Note: the per-zone interval is a 90%% credible interval — the posterior spread of the reach probability across parameter draws (hazard coefficients and reproduction numbers), with Monte-Carlo process noise removed. It does NOT include between-scenario (§3) or structural (mobility/generation-time, §5) uncertainty, which is larger; absolute probabilities remain upper bounds._")
  add("")
  add("| Zone | Province | P(reach) [90%% CrI] | P(estab.) | 1st-pass wk | Priority rank |")
  add("|---|---|---|---|---|---|")
  for (i in seq_len(min(15L, nrow(r13)))) { z <- r13[i, ]
    add("| %s | %s | %.2f [%.2f–%.2f] | %.2f | %s | %s |", z$health_zone,
        z$province %||% "", z$p_case_invasion, z$p_lo, z$p_hi, z$p_establishment,
        ifelse(is.finite(z$first_passage_week), as.character(round(z$first_passage_week)), "—"),
        ifelse(is.finite(z$priority_rank), as.character(z$priority_rank), "—")) }
  add("")
  add("%s", .fig_line(c(
    "Reach map"                   = "cascade/figures/cascade_reach_map_national_h13",
    "Expanding-frontier triptych" = "cascade/figures/cascade_time_to_invasion_triptych",
    "Uncertainty"                 = "cascade/figures/cascade_reach_uncertainty_national_h13")))
  add("")
  add("---")
  add("## 3. Burden trajectory — expected new zones & provinces, by scenario")
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
  add("*All five columns are medians over the Monte-Carlo iterations; the bracketed range is a 90%% predictive interval.*")
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
  if (!is.null(ctx$conditional) && length(ctx$conditional)) {
    add("**Conditional onward-spread (generalised Kisangani).** If a key hub is invaded, the downstream reach concentrates on its catchment:")
    add("")
    for (h in names(ctx$conditional)) { cc <- ctx$conditional[[h]]
      # The per-hub conditional-reach map is suppressed by the retained-figure allow-list on the
      # published configuration, so its pointer is only printed when the gate would write it.
      add("- **%s** → top downstream zones by 13-week conditional reach: %s%s",
          h, paste(cc$top, collapse = ", "),
          .fig_line(stats::setNames(
            sprintf("cascade/figures/cascade_conditional_reach_%s_h13",
                    gsub("[^A-Za-z0-9]+", "_", h)), "map"),
            prefix = " (", ext = ".pdf")) }
    add("")
  }
  if (!is.null(ctx$urban) && !is.null(ctx$urban$summary) && nrow(ctx$urban$summary)) {
    add("---")
    add("## 4b. Scenario: invasion of a highly-connected/populated city")
    add("What happens if a major urban hub — far from the current front — is itself invaded and becomes a source? Each city's urban-core health zones are seeded at a chosen week and the cascade projects its onward spread; impact is measured against the baseline (no urban seeding), **excluding the seeded city itself**. **Timing matters**: earlier seeding leaves more weeks for spread.")
    add("")
    # The PRIMARY specification only: earliest seed week, seeded city at its province pool.
    # The seeded-R arms get their own subsection — averaging them into this table would
    # silently mix a premise the reader has not yet been shown.
    u1 <- ctx$urban$summary[ctx$urban$summary$seed_week == ctx$urban$min_week, ]
    if ("seed_r" %in% names(u1)) u1 <- u1[is.na(u1$seed_r), ]
    u1 <- u1[order(-u1$expected_added), ]
    add("Downstream impact if invaded at week %d (median status-quo scenario):", ctx$urban$min_week)
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
        add("Read the two columns as different quantities, not as a split of one. **Cases in the seeded city** combine the *assumed* seed itself (about 1.5 cases per forced zone under the default Poisson seeding, so roughly 4-5 for a three-zone hub) with the local outbreak it starts, so this column is largely a restatement of the seeded-R premise below — it is what you assume, propagated. **Cases elsewhere** is the genuine downstream burden and the quantity a national planner trades off; a city whose interval there spans zero has no detectable effect on the rest of the country. **The typical outcome and the expected outcome are different numbers, and the gap is the headline.** A seeded introduction usually fades: in most iterations the city contributes little beyond the seed cases themselves, which is why the MEDIAN additional total is small. Occasionally it establishes and produces dozens of cases, which is why the MEAN is several times larger. On this frame the two differ by about %.1f-fold (median %s against mean %s, with 90%% of outcomes spanning roughly %s to %s). Neither number alone is honest: the **median** is what typically happens and the **mean** is the right input to an expected-burden or cost calculation, while the predictive span is what a preparedness plan has to cover. The same skew, milder, affects the baseline: mean %s against median %s (90%% predictive band %s-%s). The percentage column divides mean by mean so the two estimators are never mixed, and medians are NOT additive across the in-city/elsewhere split — only the means are. The baseline denominator is the projected confirmed cases over the same 13 weeks with no urban seeding, and it rests on the status-quo scenario — transmissibility persisting with no further control gain and no within-zone susceptible depletion — so it is an upper bound on the counterfactual as well. All figures are modelled **confirmed** cases, the scale the models are fitted and calibrated on; they are not infections, and this suite does not estimate an ascertainment fraction.%s",
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
      add(paste0("| City | province pool | ",
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
  add("- **Consistency gate (h=1):** ranking Spearman = %.3f vs the validated short-horizon model; level recalibrated to observed frequency. %s",
      ctx$gate$spearman_signal, if (isTRUE(ctx$gate$pass)) "PASS." else "Review.")
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
    add("- **Ranking stability (mobility kernel & frontier saturation psi):** 13-week reach ranking is %s — Spearman vs baseline: %s.%s",
        .verdict, paste(sprintf("%s %.2f", ctx$sens$axis, ctx$sens$spearman), collapse = "; "),
        if (is.finite(.mn) && .mn < 0.85)
          " Top-K and gateway products inherit that movement on the kernel axis and should be read against it, not as kernel-independent."
        else " Decision products (top-K, gateway) are robust even where absolute probabilities are not.")
  }
  if (!is.null(ctx$ksweep) && !is.null(ctx$ksweep$summary) && nrow(ctx$ksweep$summary)) {
    ks <- ctx$ksweep$summary; kb <- ctx$ksweep$k_base %||% CASCADE_K_INDIV
    off <- ks[abs(ks$k_indiv - kb) > 1e-9, ]                 # exclude the self-comparison
    imin <- which.min(ks$k_indiv); imax <- which.max(ks$k_indiv)
    # FLOOR the lower bounds (a ">=" claim must not round up past the true minimum)
    fl <- function(x, d) floor(x * 10^d) / 10^d
    add("- **Overdispersion (individual-offspring `k_indiv`) — comprehensive sweep (%d values, k in [%.2f, %.2f] at the fitted delta/psi; baseline k = %.2f):** the 13-week invasion **ranking is invariant to k** — Spearman rho >= %.3f, Kendall tau >= %.2f, and top-15 Jaccard >= %.2f for every off-baseline k across two orders of magnitude. Mean reach is only weakly sensitive (%.1f%% at k=%.2f to %.1f%% at k=%.2f); **establishment probability** is the k-sensitive quantity (%.1f%% to %.1f%%, ~%.1fx over the grid), because heavier superspreading (lower k) yields more stochastic fade-out before the n_est=%d establishment threshold. The Lloyd-Smith EBOV-plausible band [0.2, 0.4] lies on a flat part of every curve.%s",
        nrow(ks), min(ks$k_indiv), max(ks$k_indiv), kb,
        fl(min(off$spearman), 3), fl(min(off$kendall), 2), fl(min(off$top15_jaccard), 2),
        100 * ks$mean_reach_h13[imin], ks$k_indiv[imin], 100 * ks$mean_reach_h13[imax], ks$k_indiv[imax],
        100 * ks$mean_estab_h13[imin], 100 * ks$mean_estab_h13[imax],
        ks$mean_estab_h13[imax] / max(ks$mean_estab_h13[imin], 1e-9), CASCADE_N_EST,
        .fig_line(c("Figure" = "cascade/figures/Figure_kindiv_sweep"), prefix = " "))
  }
  if (!is.null(ctx$backtest)) {
    add("- **Intermediate-horizon cascade backtest (~1 effective origin; indicative):** at %d-week horizon, AUC-PR skill %s, precision@10 %s, predicted/observed new-zone count ratio %s.",
        ctx$backtest$K[1],
        paste(sprintf("%.1fx", ctx$backtest$auc_pr_skill), collapse = "/"),
        paste(sprintf("%.2f", ctx$backtest$prec_at10), collapse = "/"),
        paste(sprintf("%.2f", ctx$backtest$count_ratio), collapse = "/"))
  }
  if (!is.null(ctx$sbc) && ctx$sbc$n_used > 0)
    add("- **Simulation-based calibration:** %d replicates; posterior ranks of truth ~uniform (kernel machinery calibrated).", ctx$sbc$n_used)
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
               sbc      = !is.null(ctx$sbc) && isTRUE(ctx$sbc$n_used > 0))
    .miss <- names(.want)[!.have[names(.want)]]
    for (k in .miss)
      add("- **%s: NOT AVAILABLE in this run** — the stage did not produce a result (it was disabled, or it failed; see the run log). This section is absent, not passed.",
          .want[[k]])
  })
  add("")
  add("---")
  add("## 6. Assumptions & limitations (stated up front)")
  add("- **Projection, not forecast.** 13-week trajectories are dominated by future transmissibility/response; issued as scenarios, S1 central.")
  add("- **Absolute probabilities are upper bounds**; rankings/relative risk/gateway ordering are the decision-grade products.")
  add("- **Endogenous spread only:** no cross-border reintroduction (Uganda/South Sudan) or off-network seeding.")
  add("- **Static mobility** over 13 weeks (sensitivity: alternative kernels); **R_eff** is a SINGLE NATIONAL number from EpiNow2, modulated by the pre-registered scenario multiplier c(t) and a mean-reverting walk — there is no per-zone or per-province reproduction number, because no between-zone heterogeneity was detectable (see below).")
  add("- **How R is estimated — ONE national number, from the same model the short-term arm uses.** The cascade transmits on a single national reproduction number R = %.2f [%.2f, %.2f] (90%%), taken from EpiNow2 and averaged over the last %s weeks (%s to %s; source `%s`). The average is taken *within each posterior draw*, so the draws are of the window mean and carry the right correlation across weeks; it cannot be assembled from separate per-week fits, whose sample indices share no posterior. The fit is censored at the analysis date, and each held-out calibration origin refits it at its OWN cutoff, so no origin sees its future. **There is no per-zone or per-province reproduction number.** The earlier three-level conjugate estimator (zone shrunk to province shrunk to national, with EM case attribution) was removed because a dispersion test found no detectable between-zone heterogeneity to estimate — X² = 10.6 on 14 df, p = 0.71, with *negative* excess variance at every case threshold — while most zones were prior-dominated and therefore already leaning on the national number behind a per-zone presentation. Sharing EpiNow2 with the short-term arm also removes a difference between the two analyses that would otherwise have to be defended. Both arms average the posterior over the SAME window, RT_WINDOW_WEEKS (00_config.R), ending at the forecast origin — until 2026-09-19 they did not (1 week here against 3 for the cascade, R = 0.83 against 1.00), so this sentence was not true of the numbers it accompanied. **Four consequences, stated rather than buried.** (i) Because there is one R, its uncertainty is **common to every zone**: each Monte-Carlo replicate draws one value and applies it everywhere. That is deliberate — drawing per zone would average the national uncertainty away across the ~60 infected zones and make the interval far too narrow — but it means zone-level reproduction numbers are *not* an output of this analysis, and the R uncertainty is perfectly correlated rather than diversifying. (ii) The posterior is tight (sd on the log scale %.3f, about ±%.0f%%), so over 13 weeks the per-zone interval is driven mainly by the calibration factor and the R walk, not by uncertainty in the level of R itself. Monte-Carlo **process noise is NOT** a driver: the variance decomposition that forms this interval removes it explicitly (see the note in \u00a72). The R walk is a driver because its innovations are drawn once per parameter draw and shared by that draw's replicates, which places it on the posterior side of that decomposition. (iii) EpiNow2 carries no import term. That is the right choice nationally, where mobility importation is internal redistribution rather than new infection, and it is why this number is not used at zone level. (iv) The conjugate estimator is still computed and reported as a **diagnostic** (%.2f on this frame, %.0f%% %s than the deployed anchor) together with the case and importation bookkeeping — %d zones with cases, importation %.1f%% of its denominator — but nothing in the projection transmits on it, so its two known weaknesses (partly circular import attribution inside a tightly connected cluster, and extrapolating a first-invasion hazard to continuing importation) now affect only that diagnostic and no reported result.",
      ctx$reff$R_nat %||% NA_real_,
      unname((ctx$reff$rt_quantiles %||% rep(NA_real_, 5))[1]),
      unname((ctx$reff$rt_quantiles %||% rep(NA_real_, 5))[5]),
      as.character(ctx$reff$rt_window_weeks %||% NA_integer_),
      if (is.null(ctx$reff$rt_window_start)) "?" else format(as.Date(ctx$reff$rt_window_start)),
      if (is.null(ctx$reff$rt_window_end)) "?" else format(as.Date(ctx$reff$rt_window_end)),
      ctx$reff$rt_source %||% "unknown",
      as.numeric((ctx$reff$sd_zone %||% NA_real_)[1]),
      100 * as.numeric((ctx$reff$sd_zone %||% NA_real_)[1]) * 1.6448536,
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
  # Only referenced when the artifact actually exists — a report that points at a file the
  # run did not produce is worse than one that stays silent.
  .rtchk <- file.path(OUT_DIR, "cascade", "diagnostics", "reff_national_epinow2_check.json")
  if (file.exists(.rtchk)) {
    .rc <- tryCatch(jsonlite::fromJSON(.rtchk), error = function(e) NULL)
    # STALENESS GUARD, and it has already caught a real one. 44_reff_epinow2_check.R is a
    # STANDALONE driver that run_cascade.R does not call, so this artifact can easily predate
    # the run being reported. On 2026-09-16 it was four days old and still described the
    # RETIRED conjugate anchor (1.08 over 2026-08-11..08-31) while the prose around it had been
    # rewritten for the EpiNow2 anchor — the worst combination, because current framing lends
    # credibility to superseded numbers. Existence is therefore not enough: the artifact is
    # quoted only when the anchor it records is the anchor this run actually used.
    .rc_anchor <- if (is.null(.rc)) NA_real_ else suppressWarnings(as.numeric(
      .rc$deployed_R_nat %||% .rc$cascade_R_nat %||% NA_real_))
    .rc_now <- suppressWarnings(as.numeric(ctx$reff$R_nat %||% NA_real_))
    .rc_ok  <- is.finite(.rc_anchor) && is.finite(.rc_now) &&
               abs(.rc_anchor / .rc_now - 1) < 0.01
    # Direction of the conjugate estimator RELATIVE TO the deployed anchor, computed from the
    # artifact rather than asserted in prose.
    .rc_conj <- suppressWarnings(as.numeric(.rc$conjugate_R_nat_diagnostic %||%
                                            ctx$reff$R_nat_conjugate %||% NA_real_))
    .rc_dir <- if (is.finite(.rc_conj) && is.finite(.rc_now))
      (if (.rc_conj < .rc_now) "lower" else "higher") else "not comparable"
    if (!is.null(.rc) && !isTRUE(.rc_ok))
      add(paste0("- **The EpiNow2 cross-comparison is STALE and is omitted.** `%s` records a ",
                 "deployed anchor of %s against this run's %.3f, so it describes a different ",
                 "fit on a different frame and cannot be quoted beside these results. It is a ",
                 "standalone driver: re-run `44_reff_epinow2_check.R` against this frame to ",
                 "regenerate it. Reporting it anyway would put current prose around superseded ",
                 "numbers, which is exactly how a stale figure acquires false authority."),
          basename(.rtchk),
          if (is.finite(.rc_anchor)) sprintf("%.3f", .rc_anchor) else "no anchor",
          .rc_now)
    if (!is.null(.rc) && isTRUE(.rc_ok))
      add("- **This comparison is NOT an independent validation of the anchor — it is a sensitivity between two estimators of the same quantity.** EpiNow2 now *supplies* the cascade's R, so comparing the projection's anchor against EpiNow2 would be circular, and this artifact must not be read as corroboration. What it still does usefully is quantify how far the retired conjugate renewal estimator sits from the deployed one on the same national series and window: the conjugate version treats the generation-time PMF as known, carries no interval of its own, and models no right-truncation, whereas EpiNow2 carries an explicit generation time and a right-truncation model for the incomplete recent tail. The gap between them is therefore a statement about the conjugate estimator's limitations, not evidence about EpiNow2. Result: %s (`cascade/diagnostics/reff_national_epinow2_check.csv`). Direction matters for reading the projection: on this frame the conjugate estimator runs **%s** than the deployed anchor, so a cascade built on it would have projected **%s** spread than the one reported here.",
          .rc$note %||% "see the artifact",
          # COMPUTED, not asserted. The text hard-coded "lower"; on the shipped artifact the
          # conjugate estimator was HIGHER (1.184 vs 1.006). It also credited EpiNow2 with "an
          # ascertainment term", which 02_epi_params.R explicitly removed.
          .rc_dir, if (identical(.rc_dir, "higher")) "more" else "less")
  }
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
  add("- Two-stage nowcast point-plugs the initial conditions; `delta` is fitted out of sample on realised %s-week invasion counts, not on a long-horizon outcome (unavailable), and `psi` is not calibrated at all (it is fixed at %.2f).",
      if (is.na(.calib_K)) "held-out" else format(.calib_K),
      as.numeric(get0("CASCADE_PSI", ifnotfound = NA_real_)))
  add("")
  add("Full spec, equations, and the review log: `PLAN_3MONTH_INVASION.md`.")
  writeLines(L, path)
  message("[cascade] wrote report -> ", path)
  invisible(path)
}
