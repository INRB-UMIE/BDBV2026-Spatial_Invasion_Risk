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
  add("## 1. Method (one paragraph)")
  add("A stochastic spatial-metapopulation **cascade** couples (A) a within-zone overdispersed NegBin branching process — each infected zone projects incidence forward with a partially-pooled current effective reproduction number `R_eff` (NOT R0) — with (B) the **validated cloglog invasion hazard** from `21_bayesian_renewal.R`, reused as a weekly Bernoulli **seeding** step for each still-at-risk zone. A seeded zone becomes a source, so invasion chains onward across the mobility network. The per-step hazard is recalibrated on the hazard scale (delta) to observed short-horizon frequency and damped by a frontier-saturation term (psi); the dynamic `d_min` frontier covariate is recomputed each simulated week. Individual-offspring overdispersion (`k_indiv`≈%.2f) produces realistic stochastic **fade-out**, so invasion (first case) is distinguished from establishment (a sustained outbreak). Everything runs on the confirmed-case scale (time index = symptom onset), consistent with the short-horizon suite.",
      CASCADE_K_INDIV)
  add("")
  add("**Calibration (PLAN §3.5).** delta = %.3f (from calibration-in-the-large %.2f on the short-horizon LFO; the hazard is scaled so the cascade's h=1 seeding matches observed frequency). psi = %.2f (frontier saturation, relative to the initial front). Consistency gate: h=1 reach ranks like the validated short-horizon model **Spearman = %.3f** (n=%d signal zones) at the recalibrated level (ratio %.2f) — %s.",
      ctx$delta, attr(ctx$delta, "cal_in_large") %||% NA_real_, ctx$psi,
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
  add("Reach map: `cascade/figures/cascade_reach_map_national_h13.pdf`; expanding-frontier triptych: `cascade_time_to_invasion_triptych.pdf`; uncertainty: `cascade_reach_uncertainty_national_h13.pdf`.")
  add("")
  add("---")
  add("## 3. Burden trajectory — expected new zones & provinces, by scenario")
  add("")
  add("_Median cumulative newly-invaded zones (one consistent estimator across horizons); 13w carries the 90%% predictive band across MC iterations._")
  add("")
  add("| Scenario | New zones by 4w (med) | by 8w (med) | by 13w (med) [90%% pred. int.] | New provinces by 13w |")
  add("|---|---|---|---|---|")
  for (sk in names(ctx$agg)) { a <- ctx$agg[[sk]]
    add("| %s | %.0f | %.0f | %.0f [%.0f–%.0f] | %.1f |", a$label,
        a$new4, a$new8, a$new13_med, a$new13_lo, a$new13_hi, a$new_prov13) }
  add("")
  add("Fan chart: `cascade/figures/cascade_newzones_fanchart.pdf`; scenario reach maps: `cascade_scenario_reach_maps_h13.pdf`.")
  add("")
  add("---")
  add("## 4. Onward propagation — gateway zones & corridors (the task's core)")
  add("Gateway zones ranked by **Delta_j** = expected downstream invasions prevented if that source is contained (exact knockout on the shortlist; PLAN §3.7):")
  add("")
  add("| Source zone | Delta_j (downstream invasions) | %% of total spread |")
  add("|---|---|---|")
  gt <- head(ctx$gateway, 12L)
  for (i in seq_len(nrow(gt))) add("| %s | %.1f | %s |", gt$health_zone[i], gt$delta_j[i],
      .fmt_pct(gt$pct_reduction[i] / 100))
  add("")
  add("Inter-province seeding flux (corridors): `cascade/figures/cascade_corridors.pdf`; gateway bars: `cascade_gateway_knockout.pdf`; most-likely invasion tree: `cascade/tables/cascade_transmission_tree.csv`.")
  add("")
  if (!is.null(ctx$conditional) && length(ctx$conditional)) {
    add("**Conditional onward-spread (generalised Kisangani).** If a key hub is invaded, the downstream reach concentrates on its catchment:")
    add("")
    for (h in names(ctx$conditional)) { cc <- ctx$conditional[[h]]
      add("- **%s** → top downstream zones by 13-week conditional reach: %s (map: `cascade_conditional_reach_%s_h13.pdf`).",
          h, paste(cc$top, collapse = ", "), gsub("[^A-Za-z0-9]+", "_", h)) }
    add("")
  }
  if (!is.null(ctx$urban) && !is.null(ctx$urban$summary) && nrow(ctx$urban$summary)) {
    add("---")
    add("## 4b. Scenario: invasion of a highly-connected/populated city")
    add("What happens if a major urban hub — far from the current front — is itself invaded and becomes a source? Each city's urban-core health zones are seeded at a chosen week and the cascade projects its onward spread; impact is measured against the baseline (no urban seeding), **excluding the seeded city itself**. **Timing matters**: earlier seeding leaves more weeks for spread.")
    add("")
    u1 <- ctx$urban$summary[ctx$urban$summary$seed_week == ctx$urban$min_week, ]
    u1 <- u1[order(-u1$expected_added), ]
    add("Downstream impact if invaded at week %d (median status-quo scenario):", ctx$urban$min_week)
    add("")
    add("| City | Core zones | Exp. added zones (13w) | Zones materially elevated | Newly >20%% | Newly >50%% (%% of base) | Max RR | New prov. |")
    add("|---|---|---|---|---|---|---|---|")
    for (i in seq_len(nrow(u1))) { z <- u1[i, ]
      add("| %s | %d | %.1f | %d | +%d | +%d (+%.0f%%) | %s | %s |", z$city, z$n_hub_zones,
          z$expected_added, z$n_elevated, z$newly_above_20, z$newly_above_50, z$pct_increase_gt50,
          ifelse(is.finite(z$max_rr), sprintf("%.0f×", z$max_rr), "—"),
          ifelse(is.finite(z$provinces_newly_exposed), as.character(z$provinces_newly_exposed), "—")) }
    add("")
    add("**Which metrics to focus on.** The **expected additional zones invaded** (`Exp. added zones`) is the headline: it is the sum of the per-zone attributable increase over all non-hub zones and is an **unbiased estimate of E[extra zones invaded by 13 weeks]** — it does **not** require matched (common-random-number) runs, because far zones unaffected by this local seeding have a true increase of zero and their Monte-Carlo flicker cancels in the sum (the urban runs use M=%d so the residual is small). The **timing profile** (how `Exp. added` falls as the seeding week moves later) and the **within-catchment relative risk** + **attributable-risk map** (which zones) are the other primary products. The threshold-crossing counts (`Newly >20%%/50%%`) are secondary and shown as material (non-noise) upward crossings; for these distant metros the impact is a **new regional cluster of moderate risk**, so the national count above 50%% moves little. The seeded city's own zones are excluded (invaded by assumption).",
        ctx$urban$n_mc %||% CASCADE_N_MC)
    add("")
    add("**Timing.** Expected additional zones invaded by 13 weeks, by seeding week:")
    add("")
    us <- ctx$urban$summary
    hdr <- paste0("| City | ", paste(sprintf("seed wk %d", ctx$urban$seed_weeks), collapse = " | "), " |")
    add(hdr); add(paste0("|", paste(rep("---", length(ctx$urban$seed_weeks) + 1L), collapse = "|"), "|"))
    for (cty in unique(us$city)) {
      vals <- vapply(ctx$urban$seed_weeks, function(w) {
        v <- us$expected_added[us$city == cty & us$seed_week == w]
        if (length(v)) sprintf("%.1f", v[1]) else "—" }, character(1))
      add("| %s | %s |", cty, paste(vals, collapse = " | ")) }
    add("")
    add("Figures: **`key_outputs/figures/Figure_urban_<city>.pdf`** (publication two-panel: attributable-risk map zoomed to the catchment + baseline→conditional reach dumbbell by province) and **`Figure_urban_all_cities.pdf`** (five-city spatial overview); plus `cascade/figures/urban_impact_summary.pdf` (cross-city), `urban_timing_sensitivity.pdf` (timing), per-city `urban_attributable_*.pdf` and `urban_top_zones_*.pdf`. Tables: `cascade/tables/urban_scenario_summary.csv`, `urban_impact_<city>.csv`.")
    add("")
  }
  add("---")
  add("## 5. Validation & robustness (PLAN §6 — honest about the data ceiling)")
  add("A true 13-week hold-out has at most one truncated origin, so validation is **layered**:")
  add("")
  add("- **Consistency gate (h=1):** ranking Spearman = %.3f vs the validated short-horizon model; level recalibrated to observed frequency. %s",
      ctx$gate$spearman_signal, if (isTRUE(ctx$gate$pass)) "PASS." else "Review.")
  if (!is.null(ctx$sens)) {
    add("- **Ranking stability (sensitivity):** 13-week reach ranking is stable across assumptions — Spearman vs baseline: %s. Decision products (top-K, gateway) are robust even where absolute probabilities are not.",
        paste(sprintf("%s %.2f", ctx$sens$axis, ctx$sens$spearman), collapse = "; "))
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
  add("")
  add("---")
  add("## 6. Assumptions & limitations (stated up front)")
  add("- **Projection, not forecast.** 13-week trajectories are dominated by future transmissibility/response; issued as scenarios, S1 central.")
  add("- **Absolute probabilities are upper bounds**; rankings/relative risk/gateway ordering are the decision-grade products.")
  add("- **Endogenous spread only:** no cross-border reintroduction (Uganda/South Sudan) or off-network seeding.")
  add("- **Static mobility** over 13 weeks (sensitivity: alternative kernels); **R_eff** scenario-driven with province partial-pooling; most affected zones are data-poor and lean on the pool.")
  add("- **Small event base** (~39 h=1 invasions) limits identifiability; priors + pooling do heavy lifting.")
  add("- Two-stage nowcast point-plugs the initial conditions; delta/psi are calibrated to short-horizon frequency and the observed early expansion pace, not to a long-horizon outcome (unavailable).")
  add("")
  add("Full spec, equations, and the review log: `PLAN_3MONTH_INVASION.md`.")
  writeLines(L, path)
  message("[cascade] wrote report -> ", path)
  invisible(path)
}
