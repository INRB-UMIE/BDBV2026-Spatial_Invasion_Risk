# =============================================================================
# 33_cascade_eval.R — 3-MONTH CASCADE: validation & calibration (PLAN §6, §3.5)
# BDBV 2026 DRC
#
# Layered validation, honest about the ~15-week data ceiling (a true 13-week
# hold-out has at most one truncated origin):
#   cascade_consistency_gate()  h=1 RANKING vs the validated short-horizon model (§6.4),
#                               plus a level ratio that is ~delta by construction and is
#                               therefore only a plumbing check, not calibration.
#   cascade_fit_psi()           RETIRED naive-pace psi diagnostic (CASCADE_RUN_PACE_DIAG=1).
#   cascade_backtest()          intermediate-horizon (<=8w) leave-future-out cascade
#                               backtest at the earliest origins (~1 effective origin;
#                               indicative), scored on reach discrimination + the new-
#                               invasion count — the COUNT is in-sample whenever delta was
#                               fitted on the same origins and K (run_cascade.R reuses them).
#   cascade_sensitivity()       ranking-stability across GT/kernel/psi/k/scenario (§6.6).
#   cascade_sbc()               light simulation-based calibration of the kernel (§6.5).
# Small self-contained metric helpers (no coupling to 16_invasion_eval internals).
# =============================================================================

# ---- metric helpers --------------------------------------------------------
# The published top-K budgets, shared with evaluate_invasion() (16_invasion_eval.R).
.TOPK <- get0("TOPK_VALUES", ifnotfound = c(5L, 15L, 20L))

# ONE average-precision implementation for the whole suite. This file used to carry its own
# tie-aware copy under a comment saying it had "no coupling to 16_invasion_eval internals";
# the two were numerically identical, so the independence bought nothing and cost a second
# definition that could drift. 16_invasion_eval.R is loaded before this file by every entry
# point (run_cascade.R, tests/run_tests.R, 41_kindiv_sweep.R); the fallback keeps an isolated
# source() working and warns rather than silently scoring with a different estimator.
.casc_auc_pr <- function(p, y) {
  if (exists(".auc_pr", mode = "function")) return(.auc_pr(p, y))
  warning("[cascade] .auc_pr() unavailable (16_invasion_eval.R not loaded); the cascade ",
          "backtest is scored with a local fallback.", call. = FALSE)
  ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- as.integer(y[ok])
  P <- sum(y == 1); if (!P || !length(p)) return(NA_real_)
  o <- order(p, decreasing = TRUE); ps <- p[o]; y <- y[o]
  tp <- cumsum(y == 1); fp <- cumsum(y == 0)
  keep <- c(ps[-length(ps)] != ps[-1], TRUE)
  sum((tp / (tp + fp))[keep] * diff(c(0, (tp / P)[keep])), na.rm = TRUE)
}
.casc_topk <- function(p, y, K) {             # precision@K (fraction of top-K truly invaded)
  ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  if (!length(p)) return(NA_real_)
  # TIE-AWARE. This used head(order(p), K), which breaks ties by ROW ORDER — so the
  # cascade's headline precision@K depended on the order zones happened to arrive in.
  # .casc_auc_pr above is already tie-aware (it collapses equal-score runs), so the two
  # metrics beside each other disagreed about what a tie means. rank(-p, "max") <= K
  # credits a zone only if monitoring K zones necessarily includes it — the same
  # convention as .ranking_metrics (16), spatiotemporal_skill (19),
  # prospective_invasion_check (23) and lead_time_analysis (19). Denominator stays the
  # REALISED top-K size; an all-tied set yields 0, correctly scoring a non-discriminating
  # forecast rather than handing it a row-order-dependent number.
  rk <- rank(-p, ties.method = "max")
  inK <- rk <= K
  nk <- sum(inK)
  if (!nk) return(0)
  sum(y[inK] == 1) / nk
}
.casc_spearman <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman",
                                                             use = "pairwise.complete.obs"))

# ---------------------------------------------------------------------------
# CONSISTENCY GATE (PLAN §6.4). At h=1 the cascade uses no source projection, so
# with sat==1 at week 1 its per-zone reach must be PROPORTIONAL to the validated
# short-horizon model's week-1 probability from the same fit (the factor being delta).
# The h=2 leg is compared to the DETERMINISTIC mean-field to avoid the Jensen
# offset (PLAN §6.4); here we anchor on the h=1 leg.
#
# THE PASS RULE IS A NOISE-CALIBRATED DISPERSION TEST, not a Spearman threshold. A threshold on
# the rank correlation of the few zones above `signal_thresh` cannot be met reliably: as the
# outbreak declines those zones are near-tied (e.g. 6 zones at 2.2-3.7%), and the cascade's
# Monte-Carlo noise at n_mc = 3000 (~15-20% of each zone's probability) reorders them. A
# perfectly implemented cascade then fails rho >= 0.9 most of the time (measured on the
# 2026-10-16 frame: 84%). The test below asks instead whether every zone's cascade reach
# departs from level x prediction by more than its OWN Monte-Carlo error allows: Pearson
# X^2 = sum (p_casc - level*p_pred)^2 / V over the zones with >= 5 expected invasions, on k - 1
# df (one fitted level). V is the MODEL-BASED variance of the cascade's mean over D parameter
# groups of R replicates each, V = [(q(1-q) - s^2)/R + s^2] / D, with q = level*p_pred and
# s = level*p_sd (the forecast's posterior SD, scaled like the probability): process noise
# within a group plus posterior-sampling noise between groups. (The empirical between-group
# variance was tried first: for the many low-probability zones only a handful of groups record
# an invasion, the estimate is noisy, and dividing by it inflates X^2.) The pass rule is
# p >= gate_alpha. Spearman rho on the signal zones is still reported, as a diagnostic.
#
# THE GATE SIMULATES ON THE SHORT-HORIZON MODEL'S OWN GENERATION TIME. When the cascade
# marginalises the generation time over a grid (CASCADE_GT_MARGINALISE = 1, a sensitivity; OFF
# by default, where prep$G_list is NULL and nothing below changes), the grid point a group draws
# sets the week-1 import force (own = sum_k G[k] I[, t-k]) while the short-horizon model uses the
# single CASCADE_GT kernel, so the two would differ BY DESIGN, zone by zone. The week-1
# simulation here therefore always holds the generation time at that kernel: what is tested is
# that the simulator reproduces the model it is built from.
#
# WHAT THE GATE CANNOT TEST: THE LEVEL. Both sides are built from the SAME fit, so the
# cascade's week-1 probability is 1 - exp(-delta * mu) against the short-horizon model's
# 1 - exp(-mu), and their ratio is ~delta for ANY delta — dividing it by delta returns ~1 by
# construction. It used to be part of `pass` and was described as a calibration test ("level
# recalibrated to observed frequency"); it is an identity, and can fail only if the plumbing is
# broken (delta not applied, saturation not 1 at week 1). It is kept as exactly that — an
# implementation check, reported and warned on, but not part of `pass`. The one-week LEVEL is
# tested against OBSERVED week-1 invasions by cascade_delta_h1_refit() (33b), at held-out
# origins.
# ---------------------------------------------------------------------------
cascade_consistency_gate <- function(prep, fit, design, delta,
                                     zone_week_nc, mobility_matrices, gt_pmfs,
                                     covariates, osrm_mat, zones_all,
                                     n_mc = 3000L, signal_thresh = 0.02,
                                     gate_alpha = 0.01, min_expected = 5,
                                     # Tolerance on the PLUMBING check level_ratio / delta (see
                                     # the header): it is ~1 by construction for any delta, so a
                                     # departure means the factor is not reaching the hazard.
                                     level_tol = 0.25) {
  scP <- CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]
  # No pop_vec: within-zone depletion is not modelled (removed 2026-09-19), so simulate_cascade()
  # accepts the argument but uses it in no simulation of the suite.
  prep1 <- prep; prep1$G_list <- NULL; prep1$gt_weights <- NULL; prep1$gt_keys <- NULL
  sim1 <- simulate_cascade(prep1, scP, n_mc = n_mc, horizon = 1L, delta = delta,
                           psi = CASCADE_PSI)
  reach1 <- cascade_reach_table(sim1, horizons = 1L)
  offs <- bayes_forecast_offsets(zone_week_nc, mobility_matrices, gt_pmfs, covariates,
            osrm_mat, zones_all, mob = prep$kernel, gt = CASCADE_GT,
            cov_spec = CASCADE_COV_SPEC, horizons = 1L)
  affected <- zones_all[prep$affected0]
  # The short-horizon model's own week-1 prediction from the SAME fit, by its own engine: a
  # featured time-varying form is predicted by 21c's forecast (drift and frailty included).
  pred <- if (inherits(fit, "stan_invasion_fit"))
            .finalise_invasion_summary(stan_invasion_forecast_draws(fit, offs, 1L), affected)
          else predict_bayes_invasion(fit, offs, design, horizons = 1L, affected_zones = affected)
  if (!"p_sd" %in% names(pred))
    stop("[cascade_gate] the short-horizon forecast carries no posterior SD (p_sd); the gate's ",
         "Monte-Carlo variance cannot be formed.", call. = FALSE)
  m <- merge(reach1[, c("health_zone", "p_invasion")],
             pred[, c("health_zone", "p_invasion", "p_sd")], by = "health_zone",
             suffixes = c("_casc", "_pred"))
  m <- m[is.finite(m$p_invasion_casc) & is.finite(m$p_invasion_pred), ]
  sig <- m[m$p_invasion_pred > signal_thresh, ]
  rho <- .casc_spearman(sig$p_invasion_casc, sig$p_invasion_pred)
  lvl <- sum(m$p_invasion_casc) / sum(m$p_invasion_pred)
  # Monte-Carlo variance of each zone's week-1 reach, model-based (see the header).
  .D <- as.numeric(sim1$n_draws_used); .R <- as.numeric(sim1$n_rep)
  .q <- lvl * m$p_invasion_pred; .s2 <- (lvl * m$p_sd)^2
  m$var_casc <- (pmax(.q * (1 - .q) - .s2, 0) / .R + .s2) / .D
  .k <- is.finite(m$var_casc) & m$var_casc > 0 & lvl * m$p_invasion_pred * sim1$n_mc >= min_expected
  x2 <- sum((m$p_invasion_casc[.k] - lvl * m$p_invasion_pred[.k])^2 / m$var_casc[.k])
  df <- sum(.k) - 1L
  p_disp <- if (df >= 1L) stats::pchisq(x2, df, lower.tail = FALSE) else NA_real_
  pass_rank  <- is.finite(p_disp) && p_disp >= gate_alpha
  lvl_rel    <- if (is.finite(lvl) && is.finite(delta) && delta > 0) lvl / delta else NA_real_
  # Plumbing check only (see the header): NOT part of `pass`, because it cannot fail for a
  # wrong delta — only for a delta that never reached the hazard.
  pass_level <- is.finite(lvl_rel) && abs(lvl_rel - 1) <= level_tol
  pass <- pass_rank
  if (!pass_level)
    warning(sprintf(paste0("[cascade_gate] h=1 PLUMBING check failed: the cascade seeds at %.3f x the ",
                           "raw short-horizon model from the same fit, where applying delta=%.3f ",
                           "should give ~delta; ratio/delta = %.2f (tolerance +/-%.0f%%). The factor ",
                           "is not reaching the simulated hazard as intended. (This is not a ",
                           "calibration test: it is ~1 for any delta that IS applied.)"),
                   lvl, delta, lvl_rel, 100 * level_tol), call. = FALSE)
  # keep the paired point cloud (for the evaluation figure's week-1 anchor panel):
  # signal = zones the validated model flags above threshold (drive the ranking test).
  gate_data <- data.frame(health_zone = m$health_zone,
                          p_invasion_pred = m$p_invasion_pred,
                          p_invasion_casc = m$p_invasion_casc,
                          signal = m$p_invasion_pred > signal_thresh)
  list(spearman_signal = rho, n_signal = nrow(sig), level_ratio = lvl,
       x2 = x2, x2_df = df, x2_p = p_disp, gate_alpha = gate_alpha, n_tested = sum(.k),
       level_vs_delta = lvl_rel, delta = delta,
       # pass_level is the PLUMBING check (renamed in the note); `pass` is the ranking gate.
       pass = pass, pass_rank = pass_rank, pass_level = pass_level,
       level_is_identity = TRUE,
       level_tol = level_tol,
       signal_thresh = signal_thresh, data = gate_data,
       note = sprintf(paste0("h1 proportionality X2=%.1f on %d df (p=%.3f, %d zones, gate p>=%.2f) %s ",
                             "[rank diagnostic: rho=%.3f on %d signal zones]; plumbing: level=%.3f vs ",
                             "delta=%.3f (ratio %.2f, ~1 by construction, tol +/-%.0f%%) %s -> %s"),
                      x2, df, p_disp, sum(.k), gate_alpha, if (pass_rank) "PASS" else "FAIL",
                      rho, nrow(sig),
                      lvl, delta, lvl_rel, 100 * level_tol,
                      if (pass_level) "ok" else "BROKEN",
                      if (pass) "PASS" else "FAIL"))
}

# ---------------------------------------------------------------------------
# RETIRED AS AN ESTIMATOR (2026-09-11) — retained as a naive-pace DIAGNOSTIC.
#
# This fits psi so the cascade's expected cumulative new-invasion count at
# `match_week` matches "mean weekly new zones over the last `ref_weeks`" x
# match_week. That target is a CONSTANT-RATE EXTRAPOLATION of the recent observed
# pace, i.e. an assumption about the very dynamics the 13-week layer exists to
# project, and it is biased LOW because `trunc_drop = 1` removes only one week for
# reporting truncation while the onset->sample delay alone has a mean near 8 days
# and a 99th percentile near 6 weeks. On the 2026-09-07 run it returned psi = 0,
# CENSORED at the lower bound with converged = FALSE — psi was not identified at
# all — while the genuinely out-of-sample backtest showed the cascade over-running
# the realised counts by ~1.8x at 6 weeks.
#
# psi is now fitted by cascade_fit_psi_oos() (33b_cascade_calibration.R) against
# the invasions that ACTUALLY occurred in held-out windows. This function is kept
# because the modelled-vs-naive-pace comparison is still a useful sanity number,
# and it is called only when CASCADE_RUN_PACE_DIAG=1. Do not restore it as the
# estimator.
# ---------------------------------------------------------------------------
cascade_fit_psi <- function(prep, delta, zone_week_nc, zones_all,
                            psi_lo = 0, psi_hi = CASCADE_PSI_MAX, n_coarse = 6L,
                            ref_weeks = 6L, n_mc = 600L, match_week = 8L,
                            trunc_drop = 1L, tol = CASCADE_PSI_TOL, max_refine = 12L,
                            report_horizon = CASCADE_HORIZON_WEEKS, pop_vec = NULL) {
  # observed weekly new-zone rate over the reference window (first-case weeks). Drop
  # the final `trunc_drop` week(s): even after nowcast correction the most recent week
  # is the most right-truncated, so newly-invaded zones there are under-observed, which
  # would bias obs_rate low -> target_cum low -> psi_hat high (over-damping the front).
  Y <- .count_wide(zone_week_nc, zones_all, "confirmed_nc")
  first_wk <- apply(Y, 1, function(x) { w <- which(x > 0); if (length(w)) min(w) else NA_integer_ })
  nT <- ncol(Y) - trunc_drop; recent <- (nT - ref_weeks + 1L):nT
  obs_rate <- mean(tabulate(first_wk[first_wk %in% recent] - (nT - ref_weeks), nbins = ref_weeks))
  target_cum <- obs_rate * match_week          # naive observed-pace expectation at match_week
  scP <- CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]

  # Expected cumulative new zones by `week` at a given psi. simulate_cascade() calls
  # set.seed(seed) on entry, so every evaluation shares one random stream (common random
  # numbers) and the curve is a smooth, monotone-decreasing function of psi rather than a
  # noisy one — which is what makes bracketing and bisection valid here.
  trace <- list()
  f_cum <- function(ps, week = match_week) {
    s <- simulate_cascade(prep, scP, n_mc = n_mc, horizon = week, delta = delta, psi = ps,
                          pop_vec = pop_vec)
    v <- sum(rowMeans(s$new_by_week))
    if (week == match_week) trace[[length(trace) + 1L]] <<- data.frame(psi = ps, modelled_cum = v)
    v
  }

  # ---- stage 1: coarse bracket (0 plus a log-spaced ladder to psi_hi) -------
  ladder <- if (n_coarse > 1L)
    exp(seq(log(0.25), log(max(psi_hi, 0.5)), length.out = n_coarse - 1L)) else numeric(0)
  grid   <- sort(unique(c(psi_lo, ladder, psi_hi)))
  curve  <- vapply(grid, f_cum, numeric(1))
  gap    <- curve - target_cum                 # > 0 = cascade still faster than observed pace

  boundary_hit <- FALSE; boundary_side <- NA_character_
  converged <- TRUE; psi_hat <- NA_real_
  if (all(gap <= 0)) {
    # Even with NO saturation the cascade is at or below the observed pace. Saturation is
    # one-sided (it can only slow the front), so psi_lo is the constrained optimum — but the
    # target is still unmet, and the mismatch lives in delta / the hazard, not in psi. Report
    # it as a LOWER-boundary solution so it is not confused with the over-running case.
    psi_hat <- grid[1L]
    if (abs(gap[1L]) > tol * max(target_cum, 1e-9)) {
      converged <- FALSE; boundary_hit <- TRUE; boundary_side <- "lower"
      warning(sprintf(paste0("[cascade_fit_psi] at psi=%.2f (no saturation) the cascade already ",
                             "projects only %.1f new zones by week %d vs an observed-pace target of ",
                             "%.1f (%.0f%% under). psi cannot correct an UNDER-running front — the ",
                             "miss is in delta or the hazard, not the saturation term."),
                     psi_hat, curve[1L], match_week, target_cum,
                     100 * (1 - curve[1L] / max(target_cum, 1e-9))), call. = FALSE)
    }
  } else if (all(gap > 0)) {
    # still over-running at the top of the search interval -> CENSORED, not fitted.
    psi_hat <- grid[length(grid)]
    boundary_hit <- TRUE; boundary_side <- "upper"; converged <- FALSE
    warning(sprintf(paste0("[cascade_fit_psi] frontier saturation did not converge: at psi=%.1f ",
                           "(search maximum) the cascade still projects %.1f new zones by week %d ",
                           "vs an observed-pace target of %.1f (%.0f%% over). psi is CENSORED at the ",
                           "boundary and the projection will over-run. Check the mobility kernel's ",
                           "dispersion before using the long-horizon layer."),
                   psi_hat, curve[length(curve)], match_week, target_cum,
                   100 * (curve[length(curve)] / max(target_cum, 1e-9) - 1)), call. = FALSE)
  } else {
    # ---- stage 2: bisect the bracketing interval ---------------------------
    i_hi <- max(which(gap > 0))                # last psi still above target
    lo <- grid[i_hi]; hi <- grid[i_hi + 1L]
    g_lo <- gap[i_hi]; g_hi <- gap[i_hi + 1L]
    for (it in seq_len(max_refine)) {
      if (abs(g_lo) <= tol * max(target_cum, 1e-9)) { hi <- lo; g_hi <- g_lo; break }
      if (abs(g_hi) <= tol * max(target_cum, 1e-9)) break
      mid <- 0.5 * (lo + hi)
      g_mid <- f_cum(mid) - target_cum
      if (g_mid > 0) { lo <- mid; g_lo <- g_mid } else { hi <- mid; g_hi <- g_mid }
    }
    # take the endpoint closer to the target
    psi_hat <- if (abs(g_lo) <= abs(g_hi)) lo else hi
    converged <- min(abs(g_lo), abs(g_hi)) <= tol * max(target_cum, 1e-9)
    if (!converged)
      warning(sprintf(paste0("[cascade_fit_psi] bisection exhausted %d refinements without ",
                             "reaching the %.0f%% tolerance; residual %.1f zones on a target of %.1f."),
                      max_refine, 100 * tol, min(abs(g_lo), abs(g_hi)), target_cum), call. = FALSE)
  }

  # ---- runaway diagnostic ---------------------------------------------------
  # psi is anchored at `match_week`; NOTHING constrains the trajectory beyond it. Report
  # how far the fitted cascade runs past a naive continuation of the observed frontier
  # pace over the full reporting horizon, so a compounding blow-up is visible in the log
  # and the provenance sidecar rather than only in the final projected counts.
  horizon_cum   <- f_cum(psi_hat, week = report_horizon)
  horizon_naive <- obs_rate * report_horizon
  horizon_ratio <- horizon_cum / max(horizon_naive, 1e-9)

  curve_df <- do.call(rbind, trace)
  curve_df <- curve_df[order(curve_df$psi), , drop = FALSE]
  curve_df <- curve_df[!duplicated(curve_df$psi), , drop = FALSE]
  curve_df$target_cum <- target_cum
  curve_df$gap        <- curve_df$modelled_cum - target_cum
  # Evaluate the chosen psi ONCE: every f_cum() is a full simulation, and calling it
  # separately for the returned scalar and again for the note both doubled the cost and
  # reported two different numbers for the same quantity (the streams desynchronise).
  .near <- which.min(abs(curve_df$psi - psi_hat))
  modelled_at_psi <- if (abs(curve_df$psi[.near] - psi_hat) < 1e-9)
    curve_df$modelled_cum[.near] else f_cum(psi_hat)

  list(psi = psi_hat, converged = converged, boundary_hit = boundary_hit,
       boundary_side = boundary_side, psi_lo = psi_lo, psi_hi = psi_hi, tol = tol,
       curve = curve_df, modelled_cum_at_psi = modelled_at_psi,
       observed_pace_cum = target_cum, obs_weekly_rate = obs_rate, match_week = match_week,
       report_horizon = report_horizon, horizon_cum = horizon_cum,
       horizon_naive_cum = horizon_naive, horizon_pace_ratio = horizon_ratio,
       note = sprintf("psi=%.2f (%s%s); modelled %.1f vs observed-pace %.1f new zones by wk %d; at wk %d the cascade projects %.1f vs %.1f naive-pace (%.1fx)",
                      psi_hat, if (converged) "converged" else "NOT converged",
                      if (boundary_hit) sprintf(", CENSORED at the %s search bound", boundary_side) else "",
                      modelled_at_psi, target_cum, match_week,
                      report_horizon, horizon_cum, horizon_naive, horizon_ratio))
}
# ---------------------------------------------------------------------------
# Out-of-sample backtest at the cascade's own horizon
# ---------------------------------------------------------------------------

#' Score the cascade against realised multi-week outcomes from held-out origins.
#'
#' Origin construction (as-of training counts, per-origin hazard/R_eff refit, the
#' evaluation zone set and the realised outcomes) is delegated to
#' cascade_calibration_origins() (33b_cascade_calibration.R) so the backtest and the
#' CALIBRATION score the same object. It used to build its own origins with a
#' near-identical block; the two copies could drift, and one of them was wrong: the
#' evaluation set came from the model's as-of at-risk mask while the outcome came
#' from the final record, so a zone already invaded at the cutoff but not yet
#' REPORTED was scored as a false alarm. That charges the model for a reporting lag
#' and inflates the apparent over-prediction. The shared builder intersects the two
#' at-risk definitions, as run_invasion_lfo() does.
#'
#' @param origins optional pre-built origins (e.g. from cascade_calibrate()); pass
#'   them to avoid refitting the same per-origin hazards a second time.
#' @return tibble, one row per origin, with attr("detail") carrying the per-zone rows.
cascade_backtest <- function(layer, cutoffs_from_end = c(10L, 9L, 8L), K = 6L,
                             n_mc = 300L, delta = NULL, psi = CASCADE_PSI,
                             kernel = CASCADE_KERNEL, zone_province = NULL,
                             # SETTLED outcomes only (LFO_RELIABLE_EVAL_AGE_DAYS), the same
                             # threshold the calibration origins use: was 10 d, below the p90
                             # onset->confirmation delay.
                             linelist = NULL,
                             min_eval_age_days = get0("LFO_RELIABLE_EVAL_AGE_DAYS", ifnotfound = 25L),
                             analysis_date = get0("ANALYSIS_DATE", ifnotfound = NA),
                             origins = NULL, seed = CASCADE_SEED) {
  K <- as.integer(K)
  if (is.null(origins)) {
    if (!exists("cascade_calibration_origins"))
      stop("[cascade_backtest] cascade_calibration_origins() not found — source 33b_cascade_calibration.R.")
    origins <- cascade_calibration_origins(
      layer, zone_province = zone_province, cutoffs_from_end = cutoffs_from_end,
      truth_horizons = K, min_eval_age_days = min_eval_age_days, kernel = kernel,
      gt = CASCADE_GT, linelist = linelist, analysis_date = analysis_date)
  }
  if (!length(origins)) return(NULL)
  key <- as.character(K)
  # Origins are gated PER HORIZON, so a recent one may carry its 1-week outcome but not
  # its K-week one. Score the origins that do, rather than refusing to score any.
  origins <- Filter(function(o) key %in% names(o$truth), origins)
  if (!length(origins)) {
    warning(sprintf("[cascade_backtest] no origin carries a settled %d-week outcome (min_eval_age_days).", K),
            call. = FALSE)
    return(NULL)
  }

  del <- delta %||% cascade_fit_delta(cov_model = cascade_hazard_label(kernel))
  scP <- CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]
  res <- list(); det <- list()
  for (o in origins) {
    # Same model as production (within-zone depletion is modelled in neither: removed
    # 2026-09-19).
    sim   <- simulate_cascade(o$prep, scP, n_mc = n_mc, horizon = K,
                              delta = as.numeric(del), psi = psi, seed = seed,
                              pop_vec = o$pop_vec)
    reach <- cascade_reach_table(sim, horizons = K)
    p <- reach$p_invasion[match(o$eval_zones, reach$health_zone)]
    y <- o$truth[[key]]
    ok <- is.finite(p) & !is.na(y)
    p <- p[ok]; y <- y[ok]; z <- o$eval_zones[ok]
    if (!length(p)) next
    br <- mean(y); ap <- .casc_auc_pr(p, y)
    rk_avg <- rank(-p, ties.method = "average")
    .row <- tibble::tibble(
      cutoff = as.character(o$cutoff_date), K = K, n_atrisk = length(p), n_events = sum(y),
      base_rate = br, auc_pr = ap, auc_pr_skill = ap / max(br, 1e-9),
      mean_rank_of_truth = if (any(y == 1L)) mean(rk_avg[y == 1L]) else NA_real_,
      pred_new = sum(p), obs_new = sum(y), count_ratio = sum(p) / max(sum(y), 1e-9))
    # SAME top-K budgets as the short-horizon evaluation (TOPK_VALUES, 00_config.R), so a
    # cascade precision@K and an invasion precision@K mean the same watch-list size. The
    # column NAME is generated from the budget it holds: the three names used to be written
    # out by hand beside values taken from .TOPK[1:3], so changing TOPK_VALUES left a column
    # called `prec_at15` carrying precision at some other K -- and .TOPK[3] silently became
    # NA if the vector was shortened.
    for (.K in .TOPK) .row[[paste0("prec_at", .K)]] <- .casc_topk(p, y, .K)
    res[[length(res) + 1L]] <- .row
    det[[length(det) + 1L]] <- tibble::tibble(
      cutoff = as.character(o$cutoff_date), K = K,
      health_zone = z, p_invasion = p, y = as.integer(y))
  }
  if (!length(res)) return(NULL)
  out <- dplyr::bind_rows(res)
  attr(out, "detail") <- dplyr::bind_rows(det)   # survives saveRDS / direct passing
  out
}


# ---------------------------------------------------------------------------
# Sensitivity / ranking-stability (PLAN §6.6). Re-run the primary scenario under
# perturbed assumptions and report the Spearman of the 13-week reach ranking vs
# the baseline. Acceptance: decision products (ranking) stable even where absolute
# probabilities are not. Uses a modest M for speed.
# ---------------------------------------------------------------------------
cascade_sensitivity <- function(layer, design_by_kernel, fit_by_kernel, reff, delta,
                                base_reach = NULL, n_mc = 400L, psi = CASCADE_PSI) {
  # MATCHED BASELINE. Each alternative below is simulated at `n_mc` with a fixed delta, so
  # the comparison must be against a baseline run the same way. Passing the PRODUCTION
  # reach table here (M = 2000, and since 2026-09-11 with delta drawn per parameter group)
  # contrasted a large, calibration-averaged baseline against small, fixed-delta
  # alternatives and charged the difference to the sensitivity axis. The baseline is now
  # simulated here under identical settings; `base_reach` is accepted only as a fallback.
  scP <- CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]
  prep_base <- cascade_prepare(layer, fit_by_kernel[[CASCADE_KERNEL]],
                               design_by_kernel[[CASCADE_KERNEL]], reff, kernel = CASCADE_KERNEL)
  base_matched <- tryCatch(
    cascade_reach_table(simulate_cascade(prep_base, scP, n_mc = n_mc, delta = delta, psi = psi,
                                        pop_vec = layer$pop)),
    error = function(e) NULL)
  if (is.null(base_matched)) {
    if (is.null(base_reach))
      stop("[cascade_sensitivity] could not simulate a matched baseline and none was supplied.")
    warning("[cascade_sensitivity] matched baseline failed; falling back to the supplied ",
            "production reach, whose Monte-Carlo size and delta treatment differ.", call. = FALSE)
    base_matched <- base_reach
  }
  base_r <- base_matched[base_matched$horizon == max(CASCADE_REPORT_HORIZONS) &
                         !base_matched$was_active_before, c("health_zone", "p_invasion")]
  runs <- list()
  # `level_comparable` says whether the row's mean reach can be read against the baseline's.
  # It cannot for the KERNEL rows: each alternative kernel's hazard is refit on its own design,
  # so its intercept — and the scale delta corrects — is that kernel's, while delta was fitted
  # (on the held-out-origin calibration windows) for the FEATURED kernel only. Running the alternative at the featured delta
  # puts its level on an uncalibrated scale (the shipped table showed +80-105% mean reach on
  # the kernel axis, an artefact of exactly this). Refitting delta per kernel would need the
  # whole held-out-origin calibration once per kernel; instead those rows report the RANKING
  # comparison only and their mean reach is withheld (NA), with the reason recorded.
  add <- function(tag, reach, level_comparable = TRUE, note = NA_character_) {
    r <- reach[reach$horizon == max(CASCADE_REPORT_HORIZONS) & !reach$was_active_before,
               c("health_zone", "p_invasion")]
    m <- merge(base_r, r, by = "health_zone", suffixes = c("_base", "_alt"))
    runs[[tag]] <<- tibble::tibble(axis = tag,
      spearman = .casc_spearman(m$p_invasion_base, m$p_invasion_alt),
      mean_reach = if (isTRUE(level_comparable)) mean(r$p_invasion, na.rm = TRUE) else NA_real_,
      level_comparable = isTRUE(level_comparable), note = note)
  }
  # The matched baseline itself is reported as a row so the mean reach every alternative
  # is compared against is visible in the table rather than implicit (its Spearman is 1
  # by construction).
  add(sprintf("baseline (psi=%.2f)", psi), base_matched)
  # NOTE: the individual-offspring overdispersion (k_indiv) axis is handled by the
  # dedicated cascade_kindiv_sweep() (dense grid, at the deployed psi), not here, to avoid
  # duplicate work and keep a single canonical k result.
  # psi sweep. prep for the production kernel is built ONCE above; it was previously
  # rebuilt inside this loop on every psi, re-extracting the posterior draws and the
  # 519x519 travel-time matrix for no reason.
  for (ps in CASCADE_PSI_SWEEP) {
    add(sprintf("psi=%.1f", ps),
        cascade_reach_table(simulate_cascade(prep_base, scP, n_mc = n_mc, delta = delta, psi = ps,
                                             pop_vec = layer$pop)))
  }
  # R-WALK VOLATILITY axis. CASCADE_R_RW_SWEEP arrived with the temporal R model; without this
  # loop it would be a constant nothing reads — a sweep that looks available and is not. The
  # baseline row already covers the configured sigma, so it is excluded here; the sigma = 0 arm
  # is the constant-R model and therefore doubles as the backward-compatibility check.
  .rw_base <- get0("CASCADE_R_RW_SIGMA", ifnotfound = 0)
  for (sg in setdiff(get0("CASCADE_R_RW_SWEEP", ifnotfound = numeric(0)), .rw_base)) {
    add(sprintf("R-walk sigma=%.2f", sg),
        cascade_reach_table(simulate_cascade(prep_base, scP, n_mc = n_mc, delta = delta,
                                             psi = psi, rw_sigma = sg, pop_vec = layer$pop)))
  }
  # kernel sweep (needs its own fit/design)
  for (kern in setdiff(CASCADE_SENS_KERNELS, CASCADE_KERNEL)) {
    if (is.null(fit_by_kernel[[kern]])) next
    prepk <- cascade_prepare(layer, fit_by_kernel[[kern]], design_by_kernel[[kern]],
                             reff, kernel = kern)
    add(sprintf("kernel=%s", kern),
        cascade_reach_table(simulate_cascade(prepk, scP, n_mc = n_mc, delta = delta, psi = psi,
                                             pop_vec = layer$pop)),
        level_comparable = FALSE,
        note = sprintf(paste0("rank comparison only: delta was fitted for %s, not for %s, so ",
                              "this kernel's absolute reach is uncalibrated"), CASCADE_KERNEL, kern))
  }
  dplyr::bind_rows(runs)
}

# ---------------------------------------------------------------------------
# COMPREHENSIVE individual-offspring overdispersion (k_indiv) sweep (PLAN §6.6).
# Re-runs the PRIMARY scenario cascade over a dense k grid spanning heavy
# superspreading (k~0.05) to the near-Poisson limit (k~4), densely sampling the
# Lloyd-Smith (2005) EBOV-plausible band [0.2, 0.4]. Reports, per k: ranking
# stability vs the k_base baseline (Spearman, Kendall, top-15 Jaccard) and absolute
# magnitude (mean reach at each report horizon + mean establishment at Hmax).
#
# All runs — baseline AND alternatives — use the SAME (delta, psi) as production, so
# the effect attributed to k is purely k. (Contrast cascade_sensitivity(), whose k
# rows historically ran at the default psi while the baseline used the then-fitted psi;
# that conflation is avoided here and the call site passes the deployed psi, now fixed at 0.)
# Returns list(summary, byzone). `province` is a named vector zone -> province.
# ---------------------------------------------------------------------------
cascade_kindiv_sweep <- function(prep, scenario, delta, psi, pop_vec = NULL,
                                 province = NULL,
                                 k_grid = CASCADE_K_INDIV_SWEEP, k_base = CASCADE_K_INDIV,
                                 n_mc = CASCADE_N_MC, horizons = CASCADE_REPORT_HORIZONS,
                                 topN = 15L) {
  Hmax <- max(horizons)
  key  <- function(kk) sprintf("%.4g", kk)
  grid <- k_grid; if (!any(abs(grid - k_base) < 1e-9)) grid <- sort(c(grid, k_base))
  reach_by_k <- list()
  for (kk in grid) {
    sim <- simulate_cascade(prep, scenario, n_mc = n_mc, delta = delta, psi = psi,
                            k_indiv = kk, pop_vec = pop_vec)
    reach_by_k[[key(kk)]] <- cascade_reach_table(sim, horizons = horizons)
  }
  atrisk <- function(rt, H) {
    d <- rt[rt$horizon == H & !rt$was_active_before,
            c("health_zone", "p_invasion", "p_establishment")]
    d[is.finite(d$p_invasion), ]
  }
  # TIE-AWARE top-N sets. Both sets were taken with order(-p), and both data frames come
  # from atrisk() in the SAME zones_all row order — so ties were broken IDENTICALLY in base
  # and alt, and the Jaccard counted agreement that is an artefact of shared row ordering
  # rather than shared signal. That matters here more than anywhere else: p_invasion is a
  # Monte-Carlo reach probability with granularity 1/n_mc (0.0025 at the sensitivity default
  # n_mc = 400), so exact ties are COMMON and can easily straddle the top-15 boundary — and
  # top15_jaccard is the headline "the ranking is invariant to k_indiv" statistic. Using
  # rank(-p, "max") <= topN excludes a straddling tie group from BOTH sets, so the Jaccard
  # measures only genuinely separable agreement. (Spearman/Kendall beside it already handle
  # ties correctly via average ranks, which is why they were never affected.)
  .topset <- function(d, N) d$health_zone[rank(-d$p_invasion, ties.method = "max") <= N]
  # THE COLUMN NAMES BELOW ENCODE THEIR OWN ARGUMENTS. `mean_reach_h4/h8/h13`,
  # `mean_estab_h13` and `top15_jaccard` are written out as literals while their VALUES come
  # from `horizons` and `topN`, and 36_report.R and 42_kindiv_sweep_figure.R read those
  # literals and label them "13 weeks" / "Top-15 Jaccard" in published prose and figures. A
  # caller that changed either argument would therefore publish a number under the wrong K or
  # the wrong horizon, silently. Stop instead: the names are part of the published contract,
  # so a genuine change to them has to be made deliberately in all three files at once.
  .hz <- as.integer(horizons)
  if (!identical(.hz[c(1L, min(2L, length(.hz)))], c(4L, 8L)) || Hmax != 13L || topN != 15L)
    stop(sprintf(paste0("[kindiv] cascade_kindiv_sweep() writes the fixed column names ",
                        "mean_reach_h4/h8/h13, mean_estab_h13 and top15_jaccard, but was ",
                        "called with horizons = %s and topN = %d. Rename the columns here ",
                        "and in 36_report.R / 42_kindiv_sweep_figure.R before changing these."),
                 paste(.hz, collapse = "/"), topN), call. = FALSE)
  base_r   <- atrisk(reach_by_k[[key(k_base)]], Hmax)
  base_top <- .topset(base_r, topN)
  jacc <- function(a, b) { u <- length(union(a, b)); if (!u) return(NA_real_)
                           length(intersect(a, b)) / u }
  summary <- dplyr::bind_rows(lapply(grid, function(kk) {
    alt <- atrisk(reach_by_k[[key(kk)]], Hmax)
    m <- merge(base_r[c("health_zone", "p_invasion")], alt[c("health_zone", "p_invasion")],
               by = "health_zone", suffixes = c("_base", "_alt"))
    alt_top  <- .topset(alt, topN)
    reach_at <- function(H) mean(atrisk(reach_by_k[[key(kk)]], H)$p_invasion, na.rm = TRUE)
    tibble::tibble(
      k_indiv       = kk,
      spearman      = suppressWarnings(stats::cor(m$p_invasion_base, m$p_invasion_alt,
                                                  method = "spearman", use = "complete.obs")),
      kendall       = suppressWarnings(stats::cor(m$p_invasion_base, m$p_invasion_alt,
                                                  method = "kendall", use = "complete.obs")),
      top15_jaccard = jacc(base_top, alt_top),
      mean_reach_h4  = reach_at(horizons[1]),
      mean_reach_h8  = reach_at(horizons[min(2L, length(horizons))]),
      mean_reach_h13 = reach_at(Hmax),
      mean_estab_h13 = mean(alt$p_establishment, na.rm = TRUE),
      n_atrisk       = nrow(alt))
  }))
  byzone <- dplyr::bind_rows(lapply(grid, function(kk) {
    d <- atrisk(reach_by_k[[key(kk)]], Hmax); d$k_indiv <- kk; d }))
  if (!is.null(province))
    byzone <- dplyr::left_join(byzone,
      tibble::tibble(health_zone = names(province), province = unname(province)),
      by = "health_zone")
  list(summary = summary, byzone = byzone, k_base = k_base, n_mc = n_mc, horizons = horizons)
}

# ---------------------------------------------------------------------------
# Light simulation-based calibration of the invasion kernel (PLAN §6.5).
# Simulate invasion outcomes from the prior generative hazard on the real design
# rows, refit, and check the rank of each true parameter within its posterior is
# ~uniform. Few replicates by default (heavy); returns the rank statistics.
# ---------------------------------------------------------------------------
cascade_sbc <- function(design, cov_spec = CASCADE_COV_SPEC, n_sbc = 20L,
                        iter = 400L, seed = CASCADE_SEED, form = CASCADE_FORM) {
  # The check simulates from, and refits, the FIXED-BETA brms hazard. When the cascade deploys a
  # time-varying form its fitter is the Stan engine, which this does not exercise; it says so
  # rather than report a check of a different fitter as the deployed one's.
  if (!identical(form, "fixed"))
    return(list(rank_stats = NULL, n_used = 0L, skipped = TRUE,
                note = sprintf(paste0("SBC not run: the cascade hazard is the '%s' form (Stan engine); ",
                                      "this check covers the fixed-beta brms fitter only."), form)))
  # Seeded; the caller's RNG stream is restored (as simulate_cascade() does).
  .had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  .old <- if (.had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({ if (.had) assign(".Random.seed", .old, envir = globalenv())
            else suppressWarnings(rm(".Random.seed", envir = globalenv())) }, add = TRUE)
  set_fixed_seed(seed)
  d0 <- design$d; feat <- intersect(cov_spec, design$feat)
  # The truths are drawn from the SAME priors the fitter uses (00_config.R), so the ranks are
  # uniform under a correct fitter; typed values here would silently break that the day a prior
  # changes.
  .p0 <- .normal_prior_params(BAYES_PRIOR_INTERCEPT, "BAYES_PRIOR_INTERCEPT")
  .pb <- .normal_prior_params(BAYES_PRIOR_COEF, "BAYES_PRIOR_COEF")
  Z <- as.matrix(cbind(1, d0[, feat, drop = FALSE])); off <- d0$logLam
  ranks <- vector("list", n_sbc)
  for (s in seq_len(n_sbc)) {
    b0 <- stats::rnorm(1, .p0[1], .p0[2]); g <- stats::rnorm(length(feat), .pb[1], .pb[2])
    eta <- as.numeric(Z %*% c(b0, g)) + off
    y <- stats::rbinom(nrow(Z), 1, 1 - exp(-exp(eta)))       # cloglog
    if (sum(y) < 2L) next
    dd <- d0; dd$invaded <- y
    des <- design; des$d <- dd
    fit <- tryCatch(fit_bayes_renewal(des, cov_spec = feat, iter = iter, chains = 2L),
                    error = function(e) NULL)
    if (is.null(fit)) next
    dr <- posterior::as_draws_df(fit)
    truth <- setNames(c(b0, g), c("b_Intercept", paste0("b_", feat)))
    ranks[[s]] <- sapply(names(truth), function(nm)
      if (nm %in% names(dr)) mean(dr[[nm]] < truth[[nm]]) else NA_real_)
  }
  ranks <- do.call(rbind, ranks[!vapply(ranks, is.null, logical(1))])
  list(rank_stats = ranks, n_used = if (is.null(ranks)) 0L else nrow(ranks), skipped = FALSE,
       note = "SBC rank of truth within posterior should be ~Uniform(0,1) per parameter")
}
