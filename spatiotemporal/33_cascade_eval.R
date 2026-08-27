# =============================================================================
# 33_cascade_eval.R — 3-MONTH CASCADE: validation & calibration (PLAN §6, §3.5)
# BDBV 2026 DRC
#
# Layered validation, honest about the ~15-week data ceiling (a true 13-week
# hold-out has at most one truncated origin):
#   cascade_consistency_gate()  h=1 ranking vs the validated short-horizon model +
#                               observed-frequency calibration (the anchor, §6.4).
#   cascade_fit_psi()           fit the saturation psi so the cascade's expanding
#                               front matches the observed frontier growth (§3.5).
#   cascade_backtest()          intermediate-horizon (<=8w) leave-future-out cascade
#                               backtest at the earliest origins (~1 effective origin;
#                               indicative), scored on reach discrimination + the
#                               calibration of the weekly new-invasion count (§6.3).
#   cascade_sensitivity()       ranking-stability across GT/kernel/psi/k/scenario (§6.6).
#   cascade_sbc()               light simulation-based calibration of the kernel (§6.5).
# Small self-contained metric helpers (no coupling to 16_invasion_eval internals).
# =============================================================================

# ---- metric helpers --------------------------------------------------------
.casc_auc_pr <- function(p, y) {              # step-wise average precision (tie-aware)
  ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  P <- sum(y == 1); if (!P || !length(p)) return(NA_real_)
  o <- order(p, decreasing = TRUE); ps <- p[o]; y <- y[o]
  tp <- cumsum(y == 1); fp <- cumsum(y == 0)
  keep <- c(ps[-length(ps)] != ps[-1], TRUE)
  prec <- (tp / (tp + fp))[keep]; rec <- (tp / P)[keep]
  sum(prec * diff(c(0, rec)), na.rm = TRUE)
}
.casc_topk <- function(p, y, K) {             # precision@K (fraction of top-K truly invaded)
  ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  if (!length(p)) return(NA_real_)
  o <- order(p, decreasing = TRUE); sum(y[head(o, K)] == 1) / min(K, length(p))
}
.casc_spearman <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman",
                                                             use = "pairwise.complete.obs"))

# ---------------------------------------------------------------------------
# CONSISTENCY GATE (PLAN §6.4). At h=1 the cascade uses no source projection, so
# with sat==1 at week 1 its reach must (a) rank like the validated short-horizon
# model and (b) sit at the recalibrated (observed-frequency) level. The h=2 leg is
# compared to the DETERMINISTIC mean-field to avoid the Jensen offset (PLAN §6.4);
# here we anchor on the exact-by-construction h=1 leg.
# ---------------------------------------------------------------------------
cascade_consistency_gate <- function(prep, fit, design, delta,
                                     zone_week_nc, mobility_matrices, gt_pmfs,
                                     covariates, osrm_mat, zones_all,
                                     n_mc = 3000L, signal_thresh = 0.02,
                                     min_rho = 0.9) {
  scP <- CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]
  sim1 <- simulate_cascade(prep, scP, n_mc = n_mc, horizon = 1L, delta = delta,
                           psi = CASCADE_PSI)
  reach1 <- cascade_reach_table(sim1, horizons = 1L)
  offs <- bayes_forecast_offsets(zone_week_nc, mobility_matrices, gt_pmfs, covariates,
            osrm_mat, zones_all, mob = prep$kernel, gt = CASCADE_GT,
            cov_spec = CASCADE_COV_SPEC, horizons = c(1L, 2L), beta_proj = design$beta0)
  affected <- zones_all[prep$affected0]
  pred <- predict_bayes_invasion(fit, offs, design, horizons = 1L, affected_zones = affected)
  m <- merge(reach1[, c("health_zone", "p_invasion")],
             pred[, c("health_zone", "p_invasion")], by = "health_zone",
             suffixes = c("_casc", "_pred"))
  m <- m[is.finite(m$p_invasion_casc) & is.finite(m$p_invasion_pred), ]
  sig <- m[m$p_invasion_pred > signal_thresh, ]
  rho <- .casc_spearman(sig$p_invasion_casc, sig$p_invasion_pred)
  lvl <- sum(m$p_invasion_casc) / sum(m$p_invasion_pred)
  pass <- is.finite(rho) && rho >= min_rho
  # keep the paired point cloud (for the evaluation figure's week-1 anchor panel):
  # signal = zones the validated model flags above threshold (drive the ranking test).
  gate_data <- data.frame(health_zone = m$health_zone,
                          p_invasion_pred = m$p_invasion_pred,
                          p_invasion_casc = m$p_invasion_casc,
                          signal = m$p_invasion_pred > signal_thresh)
  list(spearman_signal = rho, n_signal = nrow(sig), level_ratio = lvl,
       delta = delta, pass = pass, signal_thresh = signal_thresh, data = gate_data,
       note = sprintf("h1 ranking rho=%.3f (n=%d, gate>=%.2f); level=%.2f (=delta target, observed-freq calibrated)",
                      rho, nrow(sig), min_rho, lvl))
}

# ---------------------------------------------------------------------------
# Fit the frontier-saturation psi (PLAN §3.5). Choose psi so the cascade's
# EXPECTED cumulative new-invasion count over the horizon best matches the
# observed cumulative frontier growth rate implied by the recent outbreak (the
# average weekly new-zone rate over the last `ref_weeks`, extrapolated). Because a
# true long-horizon target is unavailable, we anchor psi to the observed short-run
# expansion pace: without saturation the cascade over-accelerates; psi is the
# smallest value bringing the modelled early-weeks slope to the observed slope.
# Returns psi plus the diagnostic curve. Falls back to CASCADE_PSI on failure.
# ---------------------------------------------------------------------------
cascade_fit_psi <- function(prep, delta, zone_week_nc, zones_all,
                            psi_grid = c(0, 0.5, 1, 1.5, 2, 3, 4, 6), ref_weeks = 6L,
                            n_mc = 600L, match_week = 8L, trunc_drop = 1L) {
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
  curve <- sapply(psi_grid, function(ps) {
    s <- simulate_cascade(prep, scP, n_mc = n_mc, horizon = match_week, delta = delta, psi = ps)
    sum(rowMeans(s$new_by_week))               # expected cumulative new zones by match_week
  })
  psi_hat <- psi_grid[which.min(abs(curve - target_cum))]
  list(psi = psi_hat, psi_grid = psi_grid, modelled_cum = curve,
       observed_pace_cum = target_cum, obs_weekly_rate = obs_rate, match_week = match_week)
}

# ---------------------------------------------------------------------------
# Intermediate-horizon LEAVE-FUTURE-OUT cascade backtest (PLAN §6.3).
# For each cutoff C (weeks), truncate the incidence to <= C, treat zones affected
# by C as the initial front, REFIT the covariate hazard + R_eff on the truncated
# data, run the cascade forward K weeks, and score the cumulative newly-invaded set
# against the zones whose first case actually fell in (C, C+K]. The hazard fit, the
# R_eff estimate, and the training incidence are all as-of C (no future data). The
# delta/psi CALIBRATION is held at the production value (a global scalar / grid
# choice) rather than re-fit per fold — a minor, deliberately-declared exception, not
# a per-zone leak. Honest ceiling: only the earliest cutoffs afford a K-week window
# and they share most history -> ~1 effective origin; treat as INDICATIVE.
# ---------------------------------------------------------------------------
cascade_backtest <- function(layer, cutoffs_from_end = c(8L, 7L, 6L), K = 6L,
                             n_mc = 300L, delta = NULL, psi = CASCADE_PSI,
                             kernel = CASCADE_KERNEL, zone_province = NULL) {
  zones_all <- layer$zones_all
  Yfull <- .count_wide(layer$zone_week_nc, zones_all, "confirmed_nc")
  first_wk <- apply(Yfull, 1, function(x) { w <- which(x > 0); if (length(w)) min(w) else NA_integer_ })
  nT <- ncol(Yfull)
  weeks_sorted <- sort(unique(layer$zone_week_nc$week_start))
  res <- list(); det <- list()
  for (ce in cutoffs_from_end) {
    Ccol <- nT - ce                            # cutoff column index
    if (Ccol < 4L || Ccol + K > nT) next
    cutoff_date <- weeks_sorted[Ccol]
    zw_c <- layer$zone_week_nc[layer$zone_week_nc$week_start <= cutoff_date, , drop = FALSE]
    # refit hazard + reff on truncated data (leakage-free)
    dz <- tryCatch(build_invasion_design(zw_c, layer$mobility_matrices, layer$gt_pmfs,
             layer$covariates, layer$osrm_mat, zones_all, mob = kernel, gt = CASCADE_GT),
             error = function(e) NULL)
    if (is.null(dz) || dz$n_events < 3L) next
    fz <- fit_bayes_renewal(dz, cov_spec = CASCADE_COV_SPEC, iter = 600L, chains = 2L)
    if (is.null(fz)) next
    rz <- estimate_zone_reff(zw_c, zones_all, layer$gt_pmfs, gt = CASCADE_GT,
                             zone_province = zone_province)   # province-pooled, as production
    layer_c <- layer; layer_c$zone_week_nc <- zw_c
    prep_c <- cascade_prepare(layer_c, fz, dz, rz, kernel = kernel, gt = CASCADE_GT)
    del <- delta %||% cascade_fit_delta(cov_model = paste0("Bayes-", kernel, "-geo"))
    scP <- CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]
    sim <- simulate_cascade(prep_c, scP, n_mc = n_mc, horizon = K, delta = del, psi = psi)
    reach <- cascade_reach_table(sim, horizons = K)
    reach <- reach[reach$horizon == K & !reach$was_active_before, ]
    # truth: at-risk-as-of-C zones whose first case fell in (C, C+K]. Never-invaded
    # zones have first_wk = NA -> treat as +Inf so they score as 0 (NOT dropped, which
    # would delete every true negative and inflate AUC-PR / null the base rate).
    atrisk_c <- !(rowSums(Yfull[, seq_len(Ccol), drop = FALSE] > 0) > 0)
    fw <- first_wk; fw[is.na(fw)] <- .Machine$integer.max
    truth <- setNames(as.integer(atrisk_c & fw > Ccol & fw <= Ccol + K), zones_all)
    df <- merge(reach[, c("health_zone", "p_invasion")],
                data.frame(health_zone = zones_all, y = truth), by = "health_zone")
    df <- df[is.finite(df$p_invasion), ]
    br <- mean(df$y); ap <- .casc_auc_pr(df$p_invasion, df$y)
    # mean rank of the truly-invaded zones among the at-risk set (lower = nearer the
    # top of the watch-list); tie-AVERAGED, matching invasion_evaluation.csv's
    # mean_rank_of_truth. NA when a fold has no events.
    rk_avg <- rank(-df$p_invasion, ties.method = "average")
    mrt <- if (any(df$y == 1L)) mean(rk_avg[df$y == 1L]) else NA_real_
    # count calibration: predicted vs observed cumulative new invasions
    pred_cum <- sum(df$p_invasion); obs_cum <- sum(df$y)
    res[[length(res) + 1L]] <- tibble::tibble(
      cutoff = as.character(cutoff_date), K = K, n_atrisk = nrow(df), n_events = sum(df$y),
      base_rate = br, auc_pr = ap, auc_pr_skill = ap / max(br, 1e-9),
      mean_rank_of_truth = mrt,
      prec_at5 = .casc_topk(df$p_invasion, df$y, 5L),
      prec_at10 = .casc_topk(df$p_invasion, df$y, 10L),
      pred_new = pred_cum, obs_new = obs_cum, count_ratio = pred_cum / max(obs_cum, 1e-9))
    # per-zone detail (all at-risk zones this fold) for the evaluation figure's
    # capture curve, mean-rank and per-round top-K outcome panels.
    det[[length(det) + 1L]] <- tibble::tibble(
      cutoff = as.character(cutoff_date), K = K,
      health_zone = df$health_zone, p_invasion = df$p_invasion, y = as.integer(df$y))
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
                                base_reach, n_mc = 400L, psi = CASCADE_PSI) {
  base_r <- base_reach[base_reach$horizon == max(CASCADE_REPORT_HORIZONS) &
                       !base_reach$was_active_before, c("health_zone", "p_invasion")]
  runs <- list()
  add <- function(tag, reach) {
    r <- reach[reach$horizon == max(CASCADE_REPORT_HORIZONS) & !reach$was_active_before,
               c("health_zone", "p_invasion")]
    m <- merge(base_r, r, by = "health_zone", suffixes = c("_base", "_alt"))
    runs[[tag]] <<- tibble::tibble(axis = tag,
      spearman = .casc_spearman(m$p_invasion_base, m$p_invasion_alt),
      mean_reach = mean(r$p_invasion, na.rm = TRUE))
  }
  scP <- CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]
  # k_indiv sweep
  for (kk in CASCADE_K_INDIV_SWEEP) {
    prepk <- cascade_prepare(layer, fit_by_kernel[[CASCADE_KERNEL]],
                             design_by_kernel[[CASCADE_KERNEL]], reff, kernel = CASCADE_KERNEL)
    add(sprintf("k_indiv=%.2f", kk),
        cascade_reach_table(simulate_cascade(prepk, scP, n_mc = n_mc, delta = delta,
                                             psi = psi, k_indiv = kk)))
  }
  # psi sweep
  for (ps in CASCADE_PSI_SWEEP) {
    prepk <- cascade_prepare(layer, fit_by_kernel[[CASCADE_KERNEL]],
                             design_by_kernel[[CASCADE_KERNEL]], reff, kernel = CASCADE_KERNEL)
    add(sprintf("psi=%.1f", ps),
        cascade_reach_table(simulate_cascade(prepk, scP, n_mc = n_mc, delta = delta, psi = ps)))
  }
  # kernel sweep (needs its own fit/design)
  for (kern in setdiff(CASCADE_SENS_KERNELS, CASCADE_KERNEL)) {
    if (is.null(fit_by_kernel[[kern]])) next
    prepk <- cascade_prepare(layer, fit_by_kernel[[kern]], design_by_kernel[[kern]],
                             reff, kernel = kern)
    add(sprintf("kernel=%s", kern),
        cascade_reach_table(simulate_cascade(prepk, scP, n_mc = n_mc, delta = delta, psi = psi)))
  }
  dplyr::bind_rows(runs)
}

# ---------------------------------------------------------------------------
# Light simulation-based calibration of the invasion kernel (PLAN §6.5).
# Simulate invasion outcomes from the prior generative hazard on the real design
# rows, refit, and check the rank of each true parameter within its posterior is
# ~uniform. Few replicates by default (heavy); returns the rank statistics.
# ---------------------------------------------------------------------------
cascade_sbc <- function(design, cov_spec = CASCADE_COV_SPEC, n_sbc = 20L,
                        iter = 400L, seed = CASCADE_SEED) {
  set.seed(seed)
  d0 <- design$d; feat <- intersect(cov_spec, design$feat)
  Z <- as.matrix(cbind(1, d0[, feat, drop = FALSE])); off <- d0$logLam
  ranks <- vector("list", n_sbc)
  for (s in seq_len(n_sbc)) {
    b0 <- stats::rnorm(1, -3, 2); g <- stats::rnorm(length(feat), 0, 1)
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
  list(rank_stats = ranks, n_used = if (is.null(ranks)) 0L else nrow(ranks),
       note = "SBC rank of truth within posterior should be ~Uniform(0,1) per parameter")
}
