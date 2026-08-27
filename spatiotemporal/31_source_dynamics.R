# =============================================================================
# 31_source_dynamics.R — 3-MONTH CASCADE: Layer A (within-zone transmission)
# BDBV 2026 DRC · implements PLAN_3MONTH_INVASION.md §3.2
#
# Provides:
#   estimate_zone_reff()   per-zone current R_eff with province/national
#                          empirical-Bayes shrinkage (route ii of PLAN §3.2 —
#                          per-zone renewal ratio + shrinkage; the fast, robust
#                          fallback to a bespoke joint hierarchical Stan fit).
#   cascade_draw_R0()      one posterior/uncertainty draw of the R_eff vector.
#   cascade_reff_week()    R_j(t) = R_j0 * s_j(t) * c(t) for a projection week.
#   cascade_branch_step()  NegBin individual-offspring branching increment.
#
# Sourced AFTER 06_simple_models.R (compute_foi/daily_to_weekly_gt),
# 15_workhorse.R (.gweighted_own, estimate_R_local, .count_wide) and
# 30_projection_config.R. ASSUMPTIONS documented inline (Task 3).
# =============================================================================

# ---------------------------------------------------------------------------
# Per-zone current effective reproduction number, R_eff (NOT R0).
#
# For each zone we form the renewal ratio over the last `window` reporting-
# reasonably-complete weeks:  R_raw_j = sum_w I_j(w) / sum_w (g ⊛ I_j)(w).
# We then EMPIRICAL-BAYES SHRINK log R toward the zone's province mean and the
# province mean toward the national estimate, with weight w_j = n_j/(n_j+tau0)
# (data-rich zones move little; data-poor zones borrow strength). This is the
# route-(ii) approximation of PLAN §3.2; a joint hierarchical renewal (Stan) is
# the documented upgrade. Returns everything the simulator needs.
#
# @param zone_week_nc  nowcast-corrected zone-week tibble (needs confirmed_nc)
# @param zones_all     519-zone spine
# @param gt_pmfs,gt    generation-time PMFs + profile key
# @param zone_province named char vector zone -> province (optional; NULL = national only)
# @param window        # of trailing weeks used for the ratio (default 3)
# ---------------------------------------------------------------------------
estimate_zone_reff <- function(zone_week_nc, zones_all, gt_pmfs, gt = CASCADE_GT,
                               zone_province = NULL, window = 3L,
                               tau0 = CASCADE_R_MIN_CASES,
                               r_min = CASCADE_R_MIN, r_max = CASCADE_R_MAX) {
  Y <- .count_wide(zone_week_nc, zones_all, "confirmed_nc")
  G <- daily_to_weekly_gt(gt_pmfs[[gt]])
  nT <- ncol(Y); kmax <- length(G)
  R_nat <- as.numeric(estimate_R_local(Y, G, nT))       # national scalar, clamped [0.5,5]

  # Per-zone raw renewal ratio over the trailing `window` weeks that have a
  # defined g-weighted denominator (need at least one prior week of history).
  wks <- seq.int(max(2L, nT - window + 1L), nT)
  num <- den <- setNames(numeric(length(zones_all)), zones_all)
  for (w in wks) {
    lam_own <- .gweighted_own(Y, G, w)                  # sum_k G[k] Y[,w-k], vector over zones
    num  <- num  + Y[, w]
    den  <- den  + lam_own
  }
  # total recent cases (for the shrinkage weight): use full trailing window incidence
  ncase_all <- rowSums(Y[, wks, drop = FALSE])
  R_raw <- ifelse(den > 0 & num > 0, num / den, NA_real_)
  R_raw <- pmin(pmax(R_raw, r_min), r_max)
  logR_raw <- log(R_raw)

  # Province means of log R (case-weighted), each shrunk toward national.
  prov <- if (is.null(zone_province)) rep("ALL", length(zones_all)) else
            unname(zone_province[match(zones_all, names(zone_province))])
  prov[is.na(prov)] <- "ALL"
  logR_nat <- log(R_nat)
  prov_meanlog <- tapply(seq_along(zones_all), prov, function(idx) {
    lr <- logR_raw[idx]; wt <- ncase_all[idx]
    ok <- is.finite(lr) & wt > 0
    if (!any(ok)) return(logR_nat)
    m  <- sum(lr[ok] * wt[ok]) / sum(wt[ok])            # case-weighted province logR
    Wp <- sum(wt[ok]); wp <- Wp / (Wp + 5 * tau0)       # shrink province -> national
    wp * m + (1 - wp) * logR_nat
  })
  prov_pool <- as.numeric(prov_meanlog[prov])           # each zone's province pool mean (logR)

  # Zone-level shrinkage: data-rich -> own ratio, data-poor -> province pool.
  wj <- ncase_all / (ncase_all + tau0)
  logR_zone <- ifelse(is.finite(logR_raw), wj * logR_raw + (1 - wj) * prov_pool, prov_pool)
  R_zone <- pmin(pmax(exp(logR_zone), r_min), r_max)
  names(R_zone) <- zones_all

  # Estimation-uncertainty SD on log R per zone (wider when data-poor): used to
  # propagate R uncertainty across MC iterations for the initially-infected zones.
  sd_zone <- CASCADE_R_POOL_SD * (1 + tau0 / (ncase_all + tau0))
  names(sd_zone) <- zones_all
  # Seeding pool mean (log R) for a NEW zone = its province pool mean (already
  # shrunk to national); seeded zones inherit R_eff (~province level), not R0.
  seed_meanlog <- setNames(prov_pool, zones_all)

  list(R_zone = R_zone, logR_zone = logR_zone, sd_zone = sd_zone,
       seed_meanlog = seed_meanlog, seed_sdlog = CASCADE_R_POOL_SD,
       R_nat = R_nat, prov = prov, ncase = ncase_all, r_min = r_min, r_max = r_max)
}

# One MC draw of the initial R_eff vector for the currently-infected zones
# (propagates estimation uncertainty). LogNormal(logR_zone, sd_zone), clamped.
cascade_draw_R0 <- function(reff) {
  z <- stats::rnorm(length(reff$logR_zone))
  R <- exp(reff$logR_zone + reff$sd_zone * z)
  pmin(pmax(R, reff$r_min), reff$r_max)
}

# R for a freshly-seeded zone i (drawn at seeding time from the province R_eff pool).
cascade_draw_R_seed <- function(reff, zone_idx) {
  R <- exp(reff$seed_meanlog[zone_idx] + reff$seed_sdlog * stats::rnorm(length(zone_idx)))
  pmin(pmax(R, reff$r_min), reff$r_max)
}

# Effective R at projection week t (1-based into the horizon) for a scenario.
#   R_j(t) = R_j0 * s_j(t) * c(t)
# s_j(t) = susceptible fraction (1 unless within-zone depletion is on).
cascade_reff_week <- function(R0, t, c_mult, s_frac = NULL) {
  Rt <- R0 * c_mult
  if (!is.null(s_frac)) Rt <- Rt * s_frac
  Rt
}

# NegBin individual-offspring branching increment for one week (PLAN §3.2).
# mean = R_j * own_j ; size = own_j * k_indiv  (sum-of-branching identity), so
# extinction is possible when own_j is small and the draw concentrates when large.
# Zeros where own_j == 0 (no active infectors -> no incidence). Returns integer vec.
cascade_branch_step <- function(R_vec, own_vec, k_indiv = CASCADE_K_INDIV) {
  mu   <- R_vec * own_vec
  size <- own_vec * k_indiv
  out  <- integer(length(mu))
  act  <- own_vec > 0 & mu > 0
  if (any(act))
    out[act] <- stats::rnbinom(sum(act), mu = mu[act], size = pmax(size[act], 1e-6))
  out
}
