# =============================================================================
# 32_cascade_simulator.R — 3-MONTH CASCADE: coupled Layer A/B Monte-Carlo engine
# BDBV 2026 DRC · implements PLAN_3MONTH_INVASION.md §3.1–3.4, §3.5(sat)
#
# cascade_prepare()   precompute all static objects (transpose FOI, standardised
#                     static covariates, d_min osrm, frontier denominators, draws).
# simulate_cascade()  run M iterations of the coupled cascade for one scenario;
#                     returns first-passage times, reach, establishment, weekly
#                     new-invasion counts, and the seeding flux network.
# cascade_reach_table() summarise first-passage into per-zone reach at each horizon.
#
# The inner weekly recurrence (one MC iteration), t_for = W0 + w:
#   own_j   = sum_k G[k] I_j(t_for-k)             (g-weighted infector pressure)
#   I_j(w)  ~ NegBin(R_j(w)*own_j, own_j*k_indiv) (Layer A branching, infected j)
#   Lam_i   = (Wᵀ own)_i                          (import force = same own routed out)
#   mu_i    = delta * exp(eta0_i + gamma_dmin*z_dmin_i(w) + log Lam_i) * sat_i(w)
#   seed_i  ~ Bernoulli(1 - exp(-mu_i))           (Layer B, at-risk i)
# Verified against predict_bayes_invasion: the (eta0 + gamma*z + logLam) hazard is
# byte-identical (recon smoke test, cor=1, max|diff|=0).
#
# Sourced AFTER 30/31 and the pipeline core (06/15/20/21). ASSUMPTIONS inline.
# =============================================================================

# ---------------------------------------------------------------------------
# Precompute the static simulation context (does NOT depend on scenario / draw).
# @param layer  list(zone_week_nc, zones_all, covariates, gt_pmfs, osrm_mat,
#                     mobility_matrices)
# @param fit    covariate-bearing brms invasion fit (draws source)
# @param design build_invasion_design() output (center/scale/feat/beta0)
# @param reff   estimate_zone_reff() output
# ---------------------------------------------------------------------------
cascade_prepare <- function(layer, fit, design, reff,
                            kernel = CASCADE_KERNEL, gt = CASCADE_GT,
                            cov_spec = CASCADE_COV_SPEC) {
  zones_all <- layer$zones_all; nz <- length(zones_all)
  W  <- layer$mobility_matrices[[kernel]]
  if (is.null(W)) stop("cascade_prepare: mobility kernel '", kernel, "' not found")
  W  <- W[zones_all, zones_all]                       # align to spine order
  Wt <- t(W)                                          # inflow: Lam = Wt %*% own
  denom_fi <- colSums(W)                              # sum_j W[j,i] (frontier-sat denominator)
  G  <- cascade_weekly_gt(gt, layer$gt_pmfs); kmax <- length(G)

  # ---- per-GT weekly kernels, when the projection marginalises over GT_PRIOR -----------
  # Built ONLY when reff carries an R posterior per grid point, so the kernel and the
  # reproduction number always come from the same grid point. Each kernel is built with
  # weekly_censored_gt(mean, sd) DIRECTLY rather than through cascade_weekly_gt(), because that
  # helper looks up GT_PROFILES first and a grid key has no profile entry — it would silently
  # fall through to BINNING the daily pmf, which is the older mis-censored construction and a
  # different kernel from the one every other part of this suite uses.
  G_list <- NULL; gt_weights <- NULL; gt_keys <- NULL
  if (!is.null(reff$R_draws_by_gt) && length(reff$R_draws_by_gt) > 1L) {
    gt_keys    <- reff$gt_keys
    gt_weights <- as.numeric(reff$gt_weights)
    gg <- reff$gt_grid
    G_list <- lapply(seq_along(gt_keys), function(i)
      weekly_censored_gt(gg$gt_mean[i], gg$gt_sd[i]))
    names(G_list) <- gt_keys
    message(sprintf(paste0("[cascade_prepare] GT marginalisation ON: %d weekly kernels, mean ",
                           "lag %.2f-%.2f weeks (single-GT kernel %.2f)"),
                    length(G_list),
                    min(vapply(G_list, function(x) sum(x * seq_along(x)), numeric(1))),
                    max(vapply(G_list, function(x) sum(x * seq_along(x)), numeric(1))),
                    sum(G * seq_along(G))))
  }

  Y0 <- .count_wide(layer$zone_week_nc, zones_all, "confirmed_nc")
  Y0 <- Y0[zones_all, , drop = FALSE]
  W0 <- ncol(Y0)
  affected0 <- rowSums(Y0 > 0) > 0                    # initial infected set (absorbing origins)

  # Static covariate raw vectors, standardised with the FIT design center/scale.
  static <- .static_features(layer$covariates, zones_all)
  z_of <- function(nm, raw) {
    c0 <- design$center[[nm]]; s0 <- design$scale[[nm]]
    if (is.null(c0)) return(rep(0, nz))
    v <- (raw - c0) / s0; v[!is.finite(v)] <- 0; v
  }
  z_logpop <- if ("log_pop" %in% cov_spec) z_of("log_pop", static$log_pop) else rep(0, nz)
  z_ccvi   <- if ("ccvi"    %in% cov_spec) z_of("ccvi",    static$ccvi)    else rep(0, nz)

  # OSRM travel-time matrix aligned to the spine (missing -> Inf); for dynamic d_min.
  osrm <- layer$osrm_mat
  OS <- matrix(Inf, nz, nz, dimnames = list(zones_all, zones_all))
  ii <- intersect(zones_all, rownames(osrm)); jj <- intersect(zones_all, colnames(osrm))
  OS[ii, jj] <- osrm[ii, jj]
  OS[!is.finite(OS) | OS <= 0] <- Inf                 # self / unroutable -> Inf
  dmin_center <- design$center[["d_min"]]; dmin_scale <- design$scale[["d_min"]]
  has_dmin <- ("d_min" %in% cov_spec) && !is.null(dmin_center)
  DMIN_SENTINEL <- 1e4                                # matches .dmin_vec

  # Posterior draws of (beta0, gamma) — the covariate-bearing hazard.
  dr <- posterior::as_draws_df(fit)
  b_int <- dr[["b_Intercept"]]
  gvec <- function(nm) if (paste0("b_", nm) %in% names(dr)) dr[[paste0("b_", nm)]] else rep(0, length(b_int))
  g_logpop <- gvec("log_pop"); g_ccvi <- gvec("ccvi"); g_dmin <- gvec("d_min")
  n_draws <- length(b_int)

  list(zones_all = zones_all, nz = nz, W = W, Wt = Wt, denom_fi = denom_fi,
       G = G, kmax = kmax, G_list = G_list, gt_weights = gt_weights, gt_keys = gt_keys,
       Y0 = Y0, W0 = W0, affected0 = affected0,
       z_logpop = z_logpop, z_ccvi = z_ccvi, OS = OS, has_dmin = has_dmin,
       dmin_center = dmin_center, dmin_scale = dmin_scale, DMIN_SENTINEL = DMIN_SENTINEL,
       b_int = b_int, g_logpop = g_logpop, g_ccvi = g_ccvi, g_dmin = g_dmin,
       n_draws = n_draws, reff = reff, prov = reff$prov, kernel = kernel, gt = gt)
}

# Initial d_min travel-time vector from an active set (min over active columns of OS).
.cascade_dmin0 <- function(OS, active_idx, sentinel) {
  nz <- nrow(OS)
  if (!length(active_idx)) return(rep(sentinel, nz))
  sub <- OS[, active_idx, drop = FALSE]
  d <- apply(sub, 1, function(x) { x <- x[is.finite(x)]; if (length(x)) min(x) else sentinel })
  d[!is.finite(d)] <- sentinel; d
}

# ---------------------------------------------------------------------------
# Run the cascade for one scenario. Returns first-passage & summaries.
# @param prep         cascade_prepare() output
# @param scenario     entry of CASCADE_SCENARIOS
# @param n_mc         iterations
# @param delta,psi,k_indiv,n_seed_fun  calibration/branching knobs
# @param bar_sources  character zone names barred from EVER seeding others (knockout)
# @param force_seed_week  horizon week at which force_seed zones are seeded (timing axis)
# @param force_seed   character zone names force-seeded at `force_seed_week` (designed
#                     conditional experiment)
# ---------------------------------------------------------------------------
simulate_cascade <- function(prep, scenario, n_mc = CASCADE_N_MC,
                             horizon = CASCADE_HORIZON_WEEKS,
                             delta = CASCADE_DELTA_DEFAULT, psi = CASCADE_PSI,
                             delta_sd_log = 0,
                             delta_band = get0("RECAL_BAND", ifnotfound = c(0.01, 5.0)),
                             k_indiv = CASCADE_K_INDIV,
                             n_seed_fun = cascade_draw_seed,
                             n_est = CASCADE_N_EST, n_rep = CASCADE_N_REP,
                             pop_vec = NULL,   # ACCEPTED BUT UNUSED since depletion was removed; see the note below
                             bar_sources = character(0), force_seed = character(0),
                             force_seed_week = 1L,
                             seed_r = get0("CASCADE_SEED_R", ifnotfound = NA_real_),
                             seed_r_scope = c("forced", "all"),
                             rw_sigma = get0("CASCADE_R_RW_SIGMA", ifnotfound = 0),
                             # Fallback matches the config default: a PURE random walk. It was
                             # 0.30, so a caller reaching this fallback silently got the
                             # retired mean-reverting walk instead of the configured one.
                             rw_rho   = get0("CASCADE_R_RW_RHO",   ifnotfound = 0),
                             crn = get0("CASCADE_CRN", ifnotfound = TRUE),
                             seed = CASCADE_SEED) {
  # RNG hygiene: leave the caller's stream exactly as it was found. Without this, simulating
  # a cascade silently changes every subsequent random draw in the session.
  .old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(.old_seed)) assign(".Random.seed", .old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  seed_r_scope <- match.arg(seed_r_scope)
  if (!is.null(seed_r) && length(seed_r) == 1L && is.na(seed_r)) seed_r <- NULL
  nz <- prep$nz; zones_all <- prep$zones_all; W0 <- prep$W0; G <- prep$G; kmax <- prep$kmax
  Wt <- prep$Wt; denom_fi <- prep$denom_fi; OS <- prep$OS
  c_fun <- scenario$c_fun
  bar_idx0 <- which(zones_all %in% bar_sources)
  force_idx <- which(zones_all %in% force_seed)
  aff0_idx <- which(prep$affected0)
  # For a clean knockout, a barred zone contributes to NEITHER the import force NOR
  # the frontier (d_min / saturation): exclude it from the initial front too.
  aff0_src_idx <- setdiff(aff0_idx, bar_idx0)
  # frontier numerator for the initial infected set: sum over infected origins j of W[j,i]
  inf_num0 <- if (length(aff0_src_idx)) as.numeric(Wt[, aff0_src_idx, drop = FALSE] %*% rep(1, length(aff0_src_idx))) else numeric(nz)
  dmin0 <- if (prep$has_dmin) .cascade_dmin0(OS, aff0_src_idx, prep$DMIN_SENTINEL) else NULL
  # Frontier saturation is measured RELATIVE to the initial front (PLAN §3.5: the
  # hazard falls as neighbours *become* invaded). f_inv0 is the baseline invaded
  # fraction of each zone's source neighbourhood, so sat = 1 at week 1 (h1 ranking
  # is preserved, consistency gate) and only bites as NEW zones are seeded.
  f_inv0 <- ifelse(denom_fi > 0, inf_num0 / denom_fi, 0)
  # Nested design for a proper CREDIBLE interval on the reach probability: draw D posterior
  # PARAMETER sets (hazard coefficients (beta0, gamma) AND the currently-infected zones'
  # effective reproduction numbers), and run n_rep stochastic process replicates within each.
  # Parameters are held FIXED across the replicates of a draw, so between-draw variation is
  # posterior (parameter) uncertainty and within-draw variation is process noise — separated
  # in cascade_reach_table to form the credible interval (not the M-shrinking sampling band).
  n_rep <- max(2L, as.integer(n_rep))
  D <- max(1L, n_mc %/% n_rep); n_mc <- D * n_rep
  # COMMON RANDOM NUMBERS. Every contrast in this layer — urban seeding, gateway knockout,
  # conditional hub seeding, the psi search — is a DIFFERENCE between two simulations, and
  # the draws inside the week loop have lengths that depend on how many zones are active or
  # newly seeded. Two runs therefore desynchronise the moment an intervention changes
  # anything, and the difference measures the random stream rather than the intervention. In
  # the 2026-09-07 urban run that produced a NEGATIVE expected number of added invasions for
  # most cities, which is structurally impossible, and credited the largest increases from
  # seeding Kinshasa to zones beside the actual outbreak 1,500 km away.
  #
  # The fix is to make the draws ADDRESSABLE rather than sequential: one stream per
  # iteration, fixed-shape uniform matrices generated up front, and every draw an inverse
  # CDF indexed by ZONE. The iteration seeds are drawn BEFORE anything else so they depend
  # only on (seed, n_mc) and not on how much randomness the pre-loop draws happen to use.
  #
  # HOW FAR THE COUPLING GOES — stated with its conditions, because two earlier versions of
  # this comment overstated it and a regression test caught both.
  #
  # THE COUPLING ITSELF IS EXACT and always holds: an unreachable zone's difference is
  # exactly 0 (329 of 349 far zones, against 30 unpaired), because every draw is an inverse
  # CDF keyed by zone and both runs read the same key. That is a property of the RNG scheme.
  #
  # MONOTONICITY is a property of the MODEL, and needs two conditions:
  #   (i)  psi = 0, and
  #   (ii) gamma_dmin <= 0 for the parameter draw in use.
  # Under both, adding a source can only raise the import force, the seeding indicator can
  # only flip 0 -> 1, and scaling `own` scales the negative binomial's mean and size
  # together (stochastically ordered), so no per-zone difference can be negative.
  #
  # TWO SEPARATE CHANNELS BREAK IT. Both are model features, not coupling failures:
  #   * FRONTIER SATURATION (psi > 0). sat = exp(-psi * max(f_inv - f_inv0, 0)) is a
  #     negative feedback: seeding a city raises its neighbours' invaded fraction and lowers
  #     a third zone's hazard. This was DOMINANT while psi was fitted (it reached 11.25 on
  #     one frame). psi is now fixed at 0 (30_projection_config.R), so sat == 1 and this
  #     channel is INACTIVE in the default configuration — it revives only if psi is swept
  #     above 0. With it gone, a negative per-zone difference is noise or the re-seeding
  #     channel below, and must not be explained as saturation.
  #   * THE FRONTIER COVARIATE'S POSTERIOR TAIL. eta carries gamma_dmin * z_dmin, and
  #     seeding LOWERS d_min for neighbours — which raises their hazard only when gamma_dmin
  #     is negative. gamma_dmin is a posterior DRAW, not the mean, and its posterior straddles
  #     zero (mean -0.11, P(>0) = 0.27 at iter=600; tighter but still non-zero at iter=2000),
  #     so a minority of parameter groups invert the sign. MEASURED at psi = 0 on one fit:
  #     as-is gave 1-2 negatives per 466 zones at horizons 2, 3, 6 and 13; constraining
  #     gamma_dmin <= 0 gave EXACTLY ZERO at every one of those horizons, as did dropping the
  #     d_min term. Note the honest limit of that evidence: on a SECOND fit the same
  #     constraint still left one negative, so this channel is demonstrated and dominant but
  #     is not proven to be the only one. What was ruled out is timing desynchronisation —
  #     the worry that a zone seeded in different weeks reads different branching uniforms
  #     and so follows a different path. If that were a live channel the negatives would grow
  #     with horizon; constrained, they stay at zero from horizon 1 to 13.
  # Do not "fix" either by clipping; both are consequences of the fitted model. NOTHING in
  # the reported analysis depends on exact pathwise monotonicity — the variance reduction
  # does not require it, and it is measured directly (per-zone SD 0.017 -> 0.0006).
  #
  # PAIRING PRECONDITION. Two runs pair only if they share `seed`, `n_mc`, `n_rep` and
  # `horizon`: the nested design splits iterations into parameter groups, and a different
  # iteration count changes both the grouping and the per-iteration seeds. The contrast
  # helpers CHECK this and warn — they do not stop, because an unpaired contrast is still an
  # estimate, just a far noisier one than its caller is likely to assume.
  iter_seed <- if (isTRUE(crn)) sample.int(.Machine$integer.max, n_mc) else NULL
  param_idx  <- sample.int(prep$n_draws, D, replace = TRUE)         # (beta0,gamma) draw per group
  # GENERATION TIME, drawn per parameter group with its PRIOR WEIGHT. Drawn HERE, inside the
  # seeded region and alongside param_idx, so two runs sharing (seed, n_mc, n_rep, horizon) draw
  # the same grid points and a scenario contrast stays paired — drawn anywhere else it would
  # silently unpair every contrast in this layer.
  .has_gtgrid <- !is.null(prep$G_list) && length(prep$G_list) > 1L
  gt_idx <- if (.has_gtgrid)
    sample.int(length(prep$gt_weights), D, replace = TRUE, prob = prep$gt_weights)
  else rep(1L, D)
  # R is drawn from the posterior fitted AT this group's generation time, so the kernel used
  # below and the reproduction number used here are the same grid point.
  R0_by_grp  <- vapply(seq_len(D), function(d)
    cascade_draw_R0(prep$reff,
                    gt_key = if (.has_gtgrid) prep$gt_keys[gt_idx[d]] else NULL),
    numeric(nz))                                                    # nz x D
  # CALIBRATION UNCERTAINTY. delta is a FITTED nuisance with a zone-clustered interval
  # (33b_cascade_calibration.R), and the reach LEVEL is more sensitive to it than to any
  # other single quantity. Holding it fixed while (beta0, gamma) and R_eff are drawn made
  # the credible intervals understate exactly the uncertainty that matters most, so it is
  # drawn ONCE PER PARAMETER GROUP on the log scale — the scale its bootstrap SE is
  # reported on — and clamped to the same band the estimator searched. delta_sd_log = 0
  # (the default) reproduces the previous fixed-delta behaviour exactly.
  delta_by_grp <- if (is.finite(delta_sd_log) && delta_sd_log > 0)
    pmin(pmax(as.numeric(delta) * exp(delta_sd_log * stats::rnorm(D)), delta_band[1]), delta_band[2])
  else rep(as.numeric(delta), D)
  draw_group <- rep(seq_len(D), each = n_rep)                        # length n_mc, iter -> draw id
  # R-WALK INNOVATIONS, DRAWN PER PARAMETER GROUP (not per iteration).
  #
  # This walk exists to make the credible interval widen with horizon — 30_projection_config.R
  # introduces it because "the week-13 credible interval [was] exactly as tight as the week-4
  # one, which cannot be right". Drawn per ITERATION it could not do that, because iterations
  # inside one parameter group are the PROCESS replicates and cascade_reach_table() treats all
  # within-group variation as process noise and subtracts it (s2_param = var_total - within).
  #
  # Measured directly on that decomposition (nz=300, D=150, n_rep=8, mean interval width),
  # varying only WHERE the same shock is drawn:
  #     sigma      per-iteration   per-group
  #     0.00         0.0615          0.0615     (identical, as they must be)
  #     0.15         0.0666          0.1259
  #     0.30         0.0662          0.2972
  #     0.60         0.0746          0.6387
  # A per-group shock passes through in full; a per-iteration shock is LARGELY absorbed and
  # moves the interval only ~20% even at sigma 0.6 — and not monotonically. (It is not removed
  # to exactly zero: `within` estimates the process term as p(1-p)/(n_rep-1), i.e. assumes the
  # within-group variation is pure Bernoulli, so extra within-group variance from a
  # per-iteration shock partly leaks into s2_param. On the full simulator the net effect at the
  # shipped settings was that the interval did NOT widen, and at large sigma narrowed.)
  #
  # The future trajectory of R is EPISTEMIC uncertainty about the scenario, not observation
  # noise within a trajectory, so it belongs at the parameter-group level — the same level as
  # delta_by_grp, which does widen the interval as intended. Drawn LAST in the seeded region so
  # iter_seed, param_idx, gt_idx, R0_by_grp and delta_by_grp keep their exact streams, and
  # keyed on (zone, week, group) so a baseline and its intervention still see identical shocks
  # and every paired contrast stays paired.
  Z_Rw_grp <- if (isTRUE(crn) && is.finite(rw_sigma) && rw_sigma > 0)
    array(stats::rnorm(nz * horizon * D), dim = c(nz, horizon, D)) else NULL
  # WITHIN-ZONE SUSCEPTIBLE DEPLETION WAS REMOVED (2026-09-19). It was controlled by
  # CASCADE_WITHIN_ZONE_DEPLETION, which was FALSE and was never overridden by any caller in
  # this suite, so `s_frac` was NULL and the whole path was dead. Its only implementation also
  # divided cumulative CONFIRMED cases by a fixed, guessed ascertainment constant
  # (CASCADE_RHO = 0.45) to convert them to infections — a number this pipeline does not
  # estimate and no longer carries anywhere. The modelling assumption it encoded is unchanged
  # and is the documented one: over a 13-week horizon EVD does not materially deplete a health
  # zone's susceptibles, and the R decline that matters is carried by the pre-registered
  # scenario multiplier c(t). `pop_vec` is still ACCEPTED because 15 call sites thread it, but
  # it is no longer read; it is listed in the dead-argument inventory for removal.

  # accumulators
  tau <- matrix(NA_integer_, nz, n_mc, dimnames = list(zones_all, NULL))  # first-passage (horizon week)
  # week at which cumulative incidence first reaches n_est (establishment first-passage);
  # tracked like `tau` so establishment is HORIZON-CONSISTENT (established-by-H <= reached-by-H
  # for every H, since a zone must be seeded before it can accumulate n_est cases).
  est_tau <- matrix(NA_integer_, nz, n_mc, dimnames = list(zones_all, NULL))
  new_by_week <- matrix(0L, horizon, n_mc)     # # new invasions each horizon week
  # PROJECTED CASES. The incidence matrix I was already being built to drive the renewal
  # term, the establishment threshold and depletion — and then discarded, so the layer could
  # say how many ZONES an intervention adds but not how many CASES. Both are kept now.
  #   cases_zone    nz x n_mc, new confirmed cases per zone over the projection window ONLY
  #                 (columns W0+1 .. W0+horizon; the observed history in Y0 is excluded, so
  #                 a contrast between two runs is not diluted by shared past incidence).
  #   cases_by_week horizon x n_mc, the national weekly total, for trajectories.
  # Both are DOUBLE, not integer: at a swept seeded R these sums can run above the integer
  # range that `new_by_week` safely uses for zone counts, and rowSums() returns double
  # regardless. COST, stated because it is not free: cases_zone is nz x n_mc doubles, so at
  # M = 4000 it adds ~17 MB per simulation and roughly doubles the object (tau and est_tau
  # are nz x n_mc integers). The urban grid holds a baseline plus one conditional run at a
  # time, so the peak is ~70 MB — modest, but worth knowing before raising M.
  cases_zone    <- matrix(0, nz, n_mc, dimnames = list(zones_all, NULL))
  cases_by_week <- matrix(0, horizon, n_mc)
  flux <- matrix(0, nz, nz, dimnames = list(zones_all, zones_all))  # flux[source, dest] expected seedings
  ncol_tot <- W0 + horizon

  # Per-iteration random-number buffers, allocated once and refilled inside the loop.
  # Allocated SEPARATELY, not chained: `a <- b <- matrix(...)` binds one object to both
  # names, and the in-place refills below would then rely on copy-on-write to separate them.
  .mk <- function() if (isTRUE(crn)) matrix(0, nz, horizon) else NULL
  U_seed <- .mk(); U_br <- .mk(); U_flux <- .mk()
  U_ns <- if (isTRUE(crn)) numeric(nz) else NULL
  Z_Rs <- if (isTRUE(crn)) numeric(nz) else NULL
  # Innovations for the R walk, keyed by (zone, week) like every other CRN stream so a
  # baseline and its intervention see the SAME shocks and the contrast stays paired.
  Z_Rw <- .mk()

  # First TRUE of a cumulative-weight comparison: the paired analogue of sample.int(prob=).
  .pick_src <- function(wts, u) {
    cw <- cumsum(wts); tot <- cw[length(cw)]
    if (!is.finite(tot) || tot <= 0) return(NA_integer_)
    which.max(cw >= u * tot)
  }

  for (m in seq_len(n_mc)) {
    # Fixed-shape draws for the whole iteration. Seeding, branching and flux are keyed by
    # (zone, week); the initial seed size and a seeded zone's reproduction number are keyed
    # by ZONE ALONE, so a zone seeded in different weeks in the two runs still starts from
    # the same size and the same R — which is what keeps the coupling monotone.
    if (isTRUE(crn)) {
      set.seed(iter_seed[m])
      # Refilled IN PLACE rather than reallocated: three nz x horizon matrices per iteration
      # is ~160 KB, which at M = 4000 is 650 MB of allocation churn and the GC pressure that
      # goes with it, for no benefit. The number and ORDER of runif/rnorm calls is unchanged,
      # so the streams — and therefore every result — are identical either way.
      U_seed[] <- stats::runif(nz * horizon)
      U_br[]   <- stats::runif(nz * horizon)
      U_flux[] <- stats::runif(nz * horizon)
      U_ns[]   <- stats::runif(nz)
      Z_Rs[]   <- stats::rnorm(nz)
    }
    g <- draw_group[m]; di <- param_idx[g]                # parameter draw for this replicate
    # This GROUP's R-walk innovations: shared by every process replicate in the group, so the
    # walk is posterior (between-group) variation and survives the variance decomposition.
    # Refilled IN PLACE (like the other CRN streams): reuses the one nz x horizon allocation
    # instead of materialising a fresh slice every iteration. Z_Rw was allocated by .mk() above.
    if (!is.null(Z_Rw_grp)) Z_Rw[] <- Z_Rw_grp[, , g]
    delta_m <- delta_by_grp[g]                            # calibration factor for this group
    # This group's weekly generation-time kernel. Falls back to the single kernel when the
    # projection is not marginalising over GT. kmax varies with the grid point (a longer GT
    # supports more lags); the lag loop already guards tp >= 1, so a longer kernel is safe.
    Gm <- if (.has_gtgrid) prep$G_list[[gt_idx[g]]] else G
    kmax_m <- length(Gm)
    # Weekly lag-0 mass of THIS draw's kernel. weekly_censored_gt() drops the lag-0 term and
    # renormalises (5.8% at the medium profile) because an explicit weekly recursion cannot
    # carry a self-referential term — but the R this simulator deploys is EpiNow2's DAILY-scale
    # posterior (31_source_dynamics.R sets R_draws from the EpiNow2 fit), which retains it.
    # Applying a daily-scale R to a lag-0-dropped kernel lengthens the effective generation
    # interval and so understates growth for a given R. .R_eff_m() is the exact algebraic
    # correction; it is the identity at R = 1 and when the kernel carries no lag-0 mass, and it
    # is per-draw because the GT grid point (hence G_0) varies with the draw.
    .g0_m <- { a <- attr(Gm, "lag0_mass"); if (is.null(a) || !is.finite(a)) 0 else as.numeric(a) }
    .R_eff_m <- if (.g0_m > 0) function(R) weekly_renewal_R_eff(R, .g0_m) else function(R) R
    eta0 <- prep$b_int[di] + prep$g_logpop[di] * prep$z_logpop + prep$g_ccvi[di] * prep$z_ccvi
    gdmin <- prep$g_dmin[di]
    # incidence matrix (preallocated); first W0 cols = observed nowcast counts
    I <- matrix(0, nz, ncol_tot); I[, seq_len(W0)] <- prep$Y0
    infected <- prep$affected0
    # currently-infected zones' R_eff is the draw's fixed value (held across replicates)
    R0v <- numeric(nz); R0v[infected] <- R0_by_grp[infected, g]
    # The walk reverts to each zone's OWN starting R, so it adds spread without adding trend.
    logR_anchor <- numeric(nz)
    logR_anchor[infected] <- log(pmax(R0v[infected], 1e-12))
    inf_num <- inf_num0
    dmin <- dmin0
    cumI <- rowSums(prep$Y0)                    # cumulative incidence (for establishment)
    cz <- numeric(nz)                           # projection-window cases, this iteration
    est_week <- rep(NA_integer_, nz)            # per-iteration establishment week
    barred <- logical(nz); barred[bar_idx0] <- TRUE

    for (w in seq_len(horizon)) {
      cw <- W0 + w
      c_mult <- c_fun(w)
      # g-weighted own/infector pressure at cw (uses cols cw-1 .. cw-kmax)
      own <- numeric(nz)
      for (k in seq_len(kmax_m)) { tp <- cw - k; if (tp >= 1) own <- own + Gm[k] * I[, tp] }
      # susceptible fraction (optional within-zone depletion)
      # Layer A: advance infected zones' incidence for week cw
      # R_eff EVOLVES. Mean-reverting walk on log R toward each zone's anchor: adds
      # uncertainty that grows with horizon and then plateaus, without adding a trend (the
      # scenario multiplier c(t) owns policy). rw_sigma = 0 reproduces constant-R exactly.
      if (rw_sigma > 0 && any(infected)) {
        .ii  <- which(infected)
        .lr  <- log(pmax(R0v[.ii], 1e-12))
        .eps <- if (isTRUE(crn)) Z_Rw[.ii, w] else stats::rnorm(length(.ii))
        .lr  <- .lr + rw_rho * (logR_anchor[.ii] - .lr) + rw_sigma * .eps
        R0v[.ii] <- pmin(pmax(exp(.lr), prep$reff$r_min), prep$reff$r_max)
      }
      Rt <- R0v * c_mult
      # Convert to the weekly-kernel multiplier AFTER the scenario multiplier and the
      # susceptible-depletion factor, both of which act on the daily-scale R.
      inc <- cascade_branch_step(.R_eff_m(Rt), own, k_indiv,
                                 u = if (isTRUE(crn)) U_br[, w] else NULL)  # 0 where own==0
      I[, cw] <- inc
      cumI <- cumI + inc

      # Layer B: seeding hazard for at-risk zones (import force = own routed through Wt)
      # own for the import force must EXCLUDE barred zones' contribution (knockout).
      own_src <- own
      if (any(barred)) own_src[barred] <- 0
      Lam <- as.numeric(Wt %*% own_src)
      # AT RISK OF SEEDING is "has no transmission of its own right now", NOT "has never been
      # invaded". A zone whose outbreak burned out (own == 0: no incidence inside the
      # generation-time window) can be re-seeded from elsewhere, which is what actually happens
      # in this outbreak — zones go quiet and reappear. The old test, !infected, made the first
      # invasion absorbing and allowed each zone exactly one introduction for all 13 weeks.
      #
      # `infected` stays the EVER-invaded flag and is what first-passage accounting uses below:
      # tau, new_by_week and the flux attribution must keep counting first invasions only, or a
      # re-seeding would overwrite a zone's invasion week with a later one and be reported as a
      # new invasion. Those are separated explicitly where they are written.
      atrisk <- (own == 0) & Lam > 0
      # forced conditional seeding at week 1 (designed experiment)
      forced_now <- if (w == force_seed_week && length(force_idx))
                      force_idx[!infected[force_idx]] else integer(0)
      seeded_idx <- integer(0)
      if (any(atrisk)) {
        z_dmin <- if (prep$has_dmin) (log1p(dmin) - prep$dmin_center) / prep$dmin_scale else rep(0, nz)
        if (prep$has_dmin) z_dmin[!is.finite(z_dmin)] <- 0
        f_inv <- ifelse(denom_fi > 0, inf_num / denom_fi, 0)
        sat <- exp(-psi * pmax(f_inv - f_inv0, 0))   # increment over the initial front
        eta <- eta0 + gdmin * z_dmin + log(pmax(Lam, 1e-12))
        mu  <- delta_m * exp(eta) * sat
        p   <- 1 - exp(-mu)
        u   <- if (isTRUE(crn)) U_seed[, w] else stats::runif(nz)
        hit <- atrisk & (u < p)
        seeded_idx <- which(hit)
      }
      seeded_idx <- union(seeded_idx, forced_now)
      if (length(seeded_idx)) {
        # FIRST invasions only, separated here once and used for every piece of first-passage
        # accounting below. Zones can now be re-seeded after burning out (see `atrisk` above),
        # and a re-seeding is a new INTRODUCTION but not a new INVASION: it must not overwrite
        # the zone's tau, must not be counted in new_by_week, and must not be attributed a
        # source in `flux`, all three of which describe when and whence a zone was FIRST
        # reached. Everything that is genuinely per-introduction (seed size, the R reset, the
        # walk anchor) still applies to the whole seeded set.
        newly <- seeded_idx[!infected[seeded_idx]]
        # source attribution: multinomial over infected origins by W[j,i]*own_j
        for (i in newly) {
          wsrc <- Wt[i, ] * own_src               # W[j,i]*own_j over origins j
          tot <- sum(wsrc)
          if (tot > 0) {
            src <- if (isTRUE(crn)) .pick_src(wsrc, U_flux[i, w]) else sample.int(nz, 1L, prob = wsrc)
            if (!is.na(src)) flux[src, i] <- flux[src, i] + 1
          }
        }
        # n_seed_fun now converts UNIFORMS to counts (30_projection_config.R), so the seed
        # size is paired by zone rather than drawn from a stream whose position has shifted.
        ns <- n_seed_fun(if (isTRUE(crn)) U_ns[seeded_idx] else stats::runif(length(seeded_idx)))
        I[cbind(seeded_idx, cw)] <- pmax(ns, 1L)
        cumI[seeded_idx] <- cumI[seeded_idx] + pmax(ns, 1L)
        # The seeded-R override describes the SCENARIO'S city, not every zone the front
        # happens to reach afterwards. Applying it to all of them would quietly re-parameterise
        # the whole projection, so by default it touches only the zones the scenario forced;
        # organically seeded zones take the replicate's national R. seed_r_scope = "all" is
        # available for a deliberate whole-projection sensitivity.
        #
        # R_grp is THIS parameter group's national R — the same number every already-infected
        # zone in this replicate started from. Passing it is what keeps a seeded zone on the
        # shared national draw instead of drawing independently, which would reintroduce the
        # spurious cross-zone independence cascade_draw_R0() exists to remove.
        .zr <- if (isTRUE(crn)) Z_Rs[seeded_idx] else NULL
        .Rg <- R0_by_grp[1L, g]
        # This group's GT grid point, passed alongside R_grp so the seed_r normaliser and any
        # fallback draw come from the SAME posterior the group's kernel belongs to.
        .Gk <- if (.has_gtgrid) prep$gt_keys[gt_idx[g]] else NULL
        if (is.null(seed_r) || identical(seed_r_scope, "all")) {
          R0v[seeded_idx] <- cascade_draw_R_seed(prep$reff, seeded_idx, z = .zr,
                                                 seed_r = seed_r, R_grp = .Rg, gt_key = .Gk)
        } else {
          .is_forced <- seeded_idx %in% force_idx
          R0v[seeded_idx] <- cascade_draw_R_seed(prep$reff, seeded_idx, z = .zr,
                                                 seed_r = NULL, R_grp = .Rg, gt_key = .Gk)
          if (any(.is_forced))
            R0v[seeded_idx[.is_forced]] <- cascade_draw_R_seed(
              prep$reff, seeded_idx[.is_forced],
              z = if (isTRUE(crn)) Z_Rs[seeded_idx[.is_forced]] else NULL,
              seed_r = seed_r, R_grp = .Rg, gt_key = .Gk)
        }
        # A zone seeded mid-projection anchors on the R it was seeded with, not on a value
        # it never had, so its walk starts from its own draw.
        logR_anchor[seeded_idx] <- log(pmax(R0v[seeded_idx], 1e-12))
        infected[seeded_idx] <- TRUE
        if (length(newly)) tau[cbind(newly, m)] <- w
        new_by_week[w, m] <- length(newly)
        # update frontier numerator + d_min for the new front (unless barred as a source):
        # both use add_idx (barred zones do not advance the front) for knockout consistency.
        # `newly`, NOT `seeded_idx`. seeded_idx includes RE-SEEDINGS of already-invaded zones
    # (atrisk deliberately re-admits burnt-out zones), so every re-introduction added that
    # zone's column to inf_num again and f_inv — documented as "the mobility-weighted invaded
    # FRACTION of i's source neighbourhood", a number in [0,1] — exceeded 1 (measured max 1.35).
    # First invasions are already separated for tau/new_by_week/flux one line above; the
    # frontier update was the one place that used the wrong set. Inert while CASCADE_PSI = 0
    # (sat == 1), live in the psi sweep that feeds cascade_sensitivity.csv.
    add_idx <- newly[!barred[newly]]
        if (length(add_idx)) {
          inf_num <- inf_num + as.numeric(Wt[, add_idx, drop = FALSE] %*% rep(1, length(add_idx)))
          if (prep$has_dmin) for (z in add_idx) {
            col <- OS[, z]; dmin <- pmin(dmin, ifelse(is.finite(col), col, Inf))
          }
        }
      }
      # Cases for this week, read HERE and not immediately after the branching step:
      # `I[, cw] <- inc` is overwritten a few lines above by `I[cbind(seeded_idx, cw)] <- ns`
      # for every zone seeded this week, so reading it earlier would miss each zone's seed
      # cases. cumI is incremented separately for those zones and the branching step returns
      # 0 for them (own == 0 before seeding), so there is no double count either way.
      #
      # Accumulated into `cz` rather than re-derived after the loop as
      # rowSums(I[, W0 + seq_len(horizon)]): that subset COPIES a nz x horizon matrix on
      # every iteration, which measured ~3x slower on the M = 2000 scenario stage. The column
      # is already extracted here for the weekly total, so the running sum is free.
      .inc_w <- I[, cw]
      cases_by_week[w, m] <- sum(.inc_w)
      cz <- cz + .inc_w

      # record establishment week (first week cumI crosses n_est), at-risk zones only
      je <- is.na(est_week) & cumI >= n_est & !prep$affected0
      if (any(je)) est_week[je] <- w
    }
    est_tau[, m] <- est_week
    # Projection-window cases only (weeks 1..horizon); the observed history in Y0 is never
    # added, so a contrast between two runs is not diluted by shared past incidence. `cz`
    # accumulated exactly the columns the weekly series summed, so the two cannot drift —
    # asserted by the colSums invariant in the regression tests. A zero-length horizon
    # leaves cz at its initial zeros, which is the right answer.
    cases_zone[, m] <- cz
  }
  flux <- flux / n_mc
  list(tau = tau, est_tau = est_tau, new_by_week = new_by_week, flux = flux,
       cases_zone = cases_zone, cases_by_week = cases_by_week,
       n_mc = n_mc, horizon = horizon, zones_all = zones_all,
       affected0 = prep$affected0, prov = prep$prov,
       draw_group = draw_group, n_rep = n_rep, n_draws_used = D,
       knobs = list(delta = as.numeric(delta), delta_sd_log = delta_sd_log,
                    delta_drawn_mean = mean(delta_by_grp),
                    psi = psi, k_indiv = k_indiv, crn = isTRUE(crn),
                    seed_r = if (is.null(seed_r)) NA_real_ else as.numeric(seed_r),
                    seed_r_scope = seed_r_scope,
                    rw_sigma = rw_sigma, rw_rho = rw_rho,
                    seed = seed, n_rep = n_rep, scenario = scenario$label))
}

# ---------------------------------------------------------------------------
# Summarise first-passage into a per-zone reach table at each reporting horizon.
# Returns a tibble mirroring the short-horizon risk-score columns so it can be
# fed through compute_risk_scores() / add_risk_indices() / .write_risk_csv().
#   p_invasion = reach R_i(H); mu_forecast = cumulative hazard proxy -log(1-R);
#   p_lo/p_hi  = 90% CREDIBLE interval on the reach probability (posterior spread
#               across parameter draws; process noise removed by variance decomposition);
#   p_sd = posterior SD; est = establishment prob.
# ---------------------------------------------------------------------------
#' Row-wise quantiles, using matrixStats when available and a base fallback otherwise
#' (mirrors matrixStats_rowVars, which this file already defines for the same reason).
matrixStats_rowQuantiles <- function(x, probs) {
  if (requireNamespace("matrixStats", quietly = TRUE))
    return(matrixStats::rowQuantiles(as.matrix(x), probs = probs, na.rm = TRUE))
  t(apply(as.matrix(x), 1, stats::quantile, probs = probs, na.rm = TRUE, names = FALSE))
}

cascade_reach_table <- function(sim, horizons = CASCADE_REPORT_HORIZONS) {
  zones <- sim$zones_all; M <- sim$n_mc; R <- sim$n_rep
  # M x D indicator (1/R in each parameter-draw group) to average replicates within a draw
  dg <- sim$draw_group; D <- max(dg)
  Gnorm <- Matrix_or_matrix_group(dg, D, R, M)
  z90 <- 1.6448536
  out <- lapply(horizons, function(H) {
    reached_iter <- (!is.na(sim$tau) & sim$tau <= H) * 1.0     # nz x M (0/1)
    p <- rowMeans(reached_iter)                                # marginal reach = overall fraction
    # per-parameter-draw reach probability (mean over that draw's R process replicates)
    p_draw <- reached_iter %*% Gnorm                           # nz x D
    if (D >= 2L && R >= 2L) {
      var_total <- matrixStats_rowVars(p_draw)                 # between-draw variance of estimates
      within    <- rowMeans(p_draw * (1 - p_draw)) / (R - 1L)  # process variance of each draw's mean
      s2_param  <- pmax(var_total - within, 0)                 # posterior (parameter) variance of reach
    } else {
      # var_total must exist on EVERY path: the interval block below reads it, but it was only
      # assigned inside `D >= 2 && R >= 2` while that block is guarded on `D >= 2` alone. With
      # D >= 2 and n_rep < 2 the function aborted with "object 'var_total' not found".
      # simulate_cascade() forces n_rep >= 2, but cascade_reach_table() is called directly on
      # arbitrary sim objects (33b, 38, 41), so the guard has to hold on its own terms.
      # Zero parameter variance also makes .shrink 0, which is the correct degenerate answer.
      s2_param  <- rep(0, length(p))
      var_total <- rep(0, length(p))
    }
    psd <- sqrt(s2_param)                                      # posterior SD of the reach probability
    # EMPIRICAL QUANTILES of the per-draw reach, not a Gaussian +-1.645*SD on a bounded
    # probability. The moment-matched interval collapsed to ZERO WIDTH for 288 of 458 at-risk
    # zones (63%) at h=13 — a "90% credible interval" asserting the probability is known exactly
    # — because pmax(var_total - within, 0) clips a noisy negative variance estimate to a point
    # mass. It also produced intervals clipped at 0/1 for another 102 zones. Those rows are drawn
    # as invisible bars in the retained Figure4_cascade forest (8 of 20 rows at h=4).
    # p_draw already holds the per-parameter-draw reach, so its 5th/95th percentiles are the
    # interval the docstring describes, are bounded in [0,1] by construction, and need no
    # normal approximation.
    if (D >= 2L) {
      # POSTERIOR-ONLY, as every consumer's label claims ("posterior spread across parameter
      # draws, process noise removed": 36_report.R's SS2 note, 35_cascade_viz.R's uncertainty-map
      # subtitle, 37_cascade_figure4.R's forest bars and Figure4_cascade_h*_data.csv).
      #
      # The raw quantiles of p_draw do NOT satisfy that: p_draw is the mean of only R process
      # replicates, so Var(p_draw) = s2_param + p(1-p)/R, and at the shipped R = 8 the process
      # term is the MAJORITY of the spread for a mid-range zone. Taking them unshrunk made the
      # published interval ~29% wider than the quantity it is labelled as.
      #
      # Shrink the empirical quantiles toward p by sqrt(s2_param / Var(p_draw)). That keeps the
      # SHAPE of the posterior (asymmetry near 0/1, which the old +-1.645*sd Gaussian destroyed)
      # while scaling it to the posterior-only width, and it cannot collapse to a point the way
      # pmax(var_total - within, 0) could, because of the resolution floor below.
      qs <- matrixStats_rowQuantiles(p_draw, probs = c(0.05, 0.95))
      .shrink <- ifelse(is.finite(var_total) & var_total > 0,
                        sqrt(pmin(s2_param / pmax(var_total, .Machine$double.eps), 1)), 0)
      lo <- pmax(p + (qs[, 1] - p) * .shrink, 0)
      hi <- pmin(p + (qs[, 2] - p) * .shrink, 1)
      # Floor the width at the Monte-Carlo resolution: with D draws the quantiles cannot
      # resolve below 1/D, and a genuinely degenerate row should read as "unresolved", not
      # "certain".
      .res <- 1 / D
      .flat <- (hi - lo) < .res & p > 0 & p < 1
      if (any(.flat)) { lo[.flat] <- pmax(p[.flat] - .res / 2, 0)
                        hi[.flat] <- pmin(p[.flat] + .res / 2, 1) }
    } else {
      lo <- pmax(p - z90 * psd, 0); hi <- pmin(p + z90 * psd, 1)
    }
    est <- rowMeans(!is.na(sim$est_tau) & sim$est_tau <= H)    # established BY H (horizon-consistent)
    fp_med <- apply(sim$tau, 1, function(x) { x <- x[!is.na(x) & x <= H]; if (length(x)) stats::median(x) else NA_real_ })
    masked <- sim$affected0
    tibble::tibble(
      health_zone = zones, horizon = as.integer(H),
      mu_forecast = ifelse(masked, NA_real_, -log(pmax(1 - p, 1e-12))),
      p_invasion  = ifelse(masked, NA_real_, p),
      p_case_invasion = ifelse(masked, NA_real_, p),  # confirmed-case scale (= reach)
      p_lo = ifelse(masked, NA_real_, lo), p_hi = ifelse(masked, NA_real_, hi),
      p_sd = ifelse(masked, NA_real_, psd),
      p_establishment = ifelse(masked, NA_real_, est),
      first_passage_week = ifelse(masked, NA_real_, fp_med),
      was_active_before = masked)
  })
  dplyr::bind_rows(out)
}

# ---------------------------------------------------------------------------
# PAIRED contrast between two simulations sharing a random-number stream
# ---------------------------------------------------------------------------

#' Per-zone and total effect of an intervention, differenced ITERATION BY ITERATION.
#'
#' Every scenario question in this layer is a difference. Differencing two separately
#' summarised runs throws away the pairing and leaves an estimate whose noise floor swamped
#' the signal: at M = 4000 the 2026-09-07 urban run reported a NEGATIVE expected number of
#' added invasions for most cities and ranked zones 1,500 km from the seeded city as the
#' most affected. (A small negative is not by itself impossible — at psi > 0 frontier
#' saturation is a genuine negative feedback — but negatives of that MAGNITUDE, spread over
#' zones the seeding cannot reach, are noise and nothing else.) Under common random numbers the same contrast returns
#' exact zeros for 94% of far zones (against 9% unpaired) and a per-zone standard error
#' 25-30x smaller (0.017 -> 0.0006 over the zones the seeding cannot reach).
#'
#' THE UNIT OF INDEPENDENCE IS THE PARAMETER DRAW, NOT THE ITERATION. The nested design runs
#' n_rep process replicates inside each posterior draw, so iterations within a group share
#' (beta0, gamma, delta, R_eff) and are positively correlated. Dividing by sqrt(M) would
#' understate the standard error by roughly sqrt(n_rep); the differences are therefore
#' averaged within each draw group first, and the spread is taken across the D groups.
#'
#' MATERIALITY IS A TEST, NOT A THRESHOLD. A zone counts as elevated when its paired 90%
#' interval excludes zero — replacing a fixed 0.03 cut that was neither a noise floor nor an
#' effect size, and that under pairing would discard real small effects while admitting
#' nothing it used to exclude.
#'
#' @param exclude zones to drop (the force-seeded hubs: invaded by assumption, not
#'   prediction). Zones already infected at t0 are dropped automatically.
#' @return list(per_zone, total, se, lo, hi, n_groups, paired)
#' Are two simulations genuinely paired? Shared by every contrast helper so they cannot
#' drift apart. `horizon` is part of the precondition, not just of the question: the
#' per-iteration uniform matrices are nz x horizon, so two runs with different horizons draw
#' from differently shaped streams and are not coupled at all.
.cascade_check_paired <- function(sim_base, sim_alt, label = "contrast") {
  kb <- sim_base$knobs; ka <- sim_alt$knobs
  .same <- function(a, b) isTRUE(all.equal(a, b))
  ok <- isTRUE(kb$crn) && isTRUE(ka$crn) &&
        identical(kb$seed, ka$seed) && identical(sim_base$n_mc, sim_alt$n_mc) &&
        identical(kb$n_rep, ka$n_rep) &&
        identical(sim_base$draw_group, sim_alt$draw_group) &&
        identical(sim_base$horizon, sim_alt$horizon) &&
        # THE STOCHASTIC-STRUCTURE KNOBS MATTER TOO, and were not checked. Two runs identical
        # on everything above but differing in delta_sd_log carry systematically different
        # per-group calibration factors, so the "paired" difference partly measures the delta
        # randomisation rather than the intervention; and because delta_by_grp only consumes
        # rnorm(D) when delta_sd_log > 0, the two runs also diverge in the seeded region from
        # that point on. rw_sigma/rw_rho are the same argument for the R walk (drawn per
        # parameter group). Nothing published differs on these today — every contrast uses the
        # defaults on both sides — but the guard has to enforce what it claims.
        .same(kb$delta_sd_log, ka$delta_sd_log) &&
        .same(kb$rw_sigma, ka$rw_sigma) && .same(kb$rw_rho, ka$rw_rho)
  if (!ok)
    warning(sprintf(paste0("[%s] these two runs are NOT paired (crn/seed/n_mc/n_rep/horizon/",
                           "draw_group/delta_sd_log/rw_sigma/rw_rho differ). The difference is ",
                           "an unpaired estimate; its standard error is far larger than the ",
                           "paired one and small effects are not interpretable."), label),
            call. = FALSE)
  ok
}

#' Mean and interval of a per-iteration difference, aggregated over PARAMETER DRAWS.
#'
#' The nested design runs n_rep correlated process replicates inside each posterior draw, so
#' the iteration is not the unit of independence: dividing by sqrt(M) would understate the
#' standard error by roughly sqrt(n_rep). Differences are averaged within each draw group
#' first and the spread taken across the D groups.
#'
#' @param d_iter numeric vector of length n_mc, or a matrix with n_mc COLUMNS.
#' @return for a vector, list(estimate, se, lo, hi); for a matrix, a list of equal-length
#'   vectors, one entry per row.
.cascade_group_summary <- function(d_iter, draw_group, conf = 0.90) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  D <- max(draw_group)
  if (is.matrix(d_iter)) {
    Gn <- matrix(0, ncol(d_iter), D)
    Gn[cbind(seq_len(ncol(d_iter)), draw_group)] <- 1
    Gn <- sweep(Gn, 2, pmax(colSums(Gn), 1), "/")
    grp <- d_iter %*% Gn                                   # rows x D group means
    est <- rowMeans(grp)
    se  <- if (D >= 2L) apply(grp, 1, stats::sd) / sqrt(D) else rep(NA_real_, nrow(grp))
    return(list(estimate = est, se = se, lo = est - z * se, hi = est + z * se, groups = grp))
  }
  grp <- vapply(seq_len(D), function(k) mean(d_iter[draw_group == k]), numeric(1))
  est <- mean(grp)
  se  <- if (D >= 2L) stats::sd(grp) / sqrt(D) else NA_real_
  list(estimate = est, se = se, lo = est - z * se, hi = est + z * se, groups = grp)
}

cascade_paired_contrast <- function(sim_base, sim_cond, horizon, exclude = character(),
                                    conf = 0.90) {
  stopifnot(identical(sim_base$zones_all, sim_cond$zones_all))
  paired <- .cascade_check_paired(sim_base, sim_cond, "contrast")

  H  <- as.integer(horizon)
  if (H > sim_base$horizon || H > sim_cond$horizon)
    stop(sprintf("[contrast] horizon %d exceeds a simulation's horizon (%d / %d).",
                 H, sim_base$horizon, sim_cond$horizon), call. = FALSE)
  rb <- (!is.na(sim_base$tau) & sim_base$tau <= H) * 1.0      # nz x M
  rc <- (!is.na(sim_cond$tau) & sim_cond$tau <= H) * 1.0
  D_iter <- rc - rb

  keep <- !sim_base$affected0 & !(sim_base$zones_all %in% exclude)
  dg <- sim_base$draw_group; D <- max(dg)
  # Per zone, aggregated over draw groups by the shared helper so this and the case contrast
  # cannot drift apart in how they treat the nested design.
  pz <- .cascade_group_summary(D_iter, dg, conf)
  # The TOTAL is summed WITHIN each draw group before its spread is taken, so the interval
  # carries the correlation between zones inside a group. Combining the per-zone standard
  # errors instead would treat the zones as independent and understate it badly.
  tot <- .cascade_group_summary(colSums(pz$groups[keep, , drop = FALSE]), seq_len(D), conf)

  per_zone <- tibble::tibble(
    health_zone = sim_base$zones_all, horizon = H,
    p_base = rowMeans(rb), p_cond = rowMeans(rc),
    delta = pz$estimate, delta_se = pz$se,
    delta_lo = pz$lo, delta_hi = pz$hi,
    # A zone that never differed in ANY iteration is exactly zero, not merely
    # indistinguishable from zero; keeping the two apart is what makes the null
    # intervention check meaningful.
    exact_zero = rowSums(abs(D_iter)) == 0,
    eligible = keep)
  per_zone$elevated <- per_zone$eligible & is.finite(per_zone$delta_lo) & per_zone$delta_lo > 0
  per_zone$reduced  <- per_zone$eligible & is.finite(per_zone$delta_hi) & per_zone$delta_hi < 0

  list(per_zone = per_zone[order(-per_zone$delta), ],
       total = tot$estimate, se = tot$se, lo = tot$lo, hi = tot$hi,
       n_groups = D, n_mc = sim_base$n_mc, paired = paired, horizon = H, conf = conf)
}

#' Total new invasions per ITERATION, over the whole horizon.
#'
#' The per-iteration vector rather than its mean, so a contrast between two runs can be
#' differenced iteration by iteration instead of as a difference of two noisy means.
cascade_new_iter <- function(sim) colSums(sim$new_by_week)

#' Paired difference in expected NEW invasions between two runs (base minus intervention).
#'
#' Positive = the intervention averts invasions. As in cascade_paired_contrast(), the unit
#' of independence is the parameter draw, not the iteration, because the nested design runs
#' n_rep correlated process replicates inside each draw.
cascade_paired_total <- function(sim_base, sim_alt, conf = 0.90) {
  stopifnot(identical(sim_base$n_mc, sim_alt$n_mc),
            identical(sim_base$draw_group, sim_alt$draw_group),
            # new_by_week has one row per horizon week; summing over different numbers of
            # rows would compare invasions counted over different lengths of time.
            identical(sim_base$horizon, sim_alt$horizon))
  d <- cascade_new_iter(sim_base) - cascade_new_iter(sim_alt)
  r <- .cascade_group_summary(d, sim_base$draw_group, conf)
  list(estimate = r$estimate, se = r$se, lo = r$lo, hi = r$hi,
       paired = .cascade_check_paired(sim_base, sim_alt, "knockout"))
}

#' Paired difference in projected CASES between two runs, split in-city vs elsewhere.
#'
#' WHY THIS IS A SEPARATE FUNCTION FROM cascade_paired_contrast(). The invasion contrast asks
#' how many additional ZONES are invaded, so it excludes two sets by construction: the
#' force-seeded hub (invaded by assumption) and zones already infected at t0 (they cannot be
#' newly invaded). NEITHER exclusion is right for cases:
#'   * the seeded city is where most of the case burden lands, and dropping it would discard
#'     the largest part of the answer;
#'   * an already-infected zone CAN receive more cases when another city is seeded, but NOT
#'     through the import force directly. In this simulator `Lam` enters only the Bernoulli
#'     seeding of zones with no current transmission (`atrisk <- (own == 0) & Lam > 0`); a
#'     zone that is already transmitting draws NegBin(R_eff*own, own*k) with no import term,
#'     so its incidence is invariant to `Lam` that week. The only channel is burn-out (own
#'     falls to 0 across the whole generation-time window) followed by re-seeding.
#'     LIMITATION, stated rather than implied: `cases_elsewhere` therefore omits importation
#'     into zones that are transmitting continuously, and is an UNDER-statement of the
#'     attributable downstream burden. (This note used to assert the opposite mechanism.)
#' So every zone contributes here, and the split is by LOCATION (seeded city vs the rest)
#' rather than by eligibility. Reporting only the sum would hide that the two behave very
#' differently across a seeded-R sweep.
#'
#' WHAT "CASES" MEANS. Confirmed cases on the projection window only (weeks 1..horizon), on
#' the same confirmed-case scale as the rest of the suite — so ascertainment (rho ~ 0.45) is
#' already inside the calibration and these are modelled CONFIRMED cases, not infections.
#' The in-city total also contains the SEED CASES THEMSELVES, which are an assumption rather
#' than a prediction: E[n_seed] cases per forced zone (about 1.5 under the default Poisson
#' seeding), so roughly 4-5 cases for a three-zone metro hub, not 1-2. They are deliberately
#' left in — netting them out would report a downstream-only figure under an in-city label —
#' but the distinction matters when reading the pooled arm, where most introductions fade and
#' the assumed seed is a large share of the total, against a swept R of 2.5 where the onward
#' epidemic dominates. `n_in_city` is returned so the assumed component can be recovered.
#'
#' @param in_city zones treated as the seeded city (normally the force_seed hub zones).
#' @return list(per_zone, in_city, elsewhere, total, paired, n_groups, n_mc, conf) where the
#'   three totals are each list(estimate, se, lo, hi).
cascade_paired_cases <- function(sim_base, sim_cond, in_city = character(), conf = 0.90) {
  stopifnot(identical(sim_base$zones_all, sim_cond$zones_all))
  if (is.null(sim_base$cases_zone) || is.null(sim_cond$cases_zone))
    stop("[cases] these simulations carry no `cases_zone`; they predate case tracking and ",
         "cannot be contrasted on cases. Re-run them.", call. = FALSE)
  paired <- .cascade_check_paired(sim_base, sim_cond, "cases")

  D_iter <- sim_cond$cases_zone - sim_base$cases_zone         # nz x M, cond MINUS base
  dg <- sim_base$draw_group; D <- max(dg)
  zones <- sim_base$zones_all
  is_city <- zones %in% in_city
  if (length(in_city) && !any(is_city))
    warning("[cases] none of the supplied `in_city` zones is on the spine; the in-city total ",
            "will be zero and the split is meaningless.", call. = FALSE)

  pz <- .cascade_group_summary(D_iter, dg, conf)
  # Each total is summed WITHIN a draw group first, so its interval keeps the correlation
  # between zones in that group; combining per-zone SEs would treat zones as independent.
  #
  # TWO KINDS OF INTERVAL, answering different questions. `lo`/`hi` are Monte-Carlo
  # uncertainty on the EXPECTATION — how precisely E[additional cases] is known. `median`
  # with `pred_lo`/`pred_hi` describes the SPREAD OF OUTCOMES across iterations — what might
  # actually happen. They diverge sharply here because projected case counts are
  # right-skewed: parameter uncertainty enters an exponential, so E[R^t] far exceeds
  # median(R)^t. On this frame the baseline projection has mean 15.7k against median 13.2k,
  # with the top 5% of iterations carrying 15% of the mean. Reporting only the mean would
  # overstate the typical outcome; only the median would understate the expected burden. The
  # suite's scenario burden table already uses medians with predictive bands, so this follows
  # that convention rather than inventing a second one.
  qp <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  .tot <- function(sel) {
    if (!any(sel))
      return(list(estimate = 0, se = 0, lo = 0, hi = 0, median = 0, pred_lo = 0, pred_hi = 0))
    r <- .cascade_group_summary(colSums(pz$groups[sel, , drop = FALSE]), seq_len(D), conf)
    per_iter <- colSums(D_iter[sel, , drop = FALSE])          # length n_mc, not n_groups
    q <- unname(stats::quantile(per_iter, qp))
    c(r, list(median = unname(stats::median(per_iter)), pred_lo = q[1], pred_hi = q[2]))
  }
  t_city <- .tot(is_city); t_else <- .tot(!is_city); t_all <- .tot(rep(TRUE, length(zones)))

  per_zone <- tibble::tibble(
    health_zone = zones,
    cases_base = rowMeans(sim_base$cases_zone), cases_cond = rowMeans(sim_cond$cases_zone),
    delta = pz$estimate, delta_se = pz$se, delta_lo = pz$lo, delta_hi = pz$hi,
    exact_zero = rowSums(abs(D_iter)) == 0,
    in_city = is_city, was_active_before = sim_base$affected0)
  per_zone$elevated <- is.finite(per_zone$delta_lo) & per_zone$delta_lo > 0
  per_zone$reduced  <- is.finite(per_zone$delta_hi) & per_zone$delta_hi < 0

  # NOTE for anyone reading the three totals together: `estimate` IS additive
  # (in_city + elsewhere == total, asserted in the tests) because expectations add. `median`
  # and the predictive bounds are NOT — quantiles of a sum are not the sum of quantiles — so
  # do not present the medians as a decomposition.
  list(per_zone = per_zone[order(-per_zone$delta), ],
       in_city = t_city, elsewhere = t_else, total = t_all,
       n_in_city = sum(is_city), paired = paired, n_groups = D,
       n_mc = sim_base$n_mc, conf = conf)
}

# Build an M x D group-averaging matrix (column d has 1/R in the rows of draw group d).
Matrix_or_matrix_group <- function(dg, D, R, M) {
  G <- matrix(0, M, D)
  G[cbind(seq_len(M), dg)] <- 1 / R
  G
}
# Row-wise variance (population-style / (D-1)) without a matrixStats dependency.
matrixStats_rowVars <- function(x) {
  n <- ncol(x); mu <- rowMeans(x)
  rowSums((x - mu)^2) / (n - 1)
}
