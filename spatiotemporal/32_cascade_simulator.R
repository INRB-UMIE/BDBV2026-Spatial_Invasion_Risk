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
  G  <- daily_to_weekly_gt(layer$gt_pmfs[[gt]]); kmax <- length(G)

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
       G = G, kmax = kmax, Y0 = Y0, W0 = W0, affected0 = affected0,
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
                             k_indiv = CASCADE_K_INDIV,
                             n_seed_fun = cascade_draw_seed,
                             n_est = CASCADE_N_EST, n_rep = CASCADE_N_REP,
                             within_depletion = CASCADE_WITHIN_ZONE_DEPLETION,
                             pop_vec = NULL,
                             bar_sources = character(0), force_seed = character(0),
                             force_seed_week = 1L,
                             seed = CASCADE_SEED) {
  set.seed(seed)
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
  param_idx  <- sample.int(prep$n_draws, D, replace = TRUE)         # (beta0,gamma) draw per group
  R0_by_grp  <- vapply(seq_len(D), function(d) cascade_draw_R0(prep$reff), numeric(nz))  # nz x D
  draw_group <- rep(seq_len(D), each = n_rep)                        # length n_mc, iter -> draw id
  Nj <- if (!is.null(pop_vec)) pmax(as.numeric(pop_vec[zones_all]), 1) else rep(Inf, nz)

  # accumulators
  tau <- matrix(NA_integer_, nz, n_mc, dimnames = list(zones_all, NULL))  # first-passage (horizon week)
  # week at which cumulative incidence first reaches n_est (establishment first-passage);
  # tracked like `tau` so establishment is HORIZON-CONSISTENT (established-by-H <= reached-by-H
  # for every H, since a zone must be seeded before it can accumulate n_est cases).
  est_tau <- matrix(NA_integer_, nz, n_mc, dimnames = list(zones_all, NULL))
  new_by_week <- matrix(0L, horizon, n_mc)     # # new invasions each horizon week
  flux <- matrix(0, nz, nz, dimnames = list(zones_all, zones_all))  # flux[source, dest] expected seedings
  ncol_tot <- W0 + horizon

  for (m in seq_len(n_mc)) {
    g <- draw_group[m]; di <- param_idx[g]                # parameter draw for this replicate
    eta0 <- prep$b_int[di] + prep$g_logpop[di] * prep$z_logpop + prep$g_ccvi[di] * prep$z_ccvi
    gdmin <- prep$g_dmin[di]
    # incidence matrix (preallocated); first W0 cols = observed nowcast counts
    I <- matrix(0, nz, ncol_tot); I[, seq_len(W0)] <- prep$Y0
    infected <- prep$affected0
    # currently-infected zones' R_eff is the draw's fixed value (held across replicates)
    R0v <- numeric(nz); R0v[infected] <- R0_by_grp[infected, g]
    inf_num <- inf_num0
    dmin <- dmin0
    cumI <- rowSums(prep$Y0)                    # cumulative incidence (for establishment)
    est_week <- rep(NA_integer_, nz)            # per-iteration establishment week
    barred <- logical(nz); barred[bar_idx0] <- TRUE

    for (w in seq_len(horizon)) {
      cw <- W0 + w
      c_mult <- c_fun(w)
      # g-weighted own/infector pressure at cw (uses cols cw-1 .. cw-kmax)
      own <- numeric(nz)
      for (k in seq_len(kmax)) { tp <- cw - k; if (tp >= 1) own <- own + G[k] * I[, tp] }
      # susceptible fraction (optional within-zone depletion)
      s_frac <- if (within_depletion) pmin(pmax((Nj - cumI / max(CASCADE_RHO, 1e-3)) / Nj, 0), 1) else NULL
      # Layer A: advance infected zones' incidence for week cw
      Rt <- R0v * c_mult; if (!is.null(s_frac)) Rt <- Rt * s_frac
      inc <- cascade_branch_step(Rt, own, k_indiv)   # 0 where own==0 (at-risk zones)
      I[, cw] <- inc
      cumI <- cumI + inc

      # Layer B: seeding hazard for at-risk zones (import force = own routed through Wt)
      # own for the import force must EXCLUDE barred zones' contribution (knockout).
      own_src <- own
      if (any(barred)) own_src[barred] <- 0
      Lam <- as.numeric(Wt %*% own_src)
      atrisk <- !infected & Lam > 0
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
        mu  <- delta * exp(eta) * sat
        p   <- 1 - exp(-mu)
        u   <- stats::runif(nz)
        hit <- atrisk & (u < p)
        seeded_idx <- which(hit)
      }
      seeded_idx <- union(seeded_idx, forced_now)
      if (length(seeded_idx)) {
        # source attribution: multinomial over infected origins by W[j,i]*own_j
        for (i in seeded_idx) {
          wsrc <- Wt[i, ] * own_src               # W[j,i]*own_j over origins j
          tot <- sum(wsrc)
          if (tot > 0) { src <- sample.int(nz, 1L, prob = wsrc); flux[src, i] <- flux[src, i] + 1 }
        }
        ns <- n_seed_fun(length(seeded_idx))
        I[cbind(seeded_idx, cw)] <- pmax(ns, 1L)
        cumI[seeded_idx] <- cumI[seeded_idx] + pmax(ns, 1L)
        R0v[seeded_idx] <- cascade_draw_R_seed(prep$reff, seeded_idx)
        infected[seeded_idx] <- TRUE
        tau[cbind(seeded_idx, m)] <- w
        new_by_week[w, m] <- length(seeded_idx)
        # update frontier numerator + d_min for the new front (unless barred as a source):
        # both use add_idx (barred zones do not advance the front) for knockout consistency.
        add_idx <- seeded_idx[!barred[seeded_idx]]
        if (length(add_idx)) {
          inf_num <- inf_num + as.numeric(Wt[, add_idx, drop = FALSE] %*% rep(1, length(add_idx)))
          if (prep$has_dmin) for (z in add_idx) {
            col <- OS[, z]; dmin <- pmin(dmin, ifelse(is.finite(col), col, Inf))
          }
        }
      }
      # record establishment week (first week cumI crosses n_est), at-risk zones only
      je <- is.na(est_week) & cumI >= n_est & !prep$affected0
      if (any(je)) est_week[je] <- w
    }
    est_tau[, m] <- est_week
  }
  flux <- flux / n_mc
  list(tau = tau, est_tau = est_tau, new_by_week = new_by_week, flux = flux,
       n_mc = n_mc, horizon = horizon, zones_all = zones_all,
       affected0 = prep$affected0, prov = prep$prov,
       draw_group = draw_group, n_rep = n_rep, n_draws_used = D,
       knobs = list(delta = delta, psi = psi, k_indiv = k_indiv, scenario = scenario$label))
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
    } else s2_param <- rep(0, length(p))
    psd <- sqrt(s2_param)                                      # posterior SD of the reach probability
    lo <- pmax(p - z90 * psd, 0); hi <- pmin(p + z90 * psd, 1)
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
