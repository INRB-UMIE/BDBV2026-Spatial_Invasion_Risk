# =============================================================================
# 30_projection_config.R — 3-MONTH SPATIAL INVASION CASCADE: configuration
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Implements the configuration layer of PLAN_3MONTH_INVASION.md. This module
# defines EVERY tunable of the 13-week stochastic metapopulation cascade and
# documents the modelling ASSUMPTION behind each one (Task 3: assumptions must be
# explicit in the code). Nothing here is estimated; estimation lives in 31/33.
#
# Sourced AFTER 00_config.R (uses OUT_DIR, ANALYSIS_DATE, RANDOM_SEED,
# ASCERTAINMENT_NOMINAL, GT_PRIMARY). Depends on nothing else.
# =============================================================================

# ---------------------------------------------------------------------------
# Horizon grid (PLAN §1)
# ---------------------------------------------------------------------------
# 13 weekly steps ≈ 3 months. Weekly grid anchored exactly as the short-horizon
# suite (the current week ends on ANALYSIS_DATE). Reporting horizons ≈ 1/2/3 mo.
CASCADE_HORIZON_WEEKS  <- 13L
CASCADE_REPORT_HORIZONS <- c(4L, 8L, 13L)

# ---------------------------------------------------------------------------
# Monte-Carlo settings (PLAN §3.4, §6.7)
# ---------------------------------------------------------------------------
# M iterations = draws over BOTH parameter uncertainty (posterior of (beta0,gamma),
# R_eff, k_indiv) AND process stochasticity (branching + Bernoulli seeding). Default
# 1000 gives Monte-Carlo SE < ~0.005 on mid-range reach probabilities; production
# runs use 2000+. Set via env var CASCADE_N_MC to override without editing.
CASCADE_N_MC <- as.integer(Sys.getenv("CASCADE_N_MC",
                  unset = if (identical(Sys.getenv("CASCADE_SMOKE"), "1")) "40" else "2000"))
CASCADE_SEED <- get0("RANDOM_SEED", ifnotfound = 20260704L)
# Nested Monte-Carlo design for CREDIBLE intervals (PLAN §3.4). The M iterations are
# organised as D posterior-parameter draws x n_rep process replicates (M = D x n_rep):
# parameters (hazard coefficients + R_eff) are held fixed within a draw so the credible
# interval reflects POSTERIOR uncertainty in the reach probability (process noise is
# averaged out by variance decomposition in cascade_reach_table). n_rep >= 2 required.
CASCADE_N_REP <- as.integer(Sys.getenv("CASCADE_N_REP",
                  unset = if (identical(Sys.getenv("CASCADE_SMOKE"), "1")) "4" else "8"))

# ---------------------------------------------------------------------------
# Hazard model reused from the Bayesian suite (PLAN §3.3)
# ---------------------------------------------------------------------------
# ASSUMPTION: the cascade's per-week seeding hazard IS the validated cloglog
# invasion hazard from 21_bayesian_renewal.R, p = 1 - exp(-beta_i * Lambda_i),
# beta_i = exp(beta0 + sum_m gamma_m z_{m,i}). We MUST use a COVARIATE-BEARING
# posterior (not the intercept-only featured model) so the dynamically-updated
# d_min frontier covariate modulates the hazard as the front advances.
CASCADE_KERNEL   <- Sys.getenv("CASCADE_KERNEL", unset = "M13-dist")  # CV-featured cohort+gravity road-km
CASCADE_KERNEL_FALLBACK <- "M8"                                       # canonical short-trip composite
CASCADE_GT       <- get0("GT_PRIMARY", ifnotfound = "medium")        # 15.3 d medium; short/long swept
CASCADE_COV_SPEC <- c("log_pop", "ccvi", "d_min")                     # "geo" set: d_min is the dynamic frontier term
CASCADE_LINK     <- "cloglog"

# ---------------------------------------------------------------------------
# Layer A — within-zone transmission (PLAN §3.2)
# ---------------------------------------------------------------------------
# Individual-offspring overdispersion. Because sum of n iid NegBin(mean=R, size=k)
# = NegBin(mean=nR, size=nk), a single draw
#     I_j(w) ~ NegBin(mean = R_j * own_j, size = own_j * k_indiv)
# reproduces individual-level branching (extinction possible when own_j is small,
# the seeding/frontier regime) AND concentrates to a renewal mean once own_j is
# large (established regime) — one parameter, no regime switch. own_j is the
# generation-time-weighted recent incidence (.gweighted_own).
# ASSUMPTION: EBOV is strongly overdispersed (superspreading), k_indiv ~ 0.2-0.4
# (Lloyd-Smith 2005). This is the individual offspring dispersion, NOT a weekly
# aggregate size; it is swept as a sensitivity.
CASCADE_K_INDIV      <- 0.30
CASCADE_K_INDIV_SWEEP <- c(0.20, 0.40)

# Seed size on invasion. ASSUMPTION: an invasion event (first confirmed case) seeds
# a small confirmed index count; governs onward speed. Default 1 + Poisson(0.5).
cascade_draw_seed <- function(n = 1L) 1L + stats::rpois(n, 0.5)
CASCADE_N_SEED_SWEEP <- list(one = function(n) rep(1L, n),
                             pois05 = function(n) 1L + stats::rpois(n, 0.5),
                             pois1  = function(n) stats::rpois(n, 1))

# Establishment threshold: a seeded zone is "established" (a sustained local
# outbreak, not a single imported case) if its cumulative confirmed incidence over
# the horizon reaches n_est. PLAN §1 secondary decision layer (invasion != establishment).
CASCADE_N_EST <- 5L

# Within-zone susceptible depletion. ASSUMPTION: minor over 13 weeks (EBOV rarely
# depletes a whole zone's susceptibles); off by default, most R-decline is carried
# by the scenario control multiplier c(t) below. When on, s_j(t) = S_j/N_j with
# S_j decremented by cumulative infections / rho.
CASCADE_WITHIN_ZONE_DEPLETION <- FALSE

# ---------------------------------------------------------------------------
# Reproduction number (PLAN §3.2, §3.6) — R_eff, NOT R0
# ---------------------------------------------------------------------------
# ASSUMPTION: newly-seeded zones inherit the CURRENT effective reproduction number
# R_eff (~1 for a partially-controlled outbreak ~15 weeks in), NOT the R0-scale
# LogNormal(2.0) prior. Seeding at R0~2 would over-propagate. R_eff is estimated
# in 31_source_dynamics.R (pooled + province partial-pooling); these are bounds/priors.
CASCADE_R_MIN <- 0.30
CASCADE_R_MAX <- 4.0
CASCADE_R_POOL_SD <- 0.30      # log-scale SD of the seeding R_eff predictive prior
CASCADE_R_MIN_CASES <- 10L     # min cases for a zone-specific R_eff (else use pool)

# ---------------------------------------------------------------------------
# Scenarios (PLAN §3.6) — the projection is scenario-based
# ---------------------------------------------------------------------------
# Control multiplier c(t) on R_j(t) = R_j0 * s_j(t) * c(t), t = 1..13 (week index
# into the projection). Pre-registered; NOT fitted. S1 central.
CASCADE_SCENARIOS <- list(
  S1_status_quo = list(
    label = "Status quo",
    desc  = "Current transmissibility persists (c(t)=1); only depletion damps R.",
    c_fun = function(t) rep(1.0, length(t))),
  S2_control = list(
    label = "Strengthened control",
    desc  = "Response scale-up: c(t) decays so R_j crosses 1 by ~week 6 for R_j0~1.3.",
    # exponential decay to a 0.65 floor; halving time ~5 weeks
    c_fun = function(t) 0.65 + 0.35 * exp(-t / 5)),
  S3_deterioration = list(
    label = "Deterioration",
    desc  = "Access/funding/security shock: c(t) rises to +30% by week 13.",
    c_fun = function(t) 1.0 + 0.30 * (t / 13))
)
CASCADE_SCENARIO_PRIMARY <- "S1_status_quo"

# ---------------------------------------------------------------------------
# Calibration of the per-step hazard (PLAN §3.5) — the central correction
# ---------------------------------------------------------------------------
# (1) Frontier saturation sat_i(t) = exp(-psi * f_inv,i(t)), where f_inv,i is the
#     mobility-weighted invaded fraction of i's SOURCE neighbourhood. Self-limits the
#     front. ASSUMPTION: psi >= 0 fitted (33) or set here as a documented default and
#     swept; the mechanism is zone-level frontier saturation, NOT within-zone herd
#     immunity, and is distinct from the shrinking at-risk set.
CASCADE_PSI       <- 1.0
CASCADE_PSI_SWEEP <- c(0.0, 2.0)

# (2) Hazard-scale recalibration delta: mu_tilde = delta * mu (equivalently a shift
#     in beta0). Fitted (33) so the cascade's h=1 seeding matches OBSERVED invasion
#     frequency (calibration_in_large -> 1). Applied ON THE HAZARD so the additive
#     cumulative-hazard identity p = 1 - exp(-sum mu) is preserved. Default here is a
#     placeholder overwritten at runtime by cascade_fit_delta() from the saved
#     invasion_evaluation.csv (delta ~= 1 / calibration_in_large).
CASCADE_DELTA_DEFAULT <- 1.0

# Reads calibration_in_large for the covariate-bearing model (falls back to any
# Bayesian model, then to the report's ~2.3x) and returns delta = 1/cal_in_large,
# floored/capped to a sane band. Documented in PLAN §3.5 / §6.1.
cascade_fit_delta <- function(cov_model = "Bayes-M13-dist-geo",
                              eval_csv = file.path(OUT_DIAGNOSTICS, "invasion_evaluation.csv")) {
  cil <- tryCatch({
    ev <- utils::read.csv(eval_csv, stringsAsFactors = FALSE)
    col <- intersect(c("calibration_in_large", "cal_in_large"), names(ev))[1]
    if (is.na(col)) stop("no calibration column")
    h1 <- ev[ev$horizon == 1L, ]
    v <- h1[[col]][match(cov_model, h1$method)]
    if (!length(v) || !is.finite(v)) v <- stats::median(h1[[col]], na.rm = TRUE)  # any Bayesian model
    v
  }, error = function(e) { message("[cascade] delta: eval csv unavailable (", conditionMessage(e),
                                   "); using ~2.3x report value."); 2.3 })
  delta <- 1 / max(cil, 1e-6)
  delta <- min(max(delta, 0.2), 1.0)   # over-prediction only shrinks the hazard
  attr(delta, "cal_in_large") <- cil
  delta
}

# ---------------------------------------------------------------------------
# Ascertainment (PLAN §1) — confirmed-case scale is primary; rho is secondary
# ---------------------------------------------------------------------------
CASCADE_RHO <- get0("ASCERTAINMENT_NOMINAL", ifnotfound = 0.45)

# ---------------------------------------------------------------------------
# Conditional / attribution (PLAN §3.7)
# ---------------------------------------------------------------------------
CASCADE_KNOCKOUT_TOPN <- 15L   # gateway shortlist size for exact knockout Delta_j
CASCADE_HUBS_FOR_COND <- NULL  # NULL -> auto (current high-incidence source hubs)

# ---------------------------------------------------------------------------
# Sensitivity axes swept in 33 (PLAN §6.6)
# ---------------------------------------------------------------------------
CASCADE_SENS_GT      <- c("short", "medium", "long")
CASCADE_SENS_KERNELS <- c("M13-dist", "M8")

# ---------------------------------------------------------------------------
# Output location
# ---------------------------------------------------------------------------
OUT_CASCADE <- file.path(OUT_DIR, "cascade")
for (.d in c(OUT_CASCADE, file.path(OUT_CASCADE, "figures"),
             file.path(OUT_CASCADE, "tables"), file.path(OUT_CASCADE, "diagnostics")))
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE, showWarnings = FALSE)

message(sprintf("[30_cascade_config] horizon=%dw report=[%s] M=%d kernel=%s cov=[%s] scenarios=%d",
                CASCADE_HORIZON_WEEKS, paste(CASCADE_REPORT_HORIZONS, collapse=","),
                CASCADE_N_MC, CASCADE_KERNEL, paste(CASCADE_COV_SPEC, collapse=","),
                length(CASCADE_SCENARIOS)))
