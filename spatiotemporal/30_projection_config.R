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
# GT_PRIMARY). Depends on nothing else.
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
# M iterations = draws over BOTH parameter uncertainty (posterior of (beta0, gamma), R_eff and
# delta) AND process stochasticity (branching + Bernoulli seeding). k_indiv is a FIXED scalar
# (CASCADE_K_INDIV), swept in sensitivity (41_kindiv_sweep.R) but never drawn, so the credible
# intervals contain NO overdispersion uncertainty — the intervals used to claim otherwise. Default
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
#
# MOBILITY KERNEL: track the LFO-CV-selected model so the cascade's seeding hazard
# stays the *validated* short-horizon model rather than a hardcoded (drift-prone)
# choice. We take the kernel token of the featured Bayesian method from
# key_outputs/model_selection.json and pair it with the -geo covariate spec above
# (the selected variant may be constant-beta `-med`, but the cascade needs the
# covariate-bearing `-geo` twin for the dynamic d_min frontier). Precedence:
#   CASCADE_KERNEL env override > CV-selected kernel > "M13-dist-fill" fallback.
# Extract the kernel token: strip the "Bayes-" prefix and any trailing beta/GT spec
# suffix, KEEPING "-dist" (part of the kernel).
# The suffix set must track bayes_default_grid() (21_bayesian_renewal.R). It did not:
# the time-varying beta families it adds by default are labelled "-tvtrend" / "-tvweek",
# and those were not stripped — so a tv model winning selection yielded a token like
# "M14-tvweek", no mobility_M14-tvweek.rds exists, and the kernel silently degraded to the
# hard-coded fallback instead of the featured model's actual kernel. The parser now lives in
# 00_config.R (mobility_kernel_from_method) and is shared with run_all.R's featured-model
# refits and the report generator, so a suffix added to the grid is handled in ONE place; it
# also validates the token against MOBILITY_IDS. The file-existence guard below is the
# backstop, not the mechanism.
.cascade_kernel_from_method <- function(method) {
  if (!is.character(method) || length(method) != 1L || is.na(method)) return(NA_character_)
  if (!grepl("^Bayes-", method)) return(NA_character_)            # ensembles/renewal -> no single kernel
  mobility_kernel_from_method(method)                             # canonical parser (00_config.R)
}
.cascade_selected_kernel <- function() {
  f <- file.path(OUT_DIR, "key_outputs", "model_selection.json")
  if (!file.exists(f)) return(NA_character_)
  sel <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
  if (is.null(sel)) return(NA_character_)
  m <- sel$featured$bayesian$method
  if (is.null(m) || is.na(m)) m <- sel$featured$headline$method
  k <- .cascade_kernel_from_method(m)
  # only adopt it if its mobility kernel has actually been built on disk
  if (!is.na(k) && !file.exists(file.path(OUT_MOBILITY, sprintf("mobility_%s.rds", k))))
    return(NA_character_)
  k
}
.cascade_kernel_env <- Sys.getenv("CASCADE_KERNEL")
# Filled forms throughout: source-cell fill is the default for every composite (00_config.R).
#
# ONE FALLBACK CONSTANT, ACTUALLY USED. CASCADE_KERNEL's last resort was the LITERAL
# "M13-dist-fill" while CASCADE_KERNEL_FALLBACK — which reads as the fallback, and is what a
# reader would change — was consumed only by run_cascade.R to decide which kernels to
# PRELOAD. The two could therefore disagree, and repointing the constant (as was done when
# MOBILITY_PRIMARY moved to the cohort composite) had no effect on the kernel actually chosen.
# CASCADE_KERNEL now uses it, so there is a single place to change.
CASCADE_KERNEL_FALLBACK <- get0("MOBILITY_PRIMARY", ifnotfound = "M14-fill")
CASCADE_KERNEL <- if (nzchar(.cascade_kernel_env)) {
  .cascade_kernel_env
} else {
  .k_sel <- .cascade_selected_kernel()
  if (is.na(.k_sel)) CASCADE_KERNEL_FALLBACK else .k_sel   # selection artifact missing/unbuilt
}
message(sprintf("[30_cascade_config] mobility kernel = %s (%s)", CASCADE_KERNEL,
                if (nzchar(.cascade_kernel_env)) "env override"
                else if (!is.na(.cascade_selected_kernel())) "LFO-CV selected" else "fallback"))
CASCADE_GT       <- get0("GT_PRIMARY", ifnotfound = "medium")        # 15.3 d medium; short/long swept
# THE GRID'S OWN COVARIATE SET, not a copy of it. This restated c("log_pop", "ccvi", "d_min")
# and called itself the "geo" set. It stopped being the geo set on 2026-09-21, when commit
# 143b5a8 reset BAYES_GEO_COVARIATES to c("ccvi", "d_min") on evidence -- with 54 invasion
# events in 8,623 at-risk rows the honest ceiling is two to three covariates, and log_pop
# alone scored +0.6 AIC on M14-fill, worse than omitting it, because population already
# enters through the kernel and the offset. The cascade was not updated with it, so every
# 13-week number came from a hazard carrying a covariate the cross-validated grid had
# dropped. Reading the constant rather than restating it is what stops that recurring.
# d_min stays the dynamic frontier term: the cascade recomputes it as the front advances.
CASCADE_COV_SPEC <- get0("BAYES_GEO_COVARIATES", ifnotfound = c("ccvi", "d_min"))
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
# (Lloyd-Smith 2005; central 0.30). This is the individual offspring dispersion, NOT
# a weekly aggregate size; it is swept as a sensitivity. The COMPREHENSIVE sweep
# (cascade_kindiv_sweep, module 33/41; figure 42) runs this dense grid — spanning
# heavy superspreading (0.05) to the near-Poisson limit (4.0), densely sampling the
# plausible band — at the fitted (delta, psi). Finding: the invasion RANKING is
# invariant to k (Spearman >= 0.99 across the whole grid) and reach magnitude is
# weakly sensitive, but establishment probability rises ~3x from k=0.05 to k=4 (heavy
# superspreading -> more stochastic fade-out before n_est), which is why k is swept.
CASCADE_K_INDIV      <- 0.30
CASCADE_K_INDIV_SWEEP <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 2.00, 4.00)

# Seed size on invasion. ASSUMPTION: an invasion event (first confirmed case) seeds
# a small confirmed index count; governs onward speed. Default 1 + Poisson(0.5).
# Initial confirmed cases in a newly seeded zone, from UNIFORMS rather than from a count.
# The uniform interface is what makes the seeding paired across a baseline and its contrast:
# zone i draws the same u whether or not the scenario seeded some other city, so the two runs
# differ because of the intervention and not because the random stream moved on.
cascade_draw_seed <- function(u) 1L + stats::qpois(u, 0.5)
CASCADE_N_SEED_SWEEP <- list(one    = function(u) rep(1L, length(u)),
                             pois05 = function(u) 1L + stats::qpois(u, 0.5),
                             pois1  = function(u) stats::qpois(u, 1))

# Establishment threshold: a seeded zone is "established" (a sustained local
# outbreak, not a single imported case) if its cumulative confirmed incidence over
# the horizon reaches n_est. PLAN §1 secondary decision layer (invasion != establishment).
CASCADE_N_EST <- 5L

# Within-zone susceptible depletion was REMOVED (2026-09-19): the flag was FALSE, no caller
# ever overrode it, and its only implementation converted confirmed cases to infections by
# dividing by a guessed ascertainment constant. The modelling ASSUMPTION is unchanged and
# explicit: over 13 weeks EVD does not materially deplete a health zone's susceptibles, and
# the R decline that matters is carried by the scenario control multiplier c(t) below.

# ---------------------------------------------------------------------------
# Reproduction number (PLAN §3.2, §3.6) — R_eff, NOT R0
# ---------------------------------------------------------------------------
# ASSUMPTION: newly-seeded zones inherit the CURRENT effective reproduction number
# R_eff (~1 for a partially-controlled outbreak ~15 weeks in), NOT the R0-scale
# LogNormal(2.0) prior. Seeding at R0~2 would over-propagate. R_eff is estimated
# in 31_source_dynamics.R (pooled + province partial-pooling); these are bounds/priors.
# --- Zone reproduction number: priors, guards, and what is NOT a prior ---------------
# The estimator (31_source_dynamics.R) is a Gamma-conjugate hierarchy, zone <- province <-
# national, on the renewal denominator. It replaces a ratio num/den that was hard-CLAMPED
# to [0.30, 4.0]. That clamp was doing the work of a prior while looking like a guard: at an
# origin four weeks into the outbreak every single zone with cases sat at exactly 4.00 and
# the national estimate sat at exactly its own 5.00 ceiling, so the simulation was driven by
# a constant and nothing in any output said so. The cause is structural — only 18% of the
# weekly generation-time mass sits at a one-week lag, so a zone in its first weeks divides
# by a near-empty denominator, and its first cases are IMPORTED yet were credited to local
# transmission because the denominator had no import term.
#
# kappa is the prior strength in DENOMINATOR units, i.e. "this prior is worth kappa expected
# cases of local transmission". Posterior mean = (kappa * R_parent + cases) / (kappa + expected).
# A zone whose expected local cases are far below kappa is reported as prior-dominated
# rather than silently censored.
# --- WHERE THE CASCADE'S R ACTUALLY COMES FROM ---------------------------------------
# The cascade no longer transmits on the conjugate per-zone estimator below. It takes ONE
# NATIONAL R from the same EpiNow2 posterior the short-term arm uses (bayes_rt_week_draws(),
# 21_bayesian_renewal.R), averaged over the window set here. Three reasons, all of which were
# measured on this outbreak rather than assumed:
#   (i)   A dispersion test found NO detectable between-zone heterogeneity in R (X^2 = 10.6 on
#         14 df, p = 0.71, with negative excess variance at every case threshold). The per-zone
#         and per-province levels were therefore modelling noise, and the three-level shrinkage
#         (kappa_zone / kappa_prov / kappa_nat) was a large apparatus with nothing to estimate.
#   (ii)  Most zones were prior-dominated anyway, so the projection already leaned on the single
#         national number — while presenting itself as per-zone.
#   (iii) EpiNow2 carries an explicit generation time, a right-truncation model and an
#         ascertainment term, and it is what the short-term arm reports. Sharing it removes an
#         unjustifiable difference between the two analyses rather than having to defend one.
# The conjugate estimator is still computed, but as a DIAGNOSTIC (R_nat_conjugate) and for the
# case/import bookkeeping the hazard needs; it no longer drives transmission.
#
# WINDOW. R is averaged over this many weeks ending at the last observed week, WITHIN each
# posterior sample (see bayes_rt_week_draws(): a window mean has to be taken inside a sample).
# The window INCLUDES the final week: the retired conjugate estimator dropped it because
# it is truncation-low and it had no truncation model, whereas EpiNow2 models truncation from
# the same fitted delay the nowcast uses, so dropping it would discard the most current
# information for a reason that no longer applies.
# DERIVED from the single shared constant (00_config.R), currently 1 week: the cascade and the
# 1-2 week invasion forecast must transmit on the same reproduction number. The name is kept
# because several call sites thread it. (The paragraph above used to justify a value of 3 by
# reference to the conjugate estimator's window; this constant has not been 3 since
# RT_WINDOW_WEEKS was introduced, and the conjugate estimator is diagnostic-only.)
CASCADE_R_WINDOW_WEEKS <- get0("RT_WINDOW_WEEKS", ifnotfound = 1L)

# --- GENERATION-TIME UNCERTAINTY: marginalised, not selected -------------------------
# The short-term arm stopped SELECTING a generation time (a reviewer's point: selecting a GT by
# cross-validation is "akin to fitting" a quantity the data cannot identify) and instead places a
# prior over it, GT_PRIOR, and marginalises the posterior over a 15-point (5 mean x 3 sd) grid.
# The cascade ran on a SINGLE profile, which is the inconsistency between the two analyses that
# this flag removes: with it on, each Monte-Carlo parameter group draws a GT grid point with its
# prior weight and uses BOTH that point's weekly kernel and the R posterior fitted AT that GT.
#
# THE PAIRING IS THE POINT. R is estimated under a generation time, so grid point i's R must
# travel with grid point i's kernel; combining one point's R with another's kernel would pair a
# reproduction number with an infectiousness profile it was never estimated under. That is why
# this cannot be approximated by perturbing the kernel alone.
#
# COST: one EpiNow2 fit per grid point per estimation window (the production anchor plus each
# calibration origin), all cached in outputs/diagnostics/rt_draws. Turn OFF to run on
# CASCADE_GT alone, which is the older, narrower behaviour — it understates the intervals,
# because GT uncertainty then contributes nothing.
# DEFAULT FALSE since 2026-09-22: THE GENERATION TIME IS TREATED AS A KNOWN DISTRIBUTION.
#
# Marginalising looks like uncertainty propagation and is not. GT_PRIOR is
# mean ~ Normal(15.3, 1.82) and sd ~ Normal(9.3, 1.50); neither width is estimated from these
# data — they are chosen. Integrating over them produces intervals whose width is set by the
# analyst's prior rather than by evidence, and reports them as if the data had spoken. A
# stated limitation ("our intervals carry no generation-time uncertainty") is honest and
# checkable; an unquantified prior-driven widening is neither.
#
# This was TRUE briefly, on the reasoning that turning it off "understates the intervals
# because GT uncertainty then contributes nothing". That is true and it is the point: the
# contribution it would make is not measured, so we decline to invent it and say so instead.
# Both arms now use their single GT anchor, which also makes the 1-2 week and 13-week analyses
# consistent by construction rather than by matching two marginalisation settings.
#
# The machinery is retained, not deleted: CASCADE_GT_MARGINALISE=1 runs the sensitivity, at a
# cost of one EpiNow2 fit per grid point per estimation window.
CASCADE_GT_MARGINALISE <- {
  v <- tolower(trimws(Sys.getenv("CASCADE_GT_MARGINALISE", "")))
  if (nzchar(v)) v %in% c("1", "true", "t", "yes", "y") else FALSE
}

# --- Conjugate estimator (DIAGNOSTIC ONLY; see above) ---------------------------------
CASCADE_R_KAPPA_ZONE <- 5      # zone prior strength, in expected-case units
CASCADE_R_KAPPA_PROV <- 20     # province prior strength (shrinks province -> national)
CASCADE_R_KAPPA_NAT  <- 1      # national prior strength, centred on CASCADE_R_PRIOR_MEAN
CASCADE_R_PRIOR_MEAN <- 1.0    # neutral anchor for the national level
# NUMERICAL guard only — not a prior, and not a place to encode beliefs. It exists so a
# pathological estimate cannot make the branching process run away; if it ever binds the
# estimator says so loudly and the origin should be excluded, not censored.
CASCADE_R_MIN <- 0.05
CASCADE_R_MAX <- 10.0
# Above this the estimate is flagged as epidemiologically implausible for BDBV (published
# EVD reproduction numbers sit near 1.5-2.5). A flag, never a truncation.
CASCADE_R_PLAUSIBLE_MAX <- 4.0
CASCADE_R_POOL_SD <- 0.30      # FALLBACK log-scale SD when a posterior SD is unavailable
# Retained only as a reporting convention (a zone with fewer than this many cases in the
# window is described as data-poor in text). It is NO LONGER an estimator parameter: the
# old tau0 shrinkage weight n/(n + tau0) was replaced by the Gamma prior above, whose
# strength is kappa_zone in expected-case units. Nothing reads this in 31_source_dynamics.R.
CASCADE_R_MIN_CASES <- 10L

# --- R_eff is NOT constant over the horizon -------------------------------------------
# The projection used to hold each zone's R fixed for all 13 weeks, modulated only by the
# pre-registered scenario multiplier c(t). That made the week-13 credible interval exactly as
# tight as the week-4 one, which cannot be right: the further out you look, the less you know.
#
# WHAT THIS ADDS, AND WHAT IT DELIBERATELY DOES NOT. log R_i evolves as
#     log R_i(t) = log R_i(t-1) + rho * (log R_i(0) - log R_i(t-1)) + sigma * eps
# with rho = 0 since 2026-09-22, i.e. a PURE RANDOM WALK on log R. E[log R_i(t)] stays at the
# estimate, so this adds UNCERTAINTY that grows as sigma*sqrt(t) — it does NOT add a trend.
# That distinction is deliberate. (With rho > 0 the same expression is an Ornstein-Uhlenbeck
# walk whose spread instead plateaus at sigma/sqrt(1-(1-rho)^2); see the rho entry below for
# why that attractor was removed.) The observed national R fell steadily
# (-0.017 to -0.048 per week in log terms over the post-burn-in window), but projecting that
# decline forward is a claim about FUTURE CONTROL, which is precisely what the pre-registered
# c(t) scenarios exist to express. Fitting it here as well would contradict S1 "status quo"
# and double-count control in S2.
#
# SIGMA IS AN ASSUMPTION WITH A SWEEP, not a measurement. An earlier version of this comment
# claimed "sigma is MEASURED: an AR(1) on the weekly national log R series gives a residual
# innovation SD of 0.036-0.040 across windows beginning at weeks 8, 10 and 12". That does not
# reproduce. Refitting on the current R(t) series (v8 cache, 18 weekly points, 2026-05-04 to
# 2026-08-31) gives residual SDs of 0.069, 0.026 and 0.029 at those three windows — spanning
# the claimed range rather than sitting inside it.
#
# The reason is structural, and it applies to any estimate of this kind here: national R over
# the window is a monotone decline (3.24 -> 1.13), not fluctuation about a level, so an AR(1)
# in deviations from an anchor is being asked to separate trend from innovation on 7-11
# points. It cannot. Two further problems compound it — the series is EpiNow2's posterior
# MEAN, already smoothed by the model's own random walk on R, so an innovation SD read off it
# is biased low by construction; and 90 of the 131 historical days are "estimate based on
# partial data", so the cited windows sit entirely inside the truncation-corrected region.
#
# 0.04 is therefore retained as the CENTRAL member of the pre-specified sweep below, not as a
# fitted value. The sweep is the uncertainty statement. For scale: EpiNow2's own posterior SD
# on log R is about 0.099 per week (measured across the weekly rt_draws caches), so the walk
# is a second-order contributor to R uncertainty, not the dominant one.
CASCADE_R_RW_SIGMA <- 0.04
# REVERSION: NOW ZERO (2026-09-22). rho > 0 pulls log R back toward its value at the
# projection origin, which asserts that R returns to whatever it happened to be on the day we
# projected from. Nothing supports that. Refitting on the current series gives rho = -0.49,
# -0.91 and -0.88 at the windows the retired comment cited as "0.32-0.56" — negative, i.e. the
# opposite of reversion, because a declining series read as deviations from its own first
# value moves steadily away from it.
#
# rho = 0 makes this a pure random walk on log R: the minimal statement that week-to-week
# changes are unpredictable, with no claim about direction or attractor. It is the weaker
# assumption, and the honest one over 13 weeks. Runaway is not a risk — the simulator clamps
# R to [CASCADE_R_MIN, CASCADE_R_MAX] on every step (32_cascade_simulator.R).
#
# CONSEQUENCE, STATED: the fan no longer plateaus. At sigma = 0.04 the walk's SD on log R
# grows as sigma*sqrt(t), reaching 0.144 by week 13 against the 0.056 the rho = 0.30 walk
# plateaued at. Combined with the R posterior (0.099) that widens total R uncertainty from
# about 0.114 to 0.175. The 13-week intervals are correspondingly wider, and that is the
# point: they were previously narrowed by an attractor we cannot justify.
CASCADE_R_RW_RHO   <- 0
CASCADE_R_RW_SWEEP <- c(0, 0.02, 0.04, 0.08)   # sigma sweep; 0 reproduces constant-R exactly

# --- Nowcast uncertainty in R_eff: RETIRED ---------------------------------------------
# This constant is READ BUT IGNORED. cascade_reff() declares the argument retired
# (31_source_dynamics.R) and emits a runtime message saying so.
#
# It used to describe a bootstrap that resampled the weekly completeness factor and added the
# induced spread to sd_zone in quadrature. That is no longer how R is obtained: R comes from
# EpiNow2, which models right-truncation internally from the same fitted delay, so the
# correction's uncertainty is already inside the posterior rather than something to be added
# afterwards. Doing both would double-count it. The paragraph that used to sit here described
# the bootstrap in the present tense ("cascade_reff() NOW resamples...") — a statement of fact
# the code does not produce. It also referenced epinowcast's nowcast_cv, which the deployed
# deterministic nowcast sets to NA everywhere, so the quantity did not exist on the live path.
#
# Kept only so existing callers that pass it still resolve; it has no effect.
CASCADE_NOWCAST_BOOT <- 100L

# --- Expected imported cases in the renewal denominator ------------------------------
# Expected weekly introductions into zone i are delta * beta0 * Lambda_i, and each
# introduction brings E[n_seed] confirmed cases. Counting them in the denominator is what
# stops a newly invaded zone's imported cases being read as local transmission.
CASCADE_R_USE_IMPORTS <- TRUE
# NULL means "derive E[n_seed] from the seeding function actually in force", which is the
# correct behaviour: it was the literal 1.5 = E[1 + Poisson(0.5)], which silently became
# wrong whenever CASCADE_N_SEED_SWEEP swapped in "one" (1.0) or "pois1" (1.0). Set a number
# only to override deliberately.
CASCADE_R_IMPORT_CASES_PER_INTRO <- NULL

# --- Seeded-city reproduction number: an ASSUMPTION, made explicit --------------------
# A zone seeded during the projection inherits its province pool, and for a province with no
# case history that pool IS the national value. The urban scenarios therefore ask "what if a
# capital is seeded and then behaves like the current national average", and answer, quite
# correctly, that little happens: at a national R near 1 with k = 0.30 a one-to-two case seed
# usually fades out. That is a premise, not a finding, so it is swept rather than assumed.
# NULL/NA = use the province pool (previous behaviour). (The measured fade-out probabilities
# once quoted here were computed at R = 1.19, a value the deployed anchor has long since left.)
#
# REGIME NOTE, measured rather than assumed. The within-zone layer has no saturating
# mechanism at all (within-zone depletion is not modelled), so a zone with R above 1 compounds
# for the whole 13 weeks and it is fair to ask whether the upper arms leave the regime the
# model can represent. Measured at psi = 2, M = 400 (seed 11): expected new invasions over
# the horizon were 56.8 with no urban seeding, 62.6 seeding Kinshasa at the pooled R, and
# 67.5 at R = 2.5 — so the sweep roughly doubles the attributable effect (+5.7 -> +10.7)
# without the projection running away, and no branching mean came near mu_max. Switching
# within-zone depletion ON changed those same three numbers by about 2% (55.6 / 61.3 / 66.3),
# which is the first empirical check of this file's long-standing "minor over 13 weeks"
# claim. Within-zone depletion has since been REMOVED entirely (see the removal note earlier in
# this file and 32_cascade_simulator.R), so those "+2% with depletion ON" figures can no longer
# be reproduced and nothing is armed: pop_vec is still threaded through so every contrast calls
# simulate_cascade() identically, but the simulator accepts and ignores it. mu_max in
# cascade_branch_step() is a backstop against a future parameterisation, not a live brake.
CASCADE_SEED_R <- NA_real_
CASCADE_SEED_R_SWEEP <- c(NA, 1.5, 2.0, 2.5)   # NA = the pooled value, as one arm

# --- Common random numbers ------------------------------------------------------------
# Every scenario contrast in this layer (urban seeding, gateway knockout, conditional hub
# seeding, the out-of-sample delta search) is a DIFFERENCE between two simulations. Without
# pairing, that difference is dominated by Monte-Carlo noise: the 2026-09-07 urban run reported
# a NEGATIVE expected number of added invasions for most cities, and the zones it credited with
# the largest increase from seeding Kinshasa were beside the actual outbreak 1,500 km away.
#
# HOW TO READ A NEGATIVE NOW THAT psi = 0. Under the retired saturation term a small negative
# could be genuine — frontier saturation was a real negative feedback. With psi fixed at 0 that
# channel is gone, so seeding a city can only ADD import force, and any negative attributable
# change is Monte-Carlo noise or the (small) re-seeding channel, never protection. Negatives
# should therefore no longer be explained away as saturation; if they are large, the pairing is
# not working. Set FALSE only to reproduce pre-2026-09 outputs.
CASCADE_CRN <- {
  v <- tolower(trimws(Sys.getenv("CASCADE_CRN", "")))
  if (nzchar(v)) v %in% c("1", "true", "t", "yes", "y") else TRUE
}

# ---------------------------------------------------------------------------
# Scenarios (PLAN §3.6) — the projection is scenario-based
# ---------------------------------------------------------------------------
# Control multiplier c(t) on R_j(t) = R_j0 * s_j(t) * c(t), t = 1..13 (week index
# into the projection). Pre-registered; NOT fitted. S1 central.
CASCADE_SCENARIOS <- list(
  S1_status_quo = list(
    label = "Status quo",
    # NOTHING damps R in S1: c(t) = 1 for every week, the frontier-saturation term is off
    # (CASCADE_PSI = 0) and within-zone susceptible depletion was removed. This is the
    # authoritative in-code definition of the PRIMARY specification, so it must not name a
    # mechanism the simulator does not have.
    desc  = "Current transmissibility persists: c(t) = 1 for all 13 weeks, R undamped.",
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
# (1) Frontier saturation. The IMPLEMENTED form is INCREMENTAL:
#         sat_i(t) = exp(-psi * max(f_inv,i(t) - f_inv,i(0), 0))
#     where f_inv,i is the mobility-weighted invaded fraction of i's SOURCE neighbourhood.
#     This comment previously documented sat = exp(-psi * f_inv), WITHOUT the baseline
#     subtraction — which is a different model. The increment is deliberate and load-bearing:
#     it makes sat = 1 in week 1 for every zone, so the cascade's week-1 ranking is identical
#     to the validated short-horizon hazard and the h=1 consistency gate means what it says.
#     Saturation then bites only as NEW zones are seeded (32_cascade_simulator.R).
#     ASSUMPTION: psi >= 0, fitted out-of-sample (33b) or set here as a documented default and
#     swept; the mechanism is zone-level frontier saturation, NOT within-zone herd immunity,
#     and is distinct from the shrinking at-risk set.
#     RETIRED AS A FITTED PARAMETER, and now fixed at 0 (no frontier saturation).
#     psi was never identified on this outbreak. Across runs the fit came back censored at a
#     bound, or 30% away from its own root, and — the decisive evidence — the two probabilistic
#     criteria (deviance and AUC-PR) kept improving all the way to the search ceiling while the
#     count match sat far below it. That is the documented signature of a phenomenological term
#     absorbing misspecification it cannot represent (most likely the kernel's dispersion), not
#     of a saturation effect being measured. Carrying an unidentifiable parameter into a
#     published 13-week projection is worse than not having it, so the cascade now fits ONE
#     parameter, delta, against realised K-week counts (cascade_fit_delta_oos(), 33b).
#     Setting psi = 0 makes sat identically 1, so the hazard is mu = delta * exp(eta) and the
#     multi-week shape is whatever the mobility kernel and R imply — a stated consequence, not a
#     hidden one: the cascade no longer has a free term to bend that shape.
CASCADE_PSI       <- 0
# Kept as a SENSITIVITY only, to show what turning saturation on would do. The baseline is now 0,
# so a 0 arm would merely duplicate it; these are the two arms the earlier fits wandered between.
CASCADE_PSI_SWEEP <- c(1.0, 2.0)
# Relative width of the delta bracket, on the LOG scale, at which the out-of-sample search stops
# (log(hi/lo) <= this). delta is a multiplicative factor, so its bracket is naturally a ratio;
# 0.01 is a 1% bracket. As with psi, this is NUMERICAL width and not delta's uncertainty — the
# Monte-Carlo component dominates it — and the report must not quote it as a precision.
CASCADE_DELTA_BRACKET_TOL <- 0.01
# Search bounds for cascade_fit_psi(). psi was previously chosen from a FIXED grid whose
# maximum was 6; in the 2026-08 run the fit landed exactly on that ceiling (a censored
# boundary solution returned with no convergence flag, because which.min(|curve - target|)
# cannot distinguish "matched at 6" from "ran out of grid at 6"). The solver now brackets
# and bisects on a continuous interval and reports whether it converged INSIDE it.
# CASCADE_PSI_MAX is deliberately far above any plausible value: hitting it is a
# diagnostic that the saturation term cannot absorb the kernel's dispersion, not a fit.
CASCADE_PSI_MAX   <- 40
# Relative width the psi BRACKET must reach before the search stops. The count tolerance
# (CASCADE_PSI_TOL) alone stops at the first psi whose modelled count lands inside the band,
# which is not the estimator's value: on the 2026-09 frame it returned the first bisection
# midpoint, 7.20, when the zero crossing was near 5.5. Refining on the bracket pins psi to
# +/-2%, and the bracket is returned so the residual numerical width is visible.
CASCADE_PSI_BRACKET_TOL <- 0.02
CASCADE_PSI_TOL   <- 0.05   # relative tolerance on the modelled-vs-target cumulative

# (2) Hazard-scale recalibration delta: mu_tilde = delta * mu (equivalently a shift
#     in beta0). Applied ON THE HAZARD so the additive cumulative-hazard identity
#     p = 1 - exp(-sum mu) is preserved. Default here is a placeholder overwritten at
#     runtime by cascade_fit_delta().
CASCADE_DELTA_DEFAULT <- 1.0

# THE ESTIMATOR (changed 2026-09-11; see 33b_cascade_calibration.R for the full rationale).
# delta is the MAXIMUM-LIKELIHOOD hazard-scale factor that 16b_invasion_recalibration.R
# already fits PREQUENTIALLY for the cascade's own covariate-bearing model, read from
# invasion_recalibration.csv. It is the same transform family the cascade applies
# (mu -> delta*mu is exactly p -> 1-(1-p)^delta), so the two are directly interchangeable.
#
# It replaces the previous moment estimator delta = 1/calibration_in_large clamped to
# [0.2, 1.0], which was wrong in three separate ways:
#   * the moment estimator is the one 16b documents as WORSE on this frame — a rare-event
#     mean is dominated by the many near-zero cells, and on the 2026-09-07 run it returns
#     0.507 for Bayes-M14-geo against a maximum-likelihood 0.445, i.e. ~14% more seeding
#     hazard per week, compounded over 13 weekly steps;
#   * the [0.2, 1.0] clamp cannot represent delta > 1, so an UNDER-predicting hazard could
#     never be corrected — it would sit silently on the boundary. The search band is now
#     RECAL_BAND, which spans both sides of 1, and a boundary hit is reported;
#   * it carried no uncertainty at all, while the cascade draws (beta0, gamma) and R_eff
#     from their posteriors. delta's zone-clustered interval now travels with it
#     (attr "se_log") and simulate_cascade() draws it per parameter group.
#
# Fallbacks degrade in order and each one says so: the 16b table -> the same kernel's
# intercept-only twin -> the median over Bayesian models in that table -> the legacy moment
# estimator from invasion_evaluation.csv (band-clamped, LOUDLY warned) -> 1/2.3 from the report.
cascade_fit_delta <- function(cov_model = "Bayes-M13-dist-fill-geo",
                              eval_csv  = file.path(OUT_DIAGNOSTICS, "invasion_evaluation.csv"),
                              recal_csv = file.path(OUT_DIAGNOSTICS, "invasion_recalibration.csv"),
                              band      = get0("RECAL_BAND", ifnotfound = c(0.01, 5.0))) {
  .out <- function(d, estimator, source, cil = NA_real_, lo = NA_real_, hi = NA_real_,
                   se_log = NA_real_, boundary = NA) {
    # RECORD THE CLAMP. This line silently pins d into `band`, but `boundary` was taken only
    # from the 16b table and left NA on every fallback path — so a legacy-moment fallback with
    # 1/cal_in_large > band[2] came back as exactly band[2] with boundary_hit = NA, and
    # cascade_calibration.json recorded a censored value as an ordinary estimate. Whatever the
    # caller passed still wins (it knows about its own optimiser boundary); the clamp only ever
    # UPGRADES an otherwise-unflagged value to TRUE.
    .raw <- as.numeric(d)
    d <- min(max(.raw, band[1]), band[2])
    .clamped <- is.finite(.raw) && !isTRUE(all.equal(.raw, d))
    if (.clamped && !isTRUE(boundary)) {
      boundary <- TRUE
      warning(sprintf(paste0("[cascade] delta from '%s' was %.4g, outside the band [%.3g, %.3g]; ",
                             "clamped to %.4g and flagged boundary_hit = TRUE."),
                      estimator, .raw, band[1], band[2], d), call. = FALSE)
    }
    attr(d, "estimator")    <- estimator
    attr(d, "source")       <- source
    attr(d, "cal_in_large") <- cil
    attr(d, "lo")           <- lo
    attr(d, "hi")           <- hi
    attr(d, "se_log")       <- se_log
    attr(d, "boundary_hit") <- boundary
    d
  }

  # ---- 1. the prequential ML factor from 16b -------------------------------
  rc <- if (file.exists(recal_csv))
          tryCatch(utils::read.csv(recal_csv, stringsAsFactors = FALSE), error = function(e) NULL)
        else NULL
  if (!is.null(rc) && all(c("method", "horizon", "delta") %in% names(rc))) {
    r1 <- rc[rc$horizon == 1L, , drop = FALSE]
    pick <- function(m) { i <- match(m, r1$method)
                          if (!is.na(i) && is.finite(r1$delta[i]) && r1$delta[i] > 0) i else NA_integer_ }
    i <- pick(cov_model)
    used <- cov_model
    if (is.na(i)) {                                  # same kernel, intercept-only twin
      # TRY BOTH SPELLINGS. The intercept-only twin of "<kernel>-geo" is sometimes
      # "<kernel>-med" and sometimes the BARE kernel: the -dist-fill family is labelled
      # "Bayes-M13-dist-fill", with no "-med". Substituting "-med" unconditionally produced
      # "Bayes-M13-dist-fill-med", a label that exists in no evaluation table, so this entire
      # fallback tier was unreachable for exactly the default kernel family.
      .base <- sub("-(geo|full)$", "", cov_model)
      alt <- NA_character_
      for (.cand in c(paste0(.base, "-med"), .base)) {
        .i <- pick(.cand)
        if (!is.na(.i)) { alt <- .cand; break }
      }
      i <- if (is.na(alt)) NA_integer_ else pick(alt)
      if (!is.na(i)) {
        used <- alt
        warning(sprintf(paste0("[cascade] delta: '%s' is absent from %s; using the same kernel's ",
                               "intercept-only twin '%s'. The covariate model's own level may differ."),
                        cov_model, basename(recal_csv), alt), call. = FALSE)
      }
    }
    if (!is.na(i)) {
      g <- function(nm) if (nm %in% names(r1) && is.finite(r1[[nm]][i])) r1[[nm]][i] else NA_real_
      return(.out(r1$delta[i], "maximum_likelihood_16b", used,
                  cil = g("cal_in_large"), lo = g("delta_lo"), hi = g("delta_hi"),
                  se_log = g("se_log_delta"),
                  boundary = if ("boundary_hit" %in% names(r1)) isTRUE(as.logical(r1$boundary_hit[i])) else NA))
    }
    bay <- r1[grepl("^Bayes", r1$method) & !grepl("-ens-", r1$method) & is.finite(r1$delta), , drop = FALSE]
    if (nrow(bay)) {
      warning(sprintf(paste0("[cascade] delta: neither '%s' nor its intercept-only twin is in %s; ",
                             "using the MEDIAN factor over %d Bayesian models. This is a fallback, ",
                             "not a fit for this kernel."),
                      cov_model, basename(recal_csv), nrow(bay)), call. = FALSE)
      return(.out(stats::median(bay$delta), "median_over_models_16b", "median(Bayesian)",
                  cil = NA_real_))
    }
  }

  # ---- 2. legacy moment estimator (loudly) ---------------------------------
  cil <- tryCatch({
    ev <- utils::read.csv(eval_csv, stringsAsFactors = FALSE)
    col <- intersect(c("calibration_in_large", "cal_in_large"), names(ev))[1]
    if (is.na(col)) stop("no calibration column")
    h1 <- ev[ev$horizon == 1L, ]
    v <- h1[[col]][match(cov_model, h1$method)]
    if (!length(v) || !is.finite(v)) v <- stats::median(h1[[col]], na.rm = TRUE)
    v
  }, error = function(e) {
    # NO SILENT CONSTANT. This used to fall back to a literal 2.3 described as "the ~2.3x
    # report value" — a calibration factor mirrored from a report that is regenerated every
    # run, so it went stale the moment the record grew, and it entered the projection with
    # only a message(). If neither the prequential factor nor the moment estimator is
    # available there is no evidence for any delta, and the caller must say so rather than
    # invent one.
    warning("[cascade] delta: evaluation table unavailable (", conditionMessage(e),
            ") and no prequential factor on disk, so NO calibration factor can be derived. ",
            "Re-run 16_invasion_eval.R / 16b_invasion_recalibration.R before the cascade.",
            call. = FALSE)
    NA_real_
  })
  # cil is NA when even the evaluation table was unreadable (see the handler above), and
  # 1/NA is NA — so stop rather than hand simulate_cascade() an NA delta it would silently
  # propagate into every projected probability.
  if (!is.finite(cil) || cil <= 0)
    stop("[cascade] delta: neither the prequential maximum-likelihood factor nor the moment ",
         "estimator could be obtained, so the cascade has no calibration factor. Run ",
         "16_invasion_eval.R and 16b_invasion_recalibration.R against this frame first.",
         call. = FALSE)
  warning(sprintf(paste0("[cascade] delta: no maximum-likelihood factor available (%s missing or ",
                         "empty) — falling back to the MOMENT estimator 1/calibration_in_large = ",
                         "%.3f. That estimator over-states the hazard on this frame and carries no ",
                         "interval; run the main pipeline with INVASION_RECALIBRATE=TRUE first."),
                 basename(recal_csv), 1 / max(cil, 1e-6)), call. = FALSE)
  .out(1 / max(cil, 1e-6), "moment_1_over_cal_in_large", "invasion_evaluation.csv", cil = cil)
}

# ---------------------------------------------------------------------------
# Conditional / attribution (PLAN §3.7)
# ---------------------------------------------------------------------------
CASCADE_KNOCKOUT_TOPN <- 15L   # gateway shortlist size for exact knockout Delta_j
CASCADE_HUBS_FOR_COND <- NULL  # NULL -> auto (current high-incidence source hubs)

# ---------------------------------------------------------------------------
# Sensitivity axes swept in 33 (PLAN §6.6)
# ---------------------------------------------------------------------------
CASCADE_SENS_GT      <- c("short", "medium", "long")
# Alternative kernels for the sensitivity axis: the two other principal cohort/composite
# kernels (cascade_sensitivity excludes the selected CASCADE_KERNEL via setdiff, so
# whichever is selected the other two serve as the robustness comparators).
# Filled forms, since source-cell fill is the default for every composite (00_config.R).
CASCADE_SENS_KERNELS <- c("M14-fill", "M13-dist-fill", "M8-fill")

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
