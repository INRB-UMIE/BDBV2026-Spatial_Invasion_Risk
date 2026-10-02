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
# Sourced AFTER 06_simple_models.R (compute_foi/cascade_weekly_gt),
# 15_workhorse.R (.gweighted_own, estimate_R_local, .count_wide) and
# 30_projection_config.R. ASSUMPTIONS documented inline (Task 3).
# =============================================================================

# ---------------------------------------------------------------------------
# Per-zone current effective reproduction number, R_eff (NOT R0).
#
# THE ESTIMATOR, in one paragraph. Over the last `window` reporting-reasonably-complete
# weeks (the final week is excluded — truncation-low even after nowcasting), zone i's
# confirmed cases are modelled as Poisson with mean R_i * loc_i + imp_i: a local renewal
# term and an imported term built from the same mobility force and deployed hazard scale the
# cascade itself simulates on. Cases are attributed between the two by EM and R_i is the
# Gamma-conjugate posterior mean under a prior centred on the zone's province pool, the
# province shrunk in turn to the national estimate. Prior strengths are in EXPECTED-CASE
# units, so "kappa_zone = 5" means "this prior is worth five expected local cases".
#
# WHAT IT REPLACED. A renewal ratio, empirically shrunk toward a province mean with weight
# n/(n + tau0), and then HARD-CLAMPED to [0.30, 4.0]. The clamp was doing a prior's job
# while looking like a guard: four weeks into this outbreak every zone with cases sat at
# exactly 4.00 and the national estimate at its own ceiling, so the projection ran on a
# constant and nothing in the output said so. The shrinkage is now an explicit prior that
# is reported (`prior_dominated`, `import_share`), and the only remaining bound is a wide
# numerical guard that warns when it binds.
#
# See estimate_zone_reff()'s own header for the model, its two stated limitations
# (attribution circularity inside a connected cluster; the import-conversion
# extrapolation), and the full return contract. cascade_reff() is the supported entry point
# — call that, not estimate_zone_reff() directly, so that production and every calibration
# origin build R identically.
#
# @param zone_week_nc  nowcast-corrected zone-week tibble (needs confirmed_nc)
# @param zones_all     519-zone spine
# @param gt_pmfs,gt    generation-time PMFs + profile key
# @param zone_province named char vector zone -> province (optional; NULL = national only)
# @param window        # of trailing weeks in the estimation window (default 3)
# ---------------------------------------------------------------------------
# Expected imported cases, for the renewal denominator
# ---------------------------------------------------------------------------

#' Expected imported CONFIRMED cases per zone-week.
#'
#' Expected weekly introductions into zone i are delta * beta0 * Lambda_i(w), where Lambda is
#' the mobility import force the invasion model already uses, delta is the deployed hazard
#' recalibration factor, and each introduction brings E[n_seed] confirmed cases. Counting
#' these in the renewal decomposition is what stops a newly invaded zone's IMPORTED cases
#' being read as local transmission — a large part of the early-window inflation (Bunia's raw
#' ratio ran 6.1 and 7.6 in weeks 3 and 4 against an eventual 1.0-1.1).
#'
#' compute_foi() reads strictly past weeks, so column w is the force acting INTO week w and
#' carries no contemporaneous information.
#'
#' THE KNOWN LIMITATION, stated because it is a modelling decision and not an oversight.
#' Within a tightly connected cluster the attribution is partly circular: Rwampara's cases
#' generate import force back into Bunia, so part of what is really one epidemic gets booked
#' as importation and discounted from Bunia's local R. The capping at `num` bounds this, and
#' `import_share` reports its size (20.2% of the national denominator on this frame, moving
#' Bunia 1.11 -> 0.89), but it is not eliminated. Read the zone R's of the epicentre cluster
#' as local-transmission-net-of-modelled-importation, not as closed-population R's.
#'
#' @param delta deployed hazard recalibration factor; beta0 is the RAW coefficient and the
#'   cascade deploys delta*beta0, so introductions must carry the same factor.
#' @param n_seed_fun the seeding function actually in force; E[n_seed] is derived from it.
#' @return zones x weeks matrix of expected imported cases, or NULL when the inputs are
#'   unavailable (the estimator then degrades to local-only and says so).
#' Per-zone import-to-case coefficient beta_i, at the posterior mean.
#'
#' The invasion hazard is exp(eta0_i + gamma_dmin z_dmin_i(w) + log Lambda_i), so a
#' zone-specific conversion from import force to expected introductions is exp(eta0_i).
#'
#' NOT THE DEFAULT, and the reason is a measured one. On this frame the fitted log-population
#' coefficient is NEGATIVE (-0.29, 90% CI -0.57 to -0.02): conditional on the mobility force,
#' a larger zone has a LOWER per-unit invasion hazard. That is what a cloglog link does when
#' the offset is large — p = 1 - exp(-beta*Lambda) saturates towards 1 for the high-Lambda
#' zones, which are exactly the big cities, so the fit compensates with a negative population
#' slope. Carrying that compensation into a LINEAR count model would systematically
#' under-credit importation into large cities (Bunia's beta_i is 0.059 against a
#' population-average 0.216, a factor of 3.7), and it would do so because of a link artefact
#' rather than because of anything epidemiological. The scalar is therefore primary and this
#' is the documented sensitivity — note it moves the answer the FAVOURABLE way (import share
#' 6.0% instead of 11.4%, Bunia 1.08 instead of 1.01), which is a further reason not to
#' adopt it as the default.
#'
#' The dynamic d_min term is omitted regardless: it is defined relative to the invasion FRONT
#' among at-risk zones, and these are zones that have already been invaded. Its covariate is
#' centred, so dropping it means "at the average frontier position".
#'
#' EXTRAPOLATION, stated plainly and true of BOTH specifications: the hazard was fitted on
#' FIRST invasions of susceptible zones, and is applied here to continuing importation into
#' already-invaded zones. It is the best available estimate of the mobility-to-case
#' conversion, not a fitted quantity for this use; `import_share` reports how much work it
#' is doing (6-11% of the national denominator on this frame).
cascade_import_beta <- function(fit, design, covariates, zones_all,
                                cov_spec = CASCADE_COV_SPEC, covariates_on = TRUE) {
  dr <- tryCatch(posterior::as_draws_df(fit), error = function(e) NULL)
  if (is.null(dr) || !("b_Intercept" %in% names(dr))) return(NULL)
  if (!isTRUE(covariates_on))       # intercept only: the conversion at the covariate mean
    return(stats::setNames(rep(exp(mean(dr[["b_Intercept"]])), length(zones_all)), zones_all))
  gm <- function(nm) if (paste0("b_", nm) %in% names(dr) && nm %in% cov_spec)
                       mean(dr[[paste0("b_", nm)]]) else 0
  # Standardise with the FIT design's center/scale, exactly as cascade_prepare() does —
  # design$feat is a vector of NAMES, not values, so the raw covariates come from
  # .static_features() and the two paths must not drift apart.
  static <- .static_features(covariates, zones_all)
  z_of <- function(nm) {
    c0 <- design$center[[nm]]; s0 <- design$scale[[nm]]; raw <- static[[nm]]
    # Each guard is checked for LENGTH before value: `!is.finite(NULL)` is logical(0), and
    # `||` on a zero-length operand is an error in R >= 4.3, so an absent scale would abort
    # rather than fall back.
    if (length(c0) != 1L || length(s0) != 1L || is.null(raw) ||
        !is.finite(s0) || s0 <= 0) return(rep(0, length(zones_all)))
    v <- (as.numeric(raw) - c0) / s0; v[!is.finite(v)] <- 0; v
  }
  b <- exp(mean(dr[["b_Intercept"]]) + gm("log_pop") * z_of("log_pop") + gm("ccvi") * z_of("ccvi"))
  b[!is.finite(b)] <- 0
  stats::setNames(pmax(b, 0), zones_all)
}

cascade_import_cases <- function(zone_week_nc, zones_all, mobility_matrices, gt_pmfs,
                                 kernel = CASCADE_KERNEL, gt = CASCADE_GT, beta0 = NULL,
                                 delta = 1, n_seed_fun = cascade_draw_seed,
                                 cases_per_intro = NULL) {
  W <- mobility_matrices[[kernel]]
  # beta0 may be a scalar (the intercept) or a per-zone vector from cascade_import_beta().
  if (is.null(W) || is.null(beta0) || !any(is.finite(beta0)) || max(beta0, na.rm = TRUE) <= 0)
    return(NULL)
  b_i <- as.numeric(beta0)
  if (length(b_i) == 1L) b_i <- rep(b_i, length(zones_all))
  if (length(b_i) != length(zones_all))
    stop("[import] beta0 must be scalar or one value per zone.", call. = FALSE)
  b_i[!is.finite(b_i) | b_i < 0] <- 0
  # E[n_seed] is DERIVED from the seeding function actually in force rather than hard-coded.
  # It was the literal 1.5 = E[1 + Poisson(0.5)], which silently becomes wrong the moment
  # CASCADE_N_SEED_SWEEP swaps in "one" or "pois1".
  if (is.null(cases_per_intro)) {
    cases_per_intro <- tryCatch(mean(n_seed_fun(seq(0.0005, 0.9995, length.out = 2000))),
                                error = function(e) NA_real_)
    if (!is.finite(cases_per_intro) || cases_per_intro <= 0) cases_per_intro <- 1.5
  }
  # beta0 is the RAW invasion-hazard coefficient. The cascade deploys it multiplied by the
  # recalibration factor delta (~0.45 on this snapshot, i.e. the raw hazard over-predicts
  # introductions about two-fold), so the expected number of introductions must carry the
  # same factor. Leaving it out inflated modelled importation by 1/delta and therefore
  # over-subtracted from every zone's local transmission.
  d <- as.numeric(delta)
  if (!length(d) || !is.finite(d) || d <= 0) d <- 1
  Y <- .count_wide(zone_week_nc, zones_all, "confirmed_nc")
  G <- cascade_weekly_gt(gt, gt_pmfs)
  imp <- matrix(0, nrow(Y), ncol(Y), dimnames = dimnames(Y))
  for (w in seq_len(ncol(Y)))
    imp[, w] <- pmax(d * b_i * compute_foi(Y, W, G, t_idx = w, zones_all) * cases_per_intro, 0)
  imp
}

# ---------------------------------------------------------------------------
# Zone-level effective reproduction number (Gamma-conjugate hierarchy)
# ---------------------------------------------------------------------------

#' Partially pooled R_eff per zone, with importation in the decomposition.
#'
#' MODEL. For zone i in week w, observed confirmed cases are
#'     Y_iw ~ Poisson( R_i * loc_iw  +  imp_iw ),
#'       loc_iw = sum_k G[k] Y_i,w-k       (local renewal force)
#'       imp_iw = expected imported cases  (cascade_import_cases(); 0 when unavailable)
#' with a Gamma prior on R_i whose MEAN is the zone's province pool and whose strength
#' `kappa_zone` is in denominator units — "worth kappa expected local cases". The posterior
#' mean is obtained by EM: cases are attributed between the local and imported components in
#' proportion to their current intensities (E-step), then the conjugate Gamma-Poisson
#' posterior mean is taken over the local part (M-step). With imp = 0 this collapses exactly
#' to the plain conjugate posterior mean; with loc = 0 it returns the prior, which is the
#' honest answer for a zone with no local history.
#'
#' WHAT THIS REPLACES, and why it matters for a published projection. The previous estimator
#' was num/den hard-CLAMPED to [0.30, 4.0]. The clamp was doing a prior's job while looking
#' like a guard: at an origin four weeks into the outbreak all ten zones with cases sat at
#' exactly 4.00 and the national estimate at exactly its own 5.00 ceiling, so the simulation
#' ran on a constant and no output said so. Here there is no clamp — only a wide numerical
#' guard that warns if it ever binds — and the shrinkage is an explicit, reportable prior.
#'
#' SCALE OF THE IMPORT CORRECTION. Imports enter at the ZONE level only. Nationally, mobility
#' importation is internal redistribution rather than new infection, so the national and
#' province aggregates are formed from raw totals; at province level that is an approximation
#' (it ignores importation across provincial borders) and is recorded as such.
#'
#' WHAT IT STILL CANNOT DO. In the first weeks of a newly detected outbreak, apparent
#' transmission is dominated by ASCERTAINMENT ramp-up — a response team arrives and a backlog
#' appears at once — which no renewal estimator can separate from transmission without an
#' explicit reporting model. That is handled by refusing to calibrate on origins inside the
#' ramp (`min_weeks_from_start`, 33b_cascade_calibration.R), not by censoring here.
#'
#' Nor can it undo the CIRCULARITY inside a tightly connected cluster: a neighbour's cases
#' generate import force back into the index zone, so part of one epidemic is booked as
#' importation. Capping `imp` at the observed count bounds it and `import_share` reports its
#' size, but the epicentre's values are local transmission NET OF MODELLED IMPORTATION, not
#' closed-population reproduction numbers, and must be described that way.
#'
#' @param import_cases zones x weeks expected imported cases, or NULL for local-only.
#' @param kappa_zone,kappa_prov,kappa_nat prior strengths in expected-case units.
#' @param window trailing weeks used for the estimate (the final week is excluded).
#' @return list with the previous contract (R_zone, logR_zone, sd_zone, seed_meanlog,
#'   seed_sdlog, R_nat, prov, ncase, r_min, r_max) plus the Gamma posterior `a_post`/`b_post`,
#'   the shrinkage diagnostics `prior_dominated` / `n_prior_dominated` / `n_with_cases`
#'   (the counts are among zones WITH cases — over all 519 they are dominated by the ~460
#'   never invaded and unreadable), `implausible` / `n_implausible`, `R_prov`, the local and
#'   imported denominators `loc` / `imp`, `import_share` and `window_weeks`.
#'   NOTE `seed_sdlog` and `sd_zone` are per-zone VECTORS; `seed_sdlog` was a scalar before.
estimate_zone_reff <- function(zone_week_nc, zones_all, gt_pmfs, gt = CASCADE_GT,
                               zone_province = NULL, window = 3L,
                               import_cases = NULL,
                               kappa_zone = get0("CASCADE_R_KAPPA_ZONE", ifnotfound = 5),
                               kappa_prov = get0("CASCADE_R_KAPPA_PROV", ifnotfound = 20),
                               kappa_nat  = get0("CASCADE_R_KAPPA_NAT",  ifnotfound = 1),
                               prior_mean = get0("CASCADE_R_PRIOR_MEAN", ifnotfound = 1),
                               r_min = CASCADE_R_MIN, r_max = CASCADE_R_MAX,
                               plausible_max = get0("CASCADE_R_PLAUSIBLE_MAX", ifnotfound = 4),
                               k_indiv = get0("CASCADE_K_INDIV", ifnotfound = 0.30),
                               em_iter = 200L, em_tol = 1e-8) {
  # em_iter is generous on purpose: one iteration is a handful of vectorised operations on a
  # 519-length vector, and the adversarial case (a very large import matrix) needed ~50+ to
  # reach 1e-8. There is no reason to trade accuracy for a cost this small.
  Y <- .count_wide(zone_week_nc, zones_all, "confirmed_nc")
  G <- cascade_weekly_gt(gt, gt_pmfs)
  nT <- ncol(Y)

  # The FINAL week is excluded, matching estimate_R_local() (15_workhorse.R): the last onset
  # week is still truncation-low even after nowcasting, and including it suppresses R.
  nT_use <- if (nT >= 3L) nT - 1L else nT
  # With a single usable week there is no renewal denominator to form, and the seq.int()
  # below would run BACKWARDS (seq.int(2, 1) is c(2, 1)) and index a column that does not
  # exist. Refuse explicitly rather than fail with a subscript error.
  if (nT_use < 2L)
    stop(sprintf(paste0("[reff] the record has %d usable week(s); the renewal estimator needs ",
                        "at least 2 (the final week is always excluded as truncation-low)."),
                 nT_use), call. = FALSE)
  wks <- seq.int(max(2L, nT_use - window + 1L), nT_use)

  nz  <- length(zones_all)
  num <- loc <- imp <- stats::setNames(numeric(nz), zones_all)
  for (w in wks) {
    num <- num + Y[, w]
    loc <- loc + .gweighted_own(Y, G, w)
    if (!is.null(import_cases) && w <= ncol(import_cases))
      imp <- imp + pmax(import_cases[, w], 0)
  }
  # Importation can only ever explain cases that were actually observed. Capping at the
  # observed count keeps the E-step a genuine split rather than letting a large modelled
  # force erase a zone's entire numerator.
  imp <- pmin(imp, num)

  # One EM run: attribute cases to the local component, then take the conjugate mean.
  # OVERDISPERSION. The simulator generates incidence as NegBin(mu = R*own, size = own*k),
  # for which Var/Mean = 1 + R/k EXACTLY (own cancels; verified against simulation). Fitting R
  # under a Poisson likelihood, as this estimator did, therefore assumed Var/Mean = 1 and
  # produced a posterior about sqrt(1 + R/k) too narrow — 2.1x at R = 1.08, k = 0.30. The
  # estimator contradicted the model's own stated assumption (Lloyd-Smith overdispersion,
  # k ~ 0.3), and every credible interval downstream inherited the over-confidence.
  #
  # The fix is quasi-likelihood precision scaling: a count contributes y/phi pseudo-events
  # against loc/phi pseudo-exposure, which leaves the posterior MEAN essentially unchanged
  # while inflating Var(log R) by phi. phi depends on R, so it is iterated to a fixed point
  # alongside the case-attribution EM rather than pinned at a starting value.
  .phi <- if (is.finite(k_indiv) && k_indiv > 0) function(R) pmax(1 + R / k_indiv, 1) else {
    warning("[reff] k_indiv is not a positive number; the likelihood falls back to POISSON ",
            "and the posterior will be too narrow. Set CASCADE_K_INDIV.", call. = FALSE)
    function(R) rep(1, length(R))
  }

  .em <- function(num, loc, imp, a0, b0, start) {
    R <- pmax(start, r_min); converged <- FALSE; step <- NA_real_; n_bad <- 0L
    for (it in seq_len(em_iter)) {
      lam_loc <- R * loc
      tot     <- lam_loc + imp
      y_loc   <- ifelse(tot > 0, num * lam_loc / tot, 0)
      ph      <- .phi(R)                      # overdispersion at the CURRENT iterate
      Rn      <- (a0 + y_loc / ph) / (b0 + loc / ph)
      # A non-finite iterate is a different failure from slow convergence and must not be
      # carried forward as if it were an estimate: hold those zones at their last finite
      # value and report how many there were.
      nf <- !is.finite(Rn)
      if (any(nf)) { n_bad <- max(n_bad, sum(nf)); Rn[nf] <- R[nf] }
      # The step must be measured BEFORE the update. Measuring it after assigning R <- Rn
      # makes it identically zero, so the warning reported "largest step 0.00e+00" while
      # claiming non-convergence.
      step <- max(abs(Rn - R))
      R <- Rn
      if (is.finite(step) && step < em_tol) { converged <- TRUE; break }
    }
    if (n_bad > 0L)
      warning(sprintf(paste0("[reff] the case-attribution EM produced a non-finite iterate for up ",
                             "to %d zone(s); those were held at their previous value. Check the ",
                             "import matrix for zones with no local history."), n_bad),
              call. = FALSE)
    if (!converged)
      warning(sprintf(paste0("[reff] the case-attribution EM did not converge in %d iterations ",
                             "(largest remaining step %.2e, tolerance %.1e); the estimate is the ",
                             "last iterate."), em_iter, step, em_tol), call. = FALSE)
    R
  }

  # ---- national: raw totals (mobility importation is internal at this scale) ----------
  num_nat <- sum(num); loc_nat <- sum(loc)
  # Same overdispersion scaling as the zone level, solved by fixed point because phi depends
  # on the very R being estimated. Two or three iterations suffice; the loop is bounded.
  R_nat <- prior_mean
  for (.it in seq_len(em_iter)) {
    .ph <- .phi(R_nat)
    .Rn <- as.numeric((kappa_nat * prior_mean + num_nat / .ph) / (kappa_nat + loc_nat / .ph))
    if (!is.finite(.Rn)) { .Rn <- R_nat; break }
    if (abs(.Rn - R_nat) < em_tol) { R_nat <- .Rn; break }
    R_nat <- .Rn
  }
  phi_nat <- .phi(R_nat)

  # ---- province: shrunk toward national, again on raw totals -------------------------
  prov <- if (is.null(zone_province)) rep("ALL", nz) else
            unname(zone_province[match(zones_all, names(zone_province))])
  prov[is.na(prov)] <- "ALL"
  pn <- tapply(num, prov, sum); pl <- tapply(loc, prov, sum)
  R_prov_by <- rep(R_nat, length(pn)); names(R_prov_by) <- names(pn)
  for (.it in seq_len(em_iter)) {
    .ph <- .phi(R_prov_by)
    .Rn <- (kappa_prov * R_nat + pn / .ph) / (kappa_prov + pl / .ph)
    .Rn[!is.finite(.Rn)] <- R_prov_by[!is.finite(.Rn)]
    if (max(abs(.Rn - R_prov_by)) < em_tol) { R_prov_by <- .Rn; break }
    R_prov_by <- .Rn
  }
  phi_prov_by <- .phi(R_prov_by)
  R_prov <- as.numeric(R_prov_by[prov]); names(R_prov) <- zones_all

  # ---- zone: shrunk toward its province pool, imports removed from the numerator -----
  a0 <- kappa_zone * R_prov; b0 <- rep(kappa_zone, nz)
  R_zone <- .em(num, loc, imp, a0, b0, start = R_prov)
  y_loc  <- ifelse(R_zone * loc + imp > 0, num * (R_zone * loc) / (R_zone * loc + imp), 0)
  phi_zone <- .phi(R_zone)
  a_post <- a0 + y_loc / phi_zone
  b_post <- b0 + loc   / phi_zone

  # Wide NUMERICAL guard. Not a prior: if it binds, the estimate is pathological and the
  # origin should be excluded rather than censored, so say so.
  bound <- R_zone < r_min | R_zone > r_max
  if (any(bound))
    warning(sprintf(paste0("[reff] %d zone(s) hit the numerical guard [%.2f, %.2f]; these are ",
                           "pathological estimates, not shrinkage — exclude the origin rather ",
                           "than reading them: %s"),
                    sum(bound), r_min, r_max,
                    paste(utils::head(zones_all[bound], 5), collapse = ", ")), call. = FALSE)
  R_zone <- pmin(pmax(R_zone, r_min), r_max)
  names(R_zone) <- zones_all

  # Posterior SD on the log scale, delta method: Var(log R) ~ 1/a_post. This REPLACES an
  # assumed constant, so a data-poor zone now carries wide uncertainty because it IS
  # data-poor rather than because a constant said so. Bounded only to keep the lognormal
  # draw numerically sane.
  # EXACT, not the delta-method approximation. For R ~ Gamma(a, b), Var(log R) = trigamma(a)
  # exactly, and it does not depend on the rate. The 1/sqrt(a) approximation understates that
  # by 12% at a = 2 and 4% at a = 6 — i.e. precisely for the data-poor zones whose
  # uncertainty this replacement of a fixed 0.30 existed to represent honestly. The floor is
  # a limit on CLAIMED precision, not a belief about spread; the ceiling keeps the lognormal
  # draw numerically sane.
  .sd_log_gamma <- function(a) sqrt(trigamma(pmax(a, 1e-6)))
  sd_zone <- pmin(pmax(.sd_log_gamma(a_post), 0.05), 1.5)
  names(sd_zone) <- zones_all

  # Effective, not raw, exposure: kappa_zone is in expected-case units, and overdispersion
  # means loc real cases carry only loc/phi cases' worth of information.
  prior_dominated <- (loc / phi_zone) < kappa_zone
  implausible     <- R_zone > plausible_max
  if (any(implausible))
    message(sprintf(paste0("[reff] %d zone(s) above the plausibility flag (R > %.1f): %s. ",
                           "Published EVD reproduction numbers sit near 1.5-2.5; early-outbreak ",
                           "values this high usually reflect ascertainment ramp-up, not transmission."),
                    sum(implausible), plausible_max,
                    paste(utils::head(zones_all[implausible], 5), collapse = ", ")))
  # Report prior-domination among zones that actually HAVE cases. Counting it over all 519
  # zones is dominated by the ~460 that have never been invaded and were never going to
  # carry local information, which makes the diagnostic unreadable.
  with_cases <- num > 0
  message(sprintf(paste0("[reff] R_nat=%.2f | weeks %s | zones with cases %d | of those, ",
                         "prior-dominated %d (%.0f%%) | importation = %.1f%% of the denominator"),
                  R_nat, paste(range(wks), collapse = "-"), sum(with_cases),
                  sum(prior_dominated & with_cases),
                  100 * sum(prior_dominated & with_cases) / max(sum(with_cases), 1L),
                  100 * sum(imp) / max(sum(loc) + sum(imp), 1e-9)))

  # A NEWLY seeded zone inherits its province pool (already shrunk to national), carrying
  # that pool's own posterior uncertainty rather than an assumed spread.
  a_prov <- as.numeric((kappa_prov * R_nat + (pn / phi_prov_by))[prov])
  seed_sdlog <- pmin(pmax(.sd_log_gamma(a_prov), 0.05), 1.5)

  list(R_zone = R_zone, logR_zone = log(R_zone), sd_zone = sd_zone,
       seed_meanlog = stats::setNames(log(pmax(R_prov, r_min)), zones_all),
       seed_sdlog = stats::setNames(seed_sdlog, zones_all),
       R_nat = R_nat, R_prov = R_prov, prov = prov, ncase = num,
       loc = loc, imp = imp, a_post = a_post, b_post = b_post,
       prior_dominated = prior_dominated,
       n_prior_dominated = sum(prior_dominated & num > 0), n_with_cases = sum(num > 0),
       implausible = implausible, n_implausible = sum(implausible),
       import_share = sum(imp) / max(sum(loc) + sum(imp), 1e-9),
       kappa_zone = kappa_zone, window_weeks = wks,
       r_min = r_min, r_max = r_max)
}

#' The cascade's NATIONAL reproduction number, from the shared EpiNow2 posterior.
#'
#' This is the single point at which transmission enters the 13-week projection. It calls the
#' same bayes_rt_week_draws() the short-term arm uses — same model, same generation time, same
#' as-of censoring rule, same cache — and asks for the average over the last
#' `window_weeks` weeks WITHIN each posterior sample (see that function for why a window mean
#' cannot be assembled from separate per-week fits).
#'
#' AS-OF SAFETY. `issue_date` is what makes this usable inside leave-future-out calibration:
#' bayes_rt_week_draws() censors the line list on it, so an origin fitted at its own cutoff
#' cannot see a case reported afterwards. The caller must pass it deliberately; there is no
#' default, because a default would silently be the production analysis date and every held-out
#' origin would train on the future.
#'
#' CLAMPING is numerical, not a prior: r_min/r_max exist so a pathological draw cannot make the
#' branching process run away. If it ever binds, that is reported rather than absorbed.
#'
#' @return list(draws, source, window_weeks, window_start, week_start, end_date, n_clamped).
cascade_rt_draws <- function(linelist, gt_pmfs, week_start, issue_date, gt = CASCADE_GT,
                             # FALLBACK 1L, matching RT_WINDOW_WEEKS. It was 3L, so any entry
                             # point that sourced this file WITHOUT 30_projection_config.R (which
                             # is what defines CASCADE_R_WINDOW_WEEKS, as RT_WINDOW_WEEKS) would
                             # silently average R over three weeks while the invasion arm used
                             # one — re-opening the two-arms-on-different-R split that
                             # RT_WINDOW_WEEKS exists to close, with nothing in the output to show it.
                             window_weeks = get0("CASCADE_R_WINDOW_WEEKS",
                                                 ifnotfound = get0("RT_WINDOW_WEEKS", ifnotfound = 1L)),
                             r_min = CASCADE_R_MIN, r_max = CASCADE_R_MAX) {
  if (is.null(linelist) || !is.data.frame(linelist))
    stop("[cascade-rt] the cascade's R comes from the line list via EpiNow2; none was supplied. ",
         "Pass linelist = layer$ll.", call. = FALSE)
  # NAME THE MISSING MODULE rather than letting R report "could not find function" from three
  # frames down. Taking R from EpiNow2 gave cascade_reff() a dependency it never had, and every
  # entry point that builds its own source list must now load 22_daily_reissue.R. When
  # 44_reff_epinow2_check.R did not, it failed four minutes into a run with an error that named
  # neither the cause nor the fix.
  if (!exists("linelist_observation_date", mode = "function"))
    stop("[cascade-rt] linelist_observation_date() is not loaded — source 22_daily_reissue.R. ",
         "cascade_reff() reaches it through cascade_rt_draws() -> bayes_rt_week_draws() to ",
         "censor the line list as of the issue date.", call. = FALSE)
  issue_date <- as.Date(issue_date)
  if (length(issue_date) != 1L || is.na(issue_date))
    stop("[cascade-rt] a single, non-missing issue_date is required: it is the as-of censoring ",
         "date, and getting it wrong is how a held-out origin sees the future.", call. = FALSE)
  pmf <- gt_pmfs[[gt]]
  if (is.null(pmf))
    stop("[cascade-rt] no generation-time pmf for profile '", gt, "'.", call. = FALSE)
  d <- bayes_rt_week_draws(linelist, pmf, week_start = as.Date(week_start),
                           issue_date = issue_date, window_weeks = window_weeks)
  src <- attr(d, "source"); wst <- attr(d, "window_start"); edt <- attr(d, "end_date")
  ww  <- attr(d, "window_weeks")
  v <- as.numeric(d); v <- v[is.finite(v) & v >= 0]
  if (!length(v))
    stop("[cascade-rt] EpiNow2 returned no usable R draws for the window ending ",
         format(week_start), ".", call. = FALSE)
  nclamp <- sum(v < r_min | v > r_max)
  if (nclamp > 0L)
    warning(sprintf(paste0("[cascade-rt] %d of %d R draws (%.1f%%) fell outside the numerical ",
                           "guard [%.2f, %.2f] and were clamped. This is a guard, not a prior — ",
                           "if it binds on a material share of draws the fit should be examined."),
                    nclamp, length(v), 100 * nclamp / length(v), r_min, r_max), call. = FALSE)
  v <- pmin(pmax(v, r_min), r_max)
  list(draws = v, source = src, window_weeks = ww, window_start = wst,
       week_start = as.Date(week_start), end_date = edt, n_clamped = nclamp)
}

# One MC draw of the initial R_eff vector for the currently-infected zones. Every zone in a
# replicate shares ONE value drawn from the national EpiNow2 posterior — see the body for why
# that sharing is the point rather than a simplification.
#
# `z` lets the caller SUPPLY the standard normals instead of drawing them. That is what
# makes a scenario contrast paired: the same z gives the same R whether or not the scenario
# seeded some other city, so the difference between two runs is the intervention rather
# than the random stream having moved on.
cascade_draw_R0 <- function(reff, z = NULL, gt_key = NULL) {
  nz <- length(reff$R_zone)
  # GT MARGINALISATION. `gt_key` selects the R posterior fitted AT that generation-time grid
  # point, so a parameter group using grid point i's weekly kernel also uses grid point i's R.
  # NULL keeps the single-GT behaviour. A key that is absent is an ERROR rather than a silent
  # fall-back to the central posterior: quietly pairing one GT's R with another's kernel is
  # exactly the incoherence this argument exists to prevent, and it would leave no trace.
  draws <- if (is.null(gt_key)) reff$R_draws else {
    byk <- reff$R_draws_by_gt
    if (is.null(byk) || is.null(byk[[gt_key]]))
      stop("[cascade] no R draws for generation-time grid point '", gt_key, "'. Build reff with ",
           "cascade_reff(gt_marginalise = TRUE) so every grid point carries its own posterior.",
           call. = FALSE)
    byk[[gt_key]]
  }
  if (is.null(draws) || !length(draws))
    stop("[cascade] this reff carries no R_draws; build it with cascade_reff(), which takes the ",
         "cascade's R from the shared EpiNow2 posterior.", call. = FALSE)
  # ONE DRAW, SHARED BY EVERY ZONE. This is the substantive consequence of anchoring the cascade
  # on a single national R: the posterior uncertainty in that number is COMMON to all zones, not
  # independent across them. Drawing per zone — as this function did when each zone had its own
  # Gamma posterior — would average the national uncertainty away across the ~60 infected zones
  # and make the projection interval far too narrow, which is the same error the nowcast
  # bootstrap comment in cascade_reff() warns about for the weekly correction factor. Perfectly
  # correlated R across zones is the honest representation and it WIDENS the national totals.
  #
  # Drawn by INVERSE CDF from the empirical posterior rather than from a fitted lognormal: the
  # EpiNow2 posterior is what we have, and resampling it keeps its actual shape and spread. It
  # also needs no moment-matching correction — the mean of the draws is the posterior mean by
  # construction, where the old lognormal parameterisation had to subtract sd^2/2 by hand.
  #
  # `z` keeps paired contrasts paired: only z[1] is used, because one number is being drawn.
  u <- if (is.null(z)) stats::runif(1L) else stats::pnorm(as.numeric(z)[1])
  # `draws`, NOT reff$R_draws. This line read reff$R_draws until 2026-09-16, which silently
  # discarded the gt_key selection above: every grid point returned the CENTRAL posterior's R
  # while the simulator used that grid point's own kernel — exactly the R/kernel mis-pairing
  # the gt_key argument exists to prevent, and invisible because the projection still moved
  # (the kernels varied). Any expression here must read the selected posterior.
  R <- as.numeric(stats::quantile(draws, probs = u, names = FALSE, type = 7))
  rep(pmin(pmax(R, reff$r_min), reff$r_max), nz)
}

# R for a freshly-seeded zone i, drawn at seeding time from the province R_eff pool.
#
# `seed_r` OVERRIDES that pool with an explicit reproduction number. A zone seeded during
# the projection inherits its province pool, and for a province with no case history that
# pool IS the national value — so the urban scenarios otherwise ask "what if a capital is
# seeded and then behaves like the current national average". That is a premise, not a
# finding, and CASCADE_SEED_R / CASCADE_SEED_R_SWEEP make it explicit. The pool's own
# posterior spread is retained around whatever mean is used.
#
# `z` supplies the standard normals for pairing, as in cascade_draw_R0(). Note seed_sdlog
# is a per-zone VECTOR (posterior SD of the province pool), so it is indexed alongside.
cascade_draw_R_seed <- function(reff, zone_idx, z = NULL, seed_r = NULL, R_grp = NULL,
                                gt_key = NULL) {
  n <- length(zone_idx)
  if (!n) return(numeric(0))
  # `R_grp` is this parameter group's national R — the value every already-infected zone in the
  # same replicate started from. A newly seeded zone gets THAT number, not an independent draw:
  # under one national R there is nothing left to distinguish it, and drawing again would
  # reintroduce exactly the spurious independence cascade_draw_R0() exists to avoid. It anchors
  # on the group's starting value rather than the walked one because the walk is a property of
  # the weeks a zone has already been transmitting through, which a fresh zone has not.
  # gt_key is threaded through so the FALLBACK draw comes from this group's generation-time grid
  # point too. Without it a seeded zone reached for the central posterior while every other zone
  # in the same replicate used its grid point's R.
  base <- if (!is.null(R_grp) && length(R_grp) == 1L && is.finite(R_grp)) as.numeric(R_grp)
          else cascade_draw_R0(reff, z = z, gt_key = gt_key)[1]
  if (!is.null(seed_r) && length(seed_r) == 1L && is.finite(seed_r) && seed_r > 0) {
    # The SWEPT premise: "this city transmits at R = seed_r". Read as a MEAN, as the sweep's
    # column headings claim, so the group's draw is rescaled to put E[R] exactly at seed_r while
    # keeping the posterior's relative spread: E[base / m] = 1, hence E[seed_r * base / m] =
    # seed_r. A bare seed_r would instead assert the city's R is known exactly, which is a
    # stronger claim than the sweep intends to make.
    #
    # The normaliser MUST come from the same posterior `base` was drawn from. Using the central
    # mean while base came from a grid point would leave E[R] off seed_r by the ratio between
    # the two — small, but wrong, and the same mis-pairing as the quantile bug above.
    .src <- if (!is.null(gt_key) && !is.null(reff$R_draws_by_gt) &&
                !is.null(reff$R_draws_by_gt[[gt_key]]))
      reff$R_draws_by_gt[[gt_key]] else reff$R_draws
    m <- mean(.src)
    base <- if (is.finite(m) && m > 0) seed_r * base / m else seed_r
  }
  rep(pmin(pmax(base, reff$r_min), reff$r_max), n)
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
#
# `u` supplies per-ZONE uniforms and the draw becomes an inverse CDF. Indexing by zone (not
# by position among the active zones) is what keeps the draw paired when an intervention
# changes WHICH zones are active. It also makes the coupling monotone: scaling `own` scales
# the mean and the size together, which is the sum-of-independent-draws representation, so
# the family is stochastically ordered and its quantiles are ordered pointwise. More
# infection upstream therefore can never produce fewer offspring downstream in a paired run.
cascade_branch_step <- function(R_vec, own_vec, k_indiv = CASCADE_K_INDIV, u = NULL,
                                mu_max = get0("CASCADE_MU_MAX", ifnotfound = 1e7)) {
  mu   <- R_vec * own_vec
  size <- own_vec * k_indiv
  out  <- integer(length(mu))
  act  <- own_vec > 0 & mu > 0
  if (any(act)) {
    # NUMERICAL BACKSTOP, not a modelling assumption and not currently live. The only brake
    # on growth in this layer is within-zone susceptible depletion, which needs BOTH
    # within-zone depletion, which this suite does not model (removed 2026-09-19); an R above 1
    # compounds unchecked for the whole horizon. Measured, that does not bite at the
    # reproduction numbers this analysis sweeps (up to 2.5 over 13 weeks) — no mean came
    # near mu_max — but qnbinom's search degrades badly for very large means, so the cap
    # exists to make a future parameterisation fail loudly instead of hanging. It sits far
    # above any DRC health zone's population (~2e6), so it cannot bind on an
    # epidemiologically sensible run.
    big <- act & mu > mu_max
    if (any(big)) {
      warning(sprintf(paste0("[branch] %d zone-week mean(s) exceeded mu_max = %.3g and were ",
                             "capped. The within-zone process has left the regime this model ",
                             "represents: enable within-zone depletion with a pop_vec, or lower ",
                             "the reproduction number. The projection is NOT usable as it ",
                             "stands."), sum(big), mu_max), call. = FALSE)
      mu[big] <- mu_max
    }
    out[act] <- if (is.null(u))
      stats::rnbinom(sum(act), mu = mu[act], size = pmax(size[act], 1e-6))
    else
      stats::qnbinom(u[act], mu = mu[act], size = pmax(size[act], 1e-6))
  }
  out
}

# ---------------------------------------------------------------------------
# The one way to build R_eff for the cascade
# ---------------------------------------------------------------------------

#' Zone R_eff with the import term, from a layer and its fitted hazard.
#'
#' Both the production run and every calibration origin must build R the SAME way, or the
#' fitted (delta, psi) describe a simulator that is not the one deployed. This wrapper is
#' the only supported path, and it now does TWO things:
#'
#'   1. It derives the per-zone import coefficient from the fit, scales it by the deployed
#'      delta, and hands that to estimate_zone_reff() — which is kept as a DIAGNOSTIC and for
#'      the case/importation bookkeeping (ncase, loc, imp, import_share, prov), not as the
#'      source of the reproduction number.
#'   2. It takes the reproduction number the cascade actually transmits on from the shared
#'      EpiNow2 posterior via cascade_rt_draws(): ONE national R, averaged over
#'      `rt_window_weeks` within each posterior sample, censored at `issue_date`.
#'
#' The returned object therefore carries BOTH national estimates — `R_nat` (EpiNow2, what the
#' simulator uses) and `R_nat_conjugate` (the conjugate estimator, for comparison) — so the
#' difference between them is reportable rather than one silently replacing the other.
#'
#' @param zone_week the count frame for THIS origin (as-of, not the final record).
#' @param delta deployed hazard recalibration factor for this kernel's covariate model.
#' @param use_imports FALSE reproduces the local-only estimator, for sensitivity.
#' @param n_seed_fun the seeding function the SIMULATION will use; passed through so that
#'   E[n_seed] in the import term can never silently disagree with what is simulated (it is
#'   derived from this function, not hard-coded).
#' @param nowcast_boot RETIRED and ignored; see the body. EpiNow2 carries its own
#'   right-truncation model, so resampling the nowcast factor would double-count it.
#' @param linelist the line list the national R is fitted from (defaults to layer$ll).
#' @param issue_date REQUIRED, no default. The as-of date the R fit censors on. Production
#'   passes ANALYSIS_DATE; a calibration origin passes its OWN cutoff + 6 (the last day of its
#'   last training week — the convention the LFO, the R(t) issue date and the deployed weekly
#'   anchor all share). There is no default
#'   because the only sensible one would be the production date, and an origin that inherited
#'   it would train on data from after its own cutoff — a leak the outputs would not reveal.
#' @param rt_window_weeks weeks of R to average over, ending at the frame's last week.
cascade_reff <- function(zone_week, zones_all, layer, fit, design, delta = 1,
                         zone_province = NULL, gt = CASCADE_GT,
                         kernel = CASCADE_KERNEL,
                         use_imports = get0("CASCADE_R_USE_IMPORTS", ifnotfound = TRUE),
                         import_beta = c("intercept", "design", "covariate"),
                         n_seed_fun = cascade_draw_seed,
                         nowcast_boot = get0("CASCADE_NOWCAST_BOOT", ifnotfound = 100L),
                         # The line list and the as-of date the NATIONAL R is fitted on. There is
                         # deliberately NO default for issue_date: it would be the production
                         # analysis date, and every calibration origin that forgot to override it
                         # would quietly fit R on data from after its own cutoff. Missing it is an
                         # error the caller sees, not a leak the reader never finds.
                         linelist = layer$ll,
                         issue_date,
                         # see cascade_rt_draws(): the fallback must be the invasion arm's window.
                         rt_window_weeks = get0("CASCADE_R_WINDOW_WEEKS",
                                                ifnotfound = get0("RT_WINDOW_WEEKS", ifnotfound = 1L)),
                         # Fetch an R posterior for EVERY GT_PRIOR grid point, so the simulator
                         # can marginalise over generation time. Off by default; see
                         # CASCADE_GT_MARGINALISE in 30_projection_config.R.
                         gt_marginalise = get0("CASCADE_GT_MARGINALISE", ifnotfound = FALSE),
                         gt_prior = get0("GT_PRIOR", ifnotfound = NULL),
                         ...) {
  import_beta <- match.arg(import_beta)
  # RNG HYGIENE, at the TOP of the function and matching simulate_cascade(). Registered here
  # rather than beside the bootstrap because anything consuming randomness BEFORE the guard
  # escapes restoration — which is exactly why the first placement failed its own test. Seed
  # from a fixed derived constant so the bootstrap is reproducible for identical inputs, and
  # restore the caller's stream on exit so nothing downstream is perturbed (run_cascade calls
  # this before calibration and before the scenario simulations).
  .old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(.old_seed)) assign(".Random.seed", .old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(get0("CASCADE_SEED", ifnotfound = 20260704L) + 7919L)   # distinct from the sim seed
  imp <- NULL
  if (isTRUE(use_imports)) {
    # WHICH CONVERSION, and why the default changed. It was design$beta0, the intercept-only
    # cloglog glm fitted alongside the design. That glm UNDERFLOWS on short training frames:
    # on the 2026-09-01 record it returned exactly 0 at all five calibration origins while
    # giving 0.138 on the full frame, so every calibration origin silently fell back to the
    # local-only estimator and (delta, psi) were fitted against an R built differently from
    # the one deployed. The default is now the fitted hazard's own intercept, exp(E[b0]) —
    # always available, regularised by its prior, and the conversion the SIMULATOR itself
    # applies to an average zone, which is the consistency that matters here.
    #   "design"    the legacy glm scalar (kept for reproducing earlier outputs)
    #   "covariate" the per-zone version; a documented sensitivity, NOT a default — see
    #               cascade_import_beta() for why its negative population coefficient is a
    #               link-saturation artefact.
    beta_i <- switch(import_beta,
      intercept = tryCatch(cascade_import_beta(fit, design, layer$covariates, zones_all,
                             covariates_on = FALSE), error = function(e) NULL),
      covariate = tryCatch(cascade_import_beta(fit, design, layer$covariates, zones_all),
                           error = function(e) NULL),
      design    = design$beta0)
    if (is.null(beta_i) || !any(is.finite(beta_i)) || max(beta_i, na.rm = TRUE) <= 0) {
      warning(sprintf(paste0("[reff] the import conversion ('%s') is unusable (%s); falling back ",
                             "to the fitted hazard intercept."), import_beta,
                      if (is.null(beta_i)) "NULL" else sprintf("max = %.3g", max(beta_i, na.rm = TRUE))),
              call. = FALSE)
      beta_i <- tryCatch(cascade_import_beta(fit, design, layer$covariates, zones_all,
                           covariates_on = FALSE), error = function(e) NULL)
    }
    imp <- tryCatch(cascade_import_cases(zone_week, zones_all, layer$mobility_matrices,
                      layer$gt_pmfs, kernel = kernel, gt = gt, beta0 = beta_i, delta = delta,
                      n_seed_fun = n_seed_fun),
                    error = function(e) {
                      warning("[reff] the import term could not be built (", conditionMessage(e),
                              "); falling back to the local-only estimator.", call. = FALSE)
                      NULL })
    # A SILENT degradation to local-only is the failure this wrapper exists to prevent: it
    # puts the calibration origins and the production run on different estimators, and
    # nothing in the output would say so. Imports were asked for; if they are not there, say
    # it loudly enough to stop a publication run.
    if (is.null(imp))
      warning("[reff] importation was requested but the import matrix could not be built, so ",
              "this R_eff is LOCAL-ONLY. If other calls in the same run did build it, their ",
              "estimates are not comparable with this one — do not calibrate across the two.",
              call. = FALSE)
  }
  # The conjugate estimator is RETAINED, but demoted to a diagnostic. It still supplies the
  # case and importation bookkeeping the report and the hazard scale reference (ncase, loc,
  # imp, import_share, prov) and a second opinion on the national level. It no longer supplies
  # the R the cascade transmits on; see 30_projection_config.R for the three measured reasons.
  conj <- estimate_zone_reff(zone_week, zones_all, layer$gt_pmfs, gt = gt,
                             zone_province = zone_province, import_cases = imp, ...)

  # ---- the anchor the projection actually runs on -------------------------------------
  # ONE national R from the shared EpiNow2 posterior, averaged over rt_window_weeks WITHIN each
  # sample, censored at `issue_date` so a held-out calibration origin cannot see its own future.
  # week_start is taken from the frame itself, so an as-of frame automatically anchors on its
  # own last week rather than on production's.
  wk_last <- max(as.Date(zone_week$week_start))
  rt <- cascade_rt_draws(linelist, layer$gt_pmfs, week_start = wk_last,
                         issue_date = issue_date, gt = gt,
                         window_weeks = rt_window_weeks,
                         r_min = conj$r_min, r_max = conj$r_max)
  Rd     <- rt$draws
  nz     <- length(zones_all)
  R_nat  <- mean(Rd)
  # Reported for the CSV and the narrative only. The simulator never uses it: it resamples
  # Rd directly (cascade_draw_R0), so the posterior's real shape is preserved rather than
  # being forced through a lognormal.
  sd_log <- stats::sd(log(pmax(Rd, 1e-12)))
  if (!is.finite(sd_log)) sd_log <- 0

  # ---- NOWCAST UNCERTAINTY is deliberately NOT propagated here any more ----------------
  # The retired block resampled `confirmed_nc` by its fitted completeness factor and re-ran the
  # conjugate estimator B times, widening sd_zone. That has no place in the new path. EpiNow2 is
  # fitted to RAW daily onset counts and carries its own right-truncation model (trunc_opts,
  # built from the same fitted reporting delay the nowcast uses), so the incompleteness of the
  # recent tail is already inside the posterior. Resampling the nowcast factor on top would
  # correct for the same thing twice — exactly the double-counting that feeding EpiNow2 raw
  # rather than nowcast-corrected counts exists to avoid.
  .nb <- suppressWarnings(as.integer(nowcast_boot))
  if (length(.nb) == 1L && !is.na(.nb) && .nb > 1L)
    message("[reff] nowcast_boot (", .nb, ") is not used: R now comes from EpiNow2, which ",
            "models the incomplete recent tail itself, so resampling the nowcast factor ",
            "would double-count it.")

  out <- conj
  # Captured BEFORE the national and per-zone fields are overwritten below, so the conjugate
  # picture stays inspectable rather than one estimate silently replacing the other. These are
  # DIAGNOSTICS with no path into the simulation: nothing downstream transmits on them, and
  # cascade_reff_by_zone.csv labels every one of them `diag_`. Without this the per-zone CSV
  # would fall back to an all-NA column, which reads as "not estimated" rather than "estimated
  # and not used".
  out$R_nat_conjugate   <- conj$R_nat
  out$R_zone_conjugate  <- conj$R_zone
  out$sd_zone_conjugate <- conj$sd_zone
  out$R_prov_conjugate  <- conj$R_prov
  out$R_draws   <- Rd
  out$R_nat     <- R_nat
  out$R_zone    <- stats::setNames(rep(R_nat, nz), zones_all)
  out$logR_zone <- log(out$R_zone)
  out$sd_zone   <- stats::setNames(rep(sd_log, nz), zones_all)
  out$R_prov    <- stats::setNames(rep(R_nat, nz), zones_all)
  # A newly seeded zone is drawn from the same national posterior as everyone else, so the
  # "seed pool" is that posterior. Kept under the old names so consumers do not have to branch.
  out$seed_meanlog <- stats::setNames(rep(log(max(R_nat, 1e-12)), nz), zones_all)
  out$seed_sdlog   <- stats::setNames(rep(sd_log, nz), zones_all)
  out$rt_source       <- rt$source
  out$rt_window_weeks <- rt$window_weeks
  out$rt_window_start <- rt$window_start
  out$rt_window_end   <- rt$end_date
  out$rt_n_draws      <- length(Rd)
  out$rt_n_clamped    <- rt$n_clamped
  out$rt_issue_date   <- as.Date(issue_date)
  out$rt_quantiles    <- stats::quantile(Rd, c(0.05, 0.25, 0.5, 0.75, 0.95), names = TRUE)

  # ---- GENERATION-TIME MARGINALISATION: one R posterior per GT_PRIOR grid point --------
  # R is estimated AT a generation time, so marginalising the projection over GT_PRIOR needs a
  # posterior for each grid point; the simulator then pairs grid point i's R with grid point i's
  # weekly kernel. This is what makes the cascade's treatment of GT uncertainty the same as the
  # short-term arm's, which stopped selecting a GT and marginalises over this same grid.
  #
  # Every fit is cached, so this loop is cheap on a warm cache and expensive exactly once.
  # A grid point whose fit fails is DROPPED and its prior weight redistributed over the rest,
  # with a warning — rather than substituting the central posterior, which would silently
  # mis-pair an R with a kernel it was not estimated under.
  if (isTRUE(gt_marginalise)) {
    if (is.null(gt_prior))
      stop("[reff] gt_marginalise = TRUE needs a GT prior; GT_PRIOR is not available.",
           call. = FALSE)
    gp   <- make_gt_prior_pmfs(gt_prior)
    keys <- names(gp$gt_pmfs)
    byk  <- stats::setNames(vector("list", length(keys)), keys)
    for (k in keys) {
      rk <- tryCatch(cascade_rt_draws(linelist, gp$gt_pmfs, week_start = wk_last,
                                      issue_date = issue_date, gt = k,
                                      window_weeks = rt_window_weeks,
                                      r_min = conj$r_min, r_max = conj$r_max),
                     error = function(e) { warning(sprintf(
                       "[reff] GT grid point '%s' (%.2f/%.2f d) has no usable R and is dropped: %s",
                       k, attr(gp$gt_pmfs[[k]], "gt_mean"), attr(gp$gt_pmfs[[k]], "gt_sd"),
                       conditionMessage(e)), call. = FALSE); NULL })
      if (!is.null(rk)) byk[[k]] <- rk$draws
    }
    keep <- !vapply(byk, is.null, logical(1))
    if (!any(keep))
      stop("[reff] gt_marginalise = TRUE but no GT grid point yielded an R posterior.",
           call. = FALSE)
    if (!all(keep))
      warning(sprintf(paste0("[reff] %d of %d GT grid points were dropped; the remaining prior ",
                             "weight (%.3f of 1) is renormalised over the rest, which narrows the ",
                             "GT range actually marginalised over."),
                      sum(!keep), length(keep), sum(gp$weights[keep])), call. = FALSE)
    out$R_draws_by_gt <- byk[keep]
    out$gt_keys       <- keys[keep]
    out$gt_weights    <- gp$weights[keep] / sum(gp$weights[keep])
    out$gt_grid       <- gp$grid[keep, , drop = FALSE]
    # The prior-weighted mean across grid points: the quantity a GT-marginalised projection is
    # centred on, reported beside the single-GT anchor so the shift from marginalising is visible.
    out$R_nat_gt_marginal <- sum(vapply(out$R_draws_by_gt, mean, numeric(1)) * out$gt_weights)
    message(sprintf(paste0("[reff] GT marginalisation ON: %d grid points, prior-weighted R = ",
                           "%.3f (single-GT anchor %.3f, %.1f%% apart); R by GT spans %.3f-%.3f"),
                    length(out$gt_keys), out$R_nat_gt_marginal, R_nat,
                    100 * abs(out$R_nat_gt_marginal / max(R_nat, 1e-12) - 1),
                    min(vapply(out$R_draws_by_gt, mean, numeric(1))),
                    max(vapply(out$R_draws_by_gt, mean, numeric(1)))))
  }
  .rr <- R_nat / max(conj$R_nat, 1e-12)
  message(sprintf(paste0("[reff] cascade anchor: EpiNow2 R = %.3f [%.3f, %.3f] (90%%), %d draws ",
                         "over %s..%s (source %s, issue %s) | conjugate estimator gives %.3f ",
                         "on the same frame (EpiNow2 is %.0f%% %s)"),
                  R_nat, out$rt_quantiles[1], out$rt_quantiles[5], length(Rd),
                  format(as.Date(rt$window_start)), format(as.Date(rt$end_date)),
                  rt$source, format(as.Date(issue_date)), conj$R_nat,
                  100 * abs(.rr - 1), if (.rr >= 1) "higher" else "lower"))
  out
}
