# =============================================================================
# 21c_frailty_models.R — the TIME-VARYING FORMS of the invasion model (Stan engine)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Stage 2 of the featured-model selection (00_config.R, TV_FORMS; 16_invasion_eval.R,
# select_featured_two_stage()) refits the stage-1 specification (kernel, covariates, generation
# time) under three forms of the import coefficient, all fitted by ONE Stan model
# (stan/invasion_frailty_rw.stan) with two switches:
#   use_rw       RW1 drift u_t in the log import coefficient (a COMMON drift: response,
#                reporting, behaviour), persisting into the forecast with one innovation a week;
#   use_frailty  a zone gamma frailty Z ~ Gamma(1/theta, 1/theta), integrated out. With this
#                week's hazard h = exp(b0 + x'gamma + u) * Lam and the zone's own cumulative
#                hazard so far C = sum over past weeks of h (each week at ITS covariates and
#                drift), a zone still uninvaded has
#                    P(invaded this week) = 1 - (1 + theta * h / (1 + theta * C))^(-1/theta).
#                The factor 1/(1 + theta*C) falls with exposure, so it CHANGES THE RANKING.
# CAVEAT, stated: the frailty cannot tell a zone that resists invasion from one that was invaded
# but not detected; it down-weights exposed, silent zones either way.
#
# WHY STAN, NOT BRMS: past exposure must be weighted by the drift level and covariates in force
# when it occurred, which a brms non-linear formula cannot express. With both switches off the
# Stan model is the brms fixed-beta model (same likelihood, same priors); the 2026-10-14 check
# found forecast correlation 0.9999 and totals within 0.3% between the two fits.
#
# Everything here follows the brms models' conventions so the forms are scored, deployed and
# reported exactly like every other model: the same rolling as-of design (build_invasion_design),
# the same forecast offsets and EpiNow2 source projection (bayes_forecast_offsets,
# .bayes_project_mu), the same draws-list format, recalibration (.scale_invasion_draws) and
# summaries (.finalise_invasion_summary), and the same parameter / diagnostic / posterior-draw /
# beta-trajectory schemas as bayes_param_table(), bayes_fit_diagnostics(),
# bayes_posterior_draws() and bayes_beta_trajectory().
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))
suppressPackageStartupMessages({ library(dplyr) })

FRAILTY_STAN_FILE <- file.path(here::here(), "spatiotemporal", "stan", "invasion_frailty_rw.stan")

#' Mean and SD of a "normal(m, s)" prior string, the form every prior in 00_config.R takes.
#' Anything else stops: the Stan model implements Normal priors only, and a prior it cannot
#' represent must not be replaced by a different one silently.
.normal_prior_params <- function(spec, what) {
  m <- regmatches(spec, regexec(
    "^\\s*normal\\s*\\(\\s*([^,()]+?)\\s*,\\s*([^,()]+?)\\s*\\)\\s*$", spec))[[1]]
  v <- if (length(m) == 3L) suppressWarnings(as.numeric(m[2:3])) else c(NA_real_, NA_real_)
  if (anyNA(v) || v[2] <= 0)
    stop(sprintf("[stan-form] %s = '%s' is not a normal(mean, sd) prior with sd > 0.", what, spec),
         call. = FALSE)
  v
}

#' The effect_scale of an RW1 form's intercept row (parameter table and posterior draws): exp(b0)
#' is the walk's level BEFORE the first training week, not the import coefficient in force.
RW1_START_LEVEL_SCALE <- "import coefficient before the first training week (RW1 start level)"

#' The import force AND the standardised covariates of EVERY zone at EVERY training transition.
#'
#' Transition s is week s -> s+1. L[z, s] is the import force into zone z, computed from the
#' counts KNOWN at the end of week s, and XA[[s]] the zone covariates for that transition, both
#' exactly as build_invasion_design() builds its row (zone, .week = s + 1): the same rolling as-of
#' reconstruction (.asof_counts_wide), compute_foi() and .feature_matrix() with the d_min frontier
#' as known that week, standardised with the design's own centre and scale. They are needed for
#' every zone, not only the design rows, because a zone's cumulative exposure includes weeks
#' before it entered the at-risk design. stan_invasion_fit() checks the two agree on every design
#' row.
#' @return list(L = zones x S matrix, XA = list of S zones x K matrices).
invasion_exposure <- function(zw, linelist, zones_all, mobility_matrices, gt_pmfs, covariates,
                              osrm_mat, mob, gt, feat, center, scale, trunc_delay) {
  W <- mobility_matrices[[mob]]; G <- daily_to_weekly_gt(gt_pmfs[[gt]])
  Y_wide <- .count_wide(zw, zones_all, "confirmed_nc")
  A_wide <- .count_wide(zw, zones_all, "total_alerts")
  S_wide <- .susp_wide(zw, zones_all)
  static <- .static_features(covariates, zones_all)
  wk <- suppressWarnings(as.Date(colnames(Y_wide))); S <- ncol(Y_wide) - 1L
  if (S < 1L || anyNA(wk)) stop("[stan-form] the training counts have no dated transitions.", call. = FALSE)
  L <- matrix(0, length(zones_all), S, dimnames = list(zones_all, NULL))
  XA <- vector("list", S)
  for (s in seq_len(S)) {
    XA[[s]] <- matrix(0, length(zones_all), length(feat), dimnames = list(zones_all, feat))
    Yl <- .asof_counts_wide(linelist, zones_all, wk[s] + 6L, trunc_delay, week_spine = wk[seq_len(s)])
    # The design skips a transition whose as-of reconstruction fails; it then carries no
    # exposure either (zero), and no row of the design refers to it.
    if (is.null(Yl) || ncol(Yl) < 1L) next
    v <- compute_foi(Yl, W, G, t_idx = ncol(Yl) + 1L, zones_all); v[!is.finite(v)] <- 0
    L[, s] <- as.numeric(v)
    if (length(feat)) {
      X <- .feature_matrix(s + 1L, Y_wide, A_wide, W, G, static, osrm_mat, zones_all, feat,
                           Y_active = Yl, S_wide = S_wide)
      if (is.null(X) || !all(feat %in% colnames(X)))
        stop(sprintf("[stan-form] covariate(s) %s unavailable at transition %d.",
                     paste(setdiff(feat, colnames(X)), collapse = ", "), s), call. = FALSE)
      X <- X[, feat, drop = FALSE]
      for (f in feat) X[, f] <- (X[, f] - center[[f]]) / scale[[f]]
      if (any(!is.finite(X)))
        stop(sprintf("[stan-form] non-finite covariate at transition %d; past exposure cannot be weighted.", s),
             call. = FALSE)
      XA[[s]][] <- X
    }
  }
  list(L = L, XA = XA)
}

#' Cumulative hazards for the 1- and 2-week windows under the (frailty) model.
#'
#' All arguments are draws x zones matrices (or recycle to them). h1, h2 = this and next week's
#' hazard exp(b0 + x'gamma + u) * Lam; C = the zone's cumulative hazard over all past weeks;
#' theta = frailty variance per draw (NULL = no frailty). Without frailty H1 = h1 and
#' H2 = h1 + h2 -- exactly the pipeline's cumulative hazard. With frailty the 2-week hazard is the
#' joint survival over both weeks, not a sum of weekly hazards. H2 >= H1 always, so the windows
#' stay nested.
#' @return list(H1, H2)
frailty_cum_hazard <- function(h1, h2, C, theta = NULL) {
  if (is.null(theta)) return(list(H1 = h1, H2 = h1 + h2))
  den <- 1 + theta * C
  list(H1 = log1p(theta * h1 / den) / theta,
       H2 = log1p(theta * (h1 + h2) / den) / theta)
}

#' The model's log-likelihood in R: a reference implementation. stan_invasion_loglik() is tested
#' against it; Stan's own density is not (that needs a compiled model).
#' @param par list(b0, beta, lth = NULL (no frailty), u = NULL (no drift))
#' @param X design-row covariates (N x K); XA list of S zones x K matrices (all zones).
frailty_loglik_R <- function(par, y, X, lam, zone, tr, L, XA) {
  S <- ncol(L); u <- if (is.null(par$u)) rep(0, S) else par$u
  xb <- function(M) if (length(par$beta)) drop(M %*% par$beta) else 0
  h <- exp(par$b0 + u[tr] + xb(X)) * lam
  H <- if (is.null(par$lth)) h else {
    th <- exp(par$lth)
    E <- matrix(vapply(seq_len(S), function(s) exp(par$b0 + u[s] + xb(XA[[s]])) * L[, s],
                       numeric(nrow(L))), nrow(L))
    Cz <- t(apply(cbind(0, E), 1, cumsum))[, seq_len(S), drop = FALSE]
    log1p(th * h / (1 + th * Cz[cbind(zone, tr)])) / th
  }
  sum(ifelse(y == 1, log(-expm1(-H)), -H))
}

.FRAILTY_MODEL <- new.env(parent = emptyenv())
#' The compiled Stan model, cached per R process. run_all.R compiles it once in the main process
#' before any parallel back-test, so the workers find an up-to-date executable and never compile
#' concurrently into the same directory.
frailty_stan_model <- function(stan_file = FRAILTY_STAN_FILE,
                               bin_dir = file.path(get0("OUT_DIAGNOSTICS", ifnotfound = tempdir()),
                                                   "stan_bin")) {
  if (!is.null(.FRAILTY_MODEL$m)) return(.FRAILTY_MODEL$m)
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("[stan-form] cmdstanr is required for the time-varying forms.", call. = FALSE)
  dir.create(bin_dir, recursive = TRUE, showWarnings = FALSE)
  .FRAILTY_MODEL$m <- cmdstanr::cmdstan_model(stan_file, dir = bin_dir, quiet = TRUE)
  .FRAILTY_MODEL$m
}

#' Fit ONE time-varying form on an invasion design.
#'
#' @param des build_invasion_design() output (rolling as-of design).
#' @param ex invasion_exposure() for the same training counts and specification.
#' @param cov_spec the covariates of the specification (the fit uses those present in des$feat).
#' @param use_rw,use_frailty the form's switches (TV_FORMS, 00_config.R).
#' @param iter total iterations per chain (half warm-up), as for the brms fitters.
#' @return an object of class "stan_invasion_fit": list(draws = draws x variables matrix, feat,
#'   center, scale, S, L, XA, zones, use_rw, use_frailty, n_obs, n_events, diag, label,
#'   priors), or NULL when sampling fails at both settings (the failure is logged).
stan_invasion_fit <- function(des, ex, zones_all, cov_spec, use_rw, use_frailty, iter,
                              label = NA_character_,
                              seed = get0("RANDOM_SEED", ifnotfound = 20260704L)) {
  d <- des$d; L <- ex$L; S <- ncol(L); feat <- intersect(cov_spec, des$feat)
  zi <- match(d$.zone, zones_all); tr <- as.integer(d$.week - 1L); lam <- exp(d$logLam)
  if (anyNA(zi) || any(tr < 1L | tr > S))
    stop("[stan-form] design rows fall outside the exposure matrix (zones or transitions).", call. = FALSE)
  X <- if (length(feat)) as.matrix(d[, feat, drop = FALSE]) else matrix(0, nrow(d), 0)
  if (any(!is.finite(X)))
    stop("[stan-form] the design has non-finite covariates; brms would drop those rows and the two ",
         "engines would no longer fit the same data.", call. = FALSE)
  # The design and the exposure are built by the same as-of reconstruction; any difference means
  # they disagree about the data, so the fit is refused rather than run on inconsistent inputs.
  chk <- max(abs(lam / pmax(L[cbind(zi, tr)], 1e-300) - 1))
  if (!is.finite(chk) || chk > 1e-6)
    stop(sprintf("[stan-form] design import force differs from the exposure matrix (max rel diff %.3g)",
                 chk), call. = FALSE)
  if (length(feat)) {
    xa_rows <- matrix(vapply(seq_len(nrow(d)), function(n) ex$XA[[tr[n]]][zi[n], feat],
                             numeric(length(feat))), ncol = length(feat), byrow = TRUE)
    dx <- max(abs(xa_rows - X))
    if (!is.finite(dx) || dx > 1e-8)
      stop(sprintf("[stan-form] design covariates differ from the all-zone covariates (max abs diff %.3g)",
                   dx), call. = FALSE)
  }
  pri <- list(
    b0_prior    = .normal_prior_params(BAYES_PRIOR_INTERCEPT, "BAYES_PRIOR_INTERCEPT"),
    beta_prior  = .normal_prior_params(BAYES_PRIOR_COEF, "BAYES_PRIOR_COEF"),
    sigma_prior = .normal_prior_params(TV_FORM_PRIOR_RW_SD, "TV_FORM_PRIOR_RW_SD"),
    lth_prior   = .normal_prior_params(TV_FORM_PRIOR_LOG_THETA, "TV_FORM_PRIOR_LOG_THETA"))
  # XA as an S x Z x K array (Stan's array[S] matrix[Z, K]).
  xa <- array(0, dim = c(S, length(zones_all), length(feat)))
  if (length(feat)) for (s in seq_len(S)) xa[s, , ] <- ex$XA[[s]][, feat, drop = FALSE]
  sdat <- c(list(N = nrow(d), y = as.integer(d$invaded), K = ncol(X), X = X, lam = lam,
                 Z = length(zones_all), S = S, zone = zi, tr = tr, L = unname(L), XA = xa,
                 use_rw = as.integer(isTRUE(use_rw)), use_frailty = as.integer(isTRUE(use_frailty))),
            pri)
  # Same bounded intercept init as the brms fitters (21_bayesian_renewal.R), and the brms
  # time-varying fitter's retry settings (adapt_delta 0.95, then 0.99 from the PRIOR MEAN of the
  # intercept, BAYES_PRIOR_INTERCEPT), so an init-caused failure does not simply recur.
  lb0 <- suppressWarnings(log(des$beta0)); if (!is.finite(lb0)) lb0 <- pri$b0_prior[1]
  lb0 <- min(max(lb0, -20), -max(d$logLam))
  # Every parameter gets an init, the switched-off ones as zero-length arrays (cmdstanr otherwise
  # reports them as missing).
  init_at <- function(v) function() list(b0 = v, beta = as.array(rep(0, ncol(X))),
                                         lth = if (isTRUE(use_frailty)) as.array(0) else array(numeric(0)),
                                         sigma = if (isTRUE(use_rw)) as.array(0.3) else array(numeric(0)),
                                         zr = if (isTRUE(use_rw)) rep(0, S) else numeric(0))
  fitf <- function(ad, v) frailty_stan_model()$sample(
    data = sdat, chains = 2, parallel_chains = 2, iter_warmup = as.integer(iter) %/% 2L,
    iter_sampling = as.integer(iter) %/% 2L, seed = as.integer(seed), init = init_at(v),
    adapt_delta = ad, max_treedepth = 12, refresh = 0, show_messages = FALSE,
    show_exceptions = FALSE)
  .ok <- function(f) !is.null(f) && isTRUE(tryCatch(posterior::ndraws(f$draws()) > 0,
                                                    error = function(e) FALSE))
  fit <- tryCatch(suppressWarnings(fitf(0.95, lb0)), error = function(e) NULL)
  if (!.ok(fit)) fit <- tryCatch(suppressWarnings(fitf(0.99, pri$b0_prior[1])), error = function(e) NULL)
  if (!.ok(fit)) {
    .bayes_log_fit_failure(des, feat, "cloglog",
                           sprintf("[stan-form %s] sampling failed at adapt_delta 0.95 and 0.99", label))
    return(NULL)
  }
  dr0 <- fit$draws(format = "draws_matrix")
  dr  <- matrix(as.numeric(dr0), nrow(dr0), dimnames = list(NULL, colnames(dr0)))
  # CONVERGENCE over every SAMPLED parameter: b0, beta, log theta, sigma and the walk's
  # standardised innovations zr (lp__ and the derived u are not parameters). This is stricter
  # than bayes_fit_diagnostics(), which reads the population-level and SD terms of a brms fit.
  pars <- c("b0", if (ncol(X)) paste0("beta[", seq_len(ncol(X)), "]"),
            if (isTRUE(use_frailty)) "lth[1]",
            if (isTRUE(use_rw)) c("sigma[1]", paste0("zr[", seq_len(S), "]")))
  sm <- posterior::summarise_draws(posterior::subset_draws(fit$draws(), variable = pars),
                                   "rhat", "ess_bulk", "ess_tail")
  nd_div <- tryCatch(as.integer(sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)),
                     error = function(e) NA_integer_)
  fin <- function(x) if (is.finite(x)) x else NA_real_
  diag <- data.frame(model = label, n_divergent = nd_div,
                     pct_divergent = if (is.finite(nd_div)) 100 * nd_div / nrow(dr) else NA_real_,
                     rhat_max = fin(suppressWarnings(max(sm$rhat, na.rm = TRUE))),
                     ess_bulk_min = fin(suppressWarnings(min(sm$ess_bulk, na.rm = TRUE))),
                     ess_tail_min = fin(suppressWarnings(min(sm$ess_tail, na.rm = TRUE))),
                     n_draws = nrow(dr), stringsAsFactors = FALSE)
  structure(list(draws = dr, feat = feat, center = des$center[feat], scale = des$scale[feat],
                 S = S, L = L, XA = lapply(ex$XA, function(m) m[, feat, drop = FALSE]),
                 zones = zones_all, use_rw = isTRUE(use_rw), use_frailty = isTRUE(use_frailty),
                 n_obs = nrow(d), n_events = sum(d$invaded), diag = diag, label = label,
                 priors = pri),
            class = "stan_invasion_fit")
}

# ---- draw helpers ----------------------------------------------------------------------------
.sf_beta <- function(fit) {
  if (length(fit$feat)) fit$draws[, paste0("beta[", seq_along(fit$feat), "]"), drop = FALSE] else NULL
}
.sf_U <- function(fit) fit$draws[, paste0("u[", seq_len(fit$S), "]"), drop = FALSE]
# draws x zones linear predictor b0 + x'beta for a zones x K covariate matrix (no drift).
.sf_xb <- function(fit, Xz) {
  ND <- nrow(fit$draws); Z <- length(fit$zones)
  out <- matrix(fit$draws[, "b0"], ND, Z)
  B <- .sf_beta(fit)
  if (!is.null(B)) {
    if (is.null(Xz) || nrow(Xz) != Z || ncol(Xz) != ncol(B))
      stop("[stan-form] covariate matrix does not match the fit (zones x covariates).", call. = FALSE)
    out <- out + B %*% t(Xz)
  }
  out
}
# draws x zones cumulative hazard over the first `upto` transitions (sum_{s <= upto} h_s).
.sf_cum_exposure <- function(fit, upto = fit$S) {
  ND <- nrow(fit$draws); U <- .sf_U(fit); Z <- length(fit$zones)
  C <- matrix(0, ND, Z)
  for (s in seq_len(upto))
    C <- C + exp(.sf_xb(fit, fit$XA[[s]]) + U[, s]) * matrix(fit$L[, s], ND, Z, byrow = TRUE)
  C
}
# Forecast-week drift (draws x n_weeks): RW1 persists u_S and adds one innovation per week, the
# innovations drawn reproducibly from RANDOM_SEED with the caller's RNG stream restored;
# otherwise 0. The forecast and the beta trajectory call this with the same seed, so the
# trajectory's forecast weeks are the forecast's own drift.
.sf_future_u <- function(fit, n_weeks) {
  ND <- nrow(fit$draws)
  if (!fit$use_rw) return(matrix(0, ND, n_weeks))
  .had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  .old <- if (.had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({ if (.had) assign(".Random.seed", .old, envir = globalenv())
            else suppressWarnings(rm(".Random.seed", envir = globalenv())) }, add = TRUE)
  set_fixed_seed(as.integer(get0("RANDOM_SEED", ifnotfound = 20260704L)) + 7411L)
  e <- matrix(stats::rnorm(ND * n_weeks), ND, n_weeks)
  cs <- if (n_weeks == 1L) e else t(apply(e, 1, cumsum))
  .sf_U(fit)[, fit$S] + fit$draws[, "sigma[1]"] * cs
}

#' Pointwise log-likelihood (draws x rows) of a stan_invasion_fit on design rows. Not used for
#' stacking (the forms are not stacked members, fit_bayes_suite()); kept as the tested reference
#' for frailty_loglik_R(). `rows` are rows of the SAME design the model was fitted on
#' (standardised covariates, .zone, .week, logLam, invaded).
stan_invasion_loglik <- function(fit, rows) {
  ND <- nrow(fit$draws); Z <- length(fit$zones)
  zi <- match(rows$.zone, fit$zones); tr <- as.integer(rows$.week - 1L)
  if (anyNA(zi) || any(tr < 1L | tr > fit$S)) stop("[stan-form] log-lik rows outside the fit.", call. = FALSE)
  B <- .sf_beta(fit); U <- .sf_U(fit)
  lin <- matrix(fit$draws[, "b0"], ND, nrow(rows)) + U[, tr, drop = FALSE]
  if (!is.null(B)) lin <- lin + B %*% t(as.matrix(rows[, fit$feat, drop = FALSE]))
  h <- exp(sweep(lin, 2, rows$logLam, "+"))
  H <- if (!fit$use_frailty) h else {
    th <- exp(fit$draws[, "lth[1]"])
    Cn <- matrix(0, ND, nrow(rows)); Crun <- matrix(0, ND, Z)
    for (s in seq_len(fit$S)) {
      ix <- which(tr == s)
      if (length(ix)) Cn[, ix] <- Crun[, zi[ix], drop = FALSE]
      Crun <- Crun + exp(.sf_xb(fit, fit$XA[[s]]) + U[, s]) * matrix(fit$L[, s], ND, Z, byrow = TRUE)
    }
    log1p(th * h / (1 + th * Cn)) / th
  }
  y <- matrix(as.integer(rows$invaded), ND, nrow(rows), byrow = TRUE)
  ifelse(y == 1L, log(-expm1(-H)), -H)
}

#' Posterior draws of the cumulative hazard for the forecast weeks of a stan_invasion_fit.
#'
#' Week 1's import force is fixed by observed counts (offsets row h = 1). Week 2's comes from the
#' pipeline's per-draw source projection (.bayes_project_mu, EpiNow2 R draws), into which week-1
#' import-seeded cases enter as new sources: their expected number in zone z is h1 * E[Z_z | its
#' history], which is h1 / (1 + theta * C) for a zone with no case yet and h1 (the prior mean
#' of the frailty) for a zone that already has cases, whose frailty the model does not track
#' after invasion. Without frailty both are h1, exactly the brms projection. Covariates are those
#' of each forecast week's own rows; the RW1 walk persists its last training level and adds one
#' innovation per forecast week.
#' @return a draws list in the format of .bayes_invasion_draws() (list per horizon of
#'   list(zones, horizon, cum = draws x zones, p)), RAW (no recalibration).
stan_invasion_forecast_draws <- function(fit, off, horizons) {
  stopifnot(inherits(fit, "stan_invasion_fit"))
  hz <- sort(unique(as.integer(horizons)))
  if (length(setdiff(hz, c(1L, 2L))))
    stop("[stan-form] only horizons 1 and 2 are supported.", call. = FALSE)
  ND <- nrow(fit$draws); zones_all <- fit$zones; Z <- length(zones_all); feat <- fit$feat
  for (f in feat) {
    v <- (off[[f]] - fit$center[[f]]) / fit$scale[[f]]; v[!is.finite(v)] <- 0; off[[f]] <- v
  }
  rows_h <- function(h) {
    i <- which(off$horizon == h)
    if (!identical(as.character(off$health_zone[i]), zones_all))
      stop(sprintf("[stan-form] the horizon-%d forecast rows are not the zone spine in order.", h),
           call. = FALSE)
    i
  }
  Xh <- function(i) if (length(feat)) as.matrix(off[i, feat, drop = FALSE]) else NULL
  uf <- .sf_future_u(fit, max(hz))
  i1 <- rows_h(1L)
  b1 <- exp(.sf_xb(fit, Xh(i1)) + uf[, 1])
  h1 <- b1 * matrix(exp(off$logLam[i1]), ND, Z, byrow = TRUE)
  th <- if (fit$use_frailty) exp(fit$draws[, "lth[1]"]) else NULL
  Cf <- if (fit$use_frailty) .sf_cum_exposure(fit) else matrix(0, ND, Z)
  h2 <- matrix(0, ND, Z)
  if (max(hz) >= 2L) {
    i2 <- rows_h(2L)
    proj <- attr(off, "projection")
    affected <- rowSums(proj$Y[zones_all, , drop = FALSE] > 0) > 0
    zcond <- if (fit$use_frailty) 1 / (1 + th * Cf) else matrix(1, ND, Z)
    zcond[, affected] <- 1
    b2 <- exp(.sf_xb(fit, Xh(i2)) + uf[, 2])
    eta0 <- matrix(0, ND, nrow(off))
    eta0[, i1] <- log(b1 * zcond)
    eta0[, i2] <- log(b2)
    nd <- off; nd$.off <- nd$logLam
    mu <- .bayes_project_mu(eta0, nd, proj, "cloglog", 2L)
    h2 <- mu[, i2, drop = FALSE]                      # = b2 * projected week-2 import force
  }
  H <- frailty_cum_hazard(h1, h2, Cf, theta = th)
  lapply(stats::setNames(hz, paste0("h", hz)), function(h) {
    cum <- matrix(if (h == 1L) H$H1 else H$H2, ND, Z, dimnames = list(NULL, zones_all))
    list(zones = zones_all, horizon = h, cum = cum, p = 1 - exp(-cum))
  })
}

#' Design, exposure and fit of one form on a zone-week frame (the rolling as-of design).
#' @return list(des, fit, status): status "ok", "no_design" (the design could not be built),
#'   "no_events" (no invasion to fit -- the same case the brms closure returns NULL for) or
#'   "fit_failed" (sampling failed at both settings); des/fit are NULL unless status is "ok".
.stan_form_design_fit <- function(spec, zw, mobility_matrices, gt_pmfs, covariates, osrm_mat,
                                  zones_all, linelist, trunc_delay, iter, label) {
  des <- tryCatch(build_invasion_design(zw, mobility_matrices, gt_pmfs, covariates, osrm_mat,
           zones_all, mob = spec$mob, gt = spec$gt,
           candidates = unique(c(spec$cov, "log_pop", "ccvi", "d_min", "healthsite_density")),
           rolling_ll = linelist, rolling_delay = trunc_delay), error = function(e) NULL)
  if (is.null(des)) return(list(des = NULL, fit = NULL, status = "no_design"))
  if (des$n_events < 1L) return(list(des = NULL, fit = NULL, status = "no_events"))
  feat <- intersect(spec$cov, des$feat)
  ex <- invasion_exposure(zw, linelist, zones_all, mobility_matrices, gt_pmfs, covariates, osrm_mat,
                          spec$mob, spec$gt, feat, des$center, des$scale, trunc_delay)
  fit <- stan_invasion_fit(des, ex, zones_all, spec$cov, spec$use_rw, spec$use_frailty,
                           iter = iter, label = label)
  if (is.null(fit)) return(list(des = NULL, fit = NULL, status = "fit_failed"))
  list(des = des, fit = fit, status = "ok")
}

#' The forms need the ROLLING as-of design: their exposure IS the as-of import force. Refuse to
#' run otherwise rather than fit them on a different design from the fixed-beta model they are
#' compared with.
.stan_form_require_rolling <- function(linelist, trunc_delay) {
  if (!isTRUE(get0("ROLLING_PREDICTORS", ifnotfound = TRUE)) || is.null(linelist) || is.null(trunc_delay))
    stop("[stan-form] the time-varying forms need the rolling as-of design (ROLLING_PREDICTORS on, ",
         "a line list and the as-of truncation); set INCLUDE_TV_FORMS = FALSE otherwise.", call. = FALSE)
}

#' LFO closure for a time-varying form, with the signature run_invasion_lfo() expects.
#'
#' Mirrors make_bayes_lfo_model(): the same rolling as-of design, the same EpiNow2 R draws for the
#' 2-week source projection, the same forecast offsets and summary. Only the hazard model differs.
#' @param spec a form spec (engine "stan"; bayes_tv_form_specs(), 21_bayesian_renewal.R).
#' @param diag_dir optional directory; each fold writes ONE row to its own file
#'   "<label>__<cutoff>.csv" there (the back-test runs a fold's models in parallel workers, so a
#'   shared file appended to concurrently could interleave rows): the convergence diagnostics and
#'   a `status` -- "ok", "no_design", "no_events", "fit_failed", "offsets_failed" (fitted, but the
#'   forecast offsets failed, so the fold has no forecast) or "error: <message>" (the row is
#'   written, then the error is re-raised). The caller gathers them
#'   (collect_stan_form_diagnostics()).
make_stan_form_lfo_model <- function(spec, mobility_matrices, gt_pmfs, covariates, osrm_mat,
                                     zones_all, linelist, trunc_delay, iter, diag_dir = NULL) {
  force(spec); force(linelist); force(trunc_delay); force(iter); force(diag_dir)
  .stan_form_require_rolling(linelist, trunc_delay)
  function(zw, ti, hz, cut) {
    rt <- if (max(hz) >= 2L)
            bayes_rt_week_draws(linelist, gt_pmfs[[spec$gt]], week_start = max(zw$week_start),
                                issue_date = as.Date(cut) + 6L) else NULL
    .write_diag <- function(df, status) {
      if (is.null(diag_dir)) return(invisible(NULL))
      fitted <- !is.null(df) && !is.null(df$fit)
      row <- if (!fitted) data.frame(model = spec$label, n_divergent = NA_integer_,
                                     pct_divergent = NA_real_, rhat_max = NA_real_,
                                     ess_bulk_min = NA_real_, ess_tail_min = NA_real_,
                                     n_draws = NA_integer_)
             else df$fit$diag
      row$cutoff <- as.character(cut); row$ok <- fitted; row$status <- status
      row$theta_median <- if (fitted && df$fit$use_frailty)
        stats::median(exp(df$fit$draws[, "lth[1]"])) else NA_real_
      row$sigma_median <- if (fitted && df$fit$use_rw)
        stats::median(df$fit$draws[, "sigma[1]"]) else NA_real_
      dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(row, file.path(diag_dir, sprintf("%s__%s.csv", spec$label, as.character(cut))),
                       row.names = FALSE)
    }
    df <- tryCatch(.stan_form_design_fit(spec, zw, mobility_matrices, gt_pmfs, covariates, osrm_mat,
                                         zones_all, linelist, trunc_delay, iter, spec$label),
                   error = function(e) {
                     .write_diag(NULL, paste0("error: ", conditionMessage(e)))
                     stop(e)
                   })
    if (!identical(df$status, "ok")) { .write_diag(df, df$status); return(NULL) }
    off <- tryCatch(bayes_forecast_offsets(zw, mobility_matrices, gt_pmfs, covariates, osrm_mat,
             zones_all, spec$mob, spec$gt, spec$cov, hz, rt_draws = rt),
             error = function(e) NULL)
    .write_diag(df, if (is.null(off)) "offsets_failed" else "ok")
    if (is.null(off)) return(NULL)
    .finalise_invasion_summary(stan_invasion_forecast_draws(df$fit, off, hz),
                               affected_zones = character(0))
  }
}

#' Gather the per-fit diagnostics written by make_stan_form_lfo_model() into one table (and
#' remove the per-fit files). NULL when there are none.
collect_stan_form_diagnostics <- function(diag_dir) {
  fs <- if (dir.exists(diag_dir)) list.files(diag_dir, pattern = "__.*\\.csv$", full.names = TRUE)
        else character(0)
  if (!length(fs)) return(NULL)
  out <- dplyr::bind_rows(lapply(fs, function(f) utils::read.csv(f, stringsAsFactors = FALSE)))
  unlink(fs)
  out[order(out$model, out$cutoff), , drop = FALSE]
}

#' Posterior parameter summary in the schema of bayes_param_table(): the Intercept row is
#' exp(b0) and covariate rows per-SD hazard ratios exp(gamma), as for the brms fits. The forms
#' add `frailty_theta` (the frailty variance theta; 0 = no frailty) and `rw1_sigma_week`
#' (exp(sigma): the multiplicative spread of beta from a one-SD weekly innovation; 1 = no drift).
#' Under RW1 the intercept is the level BEFORE the first transition (u starts there; its
#' effect_scale is RW1_START_LEVEL_SCALE and beta_last_training_week carries exp(b0 + u_S)), not the
#' current import coefficient; the trajectory (stan_invasion_beta_trajectory()) gives that.
stan_invasion_param_table <- function(fit, label = "") {
  if (is.null(fit)) return(NULL)
  dr <- fit$draws; rmx <- fit$diag$rhat_max
  row <- function(term, x, is_int, scale, pdir)
    tibble::tibble(model = label, term = term, is_intercept = is_int, effect_scale = scale,
                   hr = stats::median(x), lo = stats::quantile(x, 0.05, names = FALSE),
                   hi = stats::quantile(x, 0.95, names = FALSE),
                   p_dir = if (pdir) mean(x > 1) else NA_real_, rhat = rmx)
  # Under RW1 exp(b0) is NOT the import coefficient in force: it is the walk's starting level,
  # so its row says so (effect_scale), and the level at the last training week, exp(b0 + u_S) --
  # the value bayes_beta_over_folds.csv reports as beta0 -- is a row of its own.
  out <- list(row("Intercept", exp(dr[, "b0"]), TRUE,
                  if (fit$use_rw) RW1_START_LEVEL_SCALE else "hazard ratio", FALSE))
  if (fit$use_rw)
    out[[length(out) + 1L]] <- row("beta_last_training_week",
                                   exp(dr[, "b0"] + .sf_U(fit)[, fit$S]), FALSE,
                                   "import coefficient at the last training week (RW1 level)", FALSE)
  for (k in seq_along(fit$feat))
    out[[length(out) + 1L]] <- row(fit$feat[k], exp(dr[, sprintf("beta[%d]", k)]), FALSE,
                                   "hazard ratio", TRUE)
  if (fit$use_frailty)
    out[[length(out) + 1L]] <- row("frailty_theta", exp(dr[, "lth[1]"]), FALSE,
                                   "frailty variance theta (0 = no frailty)", FALSE)
  if (fit$use_rw)
    out[[length(out) + 1L]] <- row("rw1_sigma_week", exp(dr[, "sigma[1]"]), FALSE,
                                   "exp(RW1 innovation SD) per week (1 = no drift)", FALSE)
  dplyr::bind_rows(out)
}

#' Convergence diagnostics in the schema of bayes_fit_diagnostics().
stan_invasion_diagnostics <- function(fit, label = "") {
  if (is.null(fit)) return(NULL)
  d <- fit$diag; d$model <- label
  tibble::as_tibble(d[, c("model", "n_divergent", "pct_divergent", "rhat_max", "ess_bulk_min",
                          "ess_tail_min", "n_draws")])
}

#' Posterior draws in the schema of bayes_posterior_draws(), each row stamped with the prior it
#' was fitted under (on the log scale). The Stan model puts BAYES_PRIOR_INTERCEPT on b0 itself,
#' so the intercept's stamp is always valid (bayes_posterior_draws() has to check this for brms,
#' which places it on the centred intercept). The process parameters (theta, sigma) are not
#' hazard ratios and are reported in the parameter table, not here.
stan_invasion_posterior_draws <- function(fit, label = "") {
  if (is.null(fit)) return(NULL)
  dr <- fit$draws
  one <- function(term, x, is_int, prior, scale = "hazard ratio")
    tibble::tibble(model = label, term = term, is_intercept = is_int, effect_scale = scale,
                   prior = prior, .draw = seq_along(x), hr = exp(x))
  # exp(b0) is the parameter BAYES_PRIOR_INTERCEPT describes; under RW1 it is the walk's
  # starting level, and the row says so.
  out <- list(one("Intercept", dr[, "b0"], TRUE, BAYES_PRIOR_INTERCEPT,
                  if (fit$use_rw) RW1_START_LEVEL_SCALE else "hazard ratio"))
  for (k in seq_along(fit$feat))
    out[[length(out) + 1L]] <- one(fit$feat[k], dr[, sprintf("beta[%d]", k)], FALSE, BAYES_PRIOR_COEF)
  # The form's PROCESS parameters, on the scales of stan_invasion_param_table(), so the
  # prior-against-posterior figure can show the parameter that defines the form: theta = exp(lth)
  # with lth ~ TV_FORM_PRIOR_LOG_THETA, and exp(sigma) with sigma ~ TV_FORM_PRIOR_RW_SD restricted
  # to sigma > 0 (Stan's <lower=0>), stamped "half-" so the figure draws the truncated prior.
  if (fit$use_frailty)
    out[[length(out) + 1L]] <- one("frailty_theta", dr[, "lth[1]"], FALSE, TV_FORM_PRIOR_LOG_THETA,
                                   "frailty variance theta (0 = no frailty)")
  if (fit$use_rw)
    out[[length(out) + 1L]] <- one("rw1_sigma_week", dr[, "sigma[1]"], FALSE,
                                   paste0("half-", trimws(TV_FORM_PRIOR_RW_SD)),
                                   "exp(RW1 innovation SD) per week (1 = no drift)")
  dplyr::bind_rows(out)
}

#' Posterior trajectory of beta_t in the schema of bayes_beta_trajectory(): exp(b0 + u_t) with
#' covariates at their standardised mean (0), for every training week and the forecast weeks
#' (RW1: the walk propagated one innovation per week, the same draws the forecast uses;
#' otherwise constant). Week index w is the design's `.week` (transition w - 1). The frailty
#' does not enter beta_t: it acts on each zone through its own exposure.
stan_invasion_beta_trajectory <- function(fit, design, label = "", horizons = NULL, week_dates = NULL) {
  if (is.null(fit) || is.null(design) || !".week" %in% names(design$d)) return(NULL)
  dr <- fit$draws; U <- .sf_U(fit)
  weeks <- sort(unique(as.integer(design$d$.week)))
  b <- exp(dr[, "b0"] + U[, weeks - 1L, drop = FALSE])
  fut <- if (length(horizons)) max(weeks) + sort(unique(as.integer(horizons))) else integer(0)
  if (length(fut)) {
    # Forecast week max(weeks) + k is transition S + k only when the last design week is the last
    # training week (S + 1); otherwise the forecast weeks are not those of .sf_future_u().
    if (max(weeks) != fit$S + 1L)
      stop("[stan-form] the design's last week is not the last training week.", call. = FALSE)
    uf <- .sf_future_u(fit, max(fut) - max(weeks))
    b <- cbind(b, exp(dr[, "b0"] + uf[, fut - max(weeks), drop = FALSE]))
  }
  all_w <- c(weeks, fut)
  out <- tibble::tibble(model = label, week = all_w, is_forecast = !(all_w %in% weeks),
                        beta = apply(b, 2, stats::median),
                        beta_lo = apply(b, 2, stats::quantile, 0.05, names = FALSE),
                        beta_hi = apply(b, 2, stats::quantile, 0.95, names = FALSE),
                        tv = if (fit$use_rw && fit$use_frailty) "rw1_frailty"
                             else if (fit$use_rw) "rw1" else if (fit$use_frailty) "frailty" else "none")
  out$week_date <- if (!is.null(week_dates)) unname(as.Date(week_dates[as.character(out$week)]))
                   else as.Date(NA)
  out
}

message("[21c] time-varying forms loaded (Stan: ", basename(FRAILTY_STAN_FILE), ")")
