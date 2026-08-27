# =============================================================================
# 21_bayesian_renewal.R — Bayesian mobility-informed renewal INVASION models
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# A Bayesian (brms / Stan) re-expression of the frequentist mobility-informed
# renewal invasion model, EXTENDING (not replacing) the suite. The invasion
# hazard of an at-risk zone i is
#     P(first case)_i = 1 - exp(-beta_i * Lambda_i),   beta_i = exp(eta_i),
# with the mobility import force Lambda_i entering as a FIXED OFFSET (log Lambda)
# and eta_i = intercept (+ optional covariates). This is EXACTLY a Bernoulli GLM
# with a complementary-log-log link and an offset, so it is fit directly with
# brms as
#     brm(invaded ~ <covariates> + offset(log Lambda),
#         family = bernoulli(link = "cloglog"), prior = <weakly-informative>).
#
# Why Bayesian here (beyond taste):
#   * proper posterior uncertainty on every predicted invasion probability and on
#     every parameter (credible intervals), rather than a delta-method / bootstrap
#     bolt-on;
#   * the priors REGULARISE the quasi-separation that made the frequentist JOINT
#     MLE of the covariate model diverge (coefficients -> +/-1e15): the posterior
#     stays finite and honest, with wide intervals where the data are silent;
#   * a principled combination ACROSS structural assumptions (mobility kernel /
#     covariates) via loo predictive stacking (Yao et al. 2018 Bayesian Analysis
#     13:917-1007, doi:10.1214/17-BA1091).
#
# Priors are WEAKLY-INFORMATIVE (documented below): not tight enough to drive the
# posterior, not so flat that separation reappears.
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))
suppressPackageStartupMessages({
  library(tidyverse); library(brms); library(posterior); library(loo)
})
if (!exists("%||%"))
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 ||
                              (length(a) == 1 && is.na(a))) b else a

# --- Priors (see METHODS §Bayesian) ------------------------------------------
# Intercept = log(beta0), the log baseline import->first-case coefficient.
# Empirically beta0 ~ 0.05-0.10, so normal(-3, 2) gives a 90% prior interval on
# beta0 of roughly [0.002, 1.3]: wide, gently favouring the small import->seeding
# conversion expected for a ~0.3-0.6% base-rate event, NOT a tight constraint.
# Coefficients act on STANDARDISED covariates (per-SD log hazard ratios).
# normal(0, 1) is the Gelman-style weakly-informative default: a 1-SD move in a
# covariate multiplies the hazard by up to ~e^2 within two prior SDs, shrinks
# toward "no effect", and regularises separation without imposing a direction.
BAYES_PRIOR_INTERCEPT <- "normal(-3, 2)"
BAYES_PRIOR_COEF      <- "normal(0, 1)"

.bayes_prior <- function(feat) {
  pr <- brms::prior_string(BAYES_PRIOR_INTERCEPT, class = "Intercept")
  if (length(feat)) pr <- pr + brms::prior_string(BAYES_PRIOR_COEF, class = "b")
  pr
}

#' Record a Bayesian fit failure with a design fingerprint so any residual failure is
#' diagnosable rather than silently dropped. Appends one line to a diagnostics log.
.bayes_log_fit_failure <- function(design, feat, link, msg) {
  d <- design$d
  line <- sprintf("mob=%s gt=%s link=%s feat=%s n_obs=%d n_events=%s beta0=%s logLam=[%s,%s] :: %s",
    design$mob %||% "?", design$gt %||% "?", link,
    if (length(feat)) paste(feat, collapse = "+") else "none",
    nrow(d), sum(d$invaded), signif(design$beta0, 3),
    signif(min(d$logLam, na.rm = TRUE), 3), signif(max(d$logLam, na.rm = TRUE), 3), msg)
  dir <- if (exists("OUT_DIAGNOSTICS")) OUT_DIAGNOSTICS else tempdir()
  try(cat(line, "\n", file = file.path(dir, "bayes_fit_failures.log"), append = TRUE), silent = TRUE)
  warning("[bayes] fit failed: ", msg)
}

#' Fit ONE Bayesian cloglog renewal invasion model on an at-risk design.
#' @param design list from build_invasion_design() (d has invaded, logLam,
#'   standardised covariates; plus center/scale/beta0).
#' @param cov_spec covariates to include (subset of design$feat).
#' @param link observation-process link for the binary invasion outcome. cloglog
#'   (DEFAULT, principled) makes the mobility import force a proper log-cumulative-
#'   hazard offset and yields exactly p = 1 - exp(-beta*Lambda) — the renewal
#'   invasion probability. logit / probit are supported as robustness sensitivities
#'   (their offset is a valid GLM shift but loses the clean hazard interpretation);
#'   prediction reconciles all three on the per-week hazard scale.
fit_bayes_renewal <- function(design, cov_spec = character(0),
                              iter = 2000L, chains = 2L, seed = get0("RANDOM_SEED", ifnotfound = 20260704L),
                              adapt_delta = 0.9, link = "cloglog") {
  d <- design$d; if (is.null(d) || !nrow(d)) return(NULL)
  feat <- intersect(cov_spec, design$feat)
  rhs  <- if (length(feat)) paste(feat, collapse = " + ") else "1"
  fml  <- stats::as.formula(sprintf("invaded ~ %s + offset(logLam)", rhs))
  # Deterministic initialisation at the EMPIRICAL null-model coefficient. brms/cmdstanr
  # otherwise inits the intercept randomly in [-2,2] on the link scale (NOT from the
  # normal(-3,2) prior); with a large fixed log(Lambda) offset (e.g. the inward-focused
  # M11 mobility) a non-invaded zone then gets eta = intercept + logLam large positive,
  # so under cloglog p = 1-exp(-exp(eta)) -> 1 and the Bernoulli(0) log-density is -Inf.
  # Every random start is rejected and cmdstanr returns no draws ("Unable to retrieve the
  # metadata"). Starting the intercept at log(beta0) (the null-model MLE, which already
  # absorbs the offset) and slopes at 0 makes the initial log-density finite by construction.
  lb0 <- suppressWarnings(log(design$beta0))
  if (!is.finite(lb0)) lb0 <- -3
  init_list <- lapply(seq_len(chains), function(i)
    if (length(feat)) list(Intercept = lb0, b = rep(0, length(feat))) else list(Intercept = lb0))
  .brm <- function(ad) suppressMessages(suppressWarnings(brms::brm(
    fml, data = d, family = brms::bernoulli(link = link),
    prior = .bayes_prior(feat), iter = iter, warmup = iter %/% 2L,
    chains = chains, cores = min(chains, 4L), seed = seed, init = init_list,
    backend = "cmdstanr", refresh = 0, silent = 2,
    control = list(adapt_delta = ad, max_treedepth = 12L))))
  fit <- tryCatch(.brm(adapt_delta), error = function(e) {
    # One retry at a higher adapt_delta before giving up — covers transient divergence
    # or a marginal init on a sparse fold. The real error (not the swallowed metadata
    # message) is recorded so any residual failure is diagnosable, not silent.
    tryCatch(.brm(max(adapt_delta, 0.97)), error = function(e2) {
      .bayes_log_fit_failure(design, feat, link, conditionMessage(e2)); NULL }) })
  if (!is.null(fit)) { attr(fit, "feat") <- feat; attr(fit, "cov_spec") <- cov_spec
                       attr(fit, "link") <- link }
  fit
}

#' Time-VARYING-beta Bayesian renewal fit (review §4.2 — reviewer point E1).
#'
#' Generalises fit_bayes_renewal by letting the import->invasion conversion scale
#' evolve over weeks. A smooth weekly intercept u(week) is added on the cloglog
#' (log-hazard) scale, so beta_t = exp(beta0 + u(week_t)) and
#'   cloglog(P(invade_{i,t})) = beta0 + u(week_t) + <covariates> + log Lambda_{i,t}.
#' u(week) is a Gaussian process on the integer week index (brms gp(.week)) — the
#' native, forecast-friendly analogue of a random walk on the intercept:
#'   * STRICT GENERALISATION: its marginal SD sdgp -> 0 recovers the fixed-beta
#'     model exactly, so the data (sdgp posterior) decide HOW MUCH time-variation
#'     is warranted — replacing the manuscript's "time-variation is unhelpful"
#'     assertion with an estimate.
#'   * FORECAST BEHAVIOUR: the GP extrapolates to the forecast week (see
#'     predict_bayes_invasion_tv) with GROWING uncertainty, the coverage-increasing
#'     behaviour a time-varying term should provide out-of-sample.
#'   * SMALL-SAMPLE HONESTY: with ~39 events sdgp is weakly identified; a
#'     data-driven sdgp ~ 0 ("time-variation not warranted") is a legitimate,
#'     informative result — report the sdgp posterior, do not force time-variation.
#' Falls back to the fixed-beta fit when there is no .week column or < 3 training
#' weeks. sdgp carries a half-normal (normal(0,1) on the positive-constrained
#' scale) weakly-informative prior; the GP length-scale uses the brms default.
#'
#' @inheritParams fit_bayes_renewal
#' @param k_gp GP basis dimension (default 5; small — matches the handful of
#'   training weeks, keeps the fit identifiable and cheap).
#' @return a brmsfit with attr("tv")=TRUE and attr("max_train_week"), or NULL.
fit_bayes_renewal_tv <- function(design, cov_spec = character(0),
                                 iter = 2000L, chains = 2L,
                                 seed = get0("RANDOM_SEED", ifnotfound = 20260704L),
                                 adapt_delta = 0.95, link = "cloglog", k_gp = 5L) {
  d <- design$d; if (is.null(d) || !nrow(d)) return(NULL)
  # Guard: the time-varying term needs a weekly index and >= 3 distinct weeks.
  if (!".week" %in% names(d) || length(unique(d$.week)) < 3L) {
    warning("[tv] design lacks .week or has < 3 training weeks; using fixed-beta fit")
    return(fit_bayes_renewal(design, cov_spec, iter, chains, seed, adapt_delta, link))
  }
  nwk  <- length(unique(d$.week))
  feat <- intersect(cov_spec, design$feat)
  rhs  <- if (length(feat)) paste(feat, collapse = " + ") else NULL
  kk   <- as.integer(min(k_gp, max(3L, nwk)))               # cap basis at #weeks
  terms <- c("1", sprintf("gp(.week, k = %d, c = 5/4)", kk), rhs)
  fml   <- stats::as.formula(sprintf("invaded ~ %s + offset(logLam)",
                                     paste(terms[!is.na(terms)], collapse = " + ")))
  tv_prior <- .bayes_prior(feat) +
    brms::prior_string("normal(0, 1)", class = "sdgp")       # half-normal on the time-variation scale
  lb0 <- suppressWarnings(log(design$beta0)); if (!is.finite(lb0)) lb0 <- -3
  init_list <- lapply(seq_len(chains), function(i)
    if (length(feat)) list(Intercept = lb0, b = rep(0, length(feat))) else list(Intercept = lb0))
  .brm <- function(ad) suppressMessages(suppressWarnings(brms::brm(
    fml, data = d, family = brms::bernoulli(link = link),
    prior = tv_prior, iter = iter, warmup = iter %/% 2L,
    chains = chains, cores = min(chains, 4L), seed = seed, init = init_list,
    backend = "cmdstanr", refresh = 0, silent = 2,
    control = list(adapt_delta = ad, max_treedepth = 12L))))
  fit <- tryCatch(.brm(adapt_delta), error = function(e)
    tryCatch(.brm(max(adapt_delta, 0.99)), error = function(e2) {
      .bayes_log_fit_failure(design, feat, link, paste0("[tv] ", conditionMessage(e2))); NULL }))
  if (!is.null(fit)) { attr(fit, "feat") <- feat; attr(fit, "cov_spec") <- cov_spec
                       attr(fit, "link") <- link; attr(fit, "tv") <- TRUE
                       attr(fit, "max_train_week") <- max(d$.week, na.rm = TRUE) }
  fit
}

#' Prediction for a time-varying-beta fit (review §4.2). Evaluates the GP weekly
#' intercept at the FORECAST week (max training week + horizon) so the
#' time-varying component extrapolates forward with growing uncertainty. Adds the
#' required .week column to the forecast offsets, then delegates to the shared
#' predict_bayes_invasion (identical summarisation incl. rank CrI).
predict_bayes_invasion_tv <- function(fit, offsets_df, design, horizons,
                                      affected_zones = character(0),
                                      rho = ASCERTAINMENT_NOMINAL) {
  if (is.null(fit)) return(NULL)
  mw <- attr(fit, "max_train_week"); if (is.null(mw) || !is.finite(mw)) mw <- max(design$d$.week, na.rm = TRUE)
  od <- offsets_df
  od$.week <- mw + od$horizon                # forecast week index per horizon (extrapolation point)
  predict_bayes_invasion(fit, od, design, horizons, affected_zones = affected_zones, rho = rho)
}

#' Per-zone import-force offsets log Lambda(h) at the FORECAST week for each
#' horizon, mirroring the frequentist workhorse's one-week-ahead projection: the
#' SOURCE-zone incidence is rolled forward with the local renewal R and a POINT
#' beta_proj (a nuisance projection), while the DESTINATION hazard beta is left
#' to the Bayesian posterior.
bayes_forecast_offsets <- function(zone_week_nc, mobility_matrices, gt_pmfs,
                                   covariates, osrm_mat, zones_all, mob, gt,
                                   cov_spec, horizons, beta_proj) {
  Y <- .count_wide(zone_week_nc, zones_all, "confirmed_nc")
  G <- daily_to_weekly_gt(gt_pmfs[[gt]]); W <- mobility_matrices[[mob]]
  static <- .static_features(covariates, zones_all)
  A_wide <- .count_wide(zone_week_nc, zones_all, "total_alerts")
  # Pad alerts with max(horizons) zero columns so an alert covariate evaluated at the forecast
  # week t_for = ncol(Yw)+h never indexes past A_wide (compute_foi has no upper bound guard).
  # Mirrors the frequentist workhorse's A_fore padding. Dormant today (the Bayesian grid omits
  # alert covariates) but prevents a silent out-of-bounds if the grid ever adds them.
  A_wide <- cbind(A_wide, matrix(0, nrow(A_wide), max(horizons),
                                 dimnames = list(rownames(A_wide), NULL)))
  # Suspected-but-not-confirmed leading-indicator matrix for the susp_* covariates,
  # padded with max(horizons) zero columns exactly like the alert matrix (susp is not
  # projected forward, so future weeks are zero). Active when the grid uses susp covariates.
  S_wide <- .susp_wide(zone_week_nc, zones_all)
  S_wide <- cbind(S_wide, matrix(0, nrow(S_wide), max(horizons),
                                 dimnames = list(rownames(S_wide), NULL)))
  R_local <- estimate_R_local(Y, G, ncol(Y)); Yw <- Y
  rows <- list()
  # Record EVERY intermediate week up to max(horizons), not just the requested horizons:
  # predict_bayes_invasion forms the cumulative hazard p(<=h)=1-exp(-sum_{h'<=h} mu(h')) by
  # summing the recorded per-week rows with horizon <= h. Recording only {1,2} happens to be
  # correct because they are contiguous; a non-contiguous request (e.g. c(1,3)) would silently
  # drop week 2 and under-count. Recording all weeks makes the cumulation correct for any set
  # (predict still emits only the requested horizons).
  for (h in seq_len(max(horizons))) {
    t_for <- ncol(Yw) + 1L
    Lam <- compute_foi(Yw, W, G, t_idx = t_for, zones_all)
    X <- .feature_matrix(t_for, Y, A_wide, W, G, static, osrm_mat, zones_all,
                         cov_spec, Y_active = Y, S_wide = S_wide)
    df <- data.frame(health_zone = zones_all, horizon = h,
                     Lambda = as.numeric(Lam),
                     logLam = log(pmax(as.numeric(Lam), 1e-12)))
    if (!is.null(X)) df <- cbind(df, as.data.frame(X[, , drop = FALSE]))
    rows[[length(rows) + 1L]] <- df
    own <- .gweighted_own(Yw, G, t_for)
    Yw  <- cbind(Yw, R_local * own + max(beta_proj, 0) * Lam)
  }
  dplyr::bind_rows(rows)
}

# --------------------------------------------------------------------------
# Shared draws head + summarise tail for the invasion predictors (review §2.1/§5.1).
# Factored out so predict_bayes_invasion (single GT) and predict_bayes_gt_marginal
# (GT-marginalised mixture, §2.1) summarise IDENTICALLY from posterior draws.
# --------------------------------------------------------------------------

#' Posterior DRAWS of the per-(zone,horizon) cumulative hazard and invasion
#' probability for one fitted model at one generation time. Returns a list keyed
#' by horizon ("h1","h2",...), each a list(zones, horizon, cum, p, p_inf) of
#' draws x zones matrices. This is the computational head previously inside
#' predict_bayes_invasion (see there for the link-reconciliation rationale).
.bayes_invasion_draws <- function(fit, offsets_df, design, horizons,
                                  rho = ASCERTAINMENT_NOMINAL) {
  if (is.null(fit)) return(NULL)
  feat <- design$feat
  nd <- offsets_df
  for (f in feat) if (f %in% names(nd)) {
    nd[[f]] <- (nd[[f]] - design$center[[f]]) / design$scale[[f]]
    nd[[f]][!is.finite(nd[[f]])] <- 0     # non-finite standardised covariate -> standardised mean
  }
  nd$.off <- nd$logLam; nd$logLam <- 0    # zero model offset; add per-horizon logLam back
  eta0 <- suppressWarnings(brms::posterior_linpred(
    fit, newdata = nd, allow_new_levels = TRUE))            # draws x rows (no offset)
  eta  <- sweep(eta0, 2, nd$.off, "+")                       # eta0 + logLam (link scale)
  # Reconcile observation links on the per-week HAZARD scale (see predict_bayes_invasion):
  #   cloglog: mu = exp(eta) (exact p = 1 - exp(-mu));  logit/probit: mu = -log(1 - link^-1(eta)).
  link <- tryCatch(attr(fit, "link") %||% fit$family$link, error = function(e) "cloglog")
  mu <- if (identical(link, "logit"))  -log1p(-stats::plogis(eta))
        else if (identical(link, "probit")) -log1p(-stats::pnorm(eta))
        else exp(eta)                                         # cloglog (default)
  zones <- unique(nd$health_zone)
  idx <- split(seq_len(nrow(nd)), nd$health_zone)
  lapply(stats::setNames(horizons, paste0("h", horizons)), function(h) {
    cum <- vapply(zones, function(z) {
      cols <- idx[[z]]; cols <- cols[nd$horizon[cols] <= h]
      if (!length(cols)) return(rep(0, nrow(mu)))
      rowSums(mu[, cols, drop = FALSE]) }, numeric(nrow(mu)))     # draws x zones
    list(zones = zones, horizon = h, cum = cum,
         p = 1 - exp(-cum), p_inf = 1 - exp(-cum / max(rho, 1e-3)))
  })
}

#' Summarise one horizon's posterior draws (draws x zones cum/p/p_inf) into the
#' per-zone tibble: posterior mean + median invasion probability, 90% CrI, sd,
#' infection-scale probability, and the RANK credible interval (review §5.1).
#' Affected zones are excluded from the ranking (left NA); the caller masks their
#' probabilities to NA.
.invasion_horizon_summary <- function(dr, affected_zones = character(0)) {
  zones <- dr$zones; p <- dr$p; cum <- dr$cum; p_inf <- dr$p_inf
  atrisk_cols <- which(!(zones %in% affected_zones))
  rank_med <- rep(NA_real_, length(zones))
  rank_lo  <- rep(NA_real_, length(zones))
  rank_hi  <- rep(NA_real_, length(zones))
  if (length(atrisk_cols) >= 1L) {
    rk <- apply(p[, atrisk_cols, drop = FALSE], 1,
                function(v) rank(-v, ties.method = "average"))    # (n_atrisk x draws)
    rk <- if (is.matrix(rk)) t(rk) else matrix(rk, nrow = nrow(p))  # -> draws x n_atrisk
    rank_med[atrisk_cols] <- apply(rk, 2, stats::median)
    rank_lo[atrisk_cols]  <- apply(rk, 2, stats::quantile, 0.05, names = FALSE)
    rank_hi[atrisk_cols]  <- apply(rk, 2, stats::quantile, 0.95, names = FALSE)
  }
  tibble::tibble(
    health_zone = zones, horizon = dr$horizon,
    mu_forecast = colMeans(cum), p_invasion = colMeans(p),
    p_median = apply(p, 2, stats::median),   # review §5.2: posterior median alongside the mean
    p_lo = apply(p, 2, stats::quantile, 0.05, names = FALSE),
    p_hi = apply(p, 2, stats::quantile, 0.95, names = FALSE),
    p_sd = apply(p, 2, stats::sd),
    p_infection_invasion = colMeans(p_inf),
    rank_med = rank_med, rank_lo = rank_lo, rank_hi = rank_hi)
}

#' Posterior invasion probability per (zone, horizon) with a 90% credible
#' interval. Cumulative hazard across horizons (p(<=h) = 1 - exp(-sum mu(h')));
#' the model offset is zeroed and log Lambda(h) added back explicitly so the
#' horizon-specific offset is unambiguous. Forecast covariates are standardised
#' with the fit design's center/scale; affected zones are masked to NA.
predict_bayes_invasion <- function(fit, offsets_df, design, horizons,
                                   affected_zones = character(0),
                                   rho = ASCERTAINMENT_NOMINAL) {
  if (is.null(fit)) return(NULL)
  # Draws head + per-horizon summarise tail are factored into shared helpers so the
  # GT-marginalised predictor (predict_bayes_gt_marginal, review §2.1) reuses the
  # IDENTICAL summarisation (mean/median/90% CrI + rank CrI) on pooled mixture draws.
  drl <- .bayes_invasion_draws(fit, offsets_df, design, horizons, rho = rho)
  if (is.null(drl)) return(NULL)
  res <- lapply(drl, .invasion_horizon_summary, affected_zones = affected_zones)
  dplyr::bind_rows(res) %>%
    # p_case_invasion mirrors p_invasion so EVERY Bayesian prediction (per-model, not
    # only the stacked one) carries the column the risk-score / map products key on.
    dplyr::mutate(p_case_invasion = p_invasion,
                  was_active_before = health_zone %in% affected_zones) %>%
    dplyr::mutate(dplyr::across(
      c(mu_forecast, p_invasion, p_case_invasion, p_median, p_lo, p_hi, p_sd, p_infection_invasion),
      ~ ifelse(was_active_before, NA_real_, .x)))
}

#' Pool a set of posterior-draw COMPONENTS into one mixture predictive and
#' summarise it (review §2.1/§2.4). Each component is a `.bayes_invasion_draws`
#' output (named-by-horizon list of list(zones, horizon, cum, p, p_inf)); the
#' mixture draws each component in proportion to its (renormalised) weight by
#' resampling JOINT draws — every sampled row keeps its whole zone vector, so
#' cross-zone correlation is preserved and the rank credible intervals are
#' correct. Shared by the GT-marginalised predictor and the ensemble.
#'
#' @param comp_draws named list of components (NULL components are dropped).
#' @param weights parallel numeric weights (renormalised internally).
#' @param affected_zones zones masked to NA (and excluded from ranking).
#' @param total_draws pooled draw count (default = one component's draw count).
#' @param seed RNG seed for reproducible resampling (snapshot/restore).
#' @return summarised prediction tibble with attr("weights_used").
.mix_and_summarise_draws <- function(comp_draws, weights, affected_zones = character(0),
                                     total_draws = NULL,
                                     seed = get0("RANDOM_SEED", ifnotfound = 20260704L)) {
  keep <- !vapply(comp_draws, is.null, logical(1))
  comp_draws <- comp_draws[keep]; weights <- as.numeric(weights)[keep]
  if (!length(comp_draws)) return(NULL)
  w <- weights; w[!is.finite(w) | w < 0] <- 0
  if (sum(w) <= 0) w <- rep(1, length(w))
  w <- w / sum(w)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  N <- nrow(comp_draws[[1]][[1]]$p)
  if (is.null(total_draws)) total_draws <- N
  m_k   <- pmax(1L, round(w * total_draws))
  hkeys <- names(comp_draws[[1]])
  pooled <- lapply(hkeys, function(hk) {
    parts <- lapply(seq_along(comp_draws), function(k) {
      dr <- comp_draws[[k]][[hk]]
      ix <- sample.int(nrow(dr$p), m_k[k], replace = m_k[k] > nrow(dr$p))
      list(cum = dr$cum[ix, , drop = FALSE], p = dr$p[ix, , drop = FALSE],
           p_inf = dr$p_inf[ix, , drop = FALSE])
    })
    list(zones   = comp_draws[[1]][[hk]]$zones,
         horizon = comp_draws[[1]][[hk]]$horizon,
         cum   = do.call(rbind, lapply(parts, `[[`, "cum")),
         p     = do.call(rbind, lapply(parts, `[[`, "p")),
         p_inf = do.call(rbind, lapply(parts, `[[`, "p_inf")))
  })
  out <- dplyr::bind_rows(lapply(pooled, .invasion_horizon_summary,
                                 affected_zones = affected_zones)) %>%
    dplyr::mutate(p_case_invasion = p_invasion,
                  was_active_before = health_zone %in% affected_zones) %>%
    dplyr::mutate(dplyr::across(
      c(mu_forecast, p_invasion, p_case_invasion, p_median, p_lo, p_hi, p_sd, p_infection_invasion),
      ~ ifelse(was_active_before, NA_real_, .x)))
  attr(out, "weights_used") <- stats::setNames(w, names(comp_draws))
  out
}

#' GT-MARGINALISED featured forecast (review §2.1 — marginalise, don't select).
#'
#' Instead of selecting a single generation time, this fits the given model at
#' EACH grid point of the GT PRIOR (make_gt_prior_pmfs) and forms the MIXTURE
#' predictive by resampling posterior draws from each grid point in proportion to
#' its prior weight, then summarises (mean/median/90% CrI + rank CrI) from the
#' pooled draws via the same tail as predict_bayes_invasion. Generation-time
#' uncertainty is thereby PROPAGATED into the featured forecast: intervals widen
#' relative to a fixed GT and no GT scenario is selected. The per-grid fit/predict
#' reuses build_invasion_design / fit_bayes_renewal / bayes_forecast_offsets
#' EXACTLY as fit_bayes_suite does (same argument pattern), at each grid GT key.
#'
#' Cost: one brms fit per GT grid point (default 15). Intended for the FEATURED
#' model's final forecast, not for every LFO fold (selection uses the medium
#' anchor; ranking is GT-robust). Resampling is seeded (RNG snapshot/restore) for
#' reproducibility. Returns NULL if no grid point fits.
#'
#' @param mob        mobility-kernel id.
#' @param cov_spec   covariate vector (character(0) for the intercept-only model).
#' @param gt_prior   a GT prior list (GT_PRIOR, or a GT_PRIOR_ALTS entry for the
#'                   sensitivity re-runs).
#' @param total_draws pooled draw count (default = one component's draw count).
#' @return prediction tibble (same columns as predict_bayes_invasion) with attrs
#'   `gt_grid` (grid + prior weights) and `gt_weights_used` (post-fit renormalised).
predict_bayes_gt_marginal <- function(zone_week_nc, mobility_matrices, covariates,
                                      osrm_mat, zones_all, mob, cov_spec = character(0),
                                      horizons = c(1L, 2L), affected_zones = character(0),
                                      gt_prior = GT_PRIOR, beta_proj = 0.05,
                                      rho = ASCERTAINMENT_NOMINAL, iter = 2000L,
                                      chains = 2L, link = "cloglog", total_draws = NULL,
                                      seed = get0("RANDOM_SEED", ifnotfound = 20260704L)) {
  gp   <- make_gt_prior_pmfs(gt_prior)                 # named daily PMFs + prior weights
  keys <- names(gp$gt_pmfs)
  cand <- unique(c(cov_spec, "log_pop", "ccvi", "d_min"))   # mirror fit_bayes_suite's candidates
  # --- fit + extract posterior draws at each GT grid point ---
  comp <- lapply(keys, function(k) {
    des <- tryCatch(build_invasion_design(zone_week_nc, mobility_matrices, gp$gt_pmfs,
             covariates, osrm_mat, zones_all, mob = mob, gt = k, candidates = cand),
             error = function(e) NULL)
    if (is.null(des)) return(NULL)
    fit <- fit_bayes_renewal(des, cov_spec = cov_spec, iter = iter, chains = chains,
                             seed = seed, link = link)
    if (is.null(fit)) return(NULL)
    off <- tryCatch(bayes_forecast_offsets(zone_week_nc, mobility_matrices, gp$gt_pmfs,
             covariates, osrm_mat, zones_all, mob = mob, gt = k, cov_spec = cov_spec,
             horizons = horizons, beta_proj = beta_proj), error = function(e) NULL)
    if (is.null(off)) return(NULL)
    drl <- .bayes_invasion_draws(fit, off, des, horizons, rho = rho)
    if (is.null(drl)) return(NULL)
    list(draws = drl, key = k, weight = as.numeric(gp$weights[[k]]))
  })
  comp <- comp[!vapply(comp, is.null, logical(1))]
  if (!length(comp)) { warning("[gt-marginal] no GT grid point fit; returning NULL"); return(NULL) }
  names(comp) <- vapply(comp, function(cc) cc$key, character(1))
  # Pool the per-GT posterior draws into the mixture predictive (shared machinery).
  out <- .mix_and_summarise_draws(
    lapply(comp, `[[`, "draws"), vapply(comp, function(cc) cc$weight, numeric(1)),
    affected_zones = affected_zones, total_draws = total_draws, seed = seed)
  attr(out, "gt_grid")         <- gp$grid
  attr(out, "gt_weights_used") <- attr(out, "weights_used")
  out
}

#' Kernel-diverse Bayesian ENSEMBLE as a proper MIXTURE PREDICTIVE (review §2.4).
#'
#' Pools the posterior draws of a small, pre-specified, structurally-diverse set of
#' member models (different mobility kernels, GT handled WITHIN each member by the
#' prior — never as an ensemble axis) into one coherent predictive, so the ensemble
#' carries a mean, MEDIAN, 90% CrI AND rank CrI from the mixture — not a
#' normal-approx of member summaries (cf. the legacy bayes_stacked_prediction). This
#' is the plan's response to the winner's-curse of selecting one of ~40 kernels:
#' combine kernel-diverse members instead. Weights are equal (robust default) or
#' loo-stacking; leakage-honest weights are the caller's responsibility (estimate on
#' training/inner folds, apply frozen).
#'
#' @param member_draws named list of `.bayes_invasion_draws` outputs (one per member).
#' @param weights named numeric member weights; default equal-weight over members.
#' @return mixture-predictive tibble (columns as predict_bayes_invasion) with
#'   attr("weights_used").
bayes_ensemble_mixture <- function(member_draws, weights = NULL,
                                   affected_zones = character(0), total_draws = NULL,
                                   seed = get0("RANDOM_SEED", ifnotfound = 20260704L)) {
  member_draws <- member_draws[!vapply(member_draws, is.null, logical(1))]
  if (!length(member_draws)) return(NULL)
  if (is.null(weights))
    weights <- stats::setNames(rep(1 / length(member_draws), length(member_draws)),
                               names(member_draws))
  weights <- weights[names(member_draws)]
  .mix_and_summarise_draws(member_draws, weights, affected_zones = affected_zones,
                           total_draws = total_draws, seed = seed)
}

#' Pairwise directed importation pressure / force of infection between health
#' zones, decomposed from the FEATURED Bayesian model's renewal-equation
#' structure at its posterior-median import coefficient.
#'
#' The model's per-week invasion hazard on destination i is
#'   mu_i(h) = beta * Lambda_i(h),   Lambda_i(h) = sum_j W[j,i] * source_j(h),
#'   source_j(h) = sum_k g(k) * Y_nc[j, t_for(h) - k]
#' (intercept-only cloglog renewal: p(<=h) = 1 - exp(-sum_{h'<=h} mu(h'))). The
#' DIRECTED contribution of origin j to at-risk destination i is therefore
#'   FOI_{j->i}(h) = beta * W[j,i] * source_j(h),
#' and by construction sum_j FOI_{j->i}(h) = mu_i(h) EXACTLY. `beta` is the
#' posterior-median import coefficient beta0 = median(exp(Intercept)) (the
#' Intercept `hr` in bayes_parameters), so the pairwise matrix IS the featured
#' model's posterior-median force of infection. The source projection for
#' horizons >= 2 reproduces bayes_forecast_offsets() verbatim (local renewal +
#' beta_proj * import), so Lambda_i(h) matches the model's own offsets.
#'
#' NOTE: exact for a covariate-free hazard (the featured M13-dist). If the
#' featured model carries covariates, beta is zone-varying and this scalar-beta
#' decomposition drops the per-zone covariate modulation of the DESTINATION
#' hazard — the caller warns and the renewal-share column stays valid, but `foi`
#' should then be read as the intercept-only import force.
#'
#' @return long tibble: one row per non-negligible directed (origin -> dest) pair
#'   per requested horizon, with the raw import force (beta-free), the force of
#'   infection (beta * import force), the destination totals and the renewal
#'   share of the destination hazard.
bayes_pairwise_import_force <- function(zone_week_nc, mobility_matrices, gt_pmfs,
                                        zones_all, mob, gt, beta_med, horizons,
                                        beta_proj = 0.05, affected_zones = character(0),
                                        province_map = NULL, tol = 1e-12) {
  if (is.null(beta_med) || !is.finite(beta_med) || beta_med < 0)
    stop("beta_med must be a finite non-negative scalar (posterior-median import coefficient)")
  W <- mobility_matrices[[mob]]
  if (is.null(W)) stop(sprintf("mobility matrix '%s' not found for pairwise FOI", mob))
  W <- W[zones_all, zones_all, drop = FALSE]
  G <- daily_to_weekly_gt(gt_pmfs[[gt]])
  Y <- .count_wide(zone_week_nc, zones_all, "confirmed_nc")
  R_local <- estimate_R_local(Y, G, ncol(Y)); Yw <- Y
  atrisk <- setdiff(zones_all, affected_zones)
  rows <- list()
  for (h in seq_len(max(horizons))) {
    t_for <- ncol(Yw) + 1L
    # GT-weighted past incidence of every ORIGIN j at this forecast week (== compute_foi's
    # internal Y_weighted); NA counts treated as 0 exactly as compute_foi does.
    src <- numeric(length(zones_all)); names(src) <- zones_all
    for (k in seq_along(G)) {
      tp <- t_for - k; if (tp < 1) break
      yv <- Yw[zones_all, tp]; yv[is.na(yv)] <- 0
      src <- src + G[k] * yv
    }
    Lambda <- as.numeric(t(W) %*% src); names(Lambda) <- zones_all   # dest total import force
    if (h %in% horizons) {
      origins <- zones_all[src > tol]
      if (length(origins)) {
        Wsub <- W[origins, atrisk, drop = FALSE]
        # N[j,i] = W[j,i] * src_j: recycle scales each column by the origin (row) source.
        Nsub <- Wsub * src[origins]
        keep <- which(Nsub > tol, arr.ind = TRUE)
        if (nrow(keep)) {
          oj <- origins[keep[, 1]]; di <- atrisk[keep[, 2]]
          impf <- Nsub[keep]; Ld <- Lambda[di]
          rows[[length(rows) + 1L]] <- tibble::tibble(
            origin_zone = oj, dest_zone = di, horizon = h,
            w_ji = Wsub[keep], source_origin = src[oj],
            import_force = impf, foi = beta_med * impf,
            dest_import_force_total = Ld,
            dest_hazard_week = beta_med * Ld,
            share_of_dest = ifelse(Ld > 0, impf / Ld, NA_real_))
        }
      }
    }
    # Project sources forward for the next horizon EXACTLY as bayes_forecast_offsets().
    own <- .gweighted_own(Yw, G, t_for)
    Yw  <- cbind(Yw, R_local * own + max(beta_proj, 0) * Lambda)
  }
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(out)
  if (!is.null(province_map) && all(c("nom", "province") %in% names(province_map))) {
    pm <- province_map %>% dplyr::distinct(.data$nom, .keep_all = TRUE)
    out <- out %>%
      dplyr::left_join(dplyr::transmute(pm, origin_zone = .data$nom, origin_province = .data$province),
                       by = "origin_zone") %>%
      dplyr::left_join(dplyr::transmute(pm, dest_zone = .data$nom, dest_province = .data$province),
                       by = "dest_zone")
  }
  out %>% dplyr::arrange(.data$horizon, .data$dest_zone, dplyr::desc(.data$foi))
}

#' Posterior parameter summary (hazard ratios, 90% CrI) for a fitted model.
#' The Intercept row is the baseline import coefficient beta0 = exp(Intercept);
#' covariate rows are per-SD hazard ratios exp(beta_m).
bayes_param_table <- function(fit, label = "") {
  if (is.null(fit)) return(NULL)
  dr <- posterior::as_draws_df(fit)
  vars <- grep("^b_", names(dr), value = TRUE)
  # exp(coef) is a HAZARD ratio under cloglog (the renewal default) but an ODDS ratio under
  # logit; record the scale so the combined bayes_parameters.csv is self-describing and the
  # one logit model's ORs are not silently read as HRs.
  lk <- tryCatch(attr(fit, "link") %||% fit$family$link, error = function(e) "cloglog")
  escale <- if (identical(lk, "cloglog")) "hazard ratio"
            else if (identical(lk, "logit")) "odds ratio"
            else if (identical(lk, "probit")) "probit exp(coef)" else lk
  rh  <- tryCatch(brms::rhat(fit), error = function(e) NA_real_)   # hoist: compute rhat once
  rmx <- suppressWarnings(max(rh[grep("Intercept|b_", names(rh))], na.rm = TRUE))
  if (!is.finite(rmx)) rmx <- NA_real_
  purrr::map_dfr(vars, function(v) {
    x <- exp(dr[[v]]); term <- sub("^b_", "", v)
    tibble::tibble(model = label, term = term, is_intercept = term == "Intercept",
      effect_scale = escale,
      hr = stats::median(x),
      lo = stats::quantile(x, 0.05, names = FALSE),
      hi = stats::quantile(x, 0.95, names = FALSE),
      p_dir = if (term == "Intercept") NA_real_ else mean(x > 1),  # posterior P(effect>1)
      rhat = rmx)
  })
}

#' NUTS/HMC convergence diagnostics for a fitted brms model (review §4.4).
#'
#' The reviewer asked "what Stan settings and what convergence diagnostics?" — the
#' pipeline previously extracted only Rhat. This returns a one-row-per-model tibble
#' with divergent-transition count and %, minimum bulk/tail effective sample size,
#' maximum Rhat over the population parameters (Intercept + b_*), and the post-warmup
#' draw count. Every extraction is guarded (backend-agnostic: cmdstanr or rstan) so a
#' diagnostics failure never aborts a fit; unavailable quantities return NA.
#'
#' @param fit   a brmsfit, or NULL.
#' @param label model label carried into the row.
#' @return one-row tibble (model, n_divergent, pct_divergent, rhat_max,
#'   ess_bulk_min, ess_tail_min, n_draws), or NULL when `fit` is NULL.
bayes_fit_diagnostics <- function(fit, label = "") {
  if (is.null(fit)) return(NULL)
  # Divergent transitions via brms::nuts_params (works for both backends).
  np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
  n_div <- NA_integer_; n_iter <- NA_integer_
  if (!is.null(np) && all(c("Parameter", "Value") %in% names(np))) {
    dv <- np[np$Parameter == "divergent__", , drop = FALSE]
    if (nrow(dv)) { n_div <- as.integer(sum(dv$Value, na.rm = TRUE)); n_iter <- nrow(dv) }
  }
  pct_div <- if (is.finite(n_div) && is.finite(n_iter) && n_iter > 0L)
    100 * n_div / n_iter else NA_real_
  # Rhat / bulk+tail ESS over the population coefficients only (Intercept + b_*).
  dr <- tryCatch(posterior::as_draws_df(fit), error = function(e) NULL)
  rhat_max <- NA_real_; ess_bulk_min <- NA_real_; ess_tail_min <- NA_real_; n_draws <- NA_integer_
  if (!is.null(dr)) {
    n_draws <- tryCatch(as.integer(posterior::ndraws(dr)), error = function(e) NA_integer_)
    sm <- tryCatch(posterior::summarise_draws(dr, "rhat", "ess_bulk", "ess_tail"),
                   error = function(e) NULL)
    if (!is.null(sm)) {
      sm <- sm[grepl("^b_", sm$variable), , drop = FALSE]
      if (nrow(sm)) {
        rhat_max     <- suppressWarnings(max(sm$rhat, na.rm = TRUE))
        ess_bulk_min <- suppressWarnings(min(sm$ess_bulk, na.rm = TRUE))
        ess_tail_min <- suppressWarnings(min(sm$ess_tail, na.rm = TRUE))
      }
    }
  }
  fix <- function(x) if (is.finite(x)) x else NA_real_
  tibble::tibble(
    model = label, n_divergent = n_div, pct_divergent = pct_div,
    rhat_max = fix(rhat_max), ess_bulk_min = fix(ess_bulk_min),
    ess_tail_min = fix(ess_tail_min), n_draws = n_draws)
}

#' FULL posterior DRAWS of each fitted model's parameters on the exp(coef) scale (beta0 =
#' exp(Intercept) for the intercept row; hazard/odds ratio for covariate rows), in long form
#' (model, term, is_intercept, .draw, hr). Used to visualise the posterior DISTRIBUTIONS (not
#' just the median + interval). Thinned to keep the tibble compact.
bayes_posterior_draws <- function(fits, max_draws = 1500L) {
  if (is.null(fits) || !length(fits)) return(NULL)
  purrr::map_dfr(names(fits), function(nm) {
    fit <- fits[[nm]]; if (is.null(fit)) return(NULL)
    dr <- tryCatch(posterior::as_draws_df(fit), error = function(e) NULL)
    if (is.null(dr)) return(NULL)
    vars <- grep("^b_", names(dr), value = TRUE)
    if (!length(vars)) return(NULL)
    lk <- tryCatch(attr(fit, "link") %||% fit$family$link, error = function(e) "cloglog")
    escale <- if (identical(lk, "cloglog")) "hazard ratio"
              else if (identical(lk, "logit")) "odds ratio" else lk
    idx <- if (nrow(dr) > max_draws) round(seq(1, nrow(dr), length.out = max_draws)) else seq_len(nrow(dr))
    purrr::map_dfr(vars, function(v) {
      term <- sub("^b_", "", v)
      tibble::tibble(model = nm, term = term, is_intercept = term == "Intercept",
                     effect_scale = escale, .draw = seq_along(idx), hr = exp(dr[[v]][idx]))
    })
  })
}

#' The comprehensive Bayesian model grid — the SINGLE SOURCE OF TRUTH used for BOTH
#' the current-forecast suite and the leave-future-out cross-validation, so every
#' Bayesian model that is fit is also cross-validated. It spans all four modelling
#' axes the frequentist suite varies:
#' The suite is COMPOSED from a CORE that is always fitted plus OPTIONAL families gated
#' by the 00_config.R toggles (INCLUDE_OSRM_DIST_MODELS / INCLUDE_GEO_COV_MODELS /
#' INCLUDE_M11_MODELS / INCLUDE_LOGIT_SENS_MODELS / INCLUDE_SUSPECTED_COV_MODELS), so the
#' default grid is focused while any family can be switched back on without code change.
#'   * MOBILITY kernel: the CORE four distinct families — gravity (M4), composite-gravity
#'     (M8), multi-kernel ensemble (M9), radiation-composite (M10). (M1/M2/M3/M5/M6/M7/M4b
#'     are components or minor variants folded into these composites.) OPTIONAL: the inward
#'     meeting-location FOI kernel (M11; INCLUDE_M11_MODELS) and, for each family, an OSRM
#'     ROAD-DISTANCE (km) deterrence variant (M4-dist / ... / M11-dist; INCLUDE_OSRM_DIST_MODELS)
#'     — the same kernels keyed on distance instead of travel time. -dist entries are also
#'     auto-filtered out when the road-distance matrices were not built.
#'   * GENERATION TIME: a SINGLE medium anchor for the selection grid (review §2.1 — the GT
#'     is no longer a selectable axis; it is a PRIOR marginalised into the featured forecast).
#'   * COVARIATES: none (CORE) and FULL-exogenous (log_pop, CCVI, d_min, healthsite_density;
#'     always kept). OPTIONAL: the reduced "geo" subset (log_pop, CCVI, d_min;
#'     INCLUDE_GEO_COV_MODELS). NEW (default ON, INCLUDE_SUSPECTED_COV_MODELS): suspected-
#'     but-not-confirmed leading-indicator covariates — susp_local (own preceding-week
#'     suspected cases) and susp_import (mobility-weighted import of other zones' suspected
#'     cases), as a susp-only model and a wider FULL+susp model. positivity and the
#'     total-alert covariates remain excluded: positivity = confirmed/tests is circular with
#'     the invasion outcome; the suspected covariate is the line-list classification signal,
#'     distinct from the raw alert count and used past-only (causal in LFO).
#'   * OBSERVATION PROCESS: cloglog (default, principled — makes the import force a
#'     proper log-cumulative-hazard offset, giving p = 1 - exp(-beta*Lambda)). OPTIONAL:
#'     a logit-link SENSITIVITY on the FULL covariate set (INCLUDE_LOGIT_SENS_MODELS);
#'     prediction reconciles both on the per-week hazard scale. (Poisson vs NB is the
#'     COUNT observation for the source-incidence projection, shared with and varied in
#'     the frequentist suite, not an axis of the binary invasion likelihood.)
bayes_default_grid <- function(mobility_matrices = NULL) {
  # Optional model families (00_config.R toggles; get0 fallbacks so a config-less caller
  # — e.g. a unit test — still gets the documented defaults).
  # inc_dist is the MASTER switch for EVERY -dist variant (generic M4/M8/M9/M10-dist AND cohort
  # M13/M14-dist AND consensus M17-dist): OFF => the grid carries a single distance measure
  # (travel time). inc_full toggles the FULL-exogenous covariate models (geo + healthsite_density);
  # OFF => the covariate sweep is {none, geo} only.
  inc_dist  <- isTRUE(get0("INCLUDE_OSRM_DIST_MODELS",     ifnotfound = FALSE))
  inc_geo   <- isTRUE(get0("INCLUDE_GEO_COV_MODELS",       ifnotfound = FALSE))
  inc_full  <- isTRUE(get0("INCLUDE_FULL_COV_MODELS",      ifnotfound = FALSE))
  inc_m11   <- isTRUE(get0("INCLUDE_M11_MODELS",           ifnotfound = FALSE))
  inc_logit <- isTRUE(get0("INCLUDE_LOGIT_SENS_MODELS",    ifnotfound = FALSE))
  inc_susp  <- isTRUE(get0("INCLUDE_SUSPECTED_COV_MODELS", ifnotfound = TRUE))
  inc_flowstat <- isTRUE(get0("INCLUDE_FLOWSTATIC_MODELS", ifnotfound = TRUE))
  # M9 (short-trip + M4/M5/M6a ensemble) and M15 (symmetrised inflow+outflow static) are OFF by
  # default — their matrices are not built (03_mobility_matrices.R), so the trailing availability
  # filter would drop them anyway; gating here keeps the grid self-documenting.
  inc_m9    <- isTRUE(get0("INCLUDE_M9_MODELS",            ifnotfound = FALSE))
  inc_m15   <- isTRUE(get0("INCLUDE_M15_MODELS",           ifnotfound = FALSE))
  geo  <- c("log_pop", "ccvi", "d_min")
  full <- c("log_pop", "ccvi", "d_min", "healthsite_density")
  susp <- c("susp_import", "susp_local")
  g <- list()
  add <- function(mob, gt, cov, link, label)
    g[[length(g) + 1L]] <<- list(mob = mob, gt = gt, cov = cov, link = link, label = label)

  # --- CORE (always): the distinct FOI families at medium GT, no covariates ---
  add("M4",  "medium", character(0), "cloglog", "Bayes-M4-med")
  add("M8",  "medium", character(0), "cloglog", "Bayes-M8-med")
  if (inc_m9) add("M9", "medium", character(0), "cloglog", "Bayes-M9-med")
  add("M10", "medium", character(0), "cloglog", "Bayes-M10-med")
  # --- Flowminder cohort composites (M13/M14 + geographic-distance -dist), default ON ---
  # Cohort presence source rows + gravity/radiation, on travel-time (M13/M14) or road-km
  # geographic distance (M13/M14-dist). Auto-filtered below when a matrix is absent
  # (INCLUDE_COHORT_MODELS off, or no road-distance matrix). The -dist pair does NOT depend
  # on INCLUDE_OSRM_DIST_MODELS (its matrices build under the cohort family).
  if (isTRUE(get0("INCLUDE_COHORT_MODELS", ifnotfound = TRUE))) {
    add("M13",      "medium", character(0), "cloglog", "Bayes-M13-med")
    add("M14",      "medium", character(0), "cloglog", "Bayes-M14-med")
    # Road-distance cohort twins only when the -dist axis is on (inc_dist).
    if (inc_dist) {
      add("M13-dist", "medium", character(0), "cloglog", "Bayes-M13-dist")
      add("M14-dist", "medium", character(0), "cloglog", "Bayes-M14-dist")
    }
    # Covariate variants: FULL-exogenous when inc_full, the reduced geo subset when inc_geo — on the
    # travel-time (M13/M14) and, when inc_dist, the road-km (M13/M14-dist) cohort composites.
    # Auto-filtered below when a matrix is absent.
    for (cm in c("M13", "M14", if (inc_dist) c("M13-dist", "M14-dist"))) {
      if (inc_full) add(cm, "medium", full, "cloglog", paste0("Bayes-", cm, "-full"))
      if (inc_geo)  add(cm, "medium", geo,  "cloglog", paste0("Bayes-", cm, "-geo"))
    }
  }
  # --- Flowminder combined inflow+outflow static family (M15/M16/M17 + M17-dist), default ON ---
  # M15 = combined-Flowminder static (symmetrised inflow+outflow); M16 = cohort + static (cohort
  # source rows where available, M15 elsewhere); M17/M17-dist = grand all-kernel consensus ensemble.
  # Swept EXHAUSTIVELY: each kernel x generation time {short, medium, long} x covariates {none,
  # geo (optional), full}. Auto-filtered below when a matrix is absent (M16 needs the cohort kernel;
  # M17-dist needs a road-distance matrix), so a cohort-off / distance-off run simply omits them.
  if (inc_flowstat) {
    for (mm in c(if (inc_m15) "M15", "M16", "M17", if (inc_dist) "M17-dist")) {
      # Single medium GT ANCHOR (GT is marginalised into the featured forecast, not a grid axis).
      # none + (full when inc_full) + (geo when inc_geo) covariate variants.
      add(mm, "medium", character(0), "cloglog", paste0("Bayes-", mm, "-med"))
      if (inc_full) add(mm, "medium", full, "cloglog", paste0("Bayes-", mm, "-full"))
      if (inc_geo)  add(mm, "medium", geo,  "cloglog", paste0("Bayes-", mm, "-geo"))
    }
  }
  # NOTE (review §2.1): the former Bayes-M8-short / Bayes-M8-long and the M15/M16/M17
  # short/long GT-scenario models are REMOVED. Selecting the generation time by
  # cross-validation is "akin to fitting" a quantity the data cannot identify; instead a
  # single medium anchor is used for the SELECTION grid and the generation time is
  # marginalised over its prior (GT_PRIOR) in the featured forecast. This also shrinks the
  # candidate space (fewer models selected among ~39 events; §2.2). GT SENSITIVITY is
  # reported by re-running the featured forecast under GT_PRIOR_ALTS, not as grid members.
  # --- FULL-exogenous covariate models (OPTIONAL; inc_full, OFF by default). The reduced geo
  #     subset is added separately under inc_geo below. ---
  if (inc_full) {
    add("M8",  "medium", full, "cloglog", "Bayes-M8-full")
    add("M4",  "medium", full, "cloglog", "Bayes-M4-full")
    if (inc_m9) add("M9", "medium", full, "cloglog", "Bayes-M9-full")
    add("M10", "medium", full, "cloglog", "Bayes-M10-full")
  }
  # --- NEW suspected-but-not-confirmed leading-indicator covariate models (default ON) ---
  # susp_import = mobility-weighted import of OTHER zones' preceding suspected (not yet
  # confirmed) cases; susp_local = own preceding-week suspected cases. Both past-only
  # (causal in LFO). One susp-only model and one FULL+susp "wider covariate" model.
  if (inc_susp) {
    add("M8", "medium", susp,                   "cloglog", "Bayes-M8-susp")
    add("M8", "medium", unique(c(full, susp)),  "cloglog", "Bayes-M8-full-susp")
  }
  # --- OPTIONAL: reduced "geo" covariate models (M4/M8/M9/M10 base kernels) ---
  if (inc_geo) {
    add("M8",  "medium", geo, "cloglog", "Bayes-M8-geo")
    add("M4",  "medium", geo, "cloglog", "Bayes-M4-geo")
    if (inc_m9) add("M9", "medium", geo, "cloglog", "Bayes-M9-geo")
    add("M10", "medium", geo, "cloglog", "Bayes-M10-geo")
  }
  # --- OPTIONAL: inward / meeting-location FOI kernel (M11) ---
  if (inc_m11) {
    add("M11", "medium", character(0), "cloglog", "Bayes-M11-inward")
    add("M11", "medium", full,         "cloglog", "Bayes-M11-full")
    if (inc_geo) add("M11", "medium", geo, "cloglog", "Bayes-M11-geo")
  }
  # --- OPTIONAL: OSRM ROAD-DISTANCE (km) deterrence kernels — same families keyed on
  #     distance instead of travel time. Auto-filtered below if the -dist matrices are absent. ---
  if (inc_dist) {
    for (mm in c("M4-dist", "M8-dist", if (inc_m9) "M9-dist", "M10-dist")) {
      add(mm, "medium", character(0), "cloglog", paste0("Bayes-", mm))
      if (inc_full) add(mm, "medium", full, "cloglog", paste0("Bayes-", mm, "-full"))
      if (inc_geo)  add(mm, "medium", geo,  "cloglog", paste0("Bayes-", mm, "-geo"))
    }
    if (inc_m11) {
      add("M11-dist", "medium", character(0), "cloglog", "Bayes-M11-dist")
      add("M11-dist", "medium", full,         "cloglog", "Bayes-M11-dist-full")
      if (inc_geo) add("M11-dist", "medium", geo, "cloglog", "Bayes-M11-dist-geo")
    }
  }
  # --- OPTIONAL: logit-link observation-process SENSITIVITY (slow to fit). Rides on the
  #     retained FULL covariate set so it is directly comparable to Bayes-M8-full; the
  #     geo-logit variant is added too when the geo set is enabled. Prediction reconciles
  #     the logit hazard with the cloglog default on the per-week hazard scale. ---
  if (inc_logit) {
    add("M8", "medium", full, "logit", "Bayes-M8-full-logit")
    if (inc_geo) add("M8", "medium", geo, "logit", "Bayes-M8-geo-logit")
  }

  if (!is.null(mobility_matrices)) g <- Filter(function(x) !is.null(mobility_matrices[[x$mob]]), g)
  g
}

#' Fit the Bayesian SUITE across the comprehensive grid (mobility x GT x covariates x
#' observation link), predict the current forecast for each, collect posterior
#' parameters, and compute loo predictive-stacking weights.
fit_bayes_suite <- function(zone_week_nc, mobility_matrices, gt_pmfs, covariates,
                            osrm_mat, zones_all, affected_zones, horizons,
                            grid = NULL, iter = 2000L) {
  if (is.null(grid)) grid <- bayes_default_grid(mobility_matrices)
  grid <- Filter(function(g) !is.null(mobility_matrices[[g$mob]]), grid)

  # ---- Pass 1: build every design (cheap; no MCMC yet) ----------------------
  designs <- list()
  for (g in grid) {
    des <- tryCatch(build_invasion_design(zone_week_nc, mobility_matrices, gt_pmfs,
             covariates, osrm_mat, zones_all, mob = g$mob, gt = g$gt,
             candidates = unique(c(g$cov, "log_pop", "ccvi", "d_min"))),
             error = function(e) NULL)
    if (!is.null(des)) designs[[g$label]] <- list(g = g, des = des)
  }
  if (!length(designs))
    return(list(fits = list(), preds = tibble::tibble(), params = tibble::tibble(),
                weights = NULL, stacked = FALSE, grid = grid))

  # loo predictive stacking requires every model's pointwise log-likelihood to be
  # over the SAME observations. Each design's at-risk set (!affected & Lambda>0)
  # depends on the mobility kernel (M8 short-trip rows are sparse; M4 gravity is
  # dense), so the raw designs have different row counts and loo_model_weights()
  # would abort. We therefore fit every suite model on the COMMON at-risk zones
  # (intersection of the designs' zones), keeping the full-design center/scale so
  # standardisation is unchanged. Predictions are still made on each model's full
  # forecast grid (bayes_forecast_offsets), so forecast coverage is not reduced —
  # only the rows used to estimate the stacking weights are aligned.
  # Key each at-risk observation by (zone, week). loo_model_weights matches the
  # per-model pointwise log-likelihoods BY POSITION, so we both intersect to the
  # common observations AND sort every model's fit data into the same canonical key
  # order below, otherwise the weights would be computed on mismatched rows.
  key_of <- function(dd) paste(dd$.zone, dd$.week, sep = "@")
  common_keys <- Reduce(intersect, lapply(designs, function(x) key_of(x$des$d)))
  # Also require the grid covariates to be non-NA on the shared rows: brms drops
  # rows with any NA in a formula term, so a covariate (geo) model would silently
  # fit on FEWER complete cases than an intercept-only model even on identical
  # (zone,week) keys — which breaks loo's equal-observations requirement. Restrict
  # to rows that are complete for every covariate any suite model uses.
  grid_covs <- unique(unlist(lapply(designs, function(x) x$g$cov)))
  if (length(grid_covs) && length(common_keys)) {
    for (x in designs) {
      dd <- x$des$d; present <- intersect(grid_covs, names(dd))
      if (length(present))
        common_keys <- intersect(common_keys, key_of(dd)[stats::complete.cases(dd[, present, drop = FALSE])])
    }
  }
  message(sprintf("[bayes] suite: %d models; %d common complete-case (zone,week) rows for stacking alignment",
                  length(designs), length(common_keys)))

  # ---- Pass 2: fit each model on the aligned row-set, predict on the full grid --
  # Every model is independent. Reuse the future pool configured by run_all.R;
  # PARALLEL_JOBS=1 retains the historical sequential execution and ordering.
  fit_one <- function(nm) {
    g <- designs[[nm]]$g; des <- designs[[nm]]$des
    message(sprintf("[bayes] fitting %s (mobility %s, %d covariates) ...",
                    g$label, g$mob, length(g$cov)))
    des_fit <- des
    if (length(common_keys) >= 10L) {           # align rows when a usable common set exists
      keep <- key_of(des$d) %in% common_keys
      des_fit$d <- des$d[keep, , drop = FALSE]
      des_fit$d <- des_fit$d[order(key_of(des_fit$d)), , drop = FALSE]   # canonical order for loo
      des_fit$n_events <- sum(des_fit$d$invaded); des_fit$n_obs <- nrow(des_fit$d)
    }
    fit <- fit_bayes_renewal(des_fit, cov_spec = g$cov, iter = iter,
                             link = g$link %||% "cloglog")
    if (is.null(fit)) return(NULL)
    pr <- NULL
    off <- tryCatch(bayes_forecast_offsets(zone_week_nc, mobility_matrices, gt_pmfs,
             covariates, osrm_mat, zones_all, g$mob, g$gt, g$cov, horizons,
             beta_proj = des$beta0 %||% 0.05), error = function(e) NULL)
    if (!is.null(off)) {
      pr <- predict_bayes_invasion(fit, off, des, horizons, affected_zones)
      if (!is.null(pr)) pr$method <- g$label
    }
    # Compute the model's aligned pointwise LOO here, while the suite workers are
    # already fanned out. Previously this was a second sequential pass over every
    # fit and took several minutes after the parallel sampling had finished.
    lo <- NULL
    if (length(common_keys) >= 10L) {
      dd <- des$d
      dd <- dd[key_of(dd) %in% common_keys, , drop = FALSE]
      dd <- dd[order(key_of(dd)), , drop = FALSE]
      lo <- tryCatch({
        llm <- brms::log_lik(fit, newdata = dd)
        suppressWarnings(loo::loo(llm))
      }, error = function(e) NULL)
    }
    list(fit = fit, pred = pr, params = bayes_param_table(fit, g$label), loo = lo)
  }
  .par <- get0("PARALLEL_JOBS", ifnotfound = 1L) > 1L &&
          length(designs) > 1L && requireNamespace("furrr", quietly = TRUE)
  .tf <- Sys.time()
  fit_res <- if (.par)
    furrr::future_map(names(designs), fit_one,
                      .options = furrr::furrr_options(seed = TRUE))
  else lapply(names(designs), fit_one)
  names(fit_res) <- names(designs)
  message(sprintf("[bayes-timing] current suite: %d models %s in %.1fs",
                  length(designs), if (.par) "PARALLEL" else "seq",
                  as.numeric(difftime(Sys.time(), .tf, units = "secs"))))
  fit_res <- Filter(Negate(is.null), fit_res)
  fits  <- lapply(fit_res, `[[`, "fit")
  preds <- lapply(fit_res, `[[`, "pred")
  preds <- Filter(Negate(is.null), preds)
  params <- lapply(fit_res, `[[`, "params")

  # ---- loo predictive stacking (Yao et al. 2018); honest fallback -----------
  # Each worker computed its model's pointwise log-likelihood on its own aligned rows
  # (the common complete-case (zone,week) set, in a canonical order, carrying that
  # model's own mobility offset and covariates). Relying on each brms fit's INTERNAL
  # loo is fragile — the fits use different mobility offsets and one may transiently
  # fail — and yields loo objects of differing dimension that loo_model_weights
  # rejects. Extracting log_lik on identical observations guarantees aligned
  # dimensions, so stacking on the shared evaluation set is well defined. The expensive
  # extraction is part of the parallel pass above; only the cheap weight calculation
  # remains here. Falls back to an equal-weight average if any LOO calculation failed.
  wts <- NULL; stacked <- FALSE
  if (length(fits) >= 2 && length(common_keys) >= 10L) {
    wts <- tryCatch({
      ll_list <- lapply(fit_res, `[[`, "loo")
      if (length(ll_list) != length(fits) || any(vapply(ll_list, is.null, logical(1))))
        stop("one or more aligned model LOO calculations failed")
      # Collapse each model to its per-observation ELPD (draws integrated out) and
      # stack on that [obs x model] matrix. loo_model_weights() compares full
      # [draws x obs] matrices and rejects models with DIFFERING DRAW COUNTS — which
      # happens when a cmdstanr chain partially fails; stacking_weights() on the
      # pointwise ELPD is immune because the draws are already integrated out.
      lpd <- sapply(ll_list, function(l) l$pointwise[, "elpd_loo"])
      if (!is.matrix(lpd) || ncol(lpd) < 2L) stop("insufficient aligned models for stacking")
      w <- loo::stacking_weights(lpd)
      stacked <- TRUE
      stats::setNames(as.numeric(w), names(fits)) },
      error = function(e) { warning("[bayes] loo stacking failed (", conditionMessage(e),
                                    "); reporting equal-weight average."); NULL })
  }
  list(fits = fits, preds = dplyr::bind_rows(preds),
       params = dplyr::bind_rows(params), weights = wts, stacked = stacked, grid = grid)
}

#' Data-informed posterior over the generation-time MEAN for one model spec. The GT is a fixed
#' ASSUMPTION in the renewal model, so we PROFILE it: refit the model over a grid of GT means
#' (fixed coefficient of variation), score each fit by its leave-one-out predictive density
#' (loo elpd), and combine with a weakly-informative literature prior on the GT mean to obtain
#' a posterior over it (a pseudo-Bayesian / loo-weighted model average — NOT a fully joint GT
#' estimate, which would require the renewal convolution inside Stan). Returns a tibble
#' (gt_mean, elpd, elpd_se, log_prior, weight) plus posterior summaries as attributes.
bayes_gt_posterior <- function(zone_week_nc, mobility_matrices, gt_pmfs, covariates, osrm_mat,
                               zones_all, mob = "M8", cov = character(0), link = "cloglog",
                               gt_means = seq(12, 20, by = 1), gt_cv = 0.6, max_tau = 50L,
                               prior_mean = 15.3, prior_sd = 1.5, iter = 1000L) {
  # gt_cv and max_tau are held FIXED across the grid ON PURPOSE: that keeps the weekly GT support
  # (all lags 1..ceil(max_tau/7) strictly positive) constant, so the at-risk set !aff & Lam>0 is
  # GT-INVARIANT and every fit's loo is over the SAME observations — a precondition for comparing
  # them. This is asserted below.
  fit_one <- function(m) {
    pmf <- tryCatch(make_gt_pmf(mean = m, sd = m * gt_cv, max_tau = max_tau), error = function(e) NULL)
    if (is.null(pmf)) return(NULL)
    des <- tryCatch(build_invasion_design(zone_week_nc, mobility_matrices, list(g = pmf),
             covariates, osrm_mat, zones_all, mob = mob, gt = "g",
             candidates = unique(c(cov, "log_pop", "ccvi", "d_min"))), error = function(e) NULL)
    if (is.null(des) || is.null(des$d) || !nrow(des$d)) return(NULL)
    fit <- fit_bayes_renewal(des, cov_spec = cov, iter = iter, link = link)
    if (is.null(fit)) return(NULL)
    lo <- tryCatch({ suppressWarnings(loo::loo(brms::log_lik(fit))) }, error = function(e) NULL)
    if (is.null(lo)) return(NULL)
    list(loo = lo, meta = tibble::tibble(
      gt_mean = m, elpd = lo$estimates["elpd_loo", "Estimate"],
      elpd_se = lo$estimates["elpd_loo", "SE"], n_obs = nrow(lo$pointwise),
      max_pareto_k = suppressWarnings(max(lo$diagnostics$pareto_k, na.rm = TRUE)),
      log_prior = stats::dnorm(m, prior_mean, prior_sd, log = TRUE)))
  }
  .par <- get0("PARALLEL_JOBS", ifnotfound = 1L) > 1L &&
          length(gt_means) > 1L && requireNamespace("furrr", quietly = TRUE)
  .tf <- Sys.time()
  res <- if (.par)
    furrr::future_map(gt_means, fit_one,
                      .options = furrr::furrr_options(seed = TRUE))
  else lapply(gt_means, fit_one)
  message(sprintf("[bayes-timing] GT profile: %d fits %s in %.1fs",
                  length(gt_means), if (.par) "PARALLEL" else "seq",
                  as.numeric(difftime(Sys.time(), .tf, units = "secs"))))
  res <- Filter(Negate(is.null), res)
  loos <- lapply(res, `[[`, "loo")
  meta <- lapply(res, `[[`, "meta")
  if (length(loos) < 2L) return(NULL)
  out <- dplyr::bind_rows(meta)
  # ENFORCE the identical-observations invariant (else the elpd comparison is invalid).
  if (length(unique(out$n_obs)) != 1L) {
    warning("[bayes] GT profile fits have differing observation counts (",
            paste(unique(out$n_obs), collapse = "/"), "); loo comparison not well-posed — returning NULL.")
    return(NULL)
  }
  if (any(out$max_pareto_k > 0.7, na.rm = TRUE))
    warning("[bayes] GT profile: some loo Pareto-k > 0.7 (rare-event cloglog LOO unreliable); interpret the GT preference with caution.")
  # pseudo-BMA+ : Bayesian-bootstrap the POINTWISE elpd so the elpd standard error is respected
  # (plain pseudo-BMA on the summed elpd is over-confident and collapses to one grid point — Yao
  # et al. 2018 §3.2). Then reweight by the GT-mean prior (changing the model prior from uniform
  # to Normal(prior_mean, prior_sd)). Falls back to summed-elpd softmax if the bootstrap fails.
  w_data <- tryCatch(as.numeric(loo::loo_model_weights(loos, method = "pseudobma", BB = TRUE)),
    error = function(e) { s <- exp(out$elpd - max(out$elpd)); s / sum(s) })
  w_prior <- exp(out$log_prior - max(out$log_prior))
  out$weight <- w_data * w_prior; out$weight <- out$weight / sum(out$weight)
  attr(out, "pref_mean") <- sum(out$gt_mean * out$weight)   # loo-preference mean, NOT a true posterior mean
  attr(out, "prior")     <- c(mean = prior_mean, sd = prior_sd)
  out
}

#' SENSITIVITY of the featured Bayesian model to the two-stage nowcast INPUT. The suite feeds
#' epinowcast-corrected training counts into the renewal fit as a fixed second stage, which
#' ignores nowcast uncertainty; a fully joint model would estimate the nowcast and the invasion
#' hazard together. As a pragmatic check we REFIT the same model spec on several nowcast
#' treatments of the training counts (each `zw_variants[[nm]]` must carry a `confirmed_nc`
#' column) and return the posterior draws of beta0 per scenario, so a reader can see whether the
#' nowcast choice materially shifts the inference. Returns tibble(scenario, .draw, beta0, ...).
bayes_nowcast_sensitivity <- function(zw_variants, mobility_matrices, gt_pmfs, covariates,
                                      osrm_mat, zones_all, mob = "M8", gt = "medium",
                                      cov = character(0), link = "cloglog", iter = 1000L) {
  fit_one <- function(nm) {
    zw <- zw_variants[[nm]]; if (is.null(zw) || !nrow(zw)) return(NULL)
    des <- tryCatch(build_invasion_design(zw, mobility_matrices, gt_pmfs, covariates, osrm_mat,
             zones_all, mob = mob, gt = gt,
             candidates = unique(c(cov, "log_pop", "ccvi", "d_min"))), error = function(e) NULL)
    if (is.null(des) || is.null(des$d) || !nrow(des$d)) return(NULL)
    fit <- fit_bayes_renewal(des, cov_spec = cov, iter = iter, link = link)
    if (is.null(fit)) return(NULL)
    dr <- tryCatch(posterior::as_draws_df(fit), error = function(e) NULL)
    if (is.null(dr) || !"b_Intercept" %in% names(dr)) return(NULL)
    tibble::tibble(scenario = nm, .draw = seq_len(nrow(dr)),
                   beta0 = exp(dr$b_Intercept),
                   n_events = des$n_events %||% sum(des$d$invaded),
                   n_obs = des$n_obs %||% nrow(des$d))
  }
  .par <- get0("PARALLEL_JOBS", ifnotfound = 1L) > 1L &&
          length(zw_variants) > 1L && requireNamespace("furrr", quietly = TRUE)
  .tf <- Sys.time()
  out <- if (.par)
    furrr::future_map(names(zw_variants), fit_one,
                      .options = furrr::furrr_options(seed = TRUE))
  else lapply(names(zw_variants), fit_one)
  message(sprintf("[bayes-timing] nowcast sensitivity: %d fits %s in %.1fs",
                  length(zw_variants), if (.par) "PARALLEL" else "seq",
                  as.numeric(difftime(Sys.time(), .tf, units = "secs"))))
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(NULL)
  dplyr::bind_rows(out)
}

#' loo-stacked (or equal-weight fallback) Bayesian ensemble of the suite's
#' current-forecast predictions.
bayes_stacked_prediction <- function(suite) {
  if (is.null(suite$preds) || !nrow(suite$preds)) return(NULL)
  w <- suite$weights
  if (is.null(w)) { nm <- unique(suite$preds$method); w <- stats::setNames(rep(1 / length(nm), length(nm)), nm) }
  lab <- if (isTRUE(suite$stacked) && !is.null(suite$weights)) "Bayes-stack (loo)" else "Bayes-stack (equal-weight)"
  # loo stacking (Yao 2018) is a LINEAR POOL of the member predictive distributions. The
  # pooled MEAN is the weighted mean of member means, but the pooled INTERVAL is NOT the
  # weighted mean of member 5/95% bounds (that ignores between-model disagreement and
  # under-covers when models differ). Reconstruct the pool's spread via the law of total
  # variance from each member's posterior mean + SD (p_invasion, p_sd), so the stacked 90%
  # CrI widens when the members disagree, and take a normal-approx 90% interval on [0,1].
  suite$preds %>% dplyr::filter(method %in% names(w)) %>%
    dplyr::mutate(wt = w[method]) %>%
    dplyr::group_by(health_zone, horizon, was_active_before) %>%
    dplyr::summarise(
      # .var_pool MUST be computed BEFORE p_invasion is redefined: dplyr::summarise()
      # exposes a just-created summary column to later expressions in the same call, so
      # if p_invasion were summarised first, the name would resolve to the length-1
      # pooled scalar here (not the per-member vector) and p_invasion[ok] would return
      # NAs for every group with >=2 members — silently NA-ing the stacked CrI.
      .var_pool = { ok <- is.finite(p_invasion) & is.finite(wt) &
                          (if ("p_sd" %in% names(suite$preds)) is.finite(p_sd) else TRUE)
                    if (!any(ok)) NA_real_ else {
                      wn <- wt[ok] / sum(wt[ok]); pm <- sum(wn * p_invasion[ok])
                      within <- if ("p_sd" %in% names(suite$preds)) sum(wn * p_sd[ok]^2) else 0
                      within + sum(wn * (p_invasion[ok] - pm)^2) } },
      p_invasion = stats::weighted.mean(p_invasion, wt, na.rm = TRUE),
      mu_forecast = stats::weighted.mean(mu_forecast, wt, na.rm = TRUE),
      p_infection_invasion = stats::weighted.mean(p_infection_invasion, wt, na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::mutate(
      p_lo = pmax(0, p_invasion - 1.645 * sqrt(.var_pool)),   # 90% CrI (normal approx to the pool)
      p_hi = pmin(1, p_invasion + 1.645 * sqrt(.var_pool))) %>%
    dplyr::select(-.var_pool) %>%
    # p_case_invasion mirrors p_invasion so the Bayesian stacked tibble carries the
    # SAME column the risk-map / bars / priority products key on (compute_risk_scores,
    # plot_invasion_risk_map, plot_risk_scores_bars) — otherwise the Bayesian decision
    # maps are silently skipped.
    dplyr::mutate(method = lab, p_case_invasion = p_invasion)
}

#' LFO closure: refit the Bayesian model on each fold's training data and predict
#' — plugs into run_invasion_lfo exactly like a frequentist model. cmdstanr
#' reuses the compiled Stan binary across folds (same formula), so only sampling
#' is repeated. Kept lightweight (fewer iterations) for cross-validation.
make_bayes_lfo_model <- function(mob, gt, cov_spec, mobility_matrices, gt_pmfs,
                                 covariates, osrm_mat, zones_all, iter = 1000L,
                                 link = "cloglog") {
  # FORCE the per-spec arguments so each closure captures ITS OWN mob/gt/cov/link.
  # Without this, R's lazy promises defer evaluation until the LFO runs the closure —
  # by which time the build loop has finished and `sp` holds the LAST spec, so every
  # model silently computes identical (last-spec) predictions.
  force(mob); force(gt); force(cov_spec); force(link); force(iter)
  function(zw, ti, hz, cut) {
    des <- tryCatch(build_invasion_design(zw, mobility_matrices, gt_pmfs, covariates,
             osrm_mat, zones_all, mob = mob, gt = gt,
             candidates = unique(c(cov_spec, "log_pop", "ccvi", "d_min", "healthsite_density"))),
             error = function(e) NULL)
    if (is.null(des) || des$n_events < 1L) return(NULL)
    fit <- fit_bayes_renewal(des, cov_spec = cov_spec, iter = iter, chains = 2L, link = link)
    if (is.null(fit)) return(NULL)
    off <- tryCatch(bayes_forecast_offsets(zw, mobility_matrices, gt_pmfs, covariates,
             osrm_mat, zones_all, mob, gt, cov_spec, hz, beta_proj = des$beta0 %||% 0.05),
             error = function(e) NULL)
    if (is.null(off)) return(NULL)
    predict_bayes_invasion(fit, off, des, hz, affected_zones = character(0))
  }
}

#' Refit a Bayesian renewal model at each fold cutoff (causally, training only on
#' week <= cutoff) to trace the POSTERIOR import coefficient beta0 and covariate
#' hazard ratios over time — the Bayesian analogue of compute_beta_over_folds /
#' compute_params_over_time, carrying proper posterior 90% credible intervals. A
#' covariate spec is used so BOTH beta0 and covariate-HR traces are produced.
compute_bayes_params_over_time <- function(zone_week_outbreak, cutoffs, mobility_matrices,
                                           gt_pmfs, covariates, osrm_mat, zones_all,
                                           mob = "M8", gt = "medium",
                                           cov_spec = c("log_pop", "ccvi", "d_min"),
                                           iter = 800L, nowcast_fn = NULL) {
  if (is.null(nowcast_fn)) nowcast_fn <- get0("apply_nowcast_correction")
  fit_one <- function(ci) {
    cut <- cutoffs[ci]            # index, NOT `for (cut in cutoffs)` — iterating a
                                  # Date vector unclasses cut to numeric (throws on R<4.3)
    zw  <- zone_week_outbreak %>% dplyr::filter(week_start <= cut)
    zwn <- tryCatch(nowcast_fn(zw, analysis_date = cut + 7), error = function(e) zw)
    des <- tryCatch(build_invasion_design(zwn, mobility_matrices, gt_pmfs, covariates,
             osrm_mat, zones_all, mob = mob, gt = gt,
             candidates = unique(c(cov_spec, "log_pop", "ccvi", "d_min"))),
             error = function(e) NULL)
    if (is.null(des) || des$n_events < 1L) return(NULL)
    fit <- tryCatch(fit_bayes_renewal(des, cov_spec = cov_spec, iter = iter, chains = 2L),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    tb <- bayes_param_table(fit, "")
    ic <- tb %>% dplyr::filter(is_intercept)
    bl <- if (nrow(ic)) tibble::tibble(cutoff = as.Date(cut),
      beta0 = ic$hr[1], beta0_lo = ic$lo[1], beta0_hi = ic$hi[1]) else NULL
    cv <- tb %>% dplyr::filter(!is_intercept)
    pl <- if (nrow(cv)) cv %>%
      dplyr::transmute(cutoff = as.Date(cut), term, hr, lo, hi)
      else NULL
    list(params = pl, beta = bl)
  }
  .par <- get0("PARALLEL_JOBS", ifnotfound = 1L) > 1L &&
          length(cutoffs) > 1L && requireNamespace("furrr", quietly = TRUE)
  .tf <- Sys.time()
  res <- if (.par)
    furrr::future_map(seq_along(cutoffs), fit_one,
                      .options = furrr::furrr_options(seed = TRUE))
  else lapply(seq_along(cutoffs), fit_one)
  message(sprintf("[bayes-timing] params over time: %d folds %s in %.1fs",
                  length(cutoffs), if (.par) "PARALLEL" else "seq",
                  as.numeric(difftime(Sys.time(), .tf, units = "secs"))))
  res <- Filter(Negate(is.null), res)
  pl <- Filter(Negate(is.null), lapply(res, `[[`, "params"))
  bl <- Filter(Negate(is.null), lapply(res, `[[`, "beta"))
  list(params = if (length(pl)) dplyr::bind_rows(pl) else NULL,
       beta   = if (length(bl)) dplyr::bind_rows(bl) else NULL)
}

message("[bayes] 21_bayesian_renewal.R loaded — Bayesian renewal invasion suite (brms; cloglog/logit links, full grid).")
