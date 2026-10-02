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
  # BOUND THE INIT (2026-09-17). The guard below used to catch only a NON-FINITE log(beta0);
  # a large FINITE one was passed straight through, and that is the dominant failure mode.
  # beta0 is the unpenalised null-model MLE, which for kernels whose mobility mass is
  # concentrated runs to 10^2-10^125. At init the largest-offset row then has
  # mu = beta0 * exp(max logLam); once mu exceeds ~37, 1 - exp(-mu) is EXACTLY 1 in double
  # precision, so log(1 - p) is -Inf for every non-invaded row and Stan returns no draws
  # ("Unable to retrieve the metadata"). bayes_fit_failures.log records this 28 times across
  # 8 specs, and EVERY logged failure has beta0 >= 32.7.
  #
  # Cap the intercept so the largest initial hazard is O(1) (max eta <= 0 => mu <= 1 =>
  # p <= 0.632), and floor it so it cannot be absurdly small. The posterior is nowhere near
  # this value anyway — the prior shrinks beta0 by 3-4 orders of magnitude from the MLE.
  # ORDER MATTERS: floor FIRST, then cap. Applying the floor last meant max(lb0, -20)
  # OVERRODE the cap whenever max(logLam) > 20 — it would hand back -20 where the cap
  # required -max(logLam) — and the stated guarantee "max eta <= 0 => mu <= 1" silently
  # failed in exactly the regime the cap exists for. (max(logLam) is ~1.74 on the shipped
  # designs, so this is latent, not live.) When the two conflict the CAP must win: an init
  # that is merely small costs a few warmup steps, an init that overflows the hazard returns
  # no draws at all.
  lb0 <- suppressWarnings(log(design$beta0))
  if (!is.finite(lb0)) lb0 <- -3
  lb0 <- max(lb0, -20)
  .mx_off <- suppressWarnings(max(design$d$logLam[is.finite(design$d$logLam)]))
  if (is.finite(.mx_off)) lb0 <- min(lb0, -.mx_off)
  .init_at <- function(v) lapply(seq_len(chains), function(i)
    if (length(feat)) list(Intercept = v, b = rep(0, length(feat))) else list(Intercept = v))
  .brm <- function(ad, v) suppressMessages(suppressWarnings(brms::brm(
    fml, data = d, family = brms::bernoulli(link = link),
    prior = .bayes_prior(feat), iter = iter, warmup = iter %/% 2L,
    chains = chains, cores = min(chains, 4L), seed = seed, init = .init_at(v),
    backend = "cmdstanr", refresh = 0, silent = 2,
    control = list(adapt_delta = ad, max_treedepth = 12L))))
  fit <- tryCatch(.brm(adapt_delta, lb0), error = function(e) {
    # The retry now VARIES THE INIT as well as adapt_delta. Reusing the same init_list meant an
    # init-caused failure was guaranteed to recur, so the "covers a marginal init" claim was
    # false for exactly the failures that dominate the log. -3 is the prior mean.
    tryCatch(.brm(max(adapt_delta, 0.97), -3), error = function(e2) {
      .bayes_log_fit_failure(design, feat, link, conditionMessage(e2)); NULL }) })
  if (!is.null(fit)) { attr(fit, "feat") <- feat; attr(fit, "cov_spec") <- cov_spec
                       attr(fit, "link") <- link }
  fit
}

# ===========================================================================
# TIME-VARYING IMPORT COEFFICIENT beta_t  (reviewer point E1)
# ===========================================================================
# The fixed-beta model asserts that one import->invasion conversion rate held from the
# index case to the present. That is the assumption most at risk as an outbreak matures:
# the epicentre saturates, control changes, and the zones still at risk are progressively
# further from the source. Every variant below is a STRICT GENERALISATION — each collapses
# exactly onto fit_bayes_renewal() as its extra parameter goes to zero — so the posterior,
# not an assertion, decides whether beta varies.
#
#   cloglog P(invade_{i,t}) = beta0 + u_t + x_i' gamma + log Lambda_{i,t},
#   beta_t = exp(beta0 + u_t).
#
# FIVE PROCESSES FOR u_t, AND WHAT EACH DOES AT THE FORECAST WEEK. The forecast behaviour
# is the whole point of the comparison: these models are used to predict week T+h, so what
# u_{T+h} does when it has never been observed is not a technicality.
#
#   tv     process                                  u_{T+h} | data            interval
#   -----  ---------------------------------------  ----------------------    --------------
#   trend  u_t = gamma_w * z(t)   (a covariate)      extrapolates the slope    no extra width
#   week   u_t ~ N(0, sigma) iid                     reverts to 0 (= beta0)    +sigma
#   rw1    u_t = u_{t-1} + eps_t, eps ~ N(0, sigma)  PERSISTS u_T              +sqrt(h)*sigma
#   ar1    Cov(u_t,u_s) = sigma^2 exp(-|t-s|/l)      decays u_T toward 0       grows to sigma
#   gp     Cov(u_t,u_s) = sigma^2 exp(-(t-s)^2/2l^2) smooth local extrapolation grows to sigma
#
# WHY ALL FIVE RATHER THAN ONE. They differ only in how much of the CURRENT level they
# carry into the forecast, and that is exactly the unknown: `week` carries none, `rw1`
# carries all of it, `ar1` carries a fitted fraction phi = exp(-1/l) per week and nests
# both (l -> Inf is rw1, l -> 0 is week), `gp` is `ar1` made smooth (differentiable, so it
# extrapolates the local slope a short way), and `trend` is the parametric limit in which
# the drift is deterministic and persists for ever. Cross-validation prices them.
#
# HOW THEY ARE FITTED. All in brms, none in hand-written Stan:
#   * rw1 is a group-level effect with a KNOWN covariance, (1 | gr(.week_f, cov = K)) with
#     K[i,j] = min(tau_i, tau_j). That is the Brownian/RW1 covariance exactly, and brms
#     scales it by sd^2, so the `sd` parameter IS the innovation SD sigma.
#   * ar1 and gp are EXACT Gaussian processes, gp(.week, cov = "exponential" | "exp_quad")
#     with k = NA. The exponential kernel is the Ornstein-Uhlenbeck process, i.e. AR(1) in
#     continuous time, so phi = exp(-1/l) per week. `gr = TRUE` evaluates the GP on the
#     ~12-18 UNIQUE weeks, not on the ~8,000 rows, so an exact GP is cheap here.
#
# WHY NOT THE APPROXIMATE GP. An earlier version of this file fitted `gp(.week, k = 5,
# c = 5/4)`. brms' basis-function approximation is built on a bounded domain
# L = c * max|x - mean(x)|, and the FORECAST weeks fall OUTSIDE it, where the basis is
# being evaluated off the domain it was constructed for. The exact GP has no domain: brms
# predicts a new week from the proper conditional, which is why `k` is left at NA here.
# (Verified: predictive SD at week T+1 and T+2 exceeds that at week T, as a GP must.)
#
# WHY THE rw1 FORECAST IS PROPAGATED BY HAND. brms can draw a new level of a group-level
# term, but only from the POPULATION distribution N(0, sd) — correct for `week`, and wrong
# for `rw1`, whose forecast must start from u_T. .tv_linpred_fn() therefore computes the
# population part with re_formula = NA and adds the propagated walk itself; see there.

#' Covariance matrix of a first-order random walk on the supplied week indices.
#'
#' Returns K with K[i,j] = min(tau_i, tau_j), tau = week - min(week) + 1, and the week
#' labels as dimnames (brms matches `cov` to factor levels BY NAME). Under brms'
#' parameterisation Cov(u) = sd^2 * K, so
#'     Var(u_t) = sd^2 * tau_t,   Var(u_t - u_s) = sd^2 * |tau_t - tau_s|,
#' which is a Brownian motion with per-week innovation SD = sd, started from u ~ N(0, sd^2)
#' at the FIRST training week.
#'
#' ANCHORING AT THE FIRST OBSERVED WEEK, not at week 0, is deliberate. `.week` counts from
#' the start of the epidemic grid and each cross-validation fold starts at a different
#' index, so min(tau_i, tau_j) on the RAW index would make the prior variance of the walk
#' depend on where the fold happens to begin — a fold starting at week 8 would get eight
#' innovations' worth of prior spread before its first observation. Anchoring makes the
#' model invariant to the index origin, which is what an RW1 prior should be.
#'
#' Gaps are handled correctly by construction: a two-week gap contributes 2*sd^2 of
#' variance, exactly as a random walk should.
.tv_rw1_cov <- function(weeks) {
  w <- sort(unique(as.integer(weeks)))
  tau <- w - min(w) + 1L
  K <- outer(tau, tau, pmin)
  storage.mode(K) <- "double"
  dimnames(K) <- list(as.character(w), as.character(w))
  K
}

#' Time-VARYING-beta Bayesian renewal fit.
#'
#' @param tv one of "week" (iid weekly intercept), "rw1" (random walk), "ar1"
#'   (Ornstein-Uhlenbeck / continuous-time AR(1)) or "gp" (squared-exponential GP).
#'   "trend" is NOT handled here: it rides the ordinary covariate machinery as `week_idx`
#'   and must go through fit_bayes_renewal(), or the trend would be fitted twice over.
#' @inheritParams fit_bayes_renewal
#' @return a brmsfit with attr("tv") and attr("max_train_week"), or NULL.
#'
#' PRIORS. sigma (the group-level `sd` for rw1/week, `sdgp` for ar1/gp) carries
#' half-normal(0, 1) on the log-hazard scale — the same scale as BAYES_PRIOR_COEF, so a
#' ~1 log-unit week-to-week swing is well inside the prior and a wild one is not. The GP
#' lengthscale keeps brms' default inverse-gamma, which is tuned to the observed spacing
#' and range of `.week`: it puts negligible mass on lengthscales shorter than the sampling
#' interval (unidentifiable) or far longer than the record (indistinguishable from an
#' intercept shift). With ~12-18 weekly points that admits phi = exp(-1/l) roughly in
#' [0.4, 0.95], which is the range the data can speak to.
#'
#' SMALL-SAMPLE HONESTY. There are ~54 invasion events over ~18 weeks. A posterior sigma
#' concentrating near 0 ("week-to-week variation is not warranted") is a legitimate and
#' informative result, and is reported as such rather than being read as a failure.
#'
#' Falls back to the fixed-beta fit when there is no .week column or < 4 distinct training
#' weeks: below that the process parameters are estimated from too few points to mean
#' anything, and for the GP variants the lengthscale is not identified at all.
fit_bayes_renewal_tv <- function(design, cov_spec = character(0),
                                 iter = 2000L, chains = 2L,
                                 seed = get0("RANDOM_SEED", ifnotfound = 20260704L),
                                 adapt_delta = 0.95, link = "cloglog",
                                 tv = "week") {
  d <- design$d; if (is.null(d) || !nrow(d)) return(NULL)
  tv <- match.arg(tv, c("week", "rw1", "ar1", "gp"))
  # Guard: every variant needs a weekly index and enough distinct weeks for its process
  # parameter to be estimable. 4 rather than the previous 3 because ar1/gp additionally
  # estimate a lengthscale.
  if (!".week" %in% names(d) || length(unique(d$.week)) < 4L) {
    warning(sprintf("[tv] design lacks .week or has < 4 training weeks; %s falls back to the fixed-beta fit",
                    tv), call. = FALSE)
    return(fit_bayes_renewal(design, cov_spec, iter, chains, seed, adapt_delta, link))
  }
  # brms needs the grouping variable to be a factor; keep the ORIGINAL integer in `.week`
  # (build_invasion_design's column, used by the GP terms and by .tv_offsets) and group on
  # a factor copy. Levels are set from the SORTED unique weeks so they line up with the
  # row/column order of the RW1 covariance matrix.
  wk_levels <- as.character(sort(unique(as.integer(d$.week))))
  d$.week_f <- factor(as.integer(d$.week), levels = wk_levels)
  feat <- intersect(cov_spec, design$feat)
  rhs  <- if (length(feat)) paste(feat, collapse = " + ") else NULL
  data2 <- NULL
  tv_term <- switch(tv,
    week = "(1 | .week_f)",
    rw1  = "(1 | gr(.week_f, cov = Krw))",
    ar1  = "gp(.week, cov = \"exponential\", gr = TRUE, k = NA, scale = FALSE)",
    gp   = "gp(.week, cov = \"exp_quad\",    gr = TRUE, k = NA, scale = FALSE)")
  if (identical(tv, "rw1")) data2 <- list(Krw = .tv_rw1_cov(d$.week))
  terms <- c("1", tv_term, rhs)
  fml   <- stats::as.formula(sprintf("invaded ~ %s + offset(logLam)",
                                     paste(terms[!is.na(terms)], collapse = " + ")))
  # The scale prior is on `sd` for the group-level variants and on `sdgp` for the GPs;
  # brms rejects a prior on a class the model does not have, so it is selected here.
  tv_prior <- .bayes_prior(feat) +
    brms::prior_string("normal(0, 1)",
                       class = if (tv %in% c("ar1", "gp")) "sdgp" else "sd")
  # SAME BOUNDED INIT as fit_bayes_renewal(); see there for the full argument. beta0 is the
  # unpenalised null MLE and runs to 1e2-1e125 on concentrated kernels; once
  # mu = beta0 * exp(max logLam) exceeds ~37, 1 - exp(-mu) is EXACTLY 1 in double precision,
  # log(1 - p) is -Inf for every non-invaded row, and Stan returns no draws at all.
  # ORDER MATTERS: floor first, then cap, so the cap wins when the two conflict.
  lb0 <- suppressWarnings(log(design$beta0)); if (!is.finite(lb0)) lb0 <- -3
  lb0 <- max(lb0, -20)
  .mx_off <- suppressWarnings(max(design$d$logLam[is.finite(design$d$logLam)]))
  if (is.finite(.mx_off)) lb0 <- min(lb0, -.mx_off)
  .init_at <- function(v) lapply(seq_len(chains), function(i)
    if (length(feat)) list(Intercept = v, b = rep(0, length(feat))) else list(Intercept = v))
  .brm <- function(ad, v) suppressMessages(suppressWarnings(brms::brm(
    fml, data = d, data2 = data2, family = brms::bernoulli(link = link),
    prior = tv_prior, iter = iter, warmup = iter %/% 2L,
    chains = chains, cores = min(chains, 4L), seed = seed, init = .init_at(v),
    backend = "cmdstanr", refresh = 0, silent = 2,
    control = list(adapt_delta = ad, max_treedepth = 12L))))
  # The retry VARIES THE INIT as well as adapt_delta; reusing the same init guaranteed that an
  # init-caused failure recurred identically. -3 is the prior mean.
  fit <- tryCatch(.brm(adapt_delta, lb0), error = function(e)
    tryCatch(.brm(max(adapt_delta, 0.99), -3), error = function(e2) {
      .bayes_log_fit_failure(design, feat, link,
                             paste0("[tv:", tv, "] ", conditionMessage(e2))); NULL }))
  if (!is.null(fit)) { attr(fit, "feat") <- feat; attr(fit, "cov_spec") <- cov_spec
                       attr(fit, "link") <- link; attr(fit, "tv") <- tv
                       attr(fit, "max_train_week") <- max(d$.week, na.rm = TRUE)
                       attr(fit, "train_weeks") <- as.integer(wk_levels) }
  fit
}

#' Is this tv type a STRUCTURAL term handled by fit_bayes_renewal_tv()?
#'
#' "trend" is NOT: it is the ordinary `week_idx` covariate, so it must take the fixed-beta
#' fitter or the trend would be fitted twice over. "none" is not. Everything else is.
#' ONE predicate, used at every dispatch site (the suite, the LFO closure, the
#' GT-marginalised featured path), because three hand-written `identical(tv, "week")` tests
#' is exactly how a newly added process silently gets refit as fixed-beta.
.tv_is_structural <- function(tv) (tv %||% "none") %in% c("week", "rw1", "ar1", "gp")

#' Does this tv type need the hand-propagated forecast path?
#' Only rw1: `week` is handled by brms' population draw for a new level, and ar1/gp are
#' exact GPs whose conditional at a new week brms computes itself.
.tv_needs_manual_forecast <- function(tv) identical(tv %||% "none", "rw1")

#' Forecast-week grouping/index columns for a tv fit's offsets.
#'
#' Factored out so the single-scale and both-scales predictors construct the forecast week
#' the SAME way; a second copy of this arithmetic could drift and would then give a tv model
#' two different forecast weeks depending on which predictor was called.
#'
#' `.week` (integer) is what the GP terms read; `.week_f` is what the group-level terms read.
#' BOTH are set for every tv type, so brms' newdata validation never fails on a missing
#' column regardless of which term the fitted formula happens to carry.
#'
#' For rw1 the factor level is pinned to the LAST TRAINING WEEK rather than to the forecast
#' week. That is not a fudge: the rw1 prediction path calls posterior_linpred(re_formula = NA),
#' so the level is never read — but brms still validates newdata against the grouping factor,
#' and a genuinely new level there would (a) require the covariance matrix to be extended and
#' (b) invite brms to sample it from N(0, sd), which is exactly the wrong forecast. The
#' propagated walk is added explicitly in .tv_linpred_fn().
.tv_offsets <- function(fit, offsets_df, design) {
  mw <- attr(fit, "max_train_week")
  if (is.null(mw) || !is.finite(mw)) mw <- max(design$d$.week, na.rm = TRUE)
  tv <- attr(fit, "tv") %||% "none"
  od <- offsets_df
  od$.week <- mw + od$horizon              # forecast week index per horizon
  if (.tv_needs_manual_forecast(tv)) {
    tw <- attr(fit, "train_weeks")
    if (is.null(tw)) tw <- sort(unique(as.integer(design$d$.week)))
    od$.week_f <- factor(as.integer(mw), levels = as.character(tw))
  } else {
    # `week`: the forecast week IS a new level, which is the point — brms draws its u from
    # N(0, sigma_w) under sample_new_levels = "gaussian". ar1/gp carry no grouping term, so
    # the column is inert for them; it is set anyway so one code path serves all types.
    od$.week_f <- factor(od$.week)
  }
  attr(od, "projection") <- attr(offsets_df, "projection")
  od
}

#' Linear-predictor function for a fitted model, or NULL to use the default brms call.
#'
#' Returns a function(fit, nd) -> draws x rows matrix of the linear predictor EXCLUDING the
#' offset (the caller zeroes `logLam` and adds it back per horizon).
#'
#' ONLY rw1 needs one. Its forecast is
#'     u_{T+k} = u_T + sum_{m=1..k} eps_m,   eps_m ~ N(0, sigma) iid,
#' so the level PERSISTS and the variance grows linearly in k. Three properties are load-
#' bearing and are why this is not done with a per-row rnorm():
#'   * WITHIN A DRAW, the innovations are shared across rows, so u_{T+1} and u_{T+2} are
#'     correctly correlated. The cumulative hazard summed over horizons in
#'     .bayes_invasion_draws() would otherwise average out the very uncertainty this model
#'     exists to represent.
#'   * ACROSS ZONES, the SAME u is added to every zone at a given horizon, because u_t is a
#'     time effect, not a zone effect. Drawing it per row would turn a common shock into
#'     independent noise and shrink the between-zone rank uncertainty.
#'   * The draw is REPRODUCIBLE and does not disturb the caller's RNG stream: the seed is
#'     derived from RANDOM_SEED and the previous .Random.seed is restored on exit, matching
#'     the convention used elsewhere in this suite.
.tv_linpred_fn <- function(fit) {
  tv <- attr(fit, "tv") %||% "none"
  if (!.tv_needs_manual_forecast(tv)) return(NULL)
  function(fit, nd) {
    eta <- suppressWarnings(brms::posterior_linpred(fit, newdata = nd, re_formula = NA))
    mw  <- attr(fit, "max_train_week")
    dr  <- posterior::as_draws_df(fit)
    v_u <- sprintf("r_.week_f[%s,Intercept]", mw)
    v_s <- "sd_.week_f__Intercept"
    if (!all(c(v_u, v_s) %in% names(dr))) {
      warning("[tv:rw1] the fitted walk's last level or its innovation SD is missing from the ",
              "posterior; the forecast REVERTS to beta0 instead of persisting the walk.",
              call. = FALSE)
      return(eta)
    }
    u_T   <- as.numeric(dr[[v_u]])
    sigma <- as.numeric(dr[[v_s]])
    nd_wk <- as.integer(nd$.week)
    k     <- nd_wk - as.integer(mw)          # steps ahead per row (>= 1 on forecast rows)
    kmax  <- max(k, 0L)
    ndr   <- length(u_T)
    if (nrow(eta) != ndr)
      stop(sprintf("[tv:rw1] %d linear-predictor draws but %d parameter draws; cannot pair them.",
                   nrow(eta), ndr), call. = FALSE)
    # Reproducible innovations, caller's RNG stream preserved.
    .seed <- get0("RANDOM_SEED", ifnotfound = 20260704L)
    .had  <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    .old  <- if (.had) get(".Random.seed", envir = globalenv()) else NULL
    on.exit({ if (.had) assign(".Random.seed", .old, envir = globalenv())
              else suppressWarnings(rm(".Random.seed", envir = globalenv())) }, add = TRUE)
    set.seed(.seed + 7411L)
    # cumsum of iid N(0, sigma) innovations: draws x kmax, then u_{T+k} = u_T + walk[, k].
    walk <- if (kmax >= 1L)
      t(apply(matrix(stats::rnorm(ndr * kmax), nrow = ndr) * sigma, 1, cumsum)) else NULL
    if (!is.null(walk) && kmax == 1L) walk <- matrix(walk, ncol = 1L)
    add <- vapply(seq_along(k), function(j) {
      if (k[j] <= 0L) u_T else u_T + walk[, k[j]]
    }, numeric(ndr))
    if (!is.matrix(add)) add <- matrix(add, nrow = ndr)
    eta + add
  }
}

predict_bayes_invasion_tv <- function(fit, offsets_df, design, horizons,
                                      affected_zones = character(0),
                                      delta = NULL) {
  if (is.null(fit)) return(NULL)
  predict_bayes_invasion(fit, .tv_offsets(fit, offsets_df, design), design, horizons,
                         affected_zones = affected_zones, delta = delta,
                         linpred_fn = .tv_linpred_fn(fit))
}

predict_bayes_invasion_tv_scales <- function(fit, offsets_df, design, horizons,
                                             affected_zones = character(0),
                                             delta = NULL) {
  if (is.null(fit)) return(NULL)
  predict_bayes_invasion_scales(fit, .tv_offsets(fit, offsets_df, design), design, horizons,
                                affected_zones = affected_zones, delta = delta,
                                linpred_fn = .tv_linpred_fn(fit))
}

# --------------------------------------------------------------------------
# Source projection for horizons >= 2: EpiNow2 R draws x the model's own hazard draws
# --------------------------------------------------------------------------

#' Posterior draws of the national reproduction number over one week, from EpiNow2.
#'
#' Drives the local-renewal term of the h >= 2 source projection (.bayes_project_mu()).
#' Fitted with .epinow2_rt_fit() (02_epi_params.R) to the national daily confirmed series
#' as KNOWN on `issue_date`, built by epinow2_daily_confirmed() (02_epi_params.R): dated
#' EXACTLY as the weekly counts (`date_index`: the usable onset, or an onset imputed from the
#' sample date) — the only dating consistent with EpiNow2's onset->sample truncation model — and
#' censored on observation date with the same rule as the weekly training counts
#' (linelist_observation_date(), 22_daily_reissue.R). The truncation delay is
#' effective_onset_sample_delay() (the shared full-snapshot estimate), and the generation-time
#' parameters (Gamma mean, SD and support) are those that generated `gt_pmf`. EpiNow2
#' discretises that Gamma itself, so R and the weekly renewal kernel g share parameters, not an
#' identical discretisation.
#'
#' Each returned value is one posterior sample's MEAN daily R over
#' [week_start - 7*(window_weeks - 1), min(week_start + 6, issue_date)]: EpiNow2's estimates
#' on days with data, never its forecast tail. With the default window_weeks = 1 that is the
#' single week [week_start, end_date], which is what every short-term caller wants.
#'
#' WHY A MULTI-WEEK WINDOW IS AVERAGED INSIDE ONE FIT, not across several. The cascade anchors
#' on a multi-week average of R rather than the final week alone (the final week is the least
#' constrained point of the curve). That average has to be taken WITHIN a posterior sample:
#' sample s's window mean is a draw from the posterior of the window mean, and the draws then
#' carry the right correlation across weeks. Fitting each week separately and combining them
#' afterwards cannot reproduce this — sample indices from different fits share no posterior, so
#' averaging draw s across them pairs unrelated numbers, and concatenating them yields a
#' MIXTURE over weeks (too wide, and centred on the wrong thing) rather than the average. One
#' fit spanning the window is therefore not an optimisation, it is the only correct construction.
#' It is also free here: the fit already runs to `end_date`, so the earlier weeks of the window
#' are estimated WITH the subsequent data rather than as right-edge real-time estimates.
#'
#' AN R IS ALWAYS RETURNED, so every forecast carries its h >= 2 horizon:
#'   1. EpiNow2 is fitted whenever at least one confirmed case is observable. There is no
#'      minimum-data threshold: with sparse data the posterior is prior-dominated, which is the
#'      Bayesian answer, rather than an arbitrary cut-off.
#'   2. A failed fit (an EpiNow2/Stan error, or no usable R samples) is retried once with the
#'      next seed.
#'   3. After two failed fits, the most recent cached posterior for an EARLIER week (at most
#'      RT_CARRY_MAX_WEEKS back) with the same generation time, the same averaging window and
#'      the same model/data specification is carried forward (.rt_recent_posterior()).
#'   4. With no observable case, or when nothing can be carried forward, the draws come from
#'      the model's own PRIOR for R, EPINOW2_R_PRIOR (00_config.R; the rt_opts() prior of
#'      .epinow2_rt_fit(), which EpiNow2 applies to the initial reproduction number): with no
#'      information the model's R is its prior. Each prior draw is held constant over the week.
#' The source is attached as attr(, "source") ("epinow2", "epinow2-retry", "carried-forward"
#' or "prior"); every
#' freshly computed R is appended to <cache_dir>/rt_draws_source_log.csv; a prior fallback also
#' warns. Only EpiNow2 posteriors are cached, so a fallback is re-attempted on the next call. A
#' missing line list, a GT pmf without its generating parameters, or EpiNow2 not being installed
#' are configuration errors and stop.
#'
#' Cached per (case series, generation time, week, delay, ascertainment, R prior, seed, model
#' version) and, for a multi-week window, the window length. A lock directory serialises
#' concurrent callers — parallel LFO workers ask for the same fold's R — so one fits and the
#' others wait for its cache file.
#'
#' CACHE COMPATIBILITY, deliberate. A window_weeks = 1 request hashes and names its cache file
#' exactly as before, so the existing entries stay valid and no refit is forced by this change.
#' A window of 2+ gets its own file-name prefix AND an extra hash component, so a window average
#' can never be served from (or mistaken for) a single-week entry: the two are different
#' estimands computed from the same fit.
#'
#' @param linelist line list (dat$ll).
#' @param gt_pmf daily GT pmf from make_gt_pmf() (carries gt_mean / gt_sd).
#' @param week_start Date; the last observed week of the forecast's training counts.
#' @param issue_date Date; the forecast moment (cases observed after it are invisible).
#' @param window_weeks integer >= 1; how many weeks, ending at week_start, to average R over
#'   within each posterior sample. 1 (the default) is the single-week behaviour.
#' @param lock_stale_min minutes after which a lock whose owner cannot be identified is broken.
#' @param cache_dir directory holding the cache, the locks and the source log.
#' @param n_prior_draws Monte Carlo size of a prior fallback.
#' @return numeric vector of draws of the weekly-mean R, with attr "source".
bayes_rt_week_draws <- function(linelist, gt_pmf, week_start, issue_date,
                                window_weeks = get0("RT_WINDOW_WEEKS", ifnotfound = 1L),
                                lock_stale_min = 180,
                                cache_dir = file.path(OUT_DIAGNOSTICS, "rt_draws"),
                                n_prior_draws = 4000L) {
  if (is.null(linelist) || !is.data.frame(linelist))
    stop("[rt-draws] a line list is required: horizons >= 2 project sources with EpiNow2 R.",
         call. = FALSE)
  week_start <- as.Date(week_start); issue_date <- as.Date(issue_date)
  stopifnot(length(week_start) == 1L, !is.na(week_start),
            length(issue_date) == 1L, !is.na(issue_date))
  window_weeks <- suppressWarnings(as.integer(window_weeks))
  if (length(window_weeks) != 1L || is.na(window_weeks) || window_weeks < 1L)
    stop("[rt-draws] window_weeks must be a single integer >= 1.", call. = FALSE)
  # The first day of the averaging window. week_start is the LAST week of the window, so a
  # 3-week window reaches back 14 days from it and ends at end_date below.
  window_start <- week_start - 7L * (window_weeks - 1L)
  if (!isTRUE(get0(".HAVE_EPINOW2", ifnotfound = FALSE)))
    stop("[rt-draws] EpiNow2 is not installed; it is required for the h >= 2 projection.",
         call. = FALSE)
  gt_mean <- attr(gt_pmf, "gt_mean"); gt_sd <- attr(gt_pmf, "gt_sd")
  if (!is.numeric(gt_mean) || !is.numeric(gt_sd))
    stop("[rt-draws] the GT pmf carries no gt_mean/gt_sd, so EpiNow2 cannot be given the ",
         "same generation time as the invasion model.", call. = FALSE)
  gt_max   <- length(gt_pmf)
  end_date <- min(week_start + 6L, issue_date)
  # Prefixed onto every source-log note so the window is auditable after a run. It goes in the
  # NOTE rather than a new column on purpose: the log is appended to with col.names = FALSE, so
  # widening the schema would write rows with more fields than the existing header declares.
  .wtag <- if (window_weeks > 1L)
    sprintf("[%d-week window %s..%s] ", window_weeks, format(window_start), format(end_date))
  else ""

  cases   <- epinow2_daily_confirmed(linelist, end_date, issue_date = issue_date,
                                     caller = "rt-draws")
  n_cases <- sum(cases$confirm)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (!n_cases)
    return(.rt_prior_fallback(n_prior_draws, week_start, end_date, issue_date,
                              c(gt_mean, gt_sd, gt_max), 0L, cache_dir,
                              paste0(.wtag, "no confirmed case observable by the issue date"),
                              window_weeks, window_start))

  # FINGERPRINT THE TRUNCATION THAT WILL ACTUALLY BE USED, not the onset->sample delay.
  # This fit's right-truncation is no longer derived from effective_onset_sample_delay(): it is
  # ESTIMATED per regime (02_epi_params.R section 4b). Keying on the old delay would let two
  # fits with DIFFERENT truncations collide on one cache entry — and because the as-of estimate
  # is measured from the line list, that can change without any config or code change. Resolving
  # it here is cheap: the as-of regime is a sub-second cohort fit with no Stan, and the extract
  # regime's Stan fit is cached to disk, so this is arithmetic or a file read.
  .trunc_dist <- tryCatch(epinow2_truncation(issue_date = issue_date, ll = linelist),
                          error = function(e) NULL)
  # Hash the OBJECT, not its str() — str() prints 3 significant digits and follows the
  # session's strOptions(), so it both collides on materially different fits and changes
  # between environments. Same key construction as estimate_rt_epinow2(), where the reasoning
  # is written out in full.
  .trunc_k <- if (is.null(.trunc_dist)) "none" else substr(rlang::hash(.trunc_dist), 1L, 12L)
  .hkey <- list(cases = cases, gt = c(gt_mean, gt_sd, gt_max),
                week = week_start, end = end_date, trunc = .trunc_k,
                prior = unlist(EPINOW2_R_PRIOR),
                seed = get0("RANDOM_SEED", ifnotfound = 20260704L),
                version = get0("RT_CACHE_VERSION", ifnotfound = 2L))
  # The model/data SPECIFICATION of this fit, without the data or the week: what a posterior
  # carried forward from an earlier week must share with this one (.rt_recent_posterior()).
  # Stamped on every cached posterior. RT_CACHE_VERSION is part of it, so entries fitted under the
  # retired sample-date dating (version <= 4, and unstamped) can never be carried forward.
  # "asc" is NOT an element of .hkey — ascertainment is no longer a model input at all — so
  # selecting it appended a constant NULL to the hashed list. Harmless to the hash (it is
  # deterministic), but the docstring claimed the carry-forward spec pinned an "ascertainment
  # prior", which it neither did nor could. Naming a field that does not exist in a
  # cache-invalidation key is exactly the kind of thing that later gets "fixed" by adding the
  # field back. NOTE: dropping it CHANGES the hash, so posteriors stamped under the old spec
  # are not carried forward — which is correct, they were stamped against a different list.
  rt_spec <- rlang::hash(.hkey[c("gt", "trunc", "prior", "seed", "version")])
  # Added ONLY for a multi-week window, so a single-week key is byte-identical to the keys the
  # existing cache entries were written under and none of them is invalidated by this change.
  if (window_weeks > 1L) .hkey$window <- window_weeks
  key  <- rlang::hash(.hkey)
  .pfx <- if (window_weeks > 1L) sprintf("epinow2_rtwin%dw_", window_weeks) else "epinow2_rtweek_"
  cache <- file.path(cache_dir,
                     sprintf("%s%s_%s.rds", .pfx, format(week_start, "%Y%m%d"), key))
  lock <- paste0(cache, ".lock")
  waited <- FALSE
  repeat {
    if (file.exists(cache)) return(readRDS(cache))
    if (dir.create(lock, showWarnings = FALSE)) {
      # Record the HOST as well as the pid: a pid from another machine (or from before a
      # reboot) means nothing here, and trusting it is how a lock becomes immortal.
      writeLines(c(as.character(Sys.getpid()), Sys.info()[["nodename"]]),
                 file.path(lock, "pid"))
      break
    }
    # Another process holds the lock. Break it only when that process no longer exists (a
    # killed run never reaches on.exit) or, if its pid cannot be read yet, when the lock is
    # older than lock_stale_min. Breaking is an atomic rename, so two waiters can never both
    # break — and then both fit under — the same lock.
    .pf  <- tryCatch(readLines(file.path(lock, "pid"), n = 2L), error = function(e) character(0))
    pid  <- suppressWarnings(as.integer(if (length(.pf) >= 1L) .pf[1] else NA_character_))
    host <- if (length(.pf) >= 2L) .pf[2] else NA_character_
    age <- as.numeric(difftime(Sys.time(), file.info(lock)$mtime, units = "mins"))
    # STALENESS IS A HARD CEILING, checked FIRST and regardless of the pid.
    #
    # Previously `age > lock_stale_min` was only consulted when the pid file could NOT be read,
    # so a readable pid file made the age bound unreachable and the ONLY exits from this repeat
    # were the cache appearing or that pid dying. A SIGKILLed run never reaches its on.exit, so
    # its lock survives; after a reboot or enough pid churn that integer belongs to some
    # unrelated live process, pskill(pid, 0) succeeds forever, and because run_invasion_lfo()
    # fans out over models within a fold, every Bayesian worker blocks on the same key and the
    # whole LFO deadlocks with no error. A pid from a DIFFERENT host is likewise meaningless.
    .same_host <- is.na(host) || identical(host, Sys.info()[["nodename"]])
    .alive <- .same_host && length(pid) == 1L && !is.na(pid) &&
              isTRUE(tools::pskill(pid, signal = 0L))
    dead <- (is.finite(age) && age > lock_stale_min) || !.alive
    if (isTRUE(dead)) {
      gone <- paste0(lock, ".stale.", Sys.getpid())
      if (suppressWarnings(file.rename(lock, gone))) unlink(gone, recursive = TRUE)
      next
    }
    if (!waited) {
      message(sprintf("[rt-draws] waiting for process %s, which is fitting %s",
                      if (length(pid) == 1L && !is.na(pid)) pid else "?", basename(cache)))
      waited <- TRUE
    }
    Sys.sleep(5)
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (file.exists(cache)) return(readRDS(cache))

  message(sprintf("[rt-draws] EpiNow2 R for %s (issue %s; GT %.2f/%.2f d; %d cases)",
                  if (window_weeks > 1L)
                    sprintf("the %d-week window %s..%s", window_weeks,
                            format(window_start), format(end_date))
                  else sprintf("week %s", format(week_start)),
                  format(issue_date), gt_mean, gt_sd, n_cases))
  # EpiNow2 draws its MCMC initial values from R's global RNG. Seed each attempt so the fit does
  # not depend on whatever ran before, and restore the caller's stream afterwards so a cache miss
  # and a cache hit leave every downstream random number identical.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  seed0 <- get0("RANDOM_SEED", ifnotfound = 20260704L)
  attempt <- function(seed) {
    set.seed(seed)
    # THE FOLD REGIME. issue_date is this fold's own origin, so epinow2_truncation() returns
    # the onset->SAMPLE truncation — correct for a series reaggregate_asof() censored on
    # linelist_observation_date(). Handing it the deployed (onset->appearance) truncation would
    # over-correct the folds and break `delta`, which is fitted here and applied to the live
    # forecast. `ll` is required: the as-of estimate is measured from the line list, not the
    # archive.
    fit <- .epinow2_rt_fit(cases, gt_mean, gt_sd, gt_max, truncation = .trunc_dist)
    smp <- as.data.frame(EpiNow2::get_samples(fit))
    # `!(type %in% "forecast")`, NOT `type != "forecast"`. This file's sibling documents the
    # hazard at 20_forecast_detail.R: base-R logical indexing with an NA returns an all-NA
    # PHANTOM ROW rather than dropping it. The one place it actually matters used the
    # forbidden form. If get_samples() ever stops emitting `type`, `!=` yields logical(0),
    # zero rows survive, both fit attempts "fail", and the run silently takes the R0 prior.
    smp <- smp[smp$variable == "R" & !is.na(smp$date) & !(smp$type %in% "forecast") &
                 smp$date >= window_start & smp$date <= end_date, c("sample", "value", "date")]
    if (!nrow(smp)) stop("EpiNow2 returned no R samples for the window")
    # Per-sample mean over the window's days: a draw from the posterior of the window mean.
    wk <- stats::aggregate(value ~ sample, data = smp, FUN = mean)
    d  <- wk$value[order(wk$sample)]
    if (!length(d) || any(!is.finite(d)) || any(d < 0))
      stop("EpiNow2 R draws are not finite and non-negative")
    # How many distinct days actually entered the average. A window reaching back before the
    # start of the fitted series would silently average over fewer days than asked for, and the
    # result would be labelled a 3-week mean while being something else.
    nd_have <- length(unique(smp$date)); nd_want <- as.integer(end_date - window_start) + 1L
    if (nd_have < nd_want)
      warning(sprintf(paste0("[rt-draws] the %d-week window %s..%s covers %d day(s) but only %d ",
                             "carry R estimates; the returned average is over those %d days."),
                      window_weeks, format(window_start), format(end_date), nd_want, nd_have,
                      nd_have), call. = FALSE)
    attr(d, "window_days") <- nd_have
    d
  }
  errs <- character(0); draws <- NULL; rt_source <- NA_character_
  for (i in 1:2) {
    draws <- tryCatch(attempt(seed0 + i - 1L),
                      error = function(e) { errs <<- c(errs, conditionMessage(e)); NULL })
    if (!is.null(draws)) { rt_source <- if (i == 1L) "epinow2" else "epinow2-retry"; break }
  }
  if (is.null(draws)) {
    # Before the prior: carry the last successful posterior for this GT forward. The prior is a
    # prior on the INITIAL reproduction number (R0-scale, correct at the start of the series);
    # using it for a RECENT week would centre that week near 2.0 where the data say ~1.0-1.2.
    # The prior remains the answer when no fit has EVER succeeded — there, the model's R really
    # is its prior.
    # The prefix keeps the carry-forward WITHIN the same estimand: a single-week posterior and a
    # 3-week window average are different quantities computed from the same fit, and carrying one
    # into the other's slot would silently relabel it.
    carried <- .rt_recent_posterior(cache_dir, week_start, c(gt_mean, gt_sd, gt_max),
                                    prefix = .pfx, spec = rt_spec,
                                    max_issue_date = issue_date)
    if (!is.null(carried)) {
      note <- sprintf(paste0("%sEpiNow2 failed twice (%s); carried forward the posterior for week ",
                             "%s (%d week(s) earlier, same GT and same window), which assumes R ",
                             "has not moved since."),
                      .wtag, paste(errs, collapse = " | "),
                      format(attr(carried, "carried_from")),
                      attr(carried, "carried_gap_weeks"))
      warning(sprintf("[rt-draws] week %s (issue %s): %s", format(week_start),
                      format(issue_date), note), call. = FALSE)
      .rt_log_source(cache_dir, week_start, end_date, issue_date,
                     c(gt_mean, gt_sd, gt_max), n_cases, "carried-forward", note)
      # Re-stamp to THIS week so a consumer reading week_start/end_date is not misled about
      # which week the value is being used for; carried_from keeps the provenance.
      attr(carried, "week_start") <- week_start
      attr(carried, "end_date")   <- end_date
      # NOT cached: a carried value is a stand-in, so the next call retries the real fit.
      return(carried)
    }
    return(.rt_prior_fallback(n_prior_draws, week_start, end_date, issue_date,
                              c(gt_mean, gt_sd, gt_max), n_cases, cache_dir,
                              paste0(.wtag,
                                     paste("EpiNow2 failed twice and no earlier posterior exists",
                                           "for this generation time and window:",
                                           paste(errs, collapse = " | "))),
                              window_weeks, window_start))
  }
  attr(draws, "week_start")   <- week_start
  attr(draws, "end_date")     <- end_date
  attr(draws, "gt")           <- c(mean = gt_mean, sd = gt_sd, max = gt_max)
  attr(draws, "source")       <- rt_source
  attr(draws, "window_weeks") <- window_weeks
  attr(draws, "window_start") <- window_start
  attr(draws, "rt_spec")      <- rt_spec
  # ISSUE DATE, stamped so the carry-forward path can enforce causality. One week_start can
  # legitimately hold several cached entries fitted at DIFFERENT issue dates (a rolling-origin
  # fold vs the production run), and without this attribute .rt_recent_posterior() had no way
  # to tell them apart: an LFO fold whose own fit failed could carry forward a posterior that
  # was fitted at the production ANALYSIS_DATE, i.e. from data observed after that fold's
  # cutoff. The ranking rule ("nearest week, then most end_date") cannot separate them, since
  # both have end_date = week + 6.
  attr(draws, "issue_date")   <- issue_date
  .rt_log_source(cache_dir, week_start, end_date, issue_date, c(gt_mean, gt_sd, gt_max),
                 n_cases, rt_source,
                 paste0(.wtag, if (length(errs)) errs[1] else ""))
  tmp <- tempfile(tmpdir = cache_dir, fileext = ".rds")
  saveRDS(draws, tmp)
  file.rename(tmp, cache)
  draws
}

#' Draws of R from the EpiNow2 model's own prior, EPINOW2_R_PRIOR: a LogNormal with the given
#' NATURAL-scale mean and SD — the parameterisation of EpiNow2::LogNormal(mean, sd), so
#' sdlog = sqrt(log(1 + (sd/mean)^2)) and meanlog = log(mean) - sdlog^2 / 2. Seeded; the
#' caller's RNG stream is restored.
#' Most recent cached EpiNow2 posterior for an EARLIER week at the SAME generation time.
#'
#' The last resort before the prior. EPINOW2_R_PRIOR is a prior on the INITIAL reproduction
#' number — EpiNow2's rt_opts() applies it at the first day of the fitted series, which starts
#' at OUTBREAK_START, so it is an R0-scale quantity and correct in that role. Using it as a
#' RECENT week's R is a different claim: it centres that week on ~2.0 where the data say ~1.0-1.2.
#' Carrying the last successful posterior forward instead keeps the value in the regime the data
#' support, and keeps its width, at the cost of assuming R has not moved since that week — which
#' is stated, logged and bounded by RT_CARRY_MAX_WEEKS (00_config.R) rather than assumed silently.
#'
#' MATCHING IS ON THE GENERATION TIME, not just the week. The GT-marginalised predictor fits one
#' R per GT grid point, so a single week can hold several posteriors that are NOT interchangeable;
#' carrying one across GT variants would silently mix generation times.
#'
#' THE SPECIFICATION MUST MATCH. A candidate must carry attr "rt_spec" identical to `spec` (the
#' generation time, delay, ascertainment prior, R prior, seed and RT_CACHE_VERSION of the calling
#' fit). Entries without that stamp — every posterior written before it existed, including all
#' fitted under the retired sample-date dating, which biased R upward — are never carried forward.
#'
#' A cached file is a posterior BY CONSTRUCTION — prior fallbacks are never written to the cache —
#' and an explicit "prior" source is rejected on top of the specification check.
#'
#' Ranking is deterministic — nearest week, then most data (`end_date`), then file name — so two
#' runs on the same cache carry the same posterior forward rather than whatever the filesystem
#' happens to list first.
#'
#' @return draws with attrs (week_start, end_date, gt, source = "carried-forward",
#'   carried_from, carried_gap_weeks), or NULL when no earlier same-GT posterior exists.
.rt_recent_posterior <- function(cache_dir, week_start, gt,
                                 max_back_weeks = get0("RT_CARRY_MAX_WEEKS", ifnotfound = 8L),
                                 prefix = "epinow2_rtweek_", spec,
                                 max_issue_date = NULL) {
  if (missing(spec) || !is.character(spec) || length(spec) != 1L || is.na(spec))
    stop("[rt-draws] .rt_recent_posterior() needs the calling fit's specification (spec).",
         call. = FALSE)
  if (!dir.exists(cache_dir)) return(NULL)
  # Matching on the PREFIX as well as the generation time keeps the estimand fixed: entries for a
  # different averaging window live under a different prefix and are not interchangeable with
  # these.
  #
  # Matched with startsWith()/substr(), NOT by interpolating the prefix into a regex. The prefix
  # is data, and escaping data into a pattern is both unnecessary here and easy to get wrong —
  # the first version of this line built a character class containing "{}", which R's default
  # (TRE) engine rejects outright as "Invalid contents of {}". startsWith() has no such failure
  # mode. The trailing check that what follows the prefix is 8 digits and an underscore also
  # stops one window prefix from matching another's files (e.g. ...win3w_ against ...win30w_).
  bn <- list.files(cache_dir, pattern = "\\.rds$")
  np <- nchar(prefix)
  sel <- startsWith(bn, prefix) & grepl("^[0-9]{8}_", substring(bn, np + 1L))
  if (!any(sel)) return(NULL)
  bn <- bn[sel]
  fs <- file.path(cache_dir, bn)
  wk <- suppressWarnings(as.Date(substr(bn, np + 1L, np + 8L), format = "%Y%m%d"))
  keep <- !is.na(wk) & wk < week_start & wk >= week_start - 7L * as.integer(max_back_weeks)
  if (!any(keep)) return(NULL)
  fs <- fs[keep]; wk <- wk[keep]
  # Read each candidate ONCE: the end_date used for ranking lives on the object, so a second
  # readRDS at selection time would both double the I/O and open a window for the file to change
  # between the two reads.
  obj <- lapply(fs, function(f) tryCatch(readRDS(f), error = function(e) NULL))
  ed  <- as.Date(vapply(obj, function(d) {
    e <- if (is.null(d)) NULL else attr(d, "end_date")
    if (is.null(e)) NA_real_ else as.numeric(as.Date(e))
  }, numeric(1)), origin = "1970-01-01")
  # Nearest week first; within a week, the entry built from the most data; then file name, so two
  # runs over the same cache carry the same posterior rather than whatever the filesystem lists
  # first. order() puts NA end_dates last within their week.
  ord <- order(-as.numeric(wk), -as.numeric(ed), basename(fs))
  for (i in ord) {
    d <- obj[[i]]
    if (is.null(d) || !length(d) || any(!is.finite(d)) || any(d < 0)) next
    if (identical(attr(d, "source"), "prior")) next   # never carry a prior draw forward
    if (!identical(attr(d, "rt_spec"), spec)) next     # different or unstamped specification
    # CAUSALITY. Never carry forward a posterior fitted at a LATER issue date than the caller's:
    # that fit saw data the caller's fold has not observed, and using it would leak the future
    # into a backtest. An UNSTAMPED entry predates issue-date stamping and cannot be shown to be
    # causal, so it is refused too — the carry-forward path is a rare stand-in (0 of 193 entries
    # in the current source log), so refusing an unverifiable one costs nothing and the fallback
    # (the prior) is explicit rather than silently wrong.
    if (!is.null(max_issue_date)) {
      .idt <- suppressWarnings(as.Date(attr(d, "issue_date")))
      if (length(.idt) != 1L || is.na(.idt) || .idt > as.Date(max_issue_date)) next
    }
    g <- attr(d, "gt")
    if (is.null(g) || length(g) != length(gt)) next
    if (!isTRUE(all.equal(unname(as.numeric(g)), unname(as.numeric(gt)), tolerance = 1e-8))) next
    attr(d, "source")            <- "carried-forward"
    attr(d, "carried_from")      <- wk[i]
    attr(d, "carried_gap_weeks") <- as.integer(round(as.numeric(week_start - wk[i]) / 7))
    return(d)
  }
  NULL
}

.rt_prior_draws <- function(n, seed = get0("RANDOM_SEED", ifnotfound = 20260704L)) {
  pr <- EPINOW2_R_PRIOR
  stopifnot(is.numeric(pr$mean), is.numeric(pr$sd), pr$mean > 0, pr$sd > 0, n >= 1)
  sdlog   <- sqrt(log1p((pr$sd / pr$mean)^2))
  meanlog <- log(pr$mean) - sdlog^2 / 2
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  stats::rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}

#' Prior-draw R for bayes_rt_week_draws() when no EpiNow2 posterior can be obtained: warns,
#' logs the reason, and returns the draws with the same attributes as a fitted R.
.rt_prior_fallback <- function(n, week_start, end_date, issue_date, gt, n_cases, cache_dir,
                               reason, window_weeks = 1L, window_start = week_start) {
  # SEED PER GRID POINT, not once per session. .rt_prior_draws() defaulted to the constant
  # RANDOM_SEED, so if the fallback fired for more than one generation-time grid point in the
  # GT-marginalised predictive, every grid point received the IDENTICAL n draws: the
  # marginalisation over generation time would silently collapse to a single value while
  # continuing to report 15 points. Deriving the seed from the GT and the week keeps the
  # fallback reproducible while making the grid points genuinely distinct.
  .fb_seed <- {
    base <- as.numeric(get0("RANDOM_SEED", ifnotfound = 20260704L))
    k    <- paste(c(as.numeric(gt), as.numeric(week_start)), collapse = "|")
    h <- 0
    for (cc in utils::head(as.integer(charToRaw(k)), 256L)) h <- (h * 31 + cc) %% 2147483647
    as.integer((base + h) %% 2147483647)
  }
  draws <- .rt_prior_draws(n, seed = .fb_seed)
  attr(draws, "week_start") <- week_start
  attr(draws, "end_date")   <- end_date
  attr(draws, "gt")         <- c(mean = gt[1], sd = gt[2], max = gt[3])
  attr(draws, "source")     <- "prior"
  # Stamped here too so all three sources — a fit, a carried-forward posterior and this prior —
  # return the SAME attribute contract, and a consumer never has to test for their presence.
  # A prior draw is constant over the window, so the window average of it is itself.
  attr(draws, "window_weeks") <- as.integer(window_weeks)
  attr(draws, "window_start") <- window_start
  warning(sprintf(paste0("[rt-draws] week %s (issue %s, GT %.2f/%.2f d): %s. R is drawn from the ",
                         "EpiNow2 model's prior, LogNormal(mean %.2f, SD %.2f)."),
                  format(week_start), format(issue_date), gt[1], gt[2], reason,
                  EPINOW2_R_PRIOR$mean, EPINOW2_R_PRIOR$sd), call. = FALSE)
  .rt_log_source(cache_dir, week_start, end_date, issue_date, gt, n_cases, "prior", reason)
  draws
}

#' Append one row to <cache_dir>/rt_draws_source_log.csv recording where a freshly computed R
#' came from, so every prior fallback or retried fit is traceable after a run. Never fatal.
.rt_log_source <- function(cache_dir, week_start, end_date, issue_date, gt, n_cases, source,
                           note) {
  f <- file.path(cache_dir, "rt_draws_source_log.csv")
  row <- data.frame(logged_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                    week_start = format(week_start), end_date = format(end_date),
                    issue_date = format(issue_date), gt_mean = gt[1], gt_sd = gt[2],
                    gt_max = gt[3], n_cases = n_cases, source = source,
                    note = gsub("[\r\n]+", " ", note), stringsAsFactors = FALSE)
  try(utils::write.table(row, f, sep = ",", row.names = FALSE, qmethod = "double",
                         col.names = !file.exists(f), append = file.exists(f)), silent = TRUE)
  invisible(NULL)
}

#' Pair EpiNow2 R draws with S invasion-model draws by seeded random resampling (without
#' replacement when S <= the number of R draws, with replacement otherwise). This samples the
#' PRODUCT of the two separately fitted posteriors, i.e. treats R and the invasion-model
#' parameters as a posteriori independent — a two-stage ("cut") approximation, since both are
#' fitted to the same line list. The caller's RNG stream is restored.
.pair_rt_draws <- function(rt, S, seed) {
  rt <- as.numeric(rt)
  if (!length(rt) || any(!is.finite(rt)) || any(rt < 0))
    stop("[bayes draws] horizons >= 2 need finite, non-negative EpiNow2 posterior draws of R ",
         "(bayes_rt_week_draws()).", call. = FALSE)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  rt[sample.int(length(rt), S, replace = S > length(rt))]
}

#' Per-week invasion hazard mu for every (posterior draw, forecast row).
#'
#' h = 1: mu = g^-1(eta0 + log Lambda(T+1)), with Lambda(T+1) fixed by observed incidence
#' (the stored offset). h >= 2: the import force depends on incidence in the projected
#' weeks, so it is rebuilt INSIDE each draw s:
#'   Y_j^(s)(T+m)        = R^(s) * sum_k g(k) Y_j^(s)(T+m-k)  +  mu_j^(s)(T+m),  m < h,
#'   Lambda_i^(s)(T+h)   = sum_j W[j,i] * sum_k g(k) Y_j^(s)(T+h-k),
#' where R^(s) is an EpiNow2 posterior draw and mu_j^(s)(T+m) is draw s's OWN expected
#' introductions into zone j that week (covariates and link included). The only thing that
#' converts import force into new cases is the fitted model.
#'
#' @param eta0 draws x rows linear predictor WITHOUT the offset.
#' @param nd forecast rows (health_zone, horizon, .off = stored log Lambda for h = 1).
#' @param proj attr "projection" from bayes_forecast_offsets().
#' @param link observation link of the fit.
#' @param hmax largest horizon to fill (rows beyond it stay NA).
#' @return draws x rows matrix of per-week hazards.
.bayes_project_mu <- function(eta0, nd, proj, link, hmax,
                              seed = get0("RANDOM_SEED", ifnotfound = 20260704L)) {
  # Reconcile observation links on the per-week HAZARD scale (see predict_bayes_invasion):
  #   cloglog: mu = exp(eta) (exact p = 1 - exp(-mu));  logit/probit: mu = -log(1 - link^-1(eta)).
  inv <- function(eta) if (identical(link, "logit")) -log1p(-stats::plogis(eta))
                       else if (identical(link, "probit")) -log1p(-stats::pnorm(eta))
                       else exp(eta)                           # cloglog (default)
  S  <- nrow(eta0)
  mu <- matrix(NA_real_, S, ncol(eta0))
  h1 <- which(nd$horizon == 1L)
  mu[, h1] <- inv(sweep(eta0[, h1, drop = FALSE], 2, nd$.off[h1], "+"))
  if (hmax < 2L) return(mu)
  if (is.null(proj))
    stop("[bayes draws] horizons >= 2 need the source projection attached by ",
         "bayes_forecast_offsets().", call. = FALSE)
  Z <- proj$zones; W <- proj$W; G <- proj$G
  # Weekly lag-0 mass of THIS projection's kernel; 0 (identity) when the kernel does not
  # carry the attribute, so a hand-built pmf degrades to the previous behaviour rather than
  # silently taking another profile's correction.
  .g0_w   <- { a <- attr(G, "lag0_mass"); if (is.null(a) || !is.finite(a)) 0 else as.numeric(a) }
  .R_eff_w <- function(R) if (.g0_w > 0) weekly_renewal_R_eff(R, .g0_w) else R
  Yobs <- proj$Y[Z, , drop = FALSE]; Yobs[is.na(Yobs)] <- 0
  T0 <- ncol(Yobs)
  R  <- .pair_rt_draws(proj$rt_draws, S, seed)
  cols_for <- function(h) {
    r  <- which(nd$horizon == h)
    ix <- r[match(Z, nd$health_zone[r])]
    if (anyNA(ix))
      stop(sprintf("[bayes draws] forecast rows for horizon %d do not cover every zone.", h),
           call. = FALSE)
    ix
  }
  Yproj <- vector("list", hmax - 1L)   # Yproj[[m]]: draws x zones incidence in week T0 + m
  # GT-weighted past incidence per draw (== compute_foi's Y_weighted; observed weeks are
  # shared by every draw, projected weeks are draw-specific).
  gweighted <- function(t_for) {
    s <- matrix(0, S, length(Z))
    for (k in seq_along(G)) {
      tp <- t_for - k
      if (tp < 1L) break
      s <- s + G[k] * (if (tp <= T0) rep(Yobs[, tp], each = S) else Yproj[[tp - T0]])
    }
    s
  }
  for (h in seq_len(hmax)) {
    src <- gweighted(T0 + h)
    ch  <- cols_for(h)
    if (h >= 2L) {
      Lam <- src %*% W                 # Lambda_i = sum_j W[j, i] * src_j, per draw
      Lam[Lam < 0] <- 0
      mu[, ch] <- inv(eta0[, ch, drop = FALSE] + log(pmax(Lam, 1e-12)))
    }
    if (h < hmax) {
      # R_eff, NOT R. `G` here is weekly_censored_gt()'s kernel, which drops the weekly lag-0
      # mass and renormalises (5.8% at the medium profile), while `R` is EpiNow2's DAILY-renewal
      # posterior, which retains it. Multiplying a daily-scale R into a lag-0-dropped kernel
      # under-counts the projected week's own contribution and so understates h>=2 growth.
      # weekly_renewal_R_eff() (06_simple_models.R) is the exact algebraic correction; it is the
      # identity at R = 1 and when the kernel carries no lag-0 mass.
      Yh <- .R_eff_w(R) * src + mu[, ch, drop = FALSE]
      # NA contributes nothing rather than poisoning every downstream zone's hazard through
      # src %*% W. But an INFINITE cell must NOT be zeroed: mu comes from inv(eta), and under
      # logit -log1p(-plogis(eta)) is +Inf once eta saturates (~37) while under cloglog exp(eta)
      # overflows at ~709. Blanket-zeroing turned "infinitely infectious" into "not infectious at
      # all", silently UNDER-stating every downstream zone's h>=2 hazard — the opposite of the
      # intended direction. Cap +Inf at the largest finite cell instead, and say so.
      if (any(!is.finite(Yh))) {
        .fin <- Yh[is.finite(Yh)]
        .cap <- if (length(.fin)) max(.fin) else 0
        .n_inf <- sum(is.infinite(Yh))
        if (.n_inf > 0L)
          warning(sprintf("[bayes-project] %d non-finite projected-incidence cell(s) capped at %.3g (hazard overflow).",
                          .n_inf, .cap), call. = FALSE)
        Yh[is.infinite(Yh)] <- .cap
        Yh[is.na(Yh)] <- 0
      }
      Yproj[[h]] <- Yh
    }
  }
  mu
}

#' Per-zone forecast rows for each horizon at the FORECAST week, for .bayes_invasion_draws().
#'
#' h = 1: the import force into week T+1 is fixed by observed incidence, so log Lambda is
#' stored here. h >= 2: the import force depends on incidence in projected weeks, which
#' depends on the posterior itself, so it is NOT stored (logLam = NA); the rows instead carry
#' attr "projection" (observed counts, kernel, weekly GT, EpiNow2 R draws), from which
#' .bayes_project_mu() rolls every source zone forward inside each posterior draw.
#'
#' @param rt_draws EpiNow2 posterior draws of R for the last observed week
#'   (bayes_rt_week_draws()); required when max(horizons) >= 2.
bayes_forecast_offsets <- function(zone_week_nc, mobility_matrices, gt_pmfs,
                                   covariates, osrm_mat, zones_all, mob, gt,
                                   cov_spec, horizons, rt_draws = NULL) {
  if (max(horizons) >= 2L &&
      (is.null(rt_draws) || !length(rt_draws) || any(!is.finite(rt_draws)) || any(rt_draws < 0)))
    stop("[bayes offsets] horizons >= 2 need finite, non-negative EpiNow2 posterior draws of R ",
         "(bayes_rt_week_draws()); none were supplied.", call. = FALSE)
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
  rows <- list()
  # Record EVERY intermediate week up to max(horizons), not just the requested horizons:
  # predict_bayes_invasion forms the cumulative hazard p(<=h)=1-exp(-sum_{h'<=h} mu(h')) by
  # summing the recorded per-week rows with horizon <= h. Recording only {1,2} happens to be
  # correct because they are contiguous; a non-contiguous request (e.g. c(1,3)) would silently
  # drop week 2 and under-count. Recording all weeks makes the cumulation correct for any set
  # (predict still emits only the requested horizons).
  for (h in seq_len(max(horizons))) {
    t_for <- ncol(Y) + h
    Lam <- if (h == 1L) compute_foi(Y, W, G, t_idx = t_for, zones_all)
           else rep(NA_real_, length(zones_all))       # per draw, in .bayes_project_mu()
    X <- .feature_matrix(t_for, Y, A_wide, W, G, static, osrm_mat, zones_all,
                         cov_spec, Y_active = Y, S_wide = S_wide)
    df <- data.frame(health_zone = zones_all, horizon = h,
                     Lambda = as.numeric(Lam),
                     logLam = log(pmax(as.numeric(Lam), 1e-12)))
    if (!is.null(X)) df <- cbind(df, as.data.frame(X[, , drop = FALSE]))
    rows[[length(rows) + 1L]] <- df
  }
  out <- dplyr::bind_rows(rows)
  attr(out, "projection") <- list(Y = Y, W = W[zones_all, zones_all, drop = FALSE], G = G,
                                  zones = zones_all,
                                  rt_draws = if (max(horizons) >= 2L) as.numeric(rt_draws))
  out
}

# --------------------------------------------------------------------------
# Shared draws head + summarise tail for the invasion predictors (review §2.1/§5.1).
# Factored out so predict_bayes_invasion (single GT) and predict_bayes_gt_marginal
# (GT-marginalised mixture, §2.1) summarise IDENTICALLY from posterior draws.
# --------------------------------------------------------------------------

#' Posterior DRAWS of the per-(zone,horizon) cumulative hazard and invasion
#' probability for one fitted model at one generation time. Returns a list keyed
#' by horizon ("h1","h2",...), each a list(zones, horizon, cum, p) of
#' draws x zones matrices. This is the computational head previously inside
#' predict_bayes_invasion (see there for the link-reconciliation rationale).
#' @param linpred_fn optional function(fit, nd) -> draws x rows linear predictor EXCLUDING
#'   the offset, replacing the default brms call below. Supplied only by the rw1
#'   time-varying models, whose forecast week must start from the walk's last fitted level
#'   rather than from the population distribution brms would draw a new level from
#'   (.tv_linpred_fn()). NULL — every other model — takes the default path unchanged.
.bayes_invasion_draws <- function(fit, offsets_df, design, horizons,
                                  delta = NULL, linpred_fn = NULL) {
  if (is.null(fit)) return(NULL)
  feat <- design$feat
  nd <- offsets_df
  for (f in feat) if (f %in% names(nd)) {
    nd[[f]] <- (nd[[f]] - design$center[[f]]) / design$scale[[f]]
    nd[[f]][!is.finite(nd[[f]])] <- 0     # non-finite standardised covariate -> standardised mean
  }
  nd$.off <- nd$logLam; nd$logLam <- 0    # zero model offset; add per-horizon logLam back
  # sample_new_levels = "gaussian" is REQUIRED for the tv-week models: the forecast week is a
  # new level of .week_f, and its u must be drawn from Normal(0, sigma_w) so the predictive
  # interval carries the estimated week-to-week variance. brms' default for an unseen level
  # is "uncertainty", which resamples an EXISTING group's deviation — that would attach one
  # arbitrary past week's level to the forecast rather than integrating over the population.
  # The argument is inert for models without group-level terms, so it is safe to pass always.
  eta0 <- if (is.function(linpred_fn)) linpred_fn(fit, nd) else
    suppressWarnings(brms::posterior_linpred(
      fit, newdata = nd, allow_new_levels = TRUE,
      sample_new_levels = "gaussian"))                      # draws x rows (no offset)
  if (!is.matrix(eta0) || nrow(eta0) < 1L || ncol(eta0) != nrow(nd))
    stop(sprintf(paste0("[bayes] the linear predictor is %s but the forecast grid has %d rows; ",
                        "refusing to summarise a mismatched posterior."),
                 if (is.matrix(eta0)) sprintf("%d x %d", nrow(eta0), ncol(eta0))
                 else paste("not a matrix (", class(eta0)[1], ")"), nrow(nd)),
         call. = FALSE)
  link <- tryCatch(attr(fit, "link") %||% fit$family$link, error = function(e) "cloglog")
  # Per-week hazard (draws x rows): h = 1 from the stored offset, h >= 2 from the per-draw
  # source projection, which needs these same draws (.bayes_project_mu()).
  mu <- .bayes_project_mu(eta0, nd, attr(offsets_df, "projection"), link, max(horizons))
  zones <- unique(nd$health_zone)
  idx <- split(seq_len(nrow(nd)), nd$health_zone)
  lapply(stats::setNames(horizons, paste0("h", horizons)), function(h) {
    cum <- vapply(zones, function(z) {
      cols <- idx[[z]]; cols <- cols[nd$horizon[cols] <= h]
      if (!length(cols)) return(rep(0, nrow(mu)))
      rowSums(mu[, cols, drop = FALSE]) }, numeric(nrow(mu)))     # draws x zones
    # ---- optional post-hoc recalibration (16b_invasion_recalibration.R) -------
    # Applied to the cumulative hazard INSIDE the draws, before any summarisation.
    # That placement is what makes every downstream summary exact rather than
    # approximate: mu_forecast = colMeans(cum) scales by delta exactly, p_invasion
    # stays a genuine posterior MEAN of the recalibrated probability (not the
    # recalibrated posterior mean, which differs by the Jensen gap), and p_median /
    # p_lo / p_hi remain true quantiles. Because delta * cum is monotone WITHIN each
    # draw, the per-draw zone ordering is untouched, so the rank credible interval
    # computed by .invasion_horizon_summary() is bit-identical.
    list(zones = zones, horizon = h, cum = cum, p = 1 - exp(-cum))
  }) |> .scale_invasion_draws(delta)
}

#' Rescale a draws list by a hazard-scale recalibration factor.
#'
#' Multiplies the CUMULATIVE hazard inside the draws and rebuilds the probabilities from
#' it, which is what makes every downstream summary exact rather than approximate:
#' mu_forecast = colMeans(cum) scales by delta exactly, p_invasion stays a genuine
#' posterior MEAN of the recalibrated probability (not the recalibrated posterior mean,
#' which differs by the Jensen gap), and p_median / p_lo / p_hi remain true quantiles.
#' delta * cum is monotone WITHIN each draw, so the per-draw zone ordering is untouched
#' and the rank credible interval is bit-identical.
#'
#' Factored out so ONE code path applies the factor, whether it is applied to a single
#' model's draws or to pooled mixture draws, and so both scales can be produced from the
#' SAME draws without recomputing the posterior linear predictor.
.scale_invasion_draws <- function(drl, delta) {
  if (is.null(drl) || is.null(delta) || !length(delta)) return(drl)
  lapply(drl, function(el) {
    dh <- .bayes_recal_delta(delta, el$horizon)
    if (dh == 1) return(el)
    cum <- el$cum * dh
    el$cum <- cum
    el$p   <- 1 - exp(-cum)
    el
  })
}

#' Summarise a draws list into the prediction tibble (shared tail).
.finalise_invasion_summary <- function(drl, affected_zones = character(0)) {
  if (is.null(drl)) return(NULL)
  dplyr::bind_rows(lapply(drl, .invasion_horizon_summary, affected_zones = affected_zones)) %>%
    # p_case_invasion mirrors p_invasion so EVERY Bayesian prediction (per-model, not
    # only the stacked one) carries the column the risk-score / map products key on.
    dplyr::mutate(p_case_invasion = p_invasion,
                  was_active_before = health_zone %in% affected_zones) %>%
    dplyr::mutate(dplyr::across(
      c(mu_forecast, p_invasion, p_case_invasion, p_median, p_lo, p_hi, p_sd),
      ~ ifelse(was_active_before, NA_real_, .x)))
}

#' Resolve the recalibration factor for one horizon.
#'
#' `delta` may be NULL (no recalibration), a single positive scalar applied to every
#' horizon, or a vector named by horizon ("1", "2", ... or "h1", "h2", ...). A horizon
#' with no entry is NOT recalibrated, and anything unusable degrades to 1 with a
#' warning rather than silently rescaling a forecast by a wrong or missing factor.
.bayes_recal_delta <- function(delta, h) {
  if (is.null(delta) || !length(delta)) return(1)
  d <- if (length(delta) == 1L && is.null(names(delta))) {
    as.numeric(delta)
  } else {
    nms <- names(delta)
    if (is.null(nms)) {
      warning("[recal] unnamed multi-element `delta`; no recalibration applied.", call. = FALSE)
      return(1)
    }
    hit <- match(c(as.character(h), paste0("h", h)), nms)
    hit <- hit[!is.na(hit)]
    if (!length(hit)) return(1)
    as.numeric(delta[[hit[1]]])
  }
  if (length(d) != 1L || !is.finite(d) || d <= 0) {
    warning(sprintf("[recal] unusable delta for horizon %s; no recalibration applied.", h),
            call. = FALSE)
    return(1)
  }
  d
}

#' Summarise one horizon's posterior draws (draws x zones cum/p) into the
#' per-zone tibble: posterior mean + median invasion probability, 90% CrI, sd,
#' infection-scale probability, and the RANK credible interval (review §5.1).
#' Affected zones are excluded from the ranking (left NA); the caller masks their
#' probabilities to NA.
.invasion_horizon_summary <- function(dr, affected_zones = character(0)) {
  zones <- dr$zones; p <- dr$p; cum <- dr$cum
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
    rank_med = rank_med, rank_lo = rank_lo, rank_hi = rank_hi)
}

#' Posterior invasion probability per (zone, horizon) with a 90% credible
#' interval. Cumulative hazard across horizons (p(<=h) = 1 - exp(-sum mu(h')));
#' the model offset is zeroed and log Lambda(h) added back explicitly so the
#' horizon-specific offset is unambiguous. Forecast covariates are standardised
#' with the fit design's center/scale; affected zones are masked to NA.
predict_bayes_invasion <- function(fit, offsets_df, design, horizons,
                                   affected_zones = character(0),
                                   delta = NULL, linpred_fn = NULL) {
  if (is.null(fit)) return(NULL)
  # Draws head + per-horizon summarise tail are factored into shared helpers so the
  # GT-marginalised predictor (predict_bayes_gt_marginal, review §2.1) reuses the
  # IDENTICAL summarisation (mean/median/90% CrI + rank CrI) on pooled mixture draws.
  drl <- .bayes_invasion_draws(fit, offsets_df, design, horizons, delta = delta,
                               linpred_fn = linpred_fn)
  .finalise_invasion_summary(drl, affected_zones)
}

#' BOTH probability scales from ONE posterior-linear-predictor pass.
#'
#' The manuscript reports the deployed (recalibrated) forecast as the primary product and
#' keeps an uncorrected twin alongside it. Predicting twice would double the expensive part
#' — brms::posterior_linpred over draws x rows — for a factor that is a scalar multiple of
#' the cumulative hazard, so the draws are computed once, summarised raw, then rescaled and
#' summarised again. `calibrated` is IDENTICAL to predict_bayes_invasion(delta = delta) and
#' `raw` to predict_bayes_invasion(delta = NULL); with delta NULL the two are the same
#' object, which is the honest answer when no factor applies.
#'
#' @return list(raw, calibrated) of prediction tibbles, or NULL.
predict_bayes_invasion_scales <- function(fit, offsets_df, design, horizons,
                                          affected_zones = character(0),
                                          delta = NULL, linpred_fn = NULL) {
  if (is.null(fit)) return(NULL)
  drl <- .bayes_invasion_draws(fit, offsets_df, design, horizons, delta = NULL,
                               linpred_fn = linpred_fn)
  if (is.null(drl)) return(NULL)
  raw <- .finalise_invasion_summary(drl, affected_zones)
  cal <- if (is.null(delta) || !length(delta)) raw
         else .finalise_invasion_summary(.scale_invasion_draws(drl, delta), affected_zones)
  list(raw = raw, calibrated = cal)
}

#' Pool a set of posterior-draw COMPONENTS into one mixture predictive and
#' summarise it (review §2.1/§2.4). Each component is a `.bayes_invasion_draws`
#' output (named-by-horizon list of list(zones, horizon, cum, p)); the
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
                                     seed = get0("RANDOM_SEED", ifnotfound = 20260704L),
                                     return_draws = FALSE) {
  keep <- !vapply(comp_draws, is.null, logical(1))
  comp_draws <- comp_draws[keep]; weights <- as.numeric(weights)[keep]
  if (!length(comp_draws)) return(NULL)
  w <- weights; w[!is.finite(w) | w < 0] <- 0
  if (sum(w) <= 0) w <- rep(1, length(w))
  w <- w / sum(w)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  # Restore the caller's stream exactly. When there was NO stream before, REMOVE the one
  # set.seed() just created rather than leaving a seeded generator behind — otherwise
  # summarising a mixture silently changes every subsequent random draw in the session.
  # Same convention as bootstrap_invasion_delta() (16b_invasion_recalibration.R).
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  N <- nrow(comp_draws[[1]][[1]]$p)
  if (is.null(total_draws)) total_draws <- N
  m_k   <- pmax(1L, round(w * total_draws))
  hkeys <- names(comp_draws[[1]])
  # ONE resample index per component, REUSED across horizons. It used to be redrawn inside
  # the per-horizon lapply, so the pooled h=1 and h=2 samples came from DIFFERENT posterior
  # draws of the same fit — the mixture was then not a coherent joint posterior sample, and
  # the within-draw cumulative-hazard ordering p(<=1) <= p(<=2) was destroyed. The marginal
  # means are unaffected (each horizon still averages the same underlying draws), but the
  # credible intervals and the rank CrIs at the two horizons came from unmatched draw sets,
  # and a difference taken between the two summaries could come out negative purely from
  # resampling noise. Drawing once keeps every horizon on the same draws.
  ix_k <- lapply(seq_along(comp_draws), function(k) {
    n_k <- nrow(comp_draws[[k]][[hkeys[1]]]$p)
    sample.int(n_k, m_k[k], replace = m_k[k] > n_k)
  })
  pooled <- lapply(hkeys, function(hk) {
    parts <- lapply(seq_along(comp_draws), function(k) {
      dr <- comp_draws[[k]][[hk]]
      ix <- ix_k[[k]]
      list(cum = dr$cum[ix, , drop = FALSE], p = dr$p[ix, , drop = FALSE],
           p = dr$p[ix, , drop = FALSE])
    })
    list(zones   = comp_draws[[1]][[hk]]$zones,
         horizon = comp_draws[[1]][[hk]]$horizon,
         # carried so the pooled draws can be rescaled by a recalibration factor with the
         # SAME helper a single model's draws use (.scale_invasion_draws)
         cum   = do.call(rbind, lapply(parts, `[[`, "cum")),
         p     = do.call(rbind, lapply(parts, `[[`, "p")),
         p = do.call(rbind, lapply(parts, `[[`, "p")))
  })
  names(pooled) <- hkeys
  # `return_draws` hands back the POOLED draws instead of a summary, so a caller that needs
  # more than one probability scale (the deployed forecast and its uncorrected twin) can
  # summarise the same mixture twice rather than re-drawing it. Resampling then rescaling
  # and rescaling then resampling are identical here — the factor is a per-horizon scalar
  # and the resample keeps whole draw rows — so the two routes give the same object.
  if (isTRUE(return_draws))
    return(list(pooled = pooled,
                weights_used = stats::setNames(m_k / sum(m_k), names(comp_draws)),
                weights_nominal = stats::setNames(w, names(comp_draws)),
                draws_per_component = stats::setNames(m_k, names(comp_draws))))
  out <- .finalise_invasion_summary(pooled, affected_zones)
  # REALISED mixture weights, not the nominal ones. m_k = pmax(1, round(w * total_draws))
  # rounds, and forces at least one draw from every component, so the pooled sample's actual
  # weights are m_k / sum(m_k) — a component with a negligible prior weight is over-weighted
  # by the pmax(1, .). Reporting `w` here would have described a mixture that was not the one
  # summarised; both are now recorded.
  attr(out, "weights_used")     <- stats::setNames(m_k / sum(m_k), names(comp_draws))
  attr(out, "weights_nominal")  <- stats::setNames(w, names(comp_draws))
  attr(out, "draws_per_component") <- stats::setNames(m_k, names(comp_draws))
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
                                      gt_prior = GT_PRIOR, linelist = NULL,
                                      analysis_date = ANALYSIS_DATE,
                                      iter = 2000L,
                                      chains = 2L, link = "cloglog", total_draws = NULL,
                                      tv = "none", delta = NULL, scales = FALSE,
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
    # Same tv dispatch as the suite and the LFO closure. Without it the DEPLOYED featured
    # forecast would silently refit a tv-week model as fixed-beta, so the published forecast
    # would come from a different model than the one cross-validation selected.
    fit <- if (.tv_is_structural(tv))
             fit_bayes_renewal_tv(des, cov_spec = cov_spec, iter = iter, chains = chains,
                                  seed = seed, link = link, tv = tv)
           else
             fit_bayes_renewal(des, cov_spec = cov_spec, iter = iter, chains = chains,
                               seed = seed, link = link)
    if (is.null(fit)) return(NULL)
    # R for the h >= 2 source projection, fitted with THIS grid point's generation time
    # (bayes_rt_week_draws() always returns one: EpiNow2, a retried fit, or the model's prior).
    rt <- if (max(horizons) >= 2L)
            bayes_rt_week_draws(linelist, gp$gt_pmfs[[k]],
                                week_start = max(zone_week_nc$week_start),
                                issue_date = analysis_date) else NULL
    off <- tryCatch(bayes_forecast_offsets(zone_week_nc, mobility_matrices, gp$gt_pmfs,
             covariates, osrm_mat, zones_all, mob = mob, gt = k, cov_spec = cov_spec,
             horizons = horizons, rt_draws = rt), error = function(e) NULL)
    if (is.null(off)) return(NULL)
    # .bayes_invasion_draws() reads the grouping column straight off the offsets, so the
    # forecast week's .week_f level must be attached here — the tv predictor wrapper is not
    # in this path (we need the DRAWS, not the summarised tibble, for the GT mixture).
    if (.tv_is_structural(tv)) off <- .tv_offsets(fit, off, des)
    # DEPLOYMENT RECALIBRATION must reach THIS path too. fit_bayes_suite() recalibrates
    # every suite model through predict_bayes_invasion(delta = ...), but the FEATURED
    # forecast — the one that becomes the watch-list, the maps and the daily re-issue — is
    # produced here instead whenever GT_MARGINALISE_FEATURED is on. Without this argument
    # the featured product would be the ONLY deployed forecast left on the raw scale, and
    # silently so (the transform is rank-preserving, so no ordering would look wrong).
    # APPROXIMATION, stated rather than hidden: delta was estimated on the cross-validated
    # FIXED-GT model, and is reused here for the GT-MARGINALISED predictive. Marginalising
    # widens the intervals but barely moves the mean, so the level correction transfers;
    # leaving the featured forecast on the raw scale while the rest of the suite is
    # corrected would be the larger error.
    # WHERE it is applied: components are drawn RAW here (delta = NULL on the next line) and
    # the recalibration factor is applied to the POOLED mixture below. Rescaling before or
    # after pooling is numerically identical — delta scales each draw's cumulative hazard, the
    # factor is a per-horizon scalar, and the resample keeps whole draw rows — and doing it
    # once at the end is what lets both probability scales come out of a single set of fits.
    # (This block used to ALSO claim the factor was "applied per GT grid point, before
    # mixing", which is the opposite of what the call below does.)
    drl <- .bayes_invasion_draws(fit, off, des, horizons, delta = NULL,
                                 linpred_fn = if (.tv_is_structural(tv)) .tv_linpred_fn(fit)
                                              else NULL)
    if (is.null(drl)) return(NULL)
    list(draws = drl, key = k, weight = as.numeric(gp$weights[[k]]))
  })
  # REPORT THE DROPPED GRID POINTS. A per-point failure (a fit that would not start, a NULL
  # offset, no R draws) returned NULL silently, .mix_and_summarise_draws() then renormalised the
  # weights over the SURVIVORS, and the result was stamped with the FULL grid — so a mixture over
  # 12 of 15 points was published as the GT prior. The drops are not missing at random: beta0
  # scales with the GT, so it is the grid's tails — precisely what the marginalisation exists to
  # represent — that fail first. Refuse when too much prior mass is gone rather than quietly
  # reporting a truncated prior as the whole one.
  .ok      <- !vapply(comp, is.null, logical(1))
  .w_all   <- vapply(seq_along(comp), function(i) gp$weights[[i]], numeric(1))
  .mass_ok <- if (sum(.w_all) > 0) sum(.w_all[.ok]) / sum(.w_all) else 0
  if (any(!.ok))
    warning(sprintf("[gt-marginal] %d of %d GT grid point(s) failed (%s); %.1f%% of the prior mass retained.",
                    sum(!.ok), length(.ok),
                    paste(utils::head(names(gp$gt_pmfs)[!.ok], 5), collapse = ", "),
                    100 * .mass_ok), call. = FALSE)
  comp <- comp[.ok]
  if (!length(comp)) { warning("[gt-marginal] no GT grid point fit; returning NULL"); return(NULL) }
  .min_mass <- get0("GT_MARGINAL_MIN_MASS", ifnotfound = 0.9)
  if (.mass_ok < .min_mass) {
    warning(sprintf("[gt-marginal] only %.1f%% of the GT prior mass is representable (< %.0f%%); refusing to report a truncated prior as the full one.",
                    100 * .mass_ok, 100 * .min_mass), call. = FALSE)
    return(NULL)
  }
  names(comp) <- vapply(comp, function(cc) cc$key, character(1))
  # Pool the per-GT posterior draws into the mixture predictive (shared machinery), then
  # summarise it once per probability scale.
  mix <- .mix_and_summarise_draws(
    lapply(comp, `[[`, "draws"), vapply(comp, function(cc) cc$weight, numeric(1)),
    affected_zones = affected_zones, total_draws = total_draws, seed = seed,
    return_draws = TRUE)
  .dress <- function(tb) {
    attr(tb, "weights_used")        <- mix$weights_used
    attr(tb, "weights_nominal")     <- mix$weights_nominal
    attr(tb, "draws_per_component") <- mix$draws_per_component
    attr(tb, "gt_grid")             <- gp$grid
    # The grid ACTUALLY marginalised over, and the prior mass it carries, so a consumer can
    # tell a full marginalisation from a truncated one.
    attr(tb, "gt_grid_used")        <- names(comp)
    attr(tb, "gt_prior_mass_used")  <- .mass_ok
    attr(tb, "gt_weights_used")     <- mix$weights_used
    tb
  }
  raw <- .dress(.finalise_invasion_summary(mix$pooled, affected_zones))
  cal <- if (is.null(delta) || !length(delta)) raw
         else .dress(.finalise_invasion_summary(.scale_invasion_draws(mix$pooled, delta),
                                                affected_zones))
  if (isTRUE(scales)) list(raw = raw, calibrated = cal) else cal
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
#' model's posterior-median force of infection. For horizons >= 2 the sources are rolled
#' forward with the projection .bayes_project_mu() applies per draw, evaluated at point
#' values: local renewal at the posterior MEDIAN of `rt_draws`, plus new introductions
#' beta_med * Lambda (this table's own hazard) — the posterior-median analogue of the forecast.
#'
#' @param rt_draws EpiNow2 posterior draws of R for the last observed week
#'   (bayes_rt_week_draws()); required when max(horizons) >= 2.
#'
#' NOTE: exact for a covariate-free hazard. If the featured model carries covariates,
#' beta is zone-varying and this scalar-beta decomposition drops the per-zone covariate
#' modulation of the DESTINATION hazard — the caller warns and the renewal-share column
#' stays valid, but `foi` should then be read as the intercept-only import force.
#'
#' SCALE. `beta_med` is the RAW posterior-median import coefficient, so `foi` and
#' `dest_hazard_week` are on the model's own hazard scale. When the live forecast is
#' recalibrated for deployment (INVASION_RECALIBRATE_DEPLOY, 16b), the published
#' probabilities carry an additional per-horizon factor delta that is applied to the
#' CUMULATIVE hazard inside the posterior draws, so the per-week increments this function
#' emits are NOT simply delta times these values and the two are deliberately not forced
#' to match. `share_of_dest` — the column this table is used for — is a ratio within a
#' destination week and is invariant to any such global rescaling.
#'
#' @return long tibble: one row per non-negligible directed (origin -> dest) pair
#'   per requested horizon, with the raw import force (beta-free), the force of
#'   infection (beta * import force), the destination totals and the renewal
#'   share of the destination hazard.
bayes_pairwise_import_force <- function(zone_week_nc, mobility_matrices, gt_pmfs,
                                        zones_all, mob, gt, beta_med, horizons,
                                        rt_draws = NULL, affected_zones = character(0),
                                        province_map = NULL, tol = 1e-12) {
  if (is.null(beta_med) || !is.finite(beta_med) || beta_med < 0)
    stop("beta_med must be a finite non-negative scalar (posterior-median import coefficient)")
  W <- mobility_matrices[[mob]]
  if (is.null(W)) stop(sprintf("mobility matrix '%s' not found for pairwise FOI", mob))
  W <- W[zones_all, zones_all, drop = FALSE]
  G <- daily_to_weekly_gt(gt_pmfs[[gt]])
  Y <- .count_wide(zone_week_nc, zones_all, "confirmed_nc")
  Yw <- Y
  hmax <- max(horizons)
  if (hmax >= 2L && (is.null(rt_draws) || !length(rt_draws) || any(!is.finite(rt_draws))))
    stop("[pairwise FOI] horizons >= 2 need EpiNow2 posterior draws of R (bayes_rt_week_draws()).",
         call. = FALSE)
  R_med <- if (hmax >= 2L) stats::median(rt_draws) else NA_real_
  # WEEKLY-EQUIVALENT R, matching .bayes_project_mu(). daily_to_weekly_gt() drops the lag-0
  # mass and renormalises (5.81% at the medium profile), so the raw daily-scale EpiNow2 R
  # under-projects the weekly renewal; weekly_renewal_R_eff() is the algebraic correction.
  # This function's own docstring calls it "the point analogue of .bayes_project_mu()", and
  # that function applies the correction — so without it the two used DIFFERENT renewal
  # kernels at h >= 2. Effect is ~+0.1% at the deployed national R ~= 1.02 and ~+10% at
  # R = 2.5; share_of_dest is a within-destination ratio and partly protected, but foi,
  # import_force, source_origin and dest_hazard_week at h >= 2 are not.
  .g0_pw  <- attr(G, "lag0_mass")
  .g0_pw  <- if (is.numeric(.g0_pw) && length(.g0_pw) == 1L && is.finite(.g0_pw)) .g0_pw else 0
  R_med_w <- if (is.finite(R_med) && .g0_pw > 0) weekly_renewal_R_eff(R_med, .g0_pw) else R_med
  atrisk <- setdiff(zones_all, affected_zones)
  rows <- list()
  for (h in seq_len(hmax)) {
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
    # Project sources forward for the next horizon (point analogue of .bayes_project_mu()):
    # local renewal at the median EpiNow2 R plus this table's own expected introductions.
    if (h < hmax) Yw <- cbind(Yw, R_med_w * src + beta_med * Lambda)
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

#' Posterior trajectory of the import coefficient beta_t over the training weeks.
#'
#' Returns one row per training week with the posterior median and 90% credible interval of
#'     beta_t = exp(beta0 + u_t)
#' — the import->invasion conversion rate the model believes held in that week — plus the
#' forecast weeks when `horizons` is given. For a FIXED-beta model this is a flat line, which
#' is exactly the reference the time-varying variants are read against, so it is computed the
#' same way for every model rather than only for the tv ones.
#'
#' HOW beta_t IS ISOLATED, and why this is not done by reading parameters out of the draws.
#' The five processes store their time effect in four different places (a group-level `r_`
#' array for week/rw1, a latent GP for ar1/gp, a population coefficient for trend, nothing at
#' all for the fixed-beta model), so reading raw parameters would need four extraction paths
#' and would silently return a flat line for any process it did not recognise. Instead the
#' model is asked directly: evaluate its own linear predictor at each week with the offset set
#' to zero and every covariate at its standardised MEAN, which leaves exactly beta0 + u_t.
#' One path, and a process this function has never heard of is still handled correctly.
#'
#' `week_idx` is the one covariate that must NOT be held at its mean: it IS the tv-trend
#' model's time term, so it is set to its standardised value at each week — otherwise the
#' trend model would report a flat trajectory, i.e. the opposite of what it fitted.
#'
#' @param fit a brmsfit from fit_bayes_renewal() or fit_bayes_renewal_tv(), or NULL.
#' @param design the design it was fitted on (for the covariate centring/scaling).
#' @param label model label carried into the output.
#' @param horizons optionally extend the trajectory past the last training week by these many
#'   weeks, using the model's OWN forecast path — so the plotted continuation is the one the
#'   forecast actually uses (rw1 persists, week reverts, ar1/gp decay), not an illustration.
#' @param week_dates optional named/indexed map from week index to a calendar date.
#' @return tibble(model, week, week_date, is_forecast, beta, beta_lo, beta_hi), or NULL.
bayes_beta_trajectory <- function(fit, design, label = "", horizons = NULL,
                                  week_dates = NULL) {
  if (is.null(fit) || is.null(design) || is.null(design$d) || !nrow(design$d)) return(NULL)
  if (!".week" %in% names(design$d)) return(NULL)
  tv    <- attr(fit, "tv") %||% "none"
  weeks <- sort(unique(as.integer(design$d$.week)))
  fut   <- if (is.null(horizons) || !length(horizons)) integer(0) else
    sort(unique(max(weeks) + as.integer(horizons)))
  all_w <- c(weeks, fut)
  nd <- data.frame(.week = all_w, logLam = 0)
  # Every covariate at its standardised mean (0), so only beta0 + u_t is left. week_idx is
  # the exception: standardise the week's own value with the design's centre/scale.
  for (f in design$feat) {
    nd[[f]] <- if (identical(f, "week_idx"))
      (all_w - design$center[[f]]) / design$scale[[f]] else 0
  }
  if (.tv_is_structural(tv)) {
    # Training weeks keep their own level; forecast weeks take the model's forecast path.
    tw <- attr(fit, "train_weeks"); if (is.null(tw)) tw <- weeks
    # THE FACTOR'S LEVELS MATTER, NOT ONLY ITS VALUES. brms validates newdata against the
    # fitted grouping levels and rejects a factor that merely CARRIES an unseen level, even
    # on rows that never use it. For rw1 the forecast rows are pinned to the last training
    # level (the walk is added explicitly below), so the level set is exactly the training
    # one; for the iid-week model the forecast week is genuinely a new level, which is the
    # point, and the prediction call below passes allow_new_levels for it.
    nd$.week_f <- if (.tv_needs_manual_forecast(tv))
      factor(ifelse(all_w %in% tw, all_w, max(tw)), levels = as.character(tw))
    else
      factor(all_w, levels = as.character(sort(unique(c(tw, all_w)))))
    eta <- if (.tv_needs_manual_forecast(tv)) {
      # rw1: population part + the propagated walk, exactly as the forecast path builds it.
      # .tv_linpred_fn() keys the propagation off .week - max_train_week, which is <= 0 on
      # every training row and so contributes u_T there; that is wrong for a TRAINING week,
      # whose own level is fitted. So the two blocks are evaluated separately.
      e_tr <- suppressWarnings(brms::posterior_linpred(
        fit, newdata = nd[all_w %in% tw, , drop = FALSE]))
      if (length(fut)) {
        nd_f <- nd[!(all_w %in% tw), , drop = FALSE]
        e_fu <- .tv_linpred_fn(fit)(fit, nd_f)
        cbind(e_tr, e_fu)
      } else e_tr
    } else {
      suppressWarnings(brms::posterior_linpred(fit, newdata = nd, allow_new_levels = TRUE,
                                               sample_new_levels = "gaussian"))
    }
  } else {
    eta <- suppressWarnings(brms::posterior_linpred(fit, newdata = nd))
  }
  if (!is.matrix(eta) || ncol(eta) != length(all_w)) return(NULL)
  b <- exp(eta)
  out <- tibble::tibble(
    model = label, week = all_w,
    is_forecast = !(all_w %in% weeks),
    beta    = apply(b, 2, stats::median),
    beta_lo = apply(b, 2, stats::quantile, 0.05, names = FALSE),
    beta_hi = apply(b, 2, stats::quantile, 0.95, names = FALSE),
    tv = tv)
  if (!is.null(week_dates)) {
    # unname(): indexing a NAMED lookup carries the names onto the column, so week_date came
    # back as a named vector — harmless in a plot, but it is written to a published CSV and
    # it makes every equality comparison against a plain Date fail.
    out$week_date <- unname(as.Date(week_dates[as.character(out$week)]))
  } else out$week_date <- as.Date(NA)
  out
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
  # sd_/sdgp_/lscale_ are included so a time-varying fit whose ONLY poorly mixed parameter is
  # its process scale or lengthscale cannot report a clean rhat. (grep is on names, so these
  # prefixes simply match nothing in a fixed-beta fit.)
  rmx <- suppressWarnings(max(rh[grep("Intercept|b_|sd_|sdgp_|lscale_", names(rh))], na.rm = TRUE))
  if (!is.finite(rmx)) rmx <- NA_real_
  out <- purrr::map_dfr(vars, function(v) {
    x <- exp(dr[[v]]); term <- sub("^b_", "", v)
    tibble::tibble(model = label, term = term, is_intercept = term == "Intercept",
      effect_scale = escale,
      hr = stats::median(x),
      lo = stats::quantile(x, 0.05, names = FALSE),
      hi = stats::quantile(x, 0.95, names = FALSE),
      p_dir = if (term == "Intercept") NA_real_ else mean(x > 1),  # posterior P(effect>1)
      rhat = rmx)
  })
  # tv-week: sigma_w, the SD of the weekly intercept, is the whole point of that model —
  # sigma_w -> 0 means "beta does not vary week to week". Without this it would be absent
  # from bayes_parameters.csv, since the loop above takes only population-level `b_` terms.
  # Reported as exp(sigma_w): the multiplicative spread in beta produced by a one-SD week,
  # so it sits on the same ratio axis as every other row with 1 = "no time-variation".
  # p_dir is NA — a standard deviation is non-negative, so P(> 1) is 1 by construction and
  # would read as a spurious "certain effect".
  sdv <- grep("^(sd_|sdgp_)", names(dr), value = TRUE)
  if (length(sdv)) {
    out <- dplyr::bind_rows(out, purrr::map_dfr(sdv, function(v) {
      x <- exp(dr[[v]])
      tibble::tibble(model = label, term = "sd_week", is_intercept = FALSE,
        effect_scale = paste0(escale, " (exp SD of weekly intercept)"),
        hr = stats::median(x),
        lo = stats::quantile(x, 0.05, names = FALSE),
        hi = stats::quantile(x, 0.95, names = FALSE),
        p_dir = NA_real_, rhat = rmx)
    }))
  }
  # GP LENGTHSCALE -> WEEK-TO-WEEK PERSISTENCE. For the ar1 (Ornstein-Uhlenbeck / exponential
  # kernel) variant the lengthscale l maps exactly onto the AR(1) coefficient,
  #     phi = exp(-1 / l)   per one-week step,
  # which is the interpretable quantity: phi -> 1 is a random walk (the level persists),
  # phi -> 0 is independent weeks (the level reverts to beta0 immediately). It is reported on
  # the SAME column as everything else, with 1 = "no reversion", so the row reads on the same
  # axis as sd_week. For the exp_quad (gp) variant exp(-1/l) is NOT the AR coefficient — that
  # kernel has no exact AR representation — so the lengthscale is reported in WEEKS instead
  # and the effect_scale column says which is which. Never quoted as a hazard ratio: it is not
  # one, and mislabelling it would put a time constant on a ratio axis.
  lsv <- grep("^lscale_", names(dr), value = TRUE)
  if (length(lsv)) {
    .tv <- attr(fit, "tv") %||% "none"
    out <- dplyr::bind_rows(out, purrr::map_dfr(lsv, function(v) {
      l <- dr[[v]]
      is_ou <- identical(.tv, "ar1")
      x <- if (is_ou) exp(-1 / l) else l
      tibble::tibble(model = label,
        term = if (is_ou) "ar1_phi_week" else "gp_lengthscale_weeks",
        is_intercept = FALSE,
        effect_scale = if (is_ou) "AR(1) coefficient per week (1 = random walk, 0 = independent weeks)"
                       else "GP lengthscale (weeks)",
        hr = stats::median(x),
        lo = stats::quantile(x, 0.05, names = FALSE),
        hi = stats::quantile(x, 0.95, names = FALSE),
        p_dir = NA_real_, rhat = rmx)
    }))
  }
  out
}

#' NUTS/HMC convergence diagnostics for a fitted brms model (review §4.4).
#'
#' The reviewer asked "what Stan settings and what convergence diagnostics?" — the
#' pipeline previously extracted only Rhat. This returns a one-row-per-model tibble
#' with divergent-transition count and %, minimum bulk/tail effective sample size,
#' maximum Rhat over the population parameters (b_*) AND the group-level SDs (sd_*), and
#' the post-warmup
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
      # sd_* included for the SAME reason bayes_param_table() includes it: for a tv-week
      # model sigma_w IS the parameter of interest, and it is typically the worst-mixing
      # one. Filtering to ^b_ let such a fit report a clean max-Rhat / min-ESS while its
      # weekly-intercept SD had not converged — a convergence table that cannot fail on
      # the parameter most likely to fail.
      sm <- sm[grepl("^(b_|sd_|sdgp_|lscale_)", sm$variable), , drop = FALSE]
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
#' INCLUDE_UNFILLED_MODELS = FALSE (the default) keeps only the source-cell fill twins of the
#' composites that have one (M8/M10/M13/M14/M16/M17 and the -dist twins of M8/M10/M13/M14/M17),
#' dropping the unfilled originals. INCLUDE_COHORT_SPLIT_MODELS adds the origin-split cohort family.
#'   * MOBILITY kernel: the CORE four distinct families — gravity (M4), composite-gravity
#'     (M8), multi-kernel ensemble (M9), radiation-composite (M10). (M1/M2/M3/M5/M6/M7/M4b
#'     are components or minor variants folded into these composites.) OPTIONAL: the inward
#'     meeting-location FOI kernel (M11; INCLUDE_M11_MODELS) and, for each family, an OSRM
#'     ROAD-DISTANCE (km) deterrence variant (M4-dist / ... / M11-dist; INCLUDE_OSRM_DIST_MODELS)
#'     — the same kernels keyed on distance instead of travel time. -dist entries are also
#'     auto-filtered out when the road-distance matrices were not built.
#'   * GENERATION TIME: a SINGLE medium anchor for the selection grid (review §2.1 — the GT
#'     is no longer a selectable axis; it is a PRIOR marginalised into the featured forecast).
#'   * COVARIATES: none (CORE), always fitted. The reduced "geo" subset (log_pop, CCVI,
#'     d_min) is ON by default (INCLUDE_GEO_COV_MODELS). The FULL-exogenous set (geo +
#'     healthsite_density) is OFF by default (INCLUDE_FULL_COV_MODELS) — this docstring used
#'     to call it "always kept", but the code gates it and the flag defaults FALSE, so no
#'     -full model has been fitted. The suspected-but-not-confirmed leading indicators —
#'     susp_local (own preceding-week suspected cases) and susp_import (mobility-weighted
#'     import of other zones' suspected cases), as a susp-only model and a wider FULL+susp
#'     model — are likewise OFF by default (INCLUDE_SUSPECTED_COV_MODELS), not "default ON".
#'     Verified against the shipped run: of 52 cross-validated Bayesian models, none is a
#'     -full or -susp model. positivity and the
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
  # Every ifnotfound below MUST equal the flag's default in 00_config.R, or a config-less
  # caller silently gets a different grid from the pipeline's.
  # Both inc_dist and inc_m11 were left at TRUE here when 00_config.R turned them OFF
  # (143b5a8, 2026-09-21), so a config-less caller composed the old 52-spec grid while the
  # pipeline composed 26. Re-synced 2026-09-21. The test fixture (.grid_flags in
  # tests/test_bayes_projection.R) pins both, which is why no test caught the divergence.
  inc_dist  <- isTRUE(get0("INCLUDE_OSRM_DIST_MODELS",     ifnotfound = FALSE))
  inc_geo   <- isTRUE(get0("INCLUDE_GEO_COV_MODELS",       ifnotfound = TRUE))
  inc_full  <- isTRUE(get0("INCLUDE_FULL_COV_MODELS",      ifnotfound = FALSE))
  inc_m11   <- isTRUE(get0("INCLUDE_M11_MODELS",           ifnotfound = FALSE))
  inc_logit <- isTRUE(get0("INCLUDE_LOGIT_SENS_MODELS",    ifnotfound = FALSE))
  inc_susp  <- isTRUE(get0("INCLUDE_SUSPECTED_COV_MODELS", ifnotfound = FALSE))
  inc_flowstat <- isTRUE(get0("INCLUDE_FLOWSTATIC_MODELS", ifnotfound = TRUE))
  # M9 (short-trip + M4/M5/M6a ensemble) and M15 (symmetrised OD kernel O + t(O)) are OFF by
  # default — their matrices are not built (03_mobility_matrices.R), so the trailing availability
  # filter would drop them anyway; gating here keeps the grid self-documenting.
  inc_m9    <- isTRUE(get0("INCLUDE_M9_MODELS",            ifnotfound = FALSE))
  inc_m15   <- isTRUE(get0("INCLUDE_M15_MODELS",           ifnotfound = FALSE))
  # Source-cell fill variants (M8/M13/M14/M16/M17-fill): same composites with the
  # cells the empirical source could not observe taken from the base kernel rather
  # than asserted as zero (03_mobility_matrices.R, compose_epicentre(fill=)). Their
  # matrices are only built when INCLUDE_SOURCEFILL_MODELS is on and the fill level
  # is not "none", and the availability filter at the end drops any that are absent.
  inc_fill  <- isTRUE(get0("INCLUDE_SOURCEFILL_MODELS",    ifnotfound = TRUE)) &&
               !identical(get0("MOBILITY_SOURCE_FILL", ifnotfound = "unmeasured"), "none")
  # INCLUDE_UNFILLED_MODELS = FALSE (the default) drops every model on a kernel whose source-cell
  # fill twin is in this grid and built (filter just before the time-varying block).
  inc_unfilled <- isTRUE(get0("INCLUDE_UNFILLED_MODELS", ifnotfound = FALSE))
  # Origin-split cohort composites. OFF by default since 2026-09-22 — this ifnotfound MUST
  # track INCLUDE_COHORT_SPLIT_MODELS in 00_config.R (see the note on inc_dist above).
  inc_split <- isTRUE(get0("INCLUDE_COHORT_SPLIT_MODELS", ifnotfound = FALSE))
  # Cohort-CALIBRATED gravity: the M4c base kernel and the M13c composite on it (plus the
  # -dist twins). M4c is gravity with its deterrence fitted to cohort presence rather than to
  # the relocation table; it is also the member that carries cohort information into the M17
  # consensus. Auto-filtered below when the matrices are absent.
  inc_cgrav <- isTRUE(get0("INCLUDE_COHORT_GRAVITY_MODELS", ifnotfound = TRUE))
  # COVARIATE SETS (revised 2026-09-21; `log_pop` dropped from both).
  #
  # There are only 54 invasion events in the pooled at-risk design (n_obs = 8,623), so at the
  # conventional ~10 events per parameter the honest ceiling is 2-3 covariates: the retired
  # 3-covariate `geo` gave 13.5 events/parameter and the 4-covariate `full` gave 10.8.
  #
  # `log_pop` was dropped on evidence, not taste. Fitting the cloglog hazard with the mobility
  # offset on five composite kernels (M14-fill, M8-fill, M17-fill, M13-fill, M16-fill), the
  # 2-covariate set beat the retired 3-covariate set on FOUR of the five and lost by 0.8 AIC on
  # the fifth, while spending one parameter fewer:
  #        kernel      retired geo (dAIC)   ccvi + d_min (dAIC)
  #        M14-fill          -0.3                 -1.7
  #        M8-fill          -37.1                -36.3
  #        M17-fill         -11.0                -11.4
  #        M13-fill         -21.9                -22.4
  #        M16-fill          +0.6                 -1.0
  # (dAIC against the intercept-only model at the same offset; negative is better.) Alone,
  # log_pop scores +0.6 on M14-fill — worse than omitting it. Population already enters through
  # the gravity/radiation kernels and the FOI offset, so as a covariate it partly re-estimates
  # what the offset already carries.
  #
  # `ccvi` is kept over `healthsite_density`: the two are correlated at r = -0.553 (so only one
  # belongs), and `healthsite_density + d_min` is worse than `ccvi + d_min` on all five kernels.
  #
  # `d_min` is KEPT despite being 62% correlated with the offset (r = -0.624). It is not a
  # nuisance term: it is a SUBSTITUTE for kernel quality. Added to `ccvi` it gains 38 AIC on
  # M8-fill, 24 on M13-fill and 12 on M17-fill, and COSTS ~2 on M14-fill and M16-fill — the
  # two kernels that already encode spatial proximity. That asymmetry is itself a finding — and
  # a caution,
  # since a covariate that can compensate for a weak kernel COMPRESSES differences between
  # kernels. The non-`geo` half of the grid is the arm that reads kernel quality directly.
  # Read from 00_config.R so the grid and the published parameter trace cannot diverge.
  geo  <- get0("BAYES_GEO_COVARIATES",  ifnotfound = c("ccvi", "d_min"))
  full <- get0("BAYES_FULL_COVARIATES", ifnotfound = c("ccvi", "d_min", "healthsite_density"))
  susp <- c("susp_import", "susp_local")
  # TIME-VARYING beta families (INCLUDE_TV_BETA_MODELS). Five strict generalisations of the
  # fixed-beta model, each collapsing back to it as its extra parameter -> 0:
  #   tv-trend : beta_t = exp(beta0 + gamma_w * z(week))         log-linear, EXTRAPOLATES
  #   tv-week  : u_t ~ N(0, sigma) iid                           reverts to beta0, wider CrI
  #   tv-rw1   : u_t = u_{t-1} + eps_t                           PERSISTS the current level
  #   tv-ar1   : Ornstein-Uhlenbeck, phi = exp(-1/lengthscale)   decays toward beta0
  #   tv-gp    : squared-exponential GP                          smooth local extrapolation
  # See the TIME-VARYING block above fit_bayes_renewal_tv() for the full argument.
  #
  # THEY RIDE ONE BASE MODEL, NOT THE WHOLE GRID. These are variants OF THE MODEL THE
  # CROSS-VALIDATION SELECTED, which is what the question "has beta changed since the
  # epidemic left Ituri?" is about. Sweeping them across every kernel would multiply the
  # grid sixfold, refit in every fold, and answer a question nobody asked.
  inc_tv <- isTRUE(get0("INCLUDE_TV_BETA_MODELS", ifnotfound = FALSE))
  tv_types <- get0("TV_BETA_TYPES", ifnotfound = c("trend", "week", "rw1", "ar1", "gp"))
  g <- list()
  # `tv` is carried on EVERY spec (default "none") so downstream dispatch never has to
  # test for the field's existence — a missing field would silently take the fixed-beta
  # branch and the tv models would be indistinguishable from their base kernels.
  # `tv` and `sensitivity` are carried on EVERY spec (defaults "none"/FALSE) so downstream
  # dispatch never has to test for a field's existence — a missing field would silently take
  # the fixed-beta / selectable branch.
  add <- function(mob, gt, cov, link, label, tv = "none", sensitivity = FALSE)
    g[[length(g) + 1L]] <<- list(mob = mob, gt = gt, cov = cov, link = link,
                                 label = label, tv = tv, sensitivity = isTRUE(sensitivity))

  # --- CORE (always): the distinct FOI families at medium GT, no covariates ---
  add("M4",  "medium", character(0), "cloglog", "Bayes-M4-med")
  # M4c: the same gravity form with the deterrence calibrated on cohort presence instead of
  # relocations. Entered standalone (not only inside the M17 consensus) ON PURPOSE, so the
  # cross-validation prices the steeper deterrence directly against M4's.
  if (inc_cgrav) add("M4c", "medium", character(0), "cloglog", "Bayes-M4c-med")
  add("M8",  "medium", character(0), "cloglog", "Bayes-M8-med")
  if (inc_m9) add("M9", "medium", character(0), "cloglog", "Bayes-M9-med")
  add("M10", "medium", character(0), "cloglog", "Bayes-M10-med")
  # --- Flowminder cohort composites (M13/M14 + geographic-distance -dist), default ON ---
  # Cohort presence source rows + gravity/radiation, on travel-time (M13/M14) or road-km
  # geographic distance (M13/M14-dist). Auto-filtered below when a matrix is absent
  # (INCLUDE_COHORT_MODELS off, or no road-distance matrix). The cohort -dist MATRICES are
  # built under the cohort family (no INCLUDE_OSRM_DIST_MODELS needed), but the -dist MODELS
  # below are gated on inc_dist like every other road-distance variant.
  if (isTRUE(get0("INCLUDE_COHORT_MODELS", ifnotfound = TRUE))) {
    add("M13",      "medium", character(0), "cloglog", "Bayes-M13-med")
    # M13c: the same composite over the COHORT-CALIBRATED gravity base (M4c), so both the
    # source rows and the deterrence elsewhere come from outbreak-period mobility.
    if (inc_cgrav) add("M13c", "medium", character(0), "cloglog", "Bayes-M13c-med")
    add("M14",      "medium", character(0), "cloglog", "Bayes-M14-med")
    # Road-distance cohort twins only when the -dist axis is on (inc_dist).
    if (inc_dist) {
      add("M13-dist", "medium", character(0), "cloglog", "Bayes-M13-dist")
      if (inc_cgrav) add("M13c-dist", "medium", character(0), "cloglog", "Bayes-M13c-dist")
      add("M14-dist", "medium", character(0), "cloglog", "Bayes-M14-dist")
    }
    # Covariate variants: FULL-exogenous when inc_full, the reduced geo subset when inc_geo — on the
    # travel-time (M13/M14) and, when inc_dist, the road-km (M13/M14-dist) cohort composites.
    # Auto-filtered below when a matrix is absent.
    for (cm in c("M13", if (inc_cgrav) "M13c", "M14",
                 if (inc_dist) c("M13-dist", if (inc_cgrav) "M13c-dist", "M14-dist"))) {
      if (inc_full) add(cm, "medium", full, "cloglog", paste0("Bayes-", cm, "-full"))
      if (inc_geo)  add(cm, "medium", geo,  "cloglog", paste0("Bayes-", cm, "-geo"))
    }
  }
  # --- Flowminder OD family (M15/M16/M17 + M17-dist), default ON ---
  # M15 = symmetrised Flowminder OD, S = O + t(O) (NOT inflow-informed: no independent
  # inflow table exists); M16 = cohort + DIRECTED OD
  # (cohort source rows where available, the M3 directed Flowminder OD elsewhere — NOT M15, which
  # is not built by default; see 03_mobility_matrices.R); M17/M17-dist = all-kernel consensus.
  # Each kernel is entered at the single medium-GT anchor (generation time is marginalised into
  # the featured forecast, not a grid axis) x covariates {none, geo (inc_geo), full (inc_full)}.
  # Auto-filtered below when a matrix is absent (M16 needs the cohort kernel; M17-dist needs a
  # road-distance matrix), so a cohort-off / distance-off run simply omits them.
  if (inc_flowstat) {
    for (mm in c(if (inc_m15) "M15", "M16", "M17", if (inc_dist) "M17-dist")) {
      # Single medium GT ANCHOR. GT is neither a grid axis nor marginalised: it is treated as
      # a known distribution (GT_MARGINALISE_FEATURED = FALSE, 00_config.R).
      # none + (full when inc_full) + (geo when inc_geo) covariate variants.
      add(mm, "medium", character(0), "cloglog", paste0("Bayes-", mm, "-med"))
      if (inc_full) add(mm, "medium", full, "cloglog", paste0("Bayes-", mm, "-full"))
      if (inc_geo)  add(mm, "medium", geo,  "cloglog", paste0("Bayes-", mm, "-geo"))
    }
  }
  # --- SOURCE-CELL FILL twins of the composites already in the grid -----------------
  # One entry per fill kernel, mirroring the covariate treatment of its parent, so the
  # comparison is like-for-like and the cross-validation decides whether refilling the
  # unobservable cells helps. M16-fill/M17-fill additionally ride on inc_flowstat.
  if (inc_fill) {
    for (mm in c("M8-fill", "M10-fill", "M13-fill", if (inc_cgrav) "M13c-fill", "M14-fill",
                 if (inc_flowstat) c("M16-fill", "M17-fill"))) {
      add(mm, "medium", character(0), "cloglog", paste0("Bayes-", mm, "-med"))
      if (inc_full) add(mm, "medium", full, "cloglog", paste0("Bayes-", mm, "-full"))
      if (inc_geo)  add(mm, "medium", geo,  "cloglog", paste0("Bayes-", mm, "-geo"))
    }
    # Road-distance fill twins. Labels follow their unfilled parents: M8/M10/M13/M14-dist carry
    # no "-med" token (Bayes-M13-dist), M17-dist does (Bayes-M17-dist-med).
    if (inc_dist) {
      for (mm in c("M8-dist-fill", "M10-dist-fill", "M13-dist-fill",
                   if (inc_cgrav) "M13c-dist-fill", "M14-dist-fill",
                   if (inc_flowstat) "M17-dist-fill")) {
        add(mm, "medium", character(0), "cloglog",
            paste0("Bayes-", mm, if (identical(mm, "M17-dist-fill")) "-med" else ""))
        if (inc_full) add(mm, "medium", full, "cloglog", paste0("Bayes-", mm, "-full"))
        if (inc_geo)  add(mm, "medium", geo,  "cloglog", paste0("Bayes-", mm, "-geo"))
      }
    }
  }
  # --- ORIGIN-SPLIT cohort composites (INCLUDE_COHORT_SPLIT_MODELS, default ON) -------------
  # The pooled cohort profile disaggregated per origin (split_cohort_rows(), 03), then filled.
  # Same label convention as above. Auto-filtered below when a matrix is absent.
  # The SHORT-TRIP splits (M8/M10-split) ride on inc_split alone — they need the short-trip
  # annex, not the cohort tables, so they stay in the grid when INCLUDE_COHORT_MODELS is off.
  # Keep only the split families SPLIT_FAMILIES names (00_config.R carries the fold evidence).
  # Strip "-split" first, then "-dist", so "M13-dist-split" resolves to the family "M13".
  .split_fams <- get0("SPLIT_FAMILIES", ifnotfound = c("M8", "M14"))
  .split_keep <- function(v) v[sub("-dist$", "", sub("-split$", "", v)) %in% .split_fams]
  if (inc_split) {
    for (mm in .split_keep(c("M8-split", "M10-split",
                 if (inc_dist) c("M8-dist-split", "M10-dist-split")))) {
      add(mm, "medium", character(0), "cloglog", paste0("Bayes-", mm, "-med"))
      if (inc_full) add(mm, "medium", full, "cloglog", paste0("Bayes-", mm, "-full"))
      if (inc_geo)  add(mm, "medium", geo,  "cloglog", paste0("Bayes-", mm, "-geo"))
    }
  }
  if (inc_split && isTRUE(get0("INCLUDE_COHORT_MODELS", ifnotfound = TRUE))) {
    for (mm in .split_keep(c("M13-split", if (inc_cgrav) "M13c-split", "M14-split",
                 if (inc_flowstat) c("M16-split", "M17-split"),
                 if (inc_dist) c("M13-dist-split", if (inc_cgrav) "M13c-dist-split",
                                 "M14-dist-split",
                                 if (inc_flowstat) "M17-dist-split")))) {
      no_med <- mm %in% c("M13-dist-split", "M13c-dist-split", "M14-dist-split")
      add(mm, "medium", character(0), "cloglog",
          paste0("Bayes-", mm, if (no_med) "" else "-med"))
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
    if (inc_cgrav) add("M4c", "medium", geo, "cloglog", "Bayes-M4c-geo")
    if (inc_m9) add("M9", "medium", geo, "cloglog", "Bayes-M9-geo")
    add("M10", "medium", geo, "cloglog", "Bayes-M10-geo")
  }
  # --- OPTIONAL: inward / meeting-location FOI kernel (M11) ---
  if (inc_m11) {
    add("M11", "medium", character(0), "cloglog", "Bayes-M11-inward")
    # The FULL covariate set is gated on inc_full like every other family: it was added
    # unconditionally here, so switching INCLUDE_FULL_COV_MODELS off still fitted M11-full.
    if (inc_full) add("M11", "medium", full, "cloglog", "Bayes-M11-full")
    if (inc_geo) add("M11", "medium", geo, "cloglog", "Bayes-M11-geo")
  }
  # --- OPTIONAL: OSRM ROAD-DISTANCE (km) deterrence kernels — same families keyed on
  #     distance instead of travel time. Auto-filtered below if the -dist matrices are absent. ---
  if (inc_dist) {
    for (mm in c("M4-dist", if (inc_cgrav) "M4c-dist", "M8-dist",
                 if (inc_m9) "M9-dist", "M10-dist")) {
      add(mm, "medium", character(0), "cloglog", paste0("Bayes-", mm))
      if (inc_full) add(mm, "medium", full, "cloglog", paste0("Bayes-", mm, "-full"))
      if (inc_geo)  add(mm, "medium", geo,  "cloglog", paste0("Bayes-", mm, "-geo"))
    }
    if (inc_m11) {
      add("M11-dist", "medium", character(0), "cloglog", "Bayes-M11-dist")
      if (inc_full) add("M11-dist", "medium", full, "cloglog", "Bayes-M11-dist-full")
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

  # --- OPTIONAL: fill-only composites (INCLUDE_UNFILLED_MODELS = FALSE) ---------------------
  # A kernel is "unfilled" when its source-cell fill twin ("<kernel>-fill") is in this grid AND
  # built. Every model on such a kernel is dropped — all covariate, suspected-case and link
  # variants — and the time-varying variants below are derived from what remains, so they
  # follow. Kernels with no fill twin (M4, M4-dist, M9, M9-dist, M11, M15, and the origin-split
  # family) are untouched. With no twin available nothing is dropped, so turning the fill
  # family off can never silently remove the composites altogether.
  if (!inc_unfilled) {
    built <- function(m) is.null(mobility_matrices) || !is.null(mobility_matrices[[m]])
    mobs  <- vapply(g, function(x) x$mob, character(1))
    twins <- unique(mobs[grepl("-fill$", mobs)])
    twins <- twins[vapply(twins, built, logical(1))]
    if (length(twins)) {
      g <- Filter(function(x) !(paste0(x$mob, "-fill") %in% twins), g)
    } else {
      warning("[bayes grid] INCLUDE_UNFILLED_MODELS = FALSE, but no source-cell fill kernel ",
              "is in the grid and built; no model was dropped.", call. = FALSE)
    }
  }

  # --- TIME-VARYING beta variants OF THE SELECTED MODEL -----------------------------------
  # The base spec is the model the LAST cross-validation selected (model_selection.json, the
  # pipeline's own record — the same file the cascade and the figure suites read), so the tv
  # arm is a set of variants of THE featured model rather than of an arbitrary kernel. The
  # variants inherit its mobility kernel, generation time, link AND covariates, so the only
  # thing that differs between the base model and each variant is the beta_t process — which
  # is what makes the cross-validated comparison between them interpretable.
  #
  # RESOLUTION ORDER, all overridable:
  #   1. TV_BETA_BASE_MODEL (00_config.R / environment) — an explicit label.
  #   2. the featured Bayesian model in model_selection.json.
  #   3. the richest spec on MOBILITY_PRIMARY, so a first-ever run with no selection file
  #      still produces a coherent tv arm rather than silently producing none.
  # A label that resolves to no spec in this grid is a hard error, not a warning: silently
  # attaching the tv variants to a different model than the one being generalised would make
  # every comparison in the resulting figure false.
  # sensitivity = TRUE on every arm: these are REFITS of a model already in the grid, not
  # competing hypotheses, so best_invasion_model() must never be able to feature one. The
  # deployed gate is INVASION_SELECTION_EXCLUDE (00_config.R), which matches on the label;
  # the flag here is the grid's own record of the same fact and is what a future change
  # should read instead of a regex.
  if (inc_tv && length(tv_types)) {
    base <- .bayes_featured_base_spec(g)
    for (ty in tv_types) {
      lab <- paste0(base$label, "-tv", ty)
      if (identical(ty, "trend")) {
        # tv-trend rides the ordinary covariate machinery: `week_idx` IS the trend term, so it
        # is ADDED to the base model's covariates rather than replacing them.
        add(base$mob, base$gt, unique(c(base$cov, "week_idx")), base$link, lab, tv = "trend",
            sensitivity = TRUE)
      } else {
        # week / rw1 / ar1 / gp are structural terms handled by fit_bayes_renewal_tv();
        # the base model's covariates are carried through unchanged.
        add(base$mob, base$gt, base$cov, base$link, lab, tv = ty, sensitivity = TRUE)
      }
    }
    message(sprintf("[bayes grid] time-varying beta variants of %s: %s",
                    base$label, paste(paste0("tv", tv_types), collapse = ", ")))
  }

  # --- GENERATION-TIME SENSITIVITY OF THE SELECTED MODEL ----------------------------------
  # The generation time is an ASSUMPTION, not a fitted quantity: it enters the renewal
  # convolution that builds the import force Lambda, so a different GT gives a different
  # offset for every row, and beta0 absorbs the difference. The suite is composed at ONE
  # anchor (medium) precisely so that GT is not a selection axis — but that makes "how much
  # does the GT assumption matter?" a question the main grid cannot answer.
  #
  # This adds the SELECTED model refit at the short and long GT profiles, nothing else. They
  # are ordinary suite members, so they are cross-validated on the same folds (answering
  # part i: does the GT assumption change out-of-sample skill?) and they produce a current
  # forecast (answering part ii: does it change today's invasion probabilities?).
  #
  # These are SENSITIVITY arms, not candidates. They share the featured model's kernel and
  # covariates and differ only in an assumption, so letting one of them win selection would
  # be selecting on the GT — exactly what the design excludes. They are marked here and
  # best_invasion_model() excludes any model carrying sensitivity_arm = TRUE.
  if (isTRUE(get0("INCLUDE_GT_SENSITIVITY_MODELS", ifnotfound = TRUE))) {
    base <- .bayes_featured_base_spec(g)
    alt_gt <- setdiff(get0("GT_SENSITIVITY_PROFILES", ifnotfound = c("short", "long")), base$gt)
    for (gk in alt_gt) {
      # Strip the base label's own GT token before appending the alternative one, so
      # "Bayes-M14-fill-med" becomes "Bayes-M14-fill-short", not "...-med-short".
      stem <- sub("-(med|short|long)$", "", base$label)
      add(base$mob, gk, base$cov, base$link, paste0(stem, "-gt", gk), sensitivity = TRUE)
    }
    if (length(alt_gt))
      message(sprintf("[bayes grid] generation-time sensitivity arms of %s: %s",
                      base$label, paste(alt_gt, collapse = ", ")))
  }

  if (!is.null(mobility_matrices)) g <- Filter(function(x) !is.null(mobility_matrices[[x$mob]]), g)
  g
}

#' The spec in `g` that the time-varying-beta variants generalise.
#'
#' Kept out of bayes_default_grid() so the resolution rule is testable on its own and so the
#' three sources it consults are visible in one place. Errors rather than guessing: see the
#' call site for why a wrong base silently invalidates the whole tv comparison.
.bayes_featured_base_spec <- function(g) {
  cand <- Filter(function(x) identical(x$tv %||% "none", "none"), g)
  if (!length(cand))
    stop("[bayes grid] INCLUDE_TV_BETA_MODELS is on but the grid holds no fixed-beta spec ",
         "to generalise.", call. = FALSE)
  pick <- function(lab) {
    if (is.null(lab) || !length(lab) || is.na(lab[1]) || !nzchar(lab[1])) return(NULL)
    hit <- Filter(function(x) identical(x$label, as.character(lab[1])), cand)
    if (length(hit)) hit[[1]] else NULL
  }
  # 1. explicit override
  lab_cfg <- get0("TV_BETA_BASE_MODEL", ifnotfound = NULL)
  b <- pick(lab_cfg)
  if (!is.null(b)) return(b)
  if (!is.null(lab_cfg) && length(lab_cfg) && !is.na(lab_cfg[1]) && nzchar(lab_cfg[1]))
    stop(sprintf(paste0("[bayes grid] TV_BETA_BASE_MODEL = '%s' matches no fixed-beta spec in ",
                        "the grid. Available: %s"), lab_cfg[1],
                 paste(vapply(cand, function(x) x$label, character(1)), collapse = ", ")),
         call. = FALSE)
  # 2. the pipeline's own selection record. STRIP any variant suffix the previous run's
  #    winner may carry: if a tv variant or a GT arm was featured last time, the label is
  #    e.g. "Bayes-M14-fill-geo-tvrw1", and building variants on THAT would produce
  #    "...-tvrw1-tvrw1" — a label no grid generates and a model nothing can fit. The base
  #    of a variant is its fixed-beta, anchor-GT parent, which is what belongs here.
  lab_sel <- tryCatch({
    f <- file.path(get0("OUT_DIR", ifnotfound = file.path(here::here(), "spatiotemporal", "outputs")),
                   "key_outputs", "model_selection.json")
    if (!file.exists(f)) NULL else {
      sel <- jsonlite::fromJSON(f, simplifyVector = TRUE)
      m <- sel$featured$bayesian$method
      if (is.null(m) || !length(m) || is.na(m[1])) m <- sel$featured$headline$method
      if (is.null(m) || !length(m) || is.na(m[1])) NULL
      else sub("-(tvtrend|tvweek|tvrw1|tvar1|tvgp|gtshort|gtlong)$", "", as.character(m[1]))
    }
  }, error = function(e) NULL)
  b <- pick(lab_sel)
  if (!is.null(b)) return(b)
  if (!is.null(lab_sel) && length(lab_sel) && !is.na(lab_sel[1]))
    warning(sprintf(paste0("[bayes grid] the selected model '%s' is not in the current grid, so ",
                           "the time-varying variants cannot be built on it; falling back to the ",
                           "primary kernel. Re-point TV_BETA_BASE_MODEL, or re-run selection."),
                    lab_sel[1]), call. = FALSE)
  # 3. fallback: the richest (most covariates) spec on the primary kernel, then the richest
  #    spec of all. Deterministic ties: first by covariate count, then by label order.
  prim <- as.character(get0("MOBILITY_PRIMARY", ifnotfound = NA_character_))
  on_prim <- Filter(function(x) identical(x$mob, prim) && identical(x$link, "cloglog"), cand)
  pool <- if (length(on_prim)) on_prim else cand
  ord <- order(-vapply(pool, function(x) length(x$cov), integer(1)),
               vapply(pool, function(x) x$label, character(1)), method = "radix")
  b <- pool[[ord[1]]]
  warning(sprintf(paste0("[bayes grid] no usable model selection was found; the time-varying ",
                         "variants are built on '%s'. They generalise THAT model, not ",
                         "necessarily the one this run selects."), b$label), call. = FALSE)
  b
}

#' Label of a model's source-cell fill twin, as bayes_default_grid() names it.
#'
#' "Bayes-M8-med" -> "Bayes-M8-fill-med", "Bayes-M13-dist-geo" -> "Bayes-M13-dist-fill-geo".
#' The twin id is built from the model's OWN kernel (mobility_kernel_from_method) and returned
#' only when that id exists in MOBILITY_FILL_IDS, so:
#'   * -dist models map to their -dist-fill twin where one exists (M8/M10/M13/M14/M17), NA otherwise;
#'   * -split and already-filled models return NA (they are filled already — inserting "-fill"
#'     produced labels such as "Bayes-M13-fill-split-med" that no grid ever generates);
#'   * a label whose kernel is unknown returns NA.
#' The suffix (-med/-geo/-susp/-logit/...) is preserved exactly. Vectorised.
bayes_fill_twin_label <- function(label) {
  fill_ids <- get0("MOBILITY_FILL_IDS", ifnotfound = character(0))
  vapply(label, function(lb) {
    if (is.na(lb)) return(NA_character_)
    kern <- mobility_kernel_from_method(lb)
    if (is.na(kern) || grepl("-(fill|split)$", kern)) return(NA_character_)
    twin <- paste0(kern, "-fill")
    if (!twin %in% fill_ids) return(NA_character_)
    sfx <- sub(sprintf("^Bayes-%s", kern), "", lb)   # "" or "-geo", "-full-susp", ...
    if (!nzchar(sfx) && !identical(lb, paste0("Bayes-", kern))) return(NA_character_)
    paste0("Bayes-", twin, sfx)
  }, character(1), USE.NAMES = FALSE)
}

#' Identifier of the h >= 2 source projection used by the Bayesian suite. Stamped on the
#' cached cross-validation so a cache produced under a different projection is never reused.
BAYES_H2_PROJECTION <- "epinow2-R-draws x per-draw-posterior-hazard v1"

#' Named-by-horizon recalibration factors for one suite model, or NULL if none apply.
#'
#' Returns NULL (rather than a vector of 1s) when the method has no usable row, so the
#' prediction path takes the untouched code path exactly as before.
.suite_delta_for <- function(delta_tbl, label, horizons) {
  if (is.null(delta_tbl) || !nrow(delta_tbl)) return(NULL)
  if (!all(c("method", "horizon", "delta") %in% names(delta_tbl))) {
    warning("[recal] delta_tbl lacks method/horizon/delta; suite forecasts left unrecalibrated.",
            call. = FALSE)
    return(NULL)
  }
  r <- delta_tbl[delta_tbl$method == label, , drop = FALSE]
  if (!nrow(r)) return(NULL)
  d <- stats::setNames(as.numeric(r$delta), as.character(r$horizon))
  d <- d[is.finite(d) & d > 0]
  d <- d[names(d) %in% as.character(horizons)]
  if (!length(d)) return(NULL)
  if (length(d) < length(horizons))
    warning(sprintf(paste0("[recal] %s: a delta is available for horizon(s) %s but not %s; ",
                           "the unmatched horizon(s) are left UNRECALIBRATED, so probabilities ",
                           "are not on a common footing across horizons."),
                    label, paste(names(d), collapse = ","),
                    paste(setdiff(as.character(horizons), names(d)), collapse = ",")),
            call. = FALSE)
  d
}

#' Fit the Bayesian SUITE across the comprehensive grid (mobility x GT x covariates x
#' observation link), predict the current forecast for each, collect posterior
#' parameters, and compute loo predictive-stacking weights.
#'
#' @param delta_tbl optional data.frame(method, horizon, delta) of POOLED post-hoc
#'   recalibration factors (16b_invasion_recalibration.R). When supplied, each suite
#'   model's current forecast is recalibrated with its OWN delta, applied to the
#'   cumulative hazard inside the posterior draws so every summary — mean, median,
#'   90% CrI, sd, infection scale — and the rank credible interval stay exact. A model
#'   with no matching row is left unrecalibrated. NULL (default) = no recalibration.
fit_bayes_suite <- function(zone_week_nc, mobility_matrices, gt_pmfs, covariates,
                            osrm_mat, zones_all, affected_zones, horizons,
                            grid = NULL, iter = 2000L, delta_tbl = NULL,
                            linelist = NULL, analysis_date = ANALYSIS_DATE,
                            # MUST match what the LFO uses. The recalibration delta is fitted
                            # on folds and applied to THIS suite's forecast, so if the folds
                            # use rolling as-of predictors and the deployed fit does not, beta0
                            # here is calibrated against a lambda about 2x smaller than the one
                            # it is applied to and delta -- now ~1 because the folds are
                            # calibrated -- no longer corrects it. Leaving this NULL while the
                            # LFO rolls is strictly worse than neither rolling.
                            trunc_delay = NULL) {
  if (is.null(grid)) grid <- bayes_default_grid(mobility_matrices)
  grid <- Filter(function(g) !is.null(mobility_matrices[[g$mob]]), grid)

  # ---- Pass 1: build every design (cheap; no MCMC yet) ----------------------
  designs <- list()
  for (g in grid) {
    .roll_on <- isTRUE(get0("ROLLING_PREDICTORS", ifnotfound = TRUE)) &&
                !is.null(linelist) && !is.null(trunc_delay)
    des <- tryCatch(build_invasion_design(zone_week_nc, mobility_matrices, gt_pmfs,
             covariates, osrm_mat, zones_all, mob = g$mob, gt = g$gt,
             candidates = unique(c(g$cov, "log_pop", "ccvi", "d_min")),
             rolling_ll    = if (.roll_on) linelist else NULL,
             rolling_delay = if (.roll_on) trunc_delay else NULL),
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

  # Week index -> calendar date, so every beta_t trajectory is dated rather than carrying a
  # bare column index. build_invasion_design() numbers `.week` by COLUMN POSITION in the count
  # matrix, so index w is the w-th week of the grid; forecast weeks continue the same spacing.
  .week_dates <- local({
    wk <- sort(unique(as.Date(zone_week_nc$week_start)))
    if (!length(wk)) return(NULL)
    idx <- seq_len(length(wk) + max(horizons))
    stats::setNames(as.character(wk[1] + 7L * (idx - 1L)), as.character(idx))
  })

  # ---- EpiNow2 R draws for the h >= 2 source projection, one fit per generation time --
  # Fitted HERE, sequentially, before the (possibly parallel) model fits: every model on a
  # given GT shares its R, and no worker waits on another's EpiNow2 fit. bayes_rt_week_draws()
  # always returns an R (EpiNow2, a retried fit, or the model's prior).
  rt_by_gt <- list()
  if (max(horizons) >= 2L)
    for (gk in unique(vapply(designs, function(x) x$g$gt, character(1))))
      rt_by_gt[[gk]] <- bayes_rt_week_draws(linelist, gt_pmfs[[gk]],
                                            week_start = max(zone_week_nc$week_start),
                                            issue_date = analysis_date)

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
    # tv = "week" needs the group-level fitter and the matching predictor (which injects the
    # forecast week's .week_f level). tv = "trend" is an ordinary covariate model — it must NOT
    # take the tv branch, or the trend term would be fitted twice over.
    .tv  <- g$tv %||% "none"
    fit <- if (.tv_is_structural(.tv))
             fit_bayes_renewal_tv(des_fit, cov_spec = g$cov, iter = iter,
                                  link = g$link %||% "cloglog", tv = .tv)
           else
             fit_bayes_renewal(des_fit, cov_spec = g$cov, iter = iter,
                               link = g$link %||% "cloglog")
    if (is.null(fit)) return(NULL)
    pr <- NULL; pr_raw <- NULL
    off <- tryCatch(bayes_forecast_offsets(zone_week_nc, mobility_matrices, gt_pmfs,
             covariates, osrm_mat, zones_all, g$mob, g$gt, g$cov, horizons,
             rt_draws = rt_by_gt[[g$gt]]), error = function(e) NULL)
    if (!is.null(off)) {
      .dl <- .suite_delta_for(delta_tbl, g$label, horizons)
      # BOTH probability scales from one posterior-linear-predictor pass: `pred` is the
      # deployed (recalibrated) forecast that every primary product uses, `pred_raw` its
      # uncorrected twin for the parallel raw artifact set. With no factor for this model
      # the two are the same object, which is the honest answer.
      .ps <- if (.tv_is_structural(.tv))
               predict_bayes_invasion_tv_scales(fit, off, des, horizons, affected_zones, delta = .dl)
             else
               predict_bayes_invasion_scales(fit, off, des, horizons, affected_zones, delta = .dl)
      if (!is.null(.ps)) {
        pr <- .ps$calibrated; pr_raw <- .ps$raw
        if (!is.null(pr))     pr$method     <- g$label
        if (!is.null(pr_raw)) pr_raw$method <- g$label
      }
    }
    # Compute the model's aligned pointwise LOO here, while the suite workers are
    # already fanned out. Previously this was a second sequential pass over every
    # fit and took several minutes after the parallel sampling had finished.
    lo <- NULL
    if (length(common_keys) >= 10L) {
      dd <- des$d
      dd <- dd[key_of(dd) %in% common_keys, , drop = FALSE]
      dd <- dd[order(key_of(dd)), , drop = FALSE]
      # The tv-week fit's formula references .week_f, which fit_bayes_renewal_tv() added to
      # its own local copy of the design — NOT to des$d. Without recreating it here brms
      # errors on the missing grouping variable, this model's LOO comes back NULL, and
      # stacking then fails for the ENTIRE suite (the weights block requires every model's
      # LOO). These rows are all training weeks, so no new levels arise.
      # Recreate the columns the tv formulas reference. fit_bayes_renewal_tv() adds .week_f to
      # its OWN local copy of the design, not to des$d, so without this brms errors on a
      # missing grouping variable, this model's LOO comes back NULL, and stacking then fails
      # for the ENTIRE suite. The factor LEVELS must be the fitted ones (attr "train_weeks"),
      # not the levels of this subset, or brms maps the same week to a different group.
      # `.week` itself is already present (build_invasion_design writes it) and is what the
      # ar1/gp terms read. These rows are all training weeks, so no new levels arise.
      if (.tv_is_structural(.tv)) {
        .twk <- attr(fit, "train_weeks")
        if (is.null(.twk)) .twk <- sort(unique(as.integer(dd$.week)))
        dd$.week_f <- factor(as.integer(dd$.week), levels = as.character(.twk))
      }
      lo <- tryCatch({
        llm <- brms::log_lik(fit, newdata = dd, allow_new_levels = TRUE,
                             sample_new_levels = "gaussian")
        suppressWarnings(loo::loo(llm))
      }, error = function(e) NULL)
    }
    # Posterior trajectory of beta_t, computed HERE because this is the only place the fit and
    # its design are both in scope, and computed for EVERY model (a fixed-beta model's flat
    # line is the reference the time-varying variants are read against).
    btj <- tryCatch(bayes_beta_trajectory(fit, des, g$label, horizons = horizons,
                                          week_dates = .week_dates),
                    error = function(e) NULL)
    list(fit = fit, pred = pr, pred_raw = pr_raw,
         params = bayes_param_table(fit, g$label), loo = lo, beta_traj = btj)
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
  preds_raw <- Filter(Negate(is.null), lapply(fit_res, `[[`, "pred_raw"))
  params <- lapply(fit_res, `[[`, "params")
  beta_traj <- Filter(Negate(is.null), lapply(fit_res, `[[`, "beta_traj"))

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
       preds_raw = dplyr::bind_rows(preds_raw),
       params = dplyr::bind_rows(params),
       beta_traj = if (length(beta_traj)) dplyr::bind_rows(beta_traj) else NULL,
       weights = wts, stacked = stacked, grid = grid)
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
      .groups = "drop") %>%
    # NaN -> NA. For an already-invaded zone every member carries NA (.finalise_invasion_summary
    # masks them), and weighted.mean(all-NA, na.rm = TRUE) returns NaN, not NA — so the stacked
    # product used a different missing-value sentinel from the single-model product it is meant
    # to differ from "only in the combination rule". The shipped bayes_stacked_current.rds had
    # NaN in 122 rows where the featured model had NA_real_.
    dplyr::mutate(dplyr::across(c(p_invasion, mu_forecast, .var_pool),
                                ~ ifelse(is.nan(.x), NA_real_, .x))) %>%
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
                                 link = "cloglog", tv = "none", linelist = NULL,
                                 # The regime's fitted truncation. Supplying it (with a line
                                 # list) switches the design to ROLLING as-of predictors;
                                 # NULL keeps the historic shared-matrix behaviour.
                                 trunc_delay = NULL) {
  # FORCE the per-spec arguments so each closure captures ITS OWN mob/gt/cov/link/tv.
  # Without this, R's lazy promises defer evaluation until the LFO runs the closure —
  # by which time the build loop has finished and `sp` holds the LAST spec, so every
  # model silently computes identical (last-spec) predictions. `tv` MUST be forced too,
  # or every model would inherit the last spec's time-varying setting.
  force(mob); force(gt); force(cov_spec); force(link); force(iter); force(tv); force(linelist)
  force(trunc_delay)   # same lazy-promise trap as the others: without this every closure
                       # would capture the LAST spec's delay once the build loop finished.
  function(zw, ti, hz, cut) {
    # R for the h >= 2 source projection, as known at this fold's forecast moment (cut + 6,
    # the issue date of the fold's training counts). bayes_rt_week_draws() always returns an R
    # (EpiNow2, a retried fit, or the model's prior), so every fold carries its h >= 2 forecast.
    rt <- if (max(hz) >= 2L)
            bayes_rt_week_draws(linelist, gt_pmfs[[gt]], week_start = max(zw$week_start),
                                issue_date = as.Date(cut) + 6L) else NULL
    .roll_on <- isTRUE(get0("ROLLING_PREDICTORS", ifnotfound = TRUE)) &&
                !is.null(linelist) && !is.null(trunc_delay)
    des <- tryCatch(build_invasion_design(zw, mobility_matrices, gt_pmfs, covariates,
             osrm_mat, zones_all, mob = mob, gt = gt,
             candidates = unique(c(cov_spec, "log_pop", "ccvi", "d_min", "healthsite_density")),
             rolling_ll    = if (.roll_on) linelist else NULL,
             rolling_delay = if (.roll_on) trunc_delay else NULL),
             error = function(e) NULL)
    if (is.null(des) || des$n_events < 1L) return(NULL)
    # Same dispatch as fit_bayes_suite: only tv = "week" needs the group-level fitter;
    # tv = "trend" is carried by the `week_idx` covariate in cov_spec.
    fit <- if (.tv_is_structural(tv))
             fit_bayes_renewal_tv(des, cov_spec = cov_spec, iter = iter, chains = 2L,
                                  link = link, tv = tv)
           else
             fit_bayes_renewal(des, cov_spec = cov_spec, iter = iter, chains = 2L, link = link)
    if (is.null(fit)) return(NULL)
    off <- tryCatch(bayes_forecast_offsets(zw, mobility_matrices, gt_pmfs, covariates,
             osrm_mat, zones_all, mob, gt, cov_spec, hz, rt_draws = rt),
             error = function(e) NULL)
    if (is.null(off)) return(NULL)
    if (.tv_is_structural(tv))
      predict_bayes_invasion_tv(fit, off, des, hz, affected_zones = character(0))
    else
      predict_bayes_invasion(fit, off, des, hz, affected_zones = character(0))
  }
}

#' Refit a Bayesian renewal model at each fold cutoff (causally, training only on
#' week <= cutoff) to trace the POSTERIOR import coefficient beta0 and covariate
#' hazard ratios over time, carrying proper posterior 90% credible intervals. (This replaced
#' the retired frequentist compute_beta_over_folds / compute_params_over_time traces, whose
#' per-fold Firth GLM refits drew figures FIGURE_KEEP then gated out.) A
#' covariate spec is used so BOTH beta0 and covariate-HR traces are produced.
compute_bayes_params_over_time <- function(zone_week_outbreak, cutoffs, mobility_matrices,
                                           gt_pmfs, covariates, osrm_mat, zones_all,
                                           mob = "M8", gt = "medium",
                                           cov_spec = c("log_pop", "ccvi", "d_min"),
                                           iter = 800L, nowcast_fn = NULL,
                                           linelist = NULL,
                                           # The AS-OF truncation. These origins are built by
                                           # .asof_train_slice()/reaggregate_asof(), so they
                                           # are the asof regime and must be nowcast with it —
                                           # the same object run_invasion_lfo() and
                                           # 33b_cascade_calibration.R receive. Without it this
                                           # path fell through to the onset->sample delay, so
                                           # the published beta0 trace was produced under a
                                           # DIFFERENT nowcast from the LFO that selected and
                                           # calibrated the models — the one cross-read a
                                           # reviewer would actually make. NULL warns.
                                           delay = NULL) {
  if (is.null(nowcast_fn)) nowcast_fn <- get0("apply_nowcast_correction")
  # As-of training slice, matching run_invasion_lfo(): see .asof_train_slice()
  # (20_forecast_detail.R) for why the final-count slice was a revision leak.
  .slice <- get0(".asof_train_slice")
  fit_one <- function(ci) {
    cut <- cutoffs[ci]            # index, NOT `for (cut in cutoffs)` — iterating a
                                  # Date vector unclasses cut to numeric (throws on R<4.3)
    zw  <- if (is.function(.slice)) .slice(zone_week_outbreak, cut, linelist, zones_all)
           else zone_week_outbreak %>% dplyr::filter(week_start <= cut)
    zwn <- tryCatch(nowcast_fn(zw, analysis_date = cut + 6, delay = delay),
                    error = function(e) {
                      # Named handler with a WARNING and a column guard, matching
                      # run_invasion_lfo() (16_invasion_eval.R). Returning `zw` bare was
                      # silent AND unsafe: if zw carries no `confirmed_nc`, .count_wide()
                      # (15_workhorse.R) returns an ALL-ZERO matrix, the design then has no
                      # at-risk rows, and the cutoff drops out of the published
                      # bayes_beta_over_folds / bayes_params_over_time traces with no
                      # diagnostic anywhere.
                      warning(sprintf("[bayes-folds] fold %s: nowcast failed (%s); using UNCORRECTED counts.",
                                      format(cut), conditionMessage(e)), call. = FALSE)
                      if (!"confirmed_nc" %in% names(zw) && "confirmed" %in% names(zw))
                        zw$confirmed_nc <- zw$confirmed
                      zw
                    })
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
