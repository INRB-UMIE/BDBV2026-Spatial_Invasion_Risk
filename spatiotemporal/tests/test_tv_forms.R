# =============================================================================
# test_tv_forms.R — the time-varying FORMS of the featured specification (2026-10-15)
#   * the Stan engine's R-side helpers (21c_frailty_models.R) against closed forms: the pointwise
#     log-likelihood, the training exposure, the week-1 forecast hazard, the drift draws
#   * the two-stage selection (16_invasion_eval.R, select_featured_two_stage())
#   * the forms' comparison with fixed beta (invasion_tv_form_comparison()): every form reported
#     with its bounds, burn-in excluded, the RNG stream restored
#   * the cascade (30-32): the hazard reader for both engines, and the frailty / drift limits of
#     the simulator
# Synthetic data only; no Stan or brms fit is run (fits are mocked as their draw matrices).
# =============================================================================
library(testthat)
.ST <- file.path(here::here(), "spatiotemporal")
suppressMessages(suppressWarnings({
  source(file.path(.ST, "00_config.R"))
  load_pipeline_modules("cascade", dir = .ST)        # data + model (incl. 21c) + cascade 30-32
}))

# A mocked stan_invasion_fit: the draw matrix and the stored inputs, as stan_invasion_fit() returns.
.fake_fit <- function(use_rw, use_frailty, nd = 3L, Z = 5L, S = 4L, K = 2L, sigma = 0.3, seed = 1L) {
  set.seed(seed)
  zones <- paste0("z", seq_len(Z)); feat <- c("ccvi", "d_min")[seq_len(K)]
  L  <- matrix(rexp(Z * S), Z, S, dimnames = list(zones, NULL))
  XA <- lapply(seq_len(S), function(s) matrix(rnorm(Z * K), Z, K, dimnames = list(zones, feat)))
  cols <- c("b0", if (K) paste0("beta[", seq_len(K), "]"), if (use_frailty) "lth[1]",
            if (use_rw) "sigma[1]", paste0("u[", seq_len(S), "]"))
  dr <- matrix(0, nd, length(cols), dimnames = list(NULL, cols))
  dr[, "b0"] <- rnorm(nd, -2, 0.3)
  if (K) dr[, paste0("beta[", seq_len(K), "]")] <- rnorm(nd * K, 0, 0.5)
  if (use_frailty) dr[, "lth[1]"] <- log(c(0.5, 1.2, 2.5)[seq_len(nd)])
  if (use_rw) {
    dr[, "sigma[1]"] <- sigma
    dr[, paste0("u[", seq_len(S), "]")] <- t(apply(matrix(rnorm(nd * S, 0, sigma), nd, S), 1, cumsum))
  }
  structure(list(draws = dr, feat = feat, center = stats::setNames(rep(0, K), feat),
                 scale = stats::setNames(rep(1, K), feat), S = S, L = L, XA = XA, zones = zones,
                 use_rw = use_rw, use_frailty = use_frailty, label = "fake"),
            class = "stan_invasion_fit")
}
.par_of <- function(fit, d) list(
  b0 = fit$draws[d, "b0"],
  beta = if (length(fit$feat)) fit$draws[d, paste0("beta[", seq_along(fit$feat), "]")] else numeric(0),
  lth = if (fit$use_frailty) fit$draws[d, "lth[1]"] else NULL,
  u = if (fit$use_rw) fit$draws[d, paste0("u[", seq_len(fit$S), "]")] else NULL)

test_that("the pointwise log-likelihood equals the reference likelihood, draw by draw", {
  for (cfg in list(c(FALSE, FALSE), c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE))) {
    fit <- .fake_fit(cfg[1], cfg[2])
    rows <- expand.grid(zone = seq_along(fit$zones), tr = seq_len(fit$S))
    X <- t(vapply(seq_len(nrow(rows)), function(n) fit$XA[[rows$tr[n]]][rows$zone[n], ], numeric(2)))
    y <- rbinom(nrow(rows), 1, 0.3)
    dd <- data.frame(.zone = fit$zones[rows$zone], .week = rows$tr + 1L,
                     logLam = log(fit$L[cbind(rows$zone, rows$tr)]), invaded = y, X)
    names(dd)[5:6] <- fit$feat
    ll <- stan_invasion_loglik(fit, dd)
    for (d in seq_len(nrow(fit$draws)))
      expect_equal(sum(ll[d, ]), frailty_loglik_R(.par_of(fit, d), y, X, exp(dd$logLam), rows$zone,
                                                  rows$tr, fit$L, fit$XA),
                   info = sprintf("rw=%s frailty=%s draw %d", cfg[1], cfg[2], d))
  }
})

test_that("the training exposure is each zone's cumulative hazard at its own weekly covariates", {
  fit <- .fake_fit(TRUE, TRUE)
  C <- .sf_cum_exposure(fit)
  for (d in seq_len(nrow(fit$draws))) {
    p <- .par_of(fit, d)
    ref <- Reduce(`+`, lapply(seq_len(fit$S), function(s)
      exp(p$b0 + p$u[s] + drop(fit$XA[[s]] %*% p$beta)) * fit$L[, s]))
    expect_equal(as.numeric(C[d, ]), as.numeric(ref))
  }
})

test_that("the week-1 forecast hazard is the frailty-conditional hazard; theta -> 0 gives fixed beta", {
  fit <- .fake_fit(TRUE, TRUE, sigma = 0)          # sigma 0: the forecast drift is exactly u_S
  # distinct, non-zero levels, so a forecast that used another week's level (or none) would fail
  fit$draws[, paste0("u[", seq_len(fit$S), "]")] <- matrix(c(0.1, -0.3, 0.2, 0.4, 0.5, -0.1,
                                                             -0.2, 0.3, 0.6, 0.7, -0.2, 0.9),
                                                           nrow(fit$draws), fit$S)
  off <- data.frame(health_zone = fit$zones, horizon = 1L, logLam = log(rexp(5)),
                    ccvi = rnorm(5), d_min = rnorm(5))
  dl <- stan_invasion_forecast_draws(fit, off, 1L)
  C <- .sf_cum_exposure(fit)
  for (d in seq_len(nrow(fit$draws))) {
    p <- .par_of(fit, d); th <- exp(p$lth)
    h <- exp(p$b0 + p$u[fit$S] + drop(as.matrix(off[, fit$feat]) %*% p$beta) + off$logLam)
    expect_equal(as.numeric(dl$h1$cum[d, ]), as.numeric(log1p(th * h / (1 + th * C[d, ])) / th))
  }
  fit0 <- fit; fit0$draws[, "lth[1]"] <- log(1e-10)
  fx <- fit; fx$use_frailty <- FALSE
  expect_equal(stan_invasion_forecast_draws(fit0, off, 1L)$h1$cum,
               stan_invasion_forecast_draws(fx, off, 1L)$h1$cum, tolerance = 1e-6)
  expect_identical(colnames(dl$h1$cum), fit$zones)
  expect_equal(dl$h1$p, 1 - exp(-dl$h1$cum))
})

test_that("the 2-week window: week-1 seeding at h1 * E[Z | history], joint survival, nesting", {
  fit <- .fake_fit(TRUE, TRUE, sigma = 0.4)
  Z <- length(fit$zones); ND <- nrow(fit$draws)
  set.seed(4)
  off <- rbind(data.frame(health_zone = fit$zones, horizon = 1L, logLam = log(rexp(Z)),
                          ccvi = rnorm(Z), d_min = rnorm(Z)),
               data.frame(health_zone = fit$zones, horizon = 2L, logLam = log(rexp(Z)),
                          ccvi = rnorm(Z), d_min = rnorm(Z)))
  Wm <- matrix(rexp(Z * Z), Z, Z, dimnames = list(fit$zones, fit$zones)); diag(Wm) <- 0
  Wm <- Wm / rowSums(Wm)
  Y <- matrix(0, Z, 3L, dimnames = list(fit$zones, NULL)); Y[1:2, ] <- c(2, 1, 3, 4, 5, 2)
  G <- c(0.5, 0.3, 0.2)
  attr(off, "projection") <- list(zones = fit$zones, W = Wm, G = G, Y = Y,
                                  rt_draws = stats::rgamma(400, 20, 16))
  dl <- stan_invasion_forecast_draws(fit, off, 1:2)
  # by hand
  th <- exp(fit$draws[, "lth[1]"]); Cf <- .sf_cum_exposure(fit); uf <- .sf_future_u(fit, 2L)
  xb <- function(h) { X <- as.matrix(off[off$horizon == h, fit$feat])
                      matrix(fit$draws[, "b0"], ND, Z) + fit$draws[, c("beta[1]", "beta[2]")] %*% t(X) }
  h1 <- exp(xb(1) + uf[, 1]) * matrix(exp(off$logLam[off$horizon == 1]), ND, Z, byrow = TRUE)
  zc <- 1 / (1 + th * Cf); zc[, rowSums(Y > 0) > 0] <- 1    # invaded zones: no selection on Z
  R  <- weekly_renewal_R_eff(.pair_rt_draws(attr(off, "projection")$rt_draws, ND, RANDOM_SEED), G)
  past <- function(t) Reduce(`+`, lapply(seq_along(G), function(k)
    if (t - k >= 1 && t - k <= ncol(Y)) G[k] * matrix(Y[, t - k], ND, Z, byrow = TRUE) else 0))
  Y1 <- R * past(ncol(Y) + 1) + h1 * zc
  src2 <- G[1] * Y1 + past(ncol(Y) + 2)
  h2 <- exp(xb(2) + uf[, 2]) * pmax(src2 %*% Wm, 1e-12)
  expect_equal(unname(dl$h1$cum), unname(log1p(th * h1 / (1 + th * Cf)) / th))
  expect_equal(unname(dl$h2$cum), unname(log1p(th * (h1 + h2) / (1 + th * Cf)) / th))
  expect_true(all(dl$h2$cum >= dl$h1$cum))
  # theta -> 0: the frailty path is the no-frailty path with the same draws
  f0 <- fit; f0$draws[, "lth[1]"] <- log(1e-10); fx <- fit; fx$use_frailty <- FALSE
  expect_equal(stan_invasion_forecast_draws(f0, off, 1:2)$h2$cum,
               stan_invasion_forecast_draws(fx, off, 1:2)$h2$cum, tolerance = 1e-6)
  # the 1-week innovations are a prefix of the 2-week ones
  expect_equal(as.numeric(.sf_future_u(fit, 1L)), uf[, 1])
})

test_that("forecast drift draws are reproducible and leave the caller's RNG stream untouched", {
  fit <- .fake_fit(TRUE, FALSE)
  set.seed(9); a <- runif(1); set.seed(9); u1 <- .sf_future_u(fit, 3L); b <- runif(1)
  expect_identical(a, b)
  expect_identical(.sf_future_u(fit, 3L), u1)
  expect_equal(dim(u1), c(nrow(fit$draws), 3L))
  # a one-week path is the last level plus one innovation of the draw's own SD
  set.seed(as.integer(RANDOM_SEED) + 7411L); e <- rnorm(nrow(fit$draws))
  expect_equal(as.numeric(.sf_future_u(fit, 1L)),
               as.numeric(fit$draws[, sprintf("u[%d]", fit$S)] + fit$draws[, "sigma[1]"] * e))
  expect_true(all(.sf_future_u(.fake_fit(FALSE, TRUE), 2L) == 0))
})

test_that("the parameter table reports the form's process parameters on their own scales", {
  fit <- .fake_fit(TRUE, TRUE); fit$diag <- data.frame(rhat_max = 1.001)
  tb <- stan_invasion_param_table(fit, "m")
  expect_setequal(tb$term, c("Intercept", "beta_last_training_week", "ccvi", "d_min",
                             "frailty_theta", "rw1_sigma_week"))
  expect_equal(tb$hr[tb$term == "frailty_theta"], stats::median(exp(fit$draws[, "lth[1]"])))
  # under RW1 the intercept is the walk's starting level, and says so; the level in force at the
  # last training week is its own row
  expect_identical(tb$effect_scale[tb$term == "Intercept"], RW1_START_LEVEL_SCALE)
  expect_equal(tb$hr[tb$term == "beta_last_training_week"],
               stats::median(exp(fit$draws[, "b0"] + fit$draws[, sprintf("u[%d]", fit$S)])))
  expect_identical(stan_invasion_param_table(.fake_fit(FALSE, TRUE), "m")$effect_scale[1], "hazard ratio")
  expect_true(all(is.na(tb$p_dir[tb$term %in% c("Intercept", "beta_last_training_week",
                                                "frailty_theta", "rw1_sigma_week")])))
  expect_true(all(tb$effect_scale[tb$term %in% c("ccvi", "d_min")] == "hazard ratio"))
  dr <- stan_invasion_posterior_draws(fit, "m")
  # The coefficients AND the form's process parameters (FigureS13 draws one panel per parameter
  # of the featured model), each with the prior it was fitted under: theta = exp(lth) with
  # lth ~ TV_FORM_PRIOR_LOG_THETA, exp(sigma) with sigma ~ the half-normal TV_FORM_PRIOR_RW_SD.
  expect_setequal(unique(dr$term), c("Intercept", "ccvi", "d_min", "frailty_theta", "rw1_sigma_week"))
  expect_identical(unique(dr$prior[dr$is_intercept]), BAYES_PRIOR_INTERCEPT)
  expect_identical(unique(dr$prior[dr$term == "frailty_theta"]), TV_FORM_PRIOR_LOG_THETA)
  expect_identical(unique(dr$prior[dr$term == "rw1_sigma_week"]), paste0("half-", trimws(TV_FORM_PRIOR_RW_SD)))
  expect_equal(stats::median(dr$hr[dr$term == "frailty_theta"]),
               tb$hr[tb$term == "frailty_theta"])
})

# ---- the two-stage selection -------------------------------------------------------------------
.eval_row <- function(m, h, sk, rot, ls) data.frame(
  method = m, horizon = h, auc_pr_skill = sk, mean_rank_of_truth = rot, log_score = ls,
  log_score_recal = ls, sel_auc_pr_skill_date = sk, sel_mean_rank_of_truth = rot,
  sel_log_score = ls, sel_log_score_recal = ls, partial_cv = FALSE, prob_calibrated = TRUE,
  stringsAsFactors = FALSE)
.mk_eval <- function(form_better, best = 2L) {
  base <- list(c("Bayes-M14-fill-geo", 30, 10, 0.020), c("Bayes-M4-geo", 25, 12, 0.021),
               c("Bayes-M8-fill-geo", 20, 14, 0.022))
  fs <- paste0("Bayes-M14-fill-geo-", names(TV_FORMS))
  rows <- list()
  for (h in 1:2) {
    for (b in base) rows[[length(rows) + 1L]] <- .eval_row(b[1], h, as.numeric(b[2]), as.numeric(b[3]), as.numeric(b[4]))
    for (i in seq_along(fs))
      rows[[length(rows) + 1L]] <- .eval_row(fs[i], h, if (form_better && i == best) 40 else 29 - i,
                                             if (form_better && i == best) 8 else 11 + i,
                                             if (form_better && i == best) 0.018 else 0.0205 + i / 1000)
    # a sensitivity arm that would win everything if it were eligible
    rows[[length(rows) + 1L]] <- .eval_row("Bayes-M14-fill-geo-gtshort", h, 99, 1, 0.001)
  }
  do.call(rbind, rows)
}

test_that("stage 1 ignores the forms and the arms; stage 2 features the best of the pool", {
  # the frailty forms are secondary: scored, never candidates (TV_FORMS `selectable`)
  expect_identical(TV_FORM_SELECTABLE_SUFFIXES, "tvrw1")
  expect_setequal(TV_FORM_SECONDARY_SUFFIXES, c("tvfrailty", "tvrw1frailty"))
  i_rw1 <- match("tvrw1", names(TV_FORMS))
  e <- .mk_eval(form_better = TRUE, best = i_rw1)
  tw <- suppressMessages(select_featured_two_stage(e))
  expect_identical(tw$stage1, "Bayes-M14-fill-geo")
  expect_identical(tw$forms, "Bayes-M14-fill-geo-tvrw1")
  expect_setequal(tw$secondary, paste0("Bayes-M14-fill-geo-", TV_FORM_SECONDARY_SUFFIXES))
  expect_identical(tw$featured, "Bayes-M14-fill-geo-tvrw1")
  expect_false(any(grepl("gtshort|frailty", tw$table$method)))
  # a frailty form that beats everything on every axis is still never featured
  for (sfx in TV_FORM_SECONDARY_SUFFIXES) {
    ef <- .mk_eval(form_better = TRUE, best = match(sfx, names(TV_FORMS)))
    twf <- suppressMessages(select_featured_two_stage(ef))
    expect_identical(twf$featured, "Bayes-M14-fill-geo")
    expect_false(sfx %in% tv_form_of(twf$pool))
  }
  e2 <- .mk_eval(form_better = FALSE)
  tw2 <- suppressMessages(select_featured_two_stage(e2))
  expect_identical(tw2$featured, "Bayes-M14-fill-geo")     # fixed beta kept when no form is better
  # without cross-validated forms the stage-1 pick is featured
  e3 <- e[is.na(tv_form_of(e$method)), ]
  tw3 <- suppressMessages(select_featured_two_stage(e3))
  expect_identical(tw3$featured, tw3$stage1); expect_length(tw3$forms, 0L)
  # stage 1 is exactly the single-stage selector with the default exclusion
  expect_identical(tw$stage1, suppressMessages(best_invasion_model(e)))
})

# ---- the forms against fixed beta -------------------------------------------------------------
test_that("the comparison reports every form with its bounds, excludes burn-in, restores the RNG", {
  set.seed(21)
  cut <- as.Date("2026-07-04") + 7L * (0:9)
  base <- expand.grid(health_zone = paste0("z", 1:60), cutoff = cut, horizon = 1:2,
                      stringsAsFactors = FALSE)
  base$fold_id <- as.integer(factor(base$cutoff))
  base$is_new_invasion <- rbinom(nrow(base), 1, 0.03)
  base$was_active_before <- FALSE; base$delta_preq <- 1
  meths <- c("Bayes-M14-fill-geo", paste0("Bayes-M14-fill-geo-", names(TV_FORMS)))
  lfo <- do.call(rbind, lapply(meths, function(m) {
    d <- base; d$method <- m
    d$p_invasion <- plogis(qlogis(0.03) + rnorm(nrow(d), 0, 0.4) + 1.2 * d$is_new_invasion)
    d$p_recal <- d$p_invasion; d }))
  cmp <- invasion_tv_form_comparison(lfo, "Bayes-M14-fill-geo", B = 200L)
  expect_setequal(unique(cmp$guard$form), vapply(TV_FORMS, `[[`, "", "form"))
  expect_setequal(unique(cmp$metrics$form), c("fixed", vapply(TV_FORMS, `[[`, "", "form")))
  expect_true(all(c("zone_lo_bonferroni", "origin_lo_bonferroni", "zone_lo90", "zone_hi90",
                    "origin_lo90", "origin_hi90", "guard_passed", "dm_lag") %in% names(cmp$guard)))
  expect_equal(cmp$bonferroni_level, 0.05 / length(TV_FORMS))
  # burn-in rounds are excluded from the post-burn-in statistics
  expect_equal(cmp$n_rounds, sum(!invasion_is_burnin(cut)))
  # the RNG stream is restored
  set.seed(5); a <- runif(1); set.seed(5); invisible(invasion_tv_form_comparison(lfo, "Bayes-M14-fill-geo", B = 50L))
  expect_identical(runif(1), a)
})

# ---- the cascade ------------------------------------------------------------------------------
test_that("the cascade's hazard reader returns the form's quantities from a Stan fit", {
  skip_if_not(exists(".cascade_hazard_draws", mode = "function"))
  fit <- .fake_fit(TRUE, TRUE)
  hz <- .cascade_hazard_draws(fit, fit$zones)
  expect_true(hz$use_rw && hz$use_frailty)
  expect_equal(hz$b_int, fit$draws[, "b0"])
  expect_identical(colnames(hz$slopes), fit$feat)
  expect_equal(hz$uS, fit$draws[, sprintf("u[%d]", fit$S)])
  expect_equal(hz$theta, exp(fit$draws[, "lth[1]"]))
  expect_equal(unname(hz$C0), unname(.sf_cum_exposure(fit)))
  expect_error(.cascade_hazard_draws(fit, rev(fit$zones)), "zone spine")
})

# The 4-zone metapopulation of test_cascade_fixes.R: A and B infected, C and D at risk.
.mk_prep <- function(b0 = -3) {
  z <- c("A", "B", "C", "D"); nz <- 4L; nd <- 10L
  W <- matrix(0, nz, nz, dimnames = list(z, z)); W["B", "C"] <- 1; W["B", "D"] <- 1
  OS <- matrix(Inf, nz, nz, dimnames = list(z, z))
  Y0 <- matrix(c(0, 2, 3,  0, 2, 3,  0, 0, 0,  0, 0, 0), nz, 3L, byrow = TRUE, dimnames = list(z, NULL))
  list(zones_all = z, nz = nz, W = W, Wt = t(W), denom_fi = colSums(W),
       G = c(0.6, 0.4), kmax = 2L, G_list = NULL, gt_weights = NULL, gt_keys = NULL,
       Y0 = Y0, W0 = 3L, affected0 = c(TRUE, TRUE, FALSE, FALSE),
       Z_static = matrix(0, nz, 0, dimnames = list(z, NULL)), OS = OS, has_dmin = FALSE,
       dmin_center = NULL, dmin_scale = NULL, DMIN_SENTINEL = 1e4,
       b_int = rep(b0, nd), G_static = matrix(0, nd, 0), g_dmin = rep(0, nd), n_draws = nd,
       use_rw = FALSE, use_frailty = FALSE, uS = rep(0, nd), sigma_rw = NULL, theta = NULL, C0 = NULL,
       reff = list(R_zone = stats::setNames(rep(1.2, nz), z), R_draws = rep(1.2, 50),
                   r_min = 0.05, r_max = 10, prov = rep("P", nz)),
       prov = rep("P", nz), kernel = "K", gt = "medium")
}
.scen <- list(label = "test", c_fun = function(t) rep(1, length(t)))
.sim <- function(p) simulate_cascade(p, .scen, n_mc = 80L, horizon = 5L, delta = 1, psi = 0,
                                     n_rep = 4L, rw_sigma = 0, seed = 17L)

test_that("frailty: theta -> 0 with no past exposure reproduces fixed beta; heavy exposure stops invasion", {
  skip_if_not(exists("simulate_cascade", mode = "function"))
  p0 <- .mk_prep(b0 = -1)
  base <- .sim(p0)
  pf <- p0; pf$use_frailty <- TRUE; pf$theta <- rep(1e-12, pf$n_draws)
  pf$C0 <- matrix(0, pf$n_draws, pf$nz, dimnames = list(NULL, pf$zones_all))
  expect_identical(.sim(pf)$tau, base$tau)
  expect_identical(.sim(pf)$knobs$form, "frailty")
  # Exposure without invasion is evidence of resistance: with a large frailty variance and a huge
  # training exposure the conditional hazard of the never-invaded zones collapses.
  ph <- pf; ph$theta <- rep(5, ph$n_draws); ph$C0[, c("C", "D")] <- 1e8
  expect_true(all(is.na(.sim(ph)$tau[c("C", "D"), ])))
  expect_gt(sum(!is.na(base$tau[c("C", "D"), ])), 0L)
})

test_that("drift: a zero innovation SD at level u_S is fixed beta with b0 + u_S", {
  skip_if_not(exists("simulate_cascade", mode = "function"))
  p0 <- .mk_prep(b0 = -1)
  pr <- p0; pr$use_rw <- TRUE; pr$uS <- rep(0.7, pr$n_draws); pr$sigma_rw <- rep(0, pr$n_draws)
  p1 <- .mk_prep(b0 = -1 + 0.7)
  expect_identical(.sim(pr)$tau, .sim(p1)$tau)
  expect_identical(.sim(pr)$knobs$form, "rw1")
  # with a positive SD the run is still reproducible from its seed (paired contrasts stay paired)
  pr$sigma_rw <- rep(0.5, pr$n_draws)
  expect_identical(.sim(pr)$tau, .sim(pr)$tau)
})
