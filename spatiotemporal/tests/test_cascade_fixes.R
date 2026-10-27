# =============================================================================
# test_cascade_fixes.R - regression tests for the 13-week cascade defects fixed on 2026-10-15
# =============================================================================
# Each block pins a PROPERTY that was violated, on a synthetic metapopulation small enough to
# reason about by hand, so a later edit that re-introduces the defect fails here rather than in a
# published table:
#   * a gateway knockout bars the IMPORT FORCE only, never the frontier (d_min)        [F14]
#   * a conditional "hub" that is already infected is refused, not silently a baseline  [F4]
#   * every covariate in the hazard is simulated, or cascade_prepare() stops           [F15]
#   * the R walk is ONE national walk per (week, parameter group), not per zone        [F5]
#   * a reach interval at the Monte-Carlo floor is flagged as such                     [floor]
#   * the deliberately-uncalibrated path keeps delta's SE                              [F2]
#   * settled-label age and the gate's ranking-only pass                                [F9/F10]
#   * the week-1 refit is stored and reported against EACH reference, by name          [R1]
#   * downstream readers use only the canonical, run-stamped layer cache               [R2]
#   * a floored interval is at most 1/D wide (narrower where clipped)                  [R3]
#   * the calibration JSON key, the rolling-scale and origin-training wording           [R4]
# The cascade layer (30-32) and 33/33b/34/36 are sourced here when the runner has not loaded them
# (a runner that loads only the "model" layer has none of them); every block skips cleanly when
# its function is unavailable. None of these files writes anything when sourced.
# =============================================================================

.ST <- file.path(here::here(), "spatiotemporal")
.need <- function(fn) if (!exists(fn, mode = "function"))
  skip(paste("Function", fn, "not available - source the cascade modules first"))
.src_if_missing <- function(fn, file)
  if (!exists(fn, mode = "function") && file.exists(file.path(.ST, file)))
    suppressMessages(suppressWarnings(sys.source(file.path(.ST, file), envir = globalenv())))
# The cascade layer, in its load order (30 configures 31/32), only if it is absent and the model
# layer it builds on is present.
if (!exists("simulate_cascade", mode = "function") && exists("fit_bayes_renewal", mode = "function"))
  for (.f in c("30_projection_config.R", "31_source_dynamics.R", "32_cascade_simulator.R"))
    if (file.exists(file.path(.ST, .f)))
      suppressMessages(suppressWarnings(sys.source(file.path(.ST, .f), envir = globalenv())))
.src_if_missing("cascade_conditional_seed", "34_conditional_queries.R")
.src_if_missing("cascade_consistency_gate", "33_cascade_eval.R")
.src_if_missing("cascade_calibration_origins", "33b_cascade_calibration.R")
.src_if_missing(".cascade_h1_refit_ratios", "36_report.R")

# ---- a 4-zone synthetic metapopulation ------------------------------------------------------
#   A  infected, exports to NOBODY (its row of W is zero), 10 min from C
#   B  infected, exports to C and D, 1000 min from both
#   C  at risk, close to A           D  at risk, far from everything
# With a strongly negative d_min slope, C's hazard depends heavily on whether A counts as part of
# the infected front.
.mk_prep <- function(g_dmin = -2, b0 = -3) {
  z <- c("A", "B", "C", "D"); nz <- 4L; nd <- 10L
  W <- matrix(0, nz, nz, dimnames = list(z, z)); W["B", "C"] <- 1; W["B", "D"] <- 1
  OS <- matrix(Inf, nz, nz, dimnames = list(z, z))
  OS["C", "A"] <- OS["A", "C"] <- 10
  OS["C", "B"] <- OS["B", "C"] <- OS["D", "A"] <- OS["A", "D"] <- OS["D", "B"] <- OS["B", "D"] <- 1000
  OS["A", "B"] <- OS["B", "A"] <- 500; OS["C", "D"] <- OS["D", "C"] <- 1000
  Y0 <- matrix(c(0, 2, 3,  0, 2, 3,  0, 0, 0,  0, 0, 0), nz, 3L, byrow = TRUE,
               dimnames = list(z, NULL))
  list(zones_all = z, nz = nz, W = W, Wt = t(W), denom_fi = colSums(W),
       G = c(0.6, 0.4), kmax = 2L, G_list = NULL, gt_weights = NULL, gt_keys = NULL,
       Y0 = Y0, W0 = 3L, affected0 = c(TRUE, TRUE, FALSE, FALSE),
       Z_static = matrix(0, nz, 0, dimnames = list(z, NULL)), OS = OS, has_dmin = TRUE,
       dmin_center = log1p(100), dmin_scale = 1, DMIN_SENTINEL = 1e4,
       b_int = rep(b0, nd), G_static = matrix(0, nd, 0), g_dmin = rep(g_dmin, nd),
       n_draws = nd,
       reff = list(R_zone = stats::setNames(rep(1.2, nz), z), R_draws = rep(1.2, 50),
                   r_min = 0.05, r_max = 10, prov = rep("P", nz)),
       prov = rep("P", nz), kernel = "K", gt = "medium")
}
.scen <- list(label = "test", c_fun = function(t) rep(1, length(t)))

test_that("a knockout of a zone that exports nothing changes nothing (frontier kept) [F14]", {
  .need("simulate_cascade")
  p <- .mk_prep()
  base <- simulate_cascade(p, .scen, n_mc = 80L, horizon = 4L, delta = 1, psi = 0, n_rep = 4L,
                           rw_sigma = 0, seed = 11L)
  ko   <- simulate_cascade(p, .scen, n_mc = 80L, horizon = 4L, delta = 1, psi = 0, n_rep = 4L,
                           rw_sigma = 0, seed = 11L, bar_sources = "A")
  # A seeds nobody, so barring its import force is a no-op. Under common random numbers the two
  # runs must therefore be IDENTICAL. Removing A from the d_min front (the old behaviour) moved
  # C from 10 to 1000 minutes from the nearest infected zone and cut its hazard ~1e4-fold.
  expect_identical(ko$tau, base$tau)
  expect_identical(ko$new_by_week, base$new_by_week)
  # C is near-certain to be reached in week 1 because A is 10 min away: the front matters.
  expect_gt(mean(!is.na(base$tau["C", ]) & base$tau["C", ] == 1L), 0.9)
})

test_that("a knockout still bars the import force of a zone that DOES export [F14]", {
  .need("simulate_cascade")
  p <- .mk_prep(g_dmin = 0)                 # no frontier term: only the import force matters
  base <- simulate_cascade(p, .scen, n_mc = 80L, horizon = 4L, delta = 1, psi = 0, n_rep = 4L,
                           rw_sigma = 0, seed = 11L)
  ko   <- simulate_cascade(p, .scen, n_mc = 80L, horizon = 4L, delta = 1, psi = 0, n_rep = 4L,
                           rw_sigma = 0, seed = 11L, bar_sources = "B")
  # B is the only exporter: barred, no at-risk zone has any import force, so nothing is seeded.
  expect_true(all(is.na(ko$tau[c("C", "D"), ])))
  expect_gt(sum(!is.na(base$tau[c("C", "D"), ])), 0L)
})

test_that("conditional seeding refuses a hub that is already infected [F4]", {
  .need("cascade_conditional_seed")
  p <- .mk_prep()
  # A is infected at the origin: forcing it is a no-op, so the experiment would be the baseline.
  expect_error(cascade_conditional_seed(p, .scen, "A", delta = 1, psi = 0, n_mc = 8L,
                                        horizons = 1:2),
               "already infected")
  expect_error(cascade_conditional_seed(p, .scen, "NOT_A_ZONE", delta = 1, psi = 0, n_mc = 8L,
                                        horizons = 1:2),
               "not on the health-zone spine")
  # An at-risk hub IS forced: reached in week 1 in every iteration.
  cs <- cascade_conditional_seed(p, .scen, "D", delta = 1, psi = 0, n_mc = 8L, horizons = 1:2)
  expect_true(all(cs$sim$tau["D", ] == 1L))
})

test_that("the R walk is one national walk per (week, parameter group) [F5]", {
  .need("simulate_cascade")
  body_txt <- paste(deparse(body(simulate_cascade)), collapse = "\n")
  # One innovation per (week, group): a horizon x D matrix, indexed [w, g] ...
  expect_true(grepl("matrix(stats::rnorm(horizon * D), horizon, D)", body_txt, fixed = TRUE))
  expect_true(grepl("Z_Rw_grp[w, g]", body_txt, fixed = TRUE))
  # ... and NOT the per-zone nz x horizon x D array it replaced.
  expect_false(grepl("rnorm(nz * horizon * D)", body_txt, fixed = TRUE))
  # Behaviour: with a walk switched on the simulator runs, and identical seeds stay paired.
  p <- .mk_prep(g_dmin = 0)
  a <- simulate_cascade(p, .scen, n_mc = 40L, horizon = 6L, delta = 1, psi = 0, n_rep = 4L,
                        rw_sigma = 0.3, seed = 5L)
  b <- simulate_cascade(p, .scen, n_mc = 40L, horizon = 6L, delta = 1, psi = 0, n_rep = 4L,
                        rw_sigma = 0.3, seed = 5L)
  expect_identical(a$tau, b$tau)
  expect_identical(a$cases_by_week, b$cases_by_week)
})

test_that("cascade_prepare simulates every hazard covariate or stops [F15]", {
  .need("cascade_prepare")
  z <- c("A", "B", "C", "D")
  W <- matrix(0.1, 4, 4, dimnames = list(z, z))
  layer <- list(zones_all = z, mobility_matrices = list(K = W),
                zone_week_nc = data.frame(health_zone = rep(z, 2),
                                          week_start = rep(as.Date(c("2026-01-02", "2026-01-09")), each = 4),
                                          confirmed_nc = c(1, 0, 0, 0, 2, 1, 0, 0)),
                covariates = data.frame(nom = z, ccvi = c(0.1, 0.4, 0.6, 0.9),
                                        log_pop = log(c(1e5, 2e5, 3e5, 4e5))),
                osrm_mat = matrix(60, 4, 4, dimnames = list(z, z)), gt_pmfs = NULL)
  design <- list(center = c(ccvi = 0.5, d_min = 1), scale = c(ccvi = 0.3, d_min = 1))
  fit <- data.frame(b_Intercept = rep(-3, 6), b_ccvi = rep(0.5, 6), b_d_min = rep(-0.2, 6))
  reff <- list(prov = rep("P", 4))
  pr <- cascade_prepare(layer, fit, design, reff, kernel = "K", gt = "medium",
                        cov_spec = c("ccvi", "d_min"))
  expect_equal(colnames(pr$Z_static), "ccvi")
  expect_equal(unname(pr$Z_static[, "ccvi"]), (c(0.1, 0.4, 0.6, 0.9) - 0.5) / 0.3)
  expect_equal(unname(pr$G_static[, "ccvi"]), rep(0.5, 6))
  # A covariate the simulator has no trajectory for is refused, not silently dropped.
  expect_error(cascade_prepare(layer, fit, design, reff, kernel = "K", gt = "medium",
                               cov_spec = c("ccvi", "d_min", "alert_import")),
               "cannot carry them forward")
  # A fitted slope the cov_spec would not simulate is refused too.
  fit2 <- cbind(fit, b_log_pop = rep(-0.3, 6))
  expect_error(cascade_prepare(layer, fit2, design, reff, kernel = "K", gt = "medium",
                               cov_spec = c("ccvi", "d_min")),
               "does not simulate")
})

test_that("reach intervals at the Monte-Carlo floor are flagged", {
  .need("cascade_reach_table")
  # D = 4 draw groups x R = 4 replicates; zone 1 identical in every iteration of every group
  # except one hit -> no resolvable posterior spread -> floored. Zone 2 varies strongly by group.
  D <- 4L; R <- 4L; M <- D * R
  tau <- matrix(NA_integer_, 2L, M)
  tau[1, 1] <- 1L                                    # one hit, in group 1
  tau[2, rep(c(TRUE, FALSE), each = R * 2)] <- 1L    # groups 1-2 always, 3-4 never
  sim <- list(zones_all = c("Z1", "Z2"), n_mc = M, n_rep = R, draw_group = rep(seq_len(D), each = R),
              tau = tau, est_tau = tau, affected0 = c(FALSE, FALSE))
  rt <- cascade_reach_table(sim, horizons = 1L)
  expect_true("p_ci_at_mc_floor" %in% names(rt))
  # Zone 1: the single hit is exactly Bernoulli within its group, so the posterior variance is 0
  # and the interval is the floor (p +- 1/(2D), clipped at 0) — flagged, not an estimate.
  expect_true(rt$p_ci_at_mc_floor[1])
  expect_equal(rt$p_hi[1], rt$p_invasion[1] + 1 / (2 * D))
  # [R3] p = 1/16 < 1/(2D), so the floor is clipped at 0: its width is p + 1/(2D) = 0.1875, LESS
  # than 1/D = 0.25. "Width exactly 1/D" is false for such a row; "at most 1/D" is the property.
  expect_equal(rt$p_lo[1], 0)
  expect_lt(rt$p_hi[1] - rt$p_lo[1], 1 / D)
  expect_equal(rt$p_hi[1] - rt$p_lo[1], rt$p_invasion[1] + 1 / (2 * D))
  # Zone 2: a genuine between-group spread; its interval is estimated, not floored.
  expect_false(rt$p_ci_at_mc_floor[2])
  expect_gt(rt$p_hi[2] - rt$p_lo[2], 1 / D)
})

test_that("the deliberately-uncalibrated path keeps delta's SE [F2]", {
  .need("cascade_short_horizon_delta_info")
  d <- structure(0.6, estimator = "maximum_likelihood_16b", se_log = 0.23, lo = 0.42, hi = 0.91,
                 source = "Bayes-M14-fill-geo", cal_in_large = 1.5, boundary_hit = FALSE)
  di <- cascade_short_horizon_delta_info(d)
  expect_equal(di$se_log, 0.23)
  expect_equal(di$delta, 0.6)
  expect_identical(di$estimator, "maximum_likelihood_16b")
})

test_that("calibration and backtest use SETTLED labels; the gate passes on the per-zone pattern only [F9, F10]", {
  .need("cascade_calibration_origins"); .need("cascade_backtest"); .need("cascade_consistency_gate")
  settled <- get0("LFO_RELIABLE_EVAL_AGE_DAYS", ifnotfound = 25L)
  expect_equal(eval(formals(cascade_calibration_origins)$min_eval_age_days), settled)
  expect_equal(eval(formals(cascade_backtest)$min_eval_age_days), settled)
  # The level ratio is ~delta for any delta (both sides share one fit), so it must not gate.
  gb <- paste(deparse(body(cascade_consistency_gate)), collapse = "\n")
  expect_true(grepl("pass <- pass_rank", gb, fixed = TRUE))
  expect_false(grepl("pass <- pass_rank && pass_level", gb, fixed = TRUE))
  # The pass rule is the Monte-Carlo-calibrated dispersion test, not a Spearman threshold.
  expect_true(grepl("is.finite(p_disp) && p_disp >= gate_alpha", gb, fixed = TRUE))
})

test_that("run_cascade.R makes a failed calibration fatal and stamps its outputs [F1, F8]", {
  src <- paste(readLines(file.path(.ST, "run_cascade.R"), warn = FALSE), collapse = "\n")
  code <- paste(sub("#.*$", "", readLines(file.path(.ST, "run_cascade.R"), warn = FALSE)),
                collapse = "\n")
  # the old warn-and-continue handler is gone; failure routes to .calib_fail(), which stops
  expect_false(grepl("falling back to delta-only", code, fixed = TRUE))
  expect_true(grepl("error = function(e) .calib_fail(conditionMessage(e))", code, fixed = TRUE))
  expect_true(grepl("calibrated = CALIBRATED", code, fixed = TRUE))
  # every run is stamped, and the layer cache is always rewritten
  expect_true(grepl(".write_stamp(\"complete\")", code, fixed = TRUE))
  expect_true(grepl("layer_cache_v2.rds", code, fixed = TRUE))
  # 40 is launched as its own Rscript, not source()d (a guarded script is a no-op when sourced)
  expect_false(grepl("source(file.path(ST_DIR, \"40_cascade_next_dominoes.R\"))", code, fixed = TRUE))
})

# ---- [R1] the week-1 refit against EACH reference ------------------------------------------
# cascade_calibrate() with its collaborators stubbed (the stubs live in an environment that sits
# between the function and the global one, so nothing global is touched). 16b's short-horizon
# factor is 0.5, the week-1 refit 0.8, the K-week fit 0.25: the two ratios are 1.6 and 3.2, and
# the old code stored only the first while the report called it the ratio to the K-week factor.
.calib_stubbed <- function(k_week_delta) {
  e <- new.env(parent = globalenv())
  e$cascade_fit_delta <- function(...) structure(0.5, estimator = "maximum_likelihood_16b",
                                                 se_log = 0.2, lo = 0.4, hi = 0.62)
  e$cascade_calibration_origins <- function(...) list(list(
    cutoff_date = as.Date("2026-07-03"), Ccol = 9L, R_nat = 1.1, n_with_cases = 4L,
    n_prior_dominated = 0L, import_share = 0.1, n_atrisk = 20L, complete_horizons = c(1L, 6L),
    n_events = c(`1` = 3L, `6` = 7L), training_source = "recensor", vintage_date = as.Date(NA),
    rolling_predictors = TRUE, truth = list(`1` = integer(20), `6` = integer(20))))
  e$cascade_delta_h1_refit <- function(...) list(delta = 0.8, lo = 0.6, hi = 1.1, n_events = 12L)
  e$cascade_fit_delta_oos <- function(...)
    if (is.null(k_week_delta)) NULL else
      list(delta = k_week_delta, converged = TRUE, boundary_hit = FALSE,
           bracket_lo = k_week_delta * 0.99, bracket_hi = k_week_delta * 1.01, K = 6L)
  e$cascade_calibration_report <- function(...) NULL
  f <- cascade_calibrate; environment(f) <- e
  w <- character(0)
  out <- withCallingHandlers(
    f(layer = list(ll = NULL), K = 6L, kernel = "K", gt = "medium", seed = 1L),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  list(out = out, warnings = w)
}

test_that("33b stores the week-1 refit ratio to the short-horizon AND to the deployed factor [R1]", {
  .need("cascade_calibrate"); .need("cascade_short_horizon_delta_info")
  r <- .calib_stubbed(k_week_delta = 0.25)
  di <- r$out$delta_info
  expect_equal(as.numeric(r$out$delta), 0.25)                 # the K-week factor is deployed
  expect_equal(di$h1_refit_delta, 0.8)
  expect_equal(di[["h1_refit_ratio_short"]], 0.8 / 0.5)        # against 16b's factor
  expect_equal(di[["h1_refit_ratio_deployed"]], 0.8 / 0.25)    # against the K-week factor
  # the ambiguous single ratio is no longer written
  expect_null(di[["h1_refit_ratio"]])
  # the transfer warning no longer claims the short-horizon factor is the one used
  expect_false(any(grepl("short-horizon factor is used", r$warnings, fixed = TRUE)))
  expect_true(any(grepl("Neither value is the one projected with", r$warnings, fixed = TRUE)))
})

test_that("when the K-week fit fails the deployed ratio is the short-horizon one [R1]", {
  .need("cascade_calibrate"); .need("cascade_short_horizon_delta_info")
  r <- .calib_stubbed(k_week_delta = NULL)
  di <- r$out$delta_info
  expect_null(r$out$delta_oos)
  expect_equal(as.numeric(r$out$delta), 0.5)                  # short-horizon factor carried
  expect_equal(di[["h1_refit_ratio_short"]], 1.6)
  expect_equal(di[["h1_refit_ratio_deployed"]], 1.6)
})

test_that("the report reads each week-1 ratio for its own reference, old contexts included [R1]", {
  .need(".cascade_h1_refit_ratios"); .need(".cascade_h1_ratio_phrase")
  # new context
  r <- .cascade_h1_refit_ratios(list(h1_refit_delta = 0.8, h1_refit_ratio_short = 1.6,
                                     h1_refit_ratio_deployed = 3.2), delta_deployed = 0.25)
  expect_equal(r$ratio_short, 1.6); expect_equal(r$ratio_deployed, 3.2)
  ph <- .cascade_h1_ratio_phrase(r, "the deployed K-week factor")
  expect_identical(ph, "ratio 3.20 to the deployed K-week factor; 1.60 to 16b's short-horizon factor")
  # LEGACY context: only `h1_refit_ratio`, which 33b always computed against the short-horizon
  # factor. It must be labelled as that, and the deployed ratio recomputed — never 1.6 "to the
  # K-week factor", which is what the report printed.
  r0 <- .cascade_h1_refit_ratios(list(h1_refit_delta = 0.8, h1_refit_ratio = 1.6),
                                 delta_deployed = 0.25)
  expect_equal(r0$ratio_short, 1.6); expect_equal(r0$ratio_deployed, 3.2)
  ph0 <- .cascade_h1_ratio_phrase(r0, "the deployed delta")
  expect_false(grepl("1.60 to the deployed", ph0, fixed = TRUE))
  expect_true(grepl("1.60 to 16b's short-horizon factor", ph0, fixed = TRUE))
  # Partial-matching trap: a context holding only the DEPLOYED ratio must not have it read back
  # as the short one (`ci$h1_refit_ratio` would partial-match it).
  r1 <- .cascade_h1_refit_ratios(list(h1_refit_delta = 0.8, h1_refit_ratio_deployed = 3.2),
                                 delta_deployed = 0.25)
  expect_true(is.na(r1$ratio_short)); expect_equal(r1$ratio_deployed, 3.2)
  # no refit at all -> nothing to quote
  rn <- .cascade_h1_refit_ratios(list(), delta_deployed = 0.25)
  expect_true(is.na(rn$delta)); expect_identical(.cascade_h1_ratio_phrase(rn), "")
  # both report sentences go through the helper, not the raw field
  wb <- paste(deparse(body(cascade_write_report)), collapse = "\n")
  expect_false(grepl("h1_refit_ratio)", wb, fixed = TRUE))
  expect_true(grepl("deployed_label = \"the deployed K-week factor\"", wb, fixed = TRUE))
  expect_true(grepl("deployed_label = \"the deployed delta\"", wb, fixed = TRUE))
})

# ---- [R4] origin training source, rolling scale, JSON key -------------------------------------
test_that("the report states what each calibration origin was trained on [R4f]", {
  .need(".cascade_origin_training_txt")
  od <- function(src, hz = rep("1/6", length(src)))
    list(origin_diagnostics = data.frame(cutoff_date = as.Date("2026-07-03") + 7 * seq_along(src),
                                         horizons = hz, training_source = src,
                                         stringsAsFactors = FALSE))
  t_ok <- .cascade_origin_training_txt(od(c("recensor", "vintage", "recensor")), K = 6L)
  expect_true(grepl("trained only on counts reconstructed as of that origin", t_ok, fixed = TRUE))
  expect_true(grepl("1 from the archived snapshot", t_ok, fixed = TRUE))
  expect_true(grepl("2 from the final line list re-censored", t_ok, fixed = TRUE))
  t_fs <- .cascade_origin_training_txt(od(c("recensor", "final_slice", "vintage")), K = 6L)
  expect_false(grepl("trained only on", t_fs, fixed = TRUE))
  expect_true(grepl("**1 (2026-07-17) was trained on the final-count slice", t_fs, fixed = TRUE))
  # an origin without a settled K-week outcome did not enter delta's fit and is not counted
  t_k <- .cascade_origin_training_txt(od(c("recensor", "final_slice"), hz = c("1/6", "1")), K = 6L)
  expect_true(grepl("trained only on counts reconstructed", t_k, fixed = TRUE))
  # nothing recorded -> say so, never assert
  expect_true(grepl("not recorded", .cascade_origin_training_txt(list(), K = 6L), fixed = TRUE))
  expect_true(grepl("not recorded",
                    .cascade_origin_training_txt(od(c(NA_character_, NA_character_)), K = 6L),
                    fixed = TRUE))
  # the origins themselves are the fallback when the diagnostics table is absent
  os <- list(origins = list(list(truth = list(`1` = 0L, `6` = 0L), training_source = "final_slice",
                                 cutoff_date = as.Date("2026-07-10"))))
  expect_true(grepl("**1 (2026-07-10) was trained on the final-count slice",
                    .cascade_origin_training_txt(os, K = 6L), fixed = TRUE))
})

test_that("the rolling-scale sentence branches on the cross-validated model's ROLLING_PREDICTORS [R4e]", {
  .need("cascade_write_report")
  wb <- paste(deparse(body(cascade_write_report)), collapse = "\n")
  expect_true(grepl("rolling_predictors_requested", wb, fixed = TRUE))
  expect_true(grepl("get0(\"ROLLING_PREDICTORS\"", wb, fixed = TRUE))
  expect_true(grepl("so its intercept is on that model's scale", wb, fixed = TRUE))
  code <- paste(sub("#.*$", "", readLines(file.path(.ST, "run_cascade.R"), warn = FALSE)),
                collapse = "\n")
  expect_true(grepl("rolling_predictors_requested = isTRUE(get0(\"ROLLING_PREDICTORS\"", code,
                    fixed = TRUE))
})

test_that("the in-sample pooled count is stored as count_residual_in_sample, not verification [R4d]", {
  code <- paste(sub("#.*$", "", readLines(file.path(.ST, "run_cascade.R"), warn = FALSE)),
                collapse = "\n")
  expect_true(grepl("count_residual_in_sample = if (is.null(calib$report))", code, fixed = TRUE))
  expect_false(grepl("verification = ", code, fixed = TRUE))
})

# ---- [R2] downstream readers use only the canonical layer cache -------------------------------
test_that("41 and the Figure-2 cascade script read only the canonical, run-stamped layer [R2]", {
  .code <- function(f) paste(sub("#.*$", "", readLines(file.path(.ST, f), warn = FALSE)),
                             collapse = "\n")
  for (f in c("41_kindiv_sweep.R", "make_manuscript_figure2_cascade.R")) {
    cf <- .code(f)
    expect_false(grepl("CASCADE_LAYER_CACHE", cf, fixed = TRUE), info = f)
    expect_true(grepl("\"layer_cache_v2.rds\"", cf, fixed = TRUE), info = f)
  }
  # run_cascade.R stamps the run id BEFORE it writes the canonical copy
  rc <- strsplit(.code("run_cascade.R"), "\n", fixed = TRUE)[[1]]
  i_id  <- grep("layer$cascade_run_id <- CASCADE_RUN_ID", rc, fixed = TRUE)
  i_can <- grep("saveRDS(layer, .canon_cache)", rc, fixed = TRUE)
  expect_length(i_id, 1L); expect_length(i_can, 1L)
  expect_lt(i_id, i_can)
})
