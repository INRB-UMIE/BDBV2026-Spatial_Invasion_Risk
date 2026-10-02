# =============================================================================
# tests/test_tv_beta.R — time-varying import coefficient beta_t (21_bayesian_renewal.R)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# These models are the answer to "has the import->invasion conversion rate changed since the
# epidemic left Ituri?", so the thing that has to be right is not only that they FIT but that
# they FORECAST the way each process says they should. The tests are therefore split:
#
#   A. .tv_rw1_cov()        the RW1 covariance is exactly Brownian motion, is invariant to the
#                           week-index origin, and handles gaps.
#   B. .tv_offsets()        the forecast week is max training week + horizon, and the rw1
#                           grouping level is pinned to a TRAINING level (never a new one).
#   C. dispatch             .tv_is_structural() sends trend to the covariate fitter and every
#                           other process to the group-level/GP fitter.
#   D. forecast behaviour   on a simulated series with a KNOWN drifting beta: rw1 PERSISTS the
#                           last level, week REVERTS to beta0, and rw1's predictive SD grows
#                           with the horizon. This is the property that distinguishes the
#                           models, and it is the one a silent brms change would break.
#   E. recovery             the rw1 fit recovers the direction and rough size of a simulated
#                           drift in beta.
#
# D and E fit real Stan models and are SKIPPED unless cmdstanr is installed and
# RUN_SLOW_TESTS=1, matching the convention in the other model tests.
#
# Run:  Rscript spatiotemporal/tests/run_tests.R
#   or: RUN_SLOW_TESTS=1 Rscript -e 'testthat::test_file("spatiotemporal/tests/test_tv_beta.R")'
# =============================================================================

suppressPackageStartupMessages(library(testthat))

.st <- if (requireNamespace("here", quietly = TRUE)) file.path(here::here(), "spatiotemporal") else "."
if (!exists("OUT_DIR")) try(suppressMessages(source(file.path(.st, "00_config.R"))), silent = TRUE)
if (!exists(".tv_rw1_cov", mode = "function"))
  try(suppressWarnings(suppressMessages(source(file.path(.st, "21_bayesian_renewal.R")))),
      silent = TRUE)
.skip_unless <- function(fn) if (!exists(fn, mode = "function")) skip(paste(fn, "not available"))
.slow <- identical(trimws(Sys.getenv("RUN_SLOW_TESTS", "")), "1") &&
         requireNamespace("cmdstanr", quietly = TRUE) &&
         requireNamespace("brms", quietly = TRUE)

# ---------------------------------------------------------------------------------------
# A. the RW1 covariance
# ---------------------------------------------------------------------------------------
test_that(".tv_rw1_cov is the Brownian covariance, anchored at the first observed week", {
  .skip_unless(".tv_rw1_cov")
  K <- .tv_rw1_cov(c(2, 3, 4, 5))
  expect_equal(dim(K), c(4L, 4L))
  expect_equal(dimnames(K)[[1]], c("2", "3", "4", "5"))
  # tau = 1..4, K[i,j] = min(tau_i, tau_j)
  expect_equal(unname(K), outer(1:4, 1:4, pmin) * 1.0)
  # Under brms' parameterisation Cov(u) = sd^2 K, so Var(u_t - u_s) = sd^2 |tau_t - tau_s|:
  # a unit-variance random walk. Check the increment variance for every pair.
  for (i in 1:4) for (j in 1:4)
    expect_equal(K[i, i] + K[j, j] - 2 * K[i, j], abs(i - j) * 1.0)
  expect_true(all(eigen(K, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that(".tv_rw1_cov does not depend on where the week index happens to start", {
  .skip_unless(".tv_rw1_cov")
  # Folds begin at different absolute week indices; the prior must not change because of it.
  expect_equal(unname(.tv_rw1_cov(2:7)), unname(.tv_rw1_cov(20:25)))
})

test_that(".tv_rw1_cov charges a gap the right amount of variance", {
  .skip_unless(".tv_rw1_cov")
  K <- .tv_rw1_cov(c(1, 2, 4))          # a missing week 3
  # Var(u_4 - u_2) must be 2 (two innovations), not 1.
  expect_equal(K[3, 3] + K[2, 2] - 2 * K[3, 2], 2)
})

# ---------------------------------------------------------------------------------------
# B. forecast-week construction
# ---------------------------------------------------------------------------------------
test_that(".tv_offsets puts the forecast at max training week + horizon", {
  .skip_unless(".tv_offsets")
  fit <- structure(list(), class = "list")
  attr(fit, "max_train_week") <- 12L
  attr(fit, "train_weeks")    <- 2:12
  attr(fit, "tv")             <- "gp"
  off <- data.frame(health_zone = c("A", "A"), horizon = c(1L, 2L), logLam = c(-1, -1))
  attr(off, "projection") <- list(marker = TRUE)
  od <- .tv_offsets(fit, off, list(d = data.frame(.week = 2:12)))
  expect_equal(od$.week, c(13L, 14L))
  expect_identical(attr(od, "projection"), list(marker = TRUE))   # projection must survive
})

test_that(".tv_offsets pins the rw1 grouping level to a TRAINING level, not a new one", {
  .skip_unless(".tv_offsets")
  fit <- structure(list(), class = "list")
  attr(fit, "max_train_week") <- 12L
  attr(fit, "train_weeks")    <- 2:12
  attr(fit, "tv")             <- "rw1"
  off <- data.frame(health_zone = "A", horizon = 1L, logLam = -1)
  od <- .tv_offsets(fit, off, list(d = data.frame(.week = 2:12)))
  # The rw1 path predicts with re_formula = NA and adds the propagated walk itself, so a NEW
  # level here would either fail brms' validation or invite it to draw u from N(0, sd) —
  # which is precisely the wrong forecast for a random walk.
  expect_true(as.character(od$.week_f) %in% as.character(2:12))
  expect_equal(as.integer(as.character(od$.week_f)), 12L)
  expect_equal(od$.week, 13L)   # the GP/index column still moves to the forecast week
})

test_that(".tv_offsets gives the iid-week model a genuinely new level", {
  .skip_unless(".tv_offsets")
  fit <- structure(list(), class = "list")
  attr(fit, "max_train_week") <- 12L
  attr(fit, "train_weeks")    <- 2:12
  attr(fit, "tv")             <- "week"
  off <- data.frame(health_zone = "A", horizon = c(1L, 2L), logLam = -1)
  od <- .tv_offsets(fit, off, list(d = data.frame(.week = 2:12)))
  expect_equal(as.integer(as.character(od$.week_f)), c(13L, 14L))
})

# ---------------------------------------------------------------------------------------
# C. dispatch
# ---------------------------------------------------------------------------------------
test_that(".tv_is_structural routes trend to the covariate fitter and the rest to the tv fitter", {
  .skip_unless(".tv_is_structural")
  expect_false(.tv_is_structural("none"))
  expect_false(.tv_is_structural(NULL))
  # trend IS the week_idx covariate; sending it to the group-level fitter would fit the trend
  # twice over (once as a covariate, once as a process).
  expect_false(.tv_is_structural("trend"))
  for (ty in c("week", "rw1", "ar1", "gp")) expect_true(.tv_is_structural(ty))
})

test_that("mobility_kernel_from_method strips every variant suffix the grid can produce", {
  .skip_unless("mobility_kernel_from_method")
  for (sfx in c("", "-med", "-geo", "-tvtrend", "-tvweek", "-tvrw1", "-tvar1", "-tvgp",
                "-gtshort", "-gtlong", "-geo-tvrw1", "-geo-tvgp"))
    expect_equal(mobility_kernel_from_method(paste0("Bayes-M14-fill", sfx)), "M14-fill",
                 info = sfx)
})

# ---------------------------------------------------------------------------------------
# D/E. fitted behaviour (slow: real Stan fits)
# ---------------------------------------------------------------------------------------
# One simulated system shared by the slow tests: 16 weeks x 250 at-risk zones, with the import
# coefficient DRIFTING DOWNWARD as a random walk — a fall of about 1.9 log units (beta roughly
# sevenfold lower by week 16), the qualitative pattern an outbreak being brought under control
# would show, and the one these models exist to detect.
#
# WHY THE SIMULATION IS DELIBERATELY STRONGER THAN THE REAL DATA. These tests check the
# IMPLEMENTATION (does the walk propagate, is the level persisted, is the shape recovered),
# not the statistical power of the real study. At the real signal level the walk is barely
# identified: on a 60-zone version of this simulation the ORACLE per-week maximum-likelihood
# estimator — which knows the true offsets and has no shrinkage — correlates only ~0.26 with
# the truth, so no estimator could pass a recovery test there and passing one would mean the
# test was measuring noise. At 250 zones the oracle correlation is 0.82-0.94, so a failure
# here is an implementation failure.
.sim_design <- function(seed = 42L, T = 16L, nz = 250L) {
  set.seed(seed)
  u  <- cumsum(stats::rnorm(T, mean = -0.15, sd = 0.06))   # drifting walk on log beta
  d  <- expand.grid(.z = seq_len(nz), .week = seq_len(T))
  d$logLam <- stats::rnorm(nrow(d), mean = 0.2, sd = 0.9)
  eta <- -2.5 + u[d$.week] + d$logLam
  d$invaded <- stats::rbinom(nrow(d), 1L, 1 - exp(-exp(eta)))
  d$.zone   <- paste0("Z", d$.z); d$.z <- NULL
  list(d = d, feat = character(0), n_events = sum(d$invaded), n_obs = nrow(d),
       beta0 = exp(-2.5), center = numeric(0), scale = numeric(0),
       mob = "sim", gt = "medium", .u_true = u)
}

test_that("every time-varying variant fits and carries its process attribute", {
  if (!.slow) skip("slow: set RUN_SLOW_TESTS=1 (needs cmdstanr)")
  .skip_unless("fit_bayes_renewal_tv")
  des <- .sim_design()
  for (ty in c("week", "rw1", "ar1", "gp")) {
    f <- fit_bayes_renewal_tv(des, iter = 600L, chains = 2L, tv = ty)
    expect_false(is.null(f), info = ty)
    expect_equal(attr(f, "tv"), ty)
    expect_equal(attr(f, "max_train_week"), 16L, info = ty)
    expect_equal(attr(f, "train_weeks"), 1:16, info = ty)
  }
})

test_that("the rw1 forecast PERSISTS the last level while the iid-week forecast REVERTS", {
  if (!.slow) skip("slow: set RUN_SLOW_TESTS=1 (needs cmdstanr)")
  .skip_unless("fit_bayes_renewal_tv")
  des <- .sim_design()
  # One forecast row per horizon, offset held at 0 so the linear predictor IS beta0 + u.
  off <- data.frame(health_zone = "Z1", horizon = 1:3, logLam = 0)

  eta_of <- function(ty) {
    f  <- fit_bayes_renewal_tv(des, iter = 800L, chains = 2L, tv = ty)
    expect_false(is.null(f), info = ty)
    od <- .tv_offsets(f, off, des)
    lp <- .tv_linpred_fn(f)
    e  <- if (is.function(lp)) lp(f, od) else
      brms::posterior_linpred(f, newdata = od, allow_new_levels = TRUE,
                              sample_new_levels = "gaussian")
    list(fit = f, eta = e)
  }
  r_rw <- eta_of("rw1"); r_wk <- eta_of("week")

  # The simulated walk ENDS about 1.9 log-units below where it started, so a model that
  # persists the level must forecast BELOW the fitted intercept, and one that reverts must
  # forecast AT the intercept. Compare each forecast mean with that model's OWN posterior
  # intercept — the intercept and the walk are only jointly identified, so comparing across
  # models, or with the simulation's -2.5, would not be a like-for-like test.
  b0_rw <- mean(posterior::as_draws_df(r_rw$fit)$b_Intercept)
  b0_wk <- mean(posterior::as_draws_df(r_wk$fit)$b_Intercept)
  m_rw  <- colMeans(r_rw$eta); m_wk <- colMeans(r_wk$eta)

  # week: reverts to the intercept (E[u] = 0 for a new level).
  expect_lt(abs(mean(m_wk) - b0_wk), 0.25)
  # rw1: carries the walk's last level, which is materially below the intercept here.
  expect_lt(mean(m_rw), b0_rw - 0.2)
  # ... and the persisted level is nearly constant across horizons (a walk has no drift),
  # whereas the *uncertainty* about it grows.
  expect_lt(max(abs(diff(m_rw))), 0.25)
  sd_rw <- apply(r_rw$eta, 2, stats::sd)
  expect_true(all(diff(sd_rw) > 0),
              info = paste("rw1 predictive SD by horizon:", paste(round(sd_rw, 3), collapse = ", ")))
  # The growth is the random walk's: Var(u_{T+k}) - Var(u_T) = k * sigma^2, so the SD at
  # horizon 3 must exceed that at horizon 1. (A weak form of the exact law, because the
  # intercept's own posterior variance is included in both.)
  expect_gt(sd_rw[3], sd_rw[1])
})

test_that("the rw1 innovations are shared within a draw and across zones", {
  if (!.slow) skip("slow: set RUN_SLOW_TESTS=1 (needs cmdstanr)")
  .skip_unless("fit_bayes_renewal_tv")
  des <- .sim_design()
  f <- fit_bayes_renewal_tv(des, iter = 600L, chains = 2L, tv = "rw1")
  expect_false(is.null(f))
  # Two zones at the same horizon: u_t is a TIME effect, so within a draw the two rows must
  # differ ONLY by their offsets — here identical — i.e. be exactly equal. Drawing the
  # innovation per row would turn a common shock into independent noise.
  off <- data.frame(health_zone = c("Z1", "Z2"), horizon = c(2L, 2L), logLam = c(0, 0))
  e <- .tv_linpred_fn(f)(f, .tv_offsets(f, off, des))
  expect_equal(e[, 1], e[, 2], tolerance = 1e-10)
})

test_that("the rw1 fit recovers the direction of a simulated drift in beta", {
  if (!.slow) skip("slow: set RUN_SLOW_TESTS=1 (needs cmdstanr)")
  .skip_unless("fit_bayes_renewal_tv")
  des <- .sim_design()
  f <- fit_bayes_renewal_tv(des, iter = 1000L, chains = 2L, tv = "rw1")
  expect_false(is.null(f))
  dr <- posterior::as_draws_df(f)
  u  <- vapply(1:16, function(t) mean(dr[[sprintf("r_.week_f[%d,Intercept]", t)]]), numeric(1))
  # The truth falls by ~1.9 log units from week 1 to week 16. Require the posterior mean walk
  # to fall substantially too, and to correlate strongly with the simulated one. (The absolute
  # level is not identified separately from the intercept, so only the SHAPE is testable —
  # hence a correlation and a difference, never an absolute value.)
  expect_lt(u[16] - u[1], -0.8)
  expect_gt(stats::cor(u, des$.u_true), 0.7)
  # sigma must be bounded away from 0: the data really do carry week-to-week variation.
  expect_gt(stats::quantile(dr[["sd_.week_f__Intercept"]], 0.05, names = FALSE), 0.02)
})

test_that("a fixed-beta series does NOT produce a spurious drift", {
  if (!.slow) skip("slow: set RUN_SLOW_TESTS=1 (needs cmdstanr)")
  .skip_unless("fit_bayes_renewal_tv")
  # The converse of the recovery test, and the more important one for a publication: fitted to
  # data generated with a CONSTANT beta, the walk must stay flat rather than inventing a trend.
  set.seed(7L)
  T <- 16L; nz <- 250L                       # same size as .sim_design(), constant beta
  d <- expand.grid(.z = seq_len(nz), .week = seq_len(T))
  d$logLam  <- stats::rnorm(nrow(d), 0.2, 0.9)
  d$invaded <- stats::rbinom(nrow(d), 1L, 1 - exp(-exp(-2.5 + d$logLam)))
  d$.zone <- paste0("Z", d$.z); d$.z <- NULL
  des <- list(d = d, feat = character(0), n_events = sum(d$invaded), n_obs = nrow(d),
              beta0 = exp(-2.5), center = numeric(0), scale = numeric(0),
              mob = "sim", gt = "medium")
  f <- fit_bayes_renewal_tv(des, iter = 1000L, chains = 2L, tv = "rw1")
  expect_false(is.null(f))
  dr <- posterior::as_draws_df(f)
  u  <- vapply(1:T, function(t) mean(dr[[sprintf("r_.week_f[%d,Intercept]", t)]]), numeric(1))
  expect_lt(max(u) - min(u), 0.6)
  # sigma -> 0 is the honest answer here; require the posterior to allow it.
  expect_lt(stats::quantile(dr[["sd_.week_f__Intercept"]], 0.05, names = FALSE), 0.2)
})

test_that("bayes_beta_trajectory reproduces each process's forecast behaviour", {
  if (!.slow) skip("slow: set RUN_SLOW_TESTS=1 (needs cmdstanr)")
  .skip_unless("bayes_beta_trajectory")
  des <- .sim_design()
  T <- max(des$d$.week)
  wd <- stats::setNames(as.character(as.Date("2026-04-28") + 7L * (seq_len(T + 2L) - 1L)),
                        as.character(seq_len(T + 2L)))
  get_tj <- function(ty) {
    f <- if (identical(ty, "none")) fit_bayes_renewal(des, iter = 700L, chains = 2L)
         else fit_bayes_renewal_tv(des, iter = 700L, chains = 2L, tv = ty)
    expect_false(is.null(f), info = ty)
    tj <- bayes_beta_trajectory(f, des, paste0("sim-", ty), horizons = 1:2, week_dates = wd)
    expect_false(is.null(tj), info = ty)
    expect_equal(nrow(tj), T + 2L, info = ty)
    expect_equal(sum(tj$is_forecast), 2L, info = ty)
    expect_equal(tj$week_date[1], as.Date("2026-04-28"), info = ty)
    expect_true(all(tj$beta_lo <= tj$beta & tj$beta <= tj$beta_hi), info = ty)
    tj
  }
  # Fixed beta: a flat line, and no extra forecast uncertainty. It is the reference the
  # time-varying variants are read against, so it must be flat to numerical precision.
  t0 <- get_tj("none")
  expect_lt(max(t0$beta) / min(t0$beta) - 1, 1e-8)

  .ciw <- function(tj, i) tj$beta_hi[i] / tj$beta_lo[i]
  # The simulation's beta FALLS by ~1.9 log units, so every time-varying process must end the
  # training period well below where it started.
  for (ty in c("week", "rw1", "ar1", "gp")) {
    tj <- get_tj(ty)
    expect_lt(tj$beta[T], tj$beta[1], label = paste(ty, "beta[T] < beta[1]"))
    # Forecast uncertainty must GROW: none of these processes knows next week's level.
    expect_gt(.ciw(tj, T + 2L), .ciw(tj, T), label = paste(ty, "forecast CrI widens"))
    assign(paste0(".tj_", ty), tj, inherits = FALSE)
  }
  # THE DISCRIMINATING TEST. All four fit the same falling series; they differ only in what
  # they carry into the forecast, and that is the whole reason for fitting more than one:
  #   rw1  persists the final level  -> forecast stays at beta[T]
  #   week reverts to beta0          -> forecast jumps back UP, well above beta[T]
  #   ar1  decays toward beta0       -> strictly between the two
  # Expressed as the forecast's displacement from the last fitted level, on the log scale the
  # model works on.
  lift <- function(tj) log(tj$beta[T + 2L]) - log(tj$beta[T])
  expect_lt(abs(lift(.tj_rw1)), 0.35)              # rw1: persists
  expect_gt(lift(.tj_week), 0.35)                  # week: reverts upward
  expect_gt(lift(.tj_ar1), lift(.tj_rw1))          # ar1: decays, so above rw1 ...
  expect_lt(lift(.tj_ar1), lift(.tj_week))         # ... but not all the way back
})

test_that("the ar1 lengthscale maps onto an AR coefficient in (0, 1)", {
  if (!.slow) skip("slow: set RUN_SLOW_TESTS=1 (needs cmdstanr)")
  .skip_unless("fit_bayes_renewal_tv")
  des <- .sim_design()
  f <- fit_bayes_renewal_tv(des, iter = 800L, chains = 2L, tv = "ar1")
  expect_false(is.null(f))
  tb <- bayes_param_table(f, "sim-ar1")
  row <- tb[tb$term == "ar1_phi_week", ]
  expect_equal(nrow(row), 1L)
  expect_true(row$hr > 0 && row$hr < 1)
  expect_true(row$lo <= row$hr && row$hr <= row$hi)
  # It must NOT be labelled a hazard ratio: it is a time constant, not an effect size.
  expect_false(grepl("hazard ratio", row$effect_scale))
})
