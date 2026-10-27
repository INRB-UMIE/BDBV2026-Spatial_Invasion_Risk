# =============================================================================
# tests/test_epi_params.R — Unit tests for 02_epi_params.R functions
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Tests cover:
#   - make_gt_pmf()               generation time PMF correctness
#   - compute_truncation_weights() Exponential-CDF weight computation
#   - the right-truncation: support (trunc_max), as-of series basis, cached-fit gates,
#     and the R(t) cache key
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: skip gracefully when a function is unavailable
# ---------------------------------------------------------------------------
skip_if_missing <- function(fn_name) {
  if (!existsFunction(fn_name)) {
    skip(paste("Function", fn_name, "not available — source 02_epi_params.R first"))
  }
}

existsFunction <- function(fn_name) {
  exists(fn_name, mode = "function", envir = .GlobalEnv, inherits = TRUE)
}

# =============================================================================
# 1. make_gt_pmf() — generation time PMF
# =============================================================================

test_that("GT PMF sums to exactly 1 for all three profiles", {
  skip_if_missing("make_gt_pmf")

  for (profile_name in c("short", "medium", "long")) {
    p   <- GT_PROFILES[[profile_name]]
    pmf <- make_gt_pmf(p$mean, p$sd, p$max_tau)

    expect_true(
      is.numeric(pmf),
      info = paste(profile_name, ": GT PMF must be numeric")
    )
    expect_equal(
      length(pmf), p$max_tau,
      info = paste(profile_name, ": GT PMF length must equal max_tau")
    )
    # make_gt_pmf() internally calls stopifnot(abs(sum - 1) < 1e-9) so this
    # is guaranteed; we verify it from the test side too.
    expect_equal(
      sum(pmf), 1,
      tolerance = 1e-9,
      info = paste(profile_name, ": GT PMF must sum to 1")
    )
    expect_true(
      all(pmf >= 0),
      info = paste(profile_name, ": GT PMF must be non-negative")
    )
  }
})

test_that("GT PMF first element is near-zero for slow Gamma (medium profile)", {
  skip_if_missing("make_gt_pmf")

  # medium: mean=9d, sd=4.5d → shape=4, rate=0.444
  # P(Gamma(4, 0.444) <= 1) = pgamma(1, 4, 0.444) ≈ 0.00013; should be tiny
  pmf <- make_gt_pmf(9, 4.5, 30)
  expect_true(
    pmf[1] < 0.01,
    info = "P(GT=1 day) for mean-9d Gamma should be < 1%"
  )
})

test_that("GT PMF first element is far smaller for a long than a short generation time", {
  skip_if_missing("make_gt_pmf")

  # THRESHOLD CORRECTED (was: pmf_long[1] < 1e-4, which had failed since make_gt_pmf()
  # switched to the DOUBLE-INTERVAL-CENSORED discretisation). 1e-4 was calibrated for the
  # naive single-interval form, where P(GT = 1) = F(1) - F(0) = 1.85e-05 for this Gamma.
  # Under primary censoring the recorded day-difference is floor(U + D) with U ~ Unif(0,1),
  # so tau = 1 collects the mass with D in (1 - U, 2 - U) — an order of magnitude more,
  # 2.02e-04. The code is right (primarycensored and the base-R integral fallback agree to
  # 5e-15); the constant was stale. Asserting a RATIO rather than an absolute floor also
  # makes the test scale-free, so it keeps testing the intended property — a long GT puts
  # far less mass on a 1-day interval than a short one — instead of pinning a magic number
  # that any future change to the discretisation convention would silently invalidate.
  pmf_long   <- make_gt_pmf(13, 5.5, 40)   # shape ~ 5.6
  pmf_medium <- make_gt_pmf(9,  4.5, 40)   # shape = 4

  expect_lt(pmf_long[1], pmf_medium[1] / 10)   # observed ratio ~0.038
  expect_lt(pmf_long[1], 1e-3)
  expect_gt(pmf_long[1], 0)
})

test_that("GT PMF mode is near expected value for medium profile", {
  skip_if_missing("make_gt_pmf")

  # medium: mean=9, sd=4.5 → shape=4, rate=1/4.5*9 = 0.444
  # Continuous Gamma mode = (shape-1)/rate = 3/0.444 ≈ 6.75 d
  # Discretised mode should be tau=7 (±2 allowed for discretisation rounding)
  pmf      <- make_gt_pmf(9, 4.5, 30)
  mode_tau <- which.max(pmf)
  expect_true(
    mode_tau >= 5 & mode_tau <= 10,
    info = sprintf(
      "Medium GT PMF mode should be 5-10 days; got tau=%d (pmf=%.4f)",
      mode_tau, pmf[mode_tau]
    )
  )
})

test_that("GT PMF mode shifts right as mean increases across profiles", {
  skip_if_missing("make_gt_pmf")

  mode_short  <- which.max(make_gt_pmf(5.5,  2.0, 20))
  mode_medium <- which.max(make_gt_pmf(9.0,  4.5, 30))
  mode_long   <- which.max(make_gt_pmf(13.0, 5.5, 40))

  expect_true(
    mode_short < mode_medium,
    info = "Short GT mode must be earlier than medium GT mode"
  )
  expect_true(
    mode_medium < mode_long,
    info = "Medium GT mode must be earlier than long GT mode"
  )
})

test_that("make_gt_pmf rejects invalid inputs", {
  skip_if_missing("make_gt_pmf")

  expect_error(make_gt_pmf(mean = -1, sd = 4, max_tau = 30),
               info = "Negative mean must error")
  expect_error(make_gt_pmf(mean = 9, sd = 0, max_tau = 30),
               info = "Zero sd must error")
  expect_error(make_gt_pmf(mean = 9, sd = 4, max_tau = 0),
               info = "Zero max_tau must error")
})

test_that("GT PMF is correctly normalised when truncation removes tail mass", {
  skip_if_missing("make_gt_pmf")

  # With max_tau = 5 on a mean-9d Gamma, a lot of mass is truncated.
  # The result must still sum exactly to 1 (renormalisation).
  pmf_truncated <- make_gt_pmf(9, 4.5, 5)
  expect_equal(sum(pmf_truncated), 1, tolerance = 1e-9,
               info = "Hard-truncated GT PMF must renormalise to 1")
  expect_equal(length(pmf_truncated), 5L,
               info = "Hard-truncated GT PMF must have length = max_tau")
})

# =============================================================================
# 3. compute_truncation_weights() — fitted-delay completeness weights
#
# Signature: compute_truncation_weights(week_start_dates, analysis_date, delay, rate)
#
# The weight is E_o[P(delay <= analysis_date - o)] over the seven onset days o of the
# bucket, i.e. the MEAN of the delay CDF across the week — NOT the CDF evaluated once at
# the week midpoint (which, F being concave, overstates completeness). P(Delta <= L) is
# F(L + 0.5) under the daily-rounding convention the delay is fit under. Passing `rate`
# pins a plain Exp(rate); the default resolves the fitted delay in force.
#
# Reference implementation used by the analytic tests below:
.wk_weight_ref <- function(week_start, analysis_date, rate, week_days = 7L) {
  lags <- as.numeric(analysis_date - week_start) - seq.int(0L, week_days - 1L)
  cdf  <- ifelse(lags < 0, 0, 1 - exp(-rate * (lags + 0.5)))
  mean(cdf)
}
# =============================================================================

test_that("Truncation weight equals the mean Exp CDF across the bucket's onset days", {
  skip_if_missing("compute_truncation_weights")

  analysis_date <- as.Date("2026-05-15")
  week_start_7  <- analysis_date - 11L   # onset-day lags 11..5
  week_start_14 <- analysis_date - 18L   # onset-day lags 18..12

  w <- compute_truncation_weights(
    week_start_dates = c(week_start_7, week_start_14),
    analysis_date    = analysis_date,
    rate             = 0.2286
  )

  expect_equal(as.numeric(w[1]), .wk_weight_ref(week_start_7,  analysis_date, 0.2286),
               tolerance = 1e-12,
               info = "week-averaged weight must equal mean_j (1 - exp(-rate*(lag_j+0.5)))")
  expect_equal(as.numeric(w[2]), .wk_weight_ref(week_start_14, analysis_date, 0.2286),
               tolerance = 1e-12)

  # Jensen: F is concave, so the week-AVERAGED weight must sit strictly below F evaluated
  # at the week's MEAN lag. (Note this is NOT a comparison against the retired week-midpoint
  # form: that used week_start + 3.5, i.e. a lag a full day SHORTER than the true mean lag
  # of the seven integer onset days, so the two errors partly offset and the retired value
  # can land on either side.)
  mean_lag <- as.numeric(analysis_date - week_start_7) - 3 + 0.5
  expect_lt(as.numeric(w[1]), 1 - exp(-0.2286 * mean_lag))
})

test_that("Truncation weight increases monotonically with lag (older weeks = more complete)", {
  skip_if_missing("compute_truncation_weights")

  analysis_date <- as.Date("2026-06-01")
  # Five successive weeks, oldest to most recent
  week_starts <- analysis_date - c(35L, 28L, 21L, 14L, 7L)

  w <- compute_truncation_weights(week_starts, analysis_date, rate = 0.2286)

  expect_true(
    all(diff(w) < 0),   # as lag decreases, weight decreases
    info = "Weights must be monotonically decreasing as weeks get more recent"
  )
})

test_that("Truncation weights are in (0, 1] for all past weeks", {
  skip_if_missing("compute_truncation_weights")

  analysis_date <- as.Date("2026-06-15")
  past_weeks <- seq(as.Date("2026-05-01"), analysis_date - 7L, by = "week")

  w <- compute_truncation_weights(past_weeks, analysis_date, rate = 0.2286)

  valid <- !is.na(w)
  expect_true(
    all(w[valid] > 0 & w[valid] <= 1),
    info = "All observed-week weights must be in (0, 1]"
  )
})

test_that("Future weeks receive NA truncation weight", {
  skip_if_missing("compute_truncation_weights")

  analysis_date  <- as.Date("2026-06-01")
  future_week    <- analysis_date + 7L   # clearly in the future (midpoint > analysis_date)
  past_week      <- analysis_date - 14L

  w <- suppressWarnings(
    compute_truncation_weights(c(past_week, future_week), analysis_date, rate = 0.2286)
  )

  expect_false(is.na(w[1]), info = "Past week must have a non-NA weight")
  expect_true(is.na(w[2]),  info = "Future week must receive NA weight")
})

test_that("Truncation weight formula exactly matches the week-averaged Exponential CDF", {
  skip_if_missing("compute_truncation_weights")

  rate          <- 0.2286
  analysis_date <- as.Date("2026-07-01")

  # week_start = analysis_date - 14 -> onset-day lags 14, 13, ..., 8
  week_start_a <- analysis_date - 14L
  w_a <- compute_truncation_weights(week_start_a, analysis_date, rate = rate)
  expect_equal(as.numeric(w_a), .wk_weight_ref(week_start_a, analysis_date, rate),
               tolerance = 1e-12,
               info = "Weight must equal mean_j (1 - exp(-rate*(lag_j + 0.5))), j = 0..6")

  # week_start = analysis_date - 21 -> onset-day lags 21, 20, ..., 15
  week_start_b <- analysis_date - 21L
  w_b <- compute_truncation_weights(week_start_b, analysis_date, rate = rate)
  expect_equal(as.numeric(w_b), .wk_weight_ref(week_start_b, analysis_date, rate),
               tolerance = 1e-12)
})

test_that("Truncation weights default to the fitted delay, not the fixed lab Exp rate", {
  skip_if_missing("compute_truncation_weights")
  skip_if_missing("effective_onset_sample_delay")

  analysis_date <- as.Date("2026-07-31")
  week_start    <- analysis_date - 6L        # the current (most truncated) week

  dly <- effective_onset_sample_delay()
  w_default  <- compute_truncation_weights(week_start, analysis_date)
  w_explicit <- compute_truncation_weights(week_start, analysis_date, delay = dly)

  expect_equal(as.numeric(w_default), as.numeric(w_explicit), tolerance = 1e-12,
               info = "the default must BE effective_onset_sample_delay(), not a lab constant")

  # REGRESSION GUARD for the bug this replaced: whenever the fitted delay is slower than
  # the fixed lab reference, the default weight must be strictly below the lab weight —
  # i.e. the fold correction must be LARGER, not silently pinned to the fast lab delay.
  if (isTRUE(dly$source == "data") && dly$mean > 1 / DELAY_ONSET_SAMPLE_RATE) {
    w_lab <- compute_truncation_weights(week_start, analysis_date,
                                        rate = DELAY_ONSET_SAMPLE_RATE)
    expect_lt(as.numeric(w_default), as.numeric(w_lab))
  }
})

test_that("compute_truncation_weights rejects an ambiguous delay specification", {
  skip_if_missing("compute_truncation_weights")
  skip_if_missing("effective_onset_sample_delay")
  expect_error(
    compute_truncation_weights(as.Date("2026-07-01") - 14L, as.Date("2026-07-01"),
                               delay = effective_onset_sample_delay(), rate = 0.2),
    "not both"
  )
})

test_that("Higher rate gives higher weight for same lag (faster reporting)", {
  skip_if_missing("compute_truncation_weights")

  analysis_date <- as.Date("2026-06-15")
  week_start    <- analysis_date - 14L   # lag ≈ 10.5 d

  w_slow <- compute_truncation_weights(week_start, analysis_date, rate = 0.10)
  w_fast <- compute_truncation_weights(week_start, analysis_date, rate = 0.50)

  expect_true(
    as.numeric(w_fast) > as.numeric(w_slow),
    info = "Higher reporting rate → higher truncation weight for same lag"
  )
})

# =============================================================================
# 3. The right-truncation's SUPPORT, the as-of fit's series basis, the cached-fit gates and the
#    R(t) cache key (02_epi_params.R section 4b and estimate_rt_epinow2())
# =============================================================================

.skip_unless_epinow2 <- function() {
  if (!requireNamespace("EpiNow2", quietly = TRUE)) skip("EpiNow2 not installed")
}
# An UNCERTAIN lognormal (the shape estimate_truncation() returns) and a FIXED gamma (the as-of fit).
.unc_lnorm <- function(ml = 2.10, sl = 0.80, max = 45)
  EpiNow2::LogNormal(meanlog = EpiNow2::Normal(ml, 0.02), sdlog = EpiNow2::Normal(sl, 0.02), max = max)

test_that(".trunc_support_max is the support EpiNow2 itself uses, on both discretisation paths", {
  skip_if_missing(".trunc_support_max"); .skip_unless_epinow2()
  csd <- get0("create_stan_delays", envir = asNamespace("EpiNow2"), inherits = FALSE)
  if (!is.function(csd)) skip("EpiNow2:::create_stan_delays unavailable")
  unc <- .unc_lnorm()
  fix <- EpiNow2::Gamma(shape = 1.3, rate = 0.21, max = 45)
  su <- csd(trunc = EpiNow2::trunc_opts(unc))
  sf <- csd(trunc = EpiNow2::trunc_opts(fix))
  # Uncertain parameters go to Stan parametrically, evaluated on 0..delay_max.
  expect_equal(su$delay_n_p, 1)
  expect_equal(.trunc_support_max(unc), as.numeric(su$delay_max))
  # Fixed parameters go as their discretised pmf, on 0..(length - 1): one day shorter.
  expect_equal(sf$delay_n_np, 1)
  expect_equal(.trunc_support_max(fix), as.numeric(sf$delay_np_pmf_length) - 1)
  expect_equal(.trunc_support_max(unc) - .trunc_support_max(fix), 1)
})

test_that(".trunc_as_delay_spec carries the support, the moments ON it, and the right label", {
  skip_if_missing(".trunc_as_delay_spec"); .skip_unless_epinow2()
  g  <- EpiNow2::Gamma(shape = 1.3, rate = 0.21, max = 46)
  sp <- .trunc_as_delay_spec(g, "asof")
  expect_equal(sp$trunc_max, .trunc_support_max(g))
  b  <- sp$trunc_max + 0.5; Fb <- stats::pgamma(b, 1.3, 0.21)
  m1 <- stats::integrate(function(x) x * stats::dgamma(x, 1.3, 0.21), 0, b)$value / Fb
  m2 <- stats::integrate(function(x) x^2 * stats::dgamma(x, 1.3, 0.21), 0, b)$value / Fb
  expect_equal(sp$mean, m1, tolerance = 1e-7)
  expect_equal(sp$sd, sqrt(m2 - m1^2), tolerance = 1e-7)
  expect_equal(sp$rate, 1 / sp$mean)
  expect_identical(sp$estimator, "lynden_bell_cohort_gamma_asof")

  sl <- .trunc_as_delay_spec(.unc_lnorm(), "extract")
  expect_identical(sl$family, "lnorm")
  expect_identical(sl$estimator, "estimate_truncation_extract")
  pr <- sl$params; b <- sl$trunc_max + 0.5
  Fb <- stats::plnorm(b, pr[["meanlog"]], pr[["sdlog"]])
  m1 <- stats::integrate(function(x) x * stats::dlnorm(x, pr[["meanlog"]], pr[["sdlog"]]),
                         0, b)$value / Fb
  expect_equal(sl$mean, m1, tolerance = 1e-6)
  expect_lt(sl$mean, exp(pr[["meanlog"]] + pr[["sdlog"]]^2 / 2))   # conditioned, so shorter
})

test_that("the nowcast's completeness is EpiNow2's: both renormalise on the same support", {
  skip_if_missing(".trunc_as_delay_spec"); skip_if_missing("compute_truncation_weights")
  .skip_unless_epinow2()
  if (!requireNamespace("primarycensored", quietly = TRUE)) skip("primarycensored not installed")
  ad <- as.Date("2026-09-14")
  newest <- function(sp) as.numeric(suppressMessages(compute_truncation_weights(ad - 6, ad, delay = sp)))
  # Uncertain lognormal: Stan's pmf is primary-censored on 0..delay_max, normalised at delay_max + 1
  # (EpiNow2's discretised_pmf(), stan/functions/pmfs.stan).
  sp <- .trunc_as_delay_spec(.unc_lnorm(), "extract"); n <- sp$trunc_max
  pm <- primarycensored::dprimarycensored(0:n, stats::plnorm, pwindow = 1, swindow = 1, D = n + 1,
                                          meanlog = sp$params[["meanlog"]], sdlog = sp$params[["sdlog"]])
  ref <- mean(cumsum(pm)[1:7])
  unnorm <- mean(delay_cdf(sp, 0:6 + 0.5))            # the parametric CDF without its support
  expect_lt(abs(newest(sp) / ref - 1), 0.005)          # left: the F(L + 0.5) midpoint convention
  expect_lt(abs(newest(sp) / ref - 1), abs(unnorm / ref - 1))
  # Fixed gamma: EpiNow2 uses the discretised pmf itself.
  g  <- EpiNow2::Gamma(shape = 1.3, rate = 0.21, max = 46)
  pg <- EpiNow2::get_pmf(EpiNow2::discretise(g))
  expect_lt(abs(newest(.trunc_as_delay_spec(g, "asof")) / mean(cumsum(pg)[1:7]) - 1), 0.002)
})

test_that("the as-of fit measures from date_index of USABLE recorded onsets, not the raw onset", {
  skip_if_missing(".trunc_fit_asof_cohort"); .skip_unless_epinow2()
  set.seed(7L)
  end <- as.Date("2026-09-07"); start <- end - 120L
  n <- 4000L
  ons <- start + sample(0:120, n, replace = TRUE, prob = exp(0.02 * (0:120)))
  obs <- ons + round(stats::rgamma(n, shape = 1.2, rate = 0.15))
  vis <- obs <= end
  ll <- data.frame(confirmed = TRUE,
                   date_of_symptom_onset = as.character(ons[vis]),
                   date_of_sample_collection = as.character(obs[vis]),
                   date_index = ons[vis], onset_usable = TRUE, stringsAsFactors = FALSE)
  # .trunc_fit_asof_cohort() reads the observation date through linelist_observation_date();
  # stub it to the sample column.
  old_fn <- if (exists("linelist_observation_date", envir = globalenv()))
    get("linelist_observation_date", envir = globalenv()) else NULL
  on.exit({ if (!is.null(old_fn)) assign("linelist_observation_date", old_fn, envir = globalenv())
            else rm("linelist_observation_date", envir = globalenv()) }, add = TRUE)
  assign("linelist_observation_date",
         function(ll, end, caller = "test") as.Date(ll$date_of_sample_collection),
         envir = globalenv())
  fit <- function(x) .trunc_as_delay_spec(
    suppressMessages(suppressWarnings(.trunc_fit_asof_cohort(x, end, start))), "asof")
  base <- fit(ll)
  expect_false(is.null(base))
  # EpiNow2 is handed the support the curve was estimated on, 0..TRUNC_MAX_DAYS.
  expect_equal(base$trunc_max, as.numeric(get0("TRUNC_MAX_DAYS", ifnotfound = 45L)))
  # The raw onset column is NOT read: moving it changes nothing.
  ll2 <- ll; ll2$date_of_symptom_onset <- as.character(as.Date(ll$date_of_symptom_onset) + 2L)
  expect_identical(fit(ll2)$params, base$params)
  # Records without a usable recorded onset (an imputed date_index) are left out.
  junk <- ll[seq_len(300L), ]
  junk$onset_usable <- FALSE
  junk$date_index <- as.Date(junk$date_of_sample_collection) - 40L
  expect_identical(fit(rbind(ll, junk))$params, base$params)
  # Without the loader's columns the fit refuses rather than reading the raw onset.
  expect_null(suppressMessages(.trunc_fit_asof_cohort(ll[, c("confirmed", "date_of_symptom_onset",
                                                             "date_of_sample_collection")], end, start)))
})

test_that("a CACHED extract fit is re-gated by .trunc_sanity(), not served unchecked", {
  skip_if_missing(".trunc_fit"); skip_if_missing(".trunc_fingerprint"); .skip_unless_epinow2()
  tmp <- tempfile("trunc_cache_"); dir.create(tmp)
  had_od  <- exists("OUT_DIAGNOSTICS", envir = globalenv(), inherits = FALSE)
  old_od  <- if (had_od) get("OUT_DIAGNOSTICS", envir = globalenv()) else NULL
  had_rng <- exists("TRUNC_MEAN_RANGE", envir = globalenv(), inherits = FALSE)
  old_rng <- if (had_rng) get("TRUNC_MEAN_RANGE", envir = globalenv()) else NULL
  on.exit({
    if (had_od) assign("OUT_DIAGNOSTICS", old_od, envir = globalenv())
    else rm("OUT_DIAGNOSTICS", envir = globalenv())
    if (had_rng) assign("TRUNC_MEAN_RANGE", old_rng, envir = globalenv())
    else if (exists("TRUNC_MEAN_RANGE", envir = globalenv(), inherits = FALSE))
      rm("TRUNC_MEAN_RANGE", envir = globalenv())
    unlink(tmp, recursive = TRUE)
  }, add = TRUE)
  assign("OUT_DIAGNOSTICS", tmp, envir = globalenv())
  panel <- list(
    data.frame(date = seq(as.Date("2026-07-01"), as.Date("2026-08-01"), by = "day"), confirm = 3L),
    data.frame(date = seq(as.Date("2026-07-01"), as.Date("2026-08-08"), by = "day"), confirm = 3L))
  # The cache path exactly as .trunc_fit() derives it. Checked before .trunc_fit() is called, so
  # a mismatch fails here instead of letting a cache miss start a Stan fit.
  key  <- substr(rlang::hash(.trunc_fingerprint(panel, "extract")), 1L, 16L)
  path <- file.path(tmp, sprintf("epinow2_truncation_extract_%s.rds", key))
  d <- .unc_lnorm(2.0, 0.7)
  saveRDS(d, path)
  stopifnot(file.exists(path))
  assign("TRUNC_MEAN_RANGE", c(4, 25), envir = globalenv())
  expect_identical(suppressMessages(.trunc_fit(panel, "extract")), d)
  assign("TRUNC_MEAN_RANGE", c(20, 25), envir = globalenv())      # gate tightened after caching
  expect_warning(r <- suppressMessages(.trunc_fit(panel, "extract")), "REFUSED")
  expect_null(r)
})

test_that("the R(t) cache key carries the Stan seed, and ONE GT support for key and fit", {
  skip_if_missing("estimate_rt_epinow2")
  if (!isTRUE(get0(".HAVE_EPINOW2", ifnotfound = FALSE))) skip("EpiNow2 not installed")
  # No fit runs: the truncation is stubbed, the fit refuses, and file.exists() (which the cache
  # lookup calls once) records the path it is asked about. Every stub is restored on exit.
  nms <- c("epinow2_truncation", ".epinow2_rt_fit", "RANDOM_SEED", "file.exists")
  had <- vapply(nms, function(nm) exists(nm, envir = globalenv(), inherits = FALSE), logical(1))
  old <- lapply(nms[had], function(nm) get(nm, envir = globalenv())); names(old) <- nms[had]
  on.exit(for (nm in nms) {
    if (had[[nm]]) assign(nm, old[[nm]], envir = globalenv())
    else if (exists(nm, envir = globalenv(), inherits = FALSE)) rm(list = nm, envir = globalenv())
  }, add = TRUE)
  seen <- character(0)
  assign("epinow2_truncation", function(...) EpiNow2::Gamma(shape = 1.3, rate = 0.21, max = 46),
         envir = globalenv())
  assign(".epinow2_rt_fit", function(...) stop("no fit in tests"), envir = globalenv())
  assign("file.exists", function(...) { seen <<- c(seen, ..1); FALSE }, envir = globalenv())
  cases <- data.frame(date = seq(as.Date("2026-08-01"), by = "day", length.out = 30L), confirm = 5L)
  p <- GT_PROFILES$medium
  run <- function(seed, gtp) {
    assign("RANDOM_SEED", seed, envir = globalenv())
    suppressWarnings(suppressMessages(estimate_rt_epinow2(cases, "medium", gtp, as.Date("2026-09-01"))))
  }
  expect_null(run(1L, p))
  expect_null(run(2L, p))
  expect_length(seen, 2L)
  expect_true(grepl("_s1_", seen[1], fixed = TRUE))
  expect_true(grepl("_s2_", seen[2], fixed = TRUE))
  expect_false(identical(seen[1], seen[2]))
  # Omitting max_tau keys the profile's support, which is the support the fit uses.
  run(1L, p[c("mean", "sd")])
  expect_identical(seen[3], seen[1])
})
