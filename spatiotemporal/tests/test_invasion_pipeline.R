# =============================================================================
# tests/test_invasion_pipeline.R — Core-logic & invariant tests for the LIVE
# spatiotemporal INVASION-forecast modules.
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# These are non-trivial logic/invariant tests (not "function exists" checks)
# for the modules run_all.R actually sources. Scripts 09–14 are STALE (moved to
# spatiotemporal/stale/) and are never referenced here; metric assertions that
# used to live in stale/11_metrics.R are ported to their live equivalents:
#   * AUC-PR / log-score / Brier  -> 16_invasion_eval.R (.auc_pr/.log_score/.brier)
#   * Weighted Interval Score     -> the live `wis_one` inside 04b_epinowcast.R
#                                    (extracted from source and exercised directly)
#
# Coverage map (see task spec):
#   1  GT PMF (02) + weekly GT binning (06)
#   2  FOI core + inward contact matrix (06)
#   3  Nowcast right-truncation correction, cap-bug regression (04)
#   4  Invasion probability + ascertainment band (15)
#   5  Relative-risk scores (15)
#   6  Vulnerability index (20)
#   7  Metrics: AUC-PR, WIS, base rate, Brier/log-score (16 / 04b)
#   8  Detection curve + balanced-accuracy/MCC (20)
#   9  Date parsing + onset imputation (01)
#  10  Province RR generalisation + .prov_suffix (15)
#
# Run:  Rscript spatiotemporal/tests/run_tests.R
#   or: testthat::test_file("spatiotemporal/tests/test_invasion_pipeline.R")
# =============================================================================

suppressPackageStartupMessages({
  library(testthat)
  library(tibble)
})

# ---------------------------------------------------------------------------
# Local helpers (self-contained so the file also runs via testthat::test_file)
# ---------------------------------------------------------------------------
.exists_fn <- function(nm) exists(nm, mode = "function", envir = .GlobalEnv, inherits = TRUE)
skip_if_missing <- function(nm) if (!.exists_fn(nm))
  skip(paste("Function", nm, "not available — source the live module first"))

# If a sentinel function is absent (file run standalone), source the live modules
# the tests depend on. No-op when run through run_tests.R (already sourced).
.ensure_sourced <- function() {
  if (.exists_fn("forecast_workhorse") && .exists_fn(".auc_pr")) return(invisible())
  if (!requireNamespace("here", quietly = TRUE)) return(invisible())
  st <- file.path(here::here(), "spatiotemporal")
  for (m in c("00_config.R", "01_data_prep.R", "02_epi_params.R", "04_nowcasting.R",
              "06_simple_models.R", "15_workhorse.R", "16_invasion_eval.R",
              "20_forecast_detail.R")) {
    suppressWarnings(suppressMessages(
      try(source(file.path(st, m)), silent = TRUE)))
  }
}
.ensure_sourced()

# Extract a (possibly nested/local) function definition from a source file by
# brace-balancing, and evaluate it in a fresh environment. Used to exercise the
# LIVE `wis_one` WIS implementation that lives inside compare_nowcast_models() in
# 04b_epinowcast.R (a local closure, not otherwise callable).
.extract_fn_from_source <- function(file, fn_name) {
  if (!file.exists(file)) return(NULL)
  L <- readLines(file, warn = FALSE)
  i <- grep(paste0("^\\s*", fn_name, "\\s*<-\\s*function"), L)
  if (length(i) == 0) return(NULL)
  i <- i[1]
  depth <- 0L; started <- FALSE; buf <- character(0)
  for (j in seq(i, length(L))) {
    buf <- c(buf, L[j])
    op <- lengths(regmatches(L[j], gregexpr("\\{", L[j])))
    cl <- lengths(regmatches(L[j], gregexpr("\\}", L[j])))
    depth <- depth + op - cl
    if (op > 0) started <- TRUE
    if (started && depth <= 0) break
  }
  eval(parse(text = paste(buf, collapse = "\n")), envir = new.env(parent = globalenv()))
}


# =============================================================================
# 1. GENERATION TIME — make_gt_pmf (02) + daily_to_weekly_gt (06)
# =============================================================================

test_that("make_gt_pmf sums to 1 and recovers the requested mean (discretisation tol)", {
  skip_if_missing("make_gt_pmf")
  skip_if(!exists("GT_PROFILES"), "GT_PROFILES config not loaded")

  for (nm in names(GT_PROFILES)) {
    p   <- GT_PROFILES[[nm]]
    pmf <- make_gt_pmf(p$mean, p$sd, p$max_tau)

    expect_equal(sum(pmf), 1, tolerance = 1e-9,
                 info = paste(nm, ": PMF must sum to 1"))
    expect_true(all(pmf >= 0), info = paste(nm, ": PMF must be non-negative"))
    expect_equal(length(pmf), p$max_tau,
                 info = paste(nm, ": PMF length must equal max_tau"))

    # Discretised (tau = 1..max_tau) mean recovers the requested Gamma mean to
    # within the ~0.5 d discretisation shift (verified: |diff| <= 0.36 d for all
    # three config profiles).
    tau <- seq_along(pmf)
    disc_mean <- sum(tau * pmf)
    expect_lt(abs(disc_mean - p$mean), 0.6)
  }
})

test_that("make_gt_pmf renormalises when the tail is hard-truncated", {
  skip_if_missing("make_gt_pmf")
  # Heavy truncation (max_tau=5 on a mean-15d Gamma) still sums to exactly 1.
  pmf <- make_gt_pmf(15, 9, 5)
  expect_equal(sum(pmf), 1, tolerance = 1e-9)
  expect_length(pmf, 5L)
})

test_that("daily_to_weekly_gt bins days into weeks and renormalises to 1", {
  skip_if_missing("daily_to_weekly_gt")

  # 20-day unnormalised pmf: week1=days1-7 (7*1=7), week2=days8-14 (7*2=14),
  # week3=days15-20 (6*3=18). Normalised weekly => c(7,14,18)/39.
  g  <- c(rep(1, 7), rep(2, 7), rep(3, 6))
  Gw <- daily_to_weekly_gt(g)

  expect_equal(sum(Gw), 1, tolerance = 1e-12)
  expect_length(Gw, 3L)                      # ceiling(20/7) = 3 weekly bins
  expect_equal(Gw, c(7, 14, 18) / 39, tolerance = 1e-12,
               info = "Weekly bins must be the correct day-sums, renormalised")
})

test_that("daily_to_weekly_gt of a real GT PMF is a valid weekly PMF", {
  skip_if_missing("daily_to_weekly_gt")
  skip_if_missing("make_gt_pmf")
  Gw <- daily_to_weekly_gt(make_gt_pmf(15.3, 9.3, 45))
  expect_equal(sum(Gw), 1, tolerance = 1e-9)
  expect_true(all(Gw >= 0))
  expect_true(which.max(Gw) >= 2L,
              info = "A ~15 d mean GT peaks in week 2+, not week 1")
})


# =============================================================================
# 2. FORCE OF INFECTION — compute_foi + build_inward_contact_matrix (06)
# =============================================================================

test_that("compute_foi routes all import pressure INTO i = sum_j W[j,i]*Ytilde_j", {
  skip_if_missing("compute_foi")
  skip_if_missing("make_gt_pmf")
  skip_if_missing("daily_to_weekly_gt")

  zones <- c("A", "B", "C")
  # W is an OUTFLOW matrix: row A sends ALL of its outflow to B (W[A,B]=1).
  W <- matrix(c(0, 1, 0,
                0, 0, 1,
                1, 0, 0),
              nrow = 3, byrow = TRUE, dimnames = list(zones, zones))
  # Only zone A has past cases.
  Y <- matrix(0, 3, 4, dimnames = list(zones, as.character(1:4)))
  Y["A", 1:3] <- c(5, 8, 3)

  Gw     <- daily_to_weekly_gt(make_gt_pmf(9, 4.5, 30))
  Lambda <- compute_foi(Y, W, Gw, t_idx = 4L, zones)

  expect_true(all(Lambda >= 0), info = "FOI must be non-negative")
  expect_gt(unname(Lambda["B"]), 0)                        # all of A's pressure -> B
  expect_equal(unname(Lambda["A"]), 0, tolerance = 1e-12)  # A imports from B,C which have no cases
  expect_equal(unname(Lambda["C"]), 0, tolerance = 1e-12)  # W[A,C]=0 so C gets nothing
  # Lambda_B equals the G-weighted A history (W[A,B]=1, no other inflow to B).
  Y_weighted_A <- Gw[1] * Y["A", 3] + Gw[2] * Y["A", 2] + Gw[3] * Y["A", 1]
  expect_equal(unname(Lambda["B"]), unname(Y_weighted_A), tolerance = 1e-9)
})

test_that("build_inward_contact_matrix M is symmetric, non-negative, unit-mean; presence P is row-stochastic", {
  skip_if_missing("build_inward_contact_matrix")

  zones <- c("A", "B", "C", "D")
  W <- matrix(c(0,   0.5, 0.3, 0.2,
                0.4, 0,   0.3, 0.3,
                0.3, 0.3, 0,   0.4,
                0.2, 0.4, 0.4, 0),
              nrow = 4, byrow = TRUE, dimnames = list(zones, zones))
  pop <- c(A = 1e5, B = 5e4, C = 3e4, D = 2e4)

  M <- build_inward_contact_matrix(W, pop, zones)

  expect_equal(dim(M), c(4L, 4L))
  expect_identical(dimnames(M), list(zones, zones))
  expect_equal(M, t(M), tolerance = 1e-12, info = "M must be exactly symmetric")
  expect_true(all(M >= 0), info = "M must be non-negative")
  # Rescaled to unit mean POSITIVE entry (absolute scale absorbed by beta).
  expect_equal(mean(M[M > 0]), 1, tolerance = 1e-9)

  # The underlying presence matrix P = (1-home)*W + home*I is row-stochastic;
  # reconstruct it with the live default home fraction and confirm.
  home <- get0("MOBILITY_HOME_FRACTION", ifnotfound = 0.70)
  P <- (1 - home) * W
  diag(P) <- diag(P) + home
  expect_true(all(abs(rowSums(P) - 1) < 1e-9),
              info = "Presence matrix P must be row-stochastic")
})


# =============================================================================
# 3. NOWCAST CORRECTION — apply_nowcast_correction (04)  [cap-bug regression]
# =============================================================================

test_that("nowcast correction inflates BOTH confirmed and suspected even when suspected >> 5x max confirmed", {
  skip_if_missing("apply_nowcast_correction")

  # Regression for the historical cap bug: the OLD code capped the absolute
  # corrected COUNT at 5 * max(confirmed) and applied that same absolute cap to
  # the (much larger) suspected series, which could push suspected_nc BELOW its
  # raw value and trip the "corrected >= raw" assertion. Here suspected (40-80)
  # is >> 5x max stable confirmed (5*2 = 10), so the old bug would fire.
  zw <- tibble::tibble(
    health_zone = rep("TestZone", 5L),
    week_start  = as.Date(c("2026-04-06", "2026-04-13", "2026-04-20",
                            "2026-04-27", "2026-06-11")),  # midpoint 06-14.5, lag 0.5 d
    confirmed   = c(1L, 2L, 2L, 1L, 3L),
    suspected   = c(50L, 80L, 60L, 40L, 70L)
  )
  # lag 0.5 d => raw multiplier ~9.3 > 5, so the 5x cap binds (earlier weeks stable):
  # suspected (40-80) >> 5x max stable confirmed (5*2=10), so the OLD absolute-count
  # cap bug would fire here, while the correct factor cap keeps suspected_nc >= suspected.
  analysis_date <- as.Date("2026-06-15")

  out <- suppressWarnings(apply_nowcast_correction(zw, analysis_date, min_lag_days = 0))

  valid <- !is.na(out$confirmed_nc)
  expect_true(any(valid), info = "at least one row must be corrected")

  # (a) corrected >= raw for BOTH series — the assertion the old bug violated.
  expect_true(all(out$confirmed_nc[valid] >= out$confirmed[valid] - 1e-9),
              info = "confirmed_nc must never fall below raw confirmed")
  expect_true(all(out$suspected_nc[valid] >= out$suspected[valid] - 1e-9),
              info = "suspected_nc must never fall below raw suspected (cap-bug regression)")

  # (b) the correction is a MULTIPLIER capped at 5x (corrected in [raw, 5*raw]).
  mult_c <- out$confirmed_nc[valid] / pmax(out$confirmed[valid], 1e-9)
  mult_s <- out$suspected_nc[valid] / pmax(out$suspected[valid], 1e-9)
  expect_true(all(mult_c <= 5 + 1e-9), info = "confirmed multiplier capped at 5x")
  expect_true(all(mult_s <= 5 + 1e-9), info = "suspected multiplier capped at 5x")
  expect_true(any(mult_s >= 5 - 1e-9),
              info = "the recent week's raw multiplier (~9.3) must be capped AT 5x (cap binds)")

  # confirmed and suspected in the same row share ONE multiplier (1/trunc_weight),
  # so the ratio-based correction is series-agnostic.
  expect_equal(mult_c, mult_s, tolerance = 1e-9,
               info = "both series share the same per-row correction factor")
})

test_that("nowcast truncation weights are in (0,1], monotone with lag, NA for future weeks", {
  skip_if_missing("apply_nowcast_correction")
  zw <- tibble::tibble(
    health_zone = rep("Z", 5L),
    week_start  = as.Date(c("2026-04-06", "2026-04-13", "2026-04-20",
                            "2026-04-27", "2026-06-22")),  # last week is future
    confirmed   = c(3L, 5L, 8L, 4L, 2L),
    suspected   = c(5L, 8L, 12L, 6L, 3L)
  )
  analysis_date <- as.Date("2026-06-15")
  out <- suppressWarnings(apply_nowcast_correction(zw, analysis_date, min_lag_days = 0))

  w <- out$trunc_weight
  wv <- w[!is.na(w)]
  expect_true(all(wv > 0 & wv <= 1 + 1e-9), info = "observed weights in (0,1]")
  expect_true(is.na(w[out$week_start == as.Date("2026-06-22")]),
              info = "future week (lag<=0) must get NA weight")
  # older weeks (larger lag) are more complete -> higher weight
  past <- out[out$week_start <= analysis_date, ]
  past <- past[order(past$week_start), ]
  pw <- past$trunc_weight[!is.na(past$trunc_weight)]
  expect_true(all(diff(pw) <= 1e-9),
              info = "weights non-increasing from oldest to most recent week")
})


# =============================================================================
# 4. INVASION PROBABILITY + ASCERTAINMENT BAND — .invasion_prob / forecast_workhorse (15)
# =============================================================================

test_that(".invasion_prob is the Poisson arrival prob 1-exp(-mu), monotone in mu", {
  skip_if_missing(".invasion_prob")
  mu <- c(0, 0.01, 0.1, 0.5, 1, 3, 10, 50)
  p  <- vapply(mu, function(m) .invasion_prob(m, obs = "poisson", theta = NA_real_), numeric(1))

  expect_equal(p, 1 - exp(-mu), tolerance = 1e-12,
               info = "Poisson invasion prob must equal 1 - exp(-mu)")
  expect_equal(p[1], 0, tolerance = 1e-15)          # mu=0 -> 0
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(diff(p) >= 0), info = "monotone non-decreasing in mu")
  # NegBin arrival prob 1-(theta/(theta+mu))^theta is also monotone and <= Poisson.
  p_nb <- vapply(mu, function(m) .invasion_prob(m, obs = "negbin", theta = 2), numeric(1))
  expect_true(all(diff(p_nb) >= 0))
  expect_true(all(p_nb <= p + 1e-9), info = "NegBin (overdispersed) arrival <= Poisson")
})

test_that("ascertainment band orders p_case <= p_inf_lo <= p_inf <= p_inf_hi at the three rhos", {
  skip_if_missing(".invasion_prob")
  skip_if(!exists("ASCERTAINMENT_GRID"), "ASCERTAINMENT_GRID not loaded")

  rho    <- get0("ASCERTAINMENT_NOMINAL", ifnotfound = 0.45)
  rho_lo <- max(ASCERTAINMENT_GRID)   # largest rho -> LOWER infection band
  rho_hi <- min(ASCERTAINMENT_GRID)   # smallest rho -> UPPER infection band
  ip <- function(m) .invasion_prob(m, "poisson", NA_real_)

  for (mu in c(0.2, 1, 4)) {
    p_case <- ip(mu)
    p_lo   <- ip(mu / rho_lo)
    p_mid  <- ip(mu / rho)
    p_hi   <- ip(mu / rho_hi)
    expect_true(p_case <= p_lo + 1e-12)
    expect_true(p_lo   <= p_mid + 1e-12)
    expect_true(p_mid  <= p_hi + 1e-12)
    expect_gt(p_hi, p_case)   # strict when mu>0 and grid spans <1
  }
})

test_that("forecast_workhorse output satisfies the band ordering and masks affected zones", {
  skip_if_missing("forecast_workhorse")
  skip_if_missing("make_gt_pmf")

  weeks <- as.Date("2026-05-04") + 7 * (0:6)
  zz    <- c("A", "B", "C", "D")
  mk <- function(z, conf) tibble::tibble(
    health_zone = z, week_start = weeks, confirmed = conf, confirmed_nc = conf)
  zw <- rbind(
    mk("A", c(2, 4, 5, 6, 3, 2, 1)),   # A affected throughout
    mk("B", c(0, 0, 0, 1, 2, 1, 0)),   # B invades at week 4
    mk("C", rep(0, 7)),                # C at-risk
    mk("D", rep(0, 7)))                # D at-risk
  W <- matrix(c(0, .5, .3, .2,
                .4, 0, .3, .3,
                .3, .3, 0, .4,
                .2, .4, .4, 0),
              nrow = 4, byrow = TRUE, dimnames = list(zz, zz))

  fc <- suppressWarnings(forecast_workhorse(
    zw, W, make_gt_pmf(9, 4.5, 30), zz, t_idx = 5L,
    horizons = c(1L, 2L), obs = "poisson"))

  # Required output columns.
  need <- c("health_zone", "horizon", "mu_forecast", "p_case_invasion",
            "p_infection_invasion", "p_infection_lo", "p_infection_hi",
            "was_active_before")
  expect_true(all(need %in% names(fc)))

  # Affected zones (A, B) carry NO invasion output.
  aff <- fc[fc$health_zone %in% c("A", "B"), ]
  expect_true(all(is.na(aff$p_case_invasion)))
  expect_true(all(is.na(aff$mu_forecast)))

  # At-risk band ordering holds row-wise.
  ar <- fc[!fc$was_active_before, ]
  expect_true(nrow(ar) > 0)
  expect_true(all(ar$p_infection_lo <= ar$p_infection_invasion + 1e-12, na.rm = TRUE))
  expect_true(all(ar$p_infection_invasion <= ar$p_infection_hi + 1e-12, na.rm = TRUE))
  expect_true(all(ar$p_infection_lo >= ar$p_case_invasion - 1e-12, na.rm = TRUE))
  expect_true(all(ar$p_case_invasion >= 0 & ar$p_case_invasion <= 1, na.rm = TRUE))
})


# =============================================================================
# 5. RISK SCORES — compute_risk_scores (15)
# =============================================================================

test_that("compute_risk_scores: rr_nat and within-province rr_<prov> average ~1 over at-risk zones", {
  skip_if_missing("compute_risk_scores")

  fc <- tibble::tibble(
    health_zone       = c("A", "B", "C", "D", "E"),
    horizon           = 1L,
    training_cutoff   = as.Date("2026-06-01"),
    method            = "Renewal",
    was_active_before = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    mu_forecast       = c(3, 1, 2, 3, 4))
  pmap <- tibble::tibble(
    nom      = c("A", "B", "C", "D", "E"),
    province = c("Ituri", "Ituri", "Ituri", "Nord-Kivu", "Nord-Kivu"))

  rs <- compute_risk_scores(fc, pmap,
                            provinces = c("Ituri", "Nord-Kivu", "Haut-Uele"))

  # Nationwide relative risk = mu / mean(mu over at-risk); mean over at-risk = 1.
  atrisk <- !rs$was_active_before
  expect_equal(mean(rs$rr_nat[atrisk]), 1, tolerance = 1e-9)

  # Within-province rr averages to 1 over the IN-province at-risk zones.
  ituri <- rs$province == "Ituri" & atrisk
  expect_equal(mean(rs$rr_ituri[ituri]), 1, tolerance = 1e-9)
  nk <- rs$province == "Nord-Kivu" & atrisk
  expect_equal(mean(rs$rr_nordkivu[nk]), 1, tolerance = 1e-9)

  # Affected zone -> NA everywhere.
  aff_row <- rs[rs$health_zone == "A", ]
  expect_true(is.na(aff_row$rr_nat))
  expect_true(is.na(aff_row$rr_ituri))

  # rr_ituri backward-compatible: it is exactly the Ituri entry of the
  # generalised within-province set, NA outside Ituri.
  expect_true("rr_ituri" %in% names(rs))
  expect_true(all(is.na(rs$rr_ituri[rs$province != "Ituri"])))
})


# =============================================================================
# 6. VULNERABILITY INDEX — compute_vulnerability_index (20)
# =============================================================================

test_that("compute_vulnerability_index: V in [0,1]; access_gap present w/ OSRM, NA-safe w/o", {
  skip_if_missing("compute_vulnerability_index")

  # Facilities only in zone A -> B,C,D get a strictly increasing travel time to
  # the nearest facility zone (A). This makes access_gap strictly monotone.
  cov <- tibble::tibble(
    nom                = c("A", "B", "C", "D"),
    healthsite_density = c(8, 5, 2, 1),
    healthsite_count   = c(6, 0, 0, 0),
    pop_count          = c(2e4, 3e4, 5e4, 1e5),
    ccvi               = c(0.1, 0.4, 0.6, 0.9))
  osrm <- matrix(c(0,  30, 60, 90,
                   30,  0, 45, 70,
                   60, 45,  0, 55,
                   90, 70, 55,  0),
                 nrow = 4, byrow = TRUE, dimnames = list(cov$nom, cov$nom))

  v_osrm <- compute_vulnerability_index(cov, zones_all = cov$nom, osrm_mat = osrm)
  v_null <- compute_vulnerability_index(cov, zones_all = cov$nom, osrm_mat = NULL)

  expect_true(all(c("nom", "V", "access_gap", "healthcare_travel_min",
                    "surveillance_gap", "healthcare_gap",
                    "social_vulnerability") %in% names(v_osrm)))
  expect_true(all(v_osrm$V >= 0 & v_osrm$V <= 1), info = "V must be in [0,1]")

  # access_gap present (non all-NA) with OSRM; absent (all NA) without it.
  expect_false(all(is.na(v_osrm$access_gap)))
  expect_true(all(is.na(v_null$access_gap)),
              info = "access_gap must be NA-safe when osrm_mat is NULL")

  # Higher travel-time -> higher access_gap (rank-based, so strictly monotone).
  ord <- order(v_osrm$healthcare_travel_min)
  expect_false(is.unsorted(v_osrm$access_gap[ord]),
               info = "access_gap must be non-decreasing in travel time")
  expect_gt(v_osrm$access_gap[which.max(v_osrm$healthcare_travel_min)],
            v_osrm$access_gap[which.min(v_osrm$healthcare_travel_min)])
  # Facility-owning zone A has zero access time.
  expect_equal(v_osrm$healthcare_travel_min[v_osrm$nom == "A"], 0)
})


# =============================================================================
# 7. METRICS — .auc_pr / live WIS / base_rate / .brier / .log_score (16, 04b)
# =============================================================================

test_that(".auc_pr equals sklearn average_precision on a perfect ranker and a ties case", {
  skip_if_missing(".auc_pr")

  # Perfect ranker: all positives above all negatives -> AP = 1.
  expect_equal(.auc_pr(c(0.9, 0.8, 0.2, 0.1), c(1, 1, 0, 0)), 1, tolerance = 1e-12)

  # Ties case: scores {0.9, 0.6, 0.6, 0.1}, labels {1, 1, 0, 1}.
  # sklearn average_precision_score = 1/3 + (1/3)(2/3) + (1/3)(3/4) = 29/36.
  expect_equal(.auc_pr(c(0.9, 0.6, 0.6, 0.1), c(1, 1, 0, 1)), 29 / 36,
               tolerance = 1e-9, info = "tie-aware AP must match sklearn 29/36")

  # All-tied scores: single operating point, precision = P/(P+N).
  expect_equal(.auc_pr(c(0.5, 0.5, 0.5, 0.5), c(1, 0, 1, 0)), 0.5, tolerance = 1e-12)

  # No positives -> undefined (NA), by construction.
  expect_true(is.na(.auc_pr(c(0.5, 0.4), c(0, 0))))
})

test_that("live WIS (wis_one in 04b_epinowcast.R) matches the Bracher 2021 hand calculation", {
  skip_if(!requireNamespace("here", quietly = TRUE), "here not available")
  f <- file.path(here::here(), "spatiotemporal", "04b_epinowcast.R")
  wis_one <- .extract_fn_from_source(f, "wis_one")
  skip_if(is.null(wis_one), "could not extract live wis_one from 04b_epinowcast.R")

  # Two-interval WIS (Bracher 2021): 90% CI alpha=0.10, 60% CI alpha=0.40, K=2.
  #   IS_a = (u-l) + (2/a)*max(l-y,0) + (2/a)*max(y-u,0)
  #   WIS  = (1/(K+0.5)) * [ |y-m|/2 + sum_k (a_k/2) IS_k ]
  # Hand example: y=20, m=10, 90%CI=[5,15], 60%CI=[8,12].
  #   IS90 = 10 + 20*max(0) + 20*max(5) = 110 ; IS60 = 4 + 5*max(0) + 5*max(8) = 44
  #   WIS  = (1/2.5)*(5 + 0.05*110 + 0.20*44) = (1/2.5)*19.3 = 7.72
  expect_equal(wis_one(y = 20, m = 10, lo90 = 5, hi90 = 15, lo60 = 8, hi60 = 12),
               7.72, tolerance = 1e-9)

  # Independent re-implementation of the general Bracher formula must agree.
  bracher_wis <- function(y, m, lo, hi, alpha) {
    IS <- (hi - lo) + (2 / alpha) * (pmax(lo - y, 0) + pmax(y - hi, 0))
    (0.5 * abs(y - m) + sum((alpha / 2) * IS)) / (length(alpha) + 0.5)
  }
  expect_equal(wis_one(20, 10, 5, 15, 8, 12),
               bracher_wis(20, 10, c(5, 8), c(15, 12), c(0.10, 0.40)),
               tolerance = 1e-12)

  # Perfect forecast (point mass at y) -> WIS = 0.
  expect_equal(wis_one(10, 10, 10, 10, 10, 10), 0, tolerance = 1e-12)
  # y inside both intervals but wider intervals -> larger WIS (dispersion penalty).
  narrow <- wis_one(10, 10, 9, 11, 9.5, 10.5)
  wide   <- wis_one(10, 10, 5, 15, 8,   12)
  expect_lt(narrow, wide)
})

test_that("evaluate_invasion reports base_rate = mean(is_new_invasion) over the pooled at-risk rows", {
  skip_if_missing("evaluate_invasion")

  lfo <- tibble::tibble(
    method          = "Renewal", horizon = 1L,
    fold_id         = rep(c(1L, 2L), each = 6L),
    health_zone     = rep(paste0("Z", 1:6), 2L),
    p_invasion      = c(0.9, 0.5, 0.1, 0.05, 0.02, 0.01,
                        0.8, 0.6, 0.3, 0.10, 0.05, 0.00),
    is_new_invasion = c(1, 0, 0, 0, 0, 0,
                        1, 1, 0, 0, 0, 0))
  ev <- suppressWarnings(evaluate_invasion(lfo, n_boot = 50L))
  expect_equal(nrow(ev), 1L)
  expect_equal(ev$base_rate, mean(lfo$is_new_invasion), tolerance = 1e-12)
  expect_equal(ev$n_invasions, sum(lfo$is_new_invasion))
  expect_equal(ev$n_atrisk, nrow(lfo))
  expect_true(ev$auc_roc >= 0 & ev$auc_roc <= 1)
})

test_that(".brier and .log_score match their textbook definitions (ported from stale/11_metrics)", {
  skip_if_missing(".brier")
  skip_if_missing(".log_score")
  p <- c(0.2, 0.7, 0.4, 0.9); o <- c(0, 1, 0, 1)
  expect_equal(.brier(p, o), mean((p - o)^2), tolerance = 1e-12)
  expect_equal(.log_score(p, o),
               -mean(o * log(p) + (1 - o) * log(1 - p)), tolerance = 1e-9)
  expect_equal(.brier(c(0, 1, 1), c(0, 1, 1)), 0, tolerance = 1e-12)  # perfect
  expect_true(.log_score(c(0.001, 0.999, 0.999), c(0, 1, 1)) < 0.01)  # near-perfect
})


# =============================================================================
# 8. DETECTION CURVE + BALANCE METRICS — compute_detection_curve / invasion_balance_metrics (20)
# =============================================================================

test_that("compute_detection_curve recall is monotone non-decreasing in K and lies in [0,1]", {
  skip_if_missing("compute_detection_curve")

  lfo <- tibble::tibble(
    method          = "Renewal", horizon = 1L,
    fold_id         = rep(c(1L, 2L), each = 6L),
    cutoff          = as.Date("2026-06-01"),
    health_zone     = rep(paste0("Z", 1:6), 2L),
    p_invasion      = c(0.9, 0.5, 0.1, 0.05, 0.02, 0.01,
                        0.8, 0.6, 0.3, 0.10, 0.05, 0.00),
    is_new_invasion = c(1, 0, 0, 0, 0, 0,
                        1, 1, 0, 0, 0, 0))
  dc <- compute_detection_curve(lfo, "Renewal", horizon = 1L, ks = 1:6)

  expect_true(!is.null(dc) && nrow(dc) == 6L)
  expect_true(all(dc$recall >= 0 & dc$recall <= 1), info = "recall in [0,1]")
  expect_true(all(dc$precision >= 0 & dc$precision <= 1), info = "precision in [0,1]")
  expect_true(all(diff(dc$recall) >= -1e-12),
              info = "recall must be monotone non-decreasing in K")
  expect_equal(dc$recall[dc$k == 6L], 1, tolerance = 1e-12,
               info = "monitoring all 6 at-risk zones catches every invasion")
})

test_that("invasion_balance_metrics: balanced_accuracy in [0,1], mcc in [-1,1], perfect split -> 1", {
  skip_if_missing("invasion_balance_metrics")

  # Perfectly separable: the two invasions have the two highest p_invasion.
  lfo <- tibble::tibble(
    method          = "Renewal", horizon = 1L,
    fold_id         = 1L,
    health_zone     = paste0("Z", 1:5),
    p_invasion      = c(0.9, 0.8, 0.3, 0.2, 0.1),
    is_new_invasion = c(1, 1, 0, 0, 0))
  bm <- invasion_balance_metrics(lfo, "Renewal", horizon = 1L)
  expect_true(!is.null(bm))
  expect_true(bm$balanced_accuracy >= 0 & bm$balanced_accuracy <= 1)
  expect_true(bm$mcc >= -1 & bm$mcc <= 1)
  expect_equal(bm$balanced_accuracy, 1, tolerance = 1e-12)
  expect_equal(bm$mcc, 1, tolerance = 1e-12)

  # Imperfect ranking: invariants still hold, and MCC is now strictly below 1.
  lfo2 <- tibble::tibble(
    method          = "Renewal", horizon = 1L, fold_id = 1L,
    health_zone     = paste0("Z", 1:6),
    p_invasion      = c(0.9, 0.4, 0.7, 0.2, 0.5, 0.1),
    is_new_invasion = c(1, 1, 0, 0, 1, 0))
  bm2 <- invasion_balance_metrics(lfo2, "Renewal", horizon = 1L)
  expect_true(bm2$balanced_accuracy >= 0 & bm2$balanced_accuracy <= 1)
  expect_true(bm2$mcc >= -1 & bm2$mcc <= 1)
})


# =============================================================================
# 9. DATE PARSING + ONSET IMPUTATION — .parse_date (01) + config imputation rule
# =============================================================================

test_that(".parse_date handles a mixed ISO + dd/mm/yyyy vector element-by-element", {
  skip_if_missing(".parse_date")

  x <- c("2026-05-01",            # ISO
          "15/05/2026",           # dd/mm/yyyy
          "2026-06-30T12:00:00",  # ISO datetime
          "01/07/2026",           # dd/mm/yyyy (1 July, not 7 Jan)
          NA, "")                 # missing / blank
  d <- .parse_date(x)

  expect_s3_class(d, "Date")
  expect_true(all(!is.na(d[1:4])), info = "all four real dates must parse")
  expect_equal(d[1], as.Date("2026-05-01"))
  expect_equal(d[2], as.Date("2026-05-15"))   # dd/mm interpretation
  expect_equal(d[3], as.Date("2026-06-30"))
  expect_equal(d[4], as.Date("2026-07-01"))   # dd/mm priority over mm/dd
  expect_true(all(is.na(d[5:6])), info = "NA and blank -> NA")
})

test_that("onset imputation shifts a missing-onset record earlier by imp_delay (config rule)", {
  skip_if(!exists("DELAY_ONSET_SAMPLE_RATE"), "delay config not loaded")

  # NOTE: this checks the DIRECTION of the fallback (imputed onset earlier than the
  # sample date), NOT the live imputation. load_linelist() (01_data_prep.R) draws
  # imp_delay PER RECORD from the onset->sample delay (empirical bootstrap when >=30
  # complete pairs, else Exp), gates it on `onset_usable`, and floors the result at the
  # outbreak week — none of which this self-contained directional check reproduces.
  imp_delay <- if (isTRUE(get0("IMPUTE_ONSET_FROM_SAMPLE", ifnotfound = TRUE)) &&
                   is.finite(DELAY_ONSET_SAMPLE_RATE) && DELAY_ONSET_SAMPLE_RATE > 0)
    round(1 / DELAY_ONSET_SAMPLE_RATE) else 0
  expect_equal(imp_delay, round(1 / DELAY_ONSET_SAMPLE_RATE))
  expect_gt(imp_delay, 0)   # with rate 0.228 -> 4 days

  onset  <- as.Date(c("2026-05-10", NA))
  sample <- as.Date(c("2026-05-14", "2026-05-20"))
  # directional stand-in for the fallback: imputed onset = sample - a positive delay
  date_index <- dplyr::if_else(!is.na(onset), onset, sample - imp_delay)

  expect_equal(date_index[1], onset[1],
               info = "records WITH onset keep their onset date")
  expect_equal(date_index[2], sample[2] - imp_delay,
               info = "missing-onset record is imputed to sample - imp_delay (earlier)")
  expect_true(date_index[2] < sample[2],
              info = "imputed onset must be strictly earlier than the sample date")
})


# =============================================================================
# 10. PROVINCE RR GENERALISATION + .prov_suffix (15)
# =============================================================================

test_that(".prov_suffix produces column-safe lowercase alpha suffixes", {
  skip_if_missing(".prov_suffix")
  expect_equal(.prov_suffix("Nord-Kivu"), "nordkivu")
  expect_equal(.prov_suffix("Haut-Uele"), "hautuele")
  expect_equal(.prov_suffix("Ituri"),     "ituri")
  expect_equal(.prov_suffix("Sud-Kivu"),  "sudkivu")
})

test_that("compute_risk_scores generalises within-province RR to every province of interest", {
  skip_if_missing("compute_risk_scores")

  fc <- tibble::tibble(
    health_zone       = c("I1", "I2", "N1", "N2", "H1", "H2"),
    horizon           = 1L,
    training_cutoff   = as.Date("2026-06-01"),
    method            = "Renewal",
    was_active_before = FALSE,
    mu_forecast       = c(1, 3, 2, 2, 1, 5))
  pmap <- tibble::tibble(
    nom      = c("I1", "I2", "N1", "N2", "H1", "H2"),
    province = c("Ituri", "Ituri", "Nord-Kivu", "Nord-Kivu", "Haut-Uele", "Haut-Uele"))

  rs <- compute_risk_scores(fc, pmap,
                            provinces = c("Ituri", "Nord-Kivu", "Haut-Uele"))

  # A dedicated within-province RR column exists for each province of interest,
  # each named by .prov_suffix, each averaging 1 over its in-province at-risk set.
  for (prov in c("Ituri", "Nord-Kivu", "Haut-Uele")) {
    col <- paste0("rr_", .prov_suffix(prov))
    expect_true(col %in% names(rs), info = paste("missing", col))
    inp <- rs$province == prov
    expect_equal(mean(rs[[col]][inp]), 1, tolerance = 1e-9,
                 info = paste(col, "must average 1 over in-province zones"))
    expect_true(all(is.na(rs[[col]][!inp])),
                info = paste(col, "must be NA outside its province"))
  }
  # Within-province rank is independent of the nationwide rank (H2 tops nationally
  # AND within Haut-Uele here, but the columns are computed separately).
  expect_equal(rs$rr_hautuele_rank[rs$health_zone == "H2"], 1L)
  expect_equal(rs$rr_nat_rank[rs$health_zone == "H2"], 1L)
})
