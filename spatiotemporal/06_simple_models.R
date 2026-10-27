# =============================================================================
# 06_simple_models.R — Renewal force-of-infection CORE (+ inward-FOI kernel M11)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Provides the shared mobility-informed renewal primitives used by the workhorse
# (15_workhorse.R), the Bayesian suite (21_bayesian_renewal.R) and the cascade (31/32):
#   * weekly_censored_gt()          — WEEKLY generation-time pmf, double-interval censored at
#                                     the 7-day scale (lag 0 dropped, its mass carried)
#   * weekly_renewal_R_eff()        — daily-scale R -> the multiplier for that weekly kernel
#   * cascade_weekly_gt()           — the one weekly kernel the cascade uses
#   * daily_to_weekly_gt()          — weekly kernel from a daily pmf (weekly_censored_gt() from
#                                     the pmf's generating parameters; bins a bare pmf)
#   * compute_foi()                 — the import force Lambda = t(W) %*% Ytilde
#   * build_inward_contact_matrix() — the manuscript-motivated inward /
#                                     meeting-location FOI kernel (M11), Mills 2026
# (The exploratory S1-S5 renewal/logistic variants this file once held are retired;
#  S2/S3/S4 are documented as removed in METHODS §7.5.)
#
# Core principle: Lambda_i,t = sum_j W[j,i] * sum_k G_w[k] * Y_nc[j, t-k]
# where W is the outflow matrix (W[j,i] = fraction of j's outflow to i),
# and G_w is the weekly generation-time pmf (weekly_censored_gt()).
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
})


# ---------------------------------------------------------------------------
# Helper: weekly generation-time PMF
# ---------------------------------------------------------------------------

#' WEEKLY generation-time pmf under the SAME double-interval censoring as the daily pmf.
#'
#' WHY THIS REPLACES BINNING A DAILY PMF. Binning sums a daily pmf into 7-day blocks (what
#' daily_to_weekly_gt() did, and still does for a bare pmf). But the renewal equation here runs
#' on WEEKLY aggregates, so the infector's time is known only to the week, and both the primary
#' (infector) and secondary (infectee) event times are interval-censored at a 7-day window.
#' Binning a daily-censored pmf applies a 1-day primary window and then aggregates, which
#' measures the lag from the infector's exact DAY rather than from the start of its week.
#'
#' DIRECTION OF THE CORRECTION, measured rather than assumed. The weekly index is counted
#' from the infector's WEEK, which begins on average 3.5 days before the infector's own
#' infection time, so the correct weekly pmf is SHORTER than the naive binning, not longer.
#' On the medium profile (15.3 d): lag-1 mass 0.313 against 0.202, and mean lag 2.32 weeks
#' against 2.615. The renewal denominator is therefore weighted toward MORE RECENT weeks
#' than before, which raises R where incidence is falling and lowers it where it is rising.
#'
#' THE ESTIMAND. With the infector's infection time u ~ Uniform(0, 7) inside its week and a
#' generation interval D ~ Gamma(shape, rate), the infectee falls in week k when
#' floor((u + D)/7) = k, so
#'     P(k) = (1/7) * integral_0^7 [ F(7(k+1) - u) - F(7k - u) ] du,
#' which is exactly `primarycensored`'s primary-censored pmf with pwindow = swindow = 7,
#' evaluated on the day grid 0, 7, 14, ... This is the same construction
#' `.gt_daily_censored_pmf()` (02_epi_params.R) applies at pwindow = swindow = 1.
#'
#' Lag 0 is DROPPED and the pmf renormalised: the renewal sum runs over strictly past weeks,
#' because a lag-0 term would make a week's incidence depend on itself (see lag0_mass below),
#' matching .gweighted_own() (15_workhorse.R).
#'
#' @param mean,sd generation-interval mean and SD in DAYS (moment-matched Gamma).
#' @param n_weeks weekly lags to retain; default covers the 0.999 quantile.
#' @return numeric pmf over lags 1..n_weeks, summing to 1.
weekly_censored_gt <- function(mean, sd, n_weeks = NULL, week_len = 7L) {
  stopifnot(mean > 0, sd > 0, week_len >= 1)
  shape <- (mean / sd)^2
  rate  <- mean / sd^2
  if (is.null(n_weeks)) {
    q999   <- stats::qgamma(0.999, shape = shape, rate = rate)
    n_weeks <- max(2L, as.integer(ceiling((q999 + week_len) / week_len)))
  }
  ks <- 0:n_weeks                                   # include lag 0, dropped below
  pk <- NULL
  if (requireNamespace("primarycensored", quietly = TRUE)) {
    pk <- tryCatch(
      as.numeric(primarycensored::dprimarycensored(
        ks * week_len, pdist = stats::pgamma, shape = shape, rate = rate,
        pwindow = week_len, swindow = week_len, D = Inf)),
      error = function(e) NULL)
    if (!is.null(pk) && (!all(is.finite(pk)) || any(pk < -1e-12))) pk <- NULL
  }
  if (is.null(pk)) {
    # Exact base-R fallback, verified against primarycensored below in the tests.
    Fc <- function(q) stats::pgamma(pmax(q, 0), shape = shape, rate = rate)
    pk <- vapply(ks, function(k) {
      f <- function(u) Fc(week_len * (k + 1) - u) - Fc(week_len * k - u)
      stats::integrate(f, 0, week_len, subdivisions = 1000L,
                       rel.tol = 1e-10)$value / week_len
    }, numeric(1))
  }
  pk0 <- pmax(pk, 0); tot0 <- sum(pk0)
  pmf <- pk0[-1L]                                   # drop lag 0 (see lag0_mass below)
  tot <- sum(pmf)
  if (!is.finite(tot) || tot <= 0)
    stop("[weekly_censored_gt] degenerate pmf; check mean/sd.", call. = FALSE)
  pmf <- pmf / tot
  stopifnot(abs(sum(pmf) - 1) < 1e-9, all(pmf >= 0))
  # THE DROPPED LAG-0 MASS TRAVELS WITH THE KERNEL.
  #
  # A weekly lag of 0 means infector and infectee fall in the SAME calendar week, which happens
  # whenever the daily generation interval is short relative to where in the week the infector
  # sat. At the medium profile (15.3 d / 9.3 d) that is 5.81% of the mass (measured; a 4e6-draw
  # simulation of floor((U + D)/7) agrees to 3 decimal places). Dropping it and
  # renormalising is unavoidable for an explicit weekly recursion — the lag-0 term is
  # self-referential — but it is NOT free: it lengthens the effective generation interval, so
  # applying a DAILY-scale R (which EpiNow2 estimates, retaining that mass) to this kernel
  # understates growth. weekly_renewal_R_eff() below undoes exactly that, and needs this number.
  attr(pmf, "lag0_mass") <- if (is.finite(tot0) && tot0 > 0) pk0[1L] / tot0 else 0
  # The EXACT daily-R -> weekly-multiplier map for this kernel (weekly_renewal_R_eff()):
  # tabulated over the growth rate r, R(r) = (1 + r/rate)^shape (the Euler-Lotka relation of the
  # continuous Gamma generation interval) and M(r) = 1 / sum_{L>=1} pmf_L exp(-7 r L), the
  # multiplier with which this lag-0-free weekly kernel reproduces growth at rate r.
  attr(pmf, "R_eff_grid") <- local({
    r <- seq(-0.999 * rate, 2, length.out = 40001L)
    L <- seq_along(pmf)
    list(R = (1 + r / rate)^shape,
         M = 1 / as.numeric(exp(-week_len * outer(r, L)) %*% as.numeric(pmf)))
  })
  pmf
}

#' Daily-scale R -> the multiplier for a lag-0-dropped weekly kernel.
#'
#' The weekly kernel (weekly_censored_gt()) drops the weekly lag-0 mass G_0 (infector and infectee
#' in the same week; 5.8% at the medium profile) and renormalises, while R is EpiNow2's
#' DAILY-renewal reproduction number. The multiplier M(R) applied to that kernel must therefore
#' reproduce the growth R implies, not R itself.
#'
#' EXACT FORM (2026-10-15). Under a constant R the daily renewal grows at the rate r solving the
#' Euler-Lotka equation 1 = R * E[exp(-r T)], T the generation interval; for the Gamma generation
#' interval E[exp(-r T)] = (1 + r/rate)^(-shape). A weekly kernel reproduces that growth if
#' 1 = M * sum_{L>=1} Gt_L exp(-7 r L), so M = 1 / sum_{L>=1} Gt_L exp(-7 r L). The pair (R(r),
#' M(r)) is tabulated on the kernel (attr "R_eff_grid") and interpolated here. M is finite for
#' every R, equals R at R = 1, and agrees with the discretised daily generation-time pmf's
#' Euler-Lotka relation to within 0.2% for R up to 5.
#'
#' WHY NOT THE CLOSED FORM R (1 - G_0) / (1 - R G_0) that this used to be. It solves the weekly
#' recursion with the self-term Y_W = R (G_0 Y_W + sum_{L>=1} G_L Y_{W-L}), which treats every
#' same-week infection as able to seed further same-week infections with probability G_0 again,
#' without limit. Against the exact M it is +0.4% at R = 1.5, +4% at R = 3, +13% at R = 5, +64%
#' at R = 10, and it DIVERGES at R = 1/G_0 (17.2 at the medium profile), where the code then
#' capped its denominator at 1e-6 and returned R_eff of order 1e7. Discovery-phase back-test
#' rounds produced R draws in that range, so their two-week hazards were inflated a millionfold.
#'
#' @param R daily-scale reproduction number (scalar or vector), e.g. EpiNow2 posterior draws.
#' @param kernel the weekly kernel from weekly_censored_gt() (via daily_to_weekly_gt() or
#'   cascade_weekly_gt()), carrying attr "R_eff_grid". A kernel without the grid but with attr
#'   "lag0_mass" falls back to the closed form, returned as NA (not capped) where it is
#'   undefined; a kernel with neither (a plain binned pmf) returns R unchanged.
#' @return M(R), the multiplier to apply to the renormalised weekly kernel.
weekly_renewal_R_eff <- function(R, kernel) weekly_renewal_R_eff_fun(kernel)(R)

#' The map R -> M(R) of weekly_renewal_R_eff() as a FUNCTION of R, built once per kernel.
#'
#' For a caller that applies the same kernel many times (the cascade simulator, every simulated
#' week): stats::approx() re-takes the log of the 40,001-point grid and re-regularises it on every
#' call, which was ~75% of simulate_cascade()'s time. stats::approxfun() does that once; its
#' values are identical to approx()'s (same interpolation, same rule = 2).
weekly_renewal_R_eff_fun <- function(kernel) {
  grid <- attr(kernel, "R_eff_grid")
  if (!is.null(grid)) {
    f <- stats::approxfun(log(grid$R), log(grid$M), rule = 2)
    return(function(R) {
      R <- as.numeric(R)
      out <- rep(0, length(R))
      pos <- is.finite(R) & R > 0
      out[pos] <- exp(f(log(R[pos])))
      out[!is.finite(R)] <- NA_real_
      out
    })
  }
  function(R) .weekly_renewal_R_eff_closed(as.numeric(R), kernel)
}

# The closed form, for a kernel without the tabulated grid (see weekly_renewal_R_eff()).
.weekly_renewal_R_eff_closed <- function(R, kernel) {
  g0 <- attr(kernel, "lag0_mass")
  if (is.null(g0) || !length(g0) || !is.finite(g0) || g0 <= 0) return(R)
  den <- 1 - R * g0
  bad <- !is.finite(den) | den <= 0
  if (any(bad, na.rm = TRUE))
    warning(sprintf(paste0("[gt] %d reproduction number(s) at or above 1/G0 = %.2f on a kernel without ",
                           "generating parameters: the closed form is undefined there and returns NA."),
                    sum(bad, na.rm = TRUE), 1 / g0), call. = FALSE)
  out <- R * (1 - g0) / den
  out[bad] <- NA_real_
  out
}

#' The ONE weekly generation-time pmf the cascade uses.
#'
#' Both the R_eff estimator (31_source_dynamics.R) and the simulator (32_cascade_simulator.R)
#' must weight the renewal denominator with the SAME kernel, or R is estimated against one
#' infectiousness profile and deployed against another. This is the single entry point.
#'
#' Prefers the weekly double-interval-censored pmf built from the profile's mean/sd, which
#' matches the censoring of the weekly aggregates the renewal actually runs on. When the
#' profile is unavailable it falls back to the daily pmf through daily_to_weekly_gt(): a pmf
#' carrying its generating Gamma parameters (make_gt_pmf() attaches them) yields the SAME
#' censored kernel, so that is silent; a pmf without one (bare, or discretised "naive") is
#' BINNED, the older mis-censored construction and not an equivalent, and that case warns.
cascade_weekly_gt <- function(gt = get0("CASCADE_GT", ifnotfound = "medium"),
                              gt_pmfs = NULL,
                              profiles = get0("GT_PROFILES", ifnotfound = NULL)) {
  pr <- if (!is.null(profiles)) profiles[[gt]] else NULL
  if (!is.null(pr) && is.finite(pr$mean %||% NA_real_) && is.finite(pr$sd %||% NA_real_))
    return(weekly_censored_gt(pr$mean, pr$sd))
  if (!is.null(gt_pmfs) && !is.null(gt_pmfs[[gt]])) {
    if (!.gt_pmf_has_generator(gt_pmfs[[gt]]))
      warning(sprintf(paste0("[gt] no GT_PROFILES entry for '%s', and its daily pmf carries no ",
                             "censored generator (gt_mean/gt_sd, method 'censored'); falling back ",
                             "to BINNING it, which measures ",
                             "the lag from the infector's exact day rather than from its week ",
                             "and is not the same kernel."), gt),
              call. = FALSE)
    return(daily_to_weekly_gt(gt_pmfs[[gt]]))
  }
  stop(sprintf("[gt] cannot build a weekly generation time for '%s'.", gt), call. = FALSE)
}

#' Does a daily pmf carry the generating Gamma of the censored construction (make_gt_pmf())?
#' The ONE test daily_to_weekly_gt() dispatches on, shared with cascade_weekly_gt()'s warning so
#' the warning fires exactly when the kernel is binned.
.gt_pmf_has_generator <- function(gt_pmf_daily) {
  .m <- attr(gt_pmf_daily, "gt_mean"); .s <- attr(gt_pmf_daily, "gt_sd")
  .meth <- attr(gt_pmf_daily, "gt_method")
  !is.null(.m) && !is.null(.s) && length(.m) == 1L && length(.s) == 1L &&
    is.finite(.m) && is.finite(.s) && (is.null(.meth) || identical(.meth, "censored"))
}

#' Weekly generation-time pmf from a DAILY pmf.
#'
#' When the pmf carries its generating parameters (.gt_pmf_has_generator(); make_gt_pmf()
#' attaches them) this is weekly_censored_gt() on them, NOT a function of the daily values. A bare
#' pmf (no parameters, or one discretised "naive") is binned into 7-day blocks and renormalised.
#' @param gt_pmf_daily numeric vector, g[1..tau_max] (daily, sums to ~1)
#' @param n_weeks number of weekly lags to return (NULL: weekly_censored_gt()'s default, or
#'   ceiling(tau_max / 7) when binning)
#' @return normalised numeric vector G_weekly[1..n_weeks]
daily_to_weekly_gt <- function(gt_pmf_daily, n_weeks = NULL) {
  # CENSOR ONCE, AT THE WEEKLY SCALE. When the pmf carries its generating parameters, the weekly
  # kernel is rebuilt straight from the Gamma with a 7-day PRIMARY and SECONDARY window. That is
  # the estimand a weekly renewal needs: P(infectee's week - infector's week = k), with the
  # infector uniform within its own week. Summing a daily pmf into 7-day blocks answers a
  # different question — whether the DAY-lag falls in a 7-day band — and so discards where in its
  # week the infector sat. Measured on the shipped GT_PROFILES against 4e6 draws of
  # floor((u + D)/7): the weekly form matches to within 4e-4 (Monte-Carlo noise), while binning
  # understates lag-1 transmission by 36-43% relative (short 37.5%, medium 35.7%, long 42.5%)
  # and inflates the mean generation lag by 0.30-0.34 weeks.
  #
  # This is NOT censoring twice. It BYPASSES the daily pmf and returns to the continuous
  # distribution, so primary-censoring uncertainty is applied exactly once, at the scale the
  # model actually runs at. Chaining daily then weekly censoring would double-count it.
  if (.gt_pmf_has_generator(gt_pmf_daily) && exists("weekly_censored_gt", mode = "function"))
    return(weekly_censored_gt(attr(gt_pmf_daily, "gt_mean"), attr(gt_pmf_daily, "gt_sd"),
                              n_weeks = n_weeks))
  # Fallback: a pmf without a censored generator — no generating parameters (synthetic vectors
  # in the tests, a pmf built by hand or subset) or an explicitly naive one. Binning is retained
  # rather than refused so those callers keep working.
  tau_max <- length(gt_pmf_daily)
  weeks_back <- ceiling(tau_max / 7)
  if (!is.null(n_weeks)) weeks_back <- n_weeks

  G_weekly <- numeric(weeks_back)
  for (k in seq_len(weeks_back)) {
    d_lo <- 7 * (k - 1) + 1
    d_hi <- min(7 * k, tau_max)
    if (d_lo > tau_max) break
    G_weekly[k] <- sum(gt_pmf_daily[d_lo:d_hi])
  }
  G_total <- sum(G_weekly)
  if (G_total > 0) G_weekly <- G_weekly / G_total
  G_weekly
}

# ---------------------------------------------------------------------------
# Core: Force-of-Infection computation
# ---------------------------------------------------------------------------

#' Compute spatial force-of-infection for all zones at time t
#'
#' Lambda_i,t = sum_j W[j,i] * sum_k G_w[k] * Y_nc[j, t-k]
#'
#' W is an OUTFLOW matrix: W[j,i] is the fraction of j's outflow going to i.
#' Therefore, inflow TO zone i FROM zone j uses column j of t(W).
#' FOI = t(W) %*% Y_weighted   (each row i of t(W) gives weights from each j to i)
#'
#' @param Y_wide   numeric matrix (n_zones × n_weeks), zones in rows
#' @param W        outflow mobility matrix (n_zones × n_zones), rows = origins
#' @param G_weekly weekly GT PMF (normalised)
#' @param t_idx    week index to compute FOI for
#' @param zones_all canonical zone names (must match rownames of Y_wide and W)
#' @return named numeric vector of Lambda values, one per zone
compute_foi <- function(Y_wide, W, G_weekly, t_idx, zones_all) {
  n_zones    <- length(zones_all)
  weeks_back <- length(G_weekly)

  # Weighted sum of past case counts across all zones
  Y_weighted <- numeric(n_zones)
  for (k in seq_len(weeks_back)) {
    t_past <- t_idx - k
    if (t_past < 1) break
    Y_past       <- Y_wide[zones_all, t_past]
    Y_past[is.na(Y_past)] <- 0
    Y_weighted   <- Y_weighted + G_weekly[k] * Y_past
  }

  # Spatial spread: Lambda_i = sum_j W[j,i] * Y_weighted_j
  # t(W)[i,j] = W[j,i]: fraction of j's outflow that reaches i
  Lambda <- as.numeric(t(W[zones_all, zones_all]) %*% Y_weighted)
  names(Lambda) <- zones_all
  Lambda[Lambda < 0] <- 0
  Lambda
}

#' Manuscript-motivated inward / meeting-location effective-contact matrix.
#'
#' Builds an effective-contact matrix M such that the mobility-informed
#' inward force of infection is
#'     Lambda_i = sum_k M[k,i] * Ytilde_k = ( t(M) %*% Ytilde )_i ,
#' i.e. M can be dropped straight into compute_foi() as a "mobility" matrix with
#' NO other code change. This implements the frequency-dependent, TWO-SIDED
#' inward FOI of Mills (2026, "Multi-scale measures of time-varying epidemic
#' spread on human mobility networks"): susceptible residents of zone i and
#' infectious residents of zone k both move and can meet at a shared activity
#' ("meeting") location l, with transmission frequency-normalised by 1/N_eff(l).
#' The pipeline's default import force (Lambda = t(W) %*% Ytilde) instead moves
#' ONLY infecteds (one-sided) with no frequency normalisation.
#'
#' Presence matrix  P[j,l] = home*1{l==j} + (1-home)*W[j,l]  (row-stochastic:
#' fraction of zone-j residents' time at l). Effective population present at l is
#' N_eff(l) = sum_k P[k,l]*pop_k. Then
#'     C[i,k] = sum_l P[i,l]*P[k,l]/N_eff(l) = (P diag(1/N_eff) P^T)[i,k],
#' symmetric and non-negative, is the PER-SUSCEPTIBLE-RESIDENT contact rate.
#'
#' EXTENSIVITY. C is per person, so the expected number of introductions into zone i is
#' S_i * sum_k C[i,k]*Ytilde_k, and for a fully susceptible at-risk zone S_i ~ N_i. The
#' N_i factor must therefore be KEPT (it is a depletion factor that is unnecessary, not a
#' population factor that cancels): without it a small zone and a large zone with the same
#' per-person exposure would carry the same invasion hazard. The default one-sided import
#' force needs no such factor because t(W) %*% Ytilde already counts arriving infected
#' PERSONS. The returned matrix is M = C %*% diag(N), so that compute_foi()'s
#' Lambda = t(M) %*% Ytilde gives Lambda_i = N_i * sum_k C[i,k]*Ytilde_k as intended. M is
#' consequently NOT symmetric (C is); the k=i self term still vanishes for at-risk zones,
#' which have Ytilde_i = 0.
#'
#' @param W row-stochastic outflow mobility matrix (zero diagonal), zones x zones.
#' @param pop_vec named population vector.
#' @param zones_all canonical zone order.
#' @param home_fraction fraction of residents' time in the home zone
#'   (default MOBILITY_HOME_FRACTION; a documented, sensitivity-analysable
#'   assumption).
#' @return zones x zones effective-contact matrix (rows = source k, columns = receiving
#'   zone i, matching compute_foi()'s t(W) convention), rescaled to unit mean positive
#'   entry so the calibrated import coefficient beta stays O(1) (the absolute scale is
#'   not separately identifiable and is absorbed by beta).
build_inward_contact_matrix <- function(W, pop_vec, zones_all,
                                        home_fraction = get0("MOBILITY_HOME_FRACTION",
                                                             ifnotfound = 0.70)) {
  W <- W[zones_all, zones_all, drop = FALSE]
  W[!is.finite(W)] <- 0
  N <- as.numeric(pop_vec[zones_all])
  N[!is.finite(N) | N <= 0] <- stats::median(N[is.finite(N) & N > 0], na.rm = TRUE)
  # Presence matrix: home retention + between-zone activity, then renormalise
  # rows to exactly 1 (guards numerical drift / any nonzero W diagonal).
  P <- (1 - home_fraction) * W
  diag(P) <- diag(P) + home_fraction
  P <- P / pmax(rowSums(P), 1e-12)
  Neff <- as.numeric(t(P) %*% N)
  Neff[!is.finite(Neff) | Neff <= 0] <- stats::median(Neff[is.finite(Neff) & Neff > 0], na.rm = TRUE)
  C <- P %*% diag(1 / Neff) %*% t(P)
  C <- (C + t(C)) / 2                      # enforce exact symmetry (numerical)
  # Scale by the RESIDENT population of the receiving zone (see EXTENSIVITY above).
  # compute_foi() forms t(M) %*% Ytilde, so row i must be scaled in t(M), i.e. column i
  # of M: M = C %*% diag(N) has M[k,i] = C[k,i]*N_i = C[i,k]*N_i by symmetry of C.
  M <- C %*% diag(N)
  dimnames(M) <- list(zones_all, zones_all)
  mpos <- M[M > 0]; s <- if (length(mpos)) mean(mpos) else 1
  if (is.finite(s) && s > 0) M <- M / s    # unit mean positive entry
  M
}
