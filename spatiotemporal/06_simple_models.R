# =============================================================================
# 06_simple_models.R — Renewal force-of-infection CORE (+ inward-FOI kernel M11)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Provides the shared mobility-informed renewal primitives used by the workhorse
# (15_workhorse.R) and the Bayesian suite (21_bayesian_renewal.R):
#   * daily_to_weekly_gt()          — weekly-aggregated generation-time PMF
#   * compute_foi()                 — the import force Lambda = t(W) %*% Ytilde
#   * build_inward_contact_matrix() — the manuscript-motivated inward /
#                                     meeting-location FOI kernel (M11), Mills 2026
# (The exploratory S1-S5 renewal/logistic variants this file once held are retired;
#  S2/S3/S4 are documented as removed in METHODS §7.5.)
#
# Core principle: Lambda_i,t = sum_j W[j,i] * sum_k G_w[k] * Y_nc[j, t-k]
# where W is the outflow matrix (W[j,i] = fraction of j's outflow to i),
# and G_w is the weekly-aggregated generation-time PMF.
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
})

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)[1]) a else b

# ---------------------------------------------------------------------------
# Helper: weekly generation-time PMF
# ---------------------------------------------------------------------------

#' Aggregate a daily GT PMF to weekly intervals
#' @param gt_pmf_daily numeric vector, g[1..tau_max] (daily, sums to ~1)
#' @param n_weeks number of weekly lags to return
#' @return normalised numeric vector G_weekly[1..n_weeks]
daily_to_weekly_gt <- function(gt_pmf_daily, n_weeks = NULL) {
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
#' Builds a SYMMETRIC effective-contact matrix M such that the mobility-informed
#' inward force of infection is
#'     Lambda_i = sum_k M[i,k] * Ytilde_k = ( t(M) %*% Ytilde )_i ,
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
#'     M[i,k] = sum_l P[i,l]*P[k,l]/N_eff(l) = (P diag(1/N_eff) P^T)[i,k],
#' symmetric and non-negative. For at-risk (fully susceptible) zones S_i ~ N_i so
#' no susceptible-depletion factor is needed (valid for the invasion task, and
#' the k=i self term vanishes because at-risk zones have Ytilde_i = 0).
#'
#' @param W row-stochastic outflow mobility matrix (zero diagonal), zones x zones.
#' @param pop_vec named population vector.
#' @param zones_all canonical zone order.
#' @param home_fraction fraction of residents' time in the home zone
#'   (default MOBILITY_HOME_FRACTION; a documented, sensitivity-analysable
#'   assumption).
#' @return symmetric zones x zones effective-contact matrix, rescaled to unit
#'   mean positive entry so the calibrated import coefficient beta stays O(1)
#'   (the absolute scale is not separately identifiable and is absorbed by beta).
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
  M <- P %*% diag(1 / Neff) %*% t(P)
  M <- (M + t(M)) / 2                      # enforce exact symmetry (numerical)
  dimnames(M) <- list(zones_all, zones_all)
  mpos <- M[M > 0]; s <- if (length(mpos)) mean(mpos) else 1
  if (is.finite(s) && s > 0) M <- M / s    # unit mean positive entry
  M
}
