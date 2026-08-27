# =============================================================================
# 08_stochastic_seir.R — Stochastic Spatial SEIR Model (C6)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Compartment structure (discrete-time, daily or weekly):
#   S_{i,t+1} = S_it - new_infections_it
#   E_{i,t+1} = E_it + new_infections_it - new_E_to_I_it
#   I_{i,t+1} = I_it + new_E_to_I_it - new_I_out_it
#   R_{i,t+1} = R_it + new_I_to_R_it
#   D_{i,t+1} = D_it + new_I_to_D_it
#
# Spatial coupling: Phi_{i,t} = (1-m)*I_it + m * sum_j W[j,i]*I_jt
#   (local + imported infectious pressure; W is outflow, so W[j,i] = j→i flow)
#
# Deterministic: use expected values (continuous compartments)
# Stochastic: Binomial draws at each transition (integer compartments)
#
# Literature: CFR=26.4% (cfr_summary.csv); incubation mean=5.3 d (modelling choice near the
# lower BDBV survivor estimate; clinical features from Wamala/MacNeil 2010, EID 16(7)/16(12))
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---------------------------------------------------------------------------
# SEIR parameters
# ---------------------------------------------------------------------------

#' Derive weekly SEIR rate parameters from daily rates
#'
#' @param cfr             case-fatality ratio (0–1)
#' @param incubation_days mean incubation period (days)
#' @param mean_illness_days mean time from symptom onset to outcome (days)
#' @param ascertainment  fraction of true infections that are observed
#' @return named list of weekly transition probabilities
make_seir_params <- function(cfr               = CFR_POINT_ESTIMATE,
                             incubation_days   = INCUBATION_MEAN_DAYS,
                             mean_illness_days = 8.0,
                             ascertainment     = ASCERTAINMENT_NOMINAL) {
  # Daily rates (Exponential model: P(leave in 1 day) = 1 - exp(-rate))
  sigma_daily        <- 1 / incubation_days              # E → I
  gamma_d_daily      <- cfr / mean_illness_days           # I → D
  gamma_r_daily      <- (1 - cfr) / mean_illness_days     # I → R
  gamma_total_daily  <- gamma_d_daily + gamma_r_daily     # = 1/mean_illness_days

  # Weekly transition probabilities (discretised exponential, competing risks)
  # gamma_out_w is derived from the TOTAL daily rate, not the sum of marginal probs
  list(
    sigma_w           = 1 - exp(-sigma_daily        * 7),  # P(E→I in a week)
    gamma_d_w         = 1 - exp(-gamma_d_daily      * 7),  # P(I→D in a week, marginal)
    gamma_r_w         = 1 - exp(-gamma_r_daily      * 7),  # P(I→R in a week, marginal)
    gamma_out_w       = 1 - exp(-gamma_total_daily  * 7),  # P(leave I in a week, correct)
    gamma_total_daily = gamma_total_daily,                  # for FOI calculation
    cfr               = cfr,
    ascertainment     = ascertainment,
    incubation_days   = incubation_days,
    mean_illness_days = mean_illness_days
  )
}

# ---------------------------------------------------------------------------
# State initialisation from observed data
# ---------------------------------------------------------------------------

#' Initialise SEIR compartments from observed zone-week case counts
#'
#' Back-calculates I and E from observed Y_nc under the ascertainment model.
#' Y_{i,t} ≈ rho_i * new_I_{i,t} where new_I = sigma * E (new infectious)
#' Hence I ≈ Y / (rho * sigma_weekly) and E ≈ I (approximately)
#'
#' @param Y_nc_vec  named numeric vector of confirmed_nc for the latest week (zones)
#' @param Y_cumsum  named numeric vector of cumulative confirmed_nc per zone
#' @param pop_vec   named population vector
#' @param params    list from make_seir_params()
#' @return named list: S, E, I, R, D (all named numeric vectors, length = n_zones)
initialise_seir_state <- function(Y_nc_vec, Y_cumsum, pop_vec, params) {
  zones <- names(Y_nc_vec)
  N <- pop_vec[zones]; N[is.na(N)] <- median(pop_vec, na.rm = TRUE)
  rho <- params$ascertainment

  # Observed cases ≈ rho * sigma_w * E  →  E = Y / (rho * sigma_w)
  new_I_weekly <- pmax(Y_nc_vec / rho, 0)
  E <- pmax(new_I_weekly / params$sigma_w, 0)
  # At steady state: sigma_w * E = gamma_out * I  →  I = Y / (rho * gamma_out).
  # Use the competing-risks weekly outflow gamma_out_w (= 1 - exp(-(gamma_d+gamma_r)*7)),
  # NOT the sum of the marginal weekly probabilities gamma_d_w + gamma_r_w, which
  # over-counts the true outflow and understates initial I.
  I <- pmax(new_I_weekly / params$gamma_out_w, 0)
  # recovered + dead approximated from cumulative observed (rough)
  cum_infected <- pmax(Y_cumsum / rho, I + E)
  R <- pmax(cum_infected * (1 - params$cfr) - I - E, 0)
  D <- pmax(cum_infected * params$cfr, 0)
  # susceptible
  S <- pmax(N - E - I - R - D, 0)

  # Enforce S + E + I + R + D <= N. When de-ascertained cumulative infections
  # (E+I+R+D) already exceed N, S is 0 and the OTHER compartments must be scaled
  # down proportionally — recomputing S alone (the previous behaviour) left the
  # invariant violated and let an over-stated I export excess infectious pressure.
  total <- S + E + I + R + D
  overshoot <- total > N
  if (any(overshoot)) {
    sc <- N[overshoot] / total[overshoot]      # in (0,1]; total>N here so <1
    E[overshoot] <- E[overshoot] * sc
    I[overshoot] <- I[overshoot] * sc
    R[overshoot] <- R[overshoot] * sc
    D[overshoot] <- D[overshoot] * sc
    S[overshoot] <- pmax(N[overshoot] - E[overshoot] - I[overshoot] -
                           R[overshoot] - D[overshoot], 0)
  }

  list(S = S, E = E, I = I, R = R, D = D, N = N)
}

# ---------------------------------------------------------------------------
# Deterministic SEIR
# ---------------------------------------------------------------------------

#' Run deterministic spatial SEIR for n_weeks
#'
#' @param state   named list(S, E, I, R, D, N): named numeric vectors
#' @param R_by_zone named numeric vector of zone-level R(t)
#' @param W       outflow mobility matrix: W[j,i] = fraction of j to i
#' @param params  list from make_seir_params()
#' @param m       mixing fraction (0 = no spatial coupling, 1 = fully mixed)
#' @param n_weeks integer, number of weeks to project
#' @return named list of matrices: S, E, I, R, D, Y_obs (zones × weeks)
run_seir_deterministic <- function(state, R_by_zone, W, params,
                                   m = 0.1, n_weeks = 2) {
  zones  <- names(state$S)
  n      <- length(zones)
  N      <- state$N[zones]

  # Transition probabilities
  sigma_w   <- params$sigma_w
  gamma_d_w <- params$gamma_d_w
  gamma_r_w <- params$gamma_r_w
  rho       <- params$ascertainment

  # Output matrices
  S_mat <- E_mat <- I_mat <- R_mat <- D_mat <- Y_mat <-
    matrix(0, nrow = n, ncol = n_weeks, dimnames = list(zones, NULL))

  S_t <- state$S[zones]; E_t <- state$E[zones]
  I_t <- state$I[zones]; R_t <- state$R[zones]; D_t <- state$D[zones]

  # W_inflow[i,j] = W[j,i] = outflow from j to i
  W_in <- t(W[zones, zones])

  for (t in seq_len(n_weeks)) {
    # Spatial infectious pressure at zone i: Phi_i = (1-m)*I_i + m * sum_j W[j,i]*I_j
    Phi <- (1 - m) * I_t + m * as.numeric(W_in %*% I_t)

    # Force of infection: 1 - exp(-R * gamma_total_daily * Phi/N * 7)
    # Uses the correct continuous-time beta = R * gamma_total_daily
    foi_rate <- pmin(1 - exp(-R_by_zone[zones] * params$gamma_total_daily *
                               7 * Phi / pmax(N, 1)), 1)

    # Deterministic transitions (expected values). Split the I outflow by the
    # competing-risks weekly rate gamma_out_w, then partition into death/recovery
    # by the CFR — mirroring the stochastic path. Using the marginal weekly probs
    # gamma_r_w + gamma_d_w directly would over-count the outflow and inflate the
    # implied CFR to gamma_d_w/(gamma_d_w+gamma_r_w) != cfr.
    new_inf   <- S_t * foi_rate
    new_EtoI  <- E_t * sigma_w
    new_Iout  <- I_t * params$gamma_out_w
    new_ItoD  <- new_Iout * params$cfr
    new_ItoR  <- new_Iout - new_ItoD
    new_obs   <- new_EtoI * rho  # observed new cases ≈ rho * new infectious

    S_t <- pmax(S_t - new_inf, 0)
    E_t <- pmax(E_t + new_inf - new_EtoI, 0)
    I_t <- pmax(I_t + new_EtoI - new_ItoR - new_ItoD, 0)
    R_t <- R_t + new_ItoR
    D_t <- D_t + new_ItoD

    S_mat[, t] <- S_t; E_mat[, t] <- E_t; I_mat[, t] <- I_t
    R_mat[, t] <- R_t; D_mat[, t] <- D_t; Y_mat[, t] <- new_obs
  }

  list(S = S_mat, E = E_mat, I = I_mat, R = R_mat, D = D_mat, Y_obs = Y_mat)
}

# ---------------------------------------------------------------------------
# Stochastic SEIR (Binomial draws)
# ---------------------------------------------------------------------------

#' Run ONE stochastic spatial SEIR trajectory for n_weeks
#'
#' Uses Binomial draws for each transition to preserve integer compartments
#' and produce calibrated stochastic uncertainty.
#'
#' @return named list of matrices: S, E, I, R, D, Y_obs (zones × weeks)
run_seir_stochastic_single <- function(state, R_by_zone, W, params,
                                       m = 0.1, n_weeks = 2) {
  zones  <- names(state$S)
  n      <- length(zones)
  N      <- state$N[zones]

  sigma_w   <- params$sigma_w
  gamma_d_w <- params$gamma_d_w
  gamma_r_w <- params$gamma_r_w
  rho       <- params$ascertainment
  # Correct total weekly outflow probability (competing exponential risks)
  gamma_out_w <- params$gamma_out_w
  # P(death | leaving I) = cfr under exponential competing risks
  prob_death  <- params$cfr

  S_t <- round(state$S[zones]); E_t <- round(state$E[zones])
  I_t <- round(state$I[zones]); R_t <- round(state$R[zones])
  D_t <- round(state$D[zones])

  S_mat <- E_mat <- I_mat <- R_mat <- D_mat <- Y_mat <- Inf_mat <-
    matrix(0L, nrow = n, ncol = n_weeks, dimnames = list(zones, NULL))

  W_in <- t(W[zones, zones])

  for (t in seq_len(n_weeks)) {
    Phi <- (1 - m) * I_t + m * as.numeric(W_in %*% I_t)
    foi_rate <- pmin(1 - exp(-R_by_zone[zones] * params$gamma_total_daily *
                               7 * Phi / pmax(N, 1)), 1)

    # Stochastic transitions (Binomial draws)
    new_inf  <- rbinom(n, size = pmax(S_t, 0L), prob = foi_rate)
    new_EtoI <- rbinom(n, size = pmax(E_t, 0L), prob = sigma_w)
    new_Iout <- rbinom(n, size = pmax(I_t, 0L), prob = gamma_out_w)
    new_ItoD <- rbinom(n, size = new_Iout, prob = prob_death)
    new_ItoR <- new_Iout - new_ItoD
    # Observed: Binomial(new_EtoI, rho) — ascertainment of new infectious
    new_obs  <- rbinom(n, size = pmax(new_EtoI, 0L), prob = rho)

    S_t <- pmax(S_t - new_inf, 0L)
    E_t <- pmax(E_t + new_inf - new_EtoI, 0L)
    I_t <- pmax(I_t + new_EtoI - new_Iout, 0L)
    R_t <- R_t + new_ItoR
    D_t <- D_t + new_ItoD

    S_mat[, t] <- S_t; E_mat[, t] <- E_t; I_mat[, t] <- I_t
    R_mat[, t] <- R_t; D_mat[, t] <- D_t; Y_mat[, t] <- new_obs
    Inf_mat[, t] <- new_inf   # new infections (imports land here immediately)
  }

  # Inf_obs: new infections per week. Invasion (a zone becoming infected) is
  # observable in this compartment at h=1, whereas Y_obs (ascertained cases)
  # lags by the latent E->I->reporting chain and is structurally 0 in week 1 for
  # a freshly-seeded zone. Invasion probability must be read from Inf_obs.
  list(S = S_mat, E = E_mat, I = I_mat, R = R_mat, D = D_mat,
       Y_obs = Y_mat, Inf_obs = Inf_mat)
}

#' Run n_sim stochastic SEIR trajectories in parallel and summarise
#'
#' @param n_sim   number of simulations
#' @param seed    random seed
#' @return list(Y_sim = array [zones × weeks × n_sim], summary_tibble)
run_seir_stochastic <- function(state, R_by_zone, W, params,
                                m = 0.1, n_weeks = 2,
                                n_sim = N_SIMULATIONS, seed = RANDOM_SEED) {
  zones  <- names(state$S)
  n      <- length(zones)
  set.seed(seed)

  Y_arr <- array(0L, dim = c(n, n_weeks, n_sim),
                 dimnames = list(zones, NULL, NULL))
  Inf_arr <- array(0L, dim = c(n, n_weeks, n_sim),
                   dimnames = list(zones, NULL, NULL))

  for (s in seq_len(n_sim)) {
    traj <- run_seir_stochastic_single(state, R_by_zone, W, params,
                                       m = m, n_weeks = n_weeks)
    Y_arr[, , s]   <- traj$Y_obs
    Inf_arr[, , s] <- traj$Inf_obs
  }

  # Summarise. Count forecast (mu_forecast, quantiles) uses ascertained cases;
  # invasion probability uses *cumulative new infections* through horizon h, so a
  # zone seeded by imported infection registers invasion risk at h=1 rather than
  # waiting out the latent+reporting lag.
  results <- vector("list", n_weeks)
  for (h in seq_len(n_weeks)) {
    # CUMULATIVE ascertained cases over weeks 1..h (n_zones × n_sim), matching the
    # cumulative-window count target and every other model; SEIR's own p_invasion
    # (from inf_cum_h) is likewise cumulative. Marginal week-h counts would drop the
    # earlier weeks at h>=2 and be incoherent with the invasion probability.
    Y_h <- if (h == 1L) Y_arr[, 1L, ] else
      apply(Y_arr[, seq_len(h), , drop = FALSE], c(1, 3), sum)
    # n_sim == 1 drops the 3rd dim to a vector; the apply(Y_h, 1, quantile) calls
    # below need a matrix, so coerce (mirrors the inf_cum_h guard).
    if (!is.matrix(Y_h)) Y_h <- matrix(Y_h, nrow = n)
    inf_cum_h <- if (h == 1L) Inf_arr[, 1L, ] else
      apply(Inf_arr[, seq_len(h), , drop = FALSE], c(1, 3), sum)
    if (!is.matrix(inf_cum_h)) inf_cum_h <- matrix(inf_cum_h, nrow = n)
    mu_h  <- rowMeans(Y_h, na.rm = TRUE)
    p_inf <- rowMeans(inf_cum_h > 0, na.rm = TRUE)
    # Consistency: a zone with zero expected infections cannot have positive
    # invasion probability (Monte-Carlo noise otherwise leaves p>0 at mu=0).
    p_inf[rowMeans(inf_cum_h, na.rm = TRUE) < 1e-9] <- 0
    results[[h]] <- tibble(
      health_zone = zones,
      horizon     = h,
      mu_forecast = mu_h,
      # p_infection (SEIR models infection, not ascertained cases) — named to
      # match the case/infection distinction used across the suite.
      p_invasion  = p_inf,
      p_infection_invasion = p_inf,
      q05  = apply(Y_h, 1, quantile, 0.05,  na.rm = TRUE),
      q20  = apply(Y_h, 1, quantile, 0.20,  na.rm = TRUE),
      q25  = apply(Y_h, 1, quantile, 0.25,  na.rm = TRUE),
      q75  = apply(Y_h, 1, quantile, 0.75,  na.rm = TRUE),
      q80  = apply(Y_h, 1, quantile, 0.80,  na.rm = TRUE),
      q95  = apply(Y_h, 1, quantile, 0.95,  na.rm = TRUE),
      method = "C6_stochastic"
    )
  }

  list(
    summary = bind_rows(results),
    Y_sim   = Y_arr
  )
}

# ---------------------------------------------------------------------------
# Estimate mixing fraction m from training data
# ---------------------------------------------------------------------------

#' Estimate SEIR spatial mixing fraction by minimising forecast RMSE
#'
#' @param zone_week_nc tibble with health_zone, week_start, confirmed_nc
#' @param R_by_zone    named vector of zone R(t) estimates
#' @param W            outflow mobility matrix
#' @param pop_vec      named population vector
#' @param params       list from make_seir_params()
#' @param zones_all    canonical zone names
#' @param t_train_max  training week index (optimise over folds 2..t_train_max-1)
#' @return scalar estimate of m in [0.001, 0.5]
estimate_mixing_fraction <- function(zone_week_nc, R_by_zone, W, pop_vec,
                                     params, zones_all,
                                     t_train_max = NULL) {
  weeks <- sort(unique(zone_week_nc$week_start))
  if (is.null(t_train_max)) t_train_max <- length(weeks)
  if (t_train_max < 4) return(0.05)  # not enough data

  Y_mat <- matrix(0, length(zones_all), length(weeks),
                  dimnames = list(zones_all, as.character(weeks)))
  for (k in seq_along(weeks)) {
    sub <- zone_week_nc |> filter(week_start == weeks[k], health_zone %in% zones_all)
    Y_mat[sub$health_zone, k] <- ifelse(is.na(sub$confirmed_nc), 0, sub$confirmed_nc)
  }

  # Score the TWO-week-ahead prediction, not one week. At h=1 the deterministic
  # observed count is E*sigma_w*rho = Y exactly (by the initialisation), i.e.
  # independent of the mixing fraction m, so a one-week objective is flat in m and
  # optimise() just returns the interval end-point. The mobility-driven new
  # infections only reach the observed compartment at h=2, so m is identifiable
  # only from the two-week-ahead error.
  obj <- function(m) {
    sse <- 0
    taus <- if (t_train_max - 2 >= 2) 2:(t_train_max - 2) else integer(0)
    for (tau in taus) {
      Y_nc_vec <- Y_mat[, tau]
      Y_cumsum <- rowSums(Y_mat[, seq_len(tau), drop = FALSE])
      state_tau <- initialise_seir_state(Y_nc_vec, Y_cumsum, pop_vec, params)
      det_result <- tryCatch(
        run_seir_deterministic(state_tau, R_by_zone, W, params, m = m, n_weeks = 2),
        error = function(e) NULL
      )
      if (is.null(det_result) || ncol(det_result$Y_obs) < 2L) next
      Y_pred <- det_result$Y_obs[, 2]
      Y_true <- Y_mat[zones_all, min(tau + 2, ncol(Y_mat))]
      sse <- sse + sum((Y_true - Y_pred)^2, na.rm = TRUE)
    }
    sse
  }

  opt <- tryCatch(
    optimise(obj, interval = c(0.001, 0.5)),
    error = function(e) list(minimum = 0.05)
  )
  message(sprintf("[SEIR] Estimated mixing fraction m = %.4f", opt$minimum))
  opt$minimum
}

# ---------------------------------------------------------------------------
# Master SEIR forecast function
# ---------------------------------------------------------------------------

#' Run C6 SEIR model (deterministic + stochastic) for forecast
#'
#' @param zone_week_nc tibble (health_zone, week_start, confirmed_nc)
#' @param W            outflow mobility matrix
#' @param pop_vec      named population vector
#' @param Rt_national  national R(t) (scalar, posterior mean)
#' @param Rt_zone      named vector of zone R(t) (NAs allowed)
#' @param zones_all    canonical zone names
#' @param t_idx        training cutoff week index
#' @param horizons     forecast horizons (c(1L, 2L))
#' @param n_sim        stochastic simulations
#' @param seed         random seed
#' @param m            mixing fraction (if NULL, estimated)
#' @param mobility_id  label for output
#' @param gt_profile   label for output
#' @param training_cutoff Date of training cutoff
forecast_seir <- function(zone_week_nc, W, pop_vec,
                          Rt_national, Rt_zone = NULL,
                          zones_all, t_idx, horizons = LFO_HORIZONS,
                          n_sim = N_SIMULATIONS, seed = RANDOM_SEED,
                          m = NULL, params = NULL,
                          mobility_id = MOBILITY_PRIMARY,
                          gt_profile = GT_PRIMARY,
                          training_cutoff = NULL) {
  if (is.null(params)) params <- make_seir_params()

  # Fill zone R(t): use zone-specific where available, fallback to national
  R_vec <- setNames(rep(Rt_national, length(zones_all)), zones_all)
  if (!is.null(Rt_zone)) {
    matching <- intersect(names(Rt_zone), zones_all)
    R_vec[matching] <- ifelse(is.na(Rt_zone[matching]),
                              Rt_national, Rt_zone[matching])
  }
  R_vec <- pmax(R_vec, 0.01)

  weeks    <- sort(unique(zone_week_nc$week_start))
  t_use    <- min(t_idx, length(weeks))
  Y_mat    <- matrix(0, length(zones_all), length(weeks),
                     dimnames = list(zones_all, as.character(weeks)))
  for (k in seq_along(weeks)) {
    sub <- zone_week_nc |> filter(week_start == weeks[k], health_zone %in% zones_all)
    Y_mat[sub$health_zone, k] <- ifelse(is.na(sub$confirmed_nc), 0, sub$confirmed_nc)
  }

  Y_nc_now <- Y_mat[, t_use]
  Y_cumsum <- rowSums(Y_mat[, seq_len(t_use), drop = FALSE])

  state <- initialise_seir_state(Y_nc_now, Y_cumsum, pop_vec, params)

  # Estimate mixing fraction if not given
  if (is.null(m)) {
    m <- estimate_mixing_fraction(zone_week_nc, R_vec, W, pop_vec, params,
                                  zones_all, t_train_max = t_use)
  }

  message(sprintf("[SEIR] Forecasting with m=%.3f, R_mean=%.2f, n_sim=%d",
                  m, Rt_national, n_sim))

  # Run stochastic simulations
  n_ahead <- max(horizons)
  stoch_result <- run_seir_stochastic(state, R_vec, W, params,
                                      m = m, n_weeks = n_ahead,
                                      n_sim = n_sim, seed = seed)

  stoch_result$summary |>
    filter(horizon %in% horizons) |>
    mutate(
      mobility_id     = mobility_id,
      gt_profile      = gt_profile,
      training_cutoff = training_cutoff %||% weeks[t_use]
    )
}

# Null-coalescing operator (length-safe: a length-0 or length->1 LHS must not error inside
# `&&` under R >= 4.3, and a multi-element LHS coalesces only when fully present).
`%||%` <- function(a, b)
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
