# =============================================================================
# 44_reff_epinow2_check.R — is the cascade's national R anchor defensible?
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# WHAT THIS IS FOR — AND WHAT IT IS NOT. It is NO LONGER an independent validation of the
# cascade's anchor, and must not be quoted as one. The cascade now TAKES its national R from
# EpiNow2 (cascade_rt_draws(), 31_source_dynamics.R), so "EpiNow2 agrees with the anchor" would
# be a statement about the same model twice.
#
# What remains genuinely informative is a THREE-WAY sensitivity, because the numbers compared
# here are not all the same construction:
#
#   1. THE DEPLOYED ANCHOR. EpiNow2 via bayes_rt_week_draws(): dated by date_index exactly as the
#      weekly counts (epinow2_daily_confirmed(), 02_epi_params.R), censored as-of the analysis
#      date by linelist_observation_date(), and averaged over the estimation window WITHIN each
#      posterior draw. Its interval comes from those draws.
#   2. EpiNow2 ON THE FULL DAILY SERIES (estimate_rt_epinow2(), 02_epi_params.R), averaged over
#      the same calendar days. Same model and the SAME dating (epinow2_daily_confirmed()), but no
#      as-of censoring and a summary of daily posterior means, so a gap between 1 and 2 measures
#      how much the censoring and the summary move R. (Before RT_CACHE_VERSION 5 this series
#      dated onset-less cases at their sample date, which biased its R upward.)
#   3. THE RETIRED CONJUGATE ESTIMATOR, a Gamma-Poisson renewal ratio over ITS OWN window,
#      which drops the final week and is therefore NOT the deployed window.
#      Transparent and fast, but it treats the weekly generation-time PMF as known, models no
#      right-truncation, and allows no overdispersion. A gap between 3 and 1 is a statement
#      about the conjugate estimator's limitations, NOT evidence about EpiNow2.
#
# EpiNow2 has no import term. That is the right choice nationally, where mobility importation is
# internal redistribution rather than new infection, and it is precisely why this number is used
# only at national scale and never per zone.
#
# ALIGNMENT IS STILL THE WHOLE DIFFICULTY: all three are put on the same calendar window before
# being compared, because a window average and a final-day value are different estimands.
#
# COST. estimate_rt_epinow2() (02_epi_params.R) caches to
# outputs/diagnostics/epinow2_rt_*.rds keyed on the GT profile, analysis date and
# delay parameters, so the first run per frame pays one Stan fit and every later
# run is free.
#
# RUN:  Rscript 44_reff_epinow2_check.R
# Standalone, after run_cascade.R. Writes:
#   outputs/cascade/diagnostics/reff_national_epinow2_check.csv
#   outputs/cascade/diagnostics/reff_national_epinow2_check.json
# =============================================================================

ST_DIR <- if (nzchar(Sys.getenv("ST_DIR"))) Sys.getenv("ST_DIR") else getwd()
suppressPackageStartupMessages({ library(dplyr) })
# 22_daily_reissue.R is REQUIRED even though nothing here calls it directly: cascade_reff() now
# takes the cascade's R from the shared EpiNow2 posterior, and that path runs
# cascade_rt_draws() -> bayes_rt_week_draws() -> linelist_observation_date(), which lives there.
# It was absent until 2026-09-16 and this driver died four minutes in with "could not find
# function". Any entry point that builds its OWN source list — rather than going through
# run_cascade.R — has to carry it.
for (f in c("00_config.R", "01_data_prep.R", "02_epi_params.R", "06_simple_models.R",
            "15_workhorse.R", "20_forecast_detail.R", "21_bayesian_renewal.R",
            "22_daily_reissue.R",
            "30_projection_config.R", "31_source_dynamics.R", "32_cascade_simulator.R"))
  source(file.path(ST_DIR, f))

OUT_CASCADE <- file.path(OUT_DIR, "cascade")
DIAG_DIR    <- file.path(OUT_CASCADE, "diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# The national daily confirmed series EpiNow2 is fitted to
# ---------------------------------------------------------------------------

#' Daily confirmed cases nationally, dated exactly as the weekly case counts.
#'
#' Built by the SAME function run_all.R uses for its own R(t) panel and the Bayesian suite uses
#' for its R draws — epinow2_daily_confirmed() (02_epi_params.R): confirmed cases only, dated by
#' `date_index` (onset, or an onset imputed from the sample date), gapless, up to the analysis
#' date. Using a different series here would make the comparison a comparison of two data
#' definitions rather than of two estimators.
.rtc_daily_cases <- function(ll, analysis_date = ANALYSIS_DATE) {
  stopifnot(!is.null(ll))
  out <- epinow2_daily_confirmed(ll, analysis_date, caller = "rt-check")
  if (!nrow(out) || !sum(out$confirm))
    stop("[rt-check] no confirmed cases in the analysis window.", call. = FALSE)
  out
}

# ---------------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------------

#' Three-way sensitivity on the cascade's national R (see the file header: this is NOT an
#' independent validation, because EpiNow2 now supplies the deployed anchor).
#'
#' ALIGNMENT IS THE WHOLE DIFFICULTY, and getting it wrong would manufacture a
#' disagreement. The deployed anchor is a single number over a WINDOW; EpiNow2's diagnostic
#' R(t) is a daily curve. The like-for-like comparison is that curve averaged over exactly
#' the same calendar days — not its final day, which is both the least certain point and a
#' different estimand. Both are reported, with the window average as the headline.
#'
#' INTERVALS ARE NOT INTERCHANGEABLE. Each row carries the interval its OWN estimator
#' produces: the deployed anchor's comes from its posterior draws (`reff$rt_quantiles`), and
#' the conjugate diagnostic's from its Gamma posterior. Attaching the conjugate interval to
#' the EpiNow2 point estimate — which this function did while the cascade still used the
#' conjugate anchor and `reff$R_nat` meant something else — would quote an interval that no
#' model produced.
#'
#' @param reff the cascade_reff() object from the production run. Uses `R_nat` (deployed,
#'   EpiNow2), `rt_quantiles`, `R_nat_conjugate` (diagnostic) and `ncase`.
#' @param weeks the week_start dates of the count frame, so `reff$window_weeks` (column
#'   indices) can be turned into calendar dates.
cascade_reff_epinow2_check <- function(reff, weeks, daily_cases,
                                       gt_profile = CASCADE_GT,
                                       analysis_date = ANALYSIS_DATE,
                                       # OPTIONAL, and only for the truncation basis guard
                                       # (02_epi_params.R): it compares the recorded-onset
                                       # basis the truncation is fitted on against the
                                       # date_index basis this series uses. NULL simply
                                       # skips that check; the fit is unaffected.
                                       ll = NULL) {
  gt_p <- GT_PROFILES[[gt_profile]]
  if (is.null(gt_p)) stop("[rt-check] unknown GT profile '", gt_profile, "'.", call. = FALSE)
  rt <- tryCatch(estimate_rt_epinow2(daily_cases, gt_profile, gt_p, analysis_date, ll = ll),
                 error = function(e) NULL)
  if (is.null(rt) || !nrow(rt)) {
    warning("[rt-check] EpiNow2 produced no R(t); the cascade's anchor is reported without ",
            "an independent cross-check. This is a gap in the validation, not a failure of ",
            "the projection.", call. = FALSE)
    return(NULL)
  }
  rt$date <- as.Date(rt$date)
  # EpiNow2 PROJECTS PAST THE DATA by default, so the tail of this curve is a forecast, not
  # an estimate. Reporting its last row as "R(t), final day" would present a 7-day-ahead
  # projection as an independent measurement of the current anchor — which is exactly the
  # kind of quiet category error this cross-check exists to catch. Estimates are the rows at
  # or before the analysis date; the forecast tail is kept only to be labelled as such.
  .asof <- as.Date(analysis_date)
  rt_est <- rt[rt$date <= .asof, , drop = FALSE]
  if (!nrow(rt_est)) {
    warning("[rt-check] EpiNow2 returned no dates at or before the analysis date; nothing to ",
            "compare against.", call. = FALSE)
    return(NULL)
  }

  # THE DEPLOYED ANCHOR'S WINDOW, not the conjugate diagnostic's.
  # `reff$window_weeks` holds the CONJUGATE estimator's column indices, and that estimator
  # DROPS the final week (estimate_zone_reff(): nT_use <- nT - 1L) while the deployed EpiNow2
  # anchor deliberately INCLUDES it. Using the conjugate's indices here compared two windows
  # one week apart and reported the difference as a disagreement between estimators, under a
  # row literally labelled "over the same window" — on the shipped snapshot that manufactured
  # a ~12% gap (conjugate window 08-11..08-31, deployed 08-18..09-07; 1.116 vs 0.997).
  # The deployed window is defined by the anchor itself: RT_WINDOW_WEEKS ending at the last
  # week of the count frame, which is the forecast origin.
  .nw <- as.integer(get0("RT_WINDOW_WEEKS", ifnotfound = 1L))
  if (!is.finite(.nw) || .nw < 1L) .nw <- 1L
  .last_i <- length(weeks)
  wi <- seq.int(max(1L, .last_i - .nw + 1L), .last_i)
  win_start <- min(as.Date(weeks[wi]))
  win_end   <- max(as.Date(weeks[wi])) + 6L
  # The conjugate row is reported over ITS OWN window, which is a different span; the label
  # below says so rather than implying the two are matched.
  .cj <- as.integer(reff$window_weeks)
  .cj <- .cj[.cj >= 1L & .cj <= length(weeks)]
  conj_start <- if (length(.cj)) min(as.Date(weeks[.cj])) else as.Date(NA)
  conj_end   <- if (length(.cj)) max(as.Date(weeks[.cj])) + 6L else as.Date(NA)

  inwin <- rt_est$date >= win_start & rt_est$date <= win_end
  if (!any(inwin)) {
    warning(sprintf(paste0("[rt-check] EpiNow2's R(t) curve (%s to %s) does not overlap the ",
                           "cascade's estimation window (%s to %s), so the two cannot be ",
                           "compared on this frame. Most often this means the cached EpiNow2 ",
                           "fit is from an older analysis date."),
                    format(min(rt_est$date)), format(max(rt_est$date)),
                    format(win_start), format(win_end)), call. = FALSE)
    return(NULL)
  }
  w <- rt_est[inwin, , drop = FALSE]
  # Interval for the window average. The daily posteriors are strongly autocorrelated, so
  # averaging their half-widths does NOT shrink the interval by sqrt(n) and pretending it
  # does would overstate the agreement. The mean half-width is used instead — the honest
  # conservative choice without the full posterior draws to hand.
  hw90 <- mean((w$R_hi_90 - w$R_lo_90) / 2, na.rm = TRUE)
  en_win <- mean(w$R_mean, na.rm = TRUE)
  last <- rt_est[which.max(rt_est$date), , drop = FALSE]   # last OBSERVED day, not the forecast tail

  # THE DEPLOYED ANCHOR, with the interval its OWN posterior produces. reff$R_nat is now an
  # EpiNow2 posterior mean, so its interval must come from those draws (rt_quantiles, 5%/95%)
  # and NOT from the conjugate Gamma below. While the cascade used the conjugate anchor the two
  # were the same object; they are not any more, and splicing them would quote an interval no
  # model produced.
  r_nat  <- as.numeric(reff$R_nat)
  .rq    <- reff$rt_quantiles
  dep_lo <- if (is.null(.rq) || length(.rq) < 5L) NA_real_ else unname(as.numeric(.rq)[1])
  dep_hi <- if (is.null(.rq) || length(.rq) < 5L) NA_real_ else unname(as.numeric(.rq)[5])
  if (!is.finite(dep_lo) || !is.finite(dep_hi))
    warning("[rt-check] this reff carries no rt_quantiles, so the deployed anchor is reported ",
            "without its posterior interval. Rebuild it with the current cascade_reff().",
            call. = FALSE)

  # THE RETIRED CONJUGATE ESTIMATOR and ITS OWN interval: shape = kappa_nat * prior_mean +
  # total cases. This describes the conjugate number only; nothing in the projection uses it.
  r_conj <- suppressWarnings(as.numeric(reff$R_nat_conjugate))
  if (!length(r_conj)) r_conj <- NA_real_
  a_nat <- get0("CASCADE_R_KAPPA_NAT", ifnotfound = 1) *
             get0("CASCADE_R_PRIOR_MEAN", ifnotfound = 1) + sum(reff$ncase)
  sd_log_nat <- sqrt(trigamma(max(a_nat, 1e-6)))
  cas_lo <- r_conj * exp(-1.6448536 * sd_log_nat)
  cas_hi <- r_conj * exp( 1.6448536 * sd_log_nat)

  # AGREEMENT is judged by interval overlap, not by a tolerance on the ratio: two estimates
  # with honest intervals agree when those intervals are compatible, and a fixed percentage
  # band would call a precise 10% gap a disagreement while passing a vague 40% one. The
  # comparison that matters is now the DEPLOYED anchor against the full-series EpiNow2 fit —
  # same model, different input series — because that is the part which is not circular.
  overlap <- is.finite(dep_lo) && is.finite(dep_hi) &&
    max(dep_lo, en_win - hw90) <= min(dep_hi, en_win + hw90)

  tbl <- tibble::tibble(
    quantity = c("DEPLOYED anchor: EpiNow2, as-of censored, per-draw window mean",
                 "EpiNow2 on the full daily series, averaged over the DEPLOYED window",
                 "EpiNow2 on the full daily series, last observed day (not the forecast tail)",
                 "conjugate renewal estimator over ITS OWN window, which drops the final week (DIAGNOSTIC; unused)"),
    estimate = c(r_nat, en_win, last$R_mean, r_conj),
    lo_90    = c(dep_lo, en_win - hw90, last$R_lo_90, cas_lo),
    hi_90    = c(dep_hi, en_win + hw90, last$R_hi_90, cas_hi),
    window_start = as.Date(c(win_start, win_start, last$date, conj_start)),
    window_end   = as.Date(c(win_end, win_end, last$date, conj_end)),
    gt_profile = gt_profile, n_days = c(NA_integer_, nrow(w), 1L, NA_integer_),
    drives_projection = c(TRUE, FALSE, FALSE, FALSE))

  # The DIRECTION and SIZE of each gap, not just a pass/fail. Two intervals that overlap at a
  # single point are not the same finding as two that agree comfortably, and a reader needs to
  # know which way each estimate errs.
  #
  # BOTH comparisons are measured against the DEPLOYED anchor, and both use the deployed
  # anchor's OWN interval. `margin` previously mixed the conjugate interval into an overlap
  # test that is now about the deployed one — a genuine inconsistency once reff$R_nat stopped
  # meaning the conjugate number.
  ratio  <- en_win / r_nat                               # full-series EpiNow2 vs deployed
  ratio_conj <- if (is.finite(r_conj) && r_conj > 0) r_conj / r_nat else NA_real_
  margin <- if (is.finite(dep_lo) && is.finite(dep_hi))
    min(dep_hi, en_win + hw90) - max(dep_lo, en_win - hw90) else NA_real_   # <0 = disjoint
  side <- if (ratio > 1.02) "LOWER than" else if (ratio < 0.98) "HIGHER than" else "level with"
  note <- sprintf(paste0("DEPLOYED anchor (EpiNow2, as-of censored, per-draw window mean) = ",
                         "%.2f [%.2f, %.2f]; EpiNow2 on the FULL daily series over the same ",
                         "window = %.2f [%.2f, %.2f] (%d days, GT '%s'). Intervals %s (overlap ",
                         "width %.3f); the deployed anchor is %s that fit (ratio %.2f). These ",
                         "are the SAME model on DIFFERENT input series — onset-or-sample dating ",
                         "with as-of censoring versus the full daily series — so this gap ",
                         "measures the dating and censoring choices, NOT model agreement, and ",
                         "is not an independent validation of the anchor. EpiNow2 final day = ",
                         "%.2f [%.2f, %.2f] (the last day WITH DATA; the curve runs past the ",
                         "analysis date and that tail is a forecast, not an estimate). The ",
                         "RETIRED conjugate estimator gives %.2f [%.2f, %.2f] on the same ",
                         "window (ratio %.2f vs deployed) and drives nothing: its interval is a ",
                         "within-model Gamma-Poisson posterior on the case counts, treating the ",
                         "generation-time PMF as known and allowing no overdispersion, so it is ",
                         "NARROWER than any uncertainty a reader should attach to it."),
                  r_nat, dep_lo, dep_hi, en_win, en_win - hw90, en_win + hw90,
                  nrow(w), gt_profile,
                  if (isTRUE(overlap)) "overlap" else "DO NOT OVERLAP", margin,
                  side, ratio,
                  last$R_mean, last$R_lo_90, last$R_hi_90,
                  r_conj, cas_lo, cas_hi, ratio_conj)
  if (max(rt$date) > .asof)
    message(sprintf("[rt-check] dropped %d forecast day(s) after %s from the comparison.",
                    sum(rt$date > .asof), format(.asof)))
  if (!isTRUE(overlap))
    warning(sprintf(paste0("[rt-check] the DEPLOYED anchor and the full-series EpiNow2 fit ",
                           "DISAGREE: %s Since both are EpiNow2, a disagreement here points at ",
                           "the input series — the dating rule or the as-of censoring — rather ",
                           "than at the model, and should be resolved or reported before ",
                           "publication."), note), call. = FALSE)
  else
    message("[rt-check] ", note)

  # A bare overlap is reported as marginal when the overlap is a small fraction of the
  # narrower interval — "agree" should not cover two intervals that meet at a point.
  marginal <- isTRUE(overlap) && is.finite(margin) &&
    margin < 0.25 * min(dep_hi - dep_lo, 2 * hw90)
  if (marginal)
    # The width REPORTED must be the width the test above used — the deployed anchor's, not the
    # conjugate's. Printing one while testing the other is how a "marginal" flag ends up
    # justified by an interval that played no part in raising it.
    message(sprintf(paste0("[rt-check] the agreement is MARGINAL: the intervals overlap by only ",
                           "%.3f, against a narrower interval width of %.3f. Report the ratio, ",
                           "not the pass."), margin, min(dep_hi - dep_lo, 2 * hw90)))
  # ratio_conj, deployed_lo/hi and r_conj are part of the CONTRACT, not locals: the driver
  # serialises them, and a field referenced there but missing here would serialise as null
  # without any error to notice.
  list(table = tbl, agree = overlap, marginal = marginal, ratio = ratio,
       ratio_conj = ratio_conj, margin = margin, note = note, rt = rt,
       deployed = c(estimate = r_nat, lo_90 = dep_lo, hi_90 = dep_hi),
       conjugate = c(estimate = r_conj, lo_90 = cas_lo, hi_90 = cas_hi),
       window = c(win_start, win_end), gt_profile = gt_profile)
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
if (identical(environment(), globalenv()) && !interactive()) {
  CACHE <- Sys.getenv("CASCADE_LAYER_CACHE")
  if (!nzchar(CACHE) || !file.exists(CACHE))
    stop("[rt-check] set CASCADE_LAYER_CACHE to the layer written by run_cascade.R; ",
         "rebuilding the layer here would risk comparing against a different frame.",
         call. = FALSE)
  layer <- readRDS(CACHE)
  zones_all <- layer$zones_all
  province_map <- load_province_map()
  zone_province <- stats::setNames(province_map$province, province_map$nom)

  des <- build_invasion_design(layer$zone_week_nc, layer$mobility_matrices, layer$gt_pmfs,
           layer$covariates, layer$osrm_mat, zones_all, mob = CASCADE_KERNEL, gt = CASCADE_GT)
  fit <- fit_bayes_renewal(des, cov_spec = CASCADE_COV_SPEC, iter = 2000L, chains = 2L)
  delta_imp <- tryCatch(as.numeric(cascade_fit_delta(
                 cov_model = paste0("Bayes-", CASCADE_KERNEL, "-geo"))), error = function(e) 1)
  if (!is.finite(delta_imp) || delta_imp <= 0) delta_imp <- 1
  reff <- cascade_reff(layer$zone_week_nc, zones_all, layer, fit, des, delta = delta_imp,
                       zone_province = zone_province, gt = CASCADE_GT, kernel = CASCADE_KERNEL,
                       linelist = layer$ll, issue_date = ANALYSIS_DATE)

  weeks <- sort(unique(layer$zone_week_nc$week_start))
  dc <- .rtc_daily_cases(layer$ll)
  message(sprintf("[rt-check] national daily series: %d days, %d confirmed cases, %s to %s",
                  nrow(dc), sum(dc$confirm), format(min(dc$date)), format(max(dc$date))))

  # ll is passed for the truncation BASIS GUARD only (02_epi_params.R); the fit itself does
  # not need it here, since issue_date = NULL selects the extract regime and its archive panel.
  res <- cascade_reff_epinow2_check(reff, weeks, dc, ll = layer$ll)
  if (is.null(res)) {
    message("[rt-check] no comparison written.")
  } else {
    readr::write_csv(res$table, file.path(DIAG_DIR, "reff_national_epinow2_check.csv"))
    jsonlite::write_json(list(
        # Named so a reader cannot mistake this for an independent validation: it is a
        # sensitivity between two constructions of the same model, plus the retired estimator.
        is_independent_validation = FALSE,
        agree = res$agree, marginal = res$marginal,
        ratio_fullseries_vs_deployed = res$ratio,
        ratio_conjugate_vs_deployed = res$ratio_conj,
        overlap_width = res$margin, note = res$note,
        gt_profile = res$gt_profile,
        window_start = format(res$window[1]),
        window_end = format(res$window[2]),
        deployed_R_nat = as.numeric(reff$R_nat),
        deployed_lo_90 = if (is.null(reff$rt_quantiles)) NA_real_
                         else unname(as.numeric(reff$rt_quantiles)[1]),
        deployed_hi_90 = if (is.null(reff$rt_quantiles)) NA_real_
                         else unname(as.numeric(reff$rt_quantiles)[5]),
        deployed_source = reff$rt_source,
        conjugate_R_nat_diagnostic = as.numeric(reff$R_nat_conjugate),
        n_zones_with_cases = reff$n_with_cases,
        n_prior_dominated = reff$n_prior_dominated),
      file.path(DIAG_DIR, "reff_national_epinow2_check.json"),
      auto_unbox = TRUE, pretty = TRUE, digits = 6, null = "null")
    message("[rt-check] wrote reff_national_epinow2_check.{csv,json}")
  }
}
