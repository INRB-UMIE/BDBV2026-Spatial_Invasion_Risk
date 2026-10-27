# =============================================================================
# 45_ascertainment_figures.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
#
# DETECTED INVASION vs TRUE INVASION: spatial variation in ascertainment.
#
# The modelled outcome is the onset week of a health zone's FIRST CONFIRMED CASE.
# That is a *detected* invasion, not a latent epidemiological one. This module
# quantifies, from the study's own data, how far the two can diverge and whether
# the divergence is spatially structured in a way that could bias either the
# invasion labels or the out-of-sample evaluation.
#
#   Figure A1  The ascertainment landscape
#                A) national map of the surveillance gap, invaded zones outlined
#                B) national map of discarded (not-a-case) alerts per 100k, invaded zones
#                   outlined; zones with no laboratory-resolved alert, and zones whose
#                   resolved alerts were all confirmed, drawn as two separate fixed fills
#                C) province confirmed share of laboratory-resolved alerts, with Wilson
#                   intervals at CONF_LEVEL
#                D) Spearman correlation matrix of the exogenous proxies and epicentre
#                   travel time
#
#   Figure A2  Does thin surveillance delay the first confirmed case?
#                A) arrival time vs the zone's median onset -> confirmation delay
#                B) arrival-time ECDF by surveillance-gap tercile
#                (The added-variable fits of arrival time on each proxy over the
#                travel-time baseline are still computed and written to CSV and the
#                report, but are no longer plotted; see section 6.)
#
#   Figure A3  Are invaded and never-invaded zones comparable?
#                A) surveillance gap by invasion status (violin + box + points)
#                B) Cliff's delta forest, invaded vs never-invaded, every exogenous
#                   proxy and epicentre travel time
#                C) the same contrast WITHIN quartiles of epicentre travel time
#                D) logistic odds ratios per SD, unadjusted and adjusted for
#                   epicentre travel time
#
#   Figure A4  Does the forecast rank by surveillance?
#                A) leave-future-out rank vs surveillance gap (at-risk zone-folds)
#                B) rank of the truly-invaded zones vs their surveillance gap
#                C) within-stratum discrimination by surveillance stratum, both horizons,
#                   zone-cluster bootstrap, on the basis of the headline performance:
#                   scored WITHIN each forecast round over all rounds (AUC-PR
#                   skill and average precision, each with its random-ordering chance
#                   level, and the within-round AUC-ROC)
#                D) watch-list recall vs budget K by surveillance stratum
#
#   Figure A5  Province test positivity, detection opportunity and the ranking
#                A) share of at-risk zone-forecasts whose zone already had a discarded
#                   alert notified at the forecast origin, by rank band of the
#                   H_MAIN-week watch-list
#                B) province confirmed share of laboratory-resolved alerts AS KNOWN AT EACH
#                   FORECAST ORIGIN vs the mean rank of the province's at-risk zone-forecasts
#                C) rank distributions by province (ordered by that as-known share),
#                   not-invaded vs invaded zone-forecasts
#                D) the live H_MAIN-week watch-list (top TOP_N), coloured by whether the
#                   zone has ever had an alert laboratory-resolved
#   Figure A6  Does the (secondary) frailty form mark down zones surveillance was not watching?
#                A) share with a discarded alert notified by the origin, among the
#                   zone-forecasts the frailty moves down vs up, by fixed-beta rank band
#                B) the same share for the zone-forecasts the frailty pushes out of / brings
#                   into the fixed-beta top-K watch-list
#
# NOTHING IS REFITTED. Every quantity is derived from artifacts the pipeline has
# already written (the LFO-CV frame, the risk table, the arrival-predictor table,
# the covariate layers) plus the line list itself. The module writes figures and
# the CSV of every number plotted; it writes no pipeline artifact and overwrites
# no published output.
#
# Outputs -> outputs/key_outputs/ascertainment/{FigureA1..A6}.{pdf,png}
#            outputs/key_outputs/ascertainment/*.csv, ASCERTAINMENT_REPORT.md,
#            ascertainment_run_info.json
#            (outputs/key_outputs/ascertainment/panels/ is created, but the lettered
#            sub-panels reach it only when FIGURE_KEEP_ALL=1; see section 1)
#            Under FORECAST_SCALE=raw every path above moves to its raw/ sibling
#            (fs_out_dir()), where the figure gate writes no figure unless
#            FIGURE_KEEP_ALL=1.
# Run:  Rscript 45_ascertainment_figures.R   (from spatiotemporal/)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(sf); library(scales)
  library(forcats); library(purrr)
})
sf::sf_use_s2(TRUE)
options(dplyr.summarise.inform = FALSE)
for (pk in c("ggrepel", "here", "jsonlite")) {
  if (!requireNamespace(pk, quietly = TRUE)) stop("45_ascertainment_figures.R needs ", pk)
}

HERE <- file.path(here::here(), "spatiotemporal")
OUT  <- file.path(HERE, "outputs")

# 01_data_prep.R sources 00_config.R; 16_invasion_eval.R carries the SCORING
# functions (.auc_pr, .auc_roc) so the stratified skill in Figure A4 is computed
# by exactly the code that produced the headline numbers, not a second copy.
# 43_spread_kinematics.R is source-guarded (is_script_run) and contributes the
# geography helpers (.load_osrm_named, build_geo) that define the epicentre
# travel time used in Figure 1C, so the covariate here is the same covariate.
source(file.path(HERE, "01_data_prep.R"))
stop_if_data_definition_sensitivity("45_ascertainment_figures.R")   # hard-coded published paths
source(file.path(HERE, "16_invasion_eval.R"))
source(file.path(HERE, "43_spread_kinematics.R"))
source(file.path(HERE, "forecast_scale.R"))
stopifnot(exists("load_linelist"), exists("load_population"), exists("load_static_covariates"),
          exists(".auc_pr"), exists(".auc_roc"), exists("build_geo"),
          exists(".auc_pr_by_date", mode = "function"),
          exists("invasion_common_cells", mode = "function"),
          exists("invasion_is_burnin", mode = "function"),
          exists("EPICENTRE_ZONES"), exists("ANALYSIS_DATE"), exists("OUTBREAK_START"),
          exists("CONFIRMED_STATUS"), exists("DELAY_CAP_DAYS"), exists("RANDOM_SEED"),
          exists("lfo_origin", mode = "function"))

# fs_out_dir() (forecast_scale.R): a FORECAST_SCALE=raw run writes to the raw/ sibling instead of
# overwriting the primary outputs with a table and report built on the other scale.
FIG_DIR   <- fs_out_dir(file.path(OUT, "key_outputs", "ascertainment"))
PANEL_DIR <- file.path(FIG_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

CONF_LEVEL <- 0.95          # every interval in this module, stated once
# 95%, NOT the suite's 90% band convention (onset_mi_pooling.R, the manuscript figures). Every
# label and the report derive the level from CONF_LEVEL rather than typing it, so changing it
# here changes every interval AND every caption that states it.
N_BOOT     <- 2000L         # bootstrap replicates (zone-cluster where applicable)
ALPHA      <- 1 - CONF_LEVEL
Q_LO <- ALPHA / 2; Q_HI <- 1 - ALPHA / 2

# -----------------------------------------------------------------------------
# 1. DESIGN SYSTEM — identical tokens to make_manuscript_figures.R (no titles)
# -----------------------------------------------------------------------------
# House style (palette, province scale, themes, save_dual) from figure_style.R — the
# single definition for the whole suite. That save_dual() applies the FIGURE_KEEP gate
# (00_config.R), which lists only the five composites FigureA1-A5. The lettered sub-panels
# below (FigureA1a_..., FigureA4c_..., etc.) are working artefacts and are deliberately NOT
# listed, so the gate DROPS them: their save_dual() calls write nothing on a normal run.
# Only with FIGURE_KEEP_ALL=1 (gate disabled) are they written, to panels/ (PANEL_DIR,
# save_dual()'s default directory).
source(file.path(HERE, "figure_style.R"))

# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("45_ascertainment_figures.R")) {
  message("[ascert] sourced, not run: nothing is written. Use `Rscript 45_ascertainment_figures.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

INV_COL <- c("Invaded" = "#C0392B", "Never invaded" = "#4C78C8")

#' In-panel annotation anchored to a panel CORNER via +/-Inf, so the text never has to be
#' positioned in data coordinates (which silently moves when the data change) and works
#' identically on continuous and discrete scales.
#' Drawn as a borderless label on a translucent white ground rather than bare text, so an
#' annotation can never be rendered illegible by points or labels that happen to sit under it.
# The box is anchored flush (hjust 0 = left edge, 1 = right edge) and the inset from the panel
# edge is padding on that side: a fractional hjust such as -0.04 shifts each LINE of a multi-line
# label by that fraction of its own width, which left the lines of a block ragged.
ann_in <- function(label, x = -Inf, y = Inf, hjust = 0, vjust = 1.12, size = 2.45, inset_pt = 6) {
  pad <- if (hjust >= 0.5) margin(1.6, 1.6 + inset_pt, 1.6, 1.6, "pt")
         else              margin(1.6, 1.6, 1.6, 1.6 + inset_pt, "pt")
  annotate("label", x = x, y = y, label = label, hjust = hjust, vjust = vjust,
           size = size, colour = MUTED, family = base_family, lineheight = 1.08,
           fill = grDevices::adjustcolor("white", alpha.f = 0.82), linewidth = 0,
           label.r = unit(0, "pt"), label.padding = pad)
}


# -----------------------------------------------------------------------------
# 2. STATISTICAL HELPERS
#    Every estimator used anywhere below is defined once, here, and each returns
#    its own sample size so no panel can quietly plot a statistic computed on a
#    different support from the one its neighbours use.
# -----------------------------------------------------------------------------

#' Wilson score interval for a binomial proportion. Undefined (NA) at n = 0 —
#' never 0/0 silently rendered as zero.
# A VALIDATING, tibble-returning Wilson interval — NOT the canonical wilson_ci() in 00_config.R,
# which takes (k, n) and returns data.frame(p, lo, hi). Named differently so neither can shadow
# the other: this one screens out n <= 0 and x outside [0, n] row-wise, which the ascertainment
# panels need because their denominators are per-stratum counts.
.wilson_ci_tbl <- function(x, n, conf = CONF_LEVEL) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  out <- tibble::tibble(x = as.numeric(x), n = as.numeric(n),
                        est = NA_real_, lo = NA_real_, hi = NA_real_)
  ok <- is.finite(out$n) & out$n > 0 & is.finite(out$x) & out$x >= 0 & out$x <= out$n
  if (any(ok)) {
    p <- out$x[ok] / out$n[ok]; nn <- out$n[ok]
    den <- 1 + z^2 / nn
    ctr <- (p + z^2 / (2 * nn)) / den
    hw  <- z * sqrt(p * (1 - p) / nn + z^2 / (4 * nn^2)) / den
    out$est[ok] <- p
    out$lo[ok]  <- pmax(0, ctr - hw)
    out$hi[ok]  <- pmin(1, ctr + hw)
  }
  out
}

#' Spearman rank correlation with a Bonett-Wright confidence interval and the
#' asymptotic t-approximation p-value (exact = FALSE, so ties are permitted).
#' Complete cases only; the n returned is the number of pairs actually used.
spearman_ci <- function(x, y, conf = CONF_LEVEL) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 4L || stats::sd(x) == 0 || stats::sd(y) == 0)
    return(tibble::tibble(rho = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, n = n))
  rho <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  p   <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  # Bonett & Wright (2000): Fisher z on the rank correlation with an inflated SE.
  se <- sqrt((1 + rho^2 / 2) / (n - 3))
  z  <- atanh(pmin(pmax(rho, -0.999999), 0.999999))
  zc <- stats::qnorm(1 - (1 - conf) / 2)
  tibble::tibble(rho = rho, lo = tanh(z - zc * se), hi = tanh(z + zc * se), p = p, n = n)
}

#' Cliff's delta = P(a > b) - P(a < b), computed exactly from mid-ranks (so ties
#' contribute zero), with a percentile bootstrap interval that resamples the two
#' groups independently — the zones in the two groups are distinct units.
cliffs_delta <- function(a, b, n_boot = N_BOOT, conf = CONF_LEVEL, seed = RANDOM_SEED) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  na <- length(a); nb <- length(b)
  if (na < 2L || nb < 2L)
    return(tibble::tibble(delta = NA_real_, lo = NA_real_, hi = NA_real_, n_a = na, n_b = nb))
  .d <- function(u, v) {
    r <- rank(c(u, v), ties.method = "average")            # mid-ranks: ties weigh 0.5 each way
    Ua <- sum(r[seq_along(u)]) - length(u) * (length(u) + 1) / 2   # Mann-Whitney U for `u`
    2 * (Ua / (length(u) * length(v))) - 1
  }
  d0 <- .d(a, b)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set_fixed_seed(seed)
  bs <- vapply(seq_len(n_boot), function(i)
    .d(sample(a, na, replace = TRUE), sample(b, nb, replace = TRUE)), numeric(1))
  if (is.null(old)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
  else assign(".Random.seed", old, envir = globalenv())
  q <- stats::quantile(bs, c((1 - conf) / 2, 1 - (1 - conf) / 2), na.rm = TRUE, names = FALSE)
  tibble::tibble(delta = d0, lo = q[1], hi = q[2], n_a = na, n_b = nb)
}

#' Added-variable linear fit: what does `prox` add to a baseline regression of
#' `y` on `base`? Both models are fit on the SAME complete-case rows, so the
#' reported delta-R-squared is a like-for-like nested comparison rather than two
#' models fit on different supports.
#' Slope is reported PER SD of `prox` on those rows, so proxies in different
#' units are comparable on one axis.
av_slope <- function(y, base, prox, conf = CONF_LEVEL) {
  ok <- is.finite(y) & is.finite(base) & is.finite(prox)
  y <- y[ok]; base <- base[ok]; prox <- prox[ok]; n <- length(y)
  s <- stats::sd(prox)
  if (n < 10L || !is.finite(s) || s == 0)
    return(tibble::tibble(slope_sd = NA_real_, lo = NA_real_, hi = NA_real_,
                          p = NA_real_, dR2 = NA_real_, n = n))
  m0 <- stats::lm(y ~ base)
  m1 <- stats::lm(y ~ base + prox)
  cf <- summary(m1)$coefficients
  ci <- suppressMessages(stats::confint(m1, "prox", level = conf))
  tibble::tibble(slope_sd = unname(cf["prox", "Estimate"]) * s,
                 lo = unname(ci[1, 1]) * s, hi = unname(ci[1, 2]) * s,
                 p  = unname(cf["prox", "Pr(>|t|)"]),
                 dR2 = summary(m1)$r.squared - summary(m0)$r.squared,
                 n = n)
}

#' Logistic odds ratio per SD of `prox` for a binary outcome, either unadjusted
#' or adjusted for `base`. Profile-likelihood interval where it can be computed,
#' Wald otherwise; which one was used is returned, never hidden.
logit_or <- function(yb, prox, base = NULL, conf = CONF_LEVEL) {
  ok <- is.finite(prox) & !is.na(yb)
  if (!is.null(base)) ok <- ok & is.finite(base)
  yb <- as.integer(yb[ok]); prox <- prox[ok]
  base_v <- if (is.null(base)) NULL else base[ok]
  n <- length(yb); s <- stats::sd(prox)
  bad <- tibble::tibble(or_sd = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_,
                        n = n, n_pos = sum(yb), ci_type = NA_character_)
  if (n < 20L || !is.finite(s) || s == 0 || sum(yb) < 5L || sum(yb) == n) return(bad)
  dd <- if (is.null(base_v)) data.frame(yb = yb, prox = prox)
        else data.frame(yb = yb, prox = prox, base = base_v)
  fm <- if (is.null(base_v)) yb ~ prox else yb ~ base + prox
  m <- tryCatch(stats::glm(fm, data = dd, family = stats::binomial()), error = function(e) NULL)
  if (is.null(m) || !isTRUE(m$converged)) return(bad)
  cf <- summary(m)$coefficients
  ci <- tryCatch(suppressMessages(stats::confint(m, "prox", level = conf)), error = function(e) NULL)
  ci_type <- "profile"
  if (is.null(ci) || any(!is.finite(ci))) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    ci <- c(cf["prox", "Estimate"] - z * cf["prox", "Std. Error"],
            cf["prox", "Estimate"] + z * cf["prox", "Std. Error"])
    ci_type <- "wald"
  }
  tibble::tibble(or_sd = exp(unname(cf["prox", "Estimate"]) * s),
                 lo = exp(unname(ci[1]) * s), hi = exp(unname(ci[2]) * s),
                 p  = unname(cf["prox", "Pr(>|z|)"]), n = n, n_pos = sum(yb),
                 ci_type = ci_type)
}

#' Within-round AUC-ROC: the AUC-ROC of each forecast round (`date`), averaged over the rounds
#' that hold an invasion with each round weighted by its number of invasions E_f, as
#' .auc_pr_by_date() (16_invasion_eval.R) weights its per-round average precision. It is
#' therefore the mean, over the invasions, of the share of that round's non-invaded at-risk zones
#' ranked below the invaded one (ties counting half). Its random-ordering expectation is 0.5 in
#' every round, so it is 0.5 in every subset. A round whose rows are all invasions has no
#' ordering to score and drops out, as a round without an invasion does.
.auc_roc_by_date <- function(p, y, date) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y) & !is.na(date)
  p <- p[ok]; y <- y[ok]; date <- date[ok]
  if (!length(y)) return(NA_real_)
  s <- vapply(split(seq_along(y), date), function(i) {
    E <- sum(y[i]); if (E == 0 || E == length(i)) return(c(0, 0))
    c(E * .auc_roc(p[i], y[i]), E)
  }, numeric(2))
  if (sum(s[2, ]) == 0) NA_real_ else sum(s[1, ]) / sum(s[2, ])
}

#' Discrimination on ONE subset of at-risk zone-forecasts, on the basis of the headline
#' performance, with a zone-cluster bootstrap interval.
#'
#' THE HEADLINE BASIS (evaluate_invasion()'s auc_pr_skill_date, 16_invasion_eval.R): EVERY round,
#' as performance is reported (the selection reads the post-burn-in rounds only, a choice about
#' which model to feature, not about how it is scored), and every metric scored WITHIN each
#' forecast round and then combined by .auc_pr_by_date() (event-weighted). A metric pooled
#' across rounds ranks every round's rows in one list and so rewards a probability LEVEL that
#' happens to be high on the rounds with invasions, which is not the ordering an operator acts on;
#' that is why the pooled AUC-PR was demoted from the selection.
#'
#' FIVE quantities are returned, because skill alone is a trap. AUC-PR SKILL is the per-round
#' average precision over the per-round base rate (event-weighted). Its RANDOM-ORDERING value
#' (`skill_chance`, Bestgen's exact expectation, from .auc_pr_by_date()) is NOT 1 — with few
#' invasions per round a random list scores well above the base rate — and it depends on each
#' round's size and event count, so it differs between subsets. Its ATTAINABLE MAXIMUM
#' (`skill_ceiling`, every round's average precision equal to 1) is the inverse of the
#' event-weighted base rate, and it too differs between subsets. Average precision (`auc_pr_date`)
#' has the common ceiling 1 but a subset-specific chance level (`ap_chance`); the within-round
#' AUC-ROC (`auc_roc_date`) is the one with a common null (0.5) and a common ceiling (1).
#'
#' SCORED ON THE RAW PROBABILITY (`p_invasion_raw`, kept by fs_apply_lfo_scale()), the column
#' evaluate_invasion() scores. A within-round metric is unchanged by the recalibration, whose
#' prequential delta is ONE monotone map per round, so this is the same number on either scale;
#' the raw column is used so that it is the evaluation's own number, ties included.
strat_discrimination <- function(d, n_boot = N_BOOT, conf = CONF_LEVEL, seed = RANDOM_SEED) {
  pcol <- "p_invasion_raw"
  stopifnot(all(c(pcol, "fold_id", "cutoff", "health_zone", "is_new_invasion") %in% names(d)))
  d <- d[is.finite(d[[pcol]]) & !is.na(d$is_new_invasion), , drop = FALSE]
  y <- as.integer(d$is_new_invasion); p <- d[[pcol]]; f <- d$fold_id
  n <- length(y); n_pos <- sum(y)
  # n_pos counts positive ROWS (the pipeline's convention, and the denominator of the
  # pooled recall); n_pos_zones counts the DISTINCT invaded zones behind them. At the
  # 2-week horizon the two differ because consecutive outcome windows overlap, and only
  # the second is an effective sample size.
  n_pos_zones <- dplyr::n_distinct(d$health_zone[y == 1])
  # `prob_column` records which probability column was ACTUALLY scored, so the published table
  # carries its own scale rather than inheriting the LFO frame's (lfo_probability_column).
  out <- tibble::tibble(prob_column = pcol, rounds = "all",
                        n_rounds = dplyr::n_distinct(f), n_rounds_with_invasion = 0L,
                        n = n, n_pos = n_pos, n_pos_zones = n_pos_zones,
                        base_rate = if (n > 0) mean(y) else NA_real_,
                        auc_pr_date = NA_real_, ap_lo = NA_real_, ap_hi = NA_real_,
                        ap_chance = NA_real_,
                        auc_pr_skill_date = NA_real_, skill_lo = NA_real_, skill_hi = NA_real_,
                        skill_chance = NA_real_, skill_ceiling = NA_real_,
                        auc_roc_date = NA_real_, roc_lo = NA_real_, roc_hi = NA_real_,
                        n_boot_used = 0L, n_zones = dplyr::n_distinct(d$health_zone))
  if (n == 0 || n_pos == 0) return(out)
  a0 <- .auc_pr_by_date(p, y, f)
  out$n_rounds_with_invasion <- as.integer(a0[["n_dates"]])
  out$auc_pr_date       <- a0[["ap"]]
  out$auc_pr_skill_date <- a0[["skill"]]
  out$skill_chance      <- a0[["chance"]]
  # ap / skill is the event-weighted base rate sum_f E_f b_f / sum_f E_f (16_invasion_eval.R), so
  # the chance level and the ceiling convert between the two scales exactly. ap > 0 whenever a
  # round holds an invasion (the precision at that invasion's operating point is positive).
  out$ap_chance         <- a0[["chance"]] * a0[["ap"]] / a0[["skill"]]
  out$skill_ceiling     <- a0[["skill"]] / a0[["ap"]]
  out$auc_roc_date      <- .auc_roc_by_date(p, y, f)
  zones <- unique(d$health_zone)
  if (n_pos >= 2L && length(zones) > 5L) {
    idx_by_zone <- split(seq_len(nrow(d)), factor(d$health_zone, levels = zones))
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    set_fixed_seed(seed)
    bs <- matrix(NA_real_, nrow = n_boot, ncol = 3)
    for (b in seq_len(n_boot)) {
      zs  <- sample.int(length(zones), length(zones), replace = TRUE)
      idx <- unlist(idx_by_zone[zs], use.names = FALSE)
      yb <- y[idx]
      if (sum(yb) >= 1L && sum(yb) < length(yb)) {
        # The replicate is scored exactly as the point estimate is: within its rounds, SKILL
        # against the replicate's own per-round base rates. Average precision is recorded
        # separately, because skill quantiles times a base rate are NOT the average-precision
        # quantiles (the divisor varies replicate to replicate).
        ab <- .auc_pr_by_date(p[idx], yb, f[idx])
        bs[b, 1] <- ab[["skill"]]
        bs[b, 2] <- .auc_roc_by_date(p[idx], yb, f[idx])
        bs[b, 3] <- ab[["ap"]]
      }
    }
    if (is.null(old)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
    else assign(".Random.seed", old, envir = globalenv())
    used <- sum(is.finite(bs[, 1]))
    if (used >= 50L) {
      .q <- function(v) stats::quantile(v, c((1 - conf) / 2, 1 - (1 - conf) / 2),
                                        na.rm = TRUE, names = FALSE)
      qs <- .q(bs[, 1]); qr <- .q(bs[, 2]); qa <- .q(bs[, 3])
      out$skill_lo <- qs[1]; out$skill_hi <- qs[2]
      out$roc_lo   <- qr[1]; out$roc_hi   <- qr[2]
      out$ap_lo    <- qa[1]; out$ap_hi    <- qa[2]
    }
    out$n_boot_used <- as.integer(used)
  }
  out
}

#' Pooled top-K recall with a ZONE-CLUSTER bootstrap interval.
#'
#' The point estimate is the pipeline's pooled estimator (hits / invasion rows, which
#' compute_detection_curve reports as `recall_pooled`). The INTERVAL may not be a Wilson
#' interval on those rows: at the 2-week horizon the outcome windows of consecutive
#' forecast origins overlap, so one invaded zone contributes up to two positive rows and
#' the rows are not independent Bernoulli trials. Resampling whole ZONES is the same
#' correlated unit evaluate_invasion() bootstraps over.
recall_curve_boot <- function(ev_rows, ks, n_boot = N_BOOT, conf = CONF_LEVEL, seed = RANDOM_SEED) {
  rk <- ev_rows$rank_max
  zones <- unique(ev_rows$health_zone)
  out <- tibble::tibble(k = ks,
                        caught = vapply(ks, function(k) sum(rk <= k), numeric(1)),
                        events = nrow(ev_rows), n_zones = length(zones),
                        recall = vapply(ks, function(k) mean(rk <= k), numeric(1)),
                        lo = NA_real_, hi = NA_real_)
  if (length(zones) < 3L || nrow(ev_rows) < 3L) return(out)
  idx_by_zone <- split(seq_len(nrow(ev_rows)), factor(ev_rows$health_zone, levels = zones))
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set_fixed_seed(seed)
  bs <- matrix(NA_real_, nrow = n_boot, ncol = length(ks))
  for (b in seq_len(n_boot)) {
    idx <- unlist(idx_by_zone[sample.int(length(zones), length(zones), replace = TRUE)],
                  use.names = FALSE)
    bs[b, ] <- colMeans(outer(rk[idx], ks, "<="))
  }
  if (is.null(old)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
  else assign(".Random.seed", old, envir = globalenv())
  out$lo <- apply(bs, 2L, stats::quantile, probs = (1 - conf) / 2, na.rm = TRUE, names = FALSE)
  out$hi <- apply(bs, 2L, stats::quantile, probs = 1 - (1 - conf) / 2, na.rm = TRUE, names = FALSE)
  out
}

#' Share of rows carrying `flag` within each level of `group`, with a ZONE-CLUSTER bootstrap
#' interval.
#'
#' The rows are at-risk ZONE-FORECASTS, and one zone contributes one row per forecast round in
#' which it is at risk, with a flag that is persistent across its rounds — a Wilson interval on
#' the rows treats them as independent Bernoulli trials and is too narrow by roughly the square
#' root of that multiplicity. Whole zones are resampled, as recall_curve_boot() and
#' strat_discrimination() do. All groups share each resample, and a group's share in a resample
#' is its flagged rows over its rows among the drawn zones.
#'
#' AT A BOUNDARY (share 0 or 1) every resample returns the same value, so the percentile
#' interval has zero width, which is no statement of uncertainty. There the interval is the
#' Wilson interval on the group's number of DISTINCT ZONES — the effective sample size if each
#' zone's rounds were perfectly correlated, i.e. the conservative end of the clustering — and
#' `interval` names which construction each row carries.
share_by_group_boot <- function(d, group, flag, n_boot = N_BOOT, conf = CONF_LEVEL,
                                seed = RANDOM_SEED) {
  g  <- d[[group]]
  lv <- if (is.factor(g)) levels(g) else sort(unique(g))
  zones <- unique(d$health_zone)
  zi <- factor(d$health_zone, levels = zones); gi <- factor(as.character(g), levels = as.character(lv))
  f  <- as.numeric(d[[flag]] %in% TRUE)
  Nz <- unclass(table(zi, gi))                           # zones x groups: rows
  Xz <- unclass(tapply(f, list(zi, gi), sum)); Xz[is.na(Xz)] <- 0   # ... flagged rows
  n  <- colSums(Nz); x <- colSums(Xz)
  n_zones <- colSums(Nz > 0)
  out <- tibble::tibble(group = lv, n = as.integer(n), flagged = as.integer(x),
                        n_zones = as.integer(n_zones),
                        est = ifelse(n > 0, x / n, NA_real_), lo = NA_real_, hi = NA_real_,
                        interval = NA_character_)
  if (length(zones) >= 3L) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    set_fixed_seed(seed)
    Mz <- vapply(seq_len(n_boot), function(b)
      tabulate(sample.int(length(zones), length(zones), replace = TRUE), nbins = length(zones)),
      numeric(length(zones)))                            # zones x resamples: multiplicities
    if (is.null(old)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
    else assign(".Random.seed", old, envir = globalenv())
    S <- crossprod(Xz, Mz) / crossprod(Nz, Mz)           # groups x resamples (NaN: group not drawn)
    S[!is.finite(S)] <- NA_real_
    out$lo <- apply(S, 1L, stats::quantile, probs = (1 - conf) / 2, na.rm = TRUE, names = FALSE)
    out$hi <- apply(S, 1L, stats::quantile, probs = 1 - (1 - conf) / 2, na.rm = TRUE, names = FALSE)
    out$interval <- ifelse(out$n > 0, "zone-cluster bootstrap", NA_character_)
    # The replicate shares (groups x resamples, rows in `group` order), for contrasts between
    # groups taken on the SAME resamples (share_diff_boot()).
    attr(out, "draws") <- S
  }
  bnd <- is.finite(out$est) & (out$est == 0 | out$est == 1) & out$n_zones > 0
  if (any(bnd)) {
    w <- .wilson_ci_tbl(out$est[bnd] * out$n_zones[bnd], out$n_zones[bnd], conf = conf)
    out$lo[bnd] <- w$lo; out$hi[bnd] <- w$hi
    out$interval[bnd] <- "Wilson on distinct zones (boundary share)"
  }
  out
}

#' The difference in share between two groups of share_by_group_boot()'s output (a minus b),
#' with the interval taken on the SAME zone resamples, so the zones the two groups share are
#' resampled together. NA interval when the bootstrap was not run.
share_diff_boot <- function(bw, a, b, conf = CONF_LEVEL) {
  ia <- match(a, bw$group); ib <- match(b, bw$group)
  stopifnot(!is.na(ia), !is.na(ib))
  S <- attr(bw, "draws")
  est <- bw$est[ia] - bw$est[ib]
  if (is.null(S)) return(tibble::tibble(est = est, lo = NA_real_, hi = NA_real_))
  d <- S[ia, ] - S[ib, ]
  tibble::tibble(est = est,
                 lo = stats::quantile(d, (1 - conf) / 2, na.rm = TRUE, names = FALSE),
                 hi = stats::quantile(d, 1 - (1 - conf) / 2, na.rm = TRUE, names = FALSE))
}

#' Expand a table of (x, y) into the vertex sequence of a `geom_step(direction = "hv")`
#' path, so a confidence RIBBON can be drawn along the same staircase as the line it
#' bounds. Drawing the ribbon with plain `geom_ribbon` interpolates diagonally between
#' consecutive x values, i.e. it draws uncertainty at watch-list sizes that do not exist.
step_path <- function(d, xcol = "k") {
  d <- d[order(d[[xcol]]), , drop = FALSE]
  n <- nrow(d)
  if (n < 2L) return(d)
  x <- d[[xcol]]
  out <- d[rep(seq_len(n), each = 2L), , drop = FALSE]
  # Row pair i holds y_i from x_i to x_{i+1}. The construction leaves a duplicate of the
  # final vertex, which is dropped so the result is the 2n-1 vertices ggplot2's own
  # stairstep() produces rather than 2n with a degenerate last segment.
  out[[xcol]] <- as.vector(rbind(x, c(x[-1], x[n])))
  out[-nrow(out), , drop = FALSE]
}

#' Format a p-value for a panel annotation without ever printing "p = 0".
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
fmt_ci <- function(e, l, h, d = 2) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), e, l, h)

# -----------------------------------------------------------------------------
# 3. DATA ASSEMBLY
# -----------------------------------------------------------------------------
message("\n", strrep("=", 74))
message("[ascert] 45_ascertainment_figures.R — detection vs onset, spatial ascertainment")
message(strrep("=", 74))

# ---- 3.1 the modelling spine, the invasion labels, and the vulnerability pillars ----
rs_all <- readr::read_csv(fs_risk_csv(OUT), show_col_types = FALSE)
stopifnot(all(c("health_zone", "province", "horizon", "was_active_before",
                "p_case_invasion", "surveillance_gap", "healthcare_gap",
                "social_vulnerability", "access_gap", "healthcare_travel_min") %in% names(rs_all)))
RISK_METHOD <- unique(rs_all$method)
stopifnot(length(RISK_METHOD) == 1L)
PROB_SCALE  <- unique(rs_all$prob_scale)
stopifnot(length(PROB_SCALE) == 1L)   # it is printed on a panel; two scales in one table
                                      # would be a caption naming the wrong one

spine <- rs_all %>%
  dplyr::filter(horizon == 1L) %>%
  dplyr::transmute(health_zone, province,
                   invaded = as.logical(was_active_before),
                   surveillance_gap, healthcare_gap, social_vulnerability,
                   access_gap, healthcare_travel_min)
stopifnot(nrow(spine) == dplyr::n_distinct(spine$health_zone), !anyNA(spine$invaded))
message(sprintf("[ascert] spine: %d zones, %d invaded, %d never invaded (risk method %s, scale %s)",
                nrow(spine), sum(spine$invaded), sum(!spine$invaded), RISK_METHOD, PROB_SCALE))

# The two access pillars are DEGENERATE on this covariate set: every one of the 519
# zones hosts at least one health facility, so the "own facility -> 0 minutes" rule in
# compute_vulnerability_index() zeroes travel time and hence access_gap for all of them
# (20_forecast_detail.R documents this and drops the pillar from V for the same reason).
# A constant cannot be an ascertainment axis, so it is dropped here too — but by TEST,
# not by assumption, so a future run in which the pillar becomes informative keeps it.
.degenerate <- function(v) dplyr::n_distinct(v[is.finite(v)]) <= 1L
DROPPED_PILLARS <- c("access_gap", "healthcare_travel_min")[
  c(.degenerate(spine$access_gap), .degenerate(spine$healthcare_travel_min))]
if (length(DROPPED_PILLARS))
  message("[ascert] degenerate pillar(s) dropped from the proxy set: ",
          paste(DROPPED_PILLARS, collapse = ", "))

# ---- 3.2 raw covariate layers (facility density, population, deprivation) -----------
cov <- load_static_covariates() %>%
  dplyr::transmute(health_zone = nom,
                   pop_count = as.numeric(pop_count),
                   healthsite_count = as.numeric(healthsite_count),
                   healthsite_density = as.numeric(healthsite_density),
                   ccvi = as.numeric(ccvi))
stopifnot(setequal(cov$health_zone, spine$health_zone))

# ---- 3.3 travel time from the epicentre region, for ALL zones -----------------------
# arrival_predictors.csv carries this only for the INVADED zones; the invaded-vs-not
# comparison needs it for the whole spine. It is recomputed here from the same OSRM matrix
# and the same min-over-EPICENTRE_ZONES rule that 43_spread_kinematics.R uses, then
# ASSERTED equal to the published column on every zone where both exist — so the new
# values come from a rule checked against the published one, not assumed to match it.
geo <- build_geo()
time_h <- geo$time                                   # hours, symmetrised
epi_set <- intersect(EPICENTRE_ZONES, rownames(time_h))
stopifnot(length(epi_set) == length(EPICENTRE_ZONES))
.epi_h <- apply(time_h[epi_set, , drop = FALSE], 2L,
                function(z) { z <- z[is.finite(z)]; if (length(z)) min(z) else NA_real_ })
epi_travel <- tibble::tibble(health_zone = names(.epi_h), epi_travel_h = as.numeric(.epi_h))
stopifnot(setequal(epi_travel$health_zone, spine$health_zone))

arr <- readr::read_csv(file.path(OUT, "key_outputs", "arrival_predictors.csv"),
                       show_col_types = FALSE) %>%
  dplyr::select(health_zone, first_onset, arrival_days, arrival_origin, in_fit_window,
                travel_time_h, dist_greatcircle_km, neglog_share_shorttrip,
                log10_population, cases_total)
# 43's origin of arrival_days (one date), named on the arrival axes rather than described.
ARRIVAL_ORIGIN <- unique(stats::na.omit(as.Date(arr$arrival_origin)))
if (length(ARRIVAL_ORIGIN) != 1L)
  stop("[ascert] arrival_predictors.csv carries no single arrival_origin; re-run 43.", call. = FALSE)
.arr_axis <- sprintf("Arrival time (days from %s)", format(ARRIVAL_ORIGIN, "%d %b %Y"))
.chk <- dplyr::inner_join(arr, epi_travel, by = "health_zone")
stopifnot(nrow(.chk) == nrow(arr))
.dev <- max(abs(.chk$travel_time_h - .chk$epi_travel_h), na.rm = TRUE)
if (!is.finite(.dev) || .dev > 1e-8)
  stop(sprintf(paste0("[ascert] recomputed epicentre travel time disagrees with ",
                      "arrival_predictors.csv by up to %.6g h — the two are NOT the same ",
                      "covariate and nothing downstream may treat them as one."), .dev))
message(sprintf("[ascert] epicentre travel time reproduced for all %d zones (max deviation %.2g h on the %d published)",
                nrow(epi_travel), .dev, nrow(arr)))

# ---- 3.4 line-list ascertainment measures -------------------------------------------
ll <- load_linelist()
# THE CONFIRMING TEST'S DATES. load_linelist() (apply_earliest_sample_dates(), 00_config.R) gives
# every record its earliest diagnostic sample date and every currently-confirmed record the lab
# date of its CONFIRMING test (the earliest lab_analysis_date recorded while it was classified
# confirmed), keeping the extract's values in *_extract columns. The extract's own lab date
# belongs to the record's LATEST (repeat or discharge) test, so without this step the
# onset -> confirmation delay below (3.4) and the as-known province share (section 4) would be
# dated by the wrong test. The step is skipped, with only a warning, when the snapshot archive is
# missing or unreadable, so its presence is asserted rather than assumed.
if (!"lab_analysis_date_extract" %in% names(ll))
  stop("[ascert] load_linelist() did not restore the confirming test's lab date (no ",
       "lab_analysis_date_extract column: the snapshot archive was missing or unreadable). The ",
       "delays and as-known shares here would be dated by each record's latest test.", call. = FALSE)
# The sitrep top-up rows are a COUNT reconciliation against the official cumulative, not
# laboratory records: they carry a classification but no test, no lab date and no alert.
# Including them would inflate the confirmation ratio in exactly the zones that received a
# top-up — i.e. it would manufacture the spatial pattern being measured.
n_synth <- sum(grepl("^SITREP-CONF-", ll$alert_id))
ll_real <- ll %>% dplyr::filter(is.na(alert_id) | !grepl("^SITREP-CONF-", alert_id))
stopifnot(nrow(ll) - nrow(ll_real) == n_synth)
message(sprintf("[ascert] line list: %d rows, %d synthetic sitrep top-up rows excluded, %d retained",
                nrow(ll), n_synth, nrow(ll_real)))

# Alerts are counted on the SAME window the zone-week grid uses: on or after
# OUTBREAK_START (load_linelist already enforces this) and on or before ANALYSIS_DATE.
ll_win <- ll_real %>%
  dplyr::filter(!is.na(health_zone), !is.na(date_index),
                date_index >= OUTBREAK_START, date_index <= ANALYSIS_DATE)
message(sprintf("[ascert] dated records inside [%s, %s]: %d (dropped %d outside the window or without a date_index)",
                format(OUTBREAK_START), format(ANALYSIS_DATE), nrow(ll_win), nrow(ll_real) - nrow(ll_win)))

# Sitrep top-up rows per zone, on the SAME window as ll_win. They are excluded from every
# alert count, so a zone invaded ONLY through them is invaded with no laboratory-resolved
# alert. This count is the evidence for that statement wherever the report or Figure A5C
# makes it: the claim is tested zone by zone, never inferred from the absence of a
# line-list confirmed case (which would be true of such a zone by construction).
zone_sitrep <- ll %>%
  dplyr::filter(!is.na(alert_id), grepl("^SITREP-CONF-", alert_id), !is.na(health_zone),
                !is.na(date_index), date_index >= OUTBREAK_START, date_index <= ANALYSIS_DATE) %>%
  dplyr::count(health_zone, name = "n_sitrep_conf")

# Counted by the SHARED definition (01_data_prep.R), which the manuscript watch-list figure
# (Figure3_rank_surveillance, panel D) also uses: a record is an alert in the window when its
# case date (date_index) OR its notification date falls in it, and that includes the records
# load_linelist() leaves out for lacking an onset or specimen date on or after OUTBREAK_START
# (load_unplaced_alerts(), placed by notification). Those are all unclassified, so they add to
# alerts_total and to nothing else; before 2026-10-03 they were missing from it. Every record
# of ll_win is in the count (the rule can only add line-list rows, not remove them); ll_win
# itself still feeds the delay and effort measures below.
unplaced_alerts <- load_unplaced_alerts()
zone_alerts <- zone_alert_counts(ll, unplaced_alerts, from = OUTBREAK_START, to = ANALYSIS_DATE)
n_alerts_ll       <- attr(zone_alerts, "n_line_list")
n_alerts_unplaced <- attr(zone_alerts, "n_unplaced")
n_alerts_window   <- sum(zone_alerts$alerts_total)
stopifnot(n_alerts_ll >= nrow(ll_win), n_alerts_window == n_alerts_ll + n_alerts_unplaced)
message(sprintf("[ascert] alerts in the window: %d (%d line-list records + %d with no onset or specimen date on or after %s, placed by notification)",
                n_alerts_window, n_alerts_ll, n_alerts_unplaced, format(OUTBREAK_START)))

# ALERTS WHOSE ZONE IS NOT IN THE SPINE. zone_alert_counts() keys on the line list's zone
# name, and a name that maps to no spine zone would be dropped by the left-join into
# zone_ascert below. (The line-list reader behind both load_linelist() and
# load_unplaced_alerts() resolves the spine's same-named zones, such as the two "Bili", by
# province, so this is a guard, not an expected path.) Such alerts would count towards NO zone,
# province or figure here, so they are counted, warned about, and left out of the report totals,
# which are taken from the mapped (joined) zones only. The per-source split is recomputed
# through the same shared definition, so the report's line-list / unplaced breakdown sums to the
# mapped total.
unmapped_alert_zones <- sort(setdiff(zone_alerts$health_zone, spine$health_zone))
n_alerts_unmapped <- sum(zone_alerts$alerts_total[zone_alerts$health_zone %in% unmapped_alert_zones])
.mapped_alerts <- function(a, u) {
  z <- zone_alert_counts(a, u, from = OUTBREAK_START, to = ANALYSIS_DATE)
  sum(z$alerts_total[z$health_zone %in% spine$health_zone])
}
n_alerts_ll_mapped       <- .mapped_alerts(ll, unplaced_alerts[0, , drop = FALSE])
n_alerts_unplaced_mapped <- .mapped_alerts(ll[0, , drop = FALSE], unplaced_alerts)
n_alerts_window_mapped   <- n_alerts_ll_mapped + n_alerts_unplaced_mapped
stopifnot(n_alerts_window_mapped + n_alerts_unmapped == n_alerts_window)
if (length(unmapped_alert_zones))
  warning(sprintf(paste0("[ascert] %d of %d alert(s) in the window carry a zone name that is not in ",
                         "the %d-zone spine (%d name(s): %s). They are in no zone-level count, ",
                         "panel or report total."),
                  n_alerts_unmapped, n_alerts_window, nrow(spine), length(unmapped_alert_zones),
                  paste(unmapped_alert_zones, collapse = ", ")), call. = FALSE)

# The confirmation ratio is an ALERT-level positivity: the share of laboratory-RESOLVED
# alerts that were confirmed. It is not a test-level positivity — a person may be swabbed
# more than once and DHIS2 carries no per-record test count (samples_analyzed is NA on
# every real record, which is why the pipeline's own `positivity` covariate is NA
# throughout). Alerts still under investigation (suspected/probable) and alerts never
# classified are in neither the numerator nor the denominator.
.w <- .wilson_ci_tbl(zone_alerts$n_confirmed, zone_alerts$alerts_resolved)
zone_alerts <- zone_alerts %>%
  dplyr::mutate(confirm_ratio = .w$est, confirm_lo = .w$lo, confirm_hi = .w$hi,
                resolution_frac = dplyr::if_else(alerts_total > 0,
                                                 alerts_resolved / alerts_total, NA_real_))

# Onset -> laboratory confirmation delay, per zone. Three restrictions, each needed:
#   (a) NATIVE onsets only. load_linelist() imputes a large minority of onsets backwards
#       from the sample date by an interval drawn from same-class, same-week donors'
#       recorded sample-onset intervals, so a delay computed on those records is a function
#       of the donors' intervals and would be circular. The `onset_imputed` flag
#       identifies them.
#   (b) EQUAL FOLLOW-UP. Only cases whose onset is at least DELAY_CAP_DAYS before the
#       analysis date are used, so every zone's cases have had the same opportunity to be
#       confirmed; without this, late-invaded zones look faster purely by truncation.
#   (c) Delays are capped at DELAY_CAP_DAYS, matching week_completeness(). A zone in which
#       more than half of the eligible cases exceed the cap has a censored median, flagged
#       rather than silently reported as exactly the cap.
# The confirmation date is `lab_analysis_date`, which load_linelist() has set to the CONFIRMING
# test's lab date (asserted above); the extract's value, the latest test's, would add the
# interval to a repeat or discharge test to every re-tested case. Reporting and notification
# dates are fall-backs for a confirmed record with no lab date at all.
.conf_date <- dplyr::coalesce(as.Date(ll_win$lab_analysis_date),
                              as.Date(ll_win$reporting_date),
                              as.Date(ll_win$date_of_notification))
.delay_d <- as.numeric(.conf_date - ll_win$date_index)
.delay_ok <- ll_win$confirmed %in% TRUE &
             !(ll_win$onset_imputed %in% TRUE) &
             is.finite(.delay_d) & .delay_d >= 0 &
             ll_win$date_index <= (ANALYSIS_DATE - DELAY_CAP_DAYS)
zone_delay <- tibble::tibble(health_zone = ll_win$health_zone,
                             delay = pmin(.delay_d, DELAY_CAP_DAYS))[.delay_ok, ] %>%
  dplyr::group_by(health_zone) %>%
  dplyr::summarise(n_delay = dplyr::n(),
                   delay_median = stats::median(delay),
                   delay_p75 = stats::quantile(delay, 0.75, names = FALSE),
                   delay_over_cap = mean(delay >= DELAY_CAP_DAYS),
                   .groups = "drop") %>%
  dplyr::mutate(delay_censored = delay_over_cap >= 0.5) %>%
  dplyr::filter(n_delay >= 5L)
message(sprintf("[ascert] onset->confirmation delay estimable in %d zone(s) (>=5 native-onset confirmed cases with >=%d d follow-up); %d have a censored median",
                nrow(zone_delay), DELAY_CAP_DAYS, sum(zone_delay$delay_censored)))

# (The line-list / sitrep COVERAGE RATIO proxy that used to be read here, from
# key_outputs/coverage_ratio_zone.csv, is gone: the careseeking module that wrote that file was
# deleted on 2026-09-21 and nothing in the suite writes it any more, so every run published a
# frozen copy of an older line list beside measures recomputed from the current one.)

# ---- 3.5 the assembled zone table ----------------------------------------------------
zone_ascert <- spine %>%
  dplyr::left_join(cov,         by = "health_zone") %>%
  dplyr::left_join(epi_travel,  by = "health_zone") %>%
  dplyr::left_join(zone_alerts, by = "health_zone") %>%
  dplyr::left_join(zone_delay,  by = "health_zone") %>%
  dplyr::left_join(arr,         by = "health_zone") %>%
  dplyr::mutate(
    # Zones with no alert at all have no alert counts to be missing — they are zero, and
    # a zero alert count is an observation about surveillance, not a gap in the data.
    # A confirmation RATIO, by contrast, stays NA there: 0/0 is undefined, not 0.
    alerts_total     = dplyr::coalesce(alerts_total, 0L),
    n_confirmed      = dplyr::coalesce(n_confirmed, 0L),
    n_not_a_case     = dplyr::coalesce(n_not_a_case, 0L),
    alerts_resolved  = dplyr::coalesce(alerts_resolved, 0L),
    resolved_per_100k = dplyr::if_else(is.finite(pop_count) & pop_count > 0,
                                       alerts_resolved / pop_count * 1e5, NA_real_),
    log_pop          = log10(pmax(pop_count, 1)),
    status           = factor(dplyr::if_else(invaded, "Invaded", "Never invaded"),
                              levels = c("Never invaded", "Invaded")))
# STABLE ORDER, set once and inherited by every artifact derived from this table. Without
# it the row order came from whatever order `bayes_risk_scores_all_zones.csv` happened to
# be written in, so a re-run of the pipeline reshuffled every published CSV and made the
# two versions impossible to diff even where no value had changed.
zone_ascert <- dplyr::arrange(zone_ascert, health_zone)
stopifnot(nrow(zone_ascert) == nrow(spine),
          dplyr::n_distinct(zone_ascert$health_zone) == nrow(spine),
          !is.unsorted(zone_ascert$health_zone),
          # the report totals are taken from this table: it must hold exactly the mapped alerts
          sum(zone_ascert$alerts_total) == n_alerts_window_mapped)

# Invaded zones with NO laboratory-resolved alert, checked against the sitrep top-up rows
# rather than assumed to be explained by them (report section 2).
inv_no_lab <- zone_ascert %>%
  dplyr::filter(invaded, alerts_resolved == 0L) %>%
  dplyr::select(health_zone) %>%
  dplyr::left_join(zone_sitrep, by = "health_zone") %>%
  dplyr::mutate(n_sitrep_conf = dplyr::coalesce(n_sitrep_conf, 0L))
INV_NO_LAB_SITREP_ONLY <- all(inv_no_lab$n_sitrep_conf > 0L)
message(sprintf("[ascert] invaded zones with no laboratory-resolved alert: %d (%d of them with a sitrep top-up row in the window)",
                nrow(inv_no_lab), sum(inv_no_lab$n_sitrep_conf > 0L)))
if (!INV_NO_LAB_SITREP_ONLY)
  warning(sprintf(paste0("[ascert] %d invaded zone(s) have neither a laboratory-resolved alert nor a ",
                         "sitrep top-up row in the window: %s"),
                  sum(inv_no_lab$n_sitrep_conf == 0L),
                  paste(inv_no_lab$health_zone[inv_no_lab$n_sitrep_conf == 0L], collapse = ", ")),
          call. = FALSE)

# ---- 3.6 SURVEILLANCE EFFORT, measured so it is not circular ------------------------
# "Did this zone ever have a laboratory-resolved alert?" is NOT a usable covariate for a
# comparison whose outcome is "did this zone record a confirmed case?": a line-list
# confirmed case IS a resolved alert, so every zone invaded through the line list has one.
# (Not every invaded zone: a zone confirmed only through sitrep top-up rows, which carry no
# test and are excluded from the alert counts, is invaded with no resolved alert —
# `inv_no_lab` above counts and checks them.) The link is mechanical either way. The DISCARDED alerts
# (`not_a_case`) carry no such mechanical link — they are investigations that found
# nothing — so they are what the surveillance-effort measures below count. Where the
# quantity is needed AT A FORECAST ORIGIN it is recomputed from the alerts dated on or
# before that origin, never from the end-of-outbreak total.
# Each discarded alert is dated by its NOTIFICATION date — when it entered the system and
# an investigation began — falling back to the case's own index date if that is missing.
# The laboratory result came later, so this is "an alert was being investigated by then",
# not "an alert had been resolved by then"; the panels that use it say so.
# An IMPLAUSIBLE date counts as missing: one after the analysis date (impossible in this
# extract) or before OUTBREAK_START less 30 days, the floor linelist_observation_date() and the
# earliest-sample rule (00_config.R) apply. A mistyped year would otherwise date an alert to
# before the outbreak (counting it as "already notified" at every origin) or after the last one.
.plaus_date <- function(v) {
  v <- .parse_date(v)
  v[!is.na(v) & (v < as.Date(OUTBREAK_START) - 30L | v > as.Date(ANALYSIS_DATE))] <- NA
  v
}
neg_alerts_long <- ll_win %>%
  dplyr::filter(final_mve_case_classification %in% NOT_A_CASE_STATUS) %>%
  dplyr::transmute(health_zone,
                   notif_implausible = !is.na(.parse_date(date_of_notification)) &
                                       is.na(.plaus_date(date_of_notification)),
                   alert_date = dplyr::coalesce(.plaus_date(date_of_notification), date_index))
stopifnot(!anyNA(neg_alerts_long$alert_date))
message(sprintf("[ascert] discarded (not-a-case) alerts available for the effort measure: %d across %d zone(s); %d dated by their index date because the notification date is implausible",
                nrow(neg_alerts_long), dplyr::n_distinct(neg_alerts_long$health_zone),
                sum(neg_alerts_long$notif_implausible)))

zone_ascert <- zone_ascert %>%
  dplyr::mutate(neg_alerts = n_not_a_case,
                neg_per_100k = dplyr::if_else(is.finite(pop_count) & pop_count > 0,
                                              n_not_a_case / pop_count * 1e5, NA_real_),
                ever_investigated = alerts_resolved > 0L)
stopifnot(nrow(inv_no_lab) == sum(!zone_ascert$ever_investigated & zone_ascert$invaded))

.n_na_travel <- sum(!is.finite(zone_ascert$epi_travel_h))
if (.n_na_travel > 0L)
  message(sprintf("[ascert] %d zone(s) have no road route to the epicentre region and are absent from every travel-time-conditioned panel: %s",
                  .n_na_travel, paste(zone_ascert$health_zone[!is.finite(zone_ascert$epi_travel_h)], collapse = ", ")))

# The proxy panel, split by whether a quantity can legitimately enter an INVADED-vs-NOT
# comparison. The exogenous set is fixed by the health system and geography and is
# measurable in a zone that never reported anything; the outbreak-conditioned set is
# partly generated BY the outbreak and enters only Figure A2, where every zone compared has
# already been invaded. (Figure A5's outbreak-generated quantities, the discarded alerts and
# the province confirmed share, are re-measured at each forecast origin in section 4.)
PROXY_EXO <- tibble::tribble(
  ~var,                   ~label,                                          ~worse_high,
  "surveillance_gap",     "Surveillance gap\n(facility density pctile)",    TRUE,
  "healthcare_gap",       "Health-system gap\n(facilities per capita pctile)", TRUE,
  "social_vulnerability", "Deprivation\n(CCVI percentile)",                 TRUE,
  "log_hs_density",       "Health-facility density\n(log10 per km2)",       FALSE,
  "log_pop",              "Population\n(log10)",                            FALSE)
# The access pillars join the exogenous set IF the degeneracy test above clears them, so
# the test governs the analysis rather than only a log line. On the current covariate set
# it does not clear them and they are absent from every panel.
.access_kept <- setdiff(c("access_gap", "healthcare_travel_min"), DROPPED_PILLARS)
if (length(.access_kept))
  PROXY_EXO <- dplyr::bind_rows(PROXY_EXO, tibble::tibble(
    var = .access_kept,
    label = c(access_gap = "Healthcare-access gap\n(travel-time percentile)",
              healthcare_travel_min = "Travel time to care\n(minutes)")[.access_kept],
    worse_high = TRUE))
PROXY_ENDO <- tibble::tribble(
  ~var,              ~label,                                        ~worse_high,
  "confirm_ratio",   "Confirmation ratio\n(confirmed / resolved alerts)", TRUE,
  "log_neg_100k",    "Discarded alerts per 100k\n(log10)",           FALSE,
  "delay_median",    "Onset to confirmation\n(median days)",         TRUE)

zone_ascert <- zone_ascert %>%
  dplyr::mutate(
    log_hs_density = dplyr::if_else(is.finite(healthsite_density) & healthsite_density > 0,
                                    log10(healthsite_density), NA_real_),
    # log10 of a rate that is legitimately zero: the zero is kept as its own state by the
    # `ever_investigated` flag and the rate is NA here, rather than being shifted by an
    # arbitrary pseudo-count that would place "no surveillance at all" on the same
    # continuum as "a little surveillance".
    log_neg_100k = dplyr::if_else(is.finite(neg_per_100k) & neg_per_100k > 0,
                                  log10(neg_per_100k), NA_real_))
PROXY_ALL <- dplyr::bind_rows(dplyr::mutate(PROXY_EXO, kind = "exogenous"),
                              dplyr::mutate(PROXY_ENDO, kind = "outbreak-conditioned"))
# Labels are the row key of every forest below, so a duplicate would silently merge two
# proxies into one axis position; vars must exist or a panel would plot NA.
stopifnot(all(PROXY_ALL$var %in% names(zone_ascert)),
          anyDuplicated(PROXY_ALL$var) == 0L,
          anyDuplicated(gsub("\n", " ", PROXY_ALL$label)) == 0L)

# -----------------------------------------------------------------------------
# 4. FORECAST FRAMES
# -----------------------------------------------------------------------------
# The suite's ONE reader of the featured model (figure_style.R, sourced above): `bayesian`
# first, `headline` as the fallback, and an error rather than a guess. This module used to carry
# its own copy of that rule.
FEATURED <- read_featured_model(file.path(OUT, "key_outputs", "model_selection.json"))
message("[ascert] featured model (read from model_selection.json): ", FEATURED)

lfo_all <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))
stopifnot(FEATURED %in% lfo_all$method)
# SUBSET TO THE FEATURED MODEL BEFORE APPLYING THE SCALE. fs_lfo_col() inspects every
# scored row in the frame it is handed and falls back to the raw column if ANY of them
# lacks a recalibrated probability. Two rank-only baselines (Adjacency-B7, Distance-B1)
# deliberately carry no `p_recal` — recalibrating an ordering is meaningless — so handing
# it the whole multi-method frame silently demotes the featured model to the raw scale too.
# The featured model's `p_recal` is complete, so filtering first is what actually puts
# this module on the scale its outputs claim. (Every RANK in this module is taken within
# a fold and is therefore invariant to the choice either way, and so is the within-round
# discrimination of Figure A4C, which strat_discrimination() scores on the RAW column
# `p_invasion_raw`, kept by fs_apply_lfo_scale(), as the headline evaluation does.)
lfo <- dplyr::filter(lfo_all, method == FEATURED)
PCOL <- fs_apply_lfo_scale(lfo)     # puts the selected scale into `p_invasion`
message("[ascert] LFO probability column in use: ", PCOL)

# At-risk rows of the featured model ON THE SHARED (fold x zone) SUPPORT, exactly the rows
# evaluate_invasion() scores: not already affected, finite probability, and a cell that EVERY
# scored method covers. The support is resolved by invasion_common_cells() (16_invasion_eval.R)
# on the whole multi-method frame -- the same code evaluate_invasion() and
# compute_detection_curve() use, fallback included -- and applied BEFORE ranking, so a round's
# watch-list is the one the evaluation ranks. (A horizon with no common cell is kept whole, as
# evaluate_invasion() keeps it.) Ranks are taken WITHIN (horizon, fold) — the whole watch-list
# of that round — with the pipeline's two tie conventions:
#   rank_avg  ties averaged, matching .ranking_metrics()'s mean rank of truth;
#   rank_max  ties to the worst position, matching compute_detection_curve()'s top-K,
#             so a zone counts as monitored only if monitoring K zones must include it.
# Both are monotone in p within a fold, so they are invariant to the probability scale.
lfo_at_risk <- lfo %>%
  dplyr::filter(!(as.logical(was_active_before) %in% TRUE), is.finite(p_invasion))
.cells_by_h <- lapply(stats::setNames(nm = sort(unique(lfo_at_risk$horizon))),
                      function(h) invasion_common_cells(lfo_all, h))
.on_support <- unlist(lapply(names(.cells_by_h), function(h) {
  i <- which(lfo_at_risk$horizon == as.integer(h)); cl <- .cells_by_h[[h]]
  if (is.null(cl)) i else i[paste(lfo_at_risk$fold_id[i], lfo_at_risk$health_zone[i], sep = "\r") %in% cl]
}), use.names = FALSE)
message(sprintf("[ascert] shared (fold x zone) support: %s featured at-risk row(s) outside it dropped",
                paste(sprintf("h=%s %d", names(.cells_by_h),
                              vapply(names(.cells_by_h), function(h)
                                sum(lfo_at_risk$horizon == as.integer(h)) -
                                  sum(lfo_at_risk$horizon[.on_support] == as.integer(h)), integer(1))),
                      collapse = ", ")))
lfo_f <- lfo_at_risk[sort(.on_support), , drop = FALSE] %>%
  dplyr::group_by(horizon, fold_id) %>%
  dplyr::mutate(rank_avg = rank(-p_invasion, ties.method = "average"),
                rank_max = rank(-p_invasion, ties.method = "max"),
                n_atrisk_fold = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(zone_ascert, health_zone, province, surveillance_gap,
                                 healthcare_gap, social_vulnerability, log_hs_density,
                                 epi_travel_h, ever_investigated, neg_alerts),
                   by = "health_zone")
stopifnot(nrow(lfo_f) > 0, !anyNA(lfo_f$surveillance_gap), !anyNA(lfo_f$rank_avg))

# Prior surveillance effort AT EACH FORECAST ORIGIN: discarded alerts notified on or
# before that fold's forecast origin, lfo_origin(cutoff) = cutoff + LFO_ORIGIN_OFFSET_DAYS
# (00_config.R). `cutoff` is the START of the last training week; the forecast is issued
# when that week closes, so the origin, not the cutoff, is what "already notified at the
# forecast origin" (Figure A5A, report section 6) refers to. Built once per distinct cutoff
# (one per back-test round) so that "before the origin" has exactly ONE definition across
# every panel that uses it. The join back to lfo_f stays keyed on `cutoff`.
.cut <- sort(unique(lfo_f$cutoff))
.neg_by_cut <- purrr::map_dfr(.cut, function(cu) {
  neg_alerts_long %>%
    dplyr::filter(alert_date <= lfo_origin(cu)) %>%
    dplyr::count(health_zone, name = "neg_prior") %>%
    dplyr::mutate(cutoff = cu)
})
lfo_f <- lfo_f %>%
  dplyr::left_join(.neg_by_cut, by = c("health_zone", "cutoff")) %>%
  dplyr::mutate(neg_prior = dplyr::coalesce(neg_prior, 0L),
                watched_prior = neg_prior > 0L)

# THE PROVINCE CONFIRMED SHARE AS KNOWN AT EACH FORECAST ORIGIN (Figure A5B/C). The
# end-of-outbreak share counts confirmations made after the origin -- including those of the very
# zones a forecast ranks -- so a correlation between it and the ranking is partly mechanical.
# Re-measured here instead: at the origin of round `cutoff`, a province's share is its confirmed
# alerts over its laboratory-resolved alerts (confirmed + not-a-case) among those both NOTIFIED
# and RESOLVED by that origin. A record is resolved on its lab date -- for a confirmed record the
# confirming test's (load_linelist(), asserted in 3.4) -- or, when that is implausible
# (.plaus_date()), on its notification date; it is known on the later of the two. The alerts are
# exactly those zone_alert_counts() counts (same window, sitrep top-ups excluded, mapped zones
# only), which is checked below: re-measured at the analysis date, every province's counts must
# equal Figure A1C's. Residual, stated: the classification is the final one, so a record
# reclassified after the origin is counted in its final class from its own date.
.in_win <- function(v) { v <- .parse_date(v); !is.na(v) & v >= OUTBREAK_START & v <= ANALYSIS_DATE }
res_long <- dplyr::bind_rows(
  ll_real %>%
    dplyr::filter(!is.na(health_zone), .in_win(date_index) | .in_win(date_of_notification)) %>%
    dplyr::transmute(health_zone, cls = final_mve_case_classification,
                     lab = .plaus_date(lab_analysis_date), notif = .plaus_date(date_of_notification)),
  unplaced_alerts %>%
    dplyr::filter(is.na(alert_id) | !grepl("^SITREP-CONF-", alert_id), !is.na(health_zone),
                  .in_win(date_of_notification)) %>%
    dplyr::transmute(health_zone, cls = final_mve_case_classification,
                     lab = .plaus_date(if ("lab_analysis_date" %in% names(unplaced_alerts))
                                         lab_analysis_date else NA_character_),
                     notif = .plaus_date(date_of_notification))) %>%
  dplyr::filter(cls %in% c(CONFIRMED_STATUS, NOT_A_CASE_STATUS)) %>%
  dplyr::inner_join(dplyr::select(zone_ascert, health_zone, province), by = "health_zone") %>%
  dplyr::mutate(confirmed = cls %in% CONFIRMED_STATUS,
                resolved_on = dplyr::coalesce(lab, notif),
                known = pmax(resolved_on, notif, na.rm = TRUE))
.res_undated <- res_long %>% dplyr::filter(is.na(known)) %>%
  dplyr::group_by(province) %>%
  dplyr::summarise(n_undated = dplyr::n(), n_undated_conf = sum(confirmed), .groups = "drop")
message(sprintf(paste0("[ascert] as-known province share: %d resolved alerts; %d with no plausible lab ",
                       "date, resolved on their notification date; %d with neither, in no as-known count"),
                nrow(res_long), sum(is.na(res_long$lab) & !is.na(res_long$notif)),
                sum(.res_undated$n_undated)))
.known_by <- function(o) res_long %>%
  dplyr::filter(!is.na(known), known <= o) %>%
  dplyr::group_by(province) %>%
  dplyr::summarise(confirmed_o = sum(confirmed), resolved_o = dplyr::n(), .groups = "drop")
prov_origin <- purrr::map_dfr(.cut, function(cu) dplyr::mutate(.known_by(lfo_origin(cu)), cutoff = cu)) %>%
  dplyr::mutate(origin = lfo_origin(cutoff), share_o = confirmed_o / resolved_o) %>%
  dplyr::relocate(cutoff, origin)
.recon <- zone_ascert %>%
  dplyr::group_by(province) %>%
  dplyr::summarise(confirmed = sum(n_confirmed), resolved = sum(alerts_resolved), .groups = "drop") %>%
  dplyr::left_join(.known_by(ANALYSIS_DATE), by = "province") %>%
  dplyr::left_join(.res_undated, by = "province") %>%
  dplyr::mutate(dplyr::across(c(confirmed_o, resolved_o, n_undated, n_undated_conf),
                              ~ dplyr::coalesce(.x, 0L)))
stopifnot(all(.recon$resolved  == .recon$resolved_o  + .recon$n_undated),
          all(.recon$confirmed == .recon$confirmed_o + .recon$n_undated_conf))
lfo_f <- lfo_f %>%
  dplyr::left_join(dplyr::select(prov_origin, cutoff, province, prov_confirmed_o = confirmed_o,
                                 prov_resolved_o = resolved_o),
                   by = c("cutoff", "province")) %>%
  dplyr::mutate(prov_confirmed_o = dplyr::coalesce(prov_confirmed_o, 0L),
                prov_resolved_o  = dplyr::coalesce(prov_resolved_o, 0L),
                prov_share_o     = dplyr::if_else(prov_resolved_o > 0L,
                                                  prov_confirmed_o / prov_resolved_o, NA_real_))

# THE MAIN HORIZON of every single-horizon panel (A4A, A5) and of the counts quoted from them.
# lfo_f holds one row per (horizon, fold, zone), so a count over all of it counts each
# zone-forecast once PER HORIZON; the log and the report count within a horizon.
H_MAIN <- 2L
stopifnot(H_MAIN %in% lfo_f$horizon)
.n_by_h <- lfo_f %>% dplyr::count(horizon)
message(sprintf("[ascert] LFO at-risk zone-forecasts: %s (%d zones, %d rounds); at h=%d, %d of %d (%.1f%%) had a discarded alert notified by the origin",
                paste(sprintf("h=%d %d", .n_by_h$horizon, .n_by_h$n), collapse = ", "),
                dplyr::n_distinct(lfo_f$health_zone), dplyr::n_distinct(lfo_f$fold_id), H_MAIN,
                sum(lfo_f$watched_prior[lfo_f$horizon == H_MAIN]), sum(lfo_f$horizon == H_MAIN),
                100 * mean(lfo_f$watched_prior[lfo_f$horizon == H_MAIN])))

# Surveillance stratum: a zone-level split at the MEDIAN surveillance gap of the zones the
# LFO actually scores, so a zone belongs to one stratum in every fold and the zone-cluster
# bootstrap in strat_discrimination() resamples whole zones within a stratum.
.sg_zone <- lfo_f %>% dplyr::distinct(health_zone, surveillance_gap)
SG_SPLIT <- stats::median(.sg_zone$surveillance_gap)
lfo_f <- lfo_f %>%
  dplyr::mutate(sg_stratum = factor(dplyr::if_else(surveillance_gap <= SG_SPLIT,
                                                   "Denser facilities", "Sparser facilities"),
                                    levels = c("Denser facilities", "Sparser facilities")))
message(sprintf("[ascert] surveillance stratum split at gap = %.4f (%d vs %d zones)",
                SG_SPLIT, sum(.sg_zone$surveillance_gap <= SG_SPLIT),
                sum(.sg_zone$surveillance_gap > SG_SPLIT)))

# Live (current) forecast, at-risk zones only, for the operational panels.
live <- rs_all %>%
  dplyr::filter(!(as.logical(was_active_before) %in% TRUE), is.finite(p_case_invasion)) %>%
  dplyr::select(health_zone, province, horizon, p_case_invasion) %>%
  dplyr::left_join(dplyr::select(zone_ascert, health_zone, surveillance_gap,
                                 ever_investigated, neg_alerts, epi_travel_h),
                   by = "health_zone")

# -----------------------------------------------------------------------------
# 5. FIGURE A1 — the ascertainment landscape
# -----------------------------------------------------------------------------
message("\n[ascert] Figure A1 — ascertainment landscape")

shp <- geo$shp
stopifnot(setequal(as.character(shp$Nom), zone_ascert$health_zone))
shp_a <- shp %>%
  dplyr::mutate(health_zone = as.character(Nom)) %>%
  dplyr::left_join(zone_ascert, by = "health_zone")
prov_outline <- shp_a %>% dplyr::group_by(.prov = as.character(PROVINCE)) %>%
  dplyr::summarise(.groups = "drop") %>% sf::st_geometry()
inv_outline <- shp_a %>% dplyr::filter(invaded) %>% sf::st_geometry()

pA1a <- ggplot() +
  geom_sf(data = shp_a, aes(fill = surveillance_gap), colour = "white", linewidth = 0.06) +
  geom_sf(data = prov_outline, fill = NA, colour = "grey45", linewidth = 0.18) +
  geom_sf(data = inv_outline, fill = NA, colour = "white", linewidth = 0.46) +
  geom_sf(data = inv_outline, fill = NA, colour = "#C0392B", linewidth = 0.24) +
  scale_fill_viridis_c(option = "G", direction = -1, limits = c(0, 1),
                       breaks = c(0, 0.5, 1), labels = c("0\ndensest", "0.5", "1\nsparsest"),
                       name = "Surveillance gap") +
  theme_map()

# THREE distinct states, drawn as three layers rather than forced onto one continuous
# scale. A log fill cannot represent a zero, and flooring the zeros at an arbitrary small
# value (the first version of this panel used pmax(rate, 0.05)) would draw three zones
# that investigated NO alert as though they had investigated a few. The two zero states
# are also not the same state: "no laboratory-resolved alert at all" and "resolved alerts,
# all of them confirmed cases" say different things about the surveillance system.
ZERO_NEG_FILL <- "#E7DCC1"
.no_lab   <- dplyr::filter(shp_a, !ever_investigated)
.zero_neg <- dplyr::filter(shp_a, ever_investigated, neg_alerts == 0)
.some_neg <- dplyr::filter(shp_a, ever_investigated, neg_alerts > 0)
stopifnot(nrow(.no_lab) + nrow(.zero_neg) + nrow(.some_neg) == nrow(shp_a))
pA1b <- ggplot() +
  geom_sf(data = .no_lab,   fill = NA_FILL,        colour = "white", linewidth = 0.06) +
  geom_sf(data = .zero_neg, fill = ZERO_NEG_FILL,  colour = "white", linewidth = 0.06) +
  geom_sf(data = .some_neg, aes(fill = neg_per_100k), colour = "white", linewidth = 0.06) +
  geom_sf(data = prov_outline, fill = NA, colour = "grey45", linewidth = 0.18) +
  geom_sf(data = inv_outline, fill = NA, colour = "white", linewidth = 0.46) +
  geom_sf(data = inv_outline, fill = NA, colour = "#C0392B", linewidth = 0.24) +
  scale_fill_viridis_c(option = "F", direction = -1, transform = "log10",
                       labels = scales::label_number(accuracy = 0.1),
                       name = "Discarded alerts\nper 100k") +
  ann_in(sprintf("grey: no laboratory-resolved alert (%d zones, %d of the %d never invaded)\ncream: resolved alerts, none discarded (%d zones)",
                 nrow(.no_lab),
                 sum(!zone_ascert$ever_investigated & !zone_ascert$invaded),
                 sum(!zone_ascert$invaded), nrow(.zero_neg)),
         x = -Inf, y = -Inf, hjust = 0, vjust = -0.15, size = 2.2) +
  theme_map()

prov_pos <- zone_ascert %>%
  dplyr::group_by(province) %>%
  dplyr::summarise(n_zones = dplyr::n(), n_invaded = sum(invaded),
                   confirmed = sum(n_confirmed), resolved = sum(alerts_resolved),
                   .groups = "drop")
.pw <- .wilson_ci_tbl(prov_pos$confirmed, prov_pos$resolved)
prov_pos <- prov_pos %>% dplyr::mutate(ratio = .pw$est, lo = .pw$lo, hi = .pw$hi)
N_PROV_SILENT <- sum(prov_pos$resolved == 0)
prov_active <- prov_pos %>% dplyr::filter(resolved > 0) %>%
  dplyr::mutate(province = forcats::fct_reorder(province, ratio))

pA1c <- ggplot(prov_active, aes(ratio, province)) +
  geom_vline(xintercept = 0, colour = FAINT, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                colour = MUTED, linewidth = 0.45) +
  geom_point(aes(size = resolved), colour = PT_BLUE) +
  geom_text(aes(x = hi, label = sprintf("  %d/%d", confirmed, resolved)),
            hjust = 0, size = 2.3, colour = MUTED, family = base_family) +
  scale_size_area(max_size = 3.6, name = "Resolved alerts",
                  breaks = c(10, 100, 1000, 10000), labels = scales::comma) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.30))) +
  # Rows are ordered by the point estimate, which for a proportion estimated from four
  # alerts is mostly noise, so the denominator is printed beside every row and the count
  # of thinly-supported provinces is stated on the axis. A5C carries the same disclosure
  # on its tick labels. The bars are Wilson intervals on the province's alerts, which treat
  # every alert as an independent trial: positivity that clusters by zone makes them too
  # narrow, so the axis says what they assume.
  labs(x = sprintf("Confirmed share of laboratory-resolved alerts (labels: confirmed/resolved)\nordered by the point estimate; %d of the %d rest on fewer than %d resolved alerts\nbars: %.0f%% Wilson intervals, every alert treated as independent",
                   sum(prov_active$resolved < 100), nrow(prov_active), 100L, 100 * CONF_LEVEL),
       y = NULL) +
  ann_in(sprintf("%d of %d provinces had no\nlaboratory-resolved alert",
                 N_PROV_SILENT, nrow(prov_pos)), x = -Inf, y = Inf, vjust = 1.1) +
  theme_pub() + theme(panel.grid.major.y = element_blank())

# Correlation matrix on ONE support: the exogenous proxies plus epicentre travel time,
# all defined for every zone that has a road route (Spearman, since several are skewed).
cm_vars <- c(PROXY_EXO$var, "epi_travel_h")
cm_lab  <- c(gsub("\n", " ", PROXY_EXO$label), "Epicentre travel time (h)")
cm_lab  <- stringr::str_wrap(gsub("\\s*\\(.*\\)$", "", cm_lab), 12)
cm_dat  <- zone_ascert %>% dplyr::select(dplyr::all_of(cm_vars)) %>% as.data.frame()
cm_ok   <- stats::complete.cases(cm_dat)
CM_N    <- sum(cm_ok)
cm <- stats::cor(cm_dat[cm_ok, , drop = FALSE], method = "spearman")
cm_df <- as.data.frame(as.table(cm)) %>%
  setNames(c("v1", "v2", "rho")) %>%
  dplyr::mutate(i = match(v1, cm_vars), j = match(v2, cm_vars)) %>%
  dplyr::filter(i >= j) %>%
  dplyr::mutate(v1 = factor(cm_lab[i], levels = cm_lab),
                v2 = factor(cm_lab[j], levels = rev(cm_lab)))

pA1d <- ggplot(cm_df, aes(v1, v2, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rho),
                colour = abs(rho) > 0.6), size = 2.35, family = base_family) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), breaks = c(-1, 0, 1),
                       name = sprintf("Spearman rho\n(n = %d zones)", CM_N)) +
  scale_x_discrete(position = "top") + labs(x = NULL, y = NULL) +
  # The -1 is definitional, not a finding: the surveillance gap IS 1 minus the percentile
  # rank of facility density. Both are nonetheless carried through the regressions, where
  # they are different linear predictors (a percentile and a log are monotone-related but
  # not linearly), so a result that held on only one of them would be a rank artefact.
  ann_in("rho = -1 by construction: the surveillance gap\nis 1 - the percentile rank of facility density",
         x = -Inf, y = -Inf, hjust = 0, vjust = -0.4, size = 2.1) +
  theme_pub() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 0, vjust = 0, size = 6.0),
        axis.text.y = element_text(size = 6.0), legend.position = "right",
        legend.key.height = unit(16, "pt"), legend.key.width = unit(7, "pt"))

FigA1 <- (pA1a | pA1b) / (pA1c | pA1d) +
  patchwork::plot_layout(heights = c(1.05, 0.95)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(pA1a, "FigureA1a_map_surveillance_gap", 4.2, 3.8)
save_dual(pA1b, "FigureA1b_map_alert_effort",     4.2, 3.8)
save_dual(pA1c, "FigureA1c_province_positivity",  4.2, 3.2)
save_dual(pA1d, "FigureA1d_proxy_correlations",   4.2, 3.2)
save_dual(FigA1, "FigureA1", 9.6, 8.0, dir = FIG_DIR)

# -----------------------------------------------------------------------------
# 6. FIGURE A2 — does thin surveillance delay the first confirmed case?
#    The support is the INVADED zones: every one of them has an arrival time, so the
#    outbreak-conditioned proxies are comparable across them here.
#
#    TWO PANELS. The added-variable plot and the added-variable forest were dropped
#    from the figure by editorial decision; the regressions behind them are still
#    fitted below and still published, as A2_arrival_avp_surveillance.csv and
#    A2_arrival_added_variable.csv and in the report, so no number is lost with the
#    panels. `avp_frame()` is retained for the same reason: it defines the residual
#    pair the published slope refers to.
# -----------------------------------------------------------------------------
message("\n[ascert] Figure A2 — surveillance and the timing of the first confirmed case")

# THE ARRIVAL FITS' SUPPORT IS 43's: an invaded zone whose first onset falls after the arrival
# fit window (in_fit_window FALSE, the right-truncated tail) enters no arrival fit there, nor in
# Figure 1C. Here it matters doubly: a tail zone already confirmed is the quickly reported subset,
# the very thing this panel tests.
inv_all <- zone_ascert %>% dplyr::filter(invaded, is.finite(arrival_days))
stopifnot(nrow(inv_all) == sum(zone_ascert$invaded))
inv <- inv_all %>% dplyr::filter(in_fit_window %in% TRUE)
message(sprintf("[ascert] %d invaded zone(s) in the right-truncated tail are left out of the arrival fits.",
                nrow(inv_all) - nrow(inv)))
inv_ne <- inv %>% dplyr::filter(!health_zone %in% EPICENTRE_ZONES)
message(sprintf("[ascert] arrival-time support: %d invaded zone(s); %d excluding the epicentre region",
                nrow(inv), nrow(inv_ne)))

# Added-variable plot: BOTH axes residualised on epicentre travel time, so the fitted
# slope IS the coefficient of the proxy in lm(arrival ~ travel + proxy) — a partial
# residual plot (raw x against residual y) would draw a line of a different slope from
# the statistic annotated beside it.
avp_frame <- function(d, prox) {
  ok <- is.finite(d$arrival_days) & is.finite(d$epi_travel_h) & is.finite(d[[prox]])
  dd <- d[ok, , drop = FALSE]
  tibble::tibble(zone = dd$health_zone,
                 ry = stats::resid(stats::lm(dd$arrival_days ~ dd$epi_travel_h)),
                 rx = stats::resid(stats::lm(dd[[prox]] ~ dd$epi_travel_h)))
}
avp_sg <- avp_frame(inv, "surveillance_gap")
fit_sg <- av_slope(inv$arrival_days, inv$epi_travel_h, inv$surveillance_gap)
rho_sg <- spearman_ci(avp_sg$rx, avp_sg$ry)

# Added-variable forest over every proxy, with and without the three epicentre zones.
av_forest <- purrr::pmap_dfr(
  list(PROXY_ALL$var, PROXY_ALL$label, PROXY_ALL$kind),
  function(v, lab, kd) dplyr::bind_rows(
    dplyr::mutate(av_slope(inv$arrival_days,    inv$epi_travel_h,    inv[[v]]),
                  set = "All invaded zones"),
    dplyr::mutate(av_slope(inv_ne$arrival_days, inv_ne$epi_travel_h, inv_ne[[v]]),
                  set = "Excluding the epicentre region")) %>%
    dplyr::mutate(var = v, label = gsub("\n", " ", lab), kind = kd))
av_forest <- av_forest %>%
  dplyr::mutate(set = factor(set, levels = c("All invaded zones", "Excluding the epicentre region")),
                estimable = is.finite(slope_sd))
.ord <- av_forest %>% dplyr::filter(set == "All invaded zones") %>%
  dplyr::arrange(dplyr::desc(dplyr::coalesce(slope_sd, -Inf))) %>% dplyr::pull(label)
av_forest <- av_forest %>%
  dplyr::mutate(label_f = factor(label, levels = rev(.ord)),
                lab_n = sprintf("n = %d", n))

# Arrival time against the zone's own onset -> confirmation delay.
del <- inv %>% dplyr::filter(is.finite(delay_median))
rho_del <- spearman_ci(del$delay_median, del$arrival_days)
fit_del <- av_slope(del$arrival_days, del$epi_travel_h, del$delay_median)
# The drawn line is an UNADJUSTED ordinary least-squares fit, which is neither of the two
# statistics annotated beside it, so its own slope test is computed and annotated as well.
# It is then drawn dashed and muted when that slope is not distinguishable from zero: a
# confident red trend through a null relationship is the one thing this panel must not
# suggest, since its point is that detection delay does NOT order arrival.
.ols_del <- stats::lm(arrival_days ~ delay_median, data = del)
.ols_p   <- summary(.ols_del)$coefficients["delay_median", "Pr(>|t|)"]
.ols_sig <- is.finite(.ols_p) && .ols_p < 0.05
.ols_col <- if (.ols_sig) FIT_RED else MUTED
.ols_lty <- if (.ols_sig) "solid" else "22"

# The annotation sits in a band reserved above the data (the generous top expansion), and the
# zone labels are kept below that band (ggrepel's ylim), so neither a point nor a label can sit
# under the annotation whatever the data.
.a2a_top <- max(del$arrival_days) + 0.06 * diff(range(del$arrival_days))
pA2a <- ggplot(del, aes(delay_median, arrival_days)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = .ols_col, fill = .ols_col,
              alpha = 0.10, linewidth = 0.6, linetype = .ols_lty) +
  geom_point(aes(size = n_delay), colour = PT_BLUE, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = health_zone), size = 2.05, colour = MUTED,
                           family = base_family, max.overlaps = 8, seed = RANDOM_SEED,
                           min.segment.length = 0.3, segment.colour = FAINT,
                           segment.size = 0.25, ylim = c(NA, .a2a_top)) +
  scale_size_area(max_size = 3.4, name = "Cases used") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.48))) +
  labs(x = sprintf("Median onset to laboratory confirmation (days, capped at %d)", DELAY_CAP_DAYS),
       y = .arr_axis) +
  ann_in(sprintf("Spearman rho %s   %s\nfitted line (unadjusted OLS): %s%s\nadjusted for travel time: %s days per SD, %s\nn = %d zones with >=5 native-onset confirmed cases,\neach with >=%d d of follow-up",
                 fmt_ci(rho_del$rho, rho_del$lo, rho_del$hi), fmt_p(rho_del$p),
                 fmt_p(.ols_p), if (.ols_sig) "" else " (dashed: not distinguishable from zero)",
                 fmt_ci(fit_del$slope_sd, fit_del$lo, fit_del$hi, 1), fmt_p(fit_del$p),
                 nrow(del), DELAY_CAP_DAYS)) +
  theme_pub()

# Arrival-time ECDF by surveillance-gap tercile.
ter_lab <- c("Lower third\n(densest facilities)", "Middle third", "Upper third\n(sparsest facilities)")
inv_t <- inv %>% dplyr::mutate(ter = factor(ter_lab[dplyr::ntile(surveillance_gap, 3)],
                                            levels = ter_lab))
ter_med <- inv_t %>% dplyr::group_by(ter) %>%
  dplyr::summarise(n = dplyr::n(), med = stats::median(arrival_days), .groups = "drop")

pA2b <- ggplot(inv_t, aes(arrival_days, colour = ter)) +
  stat_ecdf(geom = "step", linewidth = 0.7, pad = FALSE) +
  geom_vline(data = ter_med, aes(xintercept = med, colour = ter),
             linetype = "22", linewidth = 0.4, show.legend = FALSE) +
  scale_colour_manual(values = unname(OKABE[c(1, 5, 2)]), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = .arr_axis,
       y = "Share of invaded zones arrived") +
  ann_in(paste(sprintf("%s: median %.0f d (n = %d)",
                       gsub("\n", " ", ter_med$ter), ter_med$med, ter_med$n), collapse = "\n"),
         x = Inf, y = -Inf, hjust = 1, vjust = -0.3) +
  theme_pub() + theme(legend.position = "top", legend.direction = "horizontal")

FigA2 <- (pA2a | pA2b) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(pA2a, "FigureA2a_arrival_vs_delay",      4.6, 3.6)
save_dual(pA2b, "FigureA2b_arrival_ecdf_tercile",  4.6, 3.6)
save_dual(FigA2, "FigureA2", 9.4, 3.9, dir = FIG_DIR)

# -----------------------------------------------------------------------------
# 7. FIGURE A3 — are invaded and never-invaded zones comparable?
#    Only the EXOGENOUS proxies enter this figure. An invaded zone has a confirmed
#    case by definition, so every outbreak-conditioned quantity (confirmation
#    ratio, alert counts, confirmation delay) is mechanically tied to the outcome
#    and a contrast on it would measure the definition, not the surveillance system.
# -----------------------------------------------------------------------------
message("\n[ascert] Figure A3 — comparability of invaded and never-invaded zones")

za_t <- zone_ascert %>% dplyr::filter(is.finite(epi_travel_h))
message(sprintf("[ascert] comparability support: %d zones with a road route (%d invaded)",
                nrow(za_t), sum(za_t$invaded)))

cd_sg <- cliffs_delta(zone_ascert$surveillance_gap[zone_ascert$invaded],
                      zone_ascert$surveillance_gap[!zone_ascert$invaded])
mw_sg <- suppressWarnings(stats::wilcox.test(surveillance_gap ~ invaded, data = zone_ascert))

pA3a <- ggplot(zone_ascert, aes(status, surveillance_gap, fill = status, colour = status)) +
  geom_violin(alpha = 0.18, linewidth = 0.3, width = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  geom_jitter(width = 0.11, height = 0, size = 0.7, alpha = 0.38, stroke = 0) +
  scale_fill_manual(values = INV_COL, guide = "none") +
  scale_colour_manual(values = INV_COL, guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.03, 0.32))) +
  labs(x = NULL, y = "Surveillance gap (facility density percentile)") +
  ann_in(sprintf("Cliff's delta %s\nMann-Whitney %s\nn = %d invaded, %d never invaded",
                 fmt_ci(cd_sg$delta, cd_sg$lo, cd_sg$hi), fmt_p(mw_sg$p.value),
                 cd_sg$n_a, cd_sg$n_b)) +
  theme_pub()

cd_forest <- purrr::pmap_dfr(
  list(c(PROXY_EXO$var, "epi_travel_h"),
       c(PROXY_EXO$label, "Epicentre travel time\n(hours)")),
  function(v, lab) cliffs_delta(zone_ascert[[v]][zone_ascert$invaded],
                                zone_ascert[[v]][!zone_ascert$invaded]) %>%
    dplyr::mutate(var = v, label = gsub("\n", " ", lab)))
cd_forest <- cd_forest %>%
  dplyr::arrange(delta) %>%
  dplyr::mutate(label_f = factor(label, levels = label))

pA3b <- ggplot(cd_forest, aes(delta, label_f)) +
  geom_vline(xintercept = 0, colour = FAINT, linewidth = 0.35, linetype = "22") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                colour = MUTED, linewidth = 0.45) +
  geom_point(colour = FIT_RED, size = 1.9) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5),
                     expand = expansion(mult = c(0.03, 0.03))) +
  # Every row is estimated on the same two groups bar the travel-time row, where the two
  # zones with no road route drop out, so the support is stated once rather than per row.
  # Right-aligned and on short lines: the long proxy labels leave the plotting area narrow,
  # and a centred title wider than it ran off the right edge of the composite.
  labs(x = sprintf("Cliff's delta, invaded vs never invaded\n(negative = lower in invaded zones)\n%d invaded vs %d never-invaded zones\n(%d vs %d for travel time)",
                   sum(zone_ascert$invaded), sum(!zone_ascert$invaded),
                   sum(za_t$invaded), sum(!za_t$invaded)),
       y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank(),
                      axis.title.x = element_text(hjust = 1))

# The same contrast INSIDE quartiles of epicentre travel time. Invaded zones are not spread
# evenly across the country, so this shows in which quartiles the contrast can be estimated at
# all, and what it is there; the report states both from q_stat.
TT_Q <- 4L
Q_MIN_INV <- 5L      # invaded zones a quartile needs before its delta is estimated
za_q <- za_t %>%
  dplyr::mutate(qi = dplyr::ntile(epi_travel_h, TT_Q)) %>%
  dplyr::group_by(qi) %>%
  dplyr::mutate(q_lab = sprintf("Q%d: %.0f-%.0f h", qi[1], min(epi_travel_h), max(epi_travel_h))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(q_lab = factor(q_lab, levels = unique(q_lab[order(qi)])))
q_stat <- za_q %>%
  dplyr::group_by(q_lab) %>%
  dplyr::group_modify(function(g, k) {
    a <- g$surveillance_gap[g$invaded]; b <- g$surveillance_gap[!g$invaded]
    cd <- cliffs_delta(a, b)
    # Below Q_MIN_INV the delta is NOT ESTIMATED, in the CSV as on the panel: cliffs_delta()
    # itself only refuses fewer than 2, so the value is blanked here rather than published
    # beside a panel that says it does not exist.
    if (length(a) < Q_MIN_INV) cd <- dplyr::mutate(cd, delta = NA_real_, lo = NA_real_, hi = NA_real_)
    dplyr::mutate(cd, qi = g$qi[1], note = if (length(a) < Q_MIN_INV)
      sprintf("%d invaded zone%s\n(delta not estimated)", length(a), if (length(a) == 1L) "" else "s")
      else sprintf("delta %s\nn = %d / %d", fmt_ci(cd$delta, cd$lo, cd$hi), cd$n_a, cd$n_b))
  }) %>% dplyr::ungroup()

pA3c <- ggplot(za_q, aes(status, surveillance_gap, colour = status, fill = status)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.16, linewidth = 0.4) +
  geom_jitter(width = 0.13, height = 0, size = 0.6, alpha = 0.4, stroke = 0) +
  geom_text(data = q_stat, aes(x = 1.5, y = 1.03, label = note), inherit.aes = FALSE,
            size = 2.05, colour = MUTED, family = base_family, vjust = 0, lineheight = 1.05) +
  facet_wrap(~ q_lab, nrow = 1) +
  scale_fill_manual(values = INV_COL, guide = "none") +
  scale_colour_manual(values = INV_COL, guide = "none") +
  scale_y_continuous(limits = c(0, 1.24), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0))) +
  scale_x_discrete(labels = c("Never\ninvaded", "Invaded")) +
  labs(x = "Quartile of road travel time from the epicentre region",
       y = "Surveillance gap") +
  theme_pub() + theme(axis.text.x = element_text(size = 6.6))

# BOTH models are fit on `za_t`, the zones with a road route. Fitting the unadjusted model
# on the full spine and the adjusted one on that subset would compare coefficients
# estimated on different samples AND standardised by SDs computed on different samples, so
# the pair would not be the adjustment contrast it is drawn as.
or_forest <- purrr::pmap_dfr(
  list(PROXY_EXO$var, PROXY_EXO$label),
  function(v, lab) dplyr::bind_rows(
    dplyr::mutate(logit_or(za_t$invaded, za_t[[v]]), model = "Unadjusted"),
    dplyr::mutate(logit_or(za_t$invaded, za_t[[v]], base = za_t$epi_travel_h),
                  model = "Adjusted for epicentre travel time")) %>%
    dplyr::mutate(var = v, label = gsub("\n", " ", lab)))
or_forest <- or_forest %>%
  dplyr::mutate(model = factor(model, levels = c("Unadjusted", "Adjusted for epicentre travel time")))
.ord_or <- or_forest %>% dplyr::filter(model == "Unadjusted") %>%
  dplyr::arrange(dplyr::coalesce(or_sd, Inf)) %>% dplyr::pull(label)
or_forest <- or_forest %>% dplyr::mutate(label_f = factor(label, levels = .ord_or))
OR_CI_TYPES <- sort(unique(stats::na.omit(or_forest$ci_type)))

pA3d <- ggplot(dplyr::filter(or_forest, is.finite(or_sd)),
               aes(or_sd, label_f, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = FAINT, linewidth = 0.35, linetype = "22") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                linewidth = 0.45, position = position_dodge(width = 0.55)) +
  geom_point(size = 1.7, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = unname(OKABE[c(1, 2)]), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_continuous(transform = "log10", breaks = c(0.25, 0.5, 1, 2, 4),
                     labels = c("0.25", "0.5", "1", "2", "4"),
                     expand = expansion(mult = c(0.06, 0.06))) +
  # The support note goes in the axis title, not an in-panel annotation: this axis is
  # log-transformed, and an annotation anchored at -Inf has no finite position on it.
  labs(x = sprintf("Odds ratio for being invaded, per SD of the proxy\n(%s intervals; %d invaded of %d zones with a road route)",
                   paste(OR_CI_TYPES, collapse = "/"), sum(za_t$invaded), nrow(za_t)),
       y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank(),
                      legend.position = "top", legend.direction = "horizontal")

FigA3 <- (pA3a | pA3b) / pA3c / pA3d +
  patchwork::plot_layout(heights = c(1, 0.85, 0.85)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(pA3a, "FigureA3a_surveillance_by_status",   4.0, 3.4)
save_dual(pA3b, "FigureA3b_cliffs_delta_forest",      4.6, 3.4)
save_dual(pA3c, "FigureA3c_contrast_within_travel_q", 8.2, 3.0)
save_dual(pA3d, "FigureA3d_logistic_or",              5.6, 3.0)
save_dual(FigA3, "FigureA3", 9.4, 9.6, dir = FIG_DIR)

# -----------------------------------------------------------------------------
# 8. FIGURE A4 — does the forecast rank by surveillance?
#    All ranks are taken within a (horizon, fold) watch-list, so they are exactly
#    invariant to the probability scale and to any monotone recalibration.
# -----------------------------------------------------------------------------
message("\n[ascert] Figure A4 — the forecast and the surveillance gradient")

# (H_MAIN, the main horizon, is set in section 4.)
# Zone-level correlation, NOT zone-fold level: the at-risk zone-folds are repeated
# measurements of the same zones, roughly one per fold per horizon, so a correlation over
# them would treat that replication as independent information and report an interval far
# narrower than the data support.
rank_zone <- lfo_f %>%
  dplyr::group_by(horizon, health_zone, surveillance_gap) %>%
  dplyr::summarise(mean_rank = mean(rank_avg), n_folds = dplyr::n(), .groups = "drop")
rho_rank <- rank_zone %>%
  dplyr::group_by(horizon) %>%
  dplyr::group_modify(function(g, k) spearman_ci(g$surveillance_gap, g$mean_rank)) %>%
  dplyr::ungroup()

d_main <- lfo_f %>% dplyr::filter(horizon == H_MAIN)
dec_med <- d_main %>%
  dplyr::mutate(dec = dplyr::ntile(surveillance_gap, 10)) %>%
  dplyr::group_by(dec) %>%
  dplyr::summarise(x = stats::median(surveillance_gap), y = stats::median(rank_avg),
                   .groups = "drop")
.rr <- rho_rank %>% dplyr::filter(horizon == H_MAIN)

pA4a <- ggplot(d_main, aes(surveillance_gap, rank_avg)) +
  geom_bin2d(bins = 34) +
  geom_line(data = dec_med, aes(x, y), colour = "white", linewidth = 1.1, inherit.aes = FALSE) +
  geom_line(data = dec_med, aes(x, y), colour = INK, linewidth = 0.55, inherit.aes = FALSE) +
  geom_point(data = dplyr::filter(d_main, is_new_invasion == 1),
             colour = FIT_RED, size = 1.5, shape = 21, fill = "white", stroke = 0.7) +
  scale_fill_viridis_c(option = "G", direction = -1, transform = "log10",
                       name = "At-risk\nzone-forecasts", labels = scales::comma) +
  scale_y_reverse(expand = expansion(mult = c(0.04, 0.03))) +
  # The panel is filled edge to edge, and an annotation anchored at +/-Inf has no finite
  # position on a REVERSED scale, so the legend text lives in the axis title instead.
  labs(x = sprintf("Surveillance gap (facility density percentile)\nblack line: median rank by decile;  red rings: %d realised invasion zone-forecasts (%d distinct zones)\nzone-level Spearman rho %s, %s (n = %d zones)",
                   sum(d_main$is_new_invasion),
                   dplyr::n_distinct(d_main$health_zone[d_main$is_new_invasion == 1]),
                   fmt_ci(.rr$rho, .rr$lo, .rr$hi), fmt_p(.rr$p), .rr$n),
       y = sprintf("Rank in the %d-week watch-list (1 = highest risk)", H_MAIN)) +
  theme_pub() + theme(legend.position = "right",
                      legend.key.height = unit(16, "pt"), legend.key.width = unit(7, "pt"))

# ONE POINT PER INVADED ZONE. At h = 2 the outcome windows of consecutive origins overlap,
# so a single invasion appears as a positive row in up to two folds — the gap between
# `n_invasions` and `n_invasion_events` in invasion_evaluation.csv. Correlating over the
# rows would treat that duplication as independent information and report an interval the
# data do not support, so each zone enters once, at its mean rank over the rounds in which
# it was the event. Both counts are printed on the panel.
ev_rows <- lfo_f %>% dplyr::filter(is_new_invasion == 1)
ev <- ev_rows %>%
  dplyr::group_by(horizon, health_zone, surveillance_gap) %>%
  dplyr::summarise(rank_avg = mean(rank_avg), n_rows = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(h_lab = factor(sprintf("%d-week horizon", horizon),
                               levels = sprintf("%d-week horizon", sort(unique(horizon)))))
.ev_rows_n <- ev_rows %>% dplyr::count(horizon, name = "n_rows_total")
# The drawn line is an OLS fit on the TRANSFORMED y axis (the panel is log10 in rank), so
# its slope test is computed the same way, and the line is dashed and muted wherever that
# slope is not distinguishable from zero. Drawing a confident trend through a null
# relationship is exactly the misreading this panel exists to prevent.
.ols_ev <- ev %>% dplyr::group_by(horizon) %>%
  dplyr::group_modify(function(g, k) {
    m <- stats::lm(log10(rank_avg) ~ surveillance_gap, data = g)
    tibble::tibble(ols_p = summary(m)$coefficients["surveillance_gap", "Pr(>|t|)"])
  }) %>% dplyr::ungroup() %>%
  dplyr::mutate(sig = is.finite(ols_p) & ols_p < 0.05)
ev <- ev %>% dplyr::left_join(.ols_ev, by = "horizon")
rho_ev <- ev %>% dplyr::group_by(h_lab, horizon) %>%
  dplyr::group_modify(function(g, k) spearman_ci(g$surveillance_gap, g$rank_avg)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(.ev_rows_n, by = "horizon") %>%
  dplyr::left_join(.ols_ev, by = "horizon") %>%
  # Short lines: each facet is about half a panel wide, and the one-line-per-statistic
  # version ran past the facet's right edge.
  dplyr::mutate(lab = sprintf("Spearman rho %s\n%s\nOLS line (log rank): %s%s\n%d invaded zones (%d zone-forecasts)",
                              fmt_ci(rho, lo, hi), fmt_p(p), fmt_p(ols_p),
                              ifelse(sig, "", ", dashed"), n, n_rows_total))

pA4b <- ggplot(ev, aes(surveillance_gap, rank_avg)) +
  geom_smooth(aes(colour = sig, fill = sig, linetype = sig),
              method = "lm", formula = y ~ x, alpha = 0.10, linewidth = 0.6) +
  scale_colour_manual(values = c(`TRUE` = FIT_RED, `FALSE` = MUTED), guide = "none") +
  scale_fill_manual(values = c(`TRUE` = FIT_RED, `FALSE` = MUTED), guide = "none") +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), guide = "none") +
  geom_point(colour = PT_BLUE, size = 1.5, alpha = 0.85) +
  # A borderless label justified at 0, so its lines align on the left; the padding, not an
  # out-of-range hjust, keeps it off the panel edge.
  geom_label(data = rho_ev, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
             hjust = 0, vjust = 1, size = 2.1, colour = MUTED, family = base_family,
             lineheight = 1.05, fill = NA, linewidth = 0, label.padding = unit(3, "pt")) +
  facet_wrap(~ h_lab, nrow = 1) +
  # the top expansion is the band the four-line annotation sits in, above every point
  scale_y_continuous(transform = "log10", breaks = c(1, 3, 10, 30, 100, 300),
                     expand = expansion(mult = c(0.06, 0.42))) +
  labs(x = "Surveillance gap of the invaded zone",
       y = "Mean rank of the invaded zone\nin the round(s) it was the event") +
  theme_pub()

# Within-stratum discrimination, on the headline basis (strat_discrimination(): all rounds,
# every metric scored within a round). Each stratum is scored against ITS OWN per-round
# base rates and chance level, and the bootstrap resamples whole zones inside the stratum.
strat_tab <- lfo_f %>%
  dplyr::group_by(horizon, sg_stratum) %>%
  dplyr::group_modify(function(g, k) strat_discrimination(g)) %>%
  dplyr::ungroup()
message("[ascert] within-stratum discrimination (all rounds, scored within each round):")
print(as.data.frame(strat_tab %>% dplyr::select(horizon, sg_stratum, n_rounds, n, n_pos, n_pos_zones,
                                                auc_pr_date, ap_lo, ap_hi, ap_chance,
                                                auc_pr_skill_date, skill_lo, skill_hi, skill_chance,
                                                skill_ceiling, auc_roc_date, roc_lo, roc_hi)),
      digits = 3)
# Figure A4C is scored on the RAW probability (see strat_discrimination()), unlike the rank
# panels, which use `PCOL`; a within-round score is the same on either scale. The basis the
# panel states is checked against the table rather than assumed.
A4C_PROB_COL <- unique(strat_tab$prob_column)
stopifnot(identical(A4C_PROB_COL, "p_invasion_raw"), identical(unique(strat_tab$rounds), "all"))
A4C_N_ROUNDS <- max(strat_tab$n_rounds)
A4C_SCALE_NOTE <- sprintf("scored within each round, over all %d rounds", A4C_N_ROUNDS)

M_SKILL <- "AUC-PR skill"; M_AP <- "Average precision"; M_ROC <- "AUC-ROC"
strat_long <- dplyr::bind_rows(
  strat_tab %>% dplyr::transmute(horizon, sg_stratum, n_pos, n_pos_zones, base_rate,
                                 metric = M_SKILL, ref = skill_chance, cap = skill_ceiling,
                                 est = auc_pr_skill_date, lo = skill_lo, hi = skill_hi),
  strat_tab %>% dplyr::transmute(horizon, sg_stratum, n_pos, n_pos_zones, base_rate,
                                 metric = M_AP, ref = ap_chance, cap = NA_real_,
                                 est = auc_pr_date, lo = ap_lo, hi = ap_hi),
  strat_tab %>% dplyr::transmute(horizon, sg_stratum, n_pos, n_pos_zones, base_rate,
                                 metric = M_ROC, ref = 0.5, cap = NA_real_,
                                 est = auc_roc_date, lo = roc_lo, hi = roc_hi)) %>%
  dplyr::mutate(h_lab = factor(sprintf("%d wk", horizon),
                               levels = sprintf("%d wk", sort(unique(horizon)))),
                metric = factor(metric, levels = c(M_SKILL, M_AP, M_ROC)))
#' One metric per plot rather than one facet per metric: AUC-PR skill is a RATIO spanning
#' an order of magnitude at this event count and is only readable on a log axis, while the
#' other two are probabilities on a linear one. `scales = "free_y"` frees the limits but
#' not the transform, so facetting them together would force one or the other.
#'
#' The CHANCE reference (what a uniformly random ordering of each round would score) is drawn
#' per group, not as one line, for AUC-PR skill and average precision: it depends on each
#' round's size and invasion count, so it differs between strata and horizons. Only AUC-ROC has
#' a common one (0.5). The skill panel additionally marks each stratum's ATTAINABLE MAXIMUM
#' (every round's average precision equal to 1): without it, a higher skill in the stratum with
#' the lower base rate would read as better discrimination when it is partly a higher ceiling.
.strat_panel <- function(metric_name, ylab, log_y, show_cap = FALSE, common_ref = TRUE,
                         xlab = "Forecast horizon") {
  d <- dplyr::filter(strat_long, metric == metric_name)
  dodge <- position_dodge(width = 0.5)
  p <- ggplot(d, aes(h_lab, est, colour = sg_stratum, shape = sg_stratum))
  p <- if (common_ref) {
    # A "common" reference must actually be common; two values here would draw two lines
    # and silently imply the null differs between strata when the caller said it does not.
    stopifnot(length(unique(d$ref)) == 1L)
    p + geom_hline(yintercept = d$ref[1], colour = FAINT, linewidth = 0.35, linetype = "22")
  }
  else
    p + geom_point(aes(y = ref), shape = 45, size = 5, alpha = 0.75,
                   position = dodge, show.legend = FALSE)
  if (show_cap)
    p <- p + geom_point(aes(y = cap), shape = 45, size = 5, alpha = 0.55,
                        position = dodge, show.legend = FALSE)
  p +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.45, position = dodge) +
    geom_point(size = 2, position = dodge) +
    geom_text(aes(y = hi, label = sprintf("%d", n_pos_zones)), vjust = -0.8, size = 2.1,
              position = dodge, show.legend = FALSE) +
    scale_colour_manual(values = unname(OKABE[c(1, 2)]), name = NULL) +
    scale_shape_manual(values = c(16, 17), name = NULL) +
    labs(x = xlab, y = ylab) + theme_pub() +
    theme(legend.position = "top", legend.direction = "horizontal") +
    (if (log_y) scale_y_continuous(transform = "log10", breaks = c(1, 3, 10, 30, 100, 300),
                                   expand = expansion(mult = c(0.08, 0.16)))
     else scale_y_continuous(expand = expansion(mult = c(0.12, 0.16))))
}
# The "numbers above" convention and the probability scale apply to all three panels, so
# they sit under the middle one, where they read as a note on the row rather than on the
# panel they happen to touch.
pA4c <- (.strat_panel(M_SKILL, "AUC-PR skill (x base rate)\ndashes: random-ordering chance (lower)\nand attainable maximum (upper)",
                      TRUE, show_cap = TRUE, common_ref = FALSE) |
         .strat_panel(M_AP, "Average precision\ndashes: random-ordering chance\n(the stratum's own)", FALSE,
                      common_ref = FALSE,
                      # the basis note is no wider than the line above it (48 characters),
                      # so it cannot widen the panel
                      xlab = paste0("Forecast horizon\nnumbers above each point: distinct invaded zones\n",
                                    A4C_SCALE_NOTE)) |
         .strat_panel(M_ROC, "AUC-ROC\ndashes: 0.5 (chance)", FALSE)) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.direction = "horizontal")

# Watch-list recall by stratum. The top-K list is NATIONAL (ranks come from the whole
# at-risk set of that round, tie-broken to the worst position exactly as
# compute_detection_curve does); the recall is then the share of THAT stratum's invasions
# the national list caught. Pooled over folds rather than averaged per fold, because most
# folds contribute no event at all to a stratum and the per-fold share is undefined there;
# the interval is the zone-cluster bootstrap in recall_curve_boot(), which is what the
# overlapping 2-week outcome windows require.
KS <- 1:25
rec_strat <- ev_rows %>%
  dplyr::group_by(horizon, sg_stratum) %>%
  dplyr::group_modify(function(g, key) recall_curve_boot(g, KS)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(h_lab = factor(sprintf("%d-week horizon", horizon),
                               levels = sprintf("%d-week horizon", sort(unique(horizon)))))
rand_ref <- lfo_f %>% dplyr::group_by(horizon, fold_id) %>%
  dplyr::summarise(n_atrisk = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(horizon) %>%
  dplyr::summarise(n_atrisk = mean(n_atrisk), .groups = "drop") %>%
  tidyr::expand_grid(k = KS) %>%
  dplyr::mutate(recall = pmin(k / n_atrisk, 1),
                h_lab = factor(sprintf("%d-week horizon", horizon),
                               levels = sprintf("%d-week horizon", sort(unique(horizon)))))

rec_band <- rec_strat %>%
  dplyr::group_by(horizon, sg_stratum, h_lab) %>%
  dplyr::group_modify(function(g, key) step_path(g)) %>%
  dplyr::ungroup()
pA4d <- ggplot(rec_strat, aes(k, recall, colour = sg_stratum, fill = sg_stratum)) +
  geom_ribbon(data = rec_band, aes(ymin = lo, ymax = hi), alpha = 0.14, colour = NA) +
  geom_step(linewidth = 0.7, direction = "hv") +
  geom_line(data = rand_ref, aes(k, recall), colour = MUTED, linetype = "22",
            linewidth = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ h_lab, nrow = 1) +
  scale_colour_manual(values = unname(OKABE[c(1, 2)]), name = NULL) +
  scale_fill_manual(values = unname(OKABE[c(1, 2)]), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = sprintf("Size K of the national weekly watch-list\ndashed: random targeting; bands are %.0f%% zone-cluster bootstrap intervals",
                   100 * CONF_LEVEL),
       y = "Share of that stratum's\ninvasions caught") +
  theme_pub() + theme(legend.position = "top", legend.direction = "horizontal")

# pA4c is itself a patchwork; wrapped so the composite tags it as ONE panel (C) instead
# of handing separate tags to its three children. It gets a full row: three metric panels
# squeezed into half a row would make the interval whiskers unreadable.
FigA4 <- (pA4a | pA4b) / patchwork::wrap_elements(pA4c) / pA4d +
  patchwork::plot_layout(heights = c(1.15, 0.95, 0.85)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(pA4a, "FigureA4a_rank_vs_surveillance",   4.8, 3.6)
save_dual(pA4b, "FigureA4b_rank_of_truth",          4.6, 3.6)
save_dual(pA4c, "FigureA4c_stratified_skill",       8.6, 3.6)
save_dual(pA4d, "FigureA4d_recall_by_stratum",      8.6, 3.2)
save_dual(FigA4, "FigureA4", 10.0, 11.4, dir = FIG_DIR)

# -----------------------------------------------------------------------------
# 9. FIGURE A5 — province test positivity, detection opportunity, and the ranking
# -----------------------------------------------------------------------------
message("\n[ascert] Figure A5 — province positivity and the ranking")

# A. Was the system looking where the model pointed? For every at-risk zone-forecast,
#    whether that zone had ALREADY had a discarded (not-a-case) alert investigated by the
#    forecast origin — a surveillance-effort measure that cannot be created by the zone's
#    own later invasion, unlike "has a resolved alert" measured at the end of the outbreak.
BANDS <- c(0, 5, 10, 15, 25, 50, 100, 250, Inf)
BAND_LAB <- c("1-5", "6-10", "11-15", "16-25", "26-50", "51-100", "101-250", "251+")
# THE INTERVAL IS A ZONE-CLUSTER BOOTSTRAP, not a Wilson interval on zone-forecasts: each zone
# sits in every round in which it is at risk, with a flag that persists across them, so the
# rows are not independent trials (share_by_group_boot()).
.wb_rows <- lfo_f %>%
  dplyr::filter(horizon == H_MAIN) %>%
  dplyr::mutate(band = cut(rank_avg, breaks = BANDS, labels = BAND_LAB, right = TRUE))
.bw <- share_by_group_boot(.wb_rows, "band", "watched_prior")
watch_band <- .wb_rows %>%
  dplyr::group_by(band) %>%
  dplyr::summarise(n = dplyr::n(), watched = sum(watched_prior),
                   events = sum(is_new_invasion), .groups = "drop") %>%
  dplyr::left_join(.bw %>% dplyr::transmute(band = factor(group, levels = levels(.wb_rows$band)),
                                            n_zones, share = est, lo, hi, interval),
                   by = "band")
stopifnot(nrow(watch_band) > 0L, !anyNA(watch_band$share),
          isTRUE(all.equal(watch_band$share, watch_band$watched / watch_band$n)))

pA5a <- ggplot(watch_band, aes(band, share)) +
  geom_col(fill = PT_BLUE, alpha = 0.78, width = 0.68) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = INK, linewidth = 0.4) +
  geom_text(aes(y = hi, label = scales::comma(n)), vjust = -0.7, size = 2.1,
            colour = MUTED, family = base_family) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.16))) +
  # Three lines: the interval construction has to be on the panel (the figures carry no
  # caption), and one line would run past the right edge at this width.
  labs(x = sprintf(paste0("Rank band in the %d-week watch-list\nnumbers above: at-risk zone-forecasts in the band\n",
                          "bars: %.0f%% zone-cluster bootstrap%s"),
                   H_MAIN, 100 * CONF_LEVEL,
                   if (any(grepl("^Wilson", watch_band$interval))) "; at 0% or 100%, Wilson on distinct zones" else ""),
       y = "Share with a discarded alert already\nnotified at the forecast origin") +
  theme_pub()

# B. The province confirmed share AS KNOWN AT EACH FORECAST ORIGIN (prov_share_o, section 4)
#    against how the model ranked that province's at-risk zones. One point per province, over
#    the zone-forecasts issued when the province already had a laboratory-resolved alert: x is
#    the mean of the as-known share over them (each forecast weighted once, as in y), y their
#    mean rank, and the bar the RANGE of the as-known share across those origins. The share is
#    not the end-of-outbreak one (Figure A1C), which counts confirmations made after the origin,
#    including those of the zones being ranked.
prov_rank <- lfo_f %>%
  dplyr::filter(horizon == H_MAIN, prov_resolved_o > 0L) %>%
  dplyr::group_by(province) %>%
  dplyr::summarise(share = mean(prov_share_o), share_min = min(prov_share_o),
                   share_max = max(prov_share_o),
                   resolved_min = min(prov_resolved_o), resolved_max = max(prov_resolved_o),
                   n_rounds = dplyr::n_distinct(cutoff), mean_rank = mean(rank_avg),
                   n_zonefolds = dplyr::n(), n_events = sum(is_new_invasion), .groups = "drop")
rho_prov <- spearman_ci(prov_rank$share, prov_rank$mean_rank)
# The panel table: the end-of-outbreak counts of Figure A1C beside the as-known quantities
# plotted here, for every province with an at-risk zone-forecast at the main horizon.
prov_join <- prov_pos %>%
  dplyr::filter(province %in% lfo_f$province[lfo_f$horizon == H_MAIN]) %>%
  dplyr::select(province, n_zones, n_invaded, confirmed_end = confirmed, resolved_end = resolved,
                share_end = ratio, share_end_lo = lo, share_end_hi = hi) %>%
  dplyr::left_join(prov_rank, by = "province")

pA5b <- ggplot(prov_rank, aes(share, mean_rank)) +
  geom_errorbar(aes(xmin = share_min, xmax = share_max), orientation = "y", width = 0,
                colour = FAINT, linewidth = 0.4) +
  geom_point(aes(size = n_zonefolds), colour = PT_BLUE, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = province), size = 2.2, colour = MUTED,
                           family = base_family, seed = RANDOM_SEED, max.overlaps = 20,
                           min.segment.length = 0.2, segment.colour = FAINT, segment.size = 0.25) +
  scale_size_area(max_size = 4, name = "At-risk\nzone-forecasts", labels = scales::comma) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_reverse() +
  # Short lines: at this panel width a long line is clipped at the right edge.
  labs(x = sprintf(paste0("Province confirmed share of laboratory-resolved alerts,\n",
                          "as known at each forecast origin (mean over its zone-forecasts;\n",
                          "bars: range across those origins)\nSpearman rho %s, %s (%d provinces)"),
                   fmt_ci(rho_prov$rho, rho_prov$lo, rho_prov$hi), fmt_p(rho_prov$p), rho_prov$n),
       y = sprintf("Mean rank over the province's at-risk\nzone-forecasts (%d-week watch-list)", H_MAIN)) +
  theme_pub()

# C. The ranking distribution province by province, ordered by the as-known share of panel B,
#    with every zone-forecast issued when its province had NO laboratory-resolved alert yet
#    pooled into one group — there is no share to order them by, and showing many near-identical
#    boxes would imply a precision the data do not carry. A province therefore appears in its
#    own group for the rounds after its first resolved alert and in the pooled group before it.
SILENT_LAB <- "No laboratory-\nresolved alert\nby the origin"
prov_levels <- prov_rank %>% dplyr::arrange(dplyr::desc(share)) %>% dplyr::pull(province)
# The tick label carries the share AND its denominator. Ordering provinces by a proportion
# estimated, at the early origins, from a handful of resolved alerts would otherwise read as a
# ranking, which the printed range of denominators makes impossible to miss.
.res_rng <- function(a, b) ifelse(a == b, scales::comma(a),
                                  paste0(scales::comma(a), "-", scales::comma(b)))
.prov_tick <- prov_rank %>%
  dplyr::transmute(province,
                   tick = sprintf("%s\n%.0f%%\n(%s resolved)", province, 100 * share,
                                  .res_rng(resolved_min, resolved_max)))
prov_tick_lab <- stats::setNames(c(.prov_tick$tick[match(prov_levels, .prov_tick$province)],
                                   SILENT_LAB), c(prov_levels, SILENT_LAB))
rank_prov <- lfo_f %>%
  dplyr::filter(horizon == H_MAIN) %>%
  dplyr::mutate(prov_grp = factor(dplyr::if_else(prov_resolved_o > 0L, province, SILENT_LAB),
                                  levels = c(prov_levels, SILENT_LAB)),
                outcome = factor(dplyr::if_else(is_new_invasion == 1,
                                                "Invaded within the window", "Not invaded"),
                                 levels = c("Not invaded", "Invaded within the window")))
stopifnot(!anyNA(rank_prov$prov_grp))
prov_n <- rank_prov %>% dplyr::group_by(prov_grp) %>%
  dplyr::summarise(n = dplyr::n(), ev = sum(is_new_invasion),
                   ev_zones = dplyr::n_distinct(health_zone[is_new_invasion == 1]),
                   .groups = "drop") %>%
  dplyr::mutate(lab = sprintf("%s zone-forecasts\n%d invaded zone%s", scales::comma(n),
                              ev_zones, ifelse(ev_zones == 1L, "", "s")))
# THE LAST GROUP CAN HOLD AN INVASION, and the axis text and the report say how many, and why,
# from the data. An invasion is a confirmed case AFTER the origin, so a zone in a province with no
# laboratory-resolved alert by the origin can be invaded in two ways, counted separately: its
# province's first resolved alert came after the origin (the zone has a line-list confirmed
# case), or the zone was confirmed only through sitrep top-up rows, which carry no test and are
# excluded from every alert count. The second is CHECKED zone by zone against the sitrep rows
# (zone_sitrep, 3.4), never inferred from the absence of a line-list confirmed case.
.silent_inv <- rank_prov %>%
  dplyr::filter(prov_grp == SILENT_LAB, is_new_invasion == 1) %>%
  dplyr::distinct(health_zone, province) %>%
  dplyr::left_join(dplyr::select(zone_ascert, health_zone, n_confirmed), by = "health_zone") %>%
  dplyr::left_join(zone_sitrep, by = "health_zone") %>%
  dplyr::mutate(n_sitrep_conf = dplyr::coalesce(n_sitrep_conf, 0L))
.silent_n_rows    <- sum(rank_prov$prov_grp == SILENT_LAB & rank_prov$is_new_invasion == 1)
.silent_n_ll      <- sum(.silent_inv$n_confirmed > 0L)
.silent_n_sitrep  <- sum(.silent_inv$n_confirmed == 0L & .silent_inv$n_sitrep_conf > 0L)
.silent_n_neither <- sum(.silent_inv$n_confirmed == 0L & .silent_inv$n_sitrep_conf == 0L)
.silent_note <- if (!nrow(.silent_inv)) {
  "the last group holds no invasion"
} else {
  paste0(sprintf("the last group holds %d invaded zone%s (%d zone-forecast%s): ",
                 nrow(.silent_inv), ifelse(nrow(.silent_inv) == 1L, "", "s"),
                 .silent_n_rows, ifelse(.silent_n_rows == 1L, "", "s")),
         paste(c(if (.silent_n_ll > 0L)
                   sprintf("%d in a province whose first laboratory-resolved alert came after the origin",
                           .silent_n_ll),
                 if (.silent_n_sitrep > 0L)
                   sprintf("%d confirmed only through sitrep top-up rows, which carry no laboratory record",
                           .silent_n_sitrep),
                 if (.silent_n_neither > 0L)
                   sprintf("%d with neither a line-list confirmed case nor a sitrep top-up row in the window",
                           .silent_n_neither)),
               collapse = ";\n"))
}

pA5c <- ggplot(rank_prov, aes(prov_grp, rank_avg)) +
  geom_boxplot(data = dplyr::filter(rank_prov, outcome == "Not invaded"),
               width = 0.55, outlier.shape = NA, fill = "grey88", colour = MUTED,
               linewidth = 0.35) +
  geom_jitter(data = dplyr::filter(rank_prov, outcome == "Invaded within the window"),
              width = 0.13, height = 0, colour = FIT_RED, size = 1.5, shape = 21,
              fill = "white", stroke = 0.7) +
  geom_text(data = prov_n, aes(x = prov_grp, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.1, size = 1.95, colour = MUTED, family = base_family, lineheight = 1.05) +
  scale_y_continuous(transform = "log10", breaks = c(1, 3, 10, 30, 100, 300),
                     expand = expansion(mult = c(0.05, 0.22))) +
  scale_x_discrete(labels = prov_tick_lab) +
  labs(x = paste0("Province, ordered by its confirmed share of laboratory-resolved alerts as known at the origin (panel B; highest first)\n",
                  "grey boxes: at-risk zone-forecasts that were not invaded;  red rings: the realised invasions\n",
                  .silent_note),
       y = sprintf("Rank in the %d-week\nwatch-list", H_MAIN)) +
  theme_pub() + theme(axis.text.x = element_text(size = 5.9, lineheight = 0.95))

# D. The live watch-list against laboratory activity: where an operational list points.
live_h <- live %>% dplyr::filter(horizon == H_MAIN) %>%
  dplyr::mutate(rank_live = rank(-p_case_invasion, ties.method = "min"))
TOP_N <- 25L
# slice_min(with_ties = FALSE) would break a tie at the K boundary by row order, which is
# zone name after the joins. Assert there is no tie to break rather than hide one.
.p_sorted <- sort(live_h$p_case_invasion, decreasing = TRUE)
stopifnot(length(.p_sorted) > TOP_N, .p_sorted[TOP_N] > .p_sorted[TOP_N + 1L])
# The bars are coloured by whether the ZONE ITSELF has ever had an alert investigated to a
# laboratory outcome, which is the granularity an operational list is acted on; the province
# is in the tick label. Colouring by the PROVINCE's laboratory activity would say little about
# the zone: a province with laboratory activity can hold zones that have had none. (The log line
# below counts how many of the live top TOP_N sit in a province with any.) Both are measured at
# the analysis date, which is the live forecast's own origin.
live_top <- live_h %>% dplyr::slice_min(rank_live, n = TOP_N, with_ties = FALSE) %>%
  dplyr::left_join(dplyr::select(prov_pos, province, prov_ratio = ratio,
                                 prov_resolved = resolved), by = "province") %>%
  dplyr::mutate(lab_status = factor(dplyr::if_else(ever_investigated,
                                                   "Zone has had an alert laboratory-resolved",
                                                   "Zone has had none"),
                                    levels = c("Zone has had an alert laboratory-resolved",
                                               "Zone has had none")),
                zone_f = forcats::fct_reorder(sprintf("%s (%s)", health_zone, province),
                                              p_case_invasion))
message(sprintf("[ascert] live top-%d at h=%d: %d zone(s) have had an alert laboratory-resolved, %d have not; %d of %d sit in a province with any laboratory activity",
                TOP_N, H_MAIN, sum(live_top$ever_investigated), sum(!live_top$ever_investigated),
                sum(live_top$prov_resolved > 0), nrow(live_top)))

pA5d <- ggplot(live_top, aes(p_case_invasion, zone_f, fill = lab_status)) +
  geom_col(width = 0.7, alpha = 0.88) +
  scale_fill_manual(values = unname(c(OKABE[1], OKABE[2])), name = NULL, drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = sprintf("Live %d-week invasion probability (%s scale)", H_MAIN, PROB_SCALE),
       y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank(),
                      axis.text.y = element_text(size = 6.2),
                      legend.position = "top", legend.direction = "horizontal")

FigA5 <- (pA5a | pA5b) / pA5c / pA5d +
  patchwork::plot_layout(heights = c(1, 0.9, 1.15)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(pA5a, "FigureA5a_investigated_by_rank_band", 4.4, 3.4)
save_dual(pA5b, "FigureA5b_province_ratio_vs_rank",    4.6, 3.6)
save_dual(pA5c, "FigureA5c_rank_by_province",          8.4, 3.2)
save_dual(pA5d, "FigureA5d_live_watchlist_lab_status", 5.4, 4.0)
save_dual(FigA5, "FigureA5", 9.6, 11.0, dir = FIG_DIR)

# -----------------------------------------------------------------------------
# 9b. FIGURE A6 — does the frailty mark down the zones surveillance was not watching?
# -----------------------------------------------------------------------------
# The gamma-frailty form (TV_FORMS, 00_config.R) divides a still-uninvaded zone's hazard by
# (1 + theta*C), C the zone's own cumulative hazard, so it marks down zones exposed for weeks
# with no confirmed case. It cannot tell a zone that resists invasion from one invaded but not yet
# detected, and the back-test cannot either: its outcome is the first CONFIRMED case, so a
# silently invaded zone the frailty marks down is scored as a correct call. That is why the frailty
# is a secondary analysis and never featured. This check uses evidence from OUTSIDE the outcome.
# Were the zones it marks down, relative to the fixed-beta model it generalises (the stage-1
# pick), zones surveillance was already investigating (a discarded alert notified by the forecast
# origin: neg_prior, section 4)? Or zones nobody was looking at?
#   A  within bands of the fixed-beta model's rank, the share with a discarded alert by the origin
#      among the zone-forecasts the frailty moves DOWN vs UP;
#   B  the fixed-beta model's top-K watch-list, K = max(KS) (the largest list of Figure A4D): the
#      zone-forecasts the frailty pushes OUT of it vs those it brings IN.
# Ranks as in section 4: within (horizon, fold) on the shared support, ties averaged for the
# bands (rank_avg) and to the worst position for top-K membership (rank_max). On H_MAIN and the
# post-burn-in rounds (the rounds the selection and tv_comparison_* score). Intervals are
# zone-cluster bootstraps; a difference is taken on the same resamples (share_diff_boot()).
# Descriptive, not an adjustment: discarded alerts accumulate with time near the outbreak, and so
# does C, so the zones the frailty moves down have had longer to collect them.
message("\n[ascert] Figure A6 — the frailty form and surveillance effort")
.a6_files <- c(file.path(FIG_DIR, paste0("FigureA6.", c("pdf", "png"))),
               file.path(FIG_DIR, c("A6_frailty_by_rank_band.csv", "A6_frailty_topk_churn.csv",
                                    "A6_frailty_zone_forecasts.csv")))
A6 <- NULL
.a6_stage1 <- tryCatch(jsonlite::fromJSON(file.path(OUT, "key_outputs", "model_selection.json"))$featured$bayesian$stage1_method,
                       error = function(e) NULL)
.a6_sfx <- names(Filter(function(f) identical(f$form, "frailty"), get0("TV_FORMS", ifnotfound = list())))
if (length(.a6_stage1) == 1L && !is.na(.a6_stage1) && length(.a6_sfx) == 1L &&
    all(c(.a6_stage1, paste0(.a6_stage1, "-", .a6_sfx)) %in% lfo_all$method)) {
  A6_FIXED <- .a6_stage1; A6_FRAIL <- paste0(.a6_stage1, "-", .a6_sfx)
  A6_K <- max(KS)
  .a6_rows <- lfo_all %>%
    dplyr::filter(method %in% c(A6_FIXED, A6_FRAIL), horizon == H_MAIN,
                  !(as.logical(was_active_before) %in% TRUE), is.finite(p_invasion),
                  !invasion_is_burnin(cutoff))
  .cl <- .cells_by_h[[as.character(H_MAIN)]]
  if (!is.null(.cl))
    .a6_rows <- .a6_rows[paste(.a6_rows$fold_id, .a6_rows$health_zone, sep = "\r") %in% .cl, , drop = FALSE]
  # Each model ranked on its own rows; both must cover the same zone-forecasts.
  .a6_rank <- function(m) .a6_rows %>% dplyr::filter(method == m) %>%
    dplyr::group_by(fold_id) %>%
    dplyr::mutate(r_avg = rank(-p_invasion, ties.method = "average"),
                  r_max = rank(-p_invasion, ties.method = "max")) %>% dplyr::ungroup() %>%
    dplyr::select(health_zone, fold_id, cutoff, is_new_invasion, r_avg, r_max)
  .fx <- .a6_rank(A6_FIXED); .fr <- .a6_rank(A6_FRAIL)
  a6 <- dplyr::inner_join(.fx, dplyr::select(.fr, health_zone, fold_id, r_avg_frail = r_avg,
                                             r_max_frail = r_max), by = c("health_zone", "fold_id"))
  stopifnot(nrow(a6) == nrow(.fx), nrow(a6) == nrow(.fr), nrow(a6) > 0L,
            all(unique(a6$cutoff) %in% .cut))   # the ONE "before the origin" definition (.neg_by_cut)
  A6_BANDS <- c(0, 10, 25, 50, 100, Inf)
  A6_BAND_LAB <- c("1-10", "11-25", "26-50", "51-100", "101+")
  A6_DIR <- c("Moved down by the frailty", "Moved up by the frailty")
  a6 <- a6 %>%
    dplyr::left_join(.neg_by_cut, by = c("health_zone", "cutoff")) %>%
    dplyr::mutate(neg_prior = dplyr::coalesce(neg_prior, 0L), watched_prior = neg_prior > 0L,
                  d_rank = r_avg_frail - r_avg,
                  direction = dplyr::case_when(d_rank > 0 ~ A6_DIR[1], d_rank < 0 ~ A6_DIR[2],
                                               TRUE ~ "Unchanged"),
                  band = cut(r_avg, breaks = A6_BANDS, labels = A6_BAND_LAB, right = TRUE),
                  topk = dplyr::case_when(r_max <= A6_K & r_max_frail > A6_K ~ "out",
                                          r_max > A6_K & r_max_frail <= A6_K ~ "in",
                                          r_max <= A6_K & r_max_frail <= A6_K ~ "kept",
                                          TRUE ~ "neither"))
  A6_ROUNDS <- dplyr::n_distinct(a6$cutoff)

  # A. Share watched by fixed-beta rank band x direction (zone-forecasts with an unchanged rank
  #    are counted in the CSV, not plotted).
  .a6_bd <- a6 %>% dplyr::filter(direction %in% A6_DIR) %>%
    dplyr::mutate(grp = factor(paste(band, direction, sep = " | "),
                               levels = as.vector(t(outer(A6_BAND_LAB, A6_DIR, paste, sep = " | ")))))
  .bw6 <- share_by_group_boot(.a6_bd, "grp", "watched_prior")
  a6_band <- .bw6 %>%
    dplyr::mutate(band = factor(sub(" \\| .*$", "", group), levels = A6_BAND_LAB),
                  direction = factor(sub("^.* \\| ", "", group), levels = A6_DIR)) %>%
    dplyr::left_join(.a6_bd %>% dplyr::group_by(grp) %>%
                       dplyr::summarise(events = sum(is_new_invasion), mean_neg_prior = mean(neg_prior),
                                        mean_rank_shift = mean(d_rank), .groups = "drop") %>%
                       dplyr::transmute(group = as.character(grp), events, mean_neg_prior, mean_rank_shift),
                     by = "group") %>%
    dplyr::select(band, direction, n, n_zones, watched = flagged, share = est, lo, hi, interval,
                  events, mean_neg_prior, mean_rank_shift)
  a6_band_diff <- purrr::map_dfr(A6_BAND_LAB, function(b) {
    g <- paste(b, A6_DIR, sep = " | ")
    if (!all(g %in% .bw6$group[.bw6$n > 0])) return(NULL)
    dplyr::mutate(share_diff_boot(.bw6, g[1], g[2]), band = b, .before = 1)
  })
  if (!nrow(a6_band_diff))
    a6_band_diff <- tibble::tibble(band = character(0), est = numeric(0), lo = numeric(0), hi = numeric(0))
  a6_band <- a6_band %>%
    dplyr::left_join(dplyr::rename(a6_band_diff, diff_down_minus_up = est, diff_lo = lo, diff_hi = hi) %>%
                       dplyr::mutate(band = factor(band, levels = A6_BAND_LAB)), by = "band") %>%
    dplyr::bind_rows(a6 %>% dplyr::filter(direction == "Unchanged") %>%
                       dplyr::group_by(band) %>%
                       dplyr::summarise(n = dplyr::n(), n_zones = dplyr::n_distinct(health_zone),
                                        watched = sum(watched_prior), share = mean(watched_prior),
                                        events = sum(is_new_invasion), mean_neg_prior = mean(neg_prior),
                                        mean_rank_shift = 0, .groups = "drop") %>%
                       dplyr::mutate(direction = "Unchanged")) %>%
    dplyr::mutate(direction = factor(as.character(direction), levels = c(A6_DIR, "Unchanged"))) %>%
    dplyr::arrange(band, direction)

  # B. Churn of the fixed-beta top-K.
  A6_TK_LAB <- c(out = sprintf("Pushed out of the\nfixed-beta top %d", A6_K),
                 `in` = sprintf("Brought into the\nfixed-beta top %d", A6_K),
                 kept = sprintf("Kept in the\nfixed-beta top %d", A6_K))
  .a6_tk <- a6 %>% dplyr::filter(topk != "neither") %>%
    dplyr::mutate(grp = factor(topk, levels = names(A6_TK_LAB)))
  .bwk <- share_by_group_boot(.a6_tk, "grp", "watched_prior")
  a6_topk <- .bwk %>%
    dplyr::left_join(.a6_tk %>% dplyr::group_by(group = as.character(grp)) %>%
                       dplyr::summarise(events = sum(is_new_invasion), mean_neg_prior = mean(neg_prior),
                                        unwatched = sum(!watched_prior),
                                        unwatched_zones = dplyr::n_distinct(health_zone[!watched_prior]),
                                        .groups = "drop"), by = "group") %>%
    dplyr::mutate(label = unname(A6_TK_LAB[group]), top_k = A6_K) %>%
    dplyr::select(group, label, top_k, n, n_zones, watched = flagged, share = est, lo, hi, interval,
                  unwatched, unwatched_zones, events, mean_neg_prior)
  a6_topk_diff <- if (all(c("out", "in") %in% .bwk$group[.bwk$n > 0]))
    share_diff_boot(.bwk, "out", "in") else tibble::tibble(est = NA_real_, lo = NA_real_, hi = NA_real_)
  message(sprintf(paste0("[ascert] A6 (%s vs %s, h=%d, %d post-burn-in rounds): top-%d churn %d out / %d in; ",
                         "share with a discarded alert by the origin %.1f%% out vs %.1f%% in"),
                  A6_FRAIL, A6_FIXED, H_MAIN, A6_ROUNDS, A6_K,
                  a6_topk$n[a6_topk$group == "out"], a6_topk$n[a6_topk$group == "in"],
                  100 * a6_topk$share[a6_topk$group == "out"], 100 * a6_topk$share[a6_topk$group == "in"]))

  A6_COL <- stats::setNames(unname(c(OKABE[2], OKABE[1])), A6_DIR)
  .a6_pd <- a6_band %>% dplyr::filter(direction %in% A6_DIR) %>%
    dplyr::mutate(direction = factor(as.character(direction), levels = A6_DIR))
  .dodge <- ggplot2::position_dodge(width = 0.74)
  pA6a <- ggplot(.a6_pd, aes(band, share, fill = direction)) +
    geom_col(position = .dodge, width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = .dodge, width = 0, colour = INK,
                  linewidth = 0.4) +
    geom_text(aes(y = hi, label = scales::comma(n)), position = .dodge, vjust = -0.7, size = 2.0,
              colour = MUTED, family = base_family) +
    scale_fill_manual(values = A6_COL, name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.02, 0.16))) +
    labs(x = sprintf(paste0("Rank band in the fixed-beta %d-week watch-list\nnumbers above: zone-forecasts; ",
                            "bars: %.0f%% zone-cluster bootstrap"), H_MAIN, 100 * CONF_LEVEL),
         y = "Share with a discarded alert already\nnotified at the forecast origin") +
    theme_pub() + theme(legend.position = "top", legend.direction = "horizontal")
  .a6_kd <- a6_topk %>% dplyr::mutate(label = factor(label, levels = unname(A6_TK_LAB)))
  pA6b <- ggplot(.a6_kd, aes(label, share)) +
    geom_col(fill = PT_BLUE, alpha = 0.78, width = 0.62) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = INK, linewidth = 0.4) +
    geom_text(aes(y = hi, label = sprintf("%s (%d zones)", scales::comma(n), n_zones)), vjust = -0.7,
              size = 2.0, colour = MUTED, family = base_family) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.02, 0.16))) +
    labs(x = sprintf(paste0("Change the frailty makes to the fixed-beta top %d\nnumbers above: zone-forecasts; ",
                            "bars: %.0f%% zone-cluster bootstrap"), A6_K, 100 * CONF_LEVEL),
         y = "Share with a discarded alert already\nnotified at the forecast origin") +
    theme_pub()
  FigA6 <- (pA6a | pA6b) + patchwork::plot_layout(widths = c(1.45, 1)) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
  save_dual(pA6a, "FigureA6a_frailty_by_rank_band", 5.6, 3.6)
  save_dual(pA6b, "FigureA6b_frailty_topk_churn",   4.0, 3.6)
  save_dual(FigA6, "FigureA6", 10.0, 4.2, dir = FIG_DIR)
  A6 <- list(fixed = A6_FIXED, frail = A6_FRAIL, k = A6_K, rounds = A6_ROUNDS, band = a6_band,
             band_diff = a6_band_diff, topk = a6_topk, topk_diff = a6_topk_diff,
             n = nrow(a6), zones = dplyr::n_distinct(a6$health_zone))
} else {
  # No frailty back-test to compare (forms not cross-validated, or a stage-1 pick without forms):
  # remove the previous run's outputs so no stale Figure A6 survives beside this run's report.
  message("[ascert] Figure A6 skipped: the back-test holds no frailty form of the stage-1 model ",
          "(TV_FORMS not cross-validated this run); the previous run's A6 outputs are removed.")
  invisible(suppressWarnings(file.remove(.a6_files[file.exists(.a6_files)])))
}

# -----------------------------------------------------------------------------
# 10. WRITE EVERY NUMBER THAT APPEARS IN A PANEL
#     A figure is not a result until the values behind it can be quoted, so each
#     panel's underlying table is written out beside the figures.
# -----------------------------------------------------------------------------
message("\n[ascert] writing panel data")
.write_panel <- function(x, name) {
  readr::write_csv(x, file.path(FIG_DIR, name), na = "")
  message(sprintf("  wrote %-42s %d row(s)", name, nrow(x)))
}

.write_panel(zone_ascert %>%
     dplyr::select(health_zone, province, invaded, epi_travel_h, pop_count,
                   healthsite_count, healthsite_density, log_hs_density, log_pop,
                   surveillance_gap, healthcare_gap, social_vulnerability, ccvi,
                   alerts_total, n_confirmed, n_not_a_case, alerts_resolved,
                   confirm_ratio, confirm_lo, confirm_hi, resolution_frac,
                   neg_per_100k, log_neg_100k, resolved_per_100k, ever_investigated,
                   n_delay, delay_median, delay_p75, delay_over_cap, delay_censored,
                   first_onset, arrival_days, travel_time_h, cases_total),
   "zone_ascertainment.csv")
.write_panel(av_forest %>% dplyr::select(var, label, kind, set, slope_sd, lo, hi, p, dR2, n),
   "A2_arrival_added_variable.csv")
.write_panel(tibble::tibble(panel = "A2 (not plotted): added-variable fit, surveillance gap", stat = c("slope_per_sd", "ci_lo", "ci_hi", "p", "delta_R2",
                                          "spearman_rho", "rho_lo", "rho_hi", "rho_p", "n"),
                  value = c(fit_sg$slope_sd, fit_sg$lo, fit_sg$hi, fit_sg$p, fit_sg$dR2,
                            rho_sg$rho, rho_sg$lo, rho_sg$hi, rho_sg$p, fit_sg$n)),
   "A2_arrival_avp_surveillance.csv")
.write_panel(inv_t %>% dplyr::select(health_zone, province, arrival_days, surveillance_gap,
                           delay_median, n_delay, tercile = ter),
   "A2_arrival_by_zone.csv")
.write_panel(cd_forest %>% dplyr::arrange(var) %>%
               dplyr::select(var, label, delta, lo, hi, n_a, n_b),
   "A3_cliffs_delta.csv")
.write_panel(q_stat %>% dplyr::select(q_lab, delta, lo, hi, n_a, n_b),   # delta NA below Q_MIN_INV
   "A3_cliffs_delta_within_travel_quartile.csv")
.write_panel(or_forest %>% dplyr::select(var, label, model, or_sd, lo, hi, p, n, n_pos, ci_type),
   "A3_logistic_or.csv")
.write_panel(strat_tab, "A4_stratified_discrimination.csv")
.write_panel(rec_strat %>% dplyr::select(horizon, sg_stratum, k, caught, events, n_zones,
                                        recall, lo, hi),
   "A4_watchlist_recall_by_stratum.csv")
.write_panel(rho_rank, "A4_rank_vs_surveillance_zone_level.csv")
.write_panel(rho_ev %>% dplyr::select(h_lab, rho, lo, hi, p, n_zones = n, n_rows_total),
   "A4_rank_of_truth_vs_surveillance.csv")
.write_panel(dec_med %>% dplyr::rename(surveillance_gap_median = x, rank_median = y),
   "A4_median_rank_by_surveillance_decile.csv")
.write_panel(watch_band, "A5_investigation_by_rank_band.csv")
# *_end: the end-of-outbreak counts (Figure A1C); the rest: the as-known share and the ranks of
# Figure A5B, over the zone-forecasts issued after the province's first resolved alert (NA for
# a province that had none by any origin).
.write_panel(prov_join, "A5_province_positivity_and_rank.csv")
.write_panel(prov_origin, "A5_province_share_by_origin.csv")
.write_panel(prov_n %>% dplyr::rename(province_group = prov_grp), "A5_rank_by_province_counts.csv")
.write_panel(ev %>% dplyr::select(horizon, health_zone, surveillance_gap,
                                  mean_rank_of_truth = rank_avg, n_rows),
   "A4_rank_of_truth_by_zone.csv")
.write_panel(live_top %>% dplyr::select(health_zone, province, rank_live, p_case_invasion,
                              surveillance_gap, ever_investigated, neg_alerts,
                              prov_ratio, prov_resolved),
   "A5_live_top_watchlist.csv")
if (!is.null(A6)) {
  .write_panel(dplyr::mutate(A6$band, fixed_model = A6$fixed, frailty_model = A6$frail,
                             horizon = H_MAIN, .before = 1),
     "A6_frailty_by_rank_band.csv")
  .write_panel(dplyr::mutate(A6$topk, fixed_model = A6$fixed, frailty_model = A6$frail,
                             horizon = H_MAIN, diff_out_minus_in = A6$topk_diff$est,
                             diff_lo = A6$topk_diff$lo, diff_hi = A6$topk_diff$hi, .before = 1) %>%
                 dplyr::mutate(label = gsub("\n", " ", label)),
     "A6_frailty_topk_churn.csv")
  # Every zone-forecast behind both panels, so any number in them can be recounted.
  .write_panel(a6 %>% dplyr::transmute(health_zone, cutoff, origin = lfo_origin(cutoff), horizon = H_MAIN,
                                       rank_fixed = r_avg, rank_frailty = r_avg_frail,
                                       rank_fixed_max = r_max, rank_frailty_max = r_max_frail,
                                       direction, band, topk, neg_prior, watched_prior,
                                       is_new_invasion) %>%
                 dplyr::arrange(cutoff, rank_fixed, health_zone),
     "A6_frailty_zone_forecasts.csv")
}

.input_stamp <- function(paths) {
  paths <- paths[file.exists(paths)]
  stats::setNames(lapply(paths, function(f) list(
    mtime = format(file.info(f)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    bytes = as.numeric(file.info(f)$size),
    md5   = unname(tools::md5sum(f)))), basename(paths))
}
# The pipeline artifacts this module reads, and the line-list extract and its pointer. NOT
# stamped: the snapshot archive behind the earliest-sample / confirming-lab dates, the shapefile,
# the OSRM matrix and the covariate layers (the report says so).
INPUT_FILES <- c(fs_risk_csv(OUT),
                 file.path(OUT, "key_outputs", "arrival_predictors.csv"),
                 file.path(OUT, "key_outputs", "model_selection.json"),
                 file.path(OUT, "forecasts", "lfo_results.rds"),
                 LINELIST_JSON,
                 file.path(LINELIST_DIR, jsonlite::fromJSON(LINELIST_JSON)$folder,
                           "dhis2_processed_linelist.csv"))
.lfo_stamp <- attr(lfo_all, "lfo_stamp")

provenance <- list(
  module = "45_ascertainment_figures.R",
  inputs = .input_stamp(INPUT_FILES),
  lfo_stamp = if (is.null(.lfo_stamp)) NULL else
    .lfo_stamp[vapply(.lfo_stamp, function(v) length(v) == 1L && !is.list(v), logical(1))],
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  analysis_date = format(ANALYSIS_DATE), outbreak_start = format(OUTBREAK_START),
  linelist_folder = jsonlite::fromJSON(LINELIST_JSON)$folder,
  featured_model = FEATURED, risk_table_method = RISK_METHOD,
  # lfo_probability_column is the column the RANK panels use; Figure A4C's discrimination is
  # scored on a4c_prob_column (the raw probability, as the headline evaluation is), within each
  # of the a4c_rounds rounds (all of them).
  forecast_scale = FORECAST_SCALE, lfo_probability_column = PCOL,
  a4c_prob_column = A4C_PROB_COL, a4c_rounds = as.integer(A4C_N_ROUNDS),
  ascertainment_burnin_weeks = as.integer(ASCERTAINMENT_BURNIN_WEEKS),
  risk_prob_scale = PROB_SCALE,
  conf_level = CONF_LEVEL, n_boot = N_BOOT, random_seed = RANDOM_SEED,
  delay_cap_days = DELAY_CAP_DAYS,
  dropped_degenerate_pillars = DROPPED_PILLARS,
  n_zones = nrow(zone_ascert), n_invaded = sum(zone_ascert$invaded),
  n_zones_no_road_route = as.integer(.n_na_travel),
  n_synthetic_sitrep_rows_excluded = as.integer(n_synth),
  # Alert totals over the MAPPED zones (the report's totals); alerts whose zone name is not in
  # the spine are counted separately and are in no other field.
  n_alerts_in_window = n_alerts_window_mapped,
  n_alerts_placed_by_notification = n_alerts_unplaced_mapped,
  n_alerts_unmapped_zone = n_alerts_unmapped,
  unmapped_alert_zones = I(unmapped_alert_zones),
  n_resolved_alerts = sum(zone_ascert$alerts_resolved),
  n_zones_ever_investigated = sum(zone_ascert$ever_investigated),
  n_never_invaded_never_investigated =
    sum(!zone_ascert$ever_investigated & !zone_ascert$invaded),
  n_invaded_never_investigated = nrow(inv_no_lab),
  n_invaded_never_investigated_with_sitrep_rows = sum(inv_no_lab$n_sitrep_conf > 0L),
  n_provinces_with_lab_activity = as.integer(nrow(prov_pos) - N_PROV_SILENT),
  n_provinces_silent = as.integer(N_PROV_SILENT),
  surveillance_stratum_split = SG_SPLIT,
  lfo_folds = as.integer(dplyr::n_distinct(lfo_f$fold_id)),
  # at-risk zone-forecasts on the shared support, PER HORIZON (each zone-fold appears once per
  # horizon in lfo_f, so a count over both horizons would count it twice)
  lfo_zone_forecasts_by_horizon = as.list(stats::setNames(.n_by_h$n, paste0("h", .n_by_h$horizon))),
  lfo_horizon_main = as.integer(H_MAIN),
  # "watched prior" = a discarded alert notified on or before the forecast ORIGIN,
  # lfo_origin(cutoff) = cutoff + lfo_origin_offset_days; counted at the main horizon (A5A)
  lfo_zone_forecasts_h_main_watched_prior = sum(lfo_f$watched_prior[lfo_f$horizon == H_MAIN]),
  lfo_origin_offset_days = as.integer(LFO_ORIGIN_OFFSET_DAYS),
  live_top_n = as.integer(TOP_N),
  live_top_never_investigated = sum(!live_top$ever_investigated),
  # Figure A6 (NULL when the back-test held no frailty form): the secondary frailty form against
  # the fixed-beta model it generalises, at lfo_horizon_main, post-burn-in rounds.
  a6 = if (is.null(A6)) NULL else list(
    fixed_model = A6$fixed, frailty_model = A6$frail, top_k = as.integer(A6$k),
    rounds_post_burnin = as.integer(A6$rounds), zone_forecasts = as.integer(A6$n),
    zones = as.integer(A6$zones)))
writeLines(jsonlite::toJSON(provenance, auto_unbox = TRUE, pretty = TRUE, null = "null"),
           file.path(FIG_DIR, "ascertainment_run_info.json"))
message("  wrote ascertainment_run_info.json")

message("\n[ascert] done — figures and panel data in ", FIG_DIR)

# -----------------------------------------------------------------------------
# 11. REPORT — the headline numbers, composed from the objects above.
#     Every value is interpolated from a computed object; nothing here is typed
#     by hand, so the report cannot drift from the figures.
# -----------------------------------------------------------------------------
# EVERY DIRECTIONAL AND QUALITATIVE WORD BELOW IS DERIVED, not typed. An earlier version
# interpolated the numbers but hard-coded the adjectives ("LOWER", "HIGHER", "falls
# steeply"), and when the pipeline was re-run one of them became the opposite of what the
# data said while the sentence still read as generated output. A report that can state the
# wrong direction with the right number is worse than no report.
.dir <- function(a, b, hi = "higher", lo = "lower", eq = "the same") {
  if (!is.finite(a) || !is.finite(b)) return("not comparable")
  if (isTRUE(all.equal(a, b))) eq else if (a > b) hi else lo
}
# NA when either interval was not estimated (too few events or usable replicates): a missing
# interval is neither overlap nor separation, and the report says "not estimated".
.ci_overlap <- function(lo1, hi1, lo2, hi2)
  if (!all(is.finite(c(lo1, hi1, lo2, hi2)))) NA else lo1 <= hi2 && lo2 <= hi1
.overlap_word <- function(o) ifelse(is.na(o), "were not both estimated",
                                    ifelse(o, "overlap", "do not overlap"))
.rho_word <- function(rho, p, alpha = 0.05) {
  if (!is.finite(rho)) return("not estimated")
  if (!is.finite(p) || p >= alpha) return("not distinguishable from zero")
  if (rho > 0) "positive" else "negative"
}
.excl <- function(lo, hi, null_val) if (!all(is.finite(c(lo, hi)))) "not estimated" else
  if (lo > null_val || hi < null_val) "excludes" else "includes"

.rep_sg  <- rho_rank %>% dplyr::filter(horizon == H_MAIN)
.or_u    <- or_forest %>% dplyr::filter(var == "surveillance_gap", model == "Unadjusted")
.or_a    <- or_forest %>% dplyr::filter(var == "surveillance_gap", model != "Unadjusted")
.st      <- strat_tab %>% dplyr::arrange(horizon, sg_stratum)
.ter     <- inv_t %>% dplyr::group_by(ter) %>%
  dplyr::summarise(n = dplyr::n(), med = stats::median(arrival_days), .groups = "drop")
.band1   <- watch_band %>% dplyr::slice(1)
.bandlast<- watch_band %>% dplyr::slice(dplyr::n())
.q_last_inv <- max(za_q$qi[za_q$invaded])
.q_est   <- q_stat %>% dplyr::filter(is.finite(delta))
.cmp <- .st %>%
  dplyr::mutate(k = dplyr::if_else(grepl("^Sparser", as.character(sg_stratum)), "sp", "dn")) %>%
  dplyr::select(horizon, k, auc_pr_date, ap_lo, ap_hi, ap_chance, auc_roc_date, roc_lo, roc_hi,
                skill_chance, skill_ceiling, n_pos_zones) %>%
  tidyr::pivot_wider(names_from = k,
                     values_from = c(auc_pr_date, ap_lo, ap_hi, ap_chance, auc_roc_date, roc_lo,
                                     roc_hi, skill_chance, skill_ceiling, n_pos_zones)) %>%
  dplyr::arrange(horizon)
.cmp_sent <- function(est_sp, est_dn, lo_sp, hi_sp, lo_dn, hi_dn, d = 3, extra = NULL) {
  paste(sprintf(paste0("h=%d %s in the sparser stratum (%.", d, "f vs %.", d, "f%s; the two intervals %s)"),
                .cmp$horizon,
                mapply(.dir, est_sp, est_dn),
                est_sp, est_dn,
                if (is.null(extra)) "" else extra,
                .overlap_word(mapply(.ci_overlap, lo_sp, hi_sp, lo_dn, hi_dn))),
        collapse = "; ")
}
.ov_all <- c(mapply(.ci_overlap, .cmp$ap_lo_sp, .cmp$ap_hi_sp, .cmp$ap_lo_dn, .cmp$ap_hi_dn),
             mapply(.ci_overlap, .cmp$roc_lo_sp, .cmp$roc_hi_sp, .cmp$roc_lo_dn, .cmp$roc_hi_dn))
.any_separation <- any(!.ov_all, na.rm = TRUE)
.n_ov_missing   <- sum(is.na(.ov_all))
# Section 8's summary of Figure A3D, derived from the two odds-ratio intervals.
.or_u_x <- .excl(.or_u$lo, .or_u$hi, 1); .or_a_x <- .excl(.or_a$lo, .or_a$hi, 1)
.or_story <- if ("not estimated" %in% c(.or_u_x, .or_a_x)) {
  "the odds ratio of invasion per SD of surveillance gap could not be estimated both unadjusted and adjusted for epicentre travel time"
} else if (.or_u_x == "excludes" && .or_a_x == "includes") {
  "invaded and never-invaded zones differ in surveillance gap unadjusted (odds-ratio interval excludes 1) but not detectably once epicentre travel time is adjusted for (interval includes 1), which is what a difference carried by proximity to the epicentre would produce"
} else if (.or_u_x == "excludes") {
  "invaded and never-invaded zones differ in surveillance gap both unadjusted and after adjustment for epicentre travel time (both odds-ratio intervals exclude 1)"
} else if (.or_a_x == "excludes") {
  "invaded and never-invaded zones do not detectably differ in surveillance gap unadjusted (odds-ratio interval includes 1), but do once epicentre travel time is adjusted for (interval excludes 1)"
} else {
  "invaded and never-invaded zones do not detectably differ in surveillance gap, unadjusted or adjusted for epicentre travel time (both odds-ratio intervals include 1)"
}
# Section 4's quartile sentence: every quartile, estimated or not.
.q_words <- paste(sprintf("%s: %s", sub(":.*$", "", as.character(q_stat$q_lab)),
                          ifelse(is.finite(q_stat$delta),
                                 sprintf("delta %s, which %s zero (n = %d invaded vs %d never invaded)",
                                         fmt_ci(q_stat$delta, q_stat$lo, q_stat$hi),
                                         mapply(.excl, q_stat$lo, q_stat$hi, 0), q_stat$n_a, q_stat$n_b),
                                 sprintf("not estimated (%d invaded zone%s)", q_stat$n_a,
                                         ifelse(q_stat$n_a == 1L, "", "s")))),
                  collapse = "; ")
.q_sign_pos <- any(.q_est$lo > 0); .q_sign_neg <- any(.q_est$hi < 0)

# Section 7 (Figure A6), composed here so every word in it is derived from A6's intervals.
.pct_ci <- function(e, l, h) sprintf("%.1f%% [%.1f, %.1f]", 100 * e, 100 * l, 100 * h)
.pts_ci <- function(e, l, h) sprintf("%+.1f points [%+.1f, %+.1f]", 100 * e, 100 * l, 100 * h)
.a6_report <- if (is.null(A6)) {
  "- Not produced: this run's back-test holds no frailty form of the stage-1 model."
} else local({
  b <- A6$band; k <- A6$topk; kd <- A6$topk_diff
  .row <- function(t, ...) { r <- dplyr::filter(t, ...); if (nrow(r) == 1L) r else NULL }
  bands <- purrr::map_chr(seq_len(nrow(A6$band_diff)), function(i) {
    d <- A6$band_diff[i, ]
    dn <- .row(b, band == d$band, direction == A6_DIR[1]); up <- .row(b, band == d$band, direction == A6_DIR[2])
    sprintf(paste0("  - rank band %s: %s of the %s zone-forecasts moved down against %s of the %s ",
                   "moved up; the difference, %s, %s zero."),
            d$band, .pct_ci(dn$share, dn$lo, dn$hi), scales::comma(dn$n),
            .pct_ci(up$share, up$lo, up$hi), scales::comma(up$n),
            .pts_ci(d$est, d$lo, d$hi), .excl(d$lo, d$hi, 0))
  })
  ko <- .row(k, group == "out"); ki <- .row(k, group == "in")
  verdict <- if (is.null(ko) || is.null(ki) || !all(is.finite(c(kd$lo, kd$hi)))) {
    "The top-K contrast could not be estimated, so the check neither clears nor convicts the frailty."
  } else if (kd$lo > 0) {
    sprintf(paste0("Over this back-test the zones the frailty pushes out of the watch-list are MORE often ",
                   "ones surveillance had already investigated than the zones it brings in, so its ",
                   "re-ranking is not mainly a demotion of unwatched zones. That does not clear the %s ",
                   "unwatched zone-forecasts it did push out."), scales::comma(ko$unwatched))
  } else if (kd$hi < 0) {
    paste0("The zones the frailty pushes out of the watch-list are LESS often ones surveillance had ",
           "investigated than the zones it brings in: the confound between resisting invasion and ",
           "escaping detection is visible in its re-ranking.")
  } else {
    paste0("The zones the frailty pushes out of and brings into the watch-list do not detectably ",
           "differ in prior surveillance at this sample size, so the check neither clears nor convicts it.")
  }
  c(sprintf(paste0("The gamma-frailty form `%s` is a secondary analysis and is never featured (TV_FORMS, ",
                   "00_config.R). It divides a still-uninvaded zone's hazard by (1 + theta C), C the zone's ",
                   "own cumulative hazard, so it cannot tell a zone that resists invasion from one invaded ",
                   "but not yet detected. Neither can the back-test: its outcome is the first confirmed ",
                   "case, so a silently invaded zone the frailty marks down is scored as a correct call. ",
                   "This section checks the frailty against evidence from outside the outcome, discarded ",
                   "alerts notified by the forecast origin, relative to the fixed-beta model it generalises ",
                   "(`%s`), at h=%d on the %d post-burn-in rounds (%s zone-forecasts, %d zones)."),
            A6$frail, A6$fixed, H_MAIN, A6$rounds, scales::comma(A6$n), A6$zones),
    "",
    "- Share with a discarded alert by the origin, by the fixed-beta model's rank band (Figure A6A):",
    bands,
    if (!is.null(ko) && !is.null(ki))
      sprintf(paste0("- The frailty pushed %s zone-forecasts (%d zones) out of the fixed-beta top %d and ",
                     "brought %s (%d zones) in (Figure A6B). %s of those pushed out had a discarded alert ",
                     "by the origin, against %s of those brought in; the difference, %s, %s zero. %s of ",
                     "those pushed out (%d zones) had none: exposed zones nobody had investigated, which ",
                     "the frailty marked down. %d of the pushed-out and %d of the brought-in ",
                     "zone-forecasts were followed by a confirmed invasion within the window."),
              scales::comma(ko$n), ko$n_zones, A6$k, scales::comma(ki$n), ki$n_zones,
              .pct_ci(ko$share, ko$lo, ko$hi), .pct_ci(ki$share, ki$lo, ki$hi),
              .pts_ci(kd$est, kd$lo, kd$hi), .excl(kd$lo, kd$hi, 0),
              scales::comma(ko$unwatched), ko$unwatched_zones, ko$events, ki$events)
    else sprintf("- The frailty changed no zone's membership of the fixed-beta top %d.", A6$k),
    paste0("- ", verdict),
    paste0("- Descriptive, not an adjustment: discarded alerts accumulate with time near the outbreak, ",
           "and so does the frailty's cumulative exposure, so the zones it moves down have had longer ",
           "to collect them. A discarded alert notified by the origin means an alert was being ",
           "investigated by then, not that it had been resolved."))
})

report <- c(
"# Detection versus onset: spatial ascertainment and the invasion forecast",
"",
sprintf("Generated %s from the %s line-list snapshot, analysis date %s. Featured model `%s`; LFO probabilities on the `%s` column for the rank panels (a rank within a round does not depend on the scale). The discrimination of Figure A4C is scored as the headline performance is: within each forecast round, over all %d rounds, on the raw `p_invasion_raw` column (a within-round score is the same on the recalibrated one). All intervals %.0f%%. The pipeline artifacts this report reads, and the line-list extract with its `latest.json` pointer, are recorded with modification times and MD5 sums in `ascertainment_run_info.json`; the snapshot archive behind the earliest sample and confirming-test laboratory dates, the shapefile, the OSRM matrix and the covariate layers are not stamped. The pipeline artifacts are rewritten by `run_all.R` and `43_spread_kinematics.R`, so a number quoted from here must be reconciled against that stamp.",
        format(Sys.Date()), jsonlite::fromJSON(LINELIST_JSON)$folder, format(ANALYSIS_DATE),
        FEATURED, PCOL, A4C_N_ROUNDS, 100 * CONF_LEVEL),
"",
"## 1. The estimand",
"",
sprintf("The modelled outcome is the onset week of a zone's first **confirmed** case. Of %s alerts in the outbreak window attributed to one of the %d health zones (%s line-list records, plus %s with no onset or specimen date inside the outbreak, placed by their notification date), %s reached a laboratory outcome (%s confirmed, %s discarded); %s did not reach one (still suspected or probable, or never classified, including alerts sampled without a recorded outcome). %s Confirmation ratios below are computed on the resolved alerts only, excluding the %d synthetic sitrep top-up rows, which are a count reconciliation against the official cumulative and carry no test.",
        scales::comma(n_alerts_window_mapped), nrow(zone_ascert),
        scales::comma(n_alerts_ll_mapped),
        scales::comma(n_alerts_unplaced_mapped),
        scales::comma(sum(zone_ascert$alerts_resolved)),
        scales::comma(sum(zone_ascert$n_confirmed)),
        scales::comma(sum(zone_ascert$n_not_a_case)),
        scales::comma(n_alerts_window_mapped - sum(zone_ascert$alerts_resolved)),
        if (n_alerts_unmapped > 0L)
          sprintf("A further %s alert%s in the window %s a zone name that matches none of the %d zones (%s) and %s excluded from every count here.",
                  scales::comma(n_alerts_unmapped), if (n_alerts_unmapped == 1L) "" else "s",
                  if (n_alerts_unmapped == 1L) "carries" else "carry",
                  nrow(zone_ascert), paste(sprintf("\"%s\"", unmapped_alert_zones), collapse = ", "),
                  if (n_alerts_unmapped == 1L) "is" else "are")
        else "Every alert in the window is attributed to one of these zones.",
        n_synth),
"",
"## 2. Ascertainment in space (Figure A1)",
"",
# NOT "in those zones an invasion could not have been detected": invaded zones confirmed only
# through sitrep top-up rows (no test, so no resolved alert) are among them. Their number is
# stated, and the sitrep explanation is given only where inv_no_lab (3.5) confirms it.
sprintf("- %d of %d health zones never had a single laboratory-resolved alert in the line list during the outbreak: %d of the %d never-invaded zones and %d invaded zone%s. %s",
        sum(!zone_ascert$ever_investigated), nrow(zone_ascert),
        sum(!zone_ascert$ever_investigated & !zone_ascert$invaded), sum(!zone_ascert$invaded),
        sum(!zone_ascert$ever_investigated & zone_ascert$invaded),
        if (sum(!zone_ascert$ever_investigated & zone_ascert$invaded) == 1L) "" else "s",
        if (!nrow(inv_no_lab))
          "Every invaded zone had at least one."
        else if (INV_NO_LAB_SITREP_ONLY)
          sprintf("The invaded one%s (%s) %s detected only through sitrep top-up rows, which carry no laboratory record, not through a laboratory-resolved alert.",
                  if (nrow(inv_no_lab) == 1L) "" else "s", paste(inv_no_lab$health_zone, collapse = ", "),
                  if (nrow(inv_no_lab) == 1L) "was" else "were")
        else
          sprintf("Of the invaded ones, %d (%s) %s neither a laboratory-resolved alert nor a sitrep top-up row in the window, which these data do not explain; the other %d %s detected only through sitrep top-up rows, which carry no laboratory record.",
                  sum(inv_no_lab$n_sitrep_conf == 0L),
                  paste(inv_no_lab$health_zone[inv_no_lab$n_sitrep_conf == 0L], collapse = ", "),
                  if (sum(inv_no_lab$n_sitrep_conf == 0L) == 1L) "has" else "have",
                  sum(inv_no_lab$n_sitrep_conf > 0L),
                  if (sum(inv_no_lab$n_sitrep_conf > 0L) == 1L) "was" else "were")),
sprintf("- %d of %d provinces resolved no alert at all. Among the %d that did, the confirmed share of resolved alerts runs from %.0f%% (%s, %d/%d) to %.0f%% (%s, %d/%d); the denominators matter more than the ordering. (Figure A1C's Wilson intervals treat every alert as independent.)",
        N_PROV_SILENT, nrow(prov_pos), nrow(prov_active),
        100 * min(prov_active$ratio), prov_active$province[which.min(prov_active$ratio)],
        prov_active$confirmed[which.min(prov_active$ratio)], prov_active$resolved[which.min(prov_active$ratio)],
        100 * max(prov_active$ratio), prov_active$province[which.max(prov_active$ratio)],
        prov_active$confirmed[which.max(prov_active$ratio)], prov_active$resolved[which.max(prov_active$ratio)]),
"",
"## 3. Surveillance density and the timing of the first confirmed case (Figure A2)",
"",
sprintf("- Added to a baseline of road travel time from the epicentre region, the surveillance gap changes arrival time by %s days per SD (%s), changing R-squared by %+.3f over %d invaded zones. Not plotted: the added-variable panels were dropped from Figure A2; the fits are in `A2_arrival_avp_surveillance.csv` and `A2_arrival_added_variable.csv`.",
        fmt_ci(fit_sg$slope_sd, fit_sg$lo, fit_sg$hi, 1), fmt_p(fit_sg$p), fit_sg$dR2, fit_sg$n),
sprintf("- Median arrival time by surveillance-gap tercile: %s.",
        paste(sprintf("%s = %.0f d (n = %d)", gsub("\n", " ", .ter$ter), .ter$med, .ter$n), collapse = "; ")),
sprintf("- The zone's own median onset-to-confirmation delay (confirmation dated by the confirming test's laboratory date), as an ordering of arrival time, is %s (Spearman rho %s, %s, n = %d zones)%s",
        .rho_word(rho_del$rho, rho_del$p),
        fmt_ci(rho_del$rho, rho_del$lo, rho_del$hi), fmt_p(rho_del$p), nrow(del),
        if (identical(.rho_word(rho_del$rho, rho_del$p), "not distinguishable from zero"))
          " — on these data the arrival ordering above is therefore not explained by a zone-level reporting delay."
        else " — this is NOT null, so the arrival ordering above must not be presented as free of reporting delay."),
"- CAVEAT, stated as a limit of the design: facility density also proxies remoteness and urbanicity, so this is consistent with ascertainment delay but is not identified against genuinely later epidemiological arrival.",
"",
"## 4. Comparability of invaded and never-invaded zones (Figure A3)",
"",
sprintf("- Unadjusted, the surveillance gap of invaded zones is %s than that of never-invaded zones (Cliff's delta %s, which %s zero; Mann-Whitney %s).",
        .dir(cd_sg$delta, 0, "higher", "lower", "no different"),
        fmt_ci(cd_sg$delta, cd_sg$lo, cd_sg$hi), .excl(cd_sg$lo, cd_sg$hi, 0), fmt_p(mw_sg$p.value)),
sprintf("- Invaded zones by quartile of epicentre travel time: %s%s. Within each quartile (delta estimated where it holds at least %d invaded zones): %s.%s",
        paste(sprintf("%s %d", sub(":.*$", "", as.character(q_stat$q_lab)), q_stat$n_a), collapse = ", "),
        if (.q_last_inv < TT_Q)
          sprintf(", so no invaded zone lies beyond quartile %d of %d and the contrast cannot be made beyond it", .q_last_inv, TT_Q)
        else "",
        Q_MIN_INV, .q_words,
        if (.q_sign_pos && .q_sign_neg)
          " The quartiles whose interval excludes zero disagree in sign, so the unadjusted contrast has no single direction across travel time."
        else ""),
sprintf("- Odds ratio for being invaded, per SD of surveillance gap: %s unadjusted (interval %s 1), %s adjusted for epicentre travel time (interval %s 1).",
        fmt_ci(.or_u$or_sd, .or_u$lo, .or_u$hi), .or_u_x,
        fmt_ci(.or_a$or_sd, .or_a$lo, .or_a$hi), .or_a_x),
"",
"## 5. The forecast and the surveillance gradient (Figure A4)",
"",
sprintf("- The association between a zone's surveillance gap and its mean rank in the %d-week watch-list is %s (zone-level Spearman rho %s, %s, n = %d zones).%s",
        H_MAIN, .rho_word(.rep_sg$rho, .rep_sg$p),
        fmt_ci(.rep_sg$rho, .rep_sg$lo, .rep_sg$hi), fmt_p(.rep_sg$p), .rep_sg$n,
        if (identical(.rho_word(.rep_sg$rho, .rep_sg$p), "positive"))
          " A larger rank number is a lower place on the list, so zones with sparser facilities are ranked lower."
        else if (identical(.rho_word(.rep_sg$rho, .rep_sg$p), "negative"))
          " A larger rank number is a lower place on the list, so zones with sparser facilities are ranked higher."
        else ""),
sprintf("- Among the zones that were invaded, the association between the surveillance gap and the rank of the truth is: %s.",
        paste(sprintf("%s %s (rho %s, %s, %d zones)", rho_ev$h_lab,
                      mapply(.rho_word, rho_ev$rho, rho_ev$p),
                      fmt_ci(rho_ev$rho, rho_ev$lo, rho_ev$hi), fmt_p(rho_ev$p), rho_ev$n),
              collapse = "; ")),
sprintf("- Within-stratum discrimination (Figure A4C) is scored as the headline performance is: within each forecast round, over all %d rounds. AUC-PR skill is average precision over the base rate, but neither of its reference values is common to the strata. A uniformly random ordering does not score 1: its expected skill (the chance level) depends on each round's size and invasion count. The attainable maximum is the inverse of the event-weighted base rate. %s. The skill numbers must therefore not be read as a ratio of discrimination.",
        A4C_N_ROUNDS,
        paste0("At ", paste(sprintf("h=%d the chance level is %.2f in the sparser stratum vs %.2f in the denser, and the maximum is %s there (%.0f vs %.0f)",
                                    .cmp$horizon, .cmp$skill_chance_sp, .cmp$skill_chance_dn,
                                    mapply(.dir, .cmp$skill_ceiling_sp, .cmp$skill_ceiling_dn),
                                    .cmp$skill_ceiling_sp, .cmp$skill_ceiling_dn),
                            collapse = "; at "))),
sprintf("  - AUC-PR skill (interval in brackets; chance level and attainable maximum in parentheses): %s.",
        paste(sprintf("h=%d %s %s (chance %.2f, max %.0f)", .st$horizon, .st$sg_stratum,
                      fmt_ci(.st$auc_pr_skill_date, .st$skill_lo, .st$skill_hi, 1),
                      .st$skill_chance, .st$skill_ceiling), collapse = "; ")),
sprintf("  - Average precision, whose ceiling (1) is common to the strata but whose chance level is not: %s.",
        .cmp_sent(.cmp$auc_pr_date_sp, .cmp$auc_pr_date_dn, .cmp$ap_lo_sp, .cmp$ap_hi_sp,
                  .cmp$ap_lo_dn, .cmp$ap_hi_dn, 3,
                  extra = sprintf("; chance %.3f vs %.3f", .cmp$ap_chance_sp, .cmp$ap_chance_dn))),
sprintf("  - Within-round AUC-ROC, which shares a [0, 1] scale and a 0.5 chance level across strata: %s.",
        .cmp_sent(.cmp$auc_roc_date_sp, .cmp$auc_roc_date_dn, .cmp$roc_lo_sp, .cmp$roc_hi_sp,
                  .cmp$roc_lo_dn, .cmp$roc_hi_dn, 3)),
sprintf("  - Distinct invaded zones over all rounds (the support of every comparison above): %s. %s",
        paste(sprintf("h=%d sparser %d, denser %d", .cmp$horizon, .cmp$n_pos_zones_sp,
                      .cmp$n_pos_zones_dn), collapse = "; "),
        if (.any_separation)
          "At least one comparison of average precision or AUC-ROC separates the strata; read the per-horizon intervals above before drawing a conclusion."
        else if (.n_ov_missing == length(.ov_all))
          "No interval comparison could be made: the intervals were not estimated."
        else
          sprintf("No comparison of average precision or AUC-ROC separates the strata at this event count: every %sinterval overlaps%s, so the defensible statement is that discrimination is not detectably degraded where surveillance is thin — not that it differs in either direction.",
                  if (.n_ov_missing > 0L) "estimated " else "",
                  if (.n_ov_missing > 0L) sprintf(" (%d comparison%s not estimated)", .n_ov_missing,
                                                  if (.n_ov_missing == 1L) "" else "s") else "")),
"",
"## 6. Where the surveillance system was looking (Figure A5)",
"",
sprintf("- Of %s at-risk zone-forecasts in the %d-week watch-list, %s (%.1f%%) were for a zone that already had a discarded alert notified at the forecast origin.",
        scales::comma(sum(lfo_f$horizon == H_MAIN)), H_MAIN,
        scales::comma(sum(lfo_f$watched_prior[lfo_f$horizon == H_MAIN])),
        100 * mean(lfo_f$watched_prior[lfo_f$horizon == H_MAIN])),
sprintf("- Moving down the watch-list, that share is %s at the bottom than at the top: %.1f%% [%.1f, %.1f] in rank band %s against %.1f%% [%.1f, %.1f] in band %s.",
        .dir(.bandlast$share, .band1$share, "higher", "lower", "unchanged"),
        100 * .band1$share, 100 * .band1$lo, 100 * .band1$hi, .band1$band,
        100 * .bandlast$share, 100 * .bandlast$lo, 100 * .bandlast$hi, .bandlast$band),
sprintf("- The association between the province confirmed share of laboratory-resolved alerts, as known at each forecast origin, and the mean rank of the province's at-risk zone-forecasts is %s (Spearman rho %s, %s, %d provinces; Figure A5B).",
        .rho_word(rho_prov$rho, rho_prov$p), fmt_ci(rho_prov$rho, rho_prov$lo, rho_prov$hi),
        fmt_p(rho_prov$p), rho_prov$n),
sprintf("- In the live %d-week watch-list, %d of the top %d zones have never had a single alert laboratory-resolved.",
        H_MAIN, sum(!live_top$ever_investigated), TOP_N),
"",
"## 7. The frailty form and surveillance effort (Figure A6)",
"",
.a6_report,
"",
"## 8. What these figures do and do not establish",
"",
sprintf("- They describe how laboratory ascertainment varies in space (section 2), and show that %s.", .or_story),
"- They do NOT identify ascertainment delay against genuinely later epidemiological arrival: no measurement of unobserved transmission exists in these data, and every proxy used here is also a proxy for remoteness.",
"- Quantities generated BY the outbreak (confirmation ratio, alert counts, confirmation delay) are never used in an invaded-versus-not comparison, where they would be mechanically tied to the outcome. Figure A2 uses them only among already-invaded zones. Figure A5 re-measures them at each forecast origin (the discarded alerts notified, and the province confirmed share of the alerts resolved, by that origin), and the live list of Figure A5D at its own origin, the analysis date. Figure A1 maps their end-of-outbreak values descriptively; those are not the shares Figure A5B/C rank against.",
sprintf("- In Figure A5C the final group pools the zone-forecasts issued while their province had no laboratory-resolved alert known at the origin; %s.",
        gsub("\n", " ", .silent_note)))

writeLines(report, file.path(FIG_DIR, "ASCERTAINMENT_REPORT.md"))
message("  wrote ASCERTAINMENT_REPORT.md")

}   # end of the script-run guard

