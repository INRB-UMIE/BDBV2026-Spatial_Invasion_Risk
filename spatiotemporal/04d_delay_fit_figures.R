# =============================================================================
# 04d_delay_fit_figures.R — delay-distribution fit diagnostics (EpiDist vs data)
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# WHY THIS FILE EXISTS
# --------------------
# run_all.R step 1a fits the DHIS2 reporting delays (interval-censored MLE, and —
# when RUN_EPIDIST is on — the Bayesian EpiDist naive + marginal models) and writes
# the selected parameters to CSV. Nothing in the pipeline ever LOOKED at those fits:
# the naive-vs-censored QA figure lives inside 04c's Rscript-gated MAIN block, which
# run_all.R never triggers because it source()s the file. The onset->sample delay is
# routed into the onset imputation (01_data_prep.R), the nowcast rate (04/04b) and the
# R(t) truncation, so a silently bad fit propagates into every downstream number with
# no visual check. This module renders that check on EVERY run.
#
# WHAT IT DRAWS (all on the SAME windowed records the fits used)
#   A. dhis2_delay_epidist_fits_<date>       2x2: empirical histogram + fitted densities
#                                            (EpiDist naive & marginal w/ 95% credible
#                                            band, interval-censored MLE) per delay.
#   B. dhis2_delay_fit_cdf_<date>            2x2: empirical CDF vs fitted CDFs — the
#                                            sharper goodness-of-fit view (no binwidth
#                                            artefacts), with the max |ECDF-CDF| gap.
#   C. dhis2_delay_mean_estimators_<date>    forest: mean delay by estimator with 95%
#                                            CrI, empirical mean, and the lab-reference
#                                            delay the pipeline used to assume.
#   D. dhis2_onset_sample_fit_<date>         single large panel for the ONE delay that
#                                            drives the pipeline (density + CDF inset).
#   E. dhis2_onset_sample_by_class_<date>    onset->sample PER CASE CLASSIFICATION: fitted
#                                            completeness curves for confirmed /
#                                            test-negative / suspected, plus a mean-with-CrI
#                                            forest against the pooled fit. Read from the
#                                            published parameters via the pipeline's own
#                                            loader -- this panel fits nothing.
#   +  dhis2_delay_fit_summary.csv           one row per delay x estimator: parameters,
#                                            mean/sd, CrI, fit statistic, and a `plotted`
#                                            flag. Includes the EpiDist families that were
#                                            NOT plotted, so a diverged runner-up stays on
#                                            the record instead of vanishing behind the
#                                            gamma-preferred pick.
#
# READING THE FIGURE — the marginal SHOULD sit right of the histogram
# ------------------------------------------------------------------
# The empirical histogram is right-truncated: cases whose onset is recent can only
# appear in the data if their delay was SHORT enough to be sampled before the extract,
# so the observed delays are biased short. The EpiDist MARGINAL model corrects for that
# (and for double interval censoring), so its density is EXPECTED to sit to the right of
# the bars. That gap is the correction working, not a misfit. The NAIVE curve is drawn
# precisely so the two can be compared: naive ≈ histogram, marginal shifted right by the
# truncation correction. The censored MLE mitigates truncation only by dropping the last
# TEST_DAYS of onsets, so it normally lands between the two.
#
# COST. Free by default: it reuses run_all.R's already-computed onset->sample EpiDist
# posterior (via 04c's draw registry — no Stan refit) and only re-runs the cheap
# interval-censored MLE (~seconds) for the other three delays. Set
# DELAY_FIG_EPIDIST_ALL=TRUE to also fit EpiDist for those three (12 extra Stan fits,
# several minutes) when the descriptive delays need the same treatment.
#
# Every entry point is guarded and non-fatal: a missing package, a failed fit or an
# absent date column degrades the figure (fewer curves / fewer panels) but never
# breaks the run.
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

# 04c supplies build_dhis2_delay_populations(), .fit_all_censored(), .implied_mean(),
# epidist_draws_get() and friends. run_all.R sources it first; a standalone/interactive
# use of this file sources it here.
if (!exists("build_dhis2_delay_populations", mode = "function"))
  source(file.path(here::here(), "spatiotemporal", "04c_dhis2_delay_windows.R"))

.DFG_HAVE_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)

# ---- Design system (matches 24_review_figures.R / make_manuscript_figures.R) --
.DFG_INK <- "grey15"; .DFG_MUTED <- "grey38"; .DFG_GRID <- "grey92"; .DFG_FAM <- "sans"
.DFG_BAR <- "grey86"; .DFG_BAR_EDGE <- "grey62"
# Okabe-Ito. One colour per ESTIMATOR, fixed across all panels and figures so the eye
# can carry a curve from one delay to the next.
.DFG_COL <- c(
  "EpiDist marginal (truncation-corrected)" = "#0072B2",  # blue  — the pipeline default
  "EpiDist naive (no truncation correction)" = "#E69F00", # amber — the comparison
  "Interval-censored MLE"                    = "#009E73",  # green — the fallback
  "Empirical mean"                           = "#CC79A7",  # pink
  "Lab-linelist reference"                   = "#999999")
.DFG_LTY <- c(
  "EpiDist marginal (truncation-corrected)"  = "solid",
  "EpiDist naive (no truncation correction)" = "22",
  "Interval-censored MLE"                    = "42")

.dfg_theme <- function(base = 10.5) {
  ggplot2::theme_minimal(base_size = base, base_family = .DFG_FAM) %+replace% ggplot2::theme(
    plot.title    = ggplot2::element_text(size = base + 0.5, colour = .DFG_INK,
                                          face = "bold", hjust = 0,
                                          margin = ggplot2::margin(b = 2)),
    plot.subtitle = ggplot2::element_text(size = base - 1.8, colour = .DFG_MUTED, hjust = 0,
                                          margin = ggplot2::margin(b = 5)),
    axis.title    = ggplot2::element_text(size = base - 0.8, colour = .DFG_MUTED),
    axis.title.x  = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
    axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(r = 4), angle = 90),
    axis.text     = ggplot2::element_text(size = base - 1.6, colour = .DFG_MUTED),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = .DFG_GRID, linewidth = 0.3),
    legend.position = "top", legend.justification = "left",
    legend.title  = ggplot2::element_blank(),
    legend.text   = ggplot2::element_text(size = base - 2, colour = .DFG_INK),
    legend.key.height = ggplot2::unit(9, "pt"), legend.key.width = ggplot2::unit(18, "pt"),
    plot.margin   = ggplot2::margin(7, 9, 7, 7))
}

# cairo_pdf where it actually works (see .dfg_probe_cairo below), else the base device with
# ASCII-downgraded text so the PDF and the PNG say the same thing.
.dfg_pdf_device <- function() if (.DFG_UNICODE) grDevices::cairo_pdf else "pdf"

.dfg_save <- function(p, dir, stem, w, h, verbose = TRUE) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  # Gate on the FULL destination path, not the bare stem: FIGURE_DROP entries are
  # "<directory>/<stem>" and the raw/ exclusion inspects path components, neither of
  # which a basename can match.
  path <- file.path(dir, stem)
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !.fk(path)) return(invisible(p))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ok <- tryCatch({
    ggplot2::ggsave(paste0(path, ".pdf"), p, width = w, height = h,
                    device = .dfg_pdf_device(), bg = "white")
    TRUE
  }, error = function(e) { message("  [delay-fig] PDF save failed: ", conditionMessage(e)); FALSE })
  tryCatch(ggplot2::ggsave(paste0(path, ".png"), p, width = w, height = h, dpi = 400, bg = "white"),
           error = function(e) NULL)
  if (ok && verbose) message(sprintf("  [delay-fig] saved %-42s %.1f x %.1f in", paste0(stem, ".pdf"), w, h))
  invisible(path)
}

# ---- Text encoding: only emit non-ASCII glyphs the PDF device can actually draw --------
# capabilities("cairo") is NOT sufficient — it reports TRUE on builds whose cairo DLL fails
# to load at first use, after which ggsave silently falls back to the base pdf() device.
# That device's single-byte encoding transliterates every "→" and "—" and emits an
# mbcsToSbcs warning per string (28 per render on this machine), so the PDF said
# "Onset -> sample" while the PNG said "Onset → sample". Probe by actually opening the
# device once, then route all figure text through .dfg_txt(), which downgrades to ASCII
# when the probe fails. Cheap (one temp file at load) and makes the two formats agree.
.dfg_probe_cairo <- function() {
  if (!isTRUE(capabilities("cairo"))) return(FALSE)
  tf <- tempfile(fileext = ".pdf"); n0 <- length(grDevices::dev.list()); bad <- FALSE
  ok <- tryCatch(
    withCallingHandlers({ grDevices::cairo_pdf(tf, width = 1, height = 1); TRUE },
                        warning = function(w) { bad <<- TRUE; invokeRestart("muffleWarning") }),
    error = function(e) FALSE)
  while (length(grDevices::dev.list()) > n0) try(grDevices::dev.off(), silent = TRUE)
  unlink(tf)
  isTRUE(ok) && !bad
}
.DFG_UNICODE <- .dfg_probe_cairo()

.dfg_txt <- function(s) {
  if (is.null(s) || !length(s) || .DFG_UNICODE) return(s)
  # No padding on the dash replacements: the source strings already carry their own spaces
  # (" — "), so adding more yields "delay  -  the".
  for (r in list(c("→", "->"), c("—", "-"), c("½", "1/2"), c("·", "|"), c("×", "x")))
    s <- gsub(r[1], r[2], s, fixed = TRUE)
  s
}
.dfg_pretty <- function(s) .dfg_txt(gsub("->", "→", s, fixed = TRUE))

# ---- Distribution helpers ----------------------------------------------------
# Moment-match a mean/sd pair onto the family's native parameters. EpiDist reports the
# posterior of (mean, sd) rather than the raw Stan parameters, so this is how a posterior
# draw becomes a drawable curve.
.dfg_par_from_mean_sd <- function(family, m, s) {
  if (!is.finite(m) || !is.finite(s) || m <= 0 || s <= 0) return(NULL)
  switch(tolower(family),
    lognormal = ,
    lnorm = c(meanlog = log(m^2 / sqrt(s^2 + m^2)), sdlog = sqrt(log(1 + (s / m)^2))),
    gamma = c(shape = (m / s)^2, rate = m / s^2),
    NULL)
}

.dfg_dens <- function(family, p, x) {
  if (is.null(p)) return(rep(NA_real_, length(x)))
  tryCatch(switch(tolower(family),
    gamma     = stats::dgamma(x, shape = p[["shape"]], rate = p[["rate"]]),
    lognormal = ,
    lnorm     = stats::dlnorm(x, meanlog = p[["meanlog"]], sdlog = p[["sdlog"]]),
    weibull   = stats::dweibull(x, shape = p[["shape"]], scale = p[["scale"]]),
    exp       = stats::dexp(x, rate = p[["rate"]]),
    rep(NA_real_, length(x))), error = function(e) rep(NA_real_, length(x)))
}

.dfg_cdf <- function(family, p, x) {
  if (is.null(p)) return(rep(NA_real_, length(x)))
  tryCatch(switch(tolower(family),
    gamma     = stats::pgamma(x, shape = p[["shape"]], rate = p[["rate"]]),
    lognormal = ,
    lnorm     = stats::plnorm(x, meanlog = p[["meanlog"]], sdlog = p[["sdlog"]]),
    weibull   = stats::pweibull(x, shape = p[["shape"]], scale = p[["scale"]]),
    exp       = stats::pexp(x, rate = p[["rate"]]),
    rep(NA_real_, length(x))), error = function(e) rep(NA_real_, length(x)))
}

.dfg_qf <- function(family, p, q) {
  if (is.null(p)) return(NA_real_)
  tryCatch(switch(tolower(family),
    gamma     = stats::qgamma(q, shape = p[["shape"]], rate = p[["rate"]]),
    lognormal = ,
    lnorm     = stats::qlnorm(q, meanlog = p[["meanlog"]], sdlog = p[["sdlog"]]),
    weibull   = stats::qweibull(q, shape = p[["shape"]], scale = p[["scale"]]),
    exp       = stats::qexp(q, rate = p[["rate"]]),
    NA_real_), error = function(e) NA_real_)
}

# Shared x-range for a delay's density and CDF panels.
# QUANTILE-based, not max-based: these delays are heavy-tailed and a handful of 50-day
# outliers stretched every panel to the MAX_DELAY cap, squeezing the bulk of the mass —
# and the naive/marginal separation the figure exists to show — into the leftmost tenth of
# the axis. The 99.5th empirical percentile and the fitted 99th percentile together keep
# the visible range honest about where the mass is while still showing the corrected
# curve's tail. Both panels call this so their axes agree.
.dfg_xmax <- function(x, curves) {
  x <- x[is.finite(x) & x >= 0]
  cap <- get0("MAX_DELAY", ifnotfound = 60L)
  if (!length(x)) return(cap)
  q_emp <- unname(stats::quantile(x, 0.995, na.rm = TRUE))
  q_fit <- suppressWarnings(max(vapply(curves, function(cv) .dfg_qf(cv$family, cv$params, 0.99),
                                       numeric(1)), na.rm = TRUE))
  if (!is.finite(q_fit)) q_fit <- 0
  min(max(q_emp, q_fit, 5, na.rm = TRUE), cap, max(x) + 2)
}

.dfg_mean_sd <- function(family, p) {
  if (is.null(p)) return(c(mean = NA_real_, sd = NA_real_))
  out <- tryCatch(switch(tolower(family),
    gamma     = c(p[["shape"]] / p[["rate"]], sqrt(p[["shape"]]) / p[["rate"]]),
    lognormal = ,
    lnorm     = { m <- exp(p[["meanlog"]] + 0.5 * p[["sdlog"]]^2)
                  c(m, m * sqrt(exp(p[["sdlog"]]^2) - 1)) },
    weibull   = { g1 <- gamma(1 + 1 / p[["shape"]]); g2 <- gamma(1 + 2 / p[["shape"]])
                  c(p[["scale"]] * g1, p[["scale"]] * sqrt(g2 - g1^2)) },
    exp       = c(1 / p[["rate"]], 1 / p[["rate"]]),
    c(NA_real_, NA_real_)), error = function(e) c(NA_real_, NA_real_))
  c(mean = unname(out[1]), sd = unname(out[2]))
}

# Largest vertical gap between the empirical CDF of the (right-truncated) observations
# and a fitted CDF. Reported as a DESCRIPTIVE distance, never as a test statistic: under
# interval censoring and right truncation the KS null distribution does not apply, and for
# the truncation-corrected marginal a LARGE gap is expected by construction.
#
# CENSORING-ALIGNED. Delays are recorded as whole days, and every fit here treats an integer
# d as the interval [max(0, d-0.5), d+0.5] (.make_cens_df in 04c). So P(recorded <= d) is
# F(d + 0.5), NOT F(d), and the comparison must be made at the interval's upper bound.
# Comparing against F(d) instead is dominated by the zero atom — 14% of onset->sample delays
# are same-day, and every continuous fit has F(0) = 0 — which returned an identical 0.142 for
# fits whose means differ by 1.4 days, i.e. a statistic that could not tell them apart.
.dfg_ecdf_gap <- function(x, family, p) {
  x <- x[is.finite(x) & x >= 0]
  if (!length(x) || is.null(p)) return(NA_real_)
  u  <- sort(unique(x))
  Fn <- stats::ecdf(x)(u)                  # P(recorded <= d)
  Ff <- .dfg_cdf(family, p, u + 0.5)       # model's P(recorded <= d)
  if (all(is.na(Ff))) return(NA_real_)
  max(pmax(abs(Fn - Ff), abs(c(0, utils::head(Fn, -1)) - Ff)), na.rm = TRUE)
}

# ---- Assemble the drawable estimator set for one delay -----------------------
# Returns a list of curve specs: label, family, params, mean, sd, mean CrI, n, plus the
# posterior draws (marginal/naive) needed for the credible band. Ordered so the pipeline's
# own estimator is first.
.dfg_curves_for <- function(delay_name, cens_fits, draws_all) {
  curves <- list()
  dr <- Filter(function(d) identical(d$delay, delay_name), draws_all)
  # EpiDist: prefer the gamma parameterisation (same rule as 04c's §1.2 pick, so the drawn
  # curve is the one that actually reached the imputation), else whatever converged.
  for (mt in c("marginal", "naive")) {
    cand <- Filter(function(d) identical(d$model_type, mt), dr)
    if (!length(cand)) next
    fams <- vapply(cand, function(d) d$family, character(1))
    pick <- cand[[if ("gamma" %in% fams) match("gamma", fams) else 1L]]
    m <- stats::median(pick$mean, na.rm = TRUE); s <- stats::median(pick$sd, na.rm = TRUE)
    p <- .dfg_par_from_mean_sd(pick$family, m, s)
    if (is.null(p)) next
    curves[[length(curves) + 1L]] <- list(
      label = if (mt == "marginal") "EpiDist marginal (truncation-corrected)"
              else "EpiDist naive (no truncation correction)",
      estimator = paste0("epidist_", mt), family = pick$family, params = p,
      mean = m, sd = s,
      mean_lo = unname(stats::quantile(pick$mean, 0.025, na.rm = TRUE)),
      mean_hi = unname(stats::quantile(pick$mean, 0.975, na.rm = TRUE)),
      n = pick$n, frac_complete = pick$frac_complete,
      draws = if (mt == "marginal") pick else NULL)   # band on the headline curve only
  }
  # Interval-censored MLE, AIC-best family.
  if (!is.null(cens_fits) && length(cens_fits)) {
    fam <- .cens_best_fam(cens_fits)
    if (!is.na(fam)) {
      p  <- cens_fits[[fam]]$estimate
      ms <- .dfg_mean_sd(fam, p)
      curves[[length(curves) + 1L]] <- list(
        label = "Interval-censored MLE", estimator = "interval_censored_mle",
        family = fam, params = p, mean = unname(ms["mean"]), sd = unname(ms["sd"]),
        mean_lo = NA_real_, mean_hi = NA_real_, n = NA_integer_,
        frac_complete = NA_real_, draws = NULL,
        aic = cens_fits[[fam]]$aic)
    }
  }
  curves
}

# 95% pointwise credible band for a fitted density, from the posterior draws. Capped at
# `max_draws` curves (deterministic thinning) — 300 is visually indistinguishable from
# 2000 and keeps the grid evaluation cheap.
.dfg_density_band <- function(draws, xg, max_draws = 300L) {
  if (is.null(draws) || !length(draws$mean)) return(NULL)
  k <- length(draws$mean)
  idx <- if (k > max_draws) unique(round(seq(1, k, length.out = max_draws))) else seq_len(k)
  mat <- vapply(idx, function(i) {
    p <- .dfg_par_from_mean_sd(draws$family, draws$mean[i], draws$sd[i])
    .dfg_dens(draws$family, p, xg)
  }, numeric(length(xg)))
  if (!is.matrix(mat) || !ncol(mat)) return(NULL)
  data.frame(x = xg,
             lo = apply(mat, 1, stats::quantile, 0.025, na.rm = TRUE),
             hi = apply(mat, 1, stats::quantile, 0.975, na.rm = TRUE))
}

.dfg_legend_lab <- function(cv) {
  ci <- if (is.finite(cv$mean_lo) && is.finite(cv$mean_hi))
    sprintf(" [%.1f-%.1f]", cv$mean_lo, cv$mean_hi) else ""
  .dfg_txt(sprintf("%s — %s, mean %.2f d%s", cv$label, toupper(cv$family), cv$mean, ci))
}

# ---- Panel A: empirical histogram + fitted densities -------------------------
.dfg_panel_density <- function(x, curves, title, subtitle_extra = NULL, xmax = NULL,
                               show_legend = TRUE) {
  x <- x[is.finite(x) & x >= 0]
  if (!length(x)) return(NULL)
  if (is.null(xmax)) xmax <- .dfg_xmax(x, curves)
  xg <- seq(0.01, xmax, length.out = 400)

  p <- ggplot2::ggplot(data.frame(delay = x), ggplot2::aes(x = .data$delay)) +
    # boundary = -0.5, NOT 0. Every fit in this pipeline treats an integer delay d as the
    # interval [max(0, d-0.5), d+0.5] (.make_cens_df, 04c), and this file relies on that
    # convention everywhere else (.dfg_ecdf_gap evaluates at u + 0.5; the CDF panel steps at
    # x + 0.5). boundary = 0 gave bins [0,1), [1,2), ... which (a) pooled the d = 0 same-day
    # atom — ~14% of records, and quoted separately in this panel's own subtitle — into the
    # 1-day bar, and (b) shifted every bar half a day right of the fitted curve, visually
    # cancelling roughly half of the naive-vs-marginal truncation gap this panel exists to
    # show. 04c_dhis2_delay_windows.R plots the same data with aligned bins; the two QA
    # figures disagreed.
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = 1,
                            boundary = -0.5, fill = .DFG_BAR, colour = .DFG_BAR_EDGE,
                            linewidth = 0.25)

  # Credible band under the curves.
  band <- NULL
  for (cv in curves) if (!is.null(cv$draws)) { band <- .dfg_density_band(cv$draws, xg); break }
  if (!is.null(band))
    p <- p + ggplot2::geom_ribbon(data = band,
                                  ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
                                  inherit.aes = FALSE,
                                  fill = .DFG_COL[["EpiDist marginal (truncation-corrected)"]],
                                  alpha = 0.16)

  ld <- dplyr::bind_rows(lapply(curves, function(cv)
    data.frame(x = xg, y = .dfg_dens(cv$family, cv$params, xg),
               key = cv$label, lab = .dfg_legend_lab(cv), stringsAsFactors = FALSE)))
  if (!is.null(ld) && nrow(ld)) {
    labs_by_key <- ld[!duplicated(ld$key), c("key", "lab")]
    cols <- stats::setNames(unname(.DFG_COL[labs_by_key$key]), labs_by_key$lab)
    ltys <- stats::setNames(unname(.DFG_LTY[labs_by_key$key]), labs_by_key$lab)
    ld$lab <- factor(ld$lab, levels = labs_by_key$lab)
    p <- p +
      ggplot2::geom_line(data = ld,
                         ggplot2::aes(x = .data$x, y = .data$y, colour = .data$lab,
                                      linetype = .data$lab),
                         inherit.aes = FALSE, linewidth = 0.85, na.rm = TRUE) +
      ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
      ggplot2::scale_linetype_manual(values = ltys, drop = FALSE) +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1),
                      linetype = ggplot2::guide_legend(ncol = 1))
  }

  # Y-CAP. A gamma/lognormal with shape < 1 — which every same-day-heavy delay here fits —
  # has an unbounded density at 0. Left to autoscale, that vertical asymptote sets the axis
  # and flattens the histogram and all three curves into the bottom fifth of the panel.
  # Cap at 1.3x the tallest visible bar and CLIP (coord_cartesian, not ylim/scale limits, so
  # no data is dropped and no "removed rows" warning is emitted); note it when a curve is cut.
  # Same [d-0.5, d+0.5) bins as the geom_histogram above, so the y-axis cap is computed from
  # the bars actually drawn. With breaks at 0,1,2,... it was measuring a different histogram.
  bar_h <- graphics::hist(x[x <= xmax],
                          breaks = seq(-0.5, ceiling(xmax) + 0.5, by = 1),
                          right = FALSE, plot = FALSE)$density
  ymax  <- if (length(bar_h) && any(is.finite(bar_h))) 1.3 * max(bar_h, na.rm = TRUE) else NA_real_
  clipped <- !is.null(ld) && nrow(ld) && is.finite(ymax) && any(ld$y > ymax, na.rm = TRUE)

  sub <- .dfg_txt(sprintf("n = %d complete pairs · %.0f%% same-day · empirical mean %.2f d, median %.0f d",
                          length(x), 100 * mean(x == 0), mean(x), stats::median(x)))
  if (!is.null(subtitle_extra)) sub <- paste0(sub, "\n", subtitle_extra)
  # The naive EpiDist model cannot take exact zeros under lognormal/gamma, so 04c fits it on
  # the strictly-positive delays only. Its curve therefore has no mass at 0 and sits visibly
  # right of the first bar — a data-exclusion artefact, not evidence that the naive model is
  # worse. Say so wherever the sample sizes actually differ, so the panel is not misread.
  nv <- Filter(function(cv) identical(cv$estimator, "epidist_naive"), curves)
  if (length(nv) && is.finite(nv[[1]]$n) && nv[[1]]$n < length(x))
    sub <- paste0(sub, sprintf("\nThe naive curve is fitted on the %d non-zero delays only (lognormal/gamma cannot take exact zeros).",
                               as.integer(nv[[1]]$n)))
  if (clipped) sub <- paste0(sub, "\nDensity axis clipped: the fitted densities are unbounded at delay 0.")
  # xlim starts at -0.5, not 0: the zero bin spans [-0.5, 0.5), so clipping the panel at 0
  # rendered the same-day bar at HALF WIDTH — in the very panel the bin alignment was fixed for.
  p + ggplot2::coord_cartesian(xlim = c(-0.5, xmax),
                               ylim = if (is.finite(ymax)) c(0, ymax) else NULL,
                               expand = FALSE) +
    # 2026-09-17 streamlining: the published dhis2_delay_epidist_fits panel carries its
    # explanation in the manuscript caption, so the verbose per-panel SUBTITLE is dropped
    # (DELAY_FIG_PANEL_SUBTITLES=1 restores it for interactive QA). The short panel TITLE is
    # kept: it names WHICH delay the panel shows and is the only thing identifying it.
    ggplot2::labs(title = .dfg_pretty(title),
                  subtitle = if (identical(Sys.getenv("DELAY_FIG_PANEL_SUBTITLES", "0"), "1"))
                               sub else NULL,
                  x = "Delay (days)", y = "Density") +
    .dfg_theme() +
    (if (show_legend) ggplot2::theme() else ggplot2::theme(legend.position = "none"))
}

# ---- Panel B: empirical CDF vs fitted CDFs -----------------------------------
.dfg_panel_cdf <- function(x, curves, title, show_legend = TRUE) {
  x <- x[is.finite(x) & x >= 0]
  if (!length(x)) return(NULL)
  xmax <- .dfg_xmax(x, curves)   # SAME range as the density panel for this delay
  xg <- seq(0, xmax, length.out = 400)

  # Step plotted at d + 0.5, the upper bound of the censoring interval the fits use, so the
  # empirical and fitted curves are like-for-like: P(recorded <= d) = F(d + 0.5). Plotting it
  # at d instead puts the whole same-day atom half a day left of where every fit places it,
  # which reads as systematic misfit when it is purely a rounding convention.
  p <- ggplot2::ggplot() +
    ggplot2::stat_ecdf(data = data.frame(delay = x + 0.5), ggplot2::aes(x = .data$delay),
                       geom = "step", colour = .DFG_INK, linewidth = 0.6, pad = FALSE)
  ld <- dplyr::bind_rows(lapply(curves, function(cv)
    data.frame(x = xg, y = .dfg_cdf(cv$family, cv$params, xg),
               key = cv$label,
               lab = sprintf("%s (gap %.2f)", cv$label, .dfg_ecdf_gap(x, cv$family, cv$params)),
               stringsAsFactors = FALSE)))
  if (!is.null(ld) && nrow(ld)) {
    lk <- ld[!duplicated(ld$key), c("key", "lab")]
    ld$lab <- factor(ld$lab, levels = lk$lab)
    p <- p + ggplot2::geom_line(data = ld,
                                ggplot2::aes(x = .data$x, y = .data$y, colour = .data$lab,
                                             linetype = .data$lab),
                                linewidth = 0.85, na.rm = TRUE) +
      ggplot2::scale_colour_manual(values = stats::setNames(unname(.DFG_COL[lk$key]), lk$lab)) +
      ggplot2::scale_linetype_manual(values = stats::setNames(unname(.DFG_LTY[lk$key]), lk$lab)) +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1),
                      linetype = ggplot2::guide_legend(ncol = 1))
  }
  p + ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(title = .dfg_pretty(title),
                  subtitle = .dfg_txt(paste0("Black step = empirical CDF of the (right-truncated) observations, plotted at the d + ½ d censoring bound.\n",
                                             "Bracketed value = max |empirical - fitted| at that bound; descriptive only, not a KS test.")),
                  x = "Delay (days)", y = "Cumulative probability") +
    .dfg_theme() +
    (if (show_legend) ggplot2::theme() else ggplot2::theme(legend.position = "none"))
}

# ---- Figure C: estimator forest ----------------------------------------------
.dfg_panel_forest <- function(summary_df, lab_ref_mean = NULL) {
  # `plotted` only: a diverged runner-up family (the lognormal marginal reached 22 d on the
  # 2026-07-25 snapshot) would set the axis and compress every estimator that matters into
  # the left third. Runner-ups stay in the CSV, where they are the audit trail rather than
  # the headline.
  d <- summary_df[is.finite(summary_df$mean_d) & summary_df$plotted, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$delay_lab <- factor(.dfg_pretty(d$label), levels = rev(unique(.dfg_pretty(d$label))))
  d$key <- ifelse(d$estimator == "epidist_marginal", "EpiDist marginal (truncation-corrected)",
           ifelse(d$estimator == "epidist_naive",    "EpiDist naive (no truncation correction)",
           ifelse(d$estimator == "empirical",        "Empirical mean", "Interval-censored MLE")))
  d$key <- factor(d$key, levels = names(.DFG_COL)[names(.DFG_COL) %in% d$key])
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$mean_d, y = .data$delay_lab,
                                       colour = .data$key)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$mean_lo, xmax = .data$mean_hi),
                            position = ggplot2::position_dodge(width = 0.62),
                            linewidth = 0.7, na.rm = TRUE) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.62),
                        size = 2.1, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = .DFG_COL[levels(d$key)], drop = FALSE)
  if (!is.null(lab_ref_mean) && is.finite(lab_ref_mean))
    p <- p + ggplot2::geom_vline(xintercept = lab_ref_mean, linetype = "dotted",
                                 colour = .DFG_COL[["Lab-linelist reference"]], linewidth = 0.6)
  p + ggplot2::labs(
        title = "Mean delay by estimator",
        subtitle = .dfg_txt(paste0("Bars = 95% credible interval (EpiDist only; the censored MLE and the empirical mean are point estimates).",
                                   if (!is.null(lab_ref_mean) && is.finite(lab_ref_mean))
                                     sprintf("\nDotted line = the %.2f d lab-linelist onset→sample reference the pipeline assumed before the DHIS2 refit.", lab_ref_mean) else "")),
        x = "Mean delay (days)", y = NULL) +
    .dfg_theme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# =============================================================================
# DRIVER
# =============================================================================

#' Render the delay-distribution fit diagnostics.
#'
#' @param ll raw DHIS2 line list (as read by run_all.R step 1a).
#' @param analysis_date the run's as-of date; passed to build_dhis2_delay_populations()
#'   so the figure describes EXACTLY the records the fits used.
#' @param os_fit optional value of estimate_dhis2_onset_sample_delay() — its `cens_fits`
#'   are reused so the onset->sample MLE is not refitted.
#' @param epidist_tbl optional EpiDist summary tibble (attr(.osepi, "epidist_table")).
#'   Only used for reporting; the CURVES come from 04c's posterior-draw registry.
#' @param out_dir directory for the figures (default <OUT_DIAGNOSTICS>/delay_fits).
#' @param epidist_all fit EpiDist for the other three delays too (slow). Defaults to the
#'   DELAY_FIG_EPIDIST_ALL env flag / global, else FALSE.
#' @return invisibly, the summary data frame (also written as CSV), or NULL on a no-op.
make_delay_fit_figures <- function(ll,
                                   analysis_date = get0("ANALYSIS_DATE", ifnotfound = NULL),
                                   os_fit = NULL, epidist_tbl = NULL,
                                   out_dir = NULL, epidist_all = NULL,
                                   verbose = TRUE) {
  if (!.DFG_HAVE_PATCHWORK) {
    message("[delay-fig] patchwork not installed — skipping delay-fit figures."); return(invisible(NULL))
  }
  if (is.null(out_dir))
    out_dir <- file.path(get0("OUT_DIAGNOSTICS", ifnotfound = file.path("outputs", "diagnostics")),
                         "delay_fits")
  if (is.null(epidist_all)) {
    .e <- tolower(trimws(Sys.getenv("DELAY_FIG_EPIDIST_ALL", "")))
    epidist_all <- if (nzchar(.e)) .e %in% c("1", "true", "t", "yes", "y")
                   else isTRUE(get0("DELAY_FIG_EPIDIST_ALL", ifnotfound = FALSE))
  }

  pops <- build_dhis2_delay_populations(ll, analysis_date = analysis_date)
  if (is.null(pops) || !length(pops$delays)) {
    message("[delay-fig] no usable delay populations — skipping."); return(invisible(NULL))
  }
  if (verbose)
    message(sprintf("[delay-fig] window %s to %s (obs date %s) | %d delays",
                    pops$analysis_start, pops$trunc_date, pops$obs_date, length(pops$delays)))

  # Optionally extend the Bayesian fit to the descriptive delays. onset->sample is NOT
  # refitted here: run_all.R already fitted it in step 1a and its draws are in the registry,
  # and refitting could land on a different posterior from the one the pipeline used.
  if (epidist_all && isTRUE(get0("RUN_EPIDIST", ifnotfound = FALSE)) &&
      exists("fit_epidist_both", mode = "function")) {
    already <- unique(vapply(epidist_draws_get(), function(d) d$delay, character(1)))
    todo    <- setdiff(names(pops$delays), already)
    if (length(todo) && verbose)
      message(sprintf("[delay-fig] DELAY_FIG_EPIDIST_ALL — fitting EpiDist for: %s",
                      paste(todo, collapse = ", ")))
    for (nm in todo) tryCatch({
      d <- pops$delays[[nm]]
      fit_epidist_both(tibble::tibble(onset = d$onset, sample = d$sample), nm, pops$obs_date)
    }, error = function(e)
      message(sprintf("[delay-fig] EpiDist for %s failed (non-fatal): %s", nm, conditionMessage(e))))
  }

  draws_all <- tryCatch(epidist_draws_get(), error = function(e) list())
  if (verbose)
    message(sprintf("[delay-fig] EpiDist posteriors available: %s",
                    if (!length(draws_all)) "none (censored-MLE curves only)"
                    else paste(vapply(draws_all, function(d)
                      sprintf("%s/%s/%s", d$delay, d$model_type, d$family), character(1)),
                      collapse = ", ")))

  # ---- Per-delay curves + summary rows ---------------------------------------
  panels_dens <- list(); panels_cdf <- list(); summary_rows <- list()
  for (nm in names(pops$delays)) {
    d   <- pops$delays[[nm]]
    lbl <- pops$specs[[nm]]$lbl
    if (length(d$delay) < 5L) next
    # Reuse run_all.R's onset->sample censored fit when supplied; refit the rest (cheap).
    cens <- if (identical(nm, "onset_sample") && !is.null(os_fit$cens_fits)) os_fit$cens_fits
            else tryCatch(.fit_all_censored(d$delay, lbl), error = function(e) NULL)
    curves <- .dfg_curves_for(nm, cens, draws_all)
    if (!length(curves) && verbose)
      message(sprintf("[delay-fig] %s: no fitted curve — histogram only", lbl))

    extra <- if (identical(nm, "onset_sample"))
      "Drives the onset imputation (01_data_prep.R), the nowcast rate and the R(t) truncation." else NULL
    pd <- .dfg_panel_density(d$delay, curves, lbl, subtitle_extra = extra)
    pc <- .dfg_panel_cdf(d$delay, curves, lbl)
    if (!is.null(pd)) panels_dens[[nm]] <- pd
    if (!is.null(pc)) panels_cdf[[nm]]  <- pc

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      delay = nm, label = lbl, window = pops$window_name,
      estimator = "empirical", family = NA_character_, n = length(d$delay),
      mean_d = mean(d$delay), sd_d = stats::sd(d$delay),
      mean_lo = NA_real_, mean_hi = NA_real_, ecdf_gap = NA_real_,
      params = NA_character_, plotted = TRUE,
      note = "right-truncated; shown for reference only",
      stringsAsFactors = FALSE)
    # EVERY EpiDist sub-fit goes in the CSV, not just the plotted (gamma-preferred) pick.
    # On the 2026-07-25 snapshot the lognormal MARGINAL diverged to mean 22.1 d / sd 88.6 d
    # while the gamma marginal sat at 7.85 d — a 3x disagreement between two fits of the same
    # data. The gamma is what the pipeline uses and what the panel draws, so a reader would
    # never see that instability unless the runner-up families are recorded somewhere.
    plotted <- vapply(curves, function(cv) paste0(cv$estimator, "|", cv$family), character(1))
    for (dw in Filter(function(z) identical(z$delay, nm), draws_all)) {
      est <- paste0("epidist_", dw$model_type)
      if (paste0(est, "|", dw$family) %in% plotted) next
      pp <- .dfg_par_from_mean_sd(dw$family, stats::median(dw$mean, na.rm = TRUE),
                                  stats::median(dw$sd, na.rm = TRUE))
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        delay = nm, label = lbl, window = pops$window_name,
        estimator = est, family = dw$family, n = as.integer(dw$n),
        mean_d = stats::median(dw$mean, na.rm = TRUE), sd_d = stats::median(dw$sd, na.rm = TRUE),
        mean_lo = unname(stats::quantile(dw$mean, 0.025, na.rm = TRUE)),
        mean_hi = unname(stats::quantile(dw$mean, 0.975, na.rm = TRUE)),
        ecdf_gap = .dfg_ecdf_gap(d$delay, dw$family, pp),
        params = if (is.null(pp)) NA_character_ else
          paste(sprintf("%s=%.4f", names(pp), unname(pp)), collapse = "; "),
        plotted = FALSE,
        note = "not plotted — runner-up family for this model type",
        stringsAsFactors = FALSE)
    }
    for (cv in curves) summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      delay = nm, label = lbl, window = pops$window_name,
      estimator = cv$estimator, family = cv$family,
      n = if (is.null(cv$n) || is.na(cv$n)) length(d$delay) else as.integer(cv$n),
      mean_d = cv$mean, sd_d = cv$sd, mean_lo = cv$mean_lo, mean_hi = cv$mean_hi,
      ecdf_gap = .dfg_ecdf_gap(d$delay, cv$family, cv$params),
      params = paste(sprintf("%s=%.4f", names(cv$params), unname(cv$params)), collapse = "; "),
      plotted = TRUE,
      note = if (identical(cv$estimator, "epidist_marginal"))
        "truncation- and censoring-corrected; expected to exceed the empirical mean" else NA_character_,
      stringsAsFactors = FALSE)
  }
  if (!length(panels_dens)) {
    message("[delay-fig] nothing to draw — skipping."); return(invisible(NULL))
  }
  summary_df <- dplyr::bind_rows(summary_rows)

  # ---- Provenance caption shared by every figure ------------------------------
  # A JSON that PARSES but has no `folder` key yields NULL, and is.na(NULL) is logical(0), so
  # the `if (is.na(src))` below raised "argument is of length zero" OUTSIDE any tryCatch —
  # contradicting this file's guarantee that every entry point is guarded and non-fatal.
  # Normalise to a length-1 character here so the caption logic cannot fail.
  src <- tryCatch({
    v <- jsonlite::fromJSON(get0("LINELIST_JSON", ifnotfound = ""))$folder
    if (length(v) != 1L) NA_character_ else as.character(v)
  }, error = function(e) NA_character_)
  cap <- .dfg_txt(sprintf(
    "DHIS2 line list %s · onset window %s to %s (truncation buffer %d d) · delays capped at %d d · as-of %s",
    if (is.na(src)) "(unknown snapshot)" else src, pops$analysis_start, pops$trunc_date,
    get0("TEST_DAYS", ifnotfound = 5L), get0("MAX_DELAY", ifnotfound = 60L),
    if (is.null(analysis_date)) "n/a" else format(analysis_date)))
  stamp <- format(pops$trunc_date, "%Y%m%d")

  # Figure-level title/subtitle suppressed for the same reason (caption stays: it carries the
  # snapshot, window and as-of provenance, which a reader needs to interpret the panel).
  .show_fig_title <- identical(Sys.getenv("DELAY_FIG_TITLES", "0"), "1")
  wrap <- function(panels, title, subtitle) {
    patchwork::wrap_plots(panels, ncol = 2) +
      patchwork::plot_annotation(
        title = if (.show_fig_title) title else NULL,
        subtitle = if (.show_fig_title) subtitle else NULL, caption = cap,
        theme = ggplot2::theme(
          plot.title    = ggplot2::element_text(size = 13, face = "bold", colour = .DFG_INK),
          plot.subtitle = ggplot2::element_text(size = 9.5, colour = .DFG_MUTED),
          plot.caption  = ggplot2::element_text(size = 7.5, colour = .DFG_MUTED, hjust = 0)))
  }

  nrow_grid <- ceiling(length(panels_dens) / 2)

  # ---- A. densities -----------------------------------------------------------
  .dfg_save(
    wrap(panels_dens, .dfg_txt("DHIS2 reporting delays — fitted distributions vs observed delays"),
         .dfg_txt(paste0("Bars = observed (right-truncated) delays. Shaded band = 95% pointwise credible interval of the EpiDist marginal density.\n",
                "The marginal curve is EXPECTED to sit right of the bars: it corrects the right-truncation the histogram still carries."))),
    out_dir, sprintf("dhis2_delay_epidist_fits_%s", stamp),
    w = 13, h = 4.4 * nrow_grid, verbose = verbose)

  # ---- B. CDFs ----------------------------------------------------------------
  .dfg_save(
    wrap(panels_cdf, .dfg_txt("DHIS2 reporting delays — empirical vs fitted CDFs"),
         .dfg_txt("Binwidth-free goodness-of-fit view. The gap statistic is descriptive, not a KS test: censoring and truncation invalidate the KS null.")),
    out_dir, sprintf("dhis2_delay_fit_cdf_%s", stamp),
    w = 13, h = 4.2 * nrow_grid, verbose = verbose)

  # ---- C. estimator forest ----------------------------------------------------
  lab_ref <- tryCatch(1 / get0("DELAY_ONSET_SAMPLE_RATE", ifnotfound = NA_real_),
                      error = function(e) NA_real_)
  fp <- .dfg_panel_forest(summary_df, lab_ref_mean = lab_ref)
  if (!is.null(fp))
    .dfg_save(fp + ggplot2::labs(caption = cap) +
                ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.5,
                                                                    colour = .DFG_MUTED, hjust = 0)),
              out_dir, sprintf("dhis2_delay_mean_estimators_%s", stamp),
              w = 9.5, h = 1.1 + 0.85 * length(unique(summary_df$delay)), verbose = verbose)

  # ---- D. onset->sample focus -------------------------------------------------
  # The one delay that feeds the models gets a full-width panel of its own: it is the
  # figure that belongs in the supplement, and at 2x2 grid size the credible band and
  # the naive-vs-marginal separation are too small to judge.
  if (!is.null(panels_dens$onset_sample) && !is.null(panels_cdf$onset_sample)) {
    focus <- (panels_dens$onset_sample | panels_cdf$onset_sample) +
      patchwork::plot_annotation(
        title = .dfg_txt("Onset → sample delay — the distribution the pipeline uses"),
        subtitle = "Left: fitted densities over the observed delays. Right: fitted vs empirical CDF.",
        caption = cap,
        theme = ggplot2::theme(
          plot.title    = ggplot2::element_text(size = 13, face = "bold", colour = .DFG_INK),
          plot.subtitle = ggplot2::element_text(size = 9.5, colour = .DFG_MUTED),
          plot.caption  = ggplot2::element_text(size = 7.5, colour = .DFG_MUTED, hjust = 0)))
    .dfg_save(focus, out_dir, sprintf("dhis2_onset_sample_fit_%s", stamp),
              w = 13, h = 5.6, verbose = verbose)
  }

  # ---- E. Onset -> sample BY CASE CLASSIFICATION -------------------------------
  # The onset imputation imputes onsets for CONFIRMED records, so it must draw from the
  # confirmed delay; the pooled fit is majority test-negative on this line list and runs
  # short. run_all.R step 1a fits each stratum independently (estimate_onset_sample_strata(),
  # 04c) and write_onset_sample_long() publishes them as `<quantity>__<stratum>` rows.
  # Nothing looked at them, so the contrast that motivated the fix was unverifiable by eye.
  #
  # THIS PANEL COMPUTES NOTHING. It resolves each stratum through
  # .load_dhis2_delay_params(stratum=) -- the SAME loader effective_onset_sample_delay() uses
  # -- so the curve drawn for `confirmed` is by construction the distribution the imputation
  # draws from. Turning the published (mean, sd) into the family's native parameters is the
  # loader's own method-of-moments step, not a fit performed here.
  tryCatch({
    .ld <- get0(".load_dhis2_delay_params", mode = "function")
    .pp <- get0("DELAY_PARAMS_PATH", ifnotfound = NULL)
    if (is.function(.ld) && !is.null(.pp) && file.exists(.pp)) {
      .want <- c(confirmed = "Confirmed", not_a_case = "Test-negative",
                 suspected = "Suspected / probable")
      .specs <- lapply(names(.want), function(st)
        tryCatch(suppressWarnings(.ld(.pp, stratum = st)), error = function(e) NULL))
      names(.specs) <- names(.want)
      # A stratum resolves to the POOLED fit when its rows are absent; the loader warns and
      # sets $stratum to NA. Drawing that would put the pooled curve on the panel three times
      # under three stratum names, so drop anything that did not resolve to its own fit.
      .specs <- .specs[vapply(.specs, function(x)
        !is.null(x) && isTRUE(!is.na(x$stratum)) && is.finite(x$mean) && x$mean > 0, logical(1))]
      .pool <- tryCatch(suppressWarnings(.ld(.pp)), error = function(e) NULL)
      if (length(.specs)) {
        .lab <- function(st) unname(.want[st])
        # CDF, NOT density. Every fitted gamma here has shape < 1 (confirmed: 0.95), so its
        # density is unbounded at zero: three spikes at the origin, and the separation that
        # matters -- where the mass actually sits -- squeezed flat beneath them. The CDF is
        # bounded, and "the share of cases sampled within d days" is the reading the nowcast
        # and the imputation both act on.
        .cdf_of <- function(spec, xs) {
          pr <- spec$params
          if (identical(spec$family, "gamma") && all(c("shape", "rate") %in% names(pr)))
            stats::pgamma(xs, shape = pr[["shape"]], rate = pr[["rate"]])
          else if (identical(spec$family, "lnorm") && all(c("meanlog", "sdlog") %in% names(pr)))
            stats::plnorm(xs, meanlog = pr[["meanlog"]], sdlog = pr[["sdlog"]])
          else stats::pexp(xs, rate = 1 / spec$mean)   # Exponential summary is the last resort
        }
        .xmax <- max(25, ceiling(1.9 * max(vapply(.specs, function(x) x$mean, 0))))
        .xs   <- seq(0, .xmax, length.out = 512)
        .dd <- dplyr::bind_rows(lapply(names(.specs), function(st) tibble::tibble(
          stratum = .lab(st), x = .xs, d = .cdf_of(.specs[[st]], .xs))))
        .ord <- vapply(names(.specs), function(st) .specs[[st]]$mean, 0)
        .lev <- .lab(names(sort(.ord, decreasing = TRUE)))
        .dd$stratum <- factor(.dd$stratum, levels = .lev)

        pE1 <- ggplot2::ggplot(.dd, ggplot2::aes(x = .data$x, y = .data$d,
                                                 colour = .data$stratum)) +
          ggplot2::geom_line(linewidth = 0.75) +
          { if (!is.null(.pool) && is.finite(.pool$mean))
              ggplot2::geom_vline(xintercept = .pool$mean, linetype = "22",
                                  colour = .DFG_MUTED, linewidth = 0.45) } +
          ggplot2::scale_colour_brewer(palette = "Dark2", name = NULL) +
          ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                      limits = c(0, 1)) +
          ggplot2::labs(x = "Days from onset to sample", y = "Sampled within d days",
                        subtitle = if (!is.null(.pool) && is.finite(.pool$mean))
                          sprintf("Dashed line: the POOLED mean (%.2f d) the strata replace", .pool$mean)
                        else NULL) +
          .dfg_theme() + ggplot2::theme(legend.position = "top")

        .fr <- dplyr::bind_rows(lapply(names(.specs), function(st) {
          sp <- .specs[[st]]
          tibble::tibble(stratum = .lab(st), mean = sp$mean,
                         lo = sp$mean_lo, hi = sp$mean_hi, n = sp$n_fit)
        }))
        if (!is.null(.pool) && is.finite(.pool$mean))
          .fr <- dplyr::bind_rows(.fr, tibble::tibble(
            stratum = "Pooled (all classes)", mean = .pool$mean,
            lo = .pool$mean_lo, hi = .pool$mean_hi, n = .pool$n_fit))
        .fr$stratum <- factor(.fr$stratum, levels = rev(c(.lev, "Pooled (all classes)")))
        # n is printed because these means are not equally informed: the suspected stratum is
        # two orders of magnitude smaller than the other two and its interval shows it.
        .fr$tag <- sprintf("%.2f d  (n = %s)", .fr$mean,
                           formatC(.fr$n, format = "d", big.mark = ","))

        pE2 <- ggplot2::ggplot(.fr, ggplot2::aes(x = .data$mean, y = .data$stratum)) +
          ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                                  colour = .DFG_MUTED, linewidth = 0.7, na.rm = TRUE) +
          ggplot2::geom_point(size = 2.4, colour = .DFG_INK) +
          ggplot2::geom_text(ggplot2::aes(label = .data$tag), hjust = 0, nudge_y = 0.24,
                             size = 3, colour = .DFG_MUTED) +
          ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.30))) +
          ggplot2::labs(x = "Mean delay (days), 95% CrI", y = NULL) +
          .dfg_theme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

        .cf <- .specs[["confirmed"]]
        .sub <- if (!is.null(.cf) && !is.null(.pool) && is.finite(.pool$mean))
          sprintf(paste0("Each stratum is fitted independently, with its own interval ",
                         "censoring and truncation correction. Confirmed cases are sampled ",
                         "%.2f d later than the pooled fit implies (%.2f vs %.2f d), which is ",
                         "why the onset imputation draws from the confirmed fit."),
                  .cf$mean - .pool$mean, .cf$mean, .pool$mean)
        else "Each stratum is fitted independently, with its own censoring and truncation correction."
        .sub <- paste(strwrap(.sub, width = 116), collapse = "\n")
        figE <- (pE1 | pE2) + patchwork::plot_layout(widths = c(1.25, 1)) +
          patchwork::plot_annotation(
            title = .dfg_txt("Onset → sample delay by case classification"),
            subtitle = .sub, caption = cap,
            theme = ggplot2::theme(
              plot.title    = ggplot2::element_text(size = 13, face = "bold", colour = .DFG_INK),
              plot.subtitle = ggplot2::element_text(size = 9, colour = .DFG_MUTED),
              plot.caption  = ggplot2::element_text(size = 7.5, colour = .DFG_MUTED, hjust = 0)))
        .dfg_save(figE, out_dir, sprintf("dhis2_onset_sample_by_class_%s", stamp),
                  w = 12.4, h = 5.2, verbose = verbose)
      } else if (verbose) {
        message("  [delay-fig] no per-classification delay rows in ", basename(.pp),
                " — panel E skipped (run step 1a with RUN_EPIDIST=TRUE to write them).")
      }
    }
  }, error = function(e) message("  [delay-fig] by-classification panel skipped: ",
                                 conditionMessage(e)))

  # ---- Machine-readable companion ---------------------------------------------
  csv <- file.path(out_dir, "dhis2_delay_fit_summary.csv")
  tryCatch({
    readr::write_csv(summary_df, csv)
    if (verbose) message(sprintf("  [delay-fig] saved %-42s (%d rows)",
                                 basename(csv), nrow(summary_df)))
  }, error = function(e) message("  [delay-fig] summary CSV failed: ", conditionMessage(e)))

  # ---- CONSISTENCY CHECK against the delay the pipeline actually uses -----------
  # This module used to WRITE dhis2_delay_selected.csv (into outputs/ AND into
  # data/cfr_reference/), applying its OWN selection rule to its OWN refits. That made a
  # FIGURE module the second producer of a table METHODS.md calls "the headline estimate per
  # delay" — and the two producers disagreed: 04c/run_all step 1a wrote
  # dhis2_onset_sample_delay_params.csv with selected_estimator = epidist_marginal,
  # selected_mean = 7.665 d, while this block, running later in the same run, wrote
  # interval_censored_mle, 6.832 d. A 12% disagreement on the delay that drives the onset
  # imputation, the nowcast rate and the R(t) truncation, with the wrong one bearing the
  # name METHODS.md cites.
  #
  # A figure module does not get to publish the headline estimate. The single producer is
  # write_onset_sample_long() (04c), called from run_all.R step 1a, and the single reader is
  # effective_onset_sample_delay() (00_config.R). All this block does now is CHECK that the
  # fits it just plotted agree with that resolver, and say so loudly if they do not — which
  # is exactly what a diagnostic module is for.
  tryCatch({
    .res <- get0("effective_onset_sample_delay", mode = "function")
    if (is.function(.res)) {
      .live <- .res()
      .os <- summary_df[summary_df$delay == "onset_sample" &
                          summary_df$estimator == "epidist_marginal" &
                          summary_df$plotted %in% TRUE, , drop = FALSE]
      if (!nrow(.os))
        .os <- summary_df[summary_df$delay == "onset_sample" &
                            summary_df$estimator == "epidist_marginal", , drop = FALSE]
      if (nrow(.os) && is.finite(.os$mean_d[1]) && is.finite(.live$mean)) {
        .rel <- abs(.os$mean_d[1] - .live$mean) / .live$mean
        if (.rel > 0.02)
          warning(sprintf(paste0("[delay-fig] the onset->sample EpiDist marginal plotted here (%.3f d) ",
                                 "differs from the delay the pipeline uses (%.3f d, estimator '%s') by ",
                                 "%.1f%%. These figures then describe a delay no downstream step applies. ",
                                 "Check that run_all.R step 1a refreshed ",
                                 "data/cfr_reference/dhis2_onset_sample_delay_params.csv for THIS line list."),
                          .os$mean_d[1], .live$mean,
                          if (is.null(.live$estimator)) NA_character_ else .live$estimator, 100 * .rel),
                  call. = FALSE)
        else if (verbose)
          message(sprintf("  [delay-fig] onset->sample agrees with the pipeline resolver: %.3f d vs %.3f d",
                          .os$mean_d[1], .live$mean))
      }
    }
  }, error = function(e) message("  [delay-fig] resolver consistency check skipped: ", conditionMessage(e)))

  # One-line console QA: the number that actually reaches the imputation, next to the
  # empirical mean it corrects.
  if (verbose) {
    # `plotted` filter matters: the summary also carries runner-up EpiDist families, and
    # without it this line would quote whichever family sorted first rather than the one
    # routed into the imputation.
    os <- summary_df[summary_df$delay == "onset_sample" & summary_df$plotted, , drop = FALSE]
    em <- os$mean_d[os$estimator == "empirical"]
    mg <- os$mean_d[os$estimator == "epidist_marginal"]
    ml <- os$mean_d[os$estimator == "interval_censored_mle"]
    message(sprintf("[delay-fig] onset->sample mean (d): empirical %.2f | censored MLE %s | EpiDist marginal %s",
                    if (length(em)) em[1] else NA_real_,
                    if (length(ml)) sprintf("%.2f", ml[1]) else "n/a",
                    if (length(mg)) sprintf("%.2f", mg[1]) else "n/a"))
    # Cross-family instability alarm. Two EpiDist families fitted to the SAME data should
    # not disagree by much; when they do, one of them has not identified its tail and the
    # gamma-preferred pick is carrying a silent modelling choice into the imputation. Loud
    # in the log rather than buried in the CSV, because nobody opens the CSV on a clean run.
    for (dl in unique(summary_df$delay)) {
      mm <- summary_df[summary_df$delay == dl & summary_df$estimator == "epidist_marginal" &
                       is.finite(summary_df$mean_d), , drop = FALSE]
      if (nrow(mm) < 2) next
      if (max(mm$mean_d) > 1.5 * min(mm$mean_d)) {
        used <- if ("gamma" %in% mm$family) "gamma" else mm$family[1]  # 04c §1.2 pick rule
        message(sprintf("[delay-fig] WARNING %s: EpiDist marginal families disagree — %s. Plotted and routed into the pipeline: %s.",
                        dl, paste(sprintf("%s %.2f d", mm$family, mm$mean_d), collapse = ", "), used))
      }
    }
  }
  invisible(summary_df)
}

# Standalone: Rscript spatiotemporal/04d_delay_fit_figures.R
# (fits everything itself; env DELAY_FIG_EPIDIST_ALL=TRUE for all four Bayesian fits)
if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({ library(readr); library(jsonlite) })
  .meta <- jsonlite::fromJSON(LINELIST_JSON)
  .csv  <- file.path(LINELIST_DIR, .meta$folder, "dhis2_processed_linelist.csv")
  .ll   <- readr::read_csv(.csv, col_types = readr::cols(.default = "c"),
                           show_col_types = FALSE, na = c("", "NA", "N/A"))
  .osf  <- estimate_dhis2_onset_sample_delay(.ll, analysis_date = ANALYSIS_DATE)
  .ose  <- if (exists("estimate_dhis2_onset_sample_epidist", mode = "function"))
    estimate_dhis2_onset_sample_epidist(.ll, analysis_date = ANALYSIS_DATE) else NULL
  make_delay_fit_figures(.ll, analysis_date = ANALYSIS_DATE, os_fit = .osf,
                         epidist_tbl = attr(.ose, "epidist_table"))
}
