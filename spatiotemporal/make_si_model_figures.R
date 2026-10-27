# =============================================================================
# make_si_model_figures.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
# Supplementary MODEL figures, in the main-text house style (figure_style.R):
# no embedded titles or subtitles, panel tags only, captions live in the manuscript.
#
#   FigureS5  Cross-validated discrimination   A) ROC curves, B) precision-recall curves,
#                                              featured model vs Figure 2A's comparators
#                                              (FIGURE_BASELINE_METHODS), both horizons, on
#                                              the RAW probabilities and the shared scoring
#                                              support, pooled across the leave-future-out
#                                              folds with a band for round-to-round spread.
#                                              (The all-rounds companion of Figure 2B is
#                                              FigureS14, make_manuscript_figures.R.)
#   FigureS6  MCMC convergence                 A) max R-hat per model, B) minimum bulk and
#                                              tail effective sample size, C) divergent
#                                              transitions — every fitted model, one row each.
#   FigureS7  Calibration over time            A) expected vs observed invasions per round,
#                                              B) the O/E ratio over the rounds, C) the
#                                              pooled reliability curve, for the featured
#                                              model, one figure per horizon (_h1, _h2).
#   FigureS8  Generation-time sensitivity      A) cross-validated skill of the featured model
#                                              under the short / medium / long generation
#                                              time, B) today's invasion probabilities under
#                                              each, for the current top-ranked zones.
#   FigureS9  Time-varying beta_t              the time-varying forms of the stage-1
#                                              specification against fixed beta (run_all.R's
#                                              selection stage, tv_comparison_*):
#                                              A) cumulative expected (recalibrated) vs
#                                              observed invasions per form and horizon,
#                                              B) log-score gain over fixed beta, its 90%
#                                              zone-bootstrap interval and the smaller
#                                              Bonferroni lower bound (reported; the featured
#                                              model is chosen by the selection composite).
#   FigureS13 Prior against posterior          one panel per parameter of the featured model:
#                                              the weakly-informative prior it was fitted
#                                              under, drawn against its posterior, so a
#                                              reader can see which parameters the data
#                                              identify and which the prior carries.
#
# EVERYTHING IS READ, NOTHING IS REFITTED. The curves come from the saved cross-validation
# (lfo_results.rds), the diagnostics from bayes_convergence_diagnostics.csv, the time-varying
# comparison from the selection stage's tv_comparison_by_date.csv and tv_comparison_guard.csv,
# the posterior draws and the priors they were fitted under from bayes_posterior_draws.csv,
# and the probabilities from the published risk tables. A figure that recomputed its own
# model would be free to disagree with the pipeline.
#
# Outputs -> outputs/key_outputs/manuscript_figures/{FigureS5..S9,FigureS13}.{pdf,png} + panels/
# Run:  Rscript make_si_model_figures.R   (from anywhere; paths are anchored with here())
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(scales); library(forcats); library(purrr)
})
options(dplyr.summarise.inform = FALSE)

HERE <- file.path(here::here(), "spatiotemporal")
OUT  <- file.path(HERE, "outputs")
source(file.path(HERE, "00_config.R"))
stop_if_data_definition_sensitivity("make_si_model_figures.R")   # hard-coded published paths
source(file.path(HERE, "forecast_scale.R"))
# invasion_common_cells(): the (fold x zone) support every method is scored on in
# invasion_evaluation.csv. FigureS5 draws its curves on that support so they describe the same
# rows as the AUCs it prints. Sourced before 23/20, the order load_pipeline_modules() uses; it
# defines functions only.
source(file.path(HERE, "16_invasion_eval.R"))
# reliability_breaks()/reliability_curve(): the calibration panel and the published
# diagnostics table must come from ONE estimator and ONE binning, so the figure calls the
# same pure functions the pipeline writes the table with rather than re-deriving bins here.
source(file.path(HERE, "23_prospective_eval.R"))
# .PARAM_LABEL / .lk(): the parameter vocabulary the pipeline's own parameter figures use.
# FigureS13 names the same coefficients, so it reads them from there rather than carrying a
# second list that would drift the moment a covariate is renamed. This module defines
# functions only; it runs nothing on source().
source(file.path(HERE, "20_forecast_detail.R"))
FIG_DIR   <- fs_out_dir(file.path(OUT, "key_outputs", "manuscript_figures"))
PANEL_DIR <- fs_out_dir(file.path(OUT, "key_outputs", "manuscript_figures", "panels"))
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
source(file.path(HERE, "figure_style.R"))


# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("make_si_model_figures.R")) {
  message("[si-model] sourced, not run: nothing is written. Use `Rscript make_si_model_figures.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

# -----------------------------------------------------------------------------
# LOAD (saved artefacts only)
# -----------------------------------------------------------------------------
.read_or_null <- function(f, what) {
  if (!file.exists(f)) { message(sprintf("[si] %s not found (%s) — its figure is skipped.",
                                         what, basename(f))); return(NULL) }
  suppressWarnings(read_csv(f, show_col_types = FALSE))
}

lfo <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))
ev  <- read_csv(file.path(OUT, "diagnostics", "invasion_evaluation.csv"), show_col_types = FALSE)
# Put the selected probability scale into `p_invasion` exactly as the main-text suite does, so
# a calibration panel here and a Brier score there describe the same numbers.
PSCALE_COL <- fs_apply_lfo_scale(lfo)

# FEATURED MODEL: READ, DO NOT RECOMPUTE — the pipeline's own record of what it picked.
FEATURED <- read_featured_model(file.path(OUT, "key_outputs", "model_selection.json"))
message(sprintf("[si] featured model: %s (%s)", FEATURED, model_pretty_label(FEATURED)))

# FigureS5's comparators ARE Figure 2A's: the published FIGURE_BASELINE_METHODS, in the
# colours and under the names fs_panelA_scales() gives them (figure_style.R). This used to be a
# private set — Gravity-B4 / Distance-B1 / Adjacency-B7, the comparators Figure 2 retired — in
# a private palette in which #D55E00 was Gravity-B4 here and the Flowminder cohort inflow in
# Figure 2A, and #009E73 Distance-B1 here and Adjacency-B7 there. The membership check is
# inside build_figS5(), so a missing comparator fails that figure rather than the whole suite.
BASELINES <- fs_baseline_methods()
HZ_LAB    <- c(`1` = "1 week ahead", `2` = "2 weeks ahead")

# Date-tick thinning, shared by the panels that run over forecast rounds. Same rule as the
# publication suite: at most `max_n` labels, endpoints always kept, so a growing fold sequence
# never overprints its own axis.
.thin_si_breaks <- function(x, max_n = 6L) {
  cs <- sort(unique(as.Date(x)))
  if (length(cs) <= max_n) return(cs)
  cs[sort(unique(c(1L, round(seq(1, length(cs), length.out = max_n)), length(cs))))]
}

# `ok` records which figures actually built, so the run log says what is missing rather than
# the absence being discovered when the manuscript is assembled.
#
# A FAILED FIGURE MUST NOT LEAVE THE PREVIOUS RUN'S FILE BEHIND, AND MUST FAIL THE SCRIPT. This
# used to catch the error, warn, and carry on to an exit status of 0 — and run_all.R judges a
# suite by its exit status — while the earlier run's PDF and PNG stayed on disk under the same
# name. FigureS8 shipped that way: its builder had been failing for weeks and the September
# figure, built on a different featured model, was still in manuscript_figures/. So every file
# a figure writes is deleted BEFORE it is rebuilt (by the suite's naming convention: the
# composite "<nm>_*" in FIG_DIR and its lettered panels "FS<n><letter>_*" in PANEL_DIR), and the
# script stops with a non-zero status at the end if any figure failed.
# The stems are listed at the call site (bottom of the file) rather than matched by a pattern
# such as "FigureS5_*", which would also sweep up another suite's figure that happens to share
# the number; a "_h<N>" horizon suffix is allowed, as FIGURE_KEEP allows it.
.clear_figure_files <- function(nm, stems) {
  gone <- unlist(lapply(c(FIG_DIR, PANEL_DIR), function(d)
    unlist(lapply(stems, function(st)
      list.files(d, pattern = sprintf("^%s(_h[0-9]+)?\\.(pdf|png)$", st), full.names = TRUE)))))
  if (length(gone)) {
    unlink(gone)
    message(sprintf("[si] %s: removed %d previous file(s) before rebuilding", nm, length(gone)))
  }
  invisible(gone)
}
ok <- list()
.build <- function(nm, f, stems) {
  .clear_figure_files(nm, stems)
  ok[[nm]] <<- tryCatch({ f(); TRUE },
    error = function(e) { warning(sprintf("[si] %s FAILED: %s", nm, conditionMessage(e)),
                                  call. = FALSE); FALSE })
}

# -----------------------------------------------------------------------------
# FIGURE S5 — cross-validated ROC and precision-recall curves
# -----------------------------------------------------------------------------
# WHY BOTH CURVES, AND WHY PR IS THE ONE TO READ. Invasion is a rare event: the pooled base
# rate is ~0.8% at h=1. ROC is computed against the huge negative class, so it is close to 1
# for any model that puts the invaded zones near the top of ~500 — informative about ranking,
# but flattering in absolute terms. The precision-recall curve is scaled by the base rate and
# is the operationally meaningful one: at a given recall it says what fraction of the zones
# you would have monitored were actually invaded. Both are drawn, side by side, so neither
# can be quoted without the other.
#
# POOLED, WITH A FOLD-VARIATION BAND. The heavy line is the curve over all scored rows — the
# same pool the published AUCs in invasion_evaluation.csv are computed on. The shaded band is
# the 5th-95th percentile ACROSS ROUNDS of the featured model's own curve, so the reader can
# see how much of the pooled curve is round-to-round variation rather than a stable property.
# It is a band and not a spaghetti of per-round curves because with 1-8 invasions in a round
# the individual curves are step functions of a handful of points: drawn raw they are visual
# noise that hides the pooled curve, and the question they are there to answer is about the
# spread, not about any one round. Rounds with no invasion contribute no curve at all
# (precision and recall are undefined without a positive) and are absent by construction.

#' ROC points for one (score, outcome) vector, tie-aware.
#' Ties are collapsed to a single operating point — tied zones cannot be separated by any
#' threshold, so drawing them as separate points would claim a discrimination that does not
#' exist. Returns FPR/TPR including the (0,0) and (1,1) endpoints.
.roc_points <- function(p, y) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NULL)
  o <- order(p, decreasing = TRUE); ps <- p[o]; ys <- y[o]
  keep <- c(ps[-length(ps)] != ps[-1], TRUE)
  tibble(x = c(0, (cumsum(ys == 0) / n0)[keep]),
         y = c(0, (cumsum(ys == 1) / n1)[keep]))
}

#' Precision-recall points, tie-aware, matching .auc_pr() in 16_invasion_eval.R.
#' x = recall, y = precision.
.pr_points <- function(p, y) {
  y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  P <- sum(y == 1); if (P == 0) return(NULL)
  o <- order(p, decreasing = TRUE); ps <- p[o]; ys <- y[o]
  tp <- cumsum(ys == 1); fp <- cumsum(ys == 0)
  keep <- c(ps[-length(ps)] != ps[-1], TRUE)
  tibble(x = (tp / P)[keep], y = (tp / (tp + fp))[keep])
}

#' The vertices FigureS5 hands to its geoms, in drawing order — ONE convention per curve, the
#' one whose area is the AUC the figure prints.
#'
#' ROC: the tie-collapsed operating points in threshold order (FPR and TPR both
#'   non-decreasing; sorted on (x, y) so a vertical run climbs). Joined by straight lines
#'   (geom_path), the area under them is the TRAPEZOID sum, which is exactly AUC-ROC with tied
#'   scores counted one half (.auc_roc(), 16_invasion_eval.R): a tie group holding positives
#'   and negatives is a diagonal segment, not a staircase.
#' PR: the tie-collapsed operating points in threshold order (recall non-decreasing; at a
#'   fixed recall precision falls as negatives are added, so sorted on (x, -y)), with (0, the
#'   first precision) PREPENDED, drawn by geom_step(direction = "vh"). Each recall interval
#'   (R_{k-1}, R_k] is then drawn at the precision of the operating point that first reaches
#'   R_k, and the area under the steps is sum((R_k - R_{k-1}) * P_k) -- the tie-aware average
#'   precision .auc_pr() publishes. Straight lines would be wrong here: the achievable region
#'   between two PR operating points is not the line between them (Davis & Goadrich 2006, ICML).
#'
#' These used to be drawn by geom_step()'s default "hv", which holds the PREVIOUS point's height
#' across each interval. For PR that drew every interval at the lowest precision of the recall
#' level before it, and for ROC it drew a staircase under every diagonal tie segment, so on both
#' panels the area under the drawn curve was not the AUC printed beside it.
.s5_drawn <- function(cur, cv) {
  if (is.null(cur) || !nrow(cur)) return(cur)
  if (identical(cv, "ROC")) return(cur[order(cur$x, cur$y), , drop = FALSE])
  cur <- cur[order(cur$x, -cur$y), , drop = FALSE]
  dplyr::bind_rows(tibble(x = 0, y = cur$y[1]), cur)
}

#' Area under the curve AS DRAWN from .s5_drawn()'s vertices, by each geom's own rule:
#' the polyline (trapezoids) for ROC; geom_step(direction = "vh"), which draws the interval
#' (x_i, x_{i+1}) at height y_{i+1}, for PR.
.s5_drawn_area <- function(v, cv) {
  if (identical(cv, "ROC")) sum(diff(v$x) * (head(v$y, -1) + tail(v$y, -1)) / 2)
  else sum(diff(v$x) * v$y[-1])
}

#' Evaluate one round's DRAWN curve (.s5_drawn() vertices) on a common grid, so curves from
#' different rounds can be compared pointwise -- by the SAME convention the curve is drawn
#' with, so the band describes the curves the reader sees.
#'
#' ROC: linear interpolation along the polyline. g falls between the last vertex with x <= g
#'   (the TOP of any vertical run there) and the next one, whose x is strictly greater; at
#'   g >= the last x the curve's final value is held.
#' PR: the height of the step over g: the precision of the first vertex with x >= g (in the
#'   drawn (x, -y) order, the operating point that first reaches that recall). Beyond the last
#'   recall reached no value is defined and NA is returned, which keeps those grid points out
#'   of the quantiles rather than inventing a value at a recall no round reached.
.on_grid <- function(v, grid, cv) {
  out <- rep(NA_real_, length(grid))
  if (is.null(v) || !nrow(v)) return(out)
  n <- nrow(v)
  if (identical(cv, "ROC")) {
    j  <- findInterval(grid, v$x)                 # last vertex with x <= g
    ok <- j >= 1L
    end <- ok & j == n
    out[end] <- v$y[n]
    mid <- which(ok & j < n)
    if (length(mid)) {
      a <- j[mid]; b <- a + 1L                    # v$x[b] > grid >= v$x[a], so no zero width
      out[mid] <- v$y[a] + (grid[mid] - v$x[a]) / (v$x[b] - v$x[a]) * (v$y[b] - v$y[a])
    }
  } else {
    k  <- findInterval(grid, v$x, left.open = TRUE) + 1L   # first vertex with x >= g
    ok <- k <= n
    out[ok] <- v$y[k[ok]]
  }
  out
}

#' The rows one (method, horizon) is actually SCORED on: at-risk only, finite probability,
#' finite outcome. Every panel in this file filters identically, so a curve, a calibration bin
#' and a quoted AUC all describe the same set of rows.
.scored <- function(m, h) {
  lfo %>% filter(method == m, horizon == h, is.finite(p_invasion),
                 is.finite(is_new_invasion), !(was_active_before %in% TRUE))
}

#' The rows FigureS5 draws for one (method, horizon): .scored() restricted to the shared
#' (fold x zone) support invasion_evaluation.csv scores every method on, so a curve and the AUC
#' printed beside it are computed from the same rows. Without the restriction a method that
#' covers cells some other method lacks is drawn on more rows than it was scored on.
.s5_cells <- local({
  memo <- list()
  function(h) {
    k <- as.character(h)
    if (is.null(memo[[k]])) {
      # The support the evaluation resolved: on the RAW probabilities it was computed from.
      memo[[k]] <<- invasion_common_cells(dplyr::mutate(lfo, p_invasion = p_invasion_raw), h)
      if (is.null(memo[[k]])) stop("no shared scoring support at h=", h)
    }
    memo[[k]]
  }
})
.s5_scored <- function(m, h) {
  d <- .scored(m, h)
  d[paste(d$fold_id, d$health_zone, sep = "\r") %in% .s5_cells(h), , drop = FALSE]
}

build_figS5 <- function() {
  # THE COMPARATORS OF FIGURE 2A, in its colours and under its names (see BASELINES above).
  miss <- setdiff(c(FEATURED, BASELINES), unique(lfo$method))
  if (length(miss))
    stop("the cross-validation holds no rows for ", paste(miss, collapse = ", "),
         " — FigureS5 draws the featured model against FIGURE_BASELINE_METHODS")
  sc     <- fs_panelA_scales(FEATURED, FEATURED_LBL, BASELINES, who = "figS5")
  models <- c(FEATURED, BASELINES)
  hz     <- sort(unique(lfo$horizon))
  GRID   <- seq(0, 1, by = 0.005)
  pooled <- list(); band <- list()
  # DISCRIMINATION IS DRAWN ON THE RAW PROBABILITIES — the suite's convention, and the scale of
  # the AUCs this figure prints (auc_roc / auc_pr, the raw columns). Within a round the
  # recalibration is strictly monotone, so the per-round band is identical on either scale; but
  # the POOLED curve mixes rounds that carry different prequential factors, so on the
  # recalibrated scale the pooled ranking — and with it the pooled curve — changes. Drawing the
  # recalibrated curve beside the raw AUC put two different pooled rankings on one panel.
  for (cv in c("ROC", "PR")) {
    pts <- if (identical(cv, "ROC")) .roc_points else .pr_points
    for (h in hz) {
      for (m in models) {
        d <- .s5_scored(m, h); if (!nrow(d)) next
        cur <- pts(d$p_invasion_raw, d$is_new_invasion); if (is.null(cur)) next
        # The DRAWN vertices (.s5_drawn()): the area check below and the geoms read the same rows.
        pooled[[length(pooled) + 1L]] <- mutate(.s5_drawn(cur, cv), method = m, horizon = h,
                                                curve = cv)
      }
      # Fold-variation band, featured model only.
      dF <- .s5_scored(FEATURED, h); if (!nrow(dF)) next
      per <- lapply(sort(unique(dF$fold_id)), function(fd) {
        df <- dF %>% filter(fold_id == fd)
        .on_grid(.s5_drawn(pts(df$p_invasion_raw, df$is_new_invasion), cv), GRID, cv)
      })
      per <- per[!vapply(per, function(v) all(is.na(v)), logical(1))]
      if (length(per) < 3L) next          # a band over fewer than three rounds is not one
      M <- do.call(cbind, per)
      band[[length(band) + 1L]] <- tibble(
        x  = GRID, curve = cv, horizon = h,
        lo = apply(M, 1, stats::quantile, 0.05, na.rm = TRUE, names = FALSE),
        hi = apply(M, 1, stats::quantile, 0.95, na.rm = TRUE, names = FALSE),
        n_folds = apply(M, 1, function(v) sum(is.finite(v)))) %>%
        # Drop grid points reached by fewer than three rounds: a "5th-95th percentile" over
        # one or two rounds is a range, not a percentile, and it is exactly at the extreme
        # recall end that only one lucky round reaches.
        filter(n_folds >= 3L, is.finite(lo), is.finite(hi))
    }
  }
  if (!length(pooled)) stop("no scorable rows for the ROC/PR curves")
  pooled <- bind_rows(pooled); band <- if (length(band)) bind_rows(band) else NULL

  # THE PRINTED AUC MUST BE THE AREA UNDER THE DRAWN CURVE. Recomputed from the vertices the
  # geoms are handed, IN THE ORDER THEY ARE HANDED THEM, by each geom's own rule
  # (.s5_drawn_area(): trapezoids under the ROC polyline; geom_step "vh" heights under the PR
  # steps), and compared with the published value. A mismatch means the curve and the number
  # describe different rows, a different scale or a different drawing convention -- each a
  # defect this figure has carried -- so it is fatal.
  chk <- pooled %>% group_by(method, horizon, curve) %>%
    group_modify(~ tibble(area = .s5_drawn_area(.x, .y$curve))) %>% ungroup() %>%
    mutate(col = ifelse(curve == "ROC", "auc_roc", "auc_pr")) %>%
    rowwise() %>%
    mutate(pub = { v <- ev[[col]][ev$method == method & ev$horizon == horizon]
                   if (length(v) == 1L) v else NA_real_ }) %>% ungroup()
  bad <- chk %>% filter(is.finite(pub), abs(area - pub) > 1e-6)
  if (nrow(bad))
    stop("FigureS5's drawn curves do not reproduce the published AUCs: ",
         paste(sprintf("%s h=%d %s drawn %.4f vs published %.4f", bad$method, bad$horizon,
                       bad$curve, bad$area, bad$pub), collapse = "; "))

  # Reference lines: the ROC chance diagonal, and for PR the pooled base rate — the precision
  # a random watch-list achieves. Computed PER HORIZON from the rows actually scored, not
  # assumed, because the h=2 base rate is not the h=1 one (overlapping outcome windows).
  base_rate <- bind_rows(lapply(hz, function(h)
    tibble(horizon = h, br = mean(.s5_scored(FEATURED, h)$is_new_invasion))))

  lab_of <- function(v) unname(HZ_LAB[as.character(v)])
  # The AUCs are QUOTED FROM invasion_evaluation.csv (checked against the drawn curves above),
  # under the same names as the legend. A model with no row is omitted rather than printed as NA.
  .auc_lab <- function(h, col) {
    e <- ev %>% filter(horizon == h) %>% slice(match(models, method))
    keep <- !is.na(e$method) & is.finite(e[[col]])
    if (!any(keep)) return("")
    # THE HEADER SAYS WHAT IS POOLED. These curves rank every round's rows in ONE list, so the
    # average precision rewards a probability LEVEL on the rounds that had invasions; a case-blind
    # fixed list can lead here while trailing on the per-round skill the selection and Figure 2
    # use. Unlabelled, the two figures read as contradicting each other.
    hdr <- if (identical(col, "auc_roc")) "AUC-ROC, all rounds pooled"
           else "Average precision, all rounds pooled\n(per-round skill: Figure 2C)"
    paste(c(hdr, sprintf("%s  %.3f", unname(sc$labels[e$method[keep]]), e[[col]][keep])),
          collapse = "\n")
  }

  mk <- function(cv) {
    is_roc <- identical(cv, "ROC")
    d  <- pooled %>% filter(curve == cv) %>%
      mutate(method = factor(method, levels = models), hz = lab_of(horizon))
    bd <- if (is.null(band)) NULL else band %>% filter(curve == cv) %>% mutate(hz = lab_of(horizon))
    p <- ggplot(d, aes(x, y, colour = method))
    if (!is.null(bd) && nrow(bd))
      p <- p + geom_ribbon(data = bd, mapping = aes(x = x, ymin = lo, ymax = hi),
                           fill = sc$values[[FEATURED]], alpha = 0.13, colour = NA,
                           inherit.aes = FALSE)
    p <- p + if (is_roc)
      geom_abline(slope = 1, intercept = 0, linetype = "22", colour = FAINT, linewidth = 0.4)
    else
      geom_hline(data = base_rate %>% mutate(hz = lab_of(horizon)),
                 aes(yintercept = br), linetype = "22", colour = FAINT, linewidth = 0.4)
    ann <- bind_rows(lapply(hz, function(h)
      tibble(hz = lab_of(h), lab = .auc_lab(h, if (is_roc) "auc_roc" else "auc_pr"))))
    # ROC: straight segments in the given (threshold) order, so the drawn area is the trapezoid
    # sum. geom_path, not geom_line: geom_line re-sorts on x, which leaves a vertical run's order
    # to the sort's tie handling. PR: "vh" steps, so each recall interval is drawn at the
    # precision of the point that reaches it -- the average-precision sum (see .s5_drawn()).
    p + (if (is_roc) geom_path(linewidth = 0.8) else geom_step(direction = "vh", linewidth = 0.8)) +
      facet_wrap(~ hz, nrow = 1) +
      geom_text(data = ann, aes(x = 0.98, y = if (is_roc) 0.04 else 0.98, label = lab),
                hjust = 1, vjust = if (is_roc) 0 else 1, size = 2.5, colour = INK,
                lineheight = 1.05, inherit.aes = FALSE) +
      scale_colour_manual(values = sc$values, labels = sc$labels, breaks = models, name = NULL) +
      # Three x labels, not five: the two facets sit side by side and "100%" on the left panel
      # was printing into "0%" on the right.
      scale_x_continuous(labels = percent_format(1), limits = c(0, 1), breaks = c(0, 0.5, 1),
                         expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_continuous(labels = percent_format(1), limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                         expand = expansion(mult = c(0.02, 0.02))) +
      labs(x = if (is_roc) "False-positive rate" else "Recall (share of invasions caught)",
           y = if (is_roc) "True-positive rate" else "Precision (share of monitored zones invaded)") +
      # Two columns: five entries in one row ran ~9 in, past the edge of an 8.6 in figure, and
      # the last keys were cut off. Measured below.
      guides(colour = guide_legend(ncol = 2, byrow = TRUE)) +
      coord_equal() +
      theme_pub(9.5) + theme(legend.position = "top",
                             panel.spacing.x = unit(14, "pt"))
  }

  pA <- mk("ROC"); pB <- mk("PR")
  FIG_W <- 8.6
  .assert_legend_fits(list(pA), 1, FIG_W, labels = "S5", who = "figS5")
  save_dual(pA, "FS5A_roc", 7.2, 4.2)
  save_dual(pB, "FS5B_pr",  7.2, 4.2)
  fig <- (pA / pB) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(legend.position = "top")
  save_dual(fig, "FigureS5_roc_pr_curves", FIG_W, 9.0, dir = FIG_DIR)
  invisible(fig)
}

# -----------------------------------------------------------------------------
# FIGURE S6 — MCMC convergence diagnostics
# -----------------------------------------------------------------------------
# Every model in the suite is a Stan fit, and a published posterior is only worth reading if
# its sampler converged. This is the whole diagnostic set for every model in one panel, rather
# than a maximum quoted in prose: R-hat against the 1.01 and 1.05 conventions, the minimum
# bulk and tail effective sample size against the 400-draw rule of thumb (Vehtari et al. 2021,
# Bayesian Analysis 16:667-718), and the divergent-transition count. Models are ordered by
# R-hat so the worst case is the first row a reader's eye lands on.
build_figS6 <- function() {
  dg <- .read_or_null(file.path(OUT, "reports", "bayes_convergence_diagnostics.csv"),
                      "convergence diagnostics")
  if (is.null(dg) || !nrow(dg)) stop("no convergence diagnostics to draw")
  dg <- dg %>%
    mutate(label = model_pretty_label(model, warn_unknown = FALSE),
           is_featured = model == FEATURED) %>%
    arrange(desc(rhat_max)) %>%
    mutate(label = factor(label, levels = rev(unique(label))))
  ptcol <- c(`TRUE` = "#D55E00", `FALSE` = "#0072B2")
  # The threshold labels sit below the first row (y = 0.3), so the lower expansion makes room for
  # them; the SAME y scale on all three panels keeps their rows aligned.
  .s6_y <- scale_y_discrete(expand = expansion(add = c(1.2, 0.6)))

  pA <- ggplot(dg, aes(rhat_max, label, colour = is_featured)) +
    geom_vline(xintercept = 1.01, linetype = "22", colour = FAINT, linewidth = 0.4) +
    geom_vline(xintercept = 1.05, linetype = "12", colour = "#C0392B", linewidth = 0.4) +
    geom_point(size = 1.8) +
    annotate("text", x = 1.01, y = 0.3, label = "1.01", size = 2.4, colour = MUTED, hjust = -0.15) +
    scale_colour_manual(values = ptcol, guide = "none") +
    scale_x_continuous(breaks = pretty_breaks(4)) + .s6_y +
    labs(x = "Maximum R-hat over the population and process parameters", y = NULL) +
    theme_pub(9) + theme(panel.grid.major.y = element_blank(),
                         axis.text.y = element_text(size = 6.8, colour = INK))

  ess <- dg %>%
    select(label, is_featured, `Bulk` = ess_bulk_min, `Tail` = ess_tail_min) %>%
    pivot_longer(c(Bulk, Tail), names_to = "kind", values_to = "ess")
  pB <- ggplot(ess, aes(ess, label, colour = kind)) +
    geom_vline(xintercept = 400, linetype = "22", colour = "#C0392B", linewidth = 0.4) +
    geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
    annotate("text", x = 400, y = 0.3, label = "400", size = 2.4, colour = MUTED, hjust = -0.15) +
    scale_colour_manual(values = c(Bulk = "#0072B2", Tail = "#E69F00"), name = NULL) +
    scale_x_continuous(breaks = pretty_breaks(4)) + .s6_y +
    labs(x = "Minimum effective sample size", y = NULL) +
    theme_pub(9) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank(),
                         legend.position = "top")

  pC <- ggplot(dg, aes(pmax(pct_divergent, 0), label, colour = is_featured)) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = ptcol, guide = "none") +
    # A floor of 0.5% on the axis: with no divergences anywhere the breaks would all fall inside
    # one 0.1% rounding step and print "0.0% 0.0% 0.0%".
    scale_x_continuous(labels = percent_format(accuracy = 0.1, scale = 1),
                       limits = c(0, max(0.5, max(dg$pct_divergent, na.rm = TRUE))),
                       breaks = pretty_breaks(3)) + .s6_y +
    labs(x = "Divergent transitions\n(% of post-warmup draws)", y = NULL) +
    theme_pub(9) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank())

  save_dual(pA, "FS6A_rhat", 4.6, 5.4)
  fig <- (pA | pB | pC) + plot_layout(widths = c(1.55, 1, 1)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureS6_mcmc_diagnostics", 12.0, 5.6, dir = FIG_DIR)
  invisible(fig)
}

# -----------------------------------------------------------------------------
# FIGURE S7 — calibration of the featured model over the forecast rounds
# -----------------------------------------------------------------------------
# Discrimination says the model ranks the right zones; calibration says its probabilities mean
# what they claim. A model can rank perfectly and still be useless for planning if 5% means 15%.
#
# THE PANELS ARE ORDERED BY HOW MUCH THEY CAN CARRY. Invasion is rare and the predicted
# probabilities are hugely skewed, so a POOLED reliability diagram — the conventional lead for a
# calibration figure — rests on very little here: it is panel C, last, and its bins had to be
# balanced on EXPECTED events, sum(p), to be worth drawing at all (RELIABILITY_N_BINS and
# RELIABILITY_BINNING, 00_config.R; balancing on OBSERVED events made the bin edges depend on
# the outcome and biased observed/predicted upward). What the record
# does support is the per-round view, because each round is a genuine out-of-sample forecast:
#   A  COUNTS per round: expected invasions (the sum of that round's probabilities) against the
#      number that happened. This is the quantity an operational reader cares about, it is what
#      the recalibration factor is fitted on, and the scatter of its points IS cross-validated
#      uncertainty — one point per held-out round, not one binomial interval.
#   B  The same comparison as a ratio, on a log scale so over- and under-prediction are
#      symmetric. Flat at 1 = calibrated.
#   C  RELIABILITY, pooled: within each predicted-probability band, what fraction was invaded?
#      Each bin carries BOTH a binomial (Wilson) interval and a between-origin one.
#
# LOG y ON A, AND THE WARM-UP SHADED. Two things used to make this figure mislead:
#   * the first origin's expectation is an order of magnitude above every other round's, so on a
#     linear axis the twelve rounds that matter were squeezed into the bottom quarter;
#   * the recalibration factor is fitted PREQUENTIALLY, on folds whose outcome has already
#     closed, so the earliest origins have no correction at all (delta = 1, not estimable) and
#     the first origin that does have one fits it on a single closed fold. The corrected
#     expectation therefore collapses and recovers over the first few rounds, and the O/E hump
#     that produces is the correction finding its level, NOT the model's calibration drifting.
#     Both scales are drawn and the not-yet-estimable rounds are shaded, so the two cannot be
#     confused.
#
# The interval on A and B is a Poisson interval on the OBSERVED count given the expected one:
# with a handful of invasions per round the sampling noise dominates, and a series of bare
# points invites reading noise as miscalibration.
build_figS7_one <- function(h) {
  d <- .scored(FEATURED, h)
  if (!nrow(d)) stop("no scored rows for the featured model")
  hh <- as.integer(h)

  # --- A/B: expected vs observed per round -----------------------------------
  # `eval_reliable` is written by run_invasion_lfo(); an LFO frame produced before that column
  # existed has no reliability label, and every round is then drawn as settled — which is what
  # it meant at the time, since that fold window admitted only settled rounds.
  if (!"eval_reliable" %in% names(d)) d$eval_reliable <- TRUE
  # p_invasion is the SELECTED scale (fs_apply_lfo_scale); p_invasion_raw is its uncorrected
  # twin. When no correction is in force the two coincide and the raw series is not drawn.
  if (!"p_invasion_raw" %in% names(d)) d$p_invasion_raw <- d$p_invasion
  if (!"delta_estimable" %in% names(d)) d$delta_estimable <- TRUE
  per <- d %>% group_by(cutoff = as.Date(cutoff)) %>%
    summarise(expected = sum(p_invasion), expected_raw = sum(p_invasion_raw),
              observed = sum(is_new_invasion), reliable = all(eval_reliable %in% TRUE),
              corrected = any(delta_estimable %in% TRUE), .groups = "drop") %>%
    arrange(cutoff)
  # Poisson interval on the OBSERVED count (exact, via the chi-square relation); a round with
  # zero invasions still has an upper limit, which is why the interval is drawn rather than a
  # bare point.
  per <- per %>% mutate(
    obs_lo = ifelse(observed == 0, NA_real_, stats::qchisq(0.05, 2 * observed) / 2),
    obs_hi = stats::qchisq(0.95, 2 * (observed + 1)) / 2,
    oe     = ifelse(expected > 0, observed / expected, NA_real_),
    oe_lo  = ifelse(expected > 0, obs_lo / expected, NA_real_),
    oe_hi  = ifelse(expected > 0, obs_hi / expected, NA_real_))
  # Rounds whose outcome window has not settled are drawn hollow: their observed count can
  # still rise, so a low O/E there is not yet evidence of over-prediction.
  SHP <- c(`TRUE` = 16, `FALSE` = 21)
  # Label rounds by the FORECAST ORIGIN (cutoff + 6), not the training week's start.
  per$origin <- lfo_origin(per$cutoff)
  raw_differs <- any(abs(per$expected_raw - per$expected) > 1e-8)
  # The stretch of origins with no fitted correction yet, as a shaded band behind both panels.
  warm <- per$origin[!per$corrected]
  # ymin = 0, NOT -Inf: both panels have a log10 y axis, which maps 0 to -Inf (drawn to the panel
  # edge) but turns -Inf into NaN, and ggplot then drops the band without a word.
  .warmup <- function() if (!length(warm)) NULL else list(
    ggplot2::annotate("rect", xmin = min(per$origin) - 3, xmax = max(warm) + 3,
                      ymin = 0, ymax = Inf, fill = FAINT, alpha = 0.18),
    ggplot2::annotate("text", x = min(per$origin) - 1, y = Inf, hjust = 0, vjust = 1.5,
                      label = "no correction\nfitted yet", size = 2.3, colour = MUTED,
                      lineheight = 0.95))
  LAB_E <- "Expected (sum of forecast probabilities)"
  LAB_R <- "Expected, before recalibration"
  LAB_O <- "Observed invasions"
  cols  <- c("#0072B2", FAINT, "#B33005"); names(cols) <- c(LAB_E, LAB_R, LAB_O)
  # LOG y: the expected count spans more than an order of magnitude across the rounds while the
  # observed count stays within a narrow range, so a linear axis would be set by one point. A zero observed count has no position on a log axis and is
  # drawn at the axis floor carrying its Poisson upper limit.
  y_lo <- min(c(per$expected, per$expected_raw, per$observed[per$observed > 0],
                per$obs_lo), na.rm = TRUE) / 1.6
  pA <- ggplot2::ggplot(per, ggplot2::aes(origin)) +
    .warmup() +
    ggplot2::geom_linerange(ggplot2::aes(ymin = pmax(obs_lo, y_lo), ymax = obs_hi),
                            colour = "#B33005", alpha = 0.35, linewidth = 0.5, na.rm = TRUE) +
    (if (raw_differs) ggplot2::geom_line(ggplot2::aes(y = expected_raw, colour = LAB_R),
                                         linewidth = 0.6, linetype = "22")) +
    ggplot2::geom_line(ggplot2::aes(y = expected, colour = LAB_E), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = pmax(observed, y_lo), colour = LAB_O, shape = reliable),
                        size = 2.1, fill = "white", stroke = 0.7) +
    ggplot2::scale_colour_manual(values = cols, name = NULL,
                                 breaks = c(LAB_E, LAB_R, LAB_O)[c(TRUE, raw_differs, TRUE)],
                                 guide = ggplot2::guide_legend(nrow = 3)) +
    ggplot2::scale_shape_manual(values = SHP, guide = "none") +
    ggplot2::scale_x_date(breaks = .thin_si_breaks(per$origin), date_labels = "%d %b",
                          expand = ggplot2::expansion(mult = c(0.06, 0.04))) +
    ggplot2::scale_y_continuous(trans = "log10", limits = c(y_lo, NA)) +
    ggplot2::labs(x = "Forecast origin (as-of date)",
                  y = sprintf("Invasions in the %d-week window", hh)) +
    theme_pub(9.5) + ggplot2::theme(legend.position = "top", legend.justification = "left",
                                    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  pB <- ggplot2::ggplot(per, ggplot2::aes(origin, oe)) +
    .warmup() +
    ggplot2::geom_hline(yintercept = 1, linetype = "22", colour = FAINT, linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = oe_lo, ymax = oe_hi), colour = PT_BLUE,
                            alpha = 0.45, linewidth = 0.5, na.rm = TRUE) +
    ggplot2::geom_line(colour = PT_BLUE, linewidth = 0.6, alpha = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = reliable), colour = PT_BLUE, fill = "white",
                        size = 2.1, stroke = 0.7) +
    ggplot2::scale_shape_manual(values = SHP, guide = "none") +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::scale_x_date(breaks = .thin_si_breaks(per$origin), date_labels = "%d %b",
                          expand = ggplot2::expansion(mult = c(0.06, 0.04))) +
    ggplot2::labs(x = "Forecast origin (as-of date)", y = "Observed / expected invasions") +
    theme_pub(9.5) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  # --- C: reliability, bins balanced on expected events, log-log -------------
  # Same estimator and same binning as the published table (23_prospective_eval.R), so the
  # panel and diagnostics/forecast_reliability_h*.csv cannot describe different curves.
  brk <- reliability_breaks(d$p_invasion, d$is_new_invasion,
                            RELIABILITY_N_BINS, RELIABILITY_BINNING)
  rel <- reliability_curve(d$p_invasion, as.integer(d$is_new_invasion),
                           fold = d$fold_id, breaks = brk, conf = RELIABILITY_CONF)
  # A BIN WITH NO INVASION has an observed fraction of exactly 0, which has no position on a
  # log axis. Those bins are not dropped — dropping them would hide precisely the bins where
  # the model predicted something and nothing happened. They are drawn at the axis floor as a
  # downward triangle carrying their upper confidence limit, i.e. "at most this".
  flo <- min(c(rel$mean_pred, rel$obs_freq[rel$k > 0], rel$obs_lo_cv[rel$obs_lo_cv > 0]),
             na.rm = TRUE) / 3
  top <- max(c(rel$mean_pred, rel$obs_hi, rel$obs_hi_cv), na.rm = TRUE) * 1.6
  rel <- rel %>% mutate(zero = k == 0L, obs_plot = ifelse(zero, flo, obs_freq),
                        lo_p = pmax(obs_lo, flo), lo_cv_p = pmax(obs_lo_cv, flo))
  pC <- ggplot2::ggplot(rel, ggplot2::aes(mean_pred, obs_plot)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22", colour = FAINT,
                         linewidth = 0.4) +
    # Thick = between-origin, thin drawn over it = binomial, both at RELIABILITY_CONF. The thick
    # bar is not always the
    # wider one: where origins agree more closely than independent sampling would predict the
    # design effect is below 1, which is worth seeing rather than hiding.
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo_cv_p, ymax = obs_hi_cv), colour = PT_BLUE,
                            alpha = 0.30, linewidth = 1.9, na.rm = TRUE) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo_p, ymax = obs_hi), colour = MUTED,
                            linewidth = 0.45, na.rm = TRUE) +
    ggplot2::geom_point(data = ~ dplyr::filter(.x, !zero), colour = PT_BLUE, size = 2) +
    ggplot2::geom_point(data = ~ dplyr::filter(.x, zero), colour = PT_BLUE, size = 2,
                        shape = 25, fill = "white", stroke = 0.7) +
    ggplot2::scale_x_continuous(trans = "log10", limits = c(flo, top),
                                labels = label_percent(accuracy = 0.001, drop0trailing = TRUE)) +
    ggplot2::scale_y_continuous(trans = "log10", limits = c(flo, top),
                                labels = label_percent(accuracy = 0.001, drop0trailing = TRUE)) +
    annotation_logticks(sides = "bl", colour = FAINT, size = 0.25,
                        short = unit(0.04, "cm"), mid = unit(0.07, "cm"), long = unit(0.11, "cm")) +
    ggplot2::labs(x = sprintf("Mean predicted invasion probability\n(%d bins balanced on expected events)",
                              nrow(rel)),
                  y = "Observed share invaded") +
    ggplot2::coord_equal() + theme_pub(9.5)

  save_dual(pA, sprintf("FS7A_counts_over_time_h%d", hh), 5.6, 4.6)
  save_dual(pB, sprintf("FS7B_oe_over_time_h%d", hh), 5.6, 4.2)
  save_dual(pC, sprintf("FS7C_reliability_h%d", hh), 4.4, 4.2)
  fig <- (pA | pB | pC) + plot_layout(widths = c(1.3, 1.25, 1)) +
    plot_annotation(tag_levels = "A")
  save_dual(fig, sprintf("FigureS7_calibration_h%d", hh), 13.6, 5.0, dir = FIG_DIR)
  invisible(fig)
}

#' Both horizons. The submitted figure was h=2 only, with nothing in its name saying so, while
#' the reliability diagram it now absorbs was published at both.
build_figS7 <- function() {
  for (.h in get0("LFO_HORIZONS", ifnotfound = c(1L, 2L))) build_figS7_one(.h)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# FIGURE S8 — generation-time sensitivity of the featured model
# -----------------------------------------------------------------------------
# The generation time is an ASSUMPTION: it enters the renewal convolution that builds the
# import force, so changing it changes the offset on every row and beta_0 absorbs the
# difference. The suite is composed at a single anchor precisely so the GT is not selected on,
# which leaves the honest question — how much does the assumption matter? — for these two
# panels: A, does it change out-of-sample skill; B, does it change today's probabilities.
GT_ORDER <- c("short", "medium", "long")
GT_COL   <- c(short = "#56B4E9", medium = "#0072B2", long = "#CC79A7")
#' The generation-time arms of THE FEATURED MODEL — or a refusal.
#'
#' WHY IT CAN REFUSE. The arms used to be refits of the PREVIOUS run's winner, whatever its
#' covariate set: on LINELIST_09102026 the featured model is Bayes-M14-fill-geo (ccvi + d_min)
#' and the arms on disk, Bayes-M14-fill-gtshort / -gtlong, are refits of the covariate-free
#' Bayes-M14-fill-med. Drawn as "short / medium / long" beside the featured model, the figure
#' would attribute a covariate difference to the generation time. Two checks, both required:
#'   1. THE LABEL. bayes_gt_sensitivity_specs() (21_bayesian_renewal.R) names an arm by
#'      appending "-gt<profile>" to bayes_gt_arm_stem() of the base label (00_config.R: a
#'      trailing "-med|-short|-long" token is stripped, before the form suffix of a
#'      time-varying form), so Bayes-M14-fill-geo -> Bayes-M14-fill-geo-gtshort,
#'      Bayes-M14-fill-med -> Bayes-M14-fill-gtshort and Bayes-M14-fill-med-tvfrailty ->
#'      Bayes-M14-fill-tvfrailty-gtshort. Both forms — the stem's, and the featured label with
#'      the suffix appended — are accepted (they coincide unless the label carries one of
#'      those tokens); finding both for one profile is ambiguous and refused.
#'   2. THE SPECIFICATION. A label is a promise, not evidence. The arm must carry the featured
#'      model's mobility kernel (mobility_kernel_from_method()) and EXACTLY its covariate terms,
#'      read off the fitted coefficients in reports/bayes_parameters.csv. A model whose
#'      coefficients are not on record cannot be verified and is refused.
#' @param methods the methods present in the evaluation table.
#' @return named c(medium =, short =, long =) of the arms present (medium = FEATURED).
.gt_arms <- function(methods) {
  stem <- bayes_gt_arm_stem(FEATURED)   # 00_config.R: the rule bayes_gt_sensitivity_specs() uses
  arms <- c(medium = FEATURED)
  for (pr in c("short", "long")) {
    cand <- unique(c(paste0(stem, "-gt", pr), paste0(FEATURED, "-gt", pr)))
    hit  <- intersect(cand, methods)
    if (length(hit) > 1L)
      stop("two ", pr, "-generation-time arms of ", FEATURED, " are scored (",
           paste(hit, collapse = ", "), "); which one refits the featured model is ambiguous")
    if (length(hit)) arms[[pr]] <- hit
  }
  if (length(arms) < 2L) {
    other <- grep("-gt(short|long)$", methods, value = TRUE)
    stop("no generation-time arm of the featured model ", FEATURED, " is scored",
         if (length(other)) paste0("; the arms present (", paste(other, collapse = ", "),
                                   ") are refits of a different specification") else
           " (INCLUDE_GT_SENSITIVITY_MODELS off, or the pipeline has not been re-run)")
  }
  # 2. the specification: kernel and covariate terms must be the featured model's own.
  pf <- file.path(OUT, "reports", "bayes_parameters.csv")
  if (!file.exists(pf))
    stop("reports/bayes_parameters.csv is absent, so the generation-time arms cannot be ",
         "verified as refits of ", FEATURED)
  pp <- suppressWarnings(read_csv(pf, show_col_types = FALSE))
  terms_of <- function(m) {
    if (!any(pp$model == m))
      stop("bayes_parameters.csv holds no coefficients for ", m, "; it cannot be verified ",
           "as a refit of ", FEATURED)
    sort(unique(pp$term[pp$model == m & !(pp$is_intercept %in% TRUE)]))
  }
  k0 <- mobility_kernel_from_method(FEATURED); t0 <- terms_of(FEATURED)
  for (a in arms[names(arms) != "medium"]) {
    ka <- mobility_kernel_from_method(a); ta <- terms_of(a)
    if (!identical(ka, k0) || !identical(ta, t0))
      stop(sprintf(paste0("%s is not a refit of the featured model %s: kernel %s vs %s, ",
                          "covariates {%s} vs {%s}. Re-run the pipeline so the arms are built ",
                          "on this run's selection."),
                   a, FEATURED, ka, k0, paste(ta, collapse = ", "), paste(t0, collapse = ", ")))
  }
  arms
}

build_figS8 <- function() {
  arms <- .gt_arms(unique(ev$method))

  # --- A: cross-validated performance under each GT ---
  # Four metrics on one panel, each on its own scale, so they are faceted with a free x rather
  # than forced onto a common axis that would mean nothing.
  # THE SCALE OF EACH METRIC IS STATED, AND MATCHES FIGURE 2A. Discrimination (AUC-PR skill,
  # AUC-ROC) is quoted from the raw columns, as everywhere in the suite: it is invariant to the
  # within-round recalibration. The log score and calibration-in-the-large DO depend on the
  # scale, and this panel used to quote them RAW while Figure 2A quotes the featured model's
  # calibration RECALIBRATED — so one model carried two different calibration numbers in two
  # figures. They now follow the column this pass draws on (PSCALE_COL), as Figure 2A's do,
  # and the facet says which.
  .recal <- identical(PSCALE_COL, "p_recal")
  .ls  <- if (.recal) "log_score_recal" else "log_score"
  .cal <- if (.recal) "calibration_in_large_recal" else "calibration_in_large"
  .sw  <- if (.recal) "recalibrated" else "raw"
  if (length(setdiff(c(.ls, .cal), names(ev))))
    stop("invasion_evaluation.csv lacks ", paste(setdiff(c(.ls, .cal), names(ev)), collapse = ", "))
  # AUC-PR skill PER ROUND over all rounds (auc_pr_skill_date), as every other published AUC-PR
  # skill; all four metrics are over all rounds, and their labels say so.
  if (!"auc_pr_skill_date" %in% names(ev))
    stop("invasion_evaluation.csv lacks auc_pr_skill_date; re-run the pipeline.")
  MET <- setNames(c("AUC-PR skill per round, all rounds (x base rate)", "AUC-ROC, all rounds",
                    sprintf("Log score, %s, all rounds (lower is better)", .sw),
                    sprintf("Calibration, %s, all rounds: predicted / observed", .sw)),
                  c("auc_pr_skill_date", "auc_roc", .ls, .cal))
  a <- ev %>% filter(method %in% arms) %>%
    mutate(gt = factor(names(arms)[match(method, arms)], levels = GT_ORDER)) %>%
    select(gt, horizon, all_of(names(MET))) %>%
    pivot_longer(all_of(names(MET)), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(MET[metric], levels = unname(MET)),
           hz = unname(HZ_LAB[as.character(horizon)]))
  pA <- ggplot(a, aes(value, gt, colour = gt, shape = hz)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
    facet_wrap(~ metric, scales = "free_x", nrow = 1,
               labeller = labeller(metric = label_wrap_gen(24))) +
    scale_colour_manual(values = GT_COL, guide = "none") +
    scale_shape_manual(values = c(16, 1), name = NULL) +
    scale_x_continuous(breaks = pretty_breaks(3), expand = expansion(mult = c(0.12, 0.12))) +
    labs(x = NULL, y = "Generation-time profile") +
    theme_pub(9.5) + theme(legend.position = "top")

  # --- B: today's invasion probabilities under each GT ---
  rs <- read_csv(fs_risk_csv(OUT), show_col_types = FALSE)
  pr <- readRDS(file.path(OUT, "forecasts", "bayes_current_predictions.rds"))
  H  <- 2L
  cur <- pr %>% filter(method %in% arms, horizon == H, !(was_active_before %in% TRUE),
                       is.finite(p_invasion)) %>%
    mutate(gt = factor(names(arms)[match(method, arms)], levels = GT_ORDER))
  if (!nrow(cur)) stop("no current predictions for the generation-time arms")
  top <- cur %>% filter(gt == "medium") %>% slice_max(p_invasion, n = 20, with_ties = FALSE) %>%
    pull(health_zone)
  cur <- cur %>% filter(health_zone %in% top) %>%
    mutate(health_zone = factor(health_zone,
      levels = cur %>% filter(gt == "medium", health_zone %in% top) %>%
        arrange(p_invasion) %>% pull(health_zone)))
  # THE AXIS NAMES ITS SCALE. bayes_current_predictions.rds carries no stamp, so it is derived:
  # the risk table's prob_scale says whether this run deployed recalibration at all, and the
  # recalibration table whether EVERY arm has a fitted factor at this horizon (the suite
  # rescales a member only by its own delta). Arms on different scales are refused, not labelled.
  .ps_rs <- unique(as.character(rs$prob_scale))
  .rc <- .read_or_null(file.path(OUT, "diagnostics", "invasion_recalibration.csv"),
                       "recalibration factors")
  .has_d <- if (is.null(.rc)) rep(FALSE, length(arms))
            else unname(arms) %in% .rc$method[.rc$horizon == H & is.finite(.rc$delta)]
  .scaleB <- if (identical(.ps_rs, "raw") || !any(.has_d)) "raw"
             else if (all(.has_d) && length(.ps_rs) == 1L && grepl("^recalibrated", .ps_rs)) "recalibrated"
             else stop("[si] FigureS8B: the generation-time arms are not all on one probability ",
                       "scale (risk table: ", paste(.ps_rs, collapse = "/"), "; arms with a factor: ",
                       sum(.has_d), " of ", length(arms), ").", call. = FALSE)
  pB <- ggplot(cur, aes(p_invasion, health_zone, colour = gt)) +
    geom_linerange(aes(xmin = p_lo, xmax = p_hi), position = position_dodge(width = 0.6),
                   linewidth = 0.45, alpha = 0.55) +
    geom_point(size = 1.9, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = GT_COL, name = "Generation time") +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.04))) +
    labs(x = sprintf("P(first confirmed case within %d weeks), %s", H, .scaleB), y = NULL) +
    theme_pub(9.5) + theme(panel.grid.major.y = element_blank(), legend.position = "top",
                           axis.text.y = element_text(size = 7.6, colour = INK))

  save_dual(pA, "FS8A_gt_skill", 9.6, 3.0)
  save_dual(pB, "FS8B_gt_current", 5.2, 5.6)
  fig <- (pA / pB) + plot_layout(heights = c(0.62, 1)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureS8_generation_time_sensitivity", 10.0, 9.4, dir = FIG_DIR)
  invisible(fig)
}

# -----------------------------------------------------------------------------
# FIGURE S9 — time-varying import coefficient: the forms of the stage-1 specification
# -----------------------------------------------------------------------------
# DRAWN FROM THE SELECTION STAGE'S OUTPUTS (run_all.R, invasion_tv_form_comparison() in
# 16_invasion_eval.R); nothing is computed here. Stage 2 of the selection refits the stage-1
# specification under three forms (RW1 drift, gamma frailty, RW1 + frailty), cross-validates
# them on the same rounds and recalibrates each by the single-delta rule. Fixed beta and the
# SELECTABLE forms compete on the selection composite; the frailty forms are a secondary analysis,
# never candidates (TV_FORMS `selectable`, 00_config.R), and are labelled so. This figure shows
# how all four compare.
#
# Panel A is calibration over time: cumulative expected events (recalibrated where a factor was
# fitted; the earliest post-burn-in rounds have none yet and are drawn as open points, raw)
# against cumulative observed, post-burn-in rounds, per form and horizon. An event is a positive
# zone-round: at 2 weeks the windows of consecutive rounds overlap, so one invasion can count
# twice. A calibrated form tracks the observed step; a fixed beta drifts above it as the
# outbreak matures.
# Panel B is the horizon-mean paired log-score gain over fixed beta with its 90% zone-bootstrap
# interval, and the smaller of the two Bonferroni lower bounds (zone bootstrap; origin-level
# Diebold-Mariano with the HLN correction) -- the former guard, REPORTED: it does not select.
TV_LAB <- c(fixed = "Fixed beta", rw1 = "Random-walk drift", frailty = "Gamma frailty",
            rw1_frailty = "Random walk + frailty")
TV_COL <- c(fixed = "grey35", rw1 = "#D55E00", frailty = "#0072B2", rw1_frailty = "#009E73")

build_figS9 <- function() {
  bd <- .read_or_null(file.path(OUT, "key_outputs", "tv_comparison_by_date.csv"), "tv comparison by date")
  gd <- .read_or_null(file.path(OUT, "key_outputs", "tv_comparison_guard.csv"), "tv comparison guard")
  if (is.null(bd) || is.null(gd))
    stop("no time-varying comparison for this run (INCLUDE_TV_FORMS off, or the forms' ",
         "back-test failed or did not cover the shared rounds)")
  if (!identical(unique(bd$featured_model), FEATURED) || !identical(unique(gd$featured_model), FEATURED))
    stop(sprintf("the tv_comparison tables describe %s, not the featured model %s; re-run run_all.R",
                 paste(unique(c(bd$featured_model, gd$featured_model)), collapse = ", "), FEATURED))
  if (!"delta_estimable" %in% names(bd))
    stop("tv_comparison_by_date.csv has no delta_estimable column, so the rounds without a ",
         "fitted recalibration factor cannot be marked; re-run run_all.R")
  if (length(setdiff(unique(c(bd$form, gd$form)), names(TV_LAB))))
    stop("the tv_comparison tables carry form(s) with no label: ",
         paste(setdiff(unique(c(bd$form, gd$form)), names(TV_LAB)), collapse = ", "))
  a <- bd %>% filter(!burnin) %>% arrange(form, horizon, cutoff) %>%
    group_by(form, horizon) %>%
    mutate(cum_obs = cumsum(observed), cum_exp = cumsum(expected_recal)) %>% ungroup() %>%
    mutate(form = factor(form, levels = names(TV_LAB)),
           hz = factor(unname(HZ_LAB[as.character(horizon)]), levels = unname(HZ_LAB)),
           # FORECAST ORIGIN, not the fold cutoff: the comparison keys its rows by the training
           # week's start, and every other round axis in the suite is labelled by the date the
           # forecast was issued (lfo_origin() = cutoff + LFO_ORIGIN_OFFSET_DAYS). Plotting the
           # cutoff put this panel six days to the left of the same rounds elsewhere.
           origin = lfo_origin(as.Date(cutoff)))
  obs <- a %>% filter(form == "fixed")
  # The featured model is one of the four (fixed beta when stage 2 kept the stage-1 model).
  .fkey <- tv_form_of(FEATURED)
  .fkey <- if (is.na(.fkey)) "fixed" else TV_FORMS[[.fkey]]$form
  tv_lab <- TV_LAB; tv_lab[[.fkey]] <- paste0(tv_lab[[.fkey]], " (featured)")
  # The forms that can never be featured (the frailty ones) say so.
  .sec <- vapply(Filter(function(f) !isTRUE(f$selectable), TV_FORMS), `[[`, "", "form")
  tv_lab[.sec] <- paste0(tv_lab[.sec], " (secondary)")
  pA <- ggplot(a, aes(origin, cum_exp, colour = form)) +
    geom_step(data = obs, aes(origin, cum_obs), inherit.aes = FALSE, colour = INK, linewidth = 0.9) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = ifelse(delta_estimable %in% TRUE, "rec", "unc")), size = 1.4) +
    facet_wrap(~ hz, scales = "free_y") +
    scale_colour_manual(values = TV_COL, labels = tv_lab, name = NULL) +
    scale_shape_manual(values = c(rec = 16, unc = 1),
                       labels = c(rec = "recalibrated", unc = "uncorrected (no factor fitted yet)"),
                       name = NULL) +
    scale_x_date(labels = function(d) paste(format(d, "%d"), month.abb[as.integer(format(d, "%m"))]),
                 breaks = pretty_breaks(5)) +
    labs(x = "Forecast origin (as-of date)",
         y = paste0("Cumulative positive zone-rounds\n(black step = observed; lines = expected;\n",
                    "2-week windows overlap)")) +
    # The form key and the recalibration key on rows of their own: on one row the second ran off
    # the right edge ("uncorrected (no factor fitted ...").
    guides(colour = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 2, nrow = 1)) +
    theme_pub(9.5) + theme(legend.position = "top", legend.box = "vertical",
                           legend.box.just = "left")
  # THE ORIGIN-LEVEL BOUND IS THE DIEBOLD-MARIANO ONE, AND ONLY THAT (origin_lo_bonferroni,
  # origin_lo90 / origin_hi90). The circular origin-block bootstrap it replaced was
  # ANTI-CONSERVATIVE -- by simulation, > 0 under the null in 8-13% of replicates against 1.67%
  # nominal -- so an old table that still carries block_lo_bonferroni is not drawn. No fallback.
  .need_g <- c("zone_lo_bonferroni", "origin_lo_bonferroni", "zone_lo90", "zone_hi90")
  if (length(setdiff(.need_g, names(gd))))
    stop("tv_comparison_guard.csv lacks ", paste(setdiff(.need_g, names(gd)), collapse = ", "),
         ". A table without origin_lo_bonferroni predates the origin-level Diebold-Mariano ",
         "bound (its block_* columns are the retired, anti-conservative block bootstrap); ",
         "re-run run_all.R")
  g <- gd %>% filter(forecasts == "p_recal") %>%
    mutate(form = factor(form, levels = names(TV_LAB)),
           # the smaller of the two Bonferroni bounds (the former guard; reported only)
           bonf = pmin(zone_lo_bonferroni, origin_lo_bonferroni))
  pB <- ggplot(g, aes(gain_x1000, form, colour = form)) +
    geom_vline(xintercept = 0, colour = FAINT) +
    geom_errorbar(aes(xmin = zone_lo90, xmax = zone_hi90), width = 0.18, linewidth = 0.7,
                  orientation = "y") +
    geom_point(size = 2.6) +
    geom_point(aes(x = bonf), shape = 124, size = 4) +
    scale_colour_manual(values = TV_COL, guide = "none") +
    scale_y_discrete(labels = tv_lab, limits = rev) +
    labs(x = paste0("Log-score gain over fixed beta, recalibrated forecasts (x1000; > 0 = better)\n",
                    "Bar: 90% zone-bootstrap interval; tick: Bonferroni lower bound, the smaller of ",
                    "the zone-bootstrap and\norigin-level Diebold-Mariano bounds"),
         y = NULL) +
    theme_pub(9.5) + theme(panel.grid.major.y = element_blank())
  save_dual(pA, "FS9A_tv_calibration", 9.6, 3.6)
  save_dual(pB, "FS9B_tv_guard", 9.6, 2.6)
  fig <- (pA / pB) + plot_layout(heights = c(1, 0.6)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureS9_time_varying_beta", 10.0, 7.2, dir = FIG_DIR)
  invisible(fig)
}

# -----------------------------------------------------------------------------
# FIGURE S13 — prior against posterior, featured model
# -----------------------------------------------------------------------------
# WHAT THE PRIORS ARE DOING HERE. The priors are not decoration: 20_forecast_detail.R records
# that they "regularise the separation that diverged the frequentist joint fit", and the
# featured model is estimated from a few tens of invasion events over several thousand at-risk
# zone-weeks against a strong mobility offset. A reader is then entitled to ask the obvious question — is the
# posterior the data speaking, or the prior? — and only a picture of BOTH densities answers
# it. The published parameter table gives a median and a 90% interval, which cannot: an
# interval that is narrow because the prior is narrow looks exactly like one that is narrow
# because the data are informative.
#
# READ IT AS CONTRACTION, PANEL BY PANEL. Where the posterior is much narrower than the prior
# and sits somewhere the prior did not especially favour, the data have moved it. Where the
# two curves are close, the parameter is prior-dominated and should be quoted as such.
#
# ONE MODEL, THE FEATURED ONE. The suite-wide version of this plot (bayes_posterior_densities,
# 20_forecast_detail.R) stacks every member as ridgelines and answers a different question —
# how the estimate moves across kernels. This figure is about whether the published model's
# own parameters are identified.
#
# THE PRIOR IS READ FROM THE ARTEFACT, NOT RE-DECLARED. bayes_posterior_draws() stamps each
# row with the brms prior_string it was fitted under, so this panel cannot go on drawing
# Normal(-3,2) against a posterior fitted under something else. An un-stamped or unparseable
# prior draws the posterior alone rather than inventing one.

#' A prior_string ("normal(-3, 2)") as a density with respect to log10(value).
#'
#' BOTH CURVES MUST BE DENSITIES OF THE SAME VARIABLE. The panel's x axis is log10, and the
#' posterior is drawn as a kernel density of log10(hr) (below), so the prior has to be the
#' density of log10 of the parameter too — not the density of the parameter itself. They are
#' not the same curve: they differ by the Jacobian dx/dl = x*ln(10), which is not constant, so
#' plotting one against the other misstates the prior's SHAPE and not merely its height. For
#' the featured model's intercept the difference is stark — the correct curve is a symmetric
#' bell centred on exp(-3) = 0.050, the incorrect one decays monotonically across the whole
#' panel, because a LogNormal(-3, 2) density in x has its mode away at exp(-3 - 2^2) = 0.0009.
#'
#' The priors are stated on the coefficient (log) scale, so if log(value) ~ Normal(mu, sigma)
#' then log10(value) ~ Normal(mu/ln10, sigma/ln10). Being normal in the plotted variable also
#' guarantees mode = median, which the window rule below relies on.
#'
#' Returns NULL for anything not recognised — no prior stamped, several conflicting ones, or a
#' family this parser does not implement — and the caller then draws the posterior alone. It
#' never guesses: a figure captioned "prior" that shows the wrong prior is worse than one that
#' shows none.
.prior_log10_density <- function(spec) {
  if (is.null(spec)) return(NULL)
  spec <- unique(spec[!is.na(spec) & nzchar(trimws(spec))])
  if (length(spec) != 1L) return(NULL)
  # "half-normal(0, s)": the RW1 innovation SD, a normal prior restricted to sigma > 0 (Stan's
  # <lower=0>); its density on log10(exp(sigma)) is twice the normal's on l >= 0 and 0 below.
  m <- regmatches(spec, regexec(
    "^\\s*(half-)?normal\\s*\\(\\s*([^,()]+?)\\s*,\\s*([^,()]+?)\\s*\\)\\s*$", spec))[[1]]
  if (length(m) != 4L) return(NULL)
  half <- nzchar(m[2])
  mu <- suppressWarnings(as.numeric(m[3])); sg <- suppressWarnings(as.numeric(m[4]))
  if (!is.finite(mu) || !is.finite(sg) || sg <= 0) return(NULL)
  if (half && mu != 0) return(NULL)        # only the zero-centred truncation is a half-normal
  L <- log(10)
  if (half)
    return(list(d = function(l) ifelse(l >= 0, 2 * stats::dnorm(l, mean = 0, sd = sg / L), 0),
                q = function(p) stats::qnorm((1 + p) / 2, mean = 0, sd = sg / L)))
  list(d = function(l) stats::dnorm(l, mean = mu / L, sd = sg / L),
       q = function(p) stats::qnorm(p, mean = mu / L, sd = sg / L))
}

build_figS13 <- function() {
  dr <- .read_or_null(file.path(OUT, "key_outputs", "bayes_posterior_draws.csv"),
                      "Bayesian posterior draws")
  if (is.null(dr) || !nrow(dr))
    stop("no posterior draws on disk — run the pipeline, or 46_bayes_posterior_draws.R")
  if (!"prior" %in% names(dr))
    stop("bayes_posterior_draws.csv predates the prior stamp (no `prior` column); rebuild it ",
         "with 46_bayes_posterior_draws.R so the prior drawn is the prior fitted")
  d <- dr %>% filter(model == FEATURED)
  if (!nrow(d))
    stop(sprintf("bayes_posterior_draws.csv holds no draws for the featured model (%s); it has %s",
                 FEATURED, paste(unique(dr$model), collapse = ", ")))

  # PANEL ORDER: the import coefficient first — it is the model's scale parameter, and the one
  # reference line below does not apply to it — then the covariates in the order the draws
  # carry them, which is the brms coefficient order and therefore the order
  # bayes_parameters.csv lists them in too (both are built from the same fit).
  terms <- c(unique(d$term[d$is_intercept]), unique(d$term[!d$is_intercept]))
  # make.unique() so two terms that happen to share a .PARAM_LABEL entry cannot collide into
  # duplicate factor levels, which is an error rather than a cosmetic problem.
  # A row whose effect_scale is not a coefficient scale says what it is (an RW1 form's intercept
  # is the walk's STARTING level, not the import coefficient in force).
  .esc <- function(t) { e <- unique(d$effect_scale[d$term == t]); e[1] }
  .rw1_start <- get0("RW1_START_LEVEL_SCALE", ifnotfound = NA_character_)
  lab_of <- setNames(make.unique(vapply(terms, function(t) {
    l <- .lk(.PARAM_LABEL, t, t); e <- .esc(t)
    if (!is.na(e) && identical(e, .rw1_start)) "RW1 start level, exp(b0)"
    else if (t %in% c("frailty_theta", "rw1_sigma_week")) l   # the label already says what it is
    else if (!is.na(e) && !e %in% c("hazard ratio", "odds ratio")) paste0(l, " (", e, ")") else l
  }, character(1))), terms)
  d <- d %>% mutate(param = factor(unname(lab_of[term]), levels = unname(lab_of[terms])))

  # -- one panel at a time: window, then both densities ON THAT WINDOW ----------------------
  #
  # THE DENSITIES ARE COMPUTED HERE, NOT BY stat_density(). Two reasons, both correctness
  # rather than taste. (1) Under scale_x_log10() ggplot transforms x BEFORE the stat runs, so
  # stat_density() would silently return a density of log10(hr) while the prior curve was a
  # density of hr — the mismatch described above. Computing both explicitly makes the shared
  # measure visible in the code. (2) stat_density() extends ~3 bandwidths past the data and
  # scales by the maximum over that extended range, and with free x scales there is no way to
  # clip it per panel; the panel would then be sized by the kernel's tail rather than by the
  # window, and `scaled` would be relative to a peak that may sit off-panel.
  #
  # WINDOW: the posterior's 0.1-99.9 percentile range, extended to include the prior's
  # interquartile range so the prior's peak is always on the panel (the prior is normal in
  # log10, so its mode IS its median and the IQR brackets it), then padded by 6% of the span.
  # It is deliberately NOT the prior's full support: the intercept prior spans four orders of
  # magnitude against a posterior spanning less than one, and an axis sized to it would render
  # the posterior as a featureless spike — which is the half of the comparison a reader has to
  # judge. A prior whose tails run off the panel is itself the finding, and the axis says so.
  curves <- purrr::map_dfr(terms, function(tm) {
    v <- d$hr[d$term == tm]
    v <- v[is.finite(v) & v > 0]           # hr = exp(coef) > 0; guards a signif() underflow
    if (length(v) < 2L) return(NULL)
    lv <- log10(v)
    pd <- .prior_log10_density(d$prior[d$term == tm])
    rng <- as.numeric(stats::quantile(lv, c(0.001, 0.999), names = FALSE))
    if (!is.null(pd)) rng <- range(c(rng, pd$q(c(0.25, 0.75))))
    pad <- max(diff(rng) * 0.06, 1e-3)     # 1e-3 decades keeps a degenerate posterior drawable
    lo <- rng[1] - pad; hi <- rng[2] + pad
    # bw.SJ is the better bandwidth here (the posteriors are smooth and unimodal) but it can
    # fail outright on a degenerate sample, so it degrades to Silverman's rule rather than
    # taking the figure down with it.
    bw <- tryCatch(stats::bw.SJ(lv), error = function(e) NA_real_)
    if (!is.finite(bw) || bw <= 0) bw <- stats::bw.nrd0(lv)
    if (!is.finite(bw) || bw <= 0) return(NULL)
    kd <- stats::density(lv, bw = bw, from = lo, to = hi, n = 512L)
    out <- tibble(term = tm, curve = "Posterior", l = kd$x, dens = kd$y)
    if (!is.null(pd)) {
      gl <- seq(lo, hi, length.out = 512L)
      out <- bind_rows(out, tibble(term = tm, curve = "Prior", l = gl, dens = pd$d(gl)))
    }
    # Each curve scaled to ITS OWN maximum over the window: the comparison the figure is for
    # is shape and width, not integrated mass, and on a common absolute scale a prior spread
    # over four orders of magnitude is numerically invisible beside the posterior. A curve
    # whose mass is entirely outside the window has max 0 and is dropped rather than divided
    # by zero.
    # No helper column, so no select() call: MASS::select() masks dplyr's wherever MASS gets
    # attached (it already bites sibling scripts in this tree), and it would not fail loudly
    # here — it would fail inside a figure that had otherwise built.
    out %>% group_by(curve) %>%
      mutate(dens = { .mx <- max(dens)
                      if (is.finite(.mx) && .mx > 0) dens / .mx else NA_real_ }) %>%
      ungroup() %>% filter(is.finite(dens))
  })
  if (!nrow(curves)) stop("no drawable density for any parameter of ", FEATURED)
  # bayes_posterior_draws() carries the coefficients only, so a featured model with process
  # parameters — a time-varying form contributes frailty_theta and/or rw1_sigma_week to
  # bayes_parameters.csv — has parameters this figure cannot show. Say so, rather than let the
  # figure quietly present a subset of the model as though it were the model.
  .published_terms <- local({
    pf <- file.path(OUT, "reports", "bayes_parameters.csv")
    if (!file.exists(pf)) return(character(0))
    pp <- suppressWarnings(read_csv(pf, show_col_types = FALSE))
    pp$term[pp$model == FEATURED]
  })
  # beta_last_training_week is DERIVED (exp(b0 + u_S), an RW1 form's level in force), not a
  # parameter with a prior, so it has nothing to show here; any other published parameter
  # without a curve would make the figure present a subset of the model as the model.
  .absent <- setdiff(.published_terms, c(unique(curves$term), "beta_last_training_week"))
  if (length(.absent))
    stop("[si] FigureS13 would omit ", paste(.absent, collapse = ", "), " of ", FEATURED,
         ": bayes_posterior_draws.csv has no drawable posterior for it. Rebuild the draws ",
         "(46_bayes_posterior_draws.R, or the pipeline).", call. = FALSE)
  curves <- curves %>%
    mutate(param = factor(unname(lab_of[term]), levels = unname(lab_of[terms])),
           x = 10^l)

  # Posterior medians, and the HR = 1 reference. Keyed on is_intercept rather than on the
  # first factor level: "no effect" is a statement about a RATIO, so it means nothing for
  # beta0, and inferring which panel is the intercept from panel order would silently put the
  # line in the wrong place the day the ordering rule changes.
  meta <- d %>% group_by(term, param, is_intercept) %>%
    summarise(med = stats::median(hr), .groups = "drop") %>%
    filter(param %in% unique(curves$param))

  # The ratio panels: the covariates (hazard / odds ratios). Process parameters and the import
  # coefficient are not ratios, so "no effect = 1" does not apply to them.
  .ratio_terms <- unique(d$term[d$effect_scale %in% c("hazard ratio", "odds ratio") & !d$is_intercept])
  .ratio_meta <- meta %>% filter(term %in% .ratio_terms)
  .lt_keys <- c("Posterior median", if (nrow(.ratio_meta)) "No effect (ratio = 1)")
  PRIOR_COL <- "grey55"; POST_COL <- "#0072B2"
  .scale_word <- .lk(c(`hazard ratio` = "hazard ratio", `odds ratio` = "odds ratio"),
                     (d$effect_scale[!d$is_intercept][1]) %|NA|% "hazard ratio", "effect")
  p <- ggplot() +
    geom_area(data = curves %>% filter(curve == "Prior"),
              aes(x, dens, fill = "Prior"), colour = NA, alpha = 0.30) +
    geom_line(data = curves %>% filter(curve == "Prior"),
              aes(x, dens, colour = "Prior"), linewidth = 0.5, linetype = "22") +
    geom_area(data = curves %>% filter(curve == "Posterior"),
              aes(x, dens, fill = "Posterior"), colour = NA, alpha = 0.55) +
    geom_line(data = curves %>% filter(curve == "Posterior"),
              aes(x, dens, colour = "Posterior"), linewidth = 0.8) +
    # BOTH REFERENCE LINES ARE IN THE LEGEND. They were drawn bare, which left two unexplained
    # verticals in a figure whose caption lives in the manuscript — and the "= 1" line is
    # absent from the beta0 panel, so its meaning could not even be inferred by analogy.
    # The ratio = 1 line only on the coefficient (ratio) panels, and only when one exists: a
    # zero-row layer would carry no facet variable and be drawn in EVERY panel, beta0's included.
    (if (nrow(.ratio_meta)) geom_vline(data = .ratio_meta,
               aes(xintercept = 1, linetype = "No effect (ratio = 1)"),
               colour = MUTED, linewidth = 0.35)) +
    geom_vline(data = meta, aes(xintercept = med, linetype = "Posterior median"),
               colour = POST_COL, linewidth = 0.4) +
    scale_linetype_manual(values = c(`Posterior median` = "11",
                                     `No effect (ratio = 1)` = "42")[.lt_keys], name = NULL,
                          breaks = .lt_keys) +
    # order: the densities are what the figure is about, so they lead; the reference lines
    # follow. Without it ggplot orders the guides by scale type and puts linetype first.
    guides(fill     = guide_legend(order = 1),
           colour   = guide_legend(order = 1),
           linetype = guide_legend(order = 2,
                                   override.aes = list(colour = c(`Posterior median` = POST_COL,
                                                                  `No effect (ratio = 1)` = MUTED)[.lt_keys]))) +
    facet_wrap(~ param, scales = "free_x", nrow = 1,
               labeller = labeller(param = label_wrap_gen(width = 30))) +
    scale_x_log10(labels = label_number(drop0trailing = TRUE),
                  expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0, 1.04), breaks = c(0, 0.5, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    # Legend keys only for curves actually drawn: with no parsable prior anywhere, a forced
    # "Prior" break would advertise a curve the figure does not contain.
    scale_fill_manual(values = c(Prior = PRIOR_COL, Posterior = POST_COL), name = NULL,
                      breaks = intersect(c("Prior", "Posterior"), unique(curves$curve))) +
    scale_colour_manual(values = c(Prior = PRIOR_COL, Posterior = POST_COL), name = NULL,
                        breaks = intersect(c("Prior", "Posterior"), unique(curves$curve))) +
    labs(x = paste0("Log scale: ",
                    if (any(d$effect_scale %in% .rw1_start)) "exp(b0), the RW1 walk's start level"
                    else "beta0 = exp(intercept), the import coefficient",
                    if (length(.ratio_terms)) sprintf("; %s per +1 SD for covariates", .scale_word) else "",
                    if (any(d$term == "frailty_theta")) "; frailty variance theta" else "",
                    if (any(d$term == "rw1_sigma_week")) "; exp(RW1 innovation SD)" else ""),
         y = "Relative density (peak = 1)") +
    theme_pub(9.5) +
    theme(legend.position = "top",
          panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
          panel.spacing.x = unit(10, "pt"))

  # Width from the panels actually DRAWN, not from the factor's levels: facet_wrap drops unused
  # levels, so a parameter whose density could not be computed would otherwise leave the figure
  # sized for a panel that is not there.
  .w <- max(6.4, 2.9 * dplyr::n_distinct(curves$param))
  save_dual(p, "FigureS13_prior_posterior", .w, 3.1, dir = FIG_DIR)
  invisible(p)
}

# -----------------------------------------------------------------------------
if (!isTRUE(get0(".SI_MODEL_FIGURES_NO_RUN", ifnotfound = FALSE))) {
  # Every file each builder writes, composite first (see .clear_figure_files()).
  .build("FigureS5",  build_figS5,  c("FigureS5_roc_pr_curves", "FS5A_roc", "FS5B_pr"))
  .build("FigureS6",  build_figS6,  c("FigureS6_mcmc_diagnostics", "FS6A_rhat"))
  .build("FigureS7",  build_figS7,  c("FigureS7_calibration", "FS7A_counts_over_time",
                                      "FS7B_oe_over_time", "FS7C_reliability"))
  .build("FigureS8",  build_figS8,  c("FigureS8_generation_time_sensitivity", "FS8A_gt_skill",
                                      "FS8B_gt_current"))
  .build("FigureS9",  build_figS9,  c("FigureS9_time_varying_beta", "FS9A_tv_calibration",
                                      "FS9B_tv_guard"))
  .build("FigureS13", build_figS13, c("FigureS13_prior_posterior"))
  message(sprintf("\n[si-model] %s  ->  %s",
                  paste(sprintf("%s:%s", names(ok), ifelse(unlist(ok), "ok", "FAILED")),
                        collapse = "  "), FIG_DIR))
  if (any(!unlist(ok)))
    stop(sum(!unlist(ok)), " of ", length(ok), " SI model figure(s) failed: ",
         paste(names(ok)[!unlist(ok)], collapse = ", "),
         ". Their previous files were removed, so nothing stale is left in their place; see ",
         "the warnings above for each cause.", call. = FALSE)
}

}   # end of the script-run guard

