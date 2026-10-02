# =============================================================================
# make_si_model_figures.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
# Supplementary MODEL figures, in the main-text house style (figure_style.R):
# no embedded titles or subtitles, panel tags only, captions live in the manuscript.
#
#   FigureS5  Cross-validated discrimination   A) ROC curves, B) precision-recall curves,
#                                              featured model vs the three structural
#                                              baselines, both horizons, pooled across the
#                                              leave-future-out folds with the per-fold
#                                              curves drawn faintly behind.
#   FigureS6  MCMC convergence                 A) max R-hat per model, B) minimum bulk and
#                                              tail effective sample size, C) divergent
#                                              transitions — every fitted model, one row each.
#   FigureS7  Calibration over time            A) reliability curve, B) observed vs expected
#                                              invasions per round, C) the O/E ratio over
#                                              the rounds, for the featured model.
#   FigureS8  Generation-time sensitivity      A) cross-validated skill of the featured model
#                                              under the short / medium / long generation
#                                              time, B) today's invasion probabilities under
#                                              each, for the current top-ranked zones.
#   FigureS9  Time-varying beta_t              A) posterior trajectory of the import
#                                              coefficient under each beta_t process, with
#                                              the fixed-beta model as the reference,
#                                              B) what each process costs or buys in
#                                              cross-validated skill.
#
# EVERYTHING IS READ, NOTHING IS REFITTED. The curves come from the saved cross-validation
# (lfo_results.rds), the diagnostics from bayes_convergence_diagnostics.csv, the beta_t
# trajectories from bayes_beta_trajectory.csv and the probabilities from the published risk
# tables. A figure that recomputed its own model would be free to disagree with the pipeline.
#
# Outputs -> outputs/key_outputs/manuscript_figures/{FigureS5..S9}.{pdf,png} + panels/
# Run:  Rscript make_si_model_figures.R   (from anywhere; paths are anchored with here())
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(scales); library(forcats)
})
options(dplyr.summarise.inform = FALSE)

HERE <- file.path(here::here(), "spatiotemporal")
OUT  <- file.path(HERE, "outputs")
source(file.path(HERE, "00_config.R"))
source(file.path(HERE, "forecast_scale.R"))
FIG_DIR   <- fs_out_dir(file.path(OUT, "key_outputs", "manuscript_figures"))
PANEL_DIR <- fs_out_dir(file.path(OUT, "key_outputs", "manuscript_figures", "panels"))
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
source(file.path(HERE, "figure_style.R"))

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
fs_apply_lfo_scale(lfo)

# FEATURED MODEL: READ, DO NOT RECOMPUTE — the pipeline's own record of what it picked.
FEATURED <- local({
  f <- file.path(OUT, "key_outputs", "model_selection.json")
  if (!file.exists(f))
    stop("[si] model_selection.json is missing; these figures name the featured model and ",
         "must not guess it. Run the pipeline first.", call. = FALSE)
  sel <- jsonlite::fromJSON(f, simplifyVector = TRUE)
  m <- sel$featured$bayesian$method
  if (is.null(m) || !length(m) || is.na(m[1])) m <- sel$featured$headline$method
  if (is.null(m) || !length(m) || is.na(m[1]))
    stop("[si] model_selection.json names no featured model.", call. = FALSE)
  as.character(m[1])
})
message(sprintf("[si] featured model: %s (%s)", FEATURED, model_pretty_label(FEATURED)))

BASELINES <- intersect(c("Gravity-B4", "Distance-B1", "Adjacency-B7"), unique(lfo$method))
HZ_LAB    <- c(`1` = "1 week ahead", `2` = "2 weeks ahead")
# Featured model first (blue), then the structural baselines in the order they are listed.
MODEL_COL <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_len(1L + length(BASELINES))],
                      c(FEATURED, BASELINES))

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
ok <- list()
.build <- function(nm, f) {
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

#' Evaluate a monotone-x step curve on a common grid, so curves from different rounds can be
#' compared pointwise.
#'
#' Both curves are STEP functions of their x axis, so the value at grid point g is the curve's
#' value at the LAST operating point with x <= g (constant interpolation, held from the left).
#' Linear interpolation would be wrong for a precision-recall curve in particular: the
#' achievable region between two operating points is not the straight line between them
#' (Davis & Goadrich 2006, ICML), and the constant hold is the same convention the tie-aware
#' average precision above already uses.
#' Grid points before the curve's first x have no defined value and return NA, which keeps
#' them out of the quantiles rather than inventing a value at recall no round reached.
.on_grid <- function(cur, grid) {
  if (is.null(cur) || !nrow(cur)) return(rep(NA_real_, length(grid)))
  cur <- cur[order(cur$x), , drop = FALSE]
  idx <- findInterval(grid, cur$x)
  out <- rep(NA_real_, length(grid))
  ok <- idx >= 1L
  out[ok] <- cur$y[idx[ok]]
  out
}

#' The rows one (method, horizon) is actually SCORED on: at-risk only, finite probability,
#' finite outcome. Every panel in this file filters identically, so a curve, a calibration bin
#' and a quoted AUC all describe the same set of rows.
.scored <- function(m, h) {
  lfo %>% filter(method == m, horizon == h, is.finite(p_invasion),
                 is.finite(is_new_invasion), !(was_active_before %in% TRUE))
}

build_figS5 <- function() {
  models <- c(FEATURED, BASELINES)
  hz     <- sort(unique(lfo$horizon))
  GRID   <- seq(0, 1, by = 0.005)
  pooled <- list(); band <- list()
  for (cv in c("ROC", "PR")) {
    pts <- if (identical(cv, "ROC")) .roc_points else .pr_points
    for (h in hz) {
      for (m in models) {
        d <- .scored(m, h); if (!nrow(d)) next
        cur <- pts(d$p_invasion, d$is_new_invasion); if (is.null(cur)) next
        pooled[[length(pooled) + 1L]] <- mutate(cur, method = m, horizon = h, curve = cv)
      }
      # Fold-variation band, featured model only.
      dF <- .scored(FEATURED, h); if (!nrow(dF)) next
      per <- lapply(sort(unique(dF$fold_id)), function(fd) {
        df <- dF %>% filter(fold_id == fd)
        .on_grid(pts(df$p_invasion, df$is_new_invasion), GRID)
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

  # Reference lines: the ROC chance diagonal, and for PR the pooled base rate — the precision
  # a random watch-list achieves. Computed PER HORIZON from the rows actually scored, not
  # assumed, because the h=2 base rate is not the h=1 one (overlapping outcome windows).
  base_rate <- bind_rows(lapply(hz, function(h)
    tibble(horizon = h, br = mean(.scored(FEATURED, h)$is_new_invasion))))

  lab_of <- function(v) unname(HZ_LAB[as.character(v)])
  # The AUCs are QUOTED FROM invasion_evaluation.csv, not recomputed from the drawn curve, so
  # the number on the panel is the published one. A model with no row is omitted rather than
  # printed as NA.
  .auc_lab <- function(h, col) {
    e <- ev %>% filter(horizon == h) %>% slice(match(models, method))
    keep <- !is.na(e$method) & is.finite(e[[col]])
    if (!any(keep)) return("")
    paste(sprintf("%s  %.3f", model_pretty_label(e$method[keep], warn_unknown = FALSE),
                  e[[col]][keep]), collapse = "\n")
  }

  mk <- function(cv) {
    is_roc <- identical(cv, "ROC")
    d  <- pooled %>% filter(curve == cv) %>%
      mutate(method = factor(method, levels = models), hz = lab_of(horizon))
    bd <- if (is.null(band)) NULL else band %>% filter(curve == cv) %>% mutate(hz = lab_of(horizon))
    p <- ggplot(d, aes(x, y, colour = method))
    if (!is.null(bd) && nrow(bd))
      p <- p + geom_ribbon(data = bd, mapping = aes(x = x, ymin = lo, ymax = hi),
                           fill = MODEL_COL[[FEATURED]], alpha = 0.13, colour = NA,
                           inherit.aes = FALSE)
    p <- p + if (is_roc)
      geom_abline(slope = 1, intercept = 0, linetype = "22", colour = FAINT, linewidth = 0.4)
    else
      geom_hline(data = base_rate %>% mutate(hz = lab_of(horizon)),
                 aes(yintercept = br), linetype = "22", colour = FAINT, linewidth = 0.4)
    ann <- bind_rows(lapply(hz, function(h)
      tibble(hz = lab_of(h), lab = .auc_lab(h, if (is_roc) "auc_roc" else "auc_pr"))))
    p + geom_step(linewidth = 0.8) +
      facet_wrap(~ hz, nrow = 1) +
      geom_text(data = ann, aes(x = 0.98, y = if (is_roc) 0.04 else 0.98, label = lab),
                hjust = 1, vjust = if (is_roc) 0 else 1, size = 2.5, colour = INK,
                lineheight = 1.05, inherit.aes = FALSE) +
      scale_colour_manual(values = MODEL_COL, name = NULL,
                          labels = function(v) model_pretty_label(v, warn_unknown = FALSE)) +
      # Three x labels, not five: the two facets sit side by side and "100%" on the left panel
      # was printing into "0%" on the right.
      scale_x_continuous(labels = percent_format(1), limits = c(0, 1), breaks = c(0, 0.5, 1),
                         expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_continuous(labels = percent_format(1), limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                         expand = expansion(mult = c(0.02, 0.02))) +
      labs(x = if (is_roc) "False-positive rate" else "Recall (share of invasions caught)",
           y = if (is_roc) "True-positive rate" else "Precision (share of monitored zones invaded)") +
      coord_equal() +
      theme_pub(9.5) + theme(legend.position = "top",
                             panel.spacing.x = unit(14, "pt"))
  }

  pA <- mk("ROC"); pB <- mk("PR")
  save_dual(pA, "FS5A_roc", 7.2, 4.2)
  save_dual(pB, "FS5B_pr",  7.2, 4.2)
  fig <- (pA / pB) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(legend.position = "top")
  save_dual(fig, "FigureS5_roc_pr_curves", 8.6, 9.0, dir = FIG_DIR)
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

  pA <- ggplot(dg, aes(rhat_max, label, colour = is_featured)) +
    geom_vline(xintercept = 1.01, linetype = "22", colour = FAINT, linewidth = 0.4) +
    geom_vline(xintercept = 1.05, linetype = "12", colour = "#C0392B", linewidth = 0.4) +
    geom_point(size = 1.8) +
    annotate("text", x = 1.01, y = 0.3, label = "1.01", size = 2.4, colour = MUTED, hjust = -0.15) +
    scale_colour_manual(values = ptcol, guide = "none") +
    scale_x_continuous(breaks = pretty_breaks(4)) +
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
    scale_x_continuous(breaks = pretty_breaks(4)) +
    labs(x = "Minimum effective sample size", y = NULL) +
    theme_pub(9) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank(),
                         legend.position = "top")

  pC <- ggplot(dg, aes(pmax(pct_divergent, 0), label, colour = is_featured)) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = ptcol, guide = "none") +
    scale_x_continuous(labels = percent_format(accuracy = 0.1, scale = 1),
                       breaks = pretty_breaks(3)) +
    labs(x = "Divergent transitions (% of post-warmup draws)", y = NULL) +
    theme_pub(9) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank())

  save_dual(pA, "FS6A_rhat", 4.6, 5.4)
  fig <- (pA | pB | pC) + plot_layout(widths = c(1.55, 1, 1)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureS6_mcmc_diagnostics", 12.0, 5.6, dir = FIG_DIR)
  invisible(fig)
}

# -----------------------------------------------------------------------------
# FIGURE S7 — calibration of the featured model, and how it moves over the rounds
# -----------------------------------------------------------------------------
# Discrimination says the model ranks the right zones; calibration says its probabilities mean
# what they claim. A model can rank perfectly and still be useless for planning if 5% means
# 15%. The three panels answer three different versions of the question:
#   A  RELIABILITY, pooled: within each predicted-probability bin, what fraction was invaded?
#      Bins are equal-COUNT (quantiles of the predicted probability), not equal-width: the
#      predictions are heavily right-skewed, so equal-width bins put almost every row in the
#      first bin and leave the interesting range with a handful of observations each.
#   B  COUNTS per round: the expected number of invasions (the sum of the probabilities) beside
#      the number that happened. This is the quantity an operational reader cares about, and
#      it is the one the recalibration factor is fitted on.
#   C  The O/E ratio per round, on a log scale so over- and under-prediction are symmetric.
#      A flat series at 1 is a calibrated model; a trend is a beta_0 that is drifting, which is
#      exactly what the time-varying models in Figure S9 are built to absorb.
#
# The interval on B and C is a Poisson interval on the OBSERVED count given the expected one:
# with a handful of invasions per round the sampling noise dominates, and a series of bare
# points invites reading noise as miscalibration.
build_figS7 <- function(h = 2L) {
  d <- .scored(FEATURED, h)
  if (!nrow(d)) stop("no scored rows for the featured model")

  # --- A: reliability, equal-count bins, LOG-LOG ---
  # WHY LOG-LOG. The predicted probabilities span four orders of magnitude and are extremely
  # right-skewed: on linear axes eight of ten equal-count bins pile onto the origin and the
  # panel shows two distinguishable points. Log axes spread the bins over the range the model
  # actually uses, and the 1:1 line is still a straight line, so the reading ("points on the
  # line = calibrated") is unchanged.
  # A BIN WITH NO INVASION has an observed fraction of exactly 0, which has no position on a
  # log axis. Those bins are not dropped — dropping them would hide precisely the bins where
  # the model predicted something and nothing happened. They are drawn at the axis floor as a
  # downward triangle carrying their upper confidence limit, i.e. "at most this".
  NB <- 10L
  rel <- d %>%
    mutate(bin = cut(rank(p_invasion, ties.method = "first"), breaks = NB, labels = FALSE)) %>%
    group_by(bin) %>%
    summarise(p_hat = mean(p_invasion), obs = mean(is_new_invasion), n = n(),
              k = sum(is_new_invasion), .groups = "drop") %>%
    # Clopper-Pearson (exact) interval on the observed fraction: the bins are near zero, where
    # a Wald interval runs below it.
    mutate(lo = mapply(function(k, n) stats::binom.test(k, n)$conf.int[1], k, n),
           hi = mapply(function(k, n) stats::binom.test(k, n)$conf.int[2], k, n),
           zero = k == 0L)
  # Axis floor: one decade below the smallest quantity that has to be shown.
  flo <- min(c(rel$p_hat, rel$obs[!rel$zero], rel$hi[rel$zero]), na.rm = TRUE) / 3
  top <- max(c(rel$p_hat, rel$hi), na.rm = TRUE) * 1.6
  rel <- rel %>% mutate(obs_plot = ifelse(zero, flo, obs),
                        lo_plot  = pmax(lo, flo))
  pA <- ggplot(rel, aes(p_hat, obs_plot)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = FAINT, linewidth = 0.4) +
    geom_linerange(aes(ymin = lo_plot, ymax = hi), colour = PT_BLUE, alpha = 0.5, linewidth = 0.5) +
    geom_point(data = ~ dplyr::filter(.x, !zero), colour = PT_BLUE, size = 2) +
    geom_point(data = ~ dplyr::filter(.x, zero), colour = PT_BLUE, size = 2, shape = 25,
               fill = "white", stroke = 0.7) +
    scale_x_continuous(trans = "log10", limits = c(flo, top),
                       labels = label_percent(accuracy = 0.001, drop0trailing = TRUE)) +
    scale_y_continuous(trans = "log10", limits = c(flo, top),
                       labels = label_percent(accuracy = 0.001, drop0trailing = TRUE)) +
    annotation_logticks(sides = "bl", colour = FAINT, size = 0.25,
                        short = unit(0.04, "cm"), mid = unit(0.07, "cm"), long = unit(0.11, "cm")) +
    labs(x = sprintf("Mean predicted invasion probability (%d equal-count bins)", NB),
         y = "Observed share invaded") +
    coord_equal() +
    theme_pub(9.5)

  # --- B/C: expected vs observed per round ---
  # `eval_reliable` is written by run_invasion_lfo(); an LFO frame produced before that column
  # existed has no reliability label, and every round is then drawn as settled — which is what
  # it meant at the time, since that fold window admitted only settled rounds.
  if (!"eval_reliable" %in% names(d)) d$eval_reliable <- TRUE
  per <- d %>% group_by(cutoff = as.Date(cutoff)) %>%
    summarise(expected = sum(p_invasion), observed = sum(is_new_invasion),
              reliable = all(eval_reliable %in% TRUE), .groups = "drop") %>%
    arrange(cutoff)
  # Poisson interval on the OBSERVED count (exact, via the chi-square relation); a round with
  # zero invasions still has an upper limit, which is why the interval is drawn rather than a
  # bare point.
  per <- per %>% mutate(
    obs_lo = ifelse(observed == 0, 0, stats::qchisq(0.05, 2 * observed) / 2),
    obs_hi = stats::qchisq(0.95, 2 * (observed + 1)) / 2,
    oe     = ifelse(expected > 0, observed / expected, NA_real_),
    oe_lo  = ifelse(expected > 0, obs_lo / expected, NA_real_),
    oe_hi  = ifelse(expected > 0, obs_hi / expected, NA_real_))
  # Rounds whose outcome window has not settled are drawn hollow: their observed count can
  # still rise, so a low O/E there is not yet evidence of over-prediction.
  SHP <- c(`TRUE` = 16, `FALSE` = 21)

  # Label rounds by the FORECAST ORIGIN (cutoff + 6), not the training week's start.
  per$origin <- lfo_origin(per$cutoff)
  pB <- ggplot(per, aes(origin)) +
    geom_linerange(aes(ymin = obs_lo, ymax = obs_hi), colour = "#B33005", alpha = 0.35,
                   linewidth = 0.5) +
    geom_line(aes(y = expected, colour = "Expected (sum of forecast probabilities)"), linewidth = 0.8) +
    geom_point(aes(y = observed, colour = "Observed invasions", shape = reliable),
               size = 2.1, fill = "white", stroke = 0.7) +
    scale_colour_manual(values = c("Expected (sum of forecast probabilities)" = "#0072B2",
                                   "Observed invasions" = "#B33005"), name = NULL) +
    scale_shape_manual(values = SHP, guide = "none") +
    scale_x_date(breaks = .thin_si_breaks(per$origin), date_labels = "%d %b",
                 expand = expansion(mult = c(0.03, 0.03))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = "Forecast origin (as-of date)", y = sprintf("Invasions in the %d-week window", h)) +
    theme_pub(9.5) + theme(legend.position = "top",
                           axis.text.x = element_text(angle = 30, hjust = 1))

  pC <- ggplot(per, aes(origin, oe)) +
    geom_hline(yintercept = 1, linetype = "22", colour = FAINT, linewidth = 0.4) +
    geom_linerange(aes(ymin = oe_lo, ymax = oe_hi), colour = PT_BLUE, alpha = 0.45, linewidth = 0.5) +
    geom_line(colour = PT_BLUE, linewidth = 0.6, alpha = 0.8) +
    geom_point(aes(shape = reliable), colour = PT_BLUE, fill = "white", size = 2.1, stroke = 0.7) +
    scale_shape_manual(values = SHP, guide = "none") +
    scale_y_continuous(trans = "log10", breaks = c(0.25, 0.5, 1, 2, 4),
                       labels = c("0.25", "0.5", "1", "2", "4")) +
    scale_x_date(breaks = .thin_si_breaks(per$origin), date_labels = "%d %b",
                 expand = expansion(mult = c(0.03, 0.03))) +
    labs(x = "Forecast origin (as-of date)",
         y = "Observed / expected invasions") +
    theme_pub(9.5) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  save_dual(pA, "FS7A_reliability", 4.4, 4.2)
  save_dual(pB, "FS7B_counts_over_time", 5.6, 4.2)
  save_dual(pC, "FS7C_oe_over_time", 5.6, 4.2)
  fig <- (pA | pB | pC) + plot_layout(widths = c(1, 1.25, 1.25)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureS7_calibration_over_time", 13.2, 4.6, dir = FIG_DIR)
  invisible(fig)
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
.gt_arms <- function(methods) {
  stem <- sub("-(med|short|long|gtshort|gtlong)$", "", FEATURED)
  arms <- c(medium = FEATURED,
            short  = paste0(stem, "-gtshort"),
            long   = paste0(stem, "-gtlong"))
  arms[arms %in% methods]
}

build_figS8 <- function() {
  arms <- .gt_arms(unique(ev$method))
  if (length(arms) < 2L)
    stop("the generation-time sensitivity arms are not in the evaluation table ",
         "(INCLUDE_GT_SENSITIVITY_MODELS off, or the pipeline has not been re-run)")

  # --- A: cross-validated performance under each GT ---
  # Four metrics on one panel, each on its own scale, so they are faceted with a free x rather
  # than forced onto a common axis that would mean nothing.
  MET <- c(auc_pr_skill        = "AUC-PR skill (x base rate)",
           auc_roc             = "AUC-ROC",
           log_score           = "Log score (lower is better)",
           calibration_in_large = "Calibration: predicted / observed")
  a <- ev %>% filter(method %in% arms) %>%
    mutate(gt = factor(names(arms)[match(method, arms)], levels = GT_ORDER)) %>%
    select(gt, horizon, all_of(names(MET))) %>%
    pivot_longer(all_of(names(MET)), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(MET[metric], levels = unname(MET)),
           hz = unname(HZ_LAB[as.character(horizon)]))
  pA <- ggplot(a, aes(value, gt, colour = gt, shape = hz)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
    facet_wrap(~ metric, scales = "free_x", nrow = 1) +
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
  pB <- ggplot(cur, aes(p_invasion, health_zone, colour = gt)) +
    geom_linerange(aes(xmin = p_lo, xmax = p_hi), position = position_dodge(width = 0.6),
                   linewidth = 0.45, alpha = 0.55) +
    geom_point(size = 1.9, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = GT_COL, name = "Generation time") +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.04))) +
    labs(x = sprintf("P(first confirmed case within %d weeks)", H), y = NULL) +
    theme_pub(9.5) + theme(panel.grid.major.y = element_blank(), legend.position = "top",
                           axis.text.y = element_text(size = 7.6, colour = INK))

  save_dual(pA, "FS8A_gt_skill", 9.6, 3.0)
  save_dual(pB, "FS8B_gt_current", 5.2, 5.6)
  fig <- (pA / pB) + plot_layout(heights = c(0.62, 1)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureS8_generation_time_sensitivity", 10.0, 9.4, dir = FIG_DIR)
  invisible(fig)
}

# -----------------------------------------------------------------------------
# FIGURE S9 — time-varying import coefficient beta_t
# -----------------------------------------------------------------------------
# Panel A is the posterior trajectory of beta_t under each process, with the fixed-beta model
# drawn as the flat reference. The dashed segment past the last training week is the FORECAST
# continuation each process actually uses, which is where the processes differ most: the
# random walk persists the current level, the weekly random effect reverts to beta_0, the
# AR(1) and GP decay toward it, the log-linear trend extrapolates its slope.
#
# Panel B prices them. A time-varying beta is a strictly more flexible model, so it cannot fit
# the training data worse; the only question that matters is whether it FORECASTS better, and
# that is the cross-validated skill on the identical folds. A variant that does not beat the
# fixed-beta model is evidence that beta has not drifted detectably — a result, not a failure.
# The processes actually fitted, from the config rather than a literal list here, so cutting
# TV_BETA_TYPES cannot leave this figure looking for arms that no longer exist. The labels
# and colours below stay complete: an entry for a retired process is harmless (the factor
# drops unused levels), an entry MISSING for a fitted one is not.
.tv_types <- function() as.character(get0("TV_BETA_TYPES",
                                          ifnotfound = c("trend", "week", "rw1", "ar1", "gp")))
TV_LAB <- c(none = "Fixed beta", trend = "Log-linear trend", week = "Weekly random effect",
            rw1 = "Random walk", ar1 = "AR(1)", gp = "Gaussian process")
TV_COL <- c(none = "grey35", trend = "#E69F00", week = "#56B4E9",
            rw1 = "#D55E00", ar1 = "#009E73", gp = "#CC79A7")

build_figS9 <- function() {
  tj <- .read_or_null(file.path(OUT, "key_outputs", "bayes_beta_trajectory.csv"),
                      "beta_t trajectories")
  if (is.null(tj) || !nrow(tj)) stop("no beta_t trajectories to draw")
  stem <- sub(sprintf("-(%s)$", paste0("tv", .tv_types(), collapse = "|")), "", FEATURED)
  fam  <- tj %>% filter(model == stem | startsWith(model, paste0(stem, "-tv")))
  if (!nrow(fam))
    stop(sprintf("bayes_beta_trajectory.csv holds no rows for %s or its time-varying variants", stem))
  fam <- fam %>% mutate(tv = factor(tv, levels = names(TV_LAB)),
                        week_date = as.Date(week_date))
  xcol <- if (all(is.na(fam$week_date))) "week" else "week_date"

  pA <- ggplot(fam, aes(.data[[xcol]], beta, colour = tv, fill = tv)) +
    geom_ribbon(data = fam %>% filter(!is_forecast), aes(ymin = beta_lo, ymax = beta_hi),
                alpha = 0.10, colour = NA) +
    geom_line(data = fam %>% filter(!is_forecast), linewidth = 0.85) +
    geom_line(data = fam %>% filter(is_forecast), linewidth = 0.85, linetype = "22") +
    geom_point(data = fam %>% filter(is_forecast), size = 1.4) +
    scale_colour_manual(values = TV_COL, labels = TV_LAB, name = NULL, drop = TRUE) +
    scale_fill_manual(values = TV_COL, guide = "none") +
    scale_y_continuous(trans = "log10", labels = label_number(accuracy = 0.001)) +
    labs(x = if (identical(xcol, "week_date")) "Week (start date)" else "Week index",
         y = expression(paste("Import coefficient  ", beta[t], "  (log scale)"))) +
    theme_pub(9.5) + theme(legend.position = "top")
  if (identical(xcol, "week_date"))
    pA <- pA + scale_x_date(date_labels = "%d %b", breaks = pretty_breaks(6)) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # --- B: cross-validated skill of each variant against the fixed-beta model ---
  fam_m <- c(stem, paste0(stem, "-tv", .tv_types()))
  e <- ev %>% filter(method %in% fam_m) %>%
    mutate(tv = factor(ifelse(method == stem, "none",
                              sub("^.*-tv", "", method)), levels = names(TV_LAB)),
           hz = unname(HZ_LAB[as.character(horizon)]))
  if (!nrow(e)) stop("the time-varying variants are not in the evaluation table")
  MET <- c(auc_pr_skill = "AUC-PR skill (x base rate)",
           mean_rank_of_truth = "Mean rank of invaded zones (lower is better)",
           log_score = "Log score (lower is better)")
  b <- e %>% select(tv, hz, all_of(names(MET))) %>%
    pivot_longer(all_of(names(MET)), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(MET[metric], levels = unname(MET)))
  pB <- ggplot(b, aes(value, tv, colour = tv, shape = hz)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
    facet_wrap(~ metric, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = TV_COL, guide = "none") +
    scale_y_discrete(labels = TV_LAB) +
    scale_shape_manual(values = c(16, 1), name = NULL) +
    scale_x_continuous(breaks = pretty_breaks(3), expand = expansion(mult = c(0.12, 0.12))) +
    labs(x = NULL, y = expression(paste(beta[t], " process"))) +
    theme_pub(9.5) + theme(legend.position = "top", panel.grid.major.y = element_blank())

  save_dual(pA, "FS9A_beta_trajectory", 6.4, 4.2)
  save_dual(pB, "FS9B_tv_skill", 9.6, 3.2)
  fig <- (pA / pB) + plot_layout(heights = c(1, 0.8)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureS9_time_varying_beta", 10.0, 8.4, dir = FIG_DIR)
  invisible(fig)
}

# -----------------------------------------------------------------------------
if (!isTRUE(get0(".SI_MODEL_FIGURES_NO_RUN", ifnotfound = FALSE))) {
  .build("FigureS5", build_figS5)
  .build("FigureS6", build_figS6)
  .build("FigureS7", build_figS7)
  .build("FigureS8", build_figS8)
  .build("FigureS9", build_figS9)
  message(sprintf("\n[si-model] %s  ->  %s",
                  paste(sprintf("%s:%s", names(ok), ifelse(unlist(ok), "ok", "FAILED")),
                        collapse = "  "), FIG_DIR))
}
