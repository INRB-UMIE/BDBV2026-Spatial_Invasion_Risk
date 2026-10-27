# =============================================================================
# 39_cascade_eval_figure.R — CASCADE evaluation figure, in the Figure 2 style
# BDBV 2026 DRC · 3-month spatial invasion cascade
#
# A faithful cascade analogue of make_publication_figures.R's Figure 2 (the
# short-horizon backtest evaluation), for the LONGER-TERM invasion cascade. The
# cascade is a single model, so where Figure 2 contrasts many models / two short
# horizons, this figure contrasts the leave-future-out backtest ORIGINS (folds)
# over the backtest window K (CASCADE_CALIB_K weeks; read from the backtest's own K
# column). Same design system, same panels, honest about the data ceiling:
#
#   A  Discrimination          AUC-PR skill (× base rate) per backtest origin.
#   B  Real-time prioritisation share of true invasions caught vs a random
#                              watch-list of the same size (capture curve).
#   C  Ranking accuracy        mean rank of the truly-invaded zones (lower better).
#   D  Per-round forecasts      top-K predicted zones per origin, coloured by the
#                              realised outcome (invaded within the K-week window, or not).
#
# All four panels come from ONE source — the cascade leave-future-out backtest
# (cascade_backtest → summary + per-zone `detail`) — so they are mutually
# consistent. NOT "fully out-of-sample": each origin's model is trained only on data known at
# that origin, but when run_cascade.R calibrates, delta is fitted to the pooled new-invasion
# COUNT on these same origins and K, so the probability LEVEL (the bar lengths in D, any
# predicted/observed ratio) is in-sample; the ranking metrics in A-C are what delta was not
# fitted to. Entry: cascade_fig_evaluation(ctx), called
# from run_cascade.R (needs ctx$backtest carrying its `detail` attribute).
# Standalone regen: Rscript spatiotemporal/39_cascade_eval_figure.R (from the RDS).
# Saves Figure2_cascade.{pdf,png} to cascade/figures ONLY (the published Figure2_cascade is
# make_manuscript_figure2_cascade.R's), plus key_outputs/cascade_backtest_per_origin.csv and
# the inputs RDS.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2)
  library(patchwork); library(scales)
})

# ---- design system (verbatim from make_publication_figures.R / Figure 2) ----
if (!exists("INK"))   INK   <- "grey15"
if (!exists("MUTED")) MUTED <- "grey38"
if (!exists("FAINT")) FAINT <- "grey72"
if (!exists("GRID"))  GRID  <- "grey92"
if (!exists("OKABE")) OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7",
                                 "#E69F00","#56B4E9","#F0E442","#000000")
.CASC_FAMILY <- "sans"
CASC_ACCENT  <- "#0072B2"                       # featured-cascade blue (Okabe)
OUT_COL <- c("Invaded (first case)" = "#D55E00", "Not invaded" = "grey75")

# theme_pub, plus concise visible titles (this figure stands alone — Figure 2
# suppresses titles because its captions live in the manuscript, not in the figure).
.theme_casc <- function(base = 12) {
  theme_minimal(base_size = base, base_family = .CASC_FAMILY) %+replace% theme(
    plot.title = element_text(size = base + 0.4, face = "bold", colour = INK,
                              hjust = 0, margin = margin(b = 1)),
    plot.subtitle = element_text(size = base - 2.2, colour = MUTED, hjust = 0,
                                 margin = margin(b = 5)),
    plot.caption = element_blank(),
    axis.title = element_text(size = base - 1.6, colour = MUTED),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4), angle = 90),
    axis.text = element_text(size = base - 2.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
    legend.position = "top", legend.justification = "left",
    legend.title = element_text(size = base - 2.2, colour = MUTED),
    legend.text = element_text(size = base - 2.4, colour = INK),
    legend.key.height = unit(9, "pt"), legend.key.width = unit(15, "pt"),
    strip.text = element_text(size = base - 2, colour = INK, face = "bold",
                              margin = margin(3, 3, 3, 3)),
    plot.tag = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin = margin(6, 8, 6, 6))
}

# within-facet reordering (verbatim from make_publication_figures.R)
.reorder_within <- function(x, by, within) {
  factor(paste(x, within, sep = "___"),
         levels = unique(paste(x, within, sep = "___"))[order(within, by)])
}
.scale_y_reordered <- function() scale_y_discrete(labels = function(z) sub("___.*$", "", z))

.casc_blank <- function(msg) {
  ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 3.4,
                      colour = MUTED, lineheight = 1.05) +
    theme_void(base_family = .CASC_FAMILY)
}

# nice per-fold label: "Round n — dd Mon" from the cutoff date string
.casc_fold_labels <- function(cutoffs) {
  d <- as.Date(cutoffs)
  ord <- order(d)
  # Label by the FORECAST ORIGIN (cutoff + 6) -- 33b issues each origin at cutoff_date + 6.
  # English month abbreviations from month.abb, not format("%b"), which follows the machine's
  # LC_TIME locale (a French one prints "juin").
  o <- lfo_origin(d[ord])
  lab <- sprintf("Round %d — %02d %s", seq_along(d[ord]), as.integer(format(o, "%d")),
                 month.abb[as.integer(format(o, "%m"))])
  setNames(lab, as.character(cutoffs[ord]))
}

# recall / random capture curve: the SAME ESTIMAND as the published panel's det_curve()
# (make_manuscript_figure2_cascade.R) and compute_detection_curve() (20_forecast_detail.R) —
# rank WITHIN origin, recall = each origin's share of its own invasions caught, AVERAGED over
# the origins that had any; random = min(k / mean-at-risk, 1). This used to plot the POOLED
# recall (all hits / all invasions) under the same label, so two "Figure2_cascade" files showed
# different numbers. The pooled value is kept beside it (recall_pooled) and is not drawn.
.casc_capture_curve <- function(detail, ks = 1:25) {
  if (is.null(detail) || !nrow(detail)) return(NULL)
  pf <- detail %>% dplyr::group_by(cutoff) %>%
    dplyr::mutate(rk = rank(-p_invasion, ties.method = "max")) %>% dplyr::ungroup()
  tot <- sum(pf$y, na.rm = TRUE); if (tot == 0) return(NULL)
  natr <- pf %>% dplyr::count(cutoff) %>% dplyr::pull(n) %>% mean()
  idx_by_origin <- split(seq_len(nrow(pf)), pf$cutoff)
  purrr::map_dfr(ks, function(k) {
    inK <- pf$rk <= k
    per <- vapply(idx_by_origin, function(ix) {
      np <- sum(pf$y[ix], na.rm = TRUE)
      if (np == 0) NA_real_ else sum(pf$y[ix][inK[ix]], na.rm = TRUE) / np
    }, numeric(1))
    tibble(k = k, recall = mean(per[is.finite(per)]),
           recall_pooled = sum(pf$y[inK], na.rm = TRUE) / tot,
           recall_random = pmin(k / natr, 1))
  })
}

# ---- Panel A: discrimination (AUC-PR skill per backtest origin) -------------
.casc_panel_discrimination <- function(bt, fold_lab) {
  d <- bt %>% dplyr::mutate(fold = factor(fold_lab[as.character(cutoff)],
                                          levels = rev(unname(fold_lab))))
  ggplot(d, aes(auc_pr_skill, fold)) +
    geom_vline(xintercept = 1, linetype = "22", colour = FAINT) +
    geom_segment(aes(x = 1, xend = auc_pr_skill, yend = fold),
                 colour = GRID, linewidth = 1) +
    geom_point(colour = CASC_ACCENT, size = 2.6) +
    geom_text(aes(label = sprintf("%.0f×", auc_pr_skill)), hjust = -0.35,
              size = 2.9, colour = INK) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.16))) +
    labs(title = "A · Discrimination",
         subtitle = "AUC-PR ÷ base rate  (1× = base rate; a random ordering scores above 1×)",
         x = "AUC-PR skill (× base rate)", y = NULL) +
    .theme_casc() + theme(panel.grid.major.y = element_blank())
}

# ---- Panel B: real-time prioritisation (capture curve) ----------------------
.casc_panel_prioritisation <- function(detail) {
  cc <- .casc_capture_curve(detail)
  if (is.null(cc)) return(.casc_blank("Prioritisation curve\nunavailable"))
  ggplot(cc, aes(k, recall)) +
    geom_line(aes(y = recall_random), linetype = "22", colour = FAINT, linewidth = 0.7) +
    geom_line(colour = CASC_ACCENT, linewidth = 1) +
    geom_point(colour = CASC_ACCENT, size = 1.1) +
    annotate("text", x = max(cc$k) * 0.62, y = max(cc$recall_random) * 0.7 + 0.05,
             label = "random targeting", colour = MUTED, size = 3, angle = 4) +
    scale_y_continuous(labels = percent_format(1), limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    labs(title = "B · Real-time prioritisation",
         subtitle = "vs a same-size random watch-list",
         x = "Zones actively monitored each week (K)",
         y = "Share of true invasions caught") +
    .theme_casc()
}

# ---- Panel C: ranking accuracy (mean rank of invaded zones) -----------------
.casc_panel_meanrank <- function(bt, fold_lab) {
  d <- bt %>% dplyr::filter(is.finite(mean_rank_of_truth)) %>%
    dplyr::mutate(fold = factor(fold_lab[as.character(cutoff)],
                                levels = rev(unname(fold_lab))))
  if (!nrow(d)) return(.casc_blank("Mean-rank panel\nunavailable"))
  natr <- round(mean(bt$n_atrisk))
  ggplot(d, aes(mean_rank_of_truth, fold)) +
    geom_segment(aes(x = 0, xend = mean_rank_of_truth, yend = fold),
                 colour = GRID, linewidth = 1) +
    geom_point(colour = CASC_ACCENT, size = 2.6) +
    geom_text(aes(label = sprintf("%.0f", mean_rank_of_truth)), hjust = -0.4,
              size = 2.9, colour = INK) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.14))) +
    labs(title = "C · Ranking accuracy",
         subtitle = sprintf("among ~%d at-risk zones per origin", natr),
         x = "Mean rank of invaded zones (lower is better)", y = NULL) +
    .theme_casc() + theme(panel.grid.major.y = element_blank())
}

# ---- Panel D: per-round top-K forecasts vs realised outcome -----------------
.casc_panel_topk <- function(detail, fold_lab, top_k = 12L) {
  if (is.null(detail) || !nrow(detail)) return(.casc_blank("Per-round forecasts\nunavailable"))
  # Read the window length from the data. It was hard-coded as "+6 weeks" in the panel
  # title while the backtest horizon is configurable (CASCADE_CALIB_K), so a different K
  # would have published a figure whose caption named the wrong window.
  .K <- if ("K" %in% names(detail)) as.integer(stats::median(detail$K, na.rm = TRUE)) else NA_integer_
  d <- detail %>%
    dplyr::mutate(round_lab = factor(fold_lab[as.character(cutoff)], levels = unname(fold_lab))) %>%
    dplyr::group_by(cutoff) %>%
    dplyr::slice_max(p_invasion, n = top_k, with_ties = FALSE) %>% dplyr::ungroup() %>%
    dplyr::mutate(outcome = factor(ifelse(y == 1L, "Invaded (first case)", "Not invaded"),
                                   levels = names(OUT_COL)),
                  zone_w = .reorder_within(health_zone, p_invasion, round_lab))
  ggplot(d, aes(p_invasion, zone_w, fill = outcome)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    facet_wrap(~ round_lab, scales = "free_y", nrow = 1) +
    .scale_y_reordered() +
    scale_fill_manual(values = OUT_COL, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.08)), breaks = pretty_breaks(3)) +
    labs(title = sprintf("D · Per-round forecasts vs realised outcome (top %d%s)", top_k,
                         if (is.na(.K)) "" else sprintf(", +%d weeks", .K)),
         subtitle = "Bar = predicted P(first case) at each backtest origin; highlighted where the zone was actually invaded",
         x = "Predicted invasion probability, P(first case)", y = NULL) +
    .theme_casc() +
    theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(9, "pt"),
          strip.clip = "off", axis.text.y = element_text(size = 8, colour = INK))
}

# ---- assemble + save --------------------------------------------------------
cascade_fig_evaluation <- function(ctx, file = "Figure2_cascade") {
  bt <- ctx$backtest
  if (is.null(bt) || !nrow(bt)) {
    message("[cascade] Figure2_cascade skipped: no backtest available."); return(invisible(NULL))
  }
  detail <- attr(bt, "detail") %||% ctx$backtest_detail
  if (is.null(detail)) {   # fall back to the on-disk detail CSV
    dcsv <- file.path(OUT_CASCADE, "diagnostics", "cascade_backtest_detail.csv")
    if (file.exists(dcsv)) detail <- readr::read_csv(dcsv, show_col_types = FALSE)
  }
  fold_lab <- .casc_fold_labels(unique(bt$cutoff))

  pA <- .casc_panel_discrimination(bt, fold_lab)
  pB <- .casc_panel_prioritisation(detail)
  pC <- .casc_panel_meanrank(bt, fold_lab)
  pD <- .casc_panel_topk(detail, fold_lab)

  # The window length and the number of origins come from the data, not from prose: the
  # backtest horizon is configurable (CASCADE_CALIB_K) and the number of usable origins
  # grows with the record, so both were wrong as soon as either changed.
  # NA (not a typed 6) when the backtest predates the K column: the caption then names no window.
  .Kcap <- if ("K" %in% names(bt)) as.integer(stats::median(bt$K, na.rm = TRUE)) else NA_integer_
  # Whether the level was fitted on these origins: recorded by run_cascade.R; when absent (a
  # standalone regeneration from an older RDS) the caption says it cannot tell rather than guess.
  .ins <- ctx$backtest_level_in_sample
  cap <- paste0(
    sprintf("Leave-future-out evaluation of the 3-month invasion cascade over its backtest window (%s%d backtest origin(s) on this record). ",
            if (is.na(.Kcap)) "" else sprintf("+%d weeks; ", .Kcap), nrow(bt)),
    "A: discrimination, AUC-PR relative to the base rate (1× = the base rate; a random ordering scores somewhat above it). ",
    "B: share of true invasions caught when the top-K highest-risk zones are watched each week (per origin, averaged over origins), vs a same-size random list. ",
    "C: mean rank of the truly-invaded zones among all at-risk zones. D: the top predicted zones at each origin, coloured by ",
    if (is.na(.Kcap)) "whether they were actually invaded within the backtest window. "
    else sprintf("whether they were actually invaded within the next %d weeks. ", .Kcap),
    if (isTRUE(.ins))
      "The probability LEVEL (bar lengths in D) is in-sample: the cascade's scale factor was fitted to the pooled count on these origins. "
    else if (isFALSE(.ins)) "The scale factor was not fitted on these origins. "
    else "Whether the scale factor was fitted on these origins is not recorded for this figure. ",
    "Discrimination & ranking are the trustworthy products; absolute probabilities are conditional on the scenario.")

  fig <- (pA | pB | pC) / pD +
    plot_layout(heights = c(1, 1.15)) +
    plot_annotation(
      title = "Evaluation of the 3-month spatial-invasion cascade (leave-future-out backtest)",
      caption = paste(strwrap(cap, width = 172), collapse = "\n"),
      theme = theme(
        plot.title = element_text(size = 13.5, face = "bold", colour = INK, margin = margin(b = 2)),
        plot.caption = element_text(size = 7, colour = MUTED, hjust = 0, margin = margin(t = 8)),
        plot.margin = margin(10, 12, 8, 10)))

  # Unicode-safe vector PDF: figure_pdf_device() (00_config.R) is a device that opened cleanly
  # in this session (cairo_pdf, else Quartz); the base pdf device turns dashes into dots.
  save_pdf <- function(path, w, h) {
    ggsave(path, fig, width = w, height = h, device = figure_pdf_device(), bg = "white")
    invisible(path)
  }

  W <- 13.2; H <- 9.6
  # 2026-09-17 streamlining: the PUBLISHED Figure2_cascade is the manuscript-styled
  # one from make_manuscript_figure2_cascade.R (key_outputs/manuscript_figures/).
  # This module keeps only its internal cascade/figures/ copy — the key_outputs/figures/
  # twin was a second, differently-styled figure under the same name.
  casc_dir <- file.path(OUT_CASCADE, "figures")
  for (dir in c(casc_dir)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    save_pdf(file.path(dir, paste0(file, ".pdf")), W, H)
    ggsave(file.path(dir, paste0(file, ".png")), fig, width = W, height = H, dpi = 400, bg = "white")
  }

  # figure data (per-origin backtest summary) + inputs RDS for fast standalone regen
  # CASCADE_OUT_DIR, not OUT_DIR: the bundle follows the run's output root (30_projection_config.R).
  key_dir <- file.path(get0("CASCADE_OUT_DIR", ifnotfound = OUT_DIR), "key_outputs")
  if (!dir.exists(key_dir)) dir.create(key_dir, recursive = TRUE, showWarnings = FALSE)
  # RENAMED. This frame is the PER-ORIGIN BACKTEST SUMMARY (cutoff, K, n_atrisk, n_events,
  # base_rate, auc_pr, auc_pr_skill, mean_rank_of_truth, prec_at<K> for TOPK_VALUES,
  # pred_new, obs_new,
  # count_ratio). None of those columns appears in the published Figure 2 cascade panel, whose
  # data are the pooled prioritisation curve and the per-origin top-12 outcomes — so a reader
  # matching "Figure2_cascade_data.csv" to the figure found neither panel. The published figure
  # is the manuscript-styled one, and make_manuscript_figure2_cascade.R now writes its own data
  # under that name.
  # Stamped with the run that produced it, so a stale copy cannot pass as current.
  # THE BACKTEST'S OWN ANALYSIS DATE. From run_cascade.R (ctx without one) that is this
  # session's ANALYSIS_DATE; a standalone regeneration carries the date saved in the inputs RDS
  # (NA for an RDS written before the date was saved), never the date of the session that
  # happens to redraw it.
  .adate <- if ("analysis_date" %in% names(ctx)) ctx[["analysis_date"]]
            else get0("ANALYSIS_DATE", ifnotfound = NA)
  .adate <- if (length(.adate) == 1L && !is.na(.adate)) format(as.Date(.adate)) else NA_character_
  .bt_csv <- bt
  if (!is.null(ctx$cascade_run_id)) .bt_csv$cascade_run_id <- ctx$cascade_run_id
  .bt_csv$analysis_date <- .adate
  if (!is.null(ctx$backtest_level_in_sample)) .bt_csv$level_in_sample <- isTRUE(ctx$backtest_level_in_sample)
  readr::write_csv(.bt_csv, file.path(key_dir, "cascade_backtest_per_origin.csv"))
  saveRDS(list(backtest = bt, cascade_run_id = ctx$cascade_run_id,
               backtest_level_in_sample = ctx$backtest_level_in_sample,
               analysis_date = .adate),
          file.path(OUT_CASCADE, "diagnostics", "cascade_eval_figure_inputs.rds"))
  message(sprintf("[cascade] wrote %s.{pdf,png} + per-origin data CSV + inputs RDS", file))
  invisible(fig)
}

# ---- standalone regeneration (Rscript 39_cascade_eval_figure.R) -------------
# Rebuilds the figure from the saved inputs RDS (seconds, no model re-run). Sourcing
# this file from run_cascade.R does NOT trigger it (its --file= is run_cascade.R).
if (any(grepl("39_cascade_eval_figure.R", commandArgs(trailingOnly = FALSE), fixed = TRUE))) {
  suppressPackageStartupMessages({ library(here) })
  .ST <- Sys.getenv("CASCADE_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
  source(file.path(.ST, "00_config.R"))
  source(file.path(.ST, "30_projection_config.R"))
  .rds <- file.path(OUT_CASCADE, "diagnostics", "cascade_eval_figure_inputs.rds")
  if (!file.exists(.rds)) stop("eval-figure inputs not found: ", .rds,
                               " — run run_cascade.R once to create them.")
  .inp <- readRDS(.rds)
  # An RDS written before the analysis date was saved: unknown, recorded as NA rather than
  # stamped with this session's date.
  if (!"analysis_date" %in% names(.inp)) .inp["analysis_date"] <- list(NA_character_)
  cascade_fig_evaluation(.inp)
  message("[done] Figure2_cascade regenerated from ", basename(.rds))
}
