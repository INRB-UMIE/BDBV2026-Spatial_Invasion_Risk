# =============================================================================
# 39_cascade_eval_figure.R — CASCADE evaluation figure, in the Figure 2 style
# BDBV 2026 DRC · 3-month spatial invasion cascade
#
# A faithful cascade analogue of make_publication_figures.R's Figure 2 (the
# short-horizon backtest evaluation), for the LONGER-TERM invasion cascade. The
# cascade is a single model, so where Figure 2 contrasts many models / two short
# horizons, this figure contrasts the leave-future-out backtest ORIGINS (folds)
# at the deepest horizon the ~15-week record supports (+6 weeks). Same design
# system, same panels, honest about the data ceiling:
#
#   A  Discrimination          AUC-PR skill (× base rate) per backtest origin.
#   B  Real-time prioritisation share of true invasions caught vs a random
#                              watch-list of the same size (capture curve).
#   C  Ranking accuracy        mean rank of the truly-invaded zones (lower better).
#   D  Per-round forecasts      top-K predicted zones per origin, coloured by the
#                              realised outcome (invaded within +6 wk, or not).
#
# All four panels come from ONE source — the cascade leave-future-out backtest
# (cascade_backtest → summary + per-zone `detail`) — so they are mutually
# consistent and fully out-of-sample. Entry: cascade_fig_evaluation(ctx), called
# from run_cascade.R (needs ctx$backtest carrying its `detail` attribute).
# Standalone regen: Rscript spatiotemporal/39_cascade_eval_figure.R (from the RDS).
# Saves Figure2_cascade.{pdf,png} to key_outputs/figures (+ mirror in cascade/figures).
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
# suppresses titles because its captions live in FIGURE_CAPTIONS.md).
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
  lab <- sprintf("Round %d — %s", seq_along(d[ord]), format(d[ord], "%d %b"))
  setNames(lab, as.character(cutoffs[ord]))
}

# recall / random capture curve, pooled over folds (mirrors .det_curve exactly:
# rank WITHIN fold, total events over all folds, random = min(k / mean-at-risk, 1)).
.casc_capture_curve <- function(detail, ks = 1:25) {
  if (is.null(detail) || !nrow(detail)) return(NULL)
  pf <- detail %>% dplyr::group_by(cutoff) %>%
    dplyr::mutate(rk = rank(-p_invasion, ties.method = "max")) %>% dplyr::ungroup()
  tot <- sum(pf$y, na.rm = TRUE); if (tot == 0) return(NULL)
  natr <- pf %>% dplyr::count(cutoff) %>% dplyr::pull(n) %>% mean()
  purrr::map_dfr(ks, function(k) {
    tp <- sum(pf$y[pf$rk <= k], na.rm = TRUE)
    tibble(k = k, recall = tp / tot, recall_random = pmin(k / natr, 1))
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
         subtitle = "AUC-PR ÷ base rate  (1× = random)",
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
    labs(title = sprintf("D · Per-round forecasts vs realised outcome (top %d, +6 weeks)", top_k),
         subtitle = "Bar = predicted P(first case) at each backtest origin; highlighted where the zone was actually invaded",
         x = "Predicted invasion probability, P(first case)", y = NULL) +
    .theme_casc() +
    theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(9, "pt"),
          strip.clip = "off", axis.text.y = element_text(size = 8, colour = INK))
}

# ---- assemble + save --------------------------------------------------------
cascade_fig_evaluation <- function(ctx, file = "Figure2_cascade") {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
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

  cap <- paste0(
    "Leave-future-out evaluation of the 3-month invasion cascade at its deepest supported horizon (+6 weeks; a ~15-week ",
    "record affords ~3 backtest origins). A: out-of-sample discrimination, AUC-PR relative to the base rate (1× = random). ",
    "B: share of true invasions caught when the top-K highest-risk zones are watched each week, vs a same-size random list. ",
    "C: mean rank of the truly-invaded zones among all at-risk zones. D: the top predicted zones at each origin, coloured by ",
    "whether they were actually invaded within the next 6 weeks. Discrimination & ranking are the trustworthy products; ",
    "absolute probabilities are upper bounds.")

  fig <- (pA | pB | pC) / pD +
    plot_layout(heights = c(1, 1.15)) +
    plot_annotation(
      title = "Evaluation of the 3-month spatial-invasion cascade (leave-future-out backtest)",
      caption = paste(strwrap(cap, width = 172), collapse = "\n"),
      theme = theme(
        plot.title = element_text(size = 13.5, face = "bold", colour = INK, margin = margin(b = 2)),
        plot.caption = element_text(size = 7, colour = MUTED, hjust = 0, margin = margin(t = 8)),
        plot.margin = margin(10, 12, 8, 10)))

  # unicode-safe vector PDF via quartz (base "pdf" drops ×/÷; cairo_pdf is dead on
  # this box — capabilities("cairo")==TRUE but the DLL fails to load, silent no-write).
  save_pdf <- function(path, w, h) {
    if (isTRUE(capabilities("aqua"))) {
      ok <- tryCatch({ grDevices::quartz(type = "pdf", file = path, width = w, height = h,
                                         bg = "white"); print(fig); grDevices::dev.off(); TRUE },
                     error = function(e) FALSE)
      if (ok && file.exists(path)) return(invisible(path))
    }
    ggsave(path, fig, width = w, height = h, device = "pdf", bg = "white")
    invisible(path)
  }

  W <- 13.2; H <- 9.6
  fig_dir <- file.path(OUT_DIR, "key_outputs", "figures")
  casc_dir <- file.path(OUT_CASCADE, "figures")
  for (dir in c(fig_dir, casc_dir)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    save_pdf(file.path(dir, paste0(file, ".pdf")), W, H)
    ggsave(file.path(dir, paste0(file, ".png")), fig, width = W, height = H, dpi = 400, bg = "white")
  }

  # figure data (per-origin backtest summary) + inputs RDS for fast standalone regen
  key_dir <- file.path(OUT_DIR, "key_outputs")
  if (!dir.exists(key_dir)) dir.create(key_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(bt, file.path(key_dir, "Figure2_cascade_data.csv"))
  saveRDS(list(backtest = bt), file.path(OUT_CASCADE, "diagnostics", "cascade_eval_figure_inputs.rds"))
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
  cascade_fig_evaluation(readRDS(.rds))
  message("[done] Figure2_cascade regenerated from ", basename(.rds))
}
