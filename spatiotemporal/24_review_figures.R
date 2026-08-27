# =============================================================================
# 24_review_figures.R - Publication figures for the review-response analyses
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Renders the NEW review analyses (2026-08-06 review) as figures in the same
# modern, colourblind-safe, Nature-ready house style as the key_outputs
# manuscript figures (OKABE palette, theme_pub, PDF + 600-dpi PNG):
#   §3.5  reliability curve            <- forecast_reliability_h1.csv
#   §3.5  count-calibration over time  <- forecast_count_calibration_h1.csv
#   §3.4  prospective pre-invasion ranks <- prospective_invasion_per_zone.csv
#   §2.2  held-out optimism gap        <- invasion_heldout_optimism.csv
#
# Self-contained (defines its own copy of the design system so it does not depend
# on make_manuscript_figures.R internals). Driver make_review_figures() is guarded
# and a no-op when an input CSV is absent, so it never breaks the pipeline.
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(readr); library(dplyr) })

# ---- Design system (matches make_manuscript_figures.R) ----------------------
.RV_INK <- "grey15"; .RV_MUTED <- "grey38"; .RV_FAINT <- "grey72"; .RV_GRID <- "grey92"
.RV_OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
.RV_BLUE <- "#0072B2"; .RV_ORANGE <- "#D55E00"; .RV_GREEN <- "#009E73"; .RV_FAM <- "sans"

.rv_theme <- function(base = 11.5) {
  ggplot2::theme_minimal(base_size = base, base_family = .RV_FAM) %+replace% ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    axis.title    = ggplot2::element_text(size = base - 0.4, colour = .RV_MUTED),
    axis.title.x  = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
    axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(r = 4), angle = 90),
    axis.text     = ggplot2::element_text(size = base - 1.2, colour = .RV_MUTED),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = .RV_GRID, linewidth = 0.3),
    legend.position = "top", legend.justification = "left",
    legend.title  = ggplot2::element_text(size = base - 1.2, colour = .RV_MUTED),
    legend.text   = ggplot2::element_text(size = base - 1.4, colour = .RV_INK),
    legend.key.height = ggplot2::unit(9, "pt"), legend.key.width = ggplot2::unit(15, "pt"),
    plot.margin   = ggplot2::margin(7, 9, 7, 7))
}
.rv_save <- function(p, path, w, h) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(path, ".pdf"), p, width = w, height = h, device = "pdf", bg = "white")
  tryCatch(ggplot2::ggsave(paste0(path, ".png"), p, width = w, height = h, dpi = 600, bg = "white"),
           error = function(e) NULL)
  message(sprintf("  [review-fig] saved %-34s %.1f x %.1f in", basename(path), w, h)); invisible(p)
}
# Wrap a subtitle to `n` characters so it never overruns the (narrow) panel width.
.rv_wrap <- function(s, n) paste(strwrap(s, width = n), collapse = "\n")

# ---- §3.5 Reliability / calibration curve -----------------------------------
#' Binned predicted probability vs observed invasion frequency (Wilson intervals),
#' with the 45-degree perfect-calibration reference. Points sized by bin count.
plot_reliability_review <- function(csv, out, horizon = 1L, cal_in_large = NULL) {
  if (!file.exists(csv)) return(invisible(NULL))
  d <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE))
  d <- d[is.finite(d$mean_pred) & d$n > 0, , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  lim <- max(d$mean_pred, d$obs_hi, na.rm = TRUE) * 1.05
  sub <- .rv_wrap(sprintf("Featured model, %d-week horizon. Points below the diagonal indicate over-prediction%s.",
                 horizon, if (!is.null(cal_in_large)) sprintf(" (%.2f× overall)", cal_in_large) else ""), 40)
  # Perfect-calibration diagonal; Wilson intervals convey per-bin precision (wide = few
  # zone-weeks), so a fixed, legible point size reads more cleanly than sizing by count
  # (bin 1 alone holds ~4,000 zone-weeks and would dwarf the informative mid-range bins).
  p <- ggplot2::ggplot(d, ggplot2::aes(mean_pred, obs_freq)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22", colour = .RV_FAINT, linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = obs_lo, ymax = obs_hi),
                            colour = .RV_MUTED, linewidth = 0.45, na.rm = TRUE) +
    ggplot2::geom_point(colour = .RV_BLUE, size = 2.4, alpha = 0.95) +
    ggplot2::scale_x_continuous(limits = c(0, lim), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_y_continuous(limits = c(0, lim), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Predicted probability of first case",
                  y = "Observed invasion frequency") +
    .rv_theme() + ggplot2::coord_equal()
  .rv_save(p, out, 4.0, 4.2)
}

# ---- §3.5 Aggregate-count calibration over time -----------------------------
#' Predicted expected number of new invasions per forecast origin (sum of at-risk
#' probabilities) vs the observed number of new invasions - over-prediction and its
#' drift are read directly.
plot_count_calibration_review <- function(csv, out, horizon = 1L) {
  if (!file.exists(csv)) return(invisible(NULL))
  d <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE))
  if (!nrow(d)) return(invisible(NULL))
  # Order origins chronologically: numeric fold indices / ISO dates sort correctly numerically;
  # fall back to lexical order only if `time` is non-numeric and non-date.
  .tnum <- suppressWarnings(as.numeric(d$time))
  ord   <- if (!any(is.na(.tnum))) order(.tnum) else order(d$time)
  d$origin <- factor(d$time, levels = d$time[ord])
  dl <- rbind(
    data.frame(origin = d$origin, count = d$expected, kind = "Predicted (expected)"),
    data.frame(origin = d$origin, count = d$observed, kind = "Observed"))
  p <- ggplot2::ggplot() +
    ggplot2::geom_col(data = subset(dl, kind == "Observed"),
                      ggplot2::aes(origin, count, fill = kind), width = 0.62) +
    ggplot2::geom_line(data = subset(dl, kind == "Predicted (expected)"),
                       ggplot2::aes(origin, count, colour = kind, group = 1), linewidth = 0.7) +
    ggplot2::geom_point(data = subset(dl, kind == "Predicted (expected)"),
                        ggplot2::aes(origin, count, colour = kind), size = 1.9) +
    ggplot2::scale_fill_manual(values = c("Observed" = .RV_FAINT), name = NULL) +
    ggplot2::scale_colour_manual(values = c("Predicted (expected)" = .RV_ORANGE), name = NULL) +
    ggplot2::labs(x = "Forecast origin", y = "New invasions per origin") +
    .rv_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  .rv_save(p, out, 4.9, 4.0)
}

# ---- §3.4 Prospective pre-invasion ranks ------------------------------------
#' For zones invaded AFTER the last CV origin, the predicted-risk rank they held
#' BEFORE their first case (lollipop). Coloured by whether they fell inside the
#' top-K watch-list; the K cut-off is drawn as a reference line.
plot_prospective_ranks_review <- function(csv, out, top_k = 15L) {
  if (!file.exists(csv)) return(invisible(NULL))
  d <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE))
  if (!nrow(d)) return(invisible(NULL))
  d <- d[order(d$pre_invasion_rank), , drop = FALSE]
  d$health_zone <- factor(d$health_zone, levels = rev(d$health_zone))
  d$flag <- ifelse(d$in_topk, sprintf("In top-%d watch-list", top_k), "Below watch-list")
  cols <- setNames(c(.RV_GREEN, .RV_ORANGE),
                   c(sprintf("In top-%d watch-list", top_k), "Below watch-list"))
  n_hit <- sum(d$in_topk); med <- stats::median(d$pre_invasion_rank)
  p <- ggplot2::ggplot(d, ggplot2::aes(pre_invasion_rank, health_zone, colour = flag)) +
    ggplot2::geom_vline(xintercept = top_k + 0.5, linetype = "22", colour = .RV_FAINT, linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = pre_invasion_rank,
                                       y = health_zone, yend = health_zone), linewidth = 0.5) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_text(ggplot2::aes(label = pre_invasion_rank), hjust = -0.5, size = 3.3,
                       colour = .RV_INK, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = "Pre-invasion predicted-risk rank", y = NULL) +
    .rv_theme()
  .rv_save(p, out, 4.4, max(2.7, 0.40 * nrow(d) + 1.6))
}

# ---- §2.2 Held-out optimism gap ---------------------------------------------
#' Inner-CV (selection-optimistic) vs held-out (last two origins, never seen by
#' selection) AUC-PR skill for the featured model, with the optimism gap annotated.
plot_heldout_optimism_review <- function(csv, out) {
  if (!file.exists(csv)) return(invisible(NULL))
  d <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE))
  if (!nrow(d)) return(invisible(NULL))
  r <- d[1, ]
  bars <- data.frame(
    kind = factor(c("All-fold CV\n(selection-optimistic)", "Held-out\n(last 2 origins)"),
                  levels = c("All-fold CV\n(selection-optimistic)", "Held-out\n(last 2 origins)")),
    skill = c(r$inner_cv_skill, r$heldout_skill))
  ymax <- max(bars$skill, na.rm = TRUE) * 1.18
  p <- ggplot2::ggplot(bars, ggplot2::aes(kind, skill, fill = kind)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f×", skill)), vjust = -0.5,
                       size = 4.1, colour = .RV_INK) +
    ggplot2::scale_fill_manual(values = setNames(c(.RV_FAINT, .RV_BLUE), levels(bars$kind)), guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, ymax), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = "AUC-PR skill (× chance)") +
    .rv_theme()
  .rv_save(p, out, 3.8, 4.0)
}

#' Driver: render all four review figures from the pipeline's diagnostics CSVs.
#' @param diag_dir directory holding the review CSVs (OUT_DIAGNOSTICS).
#' @param fig_dir  output directory (default the key_outputs manuscript panels dir).
#' @param cal_in_large optional calibration-in-the-large for the reliability subtitle.
make_review_figures <- function(diag_dir, fig_dir, cal_in_large = NULL, top_k = 15L) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  message("[review-fig] rendering review-response figures (house style) ...")
  tryCatch(plot_reliability_review(file.path(diag_dir, "forecast_reliability_h1.csv"),
             file.path(fig_dir, "FigR1_reliability_h1"), 1L, cal_in_large),
           error = function(e) warning("[review-fig] reliability: ", conditionMessage(e)))
  tryCatch(plot_count_calibration_review(file.path(diag_dir, "forecast_count_calibration_h1.csv"),
             file.path(fig_dir, "FigR2_count_calibration_h1"), 1L),
           error = function(e) warning("[review-fig] count-calibration: ", conditionMessage(e)))
  tryCatch(plot_prospective_ranks_review(file.path(diag_dir, "prospective_invasion_per_zone.csv"),
             file.path(fig_dir, "FigR3_prospective_ranks"), top_k),
           error = function(e) warning("[review-fig] prospective: ", conditionMessage(e)))
  tryCatch(plot_heldout_optimism_review(file.path(diag_dir, "invasion_heldout_optimism.csv"),
             file.path(fig_dir, "FigR4_heldout_optimism")),
           error = function(e) warning("[review-fig] held-out: ", conditionMessage(e)))
  invisible(TRUE)
}
