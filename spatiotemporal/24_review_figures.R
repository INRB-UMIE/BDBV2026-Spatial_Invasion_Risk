# =============================================================================
# 24_review_figures.R - Publication figures for the review-response analyses
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Renders the NEW review analyses (2026-08-06 review) as figures in the same
# modern, colourblind-safe, Nature-ready house style as the key_outputs
# manuscript figures (OKABE palette, theme_pub, PDF + 600-dpi PNG):
#   §3.4  prospective pre-invasion ranks <- prospective_invasion_per_zone.csv
#         (zones first invaded after every scored back-test outcome window)
#   §2.2  held-out optimism gap        <- invasion_heldout_optimism.csv
#
# CALIBRATION IS NOT HERE. The reliability diagram this file used to draw (FigR1_reliability)
# is now panel C of FigureS7_calibration (make_si_model_figures.R), which leads with the
# per-round expected-vs-observed view the record can actually support and keeps the pooled
# diagram as a third panel. That figure RECOMPUTES the curve from the back-test frame with the
# same estimator (reliability_breaks()/reliability_curve(), 23_prospective_eval.R) and the same
# binning constants (RELIABILITY_*, 00_config.R); it does not read the published tables.
# forecast_reliability_h*.csv and their per-origin twins are still written by run_all.R, as the
# tabulated version of the same curve.
#
# THE HOUSE STYLE COMES FROM figure_style.R. FigR3 and FigR4 are retained deliverables, and
# this file used to define a fourth private design system at base = 11.5 while the house
# theme_pub() is 8.6 — under a header claiming "the same house style as the key_outputs
# manuscript figures". They now use the same theme as every other retained figure.
# make_review_figures() is guarded and a no-op when an input CSV is absent, so it never
# breaks the pipeline.
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(readr); library(dplyr); library(tidyr) })
source(file.path(here::here(), "spatiotemporal", "figure_style.R"))

# Aliases onto the shared tokens, so the panel code below reads unchanged.
.RV_INK <- INK; .RV_MUTED <- MUTED; .RV_FAINT <- FAINT; .RV_GRID <- GRID
.RV_OKABE <- OKABE
.RV_BLUE <- OKABE[1]; .RV_ORANGE <- OKABE[2]; .RV_GREEN <- OKABE[3]; .RV_FAM <- base_family

.rv_theme <- function(base = 8.6) theme_pub(base)
.rv_save <- function(p, path, w, h) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !.fk(path)) return(invisible(p))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(path, ".pdf"), p, width = w, height = h, device = figure_pdf_device(),
                  bg = "white")
  # The previous PNG goes first, and a failed write is reported: swallowing it left last run's
  # PNG beside this run's PDF.
  invisible(suppressWarnings(file.remove(paste0(path, ".png"))))
  tryCatch(ggplot2::ggsave(paste0(path, ".png"), p, width = w, height = h, dpi = 600, bg = "white"),
           error = function(e) warning("[review-fig] ", basename(path), ".png not written: ",
                                       conditionMessage(e), call. = FALSE))
  message(sprintf("  [review-fig] saved %-34s %.1f x %.1f in", basename(path), w, h)); invisible(p)
}

# A missing or empty input must NOT vanish silently. FigR3_prospective_ranks and
# FigR4_heldout_optimism are published deliverables; when
# run_all.R fails to write one of their input CSVs the figure simply never appeared, with no
# warning anywhere and nothing distinguishing "not produced" from "not requested". Warn, but
# only for figures the retained-figure gate would actually write, so a suppressed figure
# stays quiet.
.rv_missing <- function(csv, out, what) {
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !isTRUE(.fk(out))) return(invisible(NULL))
  # A figure this run cannot draw must not survive from an earlier run.
  suppressWarnings(file.remove(paste0(out, c(".pdf", ".png"))[file.exists(paste0(out, c(".pdf", ".png")))]))
  warning(sprintf("[review-fig] %s NOT produced: %s (%s)",
                  basename(out), what, csv), call. = FALSE)
  invisible(NULL)
}

# ---- §3.4 Prospective pre-invasion ranks ------------------------------------
#' For zones first invaded AFTER every scored back-test outcome window -- first case in a week
#' starting after the last week of any (round, horizon) outcome window in the back-test frame,
#' and no scored back-test invasion event (run_all.R; prospective_invasion_check(),
#' 23_prospective_eval.R) -- the predicted-risk rank each held in the forecast issued at the
#' last back-test origin, before its first case (lollipop). Zones invaded inside a scored
#' window are back-test events the selection already used, so they are not shown. Coloured by
#' whether they fell inside the top-K watch-list; the K cut-off is drawn as a reference line.
plot_prospective_ranks_review <- function(csv, out, top_k = 15L) {
  if (!file.exists(csv)) return(.rv_missing(csv, out, "input CSV is absent"))
  d <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE))
  if (!nrow(d)) return(.rv_missing(csv, out, "input CSV is empty"))
  d <- d[order(d$pre_invasion_rank), , drop = FALSE]
  d$health_zone <- factor(d$health_zone, levels = rev(d$health_zone))
  d$flag <- ifelse(d$in_topk, sprintf("In top-%d watch-list", top_k), "Below watch-list")
  cols <- setNames(c(.RV_GREEN, .RV_ORANGE),
                   c(sprintf("In top-%d watch-list", top_k), "Below watch-list"))
  # n_hit and med used to be computed here and never used — a hit count and a median rank
  # formatted into nothing. Removed rather than rendered: this panel's message is the per-zone
  # ranks against the watch-list line, and the pooled top-K hit rate is already published in
  # invasion_evaluation.csv (hit_at_<K> for each K in TOPK_VALUES).
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

# ---- §2.2 Held-out optimism, DECOMPOSED --------------------------------------
#' Inner vs held-out AUC-PR skill for EVERY eligible candidate, as a slope chart.
#'
#' WHAT THIS FIGURE REPLACED, AND WHY. It used to be two bars: the selected model's inner-CV
#' skill and its held-out skill, with the difference captioned as the optimism gap. That
#' difference conflates two unrelated effects. The held-out origins are a different window of
#' the series, with their own base rate and difficulty, so EVERY candidate's skill moves there
#' (a window effect). Only the amount by which the SELECTED model drops MORE than the field is
#' optimism. The decomposition is published in invasion_heldout_decomposition.csv.
#'
#' Drawing the whole field makes the decomposition visible rather than asserted: a fan of
#' near-parallel lines sloping down is a window effect; a highlighted line that is visibly
#' steeper than the fan is selection. The reader can see which one this is.
#'
#' TWO MODELS ARE HIGHLIGHTED, deliberately. The inner-fold selection need not return the
#' model the pipeline ships (is_selected vs is_featured in the CSV). Showing only the inner
#' pick would report optimism for a model nothing ships; showing only the deployed one would
#' hide it whenever the honest test selects something else. When the two coincide, the line
#' is drawn once, as the selected model.
#'
#' Skill is plotted on the COMMON (inner) base rate: both points of every line are AUC-PR
#' divided by the selection folds' base rate, so a slope is a change in AUC-PR, not a change in
#' base rate. On each set's own base rate the field would move out of sample purely because the
#' held-out window's base rate differs from the selection folds'.
#' @param csv invasion_heldout_field.csv (one row per candidate x horizon).
plot_heldout_optimism_review <- function(csv, out) {
  if (!file.exists(csv)) return(.rv_missing(csv, out, "input CSV is absent"))
  d <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE))
  if (!nrow(d)) return(.rv_missing(csv, out, "input CSV is empty"))
  need <- c("method", "horizon", "auc_pr_inner", "auc_pr_outer", "base_rate_inner",
            "is_selected", "is_featured")
  if (!all(need %in% names(d)))
    return(.rv_missing(csv, out, paste("missing column(s):",
                                       paste(setdiff(need, names(d)), collapse = ", "))))
  # COMMON BASE RATE: both ends over the inner base rate, so the slope is a difference of
  # numerators and not a base-rate shift. Computed per horizon.
  d <- d %>%
    dplyr::group_by(horizon) %>%
    dplyr::mutate(inner = auc_pr_inner / base_rate_inner[1],
                  outer = auc_pr_outer / base_rate_inner[1]) %>%
    dplyr::ungroup()
  lab_h <- function(h) paste0(h, "-week horizon")
  .ax <- if ("axis" %in% names(d)) unique(stats::na.omit(d$axis)) else character(0)
  if (length(.ax) > 1L)
    return(.rv_missing(csv, out, paste("the table mixes selection axes:", paste(.ax, collapse = "; "))))
  .rv_axis_lab <- if (length(.ax)) paste0(toupper(substr(.ax, 1, 1)), substring(.ax, 2))
                  else "AUC-PR skill on the selection rounds"
  long <- d %>%
    tidyr::pivot_longer(c(inner, outer), names_to = "set", values_to = "skill") %>%
    dplyr::mutate(set = factor(set, levels = c("inner", "outer"),
                               labels = c("Selection folds", "Held-out origins")),
                  role = dplyr::case_when(is_selected ~ "Selected on the inner folds",
                                          is_featured ~ "Deployed model",
                                          TRUE        ~ "Other candidates"),
                  panel = lab_h(horizon))
  # A model can be both selected and featured; case_when above already gives "selected"
  # priority, so the legend never shows one line under two roles.
  long$role <- factor(long$role, levels = c("Other candidates", "Deployed model",
                                            "Selected on the inner folds"))
  cols <- c("Other candidates" = .RV_FAINT, "Deployed model" = .RV_GREEN,
            "Selected on the inner folds" = .RV_ORANGE)
  # Field median drawn as a thick reference slope: this IS the window effect. Over the FIELD
  # rows only (in_field), as .heldout_field_summary() computes the published decomposition: a
  # featured form outside the field must not move the median it is compared against.
  med <- (if ("in_field" %in% names(long)) dplyr::filter(long, in_field %in% TRUE) else long) %>%
    dplyr::group_by(panel, set) %>%
    dplyr::summarise(skill = stats::median(skill, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(grp = "median")
  # NO TITLE, NO SUBTITLE, NO STRIP HEADER. .rv_theme() already blanks plot.title and
  # plot.subtitle; the facet strip was the remaining title-like band, so the horizon is set
  # as a quiet in-panel label instead. Nothing above the plotting area but the key.
  hz <- long %>% dplyr::distinct(panel)
  p <- ggplot2::ggplot(long, ggplot2::aes(set, skill, group = method)) +
    ggplot2::geom_line(ggplot2::aes(colour = role, linewidth = role, alpha = role),
                       lineend = "round") +
    ggplot2::geom_point(ggplot2::aes(colour = role, size = role)) +
    # The field median IS the window effect, so it carries a key: unlabelled it was a dashed line
    # a reader could not name.
    ggplot2::geom_line(data = med, ggplot2::aes(set, skill, group = grp,
                                                linetype = "Median of the candidates (window effect)"),
                       inherit.aes = FALSE, colour = .RV_INK, linewidth = 0.8, lineend = "round") +
    ggplot2::geom_point(data = med, ggplot2::aes(set, skill), inherit.aes = FALSE,
                        colour = .RV_INK, size = 1.5) +
    ggplot2::geom_text(data = hz, ggplot2::aes(x = -Inf, y = Inf, label = panel),
                       inherit.aes = FALSE, hjust = -0.02, vjust = 1.1,
                       size = 3.0, colour = .RV_MUTED) +
    ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = cols, name = NULL, breaks = names(cols)[c(3, 2, 1)],
                                 guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_linetype_manual(values = c("Median of the candidates (window effect)" = "12"),
                                   name = NULL,
                                   guide = ggplot2::guide_legend(order = 2,
                                     override.aes = list(colour = .RV_INK, linewidth = 0.8))) +
    ggplot2::scale_linewidth_manual(values = c(0.32, 1.05, 1.05), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(0.5, 1, 1), guide = "none") +
    ggplot2::scale_size_manual(values = c(0.7, 1.9, 1.9), guide = "none") +
    # Breathing room at both ends so no line runs into the panel edge, and headroom at the
    # top for the in-panel horizon label.
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.13, 0.13))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.07, 0.26))) +
    # NO CAPTION ON THE FIGURE. The decomposition (field median drop, each model's excess,
    # the rank correlation) is published in invasion_heldout_decomposition.csv and reaches the
    # manuscript through extract_numbers.py's OPTIMISM_WINDOW / OPTIMISM_SELECTION / OPTIMISM_RHO
    # keys, so the figure carries the pattern and the manuscript caption carries the numbers.
    # Keeping them in two places was how the old two-bar version came to assert an optimism
    # figure the picture did not support.
    # The y axis names the SELECTION'S axis, read from the table (run_all.R writes `axis` from
    # INVASION_SELECTION_BASIS): this check is the one place performance is shown on the selection's
    # rounds rather than all rounds, because the optimism of a choice is measured on what it was
    # chosen on.
    ggplot2::labs(x = NULL, y = paste0(.rv_axis_lab, " (x selection-fold base rate)")) +
    .rv_theme() +
    ggplot2::theme(
      # The strip is the title band; blanking it is the point of this revision.
      strip.text = ggplot2::element_blank(),
      # A slope chart reads along the lines, so vertical rules are noise. Keep a faint
      # horizontal rhythm only, and drop the axis ticks with it.
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = .RV_GRID, linewidth = 0.25),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 8.2, colour = .RV_INK,
                                          margin = ggplot2::margin(t = 1)),
      axis.text.y = ggplot2::element_text(size = 7.6, colour = .RV_FAINT),
      axis.title.y = ggplot2::element_text(size = 8.4, colour = .RV_MUTED, angle = 90),
      panel.spacing = ggplot2::unit(18, "pt"),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 8, colour = .RV_INK),
      legend.margin = ggplot2::margin(b = 2),
      plot.margin = ggplot2::margin(4, 10, 6, 6)) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE,
                                                   override.aes = list(linewidth = 1.1, size = 1.9)))
  .rv_save(p, out, 7.0, 3.5)
}

make_review_figures <- function(diag_dir, fig_dir, top_k = 15L) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  message("[review-fig] rendering review-response figures (house style) ...")
  tryCatch(plot_prospective_ranks_review(file.path(diag_dir, "prospective_invasion_per_zone.csv"),
             file.path(fig_dir, "FigR3_prospective_ranks"), top_k),
           error = function(e) warning("[review-fig] prospective: ", conditionMessage(e)))
  tryCatch(plot_heldout_optimism_review(file.path(diag_dir, "invasion_heldout_field.csv"),
             file.path(fig_dir, "FigR4_heldout_optimism")),
           error = function(e) warning("[review-fig] held-out: ", conditionMessage(e)))
  invisible(TRUE)
}
