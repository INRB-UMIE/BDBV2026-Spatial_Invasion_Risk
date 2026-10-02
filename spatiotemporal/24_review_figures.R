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
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !.fk(path)) return(invisible(p))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(path, ".pdf"), p, width = w, height = h, device = "pdf", bg = "white")
  tryCatch(ggplot2::ggsave(paste0(path, ".png"), p, width = w, height = h, dpi = 600, bg = "white"),
           error = function(e) NULL)
  message(sprintf("  [review-fig] saved %-34s %.1f x %.1f in", basename(path), w, h)); invisible(p)
}
# Wrap a subtitle to `n` characters so it never overruns the (narrow) panel width.

# A missing or empty input must NOT vanish silently. FigR1_reliability_h1/_h2,
# FigR3_prospective_ranks and FigR4_heldout_optimism are all published deliverables; when
# run_all.R fails to write one of their input CSVs the figure simply never appeared, with no
# warning anywhere and nothing distinguishing "not produced" from "not requested". Warn, but
# only for figures the retained-figure gate would actually write, so a suppressed figure
# stays quiet.
.rv_missing <- function(csv, out, what) {
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !isTRUE(.fk(out))) return(invisible(NULL))
  warning(sprintf("[review-fig] %s NOT produced: %s (%s)",
                  basename(out), what, csv), call. = FALSE)
  invisible(NULL)
}

# ---- §3.5 Reliability / calibration curve -----------------------------------
#' Binned predicted probability vs observed invasion frequency (Wilson intervals),
#' with the 45-degree perfect-calibration reference. Points sized by bin count.
#' Reliability of the PRIMARY (recalibrated) forecast, with the raw curve overlaid.
#'
#' The two scales belong on one panel: the raw curve is the evidence that motivates the
#' correction and the recalibrated curve is what the deployed system issues, so either
#' alone tells half the story. `csv_raw` is the twin written under diagnostics/raw/; when
#' it is absent the panel degrades to the single (primary) curve.
plot_reliability_review <- function(csv, out, horizon = 1L, cal_in_large = NULL,
                                    csv_raw = NULL, cal_in_large_raw = NULL) {
  if (!file.exists(csv)) return(.rv_missing(csv, out, "input CSV is absent"))
  d <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE))
  d <- d[is.finite(d$mean_pred) & d$n > 0, , drop = FALSE]
  if (!nrow(d)) return(.rv_missing(csv, out, "no bin has finite mean_pred and n > 0"))
  draw <- if (!is.null(csv_raw) && file.exists(csv_raw)) {
    x <- suppressMessages(readr::read_csv(csv_raw, show_col_types = FALSE))
    x <- x[is.finite(x$mean_pred) & x$n > 0, , drop = FALSE]
    if (nrow(x)) x else NULL
  } else NULL
  lim <- max(c(d$mean_pred, d$obs_hi, draw$mean_pred, draw$obs_hi), na.rm = TRUE) * 1.05
  # The series are named in the LEGEND, not in a subtitle: .rv_theme() blanks plot.title and
  # plot.subtitle by design (captions live in the manuscript text), so the `sub` string this
  # function used to build was computed and then silently discarded — and an unlabelled
  # second series would be worse than none. Calibration-in-the-large rides in the label, so
  # the panel states how far each scale is from perfect without extra furniture.
  .lab_cal <- function(base, cil) if (is.null(cil) || !is.finite(cil)) base else
    sprintf("%s (%.2f\u00d7)", base, cil)
  LAB_CAL <- .lab_cal("Recalibrated (prequential)", cal_in_large)
  LAB_RAW <- .lab_cal("Raw", cal_in_large_raw)
  # Perfect-calibration diagonal; Wilson intervals convey per-bin precision (wide = few
  # zone-weeks), so a fixed, legible point size reads more cleanly than sizing by count
  # (bin 1 alone holds ~4,000 zone-weeks and would dwarf the informative mid-range bins).
  p <- ggplot2::ggplot(d, ggplot2::aes(mean_pred, obs_freq)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22", colour = .RV_FAINT, linewidth = 0.4) +
    (if (!is.null(draw)) ggplot2::geom_line(data = draw,
        ggplot2::aes(mean_pred, obs_freq, colour = LAB_RAW), linewidth = 0.5)) +
    (if (!is.null(draw)) ggplot2::geom_point(data = draw,
        ggplot2::aes(mean_pred, obs_freq, colour = LAB_RAW), size = 1.7, shape = 1)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = obs_lo, ymax = obs_hi),
                            colour = .RV_MUTED, linewidth = 0.45, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(colour = LAB_CAL), size = 2.4, alpha = 0.95) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(.RV_BLUE, .RV_FAINT),
                                                          c(LAB_CAL, LAB_RAW)),
                                 breaks = c(LAB_CAL, LAB_RAW), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, lim), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_y_continuous(limits = c(0, lim), labels = scales::percent_format(accuracy = 1)) +
    # STATE THE HORIZON. `horizon` was accepted, documented and passed positionally, but the
    # function body never read it — it used to reach the panel through a subtitle that was
    # deleted. With FigR1 now produced at BOTH horizons, the two deliverables were visually
    # identical in what they claimed to show.
    ggplot2::labs(x = sprintf("Predicted probability of first case, %d week%s ahead",
                              as.integer(horizon), if (as.integer(horizon) == 1L) "" else "s"),
                  y = "Observed invasion frequency") +
    .rv_theme() + ggplot2::coord_equal()
  .rv_save(p, out, 4.0, 4.2)
}


# ---- §3.4 Prospective pre-invasion ranks ------------------------------------
#' For zones invaded AFTER the last CV origin, the predicted-risk rank they held
#' BEFORE their first case (lollipop). Coloured by whether they fell inside the
#' top-K watch-list; the K cut-off is drawn as a reference line.
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
  # invasion_evaluation.csv (hit_at_5/10/15).
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
  if (!file.exists(csv)) return(.rv_missing(csv, out, "input CSV is absent"))
  d <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE))
  if (!nrow(d)) return(.rv_missing(csv, out, "input CSV is empty"))
  r <- d[1, ]
  bars <- data.frame(
    # "Inner CV", not "All-fold CV": evaluate_invasion_heldout() computes inner_cv_skill on
  # setdiff(cutoffs, outer_cut) — the first 12 of 14 origins — NOT on all folds. The all-fold
  # h=1 skill for the same model is 47.1x (the number in the leaderboard, in Figure 2 and in the
  # report prose), while this bar reads 49.1x, so the manuscript carried two different
  # "all-fold CV" numbers for one model. The CSV column is literally named inner_cv_skill.
  kind = factor(c("Inner CV\n(selection folds, optimistic)", "Held-out\n(last 2 origins)"),
                  levels = c("Inner CV\n(selection folds, optimistic)", "Held-out\n(last 2 origins)")),
    skill = c(r$inner_cv_skill, r$heldout_skill))
  ymax <- max(bars$skill, na.rm = TRUE) * 1.18
  p <- ggplot2::ggplot(bars, ggplot2::aes(kind, skill, fill = kind)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f×", skill)), vjust = -0.5,
                       size = 4.1, colour = .RV_INK) +
    ggplot2::scale_fill_manual(values = setNames(c(.RV_FAINT, .RV_BLUE), levels(bars$kind)), guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, ymax), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = "AUC-PR skill (× chance)",
                  # The two bars are scored against DIFFERENT base rates (each set against its
                  # own), so their difference is not purely selection optimism — part of it is
                  # the outer window being an easier or harder problem. State that on the
                  # figure, with the common-denominator gap when the CSV carries it, rather
                  # than letting the bar difference be read as the optimism outright.
                  caption = local({
                    .g  <- suppressWarnings(as.numeric(r$optimism_gap))
                    .gc <- if ("optimism_gap_common_base" %in% names(r))
                             suppressWarnings(as.numeric(r$optimism_gap_common_base)) else NA_real_
                    .br <- if ("base_rate_ratio" %in% names(r))
                             suppressWarnings(as.numeric(r$base_rate_ratio)) else NA_real_
                    .txt <- if (is.finite(.g)) sprintf("Gap %.1f×", .g) else "Gap n/a"
                    if (is.finite(.gc)) .txt <- sprintf("%s (%.1f× on a common base rate)", .txt, .gc)
                    if (is.finite(.br)) .txt <- sprintf(
                      "%s. Each bar is scored against its own base rate; the held-out base rate is %.2f× the inner one.",
                      .txt, .br) else .txt <- paste0(.txt, ".")
                    stringr::str_wrap(.txt, width = 58)
                  })) +
    .rv_theme() +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, colour = .RV_MUTED, hjust = 0))
  .rv_save(p, out, 3.8, 4.4)
}

#' Driver: render all four review figures from the pipeline's diagnostics CSVs.
#' @param diag_dir directory holding the review CSVs (OUT_DIAGNOSTICS).
#' @param fig_dir  output directory (default the key_outputs manuscript panels dir).
#' @param cal_in_large optional calibration-in-the-large for the reliability subtitle.
#' @param horizons forecast horizons to render the reliability panel for. FigR1 is a
#'   retained deliverable at BOTH h=1 and h=2; the horizon used to be hard-coded to 1 on
#'   the input basename AND the output basename, so the h=2 panel was unreachable.
#' @param cal_in_large,cal_in_large_raw calibration-in-the-large for the reliability
#'   subtitle: either a scalar applied to every horizon, or a vector NAMED by horizon,
#'   e.g. c(`1` = 0.98, `2` = 1.07).
make_review_figures <- function(diag_dir, fig_dir, cal_in_large = NULL, top_k = 15L,
                                cal_in_large_raw = NULL,
                                horizons = get0("LFO_HORIZONS", ifnotfound = c(1L, 2L))) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  message("[review-fig] rendering review-response figures (house style) ...")
  # The calibration panels take BOTH scales: the primary file (recalibrated) and its raw
  # twin under diagnostics/raw/, written by the same block in run_all.R.
  .raw <- function(f) file.path(diag_dir, "raw", f)
  # This horizon's calibration-in-the-large from a horizon-named vector, else the scalar.
  # NULL for absent/non-finite, which plot_reliability_review() treats as "not supplied".
  .cal_for <- function(x, h) {
    if (is.null(x) || !length(x)) return(NULL)
    # A NAMED vector that lacks this horizon means the value was not computed for it; falling
    # back to x[[1]] would label the h=2 panel with h=1's calibration-in-the-large.
    if (!is.null(names(x))) {
      if (!as.character(h) %in% names(x)) return(NULL)
      v <- x[[as.character(h)]]
    } else v <- x[[1]]
    if (length(v) != 1L || !is.finite(v)) NULL else v
  }
  for (.h in horizons) {
    .rel <- sprintf("forecast_reliability_h%d.csv", .h)
    tryCatch(plot_reliability_review(file.path(diag_dir, .rel),
               file.path(fig_dir, sprintf("FigR1_reliability_h%d", .h)), .h,
               .cal_for(cal_in_large, .h), csv_raw = .raw(.rel),
               cal_in_large_raw = .cal_for(cal_in_large_raw, .h)),
             # `error =` is REQUIRED. Passed positionally, tryCatch() rejects the handler with
             # "condition handlers must be specified with a condition class" BEFORE evaluating
             # the expression — which unwound make_review_figures() on the first horizon and
             # silently lost FigR1 (both horizons), FigR3 and FigR4.
             error = local({ .hh <- .h; function(e)
               warning(sprintf("[review-fig] reliability h%d: %s", .hh, conditionMessage(e))) }))
  }
  tryCatch(plot_prospective_ranks_review(file.path(diag_dir, "prospective_invasion_per_zone.csv"),
             file.path(fig_dir, "FigR3_prospective_ranks"), top_k),
           error = function(e) warning("[review-fig] prospective: ", conditionMessage(e)))
  tryCatch(plot_heldout_optimism_review(file.path(diag_dir, "invasion_heldout_optimism.csv"),
             file.path(fig_dir, "FigR4_heldout_optimism")),
           error = function(e) warning("[review-fig] held-out: ", conditionMessage(e)))
  invisible(TRUE)
}
