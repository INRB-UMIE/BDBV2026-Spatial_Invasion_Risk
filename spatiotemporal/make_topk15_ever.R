# ---------------------------------------------------------------------------
# Top-15 per-forecast-round outcome figure, THREE-WAY outcome:
#   invaded within the forecast window | invaded later | not invaded.
# Mirrors panel E of Figure2_labelled (build_topk_folds) but at top-15, and splits
# the "invaded" category into within-window vs eventually. House style (theme_pub /
# save_dual / palette / reorder_within) is sourced from figure_style.R, as in
# make_publication_figures.R.
#
# 2026-09-17: produced at BOTH horizons (F_topk15_ever_h1 / _h2). The horizon was
# hard-coded to 1 in the data filter AND in the output basename, so the h=2 panel was
# unreachable. The "31 Jul" snapshot date in the legend was hard-coded too — it is now
# read from run_info.json$training_window_end, as make_manuscript_figures.R already does,
# because a stale date on a published legend misstates what the outcome classes mean.
# ---------------------------------------------------------------------------
suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr); library(readr); library(scales)
}))

# --- data ---
# ANCHORED path, not the bare relative "outputs": this script used to assume the working
# directory was spatiotemporal/, which held only while it was run by hand from there. Now
# that run_all.R schedules it as a child process from the REPO ROOT it failed outright with
# "'outputs/reports/bayes_risk_scores_all_zones.csv' does not exist". Resolve against the
# script's own location so it works from any working directory.
OUT <- local({
  cand <- c(file.path(tryCatch(here::here(), error = function(e) "."), "spatiotemporal", "outputs"),
            file.path(".", "outputs"), "outputs")
  hit <- cand[dir.exists(cand)]
  if (!length(hit)) stop("[topk15] cannot locate the spatiotemporal outputs directory")
  hit[1]
})
# Probability-scale switch (forecast_scale.R): default "recalibrated" (the primary set),
# FORECAST_SCALE=raw writes the twin to a raw/ sibling under the same basenames.
source(file.path(dirname(OUT), "forecast_scale.R"))
# House style (palette, province scale, themes, save_dual) from figure_style.R — the
# single definition for the whole suite.
source(file.path(dirname(OUT), "figure_style.R"))


# The pipeline config is REQUIRED (00_config.R; dirname(OUT) is the spatiotemporal directory).
# It defines the script guard below (script_run_or_forced), the forecast-origin convention
# (lfo_origin) and the retained-figure allow-list (FIGURE_KEEP / figure_is_kept) that save_dual()
# gates on. It used to be loaded inside try(silent = TRUE) as "non-fatal", but nothing below can
# run without it: a failed load ended at the guard as "could not find function", and had the
# guard been defined elsewhere, save_dual() would have resolved figure_is_kept to NULL and
# skipped the gate (FIGURE_DROP and the raw/ exclusion) without a word.
if (!is.function(get0("figure_is_kept")) || !is.function(get0("lfo_origin")) ||
    !is.function(get0("script_run_or_forced")))
  suppressWarnings(suppressMessages(source(file.path(dirname(OUT), "00_config.R"))))
stopifnot(is.function(get0("figure_is_kept")), is.function(get0("lfo_origin")),
          is.function(get0("script_run_or_forced")))
# OUT is the PUBLISHED tree (above): a data-definition sensitivity must not redraw its panels.
stop_if_data_definition_sensitivity("make_topk15_ever.R")

# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("make_topk15_ever.R")) {
  message("[topk15] sourced, not run: nothing is written. Use `Rscript make_topk15_ever.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

# PANEL_DIR is derived from OUT, NOT from a path relative to the working directory. It was the
# literal "outputs/key_outputs/figures/panels", and run_all.R launches this script from the REPO
# ROOT — where an unrelated `outputs/` tree also exists — so every run wrote the panel into that
# stray tree while the copy in spatiotemporal/outputs went stale.
PANEL_DIR <- fs_out_dir(file.path(OUT, "key_outputs", "figures", "panels"))

# FEATURED MODEL read from the pipeline's own selection record rather than hard-coded — a
# literal here silently diverges from the model the pipeline actually chose (it said
# "Bayes-M17-med" while the manuscript still described Bayes-M14-med). An unreadable
# selection file stops the build.
# read_featured_model() (figure_style.R): $featured$bayesian, else $featured$headline, and NO
# FALLBACK - this used to fall back to a typed "Bayes-M14-fill-med" when the file was unreadable.
BEST_BAYES <- read_featured_model(file.path(OUT, "key_outputs", "model_selection.json"))
message("[topk15] outputs = ", OUT, " | featured model = ", BEST_BAYES)
rs  <- read_csv(fs_risk_csv(OUT, sub = "reports"), show_col_types = FALSE)
lfo <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))
# The column actually drawn is kept: a round with no fitted recalibration factor is labelled
# "(uncorrected)" only on a pass that draws the recalibrated column (build_topk_ever()).
PSCALE_COL <- fs_apply_lfo_scale(lfo)

# Last calendar day of data used for fitting. The outcome classes ("invaded later", "not
# invaded") are judged against THIS date, so the legend must quote it rather than a literal.
DATA_END_DATE <- {
  ri <- jsonlite::fromJSON(file.path(OUT, "key_outputs", "run_info.json"))
  d  <- if (is.null(ri$training_window_end)) as.Date(NA) else as.Date(ri$training_window_end)
  if (length(d) != 1L || is.na(d))
    stop("[topk15] run_info.json lacks training_window_end (last data day used for fitting)")
  d
}
DATA_END_LBL <- format(DATA_END_DATE, "%d %b")

# zones invaded AT SOME POINT (recorded a first confirmed case by the data-end date)
affected_ever <- rs %>%
  filter(horizon == 1, as.logical(was_active_before) %in% TRUE) %>%
  pull(health_zone) %>% unique()
message(sprintf("[topk15-ever] featured=%s | zones invaded by %s: %d",
                BEST_BAYES, DATA_END_LBL, length(affected_ever)))

TOP_K <- 15L
OUT_WIN   <- "#B33005"   # invaded within the forecast window (dark burnt orange)
OUT_LATER <- "#F6BB6B"   # invaded later                      (light amber)
OUT_NEVER <- "grey80"    # not invaded by the data-end date

#' Build and save the top-15 outcome panel for one forecast horizon.
#'
#' @param hz forecast horizon in weeks (1 or 2). Drives the data filter, the legend wording
#'   (a 1-week horizon resolves within a WEEK, a 2-week one within a FORTNIGHT) and the
#'   output basename. Facet rows and canvas height follow the fold count, which differs by
#'   horizon — the old fixed nrow = 3 / 12.0 in canvas assumed h = 1's nine folds.
build_topk_ever <- function(hz) {
  lv <- c(sprintf("Invaded within the forecast %s", if (hz == 1L) "week" else "fortnight"),
          sprintf("Invaded later (by %s)", DATA_END_LBL),
          sprintf("Not invaded (by %s)", DATA_END_LBL))
  out_col <- setNames(c(OUT_WIN, OUT_LATER, OUT_NEVER), lv)

  d <- lfo %>%
    filter(method == BEST_BAYES, horizon == hz, is.finite(p_invasion),
           !(as.logical(was_active_before) %in% TRUE)) %>%
    mutate(cutoff = as.Date(cutoff))
  if (!nrow(d)) {
    warning(sprintf("[topk15-ever] no rows for %s at h=%d; panel skipped.", BEST_BAYES, hz),
            call. = FALSE)
    return(invisible(NULL))
  }

  fold_ord <- sort(unique(d$cutoff))
  # WHICH ROUNDS ARE DRAWN UNCORRECTED. The bars are on the prequentially recalibrated scale,
  # but the earliest rounds have no settled training fold to fit delta on (delta_estimable
  # FALSE) and their p_recal is the RAW probability, drawn beside recalibrated rounds with
  # nothing to say so. Read per round off the featured model's own rows at this horizon (the
  # rows drawn), exactly as make_manuscript_figures.R's .est_by_cut/.uncorr() do for Figure 2B
  # and FigureS14. On a raw pass nothing is recalibrated, so nothing is "uncorrected" relative
  # to the rest; a frame without the column predates recalibration (one scale throughout).
  unc_cuts <- if (identical(PSCALE_COL, "p_recal") && "delta_estimable" %in% names(d)) {
    est <- tapply(d$delta_estimable %in% TRUE, as.character(d$cutoff), any)
    as.Date(names(est)[!est])
  } else as.Date(character(0))
  # WHICH ROUNDS ARE PROVISIONAL. A round whose outcome window closed fewer than
  # LFO_RELIABLE_EVAL_AGE_DAYS before the analysis date (eval_reliable FALSE, written by
  # run_invasion_lfo()) is scored, but its invasions are not all reported yet: a bar drawn as
  # "invaded later" or "not invaded" there may yet prove to have been invaded within the window.
  # Drawn unmarked it read as final; FigureS7 draws the same rounds hollow, and
  # build_topk_folds() labels them the same way. A frame without the column cannot say which
  # rounds have settled, so it is refused rather than drawn as if every outcome were final.
  if (!"eval_reliable" %in% names(d))
    stop("[topk15-ever] lfo_results.rds has no eval_reliable column, so unsettled rounds cannot ",
         "be marked provisional; re-run run_all.R.", call. = FALSE)
  rel <- tapply(d$eval_reliable %in% TRUE, as.character(d$cutoff), all)
  prov_cuts <- as.Date(names(rel)[!rel])
  round_tag <- vapply(seq_along(fold_ord), function(i) {
    tg <- c(if (fold_ord[i] %in% unc_cuts)  "uncorrected",
            if (fold_ord[i] %in% prov_cuts) "provisional")
    if (length(tg)) paste0(" (", paste(tg, collapse = ", "), ")") else ""
  }, character(1))
  if (any(fold_ord %in% unc_cuts))
    message(sprintf("[topk15-ever] h=%d: %d of %d round(s) carry no fitted recalibration factor and are labelled \"(uncorrected)\".",
                    hz, sum(fold_ord %in% unc_cuts), length(fold_ord)))
  if (any(fold_ord %in% prov_cuts))
    message(sprintf("[topk15-ever] h=%d: %d of %d round(s) have not settled (eval_reliable FALSE) and are labelled \"(provisional)\".",
                    hz, sum(fold_ord %in% prov_cuts), length(fold_ord)))
  # ASCII hyphen, not an em-dash: save_dual() now draws through figure_pdf_device(), which
  # encodes U+2014, but the last-resort base pdf device (ISOLatin1) cannot, and a hyphen reads
  # the same on every device.
  lab_lk   <- setNames(sprintf("Round %d - %s%s", seq_along(fold_ord),
                               format(lfo_origin(fold_ord), "%d %b"), round_tag),
                       as.character(fold_ord))

  topk <- d %>%
    group_by(cutoff) %>%
    slice_max(p_invasion, n = TOP_K, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      outcome = factor(dplyr::case_when(
                  is_new_invasion == 1           ~ lv[1],
                  health_zone %in% affected_ever ~ lv[2],
                  TRUE                           ~ lv[3]),
                levels = lv),
      fold_lab = factor(lab_lk[as.character(cutoff)], levels = lab_lk[as.character(fold_ord)]),
      zone_w   = reorder_within(health_zone, p_invasion, fold_lab))

  # --- per-round counts, LOGGED for checking against the evaluation tables (nothing is
  # asserted here: this prints them, it does not compare them) ---
  chk <- topk %>% group_by(fold_lab) %>%
    summarise(n = n(), within = sum(outcome == lv[1]), ever = sum(outcome != lv[3]),
              .groups = "drop")
  message(sprintf("[topk15-ever] h=%d per-fold (n / within-window / ever-invaded):", hz))
  print(as.data.frame(chk))
  message(sprintf("[topk15-ever] h=%d POOLED top-%d: ever=%d/%d (%.0f%%), within=%d/%d (%.0f%%)",
                  hz, TOP_K,
                  sum(topk$outcome != lv[3]), nrow(topk), 100 * mean(topk$outcome != lv[3]),
                  sum(topk$outcome == lv[1]), nrow(topk), 100 * mean(topk$outcome == lv[1])))

  # Facet rows / canvas scale with the number of folds rather than assuming h = 1's count.
  n_folds <- length(fold_ord)
  n_row   <- max(1L, ceiling(n_folds / 3L))
  fig_h   <- max(4.5, 4.0 * n_row)

  p <- ggplot(topk, aes(p_invasion, zone_w, fill = outcome)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    facet_wrap(~ fold_lab, scales = "free_y", nrow = n_row) +
    tidytext_scale_y() +
    scale_fill_manual(values = out_col, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.06)), breaks = scales::pretty_breaks(3)) +
    # The bars' scale, as Figure 2 panel E and FigureS14 state it.
    labs(x = paste0("Predicted invasion probability, P(first case)",
                    if (!identical(PSCALE_COL, "p_recal")) " (raw)"
                    else if (any(fold_ord %in% unc_cuts)) " (recalibrated; rounds marked uncorrected are raw)"
                    else " (recalibrated)"), y = NULL) +
    theme_pub(8.6) +
    theme(panel.grid.major.y = element_blank(),
          panel.spacing.x = unit(9, "pt"), panel.spacing.y = unit(8, "pt"),
          strip.clip = "off",
          axis.text.y = element_text(size = 6.6, colour = INK),
          legend.position = "top")

  save_dual(p, sprintf("F_topk15_ever_h%d", hz), 9.6, fig_h)
  invisible(p)
}

for (.hz in c(1L, 2L)) build_topk_ever(.hz)
message("[topk15-ever] done.")

}   # end of the script-run guard

