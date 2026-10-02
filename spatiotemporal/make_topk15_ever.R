# ---------------------------------------------------------------------------
# Top-15 per-forecast-round outcome figure, THREE-WAY outcome:
#   invaded within the forecast window | invaded later | not invaded.
# Mirrors panel E of Figure2_labelled (build_topk_folds) but at top-15, and splits
# the "invaded" category into within-window vs eventually. House aesthetic copied
# verbatim from make_publication_figures.R (theme_pub / save_dual / palette).
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

# --- house-style constants (verbatim from make_publication_figures.R) ---
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
base_family <- "sans"

theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    plot.title      = element_blank(),
    plot.subtitle   = element_blank(),
    plot.caption    = element_blank(),
    axis.title      = element_text(size = base - 0.4, colour = MUTED),
    axis.title.x    = element_text(margin = margin(t = 4)),
    axis.title.y    = element_text(margin = margin(r = 4), angle = 90),
    axis.text       = element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
    legend.position = "top", legend.justification = "left",
    legend.title    = element_text(size = base - 1.2, colour = MUTED),
    legend.text     = element_text(size = base - 1.4, colour = INK),
    legend.key.height = unit(9, "pt"), legend.key.width = unit(15, "pt"),
    strip.text      = element_text(size = base - 0.6, colour = INK, face = "bold",
                                   margin = margin(3, 3, 3, 3)),
    plot.tag        = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin     = margin(6, 8, 6, 6)
  )
}

# PANEL_DIR is derived from OUT (resolved below), NOT from a path relative to the working
# directory. It was the literal "outputs/key_outputs/figures/panels", and run_all.R launches
# this script from the REPO ROOT — where an unrelated `outputs/` tree also exists — so every
# run wrote the panel into that stray tree while the copy in spatiotemporal/outputs went
# stale. The same resolution problem was already fixed for the INPUT paths just below; the
# output path was missed.
save_dual <- function(p, name, w, h, dir = PANEL_DIR) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  # Gate on the FULL destination path, not the bare stem: FIGURE_DROP entries are
  # "<directory>/<stem>" and the raw/ exclusion inspects path components, neither of
  # which can match a basename.
  if (is.function(.fk) && !.fk(file.path(dir, name))) return(invisible(p))
  p <- fs_caption(p)                       # scale statement (raw pass by default)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  message(sprintf("  saved %-34s  %.1f x %.1f in", name, w, h))
  invisible(p)
}
# tidytext-style within-facet ordering (verbatim)
reorder_within <- function(x, by, within) {
  factor(paste(x, within, sep = "___"),
         levels = unique(paste(x, within, sep = "___"))[order(within, by)])
}
tidytext_scale_y <- function() scale_y_discrete(labels = function(z) sub("___.*$", "", z))

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

# Load the retained-figure allow-list (FIGURE_KEEP / figure_is_kept, 00_config.R). Without it
# the gate in save_dual() below resolves figure_is_kept to NULL and silently no-ops, so
# FIGURE_DROP entries and the raw/ exclusion could never fire here even though this script
# writes a retained deliverable (F_topk15_ever_h1/_h2). dirname(OUT) is the spatiotemporal
# directory (OUT is .../spatiotemporal/outputs). Non-fatal: it must still run standalone.
if (!is.function(get0("figure_is_kept")) || !is.function(get0("lfo_origin")))
  try(suppressWarnings(suppressMessages(source(file.path(dirname(OUT), "00_config.R")))),
      silent = TRUE)

PANEL_DIR <- fs_out_dir(file.path(OUT, "key_outputs", "figures", "panels"))

# FEATURED MODEL read from the pipeline's own selection record rather than hard-coded — a
# literal here silently diverges from the model the pipeline actually chose (it said
# "Bayes-M17-med" while the manuscript still described Bayes-M14-med). Falls back to the
# previous literal only if the selection file is unreadable.
BEST_BAYES <- local({
  f <- file.path(OUT, "key_outputs", "model_selection.json")
  # $featured$bayesian, NOT $featured$headline. `headline` is the best model over ALL families
  # (write_model_selection.R), so it can be a renewal or baseline method — this figure is a
  # BAYESIAN product and would then plot a different model from the rest of the suite while its
  # caption still said "featured Bayesian model". `headline` is kept only as the fallback, which
  # is the convention every other consumer uses (make_manuscript_figures.R, 30_projection_config.R).
  sel <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
  m <- tryCatch(sel$featured$bayesian$method, error = function(e) NULL)
  if (is.null(m) || !length(m) || is.na(m[1]))
    m <- tryCatch(sel$featured$headline$method, error = function(e) NULL)
  if (is.null(m) || !length(m) || is.na(m[1])) {
    warning("[topk15] model_selection.json unreadable; falling back to the hard-coded featured model.",
            call. = FALSE); "Bayes-M14-fill-med"
  } else m[1]
})
message("[topk15] outputs = ", OUT, " | featured model = ", BEST_BAYES)
rs  <- read_csv(fs_risk_csv(OUT, sub = "reports"), show_col_types = FALSE)
lfo <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))
fs_apply_lfo_scale(lfo)

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
  # ASCII hyphen, not an em-dash: these facet strips are drawn on the base pdf device
  # (ISOLatin1 encoding), which cannot represent U+2014 and silently substitutes a hyphen
  # with a conversion warning — so the PDF and the PNG would otherwise disagree.
  lab_lk   <- setNames(sprintf("Round %d - %s", seq_along(fold_ord), format(lfo_origin(fold_ord), "%d %b")),
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

  # --- accuracy self-check: per-fold counts must reproduce the validated table ---
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
    labs(x = "Predicted invasion probability, P(first case)", y = NULL) +
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
