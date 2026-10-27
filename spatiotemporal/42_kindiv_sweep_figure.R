# =============================================================================
# 42_kindiv_sweep_figure.R — publication figure for the comprehensive
# individual-offspring overdispersion (k_indiv) sensitivity sweep.
#
# Reads the k_indiv sweep (written by run_cascade.R's sensitivity stage, or by
# 41_kindiv_sweep.R) and renders one 3-panel figure in the house publication
# aesthetic (Okabe-Ito, theme_pub, dual PDF + 600-dpi PNG):
#   A) Ranking stability vs k     — Spearman rho, Kendall tau, top-15 Jaccard
#   B) Magnitude vs k             — mean reach at 4/8/13 wk + mean establishment
#   C) Per-zone reach vs k        — top-12 at-risk zones, coloured by province
# The sweep's own baseline k (k_base in its sidecar) is marked with a dashed reference
# line. (The Lloyd-Smith [0.2, 0.4] plausible-band box was removed by request.) Titles
# are blank (house style: captions live externally); panels are tagged A/B/C.
# Drawn ONLY from a sweep of the latest cascade run for this analysis date, at this run's
# kernel, form, scenario, seed and k_base (PROVENANCE below); anything else is skipped.
#
# Run:  Rscript spatiotemporal/42_kindiv_sweep_figure.R   (run_cascade.R source()s it)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(scales); library(forcats); library(jsonlite); library(here)
})

# CONFIGURATION FIRST, THEN THE PATHS. The paths used to be fixed from the working directory
# BEFORE 00_config.R was sourced, so OUT_DIR — and with it a data-definition sensitivity's own
# tree (DATA_DEFINITION_TAG) — was ignored, and the provenance check below compared against
# CASCADE_KERNEL / CASCADE_K_INDIV, which only 30_projection_config.R defines and which a
# standalone run therefore never had (so those comparisons were silently skipped). When
# run_cascade.R source()s this file, its configuration is already loaded and is NOT re-sourced:
# re-sourcing 30 there would reset the driver's globals mid-run.
ST <- Sys.getenv("CASCADE_ST_DIR",
                 unset = get0("ST_DIR", ifnotfound = file.path(here::here(), "spatiotemporal")))
if (!exists("OUT_DIR") || !is.function(get0("figure_is_kept")))
  source(file.path(ST, "00_config.R"))
if (!exists("CASCADE_KERNEL") || !exists("CASCADE_FORM") || !exists("OUT_CASCADE"))
  source(file.path(ST, "30_projection_config.R"))
# cascade_stamp_check() / cascade_calibration_note() (35 only defines functions).
if (!exists("cascade_stamp_check", mode = "function"))
  source(file.path(ST, "35_cascade_viz.R"))
# House style from figure_style.R, the single definition for the whole suite. This file's
# own PROV_COL permuted Sud-Kivu, Bas-Uele and Tshopo against the main-text figures'.
source(file.path(ST, "figure_style.R"))

# OUT_CASCADE and CASCADE_OUT_DIR follow CASCADE_OUT_ROOT (30_projection_config.R), as every
# other cascade product does, and OUT_DIR the data definition (00_config.R).
OUTC <- OUT_CASCADE
DIAG <- file.path(OUTC, "diagnostics")
FIGD <- file.path(OUTC, "figures")
KEYD <- file.path(get0("CASCADE_OUT_DIR", ifnotfound = OUT_DIR), "key_outputs")
# Every figure here is written into BOTH the cascade figures/ tree and key_outputs/
# (save_dual() creates the directories).
SAVE_DIRS <- c(FIGD, KEYD)
BAND_FILL <- "#0072B2"                       # Lloyd-Smith plausible band (soft, low alpha)
# plausible-band [0.2, 0.4] box removed (user request); no-op layer (NULL is ignored
# when added to a ggplot), so the panel call sites below are left unchanged.
band_layer <- function(lo = 0.20, hi = 0.40) NULL
kbase_layer <- function(kb)
  geom_vline(xintercept = kb, colour = FAINT, linetype = "22", linewidth = 0.35)
LOGX <- scale_x_log10(breaks = c(0.05,0.1,0.2,0.3,0.5,1,2,4),
                      labels = c("0.05","0.1","0.2","0.3","0.5","1","2","4"))

# A skip is an ordinary state here, not an error: the sweep is optional and expensive.
# NEVER quit() when sourced: run_cascade.R SOURCES this file, and quit() from a sourced script
# terminates the whole calling process — a missing optional input used to kill the entire
# cascade run WITH STATUS 0, silently skipping everything after it. Sourced, raise a condition
# the caller's tryCatch already handles as a skip; run as its own Rscript, exit 0, which is the
# correct status for "nothing to draw".
.skip42 <- function(msg) {
  message("[42_kindiv_sweep_figure] Skipping - ", msg)
  if (!is_script_run("42_kindiv_sweep_figure.R"))
    stop("[42_kindiv_sweep_figure] ", msg, "; figure skipped.", call. = FALSE)
  quit(save = "no", status = 0L)
}

# ---- data -------------------------------------------------------------------
# A missing input is an ordinary state (the sweep is not run on every pass), so it skips
# cleanly instead of failing on read_csv() — required since run_all.R schedules this script.
.sweep_csv  <- file.path(DIAG, "cascade_kindiv_sweep.csv")
.byzone_csv <- file.path(DIAG, "cascade_kindiv_sweep_byzone.csv")
.meta_json  <- file.path(DIAG, "cascade_kindiv_sweep_meta.json")
.missing <- c(.sweep_csv, .byzone_csv, .meta_json)[!file.exists(c(.sweep_csv, .byzone_csv, .meta_json))]
if (length(.missing))
  .skip42(paste0("no k_indiv sweep on disk (missing: ", paste(basename(.missing), collapse = ", "),
                 "). It is written by run_cascade.R's sensitivity stage or by ",
                 "`Rscript spatiotemporal/41_kindiv_sweep.R` (RUN_KINDIV_SWEEP=1 in run_all.R)"))
meta <- tryCatch(fromJSON(.meta_json), error = function(e) NULL)
if (is.null(meta)) .skip42(paste0(basename(.meta_json), " is unreadable, so the sweep's provenance cannot be checked"))

.m1 <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) NA_character_ else as.character(x)
# PROVENANCE. RUN_KINDIV_SWEEP defaults FALSE while RUN_CASCADE_FIGURES defaults TRUE, so this
# script runs every time and 41_kindiv_sweep.R usually does not. It used to skip only when the
# CSVs were ABSENT — never when they were merely OLD — and its "provenance OK" compared nothing
# that could fail standalone, so Figure_kindiv_sweep (a retained deliverable) was re-rendered
# under a current mtime from whatever sweep was on disk: one left by a cascade run that later
# failed at calibration (whose sweep reset never ran), or one of a model the suite no longer
# features (CASCADE_KERNEL / CASCADE_FORM follow model_selection.json). Now the sweep is drawn
# only when (1) its sidecar's cascade_run_id is the LATEST cascade run's for this analysis date
# and that run completed, or is the run_cascade.R process drawing it (cascade_stamp_check());
# (2) the sidecar's analysis date is this one; and (3) it was run at this configuration's
# kernel, hazard form, primary scenario, seed and baseline k. Otherwise it is skipped.
.prov <- cascade_stamp_check(.m1(meta$cascade_run_id), DIAG)
.mis <- character(0)
if (!isTRUE(.prov$ok)) .mis <- c(.mis, .prov$why)
.want_date <- format(as.Date(ANALYSIS_DATE))
if (!identical(.m1(meta$analysis_date), .want_date))
  .mis <- c(.mis, sprintf("analysis date (sweep %s vs this run %s)", .m1(meta$analysis_date), .want_date))
.cmp <- list(kernel = CASCADE_KERNEL, form = CASCADE_FORM, scenario = CASCADE_SCENARIO_PRIMARY)
for (.nm in names(.cmp))
  if (!identical(.m1(meta[[.nm]]), as.character(.cmp[[.nm]])))
    .mis <- c(.mis, sprintf("%s (sweep %s vs this run %s)", .nm, .m1(meta[[.nm]]), .cmp[[.nm]]))
if (!identical(suppressWarnings(as.integer(.m1(meta$seed))), as.integer(CASCADE_SEED)))
  .mis <- c(.mis, sprintf("seed (sweep %s vs this run %s)", .m1(meta$seed), CASCADE_SEED))
KB <- suppressWarnings(as.numeric(.m1(meta$k_base)))
if (!is.finite(KB) || abs(KB - CASCADE_K_INDIV) > 1e-9)
  .mis <- c(.mis, sprintf("baseline k (sweep %s vs this run CASCADE_K_INDIV = %s)",
                          .m1(meta$k_base), format(CASCADE_K_INDIV)))
if (length(.mis))
  .skip42(sprintf(paste0("the k_indiv sweep on disk is not this run's - %s. Re-run the cascade ",
                         "(its sensitivity stage writes the sweep) or 41_kindiv_sweep.R."),
                  paste(.mis, collapse = "; ")))
message(sprintf("[42_kindiv_sweep_figure] sweep provenance OK (run %s, kernel %s, form %s, scenario %s, seed %s, k_base %s).",
                .m1(meta$cascade_run_id), .m1(meta$kernel), .m1(meta$form), .m1(meta$scenario),
                .m1(meta$seed), format(KB)))
# NULL when the sweep's run was calibrated; otherwise the caption the figure carries (Panels B
# and C are absolute probabilities).
CAL_NOTE <- cascade_calibration_note(meta$calibrated)
sweep  <- read_csv(.sweep_csv,  show_col_types = FALSE)
byzone <- read_csv(.byzone_csv, show_col_types = FALSE)

# ---- Panel A: ranking stability --------------------------------------------
STAB_COL <- c("Spearman rank corr." = OKABE[1], "Kendall tau" = OKABE[2], "Top-15 Jaccard" = OKABE[3])
pa_df <- sweep %>%
  transmute(k_indiv,
            `Spearman rank corr.` = spearman, `Kendall tau` = kendall, `Top-15 Jaccard` = top15_jaccard) %>%
  pivot_longer(-k_indiv, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = names(STAB_COL)))
pA <- ggplot(pa_df, aes(k_indiv, value, colour = metric)) +
  band_layer() + kbase_layer(KB) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.25) +
  scale_colour_manual(values = STAB_COL, name = NULL) +
  LOGX + coord_cartesian(ylim = c(min(0.7, min(pa_df$value, na.rm = TRUE) - 0.02), 1.005)) +
  labs(x = expression("individual-offspring overdispersion  " * italic(k)),
       y = sprintf("agreement with k = %s ranking", format(KB, nsmall = 2))) +
  theme_pub()

# ---- Panel B: magnitude (reach at horizons + establishment) ----------------
MAG_COL <- c("Reach, 4 wk" = "#9ecae1", "Reach, 8 wk" = "#4292c6", "Reach, 13 wk" = "#08519c",
             "Establishment, 13 wk" = OKABE[2])
pb_df <- sweep %>%
  transmute(k_indiv,
            `Reach, 4 wk` = mean_reach_h4, `Reach, 8 wk` = mean_reach_h8,
            `Reach, 13 wk` = mean_reach_h13, `Establishment, 13 wk` = mean_estab_h13) %>%
  pivot_longer(-k_indiv, names_to = "series", values_to = "value") %>%
  mutate(series = factor(series, levels = names(MAG_COL)))
pB <- ggplot(pb_df, aes(k_indiv, value, colour = series)) +
  band_layer() + kbase_layer(KB) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.15) +
  scale_colour_manual(values = MAG_COL, name = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  LOGX + scale_y_continuous(labels = label_percent(accuracy = 1)) +
  expand_limits(y = 0) +
  labs(x = expression("individual-offspring overdispersion  " * italic(k)),
       y = "mean probability, at-risk zones") +
  theme_pub()

# ---- Panel C: per-zone reach trajectories (top-12 at-risk by baseline) ------
base_rank <- byzone %>% filter(abs(k_indiv - KB) < 1e-9) %>%
  arrange(desc(p_invasion)) %>% slice_head(n = 12) %>% pull(health_zone)
pc_df <- byzone %>% filter(health_zone %in% base_rank) %>%
  mutate(province = ifelse(is.na(province) | !province %in% names(PROV_COL), "Other", province),
         province = factor(province, levels = names(PROV_COL)))
pC <- ggplot(pc_df, aes(k_indiv, p_invasion, group = health_zone, colour = province)) +
  band_layer() + kbase_layer(KB) +
  geom_line(linewidth = 0.45, alpha = 0.9) + geom_point(size = 0.8, alpha = 0.9) +
  scale_colour_manual(values = PROV_COL, name = NULL, drop = TRUE) +
  LOGX + scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(x = expression("individual-offspring overdispersion  " * italic(k)),
       y = "13-week reach P, top-12 zones") +
  theme_pub()

# ---- compose ----------------------------------------------------------------
# An uncalibrated sweep is labelled on the figure (CAL_NOTE; no caption otherwise). theme_pub()
# blanks plot.caption, so the annotation sets it explicitly.
fig <- (pA | pB) / pC +
  plot_annotation(tag_levels = "A",
                  # wrapped: one line runs off a 7.2-inch figure
                  caption = if (is.null(CAL_NOTE)) NULL else paste(strwrap(CAL_NOTE, 130), collapse = "\n"),
                  theme = theme(plot.caption = element_text(size = 7, colour = MUTED, hjust = 0))) &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(fig, "Figure_kindiv_sweep", w = 7.2, h = 6.4)

message("[42_kindiv_sweep_figure] done.")
