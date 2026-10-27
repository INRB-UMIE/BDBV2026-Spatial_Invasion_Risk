# =============================================================================
# make_manuscript_figures.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
# Main-text manuscript figures (Nature-ready aesthetic), rebuilt from SAVED data
# (no model re-fitting), matching the house design system used across key_outputs.
#
#   Figure 1  Spatial invasion       A) bivariate choropleth: cumulative cases x
#                                        date of first onset (locator + 4x4 legend +
#                                        province outlines + Kisangani + arrival~burden rho)
#                                     B) epidemic trajectory in ONE dual-axis panel:
#                                        stacked cumulative cases (the EPICENTRE_ZONES
#                                        epicentre zones + others) on the left, health zones
#                                        invaded on the right
#                                     C) arrival time (x) vs road travel time + short-trip
#                                        epicentre outflow (R^2 + Spearman rho). All 5 predictors -> SI.
#   Figure S1 arrival time (x) vs all 5 predictors (-log short-trip and -log Flowminder
#                                        relocation share, months read from the table / config;
#                                        great-circle, road travel time,
#                                        log10 population). Panels C and S1 are drawn from
#                                        key_outputs/arrival_predictor{s,_fits}.csv, which
#                                        43_spread_kinematics.R computes; nothing is fitted here.
# All forecast-based panels are fixed to the 2-WEEK-AHEAD horizon (Figures 2 & 3).
# Figure 1 is descriptive (observed invasion), so it carries no forecast horizon.
#   Figure 2  Predictive performance  A) prioritisation (model vs naive vs random), 2 wk,
#                                        90% fold-cluster bootstrap band + discrimination and
#                                        calibration annotation
#                                     B) forecast-vs-outcome for 3 evenly-spread rounds, 2 wk
#                                        (among the rounds with an estimable recalibration
#                                        factor, so all three share one probability scale);
#                                        top-20 at-risk zones coloured invaded-within-round /
#                                        invaded-later / not-invaded (by the last data day used for fitting)
#                                     C) what the mobile-phone mobility data adds, 2 wk: the
#                                        AUC-PR skill ladder with 90% paired zone-cluster
#                                        bootstrap intervals, read from
#                                        key_outputs/mobility_value_ladder.csv (published by
#                                        50_mobility_value.R) and drawn through the SAME
#                                        builder as FigureS24 (mobility_value_panel.R).
#                                        Skipped on an MI pass; see build_fig2c.
#   Figure S14 panel B for EVERY round (built by build_fig2; was "FigureS5_topk_all_rounds",
#                                        renumbered because S5 is the ROC/PR figure)
#   Figure 3  Invasion + forecast     MAIN TEXT (2-wk horizon):
#                                     A) relative-invasion-risk map
#                                     B) top-20 forest (invasion probability + 90% CrI, by province)
#   Figure S4 Front & rank narrative  SI (split out of the former Figure 3 A/B):
#                                     A) front-approach: distance to nearest invaded zone per round,
#                                        evolving to the present (explains B; shared axis/colour/zones)
#                                     B) 2-wk invasion-risk rank evolution across rounds up to the
#                                        present live run (tracks main-text Fig 3B's top zones)
#   Figure3_rank_surveillance         Figure 3 extended (written ALONGSIDE Figure 3, not over it):
#                                     A, B) Figure 3A and 3B, unchanged
#                                     C) the Figure S4B rank evolution of the same top-20 zones
#                                     D) per zone, alerts and laboratory-resolved alerts
#                                        (confirmed + discarded) from the outbreak start to the
#                                        live run, via zone_alert_counts() — the definition the
#                                        ascertainment appendix (45) uses. Numbers also written
#                                        to Figure3_rank_surveillance_panelD.csv.
#   Figure S2 operational detail      A) prob x vulnerability bivariate  B) priority scatter
#
# Outputs -> outputs/key_outputs/manuscript_figures/{Figure1,2,3,...}.{pdf,png} (the composed
#            figures above), plus Figure2_simple_caption.md and Figure3_rank_surveillance_panelD.csv.
#            The single-panel stems (F1A_..., FS4B_...) are saved to panels/ through save_dual(),
#            whose FIGURE_KEEP gate (00_config.R) does not list them, so they are written only
#            under FIGURE_KEEP_ALL=1.
# Run:  Rscript make_manuscript_figures.R   (from any directory; paths are anchored with here())
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(sf); library(scales)
  library(viridisLite); library(lubridate); library(forcats)
})
sf::sf_use_s2(TRUE)
options(dplyr.summarise.inform = FALSE)
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("need ggrepel")

# Anchor to the repo via here::here() so the script runs from ANY working directory
# (previously normalizePath(".") assumed the wd was spatiotemporal/, so running from the
# repo root pointed DATA at <repo>/../data and every path escaped the repo).
HERE <- file.path(here::here(), "spatiotemporal")
DATA <- file.path(HERE, "..", "data")

# -----------------------------------------------------------------------------
# INPUT ROOT / OUTPUT SUBDIRECTORY — the onset multiple-imputation pass
# -----------------------------------------------------------------------------
# 48_onset_mi_pool.R writes MI-POOLED artifacts into a MIRROR of outputs/ (same relative
# paths, same schemas), so this suite can draw the pooled Figures 2 and 3 with no change
# to a single expression below: only the root it reads from moves.
#
#   MANU_IN_ROOT   read pipeline artifacts from here instead of <HERE>/outputs
#   MANU_FIG_SUB   write figures to key_outputs/manuscript_figures/<sub> instead of
#                  key_outputs/manuscript_figures
#   MANU_ONLY      comma-separated builders to run (default: all)
#
# THE GUARD IS THE POINT. Reading pooled inputs while writing to the canonical figure
# directory would silently overwrite the PUBLISHED single-draw Figures 2 and 3 with
# MI-pooled twins carrying the same basenames and captions that no longer describe them.
# That is the single worst thing this file can do, so a non-canonical input root without a
# figure subdirectory is refused outright rather than defaulted.
MANU_IN_ROOT <- trimws(Sys.getenv("MANU_IN_ROOT", ""))
MANU_FIG_SUB <- trimws(Sys.getenv("MANU_FIG_SUB", ""))
if (nzchar(MANU_IN_ROOT) && !nzchar(MANU_FIG_SUB))
  stop("[manu] MANU_IN_ROOT is set (reading non-canonical artifacts) but MANU_FIG_SUB is not, ",
       "so these figures would be written over the PUBLISHED ones under the same basenames. ",
       "Set MANU_FIG_SUB (e.g. 'mi') or unset MANU_IN_ROOT.", call. = FALSE)
OUT <- if (nzchar(MANU_IN_ROOT)) {
  if (!dir.exists(MANU_IN_ROOT))
    stop("[manu] MANU_IN_ROOT does not exist: ", MANU_IN_ROOT, call. = FALSE)
  normalizePath(MANU_IN_ROOT, mustWork = TRUE)
} else file.path(HERE, "outputs")

# MI mode = drawing from a pooled root. Several products change behaviour under it: Figure S4
# reads `zone_first`, which comes from a SINGLE load_linelist() draw, and pairing that
# single-draw front with an MI-pooled forecast would be a mixture of two different objects
# presented as one; Figure S2's pooling rule is unspecified; Figure 2 panel C, Figure2_simple
# and the all-rounds FigureS14 are single-draw ladders or round panels. All are skipped on an MI
# pass (see build_fig2 and build_fig3).
MI_MODE <- nzchar(MANU_IN_ROOT)
MANU_ONLY <- local({
  v <- trimws(Sys.getenv("MANU_ONLY", ""))
  if (!nzchar(v)) return(NULL)
  trimws(strsplit(v, ",", fixed = TRUE)[[1]])
})
# The builders drawn from a single load_linelist() draw, which an MI pass never builds. Gated
# HERE, not only when the run list is assembled at the end: .want() also decides which inputs
# are read at load time, and on an MI pass without MANU_ONLY it used to say yes to Figure 1 —
# so the arrival-predictor tables were read from the pooled mirror, which does not carry them,
# and the suite aborted before reaching the builders it was asked for.
.SINGLE_DRAW_BUILDERS <- c("Figure1", "FigureS1", "FigureS3")
.want <- function(nm) (is.null(MANU_ONLY) || nm %in% MANU_ONLY) &&
                      !(MI_MODE && nm %in% .SINGLE_DRAW_BUILDERS)
if (!is.null(MANU_ONLY)) {
  .known <- c("Figure1", "FigureS1", "FigureS3", "Figure2", "Figure3")
  .bad <- setdiff(MANU_ONLY, .known)
  if (length(.bad))
    stop("[manu] MANU_ONLY names unknown builder(s): ", paste(.bad, collapse = ", "),
         ". Known: ", paste(.known, collapse = ", "), call. = FALSE)
  message("[manu] MANU_ONLY: ", paste(MANU_ONLY, collapse = ", "))
}
if (MI_MODE)
  message(sprintf("[manu] MI MODE — inputs %s -> figures manuscript_figures/%s",
                  OUT, MANU_FIG_SUB))
# Probability-scale switch (forecast_scale.R). Default "recalibrated": the primary set.
# FORECAST_SCALE=raw re-runs the identical build onto a raw/ sibling, same basenames.
source(file.path(HERE, "forecast_scale.R"))
# FIGURES ALWAYS LAND IN THE CANONICAL TREE, under MANU_FIG_SUB when one is given. The
# figure directory deliberately does NOT follow MANU_IN_ROOT: the pooled input mirror holds
# tables, and scattering published PDFs across per-pass roots is how a figure comes to be
# cited from a directory nobody can find again.
.FIG_BASE <- file.path(HERE, "outputs", "key_outputs", "manuscript_figures")
if (nzchar(MANU_FIG_SUB)) .FIG_BASE <- file.path(.FIG_BASE, MANU_FIG_SUB)
FIG_DIR   <- fs_out_dir(.FIG_BASE)
PANEL_DIR <- fs_out_dir(file.path(.FIG_BASE, "panels"))
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# --- constants -----------------------------------------------------------------
# 01_data_prep.R sources 00_config.R, so OUTBREAK_START and EPICENTRE_ZONES come from
# the pipeline rather than being mirrored here. They used to be local copies "so the
# script is standalone"; a local copy of a constant is a second definition that can
# drift from the one the models use, and a figure that disagrees with the model it
# illustrates is worse than a figure that fails to build.
source(file.path(HERE, "01_data_prep.R"))
stop_if_data_definition_sensitivity("make_manuscript_figures.R")   # hard-coded published paths

# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("make_manuscript_figures.R")) {
  message("[manu] sourced, not run: nothing is written. Use `Rscript make_manuscript_figures.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

# invasion_common_cells(): the (fold x zone) support every method is scored on in
# invasion_evaluation.csv, resolved by the code evaluate_invasion() and
# compute_detection_curve() use. Figure 2A's Brier score and calibration, and the round counts
# its caption quotes, are taken on that support so they describe the rows the rest of the
# panel describes. The module defines functions only (it re-sources 00_config.R, which
# 01_data_prep.R has already loaded).
source(file.path(HERE, "16_invasion_eval.R"))

stopifnot(exists("OUTBREAK_START"), exists("EPICENTRE_ZONES"), exists("load_linelist"))
SHP_PATH       <- file.path(DATA, "shapefiles", "DRC_Health_zones.shp")
PROV_INT       <- c("Ituri", "Nord-Kivu", "Haut-Uele")
EPI_ZONES      <- EPICENTRE_ZONES                        # M+B+R epicentre
# ARRIVAL_REF is read from arrival_predictors.csv below (the observed index onset).

# -----------------------------------------------------------------------------
# 1. DESIGN SYSTEM  (colourblind-safe)
# -----------------------------------------------------------------------------
# The palette, theme_pub()/theme_map() and save_dual() live in figure_style.R so that the
# SI figure suites are drawn in the SAME style as the main-text figures by construction
# rather than by a second copy that has to be kept in step by hand.
source(file.path(HERE, "figure_style.R"))
# Figure 2 panel C is the mobility-value ladder that FigureS24 also draws. Its row set, its
# two encodings, its ordering rule and the panel itself live in ONE place, shared with
# 50_mobility_value.R, for the same reason the style does: a main-text panel and its appendix
# twin built by two pieces of code can disagree, and nothing in the suite would notice.
# The file defines objects only - it starts no work and does not touch PANEL_DIR.
source(file.path(HERE, "mobility_value_panel.R"))

# --- MI STAMP ----------------------------------------------------------------------------
# An MI-pooled figure shares its BASENAME with the published single-draw twin, so once the
# PDF leaves manuscript_figures/mi/ nothing on the page distinguishes the two — and their
# intervals mean different things (a posterior credible interval against a multiple-
# imputation mixture). Same hazard, and the same remedy, as fs_stamp()'s raw/ twins: stamp
# the figure itself, not only the path. save_dual() (figure_style.R) applies it to every panel
# and every composed figure on an MI pass through .fig_mi_stamp(), which reads MI_STAMP_TXT and
# APPENDS to a caption the figure already shows. (A local save_dual wrapper used to set the
# caption outright first: it deleted a figure's own caption and, with scale stamps on, printed
# the MI stamp twice.)
MI_META <- NULL
if (MI_MODE) {
  MI_META <- tryCatch({
    ri <- jsonlite::fromJSON(file.path(OUT, "key_outputs", "run_info.json"))
    if (is.null(ri$onset_mi)) stop("run_info.json carries no `onset_mi` block")
    ri$onset_mi
  }, error = function(e)
    stop("[manu] MI mode but the pooled run_info.json does not describe the imputation set (",
         conditionMessage(e), "). Refusing to stamp a figure with provenance it cannot ",
         "verify — re-run 48_onset_mi_pool.R.", call. = FALSE))
  MI_STAMP_TXT <- sprintf(paste0("Onset multiple imputation: M = %d completions (seeds %d-%d, ",
                                 "stratum '%s'). Intervals are MI mixtures, NOT single-draw ",
                                 "posterior credible intervals."),
                          as.integer(MI_META$n_imputations), as.integer(MI_META$seed_min),
                          as.integer(MI_META$seed_max), as.character(MI_META$stratum))
  message("[manu] MI stamp: ", MI_STAMP_TXT)
}

# --- bivariate-choropleth helpers (shared by Figures 1 and 3) --------------------------------
make_bivar <- function(c00, c10, c01, c11) {
  cc <- function(h) grDevices::col2rgb(h)[, 1]
  a <- cc(c00); b <- cc(c10); d <- cc(c01); e <- cc(c11)
  g <- expand.grid(bx = 1:4, by = 1:4)
  cols <- apply(g, 1, function(r) { fx <- (r[1]-1)/3; fy <- (r[2]-1)/3
    v <- (1-fx)*(1-fy)*a + fx*(1-fy)*b + (1-fx)*fy*d + fx*fy*e
    grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255) })
  setNames(cols, paste0(g$bx, "-", g$by))
}
BIV_ARRIVAL <- make_bivar("#e8e8e8", "#c85a5a", "#5a8fc8", "#3a2f4f")  # x = cases, y = later arrival
BIV_PROBVUL <- make_bivar("#e8e8e8", "#5ac8c8", "#be64ac", "#3b4994")  # x = invasion prob, y = vulnerability
bivar_legend <- function(BIP, xlab, ylab, title = NULL, base = 9) {
  leg_df <- expand.grid(bx = 1:4, by = 1:4); leg_df$bikey <- paste0(leg_df$bx, "-", leg_df$by)
  p <- ggplot(leg_df, aes(bx, by, fill = bikey)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = BIP, guide = "none") + coord_fixed() +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal(base_size = base, base_family = base_family) +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          axis.title = element_text(size = base - 0.8, colour = MUTED, lineheight = 0.9),
          plot.title = element_text(size = base - 0.6, colour = INK, face = "bold", hjust = 0.5, margin = margin(b = 2)),
          plot.margin = margin(2,2,2,2))
  p
}

# -----------------------------------------------------------------------------
# 2. LOAD SAVED DATA
# -----------------------------------------------------------------------------
message("[load] shapefile, risk scores, LFO, evaluation, line list ...")
shp <- st_read(SHP_PATH, quiet = TRUE) %>% mutate(.key = key(Nom), .prov = as.character(PROVINCE))
rs  <- read_csv(fs_risk_csv(OUT), show_col_types = FALSE)
lfo <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))
ev  <- read_csv(file.path(OUT, "diagnostics", "invasion_evaluation.csv"), show_col_types = FALSE)
# Put the selected scale into `p_invasion` ONCE, keeping the untouched values as
# `p_invasion_raw`. Every expression below reads `p_invasion`, so the switch is total
# rather than depending on a dozen call sites being edited consistently.
# Capture WHICH column was actually resolved. The scale-dependent metrics quoted on Figure 2
# (Brier, calibration) must come from the same scale as the probabilities being drawn, and
# fs_lfo_col() can fall back to raw; keying off the requested scale instead of the resolved
# one is how a panel comes to quote a recalibrated number over raw probabilities.
PSCALE_COL <- fs_apply_lfo_scale(lfo)
# Discrimination is quoted from the RAW columns in BOTH passes, deliberately. Within a
# fold the recalibration is strictly monotone, so AUC is EXACTLY invariant; the `*_recal`
# columns differ only because POOLED AUC mixes folds carrying different prequential
# factors. Publishing a "recalibrated AUC" that differs from the raw one would invite the
# reader to conclude that recalibration changed discrimination, which is false. Brier, log
# score and reliability DO change with the scale and are computed from the selected
# probabilities wherever they appear.
# AUC-PR SKILL IS PER ROUND, OVER ALL ROUNDS: AP within each round, event-weighted, over the base
# rate (auc_pr_skill_date), with its random-ordering level (auc_pr_skill_date_chance) and the
# all-round base rate. Performance is REPORTED over every round; the selection reads the
# post-burn-in rounds by default (INVASION_SELECTION_BASIS, 16_invasion_eval.R), a deliberate
# choice the caption states beside the selection-round value (sel_auc_pr_skill_date). The pooled
# auc_pr_skill ranks every round's rows in one list and rewards a model's probability LEVEL on
# the rounds that had invasions, not its ordering within a round (.auc_pr_by_date()); quoting it
# gave Figure 2 the opposite conclusion against the comparators.
EV_SKILL  <- "auc_pr_skill_date"
EV_CHANCE <- "auc_pr_skill_date_chance"
EV_BASE   <- "base_rate"
EV_ROC   <- "auc_roc"
# CALIBRATION IS SCALE-DEPENDENT, unlike AUC: a monotone recalibration leaves the ranking (and
# so AUC) untouched but moves every calibration measure, which is the whole point of applying
# it. This column therefore tracks the column actually drawn, where the discrimination columns
# above deliberately do not. It pools EVERY round, so Figure 2A does not print it: the block
# quotes calibration on the recalibrated rounds only (fig2_calibration_metrics()). It is the
# cross-check that the block's rows are the evaluation's rows — the same estimator over all
# rounds must reproduce it.
EV_CAL <- if (identical(PSCALE_COL, "p_recal")) "calibration_in_large_recal" else "calibration_in_large"
# EV_ECE is GONE from this suite. Figure 2A's annotation block quoted the expected calibration
# error beside calibration-in-the-large, and the block had grown to six lines inside the panel;
# the base rate, which is the denominator of the AUC-PR skill on the line above it, earns that
# line more. Nothing is lost: ECE is a one-number summary of a reliability curve, and the curve
# itself is FigureS7_calibration panel C (make_si_model_figures.R, which bins the LFO frame
# directly). The columns are still computed and still published in
# diagnostics/invasion_evaluation.csv (16_invasion_eval.R) for the appendix table.

#' One value out of invasion_evaluation.csv, by column / horizon / method.
#'
#' NAME THE MISSING COLUMN. `ev[[col]]` on an absent column returns NULL, which sprintf
#' renders as "character(0)" and silently drops the whole annotation string - a published
#' panel would then lose its metric block with no error anywhere. And a column that is present
#' but carries two rows for one method-horizon is a duplicated evaluation, not a value.
#'
#' Shared by Figure 2 panel A (which annotates these numbers) and panel C (which cross-checks
#' its own ladder against them), so the two cannot read the table differently.
ev_get <- function(col, h, method = FEATURED, tbl = ev) {
  if (!col %in% names(tbl))
    stop("[fig2] invasion_evaluation.csv has no column `", col, "`; re-run the pipeline.",
         call. = FALSE)
  v <- tbl[[col]][tbl$method == method & tbl$horizon == h]
  if (length(v) != 1L)
    stop(sprintf("[fig2] %s at h=%s: expected one value for %s, got %d.",
                 method, format(h), col, length(v)), call. = FALSE)
  v
}

# FEATURED MODEL: READ, DO NOT RECOMPUTE. See the identical note in make_publication_figures.R.
# A hand-written mirror of best_invasion_model() lived here; it resolved the log-score axis
# before the eligibility gates (the bug fixed in the canonical function) and broke composite
# ties differently from the pipeline, and the top two Bayesian models tie at composite 39 on
# the shipped run. model_selection.json is the pipeline's own record of what it picked.
# read_featured_model() (figure_style.R) is the one reader every figure suite uses, and it has
# NO FALLBACK: this file used to fall back to a typed "Bayes-M14-fill-med" when the selection
# file was unreadable, which would publish a Figure 2 about a model the pipeline may not have
# chosen. A figure that cannot tell which model was chosen must fail, not guess.
FEATURED <- read_featured_model(file.path(OUT, "key_outputs", "model_selection.json"))
message(sprintf("[featured] %s", FEATURED))

prov_lookup <- rs %>% distinct(health_zone, province)

# "present" as-of date for the LIVE run (the final time point appended to the rank/front
# panels, downstream of the CV folds). Read from run_info.json and NEVER invented: it used to
# fall back to "last CV round + 28 d", a date no run had, drawn and labelled as the latest run.
# NA here is not fatal for the builders that do not draw it; build_rank_traj() refuses it.
PRESENT_DATE <- tryCatch({
  ri <- jsonlite::fromJSON(file.path(OUT, "key_outputs", "run_info.json"))
  d <- as.Date(ri$analysis_date); if (is.na(d)) stop("no analysis_date"); d
}, error = function(e) { message("[present] run_info analysis_date unavailable (", conditionMessage(e), ")"); as.Date(NA) })
message(sprintf("[present] live-run as-of date: %s", ifelse(is.na(PRESENT_DATE), "NA (unavailable)", format(PRESENT_DATE))))

# Last calendar day of data used for fitting (run_info training_window_end). Figure 2B's
# "invaded later / not invaded" outcome classes are judged against this date, so its legend
# must quote it. No fallback: a wrong date on a published legend is worse than a failed build.
DATA_END_DATE <- {
  ri <- jsonlite::fromJSON(file.path(OUT, "key_outputs", "run_info.json"))
  d <- if (is.null(ri$training_window_end)) as.Date(NA) else as.Date(ri$training_window_end)
  if (length(d) != 1 || is.na(d)) stop("run_info.json lacks training_window_end (last data day used for fitting)")
  d
}
DATA_END_LBL <- format(DATA_END_DATE, "%d %b")
message(sprintf("[data-end] last data day used for fitting: %s", format(DATA_END_DATE)))

# Province -> colour map shared by every zone panel coloured by province — Figure 3B, the
# Figure S4 pair, and panels C and D of Figure3_rank_surveillance (Figure 3A is a choropleth of
# relative risk and carries no province colour) — so every province keeps one hue. THE SUITE'S ONE PROVINCE SCALE, prov_palette() (figure_style.R),
# in its list(levels, colours) form. This file used to carry its own copy whose spare pool was
# Okabe-Ito minus PROV_COL — two hues — so a third province without an identity hue recycled
# the first spare and two provinces shared a colour on one panel, and a spare province could
# take a different hue here from the one every other figure gives it.
province_palette <- function(provinces) prov_palette_levels(provinces, prov_int = PROV_INT)

# -----------------------------------------------------------------------------
# 3. FIGURE 1 DATA  (onset-dated confirmed cases; pipeline-consistent)
# -----------------------------------------------------------------------------
# The line list is the PIPELINE's: load_linelist() resolves the latest processed
# export, restores each record's diagnostic sample date (SAMPLE_DATE_SOURCE; the earliest by
# default, apply_earliest_sample_dates()), tops each zone up to the INSP sitrep's
# confirmed count, and imputes a missing onset as the sample date minus an interval drawn from
# DONORS — records of the same class sampled in the same week (WEEK_ANCHOR grid) whose onset is
# recorded — at a per-record uniform from rlang::hash(alert_id, ONSET_IMPUTE_SEED)
# (.impute_onset_delays(), 01_data_prep.R). Figure 1 therefore shows exactly the case set, the
# zone set and the onset dates the forecasting pipeline scores.
#
# This script used to rebuild all of that itself and shift every missing onset back by a
# single ROUNDED MEAN delay. A fixed shift is not the same estimator: it collapses a
# right-skewed interval distribution onto its mean, so a zone whose first confirmed case has
# no recorded onset was dated to the day rather than sampled, and the resulting arrival dates
# matched no table in the pipeline.
# WHICH BUILDERS NEED THE LINE LIST. Figures 1, S1 and S3 are drawn from it directly, and
# so are Figure 3's SI companions (S4A's front-approach reads `zone_first`) — but those
# companions are not built in MI mode (see build_fig3). Main-text Figures 2 and 3 read only
# the published tables. Gating matters for more than speed: under MI the arrival-predictor
# tables below live in the canonical tree and NOT in the pooled mirror, so an ungated read
# would abort the whole suite at load time on a file the requested figures never use.
.need_linelist <- .want("Figure1") || .want("FigureS1") || .want("FigureS3") ||
                  (.want("Figure3") && !MI_MODE)
ll <- NULL; conf <- NULL; zone_cum <- NULL; zone_first <- NULL; daily <- NULL; zone_curve <- NULL
# The first day of the pipeline's week grid (WEEK_ANCHOR, 00_config.R: weeks end on the analysis
# weekday). It starts Figures 1B and S3's date axes, so their ticks fall on the grid's weekday
# (fs_date_breaks() steps in whole weeks from the axis start). It used to be a Monday floor —
# a second week grid that matched no table in the pipeline.
wk_floor <- floor_date(OUTBREAK_START, "week", week_start = WEEK_ANCHOR)
if (.need_linelist) {
  ll <- load_linelist()
  stopifnot(all(c("health_zone", "province", "date_index", "confirmed") %in% names(ll)))

  conf <- ll %>% filter(confirmed, !is.na(date_index), date_index >= wk_floor, !is.na(health_zone))

  zone_cum   <- conf %>% group_by(health_zone, province) %>% summarise(cases = n(), .groups = "drop")
  zone_first <- conf %>% group_by(health_zone) %>% summarise(first = min(date_index), .groups = "drop")
  daily      <- conf %>% count(date_index, name = "n") %>% arrange(date_index) %>% mutate(cum_cases = cumsum(n))
  zone_curve <- zone_first %>% arrange(first) %>% mutate(n_zones = row_number())
  message(sprintf("[fig1] %d confirmed cases, %d zones, onset %s to %s", nrow(conf), nrow(zone_cum), min(daily$date_index), max(daily$date_index)))
} else {
  message("[fig1] line list not loaded: no requested builder reads it.")
}

EPI_ORIGIN <- EPI_ZONES[1]

# -----------------------------------------------------------------------------
# FIGURE 1C / S1 DATA — arrival time vs predictors
# -----------------------------------------------------------------------------
# READ ONLY. The per-zone predictors AND the per-predictor R^2 / rho / n / fitted line
# are computed by 43_spread_kinematics.R (Section 5b) and published to key_outputs.
# This panel previously rebuilt all of it here — scraping node coordinates out of a raw
# HTML file, re-deriving the epicentre, the great-circle and road distances and the
# mobility shares, and imputing onsets with a single fixed mean delay instead of the
# pipeline's fitted draw. It therefore published correlations computed on a different
# zone set and different arrival dates from every table in the pipeline, reproducible
# from none of them. A figure must not fit a model.
FIG1C_SRC  <- file.path(OUT, "key_outputs", "arrival_predictors.csv")
FIG1C_FITS <- file.path(OUT, "key_outputs", "arrival_predictor_fits.csv")
# A MISSING table is fatal, not a skipped panel. Degrading quietly would publish a
# Figure 1 that silently lost Panel C, and the stale copy of an earlier run's panel would
# stay on disk beside it.
.read_fig1c <- function(f, what) {
  if (!file.exists(f))
    stop("[fig1c] ", what, " not found (", basename(f), "). ",
         "Figure 1C and S1 are drawn from tables that 43_spread_kinematics.R publishes; ",
         "run it before this script.", call. = FALSE)
  suppressWarnings(read_csv(f, show_col_types = FALSE))
}
# Read only for the builders that draw them (Figures 1C and S1) — same reasoning as the
# line-list gate above: a fatal read for a panel nobody asked for is still fatal.
fig1c_dat <- NULL; fig1c_fits <- NULL; ARRIVAL_REF <- as.Date(NA)
if (.want("Figure1") || .want("FigureS1")) {
  fig1c_dat  <- .read_fig1c(FIG1C_SRC,  "arrival-predictor table")
  fig1c_fits <- .read_fig1c(FIG1C_FITS, "arrival-predictor fits")
  stopifnot(nrow(fig1c_dat) > 0L, nrow(fig1c_fits) > 0L,
            "arrival_origin" %in% names(fig1c_dat))
  ARRIVAL_REF <- as.Date(fig1c_dat$arrival_origin[1])
  stopifnot(!is.na(ARRIVAL_REF),
            # one outbreak, one origin: a per-row origin would make arrival_days incomparable
            length(unique(fig1c_dat$arrival_origin)) == 1L)
  message(sprintf("[fig1c] %d zone(s) from %s, arrival origin %s",
                  nrow(fig1c_dat), basename(FIG1C_SRC), format(ARRIVAL_REF, "%Y-%m-%d")))
}

# =============================================================================
# FIGURE 1
# =============================================================================
# Panel-C predictor labels (shared by the main figure and the SI variant). These
# now sit on the y-axis strips because Panel C's axes are switched: arrival time
# on x, predictor on y.
#
# The two mobility predictors name the MONTH of the data each was computed from, never typed
# here. "May" and "March" used to be literals, correct only until the next short-trip snapshot
# arrived or the April relocation arm was switched on. An unknown month is omitted, not guessed.
#   * short-trip: the snapshot 43_spread_kinematics.R actually used, which it publishes per table
#     as `shorttrip_snapshot` (chosen by build_M1's rule). This used to be the NEWEST discovered
#     snapshot, which need not be the one the predictor was computed from.
#   * relocation: FLOWMINDER_OD_MONTH (the relocation file's month), as 43 reads it.
# The NAME comes from month.name, not format(, "%B"), which follows the session locale — the
# same rule 43 applies to its own labels.
.month_lab <- function(d) if (length(d) == 1L && !is.na(d))
  paste0(month.name[as.integer(format(as.Date(d), "%m"))], " ") else ""
.st_snap <- if (!is.null(fig1c_dat) && "shorttrip_snapshot" %in% names(fig1c_dat))
  as.Date(fig1c_dat$shorttrip_snapshot[1]) else as.Date(NA)
# THE EPICENTRE IS NAMED FROM EPI_ZONES (EPICENTRE_ZONES), the set 43_spread_kinematics.R
# measures every distance and share from, not typed: the labels used to read "M+B+R", which
# would have gone on naming three zones after the configured set changed. Spelled out, on a
# line of its own. The two shares are of the set's POOLED outflow ("from A + B + C"); the three
# distances are to the NEAREST zone of the set (43's minimum over EPICENTRE_ZONES), and say so.
EPI_LBL     <- paste(EPI_ZONES, collapse = " + ")
EPI_LBL_ANY <- paste(EPI_ZONES, collapse = ", ")
METRIC_LAB <- c(
  neglog_share_shorttrip  = sprintf("Mobility distance from\n%s\n(-log %sshort-trip outflow share)",
                                    EPI_LBL, .month_lab(.st_snap)),
  neglog_share_relocation = sprintf("Relocation distance from\n%s\n(-log %sFlowminder relocation share)",
                                    EPI_LBL,
                                    .month_lab(get0("FLOWMINDER_OD_MONTH", ifnotfound = as.Date(NA)))),
  dist_greatcircle_km     = sprintf("great-circle distance to the nearest of\n%s\n(km)", EPI_LBL_ANY),
  dist_road_km            = sprintf("road distance to the nearest of\n%s\n(km)", EPI_LBL_ANY),
  travel_time_h           = sprintf("road travel time to the nearest of\n%s\n(h)", EPI_LBL_ANY),
  log10_population        = "population (log10 people)")

# arrival-time (x) vs predictor (y), one facet per predictor, R^2 + Spearman rho.
# `metrics` selects which predictors to show (2 for the main figure, all 5 for SI).
# The predictor name is placed as the per-facet y-axis label (strip on the left), and
# the plot header is dropped, so each panel reads as a standard labelled scatter.
build_fig1c <- function(metrics, nrow = 1, base = 12.5, pt_size = 1.9, lab_size = 3.7) {
  miss <- setdiff(metrics, names(fig1c_dat))
  if (length(miss)) stop("[fig1c] predictor(s) absent from ", basename(FIG1C_SRC), ": ",
                         paste(miss, collapse = ", "), call. = FALSE)
  labs_sel <- METRIC_LAB[metrics]
  stopifnot(!anyNA(labs_sel))
  # THE POINTS ARE THE FIT'S ZONES. 43_spread_kinematics.R excludes from every published fit
  # the epicentre zones (distance 0, arrival ~day 0) and the right-truncated tail (zones whose
  # first onset falls after the last fit-window week, `in_fit_window` FALSE: that week is
  # incompletely observed); drawing either would put points under a line and an n that do not
  # include them. A table without the flags predates those rules.
  if (!"epicentre" %in% names(fig1c_dat))
    stop("[fig1c] ", basename(FIG1C_SRC), " has no `epicentre` column: it predates the ",
         "epicentre-exclusion rule. Re-run 43_spread_kinematics.R.", call. = FALSE)
  if (!"in_fit_window" %in% names(fig1c_dat))
    stop("[fig1c] ", basename(FIG1C_SRC), " has no `in_fit_window` column: it predates the ",
         "right-truncation rule. Re-run 43_spread_kinematics.R.", call. = FALSE)
  # The short-trip strip names the month of the snapshot the predictor used (METRIC_LAB); one
  # table, one snapshot.
  if ("neglog_share_shorttrip" %in% metrics &&
      (!"shorttrip_snapshot" %in% names(fig1c_dat) ||
       length(unique(fig1c_dat$shorttrip_snapshot)) != 1L))
    stop("[fig1c] ", basename(FIG1C_SRC), " does not carry one `shorttrip_snapshot` date, so ",
         "the short-trip label cannot name the month its predictor was computed from. Re-run ",
         "43_spread_kinematics.R.", call. = FALSE)
  long <- fig1c_dat %>%
    filter(!(epicentre %in% TRUE), in_fit_window %in% TRUE) %>%
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "y") %>%
    filter(is.finite(y), is.finite(arrival_days)) %>%
    mutate(metric = factor(metric, levels = metrics, labels = labs_sel))

  # n IS REPORTED PER FACET, from the published fit. The facets are NOT on common
  # support: a zone drops out of a mobility facet precisely when no flow out of the
  # epicentre was released for it, which is selection on the predictor and correlated
  # with arrival. Comparing R^2 across facets as though they shared a denominator is
  # the thing this panel invites a reader to do, so the denominator is on the panel.
  #
  # The numbers are QUOTED, not bare: inside parse = TRUE, R normalises a numeric
  # literal, so "%.2f" renders -0.70 as -0.7 and 0.0018 as an exact 0 — asserting a
  # precision, and in the population facet an exact zero, the data do not support.
  # TWO LINES, NOT ONE. On the five-facet SI version each facet is ~3.3 in wide and the
  # one-line form ("R² = 0.36, ρ = 0.60, n = 53") ran past the panel's right edge, so the
  # sample size — the very number this annotation exists to report — was clipped off every
  # facet. Stacking with atop() roughly halves the width; the y expansion below is widened to
  # make room for the second line.
  stat_df <- fig1c_fits %>%
    filter(predictor %in% metrics) %>%
    mutate(metric = factor(predictor, levels = metrics, labels = labs_sel)) %>%
    transmute(metric, intercept, slope,
              lab = sprintf("atop(R^2 == '%.2f' * ',' ~~ rho == '%.2f', n == %d ~ 'zones')",
                            r2, spearman, n))
  stopifnot(nrow(stat_df) == length(metrics))

  ggplot(long, aes(arrival_days, y)) +
    geom_point(colour = PT_BLUE, size = pt_size, alpha = 0.9) +
    # the PUBLISHED line (intercept, slope from arrival_predictor_fits.csv), not a refit
    geom_abline(data = stat_df, aes(intercept = intercept, slope = slope),
                colour = FIT_RED, linewidth = 0.9) +
    # ANCHORED TO THE PANEL, NOT TO THE DATA. -Inf/Inf place the label against the panel's
    # own top-left corner, so the annotation can never be pushed off-panel by an extreme
    # observation, and the expansions below are what guarantee the clearance. (It used to be
    # anchored at min(x), max(y) of the drawn points, which is a data value and therefore
    # moves with the data.)
    geom_text(data = stat_df, aes(x = -Inf, y = Inf, label = lab), parse = TRUE,
              hjust = -0.05, vjust = 1.15, size = lab_size, colour = INK) +
    facet_wrap(~metric, scales = "free_y", nrow = nrow, strip.position = "left") +
    # Right-hand x expansion and the larger TOP y expansion are the label's clearance: the
    # two-line annotation sits inside the top-left corner of each panel and must not overlap
    # either the points or the panel edge.
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.30))) +
    labs(x = sprintf("Arrival time (days from %s)", format(ARRIVAL_REF, "%Y-%m-%d")), y = NULL) +
    theme_pub(base) +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text = element_text(size = base - 0.4, colour = INK, face = "bold", lineheight = 0.95),
          panel.spacing.x = unit(12, "pt"), panel.spacing.y = unit(12, "pt"))
}

build_fig1 <- function() {
  # --- 1A: BIVARIATE choropleth — cumulative case burden (x) x date of first onset (y) ---
  # ARRIVAL IS THE FIRST-ONSET DATE ITSELF, for the y class and for rho below alike. Both used to
  # read a week number on a Monday grid that no table in the pipeline uses; binning only adds
  # ties, and which zones it ties depends on the weekday the bins start on (5 of 67 zones
  # changed class between the Monday and the WEEK_ANCHOR grid on the 2026-10-09 line list).
  af <- zone_cum %>% left_join(zone_first %>% dplyr::select(health_zone, first), by = "health_zone") %>%
    mutate(bx = qtile(cases), by = qtile(as.numeric(first)),
           bikey = ifelse(is.na(bx) | is.na(by), NA, paste0(bx, "-", by)))
  mp <- join_map(af)
  # non-invaded zones (no bivariate class) get a distinct "No confirmed cases" fill + legend swatch,
  # so they are visually separated from the light-grey low-cases/early-onset bivariate corner.
  mp <- mp %>% mutate(fillkey = ifelse(is.na(bikey), "No confirmed cases", bikey))
  inv <- mp %>% filter(!is.na(cases))
  pts <- suppressWarnings(st_point_on_surface(st_geometry(inv)))
  sym <- cbind(st_drop_geometry(inv)[, c("health_zone","cases")], st_coordinates(pts)) %>% rename(x = X, y = Y)
  # label high-burden zones (>=100 confirmed) plus a few named zones of interest regardless of burden
  ALWAYS_LABEL <- c("Nia Nia")
  lab <- sym %>% filter(cases >= 100 | key(health_zone) %in% key(ALWAYS_LABEL))

  # correlation between arrival time (first-onset date) and cumulative burden across invaded zones.
  # QUOTED, as in panel C: inside parse = TRUE a bare numeric literal is normalised, so "%.2f"
  # would print -0.70 as -0.7.
  cor_df <- af %>% filter(!is.na(cases), !is.na(first))
  rho_ac <- suppressWarnings(cor(as.numeric(cor_df$first), cor_df$cases, method = "spearman"))
  cor_lab <- sprintf("rho[arrival~vs~burden] == '%.2f' ~~ '(n =' ~ %d ~ 'zones)'", rho_ac, nrow(cor_df))

  # province outlines (dissolve health zones by province)
  prov_bounds <- shp %>% dplyr::filter(!is.na(.prov)) %>% dplyr::group_by(.prov) %>%
    dplyr::summarise(.groups = "drop") %>% sf::st_make_valid()

  # Kisangani (Tshopo) reference city — the downstream urban hub, marked whatever its case
  # status (its fill shows that like any other zone's)
  ki  <- shp[shp$.key == "makiso kisangani", ]
  kic <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(ki))))
  ki_df <- data.frame(x = kic[1, 1], y = kic[1, 2], lab = "Kisangani (Tshopo)")

  # extent: invaded cluster + epicentre + Kisangani, so the hub sits in-frame
  foc_keys <- unique(c(key(zone_cum$health_zone), key(EPI_ORIGIN), "makiso kisangani"))
  foc <- shp[shp$.key %in% foc_keys, ]; bb <- sf::st_bbox(foc)
  px <- 0.06 * as.numeric(bb["xmax"] - bb["xmin"]); py <- 0.06 * as.numeric(bb["ymax"] - bb["ymin"])
  ZX <- c(bb["xmin"] - px, bb["xmax"] + px); ZY <- c(bb["ymin"] - py, bb["ymax"] + py)

  p1a <- ggplot(mp) +
    geom_sf(aes(fill = fillkey), colour = "white", linewidth = 0.08) +
    geom_sf(data = prov_bounds, fill = NA, colour = "grey30", linewidth = 0.32) +
    geom_point(data = ki_df, aes(x, y), shape = 21, fill = ACCENT, colour = "white", size = 2.9, stroke = 0.6) +
    ggrepel::geom_text_repel(data = lab, aes(x, y, label = health_zone), size = 3.5, fontface = "bold",
      colour = INK, min.segment.length = 0, segment.colour = FAINT, segment.size = 0.3,
      box.padding = 0.55, point.padding = 0.35, max.overlaps = 30, seed = 1) +
    ggrepel::geom_text_repel(data = ki_df, aes(x, y, label = lab), size = 3.6, fontface = "bold.italic",
      colour = ACCENT, min.segment.length = 0, segment.colour = ACCENT, segment.size = 0.35,
      box.padding = 0.7, point.padding = 0.4, nudge_y = 0.35, seed = 1) +
    annotate("text", x = ZX[1] + 0.02 * diff(ZX), y = ZY[2] - 0.02 * diff(ZY),
             label = cor_lab, parse = TRUE, hjust = 0, vjust = 1, size = 3.8, colour = INK) +
    scale_fill_manual(values = c(BIV_ARRIVAL, "No confirmed cases" = NO_CASE_FILL),
                      breaks = "No confirmed cases", na.value = NO_CASE_FILL, name = NULL,
                      guide = guide_legend(override.aes = list(colour = "grey60"))) +
    coord_sf(xlim = ZX, ylim = ZY, expand = FALSE) +
    theme_map(13.5) + theme(panel.border = element_rect(fill = NA, colour = GRID, linewidth = 0.4),
                            legend.position = "right", legend.key.width = unit(13, "pt"),
                            legend.key.height = unit(13, "pt"))

  locator <- ggplot(shp) +
    geom_sf(fill = "grey86", colour = "white", linewidth = 0.04) +
    annotate("rect", xmin = ZX[1], xmax = ZX[2], ymin = ZY[1], ymax = ZY[2], fill = NA, colour = ACCENT, linewidth = 0.55) +
    coord_sf(expand = FALSE) + theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = FAINT, linewidth = 0.3), plot.margin = margin(2,2,2,2))
  # Two-line axis titles. On one line each was longer than the 4x4 key it labels, so the x title
  # ran left past the key and the y title ran down past it, and the two crossed in the key's
  # bottom-left corner ("Cumu[l]ative confirmed cases" printed through "Later first onset").
  biv_leg <- bivar_legend(BIV_ARRIVAL, "Cumulative\nconfirmed cases ->", "Later first\nonset ->",
                          base = 11)
  p1a <- p1a +
    patchwork::inset_element(locator, left = 0.00, bottom = 0.00, right = 0.26, top = 0.26, align_to = "panel", clip = FALSE) +
    patchwork::inset_element(biv_leg, left = 0.015, bottom = 0.60, right = 0.26, top = 0.88, align_to = "panel", clip = FALSE)

  # --- 1B: epidemic trajectory in ONE panel with two axes ---
  #   left axis  : cumulative confirmed cases, stacked by the epicentre zones (EPI_ZONES)
  #                + all others (share shifts to "other" as it spreads)
  #   right axis : number of health zones with >=1 confirmed case (step line)
  seed_zones <- zone_first %>% filter(health_zone %in% EPI_ZONES) %>% arrange(first, health_zone) %>% pull(health_zone)
  all_dates  <- sort(unique(conf$date_index))
  grp_levels <- c(seed_zones, "Other zones")
  grp_cum <- conf %>%
    mutate(grp = ifelse(health_zone %in% seed_zones, health_zone, "Other zones")) %>%
    count(date_index, grp, name = "n") %>%
    tidyr::complete(date_index = all_dates, grp = grp_levels, fill = list(n = 0)) %>%
    arrange(grp, date_index) %>% group_by(grp) %>% mutate(cum = cumsum(n)) %>% ungroup() %>%
    mutate(grp = factor(grp, levels = grp_levels))
  max_cases <- max(daily$cum_cases); max_zones <- max(zone_curve$n_zones); scl <- max_cases / max_zones
  zc2 <- zone_curve %>% dplyr::select(first, n_zones)
  if (min(zc2$first) > wk_floor)       zc2 <- bind_rows(tibble(first = wk_floor, n_zones = 0L), zc2)
  if (max(zc2$first) < max(all_dates)) zc2 <- bind_rows(zc2, tibble(first = max(all_dates), n_zones = max(zc2$n_zones)))

  # distinct hue shared by the "health zones invaded" step line AND the secondary axis
  # (line/ticks/text/title), so it is unambiguous which series maps to the right axis.
  ZONE_LINE <- "#762A83"
  # colourblind-safe hues (Okabe-Ito) for the M+B+R seed zones, grey for everything else
  seed_pal  <- c("#0072B2", "#E69F00", "#009E73")
  if (length(seed_zones) > length(seed_pal))
    stop("[fig1] ", length(seed_zones), " epicentre zones but ", length(seed_pal), " seed hues: a ",
         "zone would get an NA fill. Extend seed_pal.", call. = FALSE)
  fill_cols <- setNames(c(seed_pal[seq_along(seed_zones)], "#C6C6C6"), grp_levels)
  xlim <- c(wk_floor, max(all_dates))
  p1b <- ggplot(grp_cum, aes(date_index, cum)) +
    geom_area(aes(fill = grp), position = "stack", alpha = 0.95, colour = "white", linewidth = 0.12) +
    geom_step(data = zc2, aes(first, n_zones * scl, colour = "Cumulative health zones invaded"),
              linewidth = 1.1, direction = "hv", inherit.aes = FALSE) +
    scale_fill_manual(values = fill_cols, name = "Cumulative cases by zone",
                      guide = guide_legend(order = 1, nrow = 1, byrow = TRUE)) +
    scale_colour_manual(values = c("Cumulative health zones invaded" = ZONE_LINE), name = NULL,
                        guide = guide_legend(order = 2)) +
    scale_y_continuous(name = "Cumulative confirmed cases", labels = comma,
                       expand = expansion(mult = c(0, 0.06)),
                       sec.axis = sec_axis(~ . / scl, name = "Cumulative health zones with >=1 confirmed case",
                                           breaks = scales::pretty_breaks(5))) +
    # Breaks from fs_date_breaks(): a fixed two-week step put twelve labels on this ~6 in axis
    # and they overprinted. Seven is what the axis carries at this font size.
    scale_x_date(limits = xlim, breaks = fs_date_breaks(xlim, max_n = 7L), date_labels = "%d %b",
                 expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Symptom-onset date") +
    theme_pub(13) +
    theme(legend.position = "top", legend.box = "vertical", legend.spacing.y = unit(1, "pt"),
          axis.title.y.right = element_text(colour = ZONE_LINE, angle = 90, margin = margin(l = 5)),
          axis.text.y.right  = element_text(colour = ZONE_LINE),
          axis.ticks.y.right = element_line(colour = ZONE_LINE, linewidth = 0.4),
          axis.line.y.right  = element_line(colour = ZONE_LINE, linewidth = 0.5),
          plot.margin = margin(6, 10, 6, 6))

  # --- 1C: arrival time (x) vs road travel time + short-trip epicentre outflow (main text) ---
  p1c <- build_fig1c(c("travel_time_h", "neglog_share_shorttrip"), nrow = 1, base = 12.5)

  save_dual(p1a, "F1A_bivariate_map", 7.0, 5.8)
  save_dual(p1b, "F1B_trajectory",    6.6, 5.0)
  save_dual(p1c, "F1C_arrival_vs_predictors", 7.0, 3.6)

  fig1 <- ((wrap_elements(full = p1a) | wrap_elements(full = p1b)) / wrap_elements(full = p1c)) +
    plot_layout(heights = c(1.5, 1)) + plot_annotation(tag_levels = "A") &
    # ONE TAG STYLE (theme_pub()'s): every panel is wrapped, so the tags would otherwise take
    # patchwork's small plain default while the other figures' are bold.
    theme(plot.tag = element_text(size = 13 + 4.5, face = "bold", colour = INK))
  save_dual(fig1, "Figure1", 15.0, 12.0, dir = FIG_DIR)
  invisible(fig1)
}

# --- SI: arrival time vs ALL five predictors (same switched-axis layout as Panel C) ---
build_fig1_si <- function() {
  p <- build_fig1c(c("neglog_share_shorttrip", "neglog_share_relocation", "dist_greatcircle_km", "travel_time_h", "log10_population"), nrow = 1, base = 16, lab_size = 4.7)
  save_dual(p, "FigureS1_arrival_vs_predictors", 16.5, 6.0, dir = FIG_DIR)
  invisible(p)
}

# --- SI: 5-day rolling share of confirmed cases in the epicentre (EPI_ZONES) --------------
# For each onset day, share = (cases in the EPI_ZONES epicentre) / (all confirmed cases),
# smoothed as a centred 5-day rolling proportion (ratio of 5-day rolling sums, so days with
# few/zero cases don't blow up a naive daily ratio). Faint points = raw daily share.
build_fig_si_epi_share <- function(base = 16.5) {
  epi_k <- key(EPI_ZONES)
  agg <- conf %>% mutate(is_epi = key(health_zone) %in% epi_k) %>%
    group_by(date_index) %>% summarise(n_epi = sum(is_epi), n_tot = dplyr::n(), .groups = "drop")
  grid <- tibble::tibble(date_index = seq(min(conf$date_index), max(conf$date_index), by = "day"))
  d <- grid %>% left_join(agg, by = "date_index") %>%
    mutate(n_epi = dplyr::coalesce(n_epi, 0L), n_tot = dplyr::coalesce(n_tot, 0L)) %>% arrange(date_index)
  roll_sum <- function(x, k = 5L) { h <- k %/% 2L; n <- length(x)
    vapply(seq_len(n), function(i) sum(x[max(1L, i - h):min(n, i + h)]), numeric(1)) }
  d <- d %>% mutate(re = roll_sum(n_epi), rt = roll_sum(n_tot),
                    prop5 = ifelse(rt > 0, re / rt, NA_real_),
                    prop_daily = ifelse(n_tot > 0, n_epi / n_tot, NA_real_))
  cap <- "Line = 5-day rolling share  ·  points = daily share"
  # The same axis start as Figure 1B (the week grid's first day), and breaks from
  # fs_date_breaks(): a fixed two-week step overprinted its labels on this 9 in axis.
  xlim <- c(wk_floor, max(d$date_index))
  p <- ggplot(d, aes(date_index)) +
    geom_point(aes(y = prop_daily), colour = FAINT, size = 1.1, alpha = 0.75, na.rm = TRUE) +
    geom_line(aes(y = prop5), colour = PT_BLUE, linewidth = 1.2, na.rm = TRUE) +
    annotate("text", x = min(d$date_index), y = 0.02, label = cap, hjust = 0, vjust = 0, size = 4.3, colour = MUTED) +
    scale_y_continuous(labels = percent_format(1), limits = c(0, 1), expand = expansion(mult = c(0.02, 0.04))) +
    scale_x_date(limits = xlim, breaks = fs_date_breaks(xlim, max_n = 7L), date_labels = "%d %b",
                 expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Symptom-onset date",
         y = sprintf("Share of confirmed cases in the\n%s epicentre", EPI_LBL)) +
    theme_pub(base) + theme(axis.title.y = element_text(lineheight = 0.95))
  save_dual(p, "FigureS3_epicentre_case_share", 9.0, 5.6, dir = FIG_DIR)
  invisible(p)
}

# =============================================================================
# FIGURE 2  (key_outputs evaluation figure WITHOUT the across-model-grid panels)
# =============================================================================

# THE COMPARATORS (FIGURE_BASELINE_METHODS, 00_config.R). They span (uses case data) x (uses
# phone data): two case-blind epicentre rankings — measured Flowminder cohort inflow, and road
# travel time — which read no incidence, and two case-using structural rules — Distance-B1
# (neighbours' confirmed cases weighted by 1 / travel time) and Adjacency-B7 (nearest
# already-affected zone) — which read the case data to locate the front. None touches the
# renewal machinery. Published by run_all.R into key_outputs/detection_curves.csv.
# MEMBERSHIP, LABELS AND COLOURS ALL COME FROM figure_style.R. The set is
# FIGURE_BASELINE_METHODS (00_config.R, read through fs_baseline_methods()) so this suite and
# the publication suite cannot drift apart; BASELINE_LBL / BASELINE_COL / FEATURED_COL and
# fs_panelA_scales() are in the house-style file so that the CASCADE analogue of this panel
# (make_manuscript_figure2_cascade.R) draws the identical comparators in the identical colours
# under the identical names. This file used to define all four, and the cascade figure defined
# its own three — which is how the two published panels came to disagree about both which
# comparators exist and what #009E73 means.
BASELINE_METHODS <- fs_baseline_methods()

# READ ONLY. key_outputs/detection_curves.csv is published by run_all.R from
# compute_detection_curve() (20_forecast_detail.R) — the same estimator behind the
# recall_at_<K> columns of invasion_evaluation.csv, verified to reproduce them exactly.
#
# This panel used to build the curve itself. Three things were wrong with that:
#   * it re-derived the naive epicentre-inflow score from the raw M8 kernel, WITHOUT the
#     alias harmonisation naive_epicentre_inflow_scores() applies, so a non-canonical zone
#     spelling silently changed the baseline here but not in the pipeline — and it read a
#     hard-coded "mobility_M8-fill.rds" instead of MOBILITY_PRIMARY;
#   * it POOLED recall over folds (sum(hits)/sum(invasions)) while the published tables
#     AVERAGE it per fold, so the "share of true invasions caught" on this panel was a
#     different estimand from the recall_at_K beside it;
#   * it drew a 95% WILSON band, a binomial interval that treats hundreds of zone-weeks
#     within a fold as independent. The published interval is a fold-cluster bootstrap.
DETCURVE_SRC <- file.path(OUT, "key_outputs", "detection_curves.csv")
detection_curves <- {
  if (!file.exists(DETCURVE_SRC))
    stop("[fig2] detection_curves.csv not found. Figure 2A is drawn from the curve run_all.R ",
         "publishes; run the pipeline before this script.", call. = FALSE)
  d <- suppressWarnings(read_csv(DETCURVE_SRC, show_col_types = FALSE))
  stopifnot(nrow(d) > 0,
            all(c("method", "horizon", "k", "recall", "recall_lo", "recall_hi",
                  "recall_random") %in% names(d)))
  d
}
message(sprintf("[fig2] detection curves: %d method(s), horizons %s",
                dplyr::n_distinct(detection_curves$method),
                paste(sort(unique(detection_curves$horizon)), collapse = "/")))

prioritisation_curves <- function(h = 2L) {
  cu <- detection_curves %>% filter(horizon == h)
  if (!nrow(cu)) return(NULL)
  want <- c(FEATURED, BASELINE_METHODS)
  miss <- setdiff(want, unique(cu$method))
  # A MISSING BASELINE IS FATAL, not a quietly thinner panel. Panel 2A's entire claim is set
  # against the comparators FIGURE_BASELINE_METHODS names; dropping one on a message would
  # publish a figure that looks complete and silently makes a weaker claim, and the stale copy
  # of the full panel would sit on disk beside it. This fires if run_all.R has not been re-run
  # since the baseline set changed. The message names the set from the config, not a typed
  # list: it used to promise "the three structural baselines (gravity, Flowminder inflow,
  # travel time)" after gravity had left the set and two case-using rules had joined it.
  if (length(miss))
    stop("[fig2] Panel 2A needs ", length(want), " series but ", length(miss),
         " are absent from ", basename(DETCURVE_SRC), " at h=", h, ": ",
         paste(miss, collapse = ", "),
         ". Re-run run_all.R — it publishes the detection curves for the featured model and ",
         "the ", length(BASELINE_METHODS), " comparators in FIGURE_BASELINE_METHODS (",
         paste(BASELINE_METHODS, collapse = ", "), ").",
         call. = FALSE)
  cu <- cu %>% filter(method %in% want)
  if (!nrow(cu)) return(NULL)
  # recall_random is a property of the fold design, not of a method, so it is identical
  # across rows; assert that rather than silently taking one method's copy.
  rnd <- cu %>% distinct(k, recall_random) %>% arrange(k)
  if (nrow(rnd) != dplyr::n_distinct(cu$k)) {
    # This used to be a bare stopifnot, and when it fired -- because the methods were scored
    # on different (fold x zone) supports, so k/n_atrisk differed between families -- it said
    # only "nrow(rnd) == n_distinct(cu$k) is not TRUE", which is a long way from the cause.
    # Name the disagreement instead: which k, which methods, and the spread in n_atrisk that
    # produced it. What keeps them equal is run_all.R resolving the shared support ONCE per
    # horizon from the scored table and passing it as compute_detection_curve(support_cells=).
    bad <- rnd %>% count(k, name = "n_values") %>% filter(n_values > 1L)
    det <- cu %>% filter(k %in% bad$k) %>%
      distinct(method, k, recall_random) %>% arrange(k, method)
    stop("[fig2] recall_random is a property of the fold design and must be identical across ",
         "methods, but ", nrow(bad), " of ", dplyr::n_distinct(cu$k), " k values at h=", h,
         " carry more than one value. The methods were scored on different at-risk row sets, ",
         "so the random-targeting reference k/n_atrisk differs between them and Panel 2A ",
         "would draw one method's null beside another's curve. First disagreement:\n",
         paste(utils::capture.output(print(utils::head(as.data.frame(det), 8))),
               collapse = "\n"),
         "\nCheck that run_all.R passed support_cells = invasion_common_cells(lfo_results, h) ",
         "to compute_detection_curve(), and that detection_curves.csv is from the current run.",
         call. = FALSE)
  }
  list(curves = cu, rnd = rnd)
}

#' The featured model's scored rows at horizon `h`: the rows invasion_evaluation.csv scored it on.
#'
#' At-risk, finite, and on the SHARED (fold x zone) support, resolved by invasion_common_cells()
#' — the function evaluate_invasion() and compute_detection_curve() resolve it through. Every
#' other number in panel 2A's block comes from that table, so a Brier score taken on the
#' featured model's own rows would describe a different row set whenever any method covers
#' fewer cells. The support is resolved on the RAW probabilities (`p_invasion_raw`), as the
#' evaluation resolved it: on a recalibrated pass `p_invasion` holds p_recal.
#'
#' MI-POOLED FRAME (the row-wise stack of M arms, 48_onset_mi_pool.R): each arm was evaluated on
#' ITS OWN shared support — the at-risk set moves with the imputed labels — so the support is
#' resolved arm by arm, and `mi_arm` is kept for the per-arm metrics below.
featured_scored_rows <- function(h) {
  .one <- function(fr) {
    cells <- invasion_common_cells(dplyr::mutate(fr, p_invasion = p_invasion_raw), h)
    d <- fr %>% filter(method == FEATURED, horizon == h, is.finite(p_invasion),
                       is.finite(is_new_invasion), !(was_active_before %in% TRUE))
    # No resolvable support: every row, as evaluate_invasion() falls back.
    if (length(cells)) d <- d[paste(d$fold_id, d$health_zone, sep = "\r") %in% cells, , drop = FALSE]
    d
  }
  if (!"mi_arm" %in% names(lfo)) return(.one(lfo))
  dplyr::bind_rows(lapply(sort(unique(lfo$mi_arm)), function(a)
    .one(lfo[lfo$mi_arm %in% a, , drop = FALSE])))
}

#' Figure 2A's two scale-dependent numbers — the Brier score and calibration-in-the-large —
#' computed on rounds that share ONE probability scale.
#'
#' ON THE RECALIBRATED PASS, THE RECALIBRATED ROUNDS ONLY. A round with no fitted prequential
#' factor (delta_estimable FALSE; in the pipeline, the rounds issued before any settled,
#' post-burn-in round existed to train it on) carries its RAW probability (delta = 1). The block used to pool those rounds with the
#' recalibrated ones, so it quoted a calibration of neither scale, while panel B and FigureS14
#' label the same rounds "(uncorrected)". The header now says how many rounds the two numbers
#' cover (fig2_metric_header()). On the raw pass every round is on the one (raw) scale and all
#' are used.
#'
#' Calibration-in-the-large is sum(p) / sum(y): the predicted number of positive zone-rounds over
#' the observed number (evaluate_invasion()'s mean(p) / base rate). `calib_all` is the same
#' estimator over EVERY round, which must reproduce the evaluation table's EV_CAL column; the
#' caller checks that, so the block cannot be computed on rows other than the table's.
#'
#' MI: each arm is scored on its own rows and the arms are averaged with equal weight — the
#' Rubin point estimate the pooled evaluation table carries — so an arm is not weighted by how
#' many rows it happened to score.
#'
#' @param d featured_scored_rows(h).
#' @return list(brier, calib, calib_all, n_all, n_sel, recal); n_all / n_sel are the range over
#'   arms (one value twice on a single draw) of the rounds scored and the rounds used.
fig2_calibration_metrics <- function(d) {
  recal <- identical(PSCALE_COL, "p_recal")
  if (!nrow(d))
    stop("[fig2] the featured model ", FEATURED, " has no scored row at this horizon.", call. = FALSE)
  if (recal && !"delta_estimable" %in% names(d))
    stop("[fig2] the LFO frame carries p_recal but no `delta_estimable`, so the recalibrated ",
         "rounds cannot be told from the uncorrected ones. Re-run 16b_invasion_recalibration.R.",
         call. = FALSE)
  arms <- if ("mi_arm" %in% names(d)) split(d, d$mi_arm) else list(d)
  pa <- dplyr::bind_rows(lapply(arms, function(x) {
    xs <- if (recal) x[x$delta_estimable %in% TRUE, , drop = FALSE] else x
    tibble::tibble(n_all = dplyr::n_distinct(x$fold_id), n_sel = dplyr::n_distinct(xs$fold_id),
                   n_pos = sum(xs$is_new_invasion),
                   brier = mean((xs$p_invasion - xs$is_new_invasion)^2),
                   calib = sum(xs$p_invasion) / sum(xs$is_new_invasion),
                   calib_all = mean(x$p_invasion) / max(mean(x$is_new_invasion), 1e-9))
  }))
  if (any(pa$n_sel == 0L))
    stop("[fig2] no round of ", FEATURED, " carries a fitted recalibration factor, so the panel ",
         "would quote the Brier score and calibration of uncorrected probabilities on a ",
         "recalibrated figure.", call. = FALSE)
  if (any(pa$n_pos == 0))
    stop("[fig2] the rounds used for the Brier score and calibration hold no first case; ",
         "calibration-in-the-large (predicted / observed) is undefined.", call. = FALSE)
  list(brier = mean(pa$brier), calib = mean(pa$calib), calib_all = mean(pa$calib_all),
       n_all = range(pa$n_all), n_sel = range(pa$n_sel), recal = recal)
}

#' The first lines of panel 2A's metric block (fs_metric_block()'s header): the horizon and the
#' rounds the discrimination lines cover (ALL rounds: `n_skill_rounds` is the evaluation's
#' n_folds and must equal the rounds of the rows the block is computed from), then which rounds —
#' and which scale — the Brier score and calibration lines cover.
fig2_metric_header <- function(m, n_skill_rounds, what = "2-week-ahead forecast") {
  .n <- function(r) if (r[1] == r[2]) as.character(r[1]) else paste(r, collapse = "-")
  if (!all(as.integer(m$n_all) == as.integer(n_skill_rounds)))
    stop(sprintf(paste0("[fig2] the discrimination lines are scored on %d rounds (n_folds) but the ",
                        "block's rows span %s; the header would name the wrong rounds."),
                 as.integer(n_skill_rounds), .n(m$n_all)), call. = FALSE)
  sprintf("%s, all %d rounds\nBrier, calibration: %s",
          what, as.integer(n_skill_rounds),
          if (m$recal) paste(.n(m$n_sel), "recalibrated rounds") else "raw probabilities")
}

# -----------------------------------------------------------------------------
# Shared rank trajectory (2-week horizon), used by BOTH Figure-3 narrative panels
# so they stay a matched pair. The tracked ("keep") zones are the top-N at-risk
# zones AS OF THE PRESENT LIVE RUN (rs) — the exact set shown in Figure 3B — and
# the x-axis runs across the retrospective CV rounds (lfo) PLUS the present run,
# so the ranks are shown evolving all the way up to the latest forecast.
#   • CV rounds: rank within each fold by descending p_invasion  (lfo, horizon 2)
#   • present  : rank by descending p_case_invasion              (rs,  horizon 2)
# Both quantities are the same model's P(first case) within the 2-week window, so
# the per-timepoint orderings are comparable.
FRONT_H  <- 2L      # forecast horizon for the operational panels (2 weeks ahead)
FRONT_TOP_N <- 20L  # tracked at-risk zones (matches Figure 3B)
build_rank_traj <- function(N_TOP = FRONT_TOP_N, h = FRONT_H) {
  rk_h <- lfo %>% filter(method == FEATURED, horizon == h, is.finite(p_invasion), !(was_active_before %in% TRUE)) %>%
    mutate(cutoff = as.Date(cutoff)) %>% group_by(cutoff) %>%
    mutate(rank = rank(-p_invasion, ties.method = "min")) %>% ungroup() %>%
    transmute(health_zone, cutoff = lfo_origin(cutoff), rank)
  cutoffs <- sort(unique(rk_h$cutoff))
  if (is.na(PRESENT_DATE) || PRESENT_DATE <= max(cutoffs))
    stop(sprintf(paste0("[rank-traj] the live-run date (run_info.json analysis_date) is %s; it must ",
                        "fall after the last forecast round (%s). No date is substituted for it."),
                 if (is.na(PRESENT_DATE)) "unavailable" else format(PRESENT_DATE),
                 format(max(cutoffs))), call. = FALSE)
  present_date <- PRESENT_DATE
  now <- rs %>% filter(horizon == h, !(was_active_before %in% TRUE), is.finite(p_case_invasion)) %>%
    distinct(health_zone, .keep_all = TRUE) %>%
    mutate(rank = rank(-p_case_invasion, ties.method = "min"))
  keep_z <- now %>% arrange(desc(p_case_invasion)) %>% head(N_TOP) %>% pull(health_zone) %>% unique()
  rk_now <- now %>% transmute(health_zone, cutoff = present_date, rank)
  traj <- bind_rows(rk_h, rk_now) %>% filter(health_zone %in% keep_z)
  list(keep_z = keep_z, traj = traj, cutoffs = cutoffs, present_date = present_date)
}

#' Tick positions for a round axis: at most `max_n` labels, with the FIRST round and the
#' LIVE-RUN date always kept.
#'
#' The two Figure-S4 panels used to label every round. That was legible at 8 rounds; the
#' cross-validation now admits every round whose outcome window closed at least
#' LFO_MIN_EVAL_AGE_DAYS before the analysis date, and that many dated labels, plus the live
#' run's, overprint on a 7-inch panel. The endpoints are pinned because they are the two the
#' reader is asked to compare: where the front started and where it is now.
.round_breaks <- function(times, present_date, max_n = 6L) {
  ts <- sort(unique(as.Date(times)))
  if (length(ts) <= max_n) return(ts)
  idx <- unique(c(1L, round(seq(1, length(ts), length.out = max_n)), length(ts)))
  sort(unique(c(ts[idx], as.Date(present_date))))
}

# invasion-risk rank evolution across forecast rounds up to the present live run
# (2-week horizon; by province). Extracted so it can live in Figure 3.
build_rank_evolution <- function(base = 12, pal = NULL) {
  RANK_FLOOR <- 30L
  tj <- build_rank_traj()
  present_date <- tj$present_date
  if (is.null(pal)) pal <- province_palette(prov_lookup$province[prov_lookup$health_zone %in% tj$keep_z])
  traj <- tj$traj %>% left_join(prov_lookup, by = "health_zone") %>%
    mutate(province = factor(province, levels = pal$levels), rank_disp = pmin(rank, RANK_FLOOR))
  end_lab <- traj %>% filter(cutoff == present_date); x_pad <- as.numeric(diff(range(traj$cutoff))) * 0.06 + 1
  x_breaks <- .round_breaks(traj$cutoff, present_date); x_labs <- format(x_breaks, "%d %b")
  ggplot(traj, aes(cutoff, rank_disp, group = health_zone, colour = province)) +
    annotate("rect", xmin = present_date - 2, xmax = present_date + 2, ymin = -Inf, ymax = Inf, fill = "grey95") +
    geom_vline(xintercept = present_date, linetype = "22", colour = MUTED, linewidth = 0.45) +
    geom_line(linewidth = 1.0, alpha = 0.55, lineend = "round") +
    geom_point(colour = "white", size = 2.4) + geom_point(size = 1.45) +
    annotate("text", x = present_date, y = 0.4, label = "latest run", hjust = 0.5, vjust = 0, size = 2.5, colour = MUTED, fontface = "italic") +
    ggrepel::geom_text_repel(data = end_lab, aes(label = health_zone), size = 2.7, direction = "y", hjust = 0,
      nudge_x = x_pad, box.padding = 0.1, segment.size = 0.2, segment.colour = FAINT, min.segment.length = 0, max.overlaps = 40, seed = 1, show.legend = FALSE) +
    scale_y_reverse(breaks = c(1,5,10,15,20,25,30), labels = c("1","5","10","15","20","25","30+"), expand = expansion(mult = c(0.06, 0.04))) +
    scale_x_date(breaks = x_breaks, labels = x_labs, expand = expansion(mult = c(0.03, 0.22))) +
    scale_colour_manual(values = pal$colours, name = "Province", breaks = pal$levels, drop = FALSE, na.value = FAINT) +
    labs(x = "Forecast round  ·  shaded = latest (live) run", y = "Invasion-risk rank (1 = highest)") +
    theme_pub(base) + theme(legend.position = "right", panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
                        axis.text.x = element_text(angle = 30, hjust = 1, size = base - 4))
}

# -----------------------------------------------------------------------------
# DOES A PANEL'S LEGEND FIT INSIDE THE PANEL?
# -----------------------------------------------------------------------------
# .assert_legend_fits() and its measuring helpers MOVED TO figure_style.R (2026-09-28), with
# the full account of the defect they exist to catch. They were defined here, so the cascade
# analogue of Figure 2 — a composed figure with a legend over each panel, which had the same
# overflow on the page — could not call them. Nothing about the check is Figure-2 specific;
# `who =` names the caller in its messages.

# -----------------------------------------------------------------------------
# FIGURE 2, PANEL C — what the mobile-phone mobility data adds
# -----------------------------------------------------------------------------
# WHAT IT IS. The ladder 50_mobility_value.R computes: every method's AUC-PR skill on the
# shared (origin x zone) support, with a 90% PAIRED zone-cluster bootstrap interval, ordered
# by skill, coloured by what its mobility term is built from and shaped by what kind of method
# it is. The dashed line is the best method on the WHOLE ladder that uses no mobile-network
# data at any point (mobval_ladder_frame() guarantees it is drawn); the distance from it to the
# featured model is what the phone data buys.
#
# DRAWN BY THE SAME BUILDER AS FigureS24, from the same CSV, via mobility_value_panel.R. The
# main-text call differs in exactly three ways, all of them layout:
#   * ONE HORIZON. Panels A and B are fixed to 2 weeks ahead, and so is this. Drawing the
#     1-week facet here as well would put a second "AUC-PR skill" for the featured model in
#     the same figure as panel A's annotated 2-week number - two different values of one named
#     quantity, one figure. FigureS24 keeps both horizons side by side, which is where that
#     comparison belongs.
#   * NO FACET STRIP. With one horizon the strip is a lone centred bold heading across the
#     whole panel, which reads as the plot title theme_pub() deliberately blanks. The horizon
#     moves to the axis title instead.
#   * THE FIGURE-2 BASE FONT, with every hard-coded size in the panel scaled with it.
#
# COLOUR MEANS SOMETHING DIFFERENT HERE THAN IN PANEL A, as it already does in panel B: A
# colours METHODS (five named series), B colours OUTCOMES, C colours a CLASS of mobility term.
# Each panel carries its own legend directly above it. The alternative - recolouring C to
# panel A's method colours - cannot be done: A distinguishes five methods that C groups into
# three classes, and it would also break C away from FigureS24, which is the one thing this
# panel must not do. What IS required, and checked in mobility_value_panel.R, is that no hue
# carries one meaning in A and a different one in C (C's "Geography only" used to be A's
# featured-model blue), and that a method both panels draw carries ONE name in both
# (FEATURED_LBL / BASELINE_LBL).
MOBVAL_LADDER_SRC <- file.path(OUT, "key_outputs", "mobility_value_ladder.csv")

#' Load the mobility-value ladder and run every check that keeps panel C honest beside
#' panels A and B. The two panel-C drawers below differ in what they DRAW from it, never in
#' what is VERIFIED about it.
#'
#' @param h forecast horizon, in weeks. 2 everywhere in Figure 2.
#' @return the validated ladder frame (all horizons; the checks are run at `h`).
fig2_ladder_checked <- function(h = 2L) {
  if (!file.exists(MOBVAL_LADDER_SRC))
    stop("[fig2] Figure 2 panel C is drawn from ", basename(MOBVAL_LADDER_SRC),
         ", which is not in ", file.path(OUT, "key_outputs"),
         ". 50_mobility_value.R publishes it (run_all.R runs that module in the downstream ",
         "block, before this suite); run it, or run the pipeline, before rendering Figure 2. ",
         "Panel C is NOT optional: a Figure 2 that silently drops it would sit on disk under ",
         "the same basename as the three-panel figure it replaced.", call. = FALSE)
  lad <- suppressWarnings(read_csv(MOBVAL_LADDER_SRC, show_col_types = FALSE))
  need <- c("method", "horizon", "skill", "lo", "hi", "base_rate", "n_cells",
            "mv_label", "mv_class", "mv_kind", "phone")
  miss <- setdiff(need, names(lad))
  if (length(miss))
    stop("[fig2] ", basename(MOBVAL_LADDER_SRC), " is missing column(s): ",
         paste(miss, collapse = ", "), ". Re-run 50_mobility_value.R.", call. = FALSE)
  lad <- lad %>% mutate(horizon = as.integer(horizon))
  if (!h %in% lad$horizon)
    stop("[fig2] the mobility-value ladder carries horizons ",
         paste(sort(unique(lad$horizon)), collapse = "/"), " but panel C needs h=", h, ".",
         call. = FALSE)

  # (1) THE SAME FEATURED MODEL AS PANELS A AND B, BY CONSTRUCTION. Panel C is handed FEATURED
  # (read from model_selection.json) for its "(featured)" tag and its like-for-like row set;
  # there is no second, typed name to diverge from it. That typed name (MV_FEATURED) is what
  # stopped Figure 2 building when the selection moved from -geo to -med. Check (2) below then
  # ties the ladder's number for that model to panel A's.

  # (2) THE LADDER IS FROM THIS RUN. Panel A annotates the featured model's AUC-PR skill from
  # diagnostics/invasion_evaluation.csv; panel C draws the same quantity, for the same model,
  # at the same horizon, from a CSV that 50_mobility_value.R writes in a separate pass. A
  # ladder left over from an earlier run would put two different numbers for one quantity in
  # one figure, with nothing on the page to show which is current. 50 pins its own average
  # precision to the published table at 1e-9, so anything above a relative 1e-6 here is a
  # stale file, not arithmetic.
  ladf <- lad %>% filter(method == FEATURED, horizon == h)
  if (nrow(ladf) != 1L)
    stop("[fig2] the mobility-value ladder holds ", nrow(ladf), " rows for ", FEATURED,
         " at h=", h, "; expected exactly one.", call. = FALSE)
  .rel <- function(a, b) abs(a - b) / max(abs(b), .Machine$double.eps)
  ev_skill <- ev_get(EV_SKILL, h)
  if (.rel(ladf$skill, ev_skill) > 1e-6)
    stop(sprintf(paste0("[fig2] panel C would contradict panel A: %s at h=%d has AUC-PR skill ",
                        "%.6f in %s but %.6f in invasion_evaluation.csv. The ladder is from a ",
                        "different run - re-run 50_mobility_value.R against this run's ",
                        "lfo_results.rds."),
                 FEATURED, h, ladf$skill, basename(MOBVAL_LADDER_SRC), ev_skill),
         call. = FALSE)
  # The same run, scored on the same rows. The ladder scores the whole shared support (every
  # round), which it rebuilds independently: its cell count must equal the featured model's scored
  # rows (which ARE the shared support) and the evaluation's n_atrisk, and its positive count the
  # evaluation's n_invasions, so the two panels describe one cross-validation.
  .fs_all <- featured_scored_rows(h)
  if (!identical(as.integer(ladf$n_cells), nrow(.fs_all)) ||
      !identical(as.integer(ladf$n_cells), as.integer(ev_get("n_atrisk", h))))
    stop("[fig2] the ladder scored ", ladf$n_cells, " zone-rounds at h=", h,
         " but the featured model's scored rows number ", nrow(.fs_all), " (n_atrisk ",
         ev_get("n_atrisk", h), "); panels A and C are not on the same support.", call. = FALSE)
  ev_e <- ev_get("n_invasions", h)
  if (!is.finite(ev_e) || !identical(as.integer(ladf$n_events), as.integer(ev_e)))
    stop("[fig2] the ladder counts ", ladf$n_events, " positive zone-rounds at h=", h,
         " but invasion_evaluation.csv's n_invasions is ", ev_e, ".", call. = FALSE)

  # "positive zone-rounds", not "events": at h=2 one first case is scored in the two
  # overlapping windows that contain it (n_invasion_events in invasion_evaluation.csv counts
  # the distinct ones).
  message(sprintf("[fig2] ladder: h=%d, %s zone-rounds, %s positive zone-rounds, %s bootstrap resamples",
                  h, ladf$n_cells,
                  if ("n_events" %in% names(ladf)) ladf$n_events else "?",
                  if ("n_boot" %in% names(ladf)) format(ladf$n_boot, big.mark = ",") else "?"))
  lad
}

#' Panel C as published: the full ladder (mv_ladder_rows()).
#'
#' @param base theme_pub() base size; the panel's own sizes scale with it.
#' @param h    forecast horizon, in weeks. 2 everywhere in Figure 2.
#' @param lad  the checked ladder frame; resolved once by the caller so the published and
#'   simplified panels cannot be built from different files or different runs.
build_fig2c <- function(base, h = 2L, lad = fig2_ladder_checked(h)) {
  A <- mobval_ladder_frame(lad, "skill", featured = FEATURED)
  p <- mobval_panel(A, "skill", horizons = h, show_y = TRUE, base = base, strip = FALSE,
                    # The full width of Figure 2, not a ~5-inch appendix facet: four ticks
                    # across ten inches leaves the reader interpolating.
                    n_breaks = 7)
  message(sprintf("[fig2] 2C: %d methods at h=%d",
                  dplyr::n_distinct(A$method[A$horizon == h]), h))
  # How many printed label LINES the rows take (a wrapped label takes two), for the caller to
  # size the panel: a fixed height drew the rows of a longer ladder into one another.
  attr(p, "n_label_lines") <- sum(lengths(strsplit(levels(droplevels(A$mv_label[A$horizon == h])), "\n")))
  p
}

#' Panel C for the SIMPLIFIED figure: the same ladder reduced to its four load-bearing rows.
#'
#' Same CSV, same intervals, same shared support, same validation (`lad` is the checked frame
#' the published panel draws from). mobval_summary_frame() owns the selection rules and the
#' checks that keep the four rows disjoint, comparable and consistent with the appendix
#' ladder; see the block above MV_SUMMARY_SPEC in mobility_value_panel.R for why row 2 — the
#' no-model mobility ranking — is not optional.
build_fig2c_simple <- function(base, h = 2L, lad = fig2_ladder_checked(h)) {
  S <- mobval_summary_frame(lad, h = h, metric = "skill", featured = FEATURED)
  message("[fig2] 2C (simple), h=", h, ": ",
          paste(sprintf("%s = %s (%.2fx)", S$row_key, S$method, S$x), collapse = "; "))
  mobval_summary_panel(S, h = h, base = base, n_breaks = 7)
}

build_fig2 <- function() {
  # Larger base font for Figure 2 only (all three panels): more legible at print size than
  # the 8.6 default shared by the other figures. Hardcoded annotation/axis sizes below are
  # scaled up to match, and panel C's builder scales its own by base / 9.5.
  F2_BASE <- 11.5
  # The metric block on panel A. ev_get() (top level) is the single reader of
  # invasion_evaluation.csv, shared with panel C's cross-check, and it names the column when
  # one is missing rather than letting sprintf render NULL as "character(0)" and silently
  # drop the whole annotation from a published panel.
  gv <- function(col, h) ev_get(col, h)
  # The featured model's 2-week rows as invasion_evaluation.csv scored them (shared support);
  # also the source of the round counts the Figure2_simple caption quotes.
  fs2 <- featured_scored_rows(2L)
  # TWO DISCRIMINATION MEASURES, THEN TWO CALIBRATION MEASURES. AUC-PR skill and AUC-ROC are
  # invariant to any monotone transform of the probabilities, so a model can lead on both
  # while its probabilities are wrong by a factor of two — and the probabilities are used to
  # decide how many zones to prepare, not only which ones. Calibration (pred/obs) is
  # calibration-in-the-large: 1.0 = the expected number of positive zone-rounds matched the
  # realised one.
  #
  # THE TWO GROUPS COVER DIFFERENT ROUNDS, AND THE HEADER SAYS SO. Skill, base rate and AUC-ROC
  # are the evaluation table's, over every scored round. The Brier score and calibration are
  # computed on the probabilities this pass draws, over the rounds that share ONE scale — on
  # the recalibrated pass, the rounds carrying a fitted factor (fig2_calibration_metrics()).
  # The same estimator over every round must reproduce the table's EV_CAL value; if it does
  # not, the rows here are not the table's rows, and the block would describe two row sets.
  #
  # THE BASE RATE SITS DIRECTLY UNDER THE SKILL because it is the skill's denominator: AUC-PR
  # skill is average precision divided by it, so the two lines together give the average
  # precision itself. Without it the reader cannot recover the raw quantity, cannot tell that
  # this is a rare-event problem, and has no way to see why the Brier score below is small for
  # a reason that has nothing to do with the model. It is the same denominator panel C divides
  # by, and comes from the same table and the same shared support.
  # fs_metric_block() (figure_style.R) owns the wording and the order, so this block and the
  # cascade panel's carry the same metric lines in the same sequence.
  cm2 <- fig2_calibration_metrics(fs2)
  local({
    ev_cal <- gv(EV_CAL, 2)
    if (is.finite(ev_cal) && abs(cm2$calib_all - ev_cal) / max(abs(ev_cal), .Machine$double.eps) > 1e-6)
      stop(sprintf(paste0("[fig2] over every round, calibration-in-the-large of %s at h=2 is %.6f ",
                          "on the rows panel A's block is computed from but %.6f in ",
                          "invasion_evaluation.csv (%s): the block would not describe the ",
                          "table's rows. Is the LFO frame from the same run as the table?"),
                   FEATURED, cm2$calib_all, ev_cal, EV_CAL), call. = FALSE)
  })
  message(sprintf("[fig2] 2A block: Brier %.4f, calibration %.3f on %s of %s round(s)%s",
                  cm2$brier, cm2$calib, paste(unique(cm2$n_sel), collapse = "-"),
                  paste(unique(cm2$n_all), collapse = "-"),
                  if (cm2$recal) " (recalibrated rounds only)" else " (raw scale)"))
  disc_lab <- fs_metric_block(fig2_metric_header(cm2, n_skill_rounds = gv("n_folds", 2)),
                              skill = gv(EV_SKILL, 2), chance = gv(EV_CHANCE, 2),
                              base_rate = gv(EV_BASE, 2),
                              roc = gv(EV_ROC, 2), brier = cm2$brier, calib = cm2$calib)

  # --- 2A prioritisation: featured model vs the FIGURE_BASELINE_METHODS comparators and random
  #     (2 weeks). The published 90% fold-cluster bootstrap band is drawn on the featured model
  #     only, to keep the panel legible. ---
  pc <- prioritisation_curves(2L)
  # LEGEND ORDER: featured first, then FIGURE_BASELINE_METHODS. Only the ORDER is taken here —
  # the colours and labels are resolved inside .mk_p2a(), which needs them per variant.
  MODEL_ORD <- c(FEATURED, BASELINE_METHODS)
  curves_b <- pc$curves %>% mutate(method = factor(method, levels = intersect(MODEL_ORD, unique(method))))
  feat_c   <- curves_b %>% filter(method == FEATURED)

  #' Panel 2A. ONE code path for the published panel and the simplified variant, so the two
  #' cannot end up drawing different curves: the data, the geometry, the band, the colours and
  #' the metric block are identical and only the words and one annotation differ.
  #'
  #' @param plain      plain-language legend labels and axis title (the simplified figure).
  #' @param watch_note draw the italic "random watch-list" note on the dashed reference line.
  #'   The LINE itself is drawn either way — only the note is dropped.
  .mk_p2a <- function(plain = FALSE, watch_note = TRUE) {
    # FEATURED_LBL / FEATURED_LBL_PLAIN (figure_style.R) are also what panel C calls the same
    # model, so the two panels name it identically.
    sc <- fs_panelA_scales(FEATURED, if (plain) FEATURED_LBL_PLAIN else FEATURED_LBL,
                           BASELINE_METHODS, plain = plain)
    ggplot(curves_b, aes(k, recall, colour = method)) +
      geom_line(data = pc$rnd, aes(k, recall_random), linetype = "22", colour = FAINT, linewidth = 0.7, inherit.aes = FALSE) +
      geom_ribbon(data = feat_c, aes(ymin = recall_lo, ymax = recall_hi, fill = method), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.9) + geom_point(data = feat_c, size = 1) +
      # `+ NULL` is a no-op in ggplot2, so the simplified panel is the published one with this
      # single layer absent — not a separately maintained plot.
      (if (watch_note)
         annotate("text", x = 15, y = 0.075, label = "random watch-list", colour = MUTED,
                  size = 3.5, angle = 8)
       else NULL) +
      annotate("text", x = 0.5, y = 0.99, hjust = 0, vjust = 1, label = disc_lab, size = 3.4, colour = INK, lineheight = 0.95) +
      scale_colour_manual(values = sc$values, labels = sc$labels, name = NULL, breaks = levels(curves_b$method)) +
      scale_fill_manual(values = sc$values, guide = "none") +
      scale_y_continuous(labels = percent_format(1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
      labs(x = if (plain) "Prioritisation of top K zones"
                else      "Zones actively monitored per round (K)",
           y = "Share of true invasions caught") +
      # ONE KEY PER ROW. This legend has five entries and the longest runs to 39 characters, so
      # laid out in two rows of three it is 7.46 in wide (measured) against the 4.1 in panel A
      # is given - it ran straight across panel B and interleaved with panel B's own legend,
      # which is what a reader of the published figure actually saw. ggplot does not clip or
      # wrap a guide box that is wider than its plot, so the only fix is a narrower legend.
      # Two columns is still 5.14 in and does not fit either; one column is 2.91 in and does,
      # with room to spare. .assert_legend_fits() below re-measures this at build time.
      guides(colour = guide_legend(ncol = 1)) +
      theme_pub(F2_BASE) + theme(legend.position = "top")
  }
  p2a   <- .mk_p2a(plain = FALSE, watch_note = TRUE)    # as published
  p2a_s <- .mk_p2a(plain = TRUE,  watch_note = FALSE)   # Figure2_simple

  # --- 2B forecast-vs-outcome: top-20 predicted-risk zones for 3 evenly-spread rounds (2 weeks) ---
  # THREE-WAY outcome (mirrors F_topk15_ever_h1, adapted to the 2-week horizon): each top-ranked
  # at-risk zone is coloured by whether it was invaded WITHIN this round's 2-week forecast window
  # (is_new_invasion == 1 at h=2 ⇒ its first confirmed case has onset in days (C+6, C+20] for
  # cutoff C, i.e. the 14 days after the forecast date C+6), invaded LATER but still by the last
  # data day used for fitting, DATA_END_DATE (affected_ever, i.e. recorded a first confirmed case
  # by the present run but outside this round's window), or NOT invaded by DATA_END_DATE.
  # Palette/labels match F_topk15_ever_h1 (dark burnt-orange / light amber / grey).
  OUT_COL <- setNames(c("#B33005", "#F6BB6B", "grey80"),
                      c("Invaded within round",
                        sprintf("Invaded later (by %s)", DATA_END_LBL),
                        sprintf("Not invaded (by %s)", DATA_END_LBL)))
  LV <- names(OUT_COL)
  # zones invaded at SOME point by DATA_END_DATE (present-run risk scores; was_active_before
  # marks zones already affected as of the analysis date). horizon == 1 slice is only to de-duplicate.
  # THE PANEL IS A STACKED BAR IN BOTH MODES, so one code path draws both.
  #   single draw : each zone contributes ONE segment of full length — the realised outcome
  #                 class. Geometrically identical to the un-stacked bar this replaced.
  #   MI (M arms) : the bar LENGTH is the MI-mean probability and the SEGMENTS are the share
  #                 of imputations falling in each outcome class. That is the honest drawing,
  #                 because the outcome class is itself imputed: `is_new_invasion` is read off
  #                 a zone's first confirmed ONSET, which moves with the draw. Colouring by a
  #                 single arm's class, or by the modal one, would present a label as certain
  #                 when the label is the very thing under study.
  # `seg` is the drawn width (p x share); `p_mi` the total, used for ordering and the axis.
  MI_F2B <- file.path(OUT, "key_outputs", "Figure2B_mi_data.csv")
  # WHICH ROUNDS CARRY A RECALIBRATION FACTOR. The bars are drawn on the prequentially
  # recalibrated scale, but the earliest rounds had no settled, post-burn-in round to fit delta
  # on by their issue date (delta_estimable FALSE) and are drawn UNCORRECTED. One shared x axis then mixed the two
  # scales with nothing on the page to say so — theme_pub() blanks subtitles, so the disclosure
  # the all-rounds companion (FigureS14) carried in its subtitle was never drawn either. Any such round is now
  # labelled "(uncorrected)" in its facet strip, here and in FigureS14, and the three rounds of
  # panel B are chosen among the rounds that DO carry a factor. Read per round off the
  # featured model's own rows; an LFO frame without the column predates recalibration, and
  # every round is then on the one (raw) scale.
  .est_by_cut <- local({
    d <- lfo %>% filter(method == FEATURED, horizon == 2)
    if (!"delta_estimable" %in% names(d)) d$delta_estimable <- TRUE
    d %>% group_by(cutoff = as.Date(cutoff)) %>%
      summarise(est = any(delta_estimable %in% TRUE), .groups = "drop")
  })
  .uncorr <- function(cut) {
    # On the raw pass nothing is recalibrated, so nothing is "uncorrected" relative to the rest.
    if (!identical(PSCALE_COL, "p_recal")) return(rep(FALSE, length(cut)))
    !(.est_by_cut$est[match(as.Date(cut), .est_by_cut$cutoff)] %in% TRUE)
  }
  # WHICH ROUNDS HAVE NOT SETTLED: their outcome window closed fewer than
  # LFO_RELIABLE_EVAL_AGE_DAYS before the data end (eval_reliable FALSE), so a first case in it
  # may still be confirmed. They are drawn, as the publication suite draws them, but marked
  # "(provisional)" rather than shown as final.
  .rel_by_cut <- local({
    d <- lfo %>% filter(method == FEATURED, horizon == 2)
    if (!"eval_reliable" %in% names(d))
      stop("[fig2] lfo_results.rds has no eval_reliable column, so unsettled rounds cannot be ",
           "marked; re-run the pipeline.", call. = FALSE)
    d %>% group_by(cutoff = as.Date(cutoff)) %>%
      summarise(rel = all(eval_reliable %in% TRUE), .groups = "drop")
  })
  .prov <- function(cut) !(.rel_by_cut$rel[match(as.Date(cut), .rel_by_cut$cutoff)] %in% TRUE)
  .unc_tag <- function(cut) vapply(seq_along(cut), function(i) {
    tg <- c(if (.uncorr(cut[i])) "uncorrected", if (.prov(cut[i])) "provisional")
    if (length(tg)) paste0(" (", paste(tg, collapse = ", "), ")") else ""
  }, character(1))
  # THE BAR AXIS NAMES ITS SCALE. Panel B and FigureS14 draw probabilities, which recalibration
  # moves, and nothing on the primary figures said which scale they were on (fs_caption() stamps
  # only the raw twins by default). Rounds drawn uncorrected are named as such.
  .p_axis <- function(cuts) {
    if (!identical(PSCALE_COL, "p_recal"))
      return("Predicted invasion probability (raw), P(first case)")
    if (any(.uncorr(cuts)))
      return(paste0("Predicted invasion probability, P(first case): recalibrated, ",
                    "except rounds marked (uncorrected), which are raw"))
    "Predicted invasion probability (recalibrated), P(first case)"
  }
  if (MI_MODE) {
    if (!file.exists(MI_F2B))
      stop("[fig2] MI mode but ", basename(MI_F2B), " is absent from ", OUT,
           ". 48_onset_mi_pool.R publishes it; pool the arms before rendering.", call. = FALSE)
    .mi <- suppressWarnings(read_csv(MI_F2B, show_col_types = FALSE))
    stopifnot(all(c("cutoff", "round_index", "health_zone", "p_mi", "outcome", "share",
                    "n_arms_at_risk", "n_imputations") %in% names(.mi)))
    # The class labels are written by the pool step against the SAME DATA_END_DATE this
    # script quotes; if they disagree the legend would describe a different cutoff date.
    .bad <- setdiff(unique(.mi$outcome), LV)
    if (length(.bad))
      stop("[fig2] ", basename(MI_F2B), " carries outcome class(es) this figure does not ",
           "define: ", paste(.bad, collapse = ", "), ". Expected: ", paste(LV, collapse = ", "),
           " (check that the pool step and this script agree on DATA_END_DATE).", call. = FALSE)
    .mi <- .mi %>% mutate(cutoff = as.Date(cutoff))
    fold_ord <- sort(unique(.mi$cutoff))
    lab_lk <- setNames(sprintf("Round %d: %s%s", .mi$round_index[match(fold_ord, .mi$cutoff)],
                               format(lfo_origin(fold_ord), "%d %b"), .unc_tag(fold_ord)),
                       as.character(fold_ord))
    # FLAG A CELL THAT IS NOT AT RISK IN EVERY COMPLETION. A zone already invaded before this
    # round in some completions is ranked here on the ones where it was still at risk, so its
    # bar is an average over FEWER completions than its neighbours'. Unmarked, a bar built
    # from one completion is indistinguishable from one built from twenty.
    .nmax <- max(.mi$n_imputations, na.rm = TRUE)   # M, published by the pool step
    .mi <- .mi %>% mutate(.partial = n_arms_at_risk < .nmax,
                          disp = ifelse(.partial, paste0(health_zone, " †"), health_zone))
    .n_partial <- dplyr::n_distinct(paste(.mi$cutoff, .mi$health_zone)[.mi$.partial])
    topk <- .mi %>%
      mutate(outcome  = factor(outcome, levels = LV),
             fold_lab = factor(lab_lk[as.character(cutoff)], levels = lab_lk[as.character(fold_ord)]),
             seg      = p_mi * share,
             zone_w   = reorder_within(disp, p_mi, fold_lab))
    if (.n_partial)
      message(sprintf("[fig2] %d of %d drawn cells are at risk in fewer than %d completions (marked \u2020)",
                      .n_partial, dplyr::n_distinct(paste(.mi$cutoff, .mi$health_zone)), .nmax))
    message(sprintf("[fig2] 2B from %s: %d round(s), %d zone-rows, %d imputations",
                    basename(MI_F2B), length(fold_ord), dplyr::n_distinct(.mi$health_zone),
                    .nmax))
  } else {
    # zones invaded at SOME point by DATA_END_DATE (present-run risk scores; was_active_before
    # marks zones already affected as of the analysis date). horizon == 1 slice is only to de-duplicate.
    affected_ever <- rs %>% filter(horizon == 1, as.logical(was_active_before) %in% TRUE) %>%
      pull(health_zone) %>% unique()
    dtk <- lfo %>% filter(method == FEATURED, horizon == 2, is.finite(p_invasion), !(was_active_before %in% TRUE)) %>% mutate(cutoff = as.Date(cutoff))
    all_folds <- sort(unique(dtk$cutoff))
    # First / middle / last of the rounds whose recalibration factor was estimable (see
    # .est_by_cut above), so the three facets share one probability scale. Every round of the
    # frame when none qualifies (a raw pass, or a frame predating recalibration).
    # NOTE: 48_onset_mi_pool.R mirrors this rule for the MI twin (Figure2B_mi_data.csv); the two
    # must pick the same rounds.
    elig_folds <- all_folds[!.uncorr(all_folds)]
    if (!length(elig_folds)) elig_folds <- all_folds
    sel_idx  <- unique(round(seq(1, length(elig_folds), length.out = 3)))   # first / middle / last
    fold_ord <- elig_folds[sel_idx]
    dtk <- dtk %>% filter(cutoff %in% fold_ord)
    # honest round numbering: index within the full CV fold sequence, not 1..3
    lab_lk <- setNames(sprintf("Round %d: %s%s", match(fold_ord, all_folds),
                               format(lfo_origin(fold_ord), "%d %b"), .unc_tag(fold_ord)),
                       as.character(fold_ord))
    topk <- dtk %>% group_by(cutoff) %>% slice_max(p_invasion, n = 20, with_ties = FALSE) %>% ungroup() %>%
      mutate(outcome = factor(dplyr::case_when(
                         is_new_invasion == 1            ~ LV[1],   # invaded within this round's 2-week window
                         health_zone %in% affected_ever  ~ LV[2],   # invaded later, but by DATA_END_DATE
                         TRUE                            ~ LV[3]),  # not invaded by DATA_END_DATE
                       levels = LV),
             fold_lab = factor(lab_lk[as.character(cutoff)], levels = lab_lk[as.character(fold_ord)]),
             p_mi   = p_invasion,
             share  = 1,
             seg    = p_invasion,
             zone_w = reorder_within(health_zone, p_invasion, fold_lab))
  }
  p2b <- ggplot(topk, aes(seg, zone_w, fill = outcome)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    facet_wrap(~ fold_lab, scales = "free_y", nrow = 1) + tidytext_scale_y() +
    scale_fill_manual(values = OUT_COL, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA), expand = expansion(mult = c(0, 0.06)), breaks = scales::pretty_breaks(3)) +
    labs(x = if (MI_MODE && exists(".n_partial", inherits = FALSE) && .n_partial > 0)
               paste0(.p_axis(fold_ord), "  ·  † at risk in only some completions")
             else .p_axis(fold_ord), y = NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_pub(F2_BASE) + theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(9, "pt"), strip.clip = "off",
                        axis.text.y = element_text(size = 8.4, colour = INK), axis.text.x = element_text(size = 8.8), legend.position = "top")

  save_dual(p2a, "F2A_prioritisation",     5.2, 4.8)
  save_dual(p2b, "F2B_forecast_vs_outcome", 2.9 * length(fold_ord) + 1.0, 7.0)

  # --- 2C what the mobile-phone mobility data adds, 2 weeks (see build_fig2c above) ------
  # NOT DRAWN ON AN MI PASS. The ladder is a single-draw object: 50_mobility_value.R scores it
  # on the published lfo_results.rds, whose invasion labels come from ONE onset completion.
  # Putting it beside MI-pooled panels A and B would present two different objects as one
  # figure, which is the same reason Figures S4 and S2 are main-text-pass only (build_fig3).
  # The MI twin is therefore the two-panel Figure 2 it has always been.
  # ONE validated ladder frame for BOTH panel-C drawers, so the published full-ladder panel and
  # the simplified four-row one can never be built from different files or different runs.
  .lad  <- if (MI_MODE) NULL else fig2_ladder_checked(2L)
  p2c   <- if (MI_MODE) NULL else build_fig2c(F2_BASE, h = 2L, lad = .lad)
  p2c_s <- if (MI_MODE) NULL else build_fig2c_simple(F2_BASE, h = 2L, lad = .lad)
  if (MI_MODE)
    message("[fig2] MI mode: panel C skipped (the mobility-value ladder is a single-draw ",
            "object); Figure 2 is A|B.")

  F2_W     <- 13.8
  F2_WIDTHS <- c(0.85, 1.5)   # panel A : panel B, as published
  # ROW HEIGHTS ARE RELATIVE, AND ONLY THE PANELS ARE ELASTIC. The two numbers below divide
  # the height left after panel C's own fixed parts (its two stacked legends and its axis,
  # ~1.2 in) are taken off, so the top row renders at ~7.3 in rather than a literal 8.0. It is
  # 8.0 and not the 7.4 the two-panel figure used because panel A's legend went from two rows
  # to five, which is 0.63 in taller (see the guides() note above), and without that back the
  # panels would have paid for the legend as well.
  #
  # Measured on the rendered figure: panels A and B keep their published WIDTH to within 1%
  # (panel A 4.08 in against 4.11 in; panel B's facets 1.74 against 1.76). They are ~20%
  # shorter - panel A's plotting area 4.3 in against 5.5 in - which is the unavoidable price
  # of a third panel inside a figure that is not allowed to grow without limit. Nothing in
  # either panel depends on its height: both axes rescale.
  F2_TOP_H <- 8.0
  # Panel C's height FROM ITS ROWS: ~1.2 in of fixed parts (two stacked legends and the axis at
  # the Figure-2 base font) plus 0.235 in per printed label line, the pitch of the 5.2 in that
  # the original 17-line ladder was drawn at. A fixed height let a longer ladder (more rows, or a
  # time-varying form's two-line featured label) print its labels into one another.
  F2C_H    <- if (is.null(p2c)) 0 else 1.2 + 0.235 * attr(p2c, "n_label_lines")
  # Both rows, each against the width it will actually be laid out in: the top row shares
  # F2_W in the ratio F2_WIDTHS, panel C has the row to itself.
  .assert_legend_fits(list(p2a, p2b), F2_WIDTHS, F2_W, labels = c("2A", "2B"))
  if (!is.null(p2c)) .assert_legend_fits(list(p2c), 1, F2_W, labels = "2C")
  if (!is.null(p2c)) save_dual(p2c, "F2C_mobility_value", F2_W, F2C_H)

  #' Assemble one Figure-2 layout: a panel A, panel B beside it, and a panel C underneath.
  #'
  #' THE TOP ROW IS WRAPPED, AND THAT IS THE POINT.
  #'
  #' PATCHWORK ALIGNS EVERY PANEL IN ONE GRID, however the plots are nested: `(A|B)/C`,
  #' `A+B+C` under a design, and `free(C, type = "label")` all produce a byte-identical
  #' figure (checked). One of that grid's columns is the y-axis TEXT column, sized by the
  #' widest axis text in it - panel C's method labels, 2.2 in of them. Panels A and B then pay
  #' for a column they do not use, and worse, panel A's own "100%" ticks are right-aligned
  #' inside it, so the figure came out with 1.8 in of white space between panel A's axis title
  #' and its axis. Measured on that build: panel A's plotting area 3.42 in against 4.09 in as
  #' a standalone figure, and each of panel B's three facets 1.42 in against 1.76 in.
  #'
  #' wrap_elements(full = ) makes the assembled top row ONE element of the outer layout, so it
  #' is laid out internally exactly as the standalone 13.8-in-wide figure was, and only the
  #' height is imposed from outside. `full` rather than the default `panel`: as a panel the
  #' block would be aligned with panel C's PANEL, i.e. pushed right by the same 2.2 in of row
  #' labels this is here to escape.
  #'
  #' TAGS. The wrapped row carries its own plot_annotation(tag_levels = "A") -> A, B, and
  #' panel C is tagged by hand. Auto-tagging on the OUTER layout would see two elements and
  #' label them A and B.
  .assemble_fig2 <- function(pA, pC, cH) {
    top <- (pA | p2b) + plot_layout(widths = F2_WIDTHS) + plot_annotation(tag_levels = "A")
    if (is.null(pC)) return(top)
    patchwork::wrap_elements(full = top) / (pC + labs(tag = "C")) +
      plot_layout(heights = c(F2_TOP_H, cH))
  }

  fig2 <- .assemble_fig2(p2a, p2c, F2C_H)
  save_dual(fig2, "Figure2", F2_W, F2_TOP_H + if (is.null(p2c)) 0 else F2C_H, dir = FIG_DIR)

  # ---- THE SIMPLIFIED VARIANT, alongside the published figure ----------------
  # SAME DATA, SAME PANELS, FEWER WORDS. Panel A is the same curves with plain-language legend
  # labels and without the "random watch-list" note (the dashed reference LINE stays); panel B
  # is the identical object; panel C is the four-row summary of the very ladder the published
  # panel C draws in full. It is written under its own basename so both are on disk and the
  # manuscript chooses; nothing about Figure2 itself changes.
  if (!MI_MODE) {
    F2S_C_H <- 3.0     # four rows, one legend and one axis at the Figure-2 base font
    .assert_legend_fits(list(p2a_s, p2b), F2_WIDTHS, F2_W, labels = c("2A-s", "2B"),
                        who = "fig2-simple")
    .assert_legend_fits(list(p2c_s), 1, F2_W, labels = "2C-s", who = "fig2-simple")
    save_dual(.assemble_fig2(p2a_s, p2c_s, F2S_C_H), "Figure2_simple",
              F2_W, F2_TOP_H + F2S_C_H, dir = FIG_DIR)

    # ---- THE CAPTION, WRITTEN FROM THE SAME RUN AS THE FIGURE --------------
    # EVERY NUMBER IN IT IS READ OFF THE OBJECTS THE PANELS WERE DRAWN FROM. A caption typed
    # once and kept beside a figure that is rebuilt whenever the pipeline runs is the single
    # easiest way for a manuscript to end up quoting last month's numbers under this month's
    # figure — and unlike a wrong panel, nothing about it looks wrong. Edit the PROSE here;
    # do not paste numbers into it.
    local({
      .r  <- function(k) pc$curves %>% filter(method == FEATURED, k == !!k)
      .pc <- function(x) sprintf("%.0f%%", 100 * x)
      # "07 September" reads as a file stamp, not as prose.
      .dt <- function(d) sub("^0", "", format(as.Date(d), "%d %B %Y"))
      r10 <- .r(10); r20 <- .r(20)
      S2  <- mobval_summary_frame(.lad, h = 2L, metric = "skill", featured = FEATURED)
      .row <- function(key) S2[S2$row_key == key, ]
      # TOP TO BOTTOM IS THE PANEL'S ORDER, WHICH IS BY SKILL — not the order the selection
      # rules happen to be written in. They coincide on this run; they need not, and a caption
      # that lists four rows in the wrong order is wrong in a way no one would spot.
      S2ord <- S2[order(-S2$x), ]
      ld  <- .lad %>% filter(method == FEATURED, horizon == 2L)
      # PANEL A's ROUNDS ARE THE ROUNDS ITS CURVE AVERAGES OVER. compute_detection_curve()
      # averages recall over the rounds on the shared support that hold at least one first
      # case (a round with none has no recall to average). They used to be counted as the
      # featured model's own rounds, which is the same number only while no round drops off the
      # support and every round holds a first case.
      .rounds <- fs2 %>% group_by(cutoff = as.Date(cutoff)) %>%
        summarise(n_pos = sum(is_new_invasion), .groups = "drop")
      .n_scored <- nrow(.rounds); .n_pos_rounds <- sum(.rounds$n_pos > 0)
      ori <- lfo_origin(range(.rounds$cutoff))
      .x1 <- function(v) sprintf("%.1f", v)
      .num <- function(n) if (n >= 1 && n <= 10) c("one", "two", "three", "four", "five", "six",
                                                   "seven", "eight", "nine", "ten")[n] else
                            as.character(n)
      # PANEL A's COMPARATORS, described in the order the legend lists them. The descriptions are
      # keyed by method so a change to FIGURE_BASELINE_METHODS fails here instead of leaving the
      # caption describing a comparator the panel no longer draws.
      .bl_desc <- c(
        "Baseline-flowminder-inflow" = "measured movement of people out of the outbreak area",
        "Baseline-travel-time"       = "road travel time from the outbreak area",
        "Distance-B1"                = "recent cases in other zones divided by the travel time to them",
        # forecast_B7_adjacency() is handed the OSRM travel-time matrix (run_all.R), so it ranks by
        # road travel time in minutes, not by distance.
        "Adjacency-B7"               = "road travel time to the nearest zone that already had cases",
        "Baseline-gravity"           = "modelled movement of people out of the outbreak area")
      if (length(setdiff(BASELINE_METHODS, names(.bl_desc))))
        stop("[fig2] the Figure2_simple caption has no description for comparator(s) ",
             paste(setdiff(BASELINE_METHODS, names(.bl_desc)), collapse = ", "), call. = FALSE)
      .bl_txt <- unname(.bl_desc[BASELINE_METHODS])
      .bl_txt <- if (length(.bl_txt) > 1L)
        paste0(paste(.bl_txt[-length(.bl_txt)], collapse = "; "), "; and ", .bl_txt[length(.bl_txt)])
        else .bl_txt
      # PANEL B's ROUNDS: say how they were chosen when rounds were passed over (see .uncorr()).
      .n_unc <- sum(.uncorr(all_folds))
      .b_how <- if (.n_unc > 0L)
        sprintf(paste0(" (the first, middle and last of the %d rounds whose probabilities could be ",
                       "recalibrated; the first %s could not be, because the factor is fitted only on ",
                       "rounds after the first %d weeks whose outcomes were at least %d days old, and ",
                       "none was yet available when they were issued)"),
                length(all_folds) - .n_unc, if (.n_unc == 1L) "round" else paste(.num(.n_unc), "rounds"),
                # Burn-in = week index < ASCERTAINMENT_BURNIN_WEEKS with week 1 the week holding
                # OUTBREAK_START, so a post-burn-in round comes AFTER the first (that - 1) weeks.
                as.integer(ASCERTAINMENT_BURNIN_WEEKS) - 1L, as.integer(LFO_RELIABLE_EVAL_AGE_DAYS))
        else " (the first, middle and last)"
      # PANEL C's PROSE IS DERIVED FROM THE NUMBERS. It used to say the mobility-data model "did
      # best" and that "most of the advantage lies in the data itself" as fixed text, and on
      # LINELIST_09102026 the first was false (the data-only ranking scored 24.8 against the
      # model's 21.0). Every comparative word below is chosen by comparing the point estimates,
      # and the caption says they are point estimates.
      X <- .row("model_mob")$x; Y <- .row("data_only")$x
      M <- .row("model_nomob")$x; R <- .row("rule")$x; P <- max(M, R)
      s_xy <- if (X > Y)
        paste0("The forecasting model of panels A and B, which uses human mobility data, scored ",
               .x1(X), "; ranking the zones by the movement data alone, ",
               "with nothing fitted to the case counts, scored ", .x1(Y), ".")
        else
        paste0("Ranking the zones by the movement data alone, with nothing fitted to the case ",
               "counts, scored ", if (Y >= max(X, M, R)) "highest, at " else "", .x1(Y),
               "; the forecasting ",
               "model of panels A and B, which uses the same kind of data, scored ", .x1(X), ".")
      s_attr <- if (X <= P)
        paste0(" The forecasting model did not score above the best method without movement data (",
               .x1(P), ").")
        else if (Y >= X)
        paste0(" So the whole of the forecasting model's advantage over the best method without ",
               "movement data (", .x1(P), ") is already present in the data alone, with no model.")
        else if (Y <= P)
        paste0(" The movement data alone did no better than the best method without it (", .x1(P),
               "), so the forecasting model's advantage comes from the modelling rather than from ",
               "the data on its own.")
        else if ((Y - P) / (X - P) > 0.5)
        paste0(" So most of the forecasting model's advantage over the best method without ",
               "movement data (", .x1(P), ") is already present in the data alone.")
        else
        paste0(" So most of the forecasting model's advantage over the best method without ",
               "movement data (", .x1(P), ") comes from the modelling rather than from the data alone.")
      # THE FEATURED MODEL'S FORM. A time-varying form is featured only against its own fixed-beta
      # specification (stage 2); every other fitted row of the ladder holds beta fixed, so the
      # forecasting model's difference from them combines the movement data with the form. Say so,
      # and quote the like-for-like row (its stage-1 specification, on the ladder).
      .fp <- mv_featured_parts(FEATURED)
      .twin <- .lad %>% filter(method == .fp$base, horizon == 2L)
      .fk <- tv_form_of(FEATURED)
      .form_desc <- if (is.na(.fk)) NA_character_ else switch(TV_FORMS[[.fk]]$form,
        frailty     = "a gamma frailty, which lets zones differ in susceptibility",
        rw1         = "a weekly random walk in its import coefficient",
        rw1_frailty = "a weekly random walk in its import coefficient and a gamma frailty",
        NA_character_)
      if (!is.na(.fk) && (is.na(.form_desc) || nrow(.twin) != 1L))
        stop("[fig2] the featured model ", FEATURED, " is a time-varying form, but its form is not ",
             "described or its fixed-beta specification ", .fp$base, " is not on the ladder; the ",
             "caption could not say what panel C compares.", call. = FALSE)
      s_form <- if (is.na(.fk)) "" else
        paste0(" The forecasting model carries ", .form_desc, "; every other fitted row holds its ",
               "import coefficient fixed, so its difference from them combines the movement data with ",
               "that form. The same model with the coefficient fixed scored ", .x1(.twin$skill), ".")
      # Distinct invasions vs positive zone-rounds: at h=2 the cumulative two-week windows of
      # consecutive rounds overlap, so one first case is usually scored from two origins. The
      # ladder's n_events counts the positive ROWS (91 on LINELIST_09102026); the evaluation's
      # n_invasion_events counts distinct invaded zones (50). Quoting the first as "first cases"
      # overstated the evidence nearly twofold.
      # On the ladder's rounds (ALL rounds), from the featured model's scored rows, which ARE the
      # shared support (every method covers every at-risk cell; 50 asserts it).
      n_dist <- dplyr::n_distinct(fs2$health_zone[fs2$is_new_invasion == 1L])
      if (sum(fs2$is_new_invasion == 1L) != ld$n_events || nrow(fs2) != ld$n_cells)
        stop(sprintf(paste0("[fig2] the ladder scores %d positive zone-rounds on %d cells at h=2, ",
                            "but the featured model's scored rows hold %d on %d: the ladder is not ",
                            "on this run's rounds."),
                     ld$n_events, ld$n_cells, sum(fs2$is_new_invasion == 1L), nrow(fs2)),
             call. = FALSE)
      if (n_dist != gv("n_invasion_events", 2))
        stop(sprintf("[fig2] %d distinct first cases in the featured model's rows but %d in the table.",
                     n_dist, as.integer(gv("n_invasion_events", 2))), call. = FALSE)
      .n_all_r <- dplyr::n_distinct(fs2$cutoff)
      # HOW THE MODEL WAS CHOSEN, read from the run's own evaluation (sel_burnin_applied,
      # sel_n_folds, sel_auc_pr_skill_date), not from this session's settings. Under the default
      # basis the selection reads the post-burn-in rounds only — a deliberate choice, stated —
      # while every number in the figure is over all rounds.
      .sel_burn <- isTRUE(as.logical(gv("sel_burnin_applied", 2)))
      .n_sel_r  <- as.integer(gv("sel_n_folds", 2))
      s_choice <- if (.sel_burn)
        paste0("*How the model was chosen.* The forecasting model was chosen on the ", .n_sel_r,
               " rounds after the first ", as.integer(ASCERTAINMENT_BURNIN_WEEKS) - 1L,
               " weeks of the outbreak only; on those rounds its AUC-PR skill was ",
               .x1(gv("sel_auc_pr_skill_date", 2)), " (", .x1(gv(EV_SKILL, 2)), " over all ",
               .n_all_r, "). This was deliberate: in the first ",
               as.integer(ASCERTAINMENT_BURNIN_WEEKS) - 1L, " weeks, which zones had a confirmed ",
               "first case partly reflected where surveillance looked first, as a backlog of zones ",
               "already reached was confirmed at once, so those rounds were kept out of the choice ",
               "of model. They are still forecast and scored, and every measure of ranking in this ",
               "figure (recall at K, AUC-PR skill, AUC-ROC) is over all ", .n_all_r, " rounds.\n\n")
        else paste0("*How the model was chosen.* The forecasting model was chosen on all ",
                    .n_all_r, " rounds, the rounds every measure of ranking in this figure is ",
                    "computed over.\n\n")
      .between <- sprintf("issued between %s and %s", sub("^0", "", format(ori[1], "%d %B")),
                          .dt(ori[2]))
      .avg_over <- if (.n_pos_rounds == .n_scored)
        sprintf("Averaged over the %d weekly rounds %s", .n_scored, .between)
        else sprintf("Averaged over the %d of the %d weekly rounds %s in which at least one zone had its first case",
                     .n_pos_rounds, .n_scored, .between)
      # THE PROBABILITY SCALE OF PANEL B's BARS, and of the block's Brier score and calibration
      # (see fig2_calibration_metrics() and .p_axis()).
      .bar_scale <- if (!identical(PSCALE_COL, "p_recal")) ", as raw (uncorrected) model output"
        else if (any(.uncorr(fold_ord))) ", after recalibration except in rounds marked uncorrected"
        else ", after recalibration"
      .cal_rounds <- if (cm2$recal)
        sprintf("over the %d rounds whose probabilities were recalibrated", cm2$n_sel[1])
        else sprintf("over all %d rounds, on the raw (uncorrected) probabilities", cm2$n_all[1])
      txt <- paste0(
"**Figure 2. Which health zones the outbreak reached next, and how far ahead this could be seen.**\n",
"\n",
# AT RISK AND OUTCOME ARE BOTH ONSET-DATED, ON THE FINAL DATA (run_invasion_lfo(),
# 16_invasion_eval.R): at risk = no confirmed case with onset on or before the forecast date
# (C + 6); outcome = first such onset in the 14 days after it. Neither is about REPORTING: a zone
# whose first case fell ill before the forecast date but was confirmed later was not at risk, and
# a first case confirmed weeks after the window still counts as a hit.
"**(A)** In each weekly forecast round, every health zone with no confirmed case whose symptoms ",
"had begun by the forecast date was ranked by how likely it was to have its first case, dated by ",
"the day symptoms began, in the following two weeks. Both the list of zones and what happened ",
"to them are judged on the final data: a case confirmed only after the forecast date still ",
"counts, on the day its symptoms began (estimated where no onset was recorded). The ",
"curves show how much of what actually happened would have been covered by watching only the ",
"top K zones on that list. ", .avg_over, ", watching the top ",
"10 zones each round would have covered ", .pc(r10$recall), " of the zones that were genuinely ",
"reached (90% interval ", .pc(r10$recall_lo), " to ", .pc(r10$recall_hi), "), and the top 20 ",
"would have covered ", .pc(r20$recall), " (", .pc(r20$recall_lo), " to ", .pc(r20$recall_hi),
"). Choosing zones at random would have covered ", .pc(r10$recall_random), " and ",
.pc(r20$recall_random), ", shown by the grey dashed line. The ", .num(length(BASELINE_METHODS)),
" comparison lines rank zones by simpler things: ", .bl_txt, ". The shaded band is the ",
"uncertainty in the forecasting model's line, obtained by resampling whole rounds.\n",
"\n",
"**(B)** Three of those rounds in detail", .b_how, ": the 20 zones ranked highest in each, and ",
"what actually happened to each one — a first case inside the two weeks the ranking covered ",
"(dark red), a first case later in the outbreak (amber), or no case by ", .dt(DATA_END_DATE),
", the last day of data used to fit (grey). Bar length is the probability the ranking gave ",
"that zone", .bar_scale, ".",
if (any(.prov(fold_ord))) paste0(" A round marked provisional closed fewer than ",
  as.integer(LFO_RELIABLE_EVAL_AGE_DAYS), " days before the end of the data, so a first case ",
  "inside its window may still be confirmed.") else "",
"\n",
"\n",
"**(C)** What the movement data was worth. One row is the forecasting model itself; every ",
"other row is the best method of its kind. All are scored on exactly the ",
"same zones and rounds (all ", .n_all_r, " rounds), ",
"and the number is how well each put that round's newly reached zones at the top of its list: ",
"the average precision within each round, as a multiple of the share of zones that were ",
"reached (a list in random order scores ", .x1(.row("model_mob")$ref), " on these rounds, the ",
"dotted line). ", s_xy, s_attr, s_form,
" The best model with no access to movement data scored ", .x1(M), ", and the best rule ",
"without it (no renewal model) ", .x1(R), ". These comparisons are between point estimates: ",
"bars are 90% intervals for each row on its own, not a test of the difference between rows.\n",
"\n",
"*Metrics.* Panel A plots recall at K, averaged over rounds, with 90% intervals from a cluster ",
"bootstrap over rounds; its text block reports AUC-PR skill (average precision within each ",
"round, combined over rounds weighted by their first cases and divided by the base rate, beside ",
"the score of a random ordering), the base rate and AUC-ROC over all ", .n_all_r,
" rounds, and the Brier score and ",
"calibration-in-the-large ", .cal_rounds, ". Calibration-in-the-large is the predicted number ",
"of positive zone-rounds divided by the observed number, where a zone-round is an at-risk zone ",
"at a forecast round and is positive when that zone's first case fell inside the round's ",
"window. Panel C plots the same AUC-PR skill with 90% intervals from a paired ",
"zone-cluster bootstrap (", format(ld$n_boot, big.mark = ","), " resamples) over the ",
format(ld$n_cells, big.mark = ","), " zone-rounds scored ",
"at the two-week horizon in all ", .n_all_r, " rounds, in which ", n_dist,
" zones had their first case inside a scored ",
"window; because the two-week windows of consecutive rounds overlap, those ", n_dist,
" first cases make ", ld$n_events, " positive zone-rounds. Its rows are, from top to bottom: ",
paste(sprintf("%s = %s", S2ord$row_label, S2ord$method), collapse = "; "), ".\n",
"\n", s_choice)
      f <- file.path(FIG_DIR, "Figure2_simple_caption.md")
      writeLines(txt, f)
      message("  wrote ", basename(f))
    })
  } else {
    message("[fig2] MI mode: Figure2_simple skipped (its panel C is the single-draw ladder).")
  }

  # ---- SI: the SAME panel for EVERY round, not a three-round sample ----------
  # MAIN-TEXT PASS ONLY, for the same reason Figures S4 and S2 are (see build_fig3). This
  # companion colours each bar by a SINGLE completion's outcome class, and it needs every
  # round rather than the three the pooled Figure2B_mi_data.csv carries. Its MI analogue
  # would be a stacked-share panel with one facet per round, which nothing has specified;
  # drawing it from one arm beside MI-pooled probabilities would be worse than not drawing it.
  #
  # NUMBERED S14. It was "FigureS5_topk_all_rounds", which put two different supplementary
  # figures under one number: make_si_model_figures.R publishes FigureS5_roc_pr_curves. S14 is
  # the next number after the last one the suite assigns (FigureS13_prior_posterior). The
  # FIGURE_KEEP entry in 00_config.R must name the same stem, or save_dual() drops it - which
  # is checked below rather than left to happen silently.
  S14_STEM <- "FigureS14_topk_all_rounds"
  if (MI_MODE) {
    message("[fig2] MI mode: ", S14_STEM, " skipped (single-completion outcome classes).")
    return(invisible(fig2))
  }
  # Gated on the bare STEM: on a FORECAST_SCALE=raw pass FIG_DIR is a raw/ sibling, which the
  # gate refuses by design, and that refusal is not the defect this check is looking for.
  if (exists("figure_is_kept", mode = "function") && !isTRUE(figure_is_kept(S14_STEM)))
    stop("[fig2] ", S14_STEM, " is not on the FIGURE_KEEP allow-list (00_config.R), so ",
         "save_dual() would build it and silently write nothing. Rename its entry there ",
         "(it was \"FigureS5_topk_all_rounds\").", call. = FALSE)

  # Figure 2B shows first / middle / last so the main figure stays legible, which means the
  # reader sees three of the rounds and has to take the rest on trust. This companion
  # shows all of them. Rounds whose recalibration factor is not yet estimable (the earliest
  # ones, which had no settled, post-burn-in round to fit delta on by their issue date) are
  # drawn UNCHANGED rather than omitted -- an uncorrected forecast is still a forecast, and its
  # top-20 ordering is exactly what the operational watch-list would have been. delta is a
  # monotone transform, so within a round it cannot move the ranking at all; only the x-axis
  # scale differs between a recalibrated round and an uncorrected one, and each uncorrected
  # round says so in its own facet strip (.unc_tag(), above) and the axis title says which
  # scale the rest are on (.p_axis()). That disclosure used to be a subtitle, which theme_pub()
  # blanks, so the figure carried no trace of it.
  dall <- lfo %>%
    filter(method == FEATURED, horizon == 2, is.finite(p_invasion),
           !(was_active_before %in% TRUE)) %>%
    mutate(cutoff = as.Date(cutoff))
  if (nrow(dall)) {
    af <- sort(unique(dall$cutoff))
    lab_all <- setNames(sprintf("Round %d: %s%s", seq_along(af), format(lfo_origin(af), "%d %b"),
                                .unc_tag(af)),
                        as.character(af))
    tk <- dall %>% group_by(cutoff) %>% slice_max(p_invasion, n = 20, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(outcome = factor(dplyr::case_when(
                         is_new_invasion == 1           ~ LV[1],
                         health_zone %in% affected_ever ~ LV[2],
                         TRUE                           ~ LV[3]), levels = LV),
             fold_lab = factor(lab_all[as.character(cutoff)],
                               levels = lab_all[as.character(af)]),
             zone_w = reorder_within(health_zone, p_invasion, fold_lab))
    nhit <- tk %>% group_by(fold_lab) %>%
      summarise(hits = sum(outcome == LV[1]), .groups = "drop")
    # For the caption (the figure carries none): how much of the scored window the top 20
    # caught. Positive zone-rounds, not distinct invasions — at h=2 one first case is scored in
    # the two overlapping windows that contain it.
    message(sprintf(paste0("[fig2] %s: %d rounds (%d uncorrected); %d of %d positive zone-rounds ",
                           "were inside that round's top 20."),
                    S14_STEM, length(af), sum(.uncorr(af)), sum(nhit$hits),
                    sum(dall$is_new_invasion, na.rm = TRUE)))
    .ncol <- min(4L, length(af))
    pall <- ggplot(tk, aes(p_invasion, zone_w, fill = outcome)) +
      geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
      facet_wrap(~ fold_lab, scales = "free_y", ncol = .ncol) + tidytext_scale_y() +
      scale_fill_manual(values = OUT_COL, name = NULL, drop = FALSE) +
      scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                         expand = expansion(mult = c(0, 0.06)),
                         breaks = scales::pretty_breaks(3)) +
      labs(x = .p_axis(af), y = NULL) +
      guides(fill = guide_legend(nrow = 1)) +
      theme_pub(F2_BASE - 1) +
      theme(panel.grid.major.y = element_blank(),
            panel.spacing.x = unit(9, "pt"), strip.clip = "off",
            axis.text.y = element_text(size = 7.2, colour = INK),
            axis.text.x = element_text(size = 7.6), legend.position = "top")
    save_dual(pall, S14_STEM,
              3.0 * .ncol + 1.0, 2.5 * ceiling(length(af) / .ncol) + 1.4, dir = FIG_DIR)
  }
  invisible(fig2)
}

# =============================================================================
# FIGURE 3  (two-week operational forecast)
# =============================================================================
# Narrative panel that EXPLAINS the rank evolution. It is deliberately built as the
# mechanistic mirror of build_rank_evolution() via the SAME build_rank_traj(): the
# SAME top-20 present-run at-risk zones, the SAME timeline (CV rounds + the present
# live run), and the SAME province colours. For each zone/round it plots the great-
# circle distance to the nearest ALREADY-invaded zone -- i.e. how close the advancing
# front is. As that distance collapses here, the zone's predicted-risk rank climbs in
# the paired panel: front closes in (A) -> risk rank rises (B). Every zone can be
# traced by name/colour across both, all the way up to the latest forecast.
.zone_centroids <- NULL   # cached lon/lat lookup (keyed by shp .key)
.get_centroids <- function() {
  if (is.null(.zone_centroids)) {
    ct <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(shp))))
    .zone_centroids <<- data.frame(k = shp$.key, lon = ct[, 1], lat = ct[, 2]) %>%
      dplyr::distinct(k, .keep_all = TRUE)
  }
  .zone_centroids
}
.hav_km <- function(lo1, la1, lo2, la2) {
  R <- 6371; dl <- (lo2 - lo1) * pi/180; dp <- (la2 - la1) * pi/180
  a <- sin(dp/2)^2 + cos(la1*pi/180) * cos(la2*pi/180) * sin(dl/2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}
build_front_approach <- function(base = 12, pal = NULL) {
  # tracked zones + timeline come from the SAME shared trajectory as build_rank_evolution:
  # top-20 at-risk zones as of the present run, across CV rounds + the present live run.
  tj <- build_rank_traj()
  keep_z <- tj$keep_z; present_date <- tj$present_date
  times  <- c(tj$cutoffs, present_date)   # CV rounds + the present (live) forecast
  if (is.null(pal)) pal <- province_palette(prov_lookup$province[prov_lookup$health_zone %in% keep_z])

  coord_lk <- .get_centroids()
  xy <- function(zk) { r <- coord_lk[coord_lk$k == zk, ]; if (nrow(r)) c(r$lon[1], r$lat[1]) else c(NA, NA) }
  fo <- zone_first %>% transmute(k = key(health_zone), first)   # first-onset (invasion) date per zone

  rows <- list()
  for (z in keep_z) {
    zk <- key(z); zc <- xy(zk); if (any(is.na(zc))) next
    # Indexed, not `for (t in times)`: a for-loop strips the Date class, which is what used to
    # force an as.Date(t, origin = "1970-01-01") epoch literal back in here.
    for (ti in seq_along(times)) {
      t <- times[ti]
      inv_keys <- fo$k[fo$first <= t & fo$k != zk]
      cc <- coord_lk[coord_lk$k %in% inv_keys, ]; if (!nrow(cc)) next
      rows[[length(rows) + 1L]] <- data.frame(health_zone = z, cutoff = t,
        dist = min(.hav_km(zc[1], zc[2], cc$lon, cc$lat), na.rm = TRUE))
    }
  }
  df <- dplyr::bind_rows(rows) %>% left_join(prov_lookup, by = "health_zone") %>%
    mutate(province = factor(province, levels = pal$levels))

  # key-invasion event markers: date of first confirmed onset for notable zones (dotted
  # verticals) — validation checkpoints as the front advances up to the present run.
  ev <- tibble::tibble(
    health_zone = c("Komanda", "Nia Nia", "Makiso Kisangani", "Isiro", "Rethy"),
    lab         = c("Komanda", "Nia Nia", "Kisangani", "Isiro", "Rethy")) %>%   # spelt as in Figure 1A
    mutate(k = key(health_zone)) %>% left_join(fo, by = "k") %>% filter(!is.na(first)) %>% arrange(first)
  ymax  <- max(df$dist, na.rm = TRUE)
  x_min <- min(times); x_end <- max(c(present_date, ev$first)) + 4

  ggplot(df, aes(cutoff, dist, group = health_zone, colour = province)) +
    annotate("rect", xmin = present_date - 2, xmax = present_date + 2, ymin = -Inf, ymax = Inf, fill = "grey95") +
    geom_vline(xintercept = present_date, linetype = "22", colour = MUTED, linewidth = 0.45) +
    geom_vline(data = ev, aes(xintercept = first), inherit.aes = FALSE,
               linetype = "dotted", colour = "grey40", linewidth = 0.45) +
    geom_line(linewidth = 1.0, alpha = 0.55, lineend = "round") +
    geom_point(colour = "white", size = 2.4) + geom_point(size = 1.45) +
    geom_text(data = ev, aes(x = first, y = ymax, label = lab), inherit.aes = FALSE,
              angle = 90, hjust = 1, vjust = -0.35, size = 2.55, colour = "grey20", fontface = "italic") +
    annotate("text", x = present_date, y = ymax * 0.99, label = "latest run",
             hjust = 0.5, vjust = 1, size = 2.5, colour = MUTED, fontface = "italic") +
    scale_colour_manual(values = pal$colours, name = "Province", breaks = pal$levels, drop = FALSE, na.value = FAINT) +
    scale_x_date(date_labels = "%d %b", breaks = .round_breaks(times, present_date),
                 limits = c(x_min - 2, x_end), expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.03, 0.18))) +
    labs(x = "Forecast round  ·  dotted = key first-onset dates", y = "Distance to nearest\ninvaded zone (km)") +
    theme_pub(base) + theme(legend.position = "right", panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
                            panel.grid.major.x = element_blank(),
                            axis.text.x = element_text(angle = 45, hjust = 1, size = base - 4.2))
}

# -----------------------------------------------------------------------------
# ALERTS AND LABORATORY OUTCOMES in the currently top-ranked zones
# (Figure3_rank_surveillance, panel D)
# -----------------------------------------------------------------------------
# For each of the top-N at-risk zones of the live run — the zones of Figure 3B, in the same
# order — the number of alerts between OUTBREAK_START and the live run's ANALYSIS_DATE, and
# how many of them have a laboratory OUTCOME (confirmed or discarded). An alert is in the window
# when its case date (the onset the case data use) OR its notification date falls in it, so the
# alerts load_linelist() leaves out for lacking an in-outbreak onset or specimen date are
# counted by their notification date (zone_alert_counts(), 01_data_prep.R). Counted by
# zone_alert_counts() (01_data_prep.R), the definition the ascertainment appendix (45,
# zone_ascertainment.csv) also uses, so this panel and the appendix cannot quote different
# counts for one zone. Counts are alert RECORDS, not de-duplicated people.
#
# Every guard below protects a number printed on the panel:
#   * SAME ZONES AS B AND C. B ranks `rs` itself and C goes through build_rank_traj(); the
#     panel refuses to draw if the two sets, or their live-run order, differ.
#   * SAME DATE AND SNAPSHOT AS THE FORECAST. The alerts are counted to ANALYSIS_DATE from
#     the line list on disk; the forecast is the live run in run_info.json. A line list
#     re-pulled after the run would put newer alerts beside an older forecast.
#   * NO SILENT ZERO FROM A SPELLING. A zone the line list spells differently from the spine
#     (as it does "Bili" for "Bili (Bas-Uele)") would show too few alerts rather than fail;
#     refused when an unmatched name equals a plotted zone's name, with or without its
#     "(Province)" suffix. It cannot see other spellings, nor records removed upstream by
#     health_zone_unmatched.
#   * EVERY ALERT ON RECORD IS PLACED. A plotted zone's alert that lands outside the window on
#     both dates (or has no notification date) is refused rather than silently left out: these
#     zones are drawn as complete, and on this snapshot every one of their records is in.
#   * NO REPEATED ALERT ID. A duplicated alert_id among the plotted zones is refused. This is
#     a check on the ID only: two registrations of one person under different IDs are two
#     records here (the panel counts records, as the appendix does).
#   * NOT YET INVADED. An at-risk zone cannot already hold a confirmed case, so a confirmed
#     alert here means the live run and the line list disagree. That invariant is also what
#     lets the legend call every laboratory-resolved alert here "discarded".
#
# @param top   the top-N table drawn in Figure 3B (health_zone, province, zone factor).
# @param pal   the shared province palette (Figure 3B's).
# @param base  theme_pub() base size.
# @return list(plot = ggplot, data = the per-zone counts drawn).
build_zone_alerts <- function(top, pal, base = 12) {
  tj <- build_rank_traj(N_TOP = nrow(top))
  if (!setequal(tj$keep_z, top$health_zone))
    stop("[fig3-alerts] the top-zone set differs between panel B and panel C (build_rank_traj): ",
         paste(sort(setdiff(union(tj$keep_z, top$health_zone),
                            intersect(tj$keep_z, top$health_zone))), collapse = ", "),
         call. = FALSE)
  live_rank <- tj$traj %>% filter(cutoff == tj$present_date) %>% select(health_zone, rank_live = rank)
  top <- top %>% left_join(live_rank, by = "health_zone")
  if (anyNA(top$rank_live) || !identical(as.integer(top$rank_live), seq_len(nrow(top))))
    stop("[fig3-alerts] panel B's order is not the live-run rank 1..", nrow(top),
         " that panel C ends on (ties or a different ranking).", call. = FALSE)

  ri <- jsonlite::fromJSON(file.path(OUT, "key_outputs", "run_info.json"))
  # Compared as VALUES: identical() on two Dates is FALSE when one is integer-backed.
  if (!isTRUE(PRESENT_DATE == ANALYSIS_DATE) || !isTRUE(tj$present_date == ANALYSIS_DATE))
    stop(sprintf("[fig3-alerts] the live run is dated %s but the line list is counted to %s; ",
                 format(PRESENT_DATE), format(ANALYSIS_DATE)),
         "alerts and forecast would describe different days.", call. = FALSE)
  ll_folder  <- jsonlite::fromJSON(LINELIST_JSON)$folder
  run_folder <- ri$linelist_snapshot$folder
  if (is.null(run_folder) || !identical(as.character(ll_folder), as.character(run_folder)))
    stop("[fig3-alerts] the line list on disk (", ll_folder, ") is not the snapshot the live ",
         "run used (", if (is.null(run_folder)) "unrecorded" else run_folder, ").", call. = FALSE)

  und <- load_unplaced_alerts()
  za  <- zone_alert_counts(ll, und, from = OUTBREAK_START, to = ANALYSIS_DATE)
  unmatched <- setdiff(za$health_zone, prov_lookup$health_zone)
  base_nm   <- sub("\\s*\\([^)]*\\)$", "", top$health_zone)
  clash     <- unmatched[key(unmatched) %in% c(key(top$health_zone), key(base_nm))]
  if (length(clash))
    stop("[fig3-alerts] line-list zone name(s) ", paste(clash, collapse = ", "),
         " do not key onto the spine but match a plotted zone; their alerts would be ",
         "dropped and the zone drawn with too few. Add the alias.", call. = FALSE)
  if (length(unmatched))
    message(sprintf("[fig3-alerts] %d line-list zone name(s) not on the spine (none plotted): %s",
                    length(unmatched), paste(unmatched, collapse = ", ")))
  .inw <- function(x) !is.na(x) & x >= OUTBREAK_START & x <= ANALYSIS_DATE
  .not_sit <- !grepl("^SITREP-CONF-", ll$alert_id)
  .out <- c(sum(ll$health_zone %in% top$health_zone & .not_sit &
                  !(.inw(ll$date_index) | .inw(ll$date_of_notification))),
            sum(und$health_zone %in% top$health_zone & !.inw(und$date_of_notification)))
  if (sum(.out))
    stop("[fig3-alerts] ", sum(.out), " alert record(s) in plotted zones fall outside the window ",
         "on every date; the panel would silently omit them.", call. = FALSE)
  ids <- c(ll$alert_id[ll$health_zone %in% top$health_zone & !grepl("^SITREP-CONF-", ll$alert_id)],
           und$alert_id[und$health_zone %in% top$health_zone])
  if (anyDuplicated(ids[!is.na(ids)]))
    stop("[fig3-alerts] duplicated alert IDs among the plotted zones; the counts would ",
         "overstate the number of alerts.", call. = FALSE)

  dz <- top %>% select(health_zone, province, zone, rank_live, p_case_invasion) %>%
    left_join(za, by = "health_zone") %>%
    mutate(across(c(alerts_total, n_confirmed, n_not_a_case, alerts_resolved), ~ coalesce(.x, 0L)))
  stopifnot(nrow(dz) == nrow(top), all(dz$alerts_resolved <= dz$alerts_total))
  if (any(dz$n_confirmed > 0L))
    stop("[fig3-alerts] zone(s) ranked at risk of a FIRST case already hold a confirmed alert: ",
         paste(dz$health_zone[dz$n_confirmed > 0L], collapse = ", "), call. = FALSE)
  message(sprintf(paste0("[fig3-alerts] top %d zones: %d alerts (%d outside the line list's ",
                         "case dates, placed by notification), %d laboratory-resolved ",
                         "(%d discarded, 0 confirmed); %d zone(s) with no alert"),
                  nrow(dz), sum(dz$alerts_total),
                  sum(und$health_zone %in% top$health_zone & .inw(und$date_of_notification)),
                  sum(dz$alerts_resolved), sum(dz$n_not_a_case), sum(dz$alerts_total == 0L)))

  lev <- levels(top$zone)
  ylab_cols <- unname(pal$colours[as.character(top$province)[match(lev, as.character(top$zone))]])
  dz <- dz %>% mutate(count_lab = sprintf("%d / %d", alerts_resolved, alerts_total))
  dp <- filter(dz, alerts_total > 0L)
  d0 <- filter(dz, alerts_total == 0L)
  # "to", not an en dash: the pdf() device cannot encode U+2013 and prints a hyphen.
  win_lab <- sprintf("%s to %s", format(OUTBREAK_START, "%d %b"), format(ANALYSIS_DATE, "%d %b %Y"))
  K_ALL <- "All alerts"; K_RES <- "Laboratory-resolved (all discarded)"
  # The axis ends on the first 1-3-10 tick at or above the largest count, so the top tick is
  # always drawn (limits ending at the data maximum would drop it).
  all_brk <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000)
  x_top <- all_brk[which(all_brk >= max(3L, dz$alerts_total))[1]]
  if (is.na(x_top)) stop("[fig3-alerts] alert count beyond the axis break set.", call. = FALSE)
  x_brk <- all_brk[all_brk <= x_top]
  # THE COUNT COLUMN IS A LAYER ON THE SAME y AESTHETIC AS THE POINTS, inside the panel, past
  # the last tick (no gridline reaches it). It is NOT a secondary discrete axis: dup_axis()
  # on scale_y_discrete hands its label function the breaks in an order that does not match
  # the drawn positions, and printed every zone's counts against another zone (Rethy read
  # "1 / 4", Banalia's numbers). A layer mapped to `zone` cannot be misaligned.
  tr    <- scales::pseudo_log_trans(sigma = 1, base = 10)
  x_col <- tr$inverse(tr$transform(x_top) * 1.27)   # right edge of the count column
  n_row <- length(lev)

  p <- ggplot(dz, aes(y = zone)) +
    geom_segment(data = dp, aes(x = alerts_resolved, xend = alerts_total, yend = zone),
                 colour = FAINT, linewidth = 1.1, lineend = "round") +
    geom_point(data = dp, aes(x = alerts_total, fill = K_ALL),
               shape = 21, colour = INK, size = 3.0, stroke = 0.9) +
    geom_point(data = dp, aes(x = alerts_resolved, fill = K_RES),
               shape = 21, colour = "white", size = 2.5, stroke = 0.5) +
    geom_text(data = d0, aes(x = 0, label = "no alerts"), hjust = 0,
              size = 3.1, colour = MUTED, fontface = "italic") +
    geom_text(aes(x = x_col, label = count_lab), hjust = 1, size = 9.5 / .pt, colour = INK) +
    annotate("text", x = x_col, y = n_row + 0.9, label = "resolved / alerts",
             hjust = 1, vjust = 0, size = (base - 2.6) / .pt, colour = MUTED) +
    scale_fill_manual(values = setNames(c("white", INK), c(K_ALL, K_RES)),
                      breaks = c(K_ALL, K_RES), name = NULL,
                      guide = guide_legend(override.aes = list(colour = INK, size = c(3.0, 2.5),
                                                               stroke = c(0.9, 0.5)))) +
    scale_x_continuous(trans = tr, breaks = x_brk,
                       labels = label_number(accuracy = 1, big.mark = ","),
                       limits = c(0, x_col), expand = expansion(mult = c(0.02, 0.01))) +
    scale_y_discrete(limits = lev, expand = expansion(add = c(0.6, 1.5))) +
    labs(x = sprintf("Alerts, %s (pseudo-log scale)", win_lab), y = NULL) +
    theme_pub(base) +
    theme(legend.position = "top",
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(colour = ylab_cols, size = 9.5, face = "bold"))
  list(plot = p, data = dz)
}


build_fig3 <- function() {
  BIP <- BIV_PROBVUL
  # BOTH QUARTILES OVER THE SAME ZONES — the at-risk set. rr01_nat is NA for an already-affected
  # zone, so its quartiles were cut over the 452 at-risk zones while vulnerability's were cut
  # over all 519, affected zones included: the two axes of one bivariate class meant "top
  # quarter" of two different populations. Affected zones carry no class (they are overlaid in
  # AFFECTED_FILL), so they are excluded from both rankings.
  prep <- function(h) rs %>% filter(horizon == h) %>%
    mutate(active = was_active_before %in% TRUE, rr_log = ifelse(active, NA, pmax(rr01_nat, RR_FLOOR)),
           bx = qtile(ifelse(active, NA_real_, rr01_nat)), by = qtile(ifelse(active, NA_real_, V)),
           bikey = ifelse(is.na(bx) | is.na(by), NA, paste0(bx, "-", by)))
  # THE GREY OVERLAY HAS A LEGEND KEY. Zones already affected at the live run carry no forecast
  # and are drawn in AFFECTED_FILL over the map; nothing on the page used to say what that grey
  # meant (Figure 1A gives its "No confirmed cases" fill a swatch). The key rides on the COLOUR
  # aesthetic because each map's fill is already taken by its own scale (continuous on 3A,
  # bivariate on S2A); the border stays white, and the key alone is outlined so it shows.
  AFF_LBL <- "Already affected (not forecast)"
  aff_overlay <- function(mp) list(
    geom_sf(data = mp %>% filter(was_active_before %in% TRUE), aes(colour = AFF_LBL),
            fill = AFFECTED_FILL, linewidth = 0.08),
    scale_colour_manual(values = stats::setNames("white", AFF_LBL), name = NULL,
                        guide = guide_legend(order = 2, override.aes = list(colour = "grey60"),
                                             theme = theme(legend.key.width = unit(13, "pt"),
                                                           legend.key.height = unit(13, "pt")))))

  # --- 3A relative-invasion-risk map, 2 weeks ahead ---
  mp2 <- join_map(prep(2))
  sc_rr <- scale_fill_viridis_c(option = "viridis", trans = "log10", limits = c(RR_FLOOR, 1),
             breaks = c(0.001, 0.01, 0.1, 1), labels = c("<=0.001","0.01","0.1","1"),
             na.value = NA_FILL, name = "Relative\ninvasion risk\n(log, 0-1)", direction = 1)
  p3a <- ggplot(mp2) + geom_sf(aes(fill = rr_log), colour = "white", linewidth = 0.08) + aff_overlay(mp2) + sc_rr +
    theme_map(12) + theme(legend.position = "right")

  # --- 3B top-20 zones by 2-week invasion probability + 90% CrI, coloured by province ---
  TOP_N <- 20L
  d <- rs %>% filter(horizon == 2, !(was_active_before %in% TRUE), is.finite(p_case_invasion)) %>%
    distinct(health_zone, province, .keep_all = TRUE) %>% arrange(desc(p_case_invasion)) %>% head(TOP_N)
  # shared province palette (also drives Figure S4 and Figure3_rank_surveillance C/D, so every
  # province matches across panels)
  pal <- province_palette(d$province); prov_levels <- pal$levels; prov_cols <- pal$colours
  # disambiguate any duplicate zone names across provinces (else the y factor has duplicate levels)
  d <- d %>% mutate(province = factor(province, levels = prov_levels),
                    disp = ifelse(duplicated(health_zone) | duplicated(health_zone, fromLast = TRUE),
                                  paste0(health_zone, " (", as.character(province), ")"), health_zone))
  # MI: FLAG A ZONE ALREADY INVADED IN SOME COMPLETIONS. Under the pooling rule such a zone
  # contributes p = 0 from those completions rather than being dropped (the estimand is
  # "P(first case in the next two weeks)", which is genuinely zero once invaded), so its
  # interval mixes a forecast with a structural zero. That is the right estimate and the
  # wrong thing to leave unlabelled beside zones at risk in all M.
  .n_partial3 <- 0L
  if (MI_MODE && "n_arms_at_risk" %in% names(d)) {
    # M from the sweep's own provenance, not inferred as max(n_arms_at_risk) over the table:
    # the two agree only as long as SOME zone is at risk in every completion.
    .nmax3 <- as.integer(MI_META$n_imputations)
    .p3 <- is.finite(d$n_arms_at_risk) & d$n_arms_at_risk < .nmax3
    .n_partial3 <- sum(.p3)
    d$disp[.p3] <- paste0(d$disp[.p3], " †")
    if (.n_partial3)
      message(sprintf("[fig3] %d of %d top zones are at risk in fewer than %d completions (marked \u2020)",
                      .n_partial3, nrow(d), .nmax3))
  }
  d$disp <- make.unique(d$disp)                       # final guard: guarantee unique factor levels
  d$zone <- factor(d$disp, levels = rev(d$disp))
  lev <- levels(d$zone); ylab_cols <- unname(prov_cols[as.character(d$province)[match(lev, d$disp)]])
  # THE AXIS NAMES ITS PROBABILITY SCALE, read off the risk table's own stamp (prob_scale:
  # "recalibrated-deployed", its MI twin, or "raw") rather than off the requested pass, because
  # fs_risk_csv() falls back to the primary table when a raw twin is missing.
  .ps3 <- unique(as.character(d$prob_scale))
  .scale3 <- if (length(.ps3) == 1L && grepl("^recalibrated", .ps3)) "recalibrated"
             else if (length(.ps3) == 1L && grepl("^raw", .ps3)) "raw"
             else stop("[fig3] the live risk table carries no single prob_scale stamp (",
                       paste(.ps3, collapse = ", "), "), so panel B cannot say which ",
                       "probability scale it draws. Re-run the pipeline.", call. = FALSE)
  .xlab3 <- sprintf("Probability of invasion (%s), 2-week forecast", .scale3)
  p3b <- ggplot(d, aes(p_case_invasion, zone, colour = province)) +
    geom_linerange(aes(xmin = p_case_lo, xmax = p_case_hi), linewidth = 1.4, alpha = 0.5) + geom_point(size = 2.4) +
    scale_colour_manual(values = prov_cols, name = "Province", drop = FALSE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA), expand = expansion(mult = c(0.01, 0.06))) +
    labs(x = if (MI_MODE && .n_partial3 > 0)
               paste0(.xlab3, "  \u00b7  \u2020 already invaded in some completions")
             else .xlab3, y = NULL) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2.2, size = 3))) +
    theme_pub(12) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = ylab_cols, size = 9.5, face = "bold"), legend.position = "top")

  # --- 3C prob x vulnerability bivariate choropleth, 2 weeks ahead + legend ---
  p3c_map <- ggplot(mp2) + geom_sf(aes(fill = bikey), colour = "white", linewidth = 0.08) + aff_overlay(mp2) +
    scale_fill_manual(values = BIP, na.value = NA_FILL, guide = "none") + theme_map(15.5) +
    # The overlay's key goes UNDER the map: to its right it would sit between the map and the
    # bivariate legend and squeeze both.
    theme(legend.position = "bottom")
  legend_biv <- bivar_legend(BIP, "Invasion prob ->", "Vulnerability ->", base = 12)
  p3c <- (p3c_map | wrap_elements(full = legend_biv)) + plot_layout(widths = c(1, 0.42))

  # --- 3D preparedness-priority scatter, 2 weeks ahead ---
  # Coloured on the suite's identity provinces (prov_region_scale(), figure_style.R), not on
  # PROV_INT: collapsing every province outside the three of interest to "Other" drew the
  # Tshopo and Bas-Uele zones grey here while panel B of Figure 3 draws them in their own hues.
  psc <- rs %>% filter(horizon == 2, !(was_active_before %in% TRUE), is.finite(V), is.finite(rr01_nat))
  .rg <- prov_region_scale(psc$province)
  psc$region <- factor(.rg$region, levels = .rg$levels)
  lab_ps <- psc %>% arrange(desc(priority)) %>% head(12)
  p3d <- ggplot(psc, aes(V, rr01_nat)) +
    geom_point(aes(size = priority, colour = region), alpha = 0.78) +
    ggrepel::geom_text_repel(data = lab_ps, aes(label = health_zone), size = 3.9, colour = INK, box.padding = 0.3,
      max.overlaps = 30, seed = 1, min.segment.length = 0, segment.colour = FAINT, segment.size = 0.2) +
    scale_colour_manual(values = .rg$colours, name = "Province", breaks = .rg$levels) +
    scale_size_area(max_size = 7, name = "Priority") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = label_number(accuracy = 0.01)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = "Vulnerability", y = "Relative invasion risk") +
    theme_pub(15.5) + theme(legend.position = "right")

  # MAIN TEXT Figure 3: relative-invasion-risk map (A) + top-zone uncertainty forest (B), 2-week
  # horizon. Distinct legends per panel -> wrap_elements keeps each panel's own guides/layout.
  save_dual(p3a,  "F3A_relative_risk_map_h2", 6.4, 4.6)   # main-text Fig 3A
  save_dual(p3b,  "F3B_top_zone_uncertainty", 5.6, 4.6)   # main-text Fig 3B
  fig3 <- (wrap_elements(full = p3a) | wrap_elements(full = p3b)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 12 + 4.5, face = "bold", colour = INK))
  save_dual(fig3, "Figure3", 12.8, 5.2, dir = FIG_DIR)

  # ---- SI companions: MAIN-TEXT PASS ONLY ------------------------------------------------
  # Figures S4 and S2 are NOT drawn on an MI-pooled pass, and the reason is not cost.
  # S4A (front approach) is built from `zone_first`, the per-zone first confirmed onset of a
  # SINGLE load_linelist() draw, and S4B tracks ranks across CV rounds of that same draw.
  # Pairing either with an MI-pooled forecast would put a single-draw object and a
  # multiple-imputation mixture (M completions, run_info.json) in one figure and label the
  # result with one caption. S2's pooling semantics (a bivariate choropleth of pooled
  # probability x vulnerability, and a priority scatter) are not specified anywhere, and an
  # unspecified pooling rule is how a figure comes to mean something nobody chose. Both stay single-draw, in the canonical directory.
  if (MI_MODE) {
    message("[fig3] MI mode: Figures S4 and S2 skipped (single-draw objects; see the note here).")
    return(invisible(fig3))
  }

  # SI Figure S2: prob x vulnerability bivariate (A) + preparedness-priority scatter (B).
  # SAVED FIRST of the SI set: it needs neither the live-run date nor the line list, so a guard
  # that stops Figure S4 or Figure3_rank_surveillance below must not cost it.
  save_dual(p3c,  "FS2A_prob_vuln_bivariate", 6.4, 4.6)
  save_dual(p3d,  "FS2B_priority_scatter_h2", 6.4, 4.6)
  figS2 <- (wrap_elements(full = p3c) | wrap_elements(full = p3d)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 15.5 + 4.5, face = "bold", colour = INK))
  save_dual(figS2, "FigureS2_operational_prob_priority", 12.8, 5.2, dir = FIG_DIR)

  # --- SI narrative pair: front-approach (explanation) + rank evolution (response) ---
  # Paired on a shared time axis / colour / zone set so A mechanistically explains B.
  # Both reuse the province palette `pal` derived from the top-zone forest (main-text Fig 3B),
  # so the narrative panels share identical province hues with it.
  p3fa <- build_front_approach(base = 12, pal = pal)
  p3e  <- build_rank_evolution(base = 12, pal = pal)
  save_dual(p3fa, "FS4A_front_approach",      8.4, 5.2)   # SI Fig S4A
  save_dual(p3e,  "FS4B_rank_evolution",      8.4, 5.2)   # SI Fig S4B

  # SI Figure S4: front-approach (A) + risk-rank evolution (B) narrative pair. Identical province
  # palette -> patchwork collects ONE shared legend to the right.
  figS4 <- (p3fa | p3e) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") & theme(legend.position = "right")
  save_dual(figS4, "FigureS4_front_rank_evolution", 14.5, 5.6, dir = FIG_DIR)

  # Figure3_rank_surveillance: main-text Figure 3 A and B, unchanged, over (C) the rank
  # evolution of FigureS4B and (D) alerts and laboratory outcomes in the same zones. All four
  # panels are built from the objects above, so the zone set, province colours and live-run
  # ranks are those of Figure 3B by construction (build_zone_alerts() asserts it as well).
  # Panel C drops its province legend: panel B carries the same palette for the same zones.
  # Its axis names the horizon, which FigureS4B leaves to its caption: in this figure the
  # panel sits under B's "2-week forecast" and must not read as some other horizon.
  # ONE TAG STYLE: the wrapped top row would otherwise take patchwork's plain default while
  # C and D take theme_pub()'s bold tag, so the letters differed in weight and size.
  # Single-draw, like Figure S4, so not reached on an MI pass.
  alerts <- build_zone_alerts(d, pal, base = 12)
  p3c_rank <- p3e + labs(y = sprintf("%d-week invasion-risk rank (1 = highest)", FRONT_H)) +
    theme(legend.position = "none")
  fig3rs <- ((wrap_elements(full = p3a) | wrap_elements(full = p3b)) /
             (p3c_rank | alerts$plot) + plot_layout(widths = c(1.45, 1))) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 12 + 4.5, face = "bold", colour = INK))
  save_dual(fig3rs, "Figure3_rank_surveillance", 13.6, 10.8, dir = FIG_DIR)
  readr::write_csv(alerts$data %>%
                     transmute(rank_live, health_zone, province, p_case_invasion,
                               alerts_total, n_confirmed, n_not_a_case, alerts_resolved,
                               window_start = OUTBREAK_START, window_end = ANALYSIS_DATE),
                   file.path(FIG_DIR, "Figure3_rank_surveillance_panelD.csv"))

  invisible(fig3)
}

# -----------------------------------------------------------------------------
run <- function(f, nm) tryCatch({ message("== ", nm, " =="); f(); TRUE },
                                error = function(e){ message("!! ", nm, " FAILED: ", conditionMessage(e)); FALSE })
# ONLY THE REQUESTED BUILDERS RUN, AND ONLY THOSE ARE REPORTED. A builder excluded by
# MANU_ONLY must not appear in the status line or the failure tally as FALSE: that would
# report a figure as broken when it was never asked for, and on an MI pass — where only
# Figures 2 and 3 are pooled — it would make every run look like four failures.
.BUILDERS <- list(
  Figure1    = list(f = build_fig1,             lab = "Figure 1"),
  FigureS1   = list(f = build_fig1_si,          lab = "Figure S1"),
  FigureS3   = list(f = build_fig_si_epi_share, lab = "Figure S3"),
  Figure2    = list(f = build_fig2,             lab = "Figure 2"),
  Figure3    = list(f = build_fig3,             lab = "Figure 3"))
# .want() already leaves out .SINGLE_DRAW_BUILDERS on an MI pass.
.selected <- names(.BUILDERS)[vapply(names(.BUILDERS), .want, logical(1))]
if (!length(.selected))
  stop("[manu] no builder selected (MANU_ONLY=", paste(MANU_ONLY, collapse = ","),
       if (MI_MODE) "; MI mode additionally excludes the single-draw builders" else "", ")",
       call. = FALSE)
.panels <- vapply(.selected, function(nm) run(.BUILDERS[[nm]]$f, .BUILDERS[[nm]]$lab), logical(1))
names(.panels) <- .selected
message(sprintf("\n[done] %s  (featured %s)  ->  %s",
                paste(sprintf("%s:%s", names(.panels), .panels), collapse = "  "),
                FEATURED, FIG_DIR))

# A FAILED PANEL MUST FAIL THE SCRIPT. `run()` above catches so that one broken panel does not
# cost the others, but the script then reached its end and exited 0 -- and run_all.R decides
# OK/FAIL purely on the child's exit status. A panel that threw was therefore reported as a
# successful figure suite while the PREVIOUS run's PDF stayed on disk, so the manuscript
# carried a silently stale figure. Report per-panel status first (above), then exit non-zero.
if (any(!.panels))
  stop(sum(!.panels), " of ", length(.panels), " panel(s) failed: ",
       paste(names(.panels)[!.panels], collapse = ", "),
       ". Their files on disk are from an earlier run -- see the messages above for each cause.",
       call. = FALSE)

}   # end of the script-run guard

