# =============================================================================
# make_publication_figures.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
# Main-text multi-panel figures (Nature/Science aesthetic), rebuilt from SAVED
# data (no model re-fitting) under ONE unified publication theme.
#
#   Figure 1  Outbreak to date   A) cumulative-case map with the epicentre cohort's short-trip
#                                   destination shares (national extent; the zoomed top-N
#                                   variant is built as a separate panel)
#                                B) onset time series: cumulative confirmed cases (top) and
#                                   health zones invaded (bottom)
#   Figure 2  Model performance  A) discrimination (per-date AUC-PR skill over all rounds, both
#                                   horizons; the selection reads the post-burn-in rounds)
#                                B) prioritisation (share of invasions caught vs K, 1 week)
#                                C) ranking accuracy (mean rank of the invaded zones)
#                                D) ranking evolution of the featured model's watch-list
#                                E) top-K predicted risks per forecast round, by outcome
#             Figure2_labelled   the same panels with reader-friendly model names, a
#                                decoding caption, and panel E at 2 weeks over two rows
#   Figure 3  Spatial forecast   A) relative invasion-risk map B) uncertainty (CrI width) map
#                                C) top-zone probability and 90% CrI, both horizons
#                                D) preparedness priority (risk x vulnerability) scatter
#                                E) probability x vulnerability bivariate map
#   Panels F_topk_folds_h1/_h2: every forecast round of Figure 2E's bar chart.
#
# Outputs -> outputs/key_outputs/figures/{Figure1,Figure2,Figure2_labelled,Figure3}.{pdf,png}
#            and outputs/key_outputs/figures/panels/<panel>.{pdf,png}
# EVERY FILE GOES THROUGH save_dual()'s FIGURE_KEEP / FIGURE_DROP GATE (00_config.R): only the
# stems retained there are written, and everything else above is built and then discarded
# (FIGURE_KEEP_ALL=1 writes all of it). Which of these are deliverables is set there, not here.
# Run:  Rscript make_publication_figures.R   (from spatiotemporal/)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(sf); library(scales)
  library(viridisLite); library(lubridate); library(forcats)
})
sf::sf_use_s2(FALSE)
options(dplyr.summarise.inform = FALSE)

# Anchor to the repo via here::here(), NOT normalizePath(".") — that assumed the working
# directory was the repo root, so running this from spatiotemporal/ (where every other
# script in the suite is run from) sent HERE one level too deep and every path escaped
# the repo. make_manuscript_figures.R was already fixed for the same fault.
HERE      <- here::here()

# The pipeline config and data layer are REQUIRED, not optional. Paths, the outbreak
# start, the epicentre zones, the delay resolver and load_linelist() all come from them,
# and the figure must show the same case set as the model. The old ".cfg_ok" fallback
# quietly redrew main-text panels on local literals when the config failed to load — a
# figure that disagrees with its model, with one message in the log.
source(file.path(HERE, "spatiotemporal", "00_config.R"))
source(file.path(HERE, "spatiotemporal", "01_data_prep.R"))
# The mobility layer's loaders, for Figure 1A's short-trip shares (load_short_trip_profile(),
# the reader the M1 kernel uses). 03 builds and WRITES every mobility matrix when its own
# script guard passes, and FORCE_RUN=1 passes it from a source() too — so FORCE_RUN is
# withheld for this one source(): reading the loaders must not rebuild outputs/mobility/.
local({
  .fr <- Sys.getenv("FORCE_RUN", unset = NA)
  Sys.unsetenv("FORCE_RUN")
  on.exit(if (!is.na(.fr)) Sys.setenv(FORCE_RUN = .fr), add = TRUE)
  source(file.path(HERE, "spatiotemporal", "03_mobility_matrices.R"))
})
# (effective_onset_sample_delay() is no longer required: load_linelist() now imputes a
# missing onset from donor intervals, not from the fitted delay, and nothing here calls it.)
stopifnot(exists("OUT_DIR"), exists("OUTBREAK_START"), exists("EPICENTRE_ZONES"),
          exists("WEEK_ANCHOR"), exists("FLOWMINDER_ST_TAGS"), exists("FLOWMINDER_ST_DATES"),
          exists("load_linelist", mode = "function"),
          exists("load_short_trip_profile", mode = "function"),
          exists("load_aliases", mode = "function"))
.cfg <- function(name, default) if (exists(name)) get(name) else default

OUT       <- .cfg("OUT_DIR", file.path(HERE, "outputs"))
# Probability-scale switch (forecast_scale.R). Default "recalibrated": the primary set.
# FORECAST_SCALE=raw re-runs the identical build onto a raw/ sibling, same basenames.
source(file.path(HERE, "spatiotemporal", "forecast_scale.R"))
FIG_DIR   <- fs_out_dir(file.path(OUT, "key_outputs", "figures"))
PANEL_DIR <- fs_out_dir(file.path(OUT, "key_outputs", "figures", "panels"))
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# OUTBREAK_START and ANALYSIS_DATE (the as-of date the panels run up to) come from the config
# and nowhere else. They used to be read through .cfg() with a TYPED fallback — a literal
# "2026-04-30" and Sys.Date() — so a config that failed to define either would have drawn
# Figure 1 against an invented date (the day the script happened to run) without a word.
# A date the figure cannot derive is fatal.
for (.d in c("OUTBREAK_START", "ANALYSIS_DATE")) {
  .v <- get0(.d, ifnotfound = NULL)
  if (is.null(.v) || length(.v) != 1L || !inherits(.v, "Date") || is.na(.v))
    stop("[pubfig] ", .d, " is not a Date defined by 00_config.R; no date is substituted.",
         call. = FALSE)
}
rm(.d, .v)
SHP_PATH <- .cfg("SHAPEFILE_PATH",
                 file.path(HERE, "..", "data", "shapefiles", "DRC_Health_zones.shp"))
PROV_INT  <- .cfg("PROVINCES_OF_INTEREST", c("Ituri", "Nord-Kivu", "Haut-Uele"))
# Canonical (post-2026-07 shapefile) spellings in the fallback too: "Mongbalu" is not
# in the zone spine, so a fallback carrying it silently drops that origin.
EPI_ZONES <- .cfg("EPICENTRE_ZONES", c("Bunia", "Mongbwalu", "Rwampara"))

# The line list is resolved and logged by load_linelist() (01_data_prep.R); this script
# no longer resolves a path of its own.

# -----------------------------------------------------------------------------
# 1. DESIGN SYSTEM  (colourblind-safe; sequential = single perceptual hue)
# -----------------------------------------------------------------------------
# House style (palette, province scale, themes, save_dual) from figure_style.R — the
# single definition for the whole suite.
source(file.path(HERE, "spatiotemporal", "figure_style.R"))

# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("make_publication_figures.R")) {
  message("[pub] sourced, not run: nothing is written. Use `Rscript make_publication_figures.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

# The horizon colours (HZ_COL) come from figure_style.R with the rest of the palette. This
# file used to redefine them here, and carried a SER_COL that no panel read.

# Reader-friendly model names for the labelled Figure 2 variant. These are GENERATED from the
# method code by model_pretty_label() (00_config.R), which composes the name from the same
# tokens bayes_default_grid() composes the label from — kernel family, road-km / source-fill /
# origin-split qualifiers, covariate set, beta_t process, generation-time arm.
#
# IT USED TO BE A 100-ENTRY HAND-WRITTEN LOOKUP, and it drifted: a sixth of the cross-validated
# models had no entry, so raw codes such as "Bayes-M8-dist-split-geo" were printed on the
# published figure whose entire purpose is to avoid them, and each new model family had to be
# transcribed here by hand. A derived name cannot fall behind the grid.
# There is no table: build_fig2(model_labels = TRUE) calls model_pretty_label() on whatever
# codes the panel actually holds.

#' Thin a set of dated tick positions to at most `max_n` labels, always keeping the FIRST and
#' the LAST.
#'
#' Dated forecast rounds accumulate as the outbreak runs, so any fixed rule ("every other
#' round") eventually overprints. This picks evenly spaced positions across the sequence and
#' pins the endpoints, which are the two the reader needs: the first round and the round the
#' end-of-line labels annotate. Returns the dates themselves, so it can be passed straight to
#' scale_x_date(breaks = ).
.thin_date_breaks <- function(x, max_n = 5L) {
  cs <- sort(unique(as.Date(x)))
  if (length(cs) <= max_n) return(cs)
  idx <- unique(c(1L, round(seq(1, length(cs), length.out = max_n)), length(cs)))
  cs[sort(unique(idx))]
}


# (save_dual(), which writes the vector PDF + 600-dpi PNG, is in figure_style.R.)

# -----------------------------------------------------------------------------
# 2. LOAD SAVED DATA
# -----------------------------------------------------------------------------
message("[load] shapefile, risk scores, LFO results, evaluation, linelist ...")
shp <- st_read(SHP_PATH, quiet = TRUE) %>%
  mutate(.key = tolower(trimws(Nom)), .prov = as.character(PROVINCE))

rs  <- read_csv(fs_risk_csv(OUT), show_col_types = FALSE)
lfo <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))
ev  <- read_csv(file.path(OUT, "diagnostics", "invasion_evaluation.csv"),
                show_col_types = FALSE)
# Put the selected scale into `p_invasion` ONCE, keeping the untouched values as
# `p_invasion_raw`. Every expression below reads `p_invasion`, so the switch is total
# rather than depending on each call site being edited consistently.
# The column actually drawn is kept: the per-round panels label any round that carries no
# fitted recalibration factor, and only on a pass that draws the recalibrated column.
PSCALE_COL <- fs_apply_lfo_scale(lfo)
# NOTE on the discrimination columns used below (auc_pr_skill, auc_roc): they are read from
# the RAW evaluation columns in BOTH passes, deliberately. Within a fold the recalibration
# is strictly monotone, so AUC is EXACTLY invariant; the `*_recal` columns differ only
# because POOLED AUC mixes folds carrying different prequential factors. Publishing a
# "recalibrated AUC" that differs from the raw one would invite the reader to conclude that
# recalibration changed discrimination, which is false. Brier, log score and reliability DO
# change with the scale and are computed from the selected probabilities wherever they
# appear. (make_manuscript_figures.R names the same choice as EV_SKILL/EV_ROC because it
# quotes those numbers in a panel annotation; this suite does not, so there is nothing to
# parameterise here.)

# The featured Bayesian model is the one the PIPELINE selected. It is read from
# model_selection.json rather than recomputed here (read_featured_model() below), so the figures
# cannot track a different model from the pipeline, and the log-score axis, the eligibility
# gates and the tie-break are whatever best_invasion_model() actually applied.
# read_featured_model() (figure_style.R) is the reader every suite shares, with NO FALLBACK:
# this file used to fall back to a typed "Bayes-M14-fill-med" (or the first Bayesian method in
# the LFO) when model_selection.json was unreadable, which would label a figure with a model
# the pipeline may not have chosen.
BEST_BAYES <- read_featured_model(file.path(OUT_DIR, "key_outputs", "model_selection.json"))
# Fold / event counts for captions, computed from the actual LFO rather than fixed.
LFO_N_FOLDS  <- dplyr::n_distinct(lfo$fold_id)
# Count invasion events for ONE model only: is_new_invasion is the observed outcome, identical
# across every model in the grid, so summing over all methods would multiply the true count by the
# number of models. Restrict to the featured model (the standard at-risk row set).
LFO_N_EVENTS <- lfo %>% dplyr::filter(horizon == 1, method == BEST_BAYES) %>%
  dplyr::summarise(n = sum(is_new_invasion, na.rm = TRUE)) %>% dplyr::pull(n)
message(sprintf("[figures] Featured Bayesian model: %s; LFO %d folds, %d h=1 events.",
                BEST_BAYES, LFO_N_FOLDS, LFO_N_EVENTS))

# -----------------------------------------------------------------------------
# PRIORITISATION (detection) curve — READ, not rebuilt
# -----------------------------------------------------------------------------
# key_outputs/detection_curves.csv is published by run_all.R from compute_detection_curve()
# (20_forecast_detail.R), the estimator behind the recall_at_<K> columns of
# invasion_evaluation.csv, on the shared support. Figure 2 panel A of the manuscript suite
# draws the same file.
#
# THIS PANEL USED TO REBUILD THE CURVE ITSELF, and so described a different object from the
# manuscript's Figure 2A under the same name:
#   * it POOLED recall over folds (sum of hits / sum of invasions) where the published curve
#     AVERAGES it per fold, so "share of true invasions caught" was a different estimand;
#   * its one comparator was a "Naive-epicentre-inflow" built from MOBILITY_PRIMARY (the
#     M14-fill kernel's own epicentre rows), not the Baseline-flowminder-inflow ranking the
#     manuscript publishes, and it omitted the other FIGURE_BASELINE_METHODS altogether;
#   * it drew no interval.
# It now shows the featured model against the published comparator set, from the published
# curve, with the published fold-cluster bootstrap band.
DETCURVE_SRC <- file.path(OUT, "key_outputs", "detection_curves.csv")
BASELINE_METHODS <- fs_baseline_methods()

#' The featured model's and the comparators' curves at one horizon, with the random reference.
#' A missing series is fatal — a thinner panel would make a weaker claim under the same name.
prioritisation_curves <- function(horizon = 1L) {
  if (!file.exists(DETCURVE_SRC))
    stop("[fig2] ", basename(DETCURVE_SRC), " not found; run_all.R publishes it. Run the ",
         "pipeline before this script.", call. = FALSE)
  d <- suppressWarnings(read_csv(DETCURVE_SRC, show_col_types = FALSE))
  need <- c("method", "horizon", "k", "recall", "recall_lo", "recall_hi", "recall_random")
  if (length(setdiff(need, names(d))))
    stop("[fig2] ", basename(DETCURVE_SRC), " lacks column(s) ",
         paste(setdiff(need, names(d)), collapse = ", "), call. = FALSE)
  want <- c(BEST_BAYES, BASELINE_METHODS)
  cu <- d %>% filter(horizon == !!horizon, method %in% want)
  miss <- setdiff(want, unique(cu$method))
  if (length(miss))
    stop("[fig2] ", basename(DETCURVE_SRC), " has no h=", horizon, " curve for ",
         paste(miss, collapse = ", "), "; re-run run_all.R.", call. = FALSE)
  # recall_random is a property of the fold design, not of a method: one value per k.
  rnd <- cu %>% distinct(k, recall_random) %>% arrange(k)
  if (nrow(rnd) != dplyr::n_distinct(cu$k))
    stop("[fig2] recall_random differs between methods at h=", horizon, "; the curves were ",
         "scored on different supports.", call. = FALSE)
  list(curves = cu, rnd = rnd)
}

# (The zone-table -> shapefile join, join_map(), is in figure_style.R.)

# -----------------------------------------------------------------------------
# 3. FIGURE 1 DATA  (onset-dated confirmed cases; pipeline-consistent)
# -----------------------------------------------------------------------------
# The line list is the PIPELINE's. load_linelist() resolves the latest processed export,
# reconciles it to the INSP sitrep cumulative (appending the shortfall as confirmed rows
# for zones the line list has not caught up on; these top-up rows carry no onset by
# construction), and imputes a missing onset as the sample date minus a sample-to-onset
# interval DRAWN FROM DONORS: records of the same class, sampled in the same week, whose
# onset is recorded (.impute_onset_delays(), 01_data_prep.R). Each record's draw is a uniform
# hashed from its own alert_id and ONSET_IMPUTE_SEED (which defaults to RANDOM_SEED), not a
# stream from an RNG state. The caption quotes the loader's own description of that method
# (.ONSET_IMPUTE_DESC). This figure therefore shows exactly the case set, zone set and onset
# dates the invasion model is fitted and scored on.
#
# This script used to read the processed CSV directly and shift every missing onset back
# by a single ROUNDED MEAN delay. A fixed shift is a different estimator from the
# pipeline's, so the panel's arrival dates and cumulative curve matched no table in the
# pipeline, and the caption had to carry a standing disclaimer saying so. Both are now gone.
ll <- load_linelist()
# Figure 1B's time axis starts on the first day of the zone-week grid's first week: the
# WEEK_ANCHOR floor of OUTBREAK_START, as invasion_week_index() (00_config.R) numbers it, so
# fs_date_breaks()'s whole-week ticks fall on the grid's week starts. It used to floor to a
# MONDAY, which is neither the outbreak start nor a week start of the grid.
wk_floor <- floor_date(OUTBREAK_START, "week", week_start = WEEK_ANCHOR)
stopifnot(all(c("health_zone", "province", "date_index", "confirmed", "onset_usable") %in% names(ll)))

# No lower date bound here: load_linelist() never puts date_index before OUTBREAK_START (it
# floors it there), so the plotted set is every confirmed, dated, zoned record it returns.
conf <- ll %>% filter(confirmed %in% TRUE, !is.na(date_index), !is.na(health_zone))
stopifnot(min(conf$date_index) >= OUTBREAK_START)
# THE IMPUTED SHARE IS QUOTED BY PROVENANCE, from load_linelist()'s own published globals.
# The share used to be computed here over every plotted confirmed row, which INCLUDES the
# sitrep top-up rows: those carry no onset by construction (they are a count reconciliation,
# not case records), so they are 100% "imputed" and inflated a figure that reads as the line
# list's onset missingness. The caption now states the two separately: the share of confirmed
# LINE-LIST cases whose onset was imputed (.PCT_ONSET_IMPUTED_LINELIST) and the number of
# sitrep top-up rows (.N_SITREP_APPENDED), whose onsets are imputed by the same method.
# The method itself is stated in the loader's own words (.ONSET_IMPUTE_DESC), so the caption
# cannot describe a different one. All three are published together by load_linelist(); a
# missing one means a stale loader, and no number is substituted for it.
.onset_glob <- function(nm) {
  v <- get0(nm, envir = .GlobalEnv, ifnotfound = NULL)
  if (is.null(v) || length(v) != 1L || is.na(v))
    stop("[pubfig] load_linelist() did not publish ", nm, "; the Figure 1B caption cannot ",
         "state the onset imputation without it, and no value is substituted.", call. = FALSE)
  v
}
.ONSET_IMPUTE_TEXT   <- .onset_glob(".ONSET_IMPUTE_DESC")
PCT_IMPUTED_LINELIST <- as.numeric(.onset_glob(".PCT_ONSET_IMPUTED_LINELIST"))
N_SITREP_APPENDED    <- as.integer(.onset_glob(".N_SITREP_APPENDED"))
# Under ONSET_MODE = "complete_case" nothing is imputed (onset-less records are dropped), and
# a caption saying "imputed for x%" would be false; the mode is published with the method.
.ONSET_MODE <- .onset_glob(".ONSET_IMPUTE_MODE")
# THE SITREP TOP-UP ROWS HAVE NO SAMPLE DATE EITHER. .build_sitrep_confirmed_appends()
# (01_data_prep.R) dates each one by its SITREP REPORT DATE, less SITREP_APPEND_SHIFT_DAYS, and
# that date stands in for the sample date in the imputation. The loader's method text says
# "the sample date minus ...", which is true of the line-list records only, so the top-up rows
# get a clause of their own naming the date they are actually imputed from. The shift is read
# by the same expression the loader reads it with.
.SITREP_SHIFT <- as.integer(get0("SITREP_APPEND_SHIFT_DAYS", ifnotfound = 0L))
F1B_ONSET_CAPTION <- if (identical(.ONSET_MODE, "complete_case")) {
  sprintf("Onset dates: %s.", .ONSET_IMPUTE_TEXT)
} else {
  paste0(
    sprintf("Onset imputed for the %.1f%% of confirmed line-list cases lacking a usable onset: %s.",
            PCT_IMPUTED_LINELIST, .ONSET_IMPUTE_TEXT),
    if (N_SITREP_APPENDED > 0L)
      sprintf(paste0(" The %s confirmed case%s added to reconcile zone counts with the INSP ",
                     "sitrep carr%s no onset or sample date; %s onset%s imputed the same way, ",
                     "from %s in place of the sample date."),
              format(N_SITREP_APPENDED, big.mark = ","),
              if (N_SITREP_APPENDED == 1L) "" else "s",
              if (N_SITREP_APPENDED == 1L) "ies" else "y",
              if (N_SITREP_APPENDED == 1L) "its" else "their",
              if (N_SITREP_APPENDED == 1L) " is" else "s are",
              if (.SITREP_SHIFT > 0L)
                sprintf("%d day%s before the sitrep report date", .SITREP_SHIFT,
                        if (.SITREP_SHIFT == 1L) "" else "s")
              else "the sitrep report date")
    else "")
}

# 1A cumulative confirmed per zone
zone_cum <- conf %>% group_by(health_zone, province) %>%
  summarise(cases = n(), .groups = "drop")
# 1B national cumulative + affected-zone count over onset time
daily <- conf %>% count(date_index, name = "n") %>% arrange(date_index) %>%
  mutate(cum_cases = cumsum(n))
zone_first <- conf %>% group_by(health_zone) %>%
  summarise(first = min(date_index), .groups = "drop") %>% arrange(first) %>%
  mutate(n_zones = row_number())
message(sprintf("[fig1] %d confirmed cases, %d zones, onset %s to %s",
                nrow(conf), nrow(zone_cum), min(daily$date_index), max(daily$date_index)))

# --- 1A mobility overlay: the epicentre cohort's Flowminder short-trip destination shares ---
# WHAT THE NUMBERS ARE (data/flowminder_short_trips/metadata.yaml, Annex A of Flowminder's
# short-trip report): of the subscribers seen in the epicentre zones during 3-23 April, pooled
# into ONE cohort (EPICENTRE_ZONES), the percentage later observed in each destination zone by
# the snapshot date. A subscriber can be observed in more than one zone, so the shares sum to
# more than 100%: they are shares of the COHORT, not shares of trips. The cohort has no
# per-origin split, so the fan is anchored at the first epicentre zone and labelled as the
# cohort, not as that one zone. The strongest N shares are drawn as width-encoded,
# DIRECTIONAL arcs (arrowheads point epicentre -> destination).
#
# READ THROUGH load_short_trip_profile() (03_mobility_matrices.R) — the loader build_M1()
# uses, at the snapshot build_M1() picks (the latest on or before the analysis date) — which
# resolves the RAW annex table against the zone spine. This panel used to read the PROCESSED
# outflow CSV, which upstream canonicalisation against the 437-zone OD header had stripped of
# destinations the raw table resolves (Kilo, an invaded zone and one of the cohort's largest
# destinations, among them), so it drew a fan that the mobility layer does not use.
#
# A FAILED LOAD IS RECORDED, NOT SWALLOWED. It used to sit in try(silent = TRUE), after which
# Figure 1A was published without its mobility overlay and one log line said so. The error is
# kept in FLOW_ERR and raised by the panels that need the shares, so Figure 1 fails (and the
# script with it, below) while the other figures are still built.
FLOW_ST_TAG <- local({
  .elig <- which(FLOWMINDER_ST_DATES <= ANALYSIS_DATE)
  if (length(.elig)) FLOWMINDER_ST_TAGS[.elig[which.max(FLOWMINDER_ST_DATES[.elig])]]
  else FLOWMINDER_ST_TAGS[which.min(FLOWMINDER_ST_DATES)]   # build_M1()'s rule when all post-date
})
EPI_ORIGIN  <- EPI_ZONES[1]     # anchor for the fan (the cohort is pooled over EPI_ZONES)
EPI_LABEL   <- sprintf("Epicentre cohort\n(%s)", paste(EPI_ZONES, collapse = "/"))
ITURI_PROV  <- PROV_INT[1]      # epicentre province (Ituri) — for the "share staying within" stat
N_FLOWS     <- 12L              # strongest N destinations (keeps the map legible)
FLOW_COL    <- "#EE7733"        # Tol 'vibrant' orange; warm accent vs the cool mako fill
# Kept to the width of the legend it replaced: a wider guide shrinks the national map.
FLOW_LEGEND <- "% of epicentre\ncohort seen there"
flows_df <- NULL; flows_all <- NULL; epi_pt <- NULL; latest_flow_date <- NA; ITURI_PCT <- NA_real_
FLOW_ERR <- NULL
tryCatch({
  prof <- load_short_trip_profile(FLOW_ST_TAG, zones_all = as.character(shp$Nom),
                                  aliases = load_aliases())
  # Percent, as the annex states it and as the legend reads. The cohort's own zones are not
  # destinations of it, so any that appear are not drawn.
  st <- data.frame(nom = names(prof$props), prop = 100 * unname(prof$props),
                   stringsAsFactors = FALSE)
  st <- st[!st$nom %in% EPI_ZONES & is.finite(st$prop) & st$prop > 0, , drop = FALSE]
  if (!nrow(st)) stop("snapshot ", FLOW_ST_TAG, " has no destination share")
  # zone centroids (point-on-surface = always inside the polygon) + province, for arc geometry
  # and the destination-province breakdown. The profile's names WERE resolved against shp$Nom,
  # and the zone names are unique, so each destination is exactly one polygon; anything else
  # is an error, not a nearest-instance guess.
  pts <- suppressWarnings(sf::st_point_on_surface(sf::st_geometry(shp)))
  cc  <- sf::st_coordinates(pts)
  cent <- data.frame(nom = as.character(shp$Nom), key = shp$.key, prov = as.character(shp$.prov),
                     cx = cc[, 1], cy = cc[, 2], stringsAsFactors = FALSE)
  if (anyDuplicated(cent$nom) || !all(c(st$nom, EPI_ORIGIN) %in% cent$nom))
    stop("the short-trip destinations do not map one-to-one onto the shapefile's zones")
  epi_row <- cent[cent$nom == EPI_ORIGIN, , drop = FALSE]
  fl_all <- merge(st, cent, by = "nom")
  fl_all$x <- epi_row$cx; fl_all$y <- epi_row$cy; fl_all$xend <- fl_all$cx; fl_all$yend <- fl_all$cy
  fl_all <- fl_all[order(fl_all$prop), ]            # small first so big arcs draw on top
  fl <- head(fl_all[order(-fl_all$prop), ], N_FLOWS)
  flows_df <- fl; flows_all <- fl_all; epi_pt <- epi_row
  latest_flow_date <- as.Date(FLOW_ST_TAG, "%Y%m%d")
  # COMPUTED (not hard-coded): the share of the cohort's summed destination shares that lies in
  # Ituri (a share of the sum, since the shares themselves sum to more than 100%).
  ITURI_PCT <- round(100 * sum(fl_all$prop[fl_all$prov == ITURI_PROV], na.rm = TRUE) /
                     sum(fl_all$prop, na.rm = TRUE), 0)
  message(sprintf(paste0("[fig1] epicentre-cohort shares: %d destinations (%s table), snapshot %s; ",
                         "%.0f%% of the summed share in %s; top: %s"),
                  nrow(fl_all), prof$source, latest_flow_date, ITURI_PCT, ITURI_PROV,
                  paste(head(fl$nom, 4), collapse = ", ")))
}, error = function(e) FLOW_ERR <<- conditionMessage(e))
if (!is.null(FLOW_ERR)) message("[fig1] ERROR: the epicentre short-trip shares did not load (",
                                FLOW_ERR, "); Figure 1 will fail.")
.need_flows <- function() if (!is.null(FLOW_ERR))
  stop("[fig1] the epicentre short-trip shares did not load (", FLOW_ERR, "); Figure 1A is ",
       "not drawn without its mobility overlay.", call. = FALSE)

# F1A zoom window: the affected zones + epicentre + flow destinations sit in a small
# NE cluster, so we frame the panel on that region (+ padding) and add a national
# locator inset for context, rather than showing ~500 empty zones at national scale.
.foc_keys <- unique(c(tolower(trimws(zone_cum$health_zone)), tolower(EPI_ORIGIN),
                      if (!is.null(flows_df)) flows_df$key))
.foc <- shp[shp$.key %in% .foc_keys, ]
.bb  <- sf::st_bbox(.foc)
.px  <- 0.14 * as.numeric(.bb["xmax"] - .bb["xmin"])
.py  <- 0.14 * as.numeric(.bb["ymax"] - .bb["ymin"])
ZOOM_X <- c(as.numeric(.bb["xmin"]) - .px, as.numeric(.bb["xmax"]) + .px)
ZOOM_Y <- c(as.numeric(.bb["ymin"]) - .py, as.numeric(.bb["ymax"]) + .py)

# =============================================================================
# FIGURE 1
# =============================================================================
build_fig1 <- function() {
  .need_flows()
  # --- 1A: cumulative-case choropleth (binned sequential) ---
  mp <- join_map(zone_cum) %>%
    mutate(cases = ifelse(is.na(cases), 0, cases),
           bin = cut(cases, breaks = c(-1, 0, 9, 49, 99, 499, Inf),
                     labels = c("0","1-9","10-49","50-99","100-499","500+")))
  seq6 <- viridisLite::mako(6, begin = 0.92, end = 0.12)      # light -> dark
  names(seq6) <- c("0","1-9","10-49","50-99","100-499","500+")
  seq6["0"] <- NA_FILL
  # epicentre labels for the few biggest zones
  lab <- mp %>% filter(cases >= 100)
  lab_pts <- suppressWarnings(st_point_on_surface(st_geometry(lab)))
  lab_df  <- cbind(st_drop_geometry(lab)["health_zone"], st_coordinates(lab_pts))
  # unified label set so nothing overlaps or duplicates: high-case zones (bold, ink)
  # + top flow destinations NOT already labelled (plain, muted) — one repel call.
  case_lab <- data.frame(x = lab_df$X, y = lab_df$Y, label = lab_df$health_zone,
                         face = "bold", col = INK, stringsAsFactors = FALSE)
  lab_all <- case_lab
  if (!is.null(flows_df) && nrow(flows_df) > 0) {
    fd <- head(flows_df[!(flows_df$key %in% tolower(lab_df$health_zone)), , drop = FALSE], 6)
    if (nrow(fd) > 0)
      lab_all <- rbind(case_lab,
        data.frame(x = fd$xend, y = fd$yend, label = fd$nom,
                   face = "plain", col = MUTED, stringsAsFactors = FALSE))
  }

  p1a <- ggplot(mp) +
    geom_sf(aes(fill = bin), colour = "white", linewidth = 0.08) +
    # --- epicentre outflow fan (drawn under the labels, over the choropleth) ---
    { if (!is.null(flows_df) && nrow(flows_df) > 0) list(
        geom_curve(data = flows_df,
                   aes(x = x, y = y, xend = xend, yend = yend, linewidth = prop),
                   curvature = 0.16, angle = 90, ncp = 16, colour = FLOW_COL,
                   alpha = 0.82, lineend = "round",
                   arrow = grid::arrow(length = unit(0.04, "in"), type = "closed",
                                       angle = 20)),
        geom_point(data = epi_pt, aes(cx, cy), shape = 21, size = 3.6,
                   fill = "white", colour = FLOW_COL, stroke = 1.3)
      ) } +
    { if (nrow(lab_all) > 0) ggrepel::geom_text_repel(
        data = lab_all, aes(x, y, label = label, fontface = face, colour = col),
        size = 3.3, min.segment.length = 0, segment.colour = FAINT,
        segment.size = 0.25, box.padding = 0.5, point.padding = 0.3,
        max.overlaps = 40, seed = 1) } +
    scale_colour_identity() +
    scale_fill_manual(values = seq6, na.value = NA_FILL, drop = FALSE,
                      name = "Cumulative\nconfirmed cases",
                      guide = guide_legend(reverse = TRUE, order = 1,
                                           override.aes = list(colour="white"))) +
    scale_linewidth_continuous(
      name = FLOW_LEGEND, range = c(0.3, 2.6),
      breaks = c(5, 10, 20), limits = c(0, NA),
      guide = guide_legend(order = 2, override.aes = list(colour = FLOW_COL, alpha = 0.9))) +
    coord_sf(xlim = ZOOM_X, ylim = ZOOM_Y, expand = FALSE) +
    theme_map(11.5) +
    theme(legend.position = "right",
          plot.title    = element_blank(),   # no embedded titles: the caption is in the manuscript
          plot.subtitle = element_blank(),
          plot.caption  = element_blank(),
          panel.border  = element_rect(fill = NA, colour = GRID, linewidth = 0.4))

  # national locator inset: whole DRC in grey with the zoom window boxed
  locator <- ggplot(shp) +
    geom_sf(fill = "grey86", colour = "white", linewidth = 0.04) +
    annotate("rect", xmin = ZOOM_X[1], xmax = ZOOM_X[2],
             ymin = ZOOM_Y[1], ymax = ZOOM_Y[2],
             fill = NA, colour = FLOW_COL, linewidth = 0.55) +
    coord_sf(expand = FALSE) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = FAINT, linewidth = 0.3),
          plot.margin = margin(2, 2, 2, 2))
  p1a <- p1a + patchwork::inset_element(
    locator, left = 0.00, bottom = 0.00, right = 0.30, top = 0.30,
    align_to = "panel", clip = FALSE)

  # --- 1B: onset time series (two stacked, shared x) ---
  # Run the x-axis (and the cumulative curves) up to the analysis date: append a terminal
  # point carrying the final cumulative value forward so both step curves plateau to the
  # as-of date rather than stopping at the last onset event.
  xlim <- c(wk_floor, ANALYSIS_DATE)
  daily_x <- daily
  if (max(daily_x$date_index) < ANALYSIS_DATE)
    daily_x <- bind_rows(daily_x, tibble(date_index = ANALYSIS_DATE, n = 0L,
                                         cum_cases = max(daily_x$cum_cases)))
  zone_first_x <- zone_first
  if (max(zone_first_x$first) < ANALYSIS_DATE)
    zone_first_x <- bind_rows(zone_first_x, tibble(health_zone = NA_character_,
                                                   first = ANALYSIS_DATE,
                                                   n_zones = max(zone_first_x$n_zones)))
  # Breaks from fs_date_breaks() (figure_style.R): a fixed two-week step put twelve "%d %b"
  # labels on this ~5 in axis at 13.5 pt and they overprinted. Five is what it carries.
  xsc  <- scale_x_date(limits = xlim, breaks = fs_date_breaks(xlim, max_n = 5L),
                       date_labels = "%d %b", expand = expansion(mult = c(0.01, 0.02)))
  top <- ggplot(daily_x, aes(date_index, cum_cases)) +
    geom_area(fill = "#0072B2", alpha = 0.16) +
    geom_step(colour = "#0072B2", linewidth = 0.85, direction = "hv") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06)), labels = comma) +
    xsc +
    labs(title = "Epidemic trajectory by symptom-onset date",
         y = "Cumulative\nconfirmed cases", x = NULL) +
    theme_pub(16.5) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                        plot.margin = margin(6,8,0,6))
  bot <- ggplot(zone_first_x, aes(first, n_zones)) +
    geom_step(colour = "#D55E00", linewidth = 0.95, direction = "hv") +
    geom_point(data = zone_first, aes(first, n_zones), colour = "#D55E00", size = 1.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                       limits = c(0, NA)) + xsc +
    # "Health zones / invaded", not "Health zones with / >=1 confirmed case": the longer title
    # ran up the shared left margin into the top panel's "Cumulative / confirmed cases", and
    # the two read as one run of text. "Invaded" = at least one confirmed case, as everywhere
    # else in the suite.
    labs(y = "Health zones\ninvaded",
         x = "Symptom-onset date", caption =
         # The panel and the pipeline now share one estimator, so the caption states it once.
         # Wrapped at 90 characters and placed from the PLOT's left edge: at width 110 and the
         # default panel-anchored position, the first line ran ~0.5 in past the right edge of
         # the panel in Figure 1 and was clipped mid-word.
         # The text (line-list share and sitrep top-up rows stated separately) is built once,
         # beside the data it describes: F1B_ONSET_CAPTION, section 3.
         stringr::str_wrap(F1B_ONSET_CAPTION, width = 90)) +
    # theme_pub() blanks plot.caption, so this disclosure was COMPUTED AND DISCARDED — a main-text
    # figure in which a quarter of the plotted onsets are imputed said nothing about it. Re-enable
    # the caption for this panel only.
    theme_pub(16.5) +
    # margin(t): the caption grew to four lines when the sitrep top-up rows were stated
    # separately, and sat flush against the x-axis title.
    theme(plot.caption = element_text(size = 7.5, colour = MUTED, hjust = 0,
                                      margin = margin(t = 5)),
          plot.caption.position = "plot") +
    theme(plot.margin = margin(0,8,6,6),
                            axis.text.x = element_text(size = 13.5))
  p1b <- top / bot + plot_layout(heights = c(1, 0.78))

  save_dual(p1a, "F1A_cumulative_case_map", 6.4, 5.4)
  save_dual(p1b, "F1B_onset_timeseries",    6.2, 5.4)

  # Figure 1 panel A is the national all-destinations case map (the cohort's full reach). There
  # is no flow-less fallback: the shares are required (.need_flows() above).
  p1a_main <- make_f1a_allflows_panel()

  # No embedded figure title/subtitle — only the A/B panel tags. The descriptive text is the
  # figure caption in the manuscript, not part of the figure.
  fig1 <- (wrap_elements(full = p1a_main) | wrap_elements(full = p1b)) +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(tag_levels = "A")
  save_dual(fig1, "Figure1", 11.8, 5.0, dir = FIG_DIR)
  invisible(fig1)
}

# =============================================================================
# FIGURE 1A — ALL epicentre-cohort destinations (national extent)
# Companion to the zoomed top-N panel: shows EVERY destination of the epicentre cohort's short-trip
# shares with DIRECTIONAL arcs (arrowheads point epicentre -> destination). Most of the summed share
# lies within the epicentre province (the dense near fan; ITURI_PCT, computed and logged above); the
# thin threads are the long tail reaching Kinshasa and the far provinces. Text is enlarged for
# main-figure legibility.
# =============================================================================
make_f1a_allflows_panel <- function() {
  .need_flows()
  mp <- join_map(zone_cum) %>%
    mutate(cases = ifelse(is.na(cases), 0, cases),
           bin = cut(cases, breaks = c(-1, 0, 9, 49, 99, 499, Inf),
                     labels = c("0","1-9","10-49","50-99","100-499","500+")))
  seq6 <- viridisLite::mako(6, begin = 0.92, end = 0.12)
  names(seq6) <- c("0","1-9","10-49","50-99","100-499","500+"); seq6["0"] <- NA_FILL

  # label only the epicentre cohort (destination place-names dropped for legibility; the
  # far-western Kinshasa terminus is still annotated separately below). The cohort is pooled over
  # EPI_ZONES, so the anchor is labelled as the cohort, not as the one zone it is drawn from.
  lab_all <- data.frame(x = epi_pt$cx, y = epi_pt$cy, label = EPI_LABEL,
                        face = "bold", col = INK, stringsAsFactors = FALSE)
  # KINSHASA IS DRAWN AS ONE ARC. Its zones sit close together at the far end of the country, so
  # one arc per zone drew near-coincident threads that overplotted into a band far wider than any
  # one of their shares — and the eye reads band width as flow. They are merged into ONE arc
  # carrying their SUMMED share, ending at the strongest Kinshasa destination's own centroid
  # (data-derived; no hard-coded coordinates), which is also where the "Kinshasa" label sits.
  # flows_all itself keeps one row per zone, so ITURI_PCT and the zoomed panel are unaffected.
  arcs <- flows_all
  .kin_rows <- !is.na(arcs$prov) & arcs$prov == "Kinshasa"
  kin <- NULL
  if (any(.kin_rows)) {
    kz  <- arcs[.kin_rows, , drop = FALSE]
    kin <- kz[which.max(kz$prop), , drop = FALSE]
    kin$prop <- sum(kz$prop); kin$nom <- "Kinshasa"
    arcs <- rbind(arcs[!.kin_rows, , drop = FALSE], kin)
    arcs <- arcs[order(arcs$prop), ]                 # small first so big arcs draw on top
  }

  arr <- grid::arrow(length = unit(0.05, "in"), type = "closed", angle = 18)
  # This panel blanks plot.title, plot.subtitle and plot.caption (below): the captions live in
  # the manuscript, outside the figure. ITURI_PCT (the share of the cohort's summed destination
  # shares that lies in Ituri) is COMPUTED and LOGGED where it is derived, for the caption text.
  p <- ggplot(mp) +
    geom_sf(aes(fill = bin), colour = "white", linewidth = 0.06) +
    # every destination arc (Kinshasa merged, above); width = the cohort's share, DIRECTIONAL
    # (arrow epicentre -> destination)
    geom_curve(data = arcs,
               aes(x = x, y = y, xend = xend, yend = yend, linewidth = prop),
               curvature = 0.16, angle = 90, ncp = 12, colour = FLOW_COL,
               alpha = 0.6, lineend = "round", arrow = arr) +
    geom_point(data = epi_pt, aes(cx, cy), shape = 21, size = 3.6,
               fill = "white", colour = FLOW_COL, stroke = 1.3) +
    ggrepel::geom_text_repel(
      data = lab_all, aes(x, y, label = label, fontface = face, colour = col),
      # 4.2, not the 5.2 the one-line "Bunia (epicentre)" label had: the cohort label runs to
      # two lines and at 5.2 covered most of the invaded zones around the epicentre.
      size = 4.2, lineheight = 0.9, min.segment.length = 0, segment.colour = FAINT,
      segment.size = 0.25, box.padding = 0.5, max.overlaps = 40, seed = 1) +
    { if (!is.null(kin)) annotate("text", x = kin$xend, y = kin$yend - 0.35,
             label = "Kinshasa", size = 4.8, colour = MUTED, fontface = "italic") } +
    scale_colour_identity() +
    scale_fill_manual(values = seq6, na.value = NA_FILL, drop = FALSE,
                      name = "Cumulative\nconfirmed cases",
                      guide = guide_legend(reverse = TRUE, order = 1,
                                           override.aes = list(colour = "white"))) +
    scale_linewidth_continuous(
      name = FLOW_LEGEND, range = c(0.15, 2.9),
      breaks = c(1, 5, 10, 20), limits = c(0, NA),
      guide = guide_legend(order = 2, override.aes = list(colour = FLOW_COL, alpha = 0.9))) +
    coord_sf(expand = FALSE) +
    theme_map(17) +
    theme(legend.position = "right",
          plot.title    = element_blank(),   # no embedded titles: the caption is in the manuscript
          plot.subtitle = element_blank(),
          plot.caption  = element_blank(),
          panel.border  = element_rect(fill = NA, colour = GRID, linewidth = 0.4))

  p
}

# thin wrapper: build the national all-flows panel and save it standalone
build_f1a_allflows <- function() {
  p <- make_f1a_allflows_panel()     # stops, rather than returning nothing, if the shares failed
  save_dual(p, "F1A_cumulative_case_map_allflows", 6.6, 6.2)
  invisible(p)
}

# =============================================================================
# FIGURE 2  (evaluation; from lfo_results + invasion_evaluation.csv)
# =============================================================================
build_fig2 <- function(model_labels = NULL, file_suffix = "", fig_w = 11.4, fig_h = 13.3,
                       topk_horizon = 1L, topk_nrow = 1L) {
  # Panel E is the per-fold top-K bar chart (build_topk_folds). topk_horizon selects the
  # forecast lead (1 or 2 weeks); topk_nrow wraps the fold facets over that many rows so a
  # long fold sequence stays legible (one row of 9 folds is too narrow to read).
  # Optional intuitive per-model relabelling (panels A-C carry model names). With
  # model_labels = TRUE the model codes on the y-axes / legend are swapped for reader-friendly
  # names and a decoding caption is added; FALSE/NULL keeps the raw codes.
  #
  # The names are GENERATED from the code (model_pretty_label(), 00_config.R), not looked up:
  # the hand-written table this replaced was missing a sixth of the cross-validated models, so
  # raw codes such as "Bayes-M8-dist-split-geo" were printed on the one figure whose entire
  # purpose is to avoid them. model_pretty_label() still warns, naming any code it cannot
  # decode, so a genuinely new token shows up in the run log rather than only in the PDF.
  .use_labels <- isTRUE(model_labels) ||
                 (!is.null(model_labels) && !is.logical(model_labels) && length(model_labels) > 0)
  # The published comparators take their figure-suite names (BASELINE_LBL, figure_style.R),
  # the names panel B's legend and the manuscript's Figure 2A give them, so one comparator is
  # not "Case-blind: measured cohort inflow from the epicentre" on panel A's axis and
  # "Flowminder cohort inflow from epicentre" in the legend beside it.
  relab <- function(v) {
    v <- as.character(v)
    if (!.use_labels) return(v)
    out <- model_pretty_label(v)
    bl  <- v %in% intersect(fs_baseline_methods(), names(BASELINE_LBL))
    out[bl] <- unname(BASELINE_LBL[v[bl]])
    out
  }
  y_relabel <- if (!.use_labels) NULL else scale_y_discrete(labels = function(v) relab(v))
  # ONE label size for both variants. Panels A and C hold the same ~24 rows in the same ~4.2 in
  # of height whether or not the codes are relabelled; at 11 pt the raw codes overprinted.
  yt_size   <- 8.2

  # --- 2A discrimination: per-date AUC-PR skill over all rounds, h1 & h2 ---
  # PANEL-A/C MODEL SET. Panels A and C each occupy ~4.2 in of canvas height inside the
  # composite, and a grid of more than a few dozen rows there gives well under the ~0.11 in an
  # 8.2pt label needs — the y-axis renders as a solid block of overprinted text (it did, at an
  # earlier, larger grid). Show the leading N_MODELS by 1-week skill on the panel's axis (plus
  # the featured model, wherever it ranks) at a legible row pitch; the COMPLETE grid is
  # published separately as bayes_discrimination_summary_h{1,2}.
  # A and C use the SAME model set so the two panels can be read against each other.
  N_MODELS <- 20L
  # THE CANDIDATE GRID ONLY. Panels A and C are a comparison BETWEEN models, so they must hold
  # only models that competed: the cross-validated specifications. Three kinds of row are kept
  # out (the featured model excepted).
  #   * SENSITIVITY ARMS (-gtshort/-gtlong) are refits of ONE specification at a changed
  #     assumption. Ranking them against the field invites the reader to treat "short
  #     generation time" as a rival model, which is exactly what composing the grid at a single
  #     anchor is designed to prevent, and they are barred from selection for the same reason.
  #   * ENSEMBLES (-ens-) are combinations OF the grid, not members of it, so they cannot be
  #     read as one more kernel-plus-covariate choice.
  #   * TIME-VARYING FORMS (-tvrw1/-tvfrailty/-tvrw1frailty) are stage-2 refits of the stage-1
  #     specification; only the FEATURED one (when stage 2 picked a form) is drawn, as the model
  #     the figure is about. The others are compared in FigureS9.
  # All remain fully scored and published: the arms in the sensitivity table, the forms in
  # FigureS9 and tv_comparison_*, the ensembles in the evaluation table and in
  # bayes_discrimination_summary, which is the all-methods figure.
  .sens_re <- get0("INVASION_SELECTION_EXCLUDE", ifnotfound = "-(gtshort|gtlong|tv[a-z0-9]+)$")
  .arms_re <- get0("SENSITIVITY_ARM_REGEX", ifnotfound = "-(gtshort|gtlong)$")
  # THE BASELINE RULES BELONG IN THIS PANEL. Filtering to "^Bayes" dropped them by
  # construction, so the discrimination panel showed the fitted grid ranked only against
  # itself -- the reader could not see how much of the skill is the model rather than the
  # geography. Derived from the evaluation table rather than named here, so a change to the
  # baseline set reaches the figure instead of silently leaving a hard-coded name behind.
  # The PUBLISHED set (FIGURE_BASELINE_METHODS, 00_config.R, through fs_baseline_methods()),
  # intersected with what the evaluation table actually holds -- "if present". Deriving it
  # instead as "every non-Bayes method" would sweep in Gravity-B4 and Baseline-gravity, which
  # remain scored but are not the comparators the paper displays. It used to carry its own
  # literal fallback (gravity / Flowminder inflow / travel time) — a set the config had already
  # left behind — so a config-less run would have published a third comparator set.
  .fig_baselines  <- BASELINE_METHODS
  .eval_baselines <- intersect(.fig_baselines, unique(ev$method))
  .missing_bl     <- setdiff(.fig_baselines, unique(ev$method))
  message(sprintf("[fig2] published baselines present in the evaluation table: %s%s",
                  if (length(.eval_baselines)) paste(.eval_baselines, collapse = ", ") else "NONE",
                  if (length(.missing_bl))
                    sprintf(" (absent: %s)", paste(.missing_bl, collapse = ", ")) else ""))
  .is_primary <- function(m) grepl("^Bayes", m) & !grepl("-ens-", m) &
    (!grepl(.sens_re, m) | (m == BEST_BAYES & !grepl(.arms_re, m)))
  # THE COMPARATORS' LEGEND KEY. They were called "Structural baseline", which reads as
  # "geography only" -- but two of them (Distance-B1, Adjacency-B7) read the case data to locate
  # the front. What the published set has in common is that none is a fitted model: each is a
  # fixed ranking rule (FIGURE_BASELINE_METHODS, 00_config.R). One string, used by panels A and
  # C and their captions, so the two panels cannot name the same triangles differently.
  FAM_PRIM <- "Primary specification"
  FAM_BASE <- "Baseline rule (no fitted model)"
  # PANEL A PLOTS PER-DATE AUC-PR SKILL OVER ALL ROUNDS: auc_pr_skill_date, AUC-PR scored WITHIN
  # each forecast round over that round's base rate, event-weighted across every round
  # (evaluate_invasion(), 16_invasion_eval.R). Performance is reported over all rounds; the
  # selection applies the same estimator to the post-burn-in rounds only (a deliberate choice,
  # INVASION_SELECTION_BASIS), which the selection record and Figure 2's caption state. The bars
  # are ITS 90% zone-cluster bootstrap interval
  # (auc_pr_skill_date_lo/hi: each replicate's own per-date skill). The panel used to plot
  # the POOLED all-round skill, which rewards a model whose probability LEVEL happens to be high
  # on the rounds that had invasions rather than how it orders the zones within a round, so it
  # ranked a fixed baseline rule above every fitted model on the metric the selection rejects;
  # and its bars were the AP interval over the POINT base rate, not an interval of the skill.
  # The reference is the skill a uniformly RANDOM ordering scores on the same rounds
  # (auc_pr_skill_date_chance), not 1: with few invasions per round, chance AP sits well
  # above the base rate. An evaluation table without these columns predates them, and the panel
  # is not drawn on another metric in their place.
  .need_ev <- c("auc_pr_skill_date", "auc_pr_skill_date_lo", "auc_pr_skill_date_hi",
                "auc_pr_skill_date_chance")
  .miss_ev <- setdiff(.need_ev, names(ev))
  if (length(.miss_ev))
    stop("[fig2] invasion_evaluation.csv lacks ", paste(.miss_ev, collapse = ", "),
         ": it predates the per-date AUC-PR columns. Re-run run_all.R (evaluate_invasion()) to ",
         "regenerate it; panel A is not drawn on any other metric.", call. = FALSE)
  d <- ev %>% filter(horizon %in% c(1,2),
                     .is_primary(method) | method %in% .eval_baselines) %>%
    mutate(skill = auc_pr_skill_date,
           lo = auc_pr_skill_date_lo, hi = auc_pr_skill_date_hi,
           family = factor(ifelse(method %in% .eval_baselines, FAM_BASE, FAM_PRIM),
                           levels = c(FAM_PRIM, FAM_BASE)))
  if (!nrow(d)) stop("[fig2] no candidate specifications left after excluding arms and ensembles")
  # Rank and cap the PRIMARY specifications only; every baseline present is always kept, so
  # the comparator set cannot be truncated away by a grid that grew.
  .prim <- d %>% filter(family == FAM_PRIM)
  keep_models <- .prim %>% filter(horizon == 1, is.finite(skill)) %>%
    arrange(desc(skill)) %>% pull(method)
  if (!length(keep_models))
    keep_models <- .prim %>% filter(is.finite(skill)) %>%
      arrange(desc(skill)) %>% pull(method) %>% unique()
  n_grid  <- dplyr::n_distinct(.prim$method)      # primary specifications cross-validated
  n_basel <- length(.eval_baselines)
  if (n_grid > N_MODELS)
    message(sprintf("[fig2] %d primary specifications exceed the N_MODELS=%d cap; showing the leading %d.",
                    n_grid, N_MODELS, N_MODELS))
  keep_models <- unique(c(head(keep_models, N_MODELS), intersect(BEST_BAYES, .prim$method),
                          .eval_baselines))
  d <- d %>% filter(method %in% keep_models)
  ord <- d %>% filter(horizon == 1) %>% arrange(skill) %>% pull(method)
  if (!length(ord)) ord <- d %>% arrange(skill) %>% pull(method) %>% unique()
  # COMPLETE the level set. ord is built from the h=1 rows only, so any model present at h=2
  # but not at h=1 (the featured model is added to keep_models independently of its h=1 row)
  # would become a factor NA and be dropped from the panel silently, with only a ggplot
  # "removed rows" warning that nobody reads in a batch log.
  ord <- unique(c(setdiff(unique(d$method), ord), ord))
  d <- d %>% mutate(method = factor(method, levels = ord),
                    hz = factor(horizon))
  # The random-ordering reference, one per horizon. It is a property of the rounds and zones
  # scored, not of a method, so on the shared support every method carries the same value; two
  # values at one horizon mean the rows were scored on different supports (the same check
  # prioritisation_curves() applies to recall_random).
  chance_df <- d %>% filter(is.finite(auc_pr_skill_date_chance)) %>%
    group_by(horizon) %>%
    summarise(n_val = dplyr::n_distinct(signif(auc_pr_skill_date_chance, 10)),
              chance = dplyr::first(auc_pr_skill_date_chance), .groups = "drop") %>%
    mutate(hz = factor(horizon, levels = levels(d$hz)))
  if (any(chance_df$n_val != 1L))
    stop("[fig2] the random-ordering per-date skill differs between methods at h=",
         paste(chance_df$horizon[chance_df$n_val != 1L], collapse = ","), "; the methods were ",
         "scored on different supports.", call. = FALSE)
  p2a <- ggplot(d, aes(skill, method, colour = hz)) +
    geom_vline(data = chance_df, aes(xintercept = chance, colour = hz), linetype = "22",
               linewidth = 0.45, show.legend = FALSE) +
    geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .55),
                   linewidth = .5, alpha = .55) +
    geom_point(aes(shape = family), position = position_dodge(width = .55), size = 1.9) +
    scale_colour_manual(values = HZ_COL, name = "Horizon",
                        labels = c("1 week","2 weeks")) +
    scale_shape_manual(values = setNames(c(16, 17), c(FAM_PRIM, FAM_BASE)),
                       name = NULL, drop = FALSE, guide = guide_legend(ncol = 1)) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    # theme_pub() blanks the title, subtitle and caption, so the AXIS is the only visible text:
    # it names the metric, the rounds it is scored on, and the dashed reference.
    labs(title = "Discrimination of the Bayesian invasion models",
         subtitle = paste0("Per-date AUC-PR skill = average precision within each forecast ",
                           "round / that round's base rate, event-weighted (higher = better)"),
         # Three short lines at a smaller size (below): as one line it ran into panel B's axis
         # title, and as two it was clipped at the edge of the standalone panel.
         x = paste0("Per-date AUC-PR skill (x base rate),\n",
                    "all rounds\n",
                    "dashed lines = random ordering"),
         y = NULL,
         caption = sprintf(paste0("Points = per-date AUC-PR skill over all rounds (the selection ",
                                  "applies it to the rounds after the ascertainment burn-in only, ",
                                  "deliberately); bars = its 90%% ",
                                  "zone-cluster bootstrap interval; dashed lines = the skill a random ",
                                  "ordering scores on the same rounds, per horizon. %s of the %d ",
                                  "primary specifications (circles), ordered by 1-week skill, against ",
                                  "%s (triangles). Sensitivity arms, the time-varying forms other ",
                                  "than the featured model, and the ensembles are also ",
                                  "cross-validated but are excluded here; all methods appear in ",
                                  "bayes_discrimination_summary."),
                           if (dplyr::n_distinct(.prim$method[.prim$method %in% keep_models]) >= n_grid)
                             "All" else sprintf("The leading %d",
                               dplyr::n_distinct(.prim$method[.prim$method %in% keep_models])),
                           n_grid,
                           if (n_basel) sprintf("the %d baseline rule%s (no fitted model)", n_basel,
                                                if (n_basel == 1L) "" else "s")
                           else "no baseline rule (none in the evaluation table)")) +
    theme_pub(13.5) +
    # The horizon legend and the shape legend STACK. Side by side they ran ~9 in across a
    # ~5.7 in panel and into the neighbouring panel's legend (checked by .assert_legend_fits()).
    theme(panel.grid.major.y = element_blank(),
          legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(1, "pt"),
          axis.title.x = element_text(size = 11.5, lineheight = 0.95),
          axis.text.y = element_text(size = yt_size, colour = INK)) +
    y_relabel

  # --- 2B prioritisation: invasions caught vs random targeting, 1 week ahead. The featured
  #     model against the published comparators (FIGURE_BASELINE_METHODS), from the published
  #     detection curve (see prioritisation_curves() above): recall AVERAGED per fold, with the
  #     90% fold-cluster bootstrap band on the featured model, exactly the estimand of the
  #     manuscript's Figure 2A. Colours and comparator labels from fs_panelA_scales(), so a
  #     comparator has the same name and colour here as there; the featured model keeps the
  #     name panel A's axis gives it (relab()). ---
  pc <- prioritisation_curves(horizon = 1L)
  # The featured model's legend key is WRAPPED to two lines when long: a time-varying form's
  # plain name ("..., gamma frailty") runs ~60 characters, and on one line panel B's legend was
  # wider than its half of the labelled figure (.assert_legend_fits() refused it). Wrapping at
  # 40 characters keeps it within the width of the comparators' own labels.
  sc_b <- fs_panelA_scales(BEST_BAYES, paste(strwrap(relab(BEST_BAYES), width = 40), collapse = "\n"),
                           BASELINE_METHODS, who = "pubfig2")
  curves_b <- pc$curves %>%
    mutate(method = factor(method, levels = intersect(names(sc_b$values), unique(method))))
  feat_b <- curves_b %>% filter(method == BEST_BAYES)
  p2b <- ggplot(curves_b, aes(k, recall, colour = method)) +
    geom_line(data = pc$rnd, aes(k, recall_random), linetype = "22", colour = FAINT,
              linewidth = 0.7, inherit.aes = FALSE) +
    geom_ribbon(data = feat_b, aes(ymin = recall_lo, ymax = recall_hi, fill = method),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) + geom_point(data = feat_b, size = 1) +
    annotate("text", x = max(pc$rnd$k) * 0.6, y = max(pc$rnd$recall_random) * 0.7 + 0.045,
             label = "random targeting", colour = MUTED, size = 3.8, angle = 8) +
    scale_colour_manual(values = sc_b$values, labels = sc_b$labels, name = NULL,
                        breaks = levels(curves_b$method)) +
    scale_fill_manual(values = sc_b$values, guide = "none") +
    scale_y_continuous(labels = percent_format(1), limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    labs(# STATE THE HORIZON. theme_pub() blanks titles and subtitles, so the axis is the only
         # place this panel can say it is the 1-week task — panel E's horizon is set by
         # topk_horizon (2 weeks in Figure2_labelled) and is stated on its own axis.
         x = "Zones actively monitored each round (K), 1 week ahead",
         y = "Share of true invasions caught") +
    # ONE KEY PER ROW: five entries side by side ran across panel A's half of the figure and
    # interleaved with its legend (checked by .assert_legend_fits() below).
    guides(colour = guide_legend(ncol = 1)) +
    theme_pub(13.5)

  # (Panels 2C prioritisation and 2D predicted-vs-observed have been dropped from
  #  Figure 2 at the author's request; the figure now carries panels A-E: p2a, p2b, p2e, p2f
  #  and the top-K row.)

  # --- 2E ranking accuracy: mean rank of the invaded zones, Bayesian grid, h1 & h2 ---
  # Lower = the truly-invaded zones sat nearer the TOP of the model's watch-list (a direct
  # operational readout, complementary to AUC-PR). Sourced from invasion_evaluation.csv
  # (mean_rank_of_truth); the two horizons are joined per model so the h1->h2 shift is visible.
  # Same model set as panel A (keep_models), for the same legibility reason and so the two
  # panels are read against one another rather than against two different subsets.
  # Membership is keep_models, which now carries the baseline rules as well as the
  # primary specifications, so panels A and C describe the SAME set of methods.
  # Over ALL rounds, as panel A (performance is reported over every round; the selection's own
  # rank-of-truth, on the post-burn-in rounds, is sel_mean_rank_of_truth in the table).
  if (!"mean_rank_of_truth" %in% names(ev))
    stop("[fig2] invasion_evaluation.csv has no mean_rank_of_truth; re-run the pipeline.", call. = FALSE)
  mr <- ev %>% filter(horizon %in% c(1,2), is.finite(mean_rank_of_truth),
                      method %in% keep_models) %>%
    mutate(family = factor(ifelse(method %in% .eval_baselines, FAM_BASE, FAM_PRIM),
                           levels = c(FAM_PRIM, FAM_BASE)))
  ord_mr <- mr %>% filter(horizon == 1) %>% arrange(desc(mean_rank_of_truth)) %>% pull(method)
  if (!length(ord_mr)) ord_mr <- mr %>% arrange(desc(mean_rank_of_truth)) %>% pull(method) %>% unique()
  # Complete the level set, as for panel A: mean_rank_of_truth can be finite at h=2 and missing
  # at h=1, and such a model would otherwise become a factor NA and vanish from the panel.
  ord_mr <- unique(c(setdiff(unique(mr$method), ord_mr), ord_mr))
  mr <- mr %>% mutate(method = factor(method, levels = ord_mr), hz = factor(horizon))
  mr_seg <- mr %>% dplyr::select(method, horizon, mean_rank_of_truth) %>%
    pivot_wider(names_from = horizon, values_from = mean_rank_of_truth, names_prefix = "mr")
  p2e <- ggplot(mr, aes(mean_rank_of_truth, method)) +
    { if (all(c("mr1","mr2") %in% names(mr_seg)))
        geom_segment(data = mr_seg, aes(x = mr1, xend = mr2, y = method, yend = method),
                     colour = FAINT, linewidth = 0.5, na.rm = TRUE, inherit.aes = FALSE) } +
    geom_point(aes(colour = hz, shape = family), size = 2) +
    scale_colour_manual(values = HZ_COL, name = "Horizon", labels = c("1 week","2 weeks")) +
    scale_shape_manual(values = setNames(c(16, 17), c(FAM_PRIM, FAM_BASE)),
                       name = NULL, drop = FALSE, guide = guide_legend(ncol = 1)) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.08))) +
    labs(title = "Ranking accuracy for the invaded zones",
         subtitle = "Mean rank of the truly-invaded zones among at-risk zones (lower = nearer the top)",
         x = "Mean rank of invaded zones (lower is better)", y = NULL,
         caption = sprintf(paste0("Rank of every invaded zone (tie-averaged), meaned per fold then ",
                                  "over all folds, as in panel A. ",
                                  "The same %d methods as panel A (%d primary ",
                                  "specifications, circles; %d baseline rules with no fitted model, triangles)."),
                           dplyr::n_distinct(mr$method),
                           dplyr::n_distinct(mr$method[!mr$method %in% .eval_baselines]),
                           dplyr::n_distinct(mr$method[mr$method %in% .eval_baselines]))) +
    theme_pub(13.5) +
    # The horizon legend and the shape legend STACK. Side by side they ran ~9 in across a
    # ~5.7 in panel and into the neighbouring panel's legend (checked by .assert_legend_fits()).
    theme(panel.grid.major.y = element_blank(),
          legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(1, "pt"),
          axis.text.y = element_text(size = yt_size, colour = INK)) +
    y_relabel

  # --- 2F ranking evolution: how the featured model's top watch-list churned across forecast
  # rounds. Rank-trajectory ("bump") view — each line is a health zone (coloured by province);
  # y = its 1-week invasion-risk rank that round (1 = highest). Shows which zones rose into / fell
  # out of the top of the watch-list over time. (A true alluvial needs ggalluvial, not a pipeline
  # dependency; the bump chart conveys the same ranking evolution and always renders.)
  prov_lk <- rs %>% distinct(health_zone, province)
  # DECLUTTERED 2026-09-17 (streamlining brief: "tidy up panels D and E"). 25 tracked zones
  # produced 25 overlapping trajectories plus 25 ggrepel end-labels in a panel ~7 in wide,
  # and the x-axis drew a dated break at EVERY fold cutoff. Track the current top N_TOP and
  # thin the axis labels below. (Panel E and Figure 3C show a top 12, at their own top_k.)
  N_TOP <- 15L; RANK_FLOOR <- 30L
  rk <- lfo %>% filter(method == BEST_BAYES, horizon == 1, is.finite(p_invasion),
                       !(was_active_before %in% TRUE)) %>%
    mutate(cutoff = lfo_origin(cutoff)) %>%   # label/plot rounds by FORECAST ORIGIN
    group_by(cutoff) %>% mutate(rank = rank(-p_invasion, ties.method = "min")) %>% ungroup()
  p2f <- NULL
  if (nrow(rk) > 0 && dplyr::n_distinct(rk$cutoff) >= 2) {
    last_cut <- max(rk$cutoff)
    keep_z <- rk %>% filter(cutoff == last_cut, rank <= N_TOP) %>% pull(health_zone) %>% unique()
    traj <- rk %>% filter(health_zone %in% keep_z) %>%
      left_join(prov_lk, by = "health_zone") %>%
      mutate(rank_disp = pmin(rank, RANK_FLOOR))
    # The suite's identity provinces keep their hue (prov_region_scale(), figure_style.R);
    # collapsing everything outside PROV_INT to "Other" drew Tshopo / Bas-Uele zones grey here
    # and in their own colours in the manuscript's Figure 3.
    .rg <- prov_region_scale(traj$province)
    reg_levels <- .rg$levels
    traj$region <- factor(.rg$region, levels = reg_levels)
    end_lab <- traj %>% filter(cutoff == last_cut)
    x_pad <- as.numeric(diff(range(traj$cutoff))) * 0.06 + 1
    p2f <- ggplot(traj, aes(cutoff, rank_disp, group = health_zone, colour = region)) +
      geom_line(linewidth = 1.15, alpha = 0.5, lineend = "round") +
      geom_point(size = 1.7) +
      # show.legend = FALSE: otherwise the label layer adds an "a" glyph to every key of the
      # province legend.
      ggrepel::geom_text_repel(data = end_lab, aes(label = health_zone),
        size = 3.1, direction = "y", hjust = 0, nudge_x = x_pad, box.padding = 0.12,
        segment.size = 0.2, segment.colour = FAINT, min.segment.length = 0,
        max.overlaps = 40, seed = 1, show.legend = FALSE) +
      scale_y_reverse(breaks = c(1,5,10,15,20,25,30),
                      labels = c("1","5","10","15","20","25","30+"),
                      expand = expansion(mult = c(0.04, 0.04))) +
      # DATED TICKS, THINNED TO WHAT THE PANEL CAN ACTUALLY SHOW. Labelling every other round
      # put 7 dates on a half-width panel and they overprinted ("19 May02 Jun16 Jun30 Jun...").
      # The panel is ~6 in wide and a "19 May" label is ~0.55 in at this base size, so five
      # labels is the honest ceiling; the points themselves still mark EVERY round, and the
      # first and last are always labelled because they bound the trajectory the reader is
      # asked to follow. .thin_date_breaks() keeps this correct as the fold count grows.
      scale_x_date(date_labels = "%d %b",
                   breaks = .thin_date_breaks(traj$cutoff, max_n = 5L),
                   expand = expansion(mult = c(0.03, 0.30))) +
      scale_colour_manual(values = .rg$colours, name = "Province", breaks = reg_levels) +
      labs(title = "Invasion-ranking evolution across forecast rounds",
           subtitle = sprintf("Rank of each zone's 1-week invasion risk (1 = highest); the current top %d, tracked over rounds", N_TOP),
           x = "Forecast origin (as-of date)", y = "1-week invasion-risk rank (1 = highest)",
           caption = "Each line = a zone currently in the top watch-list; ranks worse than 30 shown at the 30+ baseline.") +
      theme_pub(13.5) +
      theme(legend.position = "right",
            panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
            # A slight angle buys the width the thinning cannot: five horizontal dates still
            # sit close on this panel once the ggrepel label gutter is subtracted.
            axis.text.x = element_text(angle = 30, hjust = 1))
  }
  if (is.null(p2f))
    p2f <- ggplot() + annotate("text", x = 0, y = 0, label = "Ranking evolution unavailable\n(need >= 2 folds)",
                               colour = MUTED, size = 3) + theme_void()

  save_dual(p2a, paste0("F2A_discrimination", file_suffix),   5.6, 4.6)
  save_dual(p2b, paste0("F2B_prioritisation", file_suffix),   5.2, 4.8)
  save_dual(p2e, paste0("F2E_mean_rank", file_suffix),        5.6, 4.6)
  save_dual(p2f, paste0("F2F_ranking_evolution", file_suffix), 7.2, 4.8)

  # Bottom row: the top-K predicted invasion risks per fold, coloured by realised outcome
  # (the F_topk_folds bar chart, embedded without its standalone title). Spans full width;
  # topk_nrow wraps the fold facets over multiple rows for legibility.
  # max_facets = 6: embedded in Figure 2 this row previously carried one narrow facet per
  # fold (14 of them), at which width the per-zone labels are unreadable in print. Six evenly
  # spaced rounds, first and last always kept, keep the trajectory legible. The
  # F_topk_folds_h1/h2 panels built below show every round (written only where the FIGURE_KEEP
  # gate retains them).
  p_topk1 <- build_topk_folds(horizon = topk_horizon, facet_nrow = topk_nrow,
                              save = FALSE, embed = TRUE, max_facets = 6L)

  # Decoding caption for the reader-friendly labelled variant (added only when model_labels set).
  # The caption must decode the tokens the panels ACTUALLY carry. It used to list every token
  # the grid can generate, including "+ suspected" and "ensemble", which panels A and C no
  # longer show; a glossary for absent terms sends the reader looking for rows that are not
  # there. Each clause is now gated on its token appearing in the drawn labels.
  # relab() is the same function the y axis uses, so the glossary is gated on exactly the
  # strings the reader sees rather than on the raw method codes.
  .lab_txt <- paste(relab(unique(as.character(d$method))), collapse = " ")
  .clause <- function(tok, txt) if (grepl(tok, .lab_txt, fixed = TRUE)) txt else NULL
  cap2 <- if (!.use_labels) NULL else stringr::str_wrap(paste0(
    "Model = the spatial mobility kernel driving invasion risk: gravity, radiation and their ",
    "composites with Flowminder mobile-phone data, and a three-kernel consensus of gravity, ",
    "cohort-calibrated gravity and radiation. ",
    paste(Filter(Negate(is.null), list(
      .clause("Cohort", paste0(
        "\"Cohort\" = Flowminder presence of the epicentre's subscriber cohorts")),
      .clause("Short-trip", paste0(
        "\"Short-trip\" = Flowminder short trips out of the epicentre zones")),
      .clause("relocation OD", "\"relocation OD\" = Flowminder origin-destination relocations"),
      .clause("Cohort-calibrated gravity", paste0(
        "\"Cohort-calibrated gravity\" = gravity whose decay is fitted to the cohort presence")),
      .clause("gamma frailty", paste0(
        "\"gamma frailty\" = the featured model's zone-level gamma frailty (zones differ in ",
        "susceptibility; a zone exposed without being invaded becomes less likely to be)")),
      .clause("beta random walk", "\"beta random walk\" = a weekly random walk in the import coefficient"),
      .clause("covariates", paste0(
        "\"+ covariates\" = geographic & social covariates modulate the import-to-invasion rate ",
        "(base models use a single constant rate)")),
      .clause("road-km", "\"road-km\" = road-distance rather than travel-time deterrence"),
      # Quotes the axis wording ("source-filled", model_pretty_label()) exactly: the glossary
      # used to define "source-cell fill", a term printed nowhere on the figure.
      .clause("source-filled", paste0(
        "\"source-filled\" = destinations the mobility source could not observe are taken ",
        "from the base kernel rather than assumed zero")),
      .clause("suspected", "\"+ suspected\" = suspected-case leading indicators"),
      .clause("ensemble", "\"ensemble\" = mean / median across models"))),
      collapse = "; "),
    ". Sensitivity arms, time-varying refits other than the featured model, and ensembles are ",
    "not shown here; see bayes_discrimination_summary."),
    width = 180)

  # No embedded figure title/subtitle — only the panel tags. Figure 2 carries the
  # discrimination (A), prioritisation (B), ranking-accuracy (C), ranking-evolution (D) and
  # per-fold top-K outcome (E) panels. Tags relabel A-E in layout order.
  # Panel E grows when its fold facets are wrapped over multiple rows.
  topk_h <- if (isTRUE(topk_nrow >= 2L)) 1.05 * topk_nrow else 1.1
  fig2 <- (p2a | p2b) / (p2e | p2f) / p_topk1 +
    plot_layout(heights = c(1, 1.05, topk_h)) +
    plot_annotation(tag_levels = "A", caption = cap2,
      theme = theme(plot.caption = element_text(size = 8.2, colour = MUTED, hjust = 0,
                                                margin = margin(t = 6))))
  # Legends checked in inches before anything is written: panel A's two legends side by side
  # used to run straight across panel B's, and ggplot neither clips nor wraps a guide box
  # wider than its plot (.assert_legend_fits(), figure_style.R — the check the manuscript
  # suite's Figure 2 uses). Rows of the layout above, each at the full figure width.
  .assert_legend_fits(list(p2a, p2b), c(1, 1), fig_w, labels = c("A", "B"),
                      who = paste0("pubfig2", file_suffix))
  .assert_legend_fits(list(p2e, p2f), c(1, 1), fig_w, labels = c("C", "D"),
                      who = paste0("pubfig2", file_suffix))
  save_dual(fig2, paste0("Figure2", file_suffix), fig_w, fig_h, dir = FIG_DIR)
  invisible(fig2)
}

# =============================================================================
# FIGURE 3  (spatial forecast; from risk-score CSV + shapefile)
# =============================================================================
# bivariate 4x4 palette (Stevens) via bilinear corner interpolation
bivar_pal <- function() {
  cc <- function(h) grDevices::col2rgb(h)[,1]
  c00<-cc("#e8e8e8"); c10<-cc("#5ac8c8"); c01<-cc("#be64ac"); c11<-cc("#3b4994")
  g <- expand.grid(bx = 1:4, by = 1:4)
  cols <- apply(g, 1, function(r){ fx<-(r[1]-1)/3; fy<-(r[2]-1)/3
    v <- (1-fx)*(1-fy)*c00 + fx*(1-fy)*c10 + (1-fx)*fy*c01 + fx*fy*c11
    grDevices::rgb(v[1],v[2],v[3], maxColorValue = 255) })
  setNames(cols, paste0(g$bx, "-", g$by))
}

build_fig3 <- function() {
  BIP <- bivar_pal()
  hz_lab <- function(h) if (h==1) "1 week ahead" else "2 weeks ahead"

  # per-horizon derived fields
  prep <- function(h) {
    d <- rs %>% filter(horizon == h) %>%
      mutate(active = was_active_before %in% TRUE,
             p = ifelse(active, NA, p_case_invasion),
             width = ifelse(active, NA, p_case_hi - p_case_lo))
    n_ar <- sum(!d$active & is.finite(d$p))
    d %>% mutate(rank_p = ifelse(is.finite(p), rank(-p, ties.method="min", na.last="keep"), NA),
                 pct = 100*(1 - (rank_p - 1)/max(n_ar - 1, 1)),
                 pct01 = pct / 100,                     # same rank measure on a 0-1 scale
                 rr_log = ifelse(active, NA, pmax(rr01_nat, RR_FLOOR)),  # true RR, floored for log10
                 # BOTH quartiles over the at-risk zones: rr01_nat is NA on an affected zone,
                 # so vulnerability ranked over all 519 zones against risk over the at-risk
                 # set put the two axes of one bivariate class on different populations.
                 bx = qtile(ifelse(active, NA_real_, rr01_nat)),
                 by = qtile(ifelse(active, NA_real_, V)),
                 bikey = ifelse(is.na(bx)|is.na(by), NA, paste0(bx,"-",by)))
  }
  # THE ALREADY-AFFECTED ZONES GET A LEGEND KEY. They are drawn in a flat grey over the fill
  # scale, and with no key the reader could take the grey for a low value of it. The key is a
  # constant `linetype` mapping, which a continuous fill scale cannot carry: the legend draws the
  # layer's own fill (AFFECTED_FILL) as its swatch, and every map of a pair gives the same key,
  # so guides = "collect" shows it once (panels A and B). Panel E has no ggplot legend to join,
  # so it keys the same grey as a swatch in its bivariate legend instead (d_map, below).
  AFF_KEY <- "Already affected\n(not at risk)"
  aff_overlay <- function(mp) list(
    geom_sf(data = mp %>% filter(was_active_before %in% TRUE), aes(linetype = AFF_KEY),
            fill = AFFECTED_FILL, colour = "white", linewidth = 0.08),
    scale_linetype_manual(values = setNames("solid", AFF_KEY), name = NULL,
                          guide = guide_legend(order = 9,
                                               override.aes = list(colour = FAINT))))
  # THE HORIZON IS PRINTED OVER EACH MAP. theme_map() blanks plot.title, so ggtitle(hz_lab(h))
  # was computed and discarded: the two maps of each pair carried nothing saying which was
  # 1 week and which 2 weeks ahead. Re-enabled for these sub-panel labels only (the figure
  # itself still carries no title, only the A-E tags).
  map_title <- theme(plot.title = element_text(size = 11.5, hjust = 0.5, colour = INK,
                                               margin = margin(b = 3)))

  one_map <- function(h, fillvar, scale_fn, title) {
    mp <- join_map(prep(h))
    ggplot(mp) +
      geom_sf(aes(fill = .data[[fillvar]]), colour = "white", linewidth = 0.08) +
      aff_overlay(mp) + scale_fn + ggtitle(title) + theme_map(13.5) + map_title
  }

  # --- 3A relative-risk map (viridis; TRUE relative invasion risk on a log10 0-1 scale) ---
  # rr01_nat is heavily right-skewed (relative risk is concentrated near the epicentre), so a
  # linear 0-1 fill collapses the national map to a single dark hue. A log10 fill shows HONEST
  # magnitudes while staying legible. log10 is undefined at 0 and unbounded below, so the scale
  # is floored at RR_FLOOR: every zone at or below it renders in the floor colour and the first
  # legend tick reads "<=0.001" rather than implying an exact value.
  sc_rank <- scale_fill_viridis_c(option="viridis", trans="log10", limits=c(RR_FLOOR, 1),
              breaks=c(0.001,0.01,0.1,1), labels=c("<=0.001","0.01","0.1","1"),
              na.value=NA_FILL, name="Relative\ninvasion risk\n(log, 0-1)", direction=1)
  a1 <- one_map(1, "rr_log", sc_rank, hz_lab(1)); a2 <- one_map(2, "rr_log", sc_rank, hz_lab(2))
  p3a <- (a1 | a2) + plot_layout(guides="collect") &
    theme(legend.position="right")

  # --- 3B uncertainty map (CrI width; rocket, dark = more uncertain) ---
  # The fill is capped at the 99th percentile of the at-risk widths (both horizons), floored at
  # 0.05, so a few very wide intervals do not wash out the rest; wider ones are SQUISHED to the
  # cap colour. The cap is rounded UP to two significant digits, so its label is exact, and the
  # top key reads ">=cap" whenever any width exceeds it — panel A labels its floor "<=0.001" for
  # the same reason. Without that the darkest key read as the maximum width.
  w_all <- rs %>% filter(!(was_active_before %in% TRUE)) %>%
    mutate(w = p_case_hi - p_case_lo) %>% pull(w)
  w_q   <- max(0.05, stats::quantile(w_all, .99, na.rm = TRUE, names = FALSE))
  w_e   <- 10^(floor(log10(w_q)) - 1)
  w_cap <- ceiling(round(w_q / w_e, 8)) * w_e
  w_brk <- pretty(c(0, w_cap), n = 4)
  w_brk <- c(w_brk[w_brk < 0.8 * w_cap], w_cap)         # room for the cap's own label
  w_lab <- format(w_brk, drop0trailing = TRUE, trim = TRUE)
  if (any(w_all > w_cap, na.rm = TRUE)) w_lab[length(w_lab)] <- paste0(">=", w_lab[length(w_lab)])
  sc_w <- scale_fill_viridis_c(option="rocket", direction=-1, limits=c(0, w_cap),
            breaks = w_brk, labels = w_lab,
            oob=scales::squish, na.value=NA_FILL, name="90% CrI\nwidth")
  b1 <- one_map(1, "width", sc_w, hz_lab(1)); b2 <- one_map(2, "width", sc_w, hz_lab(2))
  p3b <- (b1 | b2) + plot_layout(guides="collect") & theme(legend.position="right")

  # --- 3C top-zone posterior probability + 90% CrI (forest) — 1- AND 2-week MERGED into one
  # panel (dodged, coloured by horizon), so the two horizons are compared directly rather than in
  # two separate facets. Zones = the top by 1-week risk (the current watch-list). ---
  ord_z <- rs %>% filter(horizon == 1, !(was_active_before %in% TRUE), is.finite(p_case_invasion)) %>%
    slice_max(p_case_invasion, n = 12) %>% arrange(p_case_invasion) %>% pull(health_zone)
  topz <- rs %>% filter(health_zone %in% ord_z, horizon %in% c(1, 2), is.finite(p_case_invasion)) %>%
    mutate(zone = factor(health_zone, levels = ord_z), hz = factor(horizon))
  pd_c <- position_dodge(width = 0.55)
  # THE AXIS NAMES ITS PROBABILITY SCALE, read off the risk table's own stamp (as the manuscript
  # Figure 3 does); a table with no single stamp is refused rather than labelled.
  .ps3 <- if ("prob_scale" %in% names(topz)) unique(as.character(topz$prob_scale)) else character(0)
  .scale3 <- if (length(.ps3) == 1L && grepl("^recalibrated", .ps3)) "recalibrated"
             else if (length(.ps3) == 1L && grepl("^raw", .ps3)) "raw"
             else stop("[fig3] the live risk table carries no single prob_scale stamp (",
                       paste(.ps3, collapse = ", "), "), so panel C cannot say which ",
                       "probability scale it draws. Re-run the pipeline.", call. = FALSE)
  p3c <- ggplot(topz, aes(p_case_invasion, zone, colour = hz)) +
    geom_linerange(aes(xmin = p_case_lo, xmax = p_case_hi), linewidth = 0.6,
                   position = pd_c, orientation = "y") +
    geom_point(size = 1.9, position = pd_c) +
    scale_colour_manual(values = HZ_COL, name = "Horizon", labels = c("1 week","2 weeks")) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA),
                       expand = expansion(mult = c(0.01, 0.08))) +
    labs(title = "Highest-risk zones: 1- vs 2-week invasion probability",
         subtitle = "Top at-risk zones by 1-week risk; posterior mean and 90% CrI at both horizons",
         x = sprintf("Invasion probability (%s)", .scale3), y = NULL) +
    theme_pub(13.5) + theme(panel.grid.major.y = element_blank(),
                        axis.text.y = element_text(colour = INK, size = 11.5))

  # --- 3D preparedness-priority scatter (1 week ahead): invasion likelihood x vulnerability,
  # sized by the composite priority, coloured by province (reconstruct of bayes_priority_scatter_h1). ---
  # Identity provinces keep their hue, as in every other figure (prov_region_scale()).
  psc <- rs %>% filter(horizon == 1, !(was_active_before %in% TRUE),
                       is.finite(V), is.finite(rr01_nat))
  .rg3 <- prov_region_scale(psc$province)
  psc$region <- factor(.rg3$region, levels = .rg3$levels)
  lab_ps <- psc %>% arrange(desc(priority)) %>% head(12)
  p3d <- ggplot(psc, aes(V, rr01_nat)) +
    geom_point(aes(size = priority, colour = region), alpha = 0.78) +
    ggrepel::geom_text_repel(data = lab_ps, aes(label = health_zone), size = 3.1, colour = INK,
      box.padding = 0.3, max.overlaps = 30, seed = 1, min.segment.length = 0,
      segment.colour = FAINT, segment.size = 0.2) +
    scale_colour_manual(values = .rg3$colours, name = "Province", breaks = .rg3$levels) +
    scale_size_area(max_size = 7, name = "Priority") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = scales::label_number(accuracy = 0.01)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = scales::label_number(accuracy = 0.01),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = "Preparedness priority (1 week ahead): invasion risk x vulnerability",
         subtitle = "Upper-right = likely to be invaded AND vulnerable / under-resourced = highest priority",
         x = "Vulnerability",
         y = "Relative invasion risk") +
    theme_pub(13.5) + theme(legend.position = "right")

  # --- 3E prob x vulnerability bivariate choropleth + 4x4 legend (was panel D) ---
  d_map <- function(h){
    mp <- join_map(prep(h))
    ggplot(mp) + geom_sf(aes(fill = bikey), colour="white", linewidth=0.08) +
      aff_overlay(mp) +
      scale_fill_manual(values = BIP, na.value = NA_FILL, guide = "none") +
      # this pair's "already affected" key is a swatch in the bivariate legend (below)
      guides(linetype = "none") +
      ggtitle(hz_lab(h)) + theme_map(13.5) + map_title
  }
  leg_df <- expand.grid(bx=1:4, by=1:4) %>% mutate(bikey = paste0(bx,"-",by))
  # The maps' "already affected" grey is keyed by a swatch ABOVE the 4x4 grid (clear of the axis
  # titles): a collected ggplot legend for it landed on top of this grid in the composite.
  legend_biv <- ggplot(leg_df, aes(bx, by, fill = bikey)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    annotate("tile", x = 1, y = 5.05, width = 1, height = 0.6, fill = AFFECTED_FILL,
             colour = "white") +
    annotate("text", x = 1.65, y = 5.05, hjust = 0, size = 2.9, colour = MUTED,
             lineheight = 0.9, label = "Already affected\n(not at risk)") +
    scale_fill_manual(values = BIP, guide = "none") + coord_fixed(clip = "off") +
    labs(x = "Invasion prob ->", y = "Vulnerability ->") +
    theme_minimal(base_size = 9, base_family = base_family) +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          axis.title = element_text(size = 8.5, colour = MUTED),
          plot.margin = margin(2,2,2,2))
  p3e <- (d_map(1) | d_map(2) | wrap_elements(full = legend_biv)) +
    plot_layout(widths = c(1, 1, 0.42))

  # individual saves
  save_dual(p3a, "F3A_rank_map",             8.0, 4.2)
  save_dual(p3b, "F3B_uncertainty_map",      8.0, 4.2)
  save_dual(p3c, "F3C_top_zone_uncertainty", 5.8, 4.8)
  save_dual(p3d, "F3D_priority_scatter",     6.6, 4.8)
  save_dual(p3e, "F3E_prob_vuln_bivariate",  8.4, 4.2)

  # No embedded figure title/subtitle — only the A-E panel tags (each map of a pair names its
  # own horizon, above). The descriptive text, including the featured model, is the figure
  # caption in the manuscript, not part of the figure.
  fig3 <- (wrap_elements(full = p3a) / wrap_elements(full = p3b) /
           (wrap_elements(full = p3c) | wrap_elements(full = p3d)) /
           wrap_elements(full = p3e)) +
    plot_layout(heights = c(1, 1, 1.1, 1)) +
    plot_annotation(tag_levels = "A")
  save_dual(fig3, "Figure3", 10.0, 14.7, dir = FIG_DIR)
  invisible(fig3)
}

# (Within-facet reordering, reorder_within() / tidytext_scale_y(), is in figure_style.R.)

# =============================================================================
# FIGURE: top predicted invasion probabilities per fold, coloured by outcome
# Companion to the diagnostic bayes_lfo_forecast_vs_outcome tile: instead of a
# zone x cutoff heatmap, this separates each leave-future-out fold into its OWN
# facet and draws a bar chart of the highest-risk zones that round, coloured by
# whether that zone actually recorded its first case (invaded) or not. It reads
# as a per-round precision/reliability view: tall bars that are highlighted =
# confident calls that came true; tall grey bars = the round's false alarms.
# =============================================================================
#' @param save    write the standalone F_topk_folds_h<h> panel (PDF+PNG, through the
#'   FIGURE_KEEP gate).
#' @param embed   drop the embedded title/subtitle and enlarge the base size so the
#'   panel sits cleanly as a tagged row inside Figure 2 (which carries no titles).
#' @param max_facets show at most this many forecast rounds, evenly spaced and ALWAYS
#'   including the first and last. NA (default) keeps every round. Embedded in Figure 2 the
#'   full set was 14 narrow facets whose per-zone labels were unreadable at print size; the
#'   standalone F_topk_folds panels keep every round, since they have the width for it.
build_topk_folds <- function(top_k = 12L, horizon = 1L, save = TRUE, embed = FALSE,
                             facet_nrow = 1L, max_facets = NA_integer_) {
  OUT_HIT  <- "#D55E00"   # invaded (first case occurred)  — Okabe-Ito vermilion
  OUT_MISS <- "grey75"    # not invaded this round
  OUT_COL  <- c("Invaded (first case)" = OUT_HIT, "Not invaded" = OUT_MISS)

  d <- lfo %>%
    filter(method == BEST_BAYES, horizon == !!horizon, is.finite(p_invasion),
           !(was_active_before %in% TRUE))
  if (!nrow(d)) { message("[topk_folds] no at-risk rows; skipped."); return(invisible(NULL)) }

  # top-K zones by predicted probability WITHIN each fold. Facet header = round + date, plus
  # "(uncorrected)" on a round drawn without a recalibration factor and "(provisional)" on a
  # round whose outcome has not settled (both below).
  d <- d %>% mutate(cutoff = as.Date(cutoff))
  # WHICH ROUNDS ARE DRAWN UNCORRECTED. The bars are on the prequentially recalibrated scale,
  # but the earliest rounds have no settled training fold to fit delta on (delta_estimable
  # FALSE) and their p_recal is the RAW probability. Facets share no x scale, so nothing on the
  # page said that one round's bars were on a different scale from its neighbours'. Read per
  # round off the featured model's own rows at this horizon -- the rows drawn -- exactly as
  # make_manuscript_figures.R's .est_by_cut/.uncorr() do for Figure 2B and FigureS14. On a raw
  # pass nothing is recalibrated, so nothing is "uncorrected" relative to the rest; an LFO frame
  # without the column predates recalibration and every round is on the one (raw) scale.
  .uncorr_cuts <- if (identical(PSCALE_COL, "p_recal") && "delta_estimable" %in% names(d)) {
    .est <- tapply(d$delta_estimable %in% TRUE, as.character(d$cutoff), any)
    as.Date(names(.est)[!.est])
  } else as.Date(character(0))
  # WHICH ROUNDS ARE PROVISIONAL. A round whose outcome window closed fewer than
  # LFO_RELIABLE_EVAL_AGE_DAYS before the analysis date (eval_reliable FALSE, written by
  # run_invasion_lfo()) is scored, but its invasions are not all reported yet: a grey "Not
  # invaded" bar there can still turn out to have been invaded. Drawn unmarked, such a round read
  # as final beside the settled ones; FigureS7 draws the same rounds hollow. Read per round off
  # the rows drawn. A frame without the column cannot say which rounds have settled, so it is
  # refused rather than drawn as if every outcome were final.
  if (!"eval_reliable" %in% names(d))
    stop("[topk_folds] lfo_results.rds has no eval_reliable column, so unsettled rounds cannot ",
         "be marked provisional; re-run run_all.R.", call. = FALSE)
  .rel <- tapply(d$eval_reliable %in% TRUE, as.character(d$cutoff), all)
  .prov_cuts <- as.Date(names(.rel)[!.rel])
  # On its own line when the embedded row is crowded: six facets across Figure 2 leave ~1.9 in
  # per strip, and "R1 - 15 May (uncorrected)" on one line ran into the next round's strip
  # (strip.clip is off). Three per row (Figure2_labelled) carry it on one line. Evaluated when
  # the labels are built, after fold_ord is known.
  .round_tag <- function(cut) {
    sep <- if (embed && ceiling(length(fold_ord) / facet_nrow) > 3L) "\n" else " "
    cut <- as.Date(cut)
    vapply(seq_along(cut), function(i) {
      tg <- c(if (cut[i] %in% .uncorr_cuts) "uncorrected",
              if (cut[i] %in% .prov_cuts)   "provisional")
      if (length(tg)) paste0(sep, "(", paste(tg, collapse = ", "), ")") else ""
    }, character(1))
  }
  fold_ord_all <- sort(unique(d$cutoff))
  # Round NUMBERS are assigned over ALL folds before any subsetting, so a displayed facet
  # keeps its true round index (dropping rounds must not renumber the survivors).
  round_no <- setNames(seq_along(fold_ord_all), as.character(fold_ord_all))
  fold_ord <- if (!is.na(max_facets) && length(fold_ord_all) > max_facets)
    fold_ord_all[unique(round(seq(1, length(fold_ord_all), length.out = max_facets)))]
    else fold_ord_all
  if (length(fold_ord) < length(fold_ord_all))
    message(sprintf("[topk_folds] showing %d of %d rounds (evenly spaced, first and last kept).",
                    length(fold_ord), length(fold_ord_all)))
  d <- d %>% filter(cutoff %in% fold_ord)
  # Round number AND date, embedded or not. The justification for dropping the date ("the dates
  # are on panel D's axis and in the caption") was false: panel D shows dates without round
  # numbers, and the only surviving caption is the model-decoding text.
  lab_lk   <- setNames(if (embed) sprintf("R%d - %s%s", round_no[as.character(fold_ord)],
                                          format(lfo_origin(fold_ord), "%d %b"),
                                          .round_tag(fold_ord))
                       else sprintf("Round %d - %s%s", round_no[as.character(fold_ord)],
                                    format(lfo_origin(fold_ord), "%d %b"), .round_tag(fold_ord)),
                       as.character(fold_ord))
  if (any(fold_ord %in% .uncorr_cuts))
    message(sprintf("[topk_folds] h=%d: %d of %d drawn round(s) carry no fitted recalibration factor and are labelled \"(uncorrected)\".",
                    horizon, sum(fold_ord %in% .uncorr_cuts), length(fold_ord)))
  if (any(fold_ord %in% .prov_cuts))
    message(sprintf("[topk_folds] h=%d: %d of %d drawn round(s) have not settled (eval_reliable FALSE) and are labelled \"(provisional)\".",
                    horizon, sum(fold_ord %in% .prov_cuts), length(fold_ord)))

  topk <- d %>% group_by(cutoff) %>%
    slice_max(p_invasion, n = top_k, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(outcome = factor(ifelse(is_new_invasion == 1, "Invaded (first case)", "Not invaded"),
                            levels = names(OUT_COL)),
           fold_lab = factor(lab_lk[as.character(cutoff)], levels = lab_lk[as.character(fold_ord)]),
           zone_w = reorder_within(health_zone, p_invasion, fold_lab))

  base_sz <- if (embed) 13 else 8.6
  ttl <- if (embed) NULL else sprintf("Top-%d predicted invasion risks per forecast round, by realised outcome", top_k)
  sub <- if (embed) NULL else sprintf("Featured model %s, %d week ahead; bar = predicted P(first case), highlighted where the zone was invaded that round",
                                      BEST_BAYES, horizon)
  p <- ggplot(topk, aes(p_invasion, zone_w, fill = outcome)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    facet_wrap(~ fold_lab, scales = "free_y", nrow = facet_nrow) +
    tidytext_scale_y() +
    scale_fill_manual(values = OUT_COL, name = NULL, drop = FALSE) +
    # Embedded, each facet's bar area is well under an inch wide: three or four percent labels
    # ran together ("0%20%40%60%"). Two in-range breaks are legible there.
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.06)),
                       breaks = if (embed) function(lim) { b <- pretty(lim, n = 2)
                                                          b[b >= lim[1] & b <= lim[2]] }
                                else scales::pretty_breaks(3)) +
    # The horizon must survive `embed = TRUE`. Embedded in Figure 2 this panel drops its title
    # AND subtitle, which were the only text naming its horizon — and panel B is the 1-week task
    # while panels A/C carry a 1-week/2-week legend, so the axis has to say which this is.
    labs(title = ttl, subtitle = sub,
         x = sprintf("Predicted invasion probability, P(first case), %d week%s ahead%s",
                     horizon, if (horizon == 1L) "" else "s",
                     # The bars' scale: recalibrated, except rounds whose facet says "uncorrected".
                     if (!identical(PSCALE_COL, "p_recal")) " (raw)"
                     else if (length(.uncorr_cuts) && any(fold_ord %in% .uncorr_cuts))
                       " (recalibrated; rounds marked uncorrected are raw)"
                     else " (recalibrated)"), y = NULL) +
    theme_pub(base_sz) +
    # strip.clip = "off" keeps the per-round facet header from being clipped to its narrow strip
    theme(panel.grid.major.y = element_blank(),
          panel.spacing.x = unit(9, "pt"),
          strip.clip = "off",
          axis.text.y = element_text(size = if (embed) 9.2 else 6.6, colour = INK),
          axis.text.x = element_text(size = if (embed) 9.2 else 6.6),
          legend.position = "top")
  if (!embed)
    p <- p + theme(plot.title    = element_text(size = 9.6, face = "bold", colour = INK,
                                                margin = margin(b = 2)),
                   plot.subtitle = element_text(size = 8, colour = MUTED, margin = margin(b = 6)))

  if (save) save_dual(p, sprintf("F_topk_folds_h%d", horizon),
                      2.5 * ceiling(length(fold_ord) / facet_nrow) + 1.2, 4.6 * facet_nrow)
  invisible(p)
}

# -----------------------------------------------------------------------------
# RUN
# -----------------------------------------------------------------------------
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("need ggrepel")
run <- function(f, nm) tryCatch({ message("== ", nm, " =="); f(); TRUE },
                                error = function(e){ message("!! ", nm, " FAILED: ",
                                conditionMessage(e)); FALSE })
ok1 <- run(build_fig1, "Figure 1")
oka <- run(build_f1a_allflows, "Figure 1A (all outflows)")
ok2 <- run(build_fig2, "Figure 2")
# Separate reader-friendly variant: same panels, intuitive model names on A-C + decoding caption.
# Reader-friendly variant: intuitive model names on A/C, plus a 2-week-horizon panel E
# with the per-fold facets wrapped over two rows (legible where a single row is too narrow).
# The taller two-row panel E needs extra canvas height.
ok2b <- run(function() build_fig2(model_labels = TRUE, file_suffix = "_labelled",
                                  fig_w = 12.6, fig_h = 17.4,
                                  topk_horizon = 2L, topk_nrow = 2L), "Figure 2 (labelled)")
ok3 <- run(build_fig3, "Figure 3")
okt <- run(function() { build_topk_folds(horizon = 1L); build_topk_folds(horizon = 2L) },
           "Top-K per-fold outcome bars")
# `oka` was built and then never reported: Figure 1A (all outflows) could fail and leave no
# trace in the [done] line at all. Every panel this script builds is listed here.
message(sprintf(paste0("\n[done] Figure1:%s  Figure1A-allflows:%s  Figure2:%s  ",
                       "Figure2-labelled:%s  Figure3:%s  TopKfolds:%s  ->  %s"),
                ok1, oka, ok2, ok2b, ok3, okt, FIG_DIR))

# A FAILED PANEL MUST FAIL THE SCRIPT -- see the same note in make_manuscript_figures.R.
# `run()` catches so one broken panel does not cost the others, but reaching the end and
# exiting 0 told run_all.R the suite was fine while the previous run's PDF stayed on disk.
.panels <- c(Figure1 = ok1, `Figure1A-allflows` = oka, Figure2 = ok2,
             `Figure2-labelled` = ok2b, Figure3 = ok3, TopKfolds = okt)
if (any(!.panels))
  stop(sum(!.panels), " of ", length(.panels), " panel(s) failed: ",
       paste(names(.panels)[!.panels], collapse = ", "),
       ". Their files on disk are from an earlier run -- see the messages above for each cause.",
       call. = FALSE)

}   # end of the script-run guard

