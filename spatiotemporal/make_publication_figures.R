# =============================================================================
# make_publication_figures.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
# Main-text multi-panel figures (Nature/Science aesthetic), rebuilt from SAVED
# data (no model re-fitting) under ONE unified publication theme.
#
#   Figure 1  Outbreak to date   A) cumulative-case map  B) onset time series
#   Figure 2  Model performance  A) discrimination B) calibration
#                                C) prioritisation  D) predicted-vs-observed
#   Figure 3  Spatial forecast   A) rank map B) uncertainty map
#                                C) top-zone uncertainty D) prob x vulnerability
#
# Outputs -> outputs/key_outputs/figures/{Figure1,2,3}.{pdf,png} and
#            outputs/key_outputs/figures/panels/<panel>.{pdf,png}
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
stopifnot(exists("OUT_DIR"), exists("OUTBREAK_START"), exists("EPICENTRE_ZONES"),
          exists("load_linelist", mode = "function"),
          exists("effective_onset_sample_delay", mode = "function"))
.cfg <- function(name, default) if (exists(name)) get(name) else default

OUT       <- .cfg("OUT_DIR", file.path(HERE, "outputs"))
# Probability-scale switch (forecast_scale.R). Default "recalibrated": the primary set.
# FORECAST_SCALE=raw re-runs the identical build onto a raw/ sibling, same basenames.
source(file.path(HERE, "spatiotemporal", "forecast_scale.R"))
FIG_DIR   <- fs_out_dir(file.path(OUT, "key_outputs", "figures"))
PANEL_DIR <- fs_out_dir(file.path(OUT, "key_outputs", "figures", "panels"))
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

BEST_BAYES_FALLBACK <- "Bayes-M14-fill-med"   # only if the featured model cannot be derived;
                                              # the fill family is the default, so the fallback must name a kernel the default grid builds
OUTBREAK_START <- .cfg("OUTBREAK_START", as.Date("2026-04-30"))
ANALYSIS_DATE  <- .cfg("ANALYSIS_DATE", Sys.Date())        # as-of date; panels run up to here
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
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
AFFECTED_FILL <- "grey78"; NA_FILL <- "grey93"
RR_FLOOR <- 1e-3   # lower limit of the log10 relative-risk scale (log10 is undefined at 0)

# Okabe-Ito (CVD-safe categorical); fixed province identity map
OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
PROV_COL <- c("Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73",
              "Sud-Kivu"="#CC79A7","Other"="#7A7A7A")
# series colours for the evaluation panels
SER_COL <- c("Mean predicted P"="#0072B2","Observed invasion fraction"="#111111",
             "Mean P at invaded zones"="#D55E00")
HZ_COL  <- c("1"="#3B4CC0","2"="#B4413C")            # horizon 1 / 2

# Reader-friendly model names for the labelled Figure 2 variant. These are GENERATED from the
# method code by model_pretty_label() (00_config.R), which composes the name from the same
# tokens bayes_default_grid() composes the label from — kernel family, road-km / source-fill /
# origin-split qualifiers, covariate set, beta_t process, generation-time arm.
#
# IT USED TO BE A 100-ENTRY HAND-WRITTEN LOOKUP, and it drifted: 12 of the 57 cross-validated
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

base_family <- "sans"

theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    # Titles / subtitles / captions are intentionally BLANK: the main-text figures carry no
    # embedded titles; all descriptive text lives in FIGURE_CAPTIONS.md (external caption).
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
                                   margin = margin(3,3,3,3)),
    plot.tag        = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin     = margin(6, 8, 6, 6)
  )
}
theme_map <- function(base = 8.6) {
  theme_void(base_size = base, base_family = base_family) %+replace% theme(
    plot.title    = element_text(size = base - 0.4, colour = INK, hjust = 0.5,
                                 margin = margin(b = 2), face = "bold"),
    legend.position = "right",
    legend.title  = element_text(size = base - 1.4, colour = MUTED),
    legend.text   = element_text(size = base - 1.8, colour = INK),
    legend.key.height = unit(20, "pt"), legend.key.width = unit(7, "pt"),
    plot.tag      = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin   = margin(2, 2, 2, 2)
  )
}

# save vector PDF + 600-dpi PNG
save_dual <- function(p, name, w, h, dir = PANEL_DIR) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  # Gate on the FULL destination path, not the bare stem: FIGURE_DROP entries are
  # "<directory>/<stem>" and the raw/ exclusion inspects path components, neither of
  # which can match a basename.
  if (is.function(.fk) && !.fk(file.path(dir, name))) return(invisible(p))
  p <- fs_caption(p)                                 # scale statement (raw pass by default)
  pdf_dev <- "pdf"                                   # base pdf: core Helvetica, ASCII-safe
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h,
         device = pdf_dev, bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 600, bg = "white")
  message(sprintf("  saved %-34s  %.1f x %.1f in", name, w, h))
  invisible(p)
}

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
fs_apply_lfo_scale(lfo)
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
# model_selection.json rather than recomputed here (see .pick_best_bayes below), so the figures
# cannot track a different model from the pipeline, and the log-score axis, the eligibility
# gates and the tie-break are whatever best_invasion_model() actually applied.
.pick_best_bayes <- function(rs, lfo, ev, fallback) {
  sel_path <- file.path(OUT_DIR, "key_outputs", "model_selection.json")
  m <- NULL
  if (file.exists(sel_path)) {
    sel <- tryCatch(jsonlite::fromJSON(sel_path, simplifyVector = TRUE), error = function(e) NULL)
    m <- tryCatch(sel$featured$bayesian$method, error = function(e) NULL)
    if (is.null(m) || !length(m) || is.na(m[1]))
      m <- tryCatch(sel$featured$headline$method, error = function(e) NULL)
  }
  if (!is.null(m) && length(m) && !is.na(m[1])) return(as.character(m[1]))
  warning("[figures] model_selection.json is unreadable; falling back to the configured ",
          "featured model. The figures may not track the pipeline's pick.", call. = FALSE)
  lfo_methods   <- unique(lfo$method)
  bayes_singles <- unique(lfo_methods[grepl("^Bayes", lfo_methods) & !grepl("-ens-", lfo_methods)])
  if (fallback %in% lfo_methods) return(fallback)
  if (length(bayes_singles)) return(bayes_singles[1])
  fallback
}
BEST_BAYES <- .pick_best_bayes(rs, lfo, ev, BEST_BAYES_FALLBACK)
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
# PRIORITISATION (detection) curve — reconstruction of the RHS of
# bayes_model_performance_figure3: share of true invasions caught when the top-K
# highest-risk zones are monitored each week, pooled over LFO folds, with a
# random-targeting reference line and the naive epicentre-inflow baseline overlaid.
# Rebuilt from the saved LFO + primary mobility matrix + WorldPop spine (no re-fitting);
# mirrors compute_detection_curve() / naive_epicentre_inflow_scores() exactly.
# -----------------------------------------------------------------------------
NAIVE_LBL <- "Naive-epicentre-inflow"
# CALL THE PIPELINE'S SCORER, do not re-derive it. This block used to hand-roll the
# population-weighted epicentre-outflow score: it took `zall` from the WorldPop file rather
# than the canonical zone spine, and it omitted the alias harmonisation that
# naive_epicentre_inflow_scores() applies — so a non-canonical epicentre spelling (the
# pre-2026-07 "Mongbalu" for "Mongbwalu") silently dropped an origin here and changed the
# population weighting, while leaving every other figure's baseline untouched. It was the
# fourth independent copy of one quantity.
#
# NOTE ON SCOPE: this panel's "Naive-epicentre-inflow" is a DIFFERENT series from the three
# Baseline-* structural nulls the manuscript figures publish (gravity M4 / Flowminder cohort /
# OSRM travel time, see run_all.R). It is the primary kernel's own epicentre inflow, kept here
# as this suite's long-standing comparator; the manuscript set is pinned to specific matrices
# on purpose and does not follow MOBILITY_PRIMARY.
NAIVE_SCORES <- tryCatch({
  if (!exists("naive_epicentre_inflow_scores", mode = "function"))
    source(file.path(HERE, "spatiotemporal", "20_forecast_detail.R"))
  stopifnot(exists("MOBILITY_PRIMARY"))
  .pk  <- MOBILITY_PRIMARY
  .pkf <- file.path(OUT, "mobility", sprintf("mobility_%s.rds", .pk))
  if (!file.exists(.pkf))
    .pkf <- file.path(OUT, "mobility", sprintf("mobility_%s.rds", sub("-fill$", "", .pk)))
  W   <- readRDS(.pkf)
  naive_epicentre_inflow_scores(W, EPI_ZONES, load_population(), rownames(W))
}, error = function(e) { message("[figures] naive baseline unavailable (", conditionMessage(e),
                                 "); prioritisation panel drawn without it."); NULL })

# Detection curve for one method's rows: recall (share of invasions caught) + the random
# reference, pooled over folds, for a weekly budget of K zones. Mirrors compute_detection_curve().
.det_curve <- function(d, ks = 1:25) {
  if (!nrow(d)) return(NULL)
  pf   <- d %>% group_by(fold_id) %>% mutate(rk = rank(-p_invasion, ties.method = "max")) %>% ungroup()
  tot  <- sum(pf$is_new_invasion, na.rm = TRUE); if (tot == 0) return(NULL)
  natr <- pf %>% count(fold_id) %>% pull(n) %>% mean()
  purrr::map_dfr(ks, function(k) {
    tp <- sum(pf$is_new_invasion[pf$rk <= k], na.rm = TRUE)
    tibble(k = k, recall = tp / tot, recall_random = pmin(k / natr, 1))
  })
}

# Combined featured-model + naive curves (+ shared random reference) at a given horizon.
prioritisation_curves <- function(horizon = 1L, ks = 1:25) {
  feat <- lfo %>% filter(method == BEST_BAYES, horizon == !!horizon, is.finite(p_invasion))
  cf <- .det_curve(feat, ks); if (is.null(cf)) return(NULL)
  cf$method <- BEST_BAYES
  curves <- cf
  if (!is.null(NAIVE_SCORES)) {
    refn <- feat %>% mutate(p_invasion = as.numeric(NAIVE_SCORES[match(health_zone, names(NAIVE_SCORES))]))
    cn <- .det_curve(refn, ks)
    if (!is.null(cn)) { cn$method <- NAIVE_LBL; curves <- bind_rows(cf, cn) }
  }
  list(curves = curves, rnd = cf %>% distinct(k, recall_random))
}

# province-aware join of per-zone data -> shapefile (disambiguates duplicate names)
join_map <- function(dat) {
  d <- dat %>% mutate(.key = tolower(trimws(health_zone)),
                      .prov = as.character(province))
  m <- shp %>% left_join(d, by = c(".key", ".prov"))
  # fall back to key-only for any rows that failed the province match
  miss <- m %>% st_drop_geometry() %>% summarise(f = mean(is.na(health_zone))) %>% pull(f)
  if (miss > 0.5) {
    d2 <- dat %>% mutate(.key = tolower(trimws(health_zone)))
    m <- shp %>% left_join(d2 %>% dplyr::select(-any_of(".prov")), by = ".key")
  }
  m
}

# -----------------------------------------------------------------------------
# 3. FIGURE 1 DATA  (onset-dated confirmed cases; pipeline-consistent)
# -----------------------------------------------------------------------------
# The line list is the PIPELINE's. load_linelist() resolves the latest processed export,
# reconciles it to the INSP sitrep cumulative (appending the shortfall as confirmed rows
# for zones the line list has not caught up on), and imputes a missing onset by DRAWING
# from the fitted onset->sample delay (epidist marginal-truncation-corrected gamma),
# seeded from RANDOM_SEED. This figure therefore shows exactly the case set, zone set and
# onset dates the invasion model is fitted and scored on.
#
# This script used to read the processed CSV directly and shift every missing onset back
# by a single ROUNDED MEAN delay. A fixed shift is a different estimator: it collapses a
# right-skewed gamma (mean 7.67 d, SD 8.50 d) onto its mean, so the panel's arrival dates
# and cumulative curve matched no table in the pipeline, and the caption had to carry a
# standing disclaimer saying so. Both are now gone.
ll <- load_linelist()
wk_floor <- floor_date(OUTBREAK_START, "week", week_start = 1)
stopifnot(all(c("health_zone", "province", "date_index", "confirmed", "onset_usable") %in% names(ll)))

conf <- ll %>% filter(confirmed %in% TRUE, !is.na(date_index), date_index >= wk_floor,
                      !is.na(health_zone))
# Share of confirmed cases whose onset was imputed (for the caption; computed, not hardcoded).
PCT_IMPUTED <- round(100 * mean(!conf$onset_usable, na.rm = TRUE), 1)
# The delay the caption quotes is READ from the shared resolver — the same object
# load_linelist() drew from — not a second estimate computed here.
ONSET_SAMPLE_MEAN <- effective_onset_sample_delay()$mean
stopifnot(is.finite(ONSET_SAMPLE_MEAN), ONSET_SAMPLE_MEAN > 0)

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

# --- 1A mobility overlay: latest Flowminder short-trip OUTFLOW from the epicentre ---
# The short-trip cohort is pooled over the epicentre zones (from EPICENTRE_ZONES); we anchor the
# outflow fan at the principal epicentre and draw the strongest N destination shares as
# width-encoded, DIRECTIONAL arcs (arrowheads point epicentre -> destination). Data are the same
# snapshots the M1/M8 mobility matrices consume (03_mobility_matrices.R).
FLOW_ST_DIR <- file.path(dirname(dirname(SHP_PATH)), "flowminder_short_trips", "processed")
EPI_ORIGIN  <- EPI_ZONES[1]     # principal epicentre; anchor for the outflow fan
ITURI_PROV  <- PROV_INT[1]      # epicentre province (Ituri) — for the "share staying within" stat
N_FLOWS     <- 12L              # strongest N destinations (keeps the map legible)
FLOW_COL    <- "#EE7733"        # Tol 'vibrant' orange; warm accent vs the cool mako fill
flows_df <- NULL; flows_all <- NULL; epi_pt <- NULL; latest_flow_date <- NA; ITURI_PCT <- NA_real_
try({
  ffs <- list.files(FLOW_ST_DIR,
                    pattern = "^flowminder_short_trips__outflow_[0-9]{8}__static\\.csv$",
                    full.names = TRUE)
  stopifnot(length(ffs) > 0)
  tags   <- as.integer(sub(".*outflow_([0-9]{8})__.*", "\\1", basename(ffs)))
  latest <- ffs[which.max(tags)]                       # the most recent snapshot
  st <- suppressWarnings(read_csv(latest, show_col_types = FALSE))
  names(st)[1:2] <- c("nom", "prop")
  st$key <- tolower(trimws(st$nom))
  st <- st[st$key != tolower(EPI_ORIGIN) & is.finite(st$prop), ]
  # zone centroids (point-on-surface = always inside the polygon) + province, for arc geometry
  # and the destination-province breakdown.
  pts <- suppressWarnings(sf::st_point_on_surface(sf::st_geometry(shp)))
  cc  <- sf::st_coordinates(pts)
  cent <- data.frame(key = shp$.key, prov = as.character(shp$.prov),
                     cx = cc[, 1], cy = cc[, 2], stringsAsFactors = FALSE)
  epi_row <- cent[cent$key == tolower(EPI_ORIGIN), ][1, ]
  # a few zone names recur across provinces; keep the instance nearest the epicentre
  cent$d <- (cent$cx - epi_row$cx)^2 + (cent$cy - epi_row$cy)^2
  cent <- cent[order(cent$key, cent$d), ]
  cent <- cent[!duplicated(cent$key), c("key", "prov", "cx", "cy")]
  fl_all <- merge(st, cent, by = "key")
  fl_all$x <- epi_row$cx; fl_all$y <- epi_row$cy; fl_all$xend <- fl_all$cx; fl_all$yend <- fl_all$cy
  fl_all <- fl_all[order(fl_all$prop), ]            # small first so big arcs draw on top
  fl <- head(fl_all[order(-fl_all$prop), ], N_FLOWS)
  flows_df <- fl; flows_all <- fl_all; epi_pt <- epi_row
  latest_flow_date <- as.Date(as.character(max(tags)), "%Y%m%d")
  # COMPUTED (not hard-coded): share of non-self short-trip outflow whose destination is in Ituri.
  ITURI_PCT <- round(100 * sum(fl_all$prop[fl_all$prov == ITURI_PROV], na.rm = TRUE) /
                     sum(fl_all$prop, na.rm = TRUE), 0)
  message(sprintf("[fig1] epicentre flows: %d arcs from %s, snapshot %s (%.0f%% -> %s; top dest: %s)",
                  nrow(fl), EPI_ORIGIN, latest_flow_date, ITURI_PCT, ITURI_PROV,
                  paste(head(fl$nom, 4), collapse = ", ")))
}, silent = TRUE)
if (is.null(flows_df)) message("[fig1] WARNING: no Flowminder outflow snapshot found; F1A drawn without flows.")

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
      name = "Epicentre outflow\n(% of short trips)", range = c(0.3, 2.6),
      breaks = c(5, 10, 20), limits = c(0, NA),
      guide = guide_legend(order = 2, override.aes = list(colour = FLOW_COL, alpha = 0.9))) +
    coord_sf(xlim = ZOOM_X, ylim = ZOOM_Y, expand = FALSE) +
    theme_map(11.5) +
    theme(legend.position = "right",
          plot.title    = element_blank(),   # no embedded titles — see FIGURE_CAPTIONS.md
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
  xsc  <- scale_x_date(limits = xlim, date_breaks = "2 weeks",
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
    labs(y = "Health zones with\n>=1 confirmed case",
         x = "Symptom-onset date", caption =
         # The panel and the pipeline now share one estimator, so the caption states it once.
         sprintf(paste0("Onset imputed for the %s%% of confirmed cases lacking a usable onset, ",
                        "by drawing the onset-to-sample delay per record from the fitted, ",
                        "truncation-corrected distribution (gamma, mean %.1f d) — the same draw ",
                        "the invasion model is fitted on."), PCT_IMPUTED, ONSET_SAMPLE_MEAN)) +
    # theme_pub() blanks plot.caption, so this disclosure was COMPUTED AND DISCARDED — a main-text
    # figure in which a quarter of the plotted onsets are imputed said nothing about it. Re-enable
    # the caption for this panel only.
    theme_pub(16.5) +
    theme(plot.caption = element_text(size = 7.5, colour = MUTED, hjust = 0)) +
    theme(plot.margin = margin(0,8,6,6),
                            axis.text.x = element_text(size = 13.5))
  p1b <- top / bot + plot_layout(heights = c(1, 0.78))

  save_dual(p1a, "F1A_cumulative_case_map", 6.4, 5.4)
  save_dual(p1b, "F1B_onset_timeseries",    6.2, 5.4)

  # Figure 1 panel A is the national all-flows case map (full epicentre mobility
  # reach); fall back to the zoomed top-N panel if no flow snapshot is available.
  p1a_main <- make_f1a_allflows_panel()
  if (is.null(p1a_main)) p1a_main <- p1a

  # No embedded figure title/subtitle — only the A/B panel tags. Descriptive text is in
  # FIGURE_CAPTIONS.md.
  fig1 <- (wrap_elements(full = p1a_main) | wrap_elements(full = p1b)) +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(tag_levels = "A")
  save_dual(fig1, "Figure1", 11.8, 5.0, dir = FIG_DIR)
  invisible(fig1)
}

# =============================================================================
# FIGURE 1A — ALL epicentre outflows (national extent)
# Companion to the zoomed top-N panel: shows the COMPLETE short-trip reach from the epicentre with
# DIRECTIONAL arcs (arrowheads point epicentre -> destination). Most trips stay within the epicentre
# province (the dense near fan; the share is COMPUTED and shown in the caption); the thin threads are
# the long tail reaching Kinshasa and the far provinces. Text is enlarged for main-figure legibility.
# =============================================================================
make_f1a_allflows_panel <- function() {
  if (is.null(flows_all) || nrow(flows_all) == 0) return(NULL)
  mp <- join_map(zone_cum) %>%
    mutate(cases = ifelse(is.na(cases), 0, cases),
           bin = cut(cases, breaks = c(-1, 0, 9, 49, 99, 499, Inf),
                     labels = c("0","1-9","10-49","50-99","100-499","500+")))
  seq6 <- viridisLite::mako(6, begin = 0.92, end = 0.12)
  names(seq6) <- c("0","1-9","10-49","50-99","100-499","500+"); seq6["0"] <- NA_FILL

  # label only the epicentre (destination place-names dropped for legibility; the far-western
  # Kinshasa terminus is still annotated separately below)
  lab_all <- data.frame(x = epi_pt$cx, y = epi_pt$cy,
                        label = paste0(EPI_ORIGIN, " (epicentre)"),
                        face = "bold", col = INK, stringsAsFactors = FALSE)
  n_reach <- nrow(flows_all)
  # DATA-DERIVED far-western terminus label: the strongest Kinshasa-province destination (if the
  # fan reaches it), placed at its OWN centroid — no hard-coded coordinates.
  kin <- flows_all[!is.na(flows_all$prov) & flows_all$prov == "Kinshasa", , drop = FALSE]
  kin <- if (nrow(kin)) kin[which.max(kin$prop), , drop = FALSE] else NULL

  arr <- grid::arrow(length = unit(0.05, "in"), type = "closed", angle = 18)
  # cap_share was built here and NEVER USED: this panel deliberately blanks plot.title,
  # plot.subtitle and plot.caption ("no embedded titles — see FIGURE_CAPTIONS.md", below), so
  # the string was formatted and discarded on every run. The convention is right — the captions
  # live outside the figure — so the dead line is removed rather than rendered. The number
  # itself (ITURI_PCT, the share of the epicentre's short-trip outflow staying within Ituri)
  # is still COMPUTED and LOGGED where it is derived, and is available for the caption text.
  p <- ggplot(mp) +
    geom_sf(aes(fill = bin), colour = "white", linewidth = 0.06) +
    # every destination arc; width = outflow share, DIRECTIONAL (arrow epicentre -> destination)
    geom_curve(data = flows_all,
               aes(x = x, y = y, xend = xend, yend = yend, linewidth = prop),
               curvature = 0.16, angle = 90, ncp = 12, colour = FLOW_COL,
               alpha = 0.6, lineend = "round", arrow = arr) +
    geom_point(data = epi_pt, aes(cx, cy), shape = 21, size = 3.6,
               fill = "white", colour = FLOW_COL, stroke = 1.3) +
    ggrepel::geom_text_repel(
      data = lab_all, aes(x, y, label = label, fontface = face, colour = col),
      size = 5.2, min.segment.length = 0, segment.colour = FAINT,
      segment.size = 0.25, box.padding = 0.5, max.overlaps = 40, seed = 1) +
    { if (!is.null(kin)) annotate("text", x = kin$xend, y = kin$yend - 0.35,
             label = "Kinshasa", size = 4.8, colour = MUTED, fontface = "italic") } +
    scale_colour_identity() +
    scale_fill_manual(values = seq6, na.value = NA_FILL, drop = FALSE,
                      name = "Cumulative\nconfirmed cases",
                      guide = guide_legend(reverse = TRUE, order = 1,
                                           override.aes = list(colour = "white"))) +
    scale_linewidth_continuous(
      name = "Epicentre outflow\n(% of short trips)", range = c(0.15, 2.9),
      breaks = c(1, 5, 10, 20), limits = c(0, NA),
      guide = guide_legend(order = 2, override.aes = list(colour = FLOW_COL, alpha = 0.9))) +
    coord_sf(expand = FALSE) +
    theme_map(17) +
    theme(legend.position = "right",
          plot.title    = element_blank(),   # no embedded titles — see FIGURE_CAPTIONS.md
          plot.subtitle = element_blank(),
          plot.caption  = element_blank(),
          panel.border  = element_rect(fill = NA, colour = GRID, linewidth = 0.4))

  p
}

# thin wrapper: build the national all-flows panel and save it standalone
build_f1a_allflows <- function() {
  p <- make_f1a_allflows_panel()
  if (is.null(p)) {
    message("[fig1] all-flows variant skipped (no snapshot)."); return(invisible(NULL))
  }
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
  # the hand-written table this replaced was missing 12 of the 57 cross-validated models, so
  # raw codes such as "Bayes-M8-dist-split-geo" were printed on the one figure whose entire
  # purpose is to avoid them. model_pretty_label() still warns, naming any code it cannot
  # decode, so a genuinely new token shows up in the run log rather than only in the PDF.
  .use_labels <- isTRUE(model_labels) ||
                 (!is.null(model_labels) && !is.logical(model_labels) && length(model_labels) > 0)
  relab <- function(v) {
    v <- as.character(v)
    if (!.use_labels) return(v)
    model_pretty_label(v)
  }
  y_relabel <- if (!.use_labels) NULL else scale_y_discrete(labels = function(v) relab(v))
  yt_size   <- if (!.use_labels) 11 else 8.2

  # --- 2A discrimination: AUC-PR skill for the Bayesian grid, h1 & h2 ---
  # PANEL-A/C MODEL SET. The cross-validated Bayesian grid is 54 models. Panels A and C each
  # occupy ~4.2 in of canvas height inside the composite, so drawing all 54 gives ~0.08 in per
  # row against 8.2pt labels (~0.11 in) — the y-axis rendered as a solid block of overprinted
  # text and neither panel could be read. Show the leading N_MODELS by 1-week AUC-PR skill
  # (plus the featured model, wherever it ranks) at a legible row pitch; the COMPLETE 54-model
  # grid, on both metrics, is published separately as bayes_discrimination_summary_h{1,2}.
  # A and C use the SAME model set so the two panels can be read against each other.
  N_MODELS <- 20L
  # THE CANDIDATE GRID ONLY. Panels A and C are a comparison BETWEEN models, so they must hold
  # only models that competed: the cross-validated specifications. Two kinds of row were being
  # drawn alongside them and should not be.
  #   * SENSITIVITY ARMS (-gtshort/-gtlong, -tv*) are refits of ONE specification at a changed
  #     assumption. Ranking them against the field invites the reader to treat "short
  #     generation time" as a rival model, which is exactly what composing the grid at a single
  #     anchor is designed to prevent, and they are barred from selection for the same reason.
  #   * ENSEMBLES (-ens-) are combinations OF the grid, not members of it, so they cannot be
  #     read as one more kernel-plus-covariate choice.
  # Both remain fully scored and published: the arms in the sensitivity table, the ensembles in
  # the evaluation table and in bayes_discrimination_summary, which is the all-methods figure.
  .sens_re <- get0("INVASION_SELECTION_EXCLUDE", ifnotfound = "-(gtshort|gtlong|tv[a-z0-9]+)$")
  # THE STRUCTURAL BASELINES BELONG IN THIS PANEL. Filtering to "^Bayes" dropped them by
  # construction, so the discrimination panel showed the fitted grid ranked only against
  # itself -- the reader could not see how much of the skill is the model rather than the
  # geography. Derived from the evaluation table rather than named here, so a change to the
  # baseline set reaches the figure instead of silently leaving a hard-coded name behind.
  .eval_methods   <- unique(ev$method)
  .eval_baselines <- sort(setdiff(.eval_methods, grep("^Bayes", .eval_methods, value = TRUE)))
  message(sprintf("[fig2] structural baselines in the evaluation table: %s",
                  if (length(.eval_baselines)) paste(.eval_baselines, collapse = ", ") else "NONE"))
  .is_primary <- function(m) grepl("^Bayes", m) & !grepl(.sens_re, m) & !grepl("-ens-", m)
  d <- ev %>% filter(horizon %in% c(1,2),
                     .is_primary(method) | method %in% .eval_baselines) %>%
    mutate(lo = pmax(auc_pr_lo,0)/base_rate, hi = auc_pr_hi/base_rate,
           family = factor(ifelse(method %in% .eval_baselines,
                                  "Structural baseline", "Primary specification"),
                           levels = c("Primary specification", "Structural baseline")))
  if (!nrow(d)) stop("[fig2] no candidate specifications left after excluding arms and ensembles")
  # Rank and cap the PRIMARY specifications only; every baseline present is always kept, so
  # the comparator set cannot be truncated away by a grid that grew.
  .prim <- d %>% filter(family == "Primary specification")
  keep_models <- .prim %>% filter(horizon == 1, is.finite(auc_pr_skill)) %>%
    arrange(desc(auc_pr_skill)) %>% pull(method)
  if (!length(keep_models))
    keep_models <- .prim %>% filter(is.finite(auc_pr_skill)) %>%
      arrange(desc(auc_pr_skill)) %>% pull(method) %>% unique()
  n_grid  <- dplyr::n_distinct(.prim$method)      # primary specifications cross-validated
  n_basel <- length(.eval_baselines)
  if (n_grid > N_MODELS)
    message(sprintf("[fig2] %d primary specifications exceed the N_MODELS=%d cap; showing the leading %d.",
                    n_grid, N_MODELS, N_MODELS))
  keep_models <- unique(c(head(keep_models, N_MODELS), intersect(BEST_BAYES, .prim$method),
                          .eval_baselines))
  d <- d %>% filter(method %in% keep_models)
  ord <- d %>% filter(horizon == 1) %>% arrange(auc_pr_skill) %>% pull(method)
  if (!length(ord)) ord <- d %>% arrange(auc_pr_skill) %>% pull(method) %>% unique()
  # COMPLETE the level set. ord is built from the h=1 rows only, so any model present at h=2
  # but not at h=1 (the featured model is added to keep_models independently of its h=1 row)
  # would become a factor NA and be dropped from the panel silently, with only a ggplot
  # "removed rows" warning that nobody reads in a batch log.
  ord <- unique(c(setdiff(unique(d$method), ord), ord))
  d <- d %>% mutate(method = factor(method, levels = ord),
                    hz = factor(horizon))
  p2a <- ggplot(d, aes(auc_pr_skill, method, colour = hz)) +
    geom_vline(xintercept = 1, linetype = "22", colour = FAINT) +
    geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .55),
                   linewidth = .5, alpha = .55) +
    geom_point(aes(shape = family), position = position_dodge(width = .55), size = 1.9) +
    scale_colour_manual(values = HZ_COL, name = "Horizon",
                        labels = c("1 week","2 weeks")) +
    scale_shape_manual(values = c("Primary specification" = 16, "Structural baseline" = 17),
                       name = NULL, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = "Discrimination of the Bayesian invasion models",
         subtitle = "AUC-PR skill = average precision / base rate  (1 = no skill; higher = better)",
         x = "AUC-PR skill (x base rate)", y = NULL,
         caption = sprintf(paste0("Points = pooled leave-future-out estimate; bars = 90%% zone-cluster ",
                                  "bootstrap CI. %s of the %d primary specifications (circles), ",
                                  "ordered by 1-week skill, against %s (triangles). Sensitivity arms ",
                                  "(generation-time and time-varying refits) and the ensembles are ",
                                  "also cross-validated but are excluded here; all methods appear in ",
                                  "bayes_discrimination_summary."),
                           if (dplyr::n_distinct(.prim$method[.prim$method %in% keep_models]) >= n_grid)
                             "All" else sprintf("The leading %d",
                               dplyr::n_distinct(.prim$method[.prim$method %in% keep_models])),
                           n_grid,
                           if (n_basel) sprintf("the %d structural baseline%s", n_basel,
                                                if (n_basel == 1L) "" else "s")
                           else "no structural baseline (none in the evaluation table)")) +
    theme_pub(13.5) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = yt_size, colour = INK)) +
    y_relabel

  # --- 2B prioritisation: invasions caught vs random targeting (RHS of the
  #     bayes_model_performance_figure3 panel). Featured Bayesian model + naive
  #     epicentre-inflow baseline + random reference, 1 week ahead. Replaces the
  #     former calibration panel; no embedded title/subtitle. ---
  pc <- prioritisation_curves(horizon = 1L)
  MODEL_COL <- setNames(c("#0072B2", "#D55E00"), c(BEST_BAYES, NAIVE_LBL))
  MODEL_LBL <- setNames(c(relab(BEST_BAYES), "Naive epicentre-inflow"), c(BEST_BAYES, NAIVE_LBL))
  curves_b  <- pc$curves %>% mutate(method = factor(method, levels = intersect(names(MODEL_COL),
                                                                               unique(method))))
  p2b <- ggplot(curves_b, aes(k, recall, colour = method)) +
    geom_line(data = pc$rnd, aes(k, recall_random), linetype = "22", colour = FAINT,
              linewidth = 0.7, inherit.aes = FALSE) +
    geom_line(linewidth = 0.9) + geom_point(size = 1) +
    annotate("text", x = max(pc$rnd$k) * 0.6, y = max(pc$rnd$recall_random) * 0.7 + 0.045,
             label = "random targeting", colour = MUTED, size = 3.8, angle = 8) +
    scale_colour_manual(values = MODEL_COL, labels = MODEL_LBL, name = NULL,
                        breaks = levels(curves_b$method)) +
    scale_y_continuous(labels = percent_format(1), limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    labs(title = "Real-time prioritisation skill",
         subtitle = "Invasions caught vs a random watch-list of the same size (1 week ahead)",
         # STATE THE HORIZON. theme_pub() blanks plot.subtitle, which was the only text naming
         # this panel as the 1-week task — while panel E is the 2-week task. A reader saw one
         # panel labelled "2 weeks ahead" and the rest unlabelled.
         x = "Zones actively monitored each round (K), 1 week ahead",
         y = "Share of true invasions caught") +
    theme_pub(13.5)

  # (Panels 2C prioritisation and 2D predicted-vs-observed have been dropped from
  #  Figure 2 at the author's request; the figure now carries A, B, E, F only.)

  # --- 2E ranking accuracy: mean rank of the invaded zones, Bayesian grid, h1 & h2 ---
  # Lower = the truly-invaded zones sat nearer the TOP of the model's watch-list (a direct
  # operational readout, complementary to AUC-PR). Sourced from invasion_evaluation.csv
  # (mean_rank_of_truth); the two horizons are joined per model so the h1->h2 shift is visible.
  # Same model set as panel A (keep_models), for the same legibility reason and so the two
  # panels are read against one another rather than against two different subsets.
  # Membership is keep_models, which now carries the structural baselines as well as the
  # primary specifications, so panels A and C describe the SAME set of methods.
  mr <- ev %>% filter(horizon %in% c(1,2), is.finite(mean_rank_of_truth),
                      method %in% keep_models) %>%
    mutate(family = factor(ifelse(method %in% .eval_baselines,
                                  "Structural baseline", "Primary specification"),
                           levels = c("Primary specification", "Structural baseline")))
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
    scale_shape_manual(values = c("Primary specification" = 16, "Structural baseline" = 17),
                       name = NULL, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.08))) +
    labs(title = "Ranking accuracy for the invaded zones",
         subtitle = "Mean rank of the truly-invaded zones among at-risk zones (lower = nearer the top)",
         x = "Mean rank of invaded zones (lower is better)", y = NULL,
         caption = sprintf(paste0("Rank of every invaded zone (tie-averaged), meaned per fold then ",
                                  "over folds. The same %d methods as panel A (%d primary ",
                                  "specifications, circles; %d structural baselines, triangles)."),
                           dplyr::n_distinct(mr$method),
                           dplyr::n_distinct(mr$method[!mr$method %in% .eval_baselines]),
                           dplyr::n_distinct(mr$method[mr$method %in% .eval_baselines]))) +
    theme_pub(13.5) +
    theme(panel.grid.major.y = element_blank(),
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
  # and the x-axis drew a dated break at EVERY fold cutoff. Track the top 15 (the top-K
  # convention used by every other panel in the suite) and thin the axis labels below.
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
      mutate(region = ifelse(!is.na(province) & province %in% PROV_INT, province, "Other"),
             rank_disp = pmin(rank, RANK_FLOOR))
    reg_levels <- c(intersect(PROV_INT, unique(traj$region)), "Other")
    traj$region <- factor(traj$region, levels = reg_levels)
    end_lab <- traj %>% filter(cutoff == last_cut)
    x_pad <- as.numeric(diff(range(traj$cutoff))) * 0.06 + 1
    p2f <- ggplot(traj, aes(cutoff, rank_disp, group = health_zone, colour = region)) +
      geom_line(linewidth = 1.15, alpha = 0.5, lineend = "round") +
      geom_point(size = 1.7) +
      ggrepel::geom_text_repel(data = end_lab, aes(label = health_zone),
        size = 3.1, direction = "y", hjust = 0, nudge_x = x_pad, box.padding = 0.12,
        segment.size = 0.2, segment.colour = FAINT, min.segment.length = 0,
        max.overlaps = 40, seed = 1) +
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
      scale_colour_manual(values = PROV_COL, name = "Province", breaks = reg_levels,
                          na.value = "#7A7A7A") +
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
  # spaced rounds, first and last always kept, keep the trajectory legible. The standalone
  # F_topk_folds_h1/h2 panels below still show every round.
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
    "Model = the spatial mobility kernel driving invasion risk: gravity / composite gravity / ",
    "radiation-composite, their Flowminder cohort and relocation-OD counterparts, and the ",
    "all-kernel consensus over these. ",
    paste(Filter(Negate(is.null), list(
      .clause("covariates", paste0(
        "\"+ covariates\" = geographic & social covariates modulate the import-to-invasion rate ",
        "(base models use a single constant rate)")),
      .clause("road-km", "\"road-km\" = road-distance rather than travel-time deterrence"),
      .clause("source-filled", paste0(
        "\"source-cell fill\" = destinations the mobility source could not observe are taken ",
        "from the base kernel rather than assumed zero")),
      .clause("suspected", "\"+ suspected\" = suspected-case leading indicators"),
      .clause("ensemble", "\"ensemble\" = mean / median across models"))),
      collapse = "; "),
    ". Sensitivity arms and ensembles are not shown here; see bayes_discrimination_summary."),
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
qtile <- function(x) {                                  # rank-based quartile 1..4
  out <- rep(NA_integer_, length(x)); ok <- is.finite(x)
  r <- rank(x[ok], ties.method = "average")/sum(ok)
  out[ok] <- as.integer(cut(r, c(-Inf,.25,.5,.75,Inf), labels = 1:4)); out
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
                 bx = qtile(rr01_nat), by = qtile(V),
                 bikey = ifelse(is.na(bx)|is.na(by), NA, paste0(bx,"-",by)))
  }
  aff_overlay <- function(mp) geom_sf(data = mp %>% filter(was_active_before %in% TRUE),
                                      fill = AFFECTED_FILL, colour = "white", linewidth = 0.08)

  one_map <- function(h, fillvar, scale_fn, title) {
    mp <- join_map(prep(h))
    ggplot(mp) +
      geom_sf(aes(fill = .data[[fillvar]]), colour = "white", linewidth = 0.08) +
      aff_overlay(mp) + scale_fn + ggtitle(title) + theme_map(13.5)
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
  wmax <- rs %>% filter(!(was_active_before %in% TRUE)) %>%
    mutate(w = p_case_hi - p_case_lo) %>% pull(w) %>% quantile(.99, na.rm=TRUE)
  sc_w <- scale_fill_viridis_c(option="rocket", direction=-1, limits=c(0, max(0.05,wmax)),
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
  p3c <- ggplot(topz, aes(p_case_invasion, zone, colour = hz)) +
    geom_linerange(aes(xmin = p_case_lo, xmax = p_case_hi), linewidth = 0.6,
                   position = pd_c, orientation = "y") +
    geom_point(size = 1.9, position = pd_c) +
    scale_colour_manual(values = HZ_COL, name = "Horizon", labels = c("1 week","2 weeks")) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(0, NA),
                       expand = expansion(mult = c(0.01, 0.08))) +
    labs(title = "Highest-risk zones: 1- vs 2-week invasion probability",
         subtitle = "Top at-risk zones by 1-week risk; posterior mean and 90% CrI at both horizons",
         x = "Invasion probability", y = NULL) +
    theme_pub(13.5) + theme(panel.grid.major.y = element_blank(),
                        axis.text.y = element_text(colour = INK, size = 11.5))

  # --- 3D preparedness-priority scatter (1 week ahead): invasion likelihood x vulnerability,
  # sized by the composite priority, coloured by province (reconstruct of bayes_priority_scatter_h1). ---
  psc <- rs %>% filter(horizon == 1, !(was_active_before %in% TRUE),
                       is.finite(V), is.finite(rr01_nat)) %>%
    mutate(region = ifelse(!is.na(province) & province %in% PROV_INT, province, "Other"))
  psc$region <- factor(psc$region, levels = c(intersect(PROV_INT, unique(psc$region)), "Other"))
  lab_ps <- psc %>% arrange(desc(priority)) %>% head(12)
  p3d <- ggplot(psc, aes(V, rr01_nat)) +
    geom_point(aes(size = priority, colour = region), alpha = 0.78) +
    ggrepel::geom_text_repel(data = lab_ps, aes(label = health_zone), size = 3.1, colour = INK,
      box.padding = 0.3, max.overlaps = 30, seed = 1, min.segment.length = 0,
      segment.colour = FAINT, segment.size = 0.2) +
    scale_colour_manual(values = PROV_COL, name = "Province", breaks = levels(psc$region),
                        na.value = "#7A7A7A") +
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
      ggtitle(hz_lab(h)) + theme_map(13.5)
  }
  leg_df <- expand.grid(bx=1:4, by=1:4) %>% mutate(bikey = paste0(bx,"-",by))
  legend_biv <- ggplot(leg_df, aes(bx, by, fill = bikey)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = BIP, guide = "none") + coord_fixed() +
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

  # No embedded figure title/subtitle — only the A-E panel tags. Descriptive text (incl. the
  # featured model and the left=1wk / right=2wk map convention) is in FIGURE_CAPTIONS.md.
  fig3 <- (wrap_elements(full = p3a) / wrap_elements(full = p3b) /
           (wrap_elements(full = p3c) | wrap_elements(full = p3d)) /
           wrap_elements(full = p3e)) +
    plot_layout(heights = c(1, 1, 1.1, 1)) +
    plot_annotation(tag_levels = "A")
  save_dual(fig3, "Figure3", 10.0, 14.7, dir = FIG_DIR)
  invisible(fig3)
}

# helpers for within-facet reordering (avoid tidytext dependency)
reorder_within <- function(x, by, within) {
  factor(paste(x, within, sep = "___"),
         levels = unique(paste(x, within, sep="___"))[order(within, by)])
}
tidytext_scale_y <- function() scale_y_discrete(labels = function(z) sub("___.*$","",z))

# =============================================================================
# FIGURE: top predicted invasion probabilities per fold, coloured by outcome
# Companion to the diagnostic bayes_lfo_forecast_vs_outcome tile: instead of a
# zone x cutoff heatmap, this separates each leave-future-out fold into its OWN
# facet and draws a bar chart of the highest-risk zones that round, coloured by
# whether that zone actually recorded its first case (invaded) or not. It reads
# as a per-round precision/reliability view: tall bars that are highlighted =
# confident calls that came true; tall grey bars = the round's false alarms.
# =============================================================================
#' @param save    write the standalone F_topk_folds_h<h> panel (PDF+PNG).
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

  # top-K zones by predicted probability WITHIN each fold. Facet header = round + date only.
  d <- d %>% mutate(cutoff = as.Date(cutoff))
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
  # Facet label: round number only when the panel is embedded (the dates are on panel D's
  # axis and in the caption); round + date on the standalone panel, which has room.
  # Round number AND date even when embedded. The justification for dropping the date ("the dates
  # are on panel D's axis and in the caption") was false: panel D shows dates without round
  # numbers, and the only surviving caption is the model-decoding text.
  lab_lk   <- setNames(if (embed) sprintf("R%d - %s", round_no[as.character(fold_ord)],
                                          format(lfo_origin(fold_ord), "%d %b"))
                       else sprintf("Round %d - %s", round_no[as.character(fold_ord)],
                                    format(lfo_origin(fold_ord), "%d %b")),
                       as.character(fold_ord))

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
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.06)), breaks = scales::pretty_breaks(3)) +
    # The horizon must survive `embed = TRUE`. Embedded in Figure 2 this panel drops its title
    # AND subtitle, which were the only text naming it as the 2-week task — while panel B's axis
    # says "each week" and panels A/C carry a 1-week/2-week legend. A reader had no way to tell.
    labs(title = ttl, subtitle = sub,
         x = sprintf("Predicted invasion probability, P(first case), %d week%s ahead",
                     horizon, if (horizon == 1L) "" else "s"), y = NULL) +
    theme_pub(base_sz) +
    # strip.clip = "off" keeps the per-round facet header from being clipped to its narrow strip
    theme(panel.grid.major.y = element_blank(),
          panel.spacing.x = unit(9, "pt"),
          strip.clip = "off",
          axis.text.y = element_text(size = if (embed) 9.2 else 6.6, colour = INK),
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
