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

HERE      <- normalizePath(".")

# Source the pipeline config so paths/constants (outbreak start, shapefile, provinces of interest,
# epicentre zones, onset->sample delay rate, linelist discovery) come from the SINGLE source of
# truth rather than being hard-coded here — otherwise they silently drift from the pipeline (e.g. a
# pinned linelist folder goes stale). Falls back to local literals only if the config cannot load.
.cfg_ok <- tryCatch({ source(file.path(HERE, "/spatiotemporal/00_config.R")); TRUE },
                    error = function(e) { message("[figures] 00_config.R not sourced (",
                                                  conditionMessage(e), "); using local fallbacks."); FALSE })
.cfg <- function(name, default) if (.cfg_ok && exists(name)) get(name) else default

OUT       <- .cfg("OUT_DIR", file.path(HERE, "outputs"))
FIG_DIR   <- file.path(OUT, "key_outputs", "figures")
PANEL_DIR <- file.path(FIG_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

BEST_BAYES_FALLBACK <- "Bayes-M10-med"   # only used if the featured model cannot be derived
OUTBREAK_START <- .cfg("OUTBREAK_START", as.Date("2026-04-30"))
ANALYSIS_DATE  <- .cfg("ANALYSIS_DATE", Sys.Date())        # as-of date; panels run up to here
IMP_DELAY      <- round(1 / .cfg("DELAY_ONSET_SAMPLE_RATE", 0.228))   # onset<-sample imputation (days)
SHP_PATH <- .cfg("SHAPEFILE_PATH",
                 file.path(HERE, "..", "data", "shapefiles", "DRC_Health_zones.shp"))
PROV_INT  <- .cfg("PROVINCES_OF_INTEREST", c("Ituri", "Nord-Kivu", "Haut-Uele"))
EPI_ZONES <- .cfg("EPICENTRE_ZONES", c("Bunia", "Mongbalu", "Rwampara"))

# Resolve the CURRENT processed linelist via latest.json (exactly as load_linelist() does), so the
# figure always tracks the latest snapshot instead of a pinned LINELIST_<date> folder.
LINELIST <- local({
  lj <- .cfg("LINELIST_JSON", NA_character_)
  ld <- .cfg("LINELIST_DIR",
             file.path(HERE, "..", "data", "processed", "dhis2_linelist_processed"))
  if (!is.na(lj) && file.exists(lj)) {
    meta <- tryCatch(jsonlite::fromJSON(lj), error = function(e) NULL)
    if (!is.null(meta$folder))
      return(file.path(ld, meta$folder, "dhis2_processed_linelist.csv"))
  }
  # fallback: newest LINELIST_* folder on disk
  cand <- list.files(ld, pattern = "^LINELIST_", full.names = TRUE)
  if (length(cand)) file.path(sort(cand, decreasing = TRUE)[1], "dhis2_processed_linelist.csv")
  else file.path(ld, "LINELIST_current", "dhis2_processed_linelist.csv")
})
message("[figures] linelist: ", LINELIST)

# -----------------------------------------------------------------------------
# 1. DESIGN SYSTEM  (colourblind-safe; sequential = single perceptual hue)
# -----------------------------------------------------------------------------
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
AFFECTED_FILL <- "grey78"; NA_FILL <- "grey93"
RR_FLOOR <- 1e-3   # pseudo-log floor for the log10 relative-risk map (absorbs exact-zero zones)

# Okabe-Ito (CVD-safe categorical); fixed province identity map
OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
PROV_COL <- c("Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73",
              "Sud-Kivu"="#CC79A7","Other"="#7A7A7A")
# series colours for the evaluation panels
SER_COL <- c("Mean predicted P"="#0072B2","Observed invasion fraction"="#111111",
             "Mean P at invaded zones"="#D55E00")
HZ_COL  <- c("1"="#3B4CC0","2"="#B4413C")            # horizon 1 / 2

# Intuitive, reader-friendly model names for the labelled Figure 2 variant. Keys are the
# pipeline method codes (see 21_bayesian_renewal.R / bayes_default_grid): the mobility kernel
# driving invasion risk (M4=gravity, M8=composite gravity, M9=multi-kernel ensemble,
# M10=radiation-composite, M13=cohort gravity, M14=cohort radiation, M15=combined static,
# M16=cohort+static, M17=all-kernel consensus) and the beta specification (-med = constant beta,
# no covariates; -geo = geographic + social covariates modulate beta; -full = + healthsite density;
# -dist = road-km deterrence; -susp = suspected-case indicators; -short/-long = generation-time
# profile; ens-mean/median = stacked ensemble). "Bayes" is dropped (all models are Bayesian).
# Every model fitted in the default grid (9 kernels x {med, geo}) has an entry so no raw code
# leaks onto the axes; optional families (dist/full/susp/gen-time/ensemble) are covered too.
MODEL_LABELS <- c(
  # core mobility kernels, constant beta (no covariates)
  "Bayes-M4-med"       = "Gravity",
  "Bayes-M8-med"       = "Composite gravity",
  "Bayes-M9-med"       = "Multi-kernel ensemble",
  "Bayes-M10-med"      = "Radiation-composite",
  "Bayes-M13-med"      = "Cohort gravity",
  "Bayes-M14-med"      = "Cohort radiation",
  "Bayes-M15-med"      = "Combined static",
  "Bayes-M16-med"      = "Cohort + directed OD",
  "Bayes-M17-med"      = "All-kernel consensus",
  # same kernels with geographic + social covariates modulating beta
  "Bayes-M4-geo"       = "Gravity + covariates",
  "Bayes-M8-geo"       = "Composite gravity + covariates",
  "Bayes-M9-geo"       = "Multi-kernel ensemble + covariates",
  "Bayes-M10-geo"      = "Radiation-composite + covariates",
  "Bayes-M13-geo"      = "Cohort gravity + covariates",
  "Bayes-M14-geo"      = "Cohort radiation + covariates",
  "Bayes-M15-geo"      = "Combined static + covariates",
  "Bayes-M16-geo"      = "Cohort + directed OD + covariates",
  "Bayes-M17-geo"      = "All-kernel consensus + covariates",
  # road-distance (-dist) twins (ON by default): generic + cohort + consensus, with & without covariates.
  # NOTE model ids: M17's twin is "Bayes-M17-dist-med"/"-geo" (carries the -med/-geo suffix), unlike the
  # M13/M14 twins which are "Bayes-M13-dist"/"Bayes-M13-dist-geo". Keys MUST match the fitted labels exactly.
  "Bayes-M4-dist"       = "Gravity (road-km)",
  "Bayes-M4-dist-geo"   = "Gravity (road-km) + covariates",
  "Bayes-M8-dist"       = "Composite gravity (road-km)",
  "Bayes-M8-dist-geo"   = "Composite gravity (road-km) + covariates",
  "Bayes-M10-dist"      = "Radiation-composite (road-km)",
  "Bayes-M10-dist-geo"  = "Radiation-composite (road-km) + covariates",
  "Bayes-M13-dist"      = "Cohort gravity (road-km)",
  "Bayes-M13-dist-geo"  = "Cohort gravity (road-km) + covariates",
  "Bayes-M14-dist"      = "Cohort radiation (road-km)",
  "Bayes-M14-dist-geo"  = "Cohort radiation (road-km) + covariates",
  "Bayes-M17-dist-med"  = "All-kernel consensus (road-km)",
  "Bayes-M17-dist-geo"  = "All-kernel consensus (road-km) + covariates",
  "Bayes-M8-short"     = "Composite gravity, short gen.",
  "Bayes-M8-long"      = "Composite gravity, long gen.",
  "Bayes-M8-full"      = "Composite gravity + full covariates",
  "Bayes-M8-susp"      = "Composite gravity + suspected",
  "Bayes-M8-full-susp" = "Composite gravity + cov. & susp.",
  "Bayes-ens-mean"     = "Model ensemble (mean)",
  "Bayes-ens-median"   = "Model ensemble (median)")

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

rs  <- read_csv(file.path(OUT, "key_outputs", "bayes_risk_scores_all_zones.csv"),
                show_col_types = FALSE)
lfo <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))
ev  <- read_csv(file.path(OUT, "diagnostics", "invasion_evaluation.csv"),
                show_col_types = FALSE)

# The featured Bayesian model is chosen by the pipeline's calibration-aware CV composite — the
# summed within-horizon ranks of AUC-PR skill (discrimination) + mean rank-of-truth (operational) +
# log-score (probabilistic accuracy), POOLED across both horizons; lower = better — over the
# cross-validated SINGLE Bayesian models (ensembles excluded). This mirrors run_all.R's
# best_invasion_model(), so the figures always track the current featured model; NOT the loo-stacking
# weight, NEVER a hardcoded label.
.pick_best_bayes <- function(rs, lfo, ev, fallback) {
  lfo_methods   <- unique(lfo$method)
  bayes_singles <- unique(lfo_methods[grepl("^Bayes", lfo_methods) & !grepl("-ens-", lfo_methods)])
  need <- c("method", "horizon", "auc_pr_skill", "mean_rank_of_truth", "log_score")
  if (length(bayes_singles) && all(need %in% names(ev))) {
    # Discrimination-led with a spiky-model gate. Mirrors best_invasion_model() exactly: among
    # methods covering the most horizons, drop those whose worst-horizon mean rank-of-truth exceeds
    # 1.5x the field median, then take the highest total AUC-PR skill (ties: lower rank-of-truth,
    # then lower log-score); gate falls back to the ungated set if it would leave nothing.
    agg <- ev %>% filter(method %in% bayes_singles) %>%
      group_by(method) %>%
      summarise(n_h = dplyr::n_distinct(horizon), aps = sum(dplyr::coalesce(auc_pr_skill, 0)),
                mrt_max = max(dplyr::coalesce(mean_rank_of_truth, Inf)),
                mrt_sum = sum(dplyr::coalesce(mean_rank_of_truth, Inf)),
                ls_sum  = sum(dplyr::coalesce(log_score, Inf)), .groups = "drop") %>%
      filter(n_h == max(n_h))
    if (nrow(agg)) {
      finite_mrt <- agg$mrt_max[is.finite(agg$mrt_max)]
      gate <- if (length(finite_mrt)) 1.5 * stats::median(finite_mrt) else Inf
      keep <- agg %>% filter(mrt_max <= gate); if (!nrow(keep)) keep <- agg
      return(keep %>% arrange(desc(aps), mrt_sum, ls_sum) %>% slice(1) %>% pull(method))
    }
  }
  # fallbacks: the risk CSV's own method (if it is a CV single model), then any Bayes single, then literal
  if ("method" %in% names(rs)) {
    mb <- intersect(unique(stats::na.omit(rs$method)), bayes_singles)
    if (length(mb) >= 1) return(mb[1])
  }
  if (length(bayes_singles)) return(bayes_singles[1])
  if (fallback %in% lfo_methods) return(fallback)
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
# Naive per-zone score = population-weighted epicentre outflow reaching each zone (case-blind).
NAIVE_SCORES <- tryCatch({
  W   <- readRDS(file.path(OUT, "mobility", sprintf("mobility_%s.rds", .cfg("MOBILITY_PRIMARY", "M8"))))
  wp  <- read_csv(file.path(.cfg("WORLDPOP_DIR",
                 file.path(HERE, "..", "data", "worldpop", "processed")),
                 "worldpop__pop_count__static.csv"), show_col_types = FALSE)
  pv  <- setNames(wp$pop_count, wp$nom); zall <- names(pv)
  ez  <- intersect(as.character(EPI_ZONES), rownames(W))
  if (!length(ez)) stop("no epicentre zones in mobility matrix")
  p <- as.numeric(pv[match(ez, names(pv))]); good <- is.finite(p) & p > 0
  src <- if (!any(good)) rep(1, length(ez)) else { p[!good] <- stats::median(p[good]); p }
  score <- as.numeric(crossprod(src, W[ez, , drop = FALSE])); names(score) <- colnames(W)
  out <- setNames(rep(0, length(zall)), zall)
  cm  <- intersect(zall, names(score)); out[cm] <- score[cm]
  out[intersect(zall, ez)] <- 0                       # epicentre = origin, never an at-risk target
  out
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
ll <- suppressWarnings(read_csv(LINELIST, show_col_types = FALSE, guess_max = 10000)) %>%
  mutate(onset  = suppressWarnings(as.Date(date_of_symptom_onset)),
         sample = suppressWarnings(as.Date(date_of_sample_collection)),
         confirmed = final_mve_case_classification == "confirmed_case")

# --- reconcile the figure's line list to the INSP sitrep (SAME top-up the model applies) ---
# load_linelist() reconciles the DHIS2 line list to the sitrep's cumulative confirmed counts in
# memory (.build_sitrep_confirmed_appends): for zones the sitrep confirms but the line list has not
# caught up on — e.g. Oicha — it appends the shortfall as confirmed rows. This figure reads the
# processed CSV DIRECTLY (it never calls load_linelist), so without this it would silently omit
# those zones and disagree with the invasion model. We reuse the pipeline's OWN helper (single
# source of truth — no duplicated logic) and backfill `province` from the shapefile, since the
# helper leaves it blank (the model derives province from the spine, not the line list). The
# appended rows carry a sample date and blank onset, so the fixed-delay date_index logic below
# dates them exactly as it already dates any other appended confirmation.
try({
  if (!exists(".build_sitrep_confirmed_appends"))
    source(file.path(HERE, "spatiotemporal", "01_data_prep.R"))
  .al <- if (exists("ALIASES_PATH") && file.exists(ALIASES_PATH))
           suppressWarnings(readr::read_csv(ALIASES_PATH, show_col_types = FALSE,
                            col_types = readr::cols(.default = "c"))) else NULL
  .app <- .build_sitrep_confirmed_appends(ll, .al)
  if (nrow(.app) > 0) {
    prov_lk <- shp %>% st_drop_geometry() %>%
      transmute(.k = .key, province = .prov) %>% distinct(.k, .keep_all = TRUE)
    .app <- .app %>%
      mutate(.k = tolower(trimws(health_zone))) %>%
      left_join(prov_lk, by = ".k") %>% dplyr::select(-.k) %>%
      mutate(onset     = suppressWarnings(as.Date(date_of_symptom_onset)),
             sample    = suppressWarnings(as.Date(date_of_sample_collection)),
             confirmed = final_mve_case_classification == "confirmed_case")
    ll <- dplyr::bind_rows(ll, .app)
    message(sprintf("[fig1] sitrep top-up: +%d confirmed row(s) across %d zone(s): %s",
                    nrow(.app), dplyr::n_distinct(.app$health_zone),
                    paste(sort(unique(.app$health_zone)), collapse = ", ")))
  }
}, silent = TRUE)

# date_index: usable onset, else sample - mean delay, clamped to outbreak-week floor
wk_floor <- floor_date(OUTBREAK_START, "week", week_start = 1)
ll <- ll %>% mutate(
  onset_usable = !is.na(onset) & onset >= wk_floor & (is.na(sample) | onset <= sample) &
                 (is.na(sample) | onset >= sample - 90),
  date_index = as.Date(ifelse(onset_usable, onset,
                              pmax(sample - IMP_DELAY, wk_floor)), origin = "1970-01-01"))
conf <- ll %>% filter(confirmed, !is.na(date_index), date_index >= wk_floor,
                      !is.na(health_zone))
# Share of confirmed cases whose onset was imputed (for the caption; computed, not
# hardcoded). NB this standalone figure uses a FIXED onset<-sample shift (IMP_DELAY)
# for legibility; the pipeline itself draws the onset->sample delay per record.
PCT_IMPUTED <- round(100 * mean(!conf$onset_usable, na.rm = TRUE), 1)

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
         paste0("Onset imputed (fixed ", IMP_DELAY, "-day shift) for the ", PCT_IMPUTED,
                "% of cases lacking a usable onset; pipeline draws the delay per record.")) +
    theme_pub(16.5) + theme(plot.margin = margin(0,8,6,6),
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
  cap_share <- if (is.finite(ITURI_PCT)) sprintf(" ~%.0f%% within %s.", ITURI_PCT, ITURI_PROV) else ""
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
  # Optional intuitive per-model relabelling (panels A-C carry model names). When
  # model_labels is supplied, the model codes on the y-axes / legend are swapped for
  # reader-friendly names and a decoding caption is added; NULL = original codes.
  relab <- function(v) if (is.null(model_labels)) v else
    unname(ifelse(as.character(v) %in% names(model_labels), model_labels[as.character(v)], as.character(v)))
  y_relabel <- if (is.null(model_labels)) NULL else scale_y_discrete(labels = function(v) relab(v))
  yt_size   <- if (is.null(model_labels)) 11 else 8.2

  # --- 2A discrimination: AUC-PR skill for the Bayesian grid, h1 & h2 ---
  d <- ev %>% filter(str_detect(method, "^Bayes"), horizon %in% c(1,2)) %>%
    mutate(lo = pmax(auc_pr_lo,0)/base_rate, hi = auc_pr_hi/base_rate)
  ord <- d %>% filter(horizon == 1) %>% arrange(auc_pr_skill) %>% pull(method)
  if (!length(ord)) ord <- d %>% arrange(auc_pr_skill) %>% pull(method) %>% unique()
  d <- d %>% mutate(method = factor(method, levels = ord),
                    hz = factor(horizon))
  p2a <- ggplot(d, aes(auc_pr_skill, method, colour = hz)) +
    geom_vline(xintercept = 1, linetype = "22", colour = FAINT) +
    geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .55),
                   linewidth = .5, alpha = .55) +
    geom_point(position = position_dodge(width = .55), size = 1.9) +
    scale_colour_manual(values = HZ_COL, name = "Horizon",
                        labels = c("1 week","2 weeks")) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = "Discrimination of the Bayesian invasion models",
         subtitle = "AUC-PR skill = average precision / base rate  (1 = no skill; higher = better)",
         x = "AUC-PR skill (x base rate)", y = NULL,
         caption = "Points = pooled leave-future-out estimate; bars = 90% zone-cluster bootstrap CI.") +
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
         x = "Zones actively monitored each week (K)",
         y = "Share of true invasions caught") +
    theme_pub(13.5)

  # (Panels 2C prioritisation and 2D predicted-vs-observed have been dropped from
  #  Figure 2 at the author's request; the figure now carries A, B, E, F only.)

  # --- 2E ranking accuracy: mean rank of the invaded zones, Bayesian grid, h1 & h2 ---
  # Lower = the truly-invaded zones sat nearer the TOP of the model's watch-list (a direct
  # operational readout, complementary to AUC-PR). Sourced from invasion_evaluation.csv
  # (mean_rank_of_truth); the two horizons are joined per model so the h1->h2 shift is visible.
  mr <- ev %>% filter(str_detect(method, "^Bayes"), horizon %in% c(1,2),
                      is.finite(mean_rank_of_truth))
  ord_mr <- mr %>% filter(horizon == 1) %>% arrange(desc(mean_rank_of_truth)) %>% pull(method)
  if (!length(ord_mr)) ord_mr <- mr %>% arrange(desc(mean_rank_of_truth)) %>% pull(method) %>% unique()
  mr <- mr %>% mutate(method = factor(method, levels = ord_mr), hz = factor(horizon))
  mr_seg <- mr %>% dplyr::select(method, horizon, mean_rank_of_truth) %>%
    pivot_wider(names_from = horizon, values_from = mean_rank_of_truth, names_prefix = "mr")
  p2e <- ggplot(mr, aes(mean_rank_of_truth, method)) +
    { if (all(c("mr1","mr2") %in% names(mr_seg)))
        geom_segment(data = mr_seg, aes(x = mr1, xend = mr2, y = method, yend = method),
                     colour = FAINT, linewidth = 0.5, na.rm = TRUE, inherit.aes = FALSE) } +
    geom_point(aes(colour = hz), size = 2) +
    scale_colour_manual(values = HZ_COL, name = "Horizon", labels = c("1 week","2 weeks")) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.08))) +
    labs(title = "Ranking accuracy for the invaded zones",
         subtitle = "Mean rank of the truly-invaded zones among at-risk zones (lower = nearer the top)",
         x = "Mean rank of invaded zones (lower is better)", y = NULL,
         caption = "Rank of every invaded zone (tie-averaged), meaned per fold then over folds.") +
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
  N_TOP <- 25L; RANK_FLOOR <- 30L
  rk <- lfo %>% filter(method == BEST_BAYES, horizon == 1, is.finite(p_invasion),
                       !(was_active_before %in% TRUE)) %>%
    mutate(cutoff = as.Date(cutoff)) %>%
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
      scale_x_date(date_labels = "%d %b", breaks = sort(unique(traj$cutoff)),
                   expand = expansion(mult = c(0.03, 0.30))) +
      scale_colour_manual(values = PROV_COL, name = "Province", breaks = reg_levels,
                          na.value = "#7A7A7A") +
      labs(title = "Invasion-ranking evolution across forecast rounds",
           subtitle = sprintf("Rank of each zone's 1-week invasion risk (1 = highest); the current top %d, tracked over folds", N_TOP),
           x = "Forecast round (fold cutoff)", y = "Invasion-risk rank",
           caption = "Each line = a zone currently in the top watch-list; ranks worse than 30 shown at the 30+ baseline.") +
      theme_pub(13.5) +
      theme(legend.position = "right",
            panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3))
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
  p_topk1 <- build_topk_folds(horizon = topk_horizon, facet_nrow = topk_nrow,
                              save = FALSE, embed = TRUE)

  # Decoding caption for the reader-friendly labelled variant (added only when model_labels set).
  cap2 <- if (is.null(model_labels)) NULL else stringr::str_wrap(paste0(
    "Model = the spatial mobility kernel driving invasion risk: gravity / composite gravity / ",
    "radiation-composite, their Flowminder cohort and combined-static counterparts, and multi-kernel ",
    "or all-kernel consensus ensembles of these. \"+ covariates\" = geographic & social covariates ",
    "modulate the import-to-invasion rate (base models use a single constant rate); \"road-km\" = ",
    "road-distance rather than travel-time deterrence; \"+ suspected\" = suspected-case leading ",
    "indicators; \"short/long gen.\" = assumed generation-time profile; \"ensemble\" = stacked average ",
    "across models."),
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
  # rr01_nat is heavily right-skewed (~78% of at-risk zones < 0.01), so a linear 0-1 fill collapses
  # the national map to a single dark hue. A log10 fill shows HONEST magnitudes (risk truly
  # concentrated near the epicentre) while staying legible; the pseudo-log floor RR_FLOOR absorbs
  # the ~0.4% of zones with exactly-zero relative risk (log(0) is undefined).
  sc_rank <- scale_fill_viridis_c(option="viridis", trans="log10", limits=c(RR_FLOOR, 1),
              breaks=c(0.001,0.01,0.1,1), labels=c("0.001","0.01","0.1","1"),
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
build_topk_folds <- function(top_k = 12L, horizon = 1L, save = TRUE, embed = FALSE,
                             facet_nrow = 1L) {
  OUT_HIT  <- "#D55E00"   # invaded (first case occurred)  — Okabe-Ito vermilion
  OUT_MISS <- "grey75"    # not invaded this round
  OUT_COL  <- c("Invaded (first case)" = OUT_HIT, "Not invaded" = OUT_MISS)

  d <- lfo %>%
    filter(method == BEST_BAYES, horizon == !!horizon, is.finite(p_invasion),
           !(was_active_before %in% TRUE))
  if (!nrow(d)) { message("[topk_folds] no at-risk rows; skipped."); return(invisible(NULL)) }

  # top-K zones by predicted probability WITHIN each fold. Facet header = round + date only.
  d <- d %>% mutate(cutoff = as.Date(cutoff))
  fold_ord <- sort(unique(d$cutoff))
  lab_lk   <- setNames(sprintf("Round %d — %s", seq_along(fold_ord), format(fold_ord, "%d %b")),
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
    labs(title = ttl, subtitle = sub,
         x = "Predicted invasion probability, P(first case)", y = NULL) +
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
ok2b <- run(function() build_fig2(model_labels = MODEL_LABELS, file_suffix = "_labelled",
                                  fig_w = 12.6, fig_h = 17.4,
                                  topk_horizon = 2L, topk_nrow = 2L), "Figure 2 (labelled)")
ok3 <- run(build_fig3, "Figure 3")
okt <- run(function() { build_topk_folds(horizon = 1L); build_topk_folds(horizon = 2L) },
           "Top-K per-fold outcome bars")
message(sprintf("\n[done] Figure1:%s  Figure2:%s  Figure2-labelled:%s  Figure3:%s  TopKfolds:%s  ->  %s",
                ok1, ok2, ok2b, ok3, okt, FIG_DIR))
