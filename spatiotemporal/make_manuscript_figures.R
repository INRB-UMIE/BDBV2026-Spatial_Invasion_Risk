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
#                                        stacked cumulative cases (3 earliest zones + others)
#                                        on the left, health zones invaded on the right
#                                     C) arrival time (x) vs road travel time + May M+B+R
#                                        outflow (R^2 + Spearman rho). All 5 predictors -> SI.
#   Figure S1 arrival time (x) vs all 5 predictors (May & March mobility effective distance,
#                                        great-circle, road travel time, log10 population)
# All forecast-based panels are fixed to the 2-WEEK-AHEAD horizon (Figures 2 & 3).
# Figure 1 is descriptive (observed invasion), so it carries no forecast horizon.
#   Figure 2  Predictive performance  A) prioritisation (model vs naive vs random), 2 wk,
#                                        95% Wilson bands + AUC/Brier annotation
#                                     B) forecast-vs-outcome for 3 evenly-spread rounds, 2 wk;
#                                        top-12 at-risk zones coloured invaded-within-round /
#                                        invaded-later (by 31 Jul) / not-invaded (by 31 Jul)
#   Figure 3  Invasion + forecast     MAIN TEXT (2-wk horizon):
#                                     A) relative-invasion-risk map
#                                     B) top-20 forest (invasion probability + 90% CrI, by province)
#   Figure S4 Front & rank narrative  SI (split out of the former Figure 3 A/B):
#                                     A) front-approach: distance to nearest invaded zone per round,
#                                        evolving to the present (explains B; shared axis/colour/zones)
#                                     B) 2-wk invasion-risk rank evolution across rounds up to the
#                                        present live run (tracks main-text Fig 3B's top zones)
#   Figure S2 operational detail      A) prob x vulnerability bivariate  B) priority scatter
#
# Outputs -> outputs/key_outputs/manuscript_figures/{Figure1,2,3}.{pdf,png} + panels/
# Run:  Rscript make_manuscript_figures.R   (from spatiotemporal/)
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
OUT  <- file.path(HERE, "outputs")
DATA <- file.path(HERE, "..", "data")
FIG_DIR   <- file.path(OUT, "key_outputs", "manuscript_figures")
PANEL_DIR <- file.path(FIG_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# --- constants (mirror 00_config.R; kept local so the script is standalone) ---
OUTBREAK_START <- as.Date("2026-04-30")
ARRIVAL_REF    <- as.Date("2026-04-05")       # Figure 1C arrival-time origin
DELAY_RATE     <- 0.228                       # onset->sample exp. rate (/day)
IMP_DELAY      <- round(1 / DELAY_RATE)        # fixed onset<-sample shift for the standalone figure
SHP_PATH       <- file.path(DATA, "shapefiles", "DRC_Health_zones.shp")
LINELIST_DIR   <- file.path(DATA, "processed", "dhis2_linelist_processed")
LINELIST_JSON  <- file.path(LINELIST_DIR, "latest.json")
PROV_INT       <- c("Ituri", "Nord-Kivu", "Haut-Uele")
EPI_ZONES      <- c("Bunia", "Mongbwalu", "Rwampara")     # M+B+R epicentre for Figure 1C
FEATURED_FALLBACK <- "Bayes-M13-dist"
key <- function(x) tolower(trimws(x))

# -----------------------------------------------------------------------------
# 1. DESIGN SYSTEM  (colourblind-safe)
# -----------------------------------------------------------------------------
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
AFFECTED_FILL <- "grey78"; NA_FILL <- "grey93"
# distinct warm neutral for zones with NO confirmed cases (Figure 1A) — deliberately OFF the
# grey→red→blue bivariate ramp so "no cases" is never read as the low-cases/early-onset corner.
NO_CASE_FILL <- "#E7DCC1"
RR_FLOOR <- 1e-3
ACCENT   <- "#EE7733"                          # locator-box accent
PT_BLUE  <- "#4C78C8"; FIT_RED <- "#C0392B"    # Figure 1C points / linear fit
OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
PROV_COL <- c("Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73",
              "Tshopo"="#CC79A7","Bas-Uele"="#E69F00","Sud-Kivu"="#56B4E9","Other"="#7A7A7A")
HZ_COL   <- c("1"="#3B4CC0","2"="#B4413C")
base_family <- "sans"

theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    plot.title    = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    axis.title    = element_text(size = base - 0.4, colour = MUTED),
    axis.title.x  = element_text(margin = margin(t = 4)),
    axis.title.y  = element_text(margin = margin(r = 4), angle = 90),
    axis.text     = element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
    legend.position = "top", legend.justification = "left",
    legend.title  = element_text(size = base - 1.2, colour = MUTED),
    legend.text   = element_text(size = base - 1.4, colour = INK),
    legend.key.height = unit(9, "pt"), legend.key.width = unit(15, "pt"),
    strip.text    = element_text(size = base - 0.6, colour = INK, face = "bold", margin = margin(3,3,3,3)),
    plot.tag      = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin   = margin(6, 8, 6, 6))
}
theme_map <- function(base = 8.6) {
  theme_void(base_size = base, base_family = base_family) %+replace% theme(
    plot.title    = element_blank(),
    legend.position = "right",
    legend.title  = element_text(size = base - 1.4, colour = MUTED),
    legend.text   = element_text(size = base - 1.8, colour = INK),
    legend.key.height = unit(20, "pt"), legend.key.width = unit(7, "pt"),
    plot.tag      = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin   = margin(2, 2, 2, 2))
}
save_dual <- function(p, name, w, h, dir = PANEL_DIR) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  message(sprintf("  saved %-32s  %.1f x %.1f in", name, w, h)); invisible(p)
}

# --- bivariate-choropleth helpers (shared by Figures 1 and 3) --------------------------------
qtile <- function(x) {
  out <- rep(NA_integer_, length(x)); ok <- is.finite(x)
  if (!any(ok)) return(out)
  r <- rank(x[ok], ties.method = "average") / sum(ok)
  out[ok] <- as.integer(cut(r, c(-Inf, .25, .5, .75, Inf), labels = 1:4)); out
}
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
          axis.title = element_text(size = base - 0.8, colour = MUTED),
          plot.title = element_text(size = base - 0.6, colour = INK, face = "bold", hjust = 0.5, margin = margin(b = 2)),
          plot.margin = margin(2,2,2,2))
  p
}

# -----------------------------------------------------------------------------
# 2. LOAD SAVED DATA
# -----------------------------------------------------------------------------
message("[load] shapefile, risk scores, LFO, evaluation, line list ...")
shp <- st_read(SHP_PATH, quiet = TRUE) %>% mutate(.key = key(Nom), .prov = as.character(PROVINCE))
rs  <- read_csv(file.path(OUT, "key_outputs", "bayes_risk_scores_all_zones.csv"), show_col_types = FALSE)
lfo <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))
ev  <- read_csv(file.path(OUT, "diagnostics", "invasion_evaluation.csv"), show_col_types = FALSE)

# Mirror of best_invasion_model() (16_invasion_eval.R): discrimination-led selection among single
# Bayesian models — highest total AUC-PR skill, AFTER a spiky-model gate that drops any model whose
# worst-horizon mean rank-of-truth exceeds MRT_GATE_MULT x the field median (excludes over-confident
# variants that nail a few top hits but rank the rest poorly). Ties break by lower total rank-of-truth
# then lower log-score. Only methods covering the most horizons are eligible.
MRT_GATE_MULT <- 1.5
pick_featured <- function(ev, fallback) {
  need <- c("method","horizon","auc_pr_skill","mean_rank_of_truth","log_score")
  singles <- unique(ev$method[grepl("^Bayes", ev$method) & !grepl("-ens-", ev$method)])
  if (!length(singles) || !all(need %in% names(ev))) return(fallback)
  agg <- ev %>% filter(method %in% singles) %>% group_by(method) %>%
    summarise(n_h = dplyr::n_distinct(horizon), aps = sum(dplyr::coalesce(auc_pr_skill, 0)),
              mrt_max = max(dplyr::coalesce(mean_rank_of_truth, Inf)),
              mrt_sum = sum(dplyr::coalesce(mean_rank_of_truth, Inf)),
              ls_sum  = sum(dplyr::coalesce(log_score, Inf)), .groups = "drop") %>%
    filter(n_h == max(n_h))
  if (!nrow(agg)) return(fallback)
  finite_mrt <- agg$mrt_max[is.finite(agg$mrt_max)]
  gate <- if (length(finite_mrt)) MRT_GATE_MULT * stats::median(finite_mrt) else Inf
  keep <- agg %>% filter(mrt_max <= gate); if (!nrow(keep)) keep <- agg
  keep %>% arrange(desc(aps), mrt_sum, ls_sum) %>% slice(1) %>% pull(method)
}
FEATURED <- pick_featured(ev, FEATURED_FALLBACK)
message(sprintf("[featured] %s", FEATURED))

join_map <- function(dat) {
  d <- dat %>% mutate(.key = key(health_zone), .prov = as.character(province))
  m <- shp %>% left_join(d, by = c(".key", ".prov"))
  miss <- m %>% st_drop_geometry() %>% summarise(f = mean(is.na(health_zone))) %>% pull(f)
  if (miss > 0.5) { d2 <- dat %>% mutate(.key = key(health_zone))
    m <- shp %>% left_join(d2 %>% dplyr::select(-any_of(".prov")), by = ".key") }
  m
}
prov_lookup <- rs %>% distinct(health_zone, province)

# "present" as-of date for the LIVE run (the final time point appended to the rank/front
# panels, downstream of the CV folds). Read from run_info.json; fall back to the last CV
# fold + 28d if unavailable so the panels still render.
PRESENT_DATE <- tryCatch({
  ri <- jsonlite::fromJSON(file.path(OUT, "key_outputs", "run_info.json"))
  d <- as.Date(ri$analysis_date); if (is.na(d)) stop("no analysis_date"); d
}, error = function(e) { message("[present] run_info analysis_date unavailable (", conditionMessage(e), ")"); as.Date(NA) })
message(sprintf("[present] live-run as-of date: %s", ifelse(is.na(PRESENT_DATE), "NA (fallback)", format(PRESENT_DATE))))

# Province -> colour map shared by ALL Figure-3 zone panels (A, B, D) so every province
# gets its own consistent hue: fixed hues from PROV_COL where defined, deterministic spare
# Okabe-Ito hues for any province not in PROV_COL. Level order = the three provinces of
# interest first, then the remainder alphabetically. Returns levels + named colour vector.
province_palette <- function(provinces) {
  provs <- unique(as.character(provinces[!is.na(provinces)]))
  levs  <- c(intersect(PROV_INT, provs), sort(setdiff(provs, PROV_INT)))
  spare <- setdiff(OKABE, unname(PROV_COL)); si <- 1L
  cols  <- setNames(character(length(levs)), levs)
  for (pv in levs) {
    if (pv %in% names(PROV_COL)) cols[pv] <- unname(PROV_COL[pv])
    else { cols[pv] <- spare[((si - 1L) %% max(length(spare), 1L)) + 1L]; si <- si + 1L }
  }
  list(levels = levs, colours = cols)
}

# -----------------------------------------------------------------------------
# 3. FIGURE 1 DATA  (onset-dated confirmed cases; pipeline-consistent)
# -----------------------------------------------------------------------------
resolve_linelist <- function() {
  if (file.exists(LINELIST_JSON)) {
    meta <- tryCatch(jsonlite::fromJSON(LINELIST_JSON), error = function(e) NULL)
    if (!is.null(meta$folder)) {
      p <- file.path(LINELIST_DIR, meta$folder, "dhis2_processed_linelist.csv")
      if (file.exists(p)) return(p)
    }
  }
  cand <- list.files(LINELIST_DIR, pattern = "^LINELIST_", full.names = TRUE)
  if (length(cand)) file.path(sort(cand, decreasing = TRUE)[1], "dhis2_processed_linelist.csv")
  else stop("cannot resolve processed line list")
}
LINELIST <- resolve_linelist(); message("[fig1] line list: ", LINELIST)

ll <- suppressWarnings(read_csv(LINELIST, show_col_types = FALSE, guess_max = 10000)) %>%
  mutate(onset  = suppressWarnings(as.Date(date_of_symptom_onset)),
         sample = suppressWarnings(as.Date(date_of_sample_collection)),
         confirmed = final_mve_case_classification == "confirmed_case")

# reconcile to the INSP sitrep exactly as load_linelist() does (single source of truth); skipped
# gracefully if 01_data_prep.R cannot be sourced.
try({
  if (!exists(".build_sitrep_confirmed_appends")) source(file.path(HERE, "01_data_prep.R"))
  al <- if (exists("ALIASES_PATH") && file.exists(ALIASES_PATH))
          suppressWarnings(readr::read_csv(ALIASES_PATH, show_col_types = FALSE, col_types = readr::cols(.default = "c"))) else NULL
  app <- .build_sitrep_confirmed_appends(ll, al)
  if (nrow(app) > 0) {
    prov_lk <- shp %>% st_drop_geometry() %>% transmute(.k = .key, province = .prov) %>% distinct(.k, .keep_all = TRUE)
    app <- app %>% mutate(.k = key(health_zone)) %>% left_join(prov_lk, by = ".k") %>% dplyr::select(-.k) %>%
      mutate(onset = suppressWarnings(as.Date(date_of_symptom_onset)), sample = suppressWarnings(as.Date(date_of_sample_collection)),
             confirmed = final_mve_case_classification == "confirmed_case")
    ll <- dplyr::bind_rows(ll, app)
    message(sprintf("[fig1] sitrep top-up: +%d row(s), %d zone(s)", nrow(app), dplyr::n_distinct(app$health_zone)))
  }
}, silent = TRUE)

wk_floor <- floor_date(OUTBREAK_START, "week", week_start = 1)
ll <- ll %>% mutate(
  onset_usable = !is.na(onset) & onset >= wk_floor & (is.na(sample) | onset <= sample) & (is.na(sample) | onset >= sample - 90),
  date_index = as.Date(ifelse(onset_usable, onset, pmax(sample - IMP_DELAY, wk_floor)), origin = "1970-01-01"))
conf <- ll %>% filter(confirmed, !is.na(date_index), date_index >= wk_floor, !is.na(health_zone))

zone_cum   <- conf %>% group_by(health_zone, province) %>% summarise(cases = n(), .groups = "drop")
zone_first <- conf %>% group_by(health_zone) %>% summarise(first = min(date_index), .groups = "drop") %>%
  mutate(arrival_wk = as.integer(floor(as.numeric(first - wk_floor) / 7)) + 1L)
daily      <- conf %>% count(date_index, name = "n") %>% arrange(date_index) %>% mutate(cum_cases = cumsum(n))
zone_curve <- zone_first %>% arrange(first) %>% mutate(n_zones = row_number())
message(sprintf("[fig1] %d confirmed cases, %d zones, onset %s to %s", nrow(conf), nrow(zone_cum), min(daily$date_index), max(daily$date_index)))

EPI_ORIGIN <- EPI_ZONES[1]

# -----------------------------------------------------------------------------
# FIGURE 1C DATA — arrival time (days from ARRIVAL_REF) vs 5 predictors
#   May M+B+R outflow effective distance   = -log(share) of the latest May short-trip outflow
#   March Flowminder effective distance    = -log(share) of the March Flowminder O-D from M+B+R
#   great-circle distance (km)             = geodesic from the M+B+R centroid
#   road travel time from M+B+R (h)        = OSRM travel time, min over M+B+R, in hours
#   log10 population                       = log10 WorldPop count
# -----------------------------------------------------------------------------
compute_fig1c <- function() {
  arr <- zone_first %>% transmute(health_zone, k = key(health_zone), arr_days = as.numeric(first - ARRIVAL_REF))

  # zone coordinates for great-circle: prefer the effective-distance analysis's own node coords
  # (main-town points), else fall back to shapefile centroids.
  coords <- tryCatch({
    hp <- file.path(DATA, "effective_distance", "raw", "flowminder_effective_distance.html")
    stopifnot(file.exists(hp))
    txt <- paste(readLines(hp, warn = FALSE), collapse = " ")
    m <- regmatches(txt, gregexpr('"name":"[^"]+","lon":[-0-9.]+,"lat":[-0-9.]+', txt))[[1]]
    stopifnot(length(m) > 50)
    data.frame(k = key(sub('"name":"([^"]+)".*', "\\1", m)),
               lon = as.numeric(sub('.*"lon":([-0-9.]+),.*', "\\1", m)),
               lat = as.numeric(sub('.*"lat":([-0-9.]+)$', "\\1", m)))
  }, error = function(e) {
    ct <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(shp))))
    data.frame(k = shp$.key, lon = ct[, 1], lat = ct[, 2])
  })
  coords <- coords[!duplicated(coords$k), ]
  ec <- coords[coords$k %in% key(EPI_ZONES), ]; o <- c(mean(ec$lon), mean(ec$lat))
  hav <- function(lo, la) { Rk <- 6371; dl <- (lo - o[1]) * pi/180; dp <- (la - o[2]) * pi/180
    a <- sin(dp/2)^2 + cos(la * pi/180) * cos(o[2] * pi/180) * sin(dl/2)^2; 2 * Rk * asin(pmin(1, sqrt(a))) }
  gc <- data.frame(k = coords$k, gc_km = hav(coords$lon, coords$lat))

  # OSRM road travel time (min -> h), min over the M+B+R epicentre
  tt <- as.matrix(read.csv(file.path(DATA, "osrm", "processed", "osrm__travel_time__static.matrix.csv"), row.names = 1, check.names = FALSE))
  rownames(tt) <- key(rownames(tt)); colnames(tt) <- key(colnames(tt))
  ei <- intersect(key(EPI_ZONES), rownames(tt))
  th <- suppressWarnings(apply(tt[ei, , drop = FALSE], 2, function(z) { z <- z[is.finite(z)]; if (length(z)) min(z) else NA })) / 60
  road <- data.frame(k = names(th), road_h = as.numeric(th))

  # March Flowminder outflow -log(share) from M+B+R
  fm <- read.csv(file.path(DATA, "flowminder", "processed", "flowminder__outflow__static.matrix.csv"), row.names = 1, check.names = FALSE)
  rownames(fm) <- key(rownames(fm)); colnames(fm) <- key(colnames(fm))
  fi <- intersect(key(EPI_ZONES), rownames(fm))
  mo <- colSums(fm[fi, , drop = FALSE], na.rm = TRUE); mo <- mo / sum(mo, na.rm = TRUE)
  march <- data.frame(k = names(mo), march_eff = -log(ifelse(mo > 0, mo, NA)))

  # May M+B+R short-trip outflow -log(share) — latest May snapshot
  stf <- sort(list.files(file.path(DATA, "flowminder_short_trips", "processed"),
                         pattern = "^flowminder_short_trips__outflow_20260[0-9]{3}__static\\.csv$", full.names = TRUE))
  st <- read_csv(stf[length(stf)], show_col_types = FALSE); names(st)[1:2] <- c("nom", "out")
  st$k <- key(st$nom); st$sh <- st$out / sum(st$out, na.rm = TRUE)
  may <- data.frame(k = st$k, may_eff = -log(ifelse(st$sh > 0, st$sh, NA)))

  wp <- read_csv(file.path(DATA, "worldpop", "processed", "worldpop__pop_count__static.csv"), show_col_types = FALSE)
  pop <- data.frame(k = key(wp$nom), log10pop = log10(wp$pop_count))

  arr %>% left_join(may, by = "k") %>% left_join(march, by = "k") %>% left_join(gc, by = "k") %>%
    left_join(road, by = "k") %>% left_join(pop, by = "k")
}
fig1c_dat <- tryCatch(compute_fig1c(), error = function(e) { message("[fig1c] data build failed: ", conditionMessage(e)); NULL })

# =============================================================================
# FIGURE 1
# =============================================================================
# Panel-C predictor labels (shared by the main figure and the SI variant). These
# now sit on the y-axis strips because Panel C's axes are switched: arrival time
# on x, predictor on y.
METRIC_LAB <- c(
  may_eff   = "Mobility effective distance from M+B+R\n(May outflow; higher = smaller travel share)",
  march_eff = "Mobility effective distance from M+B+R\n(March Flowminder; higher = smaller travel share)",
  gc_km     = "great-circle distance from M+B+R\n(km)",
  road_h    = "road travel time from M+B+R\n(h)",
  log10pop  = "population (log10 people)")

# arrival-time (x) vs predictor (y), one facet per predictor, R^2 + Spearman rho.
# `metrics` selects which predictors to show (2 for the main figure, all 5 for SI).
# The predictor name is placed as the per-facet y-axis label (strip on the left), and
# the plot header is dropped, so each panel reads as a standard labelled scatter.
build_fig1c <- function(metrics, nrow = 1, base = 12.5, pt_size = 1.9, lab_size = 3.7) {
  labs_sel <- METRIC_LAB[metrics]
  long <- fig1c_dat %>%
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "y") %>%
    filter(is.finite(y), is.finite(arr_days)) %>%
    mutate(metric = factor(metric, levels = metrics, labels = labs_sel))
  stat_df <- long %>% group_by(metric) %>%
    summarise(r2  = summary(lm(y ~ arr_days))$r.squared,
              rho = suppressWarnings(cor(arr_days, y, method = "spearman")),
              xx  = min(arr_days), yy = max(y), .groups = "drop") %>%
    mutate(lab = sprintf("R^2 == %.2f * ',' ~~ rho == %.2f", r2, rho))
  ggplot(long, aes(arr_days, y)) +
    geom_point(colour = PT_BLUE, size = pt_size, alpha = 0.9) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = FIT_RED, linewidth = 0.9) +
    geom_text(data = stat_df, aes(xx, yy, label = lab), parse = TRUE, hjust = 0, vjust = 1, size = lab_size, colour = INK) +
    facet_wrap(~metric, scales = "free_y", nrow = nrow, strip.position = "left") +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.14))) +
    labs(x = sprintf("Arrival time (days from %s)", format(ARRIVAL_REF, "%Y-%m-%d")), y = NULL) +
    theme_pub(base) +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text = element_text(size = base - 0.4, colour = INK, face = "bold", lineheight = 0.95),
          panel.spacing.x = unit(12, "pt"), panel.spacing.y = unit(12, "pt"))
}

build_fig1 <- function() {
  # --- 1A: BIVARIATE choropleth — cumulative case burden (x) x date of first onset (y) ---
  af <- zone_cum %>% left_join(zone_first %>% dplyr::select(health_zone, first, arrival_wk), by = "health_zone") %>%
    mutate(bx = qtile(cases), by = qtile(arrival_wk),
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

  # correlation between arrival time and cumulative burden across invaded zones
  cor_df <- af %>% filter(!is.na(cases), !is.na(arrival_wk))
  rho_ac <- suppressWarnings(cor(cor_df$arrival_wk, cor_df$cases, method = "spearman"))
  cor_lab <- sprintf("rho[arrival~vs~burden] == %.2f ~~ '(n =' ~ %d ~ 'zones)'", rho_ac, nrow(cor_df))

  # province outlines (dissolve health zones by province)
  prov_bounds <- shp %>% dplyr::filter(!is.na(.prov)) %>% dplyr::group_by(.prov) %>%
    dplyr::summarise(.groups = "drop") %>% sf::st_make_valid()

  # Kisangani (Tshopo) reference city — an at-risk downstream hub, not yet invaded
  ki  <- shp[shp$.key == "makiso kisangani", ]
  kic <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(ki))))
  ki_df <- data.frame(x = kic[1, 1], y = kic[1, 2], lab = "Kisangani (Tshopo)")

  # extent: invaded cluster + epicentre + Kisangani, so the at-risk hub sits in-frame
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
  biv_leg <- bivar_legend(BIV_ARRIVAL, "Cumulative confirmed cases ->", "Later first onset ->",
                          base = 11)
  p1a <- p1a +
    patchwork::inset_element(locator, left = 0.00, bottom = 0.00, right = 0.26, top = 0.26, align_to = "panel", clip = FALSE) +
    patchwork::inset_element(biv_leg, left = 0.015, bottom = 0.60, right = 0.26, top = 0.88, align_to = "panel", clip = FALSE)

  # --- 1B: epidemic trajectory in ONE panel with two axes ---
  #   left axis  : cumulative confirmed cases, stacked by the M+B+R epicentre zones
  #                (Bunia, Mongbwalu, Rwampara) + all others (share shifts to "other" as it spreads)
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
    scale_x_date(limits = xlim, date_breaks = "2 weeks", date_labels = "%d %b",
                 expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Symptom-onset date") +
    theme_pub(13) +
    theme(legend.position = "top", legend.box = "vertical", legend.spacing.y = unit(1, "pt"),
          axis.title.y.right = element_text(colour = ZONE_LINE, angle = 90, margin = margin(l = 5)),
          axis.text.y.right  = element_text(colour = ZONE_LINE),
          axis.ticks.y.right = element_line(colour = ZONE_LINE, linewidth = 0.4),
          axis.line.y.right  = element_line(colour = ZONE_LINE, linewidth = 0.5),
          plot.margin = margin(6, 10, 6, 6))

  # --- 1C: arrival time (x) vs road travel time + May M+B+R outflow (main text) ---
  p1c <- build_fig1c(c("road_h", "may_eff"), nrow = 1, base = 12.5)

  save_dual(p1a, "F1A_bivariate_map", 7.0, 5.8)
  save_dual(p1b, "F1B_trajectory",    6.6, 5.0)
  save_dual(p1c, "F1C_arrival_vs_predictors", 7.0, 3.6)

  fig1 <- ((wrap_elements(full = p1a) | wrap_elements(full = p1b)) / wrap_elements(full = p1c)) +
    plot_layout(heights = c(1.5, 1)) + plot_annotation(tag_levels = "A")
  save_dual(fig1, "Figure1", 15.0, 12.0, dir = FIG_DIR)
  invisible(fig1)
}

# --- SI: arrival time vs ALL five predictors (same switched-axis layout as Panel C) ---
build_fig1_si <- function() {
  p <- build_fig1c(c("may_eff", "march_eff", "gc_km", "road_h", "log10pop"), nrow = 1, base = 16, lab_size = 4.7)
  save_dual(p, "FigureS1_arrival_vs_predictors", 16.5, 6.0, dir = FIG_DIR)
  invisible(p)
}

# --- SI: 5-day rolling share of confirmed cases in the M+B+R epicentre --------------------
# For each onset day, share = (cases in Bunia/Mongbwalu/Rwampara) / (all confirmed cases),
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
  p <- ggplot(d, aes(date_index)) +
    geom_point(aes(y = prop_daily), colour = FAINT, size = 1.1, alpha = 0.75, na.rm = TRUE) +
    geom_line(aes(y = prop5), colour = PT_BLUE, linewidth = 1.2, na.rm = TRUE) +
    annotate("text", x = min(d$date_index), y = 0.02, label = cap, hjust = 0, vjust = 0, size = 4.3, colour = MUTED) +
    scale_y_continuous(labels = percent_format(1), limits = c(0, 1), expand = expansion(mult = c(0.02, 0.04))) +
    scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b", expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Symptom-onset date",
         y = "Share of confirmed cases in the\nBunia + Mongbwalu + Rwampara epicentre") +
    theme_pub(base) + theme(axis.title.y = element_text(lineheight = 0.95))
  save_dual(p, "FigureS3_epicentre_case_share", 9.0, 5.6, dir = FIG_DIR)
  invisible(p)
}

# =============================================================================
# FIGURE 2  (key_outputs evaluation figure WITHOUT the across-model-grid panels)
# =============================================================================
reorder_within <- function(x, by, within)
  factor(paste(x, within, sep = "___"), levels = unique(paste(x, within, sep = "___"))[order(within, by)])
tidytext_scale_y <- function() scale_y_discrete(labels = function(z) sub("___.*$", "", z))

NAIVE_LBL <- "Naive-epicentre-inflow"
# Additional structural baselines drawn from the LFO grid (own predictions, same at-risk set
# and ground truth as the featured model), shown alongside the naive epicentre-inflow reference.
BASELINE_METHODS <- c("Gravity-B4", "Adjacency-B7")
BASELINE_LBL <- c("Gravity-B4" = "Gravity baseline", "Adjacency-B7" = "Adjacency baseline")
BASELINE_COL <- c("Gravity-B4" = "#009E73", "Adjacency-B7" = "#CC79A7")
NAIVE_SCORES <- tryCatch({
  W  <- readRDS(file.path(OUT, "mobility", "mobility_M8.rds"))
  wp <- read_csv(file.path(DATA, "worldpop", "processed", "worldpop__pop_count__static.csv"), show_col_types = FALSE)
  pv <- setNames(wp$pop_count, wp$nom); zall <- names(pv)
  ez <- intersect(as.character(EPI_ZONES), rownames(W)); stopifnot(length(ez) > 0)
  p <- as.numeric(pv[match(ez, names(pv))]); good <- is.finite(p) & p > 0
  src <- if (!any(good)) rep(1, length(ez)) else { p[!good] <- stats::median(p[good]); p }
  score <- as.numeric(crossprod(src, W[ez, , drop = FALSE])); names(score) <- colnames(W)
  out <- setNames(rep(0, length(zall)), zall); cm <- intersect(zall, names(score)); out[cm] <- score[cm]
  out[intersect(zall, ez)] <- 0; out
}, error = function(e) { message("[fig2] naive baseline unavailable (", conditionMessage(e), ")"); NULL })

# Wilson score interval for a binomial proportion x/n (clamped to [0,1]); used for
# uncertainty on the recall ("share of true invasions caught") curves.
.wilson <- function(x, n, z = 1.96) {
  if (!is.finite(n) || n <= 0) return(c(NA_real_, NA_real_))
  p <- x / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hlf <- z / d * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(max(0, ctr - hlf), min(1, ctr + hlf))
}
.det_curve <- function(d, ks = 1:25) {
  if (!nrow(d)) return(NULL)
  pf <- d %>% group_by(fold_id) %>% mutate(rk = rank(-p_invasion, ties.method = "max")) %>% ungroup()
  tot <- sum(pf$is_new_invasion, na.rm = TRUE); if (tot == 0) return(NULL)
  natr <- pf %>% count(fold_id) %>% pull(n) %>% mean()
  do.call(rbind, lapply(ks, function(k) {
    caught <- sum(pf$is_new_invasion[pf$rk <= k], na.rm = TRUE)
    ci <- .wilson(caught, tot)
    data.frame(k = k, recall = caught / tot, recall_lo = ci[1], recall_hi = ci[2],
               recall_random = pmin(k / natr, 1))
  }))
}
prioritisation_curves <- function(h = 2L, ks = 1:25) {
  feat <- lfo %>% filter(method == FEATURED, horizon == h, is.finite(p_invasion), !(was_active_before %in% TRUE))
  cf <- .det_curve(feat, ks); if (is.null(cf)) return(NULL); cf$method <- FEATURED
  curves <- cf
  # naive epicentre-inflow: a static mobility score applied to the SAME at-risk rows
  if (!is.null(NAIVE_SCORES)) {
    refn <- feat %>% mutate(p_invasion = as.numeric(NAIVE_SCORES[match(health_zone, names(NAIVE_SCORES))]))
    cn <- .det_curve(refn, ks); if (!is.null(cn)) { cn$method <- NAIVE_LBL; curves <- rbind(curves, cn) }
  }
  # structural baselines (gravity, adjacency): each ranks the same at-risk zones by its OWN prediction
  for (m in BASELINE_METHODS) {
    bm <- lfo %>% filter(method == m, horizon == h, is.finite(p_invasion), !(was_active_before %in% TRUE))
    cb <- .det_curve(bm, ks); if (!is.null(cb)) { cb$method <- m; curves <- rbind(curves, cb) }
  }
  list(curves = curves, rnd = cf[, c("k", "recall_random")])
}

# raw Brier score for predicted invasion probabilities (featured model, at-risk zones only)
brier_score <- function(h) {
  d <- lfo %>% filter(method == FEATURED, horizon == h, is.finite(p_invasion),
                      is.finite(is_new_invasion), !(was_active_before %in% TRUE))
  if (!nrow(d)) return(NA_real_)
  mean((d$p_invasion - d$is_new_invasion)^2)
}

# -----------------------------------------------------------------------------
# Shared rank trajectory (2-week horizon), used by BOTH Figure-3 narrative panels
# so they stay a matched pair. The tracked ("keep") zones are the top-N at-risk
# zones AS OF THE PRESENT LIVE RUN (rs) — the exact set shown in Figure 3D — and
# the x-axis runs across the 8 retrospective CV rounds (lfo) PLUS the present run,
# so the ranks are shown evolving all the way up to the latest forecast.
#   • CV rounds: rank within each fold by descending p_invasion  (lfo, horizon 2)
#   • present  : rank by descending p_case_invasion              (rs,  horizon 2)
# Both quantities are the same model's P(first case) within the 2-week window, so
# the per-timepoint orderings are comparable.
FRONT_H  <- 2L      # forecast horizon for the operational panels (2 weeks ahead)
FRONT_TOP_N <- 20L  # tracked at-risk zones (matches Figure 3D)
build_rank_traj <- function(N_TOP = FRONT_TOP_N, h = FRONT_H) {
  rk_h <- lfo %>% filter(method == FEATURED, horizon == h, is.finite(p_invasion), !(was_active_before %in% TRUE)) %>%
    mutate(cutoff = as.Date(cutoff)) %>% group_by(cutoff) %>%
    mutate(rank = rank(-p_invasion, ties.method = "min")) %>% ungroup() %>%
    transmute(health_zone, cutoff, rank)
  cutoffs <- sort(unique(rk_h$cutoff))
  present_date <- if (!is.na(PRESENT_DATE) && PRESENT_DATE > max(cutoffs)) PRESENT_DATE else max(cutoffs) + 28L
  now <- rs %>% filter(horizon == h, !(was_active_before %in% TRUE), is.finite(p_case_invasion)) %>%
    distinct(health_zone, .keep_all = TRUE) %>%
    mutate(rank = rank(-p_case_invasion, ties.method = "min"))
  keep_z <- now %>% arrange(desc(p_case_invasion)) %>% head(N_TOP) %>% pull(health_zone) %>% unique()
  rk_now <- now %>% transmute(health_zone, cutoff = present_date, rank)
  traj <- bind_rows(rk_h, rk_now) %>% filter(health_zone %in% keep_z)
  list(keep_z = keep_z, traj = traj, cutoffs = cutoffs, present_date = present_date)
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
  x_breaks <- sort(unique(traj$cutoff)); x_labs <- format(x_breaks, "%d %b")
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

build_fig2 <- function() {
  # Larger base font for Figure 2 only (both panels): more legible at print size than the
  # 8.6 default shared by the other figures. Hardcoded annotation/axis sizes below are
  # scaled up to match.
  F2_BASE <- 11.5
  disc <- ev %>% filter(method == FEATURED, horizon == 2) %>% arrange(horizon)
  gv <- function(col, h) disc[[col]][disc$horizon == h]
  bs2 <- brier_score(2)
  disc_lab <- sprintf("2-week-ahead forecast\nAUC-PR skill  %.2fx\nAUC-ROC  %.2f\nBrier score  %.3f",
                      gv("auc_pr_skill",2), gv("auc_roc",2), bs2)

  # --- 2A prioritisation: featured model vs naive epicentre-inflow, gravity & adjacency baselines,
  #     and random (2 weeks). 95% Wilson band is drawn on the featured model only, to keep it legible. ---
  pc <- prioritisation_curves(2L)
  MODEL_COL <- c(setNames(c("#0072B2", "#D55E00"), c(FEATURED, NAIVE_LBL)), BASELINE_COL)
  MODEL_LBL <- c(setNames(c("Mobility-informed model", "Naive epicentre-inflow"), c(FEATURED, NAIVE_LBL)), BASELINE_LBL)
  curves_b <- pc$curves %>% mutate(method = factor(method, levels = intersect(names(MODEL_COL), unique(method))))
  feat_c   <- curves_b %>% filter(method == FEATURED)
  p2a <- ggplot(curves_b, aes(k, recall, colour = method)) +
    geom_line(data = pc$rnd, aes(k, recall_random), linetype = "22", colour = FAINT, linewidth = 0.7, inherit.aes = FALSE) +
    geom_ribbon(data = feat_c, aes(ymin = recall_lo, ymax = recall_hi, fill = method), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) + geom_point(data = feat_c, size = 1) +
    annotate("text", x = 15, y = 0.075, label = "random watch-list", colour = MUTED, size = 3.5, angle = 8) +
    annotate("text", x = 0.5, y = 0.99, hjust = 0, vjust = 1, label = disc_lab, size = 3.4, colour = INK, lineheight = 0.95) +
    scale_colour_manual(values = MODEL_COL, labels = MODEL_LBL, name = NULL, breaks = levels(curves_b$method)) +
    scale_fill_manual(values = MODEL_COL, guide = "none") +
    scale_y_continuous(labels = percent_format(1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Zones actively monitored per round (K)", y = "Share of true invasions caught") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_pub(F2_BASE) + theme(legend.position = "top")

  # --- 2B forecast-vs-outcome: top-12 predicted-risk zones for 3 evenly-spread rounds (2 weeks) ---
  # THREE-WAY outcome (mirrors F_topk15_ever_h1, adapted to the 2-week horizon): each top-ranked
  # at-risk zone is coloured by whether it was invaded WITHIN this round's 2-week forecast window
  # (is_new_invasion == 1 at h=2 ⇒ first case in (cutoff, cutoff+14d]), invaded LATER but still by
  # the 31 Jul analysis snapshot (affected_ever, i.e. recorded a first confirmed case by the present
  # run but outside this round's window), or NOT invaded by 31 Jul. Palette/labels match
  # F_topk15_ever_h1 (dark burnt-orange / light amber / grey).
  OUT_COL <- c("Invaded within round"        = "#B33005",
               "Invaded later (by 31 Jul)"   = "#F6BB6B",
               "Not invaded (by 31 Jul)"     = "grey80")
  LV <- names(OUT_COL)
  # zones invaded at SOME point by the 31 Jul snapshot (present-run risk scores; was_active_before
  # marks zones already affected as of the analysis date). horizon == 1 slice is only to de-duplicate.
  affected_ever <- rs %>% filter(horizon == 1, as.logical(was_active_before) %in% TRUE) %>%
    pull(health_zone) %>% unique()
  dtk <- lfo %>% filter(method == FEATURED, horizon == 2, is.finite(p_invasion), !(was_active_before %in% TRUE)) %>% mutate(cutoff = as.Date(cutoff))
  all_folds <- sort(unique(dtk$cutoff))
  sel_idx  <- unique(round(seq(1, length(all_folds), length.out = 3)))   # first / middle / last
  fold_ord <- all_folds[sel_idx]
  dtk <- dtk %>% filter(cutoff %in% fold_ord)
  # honest round numbering: index within the full CV fold sequence, not 1..3
  lab_lk <- setNames(sprintf("Round %d: %s", match(fold_ord, all_folds), format(fold_ord, "%d %b")), as.character(fold_ord))
  topk <- dtk %>% group_by(cutoff) %>% slice_max(p_invasion, n = 12, with_ties = FALSE) %>% ungroup() %>%
    mutate(outcome = factor(dplyr::case_when(
                       is_new_invasion == 1            ~ LV[1],   # invaded within this round's 2-week window
                       health_zone %in% affected_ever  ~ LV[2],   # invaded later, but by 31 Jul
                       TRUE                            ~ LV[3]),  # not invaded by 31 Jul
                     levels = LV),
           fold_lab = factor(lab_lk[as.character(cutoff)], levels = lab_lk[as.character(fold_ord)]),
           zone_w = reorder_within(health_zone, p_invasion, fold_lab))
  p2b <- ggplot(topk, aes(p_invasion, zone_w, fill = outcome)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    facet_wrap(~ fold_lab, scales = "free_y", nrow = 1) + tidytext_scale_y() +
    scale_fill_manual(values = OUT_COL, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA), expand = expansion(mult = c(0, 0.06)), breaks = scales::pretty_breaks(3)) +
    labs(x = "Predicted invasion probability, P(first case)", y = NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_pub(F2_BASE) + theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(9, "pt"), strip.clip = "off",
                        axis.text.y = element_text(size = 9.6, colour = INK), axis.text.x = element_text(size = 8.8), legend.position = "top")

  save_dual(p2a, "F2A_prioritisation",     5.2, 4.8)
  save_dual(p2b, "F2B_forecast_vs_outcome", 2.9 * length(fold_ord) + 1.0, 4.8)

  fig2 <- (p2a | p2b) + plot_layout(widths = c(0.85, 1.5)) + plot_annotation(tag_levels = "A")
  save_dual(fig2, "Figure2", 13.8, 5.2, dir = FIG_DIR)
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
    for (t in times) {
      inv_keys <- fo$k[fo$first <= t & fo$k != zk]
      cc <- coord_lk[coord_lk$k %in% inv_keys, ]; if (!nrow(cc)) next
      rows[[length(rows) + 1L]] <- data.frame(health_zone = z, cutoff = as.Date(t, origin = "1970-01-01"),
        dist = min(.hav_km(zc[1], zc[2], cc$lon, cc$lat), na.rm = TRUE))
    }
  }
  df <- dplyr::bind_rows(rows) %>% left_join(prov_lookup, by = "health_zone") %>%
    mutate(province = factor(province, levels = pal$levels))

  # key-invasion event markers: date of first confirmed onset for notable zones (dotted
  # verticals) — validation checkpoints as the front advances up to the present run.
  ev <- tibble::tibble(
    health_zone = c("Komanda", "Nia Nia", "Makiso Kisangani", "Isiro", "Rethy"),
    lab         = c("Komanda", "Nia-Nia", "Kisangani", "Isiro", "Rethy")) %>%
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
    scale_x_date(date_labels = "%d %b", breaks = times, limits = c(x_min - 2, x_end),
                 expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.03, 0.18))) +
    labs(x = "Forecast round  ·  dotted = key first-onset dates", y = "Distance to nearest\ninvaded zone (km)") +
    theme_pub(base) + theme(legend.position = "right", panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
                            panel.grid.major.x = element_blank(),
                            axis.text.x = element_text(angle = 45, hjust = 1, size = base - 4.2))
}

# =============================================================================
# FIGURE 3 phase-space — the A<->B relationship as a single portrait
# =============================================================================
# Merges the two Figure-3 narrative panels: y = invasion-risk rank, x = distance to the
# nearest already-invaded zone. Each convex line is one zone's trajectory across the 2-week
# CV rounds (as the front closes in, distance falls and model rank climbs -> up-and-left).
# Lines are coloured by eventual outcome; the arrowhead marks the latest round. Rounds where
# a zone was actually invaded within the 2-week horizon (is_new_invasion) are crossed, so the
# distance threshold below which a high model rank converts into a real invasion is visible.
build_phase_space <- function(base = 12, RANK_KEEP = 25L, RANK_CAP = 40L) {
  h <- FRONT_H
  rk <- lfo %>% filter(method == FEATURED, horizon == h, is.finite(p_invasion), !(was_active_before %in% TRUE)) %>%
    mutate(cutoff = as.Date(cutoff)) %>% group_by(cutoff) %>%
    mutate(rank = rank(-p_invasion, ties.method = "min")) %>% ungroup() %>%
    transmute(health_zone, cutoff, rank, is_new = is_new_invasion, k = key(health_zone))
  coord_lk <- .get_centroids()
  fo <- zone_first %>% transmute(k = key(health_zone), first)   # first-onset date per invaded zone
  # keep only zones the model ever flagged (best rank <= RANK_KEEP) so the portrait stays legible
  keep <- rk %>% group_by(health_zone) %>% summarise(best = min(rank), .groups = "drop") %>%
    filter(best <= RANK_KEEP) %>% pull(health_zone)
  rk <- rk %>% filter(health_zone %in% keep) %>%
    left_join(coord_lk %>% rename(zlon = lon, zlat = lat), by = "k")
  # distance to nearest zone invaded by each round's cutoff (excluding the zone itself)
  dist_of <- function(zk, zlon, zlat, t) {
    if (is.na(zlon)) return(NA_real_)
    ik <- fo$k[fo$first <= t & fo$k != zk]; cc <- coord_lk[coord_lk$k %in% ik, ]
    if (!nrow(cc)) return(NA_real_)
    min(.hav_km(zlon, zlat, cc$lon, cc$lat), na.rm = TRUE)
  }
  rk$dist <- mapply(dist_of, rk$k, rk$zlon, rk$zlat, rk$cutoff)
  df <- rk %>% filter(is.finite(dist)) %>% arrange(health_zone, cutoff) %>%
    mutate(outcome = factor(ifelse(k %in% fo$k, "Eventually invaded", "Not (yet) invaded"),
                            levels = c("Eventually invaded", "Not (yet) invaded")),
           rank_c = pmin(rank, RANK_CAP))
  # imminent-invasion observations (zone invaded within the 2-week horizon of that round)
  events <- df %>% filter(is_new %in% 1)
  thr <- if (nrow(events)) as.numeric(stats::quantile(events$dist, 0.9, na.rm = TRUE)) else NA_real_
  path_df <- df %>% group_by(health_zone) %>% filter(dplyr::n() >= 2) %>% ungroup()
  ends <- df %>% group_by(health_zone) %>% slice_max(cutoff, n = 1, with_ties = FALSE) %>% ungroup()
  OUT_COL <- c("Eventually invaded" = "#D55E00", "Not (yet) invaded" = "#0072B2")
  # focus the x-axis on the informative front-proximity range: a few zones the model flagged via
  # long-range mobility sit >1000 km from any invaded zone and would otherwise squash the plot.
  # coord_cartesian zooms (keeps whole trajectories) rather than dropping points mid-path.
  inv_dist <- df$dist[df$outcome == "Eventually invaded"]
  XMAX <- max(thr * 1.5, as.numeric(stats::quantile(inv_dist, 0.98, na.rm = TRUE)), 150, na.rm = TRUE)
  XMAX <- ceiling(XMAX / 50) * 50
  n_far <- df %>% group_by(health_zone) %>% summarise(mn = min(dist), .groups = "drop") %>% filter(mn > XMAX) %>% nrow()

  p <- ggplot(df, aes(dist, rank_c, group = health_zone, colour = outcome))
  if (is.finite(thr)) p <- p +
    annotate("rect", xmin = -Inf, xmax = thr, ymin = -Inf, ymax = Inf, fill = "#E9895E", alpha = 0.12) +
    geom_vline(xintercept = thr, linetype = "22", colour = MUTED, linewidth = 0.5)
  p <- p +
    geom_path(data = path_df, linewidth = 0.7, alpha = 0.5, lineend = "round",
              arrow = grid::arrow(length = unit(5, "pt"), type = "closed", ends = "last")) +
    geom_point(data = events, shape = 4, size = 2.6, stroke = 1.0, colour = "#8A2E12", show.legend = FALSE) +
    geom_point(data = ends, size = 1.5) +
    scale_y_reverse(breaks = c(1, 5, 10, 20, 30, 40), labels = c("1","5","10","20","30", paste0(RANK_CAP, "+")),
                    expand = expansion(mult = c(0.05, 0.05))) +
    scale_x_continuous(breaks = pretty(c(0, XMAX), n = 5), expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(xlim = c(0, XMAX)) +
    scale_colour_manual(values = OUT_COL, name = NULL) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.4, alpha = 1))) +
    labs(x = "Distance to nearest invaded zone (km)  -  front closes in as this falls",
         y = "Invasion-risk rank (1 = highest)") +
    theme_pub(base) + theme(legend.position = "top", panel.grid.minor = element_blank())
  if (is.finite(thr)) {
    lab <- sprintf("crosses = invaded within 2 wk of that round\n90%% of these had the front within %.0f km", thr)
    if (n_far > 0) lab <- paste0(lab, sprintf("\n(%d far zone%s >%.0f km not shown)", n_far, ifelse(n_far == 1, "", "s"), XMAX))
    p <- p + annotate("text", x = thr, y = 1.2, hjust = -0.06, vjust = 1, size = 3.1, colour = MUTED,
                      lineheight = 0.95, fontface = "italic", label = lab)
  }
  save_dual(p, "Figure3_phase_space", 8.8, 6.6, dir = FIG_DIR)
  invisible(p)
}

build_fig3 <- function() {
  BIP <- BIV_PROBVUL
  prep <- function(h) rs %>% filter(horizon == h) %>%
    mutate(active = was_active_before %in% TRUE, rr_log = ifelse(active, NA, pmax(rr01_nat, RR_FLOOR)),
           bx = qtile(rr01_nat), by = qtile(V), bikey = ifelse(is.na(bx) | is.na(by), NA, paste0(bx, "-", by)))
  aff_overlay <- function(mp) geom_sf(data = mp %>% filter(was_active_before %in% TRUE), fill = AFFECTED_FILL, colour = "white", linewidth = 0.08)

  # --- 3A relative-invasion-risk map, 2 weeks ahead ---
  mp2 <- join_map(prep(2))
  sc_rr <- scale_fill_viridis_c(option = "viridis", trans = "log10", limits = c(RR_FLOOR, 1),
             breaks = c(0.001, 0.01, 0.1, 1), labels = c("0.001","0.01","0.1","1"),
             na.value = NA_FILL, name = "Relative\ninvasion risk\n(log, 0-1)", direction = 1)
  p3a <- ggplot(mp2) + geom_sf(aes(fill = rr_log), colour = "white", linewidth = 0.08) + aff_overlay(mp2) + sc_rr +
    theme_map(12) + theme(legend.position = "right")

  # --- 3B top-20 zones by 2-week invasion probability + 90% CrI, coloured by province ---
  TOP_N <- 20L
  d <- rs %>% filter(horizon == 2, !(was_active_before %in% TRUE), is.finite(p_case_invasion)) %>%
    distinct(health_zone, province, .keep_all = TRUE) %>% arrange(desc(p_case_invasion)) %>% head(TOP_N)
  # shared province palette (also drives Fig 3 A & B, so every province matches across panels)
  pal <- province_palette(d$province); prov_levels <- pal$levels; prov_cols <- pal$colours
  # disambiguate any duplicate zone names across provinces (else the y factor has duplicate levels)
  d <- d %>% mutate(province = factor(province, levels = prov_levels),
                    disp = ifelse(duplicated(health_zone) | duplicated(health_zone, fromLast = TRUE),
                                  paste0(health_zone, " (", as.character(province), ")"), health_zone))
  d$disp <- make.unique(d$disp)                       # final guard: guarantee unique factor levels
  d$zone <- factor(d$disp, levels = rev(d$disp))
  lev <- levels(d$zone); ylab_cols <- unname(prov_cols[as.character(d$province)[match(lev, d$disp)]])
  p3b <- ggplot(d, aes(p_case_invasion, zone, colour = province)) +
    geom_linerange(aes(xmin = p_case_lo, xmax = p_case_hi), linewidth = 1.4, alpha = 0.5) + geom_point(size = 2.4) +
    scale_colour_manual(values = prov_cols, name = "Province", drop = FALSE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA), expand = expansion(mult = c(0.01, 0.06))) +
    labs(x = "Probability of invasion, 2-week forecast", y = NULL) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2.2, size = 3))) +
    theme_pub(12) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = ylab_cols, size = 9.5, face = "bold"), legend.position = "top")

  # --- 3C prob x vulnerability bivariate choropleth, 2 weeks ahead + legend ---
  p3c_map <- ggplot(mp2) + geom_sf(aes(fill = bikey), colour = "white", linewidth = 0.08) + aff_overlay(mp2) +
    scale_fill_manual(values = BIP, na.value = NA_FILL, guide = "none") + theme_map(15.5)
  legend_biv <- bivar_legend(BIP, "Invasion prob ->", "Vulnerability ->", base = 12)
  p3c <- (p3c_map | wrap_elements(full = legend_biv)) + plot_layout(widths = c(1, 0.42))

  # --- 3D preparedness-priority scatter, 2 weeks ahead ---
  psc <- rs %>% filter(horizon == 2, !(was_active_before %in% TRUE), is.finite(V), is.finite(rr01_nat)) %>%
    mutate(region = ifelse(!is.na(province) & province %in% PROV_INT, province, "Other"))
  psc$region <- factor(psc$region, levels = c(intersect(PROV_INT, unique(psc$region)), "Other"))
  lab_ps <- psc %>% arrange(desc(priority)) %>% head(12)
  p3d <- ggplot(psc, aes(V, rr01_nat)) +
    geom_point(aes(size = priority, colour = region), alpha = 0.78) +
    ggrepel::geom_text_repel(data = lab_ps, aes(label = health_zone), size = 3.9, colour = INK, box.padding = 0.3,
      max.overlaps = 30, seed = 1, min.segment.length = 0, segment.colour = FAINT, segment.size = 0.2) +
    scale_colour_manual(values = PROV_COL, name = "Province", breaks = levels(psc$region), na.value = "#7A7A7A") +
    scale_size_area(max_size = 7, name = "Priority") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = label_number(accuracy = 0.01)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = "Vulnerability", y = "Relative invasion risk") +
    theme_pub(15.5) + theme(legend.position = "right")

  # --- SI narrative pair: front-approach (explanation) + rank evolution (response) ---
  # Paired on a shared time axis / colour / zone set so A mechanistically explains B.
  # Both reuse the province palette `pal` derived from the top-zone forest (main-text Fig 3B),
  # so the narrative panels share identical province hues with it.
  p3fa <- build_front_approach(base = 12, pal = pal)
  p3e  <- build_rank_evolution(base = 12, pal = pal)

  # Individual panels. The former 4-panel Figure 3 is split: the relative-risk map + top-zone
  # forest stay in the MAIN TEXT (now Figure 3 A/B); the front-approach + rank-evolution narrative
  # pair moves to the SI (now Figure S4 A/B). Basenames follow the new figure/panel roles.
  save_dual(p3a,  "F3A_relative_risk_map_h2", 6.4, 4.6)   # main-text Fig 3A
  save_dual(p3b,  "F3B_top_zone_uncertainty", 5.6, 4.6)   # main-text Fig 3B
  save_dual(p3fa, "FS4A_front_approach",      8.4, 5.2)   # SI Fig S4A
  save_dual(p3e,  "FS4B_rank_evolution",      8.4, 5.2)   # SI Fig S4B
  # operational SI panels (prob x vulnerability bivariate + preparedness-priority scatter)
  save_dual(p3c,  "FS2A_prob_vuln_bivariate", 6.4, 4.6)
  save_dual(p3d,  "FS2B_priority_scatter_h2", 6.4, 4.6)

  # MAIN TEXT Figure 3: relative-invasion-risk map (A) + top-zone uncertainty forest (B), 2-week
  # horizon. Distinct legends per panel -> wrap_elements keeps each panel's own guides/layout.
  fig3 <- (wrap_elements(full = p3a) | wrap_elements(full = p3b)) +
    plot_annotation(tag_levels = "A")
  save_dual(fig3, "Figure3", 12.8, 5.2, dir = FIG_DIR)

  # SI Figure S4: front-approach (A) + risk-rank evolution (B) narrative pair. Identical province
  # palette -> patchwork collects ONE shared legend to the right.
  figS4 <- (p3fa | p3e) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") & theme(legend.position = "right")
  save_dual(figS4, "FigureS4_front_rank_evolution", 14.5, 5.6, dir = FIG_DIR)

  # SI Figure S2: prob x vulnerability bivariate (A) + preparedness-priority scatter (B)
  figS2 <- (wrap_elements(full = p3c) | wrap_elements(full = p3d)) +
    plot_annotation(tag_levels = "A")
  save_dual(figS2, "FigureS2_operational_prob_priority", 12.8, 5.2, dir = FIG_DIR)
  invisible(fig3)
}

# -----------------------------------------------------------------------------
run <- function(f, nm) tryCatch({ message("== ", nm, " =="); f(); TRUE },
                                error = function(e){ message("!! ", nm, " FAILED: ", conditionMessage(e)); FALSE })
ok1 <- run(build_fig1, "Figure 1"); okS1 <- run(build_fig1_si, "Figure S1")
okS3 <- run(build_fig_si_epi_share, "Figure S3")
ok2 <- run(build_fig2, "Figure 2"); ok3 <- run(build_fig3, "Figure 3")
okPS <- run(build_phase_space, "Figure 3 phase-space")
message(sprintf("\n[done] Figure1:%s  FigureS1:%s  FigureS3:%s  Figure2:%s  Figure3:%s  PhaseSpace:%s  (featured %s)  ->  %s",
                ok1, okS1, okS3, ok2, ok3, okPS, FEATURED, FIG_DIR))
