# =============================================================================
# make_si_data_figures.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
# Supplementary SURVEILLANCE-STREAM and MOBILITY figures, in the main-text house style
# (figure_style.R): no embedded titles or subtitles, panel tags only.
#
#   FigureS10  Invasion over time, by source   A) the invaded set mapped at four dates, each
#                                              zone coloured by WHICH stream records it
#                                              (both / line list only / sitrep only),
#                                              B) cumulative invaded zones per stream.
#   FigureS11  Line list vs sitrep             A) national cumulative confirmed cases from each
#                                              stream, with the gap shaded, B) the zone-level
#                                              comparison at the snapshot, C) the per-zone gap
#                                              over time for the zones that disagree most.
#   FigureS12  Mobility kernels                A) the import-weight matrices actually used,
#                                              restricted to the zones that carry the mass,
#                                              B) how much the kernels agree with one another,
#                                              C) where each kernel sends the epicentre's
#                                              outflow.
#
# WHY THESE EXIST. The invasion outcome the models are scored against is the UNION of two
# surveillance streams: the DHIS2 line list, and the INSP sitrep, which is used as a FLOOR
# (01_data_prep.R tops the line list up to the sitrep's cumulative confirmed count per zone).
# Four zones on the 2026-09-07 snapshot are invaded on the sitrep alone. That is a modelling
# decision with consequences for every number in the manuscript, and these two figures are
# where it is shown rather than described.
#
# Outputs -> outputs/key_outputs/manuscript_figures/{FigureS10..S12}.{pdf,png} + panels/
# Run:  Rscript make_si_data_figures.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(sf); library(scales); library(forcats)
})
sf::sf_use_s2(TRUE)
options(dplyr.summarise.inform = FALSE)

HERE <- file.path(here::here(), "spatiotemporal")
OUT  <- file.path(HERE, "outputs")
DATA <- file.path(HERE, "..", "data")
source(file.path(HERE, "01_data_prep.R"))          # sources 00_config.R
source(file.path(HERE, "forecast_scale.R"))
FIG_DIR   <- fs_out_dir(file.path(OUT, "key_outputs", "manuscript_figures"))
PANEL_DIR <- fs_out_dir(file.path(OUT, "key_outputs", "manuscript_figures", "panels"))
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
source(file.path(HERE, "figure_style.R"))

SHP_PATH <- file.path(DATA, "shapefiles", "DRC_Health_zones.shp")
key <- function(x) tolower(trimws(x))
shp <- st_read(SHP_PATH, quiet = TRUE) %>% mutate(.key = key(Nom), .prov = as.character(PROVINCE))

ok <- list()
.build <- function(nm, f) {
  ok[[nm]] <<- tryCatch({ f(); TRUE },
    error = function(e) { warning(sprintf("[si] %s FAILED: %s", nm, conditionMessage(e)),
                                  call. = FALSE); FALSE })
}

# -----------------------------------------------------------------------------
# THE TWO STREAMS, ON ONE FOOTING
# -----------------------------------------------------------------------------
# Both are reduced to the same shape — one cumulative confirmed count per (canonical zone,
# date) — so every comparison below is like-for-like.
#
# THE LINE LIST IS TAKEN BEFORE RECONCILIATION. load_linelist() returns the line list with the
# sitrep shortfall already appended (APPEND_SITREP_CONFIRMED), which is correct for the models
# but would make this comparison circular: the reconciled line list agrees with the sitrep by
# construction wherever the sitrep is higher, and the figure would show a gap of zero on
# exactly the zones the reconciliation exists to fix. The appended rows carry alert_ids of the
# form "SITREP-CONF-<zone>-NN", so they are identifiable and are removed here.
ALIASES <- if (file.exists(ALIASES_PATH))
  suppressWarnings(read_csv(ALIASES_PATH, col_types = cols(.default = "c"), show_col_types = FALSE)) else NULL

.streams <- local({
  ll <- load_linelist()
  n_app <- sum(grepl("^SITREP-CONF-", ll$alert_id %||% character(0)))
  ll_raw <- ll %>% filter(!grepl("^SITREP-CONF-", alert_id),
                          (final_mve_case_classification == CONFIRMED_STATUS) %in% TRUE)
  message(sprintf("[si] line list: %d confirmed rows after removing %d sitrep-appended rows",
                  nrow(ll_raw), n_app))
  # The line list is dated by SPECIMEN COLLECTION here, not by imputed onset. The sitrep is a
  # REPORTING series — it records when a confirmation was reported — so comparing it against
  # an onset-dated line list would compare two different clocks and manufacture a lag that is
  # an artefact of the dating convention. The specimen date is the closest thing the line list
  # has to the sitrep's own clock, and it is observed rather than imputed.
  ll_day <- ll_raw %>%
    transmute(zone = health_zone, date = as.Date(date_of_sample_collection)) %>%
    filter(!is.na(zone), !is.na(date), date >= OUTBREAK_START) %>%
    count(zone, date, name = "inc")
  sit_day <- sitrep_cumulative_confirmed(ALIASES) %>%
    filter(inc > 0L) %>% transmute(zone = nom, date, inc = as.integer(inc))
  list(linelist = ll_day, sitrep = sit_day)
})

#' Daily increments -> a complete (zone x date) cumulative series on a shared date spine, so
#' the two streams can be differenced without a join dropping the days on which one of them
#' reported nothing.
.cumulate <- function(d, dates, zones) {
  tidyr::expand_grid(zone = zones, date = dates) %>%
    left_join(d, by = c("zone", "date")) %>%
    mutate(inc = coalesce(inc, 0L)) %>%
    arrange(zone, date) %>% group_by(zone) %>% mutate(cum = cumsum(inc)) %>% ungroup()
}

DATES <- local({
  rg <- range(c(.streams$linelist$date, .streams$sitrep$date))
  seq(rg[1], rg[2], by = "day")
})
ZONES <- sort(unique(c(.streams$linelist$zone, .streams$sitrep$zone)))
CUM   <- bind_rows(
  .cumulate(.streams$linelist, DATES, ZONES) %>% mutate(source = "Line list"),
  .cumulate(.streams$sitrep,   DATES, ZONES) %>% mutate(source = "Sitrep"))
SRC_COL <- c(`Line list` = "#0072B2", `Sitrep` = "#D55E00")

# -----------------------------------------------------------------------------
# FIGURE S10 — where and when each stream says a zone was invaded
# -----------------------------------------------------------------------------
# "Invaded" here means the same thing it means to the models: the zone has at least one
# CONFIRMED case by the date shown. The panel is about WHICH STREAM SUPPORTS THAT CLAIM, so a
# zone is coloured by the pair (line list yes/no, sitrep yes/no) rather than by a single
# status. The two disagreement classes are the interesting ones and are given the two strong
# hues; agreement is the muted one.
INV_LV  <- c("Both streams", "Line list only", "Sitrep only", "Not invaded")
INV_COL <- setNames(c("#5B5B7A", "#0072B2", "#D55E00", "grey90"), INV_LV)

build_figS10 <- function() {
  first_by <- CUM %>% filter(cum > 0) %>% group_by(source, zone) %>%
    summarise(first = min(date), .groups = "drop")

  # Four evenly spaced snapshots across the epidemic, the last pinned to the analysis date so
  # the final map is the state the manuscript actually reports.
  snaps <- as.Date(c(OUTBREAK_START + 30, OUTBREAK_START + 60, OUTBREAK_START + 90,
                     ANALYSIS_DATE))
  snaps <- sort(unique(snaps[snaps <= ANALYSIS_DATE]))

  # Map extent: the affected provinces, not the whole country — at national extent the
  # invaded zones are a few pixels and the panel says nothing.
  PROV_SHOW <- c("Ituri", "Nord-Kivu", "Haut-Uele", "Tshopo", "Bas-Uele")
  reg <- shp %>% filter(.prov %in% PROV_SHOW)
  if (!nrow(reg)) stop("no shapefile rows for the affected provinces")
  # EXTENT: the invaded set plus a margin, not the provinces' full footprint. Tshopo and
  # Bas-Uele run several hundred kilometres west of anything that has been invaded, and
  # framing on them shrinks the epidemic to a few pixels. The margin keeps the uninvaded
  # frontier — the zones the forecast is about — inside the frame.
  bb <- local({
    inv_keys <- unique(key(first_by$zone))
    g <- reg[reg$.key %in% inv_keys, ]
    b <- st_bbox(if (nrow(g)) g else reg)
    mx <- 0.55 * as.numeric(b["xmax"] - b["xmin"]); my <- 0.25 * as.numeric(b["ymax"] - b["ymin"])
    rb <- st_bbox(reg)
    c(xmin = max(b[["xmin"]] - mx, rb[["xmin"]]), xmax = min(b[["xmax"]] + mx, rb[["xmax"]]),
      ymin = max(b[["ymin"]] - my, rb[["ymin"]]), ymax = min(b[["ymax"]] + my, rb[["ymax"]]))
  })

  # THE STATUS TABLE COVERS EVERY ZONE ON THE MAP, not only the zones one of the streams
  # mentions. Building it over the reported zones alone and joining left the rest with no
  # status; they dropped out of the facet and were drawn as blank paper, so "not invaded"
  # was indistinguishable from "outside the map". Keyed on the shapefile's own `.key`, which
  # is what the join uses, so a zone cannot be absent from one snapshot and present in another.
  reg_keys <- unique(reg$.key)
  lk_key   <- first_by %>% mutate(.key = key(zone))
  status_at <- function(t) {
    l <- lk_key$.key[lk_key$source == "Line list" & lk_key$first <= t]
    s <- lk_key$.key[lk_key$source == "Sitrep"    & lk_key$first <= t]
    tibble(.key = reg_keys,
           status = factor(case_when(.key %in% l & .key %in% s ~ INV_LV[1],
                                     .key %in% l               ~ INV_LV[2],
                                     .key %in% s               ~ INV_LV[3],
                                     TRUE                      ~ INV_LV[4]),
                           levels = INV_LV),
           snap = t)
  }
  st <- bind_rows(lapply(snaps, status_at))

  mp <- reg %>%
    left_join(st, by = ".key", relationship = "many-to-many") %>%
    mutate(snap_lab = factor(format(snap, "%d %b %Y"),
                             levels = format(snaps, "%d %b %Y")))
  if (anyNA(mp$status))
    stop("some mapped zones have no invasion status; the status table and the shapefile ",
         "disagree on the zone key set.")
  prov_bounds <- reg %>% group_by(.prov) %>% summarise(.groups = "drop") %>% st_make_valid()

  pA <- ggplot(mp) +
    geom_sf(aes(fill = status), colour = "white", linewidth = 0.06) +
    geom_sf(data = prov_bounds, fill = NA, colour = "grey35", linewidth = 0.25) +
    facet_wrap(~ snap_lab, nrow = 1) +
    scale_fill_manual(values = INV_COL, name = NULL, drop = FALSE) +
    coord_sf(xlim = c(bb[["xmin"]], bb[["xmax"]]), ylim = c(bb[["ymin"]], bb[["ymax"]]),
             expand = FALSE) +
    theme_map(9.5) +
    theme(legend.position = "top", legend.key.width = unit(11, "pt"),
          legend.key.height = unit(11, "pt"),
          strip.text = element_text(size = 9, colour = INK, face = "bold",
                                    margin = margin(2, 2, 3, 2)),
          panel.border = element_rect(fill = NA, colour = GRID, linewidth = 0.3))

  # --- B: cumulative invaded-zone count per stream, plus the union the models are scored on ---
  curve <- bind_rows(lapply(DATES, function(t) {
    l <- first_by$zone[first_by$source == "Line list" & first_by$first <= t]
    s <- first_by$zone[first_by$source == "Sitrep"    & first_by$first <= t]
    tibble(date = t,
           `Line list` = length(l), `Sitrep` = length(s),
           `Union (used by the models)` = length(union(l, s)))
  })) %>% pivot_longer(-date, names_to = "series", values_to = "n")
  SER_COL <- c(`Line list` = "#0072B2", `Sitrep` = "#D55E00",
               `Union (used by the models)` = "grey20")
  pB <- ggplot(curve, aes(date, n, colour = series, linetype = series)) +
    geom_vline(xintercept = snaps, colour = FAINT, linetype = "22", linewidth = 0.3) +
    geom_step(linewidth = 0.85) +
    scale_colour_manual(values = SER_COL, name = NULL) +
    scale_linetype_manual(values = c(`Line list` = "solid", `Sitrep` = "solid",
                                     `Union (used by the models)` = "22"), name = NULL) +
    scale_x_date(date_breaks = "3 weeks", date_labels = "%d %b",
                 expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = "Date of first confirmed case (specimen date; sitrep report date)",
         y = "Health zones invaded") +
    theme_pub(9.5) + theme(legend.position = "top")

  save_dual(pA, "FS10A_invaded_maps_by_source", 11.0, 3.6)
  save_dual(pB, "FS10B_invaded_count_by_source", 6.4, 3.6)
  fig <- (pA / pB) + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureS10_invasion_by_source", 11.6, 7.4, dir = FIG_DIR)
  invisible(fig)
}

# -----------------------------------------------------------------------------
# FIGURE S11 — how far apart the two streams are, and where
# -----------------------------------------------------------------------------
# The reconciliation takes the sitrep as a FLOOR: wherever the sitrep's cumulative confirmed
# count for a zone exceeds the line list's, the difference is appended as confirmed cases. The
# size of that operation is the subject of this figure. Panel A is the national total, panel B
# is the per-zone position at the snapshot, panel C is where the gap sits over time.
build_figS11 <- function() {
  nat <- CUM %>% group_by(source, date) %>% summarise(cum = sum(cum), .groups = "drop")
  gap <- nat %>% pivot_wider(names_from = source, values_from = cum) %>%
    mutate(lo = pmin(`Line list`, Sitrep), hi = pmax(`Line list`, Sitrep))
  pA <- ggplot(nat, aes(date, cum, colour = source)) +
    geom_ribbon(data = gap, aes(x = date, ymin = lo, ymax = hi), fill = "grey65", alpha = 0.22,
                colour = NA, inherit.aes = FALSE) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = SRC_COL, name = NULL) +
    scale_x_date(date_breaks = "3 weeks", date_labels = "%d %b",
                 expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = comma, limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = NULL, y = "Cumulative confirmed cases (national)") +
    theme_pub(9.5) + theme(legend.position = "top")

  # --- C: the per-zone gap over time, for the zones that disagree most ---
  # WEEKLY, not daily: the sitrep is published irregularly, so a daily heatmap is mostly an
  # image of its publication calendar. The value shown is the gap at the END of each week.
  wk <- CUM %>%
    mutate(week = as.Date(cut(date, "week", start.on.monday = TRUE)) + 6L) %>%
    group_by(source, zone, week) %>% summarise(cum = max(cum), .groups = "drop") %>%
    pivot_wider(names_from = source, values_from = cum) %>%
    mutate(gap = Sitrep - `Line list`)
  topz <- wk %>% group_by(zone) %>% summarise(m = max(abs(gap), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>% head(20) %>% pull(zone)
  hb <- wk %>% filter(zone %in% topz) %>%
    mutate(zone = factor(zone, levels = rev(topz)))
  lim <- max(abs(hb$gap), na.rm = TRUE)
  pC <- ggplot(hb, aes(week, zone, fill = gap)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = SRC_COL[["Line list"]], mid = "grey96",
                         high = SRC_COL[["Sitrep"]], midpoint = 0, limits = c(-lim, lim),
                         name = "Sitrep minus line list\n(cumulative confirmed)") +
    scale_x_date(date_breaks = "3 weeks", date_labels = "%d %b", expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_pub(9.5) + theme(panel.grid = element_blank(), legend.position = "right",
                           legend.key.width = unit(8, "pt"), legend.key.height = unit(26, "pt"),
                           axis.text.y = element_text(size = 7.2, colour = INK))

  # --- B: the zone-level comparison at the snapshot ---
  fin <- CUM %>% filter(date == max(DATES)) %>%
    select(zone, source, cum) %>% pivot_wider(names_from = source, values_from = cum) %>%
    filter(`Line list` > 0 | Sitrep > 0)
  # log1p on both axes: the counts span three orders of magnitude and the zeros must stay
  # visible (a zone confirmed by one stream and not the other is the whole point).
  lab <- fin %>% mutate(d = abs(Sitrep - `Line list`)) %>% arrange(desc(d)) %>% head(12)
  pB <- ggplot(fin, aes(`Line list` + 1, Sitrep + 1)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = FAINT, linewidth = 0.4) +
    geom_point(colour = PT_BLUE, size = 1.9, alpha = 0.85) +
    ggrepel::geom_text_repel(data = lab, aes(label = zone), size = 2.5, colour = INK,
                             min.segment.length = 0, segment.colour = FAINT,
                             segment.size = 0.2, max.overlaps = 30, seed = 1) +
    # Explicit breaks on the shifted log axis. The default decade breaks land at 1, 10, 100,
    # 1000 and therefore LABEL 0, 9, 99, 999 once the +1 shift is undone, which reads as a
    # mistake; these are the round numbers the reader expects.
    scale_x_continuous(trans = "log10", breaks = 1 + c(0, 1, 3, 10, 30, 100, 300, 1000),
                       labels = function(v) comma(v - 1)) +
    scale_y_continuous(trans = "log10", breaks = 1 + c(0, 1, 3, 10, 30, 100, 300, 1000),
                       labels = function(v) comma(v - 1)) +
    labs(x = "Line-list confirmed cases (unreconciled)", y = "Sitrep confirmed cases") +
    coord_equal() +
    theme_pub(9.5)

  save_dual(pA, "FS11A_national_streams", 6.0, 3.4)
  save_dual(pB, "FS11B_zone_scatter", 4.6, 4.4)
  save_dual(pC, "FS11C_zone_gap_heatmap", 7.6, 4.6)
  # Layout order IS tag order: A national totals, B the zone scatter, C the gap over time.
  fig <- (pA | pB) / pC + plot_layout(heights = c(1, 1.15)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureS11_linelist_vs_sitrep", 11.6, 8.2, dir = FIG_DIR)
  invisible(fig)
}

# -----------------------------------------------------------------------------
# FIGURE S12 — the mobility kernels the models actually use
# -----------------------------------------------------------------------------
# Every invasion model is the SAME hazard with a DIFFERENT import-weight matrix W, so W is the
# scientific content of the model comparison and it has never been shown. Each W is
# row-stochastic (row i = where zone i's outward movement goes) with a zero diagonal.
#
# A 519 x 519 heatmap is unreadable, and thresholding it would show the kernels' shared
# geography rather than what separates them. The panels are therefore chosen to answer three
# specific questions:
#   A  WHAT DOES W LOOK LIKE? The submatrix over the zones that carry the mass — every invaded
#      zone plus the highest-weight destinations — ordered by province so the block structure
#      is interpretable, on a log colour scale because the weights span six orders of magnitude.
#   B  HOW MUCH DO THE KERNELS AGREE? Spearman correlation between the off-diagonal entries of
#      each pair. Rank correlation, not Pearson: the hazard uses W through a sum, but what
#      separates the kernels' PREDICTIONS is the ordering of destinations, and Pearson on
#      six-order-of-magnitude weights is a correlation between the largest few cells.
#   C  WHERE DOES THE EPICENTRE'S OUTFLOW GO? The row of W out of the M+B+R epicentre, which is
#      the single most consequential row in the whole matrix during the seeding phase.
build_figS12 <- function() {
  lfo <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))
  kern <- sort(unique(stats::na.omit(vapply(unique(lfo$method), mobility_kernel_from_method,
                                            character(1)))))
  if (!length(kern)) stop("no mobility kernel could be parsed from the cross-validated methods")
  Ws <- setNames(lapply(kern, function(k) {
    f <- file.path(OUT, "mobility", sprintf("mobility_%s.rds", k))
    if (!file.exists(f)) { warning("[si] kernel matrix missing: ", basename(f), call. = FALSE)
                           return(NULL) }
    readRDS(f)
  }), kern)
  Ws <- Filter(Negate(is.null), Ws)
  if (!length(Ws)) stop("none of the cross-validated kernels' matrices are on disk")
  zn <- rownames(Ws[[1]])
  if (!all(vapply(Ws, function(W) identical(rownames(W), zn), logical(1))))
    stop("the kernel matrices do not share a zone ordering; they cannot be compared cell by cell")
  klab <- setNames(unname(MOBILITY_LABELS[names(Ws)]), names(Ws))
  klab[is.na(klab)] <- names(Ws)[is.na(klab)]
  # Compact facet titles: the long MOBILITY_LABELS glossary entries do not fit a facet strip.
  kshort <- setNames(model_pretty_label(paste0("Bayes-", names(Ws), "-med"), warn_unknown = FALSE),
                     names(Ws))
  kshort[kshort == paste0("Bayes-", names(Ws), "-med")] <- names(Ws)[kshort == paste0("Bayes-", names(Ws), "-med")]

  # --- zone subset: every invaded zone, plus the strongest destinations overall ---
  rs  <- read_csv(fs_risk_csv(OUT), show_col_types = FALSE)
  inv <- rs %>% filter(as.logical(was_active_before) %in% TRUE) %>% pull(health_zone) %>% unique()
  inv <- intersect(inv, zn)
  Wsum <- Reduce(`+`, Ws) / length(Ws)
  # Inflow TO each zone from the invaded set — the quantity the import force is built from.
  inflow <- if (length(inv)) colSums(Wsum[inv, , drop = FALSE]) else colSums(Wsum)
  extra <- setdiff(names(sort(inflow, decreasing = TRUE)), inv)
  sel <- unique(c(inv, head(extra, max(0L, 55L - length(inv)))))
  prov <- rs %>% distinct(health_zone, province)
  ord <- tibble(zone = sel) %>% left_join(prov, by = c("zone" = "health_zone")) %>%
    mutate(province = coalesce(province, "Other")) %>%
    arrange(factor(province, levels = c("Ituri", "Nord-Kivu", "Haut-Uele", "Tshopo", "Bas-Uele")),
            province, zone)
  sel <- ord$zone

  # --- A: the submatrices ---
  FLOOR <- 1e-6      # weights below this are drawn as the floor; they cannot seed anything
  hm <- bind_rows(lapply(names(Ws), function(k) {
    M <- Ws[[k]][sel, sel, drop = FALSE]
    as.data.frame(as.table(M)) %>%
      setNames(c("from", "to", "w")) %>%
      mutate(kernel = kshort[[k]])
  })) %>%
    mutate(from = factor(as.character(from), levels = sel),
           to   = factor(as.character(to),   levels = rev(sel)),
           w    = pmax(w, FLOOR),
           kernel = factor(kernel, levels = unname(kshort[names(Ws)])))
  pA <- ggplot(hm, aes(from, to, fill = w)) +
    geom_raster() +
    # Wrapped strip labels: the kernel names run to ~40 characters and a 2-inch facet strip
    # clipped them mid-word ("hort-trip + radiation (source-fille").
    facet_wrap(~ kernel, nrow = 2, labeller = labeller(kernel = label_wrap_gen(width = 22))) +
    scale_fill_viridis_c(trans = "log10", option = "magma", direction = -1,
                         name = "Import weight\nW[from, to]",
                         labels = label_number(accuracy = 0.000001, drop0trailing = TRUE),
                         na.value = "white") +
    labs(x = "Origin zone", y = "Destination zone") +
    coord_equal() +
    theme_pub(8.5) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "right",
          legend.key.width = unit(8, "pt"), legend.key.height = unit(26, "pt"),
          strip.text = element_text(size = 6.6, lineheight = 0.95,
                                    margin = margin(2, 2, 2, 2)))

  # --- B: agreement between kernels (Spearman on the off-diagonal cells) ---
  offd <- which(row(Ws[[1]]) != col(Ws[[1]]))
  V <- vapply(Ws, function(W) as.numeric(W)[offd], numeric(length(offd)))
  CM <- suppressWarnings(stats::cor(V, method = "spearman"))
  cdf <- as.data.frame(as.table(CM)) %>% setNames(c("a", "b", "rho")) %>%
    mutate(a = factor(kshort[as.character(a)], levels = unname(kshort[names(Ws)])),
           b = factor(kshort[as.character(b)], levels = rev(unname(kshort[names(Ws)]))))
  pB <- ggplot(cdf, aes(a, b, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 2.1,
              colour = ifelse(cdf$rho > 0.7, "white", INK)) +
    scale_fill_gradient(low = "grey95", high = "#08519C", limits = c(0, 1),
                        name = "Spearman rho\n(off-diagonal cells)") +
    labs(x = NULL, y = NULL) + coord_equal() +
    theme_pub(8.5) + theme(panel.grid = element_blank(),
                           axis.text.x = element_text(angle = 35, hjust = 1, size = 6.6),
                           axis.text.y = element_text(size = 6.6),
                           legend.position = "right", legend.key.width = unit(8, "pt"))

  # --- C: where the epicentre's outflow goes ---
  epi <- intersect(EPICENTRE_ZONES, zn)
  if (!length(epi)) stop("no epicentre zone is present in the kernel matrices")
  # Population-weighted mean of the epicentre rows, so the three seed zones are combined the
  # way the import force combines them (each row is weighted by the zone's case load, which at
  # the seeding phase is dominated by the epicentre itself). Equal weights would let the
  # smallest seed zone speak as loudly as Bunia.
  wrow <- local({
    cw <- rs %>% filter(health_zone %in% epi) %>%
      transmute(health_zone, w = 1) %>% tibble::deframe()
    cw <- cw[epi]; cw[!is.finite(cw)] <- 1; cw / sum(cw)
  })
  outfl <- bind_rows(lapply(names(Ws), function(k) {
    r <- colSums(Ws[[k]][epi, , drop = FALSE] * as.numeric(wrow))
    tibble(kernel = kshort[[k]], zone = names(r), w = as.numeric(r))
  }))
  topd <- outfl %>% group_by(zone) %>% summarise(m = max(w), .groups = "drop") %>%
    arrange(desc(m)) %>% head(22) %>% pull(zone)
  oc <- outfl %>% filter(zone %in% topd) %>%
    mutate(zone = factor(zone, levels = rev(topd)),
           kernel = factor(kernel, levels = unname(kshort[names(Ws)])))
  pC <- ggplot(oc, aes(pmax(w, FLOOR), zone, colour = kernel)) +
    geom_point(size = 1.6, alpha = 0.9) +
    scale_x_continuous(trans = "log10", labels = label_number(accuracy = 0.0001, drop0trailing = TRUE)) +
    scale_colour_manual(values = setNames(
      grDevices::colorRampPalette(OKABE[1:6])(length(Ws)), unname(kshort[names(Ws)])),
      name = NULL) +
    labs(x = "Share of epicentre outflow (log scale)", y = NULL) +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2.4))) +
    theme_pub(8.5) + theme(panel.grid.major.y = element_line(colour = GRID, linewidth = 0.25),
                           legend.position = "right",
                           axis.text.y = element_text(size = 6.8, colour = INK))

  save_dual(pA, "FS12A_kernel_matrices", 9.6, 5.6)
  save_dual(pB, "FS12B_kernel_agreement", 5.0, 4.4)
  save_dual(pC, "FS12C_epicentre_outflow", 6.4, 4.6)
  # heights: panel A is two rows of SQUARE matrices, so its natural height is set by the
  # panel width; allocating more than that leaves a band of blank paper between A and B.
  fig <- pA / (pB | pC) + plot_layout(heights = c(1, 1.05)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureS12_mobility_kernels", 12.0, 9.2, dir = FIG_DIR)
  invisible(fig)
}

# -----------------------------------------------------------------------------
if (!isTRUE(get0(".SI_DATA_FIGURES_NO_RUN", ifnotfound = FALSE))) {
  .build("FigureS10", build_figS10)
  .build("FigureS11", build_figS11)
  .build("FigureS12", build_figS12)
  message(sprintf("\n[si-data] %s  ->  %s",
                  paste(sprintf("%s:%s", names(ok), ifelse(unlist(ok), "ok", "FAILED")),
                        collapse = "  "), FIG_DIR))
}
