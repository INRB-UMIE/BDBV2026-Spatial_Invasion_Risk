# =============================================================================
# 40_cascade_next_dominoes.R — CASCADE: the "next dominoes" narrative figure
# BDBV 2026 DRC · 3-month spatial invasion cascade
#
# ONE coherent figure answering "if a major city is invaded, where do we look
# next?" — combining the eastern-epicentre frontier (Figure4_cascade anchor) with
# the urban conditional scenarios (Figure_urban_*), across all five hub cities.
#
#   Panel A  HUB FLOW MAP — each hub city's mobility to/from ANYWHERE in the DRC (the
#            same hubs as panels B & C), from the Flowminder directed origin->destination
#            table. Edges link each labelled hub to every zone it exchanges people with,
#            coloured by DIRECTION (outflow vs inflow), width = relocations. "Where does each
#            hub pull from and push to across the country?"
#   Panel B  WHICH domino matters — cities ranked by expected additional health
#            zones invaded within 13 weeks if that city falls (seeded week 1).
#   Panel C  WHERE to look next — per-city watch-lists (one facet per city). Each
#            zone is ranked by the ATTRIBUTABLE increase it receives from the city
#            (so the list is the city's OWN catchment, not the pre-saturated east),
#            and drawn as a baseline -> conditional dumbbell so BOTH quantities show:
#              • dot position  = TOTAL P(invasion | city invaded)   [where to look]
#              • segment length = ATTRIBUTABLE increase P_cond - P_base [what the
#                                 city itself adds]
#
# Rationale for showing BOTH (was: attributable-only): attributable Δ answers the
# atttribution question ("does this city matter?") but is pale wherever baseline
# is already high; the operational "where do I send teams" question needs the TOTAL
# conditional risk. The dumbbell carries both in a single mark.
#
# Standalone & regenerable: reads only saved CSVs + the Flowminder OD matrix (no
# re-simulation), mirroring 37_cascade_figure4.R. Works in both spatiotemporal/ and
# spatiotemporal_conditional/.
#
# Run:  Rscript spatiotemporal/40_cascade_next_dominoes.R
# Out:  outputs/key_outputs/manuscript_figures/Figure4.{pdf,png}
#       outputs/key_outputs/Figure_next_dominoes_watchlist.csv
# Panel titles/subtitles and the embedded caption are intentionally OMITTED — figure
# identity is carried by the A/B/C tags and axis titles; the caption lives in the manuscript.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(patchwork)
  library(sf); library(scales); library(here)
})
ST_DIR <- Sys.getenv("CASCADE_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
source(file.path(ST_DIR, "00_config.R"))

# ---- style constants (shared with 37_cascade_figure4.R) --------------------
INK <- "grey15"; MUTED <- "grey38"; GRID <- "grey92"
OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
# Province identity spanning the eastern epicentre AND the western/southern urban
# catchments (Okabe-Ito, CVD-safe). Merged from 37 (east) + 38 (urban).
PROV_COL <- c(
  "Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73","Sud-Kivu"="#CC79A7",
  "Kinshasa"="#E69F00","Haut-Katanga"="#56B4E9","Kasaï-Oriental"="#CC79A7",
  "Kasaï-Central"="#009E73","Kasaï"="#D55E00","Kwilu"="#0072B2",
  "Kwango"="#F0E442","Lomami"="#000000","Mai-Ndombe"="#7A7A7A")
# Spare palette for provinces WITHOUT an identity colour (e.g. peripheral zones a hub's
# catchment reaches — Tshopo, Maniema). PROV_COL already consumes every OKABE hue, so
# setdiff(OKABE, PROV_COL) is empty; without a real fallback such provinces map to NA and
# their Panel C dots are silently dropped. These Paul-Tol muted hues (CVD-safe) are disjoint
# from PROV_COL so every province gets a distinct, visible colour.
SPARE_COL <- c("#882255","#332288","#661100","#999933","#44AA99","#AA4499","#DDCC77","#117733")
PROV_INT <- get0("PROVINCES_OF_INTEREST", ifnotfound = c("Ituri","Nord-Kivu","Haut-Uele"))
base_family <- "sans"
HORIZON <- max(get0("CASCADE_REPORT_HORIZONS", ifnotfound = c(4L, 8L, 13L)))
TOP_N   <- 5L                                    # max zones per city watch-list (Panel C)

# ---- display labels: city vs health zone of the same name -------------------
# "Kinshasa" is BOTH a hub city and a health zone inside it, and the health zone ranks
# 4th in the Kinshasa watch-list — so Panel C's facet strip and one of its own rows both
# read "Kinshasa", which is genuinely ambiguous to a reader. Disambiguate at the DISPLAY
# layer only: `city` and `health_zone` are join keys (against the urban impact CSVs, the
# seed-zone provenance table and CITY_ORDER), so the stored values must NOT be renamed.
# city_label() is applied wherever a CITY is named (Panels A, B and the Panel C facet
# strips) so the figure reads consistently; zone_label() only to health-zone axis text.
CITY_DISPLAY <- c(Kinshasa = "Kinshasa City")
ZONE_DISPLAY <- c(Kinshasa = "Kinshasa (HZ)")
.relabel <- function(x, map) { x <- as.character(x); i <- match(x, names(map))
                               ifelse(is.na(i), x, unname(map[i])) }
city_label <- function(x) .relabel(x, CITY_DISPLAY)
zone_label <- function(x) .relabel(x, ZONE_DISPLAY)

# Urban-core health zones per hub city — used by Panel A to aggregate mobility in/out
# of each hub as a unit. These MUST be the zones the cascade actually force-seeded, so
# they are read from the seed-zone provenance CSV that 38_urban_scenarios.R writes
# alongside the scenario outputs (see urban_hub_selection.csv for the derivation rule
# and every rejected candidate). The frozen literal below is only a last-resort fallback
# for an output tree written before that CSV existed; it is flagged loudly when used.
URBAN_HUBS_FALLBACK <- list(
  Kinshasa     = c("Gombe", "Limete", "Kimbanseke"),
  Lubumbashi   = c("Lubumbashi", "Kampemba", "Ruashi"),
  `Mbuji-Mayi` = c("Diulu", "Bonzola", "Nzaba"),
  Kananga      = c("Kananga", "Katoka", "Lukonga"),
  Tshikapa     = c("Tshikapa", "Kanzala"))
read_urban_hubs <- function(tbl_dir) {
  f <- file.path(tbl_dir, "urban_hub_selection.csv")
  if (!file.exists(f)) {
    warning("[dominoes] ", basename(f), " not found — falling back to the frozen hub ",
            "list; Panel A may not match the seeded zones. Re-run run_cascade.R.",
            call. = FALSE)
    return(URBAN_HUBS_FALLBACK)
  }
  s <- readr::read_csv(f, show_col_types = FALSE)
  flag <- if ("seeded" %in% names(s)) "seeded" else "selected"
  s <- s[which(s[[flag]]), ]
  out <- split(s$Nom, s$city)
  out[order(match(names(out), names(URBAN_HUBS_FALLBACK)))]
}
# MATERIALITY. The impact tables this figure reads are now written by a PAIRED contrast
# (common random numbers, 38_urban_scenarios.R), so "materially elevated" is a test — the
# zone's 90% paired interval excludes zero — carried in the `elevated` column. That test
# replaces the fixed cut below, which was never a noise floor nor a pre-registered effect
# size: paired, the per-zone standard error is ~0.0006, so a 0.03 cut would discard real
# effects fifty times larger than the noise. The cut survives only as the fallback for
# tables written before the pairing existed, and must keep matching 38's URBAN_DELTA_MIN —
# 38 is not sourced by this standalone script, so it is read if present.
# get0() ALWAYS MISSES here — this script sources 00_config.R only, never
# 38_urban_scenarios.R, which is where URBAN_DELTA_MIN is defined — so the literal below is
# what runs, and the comment's promise that it "is read if present" was never kept. Source the
# definition so the two genuinely cannot drift; fall back to the literal only if 38 is absent.
if (!exists("URBAN_DELTA_MIN")) {
  .u38 <- file.path(ST_DIR, "38_urban_scenarios.R")
  if (file.exists(.u38)) {
    .e38 <- new.env()
    try(suppressWarnings(suppressMessages(
      eval(parse(text = grep("^URBAN_DELTA_MIN\\s*<-", readLines(.u38, warn = FALSE), value = TRUE)),
           envir = .e38))), silent = TRUE)
    if (exists("URBAN_DELTA_MIN", envir = .e38))
      URBAN_DELTA_MIN <- get("URBAN_DELTA_MIN", envir = .e38)
  }
}
DELTA_MIN <- get0("URBAN_DELTA_MIN", ifnotfound = 0.03)
if (!exists("URBAN_DELTA_MIN"))
  warning("[dominoes] URBAN_DELTA_MIN not found in 38_urban_scenarios.R; using the literal ",
          DELTA_MIN, ", which can drift from the urban tables.", call. = FALSE)

# top materially-elevated zones for one city, ranked by attributable increase (delta).
# Shared by the figure and the CSV so both stay identical.
city_watchlist <- function(df, top_n = TOP_N) {
  keep <- if ("elevated" %in% names(df)) df$elevated %in% TRUE
          else is.finite(df$delta) & df$delta > DELTA_MIN
  df[keep & is.finite(df$delta), , drop = FALSE] |>
    dplyr::slice_max(order_by = delta, n = top_n, with_ties = FALSE)
}

# assign a stable colour to every province (identity where known, spare Tol hue else)
prov_palette <- function(provs) {
  provs <- unique(as.character(provs[!is.na(provs)]))
  provs <- c(intersect(PROV_INT, provs), sort(setdiff(provs, PROV_INT)))
  # Fallback pool disjoint from the identity colours; recycled if ever exhausted.
  spare <- setdiff(SPARE_COL, unname(PROV_COL)); si <- 1L
  cols <- setNames(character(length(provs)), provs)
  for (p in provs) {
    if (p %in% names(PROV_COL)) cols[p] <- unname(PROV_COL[p])
    else { cols[p] <- spare[((si - 1L) %% length(spare)) + 1L]; si <- si + 1L }
  }
  stopifnot(!anyNA(cols), all(nzchar(cols)))   # every province must get a visible colour
  cols
}

# House themes (verbatim from 37_cascade_figure4.R / make_publication_figures.R):
# embedded titles/subtitles/captions are BLANKED — figure identity is carried by the
# A/B/C tags and axis titles; the caption lives in the manuscript, not in the figure.
theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    axis.title = element_text(size = base - 0.4, colour = MUTED),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4), angle = 90),
    axis.text = element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
    strip.text = element_text(size = base - 0.6, colour = INK, face = "bold"),
    legend.position = "top", legend.justification = "left",
    legend.title = element_text(size = base - 1.2, colour = MUTED),
    legend.text = element_text(size = base - 1.4, colour = INK),
    legend.key.height = unit(9, "pt"), legend.key.width = unit(15, "pt"),
    plot.tag = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin = margin(6, 8, 6, 6))
}
theme_map <- function(base = 8.6) {
  theme_void(base_size = base, base_family = base_family) %+replace% theme(
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    legend.position = "right", legend.justification = "centre",
    legend.title = element_text(size = base - 1.4, colour = MUTED),
    legend.text = element_text(size = base - 1.8, colour = INK),
    legend.key.height = unit(12, "pt"), legend.key.width = unit(14, "pt"),
    plot.tag = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin = margin(2, 2, 2, 2))
}
save_dual <- function(p, name, w, h, dir) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  # Gate on the FULL destination path, not the bare stem: FIGURE_DROP entries are
  # "<directory>/<stem>" and the raw/ exclusion inspects path components, neither of
  # which can match a basename.
  if (is.function(.fk) && !.fk(file.path(dir, name))) return(invisible(p))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  message(sprintf("  saved %-30s  %.1f x %.1f in", name, w, h)); invisible(p)
}

# ---- paths / inputs --------------------------------------------------------
FIG_DIR   <- file.path(OUT_DIR, "key_outputs", "figures")
MANU_DIR  <- file.path(OUT_DIR, "key_outputs", "manuscript_figures")   # main-text Figure 4 target
KEY_DIR   <- file.path(OUT_DIR, "key_outputs")
TBL_DIR   <- file.path(OUT_DIR, "cascade", "tables")
URBAN_HUBS <- read_urban_hubs(TBL_DIR)     # the zones the cascade actually seeded
reach_csv <- file.path(TBL_DIR, "cascade_reach_scores_all_zones.csv")
summ_csv  <- file.path(TBL_DIR, "urban_scenario_summary.csv")
stopifnot(file.exists(reach_csv), file.exists(summ_csv))

rs   <- read_csv(reach_csv, show_col_types = FALSE)
summ <- read_csv(summ_csv,  show_col_types = FALSE)

# per-city impact CSVs (p_base, p_cond, delta, rr, province) — one row per zone
city_files <- list.files(TBL_DIR, pattern = "^urban_impact_.*\\.csv$", full.names = TRUE)
city_files <- city_files[!grepl("summary", city_files)]
imp <- bind_rows(lapply(city_files, function(f)
  read_csv(f, show_col_types = FALSE))) |>
  filter(is.finite(p_base), is.finite(p_cond))
stopifnot(nrow(imp) > 0)

# city ordering by impact (expected additional zones, earliest seed = max impact)
min_wk <- min(summ$seed_week, na.rm = TRUE)
# The summary now carries one row per (city, seed_week, seeded-R arm). The pooled arm is
# the primary specification; without this filter every city would appear once per arm.
sumw   <- summ |>
  { \(d) if ("seed_r" %in% names(d)) filter(d, is.na(seed_r)) else d }() |>
  filter(seed_week == min_wk) |> arrange(desc(expected_added))
CITY_ORDER <- sumw$city

# shared province palette across every panel
all_provs <- unique(c(imp$province, rs$province[rs$horizon == HORIZON]))
PCOL <- prov_palette(all_provs)

# =============================================================================
# Panel A — HUB FLOW MAP: each city's mobility to/from ANYWHERE in the DRC
# =============================================================================
# For each of the five hub cities (aggregated over its constituent zones), every edge
# links the hub to another DRC health zone it exchanges people with: OUTFLOW (hub ->
# zone) and INFLOW (zone -> hub), coloured by DIRECTION (not by city), width = relocations
# (Flowminder ESTIMATED RELOCATIONS, monthly home-location changes — NOT trips).
# The two directions of a hub<->partner pair are drawn as opposite-bowing arcs so both
# stay visible. Read off the Flowminder directed OD table O[origin_row, dest_col]; its
# processed "inflow"/"outflow" exports are byte-identical (03_mobility_matrices.R::
# build_M15), so the ONE table gives outflow (hub rows) and inflow (hub columns).
load_flow_od <- function() {
  # Follow the configured export (00_config.R FLOWMINDER_OD_FILE) so the sensitivity arm
  # and the model kernels read the SAME month.
  f <- file.path(FLOWMINDER_DIR,
                 get0("FLOWMINDER_OD_FILE",
                      ifnotfound = "flowminder__outflow__static.matrix.csv"))
  stopifnot(file.exists(f))
  df <- readr::read_csv(f, col_types = readr::cols(.default = "d", nom = "c"),
                        show_col_types = FALSE)
  M <- as.matrix(df[, setdiff(colnames(df), "nom")]); rownames(M) <- df$nom
  M
}

# Manual label offsets (degrees); Kananga & Mbuji-Mayi sit close, so pull them apart.
HUB_LABEL_NUDGE <- list(
  Kinshasa = c(-1.5, 1.1), Lubumbashi = c(0.0, -1.3), `Mbuji-Mayi` = c(1.7, 0.9),
  Kananga = c(-1.8, -0.7), Tshikapa = c(-0.6, -1.2))
FLOW_COL  <- c(Outflow = "#B2182B", Inflow = "#2166AC")  # direction-only palette (not city)
FLOW_MIN  <- 5L    # drop negligible edges (< this many relocations) so the map is not a hairball

build_hubflowmap <- function() {
  M      <- load_flow_od()
  cities <- names(URBAN_HUBS)
  present <- vapply(URBAN_HUBS, function(hz)
    all(hz %in% rownames(M)) && all(hz %in% colnames(M)), logical(1))
  if (!all(present)) warning("[flowmap] hub zone(s) absent from OD table: ",
                             paste(cities[!present], collapse = ", "))

  shp <- sf::st_read(SHAPEFILE_PATH, quiet = TRUE)
  shp$.key <- tolower(trimws(shp$Nom))
  # Centroid of every health zone (WGS84), de-duplicated on the name key used to join OD.
  zc <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(shp))))
  zc <- data.frame(key = shp$.key, lon = zc[, 1], lat = zc[, 2], stringsAsFactors = FALSE)
  zc <- zc[!duplicated(zc$key), ]
  # Hub point = centroid of the UNION of its constituent zones.
  ctr <- do.call(rbind, lapply(cities, function(nm) {
    g  <- shp[shp$.key %in% tolower(trimws(URBAN_HUBS[[nm]])), ]
    xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_union(sf::st_geometry(g)))))
    data.frame(city = nm, lon = xy[1, 1], lat = xy[1, 2], stringsAsFactors = FALSE)
  }))
  rownames(ctr) <- ctr$city

  # One edge per (hub, partner-zone, direction). Outflow = colSums over hub rows; inflow =
  # rowSums over hub columns. Partner = every DRC zone EXCEPT the hub's own zones. Store the
  # hub end (hx,hy) and partner end (px,py) separately so each layer can orient its arrow to
  # the true destination: outflow hub->partner, inflow partner->hub.
  edge_dir <- function(nm, direction) {
    hz  <- URBAN_HUBS[[nm]]
    vec <- if (direction == "Outflow") colSums(M[hz, , drop = FALSE])
           else                        rowSums(M[, hz, drop = FALSE])
    vec <- vec[!(names(vec) %in% hz)]
    vec <- vec[is.finite(vec) & vec >= FLOW_MIN]
    if (!length(vec)) return(NULL)
    xy <- zc[match(tolower(trimws(names(vec))), zc$key), c("lon", "lat")]
    ok <- !is.na(xy$lon)
    if (!any(ok)) return(NULL)
    data.frame(city = nm, direction = direction, trips = as.numeric(vec[ok]),
               hx = ctr[nm, "lon"], hy = ctr[nm, "lat"],
               px = xy$lon[ok], py = xy$lat[ok], stringsAsFactors = FALSE)
  }
  ed_out <- do.call(rbind, lapply(cities, edge_dir, direction = "Outflow"))
  ed_in  <- do.call(rbind, lapply(cities, edge_dir, direction = "Inflow"))
  stopifnot(!is.null(ed_out), !is.null(ed_in), nrow(ed_out) > 0, nrow(ed_in) > 0)
  # Draw the largest flows last (on top); each direction keyed to its legend label.
  ed_out <- ed_out[order(ed_out$trips), ]; ed_out$direction <- factor("Outflow", names(FLOW_COL))
  ed_in  <- ed_in[order(ed_in$trips),  ]; ed_in$direction  <- factor("Inflow",  names(FLOW_COL))

  ctr$city <- factor(ctr$city, levels = cities)
  lab <- ctr
  # A city absent from HUB_LABEL_NUDGE yields a NULL element, which rbind() DROPS — so `nud`
  # came back shorter than `lab` and the offsets were recycled onto the wrong hubs, silently,
  # in a manuscript figure. Hub names come from urban_hub_selection.csv and are not guaranteed
  # to be keys here. Default any unknown city to no nudge and say which ones.
  .nud_key <- as.character(lab$city)
  .nud_missing <- setdiff(unique(.nud_key), names(HUB_LABEL_NUDGE))
  if (length(.nud_missing))
    message("[dominoes] no label nudge for ", paste(.nud_missing, collapse = ", "),
            "; drawn unoffset.")
  nud <- t(vapply(.nud_key, function(k) {
    v <- HUB_LABEL_NUDGE[[k]]
    if (is.null(v) || length(v) < 2L) c(0, 0) else as.numeric(v[1:2])
  }, numeric(2)))
  stopifnot(nrow(nud) == nrow(lab))
  lab$lonx <- lab$lon + nud[, 1]; lab$laty <- lab$lat + nud[, 2]

  country <- sf::st_union(sf::st_geometry(shp))          # crisp national outline
  arw <- grid::arrow(length = unit(3.2, "pt"), type = "closed")
  # Both layers use the SAME positive curvature but opposite endpoint order, so a hub<->
  # partner pair bows to opposite sides AND each arrowhead lands on its true destination.
  # sqrt width scale: intra-metro flows dwarf cross-country ones (~460x), so a linear map
  # would erase the long-range reach; sqrt keeps small edges visible without flooding.
  ggplot() +
    geom_sf(data = shp, fill = "grey91", colour = "grey80", linewidth = 0.12) +
    geom_sf(data = country, fill = NA, colour = "grey45", linewidth = 0.4) +
    geom_curve(data = ed_in,  aes(x = px, y = py, xend = hx, yend = hy,
               colour = direction, linewidth = trips),
               curvature = 0.22, alpha = 0.6, lineend = "round", arrow = arw) +
    geom_curve(data = ed_out, aes(x = hx, y = hy, xend = px, yend = py,
               colour = direction, linewidth = trips),
               curvature = 0.22, alpha = 0.6, lineend = "round", arrow = arw) +
    geom_point(data = ctr, aes(lon, lat), shape = 21, size = 2.6,
               fill = "grey20", colour = "white", stroke = 0.5) +
    geom_label(data = lab, aes(lonx, laty, label = city_label(city)),
               size = 2.9, fontface = "bold", colour = INK,
               fill = "white", alpha = 0.75, linewidth = 0, label.padding = unit(1, "pt")) +
    scale_colour_manual(values = FLOW_COL, name = NULL, limits = names(FLOW_COL),
                        guide = guide_legend(order = 1,
                          override.aes = list(linewidth = 1.6, alpha = 1))) +
    # Min width raised (Flowminder floors flows at 15, so the smallest edges cluster near
    # the low end); this keeps those faint 15-20-relocation links visible rather than hairline.
    scale_linewidth(range = c(0.3, 2.3), transform = "sqrt", name = "Relocations",
                    breaks = c(100, 1000, 4000),
                    guide = guide_legend(order = 2, override.aes = list(colour = "grey45"))) +
    coord_sf(expand = FALSE) +
    theme_map(11)
}

# =============================================================================
# Panel B — WHICH domino: cities ranked by expected additional zones
# =============================================================================
build_cityrank <- function() {
  d <- sumw |> mutate(city = factor(city, levels = rev(CITY_ORDER)))
  ggplot(d, aes(expected_added, city)) +
    # Single NEUTRAL fill (magnitude is carried by bar length alone). Deliberately non-semantic:
    # a sequential ramp here would read as a third scale and clash with Panel A's viridis
    # probability ramp and Panel C's province hues; grey keeps the bars unambiguous. Grey also
    # avoids the province palette (the old #0072B2 = Ituri), so no false province cue.
    geom_col(fill = "grey45", width = 0.66) +
    # The paired 90% Monte-Carlo interval. Shown because several of these bars are shorter
    # than the noise the unpaired version of this analysis carried, and a bar whose interval
    # touches zero should not be read as an ordering.
    {if ("expected_added_lo" %in% names(d))
       geom_linerange(aes(xmin = expected_added_lo, xmax = expected_added_hi),
                      colour = INK, linewidth = 0.45)} +
    geom_text(aes(label = sprintf("+%d newly >20%%", newly_above_20)),
              hjust = -0.08, size = 2.7, colour = MUTED) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.34))) +
    scale_y_discrete(labels = city_label) +
    labs(x = sprintf("Expected additional zones invaded within %d weeks", HORIZON), y = NULL) +
    theme_pub(11) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 10, colour = INK))
}

# =============================================================================
# Panel C — WHERE to look next: per-city watch-lists. For each zone, a PAIR of
# horizontal error bars (posterior mean + 90% CrI): baseline P(invasion) in grey
# (upper) and, below it, the conditional-on-city-invaded probability in the zone's
# province colour. Mirrors the two-horizon forest of the old Figure 3 panel C
# (make_publication_figures.R), but the pair here is baseline vs conditional.
# =============================================================================
build_watchlists <- function() {
  # rank each city's materially-elevated zones by ATTRIBUTABLE increase (its own
  # catchment); each facet is sorted independently by delta.
  d <- imp |>
    group_by(city) |>
    group_modify(~ city_watchlist(.x)) |>
    ungroup() |>
    mutate(city = factor(city, levels = CITY_ORDER),
           province = factor(province, levels = names(PCOL)),
           row = paste(city, health_zone, sep = "___"))
  ord <- d |> arrange(city, delta) |> pull(row) |> unique()
  d$row <- factor(d$row, levels = ord)
  # Draw the 90% credible-interval bars when the reach CI columns are present (emitted
  # by 38_urban_scenarios.R); degrade to points-only on older CSVs. NUDGE offsets the
  # baseline (up) and conditional (down) bars within each zone's row band.
  has_ci <- all(c("p_base_lo", "p_base_hi", "p_cond_lo", "p_cond_hi") %in% names(d)) &&
            any(is.finite(d$p_cond_lo))
  NUDGE <- 0.22
  up <- position_nudge(y = NUDGE); dn <- position_nudge(y = -NUDGE)
  xlab <- paste0(sprintf("P(invasion by %d wk): baseline (grey, upper) vs conditional on city invaded (colour, lower)", HORIZON),
                 if (has_ci) "; point = mean, bar = 90% CrI" else "")

  p <- ggplot(d)
  if (has_ci)
    p <- p +
      geom_linerange(aes(y = row, xmin = p_base_lo, xmax = p_base_hi),
                     position = up, orientation = "y", colour = "grey62", linewidth = 0.7, na.rm = TRUE) +
      geom_linerange(aes(y = row, xmin = p_cond_lo, xmax = p_cond_hi, colour = province),
                     position = dn, orientation = "y", linewidth = 0.7, na.rm = TRUE)
  p +
    geom_point(aes(x = p_base, y = row), position = up, colour = "grey45", size = 1.7) +
    geom_point(aes(x = p_cond, y = row, colour = province), position = dn, size = 2.1) +
    facet_wrap(~city, scales = "free_y", nrow = 1,
               labeller = ggplot2::as_labeller(city_label)) +
    scale_y_discrete(labels = function(x) zone_label(sub("^.*___", "", x))) +
    scale_colour_manual(values = PCOL, name = "Province", drop = TRUE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = xlab, y = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3), nrow = 1)) +
    theme_pub(11) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 9, colour = INK))
}

# =============================================================================
# Assemble + write
# =============================================================================
# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. run_all.R / run_cascade.R launch it as its own subprocess. Guarded so a source()
# — a test, an ad-hoc check, a sweep that loads "every module" — defines what it needs and
# stops, instead of silently rewriting files. DOMINOES_FORCE_RUN=1 forces a run from a source().
if (!is_script_run("40_cascade_next_dominoes.R") && !identical(Sys.getenv("DOMINOES_FORCE_RUN"), "1")) {
  message("[dominoes] sourced, not run: nothing is written. Use `Rscript 40_cascade_next_dominoes.R` to run it.")
} else {

message("== Figure: the next dominoes ==")
# A failed panel must NOT be replaced by plot_spacer(). manuscript_figures/Figure4 is a
# published deliverable and plot_annotation(tag_levels = "A") tags panels positionally, so a
# spacer ships as an EMPTY box bearing a real panel letter — a blank "B" in the manuscript,
# with only a message in the log to say why. Record each failure instead and refuse to write
# the figure unless all three panels built.
.panel_fail <- character(0)
.panel <- function(expr, nm) tryCatch(expr, error = function(e) {
  message("!! ", nm, ": ", conditionMessage(e))
  .panel_fail <<- c(.panel_fail, sprintf("%s (%s)", nm, conditionMessage(e)))
  patchwork::plot_spacer()
})
pA <- .panel(build_hubflowmap(), "flowmap")
pB <- .panel(build_cityrank(),   "cityrank")
pC <- .panel(build_watchlists(), "watchlist")

top <- (patchwork::wrap_elements(full = pA) | patchwork::wrap_elements(full = pB)) +
  patchwork::plot_layout(widths = c(0.52, 0.48))
# No embedded caption — the figure carries only the A/B/C tags. The caption is supplied
# separately in the manuscript (see the accompanying text).
fig <- (top / patchwork::wrap_elements(full = pC)) +
  patchwork::plot_layout(heights = c(0.46, 0.54)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.tag = element_text(size = 13, face = "bold", colour = INK)))

if (length(.panel_fail)) {
  # Remove any previous Figure4 so a stale copy cannot stand in for this run's missing one.
  for (ext in c("pdf", "png")) {
    .old <- file.path(MANU_DIR, paste0("Figure4.", ext))
    if (file.exists(.old)) file.remove(.old)
  }
  message(sprintf("!! [dominoes] %d of 3 panels failed; manuscript Figure4 NOT written: %s",
                  length(.panel_fail), paste(.panel_fail, collapse = "; ")))
} else {
  save_dual(fig, "Figure4", 13.5, 9.6, MANU_DIR)
}

# ---- watch-list CSV underlying Panel C -------------------------------------
watch <- imp |>
  group_by(city) |>
  group_modify(~ city_watchlist(.x)) |>
  arrange(city, desc(delta)) |>
  mutate(rank_in_city = row_number()) |>
  ungroup()
.has_ci_csv <- all(c("p_base_lo", "p_base_hi", "p_cond_lo", "p_cond_hi") %in% names(watch))
watch <- watch |>
  transmute(city = factor(city, levels = CITY_ORDER), rank_in_city, health_zone, province,
            # posterior means + 90% credible intervals mirror the Panel C paired error bars
            p_baseline = round(p_base, 4),
            p_baseline_lo = if (.has_ci_csv) round(p_base_lo, 4) else NA_real_,
            p_baseline_hi = if (.has_ci_csv) round(p_base_hi, 4) else NA_real_,
            p_conditional = round(p_cond, 4),
            p_conditional_lo = if (.has_ci_csv) round(p_cond_lo, 4) else NA_real_,
            p_conditional_hi = if (.has_ci_csv) round(p_cond_hi, 4) else NA_real_,
            attributable_increase = round(delta, 4), relative_risk = round(rr, 3),
            horizon = HORIZON, seed_week) |>
  arrange(city, rank_in_city)
if (!dir.exists(KEY_DIR)) dir.create(KEY_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(watch, file.path(KEY_DIR, "Figure_next_dominoes_watchlist.csv"))
# The watch-list CSV is an independent deliverable and is written either way; only now, with
# every other product safely on disk, does a panel failure become a hard error so the caller
# (run_cascade.R's tryCatch, or run_one()'s exit status) records it instead of the run
# reporting success with a figure missing.
if (length(.panel_fail))
  stop(sprintf("[dominoes] %d of 3 panels failed, so manuscript Figure4 was not written (a spacer would have shipped as a blank tagged panel): %s",
               length(.panel_fail), paste(.panel_fail, collapse = "; ")), call. = FALSE)
message("[done] Figure4 (next dominoes) -> ", MANU_DIR)

}   # end of the script-run guard
