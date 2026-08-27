# =============================================================================
# 40_cascade_next_dominoes.R — CASCADE: the "next dominoes" narrative figure
# BDBV 2026 DRC · 3-month spatial invasion cascade
#
# ONE coherent figure answering "if a major city is invaded, where do we look
# next?" — combining the eastern-epicentre frontier (Figure4_cascade anchor) with
# the urban conditional scenarios (Figure_urban_*), across all five hub cities.
#
#   Panel A  ANCHOR — national 13-week reach map from the CURRENT eastern
#            epicentre (baseline cascade, no city seeded). "Where it is going now."
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
# Standalone & regenerable: reads only saved CSVs (no re-simulation), mirroring
# 37_cascade_figure4.R. Works in both spatiotemporal/ and spatiotemporal_conditional/.
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
AFFECTED_FILL <- "grey78"; NA_FILL <- "grey93"
OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
# Province identity spanning the eastern epicentre AND the western/southern urban
# catchments (Okabe-Ito, CVD-safe). Merged from 37 (east) + 38 (urban).
PROV_COL <- c(
  "Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73","Sud-Kivu"="#CC79A7",
  "Kinshasa"="#E69F00","Haut-Katanga"="#56B4E9","Kasaï-Oriental"="#CC79A7",
  "Kasaï-Central"="#009E73","Kasaï"="#D55E00","Kwilu"="#0072B2",
  "Kwango"="#F0E442","Lomami"="#000000","Mai-Ndombe"="#7A7A7A")
PROV_INT <- get0("PROVINCES_OF_INTEREST", ifnotfound = c("Ituri","Nord-Kivu","Haut-Uele"))
base_family <- "sans"
HORIZON <- max(get0("CASCADE_REPORT_HORIZONS", ifnotfound = c(4L, 8L, 13L)))
TOP_N   <- 8L                                    # max zones per city watch-list
# Below this per-zone attributable increase, the baseline->conditional difference is
# Monte-Carlo flicker, not signal (mirrors URBAN_DELTA_MIN in 38_urban_scenarios.R).
# Watch-lists show only MATERIALLY-elevated zones so noise is never shown as a target.
DELTA_MIN <- 0.03

# top materially-elevated zones for one city, ranked by attributable increase (delta).
# Shared by the figure and the CSV so both stay identical.
city_watchlist <- function(df, top_n = TOP_N)
  df |> dplyr::filter(is.finite(delta), delta > DELTA_MIN) |>
    dplyr::slice_max(order_by = delta, n = top_n, with_ties = FALSE)

# assign a stable colour to every province (identity where known, spare Okabe else)
prov_palette <- function(provs) {
  provs <- unique(as.character(provs[!is.na(provs)]))
  provs <- c(intersect(PROV_INT, provs), sort(setdiff(provs, PROV_INT)))
  spare <- setdiff(OKABE, unname(PROV_COL)); si <- 1L
  cols <- setNames(character(length(provs)), provs)
  for (p in provs) {
    if (p %in% names(PROV_COL)) cols[p] <- unname(PROV_COL[p])
    else { cols[p] <- spare[((si - 1L) %% max(length(spare), 1L)) + 1L]; si <- si + 1L }
  }
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
    legend.position = "right",
    legend.title = element_text(size = base - 1.4, colour = MUTED),
    legend.text = element_text(size = base - 1.8, colour = INK),
    legend.key.height = unit(20, "pt"), legend.key.width = unit(7, "pt"),
    plot.tag = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin = margin(2, 2, 2, 2))
}
save_dual <- function(p, name, w, h, dir) {
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
sumw   <- summ |> filter(seed_week == min_wk) |> arrange(desc(expected_added))
CITY_ORDER <- sumw$city

# shared province palette across every panel
all_provs <- unique(c(imp$province, rs$province[rs$horizon == HORIZON]))
PCOL <- prov_palette(all_provs)

# ---- shapefile join (as in 37_cascade_figure4.R) ---------------------------
shp <- st_read(SHAPEFILE_PATH, quiet = TRUE) |>
  mutate(.key = tolower(trimws(Nom)))
join_map <- function(dat) {
  d <- dat |> mutate(.key = tolower(trimws(health_zone)))
  shp |> left_join(d, by = ".key")
}

# =============================================================================
# Panel A — ANCHOR: current 13-week frontier from the eastern epicentre
# =============================================================================
build_anchor <- function() {
  d <- rs |> filter(horizon == HORIZON) |>
    mutate(active = was_active_before %in% TRUE,
           fill_p = ifelse(active, NA_real_, p_case_invasion))
  mp  <- join_map(d |> dplyr::select(health_zone, fill_p, active))
  ggplot(mp) +
    geom_sf(aes(fill = fill_p), colour = "white", linewidth = 0.07) +
    geom_sf(data = mp |> filter(active %in% TRUE),
            fill = AFFECTED_FILL, colour = "white", linewidth = 0.07) +
    # SAME legend palette as manuscript_figures/Figure3.pdf panel A (continuous viridis,
    # direction = 1: dark = low, yellow = high). Kept linear/percent since this is a
    # probability P(invaded by HORIZON wk), not the log-scaled relative risk of Figure 3.
    scale_fill_viridis_c(option = "viridis", direction = 1, na.value = NA_FILL,
                         limits = c(0, NA), labels = percent_format(1),
                         name = sprintf("P(invaded\nby %d wk)", HORIZON)) +
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
    geom_text(aes(label = sprintf("+%d newly >20%%", newly_above_20)),
              hjust = -0.08, size = 2.7, colour = MUTED) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.34))) +
    labs(x = sprintf("Expected additional zones invaded within %d weeks", HORIZON), y = NULL) +
    theme_pub(11) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 10, colour = INK))
}

# =============================================================================
# Panel C — WHERE to look next: per-city watch-lists (dumbbell = both metrics)
# =============================================================================
build_watchlists <- function() {
  # rank each city's materially-elevated zones by ATTRIBUTABLE increase (its own
  # catchment); reorder_within so each facet is sorted independently by delta.
  d <- imp |>
    group_by(city) |>
    group_modify(~ city_watchlist(.x)) |>
    ungroup() |>
    mutate(city = factor(city, levels = CITY_ORDER),
           province = factor(province, levels = names(PCOL)),
           row = paste(city, health_zone, sep = "___"))
  ord <- d |> arrange(city, delta) |> pull(row) |> unique()
  d$row <- factor(d$row, levels = ord)

  ggplot(d) +
    geom_segment(aes(x = p_base, xend = p_cond, y = row, yend = row),
                 colour = "grey75", linewidth = 1) +
    geom_point(aes(p_base, row), colour = "grey55", size = 1.6) +
    geom_point(aes(p_cond, row, colour = province), size = 2.5) +
    facet_wrap(~city, scales = "free_y", nrow = 1) +
    scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_colour_manual(values = PCOL, name = "Province", drop = TRUE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = sprintf("P(invasion by %d wk): baseline (grey) to conditional on city invaded (colour)", HORIZON),
         y = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3), nrow = 1)) +
    theme_pub(11) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 9, colour = INK))
}

# =============================================================================
# Assemble + write
# =============================================================================
message("== Figure: the next dominoes ==")
pA <- tryCatch(build_anchor(),     error = function(e) { message("!! anchor: ",   conditionMessage(e)); patchwork::plot_spacer() })
pB <- tryCatch(build_cityrank(),   error = function(e) { message("!! cityrank: ", conditionMessage(e)); patchwork::plot_spacer() })
pC <- tryCatch(build_watchlists(), error = function(e) { message("!! watchlist: ",conditionMessage(e)); patchwork::plot_spacer() })

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

save_dual(fig, "Figure4", 13.5, 9.6, MANU_DIR)

# ---- watch-list CSV underlying Panel C -------------------------------------
watch <- imp |>
  group_by(city) |>
  group_modify(~ city_watchlist(.x)) |>
  arrange(city, desc(delta)) |>
  mutate(rank_in_city = row_number()) |>
  ungroup() |>
  transmute(city = factor(city, levels = CITY_ORDER), rank_in_city, health_zone, province,
            p_baseline = round(p_base, 4), p_conditional = round(p_cond, 4),
            attributable_increase = round(delta, 4), relative_risk = round(rr, 3),
            horizon = HORIZON, seed_week) |>
  arrange(city, rank_in_city)
if (!dir.exists(KEY_DIR)) dir.create(KEY_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(watch, file.path(KEY_DIR, "Figure_next_dominoes_watchlist.csv"))
message("[done] Figure4 (next dominoes) -> ", MANU_DIR)
