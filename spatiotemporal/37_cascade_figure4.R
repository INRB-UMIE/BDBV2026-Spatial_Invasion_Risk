# =============================================================================
# 37_cascade_figure4.R — CASCADE version of publication Figure 4
# BDBV 2026 DRC · 3-month spatial invasion cascade
#
# A faithful cascade analogue of make_publication_figures.R's Figure 4:
#   Panel A  log-scaled RELATIVE invasion-risk map (rr01_nat, 0-1) for horizon H
#   Panel B  top-20 at-risk zones by H-week reach probability with 90% CrI,
#            every zone coloured by its province (labels coloured to match)
# Combined A|B at 60:40, tagged A/B, no embedded titles (house convention).
# Reads THIS folder's saved cascade reach scores (primary = status-quo scenario)
# so it works standalone in both spatiotemporal/ and spatiotemporal_conditional/.
#
# Run:  Rscript spatiotemporal/37_cascade_figure4.R
# Out:  outputs/key_outputs/figures/Figure4_cascade_h{4,8,13}.{pdf,png}
#       outputs/key_outputs/figures/panels/F4A_cascade_reach_map_h*, F4B_cascade_rank_h*
#       outputs/key_outputs/Figure4_cascade_h{4,8,13}_data.csv
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(patchwork)
  library(sf); library(scales); library(here)
})
ST_DIR <- Sys.getenv("CASCADE_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
source(file.path(ST_DIR, "00_config.R"))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- style constants (copied verbatim from make_publication_figures.R) ------
INK <- "grey15"; MUTED <- "grey38"; GRID <- "grey92"
AFFECTED_FILL <- "grey78"; NA_FILL <- "grey93"; RR_FLOOR <- 1e-3
OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
PROV_COL <- c("Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73",
              "Sud-Kivu"="#CC79A7","Other"="#7A7A7A")
PROV_INT <- get0("PROVINCES_OF_INTEREST", ifnotfound = c("Ituri","Nord-Kivu","Haut-Uele"))
base_family <- "sans"

theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    axis.title = element_text(size = base - 0.4, colour = MUTED),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4), angle = 90),
    axis.text = element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
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
  message(sprintf("  saved %-34s  %.1f x %.1f in", name, w, h)); invisible(p)
}

# ---- load shapefile + cascade reach scores ---------------------------------
shp <- st_read(SHAPEFILE_PATH, quiet = TRUE) %>%
  mutate(.key = tolower(trimws(Nom)), .prov = as.character(PROVINCE))
join_map <- function(dat) {
  d <- dat %>% mutate(.key = tolower(trimws(health_zone)), .prov = as.character(province))
  m <- shp %>% left_join(d, by = c(".key", ".prov"))
  miss <- m %>% st_drop_geometry() %>% summarise(f = mean(is.na(health_zone))) %>% pull(f)
  if (miss > 0.5) { d2 <- dat %>% mutate(.key = tolower(trimws(health_zone)))
    m <- shp %>% left_join(d2 %>% dplyr::select(-any_of(".prov")), by = ".key") }
  m
}
FIG_DIR <- file.path(OUT_DIR, "key_outputs", "figures")
PANEL_DIR <- file.path(FIG_DIR, "panels")
KEY_DIR <- file.path(OUT_DIR, "key_outputs")
reach_csv <- file.path(OUT_DIR, "cascade", "tables", "cascade_reach_scores_all_zones.csv")
if (!file.exists(reach_csv)) stop("cascade reach scores not found: ", reach_csv,
                                  " — run run_cascade.R first.")
rs <- read_csv(reach_csv, show_col_types = FALSE)
stopifnot(all(c("health_zone","province","horizon","was_active_before",
                "p_case_invasion","p_case_lo","p_case_hi","rr01_nat") %in% names(rs)))

# ---- build_fig4_cascade(H) --------------------------------------------------
build_fig4_cascade <- function(H, top_n = 20L) {
  wk_label <- sprintf("%d-week", H)
  # --- Panel A: log relative-invasion-risk map for horizon H ---
  rr_df <- rs %>% filter(horizon == H) %>%
    mutate(active = was_active_before %in% TRUE,
           rr_log = ifelse(active, NA, pmax(rr01_nat, RR_FLOOR)))
  mp <- join_map(rr_df)
  p_map <- ggplot(mp) +
    geom_sf(aes(fill = rr_log), colour = "white", linewidth = 0.08) +
    geom_sf(data = mp %>% filter(was_active_before %in% TRUE),
            fill = AFFECTED_FILL, colour = "white", linewidth = 0.08) +
    scale_fill_viridis_c(option = "viridis", trans = "log10", limits = c(RR_FLOOR, 1),
                         breaks = c(0.001, 0.01, 0.1, 1), labels = c("0.001","0.01","0.1","1"),
                         na.value = NA_FILL,
                         name = sprintf("Relative\ninvasion risk\n(%s, log, 0-1)", wk_label),
                         direction = 1) +
    theme_map(11)

  # --- Panel B: top-N zones by H-week reach + 90% CrI, coloured by province ---
  d <- rs %>% filter(horizon == H, !(was_active_before %in% TRUE), is.finite(p_case_invasion)) %>%
    arrange(desc(p_case_invasion)) %>% head(top_n)
  prov_present <- unique(d$province[!is.na(d$province)])
  prov_levels  <- c(intersect(PROV_INT, prov_present), sort(setdiff(prov_present, PROV_INT)))
  spare <- setdiff(OKABE, unname(PROV_COL))
  prov_cols <- setNames(character(length(prov_levels)), prov_levels); si <- 1L
  for (pv in prov_levels) {
    if (pv %in% names(PROV_COL)) prov_cols[pv] <- unname(PROV_COL[pv])
    else { prov_cols[pv] <- spare[((si - 1L) %% max(length(spare), 1L)) + 1L]; si <- si + 1L }
  }
  d <- d %>% mutate(province = factor(province, levels = prov_levels),
                    zone = factor(health_zone, levels = rev(health_zone)))
  lev <- levels(d$zone)
  ylab_cols <- unname(prov_cols[as.character(d$province)[match(lev, d$health_zone)]])
  p_forest <- ggplot(d, aes(p_case_invasion, zone, colour = province)) +
    geom_linerange(aes(xmin = p_case_lo, xmax = p_case_hi), linewidth = 1.5, alpha = 0.5) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = prov_cols, name = "Province", drop = FALSE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0.01, 0.06))) +
    labs(x = sprintf("Probability of invasion within %d weeks (cascade)", H), y = NULL) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2.4, size = 3))) +
    theme_pub(11) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(colour = ylab_cols, size = 9.5, face = "bold"),
          legend.position = "top")

  # --- data underlying both panels -> key_outputs/Figure4_cascade_h{H}_data.csv ---
  rank_tbl <- rs %>% filter(horizon == H, !(was_active_before %in% TRUE), is.finite(p_case_invasion)) %>%
    arrange(desc(p_case_invasion)) %>%
    transmute(health_zone, province, invasion_prob_rank = row_number())
  fig4_data <- rs %>% filter(horizon == H) %>%
    left_join(rank_tbl, by = c("health_zone", "province")) %>%
    mutate(in_top20_panelB = !is.na(invasion_prob_rank) & invasion_prob_rank <= top_n,
           rr_log10_fill = ifelse(was_active_before %in% TRUE, NA_real_, pmax(rr01_nat, RR_FLOOR))) %>%
    arrange(invasion_prob_rank) %>%
    transmute(invasion_prob_rank, in_top20_panelB, health_zone, province, horizon,
              was_active_before, p_invasion = p_case_invasion, p_invasion_lo = p_case_lo,
              p_invasion_hi = p_case_hi, relative_invasion_risk = rr01_nat, rr_log10_fill)
  if (!dir.exists(KEY_DIR)) dir.create(KEY_DIR, recursive = TRUE, showWarnings = FALSE)
  write_csv(fig4_data, file.path(KEY_DIR, sprintf("Figure4_cascade_h%d_data.csv", H)))

  save_dual(p_map,    sprintf("F4A_cascade_reach_map_h%d", H), 6.4, 6.4, PANEL_DIR)
  save_dual(p_forest, sprintf("F4B_cascade_rank_h%d", H),      6.6, 6.4, PANEL_DIR)
  fig4 <- (wrap_elements(full = p_map) | wrap_elements(full = p_forest)) +
    plot_layout(widths = c(0.6, 0.4)) + plot_annotation(tag_levels = "A")
  save_dual(fig4, sprintf("Figure4_cascade_h%d", H), 12.8, 6.6, FIG_DIR)
  invisible(fig4)
}

# ---- RUN for each reporting horizon (13-week = the headline 3-month figure) --
horizons <- get0("CASCADE_REPORT_HORIZONS", ifnotfound = c(4L, 8L, 13L))
for (H in horizons) {
  message(sprintf("== Figure 4 (cascade, %d-week reach + uncertainty) ==", H))
  tryCatch(build_fig4_cascade(H), error = function(e)
    message("!! Figure4_cascade_h", H, " FAILED: ", conditionMessage(e)))
}
# copy the headline (13-week) to the canonical Figure4_cascade name
for (ext in c("pdf","png")) {
  src <- file.path(FIG_DIR, sprintf("Figure4_cascade_h%d.%s", max(horizons), ext))
  if (file.exists(src)) file.copy(src, file.path(FIG_DIR, sprintf("Figure4_cascade.%s", ext)),
                                  overwrite = TRUE)
}
message("[done] Figure4_cascade -> ", FIG_DIR)
