# =============================================================================
# figure_style.R — the manuscript figure design system
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# ONE definition of the house style, sourced by every manuscript/SI figure builder. It was
# previously a local copy inside make_manuscript_figures.R, which meant any new figure suite
# had to re-declare it and could drift from the figures it sits beside in the manuscript —
# a caption saying "same style" is not a mechanism.
#
# Contents: the colourblind-safe palette, theme_pub()/theme_map() (both title-free: the main
# text carries the captions), and save_dual(), which writes a PDF + 600-dpi PNG through the
# FIGURE_KEEP gate and the probability-scale stamp.
#
# save_dual() needs a default output directory. Each suite sets PANEL_DIR before sourcing
# this file (or passes `dir =` at every call site); the default is resolved at CALL time, so
# sourcing order does not matter.
# =============================================================================
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
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  # Gate on the FULL destination path, not the bare stem: FIGURE_DROP entries are
  # "<directory>/<stem>" and the raw/ exclusion inspects path components, neither of
  # which can match a basename.
  if (is.function(.fk) && !.fk(file.path(dir, name))) return(invisible(p))
  p <- fs_caption(p)                       # scale statement (raw pass by default)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  message(sprintf("  saved %-32s  %.1f x %.1f in", name, w, h)); invisible(p)
}

