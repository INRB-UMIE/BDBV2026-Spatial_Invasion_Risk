# =============================================================================
# make_spread_kinematics_compact.R
# BDBV 2026 DRC — Spatiotemporal invasion
#
# A three-panel MANUSCRIPT version of Figure_spread_kinematics, carrying only the
# three time-series panels of the five-panel diagnostic figure (its B, D and F,
# relabelled A-C here) and dropping the map and the dispersion panel.
#
#   A  reach from the index zone       (was Panel B)
#   B  centre-of-gravity step speed    (was Panel D)
#   C  new health zones invaded        (built here from spread_kinematics_weekly.csv)
#
# Differences from the diagnostic figure, all deliberate:
#   * NO in-panel explanatory text. The diagnostic figure annotates each panel with
#     the fitted rate, its interval and a sentence of interpretation; here that
#     belongs to the caption, so every annotate("label")/note() block is removed —
#     together with the two overlays whose meaning lived ONLY in that text (the
#     week-of-90%-extent time marker and the median-reach fit line in Panel A).
#     Reference lines that a caption can name in one clause are KEPT: the two
#     centre-of-gravity medians in B and the mean invasion rate in C.
#   * NO SHADED RECTANGLES. The right-truncation bands in all three panels and the
#     bootstrap band on the invasion rate in C are gone; right truncation is carried
#     by the open-symbol / pale-bar convention alone, which says the same thing
#     without a background block sitting behind the data.
#   * House manuscript style: no titles, no subtitles, identity carried by the A-C
#     tags and the axis titles (make_manuscript_figures.R).
#   * Large type throughout (FONT_SCALE below) and a modern colourblind-safe
#     palette (Paul Tol "vibrant"), replacing Okabe-Ito.
#
# Rebuilt from SAVED data — the CSVs written by 43_spread_kinematics.R — so it
# re-fits nothing and cannot drift from the numbers that module published.
#
# Run:  Rscript spatiotemporal/make_spread_kinematics_compact.R
# Out:  outputs/key_outputs/manuscript_figures/Figure_spread_kinematics_compact.{pdf,png}
#       outputs/key_outputs/manuscript_figures/panels/SKC_{A,B,C}_*.{pdf,png}
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(scales); library(here)
})

ST_DIR    <- Sys.getenv("SPREAD_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
KEY_DIR   <- file.path(ST_DIR, "outputs", "key_outputs")
FIG_DIR   <- file.path(KEY_DIR, "manuscript_figures")
PANEL_DIR <- file.path(FIG_DIR, "panels")
for (d in c(FIG_DIR, PANEL_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Load the retained-figure allow-list (FIGURE_KEEP / figure_is_kept, 00_config.R). Without it
# the gate in save_dual() below resolves figure_is_kept to NULL and no-ops, so this script wrote
# its SKC_* working panels despite none of them being a deliverable. Non-fatal: the script must
# still run standalone from a checkout with no config on the path.
if (!is.function(get0("figure_is_kept")))
  try(suppressWarnings(suppressMessages(source(file.path(ST_DIR, "00_config.R")))), silent = TRUE)

W <- readr::read_csv(file.path(KEY_DIR, "spread_kinematics_weekly.csv"), show_col_types = FALSE)
S <- readr::read_csv(file.path(KEY_DIR, "spread_kinematics_summary.csv"), show_col_types = FALSE)
stopifnot(nrow(W) >= 4, all(c("week", "truncated", "new_zones") %in% names(W)))

#' One published statistic, by its exact row name; errors rather than returning NA so a
#' renamed row in the source module can never reach the figure as a silently missing line.
stat <- function(nm, col = "value") {
  i <- which(S$statistic == nm)
  stopifnot(length(i) == 1L)
  v <- S[[col]][i]; stopifnot(is.finite(v)); v
}
# READ, not computed. 43_spread_kinematics.R now publishes "Mean new invasions per week"
# (value + 90% week-bootstrap interval) to spread_kinematics_summary.csv, over exactly the
# window this panel draws: weeks 1..fit_max, excluding week 0 (the seeding condition, whose
# zones already had confirmed cases when the series opens) and the right-truncated tail.
# This script used to compute the mean itself from W, and its Panel C comment asserted that
# the interval "is in spread_kinematics_summary.csv" — there was no such row in any of the 29.
NZ_MEAN <- stat("Mean new invasions per week")
NZ_LO   <- stat("Mean new invasions per week", "lo")
NZ_HI   <- stat("Mean new invasions per week", "hi")
V_COG   <- stat("Case centre-of-gravity speed")
V_INV   <- stat("New-invasion centre-of-gravity speed")
FIT_MAX <- max(W$week[!W$truncated])

# -----------------------------------------------------------------------------
# Design system — modern, colourblind-safe (Paul Tol "vibrant")
# -----------------------------------------------------------------------------
# One SEMANTIC colour per quantity, held across all three panels: a reader who learns
# "orange = cases" in A must not meet orange standing for something else in B or C.
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
CASE_COL <- "#EE7733"   # anything CASE-weighted (that week's confirmed cases)
INV_COL  <- "#009988"   # anything about the zones NEWLY invaded that week
EXT_COL  <- "#0077BB"   # the cumulative invaded set — its outer envelope
EXT_COL2 <- "#33BBEE"   # lighter shade — the middle of that same cumulative set
base_family <- "sans"
FONT_SCALE <- 1.85      # one knob for every text size in the figure

theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    axis.title = element_text(size = base - 0.4, colour = MUTED),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5), angle = 90),
    axis.text = element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.35),
    legend.position = "top", legend.justification = "left",
    legend.title = element_blank(),
    legend.text = element_text(size = base - 1.4, colour = INK),
    legend.key.height = unit(11, "pt"), legend.key.width = unit(20, "pt"),
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(4, "pt"),
    plot.margin = margin(6, 10, 6, 6))
}
BASE <- 8.6 * FONT_SCALE

save_dual <- function(p, name, w, h, dir = PANEL_DIR) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  # Gate on the FULL destination path, not the bare stem: FIGURE_DROP entries are
  # "<directory>/<stem>" and the raw/ exclusion inspects path components, neither of
  # which can match a basename.
  if (is.function(.fk) && !.fk(file.path(dir, name))) return(invisible(p))
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  message(sprintf("  saved %-38s %.1f x %.1f in", name, w, h)); invisible(p)
}

# Shared scaffolding: the week axis and the hollow/solid point convention, identical in
# all three panels so one caption clause covers them. FIT_MAX is used only to decide
# which weeks get the truncated treatment — there is no shaded band anywhere here.
WK_SCALE <- scale_x_continuous(breaks = scales::breaks_width(2),
                               expand = expansion(mult = c(0.02, 0.03)))
SHAPE_TRUNC <- scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none")

# -----------------------------------------------------------------------------
# Panel A — reach from the index zone
# -----------------------------------------------------------------------------
build_reach <- function() {
  lv <- c("Furthest zone so far", "Furthest new zone that week",
          "Median zone so far", "Weekly cases, mean")
  RCOL <- setNames(c(EXT_COL, INV_COL, EXT_COL2, CASE_COL), lv)
  d <- W %>%
    dplyr::transmute(week, truncated,
                     `Furthest zone so far`        = reach_max_km,
                     `Furthest new zone that week` = new_max_km,
                     `Median zone so far`          = reach_med_km,
                     `Weekly cases, mean`          = case_mean_km) %>%
    tidyr::pivot_longer(-c(week, truncated), names_to = "series", values_to = "km") %>%
    dplyr::filter(is.finite(km)) %>%
    dplyr::mutate(series = factor(series, levels = lv))
  ggplot(d, aes(week, km, colour = series)) +
    geom_line(aes(group = series), linewidth = 0.9) +
    geom_point(aes(shape = truncated), size = 2.6, stroke = 0.9, fill = "white") +
    SHAPE_TRUNC +
    scale_colour_manual(values = RCOL, guide = guide_legend(nrow = 2, byrow = TRUE)) +
    WK_SCALE + expand_limits(y = 0) +
    labs(x = "Epidemic week", y = "Road distance from index zone (km)") +
    theme_pub(BASE)
}

# -----------------------------------------------------------------------------
# Panel B — centre-of-gravity step speed
# -----------------------------------------------------------------------------
build_cogspeed <- function() {
  lv <- c("Weekly cases", "Newly invaded zones")
  CCOL <- setNames(c(CASE_COL, INV_COL), lv)
  d <- W %>%
    dplyr::transmute(week, truncated,
                     `Weekly cases` = cog_speed_km_wk,
                     `Newly invaded zones` = inv_cog_speed_km_wk) %>%
    tidyr::pivot_longer(-c(week, truncated), names_to = "series", values_to = "kmwk") %>%
    dplyr::filter(is.finite(kmwk)) %>%
    dplyr::mutate(series = factor(series, levels = lv))
  # Log y: the invasion centre moves one to two orders of magnitude faster than the case
  # centre, so a linear axis flattens the case series onto the baseline. coord_cartesian
  # (which clips) rather than scale limits (which drop data) pins the window.
  rng <- range(d$kmwk); ylo <- rng[1] / 1.2; yhi <- rng[2] * 1.2
  ggplot(d, aes(week, kmwk, colour = series)) +
    geom_hline(yintercept = V_COG, colour = CASE_COL, linewidth = 0.7, linetype = "22") +
    geom_hline(yintercept = V_INV, colour = INV_COL,  linewidth = 0.7, linetype = "22") +
    geom_line(aes(group = series), linewidth = 0.9) +
    geom_point(aes(shape = truncated), size = 2.6, stroke = 0.9, fill = "white") +
    SHAPE_TRUNC +
    scale_colour_manual(values = CCOL, guide = guide_legend(nrow = 2, byrow = TRUE)) +
    WK_SCALE +
    scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300),
                  labels = c("1", "3", "10", "30", "100", "300")) +
    annotation_logticks(sides = "l", colour = FAINT, linewidth = 0.25,
                        short = unit(2, "pt"), mid = unit(3.2, "pt"), long = unit(4.6, "pt")) +
    coord_cartesian(ylim = c(ylo, yhi)) +
    labs(x = "Epidemic week", y = "Centre-of-gravity speed (km per week, log scale)") +
    theme_pub(BASE)
}

# -----------------------------------------------------------------------------
# Panel C — new health zones invaded per week
# -----------------------------------------------------------------------------
build_rate <- function() {
  d <- W %>% dplyr::transmute(week, truncated, new_zones)
  # The mean is drawn ONLY across the weeks it summarises: weeks 1 to FIT_MAX. Week 0 is
  # excluded because its zones already had confirmed cases when the series opens (the
  # seeding condition, not one week's spread) and the trailing weeks because they are
  # right-truncated. Running the line to the panel edge instead would assert a reference
  # level over weeks that never entered it. Its 90% bootstrap interval was drawn as a shaded
  # band and has been removed with the other rectangles; the interval is published as the
  # lo/hi of the same summary row (NZ_LO/NZ_HI here) and belongs in the caption or the text.
  ggplot(d, aes(week, new_zones)) +
    annotate("segment", x = 0.5, xend = FIT_MAX + 0.5, y = NZ_MEAN, yend = NZ_MEAN,
             colour = INV_COL, linewidth = 0.8, linetype = "22") +
    geom_col(aes(alpha = truncated), fill = INV_COL, width = 0.7) +
    WK_SCALE +
    scale_y_continuous(breaks = scales::breaks_width(3),
                       expand = expansion(mult = c(0, 0.06))) +
    # The dashed reference level now states its own uncertainty on the panel, from the
    # published 90% week-bootstrap interval, rather than leaving a bare line to be read as
    # exact. Both numbers come from spread_kinematics_summary.csv.
    labs(x = "Epidemic week", y = "New health zones invaded",
         subtitle = sprintf("mean %.1f new zones/week over weeks 1-%d [90%% CI %.1f-%.1f]",
                            NZ_MEAN, FIT_MAX, NZ_LO, NZ_HI)) +
    theme_pub(BASE) +
    # This panel has nothing to key, but the other two carry a two-row legend above the
    # plot area. Dropping the legend outright would let this panel grow taller than its
    # neighbours and break the shared baseline, so the space is reserved with a legend
    # whose keys and labels are both invisible — NOT with override.aes alone, which
    # blanks the key but still prints the factor levels ("FALSE"/"TRUE") as labels.
    guides(alpha = guide_legend(nrow = 2, byrow = TRUE, title = NULL,
                                override.aes = list(alpha = 0, fill = NA))) +
    scale_alpha_manual(values = c(`FALSE` = 0.88, `TRUE` = 0.32),
                       labels = c(" ", " "), name = NULL) +
    theme(legend.position = "top", legend.key = element_blank(),
          legend.text = element_text(colour = NA),
          # theme_pub() BLANKS plot.subtitle, so the reference level's interval would have been
          # computed and then thrown away — the same defect as a caption that never renders.
          # Re-enable it for this panel only.
          plot.subtitle = element_text(size = BASE - 2.5, colour = MUTED, hjust = 0))
}

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
pA <- build_reach(); pB <- build_cogspeed(); pC <- build_rate()
save_dual(pA, "SKC_A_reach_from_index_zone", 6.4, 5.4)
save_dual(pB, "SKC_B_cog_speed",             6.4, 5.4)
save_dual(pC, "SKC_C_new_zones_per_week",    6.4, 5.4)

fig <- (pA | pB | pC) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = BASE + 5, face = "bold", colour = INK))
save_dual(fig, "Figure_spread_kinematics_compact", 19.2, 5.6, dir = FIG_DIR)
message("[compact] done — ", FIG_DIR)
