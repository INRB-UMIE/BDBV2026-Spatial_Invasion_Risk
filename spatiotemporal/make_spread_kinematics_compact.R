# =============================================================================
# make_spread_kinematics_compact.R
# BDBV 2026 DRC — Spatiotemporal invasion
#
# Figure_spread_kinematics_compact: three time-series panels of the outbreak's spatial
# spread, drawn from the tables 43_spread_kinematics.R publishes (43 itself draws no
# figure; the five-panel diagnostic figure it once built has been removed).
#
#   A  road distance from the index zone: the cumulative invaded set's furthest,
#      90th-percentile (the envelope the saturation verdict is computed on) and median
#      zone, the furthest zone newly invaded that week, and the case-weighted mean
#   B  centre-of-gravity step speed (weekly cases; newly invaded zones), log scale
#   C  new health zones invaded per week, with the published mean over weeks 1..fit_max
#
# Design, all deliberate:
#   * NO in-panel explanatory text, with ONE exception: Panel C's subtitle states the mean
#     invasion rate and its published 90% week-bootstrap interval, so the dashed reference
#     level is never read as exact. Everything else (fitted rates, the saturation verdict)
#     belongs to the caption. Reference lines a caption can name in one clause are drawn:
#     the two centre-of-gravity medians in B and the mean invasion rate in C.
#   * NO SHADED RECTANGLES. Right truncation is carried by the open-symbol / pale-bar
#     convention alone (weeks flagged `truncated` in the weekly CSV).
#   * A line BREAKS at a week where its series is undefined (no new invasion that week),
#     rather than being drawn straight through it.
#   * House manuscript style: no titles, identity carried by the A-C tags and the axis
#     titles; legends untitled (the axis title names the quantity).
#   * Large type throughout (FONT_SCALE below) and a colourblind-safe palette (Paul Tol
#     "vibrant").
#
# Rebuilt from SAVED data — the CSVs written by 43_spread_kinematics.R — so it re-fits
# nothing and cannot drift from the numbers that module published.
#
# Run:  Rscript spatiotemporal/make_spread_kinematics_compact.R
# In:   <OUT_DIR>/key_outputs/spread_kinematics_{weekly,summary}.csv
# Out:  <OUT_DIR>/key_outputs/manuscript_figures/Figure_spread_kinematics_compact.{pdf,png}
#       (OUT_DIR from 00_config.R: spatiotemporal/outputs/, or outputs_sensitivity/<tag>/
#       under a data-definition sensitivity, i.e. wherever 43 wrote its tables. The single-
#       panel SKC_* saves below are not on FIGURE_KEEP, so the save_dual() gate skips them.)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(scales); library(here)
})

ST_DIR <- Sys.getenv("SPREAD_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
# The config is REQUIRED, not optional: it defines OUT_DIR (where 43 wrote the tables this
# figure reads), the run guard script_run_or_forced() below, and the retained-figure allow-list
# (FIGURE_KEEP / figure_is_kept) — without which save_dual()'s gate no-ops and the script wrote
# its SKC_* working panels despite none of them being a deliverable. A session that has already
# sourced it (a test, run_all's own session) is not made to re-source it.
if (!exists("OUT_DIR") || !is.function(get0("script_run_or_forced")) ||
    !is.function(get0("figure_is_kept")))
  source(file.path(ST_DIR, "00_config.R"))
# READ FROM, AND WRITE TO, THE RUN'S OWN TREE. This was hard-coded to spatiotemporal/outputs/,
# so under a data-definition sensitivity (43 then writes to outputs_sensitivity/<tag>/) the
# figure would have drawn the PUBLISHED run's tables.
KEY_DIR   <- file.path(OUT_DIR, "key_outputs")
FIG_DIR   <- file.path(KEY_DIR, "manuscript_figures")
PANEL_DIR <- file.path(FIG_DIR, "panels")

W <- readr::read_csv(file.path(KEY_DIR, "spread_kinematics_weekly.csv"), show_col_types = FALSE)
S <- readr::read_csv(file.path(KEY_DIR, "spread_kinematics_summary.csv"), show_col_types = FALSE)
stopifnot(nrow(W) >= 4,
          all(c("week", "truncated", "new_zones", "reach_max_km", "reach_p90_km", "reach_med_km",
                "new_max_km", "case_mean_km", "cog_speed_km_wk", "inv_cog_speed_km_wk") %in% names(W)))

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
# House style (palette, province scale, themes, save_dual) from figure_style.R — the
# single definition for the whole suite.
source(file.path(ST_DIR, "figure_style.R"))

# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("make_spread_kinematics_compact.R")) {
  message("[skc] sourced, not run: nothing is written. Use `Rscript make_spread_kinematics_compact.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

CASE_COL <- "#EE7733"   # anything CASE-weighted (that week's confirmed cases)
INV_COL  <- "#009988"   # anything about the zones NEWLY invaded that week
EXT_COL  <- "#0077BB"   # the cumulative invaded set — its outer envelope
EXT_COL2 <- "#33BBEE"   # lighter shade — the middle of that same cumulative set
FONT_SCALE <- 1.85      # one knob for every text size in the figure

BASE <- 8.6 * FONT_SCALE


# Shared scaffolding: the week axis and the hollow/solid point convention, identical in
# all three panels so one caption clause covers them. Which weeks get the truncated
# treatment is read from the weekly CSV's `truncated` column; FIT_MAX (the last untruncated
# week) only bounds Panel C's reference segment. There is no shaded band anywhere here.
WK_SCALE <- scale_x_continuous(breaks = scales::breaks_width(2),
                               expand = expansion(mult = c(0.02, 0.03)))
SHAPE_TRUNC <- scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none")

# -----------------------------------------------------------------------------
# Panel A — reach from the index zone
# -----------------------------------------------------------------------------
build_reach <- function() {
  # The 90th-percentile series is DRAWN because it is the envelope the summary CSV's saturation
  # verdict and its four "Envelope (p90 reach)" rows are computed on; the running maximum alone
  # jumps whenever ONE zone sets a new record and would show a reader a different envelope
  # from the one the verdict describes. Same colour as the maximum (both are the outer edge of
  # the invaded set), dashed to tell them apart.
  lv <- c("Furthest zone so far", "90th-percentile zone so far", "Median zone so far",
          "Furthest new zone that week", "Weekly cases, mean")
  RCOL <- setNames(c(EXT_COL, EXT_COL, EXT_COL2, INV_COL, CASE_COL), lv)
  RLTY <- setNames(c("solid", "22", "solid", "solid", "solid"), lv)
  d <- W %>%
    dplyr::transmute(week, truncated,
                     `Furthest zone so far`        = reach_max_km,
                     `90th-percentile zone so far` = reach_p90_km,
                     `Median zone so far`          = reach_med_km,
                     `Furthest new zone that week` = new_max_km,
                     `Weekly cases, mean`          = case_mean_km) %>%
    tidyr::pivot_longer(-c(week, truncated), names_to = "series", values_to = "km") %>%
    # Undefined weeks (no new invasion that week) are KEPT as NA so the line breaks there;
    # filtering them out drew the line straight across a week with no new zone.
    dplyr::mutate(km = ifelse(is.finite(km), km, NA_real_),
                  series = factor(series, levels = lv))
  ggplot(d, aes(week, km, colour = series, linetype = series)) +
    geom_line(aes(group = series), linewidth = 0.9, na.rm = TRUE) +
    geom_point(aes(shape = truncated), size = 2.6, stroke = 0.9, fill = "white", na.rm = TRUE) +
    SHAPE_TRUNC +
    # untitled legend (house style: the axis title names the quantity); colour and linetype
    # share name and breaks, so they merge into ONE legend. THREE rows for five keys: on two
    # rows the third column ran across into Panel B's legend in the assembled figure.
    scale_colour_manual(name = NULL, values = RCOL,
                        guide = guide_legend(nrow = 3, byrow = TRUE)) +
    scale_linetype_manual(name = NULL, values = RLTY,
                          guide = guide_legend(nrow = 3, byrow = TRUE)) +
    WK_SCALE + expand_limits(y = 0) +
    labs(x = "Epidemic week", y = "Road distance from index zone (km)") +
    theme_pub(BASE) +
    theme(legend.key.width = unit(26, "pt"))   # wide enough to show the dash
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
    dplyr::mutate(series = factor(series, levels = lv))
  # A zero step (the centre did not move) has no place on a log axis. It is set to NA — the
  # line breaks there, as at an undefined week — and the count is reported, never hidden.
  .nonpos <- is.finite(d$kmwk) & d$kmwk <= 0
  if (any(.nonpos))
    message(sprintf("[skc] Panel B: %d zero centre-of-gravity step(s) cannot be drawn on the log axis and are left as gaps",
                    sum(.nonpos)))
  d$kmwk[!is.finite(d$kmwk) | .nonpos] <- NA_real_
  stopifnot(any(is.finite(d$kmwk)))
  # Log y: the invasion centre moves one to two orders of magnitude faster than the case
  # centre, so a linear axis flattens the case series onto the baseline. coord_cartesian
  # (which clips) rather than scale limits (which drop data) pins the window.
  rng <- range(d$kmwk, na.rm = TRUE); ylo <- rng[1] / 1.2; yhi <- rng[2] * 1.2
  # 1-3-10 breaks up to 1000: the new-invasion series reaches several hundred km/week, and a
  # top label of 300 left everything above it unreadable.
  brk <- c(1, 3, 10, 30, 100, 300, 1000)
  ggplot(d, aes(week, kmwk, colour = series)) +
    geom_hline(yintercept = V_COG, colour = CASE_COL, linewidth = 0.7, linetype = "22") +
    geom_hline(yintercept = V_INV, colour = INV_COL,  linewidth = 0.7, linetype = "22") +
    geom_line(aes(group = series), linewidth = 0.9, na.rm = TRUE) +
    geom_point(aes(shape = truncated), size = 2.6, stroke = 0.9, fill = "white", na.rm = TRUE) +
    SHAPE_TRUNC +
    scale_colour_manual(name = NULL, values = CCOL,
                        guide = guide_legend(nrow = 2, byrow = TRUE)) +
    WK_SCALE +
    scale_y_log10(breaks = brk, labels = format(brk, big.mark = ",", trim = TRUE)) +
    annotation_logticks(sides = "l", colour = FAINT, linewidth = 0.25,
                        short = unit(2, "pt"), mid = unit(3.2, "pt"), long = unit(4.6, "pt")) +
    coord_cartesian(ylim = c(ylo, yhi)) +
    # two lines: on one line the title ran into Panel A's legend in the assembled figure
    labs(x = "Epidemic week", y = "Centre-of-gravity speed\n(km per week, log scale)") +
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
  # level over weeks that never entered it. Its 90% week-bootstrap interval was drawn as a
  # shaded band and has been removed with the other rectangles; it is published as the lo/hi
  # of the same summary row (NZ_LO/NZ_HI here) and stated in the subtitle below.
  ggplot(d, aes(week, new_zones)) +
    annotate("segment", x = 0.5, xend = FIT_MAX + 0.5, y = NZ_MEAN, yend = NZ_MEAN,
             colour = INV_COL, linewidth = 0.8, linetype = "22") +
    geom_col(aes(alpha = truncated), fill = INV_COL, width = 0.7) +
    WK_SCALE +
    scale_y_continuous(breaks = scales::breaks_width(3),
                       expand = expansion(mult = c(0, 0.06))) +
    # The dashed reference level states its spread on the panel, from the published 90%
    # week-bootstrap interval, rather than leaving a bare line to be read as exact. Both numbers
    # come from spread_kinematics_summary.csv. It is NOT called a CI: weeks are resampled as
    # exchangeable although new_zones trends and is serially correlated (43's basis column).
    labs(x = "Epidemic week", y = "New health zones invaded",
         subtitle = sprintf("mean %.1f new zones/week, weeks 1-%d\n(90%% week-bootstrap interval %.1f-%.1f)",
                            NZ_MEAN, FIT_MAX, NZ_LO, NZ_HI)) +
    theme_pub(BASE) +
    # This panel has nothing to key, but the other two carry a legend above the plot area
    # (three rows in A, two in B). Dropping the legend outright would let this panel grow taller than its
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

}   # end of the script-run guard

