# =============================================================================
# 37_cascade_figure4.R — CASCADE version of publication Figure 4
# BDBV 2026 DRC · 3-month spatial invasion cascade
#
# A faithful cascade analogue of make_publication_figures.R's Figure 4:
#   Panel A  log-scaled RELATIVE invasion-risk map (rr01_nat, 0-1) for horizon H
#   Panel B  top-20 at-risk zones by H-week reach probability with 90% CrI,
#            every zone coloured by its province (labels coloured to match)
# Combined A|B at 60:40, tagged A/B, no embedded titles (house convention).
# Reads the saved cascade reach scores (primary = status-quo scenario) under OUT_DIR, and draws
# them only when they are the LATEST cascade run's for this analysis date
# (cascade_stamp_check(), 35_cascade_viz.R); an uncalibrated run's figure says so.
#
# Run:  Rscript spatiotemporal/37_cascade_figure4.R   (run_cascade.R source()s it)
# Out:  outputs/key_outputs/figures/Figure4_cascade_h<H>.{pdf,png} for each horizon in the reach
#       table; the retained-figure gate (FIGURE_KEEP / FIGURE_DROP, 00_config.R) keeps only the
#       headline horizon, which is also copied to figures/Figure4_cascade.{pdf,png}. The other
#       horizons and the panels/F4A_cascade_reach_map_h*, F4B_cascade_rank_h* parts are written
#       only with FIGURE_KEEP_ALL=1.
#       outputs/key_outputs/Figure4_cascade_h<H>_data.csv (+ _total.csv) for every horizon
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(patchwork)
  library(sf); library(scales); library(here)
})
ST_DIR <- Sys.getenv("CASCADE_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
source(file.path(ST_DIR, "00_config.R"))

# House style: palette, themes and save_dual() all come from figure_style.R, which is the
# single definition for the whole suite. This file used to carry its own copies, and its
# PROV_COL disagreed with the main-text figures' on Sud-Kivu.
source(file.path(ST_DIR, "figure_style.R"))
# The provenance helpers (cascade_stamp_check(), cascade_calibration_note()) live beside the
# reach-table writer that stamps what they check. Sourcing 35 only defines functions; when
# run_cascade.R source()s this file they are already loaded.
if (!exists("cascade_stamp_check", mode = "function"))
  source(file.path(ST_DIR, "35_cascade_viz.R"))


# ---- load shapefile + cascade reach scores ---------------------------------
shp <- st_read(SHAPEFILE_PATH, quiet = TRUE) %>%
  mutate(.key = tolower(trimws(Nom)), .prov = as.character(PROVINCE))
FIG_DIR <- file.path(OUT_DIR, "key_outputs", "figures")
PANEL_DIR <- file.path(FIG_DIR, "panels")
KEY_DIR <- file.path(OUT_DIR, "key_outputs")
reach_csv <- file.path(OUT_DIR, "cascade", "tables", "cascade_reach_scores_all_zones.csv")
if (!file.exists(reach_csv)) stop("cascade reach scores not found: ", reach_csv,
                                  " — run run_cascade.R first.")
rs <- read_csv(reach_csv, show_col_types = FALSE)
stopifnot(all(c("health_zone","province","horizon","was_active_before",
                "p_case_invasion","p_case_lo","p_case_hi","rr01_nat") %in% names(rs)))

# ---- PROVENANCE: only the LATEST, COMPLETED (or this very) cascade run for THIS date ----------
# This script used to draw whatever reach table was on disk. A cascade run whose calibration
# fails stops BEFORE it rewrites that table, and run_all.R's stage C runs this script anyway, so
# Figure4_cascade was redrawn from the previous run's table, re-stamped with a current mtime and
# exited 0. The table now carries the run id that wrote it (cascade_write_reach_csv()), and it is
# drawn only when cascade_run_stamp.json says that run is the latest for this analysis date and
# completed — or is the run_cascade.R process drawing its own figures (status "running", its id
# passed in; cascade_live_run_id()). The fan-chart table read below is reset and rewritten by the
# same run AFTER the reach table, so a reach table that passes vouches for it too.
.prov <- cascade_stamp_check(rs$cascade_run_id, file.path(OUT_DIR, "cascade", "diagnostics"))
if (!isTRUE(.prov$ok))
  stop("[fig4] refusing to draw Figure4_cascade: ", .prov$why, ". Re-run run_cascade.R.",
       call. = FALSE)
CASC_RUN_ID <- .prov$stamp$cascade_run_id
# NULL when the run was calibrated; otherwise the note every figure carries as its caption (the
# data files carry the `calibrated` flag itself).
CAL_NOTE <- cascade_calibration_note(rs$calibrated, rs$calibration_status)
if (!is.null(CAL_NOTE)) message("[fig4] ", CAL_NOTE)

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
                         # the scale is floored at RR_FLOOR, so the first tick is an upper bound
                         breaks = c(0.001, 0.01, 0.1, 1), labels = c("<=0.001","0.01","0.1","1"),
                         na.value = NA_FILL,
                         name = sprintf("Relative\ninvasion risk\n(%s, log, 0-1)", wk_label),
                         direction = 1) +
    theme_map(11)

  # --- Panel B: top-N zones by H-week reach + 90% CrI, coloured by province ---
  d <- rs %>% filter(horizon == H, !(was_active_before %in% TRUE), is.finite(p_case_invasion)) %>%
    arrange(desc(p_case_invasion)) %>% head(top_n)
  # THE SUITE'S ONE PROVINCE SCALE (prov_palette(), figure_style.R), over the provinces this
  # panel draws. The private loop this replaced took its spare hues from OKABE minus PROV_COL
  # (two colours), so a third province without an identity hue recycled the first, and its
  # spares disagreed with every other figure's.
  prov_cols   <- prov_palette(d$province)
  prov_levels <- names(prov_cols)
  d <- d %>% mutate(province = factor(province, levels = prov_levels),
                    zone = factor(health_zone, levels = rev(health_zone)))
  lev <- levels(d$zone)
  ylab_cols <- unname(prov_cols[as.character(d$province)[match(lev, d$health_zone)]])
  # A MONTE-CARLO-FLOOR interval (p_ci_at_mc_floor, cascade_reach_table()) is p +- 1/(2D):
  # its width is the resolution of D parameter draws by construction, not an estimated
  # credible interval, so it is NOT drawn as one. Its point stays; the axis title says why
  # the bar is missing whenever any is.
  .floor <- if ("p_ci_at_mc_floor" %in% names(d)) as.logical(d$p_ci_at_mc_floor) %in% TRUE
            else rep(FALSE, nrow(d))
  # Line-broken: the panel is 40% of the figure's width and one line was clipped.
  x_lab <- paste0(sprintf("Probability of invasion within %d weeks (cascade)\npoint = mean, bar = 90%% CrI", H),
                  if (any(.floor)) sprintf("\n%d without a bar: below Monte-Carlo resolution", sum(.floor)) else "")
  p_forest <- ggplot(d, aes(p_case_invasion, zone, colour = province)) +
    # EVERY ROW, with a floored interval blanked to NA, NOT a row subset: the first layer trains
    # the discrete y axis, so a subset placed its rows first and appended the others on top --
    # scrambling the descending order and misaligning the coloured zone names (ylab_cols follows
    # the factor order).
    geom_linerange(data = d %>% mutate(p_case_lo = ifelse(.floor, NA_real_, p_case_lo),
                                       p_case_hi = ifelse(.floor, NA_real_, p_case_hi)),
                   aes(xmin = p_case_lo, xmax = p_case_hi),
                   linewidth = 1.5, alpha = 0.5, na.rm = TRUE) +
    # The factor order, explicitly: the axis must follow `lev`, which ylab_cols is built on.
    scale_y_discrete(limits = lev) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = prov_cols, name = "Province", drop = FALSE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0.01, 0.06))) +
    labs(x = x_lab, y = NULL) +
    # Two rows: on one row the fifth province's key ran off the panel edge with no label.
    guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                                 override.aes = list(linewidth = 2.4, size = 3))) +
    theme_pub(11) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(colour = ylab_cols, size = 9.5, face = "bold"),
          legend.position = "top")

  # --- data underlying both panels -> key_outputs/Figure4_cascade_h{H}_data.csv ---
  rank_tbl <- rs %>% filter(horizon == H, !(was_active_before %in% TRUE), is.finite(p_case_invasion)) %>%
    arrange(desc(p_case_invasion)) %>%
    transmute(health_zone, province, invasion_prob_rank = row_number())
  # `p_interval_at_mc_floor` says which per-zone intervals are the Monte-Carlo resolution floor
  # (width at most 1/D: exactly 1/D unless clipped at 0 or 1) rather than an estimated credible
  # interval (cascade_reach_table()). NA when the reach table predates the flag.
  fig4_data <- rs %>% filter(horizon == H) %>%
    left_join(rank_tbl, by = c("health_zone", "province")) %>%
    mutate(in_top20_panelB = !is.na(invasion_prob_rank) & invasion_prob_rank <= top_n,
           rr_log10_fill = ifelse(was_active_before %in% TRUE, NA_real_, pmax(rr01_nat, RR_FLOOR)),
           p_interval_at_mc_floor = if ("p_ci_at_mc_floor" %in% names(.)) as.logical(p_ci_at_mc_floor)
                                    else NA) %>%
    arrange(invasion_prob_rank) %>%
    transmute(invasion_prob_rank, in_top20_panelB, health_zone, province, horizon,
              was_active_before, p_invasion = p_case_invasion, p_invasion_lo = p_case_lo,
              p_invasion_hi = p_case_hi, p_interval_at_mc_floor,
              relative_invasion_risk = rr01_nat, rr_log10_fill,
              # provenance: the run the figure was drawn from, and whether its level is calibrated
              cascade_run_id = CASC_RUN_ID,
              calibrated = if ("calibrated" %in% names(.)) as.logical(calibrated) else NA)
  if (!dir.exists(KEY_DIR)) dir.create(KEY_DIR, recursive = TRUE, showWarnings = FALSE)
  write_csv(fig4_data, file.path(KEY_DIR, sprintf("Figure4_cascade_h%d_data.csv", H)))

  # THE TOTAL, beside the per-zone data — because a total is what a reader builds from this
  # file, and the per-zone intervals CANNOT be summed into one (quantiles are not additive; a
  # manuscript audit found exactly that). The expected number of new zones is additive (the sum
  # of the per-zone p); its interval is read from the simulated national total across
  # Monte-Carlo iterations (cascade_newzones_fanchart.csv, primary scenario, week H). Written
  # only when that table is present; otherwise any older copy is removed rather than left stale.
  tot_path <- file.path(KEY_DIR, sprintf("Figure4_cascade_h%d_total.csv", H))
  fan_csv  <- file.path(OUT_DIR, "cascade", "tables", "cascade_newzones_fanchart.csv")
  fan <- if (file.exists(fan_csv)) read_csv(fan_csv, show_col_types = FALSE) else NULL
  # The primary scenario is flagged IN the table by run_cascade.R (`is_primary`), so this script
  # needs no scenario config of its own (it sources only 00_config.R).
  fw  <- if (!is.null(fan) && all(c("week", "scenario", "is_primary", "med", "lo", "hi") %in% names(fan)))
           fan[fan$week == H & fan$is_primary %in% TRUE, , drop = FALSE] else NULL
  if (!is.null(fw) && nrow(fw) == 1L) {
    write_csv(tibble::tibble(
        horizon = H, scenario = fw$scenario,
        expected_new_zones = sum(rs$p_case_invasion[rs$horizon == H & !(rs$was_active_before %in% TRUE)],
                                 na.rm = TRUE),
        mean_simulated_total = if ("mean" %in% names(fw)) fw$mean else NA_real_,
        median_new_zones = fw$med, lo_90 = fw$lo, hi_90 = fw$hi,
        interval = "5th-95th percentile of the SIMULATED national total across Monte-Carlo iterations (NOT a sum of per-zone intervals)",
        cascade_run_id = CASC_RUN_ID,
        calibrated = { .c <- unique(as.logical(rs$calibrated)); if (length(.c) == 1L) .c else NA }),
      tot_path)
  } else if (file.exists(tot_path)) {
    file.remove(tot_path)
    message("[fig4] ", basename(tot_path), " not written (no fan-chart table for week ", H,
            "); removed the stale copy.")
  }

  save_dual(p_map,    sprintf("F4A_cascade_reach_map_h%d", H), 6.4, 6.4, PANEL_DIR)
  save_dual(p_forest, sprintf("F4B_cascade_rank_h%d", H),      6.6, 6.4, PANEL_DIR)
  # An uncalibrated run is labelled ON the figure (CAL_NOTE; NULL, i.e. no caption, otherwise).
  # theme_pub() blanks plot.caption, so the annotation sets it explicitly.
  fig4 <- (wrap_elements(full = p_map) | wrap_elements(full = p_forest)) +
    plot_layout(widths = c(0.6, 0.4)) +
    plot_annotation(tag_levels = "A", caption = CAL_NOTE,
                    theme = theme(plot.caption = element_text(size = 8, colour = MUTED, hjust = 0)))
  save_dual(fig4, sprintf("Figure4_cascade_h%d", H), 12.8, 6.6, FIG_DIR)
  invisible(fig4)
}

# ---- RUN for each reporting horizon (the deepest = the headline 3-month figure) --
# READ FROM THE TABLE. CASCADE_REPORT_HORIZONS lives in 30_projection_config.R, which this
# script does not source, so a get0() of it always fell back to a typed literal. The reach
# table holds exactly the horizons the run reported (cascade_reach_table()'s default).
horizons <- sort(unique(as.integer(rs$horizon)))
stopifnot(length(horizons) >= 1L, !anyNA(horizons))
.t_start <- Sys.time()
.built <- stats::setNames(logical(length(horizons)), as.character(horizons))
for (H in horizons) {
  message(sprintf("== Figure 4 (cascade, %d-week reach + uncertainty) ==", H))
  .built[as.character(H)] <- tryCatch({ build_fig4_cascade(H); TRUE },
    error = function(e) { message("!! Figure4_cascade_h", H, " FAILED: ",
                                  conditionMessage(e)); FALSE })
}

# Copy the headline (13-week) figure to the canonical Figure4_cascade name.
# ONLY when this run actually produced it. The build is wrapped in tryCatch, so a failure was
# previously swallowed and this loop then found LAST run's Figure4_cascade_h13.pdf still on
# disk and republished it as the canonical headline — shipping a stale figure under a
# current-looking name, with "[done]" printed and nothing flagged anywhere. Two conditions are
# required: the builder returned without error, AND the source file was written during THIS
# run (mtime at or after .t_start), which also catches a build that "succeeded" without
# writing because the retained-figure gate suppressed it.
.H_head <- max(horizons)
for (ext in c("pdf", "png")) {
  src <- file.path(FIG_DIR, sprintf("Figure4_cascade_h%d.%s", .H_head, ext))
  dst <- file.path(FIG_DIR, sprintf("Figure4_cascade.%s", ext))
  fresh <- isTRUE(.built[as.character(.H_head)]) && file.exists(src) &&
           file.mtime(src) >= .t_start - 1
  if (fresh) {
    file.copy(src, dst, overwrite = TRUE)
  } else if (file.exists(dst)) {
    # A canonical file that this run did not regenerate is stale, and everything downstream
    # (the cascade report, the key-outputs manifest) reads it as current. Remove it so the
    # failure is loud and detectable rather than silently published. The per-horizon source
    # file is untouched, so nothing is lost.
    file.remove(dst)
    warning(sprintf("[fig4] Figure4_cascade_h%d was not produced by this run; removed the STALE %s rather than republish it.",
                    .H_head, basename(dst)), call. = FALSE, immediate. = TRUE)
  }
}
if (!isTRUE(.built[as.character(.H_head)]))
  message("[fig4] headline horizon h", .H_head, " did NOT build — canonical Figure4_cascade not written.")
message("[done] Figure4_cascade -> ", FIG_DIR)
