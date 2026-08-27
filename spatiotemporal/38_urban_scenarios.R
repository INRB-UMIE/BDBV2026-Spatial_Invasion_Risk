# =============================================================================
# 38_urban_scenarios.R — CASCADE: urban-hub invasion scenarios (timing-aware)
# BDBV 2026 DRC · implements the "what if a major city is invaded" suite.
#
# Generalises the Kisangani conditional analysis to the country's largest, most
# connected metros. For each hub city and each seeding WEEK, the hub's urban-core
# health zones are FORCE-SEEDED at that week; the cascade then projects the hub's
# onward spread through the mobility network. Impact is measured against the
# baseline (no urban seeding) cascade:
#   - per-zone attributable increase   dP_i = P_cond(i) - P_base(i)
#   - relative risk                    RR_i = P_cond(i) / P_base(i)
#   - threshold crossings              # zones the seeding pushes above P>t (material,
#                                      non-noise upward crossings), and % of baseline
#   - expected additional invasions    E[added] = sum_i dP_i  (i != hub; noise-robust)
#   - provinces newly exposed
# TIMING is a first-class axis: earlier seeding leaves more weeks for onward
# spread, so impact by week 13 is larger for earlier seeds.
#
# Sourced AFTER 30-35 (uses simulate_cascade, cascade_reach_table, and the 35
# shapefile/plot helpers .cascade_shape/.cascade_join/.cascade_ggsave/OUT_CASCADE).
# ASSUMPTIONS documented inline. All quantities are downstream of, and consistent
# with, the confirmed-case-scale cascade (time index = symptom onset).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

# Urban-core health zones per city (verified against the WorldPop spine). "[City]
# invaded" is operationalised as seeding these central/most-connected zones; the
# cascade then carries spread to the rest of the metro and its catchment.
URBAN_HUBS <- list(
  Kinshasa     = c("Gombe", "Limete", "Kimbanseke"),
  Lubumbashi   = c("Lubumbashi", "Kampemba", "Ruashi"),
  `Mbuji-Mayi` = c("Diulu", "Bonzola", "Nzaba"),
  Kananga      = c("Kananga", "Katoka", "Lukonga"),
  Tshikapa     = c("Tshikapa", "Kanzala"))
URBAN_SEED_WEEKS <- c(1L, 4L, 8L)          # timing axis (weeks into the projection)
URBAN_THRESHOLDS <- c(0.10, 0.20, 0.50)    # includes the P>50% example metric
URBAN_ATTR_FLOOR <- 5e-4                   # RR undefined below this baseline prob (0/0 guard)
# Urban runs use a HIGHER M than the main cascade: the baseline and conditional are
# separate stochastic runs whose streams diverge after the seed week, so far zones
# (unaffected by this LOCAL seeding) show a per-zone attributable-difference NOISE that
# must sit below the material threshold. Empirically the pure noise SD is ~0.008 at
# M=2000 (q99 ~0.028; a few far zones flicker above 0.03) and drops to ~0.006 at M=4000
# (0 zones above 0.03). At M=4000, delta_min=0.03 is well above the noise, so it keeps
# real moderate-catchment zones while excluding flicker. The expected-additional-
# invasions metric (sum of dP over ALL zones) is noise-robust and is the headline.
URBAN_N_MC <- as.integer(Sys.getenv("URBAN_N_MC",
                unset = if (identical(Sys.getenv("CASCADE_SMOKE"), "1")) "40" else "4000"))
URBAN_DELTA_MIN <- 0.03

# Validate hub zones against the spine; drop (and warn on) any that are absent.
resolve_urban_hubs <- function(hubs, zones_all) {
  out <- lapply(names(hubs), function(nm) {
    ex <- intersect(hubs[[nm]], zones_all); miss <- setdiff(hubs[[nm]], zones_all)
    if (length(miss)) message("[urban] ", nm, ": not in spine, dropped: ",
                              paste(miss, collapse = ", "))
    ex
  })
  names(out) <- names(hubs)
  out[vapply(out, length, integer(1)) > 0]
}

# Run one urban scenario: force-seed hub_zones at seed_week; return reach + sim.
run_urban_scenario <- function(prep, scenario, hub_zones, seed_week, delta, psi, n_mc) {
  sim <- simulate_cascade(prep, scenario, n_mc = n_mc, delta = delta, psi = psi,
                          force_seed = hub_zones, force_seed_week = seed_week)
  list(sim = sim, reach = cascade_reach_table(sim, horizons = CASCADE_REPORT_HORIZONS))
}

# Impact of a scenario vs baseline, at one horizon. The seeded hub zones and any
# already-affected (masked) zones are EXCLUDED from the impact summaries (the hub
# is invaded by assumption, not prediction), so metrics reflect DOWNSTREAM effect.
urban_impact_metrics <- function(reach_base, reach_cond, hub_zones, province_map,
                                 horizon = max(CASCADE_REPORT_HORIZONS),
                                 thresholds = URBAN_THRESHOLDS,
                                 attr_floor = URBAN_ATTR_FLOOR,
                                 delta_min = URBAN_DELTA_MIN) {
  gb <- function(rt) rt[rt$horizon == horizon, c("health_zone", "p_case_invasion")]
  b <- gb(reach_base); names(b)[2] <- "p_base"
  c2 <- gb(reach_cond); names(c2)[2] <- "p_cond"
  m <- merge(b, c2, by = "health_zone")
  m <- m[!(m$health_zone %in% hub_zones) & is.finite(m$p_base) & is.finite(m$p_cond), ]
  m$delta <- m$p_cond - m$p_base
  m$rr    <- ifelse(m$p_base > attr_floor, m$p_cond / m$p_base, NA_real_)
  if (!is.null(province_map))
    m$province <- province_map$province[match(m$health_zone, province_map$nom)]
  material <- is.finite(m$delta) & m$delta > delta_min   # zones the seeding materially elevated
  thr <- dplyr::bind_rows(lapply(thresholds, function(t) {
    # zones pushed from <=t to >t by the seeding, with a material (non-noise) change;
    # far zones unaffected by this local seeding are excluded from the count.
    newly <- sum(m$p_base <= t & m$p_cond > t & material, na.rm = TRUE)
    nb    <- sum(m$p_base > t, na.rm = TRUE)
    tibble::tibble(threshold = t, n_base = nb, newly_above = newly,
                   pct_increase = 100 * newly / max(nb, 1))
  }))
  # provinces newly exposed: a province with a materially-elevated zone above 0.20 that
  # had no baseline zone above 0.20.
  prov_new <- NA_integer_
  if (!is.null(province_map)) {
    pb <- unique(m$province[m$p_base > 0.20])
    pc <- unique(m$province[m$p_cond > 0.20 & material])
    prov_new <- length(setdiff(pc[!is.na(pc)], pb[!is.na(pb)]))
  }
  list(per_zone = m[order(-m$delta), ],
       thresholds = thr,
       expected_added = sum(m$delta, na.rm = TRUE),          # E[extra zones invaded by H]
       n_elevated = sum(material, na.rm = TRUE),             # # zones materially elevated
       max_rr = suppressWarnings(max(m$rr[material], na.rm = TRUE)),
       provinces_newly_exposed = prov_new,
       horizon = horizon)
}

# Run the full hub x seed-week grid. Returns a per-scenario summary table plus,
# for the EARLIEST seed week (max impact), the per-zone impact + sim for detailed
# figures. reach_base = the baseline (no urban seeding) reach table.
urban_scenario_grid <- function(prep, scenario, hubs, seed_weeks, delta, psi, n_mc,
                                reach_base, province_map) {
  summ <- list(); detail <- list()
  for (city in names(hubs)) {
    message("[urban]   scenario: ", city, " (", length(seed_weeks), " seed-weeks, M=", n_mc, ")")
    for (sw in seed_weeks) {
      sc <- run_urban_scenario(prep, scenario, hubs[[city]], sw, delta, psi, n_mc)
      im <- urban_impact_metrics(reach_base, sc$reach, hubs[[city]], province_map)
      t50 <- im$thresholds[im$thresholds$threshold == 0.50, ]
      t20 <- im$thresholds[im$thresholds$threshold == 0.20, ]
      summ[[paste(city, sw)]] <- tibble::tibble(
        city = city, seed_week = sw, n_hub_zones = length(hubs[[city]]),
        expected_added = im$expected_added, n_elevated = im$n_elevated,
        newly_above_20 = t20$newly_above, n_base_gt50 = t50$n_base,
        newly_above_50 = t50$newly_above, pct_increase_gt50 = t50$pct_increase,
        max_rr = im$max_rr, provinces_newly_exposed = im$provinces_newly_exposed)
      if (sw == min(seed_weeks)) {
        im$per_zone$city <- city; im$per_zone$seed_week <- sw
        detail[[city]] <- list(impact = im, reach = sc$reach)
      }
    }
  }
  list(summary = dplyr::bind_rows(summ), detail = detail,
       seed_weeks = seed_weeks, min_week = min(seed_weeks), n_mc = n_mc)
}

# ---- visualisations --------------------------------------------------------
# 1. Attributable-risk choropleth: where the hub's onward spread concentrates.
urban_map_attributable <- function(impact, hub_label, seed_week,
                                   file_prefix = "urban_attributable") {
  d <- impact$per_zone
  # Show only MATERIALLY-elevated zones (delta > URBAN_DELTA_MIN); grey out the rest,
  # which carry only Monte-Carlo flicker (this local seeding does not reach them). This
  # matches the impact metric and removes visual noise so the real catchment stands out.
  d$delta_show <- ifelse(is.finite(d$delta) & d$delta > URBAN_DELTA_MIN, d$delta, NA_real_)
  sp <- .cascade_join(d, "delta_show"); if (is.null(sp)) return(invisible(NULL))
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = "grey80", linewidth = 0.04) +
    scale_fill_viridis_c(option = "rocket", direction = -1, na.value = "grey93",
                         name = "Attributable\nincrease in\nP(invasion)", limits = c(0, NA)) +
    labs(title = sprintf("BDBV 2026 — downstream invasion risk attributable to invasion of %s",
                         hub_label),
         subtitle = sprintf("Materially-elevated zones only: increase in 13-week reach vs baseline if %s is seeded at week %d (projection; upper bound).",
                            hub_label, seed_week)) +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7.5))
  .cascade_ggsave(p, sprintf("%s_%s", file_prefix, gsub("[^A-Za-z0-9]+", "_", hub_label)))
}

# 2. Top downstream zones by attributable increase (with relative risk).
urban_fig_top_zones <- function(impact, hub_label, top_n = 15L,
                                file_prefix = "urban_top_zones") {
  d <- impact$per_zone[is.finite(impact$per_zone$delta) &
                       impact$per_zone$delta > URBAN_DELTA_MIN, ]   # material zones only
  if (!nrow(d)) return(invisible(NULL))
  d <- head(d[order(-d$delta), ], top_n)
  d$health_zone <- factor(d$health_zone, levels = rev(d$health_zone))
  p <- ggplot(d, aes(delta, health_zone)) +
    geom_col(fill = "#b2182b") +
    geom_text(aes(label = ifelse(is.finite(rr), sprintf("RR %.1f", rr), "")),
              hjust = -0.1, size = 3, colour = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.18))) +
    labs(title = sprintf("BDBV 2026 — zones most affected by invasion of %s", hub_label),
         subtitle = "Attributable increase in 13-week reach probability; RR = conditional / baseline.",
         x = "attributable increase in P(invasion)", y = NULL) +
    theme_minimal(base_size = 11)
  .cascade_ggsave(p, sprintf("%s_%s", file_prefix, gsub("[^A-Za-z0-9]+", "_", hub_label)),
                  w = 8, h = 6)
}

# 3. Cross-hub summary bar at the earliest seed week (expected additional invasions).
urban_fig_summary_bar <- function(summary, min_week, file = "urban_impact_summary") {
  d <- summary[summary$seed_week == min_week, ]
  d$city <- factor(d$city, levels = d$city[order(d$expected_added)])
  p <- ggplot(d, aes(expected_added, city)) +
    geom_col(fill = "#2166ac") +
    geom_text(aes(label = sprintf("%d zones elevated; +%d newly >20%%",
                                  n_elevated, newly_above_20)),
              hjust = -0.02, size = 2.9, colour = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.55))) +
    labs(title = sprintf("BDBV 2026 — downstream impact if a major city is invaded (seeded week %d)", min_week),
         subtitle = "Expected additional health zones invaded within 13 weeks (excludes the seeded city itself).",
         x = "expected additional zones invaded (13 weeks)", y = NULL) +
    theme_minimal(base_size = 11)
  .cascade_ggsave(p, file, w = 9, h = 5.5)
}

# 4. Timing sensitivity: impact vs seeding week, per hub.
urban_fig_timing <- function(summary, file = "urban_timing_sensitivity") {
  p <- ggplot(summary, aes(seed_week, expected_added, colour = city, group = city)) +
    geom_line(linewidth = 1) + geom_point(size = 2.4) +
    scale_x_continuous(breaks = sort(unique(summary$seed_week))) +
    labs(title = "BDBV 2026 — earlier invasion of a city drives more onward spread",
         subtitle = "Expected additional zones invaded by 13 weeks, by the week the city is seeded. Later seeding leaves fewer weeks for spread.",
         x = "week the city is invasion-seeded", y = "expected additional zones (by 13 wk)",
         colour = "City") +
    theme_minimal(base_size = 11)
  .cascade_ggsave(p, file, w = 9, h = 5.5)
}

# Write per-hub per-zone impact CSVs + the cross-scenario summary table.
urban_write_tables <- function(grid, dir) {
  readr::write_csv(grid$summary, file.path(dir, "urban_scenario_summary.csv"))
  for (city in names(grid$detail)) {
    d <- grid$detail[[city]]$impact$per_zone
    readr::write_csv(d, file.path(dir, sprintf("urban_impact_%s.csv",
                                  gsub("[^A-Za-z0-9]+", "_", city))))
  }
  invisible(NULL)
}

# ---- publication figures (Figure4-style) -----------------------------------
# Province colour identity (extends the Figure-4 eastern palette to the western/
# southern provinces that host the urban catchments), Okabe-Ito CVD-safe.
.URBAN_OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
.URBAN_PROV_COL <- c(
  "Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73","Sud-Kivu"="#CC79A7",
  "Kinshasa"="#E69F00","Haut-Katanga"="#56B4E9","Kasaï-Oriental"="#CC79A7",
  "Kasaï-Central"="#009E73","Kasaï"="#D55E00","Kwilu"="#0072B2",
  "Kwango"="#F0E442","Lomami"="#000000","Mai-Ndombe"="#7A7A7A")
.urban_prov_palette <- function(provs) {
  provs <- unique(as.character(provs[!is.na(provs)]))
  cols <- setNames(character(length(provs)), provs)
  spare <- setdiff(.URBAN_OKABE, unname(.URBAN_PROV_COL)); si <- 1L
  for (p in provs) {
    if (p %in% names(.URBAN_PROV_COL)) cols[p] <- unname(.URBAN_PROV_COL[p])
    else { cols[p] <- spare[((si - 1L) %% max(length(spare), 1L)) + 1L]; si <- si + 1L }
  }
  cols
}
.URBAN_KEYFIG_DIR <- function() file.path(OUT_DIR, "key_outputs", "figures")

# Figure-4 analogue for a city scenario: Panel A = attributable-risk choropleth
# zoomed to the catchment (seeded core in grey), Panel B = baseline -> conditional
# reach dumbbell for the top affected zones, coloured by province. Built from the
# per-zone impact only (regenerable without re-simulating).
urban_figure4 <- function(impact, hub_label, hub_zones, top_n = 15L,
                          file_prefix = "Figure_urban") {
  if (!requireNamespace("patchwork", quietly = TRUE)) return(invisible(NULL))
  sp0 <- .cascade_shape(); if (is.null(sp0)) return(invisible(NULL))
  d <- impact$per_zone
  mat <- d[is.finite(d$delta) & d$delta > URBAN_DELTA_MIN, ]
  if (!nrow(mat)) return(invisible(NULL))
  # --- Panel A: attributable map, zoomed to catchment + seeded core ---
  d$dshow <- ifelse(is.finite(d$delta) & d$delta > URBAN_DELTA_MIN, d$delta, NA_real_)
  spm <- .cascade_join(d, "dshow"); if (is.null(spm)) return(invisible(NULL))
  hubkey <- tolower(trimws(hub_zones)); matkey <- tolower(trimws(mat$health_zone))
  focus <- spm[spm$.key %in% c(matkey, hubkey), ]
  bb <- sf::st_bbox(focus)
  mx <- 0.18 * (bb[["xmax"]] - bb[["xmin"]]); my <- 0.18 * (bb[["ymax"]] - bb[["ymin"]])
  pA <- ggplot(spm) +
    geom_sf(aes(fill = .val), colour = "grey78", linewidth = 0.06) +
    geom_sf(data = spm[spm$.key %in% hubkey, ], fill = "grey35",
            colour = "white", linewidth = 0.15) +
    scale_fill_viridis_c(option = "rocket", direction = -1, na.value = "grey92",
                         name = "Attributable\nΔ P(invasion)", limits = c(0, NA)) +
    coord_sf(xlim = c(bb[["xmin"]] - mx, bb[["xmax"]] + mx),
             ylim = c(bb[["ymin"]] - my, bb[["ymax"]] + my), expand = FALSE) +
    theme_void(base_size = 11) +
    theme(legend.position = "right", plot.margin = margin(2, 2, 2, 2))
  # --- Panel B: baseline -> conditional dumbbell, top zones by attributable increase ---
  db <- head(mat[order(-mat$delta), ], top_n)
  db$zone <- factor(db$health_zone, levels = rev(db$health_zone))
  pcols <- .urban_prov_palette(db$province)
  lev <- levels(db$zone); ylab_cols <- unname(pcols[as.character(db$province)[match(lev, db$health_zone)]])
  pB <- ggplot(db) +
    geom_segment(aes(x = p_base, xend = p_cond, y = zone, yend = zone),
                 colour = "grey78", linewidth = 1) +
    geom_point(aes(p_base, zone), colour = "grey55", size = 1.9) +
    geom_point(aes(p_cond, zone, colour = province), size = 2.8) +
    scale_colour_manual(values = pcols, name = "Province") +
    scale_x_continuous(labels = scales::percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0.01, 0.06))) +
    labs(x = "P(invasion by 13 wk): baseline (grey) → if city invaded (colour)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(colour = ylab_cols, face = "bold", size = 9),
          legend.position = "top")
  fig <- (patchwork::wrap_elements(full = pA) | patchwork::wrap_elements(full = pB)) +
    patchwork::plot_layout(widths = c(0.55, 0.45)) +
    patchwork::plot_annotation(tag_levels = "A",
      title = sprintf("BDBV 2026 — downstream invasion risk if %s is invaded (13-week projection)", hub_label),
      subtitle = "A: attributable increase in reach (seeded core in dark grey). B: baseline → conditional reach for the most-affected zones, by province.",
      theme = ggplot2::theme(plot.subtitle = element_text(size = 8)))
  path <- file.path(.URBAN_KEYFIG_DIR(), sprintf("%s_%s", file_prefix,
                    gsub("[^A-Za-z0-9]+", "_", hub_label)))
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(path, ".pdf"), fig, width = 12.5, height = 6.2, bg = "white")
  ggsave(paste0(path, ".png"), fig, width = 12.5, height = 6.2, dpi = 150, bg = "white")
  message("  saved ", basename(path)); invisible(path)
}

# Multi-city small-multiples: attributable-risk maps for all cities in one figure
# (national extent), a comprehensive spatial overview. impacts = named list of
# per_zone data frames (one per city).
urban_map_facet <- function(impacts, file = "Figure_urban_all_cities") {
  sp0 <- .cascade_shape(); if (is.null(sp0)) return(invisible(NULL))
  frames <- lapply(names(impacts), function(city) {
    d <- impacts[[city]]
    d$dshow <- ifelse(is.finite(d$delta) & d$delta > URBAN_DELTA_MIN, d$delta, NA_real_)
    s <- .cascade_join(d, "dshow"); if (is.null(s)) return(NULL)
    s$city <- city; s
  })
  frames <- frames[!vapply(frames, is.null, logical(1))]
  if (!length(frames)) return(invisible(NULL))
  sp <- do.call(rbind, frames)
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = NA) +
    facet_wrap(~city, nrow = 1) +
    scale_fill_viridis_c(option = "rocket", direction = -1, na.value = "grey93",
                         name = "Attributable\nΔ P(invasion)", limits = c(0, NA)) +
    labs(title = "BDBV 2026 — downstream invasion risk attributable to invasion of each city (13-week projection)",
         subtitle = "Increase in reach vs baseline if the city is seeded at week 1; materially-elevated zones only. Projection, upper bound.") +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7.5),
                                       strip.text = element_text(face = "bold"))
  path <- file.path(.URBAN_KEYFIG_DIR(), file)
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(path, ".pdf"), p, width = 15, height = 4.6, bg = "white")
  ggsave(paste0(path, ".png"), p, width = 15, height = 4.6, dpi = 150, bg = "white")
  message("  saved ", basename(path)); invisible(path)
}
