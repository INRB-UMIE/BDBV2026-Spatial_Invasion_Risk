# =============================================================================
# 35_cascade_viz.R — 3-MONTH CASCADE: policy-focused outputs (PLAN §7)
# BDBV 2026 DRC
#
# Enriches reach tables with province / vulnerability / relative-risk / priority
# (reusing compute_risk_scores + add_risk_indices + compute_vulnerability_index),
# writes the all-zone reach CSVs, and renders the decision figures:
#   - 3-month cumulative-reach choropleth + uncertainty (national + provincial)
#   - time-to-invasion triptych (P invaded by wk 4/8/13) — the expanding frontier
#   - aggregate fan chart (cumulative new zones over 13 wk, per scenario)
#   - scenario comparison reach maps
#   - seeding corridor / flux map + gateway (knockout Delta_j) bars
#   - conditional onward-spread map for key hubs (generalised Kisangani)
#   - establishment vs reach
# Framing stamped on every product: "3-month PROJECTION (scenario-based); trust
# rankings/relative risk; absolute probabilities are UPPER BOUNDS."
#
# Sourced AFTER 30-34, 15_workhorse.R (compute_risk_scores), 20_forecast_detail.R
# (add_risk_indices, compute_vulnerability_index, map helpers) and 17_invasion_viz.R.
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

CASCADE_FRAMING <- paste0("3-month PROJECTION (scenario-based). Trust rankings & ",
  "relative risk; absolute probabilities are upper bounds. Time index = symptom onset.")

# ---- enrichment + CSV ------------------------------------------------------
cascade_enrich_reach <- function(reach, province_map, vuln, method_label) {
  rs <- compute_risk_scores(reach, province_map)
  rs <- add_risk_indices(rs, vuln)
  rs$method <- method_label
  rs
}

cascade_write_reach_csv <- function(rs, path) {
  keep <- c("health_zone", "province", "horizon", "was_active_before",
            "p_case_invasion", "p_lo", "p_hi", "p_establishment", "first_passage_week",
            "mu_forecast", "rr_nat", "rr_nat_rank", "rr01_nat",
            "rr_ituri", "rr_ituri_rank", "rr_nordkivu", "rr_nordkivu_rank",
            "rr_hautuele", "rr_hautuele_rank", "V", "priority", "priority_rank",
            "surveillance_gap", "healthcare_gap", "access_gap", "social_vulnerability",
            "healthcare_travel_min", "method")
  out <- rs |> dplyr::select(dplyr::any_of(keep)) |>
    dplyr::rename(dplyr::any_of(c(p_case_lo = "p_lo", p_case_hi = "p_hi"))) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 4)))
  if (all(c("horizon", "p_case_invasion") %in% names(out)))
    out <- out |> dplyr::arrange(horizon, dplyr::desc(dplyr::coalesce(p_case_invasion, -1)))
  readr::write_csv(out, path)
  message(sprintf("[cascade] wrote %d rows -> %s", nrow(out), basename(path)))
  invisible(out)
}

# ---- shapefile helpers -----------------------------------------------------
.cascade_shape <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    sp <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE), error = function(e) NULL)
    if (!is.null(sp)) { sp$.key <- tolower(trimws(sp$Nom)); cache <<- sp }
    sp
  }
})
.cascade_join <- function(df, valcol) {
  sp <- .cascade_shape(); if (is.null(sp)) return(NULL)
  df$.key <- tolower(trimws(df$health_zone))
  sp$.val <- df[[valcol]][match(sp$.key, df$.key)]
  sp
}
.cascade_ggsave <- function(p, file, w = 9, h = 7.5) {
  path <- file.path(OUT_CASCADE, "figures", file)
  ggsave(paste0(path, ".pdf"), p, width = w, height = h)
  ggsave(paste0(path, ".png"), p, width = w, height = h, dpi = 150)
  invisible(path)
}
.PROV_LIMS <- get0("PROVINCES_OF_INTEREST", ifnotfound = c("Ituri","Nord-Kivu","Haut-Uele"))

# ---- 1. cumulative-reach choropleth (national) -----------------------------
cascade_map_reach <- function(rs, horizon = 13L, file = "cascade_reach_map_national",
                              scenario_label = "") {
  d <- rs[rs$horizon == horizon, ]
  sp <- .cascade_join(d, "p_case_invasion"); if (is.null(sp)) return(invisible(NULL))
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = "grey70", linewidth = 0.05) +
    scale_fill_viridis_c(option = "inferno", na.value = "grey92", direction = -1,
                         name = sprintf("P(invaded\nby %d wk)", horizon), limits = c(0, NA)) +
    labs(title = sprintf("BDBV 2026 — %d-week spatial invasion reach%s", horizon,
                         if (nzchar(scenario_label)) paste0(" [", scenario_label, "]") else ""),
         subtitle = CASCADE_FRAMING, caption = "Grey = already-affected (masked) or no data.") +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7))
  .cascade_ggsave(p, paste0(file, "_h", horizon))
}

# ---- 2. uncertainty map (credible-interval width) --------------------------
cascade_map_uncertainty <- function(rs, horizon = 13L,
                                    file = "cascade_reach_uncertainty_national") {
  d <- rs[rs$horizon == horizon, ]; d$.w <- d$p_hi - d$p_lo
  sp <- .cascade_join(d, ".w"); if (is.null(sp)) return(invisible(NULL))
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = "grey70", linewidth = 0.05) +
    scale_fill_viridis_c(option = "mako", na.value = "grey92",
                         name = "90% CrI\nwidth", direction = -1) +
    labs(title = sprintf("BDBV 2026 — %d-week reach: forecast uncertainty", horizon),
         subtitle = "Width of the 90% credible interval on the reach probability (posterior parameter uncertainty; excludes between-scenario/structural uncertainty).") +
    theme_void(base_size = 11)
  .cascade_ggsave(p, paste0(file, "_h", horizon))
}

# ---- 3. time-to-invasion triptych (expanding frontier) ---------------------
cascade_map_timetoinvasion <- function(rs, horizons = CASCADE_REPORT_HORIZONS,
                                       file = "cascade_time_to_invasion_triptych") {
  sp0 <- .cascade_shape(); if (is.null(sp0)) return(invisible(NULL))
  frames <- lapply(horizons, function(H) {
    d <- rs[rs$horizon == H, ]; s <- .cascade_join(d, "p_case_invasion")
    s$panel <- factor(sprintf("by %d weeks", H),
                      levels = sprintf("by %d weeks", horizons)); s
  })
  sp <- do.call(rbind, frames)
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = NA) + facet_wrap(~panel, nrow = 1) +
    scale_fill_viridis_c(option = "inferno", na.value = "grey92", direction = -1,
                         name = "P(invaded)", limits = c(0, NA)) +
    labs(title = "BDBV 2026 — spatial invasion: the expanding frontier",
         subtitle = CASCADE_FRAMING) +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7))
  .cascade_ggsave(p, file, w = 13, h = 5.5)
}

# ---- 4. aggregate fan chart (new zones over time, per scenario) ------------
cascade_fig_fanchart <- function(sims_by_scenario,
                                 file = "cascade_newzones_fanchart") {
  rows <- lapply(names(sims_by_scenario), function(sk) {
    sim <- sims_by_scenario[[sk]]
    cum <- apply(sim$new_by_week, 2, cumsum)               # weeks x iters
    tibble::tibble(week = seq_len(nrow(cum)),
      scenario = sim$knobs$scenario,
      med = apply(cum, 1, stats::median),
      lo  = apply(cum, 1, stats::quantile, 0.05),
      hi  = apply(cum, 1, stats::quantile, 0.95))
  })
  d <- dplyr::bind_rows(rows)
  p <- ggplot(d, aes(week, med, colour = scenario, fill = scenario)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    labs(title = "BDBV 2026 — expected NEW health zones invaded over 13 weeks",
         subtitle = "Median + 90% band, by scenario. Projection, not a point forecast.",
         x = "weeks ahead", y = "cumulative newly-invaded zones",
         colour = "Scenario", fill = "Scenario") +
    theme_minimal(base_size = 12)
  .cascade_ggsave(p, file, w = 9, h = 6)
  d
}

# ---- 5. scenario comparison reach maps -------------------------------------
cascade_map_scenarios <- function(rs_by_scenario, horizon = 13L,
                                  file = "cascade_scenario_reach_maps") {
  sp0 <- .cascade_shape(); if (is.null(sp0)) return(invisible(NULL))
  frames <- lapply(names(rs_by_scenario), function(sk) {
    d <- rs_by_scenario[[sk]]; d <- d[d$horizon == horizon, ]
    s <- .cascade_join(d, "p_case_invasion"); s$scenario <- sk; s
  })
  sp <- do.call(rbind, frames)
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = NA) + facet_wrap(~scenario, nrow = 1) +
    scale_fill_viridis_c(option = "inferno", na.value = "grey92", direction = -1,
                         name = sprintf("P(invaded\nby %dw)", horizon), limits = c(0, NA)) +
    labs(title = sprintf("BDBV 2026 — %d-week reach under control scenarios", horizon),
         subtitle = CASCADE_FRAMING) +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7))
  .cascade_ggsave(p, paste0(file, "_h", horizon), w = 14, h = 5.5)
}

# ---- 6. gateway (knockout Delta_j) bars ------------------------------------
cascade_fig_gateway <- function(gateway_tbl, file = "cascade_gateway_knockout",
                                top_n = 15L) {
  d <- head(gateway_tbl, top_n)
  d$health_zone <- factor(d$health_zone, levels = rev(d$health_zone))
  p <- ggplot(d, aes(delta_j, health_zone)) +
    geom_col(fill = "#b2182b") +
    labs(title = "BDBV 2026 — gateway zones: onward-spread impact",
         subtitle = "Delta_j = expected downstream invasions prevented if this source is contained (knockout).",
         x = "expected downstream invasions attributable (Delta_j)", y = NULL) +
    theme_minimal(base_size = 12)
  .cascade_ggsave(p, file, w = 9, h = 6)
}

# ---- 7. seeding corridors (province flux heat + top zone edges) ------------
cascade_fig_corridors <- function(flux_net, file = "cascade_corridors") {
  pm <- flux_net$province_matrix
  if (is.null(pm)) return(invisible(NULL))
  dm <- as.data.frame(as.table(pm)); names(dm) <- c("source", "dest", "flux")
  p <- ggplot(dm, aes(dest, source, fill = flux)) + geom_tile(colour = "white") +
    scale_fill_viridis_c(option = "rocket", direction = -1, name = "expected\nseedings") +
    labs(title = "BDBV 2026 — inter-province seeding flux (invasion corridors)",
         subtitle = "Expected number of zone invasions seeded from source province to destination province.",
         x = "destination province", y = "source province") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  .cascade_ggsave(p, file, w = 9, h = 7)
}

# ---- 8. conditional onward-spread map (generalised Kisangani) --------------
cascade_map_conditional <- function(cond_reach, hub_label, horizon = 13L,
                                    file = "cascade_conditional_reach") {
  d <- cond_reach[cond_reach$horizon == horizon, ]
  sp <- .cascade_join(d, "p_case_invasion"); if (is.null(sp)) return(invisible(NULL))
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = "grey75", linewidth = 0.05) +
    scale_fill_viridis_c(option = "inferno", na.value = "grey92", direction = -1,
                         name = sprintf("P(invaded\nby %dw)", horizon), limits = c(0, NA)) +
    labs(title = sprintf("BDBV 2026 — downstream invasion risk GIVEN %s is invaded", hub_label),
         subtitle = paste0("Designed-seeding conditional experiment. ", CASCADE_FRAMING)) +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7))
  .cascade_ggsave(p, paste0(file, "_", gsub("[^A-Za-z0-9]+", "_", hub_label), "_h", horizon))
}

# ---- 9. establishment vs reach ---------------------------------------------
cascade_fig_establishment <- function(rs, horizon = 13L,
                                      file = "cascade_establishment_vs_reach") {
  d <- rs[rs$horizon == horizon & !rs$was_active_before, ]
  d <- d[is.finite(d$p_case_invasion) & is.finite(d$p_establishment), ]
  p <- ggplot(d, aes(p_case_invasion, p_establishment)) +
    geom_abline(slope = 1, linetype = 2, colour = "grey60") +
    geom_point(aes(colour = V), alpha = 0.7) +
    scale_colour_viridis_c(option = "viridis", name = "Vulnerability V") +
    labs(title = "BDBV 2026 — invasion vs establishment (13 weeks)",
         subtitle = "Reach = P(first case); establishment = P(>=5 cases). Gap = stochastic fade-out.",
         x = "P(reached)", y = "P(established, >=5 cases)") +
    theme_minimal(base_size = 12)
  .cascade_ggsave(p, paste0(file, "_h", horizon), w = 8, h = 6.5)
}

# ---- optional: reuse the existing rank/uncertainty/choropleth helpers ------
# Wrapped so a signature mismatch degrades gracefully rather than aborting the run.
cascade_render_reuse_maps <- function(rs, file_prefix = "cascade_") {
  sp <- .cascade_shape()
  for (H in CASCADE_REPORT_HORIZONS) {
    tryCatch(plot_invasion_rank_map(rs, horizon = H, method_label = "Cascade 3-mo",
             shapefile = sp, save = TRUE, file = paste0(file_prefix, "invasion_rank_map"),
             extent = "national"), error = function(e) message("[cascade] rank map h", H, ": ", conditionMessage(e)))
    tryCatch(plot_prob_vuln_choropleth(rs, shapefile = sp, horizon = H, province_zoom = NULL,
             save = TRUE, file = paste0(file_prefix, "prob_vuln_choropleth"),
             model_label = "Cascade 3-mo"),
             error = function(e) message("[cascade] choropleth h", H, ": ", conditionMessage(e)))
  }
}
