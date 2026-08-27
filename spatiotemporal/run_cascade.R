# =============================================================================
# run_cascade.R — 3-MONTH SPATIAL INVASION CASCADE: master driver
# BDBV 2026 DRC · runs modules 30–36 of PLAN_3MONTH_INVASION.md
#
# Usage:  Rscript spatiotemporal/run_cascade.R
# Env flags:
#   CASCADE_N_MC          MC iterations (default 1000; smoke 40)
#   CASCADE_SMOKE=1       fast path (M=40, heavy stages off)
#   CASCADE_LAYER_CACHE   optional rds path to cache/reuse the reconstructed layer
#   CASCADE_RUN_BACKTEST / CASCADE_RUN_SBC / CASCADE_RUN_SENS  (default 1; smoke 0)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(here) })
ST_DIR <- Sys.getenv("CASCADE_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
t0 <- Sys.time(); tick <- function(m) message(sprintf("[cascade t=%.0fs] %s",
                                     as.numeric(Sys.time() - t0, units = "secs"), m))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

source(file.path(ST_DIR, "00_config.R"))
for (f in c("01_data_prep.R","02_epi_params.R","03_mobility_matrices.R","04_nowcasting.R",
            "04b_epinowcast.R","05_baseline_models.R","06_simple_models.R","07_hhh4_model.R",
            "08_stochastic_seir.R","15_workhorse.R","16_invasion_eval.R","18_ensemble.R",
            "19_spacetime_eval.R","20_forecast_detail.R","21_bayesian_renewal.R",
            "17_invasion_viz.R",
            "30_projection_config.R","31_source_dynamics.R","32_cascade_simulator.R",
            "33_cascade_eval.R","34_conditional_queries.R","35_cascade_viz.R","36_report.R",
            "38_urban_scenarios.R","39_cascade_eval_figure.R"))
  source(file.path(ST_DIR, f))
SMOKE <- identical(Sys.getenv("CASCADE_SMOKE"), "1")
RUN_BT   <- !SMOKE && !identical(Sys.getenv("CASCADE_RUN_BACKTEST"), "0")
RUN_SBC  <- !SMOKE && !identical(Sys.getenv("CASCADE_RUN_SBC"), "0")
RUN_SENS <- !SMOKE && !identical(Sys.getenv("CASCADE_RUN_SENS"), "0")
RUN_URBAN <- !identical(Sys.getenv("CASCADE_RUN_URBAN"), "0")
tick(sprintf("sourced (M=%d, backtest=%s sbc=%s sens=%s)", CASCADE_N_MC, RUN_BT, RUN_SBC, RUN_SENS))

# ---- reconstruct modelling layer (cache option) ----------------------------
CACHE <- Sys.getenv("CASCADE_LAYER_CACHE")
if (nzchar(CACHE) && file.exists(CACHE)) {
  layer <- readRDS(CACHE); tick("layer loaded from cache")
} else {
  dat <- prep_all_data()
  gt_pmfs <- compute_all_gt_pmfs(); osrm_mat <- load_osrm()
  zones_all <- names(dat$pop)
  zw <- nowcast_zone_week_epinowcast(zone_week = dat$zone_week, linelist = dat$ll,
          analysis_date = ANALYSIS_DATE, outbreak_start = OUTBREAK_START)
  mob_dir <- file.path(OUT_MOBILITY)
  need <- unique(c(CASCADE_KERNEL, CASCADE_KERNEL_FALLBACK, CASCADE_SENS_KERNELS))
  mobility_matrices <- setNames(lapply(need, function(k)
    readRDS(file.path(mob_dir, sprintf("mobility_%s.rds", k)))), need)
  layer <- list(zone_week_nc = zw, zones_all = zones_all, covariates = dat$covariates,
                gt_pmfs = gt_pmfs, osrm_mat = osrm_mat, mobility_matrices = mobility_matrices,
                pop = dat$pop)
  if (nzchar(CACHE)) saveRDS(layer, CACHE)
  tick("layer reconstructed")
}
zones_all <- layer$zones_all; gt_pmfs <- layer$gt_pmfs
province_map <- load_province_map()
zone_province <- setNames(province_map$province, province_map$nom)
vuln <- compute_vulnerability_index(layer$covariates, zones_all, layer$osrm_mat)

# ---- fit covariate-bearing hazards (featured + sensitivity kernels) --------
kernels <- unique(c(CASCADE_KERNEL, if (RUN_SENS) CASCADE_SENS_KERNELS))
kernels <- kernels[kernels %in% names(layer$mobility_matrices)]
iter_fit <- if (SMOKE) 800L else 2000L
designs <- list(); fits <- list()
for (kern in kernels) {
  designs[[kern]] <- build_invasion_design(layer$zone_week_nc, layer$mobility_matrices, gt_pmfs,
                       layer$covariates, layer$osrm_mat, zones_all, mob = kern, gt = CASCADE_GT)
  fits[[kern]] <- fit_bayes_renewal(designs[[kern]], cov_spec = CASCADE_COV_SPEC,
                                    iter = iter_fit, chains = 2L)
}
tick(sprintf("fitted %d kernel hazard(s): %s", length(fits), paste(names(fits), collapse=",")))

reff  <- estimate_zone_reff(layer$zone_week_nc, zones_all, gt_pmfs, gt = CASCADE_GT,
                            zone_province = zone_province)
delta <- cascade_fit_delta(cov_model = paste0("Bayes-", CASCADE_KERNEL, "-geo"))
prep  <- cascade_prepare(layer, fits[[CASCADE_KERNEL]], designs[[CASCADE_KERNEL]], reff,
                         kernel = CASCADE_KERNEL, gt = CASCADE_GT)
psi   <- if (SMOKE) CASCADE_PSI else cascade_fit_psi(prep, delta, layer$zone_week_nc, zones_all)$psi
tick(sprintf("reff R_nat=%.2f  delta=%.3f (cal=%.2f)  psi=%.2f", reff$R_nat, delta,
             attr(delta,"cal_in_large") %||% NA, psi))

# ---- consistency gate ------------------------------------------------------
gate <- cascade_consistency_gate(prep, fits[[CASCADE_KERNEL]], designs[[CASCADE_KERNEL]], delta,
          layer$zone_week_nc, layer$mobility_matrices, gt_pmfs, layer$covariates,
          layer$osrm_mat, zones_all, n_mc = if (SMOKE) 1000L else 3000L)
tick(paste("consistency gate:", gate$note))

# ---- run all scenarios -----------------------------------------------------
sims <- list(); reaches <- list(); rs_by_scenario <- list()
for (sk in names(CASCADE_SCENARIOS)) {
  sims[[sk]] <- simulate_cascade(prep, CASCADE_SCENARIOS[[sk]], n_mc = CASCADE_N_MC,
                                 delta = delta, psi = psi, pop_vec = layer$pop)
  reaches[[sk]] <- cascade_reach_table(sims[[sk]])
  rs_by_scenario[[sk]] <- cascade_enrich_reach(reaches[[sk]], province_map, vuln,
                            method_label = paste0("Cascade-", CASCADE_KERNEL, "-", sk))
}
scP_key <- CASCADE_SCENARIO_PRIMARY
rs_primary <- rs_by_scenario[[scP_key]]; base_sim <- sims[[scP_key]]
tick(sprintf("simulated %d scenarios x M=%d", length(sims), CASCADE_N_MC))

# ---- aggregates for the report ---------------------------------------------
prov_of <- function(z) zone_province[match(z, names(zone_province))]
affected_provinces0 <- unique(prov_of(zones_all[base_sim$affected0]))
agg <- lapply(names(sims), function(sk) {
  sim <- sims[[sk]]; cum <- apply(sim$new_by_week, 2, cumsum)
  # new provinces reached by 13w (provinces with no affected zone at t0)
  newprov <- vapply(seq_len(sim$n_mc), function(m) {
    zr <- zones_all[!is.na(sim$tau[, m]) & sim$tau[, m] <= 13L]
    length(setdiff(unique(prov_of(zr)), affected_provinces0)) }, numeric(1))
  # ONE consistent estimator throughout: the MEDIAN of the cumulative-new-zone
  # distribution at each horizon (+ 90% predictive band), so 4/8/13w are comparable
  # (means and medians differ for the right-skewed count; do not mix them).
  q <- function(h, p) stats::quantile(cum[h, ], p, names = FALSE)
  list(label = CASCADE_SCENARIOS[[sk]]$label,
       new4 = stats::median(cum[4, ]), new8 = stats::median(cum[8, ]),
       new13_med = stats::median(cum[13, ]), new13_lo = q(13, .05), new13_hi = q(13, .95),
       new_prov13 = mean(newprov))
})
names(agg) <- names(sims)

# ---- outputs: CSVs, maps, figures ------------------------------------------
TAB <- file.path(OUT_CASCADE, "tables")
cascade_write_reach_csv(rs_primary, file.path(TAB, "cascade_reach_scores_all_zones.csv"))
for (sk in names(rs_by_scenario))
  cascade_write_reach_csv(rs_by_scenario[[sk]], file.path(TAB, sprintf("cascade_reach_%s.csv", sk)))
for (H in CASCADE_REPORT_HORIZONS) {
  cascade_map_reach(rs_primary, H, scenario_label = CASCADE_SCENARIOS[[scP_key]]$label)
  cascade_map_uncertainty(rs_primary, H)
}
cascade_map_timetoinvasion(rs_primary)
fan <- cascade_fig_fanchart(sims)
cascade_map_scenarios(rs_by_scenario, horizon = max(CASCADE_REPORT_HORIZONS))
cascade_fig_establishment(rs_primary)
cascade_render_reuse_maps(rs_primary)
tick("maps + figures written")

# ---- onward propagation: gateway, corridors, tree, conditional -------------
Ync <- .count_wide(layer$zone_week_nc, zones_all, "confirmed_nc")
recent_inc <- rowSums(Ync[, tail(seq_len(ncol(Ync)), 4), drop = FALSE])
cand_hubs <- names(sort(recent_inc[base_sim$affected0], decreasing = TRUE))
cand_hubs <- head(cand_hubs[recent_inc[cand_hubs] > 0], CASCADE_KNOCKOUT_TOPN)
gateway <- cascade_knockout(prep, CASCADE_SCENARIOS[[scP_key]], base_sim, delta, psi = psi,
                            n_mc = min(CASCADE_N_MC, if (SMOKE) 40L else 600L), candidates = cand_hubs)
readr::write_csv(gateway, file.path(TAB, "cascade_gateway_knockout.csv"))
cascade_fig_gateway(gateway)
flux_net <- cascade_flux_network(base_sim, zone_province)
readr::write_csv(flux_net$all_edges, file.path(TAB, "cascade_flux_edges.csv"))
cascade_fig_corridors(flux_net)
readr::write_csv(cascade_transmission_tree(base_sim), file.path(TAB, "cascade_transmission_tree.csv"))
tick("gateway + corridors + tree written")

conditional <- list()
for (hub in head(gateway$health_zone, 2L)) {
  cs <- cascade_conditional_seed(prep, CASCADE_SCENARIOS[[scP_key]], hub, delta, psi = psi,
                                 n_mc = min(CASCADE_N_MC, if (SMOKE) 40L else 800L))
  rs_c <- cascade_enrich_reach(cs$reach, province_map, vuln, method_label = paste0("cond-", hub))
  cascade_map_conditional(rs_c, hub)
  r13 <- cs$reach[cs$reach$horizon == 13L & !cs$reach$was_active_before, ]
  conditional[[hub]] <- list(top = head(r13$health_zone[order(-r13$p_invasion)], 6))
  readr::write_csv(rs_c, file.path(TAB, sprintf("cascade_conditional_%s.csv",
                                    gsub("[^A-Za-z0-9]+","_",hub))))
}
tick("conditional hubs written")

# ---- urban-hub invasion scenarios (timing-aware) ---------------------------
urban <- NULL
if (RUN_URBAN) {
  urban <- tryCatch({
    hubs_u <- resolve_urban_hubs(URBAN_HUBS, zones_all)
    sw_u   <- if (SMOKE) c(1L, 4L) else URBAN_SEED_WEEKS
    # dedicated baseline at URBAN_N_MC (matched to the urban conditional runs, so the
    # per-zone attributable difference for far zones is pure low-M noise, not a mix)
    reach_base_u <- cascade_reach_table(simulate_cascade(prep, CASCADE_SCENARIOS[[scP_key]],
                      n_mc = URBAN_N_MC, delta = delta, psi = psi),
                      horizons = CASCADE_REPORT_HORIZONS)
    grid_u <- urban_scenario_grid(prep, CASCADE_SCENARIOS[[scP_key]], hubs_u, sw_u,
                delta, psi, n_mc = URBAN_N_MC, reach_base = reach_base_u,
                province_map = province_map)
    urban_write_tables(grid_u, TAB)
    for (city in names(grid_u$detail)) {
      im <- grid_u$detail[[city]]$impact
      urban_map_attributable(im, city, grid_u$min_week)
      urban_fig_top_zones(im, city)
      urban_figure4(im, city, hubs_u[[city]])          # publication two-panel (map + dumbbell)
    }
    urban_fig_summary_bar(grid_u$summary, grid_u$min_week)
    urban_fig_timing(grid_u$summary)
    urban_map_facet(setNames(lapply(names(grid_u$detail),
                    function(c) grid_u$detail[[c]]$impact$per_zone), names(grid_u$detail)))
    grid_u
  }, error = function(e) { message("[cascade] urban scenarios: ", conditionMessage(e)); NULL })
  tick(sprintf("urban scenarios written (%d hubs x %d seed-weeks)",
               length(if (!is.null(urban)) urban$detail else 0L),
               length(if (!is.null(urban)) urban$seed_weeks else 0L)))
}

# ---- validation: sensitivity, backtest, SBC --------------------------------
sens <- if (RUN_SENS) tryCatch(cascade_sensitivity(layer, designs, fits, reff, delta,
                        base_reach = reaches[[scP_key]]),
                        error = function(e) { message("[cascade] sensitivity: ", conditionMessage(e)); NULL }) else NULL
if (!is.null(sens)) readr::write_csv(sens, file.path(OUT_CASCADE, "diagnostics", "cascade_sensitivity.csv"))
backtest <- if (RUN_BT) tryCatch(cascade_backtest(layer, delta = delta, psi = psi,
                        zone_province = zone_province),
                        error = function(e) { message("[cascade] backtest: ", conditionMessage(e)); NULL }) else NULL
if (!is.null(backtest)) {
  readr::write_csv(backtest, file.path(OUT_CASCADE, "diagnostics", "cascade_backtest.csv"))
  # per-zone detail (drives the Figure2_cascade capture-curve / mean-rank / top-K panels)
  bt_detail <- attr(backtest, "detail")
  if (!is.null(bt_detail)) readr::write_csv(bt_detail,
    file.path(OUT_CASCADE, "diagnostics", "cascade_backtest_detail.csv"))
}
sbc <- if (RUN_SBC) tryCatch(cascade_sbc(designs[[CASCADE_KERNEL]], n_sbc = 20L),
                        error = function(e) { message("[cascade] sbc: ", conditionMessage(e)); NULL }) else NULL
tick("validation done")

# ---- report + key_outputs bundle -------------------------------------------
ctx <- list(reach_primary = rs_primary, n_mc = CASCADE_N_MC, delta = delta, psi = psi,
            gate = gate, agg = agg, gateway = gateway, conditional = conditional,
            sens = sens, backtest = backtest, sbc = sbc, urban = urban)
report_path <- file.path(ST_DIR, "SPATIAL_INVASION_3MONTH_REPORT.md")
cascade_write_report(report_path, ctx)

# evaluation figure (Figure 2 style): leave-future-out backtest of the cascade —
# discrimination (A), prioritisation capture curve (B), ranking accuracy (C),
# per-round forecasts vs realised outcomes (D). Needs ctx$backtest + its detail attr.
tryCatch(cascade_fig_evaluation(ctx),
         error = function(e) message("[cascade] evaluation figure skipped: ", conditionMessage(e)))

# key_outputs bundle (mirror the short-horizon key_outputs style)
KO <- file.path(OUT_DIR, "key_outputs")
if (!dir.exists(KO)) dir.create(KO, recursive = TRUE)
figdir <- file.path(OUT_CASCADE, "figures")
ko_files <- c(file.path(TAB, "cascade_reach_scores_all_zones.csv"),
              file.path(TAB, "cascade_gateway_knockout.csv"),
              file.path(TAB, "urban_scenario_summary.csv"),
              list.files(figdir, pattern = "^cascade_(reach_map_national_h13|time_to_invasion|newzones_fanchart|scenario_reach_maps_h13|gateway_knockout|corridors)", full.names = TRUE),
              list.files(figdir, pattern = "^(urban_impact_summary|urban_timing_sensitivity)", full.names = TRUE))
for (f in ko_files) if (file.exists(f)) file.copy(f, file.path(KO, basename(f)), overwrite = TRUE)

# publication-style Figure 4 (cascade): relative-risk map + top-20 ranked reach,
# written to key_outputs/figures/Figure4_cascade_h{4,8,13}.{pdf,png} (+ data CSVs).
tryCatch(source(file.path(ST_DIR, "37_cascade_figure4.R")),
         error = function(e) message("[cascade] Figure4_cascade skipped: ", conditionMessage(e)))

# "next dominoes" narrative figure (main-text Figure 4): anchor frontier map + city ranking
# + per-city watch-lists, written to key_outputs/manuscript_figures/Figure4.{pdf,png} (+ CSV).
tryCatch(source(file.path(ST_DIR, "40_cascade_next_dominoes.R")),
         error = function(e) message("[cascade] Figure4 (next dominoes) skipped: ", conditionMessage(e)))
tick(sprintf("DONE — report + %d key_outputs files + Figure4_cascade + Figure4 (next dominoes)", length(ko_files)))
