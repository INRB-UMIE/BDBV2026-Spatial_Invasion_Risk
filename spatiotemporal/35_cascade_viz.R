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
# rankings/relative risk; absolute probabilities are conditional on the scenario."
#
# Sourced AFTER 30-34, 15_workhorse.R (compute_risk_scores), 20_forecast_detail.R
# (add_risk_indices, compute_vulnerability_index, map helpers) and 17_invasion_viz.R.
# Sourcing defines functions only (nothing is written), which is why the standalone readers
# 37_cascade_figure4.R, 40_cascade_next_dominoes.R and 42_kindiv_sweep_figure.R also source it,
# for the provenance helpers cascade_stamp_check() / cascade_calibration_note().
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

# NOT "UPPER BOUNDS". That claim was carried on every product with nothing to support it: the
# level is calibrated (delta is fitted to realised K-week counts at held-out origins), and the
# status-quo scenario holds transmissibility at its current value rather than at a ceiling, so
# no argument makes the absolute probabilities one-sided. What IS true, and is said instead, is
# that they are conditional — on the scenario's c(t) and on a level calibrated over K weeks
# rather than over the 13 the projection runs.
CASCADE_FRAMING <- paste0("3-month PROJECTION (scenario-based). Trust rankings & ",
  "relative risk; absolute probabilities are conditional on the scenario. Time index = symptom onset.")

# ---- enrichment + CSV ------------------------------------------------------
cascade_enrich_reach <- function(reach, province_map, vuln, method_label) {
  rs <- compute_risk_scores(reach, province_map)
  rs <- add_risk_indices(rs, vuln)
  rs$method <- method_label
  rs
}

#' Write a reach table, STAMPED with the run that produced it.
#'
#' The provenance columns (cascade_run_id, analysis_date, calibrated, calibration_status) are
#' what 37_cascade_figure4.R and 40_cascade_next_dominoes.R check against
#' cascade_run_stamp.json before drawing Figure 4: without them a table left by an earlier, or
#' an uncalibrated, run was redrawn as current with nothing to say so. The defaults read
#' run_cascade.R's globals (CASCADE_RUN_ID, CALIBRATED, CALIB_STATUS); outside that driver they
#' are NA, which every reader refuses or labels as unrecorded.
cascade_write_reach_csv <- function(rs, path,
                                    run_id = get0("CASCADE_RUN_ID", ifnotfound = NA_character_),
                                    calibrated = get0("CALIBRATED", ifnotfound = NA),
                                    calibration_status = get0("CALIB_STATUS",
                                                              ifnotfound = NA_character_)) {
  keep <- c("health_zone", "province", "horizon", "was_active_before",
            # p_sd is the posterior SD of the reach probability; it is computed in
            # cascade_reach_table() but was never exported, so the honest posterior width never
            # reached a consumer and p_lo/p_hi were the only signal of uncertainty.
            "p_case_invasion", "p_lo", "p_hi", "p_ci_at_mc_floor", "p_sd", "p_establishment",
            "first_passage_week",
            "mu_forecast", "rr_nat", "rr_nat_rank", "rr01_nat",
            "rr_ituri", "rr_ituri_rank", "rr_nordkivu", "rr_nordkivu_rank",
            "rr_hautuele", "rr_hautuele_rank", "V", "priority", "priority_rank",
            "surveillance_gap", "healthcare_gap", "access_gap", "social_vulnerability",
            "healthcare_travel_min", "method")
  # SIGNIFICANT DIGITS, not 4 decimal places — the identical defect that was fixed in
  # .write_risk_csv() (run_all.R) and must not survive in the cascade's twin. round(x, 4)
  # publishes every probability below 5e-5 as EXACTLY 0, and a renewal hazard gives
  # p = 1 - exp(-mu) > 0 strictly, so each of those zeros is an artefact of the writer. In
  # this file NA already means "not applicable" (an affected zone), so a printed 0 is read as
  # a real, vanishing risk.
  out <- rs |> dplyr::select(dplyr::any_of(keep)) |>
    dplyr::rename(dplyr::any_of(c(p_case_lo = "p_lo", p_case_hi = "p_hi"))) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 6)))
  .a1 <- function(x, cast) if (length(x) == 1L) cast(x) else cast(NA)
  out$cascade_run_id     <- .a1(run_id, as.character)
  out$analysis_date      <- format(as.Date(get0("ANALYSIS_DATE", ifnotfound = NA)))
  out$calibrated         <- .a1(calibrated, as.logical)
  out$calibration_status <- .a1(calibration_status, as.character)
  if (all(c("horizon", "p_case_invasion") %in% names(out)))
    out <- out |> dplyr::arrange(horizon, dplyr::desc(dplyr::coalesce(p_case_invasion, -1)))
  readr::write_csv(out, path)
  message(sprintf("[cascade] wrote %d rows -> %s", nrow(out), basename(path)))
  invisible(out)
}

# ---- provenance: is a cascade product THIS run's? ---------------------------
#' The id of the cascade run in progress in, or launching, this process ("" if none).
#'
#' run_cascade.R defines CASCADE_RUN_ID at top level, so a figure script it source()s sees it as
#' a global, and it passes it to a figure script it launches as a child process in the
#' environment variable of the same name. That run's stamp still reads "running" while it draws
#' its own figures, so this is the one case in which a "running" stamp is accepted.
cascade_live_run_id <- function() {
  v <- Sys.getenv("CASCADE_RUN_ID", "")
  if (nzchar(v)) return(v)
  g <- get0("CASCADE_RUN_ID", envir = globalenv(), inherits = FALSE, ifnotfound = "")
  if (is.character(g) && length(g) == 1L && !is.na(g)) g else ""
}

#' Do products carrying `run_ids` come from the LATEST cascade run for this analysis date?
#'
#' Checked against diagnostics/cascade_run_stamp.json (run_cascade.R writes it "running" at the
#' start, "failed: ..." when calibration fails, "complete" at the end). Accepted only when the
#' stamp's analysis_date is `analysis_date`, its status is "complete" — or "running" for the run
#' in progress here (cascade_live_run_id()) — and every one of `run_ids` is the stamp's
#' cascade_run_id. Anything else (no stamp, a failed or abandoned run, another date, a product an
#' earlier run left behind, an unstamped product) is refused: the same rule
#' make_manuscript_figure2_cascade.R applies to the backtest.
#' @param run_ids the cascade_run_id value(s) the product carries (NULL when it has none).
#' @param diag_dir the cascade diagnostics directory holding the stamp.
#' @return list(ok, why, stamp): `why` says what failed (NA when ok).
cascade_stamp_check <- function(run_ids, diag_dir,
                                analysis_date = get0("ANALYSIS_DATE", ifnotfound = NA),
                                live_run_id = cascade_live_run_id()) {
  f <- file.path(diag_dir, "cascade_run_stamp.json")
  stamp <- if (file.exists(f)) tryCatch(jsonlite::fromJSON(f), error = function(e) NULL) else NULL
  .no <- function(why) list(ok = FALSE, why = why, stamp = stamp)
  if (is.null(stamp)) return(.no(sprintf("no readable %s in %s", basename(f), diag_dir)))
  .s1 <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) NA_character_ else as.character(x)
  sid <- .s1(stamp$cascade_run_id); sst <- .s1(stamp$status); sdate <- .s1(stamp$analysis_date)
  want <- if (length(analysis_date) == 1L && !is.na(analysis_date))
    format(as.Date(analysis_date)) else NA_character_
  if (is.na(sid)) return(.no("the run stamp carries no cascade_run_id"))
  if (is.na(want) || !identical(sdate, want))
    return(.no(sprintf("the latest cascade run (%s) is for analysis date %s, not %s",
                       sid, sdate, want)))
  live <- identical(sst, "running") && nzchar(live_run_id) && identical(live_run_id, sid)
  if (!identical(sst, "complete") && !live)
    return(.no(sprintf("the latest cascade run (%s) did not complete (status: %s), so its outputs are partial",
                       sid, sst)))
  ids <- unique(as.character(run_ids))
  if (!length(ids) || anyNA(ids) || !all(ids == sid))
    return(.no(sprintf("the product was not written by the latest cascade run %s (it carries %s)",
                       sid, if (length(ids)) paste(ids, collapse = ", ") else "no run id")))
  list(ok = TRUE, why = NA_character_, stamp = stamp)
}

#' The note an UNCALIBRATED cascade product carries; NULL when the run was calibrated.
#'
#' run_cascade.R sets CALIBRATED = FALSE only when calibration was disabled deliberately
#' (CASCADE_CALIBRATE=0 / CASCADE_SMOKE=1; a FAILED calibration stops the run instead), and its
#' header promises that every product then says so. ASCII only: figure_pdf_device() encodes an
#' em-dash, but its last-resort base pdf device (ISOLatin1) cannot.
#' @param calibrated the product's `calibrated` flag(s); NA or mixed = not recorded.
#' @param status the run's calibration_status text, quoted when recorded.
cascade_calibration_note <- function(calibrated, status = NA_character_) {
  cal <- unique(as.logical(calibrated))
  if (length(cal) == 1L && isTRUE(cal)) return(NULL)
  st <- unique(as.character(status)); st <- st[!is.na(st) & nzchar(st)]
  if (length(cal) == 1L && isFALSE(cal))
    return(paste0("UNCALIBRATED cascade run", if (length(st) == 1L) paste0(" (", st, ")"),
                  ": its probability level was not fitted to realised multi-week counts, so ",
                  "absolute probabilities and expected counts are not calibrated at the ",
                  "horizons shown."))
  "Calibration status of this cascade run is not recorded: absolute probabilities and expected counts may be uncalibrated."
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
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  # Gate on the FULL destination path, not the bare stem: FIGURE_DROP entries are
  # "<directory>/<stem>" and the raw/ exclusion inspects path components, neither of
  # which a basename can match.
  path <- file.path(OUT_CASCADE, "figures", file)
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !.fk(path)) return(invisible(p))
  ggsave(paste0(path, ".pdf"), p, width = w, height = h, device = figure_pdf_device())
  ggsave(paste0(path, ".png"), p, width = w, height = h, dpi = 150)
  invisible(path)
}
.PROV_LIMS <- get0("PROVINCES_OF_INTEREST", ifnotfound = c("Ituri","Nord-Kivu","Haut-Uele"))

# ---- 1. cumulative-reach choropleth (national) -----------------------------
cascade_map_reach <- function(rs, horizon = CASCADE_HORIZON_WEEKS, file = "cascade_reach_map_national",
                              scenario_label = "") {
  d <- rs[rs$horizon == horizon, ]
  sp <- .cascade_join(d, "p_case_invasion"); if (is.null(sp)) return(invisible(NULL))
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = "grey70", linewidth = 0.05) +
    scale_fill_viridis_c(option = "inferno", na.value = "grey92", direction = -1,
                         name = sprintf("P(invaded\nby %d wk)", horizon), limits = c(0, NA)) +
    # ASCII hyphens in every rendered title here: the last-resort base pdf device of
    # figure_pdf_device() cannot encode an em-dash (see cascade_fig_corridors()).
    labs(title = sprintf("BDBV 2026 - %d-week spatial invasion reach%s", horizon,
                         if (nzchar(scenario_label)) paste0(" [", scenario_label, "]") else ""),
         subtitle = CASCADE_FRAMING, caption = "Grey = already-affected (masked) or no data.") +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7))
  .cascade_ggsave(p, paste0(file, "_h", horizon))
}

# ---- 2. uncertainty map (credible-interval width) --------------------------
cascade_map_uncertainty <- function(rs, horizon = CASCADE_HORIZON_WEEKS,
                                    file = "cascade_reach_uncertainty_national") {
  d <- rs[rs$horizon == horizon, ]; d$.w <- d$p_hi - d$p_lo
  # A Monte-Carlo-floor interval (p_ci_at_mc_floor, cascade_reach_table()) is not an estimated
  # credible interval: its width is the 1/D resolution of D parameter draws by construction.
  # Shading it as a CrI width would map an artefact of D onto the fill scale, so it is left
  # unshaded and the subtitle says so.
  if ("p_ci_at_mc_floor" %in% names(d)) d$.w[d$p_ci_at_mc_floor %in% TRUE] <- NA_real_
  sp <- .cascade_join(d, ".w"); if (is.null(sp)) return(invisible(NULL))
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = "grey70", linewidth = 0.05) +
    scale_fill_viridis_c(option = "mako", na.value = "grey92",
                         name = "90% CrI\nwidth", direction = -1) +
    labs(title = sprintf("BDBV 2026 - %d-week reach: forecast uncertainty", horizon),
         subtitle = paste("Width of the 90% credible interval on the reach probability (posterior parameter uncertainty; excludes between-scenario/structural uncertainty).",
                          "Grey = already affected, no data, or an interval below the Monte-Carlo resolution (1/D).")) +
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
    labs(title = "BDBV 2026 - spatial invasion: the expanding frontier",
         subtitle = CASCADE_FRAMING) +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7))
  .cascade_ggsave(p, file, w = 13, h = 5.5)
}

# ---- 4. aggregate fan chart (new zones over time, per scenario) ------------
#' @return the per-week median and 90% band that the ribbon draws (one row per scenario x week).
#'   The function has always returned it; run_cascade.R assigned the result to `fan` and then
#'   never used it, so the 90% credible band on a PUBLISHED figure existed only as pixels and no
#'   reader could recover a number from it. The caller now writes it.
#'
#'   THIS IS WHERE A TOTAL'S INTERVAL COMES FROM. The band is the 5th/95th percentile of the
#'   SIMULATED NATIONAL TOTAL across Monte-Carlo iterations (a 90% predictive interval for the
#'   count). It must never be rebuilt by summing the per-zone p_lo/p_hi of the reach table:
#'   quantiles are not additive, and the per-zone bounds are posterior intervals on
#'   probabilities, not on counts (a manuscript audit found a total's interval formed that
#'   way). The total's EXPECTATION, by contrast, is additive: it is the sum of the per-zone p,
#'   and `mean` carries it beside the median so the two are never confused.
cascade_fig_fanchart <- function(sims_by_scenario,
                                 file = "cascade_newzones_fanchart") {
  rows <- lapply(names(sims_by_scenario), function(sk) {
    sim <- sims_by_scenario[[sk]]
    cum <- apply(sim$new_by_week, 2, cumsum)               # weeks x iters
    tibble::tibble(week = seq_len(nrow(cum)),
      scenario = sim$knobs$scenario,
      mean = rowMeans(cum),
      med = apply(cum, 1, stats::median),
      lo  = apply(cum, 1, stats::quantile, 0.05),
      hi  = apply(cum, 1, stats::quantile, 0.95))
  })
  d <- dplyr::bind_rows(rows)
  p <- ggplot(d, aes(week, med, colour = scenario, fill = scenario)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    # MEDIAN, not "expected": the line is the median of the simulated cumulative count and the
    # ribbon its 5th-95th percentiles across iterations (the mean is in the returned data).
    labs(title = sprintf("BDBV 2026 - cumulative NEW health zones invaded over %d weeks (median)",
                         max(d$week)),
         subtitle = paste("Median + 90% predictive band across Monte-Carlo iterations, by scenario.",
                          "Projection, not a point forecast."),
         x = "weeks ahead", y = "cumulative newly-invaded zones",
         colour = "Scenario", fill = "Scenario") +
    theme_minimal(base_size = 12)
  .cascade_ggsave(p, file, w = 9, h = 6)
  d
}

# ---- 5. scenario comparison reach maps -------------------------------------
cascade_map_scenarios <- function(rs_by_scenario, horizon = CASCADE_HORIZON_WEEKS,
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
    labs(title = sprintf("BDBV 2026 - %d-week reach under control scenarios", horizon),
         subtitle = CASCADE_FRAMING) +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7))
  .cascade_ggsave(p, paste0(file, "_h", horizon), w = 14, h = 5.5)
}

# ---- 6. gateway (knockout Delta_j) bars ------------------------------------
cascade_fig_gateway <- function(gateway_tbl, file = "cascade_gateway_knockout",
                                top_n = 15L) {
  d <- head(gateway_tbl, top_n)
  d$health_zone <- factor(d$health_zone, levels = rev(d$health_zone))
  has_ci <- all(c("delta_j_lo", "delta_j_hi") %in% names(d)) && any(is.finite(d$delta_j_lo))
  p <- ggplot(d, aes(delta_j, health_zone)) +
    geom_col(fill = "#b2182b") +
    # The paired 90% interval. Drawn because the ranking of the trailing zones is not
    # separated by this analysis, and a bare bar chart invites reading an order into it.
    {if (has_ci) geom_linerange(aes(xmin = delta_j_lo, xmax = delta_j_hi),
                                colour = "grey25", linewidth = 0.5)} +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
    # The BAR is the point estimate; the LINE is the interval (cascade_paired_total(): the
    # paired 90% Monte-Carlo interval on the expectation, mean +- z * SE over parameter-draw
    # groups), not a spread of outcomes.
    labs(title = "BDBV 2026 - gateway zones: onward-spread impact",
         subtitle = paste("Delta_j = expected downstream invasions prevented if this source is",
                          "contained (knockout).\nPaired contrast against a same-seed baseline;",
                          "bar = estimate, line = 90% Monte-Carlo interval on the expected value.",
                          "Zones whose interval spans zero are not separated by this analysis."),
         x = "expected downstream invasions attributable (Delta_j)", y = NULL) +
    theme_minimal(base_size = 12)
  .cascade_ggsave(p, file, w = 9, h = 6)
}

# ---- 7. seeding corridors (province flux heat + top zone edges) ------------
#' @param horizon,scenario_label what `flux_net` covers. cascade_flux_network() keeps neither,
#'   so they are formals; the defaults are run_cascade.R's call, which builds `flux_net` from the
#'   PRIMARY scenario's simulation over the full projection horizon (simulate_cascade()'s default).
cascade_fig_corridors <- function(flux_net, file = "cascade_corridors",
                                  horizon = CASCADE_HORIZON_WEEKS,
                                  scenario_label = CASCADE_SCENARIOS[[CASCADE_SCENARIO_PRIMARY]]$label) {
  pm <- flux_net$province_matrix
  if (is.null(pm)) return(invisible(NULL))
  dm <- as.data.frame(as.table(pm)); names(dm) <- c("source", "dest", "flux")
  p <- ggplot(dm, aes(dest, source, fill = flux)) + geom_tile(colour = "white") +
    scale_fill_viridis_c(option = "rocket", direction = -1, name = "expected\nseedings") +
    # ASCII hyphen, not an em-dash: this title is rendered, and if figure_pdf_device() falls back
    # to the base pdf device (ISOLatin1) it cannot encode U+2014 and substitutes a hyphen with a
    # conversion warning — so the PDF and the 150-dpi PNG could disagree.
    # cascade_corridors is a retained deliverable, so its title is not suppressed.
    labs(title = "BDBV 2026 - inter-province seeding flux (invasion corridors)",
         # sim$flux counts FIRST invasions only (a re-seeding is not attributed a source), per
         # Monte-Carlo iteration, so the tiles are expected first-invasion seedings.
         subtitle = sprintf(paste0("Expected number of first invasions of zones in the destination province ",
                                   "seeded from the source province, over the %d-week projection ",
                                   "(scenario: %s)."), as.integer(horizon), scenario_label),
         x = "destination province", y = "source province") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  .cascade_ggsave(p, file, w = 9, h = 7)
}

# ---- 8. conditional onward-spread map (generalised Kisangani) --------------
cascade_map_conditional <- function(cond_reach, hub_label, horizon = CASCADE_HORIZON_WEEKS,
                                    file = "cascade_conditional_reach") {
  d <- cond_reach[cond_reach$horizon == horizon, ]
  sp <- .cascade_join(d, "p_case_invasion"); if (is.null(sp)) return(invisible(NULL))
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = "grey75", linewidth = 0.05) +
    scale_fill_viridis_c(option = "inferno", na.value = "grey92", direction = -1,
                         name = sprintf("P(invaded\nby %dw)", horizon), limits = c(0, NA)) +
    labs(title = sprintf("BDBV 2026 - downstream invasion risk GIVEN %s is invaded", hub_label),
         subtitle = paste0("Designed-seeding conditional experiment. ", CASCADE_FRAMING)) +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7))
  .cascade_ggsave(p, paste0(file, "_", gsub("[^A-Za-z0-9]+", "_", hub_label), "_h", horizon))
}

# ---- 9. establishment vs reach ---------------------------------------------
#' @param n_est the establishment threshold the reach table's p_establishment was computed at
#'   (simulate_cascade()'s n_est; run_cascade.R leaves it at CASCADE_N_EST).
cascade_fig_establishment <- function(rs, horizon = CASCADE_HORIZON_WEEKS,
                                      file = "cascade_establishment_vs_reach",
                                      n_est = CASCADE_N_EST) {
  d <- rs[rs$horizon == horizon & !rs$was_active_before, ]
  d <- d[is.finite(d$p_case_invasion) & is.finite(d$p_establishment), ]
  p <- ggplot(d, aes(p_case_invasion, p_establishment)) +
    geom_abline(slope = 1, linetype = 2, colour = "grey60") +
    geom_point(aes(colour = V), alpha = 0.7) +
    scale_colour_viridis_c(option = "viridis", name = "Vulnerability V") +
    # The gap is every iteration that reached the zone but had not accumulated n_est cases by
    # the horizon: a stochastic fade-out, OR a first case too late in the window to get there.
    labs(title = sprintf("BDBV 2026 - invasion vs establishment (%d weeks)", horizon),
         subtitle = sprintf(paste0("Reach = P(first case); establishment = P(>=%d cases by week %d). ",
                                   "Gap = reached but not established by then (fade-out, or reached too late)."),
                            as.integer(n_est), as.integer(horizon)),
         x = "P(reached)", y = sprintf("P(established, >=%d cases)", as.integer(n_est))) +
    theme_minimal(base_size = 12)
  .cascade_ggsave(p, paste0(file, "_h", horizon), w = 8, h = 6.5)
}

# ---- optional: reuse the existing rank/uncertainty/choropleth helpers ------
# Wrapped so a signature mismatch degrades gracefully rather than aborting the run.
# Those helpers write to OUT_MAPS, the PUBLISHED tree, and take no output directory. In a run
# sandboxed by CASCADE_OUT_ROOT (an onset-imputation arm; 30_projection_config.R) every other
# cascade product follows CASCADE_OUT_DIR, so these maps are skipped there rather than written
# over the published ones.
cascade_render_reuse_maps <- function(rs, file_prefix = "cascade_") {
  .cod <- get0("CASCADE_OUT_DIR", ifnotfound = OUT_DIR)
  if (!identical(normalizePath(.cod, mustWork = FALSE), normalizePath(OUT_DIR, mustWork = FALSE))) {
    message("[cascade] CASCADE_OUT_ROOT is set: the reused rank/choropleth maps (which write to ",
            "OUT_MAPS, the published tree) are skipped.")
    return(invisible(NULL))
  }
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
