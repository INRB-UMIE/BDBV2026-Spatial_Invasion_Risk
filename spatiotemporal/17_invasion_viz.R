# =============================================================================
# 17_invasion_viz.R — Legible, invasion-focused visualisation suite
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Replaces the illegible, model-cluttered figures with a curated, legible set for
# the invasion task. Design rules (dataviz): window every time axis to the signal;
# <= 8 categorical colours in fixed CVD-safe order (Okabe-Ito), rest -> "Other";
# sequential viridis for magnitude; mask already-affected zones out of every
# invasion figure; small multiples over overplotting; recessive grid; legible
# fonts (>= 9pt).
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))
suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(patchwork); library(scales); library(viridis)
})

# CVD-safe qualitative palette (Okabe-Ito), assigned in fixed order.
OKABE_ITO <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
               "#56B4E9", "#F0E442", "#000000")

theme_inv <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      plot.title    = element_text(face = "bold", size = base + 2),
      plot.subtitle = element_text(colour = "grey35", size = base - 1),
      plot.caption  = element_text(colour = "grey55", size = base - 3),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      axis.title    = element_text(size = base - 1),
      legend.position = "top", legend.title = element_text(size = base - 2),
      strip.text    = element_text(face = "bold", size = base - 1)
    )
}

.save <- function(p, path, w = 9, h = 6) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !.fk(path)) return(invisible(p))
  # The device is figure_pdf_device() (00_config.R): one that OPENED cleanly in this session and
  # can draw the figures' dashes and Greek letters; the base pdf device turns them into dots.
  # Delete any stale target first, then verify the write.
  if (file.exists(path)) suppressWarnings(file.remove(path))
  tryCatch(ggsave(path, p, width = w, height = h, device = figure_pdf_device(),
                  limitsize = FALSE),
           error = function(e) warning("[viz] ggsave error ", basename(path),
                                       ": ", conditionMessage(e)))
  if (file.exists(path)) message("[viz] saved -> ", path)
  else warning("[viz] FAILED to write ", basename(path))
  invisible(p)
}

# ---------------------------------------------------------------------------
# Probability-scale labels: every probability this suite draws or prints says its scale
# ---------------------------------------------------------------------------

#' Scale of a CURRENT-FORECAST risk frame's p_case_invasion, as a label: "recalibrated, as
#' deployed" or "raw". Resolved, in order, from an explicit `prob_scale`, the frame's own
#' `prob_scale` column (the published risk CSVs carry one, .write_risk_csv() in run_all.R), and
#' otherwise the deploy switch INVASION_RECALIBRATE_DEPLOY (00_config.R) under which run_all.R
#' builds the in-memory frame.
.risk_prob_scale <- function(risk_df, prob_scale = NULL) {
  v <- prob_scale
  if (is.null(v) && "prob_scale" %in% names(risk_df))
    v <- unique(stats::na.omit(as.character(risk_df$prob_scale)))
  if (is.null(v) || !length(v))
    v <- if (isTRUE(get0("INVASION_RECALIBRATE_DEPLOY", ifnotfound = FALSE)))
           "recalibrated-deployed" else "raw"
  if (length(v) != 1L) {
    warning("[viz] risk frame carries more than one probability scale: ",
            paste(v, collapse = ", "), call. = FALSE)
    return(paste(v, collapse = " + "))
  }
  switch(v, raw = "raw", `recalibrated-deployed` = "recalibrated, as deployed", v)
}

#' Is an LFO frame's p_invasion RAW? fs_apply_lfo_scale() (forecast_scale.R) swaps the selected
#' scale into p_invasion and keeps the original as p_invasion_raw; a frame without that column
#' (run_all.R passes lfo_results unswapped) is on the raw scale. A swapped frame whose values all
#' equal the raw ones (no round carried a fitted factor) is raw too.
.lfo_is_raw <- function(d) {
  if (!"p_invasion_raw" %in% names(d)) return(TRUE)
  ok <- is.finite(d$p_invasion) & is.finite(d$p_invasion_raw)
  !any(d$p_invasion[ok] != d$p_invasion_raw[ok])
}

# ---------------------------------------------------------------------------
# 1. Legible epidemic curve (windowed, top-N zones + Other, weekly)
# ---------------------------------------------------------------------------

#' @param zone_week_nc nowcast-corrected zone-week tibble.
#' @param province_map accepted for call compatibility and NOT used: the named series are the
#'   top zones NATIONALLY by total cases, whatever their province.
#' @param window_start first week drawn; NULL (default) = two weeks before OUTBREAK_START. No
#'   date is invented when OUTBREAK_START is absent: the call fails instead.
plot_epidemic_curve_legible <- function(zone_week_nc, province_map = NULL,
                                        top_n = 6L, window_start = NULL,
                                        save = TRUE) {
  if (is.null(window_start)) {
    .ob <- get0("OUTBREAK_START", ifnotfound = NULL)
    if (is.null(.ob)) stop("[viz] epidemic curve: OUTBREAK_START is not defined and no ",
                           "window_start was given.", call. = FALSE)
    window_start <- as.Date(.ob) - 14
  }
  d <- zone_week_nc %>%
    dplyr::filter(week_start >= window_start) %>%
    dplyr::mutate(conf = ifelse(is.na(confirmed_nc), 0, confirmed_nc))
  # rank named zones by total cases, nationally
  tot <- d %>% dplyr::group_by(health_zone) %>%
    dplyr::summarise(t = sum(conf), .groups = "drop") %>% dplyr::arrange(dplyr::desc(t))
  top <- head(tot$health_zone[tot$t > 0], top_n)
  d <- d %>% dplyr::mutate(zone = factor(ifelse(health_zone %in% top, health_zone, "Other"),
                                         levels = c(top, "Other")))
  agg <- d %>% dplyr::group_by(week_start, zone) %>%
    dplyr::summarise(cases = sum(conf), .groups = "drop")
  pal <- c(setNames(OKABE_ITO[seq_along(top)], top), Other = "grey75")

  p <- ggplot(agg, aes(week_start, cases, fill = zone)) +
    geom_col(width = 6, colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_x_date(date_breaks = "1 week", date_labels = "%d %b",
                 expand = expansion(mult = c(0.01, 0.02))) +
    labs(title = "BDBV 2026 — weekly confirmed cases by health zone",
         subtitle = sprintf("Nowcast-corrected onset weeks; top %d zones named, remainder pooled", top_n),
         x = NULL, y = "Confirmed cases",
         caption = paste0("Analysis date: ", ANALYSIS_DATE)) +
    theme_inv() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, "epidemic_curve.pdf"), w = 9, h = 5.5)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 2. Invasion risk map — Ituri zoom + national inset (affected zones masked)
# ---------------------------------------------------------------------------

#' @param risk_df forecast/risk tibble for ONE method & horizon with health_zone,
#'   p_case_invasion (or p_invasion), was_active_before.
plot_invasion_risk_map <- function(risk_df, horizon = 1L, method_label = "",
                                    shapefile = NULL, save = TRUE, window_txt = "",
                                    file = "invasion_risk_map") {
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile)) return(invisible(NULL))
  pcol <- if ("p_case_invasion" %in% names(risk_df)) "p_case_invasion" else "p_invasion"
  # Province-aware join: two zone NAMES are duplicated nationally (Bili in
  # Nord-Ubangi & Bas-Uele; Lubunga in Kasai-Central & Tshopo), so a name-only join
  # is many-to-many on the national panel and mis-colours one of each pair.
  has_prov <- "province" %in% names(risk_df)
  rd <- risk_df %>% dplyr::filter(horizon == !!horizon) %>%
    dplyr::transmute(.key = tolower(trimws(health_zone)),
                     # Normalise the province key the SAME way as the zone key: a raw-vs-normalised
                     # mismatch (case/whitespace) would silently fail the composite join for a whole
                     # province, rendering it as "no data" grey.
                     .prov = if (has_prov) tolower(trimws(as.character(province))) else NA_character_,
                     p = .data[[pcol]],
                     affected = was_active_before)
  shp <- shapefile %>% dplyr::mutate(.key = tolower(trimws(Nom)),
                                     .prov = tolower(trimws(as.character(PROVINCE))))
  m <- if (has_prov && !all(is.na(rd$.prov)))
         dplyr::left_join(shp, rd, by = c(".key", ".prov"))
       else dplyr::left_join(shp, dplyr::select(rd, -.prov), by = ".key")
  .lgd <- sprintf("P(first case)\n(%s)", .risk_prob_scale(risk_df))

  base_layer <- function(dat, title, legend = FALSE) {
    ggplot(dat) +
      geom_sf(aes(fill = p), colour = "grey80", linewidth = 0.12) +
      geom_sf(data = dplyr::filter(dat, affected %in% TRUE),
              fill = "grey55", colour = "grey80", linewidth = 0.12) +
      scale_fill_viridis_c(option = "plasma", name = .lgd,
                           limits = c(0, max(0.05, quantile(rd$p, 0.99, na.rm = TRUE))),
                           oob = scales::squish, na.value = "grey90") +
      labs(title = title) +
      theme_void(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 12),
            legend.position = if (legend) "right" else "none")
  }
  ituri_zones <- shapefile$Nom[shapefile$PROVINCE == "Ituri"]
  m_ituri <- m %>% dplyr::filter(Nom %in% ituri_zones)
  p_it  <- base_layer(m_ituri, sprintf("Ituri — invasion risk (%s, h=%dw)", method_label, horizon), legend = TRUE) +
    labs(subtitle = "Grey = zones already affected (no invasion probability)") +
    theme(plot.subtitle = element_text(size = 9, colour = "grey35"))
  p_nat <- base_layer(m, "National context") +
    theme(plot.title = element_text(size = 10))
  out <- p_it + p_nat + patchwork::plot_layout(widths = c(2, 1))
  if (nzchar(window_txt))
    out <- out + patchwork::plot_annotation(
      caption = window_txt,
      theme = ggplot2::theme(plot.caption = element_text(colour = "grey40", size = 8, hjust = 0)))
  if (save) .save(out, file.path(OUT_MAPS, sprintf("%s_h%d.pdf", file, horizon)), w = 12, h = 6.5)
  invisible(out)
}

#' Map of invasion-probability RANKINGS across zones (1 = highest risk). Ranks are more robust
#' than the absolute probabilities (which over-predict out-of-sample), so this is the operational
#' targeting view: at-risk zones are ranked nationally by expected introductions (the published
#' rr_nat_rank, posterior-mean mu; equivalently P(first case) up to the ordering of posterior means),
#' coloured by rank (bright = top), and the top-K are labelled with their rank number.
#' Already-affected zones are unranked (grey).
#' @param show_title FALSE drops the map title and the explanatory subtitle. The published
#'   bayes_invasion_rank_map_national panels carry their caption in the manuscript text.
plot_invasion_rank_map <- function(risk_df, horizon = 1L, method_label = "",
                                   shapefile = NULL, top_k = 12L, save = TRUE,
                                   window_txt = "", file = "invasion_rank_map",
                                   extent = c("both", "national", "ituri"),
                                   show_title = TRUE) {
  extent <- match.arg(extent)
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile)) return(invisible(NULL))
  pcol <- if ("p_case_invasion" %in% names(risk_df)) "p_case_invasion" else "p_invasion"
  has_prov <- "province" %in% names(risk_df)
  rd <- risk_df %>% dplyr::filter(horizon == !!horizon) %>%
    dplyr::transmute(.key = tolower(trimws(health_zone)),
                     # Normalise the province key like the zone key (see plot_invasion_risk_map).
                     .prov = if (has_prov) tolower(trimws(as.character(province))) else NA_character_,
                     p = .data[[pcol]], affected = was_active_before %in% TRUE,
                     .rank_pub = if ("rr_nat_rank" %in% names(risk_df)) as.numeric(rr_nat_rank)
                                 else NA_real_)
  # RANK invasion probability across all at-risk zones (1 = highest); affected zones unranked.
  # Fill CONTINUOUSLY over the FULL rank range (no cap, no binning): rank is UNIFORMLY spread
  # 1..N, unlike the highly-skewed probability (where 400+ zones pile at ~0 and look identical),
  # so the colour varies smoothly across EVERY zone and all 519 are distinguishable — including
  # the long tail beyond the top 50. Also carry a percentile for an interpretable legend.
  # THIS RANK IS THE PUBLISHED ONE: rr_nat_rank from compute_risk_scores() (15_workhorse.R),
  # min_rank(desc(mu)) over the at-risk zones, printed in bayes_risk_scores_all_zones.csv. It is
  # READ, not re-derived from p: the columns are posterior MEANS of mu and of p, and since
  # E[1 - exp(-mu)] != 1 - exp(-E[mu]) their orderings can differ (adjacent zones swapped on
  # the 2026-10-09 run). Only a frame without rr_nat_rank falls back to rank(-p, "min").
  #
  # NOT the same quantity as rank_med/rank_lo/rank_hi in that file, which are the POSTERIOR
  # MEDIAN rank and its 90% CrI (21_bayesian_renewal.R). Those describe how the ranking varies
  # across draws; these labels are the point-estimate ranking. Both are published, they answer
  # different questions, and the map deliberately shows the point estimate.
  n_ar <- sum(!(rd$affected) & is.finite(rd$p))
  rd <- rd %>% dplyr::mutate(
    rank_p = ifelse(affected | !is.finite(p), NA_real_,
                    if (all(is.na(.rank_pub)))
                      rank(ifelse(affected | !is.finite(p), NA_real_, -p), ties.method = "min",
                           na.last = "keep")       # among the at-risk zones only
                    else .rank_pub),
    pct = 100 * (1 - (rank_p - 1) / max(n_ar - 1, 1)))   # 100 = top, 0 = lowest at-risk
  shp <- shapefile %>% dplyr::mutate(.key = tolower(trimws(Nom)),
                                     .prov = tolower(trimws(as.character(PROVINCE))))
  m <- if (has_prov && !all(is.na(rd$.prov)))
         dplyr::left_join(shp, dplyr::select(rd, .key, .prov, p, affected, rank_p, pct), by = c(".key", ".prov"))
       else dplyr::left_join(shp, dplyr::select(rd, .key, p, affected, rank_p, pct), by = ".key")
  # s2-safe centroids for the top-K rank labels (some national polygons have invalid rings)
  .s2 <- suppressMessages(sf::sf_use_s2()); on.exit(suppressMessages(sf::sf_use_s2(.s2)), add = TRUE)
  suppressMessages(try(sf::sf_use_s2(FALSE), silent = TRUE))
  lab_df <- function(dat) {
    tl <- dat %>% dplyr::filter(!is.na(rank_p)) %>% dplyr::arrange(rank_p) %>% head(top_k)
    if (!nrow(tl)) return(NULL)
    ctr <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(
      sf::st_make_valid(sf::st_geometry(tl)))))
    data.frame(x = ctr[, 1], y = ctr[, 2], lab = as.integer(tl$rank_p))
  }
  base_layer <- function(dat, title, legend = FALSE, labels = FALSE) {
    g <- ggplot(dat) +
      geom_sf(aes(fill = pct), colour = "grey80", linewidth = 0.12) +
      geom_sf(data = dplyr::filter(dat, affected %in% TRUE),
              fill = "grey55", colour = "grey80", linewidth = 0.12) +
      scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 100),
                           name = "Invasion-risk\npercentile\n(100 = highest)",
                           breaks = c(0, 25, 50, 75, 100), na.value = "grey92") +
      labs(title = title) + theme_void(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 12),
            legend.position = if (legend) "right" else "none")
    if (labels) { ld <- lab_df(dat)
      if (!is.null(ld)) g <- g + geom_text(data = ld, aes(x, y, label = lab),
                                           size = 2.6, fontface = "bold", colour = "grey15") }
    g
  }
  ituri_zones <- shapefile$Nom[shapefile$PROVINCE == "Ituri"]
  m_ituri <- m %>% dplyr::filter(Nom %in% ituri_zones)
  sub <- sprintf("At-risk zones shaded by national invasion-risk percentile (rank of expected introductions, as published; bright = top); numbers label the top %d; dark grey = already affected.", top_k)
  if (extent == "national") {
    out <- base_layer(m, if (show_title)
                        sprintf("DRC national — invasion-risk RANK (%s, h=%dw)", method_label, horizon)
                      else NULL,
                      legend = TRUE, labels = TRUE) +
      labs(subtitle = if (show_title) sub else NULL, caption = window_txt) +
      theme(plot.subtitle = element_text(size = 9, colour = "grey35"),
            plot.caption = element_text(colour = "grey40", size = 8, hjust = 0))
    if (save) .save(out, file.path(OUT_MAPS, sprintf("%s_national_h%d.pdf", file, horizon)), w = 9, h = 7)
    return(invisible(out))
  }
  if (extent == "ituri") {
    out <- base_layer(m_ituri, sprintf("Ituri — invasion-risk RANK (%s, h=%dw)", method_label, horizon),
                      legend = TRUE, labels = TRUE) +
      labs(subtitle = sub, caption = window_txt) +
      theme(plot.subtitle = element_text(size = 9, colour = "grey35"),
            plot.caption = element_text(colour = "grey40", size = 8, hjust = 0))
    if (save) .save(out, file.path(OUT_MAPS, sprintf("%s_ituri_h%d.pdf", file, horizon)), w = 8, h = 7)
    return(invisible(out))
  }
  p_it  <- base_layer(m_ituri, sprintf("Ituri — invasion-risk RANK (%s, h=%dw)", method_label, horizon),
                      legend = TRUE, labels = TRUE) +
    labs(subtitle = sub) + theme(plot.subtitle = element_text(size = 9, colour = "grey35"))
  p_nat <- base_layer(m, "National context", labels = FALSE) +
    theme(plot.title = element_text(size = 10))
  out <- p_it + p_nat + patchwork::plot_layout(widths = c(2, 1))
  if (nzchar(window_txt))
    out <- out + patchwork::plot_annotation(caption = window_txt,
      theme = ggplot2::theme(plot.caption = element_text(colour = "grey40", size = 8, hjust = 0)))
  if (save) .save(out, file.path(OUT_MAPS, sprintf("%s_h%d.pdf", file, horizon)), w = 12, h = 6.5)
  invisible(out)
}

# ---------------------------------------------------------------------------
# 3. Risk-score ranked bars (top at-risk zones; Ituri + national)
# ---------------------------------------------------------------------------

plot_risk_scores_bars <- function(risk_df, horizon = 1L, top_n = 15L, save = TRUE,
                                   file = "risk_scores_bars", model_label = "") {
  d <- risk_df %>% dplyr::filter(horizon == !!horizon, !(was_active_before %in% TRUE),
                                 !is.na(p_case_invasion))
  .ml <- if (nzchar(model_label)) paste0(" — ", model_label) else ""
  .xlab <- sprintf("Probability of first case(s) within %d week%s (%s)", as.integer(horizon),
                   if (as.integer(horizon) > 1L) "s" else "", .risk_prob_scale(risk_df))
  mk <- function(dd, ttl) {
    dd <- dd %>% dplyr::arrange(dplyr::desc(p_case_invasion)) %>% head(top_n) %>%
      dplyr::mutate(health_zone = factor(health_zone, levels = rev(health_zone)))
    # P(case) ONLY. The companion "P(infection)" series was p_case / rho with rho a fixed,
    # guessed constant (0.45); it added no information the primary score did not already
    # carry, rescaled every bar by the same factor, and invited reading a guessed constant as
    # an estimate. Ascertainment has been removed from the pipeline entirely.
    long <- dd %>% dplyr::select(health_zone, `P(case)` = p_case_invasion) %>%
      tidyr::pivot_longer(-health_zone)
    ggplot(long, aes(value, health_zone, fill = name)) +
      geom_col(width = 0.65) +
      scale_fill_manual(values = c(`P(case)` = OKABE_ITO[1]), name = NULL) +
      scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
      labs(title = ttl, x = .xlab, y = NULL) +
      theme_inv(11)
  }
  p_it  <- mk(d %>% dplyr::filter(province == "Ituri"), sprintf("Ituri — highest invasion risk (h=%dw)%s", horizon, .ml))
  p_nat <- mk(d, sprintf("National — highest invasion risk%s", .ml))
  out <- p_it + p_nat
  if (save) .save(out, file.path(OUT_REPORTS, sprintf("%s_h%d.pdf", file, horizon)), w = 12, h = 6)
  invisible(out)
}

# ---------------------------------------------------------------------------
# 4. Discrimination summary — AUC-PR skill & rank-of-truth per model
# ---------------------------------------------------------------------------

#' Classify a model name into its methodological family (shared across the performance
#' plots so Bayesian renewal models are not mislabelled as naive baselines).
.invasion_family <- function(method) {
  m <- as.character(method)
  dplyr::case_when(
    # The Bayesian ensembles are labelled "Bayes-ens-mean" etc. (append_ensembles(),
    # prefix "Bayes-ens"), so they must be caught BEFORE the "^Bayes" branch.
    grepl("^Ensemble|-ens-", m)           ~ "Ensemble",
    grepl("^Bayes",    m)                 ~ "Bayesian renewal",
    grepl("^Renewal",  m)                 ~ "Renewal (mobility)",
    grepl("hhh4",      m)                 ~ "hhh4",
    grepl("SEIR",      m)                 ~ "SEIR",
    TRUE                                  ~ "Baseline (naive)")
}
.FAMILY_COLOURS <- c(`Renewal (mobility)` = OKABE_ITO[1], `Bayesian renewal` = OKABE_ITO[7],
                     hhh4 = OKABE_ITO[3], SEIR = OKABE_ITO[4],
                     `Baseline (naive)` = OKABE_ITO[2], Ensemble = OKABE_ITO[6])

#' @param restrict optional regex on `method` to show only a subset (e.g. "^Bayes" for a
#'   Bayesian-only view). @param label optional title prefix. @param file output stem.
plot_discrimination_summary <- function(eval_tbl, horizon = 1L, save = TRUE,
                                        restrict = NULL, label = NULL,
                                        file = "discrimination_summary") {
  d <- eval_tbl %>% dplyr::filter(horizon == !!horizon)
  if (!is.null(restrict)) d <- d %>% dplyr::filter(grepl(restrict, method))
  # In a family-restricted view, drop partial-CV models so the subset is compared on a
  # comparable common support (unless that would leave too few to plot).
  if (!is.null(restrict) && "partial_cv" %in% names(d) && sum(!(d$partial_cv %in% TRUE)) >= 2L)
    d <- d %>% dplyr::filter(!(partial_cv %in% TRUE))
  if (!nrow(d)) return(invisible(NULL))
  # PER-DATE AUC-PR SKILL AND MEAN RANK OF THE INVADED ZONES, OVER ALL ROUNDS: performance is
  # reported over every round (the selection applies the same estimators to the post-burn-in
  # rounds only, a deliberate choice; INVASION_SELECTION_BASIS, 16_invasion_eval.R). The pooled
  # skill rewards a probability LEVEL, not the ordering within a round, and is not drawn.
  .sel <- all(c("auc_pr_skill_date", "mean_rank_of_truth") %in% names(d))
  if (!.sel) stop("[viz] the evaluation table has no per-date AUC-PR column (auc_pr_skill_date); ",
                  "re-run the pipeline.", call. = FALSE)
  d <- d %>% dplyr::mutate(family = .invasion_family(method), sk = auc_pr_skill_date,
                           mr = mean_rank_of_truth,
                           method = forcats::fct_reorder(method, sk))
  .chance <- if ("auc_pr_skill_date_chance" %in% names(d))
    unique(round(d$auc_pr_skill_date_chance[is.finite(d$auc_pr_skill_date_chance)], 10)) else NA_real_
  .chance <- if (length(.chance) == 1L) .chance else NA_real_
  .ttl <- if (!is.null(label)) sprintf("%s invasion discrimination (h=%dw)", label, horizon)
          else sprintf("Invasion discrimination (h=%dw)", horizon)
  # The skill's interval (auc_pr_skill_date_lo/hi: zone-bootstrap replicates of the per-date
  # skill over all rounds); a table without the columns gets no bar.
  .has_ci <- all(c("auc_pr_skill_date_lo", "auc_pr_skill_date_hi") %in% names(d))
  p1 <- ggplot(d, aes(sk, method, colour = family)) +
    geom_segment(aes(x = 0, xend = sk, yend = method), colour = "grey80", linewidth = 0.5) +
    { if (.has_ci)
        geom_errorbar(aes(xmin = auc_pr_skill_date_lo, xmax = auc_pr_skill_date_hi),
                      orientation = "y", width = 0, linewidth = 0.5, na.rm = TRUE) } +
    geom_point(size = 3) +
    { if (is.finite(.chance)) geom_vline(xintercept = .chance, linetype = "dashed", colour = "grey50") } +
    scale_colour_manual(values = .FAMILY_COLOURS, name = NULL) +
    labs(title = .ttl,
         # Wrapped: each subtitle sits over part of the canvas, and unwrapped they ran into
         # each other and off the right edge.
         subtitle = paste(strwrap(paste0("Per-date AUC-PR skill over all rounds: AP within each ",
                           "round over its base rate, event-weighted. ",
                           if (is.finite(.chance)) sprintf("Dashed: a random ordering (%.1f). ", .chance) else "",
                           "Higher = better.",
                           if (.has_ci) " Bars: 90% zone-cluster bootstrap interval." else ""), 60),
                          collapse = "\n"),
         x = "Per-date AUC-PR skill (× base rate)", y = NULL) +
    theme_inv(11)
  .mrt_max <- max(d$mr, na.rm = TRUE)
  p2 <- ggplot(d, aes(mr, method, colour = family)) +
    geom_segment(aes(x = .mrt_max, xend = mr, yend = method),
                 colour = "grey80", linewidth = 0.5) +
    geom_point(size = 3) +
    scale_colour_manual(values = .FAMILY_COLOURS, name = NULL, guide = "none") +
    labs(title = "Mean rank of invaded zones",
         subtitle = paste(strwrap(paste0("All rounds: rank of EVERY invaded zone (tie-averaged), ",
                                         "meaned per round then over rounds; lower = nearer the top"), 50),
                          collapse = "\n"),
         x = "Mean rank among at-risk zones", y = NULL) +
    theme_inv(11)
  out <- p1 + p2 + patchwork::plot_layout(guides = "collect") & theme(legend.position = "top")
  # The canvas height MUST scale with the number of models on the y-axis. theme_inv(11) draws
  # axis text at ~8.8pt (~0.12 in per glyph), so a fixed 5.5-in canvas silently overprints the
  # labels once the grid exceeds ~25 models — and the Bayesian grid is currently 54, which
  # rendered the y-axis of bayes_discrimination_summary as an unreadable block of ink.
  # 0.16 in per row leaves ~0.04 in of leading; +1.8 in covers title/subtitle/legend/axis.
  n_lab <- dplyr::n_distinct(d$method)
  h_fig <- max(5.5, 1.8 + 0.16 * n_lab)
  if (save) .save(out, file.path(OUT_DIAGNOSTICS, sprintf("%s_h%d.pdf", file, horizon)),
                  w = 12, h = h_fig)
  invisible(out)
}

#' Models vs baselines/simple models (#8): does mechanistic mobility structure actually
#' buy skill over naive distance/gravity baselines and simple comparators? Two panels over
#' the SAME leave-future-out folds — discrimination (AUC-PR skill) and operational detection
#' (recall@K at the TOPK_VALUES budget invasion_topk_cols() resolves, i.e. the share of the
#' true invasions caught by a K-zone watch-list) — with every model
#' drawn, coloured by family, the naive-baseline band shaded, and the no-skill / random-watch
#' -list references marked. Makes the renewal/Bayesian-over-baseline gap explicit.
plot_model_vs_baseline <- function(eval_tbl, horizon = 1L, save = TRUE) {
  d <- eval_tbl %>% dplyr::filter(horizon == !!horizon)
  if (!nrow(d)) return(invisible(NULL))
  # partial-CV models are not comparably scored — drop them from the head-to-head.
  if ("partial_cv" %in% names(d)) d <- d %>% dplyr::filter(!(partial_cv %in% TRUE))
  d <- d %>% dplyr::mutate(family = .invasion_family(method),
                           is_baseline = family == "Baseline (naive)")
  .tk  <- invasion_topk_cols(d)
  if (is.null(.tk)) stop("[viz] the evaluation table carries no recall_at_<K> column.", call. = FALSE)
  rcol <- .tk$recall; kk <- .tk$k
  # recall@K is a PER-FOLD-averaged metric, so the random watch-list baseline is K / (at-risk
  # zones PER FOLD) = K * n_folds / n_atrisk (n_atrisk is pooled across folds); dividing by the
  # pooled n_atrisk would understate the random reference ~n_folds-fold.
  # n_folds_scored, NOT n_folds. Since the common-support fix n_folds is each method's NATIVE
  # reach while n_folds_scored is what its metrics were actually computed on; taking the median
  # of the former put the random reference on a fold count no row was scored at whenever the
  # grid's composition shifted the median.
  .nfc <- if ("n_folds_scored" %in% names(d)) "n_folds_scored" else "n_folds"
  .nf <- if (.nfc %in% names(d)) stats::median(d[[.nfc]], na.rm = TRUE) else 1
  rand <- if ("n_atrisk" %in% names(d)) kk * .nf / stats::median(d$n_atrisk, na.rm = TRUE) else NA_real_
  base_best <- suppressWarnings(max(d$auc_pr_skill[d$is_baseline], na.rm = TRUE))
  dA <- d %>% dplyr::mutate(method = forcats::fct_reorder(method, auc_pr_skill))
  pA <- ggplot(dA, aes(auc_pr_skill, method, colour = family)) +
    { if (is.finite(base_best))
        annotate("rect", xmin = 0, xmax = base_best, ymin = -Inf, ymax = Inf,
                 fill = "grey85", alpha = 0.5) } +
    geom_segment(aes(x = 1, xend = auc_pr_skill, yend = method), colour = "grey80", linewidth = 0.4) +
    geom_point(size = 2.6) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    { if (is.finite(base_best)) geom_vline(xintercept = base_best, linetype = "dotted", colour = OKABE_ITO[2]) } +
    scale_colour_manual(values = .FAMILY_COLOURS, name = NULL) +
    labs(title = sprintf("Discrimination vs baselines (h=%dw)", horizon),
         subtitle = "AUC-PR skill = AUC-PR / base rate; 1 = a constant forecast (dashed). Shaded/dotted = best naive baseline.",
         x = "AUC-PR skill (x base rate)", y = NULL) +
    theme_inv(10)
  dB <- d %>% dplyr::filter(is.finite(.data[[rcol]])) %>%
    dplyr::mutate(method = forcats::fct_reorder(method, .data[[rcol]]))
  pB <- ggplot(dB, aes(.data[[rcol]], method, colour = family)) +
    geom_segment(aes(x = 0, xend = .data[[rcol]], yend = method), colour = "grey80", linewidth = 0.4) +
    geom_point(size = 2.6) +
    { if (is.finite(rand)) geom_vline(xintercept = rand, linetype = "dashed", colour = "grey50") } +
    scale_colour_manual(values = .FAMILY_COLOURS, name = NULL, guide = "none") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = sprintf("Detection: top-%d watch-list recall", kk),
         subtitle = sprintf("Share of true invasions caught by a %d-zone list; dashed = random watch-list (%.0f%%).",
                            kk, 100 * (rand %||% NA_real_)),
         x = sprintf("recall@%d", kk), y = NULL) +
    theme_inv(10)
  out <- pA + pB + patchwork::plot_layout(guides = "collect") & theme(legend.position = "top")
  if (save) .save(out, file.path(OUT_DIAGNOSTICS, sprintf("model_vs_baseline_h%d.pdf", horizon)), w = 13, h = 6)
  invisible(out)
}

# ---------------------------------------------------------------------------
# 5. LFO forecast-vs-outcome tile (at-risk zones that were invaded)
# ---------------------------------------------------------------------------

plot_lfo_forecast_vs_outcome <- function(lfo_results, method, horizon = 1L, save = TRUE,
                                         file = "lfo_forecast_vs_outcome") {
  d <- lfo_results %>% dplyr::filter(method == !!method, horizon == !!horizon)
  invaded_zones <- d %>% dplyr::filter(is_new_invasion == 1) %>% dplyr::pull(health_zone) %>% unique()
  if (length(invaded_zones) == 0) return(invisible(NULL))
  # For each invaded zone, show predicted p at each fold cutoff + mark the invasion fold
  # THE AS-ISSUED SCALE: the prequential recalibrated probability (delta fitted only on outcomes
  # settled by each round's issue date) when every row carries one; else the raw probability.
  # The legend says which.
  .use_recal <- "p_recal" %in% names(d) && all(is.finite(d$p_recal))
  # ROUNDS WITHOUT A FITTED FACTOR ARE RAW. Before the first post-burn-in round is settled by a
  # round's issue date there is nothing to fit delta on, so that round's p_recal IS p_invasion
  # (delta_preq = 1, delta_estimable = FALSE; 16b). One colour scale would otherwise mix raw
  # early columns with recalibrated later ones unannounced, so those columns are outlined
  # (dashed) and the legend and subtitle say so. Derived from delta_estimable, never assumed.
  .unfit_origins <- if (.use_recal && "delta_estimable" %in% names(d))
    sort(unique(lfo_origin(d$cutoff[!(d$delta_estimable %in% TRUE)]))) else as.Date(character(0))
  .first_fit <- if (.use_recal && "delta_estimable" %in% names(d) && any(d$delta_estimable %in% TRUE))
    min(lfo_origin(d$cutoff[d$delta_estimable %in% TRUE])) else as.Date(NA)
  dd <- d %>% dplyr::filter(health_zone %in% invaded_zones) %>%
    dplyr::mutate(cutoff = lfo_origin(cutoff),   # label rounds by forecast origin
                  health_zone = factor(health_zone, levels = rev(sort(invaded_zones))),
                  .p_fill = if (.use_recal) p_recal else p_invasion)
  .lgd <- if (!.use_recal) "Predicted P\n(raw)"
          else if (!length(.unfit_origins)) "Predicted P\n(recalibrated,\nas issued)"
          else if (is.na(.first_fit)) "Predicted P\n(raw: no round\ncarries a fitted\nfactor)"
          else "Predicted P\n(recalibrated,\nas issued;\ndashed = raw)"
  # The outcome window is CUMULATIVE (weeks C+7 .. C+7h), so at h >= 2 consecutive rounds'
  # windows overlap and one first case is marked in every round whose window contains it.
  .sub_scale <- if (!.use_recal || !length(.unfit_origins)) ""
                else if (is.na(.first_fit))
                  ". No round carries a fitted recalibration factor: all probabilities are raw."
                else sprintf(paste0(". Dashed columns (origins before %s) carry no fitted ",
                                    "recalibration factor yet and are drawn uncorrected (raw)."),
                             format(.first_fit, "%d %b"))
  .sub <- paste0("Rows = zones that were invaded; x = FORECAST ORIGIN; X = a round whose ",
                 if (as.integer(horizon) > 1L)
                   sprintf(paste0("%d-week outcome window contains the zone's first case (windows ",
                                  "overlap, so one invasion is marked in up to %d consecutive rounds)"),
                           as.integer(horizon), as.integer(horizon))
                 else "outcome week contains the zone's first case",
                 .sub_scale)
  p <- ggplot(dd, aes(cutoff, health_zone, fill = .p_fill)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    { if (.use_recal && length(.unfit_origins) && !is.na(.first_fit))
        geom_tile(data = dplyr::filter(dd, cutoff %in% .unfit_origins), fill = NA,
                  colour = "grey20", linewidth = 0.35, linetype = "22") } +
    geom_point(data = dplyr::filter(dd, is_new_invasion == 1),
               shape = 4, size = 3, stroke = 1.1, colour = "black") +
    scale_fill_viridis_c(option = "plasma", name = .lgd,
                         labels = percent_format(accuracy = 1)) +
    # pretty_breaks rather than the default: with the cross-validation now running to the last
    # round whose outcome window closes by the analysis date, the tile has ~16 columns and the
    # default breaks print one label per column.
    scale_x_date(date_labels = "%d %b", breaks = scales::pretty_breaks(8)) +
    labs(title = sprintf("LFO: predicted invasion risk vs realised invasions (%s, h=%dw)", method, horizon),
         subtitle = paste(strwrap(.sub, width = 105), collapse = "\n"),
         x = "Forecast origin (as-of date)", y = NULL) +
    theme_inv(11)
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, sprintf("%s_h%d.pdf", file, horizon)),
                  w = 9, h = max(4, 0.35 * length(invaded_zones) + 2))
  invisible(p)
}

# ---------------------------------------------------------------------------
# 5b. Spatiotemporal skill over time (Q5) — per-fold discrimination vs date
# ---------------------------------------------------------------------------

#' The featured model a figure is drawn for: the caller's `featured`, else this session's
#' selection (run_all.R's `best_bayes_method`, set by the two-stage selection before any figure
#' is drawn; key_outputs/model_selection.json is rewritten only AFTER the figures, so reading it
#' there would return the previous run's pick), else the run's selection record through
#' read_featured_model() (figure_style.R) when that is loaded. Never a typed name or a guess:
#' with none of the three this stops.
.featured_for_figure <- function(featured = NULL) {
  .ok <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  if (.ok(featured)) return(featured)
  f <- get0("best_bayes_method", envir = globalenv(), ifnotfound = NULL)
  if (.ok(f)) return(f)
  if (exists("read_featured_model", mode = "function"))
    return(read_featured_model(file.path(OUT_DIR, "key_outputs", "model_selection.json")))
  stop("[viz] the featured model is unknown: no `featured` argument, no best_bayes_method in ",
       "this session, and read_featured_model() (figure_style.R) is not loaded. The figure is ",
       "not drawn rather than guessing which model it is about.", call. = FALSE)
}

#' @param st_skill output of spatiotemporal_skill(); @param methods optional
#'   subset of method names to draw (keeps the panel legible); the featured model is added to it.
#' @param featured the featured model, ALWAYS drawn (labelled "(featured)"); NULL resolves it
#'   through .featured_for_figure().
plot_skill_over_time <- function(st_skill, metric = "auc_pr_skill",
                                 methods = NULL, featured = NULL, save = TRUE) {
  if (is.null(st_skill) || nrow(st_skill) == 0) return(invisible(NULL))
  featured <- .featured_for_figure(featured)
  d <- st_skill
  if (!is.null(methods)) d <- d %>% dplyr::filter(method %in% union(methods, featured))
  # X AXIS = THE WEEK BEING FORECAST, not the round's cutoff.
  #
  # A round with cutoff C at horizon h predicts the window (C, C + 7h], which ENDS on
  # C + 7h + 6. Plotting against C put the last point five weeks before the end of the data
  # and made the two horizon facets disagree about which calendar week a point referred to —
  # the same cutoff means a different target week at h=1 and h=2. Plotting against the target
  # week's last day puts both facets on one calendar footing and makes the series run to the
  # analysis date, which is what "skill over time" is asking about.
  d <- d %>% dplyr::mutate(cutoff = as.Date(cutoff),
                           target_end = cutoff + 7L * as.integer(horizon) + 6L)
  # ROUNDS WITH NOTHING TO SCORE ARE MARKED, NOT DELETED. Discrimination is undefined without
  # a positive, so an event-free round has no point — but deleting it silently leaves a gap
  # the reader will read as missing data or as a model failure. They are drawn as ticks on the
  # axis and named in the subtitle. (A round can be event-free either because no zone was
  # invaded in its window or, for the newest rounds, because the invasions that did occur are
  # not yet laboratory-confirmed.)
  empty <- if ("n_events" %in% names(d))
    d %>% dplyr::filter(n_events %in% 0) %>% dplyr::distinct(horizon, target_end) else
    d[0, c("horizon", "target_end"), drop = FALSE]
  d <- d %>% dplyr::filter(is.finite(.data[[metric]]))
  if (nrow(d) == 0) return(invisible(NULL))
  if (!featured %in% d$method)
    stop(sprintf("[viz] skill over time: the featured model %s has no finite %s in the per-fold table.",
                 featured, metric), call. = FALSE)
  # <= 6 methods for legibility (Okabe-Ito): THE FEATURED MODEL, plus the five best others by mean
  # per-fold skill. The pool for those five holds only models that competed as themselves: the
  # sensitivity arms (SENSITIVITY_ARM_REGEX: refits of one specification at a changed assumption)
  # and the time-varying forms other than the featured one (tv_form_of(): stage-2 refits of the
  # stage-1 pick, compared in their own table) are left out, as from selection.
  .mm <- d %>% dplyr::group_by(method) %>%
    dplyr::summarise(m = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(m))
  .pool <- .mm$method[.mm$method != featured & !grepl(SENSITIVITY_ARM_REGEX, .mm$method) &
                      is.na(tv_form_of(.mm$method))]
  keep <- .mm$method[.mm$method %in% c(featured, head(.pool, 5L))]   # mean-skill order
  .hlab <- function(v) paste0("h = ", v, " week", ifelse(v > 1, "s", ""))
  d <- d %>% dplyr::filter(method %in% keep) %>%
    dplyr::mutate(method = factor(method, levels = keep), horizon = .hlab(horizon))
  if (nrow(empty)) empty <- empty %>% dplyr::mutate(horizon = .hlab(horizon))
  # The per-fold auc_pr_skill from 19_spacetime_eval.R is now the prevalence-NORMALISED skill
  # (ap - base)/(1 - base), bounded in [0,1], not the raw lift ap/base. Label and reference line
  # follow: 0 = no better than prevalence, not 1.
  ylab <- switch(metric, auc_pr_skill = "AUC-PR skill  (0 = prevalence, 1 = perfect)",
                 hit_at_k = "Hit-rate @K (top-K contains an invaded zone)",
                 prec_at_k = "Precision @K", metric)
  p <- ggplot(d, aes(target_end, .data[[metric]], colour = method, group = method)) +
    { if (metric == "auc_pr_skill") geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") } +
    { if (nrow(empty)) geom_rug(data = empty, mapping = aes(x = target_end), sides = "b",
                                colour = "grey55", linewidth = 0.6, length = unit(4, "pt"),
                                inherit.aes = FALSE) } +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_wrap(~ horizon) +
    scale_colour_manual(values = setNames(OKABE_ITO[seq_along(keep)], keep),
                        labels = setNames(ifelse(keep == featured, paste0(keep, " (featured)"), keep), keep),
                        name = NULL) +
    scale_x_date(date_labels = "%d %b", breaks = scales::pretty_breaks(6)) +
    labs(title = "Forecast skill over time",
         subtitle = paste(strwrap(paste0("Per-fold skill against the LAST DAY OF THE WEEK BEING FORECAST, so both horizons share a calendar. ",
                           if (nrow(empty))
                             sprintf("Ticks on the x-axis mark %d round(s) in which no invasion was recorded, where discrimination is undefined. ",
                                     nrow(empty)) else "",
                           if (identical(metric, "auc_pr_skill"))
                             "Prevalence-normalised, so rounds with different event counts are comparable." else ""),
                          140), collapse = "\n"),
         x = "Last day of the forecast window", y = ylab) +
    theme_inv(11)
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, sprintf("skill_over_time_%s.pdf", metric)),
                  w = 11, h = 5)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 5b-bis. AUC-PR decomposition (Q5) — raw AUC-PR vs its prevalence floor, and the lift
# ---------------------------------------------------------------------------

#' WHY THIS EXISTS ALONGSIDE plot_skill_over_time(). That figure plots ONE number per round,
#' the prevalence-normalised skill (ap - base)/(1 - base). This one takes the same rounds
#' apart into the three quantities the single number compresses: the raw AUC-PR, the
#' prevalence floor it is being judged against, and the ratio of the two. It is the evidence
#' behind the normalisation choice documented in 19_spacetime_eval.R — the bottom row's lift
#' peaks on the rounds with the FEWEST events (its attainable maximum is n / n_pos), so read
#' against the top row it shows the lift tracking the event count rather than the ranking.
#'
#' Four panels: top row raw AUC-PR against the floor, bottom row the lift, each by horizon.
#'
#' @param st_skill output of spatiotemporal_skill().
#' @param method the one method to draw (the featured model). REQUIRED: the figure carries no
#'   title, so nothing inside it would reveal a substituted model, and it does not fall back to
#'   the best-scoring one.
plot_aucpr_decomposition <- function(st_skill, method, save = TRUE,
                                     file = "aucpr_decomposition_by_horizon") {
  if (missing(method) || is.null(method) || length(method) != 1L || is.na(method))
    stop("[viz] plot_aucpr_decomposition() needs `method` (the featured model); it has no ",
         "title, so a substituted model would be invisible.", call. = FALSE)
  if (is.null(st_skill) || nrow(st_skill) == 0) return(invisible(NULL))
  d <- st_skill
  if (!method %in% d$method) {
    warning(sprintf("[viz] AUC-PR decomposition: %s is not in the per-fold table; not drawn.", method),
            call. = FALSE)
    return(invisible(NULL))
  }
  d <- d %>% dplyr::filter(method == !!method)
  # SAME CALENDAR AS ITS SIBLING. plot_skill_over_time() plots against the last day of the
  # week being forecast so h=1 and h=2 refer to the same calendar week; this figure is read
  # beside it, so it uses the same x.
  d <- d %>% dplyr::mutate(cutoff = as.Date(cutoff),
                           target_end = cutoff + 7L * as.integer(horizon) + 6L)
  .hlab <- function(v) paste0(v, "-week ahead")
  empty <- if ("n_events" %in% names(d))
    d %>% dplyr::filter(n_events %in% 0) %>%
      dplyr::distinct(horizon, target_end) %>% dplyr::mutate(horizon = .hlab(horizon)) else
    d[0, c("horizon", "target_end"), drop = FALSE]
  d <- d %>% dplyr::filter(is.finite(auc_pr)) %>% dplyr::mutate(horizon = .hlab(horizon))
  if (nrow(d) == 0) return(invisible(NULL))

  .rug <- function() if (nrow(empty))
    geom_rug(data = empty, mapping = aes(x = target_end), sides = "b",
             colour = "grey55", linewidth = 0.6, length = unit(4, "pt"), inherit.aes = FALSE)
  col_ap <- OKABE_ITO[1]; col_base <- "#B2182B"; col_lift <- OKABE_ITO[3]

  # TOP: the raw AUC-PR and the floor it is measured against, on one axis. The digit above
  # each point is that round's positive count -- the denominator of both the floor and the
  # lift's ceiling, and the thing the bottom row is actually tracking.
  top <- d %>%
    dplyr::select(horizon, target_end, n_events, auc_pr, base_rate) %>%
    tidyr::pivot_longer(c(auc_pr, base_rate), names_to = "series", values_to = "value") %>%
    dplyr::mutate(series = factor(series, levels = c("auc_pr", "base_rate"),
                                  labels = c("AUC-PR (raw)", "Prevalence floor (base rate)")))
  p_top <- ggplot(top, aes(target_end, value, colour = series, group = series)) +
    .rug() +
    geom_line(aes(linetype = series), linewidth = 0.8) + geom_point(size = 1.9) +
    geom_text(data = dplyr::filter(top, series == "AUC-PR (raw)"),
              aes(label = n_events), vjust = -0.9, size = 2.6, show.legend = FALSE) +
    facet_wrap(~ horizon) +
    scale_colour_manual(values = c("AUC-PR (raw)" = col_ap,
                                   "Prevalence floor (base rate)" = col_base), name = NULL) +
    scale_linetype_manual(values = c("AUC-PR (raw)" = "solid",
                                     "Prevalence floor (base rate)" = "21"), name = NULL) +
    scale_x_date(date_labels = "%d %b", breaks = scales::pretty_breaks(6)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = NULL, y = "AUC-PR vs its prevalence floor") +
    theme_inv(11)

  # BOTTOM: the ratio of the two rows above. Its scale is unbounded and round-specific, which
  # is the point of showing it under them rather than on its own.
  p_bot <- ggplot(d, aes(target_end, auc_pr_lift)) +
    .rug() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
    geom_line(colour = col_lift, linewidth = 0.8) +
    geom_point(colour = col_lift, size = 1.9) +
    facet_wrap(~ horizon) +
    scale_x_date(date_labels = "%d %b", breaks = scales::pretty_breaks(6)) +
    labs(x = "Last day of the forecast window",
         y = "Lift over floor (AUC-PR / base rate)") +
    theme_inv(11) + theme(strip.text = element_blank())

  # NO TITLE OR SUBTITLE, by request. Nothing inside the figure now records which model it
  # shows, that the digits are per-round invasion counts, or that the lift's ceiling is
  # n_at-risk / n_events and so rises when events are scarce -- the caption has to supply all
  # three. The axis titles and the legend are the only text left.
  p <- (p_top / p_bot) + patchwork::plot_layout(heights = c(1, 0.85))
  # Height trimmed from 7.5: the three subtitle lines it allowed for are gone.
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, paste0(file, ".pdf")), w = 11, h = 6.8)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 5c. Space-time risk heatmap (Q5) — top-risk zones × cutoff, invasions marked
# ---------------------------------------------------------------------------

#' Zones invaded during the back-test, then the highest peak-risk zones that never were, up to
#' `top_n` rows, across every fold, with the realised first-case folds marked — surfaces both
#' good calls and chronic false alarms in one legible tile plot.
plot_spacetime_risk <- function(lfo_results, method, horizon = 1L, top_n = 20L,
                                save = TRUE, file = "spacetime_risk") {
  d <- lfo_results %>% dplyr::filter(method == !!method, horizon == !!horizon,
                                     is.finite(p_invasion))
  if (nrow(d) == 0) return(invisible(NULL))
  .raw <- .lfo_is_raw(d)
  top <- d %>% dplyr::group_by(health_zone) %>%
    dplyr::summarise(peak = max(p_invasion, na.rm = TRUE),
                     inv = sum(is_new_invasion, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(inv > 0), dplyr::desc(peak)) %>%
    dplyr::slice_head(n = top_n) %>% dplyr::pull(health_zone)
  dd <- d %>% dplyr::filter(health_zone %in% top) %>%
    dplyr::mutate(cutoff = lfo_origin(cutoff),   # label rounds by forecast origin
                  health_zone = factor(health_zone, levels = rev(top)))
  p <- ggplot(dd, aes(cutoff, health_zone, fill = p_invasion)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_point(data = dplyr::filter(dd, is_new_invasion == 1),
               shape = 4, size = 3, stroke = 1.1, colour = "black") +
    scale_fill_viridis_c(option = "plasma",
                         name = if (.raw) "Predicted P\n(raw)"
                                else "Predicted P\n(recalibrated,\nas issued; raw\nbefore a factor\nwas fitted)",
                         labels = percent_format(accuracy = 1)) +
    scale_x_date(date_labels = "%d %b") +
    labs(title = sprintf("Space-time invasion risk (%s, h=%dw)", method, horizon),
         subtitle = paste0("Zones invaded in the back-test first, then the highest peak-risk zones; x = FORECAST ORIGIN; ",
                           "X = a round whose outcome window contains the zone's first case",
                           if (as.integer(horizon) > 1L)
                             sprintf(" (windows overlap, so up to %d consecutive rounds)", as.integer(horizon))
                           else ""),
         x = "Forecast origin (as-of date)", y = NULL) +
    theme_inv(11)
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, sprintf("%s_h%d.pdf", file, horizon)),
                  w = 9, h = max(4, 0.32 * length(top) + 2))
  invisible(p)
}

# ---------------------------------------------------------------------------
# 5d. Lead-time / anticipation (Q5) — weeks ahead each zone was flagged
# ---------------------------------------------------------------------------

#' @param lead_df output of lead_time_analysis().
plot_lead_time <- function(lead_df, save = TRUE) {
  if (is.null(lead_df) || nrow(lead_df) == 0) return(invisible(NULL))
  d <- lead_df %>%
    dplyr::mutate(health_zone = forcats::fct_reorder(health_zone, anticipation_weeks))
  p <- ggplot(d, aes(anticipation_weeks, health_zone, colour = ever_topk)) +
    geom_segment(aes(x = 0, xend = anticipation_weeks, yend = health_zone),
                 colour = "grey80", linewidth = 0.5) +
    geom_point(size = 3) +
    scale_colour_manual(values = c(`TRUE` = OKABE_ITO[3], `FALSE` = OKABE_ITO[2]),
                        labels = c(`TRUE` = "Flagged in top-K", `FALSE` = "Never top-K"),
                        name = NULL) +
    labs(title = "Lead time: anticipation before first case",
         subtitle = "Weeks ahead each newly affected zone entered the top-K risk list (0 = never)",
         x = "Anticipation (weeks)", y = NULL) +
    theme_inv(11)
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, "lead_time.pdf"),
                  w = 8, h = max(3.5, 0.32 * nrow(d) + 1.5))
  invisible(p)
}

# ---------------------------------------------------------------------------
# 5e. Spatial-error map (Q5) — where the hits, misses and false alarms sit
# ---------------------------------------------------------------------------

#' @param zone_err output of zone_spatial_error(); @param shapefile province sf.
plot_spatial_error_map <- function(zone_err, shapefile = NULL, save = TRUE,
                                   file = "spatial_error_map", model_label = "",
                                   horizon = 1L) {
  if (is.null(zone_err) || nrow(zone_err) == 0) return(invisible(NULL))
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile)) return(invisible(NULL))
  ze <- zone_err %>% dplyr::transmute(.key = tolower(trimws(health_zone)), flag)
  shp <- shapefile %>% dplyr::mutate(.key = tolower(trimws(Nom)))
  m <- dplyr::left_join(shp, ze, by = ".key")
  ituri_zones <- shapefile$Nom[shapefile$PROVINCE == "Ituri"]
  m_it <- m %>% dplyr::filter(Nom %in% ituri_zones)
  pal <- c(hit = OKABE_ITO[3], miss = OKABE_ITO[2], false_alarm = OKABE_ITO[5],
           correct_negative = "grey85")
  labs_flag <- c(hit = "Hit (invaded, top-K)", miss = "Miss (invaded, not top-K)",
                 false_alarm = "False alarm (high risk, no case)",
                 correct_negative = "Correct negative")
  mk <- function(dat, ttl) ggplot(dat) +
    geom_sf(aes(fill = flag), colour = "grey80", linewidth = 0.12) +
    scale_fill_manual(values = pal, labels = labs_flag, na.value = "grey90", name = NULL) +
    labs(title = ttl) + theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12), legend.position = "right")
  p_it <- mk(m_it, sprintf("Ituri — spatial error map (h=%dw)%s", horizon,
                           if (nzchar(model_label)) paste0(" — ", model_label) else ""))
  if (save) .save(p_it, file.path(OUT_MAPS, sprintf("%s_h%d.pdf", file, horizon)), w = 9, h = 6.5)
  invisible(p_it)
}

# ---------------------------------------------------------------------------
# 6. Nowcast fan (kept)
# ---------------------------------------------------------------------------
plot_nowcast_fan2 <- function(factors_df, save = TRUE) {
  if (is.character(factors_df)) {
    if (!file.exists(factors_df)) return(invisible(NULL))
    factors_df <- readr::read_csv(factors_df, show_col_types = FALSE)
  }
  if (is.null(factors_df) || nrow(factors_df) == 0) return(invisible(NULL))
  df <- factors_df %>% dplyr::mutate(week_start = as.Date(week_start))
  p <- ggplot(df, aes(week_start)) +
    geom_ribbon(aes(ymin = nowcast_lo, ymax = nowcast_hi), fill = OKABE_ITO[1], alpha = 0.18) +
    geom_col(aes(y = observed), width = 4, fill = "grey70", alpha = 0.7) +
    geom_line(aes(y = nowcast_mean), colour = OKABE_ITO[1], linewidth = 0.9) +
    geom_point(aes(y = nowcast_mean), colour = OKABE_ITO[1], size = 2) +
    geom_text(aes(y = nowcast_hi, label = sprintf("×%.1f", factor)), vjust = -0.6, size = 3) +
    scale_x_date(date_labels = "%d %b") +
    # nowcast_lo / nowcast_hi are the 5% / 95% quantiles of each posterior draw's WEEKLY total
    # (epinowcast_weekly_factors(), 04b_epinowcast.R), so the band is a weekly 90% interval;
    # they are NA, and no band is drawn, when the draws could not be read.
    labs(title = "epinowcast right-truncation correction",
         subtitle = paste0("Grey = observed (right-censored); blue = nowcast mean, band = 90% interval ",
                           "of the weekly total; ×N = correction factor"),
         x = "Onset week", y = "Confirmed cases", caption = paste0("Analysis date: ", ANALYSIS_DATE)) +
    theme_inv(11)
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, "nowcast_fan.pdf"), w = 8, h = 5)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Invasion analysis report (markdown)
# ---------------------------------------------------------------------------

#' Write a concise invasion-analysis report (markdown) from the evaluation table
#' and the primary-model risk scores.
#' @param prob_scale the probability scale of `risk_scores` ("recalibrated-deployed" or "raw");
#'   NULL reads the frame's `prob_scale` column, else INVASION_RECALIBRATE_DEPLOY
#'   (.risk_prob_scale()).
write_invasion_report <- function(eval_tbl, risk_scores, lfo_results = NULL,
                                  zone_week = NULL, training_cutoff = NULL,
                                  best_method = NULL, primary_method = NULL,
                                  n_models = NA, out_path = NULL,
                                  best_bayes_method = NULL, bayes_weights = NULL,
                                  prob_scale = NULL) {
  if (is.null(out_path)) out_path <- file.path(OUT_REPORTS, "invasion_report.md")
  affected <- if (!is.null(zone_week)) {
    zone_week %>% dplyr::filter(week_start <= training_cutoff, confirmed > 0) %>%
      dplyr::pull(health_zone) %>% unique()
  } else character(0)
  n_zones <- if (!is.null(zone_week)) dplyr::n_distinct(zone_week$health_zone) else NA
  outbreak_summary <- list(
    # `cutoff` is the human-facing END of the training window (last calendar day
    # trained on); the affected-zone filter above uses the week-start anchor. Under
    # the analysis-date-anchored weekly grid these differ by 6 days (week end vs start).
    analysis_date = ANALYSIS_DATE,
    cutoff = if (!is.null(training_cutoff)) training_cutoff + 6L else training_cutoff,
    total_confirmed = if (!is.null(zone_week)) sum(zone_week$confirmed, na.rm = TRUE) else NA,
    zones_affected = length(affected),
    n_atrisk = n_zones - length(affected), n_zones = n_zones)
  # primary_method named the CV-selected best FREQUENTIST RENEWAL model. That arm has been
  # removed, so run_all.R now passes NULL and the fallback below resolves to best_method — the
  # featured BAYESIAN model. Every line that used to describe primary_method as a renewal model
  # would then have applied a false label to it, so those lines are conditioned on whether a
  # genuine second (non-Bayesian) featured model exists.
  .have_primary <- !is.null(primary_method) && length(primary_method) == 1L &&
                   !is.na(primary_method) && !identical(primary_method, best_method)
  if (is.null(primary_method)) primary_method <- best_method %||% "(featured model)"
  # The scale of the CURRENT-FORECAST probabilities this report prints (P(case) tables).
  .rs_scale <- .risk_prob_scale(risk_scores, prob_scale)
  .rs_recal <- identical(.rs_scale, "recalibrated, as deployed")
  fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
  L <- c()
  L <- c(L, "# BDBV 2026 — Spatiotemporal Invasion Forecast",
         "",
         sprintf("**Analysis date:** %s  ", outbreak_summary$analysis_date),
         # The LAST DAY trained on (cutoff week start + 6), not the pipeline's cutoff C, which is
         # the first day of that week.
         sprintf("**Training data through:** %s  ", outbreak_summary$cutoff),
         sprintf("**Confirmed cases:** %d across %d affected health zones  ",
                 outbreak_summary$total_confirmed, outbreak_summary$zones_affected),
         sprintf("**At-risk zones (no cases yet):** %d of %d  ",
                 outbreak_summary$n_atrisk, outbreak_summary$n_zones),
         {
           # Fit window starts at the outbreak (training) start, not the earliest
           # historical week in the grid, which carries pre-outbreak zeros.
           .obstart <- get0("OUTBREAK_START", ifnotfound = NA)
           .datamin <- if (!is.null(zone_week)) suppressWarnings(min(zone_week$week_start, na.rm = TRUE)) else NA
           .ts <- if (!is.na(.obstart) && !is.na(.datamin)) max(as.Date(.obstart), as.Date(.datamin)) else .datamin
           .hz <- get0("LFO_HORIZONS", ifnotfound = c(1L, 2L))
           if (!is.na(.ts) && !is.null(training_cutoff))
             sprintf(paste0("**Forecast windows:** fit %s to %s (training); ",
                            "prediction %s to %s (next %s).  "),
                     format(as.Date(.ts), "%d %b %Y"), format(as.Date(training_cutoff) + 6, "%d %b %Y"),
                     format(as.Date(training_cutoff) + 7, "%d %b %Y"),
                     # The last predicted week ends on cutoff + 7*hmax + 6 (inclusive), as in
                     # .window_caption() and run_info.json; "+ 7*hmax" stopped six days short.
                     format(as.Date(training_cutoff) + max(.hz) * 7 + 6, "%d %b %Y"),
                     if (max(.hz) > 1L) sprintf("%d-%d weeks", min(.hz), max(.hz)) else "week")
           else ""
         },
         "",
         "## Task",
         "",
         "Predict, for each health zone with **no confirmed cases up to the forecast date**,",
         "the probability it records its first case(s) in the next 1-2 weeks (spatial",
         "*invasion*). Already-affected zones carry no invasion probability. Forecasts come",
         "from mobility-informed renewal models; the mobility force of infection is",
         "Lambda_i = sum_j W[j,i] sum_k g(k) Y_nc[j,t-k], and the model probability is",
         "P(first case) = 1 - exp(-beta_i * Lambda_i) (cloglog), where log beta_i is an intercept",
         "plus the zone covariates of the specification, if any, learned from realised invasions;",
         "a time-varying form adds a weekly random walk to log beta and/or a zone-level gamma",
         if (.rs_recal)
           paste("frailty. The published probability then carries one fitted recalibration factor",
                 "delta, p -> 1 - (1 - p)^delta (see *Risk scores provided*).")
         else sprintf("frailty. The published probabilities are on the %s scale.", .rs_scale),
         "",
         "The model is fitted in a **Bayesian** framework (brms, or for a time-varying form its",
         "own Stan model through cmdstanr; cloglog with offset(log Lambda)). Every map and",
         "decision product is produced for the skill-selected featured model (files prefixed",
         "`bayes_`); the published panels name it in their manuscript captions rather than in a",
         "figure title. It is scored against the comparator rules listed in the model-comparison",
         "table below, on the identical back-test rounds and zones.",
         "",
         if (!is.null(best_bayes_method))
           # The Bayesian featured model is chosen by best_invasion_model() on the
           # leave-future-out CV composite, exactly like the frequentist one — run_all.R
           # states plainly that the loo weights "NO LONGER pick the featured single model"
           # (they define the stacked ENSEMBLE). This sentence claimed the opposite and was
           # published in invasion_report.md on every run.
           # Plain, unindexed %s only. An earlier edit used "%2$s" for the model and left the
           # trailing "%s" unindexed: R accepts the mix but the unindexed specifier then
           # consumes argument ONE, so the loo-weights slot would have silently rendered
           # primary_method instead. primary_method is dropped from the arguments entirely —
           # there is no second featured model to name.
           sprintf(paste0("Featured model: **%s** (Bayesian), selected by the leave-future-out ",
                          "composite rule recorded in `model_selection.md` (its axes, rounds and ",
                          "tie-breaks are stated there). The loo predictive-stacking weights ",
                          "define the stacked ENSEMBLE only; they do not pick the featured single ",
                          "model.%s"),
                   best_bayes_method,
                   if (!is.null(bayes_weights) && length(bayes_weights))
                     sprintf(" loo-stacking weights: %s.",
                             paste(sprintf("%s=%.2f", names(bayes_weights), bayes_weights), collapse = ", "))
                   else "")
         else "",
         "")

  if (!is.null(eval_tbl)) {
    e1 <- eval_tbl %>% dplyr::filter(horizon == 1)
    # Only co-rank models scored on the COMMON support (comparable base rate); partial-CV
    # models (scored on fewer shared cells) are listed separately.
    has_partial <- "partial_cv" %in% names(e1)
    e1c <- if (has_partial) e1 %>% dplyr::filter(!(partial_cv %in% TRUE)) else e1
    # RANKED BY PER-DATE AUC-PR SKILL OVER ALL ROUNDS, the reported performance. The selection
    # applies the same estimator to its own rounds (post-burn-in by default, a deliberate choice;
    # sel_auc_pr_skill_date), shown beside it, and combines several axes, so this order is not
    # the selection's; the pooled skill is shown for comparison only.
    .sel <- all(c("auc_pr_skill_date", "sel_auc_pr_skill_date") %in% names(e1c))
    e1c <- e1c %>% dplyr::arrange(dplyr::desc(if (.sel) .data$auc_pr_skill_date else .data$auc_pr_skill))
    # INTERVALS OF THE SKILL ITSELF (evaluate_invasion(): each zone-bootstrap replicate's AUC-PR
    # over that replicate's own base rate, and the per-date selection axis on the same
    # resamples). The AUC-PR interval divided by the full-sample base rate is not an interval for
    # the skill, so a table without these columns prints the point estimate alone.
    have_ci  <- all(c("auc_pr_skill_lo", "auc_pr_skill_hi") %in% names(e1c))
    have_sci <- .sel && all(c("auc_pr_skill_date_lo", "auc_pr_skill_date_hi") %in% names(e1c))
    have_ssci <- .sel && all(c("sel_auc_pr_skill_date_lo", "sel_auc_pr_skill_date_hi") %in% names(e1c))
    .num <- function(x, fmt) ifelse(is.finite(x), sprintf(fmt, x), "n/a")
    .ci_cell <- function(x, lo, hi)
      if (is.finite(x) && is.finite(lo) && is.finite(hi)) sprintf("%.1f [%.1f-%.1f]", x, lo, hi)
      else .num(x, "%.1f")
    L <- c(L, "## Model comparison (leave-future-out, 1-week-ahead invasion)",
           "",
           sprintf("All rounds: %d at-risk zone-rounds with %d invasions on the common support (base rate %.2f%%).",
                   e1c$n_atrisk[1] %||% NA, e1c$n_invasions[1] %||% NA, 100 * (e1c$base_rate[1] %||% NA)),
           if (.sel) paste0("Ranked by per-date AUC-PR skill over all rounds, the reported performance. ",
                            "The selection applies the same estimator to its own rounds",
                            if (isTRUE(as.logical(e1c$sel_burnin_applied[1])))
                              paste0(" (the rounds after the ascertainment burn-in only: a deliberate ",
                                     "choice, to keep the discovery phase's ascertainment out of the ",
                                     "model choice)") else "",
                            ", shown in the next column, and combines several axes (see ",
                            "`model_selection.md`), so this order is not the selection's. Pooled AUC-PR ",
                            "skill, mean rank of the invaded zones and the log score are shown for ",
                            "comparison. Bracketed: 90% zone-cluster bootstrap interval of the skill.")
           else "Ranked by pooled AUC-PR skill (AUC-PR / base rate; 1 = a constant forecast, not chance: a random ordering scores above it), with a 90% zone-cluster bootstrap interval of the skill.",
           "The top models overlap within back-test noise: read the ranking as a leading *pack*, not a strict order.",
           "",
           sprintf("| Model | Per-date AUC-PR skill, all rounds%s | Per-date AUC-PR skill, selection rounds%s | Pooled AUC-PR skill, all rounds%s | Mean rank of invaded zones (all rounds) | Log score (raw, all rounds) | Expected/observed (raw) |",
                   if (have_sci) " [90% CI]" else "", if (have_ssci) " [90% CI]" else "",
                   if (have_ci) " [90% CI]" else ""),
           "|---|---|---|---|---|---|---|")
    for (i in seq_len(nrow(e1c))) {
      skill_cell <- if (have_ci) .ci_cell(e1c$auc_pr_skill[i], e1c$auc_pr_skill_lo[i], e1c$auc_pr_skill_hi[i])
                    else .num(e1c$auc_pr_skill[i], "%.1f")
      all_cell <- if (have_sci) .ci_cell(e1c$auc_pr_skill_date[i], e1c$auc_pr_skill_date_lo[i],
                                         e1c$auc_pr_skill_date_hi[i])
                  else if (.sel) .num(e1c$auc_pr_skill_date[i], "%.1f") else "n/a"
      sel_cell <- if (have_ssci) .ci_cell(e1c$sel_auc_pr_skill_date[i], e1c$sel_auc_pr_skill_date_lo[i],
                                          e1c$sel_auc_pr_skill_date_hi[i])
                  else if (.sel) .num(e1c$sel_auc_pr_skill_date[i], "%.1f") else "n/a"
      # Each axis's chance level (a random ordering scores above the base rate on it).
      if (.sel && "auc_pr_skill_date_chance" %in% names(e1c) &&
          is.finite(e1c$auc_pr_skill_date_chance[i]))
        all_cell <- sprintf("%s (chance %.1f)", all_cell, e1c$auc_pr_skill_date_chance[i])
      if (.sel && "sel_auc_pr_skill_date_chance" %in% names(e1c) &&
          is.finite(e1c$sel_auc_pr_skill_date_chance[i]))
        sel_cell <- sprintf("%s (chance %.1f)", sel_cell, e1c$sel_auc_pr_skill_date_chance[i])
      L <- c(L, sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                        e1c$method[i], all_cell, sel_cell,
                        skill_cell, .num(e1c$mean_rank_of_truth[i], "%.1f"),
                        .num(e1c$log_score[i], "%.4f"), .num(e1c$calibration_in_large[i], "%.2f")))
    }
    if (has_partial && any(e1$partial_cv %in% TRUE)) {
      pc <- e1 %>% dplyr::filter(partial_cv %in% TRUE) %>% dplyr::arrange(dplyr::desc(auc_pr_skill))
      L <- c(L, "",
             sprintf("_Partial-CV (scored on fewer shared cells, so not co-ranked): %s._",
                     paste(pc$method, collapse = "; ")))
    }
    .fr <- e1[match(best_bayes_method %||% NA_character_, e1$method), , drop = FALSE]
    L <- c(L, "", sprintf("**Featured operational model:** %s (selected by the composite rule).",
                          best_bayes_method %||% "n/a"),
           "",
           # Derived for THE FEATURED MODEL, not the field maximum: this said "best AUC-ROC" over
           # every method while describing the featured one.
           if (nrow(.fr) && is.finite(.fr$calibration_in_large[1]))
             sprintf(paste0("> **Calibration.** For the featured model the expected number of 1-week ",
                            "invasions over all back-test rounds is %.2f times the observed number on ",
                            "the raw probabilities%s; its AUC-ROC is %.3f. Absolute probabilities are ",
                            "less reliable than the ranking (see FigureS7_calibration)."),
                     .fr$calibration_in_large[1],
                     if ("calibration_in_large_recal" %in% names(.fr) && is.finite(.fr$calibration_in_large_recal[1]))
                       sprintf(" and %.2f times on the recalibrated probabilities as issued%s",
                               .fr$calibration_in_large_recal[1],
                               # ROUNDS WITHOUT A FITTED FACTOR ENTER THAT RATIO RAW (delta = 1). A
                               # burn-in round never carries one -- delta trains only on post-burn-in
                               # rounds settled by the issue date, all of which come later -- so the
                               # fitted rounds are exactly sel_n_folds_delta_fitted of n_folds_scored.
                               if (all(c("n_folds_scored", "sel_n_folds_delta_fitted") %in% names(.fr)) &&
                                   is.finite(.fr$n_folds_scored[1]) && is.finite(.fr$sel_n_folds_delta_fitted[1]))
                                 sprintf(" (%d of the %d rounds carry no fitted factor yet and enter it uncorrected)",
                                         as.integer(.fr$n_folds_scored[1] - .fr$sel_n_folds_delta_fitted[1]),
                                         as.integer(.fr$n_folds_scored[1]))
                               else "")
                     else "",
                     .fr$auc_roc[1])
           else "", "")
  }

  if (!is.null(risk_scores) && "p_case_invasion" %in% names(risk_scores)) {
    rs1 <- risk_scores %>% dplyr::filter(horizon == 1, !was_active_before, !is.na(p_case_invasion))
    tbl <- function(dd, ttl) {
      dd <- dd %>% dplyr::arrange(dplyr::desc(p_case_invasion)) %>% head(10)
      c(sprintf("### %s", ttl), "",
        sprintf("| Zone | Province | P(case), %s | Rel. risk (vs mean at-risk zone nationally) |",
                .rs_scale),
        "|---|---|---|---|",
        vapply(seq_len(nrow(dd)), function(i)
          sprintf("| %s | %s | %.3f | %.1fx |", dd$health_zone[i], dd$province[i],
                  dd$p_case_invasion[i],
                  (if ("rr_nat" %in% names(dd)) dd$rr_nat[i] else NA)), ""), "")
    }
    L <- c(L, "## Highest-risk at-risk zones (next week)", "",
           tbl(rs1 %>% dplyr::filter(province == "Ituri"), "Within Ituri"),
           tbl(rs1, "Nationwide"))
  }
  # TWO PROBABILITY SCALES ARE PUBLISHED, and the models are fitted on the raw one; say which
  # one the tables above carry.
  L <- c(L, "## Risk scores provided", "",
         if (.rs_recal)
           c("- **P(case)** — probability of >=1 CONFIRMED case in the forecast window, on the",
             "  recalibrated scale as deployed: the model's raw probability passed through the pooled",
             "  factor, p -> 1 - (1 - p)^delta, with delta fitted on the settled post-burn-in 1-week",
             "  back-test rows of the shared support (`diagnostics/invasion_recalibration.csv`). The models are fitted, and",
             "  their discrimination scored, on the raw probability; the uncorrected twin of the risk",
             "  table is `reports/raw/bayes_risk_scores_all_zones.csv` (published copy in `key_outputs/raw/`).")
         else
           sprintf("- **P(case)** — probability of >=1 CONFIRMED case in the forecast window, on the %s scale.",
                   .rs_scale),
         "- **Relative risk** — multiple of the mean at-risk zone nationwide (both tables use the",
         "  national reference).",
         "- Every zone, both horizons, with the posterior 90% interval: `key_outputs/bayes_risk_scores_all_zones.csv`.", "")

  # ---- Findings / strengths / weaknesses / improvements ------------------------------------
  # Every quantitative statement is THE FEATURED MODEL's (the composite selection's pick), read
  # from its own row; none is a maximum over the grid attributed to it.
  if (!is.null(eval_tbl)) {
    ez <- eval_tbl %>% dplyr::filter(horizon == 1)
    if ("partial_cv" %in% names(ez)) ez <- ez %>% dplyr::filter(!(partial_cv %in% TRUE))
    fz <- ez[match(best_bayes_method %||% NA_character_, ez$method), , drop = FALSE]
    fz <- fz[!is.na(fz$method), , drop = FALSE]
    .tk  <- invasion_topk_cols(ez)
    rcol <- if (is.null(.tk)) NA_character_ else .tk$recall
    kk   <- if (is.null(.tk)) NA_integer_ else .tk$k
    f_rec <- if (nrow(fz) && !is.na(rcol)) fz[[rcol]][1] else NA_real_
    # Random watch-list: K of each round's at-risk zones, i.e. K * (rounds scored) / (pooled
    # at-risk zone-rounds) of the invasions in expectation.
    .nf_r <- if (nrow(fz) && "n_folds_scored" %in% names(fz)) fz$n_folds_scored[1]
             else if (nrow(fz) && "n_folds" %in% names(fz)) fz$n_folds[1] else NA_real_
    rand_rec <- if (nrow(fz) && "n_atrisk" %in% names(fz)) kk * .nf_r / fz$n_atrisk[1] else NA_real_
    L <- c(L, "## Findings", "",
      if (nrow(fz))
        sprintf("- **Discrimination (featured model `%s`, one week, all rounds).** AUC-ROC %.3f; pooled AUC-PR skill %.1f times the base rate%s%s.",
                best_bayes_method, fz$auc_roc[1], fz$auc_pr_skill[1],
                if (all(c("auc_pr_skill_lo", "auc_pr_skill_hi") %in% names(fz)) &&
                    is.finite(fz$auc_pr_skill_lo[1]) && is.finite(fz$auc_pr_skill_hi[1]))
                  sprintf(" (90%% zone-cluster bootstrap interval %.1f-%.1f)", fz$auc_pr_skill_lo[1], fz$auc_pr_skill_hi[1])
                else "",
                if (all(c("auc_pr_skill_date", "sel_auc_pr_skill_date") %in% names(fz)) &&
                    is.finite(fz$auc_pr_skill_date[1]))
                  sprintf("; per-date AUC-PR skill over all rounds %.1f%s (%.1f on the %s the selection reads)",
                          fz$auc_pr_skill_date[1],
                          if ("auc_pr_skill_date_chance" %in% names(fz) && is.finite(fz$auc_pr_skill_date_chance[1]))
                            sprintf(", against %.1f for a random ordering", fz$auc_pr_skill_date_chance[1]) else "",
                          fz$sel_auc_pr_skill_date[1],
                          if (isTRUE(as.logical(fz$sel_burnin_applied[1]))) "post-burn-in rounds" else "rounds")
                else "")
      else "- The featured model is absent from the evaluation table.",
      # recall_at_K is the MEAN OVER ROUNDS (with an invasion) of each round's recall
      # (.ranking_metrics(), 16_invasion_eval.R), not the pooled share of all invasions, and the
      # random reference is likewise per round; the sentence says so.
      if (nrow(fz) && is.finite(f_rec) && is.finite(rand_rec))
        sprintf("- **Targeting.** Averaged over the back-test rounds with an invasion, a top-%d watch-list from the featured model contains %.0f%% of that round's newly invaded zones, against %.0f%% expected for a random list of the same size (%.1f times).",
                kk, 100 * f_rec, 100 * rand_rec, f_rec / rand_rec)
      else "",
      "- **Uncertainty.** Every zone's posterior 90% interval is in `bayes_risk_scores_all_zones.csv`; skill estimates carry zone-cluster bootstrap intervals.",
      "")
    L <- c(L, "## Strengths", "",
      "- **Out-of-sample evaluation as it would have run.** Each back-test round is trained on the data as they stood at its forecast origin: the archived line-list snapshot of that date, nowcast with the truncation estimable on that date. Rounds before the archive, or before it supports that truncation fit, use the final line list censored to the origin instead, and every round records which. Outcomes are first confirmed onset weeks on the final data; all models are compared on one shared set of zone-rounds.",
      # The back-tested ensembles are the mean/median combinations of a fixed member set
      # (append_ensembles(), "Bayes-ens-*"); the loo-stacked ensemble is a LIVE forecast only and
      # has no back-test row.
      "- **A grid of structural assumptions** (mobility kernels, covariate sets) scored on the same rounds, together with generation-time sensitivity refits and time-varying forms of the selected specification, and mean/median ensembles of a fixed member set.",
      "- **Comparators on identical rounds**, including rules that use no case data and rules that use no phone data (the model-comparison table above; Figure 2).",
      "")
    # THE ONSET-IMPUTATION SHARE, BY PROVENANCE (01_data_prep.R): the genuine line-list
    # missingness, and the share of all modelled confirmed records once the sitrep top-up rows
    # (no onset by construction) are counted. The mechanism is the loader's own description.
    .pi_ll <- get0(".PCT_ONSET_IMPUTED_LINELIST", ifnotfound = NA_real_)
    .pi_all <- get0(".PCT_ONSET_IMPUTED", ifnotfound = NA_real_)
    .n_sr <- get0(".N_SITREP_APPENDED", ifnotfound = NA_integer_)
    .desc <- get0(".ONSET_IMPUTE_DESC", ifnotfound = NA_character_)
    .imp_line <- if (!identical(get0(".ONSET_IMPUTE_MODE", ifnotfound = NA_character_), "impute"))
      sprintf("- **Onset.** %s.", if (is.na(.desc)) "The onset handling of this run is not recorded in this session" else .desc)
    else sprintf(paste0("- **Onset imputation.** %s; for these, onset = %s, widening the donor window by one ",
                        "week on each side at a time until it holds at least %d donors. One completion per run (multiple ",
                        "imputation is run separately), so imputation-dependent intervals are conditionally narrow."),
                 if (is.finite(.pi_ll) && is.finite(.pi_all) && is.finite(.n_sr))
                   sprintf(paste0("%.0f%% of line-list confirmed cases have no usable recorded onset, and the %d ",
                                  "sitrep top-up rows carry none by construction (%.0f%% of all modelled confirmed ",
                                  "records)"), .pi_ll, as.integer(.n_sr), .pi_all)
                 else "The share of confirmed cases without a usable onset is not available in this session",
                 .desc, as.integer(get0("ONSET_DONOR_MIN", ifnotfound = NA_integer_)))
    L <- c(L, "## Weaknesses & limitations", "",
      "- **Absolute calibration** is weaker than the ranking (calibration line above; FigureS7_calibration): prioritise on the ranking and relative risk.",
      sprintf("- **Few events.** %s one-week invasions on the common support; intervals are wide and the leaderboard is a pack, not a strict order.",
              if (nrow(fz)) as.character(fz$n_invasions[1]) else "Few"),
      "- **No BDBV-specific generation time exists;** a Zaire-ebolavirus serial-interval proxy is used (short/medium/long sensitivity).",
      .imp_line,
      "- **Two-stage nowcast.** Truncation-corrected counts enter the models as a fixed input (a deterministic correction from the fitted onset-to-appearance truncation); its uncertainty is not propagated.",
      "- **R(t)** drives the two-week horizon's projection of source incidence and is least certain at the most recent weeks, which is where it is used.",
      "- **Mobility kernels** are built once per run from the data available at the analysis date, so the earliest back-test rounds use mobility data released after their origin.",
      "- **Duplicated health-zone names** across provinces (Bili, Lubunga) carry an arbitrary province label in free text (maps disambiguate via a name+province join).",
      "")
    L <- c(L, "## Where to improve", "",
      "- **A joint delay + transmission model** for the recent tail, propagating the nowcast's uncertainty, instead of a two-step nowcast -> forecast.",
      "- **Mobility kernels dated to each back-test origin.**",
      "")
  }
  writeLines(L, out_path)
  message("[report] invasion report -> ", out_path)
  invisible(out_path)
}

# ---------------------------------------------------------------------------
# Reliability diagram — predicted vs observed invasion frequency
# ---------------------------------------------------------------------------

#' Rank-based AUC-ROC (Mann-Whitney): P(random invaded ranked above random not).
# NOTE: this function is ALSO defined in 16_invasion_eval.R. run_all.R sources 16 first and
# 17 last, so an UNCONDITIONAL definition here used to win for the whole session — including
# for evaluate_invasion() in 16 and spatiotemporal_skill() in 19, which both call `.auc_roc`.
# The two bodies used to DIFFER: this one omitted the finite/NA filter, so on any input with
# a non-finite probability or a missing outcome it ranked NAs last (rank()'s default
# na.last = TRUE gives them the HIGHEST ranks) while n1/n0 counted a different set — e.g.
# p = (0.9, 0.2, NA, 0.7, 0.1), y = (1,0,0,1,0) returned 0.667 against the correct 1.000.
# Both callers happen to pre-filter today, so nothing was wrong in practice, but the weaker
# definition was silently shadowing the safe one purely by source order. The bodies are now
# IDENTICAL, so the order cannot matter; keep them in step if either is edited.
# DO NOT SHADOW THE EVALUATOR. This is a VISUALISATION module, but it is sourced AFTER
# 16_invasion_eval.R in run_all.R, so an unconditional definition here wins for the whole
# session — including inside evaluate_invasion() (16) and spatiotemporal_skill() (19), which
# write invasion_evaluation.csv and skill_over_time.csv. The two bodies are identical today;
# the moment one is edited, a figure file would silently redefine the published estimator.
# Defined only as a STANDALONE fallback, so 16's version always governs when both are loaded.
if (!exists(".auc_roc", mode = "function")) {
  .auc_roc <- function(p, y) {
    y <- as.integer(y); ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
    n1 <- sum(y == 1); n0 <- sum(y == 0)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- rank(p, ties.method = "average")
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
}

#' Calibration / reliability — the Fig-3A analogue of Kraemer & Cauchemez 2017
#' (Lancet ID): predicted invasion-probability bins vs observed proportion invaded,
#' with the y=x diagonal, annotated with AUC-ROC and AUC-PR skill, plus their
#' signature high-risk callout (the small share of zone-rounds scored high risk and
#' the much higher observed invasion rate among them). Pooled over the at-risk LFO rows of the
#' shared (fold x zone) support; points below the diagonal = over-prediction.
#' @param support_cells the shared support ("fold\\rzone" ids); NULL resolves it from
#'   `lfo_results` with invasion_common_cells() (16_invasion_eval.R), so the AUC annotations are
#'   on the rows evaluate_invasion() scored.
plot_reliability <- function(lfo_results, horizon = 1L, methods = NULL, save = TRUE,
                             panel_tag = NULL, file = "reliability", support_cells = NULL) {
  # Bins balance EXPECTED events (sum of p), never the outcomes: reliability_breaks()
  # (23_prospective_eval.R) with the published binning (RELIABILITY_N_BINS / _BINNING).
  if (!exists("reliability_breaks", mode = "function"))
    stop("[reliability] reliability_breaks() (23_prospective_eval.R) is not loaded.", call. = FALSE)
  if (is.null(support_cells) && exists("invasion_common_cells", mode = "function"))
    support_cells <- tryCatch(invasion_common_cells(lfo_results, horizon), error = function(e) NULL)
  d <- lfo_results %>% dplyr::filter(horizon == !!horizon, is.finite(p_invasion))
  if ("was_active_before" %in% names(d))
    d <- d %>% dplyr::filter(!(was_active_before %in% TRUE))
  if (length(support_cells))
    d <- d[paste(d$fold_id, d$health_zone, sep = "\r") %in% support_cells, , drop = FALSE]
  else warning("[reliability] no shared support resolved; each method is drawn on its own rows.",
               call. = FALSE)
  .raw <- .lfo_is_raw(d)
  # A RELIABILITY PLOT NEEDS A PROBABILITY ON ITS X-AXIS. Rank-only comparators
  # (prob_calibrated = FALSE: Distance-B1's unscaled neighbour-count mean, Adjacency-B7's
  # 1/(1 + travel-time minutes)) carry an ordering, not a probability, so plotting them here
  # would draw a calibration curve for an arbitrary monotone transform. They are dropped
  # BEFORE the default-method selection below, which picks by AUC-ROC — a rank metric these
  # can score well on, so without this they could be selected INTO the panel precisely because
  # they discriminate. Absent column = calibrated (every Bayesian model).
  if ("prob_calibrated" %in% names(d))
    d <- d %>% dplyr::filter(!(prob_calibrated %in% FALSE))
  if (is.null(methods)) {
    # The frequentist "^Renewal" family is gone, so this now resolves to every remaining
    # (calibrated) method; the top-8 trim below keeps the panel legible.
    methods <- unique(d$method)
  } else {
    # An explicit request for a rank-only method is a caller error, not something to draw.
    .drop <- setdiff(methods, unique(d$method))
    if (length(.drop))
      warning(sprintf(paste0("[reliability] %s cannot appear on a reliability plot (no ",
                             "calibrated probability) and %s dropped."),
                      paste(.drop, collapse = ", "),
                      if (length(.drop) > 1) "were" else "was"), call. = FALSE)
    methods <- intersect(methods, unique(d$method))
    if (!length(methods)) return(invisible(NULL))
  }
  # Keep the plot legible and within the categorical palette: if many renewal variants are
  # eligible, show the top methods by pooled discrimination (AUC-ROC). The AUC annotation is
  # built over `methods` below, so it stays in sync with the series actually drawn.
  if (length(methods) > 8) {
    .roc <- vapply(methods, function(mm) { dm <- d[d$method == mm, ]
      suppressWarnings(.auc_roc(dm$p_invasion, dm$is_new_invasion)) }, numeric(1))
    methods <- names(sort(.roc, decreasing = TRUE))[seq_len(8)]
  }
  # Edges per method, from that method's own forecasts alone.
  rel <- d %>% dplyr::filter(method %in% methods) %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(bin = cut(p_invasion,
                            reliability_breaks(p_invasion, is_new_invasion, RELIABILITY_N_BINS,
                                               RELIABILITY_BINNING),
                            include.lowest = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(method, bin) %>%
    dplyr::summarise(pm = mean(p_invasion), ym = mean(is_new_invasion),
                     n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(is.finite(pm), n >= 3)
  if (nrow(rel) == 0) return(invisible(NULL))
  # AUC-ROC + AUC-PR skill per featured method (pooled, unbinned).
  ap_fn <- get0(".auc_pr")
  auc <- d %>% dplyr::filter(method %in% methods) %>% dplyr::group_by(method) %>%
    dplyr::summarise(
      roc = .auc_roc(p_invasion, is_new_invasion),
      brate = mean(is_new_invasion, na.rm = TRUE),
      prsk = if (is.function(ap_fn) && brate > 0) ap_fn(p_invasion, is_new_invasion) / brate else NA_real_,
      .groups = "drop") %>% dplyr::arrange(dplyr::desc(roc))
  auc_lab <- paste(sprintf("%s: AUC %.2f, AUC-PR skill %.1fx", auc$method, auc$roc, auc$prsk),
                   collapse = "\n")
  # High-risk callout (their "0.31% of district-weeks >25%/wk -> 77%/month" line). The
  # probability is CUMULATIVE over the horizon's window, so the threshold is per window.
  hi <- d %>% dplyr::filter(method == methods[1])
  hi_share <- mean(hi$p_invasion > 0.15, na.rm = TRUE)
  hi_obs   <- if (any(hi$p_invasion > 0.15)) mean(hi$is_new_invasion[hi$p_invasion > 0.15]) else NA_real_
  .win <- if (as.integer(horizon) > 1L) sprintf("over %d weeks", as.integer(horizon)) else "for the week"
  callout <- if (is.finite(hi_obs))
    sprintf("%s: the %.1f%% of at-risk zone-rounds scored >15%% %s had a %.0f%% observed invasion rate.",
            methods[1], 100 * hi_share, .win, 100 * hi_obs)
    else "Pooled over folds; few events, so absolute levels are indicative."
  # Shared x/y maximum so the y = x diagonal is a true 45-degree reference and
  # over-/under-confidence read correctly (predicted and observed on the SAME scale).
  .axmax <- max(c(rel$pm, rel$ym), 0.02, na.rm = TRUE) * 1.05
  p <- ggplot(rel, aes(pm, ym)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_line(aes(colour = method), linewidth = 0.8) +
    geom_point(aes(colour = method, size = n)) +
    # linewidth, not label.size: ggplot2 4 ignores label.size (with a warning) and draws a border.
    annotate("label", x = 0, y = Inf, hjust = 0, vjust = 1.05, label = auc_lab,
             size = 2.9, linewidth = 0, fill = "white", alpha = 0.7, colour = "grey20") +
    scale_colour_manual(values = setNames(
      if (length(methods) <= length(OKABE_ITO)) OKABE_ITO[seq_along(methods)]
      else grDevices::hcl.colors(length(methods), "Dark 3"), methods), name = NULL) +
    scale_size_area(max_size = 7, name = "zone-rounds") +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, .axmax)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, .axmax)) +
    coord_equal(xlim = c(0, .axmax), ylim = c(0, .axmax)) +
    labs(title = sprintf("%sCalibration: predicted vs observed invasion (h = %dw)",
                         if (!is.null(panel_tag)) paste0(panel_tag, "  ") else "", horizon),
         subtitle = paste0("On the diagonal = calibrated probabilities; below = over-prediction.\n",
                           if (identical(get0("RELIABILITY_BINNING", ifnotfound = "expected"), "width"))
                             "Equal-width bins." else "Bins balanced on expected events."),
         x = if (.raw) "Predicted invasion probability (raw)"
             else "Predicted invasion probability\n(recalibrated as issued; raw before a factor was fitted)",
         y = "Observed invasion frequency",
         caption = paste(strwrap(callout, width = 85), collapse = "\n")) +
    theme_inv(11) + theme(plot.caption = element_text(colour = "grey30", size = 8.5, hjust = 0, face = "italic"))
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, sprintf("%s_h%d.pdf", file, horizon)), 7, 6)
  invisible(p)
}

#' The combined "Figure 3" briefing panel, after Kraemer & Cauchemez 2017 (Lancet
#' Inf Dis), the published ancestor of this invasion model: (A) calibration with AUC
#' and (B) real-time prioritisation vs random targeting, side by side — the single
#' most decision-legible summary of how good the model is at predicting invasion.
#' @param support_cells the shared (fold x zone) support resolved from the SCORED table.
#'   run_all.R passes `lfo_fig3` here (it carries the appended rank-only baselines), so
#'   without this panel B would derive its own, smaller support and its random-targeting
#'   reference would not match the published detection curve.
plot_paper_figure3 <- function(lfo_results, methods, horizon = 1L, save = TRUE,
                               file = "model_performance_figure3", support_cells = NULL) {
  a <- tryCatch(plot_reliability(lfo_results, horizon = horizon, methods = methods[1],
                                 save = FALSE, panel_tag = "(A)", support_cells = support_cells),
                error = function(e) NULL)
  b <- if (exists("plot_detection_curve"))
    tryCatch(plot_detection_curve(lfo_results, methods, horizon = horizon, save = FALSE,
                                  support_cells = support_cells),
             error = function(e) NULL) else NULL
  if (is.null(a) && is.null(b)) return(invisible(NULL))
  if (is.null(a)) return(invisible(b))
  if (is.null(b)) return(invisible(a))
  b <- b + labs(title = "(B)  Prioritisation: invasions caught vs random targeting")
  out <- a + b + patchwork::plot_layout(widths = c(1, 1.25)) +
    patchwork::plot_annotation(
      title = sprintf("How good is the invasion model? (h = %dw)", horizon),
      subtitle = "(A) calibration + AUC and (B) real-time prioritisation vs random — after Kraemer & Cauchemez 2017, Lancet Inf Dis, Fig 3",
      theme = ggplot2::theme(plot.title = element_text(face = "bold", size = 13),
                             plot.subtitle = element_text(size = 9, colour = "grey35")))
  if (save) .save(out, file.path(OUT_REPORTS, sprintf("%s_h%d.pdf", file, horizon)), 14, 6)
  invisible(out)
}

message("[viz] 17_invasion_viz.R loaded — legible invasion figure suite.")

# ---------------------------------------------------------------------------
# Delta stability: is the calibration offset constant over time?
# ---------------------------------------------------------------------------

#' Visualise the delta-stability diagnostic. READS the published tables; computes nothing.
#'
#' delta is an intercept shift on the cloglog scale and beta_0 is the intercept, so the
#' per-fold delta IS a time-varying beta_0 fitted post hoc. This panel is the evidence for two
#' design choices at once: whether the deployment delta should pool all history or use a
#' recent window, and whether a time-varying beta_0 is worth fitting at all.
#'
#' Panel A forest: each fold's INDEPENDENTLY fitted delta with its normal interval on log delta
#' (at the level 16b builds it, read back from the table), against the DerSimonian-Laird
#' random-effects mean. Panel B: I-squared for every method, i.e. the share of the
#' between-fold spread that is NOT sampling noise.
#'
#' EVERY NUMBER IS READ from invasion_delta_stability{,_summary}.csv, which run_all.R writes
#' from invasion_delta_stability() (16b). Nothing here re-fits, re-pools or re-tests.
#'
#' @param per_fold,summary_tbl the two published tables (or paths to them).
#' @param method which method to draw in panel A; defaults to the featured model
#'   (.featured_for_figure()), and the panel is skipped with a warning if the per-fold table
#'   does not carry it at this horizon.
plot_delta_stability <- function(per_fold, summary_tbl, horizon = 1L, method = NULL,
                                 save = TRUE, file = "delta_stability") {
  rd <- function(x) if (!is.character(x)) x
                    else if (file.exists(x)) readr::read_csv(x, show_col_types = FALSE) else NULL
  pf <- rd(per_fold); sm <- rd(summary_tbl)
  if (is.null(pf) || !nrow(pf) || is.null(sm) || !nrow(sm)) return(invisible(NULL))
  pf <- pf %>% dplyr::filter(horizon == !!horizon, is.finite(delta), is.finite(se_log))
  sm <- sm %>% dplyr::filter(horizon == !!horizon)
  if (!nrow(pf) || !nrow(sm)) return(invisible(NULL))
  if (is.null(method)) method <- .featured_for_figure(NULL)
  if (!method %in% pf$method) {
    warning(sprintf("[delta-stability] %s is not in the per-fold table at h=%s; panel skipped.",
                    method, horizon), call. = FALSE)
    return(invisible(NULL))
  }
  a <- pf %>% dplyr::filter(method == !!method)
  srow <- sm %>% dplyr::filter(method == !!method)
  if (!nrow(a) || !nrow(srow)) return(invisible(NULL))
  a$lab <- if (all(!is.na(a$cutoff))) format(lfo_origin(a$cutoff), "%d %b") else as.character(a$fold_id)
  a$lab <- factor(a$lab, levels = rev(a$lab[order(a$fold_id)]))
  # THE INTERVAL LEVEL IS READ FROM THE TABLE, not assumed: 16b publishes lo/hi = delta x/
  # exp(z * se_log), so z = log(hi / delta) / se_log, and it builds the pooled interval with the
  # same normal multiplier. (The tables carry no level column; the label follows them if the
  # multiplier is changed there.)
  .z <- stats::median(log(a$hi / a$delta) / a$se_log, na.rm = TRUE)
  .lev <- if (is.finite(.z) && .z > 0) sprintf("%.0f%%", 100 * (2 * stats::pnorm(.z) - 1)) else "level unknown"

  pA <- ggplot(a, aes(x = .data$delta, y = .data$lab)) +
    annotate("rect", xmin = srow$pooled_lo[1], xmax = srow$pooled_hi[1],
             ymin = -Inf, ymax = Inf, fill = "grey88") +
    geom_vline(xintercept = srow$pooled_delta[1], colour = "grey35", linewidth = 0.6) +
    geom_vline(xintercept = 1, linetype = "dotted", colour = "grey55") +
    geom_linerange(aes(xmin = .data$lo, xmax = .data$hi), colour = "grey30", linewidth = 0.6) +
    geom_point(aes(size = .data$n_events), colour = "grey10") +
    scale_size_continuous(range = c(1.2, 3.4), name = "events") +
    scale_x_log10() +
    # The band is the DerSimonian-Laird RANDOM-EFFECTS MEAN of the per-fold log deltas, the
    # quantity the heterogeneity test is centred on; it is NOT the deployed factor, which is
    # fitted by penalised ML on the pooled 1-week rows (invasion_recalibration.csv).
    labs(x = sprintf("fold-specific delta with its %s interval  (log scale; 1 = already calibrated)", .lev),
         y = NULL,
         subtitle = sprintf(paste0("%s\nband = DerSimonian-Laird random-effects mean delta %.2f ",
                                   "[%s: %.2f, %.2f]\n(not the deployed factor)"),
                            method, srow$pooled_delta[1], .lev, srow$pooled_lo[1],
                            srow$pooled_hi[1])) +
    theme_minimal(base_size = 10)

  sb <- sm %>% dplyr::filter(is.finite(.data$I2)) %>%
    dplyr::mutate(method = forcats::fct_reorder(.data$method, .data$I2))
  pB <- ggplot(sb, aes(x = 100 * .data$I2, y = .data$method)) +
    geom_segment(aes(x = 0, xend = 100 * .data$I2, yend = .data$method),
                 colour = "grey80", linewidth = 0.4) +
    geom_point(colour = "grey10", size = 1.8) +
    geom_vline(xintercept = 25, linetype = "dashed", colour = "grey55") +
    scale_x_continuous(limits = c(0, 100)) +
    labs(x = expression(paste(I^2, "  (% of between-fold spread that is NOT sampling noise)")),
         y = NULL,
         subtitle = sprintf("all %d methods; dashed line = 25%%, the conventional floor for 'some' heterogeneity",
                            nrow(sb))) +
    theme_minimal(base_size = 9)

  # WHICH POOL. The deployed factor is trained on the 1-WEEK rows of the settled post-burn-in
  # rounds (shared support), so only the h = 1 panel shows its training pool; at h >= 2 the
  # rows are that horizon's own settled rounds and the panel asks whether the h-week offset is
  # itself stable.
  .pool <- if (as.integer(horizon) == 1L)
    paste0("Folds: post-burn-in rounds with settled outcomes, on the shared support -- the ",
           "1-week rows the deployed factor is trained on.")
  else sprintf(paste0("Folds: post-burn-in rounds whose %d-week outcomes are settled, on the shared ",
                      "support. The deployed factor is trained on the 1-week rows only, so this panel ",
                      "asks whether the %d-week offset is itself stable; it is not the deployment pool."),
               as.integer(horizon), as.integer(horizon))
  # The single-fold precision: exp(median se_log) is ONE standard error on the log scale; the
  # drawn per-fold interval is exp(z * se_log) wide on each side.
  cap <- sprintf(paste0("Each fold's delta is fitted on that fold ALONE, so the spread is not an ",
                        "expanding window converging. Cochran Q p = %.3f on %d df; tau^2 = %.4f; ",
                        "I^2 = %.0f%%; trend in log delta per fold = %+.3f (p = %.2f). A single ",
                        "fold's standard error on log delta (median over folds) is a factor of x/%.2f, ",
                        "so its %s interval spans x/%.2f. Heterogeneity at p<0.05 in %d of the %d ",
                        "methods with a test. %s"),
                 srow$Q_p[1], srow$Q_df[1], srow$tau2[1], 100 * srow$I2[1],
                 srow$trend_slope[1], srow$trend_p[1], exp(srow$median_se_log[1]),
                 .lev, exp(.z * srow$median_se_log[1]),
                 sum(sm$Q_p < 0.05, na.rm = TRUE), sum(is.finite(sm$Q_p)), .pool)
  p <- (pA | pB) + patchwork::plot_layout(widths = c(1, 1.15)) +
    patchwork::plot_annotation(
      title = sprintf("Is the calibration offset stable over time?  (h = %s)", horizon),
      caption = paste(strwrap(cap, width = 150), collapse = "\n"),
      theme = ggplot2::theme(
        plot.title   = ggplot2::element_text(size = 12, face = "bold"),
        plot.caption = ggplot2::element_text(size = 7.5, colour = "grey35", hjust = 0)))
  # Written beside the other diagnostics (OUT_DIAGNOSTICS); "delta_stability" is on FIGURE_KEEP,
  # so .save()'s gate keeps delta_stability_h<N>.pdf.
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, sprintf("%s_h%d.pdf", file, horizon)), w = 12, h = 5.2)
  invisible(p)
}
