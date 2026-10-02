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
  # Use the base pdf device (NOT cairo_pdf): cairo_pdf fails to load its DLL in
  # the headless/batch R environment and then silently writes nothing, leaving a
  # stale file. capabilities("cairo") reports COMPILE-time support and is TRUE
  # even when the runtime device fails, so it cannot be trusted. The base pdf
  # device always works. Delete any stale target first, then verify the write.
  if (file.exists(path)) suppressWarnings(file.remove(path))
  tryCatch(ggsave(path, p, width = w, height = h, device = "pdf",
                  limitsize = FALSE),
           error = function(e) warning("[viz] ggsave error ", basename(path),
                                       ": ", conditionMessage(e)))
  if (file.exists(path)) message("[viz] saved -> ", path)
  else warning("[viz] FAILED to write ", basename(path))
  invisible(p)
}

# ---------------------------------------------------------------------------
# 1. Legible epidemic curve (windowed, top-N zones + Other, weekly)
# ---------------------------------------------------------------------------

#' @param zone_week_nc nowcast-corrected zone-week tibble.
#' @param province_map tibble(nom, province) to restrict the named series to Ituri.
plot_epidemic_curve_legible <- function(zone_week_nc, province_map = NULL,
                                        top_n = 6L,
                                        window_start = get0("OUTBREAK_START", ifnotfound = Sys.Date() - 120) - 14,
                                        save = TRUE) {
  d <- zone_week_nc %>%
    dplyr::filter(week_start >= window_start) %>%
    dplyr::mutate(conf = ifelse(is.na(confirmed_nc), 0, confirmed_nc))
  # rank named zones by total cases (prefer Ituri context if map available)
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

  base_layer <- function(dat, title, legend = FALSE) {
    ggplot(dat) +
      geom_sf(aes(fill = p), colour = "grey80", linewidth = 0.12) +
      geom_sf(data = dplyr::filter(dat, affected %in% TRUE),
              fill = "grey55", colour = "grey80", linewidth = 0.12) +
      scale_fill_viridis_c(option = "plasma", name = "P(first case)",
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
#' targeting view: at-risk zones are ranked nationally by P(first case), coloured by rank (bright =
#' top), and the top-K are labelled with their rank number. Already-affected zones are unranked (grey).
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
                     p = .data[[pcol]], affected = was_active_before %in% TRUE)
  # RANK invasion probability across all at-risk zones (1 = highest); affected zones unranked.
  # Fill CONTINUOUSLY over the FULL rank range (no cap, no binning): rank is UNIFORMLY spread
  # 1..N, unlike the highly-skewed probability (where 400+ zones pile at ~0 and look identical),
  # so the colour varies smoothly across EVERY zone and all 519 are distinguishable — including
  # the long tail beyond the top 50. Also carry a percentile for an interpretable legend.
  # THIS RANK IS THE PUBLISHED ONE — checked, not assumed. compute_risk_scores()
  # (15_workhorse.R) builds rr_nat_rank as min_rank(desc(mu)); this uses
  # rank(-p, ties.method = "min"). p = 1 - exp(-mu) is strictly increasing in mu, so the two
  # orderings and their tie groups are identical and the integers printed on this map equal
  # rr_nat_rank in bayes_risk_scores_all_zones.csv. Do not "align" them by switching either
  # side to a different tie method: "min" is what the published column uses.
  #
  # NOT the same quantity as rank_med/rank_lo/rank_hi in that file, which are the POSTERIOR
  # MEDIAN rank and its 90% CrI (21_bayesian_renewal.R). Those describe how the ranking varies
  # across draws; these labels are the point-estimate ranking. Both are published, they answer
  # different questions, and the map deliberately shows the point estimate.
  n_ar <- sum(!(rd$affected) & is.finite(rd$p))
  rd <- rd %>% dplyr::mutate(
    rank_p = ifelse(affected | !is.finite(p), NA_real_, rank(-p, ties.method = "min", na.last = "keep")),
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
  sub <- sprintf("At-risk zones shaded by national invasion-risk percentile (rank of P(first case); bright = top); numbers label the top %d; dark grey = already affected.", top_k)
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
      labs(title = ttl, x = "Next-week probability of first case(s)", y = NULL) +
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
    grepl("^Ensemble", m)                 ~ "Ensemble",
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
  d <- d %>% dplyr::mutate(family = .invasion_family(method),
                           method = forcats::fct_reorder(method, auc_pr_skill))
  .ttl <- if (!is.null(label)) sprintf("%s invasion discrimination (h=%dw)", label, horizon)
          else sprintf("Invasion discrimination (h=%dw)", horizon)
  p1 <- ggplot(d, aes(auc_pr_skill, method, colour = family)) +
    geom_segment(aes(x = 1, xend = auc_pr_skill, yend = method), colour = "grey80", linewidth = 0.5) +
    { if (all(c("auc_pr_lo","auc_pr_hi") %in% names(d)))
        geom_errorbar(aes(xmin = pmax(auc_pr_lo,0)/base_rate, xmax = auc_pr_hi/base_rate),
                      orientation = "y", width = 0, linewidth = 0.5, na.rm = TRUE) } +
    geom_point(size = 3) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = .FAMILY_COLOURS, name = NULL) +
    labs(title = .ttl,
         subtitle = "AUC-PR skill = AUC-PR ÷ base rate; 1 = no skill (dashed). Higher = better.",
         x = "AUC-PR skill (× base rate)", y = NULL) +
    theme_inv(11)
  p2 <- ggplot(d, aes(mean_rank_of_truth, method, colour = family)) +
    geom_segment(aes(x = max(d$mean_rank_of_truth, na.rm=TRUE), xend = mean_rank_of_truth, yend = method),
                 colour = "grey80", linewidth = 0.5) +
    geom_point(size = 3) +
    scale_colour_manual(values = .FAMILY_COLOURS, name = NULL, guide = "none") +
    labs(title = "Mean rank of invaded zones",
         subtitle = "Rank of EVERY invaded zone (tie-averaged), meaned per fold then over folds; lower = nearer the top",
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
#' (recall@10, share of the true invasions caught by a 10-zone watch-list) — with every model
#' drawn, coloured by family, the naive-baseline band shaded, and the no-skill / random-watch
#' -list references marked. Makes the renewal/Bayesian-over-baseline gap explicit.
plot_model_vs_baseline <- function(eval_tbl, horizon = 1L, save = TRUE) {
  d <- eval_tbl %>% dplyr::filter(horizon == !!horizon)
  if (!nrow(d)) return(invisible(NULL))
  # partial-CV models are not comparably scored — drop them from the head-to-head.
  if ("partial_cv" %in% names(d)) d <- d %>% dplyr::filter(!(partial_cv %in% TRUE))
  d <- d %>% dplyr::mutate(family = .invasion_family(method),
                           is_baseline = family == "Baseline (naive)")
  rcol <- if ("recall_at_10" %in% names(d)) "recall_at_10" else "recall_at_5"
  kk   <- if (rcol == "recall_at_10") 10L else 5L
  # recall@K is a PER-FOLD-averaged metric, so the random watch-list baseline is K / (at-risk
  # zones PER FOLD) = K * n_folds / n_atrisk (n_atrisk is pooled across folds); dividing by the
  # pooled n_atrisk would understate the random reference ~n_folds-fold.
  # n_folds_scored, NOT n_folds. Since the common-support fix n_folds is each method's NATIVE
  # reach while n_folds_scored is what its metrics were actually computed on; taking the median
  # of the former put the random reference on a fold count no row was scored at whenever the
  # grid's composition shifted the median (it is 11 either way on this run, by majority).
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
         subtitle = "AUC-PR skill = AUC-PR / base rate; 1 = no skill (dashed). Shaded/dotted = best naive baseline.",
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
  dd <- d %>% dplyr::filter(health_zone %in% invaded_zones) %>%
    dplyr::mutate(cutoff = lfo_origin(cutoff),   # label rounds by forecast origin
                  health_zone = factor(health_zone, levels = rev(sort(invaded_zones))))
  p <- ggplot(dd, aes(cutoff, health_zone, fill = p_invasion)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_point(data = dplyr::filter(dd, is_new_invasion == 1),
               shape = 4, size = 3, stroke = 1.1, colour = "black") +
    scale_fill_viridis_c(option = "plasma", name = "Predicted P", labels = percent_format(accuracy = 1)) +
    # pretty_breaks rather than the default: with the cross-validation now running to the last
    # round whose outcome window closes by the analysis date, the tile has ~16 columns and the
    # default breaks print one label per column.
    scale_x_date(date_labels = "%d %b", breaks = scales::pretty_breaks(8)) +
    labs(title = sprintf("LFO: predicted invasion risk vs realised invasions (%s, h=%dw)", method, horizon),
         subtitle = "Rows = zones that were invaded; x = the FORECAST ORIGIN of the round in which the invasion occurred",
         x = "Forecast origin (as-of date)", y = NULL) +
    theme_inv(11)
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, sprintf("%s_h%d.pdf", file, horizon)),
                  w = 9, h = max(4, 0.35 * length(invaded_zones) + 2))
  invisible(p)
}

# ---------------------------------------------------------------------------
# 5b. Spatiotemporal skill over time (Q5) — per-fold discrimination vs date
# ---------------------------------------------------------------------------

#' @param st_skill output of spatiotemporal_skill(); @param methods optional
#'   subset of method names to draw (keeps the panel legible).
plot_skill_over_time <- function(st_skill, metric = "auc_pr_skill",
                                 methods = NULL, save = TRUE) {
  if (is.null(st_skill) || nrow(st_skill) == 0) return(invisible(NULL))
  d <- st_skill
  if (!is.null(methods)) d <- d %>% dplyr::filter(method %in% methods)
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
  # keep <= 6 methods for legibility (Okabe-Ito), rest dropped
  keep <- d %>% dplyr::group_by(method) %>%
    dplyr::summarise(m = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(m)) %>% dplyr::pull(method) %>% head(6)
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
    scale_colour_manual(values = setNames(OKABE_ITO[seq_along(keep)], keep), name = NULL) +
    scale_x_date(date_labels = "%d %b", breaks = scales::pretty_breaks(6)) +
    labs(title = "Forecast skill over time",
         subtitle = paste0("Per-fold skill against the LAST DAY OF THE WEEK BEING FORECAST, so both horizons share a calendar. ",
                           if (nrow(empty))
                             sprintf("Ticks on the x-axis mark %d round(s) in which no invasion was recorded, where discrimination is undefined. ",
                                     nrow(empty)) else "",
                           "Prevalence-normalised, so rounds with different event counts are comparable."),
         x = "Last day of the forecast window", y = ylab) +
    theme_inv(11)
  if (save) .save(p, file.path(OUT_DIAGNOSTICS, sprintf("skill_over_time_%s.pdf", metric)),
                  w = 11, h = 5)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 5c. Space-time risk heatmap (Q5) — top-risk zones × cutoff, invasions marked
# ---------------------------------------------------------------------------

#' Highest-risk at-risk zones (by peak predicted risk) across every fold, with the
#' realised first-case folds marked — surfaces both good calls and chronic false
#' alarms in one legible tile plot.
plot_spacetime_risk <- function(lfo_results, method, horizon = 1L, top_n = 20L,
                                save = TRUE, file = "spacetime_risk") {
  d <- lfo_results %>% dplyr::filter(method == !!method, horizon == !!horizon,
                                     is.finite(p_invasion))
  if (nrow(d) == 0) return(invisible(NULL))
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
    scale_fill_viridis_c(option = "plasma", name = "Predicted P",
                         labels = percent_format(accuracy = 1)) +
    scale_x_date(date_labels = "%d %b") +
    labs(title = sprintf("Space-time invasion risk (%s, h=%dw)", method, horizon),
         subtitle = "Top zones by peak risk; x = FORECAST ORIGIN of the round at which the zone recorded its first case",
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
    labs(title = "epinowcast right-truncation correction",
         subtitle = "Grey = observed (right-censored); blue = nowcast posterior (90% interval); ×N = correction factor",
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
write_invasion_report <- function(eval_tbl, risk_scores, lfo_results = NULL,
                                  zone_week = NULL, training_cutoff = NULL,
                                  best_method = NULL, primary_method = NULL,
                                  n_models = NA, out_path = NULL,
                                  best_bayes_method = NULL, bayes_weights = NULL) {
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
  fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
  L <- c()
  L <- c(L, "# BDBV 2026 — Spatiotemporal Invasion Forecast",
         "",
         sprintf("**Analysis date:** %s  ", outbreak_summary$analysis_date),
         sprintf("**Training cutoff:** %s  ", outbreak_summary$cutoff),
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
                            "prediction %s to %s (next 1-2 weeks).  "),
                     format(as.Date(.ts), "%d %b %Y"), format(as.Date(training_cutoff) + 6, "%d %b %Y"),
                     format(as.Date(training_cutoff) + 7, "%d %b %Y"),
                     format(as.Date(training_cutoff) + max(.hz) * 7, "%d %b %Y"))
           else ""
         },
         "",
         "## Task",
         "",
         "Predict, for each health zone with **no confirmed cases up to the forecast date**,",
         "the probability it records its first case(s) in the next 1-2 weeks (spatial",
         "*invasion*). Already-affected zones carry no invasion probability. Forecasts come",
         "from mobility-informed renewal models; the mobility force of infection is",
         "Lambda_i = sum_j W[j,i] sum_k g(k) Y_nc[j,t-k], the import coefficient beta is",
         "learned from realised invasions (cloglog), and P(first case) = 1 - exp(-beta*Lambda).",
         "",
         "The model is fitted in a **Bayesian** framework (brms/cmdstanr, cloglog with",
         "offset(log Lambda)). Every map and decision product is produced for the",
         "skill-selected featured model (figures prefixed `bayes_`; the loo-stacked ensemble",
         "prefixed `bayes_ensemble_`), and every figure's title names the exact model used.",
         "It is scored against three structural baselines on the identical folds: source-mass x",
         "distance-decay import pressure (B4), a distance-weighted spatial null (B1) and a",
         "nearest-affected/adjacency null (B7). The frequentist penalised-MLE renewal arm that",
         "used to run alongside it has been removed.",
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
           sprintf(paste0("Featured model: **%s** (Bayesian), selected ",
                          "by the leave-future-out CV composite (within-horizon ranks of AUC-PR ",
                          "skill + mean rank-of-truth + log score, pooled across horizons). The loo ",
                          "predictive-stacking weights below define the stacked ENSEMBLE only; they do ",
                          "not pick the featured single model.%s"),
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
    # models (fewer folds, e.g. a Bayesian fit that could not be estimated on an event-sparse
    # fold) are listed separately so a thin-subset score cannot masquerade as a full-CV win.
    has_partial <- "partial_cv" %in% names(e1)
    e1c <- (if (has_partial) e1 %>% dplyr::filter(!(partial_cv %in% TRUE)) else e1) %>%
      dplyr::arrange(dplyr::desc(auc_pr_skill))
    have_ci <- all(c("auc_pr_lo", "auc_pr_hi") %in% names(e1c))
    L <- c(L, "## Model comparison (leave-future-out, 1-week-ahead invasion)",
           "",
           sprintf("Pooled over %d at-risk zone-weeks with %d invasion events on the common support (base rate %.2f%%).",
                   e1c$n_atrisk[1] %||% NA, e1c$n_invasions[1] %||% NA, 100 * (e1c$base_rate[1] %||% NA)),
           "Ranked by AUC-PR skill (AUC-PR / base rate; 1 = no skill), with a 90% zone-cluster bootstrap",
           "interval. The top models overlap within CV noise — read the ranking as a leading *pack*, not a",
           "strict order. Count-WIS is deliberately NOT used to rank (it rewards predicting zero over the",
           "500+ never-affected zones).",
           "",
           if (have_ci) "| Model | AUC-PR skill [90% CI] | AUC-PR | Rank of truth | Log-score | Calib. |"
           else         "| Model | AUC-PR skill | AUC-PR | Rank of truth | Log-score | Calib. |",
           "|---|---|---|---|---|---|")
    for (i in seq_len(nrow(e1c))) {
      skill_cell <- if (have_ci && is.finite(e1c$auc_pr_lo[i]) && is.finite(e1c$base_rate[i]))
        sprintf("%.1f [%.1f-%.1f]", e1c$auc_pr_skill[i],
                pmax(e1c$auc_pr_lo[i], 0) / e1c$base_rate[i], e1c$auc_pr_hi[i] / e1c$base_rate[i])
        else sprintf("%.1f", e1c$auc_pr_skill[i])
      L <- c(L, sprintf("| %s | %s | %.3f | %.1f | %.3f | %.1fx |",
                        e1c$method[i], skill_cell, e1c$auc_pr[i],
                        e1c$mean_rank_of_truth[i], e1c$log_score[i], e1c$calibration_in_large[i]))
    }
    if (has_partial && any(e1$partial_cv %in% TRUE)) {
      pc <- e1 %>% dplyr::filter(partial_cv %in% TRUE) %>% dplyr::arrange(dplyr::desc(auc_pr_skill))
      L <- c(L, "",
             # "scored on <100% of folds" — the literal 100 sat where a threshold belongs, so
             # this rendered the vacuous "<100%". A partial-CV model is simply one that did
             # not cover every fold; state that, and let the per-model coverage carry the
             # detail rather than inventing a bound.
             sprintf("_Partial-CV (not co-ranked; these models did not cover every fold, so their skill denominators are not comparable): %s._",
                     paste(sprintf("%s (%.0f%% of folds)", pc$method,
                                   100 * (pc$coverage %||% NA)), collapse = "; ")))
    }
    L <- c(L, "", if (.have_primary)
             sprintf("**Primary operational model:** %s (best mobility-informed renewal).", primary_method)
           else sprintf("**Featured operational model:** %s (best cross-validated Bayesian model by the CV composite).",
                        primary_method),
           "",
           "> **Calibration caveat.** The models *discriminate* which zones are at risk very well",
           # Derived, not asserted: this said "AUC-ROC ~0.99" where the best model on this run
           # reaches 0.980 and the featured one 0.978. A published literal that rounds the
           # headline discrimination number UP is exactly the kind of claim a reviewer checks.
           sprintf("> (best AUC-ROC %s; realised invasions fall in the highest-probability bins), but their",
                   { .ar <- if (!is.null(eval_tbl) && "auc_roc" %in% names(eval_tbl))
                              suppressWarnings(max(eval_tbl$auc_roc[eval_tbl$horizon == 1], na.rm = TRUE))
                            else NA_real_
                     if (is.finite(.ar)) sprintf("%.3f", .ar) else "not available in this run" }),
           "> *absolute* probabilities over-predict out-of-sample (see the reliability figure) — the",
           "> invasion hazard per unit mobility force falls as susceptible neighbours of the front are",
           "> themselves invaded. Treat the absolute P(case) as an **upper bound** and prioritise on the",
           "> **ranking and relative-risk** scores, which are robust.", "")
  }

  if (!is.null(risk_scores) && "p_case_invasion" %in% names(risk_scores)) {
    rs1 <- risk_scores %>% dplyr::filter(horizon == 1, !was_active_before, !is.na(p_case_invasion))
    tbl <- function(dd, ttl) {
      dd <- dd %>% dplyr::arrange(dplyr::desc(p_case_invasion)) %>% head(10)
      c(sprintf("### %s", ttl), "",
        "| Zone | Province | P(case) | Rel. risk |",
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
  L <- c(L, "## Risk scores provided", "",
         "- **P(case)** — probability of >=1 CONFIRMED case. This is the only probability the",
         "  suite reports: it is the quantity the models are fitted to and evaluated against.",
         "- **Relative risk (Ituri / nationwide)** — multiple of the mean at-risk zone.",
         "- Ranked tables: `risk_table_ituri.csv`, `risk_table_national.csv`.", "")

  # ---- Findings / strengths / weaknesses / improvements (data-driven headline + audit) ----
  if (!is.null(eval_tbl)) {
    ez <- eval_tbl %>% dplyr::filter(horizon == 1)
    if ("partial_cv" %in% names(ez)) ez <- ez %>% dplyr::filter(!(partial_cv %in% TRUE))
    top <- ez %>% dplyr::arrange(dplyr::desc(auc_pr_skill)) %>% head(1)
    rcol <- if ("recall_at_10" %in% names(ez)) "recall_at_10" else "recall_at_5"
    kk   <- if (rcol == "recall_at_10") 10L else 5L
    # The targeting sentence below describes THE featured model, so read that model's
    # recall — not the maximum over the whole grid, which is a best-of-~30 statistic and
    # overstates what the published watch-list actually catches. Falls back to the field
    # maximum only when the featured model is absent, and says so.
    .feat_for_rec <- intersect(c(best_bayes_method, best_method, primary_method), ez$method)
    .rec_is_max <- !length(.feat_for_rec)
    best_rec <- if (.rec_is_max) suppressWarnings(max(ez[[rcol]], na.rm = TRUE))
                else ez[[rcol]][match(.feat_for_rec[1], ez$method)]
    # random watch-list baseline: K per-fold at-risk zones = K * n_folds / pooled n_atrisk.
    .nf_r <- if ("n_folds" %in% names(ez)) stats::median(ez$n_folds, na.rm = TRUE) else 1
    rand_rec <- if ("n_atrisk" %in% names(ez)) kk * .nf_r / stats::median(ez$n_atrisk, na.rm = TRUE) else NA_real_
    best_roc <- suppressWarnings(max(ez$auc_roc, na.rm = TRUE))
    L <- c(L, "## Findings", "",
      # "paradigm-robust" was a claim about agreement between a frequentist and a Bayesian
      # featured model. There is no frequentist arm, so with .have_primary FALSE the sentence
      # would have compared the Bayesian model with ITSELF and called the result robustness.
      if (.have_primary)
        sprintf("- **Discrimination is strong and paradigm-robust.** Best AUC-ROC %.2f; best AUC-PR skill %.0fx the base rate (chance = 1). The frequentist (%s) and Bayesian (%s) featured models rank zones near-identically.",
                best_roc, top$auc_pr_skill[1] %||% NA, primary_method %||% "n/a", best_bayes_method %||% "n/a")
      else
        sprintf("- **Discrimination is strong.** Best AUC-ROC %.2f; best AUC-PR skill %.0fx the base rate (chance = 1), for the featured Bayesian model `%s`. The structural baselines (source-mass x distance decay, distance-weighted, nearest-affected) are scored on the identical folds.",
                best_roc, top$auc_pr_skill[1] %||% NA, best_bayes_method %||% "n/a"),
      sprintf("- **Operationally useful targeting.** A top-%d watch-list from %s catches %.0f%% of the next invasions vs %.0f%% for a random list of the same size — a %.0fx lift.",
              kk, if (.rec_is_max) "the best-performing model in the grid"
                  else sprintf("the featured model (`%s`)", .feat_for_rec[1]),
              100 * best_rec, 100 * (rand_rec %||% NA), best_rec / (rand_rec %||% NA_real_)),
      "- **Absolute probabilities are upper bounds.** Discrimination >> calibration: the ranking and relative risk are trustworthy; the absolute P(case) over-predicts out-of-sample (see reliability figure).",
      "- **Uncertainty is reported everywhere** — the Bayesian posterior 90% CrI on every map (`*_uncertainty_map_*`), and a cluster-bootstrap CI on every skill estimate.",
      "")
    L <- c(L, "## Strengths", "",
      "- **Leakage-free evaluation.** Leave-future-out CV with per-fold nowcasting as of each origin (cut+6, the last day of the last training week — the deployed convention); ground truth is the raw observed first-case week; models are compared on a single common support with one shared base rate.",
      "- **Robustness across structural assumptions** — the featured model is one of a grid spanning mobility kernels, generation times and covariate sets, all scored on the same folds, with a loo-stacked posterior ensemble as the combination anchor.",
      "- **Mechanistic mobility structure beats naive baselines and simple comparators** (see `model_vs_baseline`), so the skill is from the import-force model, not the base rate.",
      "- **Real-time ready.** Dates derive from the data snapshot; the onset->sample delay can be re-estimated from the data being processed (lab or DHIS2).",
      "")
    L <- c(L, "## Weaknesses & limitations", "",
      "- **Absolute calibration.** Out-of-sample over-prediction of P(case); use ranking/relative risk for decisions, not the raw probability.",
      sprintf("- **Few events.** Only ~%s invasion events inform the h=1 fit; skill CIs are wide and the leaderboard is a pack, not a strict order.",
              top$n_invasions[1] %||% "a handful"),
      "- **No BDBV-specific generation time exists;** a Zaire-ebolavirus serial-interval proxy is used (short/med/long sensitivity brackets the published range; the Bayesian GT-preference figure shows which the data favour).",
      # The share is DERIVED. "~15%" was hard-coded and wrong: the figure code computing the same
      # quantity on the same line list gives 24.4% (the gap is largely the sitrep-appended rows,
      # which carry a sample date and no onset). The estimator description was also stale — the
      # draw is now from the fitted, truncation-corrected onset->sample delay, not the empirical pool.
      # It is NO LONGER the same distribution the nowcast and R(t) use, and saying so was wrong
      # from 2026-09-21: those correct for onset->APPEARANCE (a fitted truncation, mean ~9.8 d),
      # while imputation targets the SAMPLE date and so keeps the onset->sample fit (mean ~7.7 d).
      # Both are right for their own job; the sentence below no longer claims they are one.
      sprintf("- **Onset imputation** (%s of confirmed cases) uses a single stochastic draw from the fitted, truncation-corrected onset->sample delay (the delay to the sample date, which is what is being imputed; the nowcast and R(t) instead use a separately fitted onset-to-appearance truncation); a single draw does not fully propagate imputation uncertainty (multiple imputation is the extension), so imputation-dependent intervals are conditionally slightly narrow.",
              { .pi <- get0(".PCT_ONSET_IMPUTED", ifnotfound = NA_real_)
                if (is.finite(.pi)) sprintf("~%.0f%%", .pi) else "a material share" }),
      "- **Two-stage nowcast.** Epinowcast-corrected counts feed the models as a fixed input; the nowcast-sensitivity figure shows robustness to the nowcast choice but does not propagate its uncertainty (a joint model would).",
      "- **Diagnostic R(t)** is right-truncation-corrected but still uncertain at the most recent 1-2 weeks; it does not feed the forecasts.",
      "- **Per-fold LFO nowcast** is a fast delay-CDF approximation (delay-consistent with epinowcast); the full probabilistic epinowcast is applied only to the live current week (refitting it per fold is intractable).",
      "- **Duplicated health-zone names** across provinces (Bili, Lubunga) carry an arbitrary province label in free text (maps disambiguate via a name+province join).",
      "")
    L <- c(L, "## Where to improve", "",
      "- **Recalibrate absolute probabilities** (e.g. isotonic / Platt on held-out folds, or a decaying-susceptibility term) so P(case) is usable directly.",
      "- **Estimate a joint delay+transmission model** for the recent tail rather than a two-step nowcast->R(t).",
      "- **Tune the mobility/decay kernels by LFO** (currently fixed defaults for the minor kernels) and add a data-driven weight on IDP movement.",
      "- **Carry posterior draws** out of the Bayesian predictor so the stacked ensemble interval is an exact mixture rather than a normal approximation.",
      "- **Refit epinowcast per fold** if compute allows, to make the LFO nowcast fully probabilistic.",
      "")
  }
  writeLines(L, out_path)
  message("[report] invasion report -> ", out_path)
  invisible(out_path)
}

# ---------------------------------------------------------------------------
# Reliability diagram — predicted vs observed invasion frequency
# ---------------------------------------------------------------------------

#' Reliability of invasion probabilities on the pooled at-risk LFO rows. Few
#' bins (events are scarce). Points below the diagonal = over-confident. This is
#' the honest counterpart to the ranking metrics: the model discriminates well
#' but its absolute probabilities are upper bounds out-of-sample.
#' Rank-based AUC-ROC (Mann-Whitney): P(random invaded ranked above random not).
# NOTE: this function is ALSO defined in 16_invasion_eval.R. run_all.R sources 16 first and
# 17 last, so whichever body sits here WINS for the whole session — including for
# evaluate_invasion() in 16 and spatiotemporal_skill() in 19, which both call `.auc_roc`.
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
#' signature high-risk callout (the small share of zone-weeks scored high risk and
#' the much higher observed invasion rate among them).
plot_reliability <- function(lfo_results, horizon = 1L, methods = NULL, save = TRUE,
                             panel_tag = NULL, file = "reliability") {
  d <- lfo_results %>% dplyr::filter(horizon == !!horizon, is.finite(p_invasion))
  if ("was_active_before" %in% names(d))
    d <- d %>% dplyr::filter(!(was_active_before %in% TRUE))
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
  brks <- c(-1e-9, 0.005, 0.02, 0.05, 0.15, 1)
  rel <- d %>% dplyr::filter(method %in% methods) %>%
    dplyr::mutate(bin = cut(p_invasion, brks)) %>%
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
  # High-risk callout (their "0.31% of district-weeks >25%/wk -> 77%/month" line):
  hi <- d %>% dplyr::filter(method == methods[1])
  hi_share <- mean(hi$p_invasion > 0.15, na.rm = TRUE)
  hi_obs   <- if (any(hi$p_invasion > 0.15)) mean(hi$is_new_invasion[hi$p_invasion > 0.15]) else NA_real_
  callout <- if (is.finite(hi_obs))
    sprintf("The %.1f%% of at-risk zone-weeks we score >15%%/week had a %.0f%% observed invasion rate.",
            100 * hi_share, 100 * hi_obs) else "Pooled over folds; few events, so absolute levels are indicative."
  # Shared x/y maximum so the y = x diagonal is a true 45-degree reference and
  # over-/under-confidence read correctly (predicted and observed on the SAME scale).
  .axmax <- max(c(rel$pm, rel$ym), 0.02, na.rm = TRUE) * 1.05
  p <- ggplot(rel, aes(pm, ym)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_line(aes(colour = method), linewidth = 0.8) +
    geom_point(aes(colour = method, size = n)) +
    annotate("label", x = 0, y = Inf, hjust = 0, vjust = 1.05, label = auc_lab,
             size = 2.9, label.size = 0, fill = "white", alpha = 0.7, colour = "grey20") +
    scale_colour_manual(values = setNames(
      if (length(methods) <= length(OKABE_ITO)) OKABE_ITO[seq_along(methods)]
      else grDevices::hcl.colors(length(methods), "Dark 3"), methods), name = NULL) +
    scale_size_area(max_size = 7, name = "zone-weeks") +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, .axmax)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, .axmax)) +
    coord_equal(xlim = c(0, .axmax), ylim = c(0, .axmax)) +
    labs(title = sprintf("%sCalibration: predicted vs observed invasion (h = %dw)",
                         if (!is.null(panel_tag)) paste0(panel_tag, "  ") else "", horizon),
         subtitle = "On the diagonal = calibrated probabilities; below = over-confident.",
         x = "Predicted invasion probability", y = "Observed invasion frequency",
         caption = callout) +
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
                                 save = FALSE, panel_tag = "(A)"), error = function(e) NULL)
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
#' Panel A forest: each fold's INDEPENDENTLY fitted delta with its 95% interval, against the
#' random-effects pooled value. Panel B: I-squared for every method, i.e. the share of the
#' between-fold spread that is NOT sampling noise.
#'
#' EVERY NUMBER IS READ from invasion_delta_stability{,_summary}.csv, which run_all.R writes
#' from invasion_delta_stability() (16b). Nothing here re-fits, re-pools or re-tests.
#'
#' @param per_fold,summary_tbl the two published tables (or paths to them).
#' @param method which method to draw in panel A; defaults to the first alphabetically.
plot_delta_stability <- function(per_fold, summary_tbl, horizon = 1L, method = NULL,
                                 save = TRUE, file = "delta_stability") {
  rd <- function(x) if (is.character(x)) readr::read_csv(x, show_col_types = FALSE) else x
  pf <- rd(per_fold); sm <- rd(summary_tbl)
  if (is.null(pf) || !nrow(pf) || is.null(sm) || !nrow(sm)) return(invisible(NULL))
  pf <- pf %>% dplyr::filter(horizon == !!horizon, is.finite(delta), is.finite(se_log))
  sm <- sm %>% dplyr::filter(horizon == !!horizon)
  if (!nrow(pf) || !nrow(sm)) return(invisible(NULL))
  if (is.null(method)) {
    .pref <- c(get0("best_bayes_method", ifnotfound = NA_character_),
               "Bayes-M14-fill-med", "Bayes-M14-fill-geo")
    .pref <- .pref[!is.na(.pref) & .pref %in% pf$method]
    method <- if (length(.pref)) .pref[1] else sort(unique(pf$method))[1]
  }
  a <- pf %>% dplyr::filter(method == !!method)
  srow <- sm %>% dplyr::filter(method == !!method)
  if (!nrow(a) || !nrow(srow)) return(invisible(NULL))
  a$lab <- if (all(!is.na(a$cutoff))) format(lfo_origin(a$cutoff), "%d %b") else as.character(a$fold_id)
  a$lab <- factor(a$lab, levels = rev(a$lab[order(a$fold_id)]))

  pA <- ggplot(a, aes(x = .data$delta, y = .data$lab)) +
    annotate("rect", xmin = srow$pooled_lo[1], xmax = srow$pooled_hi[1],
             ymin = -Inf, ymax = Inf, fill = "grey88") +
    geom_vline(xintercept = srow$pooled_delta[1], colour = "grey35", linewidth = 0.6) +
    geom_vline(xintercept = 1, linetype = "dotted", colour = "grey55") +
    geom_linerange(aes(xmin = .data$lo, xmax = .data$hi), colour = "grey30", linewidth = 0.6) +
    geom_point(aes(size = .data$n_events), colour = "grey10") +
    scale_size_continuous(range = c(1.2, 3.4), name = "events") +
    scale_x_log10() +
    labs(x = expression(paste("fold-specific ", delta, "  (log scale; 1 = already calibrated)")),
         y = NULL, subtitle = sprintf("%s  -  band = pooled delta %.2f [%.2f, %.2f]",
                                      method, srow$pooled_delta[1], srow$pooled_lo[1],
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

  cap <- sprintf(paste0("Each fold's delta is fitted on that fold ALONE, so the spread is not an ",
                        "expanding window converging. Cochran Q p = %.3f on %d df; tau^2 = %.4f; ",
                        "I^2 = %.0f%%; trend in log delta per fold = %+.3f (p = %.2f); a single ",
                        "fold pins delta only to x/%.2f. Heterogeneity at p<0.05 in %d of %d ",
                        "methods. The anchor-week nowcast multiplier is constant across folds ",
                        "by construction and so cannot be tested."),
                 srow$Q_p[1], srow$Q_df[1], srow$tau2[1], 100 * srow$I2[1],
                 srow$trend_slope[1], srow$trend_p[1], exp(srow$median_se_log[1]),
                 sum(sm$Q_p < 0.05, na.rm = TRUE), nrow(sm))
  p <- (pA | pB) + patchwork::plot_layout(widths = c(1, 1.15)) +
    patchwork::plot_annotation(
      title = sprintf("Is the calibration offset stable over time?  (h = %s)", horizon),
      caption = paste(strwrap(cap, width = 150), collapse = "\n"),
      theme = ggplot2::theme(
        plot.title   = ggplot2::element_text(size = 12, face = "bold"),
        plot.caption = ggplot2::element_text(size = 7.5, colour = "grey35", hjust = 0)))
  if (save) .save(p, file.path(OUT_FIGURES, sprintf("%s_h%d.pdf", file, horizon)), w = 12, h = 5.2)
  invisible(p)
}
