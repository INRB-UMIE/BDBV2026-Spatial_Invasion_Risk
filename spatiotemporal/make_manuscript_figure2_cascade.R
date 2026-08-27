# =============================================================================
# make_manuscript_figure2_cascade.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
# A cascade analogue of manuscript_figures/Figure2, for the 3-MONTH invasion
# cascade's leave-future-out backtest (NOT the short-term 2-week forecast).
# Same house aesthetic as make_manuscript_figures.R Figure 2; TWO panels:
#
#   A  Prioritisation — share of true invasions caught when the top-K highest-risk
#      zones are watched each round, for the cascade vs structural baselines
#      (epicentre mobility inflow, proximity to epicentre) vs a random watch-list.
#      95% Wilson band drawn on the cascade only.
#   B  Forecast-vs-outcome — the top-12 predicted zones at each backtest origin,
#      coloured by the realised outcome: invaded WITHIN the round's +6-week window,
#      invaded LATER (but by the 31 Jul analysis snapshot), or not invaded by 31 Jul.
#
# Titles/subtitles/caption are intentionally OMITTED (the caption lives in the
# manuscript). Reads only saved CSVs/RDS (no re-simulation).
#
# Run:  Rscript spatiotemporal/make_manuscript_figure2_cascade.R
# Out:  outputs/key_outputs/manuscript_figures/Figure2_cascade.{pdf,png}
#       outputs/key_outputs/manuscript_figures/panels/F2{A,B}_cascade_*.{pdf,png}
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(readr); library(scales); library(sf); library(here)
})

ST_DIR <- Sys.getenv("CASCADE_ST_DIR", unset = file.path(here::here(), "spatiotemporal"))
source(file.path(ST_DIR, "00_config.R"))   # OUT_DIR, DATA_DIR, SHAPEFILE_PATH, EPICENTRE_ZONES

# ---- house style (verbatim from make_manuscript_figures.R Figure 2) ---------
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
base_family <- "sans"
F2_BASE <- 11.5
key <- function(x) tolower(trimws(x))

theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    axis.title = element_text(size = base - 0.4, colour = MUTED),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4), angle = 90),
    axis.text = element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
    legend.position = "top", legend.justification = "left",
    legend.title = element_text(size = base - 1.2, colour = MUTED),
    legend.text = element_text(size = base - 1.4, colour = INK),
    legend.key.height = unit(9, "pt"), legend.key.width = unit(15, "pt"),
    strip.text = element_text(size = base - 0.6, colour = INK, face = "bold",
                              margin = margin(3, 3, 3, 3)),
    plot.tag = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin = margin(6, 8, 6, 6))
}
reorder_within <- function(x, by, within)
  factor(paste(x, within, sep = "___"), levels = unique(paste(x, within, sep = "___"))[order(within, by)])
tidytext_scale_y <- function() scale_y_discrete(labels = function(z) sub("___.*$", "", z))

MANU_DIR  <- file.path(OUT_DIR, "key_outputs", "manuscript_figures")
PANEL_DIR <- file.path(MANU_DIR, "panels")
save_dual <- function(p, name, w, h, dir = PANEL_DIR) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  message(sprintf("  saved %-32s  %.1f x %.1f in", name, w, h)); invisible(p)
}

# Wilson score interval for a binomial proportion x/n (clamped to [0,1]).
.wilson <- function(x, n, z = 1.96) {
  if (!is.finite(n) || n <= 0) return(c(NA_real_, NA_real_))
  p <- x / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hlf <- z / d * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(max(0, ctr - hlf), min(1, ctr + hlf))
}
# step-wise average precision (tie-aware) — matches 33_cascade_eval.R::.casc_auc_pr
.auc_pr <- function(p, y) {
  ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  P <- sum(y == 1); if (!P || !length(p)) return(NA_real_)
  o <- order(p, decreasing = TRUE); ps <- p[o]; y <- y[o]
  tp <- cumsum(y == 1); fp <- cumsum(y == 0)
  keep <- c(ps[-length(ps)] != ps[-1], TRUE)
  prec <- (tp / (tp + fp))[keep]; rec <- (tp / P)[keep]
  sum(prec * diff(c(0, rec)), na.rm = TRUE)
}
# AUC-ROC via the Mann–Whitney identity (tie-averaged ranks)
.auc_roc <- function(p, y) {
  ok <- is.finite(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0); if (!n1 || !n0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# =============================================================================
# Inputs (all on disk; no re-simulation)
# =============================================================================
CASC_DIAG <- file.path(OUT_DIR, "cascade", "diagnostics")
detail <- read_csv(file.path(CASC_DIAG, "cascade_backtest_detail.csv"), show_col_types = FALSE) %>%
  mutate(cutoff = as.Date(cutoff), y = as.integer(y))
stopifnot(nrow(detail) > 0, all(c("cutoff","health_zone","p_invasion","y") %in% names(detail)))
reach  <- read_csv(file.path(OUT_DIR, "cascade", "tables", "cascade_reach_scores_all_zones.csv"),
                   show_col_types = FALSE)
K_HORIZON <- as.integer(round(mean(detail$K)))    # backtest window (weeks); 6 in production
N_ORIGINS <- dplyr::n_distinct(detail$cutoff)

# zones invaded at SOME point by the 31 Jul snapshot (present cascade run; was_active_before
# marks zones already affected as of the analysis date) — the "invaded by 31 Jul" set.
affected_ever <- reach %>% filter(as.logical(was_active_before) %in% TRUE) %>%
  distinct(health_zone) %>% pull(health_zone)

# honest per-origin labels: Round 1..N (these are the full set of backtest origins)
fold_ord <- sort(unique(detail$cutoff))
lab_lk   <- setNames(sprintf("Round %d: %s", seq_along(fold_ord), format(fold_ord, "%d %b")),
                     as.character(fold_ord))

# ---- static structural baselines, scored on the SAME backtest rows ----------
# (i) Epicentre mobility inflow: population-weighted M8 outflow FROM the epicentre zones to
#     every zone (identical construction to make_manuscript_figures.R's NAIVE_SCORES).
# (ii) Proximity to epicentre: inverse great-circle distance to the nearest epicentre centroid.
# Both are FIXED per-zone reference scores (a leave-future-out cascade is compared against a
# static naive front); every zone gets a finite score (0 where uncovered) so ranking never drops rows.
build_baselines <- function(zones) {
  W  <- readRDS(file.path(OUT_DIR, "mobility", "mobility_M8.rds"))
  wp <- read_csv(file.path(DATA_DIR, "worldpop", "processed", "worldpop__pop_count__static.csv"),
                 show_col_types = FALSE)
  pv <- setNames(wp$pop_count, wp$nom); zall <- colnames(W)
  ez <- intersect(as.character(EPICENTRE_ZONES), rownames(W)); stopifnot(length(ez) > 0)
  p <- as.numeric(pv[match(ez, names(pv))]); good <- is.finite(p) & p > 0
  src <- if (!any(good)) rep(1, length(ez)) else { p[!good] <- stats::median(p[good]); p }
  score <- as.numeric(crossprod(src, W[ez, , drop = FALSE])); names(score) <- colnames(W)
  naive <- setNames(rep(0, length(zall)), zall)
  cm <- intersect(zall, names(score)); naive[cm] <- score[cm]
  naive[intersect(zall, ez)] <- 0                       # epicentre self-rows carry no risk

  shp <- st_read(SHAPEFILE_PATH, quiet = TRUE)
  ct  <- suppressWarnings(st_coordinates(st_point_on_surface(st_geometry(shp))))
  cen <- data.frame(k = key(shp$Nom), lon = ct[, 1], lat = ct[, 2]) %>% distinct(k, .keep_all = TRUE)
  hav <- function(lo1, la1, lo2, la2) { R <- 6371; dl <- (lo2 - lo1) * pi / 180; dp <- (la2 - la1) * pi / 180
    a <- sin(dp / 2)^2 + cos(la1 * pi / 180) * cos(la2 * pi / 180) * sin(dl / 2)^2; 2 * R * asin(pmin(1, sqrt(a))) }
  epi_ct <- cen[cen$k %in% key(EPICENTRE_ZONES), ]; stopifnot(nrow(epi_ct) > 0)
  d2e <- vapply(seq_len(nrow(cen)), function(i)
    min(hav(cen$lon[i], cen$lat[i], epi_ct$lon, epi_ct$lat)), numeric(1))
  prox <- setNames(1 / (1 + d2e), cen$k)                # higher = closer to the epicentre

  # score every backtest zone (0 where uncovered)
  data.frame(health_zone = zones,
             naive = ifelse(zones %in% names(naive), naive[zones], 0),
             prox  = ifelse(key(zones) %in% names(prox), prox[key(zones)], 0))
}
bl <- build_baselines(unique(detail$health_zone))
dd <- detail %>% left_join(bl, by = "health_zone") %>%
  mutate(naive = coalesce(naive, 0), prox = coalesce(prox, 0))

# =============================================================================
# Panel A — prioritisation (capture curve): cascade vs baselines vs random
# =============================================================================
# recall pooled over origins, ranked WITHIN origin (mirrors make_manuscript_figures.R::.det_curve).
det_curve <- function(d, scorecol, ks = 1:25) {
  pf <- d %>% group_by(cutoff) %>%
    mutate(rk = rank(-.data[[scorecol]], ties.method = "max")) %>% ungroup()
  tot <- sum(pf$y, na.rm = TRUE); if (tot == 0) return(NULL)
  natr <- pf %>% count(cutoff) %>% pull(n) %>% mean()
  rows <- lapply(ks, function(k) {
    caught <- sum(pf$y[pf$rk <= k], na.rm = TRUE)
    ci <- .wilson(caught, tot)
    data.frame(k = k, recall = caught / tot, recall_lo = ci[1], recall_hi = ci[2],
               recall_random = pmin(k / natr, 1))
  })
  do.call(rbind, rows)
}

CASC_LBL <- "3-month cascade"; NAIVE_LBL <- "Epicentre mobility inflow"; PROX_LBL <- "Proximity to epicentre"
MODEL_COL <- c("#0072B2", "#D55E00", "#009E73"); names(MODEL_COL) <- c(CASC_LBL, NAIVE_LBL, PROX_LBL)

cc_casc  <- det_curve(dd, "p_invasion"); cc_casc$method  <- CASC_LBL
cc_naive <- det_curve(dd, "naive");      cc_naive$method <- NAIVE_LBL
cc_prox  <- det_curve(dd, "prox");       cc_prox$method  <- PROX_LBL
curves <- bind_rows(cc_casc, cc_naive, cc_prox) %>%
  mutate(method = factor(method, levels = names(MODEL_COL)))
feat_c <- curves %>% filter(method == CASC_LBL)             # Wilson band on the cascade only
rnd    <- cc_casc[, c("k", "recall_random")]

# pooled discrimination summary (all origins), for the in-panel annotation
pooled_ap  <- .auc_pr(dd$p_invasion, dd$y)
pooled_br  <- mean(dd$y, na.rm = TRUE)
disc_lab <- sprintf("Leave-future-out (+%d wk, %d origins)\nAUC-PR skill  %.2fx\nAUC-ROC  %.2f\nBrier score  %.3f",
                    K_HORIZON, N_ORIGINS, pooled_ap / max(pooled_br, 1e-9),
                    .auc_roc(dd$p_invasion, dd$y), mean((dd$p_invasion - dd$y)^2, na.rm = TRUE))

pA <- ggplot(curves, aes(k, recall, colour = method)) +
  geom_line(data = rnd, aes(k, recall_random), linetype = "22", colour = FAINT, linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = feat_c, aes(ymin = recall_lo, ymax = recall_hi, fill = method), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(data = feat_c, size = 1) +
  annotate("text", x = 15, y = 0.06, label = "random watch-list", colour = MUTED, size = 3.5, angle = 5) +
  annotate("text", x = 0.5, y = 0.99, hjust = 0, vjust = 1, label = disc_lab, size = 3.4, colour = INK, lineheight = 0.95) +
  scale_colour_manual(values = MODEL_COL, name = NULL, breaks = names(MODEL_COL)) +
  scale_fill_manual(values = MODEL_COL, guide = "none") +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Zones actively monitored per round (K)", y = "Share of true invasions caught") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_pub(F2_BASE) + theme(legend.position = "top")

# =============================================================================
# Panel B — forecast-vs-outcome (three-way), top-12 zones per origin
# =============================================================================
OUT_COL <- c("Invaded within round"       = "#B33005",
             "Invaded later (by 31 Jul)"  = "#F6BB6B",
             "Not invaded (by 31 Jul)"    = "grey80")
LV <- names(OUT_COL)
topk <- detail %>%
  group_by(cutoff) %>% slice_max(p_invasion, n = 12, with_ties = FALSE) %>% ungroup() %>%
  mutate(outcome = factor(case_when(
                     y == 1L                         ~ LV[1],   # invaded within this round's +K-week window
                     health_zone %in% affected_ever  ~ LV[2],   # invaded later, but by 31 Jul
                     TRUE                            ~ LV[3]),  # not invaded by 31 Jul
                   levels = LV),
         fold_lab = factor(lab_lk[as.character(cutoff)], levels = lab_lk[as.character(fold_ord)]),
         zone_w   = reorder_within(health_zone, p_invasion, fold_lab))

pB <- ggplot(topk, aes(p_invasion, zone_w, fill = outcome)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  facet_wrap(~ fold_lab, scales = "free_y", nrow = 1) + tidytext_scale_y() +
  scale_fill_manual(values = OUT_COL, name = NULL, drop = FALSE) +
  scale_x_continuous(labels = percent_format(1), limits = c(0, NA), expand = expansion(mult = c(0, 0.06)), breaks = pretty_breaks(3)) +
  labs(x = "Predicted invasion probability, P(first case)", y = NULL) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_pub(F2_BASE) + theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(9, "pt"), strip.clip = "off",
                      axis.text.y = element_text(size = 9.6, colour = INK), axis.text.x = element_text(size = 8.8), legend.position = "top")

# =============================================================================
# Assemble + save (no title, no caption — tags only)
# =============================================================================
save_dual(pA, "F2A_cascade_prioritisation",  5.2, 4.8)
save_dual(pB, "F2B_cascade_forecast_vs_outcome", 2.9 * length(fold_ord) + 1.0, 4.8)

fig <- (pA | pB) + plot_layout(widths = c(0.85, 1.5)) + plot_annotation(tag_levels = "A")
save_dual(fig, "Figure2_cascade", 13.8, 5.2, dir = MANU_DIR)
message(sprintf("[done] Figure2_cascade (%d origins, +%d wk) -> %s", N_ORIGINS, K_HORIZON, MANU_DIR))
