# =============================================================================
# make_manuscript_figure2_cascade.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
# A cascade analogue of manuscript_figures/Figure2, for the 3-MONTH invasion
# cascade's leave-future-out backtest (NOT the short-term 2-week forecast).
# Same house aesthetic as make_manuscript_figures.R Figure 2; TWO panels:
#
#   A  Prioritisation — share of true invasions caught when the top-K highest-risk
#      zones are watched each round, for the cascade vs structural baselines
#      (gravity model, Flowminder epicentre inflow, road travel time) vs a random watch-list.
#      90% origin-cluster bootstrap band drawn on the cascade only.
#   B  Forecast-vs-outcome — the top-12 predicted zones at each backtest origin,
#      coloured by the realised outcome: invaded WITHIN the round's +6-week window,
#      invaded LATER (but by the last data day used for fitting), or not invaded by then.
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
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  # Gate on the FULL destination path, not the bare stem: FIGURE_DROP entries are
  # "<directory>/<stem>" and the raw/ exclusion inspects path components, neither of
  # which can match a basename.
  if (is.function(.fk) && !.fk(file.path(dir, name))) return(invisible(p))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  message(sprintf("  saved %-32s  %.1f x %.1f in", name, w, h)); invisible(p)
}

# (.wilson() removed: the capture-curve interval is a cluster bootstrap over origins, not a
#  binomial interval on pooled rows — see det_curve() below.)
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

# Last calendar day of data used for fitting (run_info training_window_end). The
# "invaded later / not invaded" classes are judged against it, so the legend quotes it.
# No fallback: a wrong date on a published legend is worse than a failed build.
DATA_END_DATE <- {
  ri <- jsonlite::fromJSON(file.path(OUT_DIR, "key_outputs", "run_info.json"))
  d <- if (is.null(ri$training_window_end)) as.Date(NA) else as.Date(ri$training_window_end)
  if (length(d) != 1 || is.na(d)) stop("run_info.json lacks training_window_end (last data day used for fitting)")
  d
}
DATA_END_LBL <- format(DATA_END_DATE, "%d %b")

# zones invaded at SOME point by DATA_END_DATE (present cascade run; was_active_before
# marks zones already affected as of the analysis date) — the "invaded by DATA_END_DATE" set.
affected_ever <- reach %>% filter(as.logical(was_active_before) %in% TRUE) %>%
  distinct(health_zone) %>% pull(health_zone)

# honest per-origin labels: Round 1..N (these are the full set of backtest origins)
fold_ord <- sort(unique(detail$cutoff))
lab_lk   <- setNames(sprintf("Round %d: %s", seq_along(fold_ord), format(lfo_origin(fold_ord), "%d %b")),
                     as.character(fold_ord))

# ---- THREE static structural baselines, scored on the SAME backtest rows -----
# Identical set, identical construction and identical labels to manuscript Figure 2, so a
# reader comparing the two figures sees one baseline set. All three rank zones by how strongly
# they connect to the EPICENTRE, differing ONLY in the connectivity matrix; none reads
# incidence and none touches the renewal machinery.
#   (1) Gravity     — M4, the FITTED gravity kernel (both masses, all exponents estimated).
#   (2) Flowminder  — M_cohort, Flowminder cohort subscriber-day presence (during-outbreak).
#                     Observed mobility, no model.
#   (3) Travel time — OSRM road travel time, 1/(1 + minutes) to the nearest epicentre zone.
#
# THE SCORERS ARE THE PIPELINE'S. naive_epicentre_inflow_scores() and
# epicentre_travel_time_scores() (20_forecast_detail.R) are called directly rather than
# reimplemented here. The previous hand copy omitted the alias harmonisation those apply, so a
# non-canonical epicentre spelling silently changed this figure's baseline but not Figure 2's.
#
# WHAT CHANGED AND WHY. This figure previously showed only TWO comparators — an inflow score
# on MOBILITY_PRIMARY (M8-fill at the time) and "Proximity to epicentre". M8-fill's epicentre rows are
# only 61-72% Flowminder; the rest is gravity fill over the destinations Flowminder did not
# measure, which made the inflow and gravity baselines partly the same model. And "proximity"
# was inverse great-circle distance to the ORIGINAL epicentre — static, never updating, and
# not a travel-time quantity at all. There was no gravity baseline.
# 03 supplies load_osrm() + harmonise_names(); 20 supplies the two scorers. Both are
# function-definition modules (their only top-level statement re-sources 00_config.R), so
# sourcing them runs nothing and writes nothing.
if (!exists("load_osrm", mode = "function"))
  source(file.path(ST_DIR, "03_mobility_matrices.R"))
if (!exists("naive_epicentre_inflow_scores", mode = "function") ||
    !exists("epicentre_travel_time_scores", mode = "function"))
  source(file.path(ST_DIR, "20_forecast_detail.R"))

build_baselines <- function(zones) {
  MD <- file.path(OUT_DIR, "mobility")
  .k <- function(id) {
    f <- file.path(MD, sprintf("mobility_%s.rds", id))
    if (file.exists(f)) readRDS(f) else NULL
  }
  wp <- read_csv(file.path(DATA_DIR, "worldpop", "processed", "worldpop__pop_count__static.csv"),
                 show_col_types = FALSE)
  pv <- setNames(wp$pop_count, wp$nom)
  # load_osrm() is the pipeline's own loader (03_mobility_matrices.R); default kind is
  # "travel_time" in minutes — the same matrix run_all.R passes as osrm_mat.
  om <- tryCatch(load_osrm("travel_time"), error = function(e) {
    warning("[fig2-cascade] OSRM travel-time matrix unavailable: ", conditionMessage(e),
            call. = FALSE); NULL })

  # M_cohort is the Flowminder COHORT kernel (subscriber-day presence, during-outbreak) —
  # the best mobility data available and far less censored than the short-trip annex M1
  # (305 destinations vs 142), which is kept only as a fallback. Both are POOLED over their
  # origin zones, so the epicentre rows are identical within each.
  M4 <- .k("M4"); MC <- .k("M_cohort"); M1 <- .k("M1")
  if (is.null(MC) && !is.null(M1)) {
    warning("[fig2-cascade] M_cohort not on disk; the Flowminder baseline falls back to the ",
            "pooled short-trip kernel M1. Re-run the mobility build.", call. = FALSE)
    MC <- M1
  }
  zall <- if (!is.null(M4)) colnames(M4) else if (!is.null(MC)) colnames(MC) else names(pv)

  .score <- function(v, what) {
    if (is.null(v) || !any(is.finite(v) & v > 0)) {
      warning(sprintf("[fig2-cascade] the %s baseline is unavailable or all-zero; it is ",
                      what), "omitted from Panel A.", call. = FALSE)
      return(NULL)
    }
    v
  }
  grav <- .score(if (is.null(M4)) NULL else
    naive_epicentre_inflow_scores(M4, EPICENTRE_ZONES, pv, zall), "gravity (M4)")
  flow <- .score(if (is.null(MC)) NULL else
    naive_epicentre_inflow_scores(MC, EPICENTRE_ZONES, pv, zall), "Flowminder cohort inflow")
  ttim <- .score(if (is.null(om)) NULL else
    epicentre_travel_time_scores(om, EPICENTRE_ZONES, zall), "travel-time (OSRM)")

  pick <- function(v) if (is.null(v)) rep(0, length(zones)) else
    ifelse(zones %in% names(v), v[zones], 0)
  data.frame(health_zone = zones,
             grav = pick(grav), flow = pick(flow), ttim = pick(ttim))
}
bl <- build_baselines(unique(detail$health_zone))
dd <- detail %>% left_join(bl, by = "health_zone") %>%
  mutate(grav = coalesce(grav, 0), flow = coalesce(flow, 0), ttim = coalesce(ttim, 0))

# =============================================================================
# Panel A — prioritisation (capture curve): cascade vs baselines vs random
# =============================================================================
# recall pooled over origins, ranked WITHIN origin (mirrors make_manuscript_figures.R::.det_curve).
# SAME TWO ESTIMATORS AS THE REST OF THE SUITE. This used to compute recall POOLED over
# origins (sum(hits) / sum(invasions)) with a WILSON interval — both of which
# compute_detection_curve() (20_forecast_detail.R) and manuscript Figure 2 identify as
# defects and have replaced:
#   * the published recall_at_K is the per-ORIGIN share AVERAGED over origins, so a pooled
#     figure carries a different estimand under the same name ("share of invasions caught");
#   * a Wilson interval is binomial and treats the hundreds of zone rows within one origin as
#     independent, which is badly anticonservative on clustered data. The interval here is a
#     cluster bootstrap over ORIGINS, the independent unit, matching the manuscript figure.
det_curve <- function(d, scorecol, ks = 1:25, n_boot = 400L,
                      seed = get0("RANDOM_SEED", ifnotfound = 20260704L)) {
  pf <- d %>% group_by(cutoff) %>%
    mutate(rk = rank(-.data[[scorecol]], ties.method = "max")) %>% ungroup()
  tot <- sum(pf$y, na.rm = TRUE); if (tot == 0) return(NULL)
  natr <- pf %>% count(cutoff) %>% pull(n) %>% mean()
  idx_by_origin <- split(seq_len(nrow(pf)), pf$cutoff)
  rows <- lapply(ks, function(k) {
    inK <- pf$rk <= k
    # per-origin recall, then the mean over origins
    per <- vapply(idx_by_origin, function(ix) {
      yv <- pf$y[ix]; np <- sum(yv, na.rm = TRUE)
      if (np == 0) return(NA_real_)
      sum(yv[inK[ix]], na.rm = TRUE) / np
    }, numeric(1))
    per <- per[is.finite(per)]
    ci <- if (length(per) >= 2L) {
      .old <- if (exists(".Random.seed", envir = globalenv()))
        get(".Random.seed", envir = globalenv()) else NULL
      set.seed(as.integer((as.numeric(seed) + 7919 * k) %% 2147483647))
      bs <- vapply(seq_len(n_boot),
                   function(i) mean(per[sample.int(length(per), length(per), TRUE)]), numeric(1))
      if (!is.null(.old)) assign(".Random.seed", .old, envir = globalenv())
      unname(stats::quantile(bs, c(0.05, 0.95), names = FALSE, na.rm = TRUE))
    } else c(NA_real_, NA_real_)
    data.frame(k = k, recall = mean(per), recall_lo = ci[1], recall_hi = ci[2],
               recall_pooled = sum(pf$y[inK], na.rm = TRUE) / tot,
               recall_random = pmin(k / natr, 1))
  })
  do.call(rbind, rows)
}

# LABELS AND COLOURS MATCH manuscript Figure 2 exactly (BASELINE_LBL / BASELINE_COL in
# make_manuscript_figures.R), so the same comparator is the same colour and the same words in
# both figures.
CASC_LBL <- "3-month cascade"
GRAV_LBL <- "Gravity model (fitted flows)"
FLOW_LBL <- "Flowminder cohort inflow from epicentre"
TTIM_LBL <- "Road travel time from epicentre"
MODEL_COL <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7")
names(MODEL_COL) <- c(CASC_LBL, GRAV_LBL, FLOW_LBL, TTIM_LBL)

.cc <- function(col, lbl) {
  if (!any(is.finite(dd[[col]]) & dd[[col]] > 0)) {
    message("[fig2-cascade] ", lbl, " is all-zero on the backtest rows; omitted from Panel A.")
    return(NULL)
  }
  out <- det_curve(dd, col)
  if (is.null(out)) return(NULL)
  out$method <- lbl
  out
}
cc_casc <- det_curve(dd, "p_invasion"); cc_casc$method <- CASC_LBL
curves <- bind_rows(cc_casc, .cc("grav", GRAV_LBL), .cc("flow", FLOW_LBL), .cc("ttim", TTIM_LBL)) %>%
  mutate(method = factor(method, levels = names(MODEL_COL)))
feat_c <- curves %>% filter(method == CASC_LBL)             # bootstrap band on the cascade only
rnd    <- cc_casc[, c("k", "recall_random")]

# Pooled summary (all origins) for the in-panel annotation: two discrimination measures, then
# two CALIBRATION measures.
#
# WHY CALIBRATION BELONGS HERE. AUC-PR skill and AUC-ROC are invariant to any monotone
# transform of the probabilities, so a model can top both while its probabilities are wrong by
# a factor of two — and these probabilities are used to decide how many zones to prepare, not
# only which. The Brier score already on this panel is dominated by the base rate in a rare-
# event setting and is close to unreadable in absolute terms, so it does not fill the gap.
#   * CALIBRATION-IN-THE-LARGE, predicted / observed: the single number an operational reader
#     needs. 1.0 = the expected invasion count matched the realised one; 1.4 = 40% over.
#   * ECE, the mean absolute gap between predicted and observed probability across
#     equal-count bins: whether the model is calibrated ACROSS THE RANGE, not just on average
#     (a model that over-predicts the top zones and under-predicts the rest can still have a
#     calibration-in-the-large of exactly 1).
.ece <- function(p, y, nb = 10L) {
  ok <- is.finite(p) & is.finite(y); p <- p[ok]; y <- as.integer(y[ok])
  if (length(p) < nb) return(NA_real_)
  b <- cut(rank(p, ties.method = "first"), breaks = nb, labels = FALSE)
  # WEIGHTED by bin size, which is what makes it an EXPECTED calibration error. Equal-count
  # bins make the weights equal anyway, but only when n divides exactly; weighting keeps the
  # definition right when it does not.
  w <- as.numeric(table(b)) / length(p)
  sum(w * abs(tapply(p, b, mean) - tapply(y, b, mean)), na.rm = TRUE)
}
pooled_ap  <- .auc_pr(dd$p_invasion, dd$y)
pooled_br  <- mean(dd$y, na.rm = TRUE)
pooled_cil <- mean(dd$p_invasion, na.rm = TRUE) / max(pooled_br, 1e-12)
disc_lab <- sprintf(paste0("Leave-future-out (+%d wk, %d origins)\nAUC-PR skill  %.2fx\n",
                           "AUC-ROC  %.2f\nBrier score  %.3f\n",
                           "Calibration (pred/obs)  %.2f\nCalibration error (ECE)  %.3f"),
                    K_HORIZON, N_ORIGINS, pooled_ap / max(pooled_br, 1e-9),
                    .auc_roc(dd$p_invasion, dd$y), mean((dd$p_invasion - dd$y)^2, na.rm = TRUE),
                    pooled_cil, .ece(dd$p_invasion, dd$y))

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
OUT_COL <- setNames(c("#B33005", "#F6BB6B", "grey80"),
                    c("Invaded within round",
                      sprintf("Invaded later (by %s)", DATA_END_LBL),
                      sprintf("Not invaded (by %s)", DATA_END_LBL)))
LV <- names(OUT_COL)
topk <- detail %>%
  group_by(cutoff) %>% slice_max(p_invasion, n = 12, with_ties = FALSE) %>% ungroup() %>%
  mutate(outcome = factor(case_when(
                     y == 1L                         ~ LV[1],   # invaded within this round's +K-week window
                     health_zone %in% affected_ever  ~ LV[2],   # invaded later, but by DATA_END_DATE
                     TRUE                            ~ LV[3]),  # not invaded by DATA_END_DATE
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
                      axis.text.y = element_text(size = 9.6, colour = INK), axis.text.x = element_text(size = 8.8), legend.position = "top",
                      # extra right margin: the last facet's unclipped strip title overhangs the panel
                      plot.margin = margin(6, 40, 6, 6))

# =============================================================================
# Assemble + save (no title, no caption — tags only)
# =============================================================================
save_dual(pA, "F2A_cascade_prioritisation",  5.2, 4.8)
save_dual(pB, "F2B_cascade_forecast_vs_outcome", 2.9 * length(fold_ord) + 1.0, 4.8)

fig <- (pA | pB) + plot_layout(widths = c(0.85, 1.5)) + plot_annotation(tag_levels = "A")
save_dual(fig, "Figure2_cascade", 13.8, 5.2, dir = MANU_DIR)

# THE FIGURE'S OWN DATA, under the name that claims to be it. 39_cascade_eval_figure.R used to
# write "Figure2_cascade_data.csv" holding the per-origin backtest summary — a frame with none
# of this figure's columns in it — while this, the published figure, shipped no data file at all.
# Panel A is the prioritisation curves (with the bootstrap band on the cascade arm); panel B is the
# per-origin top-12 zones and their realised outcomes. Both are written, long, with a `panel` key.
local({
  key_dir <- file.path(OUT_DIR, "key_outputs")
  dir.create(key_dir, recursive = TRUE, showWarnings = FALSE)
  .a <- dplyr::mutate(curves, panel = "A_prioritisation")
  .b <- dplyr::mutate(dplyr::select(topk, dplyr::any_of(c("origin", "cutoff", "health_zone",
                                                          "p_invasion", "outcome", "rank"))),
                      panel = "B_topk_outcomes")
  readr::write_csv(dplyr::bind_rows(.a, .b),
                   file.path(key_dir, "Figure2_cascade_data.csv"))
  message("[fig2-cascade] wrote Figure2_cascade_data.csv (", nrow(.a), " curve rows + ",
          nrow(.b), " top-K rows)")
})
message(sprintf("[done] Figure2_cascade (%d origins, +%d wk) -> %s", N_ORIGINS, K_HORIZON, MANU_DIR))
