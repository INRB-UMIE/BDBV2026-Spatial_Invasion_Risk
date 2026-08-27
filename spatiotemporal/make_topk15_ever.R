# ---------------------------------------------------------------------------
# Top-15 per-forecast-round outcome figure, ONE-WEEK horizon, THREE-WAY outcome:
#   invaded within the forecast week | invaded later (by 31 Jul) | not invaded.
# Mirrors panel E of Figure2_labelled (build_topk_folds) but at top-15, h=1, and
# splits the "invaded" category into within-window vs eventually. House aesthetic
# copied verbatim from make_publication_figures.R (theme_pub / save_dual / palette).
# ---------------------------------------------------------------------------
suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr); library(readr); library(scales)
}))

# --- house-style constants (verbatim from make_publication_figures.R) ---
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
base_family <- "sans"

theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    plot.title      = element_blank(),
    plot.subtitle   = element_blank(),
    plot.caption    = element_blank(),
    axis.title      = element_text(size = base - 0.4, colour = MUTED),
    axis.title.x    = element_text(margin = margin(t = 4)),
    axis.title.y    = element_text(margin = margin(r = 4), angle = 90),
    axis.text       = element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
    legend.position = "top", legend.justification = "left",
    legend.title    = element_text(size = base - 1.2, colour = MUTED),
    legend.text     = element_text(size = base - 1.4, colour = INK),
    legend.key.height = unit(9, "pt"), legend.key.width = unit(15, "pt"),
    strip.text      = element_text(size = base - 0.6, colour = INK, face = "bold",
                                   margin = margin(3, 3, 3, 3)),
    plot.tag        = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin     = margin(6, 8, 6, 6)
  )
}

PANEL_DIR <- "outputs/key_outputs/figures/panels"
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
save_dual <- function(p, name, w, h, dir = PANEL_DIR) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  message(sprintf("  saved %-34s  %.1f x %.1f in", name, w, h))
  invisible(p)
}
# tidytext-style within-facet ordering (verbatim)
reorder_within <- function(x, by, within) {
  factor(paste(x, within, sep = "___"),
         levels = unique(paste(x, within, sep = "___"))[order(within, by)])
}
tidytext_scale_y <- function() scale_y_discrete(labels = function(z) sub("___.*$", "", z))

# --- data ---
OUT        <- "outputs"
BEST_BAYES <- "Bayes-M17-med"          # featured model (current run)
rs  <- read_csv(file.path(OUT, "reports", "bayes_risk_scores_all_zones.csv"), show_col_types = FALSE)
lfo <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))

# zones invaded AT SOME POINT (recorded a first confirmed case by the 31 Jul snapshot)
affected_ever <- rs %>%
  filter(horizon == 1, as.logical(was_active_before) %in% TRUE) %>%
  pull(health_zone) %>% unique()
message(sprintf("[topk15-ever] featured=%s | zones invaded by 31 Jul: %d", BEST_BAYES, length(affected_ever)))

TOP_K <- 15L
d <- lfo %>%
  filter(method == BEST_BAYES, horizon == 1L, is.finite(p_invasion),
         !(as.logical(was_active_before) %in% TRUE)) %>%
  mutate(cutoff = as.Date(cutoff))
stopifnot(nrow(d) > 0)

fold_ord <- sort(unique(d$cutoff))
lab_lk   <- setNames(sprintf("Round %d — %s", seq_along(fold_ord), format(fold_ord, "%d %b")),
                     as.character(fold_ord))

OUT_WIN   <- "#B33005"   # invaded within the forecast week   (dark burnt orange)
OUT_LATER <- "#F6BB6B"   # invaded later, by 31 Jul            (light amber — clear dark/light split)
OUT_NEVER <- "grey80"    # not invaded by 31 Jul
LV <- c("Invaded within the forecast week", "Invaded later (by 31 Jul)", "Not invaded (by 31 Jul)")
OUT_COL <- setNames(c(OUT_WIN, OUT_LATER, OUT_NEVER), LV)

topk <- d %>%
  group_by(cutoff) %>%
  slice_max(p_invasion, n = TOP_K, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    outcome = factor(dplyr::case_when(
                is_new_invasion == 1                 ~ LV[1],
                health_zone %in% affected_ever       ~ LV[2],
                TRUE                                  ~ LV[3]),
              levels = LV),
    fold_lab = factor(lab_lk[as.character(cutoff)], levels = lab_lk[as.character(fold_ord)]),
    zone_w   = reorder_within(health_zone, p_invasion, fold_lab))

# --- accuracy self-check: per-fold counts must reproduce the validated table ---
chk <- topk %>% group_by(fold_lab) %>%
  summarise(n = n(), within = sum(outcome == LV[1]), ever = sum(outcome != LV[3]), .groups = "drop")
message("[topk15-ever] per-fold (n / within-window / ever-invaded):")
print(as.data.frame(chk))
message(sprintf("[topk15-ever] POOLED top-15: ever=%d/%d (%.0f%%), within=%d/%d (%.0f%%)",
                sum(topk$outcome != LV[3]), nrow(topk), 100*mean(topk$outcome != LV[3]),
                sum(topk$outcome == LV[1]), nrow(topk), 100*mean(topk$outcome == LV[1])))

p <- ggplot(topk, aes(p_invasion, zone_w, fill = outcome)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  facet_wrap(~ fold_lab, scales = "free_y", nrow = 3L) +
  tidytext_scale_y() +
  scale_fill_manual(values = OUT_COL, name = NULL, drop = FALSE) +
  scale_x_continuous(labels = percent_format(1), limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.06)), breaks = scales::pretty_breaks(3)) +
  labs(x = "Predicted invasion probability, P(first case)", y = NULL) +
  theme_pub(8.6) +
  theme(panel.grid.major.y = element_blank(),
        panel.spacing.x = unit(9, "pt"), panel.spacing.y = unit(8, "pt"),
        strip.clip = "off",
        axis.text.y = element_text(size = 6.6, colour = INK),
        legend.position = "top")

save_dual(p, "F_topk15_ever_h1", 9.6, 12.0)
message("[topk15-ever] done.")
