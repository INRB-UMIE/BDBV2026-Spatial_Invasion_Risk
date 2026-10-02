# =============================================================================
# 42_kindiv_sweep_figure.R — publication figure for the comprehensive
# individual-offspring overdispersion (k_indiv) sensitivity sweep.
#
# Reads the outputs of 41_kindiv_sweep.R and renders one 3-panel figure in the
# house publication aesthetic (Okabe-Ito, theme_pub, dual PDF + 600-dpi PNG):
#   A) Ranking stability vs k     — Spearman rho, Kendall tau, top-15 Jaccard
#   B) Magnitude vs k             — mean reach at 4/8/13 wk + mean establishment
#   C) Per-zone reach vs k        — top-12 at-risk zones, coloured by province
# The production baseline k = 0.30 is marked with a dashed reference line. (The
# Lloyd-Smith [0.2, 0.4] plausible-band box was removed by request.) Titles are
# blank (house style: captions live externally); panels are tagged A/B/C.
#
# Run:  Rscript spatiotemporal/42_kindiv_sweep_figure.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(scales); library(forcats); library(jsonlite)
})

# Prefer pipeline globals when sourced from run_cascade.R; else derive from cwd.
HERE <- normalizePath(".")
ST   <- if (basename(HERE) == "spatiotemporal") HERE else file.path(HERE, "spatiotemporal")
OUTC <- if (exists("OUT_CASCADE")) OUT_CASCADE else file.path(ST, "outputs", "cascade")
DIAG <- file.path(OUTC, "diagnostics")
FIGD <- file.path(OUTC, "figures")
KEYD <- if (exists("OUT_DIR")) file.path(OUT_DIR, "key_outputs") else file.path(ST, "outputs", "key_outputs")
for (d in c(FIGD, KEYD)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Load the retained-figure allow-list (FIGURE_KEEP / figure_is_kept, 00_config.R). Without it
# the gate in save_dual() below resolves figure_is_kept to NULL and silently no-ops, so
# FIGURE_DROP entries and the raw/ exclusion could never fire in this script even though it
# writes a retained deliverable (Figure_kindiv_sweep) alongside non-retained working panels.
# Non-fatal: the script must still run standalone from a checkout with no config on the path.
if (!is.function(get0("figure_is_kept")))
  try(suppressWarnings(suppressMessages(source(file.path(ST, "00_config.R")))), silent = TRUE)

# ---- design system (mirrors make_publication_figures.R) ---------------------
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
BAND_FILL <- "#0072B2"                       # Lloyd-Smith plausible band (soft, low alpha)
OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
PROV_COL <- c("Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73",
              "Sud-Kivu"="#CC79A7","Bas-Uele"="#E69F00","Tshopo"="#56B4E9","Other"="#7A7A7A")
base_family <- "sans"
theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    axis.title   = element_text(size = base - 0.4, colour = MUTED),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4), angle = 90),
    axis.text    = element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
    legend.position = "top", legend.justification = "left",
    legend.title = element_text(size = base - 1.2, colour = MUTED),
    legend.text  = element_text(size = base - 1.4, colour = INK),
    legend.key.height = unit(9, "pt"), legend.key.width = unit(15, "pt"),
    plot.tag = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin = margin(6, 8, 6, 6))
}
save_dual <- function(p, name, w, h, dirs = c(FIGD, KEYD)) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  # Gate on each destination path (these helpers write into several directories), so
  # FIGURE_DROP's "<dir>/<stem>" entries and the raw/ rule can actually match.
  if (is.function(.fk) && !any(vapply(dirs, function(.d) .fk(file.path(.d, name)),
                                      logical(1)))) return(invisible(p))
  for (d in dirs) {
    ggsave(file.path(d, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
    ggsave(file.path(d, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  }
  message(sprintf("  saved %-30s  %.1f x %.1f in", name, w, h))
}
# plausible-band [0.2, 0.4] box removed (user request); no-op layer (NULL is ignored
# when added to a ggplot), so the panel call sites below are left unchanged.
band_layer <- function(lo = 0.20, hi = 0.40) NULL
kbase_layer <- function(kb = 0.30)
  geom_vline(xintercept = kb, colour = FAINT, linetype = "22", linewidth = 0.35)
LOGX <- scale_x_log10(breaks = c(0.05,0.1,0.2,0.3,0.5,1,2,4),
                      labels = c("0.05","0.1","0.2","0.3","0.5","1","2","4"))

# ---- data -------------------------------------------------------------------
# This figure consumes 41_kindiv_sweep.R's outputs. That sweep is expensive and is NOT
# run on every pipeline pass, so a missing input is an ordinary state, not an error:
# exit cleanly with an explanatory message instead of an unhandled read_csv() failure.
# (Required now that run_all.R schedules this script — an unguarded stop() here would
# show up as a pipeline-stage failure whenever the sweep had simply not been run.)
.sweep_csv  <- file.path(DIAG, "cascade_kindiv_sweep.csv")
.byzone_csv <- file.path(DIAG, "cascade_kindiv_sweep_byzone.csv")
.missing <- c(.sweep_csv, .byzone_csv)[!file.exists(c(.sweep_csv, .byzone_csv))]
if (length(.missing)) {
  message("[42_kindiv_sweep_figure] Skipping — 41_kindiv_sweep.R has not been run ",
          "(missing: ", paste(basename(.missing), collapse = ", "), "). ",
          "Run `Rscript spatiotemporal/41_kindiv_sweep.R` first, or set ",
          "RUN_KINDIV_SWEEP=1 in run_all.R's downstream block.")
  # NEVER quit() unconditionally here. run_cascade.R SOURCES this file, and quit() from a
  # sourced script terminates the whole calling process — so a missing optional input used to
  # kill the entire cascade run WITH STATUS 0, silently skipping the evaluation figure,
  # Figure 4 and everything after it, while reporting success. sys.nframe() is 0 only when
  # R is executing this file directly (Rscript 42_....R); any source() puts frames on the
  # stack. Sourced: raise a condition the caller's tryCatch already handles as a skip.
  # Standalone and non-interactive: exit 0, which is the correct status for "nothing to do".
  if (sys.nframe() > 0L)
    stop("[42_kindiv_sweep_figure] required sweep inputs missing; figure skipped.", call. = FALSE)
  if (!interactive()) quit(save = "no", status = 0L)
}
sweep  <- read_csv(.sweep_csv,  show_col_types = FALSE)
byzone <- read_csv(.byzone_csv, show_col_types = FALSE)
meta   <- tryCatch(fromJSON(file.path(DIAG, "cascade_kindiv_sweep_meta.json")), error = function(e) NULL)
# k_base: prefer the value the SWEEP recorded (meta), then the live config constant, and only
# then a literal. The literal 0.30 mirrored CASCADE_K_INDIV, so the figure's baseline marker
# could silently disagree with the model if that constant ever changed — in a retained
# key_outputs deliverable, with nothing to reveal the mismatch.
KB <- if (!is.null(meta) && !is.null(meta$k_base) && is.finite(suppressWarnings(as.numeric(meta$k_base))))
  as.numeric(meta$k_base) else as.numeric(get0("CASCADE_K_INDIV", ifnotfound = 0.30))
if (!is.null(meta) && !is.null(meta$k_base) &&
    is.finite(suppressWarnings(as.numeric(meta$k_base))) &&
    exists("CASCADE_K_INDIV") &&
    abs(as.numeric(meta$k_base) - CASCADE_K_INDIV) > 1e-9)
  warning(sprintf(paste0("[42_kindiv_sweep_figure] the sweep was run at k_base = %s but the ",
                         "pipeline now uses CASCADE_K_INDIV = %s; the baseline marker on this ",
                         "figure is the SWEEP's value. Re-run 41_kindiv_sweep.R to realign."),
                  format(as.numeric(meta$k_base)), format(CASCADE_K_INDIV)), call. = FALSE)

.m1 <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) "?" else as.character(x)
# PROVENANCE CHECK. RUN_KINDIV_SWEEP defaults FALSE while RUN_CASCADE_FIGURES defaults TRUE,
# so this script runs every time and 41_kindiv_sweep.R usually does not. It previously skipped
# only when the CSVs were ABSENT — never when they were merely OLD — so Figure_kindiv_sweep (a
# retained deliverable) was re-stamped with a current mtime from whatever sweep happened to be
# on disk. That matters because CASCADE_KERNEL is derived from model_selection.json and can
# change between runs: the figure would then describe a sweep of a DIFFERENT model from the one
# the rest of the suite features. The sidecar records the sweep's own kernel/scenario/seed, so
# the mismatch is checkable. Warn loudly rather than refuse: the sweep is expensive and a stale
# one is still informative, but it must never be published as if it were this run's.
if (!is.null(meta)) {
  .now_kernel   <- get0("CASCADE_KERNEL",   ifnotfound = NA_character_)
  .now_scenario <- tryCatch(get0("CASCADE_SCENARIOS", ifnotfound = NULL)[[1]]$label,
                            error = function(e) NA_character_)
  .mis <- character(0)
  if (!is.na(.now_kernel) && !is.null(meta$kernel) &&
      !identical(as.character(meta$kernel), as.character(.now_kernel)))
    .mis <- c(.mis, sprintf("kernel (sweep %s vs this run %s)", meta$kernel, .now_kernel))
  if (!is.null(meta$seed) && !identical(as.integer(meta$seed),
                                        as.integer(get0("RANDOM_SEED", ifnotfound = meta$seed))))
    .mis <- c(.mis, sprintf("seed (sweep %s vs this run %s)", meta$seed,
                            get0("RANDOM_SEED", ifnotfound = NA)))
  if (length(.mis))
    warning(sprintf(paste0("[42_kindiv_sweep_figure] the k_indiv sweep on disk was produced under a ",
                           "DIFFERENT configuration - %s. Figure_kindiv_sweep therefore describes that ",
                           "sweep, not this run. Re-run 41_kindiv_sweep.R (RUN_KINDIV_SWEEP=1) to refresh it."),
                   paste(.mis, collapse = "; ")), call. = FALSE, immediate. = TRUE)
  else
    # Explicit guard rather than %||%: this script can run standalone, and base R only gained
    # %||% in 4.4, so relying on it here would make the provenance line an error on 4.3.
    message(sprintf("[42_kindiv_sweep_figure] sweep provenance OK (kernel %s, scenario %s, seed %s).",
                    .m1(meta$kernel), .m1(meta$scenario), .m1(meta$seed)))
}

# ---- Panel A: ranking stability --------------------------------------------
STAB_COL <- c("Spearman rank corr." = OKABE[1], "Kendall tau" = OKABE[2], "Top-15 Jaccard" = OKABE[3])
pa_df <- sweep %>%
  transmute(k_indiv,
            `Spearman rank corr.` = spearman, `Kendall tau` = kendall, `Top-15 Jaccard` = top15_jaccard) %>%
  pivot_longer(-k_indiv, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = names(STAB_COL)))
pA <- ggplot(pa_df, aes(k_indiv, value, colour = metric)) +
  band_layer() + kbase_layer(KB) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.25) +
  scale_colour_manual(values = STAB_COL, name = NULL) +
  LOGX + coord_cartesian(ylim = c(min(0.7, min(pa_df$value, na.rm = TRUE) - 0.02), 1.005)) +
  labs(x = expression("individual-offspring overdispersion  " * italic(k)),
       y = "agreement with k = 0.30 ranking") +
  theme_pub()

# ---- Panel B: magnitude (reach at horizons + establishment) ----------------
MAG_COL <- c("Reach, 4 wk" = "#9ecae1", "Reach, 8 wk" = "#4292c6", "Reach, 13 wk" = "#08519c",
             "Establishment, 13 wk" = OKABE[2])
pb_df <- sweep %>%
  transmute(k_indiv,
            `Reach, 4 wk` = mean_reach_h4, `Reach, 8 wk` = mean_reach_h8,
            `Reach, 13 wk` = mean_reach_h13, `Establishment, 13 wk` = mean_estab_h13) %>%
  pivot_longer(-k_indiv, names_to = "series", values_to = "value") %>%
  mutate(series = factor(series, levels = names(MAG_COL)))
pB <- ggplot(pb_df, aes(k_indiv, value, colour = series)) +
  band_layer() + kbase_layer(KB) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.15) +
  scale_colour_manual(values = MAG_COL, name = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  LOGX + scale_y_continuous(labels = label_percent(accuracy = 1)) +
  expand_limits(y = 0) +
  labs(x = expression("individual-offspring overdispersion  " * italic(k)),
       y = "mean probability, at-risk zones") +
  theme_pub()

# ---- Panel C: per-zone reach trajectories (top-12 at-risk by baseline) ------
base_rank <- byzone %>% filter(abs(k_indiv - KB) < 1e-9) %>%
  arrange(desc(p_invasion)) %>% slice_head(n = 12) %>% pull(health_zone)
pc_df <- byzone %>% filter(health_zone %in% base_rank) %>%
  mutate(province = ifelse(is.na(province) | !province %in% names(PROV_COL), "Other", province),
         province = factor(province, levels = names(PROV_COL)))
pC <- ggplot(pc_df, aes(k_indiv, p_invasion, group = health_zone, colour = province)) +
  band_layer() + kbase_layer(KB) +
  geom_line(linewidth = 0.45, alpha = 0.9) + geom_point(size = 0.8, alpha = 0.9) +
  scale_colour_manual(values = PROV_COL, name = NULL, drop = TRUE) +
  LOGX + scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(x = expression("individual-offspring overdispersion  " * italic(k)),
       y = "13-week reach P, top-12 zones") +
  theme_pub()

# ---- compose ----------------------------------------------------------------
fig <- (pA | pB) / pC +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(fig, "Figure_kindiv_sweep", w = 7.2, h = 6.4)

# a compact 2-panel variant (A,B side-by-side) for slide/embed use
fig2 <- (pA | pB) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(fig2, "Figure_kindiv_sweep_compact", w = 7.2, h = 3.2)

message("[42_kindiv_sweep_figure] done.")
