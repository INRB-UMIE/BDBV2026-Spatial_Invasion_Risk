# =============================================================================
# 27_mobility_comparison.R — Mobility-kernel comparison visualisations
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Purpose: visualise the DIFFERENCES between every mobility kernel built by
#   03_mobility_matrices.R (M1..M17 and the -dist twins), with an emphasis on the
#   Flowminder additions — the combined inflow+outflow static kernel (M15), the
#   cohort+static composite (M16), and the grand all-kernel consensus ensemble
#   (M17/M17-dist). Reads the saved outputs/mobility/mobility_*.rds files, so it
#   runs after build_all_mobility_matrices() and needs no model fitting.
#
# Panels:
#   A  Kernel similarity heatmap — cosine similarity of the vectorised off-diagonal
#      outflow weights between every pair of kernels (how M15/M16/M17 relate to the
#      established kernels).
#   B  Origin coverage & sparsity — the share of origin zones with any outflow, per
#      kernel (highlights that M15 is materially less sparse than the directed M3).
#   C  Source-zone outflow profiles — top destinations for key invasion source zones
#      (epicentre Bunia; Kisangani hub Makiso Kisangani; Nord-Kivu hub Beni) compared
#      across a curated kernel set {M8, M13, M14, M16, M17}, so the effect of pairing
#      the cohort data with static flows (M16) vs gravity/radiation (M13/M14) is visible.
#   D  Composite divergence — per-source total-variation distance of each cohort
#      composite (M13/M14/M16) from the plain gravity kernel M4, showing how much the
#      cohort/static source rows re-route import pressure.
#
# Usage:
#   source(".../27_mobility_comparison.R"); make_mobility_comparison_figures()
#   # or run standalone:  Rscript spatiotemporal/27_mobility_comparison.R
#
# Outputs:
#   outputs/mobility/mobility_kernel_comparison.pdf (+ .png)
#   outputs/mobility/mobility_kernel_similarity.csv
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a)) a else b

# --- House style (mirrors 17_invasion_viz.R) --------------------------------
.MC_OKABE <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
               "#56B4E9", "#F0E442", "#000000")
.mc_theme <- function(base = VIZ_THEME_BASE %||% 12) {
  ggplot2::theme_minimal(base_size = base, base_family = "sans") +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base + 1),
      plot.subtitle = ggplot2::element_text(colour = "grey38", size = base - 1),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      strip.text = ggplot2::element_text(face = "bold"))
}

# Robust save: base "pdf" device (NOT cairo) for the batch environment, plus a PNG.
.mc_save <- function(p, path, w = 13, h = 10) {
  if (file.exists(path)) try(file.remove(path), silent = TRUE)
  ggplot2::ggsave(path, p, width = w, height = h, device = "pdf", limitsize = FALSE)
  png <- sub("\\.pdf$", ".png", path)
  try(ggplot2::ggsave(png, p, width = w, height = h, dpi = 300, device = "png",
                      limitsize = FALSE), silent = TRUE)
  message(sprintf("[mob-viz] wrote %s", path))
}

#' Load every saved mobility kernel from OUT_MOBILITY.
#' @return named list id -> row-stochastic matrix (aligned to a common zone set).
.mc_load_kernels <- function(out_mobility = OUT_MOBILITY) {
  fs <- list.files(out_mobility, pattern = "^mobility_.*\\.rds$", full.names = TRUE)
  if (!length(fs)) stop("[mob-viz] no mobility_*.rds files in ", out_mobility,
                        " — run 03_mobility_matrices.R / build_all_mobility_matrices() first.")
  ids <- sub("^mobility_(.*)\\.rds$", "\\1", basename(fs))
  ks  <- lapply(fs, readRDS)
  names(ks) <- ids
  # keep only square named matrices on a common zone spine
  ok <- vapply(ks, function(W) is.matrix(W) && nrow(W) == ncol(W) &&
                 !is.null(rownames(W)), logical(1))
  ks <- ks[ok]
  # order ids naturally (M1, M2a, ... then -dist), stable and readable
  ord <- order(match(names(ks), MOBILITY_IDS), nchar(names(ks)), names(ks))
  ks[ord]
}

# Vectorise off-diagonal weights of a kernel (row-major), for similarity.
.mc_vec_offdiag <- function(W) { d <- W; diag(d) <- NA_real_; as.numeric(d) }

#' Panel A — pairwise cosine similarity of the off-diagonal outflow fields.
.mc_panel_similarity <- function(kernels) {
  ids <- names(kernels)
  # common zone set across kernels (they are all 519 in practice, but be defensive)
  zc  <- Reduce(intersect, lapply(kernels, rownames))
  V <- vapply(kernels, function(W) .mc_vec_offdiag(W[zc, zc, drop = FALSE]),
              numeric(length(zc) * length(zc)))
  # cosine similarity between columns, ignoring NA (diagonal) entries pairwise-complete
  cs <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
  for (a in seq_along(ids)) for (b in seq_along(ids)) {
    x <- V[, a]; y <- V[, b]; keep <- is.finite(x) & is.finite(y)
    denom <- sqrt(sum(x[keep]^2)) * sqrt(sum(y[keep]^2))
    cs[a, b] <- if (denom > 0) sum(x[keep] * y[keep]) / denom else NA_real_
  }
  readr::write_csv(as.data.frame(cs) %>% tibble::rownames_to_column("kernel"),
                   file.path(OUT_MOBILITY, "mobility_kernel_similarity.csv"))
  df <- as.data.frame(as.table(cs)); names(df) <- c("k1", "k2", "cos")
  df$k1 <- factor(df$k1, levels = ids); df$k2 <- factor(df$k2, levels = rev(ids))
  ggplot2::ggplot(df, ggplot2::aes(k1, k2, fill = cos)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(option = "mako", limits = c(0, 1), na.value = "grey85",
                                  name = "cosine\nsimilarity") +
    ggplot2::labs(title = "A. Kernel similarity",
                  subtitle = "Cosine similarity of off-diagonal outflow weights",
                  x = NULL, y = NULL) +
    .mc_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}

#' Panel B — origin coverage (share of origins with any outflow) and sparsity.
.mc_panel_coverage <- function(kernels) {
  df <- purrr::imap_dfr(kernels, function(W, id) {
    rs <- rowSums(W)
    tibble::tibble(id = id,
                   origin_coverage = mean(rs > 1e-9),
                   sparsity        = mean(W == 0))
  })
  df$id <- factor(df$id, levels = df$id[order(df$origin_coverage)])
  ggplot2::ggplot(df, ggplot2::aes(origin_coverage, id)) +
    ggplot2::geom_col(fill = .MC_OKABE[1], width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = scales::percent(origin_coverage, accuracy = 1)),
                       hjust = -0.1, size = 3) +
    ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1.08),
                                expand = c(0, 0)) +
    ggplot2::labs(title = "B. Origin coverage",
                  subtitle = "Share of origin zones with any outflow (higher = less sparse)",
                  x = "origins with outflow", y = NULL) +
    .mc_theme()
}

#' Panel C — top-destination outflow profiles for key source zones across kernels.
.mc_panel_source_profiles <- function(kernels,
                                      sources = c("Bunia", "Makiso Kisangani", "Beni"),
                                      show    = c("M8", "M13", "M14", "M16", "M17"),
                                      topn    = 8L) {
  show <- intersect(show, names(kernels))
  zc   <- Reduce(intersect, lapply(kernels[show], rownames))
  sources <- intersect(sources, zc)
  if (!length(show) || !length(sources)) return(NULL)
  rows <- list()
  for (src in sources) for (id in show) {
    w <- kernels[[id]][src, ]
    if (sum(w) <= 0) next
    top <- sort(w, decreasing = TRUE)[seq_len(min(topn, sum(w > 0)))]
    rows[[length(rows) + 1L]] <- tibble::tibble(
      source = src, kernel = id, dest = names(top), weight = as.numeric(top))
  }
  if (!length(rows)) return(NULL)
  df <- dplyr::bind_rows(rows)
  df$kernel <- factor(df$kernel, levels = show)
  # order destinations within each (source) facet by their max weight for readability
  df <- df %>% dplyr::group_by(source, dest) %>%
    dplyr::mutate(.ord = max(weight)) %>% dplyr::ungroup()
  ggplot2::ggplot(df, ggplot2::aes(x = weight,
                                   y = tidytext_reorder(dest, .ord),
                                   fill = kernel)) +
    ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single"),
                      width = 0.75) +
    ggplot2::facet_wrap(~ source, scales = "free", ncol = 1) +
    ggplot2::scale_fill_manual(values = .MC_OKABE[seq_along(show)],
                               labels = mobility_pretty_label(show)) +
    ggplot2::labs(title = "C. Source-zone outflow profiles",
                  subtitle = "Top destinations for key invasion source zones, by kernel",
                  x = "outflow share", y = NULL, fill = "kernel") +
    .mc_theme() +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8))
}

# Small helper: reorder a factor by a numeric key (avoids a hard tidytext dependency).
tidytext_reorder <- function(x, by) factor(x, levels = unique(x[order(by)]))

#' Panel D — per-source total-variation distance of each cohort composite from M4 gravity.
.mc_panel_divergence <- function(kernels,
                                 comps = c("M13", "M14", "M16", "M17"),
                                 ref   = "M4") {
  comps <- intersect(comps, names(kernels))
  if (is.null(kernels[[ref]]) || !length(comps)) return(NULL)
  R <- kernels[[ref]]
  zc <- rownames(R)
  df <- purrr::map_dfr(comps, function(id) {
    W <- kernels[[id]][zc, zc, drop = FALSE]
    # rows where BOTH have outflow (TV distance is only meaningful there)
    both <- rowSums(W) > 1e-9 & rowSums(R) > 1e-9
    tv <- 0.5 * rowSums(abs(W[both, , drop = FALSE] - R[both, , drop = FALSE]))
    tibble::tibble(kernel = id, tv = as.numeric(tv))
  })
  df$kernel <- factor(df$kernel, levels = comps)
  ggplot2::ggplot(df, ggplot2::aes(tv, kernel, fill = kernel)) +
    ggplot2::geom_violin(colour = NA, alpha = 0.5, scale = "width") +
    ggplot2::geom_boxplot(width = 0.18, outlier.size = 0.4, fill = "white") +
    ggplot2::scale_fill_manual(values = .MC_OKABE[seq_along(comps)], guide = "none") +
    ggplot2::scale_y_discrete(labels = function(v) mobility_pretty_label(v)) +
    ggplot2::labs(title = sprintf("D. Divergence from %s (gravity)", ref),
                  subtitle = "Per-source total-variation distance of composite outflow from the gravity kernel",
                  x = "total-variation distance", y = NULL) +
    .mc_theme()
}

#' Build and save the mobility-kernel comparison figure.
#' @param out_mobility directory of mobility_*.rds and the figure output.
#' @return (invisibly) the assembled patchwork object.
make_mobility_comparison_figures <- function(out_mobility = OUT_MOBILITY) {
  message("[mob-viz] Building mobility-kernel comparison figure ...")
  kernels <- .mc_load_kernels(out_mobility)
  message(sprintf("[mob-viz] loaded %d kernels: %s",
                  length(kernels), paste(names(kernels), collapse = ", ")))

  pA <- .mc_panel_similarity(kernels)
  pB <- .mc_panel_coverage(kernels)
  pC <- .mc_panel_source_profiles(kernels)
  pD <- .mc_panel_divergence(kernels)

  # Assemble what is available (panels C/D can be NULL if their inputs are absent).
  top    <- pA + pB + patchwork::plot_layout(widths = c(1.25, 1))
  bottom <- if (!is.null(pC) && !is.null(pD))
              pC + pD + patchwork::plot_layout(widths = c(1, 1))
            else pC %||% pD
  fig <- if (!is.null(bottom)) top / bottom + patchwork::plot_layout(heights = c(1, 1.1))
         else top
  fig <- fig + patchwork::plot_annotation(
    title = "Mobility-kernel comparison — Flowminder static (M15), cohort+static (M16), all-kernel ensemble (M17)",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))

  .mc_save(fig, file.path(out_mobility, "mobility_kernel_comparison.pdf"), w = 14, h = 12)
  invisible(fig)
}

# --- MAIN (standalone run) ---------------------------------------------------
if (identical(environment(), globalenv()) && !exists(".MC_SOURCED_ONLY")) {
  if (sys.nframe() == 0L) {
    tryCatch(make_mobility_comparison_figures(),
             error = function(e) message("[mob-viz] failed: ", conditionMessage(e)))
  }
}
