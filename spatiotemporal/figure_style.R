# =============================================================================
# figure_style.R — the figure design system for the whole suite
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# ONE definition of the house style. Until 2026-09-27 ten files carried their own copies of
# theme_pub()/theme_map()/save_dual() and seven carried their own PROV_COL, and they had
# drifted: the same province was drawn in different colours in different manuscript figures,
# one save_dual() had no FIGURE_KEEP gate, and only three applied the probability-scale
# stamp. Every figure builder now sources this file and defines none of it itself.
#
# Contents
#   * the colourblind-safe palette, including ONE province scale (PROV_COL + prov_palette())
#   * theme_pub() / theme_map(), both title-free: the main text carries the captions
#   * save_dual(), which writes a PDF + 600-dpi PNG through the FIGURE_KEEP gate, the
#     probability-scale stamp and the MI stamp
#
# save_dual() needs an output directory. A suite either sets PANEL_DIR before sourcing this
# file, or passes `dir =` / `dirs =` at the call site; the default is resolved at CALL time,
# so sourcing order does not matter.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. DESIGN SYSTEM  (colourblind-safe)
# -----------------------------------------------------------------------------
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
AFFECTED_FILL <- "grey78"; NA_FILL <- "grey93"
# distinct warm neutral for zones with NO confirmed cases in the manuscript suite's Figure 1A
# bivariate map (make_manuscript_figures.R) — deliberately OFF the grey→red→blue bivariate ramp
# so "no cases" is never read as the low-cases/early-onset corner. (The publication suite's
# Figure 1A is a sequential map and draws its zero class in NA_FILL.)
NO_CASE_FILL <- "#E7DCC1"
RR_FLOOR <- 1e-3
ACCENT   <- "#EE7733"                          # locator-box accent
PT_BLUE  <- "#4C78C8"; FIT_RED <- "#C0392B"    # Figure 1C points / linear fit
OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")

# PROVINCE IDENTITY — ONE scale for every figure in the suite.
#
# A province must keep its hue across Figures 1, 2, 3 and 4 or a reader comparing panels is
# misled. Seven files used to define this vector and four of them disagreed: Sud-Kivu was
# #56B4E9 in the main-text figures and #CC79A7 in Figures 2 and 4, and Tshopo and Bas-Uele
# were permuted again in the kinematics and sweep figures. This assignment is the one the
# main-text manuscript figures already used.
#
# Provinces WITHOUT an entry here (the western/southern urban catchments that only the
# cascade's hub panels reach) are coloured by prov_palette() from SPARE_COL, which is
# deliberately disjoint from these hues — so a spare province never takes an identity
# province's colour. Two SPARE provinces collide only on a panel holding more of them than
# SPARE_COL has hues, where the pool recycles.
PROV_COL <- c("Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73",
              "Tshopo"="#CC79A7","Bas-Uele"="#E69F00","Sud-Kivu"="#56B4E9","Other"="#7A7A7A")
# Paul-Tol muted hues (CVD-safe), disjoint from PROV_COL.
SPARE_COL <- c("#882255","#332288","#661100","#999933","#44AA99","#AA4499","#DDCC77","#117733")
HZ_COL   <- c("1"="#3B4CC0","2"="#B4413C")
base_family <- "sans"

#' A named colour vector covering exactly the provinces supplied.
#'
#' Provinces of interest keep their identity hue and are ordered first; anything else is
#' assigned a SPARE_COL hue, recycled only if that pool is exhausted. NA provinces are
#' dropped rather than coloured, because an unnamed slice is not a province.
prov_palette <- function(provs,
                         prov_int = get0("PROVINCES_OF_INTEREST",
                                         ifnotfound = c("Ituri", "Nord-Kivu", "Haut-Uele"))) {
  provs <- unique(as.character(provs[!is.na(provs)]))
  if (!length(provs)) return(stats::setNames(character(0), character(0)))
  provs <- c(intersect(prov_int, provs), sort(setdiff(provs, prov_int)))
  spare <- setdiff(SPARE_COL, unname(PROV_COL)); si <- 1L
  cols  <- stats::setNames(character(length(provs)), provs)
  for (p in provs) {
    if (p %in% names(PROV_COL)) {
      cols[p] <- unname(PROV_COL[p])
    } else {
      cols[p] <- spare[((si - 1L) %% length(spare)) + 1L]; si <- si + 1L
    }
  }
  stopifnot(!anyNA(cols), all(nzchar(cols)))   # every province must get a visible colour
  cols
}

#' prov_palette() in the list(levels, colours) form a discrete scale and a factor both need.
#'
#' The manuscript suite used to carry its own province_palette() whose spare pool was OKABE
#' minus PROV_COL — two hues (#F0E442, #000000) — so a third province without an identity hue
#' recycled the first spare and two provinces met on one panel in the same colour. It also
#' disagreed with prov_palette() about which hue a given spare province received, so Figure 3
#' and the cascade figures could colour the same province differently. This is the same
#' assignment as prov_palette(), with the level order it uses.
prov_palette_levels <- function(provs, ...) {
  cols <- prov_palette(provs, ...)
  list(levels = names(cols), colours = cols)
}

#' Collapse a province vector to the IDENTITY provinces (those with a PROV_COL hue) + "Other".
#'
#' For panels that place every at-risk zone in the country on one plot (the priority scatters),
#' where giving each of ~25 provinces its own hue would recycle the spare pool. Collapsing on
#' PROVINCES_OF_INTEREST instead — as the scatters used to — drew Tshopo and Bas-Uele zones
#' grey ("Other") in Figure S2B while Figure 3B drew the same zones pink and amber, so one zone
#' changed colour between two figures that sit side by side. Every province with an identity
#' hue now keeps it here too; "Other" takes PROV_COL["Other"].
#' @return list(region = character vector, levels, colours), "Other" last.
prov_region_scale <- function(province) {
  id     <- setdiff(names(PROV_COL), "Other")
  region <- ifelse(!is.na(province) & province %in% id, as.character(province), "Other")
  pres   <- unique(region)
  cols   <- prov_palette(setdiff(pres, "Other"))
  levs   <- c(names(cols), if ("Other" %in% pres) "Other")
  cols   <- c(cols, if ("Other" %in% pres) PROV_COL["Other"])
  list(region = region, levels = levs, colours = cols[levs])
}

# -----------------------------------------------------------------------------
# 2. THEMES
# -----------------------------------------------------------------------------
theme_pub <- function(base = 8.6) {
  ggplot2::theme_minimal(base_size = base, base_family = base_family) %+replace% ggplot2::theme(
    plot.title    = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(),
    plot.caption  = ggplot2::element_blank(),
    axis.title    = ggplot2::element_text(size = base - 0.4, colour = MUTED),
    axis.title.x  = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
    axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(r = 4), angle = 90),
    axis.text     = ggplot2::element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = GRID, linewidth = 0.3),
    legend.position = "top", legend.justification = "left",
    legend.title  = ggplot2::element_text(size = base - 1.2, colour = MUTED),
    legend.text   = ggplot2::element_text(size = base - 1.4, colour = INK),
    legend.key.height = ggplot2::unit(9, "pt"), legend.key.width = ggplot2::unit(15, "pt"),
    strip.text    = ggplot2::element_text(size = base - 0.6, colour = INK, face = "bold",
                                          margin = ggplot2::margin(3, 3, 3, 3)),
    plot.tag      = ggplot2::element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin   = ggplot2::margin(6, 8, 6, 6))
}

theme_map <- function(base = 8.6) {
  ggplot2::theme_void(base_size = base, base_family = base_family) %+replace% ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    legend.position = "right",
    legend.title  = ggplot2::element_text(size = base - 1.4, colour = MUTED),
    legend.text   = ggplot2::element_text(size = base - 1.8, colour = INK),
    legend.key.height = ggplot2::unit(20, "pt"), legend.key.width = ggplot2::unit(7, "pt"),
    plot.tag      = ggplot2::element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin   = ggplot2::margin(2, 2, 2, 2))
}

# -----------------------------------------------------------------------------
# 3. SAVING
# -----------------------------------------------------------------------------

#' The caption a figure currently DISPLAYS, as one string, or NULL when it shows none.
#'
#' A patchwork's caption is its plot_annotation() caption, drawn on a blank ggplot() carrying
#' the annotation theme over the session theme (patchwork's annotate_table()); a ggplot's is its
#' labs(caption), drawn under its own theme. Either is shown only when it is non-empty AND the
#' theme that draws it does not blank plot.caption — theme_pub()/theme_map() blank it, so a
#' caption set on such a panel is computed and never drawn.
.fig_shown_caption <- function(p) {
  if (inherits(p, "patchwork")) {
    txt <- p$patches$annotation$caption; th <- p$patches$annotation$theme
  } else if (inherits(p, "ggplot")) {
    txt <- p$labels$caption; th <- p$theme
  } else return(NULL)
  if (is.null(txt) || !is.character(txt) || !length(txt) || !any(nzchar(txt))) return(NULL)
  el <- tryCatch(ggplot2::calc_element("plot.caption", ggplot2::complete_theme(th)),
                 error = function(e) NULL)
  if (is.null(el) || inherits(el, "element_blank")) return(NULL)
  paste(txt, collapse = "\n")
}

#' Apply a caption stamp (fs_caption(), .fig_mi_stamp()) so that it is APPENDED to a caption
#' the figure already shows, instead of replacing it.
#'
#' Both stampers set the caption outright — `plot_annotation(caption = )` on a patchwork,
#' `labs(caption = )` on a ggplot — so on a figure that carries its own caption (the decoding
#' caption of Figure2_labelled) a scale or MI stamp silently DELETED it. Here the shown caption
#' is read first; if the stamper changed it, the result carries the old caption, then the stamp,
#' drawn in the OLD caption's style (the stamper's small grey style would otherwise restyle the
#' figure's own caption). A caption that was not shown (a blanked theme) is not resurrected:
#' the stamp then stands alone, exactly as the stamper drew it.
.fig_stamp_append <- function(p, stamper) {
  old <- .fig_shown_caption(p)
  q   <- stamper(p)
  if (is.null(old)) return(q)
  new <- .fig_shown_caption(q)
  if (is.null(new) || identical(new, old)) return(q)        # the stamper did nothing
  both <- paste(old, new, sep = "\n")
  if (inherits(q, "patchwork")) {
    q$patches$annotation$caption <- both
    q$patches$annotation$theme   <- p$patches$annotation$theme
    return(q)
  }
  old_el <- ggplot2::calc_element("plot.caption", ggplot2::complete_theme(p$theme))
  q + ggplot2::labs(caption = both) + ggplot2::theme(plot.caption = old_el)
}

#' Stamp an MI-pooled figure so the PDF says so once it leaves its folder.
#'
#' An MI-pooled figure shares its basename with the published single-draw twin and its
#' intervals mean something different (posterior x process x imputation, Rubin-pooled).
#' No-op unless MI_STAMP_TXT is set. Sets the caption outright; save_dual() applies it through
#' .fig_stamp_append(), so a caption the figure already shows is kept.
.fig_mi_stamp <- function(p) {
  mi <- get0("MI_STAMP_TXT", ifnotfound = NULL)
  if (is.null(mi) || !nzchar(mi) || is.null(p)) return(p)
  th <- ggplot2::theme(plot.caption = ggplot2::element_text(size = 6.5, colour = MUTED,
                                                            hjust = 0))
  if (inherits(p, "patchwork") && requireNamespace("patchwork", quietly = TRUE))
    return(p + patchwork::plot_annotation(caption = mi, theme = th))
  if (inherits(p, "ggplot")) return(p + ggplot2::labs(caption = mi) + th)
  p
}

#' Write a figure as PDF + 600-dpi PNG into one or more directories.
#'
#' THE ONE implementation. Three behaviours used to be present in some copies and absent in
#' others, each with a consequence: the FIGURE_KEEP gate (without it a suppressed figure is
#' published anyway), fs_caption() (without it a raw-scale twin carries no scale statement
#' while sharing a basename with the recalibrated one), and dir.create() (without it the
#' first write into a new directory errors).
#'
#' @param p a ggplot or patchwork object.
#' @param name figure stem, WITHOUT extension — what FIGURE_KEEP matches on.
#' @param w,h inches.
#' @param dirs one or more output directories. Defaults to `dir`, then SAVE_DIRS, then
#'   PANEL_DIR — all resolved at CALL time.
#' @param dir  convenience alias for a single directory.
#' @param context passed to fs_caption(): "lfo", "current" or "mixed".
save_dual <- function(p, name, w, h, dir = NULL, dirs = NULL, context = "mixed") {
  # Directory resolution, in order: an explicit `dirs=`, an explicit `dir=`, the suite-level
  # SAVE_DIRS (set by the builders that write the SAME figure into two trees), then PANEL_DIR.
  if (is.null(dirs)) dirs <- dir
  if (is.null(dirs)) dirs <- get0("SAVE_DIRS",  ifnotfound = NULL)
  if (is.null(dirs)) dirs <- get0("PANEL_DIR",  ifnotfound = NULL)
  if (is.null(dirs) || !length(dirs))
    stop("[fig] save_dual(): no output directory (pass dir=/dirs=, or set SAVE_DIRS or ",
         "PANEL_DIR before calling).", call. = FALSE)
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): skip any figure that is not on the
  # published allow-list. Gate on the FULL destination path, not the bare stem: FIGURE_DROP
  # entries are "<directory>/<stem>" and the raw/ exclusion inspects path components, neither
  # of which can match a basename. A figure written to several directories is kept if ANY
  # destination is allowed, and is then written only to the allowed ones.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk)) {
    keep <- vapply(dirs, function(d) isTRUE(.fk(file.path(d, name))), logical(1))
    if (!any(keep)) return(invisible(p))
    dirs <- dirs[keep]
  }
  # Both stamps APPEND to a caption the figure already shows (.fig_stamp_append()): applied
  # raw, each replaced it, so a scale or MI stamp deleted Figure2_labelled's decoding caption.
  if (exists("fs_caption", mode = "function"))
    p <- .fig_stamp_append(p, function(x) fs_caption(x, context = context))
  p <- .fig_stamp_append(p, .fig_mi_stamp)
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(d, paste0(name, ".pdf")), p, width = w, height = h,
                    device = figure_pdf_device(), bg = "white")
    ggplot2::ggsave(file.path(d, paste0(name, ".png")), p, width = w, height = h,
                    dpi = 600, bg = "white")
  }
  message(sprintf("  saved %-38s %.1f x %.1f in", name, w, h))
  invisible(p)
}

# -----------------------------------------------------------------------------
# 4. SMALL PLOTTING UTILITIES (one definition; four files used to carry copies)
# -----------------------------------------------------------------------------

#' Order a discrete axis WITHIN each facet (the tidytext idiom, without the dependency).
#' Pair with tidytext_scale_y(), which strips the disambiguating suffix off the labels.
reorder_within <- function(x, by, within) {
  factor(paste(x, within, sep = "___"),
         levels = unique(paste(x, within, sep = "___"))[order(within, by)])
}
tidytext_scale_y <- function()
  ggplot2::scale_y_discrete(labels = function(z) sub("___.*$", "", z))

#' Regularly spaced date breaks for a calendar axis, at most `max_n` of them.
#'
#' A fixed `date_breaks = "2 weeks"` overprints once the outbreak runs long enough: on the
#' 2026-10 runs it put twelve "%d %b" labels on a ~6 in Figure 1B axis and they ran into one
#' another ("04 May18 May01 Jun..."). This keeps WHOLE-WEEK spacing anchored on the axis start,
#' so the ticks still fall on the same weekday as the week grid, and widens the step as the
#' span grows rather than as somebody remembers to.
#' @param lims length-2 Date vector (axis limits).
#' @param max_n the most labels the axis can carry at its width and font size.
fs_date_breaks <- function(lims, max_n = 6L) {
  lims <- sort(as.Date(lims))
  wks  <- as.numeric(diff(lims)) / 7
  # n breaks span n - 1 steps, so the step is set against max_n - 1 intervals.
  step <- max(1L, as.integer(ceiling(wks / max(1L, max_n - 1L))))
  seq(lims[1], lims[2], by = 7L * step)
}

#' Canonical zone-name join key: case- and whitespace-insensitive.
key <- function(x) tolower(trimws(x))

#' Rank-based quartile, 1..4. NA in, NA out; all-NA in, all-NA out.
qtile <- function(x) {
  out <- rep(NA_integer_, length(x)); ok <- is.finite(x)
  if (!any(ok)) return(out)
  r <- rank(x[ok], ties.method = "average") / sum(ok)
  out[ok] <- as.integer(cut(r, c(-Inf, .25, .5, .75, Inf), labels = 1:4))
  out
}

#' Join a zone-keyed table onto the shapefile: every polygon kept, each data row attached to
#' the ONE polygon it names.
#'
#' Each data row is matched on (name, province) first; a row that does not match that way —
#' its province is missing, NA, or spelled differently from the shapefile's — is matched on the
#' name alone, but only where that name identifies exactly one polygon. Province is part of the
#' key so that a zone NAME recurring across provinces can never be attached to the wrong one.
#' (The current shapefile's names are unique: the two zone pairs that share a base name, Bili
#' and Lubunga, carry a province qualifier in the name itself.)
#'
#' The match is decided PER DATA ROW. It used to be decided for the whole table from the share
#' of shapefile POLYGONS left without data, which for any sparse table (Figure 1's invaded zones
#' are a small fraction of the country) exceeded its threshold, so every such table was joined
#' on the name alone and a row carrying the wrong province was accepted without a word; and a
#' table with no province column errored instead of falling back. A row matched on its name
#' against a different province is now reported, and a row that names no polygon is reported
#' too: it is not drawn, and nothing else says so.
#'
#' @param dat a data frame with `health_zone` and (ideally) `province`.
#' @param shp the sf object to join onto, carrying the `.key` (lower-cased, trimmed name) and
#'   `.prov` columns every suite adds when it loads it; resolved from the caller's scope by
#'   default, because each figure suite loads its own shapefile once.
#' @return `shp` with `dat`'s columns joined (NA for polygons no row names).
join_map <- function(dat, shp = get0("shp", envir = parent.frame())) {
  if (is.null(shp)) stop("[fig] join_map(): no shapefile in scope (pass shp=).", call. = FALSE)
  if (!all(c(".key", ".prov") %in% names(shp)))
    stop("[fig] join_map(): the shapefile lacks its .key/.prov columns.", call. = FALSE)
  .k <- function(z) tolower(trimws(z))
  s_key  <- as.character(shp$.key); s_prov <- as.character(shp$.prov)
  d_key  <- .k(dat$health_zone)
  d_prov <- if ("province" %in% names(dat)) as.character(dat$province)
            else rep(NA_character_, nrow(dat))
  # (name, province) first; then the name alone, where it names exactly one polygon.
  by_both <- match(paste(d_key, d_prov, sep = "\r"), paste(s_key, s_prov, sep = "\r"))
  by_both[is.na(d_prov)] <- NA_integer_
  uniq    <- names(which(table(s_key) == 1L))
  by_name <- ifelse(d_key %in% uniq, match(d_key, s_key), NA_integer_)
  poly    <- ifelse(is.na(by_both), by_name, by_both)
  .zl <- function(i) paste(utils::head(unique(as.character(dat$health_zone[i])), 6),
                           collapse = ", ")
  off <- which(is.na(by_both) & !is.na(by_name) & !is.na(d_prov))
  if (length(off))
    message(sprintf(paste0("[fig] join_map(): %d row(s) matched on the zone name alone, their ",
                           "province differing from the shapefile's (%s)."), length(off), .zl(off)))
  if (anyNA(poly))
    warning(sprintf(paste0("[fig] join_map(): %d row(s) name no single shapefile zone and are ",
                           "not drawn (%s)."), sum(is.na(poly)), .zl(is.na(poly))), call. = FALSE)
  # Two rows on one polygon would draw it twice, the later fill over the earlier.
  dup <- which(!is.na(poly) & (duplicated(poly) | duplicated(poly, fromLast = TRUE)))
  if (length(dup))
    warning(sprintf(paste0("[fig] join_map(): %d row(s) name a zone another row also names, so ",
                           "its polygon is drawn once per row (%s)."), length(dup), .zl(dup)),
            call. = FALSE)
  d <- dat; d$.poly <- poly
  d <- d[!is.na(d$.poly), setdiff(names(d), c(".key", ".prov")), drop = FALSE]
  m <- shp; m$.poly <- seq_len(nrow(m))
  m <- dplyr::left_join(m, d, by = ".poly")
  m$.poly <- NULL
  m
}

# -----------------------------------------------------------------------------
# 5. PANEL-A COMPARATOR IDENTITY  (Figure 2 and Figure 2-cascade)
# -----------------------------------------------------------------------------
# ONE definition of who the prioritisation panel's comparators ARE, what they are CALLED and
# what COLOUR they are drawn in. Figure 2 (make_manuscript_figures.R) and its cascade analogue
# (make_manuscript_figure2_cascade.R) are the same panel over two horizons — a 2-week forecast
# and a 6-week leave-future-out cascade backtest — and a reader is meant to compare them
# directly. They had drifted apart: the cascade figure carried its own private set of THREE
# baselines (gravity / Flowminder inflow / travel time) against Figure 2's FOUR, with gravity
# drawn in #009E73 — the green Figure 2 gives to "Nearest affected zone". So the same hue named
# two different comparators across two figures that sit side by side in the manuscript, and the
# cascade panel silently omitted both case-using nulls.
#
# The MEMBERSHIP is FIGURE_BASELINE_METHODS (00_config.R); the map below carries labels and
# colours for every comparator the suite can publish, including ones not currently selected, so
# it survives a change to that config vector.

BASELINE_LBL <- c(
  # M4: fitted gravity, log(flow) = b0 + a log(pop_i) + b log(pop_j) - c log(d+1), with the
  # exponents estimated on the Flowminder relocation OD under a censored likelihood (build_M4(),
  # 03_mobility_matrices.R; the fitted values are in that module's log).
  "Baseline-gravity"           = "Gravity model (fitted flows)",
  # M_cohort: Flowminder COHORT subscriber-day presence from the epicentre — measured during
  # the outbreak and reaching 305 destinations against the short-trip annex's 142. Pooled over
  # the cohort's origin zones (identical rows), as the short-trip annex is.
  "Baseline-flowminder-inflow" = "Flowminder cohort inflow from epicentre",
  # OSRM road travel time, 1/(1 + minutes) to the nearest epicentre zone. Pure geography.
  "Baseline-travel-time"       = "Road travel time from epicentre",
  # THE TWO CASE-USING NULLS (added 2026-09-25 with FIGURE_BASELINE_METHODS). Both read the
  # case data to locate the front but carry no mobility network, so the panel spans
  # (uses case data) x (uses phone data) rather than only the case-blind row.
  # B1: neighbours' confirmed cases weighted by 1/(travel time), alpha = 1.
  "Distance-B1"                = "Neighbours' cases / travel time",
  # B7 (method = "distance"): 1/(1 + travel time to the nearest already-affected zone).
  "Adjacency-B7"               = "Nearest affected zone")

# Okabe-Ito, all distinguishable from the featured series (FEATURED_COL) and from the grey
# random-targeting reference. Dark red is free for gravity because it is not currently
# published; if it returns to FIGURE_BASELINE_METHODS it must NOT take the green, which
# belongs to Adjacency-B7 in both figures.
BASELINE_COL <- c("Baseline-gravity"           = "#661100",
                  "Baseline-flowminder-inflow" = "#D55E00",
                  "Baseline-travel-time"       = "#CC79A7",
                  "Distance-B1"                = "#E69F00",
                  "Adjacency-B7"               = "#009E73")

#' The featured series' colour, in every prioritisation panel in the suite.
FEATURED_COL <- "#0072B2"

#' What the featured model is CALLED in the manuscript suite's panels that name it by role —
#' Figure 2 panel A's legend, the mobility-value ladder's featured row (Figure 2 panel C,
#' FigureS24-S26) and FigureS5. (The publication suite names it by its method code, or by
#' model_pretty_label() in Figure2_labelled, not by this label.) The
#' ladder used to call the same model "Cohort + radiation + covariates (featured)" one panel
#' below a legend that called it "Mobility-informed model", so a reader comparing panels A and C
#' had nothing on the page saying they were one model. The plain twin is the simplified
#' variants' wording, used by panel A AND by panel C's first row of Figure2_simple.
FEATURED_LBL       <- "Mobility-informed model"
FEATURED_LBL_PLAIN <- "Forecasting model (uses human mobility data)"

# PLAIN-LANGUAGE LABELS, for the simplified variants of Figures 2 and 2-cascade.
#
# Same comparators, same colours, same order — only the words change. The published figures
# name each rule by its mechanism ("Flowminder cohort inflow from epicentre"), which is precise
# but assumes the reader knows that Flowminder is the mobile-network operator panel and that
# "inflow" is a directed movement volume. These say the same things in the terms the manuscript
# text uses, and in particular they say "human mobility data" rather than naming the vendor or
# the medium.
#
# EACH ONE IS THE SAME OBJECT AS ITS BASELINE_LBL TWIN, not a looser paraphrase:
#   Baseline-flowminder-inflow  ranks zones by measured movement from the epicentre with no
#     model of any kind fitted on top — hence "alone, no model". It is the SAME method the
#     simplified panel C draws as "Human mobility data alone, no model", and the two panels
#     deliberately use one wording for it.
#   Baseline-travel-time  is road travel time, in minutes, to the nearest epicentre zone.
#   Distance-B1  weights the confirmed cases in other zones by 1 / (travel time).
#   Adjacency-B7  is 1 / (1 + travel time to the nearest zone that already has cases).
#   Baseline-gravity  is the fitted gravity kernel M4 — whose deterrence exponents are
#     estimated FROM the mobility data, so it is not a mobility-free comparator and its plain
#     label must not imply that it is.
# ONLY THE JARGON CHANGES. "Epicentre", "road travel time" and "affected zone" are already
# plain English and are also SHORTER than any paraphrase of them, so they stay: this legend
# sits over a ~3 in panel and its width is set by its longest label (.assert_legend_fits()
# measures it at build time and stops the build when it does not fit). What these replace is
# the vendor name and the movement-data vocabulary — "Flowminder cohort inflow" — which a
# reader cannot decode without the methods section.
BASELINE_LBL_PLAIN <- c(
  "Baseline-gravity"           = "Modelled movement from epicentre",
  "Baseline-flowminder-inflow" = "Human mobility data alone, no model",
  "Baseline-travel-time"       = "Road travel time from epicentre",
  "Distance-B1"                = "Nearby cases, weighted by travel time",
  "Adjacency-B7"               = "Travel time to nearest affected zone")

#' The comparator set the prioritisation panels publish.
#'
#' Read from 00_config.R so the manuscript suite, the publication suite and the cascade
#' analogue cannot drift apart; the literal is the fallback for a standalone run without
#' the config loaded.
fs_baseline_methods <- function()
  get0("FIGURE_BASELINE_METHODS",
       ifnotfound = c("Baseline-flowminder-inflow", "Baseline-travel-time",
                      "Distance-B1", "Adjacency-B7"))

#' scale_colour_manual() values and labels for a prioritisation panel.
#'
#' @param featured     the featured series' KEY (a method name, or any tag the panel uses).
#' @param featured_lbl what the legend calls it ("Mobility-informed model", "3-month cascade").
#' @param methods      the comparator keys to publish; FIGURE_BASELINE_METHODS by default.
#' @param who          caller tag for the error message.
#' @param plain        TRUE for the plain-language labels the simplified figure variants use
#'   (BASELINE_LBL_PLAIN). The MEMBERSHIP and the COLOURS are identical either way, so a
#'   simplified variant and its published twin differ only in wording.
#' @return list(values, labels), named and in legend order (featured first).
fs_panelA_scales <- function(featured, featured_lbl, methods = fs_baseline_methods(),
                             who = "fig2", plain = FALSE) {
  lbl_map <- if (plain) BASELINE_LBL_PLAIN else BASELINE_LBL
  # The featured series is prepended to the comparators, so it must not BE one of them: a
  # repeated name would give scale_colour_manual() two values and two labels under one key,
  # and it silently keeps the first — so the panel would draw the comparator's colour for the
  # featured model, or the other way round, with nothing to show for it.
  if (featured %in% methods)
    stop("[", who, "] the featured series '", featured, "' is also in the comparator set, so ",
         "the legend would carry two entries under one name.", call. = FALSE)
  # Fail with the missing NAME rather than with ggplot's "`breaks` and `labels` must have the
  # same length", which says nothing about which comparator was not mapped.
  unmapped <- setdiff(methods, intersect(names(BASELINE_COL), names(lbl_map)))
  if (length(unmapped))
    stop("[", who, "] FIGURE_BASELINE_METHODS names a comparator with no label/colour: ",
         paste(unmapped, collapse = ", "),
         ". Add it to BASELINE_LBL, BASELINE_LBL_PLAIN and BASELINE_COL in figure_style.R.",
         call. = FALSE)
  # Subset to the PUBLISHED comparators on BOTH sides: BASELINE_LBL deliberately carries
  # entries that are not currently selected, and leaving them in makes `labels` longer than
  # `values` — the "`breaks` and `labels` must have the same length" abort.
  list(values = c(stats::setNames(FEATURED_COL, featured), BASELINE_COL[methods]),
       labels = c(stats::setNames(featured_lbl,  featured), lbl_map[methods]))
}

#' The in-panel metric block, identical in wording and order across the prioritisation panels.
#'
#' TWO DISCRIMINATION MEASURES, THEN TWO CALIBRATION MEASURES. AUC-PR skill and AUC-ROC are
#' invariant to any monotone transform of the probabilities, so a model can lead on both while
#' its probabilities are wrong by a factor of two — and the probabilities are used to decide how
#' many zones to prepare, not only which ones. Calibration (pred/obs) is calibration-in-the-
#' large: 1.0 = the expected invasion count matched the realised one over the cross-validation.
#'
#' AUC-PR SKILL IS PER ROUND: average precision within each forecast round, combined over rounds
#' weighted by their first cases, divided by the (event-weighted) base rate — the selection's axis
#' (16_invasion_eval.R, .auc_pr_by_date()). Its RANDOM-ORDERING level (`chance`) is printed beside
#' it: a skill of 1 is a constant forecast, and with few first cases per round a list in random
#' order scores well above 1, so the reference a reader needs is that level, not 1.
#'
#' THE BASE RATE SITS UNDER THE SKILL because it puts the multiple in context: without it the
#' reader cannot tell that this is a rare-event problem, and has no way to see why the Brier score
#' below is small for a reason that has nothing to do with the model. It is the share of scored
#' zone-rounds that were positive on the SAME rounds as the skill (the header says which).
#'
#' @param header  the panel's first line(s), naming the evaluation and the rounds each line covers
#'   ("2-week-ahead forecast, ...", "Leave-future-out (+6 wk, 3 origins)").
#' @param chance  the per-round skill of a random ordering on the same rounds.
fs_metric_block <- function(header, skill, chance, base_rate, roc, brier, calib)
  sprintf(paste0("%s\nAUC-PR skill, per round  %.2fx (random %.2fx)\nBase rate  %.2f%%\n",
                 "AUC-ROC  %.2f\nBrier score  %.3f\nCalibration (pred/obs)  %.2f"),
          header, skill, chance, 100 * base_rate, roc, brier, calib)

# -----------------------------------------------------------------------------
# 6. LEGEND-FIT CHECK  (any composed figure whose panels carry their own legends)
# -----------------------------------------------------------------------------
# Figure 2A's legend acquired a fifth entry and, laid out in two rows, measured 7.46 in against
# the ~4.1 in panel A is given. ggplot2 neither clips nor wraps a guide box that is wider than
# its plot — it simply draws past the edge — so the legend ran across panel B and interleaved,
# line for line, with panel B's own legend. The figure rendered without a warning and was
# published that way. Nothing about that is specific to Figure 2: any composed figure whose
# panels carry their own legends can acquire it from one added series or one lengthened label.
# So it is checked, in inches, at build time.
#
# THE GEOMETRY. ggplot2 (>= 3.5) places a legend relative to the PANEL area by default
# (`legend.location = "panel"`), so the width to beat is the panel's, not the plot's. In a
# patchwork row the panels are null units: the fixed parts of every plot in the row (axis
# titles, tick labels, margins) come off the total first, and what is left is divided in the
# ratio given to plot_layout(widths = ). That is why both plots in the row have to be passed
# together.

#' A plot's inelastic width, and how much of it sits LEFT of its first panel, in inches.
#'
#' Do not test a gtable width for "null"-ness and skip it: under ggplot2 4 the whole y axis is
#' ONE packed column whose width is `sum(0cm, 1null, 5.25cm, 0null, ...)`, so a null test throws
#' away the 5.25 cm of axis text with it — which is exactly the 2 in of health-zone labels this
#' check exists to account for. convertWidth() already does the right thing: a null contributes
#' 0 inches, everything else contributes its size. So sum them all.
#'
#' `left` is what the legend is pushed right by, since a top legend starts at the first panel's
#' left edge; for a faceted plot only the FIRST facet's axis column counts, because the rest are
#' inside the area the legend spans.
.gt_widths_in <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  w <- g$widths
  inch <- vapply(seq_along(w), function(i) grid::convertWidth(w[i], "in", valueOnly = TRUE),
                 numeric(1))
  pk <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  lft <- if (nrow(pk)) min(pk$l) else length(w) + 1L
  c(fixed = sum(inch), left = sum(inch[seq_len(max(0L, lft - 1L))]))
}

#' Width of a plot's legend box in inches, 0 when it has none.
.guide_box_width_in <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  i <- grep("^guide-box", g$layout$name)
  if (!length(i)) return(0)
  max(0, vapply(i, function(k) {
    gb <- g$grobs[[k]]
    if (is.null(gb) || inherits(gb, "zeroGrob")) 0
    else grid::convertWidth(grid::grobWidth(gb), "in", valueOnly = TRUE)
  }, numeric(1)))
}

#' How much room each plot's legend has, in inches, for a row laid out at `widths`.
#'
#' A plot's region is its own inelastic width plus its share of what is left over; the legend
#' starts `left` into that region and may run to its right-hand edge.
.legend_room_in <- function(plots, widths, total_in, who = "fig") {
  gw    <- vapply(plots, .gt_widths_in, numeric(2))
  free  <- total_in - sum(gw["fixed", ])
  if (free <= 0)
    stop("[", who, "] the axis labels alone are wider than the figure; nothing is left for ",
         "panels.", call. = FALSE)
  unname(gw["fixed", ] + free * widths / sum(widths) - gw["left", ])
}

#' Stop if any plot's legend is wider than the panel it will be drawn over.
#'
#' A measuring device has to be open for convertWidth() to resolve the text grobs, and under
#' Rscript there is none until ggsave opens one — hence the pdf(NULL).
.assert_legend_fits <- function(plots, widths, total_in, labels = names(plots), tol = 0.02,
                                who = "fig2") {
  grDevices::pdf(NULL, width = total_in, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  room <- .legend_room_in(plots, widths, total_in, who = who)
  leg  <- vapply(plots, .guide_box_width_in, numeric(1))
  for (i in seq_along(plots))
    message(sprintf("[%s] legend fit %-4s legend %5.2f in, room %5.2f in  %s",
                    who, labels[i], leg[i], room[i], if (leg[i] <= room[i] + tol) "ok" else "OVERFLOW"))
  bad <- which(leg > room + tol)
  if (length(bad))
    stop(sprintf(paste0("[%s] panel %s has a %.2f in legend in %.2f in of room, so it would be ",
                        "drawn across its neighbour (this is exactly how panels A and B came to ",
                        "share two legend lines). Give it more rows — guide_legend(ncol = ) — ",
                        "or shorten its labels."),
                 who, paste(labels[bad], collapse = "/"), max(leg[bad]), min(room[bad])),
         call. = FALSE)
  invisible(room)
}

# -----------------------------------------------------------------------------
# 7. THE FEATURED MODEL  (read from the pipeline's selection record, never typed)
# -----------------------------------------------------------------------------

#' The featured Bayesian model of a run, exactly as write_model_selection.R recorded it.
#'
#' `$featured$bayesian` first, `$featured$headline` only when that is absent (an older file):
#' write_model_selection.R now records the headline as the featured model itself, and every
#' figure that calls this is a Bayesian product.
#'
#' NO FALLBACK MODEL. Four figure suites used to fall back to a typed model name when the
#' selection file was unreadable, and Figure 2's mobility ladder typed its featured model
#' outright (MV_FEATURED). When the selection moved from Bayes-M14-fill-geo to
#' Bayes-M14-fill-med the typed name stayed behind and Figure 2 stopped building. A figure that
#' cannot tell which model the pipeline chose must fail, not guess.
#'
#' @param sel_path path to the run's key_outputs/model_selection.json.
#' @return the featured method label, e.g. "Bayes-M14-fill-med".
read_featured_model <- function(sel_path) {
  if (!file.exists(sel_path))
    stop("[featured] ", sel_path, " is missing. The figures name the featured model and must ",
         "not guess it; run the pipeline (write_model_selection.R) first.", call. = FALSE)
  sel <- tryCatch(jsonlite::fromJSON(sel_path, simplifyVector = TRUE),
                  error = function(e) stop("[featured] cannot parse ", sel_path, ": ",
                                           conditionMessage(e), call. = FALSE))
  pick <- function(x) if (is.null(x) || !length(x) || is.na(x[1]) || !nzchar(x[1])) NULL else as.character(x[1])
  m <- pick(sel$featured$bayesian$method)
  if (is.null(m)) m <- pick(sel$featured$headline$method)
  if (is.null(m))
    stop("[featured] ", sel_path, " names no featured model.", call. = FALSE)
  m
}
