# =============================================================================
# mobility_value_panel.R — the mobility-value ladder, as a DRAWABLE PANEL
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# WHY THIS FILE EXISTS. The ladder "what does the mobile-phone mobility data add, over and
# above population and geography?" is drawn in two places:
#
#   * FigureS24 / FigureS25 / FigureS26   (50_mobility_value.R, the appendix module that
#                                          computes the ladder and the paired bootstrap)
#   * Figure 2, panel C                   (make_manuscript_figures.R, main text, 2 weeks
#                                          ahead only — the horizon panels A and B are fixed
#                                          to)
#
# A main-text panel and its appendix twin that are built by two pieces of code can disagree
# — about which rows are shown, how they are ordered, what the colours mean, or where the
# "best phone-free" line sits — and nothing in the suite would notice. So the row set, the
# labels, the two encodings, the ordering rule and the panel geometry live HERE, once, and
# both callers draw through mobval_panel(). The only things the main-text call changes are the
# horizon it asks for, the base font size, whether the single-facet strip is replaced by a
# horizon-bearing axis title, and the number of x-axis breaks (`n_breaks`, for its full-width
# panel).
#
# WHAT IT NEEDS FROM THE CALLER. The caller must already have sourced figure_style.R (for
# theme_pub(), OKABE, INK/MUTED/FAINT/GRID); this file sources it only if that has not
# happened. It defines objects and functions and does no work at source time — in
# particular it does NOT set PANEL_DIR, which is the caller's business and which
# save_dual() resolves at call time.
#
# HOUSE STYLE (figure_style.R): no plot title, no subtitle, no in-figure caption. theme_pub()
# blanks all three, so a stray labs() cannot reintroduce one, and the manuscript carries the
# caption. Everything the panel must say for itself it says through the axis title, the facet
# strip (which horizon), the legend and the two reference lines.
#
# ORDERING. Each ladder is ordered by the metric it draws, at the LONGER horizon, and the
# same order is used for both horizons so a row sits at the same height in both. The two
# metrics therefore order the rows differently in FigureS26 — deliberately: they are two
# ladders, each tagged and labelled on its own, not one ladder drawn twice, and the
# disagreement between them (the case-blind inflow rule) is a result, not a blemish to be
# sorted away.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

# theme_pub(), the Okabe-Ito palette and the ink/grid greys. Sourced only if the caller has
# not already done so: both callers do, and figure_style.R is cheap and side-effect-free,
# but a standalone `source("mobility_value_panel.R")` should still work.
if (!exists("theme_pub", mode = "function") || !exists("FEATURED_LBL"))
  source(file.path(here::here(), "spatiotemporal", "figure_style.R"))

# Colour = what the mobility term is BUILT FROM: the two phone-free classes take two cool
# tones and the two phone classes two warm ones. There is no grey class: every row of the
# ladder has a movement or geography term, and Adjacency-B7 (nearest affected zone) is
# "Geography only" like the other travel-time rules.
# The cool/warm split is the encoding, not decoration: a reader who takes in nothing but the
# colours still sees that the top of the ladder is warm (phone) and the middle is cool (no
# phone). That is worth more than matching Figure 2 panel A's method colours, which cannot be
# done anyway - panel A distinguishes five METHODS, three of which this panel groups into one
# class - and which would also break panel C away from FigureS24. Figure 2 already carries
# three per-panel colour meanings (methods in A, outcomes in B, classes here), each with its
# own legend directly above its own panel; see the block above build_fig2c in
# make_manuscript_figures.R.
#
# BUT NO HUE MAY MEAN ONE THING IN PANEL A AND ANOTHER HERE. This palette used to be Okabe-Ito
# throughout, and panel A draws in Okabe-Ito too: #0072B2 was the FEATURED MODEL in panel A and
# "Geography only" here (the class of the dashed best-phone-free line, i.e. the method the
# featured model is measured AGAINST), and #E69F00 was Distance-B1 in panel A and
# "Phone-calibrated decay" here. The two phone-free-or-calibrated classes therefore take Paul Tol
# hues that panel A never uses (indigo #332288, wine #882255). The one hue the panels SHARE,
# vermillion #D55E00, means the same thing in both: panel A gives it to the Flowminder cohort
# inflow, which is measured phone movement, and that is the class it names here. Checked below
# against FEATURED_COL and BASELINE_COL so a palette edit cannot quietly reintroduce the clash.
MV_CLS <- c("Geography only", "Geography + population",
            "Phone-calibrated decay", "Measured phone movement")
MV_PAL <- setNames(c("#332288", OKABE[6], "#882255", OKABE[2]), MV_CLS)
local({
  shared <- intersect(toupper(MV_PAL), toupper(c(FEATURED_COL, unname(BASELINE_COL))))
  ok <- toupper(c(OKABE[2]))   # vermillion: measured phone movement in both panels (see above)
  if (length(setdiff(shared, ok)))
    stop("[mobval] MV_PAL reuses panel-A colour(s) ", paste(setdiff(shared, ok), collapse = ", "),
         " for a different meaning; Figure 2 would give one hue two meanings.", call. = FALSE)
})
# Shape = WHAT THE METHOD IS, an orthogonal axis colour cannot carry. Without it the
# case-blind epicentre rules sit among the fitted renewal models with nothing to separate
# them, and a reader could take the cohort-inflow row for a model.
MV_KIND <- c("Fitted renewal model (uses cases)", "Structural rule (uses cases)",
             "Case-blind ranking (no cases)")
MV_SHP  <- setNames(c(19, 17, 22), MV_KIND)          # solid circle / solid triangle / hollow square
# THE TWO REFERENCE LINES MUST NOT READ AS ONE. Both were grey, one dashed and one dotted, and
# at this line weight the dash pattern alone did not separate them - a reader met two similar
# grey rules and had to work out which was which. They are now told apart by colour as well,
# and the colours are not arbitrary:
#
#   random ordering (dotted) stays the house FAINT grey. It is a property of the AXIS, not of
#                            any method - the per-round skill a uniformly random ordering scores
#                            on these rounds (`chance` in the ladder CSV, ~2 rather than 1:
#                            with few invasions per round a random list's AP sits well above
#                            the base rate) - so it should read as part of the furniture.
#   best phone-free (dashed) takes THE CLASS COLOUR OF THE METHOD IT IS DRAWN AT. It marks one
#                            specific method - the best of EVERY method on the ladder that uses
#                            no mobile-phone data, of any kind (see mv_ladder_rows()) - so a
#                            reader can see which KIND of phone-free method the model is
#                            measured against without reading down the rows to find it. On
#                            LINELIST_09102026, at both horizons, that is road travel time
#                            from the epicentre on the skill ladder and travel-time decay with
#                            covariates (Bayes-M6a-geo) on the rank ladder, both "Geography
#                            only".
#
# The class colour is DERIVED in mobval_panel(), never hard-coded: if the best phone-free
# method ever changed class, a fixed colour would go on claiming the old one. Only the label
# stays MUTED grey - #56B4E9 is too light for 6-point italic text, and the label sits on its
# own line, so nothing is ambiguous about which rule it names.
MV_NOSKILL_COL <- FAINT
MV_KIND_MAP <- setNames(MV_KIND,
                        c("Fitted renewal model", "Structural rule", "Case-blind ranking"))
# A representative ladder rather than every scored row: one member of each kernel family, the
# rank-only geographic comparators, and the two case-blind epicentre rules that use no fitted
# kernel - measured cohort inflow (phone) and road travel time (no phone).
#
# BASELINE-TRAVEL-TIME IS ON THE LADDER BECAUSE IT IS THE BEST PHONE-FREE METHOD. It was left
# off, so the dashed "best phone-free" line - computed over the drawn rows only - sat at the
# best phone-free FITTED kernel (24.73 at 1 week, 17.21 at 2 weeks on LINELIST_09102026) while
# road travel time, which fits nothing and uses no phone data, scored 26.12 and 18.56. The
# figure therefore understated what can be done without the phone data. It is listed here so
# it is always drawn, and mv_ladder_rows() additionally adds whichever phone-free method is
# best on the run, so the line is the best over EVERY phone-free row of the ladder by
# construction.
#
# THE FEATURED MODEL IS NEVER TYPED HERE. It is read from the run's model_selection.json
# (read_featured_model(), figure_style.R) by each caller and passed in. And the kernel-family
# rows are taken at the FEATURED MODEL'S OWN COVARIATE SET (its "-med" / "-geo" suffix): the
# ladder exists to isolate the mobility kernel, so every fitted row must share the featured
# model's covariates, or a kernel difference is confounded with a covariate difference. Both
# used to be literals - MV_FEATURED = "Bayes-M14-fill-geo" and every row "-geo" - and when the
# selection moved to Bayes-M14-fill-med, Figure 2 stopped building on a name nobody had updated.
MV_SHOW_FIXED   <- c("Adjacency-B7", "Distance-B1", "Gravity-B4", "Baseline-flowminder-inflow",
                     "Baseline-travel-time")
MV_SHOW_KERNELS <- c("M6a", "M4", "M4c", "M5", "M8-fill", "M13-fill", "M13c-fill", "M10-fill",
                     "M17-fill", "M16-fill", "M14-fill")

#' Split a featured Bayesian label into its kernel and its covariate set:
#' "Bayes-M14-fill-med" -> kernel "M14-fill", covset "med". A time-varying FORM
#' ("Bayes-M14-fill-geo-tvfrailty") is parsed from its stage-1 specification (`base`,
#' "Bayes-M14-fill-geo"): the form is not a kernel or covariate token.
mv_featured_parts <- function(featured) {
  featured <- as.character(featured)
  # The form suffixes are read from 00_config.R, never guessed: a partial fallback list would
  # leave an unlisted form's suffix on the label, where the parse below would take it for the
  # covariate set ("Bayes-M14-fill-geo-tvfrailty" -> kernel "M14-fill-geo", covset "tvfrailty").
  sfx <- get0("TV_FORM_SUFFIXES", ifnotfound = NULL)
  if (!length(sfx))
    stop("[mobval] TV_FORM_SUFFIXES (00_config.R) is not loaded, so a time-varying form's ",
         "label cannot be told from a covariate set.", call. = FALSE)
  base <- if (length(featured) == 1L && !is.na(featured))
    sub(sprintf("-(%s)$", paste(sfx, collapse = "|")), "", featured) else featured
  if (length(featured) != 1L || is.na(featured) || !grepl("^Bayes-.+-[A-Za-z0-9]+$", base))
    stop("[mobval] the featured model '", featured, "' is not a Bayes-<kernel>-<covariates> ",
         "label, so the ladder's like-for-like rows cannot be derived from it.", call. = FALSE)
  list(kernel = sub("^Bayes-(.+)-[A-Za-z0-9]+$", "\\1", base),
       covset = sub("^Bayes-.+-([A-Za-z0-9]+)$", "\\1", base), base = base)
}

#' The ladder's row set for a featured model: the fixed comparators, every kernel family at the
#' featured model's covariate set, and the featured model itself (in case its kernel is not one
#' of the families).
mv_show <- function(featured) {
  cs <- mv_featured_parts(featured)$covset
  unique(c(MV_SHOW_FIXED, paste0("Bayes-", MV_SHOW_KERNELS, "-", cs), featured))
}

#' The two metrics, declared in one place so a panel cannot pair one metric's points with
#' another's interval. `best` is the direction of "better", which sets both the row order and
#' which phone-free method the dashed ceiling line is drawn at.
#'
#' `ref` names the ladder column holding the metric's value under a random ordering, drawn as
#' the dotted rule. For the per-round AUC-PR skill that is `chance` (one value per horizon, the
#' same for every method: it depends only on the rounds and their outcomes); 1 is NOT it — a
#' skill of 1 is a constant forecast, and a random ordering scores well above it. For the mean
#' rank it is ~(number of candidate zones + 1)/2 - in the hundreds here, an order of magnitude
#' above the featured model's mean rank - so drawing it would compress every point in the figure
#' into a narrow strip at the left of the panel. It is carried in the ladder CSV as `null_rank`
#' instead, for the manuscript caption to quote, and `ref` is NA.
#'
#' `lab_h` is the same axis title with the horizon folded into it, for the single-horizon
#' call (`strip = FALSE`) that Figure 2 panel C makes: with one facet the strip would be a
#' lone centred bold heading over the whole panel, which reads as the plot title the house
#' style deliberately blanks. The horizon has to be stated somewhere, so it is stated on the
#' axis instead. `%s` takes the panel label ("2 weeks ahead").
MV_METRIC <- list(
  skill = list(x = "skill", lo = "lo", hi = "hi", best = "max", ref = "chance",
               lab   = "AUC-PR skill per round, all rounds (dotted = random ordering)",
               lab_h = "AUC-PR skill per round at %s, all rounds (dotted = random ordering)"),
  rank  = list(x = "mean_rank", lo = "rank_lo", hi = "rank_hi", best = "min",
               ref = NA_character_,
               lab   = "Mean rank of an invaded zone, all rounds (1 = top of the watch-list; lower is better)",
               lab_h = "Mean rank of an invaded zone at %s, all rounds (1 = top; lower is better)"))

#' The rows of one ladder: the shown methods, labelled, classified and ORDERED for `metric`.
#' Built once per metric and reused by every panel that draws it, so the single-metric figure
#' and the combined figure cannot end up with different orderings.
#' `featured` is the run's featured model (read_featured_model()); it sets the row set
#' (mv_show()) and which row is tagged "(featured)".
mobval_ladder_frame <- function(lad, metric, featured) {
  spec <- MV_METRIC[[metric]]
  if (is.null(spec)) stop("[mobval] unknown metric '", metric, "'.", call. = FALSE)
  show <- mv_ladder_rows(lad, featured)
  # The featured model must be ON the ladder it is tagged in: the "(featured)" suffix is the
  # only thing that marks it, so a ladder without it would argue about a model it does not draw.
  stopifnot(featured %in% show)
  absent <- setdiff(show, lad$method)
  if (length(absent))
    stop("[mobval] the ladder figure asks for methods that were not scored: ",
         paste(absent, collapse = ", "), " (kernel rows are taken at the featured model's ",
         "covariate set, '", mv_featured_parts(featured)$covset, "').", call. = FALSE)
  need <- c(spec$x, spec$lo, spec$hi, if (!is.na(spec$ref)) spec$ref)
  if (!all(need %in% names(lad)))
    stop("[mobval] the ladder table carries no ", metric, " interval or reference (missing ",
         paste(setdiff(need, names(lad)), collapse = ", "),
         "); re-run mobval_compute() rather than drawing from an older CSV.", call. = FALSE)
  # ONE NAME PER METHOD ACROSS FIGURE 2. Panel A names the featured model FEATURED_LBL and the
  # comparators BASELINE_LBL (figure_style.R); this ladder used to carry its own wording for
  # the same objects ("Cohort inflow, 3 epicentre rows only (case-blind, static)" for panel A's
  # "Flowminder cohort inflow from epicentre", "Cases weighted by inverse travel time" for
  # "Neighbours' cases / travel time"), so a reader matching rows to curves had to guess. The
  # comparators panel A publishes now take panel A's names; whether a row is case-blind is
  # still on the page, carried by its shape. The featured row keeps its kernel description
  # after panel A's name, because which kernel it is matters on a ladder of kernels.
  .bl <- intersect(fs_baseline_methods(), names(BASELINE_LBL))
  A <- lad %>% filter(method %in% show) %>%
    mutate(x = .data[[spec$x]], xlo = .data[[spec$lo]], xhi = .data[[spec$hi]],
           mv_class = factor(mv_class, levels = MV_CLS),
           mv_kind  = factor(unname(MV_KIND_MAP[mv_kind]), levels = MV_KIND),
           panel = sprintf("%s week%s ahead", horizon, ifelse(horizon == 1, "", "s")),
           mv_label = ifelse(method %in% .bl, unname(BASELINE_LBL[method]), mv_label),
           mv_label = ifelse(method == featured,
                             paste0(FEATURED_LBL, " (featured): ", mv_label), mv_label),
           # Wrapped here rather than shortened: only the featured row is long enough to wrap,
           # and both halves of it are doing work - panel A's name, and which kernel it is.
           # AT MOST TWO LINES: the row pitch holds two, and a third prints over the neighbouring
           # rows. A time-varying form's name ("..., gamma frailty") does not fit two lines of 42
           # characters, so a label is wrapped at the narrowest width from 42 up that gives two.
           # The featured row breaks after "(featured):" when the model's own name then fits one
           # line of up to 52 characters (it reads as name over description); otherwise, and for
           # every other row, a greedy wrap at the narrowest such width.
           mv_label = unname(vapply(mv_label, function(z) {
             pre <- paste0(FEATURED_LBL, " (featured): ")
             if (startsWith(z, pre) && nchar(z) - nchar(pre) <= 52L)
               return(paste0(sub(" $", "", pre), "\n", substring(z, nchar(pre) + 1L)))
             w <- 42L
             while (length(strwrap(z, width = w)) > 2L) w <- w + 1L
             paste(strwrap(z, width = w), collapse = "\n")
           }, character(1))))
  if (anyNA(A$mv_class) || anyNA(A$mv_kind))
    stop("[mobval] a shown method carries a class or kind that is not in the legend.",
         call. = FALSE)
  if (any(!is.finite(A$x) | !is.finite(A$xlo) | !is.finite(A$xhi)))
    stop("[mobval] a non-finite ", metric, " point or interval reached the figure.",
         call. = FALSE)
  # THE DASHED LINE IS THE BEST OF EVERY PHONE-FREE ROW OF THE LADDER, not of the drawn rows.
  # mobval_panel() computes it from the rows it is given, so the drawn set must contain the
  # ladder's best phone-free method at every horizon; mv_ladder_rows() adds it, and this says
  # so out loud rather than trusting the construction (the line once sat 1.4 below the true
  # best because that method was simply not drawn).
  best_fun <- if (identical(spec$best, "max")) max else min
  for (hh in unique(A$horizon)) {
    all_pf <- lad[[spec$x]][lad$horizon == hh & lad$phone %in% FALSE & is.finite(lad[[spec$x]])]
    drw_pf <- A$x[A$horizon == hh & A$phone %in% FALSE]
    if (length(all_pf) &&
        (!length(drw_pf) || !isTRUE(all.equal(best_fun(all_pf), best_fun(drw_pf)))))
      stop(sprintf(paste0("[mobval] at h=%d the best phone-free %s on the ladder is %.3f but the ",
                          "best DRAWN phone-free row is %s; the dashed line would understate ",
                          "what can be done without the phone data."),
                   hh, metric, best_fun(all_pf),
                   if (length(drw_pf)) sprintf("%.3f", best_fun(drw_pf)) else "absent"),
           call. = FALSE)
  }
  # Order from the longest horizon, but only after checking every row has a level there: a
  # label present at one horizon and missing at the other would become NA and vanish from the
  # panel without any error.
  ord_h <- if (2L %in% A$horizon) 2L else max(A$horizon)
  o <- A %>% filter(horizon == ord_h)
  lev <- if (identical(spec$best, "max")) o %>% arrange(x) %>% pull(mv_label)
         else                             o %>% arrange(desc(x)) %>% pull(mv_label)
  if (anyDuplicated(lev) || !all(A$mv_label %in% lev))
    stop("[mobval] the ladder's rows are not one-to-one across horizons; a method appears at ",
         "one horizon but not the other, or two carry the same label.", call. = FALSE)
  A$mv_label <- factor(A$mv_label, levels = lev)
  A
}

#' One panel: a metric, one or both horizons, row labels on or off.
#'
#' `show_y` is off for the right-hand column of the combined figure, where patchwork has
#' already aligned the rows against the labelled panel beside it and a second copy of the row
#' labels would only take width away from the data.
#'
#' `base` is the theme_pub() base size. Every hard-coded size in the panel — the row labels,
#' the points, the intervals and the "best phone-free" annotation — is scaled by base / 9.5
#' with it, so the main-text call (Figure 2 panel C, which runs at the larger Figure-2 base
#' font) is the SAME panel at a larger size and not a panel whose annotations have quietly
#' shrunk relative to its text. At the default 9.5 the scale factor is exactly 1, so the three
#' appendix figures are unchanged.
#'
#' `strip` = FALSE drops the facet and folds the horizon into the axis title instead (see
#' `lab_h`). It is only meaningful for a SINGLE horizon, and refuses otherwise: dropping the
#' facet with two horizons in `a` would silently overplot 1-week and 2-week points on one set
#' of rows, which is exactly the kind of figure that looks right and is not.
#'
#' `n_breaks` is the requested x-axis break count. The default 4 suits the ~5-inch facets of
#' the appendix figures; a panel given the full width of Figure 2 can carry more.
mobval_panel <- function(A, metric, horizons = sort(unique(A$horizon)), show_y = TRUE,
                         base = 9.5, strip = TRUE, n_breaks = 4) {
  spec <- MV_METRIC[[metric]]
  if (is.null(spec)) stop("[mobval] unknown metric '", metric, "'.", call. = FALSE)
  a <- A %>% filter(horizon %in% horizons)
  if (!nrow(a))
    stop("[mobval] no rows at horizon(s) ", paste(horizons, collapse = ", "), ".", call. = FALSE)
  if (!strip && dplyr::n_distinct(a$panel) != 1L)
    stop("[mobval] strip = FALSE removes the facet, so it draws ONE horizon; ",
         dplyr::n_distinct(a$panel), " were selected (",
         paste(sort(unique(a$panel)), collapse = ", "),
         "). Pass a single value in `horizons`.", call. = FALSE)
  sc <- base / 9.5                     # every hard-coded size below scales with the base font
  # The dashed line: the BEST method here that uses no mobile-phone data at any point, drawn
  # per horizon. Everything to its right (skill) or its left (rank) is what the phone data buys.
  best_row <- if (identical(spec$best, "max")) dplyr::slice_max else dplyr::slice_min
  ceil <- a %>% filter(!phone) %>% group_by(panel) %>%
    best_row(x, n = 1, with_ties = FALSE) %>% ungroup() %>%
    dplyr::select(panel, x, mv_class)
  if (nrow(ceil) != dplyr::n_distinct(a$panel))
    stop("[mobval] a horizon has no phone-free method; the reference line would be missing ",
         "from one panel only.", call. = FALSE)
  # The dashed rule is drawn in the class colour of the method it marks (see MV_NOSKILL_COL
  # above), PER HORIZON. It used to take one colour for the whole figure and stopped when the
  # horizons disagreed; they can - with the best phone-free model taken "best overall", the
  # rank ladder's best phone-free method was radiation ("Geography + population") at 1 week and
  # travel-time decay + covariates ("Geography only") at 2 weeks on LINELIST_23092026. Each
  # facet's rule therefore carries its OWN method's class colour, mapped through the colour
  # scale on the facet's own row of `ceil` (never a loose colour vector, which a facet split
  # can misalign). Both conditions that make the colour truthful are still checked per facet.
  ceil_cls <- as.character(ceil$mv_class)
  ceil_col <- unname(MV_PAL[ceil_cls])
  if (anyNA(ceil_col))
    stop("[mobval] a best phone-free method carries class '",
         paste(ceil_cls[is.na(ceil_col)], collapse = "', '"), "', which is not in the palette.",
         call. = FALSE)
  if (any(ceil_col == MV_NOSKILL_COL))
    stop("[mobval] a best phone-free method is in the '",
         paste(ceil_cls[ceil_col == MV_NOSKILL_COL], collapse = "', '"), "' class, whose colour ",
         "is the same grey as the random-ordering rule - the two reference lines would look alike ",
         "again. Give that class a different colour, or the random-ordering rule one.", call. = FALSE)
  nlev <- nlevels(a$mv_label)

  p <- ggplot(a, aes(x, mv_label, colour = mv_class))
  if (!is.na(spec$ref)) {
    # One random-ordering value per horizon (it does not depend on the method), drawn per facet.
    ref <- a %>% group_by(panel) %>%
      summarise(x = unique(round(.data[[spec$ref]], 12)), .groups = "drop")
    if (anyDuplicated(ref$panel) || any(!is.finite(ref$x)))
      stop("[mobval] the ladder's '", spec$ref, "' is not one finite value per horizon; it is ",
           "a property of the rounds, not of the method.", call. = FALSE)
    p <- p + geom_vline(data = ref, aes(xintercept = x), colour = MV_NOSKILL_COL,
                        linetype = "dotted", linewidth = 0.4, inherit.aes = FALSE)
  }
  p +
    geom_vline(data = ceil, aes(xintercept = x, colour = mv_class), linetype = "22",
               linewidth = 0.5, show.legend = FALSE) +
    # The line is named in the plot, not in a caption: the figures carry no caption, and a
    # reference line a reader cannot name is worse than no line at all. It sits in the strip of
    # empty space the y expansion opens above the top row, so it can never cross a point.
    geom_text(data = ceil, aes(x = x, y = nlev + 0.95, label = "best phone-free"),
              hjust = -0.06, vjust = 0.4, size = 2.15 * sc, colour = MUTED,
              fontface = "italic", inherit.aes = FALSE) +
    geom_errorbar(aes(xmin = xlo, xmax = xhi), width = 0, linewidth = 0.45 * sc, alpha = 0.8,
                  orientation = "y") +
    geom_point(aes(shape = mv_kind), size = 2.1 * sc, fill = "white", stroke = 0.7 * sc) +
    (if (strip) facet_wrap(~ panel, nrow = 1, scales = "free_x") else NULL) +
    scale_colour_manual(values = MV_PAL, name = "The mobility term is built from",
                        drop = FALSE) +
    scale_shape_manual(values = MV_SHP, name = "What the method is", drop = FALSE) +
    # Anchored at zero. Both metrics are ratio-scale, and an axis that started at the smallest
    # interval bound would exaggerate every difference in the figure.
    scale_x_continuous(limits = c(0, NA), breaks = scales::pretty_breaks(n_breaks),
                       expand = expansion(mult = c(0.02, 0.06))) +
    scale_y_discrete(expand = expansion(add = c(0.7, 1.5))) +
    labs(x = if (strip) spec$lab else sprintf(spec$lab_h, unique(a$panel)), y = NULL) +
    # One row each, stacked: the class keys (one per MV_CLS entry) and the three shape keys
    # both fit the panel width, and a single row per legend keeps each axis of the encoding
    # readable as one sentence.
    guides(colour = guide_legend(nrow = 1, title.position = "top", order = 1),
           shape  = guide_legend(nrow = 1, title.position = "top", order = 2)) +
    theme_pub(base) +
    # The two legends STACK. Side by side (ggplot's default) the shape legend's last key runs
    # off the right edge of a 10-inch figure, and a key that is cut off is worse than a
    # taller legend.
    theme(legend.position = "top", legend.box = "vertical", legend.box.just = "left",
          legend.spacing.y = unit(1, "pt"), legend.margin = margin(0, 0, 2, 0),
          axis.text.y = if (show_y) element_text(size = 7.4 * sc, colour = INK, hjust = 1,
                                                 lineheight = 0.95) else element_blank(),
          # Kept, faintly: the rows are read against two panels, and the eye has to travel from
          # a label on the left to a point on the far right of the second panel.
          panel.grid.major.y = element_line(colour = GRID, linewidth = 0.25),
          panel.spacing.x = unit(11, "pt"))
}


# =============================================================================
# THE FOUR-ROW SUMMARY OF THE SAME LADDER  (Figure2_simple, panel C)
# =============================================================================
# WHAT IT IS. The full ladder above answers "what does the mobility data add?" by showing
# every rung. This answers the same question with the four rungs that carry the argument,
# for a reader who is not going to read every row:
#
#   1  the FEATURED MODEL - the forecasting model panels A and B draw, which uses human
#      mobility data. It is the featured model BY DEFINITION, not the argmax of a search: the
#      label used to say "Best model using human mobility data", which is not a claim the
#      selection makes (the featured model is chosen by the composite of three metrics, not by
#      AUC-PR skill alone, so another candidate - or a generation-time sensitivity arm, which
#      is not eligible - can score higher on this one metric).
#   2  human mobility data with NO MODEL AT ALL — the case-blind ranking straight off the
#      MEASURED movement out of the epicentre (mv_class "Measured phone movement", so the
#      phone-calibrated gravity ranking, which is a fitted model of flows, can never be
#      presented as "data alone")
#   3  the best FITTED MODEL that uses no human mobility data
#   4  the best RULE that uses no human mobility data and is not a renewal model — the
#      structural rules that use the case data to locate the front AND the case-blind
#      geographic ranking (road travel time from the epicentre)
#
# ROW 4 USED TO SEARCH THE STRUCTURAL RULES ONLY, and was labelled "Best fixed rule, nothing
# fitted". Both were wrong. Searching only mv_kind == "Structural rule" left out road travel
# time, which is phone-free, case-blind and fits nothing, and which outscores every structural
# rule (18.56 against Gravity-B4's 12.70 at 2 weeks on LINELIST_09102026); and the method it
# then reported, Gravity-B4, DOES fit something - a scale constant, by maximum likelihood
# (forecast_B4(), 05_baseline_models.R), even though a scalar cannot move its ranking. The row
# now searches every phone-free non-renewal method and its label claims only what is true of
# all of them: no mobility data, no renewal model.
#
# ROW 2 IS NOT OPTIONAL AND IS THE REASON THIS IS FOUR ROWS RATHER THAN THREE. Without it the
# panel reads "model with mobility data against model without", and a reader would reasonably
# attribute the whole gap to the modelling, when part or all of it may be in the DATA rather
# than in what the model does with it. Dropping row 2 could make the figure overstate the
# model's contribution, so it stays (the caption written by make_manuscript_figures.R says
# which way the numbers fall on the run).
#
# THE SAME LADDER CSV, THE SAME SHARED SUPPORT. Every row here is a row of
# mobility_value_ladder.csv with its published interval; nothing is recomputed. The four are
# drawn from DISJOINT groups of (mv_kind x phone), so one method cannot occupy two rows.
#
# "USES HUMAN MOBILITY DATA" IS THE LADDER'S OWN `phone` COLUMN, which is TRUE for a method
# that touches the mobile-network data at any point — including one whose kernel only had its
# distance-decay exponents CALIBRATED on it (Baseline-gravity, Bayes-M4*). That is why row 3's
# "best model without human mobility data" is never a gravity-kernel model: the phone-free
# fitted models are the radiation and travel-time-decay ones.

#' The four rows, as selection rules over the ladder. `mv_kind` is a LIST of the kinds the row
#' may draw from; NA in any rule column means "do not filter on it". Row 1 ("model_mob") is
#' not a search: it is the featured model (see mobval_summary_frame()), and its kind/phone
#' entries are what that model is REQUIRED to be for the row label to be true.
MV_SUMMARY_SPEC <- tibble::tribble(
  ~row_key,      ~row_label,                                ~mv_kind,
                 ~phone, ~mv_class,
  "model_mob",   FEATURED_LBL_PLAIN,                        list("Fitted renewal model"),
                 TRUE,   NA_character_,
  "data_only",   "Human mobility data alone, no model",     list("Case-blind ranking"),
                 TRUE,   "Measured phone movement",
  "model_nomob", "Best model without human mobility data",  list("Fitted renewal model"),
                 FALSE,  NA_character_,
  "rule",        "Best rule without human mobility data",
                 list(c("Structural rule", "Case-blind ranking")),
                 FALSE,  NA_character_)

# Two colours, one meaning: does this row touch the human mobility data? The row LABELS
# already say whether anything was fitted, so shape carries nothing here and is not used —
# one encoding is the point of the simplified panel.
#
# DELIBERATELY NOT FEATURED_COL. Panel A draws the forecasting model — the one that uses the
# mobility data — in #0072B2, so giving that blue to "no human mobility data" here would put
# one colour on opposite meanings in two panels of one figure. Vermillion is the tone MV_PAL
# already gives to measured movement (and panel A to the measured cohort inflow); grey is the
# neutral tone for the absence of the data, and no MV_PAL class takes it.
MV_USES_MOB <- c("Uses human mobility data", "No human mobility data")
MV_MOB_PAL  <- setNames(c("#D55E00", "#7A7A7A"), MV_USES_MOB)

#' One summary row's pick over a one-horizon slice `d` of the ladder: the best method of the
#' row's kind(s), phone flag and class over EVERY method in `d`, not only the rows the full
#' ladder draws. An NA kind, phone flag or class means "any". Zero rows when nothing
#' qualifies; the caller decides whether that is fatal.
.mobval_pick <- function(d, s, spec) {
  best_row <- if (identical(spec$best, "max")) dplyr::slice_max else dplyr::slice_min
  cand  <- d
  kinds <- unlist(s$mv_kind)
  if (length(kinds) && !all(is.na(kinds))) cand <- cand %>% filter(mv_kind %in% kinds)
  if (!is.na(s$phone))   cand <- cand %>% filter(phone %in% s$phone)
  cls <- if ("mv_class" %in% names(s)) s$mv_class else NA_character_
  if (!is.na(cls))       cand <- cand %>% filter(mv_class %in% cls)
  cand <- cand %>% filter(is.finite(.data[[spec$x]]))
  if (!nrow(cand)) return(cand)
  cand %>% best_row(.data[[spec$x]], n = 1, with_ties = FALSE)
}

#' The rows the full ladder draws (Figure 2 panel C, FigureS24-S26): mv_show(featured) - the
#' fixed comparators and every kernel family at the featured model's covariate set - PLUS, at
#' every horizon the ladder carries,
#'   * the best model WITHOUT human mobility data over the WHOLE ladder, any covariate set: the
#'     simplified panel's row 3, which is chosen "best overall" (decided 2026-10-03). The
#'     like-for-like rows alone leave that panel naming a method the ladder never draws
#'     whenever the best phone-free model carries a different covariate set from the featured
#'     one - on LINELIST_23092026, travel-time decay WITH covariates against a covariate-free
#'     featured model;
#'   * the best method of ANY KIND without human mobility data, under BOTH metrics: the method
#'     the dashed "best phone-free" line is drawn at in the skill ladder and in the rank ladder.
#'     mobval_panel() takes that line over the rows it is given, so a phone-free method missing
#'     from the drawn set silently lowers the line (mobval_ladder_frame() asserts it cannot).
#' Every extra row is found by a rule and never typed, and the union is taken over both metrics
#' so the skill and rank ladders draw the SAME rows.
mv_ladder_rows <- function(lad, featured) {
  s <- MV_SUMMARY_SPEC[MV_SUMMARY_SPEC$row_key == "model_nomob", ]
  stopifnot(nrow(s) == 1L)
  pf_any <- tibble::tibble(mv_kind = list(NA_character_), phone = FALSE,
                           mv_class = NA_character_)
  extra <- unlist(lapply(sort(unique(lad$horizon)), function(hh) {
    dh <- lad[lad$horizon == hh, , drop = FALSE]
    c(.mobval_pick(dh, s, MV_METRIC$skill)$method,
      unlist(lapply(MV_METRIC, function(sp)
        if (sp$x %in% names(dh)) .mobval_pick(dh, pf_any, sp)$method else character(0))))
  }))
  unique(c(mv_show(featured), extra))
}

#' Pick the four summary rows out of a ladder, at one horizon.
#'
#' @param lad    the mobility_value_ladder.csv frame.
#' @param h      forecast horizon in weeks.
#' @param metric "skill" (the only one the simplified panel draws; see MV_METRIC).
#' @param featured the run's featured model (read_featured_model()); row 1 IS it.
#' @return a tibble of four rows: row_key, row_label (a factor whose levels run worst to best,
#'   so the best row is drawn at the TOP), method, mv_label, x/xlo/xhi, phone, n_cells,
#'   uses_mob; stamped with attr "mv_metric".
mobval_summary_frame <- function(lad, h = 2L, metric = "skill", featured) {
  show <- mv_ladder_rows(lad, featured)
  spec <- MV_METRIC[[metric]]
  if (is.null(spec)) stop("[mobval] unknown metric '", metric, "'.", call. = FALSE)
  need <- c("method", "horizon", "mv_kind", "mv_class", "phone", "mv_label", "n_cells",
            spec$x, spec$lo, spec$hi)
  miss <- setdiff(need, names(lad))
  if (length(miss))
    stop("[mobval] the ladder table is missing column(s) ", paste(miss, collapse = ", "),
         "; re-run 50_mobility_value.R.", call. = FALSE)
  if (!is.logical(lad$phone))
    stop("[mobval] the ladder's `phone` column is ", class(lad$phone)[1], ", not logical; a ",
         "string \"FALSE\" would test TRUE and put phone-free methods in the mobility rows.",
         call. = FALSE)
  d <- lad %>% filter(horizon == as.integer(h))
  if (!nrow(d))
    stop("[mobval] the ladder carries no rows at h=", h, ".", call. = FALSE)

  .rule_txt <- function(s) paste0(
    "kind ", paste(unlist(s$mv_kind), collapse = " or "),
    if (is.na(s$phone)) "" else paste0(", phone = ", s$phone),
    if (is.na(s$mv_class)) "" else paste0(", class ", s$mv_class))
  picked <- lapply(seq_len(nrow(MV_SUMMARY_SPEC)), function(i) {
    s    <- MV_SUMMARY_SPEC[i, ]
    # ROW 1 IS THE FEATURED MODEL, NOT A SEARCH. Its kind / phone entries are what the featured
    # model must BE for the row's label ("uses human mobility data", a fitted model) to be true,
    # so they are checked against it rather than searched over.
    cand <- if (identical(s$row_key, "model_mob")) {
      f <- d %>% filter(method == featured)
      if (nrow(f) != 1L)
        stop("[mobval] the ladder holds ", nrow(f), " rows for the featured model '", featured,
             "' at h=", h, "; the simplified panel's first row is that model.", call. = FALSE)
      if (!(f$mv_kind %in% unlist(s$mv_kind)) || !(f$phone %in% s$phone))
        stop("[mobval] the featured model '", featured, "' is classed '", f$mv_kind,
             "', phone = ", f$phone, ", but row 1 is labelled '", s$row_label, "' (",
             .rule_txt(s), "); the label would be false.", call. = FALSE)
      f
    } else .mobval_pick(d, s, spec)
    if (!nrow(cand))
      stop("[mobval] no method at h=", h, " matches the summary row '", s$row_label,
           "' (", .rule_txt(s), "). The simplified ",
           "panel C cannot be drawn with a row missing — it would make a different claim ",
           "under the same design.", call. = FALSE)
    cand %>%
      transmute(row_key = s$row_key, row_label = s$row_label, method, mv_label,
                x = .data[[spec$x]], xlo = .data[[spec$lo]], xhi = .data[[spec$hi]],
                ref = if (!is.na(spec$ref)) .data[[spec$ref]] else NA_real_,
                phone, n_cells)
  })
  S <- bind_rows(picked)

  # --- the checks that keep the four rows honest ------------------------------
  # (a) DISJOINT BY CONSTRUCTION, asserted anyway: (mv_kind x phone) partitions the ladder, so
  #     one method cannot satisfy two of the rules above. If the classification table ever
  #     changed so that it could, the panel would draw one method twice under two labels.
  if (anyDuplicated(S$method))
    stop("[mobval] the same method (", paste(S$method[duplicated(S$method)], collapse = ", "),
         ") was picked for more than one summary row; the selection rules have stopped being ",
         "disjoint.", call. = FALSE)
  # (b) ONE SUPPORT. Every row must have been scored on the same (origin x zone) cells, or the
  #     four numbers are not comparable and the panel is a ranking of different experiments.
  if (dplyr::n_distinct(S$n_cells) != 1L)
    stop("[mobval] the four summary rows were scored on different supports (",
         paste(sort(unique(S$n_cells)), collapse = " / "), " cells); they are not comparable.",
         call. = FALSE)
  # (c) FINITE INTERVALS.
  if (any(!is.finite(c(S$x, S$xlo, S$xhi))))
    stop("[mobval] a non-finite point or interval reached the simplified panel.", call. = FALSE)
  # (d) EVERY ROW MUST ALSO BE ON THE APPENDIX LADDER. FigureS24 draws mv_show(); a main-text row
  #     naming a method that ladder omits would send the reader to an appendix that does not
  #     contain it. The argmaxes here are taken over the WHOLE ladder, not over mv_show(), so
  #     that "best" means best — this check is what makes the two consistent, by failing
  #     rather than by quietly restricting the search.
  off <- setdiff(S$method, show)
  if (length(off))
    stop("[mobval] the simplified panel calls ", paste(off, collapse = ", "),
         " best at h=", h, ", but the appendix ladder (mv_show()) omits it, so it is never ",
         "drawn there. Extend MV_SHOW_FIXED / MV_SHOW_KERNELS, or the two figures will name ",
         "different methods.", call. = FALSE)
  # (d2) ROWS 3 AND 4 TOGETHER COVER EVERY PHONE-FREE METHOD. FigureS24's dashed "best
  #      phone-free" rule is the best method of ANY KIND that uses no mobility data; rows 3 and
  #      4 split that set into fitted renewal models and everything else, so the better of the
  #      two IS that method. Before row 4 took in the case-blind rules this could fail - road
  #      travel time was the best phone-free method and neither row could pick it - so it is
  #      asserted, and which row holds it is logged for whoever writes the caption.
  .pf <- .mobval_pick(d, tibble::tibble(mv_kind = list(NA_character_), phone = FALSE,
                                        mv_class = NA_character_), spec)
  .r34 <- S[S$row_key %in% c("model_nomob", "rule"), ]
  if (nrow(.pf) == 1L && !(as.character(.pf$method) %in% as.character(.r34$method)))
    stop(sprintf(paste0("[mobval] h=%d: the best phone-free method on the ladder is %s (%.2f), ",
                        "but neither phone-free summary row picked it (%s); rows 3 and 4 no ",
                        "longer cover every phone-free method."),
                 h, .pf$method, .pf[[spec$x]], paste(.r34$method, collapse = ", ")),
         call. = FALSE)
  if (nrow(.pf) == 1L)
    message(sprintf("[mobval] h=%d: the best phone-free method is %s (%.2f), row '%s'.",
                    h, .pf$method, .pf[[spec$x]],
                    .r34$row_label[match(.pf$method, .r34$method)]))
  # (e) ROW 1 IS THE MODEL PANELS A AND B DRAW, by construction above. Asserted anyway: it is
  #     the one row whose identity the rest of the figure depends on.
  stopifnot(identical(as.character(S$method[S$row_key == "model_mob"]), as.character(featured)))

  # (f) NO NA IN `phone`. The colour is the ONLY encoding on this panel, and an NA would make
  #     ifelse() return NA, give the row an NA factor level, and let ggplot DROP the point with
  #     nothing but a warning — a four-row panel silently published with three rows. Every row
  #     is now filtered (or checked) on `phone` with %in%, which excludes NA; this stays as the
  #     guard for any rule that is later given phone = NA ("any").
  if (anyNA(S$phone))
    stop("[mobval] the ladder gives method(s) ", paste(S$method[is.na(S$phone)], collapse = ", "),
         " an NA `phone` flag, so the simplified panel cannot say whether they use human ",
         "mobility data. Fix mobval_classes() in 50_mobility_value.R.", call. = FALSE)
  S$uses_mob <- factor(ifelse(S$phone, MV_USES_MOB[1], MV_USES_MOB[2]), levels = MV_USES_MOB)
  # Best at the TOP: ggplot draws the first level of a discrete y at the bottom, so the levels
  # run worst-to-best.
  lev <- S$row_label[order(if (identical(spec$best, "max")) S$x else -S$x)]
  S$row_label <- factor(S$row_label, levels = lev)
  # LAST, after every column assignment. The panel draws one metric's semantics (a random-ordering
  # reference at the ladder's `chance`, an axis title in multiples of the base rate), so it refuses a frame
  # built from a different metric rather than mislabel it — but a custom attribute set before
  # a `$<-` is at the mercy of whatever that method does with attributes, so it is stamped
  # once nothing else will touch the frame.
  attr(S, "mv_metric") <- metric
  S
}

#' The simplified panel C: four rows, one colour encoding, the value printed at each point.
#'
#' @param S    mobval_summary_frame() output.
#' @param h    horizon, for the axis title.
#' @param base theme_pub() base size; the hard-coded sizes scale with it as in mobval_panel().
#' @param n_breaks requested x-axis break count.
mobval_summary_panel <- function(S, h = 2L, base = 9.5, n_breaks = 7) {
  # THIS PANEL IS SPECIFIC TO THE SKILL METRIC: the random-ordering reference below is the skill
  # metric's own (the ladder's `chance`), and the axis title reads the value as a multiple of
  # the base rate. Neither is true of the mean-rank metric (whose random-ordering value is
  # ~(number of candidate zones + 1)/2 and whose better direction is DOWN), so a rank frame is
  # refused rather than drawn under a skill axis.
  if (!identical(attr(S, "mv_metric"), "skill"))
    stop("[mobval] mobval_summary_panel() draws the 'skill' metric; it was given a frame ",
         "stamped '", attr(S, "mv_metric"), "'.", call. = FALSE)
  .ref <- unique(round(S$ref, 12))
  if (length(.ref) != 1L || !is.finite(.ref))
    stop("[mobval] the summary rows do not carry one finite random-ordering reference.", call. = FALSE)
  sc <- base / 9.5
  hl <- sprintf("%d week%s ahead", as.integer(h), if (as.integer(h) == 1L) "" else "s")
  ggplot(S, aes(x, row_label, colour = uses_mob)) +
    # The reference is a property of the AXIS, not of any method: the per-round skill a
    # uniformly random ordering scores on these rounds (not 1, which is a constant forecast).
    geom_vline(xintercept = .ref, colour = FAINT, linetype = "dotted", linewidth = 0.4) +
    geom_errorbar(aes(xmin = xlo, xmax = xhi), width = 0, linewidth = 0.6 * sc, alpha = 0.85,
                  orientation = "y") +
    geom_point(size = 2.6 * sc) +
    # The number, above its own point. With four rows there is vertical room for it, and a
    # reader should not have to measure a position against an axis to get the headline value.
    geom_text(aes(label = sprintf("%.1f×", x)), vjust = -1.05, size = 2.9 * sc,
              show.legend = FALSE) +
    scale_colour_manual(values = MV_MOB_PAL, name = NULL, drop = FALSE) +
    # Anchored at zero: the metric is ratio-scale, and an axis starting at the lowest interval
    # bound would exaggerate every difference in the panel. The right-hand expansion leaves
    # room for the value printed above the widest interval's point.
    scale_x_continuous(limits = c(0, NA), breaks = scales::pretty_breaks(n_breaks),
                       expand = expansion(mult = c(0.02, 0.08))) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
    labs(x = sprintf(paste0("Per-round AUC-PR skill, %s, all rounds (multiple of the base rate; ",
                            "dotted: random ordering, %.1f)"), hl, .ref),
         y = NULL) +
    guides(colour = guide_legend(nrow = 1)) +
    theme_pub(base) +
    theme(legend.position = "top",
          axis.text.y = element_text(size = 8.4 * sc, colour = INK, hjust = 1),
          panel.grid.major.y = element_line(colour = GRID, linewidth = 0.25))
}
