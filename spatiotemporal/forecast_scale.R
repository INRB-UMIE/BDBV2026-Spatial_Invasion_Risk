# =============================================================================
# forecast_scale.R — which probability scale a figure/table pass is drawn on
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# TWO SCALES ARE IN CIRCULATION, produced by DIFFERENT estimators. Keeping them
# apart is the whole point of this file, because a caption that names the wrong
# one is a false statement about how the forecast was validated.
#
#   "recalibrated"  — the PRIMARY artifact set (default).
#       * retrospective, cross-validated rows use `p_recal`: the PREQUENTIAL
#         factor from 16b_invasion_recalibration.R, fitted for each fold only on
#         folds whose outcome window had already closed. Leakage-free, and the
#         honest answer to "how would this have performed in real time".
#       * the live current forecast uses the POOLED factor — fitted on all folds,
#         which for a live forecast genuinely is the past. This is the deployed
#         scale, stamped `prob_scale = "recalibrated-deployed"` in the tables.
#       There is NO leave-future-out version of a live forecast: there is no
#       outcome to hold out. The two halves are both calibrated, by different
#       estimators, and captions must say which.
#
#   "raw"           — the uncorrected model output, written to a `raw/` sibling
#       directory with the SAME basenames, so a raw file is the same file on the
#       other scale and nothing downstream learns a second naming convention.
#
# WHAT ACTUALLY DIFFERS. The transform p -> 1-(1-p)^delta is strictly monotone,
# so a ranking taken WITHIN A FOLD is EXACTLY invariant: capture/detection curves,
# top-K precision, rank-of-truth, per-round top-K panels, watch-lists. Every
# ranking in the figure suites groups by fold_id or cutoff first, so all of them
# are invariant — but note the qualifier is doing real work. The PREQUENTIAL
# factor varies by fold, so a ranking POOLED across folds is NOT invariant (on the
# 2026-09-07 frame the pooled h=1 ordering does change), which is the same
# fold-mixing caveat 16b documents for pooled AUC. Do not add an ungrouped rank.
#
# `mu`-based relative risks scale out exactly within a fold. What moves is
# reliability, count calibration, Brier, log score, calibration-in-the-large,
# per-round probability bars, top-zone probability intervals, and the probability
# columns of the risk tables. The magnitude of the shift is the calibration-in-the-large
# reported per model and horizon in invasion_evaluation.csv, and it moves every run — the
# worked example once given here (78.6 -> 52.0 against "41 realised invasions") was frozen to
# a frame whose h=1 event count is now 37, in the one file whose whole purpose is stopping a
# caption from naming the wrong scale. The 0-1 relative-risk index divides transformed
# probabilities, so it is invariant only in the small-probability limit.
#
# Usage (both figure suites): source this file, then
#     PCOL <- fs_apply_lfo_scale(lfo)   # swaps lfo$p_invasion in place
#     rs   <- readr::read_csv(fs_risk_csv(OUT))
#     dir  <- fs_out_dir(FIG_DIR)
# and add fs_stamp() to each figure's caption.
# =============================================================================

FORECAST_SCALE <- local({
  v <- tolower(trimws(Sys.getenv("FORECAST_SCALE", "recalibrated")))
  if (v %in% c("", "recal", "recalibrated", "calibrated")) "recalibrated"
  else if (v %in% c("raw", "uncorrected", "uncalibrated")) "raw"
  else {
    warning(sprintf("[scale] unreadable FORECAST_SCALE='%s'; using 'recalibrated'.", v),
            call. = FALSE)
    "recalibrated"
  }
})

fs_is_raw <- function(scale = FORECAST_SCALE) identical(scale, "raw")

#' Which probability column of an LFO frame this pass should draw on.
#'
#' Falls back to the raw column — loudly — when the recalibrated one is absent or is
#' missing on rows that are actually scored, because a silent fallback would publish a
#' figure whose caption says "recalibrated" over raw probabilities. The rows that matter
#' are the ones every figure filters to: at-risk, finite probability.
fs_lfo_col <- function(lfo, scale = FORECAST_SCALE) {
  if (fs_is_raw(scale)) return("p_invasion")
  if (!"p_recal" %in% names(lfo)) {
    warning("[scale] FORECAST_SCALE='recalibrated' but the LFO frame carries no `p_recal` ",
            "(was INVASION_RECALIBRATE on for the run that wrote it?); drawing on the RAW ",
            "probabilities.", call. = FALSE)
    return("p_invasion")
  }
  scored <- is.finite(lfo$p_invasion)
  if ("was_active_before" %in% names(lfo))
    scored <- scored & !(as.logical(lfo$was_active_before) %in% TRUE)
  # A RANK-ONLY METHOD HAS NO p_recal BY CONSTRUCTION, NOT BY FAILURE. Distance-B1 and
  # Adjacency-B7 emit an ordering, not a probability, so 16b deliberately never fits them a
  # delta. Counting their rows as "missing" made this function return the RAW column for the
  # ENTIRE frame -- 22,112 of 234,288 scored rows were enough to demote every panel in both
  # figure suites to raw probabilities, while fs_stamp() (which keys on the REQUESTED scale,
  # not the column actually used) went on captioning them "recalibrated". Every published
  # panel was mislabelled. Judge the column on the rows that can carry it.
  if ("prob_calibrated" %in% names(lfo))
    scored <- scored & !(lfo$prob_calibrated %in% FALSE)
  bad <- sum(scored & !is.finite(lfo$p_recal))
  if (bad > 0L) {
    warning(sprintf(paste0("[scale] %d of %d scored rows have no recalibrated probability; ",
                           "drawing on the RAW probabilities instead of a mixed scale."),
                    bad, sum(scored)), call. = FALSE)
    return("p_invasion")
  }
  "p_recal"
}

#' Put the selected scale into `p_invasion`, preserving the raw column.
#'
#' Every downstream expression in the figure suites reads `p_invasion`; swapping the
#' column once here is what makes the scale switch total rather than depending on a
#' dozen call sites being edited consistently. The untouched values stay available as
#' `p_invasion_raw` for the panels that show both scales.
#'
#' @return the name of the source column, for logging.
fs_apply_lfo_scale <- function(lfo, scale = FORECAST_SCALE, envir = parent.frame(),
                               name = deparse(substitute(lfo))) {
  # FORCE `name` before `lfo` is touched. Assigning to a formal replaces its promise, after
  # which the lazily-evaluated default substitute(lfo) yields the VALUE rather than the
  # symbol — and deparsing a 300k-row data frame produced a 160 MB log line before this
  # line existed. The guard also keeps an inline expression from becoming the label.
  force(name)
  if (length(name) != 1L || is.na(name) || nchar(name) > 40L) name <- "lfo"
  col <- fs_lfo_col(lfo, scale)
  lfo$p_invasion_raw <- lfo$p_invasion
  if (!identical(col, "p_invasion")) lfo$p_invasion <- lfo[[col]]
  # `mu_forecast` is the cumulative hazard and scales EXACTLY by the same factor, so it is
  # carried onto the selected scale too; a figure mixing a recalibrated probability with a
  # raw hazard would be internally inconsistent.
  if (!identical(col, "p_invasion") && "mu_forecast" %in% names(lfo)) {
    if ("delta_preq" %in% names(lfo)) {
      lfo$mu_forecast_raw <- lfo$mu_forecast
      d <- lfo$delta_preq; d[!is.finite(d) | d <= 0] <- 1
      lfo$mu_forecast <- lfo$mu_forecast * d
    } else {
      # NO delta_preq COLUMN: the hazard cannot be moved onto the recalibrated scale. Dropping
      # it is the only safe option — leaving it in place produced exactly the mixture this
      # function exists to prevent, a recalibrated p_invasion sitting beside a RAW
      # mu_forecast in the same frame, silently. Any consumer that needs the hazard now fails
      # loudly on a missing column instead of plotting two scales as one.
      warning("[scale] delta_preq absent: mu_forecast cannot be put on the '", col,
              "' scale and has been DROPPED, so it cannot be mixed with a recalibrated ",
              "p_invasion. Re-run 16b_invasion_recalibration.R to restore it.", call. = FALSE)
      lfo$mu_forecast <- NULL
    }
  }
  assign(name, lfo, envir = envir)
  message(sprintf("[scale] %s: drawing on '%s' (FORECAST_SCALE=%s)", name, col, scale))
  invisible(col)
}

#' Path of the current-forecast risk table for this scale.
fs_risk_csv <- function(out_dir, base = "bayes_risk_scores_all_zones.csv",
                        sub = "key_outputs", scale = FORECAST_SCALE) {
  p <- file.path(out_dir, sub, base)
  if (!fs_is_raw(scale)) return(p)
  r <- file.path(out_dir, sub, "raw", base)
  if (file.exists(r)) return(r)
  warning(sprintf(paste0("[scale] FORECAST_SCALE='raw' but %s does not exist; falling back to ",
                         "the primary (RECALIBRATED) table. Re-run the pipeline to write the ",
                         "raw twin."), r), call. = FALSE)
  p
}

#' Output directory for this scale: the primary one, or its `raw/` sibling.
fs_out_dir <- function(dir, scale = FORECAST_SCALE) {
  d <- if (fs_is_raw(scale)) file.path(dir, "raw") else dir
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' One-line scale statement for a figure caption.
#'
#' Stamped INTO the figure, not only into the path: once a PDF is dragged into a
#' manuscript the folder it came from is gone, and two directories holding identical
#' basenames is precisely how the wrong panel gets published.
#' @param lfo optionally the frame being drawn. SUPPLY IT whenever the caption sits beside
#'   LFO probabilities: the stamp then reports the column fs_lfo_col() actually resolved
#'   rather than the scale that was requested. Without it a frame that fell back to raw was
#'   still captioned "recalibrated", which is how every panel in both suites came to carry a
#'   caption contradicting its own numbers.
fs_stamp <- function(scale = FORECAST_SCALE, context = c("lfo", "current", "mixed"),
                     lfo = NULL) {
  context <- match.arg(context)
  if (fs_is_raw(scale))
    return("Probabilities: RAW (uncorrected model output).")
  if (!is.null(lfo) && identical(context, "lfo") &&
      identical(fs_lfo_col(lfo, scale = scale), "p_invasion"))
    return(paste("Probabilities: RAW (uncorrected model output) \u2014 the recalibrated scale",
                 "was requested but is unavailable on these rows."))
  switch(context,
    lfo     = paste("Probabilities: recalibrated on the prequential hazard-scale factor",
                    "(fitted only on folds whose outcome had already closed)."),
    current = paste("Probabilities: recalibrated on the pooled hazard-scale factor,",
                    "as deployed."),
    mixed   = paste("Probabilities: recalibrated — prequential factor for the",
                    "cross-validated panels, pooled (deployed) factor for the current forecast."))
}

# When to stamp the scale onto the figure itself. The raw twins share basenames with the
# primary set, so once a PDF leaves its folder nothing identifies it — that is where the
# mis-citation risk lives, and those are stamped by default. The primary figures are left
# visually unchanged so manuscript layout and hand-written captions are not disturbed;
# set FS_STAMP_FIGURES=always to label both sets, or =never to label neither.
FS_STAMP_FIGURES <- local({
  v <- tolower(trimws(Sys.getenv("FS_STAMP_FIGURES", "raw-only")))
  if (v %in% c("", "raw-only", "rawonly", "raw")) "raw-only"
  else if (v %in% c("always", "both")) "always"
  else if (v %in% c("never", "none", "off")) "never"
  else { warning(sprintf("[scale] unreadable FS_STAMP_FIGURES='%s'; using 'raw-only'.", v),
                 call. = FALSE); "raw-only" }
})

#' Add the scale statement to a figure, if this pass should be stamped.
#'
#' Handles both a patchwork assembly (plot_annotation merges with the existing
#' tag_levels rather than replacing it) and a bare ggplot (labs). Anything else is
#' returned untouched, so this can be dropped into a generic save helper.
fs_caption <- function(p, context = "mixed", scale = FORECAST_SCALE) {
  stamp <- switch(FS_STAMP_FIGURES,
                  always    = TRUE,
                  never     = FALSE,
                  `raw-only` = fs_is_raw(scale),
                  fs_is_raw(scale))
  if (!isTRUE(stamp) || is.null(p)) return(p)
  txt <- fs_stamp(scale, context)
  # FORCE THE CAPTION VISIBLE. Every house theme in this suite (theme_pub in the three make_*
  # suites, .rv_theme, theme_inv) sets plot.caption = element_blank(), so `labs(caption = )`
  # alone is computed and then discarded — which silently disabled this stamp on exactly the
  # single-ggplot panels it exists to protect (F2A, F2B, F_topk15_ever, ...), leaving a raw/
  # sibling indistinguishable from its primary twin. The patchwork branch was unaffected
  # because plot_annotation carries its own theme; give the ggplot branch the same.
  .cap_theme <- ggplot2::theme(
    plot.caption = ggplot2::element_text(size = 6.5, colour = "grey38", hjust = 0))
  if (inherits(p, "patchwork") && requireNamespace("patchwork", quietly = TRUE))
    return(p + patchwork::plot_annotation(caption = txt, theme = .cap_theme))
  if (inherits(p, "ggplot") && requireNamespace("ggplot2", quietly = TRUE))
    return(p + ggplot2::labs(caption = txt) + .cap_theme)
  p
}

message(sprintf("[scale] FORECAST_SCALE = %s (figure stamp: %s)",
                FORECAST_SCALE, FS_STAMP_FIGURES))
