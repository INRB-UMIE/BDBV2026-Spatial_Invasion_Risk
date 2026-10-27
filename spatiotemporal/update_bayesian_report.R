#!/usr/bin/env Rscript
# =============================================================================
# update_bayesian_report.R — refresh the dynamic RESULTS in
#   spatiotemporal/BAYESIAN_INVASION_REPORT.md
# from the pipeline's own output objects, so the report's numbers/tables track the
# latest run instead of being hand-edited (and drifting stale).
#
# WHAT IT UPDATES. Only the machine-derivable blocks, each delimited in the Markdown
# by an HTML-comment marker pair:
#     <!-- AUTOGEN:<key> -->  ... generated content ...  <!-- /AUTOGEN:<key> -->
# The surrounding methodology / interpretation PROSE is never touched. Re-running is
# idempotent (content is replaced between the same markers). If a marker pair is
# missing the block is skipped with a warning, so the script never corrupts the file.
#
# DATA SOURCES (all written by run_all.R):
#   outputs/diagnostics/invasion_evaluation.csv   leaderboard, support counts, folds
#   outputs/forecasts/bayes_parameters.rds|.csv   posterior covariate hazard ratios
#   outputs/forecasts/bayes_stacking_weights.rds  loo predictive-stacking weights
#   outputs/forecasts/bayes_risk_scores_current.rds  featured-model current forecast
#   outputs/reports/invasion_report.md            (only) the total confirmed-case count
#   00_config.R                                   ANALYSIS_DATE
#   outputs/key_outputs/model_selection.json      the featured model (read_featured_model())
#   outputs/key_outputs/detection_curves.csv      the pooled-recall comparison (targeting)
#
# The featured single Bayesian model is READ from the pipeline's own selection record,
# never recomputed here (see the note at `featured` below), and the covariate table is drawn
# from the featured model when it carries covariates, else the best AVAILABLE covariate
# model — so the report adapts when the config toggles change which models are fitted.
#
# USAGE:   Rscript spatiotemporal/update_bayesian_report.R
#          Rscript spatiotemporal/update_bayesian_report.R --check   # dry-run: report
#                  which blocks would change; write nothing; non-zero exit if any differ
# =============================================================================

suppressWarnings(suppressMessages({
  library(here); library(dplyr); library(readr); library(stringr); library(jsonlite)
}))

# %||% and %|NA|% come from 00_config.R (sourced below). This file used to define its own
# %||% that ALSO fell back on a scalar NA — a third semantics, installed unconditionally
# over its whole session — while eleven other modules defined the NULL-only form under a
# comment asserting they were all identical. The NA-dependent call sites here now say so by
# using %|NA|% explicitly.

ST_DIR <- file.path(here::here(), "spatiotemporal")
# Load config (paths, ANALYSIS_DATE, INVASION_SELECTION_EXCLUDE, invasion_topk_cols()) and
# figure_style.R for read_featured_model(), the one reader of the pipeline's selection record
# every figure suite uses. Both define objects only.
suppressWarnings(suppressMessages({
  source(file.path(ST_DIR, "00_config.R"))
  source(file.path(ST_DIR, "figure_style.R"))
}))
# Under a data-definition sensitivity OUT_DIR is redirected to outputs_sensitivity/<tag>/, but
# REPORT_PATH below is the PUBLISHED report: running here would write a sensitivity arm's
# numbers into it. Refuse, as 45/46 and the figure suites do.
stop_if_data_definition_sensitivity("update_bayesian_report.R")

# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and rewrites the published
# report. It is launched as its own subprocess. Guarded so a source() — a test, an ad-hoc
# check, a sweep that loads "every module" — defines what it needs and stops, instead of
# silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("update_bayesian_report.R")) {
  message("[report] sourced, not run: nothing is written. Use `Rscript update_bayesian_report.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

REPORT_PATH <- file.path(ST_DIR, "BAYESIAN_INVASION_REPORT.md")
DRY_RUN     <- any(c("--check", "--dry-run") %in% commandArgs(trailingOnly = TRUE))

# ---------------------------------------------------------------------------
# Small formatting helpers (match the report's typography: en-dash, × sign)
# ---------------------------------------------------------------------------
DASH <- "–"   # – (en dash), as used in the report's CrI ranges
TIMES <- "×"  # × (multiplication sign)
f2 <- function(x) formatC(x, format = "f", digits = 2)
f3 <- function(x) formatC(x, format = "f", digits = 3)
f1 <- function(x) formatC(x, format = "f", digits = 1)
cri <- function(lo, hi, d = 3) sprintf("[%s%s%s]", formatC(lo, format = "f", digits = d),
                                       DASH, formatC(hi, format = "f", digits = d))

# ---------------------------------------------------------------------------
# Load the pipeline outputs (each guarded; a missing file disables its blocks)
# ---------------------------------------------------------------------------
O <- OUT_DIR
read_rds_safe <- function(p) if (file.exists(p)) tryCatch(readRDS(p), error = function(e) NULL) else NULL
read_csv_safe <- function(p) if (file.exists(p)) tryCatch(readr::read_csv(p, show_col_types = FALSE),
                                                          error = function(e) NULL) else NULL

eval_tbl   <- read_csv_safe(file.path(O, "diagnostics", "invasion_evaluation.csv"))
risk_sc    <- read_rds_safe(file.path(O, "forecasts", "bayes_risk_scores_current.rds"))
weights    <- read_rds_safe(file.path(O, "forecasts", "bayes_stacking_weights.rds"))
params     <- read_rds_safe(file.path(O, "forecasts", "bayes_parameters.rds"))
if (is.null(params)) params <- read_csv_safe(file.path(O, "reports", "bayes_parameters.csv"))

if (is.null(eval_tbl) && is.null(risk_sc))
  stop("[update_report] no evaluation or risk-score outputs found under ", O,
       " — run the pipeline first.", call. = FALSE)

# ---------------------------------------------------------------------------
# Derived scalars used across blocks
# ---------------------------------------------------------------------------
# THE FEATURED MODEL IS READ, NOT RECOMPUTED. This used to re-run best_invasion_model() on the
# evaluation table — without the burn-in selection view, the sensitivity-arm exclusion and the
# recalibrated log-score axis the pipeline's own pick applies — and then to fall back to
# risk_sc$method[1], or to the literal "(featured model)", when that failed. So the report could
# name a model the pipeline did not feature, or name none and still publish. model_selection.json
# is the pipeline's record of what it picked; read_featured_model() has NO fallback, and a report
# that cannot tell which model was featured stops.
featured <- read_featured_model(file.path(O, "key_outputs", "model_selection.json"))

# THE SELECTION AXES ARE READ FROM THE SAME RECORD, not assumed. The composite that picks the
# featured model ranks on three columns that depend on two recorded settings:
#   selection_basis               "per_date_post_burnin": per-date AUC-PR skill, mean rank of
#                                 truth and log score on the post-burn-in rounds (the sel_*
#                                 columns of invasion_evaluation.csv); "per_date_all_rounds":
#                                 the same sel_* columns, computed over every round;
#                                 "pooled_all_dates": the pooled all-round columns. A record
#                                 written before the field existed (2026-10-14) was selected
#                                 under the latter.
#   selection_axis$log_score_axis "log_score_recal" (recalibrated) or "log_score" (raw).
# The leaderboard used to rank and print the pooled, all-round, RAW columns beside a label
# saying "best CV composite", so the table did not show the numbers the pick was made on.
SEL <- local({
  j <- jsonlite::fromJSON(file.path(O, "key_outputs", "model_selection.json"), simplifyVector = TRUE)
  basis <- j$selection_basis
  if (is.null(basis)) {
    message("[update_report] model_selection.json carries no selection_basis (written before ",
            "2026-10-14): reading it as that rule, pooled_all_dates.")
    basis <- "pooled_all_dates"
  }
  ls_axis <- j$selection_axis$log_score_axis %||% "log_score"
  if (length(basis) != 1L || !basis %in% c("per_date_post_burnin", "per_date_all_rounds", "pooled_all_dates"))
    stop("[update_report] unrecognised selection_basis in model_selection.json: ",
         paste(basis, collapse = ", "), call. = FALSE)
  if (length(ls_axis) != 1L || !ls_axis %in% c("log_score", "log_score_recal"))
    stop("[update_report] unrecognised selection_axis$log_score_axis in model_selection.json: ",
         paste(ls_axis, collapse = ", "), call. = FALSE)
  # `post` = the selection reads the per-date sel_* columns (on the post-burn-in rounds, or on
  # every round under "per_date_all_rounds"); `burnin` = those rounds exclude the burn-in.
  post   <- basis %in% c("per_date_post_burnin", "per_date_all_rounds")
  burnin <- identical(basis, "per_date_post_burnin")
  rounds <- if (burnin) "post-burn-in rounds" else "all rounds"
  list(basis = basis, ls_axis = ls_axis, post = post, burnin = burnin, rounds = rounds,
       skill = if (post) "sel_auc_pr_skill_date" else "auc_pr_skill",
       rot   = if (post) "sel_mean_rank_of_truth" else "mean_rank_of_truth",
       ls    = paste0(if (post) "sel_" else "", ls_axis),
       skill_lab = if (post) sprintf("Per-date AUC-PR skill, %s", rounds)
                   else "Pooled AUC-PR skill, all rounds",
       rot_lab   = sprintf("Mean rank of truth, %s", rounds),
       ls_lab    = sprintf("%s log score, %s",
                           if (identical(ls_axis, "log_score_recal")) "Recalibrated" else "Raw", rounds))
})
#' The selection-axis columns the evaluation table must carry; NULL when they are all there.
sel_missing <- function(ev) {
  m <- setdiff(c(SEL$skill, SEL$rot, SEL$ls), names(ev))
  if (length(m)) m else NULL
}
#' The leaderboard's rows: the six candidates leading on the selection's AUC-PR axis at h=1,
#' PLUS the featured model wherever it ranks (it is chosen by the composite of three axes over
#' both horizons, not by this one). ONE definition, shared by the table and by the
#' featured-pick sentence's "tightest calibration among the leaders" check.
leader_set <- function(ev) {
  e <- candidate_field(ev)
  if (!nrow(e) || !is.null(sel_missing(e))) return(NULL)
  e <- e %>% dplyr::arrange(dplyr::desc(.data[[SEL$skill]]))
  dplyr::bind_rows(utils::head(e, 6), e %>% dplyr::filter(method == featured)) %>%
    dplyr::distinct(method, .keep_all = TRUE)
}

#' The CANDIDATE field: every Bayesian single model the TWO-STAGE selection chose among — the
#' stage-1 fixed-beta candidates plus the stage-2 pool (the stage-1 pick's cross-validated
#' time-varying forms, model_selection.json's stage2_pool), which INVASION_SELECTION_EXCLUDE
#' removes from stage 1. The generation-time sensitivity arms are never candidates, and
#' ensembles are combinations of the grid, not members of it; partial-CV rows are scored on
#' fewer cells. Ranking any of them beside the featured model presents a refit or an ensemble as
#' a rival the selection passed over; leaving a stage-2 rival out judges the featured form
#' against a field it did not face. One definition, used by every block below.
.sel_excl <- get0("INVASION_SELECTION_EXCLUDE", ifnotfound = NULL)
if (is.null(.sel_excl)) stop("[update_report] INVASION_SELECTION_EXCLUDE is not defined by 00_config.R",
                             call. = FALSE)
.stage2_pool <- tryCatch(unlist(jsonlite::fromJSON(file.path(O, "key_outputs", "model_selection.json"),
                                                   simplifyVector = TRUE)$featured$bayesian$stage2_pool),
                         error = function(e) NULL)
candidate_field <- function(ev) {
  e <- ev %>% dplyr::filter(grepl("^Bayes", method), !grepl("-ens-", method),
                            !grepl(.sel_excl, method) | method %in% c(featured, .stage2_pool))
  if ("partial_cv" %in% names(e)) e <- e %>% dplyr::filter(!(partial_cv %in% TRUE))
  e
}

# h=1 evaluation support (from the featured model's row, or the first full-CV row)
ev1 <- if (!is.null(eval_tbl)) eval_tbl %>% dplyr::filter(horizon == 1) else NULL
ev2 <- if (!is.null(eval_tbl)) eval_tbl %>% dplyr::filter(horizon == 2) else NULL
supp_row <- function(ev) {
  if (is.null(ev) || !nrow(ev)) return(NULL)
  e <- if ("partial_cv" %in% names(ev)) ev %>% dplyr::filter(!(partial_cv %in% TRUE)) else ev
  if (!nrow(e)) e <- ev
  e[1, ]
}
# The support counts are shared by every method (one scoring support), so the featured model's
# own row is read when present; supp_row() is the fallback for a table that lacks it.
.feat_row <- function(ev) {
  if (is.null(ev)) return(NULL)
  r <- ev[ev$method == featured, , drop = FALSE]
  if (nrow(r) == 1L) r else supp_row(ev)
}
s1 <- .feat_row(ev1); s2 <- .feat_row(ev2)

# ---------------------------------------------------------------------------
# Marker-replacement engine (string-safe, no regex backreference pitfalls)
# ---------------------------------------------------------------------------
inject <- function(txt, key, content) {
  s <- sprintf("<!-- AUTOGEN:%s -->", key); e <- sprintf("<!-- /AUTOGEN:%s -->", key)
  i <- regexpr(s, txt, fixed = TRUE); j <- regexpr(e, txt, fixed = TRUE)
  if (i[1] < 0L || j[1] < 0L) { warning(sprintf("[update_report] marker '%s' not found — skipped.", key)); return(txt) }
  if (j[1] < i[1]) { warning(sprintf("[update_report] marker '%s' end precedes start — skipped.", key)); return(txt) }
  before <- substr(txt, 1L, i[1] + attr(i, "match.length") - 1L)
  after  <- substr(txt, j[1], nchar(txt))
  paste0(before, "\n", content, "\n", after)
}

# ---------------------------------------------------------------------------
# Block generators — each returns the Markdown to sit between its markers, or NULL
# to leave the block unchanged (missing inputs).
# ---------------------------------------------------------------------------
blocks <- list()
# The set of keys to inject is taken from the MARKERS IN THE REPORT (see the Apply section),
# not from names(blocks). NEEDED because `blocks$k <- local({...})` DELETES the element when
# the generator returns NULL (R drops a NULL assigned with `$`), so iterating names(blocks)
# makes the "no data" branch unreachable: a block that silently failed to regenerate is
# indistinguishable from one that was never meant to change, and the report ships whatever
# stale text the marker still holds with no warning anywhere.

# --- data_summary: intro one-liner (dates, at-risk/total, confirmed/affected) ------
blocks$data_summary <- local({
  if (is.null(risk_sc)) return(NULL)
  h1 <- risk_sc %>% dplyr::filter(horizon == 1)
  n_zones  <- dplyr::n_distinct(risk_sc$health_zone)
  affected <- sum(h1$was_active_before, na.rm = TRUE)
  at_risk  <- n_zones - affected
  ad  <- as.character(get0("ANALYSIS_DATE", ifnotfound = ""))
  # Total confirmed cases: the pipeline's own freshly-written invasion_report.md header
  # ("**Confirmed cases:** N across M affected") is the authoritative count; parse it,
  # else leave a placeholder note rather than fabricate.
  conf <- NA_integer_
  irep <- file.path(O, "reports", "invasion_report.md")
  if (file.exists(irep)) {
    m <- stringr::str_match(paste(readLines(irep, warn = FALSE), collapse = "\n"),
                            "Confirmed cases:\\*\\*\\s*([0-9,]+)\\s+across")
    if (!is.na(m[1, 2])) conf <- as.integer(gsub(",", "", m[1, 2]))
  }
  conf_txt <- if (is.na(conf)) "the confirmed" else format(conf, big.mark = ",")
  # The CUTOFF is not the analysis date. This line filled both slots from `ad`, so it always
  # asserted they coincide; the pipeline records them separately in run_info.json
  # (linelist_cutoff_date = the last week-anchor start used for training, training_window_end =
  # the last day carrying training data), and on the shipped run they differ (2026-09-01 vs
  # 2026-09-07). They coincide only when the final week carries data, which run_all.R warns
  # about precisely because it is not guaranteed.
  .cut <- local({
    ri <- file.path(O, "key_outputs", "run_info.json")
    if (!file.exists(ri)) return(NA_character_)
    v <- tryCatch(jsonlite::fromJSON(ri, simplifyVector = TRUE)$linelist_cutoff_date,
                  error = function(e) NULL)
    if (is.null(v) || length(v) != 1L || is.na(v)) NA_character_ else as.character(v)
  })
  sprintf(paste0("The live forecast (analysis date %s, training cutoff %s) covers **%d at-risk zones** ",
                 "out of %d, with %s confirmed cases across %d affected zones."),
          ad, if (is.na(.cut)) "not recorded" else .cut, at_risk, n_zones, conf_txt, affected)
})

# --- sample_events: §1 sample-size caveat (events + base rates) --------------------
blocks$sample_events <- local({
  if (is.null(s1)) return(NULL)
  br1 <- 100 * (s1$base_rate %||% NA); br2 <- if (!is.null(s2)) 100 * (s2$base_rate %||% NA) else NA
  # THE EVENTS ARE DISTINCT ZONES. It used to say "~N events from a handful of distinct zones",
  # as though a few zones were invaded repeatedly; an invasion is a zone's FIRST confirmed case,
  # so each event is a different zone. n_invasion_events counts them (it equals the positive-row
  # count at h=1, where the windows are disjoint).
  n_ev <- as.integer(s1$n_invasion_events %||% s1$n_invasions %||% NA)
  sprintf(paste0("Pooled over the leave-future-out folds there are only %d realised invasions at ",
                 "*h*=1 — %d distinct health zones, each counted once, at its first confirmed case — ",
                 "against a per-zone-week base rate of ",
                 "~%s %% at *h*=1 (~%s %% at *h*=2). Every ranking is therefore provisional and ",
                 "differences between good models are frequently within Monte-Carlo noise."),
          as.integer(s1$n_invasions %||% NA), n_ev, f2(br1),
          if (is.na(br2)) "?" else f1(br2))
})

# --- folds_h1: §4 number of h=1 folds ---------------------------------------------
blocks$folds_h1 <- local({
  if (is.null(s1) || is.null(s1$n_folds)) return(NULL)
  sprintf("About **%d folds** survive for *h*=1.", as.integer(s1$n_folds))
})

# --- support_pool: §4 pooled at-risk zone-weeks / events / base rate ---------------
blocks$support_pool <- local({
  if (is.null(s1)) return(NULL)
  sprintf("The live *h*=1 comparison pools %s at-risk zone-weeks with **%d invasion events** (base rate %s %%).",
          format(as.integer(s1$n_atrisk %||% NA), big.mark = ","),
          as.integer(s1$n_invasions %||% NA), f2(100 * (s1$base_rate %||% NA)))
})

# --- leaderboard: §5 AUC-PR-skill table (h=1, Bayesian single models, full-CV) -----
blocks$leaderboard <- local({
  if (is.null(ev1)) return(NULL)
  if (!nrow(candidate_field(ev1))) return(NULL)
  # The model_selection.json and the evaluation table come from the same run, so a missing
  # selection-axis column means they do not; the block is left stale, loudly, rather than
  # filled with columns the selection did not use.
  .miss <- sel_missing(candidate_field(ev1))
  if (!is.null(.miss)) {
    warning(sprintf(paste0("[update_report] invasion_evaluation.csv lacks the selection-axis ",
                           "column(s) %s named by model_selection.json (selection_basis %s, ",
                           "log_score_axis %s): the leaderboard is not regenerated."),
                    paste(.miss, collapse = ", "), SEL$basis, SEL$ls_axis), call. = FALSE)
    return(NULL)
  }
  e <- leader_set(ev1)
  num <- function(x, d) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "n/a")
  # The pooled, all-round AUC-PR skill [90 % zone-cluster bootstrap CI OF THE SKILL: each
  # replicate's AP over its own base rate, auc_pr_skill_lo/hi], on the raw probabilities. The AP
  # interval divided by the point base rate is not that interval (the base rate is resampled too).
  have_ci <- all(c("auc_pr_skill_lo", "auc_pr_skill_hi") %in% names(e))
  pooled_skill <- function(i)
    if (have_ci && is.finite(e$auc_pr_skill_lo[i]) && is.finite(e$auc_pr_skill_hi[i]))
      sprintf("%s [%s%s%s]", f1(e$auc_pr_skill[i]),
              f1(e$auc_pr_skill_lo[i]), DASH, f1(e$auc_pr_skill_hi[i]))
    else num(e$auc_pr_skill[i], 1)
  all_skill <- function(i)
    if (all(c("auc_pr_skill_date_lo", "auc_pr_skill_date_hi") %in% names(e)) &&
        is.finite(e$auc_pr_skill_date_lo[i]) && is.finite(e$auc_pr_skill_date_hi[i]))
      sprintf("%s [%s%s%s]", f1(e$auc_pr_skill_date[i]),
              f1(e$auc_pr_skill_date_lo[i]), DASH, f1(e$auc_pr_skill_date_hi[i]))
    else num(e$auc_pr_skill_date[i], 1)
  # SELECTION-AXIS COLUMNS FIRST, then the comparison-only ones. Under the pooled basis the
  # selection's AUC-PR axis IS the pooled skill, so it is shown once, with its interval.
  hdr <- c("Model",
           sprintf("%s%s (selection axis)", SEL$skill_lab, if (SEL$post) "" else " [90 % CI]"),
           sprintf("%s (selection axis)", SEL$rot_lab),
           sprintf("%s (selection axis)", SEL$ls_lab),
           # The value the figures REPORT (per-date skill over all rounds) when the selection
           # read a subset of the rounds; under the other bases the selection axis already is it.
           if (SEL$burnin) "Per-date AUC-PR skill, ALL rounds [90 % CI] (as reported; not a selection axis)",
           "Calibration-in-the-large, raw, all rounds (comparison only)")
  rows <- vapply(seq_len(nrow(e)), function(i) {
    name <- if (identical(e$method[i], featured))
      sprintf("**%s** (featured, two-stage CV composite)", e$method[i]) else e$method[i]
    cells <- c(name,
               if (SEL$post) num(e[[SEL$skill]][i], 1) else pooled_skill(i),
               num(e[[SEL$rot]][i], 1),
               num(e[[SEL$ls]][i], 4),
               if (SEL$burnin) all_skill(i),
               paste0(num(e$calibration_in_large[i], 1), TIMES))
    paste0("| ", paste(cells, collapse = " | "), " |")
  }, character(1))
  # What a uniformly random ordering scores on the per-date axis (NOT 1; see .auc_pr_by_date()).
  chance <- if (SEL$post && "sel_auc_pr_skill_date_chance" %in% names(e)) {
    ch <- e$sel_auc_pr_skill_date_chance[is.finite(e$sel_auc_pr_skill_date_chance)]
    if (!length(ch)) "" else if (diff(range(ch)) < 0.05)
      sprintf(" On the per-date axis a uniformly random ordering scores about %s on these rounds, not 1.", f1(ch[1]))
    else sprintf(" On the per-date axis a uniformly random ordering scores %s%s%s on these rounds, not 1.",
                 f1(min(ch)), DASH, f1(max(ch)))
  } else ""
  note <- sprintf(paste0("Selection axes as recorded in `model_selection.json` (basis `%s`, log-score ",
                         "axis `%s`), shown at *h*=1; the composite sums the within-horizon ranks of all ",
                         "three over both horizons. Columns marked comparison only or as reported are ",
                         "not used by the selection.%s%s"),
                  SEL$basis, SEL$ls_axis, chance,
                  if (SEL$burnin) paste0(" The selection deliberately excludes the burn-in rounds, in ",
                                         "which where first cases were confirmed partly reflected ",
                                         "where surveillance looked first; performance is reported ",
                                         "over all rounds.") else "")
  paste(c(paste0("| ", paste(hdr, collapse = " | "), " |"),
          paste0("|", paste(rep("---", length(hdr)), collapse = "|"), "|"),
          rows, "", note), collapse = "\n")
})

# --- featured_pick: §5 sentence naming the featured model + its headline stats ------
blocks$featured_pick <- local({
  if (is.null(ev1)) return(NULL)
  fr <- ev1 %>% dplyr::filter(method == featured)
  fr2 <- if (!is.null(ev2)) ev2 %>% dplyr::filter(method == featured) else NULL
  if (!nrow(fr)) return(NULL)
  # THE SENTENCE QUOTES THE SELECTION AXES (SEL, the leaderboard's columns). It used to quote the
  # pooled, all-round, raw columns while describing the composite, i.e. numbers the pick was not
  # made on. Without the selection-axis columns the block is left stale, loudly (the leaderboard
  # block warns why), rather than filled with the wrong ones.
  lead <- leader_set(ev1)
  if (is.null(lead) || !is.null(sel_missing(fr))) return(NULL)
  num <- function(x, d) if (length(x) == 1L && is.finite(x)) formatC(x, format = "f", digits = d) else "n/a"
  skill_txt <- if (SEL$post)
    sprintf("per-date AUC-PR skill %s on the %s%s%s", num(fr[[SEL$skill]][1], 1), SEL$rounds,
            if ("sel_auc_pr_skill_date_chance" %in% names(fr) && is.finite(fr$sel_auc_pr_skill_date_chance[1]))
              sprintf(", against %s for a random ordering", num(fr$sel_auc_pr_skill_date_chance[1], 1)) else "",
            if (SEL$burnin && "auc_pr_skill_date" %in% names(fr))
              sprintf(" (%s over all rounds, the value the figures report)", num(fr$auc_pr_skill_date[1], 1))
            else "")
    else sprintf("pooled AUC-PR skill %s%s over all rounds", num(fr[[SEL$skill]][1], 1), TIMES)
  calib <- num(fr$calibration_in_large[1], 1)
  # CHECK the superlative instead of asserting it. "the tightest calibration among the leaders"
  # was hard-coded prose: on the shipped run Bayes-M16-fill-med was tighter (1.991 vs 2.033),
  # and the leaderboard printed immediately above rounds both to 2.0x, so the table visibly tied
  # while the prose claimed a win. Calibration-in-the-large is best at 1, so "tightest" means
  # smallest |x - 1| among the same models the leaderboard block lists (leader_set()).
  calib_word <- local({
    d <- abs(lead$calibration_in_large - 1)
    if (!nrow(lead) || !any(is.finite(d))) "competitive"
    else if (identical(lead$method[which.min(d)], featured)) "the tightest"
    else "competitive"
  })
  rank2 <- if (!is.null(fr2) && nrow(fr2) && SEL$rot %in% names(fr2))
    sprintf("%s, %s", num(fr2[[SEL$rot]][1], 1), SEL$rounds) else "?"
  axes_txt <- sprintf("%s AUC-PR skill, mean rank of truth and the %s log score, all %s",
                      if (SEL$post) "per-date" else "pooled",
                      if (identical(SEL$ls_axis, "log_score_recal")) "recalibrated" else "raw",
                      if (SEL$burnin) "on the post-burn-in rounds" else "over all rounds")
  # NEVER fabricate the number. The fallback used to be the literal "0.99", so a run whose
  # evaluation table lacked auc_roc would publish a discrimination figure nothing computed
  # (and one that rounds the real value, ~0.98, upward).
  # SAME FIELD AS THE LEADERBOARD. max() ran over EVERY h=1 row — the three structural
  # baselines and any partial_cv rows included — so "best AUC-ROC" could be a rank-only
  # comparator's number published as the Bayesian suite's discrimination. It lands on a
  # Bayesian model on the current frame by luck, not by construction.
  # The candidate field — sensitivity arms and ensembles out — as the leaderboard uses.
  .ev1_elig <- candidate_field(ev1)
  if (!nrow(.ev1_elig)) return(NULL)
  auroc <- if ("auc_roc" %in% names(.ev1_elig) && any(is.finite(.ev1_elig$auc_roc)))
             f2(max(.ev1_elig$auc_roc, na.rm = TRUE)) else "not computed in this run"
  .st <- tryCatch(jsonlite::fromJSON(file.path(O, "key_outputs", "model_selection.json"),
                                     simplifyVector = TRUE)$featured$bayesian, error = function(e) NULL)
  .s1 <- .st$stage1_method %||% featured
  # Stage 2 compared anything only when a form was cross-validated (stage2_pool holds the stage-1
  # model plus its forms); otherwise stage 1's pick is featured with no stage-2 comparison.
  .n_pool <- length(unlist(.st$stage2_pool))
  .stage_txt <- if (identical(.s1, featured) && .n_pool <= 1L)
    sprintf("stage 1 picks **%s** among the fixed-beta candidates; no time-varying form was cross-validated, so stage 1's pick is featured", .s1)
  else if (identical(.s1, featured))
    sprintf("stage 1 picks **%s** among the fixed-beta candidates and stage 2 keeps it against its time-varying forms", .s1)
  else sprintf("stage 1 picks the specification **%s** among the fixed-beta candidates and stage 2, against its time-varying forms, features **%s**", .s1, featured)
  .sec <- unlist(.st$secondary_forms)
  if (length(.sec))
    .stage_txt <- paste0(.stage_txt, sprintf(paste0(" (the frailty forms %s are a secondary analysis, ",
                                                    "never candidates)"),
                                             paste(sprintf("`%s`", .sec), collapse = ", ")))
  # NO UNCHECKED ADJECTIVES: "solid" discrimination and a "strong" mean rank were asserted; the
  # featured model's place in the candidate field on each axis is computed and stated instead.
  .place <- function(ev, col, best_low) {
    e <- candidate_field(ev); if (!nrow(e) || !col %in% names(e)) return("")
    v <- e[[col]]; me <- v[e$method == featured][1]
    if (!is.finite(me)) return("")
    sprintf(", %s of %d candidates", .ordinal(1L + sum(if (best_low) v < me else v > me, na.rm = TRUE)),
            sum(is.finite(v)))
  }
  .ordinal <- function(k) paste0(k, if (k %% 100L %in% 11:13) "th" else
                                   c("th", "st", "nd", "rd", rep("th", 6))[k %% 10L + 1L])
  sprintf(paste0("The featured model is chosen in two stages by the **CV composite** (summed within-horizon ranks of ",
    "%s, summed over both horizons; the basis is recorded in `model_selection.json`): ",
    gsub("%", "%%", .stage_txt, fixed = TRUE), ". **%s** has ",
    "*h*=1 %s%s; %s calibration among the leaders (raw calibration-in-the-large over all rounds %s%s, ",
    "not a selection axis); and a mean rank of truth at *h*=2 of %s%s. Best pooled AUC-ROC over all rounds ",
    "among the candidate models %s %s. The **loo stacking weights** are ",
    "used only to build the loo-stacked ENSEMBLE (`bayes_ensemble_*`), NOT to pick the featured single model."),
    axes_txt, featured, skill_txt, .place(ev1, SEL$skill, FALSE), calib_word, calib, TIMES, rank2,
    if (!is.null(ev2)) .place(ev2, SEL$rot, TRUE) else "", "≈", auroc)
})

# --- stacking_weights: §5 the loo-stacking weight spread ---------------------------
blocks$stacking_weights <- local({
  if (is.null(weights) || !length(weights)) return(NULL)
  w <- sort(weights[is.finite(weights)], decreasing = TRUE)
  top <- w[w >= 0.01]; if (!length(top)) top <- head(w, 3)
  lst <- paste(sprintf("%s %s", names(top), f2(top)), collapse = ", ")
  rest <- if (length(w) > length(top)) sprintf(" (others %s %s)", "≤", f2(max(w[!(names(w) %in% names(top))]))) else ""
  sprintf(paste0("In this run the loo predictive-stacking weights spread across %s%s — so the loo-stacked ",
                 "ENSEMBLE mixes several structures rather than concentrating on the single featured model."),
          lst, rest)
})

# --- covariate_table: §5 posterior covariate hazard ratios (best covariate model) --
blocks$covariate_table <- local({
  if (is.null(params) || !"is_intercept" %in% names(params)) return(NULL)
  # Covariate coefficients only: a time-varying form's process parameters (frailty_theta,
  # rw1_sigma_week) are not hazard ratios per +1 SD and carry their own effect_scale.
  if ("effect_scale" %in% names(params))
    params <- params %>% dplyr::filter(is_intercept | effect_scale %in% c("hazard ratio", "odds ratio"))
  cov_models <- params %>% dplyr::filter(!is_intercept) %>% dplyr::distinct(model) %>% dplyr::pull(model)
  if (!length(cov_models)) return(NULL)
  # Prefer the featured model if it carries covariates; else the "full" exogenous model;
  # else the first covariate model available.
  pick <- if (featured %in% cov_models) featured
          else if (any(grepl("-full$", cov_models))) grep("-full$", cov_models, value = TRUE)[1]
          else cov_models[1]
  d <- params %>% dplyr::filter(model == pick, !is_intercept, is.finite(hr), hr > 0)
  if (!nrow(d)) return(NULL)
  # THE READINGS ARE DERIVED FROM THE FIT, not written beside it. They used to be fixed text per
  # term, and on LINELIST_09102026 the table printed d_min HR 0.97 [0.73–1.26] next to "Farther
  # from the frontier ⇒ lower risk (the strongest, most interpretable driver)": an interval that
  # spans 1, described as the strongest effect. Each reading now states the DIRECTION from the
  # sign of log HR, says whether the 90 % CrI excludes 1, and calls a term the largest effect
  # only when it has the largest |log HR| of the model's covariates AND its interval excludes 1.
  # Only what "a higher value" MEANS for each term is written by hand (covariates enter
  # standardised, so the HR is per +1 SD).
  higher <- c(d_min = "longer travel time to the nearest affected zone",
              ccvi  = "more deprivation",
              log_pop = "larger population (conditional on the mobility offset)",
              healthsite_density = "more health facilities")
  label <- c(d_min = "`d_min` (travel time to nearest affected zone)", ccvi = "`ccvi` (deprivation)",
             log_pop = "`log_pop`", healthsite_density = "`healthsite_density`")
  d <- d %>% dplyr::mutate(lhr = abs(log(hr)),
                           excl = is.finite(lo) & is.finite(hi) & (lo > 1 | hi < 1)) %>%
    dplyr::arrange(dplyr::desc(lhr))
  top <- d$term[which.max(d$lhr)]
  reading_of <- function(i) {
    t   <- d$term[i]
    eff <- sprintf("%s ⇒ %s risk", higher[t] %|NA|% sprintf("higher `%s`", t),
                   if (d$hr[i] > 1) "higher" else "lower")
    if (d$excl[i])
      paste0(toupper(substr(eff, 1, 1)), substr(eff, 2, nchar(eff)), " (90 % CrI excludes 1)",
             if (identical(t, top) && nrow(d) > 1L) "; the largest effect among the covariates" else "")
    else
      sprintf("No credible effect: the 90 %% CrI spans 1 (point estimate: %s)", eff)
  }
  rows <- vapply(seq_len(nrow(d)), function(i) {
    t <- d$term[i]
    sprintf("| %s | %s %s | %s | %s |", label[t] %|NA|% sprintf("`%s`", t),
            f2(d$hr[i]), cri(d$lo[i], d$hi[i], 2), f2(d$p_dir[i]), reading_of(i))
  }, character(1))
  hdr <- sprintf("| Covariate | %s HR [90 %% CrI] | P(HR>1) | Reading |", pick)
  paste(c(hdr, "|---|---|---|---|", rows), collapse = "\n")
})

# --- top_zones: §5 highest-risk at-risk zones (featured model, h=1) -----------------
blocks$top_zones <- local({
  if (is.null(risk_sc)) return(NULL)
  d <- risk_sc %>% dplyr::filter(horizon == 1, !was_active_before, !is.na(p_case_invasion)) %>%
    dplyr::arrange(dplyr::desc(p_case_invasion)) %>% head(5)
  if (!nrow(d)) return(NULL)
  rows <- vapply(seq_len(nrow(d)), function(i)
    sprintf("| %s | %s | %s %s | %s%s |", d$health_zone[i], d$province[i] %|NA|% "",
            f3(d$p_case_invasion[i]), cri(d$p_lo[i], d$p_hi[i], 3),
            f1(d$rr_nat[i]), TIMES), character(1))
  .ps <- if ("prob_scale" %in% names(d)) unique(as.character(d$prob_scale)) else NA_character_
  .ps_txt <- if (identical(.ps, "recalibrated-deployed")) "recalibrated, as deployed"
             else if (identical(.ps, "raw")) "raw, uncorrected"
             else stop("[update_report] top_zones: the risk table's prob_scale is ",
                       paste(.ps, collapse = "/"), "; refusing to label it.", call. = FALSE)
  cap <- sprintf(paste0("**Highest-risk at-risk zones (%s, next week, posterior mean [90 %% CrI], ", .ps_txt, ").** ",
    "These are the top of the national relative-invasion-risk map (**Figure 3A** of `key_outputs/figures/Figure3`; the rank map is `bayes_invasion_rank_map_national_h*`), and their posterior spread ",
    "is shown zone-by-zone in **Figure 3C** (1- vs 2-week merged); the companion uncertainty map ",
    "(**Figure 3B**) locates where those forecasts are least certain."), featured)
  paste(c(cap, "",
          "| Zone | Province | P(case) [90 % CrI] | Rel. risk (nat.) |",
          "|---|---|---|---|", rows), collapse = "\n")
})

# --- featured_name: §3.5 the "In the current run this is <model>" sentence ---------
# Decode a Bayes-<...> label into a short human description so the sentence stays true
# to the actually-featured model (kernel family / GT / covariates / road-distance).
describe_bayes <- function(lbl) {
  # Use the canonical parser (00_config.R): a local regex dropped the -fill / -split tokens
  # and described a filled composite as its unfilled parent.
  kern <- mobility_kernel_from_method(lbl)
  base <- if (is.na(kern)) sub("^Bayes-", "", lbl)
          else sub("-(dist|fill|split)", "", gsub("-(dist|fill|split)", "", kern))
  fam  <- c(M4 = "gravity", M4c = "cohort-calibrated gravity", M8 = "short-trip + gravity composite",
            M9 = "multi-kernel ensemble", M10 = "short-trip + radiation composite",
            M11 = "inward meeting-location FOI", M13 = "cohort + gravity composite",
            M13c = "cohort + cohort-calibrated gravity composite",
            M14 = "cohort + radiation composite", M15 = "symmetrised relocation OD",
            M16 = "cohort + relocation OD composite",
            M17 = "gravity + cohort-gravity + radiation consensus")[base]
  fam  <- if (is.na(fam)) base else fam
  gt   <- if (grepl("-(gt)?short", lbl)) "short" else if (grepl("-(gt)?long", lbl)) "long" else "medium"
  cov  <- if (grepl("full-susp", lbl)) "full exogenous + suspected-case covariates"
          else if (grepl("-susp", lbl)) "suspected-case covariates"
          else if (grepl("-full", lbl)) "full exogenous covariates"
          else if (grepl("-geo", lbl))  "geo covariates"
          else "no covariates"
  qual <- paste0(
    if (!is.na(kern) && grepl("-dist", kern)) " on road distance" else "",
    if (!is.na(kern) && grepl("-split", kern)) ", cohort rows split per origin" else "",
    if (!is.na(kern) && grepl("-(fill|split)", kern)) ", with source-cell fill" else "")
  .f <- tv_form_of(lbl)
  form <- if (is.na(.f)) ""
          else sprintf(", import coefficient with %s",
                       c(rw1 = "a random-walk weekly drift", frailty = "a gamma zone frailty",
                         rw1_frailty = "a random-walk weekly drift and a gamma zone frailty")[[TV_FORMS[[.f]]$form]])
  sprintf("the %s kernel%s with the %s generation-time profile, %s%s", fam, qual, gt, cov, form)
}
blocks$featured_name <- local({
  if (is.null(featured) || is.na(featured)) return(NULL)
  sprintf(paste0("In the current run this is **%s** (%s). Ensembles are excluded from this pick so ",
                 "the featured single model always has a current forecast."),
          featured, describe_bayes(featured))
})

# --- targeting: catch-rate at the published watch-list size vs random --------------
# THE BUDGET IS RESOLVED FROM THE TABLE, not hard-coded. This block required a literal
# `recall_at_10` column and returned NULL when it was absent, so changing TOPK_VALUES
# deleted the whole paragraph from the report with nothing saying so.
blocks$targeting <- local({
  if (is.null(ev1)) return(NULL)
  fr <- ev1 %>% dplyr::filter(method == featured)
  if (!nrow(fr) || !all(c("n_atrisk", "n_folds") %in% names(fr))) return(NULL)
  tk <- invasion_topk_cols(fr)
  if (is.null(tk)) return(NULL)
  rec <- fr[[tk$recall]][1]
  per_fold <- (fr$n_atrisk[1] %|NA|% NA) / (fr$n_folds[1] %|NA|% NA)   # at-risk zones per fold
  rnd <- if (is.finite(per_fold) && per_fold > 0) min(tk$k / per_fold, 1) else NA
  if (!is.finite(rec) || !is.finite(rnd) || rnd <= 0) return(NULL)
  # THE POOLED COMPARISON IS COMPUTED, NOT ASSERTED. The parenthesis used to say "the panel's
  # pooled curve reads slightly higher" as fixed text. The panel now draws the per-fold MEAN
  # (detection_curves.csv, the estimand quoted here), and whether pooling over folds reads
  # higher or lower is a property of the run: it is read off the same file's recall_pooled.
  pooled_txt <- local({
    fp <- file.path(O, "key_outputs", "detection_curves.csv")
    dc <- if (file.exists(fp)) read_csv_safe(fp) else NULL
    if (is.null(dc) || !"recall_pooled" %in% names(dc)) return("")
    r <- dc[dc$method == featured & dc$horizon == 1 & dc$k == tk$k, , drop = FALSE]
    if (nrow(r) != 1L || !is.finite(r$recall_pooled)) return("")
    rp <- r$recall_pooled
    rel <- if (round(100 * rp) == round(100 * rec)) "the same, to the nearest per cent"
           else if (rp > rec) "higher" else "lower"
    sprintf("; pooled over folds instead it is %.0f %%, %s", 100 * rp, rel)
  })
  sprintf(paste0("**Operationally useful targeting (the prioritisation panel, Figure 2B of `key_outputs/figures/Figure2_labelled`).** A top-%d watch-list catches ~%.0f %% ",
                 "of the next invasions per round (mean over folds, as the panel draws it%s) versus ~%.0f %% for a random list of the same size (a ~%.0f%s lift)."),
          tk$k, 100 * rec, pooled_txt, 100 * rnd, rec / rnd, TIMES)
})

# --- calibration: §5 calibration-in-the-large for the FEATURED model ---------------
# Replaces a hand-written "~2-3x (2.6x for Bayes-M10-med)": both the magnitude and the model
# name were literals that no longer tracked the featured pick or the current folds.
# calibration_in_large = mean(predicted p) / observed base rate over the leave-future-out rows,
# so >1 is over-prediction; both scales are quoted (recalibrated = after the prequential delta).
blocks$calibration <- local({
  # Per horizon, from the evaluation table: invasion_recalibration.csv's cal_in_large is the RAW
  # ratio on the settled post-burn-in 1-week rows that fit delta, copied onto every horizon's row
  # (fit_horizon = 1), so quoting it per horizon reported the h=1 fitting ratio twice.
  if (is.null(featured) || is.na(featured) || is.null(eval_tbl)) return(NULL)
  need <- c("method", "horizon", "calibration_in_large", "calibration_in_large_recal")
  if (!all(need %in% names(eval_tbl))) return(NULL)
  r <- eval_tbl %>% dplyr::filter(method == featured, is.finite(calibration_in_large)) %>%
    dplyr::arrange(horizon)
  if (!nrow(r)) return(NULL)
  bits <- sprintf("%.2f%s raw and %s recalibrated at *h*=%d", r$calibration_in_large, TIMES,
                  ifelse(is.finite(r$calibration_in_large_recal),
                         sprintf("%.2f%s", r$calibration_in_large_recal, TIMES), "not available"),
                  as.integer(r$horizon))
  sprintf(paste0("Over all leave-future-out back-test rounds, the featured model **%s** has ",
                 "calibration-in-the-large (mean predicted probability / observed invasion rate; ",
                 "above 1 = over-prediction) of %s."),
          featured, paste(bits, collapse = "; "))
})

# --- priority: §5 vulnerability-adjusted preparedness priority list ----------------
blocks$priority <- local({
  if (is.null(risk_sc) || !"priority" %in% names(risk_sc)) return(NULL)
  d <- risk_sc %>% dplyr::filter(horizon == 1, !is.na(priority)) %>%
    dplyr::arrange(dplyr::desc(priority)) %>% head(5)
  if (!nrow(d)) return(NULL)
  lst <- paste(sprintf("**%s** (%s)", d$health_zone, f2(d$priority)), collapse = ", ")
  sprintf("Top priorities are %s.", lst)
})

# --- onset_imputation: §1 share + mechanism, read from run_info.json --------
# This sentence used to be hand-written ("~15 % ... onset = sample_date - Delta"). Both halves
# went stale: the share moved, and the interval is now drawn from same-class, same-week DONORS
# (run_info's `mechanism`), not from a fitted delay. Deriving it from the run's own metadata is
# the only way it cannot drift from the code again.
blocks$onset_imputation <- local({
  ri <- file.path(O, "key_outputs", "run_info.json")
  if (!file.exists(ri)) return(NULL)
  .ri <- tryCatch(jsonlite::fromJSON(ri, simplifyVector = TRUE), error = function(e) NULL)
  oi <- .ri$onset_imputation
  .dmin <- suppressWarnings(as.integer(.ri$settings$onset_donor_min))
  if (is.null(oi) || is.null(oi$mechanism) || is.na(oi$mechanism)) return(NULL)
  pct  <- suppressWarnings(as.numeric(oi$confirmed_pct))
  lpct <- suppressWarnings(as.numeric(oi$linelist_pct))
  nsit <- suppressWarnings(as.numeric(oi$n_sitrep_rows))
  share <- if (is.finite(pct)) sprintf("The %.0f %% of confirmed records that lack", pct)
           else "Confirmed records that lack"
  # Say what that percentage is made of. It is NOT all DHIS2 reporting incompleteness: the
  # sitrep-reconciliation rows are a count reconciliation with no onset by construction, so
  # quoting only the combined figure overstates line-list missingness.
  split_sentence <- if (is.finite(pct) && is.finite(lpct) && is.finite(nsit) && nsit > 0)
    sprintf(paste0(" That figure combines **%.0f %%** of DHIS2 line-list records with no usable ",
                   "recorded onset with %s sitrep-reconciliation records, which carry no onset date by ",
                   "construction and are therefore imputed in full."),
            lpct, format(nsit, big.mark = ",")) else ""
  # The REPORTED onset->sample delay, for reference only: the imputation draws donor intervals
  # and does not use this fit (it used to be called "the deployed fit" here, as if it did). Its
  # mean and family are read from the fit summary the pipeline actually wrote, not asserted.
  fit_sentence <- local({
    fp <- file.path(O, "diagnostics", "delay_fits", "dhis2_delay_fit_summary.csv")
    if (!file.exists(fp)) return("")
    f <- tryCatch(readr::read_csv(fp, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(f) || !nrow(f)) return("")
    r <- f %>% dplyr::filter(delay == "onset_sample", estimator == "epidist_marginal",
                             plotted %in% TRUE, is.finite(mean_d))
    if (!nrow(r)) return("")
    r <- r[1, ]
    ci <- if (is.finite(r$mean_lo) && is.finite(r$mean_hi))
            sprintf(" [95 %% CrI %.1f, %.1f]", r$mean_lo, r$mean_hi) else ""
    sprintf(paste0(" For reference, the fitted onset-to-sample delay (EpiDist marginal, ",
                   "truncation-corrected, all classifications) is a **%s** with mean **%.1f d**%s ",
                   "(n = %s complete pairs); the imputation does not use it."),
            r$family, r$mean_d, ci, format(r$n, big.mark = ","))
  })
  # The mechanism IS the definition of Delta (the loader's own description: the sample date
  # minus an interval drawn from same-class, same-week donors), so it is quoted as the onset,
  # with the donor-window widening the description leaves out.
  widen <- if (length(.dmin) == 1L && is.finite(.dmin))
    sprintf(", widening the donor window by one week on each side at a time until it holds at least %d donors", .dmin)
  else ""
  sprintf(paste0("  %s a usable recorded onset receive an imputed onset: %s%s (`onset = sample_date - Delta`, ",
                 "Delta drawn **per record**).%s%s"),
          share, oi$mechanism, widen, fit_sentence, split_sentence)
})

# ---------------------------------------------------------------------------
# Apply
# ---------------------------------------------------------------------------
if (!file.exists(REPORT_PATH)) stop("[update_report] report not found: ", REPORT_PATH, call. = FALSE)
orig <- paste(readLines(REPORT_PATH, warn = FALSE), collapse = "\n")
# Cross-check generators against the report's own markers, both directions: a generator with no
# marker would never be injected, and a marker with no generator would keep stale text forever.
BLOCK_KEYS <- unique(stringr::str_match_all(orig, "<!-- AUTOGEN:([A-Za-z0-9_]+) -->")[[1]][, 2])
# (A marker with no generator is the "stale" case; the loop below warns about each one.)
.orphan_gen <- setdiff(names(blocks), BLOCK_KEYS)
if (length(.orphan_gen))
  stop("[update_report] generator(s) with no matching marker in the report (their output would ",
       "be discarded): ", paste(.orphan_gen, collapse = ", "), call. = FALSE)
txt  <- orig
changed <- character(0); stale <- character(0)
for (key in BLOCK_KEYS) {
  content <- if (key %in% names(blocks)) blocks[[key]] else NULL
  if (is.null(content)) {
    message(sprintf("[update_report] %-16s : no data — marker left at its PREVIOUS content", key))
    stale <- c(stale, key); next
  }
  new <- inject(txt, key, content)
  # detect whether this block's region actually changed (compare the whole string is
  # coarse; instead re-extract the region after injection isn't trivial — compare texts)
  if (!identical(new, txt)) changed <- c(changed, key)
  txt <- new
}

# Announce every block that could NOT be regenerated. Each one leaves its marker holding the
# previously-rendered text, which is indistinguishable from current text to any reader of the
# published report — so it must be loud here.
if (length(stale))
  warning(sprintf("[update_report] %d block(s) had no data and keep their PREVIOUS content: %s",
                  length(stale), paste(stale, collapse = ", ")), call. = FALSE, immediate. = TRUE)

if (identical(txt, orig)) {
  message("[update_report] report already up to date (no marked block changed).")
  quit(status = 0)
}

if (DRY_RUN) {
  message(sprintf("[update_report] --check: %d block(s) would change: %s",
                  length(changed), paste(changed, collapse = ", ")))
  quit(status = if (length(changed)) 1L else 0L)
}

# Timestamped backup, then write.
bak <- paste0(REPORT_PATH, ".bak")
writeLines(orig, bak)
writeLines(txt, REPORT_PATH)
message(sprintf("[update_report] updated %d block(s): %s", length(changed), paste(changed, collapse = ", ")))
message(sprintf("[update_report] featured model: %s | wrote %s (backup: %s)",
                featured, basename(REPORT_PATH), basename(bak)))
if (length(stale))
  message(sprintf("[update_report] NOT regenerated (previous content retained): %s",
                  paste(stale, collapse = ", ")))

}   # end of the script-run guard

