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
#   16_invasion_eval.R                            best_invasion_model() (featured pick)
#
# The featured single Bayesian model is chosen by the SAME calibration-aware CV
# composite the pipeline uses (best_invasion_model), and the covariate table is drawn
# from the best AVAILABLE covariate model — so the report adapts when the config
# toggles change which models are fitted (e.g. the "geo" model being off by default).
#
# USAGE:   Rscript spatiotemporal/update_bayesian_report.R
#          Rscript spatiotemporal/update_bayesian_report.R --check   # dry-run: report
#                  which blocks would change; write nothing; non-zero exit if any differ
# =============================================================================

suppressWarnings(suppressMessages({
  library(here); library(dplyr); library(readr); library(stringr); library(jsonlite)
}))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 ||
                            (length(a) == 1 && is.na(a))) b else a

ST_DIR <- file.path(here::here(), "spatiotemporal")
# Load config (paths, ANALYSIS_DATE) and the evaluation helpers
# (best_invasion_model) — the SAME selector the pipeline features with, so this report
# can never disagree with run_all.R on which model is featured.
suppressWarnings(suppressMessages({
  source(file.path(ST_DIR, "00_config.R"))
  source(file.path(ST_DIR, "16_invasion_eval.R"))
}))

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
# Featured single Bayesian model, by the pipeline's calibration-aware CV composite
# (exclude the ensembles so the featured model always has a current forecast).
featured <- NULL
if (!is.null(eval_tbl) && "method" %in% names(eval_tbl)) {
  bayes_singles <- eval_tbl %>% dplyr::filter(grepl("^Bayes", method), !grepl("-ens-", method))
  featured <- tryCatch(best_invasion_model(if (nrow(bayes_singles)) bayes_singles else eval_tbl),
                       error = function(e) NA_character_)
}
if (is.null(featured) || is.na(featured))
  featured <- if (!is.null(risk_sc) && "method" %in% names(risk_sc)) risk_sc$method[1] else "(featured model)"

# h=1 evaluation support (from the featured model's row, or the first full-CV row)
ev1 <- if (!is.null(eval_tbl)) eval_tbl %>% dplyr::filter(horizon == 1) else NULL
ev2 <- if (!is.null(eval_tbl)) eval_tbl %>% dplyr::filter(horizon == 2) else NULL
supp_row <- function(ev) {
  if (is.null(ev) || !nrow(ev)) return(NULL)
  e <- if ("partial_cv" %in% names(ev)) ev %>% dplyr::filter(!(partial_cv %in% TRUE)) else ev
  if (!nrow(e)) e <- ev
  e[1, ]
}
s1 <- supp_row(ev1); s2 <- supp_row(ev2)

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
          ad, if (is.na(.cut)) ad else .cut, at_risk, n_zones, conf_txt, affected)
})

# --- sample_events: §1 sample-size caveat (events + base rates) --------------------
blocks$sample_events <- local({
  if (is.null(s1)) return(NULL)
  br1 <- 100 * (s1$base_rate %||% NA); br2 <- if (!is.null(s2)) 100 * (s2$base_rate %||% NA) else NA
  sprintf(paste0("Pooled over the leave-future-out folds there are only ~%d realised invasion events ",
                 "at *h*=1 from a handful of distinct zones, against a per-zone-week base rate of ",
                 "~%s %% at *h*=1 (~%s %% at *h*=2). Every ranking is therefore provisional and ",
                 "differences between good models are frequently within Monte-Carlo noise."),
          as.integer(s1$n_invasions %||% NA), f2(br1),
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
  e <- ev1 %>% dplyr::filter(grepl("^Bayes", method), !grepl("-ens-", method))
  if ("partial_cv" %in% names(e)) e <- e %>% dplyr::filter(!(partial_cv %in% TRUE))
  if (!nrow(e)) return(NULL)
  e <- e %>% dplyr::arrange(dplyr::desc(auc_pr_skill)) %>% head(6)
  have_ci <- all(c("auc_pr_lo", "auc_pr_hi") %in% names(e))
  rows <- vapply(seq_len(nrow(e)), function(i) {
    skill <- if (have_ci && is.finite(e$auc_pr_lo[i]) && is.finite(e$base_rate[i]))
      sprintf("%s [%s%s%s]", f1(e$auc_pr_skill[i]),
              f1(pmax(e$auc_pr_lo[i], 0) / e$base_rate[i]), DASH, f1(e$auc_pr_hi[i] / e$base_rate[i]))
      else f1(e$auc_pr_skill[i])
    name <- if (identical(e$method[i], featured))
      sprintf("**%s** (featured, best CV composite)", e$method[i]) else e$method[i]
    sprintf("| %s | %s | %s | %s | %s%s |", name, skill,
            f1(e$mean_rank_of_truth[i]), f3(e$log_score[i]), f1(e$calibration_in_large[i]), TIMES)
  }, character(1))
  paste(c("| Model | AUC-PR skill [90 % CI] | Rank of truth | Log-score | Calibration |",
          "|---|---|---|---|---|", rows), collapse = "\n")
})

# --- featured_pick: §5 sentence naming the featured model + its headline stats ------
blocks$featured_pick <- local({
  if (is.null(ev1)) return(NULL)
  fr <- ev1 %>% dplyr::filter(method == featured)
  fr2 <- if (!is.null(ev2)) ev2 %>% dplyr::filter(method == featured) else NULL
  if (!nrow(fr)) return(NULL)
  skill <- f1(fr$auc_pr_skill[1]); calib <- f1(fr$calibration_in_large[1])
  # CHECK the superlative instead of asserting it. "the tightest calibration among the leaders"
  # was hard-coded prose: on the shipped run Bayes-M16-fill-med was tighter (1.991 vs 2.033),
  # and the leaderboard printed immediately above rounds both to 2.0x, so the table visibly tied
  # while the prose claimed a win. Calibration-in-the-large is best at 1, so "tightest" means
  # smallest |x - 1| among the same six models the leaderboard block lists.
  calib_word <- local({
    lead <- ev1 %>% dplyr::filter(grepl("^Bayes", method), !grepl("-ens-", method))
    if ("partial_cv" %in% names(lead)) lead <- lead %>% dplyr::filter(!(partial_cv %in% TRUE))
    lead <- lead %>% dplyr::arrange(dplyr::desc(auc_pr_skill)) %>% utils::head(6)
    d <- abs(lead$calibration_in_large - 1)
    if (!nrow(lead) || !any(is.finite(d))) "competitive"
    else if (identical(lead$method[which.min(d)], featured)) "the tightest"
    else "competitive"
  })
  rank2 <- if (!is.null(fr2) && nrow(fr2)) f1(fr2$mean_rank_of_truth[1]) else "?"
  # NEVER fabricate the number. The fallback used to be the literal "0.99", so a run whose
  # evaluation table lacked auc_roc would publish a discrimination figure nothing computed
  # (and one that rounds the real value, ~0.98, upward).
  # SAME FIELD AS THE LEADERBOARD. max() ran over EVERY h=1 row — the three structural
  # baselines and any partial_cv rows included — so "best AUC-ROC" could be a rank-only
  # comparator's number published as the Bayesian suite's discrimination. It lands on a
  # Bayesian model on the current frame by luck, not by construction.
  .ev1_elig <- ev1[grepl("^Bayes", ev1$method) &
                     !(ev1$partial_cv %in% TRUE), , drop = FALSE]
  if (!nrow(.ev1_elig)) .ev1_elig <- ev1
  auroc <- if ("auc_roc" %in% names(.ev1_elig) && any(is.finite(.ev1_elig$auc_roc)))
             f2(max(.ev1_elig$auc_roc, na.rm = TRUE)) else "not computed in this run"
  sprintf(paste0("The featured model is chosen by the **CV composite** (summed within-horizon ranks of ",
    "AUC-PR skill + mean rank-of-truth + log-score, pooled over both horizons): **%s** wins it, pairing ",
    "solid discrimination (*h*=1 AUC-PR skill %s%s) with %s calibration (%s%s) ",
    "and a strong mean rank of truth at *h*=2 (%s). Best AUC-ROC %s %s. The **loo stacking weights** are ",
    "used only to build the loo-stacked ENSEMBLE (`bayes_ensemble_*`), NOT to pick the featured single model."),
    featured, skill, TIMES, calib_word, calib, TIMES, rank2, "≈", auroc)
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
  cov_models <- params %>% dplyr::filter(!is_intercept) %>% dplyr::distinct(model) %>% dplyr::pull(model)
  if (!length(cov_models)) return(NULL)
  # Prefer the featured model if it carries covariates; else the "full" exogenous model;
  # else the first covariate model available.
  pick <- if (featured %in% cov_models) featured
          else if (any(grepl("-full$", cov_models))) grep("-full$", cov_models, value = TRUE)[1]
          else cov_models[1]
  d <- params %>% dplyr::filter(model == pick, !is_intercept)
  # Fixed, human-authored READINGS per covariate term; numbers are filled from the fit.
  reading <- c(
    d_min = "Farther from the frontier ⇒ lower risk (the strongest, most interpretable driver)",
    ccvi  = "More deprived ⇒ higher risk",
    log_pop = "Negative *residual* conditional on the gravity offset — a density-vs-frequency correction, not “big cities are safer”",
    healthsite_density = "Health-facility density — not credibly different from 1 under these few events")
  label <- c(d_min = "`d_min` (travel time to nearest affected zone)", ccvi = "`ccvi` (deprivation)",
             log_pop = "`log_pop`", healthsite_density = "`healthsite_density`")
  ord <- c("d_min", "ccvi", "log_pop", "healthsite_density")
  d <- d %>% dplyr::arrange(match(term, ord))
  rows <- vapply(seq_len(nrow(d)), function(i) {
    t <- d$term[i]
    sprintf("| %s | %s %s | %s | %s |", label[t] %||% sprintf("`%s`", t),
            f2(d$hr[i]), cri(d$lo[i], d$hi[i], 2), f2(d$p_dir[i]), reading[t] %||% "")
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
    sprintf("| %s | %s | %s %s | %s%s |", d$health_zone[i], d$province[i] %||% "",
            f3(d$p_case_invasion[i]), cri(d$p_lo[i], d$p_hi[i], 3),
            f1(d$rr_nat[i]), TIMES), character(1))
  cap <- sprintf(paste0("**Highest-risk at-risk zones (%s, next week, posterior mean [90 %% CrI]).** ",
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
  fam  <- c(M4 = "gravity", M8 = "short-trip + gravity composite",
            M9 = "multi-kernel ensemble", M10 = "short-trip + radiation composite",
            M11 = "inward meeting-location FOI", M13 = "cohort + gravity composite",
            M14 = "cohort + radiation composite", M15 = "symmetrised relocation OD",
            M16 = "cohort + relocation OD composite",
            M17 = "all-kernel consensus ensemble")[base]
  fam  <- if (is.na(fam)) base else fam
  gt   <- if (grepl("-short", lbl)) "short" else if (grepl("-long", lbl)) "long" else "medium"
  cov  <- if (grepl("full-susp", lbl)) "full exogenous + suspected-case covariates"
          else if (grepl("-susp", lbl)) "suspected-case covariates"
          else if (grepl("-full", lbl)) "full exogenous covariates"
          else if (grepl("-geo", lbl))  "geo covariates"
          else "no covariates"
  qual <- paste0(
    if (!is.na(kern) && grepl("-dist", kern)) " on road distance" else "",
    if (!is.na(kern) && grepl("-split", kern)) ", cohort rows split per origin" else "",
    if (!is.na(kern) && grepl("-(fill|split)", kern)) ", with source-cell fill" else "")
  sprintf("the %s kernel%s with the %s generation-time profile, %s", fam, qual, gt, cov)
}
blocks$featured_name <- local({
  if (is.null(featured) || is.na(featured)) return(NULL)
  sprintf(paste0("In the current run this is **%s** (%s). Ensembles are excluded from this pick so ",
                 "the featured single model always has a current forecast."),
          featured, describe_bayes(featured))
})

# --- targeting: §5 Figure 2B top-10 catch-rate vs random (operational lift) --------
blocks$targeting <- local({
  if (is.null(ev1)) return(NULL)
  fr <- ev1 %>% dplyr::filter(method == featured)
  if (!nrow(fr) || !all(c("recall_at_10", "n_atrisk", "n_folds") %in% names(fr))) return(NULL)
  rec <- fr$recall_at_10[1]
  per_fold <- (fr$n_atrisk[1] %||% NA) / (fr$n_folds[1] %||% NA)   # at-risk zones per fold
  rnd <- if (is.finite(per_fold) && per_fold > 0) min(10 / per_fold, 1) else NA
  if (!is.finite(rec) || !is.finite(rnd) || rnd <= 0) return(NULL)
  sprintf(paste0("**Operationally useful targeting (the prioritisation panel, Figure 2B of `key_outputs/figures/Figure2_labelled`).** A top-10 watch-list catches ~%.0f %% ",
                 "of the next invasions per round (mean over folds; the panel's pooled curve reads slightly higher) versus ~%.0f %% for a random list of the same size (a ~%.0f%s lift)."),
          100 * rec, 100 * rnd, rec / rnd, TIMES)
})

# --- calibration: §5 calibration-in-the-large for the FEATURED model ---------------
# Replaces a hand-written "~2-3x (2.6x for Bayes-M10-med)": both the magnitude and the model
# name were literals that no longer tracked the featured pick or the current folds.
# cal_in_large = mean(predicted p) / observed base rate over the leave-future-out rows, so
# >1 is over-prediction and the ratio is read directly as the over-prediction factor.
blocks$calibration <- local({
  fp <- file.path(O, "diagnostics", "invasion_recalibration.csv")
  if (is.null(featured) || is.na(featured) || !file.exists(fp)) return(NULL)
  rc <- tryCatch(readr::read_csv(fp, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(rc) || !all(c("method","horizon","cal_in_large") %in% names(rc))) return(NULL)
  r <- rc %>% dplyr::filter(method == featured, is.finite(cal_in_large)) %>%
    dplyr::arrange(horizon)
  if (!nrow(r)) return(NULL)
  bits <- sprintf("%.1f%s at *h*=%d", r$cal_in_large, TIMES, as.integer(r$horizon))
  sprintf(paste0("Over the leave-future-out folds the featured model **%s** over-predicts the ",
                 "absolute invasion probability by %s (calibration-in-the-large = mean predicted ",
                 "probability / observed base rate)."),
          featured, paste(bits, collapse = " and "))
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
# went stale: the realised share is ~24 %, and the delay is now drawn from the shared
# truncation-corrected EpiDist fit rather than a fixed shift. Deriving it from the run's own
# metadata is the only way it cannot drift from the code again.
blocks$onset_imputation <- local({
  ri <- file.path(O, "key_outputs", "run_info.json")
  if (!file.exists(ri)) return(NULL)
  oi <- tryCatch(jsonlite::fromJSON(ri, simplifyVector = TRUE)$onset_imputation,
                 error = function(e) NULL)
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
    sprintf(paste0(" That figure combines **%.0f %%** genuine missingness in the DHIS2 line ",
                   "list with %s sitrep-reconciliation records, which carry no onset date by ",
                   "construction and are therefore imputed in full."),
            lpct, format(nsit, big.mark = ",")) else ""
  # The fitted delay's mean and family are read from the fit summary the pipeline actually
  # wrote, not asserted: the report previously claimed "mean ~5.9 d" where the deployed,
  # truncation-corrected fit is materially longer.
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
            sprintf(" [%.1f, %.1f]", r$mean_lo, r$mean_hi) else ""
    sprintf(paste0(" The deployed fit is a **%s** with a corrected mean of **%.1f d**%s ",
                   "(n = %s complete pairs)."),
            r$family, r$mean_d, ci, format(r$n, big.mark = ","))
  })
  sprintf(paste0("  %s an onset date receive an imputed onset `onset = sample_date - Delta`, ",
                 "where the pipeline draws Delta **per record** via %s.%s%s"),
          share, oi$mechanism, fit_sentence, split_sentence)
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
