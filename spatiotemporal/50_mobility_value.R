# =============================================================================
# 50_mobility_value.R
# BDBV 2026 DRC - Spatiotemporal invasion forecasting
#
# WHAT DOES THE MOBILE-PHONE MOBILITY DATA ADD, over and above population and geography?
#
# WHY THIS MODULE EXISTS. The model grid cannot answer that question on its own. Every kernel
# in it uses Flowminder data somewhere: the seven composites carry measured cohort or
# short-trip rows directly, and M4/M4c use the phone data to fit their gravity deterrence. The
# comparators that use none - Gravity-B4, Distance-B1, Adjacency-B7, and the case-blind road-
# travel-time rule -
# differ from the featured model in more than one way at a time (case-blind, or unfitted
# exponents, or not a renewal model at all), so a comparison against them confounds "no phone
# data" with "different model".
#
# WHAT IT ADDS. Two kernels 03_mobility_matrices.R already builds but bayes_default_grid()
# never fits as a model's own kernel, and which use no mobile-network data at any point:
#   M5   radiation (Simini et al. 2012) - WorldPop population + OSRM travel time, nothing fitted
#   M6a  exponential decay in OSRM road travel time - geography alone
# Fitting the SAME Bayesian renewal model on these, at the same generation time, covariates,
# folds, nowcast and shared support, isolates the phone data as the only thing that differs.
#
# THE NESTED ABLATIONS ARE THE STRONGEST RESULT HERE. Each composite kernel IS its background
# kernel with measured movement overlaid on the origins the phone data could observe (the
# source-cell fill supplies the rest from that same background). So M14-fill against M5 is
# literally "add the measured cohort rows to the radiation kernel and change nothing else", and
# M8-fill against M4 is "add the measured epicentre short trips to the fitted gravity kernel".
#
# THESE ARE ABLATIONS, NOT CANDIDATES, and the mechanism is structural rather than a flag: they
# are not in bayes_default_grid(), their rows are written to their own RDS and are never merged
# into the pipeline's lfo_results.rds, so no selection path can reach them. On the runs seen so
# far including them would not have changed the pick anyway (the best phone-free kernel scores
# below the featured model at both horizons). The module
# ASSERTS that appending them leaves every published method's average precision unchanged - see
# the support note below.
#
# WHY IT APPENDS INSTEAD OF RE-EVALUATING. Every model is scored on one shared support, the
# (origin x zone) cells every method covers. Re-running evaluate_invasion() over a DIFFERENT set
# of methods would move that support and silently change every published number. So the arms'
# LFO rows are appended to the published lfo_results.rds and everything is scored on the union.
# Every method covers every at-risk cell (a sparse kernel gives a zero import force, i.e. a
# finite p, not a missing one), and mobval_assert_support_unchanged() fails loudly if appending
# the arms ever moves the support.
#
# WHY THE BOOTSTRAP IS PAIRED. All methods are scored on the same zones and the same origins, so
# their errors are strongly correlated and the marginal intervals are much wider than the
# interval of the difference. Reading whether two marginal intervals overlap would understate
# the evidence for a difference. Zones are resampled once per horizon and BOTH methods are
# recomputed on each resample, so every contrast is paired and all contrasts share resamples.
#
# TWO METRICS, NOT ONE, BOTH PER ROUND. AUC-PR skill PER ROUND (average precision within each
# forecast round, event-weighted, over the base rate: the estimator the selection uses)
# answers "how well does the weekly list put that week's invaded zones at the top"; the mean
# rank of an invaded zone answers "how far down the weekly watch-list did the team have to read
# to find it". Both are scored over ALL rounds, the basis on which performance is reported; the
# selection reads the post-burn-in rounds by default, which is the selection's own deliberate
# choice (INVASION_SELECTION_BASIS, 16_invasion_eval.R) and not a reporting basis. The pooled AP
# (every round's rows ranked in one list) is NOT used: it ranks every round's rows in one list and so rewards a
# model's probability LEVEL on the rounds that had invasions (16_invasion_eval.R,
# .auc_pr_by_date()); a case-blind fixed ordering scored first pooled and twentieth per round.
# A conclusion about what the phone data adds should not rest on whichever metric happens to be
# quoted, so both are drawn, at both horizons. The skill's reference is its CHANCE level (a
# random ordering's expected per-round skill, `chance`), not 1.
#
# Usage:   Rscript spatiotemporal/50_mobility_value.R
#          (safe to source(); the work is guarded by is_script_run())
#
# Reads:   outputs/forecasts/lfo_results.rds   (the published pooled LFO, every method)
# Writes:  outputs/key_outputs/mobility_value_ladder.csv   per-method skill (+ chance), rank, CIs,
#                                                 recall@15, all over every round
#          outputs/key_outputs/mobility_value_contrasts.csv paired contrasts
#          outputs/forecasts/mobility_value_arms.rds        the ablation arms' LFO rows (cached)
#          outputs/key_outputs/manuscript_figures/FigureS24_mobility_value.pdf/.png
#                                                 (AUC-PR skill, both horizons)
#          outputs/key_outputs/manuscript_figures/FigureS25_mobility_value_rank.pdf/.png
#                                                 (mean rank of truth, both horizons)
#          outputs/key_outputs/manuscript_figures/FigureS26_mobility_value_skill_and_rank.pdf/.png
#                                                 (both metrics x both horizons, four panels)
#
# ALSO READ BY THE MAIN TEXT. mobility_value_ladder.csv is what Figure 2 panel C is drawn from
# (make_manuscript_figures.R), through the shared panel builder in mobility_value_panel.R, at
# the 2-week horizon panels A and B are fixed to. So this module is no longer only an appendix:
# run_all.R runs it in the downstream block BEFORE the figure suites, and Figure 2 refuses to
# build without its CSV rather than quietly losing a panel. Two consequences worth knowing:
#   * the ladder scores the RAW p_invasion in both passes, as the selection does (per-round
#     metrics are invariant to the prequential recalibration anyway: it is monotone within a
#     round);
#   * Figure 2 CHECKS the featured model's skill in this CSV against the same quantity in
#     diagnostics/invasion_evaluation.csv (auc_pr_skill_date), so a ladder left over from an
#     earlier run cannot be drawn beside a current panel A.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(tidyr)
})

ST_DIR     <- here::here("spatiotemporal")
if (!exists("ANALYSIS_DATE")) source(file.path(ST_DIR, "00_config.R"))
stop_if_data_definition_sensitivity("50_mobility_value.R")   # hard-coded published paths
# fs_out_dir()/fs_caption(): the raw-scale sibling directory and the probability-scale stamp,
# the same two the manuscript suites publish through.
if (!exists("fs_out_dir"))   source(file.path(ST_DIR, "forecast_scale.R"))

OUT_KEY    <- file.path(ST_DIR, "outputs", "key_outputs")
OUT_FIG    <- fs_out_dir(file.path(OUT_KEY, "manuscript_figures"))
OUT_FCST   <- file.path(ST_DIR, "outputs", "forecasts")
ARMS_RDS   <- file.path(OUT_FCST, "mobility_value_arms.rds")

# theme_pub(), the Okabe-Ito palette and save_dual() - ONE definition of the house style,
# shared with make_manuscript_figures.R and the two SI suites, so these figures cannot drift
# from the ones they sit beside. save_dual() resolves its destination at call time from
# PANEL_DIR, which must therefore be set before it is called; every call here passes `dir`
# explicitly as well.
PANEL_DIR <- OUT_FIG
source(file.path(ST_DIR, "figure_style.R"))
# The ladder's DRAWING layer - the row set, the two encodings, the ordering rule and
# mobval_panel() itself. Shared with make_manuscript_figures.R, which draws the same ladder as
# Figure 2 panel C; see the header of that file for why it is not defined here.
source(file.path(ST_DIR, "mobility_value_panel.R"))

MOBVAL_KERNELS <- c("M6a", "M5")   # phone-free: geography only, then + population
# 10,000, not 2,000. The quoted quantities are the 5th and 95th percentiles of a bootstrap
# distribution, and at 2,000 replicates their SECOND DECIMAL is not stable: re-running with a
# different resample stream moved interval bounds by up to 0.03, which is larger than the
# precision the Article quotes them to. 10,000 holds them to +/-0.01.
MOBVAL_NBOOT   <- 10000L
MOBVAL_SEED    <- 20260925L

# ---------------------------------------------------------------------------
# Metrics. Both reimplement the pipeline's own estimators exactly, so the numbers here can sit
# beside those in invasion_evaluation.csv; the module checks that they do.
# ---------------------------------------------------------------------------

#' Tie-aware average precision. Tied scores cannot be separated by any threshold, so a tied run
#' collapses to ONE operating point; treating them separately would draw a diagonal through a
#' tied block and claim discrimination that is not there.
mv_ap <- function(y, p) {
  y <- as.integer(y); P <- sum(y)
  if (P == 0L) return(NA_real_)
  o <- order(p, decreasing = TRUE); y <- y[o]; ps <- p[o]
  tp <- cumsum(y); fp <- cumsum(1 - y)
  last <- c(ps[-length(ps)] != ps[-1], TRUE)
  tpl <- tp[last]; fpl <- fp[last]
  sum(diff(c(0, tpl / P)) * (tpl / (tpl + fpl)))
}

#' PER-ROUND AUC-PR skill, replicating .auc_pr_by_date() in 16_invasion_eval.R: average precision
#' within each forecast round (fold), combined across rounds with event weights, over the
#' event-weighted base rate; and the same quantity under a uniformly random ordering (Bestgen
#' 2015), the axis's chance level. This is the selection's AUC-PR estimator (the ladder applies it
#' over all rounds): the pooled AP rewards a model's probability LEVEL on the rounds that had
#' invasions, not its ordering within a round.
#' @return c(ap, skill, chance).
mv_ap_by_date <- function(y, p, date) {
  y <- as.integer(y)
  if (!length(y) || sum(y) == 0L) return(c(ap = NA_real_, skill = NA_real_, chance = NA_real_))
  s <- vapply(split(seq_along(y), date), function(i) {
    E <- sum(y[i]); if (E == 0) return(c(0, 0, 0, 0))
    n <- length(i); Hn <- sum(1 / seq_len(n))
    eap <- (Hn + if (n > 1L) (E - 1) / (n - 1) * (n - Hn) else 0) / n
    c(E * mv_ap(y[i], p[i]), E, E * mean(y[i]), E * eap)
  }, numeric(4))
  c(ap = sum(s[1, ]) / sum(s[2, ]), skill = sum(s[1, ]) / sum(s[3, ]),
    chance = sum(s[4, ]) / sum(s[3, ]))
}

#' Recall at a weekly watch-list of K, replicating .ranking_metrics() in 16_invasion_eval.R.
#' Three details matter: the list is drawn PER FORECAST ORIGIN; ties use ties.method = "max", so
#' a zone counts as inside the top K only if monitoring K zones necessarily includes it; and the
#' per-origin recalls are AVERAGED, not pooled, with origins carrying no invasion dropped.
#' Pooling the cells and using order() instead overstates recall by 2-4 points.
mv_recall_at <- function(df, K) {
  per <- vapply(split(df, df$fold_id), function(g) {
    g <- g[is.finite(g$p_invasion), ]
    y <- as.integer(g$is_new_invasion)
    if (!nrow(g) || sum(y) == 0L) return(NA_real_)
    sum(y[rank(-g$p_invasion, ties.method = "max") <= K]) / sum(y)
  }, numeric(1))
  mean(per, na.rm = TRUE)
}

#' Area under the ROC curve, tie-aware (ties share the average rank), matching .auc_roc().
mv_auc_roc <- function(y, p) {
  y <- as.integer(y); P <- sum(y); N <- length(y) - P
  if (P == 0L || N == 0L) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - P * (P + 1) / 2) / (P * N)
}

#' Mean rank of the truly invaded zones, per origin then averaged over origins - the estimator
#' .ranking_metrics() uses. Ties share the AVERAGE rank, so a truth tied at the bottom with
#' hundreds of zones is not handed an arbitrary position.
mv_mean_rank <- function(df) {
  per <- vapply(split(df, df$fold_id), function(g) {
    g <- g[is.finite(g$p_invasion), ]
    y <- as.integer(g$is_new_invasion)
    if (!nrow(g) || sum(y) == 0L) return(NA_real_)
    mean(rank(-g$p_invasion, ties.method = "average")[y == 1])
  }, numeric(1))
  mean(per, na.rm = TRUE)
}

#' THE SAME ESTIMATOR UNDER EVERY ZONE-CLUSTER RESAMPLE AT ONCE.
#'
#' WHY NOT A LOOP. mv_mean_rank() called inside the resample loop would re-rank every scored row
#' once per resample, method and horizon. It is also unnecessary: a resample changes nothing about the scores,
#' only the MULTIPLICITY m_z with which each zone enters, and the average rank of a truth t
#' within its own forecast origin is LINEAR in those multiplicities. For a truth t in origin f
#' with score p_t, summing over the zones z present in that origin,
#'
#'     rank_t = sum_z m_z 1{p_z > p_t} + ( sum_z m_z 1{p_z = p_t} + 1 ) / 2
#'            = ( a_t . m ) + 1/2,    a_t(z) = 1{p_z > p_t} + 1{p_z = p_t}/2
#'
#' so every resample's ranks are one matrix product A %*% M: A is (truths x zones), M is
#' (zones x resamples). The tie term is exactly rank(ties.method = "average") - which is why
#' the scores are compared with `==` and not with a tolerance: that is the comparison the
#' published estimator makes, and a tolerance here would silently score a different quantity.
#'
#' The rest follows mv_mean_rank() row for row: within an origin each truth is weighted by its
#' OWN multiplicity (a zone drawn twice contributes its rank twice), an origin whose truths
#' were all dropped from the resample contributes nothing rather than a zero, and the resample
#' statistic is the mean over the origins that survive - the na.rm mv_mean_rank() does.
#'
#' mv_rank_setup() holds the parts that do not depend on the method: which cells are truths,
#' which origin and which zone each cell belongs to, and the cells making up each origin.
mv_rank_setup <- function(fold_idx, zone_idx, y, nz) {
  pos <- which(y == 1L)
  if (!length(pos))
    stop("[mobval] no invasion in the support; the rank bootstrap is undefined.", call. = FALSE)
  nf <- max(fold_idx)
  list(pos = pos, zone = zone_idx, nz = nz,
       # factor(levels = seq_len(nf)) so the list is indexable BY ORIGIN NUMBER: split() on a
       # bare integer vector drops empty origins and the positional lookup below would then
       # silently read the wrong origin's cells.
       rows_of_fold = split(seq_along(fold_idx), factor(fold_idx, levels = seq_len(nf))),
       pos_fold = fold_idx[pos], pos_zone = zone_idx[pos])
}

#' @param M zones x resamples multiplicity matrix (M[z, b] = times zone z was drawn in b).
#'   Pass a single all-ones column to recover the point estimate; the module asserts that it
#'   does, against mv_mean_rank(), for every method and horizon.
mv_boot_mean_rank <- function(p, st, M) {
  A <- matrix(0, length(st$pos), st$nz)
  for (i in seq_along(st$pos)) {
    j  <- st$rows_of_fold[[ st$pos_fold[i] ]]
    p0 <- p[st$pos[i]]
    A[i, st$zone[j]] <- (p[j] > p0) + 0.5 * (p[j] == p0)
  }
  R <- A %*% M + 0.5                      # truths x resamples: the rank of each truth
  W <- M[st$pos_zone, , drop = FALSE]     # ... and how many times it appears
  num <- rowsum(R * W, st$pos_fold, reorder = TRUE)
  den <- rowsum(W,     st$pos_fold, reorder = TRUE)
  nf  <- colSums(den > 0)
  fm  <- num / den; fm[den == 0] <- 0     # an origin with no surviving truth drops out
  out <- colSums(fm) / nf
  out[nf == 0L] <- NA_real_               # ... and a resample with no truth at all is no draw
  out
}

#' Reduce to the shared support: the (origin x zone) cells every method covers, per horizon.
#' Scoring each method on its own reach inflates whichever reaches furthest back.
mv_shared_support <- function(lfo) {
  d <- lfo %>%
    filter(!(was_active_before %in% TRUE), is.finite(p_invasion), is.finite(is_new_invasion)) %>%
    mutate(cell = paste(fold_id, health_zone, sep = "\r"))
  n_meth <- dplyr::n_distinct(d$method)
  keep <- d %>% count(horizon, cell, name = "n_m") %>%
    group_by(horizon) %>% filter(n_m == max(n_m)) %>% ungroup()
  # `max(n_m)` is the intersection ONLY IF some cell is covered by every method. If none is, this
  # would silently fall back to the most-covered cells and define a support that no longer means
  # "every method scored here" - the comparison would then be between methods that did not all
  # see the same rows. Fail instead of computing something that looks right.
  if (any(keep$n_m != n_meth))
    stop("[mobval] no cell is covered by all ", n_meth, " methods at horizon(s) ",
         paste(unique(keep$horizon[keep$n_m != n_meth]), collapse = ", "),
         "; the shared support would not be an intersection.", call. = FALSE)
  out <- semi_join(d, keep, by = c("horizon", "cell"))
  # Every method must end up on identical cells, or nothing downstream is paired.
  chk <- out %>% count(method, horizon) %>% group_by(horizon) %>%
    summarise(k = dplyr::n_distinct(n), .groups = "drop")
  if (any(chk$k != 1L))
    stop("[mobval] methods do not sit on identical cells after the support reduction.",
         call. = FALSE)
  out
}

# ---------------------------------------------------------------------------
# Step 1 - fit the ablation arms (cached; ~1 minute for four specs)
# ---------------------------------------------------------------------------
#' The settings that change what the arms ARE. A cache keyed on nothing would be reused after a
#' config edit and the contrasts would then compare freshly-fitted published models against
#' stale arms. Cell coverage and outcomes are additionally verified downstream against the
#' published run, so between the two a stale cache cannot pass silently.
#'
#' THE KEY USED TO MISS MOST OF WHAT MOVES A CROSS-VALIDATION: the h >= 2 source projection,
#' the R(t) cache version, the line-list snapshot, the as-of fold truncation, the rolling-
#' predictor floor, the fold-admission age, the sample-date rule, the onset stratum and mode,
#' and the kernel configuration. Outcome and support checks downstream catch some of those
#' changes (a different onset completion moves outcomes) but not all of them (a different fold
#' truncation or projection moves only the predictions), so they are keyed in two ways:
#'   * the PUBLISHED LFO's own stamp (`pub_stamp`, run_all.R's .lfo_stamp on lfo_results.rds),
#'     minus the entries naming the OTHER models (the grid labels and the baselines' definitions);
#'     the entries hashing the grid's specs and kernels are kept, which can only force a needless
#'     refit (the safe direction). The arms are appended to exactly that cross-validation and must have been fitted
#'     for it — this covers the code-level settings (BAYES_H2_PROJECTION is defined in 21, which
#'     is not loaded when this key is first computed) and the FITTED as-of truncation;
#'   * the config-level settings of THIS session, listed explicitly, because the arms are fitted
#'     in this session and a config edit not yet propagated to the published run must not be
#'     served arms fitted under the old value.
mobval_arms_key <- function(pub_stamp = NULL) {
  .cfg <- function(nms) lapply(stats::setNames(nm = nms), function(nm) get0(nm, ifnotfound = NA))
  # list_linelist_snapshots() (below) is defined in 01_data_prep.R, which this module otherwise
  # loads only on a cache MISS (load_pipeline_modules() in mobval_fit_arms()). run_all.R launches
  # 50 as its own Rscript with only 00_config.R loaded, so without this the key itself failed
  # ("could not find function") before the cache could even be consulted.
  if (!exists("list_linelist_snapshots", mode = "function"))
    source(file.path(ST_DIR, "01_data_prep.R"))
  k <- list(
  kernels      = MOBVAL_KERNELS,
  analysis_date = as.character(ANALYSIS_DATE),
  random_seed  = get0("RANDOM_SEED",           ifnotfound = NA),
  onset_impute_seed = get0("ONSET_IMPUTE_SEED", ifnotfound = NA),   # the onset draw's only seed
  geo_cov      = get0("BAYES_GEO_COVARIATES",  ifnotfound = c("ccvi", "d_min")),
  horizons     = LFO_HORIZONS,
  rolling      = get0("ROLLING_PREDICTORS",    ifnotfound = NA),
  fill         = get0("MOBILITY_SOURCE_FILL",  ifnotfound = NA),
  min_train_wk = get0("LFO_MIN_TRAINING_WEEKS", ifnotfound = NA),
  iter         = 600L,
  # Where each round's training data come from. ALWAYS stamped: the default changed from
  # re-censoring to the archived snapshots on 2026-10-15, so a cache without it was built under
  # the other regime and must not be reused.
  training_source = get0("LFO_TRAINING_SOURCE", ifnotfound = "vintage"),
  snapshots = if (identical(get0("LFO_TRAINING_SOURCE", ifnotfound = "vintage"), "vintage"))
    rlang::hash(list(list_linelist_snapshots()$folder, get0("TRUNC_SERIES_VERSION", ifnotfound = NA))))
  # The other data-definition switches (00_config.R), added only when off their defaults.
  k$data_definition <- nondefault_data_definition()
  k$lfo_config <- .cfg(c("RT_CACHE_VERSION", "ROLLING_PREDICTOR_FLOOR", "LFO_MIN_EVAL_AGE_DAYS",
                         "SAMPLE_DATE_SOURCE", "ONSET_MODE", "ONSET_IMPUTE_STRATUM",
                         "ONSET_DONOR_MIN", "ASCERTAINMENT_BURNIN_WEEKS",
                         "LFO_VINTAGE_MAX_GAP_DAYS", "TRUNC_SERIES_VERSION", "TRUNC_WINDOW_DAYS",
                         "TRUNC_THIN_DAYS", "TRUNC_MIN_SNAPSHOTS", "TRUNC_MAX_DAYS",
                         "OUTBREAK_START", "WEEK_ANCHOR"))
  # What the arms' kernels are built from. M5 (radiation) and M6a (travel-time decay) read only
  # the population and the OSRM travel-time matrix -- load_osrm() applies OSRM_GAP_FILL itself --
  # so their CONTENT is keyed, together with the code that builds the kernels. (No Flowminder
  # input or fill switch reaches either kernel, so none is keyed.)
  if (!exists("load_osrm", mode = "function"))
    source(file.path(ST_DIR, "03_mobility_matrices.R"))
  k$kernel_inputs <- list(osrm = rlang::hash(load_osrm()), pop = rlang::hash(load_population()),
                          code = unname(tools::md5sum(file.path(ST_DIR, "03_mobility_matrices.R"))))
  k$linelist_snapshot <- tryCatch({
    .lj <- jsonlite::fromJSON(LINELIST_JSON, simplifyVector = TRUE)
    paste(.lj$folder, .lj$processed_at, sep = "|")
  }, error = function(e) NA_character_)
  k$published_lfo_stamp <- if (is.null(pub_stamp)) NULL else
    pub_stamp[setdiff(names(pub_stamp),
                      c("bayes_grid", "baseline_models", "baseline_bodies", "invasion_models"))]
  k
}

#' Guard against silent divergence from the pipeline. The arms are only comparable to the
#' published models if they are cross-validated the way run_all.R cross-validates. If that call
#' site changes, this fails instead of quietly reporting numbers from a different
#' cross-validation.
mobval_assert_matches_run_all <- function() {
  src <- readLines(file.path(ST_DIR, "run_all.R"), warn = FALSE)
  i <- grep("^  run_invasion_lfo\\(zone_week_outbreak", src)
  if (length(i) != 1L)
    stop("[mobval] could not locate run_all.R's LFO call site; the argument guard did not run.",
         call. = FALSE)
  blk <- paste(src[i:(i + 20L)], collapse = " ")
  need <- c("delay = .trunc_asof", "horizons = LFO_HORIZONS", "analysis_date = ANALYSIS_DATE",
            "min_atrisk_events = 0L", "nowcast_fn = apply_nowcast_correction")
  miss <- need[!vapply(need, function(x) grepl(x, blk, fixed = TRUE), logical(1))]
  if (length(miss))
    stop("[mobval] run_all.R's LFO call no longer passes: ", paste(miss, collapse = "; "),
         ". Re-sync this module before trusting the comparison.", call. = FALSE)
  # The Bayesian models must also still be built with the same sampling effort. A missing call
  # site is an error, not a skip: the check used to look for a call run_all.R no longer makes,
  # and passed silently.
  j <- grep("bayes_lfo_models[[sp$label]] <- make_suite_lfo_model(sp,", src, fixed = TRUE)
  if (length(j) != 1L)
    stop("[mobval] could not locate run_all.R's back-test model construction; the sampling-",
         "effort guard cannot run. Re-sync this module.", call. = FALSE)
  blk2 <- paste(src[j:(j + 4L)], collapse = " ")
  # The effort is a named constant in run_all.R (BAYES_LFO_ITER_BRMS, stamped into the LFO cache);
  # the call passes it by name. Either form must resolve to 600, the effort the arms use.
  .iter_ok <- grepl("iter_brms = 600L", blk2, fixed = TRUE) ||
    (grepl("iter_brms = BAYES_LFO_ITER_BRMS", blk2, fixed = TRUE) &&
       any(grepl("^BAYES_LFO_ITER_BRMS <- 600L", src)))
  if (!.iter_ok)
    stop("[mobval] run_all.R no longer fits the grid at iter = 600; the arms would differ ",
         "from the published models in sampling effort.", call. = FALSE)
  invisible(TRUE)
}

mobval_fit_arms <- function(force = FALSE, pub_stamp = NULL) {
  mobval_assert_matches_run_all()
  key <- mobval_arms_key(pub_stamp)
  if (!force && file.exists(ARMS_RDS)) {
    cached <- readRDS(ARMS_RDS)
    if (is.list(cached) && !is.data.frame(cached) && identical(cached$key, key)) {
      message("[mobval] reusing cached ablation arms: ", basename(ARMS_RDS))
      return(cached$arms)
    }
    message("[mobval] cached arms do not match the current settings - refitting.")
  }
  # ONE module list, from 00_config.R — see load_pipeline_modules().
  load_pipeline_modules("model", dir = ST_DIR)
  # THE SAME EXPORT CEILING AS run_all.R. run_invasion_lfo() goes through furrr whenever
  # PARALLEL_JOBS > 1, and future refuses to ship more than 500 MiB of globals by default;
  # run_all.R raises that guard to 12 GiB (a threshold, not an allocation), but this script runs
  # as its own process and did not, so the arms' back-test aborted on the closures' size.
  .old_max <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 12 * 1024^3)
  on.exit(options(future.globals.maxSize = .old_max), add = TRUE)

  dat        <- prep_all_data()
  zones_all  <- names(dat$pop); pop_vec <- dat$pop; covariates <- dat$covariates
  osrm_mat   <- load_osrm()
  # ONLY THE TWO KERNELS THE ARMS USE, built by the same calls build_all_mobility_matrices() makes
  # (03_mobility_matrices.R: build_M5(pop_vec, osrm_mat, zones_all) and
  # build_M6(osrm_mat, zones_all, "exp", kappa = 120)). Building the whole set cost twice over:
  # every back-test closure captured all of it (531 MiB, over future's 500 MiB per-worker limit,
  # so the arms' cross-validation aborted), and its save step rewrote the published
  # outputs/mobility/ files with this session's settings.
  mobility_matrices <- list(
    M5  = build_M5(pop_vec, osrm_mat, zones_all),
    M6a = build_M6(osrm_mat, zones_all, variant = "exp", kappa = 120))
  for (k in MOBVAL_KERNELS)
    if (is.null(mobility_matrices[[k]]))
      stop("[mobval] kernel ", k, " is not built; cannot fit the ablation arm.", call. = FALSE)
  gt_pmfs <- compute_all_gt_pmfs()

  # The as-of fold truncation, fitted the way run_all.R fits it. Falling back here would mean
  # the arms were corrected by a different estimator from the published models.
  trunc_asof <- tryCatch(
    .trunc_as_delay_spec(epinow2_truncation(regime = "asof", ll = dat$ll), "asof"),
    error = function(e) NULL)
  if (is.null(trunc_asof))
    stop("[mobval] the as-of truncation did not fit; the arms would not be comparable to the ",
         "published run.", call. = FALSE)

  geo <- get0("BAYES_GEO_COVARIATES", ifnotfound = c("ccvi", "d_min"))
  models <- list()
  for (k in MOBVAL_KERNELS) for (cv in list(character(0), geo)) {
    lab <- sprintf("Bayes-%s-%s", k, if (length(cv)) "geo" else "med")
    models[[lab]] <- make_bayes_lfo_model(
      k, "medium", cv, mobility_matrices, gt_pmfs, covariates, osrm_mat, zones_all,
      iter = 600L, link = "cloglog", tv = "none", linelist = dat$ll,
      trunc_delay = trunc_asof)
  }
  zone_week_outbreak <- dat$zone_week %>%
    filter(week_start >= lubridate::floor_date(OUTBREAK_START, "week",
                                               week_start = get0("WEEK_ANCHOR", ifnotfound = 1L)))
  message("[mobval] cross-validating ", length(models), " ablation arms ...")
  arms <- run_invasion_lfo(zone_week_outbreak, models,
                           delay = trunc_asof, horizons = LFO_HORIZONS,
                           analysis_date = ANALYSIS_DATE, min_atrisk_events = 0L,
                           nowcast_fn = apply_nowcast_correction, linelist = dat$ll)
  if (is.null(arms)) stop("[mobval] run_invasion_lfo() returned NULL.", call. = FALSE)
  got <- sort(unique(as.character(arms$method)))
  want <- sort(names(models))
  if (!identical(got, want))
    stop("[mobval] the cross-validation returned ", paste(got, collapse = ", "),
         " but the arms requested were ", paste(want, collapse = ", "), ".", call. = FALSE)
  saveRDS(list(key = key, arms = arms), ARMS_RDS)
  arms
}

# ---------------------------------------------------------------------------
# Step 2 - append, then PROVE nothing published moved
# ---------------------------------------------------------------------------
mobval_assert_support_unchanged <- function(pub, arms) {
  key <- c("horizon", "fold_id", "health_zone")
  # (a) The arms must COVER the published shared support and agree with it on every outcome.
  #     Comparing against one arbitrary published method (pub$method[1]) would be fragile:
  #     methods legitimately differ in coverage, which is the reason a shared support exists at
  #     all, so that check could fail for a reason that has nothing to do with the arms.
  sup <- mv_shared_support(pub) %>% distinct(across(all_of(key)), is_new_invasion) %>%
    rename(.p = is_new_invasion)
  # Each arm's SCORED cells (at risk, finite prediction), arm by arm: a cell one arm covers must
  # not stand in for another arm's NA.
  .sc <- arms %>% filter(!(was_active_before %in% TRUE), is.finite(p_invasion))
  for (m in unique(arms$method)) {
    miss <- anti_join(sup, .sc %>% filter(method == m) %>% distinct(across(all_of(key))), by = key)
    if (nrow(miss))
      stop("[mobval] the ablation arm ", m, " does not cover ", nrow(miss),
           " cells of the published shared support; appending it would shrink the support.",
           call. = FALSE)
  }
  a <- .sc %>% distinct(across(all_of(key)), is_new_invasion) %>% rename(.a = is_new_invasion)
  j <- inner_join(sup, a, by = key)
  if (any(j$.a != j$.p))
    stop("[mobval] the ablation arms disagree with the published run on ", sum(j$.a != j$.p),
         " invasion outcomes; the onset imputation seed or the data have changed.", call. = FALSE)

  # (b) appending must not move any published number
  sc <- function(x) mv_shared_support(x) %>% group_by(method, horizon) %>%
    summarise(n = n(), ap = mv_ap(is_new_invasion, p_invasion), .groups = "drop")
  before <- sc(pub)
  cm <- intersect(names(pub), names(arms))
  after <- sc(bind_rows(pub[, cm, drop = FALSE], arms[, cm, drop = FALSE]))
  cmp <- inner_join(before, after, by = c("method", "horizon"), suffix = c("_b", "_a"))
  if (nrow(cmp) != nrow(before))
    stop("[mobval] a published method vanished from the evaluation after appending the arms.",
         call. = FALSE)
  # Compare the support size PER HORIZON; comparing the set of distinct sizes would pass if two
  # horizons happened to swap values.
  sz <- function(x) x %>% group_by(horizon) %>% summarise(n = dplyr::first(n), .groups = "drop")
  if (!isTRUE(all.equal(sz(before), sz(after))))
    stop("[mobval] appending the arms changed the shared support: the arms must cover every cell ",
         "of the published support.", call. = FALSE)
  if (any(!is.finite(cmp$ap_b) | !is.finite(cmp$ap_a)))
    stop("[mobval] a published average precision is not finite; cannot verify invariance.",
         call. = FALSE)
  worst <- max(abs(cmp$ap_a - cmp$ap_b))
  if (worst > 1e-12)
    stop("[mobval] appending the arms moved a published average precision by ", worst,
         "; the comparison is no longer on the published footing.", call. = FALSE)
  message(sprintf("[mobval] support unchanged (%s cells); %d published method-horizons moved by 0",
                  paste(unique(after$n), collapse = "/"), nrow(cmp)))
  invisible(TRUE)
}

#' Cross-check both estimators against the pipeline's own evaluation table.
mobval_assert_reproduces_published <- function(d) {
  f <- file.path(ST_DIR, "outputs", "diagnostics", "invasion_evaluation.csv")
  if (!file.exists(f)) { warning("[mobval] invasion_evaluation.csv absent; check skipped.");
                         return(invisible(FALSE)) }
  ev <- utils::read.csv(f) %>%
    transmute(method, horizon = as.integer(horizon), pub_ap = auc_pr,
              pub_recall15 = recall_at_15, pub_roc = auc_roc,
              pub_rank = mean_rank_of_truth)
  # Split explicitly rather than summarise(): mv_recall_at() needs the whole group as a frame,
  # and cur_data_all() does not hand it over intact under every dplyr version.
  grp <- split(d, list(d$method, d$horizon), drop = TRUE)
  mine <- bind_rows(lapply(grp, function(g) tibble::tibble(
    method = g$method[1], horizon = as.integer(g$horizon[1]),
    ap = mv_ap(g$is_new_invasion, g$p_invasion), recall15 = mv_recall_at(g, 15L),
    roc = mv_auc_roc(g$is_new_invasion, g$p_invasion), rank = mv_mean_rank(g))))
  chk <- inner_join(mine, ev, by = c("method", "horizon"))
  if (!nrow(chk))
    stop("[mobval] no method matched invasion_evaluation.csv; the check would be vacuous.",
         call. = FALSE)
  if (any(!is.finite(chk$ap) | !is.finite(chk$pub_ap)))
    stop("[mobval] non-finite average precision in the reproduction check.", call. = FALSE)
  d_ap <- max(abs(chk$ap - chk$pub_ap))
  d_rc <- max(abs(chk$recall15 - chk$pub_recall15), na.rm = TRUE)   # recall is NA-free in practice
  d_rr <- max(abs(chk$roc - chk$pub_roc), na.rm = TRUE)
  d_mr <- max(abs(chk$rank - chk$pub_rank), na.rm = TRUE)
  if (d_rr > 1e-6) stop("[mobval] AUC-ROC does not reproduce invasion_evaluation.csv (max |d| = ",
                        signif(d_rr, 3), ").", call. = FALSE)
  if (d_mr > 1e-6) stop("[mobval] mean rank does not reproduce invasion_evaluation.csv (max |d| = ",
                        signif(d_mr, 3), ").", call. = FALSE)
  if (d_ap > 1e-9)
    stop("[mobval] average precision does not reproduce invasion_evaluation.csv (max |d| = ",
         signif(d_ap, 3), ").", call. = FALSE)
  if (d_rc > 1e-6)
    stop("[mobval] recall@15 does not reproduce invasion_evaluation.csv (max |d| = ",
         signif(d_rc, 3), ").", call. = FALSE)
  # WHAT THE LADDER PUBLISHES: per-round AUC-PR skill and its chance level over ALL rounds must
  # equal evaluate_invasion()'s auc_pr_skill_date / auc_pr_skill_date_chance (the mean rank is
  # checked against mean_rank_of_truth above and again here).
  evs <- utils::read.csv(f) %>%
    transmute(method, horizon = as.integer(horizon), pub_sel_skill = auc_pr_skill_date,
              pub_sel_chance = auc_pr_skill_date_chance, pub_sel_rank = mean_rank_of_truth)
  ds <- d
  grs <- split(ds, list(ds$method, ds$horizon), drop = TRUE)
  mine_s <- bind_rows(lapply(grs, function(g) {
    a <- mv_ap_by_date(g$is_new_invasion, g$p_invasion, g$fold_id)
    tibble::tibble(method = g$method[1], horizon = as.integer(g$horizon[1]),
                   sel_skill = a[["skill"]], sel_chance = a[["chance"]], sel_rank = mv_mean_rank(g))
  }))
  chs <- inner_join(mine_s, evs, by = c("method", "horizon"))
  if (!nrow(chs) || any(!is.finite(chs$sel_skill) | !is.finite(chs$pub_sel_skill)))
    stop("[mobval] the selection-axis check found no finite rows to compare.", call. = FALSE)
  .rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12), na.rm = TRUE)
  if (.rel(chs$sel_skill, chs$pub_sel_skill) > 1e-9 || .rel(chs$sel_chance, chs$pub_sel_chance) > 1e-9 ||
      max(abs(chs$sel_rank - chs$pub_sel_rank), na.rm = TRUE) > 1e-6)
    stop("[mobval] the per-round skill / chance / mean rank over all rounds do not reproduce ",
         "invasion_evaluation.csv's auc_pr_skill_date / auc_pr_skill_date_chance / ",
         "mean_rank_of_truth; the ladder would publish different numbers from the table.", call. = FALSE)
  message(sprintf(paste0("[mobval] reproduced %d published method-horizons exactly (AP, recall@15, ",
                         "AUC-ROC, mean rank; and the per-round skill and its chance level on %d)"),
                  nrow(chk), nrow(chs)))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Step 3 - the ladder and the paired contrasts
# ---------------------------------------------------------------------------

#' What each method's mobility term is BUILT FROM. Stated explicitly rather than inferred from a
#' name, because this classification is the whole point of the figure, and a name can mislead:
#' "Baseline-gravity" sounds phone-free but 47_onset_mi_arm.R builds it from the M4 kernel, whose
#' deterrence exponents build_M4() fits to the Flowminder relocation table - so it is
#' phone-CALIBRATED. "Gravity-B4" is the phone-free gravity yardstick: its exponents are FIXED
#' (beta = 0.5, gamma_d = 1.0 in 05_baseline_models.R) and the only parameter it estimates is a
#' scale constant k, by maximum likelihood on the training rounds' invasions (p = 1 - exp(-k
#' lambda)), from case data, never from phone data.
#' Each entry below was checked against the code that constructs that method.
#' "Adjacency-B7" is "Geography only", like Distance-B1: no contiguity matrix exists in this
#' pipeline, so run_all.R scores it through forecast_B7_adjacency(method = "distance"), which
#' ranks each zone by 1 / (1 + OSRM travel time to the nearest affected zone) — built from the
#' same road travel-time matrix as B1. It used to be classed "No movement model".
#'
#' `mv_kind` is a SECOND, ORTHOGONAL axis: whether the method reads the case data, and whether it
#' is a fitted model at all. Colour alone could not say this, and without it the case-blind
#' epicentre rules sit in the same panel as the fitted renewal models with nothing to separate
#' them - a reader could take the cohort-inflow row for a model. Three kinds:
#'   "Fitted renewal model"  the Bayesian renewal models: importation pressure from each week's
#'                           nowcasted incidence, with beta0 (and any covariate coefficients)
#'                           estimated on the training folds.
#'   "Structural rule"       uses the case data to locate the front or weight neighbours, and
#'                           fits no case-driven parameter beyond B4's maximum-likelihood scale
#'                           constant (B1 and B7 fit nothing; B4's exponents are fixed).
#'   "Case-blind ranking"    a fixed ordering computed from connectivity to the epicentre alone;
#'                           identical at every forecast origin, never sees a case.
mobval_classes <- function() .mobval_with_forms(tibble::tribble(
 ~method,                      ~mv_label,                                             ~mv_class,                ~phone, ~mv_kind,
 "Adjacency-B7",               "Nearest affected zone",                               "Geography only",         FALSE, "Structural rule",
 "Distance-B1",                "Cases weighted by inverse travel time",               "Geography only",         FALSE, "Structural rule",
 "Baseline-travel-time",       "Travel time from epicentre (case-blind)",             "Geography only",         FALSE, "Case-blind ranking",
 "Bayes-M6a-med",              "Travel-time decay",                                   "Geography only",         FALSE, "Fitted renewal model",
 "Bayes-M6a-geo",              "Travel-time decay + covariates",                      "Geography only",         FALSE, "Fitted renewal model",
 "Gravity-B4",                 "Source mass and travel time",                         "Geography + population", FALSE, "Structural rule",
 "Baseline-gravity",           "Modelled gravity from epicentre (case-blind)",        "Phone-calibrated decay", TRUE, "Case-blind ranking",
 "Bayes-M5-med",               "Radiation",                                           "Geography + population", FALSE, "Fitted renewal model",
 "Bayes-M5-geo",               "Radiation + covariates",                              "Geography + population", FALSE, "Fitted renewal model",
 "Bayes-M4-med",               "Gravity",                                             "Phone-calibrated decay", TRUE, "Fitted renewal model",
 "Bayes-M4-geo",               "Gravity + covariates",                                "Phone-calibrated decay", TRUE, "Fitted renewal model",
 "Bayes-M4c-med",              "Cohort-calibrated gravity",                           "Phone-calibrated decay", TRUE, "Fitted renewal model",
 "Bayes-M4c-geo",              "Cohort-calibrated gravity + covariates",              "Phone-calibrated decay", TRUE, "Fitted renewal model",
 "Baseline-flowminder-inflow", "Cohort inflow, 3 epicentre rows only (case-blind, static)", "Measured phone movement", TRUE, "Case-blind ranking",
 "Bayes-M8-fill-med",          "Short-trip + gravity",                                "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M8-fill-geo",          "Short-trip + gravity + covariates",                   "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M10-fill-med",         "Short-trip + radiation",                              "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M10-fill-geo",         "Short-trip + radiation + covariates",                 "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M13-fill-med",         "Cohort + gravity",                                    "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M13-fill-geo",         "Cohort + gravity + covariates",                       "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M13c-fill-med",        "Cohort + cohort-gravity",                             "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M13c-fill-geo",        "Cohort + cohort-gravity + covariates",                "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M14-fill-med",         "Cohort + radiation",                                  "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M14-fill-geo",         "Cohort + radiation + covariates",                     "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M16-fill-med",         "Cohort + relocation OD",                              "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M16-fill-geo",         "Cohort + relocation OD + covariates",                 "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M17-fill-med",         "Three-kernel consensus",                              "Measured phone movement", TRUE, "Fitted renewal model",
 "Bayes-M17-fill-geo",         "Three-kernel consensus + covariates",                 "Measured phone movement", TRUE, "Fitted renewal model"))

#' The time-varying FORMS of every registered fitted model (TV_FORMS, 00_config.R): the same
#' kernel and covariates, so the same class and phone use, with the form added to the label. A
#' featured form is then a registered ladder row like any other.
.mobval_with_forms <- function(reg) {
  forms <- get0("TV_FORMS", ifnotfound = list())
  if (!length(forms)) return(reg)
  desc <- c(rw1 = "random-walk beta", frailty = "gamma frailty",
            rw1_frailty = "random-walk beta + gamma frailty")
  fitted <- reg[reg$mv_kind == "Fitted renewal model", , drop = FALSE]
  dplyr::bind_rows(reg, dplyr::bind_rows(lapply(names(forms), function(sfx)
    dplyr::mutate(fitted, method = paste0(method, "-", sfx),
                  mv_label = paste0(mv_label, ", ", desc[[forms[[sfx]]$form]])))))
}

#' Contrasts. Each isolates ONE difference. The nested ablations come first because they are the
#' cleanest: a composite against the very kernel it is built on top of.
#'
#' THE FEATURED MODEL IS AN ARGUMENT, NOT A LITERAL. Every contrast that is about the featured
#' model, or that is quoted beside the ladder, is built from `featured` (read from the run's
#' model_selection.json) and at ITS covariate set - the same rule mv_show() uses for the ladder
#' rows - so the contrasts and the ladder describe one like-for-like comparison. They used to be
#' typed for Bayes-M14-fill-geo and "-geo" throughout, which the selection then moved away from.
#' The nested ablations list BOTH covariate sets on purpose: they are structural, not featured.
#' A FEATURED TIME-VARYING FORM differs from the kernel rows in its form as well as its kernel, so
#' its stage-1 specification (the same kernel and covariates, fixed beta) gets the same
#' contrasts too, which isolate the kernel.
mobval_contrasts <- function(featured) {
  fp <- mv_featured_parts(featured)
  at <- function(kernel) paste0("Bayes-", kernel, "-", fp$covset)
  # Every covariate arm of `kernel` that mobval_classes() registers - read, not typed.
  other_arms <- function(kernel)
    grep(paste0("^Bayes-", kernel, "-[A-Za-z0-9]+$"), mobval_classes()$method, value = TRUE)
  list(
 list(h = "Nested ablation: a composite against its own background kernel",
      pairs = list(c("Bayes-M13-fill-geo",  "Bayes-M4-geo"),
                   c("Bayes-M13-fill-med",  "Bayes-M4-med"),
                   c("Bayes-M13c-fill-geo", "Bayes-M4c-geo"),
                   c("Bayes-M13c-fill-med", "Bayes-M4c-med"),
                   c("Bayes-M14-fill-geo",  "Bayes-M5-geo"),
                   c("Bayes-M14-fill-med",  "Bayes-M5-med"),
                   c("Bayes-M8-fill-geo",   "Bayes-M4-geo"),
                   c("Bayes-M8-fill-med",   "Bayes-M4-med"),
                   c("Bayes-M10-fill-geo",  "Bayes-M5-geo"),
                   c("Bayes-M10-fill-med",  "Bayes-M5-med"))),
 # At EVERY covariate arm the registry carries, the featured model's own first: the simplified
 # Figure 2 panel C and the ladders take the best phone-free model "best overall" (any covariate
 # set), so the ratio against it has to be published whichever arm it falls in.
 list(h = "Featured model against kernels using no phone data",
      pairs = c(lapply(c(at("M5"), at("M6a"), setdiff(other_arms("M5"), at("M5")),
                         setdiff(other_arms("M6a"), at("M6a"))),
                       function(b) c(featured, b)),
                list(c(featured, "Gravity-B4"),
                     c(featured, "Adjacency-B7")))),
 if (!identical(fp$base, featured))
   list(h = "Featured model's stage-1 specification (fixed beta) against kernels using no phone data",
        pairs = c(lapply(c(at("M5"), at("M6a"), setdiff(other_arms("M5"), at("M5")),
                           setdiff(other_arms("M6a"), at("M6a"))),
                         function(b) c(fp$base, b)),
                  list(c(fp$base, "Gravity-B4"), c(fp$base, "Adjacency-B7")))),
 list(h = "Case-blind rules: measured movement against modelled movement",
      pairs = list(c("Baseline-flowminder-inflow", "Baseline-gravity"),
                   c("Baseline-flowminder-inflow", "Baseline-travel-time"))),
 list(h = "Phone data used only to calibrate a decay, against parameter-free radiation",
      pairs = list(c(at("M4"), at("M5")), c(at("M4c"), at("M5")))),
 list(h = "For scale: what the contextual covariates add",
      pairs = list(c(paste0("Bayes-", fp$kernel, "-geo"), paste0("Bayes-", fp$kernel, "-med")),
                   c("Bayes-M5-geo",       "Bayes-M5-med"))))
}

mobval_compute <- function(d, featured) {
  set_fixed_seed(MOBVAL_SEED)
  cls <- mobval_classes()
  pairs_all <- unique(do.call(rbind, lapply(mobval_contrasts(featured), function(g)
    do.call(rbind, lapply(g$pairs, function(pr)
      data.frame(group = g$h, a = pr[1], b = pr[2], stringsAsFactors = FALSE))))))
  lad <- list(); con <- list()

  for (hh in sort(unique(d$horizon))) {
    dh <- d %>% filter(horizon == hh)
    cells <- dh %>% distinct(cell, fold_id, cutoff, health_zone, is_new_invasion) %>% arrange(cell)
    if (anyDuplicated(cells$cell))
      stop("[mobval] a cell carries more than one outcome at horizon ", hh, call. = FALSE)
    y <- as.integer(cells$is_new_invasion)
    # ALL ROUNDS. The ladder REPORTS performance, and performance is reported over every round of
    # the shared support (the selection may read the post-burn-in rounds only; that choice is
    # the selection's, see INVASION_SELECTION_BASIS, 16_invasion_eval.R). Every published metric
    # here — AUC-PR skill PER ROUND, the mean rank of the invaded zones, recall@15 — is the
    # evaluation table's all-round column (auc_pr_skill_date, mean_rank_of_truth). `sel`/`bi` are
    # kept as the row index so the estimator code is one path.
    bi  <- rep(FALSE, nrow(cells))
    sel <- which(!bi)
    ys  <- y[sel]; fs <- cells$fold_id[sel]
    if (!length(sel) || sum(ys) == 0L)
      stop("[mobval] no invasion at horizon ", hh, "; the ladder is undefined there.", call. = FALSE)
    br <- mean(ys)                       # the all-round base rate (base_rate)

    # Align every method to the SAME row order as `cells`. This is what makes the bootstrap
    # paired: resampled row indices then select the same (origin x zone) cells for all methods.
    methods <- intersect(cls$method, unique(dh$method))
    if (!length(methods))
      stop("[mobval] no classified method present at horizon ", hh,
           "; mobval_classes() and the LFO have diverged.", call. = FALSE)
    P <- list(); PMs <- list()
    for (m in methods) {
      pm <- dh %>% filter(method == m) %>% arrange(cell)
      if (!identical(pm$cell, cells$cell))
        stop("[mobval] method ", m, " does not align to the shared cells at horizon ", hh,
             call. = FALSE)
      P[[m]] <- pm$p_invasion
      # Kept, not re-filtered below: the per-origin metrics (recall@15, mean rank) must be
      # computed on the very rows the alignment above verified, not on a second filter that
      # could diverge from it.
      PMs[[m]] <- pm[sel, , drop = FALSE]
    }
    have <- pairs_all[pairs_all$a %in% methods & pairs_all$b %in% methods, , drop = FALSE]

    # Zone-cluster resamples, stored as ZONE picks rather than expanded row indices: Z zones x B
    # integers instead of N rows x B (N is many times Z). Expanding inside the replicate loop
    # also lets every method and every contrast reuse ONE expansion and ONE set of skills per
    # replicate, instead of recomputing a method's skill once for the ladder and again for each
    # contrast it appears in.
    zl <- split(seq_len(nrow(cells)), cells$health_zone)
    nz <- length(zl)
    picks <- matrix(sample.int(nz, nz * MOBVAL_NBOOT, replace = TRUE), nrow = nz)
    # The same resamples, written as MULTIPLICITIES (zones x resamples), which is the form the
    # rank bootstrap needs. Both views are built from `picks`, so the skill interval and the
    # rank interval below are the same resamples of the same zones - not two bootstraps.
    M <- apply(picks, 2L, tabulate, nbins = nz)
    storage.mode(M) <- "double"

    boot_skill <- matrix(NA_real_, nrow = length(methods), ncol = MOBVAL_NBOOT,
                         dimnames = list(methods, NULL))
    boot_dif <- matrix(NA_real_, nrow = nrow(have), ncol = MOBVAL_NBOOT)
    boot_rat <- matrix(NA_real_, nrow = nrow(have), ncol = MOBVAL_NBOOT)

    for (b in seq_len(MOBVAL_NBOOT)) {
      ii <- unlist(zl[picks[, b]], use.names = FALSE)
      ii <- ii[!bi[ii]]                  # (every cell: bi is all FALSE, see above)
      yy <- y[ii]
      # A resample containing no invasion defines no precision-recall curve; it contributes
      # nothing rather than a zero, which would drag the intervals down.
      if (sum(yy) == 0L) next
      ff <- cells$fold_id[ii]
      # Each resample's per-round skill over ITS OWN event-weighted base rate: the plug-in
      # estimate of the skill on the resampled population. Within a resample every method sees
      # the same rows and the same denominator, so differences and ratios are paired.
      sk <- vapply(methods, function(m) mv_ap_by_date(yy, P[[m]][ii], ff)[["skill"]], numeric(1))
      boot_skill[, b] <- sk
      if (nrow(have)) {
        aa <- sk[have$a]; bb <- sk[have$b]
        boot_dif[, b] <- aa - bb
        boot_rat[, b] <- ifelse(is.finite(bb) & bb > 0, aa / bb, NA_real_)
      }
    }
    used <- sum(!is.na(boot_skill[1, ]))
    if (used < 0.5 * MOBVAL_NBOOT)
      stop("[mobval] only ", used, " of ", MOBVAL_NBOOT, " resamples contained an ",
           "invasion at horizon ", hh, "; the intervals would not be trustworthy.", call. = FALSE)

    # --- the second metric: the mean rank of an invaded zone, on the SAME resamples ---------
    folds <- sort(unique(fs))
    rk <- mv_rank_setup(fold_idx = match(fs, folds),
                        zone_idx = match(cells$health_zone[sel], names(zl)), y = ys, nz = nz)
    boot_rank <- matrix(NA_real_, nrow = length(methods), ncol = MOBVAL_NBOOT,
                        dimnames = list(methods, NULL))
    M1 <- matrix(1, nz, 1)                # unit multiplicities = the observed sample
    for (m in methods) {
      boot_rank[m, ] <- mv_boot_mean_rank(P[[m]][sel], rk, M)
      # PROVE the fast form is the published estimator, method by method, rather than trusting
      # the derivation: at unit multiplicities it must return mv_mean_rank() exactly.
      fast <- mv_boot_mean_rank(P[[m]][sel], rk, M1)[1]
      slow <- mv_mean_rank(PMs[[m]])
      if (!isTRUE(all.equal(fast, slow, tolerance = 1e-10)))
        stop("[mobval] the rank bootstrap does not reproduce mv_mean_rank() for ", m,
             " at horizon ", hh, " (", signif(fast, 10), " vs ", signif(slow, 10), ").",
             call. = FALSE)
    }
    # Both bootstraps must have discarded the SAME resamples - those with no invasion left.
    if (!identical(is.na(boot_rank[1, ]), is.na(boot_skill[1, ])))
      stop("[mobval] the skill and rank bootstraps disagree about which resamples carry an ",
           "invasion at horizon ", hh, "; they are no longer paired.", call. = FALSE)

    # What a mean rank would be if the watch-list were ordered at random: within an origin of
    # n candidate zones the expected rank of any zone is (n + 1) / 2, averaged over the origins
    # the estimator averages over. It is an order of magnitude above the featured model's mean
    # rank and several times the worst drawn row's, so it is
    # reported here and in the CSV rather than drawn, where it would flatten the whole axis.
    fold_n <- table(fs)
    pos_f  <- unique(as.character(fs[ys == 1L]))
    null_rank <- mean((as.numeric(fold_n[pos_f]) + 1) / 2)

    q <- function(v, pr) if (all(is.na(v))) NA_real_ else
      stats::quantile(v, pr, na.rm = TRUE, names = FALSE)
    pt <- lapply(stats::setNames(methods, methods), function(m) mv_ap_by_date(ys, P[[m]][sel], fs))
    for (m in methods) {
      pm <- PMs[[m]]
      lad[[length(lad) + 1L]] <- tibble::tibble(
        method = m, horizon = hh, ap = pt[[m]][["ap"]], base_rate = br,
        skill = pt[[m]][["skill"]], chance = pt[[m]][["chance"]],
        lo = q(boot_skill[m, ], 0.05), hi = q(boot_skill[m, ], 0.95),
        auc_roc = mv_auc_roc(ys, P[[m]][sel]), mean_rank = mv_mean_rank(pm),
        rank_lo = q(boot_rank[m, ], 0.05), rank_hi = q(boot_rank[m, ], 0.95),
        null_rank = null_rank,
        recall15 = 100 * mv_recall_at(pm, 15L),      # percentage, unlike the 0-1 check above
        n_cells = length(sel), n_events = sum(ys), n_origins = length(fold_n), n_boot = used)
    }
    for (k in seq_len(nrow(have))) {
      sa <- pt[[have$a[k]]][["skill"]]; sb <- pt[[have$b[k]]][["skill"]]
      dif <- boot_dif[k, ]; rat <- boot_rat[k, ]
      con[[length(con) + 1L]] <- tibble::tibble(
        group = have$group[k], horizon = hh, better = have$a[k], against = have$b[k],
        skill_better = sa, skill_against = sb, d_skill = sa - sb,
        d_skill_lo = q(dif, 0.05), d_skill_hi = q(dif, 0.95),
        ratio = sa / sb, ratio_lo = q(rat, 0.05), ratio_hi = q(rat, 0.95),
        p_gt0 = mean(dif > 0, na.rm = TRUE), n_boot = sum(!is.na(dif)),
        base_rate = br, n_events = sum(ys))
    }
  }
  list(ladder = bind_rows(lad) %>% left_join(cls, by = "method") %>% arrange(horizon, desc(skill)),
       contrasts = bind_rows(con) %>% distinct(better, against, horizon, .keep_all = TRUE))
}

# ---------------------------------------------------------------------------
# Step 4 - the figures
#
# THE PANEL BUILDER IS NOT HERE. MV_CLS/MV_PAL/mv_show()/MV_METRIC, mobval_ladder_frame() and
# mobval_panel() live in mobility_value_panel.R (sourced at the head of this file). They moved
# there when Figure 2 gained panel C: the main-text panel and these three appendix figures are
# the SAME ladder, and two copies of the row set, the ordering rule and the colour meanings
# are two things that can disagree without anything noticing. Everything about how a panel
# LOOKS - the house style, the two encodings, the ordering - is documented there.
#
# THREE FIGURES, ONE PANEL BUILDER. FigureS24 is the AUC-PR skill ladder, FigureS25 the same
# ladder on the mean rank of an invaded zone, and FigureS26 puts both metrics at both horizons
# in one four-panel grid. All three call mobval_panel() at its default base font with the
# facet strip on, so none of them is changed by the main-text call, which differs only in
# asking for one horizon, at a larger base size, with the horizon on the axis title.
# ---------------------------------------------------------------------------

#' Draw and save all three figures from ONE ladder table.
mobval_save_figures <- function(lad, featured) {
  A_skill <- mobval_ladder_frame(lad, "skill", featured)
  A_rank  <- mobval_ladder_frame(lad, "rank",  featured)
  hz <- sort(unique(A_skill$horizon))
  if (length(hz) != 2L)
    stop("[mobval] expected two horizons in the ladder, found ", length(hz),
         "; the four-panel figure assumes one column per horizon.", call. = FALSE)

  save_dual(mobval_panel(A_skill, "skill"), "FigureS24_mobility_value", 10.2, 5.9, dir = OUT_FIG)
  save_dual(mobval_panel(A_rank,  "rank"),  "FigureS25_mobility_value_rank", 10.2, 5.9,
            dir = OUT_FIG)

  # A, B = AUC-PR skill at each horizon; C, D = mean rank at each horizon. Row labels on the
  # left column only, so the columns are widened 1.24 : 1 - the ratio that leaves the four
  # PANEL areas equal once the left column has paid for its row labels. One legend for the
  # four, collected by patchwork.
  quad <- mobval_panel(A_skill, "skill", hz[1], TRUE) +
          mobval_panel(A_skill, "skill", hz[2], FALSE) +
          mobval_panel(A_rank,  "rank",  hz[1], TRUE) +
          mobval_panel(A_rank,  "rank",  hz[2], FALSE) +
          patchwork::plot_layout(ncol = 2, widths = c(1.24, 1), guides = "collect") +
          # The collected legend is laid out by the ASSEMBLY's theme, not by the panels', so
          # the stacking has to be restated here: without it patchwork puts the two legends
          # side by side again and the last shape key falls off the right edge.
          patchwork::plot_annotation(
            tag_levels = "A",
            theme = ggplot2::theme(legend.box = "vertical", legend.box.just = "left",
                                   legend.justification = "left")) &
          ggplot2::theme(legend.position = "top")
  save_dual(quad, "FigureS26_mobility_value_skill_and_rank", 11.6, 11.0, dir = OUT_FIG)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
run_mobility_value <- function(force_refit = FALSE) {
  dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
  pub_f <- file.path(OUT_FCST, "lfo_results.rds")
  if (!file.exists(pub_f))
    stop("[mobval] ", pub_f, " not found; run the pipeline's LFO first.", call. = FALSE)
  pub  <- readRDS(pub_f)
  # The run's featured model, from its own selection record (no typed fallback). It sets the
  # ladder's row set and covariate set, the "(featured)" tag and the featured contrasts.
  featured <- read_featured_model(file.path(OUT_KEY, "model_selection.json"))
  if (!featured %in% pub$method)
    stop("[mobval] the featured model '", featured, "' has no rows in ", basename(pub_f),
         "; model_selection.json and the LFO are from different runs.", call. = FALSE)
  message("[mobval] featured model (model_selection.json): ", featured)
  # The published cross-validation's own stamp keys the arms' cache (mobval_arms_key()).
  pub_stamp <- attr(pub, "lfo_stamp")
  if (is.null(pub_stamp))
    warning("[mobval] ", basename(pub_f), " carries no lfo_stamp; the arms' cache is keyed on ",
            "this session's settings alone.", call. = FALSE)
  arms <- mobval_fit_arms(force = force_refit, pub_stamp = pub_stamp)

  mobval_assert_support_unchanged(pub, arms)
  cm <- intersect(names(pub), names(arms))
  d  <- mv_shared_support(bind_rows(pub[, cm, drop = FALSE], arms[, cm, drop = FALSE]))
  mobval_assert_reproduces_published(d)

  n_origins <- dplyr::n_distinct(d$fold_id)
  res <- mobval_compute(d, featured)
  utils::write.csv(res$ladder, file.path(OUT_KEY, "mobility_value_ladder.csv"), row.names = FALSE)
  utils::write.csv(res$contrasts, file.path(OUT_KEY, "mobility_value_contrasts.csv"),
                   row.names = FALSE)

  mobval_save_figures(res$ladder, featured)

  # The numbers a caption needs, now that the figures carry none. They are columns of the
  # ladder CSV too; printed here so a run log is enough to write the caption from.
  nh <- res$ladder %>% distinct(horizon, n_cells, n_events, base_rate, null_rank, n_boot)
  message(sprintf("[mobval] %d weekly forecast origins; %s; %s resamples at seed %d",
                  n_origins,
                  paste(sprintf("%s wk: %d zone-weeks, %d events, base rate %.2f%%, rank under a random ordering %.1f",
                                nh$horizon, nh$n_cells, nh$n_events, 100 * nh$base_rate,
                                nh$null_rank), collapse = "; "),
                  format(MOBVAL_NBOOT, big.mark = ","), MOBVAL_SEED))

  .pf <- grep("^Bayes-(M5|M6a)-[A-Za-z0-9]+$", mobval_classes()$method, value = TRUE)
  key <- res$contrasts %>%
    filter(better == featured, against %in% .pf) %>%
    arrange(horizon, against)
  # The featured contrasts are the ones the text quotes; a pair whose method was not scored
  # is dropped silently by mobval_compute(), so their presence is asserted here.
  if (nrow(key) != length(.pf) * dplyr::n_distinct(res$ladder$horizon))
    stop("[mobval] expected the featured model's contrast against each of ",
         paste(.pf, collapse = ", "), " at every horizon; found ", nrow(key), ".", call. = FALSE)
  for (i in seq_len(nrow(key))) with(key[i, ], message(sprintf(
    "[mobval] h%s  featured vs %-14s  skill %5.1f vs %5.1f   x%.2f (%.2f-%.2f)",
    horizon, against, skill_better, skill_against, ratio, ratio_lo, ratio_hi)))
  invisible(res)
}

if (is_script_run("50_mobility_value.R")) run_mobility_value()
