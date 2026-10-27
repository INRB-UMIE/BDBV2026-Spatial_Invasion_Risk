# =============================================================================
# make_manuscript_figure2_cascade.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
# A cascade analogue of manuscript_figures/Figure2, for the 3-MONTH invasion
# cascade's leave-future-out backtest (NOT the short-term 2-week forecast).
# Same house aesthetic as make_manuscript_figures.R Figure 2; TWO panels:
#
#   A  Prioritisation — share of true invasions caught when the top-K highest-risk
#      zones are watched each round, for the cascade vs THE SAME COMPARATOR SET manuscript
#      Figure 2A publishes (FIGURE_BASELINE_METHODS), in the same colours and under the same
#      labels, vs a random watch-list. 90% origin-cluster bootstrap band on the cascade only.
#      The in-panel metric block is fs_metric_block()'s, also shared with Figure 2A.
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

# ---- house style ------------------------------------------------------------
# House style (palette, province scale, themes, save_dual) from figure_style.R — the
# single definition for the whole suite.
source(file.path(ST_DIR, "figure_style.R"))
F2_BASE <- 11.5


MANU_DIR  <- file.path(OUT_DIR, "key_outputs", "manuscript_figures")
PANEL_DIR <- file.path(MANU_DIR, "panels")

# (.wilson() removed: the capture-curve interval is a cluster bootstrap over origins, not a
#  binomial interval on pooled rows — see det_curve() below.)
# THE SCORERS COME FROM 16_invasion_eval.R. This file used to carry private copies of
# .auc_pr and .auc_roc; they were numerically equivalent to the canonical ones, but a
# figure that scores its own panel with its own implementation can drift from the table it
# sits beside, and there is no mechanism that would notice.
source(file.path(ST_DIR, "16_invasion_eval.R"))

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

# ---- THE COMPARATOR SET, IDENTICAL TO MANUSCRIPT FIGURE 2 --------------------
# WHO THE COMPARATORS ARE is FIGURE_BASELINE_METHODS (00_config.R), read through
# fs_baseline_methods(); what they are CALLED and what COLOUR they are drawn in comes from
# BASELINE_LBL / BASELINE_COL in figure_style.R. Figure 2 reads the same three things, so the
# two panels cannot name different comparators or give one hue two meanings.
#
# WHAT CHANGED AND WHY (2026-09-28). This figure carried its own private set of THREE
# comparators — fitted gravity (M4), Flowminder cohort inflow, road travel time — while
# Figure 2 published FOUR: Flowminder inflow, travel time, and the two CASE-USING nulls
# Distance-B1 and Adjacency-B7. So the cascade panel omitted both of the comparators that
# read incidence, which are the hard ones to beat, and added a gravity arm the main figure
# does not show — drawn in #009E73, the green Figure 2 gives to "Nearest affected zone". Two
# figures meant to be read side by side disagreed about both the membership and the palette.
#
# THE TWO FAMILIES, AND HOW EACH IS SCORED ON THE BACKTEST ROWS.
#
#  (1) CASE-BLIND EPICENTRE CONNECTIVITY — Baseline-gravity (M4, the fitted gravity kernel),
#      Baseline-flowminder-inflow (M_cohort, Flowminder cohort subscriber-day presence) and
#      Baseline-travel-time (OSRM road minutes). All three rank zones by how strongly they
#      connect to the EPICENTRE, differing ONLY in the connectivity matrix; none reads
#      incidence, so one static score vector serves every origin.
#      THE SCORERS ARE THE PIPELINE'S: naive_epicentre_inflow_scores() and
#      epicentre_travel_time_scores() (20_forecast_detail.R) are called directly rather than
#      reimplemented here. The previous hand copy omitted the alias harmonisation those apply,
#      so a non-canonical epicentre spelling silently changed this figure's baseline but not
#      Figure 2's.
#
#  (2) CASE-USING NULLS — Distance-B1 (neighbours' confirmed cases weighted by 1/(travel
#      time)) and Adjacency-B7 (1/(1 + travel time to the nearest already-affected zone)).
#      These read the case data, so they must be rebuilt AT EACH ORIGIN from the counts that
#      were knowable there. forecast_B1() and forecast_B7_adjacency() (05_baseline_models.R)
#      are called with exactly the arguments run_all.R's BASELINE_MODELS uses, on the as-of
#      training slice each cascade origin was itself built on (33b_cascade_calibration.R:
#      asof_counts_source(ll, zones, cutoff + 6, week_spine = weeks <= cutoff)$raw — the archived
#      snapshot of that date when one exists, else the line list re-censored to it).
#
#      THE NOWCAST CORRECTION IS DELIBERATELY NOT APPLIED to that slice, and this panel is
#      exact without it. apply_nowcast_correction() multiplies a week by a lag-derived factor
#      that is the SAME for every zone in that week. B7 reads only the support (confirmed > 0),
#      which no positive multiplier changes. B1 is linear in the final week's counts —
#      mu_cum = (sum_k W_norm^k) Y_t, since it feeds its own output forward — so a common
#      scalar on Y_t scales every zone's intensity equally and leaves the WITHIN-ORIGIN
#      RANKING untouched, which is the only thing a recall-at-K curve depends on. Applying it
#      properly would mean fitting the as-of truncation regime (epinow2_truncation(regime =
#      "asof")), which is not on disk and costs an EpiNow2 fit — in a script whose contract is
#      that it reads saved artifacts and re-simulates nothing.
#
# 03 supplies load_osrm() + harmonise_names(); 20 supplies the two epicentre scorers; 05
# supplies the two case-using baselines; 22 supplies reaggregate_asof(). All four are
# function-definition modules (their only top-level statement re-sources 00_config.R or emits
# a load message), so sourcing them runs nothing and writes nothing.
if (!exists("load_osrm", mode = "function"))
  source(file.path(ST_DIR, "03_mobility_matrices.R"))
if (!exists("naive_epicentre_inflow_scores", mode = "function") ||
    !exists("epicentre_travel_time_scores", mode = "function"))
  source(file.path(ST_DIR, "20_forecast_detail.R"))
if (!exists("forecast_B1", mode = "function") ||
    !exists("forecast_B7_adjacency", mode = "function"))
  source(file.path(ST_DIR, "05_baseline_models.R"))
if (!exists("reaggregate_asof", mode = "function"))
  source(file.path(ST_DIR, "22_daily_reissue.R"))
# reaggregate_asof() re-aggregates through aggregate_to_zone_week(), so 01 has to be loaded
# with it. 01 is a function-definition module too (its only top-level statement re-sources
# 00_config.R); sourcing it reads no data and writes nothing.
if (!exists("aggregate_to_zone_week", mode = "function"))
  source(file.path(ST_DIR, "01_data_prep.R"))
# 02 too, for the same reason: asof_counts_source() decides whether a round reads its archived
# snapshot by asking vintage_truncation(), which calls epinow2_truncation() (02). Without 02 that
# call errors, is caught, and EVERY round silently re-censors — a different slice from the one
# 33b built inside run_cascade.R, where 02 is loaded. 02 is a function-definition module; the
# per-snapshot truncation fits are disk-cached by the LFO and cascade runs that built these rounds.
if (!exists("epinow2_truncation", mode = "function"))
  source(file.path(ST_DIR, "02_epi_params.R"))

BASELINE_METHODS <- fs_baseline_methods()
CASC_KEY <- "Cascade"                      # this panel's featured series
CASC_LBL <- "3-month cascade"

# THE SIMPLIFIED FIGURE'S LABEL FOR THE CASCADE, DERIVED RATHER THAN ASSERTED. The simplified
# variant names each series by whether it uses human mobility data, so the cascade's own label
# has to state that — and that is a claim about whichever mobility kernel the CV selection
# handed the cascade, not a constant (see [[selection-rule-cascade-coupling]]: the featured-model
# pick silently sets CASCADE_KERNEL). A hard-coded "(uses human mobility data)" would go on
# claiming it after a selection that moved the cascade onto a geography-only kernel.
#
# TWO READS, BOTH FROM THE RUN THIS FIGURE DRAWS.
#   * WHICH KERNEL: cascade_calibration.json's `kernel`, written by the cascade run that
#     produced cascade_backtest_detail.csv — not CASCADE_KERNEL from 30_projection_config.R,
#     which is what a FRESH config resolution would pick today and can differ from what this
#     backtest was actually run with.
#   * WHETHER IT IS MOBILITY DATA: the mobility-value ladder's `phone` column, which
#     50_mobility_value.R sets from a table checked against the code that builds each kernel.
#     It is TRUE for a kernel that touches the mobile-network data at any point, including one
#     whose decay exponents were merely calibrated on it.
# If either read fails the parenthetical is DROPPED rather than guessed.
CASC_LBL_PLAIN <- local({
  .kern <- tryCatch({
    j <- jsonlite::fromJSON(file.path(CASC_DIAG, "cascade_calibration.json"))
    k <- j$kernel
    # Only from the cascade run that wrote the backtest detail drawn here (same cascade_run_id).
    .rid <- if ("cascade_run_id" %in% names(detail)) unique(detail$cascade_run_id) else NA_character_
    if (!identical(j$cascade_run_id %||% NA_character_, .rid[1])) k <- NA_character_
    if (is.character(k) && length(k) == 1L && nzchar(k)) k else NA_character_
  }, error = function(e) NA_character_)
  ph <- NA
  if (!is.na(.kern)) {
    f <- file.path(OUT_DIR, "key_outputs", "mobility_value_ladder.csv")
    m <- paste0("Bayes-", .kern, "-geo")
    if (file.exists(f)) {
      z <- suppressWarnings(readr::read_csv(f, show_col_types = FALSE))
      if (all(c("method", "phone") %in% names(z)) && is.logical(z$phone)) {
        v <- unique(z$phone[z$method == m])
        if (length(v) == 1L) ph <- v
      }
    }
  }
  if (isTRUE(ph))       "3-month cascade (uses human mobility data)"
  else if (isFALSE(ph)) "3-month cascade (no human mobility data)"
  else {
    message("[fig2-cascade] could not establish whether kernel '",
            if (is.na(.kern)) "?" else .kern, "' uses human mobility data; the simplified ",
            "panel's cascade label carries no mobility-data claim.")
    "3-month cascade"
  }
})
message("[fig2-cascade] comparators: ", paste(BASELINE_METHODS, collapse = ", "))

# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("make_manuscript_figure2_cascade.R")) {
  message("[fig2-cascade] sourced, not run: nothing is written. Use `Rscript make_manuscript_figure2_cascade.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

# ---- PROVENANCE: every input must come from the LATEST, COMPLETED cascade run for THIS date ----
# This script used to read cascade_backtest_detail.csv and the layer cache with no provenance at
# all, and both outlived the runs that wrote them: a failed backtest left the previous run's
# detail on disk, and the layer cache was written only on hand-runs, so this panel could draw a
# three-week-old backtest on a different week grid beside this run's tables. run_cascade.R now
# stamps a run id into each product and records the run's status; anything that does not match
# is refused rather than drawn.
CASC_STAMP <- tryCatch(jsonlite::fromJSON(file.path(CASC_DIAG, "cascade_run_stamp.json")),
                       error = function(e) NULL)
if (is.null(CASC_STAMP))
  stop("[fig2-cascade] no cascade_run_stamp.json in ", CASC_DIAG, ": the cascade outputs carry no ",
       "provenance. Re-run run_cascade.R.", call. = FALSE)
if (!identical(CASC_STAMP$status, "complete"))
  stop(sprintf("[fig2-cascade] the latest cascade run (%s) did not complete (status: %s); its outputs are partial.",
               CASC_STAMP$cascade_run_id, CASC_STAMP$status), call. = FALSE)
if (!identical(CASC_STAMP$analysis_date, format(as.Date(ANALYSIS_DATE))))
  stop(sprintf("[fig2-cascade] the cascade outputs are for analysis date %s, this figure for %s.",
               CASC_STAMP$analysis_date, format(as.Date(ANALYSIS_DATE))), call. = FALSE)
if (!all(c("cascade_run_id", "level_in_sample") %in% names(detail)) ||
    !all(detail$cascade_run_id %in% CASC_STAMP$cascade_run_id))
  stop(sprintf(paste0("[fig2-cascade] cascade_backtest_detail.csv was not written by the latest ",
                      "cascade run (%s): it is stale or unstamped. Re-run run_cascade.R with the ",
                      "backtest on."), CASC_STAMP$cascade_run_id), call. = FALSE)
# The reach table (which zones were already invaded) must be the same run's: 35 stamps it.
if (!"cascade_run_id" %in% names(reach) || !all(reach$cascade_run_id %in% CASC_STAMP$cascade_run_id))
  stop(sprintf(paste0("[fig2-cascade] cascade_reach_scores_all_zones.csv was not written by the ",
                      "latest cascade run (%s): it is stale or unstamped. Re-run run_cascade.R."),
               CASC_STAMP$cascade_run_id), call. = FALSE)
# THE LEVEL IS IN-SAMPLE whenever delta was fitted on these same origins and K (run_cascade.R
# records it per row): calibration-in-the-large is then ~1 by construction and must say so.
LEVEL_IN_SAMPLE <- isTRUE(all(as.logical(detail$level_in_sample)))

# ONE TRAVEL-TIME MATRIX FOR THE WHOLE PANEL. load_osrm() is the pipeline's own loader
# (03_mobility_matrices.R); its default kind is "travel_time" in minutes — the same matrix
# run_all.R passes as osrm_mat, so this figure's travel-time baseline is built on the same
# geography as manuscript Figure 2's. Loaded ONCE: the loader refits the great-circle gap-fill
# on every call, and two calls would be two fits.
OSRM_FIG <- tryCatch(load_osrm("travel_time"), error = function(e) {
  warning("[fig2-cascade] OSRM travel-time matrix unavailable: ", conditionMessage(e),
          call. = FALSE); NULL })

# ---- (1) the CASE-BLIND epicentre-connectivity scores (one vector per method) ----
static_baseline_scores <- function() {
  MD <- file.path(OUT_DIR, "mobility")
  .k <- function(id) {
    f <- file.path(MD, sprintf("mobility_%s.rds", id))
    if (file.exists(f)) readRDS(f) else NULL
  }
  wp <- read_csv(file.path(DATA_DIR, "worldpop", "processed", "worldpop__pop_count__static.csv"),
                 show_col_types = FALSE)
  pv <- setNames(wp$pop_count, wp$nom)
  om <- OSRM_FIG

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

  list(`Baseline-gravity` = if (is.null(M4)) NULL else
         naive_epicentre_inflow_scores(M4, EPICENTRE_ZONES, pv, zall),
       `Baseline-flowminder-inflow` = if (is.null(MC)) NULL else
         naive_epicentre_inflow_scores(MC, EPICENTRE_ZONES, pv, zall),
       `Baseline-travel-time` = if (is.null(om)) NULL else
         epicentre_travel_time_scores(om, EPICENTRE_ZONES, zall))
}

# ---- (2) the CASE-USING nulls, rebuilt at every origin -----------------------
# One call per origin per method, from the layer cache the cascade itself was run on, so the
# counts these read are the counts that origin's cascade fit read. The cache carries `ll`
# precisely so a leakage-free reconstruction is possible without the line list; refusing
# rather than falling back to the FINAL revised counts is deliberate — that fallback would
# hand both case-using nulls cases reported after their own origin and quietly flatter them
# against the cascade arm, which gets no such help.
case_baseline_scores <- function(cuts, want) {
  if (!length(want)) return(NULL)
  # ONLY THE CANONICAL LAYER (beside CASC_DIAG), never CASCADE_LAYER_CACHE. run_cascade.R writes
  # that optional build cache only when it REBUILDS the layer, before it stamps the run id, and
  # never re-writes it when it loads it — so a layer read from there carries no cascade_run_id and
  # the provenance check below stopped this figure every time the variable was set. The canonical
  # file is re-written by every run_cascade.R run, after the run id is set, with the layer that
  # run actually used.
  cache <- file.path(OUT_DIR, "cascade", "layer_cache_v2.rds")
  if (!file.exists(cache))
    stop("[fig2-cascade] Panel A needs the case-using nulls (", paste(want, collapse = ", "),
         ") and they are rebuilt per origin from the cascade layer cache, which is not at ",
         cache, ". Run run_cascade.R (it writes the layer it used there).", call. = FALSE)
  layer <- readRDS(cache)
  if (!identical(layer$cascade_run_id %||% NA_character_, CASC_STAMP$cascade_run_id))
    stop(sprintf(paste0("[fig2-cascade] the layer cache at %s was not used by the latest cascade ",
                        "run (%s vs %s), so its as-of counts may sit on a different week grid from ",
                        "the backtest. Re-run run_cascade.R (it rewrites the cache)."),
                 cache, layer$cascade_run_id %||% "unstamped", CASC_STAMP$cascade_run_id),
         call. = FALSE)
  if (is.null(layer$ll))
    stop("[fig2-cascade] the layer cache at ", cache, " carries no line list (`ll`), so each ",
         "origin's training counts cannot be reconstructed as-of. Re-build the layer with ",
         "run_cascade.R.", call. = FALSE)
  zones_all <- layer$zones_all
  om        <- layer$osrm_mat
  weeks     <- sort(unique(layer$zone_week_nc$week_start))

  # ONE GEOGRAPHY ACROSS THE PANEL. The travel-time baseline is built on OSRM_FIG (a fresh
  # load_osrm(), matching manuscript Figure 2) and these two on the matrix the cascade run was
  # driven by. They are the same object on this run — checked here rather than assumed,
  # because if a cached layer ever went stale against the OSRM build, three comparators in one
  # panel would sit on two different road networks and nothing on the page would say so.
  if (!is.null(OSRM_FIG)) {
    .ok <- identical(dim(OSRM_FIG), dim(om)) &&
           setequal(rownames(OSRM_FIG), rownames(om)) &&
           setequal(colnames(OSRM_FIG), colnames(om))
    .d  <- if (.ok) max(abs(OSRM_FIG - om[rownames(OSRM_FIG), colnames(OSRM_FIG)]),
                        na.rm = TRUE) else NA_real_
    if (!.ok || !is.finite(.d) || .d > 1e-6)
      warning(sprintf(paste0("[fig2-cascade] the cascade layer's travel-time matrix differs ",
                             "from the current OSRM build (%s). Distance-B1 and Adjacency-B7 ",
                             "are on the layer's geography while Baseline-travel-time is on ",
                             "the current one; re-build the cascade layer."),
                      if (.ok) sprintf("max|diff| = %.4g min", .d) else "different shape"),
              call. = FALSE)
    else
      message(sprintf("[fig2-cascade] travel-time geography: layer matches the current OSRM build (max|diff| = %.3g).", .d))
  }

  purrr::map_dfr(cuts, function(cut) {
    # The SAME as-of slice 33b_cascade_calibration.R builds for this origin, through the same
    # dispatcher (asof_counts_source(): the archived snapshot of cutoff + 6 when one exists, else
    # the line list re-censored to that date), on the weeks up to the cutoff. Only the RAW counts
    # are used — see the header for why the nowcast multiplier is not applied.
    zw <- if (exists("asof_counts_source", mode = "function"))
            asof_counts_source(layer$ll, zones_all, cut + 6, week_spine = weeks[weeks <= cut])$raw
          else reaggregate_asof(layer$ll, zones_all, cut + 6, week_spine = weeks[weeks <= cut])
    zw <- zw[zw$week_start <= cut, , drop = FALSE]
    if (!"confirmed" %in% names(zw) || !nrow(zw))
      stop("[fig2-cascade] the as-of reconstruction at ", format(cut), " returned no ",
           "`confirmed` counts, so the case-using nulls have nothing to read.", call. = FALSE)
    zw$confirmed_nc <- zw$confirmed
    Yw <- zone_week_to_wide(zw, zones_all)
    ti <- length(unique(zw$week_start))
    # NEITHER NULL MAY BE SILENTLY BLIND. B1 reads only week `ti`; B7 reads the whole history.
    # If the reconstruction ever came back empty — a renamed line-list column, a censoring
    # date that excludes everything — B1's intensity would be 0 for every zone and B7 would
    # find no affected zone, both comparators would become flat ties, and their capture curves
    # would silently collapse to the random line with no error anywhere. Two comparators
    # quietly replaced by noise in a published panel is exactly the failure this figure's
    # other guards exist to prevent.
    if (!any(Yw[, min(ti, ncol(Yw))] > 0))
      stop("[fig2-cascade] the as-of counts at ", format(cut), " are zero in every zone for ",
           "the final training week, so Distance-B1 would rank every zone equally.",
           call. = FALSE)
    if (!any(rowSums(Yw[, seq_len(min(ti, ncol(Yw))), drop = FALSE] > 0) > 0))
      stop("[fig2-cascade] no zone has a case in the as-of history at ", format(cut),
           ", so Adjacency-B7 would rank every zone equally.", call. = FALSE)
    fc <- list()
    if ("Distance-B1" %in% want)
      fc[["Distance-B1"]] <- forecast_B1(Yw, om, ti, K_HORIZON, zones_all, alpha = 1)
    if ("Adjacency-B7" %in% want)
      fc[["Adjacency-B7"]] <- forecast_B7_adjacency(Yw, om, ti, K_HORIZON, zones_all,
                                                    method = "distance")
    purrr::imap_dfr(fc, function(f, nm)
      tibble::tibble(cutoff = cut, health_zone = f$health_zone, method = nm,
                     score = as.numeric(f$p_invasion)))
  })
}

# ---- assemble one long (cutoff x zone x method) score frame ------------------
# A MISSING COMPARATOR IS FATAL, exactly as it is in Figure 2A (prioritisation_curves()).
# Panel A's whole claim is that the cascade beats this set of nulls; dropping one on a
# warning publishes a figure that looks complete while making a weaker claim, with the stale
# copy of the full panel sitting on disk beside it.
.stat <- static_baseline_scores()
.want_case <- intersect(BASELINE_METHODS, c("Distance-B1", "Adjacency-B7"))
.case <- case_baseline_scores(fold_ord, .want_case)

bl_long <- dplyr::bind_rows(
  purrr::imap_dfr(.stat[intersect(BASELINE_METHODS, names(.stat))], function(v, nm) {
    if (is.null(v) || !any(is.finite(v) & v > 0))
      stop("[fig2-cascade] the ", nm, " comparator is unavailable or all-zero, so Panel A ",
           "cannot draw the comparator set manuscript Figure 2 publishes. Re-run the ",
           "mobility build.", call. = FALSE)
    tidyr::expand_grid(cutoff = fold_ord, health_zone = names(v)) %>%
      dplyr::mutate(method = nm, score = as.numeric(v[health_zone]))
  }),
  .case)

.miss <- setdiff(BASELINE_METHODS, unique(bl_long$method))
if (length(.miss))
  stop("[fig2-cascade] Panel A needs ", length(BASELINE_METHODS), " comparators but ",
       length(.miss), " produced no scores: ", paste(.miss, collapse = ", "), ".",
       call. = FALSE)
# AT MOST ONE SCORE PER (ORIGIN, ZONE, METHOD), CHECKED BEFORE THE JOIN. A duplicate key here
# — a kernel whose column names repeat a zone, say — would make the left_join below FAN OUT
# the backtest rows for that comparator alone, so it would be scored on a larger support than
# the cascade and its recall-at-K would not be comparable. The shared-support assertion on
# `rnd` further down would eventually catch it, but only as a puzzling message about the
# random-targeting reference.
if (anyDuplicated(bl_long[, c("cutoff", "health_zone", "method")]))
  stop("[fig2-cascade] a comparator supplied more than one score for the same (origin, zone). ",
       "The join onto the backtest rows would duplicate them and put that comparator on a ",
       "different support from the cascade.", call. = FALSE)

# THE SHARED SUPPORT IS THE CASCADE'S OWN ROWS. Every comparator is scored on exactly the
# (origin x zone) cells the backtest scored the cascade on, so the random-targeting reference
# k / n_atrisk is one curve rather than one per family — the property Figure 2A asserts and
# checks. A zone the comparator has no entry for scores 0 (a tie at the bottom of that
# origin) rather than dropping the row, because dropping it would move the support.
rows <- dplyr::bind_rows(
  detail %>% dplyr::transmute(cutoff, health_zone, y, method = CASC_KEY, score = p_invasion),
  detail %>% dplyr::select(cutoff, health_zone, y) %>%
    tidyr::expand_grid(method = BASELINE_METHODS) %>%
    dplyr::left_join(bl_long, by = c("cutoff", "health_zone", "method")) %>%
    # `unscored` is the join MISS, recorded before the fill. It is not the same thing as a
    # score of zero, which several comparators legitimately produce: Adjacency-B7 assigns
    # exactly 0 to every zone that is already affected, and a measured-movement kernel has
    # true zeros wherever it observed no travel. Counting `score == 0` instead — as this did
    # — reported those as missing data in the build log.
    dplyr::mutate(unscored = is.na(score), score = dplyr::coalesce(score, 0))
) %>% dplyr::mutate(unscored = dplyr::coalesce(unscored, FALSE),
                    method = factor(method, levels = c(CASC_KEY, BASELINE_METHODS)))
local({
  z <- rows %>% dplyr::filter(method != CASC_KEY) %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(n_miss = sum(unscored), n_zero = sum(score == 0), n = dplyr::n(),
                     .groups = "drop")
  for (i in seq_len(nrow(z)))
    message(sprintf("[fig2-cascade] %s: %d/%d cells unscored (filled with 0 and ranked last); %d cells score 0 in all.",
                    z$method[i], z$n_miss[i], z$n[i], z$n_zero[i]))
})

# =============================================================================
# Panel A — prioritisation (capture curve): cascade vs baselines vs random
# =============================================================================
# recall per origin (ranked WITHIN origin, ties "max"), then averaged over origins: the estimator
# of compute_detection_curve() (20_forecast_detail.R), which Figure 2A's curves are read from.
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
      set_fixed_seed(as.integer((as.numeric(seed) + 7919 * k) %% 2147483647))
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

# COLOURS, LABELS AND MEMBERSHIP ALL COME FROM figure_style.R, shared with manuscript
# Figure 2, so the same comparator is the same hue and the same words in both panels. This
# block used to hold a private four-entry palette that gave "Gravity model (fitted flows)"
# the green Figure 2 reserves for "Nearest affected zone".
#
# LEGEND ORDER: the cascade first, then FIGURE_BASELINE_METHODS. Only the ORDER is fixed here;
# the colours and labels are resolved per variant inside .mk_pA(), because the simplified
# figure uses the plain-language wording for the same comparators.
MODEL_ORD <- c(CASC_KEY, BASELINE_METHODS)
# The PUBLISHED labels, for the figure's data file: one method name per row, in the words the
# published panel uses. (The simplified variant's wording is BASELINE_LBL_PLAIN.)
MODEL_LBL <- fs_panelA_scales(CASC_KEY, CASC_LBL, BASELINE_METHODS,
                              who = "fig2-cascade")$labels

curves <- rows %>% dplyr::group_split(method) %>%
  purrr::map_dfr(function(d) {
    out <- det_curve(d, "score")
    if (is.null(out)) {
      # An empty curve here means the origin set carries no events at all, which would have
      # stopped det_curve() on the cascade arm too — so it is a broken input, not a thin panel.
      stop("[fig2-cascade] ", as.character(d$method[1]), " produced no capture curve: the ",
           "backtest rows carry no invasions.", call. = FALSE)
    }
    out$method <- as.character(d$method[1])
    out
  }) %>%
  dplyr::mutate(method = factor(method, levels = MODEL_ORD))
feat_c <- curves %>% filter(method == CASC_KEY)             # bootstrap band on the cascade only
# ONE RANDOM-TARGETING REFERENCE, NOT ONE PER METHOD. recall_random is k / n_atrisk, a property
# of the ORIGIN DESIGN and not of any method, so it must be identical across the curves — which
# it is only because every comparator above was scored on the cascade's own backtest rows. The
# same invariant is asserted in Figure 2A (prioritisation_curves()), where it fires when the
# families end up on different supports; here it would fire if a comparator ever started
# dropping rows instead of scoring them 0.
rnd <- curves %>% dplyr::distinct(k, recall_random) %>% dplyr::arrange(k)
if (nrow(rnd) != dplyr::n_distinct(curves$k))
  stop("[fig2-cascade] the comparators were scored on different row sets, so Panel A would ",
       "draw one method's random-targeting reference beside another's curve.", call. = FALSE)

# Pooled summary (all origins) for the in-panel annotation. THE SIX LINES, THEIR WORDING AND
# THEIR ORDER ARE fs_metric_block()'s (figure_style.R) — the same block manuscript Figure 2A
# carries, so a reader comparing the two panels reads one set of quantities in one sequence.
#
# The base rate is part of that block and was missing here, which mattered more on this panel
# than on Figure 2's: AUC-PR skill is average precision DIVIDED BY the base rate, and the
# cascade's +6-week window has a base rate several times the 2-week forecast's, so the two
# figures' skill multiples are not on a common denominator and the reader had no way to see it.
#
# THE EQUAL-COUNT ECE IS NO LONGER DRAWN. It is a genuine extra measure — whether the model is
# calibrated ACROSS THE RANGE rather than only on average — but Figure 2A does not carry it,
# and a seventh line under the same heading would be a second quantity a reader has to
# reconcile across two panels that are meant to be read as one. It is still computed and
# reported to the log below, and the reliability data behind it is published by
# 33b_cascade_calibration.R (cascade_calibration_reliability.csv).
#
# EQUAL-COUNT (rank-binned) ECE, NOT the equal-width .ece() in 16_invasion_eval.R. Renamed so
# it cannot silently shadow the canonical one now that this file sources 16: on a rare-event
# probability distribution the two give materially different numbers, and which is meant has
# to be visible at the call site. Equal-count bins are the right choice HERE because the
# cascade backtest pools few origins and equal-width bins would put almost every row in the
# first bin.
.ece_equalcount <- function(p, y, nb = 10L) {
  ok <- is.finite(p) & is.finite(y); p <- p[ok]; y <- as.integer(y[ok])
  if (length(p) < nb) return(NA_real_)
  b <- cut(rank(p, ties.method = "first"), breaks = nb, labels = FALSE)
  # WEIGHTED by bin size, which is what makes it an EXPECTED calibration error. Equal-count
  # bins make the weights equal anyway, but only when n divides exactly; weighting keeps the
  # definition right when it does not.
  w <- as.numeric(table(b)) / length(p)
  sum(w * abs(tapply(p, b, mean) - tapply(y, b, mean)), na.rm = TRUE)
}
# AUC-PR skill PER ORIGIN (16's .auc_pr_by_date(): AP within each origin, event-weighted, over the
# base rate) and its random-ordering level: the axis Figure 2A and the selection use. Pooling the
# origins' rows in one list would reward the probability LEVEL on the origins with more first
# cases, not the ordering within an origin.
.apd_c     <- .auc_pr_by_date(detail$p_invasion, detail$y, detail$cutoff)
pooled_br  <- mean(detail$y, na.rm = TRUE)
pooled_cil <- mean(detail$p_invasion, na.rm = TRUE) / max(pooled_br, 1e-12)
pooled_ece <- .ece_equalcount(detail$p_invasion, detail$y)
disc_lab <- fs_metric_block(sprintf("Leave-future-out (+%d wk, %d origins)", K_HORIZON, N_ORIGINS),
                            skill     = .apd_c[["skill"]],
                            chance    = .apd_c[["chance"]],
                            base_rate = pooled_br,
                            roc       = .auc_roc(detail$p_invasion, detail$y),
                            brier     = mean((detail$p_invasion - detail$y)^2, na.rm = TRUE),
                            calib     = pooled_cil)
# CALIBRATION-IN-THE-LARGE IS IN-SAMPLE HERE, and the panel says so. delta was fitted so the
# pooled predicted count matches the observed count on exactly these origins and K, so this line
# is ~1 by construction; under a "Leave-future-out" heading it read as a held-out calibration
# result. The ranking lines above it are not what delta was fitted to. The wording comes from
# fs_metric_block() (figure_style.R), so the substitution is checked rather than assumed.
if (LEVEL_IN_SAMPLE) {
  .lab_ins <- sub("Calibration (pred/obs)", "Calibration (pred/obs, in-sample)", disc_lab, fixed = TRUE)
  if (identical(.lab_ins, disc_lab))
    stop("[fig2-cascade] fs_metric_block() no longer prints 'Calibration (pred/obs)', so the ",
         "in-sample label could not be attached to the cascade's calibration line. Update the ",
         "substitution here before publishing.", call. = FALSE)
  disc_lab <- .lab_ins
}
message(sprintf("[fig2-cascade] 2A: base rate %.2f%%, calibration-in-the-large %.2f, ",
                100 * pooled_br, pooled_cil),
        sprintf("equal-count ECE %.3f (computed, not drawn — Figure 2A does not carry it)",
                pooled_ece))

#' Panel A. ONE code path for the published panel and the simplified variant, exactly as
#' Figure 2A's .mk_p2a() is: the curves, the band, the colours and the metric block are
#' identical and only the words and one annotation differ.
#'
#' @param plain      plain-language legend labels and axis title (the simplified figure).
#' @param watch_note draw the italic "random watch-list" note on the dashed reference line.
#'   The LINE itself is drawn either way — only the note is dropped.
.mk_pA <- function(plain = FALSE, watch_note = TRUE) {
  sc <- fs_panelA_scales(CASC_KEY, if (plain) CASC_LBL_PLAIN else CASC_LBL,
                         BASELINE_METHODS, who = "fig2-cascade", plain = plain)
  ggplot(curves, aes(k, recall, colour = method)) +
  geom_line(data = rnd, aes(k, recall_random), linetype = "22", colour = FAINT, linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = feat_c, aes(ymin = recall_lo, ymax = recall_hi, fill = method), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(data = feat_c, size = 1) +
  # `+ NULL` is a no-op in ggplot2, so the simplified panel is the published one with this
  # single layer absent — not a separately maintained plot.
  (if (watch_note)
     annotate("text", x = 15, y = 0.075, label = "random watch-list", colour = MUTED,
              size = 3.5, angle = 8)
   else NULL) +
  annotate("text", x = 0.5, y = 0.99, hjust = 0, vjust = 1, label = disc_lab, size = 3.4, colour = INK, lineheight = 0.95) +
  scale_colour_manual(values = sc$values, labels = sc$labels, name = NULL,
                      breaks = levels(curves$method)) +
  scale_fill_manual(values = sc$values, guide = "none") +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = if (plain) "Prioritisation of top K zones"
            else      "Zones actively monitored per round (K)",
       y = "Share of true invasions caught") +
  # ONE KEY PER ROW, as in Figure 2A and for the same reason: five entries, the longest 39
  # characters, laid out in two rows measure 7.46 in against the ~4.1 in this panel is given,
  # and ggplot neither clips nor wraps an over-wide guide box — it draws it straight across
  # panel B, on top of panel B's own legend. That is what the published version of this figure
  # did. .assert_legend_fits() below re-measures it at build time.
  guides(colour = guide_legend(ncol = 1)) +
  theme_pub(F2_BASE) + theme(legend.position = "top")
}
pA   <- .mk_pA(plain = FALSE, watch_note = TRUE)    # as published
pA_s <- .mk_pA(plain = TRUE,  watch_note = FALSE)   # Figure2_cascade_simple

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
  labs(x = paste0("Predicted invasion probability, P(first case)",
                  if (LEVEL_IN_SAMPLE) ", after the in-sample scale factor" else ""), y = NULL) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_pub(F2_BASE) + theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(9, "pt"), strip.clip = "off",
                      axis.text.y = element_text(size = 9.6, colour = INK), axis.text.x = element_text(size = 8.8), legend.position = "top",
                      # extra right margin: the last facet's unclipped strip title overhangs the panel
                      plot.margin = margin(6, 40, 6, 6))

# =============================================================================
# Assemble + save (no title, no caption — tags only)
# =============================================================================
# 14.6 in, not Figure 2's 13.8. The A:B WIDTH RATIO is Figure 2's and stays there, but this
# figure's panel B carries FIVE faceted rounds against Figure 2's three, and each facet's
# zone-label column is inelastic — so at 13.8 in panel A was left 2.75 in of room for a 2.91 in
# legend and .assert_legend_fits() stopped the build. The extra 0.8 in goes to panel A in the
# ratio below, which is the smallest change that keeps both the legend and the published ratio.
F2C_W      <- 14.6
F2C_WIDTHS <- c(0.85, 1.5)   # panel A : panel B, as published
# 5.85 in, not the 5.2 this figure used to be saved at. Panel A's legend went from two rows to
# five (guides(ncol = 1) above), which is 0.65 in taller; without giving the figure that back,
# the panels would have paid for the legend out of their own height.
F2C_H      <- 5.85

save_dual(pA, "F2A_cascade_prioritisation",  5.2, F2C_H)
save_dual(pB, "F2B_cascade_forecast_vs_outcome", 2.9 * length(fold_ord) + 1.0, F2C_H)

# BOTH LEGENDS, MEASURED AGAINST THE WIDTH THEY WILL ACTUALLY BE LAID OUT IN — the check
# Figure 2 runs, now that it lives in figure_style.R and this file can call it. The published
# version of this figure had panel A's two-row legend running straight across panel B.
.assert_legend_fits(list(pA, pB), F2C_WIDTHS, F2C_W, labels = c("2A", "2B"),
                    who = "fig2-cascade")

.assemble <- function(p) (p | pB) + plot_layout(widths = F2C_WIDTHS) +
  plot_annotation(tag_levels = "A")
save_dual(.assemble(pA), "Figure2_cascade", F2C_W, F2C_H, dir = MANU_DIR)

# ---- THE SIMPLIFIED VARIANT, alongside the published figure ------------------
# SAME DATA, SAME PANELS, FEWER WORDS: the same curves with plain-language legend labels and
# without the "random watch-list" note (the dashed reference LINE stays), beside the identical
# panel B. It has NO panel C: the mobility-value ladder Figure2_simple's panel C summarises is
# scored on the 2-week leave-future-out folds, and there is no cascade analogue of it on disk —
# building one means re-running this backtest under phone-free kernels, which is a modelling
# run and not a figure. Drawing Figure 2's ladder under cascade panels would put a 2-week
# result under a 13-week figure.
# 15.4 in, not the published figure's 14.6. Its panel A legend is WIDER, not narrower: the
# plain-language labels drop the vendor name but the featured series gains "(uses human
# mobility data)", which takes the longest label from 39 to 41 characters and the guide box
# from 2.91 in to 3.25 in — against the 3.04 in panel A has at 14.6 in, so the check below
# stopped the build. The A:B ratio is unchanged and the extra 0.8 in goes to panel A through
# it. Shortening the label instead would mean dropping the mobility-data claim from the one
# series the figure is about, which is the opposite of what the variant is for.
F2CS_W <- 15.4
.assert_legend_fits(list(pA_s, pB), F2C_WIDTHS, F2CS_W, labels = c("2A-s", "2B"),
                    who = "fig2-cascade-simple")
save_dual(.assemble(pA_s), "Figure2_cascade_simple", F2CS_W, F2C_H, dir = MANU_DIR)

# ---- THE CAPTION, WRITTEN FROM THE SAME RUN AS THE FIGURE --------------------
# EVERY NUMBER IN IT IS READ OFF THE OBJECTS THE PANELS WERE DRAWN FROM, for the reason given
# at the same point in make_manuscript_figures.R: a caption typed once and kept beside a figure
# that is rebuilt whenever the pipeline runs is the easiest way for a manuscript to quote last
# month's numbers under this month's figure, and nothing about it looks wrong. Edit the PROSE
# here; do not paste numbers into it.
local({
  .r  <- function(kk) curves %>% dplyr::filter(method == CASC_KEY, k == kk)
  .pc <- function(x) sprintf("%.0f%%", 100 * x)
  .dt <- function(d) sub("^0", "", format(as.Date(d), "%d %B %Y"))
  r10 <- .r(10); r20 <- .r(20)
  ori <- lfo_origin(range(fold_ord))
  # THE CEILING ON RECALL AT K, EXACTLY AS THE CURVE COMPUTES IT. Per-round recall cannot
  # exceed min(1, K / (invasions in that round)), and the published curve is the MEAN of the
  # per-round values — so the ceiling is the mean of the per-round ceilings, not a ratio of
  # the totals. Worth stating in the caption: a reader comparing this panel with Figure 2A's
  # would otherwise read a lower curve as a worse ranking, when most of the difference is that
  # a ten-zone budget is a smaller fraction of what there is to catch over six weeks.
  .np <- detail %>% dplyr::group_by(cutoff) %>%
    dplyr::summarise(np = sum(y, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(np > 0) %>% dplyr::pull(np)
  .ev_per_round <- mean(.np)
  .ceil10 <- mean(pmin(1, 10 / .np)); .ceil20 <- mean(pmin(1, 20 / .np))
  n_topk <- dplyr::n_distinct(topk$health_zone[topk$cutoff == fold_ord[1]])
  .settle <- get0("LFO_RELIABLE_EVAL_AGE_DAYS", ifnotfound = NA_integer_)
  txt <- paste0(
"**Figure 2 (longer horizon). The same question, scored ", K_HORIZON, " weeks ahead.**\n",
"\n",
"**(A)** As in Figure 2A, but the ranking comes from the three-month spread simulation rather ",
"than the two-week forecast. It is scored over ", K_HORIZON, " weeks — the longest window whose ",
"outcome had been settled (closed for at least ", .settle, " days) at every round — not over the full three ",
"months the simulation runs. Over the ", N_ORIGINS,
" rounds issued between ", sub("^0", "", format(ori[1], "%d %B")), " and ", .dt(ori[2]),
", watching the top 10 zones each round would have covered ", .pc(r10$recall), " of the zones ",
"that were genuinely reached in the following ", K_HORIZON, " weeks (90% interval ",
.pc(r10$recall_lo), " to ", .pc(r10$recall_hi), "), and the top 20 would have covered ",
.pc(r20$recall), " (", .pc(r20$recall_lo), " to ", .pc(r20$recall_hi), "). Choosing zones at ",
"random would have covered ", .pc(r10$recall_random), " and ", .pc(r20$recall_random),
", shown by the grey dashed line. The ", c("one", "two", "three", "four", "five", "six")[length(BASELINE_METHODS)] %|NA|% as.character(length(BASELINE_METHODS)),
" comparison lines and the shaded band are exactly ",
"as in Figure 2A, so the two figures can be read against one another.\n",
"\n",
"**(B)** All ", N_ORIGINS, " rounds: the ", n_topk, " zones ranked highest in each, and what ",
"actually happened to each one — a first case inside the ", K_HORIZON, " weeks the ranking ",
"covered (dark red), a first case later in the outbreak (amber), or no case by ",
.dt(DATA_END_DATE), ", the last day of data used to fit (grey). Bar length is the probability ",
"the ranking gave that zone.\n",
"\n",
"A fixed watch-list is a harder test over a longer window, which is part of why the top line ",
"sits lower than Figure 2A's. More zones are reached in ", K_HORIZON, " weeks than in two — an ",
"average of ", sprintf("%.1f", .ev_per_round), " per round here, a base rate of ",
sprintf("%.2f%%", 100 * mean(detail$y)), " — so a list of ten zones cannot cover more than ",
.pc(.ceil10), " of them however well it is ranked",
if (.ceil20 < 0.999) paste0(", and a list of twenty no more than ", .pc(.ceil20)) else "",
". Here ", paste(range(as.integer(table(detail$cutoff))), collapse = " to "),
" zones were still at risk in a round.\n",
"\n",
"*Metrics.* Recall at K, averaged over rounds, with 90% intervals from a cluster bootstrap ",
"over rounds",
# A bootstrap over few rounds has few distinct resamples (choose(2n - 1, n): 10 for 3 rounds), so
# its percentile interval is coarse; say so rather than let it read as a precise interval.
if (N_ORIGINS < 5L) sprintf(paste0(" (with only %d rounds the bootstrap has %d distinct resamples, ",
                                   "so these intervals are coarse)"),
                            N_ORIGINS, as.integer(choose(2L * N_ORIGINS - 1L, N_ORIGINS))) else "",
". The text block reports AUC-PR skill (average precision within each round, combined ",
"over rounds weighted by their positive zone-rounds and divided by the base rate, beside the ",
"score of a random ordering), AUC-ROC, the Brier score and calibration-in-the-large (predicted ",
"divided by observed number of positive zone-rounds), over the ",
format(nrow(detail), big.mark = ","), " zone-rounds the backtest scored, in which ",
dplyr::n_distinct(detail$health_zone[detail$y == 1]), " zones had their first case inside a ",
"scored window (", sum(detail$y), " positive zone-rounds, because the ", K_HORIZON,
"-week windows of consecutive rounds overlap). ",
if (LEVEL_IN_SAMPLE) paste0(
"Calibration-in-the-large is in-sample: the simulation's single scale factor was fitted so that ",
"the predicted number of positive zone-rounds matches the observed number over these same rounds, so it ",
"is close to 1 by construction (the Brier score depends on that level too). The ranking ",
"measures — recall at K, AUC-PR and AUC-ROC — are not what the factor was fitted to. ") else "",
"There is no panel C: the comparison of ",
"methods in Figure 2C is scored on the two-week rounds, and no three-month equivalent of it ",
"exists.\n")
  f <- file.path(MANU_DIR, "Figure2_cascade_simple_caption.md")
  writeLines(txt, f)
  message("  wrote ", basename(f))
})

# THE FIGURE'S OWN DATA, under the name that claims to be it. 39_cascade_eval_figure.R used to
# write "Figure2_cascade_data.csv" holding the per-origin backtest summary — a frame with none
# of this figure's columns in it — while this, the published figure, shipped no data file at all.
# Panel A is the prioritisation curves (with the bootstrap band on the cascade arm); panel B is the
# per-origin top-12 zones and their realised outcomes. Both are written, long, with a `panel` key.
local({
  key_dir <- file.path(OUT_DIR, "key_outputs")
  dir.create(key_dir, recursive = TRUE, showWarnings = FALSE)
  # `method` is the comparator KEY (what FIGURE_BASELINE_METHODS names it); `method_label` is
  # what the legend calls it, so the file is readable without the palette beside it.
  .a <- dplyr::mutate(curves, method_label = unname(MODEL_LBL[as.character(method)]),
                      panel = "A_prioritisation")
  .b <- dplyr::mutate(dplyr::select(topk, dplyr::any_of(c("origin", "cutoff", "health_zone",
                                                          "p_invasion", "outcome", "rank"))),
                      panel = "B_topk_outcomes")
  readr::write_csv(dplyr::bind_rows(.a, .b),
                   file.path(key_dir, "Figure2_cascade_data.csv"))
  message("[fig2-cascade] wrote Figure2_cascade_data.csv (", nrow(.a), " curve rows + ",
          nrow(.b), " top-K rows)")
})
message(sprintf("[done] Figure2_cascade (%d origins, +%d wk) -> %s", N_ORIGINS, K_HORIZON, MANU_DIR))

}   # end of the script-run guard

