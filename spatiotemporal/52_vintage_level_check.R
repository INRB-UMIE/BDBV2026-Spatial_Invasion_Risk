# =============================================================================
# 52_vintage_level_check.R — do the two "what was known then" reconstructions agree?
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# WHY. A back-test round can rebuild its information state two ways (LFO_TRAINING_SOURCE,
# 00_config.R):
#   vintage   (the default) the ARCHIVED snapshot of the issue date, nowcast with the "extract"
#             truncation (onset -> appearance) FITTED ON THAT DATE from the snapshots then
#             available -- what the live forecast of that day would have used;
#   recensor  the FINAL line list censored to the issue date, nowcast with the shared "asof"
#             truncation (onset -> sample) -- the fallback where no snapshot can be used. That
#             truncation is ONE fit on the final line list up to the analysis date (shared across
#             folds, 02_epi_params.R), so this arm's agreement with the final counts is partly
#             in-sample; the snapshot arm's extract fit uses only the snapshots of its own date.
# Under "vintage" a round's early rolling-predictor transitions (before the archive starts,
# 2026-06-09, or before it supports the truncation fit) still re-censor, so one design mixes
# the two. That is harmless only if both
# reconstructions put the import force at the same LEVEL, and the import force is what the
# invasion hazard is linear in. This module measures that at every post-burn-in origin, and
# how close each nowcast of the origin week comes to the count that week eventually reached.
#
# WHAT IT COMPARES, per origin (cutoff C, issue date C + 6):
#   * confirmed cases known, cumulative and in the origin week (raw), and how many zones had one;
#   * the nowcast of the origin week and of the last three weeks, against the FINAL counts for the
#     same weeks (the eventual truth both nowcasts estimate; final data also carry later
#     corrections, so neither is expected to hit it exactly);
#   * the import force Lambda of the FEATURED kernel and generation time, summed over the zones
#     at risk at C on the final data (the rows the round is scored on), and the Spearman
#     correlation of Lambda across those zones (does the ranking survive the change of source?).
#
# OUTPUT: diagnostics/vintage_level_check.csv and diagnostics/vintage_level_check.md.
# Runtime: one load_linelist() per snapshot plus, for each origin, a vintage_truncation() fit
# (an EpiNow2/Stan estimate_truncation run whenever its disk cache misses; cached fits are reused).
# =============================================================================

ST_DIR <- file.path(here::here(), "spatiotemporal")
source(file.path(ST_DIR, "00_config.R"))
load_pipeline_modules("model", dir = ST_DIR)
suppressPackageStartupMessages({ library(dplyr) })

if (!script_run_or_forced("52_vintage_level_check.R")) {
  message("[vintage-check] sourced, not run: nothing is written. Use `Rscript ",
          "52_vintage_level_check.R` (or FORCE_RUN=1) to run it.")
} else {

.tag <- function(...) message("[vintage-check] ", sprintf(...))

# ---- 1. featured kernel and generation time ------------------------------------------------------
sel_path <- file.path(OUT_DIR, "key_outputs", "model_selection.json")
if (!file.exists(sel_path)) stop("[vintage-check] model_selection.json not found; run run_all.R first.")
.sel     <- jsonlite::fromJSON(sel_path)
# read_featured_model()'s rule (figure_style.R, not loaded here): featured$bayesian, else headline.
featured <- .sel$featured$bayesian$method
if (is.null(featured) || !length(featured) || is.na(featured[1])) featured <- .sel$featured$headline$method
kernel   <- mobility_kernel_from_method(featured)
if (!is.character(featured) || is.na(featured) || is.na(kernel))
  stop("[vintage-check] could not resolve the featured model or its kernel from model_selection.json.")
mm <- list(); mm[[kernel]] <- readRDS(file.path(OUT_MOBILITY, sprintf("mobility_%s.rds", kernel)))
# A featured time-varying form resolves to its own spec (its stage-1 kernel and generation time;
# the form does not enter the import force this check compares).
.grid <- bayes_default_grid(mm)
.s1   <- bayes_spec_for_label(.sel$featured$bayesian$stage1_method %||% NA_character_, .grid)
spec  <- bayes_spec_for_label(featured, c(.grid, bayes_tv_forms_of(.s1)))
if (is.null(spec))
  stop(sprintf("[vintage-check] the featured model %s is not a spec of bayes_default_grid() or a form of its stage-1 pick.", featured))

# ---- 2. data and both truncations, exactly as run_all.R resolves them ----------------------------
dat       <- prep_all_data()
zones_all <- names(dat$pop)
gt_pmfs   <- compute_all_gt_pmfs()
G <- daily_to_weekly_gt(gt_pmfs[[spec$gt]]); W <- mm[[kernel]]
trunc_asof <- .trunc_as_delay_spec(epinow2_truncation(regime = "asof", ll = dat$ll), "asof")
if (is.null(trunc_asof))
  stop("[vintage-check] the as-of truncation is unavailable; the two sources cannot be compared.")
.tag("featured %s (kernel %s, GT %s); as-of truncation mean %.2f d; each snapshot is nowcast with the extract truncation fitted on its own date",
     featured, kernel, spec$gt, trunc_asof$mean)

zwo <- dat$zone_week %>%
  dplyr::filter(week_start >= lubridate::floor_date(OUTBREAK_START, "week",
                                                    week_start = get0("WEEK_ANCHOR", ifnotfound = 1L)))
weeks <- sort(unique(zwo$week_start))
cuts  <- vapply(invasion_fold_candidates(weeks, LFO_HORIZONS, ANALYSIS_DATE,
                                         min_train_weeks = LFO_MIN_TRAINING_WEEKS,
                                         min_eval_age_days = LFO_MIN_EVAL_AGE_DAYS),
                function(f) as.numeric(f$cutoff), 0)
cuts  <- as.Date(cuts, origin = "1970-01-01")
cuts  <- cuts[!invasion_is_burnin(cuts)]
final_Y <- .count_wide(zwo, zones_all, "confirmed")

lam_of <- function(nc) {
  Y <- .count_wide(nc, zones_all, "confirmed_nc")
  compute_foi(Y, W, G, t_idx = ncol(Y) + 1L, zones_all)
}
last_k <- function(Y, k) { j <- seq.int(max(1L, ncol(Y) - k + 1L), ncol(Y)); sum(Y[, j, drop = FALSE]) }

# ---- 3. per origin -------------------------------------------------------------------------------
rows <- list()
for (i in seq_along(cuts)) {
  C <- cuts[i]; issue <- C + 6L; spine <- weeks[weeks <= C]
  v <- vintage_for_issue(issue, source = "vintage")
  if (is.null(v)) { .tag("issue %s: no snapshot within %d d; skipped", format(issue), LFO_VINTAGE_MAX_GAP_DAYS); next }
  vll <- load_linelist_vintage(v$folder, v$date)
  # The truncation the live forecast of that day could have estimated (NULL before the archive
  # supports the fit): exactly what a snapshot round of the back-test is nowcast with.
  tr <- vintage_truncation(v$date)
  r_raw <- suppressMessages(reaggregate_asof(dat$ll, zones_all, issue, week_spine = spine)) %>%
    dplyr::filter(week_start <= C)
  # A snapshot is read as it stands, NOT re-censored -- exactly as the back-test reads it
  # (asof_counts_source(), 22_daily_reissue.R).
  v_raw <- suppressMessages(reaggregate_asof(vll, zones_all, issue, week_spine = spine,
                                             censor = FALSE)) %>%
    dplyr::filter(week_start <= C)
  r_nc <- suppressMessages(suppressWarnings(apply_nowcast_correction(r_raw, analysis_date = issue, delay = trunc_asof)))
  v_nc <- if (is.null(tr)) NULL else
    suppressMessages(suppressWarnings(apply_nowcast_correction(v_raw, analysis_date = v$date, delay = tr$spec)))
  Yr <- .count_wide(r_raw, zones_all, "confirmed"); Yv <- .count_wide(v_raw, zones_all, "confirmed")
  Nr <- .count_wide(r_nc, zones_all, "confirmed_nc")
  Nv <- if (is.null(v_nc)) NULL else .count_wide(v_nc, zones_all, "confirmed_nc")
  Yf <- final_Y[, colnames(final_Y) <= format(C), drop = FALSE]
  atrisk <- rowSums(Yf > 0) == 0
  Lr <- lam_of(r_nc); Lv <- if (is.null(v_nc)) rep(NA_real_, length(Lr)) else lam_of(v_nc)
  .nv <- function(k) if (is.null(Nv)) NA_real_ else last_k(Nv, k)
  rows[[length(rows) + 1L]] <- data.frame(
    cutoff = C, issue_date = issue, vintage_date = v$date, gap_days = v$gap,
    vintage_truncation_mean = if (is.null(tr)) NA_real_ else tr$spec$mean,
    conf_cum_recensor = sum(Yr), conf_cum_vintage = sum(Yv), conf_cum_final = sum(Yf),
    zones_affected_recensor = sum(rowSums(Yr) > 0), zones_affected_vintage = sum(rowSums(Yv) > 0),
    zones_affected_final = sum(!atrisk),
    origin_week_raw_recensor = last_k(Yr, 1), origin_week_raw_vintage = last_k(Yv, 1),
    origin_week_nc_recensor = last_k(Nr, 1), origin_week_nc_vintage = .nv(1),
    origin_week_final = last_k(Yf, 1),
    last3_nc_recensor = last_k(Nr, 3), last3_nc_vintage = .nv(3), last3_final = last_k(Yf, 3),
    n_atrisk = sum(atrisk),
    sum_lambda_recensor = sum(Lr[atrisk]), sum_lambda_vintage = sum(Lv[atrisk]),
    lambda_ratio_vintage_over_recensor = sum(Lv[atrisk]) / sum(Lr[atrisk]),
    lambda_spearman = suppressWarnings(stats::cor(Lr[atrisk], Lv[atrisk], method = "spearman")),
    stringsAsFactors = FALSE)
}
res <- dplyr::bind_rows(rows)
if (!nrow(res)) stop("[vintage-check] no post-burn-in origin has a snapshot; nothing to compare.")
res <- res %>% dplyr::mutate(
  # COMPLETENESS OF THE ORIGIN WEEK: what each source held at the origin as a share of what the
  # week eventually reached (observed), against the share its truncation model ASSUMED (raw over
  # nowcast). observed / assumed > 1 means the truncation over-corrects that week by that factor.
  origin_week_completeness_observed_recensor = origin_week_raw_recensor / origin_week_final,
  origin_week_completeness_assumed_recensor  = origin_week_raw_recensor / origin_week_nc_recensor,
  origin_week_completeness_observed_vintage  = origin_week_raw_vintage / origin_week_final,
  origin_week_completeness_assumed_vintage   = origin_week_raw_vintage / origin_week_nc_vintage,
  origin_week_nc_err_recensor = origin_week_nc_recensor / origin_week_final - 1,
  origin_week_nc_err_vintage  = origin_week_nc_vintage  / origin_week_final - 1,
  last3_nc_err_recensor = last3_nc_recensor / last3_final - 1,
  last3_nc_err_vintage  = last3_nc_vintage  / last3_final - 1)

dir.create(OUT_DIAGNOSTICS, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, file.path(OUT_DIAGNOSTICS, "vintage_level_check.csv"), row.names = FALSE)

# ---- 4. summary ----------------------------------------------------------------------------------
# BOTH ARMS ON ONE SET OF ORIGINS: those whose snapshot carries an extract-truncation fit (the
# only ones with a vintage nowcast). Taking the recensor medians over every origin and the
# vintage medians over the fitted subset compared the two arms on different origin sets.
gm  <- function(x) exp(mean(log(x[is.finite(x) & x > 0])))
med_abs <- function(x) sprintf("%.0f%%", 100 * stats::median(abs(x), na.rm = TRUE))
rf <- res[is.finite(res$vintage_truncation_mean), , drop = FALSE]
if (!nrow(rf)) stop("[vintage-check] no origin has a snapshot truncation fit; the two arms cannot be compared.")
s <- list(
  n = nrow(res), n_fit = nrow(rf),
  ratio_gm = gm(rf$lambda_ratio_vintage_over_recensor),
  ratio_range = range(rf$lambda_ratio_vintage_over_recensor, na.rm = TRUE),
  ratio_trend = suppressWarnings(stats::cor(seq_len(nrow(rf)), log(rf$lambda_ratio_vintage_over_recensor),
                                            method = "spearman", use = "complete.obs")),
  fit_range = range(rf$vintage_truncation_mean),
  rho_min = min(rf$lambda_spearman, na.rm = TRUE), rho_med = stats::median(rf$lambda_spearman, na.rm = TRUE))
md <- c(
  sprintf("# Vintage vs re-censored information state (%s)", format(ANALYSIS_DATE)), "",
  sprintf(paste0("_Featured_: `%s` (kernel %s, generation time %s). %d post-burn-in origins, each with a ",
                 "snapshot %s day(s) before it. Truncation: re-censored data use the shared as-of fit (mean %.2f d), ",
                 "made ONCE on the final line list up to the analysis date, so the re-censored arm's agreement ",
                 "with the final counts is partly in-sample; each snapshot uses the extract fit made on its own ",
                 "date from the snapshots then available (%d of %d origins have one, means %.2f-%.2f d; the ",
                 "others predate enough archive for the fit and re-censor in the back-test). Every summary ",
                 "line below is on those %d origins; the table lists every origin."),
          featured, kernel, spec$gt, s$n, paste(sort(unique(res$gap_days)), collapse = "/"),
          trunc_asof$mean, s$n_fit, s$n, s$fit_range[1], s$fit_range[2], s$n_fit), "",
  sprintf(paste0("**Import-force level, vintage / recensor** (summed over the zones at risk): geometric mean %.2f, ",
                 "range %.2f-%.2f; Spearman trend over origins %+.2f."),
          s$ratio_gm, s$ratio_range[1], s$ratio_range[2], s$ratio_trend),
  sprintf("**Ranking agreement** (Spearman of Lambda across at-risk zones): median %.3f, minimum %.3f.",
          s$rho_med, s$rho_min), "",
  sprintf(paste0("**Nowcast of the origin week vs its final count** (median absolute relative error, ",
                 "|nowcast / final - 1|): recensor %s, vintage %s. Last three weeks: recensor %s, vintage %s."),
          med_abs(rf$origin_week_nc_err_recensor), med_abs(rf$origin_week_nc_err_vintage),
          med_abs(rf$last3_nc_err_recensor), med_abs(rf$last3_nc_err_vintage)),
  sprintf(paste0("**Origin-week completeness, observed vs assumed by the truncation** (median over origins): ",
                 "recensor %.2f vs %.2f; snapshot %.2f vs %.2f (last four origins %.2f vs %.2f). The extract ",
                 "estimator is the one the LIVE nowcast uses; observed/assumed above 1 means it over-corrects ",
                 "the week the forecast is issued from."),
          stats::median(rf$origin_week_completeness_observed_recensor, na.rm = TRUE),
          stats::median(rf$origin_week_completeness_assumed_recensor, na.rm = TRUE),
          stats::median(rf$origin_week_completeness_observed_vintage, na.rm = TRUE),
          stats::median(rf$origin_week_completeness_assumed_vintage, na.rm = TRUE),
          stats::median(utils::tail(rf$origin_week_completeness_observed_vintage, 4), na.rm = TRUE),
          stats::median(utils::tail(rf$origin_week_completeness_assumed_vintage, 4), na.rm = TRUE)), "",
  "| Origin | Snapshot | Confirmed known (re-censored / snapshot / final) | Origin week raw (r / v) | Origin week nowcast (r / v / final) | Sum Lambda ratio v/r | Spearman |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %d / %d / %d | %d / %d | %.0f / %.0f / %d | %.2f | %.3f |",
          format(res$issue_date), format(res$vintage_date), res$conf_cum_recensor, res$conf_cum_vintage,
          res$conf_cum_final, res$origin_week_raw_recensor, res$origin_week_raw_vintage,
          res$origin_week_nc_recensor, res$origin_week_nc_vintage, res$origin_week_final,
          res$lambda_ratio_vintage_over_recensor, res$lambda_spearman), "",
  paste0("_A level ratio far from 1, or one that trends over origins, means the two sources put the import ",
         "force at different scales; under LFO_TRAINING_SOURCE = \"vintage\" the rolling predictors before ",
         "the archive starts would then sit on a different scale from the later ones. The final counts include ",
         "later corrections, so neither nowcast is expected to match them exactly._"))
writeLines(md, file.path(OUT_DIAGNOSTICS, "vintage_level_check.md"))
.tag("%d origins (%d with a snapshot truncation fit); on those, Lambda vintage/recensor geometric mean %.2f (range %.2f-%.2f, trend %+.2f); Spearman median %.3f",
     s$n, s$n_fit, s$ratio_gm, s$ratio_range[1], s$ratio_range[2], s$ratio_trend, s$rho_med)
.tag("wrote diagnostics/vintage_level_check.{csv,md}")

}   # end of the script-run guard
