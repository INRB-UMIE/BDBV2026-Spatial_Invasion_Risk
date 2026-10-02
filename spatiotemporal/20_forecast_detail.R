# =============================================================================
# 20_forecast_detail.R — Best-model documentation, ensemble maps, and
#                        space-time forecast visualisation WITH uncertainty
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Answers four follow-up requests, all additive (nothing in 15/16/17 is changed):
#   (Q2) how the best model is selected            -> describe_selection()
#   (Q3) what goes into the best model             -> describe_model_spec() and
#        (structure, covariates, GT, observation,     write_model_details_report()
#        calibration, nowcast)
#   (Q1) visualise the ensemble forecast           -> plot_forecast_map_panel()
#   (Q4) visualise forecasts across space & time   -> ensemble_member_uncertainty(),
#        WITH uncertainty                              plot_forecast_uncertainty(),
#                                                      plot_spacetime_forecast()
#
# Forecast uncertainty is taken as the DISAGREEMENT ACROSS the pre-specified
# ensemble members (min / median / max of the per-zone invasion probability).
# For a rare event with ~15 training signals this structural spread is the
# honest, leakage-free uncertainty band; it is wider where the mobility / GT /
# observation choice matters most and narrow where the members concur.
# =============================================================================

source(file.path(here::here(), "spatiotemporal", "00_config.R"))
suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(patchwork); library(scales)
})

# NOTE: base R >= 4.4 ships `%||%`, so `exists("%||%")` is TRUE and the NA-aware definition
# below NEVER installed — `NA %||% "fallback"` returned NA. That is why .lk() exists. Give the
# NA-aware helper its OWN name so it is always available, and keep %||% as base's NULL-coalesce.
if (!exists("%||%"))
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
#' NULL/NA/length-safe coalesce: returns `b` unless `a` is a single non-NA value.
`%|NA|%` <- function(a, b) if (is.null(a) || length(a) != 1L || is.na(a)) b else a
if (!exists("OKABE_ITO"))
  OKABE_ITO <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                 "#56B4E9", "#F0E442", "#000000")
if (!exists("theme_inv"))
  theme_inv <- function(base = 12) ggplot2::theme_minimal(base_size = base)

# Stable colour map for provinces (health areas) so a province reads the SAME
# across every figure. Known DRC provinces get fixed colourblind-safe hues;
# any others are filled from a fallback ramp. Returns a named vector for the
# provinces actually present.
.prov_colours <- function(provs) {
  provs <- sort(unique(as.character(provs)))
  provs[is.na(provs) | provs == ""] <- "Unknown"
  provs <- sort(unique(provs))
  known <- c("Ituri" = "#D55E00", "Nord-Kivu" = "#0072B2", "Sud-Kivu" = "#009E73",
             "Haut-Uele" = "#CC79A7", "Tshopo" = "#E69F00", "Bas-Uele" = "#56B4E9",
             "Maniema" = "#F0E442", "Haut-Katanga" = "#332288", "Tshuapa" = "#117733",
             "Unknown" = "#999999")
  pal <- stats::setNames(rep(NA_character_, length(provs)), provs)
  hit <- intersect(provs, names(known)); pal[hit] <- known[hit]
  miss <- provs[is.na(pal)]
  if (length(miss)) {
    fb <- grDevices::hcl.colors(max(length(miss), 2), "Dark3")
    pal[miss] <- fb[seq_along(miss)]
  }
  pal
}
# Province label for a row: the province, or "Unknown" if missing.
.prov_label <- function(province) ifelse(is.na(province) | province == "", "Unknown", province)
.fd_save <- function(p, path, w = 9, h = 6) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !.fk(path)) return(invisible(p))
  if (file.exists(path)) suppressWarnings(file.remove(path))
  tryCatch(ggplot2::ggsave(path, p, width = w, height = h, device = "pdf",
                           limitsize = FALSE),
           error = function(e) warning("[viz] ", basename(path), ": ",
                                       conditionMessage(e)))
  if (file.exists(path)) message("[viz] saved -> ", basename(path))
  invisible(p)
}

# Raster twin of .fd_save, for figures that also need a bitmap (slides, Word, the
# manuscript submission portal). 600 dpi on an opaque white ground, matching the
# PNG convention used by the key-output figure modules.
.fd_save_png <- function(p, path, w = 9, h = 6, dpi = 600) {
  # Retained-figure gate (FIGURE_KEEP, 00_config.R): silently skip any figure that
  # is not on the published allow-list. get0() so the helper still works standalone.
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !.fk(path)) return(invisible(p))
  if (file.exists(path)) suppressWarnings(file.remove(path))
  tryCatch(ggplot2::ggsave(path, p, width = w, height = h, dpi = dpi, bg = "white",
                           limitsize = FALSE),
           error = function(e) warning("[viz] ", basename(path), ": ",
                                       conditionMessage(e)))
  if (file.exists(path)) message("[viz] saved -> ", basename(path))
  invisible(p)
}

# ---------------------------------------------------------------------------
# Q2/Q3 — describe the model that is selected, and how it is selected
# ---------------------------------------------------------------------------

# Human-readable descriptions of each structural component the method names
# encode. New variants fall back to a name-parsed default, so nothing breaks if
# the model grid grows.
.MOB_DESC <- c(
  M8  = "M8 composite — outbreak-specific short-trip flows from the epicentre cluster (Flowminder) for epicentre origins, calibrated gravity elsewhere (recommended)",
  M4  = "M4 gravity — negative-binomial gravity GLM fitted to the national Flowminder relocation matrix (monthly home-location changes, not trips; power-law distance deterrence)",
  M4b = "M4b gravity — as M4 but with exponential distance deterrence",
  M9  = "M9 multi-kernel — average of gravity, radiation and travel-time-decay kernels",
  M10 = "M10 radiation-composite — parameter-free radiation model blended with epicentre short-trip flows",
  M3  = "M3 — raw Flowminder national relocation flows (monthly home-location changes, not trips)",
  M5  = "M5 radiation model (Simini et al. 2012)",
  M6a = "M6a travel-time exponential-decay kernel",
  M6b = "M6b travel-time power-law-decay kernel",
  M13 = "M13 cohort-composite — Flowminder cohort presence rows for the cohort origins, calibrated gravity elsewhere",
  M14 = "M14 cohort-radiation composite — Flowminder cohort presence rows, parameter-free radiation elsewhere",
  M15 = "M15 symmetrised relocation OD — O + t(O) from one directed relocation table (not inflow-informed)",
  M16 = "M16 cohort + relocation OD — cohort presence rows, directed Flowminder relocation flows (coverage gaps filled) elsewhere",
  M17 = "M17 all-kernel consensus — equal-weight mean of the relocation OD (gap-filled), gravity and radiation kernels, cohort rows overlaid")
.GT_DESC <- c(
  medium = "GT-Medium (mean 15.3 d, sd 9.3 d; Zaire EBOV serial interval, WHO Ebola Response Team 2014 NEJM)",
  short  = "GT-Short (mean 12.0 d, sd 6.5 d; low-end sensitivity, below Nash 2024 pooled 95% CI)",
  long   = "GT-Long (mean 18.0 d, sd 10.5 d; high-end sensitivity, above Nash 2024 pooled 95% CI)")
.COV_DESC <- c(
  log_pop      = "log population size",
  ccvi         = "CCVI socioeconomic vulnerability index",
  positivity   = "zone test positivity",
  d_min        = "travel time to the nearest currently-affected zone",
  week_idx     = "calendar week index (log-linear time trend in the import coefficient)",
  sd_week      = "SD of the weekly random intercept (week-to-week variation in beta)",
  alert_import = "mobility-weighted suspected-alert pressure imported from other zones",
  alert_local  = "local suspected-alert count",
  susp_import  = "mobility-weighted suspected-but-not-confirmed case pressure imported from other zones",
  susp_local   = "local preceding-week suspected-but-not-confirmed case count",
  # KEYED "healthsite_density", the name every model and the covariate screen actually use
  # (run_all.R:541/550, the -full covariate set). The key was "log_healthsite_count", which
  # matches nothing, so .lk() fell through to the raw variable name and the published model
  # specification listed three covariates in prose and one as a bare identifier.
  healthsite_density = "healthcare-site density (facilities per capita)")
# Method-name -> covariate set / special structure (mirrors the run_all grid).
.VARIANT_SPEC <- list(
  "Renewal-M8-geo"      = list(cov = c("log_pop", "ccvi", "d_min")),
  "Renewal-M8-alert"    = list(cov = c("alert_import", "alert_local")),
  "Renewal-M8-susp"     = list(cov = c("susp_import", "susp_local")),
  # FULL-exogenous covariate set (matches the run_all Renewal-M8-cov definition and the
  # Bayesian "full" grid): log population, CCVI, health-site density, travel time to the
  # nearest affected zone. positivity is excluded (circular + full-data static leak).
  "Renewal-M8-cov"      = list(cov = c("log_pop", "ccvi", "healthsite_density", "d_min")),
  "Renewal-M8-cov-susp" = list(cov = c("log_pop", "ccvi", "healthsite_density", "d_min",
                                       "susp_import", "susp_local")),
  "Renewal-M8-report" = list(note = "ascertainment/reporting-rate structure using health-site density as a completeness proxy"),
  "Renewal-M8-raw"    = list(nowcast = "raw (no nowcast correction)"),
  "Renewal-M8-cwt"    = list(note = "calibration weighted by each training week's reporting completeness"))

#' Structured specification of a model, parsed from its method name.
#' @return named list(family, mobility, gt, observation, covariates, calibration,
#'   nowcast, notes) of human-readable strings.
#' Name-safe lookup on a NAMED ATOMIC vector: `x[[key]]` ERRORS on a missing name
#' (so `%||%` can never supply the fallback); this returns `fallback` instead —
#' essential now the model grid can grow beyond the hard-coded description keys.
.lk <- function(x, key, fallback)
  if (!is.null(key) && length(key) == 1 && !is.na(key) && key %in% names(x)) x[[key]] else fallback

describe_model_spec <- function(method) {
  m <- as.character(method)
  # GUARD the empty case. as.character(NULL) is character(0), so every grepl() below returns
  # logical(0) and the first `if (is_ens)` threw "argument is of length zero". run_all.R passes
  # primary_method = NULL whenever the frequentist arm did not run — which is the DEFAULT — so
  # write_model_details_report() aborted on the default path. The abort is swallowed by vwrap(),
  # with the result that outputs/reports/model_specification.md was never regenerated (it sat
  # six weeks stale, naming a model the pipeline no longer selects) and the selection section
  # was silently missing from the current invasion report.
  if (!length(m) || is.na(m[1]) || !nzchar(m[1]))
    return(list(kernel = NA_character_, generation_time = NA_character_,
                observation = NA_character_, covariates = NA_character_,
                calibration = NA_character_, nowcast = NA_character_,
                notes = "no model of this family was fitted in this run"))
  is_renewal <- grepl("^Renewal", m)
  is_ens     <- grepl("^Ensemble", m)

  mob_key <- stringr::str_match(m, "-(M[0-9]+[a-z]?)")[, 2]
  gt_key  <- dplyr::case_when(grepl("-short", m) ~ "short",
                              grepl("-long", m)  ~ "long",
                              TRUE ~ "medium")
  obs     <- if (grepl("-NB", m)) "Negative-binomial (overdispersion from renewal residuals)"
             else "Poisson (arrival-hazard limit)"
  vs      <- .VARIANT_SPEC[[m]]

  if (is_ens) {
    return(list(
      family = "Ensemble", mobility = "several (member matrices)",
      gt = "several (member profiles)",
      observation = "combined across members",
      covariates = "n/a",
      calibration = if (grepl("mean", m))
        "linear opinion pool: p = mean over members (Vincent-averaged count quantiles)"
        else "median of member probabilities (robust combiner)",
      # The DEPLOYED nowcast is apply_nowcast_correction() (run_all.R:321), the deterministic
      # delay-CDF right-truncation correction, on BOTH the fold and the deployment paths.
      # epinowcast was retained only as the sensitivity arm on 2026-09-17; naming it here
      # published a false "Incidence input / nowcast" line for every model.
      nowcast = "right-truncation-corrected members (deterministic delay-CDF nowcast)",
      notes = "pre-specified member set; no selection on the test folds"))
  }
  # Bayesian renewal grid (Bayes-<mob>-<gt/cov/link>): same mechanistic model as the
  # frequentist renewal, fit with brms — so it must be DESCRIBED, not dumped into the
  # "Comparator / NA" bucket (the featured Bayesian model routinely wins overall).
  if (grepl("^Bayes", m)) {
    # Keyed on the model's OWN kernel id so a filled / split / road-distance composite is
    # never described as its unfilled travel-time parent.
    kern <- mobility_kernel_from_method(m)
    base <- if (is.na(kern)) NA_character_ else gsub("-(dist|fill|split)", "", kern)
    mobn <- c(M1 = "short-trip", M4 = "gravity", M8 = "composite-gravity",
              M9 = "multi-kernel ensemble", M10 = "radiation-composite",
              M11 = "inward meeting-location FOI", M13 = "cohort + gravity composite",
              M14 = "cohort + radiation composite", M15 = "symmetrised relocation OD",
              M16 = "cohort + relocation OD composite",
              M17 = "all-kernel consensus ensemble")
    cov_desc <- if (grepl("-full", m)) "full exogenous set (log-pop, CCVI, health-site density, travel-time-to-affected)"
                else if (grepl("-geo", m)) "geographic exogenous set (log-pop, CCVI, travel-time-to-affected)"
                else "none (intercept-only import hazard)"
    lnk <- if (grepl("-logit", m)) "Logit (log-odds; beta is not a hazard multiplier)"
           else "Complementary log-log with log(Lambda) offset (renewal hazard: p = 1 - exp(-beta*Lambda))"
    return(list(
      family = "Bayesian mobility-informed renewal (brms; posterior)",
      # Keyed on the model's own kernel: the base family, then every kernel qualifier
      # (road-distance deterrence, source-cell fill, cohort origin-split). Without the
      # fill/split notes a filled composite was described exactly like its unfilled parent.
      mobility = {
        .k <- if (is.na(base)) mob_key else base
        .mob_base <- .lk(mobn, .k, .k %|NA|% "n/a")
        paste0(.mob_base,
               if (!is.na(kern) && grepl("-dist", kern)) " (OSRM road-distance deterrence)" else "",
               if (!is.na(kern) && grepl("-split", kern)) "; cohort rows split per origin" else "",
               if (!is.na(kern) && grepl("-(fill|split)", kern)) "; source-cell fill" else "")
      },
      gt = gt_key, observation = lnk, covariates = cov_desc,
      calibration = "posterior; weakly-informative priors (Intercept~Normal(-3,2), coefs~Normal(0,1)); loo predictive stacking for the ensemble",
      # The DEPLOYED nowcast is apply_nowcast_correction() (run_all.R:321), the deterministic
      # delay-CDF right-truncation correction, on BOTH the fold and the deployment paths.
      # epinowcast was retained only as the sensitivity arm on 2026-09-17; naming it here
      # published a false "Incidence input / nowcast" line for every model.
      nowcast = "right-truncation-corrected training counts (deterministic delay-CDF nowcast)",
      notes = paste("Bayesian analogue of the frequentist renewal model; the featured Bayesian",
                    "model is chosen by the leave-future-out CV composite (NOT by the loo",
                    "predictive-stacking weight, which defines the stacked ensemble only)")))
  }
  if (!is_renewal) {
    cmp <- c(hhh4 = "Endemic-epidemic hhh4 (Held & Paul); neighbourhood-coupled autoregression",
             `Gravity-B4` = "gravity-decay invasion baseline (no renewal dynamics)",
             `Distance-B1` = "distance-weighted neighbour-average baseline",
             `SEIR-C6` = "stochastic spatial SEIR metapopulation")
    return(list(family = "Comparator", mobility = NA, gt = NA, observation = NA,
                covariates = NA, calibration = NA, nowcast = NA,
                notes = .lk(cmp, m, m)))
  }

  covs <- vs$cov
  list(
    family = "Mobility-informed renewal (invasion model)",
    mobility = if (!is.na(mob_key)) .lk(.MOB_DESC, mob_key, mob_key) else "M8 composite",
    gt = .lk(.GT_DESC, gt_key, gt_key),
    observation = obs,
    covariates = if (is.null(covs)) "none (import force only)"
                 else paste(vapply(covs, function(c) .lk(.COV_DESC, c, c), ""),
                            collapse = "; "),
    calibration = "import coefficient beta fit by complementary-log-log regression of realised invasions on log(import force) => p = 1 - exp(-beta * Lambda), calibrated to the rare base rate",
    # See the note above: the deployed nowcast is the deterministic delay-CDF correction.
    nowcast = vs$nowcast %||% "deterministic delay-CDF right-truncation correction",
    notes = vs$note %||% NA_character_)
}

#' Markdown paragraph describing HOW the best model is chosen.
describe_selection <- function() {
  # KEEP IN STEP WITH best_invasion_model() (16_invasion_eval.R). This text is appended to
  # the invasion report and written to model_specification.md, and it described the RETIRED
  # "highest total AUC-PR skill" rule long after selection moved to the CV composite — so the
  # published methods statement named a different criterion from the one that picked the
  # featured model (and, through it, the cascade kernel). It also quoted a fixed "~15 invasion
  # events", which the accruing record left far behind; the count is deliberately no longer
  # hard-coded here (invasion_evaluation.csv carries it per horizon).
  .axis <- if (isTRUE(get0("INVASION_SELECT_ON_RECAL", ifnotfound = FALSE)))
    "the RECALIBRATED log score (16b removes calibration-in-the-large, so this axis measures refinement)"
  else "the raw log score"
  c("## How the models are selected", "",
    "Every model is scored in the **same** leave-future-out cross-validation on the",
    "at-risk zones only (no already-affected zone ever counts), pooled across folds.",
    "Selection is by a **three-axis CV composite with a spiky-model gate** (not by any single",
    "metric, and never by count-WIS, which would reward predicting zero everywhere):",
    "",
    "1. **AUC-PR skill** — Average Precision divided by the invasion base rate:",
    "   how well the model concentrates the few invasions at the top of its ranking;",
    "2. **Mean rank of the invaded zone** — how near the top the zones that actually",
    "   invaded were placed (operational targeting);",
    "3. **Log-score** — a proper score that penalises over-confident probabilities;",
    sprintf("   this run scores %s.", .axis),
    "",
    "Among models present at BOTH horizons, any whose worst-horizon mean rank-of-truth",
    "exceeds 1.5x the field median is first dropped as **spiky** (high AUC-PR driven by a",
    "few top hits while the rest are scattered), as are models cross-validated on fewer",
    "folds than the field (an incomparable denominator). Within each horizon the survivors",
    "are then RANKED on each of the three axes; the three ranks are summed within horizon",
    "and across horizons, and the **lowest total wins**. Ties are broken by higher total",
    "AUC-PR skill, then lower total rank-of-truth, then lower log score. Requiring both",
    "horizons stops a model winning by performing well at 1-week while sitting out the",
    "2-week task.",
    "",
    "The **best-fitting model** is the best over ALL families on this composite;",
    "the **best renewal model** is the best *mobility-informed renewal* variant (so",
    "the featured risk product is an interpretable, mechanistic model).",
    "With only a few tens of invasion events, differences between the leading models are",
    "small relative to cross-validation noise — hence the pre-specified ensemble is reported",
    "alongside the single best model, and the last origins are held out for a separate",
    "optimism check.", "")
}

#' Append a "best-model specification" + "selection" section to a report file
#' (default the invasion report), and also write a standalone model_specification.md.
write_model_details_report <- function(best_method, primary_method, eval_tbl,
                                       report_path = file.path(OUT_REPORTS, "invasion_report.md"),
                                       spec_path   = file.path(OUT_REPORTS, "model_specification.md"),
                                       bayes_params = NULL, best_bayes_method = NULL,
                                       freq_beta0 = NA_real_, freq_cov = NULL) {
  # Full NUMERIC parameterisation of the two featured models (#0: "all parameters").
  .params_md <- function() {
    out <- c("## Best-model parameters (all estimated quantities)", "",
      sprintf("### Frequentist — `%s`", primary_method), "",
      if (is.finite(freq_beta0))
        sprintf("- **Import coefficient beta0 = %.4g** — intercept-only cloglog calibration; P(first case)=1-exp(-beta0*Lambda) at mean covariates.", freq_beta0)
      else "- **Import coefficient beta0:** see `best_model_parameters.pdf`.")
    if (!is.null(freq_cov) && nrow(freq_cov)) {
      out <- c(out, "",
        "Covariate hazard ratios (per +1 SD; marginal HR with 95% CI, and status once the mobility import force is adjusted for):",
        "", "| Covariate | HR [95% CI] | p | Adjusted status |", "|---|---|---|---|")
      for (i in seq_len(nrow(freq_cov)))
        out <- c(out, sprintf("| %s | %.2f [%.2f-%.2f] | %.3f | %s |",
          (freq_cov$label[i] %||% freq_cov$term[i]), freq_cov$hr[i], freq_cov$lo[i],
          freq_cov$hi[i], freq_cov$p[i], freq_cov$adj_status[i] %||% "—"))
    }
    if (!is.null(bayes_params) && !is.null(best_bayes_method)) {
      bp <- bayes_params %>% dplyr::filter(model == best_bayes_method)
      if (nrow(bp)) {
        escale <- if ("effect_scale" %in% names(bp)) unique(bp$effect_scale)[1] else "hazard ratio"
        out <- c(out, "",
          sprintf("### Bayesian — `%s` (posterior median %s, 90%% CrI, Rhat)", best_bayes_method, escale), "",
          "| Parameter | Estimate [90% CrI] | P(effect>1) | Rhat |", "|---|---|---|---|")
        for (i in seq_len(nrow(bp))) {
          nm <- if (isTRUE(bp$is_intercept[i])) "Intercept (beta0)" else bp$term[i]
          out <- c(out, sprintf("| %s | %.3g [%.3g-%.3g] | %s | %.3f |", nm, bp$hr[i], bp$lo[i], bp$hi[i],
            if (is.na(bp$p_dir[i])) "-" else sprintf("%.2f", bp$p_dir[i]), bp$rhat[i]))
        }
      }
    }
    c(out, "")
  }
  spec_block <- function(label, method) {
    s <- describe_model_spec(method)
    e <- if (!is.null(eval_tbl)) eval_tbl %>% dplyr::filter(horizon == 1, method == !!method) else NULL
    perf <- if (!is.null(e) && nrow(e) == 1)
      sprintf("AUC-PR skill %.1fx · AUC-ROC %.3f · log-score %.3f · calibration %.1fx",
              e$auc_pr_skill, e$auc_roc, e$log_score, e$calibration_in_large) else "—"
    c(sprintf("### %s — `%s`", label, method), "",
      sprintf("- **Family / structure:** %s", s$family),
      sprintf("- **Mobility kernel W:** %s", s$mobility),
      sprintf("- **Generation-time weighting g(k):** %s", s$gt),
      sprintf("- **Observation process:** %s", s$observation),
      sprintf("- **Covariates on the invasion hazard:** %s", s$covariates),
      sprintf("- **Calibration:** %s", s$calibration),
      sprintf("- **Incidence input / nowcast:** %s", s$nowcast),
      if (!is.na(s$notes)) sprintf("- **Note:** %s", s$notes) else NULL,
      sprintf("- **Cross-validated skill (h=1):** %s", perf), "")
  }
  # A NULL primary_method means "no renewal model was fitted", not "a different model". Treat
  # it as absent rather than as a second featured model, otherwise the report advertises a
  # "Best renewal model" section for a family the run did not fit.
  .have_primary <- length(as.character(primary_method)) == 1L &&
                   !is.na(primary_method) && nzchar(as.character(primary_method))
  same <- .have_primary && identical(as.character(best_method), as.character(primary_method))
  blocks <- if (!.have_primary)
    c(paste0("The headline model (best over all families) is the featured model; no ",
             "frequentist renewal model was fitted in this run."), "",
      spec_block("Featured model", best_method))
  else if (same)
    c(paste0("The headline model (best over all families) is **also** the primary ",
             "mobility-informed renewal model:"), "",
      spec_block("Best-fitting model", best_method))
  else
    c("The two featured models are:", "",
      spec_block("Headline (best over all families)", best_method),
      spec_block("Best renewal model", primary_method))
  body <- c(
    "# BDBV 2026 — Model selection & best-model specification", "",
    describe_selection(),
    "## Best-model specification", "",
    "The invasion hazard is, for every mobility-informed renewal variant,",
    "`P(first case) = 1 - exp(-beta_i * Lambda_i)` where the import force",
    "`Lambda_i = sum_j W[j,i] * sum_k g(k) * Y_nc[j, t-k]` is the generation-time-",
    "weighted, mobility-routed incidence arriving at zone i, and `beta_i` (optionally",
    "modulated by covariates) is learned from realised invasions.", "",
    blocks, .params_md(), describe_priority_index())
  writeLines(body, spec_path)
  message("[report] model specification -> ", basename(spec_path))

  # Also append to the main invasion report if it exists — but idempotently, so
  # re-running does not accumulate duplicate sections.
  if (file.exists(report_path)) {
    existing <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
    if (!grepl("How the models are selected", existing, fixed = TRUE)) {
      appended <- if (same)
        c(describe_selection(), spec_block("Best-fitting model", best_method),
          describe_priority_index())
      else
        c(describe_selection(), spec_block("Headline model", best_method),
          spec_block("Best renewal model", primary_method), describe_priority_index())
      cat("\n\n", paste(appended, collapse = "\n"),
          file = report_path, append = TRUE, sep = "")
      message("[report] appended selection + spec to ", basename(report_path))
    } else {
      message("[report] selection + spec already present in ", basename(report_path), "; not re-appending")
    }
  }
  invisible(spec_path)
}

# ---------------------------------------------------------------------------
# Q4 — forecast uncertainty from ensemble-member disagreement
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# FREQUENTIST-ENSEMBLE VIZ — REMOVED
# ---------------------------------------------------------------------------
# ensemble_member_uncertainty(), plot_forecast_uncertainty() and plot_spacetime_forecast()
# visualised the MEMBER SPREAD of the frequentist renewal ensemble (min/median/max of
# p_invasion across pre-specified Renewal-M* members). That ensemble no longer exists —
# ENSEMBLE_MEMBERS was permanently empty and its call sites in run_all.R have been removed —
# so all three were unreachable. The Bayesian posterior 90% CrI maps are the equivalent
# product and are rendered by the Bayesian suite.

# Internal masked choropleth on the same key/mask convention as 17's map.
# `affected_keys` (lower/trimmed zone names) are shown as dark grey even when
# absent from `scored`, so every panel masks the SAME already-affected zones
# (the ensemble legitimately drops them, but the map must still grey them out).
.fd_map <- function(scored, shapefile, title, prob_col = "p_invasion",
                    ituri_only = FALSE, province_zoom = NULL, lim = NULL, palette = "plasma",
                    legend_name = "P(first case)", affected_keys = NULL) {
  if (is.null(shapefile)) return(NULL)
  if (!prob_col %in% names(scored))
    prob_col <- intersect(c("p_invasion", "p_case_invasion", "p_width"), names(scored))[1]
  # Include province in the join key when available: two health-zone NAMES are
  # duplicated nationally (Bili in Nord-Ubangi & Bas-Uele; Lubunga in Kasai-Central
  # & Tshopo), so a name-only join is many-to-many on national maps and mis-colours
  # one of each pair. Joining on (name, province) resolves it.
  has_prov <- "province" %in% names(scored)
  rd <- scored %>% dplyr::transmute(.key = tolower(trimws(health_zone)),
                                    .prov = if (has_prov) as.character(province) else NA_character_,
                                    p = .data[[prob_col]],
                                    affected = if ("was_active_before" %in% names(scored))
                                      was_active_before else FALSE)
  shp <- shapefile %>% dplyr::mutate(.key = tolower(trimws(Nom)), .prov = as.character(PROVINCE))
  if (ituri_only) shp <- shp %>% dplyr::filter(PROVINCE == "Ituri")
  else if (!is.null(province_zoom)) shp <- shp %>% dplyr::filter(PROVINCE == province_zoom)
  m <- (if (has_prov && !all(is.na(rd$.prov)))
          dplyr::left_join(shp, rd, by = c(".key", ".prov"))
        else dplyr::left_join(shp, dplyr::select(rd, -.prov), by = ".key")) %>%
    dplyr::mutate(affected = (affected %in% TRUE) | (.key %in% affected_keys))
  if (is.null(lim)) lim <- c(0, max(0.05, stats::quantile(rd$p, 0.99, na.rm = TRUE)))
  ggplot(m) +
    geom_sf(aes(fill = p), colour = "grey80", linewidth = 0.12) +
    geom_sf(data = dplyr::filter(m, affected %in% TRUE),
            fill = "grey55", colour = "grey80", linewidth = 0.12) +
    scale_fill_viridis_c(option = palette, name = legend_name, limits = lim,
                         oob = scales::squish, na.value = "grey90") +
    labs(title = title) + ggplot2::theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11), legend.position = "right")
}

# plot_forecast_map_panel() — REMOVED with the frequentist ensemble above: it drew the
# best-renewal-model vs ensemble side-by-side panel plus the member-disagreement map.

#' Build the invasion design matrix (one row per at-risk zone-week transition):
#' `invaded`, the mobility import offset `logLam`, and the STANDARDISED candidate
#' covariates — exactly the rows fit_import_model uses. Returns the design plus
#' the baseline import coefficient beta0 (intercept-only cloglog).
# ---------------------------------------------------------------------------
# ROLLING AS-OF PREDICTORS (2026-09-22)
# ---------------------------------------------------------------------------
# Session memo for the as-of count matrix. The 18 model closures each call
# build_invasion_design() on the SAME fold data, so without this the 90 per-transition
# reconstructions would be repeated once per model (1,620 calls instead of 90). Keyed on the
# issue date plus a hash of the line list and the delay, so a different fold, snapshot or
# truncation can never be served a stale matrix. parent = emptyenv() so a typo'd key cannot
# resolve up the search path.
.asof_memo <- new.env(parent = emptyenv())

#' Confirmed counts as they were KNOWN at the end of a given week, nowcast-corrected.
#' @return zones x weeks matrix, or NULL if the reconstruction fails.
.asof_counts_wide <- function(ll, zones_all, issue_date, delay, week_spine) {
  # Fingerprint the line list by row count AND the sum of its onset dates: two snapshots of
  # equal size would otherwise collide on nrow alone, and a collision here serves one fold's
  # counts to another. Cheap enough to recompute per call (one pass over a date column).
  .llfp <- c(nrow(ll), sum(as.numeric(suppressWarnings(as.Date(ll$date_of_symptom_onset))),
                           na.rm = TRUE))
  key <- paste0(format(as.Date(issue_date)), "_",
                substr(rlang::hash(list(.llfp, zones_all, delay, week_spine)), 1L, 12L))
  hit <- .asof_memo[[key]]
  if (!is.null(hit)) return(hit)
  tw <- tryCatch(suppressMessages(suppressWarnings(
          reaggregate_asof(ll, zones_all, as.Date(issue_date), week_spine = week_spine))),
        error = function(e) NULL)
  if (is.null(tw) || !nrow(tw)) return(NULL)
  nc <- tryCatch(suppressMessages(suppressWarnings(
          apply_nowcast_correction(tw, analysis_date = as.Date(issue_date), delay = delay))),
        error = function(e) NULL)
  if (is.null(nc)) return(NULL)
  out <- .count_wide(nc, zones_all, "confirmed_nc")
  assign(key, out, envir = .asof_memo)
  out
}

build_invasion_design <- function(zone_week_nc, mobility_matrices, gt_pmfs,
                                  covariates, osrm_mat, zones_all,
                                  mob = "M8", gt = "medium",
                                  # `positivity` REMOVED from the default (2026-09-17). Every
                                  # caller that matters already excluded it by passing an
                                  # explicit list (run_all.R:494, :1549) because it is circular
                                  # and full-data static; leaving it in the DEFAULT meant any
                                  # other caller silently reinstated it. It is additionally
                                  # NA-everywhere on real DHIS2 data — see the positivity guard
                                  # in 01_data_prep.R::aggregate_to_zone_week().
                                  candidates = c("alert_import", "alert_local",
                                                 "susp_import", "susp_local", "log_pop",
                                                 "healthsite_density", "ccvi",
                                                 "d_min"),
                                  # ROLLING AS-OF PREDICTORS. NULL (the default) keeps the
                                  # historic behaviour exactly: the import force for every
                                  # transition is computed from ONE count matrix, the fold's
                                  # as-of counts at its own origin. Supply a line list and the
                                  # regime's truncation to compute each transition's import
                                  # force from the counts that were KNOWN at that transition.
                                  # See the block above build_invasion_design() for why.
                                  rolling_ll = NULL, rolling_delay = NULL,
                                  rolling_floor = get0("ROLLING_PREDICTOR_FLOOR",
                                                       ifnotfound = 3L)) {
  Y_wide <- .count_wide(zone_week_nc, zones_all, "confirmed_nc")
  G <- daily_to_weekly_gt(gt_pmfs[[gt]]); W <- mobility_matrices[[mob]]
  A_wide <- .count_wide(zone_week_nc, zones_all, "total_alerts")
  # Suspected-but-not-confirmed leading-indicator matrix (nowcast-corrected suspected
  # series, raw fallback), for the susp_import / susp_local covariates — like the alert matrix.
  S_wide <- .susp_wide(zone_week_nc, zones_all)
  static <- .static_features(covariates, zones_all)
  # WHY ROLLING PREDICTORS EXIST. The import force is built from nowcast-corrected counts, and
  # the generation-time kernel puts 31% of its weight on the most recent week -- exactly the
  # week the nowcast inflates most (3.45x at a fold origin, 5.92x deployed). With one count
  # matrix per fold, every TRAINING transition's predictors are complete weeks (multiplier ~1)
  # while the FORECAST's predictor carries that inflation, so beta0 is calibrated against a
  # lambda about 2x smaller than the one it is applied to. Measured on the 12 folds: pooled
  # observed/expected 0.497 with the shared matrix, 0.930 with raw counts -- the two arms
  # differing by exactly the 1.978x inflation. That artefact is what the recalibration delta
  # (~0.4) was silently absorbing.
  #
  # THE FLOOR IS LOAD-BEARING, not a tidy-up. At small t the as-of reconstruction is very
  # sparse, lambda is tiny, and the few invasions that did occur imply a huge beta0 (2.14 at
  # the first fold against 0.20 at the last). Those rows dominate the fit. Measured pooled
  # observed/expected under rolling predictors: 0.440 with no floor -- WORSE than changing
  # nothing -- against 0.741 at t >= 3 on the realised brms run. Never run this without a
  # floor. (An earlier note here said 1.071; that came from an ML harness where four early
  # folds separated, gave beta0 = 0 and predicted zero events. brms regularises them.)
  .roll <- !is.null(rolling_ll) && !is.null(rolling_delay) &&
           exists("reaggregate_asof", mode = "function")
  .wk <- suppressWarnings(as.Date(colnames(Y_wide)))
  if (.roll && anyNA(.wk)) {
    warning("[design] rolling predictors need dated columns on the count matrix; ",
            "falling back to the shared matrix.", call. = FALSE)
    .roll <- FALSE
  }
  # THE FLOOR MUST NEVER EMPTY A FOLD. It skips sparse early transitions, which costs a
  # long fold almost nothing -- but the EARLIEST fold has only two transitions, so a floor of
  # 3 skipped both and the fold produced no training rows and therefore NO FORECASTS AT ALL.
  # That silently removed the first round from every Bayesian model while the structural
  # baselines, which use no design, kept it: the earliest weeks of the outbreak went
  # uncross-validated, and the common-support alignment then (correctly) withheld that round
  # from the baselines too, so it vanished from the comparison entirely.
  #
  # Capping the floor at the number of available transitions is safe for the thing the floor
  # protects. The import force enters the hazard as `offset(logLam)` with its coefficient
  # FIXED at 1 (21_bayesian_renewal.R), so the only term a sparse early fold destabilises is
  # the INTERCEPT -- and an intercept shifts every zone identically on the cloglog scale.
  # Within-fold ranking (top-K, AUC-PR, rank-of-truth, the prioritisation panels) is therefore
  # untouched by an unstable beta_0; only the probability LEVEL moves, which is exactly what
  # the one-parameter recalibration of 16b exists to correct. One transition still supplies a
  # row per at-risk zone (~500), so the design is not thin in rows, only in import force.
  .nT <- ncol(Y_wide) - 1L
  .t0 <- if (.roll) min(max(1L, as.integer(rolling_floor)), max(1L, .nT)) else 1L
  if (.roll && .t0 < max(1L, as.integer(rolling_floor)))
    message(sprintf(paste0("[design] rolling floor %d capped to %d: this fold has only %d ",
                           "transition(s). The fold is KEPT (an unstable intercept cannot ",
                           "move the within-fold ranking); its probability level is carried ",
                           "by the recalibration factor."),
                    as.integer(rolling_floor), .t0, .nT))
  rows <- list()
  for (t in seq_len(ncol(Y_wide) - 1L)) {
    if (.roll && t < .t0) next
    # The import force for transition t -> t+1. Under rolling predictors it comes from the
    # counts known at the END of week t (week_start + 6, the same origin convention used
    # everywhere else); otherwise from the single shared matrix.
    Y_lam <- Y_wide
    if (.roll) {
      Y_lam <- .asof_counts_wide(rolling_ll, zones_all, .wk[t] + 6L, rolling_delay,
                                 week_spine = .wk[seq_len(t)])
      if (is.null(Y_lam) || ncol(Y_lam) < 1L) next
    }
    # compute_foi() reads weeks strictly BEFORE t_idx, indexed against the matrix it is given.
    # Shared matrix: the target week sits at column t + 1.
    # Rolling matrix: it holds weeks 1..t only, so the target week is the column after the end.
    .t_idx <- if (.roll) ncol(Y_lam) + 1L else t + 1L
    Lam <- compute_foi(Y_lam, W, G, t_idx = .t_idx, zones_all)
    aff <- rowSums(Y_wide[, seq_len(t), drop = FALSE] > 0) > 0
    ar  <- !aff & Lam > 0
    if (!any(ar)) next
    X <- .feature_matrix(t + 1L, Y_wide, A_wide, W, G, static, osrm_mat, zones_all,
                         candidates, S_wide = S_wide)
    df <- data.frame(invaded = as.integer(Y_wide[ar, t + 1L] > 0), logLam = log(Lam[ar]),
                     .zone = zones_all[ar], .week = t + 1L, stringsAsFactors = FALSE)
    if (!is.null(X)) df <- cbind(df, as.data.frame(X[ar, , drop = FALSE]))
    rows[[length(rows) + 1L]] <- df
  }
  d <- dplyr::bind_rows(rows)
  feat <- intersect(candidates, names(d))
  center <- stats::setNames(numeric(length(feat)), feat)
  scale  <- stats::setNames(numeric(length(feat)), feat)
  for (f in feat) { m <- mean(d[[f]], na.rm = TRUE); s <- stats::sd(d[[f]], na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- 1; center[[f]] <- m; scale[[f]] <- s
    d[[f]] <- (d[[f]] - m) / s }
  beta0 <- tryCatch(exp(unname(coef(suppressWarnings(
    stats::glm(invaded ~ 1, binomial("cloglog"), offset = d$logLam, data = d)))[1])),
    error = function(e) NA_real_)
  list(d = d, feat = feat, n_events = sum(d$invaded), n_obs = nrow(d), beta0 = beta0,
       center = center, scale = scale, mob = mob, gt = gt)
}

# covariate_associations() and plot_model_parameters() — REMOVED. They ran the Firth-penalised
# cloglog covariate-association SCREEN (marginal hazard ratios with 95% Wald CIs, separation
# flags) on build_invasion_design()'s output and drew it. Both were reachable only from the
# frequentist branch of run_all.R, which is gone; the surviving covariate evidence is the
# BAYESIAN posterior, traced over folds by compute_bayes_params_over_time() and published as
# key_outputs/bayes_params_over_time.csv. build_invasion_design() itself is retained — the
# Bayesian suite builds every one of its designs with it.

#' National effective reproduction number R(t) from the renewal estimate
#' (EpiNow2), with 60% and 90% credible bands and the three GT-profile means for
#' sensitivity. R(t) is a renewal-model output (it is what projects incidence one
#' week ahead for the 2-week invasion horizon); this plots it explicitly.
plot_rt <- function(rt_primary, rt_all = NULL, save = TRUE) {
  if (is.null(rt_primary) || nrow(rt_primary) == 0) return(invisible(NULL))
  # ESTIMATES ONLY. EpiNow2 projects 7 days past the data (forecast_opts(horizon = 7)), and
  # those rows are tagged type == "forecast". Plotting them in the same ink as the estimates
  # showed a week of projection as if it were inference — and it is the most eye-catching part
  # of the curve, being the part that moves. Drop it here and in the GT-sensitivity overlay.
  .tail_removed <- TRUE
  .est_only <- function(x) {
    if (is.null(x) || !nrow(x)) return(x)
    # !%in% "forecast", NOT != : base-R logical indexing with an NA returns an all-NA PHANTOM
    # ROW rather than dropping it.
    if ("type" %in% names(x)) return(x[!x$type %in% "forecast", , drop = FALSE])
    # No `type` column: every cached R(t) object written before RT_CACHE_VERSION 6 lacks it.
    # Fall back to the as-of date, as run_all.R and 44_reff_epinow2_check.R both do — silently
    # returning the input would leave EpiNow2's 7-day forecast tail drawn in the same ink as the
    # estimates while the caption claims it was excluded.
    .cut <- suppressWarnings(as.Date(get0("ANALYSIS_DATE", ifnotfound = NA)))
    if (length(.cut) == 1L && !is.na(.cut))
      return(x[as.Date(x$date) <= .cut, , drop = FALSE])
    .tail_removed <<- FALSE
    warning("[rt] no `type` column and no ANALYSIS_DATE: the EpiNow2 forecast tail cannot be excluded.",
            call. = FALSE)
    x
  }
  rt_primary <- .est_only(rt_primary)
  if (!nrow(rt_primary)) return(invisible(NULL))
  d <- rt_primary %>% dplyr::mutate(date = as.Date(date))
  p <- ggplot(d, aes(date, R_mean)) +
    geom_hline(yintercept = 1, linetype = "22", colour = "grey40") +
    geom_ribbon(aes(ymin = R_lo_90, ymax = R_hi_90), fill = OKABE_ITO[1], alpha = 0.15) +
    geom_ribbon(aes(ymin = R_lo_60, ymax = R_hi_60), fill = OKABE_ITO[1], alpha = 0.30)
  if (!is.null(rt_all) && length(rt_all) > 1) {
    others <- dplyr::bind_rows(lapply(names(rt_all), function(nm) {
      x <- .est_only(rt_all[[nm]]); if (is.null(x) || !nrow(x)) return(NULL)
      dplyr::mutate(x, date = as.Date(date), gt = nm) }))
    p <- p + geom_line(data = others, aes(date, R_mean, linetype = gt),
                       colour = "grey30", linewidth = 0.4, inherit.aes = FALSE) +
      scale_linetype_manual(values = c(short = "dotted", medium = "solid", long = "dashed"),
                            name = "GT profile")
  }
  p <- p + geom_line(colour = OKABE_ITO[1], linewidth = 0.9) +
    labs(title = "National effective reproduction number R(t)",
         subtitle = "Renewal (EpiNow2) estimate | shaded = 60% & 90% credible intervals | R=1 = growth/decline threshold",
         x = NULL, y = "R(t)",
         caption = paste0("A renewal-model output on the national onset series; also used to project one week of incidence for the 2-week invasion horizon. ",
                          if (.tail_removed) "Estimates only: EpiNow2's 7-day forecast tail is excluded."
                          else "WARNING: the 7-day forecast tail could NOT be identified and is included.")) +
    theme_inv(11) + theme(plot.caption = element_text(colour = "grey40", size = 8, hjust = 0))
  if (save) .fd_save(p, file.path(OUT_DIAGNOSTICS, "rt_national.pdf"), w = 9, h = 4.6)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Q2 + Q3 — 0-1 relative risk, and a vulnerability/capacity-adjusted priority
# ---------------------------------------------------------------------------

# Percentile rank in [0,1] (0 = lowest, 1 = highest), NA-safe. `ties = "average"` by default;
# use "min" for an axis with a meaningful FLOOR that is heavily tied (e.g. travel time to the
# nearest facility, which is 0 for the ~all zones that have their own facility), so that floor
# maps to 0 rather than the tie centroid (~0.5) — otherwise best-in-class maps to "median".
.rank01 <- function(x, ties = "average") {
  x[!is.finite(x)] <- NA_real_   # Inf/-Inf would rank as ordinary extremes and push percentiles past 1
  n <- sum(is.finite(x)); if (n <= 1) return(rep(0.5, length(x)))
  (rank(x, na.last = "keep", ties.method = ties) - 1) / (n - 1)
}

#' Q3 — transparent vulnerability-and-capacity index V in [0,1] (1 = most
#' vulnerable / least prepared). Improves the ebola_v23 composite by (a) using
#' data-driven percentile ranks instead of hardcoded denominators, (b) reducing
#' to the exogenous capacity/vulnerability pillars the user named — surveillance,
#' healthcare load, healthcare access, deprivation — and (c) being wired to the
#' actual invasion hazard downstream (the priority score below), which ebola_v23's
#' shipped code never did.
#'
#' Pillars (each a percentile "gap", higher = worse-off):
#'   * surveillance gap = 1 - rank(health-facility DENSITY)  [detection reach;
#'       dedicated PCR-testing data covers only ~19/519 zones so facility density
#'       is the robust per-zone surveillance proxy]
#'   * healthcare gap   = 1 - rank(health facilities per CAPITA)  [treatment load]
#'   * healthcare-access gap = rank(travel time to nearest zone with a facility;
#'       0 for own-facility zones)  [physical reach of care; needs an OSRM matrix]
#'   * social vulnerability = rank(CCVI socioeconomic deprivation)
#' V = equal-weight mean of the AVAILABLE pillars (four when an OSRM travel-time
#' matrix is supplied, as in the default run; the access axis is dropped and V
#' averages the remaining three otherwise).
compute_vulnerability_index <- function(covariates, zones_all = NULL, osrm_mat = NULL) {
  cov <- covariates
  if (!is.null(zones_all)) cov <- cov[cov$nom %in% zones_all, , drop = FALSE]
  gv <- function(nm) if (nm %in% names(cov)) suppressWarnings(as.numeric(cov[[nm]])) else rep(NA_real_, nrow(cov))
  dens   <- gv("healthsite_density")
  hcount <- gv("healthsite_count"); pop <- pmax(gv("pop_count"), 1)
  hc_percap <- hcount / pop
  ccvi   <- gv("ccvi"); if (all(is.na(ccvi))) ccvi <- gv("socioeconomic_deprivation")
  # Healthcare-ACCESS axis (Task 2): road travel time (minutes) from each zone to
  # the NEAREST zone that has any health facility (0 for zones with their own
  # facilities). Zones physically far from care are more vulnerable — this is a
  # distinct axis from facility DENSITY (surveillance) and facilities-per-capita
  # (healthcare load): a zone can have few local facilities yet sit close to a
  # well-served neighbour, or vice-versa. Derived from the OSRM travel-time matrix.
  t_access <- rep(NA_real_, nrow(cov))
  if (!is.null(osrm_mat) && !is.null(rownames(osrm_mat))) {
    zn <- as.character(cov$nom); have <- zn %in% rownames(osrm_mat)
    fac_zones <- intersect(zn[is.finite(hcount) & hcount > 0], colnames(osrm_mat))
    if (length(fac_zones) && any(have)) {
      sub <- osrm_mat[zn[have], fac_zones, drop = FALSE]
      tt  <- apply(sub, 1, function(r) { r <- r[is.finite(r)]; if (length(r)) min(r) else NA_real_ })
      ta  <- rep(NA_real_, nrow(cov)); ta[have] <- tt
      ta[is.finite(hcount) & hcount > 0] <- 0            # own facilities -> zero access time
      big <- suppressWarnings(max(ta[is.finite(ta)], na.rm = TRUE))
      ta[!is.finite(ta)] <- if (is.finite(big)) big else NA_real_   # unreachable -> worst
      t_access <- ta
    }
  }
  # ties="min" so own-facility zones (t_access = 0, the floor) map to access_gap = 0, not the
  # tie centroid — the access axis is otherwise degenerate here (nearly every zone has a facility).
  access_gap <- if (all(is.na(t_access))) rep(NA_real_, nrow(cov)) else .rank01(t_access, ties = "min")
  d <- tibble::tibble(
    nom            = cov$nom,
    surveillance_gap = 1 - .rank01(dens),
    healthcare_gap   = 1 - .rank01(hc_percap),
    access_gap       = access_gap,
    social_vulnerability = .rank01(ccvi),
    healthcare_travel_min = t_access)
  # Equal-weight mean over the AVAILABLE pillars (access axis is dropped if no OSRM
  # matrix is supplied, so behaviour is unchanged for callers that omit it).
  pillars <- c("surveillance_gap", "healthcare_gap", "access_gap", "social_vulnerability")
  # DROP DEGENERATE PILLARS, not just all-NA ones (2026-09-17).
  #
  # Every one of the 519 zones has at least one health facility (min healthsite_count = 5), so
  # the `hcount > 0 -> t_access = 0` override above zeroes the ENTIRE access axis: t_access is 0
  # for all 519 zones, hence access_gap is 0 for all of them. The old filter kept it (it is not
  # all-NA), so V was the mean of THREE informative pillars and one constant zero — i.e. exactly
  # 3/4 of the intended index. Verified against outputs/forecasts/bayes_risk_scores_current.rds:
  # max|V - 0.75 * mean(3 informative pillars)| = 5.6e-17, published V range [0.014, 0.734]
  # against a correct [0.019, 0.978].
  #
  # This deflated the published V column and compressed the x-axis of the retained
  # bayes_priority_scatter panels (a percent-formatted axis topping out at 73%). `priority` and
  # all ranks were UNAFFECTED — the constant 3/4 cancels in the priority rescale and in rank(V).
  #
  # A zero-variance axis carries no information about which zone is worse off, so averaging it in
  # can only dilute. Keep a pillar only if it is neither all-NA nor constant.
  .informative <- function(p) {
    v <- d[[p]]
    !all(is.na(v)) && is.finite(stats::sd(v, na.rm = TRUE)) && stats::sd(v, na.rm = TRUE) > 0
  }
  .dropped <- pillars[!vapply(pillars, .informative, logical(1))]
  pillars  <- pillars[vapply(pillars, .informative, logical(1))]
  if (length(.dropped))
    message(sprintf("[vulnerability] %d pillar(s) dropped as uninformative (all-NA or zero variance): %s; V is the mean of the remaining %d.",
                    length(.dropped), paste(.dropped, collapse = ", "), length(pillars)))
  if (!length(pillars))
    stop("[vulnerability] no informative pillar remains; V cannot be formed.", call. = FALSE)
  # Publish the pillar set that ACTUALLY formed V, so captions and the methods description
  # derive it instead of asserting a fixed count. The published caption used to name four
  # pillars including "healthcare access (travel time to nearest facility)" while the access
  # axis is constant on this data and is dropped just above — so the retained
  # bayes_priority_scatter panels described an index component that does not enter the index.
  .V_LABELS <- c(surveillance_gap     = "surveillance (facility density)",
                 healthcare_gap       = "healthcare (facilities per capita)",
                 access_gap           = "healthcare access (travel time to nearest facility)",
                 social_vulnerability = "CCVI deprivation")
  assign(".V_PILLARS", pillars, envir = .GlobalEnv)
  assign(".V_PILLAR_LABELS", unname(ifelse(pillars %in% names(.V_LABELS),
                                           .V_LABELS[pillars], pillars)),
         envir = .GlobalEnv)
  M <- as.matrix(d[, pillars, drop = FALSE])
  # Equal-weight mean over the available pillars; a zone with EVERY pillar NA gets V = NA (not
  # NaN, which rowMeans(na.rm=TRUE) would otherwise return and then propagate into priority).
  d$V <- ifelse(rowSums(is.finite(M)) == 0L, NA_real_, rowMeans(M, na.rm = TRUE))
  d
}

#' One sentence naming the pillars that actually formed V on THIS run, from the set
#' .attach_vulnerability() published. Falls back to the full nominal list only when the
#' index has not been built yet (standalone use), and never asserts a pillar count.
.v_pillar_sentence <- function() {
  labs <- get0(".V_PILLAR_LABELS", ifnotfound = NULL)
  if (is.null(labs) || !length(labs))
    labs <- c("surveillance (facility density)", "healthcare (facilities per capita)",
              "healthcare access (travel time to nearest facility)", "CCVI deprivation")
  joined <- if (length(labs) == 1L) labs else
    paste0(paste(labs[-length(labs)], collapse = ", "), " and ", labs[length(labs)])
  sprintf("Vulnerability = mean percentile gap in %s.", joined)
}

#' Q2 + Q3 — attach 0-1 relative-risk indices and the vulnerability-adjusted
#' preparedness-priority score to a risk table.
#'   rr01_nat / rr01_ituri : p_invasion / max(p_invasion) within the geography,
#'       so 1 = the single highest-risk at-risk zone (0-1, Q2).
#'   priority : rr01_nat * V, rescaled to [0,1] (1 = top preparedness priority):
#'       high where a zone is BOTH likely to be invaded AND vulnerable / poorly
#'       resourced (Q3). Multiplicative so a zone needs both to rank high.
add_risk_indices <- function(risk_scores, vuln,
                             provinces = get0("PROVINCES_OF_INTEREST",
                                              ifnotfound = c("Ituri", "Nord-Kivu", "Haut-Uele"))) {
  rs <- risk_scores %>% dplyr::left_join(
    vuln %>% dplyr::rename(health_zone = nom), by = "health_zone")
  if (!"province" %in% names(rs)) rs$province <- NA_character_
  rs$.active <- rs$was_active_before %in% TRUE
  rs <- rs %>% dplyr::group_by(horizon) %>%
    dplyr::mutate(
      .mx_nat = suppressWarnings(max(p_invasion[!.active], na.rm = TRUE)),
      rr01_nat = ifelse(.active, NA_real_, p_invasion / .mx_nat),
      priority_raw = ifelse(.active, NA_real_, rr01_nat * V),
      .mx_pr = suppressWarnings(max(priority_raw, na.rm = TRUE)),
      priority = ifelse(is.finite(.mx_pr) & .mx_pr > 0, priority_raw / .mx_pr, NA_real_),
      priority_rank = dplyr::min_rank(dplyr::desc(priority))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.mx_nat, -.mx_pr)
  # 0-1 relative invasion risk WITHIN each province of interest (drives the
  # province-specific choropleths); rr01_ituri retained backward-compatibly.
  psfx <- get0(".prov_suffix", ifnotfound = function(p) tolower(gsub("[^A-Za-z]", "", p)))
  for (prov in provinces) {
    s <- psfx(prov)
    rs <- rs %>% dplyr::group_by(horizon) %>%
      dplyr::mutate(
        .inp = province %in% prov,
        .mx  = suppressWarnings(max(p_invasion[!.active & .inp], na.rm = TRUE)),
        !!paste0("rr01_", s) := ifelse(!.active & .inp & is.finite(.mx),
                                       p_invasion / .mx, NA_real_)) %>%
      dplyr::ungroup() %>% dplyr::select(-.inp, -.mx)
  }
  rs %>% dplyr::select(-.active, -dplyr::any_of("priority_raw"))
}

#' Q3 — hazard-vs-vulnerability scatter: where a zone sits on likelihood (y) and
#' vulnerability/capacity (x); the diagonal is the priority score. Top-priority
#' zones (upper right) are labelled.
#' @param show_title FALSE drops the plot title AND subtitle (the caption stays). The
#'   published bayes_priority_scatter panels carry their caption in the manuscript text, so
#'   in-figure titles duplicate it; other callers keep titles for at-a-glance diagnostics.
plot_priority_scatter <- function(rs, horizon = 1L, top_n = 12L, save = TRUE, window_txt = "",
                                  file = "priority_scatter", model_label = "",
                                  show_title = TRUE) {
  d <- rs %>% dplyr::filter(horizon == !!horizon, !(was_active_before %in% TRUE),
                            is.finite(V), is.finite(rr01_nat))
  if (nrow(d) == 0) return(invisible(NULL))
  if (!"province" %in% names(d)) d$province <- NA_character_
  d <- d %>% dplyr::mutate(region = .prov_label(province))
  lab <- d %>% dplyr::arrange(dplyr::desc(priority)) %>% head(top_n)
  p <- ggplot(d, aes(V, rr01_nat)) +
    geom_point(aes(size = priority, colour = region), alpha = 0.75) +
    geom_text(data = lab, aes(label = health_zone), size = 3, vjust = -0.9,
              check_overlap = TRUE, colour = "grey15") +
    scale_colour_manual(values = .prov_colours(d$region), name = "Province") +
    scale_size_area(max_size = 8, name = "Priority") +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = if (show_title)
           sprintf("Preparedness priority: invasion likelihood x vulnerability (h=%dw)%s",
                   horizon, if (nzchar(model_label)) paste0(" — ", model_label) else "")
         else NULL,
         subtitle = if (show_title)
           "Upper-right = likely to be invaded AND vulnerable/under-resourced = highest priority"
         else NULL,
         x = "Vulnerability & capacity gap (0-1; higher = more vulnerable, less prepared)",
         y = "Relative invasion risk (0-1; 1 = highest)",
         caption = paste0(.v_pillar_sentence(), " ", window_txt)) +
    theme_inv(11) + theme(legend.position = "right",
                          plot.caption = element_text(colour = "grey40", size = 8))
  if (save) {
    .pp <- file.path(OUT_REPORTS, sprintf("%s_h%d.pdf", file, horizon))
    .fd_save(p, .pp, w = 9.5, h = 6.5)
    # Also place the priority scatter in key_outputs/ (in ADDITION to reports) so it ships
    # with the headline deliverables — copied here at save time, not only via the end-of-run
    # key_outputs manifest, so it is present even on a partial run.
    # ASK THE GATE before copying. .fd_save() returns early when figure_is_kept() refuses,
    # WITHOUT clearing a stale target, so an unconditional copy republished whatever an earlier
    # run had left in reports/ into the published key-outputs bundle — every run, with a fresh
    # mtime, and outside the reach of the end-of-run archive sweep (which does not cover
    # key_outputs/ itself or OUT_REPORTS). Copy only a file this gate would write, and clear
    # any stale copy otherwise.
    .ko  <- file.path(get0("OUT_DIR", ifnotfound = dirname(OUT_REPORTS)), "key_outputs")
    .fk  <- get0("figure_is_kept", ifnotfound = NULL)
    .keep <- !is.function(.fk) || isTRUE(.fk(.pp))
    .dst <- file.path(.ko, basename(.pp))
    if (.keep && file.exists(.pp)) {
      dir.create(.ko, showWarnings = FALSE, recursive = TRUE)
      file.copy(.pp, .dst, overwrite = TRUE)
    } else if (file.exists(.dst)) {
      file.remove(.dst)
    }
  }
  invisible(p)
}

#' Q3 — top-priority zones with the component breakdown (invasion hazard and the
#' three vulnerability pillars), so the score is fully transparent.
plot_priority_bars <- function(rs, horizon = 1L, top_n = 15L, save = TRUE,
                               file = "priority_bars", model_label = "") {
  d <- rs %>% dplyr::filter(horizon == !!horizon, !(was_active_before %in% TRUE),
                            is.finite(priority)) %>%
    dplyr::arrange(dplyr::desc(priority)) %>% head(top_n)
  if (nrow(d) == 0) return(invisible(NULL))
  # Show every vulnerability pillar that feeds V (access_gap included when present).
  comp <- c(`Invasion risk (0-1)` = "rr01_nat", `Surveillance gap` = "surveillance_gap",
            `Healthcare gap` = "healthcare_gap", `Access gap` = "access_gap",
            `Social vulnerability` = "social_vulnerability")
  comp <- comp[vapply(comp, function(cc) cc %in% names(d) && !all(is.na(d[[cc]])), logical(1))]
  long <- d %>% dplyr::select(health_zone, dplyr::all_of(unname(comp))) %>%
    dplyr::mutate(health_zone = factor(health_zone, levels = rev(health_zone))) %>%
    tidyr::pivot_longer(-health_zone) %>%
    dplyr::mutate(name = factor(names(comp)[match(name, comp)], levels = names(comp)))
  pal <- c(`Invasion risk (0-1)` = OKABE_ITO[8], `Surveillance gap` = OKABE_ITO[1],
           `Healthcare gap` = OKABE_ITO[3], `Access gap` = OKABE_ITO[5],
           `Social vulnerability` = OKABE_ITO[2])
  p <- ggplot(long, aes(value, health_zone, fill = name)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    scale_fill_manual(values = pal[names(comp)], name = NULL) +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.03))) +
    labs(title = sprintf("Top preparedness-priority zones — component breakdown (h=%dw)%s",
                         horizon, if (nzchar(model_label)) paste0(" — ", model_label) else ""),
         subtitle = "Priority = invasion risk x mean(vulnerability pillars); each component shown 0-1",
         x = "Component value (0-1)", y = NULL) +
    theme_inv(11) + theme(legend.position = "top")
  if (save) .fd_save(p, file.path(OUT_REPORTS, sprintf("%s_h%d.pdf", file, horizon)), w = 10, h = 6.5)
  invisible(p)
}

#' Q3 — priority choropleth (Ituri + national), masked for affected zones.
plot_priority_map <- function(rs, shapefile = NULL, horizon = 1L, save = TRUE, window_txt = "",
                              file = "priority_map", model_label = "") {
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile)) return(invisible(NULL))
  d <- rs %>% dplyr::filter(horizon == !!horizon)
  aff_keys <- d %>% dplyr::filter(was_active_before %in% TRUE) %>%
    dplyr::pull(health_zone) %>% trimws() %>% tolower() %>% unique()
  p_it  <- .fd_map(d, shapefile, "Ituri — preparedness priority", prob_col = "priority",
                   ituri_only = TRUE, lim = c(0, 1), palette = "inferno",
                   legend_name = "Priority (0-1)", affected_keys = aff_keys)
  p_nat <- .fd_map(d, shapefile, "National", prob_col = "priority",
                   ituri_only = FALSE, lim = c(0, 1), palette = "inferno",
                   legend_name = "Priority", affected_keys = aff_keys)
  out <- p_it + p_nat + patchwork::plot_layout(widths = c(2, 1)) +
    patchwork::plot_annotation(
      title = sprintf("Vulnerability-adjusted preparedness priority (h=%dw)%s", horizon,
                      if (nzchar(model_label)) paste0(" — ", model_label) else ""),
      subtitle = "Invasion risk x vulnerability/capacity gap | grey = already affected",
      caption = window_txt,
      theme = ggplot2::theme(plot.title = element_text(face = "bold", size = 13),
                             plot.caption = element_text(colour = "grey40", size = 8, hjust = 0)))
  if (save) .fd_save(out, file.path(OUT_MAPS, sprintf("%s_h%d.pdf", file, horizon)), w = 12, h = 6.5)
  invisible(out)
}

#' Q3 — markdown describing the priority index (for the report).
describe_priority_index <- function() {
  c("## Vulnerability-adjusted preparedness priority", "",
    "Alongside the invasion probability we report a **preparedness-priority score**",
    "that asks not just *how likely* a zone is to be invaded but *how badly it would",
    "go*. It multiplies the invasion hazard by a transparent **vulnerability &",
    "capacity index** V in [0,1] (1 = most vulnerable / least prepared), the equal-",
    # The pillar list is DERIVED from the run, not asserted. This used to promise "four when an
    # OSRM travel-time matrix is supplied, as in the default run" — which is precisely the case
    # where the access axis IS dropped: every one of the 519 zones has a facility, so t_access
    # is 0 everywhere and the zero-variance pillar is excluded by .attach_vulnerability().
    # V is therefore a THREE-pillar mean on the shipped data, and the numbered list below marks
    # any nominal pillar that did not enter it.
    sprintf("weight mean of the percentile-ranked pillars that carry information on this run (%d of 4 below):",
            length(get0(".V_PILLARS", ifnotfound = c("surveillance_gap", "healthcare_gap",
                                                     "access_gap", "social_vulnerability")))), "",
    local({
      used <- get0(".V_PILLARS", ifnotfound = c("surveillance_gap", "healthcare_gap",
                                                "access_gap", "social_vulnerability"))
      mark <- function(k) if (k %in% used) "" else " *(dropped this run: no variation across zones)*"
      c(paste0("1. **Surveillance gap** = 1 - rank(health-facility density) — detection reach",
               mark("surveillance_gap")),
        "   (dedicated PCR-testing data covers only ~4% of zones, so facility density is",
        "   the robust per-zone surveillance proxy);",
        paste0("2. **Healthcare gap** = 1 - rank(health facilities per capita) — treatment load;",
               mark("healthcare_gap")),
        "3. **Healthcare access gap** = rank(road travel time to the nearest zone with a facility;",
        paste0("   own-facility zones = 0) — physical reach of care, distinct from facility density/per-capita;",
               mark("access_gap")),
        paste0("4. **Social vulnerability** = rank(CCVI socioeconomic deprivation).",
               mark("social_vulnerability")))
    }), "",
    "`priority = (invasion risk, 0-1) x V`, rescaled so 1 is the top priority. Being",
    "multiplicative, a zone must be BOTH at material risk of invasion AND vulnerable",
    "to score highly — a well-resourced zone that is likely to be invaded, or a",
    "vulnerable zone at negligible risk, both rank lower. This improves the ebola_v23",
    "composite by using data-driven percentile ranks (not hardcoded denominators), the",
    "named exogenous pillars, and by actually wiring the score to the fitted invasion",
    "model — which the shipped ebola_v23 prioritisation never did.", "")
}

# ---------------------------------------------------------------------------
# Date-window annotation (Task 2): every forecast figure states its FIT window
# (training data used) and its PREDICTION window (weeks being forecast).
# ---------------------------------------------------------------------------
#' Caption stating the fit and prediction date windows for a cutoff + horizons.
.window_caption <- function(train_start, cutoff, horizons) {
  cutoff <- as.Date(cutoff); train_start <- as.Date(train_start)
  hmax <- max(horizons)
  # cutoff is the bucket START (needed for the prediction math); the last day actually
  # trained on is the bucket END = cutoff + 6 (== ANALYSIS_DATE under the analysis-date-
  # anchored grid), and the last predicted week ends at cutoff + hmax*7 + 6. Both window
  # endpoints are shown as inclusive first/last days.
  sprintf("Fit window: %s to %s (training). Prediction window: %s to %s (next %d week%s).",
          format(train_start, "%d %b %Y"), format(cutoff + 6, "%d %b %Y"),
          format(cutoff + 7, "%d %b %Y"), format(cutoff + hmax * 7 + 6, "%d %b %Y"),
          hmax, if (hmax > 1) "s" else "")
}

# ---------------------------------------------------------------------------
# Task 6 — Bayesian suite visualisations (posterior parameters + posterior
# invasion probabilities WITH credible intervals + stacking weights)
# ---------------------------------------------------------------------------

#' Spell out a Bayesian model code (e.g. "Bayes-M8-geo") as its modelling
#' assumptions, so plots are legible without a codebook.
.bayes_model_label <- function(code) {
  vapply(code, function(cc) {
    mob  <- sub("^Bayes-(M[0-9]+).*$", "\\1", cc)   # matches M4 in both "M4" and "M4-dist"
    mobn <- c(M1 = "short-trip", M4 = "gravity", M8 = "composite-gravity",
              M9 = "multi-kernel ensemble", M10 = "radiation-composite",
              M11 = "inward meeting-location FOI",
              M13 = "cohort+gravity", M14 = "cohort+radiation",
              # M15/M16/M17 were missing, so the Flowminder-static family rendered as bare
              # matrix ids in every legend that uses this decoder.
              # M16 composes the cohort rows with the DIRECTED OD kernel M3, not with the
              # symmetrised static M15 (03_mobility_matrices.R); "cohort+static" named the
              # wrong base kernel in every legend that uses this decoder.
              M15 = "combined-static", M16 = "cohort+relocation OD",
              M17 = "all-kernel consensus")[mob]
    # OSRM road-DISTANCE deterrence variant (vs the default travel-time deterrence).
    if (grepl("-dist", cc) && !is.na(mobn)) mobn <- paste0(mobn, ", road-distance")
    # SOURCE-CELL FILL twin: without this it decodes to exactly its unfilled parent's string.
    if (grepl("-fill", cc) && !is.na(mobn)) mobn <- paste0(mobn, ", source-cell fill")
    # ORIGIN-SPLIT cohort kernel (always filled by construction).
    if (grepl("-split", cc) && !is.na(mobn)) mobn <- paste0(mobn, ", origin-split cohort, source-cell fill")
    gt   <- if (grepl("-short", cc)) "short GT" else if (grepl("-long", cc)) "long GT" else "med GT"
    cov  <- if (grepl("-full", cc)) " · +full covariates" else
            if (grepl("-geo", cc)) " · +geo covariates" else ""
    lnk  <- if (grepl("-logit", cc)) " · logit link" else if (grepl("-probit", cc)) " · probit link" else ""
    # TIME-VARYING beta family. Without this a "-tv*" model decodes to exactly the same
    # string as its fixed-beta base kernel, so two genuinely different models would be
    # indistinguishable in every legend and in the parameter-panel key. Every process in
    # the grid must appear here: until 2026-09-23 only trend and week did, so an ar1, rw1
    # or gp arm silently carried its base kernel's label.
    .tvlab <- c(trend = "time-trend", week = "weekly-varying", rw1 = "random-walk",
                ar1 = "AR(1)", gp = "Gaussian-process")
    .tvhit <- regmatches(cc, regexpr("-tv[a-z0-9]+$", cc))
    tvv  <- if (length(.tvhit)) {
      .k <- sub("^-tv", "", .tvhit)
      sprintf(" · %s β", if (.k %in% names(.tvlab)) .tvlab[[.k]] else .k)
    } else ""
    sprintf("%s  —  %s · %s%s%s%s", cc, if (is.na(mobn)) mob else mobn, gt, cov, lnk, tvv)
  }, character(1))
}

#' Bayesian invasion-model PARAMETERS (posterior) for the BEST-FITTING subset of the
#' suite. The full suite has ~dozens of mobility/GT/covariate variants whose spelled-out
#' names (`.bayes_model_label()`) are far too long to serve as y-axis ticks — the panels
#' collapsed to invisibility. So we rank models by loo predictive-stacking weight (`weights`;
#' the same metric that drives the ensemble), keep the top `top_n`, give each a SHORT code
#' (B1, B2, …) on the axis, and spell out what every code means in a custom legend panel.
#' Two panels: (a) the import coefficient beta0 = exp(intercept) for each shown model; and
#' (b) covariate hazard ratios for the shown models that include covariates.
plot_bayes_parameters <- function(params, weights = NULL, window_txt = "", save = TRUE,
                                  top_n = 6L) {
  if (is.null(params) || !nrow(params)) return(invisible(NULL))
  lab <- function(t) .lk(.PARAM_LABEL, t, t)
  # description WITHOUT the leading "Bayes-… — " code prefix, for the legend rows.
  descr_of <- function(code) sub("^.*?  —  ", "", .bayes_model_label(code))

  # ---- pick the best-fitting subset, rank by loo stacking weight ------------------
  models <- unique(params$model)
  if (!is.null(weights) && length(weights)) {
    wv  <- weights[intersect(names(sort(weights, decreasing = TRUE)), models)]
    sel <- head(names(wv), top_n)
  } else {
    # No weights available (e.g. equal-weight fallback): can't rank, so keep the first
    # top_n and flag it in the legend rather than crushing every model onto the axis.
    sel <- head(models, top_n); wv <- stats::setNames(rep(NA_real_, length(sel)), sel)
  }
  keymap <- stats::setNames(sprintf("B%d", seq_along(sel)), sel)   # B1 = best-fitting
  # factor levels reversed so B1 sits at the TOP of each panel.
  klev   <- rev(unname(keymap[sel]))
  key_f  <- function(m) factor(unname(keymap[m]), levels = klev)

  # (a) import coefficient beta0 for each shown model
  base <- params %>% dplyr::filter(is_intercept, model %in% sel) %>%
    dplyr::mutate(key = key_f(model))
  pA <- ggplot(base, aes(hr, key)) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.9, colour = OKABE_ITO[3]) +
    geom_point(size = 2.8, colour = OKABE_ITO[3]) +
    scale_x_log10() +
    labs(title = "(a) Import coefficient beta0",
         subtitle = "exp(intercept): calibrated import->first-case hazard (posterior median, 90% CrI)",
         x = "beta0  (log scale)", y = NULL) +
    theme_inv(10)
  # (b) covariate hazard ratios (shown models that include covariates)
  cov <- params %>% dplyr::filter(!is_intercept, model %in% sel) %>%
    dplyr::mutate(label = vapply(term, lab, character(1)),
                  key = key_f(model),
                  sig = dplyr::case_when(lo > 1 ~ "raises risk (90% CrI > 1)",
                                         hi < 1 ~ "lowers risk (90% CrI < 1)",
                                         TRUE ~ "credible interval spans 1"))

  # ---- custom legend: short code -> what the model actually is --------------------
  leg <- tibble::tibble(code = unname(keymap[sel]), descr = descr_of(sel),
                        w = as.numeric(wv[sel]), y = rev(seq_along(sel)))
  leg$rhs <- ifelse(is.na(leg$w), leg$descr,
                    sprintf("%s   (stacking weight %.2f)", leg$descr, leg$w))
  pL <- ggplot(leg, aes(y = y)) +
    geom_text(aes(x = 0, label = code), fontface = "bold", hjust = 0, size = 3.2,
              colour = OKABE_ITO[3]) +
    geom_text(aes(x = 0.05, label = rhs), hjust = 0, size = 3.0) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(expand = ggplot2::expansion(add = 0.6)) +
    labs(title = paste0("Model key — ", length(sel), " best-fitting models",
                        if (all(is.na(leg$w))) "" else " (by loo stacking weight)")) +
    theme_void(10) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0),
          plot.margin = ggplot2::margin(4, 8, 4, 8))

  out <- if (nrow(cov)) {
    # Each (model x covariate) gets its OWN row, grouped into a facet per covariate — instead
    # of dodging several models onto one covariate row (which crushed every interval to a
    # sliver). space="free_y" sizes each facet to its own model count so no row is squashed.
    cov <- cov %>% dplyr::mutate(label = factor(label, levels = unique(label[order(term)])))
    pB <- ggplot(cov, aes(hr, key, colour = sig)) +
      geom_vline(xintercept = 1, linetype = "22", colour = "grey45") +
      geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.9) +
      geom_point(size = 2.8) +
      ggplot2::facet_grid(label ~ ., scales = "free_y", space = "free_y", switch = "y") +
      scale_colour_manual(values = c(`raises risk (90% CrI > 1)` = OKABE_ITO[2],
        `lowers risk (90% CrI < 1)` = OKABE_ITO[1],
        `credible interval spans 1` = "grey55"), name = NULL) +
      scale_x_log10() +
      labs(title = "(b) Covariate hazard ratios — models that include covariates",
           subtitle = "Per +1 SD (standardised); >1 raises invasion risk (posterior median, 90% CrI). One row per model, grouped by covariate.",
           x = "Hazard ratio per +1 SD  (log scale)", y = NULL) +
      theme_inv(10) +
      theme(legend.position = "right", strip.placement = "outside",
            strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 8.5),
            panel.spacing.y = grid::unit(3, "pt"))
    # Heights proportional to the ACTUAL row counts (a: one per model; b: one per
    # model-covariate; legend: one per model), so neither panel is compressed.
    pA / pB / pL + patchwork::plot_layout(
      heights = c(nrow(base) + 1, nrow(cov) + 1, length(sel) + 1))
  } else {
    pA / pL + patchwork::plot_layout(heights = c(nrow(base) + 1, length(sel) + 1))
  }
  out <- out + patchwork::plot_annotation(
    title = "Bayesian invasion-model parameters (posterior) — best-fitting models",
    caption = paste0("Showing the ", length(sel), " best-fitting models (of ", length(models),
                     " in the suite); each is labelled B1… on the axes and spelled out in the model key. ",
                     "Weakly-informative priors (intercept ~ Normal(-3,2); coefficients ~ Normal(0,1) on standardised covariates) regularise the separation that diverged the frequentist joint fit. ",
                     window_txt),
    theme = ggplot2::theme(plot.title = element_text(face = "bold", size = 13),
                           plot.caption = element_text(colour = "grey40", size = 8, hjust = 0)))
  # Height scales with the row counts across the panels plus the legend and titles.
  .h <- if (nrow(cov)) max(8, (nrow(base) + nrow(cov) + length(sel)) * 0.42 + 3)
        else max(4.6, (nrow(base) + length(sel)) * 0.5 + 1.5)
  if (save) .fd_save(out, file.path(OUT_REPORTS, "bayes_parameters.pdf"), w = 11, h = .h)
  invisible(out)
}

#' Full POSTERIOR DISTRIBUTIONS (not just point + interval) of the Bayesian parameters, from
#' posterior draws: (a) the import coefficient beta0 per model as density ridgelines, and
#' (b) the covariate hazard-ratio posteriors, grouped by covariate. Complements
#' plot_bayes_parameters (which shows only medians + 90% CrI).
plot_bayes_posterior_densities <- function(draws, window_txt = "", save = TRUE) {
  if (is.null(draws) || !nrow(draws) || !requireNamespace("ggridges", quietly = TRUE))
    return(invisible(NULL))
  lab <- function(t) .lk(.PARAM_LABEL, t, t)
  b0 <- draws %>% dplyr::filter(is_intercept) %>% dplyr::mutate(mlab = .bayes_model_label(model))
  b0$mlab <- stats::reorder(b0$mlab, b0$hr, FUN = stats::median)
  pA <- ggplot(b0, aes(x = hr, y = mlab)) +
    ggridges::geom_density_ridges(scale = 2.1, rel_min_height = 0.01, fill = OKABE_ITO[3],
      alpha = 0.72, colour = "grey30", linewidth = 0.25, quantile_lines = TRUE,
      quantiles = c(0.05, 0.5, 0.95), vline_colour = "grey25") +
    scale_x_log10() +
    labs(title = "(a) Posterior of the import coefficient beta0, per model",
         subtitle = "Full posterior density of exp(intercept) — calibrated import->first-case hazard (5/50/95% lines)",
         x = "beta0  (log scale)", y = NULL) +
    theme_inv(10)
  cov <- draws %>% dplyr::filter(!is_intercept)
  out <- if (nrow(cov)) {
    cov <- cov %>% dplyr::mutate(label = vapply(term, lab, character(1)), mlab = .bayes_model_label(model))
    cov$label <- factor(cov$label, levels = unique(cov$label[order(cov$term)]))
    pB <- ggplot(cov, aes(x = hr, y = mlab, fill = label)) +
      geom_vline(xintercept = 1, linetype = "22", colour = "grey45") +
      ggridges::geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.7,
                                    colour = "grey30", linewidth = 0.2) +
      ggplot2::facet_grid(label ~ ., scales = "free_y", space = "free_y", switch = "y") +
      scale_x_log10() + scale_fill_viridis_d(option = "viridis", guide = "none") +
      labs(title = "(b) Posterior of covariate effects per +1 SD",
           subtitle = "Mass right of 1 (dashed) raises invasion risk; one density per model, grouped by covariate (hazard ratio for cloglog models; odds ratio for any logit-link model)",
           x = "Effect per +1 SD: hazard/odds ratio  (log scale)", y = NULL) +
      theme_inv(10) +
      theme(strip.placement = "outside",
            strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 8.5),
            panel.spacing.y = grid::unit(3, "pt"))
    pA / pB + patchwork::plot_layout(heights = c(dplyr::n_distinct(b0$model),
                                                 dplyr::n_distinct(paste(cov$model, cov$term))))
  } else pA
  out <- out + patchwork::plot_annotation(
    title = "Bayesian posterior parameter distributions",
    caption = paste0("Full posteriors (not just point + interval). beta0 = exp(intercept); covariate effects are hazard ratios per standardised SD. ", window_txt),
    theme = ggplot2::theme(plot.title = element_text(face = "bold", size = 13),
                           plot.caption = element_text(colour = "grey40", size = 8, hjust = 0)))
  nrw <- dplyr::n_distinct(b0$model) + (if (nrow(cov)) dplyr::n_distinct(paste(cov$model, cov$term)) else 0)
  .h <- max(10, nrw * 0.42 + 4)
  if (save) .fd_save(out, file.path(OUT_REPORTS, "bayes_posterior_densities.pdf"), w = 11, h = .h)
  invisible(out)
}

#' Generation-time PREFERENCE (from bayes_gt_posterior): the loo-predictive pseudo-BMA+ weight
#' per GT mean (bars) against the literature prior (dashed) — which generation time the invasion
#' data actually support. This is a loo-weighted model-average preference, NOT a fully Bayesian
#' posterior (the GT is not jointly estimated), so it is labelled as a preference weight.
plot_bayes_gt_posterior <- function(gt_post, model_label = "", window_txt = "", save = TRUE) {
  if (is.null(gt_post) || nrow(gt_post) < 2) return(invisible(NULL))
  pm  <- attr(gt_post, "pref_mean") %||% sum(gt_post$gt_mean * gt_post$weight)
  pri <- attr(gt_post, "prior")
  prg <- if (!is.null(pri)) { dd <- stats::dnorm(gt_post$gt_mean, pri[["mean"]], pri[["sd"]]); dd / sum(dd) } else NULL
  d <- gt_post; d$prior_w <- if (!is.null(prg)) prg else NA_real_
  p <- ggplot(d, aes(gt_mean, weight)) +
    geom_col(fill = OKABE_ITO[6], alpha = 0.85, width = 0.8) +
    { if (!is.null(prg)) geom_line(aes(y = prior_w), colour = "grey40", linetype = "22", linewidth = 0.6) } +
    geom_vline(xintercept = pm, colour = OKABE_ITO[2], linewidth = 0.8) +
    annotate("text", x = pm, y = max(d$weight, na.rm = TRUE) * 0.96,
             label = sprintf(" preference mean %.1f d", pm), hjust = 0, size = 3.1, colour = OKABE_ITO[2]) +
    labs(title = sprintf("Generation-time preference (loo pseudo-BMA+)%s",
                         if (nzchar(model_label)) paste0(" — ", model_label) else ""),
         subtitle = "loo-predictive pseudo-BMA+ weight per GT mean (bars) vs literature prior (dashed): which GT the invasion data support (not a joint GT posterior)",
         x = "Generation-time mean (days)", y = "loo-preference weight", caption = window_txt) +
    theme_inv(11)
  if (save) .fd_save(p, file.path(OUT_REPORTS, "bayes_gt_posterior.pdf"), w = 8.5, h = 5)
  invisible(p)
}

#' Nowcast-input sensitivity (from bayes_nowcast_sensitivity): the posterior of beta0 when the
#' featured model is refit on differently nowcast training counts (raw / epinowcast / fast
#' delay-CDF). Overlapping posteriors mean the two-stage nowcast choice does not drive the
#' inference; a clear shift would flag the need for a joint nowcast+renewal model.
plot_bayes_nowcast_sensitivity <- function(sens, model_label = "", window_txt = "", save = TRUE) {
  if (is.null(sens) || !nrow(sens) || !requireNamespace("ggridges", quietly = TRUE))
    return(invisible(NULL))
  lab <- c(raw = "Raw counts (no nowcast)", epinowcast = "epinowcast-corrected (default)",
           fast = "Fast delay-CDF nowcast")
  sens <- sens %>% dplyr::mutate(scn = dplyr::coalesce(unname(lab[scenario]), scenario))
  p <- ggplot(sens, aes(x = beta0, y = scn, fill = scn)) +
    ggridges::geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.72,
      colour = "grey30", linewidth = 0.25, quantile_lines = TRUE, quantiles = c(0.05, 0.5, 0.95),
      vline_colour = "grey20") +
    scale_x_log10() + scale_fill_viridis_d(option = "cividis", guide = "none") +
    labs(title = sprintf("Nowcast-input sensitivity of beta0%s",
                         if (nzchar(model_label)) paste0(" — ", model_label) else ""),
         subtitle = paste0("beta0 posterior with the SAME model fit on differently nowcast training counts. Overlap => the ",
                           "CHOICE of point nowcast does not drive the point estimate.\nCaveat: this probes the nowcast CHOICE, ",
                           "not its uncertainty (a joint model / pooling nowcast posterior draws would); and shifts partly reflect ",
                           "at-risk-set & outcome recomposition (n_events differs by scenario)."),
         x = "beta0  (log scale)", y = NULL, caption = window_txt) +
    theme_inv(11)
  if (save) .fd_save(p, file.path(OUT_REPORTS, "bayes_nowcast_sensitivity.pdf"), w = 9, h = 4.8)
  invisible(p)
}

#' Top at-risk zones by posterior mean invasion probability WITH 90% credible
#' intervals, from a Bayesian model (or the stacked ensemble). Proper posterior
#' uncertainty, coloured by province.
plot_bayes_invasion_uncertainty <- function(preds, province_map = NULL, horizon = 1L,
                                            model = NULL, top_n = 20L, window_txt = "",
                                            save = TRUE, file = "bayes_invasion_uncertainty") {
  d <- preds %>% dplyr::filter(horizon == !!horizon, !(was_active_before %in% TRUE),
                               is.finite(p_invasion))
  if (!is.null(model) && "method" %in% names(d)) d <- d %>% dplyr::filter(method == model)
  if (!nrow(d)) return(invisible(NULL))
  if (!"province" %in% names(d) && !is.null(province_map))
    d <- d %>% dplyr::left_join(province_map %>% dplyr::transmute(health_zone = nom, province),
                                by = "health_zone")
  if (!"province" %in% names(d)) d$province <- NA_character_
  d <- d %>% dplyr::mutate(region = .prov_label(province)) %>%
    dplyr::arrange(dplyr::desc(p_invasion)) %>% head(top_n) %>%
    dplyr::mutate(health_zone = factor(health_zone, levels = rev(health_zone)))
  p <- ggplot(d, aes(p_invasion, health_zone, colour = region)) +
    geom_linerange(aes(xmin = p_lo, xmax = p_hi), linewidth = 2.4, alpha = 0.4) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = .prov_colours(d$region), name = "Province") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(title = sprintf("Bayesian next-%dw invasion probability with 90%% credible intervals", horizon),
         subtitle = if (!is.null(model)) .bayes_model_label(model) else paste(unique(d$method), collapse = " / "),
         x = "Posterior P(first confirmed case)", y = NULL,
         caption = paste0("Point = posterior mean; bar = 90% credible interval. ",
                          "Per-model bands are the true posterior 90% CrI; the stacked band ",
                          "averages member quantiles and therefore UNDERSTATES between-model spread. ",
                          "At-risk zones only. ", window_txt)) +
    theme_inv(11) + theme(legend.position = "top",
                          plot.caption = element_text(colour = "grey40", size = 8, hjust = 0))
  if (save) .fd_save(p, file.path(OUT_REPORTS, sprintf("%s_h%d.pdf", file, horizon)), w = 9, h = 7)
  invisible(p)
}

#' loo predictive-stacking weights across the Bayesian assumption set.
plot_bayes_stacking <- function(weights, save = TRUE) {
  if (is.null(weights) || !length(weights)) return(invisible(NULL))
  d <- tibble::tibble(model = names(weights), weight = as.numeric(weights)) %>%
    dplyr::arrange(weight) %>% dplyr::mutate(model = factor(model, levels = model))
  p <- ggplot(d, aes(weight, model)) +
    geom_col(fill = OKABE_ITO[1], width = 0.62) +
    geom_text(aes(label = scales::percent(weight, accuracy = 1)), hjust = -0.15, size = 3.1) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, max(d$weight) * 1.15), expand = c(0, 0)) +
    labs(title = "Bayesian model-averaging weights (loo predictive stacking)",
         subtitle = "Weight each structural assumption receives in the stacked forecast (Yao et al. 2018)",
         x = "Stacking weight", y = NULL) +
    theme_inv(11)
  if (save) .fd_save(p, file.path(OUT_REPORTS, "bayes_stacking_weights.pdf"), w = 8, h = 4.2)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Task 4 — evaluation OVER TIME across folds
# ---------------------------------------------------------------------------

#' Refit the (identifiable, exogenous) covariate cloglog GLM at each fold cutoff
#' to trace how the estimated invasion drivers move as the outbreak accrues.
# Shared as-of training slice for the over-time refits (params / beta / Bayesian params).
#
# 2026-09-17: these refits used to do `zone_week_outbreak %>% filter(week_start <= cut)`,
# i.e. slice the FINAL onset-bucketed counts. That folds in cases reported only AFTER the
# fold's forecast moment and then inflates them again with the nowcast — the training-side
# revision leak that run_invasion_lfo() (16_invasion_eval.R) re-aggregates the line list to
# avoid. The over-time traces were therefore on a different footing from the LFO folds they
# are plotted against. Mirrors the LFO's fallback behaviour exactly: warn loudly, degrade to
# the final-count slice rather than erroring.
.asof_train_slice <- function(zone_week_outbreak, cut, linelist = NULL, zones_all = NULL) {
  if (!is.null(linelist) && exists("reaggregate_asof") && !is.null(zones_all)) {
    tr <- tryCatch(
      reaggregate_asof(linelist, zones_all, cut + 6,
                       week_spine = sort(unique(zone_week_outbreak$week_start))[
                         sort(unique(zone_week_outbreak$week_start)) <= cut]),
      error = function(e) {
        warning(sprintf("[over-time] cutoff %s: reaggregate_asof failed (%s); this cutoff falls back to the final-count (revision-leaky) slice.",
                        format(cut), conditionMessage(e)), call. = FALSE); NULL })
    if (!is.null(tr)) return(dplyr::filter(tr, week_start <= cut))
  }
  dplyr::filter(zone_week_outbreak, week_start <= cut)
}

plot_params_over_time <- function(pot, save = TRUE, file = "params_over_time",
                                  model_label = "", show_title = TRUE) {
  if (is.null(pot) || !nrow(pot)) return(invisible(NULL))
  # x is the FORECAST ORIGIN (cutoff + 6), not the training week's start.
  pot <- pot %>% dplyr::mutate(label = vapply(term, function(t) .lk(.PARAM_LABEL, t, t), character(1)),
                               origin = lfo_origin(cutoff))
  p <- ggplot(pot, aes(origin, hr)) +
    geom_hline(yintercept = 1, linetype = "22", colour = "grey45") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = OKABE_ITO[1], alpha = 0.18) +
    geom_line(colour = OKABE_ITO[1], linewidth = 0.7) + geom_point(size = 1.2) +
    facet_wrap(~ label, scales = "free_y") +
    scale_y_log10() +
    labs(title = if (show_title)
           sprintf("Invasion-driver estimates over time (refit at each forecast round)%s",
                   if (nzchar(model_label)) paste0(" — ", model_label) else "")
         else NULL,
         subtitle = if (show_title)
           "Hazard ratio per +1 SD with 95% CI / 90% CrI; x = forecast date (training data grows left to right)"
         else NULL,
         x = "Forecast origin (as-of date)", y = "Hazard ratio per +1 SD (log scale)",
         caption = sprintf("Identifiable exogenous drivers only (%s); shows whether/when each effect stabilises. Note: the design is re-standardised per fold, so the +1 SD unit varies slightly across cutoffs.",
                           paste(get0("BAYES_GEO_COVARIATES", ifnotfound = c("ccvi", "d_min")), collapse = ", "))) +
    theme_inv(11) + theme(plot.caption = element_text(colour = "grey40", size = 8, hjust = 0))
  if (save) .fd_save(p, file.path(OUT_DIAGNOSTICS, sprintf("%s.pdf", file)), w = 10, h = 5)
  invisible(p)
}

#' Per-fold predicted-vs-observed invasion, over time: for a method, each fold's
#' mean predicted probability vs the realised at-risk invasion fraction, plus the
#' predicted probability of the zones that DID invade. Answers "are the forecast
#' probabilities right, fold by fold?"
plot_predobs_over_folds <- function(lfo_results, method, horizon = 1L, save = TRUE,
                                    file = "predicted_vs_observed_over_folds") {
  d <- lfo_results %>% dplyr::filter(method == !!method, horizon == !!horizon,
                                     is.finite(p_invasion))
  if (!nrow(d)) return(invisible(NULL))
  by_fold <- d %>% dplyr::group_by(cutoff) %>%
    dplyr::summarise(mean_pred = mean(p_invasion),
      obs_frac = mean(is_new_invasion),
      pred_at_invaded = mean(p_invasion[is_new_invasion == 1]),
      n_atrisk = dplyr::n(), n_inv = sum(is_new_invasion), .groups = "drop") %>%
    dplyr::mutate(cutoff = as.Date(cutoff))
  long <- by_fold %>%
    dplyr::mutate(cutoff = lfo_origin(cutoff)) %>%
    dplyr::select(cutoff, `Mean predicted P` = mean_pred,
                  `Observed invasion fraction` = obs_frac,
                  `Mean P at invaded zones` = pred_at_invaded) %>%
    tidyr::pivot_longer(-cutoff)
  p <- ggplot(long, aes(cutoff, value, colour = name, shape = name)) +
    geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
    scale_colour_manual(values = c(`Mean predicted P` = OKABE_ITO[1],
      `Observed invasion fraction` = OKABE_ITO[8],
      `Mean P at invaded zones` = OKABE_ITO[2]), name = NULL) +
    scale_shape_manual(values = c(16, 15, 17), name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(title = sprintf("Predicted vs observed invasion over folds — %s (h=%dw)", method, horizon),
         subtitle = "Each point = one fold (forecast date). Well-calibrated: mean predicted ~ observed fraction; discriminating: P at invaded zones > mean.",
         x = "Forecast origin (as-of date)", y = "Probability / fraction",
         caption = "At-risk zones only, per fold; the invaded-zone line is NA where a fold had no invasion.") +
    theme_inv(11) + theme(legend.position = "top",
                          plot.caption = element_text(colour = "grey40", size = 8, hjust = 0))
  if (save) .fd_save(p, file.path(OUT_DIAGNOSTICS, sprintf("%s_h%d.pdf", file, horizon)), w = 9.5, h = 5)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Task 0b — Frequentist vs Bayesian: separation + robustness comparison
# ---------------------------------------------------------------------------

#' Tag a method name with its inferential family.
.model_family <- function(m) ifelse(grepl("^Bayes", m), "Bayesian", "Frequentist")

# plot_freq_bayes_ranking() and plot_freq_bayes_agreement() — REMOVED. They compared the two
# inferential paradigms (AUC-PR skill by family; per-zone p_case agreement with a Spearman
# correlation). With no frequentist arm the comparison is a model against itself.

#' One masked invasion choropleth per province of interest (province zoom + a
#' national context panel), so the frontier can be read within Nord-Kivu and
#' Haut-Uele exactly as for Ituri. Affected zones are greyed.
plot_province_risk_maps <- function(rs, shapefile = NULL, horizon = 1L,
                                    provinces = get0("PROVINCES_OF_INTEREST",
                                                     ifnotfound = c("Ituri", "Nord-Kivu", "Haut-Uele")),
                                    method_label = "", window_txt = "", save = TRUE,
                                    file = "invasion_risk_map") {
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile)) return(invisible(NULL))
  d <- rs %>% dplyr::filter(horizon == !!horizon)
  pcol <- if ("p_case_invasion" %in% names(d)) "p_case_invasion" else "p_invasion"
  aff_keys <- d %>% dplyr::filter(was_active_before %in% TRUE) %>%
    dplyr::pull(health_zone) %>% trimws() %>% tolower() %>% unique()
  out <- list()
  for (prov in provinces) {
    if (!prov %in% shapefile$PROVINCE) next
    lim <- c(0, max(0.05, stats::quantile(d[[pcol]], 0.99, na.rm = TRUE)))
    p_prov <- .fd_map(d, shapefile, sprintf("%s — invasion risk (%s, h=%dw)", prov, method_label, horizon),
                      prob_col = pcol, province_zoom = prov, lim = lim, palette = "plasma",
                      legend_name = "P(first case)", affected_keys = aff_keys)
    p_nat  <- .fd_map(d, shapefile, "National context", prob_col = pcol,
                      lim = lim, palette = "plasma", legend_name = "P", affected_keys = aff_keys)
    comb <- p_prov + p_nat + patchwork::plot_layout(widths = c(2, 1)) +
      patchwork::plot_annotation(caption = window_txt,
        theme = ggplot2::theme(plot.caption = element_text(colour = "grey40", size = 8, hjust = 0)))
    if (save) .fd_save(comb, file.path(OUT_MAPS,
      sprintf("%s_%s_h%d.pdf", file, .prov_suffix(prov), horizon)), w = 12, h = 6.5)
    out[[prov]] <- comb
  }
  invisible(out)
}

# ---------------------------------------------------------------------------
# Task 1 — bivariate invasion-probability x vulnerability choropleth
# ---------------------------------------------------------------------------

#' Bivariate choropleth encoding BOTH the invasion probability AND the
#' vulnerability/capacity gap on one map via a 4x4 quartile x quartile colour grid: the
#' darker-blue a zone, the more it is BOTH likely to be invaded AND poorly resourced (the
#' operational hot-corner). This shows the two dimensions the priority score
#' multiplies, without collapsing them to a single number.
#' @param show_title FALSE drops the panel title and subtitle (caption retained).
plot_prob_vuln_choropleth <- function(rs, shapefile = NULL, horizon = 1L,
                                      province_zoom = "Ituri", window_txt = "",
                                      save = TRUE, file = NULL, model_label = "",
                                      show_title = TRUE) {
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile) || !"V" %in% names(rs)) return(invisible(NULL))
  pcol <- if ("rr01_nat" %in% names(rs)) "rr01_nat" else
          if ("p_case_invasion" %in% names(rs)) "p_case_invasion" else "p_invasion"
  d <- rs %>% dplyr::filter(horizon == !!horizon, !(was_active_before %in% TRUE),
                            is.finite(.data[[pcol]]), is.finite(V))
  if (!nrow(d)) return(invisible(NULL))
  # 4x4 bivariate classification (national QUARTILES): x = invasion prob, y = vulnerability.
  NBIV <- 4L
  qt <- function(x) {
    b <- seq(0, 1, length.out = NBIV + 1L); b[1] <- -Inf; b[length(b)] <- Inf
    cut(rank(x, ties.method = "average") / length(x), breaks = b,
        labels = as.character(seq_len(NBIV)))
  }
  d <- d %>% dplyr::mutate(bx = qt(.data[[pcol]]), by = qt(V),
                           bikey = paste0(bx, "-", by))
  # 4x4 palette by BILINEAR interpolation (in RGB) of the four corners of the classic
  # teal-magenta bivariate scheme, so the hot-corner (high prob AND high vulnerability) is
  # the darkest blue and the two single-axis edges keep their distinct hues.
  .bilin_pal <- function(n, c00, c10, c01, c11) {
    g <- function(h) grDevices::col2rgb(h)[, 1] / 255
    m00 <- g(c00); m10 <- g(c10); m01 <- g(c01); m11 <- g(c11)
    keys <- character(0); cols <- character(0)
    for (bx in seq_len(n)) for (by in seq_len(n)) {
      fx <- (bx - 1) / (n - 1); fy <- (by - 1) / (n - 1)
      v <- (1 - fx) * (1 - fy) * m00 + fx * (1 - fy) * m10 +
           (1 - fx) * fy * m01 + fx * fy * m11
      v <- pmin(pmax(v, 0), 1)
      keys <- c(keys, paste0(bx, "-", by)); cols <- c(cols, grDevices::rgb(v[1], v[2], v[3]))
    }
    stats::setNames(cols, keys)
  }
  bipal <- .bilin_pal(NBIV, "#e8e8e8", "#5ac8c8", "#be64ac", "#3b4994")
  if (!is.null(province_zoom) && !province_zoom %in% shapefile$PROVINCE) return(invisible(NULL))
  has_prov <- "province" %in% names(d)
  shp <- shapefile %>% dplyr::mutate(.key = tolower(trimws(Nom)), .prov = as.character(PROVINCE))
  if (!is.null(province_zoom)) shp <- shp %>% dplyr::filter(PROVINCE == province_zoom)
  aff_keys <- rs %>% dplyr::filter(horizon == !!horizon, was_active_before %in% TRUE) %>%
    dplyr::pull(health_zone) %>% trimws() %>% tolower() %>% unique()
  # province-aware join to avoid the duplicate-name (Bili/Lubunga) many-to-many
  rj <- d %>% dplyr::transmute(.key = tolower(trimws(health_zone)),
                               .prov = if (has_prov) as.character(province) else NA_character_, bikey)
  m <- (if (has_prov) dplyr::left_join(shp, rj, by = c(".key", ".prov"))
        else dplyr::left_join(shp, dplyr::select(rj, -.prov), by = ".key")) %>%
    dplyr::mutate(affected = .key %in% aff_keys)
  main <- ggplot(m) +
    geom_sf(aes(fill = bikey), colour = "grey80", linewidth = 0.12) +
    geom_sf(data = dplyr::filter(m, affected), fill = "grey55", colour = "grey80", linewidth = 0.12) +
    scale_fill_manual(values = bipal, na.value = "grey92", guide = "none") +
    labs(title = if (show_title)
           sprintf("%s — invasion risk x vulnerability (h=%dw)%s",
                   province_zoom %||% "National", horizon,
                   if (nzchar(model_label)) paste0(" — ", model_label) else "")
         else NULL,
         subtitle = if (show_title)
           "Darker blue = both likelier to be invaded AND more vulnerable/under-resourced (national quartiles, 4x4)"
         else NULL,
         caption = window_txt) +
    ggplot2::theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey35"),
          plot.caption = element_text(colour = "grey40", size = 8, hjust = 0))
  # 4x4 legend key
  leg_df <- expand.grid(bx = seq_len(NBIV), by = seq_len(NBIV))
  leg_df$bikey <- paste0(leg_df$bx, "-", leg_df$by)
  legend <- ggplot(leg_df, aes(bx, by, fill = bikey)) +
    geom_tile() + scale_fill_manual(values = bipal, guide = "none") +
    labs(x = "Invasion prob →", y = "Vulnerability →") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 8) +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          axis.title = element_text(size = 7.5, colour = "grey30"))
  out <- main + legend + patchwork::plot_layout(widths = c(4, 1))
  fn <- file %||% sprintf("prob_vuln_choropleth_%s_h%d", .prov_suffix(province_zoom %||% "national"), horizon)
  if (save) .fd_save(out, file.path(OUT_MAPS, paste0(fn, ".pdf")), w = 11, h = 6.5)
  invisible(out)
}

# ---------------------------------------------------------------------------
# Task 6 — an intuitive skill metric for the highly imbalanced invasion task
# ---------------------------------------------------------------------------

#' Detection-vs-budget curve, pooled over LFO folds and at-risk zone-weeks. Answers
#' the operational question directly: "if we actively monitor the top-K highest-risk
#' zones each week, what fraction of the zones that actually get invaded do we
#' catch?" — i.e. sensitivity/recall at a fixed weekly alert budget, with the
#' matching precision (share of monitored zones that were truly invaded).
# 90% fold-cluster bootstrap interval for a per-fold statistic. The seed is offset by the
# budget k so neighbouring points on the curve are not driven by one shared resample (which
# would make the ribbon artificially smooth), while the whole curve stays reproducible.
.detc_boot_ci <- function(v, n_boot, seed, k) {
  v <- v[is.finite(v)]
  if (!length(v) || !is.finite(n_boot) || n_boot < 2L || length(v) < 2L)
    return(c(NA_real_, NA_real_))
  .old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(.old)) assign(".Random.seed", .old, envir = globalenv()), add = TRUE)
  set.seed(as.integer((as.numeric(seed) + 7919 * k) %% 2147483647))
  bs <- vapply(seq_len(n_boot), function(i) mean(v[sample.int(length(v), length(v), TRUE)]), numeric(1))
  unname(stats::quantile(bs, c(0.05, 0.95), names = FALSE, na.rm = TRUE))
}

#' @param n_boot fold-cluster bootstrap replicates for the recall interval (0 = none).
#'   The resampling unit is the FOLD, matching evaluate_invasion()'s cluster bootstrap and
#'   the structure of the data: zone-weeks within a fold share the same epidemic state and
#'   the same at-risk set, so a binomial (e.g. Wilson) interval on the pooled hit count
#'   treats hundreds of dependent rows as independent and is badly anticonservative. The
#'   manuscript prioritisation panel drew exactly such a Wilson ribbon.
#' @param seed bootstrap seed, so the published interval is reproducible.
#' @param common_support restrict to the (fold x zone) cells EVERY scored method covers
#'   (invasion_common_cells(), 16_invasion_eval.R). TRUE by default and it must stay that way
#'   for anything the manuscript prints: the random-targeting reference below is `k /
#'   n_atrisk`, so a method scored on its own larger row set gets a different reference line
#'   and a different denominator from the evaluation table it is printed beside.
#' @param support_cells the shared support to use, already resolved. SUPPLY THIS whenever the
#'   table passed in is not the one evaluate_invasion() scored. run_all.R appends rank-only
#'   baseline rows to a local copy (`lfo_fig3`) before building the curves, and those rows
#'   carry NA wherever the baseline has no score for a zone. Deriving the support from that
#'   copy would drop those cells for EVERY method and hand the curves a smaller support than
#'   the evaluation table -- reintroducing, from the other side, the mismatch this argument
#'   exists to prevent. Passing the cells resolved from the SCORED table pins the two together.
compute_detection_curve <- function(lfo_results, method, horizon = 1L, ks = 1:25,
                                    n_boot = 400L, seed = get0("RANDOM_SEED", ifnotfound = 20260704L),
                                    common_support = TRUE, support_cells = NULL) {
  # Resolved from the FULL table, before the single-method filter below -- the shared support
  # is a property of the whole grid, not of this method.
  .cells <- if (!is.null(support_cells)) {
    support_cells
  } else if (isTRUE(common_support)) {
    if (!exists("invasion_common_cells", mode = "function"))
      stop("[detection] common_support = TRUE but invasion_common_cells() is not loaded ",
           "(16_invasion_eval.R). Silently skipping the restriction would publish a curve ",
           "whose k / n_atrisk reference does not match the evaluation table printed beside ",
           "it. Source 16_invasion_eval.R, or pass support_cells, or set common_support = FALSE.",
           call. = FALSE)
    invasion_common_cells(lfo_results, horizon)
  } else NULL
  d <- lfo_results %>% dplyr::filter(method == !!method, horizon == !!horizon,
                                     is.finite(p_invasion))
  if (!is.null(.cells) && length(.cells)) {
    .have <- paste(d$fold_id, d$health_zone, sep = "\r")
    d <- d[.have %in% .cells, , drop = FALSE]
    # The is.finite(p_invasion) filter above runs BEFORE this restriction, so a method carrying
    # a non-finite score on a shared cell would quietly be scored on fewer cells than it was
    # handed -- a different n_atrisk, and so a different k / n_atrisk random reference, from
    # every other series on the same panel. Say so rather than let the panel draw two nulls.
    .missing <- length(.cells) - length(unique(.have[.have %in% .cells]))
    if (.missing > 0L)
      warning(sprintf(paste0("[detection] %s h=%s: %d of %d shared-support cells carry no ",
                             "finite score, so this curve is computed on a smaller row set ",
                             "than the others on the panel and its k / n_atrisk reference ",
                             "will not match theirs."),
                      method, horizon, .missing, length(.cells)), call. = FALSE)
  }
  # AT-RISK ONLY, exactly as evaluate_invasion() scores (16_invasion_eval.R: d0 drops
  # was_active_before). Without this filter the curve was computed on a LARGER row set than
  # the published metrics: already-affected zones sat in the ranking, taking top-K slots and
  # inflating the at-risk denominator behind the random-targeting reference. The manuscript
  # panels annotate a recall from this curve beside a recall_at_K from invasion_evaluation.csv,
  # so the two must be the same estimand on the same rows.
  if ("was_active_before" %in% names(d))
    d <- d[!(as.logical(d$was_active_before) %in% TRUE), , drop = FALSE]
  if (!nrow(d) || !"fold_id" %in% names(d)) return(NULL)
  per_fold <- d %>% dplyr::group_by(fold_id) %>%
    dplyr::mutate(rk = rank(-p_invasion, ties.method = "max")) %>% dplyr::ungroup()
  tot_pos <- sum(per_fold$is_new_invasion, na.rm = TRUE); nfold <- dplyr::n_distinct(per_fold$fold_id)
  if (tot_pos == 0) return(NULL)
  # RANDOM-targeting baseline (Kraemer & Cauchemez 2017, Lancet ID, Fig 3B): if you
  # monitor K zones drawn at random each fold, the expected share of invasions caught
  # is K / (mean at-risk zones per fold) — the hypergeometric expectation. This is the
  # "no-skill" reference the model must beat, and the gap is the intuitive skill story.
  n_atrisk <- per_fold %>% dplyr::count(fold_id) %>% dplyr::pull(n) %>% mean()
  purrr::map_dfr(ks, function(k) {
    inK <- per_fold$rk <= k
    tp  <- sum(per_fold$is_new_invasion[inK], na.rm = TRUE)
    # Denominator is the REALISED number of monitored zone-folds, not k * nfold. With
    # ties.method = "max" (used above, deliberately, so a zone is credited only when
    # monitoring K zones necessarily includes it) a tie group straddling the K boundary is
    # excluded, so a fold can contribute FEWER than k zones — and a fold with fewer than k
    # at-risk zones always does. Dividing by k * nfold therefore understated precision for
    # a reason unrelated to the model. (Same defect as prec_at_k in 19_spacetime_eval.)
    n_mon <- sum(inK, na.rm = TRUE)
    # TWO recall estimators, because they are genuinely different numbers and the panels
    # print one of them next to the published one:
    #   recall_pooled = all hits / all invasions, over the pooled fold-zone rows. A fold with
    #     many invasions dominates it.
    #   recall       = the per-fold share, AVERAGED over folds — the estimator
    #     .ranking_metrics() uses, so recall at k = 5/10/15 EQUALS the published
    #     recall_at_5/10/15 in invasion_evaluation.csv. This is the one to plot and annotate.
    # The curves used to carry only the pooled version while the tables published the
    # averaged one, under the same name ("share of true invasions caught").
    .per_fold_recall <- vapply(split(seq_len(nrow(per_fold)), per_fold$fold_id), function(ix) {
      yv <- per_fold$is_new_invasion[ix]; np <- sum(yv, na.rm = TRUE)
      if (np == 0) return(NA_real_)
      sum(yv[inK[ix]], na.rm = TRUE) / np
    }, numeric(1))
    .ci <- .detc_boot_ci(.per_fold_recall, n_boot, seed, k)
    tibble::tibble(method = method, horizon = horizon, k = k,
                   recall = mean(.per_fold_recall, na.rm = TRUE),
                   recall_lo = .ci[1], recall_hi = .ci[2],
                   recall_pooled = tp / tot_pos,
                   precision = if (n_mon > 0) tp / n_mon else NA_real_,
                   n_monitored = n_mon,
                   caught = tp / nfold, n_events = tot_pos / nfold,
                   n_folds = nfold, n_atrisk_mean = n_atrisk,
                   recall_random = pmin(k / n_atrisk, 1))
  })
}

#' NAIVE COMPARATOR — rank at-risk zones purely by their mobility INFLOW FROM THE
#' EPICENTRE, ignoring the case data entirely. A structural, incidence-free baseline
#' the mobility-informed renewal model must beat on the prioritisation (detection)
#' curve: score_i = sum_{e in epicentre} src_e * W[e, i], where W[e, i] is the fraction
#' of epicentre zone e's outflow that reaches i (the SAME outflow convention compute_foi
#' uses: Lambda_i = sum_j W[j, i] * ...), and src_e weights each epicentre source by its
#' population (traveller volume) when pop_vec is supplied, else uniformly. Higher score =
#' more strongly connected to the outbreak origin = naively higher invasion risk.
#'
#' @return named numeric vector over zones_all (0 for zones with no epicentre
#'   connectivity or outside W; epicentre zones themselves set to 0 — they are the
#'   already-affected origin, never an at-risk invasion target).
#' Resolve the epicentre origin zones against a matrix's row names, harmonising spellings.
#'
#' SHARED by every epicentre-anchored structural baseline so they cannot resolve the origin
#' set differently. A raw intersect() silently drops any origin given in a non-canonical
#' spelling (the pre-2026-07 "Mongbalu" for "Mongbwalu" did exactly that), which would change
#' one baseline's origins but not another's — and the three baselines are only comparable
#' because they share an origin set.
#'
#' @return character vector of resolved origins present in `row_names`, or character(0).
.resolve_epicentre_origins <- function(row_names, epicentre_zones, what = "baseline") {
  req <- unique(trimws(as.character(epicentre_zones)))
  can <- req
  if (!all(req %in% row_names) && exists("harmonise_names", mode = "function")) {
    al <- tryCatch(load_aliases(), error = function(e) NULL)
    if (!is.null(al)) can <- suppressWarnings(harmonise_names(req, al, row_names))
  }
  ok <- can %in% row_names
  if (any(!ok))
    warning(sprintf("[%s] %d origin zone(s) absent from the matrix and excluded: %s",
                    what, sum(!ok), paste(req[!ok], collapse = ", ")), call. = FALSE)
  unique(can[ok])
}

#' Rank zones by ROAD TRAVEL TIME from the epicentre — a pure-geography structural baseline.
#'
#' score_i = 1 / (1 + min_{e in epicentre} t[e, i]), with t the OSRM travel-time matrix in
#' minutes. Higher = quicker to reach from the epicentre = naively higher invasion risk.
#' The minimum (not a population-weighted sum) is the right aggregation: a zone's exposure is
#' governed by its NEAREST epicentre seed, and travel time is a cost, not a flow to be summed.
#'
#' Renewal-free and case-free by construction: it reads no incidence at all. Epicentre zones
#' are set to 0 — they are the origin, never an at-risk target — exactly as
#' naive_epicentre_inflow_scores() does, so the three baselines share a support.
#'
#' @return named numeric vector over zones_all.
epicentre_travel_time_scores <- function(osrm_mat, epicentre_zones, zones_all = rownames(osrm_mat)) {
  if (is.null(osrm_mat) || is.null(zones_all)) return(NULL)
  ez <- .resolve_epicentre_origins(rownames(osrm_mat), epicentre_zones, "travel-time")
  if (!length(ez)) { warning("[travel-time] no epicentre zones in the OSRM matrix."); return(NULL) }
  sub <- osrm_mat[ez, , drop = FALSE]
  tmin <- apply(sub, 2L, function(z) { z <- z[is.finite(z)]; if (length(z)) min(z) else NA_real_ })
  sc <- 1 / (1 + tmin)
  sc[!is.finite(sc)] <- 0                     # unroutable = no connectivity, not NA
  names(sc) <- colnames(osrm_mat)
  out <- stats::setNames(rep(0, length(zones_all)), zones_all)
  common <- intersect(zones_all, names(sc))
  out[common] <- sc[common]
  out[intersect(zones_all, ez)] <- 0
  out
}

naive_epicentre_inflow_scores <- function(W, epicentre_zones, pop_vec = NULL,
                                          zones_all = rownames(W)) {
  if (is.null(W) || is.null(zones_all)) return(NULL)
  ez <- .resolve_epicentre_origins(rownames(W), epicentre_zones, "naive")
  if (!length(ez)) { warning("[naive] no epicentre zones present in the mobility matrix."); return(NULL) }
  # Source weight per epicentre zone = its population (relative traveller volume), median-
  # imputed for any missing/non-positive entry; uniform when no pop_vec is provided.
  src_w <- if (!is.null(pop_vec)) {
    p <- as.numeric(pop_vec[match(ez, names(pop_vec))])
    good <- is.finite(p) & p > 0
    if (!any(good)) rep(1, length(ez)) else { p[!good] <- stats::median(p[good]); p }
  } else rep(1, length(ez))
  # inflow to each destination = sum_e src_w[e] * W[e, ]  (t(src_w) %*% W[ez, ]).
  sub   <- W[ez, , drop = FALSE]
  score <- as.numeric(crossprod(src_w, sub))          # length = ncol(W)
  names(score) <- colnames(W)
  out    <- stats::setNames(rep(0, length(zones_all)), zones_all)
  common <- intersect(zones_all, names(score))
  out[common] <- score[common]
  out[intersect(zones_all, ez)] <- 0                  # epicentre = origin, not an at-risk target
  out
}

#' Inject the naive epicentre-inflow ranking into an lfo_results table as an extra
#' "method", so plot_detection_curve / plot_paper_figure3 draw its prioritisation
#' curve alongside the real models. The naive rows REUSE the exact fold structure and
#' ground-truth invasion outcomes of an existing method (the at-risk (fold, zone) set
#' and is_new_invasion are model-independent), overwriting ONLY the ranking key
#' p_invasion with the static per-zone inflow score — a fair, like-for-like comparison
#' on identical folds. Idempotent; a no-op if the inputs are empty. The naive model
#' carries a RANKING only (no calibrated probability / count), so the calibration and
#' count columns are blanked to NA — the detection/precision curves use only the rank.
append_naive_detection_curve_model <- function(lfo_results, scores,
                                               label = "Naive-epicentre-inflow") {
  if (is.null(lfo_results) || !nrow(lfo_results) || is.null(scores)) return(lfo_results)
  if (label %in% lfo_results$method) return(lfo_results)          # idempotent
  # The row SKELETON these baseline scores are hung on used to be `unique(method)[1]` -- i.e.
  # whichever method happened to sort first in the table. That silently set the baseline's
  # fold coverage from row order: pick a method short of a fold (the rolling-predictor floor
  # costs every Bayesian model the earliest one) and the baseline inherits the gap, so the
  # figure's structural nulls would be drawn on fewer cells than the model they are there to
  # beat. Take the method with the widest (fold x zone) coverage instead, tie-broken by name
  # so the choice does not move with row order. The curve is restricted to the shared support
  # afterwards, so a wider skeleton costs nothing and a narrower one cannot be recovered.
  .ref_method <- if (all(c("fold_id", "health_zone") %in% names(lfo_results))) {
    .cov <- lfo_results %>%
      dplyr::distinct(method, fold_id, health_zone) %>%
      dplyr::count(method, name = "n_cells") %>%
      dplyr::arrange(dplyr::desc(n_cells), method)
    .cov$method[1]
  } else unique(lfo_results$method)[1]
  ref <- lfo_results %>% dplyr::filter(method == .ref_method)
  if (!nrow(ref) || !"health_zone" %in% names(ref)) return(lfo_results)
  naive <- ref
  naive$method     <- label
  # ZERO-FILL UNOBSERVED DESTINATIONS. A zone the mobility source never observed has no
  # measured inflow, and zero is what this baseline would have told a user on the day: it
  # ranks such a zone at the bottom. Leaving NA instead DROPS those rows, putting the
  # baseline on a narrower support than the models -- and because the shared support is
  # anchored on every scored method, that would drag every model down to the baseline's
  # coverage. The fill is faithful to the comparator, but it is REPORTED rather than silent:
  # a large zero block is a property of the DATA SOURCE, not of the baseline's skill, and
  # zeros tie at an averaged rank, which depresses mean_rank_of_truth in proportion to the
  # censoring. Read this coverage line beside any baseline's rank metrics.
  .sc_raw <- as.numeric(scores[match(naive$health_zone, names(scores))])
  .n_zone <- dplyr::n_distinct(naive$health_zone)
  .n_obs  <- dplyr::n_distinct(naive$health_zone[is.finite(.sc_raw) & .sc_raw > 0])
  naive$p_invasion <- ifelse(is.finite(.sc_raw), .sc_raw, 0)
  message(sprintf("[baseline] %s: %d/%d zones observed (%.1f%%), %d zero-filled.",
                  label, .n_obs, .n_zone, 100 * .n_obs / max(.n_zone, 1L), .n_zone - .n_obs))
  if ("p_case_invasion" %in% names(naive)) naive$p_case_invasion <- naive$p_invasion
  for (col in intersect(c("mu_forecast", "p_lo", "p_hi", "p_sd"),
                        names(naive)))
    naive[[col]] <- NA_real_
  # The recalibration columns must be blanked for the SAME reason: this row set is a copy of
  # another method's rows with only the score replaced, so it would otherwise inherit that
  # method's prequential factor and its RECALIBRATED probability. Any panel drawn on the
  # recalibrated column would then plot the model's own corrected probabilities under the
  # baseline's name. The naive score is a mobility inflow, not a probability, so there is
  # nothing to recalibrate: p_recal is set to the score itself (identical on both scales,
  # which is correct for a rank-only comparator) and the fitted-factor bookkeeping is NA.
  if ("p_recal" %in% names(naive)) naive$p_recal <- naive$p_invasion
  # RANK-ONLY by construction: these rows carry a mobility inflow or an inverse travel time,
  # not a probability. Without this they inherit prob_calibrated from the method whose rows
  # were copied (TRUE for any Bayesian model), and evaluate_invasion() would compute a log
  # score, Brier skill and calibration-in-the-large on a connectivity score — the same defect
  # already fixed for Distance-B1 and Adjacency-B7 in 05_baseline_models.R.
  if ("prob_calibrated" %in% names(naive)) naive$prob_calibrated <- FALSE
  for (col in intersect(c("delta_preq", "n_train_events", "n_train_folds"), names(naive)))
    naive[[col]] <- NA_real_
  if ("delta_estimable" %in% names(naive)) naive$delta_estimable <- NA
  dplyr::bind_rows(lfo_results, naive)
}

#' Single-number "balanced" skill at the Youden-optimal threshold — balanced
#' accuracy (mean of sensitivity & specificity), F1, and Matthews correlation.
#' Unlike raw accuracy or AUC-ROC these are not inflated by the ~99.7% of zone-weeks
#' with no invasion, so they are honest for a rare binary event.
invasion_balance_metrics <- function(lfo_results, method, horizon = 1L) {
  d <- lfo_results %>% dplyr::filter(method == !!method, horizon == !!horizon,
                                     is.finite(p_invasion))
  if (!nrow(d) || sum(d$is_new_invasion, na.rm = TRUE) == 0) return(NULL)
  y <- as.integer(d$is_new_invasion); p <- d$p_invasion
  best <- NULL
  for (t in sort(unique(p))) {
    pred <- p >= t
    tp <- sum(pred & y == 1, na.rm = TRUE); fp <- sum(pred & y == 0, na.rm = TRUE)
    fn <- sum(!pred & y == 1, na.rm = TRUE); tn <- sum(!pred & y == 0, na.rm = TRUE)
    sens <- tp / (tp + fn); spec <- tn / (tn + fp); j <- sens + spec - 1
    if (is.null(best) || j > best$j)
      best <- list(t = t, sens = sens, spec = spec, j = j, tp = tp, fp = fp, fn = fn, tn = tn)
  }
  b <- best; prec <- b$tp / max(b$tp + b$fp, 1); rec <- b$sens
  f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  # coerce to double BEFORE multiplying: with ~10^5 pooled at-risk rows the 4-way
  # integer product overflows R's 32-bit integer (-> NA -> crashes the guard).
  tp <- as.double(b$tp); fp <- as.double(b$fp); fn <- as.double(b$fn); tn <- as.double(b$tn)
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (is.finite(den) && den > 0) (tp * tn - fp * fn) / den else 0
  tibble::tibble(method = method, horizon = horizon, threshold = b$t,
                 sensitivity = b$sens, specificity = b$spec, youden_j = b$j,
                 balanced_accuracy = (b$sens + b$spec) / 2, f1 = f1, mcc = mcc)
}

#' The detection-vs-budget curve for one or more methods — the intuitive headline
#' figure: monitor K zones (x), catch this fraction of invasions (y).
#' @param support_cells the shared (fold x zone) support, already resolved. PASS IT whenever
#'   `lfo_results` is not the table evaluate_invasion() scored -- run_all.R draws this panel
#'   from `lfo_fig3`, which carries appended rank-only baseline rows, and letting
#'   compute_detection_curve() derive the support from that copy would put this panel's
#'   random-targeting reference on a different denominator from the published curve.
plot_detection_curve <- function(lfo_results, methods, horizon = 1L, save = TRUE,
                                 annotate_k = 10L, support_cells = NULL) {
  curves <- purrr::map_dfr(methods, function(m)
    compute_detection_curve(lfo_results, m, horizon, support_cells = support_cells))
  if (is.null(curves) || !nrow(curves)) return(invisible(NULL))
  mth <- unique(curves$method)
  pal <- if (length(mth) <= length(OKABE_ITO)) OKABE_ITO[seq_along(mth)]
         else grDevices::hcl.colors(length(mth), "Dark 3")
  # the RANDOM-targeting reference line (Kraemer & Cauchemez 2017 Fig 3B) — identical
  # across methods, so take it from the first curve.
  rnd <- curves %>% dplyr::filter(method == mth[1]) %>% dplyr::distinct(k, recall_random)
  # headline callout: best model's catch-rate vs random at the chosen budget K.
  best_at_k <- curves %>% dplyr::filter(k == annotate_k) %>% dplyr::slice_max(recall, n = 1)
  cap <- if (nrow(best_at_k))
    sprintf("At a budget of K=%d zones/week, %s catches %.0f%% of the next invasions vs %.0f%% for random targeting (%.1fx better).",
            annotate_k, best_at_k$method[1], 100 * best_at_k$recall[1],
            100 * best_at_k$recall_random[1],
            best_at_k$recall[1] / max(best_at_k$recall_random[1], 1e-6))
  else "A good model reaches a high catch-rate at small K; the dashed line is random targeting."
  p <- ggplot(curves, aes(k, recall)) +
    geom_line(data = rnd, aes(k, recall_random), linetype = "22", colour = "grey45", linewidth = 0.7) +
    geom_line(aes(colour = method), linewidth = 0.9) + geom_point(aes(colour = method), size = 1) +
    { if (annotate_k %in% curves$k) geom_vline(xintercept = annotate_k, colour = "grey80", linewidth = 0.3) } +
    annotate("text", x = max(rnd$k) * 0.62, y = max(rnd$recall_random) * 0.7 + 0.03,
             label = "random targeting", colour = "grey45", size = 3, angle = 12) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_colour_manual(values = setNames(pal, mth), name = NULL) +
    labs(title = sprintf("If we monitor the top-K highest-risk zones, how many invasions do we catch? (h=%dw)", horizon),
         subtitle = "Share of true next-week invasions caught at a fixed weekly alert budget of K zones, averaged over LFO folds",
         x = "Zones actively monitored each week (K)", y = "Share of true invasions caught",
         caption = cap) +
    theme_inv(11) + theme(legend.position = "top",
                          plot.caption = element_text(colour = "grey30", size = 8.5, hjust = 0, face = "italic"))
  if (save) .fd_save(p, file.path(OUT_DIAGNOSTICS, sprintf("detection_vs_budget_h%d.pdf", horizon)), w = 9, h = 5.6)
  invisible(p)
}

#' The PRECISION counterpart of the detection curve (#3): of the top-K highest-risk
#' zones you monitor each week, what FRACTION were actually invaded? — "how many of
#' our alerts are true". Pooled over LFO folds, with the base-rate reference (a
#' random watch-list would hit invasions only at the ~invasion prevalence).
#' @param support_cells as in plot_detection_curve(): supply it when `lfo_results` is not the
#'   table evaluate_invasion() scored.
plot_topk_precision <- function(lfo_results, methods, horizon = 1L, save = TRUE,
                                file = "topk_precision", title_suffix = "",
                                support_cells = NULL) {
  curves <- purrr::map_dfr(methods, function(m)
    compute_detection_curve(lfo_results, m, horizon, support_cells = support_cells))
  if (is.null(curves) || !nrow(curves)) return(invisible(NULL))
  mth <- unique(curves$method)
  pal <- if (length(mth) <= length(OKABE_ITO)) OKABE_ITO[seq_along(mth)]
         else grDevices::hcl.colors(length(mth), "Dark 3")
  # Base rate = invasions per at-risk zone-week = a random watch-list's hit rate. Taken
  # DIRECTLY from the curve's own n_events / n_atrisk_mean. It used to be recovered
  # algebraically by inverting recall_random (= k / n_atrisk) and averaging over k, which
  # is the same quantity only while recall_random is unclipped — pmin(k / n_atrisk, 1)
  # saturates at 1 for k >= n_atrisk, so every saturated k contributed n_events/k instead
  # and dragged the averaged "base rate" DOWNWARD, flattering precision against it.
  base_rate <- suppressWarnings(mean(
    (curves$n_events / curves$n_atrisk_mean)[is.finite(curves$n_atrisk_mean) & curves$n_atrisk_mean > 0],
    na.rm = TRUE))
  p <- ggplot(curves, aes(k, precision, colour = method)) +
    { if (is.finite(base_rate)) geom_hline(yintercept = base_rate, linetype = "22", colour = "grey55") } +
    geom_line(linewidth = 0.9) + geom_point(size = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
    scale_colour_manual(values = setNames(pal, mth), name = NULL) +
    labs(title = sprintf("Of the top-K highest-risk zones, what %% were actually invaded? (h=%dw)%s",
                         horizon, if (nzchar(title_suffix)) paste0(" — ", title_suffix) else ""),
         subtitle = "Precision at a fixed weekly alert budget of K zones, pooled over LFO fold-zone rows; dashed = random-watch-list base rate",
         x = "Zones monitored each week (K)", y = "% of top-K zones that were invaded",
         caption = "Higher = fewer false alarms. Precision declines with K as the budget reaches lower-risk zones; the gap above the base rate is the model's targeting value.") +
    theme_inv(11) + theme(legend.position = "top",
                          plot.caption = element_text(colour = "grey30", size = 8.5, hjust = 0, face = "italic"))
  if (save) .fd_save(p, file.path(OUT_DIAGNOSTICS, sprintf("%s_h%d.pdf", file, horizon)), w = 9, h = 5.6)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Task 4 — import coefficient over folds (PLOTTERS ONLY)
# ---------------------------------------------------------------------------
# The frequentist COMPUTERS that used to live here — compute_params_over_time(),
# compute_beta_over_folds() and plot_reporting_rate_map() — have been removed. The first two
# refitted a Firth cloglog GLM at every fold cutoff to draw figures that FIGURE_KEEP then
# gated out of the published tree; the third mapped a "reporting-rate proxy the model uses"
# that no surviving model uses. The plotters below are retained because the BAYESIAN traces
# (compute_bayes_params_over_time) feed them, and those figures ARE published — now with
# key_outputs/bayes_{params_over_time,beta_over_folds}.csv behind them.

plot_beta_over_folds <- function(bof, save = TRUE, file = "beta_and_completeness_over_folds",
                                 model_label = "", ci_label = "95% CI",
                                 show_title = TRUE, base_size = 11, png = FALSE) {
  if (is.null(bof) || !nrow(bof)) return(invisible(NULL))
  .ml <- if (nzchar(model_label)) paste0(" — ", model_label) else ""
  .s  <- base_size / 11
  # Plot against the FORECAST ORIGIN (cutoff + 6) -- what the axis label claims to show.
  bof$cutoff <- lfo_origin(bof$cutoff)
  p1 <- ggplot(bof, aes(cutoff, beta0)) +
    geom_ribbon(aes(ymin = beta0_lo, ymax = beta0_hi), fill = OKABE_ITO[1], alpha = 0.18) +
    geom_line(colour = OKABE_ITO[1], linewidth = 0.7 * .s) + geom_point(size = 1.4 * .s) +
    scale_y_log10() +
    labs(title = if (show_title) paste0("Import coefficient beta0 over folds", .ml) else NULL,
         subtitle = if (show_title)
           sprintf("exp(intercept) of the renewal fit, refit at each cutoff (%s)", ci_label) else NULL,
         x = "Forecast origin (as-of date)", y = "beta0 = P-scale import->invasion hazard") +
    theme_inv(base_size)
  # completeness panel only when it is populated (frequentist path); the Bayesian
  # beta0-over-folds has no completeness column, so show beta0 alone.
  out <- if ("mean_completeness" %in% names(bof) && any(is.finite(bof$mean_completeness))) {
    p2 <- ggplot(bof, aes(cutoff, mean_completeness)) +
      geom_line(colour = OKABE_ITO[8], linewidth = 0.7 * .s) + geom_point(size = 1.4 * .s) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(title = if (show_title) "Mean reporting completeness over folds" else NULL,
           subtitle = if (show_title)
             "Observed / nowcast-corrected counts in the recent weeks at each cutoff" else NULL,
           x = "Forecast origin (as-of date)", y = "Reporting completeness") +
      theme_inv(base_size)
    p1 / p2
  } else p1
  if (save) {
    .h <- if (inherits(out, "patchwork")) 8 else 4.5
    .fd_save(out, file.path(OUT_DIAGNOSTICS, sprintf("%s.pdf", file)), w = 9, h = .h)
    if (png) .fd_save_png(out, file.path(OUT_DIAGNOSTICS, sprintf("%s.png", file)), w = 9, h = .h)
  }
  invisible(out)
}

#' Invasion-probability map WITH its uncertainty side by side (#2): (left) the mean
#' P(first case) and (right) the uncertainty WIDTH (p_hi - p_lo = the 90% credible /
#' ensemble interval), so where the forecast is confident vs uncertain is spatially
#' legible. Requires p_lo/p_hi (Bayesian posterior or ensemble spread). Ituri + national.
plot_invasion_uncertainty_map <- function(rs, shapefile = NULL, horizon = 1L,
                                          model_label = "", window_txt = "", save = TRUE,
                                          file = "invasion_uncertainty_map",
                                          extent = c("ituri", "national")) {
  extent <- match.arg(extent)
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile) || !all(c("p_lo", "p_hi") %in% names(rs))) return(invisible(NULL))
  d <- rs %>% dplyr::filter(horizon == !!horizon)
  pcol <- if ("p_case_invasion" %in% names(d)) "p_case_invasion" else "p_invasion"
  d <- d %>% dplyr::mutate(.width = p_hi - p_lo)
  if (!any(is.finite(d[[pcol]]))) return(invisible(NULL))
  aff <- d %>% dplyr::filter(was_active_before %in% TRUE) %>%
    dplyr::pull(health_zone) %>% trimws() %>% tolower() %>% unique()
  it <- (extent == "ituri"); reg <- if (it) "Ituri" else "DRC (national)"
  lim_p <- c(0, max(0.05, stats::quantile(d[[pcol]], 0.99, na.rm = TRUE)))
  lim_w <- c(0, max(0.05, stats::quantile(d$.width, 0.99, na.rm = TRUE)))
  p_prob <- .fd_map(d, shapefile, sprintf("%s — P(first case), h=%dw", reg, horizon),
                    prob_col = pcol, ituri_only = it, lim = lim_p, palette = "plasma",
                    legend_name = "P(first case)", affected_keys = aff)
  p_unc  <- .fd_map(d, shapefile, sprintf("%s — forecast uncertainty (90%% interval width)", reg),
                    prob_col = ".width", ituri_only = it, lim = lim_w, palette = "viridis",
                    legend_name = "Interval width", affected_keys = aff)
  out <- p_prob + p_unc + patchwork::plot_layout(widths = c(1, 1)) +
    patchwork::plot_annotation(
      title = sprintf("%s invasion probability with uncertainty%s (h=%dw)", reg,
                      if (nzchar(model_label)) paste0(" — ", model_label) else "", horizon),
      subtitle = "Left: mean P(first confirmed case). Right: width of the 90% credible/ensemble interval (darker = more uncertain).",
      caption = window_txt,
      theme = ggplot2::theme(plot.title = element_text(face = "bold", size = 13),
                             plot.subtitle = element_text(size = 9, colour = "grey35"),
                             plot.caption = element_text(colour = "grey40", size = 8, hjust = 0)))
  if (save) .fd_save(out, file.path(OUT_MAPS, sprintf("%s_h%d.pdf", file, horizon)),
                     w = if (it) 12 else 13, h = 6.5)
  invisible(out)
}

# ---------------------------------------------------------------------------
# ANIMATION VIZ — REMOVED
# ---------------------------------------------------------------------------
# .frames_to_gif(), make_horizon_animation(), animate_invasion_over_time(),
# animate_invasion_rank_over_time() and plot_mobility_flows_over_time() rendered GIF /
# flipbook animations of the invasion maps, the rank evolution and the mobility-routed import
# force. NOTHING CALLED THEM: run_all.R's only reference is the cleanup block that DELETES
# outputs/maps/animations/ so an earlier run's GIFs cannot look current, and the comment there
# records why ("they duplicated the static decision products while dominating the
# visualisation runtime"). ~200 lines of unreachable rendering, carrying a soft gifski
# dependency, in the module that also holds the detection curve and the risk-score producers.
# The deletion sweep in run_all.R stays: it is what keeps stale GIFs from a pre-removal run
# out of the published tree.

message("[forecast_detail] 20_forecast_detail.R loaded (spec/selection + ensemble/uncertainty + params/priority + bayes/over-time + province/bivariate + detection + topk + reporting/beta + uncertainty-map).")
