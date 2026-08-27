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

if (!exists("%||%"))
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 ||
                              (length(a) == 1 && is.na(a))) b else a
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
  if (file.exists(path)) suppressWarnings(file.remove(path))
  tryCatch(ggplot2::ggsave(path, p, width = w, height = h, device = "pdf",
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
  M4  = "M4 gravity — negative-binomial gravity GLM on the national origin-destination matrix (power-law distance deterrence)",
  M4b = "M4b gravity — as M4 but with exponential distance deterrence",
  M9  = "M9 multi-kernel — average of gravity, radiation and travel-time-decay kernels",
  M10 = "M10 radiation-composite — parameter-free radiation model blended with epicentre short-trip flows",
  M3  = "M3 — raw Flowminder national origin-destination flows",
  M5  = "M5 radiation model (Simini et al. 2012)",
  M6a = "M6a travel-time exponential-decay kernel",
  M6b = "M6b travel-time power-law-decay kernel")
.GT_DESC <- c(
  medium = "GT-Medium (mean 15.3 d, sd 9.3 d; Zaire EBOV serial interval, WHO Ebola Response Team 2014 NEJM)",
  short  = "GT-Short (mean 12.0 d, sd 6.5 d; low-end sensitivity, below Nash 2024 pooled 95% CI)",
  long   = "GT-Long (mean 18.0 d, sd 10.5 d; high-end sensitivity, above Nash 2024 pooled 95% CI)")
.COV_DESC <- c(
  log_pop      = "log population size",
  ccvi         = "CCVI socioeconomic vulnerability index",
  positivity   = "zone test positivity",
  d_min        = "travel time to the nearest currently-affected zone",
  alert_import = "mobility-weighted suspected-alert pressure imported from other zones",
  alert_local  = "local suspected-alert count",
  susp_import  = "mobility-weighted suspected-but-not-confirmed case pressure imported from other zones",
  susp_local   = "local preceding-week suspected-but-not-confirmed case count",
  log_healthsite_count = "healthcare-site density")
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
      nowcast = "epinowcast-corrected members",
      notes = "pre-specified member set; no selection on the test folds"))
  }
  # Bayesian renewal grid (Bayes-<mob>-<gt/cov/link>): same mechanistic model as the
  # frequentist renewal, fit with brms — so it must be DESCRIBED, not dumped into the
  # "Comparator / NA" bucket (the featured Bayesian model routinely wins overall).
  if (grepl("^Bayes", m)) {
    mobn <- c(M1 = "short-trip", M4 = "gravity", M8 = "composite-gravity",
              M9 = "multi-kernel ensemble", M10 = "radiation-composite",
              M11 = "inward meeting-location FOI")
    cov_desc <- if (grepl("-full", m)) "full exogenous set (log-pop, CCVI, health-site density, travel-time-to-affected)"
                else if (grepl("-geo", m)) "geographic exogenous set (log-pop, CCVI, travel-time-to-affected)"
                else "none (intercept-only import hazard)"
    lnk <- if (grepl("-logit", m)) "Logit (log-odds; beta is not a hazard multiplier)"
           else "Complementary log-log with log(Lambda) offset (renewal hazard: p = 1 - exp(-beta*Lambda))"
    return(list(
      family = "Bayesian mobility-informed renewal (brms; posterior)",
      # mob_key strips the "-dist" suffix for the description lookup, so re-attach a
      # road-distance note for the OSRM-km kernels (their offset uses km deterrence, not travel time).
      mobility = {
        .mob_base <- .lk(mobn, mob_key, mob_key %||% "n/a")
        if (grepl("-dist", m)) paste0(.mob_base, " (OSRM road-distance deterrence)") else .mob_base
      },
      gt = gt_key, observation = lnk, covariates = cov_desc,
      calibration = "posterior; weakly-informative priors (Intercept~Normal(-3,2), coefs~Normal(0,1)); loo predictive stacking for the ensemble",
      nowcast = "epinowcast-corrected training counts",
      notes = "Bayesian analogue of the frequentist renewal model; featured Bayesian model = highest loo predictive-stacking weight"))
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
    nowcast = vs$nowcast %||% "epinowcast right-truncation correction",
    notes = vs$note %||% NA_character_)
}

#' Markdown paragraph describing HOW the best model is chosen.
describe_selection <- function() {
  c("## How the models are selected", "",
    "Every model is scored in the **same** leave-future-out cross-validation on the",
    "at-risk zones only (no already-affected zone ever counts), pooled across folds.",
    "Selection is **discrimination-led with a spiky-model gate** (not by any single",
    "metric, and never by count-WIS, which would reward predicting zero everywhere):",
    "",
    "1. **AUC-PR skill** — Average Precision divided by the invasion base rate:",
    "   how well the model concentrates the few invasions at the top of its ranking;",
    "2. **Mean rank of the invaded zone** — how near the top the zones that actually",
    "   invaded were placed (operational targeting);",
    "3. **Log-score** — a proper score that penalises over-confident probabilities.",
    "",
    "Among models present at BOTH horizons, any whose worst-horizon mean rank-of-truth",
    "exceeds 1.5x the field median is first dropped as **spiky** (high AUC-PR driven by a",
    "few top hits while the rest are scattered). The winner is the highest **total AUC-PR",
    "skill** (pooled across horizons) among the survivors; ties are broken by lower total",
    "rank-of-truth, then lower log-score. Requiring both horizons stops a model winning by",
    "performing well at 1-week while sitting out the 2-week task.",
    "",
    "The **best-fitting model** is the best over ALL families on this composite;",
    "the **best renewal model** is the best *mobility-informed renewal* variant (so",
    "the featured risk product is an interpretable, mechanistic model).",
    "With only ~15 invasion events the differences between the leading models are",
    "within cross-validation noise — hence the pre-specified ensemble is reported",
    "alongside the single best model.", "")
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
  same <- identical(as.character(best_method), as.character(primary_method))
  blocks <- if (same)
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

#' Per-zone forecast uncertainty as the spread across ensemble members.
#'
#' For each (health_zone, horizon) returns the min / median / max / mean of the
#' invasion probability across the member models — the structural uncertainty
#' band. Already-affected zones (was_active_before) are dropped (no invasion
#' probability by construction). Works on the current forecast OR on one LFO
#' fold (any long tibble with method + prob_col + optional cutoff).
#'
#' @param fc long forecast tibble containing the members.
#' @param members member method names whose spread defines the band.
#' @param prob_col probability column (default "p_invasion").
#' @param by extra grouping columns (e.g. "cutoff" for the space-time version).
ensemble_member_uncertainty <- function(fc, members, prob_col = "p_invasion",
                                        by = character(0)) {
  if (!prob_col %in% names(fc))
    prob_col <- intersect(c("p_invasion", "p_case_invasion"), names(fc))[1]
  keep <- fc$method %in% members
  if ("was_active_before" %in% names(fc)) keep <- keep & !(fc$was_active_before %in% TRUE)
  d <- fc[keep, , drop = FALSE]
  d <- d[is.finite(d[[prob_col]]), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  grp <- c("health_zone", "horizon", by)
  grp <- intersect(grp, names(d))
  carry <- intersect(c("is_new_invasion", "province"), names(d))
  d %>% dplyr::rename(.p = dplyr::all_of(prob_col)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
    dplyr::summarise(
      p_lo = min(.p), p_med = stats::median(.p), p_hi = max(.p),
      p_mean = mean(.p), n_members = dplyr::n(),
      dplyr::across(dplyr::all_of(carry), ~ dplyr::first(.x)),
      .groups = "drop") %>%
    dplyr::mutate(p_width = p_hi - p_lo)
}

#' Q4 — current forecast for the top at-risk zones WITH an uncertainty band.
#' Each member model is a hollow dot, so the disagreement is visible as several
#' marks; the shaded range spans their min-max, the filled dot is the ensemble
#' median, and the black tick is the best renewal model. Members differ mainly in
#' the mobility kernel (M8 epicentre-flow vs M4 gravity vs M9/M10), so a zone one
#' member rates ~0 and others rate ~10% shows a wide, informative spread.
plot_forecast_uncertainty <- function(fc, members, primary_method = NULL,
                                       province_map = NULL, horizon = 1L,
                                       top_n = 20L, save = TRUE, window_txt = "") {
  unc <- ensemble_member_uncertainty(fc, members, "p_invasion")
  if (is.null(unc)) return(invisible(NULL))
  unc <- unc %>% dplyr::filter(horizon == !!horizon)
  if (!"province" %in% names(unc) && !is.null(province_map))
    unc <- unc %>% dplyr::left_join(
      province_map %>% dplyr::transmute(health_zone = nom, province), by = "health_zone")
  if (!"province" %in% names(unc)) unc$province <- NA_character_
  unc <- unc %>% dplyr::mutate(region = .prov_label(province))
  top <- unc %>% dplyr::arrange(dplyr::desc(p_med)) %>% head(top_n) %>%
    dplyr::mutate(health_zone = factor(health_zone, levels = rev(health_zone)))
  zlev <- levels(top$health_zone)

  # Individual member estimates for the displayed zones (one dot per member), so
  # the disagreement reads as several marks rather than one flat band.
  mem_pts <- fc %>% dplyr::filter(method %in% members, horizon == !!horizon,
                                  !(was_active_before %in% TRUE),
                                  is.finite(p_invasion), health_zone %in% zlev) %>%
    dplyr::transmute(health_zone = factor(health_zone, levels = zlev),
                     method, p = p_invasion)

  prim <- NULL
  if (!is.null(primary_method)) {
    prim <- fc %>% dplyr::filter(method == primary_method, horizon == !!horizon,
                                 !(was_active_before %in% TRUE)) %>%
      dplyr::transmute(health_zone, p_prim = p_invasion) %>%
      dplyr::filter(health_zone %in% zlev) %>%
      dplyr::mutate(health_zone = factor(health_zone, levels = zlev))
  }
  p <- ggplot(top, aes(p_med, health_zone)) +
    geom_linerange(aes(xmin = p_lo, xmax = p_hi, colour = region),
                   linewidth = 2.6, alpha = 0.22) +
    geom_point(data = mem_pts, aes(p, health_zone), shape = 21, size = 1.9,
               fill = "white", colour = "grey30", stroke = 0.5, alpha = 0.9) +
    geom_point(aes(colour = region), size = 3.1) +
    { if (!is.null(prim) && nrow(prim) > 0)
        geom_point(data = prim, aes(p_prim, health_zone),
                   shape = 124, size = 3.8, colour = "black") else NULL } +
    scale_colour_manual(values = .prov_colours(top$region), name = "Province") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.03, 0.05))) +
    labs(title = sprintf("Next-%dw invasion probability with ensemble uncertainty", horizon),
         subtitle = paste0("Hollow dots = ", length(members),
                           " member models | filled = ensemble median | range = min-max",
                           if (!is.null(prim)) " | tick = best renewal model" else ""),
         x = "P(first confirmed case)", y = NULL,
         caption = paste0("Uncertainty = disagreement across pre-specified members (differing mobility kernel / GT / observation); at-risk zones only. ", window_txt)) +
    theme_inv(11) + theme(legend.position = "top",
                          plot.caption = element_text(colour = "grey40", size = 8))
  if (save) .fd_save(p, file.path(OUT_REPORTS, sprintf("forecast_uncertainty_h%d.pdf", horizon)),
                     w = 9.5, h = 7)
  invisible(p)
}

#' Q4 — space-AND-time forecast trajectories WITH uncertainty.
#' For the zones that actually invaded during the CV window (the operationally
#' informative stories), plot the predicted invasion probability across fold
#' cutoffs (time) as a median line + member-spread ribbon, one panel per zone
#' (space), with the realised invasion week marked. Shows whether — and how
#' confidently — each first case was anticipated.
plot_spacetime_forecast <- function(lfo_results, members, province_map = NULL,
                                     horizon = 1L, max_zones = 12L, save = TRUE) {
  if (is.null(lfo_results) || !"cutoff" %in% names(lfo_results)) return(invisible(NULL))
  unc <- ensemble_member_uncertainty(
    lfo_results %>% dplyr::filter(horizon == !!horizon), members,
    "p_invasion", by = "cutoff")
  if (is.null(unc) || nrow(unc) == 0) return(invisible(NULL))
  unc <- unc %>% dplyr::mutate(cutoff = as.Date(cutoff))

  # Prefer zones that invaded during the window AND have a real (>=2-week)
  # forecast lead-up, so every panel shows a trajectory climbing toward the
  # invasion mark rather than a lone point. Rank by number of pre-invasion
  # observations, then peak risk; fall back to highest-risk zones if none.
  stats <- unc %>% dplyr::group_by(health_zone) %>%
    dplyr::summarise(n_cut = dplyr::n_distinct(cutoff),
                     invaded = any(is_new_invasion %in% 1),
                     peak = max(p_med, na.rm = TRUE), .groups = "drop")
  multi <- stats %>% dplyr::filter(invaded, n_cut >= 2) %>%
    dplyr::arrange(dplyr::desc(peak))
  if (nrow(multi) >= 4) {
    invaded <- multi %>% head(max_zones) %>% dplyr::pull(health_zone)
  } else {
    # too few multi-week stories: take invaded zones by lead-up length then risk,
    # or (no invasions at all) the highest-risk at-risk zones.
    pool <- stats %>% dplyr::filter(invaded)
    if (nrow(pool) == 0) pool <- stats
    invaded <- pool %>% dplyr::arrange(dplyr::desc(n_cut), dplyr::desc(peak)) %>%
      head(max_zones) %>% dplyr::pull(health_zone)
  }
  dd <- unc %>% dplyr::filter(health_zone %in% invaded)
  if (!"province" %in% names(dd) && !is.null(province_map))
    dd <- dd %>% dplyr::left_join(
      province_map %>% dplyr::transmute(health_zone = nom, province), by = "health_zone")
  if (!"province" %in% names(dd)) dd$province <- NA_character_
  dd <- dd %>% dplyr::mutate(region = .prov_label(province))
  ev <- dd %>% dplyr::filter(is_new_invasion %in% 1)
  prov_pal <- .prov_colours(dd$region)

  p <- ggplot(dd, aes(cutoff, p_med)) +
    geom_ribbon(aes(ymin = p_lo, ymax = p_hi, fill = region), alpha = 0.25) +
    geom_line(aes(colour = region), linewidth = 0.7) +
    geom_point(aes(colour = region), size = 1.1) +
    { if (nrow(ev) > 0) geom_vline(data = ev, aes(xintercept = cutoff),
             linetype = "22", colour = "grey35", linewidth = 0.4) else NULL } +
    { if (nrow(ev) > 0) geom_point(data = ev, aes(cutoff, p_med),
             shape = 4, size = 2.6, stroke = 1, colour = "black") else NULL } +
    facet_wrap(~ health_zone, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = prov_pal, name = "Province") +
    scale_fill_manual(values = prov_pal, name = "Province") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = sprintf("Forecast invasion probability across space and time — Frequentist ensemble (h=%dw)", horizon),
         subtitle = "Median line + member-spread ribbon per zone | dashed line & cross = week the first case actually occurred",
         x = "Forecast date (fold cutoff)", y = "P(first case) with ensemble uncertainty",
         caption = "Zones that invaded during cross-validation; a band rising to the mark = an anticipated invasion.") +
    theme_inv(11) + theme(legend.position = "top",
                          strip.text = element_text(face = "bold", size = 9),
                          plot.caption = element_text(colour = "grey40", size = 8))
  if (save) .fd_save(p, file.path(OUT_MAPS, sprintf("spacetime_forecast_uncertainty_h%d.pdf", horizon)),
                     w = 11, h = 8)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Q1 — visualise the ENSEMBLE forecast on the map (vs the best renewal model), and a
#      companion map of WHERE the models disagree (spatial uncertainty).
# ---------------------------------------------------------------------------

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

#' Q1 — side-by-side Ituri maps: best renewal model vs ensemble (same colour
#' scale), plus a third panel of member disagreement (spatial uncertainty).
plot_forecast_map_panel <- function(fc, members, primary_method,
                                    ensemble_method = "Ensemble-mean",
                                    shapefile = NULL, horizon = 1L, save = TRUE) {
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile)) return(invisible(NULL))
  prim <- fc %>% dplyr::filter(method == primary_method, horizon == !!horizon)
  ens  <- fc %>% dplyr::filter(method == ensemble_method, horizon == !!horizon)
  if (nrow(prim) == 0 || nrow(ens) == 0) return(invisible(NULL))
  # Common colour scale across the two forecast panels for honest comparison.
  lim <- c(0, max(0.05, stats::quantile(c(prim$p_invasion, ens$p_invasion), 0.99, na.rm = TRUE)))
  unc <- ensemble_member_uncertainty(fc %>% dplyr::filter(horizon == !!horizon),
                                     members, "p_invasion")
  # Shared already-affected mask (from any model carrying was_active_before) so
  # all panels grey out the SAME zones regardless of which rows each retains.
  aff_keys <- fc %>% dplyr::filter(horizon == !!horizon, was_active_before %in% TRUE) %>%
    dplyr::pull(health_zone) %>% trimws() %>% tolower() %>% unique()

  p1 <- .fd_map(prim, shapefile, sprintf("Best renewal model\n%s", primary_method),
                ituri_only = TRUE, lim = lim, affected_keys = aff_keys)
  p2 <- .fd_map(ens, shapefile, sprintf("Ensemble\n%s", ensemble_method),
                ituri_only = TRUE, lim = lim, affected_keys = aff_keys)
  panels <- list(p1, p2)
  if (!is.null(unc)) {
    p3 <- .fd_map(unc, shapefile, "Member disagreement\n(min-max width)",
                  prob_col = "p_width", ituri_only = TRUE, palette = "viridis",
                  legend_name = "Uncertainty", affected_keys = aff_keys)
    panels <- c(panels, list(p3))
  }
  out <- patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_annotation(
      title = sprintf("Ituri invasion forecast: best model vs FREQUENTIST ensemble, with spatial uncertainty (h=%dw)", horizon),
      subtitle = sprintf("Best model: %s | ensemble: %s | grey = already affected (Bayesian analogues: bayes_ensemble_* maps + posterior-CrI figures)",
                         primary_method %||% "renewal", ensemble_method),
      theme = ggplot2::theme(plot.title = element_text(face = "bold", size = 13)))
  if (save) .fd_save(out, file.path(OUT_MAPS, sprintf("forecast_ensemble_panel_h%d.pdf", horizon)),
                     w = 15, h = 5.6)
  invisible(out)
}

# ---------------------------------------------------------------------------
# Q1 — key parameter estimates: the mobility import coefficient + a covariate
#      association screen (the full candidate range, each estimable)
# ---------------------------------------------------------------------------

.PARAM_LABEL <- c(
  alert_import = "Imported alert pressure (mobility-weighted)",
  alert_local  = "Local suspected alerts",
  susp_import  = "Imported suspected-case pressure (mobility-weighted)",
  susp_local   = "Local suspected (not-confirmed) cases",
  log_pop      = "Population size (log)",
  ccvi         = "Socioeconomic vulnerability (CCVI)",
  positivity   = "Test positivity",
  d_min        = "Travel time to nearest affected zone",
  healthsite_density = "Health-facility density")

#' Build the invasion design matrix (one row per at-risk zone-week transition):
#' `invaded`, the mobility import offset `logLam`, and the STANDARDISED candidate
#' covariates — exactly the rows fit_import_model uses. Returns the design plus
#' the baseline import coefficient beta0 (intercept-only cloglog).
build_invasion_design <- function(zone_week_nc, mobility_matrices, gt_pmfs,
                                  covariates, osrm_mat, zones_all,
                                  mob = "M8", gt = "medium",
                                  candidates = c("alert_import", "alert_local",
                                                 "susp_import", "susp_local", "log_pop",
                                                 "healthsite_density", "ccvi", "positivity",
                                                 "d_min")) {
  Y_wide <- .count_wide(zone_week_nc, zones_all, "confirmed_nc")
  G <- daily_to_weekly_gt(gt_pmfs[[gt]]); W <- mobility_matrices[[mob]]
  A_wide <- .count_wide(zone_week_nc, zones_all, "total_alerts")
  # Suspected-but-not-confirmed leading-indicator matrix (nowcast-corrected suspected
  # series, raw fallback), for the susp_import / susp_local covariates — like the alert matrix.
  S_wide <- .susp_wide(zone_week_nc, zones_all)
  static <- .static_features(covariates, zones_all)
  rows <- list()
  for (t in seq_len(ncol(Y_wide) - 1L)) {
    Lam <- compute_foi(Y_wide, W, G, t_idx = t + 1L, zones_all)
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

#' For each candidate covariate, its UNIVARIATE association with invasion, both
#' MARGINAL (raw) and ADJUSTED for the mobility import force (offset). The
#' marginal HR is always estimable; the adjusted fit can separate (HR -> 0/Inf)
#' when a covariate has no identifiable signal beyond mobility — that is itself the
#' finding, flagged rather than hidden.
covariate_associations <- function(design) {
  d <- design$d
  fit1 <- function(f, use_off) {
    fm <- stats::as.formula(paste("invaded ~", f))
    m <- tryCatch(suppressWarnings(stats::glm(fm, binomial("cloglog"),
           offset = if (use_off) d$logLam else NULL, data = d,
           method = if (requireNamespace("brglm2", quietly = TRUE)) brglm2::brglmFit else "glm.fit")),
         error = function(e) NULL)
    if (is.null(m)) return(c(hr = NA, lo = NA, hi = NA, p = NA))
    cf <- coef(m)[f]; se <- tryCatch(summary(m)$coefficients[f, 2], error = function(e) NA)
    pv <- tryCatch(summary(m)$coefficients[f, 4], error = function(e) NA)
    c(hr = exp(unname(cf)), lo = exp(unname(cf) - 1.96 * se),
      hi = exp(unname(cf) + 1.96 * se), p = unname(pv))
  }
  purrr::map_dfr(design$feat, function(f) {
    mg <- fit1(f, FALSE); aj <- fit1(f, TRUE)
    tibble::tibble(term = f, label = unname(.PARAM_LABEL[f]) %||% f,
      hr = mg["hr"], lo = mg["lo"], hi = mg["hi"], p = mg["p"],
      hr_adj = aj["hr"], p_adj = aj["p"],
      adj_status = dplyr::case_when(
        !is.finite(aj["hr"]) | aj["hr"] > 20 | aj["hr"] < 0.05 ~ "separates (no added signal)",
        aj["p"] < 0.05 ~ "still significant",
        TRUE ~ "attenuated (mediated by mobility)"))
  })
}

#' Q1 — covariate association screen: the marginal invasion hazard-ratio per +1 SD
#' for every candidate covariate (all estimable), with 95% CIs, annotated by what
#' happens after adjusting for the mobility import force (the dominant driver).
plot_model_parameters <- function(assoc, best_method, beta0 = NA, n_events = NA,
                                  save = TRUE) {
  if (is.null(assoc) || nrow(assoc) == 0) {
    message("[params] no covariate associations to plot"); return(invisible(NULL)) }
  d <- assoc %>% dplyr::filter(is.finite(hr)) %>% dplyr::arrange(hr) %>%
    dplyr::mutate(label = factor(label, levels = label),
                  dir = ifelse(hr >= 1, "associated with higher risk", "associated with lower risk"))
  p <- ggplot(d, aes(hr, label, colour = dir)) +
    geom_vline(xintercept = 1, linetype = "22", colour = "grey45") +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 3) +
    geom_text(aes(x = hi, label = paste0("  ", adj_status)), hjust = 0, size = 2.7,
              colour = "grey35", na.rm = TRUE) +
    scale_colour_manual(values = c(`associated with higher risk` = OKABE_ITO[2],
                                   `associated with lower risk` = OKABE_ITO[1]), name = NULL) +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.35))) +
    labs(title = sprintf("Invasion drivers — %s", best_method),
         subtitle = sprintf("Mobility import force is the dominant driver (offset, beta0=%.3g); bars = marginal HR per +1 SD, text = effect after adjusting for mobility (%s events)",
                            beta0, n_events),
         x = "Marginal invasion hazard ratio per +1 SD  (log scale; >1 = higher risk)", y = NULL,
         caption = stringr::str_wrap("Univariate complementary-log-log fits (Firth-penalised, 95% CI). The import force already carries most population/distance structure, so several covariates add little (or separate) once it is conditioned on.", width = 150)) +
    theme_inv(11) + theme(legend.position = "top",
                          plot.subtitle = element_text(size = 8.5, colour = "grey30"),
                          plot.caption = element_text(colour = "grey40", size = 8, hjust = 0),
                          plot.caption.position = "plot")
  if (save) .fd_save(p, file.path(OUT_REPORTS, "best_model_parameters.pdf"), w = 10.5, h = 5)
  invisible(p)
}

#' National effective reproduction number R(t) from the renewal estimate
#' (EpiNow2), with 60% and 90% credible bands and the three GT-profile means for
#' sensitivity. R(t) is a renewal-model output (it is what projects incidence one
#' week ahead for the 2-week invasion horizon); this plots it explicitly.
plot_rt <- function(rt_primary, rt_all = NULL, save = TRUE) {
  if (is.null(rt_primary) || nrow(rt_primary) == 0) return(invisible(NULL))
  d <- rt_primary %>% dplyr::mutate(date = as.Date(date))
  p <- ggplot(d, aes(date, R_mean)) +
    geom_hline(yintercept = 1, linetype = "22", colour = "grey40") +
    geom_ribbon(aes(ymin = R_lo_90, ymax = R_hi_90), fill = OKABE_ITO[1], alpha = 0.15) +
    geom_ribbon(aes(ymin = R_lo_60, ymax = R_hi_60), fill = OKABE_ITO[1], alpha = 0.30)
  if (!is.null(rt_all) && length(rt_all) > 1) {
    others <- dplyr::bind_rows(lapply(names(rt_all), function(nm) {
      x <- rt_all[[nm]]; if (is.null(x) || !nrow(x)) return(NULL)
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
         caption = "A renewal-model output on the national onset series; also used to project one week of incidence for the 2-week invasion horizon.") +
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
  pillars <- pillars[vapply(pillars, function(p) !all(is.na(d[[p]])), logical(1))]
  M <- as.matrix(d[, pillars, drop = FALSE])
  # Equal-weight mean over the available pillars; a zone with EVERY pillar NA gets V = NA (not
  # NaN, which rowMeans(na.rm=TRUE) would otherwise return and then propagate into priority).
  d$V <- ifelse(rowSums(is.finite(M)) == 0L, NA_real_, rowMeans(M, na.rm = TRUE))
  d
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
plot_priority_scatter <- function(rs, horizon = 1L, top_n = 12L, save = TRUE, window_txt = "",
                                  file = "priority_scatter", model_label = "") {
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
    labs(title = sprintf("Preparedness priority: invasion likelihood x vulnerability (h=%dw)%s",
                         horizon, if (nzchar(model_label)) paste0(" — ", model_label) else ""),
         subtitle = "Upper-right = likely to be invaded AND vulnerable/under-resourced = highest priority",
         x = "Vulnerability & capacity gap (0-1; higher = more vulnerable, less prepared)",
         y = "Relative invasion risk (0-1; 1 = highest)",
         caption = paste0("Vulnerability = mean percentile gap in surveillance (facility density), healthcare (facilities per capita), healthcare access (travel time to nearest facility) and CCVI deprivation. ", window_txt)) +
    theme_inv(11) + theme(legend.position = "right",
                          plot.caption = element_text(colour = "grey40", size = 8))
  if (save) {
    .pp <- file.path(OUT_REPORTS, sprintf("%s_h%d.pdf", file, horizon))
    .fd_save(p, .pp, w = 9.5, h = 6.5)
    # Also place the priority scatter in key_outputs/ (in ADDITION to reports) so it ships
    # with the headline deliverables — copied here at save time, not only via the end-of-run
    # key_outputs manifest, so it is present even on a partial run.
    .ko <- file.path(get0("OUT_DIR", ifnotfound = dirname(OUT_REPORTS)), "key_outputs")
    dir.create(.ko, showWarnings = FALSE, recursive = TRUE)
    if (file.exists(.pp))
      file.copy(.pp, file.path(.ko, basename(.pp)), overwrite = TRUE)
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
    "weight mean of the available percentile-ranked pillars (four when an OSRM travel-time",
    "matrix is supplied, as in the default run; the healthcare-access axis is dropped otherwise):", "",
    "1. **Surveillance gap** = 1 - rank(health-facility density) — detection reach",
    "   (dedicated PCR-testing data covers only ~4% of zones, so facility density is",
    "   the robust per-zone surveillance proxy);",
    "2. **Healthcare gap** = 1 - rank(health facilities per capita) — treatment load;",
    "3. **Healthcare access gap** = rank(road travel time to the nearest zone with a facility;",
    "   own-facility zones = 0) — physical reach of care, distinct from facility density/per-capita;",
    "4. **Social vulnerability** = rank(CCVI socioeconomic deprivation).", "",
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
              M13 = "cohort+gravity", M14 = "cohort+radiation")[mob]
    # OSRM road-DISTANCE deterrence variant (vs the default travel-time deterrence).
    if (grepl("-dist", cc) && !is.na(mobn)) mobn <- paste0(mobn, ", road-distance")
    gt   <- if (grepl("-short", cc)) "short GT" else if (grepl("-long", cc)) "long GT" else "med GT"
    cov  <- if (grepl("-full", cc)) " · +full covariates" else
            if (grepl("-geo", cc)) " · +geo covariates" else ""
    lnk  <- if (grepl("-logit", cc)) " · logit link" else if (grepl("-probit", cc)) " · probit link" else ""
    sprintf("%s  —  %s · %s%s%s", cc, if (is.na(mobn)) mob else mobn, gt, cov, lnk)
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
  lab <- function(t) unname(.PARAM_LABEL[t]) %||% t
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
  lab <- function(t) unname(.PARAM_LABEL[t]) %||% t
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
compute_params_over_time <- function(zone_week_outbreak, cutoffs, mobility_matrices,
                                     gt_pmfs, covariates, osrm_mat, zones_all,
                                     mob = "M8", gt = "medium",
                                     cov_spec = c("log_pop", "ccvi", "d_min"),
                                     nowcast_fn = NULL) {
  if (is.null(nowcast_fn)) nowcast_fn <- get0("apply_nowcast_correction")
  purrr::map_dfr(cutoffs, function(cut) {
    zw <- zone_week_outbreak %>% dplyr::filter(week_start <= cut)
    zwn <- tryCatch(nowcast_fn(zw, analysis_date = cut + 7), error = function(e) zw)
    des <- tryCatch(build_invasion_design(zwn, mobility_matrices, gt_pmfs, covariates,
             osrm_mat, zones_all, mob = mob, gt = gt, candidates = cov_spec),
             error = function(e) NULL)
    if (is.null(des) || des$n_events < 2L) return(NULL)
    m <- tryCatch(suppressWarnings(stats::glm(
      stats::as.formula(paste("invaded ~", paste(des$feat, collapse = " + "))),
      family = binomial("cloglog"), offset = des$d$logLam, data = des$d,
      method = if (requireNamespace("brglm2", quietly = TRUE)) brglm2::brglmFit else "glm.fit")),
      error = function(e) NULL)
    if (is.null(m)) return(NULL)
    co <- summary(m)$coefficients
    terms <- intersect(des$feat, rownames(co))
    tibble::tibble(cutoff = as.Date(cut), n_events = des$n_events, term = terms,
      hr = exp(co[terms, 1]),
      lo = exp(co[terms, 1] - 1.96 * co[terms, 2]),
      hi = exp(co[terms, 1] + 1.96 * co[terms, 2]))
  })
}

#' Parameter estimates over time (fold cutoffs): HR + 95% CI per driver.
plot_params_over_time <- function(pot, save = TRUE, file = "params_over_time",
                                  model_label = "") {
  if (is.null(pot) || !nrow(pot)) return(invisible(NULL))
  pot <- pot %>% dplyr::mutate(label = vapply(term, function(t) unname(.PARAM_LABEL[t]) %||% t, character(1)))
  p <- ggplot(pot, aes(cutoff, hr)) +
    geom_hline(yintercept = 1, linetype = "22", colour = "grey45") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = OKABE_ITO[1], alpha = 0.18) +
    geom_line(colour = OKABE_ITO[1], linewidth = 0.7) + geom_point(size = 1.2) +
    facet_wrap(~ label, scales = "free_y") +
    scale_y_log10() +
    labs(title = sprintf("Invasion-driver estimates over time (refit at each fold cutoff)%s",
                         if (nzchar(model_label)) paste0(" — ", model_label) else ""),
         subtitle = "Hazard ratio per +1 SD with 95% CI / 90% CrI; x = forecast date (training data grows left to right)",
         x = "Fold cutoff (forecast date)", y = "Hazard ratio per +1 SD (log scale)",
         caption = "Identifiable exogenous drivers only (log_pop, CCVI, d_min); shows whether/when each effect stabilises. Note: the design is re-standardised per fold, so the +1 SD unit varies slightly across cutoffs.") +
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
         x = "Fold cutoff (forecast date)", y = "Probability / fraction",
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

#' Model-ranking comparison: AUC-PR skill per method coloured by inferential family,
#' so the Bayesian suite is ranked ON THE SAME leave-future-out folds as the
#' frequentist models and the two paradigms are directly comparable.
plot_freq_bayes_ranking <- function(eval_tbl, horizon = 1L, save = TRUE) {
  if (is.null(eval_tbl) || !nrow(eval_tbl) || !"auc_pr_skill" %in% names(eval_tbl))
    return(invisible(NULL))
  d <- eval_tbl %>% dplyr::filter(horizon == !!horizon, is.finite(auc_pr_skill)) %>%
    dplyr::mutate(family = .model_family(method),
                  method = stats::reorder(method, auc_pr_skill))
  if (!nrow(d)) return(invisible(NULL))
  p <- ggplot(d, aes(auc_pr_skill, method, colour = family)) +
    geom_vline(xintercept = 1, linetype = "22", colour = "grey55") +
    geom_segment(aes(x = 1, xend = auc_pr_skill, yend = method), linewidth = 0.5) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = c(Frequentist = OKABE_ITO[1], Bayesian = OKABE_ITO[2]), name = NULL) +
    labs(title = sprintf("Frequentist vs Bayesian model skill (h=%dw)", horizon),
         subtitle = "AUC-PR skill (Average Precision / base rate) on the same LFO folds; dashed = no skill (1)",
         x = "AUC-PR skill (x base rate)", y = NULL,
         caption = "Both paradigms fit the SAME cloglog renewal invasion model; only the estimation (penalised MLE vs posterior) differs.") +
    theme_inv(11) + theme(legend.position = "top",
                          plot.caption = element_text(colour = "grey40", size = 8, hjust = 0))
  if (save) .fd_save(p, file.path(OUT_DIAGNOSTICS, sprintf("freq_vs_bayes_ranking_h%d.pdf", horizon)), w = 9, h = 6)
  invisible(p)
}

#' Per-zone agreement between the frequentist featured forecast and the Bayesian
#' stacked posterior — a robustness check that the two paradigms flag the same
#' at-risk zones. Points on the diagonal = agreement.
plot_freq_bayes_agreement <- function(freq_rs, bayes_stack, province_map = NULL,
                                      horizon = 1L, save = TRUE) {
  if (is.null(freq_rs) || is.null(bayes_stack)) return(invisible(NULL))
  fcol <- if ("p_case_invasion" %in% names(freq_rs)) "p_case_invasion" else "p_invasion"
  a <- freq_rs %>% dplyr::filter(horizon == !!horizon, !(was_active_before %in% TRUE)) %>%
    dplyr::transmute(health_zone, p_freq = .data[[fcol]])
  b <- bayes_stack %>% dplyr::filter(horizon == !!horizon, !(was_active_before %in% TRUE)) %>%
    dplyr::transmute(health_zone, p_bayes = p_invasion)
  d <- dplyr::inner_join(a, b, by = "health_zone") %>%
    dplyr::filter(is.finite(p_freq), is.finite(p_bayes))
  if (!nrow(d)) return(invisible(NULL))
  if (!is.null(province_map))
    d <- d %>% dplyr::left_join(province_map %>% dplyr::transmute(health_zone = nom, province),
                                by = "health_zone")
  if (!"province" %in% names(d)) d$province <- NA_character_
  d$region <- .prov_label(d$province)
  rho <- suppressWarnings(stats::cor(d$p_freq, d$p_bayes, method = "spearman"))
  # Shared axis maximum + coord_equal so the dashed y=x line is a TRUE 45-degree reference
  # (otherwise independent x/y ranges make points look off-diagonal when they in fact agree).
  .axmax <- max(d$p_freq, d$p_bayes, 0.02, na.rm = TRUE) * 1.05
  p <- ggplot(d, aes(p_freq, p_bayes, colour = region)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_point(size = 2, alpha = 0.78) +
    scale_colour_manual(values = .prov_colours(d$region), name = "Province") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, .axmax)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, .axmax)) +
    coord_equal(xlim = c(0, .axmax), ylim = c(0, .axmax)) +
    labs(title = sprintf("Frequentist vs Bayesian invasion probability (h=%dw)", horizon),
         subtitle = sprintf("At-risk zones; Spearman rho = %.2f. On the diagonal = the paradigms agree.",
                            ifelse(is.finite(rho), rho, NA_real_)),
         x = "Frequentist featured P(first case)", y = "Bayesian stacked P(first case)") +
    theme_inv(11) + theme(legend.position = "right")
  if (save) .fd_save(p, file.path(OUT_DIAGNOSTICS, sprintf("freq_vs_bayes_agreement_h%d.pdf", horizon)), w = 8.5, h = 6.5)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Task 7 — per-province masked invasion maps (Ituri, Nord-Kivu, Haut-Uele, ...)
# ---------------------------------------------------------------------------

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
plot_prob_vuln_choropleth <- function(rs, shapefile = NULL, horizon = 1L,
                                      province_zoom = "Ituri", window_txt = "",
                                      save = TRUE, file = NULL, model_label = "") {
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
    labs(title = sprintf("%s — invasion risk x vulnerability (h=%dw)%s",
                         province_zoom %||% "National", horizon,
                         if (nzchar(model_label)) paste0(" — ", model_label) else ""),
         subtitle = "Darker blue = both likelier to be invaded AND more vulnerable/under-resourced (national quartiles, 4x4)",
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
compute_detection_curve <- function(lfo_results, method, horizon = 1L, ks = 1:25) {
  d <- lfo_results %>% dplyr::filter(method == !!method, horizon == !!horizon,
                                     is.finite(p_invasion))
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
    tp <- sum(per_fold$is_new_invasion[per_fold$rk <= k], na.rm = TRUE)
    tibble::tibble(method = method, horizon = horizon, k = k,
                   recall = tp / tot_pos, precision = tp / (k * nfold),
                   caught = tp / nfold, n_events = tot_pos / nfold,
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
naive_epicentre_inflow_scores <- function(W, epicentre_zones, pop_vec = NULL,
                                          zones_all = rownames(W)) {
  if (is.null(W) || is.null(zones_all)) return(NULL)
  ez <- intersect(as.character(epicentre_zones), rownames(W))
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
  ref <- lfo_results %>% dplyr::filter(method == unique(lfo_results$method)[1])
  if (!nrow(ref) || !"health_zone" %in% names(ref)) return(lfo_results)
  naive <- ref
  naive$method     <- label
  naive$p_invasion <- as.numeric(scores[match(naive$health_zone, names(scores))])
  if ("p_case_invasion" %in% names(naive)) naive$p_case_invasion <- naive$p_invasion
  for (col in intersect(c("mu_forecast", "p_infection_invasion", "p_lo", "p_hi", "p_sd"),
                        names(naive)))
    naive[[col]] <- NA_real_
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
plot_detection_curve <- function(lfo_results, methods, horizon = 1L, save = TRUE,
                                 annotate_k = 10L) {
  curves <- purrr::map_dfr(methods, function(m) compute_detection_curve(lfo_results, m, horizon))
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
         subtitle = "Share of true next-week invasions caught at a fixed weekly alert budget of K zones, pooled over LFO folds",
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
plot_topk_precision <- function(lfo_results, methods, horizon = 1L, save = TRUE,
                                file = "topk_precision", title_suffix = "") {
  curves <- purrr::map_dfr(methods, function(m) compute_detection_curve(lfo_results, m, horizon))
  if (is.null(curves) || !nrow(curves)) return(invisible(NULL))
  mth <- unique(curves$method)
  pal <- if (length(mth) <= length(OKABE_ITO)) OKABE_ITO[seq_along(mth)]
         else grDevices::hcl.colors(length(mth), "Dark 3")
  # base rate = invasions per at-risk zone-week = a random watch-list's hit rate.
  # recall_random = k / n_atrisk => n_atrisk = k/recall_random, so base = n_events/n_atrisk.
  base_rate <- suppressWarnings(mean(
    (curves$n_events * curves$recall_random / curves$k)[curves$k > 0 & is.finite(curves$recall_random)],
    na.rm = TRUE))
  p <- ggplot(curves, aes(k, precision, colour = method)) +
    { if (is.finite(base_rate)) geom_hline(yintercept = base_rate, linetype = "22", colour = "grey55") } +
    geom_line(linewidth = 0.9) + geom_point(size = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
    scale_colour_manual(values = setNames(pal, mth), name = NULL) +
    labs(title = sprintf("Of the top-K highest-risk zones, what %% were actually invaded? (h=%dw)%s",
                         horizon, if (nzchar(title_suffix)) paste0(" — ", title_suffix) else ""),
         subtitle = "Precision at a fixed weekly alert budget of K zones, pooled over LFO folds; dashed = random-watch-list base rate",
         x = "Zones monitored each week (K)", y = "% of top-K zones that were invaded",
         caption = "Higher = fewer false alarms. Precision declines with K as the budget reaches lower-risk zones; the gap above the base rate is the model's targeting value.") +
    theme_inv(11) + theme(legend.position = "top",
                          plot.caption = element_text(colour = "grey30", size = 8.5, hjust = 0, face = "italic"))
  if (save) .fd_save(p, file.path(OUT_DIAGNOSTICS, sprintf("%s_h%d.pdf", file, horizon)), w = 9, h = 5.6)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Task 4 — reporting rate (across space & folds) and import coefficient over folds
# ---------------------------------------------------------------------------

#' Spatial map of the per-zone RELATIVE reporting-rate proxy the model uses to
#' up-weight under-ascertained source zones (health-site density, geometric-mean
#' normalised, bounded [0.25, 4]). r < 1 = likely under-reporting (its true import
#' pressure exceeds observed counts); r > 1 = better-than-average ascertainment.
plot_reporting_rate_map <- function(covariates, shapefile = NULL, zones_all = NULL,
                                    save = TRUE) {
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile) || is.null(covariates) || !"healthsite_density" %in% names(covariates))
    return(invisible(NULL))
  # Derive the reporting-rate proxy on the covariate zones directly (the shapefile
  # join maps them by name); using match(zones_all, covariates$nom) would NA-out the
  # whole vector if the population-spine names and covariate names differ.
  za   <- as.character(covariates$nom)
  dens <- suppressWarnings(as.numeric(covariates$healthsite_density))
  if (!is.null(zones_all)) { keep <- za %in% zones_all; if (sum(keep) >= 5) { za <- za[keep]; dens <- dens[keep] } }
  rr <- get0(".reporting_rate_vec")
  r  <- if (is.function(rr)) rr(stats::setNames(dens, za), za) else NULL
  if (is.null(r)) { message("[reporting_rate_map] no valid reporting-rate proxy; skipped."); return(invisible(NULL)) }
  rd <- tibble::tibble(.key = tolower(trimws(za)), r = as.numeric(r))
  m <- shapefile %>% dplyr::mutate(.key = tolower(trimws(Nom))) %>%
    dplyr::left_join(rd, by = ".key")
  mk <- function(dat, ttl) ggplot(dat) +
    geom_sf(aes(fill = r), colour = "grey80", linewidth = 0.12) +
    scale_fill_gradient2(midpoint = 1, low = "#b2182b", mid = "grey92", high = "#2166ac",
                         name = "Relative\nreporting rate", trans = "log10",
                         na.value = "grey95") +
    labs(title = ttl) + ggplot2::theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))
  out <- mk(dplyr::filter(m, PROVINCE == "Ituri"), "Ituri — relative reporting-rate proxy") +
    mk(m, "National") + patchwork::plot_layout(widths = c(2, 1)) +
    patchwork::plot_annotation(
      caption = "Ascertainment proxy = health-site density (relative, bounded). Under-reporting source zones (red) exert MORE true import pressure than their observed counts; this is exactly the reporting-rate structure the model applies to the mobility rows.",
      theme = ggplot2::theme(plot.caption = element_text(colour = "grey40", size = 8, hjust = 0)))
  if (save) .fd_save(out, file.path(OUT_MAPS, "reporting_rate_map.pdf"), w = 12, h = 6.5)
  invisible(out)
}

#' Refit the import coefficient beta0 = exp(intercept) of the cloglog renewal model
#' at each fold cutoff (causally, training only on week <= cutoff), and the mean
#' nowcast completeness applied that fold — so one can see whether the calibrated
#' import->invasion conversion and the reporting correction are stable over time.
compute_beta_over_folds <- function(zone_week_outbreak, cutoffs, mobility_matrices,
                                    gt_pmfs, covariates, osrm_mat, zones_all,
                                    mob = "M8", gt = "medium", nowcast_fn = NULL) {
  if (is.null(nowcast_fn)) nowcast_fn <- get0("apply_nowcast_correction")
  purrr::map_dfr(cutoffs, function(cut) {
    zw  <- zone_week_outbreak %>% dplyr::filter(week_start <= cut)
    zwn <- tryCatch(nowcast_fn(zw, analysis_date = cut + 7), error = function(e) zw)
    des <- tryCatch(build_invasion_design(zwn, mobility_matrices, gt_pmfs, covariates,
             osrm_mat, zones_all, mob = mob, gt = gt, candidates = character(0)),
             error = function(e) NULL)
    if (is.null(des) || des$n_events < 1L) return(NULL)
    m <- tryCatch(suppressWarnings(stats::glm(invaded ~ 1, binomial("cloglog"),
           offset = des$d$logLam, data = des$d)), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    co <- summary(m)$coefficients
    b0 <- co[1, 1]; se <- co[1, 2]
    comp <- tryCatch({
      recent <- zwn %>% dplyr::filter(week_start > cut - 21)
      denom <- sum(recent$confirmed_nc, na.rm = TRUE)   # AGGREGATE completeness, not a
      num   <- sum(recent$confirmed,    na.rm = TRUE)   # mean of per-zone ratios (which
      if (denom > 0) min(num / denom, 1) else NA_real_  # ~500 zero zones drag to ~0)
    }, error = function(e) NA_real_)
    tibble::tibble(cutoff = as.Date(cut), n_events = des$n_events,
      beta0 = exp(b0), beta0_lo = exp(b0 - 1.96 * se), beta0_hi = exp(b0 + 1.96 * se),
      mean_completeness = comp)
  })
}

#' Two-panel over-folds figure: import coefficient beta0 (with 95% CI) and mean
#' reporting completeness, against forecast date.
plot_beta_over_folds <- function(bof, save = TRUE, file = "beta_and_completeness_over_folds",
                                 model_label = "", ci_label = "95% CI") {
  if (is.null(bof) || !nrow(bof)) return(invisible(NULL))
  .ml <- if (nzchar(model_label)) paste0(" — ", model_label) else ""
  p1 <- ggplot(bof, aes(cutoff, beta0)) +
    geom_ribbon(aes(ymin = beta0_lo, ymax = beta0_hi), fill = OKABE_ITO[1], alpha = 0.18) +
    geom_line(colour = OKABE_ITO[1], linewidth = 0.7) + geom_point(size = 1.4) +
    scale_y_log10() +
    labs(title = paste0("Import coefficient beta0 over folds", .ml),
         subtitle = sprintf("exp(intercept) of the renewal fit, refit at each cutoff (%s)", ci_label),
         x = "Fold cutoff (forecast date)", y = "beta0 = P-scale import->invasion hazard") +
    theme_inv(11)
  # completeness panel only when it is populated (frequentist path); the Bayesian
  # beta0-over-folds has no completeness column, so show beta0 alone.
  out <- if ("mean_completeness" %in% names(bof) && any(is.finite(bof$mean_completeness))) {
    p2 <- ggplot(bof, aes(cutoff, mean_completeness)) +
      geom_line(colour = OKABE_ITO[8], linewidth = 0.7) + geom_point(size = 1.4) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(title = "Mean reporting completeness over folds",
           subtitle = "Observed / nowcast-corrected counts in the recent weeks at each cutoff",
           x = "Fold cutoff (forecast date)", y = "Reporting completeness") +
      theme_inv(11)
    p1 / p2
  } else p1
  if (save) .fd_save(out, file.path(OUT_DIAGNOSTICS, sprintf("%s.pdf", file)), w = 9,
                     h = if (inherits(out, "patchwork")) 8 else 4.5)
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

#' Animate a plot across horizons into a GIF (#4). `render_fn(hz)` returns a
#' ggplot/patchwork for horizon hz; frames are written with gifski (a true GIF) or,
#' if gifski is unavailable, a multi-page "flipbook" PDF. Used for the key Bayesian
#' outputs so the 1- vs 2-week picture can be viewed as a short loop.
#' Render a list of ggplot/patchwork frames to a GIF (gifski) or, if gifski is
#' unavailable, a multi-page "flipbook" PDF. Shared encoder for all animators.
.frames_to_gif <- function(plots, out_stem, w = 11, h = 6.5, dpi = 110, delay = 1.2) {
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) < 1) return(invisible(NULL))
  if (requireNamespace("gifski", quietly = TRUE)) {
    dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
    pngs <- vapply(seq_along(plots), function(i) {
      f <- tempfile(fileext = ".png")
      suppressMessages(ggplot2::ggsave(f, plots[[i]], width = w, height = h, dpi = dpi, device = dev))
      f }, character(1))
    gif <- paste0(out_stem, ".gif")
    ok <- tryCatch({ suppressWarnings(gifski::gifski(pngs, gif, width = round(w * dpi),
            height = round(h * dpi), delay = delay, progress = FALSE)); TRUE }, error = function(e) FALSE)
    unlink(pngs)
    if (ok) { message(sprintf("[anim] wrote %s (%d frames)", basename(gif), length(plots))); return(invisible(gif)) }
  }
  pdf_out <- paste0(out_stem, "_flipbook.pdf")
  grDevices::pdf(pdf_out, width = w, height = h, onefile = TRUE)
  for (p in plots) print(p); grDevices::dev.off()
  message(sprintf("[anim] gifski unavailable; wrote flipbook %s", basename(pdf_out)))
  invisible(pdf_out)
}

make_horizon_animation <- function(render_fn, horizons, out_stem, w = 11, h = 6.5,
                                   dpi = 110, delay = 1.2) {
  plots <- lapply(horizons, function(hz) tryCatch(render_fn(hz), error = function(e) NULL))
  if (length(Filter(Negate(is.null), plots)) < length(horizons))
    warning("[anim] dropped horizon frame(s) for ", basename(out_stem))
  .frames_to_gif(plots, out_stem, w, h, dpi, delay)
}

#' TIME-EVOLUTION animation (#1/#2/#3): how the invasion-probability map evolves across the
#' leave-future-out forecast ORIGINS (successive weeks), for ONE model and ONE horizon, with
#' a FIXED colour scale across every frame (#2 — frames are directly comparable). Produced
#' separately per horizon (#1) and for the Ituri zoom or the whole DRC (extent="national", #3).
animate_invasion_over_time <- function(lfo_results, method, horizon = 1L, shapefile = NULL,
                                       extent = c("ituri", "national"), out_stem = NULL,
                                       palette = "plasma", model_label = NULL) {
  extent <- match.arg(extent)
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile) || is.null(lfo_results) || is.null(method)) return(invisible(NULL))
  pcol <- if ("p_case_invasion" %in% names(lfo_results)) "p_case_invasion" else "p_invasion"
  d <- lfo_results %>% dplyr::filter(method == !!method, horizon == !!horizon,
                                     is.finite(.data[[pcol]]))
  if (!nrow(d) || !"cutoff" %in% names(d)) return(invisible(NULL))
  folds <- sort(unique(as.Date(d$cutoff)))
  if (length(folds) < 2) return(invisible(NULL))
  # FIXED legend across frames (#2): a single shared colour limit over ALL origins.
  lim <- c(0, max(0.05, stats::quantile(d[[pcol]], 0.99, na.rm = TRUE)))
  reg <- if (extent == "ituri") "Ituri" else "DRC (national)"
  lab <- if (!is.null(model_label)) paste0(" — ", model_label) else ""
  frames <- lapply(folds, function(ct) {
    di <- d %>% dplyr::filter(as.Date(cutoff) == ct)
    aff <- di %>% dplyr::filter(was_active_before %in% TRUE) %>%
      dplyr::pull(health_zone) %>% trimws() %>% tolower() %>% unique()
    .fd_map(di, shapefile,
            sprintf("%s — P(first case) in next %dw, forecast origin %s%s",
                    reg, horizon, format(ct, "%d %b %Y"), lab),
            prob_col = pcol, ituri_only = (extent == "ituri"), lim = lim,
            palette = palette, legend_name = "P(first case)", affected_keys = aff)
  })
  if (is.null(out_stem))
    out_stem <- file.path(OUT_MAPS, "animations",
                          sprintf("invasion_over_time_%s_h%d", extent, horizon))
  dir.create(dirname(out_stem), showWarnings = FALSE, recursive = TRUE)
  .frames_to_gif(frames, out_stem, w = if (extent == "national") 9 else 8, h = 7,
                 dpi = 110, delay = 1.1)
}

#' TIME-EVOLUTION of the invasion-probability RANKING map across the LFO forecast origins, for
#' ONE model and horizon, at the given extent (national by default). Each frame ranks the at-risk
#' zones by P(first case) that week (1 = highest) and draws the rank map with a FIXED scale, so
#' the shifting priority ordering is visible as a GIF. Reuses plot_invasion_rank_map per frame.
animate_invasion_rank_over_time <- function(lfo_results, method, horizon = 1L, shapefile = NULL,
                                            extent = c("national", "ituri"), out_stem = NULL,
                                            model_label = NULL, top_k = 12L) {
  extent <- match.arg(extent)
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile) || is.null(lfo_results) || is.null(method)) return(invisible(NULL))
  pcol <- if ("p_case_invasion" %in% names(lfo_results)) "p_case_invasion" else "p_invasion"
  d <- lfo_results %>% dplyr::filter(method == !!method, horizon == !!horizon,
                                     is.finite(.data[[pcol]]))
  if (!nrow(d) || !"cutoff" %in% names(d)) return(invisible(NULL))
  folds <- sort(unique(as.Date(d$cutoff)))
  if (length(folds) < 2) return(invisible(NULL))
  reg <- if (extent == "national") "DRC (national)" else "Ituri"
  frames <- lapply(folds, function(ct) {
    di <- d %>% dplyr::filter(as.Date(cutoff) == ct)
    ml <- sprintf("%s · origin %s", model_label %||% method, format(ct, "%d %b %Y"))
    # rank map has a fixed rank scale (limits 1..cap, cap~=60 every fold) => frames comparable
    tryCatch(plot_invasion_rank_map(di, horizon = horizon, method_label = ml, shapefile = shapefile,
             top_k = top_k, save = FALSE, extent = extent), error = function(e) NULL)
  })
  if (is.null(out_stem))
    out_stem <- file.path(OUT_MAPS, "animations",
                          sprintf("invasion_rank_over_time_%s_h%d", extent, horizon))
  dir.create(dirname(out_stem), showWarnings = FALSE, recursive = TRUE)
  .frames_to_gif(frames, out_stem, w = if (extent == "national") 9 else 8, h = 7, dpi = 110, delay = 1.1)
}

#' How the ASSUMED mobility patterns route infection OVER TIME. The mobility matrix W is
#' time-invariant, but the force it routes — Lambda_i(t) = sum_j W[j,i] * sum_k g(k) Y_j(t-k)
#' — evolves as the outbreak spreads. Each weekly frame maps the import-force FIELD (fill,
#' fixed log scale across frames) plus the dominant import FLOWS (arrows: dominant source ->
#' destination) for the top at-risk zones, so the shifting mobility pathways are visible.
#' Rendered as a GIF (one per extent) and a static small-multiples panel of the last frames.
plot_mobility_flows_over_time <- function(zone_week_nc, mobility_matrices, gt_pmfs, zones_all,
                                          shapefile = NULL,
                                          mob = get0("MOBILITY_PRIMARY", ifnotfound = "M8"),
                                          gt = "medium", extent = c("ituri", "national"),
                                          top_k = 12L, out_stem = NULL, save = TRUE) {
  extent <- match.arg(extent)
  if (is.null(shapefile)) shapefile <- tryCatch(sf::st_read(SHAPEFILE_PATH, quiet = TRUE),
                                                error = function(e) NULL)
  if (is.null(shapefile) || is.null(mobility_matrices[[mob]])) return(invisible(NULL))
  W  <- mobility_matrices[[mob]]
  G  <- daily_to_weekly_gt(gt_pmfs[[gt]])
  Yw <- .count_wide(zone_week_nc, zones_all, "confirmed_nc")
  weeks <- as.Date(colnames(Yw))
  zc <- intersect(zones_all, rownames(W))
  if (length(zc) < 2 || length(weeks) < 3) return(invisible(NULL))
  key <- function(z) tolower(trimws(z))
  # zone centroids in the shapefile CRS (coord_sf overlays geom_segment in the same coords)
  shp <- shapefile %>% dplyr::mutate(.key = key(Nom))
  if (extent == "ituri") shp <- shp %>% dplyr::filter(PROVINCE == "Ituri")
  # Robust zone anchor points: some national polygons have invalid rings (duplicate vertices),
  # which makes s2 st_centroid error. Disable s2 (restored on exit), make the geometry valid,
  # and use point-on-surface (guaranteed inside the polygon, better arrow anchor than a centroid
  # that can fall outside a concave zone).
  .s2_old <- suppressMessages(sf::sf_use_s2()); on.exit(suppressMessages(sf::sf_use_s2(.s2_old)), add = TRUE)
  suppressMessages(try(sf::sf_use_s2(FALSE), silent = TRUE))
  gshp <- tryCatch(sf::st_make_valid(sf::st_geometry(shp)), error = function(e) sf::st_geometry(shp))
  ctr <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(gshp)))
  cen <- tibble::tibble(.key = shp$.key, x = ctr[, 1], y = ctr[, 2])
  # g-weighted source incidence at week index t
  .yweight <- function(t) {
    yw <- numeric(length(zc)); names(yw) <- zc
    for (k in seq_along(G)) { tp <- t - k; if (tp < 1) break
      yw <- yw + G[k] * Yw[zc, tp] }
    yw
  }
  fr_idx <- which(seq_along(weeks) >= 3L)                       # need >=2 weeks of history
  if (!length(fr_idx)) return(invisible(NULL))
  lam_max <- max(vapply(fr_idx, function(t) {
    yw <- .yweight(t); max(as.numeric(t(W[zc, zc]) %*% yw), 0) }, numeric(1)), 0.01)
  Wc <- W[zc, zc]
  frame_at <- function(t) {
    yw   <- .yweight(t)
    flow <- sweep(Wc, 1, yw, "*")                # flow[j,i] = W[j,i] * yweight[j]
    Lam  <- colSums(flow)                        # dest import force = sum_j flow[j,i]
    cum  <- rowSums(Yw[zc, seq_len(t), drop = FALSE])
    atrisk <- names(Lam)[cum == 0 & is.finite(Lam) & Lam > 0]
    dest <- if (length(atrisk)) atrisk[order(Lam[atrisk], decreasing = TRUE)][seq_len(min(top_k, length(atrisk)))] else character(0)
    seg <- NULL
    if (length(dest)) seg <- do.call(rbind, lapply(dest, function(i) {
      col <- flow[, i]; if (all(col <= 0)) return(NULL)
      j <- names(col)[which.max(col)]
      a <- cen[cen$.key == key(j), ]; b <- cen[cen$.key == key(i), ]
      if (!nrow(a) || !nrow(b)) return(NULL)
      data.frame(x = a$x[1], y = a$y[1], xend = b$x[1], yend = b$y[1], w = col[[j]])
    }))
    m <- dplyr::left_join(shp, tibble::tibble(.key = key(zc), lam = Lam), by = ".key")
    g <- ggplot(m) +
      geom_sf(aes(fill = lam), colour = "grey85", linewidth = 0.1) +
      scale_fill_viridis_c(option = "mako", trans = "log1p", name = "Import force Λ",
                           limits = c(0, lam_max), na.value = "grey92") +
      labs(title = sprintf("%s — mobility-routed import force & dominant flows (%s), week of %s",
                           if (extent == "ituri") "Ituri" else "DRC (national)", mob,
                           format(weeks[t], "%d %b %Y")),
           subtitle = "Static mobility matrix; arrows = dominant source->destination import pathway for the top at-risk zones") +
      ggplot2::theme_void(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 11), legend.position = "right")
    if (!is.null(seg) && nrow(seg))
      g <- g + geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend, linewidth = w),
                            arrow = grid::arrow(length = grid::unit(5, "pt"), type = "closed"),
                            colour = OKABE_ITO[2], alpha = 0.75, lineend = "round") +
        ggplot2::scale_linewidth(range = c(0.2, 1.6), guide = "none")
    g
  }
  frames <- lapply(fr_idx, function(t) tryCatch(frame_at(t), error = function(e) NULL))
  if (is.null(out_stem))
    out_stem <- file.path(OUT_MAPS, "animations", sprintf("mobility_flows_over_time_%s", extent))
  dir.create(dirname(out_stem), showWarnings = FALSE, recursive = TRUE)
  if (save) {
    .frames_to_gif(frames, out_stem, w = if (extent == "national") 9 else 8, h = 7, dpi = 110, delay = 1.1)
    # static small-multiples of the last (up to) 6 frames for the report
    fr <- Filter(Negate(is.null), frames); fr <- tail(fr, 6)
    if (length(fr)) {
      panel <- patchwork::wrap_plots(lapply(fr, function(g) g + theme(legend.position = "none")),
                                     ncol = 3) +
        patchwork::plot_annotation(
          title = sprintf("Mobility-routed import force over time (%s, %s)",
                          if (extent == "ituri") "Ituri" else "national", mob),
          theme = ggplot2::theme(plot.title = element_text(face = "bold", size = 13)))
      .fd_save(panel, file.path(OUT_MAPS, sprintf("mobility_flows_over_time_%s.pdf", extent)),
               w = if (extent == "national") 15 else 13, h = 8.5)
    }
  }
  invisible(frames)
}

message("[forecast_detail] 20_forecast_detail.R loaded (spec/selection + ensemble/uncertainty + params/priority + bayes/over-time + freq-vs-bayes + province/bivariate + detection + topk + reporting/beta + uncertainty-map + horizon & over-time animation viz).")
