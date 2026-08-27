# ============================================================================
# 04c_delay_four_windows.R  — MVE linelist delay analysis (v7)
#
# Comprehensive delay distribution analysis using a single analytical window:
#   Onset dates in [ANALYSIS_START, TRUNC_DATE]
#   where ANALYSIS_START = 2026-05-05 (protocol-mandated fixed start)
#   and   TRUNC_DATE     = max(sample_date) − 5 days (derived from data;
#                          matches TRAIN_CUTOFF in 15/16)
#
# Both dates are computed dynamically — no hard-coding of cutoff dates.
#
# Methods:
#   A. Naive uncensored MLE      — fitdist(); zeros excluded; lnorm uses x+0.5
#   B. Interval-censored MLE     — fitdistcens(); d=0→[0,0.5]; d>0→[d-0.5,d+0.5]
#   C. Bayesian EpiDist          — naive + marginal × lognormal + gamma (optional)
#
# Primary output for growth-rate censoring correction:
#   mle_delay_params_censored.csv  ← onset-to-sample interval-censored MLE params
#   (loaded by growth-rate notebooks to compute P(detected by TRUNC_DATE | onset))
#
# All figures: <OUTPUT_DIR>/<window_name>/
# ============================================================================

# ── renv fallback ────────────────────────────────────────────────────────────
# On a repo that activates renv (e.g. the main branch) the project library may
# not contain the CFR dependencies. renv replaces .libPaths() (and R_LIBS_USER),
# so if a core package is missing, append the user's default library that holds
# it — letting these scripts run identically in or out of renv.
if (!requireNamespace("tidyverse", quietly = TRUE)) local({
  for (l in c(Sys.glob(file.path(path.expand("~"), "Library", "R", "*", "*", "library")),
              Sys.glob(file.path(path.expand("~"), "R", "*-library", "*")),
              Sys.glob(file.path(path.expand("~"), "R", "*", "*", "library"))))
    if (dir.exists(file.path(l, "tidyverse"))) { .libPaths(c(.libPaths(), l)); break }
})

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(fitdistrplus)
  library(patchwork)
  library(scales)
  library(viridis)
  library(ggridges)
  library(gridExtra)
})
# EpiDist (Bayesian) fitting is optional — only needed when RUN_EPIDIST=TRUE. Its
# packages can be absent on a repo that only consumes the censored-MLE params
# (08/09 read those, not the EpiDist output), so load them lazily.
.HAVE_EPIDIST <- requireNamespace("epidist", quietly = TRUE) &&
                 requireNamespace("epiparameter", quietly = TRUE)
if (.HAVE_EPIDIST) suppressPackageStartupMessages({
  library(epidist); library(epiparameter)
})

select <- dplyr::select
filter <- dplyr::filter
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Configuration ──────────────────────────────────────────────────────────────
# Robustly locate config_r.R regardless of the working directory when the script
# is invoked (IDE, Rscript from project root, Rscript from scripts/, etc.).
#
# Search order for config_r.R:
#   1. Same directory as the script (works when --file= is resolved correctly)
#   2. scripts/ subdirectory of cwd (cwd = realtime_pipeline/)
#   3. realtime_pipeline/scripts/ subdirectory of cwd (cwd = project root)
.script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f    <- args[startsWith(args, "--file=")]
  if (length(f) > 0)
    normalizePath(dirname(sub("--file=", "", f)), mustWork = FALSE)
  else
    getwd()
}, error = function(e) getwd())

.config_candidates <- c(
  file.path(.script_dir,  "config_r.R"),
  file.path(getwd(),      "scripts",              "config_r.R"),
  file.path(getwd(),      "realtime_pipeline",
            "scripts",   "config_r.R")
)
.config_path <- Filter(file.exists, .config_candidates)[1]
if (length(.config_path) == 0 || is.na(.config_path))
  stop(paste(
    "Cannot locate config_r.R. Tried:\n",
    paste(.config_candidates, collapse = "\n  "),
    "\nRun from the project root or realtime_pipeline/scripts/ directory."
  ))

source(.config_path, local = TRUE, chdir = FALSE)
cat(sprintf("Sourced config: %s\n", .config_path))

MAX_DELAY <- 60L

# EpiDist toggle — Bayesian fitting (~30-60 min per window). Now the DEFAULT when the
# sourced 00_config.R sets RUN_EPIDIST=TRUE (gated on the epidist/epiparameter packages);
# a run without that config falls back to env RUN_EPIDIST (default FALSE). An explicit
# env var overrides the config value, so the real-time delay refresh can still skip it
# (the censored-MLE params 08/09 consume are written BEFORE this stage) via RUN_EPIDIST=FALSE.
.RUN_EPIDIST_DEFAULT <- isTRUE(get0("RUN_EPIDIST", ifnotfound = FALSE))
RUN_EPIDIST <- (tolower(trimws(Sys.getenv("RUN_EPIDIST",
                 if (.RUN_EPIDIST_DEFAULT) "true" else "false"))) %in%
                c("true", "t", "1", "yes", "y")) && .HAVE_EPIDIST

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

PALETTE <- list(
  confirmed = "#E74C3C", probable  = "#E67E22", negative  = "#95A5A6",
  gamma     = "#E74C3C", lognormal = "#3498DB", weibull   = "#2ECC71",
  exponential = "#F39C12", marginal  = "#8E44AD", naive     = "#E67E22",
  lit       = "#1ABC9C", bunia     = "#E74C3C",
  pos_col   = "#E74C3C", neg_col   = "#3498DB",
  inv_col   = "#F39C12", blk_col   = "#95A5A6"
)

cat(strrep("=", 80), "\n")
cat("04c_delay_four_windows.R — MVE ITURI 2026 (v7 pipeline)\n")
cat(sprintf("Naive + Censored MLE + EpiDist | %s\n", Sys.time()))
cat(strrep("=", 80), "\n\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# ── A. Naive uncensored MLE ───────────────────────────────────────────────────
fit_mle_naive <- function(x, delay_name) {
  x     <- x[!is.na(x) & x >= 0]
  x_pos <- x[x > 0]
  if (length(x_pos) < 5) {
    cat(sprintf("  [%s] n_pos=%d < 5, skipping naive MLE\n",
                delay_name, length(x_pos)))
    return(NULL)
  }
  out <- list()
  for (fam in c("gamma", "lnorm", "weibull", "exp")) {
    data_fit <- if (fam == "lnorm") x_pos + 0.5 else x_pos
    tryCatch({
      fit <- fitdist(data_fit, fam, method = "mle")
      out[[fam]] <- list(family = fam, params = fit$estimate,
                         aic = fit$aic, bic = fit$bic, fit = fit)
    }, error = function(e)
      cat(sprintf("  [%s] %s naive fit failed: %s\n",
                  delay_name, fam, e$message)))
  }
  out
}

naive_implied_mean <- function(fam, params) {
  tryCatch(switch(fam,
    gamma   = unname(params["shape"] / params["rate"]),
    lnorm   = unname(exp(params["meanlog"] + 0.5 * params["sdlog"]^2) - 0.5),
    weibull = unname(params["scale"] * gamma(1 + 1 / params["shape"])),
    exp     = unname(1 / params["rate"]),
    NA_real_), error = function(e) NA_real_)
}

naive_cdf <- function(fam, params, t) {
  tryCatch(switch(fam,
    gamma   = pgamma(t,       shape   = params["shape"],   rate  = params["rate"]),
    lnorm   = plnorm(t + 0.5, meanlog = params["meanlog"], sdlog = params["sdlog"]),
    weibull = pweibull(t,     shape   = params["shape"],   scale = params["scale"]),
    exp     = pexp(t,         rate    = params["rate"]),
    NA_real_), error = function(e) NA_real_)
}

plot_mle_fits <- function(x, mle_fits, n_total, title_str, xlim_max = NULL) {
  x_clean <- x[!is.na(x) & x >= 0]
  if (length(x_clean) == 0)
    return(ggplot() + labs(title = title_str) + theme_void())
  if (is.null(xlim_max)) xlim_max <- min(max(x_clean) + 2, MAX_DELAY)
  x_seq  <- seq(0.01, xlim_max, by = 0.1)
  colours <- c(gamma = PALETTE$gamma, lnorm = PALETTE$lognormal,
                weibull = PALETTE$weibull, exp = PALETTE$exponential)
  p <- ggplot(data.frame(delay = x_clean), aes(x = delay)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                   fill = "grey80", colour = "grey50", alpha = 0.7) +
    xlim(0, xlim_max) +
    labs(title = sprintf("%s (n=%d)", title_str, n_total),
         x = "Delay (days)", y = "Density") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
  if (is.null(mle_fits) || length(mle_fits) == 0) return(p)
  for (fam in names(mle_fits)) {
    r    <- mle_fits[[fam]]
    dens <- tryCatch(switch(fam,
      gamma   = dgamma(x_seq,        shape   = r$params["shape"],   rate  = r$params["rate"]),
      lnorm   = dlnorm(x_seq + 0.5,  meanlog = r$params["meanlog"], sdlog = r$params["sdlog"]),
      weibull = dweibull(x_seq,      shape   = r$params["shape"],   scale = r$params["scale"]),
      exp     = dexp(x_seq,          rate    = r$params["rate"])
    ), error = function(e) rep(NA_real_, length(x_seq)))
    p <- p + geom_line(
      data = data.frame(x = x_seq, y = dens,
                        family = toupper(fam)),
      aes(x = x, y = y, colour = family, linetype = family),
      linewidth = 1.1)
  }
  aic_v <- sapply(mle_fits, function(r) round(r$aic, 1))
  labs_ <- sprintf("%s (AIC=%.1f)", names(aic_v), aic_v)
  cols  <- toupper(unname(colours[names(mle_fits)]))
  p + scale_colour_manual(
    values = setNames(cols, toupper(names(mle_fits))),
    labels = setNames(labs_, toupper(names(mle_fits))), name = "Distribution") +
    scale_linetype_manual(
      values = setNames(c("solid", "dashed", "dotdash", "longdash"),
                        toupper(names(mle_fits)))[seq_along(mle_fits)],
      labels = setNames(labs_, toupper(names(mle_fits))), name = "Distribution") +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
}

.mle_row <- function(fits, delay_type, window_name, n_obs) {
  if (is.null(fits) || length(fits) == 0) return(NULL)
  purrr::map_dfr(fits, function(r) {
    p <- r$params
    tibble(window = window_name, delay_type = delay_type,
           family = toupper(r$family), n = n_obs,
           aic = round(r$aic, 4), bic = round(r$bic, 4),
           implied_mean_d = round(naive_implied_mean(r$family, p), 3),
           shape   = if ("shape"   %in% names(p)) unname(p["shape"])   else NA_real_,
           rate    = if ("rate"    %in% names(p)) unname(p["rate"])    else NA_real_,
           scale   = if ("scale"   %in% names(p)) unname(p["scale"])   else NA_real_,
           meanlog = if ("meanlog" %in% names(p)) unname(p["meanlog"]) else NA_real_,
           sdlog   = if ("sdlog"   %in% names(p)) unname(p["sdlog"])   else NA_real_)
  }) %>%
    group_by(delay_type) %>%
    mutate(delta_aic = round(aic - min(aic), 4), best = delta_aic == 0) %>%
    ungroup() %>%
    arrange(delay_type, aic)
}

# ── B. Interval-censored MLE ──────────────────────────────────────────────────
.make_cens_df <- function(x) {
  stopifnot(is.numeric(x), all(x >= 0, na.rm = TRUE))
  data.frame(left  = ifelse(x == 0L, 0, as.numeric(x) - 0.5),
             right = as.numeric(x) + 0.5)
}

.get_cens_starts <- function(x, fam) {
  m <- mean(x, na.rm = TRUE); v <- var(x, na.rm = TRUE)
  if (!is.finite(m) || m <= 0) m <- 1
  if (!is.finite(v) || v <= 0) v <- m^2
  switch(fam,
    gamma = {
      sh <- m^2 / v; ra <- m / v
      if (!is.finite(sh) || sh <= 0) sh <- 1
      if (!is.finite(ra) || ra <= 0) ra <- 1 / m
      list(shape = sh, rate = ra)
    },
    lnorm = {
      ml <- log(m^2 / sqrt(v + m^2)); sl <- sqrt(log(1 + v / m^2))
      if (!is.finite(ml) || !is.finite(sl) || sl <= 0) { ml <- log(m); sl <- 0.5 }
      list(meanlog = ml, sdlog = sl)
    },
    weibull = {
      cv <- sqrt(v) / m
      k  <- if (is.finite(cv) && cv > 0 && cv < 10) cv^(-1.086) else 1
      sc <- tryCatch(m / gamma(1 + 1 / k), error = function(e) m)
      if (!is.finite(k) || k <= 0) k <- 1
      if (!is.finite(sc) || sc <= 0) sc <- m
      list(shape = k, scale = sc)
    },
    exp = list(rate = 1 / m)
  )
}

.fit_one_censored <- function(x, fam, delay_name = "") {
  if (length(x) < 5L) return(NULL)
  cdf <- .make_cens_df(x)
  m <- mean(x)
  # FAMILY-APPROPRIATE multi-start: gamma-shaped fallbacks (shape/rate) are invalid
  # parameter names for weibull (shape/scale) and lnorm (meanlog/sdlog) and made the
  # multi-start illusory (they errored, leaving only the primary start). Give each
  # family its own alternative starts so a failed primary can still recover.
  sl <- switch(fam,
    exp     = list(list(rate = 1 / m)),
    gamma   = list(.get_cens_starts(x, fam), list(shape = 1, rate = 1 / m), list(shape = 2, rate = 2 / m)),
    lnorm   = list(.get_cens_starts(x, fam), list(meanlog = log(m), sdlog = 0.5),
                   list(meanlog = log(max(stats::median(x), 0.5)), sdlog = 1)),
    weibull = list(.get_cens_starts(x, fam), list(shape = 1, scale = m), list(shape = 1.5, scale = m)),
    list(.get_cens_starts(x, fam)))
  for (s in sl) {
    fit <- tryCatch(
      suppressWarnings(fitdistcens(cdf, fam, start = s)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      p <- fit$estimate
      # lnorm's meanlog is legitimately negative when the median delay < 1 day, so
      # test only the strictly-positive parameters (sdlog); all others must be > 0.
      pos <- if (fam == "lnorm") p[["sdlog"]] else p
      if (all(is.finite(p)) && is.finite(fit$aic) && all(pos > 0)) return(fit)
    }
  }
  message(sprintf("  [%s/%s] all censored starts failed", delay_name, fam))
  NULL
}

.fit_all_censored <- function(x, delay_name) {
  x <- x[!is.na(x) & x >= 0L]
  if (length(x) < 5L) {
    cat(sprintf("  [%s] n=%d<5, skipping censored\n", delay_name, length(x)))
    return(NULL)
  }
  cat(sprintf("  %-34s n=%-4d  n_zero=%-4d (%.1f%%)\n",
              delay_name, length(x), sum(x == 0L),
              sum(x == 0L) / length(x) * 100))
  out <- list()
  for (fam in c("gamma", "lnorm", "weibull", "exp")) {
    fit <- .fit_one_censored(x, fam, delay_name)
    if (!is.null(fit)) out[[fam]] <- fit
  }
  out
}

.implied_mean <- function(fam, params) {
  tryCatch(switch(fam,
    gamma   = unname(params["shape"] / params["rate"]),
    lnorm   = unname(exp(params["meanlog"] + 0.5 * params["sdlog"]^2)),
    weibull = unname(params["scale"] * gamma(1 + 1 / params["shape"])),
    exp     = unname(1 / params["rate"]),
    NA_real_), error = function(e) NA_real_)
}

.cens_cdf <- function(fam, params, t) {
  tryCatch(switch(fam,
    gamma   = pgamma(t,   shape   = params["shape"],   rate  = params["rate"]),
    lnorm   = plnorm(t,   meanlog = params["meanlog"], sdlog = params["sdlog"]),
    weibull = pweibull(t, shape   = params["shape"],   scale = params["scale"]),
    exp     = pexp(t,     rate    = params["rate"]),
    NA_real_), error = function(e) NA_real_)
}

.density_from_fit <- function(fam, params, x_grid) {
  tryCatch(switch(fam,
    gamma   = dgamma(x_grid,   shape   = params["shape"],   rate  = params["rate"]),
    lnorm   = dlnorm(x_grid,   meanlog = params["meanlog"], sdlog = params["sdlog"]),
    weibull = dweibull(x_grid, shape   = params["shape"],   scale = params["scale"]),
    exp     = dexp(x_grid,     rate    = params["rate"]),
    rep(NA_real_, length(x_grid))),
  error = function(e) rep(NA_real_, length(x_grid)))
}

.build_cens_tbl <- function(fits, delay_type, window_name, n_obs) {
  if (is.null(fits) || length(fits) == 0L) return(NULL)
  purrr::map_dfr(names(fits), function(fam) {
    fit <- fits[[fam]]; p <- fit$estimate
    tibble(window = window_name, delay_type = delay_type,
           family = toupper(fam), n = n_obs,
           aic = round(fit$aic, 4), bic = round(fit$bic, 4),
           implied_mean_d = round(.implied_mean(fam, p), 3),
           shape   = if ("shape"   %in% names(p)) round(unname(p["shape"]), 4)   else NA_real_,
           rate    = if ("rate"    %in% names(p)) round(unname(p["rate"]), 4)    else NA_real_,
           scale   = if ("scale"   %in% names(p)) round(unname(p["scale"]), 4)   else NA_real_,
           meanlog = if ("meanlog" %in% names(p)) round(unname(p["meanlog"]), 4) else NA_real_,
           sdlog   = if ("sdlog"   %in% names(p)) round(unname(p["sdlog"]), 4)   else NA_real_)
  }) %>%
    mutate(delta_aic = round(aic - min(aic), 4), best = delta_aic == 0) %>%
    arrange(aic)
}

.plot_censored_fits <- function(x, fits, panel_title, xlim_max = NULL) {
  x_c <- x[!is.na(x) & x >= 0L]
  if (length(x_c) == 0)                    # no usable delays -> empty panel, no crash
    return(ggplot() + labs(title = panel_title) + theme_void())
  if (is.null(xlim_max)) xlim_max <- min(max(x_c) + 2L, MAX_DELAY)
  xg  <- seq(0.01, xlim_max, by = 0.1)
  p   <- ggplot(data.frame(delay = x_c), aes(x = delay)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                   fill = "grey85", colour = "grey60", alpha = 0.75, linewidth = 0.3) +
    xlim(0, xlim_max) +
    labs(title = panel_title,
         subtitle = sprintf("n=%d  (zeros=%d, %.1f%%)",
                            length(x_c), sum(x_c == 0L),
                            sum(x_c == 0L) / length(x_c) * 100),
         x = "Delay (days)", y = "Density") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8, colour = "grey45"))
  if (is.null(fits) || length(fits) == 0L) return(p)
  aic_ord <- sort(sapply(fits, function(f) f$aic))
  top2    <- names(aic_ord)[seq_len(min(2L, length(aic_ord)))]
  cols2   <- c("#d35400", "#1a6faf")
  for (i in seq_along(top2)) {
    fam  <- top2[[i]]
    dens <- .density_from_fit(fam, fits[[fam]]$estimate, xg)
    mu   <- .implied_mean(fam, fits[[fam]]$estimate)
    lbl  <- sprintf("%s (AIC=%.1f, mean=%.2fd)",
                    toupper(fam), fits[[fam]]$aic, mu)
    p    <- p + geom_line(
      data = data.frame(x = xg, y = dens, label = lbl),
      aes(x = x, y = y, colour = label), linewidth = 1.1, na.rm = TRUE)
  }
  leg <- setNames(cols2[seq_along(top2)],
                  sapply(top2, function(fam)
                    sprintf("%s (AIC=%.1f, mean=%.2fd)", toupper(fam),
                            fits[[fam]]$aic,
                            .implied_mean(fam, fits[[fam]]$estimate))))
  p + scale_colour_manual(values = leg, name = NULL) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
}

.cens_best_fam  <- function(f)
  if (!is.null(f) && length(f) > 0)
    names(sort(sapply(f, `[[`, "aic")))[1] else NA_character_
.cens_best_mean <- function(f) {
  b <- .cens_best_fam(f)
  if (is.na(b)) NA_real_ else .implied_mean(b, f[[b]]$estimate)
}
.naive_best_fam  <- function(f)
  if (!is.null(f) && length(f) > 0)
    names(sort(sapply(f, `[[`, "aic")))[1] else NA_character_
.naive_best_mean <- function(f) {
  b <- .naive_best_fam(f)
  if (is.na(b)) NA_real_ else naive_implied_mean(b, f[[b]]$params)
}

# ── C. Bayesian EpiDist ───────────────────────────────────────────────────────
fit_epidist_model <- function(df_pairs, pdate_col, sdate_col,
                               obs_date_val, model_type, family_fn,
                               family_name, delay_name,
                               chains = STAN_CHAINS, iter = STAN_ITER,
                               warmup = STAN_WARMUP, cores = STAN_CORES) {
  tryCatch({
    df_input <- df_pairs %>%
      select(pdate_lwr = all_of(pdate_col), sdate_lwr = all_of(sdate_col)) %>%
      filter(!is.na(pdate_lwr), !is.na(sdate_lwr),
             as.numeric(sdate_lwr - pdate_lwr) >= 0) %>%
      mutate(pdate_upr = pdate_lwr + 1L, sdate_upr = sdate_lwr + 1L,
             obs_date  = pmax(obs_date_val, sdate_lwr + 1L))
    if (model_type == "naive" && family_name %in% c("lognormal", "gamma"))
      df_input <- df_input %>%
        filter(as.numeric(sdate_lwr - pdate_lwr) > 0)
    if (nrow(df_input) < 5) {
      cat(sprintf("  [%s %s %s] n=%d insufficient\n",
                  delay_name, model_type, family_name, nrow(df_input)))
      return(NULL)
    }
    n_complete <- sum(df_input$sdate_lwr <= obs_date_val, na.rm = TRUE)
    frac_c     <- n_complete / nrow(df_input)
    if (frac_c < 0.30) {
      cat(sprintf(
        "  [%s %s %s] skipped — %.0f%% right-truncated\n",
        delay_name, model_type, family_name, (1 - frac_c) * 100))
      return(NULL)
    }
    cat(sprintf("  Fitting %s %s %s (n=%d, %.0f%% complete)...\n",
                delay_name, model_type, family_name,
                nrow(df_input), frac_c * 100))
    ll_obj <- df_input %>%
      as_epidist_linelist_data(
        pdate_lwr = "pdate_lwr", pdate_upr = "pdate_upr",
        sdate_lwr = "sdate_lwr", sdate_upr = "sdate_upr",
        obs_date  = "obs_date")
    model_obj <- switch(model_type,
      naive    = as_epidist_naive_model(ll_obj),
      marginal = as_epidist_aggregate_data(ll_obj) %>%
                   as_epidist_marginal_model())
    fit <- epidist(model_obj, family = family_fn, chains = chains,
                   iter = iter, warmup = warmup, cores = cores,
                   refresh = 0, silent = 2,
                   control = list(adapt_delta = 0.95, max_treedepth = 12))
    list(model_type = model_type, family = family_name, delay = delay_name,
         n = nrow(df_input), frac_complete = frac_c, fit = fit,
         summary = summary(fit))
  }, error = function(e) {
    cat(sprintf("  [%s %s %s] FAILED: %s\n",
                delay_name, model_type, family_name, e$message))
    NULL
  })
}

extract_epidist_summary <- function(fit_obj) {
  if (is.null(fit_obj)) return(NULL)
  tryCatch({
    samps <- predict_delay_parameters(fit_obj$fit) %>% add_mean_sd()
    tibble(delay = fit_obj$delay, model_type = fit_obj$model_type,
           family = fit_obj$family, n = fit_obj$n,
           mean_post_med = round(median(samps$mean, na.rm = TRUE), 3),
           mean_post_lo  = round(quantile(samps$mean, 0.025, na.rm = TRUE), 3),
           mean_post_hi  = round(quantile(samps$mean, 0.975, na.rm = TRUE), 3),
           sd_post_med   = round(median(samps$sd, na.rm = TRUE), 3))
  }, error = function(e)
    tibble(delay = fit_obj$delay, model_type = fit_obj$model_type,
           family = fit_obj$family, n = fit_obj$n,
           mean_post_med = NA_real_, mean_post_lo = NA_real_,
           mean_post_hi  = NA_real_, sd_post_med  = NA_real_))
}

get_posterior_mean_sd <- function(fit_obj) {
  if (is.null(fit_obj)) return(NULL)
  tryCatch({
    samps <- predict_delay_parameters(fit_obj$fit) %>% add_mean_sd()
    list(mean    = median(samps$mean, na.rm = TRUE),
         sd      = median(samps$sd,   na.rm = TRUE),
         mean_lo = quantile(samps$mean, 0.025, names = FALSE),
         mean_hi = quantile(samps$mean, 0.975, names = FALSE))
  }, error = function(e) NULL)
}

fit_both_families <- function(df_ep, delay_name, trunc_date) {
  lapply(c("lognormal", "gamma"), function(fam) {
    fam_fn <- switch(fam, lognormal = lognormal(), gamma = Gamma())
    lapply(c("naive", "marginal"), function(mt)
      fit_epidist_model(df_ep, "pdate_lwr", "sdate_lwr",
                        trunc_date, mt, fam_fn, fam, delay_name))
  }) %>% unlist(recursive = FALSE)
}

# ── C2. Naive vs censored overlay panel (for D03b methods comparison) ─────────
.make_cmp_panel <- function(x_vals, naive_fits, cens_fits, panel_title, xlim_max) {
  x_c  <- x_vals[!is.na(x_vals) & x_vals >= 0L]
  xg   <- seq(0.01, xlim_max, by = 0.1)
  p    <- ggplot(data.frame(delay = x_c), aes(x = delay)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                   fill = "grey85", colour = "grey60", alpha = 0.75, linewidth = 0.3) +
    xlim(0, xlim_max) +
    labs(title = panel_title,
         subtitle = sprintf("n=%d | solid=censored | dashed=naive", length(x_c)),
         x = "Delay (days)", y = "Density") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8, colour = "grey45"),
          legend.position = "bottom", legend.text = element_text(size = 7))
  lines_df <- data.frame()
  if (!is.null(naive_fits) && length(naive_fits) > 0) {
    bn <- .naive_best_fam(naive_fits); pn <- naive_fits[[bn]]$params
    dn <- tryCatch(switch(bn,
      gamma   = dgamma(xg,     shape = pn["shape"],   rate  = pn["rate"]),
      lnorm   = dlnorm(xg+0.5, meanlog = pn["meanlog"], sdlog = pn["sdlog"]),
      weibull = dweibull(xg,   shape = pn["shape"],   scale = pn["scale"]),
      exp     = dexp(xg,       rate  = pn["rate"]),
      rep(NA_real_, length(xg))), error = function(e) rep(NA_real_, length(xg)))
    lbl <- sprintf("Naive %s (mean=%.2fd)", toupper(bn), .naive_best_mean(naive_fits))
    lines_df <- bind_rows(lines_df,
      data.frame(x = xg, y = dn, Method = lbl, ltype = "dashed", stringsAsFactors = FALSE))
  }
  if (!is.null(cens_fits) && length(cens_fits) > 0) {
    bc <- .cens_best_fam(cens_fits)
    dc_y <- .density_from_fit(bc, cens_fits[[bc]]$estimate, xg)
    lbl <- sprintf("Censored %s (mean=%.2fd)", toupper(bc), .cens_best_mean(cens_fits))
    lines_df <- bind_rows(lines_df,
      data.frame(x = xg, y = dc_y, Method = lbl, ltype = "solid", stringsAsFactors = FALSE))
  }
  if (nrow(lines_df) > 0) {
    ltypes <- setNames(lines_df$ltype, lines_df$Method)
    ltypes <- ltypes[!duplicated(names(ltypes))]
    cols   <- setNames(c("#d35400", "#1a6faf")[seq_len(length(ltypes))], names(ltypes))
    p <- p + geom_line(data = lines_df,
                       aes(x = x, y = y, colour = Method, linetype = Method),
                       linewidth = 1.1, na.rm = TRUE) +
      scale_colour_manual(values = cols, name = NULL) +
      scale_linetype_manual(values = ltypes, name = NULL)
  }
  p
}

# ── C3. EpiDist visualisation (for D07) ───────────────────────────────────────
plot_epidist_comparison <- function(fits_list, obs_data, label_str, xlim_max = 25) {
  non_null <- fits_list[!sapply(fits_list, is.null)]
  if (length(non_null) == 0) return(NULL)
  xg      <- seq(0.1, xlim_max, by = 0.1)
  df_hist <- data.frame(delay = obs_data[obs_data >= 0])
  p <- ggplot(df_hist, aes(x = delay)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                   fill = "grey85", colour = "grey60", alpha = 0.8) +
    xlim(0, xlim_max) + labs(title = label_str, x = "Delay (days)", y = "Density") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ld <- purrr::map_dfr(non_null, function(fo) {
    pm <- get_posterior_mean_sd(fo)
    if (is.null(pm) || is.na(pm$mean)) return(NULL)
    m <- pm$mean; s <- pm$sd
    yv <- tryCatch(switch(fo$family,
      lognormal = dlnorm(xg, meanlog = log(m) - 0.5*log(1 + (s/m)^2),
                         sdlog = sqrt(log(1 + (s/m)^2))),
      gamma     = dgamma(xg, shape = (m/s)^2, rate = m/s^2),
      rep(NA_real_, length(xg))), error = function(e) rep(NA_real_, length(xg)))
    data.frame(x = xg, y = yv,
               label = sprintf("%s/%s (mean=%.1fd, n=%d)",
                               fo$model_type, fo$family, pm$mean, fo$n))
  })
  if (!is.null(ld) && nrow(ld) > 0)
    p <- p + geom_line(data = ld,
                       aes(x = x, y = y, colour = label, linetype = label),
                       linewidth = 1.1)
  p
}

# ── D. Temporal analysis ──────────────────────────────────────────────────────
make_temp_plot <- function(df_plot, y_label, title_str) {
  df_plot <- df_plot %>% filter(!is.na(date), !is.na(delay))
  ggplot(df_plot, aes(x = date, y = delay, colour = result_label)) +
    geom_jitter(size = 2.0, alpha = 0.65, width = 0.2, height = 0) +
    scale_colour_manual(
      values = c("Positive" = PALETTE$pos_col, "Negative" = PALETTE$neg_col,
                 "Invalid" = PALETTE$inv_col, "Blank/unknown" = PALETTE$blk_col),
      name = NULL) +
    scale_x_date(date_breaks = "1 week", date_labels = "%d %b") +
    labs(title    = title_str,
         subtitle = sprintf("n=%d | x=symptom onset | LOESS span=0.75", nrow(df_plot)),
         x = "Symptom onset date", y = y_label) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8, colour = "grey45")) +
    { if (nrow(df_plot) >= 6)
        geom_smooth(aes(group = 1), method = "loess", formula = y ~ x,
                    span = 0.75, se = TRUE, colour = "black",
                    fill = "grey70", linewidth = 1, alpha = 0.20)
      else list() }
}

# ── E. Detection completeness ─────────────────────────────────────────────────
.p_not_processed <- function(t_vals, cens_fits) {
  fam    <- .cens_best_fam(cens_fits)
  params <- if (!is.na(fam)) cens_fits[[fam]]$estimate else NULL
  if (is.null(params)) return(rep(NA_real_, length(t_vals)))
  pmax(0, pmin(1, 1 - .cens_cdf(fam, params, t_vals)))
}

# Returns posterior CrI from EpiDist fits for the proportion of cases
# not yet processed. Uses full MCMC draws for proper uncertainty quantification.
.epidist_posterior_detect <- function(epidist_fits, t_vals,
                                       prefer_mt = c("marginal", "naive")) {
  n <- length(t_vals)
  for (mt in prefer_mt) {
    for (fo in epidist_fits) {
      if (is.null(fo) || fo$model_type != mt) next
      samps <- tryCatch(
        predict_delay_parameters(fo$fit) %>% add_mean_sd(),
        error = function(e) NULL
      )
      if (is.null(samps) || nrow(samps) < 10) next
      post_props <- vapply(seq_len(nrow(samps)), function(j) {
        m <- samps$mean[j]; s <- samps$sd[j]
        if (is.na(m) || is.na(s) || m <= 0 || s <= 0) return(NA_real_)
        p_not <- if (fo$family == "lognormal") {
          ml <- log(m^2 / sqrt(s^2 + m^2))
          sl <- sqrt(log(1 + (s/m)^2))
          1 - plnorm(t_vals, meanlog = ml, sdlog = sl)
        } else {
          sh <- (m/s)^2; ra <- m/s^2
          1 - pgamma(t_vals, shape = sh, rate = ra)
        }
        mean(pmax(0, pmin(1, p_not)))
      }, numeric(1))
      post_props <- post_props[!is.na(post_props)]
      if (length(post_props) < 10) next
      return(list(
        prop_med   = median(post_props),
        prop_lo    = quantile(post_props, 0.025, names = FALSE),
        prop_hi    = quantile(post_props, 0.975, names = FALSE),
        exp_med    = median(post_props) * n,
        exp_lo     = quantile(post_props, 0.025, names = FALSE) * n,
        exp_hi     = quantile(post_props, 0.975, names = FALSE) * n,
        n_draws    = length(post_props),
        model_type = mt,
        family     = fo$family
      ))
    }
  }
  NULL
}

# compute_detection_completeness: now accepts EpiDist fits for posterior CrIs.
compute_detection_completeness <- function(df_window, trunc_date,
                                            cens_fits_os,
                                            cens_fits_ob    = NULL,
                                            epidist_fits_os = NULL,
                                            epidist_fits_ob = NULL) {
  df <- df_window %>%
    filter(!is.na(symptom_onset), symptom_onset <= trunc_date) %>%
    mutate(
      t_available  = as.integer(trunc_date - symptom_onset),
      has_sample   = !is.na(sample_date)   & sample_date   <= trunc_date,
      has_analysis = !is.na(analysis_date) & analysis_date <= trunc_date
    )
  n_onset    <- nrow(df)
  n_sampled  <- sum(df$has_sample)
  n_analyzed <- sum(df$has_analysis)
  t_vals     <- df$t_available

  df <- df %>% mutate(
    p_not_sampled_cens  = .p_not_processed(t_available, cens_fits_os),
    p_not_analyzed_cens = if (!is.null(cens_fits_ob))
      .p_not_processed(t_available, cens_fits_ob)
    else rep(NA_real_, n_onset)
  )
  ed_os <- .epidist_posterior_detect(epidist_fits_os %||% list(), t_vals)
  ed_ob <- .epidist_posterior_detect(epidist_fits_ob %||% list(), t_vals)

  n_not_entered <- n_onset - n_analyzed

  list(
    df              = df,
    trunc_date      = trunc_date,
    n_onset         = n_onset,
    n_sampled_obs          = n_sampled,
    n_unsampled_obs        = n_onset - n_sampled,
    prop_unsampled_emp     = if (n_onset > 0) (n_onset - n_sampled) / n_onset else NA_real_,
    exp_unsampled_cens     = sum(df$p_not_sampled_cens,  na.rm = TRUE),
    prop_unsampled_cens    = mean(df$p_not_sampled_cens, na.rm = TRUE),
    ed_os                  = ed_os,
    n_analyzed_obs         = n_analyzed,
    n_not_analyzed_obs     = n_not_entered,
    prop_not_analyzed_emp  = if (n_onset > 0) n_not_entered / n_onset else NA_real_,
    exp_not_analyzed_cens  = sum(df$p_not_analyzed_cens,  na.rm = TRUE),
    prop_not_analyzed_cens = mean(df$p_not_analyzed_cens, na.rm = TRUE),
    ed_ob                  = ed_ob,
    n_not_entered_obs      = n_not_entered,
    prop_not_entered_emp   = if (n_onset > 0) n_not_entered / n_onset else NA_real_,
    best_cens_fam_os = .cens_best_fam(cens_fits_os),
    best_cens_fam_ob = if (!is.null(cens_fits_ob)) .cens_best_fam(cens_fits_ob) else NA_character_
  )
}

# ── E2. Full 5-panel detection completeness figure ────────────────────────────
plot_detection_completeness <- function(dc, window_name) {
  df <- dc$df

  .add_ed_ribbon <- function(p, df_plot, ed_result) {
    if (is.null(ed_result)) return(p)
    x_range <- range(df_plot$symptom_onset, na.rm = TRUE)
    p + annotate("rect",
                 xmin = x_range[1], xmax = x_range[2],
                 ymin = ed_result$prop_lo, ymax = ed_result$prop_hi,
                 fill = "#8E44AD", alpha = 0.12) +
      geom_hline(yintercept = ed_result$prop_med,
                 colour = "#8E44AD", lty = "dashed", linewidth = 0.9)
  }

  p_a <- NULL
  if ("p_not_sampled_cens" %in% names(df) && any(!is.na(df$p_not_sampled_cens))) {
    p_a <- ggplot(df, aes(x = symptom_onset, y = p_not_sampled_cens)) +
      geom_col(fill = "#E74C3C", alpha = 0.55, width = 0.9) +
      geom_line(colour = "#C0392B", linewidth = 0.9) +
      geom_hline(yintercept = 0.5, lty = 2, colour = "grey50", linewidth = 0.4) +
      scale_x_date(date_breaks = "1 week", date_labels = "%d %b") +
      scale_y_continuous(labels = scales::percent_format(1), limits = c(0, 1)) +
      labs(title = "A. P(not yet sampled) by onset date",
           subtitle = sprintf("Censored %s | Dashed+ribbon = EpiDist 95%% CrI",
                              toupper(dc$best_cens_fam_os)),
           x = "Symptom onset date", y = "P(not sampled by extraction date)") +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(),
            plot.subtitle = element_text(size = 8, colour = "grey45"))
    p_a <- .add_ed_ribbon(p_a, df, dc$ed_os)
  }

  p_a2 <- NULL
  if ("p_not_analyzed_cens" %in% names(df) && any(!is.na(df$p_not_analyzed_cens))) {
    p_a2 <- ggplot(df, aes(x = symptom_onset, y = p_not_analyzed_cens)) +
      geom_col(fill = "#E67E22", alpha = 0.55, width = 0.9) +
      geom_line(colour = "#CA6F1E", linewidth = 0.9) +
      geom_hline(yintercept = 0.5, lty = 2, colour = "grey50", linewidth = 0.4) +
      scale_x_date(date_breaks = "1 week", date_labels = "%d %b") +
      scale_y_continuous(labels = scales::percent_format(1), limits = c(0, 1)) +
      labs(title = "A2. P(not yet analyzed / detected) by onset date",
           subtitle = sprintf("Censored %s | Dashed+ribbon = EpiDist 95%% CrI",
                              toupper(if (!is.na(dc$best_cens_fam_ob)) dc$best_cens_fam_ob else "?")),
           x = "Symptom onset date", y = "P(no analysis result by extraction date)") +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(),
            plot.subtitle = element_text(size = 8, colour = "grey45"))
    p_a2 <- .add_ed_ribbon(p_a2, df, dc$ed_ob)
  }

  df_b <- df %>%
    mutate(status = case_when(has_sample ~ "Sampled", TRUE ~ "Not yet sampled"),
           status = factor(status, levels = c("Sampled", "Not yet sampled"))) %>%
    count(symptom_onset, status) %>%
    complete(symptom_onset = seq(min(df$symptom_onset), max(df$symptom_onset), "day"),
             status, fill = list(n = 0))
  p_b <- ggplot(df_b, aes(x = symptom_onset, y = n, fill = status)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = c("Sampled" = "#27AE60", "Not yet sampled" = "#E74C3C"),
                      name = NULL) +
    scale_x_date(date_breaks = "1 week", date_labels = "%d %b") +
    labs(title = "B. Daily onset counts: sampled vs not yet sampled",
         x = "Symptom onset date", y = "Cases") +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  df_c <- df %>% arrange(symptom_onset) %>%
    mutate(cum_unsampled_cens    = cumsum(replace_na(p_not_sampled_cens, 0)),
           cum_not_analyzed_cens = cumsum(replace_na(p_not_analyzed_cens, 0)),
           cum_unsampled_obs     = cumsum(!has_sample),
           cum_not_analyzed_obs  = cumsum(!has_analysis))
  p_c <- ggplot(df_c, aes(x = symptom_onset)) +
    geom_line(aes(y = cum_unsampled_obs,  colour = "Emp. unsampled"),  linewidth = 0.9) +
    geom_line(aes(y = cum_not_analyzed_obs, colour = "Emp. unanalyzed"), linewidth = 0.9, lty = "dotted") +
    geom_line(aes(y = cum_unsampled_cens, colour = "Model unsampled (cens)"), linewidth = 0.9, lty = "dashed") +
    geom_line(aes(y = cum_not_analyzed_cens, colour = "Model unanalyzed (cens)"), linewidth = 0.9, lty = "longdash") +
    scale_colour_manual(
      values = c("Emp. unsampled" = "#E74C3C", "Emp. unanalyzed" = "#C0392B",
                 "Model unsampled (cens)" = "#1a6faf", "Model unanalyzed (cens)" = "#2980B9"),
      name = NULL) +
    scale_x_date(date_breaks = "1 week", date_labels = "%d %b") +
    labs(title = "C. Cumulative expected unsampled and unanalyzed",
         x = "Symptom onset date", y = "Cumulative cases") +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          legend.text = element_text(size = 8))

  .fmt_ed <- function(ed_res) {
    if (is.null(ed_res)) return("  (EpiDist not run or no convergence)")
    sprintf("  EpiDist (%s/%s, n_draws=%d):\n    Prop = %.1f%%  95%% CrI [%.1f%%, %.1f%%]\n    Exp  = %.1f  95%% CrI [%.1f, %.1f]",
            ed_res$model_type, ed_res$family, ed_res$n_draws,
            100*ed_res$prop_med, 100*ed_res$prop_lo, 100*ed_res$prop_hi,
            ed_res$exp_med, ed_res$exp_lo, ed_res$exp_hi)
  }
  summ_txt <- paste0(
    sprintf("Window: %s | Extraction: %s\n", window_name, dc$trunc_date),
    sprintf("Onset cases: %d\n\n", dc$n_onset),
    "--- (1) NOT YET SAMPLED (onset→sample) ---\n",
    sprintf("  Empirical:   %d / %d  (%.1f%%)\n",
            dc$n_unsampled_obs, dc$n_onset, 100*dc$prop_unsampled_emp),
    sprintf("  Censored MLE (%s): %.1f exp  (%.1f%%)\n",
            toupper(dc$best_cens_fam_os), dc$exp_unsampled_cens, 100*dc$prop_unsampled_cens),
    .fmt_ed(dc$ed_os), "\n\n",
    "--- (2) NOT YET ANALYZED / DETECTED ---\n",
    sprintf("  Empirical:   %d / %d  (%.1f%%)\n",
            dc$n_not_analyzed_obs, dc$n_onset, 100*dc$prop_not_analyzed_emp),
    sprintf("  Censored MLE (%s): %.1f exp  (%.1f%%)\n",
            toupper(if (!is.na(dc$best_cens_fam_ob)) dc$best_cens_fam_ob else "N/A"),
            dc$exp_not_analyzed_cens, 100*dc$prop_not_analyzed_cens),
    .fmt_ed(dc$ed_ob), "\n\n",
    "--- (3) NOT YET ENTERED (no analysis_date) ---\n",
    sprintf("  Empirical:   %d / %d  (%.1f%%)",
            dc$n_not_entered_obs, dc$n_onset, 100*dc$prop_not_entered_emp)
  )
  p_d <- ggplot() +
    annotate("text", x = 0.05, y = 0.95, label = summ_txt,
             size = 2.9, hjust = 0, vjust = 1, family = "mono") +
    labs(title = "D. Detection completeness summary") +
    theme_void() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

  list(p_a = p_a, p_a2 = p_a2, p_b = p_b, p_c = p_c, p_d = p_d)
}

# ── F. Panel-B scatter+inset figures ─────────────────────────────────────────
# make_inset_cens: histogram inset with top-2 censored MLE curves (orange/blue)
# and optional EpiDist posterior (dashed purple).
make_inset_cens <- function(vals, cens_fits, epidist_fits,
                             x_label, xlim_max, base_sz = 7) {
  vals_v <- vals[!is.na(vals) & vals >= 0]
  xg     <- seq(0.01, xlim_max, by = 0.1)

  p_in <- ggplot(data.frame(delay = vals_v), aes(x = delay)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                   fill = "grey88", colour = "grey65", alpha = 0.8, linewidth = 0.3) +
    xlim(0, xlim_max) +
    labs(x = x_label, y = "Density") +
    theme_bw(base_size = base_sz) +
    theme(panel.grid        = element_blank(),
          axis.text         = element_text(size = base_sz - 1),
          legend.position   = "bottom",
          legend.key.size   = unit(0.25, "cm"),
          legend.text       = element_text(size = base_sz - 2),
          legend.spacing.y  = unit(0, "pt"),
          legend.background = element_rect(fill = "white", colour = NA),
          plot.background   = element_rect(fill = "white", colour = "grey70",
                                           linewidth = 0.4),
          panel.background  = element_rect(fill = "white"))

  line_rows <- list()

  if (!is.null(cens_fits) && length(cens_fits) > 0) {
    aic_sorted <- sort(sapply(cens_fits, function(f) f$aic))
    n_show     <- min(2L, length(aic_sorted))
    top_fams   <- names(aic_sorted)[seq_len(n_show)]
    cens_cols  <- c("#d35400", "#1a6faf")
    aic_best   <- aic_sorted[[1]]
    for (i in seq_len(n_show)) {
      fam  <- top_fams[[i]]
      dc_v <- .density_from_fit(fam, cens_fits[[fam]]$estimate, xg)
      mu   <- round(.implied_mean(fam, cens_fits[[fam]]$estimate), 1)
      daic <- round(aic_sorted[[i]] - aic_best, 1)
      lbl  <- if (i == 1L) sprintf("Censored %s (%.1fd)", toupper(fam), mu)
              else sprintf("Censored %s (%.1fd, ΔAIC=%.1f)", toupper(fam), mu, daic)
      line_rows[[sprintf("cens%d", i)]] <- data.frame(
        x = xg, y = dc_v, Fit = lbl, ltype = "solid", col = cens_cols[[i]],
        stringsAsFactors = FALSE)
    }
  }

  if (!is.null(epidist_fits) && length(epidist_fits) > 0) {
    for (mt_try in c("marginal", "naive")) {
      done <- FALSE
      for (fo in epidist_fits) {
        if (is.null(fo) || fo$model_type != mt_try) next
        pm <- get_posterior_mean_sd(fo)
        if (is.null(pm) || is.na(pm$mean) || pm$mean <= 0 || pm$sd <= 0) next
        m <- pm$mean; s <- pm$sd
        ded <- tryCatch(switch(fo$family,
          lognormal = dlnorm(xg,
                             meanlog = log(m^2 / sqrt(s^2 + m^2)),
                             sdlog   = sqrt(log(1 + (s/m)^2))),
          gamma     = dgamma(xg, shape = (m/s)^2, rate = m/s^2),
          rep(NA_real_, length(xg))),
          error = function(e) rep(NA_real_, length(xg)))
        if (all(is.na(ded))) next
        lbl <- sprintf("EpiDist %s/%s (%.1fd)", fo$model_type, fo$family, pm$mean)
        line_rows[["ed"]] <- data.frame(
          x = xg, y = ded, Fit = lbl, ltype = "dashed", col = "#7B2D8B",
          stringsAsFactors = FALSE)
        done <- TRUE; break
      }
      if (done) break
    }
  }

  if (length(line_rows) == 0) return(p_in)
  all_lines  <- bind_rows(line_rows)
  lbl_levels <- unique(all_lines$Fit)
  ltype_map  <- setNames(all_lines$ltype[!duplicated(all_lines$Fit)], lbl_levels)
  col_map    <- setNames(all_lines$col[  !duplicated(all_lines$Fit)], lbl_levels)

  p_in +
    geom_line(data = all_lines,
              aes(x = x, y = y, colour = Fit, linetype = Fit),
              linewidth = 0.9, na.rm = TRUE) +
    scale_colour_manual(values   = col_map,   name = NULL, breaks = lbl_levels) +
    scale_linetype_manual(values = ltype_map, name = NULL, breaks = lbl_levels) +
    guides(colour   = guide_legend(nrow = length(lbl_levels)),
           linetype = guide_legend(nrow = length(lbl_levels)))
}

make_scatter_panel <- function(df_delays, date_col, x_label, y_label,
                                subtitle, cens_fits, epidist_fits,
                                inset_xlbl, inset_xlim, panel_tag = NULL) {
  pal <- c("Positive" = PALETTE$pos_col, "Negative" = PALETTE$neg_col,
           "Invalid"           = PALETTE$inv_col, "Blank/unknown" = PALETTE$blk_col)
  df_plot <- df_delays %>%
    dplyr::select(date = all_of(date_col), delay, result_label) %>%
    dplyr::filter(!is.na(date), !is.na(delay))
  stopifnot(inherits(df_plot$date, "Date"), all(df_plot$delay >= 0, na.rm = TRUE))
  p <- ggplot(df_plot, aes(x = date, y = delay, colour = result_label)) +
    geom_jitter(size = 2.0, alpha = 0.65, width = 0.25, height = 0) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_x_date(date_breaks = "1 week", date_labels = "%d %b") +
    labs(subtitle = subtitle, x = x_label, y = y_label) +
    theme_bw(base_size = 10) +
    theme(legend.position  = "top",
          panel.grid.minor = element_blank(),
          plot.subtitle    = element_text(size = 8, colour = "grey45"))
  if (nrow(df_plot) >= 6)
    p <- p + geom_smooth(aes(group = 1), method = "loess", formula = y ~ x,
                          span = 0.75, se = TRUE, colour = "#2c3e50",
                          fill = "#2c3e50", linewidth = 1.0, alpha = 0.10)
  if (!is.null(panel_tag))
    p <- p + labs(tag = panel_tag) +
      theme(plot.tag = element_text(size = 12, face = "bold"))
  ins <- make_inset_cens(df_delays$delay, cens_fits, epidist_fits,
                          inset_xlbl, inset_xlim)
  p + inset_element(ins, left = 0.55, bottom = 0.46, right = 1.0, top = 1.0)
}

# ============================================================================
# DATA LOADING
# ============================================================================
cat("Loading:", DATA_PATH, "\n")
if (!file.exists(DATA_PATH))
  stop("Processed CSV not found — run 00_process_mve_linelist.py first.")

df_all <- read_csv(DATA_PATH, show_col_types = FALSE) %>%
  mutate(across(c(symptom_onset, sample_date, receipt_date, analysis_date,
                  analysis_date_bunia, analysis_date_kinshasa, death_date),
                as.Date)) %>%
  mutate(
    death_date   = if_else(!is.na(death_date) & death_date == as.Date("1970-01-01"),
                           as.Date(NA), death_date),
    age_group    = cut(age_years,
                       breaks = c(0, 5, 10, 20, 30, 40, 50, 60, Inf),
                       labels = c("0-4","5-9","10-19","20-29","30-39","40-49","50-59","60+"),
                       right = FALSE, include.lowest = TRUE),
    gender_clean = case_when(
      toupper(gender) %in% c("M", "MALE", "MASCULIN")   ~ "Male",
      toupper(gender) %in% c("F", "FEMALE", "FEMININ")  ~ "Female",
      TRUE ~ NA_character_),
    result_label = case_when(
      test_positive_either == TRUE ~ "Positive",
      is_negative_strict   == TRUE ~ "Negative",
      is_invalid_strict    == TRUE ~ "Invalid",
      TRUE                         ~ "Blank/unknown"),
    result_label = factor(result_label,
                          levels = c("Positive", "Negative",
                                     "Invalid", "Blank/unknown"))
  )

cat(sprintf("  Loaded: %d records\n\n", nrow(df_all)))

# ── Derive dynamic analysis dates ─────────────────────────────────────────────
# TRUNC_DATE = max(sample_date) − 5 days.
# The final 5 days before the extraction date are heavily right-censored
# (onset→sample delay means recent cases haven't arrived yet), so delay
# distribution fitting stops 5 days before the linelist end.  This matches
# the TRAIN_CUTOFF convention in 15/16.
#
# RIGHT-TRUNCATION CAVEAT (audit H1): the 5-day buffer is adequate for the SHORT
# delays (onset→sample mean ~4 d, sample→receipt <1 d) but is NOT sufficient for
# the LONG onset→death delay (fitted mean ~7 d, q99 ~24 d): df_od below conditions
# on death OBSERVED by the snapshot, so late-onset cases that will die at long
# delays are not yet observed and are excluded → the frequentist onset→death fit
# is right-truncated and biased SHORT. The interval-censored MLE here corrects
# daily ROUNDING only, NOT this truncation; only the EpiDist path (RUN_EPIDIST,
# off by default) applies a truncation correction. The downstream cCFR consumer
# (08_cfr.R) documents this and leans on literature-prior delays as the upward
# sensitivity. Decision: documented, NOT re-estimated.
TEST_DAYS  <- 5L
.dyn       <- get_dynamic_dates(df_all)
OBS_DATE   <- .dyn$obs_date
TRUNC_DATE <- .dyn$trunc_date - TEST_DAYS   # max(sample_date) − 5

cat(sprintf("  ANALYSIS_START (fixed)           : %s\n", ANALYSIS_START))
cat(sprintf("  OBS_DATE       (max all dates)   : %s\n", OBS_DATE))
cat(sprintf("  max(sample_date)                 : %s\n", .dyn$trunc_date))
cat(sprintf("  TRUNC_DATE     (max − %dd)        : %s\n\n",
            TEST_DAYS, TRUNC_DATE))

# ── Build single analytical window ────────────────────────────────────────────
WINDOW_NAME    <- sprintf("analytical_%s", format(TRUNC_DATE, "%Y-%m-%d"))
cur_out_dir    <- file.path(OUTPUT_DIR, WINDOW_NAME)
dir.create(cur_out_dir, showWarnings = FALSE, recursive = TRUE)

df <- df_all %>%
  filter(!is.na(symptom_onset),
         symptom_onset >= ANALYSIS_START,
         symptom_onset <= TRUNC_DATE)

cat(sprintf("%s\nWINDOW: %s  (%s – %s)\n%s\n",
            strrep("=", 70), WINDOW_NAME, ANALYSIS_START, TRUNC_DATE,
            strrep("=", 70)))
cat(sprintf("Records after filter: %d of %d\n\n", nrow(df), nrow(df_all)))

pos_either <- df %>% filter(test_positive_either == TRUE)
pos_conf   <- df %>% filter(confirmed_final == TRUE)
probable   <- df %>% filter(probable_initial == TRUE)
likely     <- df %>% filter(total_likely == TRUE)

# ── Delay populations ─────────────────────────────────────────────────────────
df_ob <- df %>%
  filter(!is.na(symptom_onset), !is.na(analysis_date)) %>%
  mutate(delay = as.integer(analysis_date - symptom_onset)) %>%
  filter(delay >= 0, delay <= MAX_DELAY)

df_os <- df %>%
  filter(!is.na(symptom_onset), !is.na(sample_date)) %>%
  mutate(delay = as.integer(sample_date - symptom_onset)) %>%
  filter(delay >= 0, delay <= MAX_DELAY)

# NB (audit M2): onset->death uses `death_date > symptom_onset`, which EXCLUDES
# same-day (delay==0) deaths, unlike the other delays which keep delay==0 (and the
# interval-censored fit maps d=0 -> [0,0.5]). Same-day recorded onset-and-death is
# biologically implausible for a filovirus and usually a data artefact (onset
# recorded = admission near death), so it is dropped here deliberately; the
# asymmetry is intentional and is now logged in the drop-accounting block below.
df_od <- df %>%
  filter(!is.na(symptom_onset), !is.na(death_date),
         died == TRUE, death_date > symptom_onset) %>%
  mutate(delay = as.integer(death_date - symptom_onset)) %>%
  filter(delay >= 0, delay <= MAX_DELAY)

df_ra <- df %>%
  filter(!is.na(receipt_date), !is.na(analysis_date)) %>%
  mutate(delay = as.integer(analysis_date - receipt_date)) %>%
  filter(delay >= 0, delay <= MAX_DELAY)

df_sr <- df %>%
  filter(!is.na(sample_date), !is.na(receipt_date)) %>%
  mutate(delay = as.integer(receipt_date - sample_date)) %>%
  filter(delay >= 0, delay <= MAX_DELAY)

df_rb <- df %>%
  filter(!is.na(receipt_date), !is.na(analysis_date_bunia)) %>%
  mutate(delay = as.integer(analysis_date_bunia - receipt_date)) %>%
  filter(delay >= 0, delay <= MAX_DELAY)

# sample -> first analysis (sample collection to first lab result). Added for the
# main-text Figure 1D; same interval-censored MLE machinery and windowing as the
# other delays. `analysis_date` is the FIRST analysis date (Bunia, else Kinshasa).
df_sa <- df %>%
  filter(!is.na(sample_date), !is.na(analysis_date)) %>%
  mutate(delay = as.integer(analysis_date - sample_date)) %>%
  filter(delay >= 0, delay <= MAX_DELAY)

for (info in list(
    list(d = df_ob, nm = "onset->analysis"),
    list(d = df_os, nm = "onset->sample"),
    list(d = df_od, nm = "onset->death"),
    list(d = df_ra, nm = "receipt->analysis"))) {
  vv <- info$d$delay
  if (length(vv) == 0)
    cat(sprintf("  %-24s  n=0\n", info$nm))
  else
    cat(sprintf("  %-24s  n=%-4d  mean=%.1fd  median=%.0fd  [%d,%d]\n",
                info$nm, length(vv), mean(vv), median(vv), min(vv), max(vv)))
}

# ── Explicit drop accounting (audit M1) ───────────────────────────────────────
# Previously every delay population silently dropped missing-date, negative, zero
# (onset->death only) and >MAX_DELAY rows. We now log the full breakdown per delay
# so the surviving N is reconciled (n_both = n_neg + n_zero_dropped + n_used +
# n_over_max). No filtering logic is changed - this is reporting only.
cat("\n  -- delay drop accounting (within-window) ----------------------------\n")
.drop_acct <- function(nm, d1, d2, drop_zero = FALSE) {
  both <- df %>% filter(!is.na(.data[[d1]]), !is.na(.data[[d2]]))
  if (nm == "onset->death") both <- both %>% filter(died == TRUE)
  dl   <- as.integer(both[[d2]] - both[[d1]])
  n_both <- length(dl)
  n_neg  <- sum(dl < 0, na.rm = TRUE)
  n_zero <- sum(dl == 0, na.rm = TRUE)
  n_big  <- sum(dl > MAX_DELAY, na.rm = TRUE)
  z_drop <- if (drop_zero) n_zero else 0L     # onset->death drops same-day
  n_used <- n_both - n_neg - n_big - z_drop
  cat(sprintf(paste0("  %-20s n_both=%-4d used=%-4d  dropped: neg=%d%s ",
                     ">%dd=%d  (reconciles=%s)\n"),
              nm, n_both, n_used, n_neg,
              if (drop_zero) sprintf(" sameday=%d", n_zero) else "",
              MAX_DELAY, n_big,
              identical(n_used + n_neg + n_big + z_drop, n_both)))
}
.drop_acct("onset->analysis",   "symptom_onset", "analysis_date")
.drop_acct("onset->sample",     "symptom_onset", "sample_date")
.drop_acct("onset->death",      "symptom_onset", "death_date", drop_zero = TRUE)
.drop_acct("receipt->analysis", "receipt_date",  "analysis_date")
.drop_acct("sample->receipt",   "sample_date",   "receipt_date")
.drop_acct("sample->analysis",  "sample_date",   "analysis_date")

onset_bunia_vals    <- df_ob$delay
onset_sample_vals   <- df_os$delay
onset_death_vals    <- df_od$delay
receipt_first_vals  <- df_ra$delay
sample_receipt_vals <- df_sr$delay

# ============================================================================
# UNIT TEST CHECKS — BATCH 1: DATA INTEGRITY
# ============================================================================
.n_pass <- 0L; .n_fail <- 0L
chk <- function(desc, condition, got = NULL, expected = NULL) {
  tag <- if (isTRUE(condition)) {
    .n_pass <<- .n_pass + 1L; "PASS"
  } else {
    .n_fail <<- .n_fail + 1L; "FAIL"
  }
  detail <- if (!is.null(got)) {
    if (!is.null(expected)) sprintf("  got=%s expected=%s", got, expected)
    else sprintf("  got=%s", got)
  } else ""
  cat(sprintf("  [%s] %s%s\n", tag, desc, detail))
  invisible(condition)
}

cat("-- CHECKS 1: DATA INTEGRITY ----------------------------------------\n")
chk("TRUNC_DATE > ANALYSIS_START", TRUNC_DATE > ANALYSIS_START,
    got = as.character(TRUNC_DATE))
chk("Records in window > 0", nrow(df) > 0, got = nrow(df))
n_part <- sum(df$test_positive_either, na.rm = TRUE) +
          sum(df$is_negative_strict,   na.rm = TRUE) +
          sum(df$is_invalid_strict,    na.rm = TRUE) +
          sum(df$is_blank,             na.rm = TRUE)
chk("Partition pos+neg+inv+blank==total", n_part == nrow(df),
    got = n_part, expected = nrow(df))
chk("No epoch death dates",
    sum(df$death_date == as.Date("1970-01-01"), na.rm = TRUE) == 0,
    got = sum(df$death_date == as.Date("1970-01-01"), na.rm = TRUE))
for (col in c("symptom_onset", "sample_date", "receipt_date",
              "analysis_date", "death_date"))
  chk(sprintf("%s is Date", col), inherits(df[[col]], "Date"),
      got = class(df[[col]])[1])
onset_v <- df$symptom_onset[!is.na(df$symptom_onset)]
if (length(onset_v) > 0) {
  chk(sprintf("All onset >= %s", ANALYSIS_START),
      all(onset_v >= ANALYSIS_START), got = as.character(min(onset_v)))
  chk(sprintf("All onset <= %s", TRUNC_DATE),
      all(onset_v <= TRUNC_DATE),     got = as.character(max(onset_v)))
}
chk("onset->sample n > 0",   nrow(df_os) > 0, got = nrow(df_os))
chk("onset->analysis n > 0", nrow(df_ob) > 0, got = nrow(df_ob))
chk("All onset->sample delays in [0,MAX_DELAY]",
    all(df_os$delay >= 0) && all(df_os$delay <= MAX_DELAY))
age_v <- df$age_years[!is.na(df$age_years)]
if (length(age_v) > 0)
  chk("age_years in [0,120]",
      all(age_v >= 0) && all(age_v <= 120),
      got = sprintf("min=%.1f max=%.1f", min(age_v), max(age_v)))
ct_raw <- suppressWarnings(as.numeric(df$altona_ct_value))
ct_v   <- ct_raw[!is.na(ct_raw) & ct_raw > 0]
if (length(ct_v) > 0)
  chk(sprintf("Altona CT in [%.0f,%.0f]", CT_MIN, CT_MAX),
      all(ct_v >= CT_MIN) && all(ct_v <= CT_MAX),
      got = sprintf("min=%.2f max=%.2f", min(ct_v), max(ct_v)))
cat(sprintf("  BATCH 1: %d PASS, %d FAIL\n\n", .n_pass, .n_fail))

# ============================================================================
# D01: EPIDEMIOLOGICAL OVERVIEW
# ============================================================================
cat("-- D01: EPIDEMIOLOGICAL OVERVIEW --\n")
onset_daily <- df %>%
  filter(!is.na(symptom_onset)) %>%
  count(symptom_onset, result_label) %>%
  complete(symptom_onset = seq(min(symptom_onset), max(symptom_onset), "day"),
           result_label, fill = list(n = 0))

p_epi <- ggplot(onset_daily, aes(x = symptom_onset, y = n, fill = result_label)) +
  geom_col(width = 0.9) +
  scale_fill_manual(
    values = c("Positive" = PALETTE$pos_col, "Negative" = PALETTE$neg_col,
               "Invalid" = PALETTE$inv_col, "Blank/unknown" = PALETTE$blk_col),
    name = NULL) +
  scale_x_date(date_breaks = "1 week", date_labels = "%d %b",
               expand = expansion(add = 0.5)) +
  geom_vline(xintercept = TRUNC_DATE,
             lty = 3, colour = "grey40", linewidth = 0.5) +
  annotate("rect", xmin = ANALYSIS_START, xmax = TRUNC_DATE,
           ymin = -Inf, ymax = Inf, alpha = 0.07, fill = "grey30") +
  labs(title = "A. Epidemic curve — all records by result category",
       subtitle = sprintf("Window: %s – %s  |  shaded = analysis window",
                          ANALYSIS_START, TRUNC_DATE),
       x = "Symptom onset date", y = "Cases") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ct_data <- pos_conf %>%
  mutate(ct = suppressWarnings(as.numeric(altona_ct_value))) %>%
  filter(!is.na(ct), ct > 0, ct < 50)
p_ct <- if (nrow(ct_data) > 0) {
  ggplot(ct_data, aes(x = analysis_date_kinshasa, y = ct)) +
    geom_point(colour = PALETTE$confirmed, size = 2.5, alpha = 0.8) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1,
                formula = y ~ x, span = 0.8, colour = "darkred") +
    geom_hline(yintercept = 35, lty = 2, colour = "grey50", linewidth = 0.5) +
    labs(title = "B. Altona CT (Kinshasa, confirmed)",
         x = "Kinshasa analysis date", y = "CT") +
    theme_bw(base_size = 11)
} else ggplot() + annotate("text", x = 0, y = 0, label = "No CT data") + theme_void()

p_pyr <- {
  pyr <- pos_conf %>%
    filter(!is.na(age_group), !is.na(gender_clean)) %>%
    count(age_group, gender_clean) %>%
    mutate(n_plot = if_else(gender_clean == "Male", -n, n))
  ggplot(pyr, aes(x = age_group, y = n_plot, fill = gender_clean)) +
    geom_col(width = 0.85) + coord_flip() +
    scale_y_continuous(labels = function(x) abs(x)) +
    scale_fill_manual(values = c("Male" = "#3498DB", "Female" = "#E74C3C"), name = NULL) +
    geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.5) +
    labs(title = "C. Age-sex pyramid (confirmed)", x = "Age group", y = "Cases") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
}
p_zone <- {
  zc <- pos_conf %>% count(health_zone) %>% arrange(desc(n)) %>%
    mutate(health_zone = factor(health_zone, levels = health_zone))
  ggplot(zc, aes(x = health_zone, y = n, fill = n)) +
    geom_col(colour = "white") +
    geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
    scale_fill_viridis_c(option = "plasma", direction = -1, guide = "none") +
    labs(title = "D. Confirmed cases by health zone", x = NULL, y = "Cases") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          panel.grid.minor = element_blank())
}
pdf(file.path(cur_out_dir, "D01_v7_epidemiological_overview.pdf"),
    width = 18, height = 14)
print((p_epi | p_ct) / (p_pyr | p_zone) + plot_annotation(
  title    = sprintf("MVE Ituri BDBV 2026 — Epidemiological Overview (%s)", WINDOW_NAME),
  subtitle = sprintf("N=%d | pos-either=%d | confirmed=%d | probable=%d",
                     nrow(df), nrow(pos_either), nrow(pos_conf), nrow(probable)),
  theme    = theme(plot.title = element_text(size = 14, face = "bold"))))
dev.off(); cat("D01_v7_epidemiological_overview.pdf\n")

# ============================================================================
# D02: DELAY OVERVIEW
# ============================================================================
cat("-- D02: DELAY OVERVIEW --\n")
make_delay_long <- function(df_d, lbl)
  df_d %>% transmute(delay_days = delay, delay_label = lbl, result_label)
delays_long <- bind_rows(
  make_delay_long(df_ob, "Onset → first analysis"),
  make_delay_long(df_os, "Onset → sample"),
  make_delay_long(df_od, "Onset → death"),
  make_delay_long(df_ra, "Receipt → first analysis"),
  make_delay_long(df_sr, "Sample → receipt"),
  make_delay_long(df_rb, "Receipt → Bunia only"))

FOUR_DELAYS <- c("Onset → first analysis", "Onset → sample",
                 "Onset → death", "Receipt → first analysis")
p_ov <- delays_long %>%
  filter(delay_label %in% FOUR_DELAYS) %>%
  mutate(delay_label = factor(delay_label, levels = FOUR_DELAYS)) %>%
  ggplot(aes(x = delay_days, fill = result_label)) +
  geom_histogram(binwidth = 1, alpha = 0.8, position = "stack",
                 colour = "white", linewidth = 0.2) +
  facet_wrap(~delay_label, scales = "free", ncol = 2) +
  scale_fill_manual(
    values = c("Positive" = PALETTE$pos_col, "Negative" = PALETTE$neg_col,
               "Invalid" = PALETTE$inv_col, "Blank/unknown" = PALETTE$blk_col),
    name = "Result") +
  labs(title = sprintf(
         "Delay distributions (%s)\nn: analysis=%d  sample=%d  death=%d  receipt->analysis=%d",
         WINDOW_NAME, nrow(df_ob), nrow(df_os), nrow(df_od), nrow(df_ra)),
       x = "Delay (days)", y = "Count") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
pdf(file.path(cur_out_dir, "D02_v7_delay_overview.pdf"), width = 14, height = 10)
print(p_ov); dev.off(); cat("D02_v7_delay_overview.pdf\n")

# D02b: Bunia vs Kinshasa lab turnaround comparison
if (nrow(df_rb) > 0 || nrow(df_ra) > 0) {
  df_bk <- bind_rows(
    df_rb %>% select(delay, result_label) %>% mutate(lab = "Bunia (receipt→analysis)"),
    df_ra %>% select(delay, result_label) %>% mutate(lab = "Kinshasa (receipt→analysis)")
  )
  p_d02b <- ggplot(df_bk, aes(x = delay, fill = lab)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1, alpha = 0.7,
                   position = "identity", colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("Bunia (receipt→analysis)"    = PALETTE$bunia,
                                 "Kinshasa (receipt→analysis)" = "#3498DB"),
                      name = NULL) +
    facet_wrap(~lab, ncol = 1) +
    labs(title  = sprintf("Bunia vs Kinshasa turnaround — %s", WINDOW_NAME),
         x = "Delay (days)", y = "Density") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
  pdf(file.path(cur_out_dir, "D02b_v7_bunia_vs_kinshasa_delays.pdf"), width = 12, height = 8)
  print(p_d02b); dev.off(); cat("D02b_v7_bunia_vs_kinshasa_delays.pdf\n")
}

# ============================================================================
# SECTION 3a: NAIVE UNCENSORED MLE
# ============================================================================
cat("-- S3a: NAIVE UNCENSORED MLE --\n")
mle_ob <- fit_mle_naive(onset_bunia_vals,    "onset->analysis")
mle_os <- fit_mle_naive(onset_sample_vals,   "onset->sample  ")
mle_od <- fit_mle_naive(onset_death_vals,    "onset->death   ")
mle_ra <- fit_mle_naive(receipt_first_vals,  "receipt->analysis")
mle_sr <- fit_mle_naive(sample_receipt_vals, "sample->receipt")

for (info in list(
    list(f = mle_ob, nm = "Onset->analysis",   n = length(onset_bunia_vals)),
    list(f = mle_os, nm = "Onset->sample",     n = length(onset_sample_vals)),
    list(f = mle_od, nm = "Onset->death",      n = length(onset_death_vals)),
    list(f = mle_ra, nm = "Receipt->analysis", n = length(receipt_first_vals)))) {
  cat(sprintf("\n  %s (n=%d):\n", info$nm, info$n))
  for (fn in names(info$f)) {
    r <- info$f[[fn]]
    cat(sprintf("    %-8s AIC=%.1f  mean=%.2fd  params: %s\n",
                fn, r$aic, naive_implied_mean(fn, r$params),
                paste(round(r$params, 3), names(r$params),
                      sep = "=", collapse = ", ")))
  }
}

.mle_tbl <- bind_rows(
  .mle_row(mle_os, "onset_sample",          WINDOW_NAME, length(onset_sample_vals)),
  .mle_row(mle_ob, "onset_first_analysis",  WINDOW_NAME, length(onset_bunia_vals)),
  .mle_row(mle_od, "onset_death",           WINDOW_NAME, length(onset_death_vals)),
  .mle_row(mle_ra, "receipt_first_analysis",WINDOW_NAME, length(receipt_first_vals)),
  .mle_row(mle_sr, "sample_receipt",        WINDOW_NAME, length(sample_receipt_vals))
)
write_csv(.mle_tbl, file.path(cur_out_dir, "mle_delay_params.csv"))
cat(sprintf("\nSaved mle_delay_params.csv (%d rows)\n", nrow(.mle_tbl)))

# Naive MLE checks
cat("-- CHECKS 2: NAIVE MLE --------------------------------------------\n")
.n_pass <- 0L; .n_fail <- 0L
for (mi in list(
    list(f = mle_ob, nm = "onset->analysis",   lo = 5,   hi = 25),
    list(f = mle_os, nm = "onset->sample",     lo = 1,   hi = 15),
    list(f = mle_od, nm = "onset->death",      lo = 2,   hi = 20),
    list(f = mle_ra, nm = "receipt->analysis", lo = 0.1, hi = 15))) {
  if (!is.null(mi$f) && length(mi$f) > 0) {
    chk(sprintf("%s: all AIC finite", mi$nm),
        all(is.finite(sapply(mi$f, `[[`, "aic"))))
    bn <- .naive_best_fam(mi$f); mu <- .naive_best_mean(mi$f)
    chk(sprintf("%s: best(%s) mean in [%.1f,%.1f]",
                mi$nm, toupper(bn), mi$lo, mi$hi),
        !is.na(mu) && mu >= mi$lo && mu <= mi$hi, got = round(mu, 2))
    chk(sprintf("%s: all params > 0", mi$nm),
        all(sapply(mi$f, function(r) all(r$params > 0, na.rm = TRUE))))
  } else {
    chk(sprintf("%s: at least 1 family fitted", mi$nm), FALSE, got = "NULL")
  }
}
cat(sprintf("  BATCH 2: %d PASS, %d FAIL\n\n", .n_pass, .n_fail))

# Naive MLE figures
if (!is.null(mle_ob)) {
  pdf(file.path(cur_out_dir, "D03_v7_mle_onset_analysis.pdf"), width = 14, height = 7)
  print(plot_mle_fits(onset_bunia_vals, mle_ob, nrow(df_ob),
                       "Onset→first analysis: Naive MLE", 30))
  dev.off()
}
if (!is.null(mle_os)) {
  pdf(file.path(cur_out_dir, "D04_v7_mle_onset_sample.pdf"), width = 14, height = 7)
  print(plot_mle_fits(onset_sample_vals, mle_os, nrow(df_os),
                       "Onset→sample: Naive MLE", 20))
  dev.off()
}
if (!is.null(mle_od)) {
  pdf(file.path(cur_out_dir, "D04b_v7_mle_onset_death.pdf"), width = 14, height = 7)
  print(plot_mle_fits(onset_death_vals, mle_od, nrow(df_od),
                       "Onset→death: Naive MLE", 25))
  dev.off()
}
if (!is.null(mle_ra)) {
  pdf(file.path(cur_out_dir, "D05_v7_mle_receipt_analysis.pdf"), width = 14, height = 7)
  print(plot_mle_fits(receipt_first_vals, mle_ra, nrow(df_ra),
                       "Receipt→analysis: Naive MLE", 20))
  dev.off()
}
cat("Naive MLE figures saved.\n\n")

# D03_combined: all four naive MLE fits in one 2×2 grid
{
  p_comb_list <- Filter(Negate(is.null), list(
    if (!is.null(mle_ob)) plot_mle_fits(onset_bunia_vals,   mle_ob, nrow(df_ob),
                                         "A. Onset→first analysis", 30) else NULL,
    if (!is.null(mle_os)) plot_mle_fits(onset_sample_vals,  mle_os, nrow(df_os),
                                         "B. Onset→sample", 20) else NULL,
    if (!is.null(mle_od)) plot_mle_fits(onset_death_vals,   mle_od, nrow(df_od),
                                         "C. Onset→death", 25) else NULL,
    if (!is.null(mle_ra)) plot_mle_fits(receipt_first_vals, mle_ra, nrow(df_ra),
                                         "D. Receipt→first analysis", 20) else NULL
  ))
  if (length(p_comb_list) >= 2) {
    pdf(file.path(cur_out_dir, "D03_combined_v7_mle_all_delays.pdf"), width = 16, height = 12)
    print(wrap_plots(p_comb_list, ncol = 2) + plot_annotation(
      title    = sprintf("Naive MLE — all delay types (%s)", WINDOW_NAME),
      subtitle = "AIC-ranked fits; x+0.5 shift for lognormal; zeros excluded",
      theme    = theme(plot.title = element_text(size = 13, face = "bold"))))
    dev.off(); cat("D03_combined_v7_mle_all_delays.pdf\n")
  }
}

# ============================================================================
# SECTION 3b: INTERVAL-CENSORED MLE
# ============================================================================
cat("-- S3b: INTERVAL-CENSORED MLE --\n")
cat("  d=0→[0,0.5];  d>0→[d-0.5,d+0.5];  zeros included\n\n")

cens_ob <- .fit_all_censored(df_ob$delay, "onset->analysis  ")
cens_os <- .fit_all_censored(df_os$delay, "onset->sample    ")
cens_od <- .fit_all_censored(df_od$delay, "onset->death     ")
cens_ra <- .fit_all_censored(df_ra$delay, "receipt->analysis")
cens_sr <- .fit_all_censored(df_sr$delay, "sample->receipt  ")
cens_sa <- .fit_all_censored(df_sa$delay, "sample->analysis ")

.cens_tbl <- bind_rows(
  .build_cens_tbl(cens_ob, "onset_first_analysis",  WINDOW_NAME, nrow(df_ob)),
  .build_cens_tbl(cens_os, "onset_sample",           WINDOW_NAME, nrow(df_os)),
  .build_cens_tbl(cens_od, "onset_death",            WINDOW_NAME, nrow(df_od)),
  .build_cens_tbl(cens_ra, "receipt_first_analysis", WINDOW_NAME, nrow(df_ra)),
  .build_cens_tbl(cens_sr, "sample_receipt",         WINDOW_NAME, nrow(df_sr)),
  .build_cens_tbl(cens_sa, "sample_first_analysis",  WINDOW_NAME, nrow(df_sa))
)

for (dt in unique(.cens_tbl$delay_type)) {
  sub <- .cens_tbl %>% dplyr::filter(delay_type == dt)
  cat(sprintf("\n  %s (n=%d):\n", dt, sub$n[1]))
  sub %>% dplyr::select(family, aic, delta_aic, implied_mean_d,
                         shape, rate, scale, meanlog, sdlog) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 3))) %>%
    print(n = Inf)
}

write_csv(.cens_tbl, file.path(cur_out_dir, "mle_delay_params_censored.csv"))
cat(sprintf("\nSaved mle_delay_params_censored.csv (%d rows)\n", nrow(.cens_tbl)))

# ── Also save to root data/processed/ for growth-rate notebooks ───────────────
# DATA_PATH is defined in config_r.R as <ROOT>/data/processed/<filename>.
# dirname(DATA_PATH) therefore gives the correct data/processed/ directory
# regardless of where OUTPUT_DIR points.
delay_params_out <- gsub("_processed\\.csv$", "_delay_params_censored.csv",
                         DATA_PATH)
# Guard: if DATA_PATH does not end in `_processed.csv` the gsub is a no-op and
# delay_params_out == DATA_PATH, so write_csv would OVERWRITE the source line list
# with this 20-row params table (silent data loss). Refuse to write onto the source.
if (identical(delay_params_out, DATA_PATH))
  stop(sprintf(paste0("[04c] Refusing to write delay params over the source data file ",
                      "(DATA_PATH '%s' does not end in '_processed.csv'). Rename the ",
                      "source or set an explicit output path."), DATA_PATH), call. = FALSE)
dir.create(dirname(delay_params_out), showWarnings = FALSE, recursive = TRUE)
write_csv(.cens_tbl, delay_params_out)
cat(sprintf("Saved (root copy): %s\n\n", delay_params_out))

# Censored MLE checks
cat("-- CHECKS 3: CENSORED MLE ------------------------------------------\n")
.n_pass <- 0L; .n_fail <- 0L
for (ci in list(
    list(f = cens_ob, nm = "onset->analysis",   lo = 3,   hi = 25),
    list(f = cens_os, nm = "onset->sample",     lo = 0.5, hi = 15),
    list(f = cens_od, nm = "onset->death",      lo = 2,   hi = 20),
    list(f = cens_ra, nm = "receipt->analysis", lo = 0.1, hi = 10))) {
  if (!is.null(ci$f) && length(ci$f) > 0) {
    chk(sprintf("%s: >=1 family converged", ci$nm), TRUE,
        got = sprintf("n_fits=%d best=%s",
                      length(ci$f), toupper(.cens_best_fam(ci$f))))
    mu <- .cens_best_mean(ci$f)
    chk(sprintf("%s: best mean in [%.1f,%.1f]", ci$nm, ci$lo, ci$hi),
        !is.na(mu) && mu >= ci$lo && mu <= ci$hi, got = round(mu, 2))
    chk(sprintf("%s: AIC finite", ci$nm),
        is.finite(.cens_tbl %>%
          dplyr::filter(delay_type == gsub("->", "_", ci$nm), best) %>%
          dplyr::pull(aic) %>% .[1] %>% as.numeric()))
  } else {
    chk(sprintf("%s: >=1 family converged", ci$nm), FALSE, got = "NULL")
  }
}
cat(sprintf("  BATCH 3: %d PASS, %d FAIL\n\n", .n_pass, .n_fail))
# D03b: naive vs censored methods comparison (4-panel) — requires cens_* already fitted
{
  p_cmp_list <- Filter(Negate(is.null), list(
    if (!is.null(mle_ob) && !is.null(cens_ob))
      .make_cmp_panel(onset_bunia_vals,   mle_ob, cens_ob, "A. Onset→first analysis", 30) else NULL,
    if (!is.null(mle_os) && !is.null(cens_os))
      .make_cmp_panel(onset_sample_vals,  mle_os, cens_os, "B. Onset→sample",         20) else NULL,
    if (!is.null(mle_od) && !is.null(cens_od))
      .make_cmp_panel(onset_death_vals,   mle_od, cens_od, "C. Onset→death",          25) else NULL,
    if (!is.null(mle_ra) && !is.null(cens_ra))
      .make_cmp_panel(receipt_first_vals, mle_ra, cens_ra, "D. Receipt→analysis",     20) else NULL
  ))
  if (length(p_cmp_list) >= 2) {
    pdf(file.path(cur_out_dir, "D03b_v7_methods_comparison.pdf"), width = 16, height = 12)
    print(wrap_plots(p_cmp_list, ncol = 2) + plot_annotation(
      title    = sprintf("Naive vs Interval-Censored MLE comparison (%s)", WINDOW_NAME),
      subtitle = "Solid = censored best; dashed = naive best | both at implied mean",
      theme    = theme(plot.title = element_text(size = 13, face = "bold"))))
    dev.off(); cat("D03b_v7_methods_comparison.pdf\n")
  }
}


# Censored MLE figures
pdf(file.path(cur_out_dir, "D_censored_mle_fits.pdf"), width = 16, height = 12)
p_cens_list <- list(
  if (!is.null(cens_ob)) .plot_censored_fits(df_ob$delay, cens_ob, "A. Onset→analysis", 30)  else NULL,
  if (!is.null(cens_os)) .plot_censored_fits(df_os$delay, cens_os, "B. Onset→sample",   20)  else NULL,
  if (!is.null(cens_od)) .plot_censored_fits(df_od$delay, cens_od, "C. Onset→death",    25)  else NULL,
  if (!is.null(cens_ra)) .plot_censored_fits(df_ra$delay, cens_ra, "D. Receipt→analysis",20) else NULL
)
p_cens_list <- p_cens_list[!sapply(p_cens_list, is.null)]
if (length(p_cens_list) >= 2)
  print(wrap_plots(p_cens_list, ncol = 2) + plot_annotation(
    title = sprintf("Interval-Censored MLE — %s", WINDOW_NAME),
    subtitle = "Top-2 families by AIC; solid lines. d=0→[0,0.5]; d>0→[d-0.5,d+0.5]",
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey40"))))
dev.off()
cat("D_censored_mle_fits.pdf\n\n")

# ============================================================================
# SECTION 3c: EpiDist (optional — set RUN_EPIDIST <- TRUE to enable)
# ============================================================================
tbl_epidist <- tibble()

if (RUN_EPIDIST) {
  cat("-- S3c: BAYESIAN EPIDIST --\n")
  options(brms.backend = "rstan")
  ed_os <- fit_both_families(
    df_os %>% rename(pdate_lwr = symptom_onset, sdate_lwr = sample_date),
    "onset_sample", TRUNC_DATE)
  ed_ob <- fit_both_families(
    df_ob %>% rename(pdate_lwr = symptom_onset, sdate_lwr = analysis_date),
    "onset_analysis", TRUNC_DATE)
  # onset→death: only fit if we have ≥ 20 deaths with known onset and death dates
  ed_od <- if (nrow(df_od) >= 20L) {
    fit_both_families(
      df_od %>% rename(pdate_lwr = symptom_onset, sdate_lwr = death_date),
      "onset_death", TRUNC_DATE)
  } else {
    cat("  [EpiDist] onset→death: fewer than 20 pairs — skipping\n")
    list()
  }
  tbl_epidist <- bind_rows(
    purrr::map_dfr(ed_os, extract_epidist_summary),
    purrr::map_dfr(ed_ob, extract_epidist_summary),
    purrr::map_dfr(ed_od, extract_epidist_summary)
  )
  if (nrow(tbl_epidist) > 0) {
    write_csv(tbl_epidist, file.path(cur_out_dir, "epidist_estimates.csv"))
    cat("Saved epidist_estimates.csv\n\n")
  }
} else {
  cat("-- S3c: EpiDist SKIPPED (RUN_EPIDIST=FALSE) --\n\n")
  ed_os <- list(); ed_ob <- list(); ed_od <- list()
}

# D07: EpiDist posterior comparison (only when RUN_EPIDIST=TRUE)
if (RUN_EPIDIST && length(ed_os) > 0) {
  cat("-- D07: EPIDIST COMPARISON --\n")
  p_ed_list <- Filter(Negate(is.null), list(
    plot_epidist_comparison(ed_os, onset_sample_vals, "A. Onset→sample (EpiDist)", 25),
    plot_epidist_comparison(ed_ob, onset_bunia_vals,  "B. Onset→first analysis (EpiDist)", 35)
  ))
  if (length(p_ed_list) > 0) {
    pdf(file.path(cur_out_dir, "D07_v7_epidist_comparison.pdf"), width = 14, height = 6)
    print(wrap_plots(p_ed_list, ncol = 2) + plot_annotation(
      title    = sprintf("EpiDist Bayesian posterior — %s", WINDOW_NAME),
      subtitle = "Posterior mean density overlaid on histogram | marginal model preferred",
      theme    = theme(plot.title = element_text(size = 13, face = "bold"))))
    dev.off(); cat("D07_v7_epidist_comparison.pdf\n")
  }
}

# ============================================================================
# D06: TEMPORAL ANALYSIS OF DELAYS
# ============================================================================
cat("-- D06: TEMPORAL DELAY ANALYSIS --\n")
# df_os and df_ob already carry result_label (inherited from df).
# Joining back onto df would create result_label.x / result_label.y suffixes,
# making the column invisible to make_temp_plot. Just rename date directly.
df_os_result <- df_os %>% rename(date = symptom_onset)
df_ob_result <- df_ob %>% rename(date = symptom_onset)

# Verify result_label exists before plotting
stopifnot("result_label" %in% names(df_os_result))
stopifnot("result_label" %in% names(df_ob_result))

p_t_os <- make_temp_plot(df_os_result, "Onset→sample (days)",
                          "Onset-to-sample delay over time")
p_t_ob <- make_temp_plot(df_ob_result, "Onset→analysis (days)",
                          "Onset-to-analysis delay over time")
pdf(file.path(cur_out_dir, "D06_v7_temporal_analysis.pdf"), width = 14, height = 10)
print(p_t_os / p_t_ob + plot_annotation(
  title    = sprintf("Temporal delay trends — %s", WINDOW_NAME),
  subtitle = "LOESS span=0.75 with 95% CI; coloured by result category",
  theme    = theme(plot.title = element_text(size = 13, face = "bold"))))
dev.off(); cat("D06_v7_temporal_analysis.pdf\n\n")

# ============================================================================
# D10: MODEL COMPARISON — extended 3-panel (ΔAIC bar + method comparison + EpiDist)
# ============================================================================
cat("-- D10: MODEL COMPARISON (extended) --\n")

# Panel 1: ΔAIC comparison across all delay types (naive and censored)
if (nrow(.mle_tbl) > 0 && nrow(.cens_tbl) > 0) {
  aic_cmp <- bind_rows(
    .mle_tbl  %>% mutate(method = "Naive MLE"),
    .cens_tbl %>% mutate(method = "Censored MLE")
  ) %>%
    group_by(method, delay_type) %>%
    arrange(aic) %>%
    mutate(delta_aic = aic - min(aic)) %>%
    ungroup()

  p_d10_aic <- ggplot(aic_cmp, aes(x = family, y = delta_aic,
                                    fill = method, group = method)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    facet_wrap(~delay_type, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Naive MLE" = "#e67e22", "Censored MLE" = "#2980b9"),
                      name = NULL) +
    labs(title   = sprintf("ΔAIC by family and method — %s", WINDOW_NAME),
         subtitle = "Lower ΔAIC = better relative fit; ΔAIC = 0 = best in group",
         x = "Distribution family", y = "ΔAIC") +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  # Panel 2: Naive vs censored overlay for onset→sample (single best families)
  p_d10_cmp <- if (!is.null(mle_os) && !is.null(cens_os))
    .make_cmp_panel(onset_sample_vals, mle_os, cens_os,
                    "Onset→sample: Naive vs Censored", 20) else NULL

  # Panel 3: EpiDist posterior (if available)
  p_d10_ed <- if (RUN_EPIDIST && length(ed_os) > 0)
    plot_epidist_comparison(ed_os, onset_sample_vals,
                            "Onset→sample: EpiDist posterior", 20) else NULL

  p_d10_panels <- Filter(Negate(is.null),
                         list(p_d10_aic, p_d10_cmp, p_d10_ed))
  if (length(p_d10_panels) >= 1) {
    pdf(file.path(cur_out_dir, "D10_v7_model_comparison.pdf"), width = 16, height = 12)
    print(wrap_plots(p_d10_panels, ncol = 1) + plot_annotation(
      title    = sprintf("Model comparison — %s", WINDOW_NAME),
      subtitle = sprintf("Top: ΔAIC across all delays | Mid: naive vs censored (onset→sample) | Bot: EpiDist posterior"),
      theme    = theme(plot.title = element_text(size = 13, face = "bold"))))
    dev.off(); cat("D10_v7_model_comparison.pdf\n")
  }
}

# ============================================================================
# D11: DETECTION COMPLETENESS (5-panel: P_sample, P_analyze, bar, cumulative, summary)
# ============================================================================
cat("-- D11: DETECTION COMPLETENESS --\n")
# Pass EpiDist fits when available — they produce posterior CrI bands on panels A and A2
dc_res <- compute_detection_completeness(
  df_window       = df,
  trunc_date      = TRUNC_DATE,
  cens_fits_os    = cens_os,
  cens_fits_ob    = cens_ob,
  epidist_fits_os = if (RUN_EPIDIST) ed_os else NULL,
  epidist_fits_ob = if (RUN_EPIDIST) ed_ob else NULL
)

cat(sprintf("\n  Cases with onset:       %d\n", dc_res$n_onset))
cat(sprintf("  NOT YET SAMPLED (onset->sample):\n"))
cat(sprintf("    Empirical:            %d / %d  (%.1f%%)\n",
            dc_res$n_unsampled_obs, dc_res$n_onset,
            100 * dc_res$prop_unsampled_emp))
cat(sprintf("    Censored MLE (%s):  %.1f exp  (%.1f%%)\n",
            toupper(dc_res$best_cens_fam_os),
            dc_res$exp_unsampled_cens,
            100 * dc_res$prop_unsampled_cens))
cat(sprintf("  NOT YET ANALYZED (onset->analysis):\n"))
cat(sprintf("    Empirical:            %d / %d  (%.1f%%)\n",
            dc_res$n_not_analyzed_obs, dc_res$n_onset,
            100 * dc_res$prop_not_analyzed_emp))
cat(sprintf("    Censored MLE (%s):  %.1f exp  (%.1f%%)\n",
            toupper(dc_res$best_cens_fam_ob %||% "N/A"),
            dc_res$exp_not_analyzed_cens,
            100 * dc_res$prop_not_analyzed_cens))

# Full 5-panel detection completeness figure
dc_panels <- plot_detection_completeness(dc_res, WINDOW_NAME)
dc_pdf <- file.path(cur_out_dir, "D11_v7_detection_completeness.pdf")
# Assemble panels: A (P_sample) and A2 (P_analyze) on top row,
# B (bar chart) and C (cumulative) on middle row, D (text summary) on bottom.
dc_plot_list <- Filter(Negate(is.null),
                       list(dc_panels$p_a, dc_panels$p_a2,
                            dc_panels$p_b, dc_panels$p_c,
                            dc_panels$p_d))
if (length(dc_plot_list) >= 1) {
  pdf(dc_pdf, width = 16, height = if (length(dc_plot_list) >= 4) 16 else 8)
  if (length(dc_plot_list) == 5) {
    print((dc_panels$p_a | dc_panels$p_a2) /
          (dc_panels$p_b | dc_panels$p_c)  /
           dc_panels$p_d +
      plot_annotation(
        title    = sprintf("Detection completeness — %s", WINDOW_NAME),
        subtitle = sprintf("TRUNC_DATE=%s | Censored MLE with EpiDist CrI bands (when available)",
                           TRUNC_DATE),
        theme    = theme(plot.title = element_text(size = 13, face = "bold"))))
  } else {
    print(wrap_plots(dc_plot_list, ncol = 2) + plot_annotation(
      title = sprintf("Detection completeness — %s", WINDOW_NAME)))
  }
  dev.off(); cat("D11_v7_detection_completeness.pdf\n")
}

# ============================================================================
# PANEL-B: SCATTER + INSET FIGURES (forward and backward views)
# ============================================================================
cat("-- PANEL-B: SCATTER+INSET FIGURES --\n")

# Reference EpiDist fits for inset purple line (empty list if not run)
epidist_fits_os <- if (RUN_EPIDIST) ed_os else list()
epidist_fits_ob <- if (RUN_EPIDIST) ed_ob else list()
epidist_fits_od <- if (RUN_EPIDIST && exists("ed_od")) ed_od else list()

# Helper: save PDF + PNG at 200 dpi
.save_panel_2x2 <- function(p1, p2, p3, p4, filepath,
                              ann_subtitle, w = 16, h = 14, dpi = 200) {
  assembled <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      subtitle = ann_subtitle,
      theme    = theme(plot.subtitle = element_text(size = 9, colour = "grey40"),
                       plot.tag      = element_text(size = 12, face = "bold")))
  pdf(filepath, width = w, height = h)
  print(assembled); dev.off()
  cat(sprintf("%s\n", basename(filepath)))
  png_path <- sub("\\.pdf$", ".png", filepath)
  ggsave(png_path, plot = assembled, width = w, height = h,
         units = "in", dpi = dpi, bg = "white")
  cat(sprintf("%s  (%d dpi)\n", basename(png_path), dpi))
}

if (nrow(df_os) >= 6) {
  # panel_B_onset_sample: onset→sample scatter with inset (forward, single panel)
  p_fwd_os <- make_scatter_panel(
    df_delays    = df_os,
    date_col     = "symptom_onset",
    x_label      = "Symptom onset date",
    y_label      = "Onset → sample delay (days)",
    subtitle     = sprintf("%s  |  n=%d  |  LOESS span=0.75", WINDOW_NAME, nrow(df_os)),
    cens_fits    = cens_os,
    epidist_fits = epidist_fits_os,
    inset_xlbl   = "Days (onset→sample)",
    inset_xlim   = 25)
  pdf(file.path(cur_out_dir, "panel_B_onset_sample.pdf"), width = 10, height = 6)
  print(p_fwd_os); dev.off(); cat("panel_B_onset_sample.pdf\n")
}

if (nrow(df_ob) >= 6) {
  # panel_B_onset_analysis: onset→analysis scatter with inset
  p_fwd_oa <- make_scatter_panel(
    df_delays    = df_ob,
    date_col     = "symptom_onset",
    x_label      = "Symptom onset date",
    y_label      = "Onset → first analysis (days)",
    subtitle     = sprintf("%s  |  n=%d  |  LOESS span=0.75", WINDOW_NAME, nrow(df_ob)),
    cens_fits    = cens_ob,
    epidist_fits = epidist_fits_ob,
    inset_xlbl   = "Days (onset→analysis)",
    inset_xlim   = 35)
  pdf(file.path(cur_out_dir, "panel_B_onset_analysis.pdf"), width = 10, height = 6)
  print(p_fwd_oa); dev.off(); cat("panel_B_onset_analysis.pdf\n")
}

# panel_B_all_delays: 2×2 forward-looking grid
# Panel layout:
#   A: onset → sample     (primary surveillance delay; used for growth-rate correction)
#   B: onset → analysis   (full detection pipeline delay)
#   C: receipt → analysis (lab processing delay)
#   D: onset → death      (disease severity marker; replaces sample→receipt)
if (nrow(df_os) >= 6 && nrow(df_ob) >= 6) {
  tryCatch({
    .g_os <- make_scatter_panel(
      df_delays = df_os, date_col = "symptom_onset",
      x_label = "Symptom onset date", y_label = "Onset → sample (days)",
      subtitle = NULL, cens_fits = cens_os, epidist_fits = epidist_fits_os,
      inset_xlbl = "Days (onset→sample)", inset_xlim = 25, panel_tag = "A")
    .g_oa <- make_scatter_panel(
      df_delays = df_ob, date_col = "symptom_onset",
      x_label = "Symptom onset date", y_label = "Onset → analysis (days)",
      subtitle = NULL, cens_fits = cens_ob, epidist_fits = epidist_fits_ob,
      inset_xlbl = "Days (onset→analysis)", inset_xlim = 35, panel_tag = "B")
    .g_ra <- make_scatter_panel(
      df_delays = df_ra, date_col = "receipt_date",
      x_label = "Receipt date", y_label = "Receipt → analysis (days)",
      subtitle = NULL, cens_fits = cens_ra, epidist_fits = list(),
      inset_xlbl = "Days (receipt→analysis)", inset_xlim = 20, panel_tag = "C")
    # Panel D: onset → death (informative for disease severity and case ascertainment)
    # Uses EpiDist fits when available (ed_od); falls back to censored MLE only.
    .g_od_panel <- if (nrow(df_od) >= 6) {
      make_scatter_panel(
        df_delays = df_od, date_col = "symptom_onset",
        x_label = "Symptom onset date", y_label = "Onset → death (days)",
        subtitle = NULL, cens_fits = cens_od, epidist_fits = epidist_fits_od,
        inset_xlbl = "Days (onset→death)", inset_xlim = 35, panel_tag = "D")
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = sprintf("Onset → death\nn = %d (< 6 — not plotted)", nrow(df_od)),
                 size = 4, colour = "grey50") +
        theme_void() + labs(tag = "D") +
        theme(plot.tag = element_text(size = 12, face = "bold"))
    }
    .save_panel_2x2(
      .g_os, .g_oa, .g_ra, .g_od_panel,
      filepath     = file.path(cur_out_dir, "panel_B_all_delays.pdf"),
      ann_subtitle = sprintf(
        "%s  |  Scatter: coloured by result category  |  Inset: censored MLE (± EpiDist purple dashed)",
        WINDOW_NAME))
  }, error = function(e)
    message("  panel_B_all_delays failed: ", e$message))
}

# ============================================================================
# D12: DELAY DISTRIBUTIONS BY TEST OUTCOME
# Compare censored MLE fits for onset→sample and onset→analysis, stratified by
# result_label (Positive / Negative / Invalid / Blank). Highlights differences
# in delay behaviour between test-result groups (e.g., positives sampled earlier
# due to clinical severity vs negatives sampled as part of contact tracing).
# ============================================================================
cat("-- D12: DELAY DISTRIBUTIONS BY TEST OUTCOME --\n")

tryCatch({

  # Helper: fit censored MLE for one stratum and return a tidy summary row
  .fit_cens_stratum <- function(delay_vals, delay_type, outcome_label) {
    v <- delay_vals[!is.na(delay_vals) & delay_vals >= 0L]
    if (length(v) < 8L) return(NULL)
    fits <- .fit_all_censored(v, sprintf("%s/%s", delay_type, outcome_label))
    best <- .cens_best_fam(fits)
    if (is.na(best)) return(NULL)
    m   <- .cens_best_mean(fits)
    med <- median(v, na.rm = TRUE)
    tibble(delay_type = delay_type, outcome = outcome_label, n = length(v),
           best_family = best,
           implied_mean = round(m, 2),
           empirical_median = round(med, 1),
           aic = round(fits[[best]]$aic, 1))
  }

  # df_os and df_ob already carry result_label (inherited from df via filter/mutate)
  df_os_out <- df_os %>% dplyr::mutate(result_label = as.character(result_label))
  df_ob_out <- df_ob %>% dplyr::mutate(result_label = as.character(result_label))

  outcomes <- c("Positive", "Negative", "Invalid", "Blank/unknown")

  # Fit MLE per outcome × delay type
  strata_rows <- list()
  for (out in outcomes) {
    v_os <- df_os_out %>% dplyr::filter(result_label == out) %>% dplyr::pull(delay)
    v_ob <- df_ob_out %>% dplyr::filter(result_label == out) %>% dplyr::pull(delay)
    strata_rows <- c(strata_rows,
      list(.fit_cens_stratum(v_os, "onset_sample",   out)),
      list(.fit_cens_stratum(v_ob, "onset_analysis",  out)))
  }
  strata_tbl <- bind_rows(Filter(Negate(is.null), strata_rows))

  if (nrow(strata_tbl) > 0) {
    write_csv(strata_tbl, file.path(cur_out_dir, "delay_by_outcome.csv"))
    cat("  Saved delay_by_outcome.csv\n")
  }

  # ── Figure: faceted ridge / overlapping density plot ──────────────────────
  # For each delay type, show all outcome-specific censored MLE densities
  # overlaid on the pooled histogram, coloured by outcome.
  outcome_cols <- c(
    "Positive"     = "#e74c3c",
    "Negative"     = "#3498db",
    "Invalid"      = "#f39c12",
    "Blank/unknown"= "#95a5a6")

  .make_outcome_panel <- function(df_delays_out, delay_label, xlim_max,
                                   all_cens_fits, panel_tag) {
    xg <- seq(0.01, xlim_max, by = 0.1)
    # Pooled histogram data
    all_vals <- df_delays_out$delay[!is.na(df_delays_out$delay) &
                                      df_delays_out$delay >= 0]
    p <- ggplot(data.frame(delay = all_vals), aes(x = delay)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                     fill = "grey88", colour = "grey65",
                     alpha = 0.7, linewidth = 0.3) +
      xlim(0, xlim_max) +
      labs(x = sprintf("%s (days)", delay_label), y = "Density",
           tag = panel_tag) +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "top",
            plot.tag = element_text(size = 12, face = "bold"))

    line_data <- purrr::map_dfr(names(all_cens_fits), function(out) {
      cf <- all_cens_fits[[out]]
      if (is.null(cf) || length(cf) == 0) return(NULL)
      fam <- .cens_best_fam(cf)
      if (is.na(fam)) return(NULL)
      params <- cf[[fam]]$estimate
      yv <- tryCatch(.density_from_fit(fam, params, xg),
                     error = function(e) rep(NA_real_, length(xg)))
      mu  <- round(.implied_mean(fam, params), 1)
      n   <- length(df_delays_out$delay[
        !is.na(df_delays_out$delay) &
        df_delays_out$result_label == out])
      data.frame(x = xg, y = yv,
                 Outcome = sprintf("%s (n=%d, μ=%.1fd)", out, n, mu),
                 outcome_key = out,
                 stringsAsFactors = FALSE)
    })

    if (!is.null(line_data) && nrow(line_data) > 0) {
      col_map <- setNames(
        outcome_cols[line_data$outcome_key[!duplicated(line_data$Outcome)]],
        line_data$Outcome[!duplicated(line_data$Outcome)])
      p <- p +
        geom_line(data = line_data,
                  aes(x = x, y = y, colour = Outcome),
                  linewidth = 1.2, na.rm = TRUE) +
        scale_colour_manual(values = col_map, name = NULL) +
        guides(colour = guide_legend(nrow = 2))
    }
    p
  }

  # Build per-outcome cens fits for both delay types
  cens_os_by_out <- setNames(lapply(outcomes, function(out) {
    v <- df_os_out %>% dplyr::filter(result_label == out) %>% dplyr::pull(delay)
    v <- v[!is.na(v) & v >= 0L]
    if (length(v) < 8L) return(NULL)
    .fit_all_censored(v, sprintf("onset_sample/%s", out))
  }), outcomes)

  cens_ob_by_out <- setNames(lapply(outcomes, function(out) {
    v <- df_ob_out %>% dplyr::filter(result_label == out) %>% dplyr::pull(delay)
    v <- v[!is.na(v) & v >= 0L]
    if (length(v) < 8L) return(NULL)
    .fit_all_censored(v, sprintf("onset_analysis/%s", out))
  }), outcomes)

  p_d12_os <- .make_outcome_panel(
    df_os_out %>% mutate(result_label = as.character(result_label)),
    "Onset → sample", xlim_max = 25,
    all_cens_fits = cens_os_by_out, panel_tag = "A")
  p_d12_ob <- .make_outcome_panel(
    df_ob_out %>% mutate(result_label = as.character(result_label)),
    "Onset → analysis", xlim_max = 35,
    all_cens_fits = cens_ob_by_out, panel_tag = "B")

  # AIC comparison facet: ΔAIC within each outcome (which family fits best)
  aic_by_outcome <- purrr::map_dfr(outcomes, function(out) {
    purrr::map_dfr(list(
      list(fits = cens_os_by_out[[out]], delay = "Onset→sample",   src = df_os_out),
      list(fits = cens_ob_by_out[[out]], delay = "Onset→analysis", src = df_ob_out)
    ), function(x) {
      if (is.null(x$fits)) return(NULL)
      aics <- sapply(x$fits, function(f) if (!is.null(f)) f$aic else NA_real_)
      aics <- aics[!is.na(aics)]
      if (length(aics) == 0) return(NULL)
      best_aic <- min(aics)
      n_out <- sum(!is.na(x$src$delay) & x$src$result_label == out)
      tibble(outcome = out, delay_type = x$delay,
             family = toupper(names(aics)),
             delta_aic = round(aics - best_aic, 1),
             n = n_out)
    })
  })

  p_d12_aic <- if (nrow(aic_by_outcome) > 0) {
    ggplot(aic_by_outcome,
           aes(x = family, y = delta_aic, fill = outcome)) +
      geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
      facet_wrap(~ delay_type, ncol = 2, scales = "free_y") +
      scale_fill_manual(values = outcome_cols, name = NULL) +
      labs(title = "ΔAIC by family and outcome",
           subtitle = "Lower = better; ΔAIC=0 = best within each outcome",
           x = "Distribution family", y = "ΔAIC") +
      theme_bw(base_size = 10) +
      theme(legend.position = "top", panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"))
  } else NULL

  panels_d12 <- Filter(Negate(is.null), list(p_d12_os, p_d12_ob, p_d12_aic))
  if (length(panels_d12) >= 2) {
    layout_d12 <- if (length(panels_d12) == 3)
      (p_d12_os | p_d12_ob) / p_d12_aic
    else
      wrap_plots(panels_d12, ncol = 2)
    pdf(file.path(cur_out_dir, "D12_delay_by_outcome.pdf"), width = 16, height = 12)
    print(layout_d12 + plot_annotation(
      title    = sprintf("Delay distributions by test outcome — %s", WINDOW_NAME),
      subtitle = paste("Censored MLE best-fit density per result category.",
                       "Grey bars = pooled histogram. Coloured lines = best-fit family per outcome.",
                       "Panel C: ΔAIC — lower = better relative fit within each outcome."),
      theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                       plot.subtitle = element_text(size = 9,  colour = "grey40"))))
    dev.off(); cat("  D12_delay_by_outcome.pdf\n")
  }

}, error = function(e)
  message("  [D12] delay-by-outcome failed: ", e$message))

# ============================================================================
# FINDINGS SUMMARY
# ============================================================================
cat("-- FINDINGS SUMMARY --\n")
lines_out <- c(
  strrep("=", 72),
  sprintf("FINDINGS SUMMARY — %s", WINDOW_NAME),
  sprintf("Window:   %s – %s  |  %s", ANALYSIS_START, TRUNC_DATE, Sys.Date()),
  sprintf("Script:   04c_delay_four_windows.R (v7 pipeline)"),
  strrep("=", 72),
  "",
  sprintf("Records in window:  %d  |  pos-either: %d  |  with onset: %d",
          nrow(df), nrow(pos_either), nrow(df[!is.na(df$symptom_onset), ])),
  "",
  "DELAY ESTIMATES (best-fit implied mean):",
  sprintf("  %-26s  %-10s %8s  %-10s %8s",
          "Delay", "Naive fam", "mean(d)", "Censored fam", "mean(d)")
)

for (dt in c("onset_sample", "onset_first_analysis",
             "onset_death", "receipt_first_analysis")) {
  nr <- .mle_tbl  %>% dplyr::filter(delay_type == dt, best) %>% dplyr::slice(1)
  cr <- .cens_tbl %>% dplyr::filter(delay_type == dt, best) %>% dplyr::slice(1)
  nf <- if (nrow(nr) > 0) nr$family[1] else "—"
  nm <- if (nrow(nr) > 0 && !is.na(nr$implied_mean_d[1]))
    sprintf("%.2f", nr$implied_mean_d[1]) else "—"
  cf <- if (nrow(cr) > 0) cr$family[1] else "—"
  cm <- if (nrow(cr) > 0 && !is.na(cr$implied_mean_d[1]))
    sprintf("%.2f", cr$implied_mean_d[1]) else "—"
  lines_out <- c(lines_out,
    sprintf("  %-26s  %-10s %8s  %-10s %8s", dt, nf, nm, cf, cm))
}

lines_out <- c(lines_out, "",
  "DETECTION COMPLETENESS:",
  sprintf("  Cases with onset:             %d", dc_res$n_onset),
  "  (1) NOT YET SAMPLED (onset→sample):",
  sprintf("    Empirical:   %d  (%.1f%%)",
          dc_res$n_unsampled_obs, 100 * dc_res$prop_unsampled_emp),
  sprintf("    Censored MLE (%s): %.1f exp  (%.1f%%)",
          toupper(dc_res$best_cens_fam_os),
          dc_res$exp_unsampled_cens, 100 * dc_res$prop_unsampled_cens),
  "  (2) NOT YET ANALYZED (onset→analysis):",
  sprintf("    Empirical:   %d  (%.1f%%)",
          dc_res$n_not_analyzed_obs, 100 * dc_res$prop_not_analyzed_emp),
  sprintf("    Censored MLE (%s): %.1f exp  (%.1f%%)",
          toupper(dc_res$best_cens_fam_ob %||% "N/A"),
          dc_res$exp_not_analyzed_cens, 100 * dc_res$prop_not_analyzed_cens),
  "",
  sprintf("CENSORED ONSET-TO-SAMPLE PARAMS (for growth-rate notebooks):"),
  sprintf("  Family: %s | mean=%.3fd",
          toupper(.cens_best_fam(cens_os)),
          .cens_best_mean(cens_os))
)

writeLines(lines_out, file.path(cur_out_dir, "findings_summary.txt"))
cat("Saved findings_summary.txt\n")
cat("\nCensored onset-to-sample params (for growth-rate notebooks):\n")
if (!is.null(cens_os) && length(cens_os) > 0) {
  bc  <- .cens_best_fam(cens_os)
  cat(sprintf("  Best family: %s | mean=%.3fd\n", toupper(bc),
              .cens_best_mean(cens_os)))
  cat("  Parameters:\n")
  print(round(cens_os[[bc]]$estimate, 4))
}
cat("\nDelay params CSV for growth-rate notebooks:\n")
cat(sprintf("  %s\n", delay_params_out))
cat(sprintf("\n%s\n", strrep("=", 80)))
cat(sprintf("04c_delay_four_windows.R COMPLETE — window: %s\n", WINDOW_NAME))
cat(sprintf("%s\n", strrep("=", 80)))
