# =============================================================================
# 45_ascertainment_figures.R
# BDBV 2026 DRC — Spatiotemporal invasion forecasting
#
# DETECTED INVASION vs TRUE INVASION: spatial variation in ascertainment.
#
# The modelled outcome is the onset week of a health zone's FIRST CONFIRMED CASE.
# That is a *detected* invasion, not a latent epidemiological one. This module
# quantifies, from the study's own data, how far the two can diverge and whether
# the divergence is spatially structured in a way that could bias either the
# invasion labels or the out-of-sample evaluation.
#
#   Figure A1  The ascertainment landscape
#                A) national map of the surveillance gap, invaded zones outlined
#                B) province case-confirmation ratio among laboratory-resolved
#                   alerts, with Wilson 95% intervals
#                C) zone confirmation ratio vs laboratory-resolved alerts per 100k
#                D) Spearman correlation matrix of the ascertainment proxies
#
#   Figure A2  Does thin surveillance delay the first confirmed case?
#                A) arrival-time residual after epicentre travel time vs surveillance gap
#                B) added-variable forest: per-SD slope of each proxy on arrival time,
#                   over the travel-time baseline, with and without the epicentre zones
#                C) arrival time vs the zone's median onset -> confirmation delay
#                D) arrival-time ECDF by surveillance-gap tercile
#
#   Figure A3  Are invaded and never-invaded zones comparable?
#                A) surveillance gap by invasion status (violin + box + points)
#                B) Cliff's delta forest, invaded vs never-invaded, every proxy
#                C) the same contrast WITHIN quartiles of epicentre travel time
#                D) logistic odds ratios per SD, unadjusted and adjusted for
#                   epicentre travel time
#
#   Figure A4  Does the forecast rank by surveillance?
#                A) leave-future-out rank vs surveillance gap (at-risk zone-folds)
#                B) rank of the truly-invaded zones vs their surveillance gap
#                C) within-stratum discrimination (AUC-PR skill, AUC-ROC) by
#                   surveillance stratum, both horizons, zone-cluster bootstrap
#                D) watch-list recall vs budget K by surveillance stratum
#
#   Figure A5  Province test positivity and the ranking
#                A) province confirmation ratio vs mean rank of its at-risk zone-folds
#                B) rank distributions by province, events vs non-events
#                C) province laboratory-resolved alerts per 100k vs mean rank
#                D) live 2-week invasion probability vs province confirmation ratio
#
# NOTHING IS REFITTED. Every quantity is derived from artifacts the pipeline has
# already written (the LFO-CV frame, the risk table, the arrival-predictor table,
# the covariate layers) plus the line list itself. The module writes figures and
# the CSV of every number plotted; it writes no pipeline artifact and overwrites
# no published output.
#
# Outputs -> outputs/key_outputs/ascertainment/{FigureA1..A5}.{pdf,png}
#            outputs/key_outputs/ascertainment/panels/
#            outputs/key_outputs/ascertainment/*.csv
# Run:  Rscript 45_ascertainment_figures.R   (from spatiotemporal/)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(sf); library(scales)
  library(forcats); library(purrr)
})
sf::sf_use_s2(TRUE)
options(dplyr.summarise.inform = FALSE)
for (pk in c("ggrepel", "here", "jsonlite")) {
  if (!requireNamespace(pk, quietly = TRUE)) stop("45_ascertainment_figures.R needs ", pk)
}

HERE <- file.path(here::here(), "spatiotemporal")
OUT  <- file.path(HERE, "outputs")

# 01_data_prep.R sources 00_config.R; 16_invasion_eval.R carries the SCORING
# functions (.auc_pr, .auc_roc) so the stratified skill in Figure A4 is computed
# by exactly the code that produced the headline numbers, not a second copy.
# 43_spread_kinematics.R is source-guarded (is_script_run) and contributes the
# geography helpers (.load_osrm_named, build_geo) that define the epicentre
# travel time used in Figure 1C, so the covariate here is the same covariate.
source(file.path(HERE, "01_data_prep.R"))
source(file.path(HERE, "16_invasion_eval.R"))
source(file.path(HERE, "43_spread_kinematics.R"))
source(file.path(HERE, "forecast_scale.R"))
stopifnot(exists("load_linelist"), exists("load_population"), exists("load_static_covariates"),
          exists(".auc_pr"), exists(".auc_roc"), exists("build_geo"),
          exists("EPICENTRE_ZONES"), exists("ANALYSIS_DATE"), exists("OUTBREAK_START"),
          exists("CONFIRMED_STATUS"), exists("DELAY_CAP_DAYS"), exists("RANDOM_SEED"))

FIG_DIR   <- file.path(OUT, "key_outputs", "ascertainment")
PANEL_DIR <- file.path(FIG_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

CONF_LEVEL <- 0.95          # every interval in this module, stated once
N_BOOT     <- 2000L         # bootstrap replicates (zone-cluster where applicable)
ALPHA      <- 1 - CONF_LEVEL
Q_LO <- ALPHA / 2; Q_HI <- 1 - ALPHA / 2

# -----------------------------------------------------------------------------
# 1. DESIGN SYSTEM — identical tokens to make_manuscript_figures.R (no titles)
# -----------------------------------------------------------------------------
INK <- "grey15"; MUTED <- "grey38"; FAINT <- "grey72"; GRID <- "grey92"
NA_FILL <- "grey93"
OKABE <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
PT_BLUE <- "#4C78C8"; FIT_RED <- "#C0392B"
INV_COL <- c("Invaded" = "#C0392B", "Never invaded" = "#4C78C8")
PROV_COL <- c("Ituri" = "#0072B2", "Nord-Kivu" = "#D55E00", "Haut-Uele" = "#009E73",
              "Tshopo" = "#CC79A7", "Bas-Uele" = "#E69F00", "Sud-Kivu" = "#56B4E9",
              "Tshuapa" = "#7A7A7A", "Kinshasa" = "#8C6D31", "Kasai" = "#6A51A3",
              "Other" = "#9E9E9E")
base_family <- "sans"

theme_pub <- function(base = 8.6) {
  theme_minimal(base_size = base, base_family = base_family) %+replace% theme(
    plot.title    = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    axis.title    = element_text(size = base - 0.4, colour = MUTED),
    axis.title.x  = element_text(margin = margin(t = 4)),
    axis.title.y  = element_text(margin = margin(r = 4), angle = 90),
    axis.text     = element_text(size = base - 1.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
    legend.position = "top", legend.justification = "left",
    legend.title  = element_text(size = base - 1.2, colour = MUTED),
    legend.text   = element_text(size = base - 1.4, colour = INK),
    legend.key.height = unit(9, "pt"), legend.key.width = unit(15, "pt"),
    strip.text    = element_text(size = base - 0.6, colour = INK, face = "bold",
                                 margin = margin(3, 3, 3, 3)),
    plot.tag      = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin   = margin(6, 8, 6, 6))
}
theme_map <- function(base = 8.6) {
  theme_void(base_size = base, base_family = base_family) %+replace% theme(
    plot.title    = element_blank(),
    legend.position = "right",
    legend.title  = element_text(size = base - 1.4, colour = MUTED),
    legend.text   = element_text(size = base - 1.8, colour = INK),
    legend.key.height = unit(20, "pt"), legend.key.width = unit(7, "pt"),
    plot.tag      = element_text(size = base + 4.5, face = "bold", colour = INK),
    plot.margin   = margin(2, 2, 2, 2))
}
#' In-panel annotation anchored to a panel CORNER via +/-Inf, so the text never has to be
#' positioned in data coordinates (which silently moves when the data change) and works
#' identically on continuous and discrete scales.
#' Drawn as a borderless label on a translucent white ground rather than bare text, so an
#' annotation can never be rendered illegible by points or labels that happen to sit under it.
ann_in <- function(label, x = -Inf, y = Inf, hjust = -0.04, vjust = 1.12, size = 2.45)
  annotate("label", x = x, y = y, label = label, hjust = hjust, vjust = vjust,
           size = size, colour = MUTED, family = base_family, lineheight = 1.08,
           fill = grDevices::adjustcolor("white", alpha.f = 0.82), linewidth = 0,
           label.r = unit(0, "pt"), label.padding = unit(1.6, "pt"))

save_dual <- function(p, name, w, h, dir = PANEL_DIR) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  message(sprintf("  saved %-42s %.1f x %.1f in", name, w, h)); invisible(p)
}

# -----------------------------------------------------------------------------
# 2. STATISTICAL HELPERS
#    Every estimator used anywhere below is defined once, here, and each returns
#    its own sample size so no panel can quietly plot a statistic computed on a
#    different support from the one its neighbours use.
# -----------------------------------------------------------------------------

#' Wilson score interval for a binomial proportion. Undefined (NA) at n = 0 —
#' never 0/0 silently rendered as zero.
wilson_ci <- function(x, n, conf = CONF_LEVEL) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  out <- tibble::tibble(x = as.numeric(x), n = as.numeric(n),
                        est = NA_real_, lo = NA_real_, hi = NA_real_)
  ok <- is.finite(out$n) & out$n > 0 & is.finite(out$x) & out$x >= 0 & out$x <= out$n
  if (any(ok)) {
    p <- out$x[ok] / out$n[ok]; nn <- out$n[ok]
    den <- 1 + z^2 / nn
    ctr <- (p + z^2 / (2 * nn)) / den
    hw  <- z * sqrt(p * (1 - p) / nn + z^2 / (4 * nn^2)) / den
    out$est[ok] <- p
    out$lo[ok]  <- pmax(0, ctr - hw)
    out$hi[ok]  <- pmin(1, ctr + hw)
  }
  out
}

#' Spearman rank correlation with a Bonett-Wright confidence interval and the
#' asymptotic t-approximation p-value (exact = FALSE, so ties are permitted).
#' Complete cases only; the n returned is the number of pairs actually used.
spearman_ci <- function(x, y, conf = CONF_LEVEL) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 4L || stats::sd(x) == 0 || stats::sd(y) == 0)
    return(tibble::tibble(rho = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, n = n))
  rho <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  p   <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  # Bonett & Wright (2000): Fisher z on the rank correlation with an inflated SE.
  se <- sqrt((1 + rho^2 / 2) / (n - 3))
  z  <- atanh(pmin(pmax(rho, -0.999999), 0.999999))
  zc <- stats::qnorm(1 - (1 - conf) / 2)
  tibble::tibble(rho = rho, lo = tanh(z - zc * se), hi = tanh(z + zc * se), p = p, n = n)
}

#' Cliff's delta = P(a > b) - P(a < b), computed exactly from mid-ranks (so ties
#' contribute zero), with a percentile bootstrap interval that resamples the two
#' groups independently — the zones in the two groups are distinct units.
cliffs_delta <- function(a, b, n_boot = N_BOOT, conf = CONF_LEVEL, seed = RANDOM_SEED) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  na <- length(a); nb <- length(b)
  if (na < 2L || nb < 2L)
    return(tibble::tibble(delta = NA_real_, lo = NA_real_, hi = NA_real_, n_a = na, n_b = nb))
  .d <- function(u, v) {
    r <- rank(c(u, v), ties.method = "average")            # mid-ranks: ties weigh 0.5 each way
    Ua <- sum(r[seq_along(u)]) - length(u) * (length(u) + 1) / 2   # Mann-Whitney U for `u`
    2 * (Ua / (length(u) * length(v))) - 1
  }
  d0 <- .d(a, b)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  bs <- vapply(seq_len(n_boot), function(i)
    .d(sample(a, na, replace = TRUE), sample(b, nb, replace = TRUE)), numeric(1))
  if (is.null(old)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
  else assign(".Random.seed", old, envir = globalenv())
  q <- stats::quantile(bs, c((1 - conf) / 2, 1 - (1 - conf) / 2), na.rm = TRUE, names = FALSE)
  tibble::tibble(delta = d0, lo = q[1], hi = q[2], n_a = na, n_b = nb)
}

#' Added-variable linear fit: what does `prox` add to a baseline regression of
#' `y` on `base`? Both models are fit on the SAME complete-case rows, so the
#' reported delta-R-squared is a like-for-like nested comparison rather than two
#' models fit on different supports.
#' Slope is reported PER SD of `prox` on those rows, so proxies in different
#' units are comparable on one axis.
av_slope <- function(y, base, prox, conf = CONF_LEVEL) {
  ok <- is.finite(y) & is.finite(base) & is.finite(prox)
  y <- y[ok]; base <- base[ok]; prox <- prox[ok]; n <- length(y)
  s <- stats::sd(prox)
  if (n < 10L || !is.finite(s) || s == 0)
    return(tibble::tibble(slope_sd = NA_real_, lo = NA_real_, hi = NA_real_,
                          p = NA_real_, dR2 = NA_real_, n = n))
  m0 <- stats::lm(y ~ base)
  m1 <- stats::lm(y ~ base + prox)
  cf <- summary(m1)$coefficients
  ci <- suppressMessages(stats::confint(m1, "prox", level = conf))
  tibble::tibble(slope_sd = unname(cf["prox", "Estimate"]) * s,
                 lo = unname(ci[1, 1]) * s, hi = unname(ci[1, 2]) * s,
                 p  = unname(cf["prox", "Pr(>|t|)"]),
                 dR2 = summary(m1)$r.squared - summary(m0)$r.squared,
                 n = n)
}

#' Logistic odds ratio per SD of `prox` for a binary outcome, either unadjusted
#' or adjusted for `base`. Profile-likelihood interval where it can be computed,
#' Wald otherwise; which one was used is returned, never hidden.
logit_or <- function(yb, prox, base = NULL, conf = CONF_LEVEL) {
  ok <- is.finite(prox) & !is.na(yb)
  if (!is.null(base)) ok <- ok & is.finite(base)
  yb <- as.integer(yb[ok]); prox <- prox[ok]
  base_v <- if (is.null(base)) NULL else base[ok]
  n <- length(yb); s <- stats::sd(prox)
  bad <- tibble::tibble(or_sd = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_,
                        n = n, n_pos = sum(yb), ci_type = NA_character_)
  if (n < 20L || !is.finite(s) || s == 0 || sum(yb) < 5L || sum(yb) == n) return(bad)
  dd <- if (is.null(base_v)) data.frame(yb = yb, prox = prox)
        else data.frame(yb = yb, prox = prox, base = base_v)
  fm <- if (is.null(base_v)) yb ~ prox else yb ~ base + prox
  m <- tryCatch(stats::glm(fm, data = dd, family = stats::binomial()), error = function(e) NULL)
  if (is.null(m) || !isTRUE(m$converged)) return(bad)
  cf <- summary(m)$coefficients
  ci <- tryCatch(suppressMessages(stats::confint(m, "prox", level = conf)), error = function(e) NULL)
  ci_type <- "profile"
  if (is.null(ci) || any(!is.finite(ci))) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    ci <- c(cf["prox", "Estimate"] - z * cf["prox", "Std. Error"],
            cf["prox", "Estimate"] + z * cf["prox", "Std. Error"])
    ci_type <- "wald"
  }
  tibble::tibble(or_sd = exp(unname(cf["prox", "Estimate"]) * s),
                 lo = exp(unname(ci[1]) * s), hi = exp(unname(ci[2]) * s),
                 p  = unname(cf["prox", "Pr(>|z|)"]), n = n, n_pos = sum(yb),
                 ci_type = ci_type)
}

#' Discrimination on ONE subset of at-risk zone-folds, with a zone-cluster
#' bootstrap interval. `.auc_pr` / `.auc_roc` come from 16_invasion_eval.R, so
#' these are the pipeline's own estimators, tie-handling included.
#'
#' THREE quantities are returned, and the reason is a trap worth naming. AUC-PR SKILL
#' is average precision divided by the subset's own base rate, so 1 is "no skill" in
#' every subset — but its ATTAINABLE MAXIMUM is 1 / base rate, which differs between
#' subsets whenever their event rates differ. Comparing raw skill across strata is
#' therefore not scale-free at the top end: a rarer stratum has more room above 1.
#' `auc_pr` (average precision itself, the ceiling-normalised quantity) and `auc_roc`
#' (a common [0, 1] scale with a common 0.5 null) are the comparable ones, and
#' `skill_ceiling` is reported so the skill column can be read correctly.
strat_discrimination <- function(d, n_boot = N_BOOT, conf = CONF_LEVEL, seed = RANDOM_SEED) {
  d <- d[is.finite(d$p_invasion) & !is.na(d$is_new_invasion), , drop = FALSE]
  y <- as.integer(d$is_new_invasion); p <- d$p_invasion
  n <- length(y); n_pos <- sum(y)
  # n_pos counts positive ROWS (the pipeline's convention, and the denominator of the
  # pooled recall); n_pos_zones counts the DISTINCT invaded zones behind them. At the
  # 2-week horizon the two differ because consecutive outcome windows overlap, and only
  # the second is an effective sample size.
  n_pos_zones <- dplyr::n_distinct(d$health_zone[y == 1])
  out <- tibble::tibble(n = n, n_pos = n_pos, n_pos_zones = n_pos_zones,
                        base_rate = if (n > 0) mean(y) else NA_real_,
                        auc_pr = NA_real_, ap_lo = NA_real_, ap_hi = NA_real_,
                        auc_pr_skill = NA_real_, skill_lo = NA_real_, skill_hi = NA_real_,
                        skill_ceiling = NA_real_,
                        auc_roc = NA_real_, roc_lo = NA_real_, roc_hi = NA_real_,
                        n_boot_used = 0L, n_zones = dplyr::n_distinct(d$health_zone))
  if (n == 0 || n_pos == 0) return(out)
  br <- mean(y)
  out$auc_pr <- .auc_pr(p, y)
  out$auc_pr_skill <- out$auc_pr / br
  out$skill_ceiling <- 1 / br          # average precision cannot exceed 1
  out$auc_roc <- .auc_roc(p, y)
  zones <- unique(d$health_zone)
  if (n_pos >= 2L && length(zones) > 5L) {
    idx_by_zone <- split(seq_len(nrow(d)), factor(d$health_zone, levels = zones))
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(seed)
    bs <- matrix(NA_real_, nrow = n_boot, ncol = 3)
    for (b in seq_len(n_boot)) {
      zs  <- sample.int(length(zones), length(zones), replace = TRUE)
      idx <- unlist(idx_by_zone[zs], use.names = FALSE)
      yb <- y[idx]; pb <- p[idx]
      if (sum(yb) >= 1L && sum(yb) < length(yb)) {
        ap <- .auc_pr(pb, yb)
        # SKILL is scored against the REPLICATE's own base rate, exactly as the point
        # estimate is scored against the observed one; average precision is recorded
        # separately, because skill quantiles times the observed base rate are NOT the
        # average-precision quantiles (the divisor varies replicate to replicate).
        bs[b, 1] <- ap / mean(yb)
        bs[b, 2] <- .auc_roc(pb, yb)
        bs[b, 3] <- ap
      }
    }
    if (is.null(old)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
    else assign(".Random.seed", old, envir = globalenv())
    used <- sum(is.finite(bs[, 1]))
    if (used >= 50L) {
      .q <- function(v) stats::quantile(v, c((1 - conf) / 2, 1 - (1 - conf) / 2),
                                        na.rm = TRUE, names = FALSE)
      qs <- .q(bs[, 1]); qr <- .q(bs[, 2]); qa <- .q(bs[, 3])
      out$skill_lo <- qs[1]; out$skill_hi <- qs[2]
      out$roc_lo   <- qr[1]; out$roc_hi   <- qr[2]
      out$ap_lo    <- qa[1]; out$ap_hi    <- qa[2]
    }
    out$n_boot_used <- as.integer(used)
  }
  out
}

#' Pooled top-K recall with a ZONE-CLUSTER bootstrap interval.
#'
#' The point estimate is the pipeline's pooled estimator (hits / invasion rows, which
#' compute_detection_curve reports as `recall_pooled`). The INTERVAL may not be a Wilson
#' interval on those rows: at the 2-week horizon the outcome windows of consecutive
#' forecast origins overlap, so one invaded zone contributes up to two positive rows and
#' the rows are not independent Bernoulli trials. Resampling whole ZONES is the same
#' correlated unit evaluate_invasion() bootstraps over.
recall_curve_boot <- function(ev_rows, ks, n_boot = N_BOOT, conf = CONF_LEVEL, seed = RANDOM_SEED) {
  rk <- ev_rows$rank_max
  zones <- unique(ev_rows$health_zone)
  out <- tibble::tibble(k = ks,
                        caught = vapply(ks, function(k) sum(rk <= k), numeric(1)),
                        events = nrow(ev_rows), n_zones = length(zones),
                        recall = vapply(ks, function(k) mean(rk <= k), numeric(1)),
                        lo = NA_real_, hi = NA_real_)
  if (length(zones) < 3L || nrow(ev_rows) < 3L) return(out)
  idx_by_zone <- split(seq_len(nrow(ev_rows)), factor(ev_rows$health_zone, levels = zones))
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  bs <- matrix(NA_real_, nrow = n_boot, ncol = length(ks))
  for (b in seq_len(n_boot)) {
    idx <- unlist(idx_by_zone[sample.int(length(zones), length(zones), replace = TRUE)],
                  use.names = FALSE)
    bs[b, ] <- colMeans(outer(rk[idx], ks, "<="))
  }
  if (is.null(old)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
  else assign(".Random.seed", old, envir = globalenv())
  out$lo <- apply(bs, 2L, stats::quantile, probs = (1 - conf) / 2, na.rm = TRUE, names = FALSE)
  out$hi <- apply(bs, 2L, stats::quantile, probs = 1 - (1 - conf) / 2, na.rm = TRUE, names = FALSE)
  out
}

#' Expand a table of (x, y) into the vertex sequence of a `geom_step(direction = "hv")`
#' path, so a confidence RIBBON can be drawn along the same staircase as the line it
#' bounds. Drawing the ribbon with plain `geom_ribbon` interpolates diagonally between
#' consecutive x values, i.e. it draws uncertainty at watch-list sizes that do not exist.
step_path <- function(d, xcol = "k") {
  d <- d[order(d[[xcol]]), , drop = FALSE]
  n <- nrow(d)
  if (n < 2L) return(d)
  x <- d[[xcol]]
  out <- d[rep(seq_len(n), each = 2L), , drop = FALSE]
  # Row pair i holds y_i from x_i to x_{i+1}. The construction leaves a duplicate of the
  # final vertex, which is dropped so the result is the 2n-1 vertices ggplot2's own
  # stairstep() produces rather than 2n with a degenerate last segment.
  out[[xcol]] <- as.vector(rbind(x, c(x[-1], x[n])))
  out[-nrow(out), , drop = FALSE]
}

#' Format a p-value for a panel annotation without ever printing "p = 0".
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
fmt_ci <- function(e, l, h, d = 2) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), e, l, h)

# -----------------------------------------------------------------------------
# 3. DATA ASSEMBLY
# -----------------------------------------------------------------------------
message("\n", strrep("=", 74))
message("[ascert] 45_ascertainment_figures.R — detection vs onset, spatial ascertainment")
message(strrep("=", 74))

# ---- 3.1 the modelling spine, the invasion labels, and the vulnerability pillars ----
rs_all <- readr::read_csv(fs_risk_csv(OUT), show_col_types = FALSE)
stopifnot(all(c("health_zone", "province", "horizon", "was_active_before",
                "p_case_invasion", "surveillance_gap", "healthcare_gap",
                "social_vulnerability", "access_gap", "healthcare_travel_min") %in% names(rs_all)))
RISK_METHOD <- unique(rs_all$method)
stopifnot(length(RISK_METHOD) == 1L)
PROB_SCALE  <- unique(rs_all$prob_scale)
stopifnot(length(PROB_SCALE) == 1L)   # it is printed on a panel; two scales in one table
                                      # would be a caption naming the wrong one

spine <- rs_all %>%
  dplyr::filter(horizon == 1L) %>%
  dplyr::transmute(health_zone, province,
                   invaded = as.logical(was_active_before),
                   surveillance_gap, healthcare_gap, social_vulnerability,
                   access_gap, healthcare_travel_min)
stopifnot(nrow(spine) == dplyr::n_distinct(spine$health_zone), !anyNA(spine$invaded))
message(sprintf("[ascert] spine: %d zones, %d invaded, %d never invaded (risk method %s, scale %s)",
                nrow(spine), sum(spine$invaded), sum(!spine$invaded), RISK_METHOD, PROB_SCALE))

# The two access pillars are DEGENERATE on this covariate set: every one of the 519
# zones hosts at least one health facility, so the "own facility -> 0 minutes" rule in
# compute_vulnerability_index() zeroes travel time and hence access_gap for all of them
# (20_forecast_detail.R documents this and drops the pillar from V for the same reason).
# A constant cannot be an ascertainment axis, so it is dropped here too — but by TEST,
# not by assumption, so a future run in which the pillar becomes informative keeps it.
.degenerate <- function(v) dplyr::n_distinct(v[is.finite(v)]) <= 1L
DROPPED_PILLARS <- c("access_gap", "healthcare_travel_min")[
  c(.degenerate(spine$access_gap), .degenerate(spine$healthcare_travel_min))]
if (length(DROPPED_PILLARS))
  message("[ascert] degenerate pillar(s) dropped from the proxy set: ",
          paste(DROPPED_PILLARS, collapse = ", "))

# ---- 3.2 raw covariate layers (facility density, population, deprivation) -----------
cov <- load_static_covariates() %>%
  dplyr::transmute(health_zone = nom,
                   pop_count = as.numeric(pop_count),
                   healthsite_count = as.numeric(healthsite_count),
                   healthsite_density = as.numeric(healthsite_density),
                   ccvi = as.numeric(ccvi))
stopifnot(setequal(cov$health_zone, spine$health_zone))

# ---- 3.3 travel time from the epicentre region, for ALL zones -----------------------
# arrival_predictors.csv carries this only for the INVADED zones; the invaded-vs-not
# comparison needs it for the whole spine. It is recomputed here from the same OSRM matrix
# and the same min-over-EPICENTRE_ZONES rule that 43_spread_kinematics.R uses, then
# ASSERTED equal to the published column on every zone where both exist — so the new
# values come from a rule checked against the published one, not assumed to match it.
geo <- build_geo()
time_h <- geo$time                                   # hours, symmetrised
epi_set <- intersect(EPICENTRE_ZONES, rownames(time_h))
stopifnot(length(epi_set) == length(EPICENTRE_ZONES))
.epi_h <- apply(time_h[epi_set, , drop = FALSE], 2L,
                function(z) { z <- z[is.finite(z)]; if (length(z)) min(z) else NA_real_ })
epi_travel <- tibble::tibble(health_zone = names(.epi_h), epi_travel_h = as.numeric(.epi_h))
stopifnot(setequal(epi_travel$health_zone, spine$health_zone))

arr <- readr::read_csv(file.path(OUT, "key_outputs", "arrival_predictors.csv"),
                       show_col_types = FALSE) %>%
  dplyr::select(health_zone, first_onset, arrival_days, arrival_origin,
                travel_time_h, dist_greatcircle_km, neglog_share_shorttrip,
                log10_population, cases_total)
.chk <- dplyr::inner_join(arr, epi_travel, by = "health_zone")
stopifnot(nrow(.chk) == nrow(arr))
.dev <- max(abs(.chk$travel_time_h - .chk$epi_travel_h), na.rm = TRUE)
if (!is.finite(.dev) || .dev > 1e-8)
  stop(sprintf(paste0("[ascert] recomputed epicentre travel time disagrees with ",
                      "arrival_predictors.csv by up to %.6g h — the two are NOT the same ",
                      "covariate and nothing downstream may treat them as one."), .dev))
message(sprintf("[ascert] epicentre travel time reproduced for all %d zones (max deviation %.2g h on the %d published)",
                nrow(epi_travel), .dev, nrow(arr)))

# ---- 3.4 line-list ascertainment measures -------------------------------------------
ll <- load_linelist()
# The sitrep top-up rows are a COUNT reconciliation against the official cumulative, not
# laboratory records: they carry a classification but no test, no lab date and no alert.
# Including them would inflate the confirmation ratio in exactly the zones that received a
# top-up — i.e. it would manufacture the spatial pattern being measured.
n_synth <- sum(grepl("^SITREP-CONF-", ll$alert_id))
ll_real <- ll %>% dplyr::filter(is.na(alert_id) | !grepl("^SITREP-CONF-", alert_id))
stopifnot(nrow(ll) - nrow(ll_real) == n_synth)
message(sprintf("[ascert] line list: %d rows, %d synthetic sitrep top-up rows excluded, %d retained",
                nrow(ll), n_synth, nrow(ll_real)))

# Alerts are counted on the SAME window the zone-week grid uses: on or after
# OUTBREAK_START (load_linelist already enforces this) and on or before ANALYSIS_DATE.
ll_win <- ll_real %>%
  dplyr::filter(!is.na(health_zone), !is.na(date_index),
                date_index >= OUTBREAK_START, date_index <= ANALYSIS_DATE)
message(sprintf("[ascert] alerts inside [%s, %s]: %d (dropped %d outside the window or undated)",
                format(OUTBREAK_START), format(ANALYSIS_DATE), nrow(ll_win), nrow(ll_real) - nrow(ll_win)))

zone_alerts <- ll_win %>%
  dplyr::group_by(health_zone) %>%
  dplyr::summarise(
    alerts_total    = dplyr::n(),
    n_confirmed     = sum(final_mve_case_classification %in% CONFIRMED_STATUS),
    n_not_a_case    = sum(final_mve_case_classification %in% NOT_A_CASE_STATUS),
    .groups = "drop") %>%
  dplyr::mutate(alerts_resolved = n_confirmed + n_not_a_case)

# The confirmation ratio is an ALERT-level positivity: the share of laboratory-RESOLVED
# alerts that were confirmed. It is not a test-level positivity — a person may be swabbed
# more than once and DHIS2 carries no per-record test count (samples_analyzed is NA on
# every real record, which is why the pipeline's own `positivity` covariate is NA
# throughout). Alerts still under investigation (suspected/probable) and alerts never
# classified are in neither the numerator nor the denominator.
.w <- wilson_ci(zone_alerts$n_confirmed, zone_alerts$alerts_resolved)
zone_alerts <- zone_alerts %>%
  dplyr::mutate(confirm_ratio = .w$est, confirm_lo = .w$lo, confirm_hi = .w$hi,
                resolution_frac = dplyr::if_else(alerts_total > 0,
                                                 alerts_resolved / alerts_total, NA_real_))

# Onset -> laboratory confirmation delay, per zone. Three restrictions, each needed:
#   (a) NATIVE onsets only. load_linelist() imputes a large minority of onsets backwards
#       from the sample date by a draw from the fitted onset->sample delay, so a delay
#       computed on those records is a function of that same distribution and would be
#       circular. The `onset_imputed` flag identifies them.
#   (b) EQUAL FOLLOW-UP. Only cases whose onset is at least DELAY_CAP_DAYS before the
#       analysis date are used, so every zone's cases have had the same opportunity to be
#       confirmed; without this, late-invaded zones look faster purely by truncation.
#   (c) Delays are capped at DELAY_CAP_DAYS, matching week_completeness(). A zone in which
#       more than half of the eligible cases exceed the cap has a censored median, flagged
#       rather than silently reported as exactly the cap.
.conf_date <- dplyr::coalesce(as.Date(ll_win$lab_analysis_date),
                              as.Date(ll_win$reporting_date),
                              as.Date(ll_win$date_of_notification))
.delay_d <- as.numeric(.conf_date - ll_win$date_index)
.delay_ok <- ll_win$confirmed %in% TRUE &
             !(ll_win$onset_imputed %in% TRUE) &
             is.finite(.delay_d) & .delay_d >= 0 &
             ll_win$date_index <= (ANALYSIS_DATE - DELAY_CAP_DAYS)
zone_delay <- tibble::tibble(health_zone = ll_win$health_zone,
                             delay = pmin(.delay_d, DELAY_CAP_DAYS))[.delay_ok, ] %>%
  dplyr::group_by(health_zone) %>%
  dplyr::summarise(n_delay = dplyr::n(),
                   delay_median = stats::median(delay),
                   delay_p75 = stats::quantile(delay, 0.75, names = FALSE),
                   delay_over_cap = mean(delay >= DELAY_CAP_DAYS),
                   .groups = "drop") %>%
  dplyr::mutate(delay_censored = delay_over_cap >= 0.5) %>%
  dplyr::filter(n_delay >= 5L)
message(sprintf("[ascert] onset->confirmation delay estimable in %d zone(s) (>=5 native-onset confirmed cases with >=%d d follow-up); %d have a censored median",
                nrow(zone_delay), DELAY_CAP_DAYS, sum(zone_delay$delay_censored)))

# ---- 3.5 line-list / sitrep coverage ratio (the careseeking module's measure) --------
# Keyed on a lower-cased zone name; the join is asserted rather than assumed, because a
# silent mismatch would drop precisely the zones with the most extreme ratios.
cover_path <- file.path(OUT, "key_outputs", "coverage_ratio_zone.csv")
zone_cover <- if (file.exists(cover_path)) {
  cr <- readr::read_csv(cover_path, show_col_types = FALSE)
  key <- tibble::tibble(health_zone = spine$health_zone, .k = tolower(trimws(spine$health_zone)))
  crj <- cr %>% dplyr::mutate(.k = tolower(trimws(zone_n))) %>%
    dplyr::inner_join(key, by = ".k") %>%
    dplyr::transmute(health_zone, coverage_ratio = as.numeric(coverage_ratio),
                     coverage_usable = as.logical(usable))
  miss <- setdiff(tolower(trimws(cr$zone_n[cr$usable %in% TRUE])), key$.k)
  message(sprintf("[ascert] coverage ratio: %d/%d rows joined to the spine; %d usable zone(s) unmatched%s",
                  nrow(crj), nrow(cr), length(miss),
                  if (length(miss)) paste0(": ", paste(miss, collapse = ", ")) else ""))
  if (length(miss) > 0L)
    warning(sprintf("[ascert] %d usable coverage-ratio zone(s) did not match the spine and are absent from the coverage panels.",
                    length(miss)), call. = FALSE)
  # One row per zone, or the left-join below would ROW-MULTIPLY the 519-zone spine and
  # every subsequent count would be silently wrong.
  stopifnot(anyDuplicated(crj$health_zone) == 0L)
  crj
} else {
  message("[ascert] coverage_ratio_zone.csv absent — the coverage-ratio proxy is skipped.")
  tibble::tibble(health_zone = character(), coverage_ratio = numeric(), coverage_usable = logical())
}

# ---- 3.6 the assembled zone table ----------------------------------------------------
zone_ascert <- spine %>%
  dplyr::left_join(cov,         by = "health_zone") %>%
  dplyr::left_join(epi_travel,  by = "health_zone") %>%
  dplyr::left_join(zone_alerts, by = "health_zone") %>%
  dplyr::left_join(zone_delay,  by = "health_zone") %>%
  dplyr::left_join(zone_cover,  by = "health_zone") %>%
  dplyr::left_join(arr,         by = "health_zone") %>%
  dplyr::mutate(
    # Zones with no alert at all have no alert counts to be missing — they are zero, and
    # a zero alert count is an observation about surveillance, not a gap in the data.
    # A confirmation RATIO, by contrast, stays NA there: 0/0 is undefined, not 0.
    alerts_total     = dplyr::coalesce(alerts_total, 0L),
    n_confirmed      = dplyr::coalesce(n_confirmed, 0L),
    n_not_a_case     = dplyr::coalesce(n_not_a_case, 0L),
    alerts_resolved  = dplyr::coalesce(alerts_resolved, 0L),
    resolved_per_100k = dplyr::if_else(is.finite(pop_count) & pop_count > 0,
                                       alerts_resolved / pop_count * 1e5, NA_real_),
    log_pop          = log10(pmax(pop_count, 1)),
    status           = factor(dplyr::if_else(invaded, "Invaded", "Never invaded"),
                              levels = c("Never invaded", "Invaded")))
# STABLE ORDER, set once and inherited by every artifact derived from this table. Without
# it the row order came from whatever order `bayes_risk_scores_all_zones.csv` happened to
# be written in, so a re-run of the pipeline reshuffled every published CSV and made the
# two versions impossible to diff even where no value had changed.
zone_ascert <- dplyr::arrange(zone_ascert, health_zone)
stopifnot(nrow(zone_ascert) == nrow(spine),
          dplyr::n_distinct(zone_ascert$health_zone) == nrow(spine),
          !is.unsorted(zone_ascert$health_zone))

# ---- 3.7 SURVEILLANCE EFFORT, measured so it is not circular ------------------------
# "Did this zone ever have a laboratory-resolved alert?" is NOT a usable covariate for a
# comparison whose outcome is "did this zone record a confirmed case?": a confirmed case
# IS a resolved alert, so every invaded zone has one by construction. The DISCARDED alerts
# (`not_a_case`) carry no such mechanical link — they are investigations that found
# nothing — so they are what the surveillance-effort measures below count. Where the
# quantity is needed AT A FORECAST ORIGIN it is recomputed from the alerts dated on or
# before that origin, never from the end-of-outbreak total.
# Each discarded alert is dated by its NOTIFICATION date — when it entered the system and
# an investigation began — falling back to the case's own index date if that is missing.
# The laboratory result came later, so this is "an alert was being investigated by then",
# not "an alert had been resolved by then"; the panels that use it say so.
neg_alerts_long <- ll_win %>%
  dplyr::filter(final_mve_case_classification %in% NOT_A_CASE_STATUS) %>%
  dplyr::transmute(health_zone,
                   alert_date = dplyr::coalesce(as.Date(date_of_notification), date_index))
stopifnot(!anyNA(neg_alerts_long$alert_date))
message(sprintf("[ascert] discarded (not-a-case) alerts available for the effort measure: %d across %d zone(s)",
                nrow(neg_alerts_long), dplyr::n_distinct(neg_alerts_long$health_zone)))

zone_ascert <- zone_ascert %>%
  dplyr::mutate(neg_alerts = n_not_a_case,
                neg_per_100k = dplyr::if_else(is.finite(pop_count) & pop_count > 0,
                                              n_not_a_case / pop_count * 1e5, NA_real_),
                ever_investigated = alerts_resolved > 0L)

.n_na_travel <- sum(!is.finite(zone_ascert$epi_travel_h))
if (.n_na_travel > 0L)
  message(sprintf("[ascert] %d zone(s) have no road route to the epicentre region and are absent from every travel-time-conditioned panel: %s",
                  .n_na_travel, paste(zone_ascert$health_zone[!is.finite(zone_ascert$epi_travel_h)], collapse = ", ")))

# The proxy panel, split by whether a quantity can legitimately enter an INVADED-vs-NOT
# comparison. The exogenous set is fixed by the health system and geography and is
# measurable in a zone that never reported anything; the outbreak-conditioned set is
# partly generated BY the outbreak and is used only where every zone compared has
# already been invaded (Figure A2) or where the comparison is across forecast origins
# with the quantity re-measured at each (Figure A5).
PROXY_EXO <- tibble::tribble(
  ~var,                   ~label,                                          ~worse_high,
  "surveillance_gap",     "Surveillance gap\n(facility density pctile)",    TRUE,
  "healthcare_gap",       "Health-system gap\n(facilities per capita pctile)", TRUE,
  "social_vulnerability", "Deprivation\n(CCVI percentile)",                 TRUE,
  "log_hs_density",       "Health-facility density\n(log10 per km2)",       FALSE,
  "log_pop",              "Population\n(log10)",                            FALSE)
# The access pillars join the exogenous set IF the degeneracy test above clears them, so
# the test governs the analysis rather than only a log line. On the current covariate set
# it does not clear them and they are absent from every panel.
.access_kept <- setdiff(c("access_gap", "healthcare_travel_min"), DROPPED_PILLARS)
if (length(.access_kept))
  PROXY_EXO <- dplyr::bind_rows(PROXY_EXO, tibble::tibble(
    var = .access_kept,
    label = c(access_gap = "Healthcare-access gap\n(travel-time percentile)",
              healthcare_travel_min = "Travel time to care\n(minutes)")[.access_kept],
    worse_high = TRUE))
PROXY_ENDO <- tibble::tribble(
  ~var,              ~label,                                        ~worse_high,
  "confirm_ratio",   "Confirmation ratio\n(confirmed / resolved alerts)", TRUE,
  "log_neg_100k",    "Discarded alerts per 100k\n(log10)",           FALSE,
  "delay_median",    "Onset to confirmation\n(median days)",         TRUE,
  "coverage_ratio",  "Line-list / sitrep\ncoverage ratio",           FALSE)

zone_ascert <- zone_ascert %>%
  dplyr::mutate(
    log_hs_density = dplyr::if_else(is.finite(healthsite_density) & healthsite_density > 0,
                                    log10(healthsite_density), NA_real_),
    # log10 of a rate that is legitimately zero: the zero is kept as its own state by the
    # `ever_investigated` flag and the rate is NA here, rather than being shifted by an
    # arbitrary pseudo-count that would place "no surveillance at all" on the same
    # continuum as "a little surveillance".
    log_neg_100k = dplyr::if_else(is.finite(neg_per_100k) & neg_per_100k > 0,
                                  log10(neg_per_100k), NA_real_),
    coverage_ratio = dplyr::if_else(coverage_usable %in% TRUE, coverage_ratio, NA_real_))
PROXY_ALL <- dplyr::bind_rows(dplyr::mutate(PROXY_EXO, kind = "exogenous"),
                              dplyr::mutate(PROXY_ENDO, kind = "outbreak-conditioned"))
# Labels are the row key of every forest below, so a duplicate would silently merge two
# proxies into one axis position; vars must exist or a panel would plot NA.
stopifnot(all(PROXY_ALL$var %in% names(zone_ascert)),
          anyDuplicated(PROXY_ALL$var) == 0L,
          anyDuplicated(gsub("\n", " ", PROXY_ALL$label)) == 0L)

# -----------------------------------------------------------------------------
# 4. FORECAST FRAMES
# -----------------------------------------------------------------------------
FEATURED <- local({
  f <- file.path(OUT, "key_outputs", "model_selection.json")
  stopifnot(file.exists(f))
  sel <- jsonlite::fromJSON(f, simplifyVector = TRUE)
  m <- tryCatch(sel$featured$bayesian$method, error = function(e) NULL)
  if (is.null(m) || !length(m) || is.na(m[1])) m <- tryCatch(sel$featured$headline$method, error = function(e) NULL)
  stopifnot(!is.null(m), length(m) >= 1L, !is.na(m[1]))
  as.character(m[1])
})
message("[ascert] featured model (read from model_selection.json): ", FEATURED)

lfo_all <- readRDS(file.path(OUT, "forecasts", "lfo_results.rds"))
stopifnot(FEATURED %in% lfo_all$method)
# SUBSET TO THE FEATURED MODEL BEFORE APPLYING THE SCALE. fs_lfo_col() inspects every
# scored row in the frame it is handed and falls back to the raw column if ANY of them
# lacks a recalibrated probability. Two rank-only baselines (Adjacency-B7, Distance-B1)
# deliberately carry no `p_recal` — recalibrating an ordering is meaningless — so handing
# it the whole multi-method frame silently demotes the featured model to the raw scale too.
# The featured model's `p_recal` is complete, so filtering first is what actually puts
# this module on the scale its outputs claim. (Every RANK in this module is taken within
# a fold and is therefore invariant to the choice either way; the pooled discrimination
# in Figure A4C is not, which is why it matters.)
lfo <- dplyr::filter(lfo_all, method == FEATURED)
PCOL <- fs_apply_lfo_scale(lfo)     # puts the selected scale into `p_invasion`
message("[ascert] LFO probability column in use: ", PCOL)

# At-risk rows of the featured model, exactly the support evaluate_invasion() scores:
# not already affected, finite probability. Ranks are taken WITHIN (horizon, fold) —
# the whole watch-list of that round — with the pipeline's two tie conventions:
#   rank_avg  ties averaged, matching .ranking_metrics()'s mean rank of truth;
#   rank_max  ties to the worst position, matching compute_detection_curve()'s top-K,
#             so a zone counts as monitored only if monitoring K zones must include it.
# Both are monotone in p within a fold, so they are invariant to the probability scale.
lfo_f <- lfo %>%
  dplyr::filter(!(as.logical(was_active_before) %in% TRUE), is.finite(p_invasion)) %>%
  dplyr::group_by(horizon, fold_id) %>%
  dplyr::mutate(rank_avg = rank(-p_invasion, ties.method = "average"),
                rank_max = rank(-p_invasion, ties.method = "max"),
                n_atrisk_fold = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(zone_ascert, health_zone, province, surveillance_gap,
                                 healthcare_gap, social_vulnerability, log_hs_density,
                                 epi_travel_h, ever_investigated, neg_alerts),
                   by = "health_zone")
stopifnot(nrow(lfo_f) > 0, !anyNA(lfo_f$surveillance_gap), !anyNA(lfo_f$rank_avg))

# Prior surveillance effort AT EACH FORECAST ORIGIN: discarded alerts notified on or
# before that fold's cutoff. Built once per distinct cutoff (a dozen of them) so that
# "before the origin" has exactly ONE definition across every panel that uses it.
.cut <- sort(unique(lfo_f$cutoff))
.neg_by_cut <- purrr::map_dfr(.cut, function(cu) {
  neg_alerts_long %>%
    dplyr::filter(alert_date <= cu) %>%
    dplyr::count(health_zone, name = "neg_prior") %>%
    dplyr::mutate(cutoff = cu)
})
lfo_f <- lfo_f %>%
  dplyr::left_join(.neg_by_cut, by = c("health_zone", "cutoff")) %>%
  dplyr::mutate(neg_prior = dplyr::coalesce(neg_prior, 0L),
                watched_prior = neg_prior > 0L)
message(sprintf("[ascert] LFO at-risk zone-folds: %d (%d zones, %d folds); with a prior discarded alert: %d (%.1f%%)",
                nrow(lfo_f), dplyr::n_distinct(lfo_f$health_zone),
                dplyr::n_distinct(lfo_f$fold_id), sum(lfo_f$watched_prior),
                100 * mean(lfo_f$watched_prior)))

# Surveillance stratum: a zone-level split at the MEDIAN surveillance gap of the zones the
# LFO actually scores, so a zone belongs to one stratum in every fold and the zone-cluster
# bootstrap in strat_discrimination() resamples whole zones within a stratum.
.sg_zone <- lfo_f %>% dplyr::distinct(health_zone, surveillance_gap)
SG_SPLIT <- stats::median(.sg_zone$surveillance_gap)
lfo_f <- lfo_f %>%
  dplyr::mutate(sg_stratum = factor(dplyr::if_else(surveillance_gap <= SG_SPLIT,
                                                   "Denser facilities", "Sparser facilities"),
                                    levels = c("Denser facilities", "Sparser facilities")))
message(sprintf("[ascert] surveillance stratum split at gap = %.4f (%d vs %d zones)",
                SG_SPLIT, sum(.sg_zone$surveillance_gap <= SG_SPLIT),
                sum(.sg_zone$surveillance_gap > SG_SPLIT)))

# Live (current) forecast, at-risk zones only, for the operational panels.
live <- rs_all %>%
  dplyr::filter(!(as.logical(was_active_before) %in% TRUE), is.finite(p_case_invasion)) %>%
  dplyr::select(health_zone, province, horizon, p_case_invasion) %>%
  dplyr::left_join(dplyr::select(zone_ascert, health_zone, surveillance_gap,
                                 ever_investigated, neg_alerts, epi_travel_h),
                   by = "health_zone")

# -----------------------------------------------------------------------------
# 5. FIGURE A1 — the ascertainment landscape
# -----------------------------------------------------------------------------
message("\n[ascert] Figure A1 — ascertainment landscape")

shp <- geo$shp
stopifnot(setequal(as.character(shp$Nom), zone_ascert$health_zone))
shp_a <- shp %>%
  dplyr::mutate(health_zone = as.character(Nom)) %>%
  dplyr::left_join(zone_ascert, by = "health_zone")
prov_outline <- shp_a %>% dplyr::group_by(.prov = as.character(PROVINCE)) %>%
  dplyr::summarise(.groups = "drop") %>% sf::st_geometry()
inv_outline <- shp_a %>% dplyr::filter(invaded) %>% sf::st_geometry()

pA1a <- ggplot() +
  geom_sf(data = shp_a, aes(fill = surveillance_gap), colour = "white", linewidth = 0.06) +
  geom_sf(data = prov_outline, fill = NA, colour = "grey45", linewidth = 0.18) +
  geom_sf(data = inv_outline, fill = NA, colour = "white", linewidth = 0.46) +
  geom_sf(data = inv_outline, fill = NA, colour = "#C0392B", linewidth = 0.24) +
  scale_fill_viridis_c(option = "G", direction = -1, limits = c(0, 1),
                       breaks = c(0, 0.5, 1), labels = c("0\ndensest", "0.5", "1\nsparsest"),
                       name = "Surveillance gap") +
  theme_map()

# THREE distinct states, drawn as three layers rather than forced onto one continuous
# scale. A log fill cannot represent a zero, and flooring the zeros at an arbitrary small
# value (the first version of this panel used pmax(rate, 0.05)) would draw three zones
# that investigated NO alert as though they had investigated a few. The two zero states
# are also not the same state: "no laboratory-resolved alert at all" and "resolved alerts,
# all of them confirmed cases" say different things about the surveillance system.
ZERO_NEG_FILL <- "#E7DCC1"
.no_lab   <- dplyr::filter(shp_a, !ever_investigated)
.zero_neg <- dplyr::filter(shp_a, ever_investigated, neg_alerts == 0)
.some_neg <- dplyr::filter(shp_a, ever_investigated, neg_alerts > 0)
stopifnot(nrow(.no_lab) + nrow(.zero_neg) + nrow(.some_neg) == nrow(shp_a))
pA1b <- ggplot() +
  geom_sf(data = .no_lab,   fill = NA_FILL,        colour = "white", linewidth = 0.06) +
  geom_sf(data = .zero_neg, fill = ZERO_NEG_FILL,  colour = "white", linewidth = 0.06) +
  geom_sf(data = .some_neg, aes(fill = neg_per_100k), colour = "white", linewidth = 0.06) +
  geom_sf(data = prov_outline, fill = NA, colour = "grey45", linewidth = 0.18) +
  geom_sf(data = inv_outline, fill = NA, colour = "white", linewidth = 0.46) +
  geom_sf(data = inv_outline, fill = NA, colour = "#C0392B", linewidth = 0.24) +
  scale_fill_viridis_c(option = "F", direction = -1, transform = "log10",
                       labels = scales::label_number(accuracy = 0.1),
                       name = "Discarded alerts\nper 100k") +
  ann_in(sprintf("grey: no laboratory-resolved alert (%d zones, %d of the %d never invaded)\ncream: resolved alerts, none discarded (%d zones)",
                 nrow(.no_lab),
                 sum(!zone_ascert$ever_investigated & !zone_ascert$invaded),
                 sum(!zone_ascert$invaded), nrow(.zero_neg)),
         x = -Inf, y = -Inf, hjust = -0.02, vjust = -0.15, size = 2.2) +
  theme_map()

prov_pos <- zone_ascert %>%
  dplyr::group_by(province) %>%
  dplyr::summarise(n_zones = dplyr::n(), n_invaded = sum(invaded),
                   confirmed = sum(n_confirmed), resolved = sum(alerts_resolved),
                   .groups = "drop")
.pw <- wilson_ci(prov_pos$confirmed, prov_pos$resolved)
prov_pos <- prov_pos %>% dplyr::mutate(ratio = .pw$est, lo = .pw$lo, hi = .pw$hi)
N_PROV_SILENT <- sum(prov_pos$resolved == 0)
prov_active <- prov_pos %>% dplyr::filter(resolved > 0) %>%
  dplyr::mutate(province = forcats::fct_reorder(province, ratio))

pA1c <- ggplot(prov_active, aes(ratio, province)) +
  geom_vline(xintercept = 0, colour = FAINT, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                colour = MUTED, linewidth = 0.45) +
  geom_point(aes(size = resolved), colour = PT_BLUE) +
  geom_text(aes(x = hi, label = sprintf("  %d/%d", confirmed, resolved)),
            hjust = 0, size = 2.3, colour = MUTED, family = base_family) +
  scale_size_area(max_size = 3.6, name = "Resolved alerts",
                  breaks = c(10, 100, 1000, 10000), labels = scales::comma) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.30))) +
  # Rows are ordered by the point estimate, which for a proportion estimated from four
  # alerts is mostly noise, so the denominator is printed beside every row and the count
  # of thinly-supported provinces is stated on the axis. A5C carries the same disclosure
  # on its tick labels.
  labs(x = sprintf("Confirmed share of laboratory-resolved alerts (labels: confirmed/resolved)\nordered by the point estimate; %d of the %d rest on fewer than %d resolved alerts",
                   sum(prov_active$resolved < 100), nrow(prov_active), 100L),
       y = NULL) +
  ann_in(sprintf("%d of %d provinces had no\nlaboratory-resolved alert",
                 N_PROV_SILENT, nrow(prov_pos)), x = -Inf, y = Inf, vjust = 1.1) +
  theme_pub() + theme(panel.grid.major.y = element_blank())

# Correlation matrix on ONE support: the exogenous proxies plus epicentre travel time,
# all defined for every zone that has a road route (Spearman, since several are skewed).
cm_vars <- c(PROXY_EXO$var, "epi_travel_h")
cm_lab  <- c(gsub("\n", " ", PROXY_EXO$label), "Epicentre travel time (h)")
cm_lab  <- stringr::str_wrap(gsub("\\s*\\(.*\\)$", "", cm_lab), 12)
cm_dat  <- zone_ascert %>% dplyr::select(dplyr::all_of(cm_vars)) %>% as.data.frame()
cm_ok   <- stats::complete.cases(cm_dat)
CM_N    <- sum(cm_ok)
cm <- stats::cor(cm_dat[cm_ok, , drop = FALSE], method = "spearman")
cm_df <- as.data.frame(as.table(cm)) %>%
  setNames(c("v1", "v2", "rho")) %>%
  dplyr::mutate(i = match(v1, cm_vars), j = match(v2, cm_vars)) %>%
  dplyr::filter(i >= j) %>%
  dplyr::mutate(v1 = factor(cm_lab[i], levels = cm_lab),
                v2 = factor(cm_lab[j], levels = rev(cm_lab)))

pA1d <- ggplot(cm_df, aes(v1, v2, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rho),
                colour = abs(rho) > 0.6), size = 2.35, family = base_family) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), breaks = c(-1, 0, 1),
                       name = sprintf("Spearman rho\n(n = %d zones)", CM_N)) +
  scale_x_discrete(position = "top") + labs(x = NULL, y = NULL) +
  # The -1 is definitional, not a finding: the surveillance gap IS 1 minus the percentile
  # rank of facility density. Both are nonetheless carried through the regressions, where
  # they are different linear predictors (a percentile and a log are monotone-related but
  # not linearly), so a result that held on only one of them would be a rank artefact.
  ann_in("rho = -1 by construction: the surveillance gap\nis 1 - the percentile rank of facility density",
         x = -Inf, y = -Inf, hjust = -0.03, vjust = -0.4, size = 2.1) +
  theme_pub() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 0, vjust = 0, size = 6.0),
        axis.text.y = element_text(size = 6.0), legend.position = "right",
        legend.key.height = unit(16, "pt"), legend.key.width = unit(7, "pt"))

FigA1 <- (pA1a | pA1b) / (pA1c | pA1d) +
  patchwork::plot_layout(heights = c(1.05, 0.95)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(pA1a, "FigureA1a_map_surveillance_gap", 4.2, 3.8)
save_dual(pA1b, "FigureA1b_map_alert_effort",     4.2, 3.8)
save_dual(pA1c, "FigureA1c_province_positivity",  4.2, 3.2)
save_dual(pA1d, "FigureA1d_proxy_correlations",   4.2, 3.2)
save_dual(FigA1, "FigureA1", 9.6, 8.0, dir = FIG_DIR)

# -----------------------------------------------------------------------------
# 6. FIGURE A2 — does thin surveillance delay the first confirmed case?
#    The support is the INVADED zones: every one of them has an arrival time, so the
#    outbreak-conditioned proxies are comparable across them here.
#
#    TWO PANELS. The added-variable plot and the added-variable forest were dropped
#    from the figure by editorial decision; the regressions behind them are still
#    fitted below and still published, as A2_arrival_avp_surveillance.csv and
#    A2_arrival_added_variable.csv and in the report, so no number is lost with the
#    panels. `avp_frame()` is retained for the same reason: it defines the residual
#    pair the published slope refers to.
# -----------------------------------------------------------------------------
message("\n[ascert] Figure A2 — surveillance and the timing of the first confirmed case")

inv <- zone_ascert %>% dplyr::filter(invaded, is.finite(arrival_days))
stopifnot(nrow(inv) == sum(zone_ascert$invaded))
inv_ne <- inv %>% dplyr::filter(!health_zone %in% EPICENTRE_ZONES)
message(sprintf("[ascert] arrival-time support: %d invaded zone(s); %d excluding the epicentre region",
                nrow(inv), nrow(inv_ne)))

# Added-variable plot: BOTH axes residualised on epicentre travel time, so the fitted
# slope IS the coefficient of the proxy in lm(arrival ~ travel + proxy) — a partial
# residual plot (raw x against residual y) would draw a line of a different slope from
# the statistic annotated beside it.
avp_frame <- function(d, prox) {
  ok <- is.finite(d$arrival_days) & is.finite(d$epi_travel_h) & is.finite(d[[prox]])
  dd <- d[ok, , drop = FALSE]
  tibble::tibble(zone = dd$health_zone,
                 ry = stats::resid(stats::lm(dd$arrival_days ~ dd$epi_travel_h)),
                 rx = stats::resid(stats::lm(dd[[prox]] ~ dd$epi_travel_h)))
}
avp_sg <- avp_frame(inv, "surveillance_gap")
fit_sg <- av_slope(inv$arrival_days, inv$epi_travel_h, inv$surveillance_gap)
rho_sg <- spearman_ci(avp_sg$rx, avp_sg$ry)

# Added-variable forest over every proxy, with and without the three epicentre zones.
av_forest <- purrr::pmap_dfr(
  list(PROXY_ALL$var, PROXY_ALL$label, PROXY_ALL$kind),
  function(v, lab, kd) dplyr::bind_rows(
    dplyr::mutate(av_slope(inv$arrival_days,    inv$epi_travel_h,    inv[[v]]),
                  set = "All invaded zones"),
    dplyr::mutate(av_slope(inv_ne$arrival_days, inv_ne$epi_travel_h, inv_ne[[v]]),
                  set = "Excluding the epicentre region")) %>%
    dplyr::mutate(var = v, label = gsub("\n", " ", lab), kind = kd))
av_forest <- av_forest %>%
  dplyr::mutate(set = factor(set, levels = c("All invaded zones", "Excluding the epicentre region")),
                estimable = is.finite(slope_sd))
.ord <- av_forest %>% dplyr::filter(set == "All invaded zones") %>%
  dplyr::arrange(dplyr::desc(dplyr::coalesce(slope_sd, -Inf))) %>% dplyr::pull(label)
av_forest <- av_forest %>%
  dplyr::mutate(label_f = factor(label, levels = rev(.ord)),
                lab_n = sprintf("n = %d", n))

# Arrival time against the zone's own onset -> confirmation delay.
del <- inv %>% dplyr::filter(is.finite(delay_median))
rho_del <- spearman_ci(del$delay_median, del$arrival_days)
fit_del <- av_slope(del$arrival_days, del$epi_travel_h, del$delay_median)
# The drawn line is an UNADJUSTED ordinary least-squares fit, which is neither of the two
# statistics annotated beside it, so its own slope test is computed and annotated as well.
# It is then drawn dashed and muted when that slope is not distinguishable from zero: a
# confident red trend through a null relationship is the one thing this panel must not
# suggest, since its point is that detection delay does NOT order arrival.
.ols_del <- stats::lm(arrival_days ~ delay_median, data = del)
.ols_p   <- summary(.ols_del)$coefficients["delay_median", "Pr(>|t|)"]
.ols_sig <- is.finite(.ols_p) && .ols_p < 0.05
.ols_col <- if (.ols_sig) FIT_RED else MUTED
.ols_lty <- if (.ols_sig) "solid" else "22"

pA2a <- ggplot(del, aes(delay_median, arrival_days)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = .ols_col, fill = .ols_col,
              alpha = 0.10, linewidth = 0.6, linetype = .ols_lty) +
  geom_point(aes(size = n_delay), colour = PT_BLUE, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = health_zone), size = 2.05, colour = MUTED,
                           family = base_family, max.overlaps = 8, seed = RANDOM_SEED,
                           min.segment.length = 0.3, segment.colour = FAINT,
                           segment.size = 0.25) +
  scale_size_area(max_size = 3.4, name = "Cases used") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = sprintf("Median onset to laboratory confirmation (days, capped at %d)", DELAY_CAP_DAYS),
       y = "Arrival time (days from the first national onset)") +
  ann_in(sprintf("Spearman rho %s   %s\nfitted line (unadjusted OLS): %s%s\nadjusted for travel time: %s days per SD, %s\nn = %d zones with >=5 native-onset cases",
                 fmt_ci(rho_del$rho, rho_del$lo, rho_del$hi), fmt_p(rho_del$p),
                 fmt_p(.ols_p), if (.ols_sig) "" else " (dashed: not distinguishable from zero)",
                 fmt_ci(fit_del$slope_sd, fit_del$lo, fit_del$hi, 1), fmt_p(fit_del$p),
                 nrow(del))) +
  theme_pub()

# Arrival-time ECDF by surveillance-gap tercile.
ter_lab <- c("Lower third\n(densest facilities)", "Middle third", "Upper third\n(sparsest facilities)")
inv_t <- inv %>% dplyr::mutate(ter = factor(ter_lab[dplyr::ntile(surveillance_gap, 3)],
                                            levels = ter_lab))
ter_med <- inv_t %>% dplyr::group_by(ter) %>%
  dplyr::summarise(n = dplyr::n(), med = stats::median(arrival_days), .groups = "drop")

pA2b <- ggplot(inv_t, aes(arrival_days, colour = ter)) +
  stat_ecdf(geom = "step", linewidth = 0.7, pad = FALSE) +
  geom_vline(data = ter_med, aes(xintercept = med, colour = ter),
             linetype = "22", linewidth = 0.4, show.legend = FALSE) +
  scale_colour_manual(values = unname(OKABE[c(1, 5, 2)]), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Arrival time (days from the first national onset)",
       y = "Share of invaded zones arrived") +
  ann_in(paste(sprintf("%s: median %.0f d (n = %d)",
                       gsub("\n", " ", ter_med$ter), ter_med$med, ter_med$n), collapse = "\n"),
         x = Inf, y = -Inf, hjust = 1.04, vjust = -0.3) +
  theme_pub() + theme(legend.position = "top", legend.direction = "horizontal")

FigA2 <- (pA2a | pA2b) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(pA2a, "FigureA2a_arrival_vs_delay",      4.6, 3.6)
save_dual(pA2b, "FigureA2b_arrival_ecdf_tercile",  4.6, 3.6)
save_dual(FigA2, "FigureA2", 9.4, 3.9, dir = FIG_DIR)

# -----------------------------------------------------------------------------
# 7. FIGURE A3 — are invaded and never-invaded zones comparable?
#    Only the EXOGENOUS proxies enter this figure. An invaded zone has a confirmed
#    case by definition, so every outbreak-conditioned quantity (confirmation
#    ratio, alert counts, confirmation delay) is mechanically tied to the outcome
#    and a contrast on it would measure the definition, not the surveillance system.
# -----------------------------------------------------------------------------
message("\n[ascert] Figure A3 — comparability of invaded and never-invaded zones")

za_t <- zone_ascert %>% dplyr::filter(is.finite(epi_travel_h))
message(sprintf("[ascert] comparability support: %d zones with a road route (%d invaded)",
                nrow(za_t), sum(za_t$invaded)))

cd_sg <- cliffs_delta(zone_ascert$surveillance_gap[zone_ascert$invaded],
                      zone_ascert$surveillance_gap[!zone_ascert$invaded])
mw_sg <- suppressWarnings(stats::wilcox.test(surveillance_gap ~ invaded, data = zone_ascert))

pA3a <- ggplot(zone_ascert, aes(status, surveillance_gap, fill = status, colour = status)) +
  geom_violin(alpha = 0.18, linewidth = 0.3, width = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  geom_jitter(width = 0.11, height = 0, size = 0.7, alpha = 0.38, stroke = 0) +
  scale_fill_manual(values = INV_COL, guide = "none") +
  scale_colour_manual(values = INV_COL, guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.03, 0.2))) +
  labs(x = NULL, y = "Surveillance gap (facility density percentile)") +
  ann_in(sprintf("Cliff's delta %s\nMann-Whitney %s\nn = %d invaded, %d never invaded",
                 fmt_ci(cd_sg$delta, cd_sg$lo, cd_sg$hi), fmt_p(mw_sg$p.value),
                 cd_sg$n_a, cd_sg$n_b)) +
  theme_pub()

cd_forest <- purrr::pmap_dfr(
  list(c(PROXY_EXO$var, "epi_travel_h"),
       c(PROXY_EXO$label, "Epicentre travel time\n(hours)")),
  function(v, lab) cliffs_delta(zone_ascert[[v]][zone_ascert$invaded],
                                zone_ascert[[v]][!zone_ascert$invaded]) %>%
    dplyr::mutate(var = v, label = gsub("\n", " ", lab)))
cd_forest <- cd_forest %>%
  dplyr::arrange(delta) %>%
  dplyr::mutate(label_f = factor(label, levels = label))

pA3b <- ggplot(cd_forest, aes(delta, label_f)) +
  geom_vline(xintercept = 0, colour = FAINT, linewidth = 0.35, linetype = "22") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                colour = MUTED, linewidth = 0.45) +
  geom_point(colour = FIT_RED, size = 1.9) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5),
                     expand = expansion(mult = c(0.03, 0.03))) +
  # Every row is estimated on the same two groups bar the travel-time row, where the two
  # zones with no road route drop out, so the support is stated once rather than per row.
  labs(x = sprintf("Cliff's delta, invaded vs never invaded (negative = lower in invaded zones)\n%d invaded vs %d never-invaded zones (%d vs %d for travel time)",
                   sum(zone_ascert$invaded), sum(!zone_ascert$invaded),
                   sum(za_t$invaded), sum(!za_t$invaded)),
       y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank())

# The same contrast INSIDE quartiles of epicentre travel time. Invaded zones are not
# spread across the country, so this is where the unadjusted contrast is shown to be
# identified only in the near field — which is itself the answer to "are they comparable".
TT_Q <- 4L
za_q <- za_t %>%
  dplyr::mutate(qi = dplyr::ntile(epi_travel_h, TT_Q)) %>%
  dplyr::group_by(qi) %>%
  dplyr::mutate(q_lab = sprintf("Q%d: %.0f-%.0f h", qi[1], min(epi_travel_h), max(epi_travel_h))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(q_lab = factor(q_lab, levels = unique(q_lab[order(qi)])))
q_stat <- za_q %>%
  dplyr::group_by(q_lab) %>%
  dplyr::group_modify(function(g, k) {
    a <- g$surveillance_gap[g$invaded]; b <- g$surveillance_gap[!g$invaded]
    cd <- cliffs_delta(a, b)
    dplyr::mutate(cd, note = if (length(a) < 5L)
      sprintf("%d invaded zone%s\n(delta not estimated)", length(a), if (length(a) == 1L) "" else "s")
      else sprintf("delta %s\nn = %d / %d", fmt_ci(cd$delta, cd$lo, cd$hi), cd$n_a, cd$n_b))
  }) %>% dplyr::ungroup()

pA3c <- ggplot(za_q, aes(status, surveillance_gap, colour = status, fill = status)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.16, linewidth = 0.4) +
  geom_jitter(width = 0.13, height = 0, size = 0.6, alpha = 0.4, stroke = 0) +
  geom_text(data = q_stat, aes(x = 1.5, y = 1.03, label = note), inherit.aes = FALSE,
            size = 2.05, colour = MUTED, family = base_family, vjust = 0, lineheight = 1.05) +
  facet_wrap(~ q_lab, nrow = 1) +
  scale_fill_manual(values = INV_COL, guide = "none") +
  scale_colour_manual(values = INV_COL, guide = "none") +
  scale_y_continuous(limits = c(0, 1.24), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0))) +
  scale_x_discrete(labels = c("Never\ninvaded", "Invaded")) +
  labs(x = "Quartile of road travel time from the epicentre region",
       y = "Surveillance gap") +
  theme_pub() + theme(axis.text.x = element_text(size = 6.6))

# BOTH models are fit on `za_t`, the zones with a road route. Fitting the unadjusted model
# on the full spine and the adjusted one on that subset would compare coefficients
# estimated on different samples AND standardised by SDs computed on different samples, so
# the pair would not be the adjustment contrast it is drawn as.
or_forest <- purrr::pmap_dfr(
  list(PROXY_EXO$var, PROXY_EXO$label),
  function(v, lab) dplyr::bind_rows(
    dplyr::mutate(logit_or(za_t$invaded, za_t[[v]]), model = "Unadjusted"),
    dplyr::mutate(logit_or(za_t$invaded, za_t[[v]], base = za_t$epi_travel_h),
                  model = "Adjusted for epicentre travel time")) %>%
    dplyr::mutate(var = v, label = gsub("\n", " ", lab)))
or_forest <- or_forest %>%
  dplyr::mutate(model = factor(model, levels = c("Unadjusted", "Adjusted for epicentre travel time")))
.ord_or <- or_forest %>% dplyr::filter(model == "Unadjusted") %>%
  dplyr::arrange(dplyr::coalesce(or_sd, Inf)) %>% dplyr::pull(label)
or_forest <- or_forest %>% dplyr::mutate(label_f = factor(label, levels = .ord_or))
OR_CI_TYPES <- sort(unique(stats::na.omit(or_forest$ci_type)))

pA3d <- ggplot(dplyr::filter(or_forest, is.finite(or_sd)),
               aes(or_sd, label_f, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = FAINT, linewidth = 0.35, linetype = "22") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                linewidth = 0.45, position = position_dodge(width = 0.55)) +
  geom_point(size = 1.7, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = unname(OKABE[c(1, 2)]), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_continuous(transform = "log10", breaks = c(0.25, 0.5, 1, 2, 4),
                     labels = c("0.25", "0.5", "1", "2", "4"),
                     expand = expansion(mult = c(0.06, 0.06))) +
  # The support note goes in the axis title, not an in-panel annotation: this axis is
  # log-transformed, and an annotation anchored at -Inf has no finite position on it.
  labs(x = sprintf("Odds ratio for being invaded, per SD of the proxy\n(%s intervals; %d invaded of %d zones with a road route)",
                   paste(OR_CI_TYPES, collapse = "/"), sum(za_t$invaded), nrow(za_t)),
       y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank(),
                      legend.position = "top", legend.direction = "horizontal")

FigA3 <- (pA3a | pA3b) / pA3c / pA3d +
  patchwork::plot_layout(heights = c(1, 0.85, 0.85)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(pA3a, "FigureA3a_surveillance_by_status",   4.0, 3.4)
save_dual(pA3b, "FigureA3b_cliffs_delta_forest",      4.6, 3.4)
save_dual(pA3c, "FigureA3c_contrast_within_travel_q", 8.2, 3.0)
save_dual(pA3d, "FigureA3d_logistic_or",              5.6, 3.0)
save_dual(FigA3, "FigureA3", 9.4, 9.6, dir = FIG_DIR)

# -----------------------------------------------------------------------------
# 8. FIGURE A4 — does the forecast rank by surveillance?
#    All ranks are taken within a (horizon, fold) watch-list, so they are exactly
#    invariant to the probability scale and to any monotone recalibration.
# -----------------------------------------------------------------------------
message("\n[ascert] Figure A4 — the forecast and the surveillance gradient")

H_MAIN <- 2L
stopifnot(H_MAIN %in% lfo_f$horizon)

# Zone-level correlation, NOT zone-fold level: the at-risk zone-folds are repeated
# measurements of the same zones, roughly one per fold per horizon, so a correlation over
# them would treat that replication as independent information and report an interval far
# narrower than the data support.
rank_zone <- lfo_f %>%
  dplyr::group_by(horizon, health_zone, surveillance_gap) %>%
  dplyr::summarise(mean_rank = mean(rank_avg), n_folds = dplyr::n(), .groups = "drop")
rho_rank <- rank_zone %>%
  dplyr::group_by(horizon) %>%
  dplyr::group_modify(function(g, k) spearman_ci(g$surveillance_gap, g$mean_rank)) %>%
  dplyr::ungroup()

d_main <- lfo_f %>% dplyr::filter(horizon == H_MAIN)
dec_med <- d_main %>%
  dplyr::mutate(dec = dplyr::ntile(surveillance_gap, 10)) %>%
  dplyr::group_by(dec) %>%
  dplyr::summarise(x = stats::median(surveillance_gap), y = stats::median(rank_avg),
                   .groups = "drop")
.rr <- rho_rank %>% dplyr::filter(horizon == H_MAIN)

pA4a <- ggplot(d_main, aes(surveillance_gap, rank_avg)) +
  geom_bin2d(bins = 34) +
  geom_line(data = dec_med, aes(x, y), colour = "white", linewidth = 1.1, inherit.aes = FALSE) +
  geom_line(data = dec_med, aes(x, y), colour = INK, linewidth = 0.55, inherit.aes = FALSE) +
  geom_point(data = dplyr::filter(d_main, is_new_invasion == 1),
             colour = FIT_RED, size = 1.5, shape = 21, fill = "white", stroke = 0.7) +
  scale_fill_viridis_c(option = "G", direction = -1, transform = "log10",
                       name = "At-risk\nzone-forecasts", labels = scales::comma) +
  scale_y_reverse(expand = expansion(mult = c(0.04, 0.03))) +
  # The panel is filled edge to edge, and an annotation anchored at +/-Inf has no finite
  # position on a REVERSED scale, so the legend text lives in the axis title instead.
  labs(x = sprintf("Surveillance gap (facility density percentile)\nblack line: median rank by decile;  red rings: %d realised invasion zone-forecasts (%d distinct zones)\nzone-level Spearman rho %s, %s (n = %d zones)",
                   sum(d_main$is_new_invasion),
                   dplyr::n_distinct(d_main$health_zone[d_main$is_new_invasion == 1]),
                   fmt_ci(.rr$rho, .rr$lo, .rr$hi), fmt_p(.rr$p), .rr$n),
       y = sprintf("Rank in the %d-week watch-list (1 = highest risk)", H_MAIN)) +
  theme_pub() + theme(legend.position = "right",
                      legend.key.height = unit(16, "pt"), legend.key.width = unit(7, "pt"))

# ONE POINT PER INVADED ZONE. At h = 2 the outcome windows of consecutive origins overlap,
# so a single invasion appears as a positive row in up to two folds — the gap between
# `n_invasions` and `n_invasion_events` in invasion_evaluation.csv. Correlating over the
# rows would treat that duplication as independent information and report an interval the
# data do not support, so each zone enters once, at its mean rank over the rounds in which
# it was the event. Both counts are printed on the panel.
ev_rows <- lfo_f %>% dplyr::filter(is_new_invasion == 1)
ev <- ev_rows %>%
  dplyr::group_by(horizon, health_zone, surveillance_gap) %>%
  dplyr::summarise(rank_avg = mean(rank_avg), n_rows = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(h_lab = factor(sprintf("%d-week horizon", horizon),
                               levels = sprintf("%d-week horizon", sort(unique(horizon)))))
.ev_rows_n <- ev_rows %>% dplyr::count(horizon, name = "n_rows_total")
# The drawn line is an OLS fit on the TRANSFORMED y axis (the panel is log10 in rank), so
# its slope test is computed the same way, and the line is dashed and muted wherever that
# slope is not distinguishable from zero. Drawing a confident trend through a null
# relationship is exactly the misreading this panel exists to prevent.
.ols_ev <- ev %>% dplyr::group_by(horizon) %>%
  dplyr::group_modify(function(g, k) {
    m <- stats::lm(log10(rank_avg) ~ surveillance_gap, data = g)
    tibble::tibble(ols_p = summary(m)$coefficients["surveillance_gap", "Pr(>|t|)"])
  }) %>% dplyr::ungroup() %>%
  dplyr::mutate(sig = is.finite(ols_p) & ols_p < 0.05)
ev <- ev %>% dplyr::left_join(.ols_ev, by = "horizon")
rho_ev <- ev %>% dplyr::group_by(h_lab, horizon) %>%
  dplyr::group_modify(function(g, k) spearman_ci(g$surveillance_gap, g$rank_avg)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(.ev_rows_n, by = "horizon") %>%
  dplyr::left_join(.ols_ev, by = "horizon") %>%
  dplyr::mutate(lab = sprintf("Spearman rho %s, %s\nfitted line (OLS on log rank): %s%s\n%d invaded zones (%d zone-forecasts)",
                              fmt_ci(rho, lo, hi), fmt_p(p), fmt_p(ols_p),
                              ifelse(sig, "", " (dashed)"), n, n_rows_total))

pA4b <- ggplot(ev, aes(surveillance_gap, rank_avg)) +
  geom_smooth(aes(colour = sig, fill = sig, linetype = sig),
              method = "lm", formula = y ~ x, alpha = 0.10, linewidth = 0.6) +
  scale_colour_manual(values = c(`TRUE` = FIT_RED, `FALSE` = MUTED), guide = "none") +
  scale_fill_manual(values = c(`TRUE` = FIT_RED, `FALSE` = MUTED), guide = "none") +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), guide = "none") +
  geom_point(colour = PT_BLUE, size = 1.5, alpha = 0.85) +
  geom_text(data = rho_ev, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.06, vjust = 1.2, size = 2.25, colour = MUTED,
            family = base_family, lineheight = 1.08) +
  facet_wrap(~ h_lab, nrow = 1) +
  scale_y_continuous(transform = "log10", breaks = c(1, 3, 10, 30, 100, 300),
                     expand = expansion(mult = c(0.06, 0.16))) +
  labs(x = "Surveillance gap of the invaded zone",
       y = "Mean rank of the invaded zone\nin the round(s) it was the event") +
  theme_pub()

# Within-stratum discrimination. Each stratum is scored against ITS OWN base rate, which
# is what AUC-PR SKILL is for, and the bootstrap resamples whole zones inside the stratum.
strat_tab <- lfo_f %>%
  dplyr::group_by(horizon, sg_stratum) %>%
  dplyr::group_modify(function(g, k) strat_discrimination(g)) %>%
  dplyr::ungroup()
message("[ascert] within-stratum discrimination:")
print(as.data.frame(strat_tab %>% dplyr::select(horizon, sg_stratum, n, n_pos, n_pos_zones,
                                                base_rate, auc_pr, ap_lo, ap_hi,
                                                auc_pr_skill, skill_lo, skill_hi, skill_ceiling,
                                                auc_roc, roc_lo, roc_hi)), digits = 3)

M_SKILL <- "AUC-PR skill (x base rate)"; M_AP <- "Average precision"; M_ROC <- "AUC-ROC"
strat_long <- dplyr::bind_rows(
  strat_tab %>% dplyr::transmute(horizon, sg_stratum, n_pos, n_pos_zones, base_rate,
                                 metric = M_SKILL, ref = 1, cap = skill_ceiling,
                                 est = auc_pr_skill, lo = skill_lo, hi = skill_hi),
  strat_tab %>% dplyr::transmute(horizon, sg_stratum, n_pos, n_pos_zones, base_rate,
                                 metric = M_AP, ref = base_rate, cap = NA_real_,
                                 est = auc_pr, lo = ap_lo, hi = ap_hi),
  strat_tab %>% dplyr::transmute(horizon, sg_stratum, n_pos, n_pos_zones, base_rate,
                                 metric = M_ROC, ref = 0.5, cap = NA_real_,
                                 est = auc_roc, lo = roc_lo, hi = roc_hi)) %>%
  dplyr::mutate(h_lab = factor(sprintf("%d wk", horizon),
                               levels = sprintf("%d wk", sort(unique(horizon)))),
                metric = factor(metric, levels = c(M_SKILL, M_AP, M_ROC)))
#' One metric per plot rather than one facet per metric: AUC-PR skill is a RATIO spanning
#' an order of magnitude at this event count and is only readable on a log axis, while the
#' other two are probabilities on a linear one. `scales = "free_y"` frees the limits but
#' not the transform, so facetting them together would force one or the other.
#'
#' The NO-SKILL reference is drawn per group, not as one line, because it differs between
#' strata for average precision (it is the stratum's base rate) even though it is common
#' for skill (1) and AUC-ROC (0.5). The skill panel additionally marks each stratum's
#' ATTAINABLE MAXIMUM, 1 / base rate: without it, the higher skill of the sparser stratum
#' reads as better discrimination when it is mostly a lower base rate.
.strat_panel <- function(metric_name, ylab, log_y, show_cap = FALSE, common_ref = TRUE,
                         xlab = "Forecast horizon") {
  d <- dplyr::filter(strat_long, metric == metric_name)
  dodge <- position_dodge(width = 0.5)
  p <- ggplot(d, aes(h_lab, est, colour = sg_stratum, shape = sg_stratum))
  p <- if (common_ref) {
    # A "common" reference must actually be common; two values here would draw two lines
    # and silently imply the null differs between strata when the caller said it does not.
    stopifnot(length(unique(d$ref)) == 1L)
    p + geom_hline(yintercept = d$ref[1], colour = FAINT, linewidth = 0.35, linetype = "22")
  }
  else
    p + geom_point(aes(y = ref), shape = 45, size = 5, alpha = 0.75,
                   position = dodge, show.legend = FALSE)
  if (show_cap)
    p <- p + geom_point(aes(y = cap), shape = 45, size = 5, alpha = 0.55,
                        position = dodge, show.legend = FALSE)
  p +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.45, position = dodge) +
    geom_point(size = 2, position = dodge) +
    geom_text(aes(y = hi, label = sprintf("%d", n_pos_zones)), vjust = -0.8, size = 2.1,
              position = dodge, show.legend = FALSE) +
    scale_colour_manual(values = unname(OKABE[c(1, 2)]), name = NULL) +
    scale_shape_manual(values = c(16, 17), name = NULL) +
    labs(x = xlab, y = ylab) + theme_pub() +
    theme(legend.position = "top", legend.direction = "horizontal") +
    (if (log_y) scale_y_continuous(transform = "log10", breaks = c(1, 3, 10, 30, 100, 300),
                                   expand = expansion(mult = c(0.08, 0.16)))
     else scale_y_continuous(expand = expansion(mult = c(0.12, 0.16))))
}
# The "numbers above" convention applies to all three panels, so it sits under the middle
# one, where it reads as a note on the row rather than on the panel it happens to touch.
pA4c <- (.strat_panel(M_SKILL, "AUC-PR skill (x base rate)\ndashes: no skill (1, lower) and the\nattainable maximum 1 / base rate (upper)",
                      TRUE, show_cap = TRUE) |
         .strat_panel(M_AP, "Average precision\ndashes: the stratum's own base rate,\nwhich is its no-skill value", FALSE,
                      common_ref = FALSE,
                      xlab = "Forecast horizon\nnumbers above each point: distinct invaded zones") |
         .strat_panel(M_ROC, "AUC-ROC\ndashes: 0.5 (no skill)", FALSE)) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.direction = "horizontal")

# Watch-list recall by stratum. The top-K list is NATIONAL (ranks come from the whole
# at-risk set of that round, tie-broken to the worst position exactly as
# compute_detection_curve does); the recall is then the share of THAT stratum's invasions
# the national list caught. Pooled over folds rather than averaged per fold, because most
# folds contribute no event at all to a stratum and the per-fold share is undefined there;
# the interval is the zone-cluster bootstrap in recall_curve_boot(), which is what the
# overlapping 2-week outcome windows require.
KS <- 1:25
rec_strat <- ev_rows %>%
  dplyr::group_by(horizon, sg_stratum) %>%
  dplyr::group_modify(function(g, key) recall_curve_boot(g, KS)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(h_lab = factor(sprintf("%d-week horizon", horizon),
                               levels = sprintf("%d-week horizon", sort(unique(horizon)))))
rand_ref <- lfo_f %>% dplyr::group_by(horizon, fold_id) %>%
  dplyr::summarise(n_atrisk = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(horizon) %>%
  dplyr::summarise(n_atrisk = mean(n_atrisk), .groups = "drop") %>%
  tidyr::expand_grid(k = KS) %>%
  dplyr::mutate(recall = pmin(k / n_atrisk, 1),
                h_lab = factor(sprintf("%d-week horizon", horizon),
                               levels = sprintf("%d-week horizon", sort(unique(horizon)))))

rec_band <- rec_strat %>%
  dplyr::group_by(horizon, sg_stratum, h_lab) %>%
  dplyr::group_modify(function(g, key) step_path(g)) %>%
  dplyr::ungroup()
pA4d <- ggplot(rec_strat, aes(k, recall, colour = sg_stratum, fill = sg_stratum)) +
  geom_ribbon(data = rec_band, aes(ymin = lo, ymax = hi), alpha = 0.14, colour = NA) +
  geom_step(linewidth = 0.7, direction = "hv") +
  geom_line(data = rand_ref, aes(k, recall), colour = MUTED, linetype = "22",
            linewidth = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ h_lab, nrow = 1) +
  scale_colour_manual(values = unname(OKABE[c(1, 2)]), name = NULL) +
  scale_fill_manual(values = unname(OKABE[c(1, 2)]), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Size K of the national weekly watch-list\ndashed: random targeting; bands are 95% zone-cluster bootstrap intervals",
       y = "Share of that stratum's\ninvasions caught") +
  theme_pub() + theme(legend.position = "top", legend.direction = "horizontal")

# pA4c is itself a patchwork; wrapped so the composite tags it as ONE panel (C) instead
# of handing separate tags to its three children. It gets a full row: three metric panels
# squeezed into half a row would make the interval whiskers unreadable.
FigA4 <- (pA4a | pA4b) / patchwork::wrap_elements(pA4c) / pA4d +
  patchwork::plot_layout(heights = c(1.15, 0.95, 0.85)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(pA4a, "FigureA4a_rank_vs_surveillance",   4.8, 3.6)
save_dual(pA4b, "FigureA4b_rank_of_truth",          4.6, 3.6)
save_dual(pA4c, "FigureA4c_stratified_skill",       8.6, 3.6)
save_dual(pA4d, "FigureA4d_recall_by_stratum",      8.6, 3.2)
save_dual(FigA4, "FigureA4", 10.0, 11.4, dir = FIG_DIR)

# -----------------------------------------------------------------------------
# 9. FIGURE A5 — province test positivity, detection opportunity, and the ranking
# -----------------------------------------------------------------------------
message("\n[ascert] Figure A5 — province positivity and the ranking")

# A. Was the system looking where the model pointed? For every at-risk zone-forecast,
#    whether that zone had ALREADY had a discarded (not-a-case) alert investigated by the
#    forecast origin — a surveillance-effort measure that cannot be created by the zone's
#    own later invasion, unlike "has a resolved alert" measured at the end of the outbreak.
BANDS <- c(0, 5, 10, 15, 25, 50, 100, 250, Inf)
BAND_LAB <- c("1-5", "6-10", "11-15", "16-25", "26-50", "51-100", "101-250", "251+")
watch_band <- lfo_f %>%
  dplyr::filter(horizon == H_MAIN) %>%
  dplyr::mutate(band = cut(rank_avg, breaks = BANDS, labels = BAND_LAB, right = TRUE)) %>%
  dplyr::group_by(band) %>%
  dplyr::summarise(n = dplyr::n(), watched = sum(watched_prior),
                   events = sum(is_new_invasion), .groups = "drop")
.bw <- wilson_ci(watch_band$watched, watch_band$n)
watch_band <- watch_band %>% dplyr::mutate(share = .bw$est, lo = .bw$lo, hi = .bw$hi)

pA5a <- ggplot(watch_band, aes(band, share)) +
  geom_col(fill = PT_BLUE, alpha = 0.78, width = 0.68) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = INK, linewidth = 0.4) +
  geom_text(aes(y = hi, label = scales::comma(n)), vjust = -0.7, size = 2.1,
            colour = MUTED, family = base_family) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.16))) +
  labs(x = sprintf("Rank band in the %d-week watch-list\nnumbers above: at-risk zone-forecasts in the band", H_MAIN),
       y = "Share with a discarded alert already\nnotified at the forecast origin") +
  theme_pub()

# B. Province confirmation ratio against how the model ranks that province's at-risk zones.
prov_rank <- lfo_f %>%
  dplyr::filter(horizon == H_MAIN) %>%
  dplyr::group_by(province) %>%
  dplyr::summarise(mean_rank = mean(rank_avg), n_zonefolds = dplyr::n(),
                   n_events = sum(is_new_invasion), .groups = "drop")
prov_join <- prov_pos %>%
  dplyr::left_join(prov_rank, by = "province") %>%
  dplyr::filter(!is.na(mean_rank))
prov_act <- prov_join %>% dplyr::filter(resolved > 0)
rho_prov <- spearman_ci(prov_act$ratio, prov_act$mean_rank)

pA5b <- ggplot(prov_act, aes(ratio, mean_rank)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                colour = FAINT, linewidth = 0.4) +
  geom_point(aes(size = resolved), colour = PT_BLUE, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = province), size = 2.2, colour = MUTED,
                           family = base_family, seed = RANDOM_SEED, max.overlaps = 20,
                           min.segment.length = 0.2, segment.colour = FAINT, segment.size = 0.25) +
  scale_size_area(max_size = 4, name = "Resolved alerts", labels = scales::comma,
                  breaks = c(100, 1000, 10000)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_reverse() +
  # Three short lines rather than two long ones: at this panel width the single-line
  # version was clipped at the right edge, losing the end of the statistic.
  labs(x = sprintf("Province confirmed share of laboratory-resolved alerts\n95%% Wilson bars\nSpearman rho %s, %s (%d provinces)",
                   fmt_ci(rho_prov$rho, rho_prov$lo, rho_prov$hi), fmt_p(rho_prov$p), rho_prov$n),
       y = sprintf("Mean rank of the province's\nat-risk zones (%d-week watch-list)", H_MAIN)) +
  theme_pub()

# C. The ranking distribution province by province, ordered by confirmation ratio, with the
#    provinces that never resolved an alert pooled into one group — there is no positivity
#    to order them by, and showing 20 near-identical boxes would imply a precision the
#    data do not carry.
SILENT_LAB <- "No laboratory-\nresolved alert"
prov_levels <- prov_act %>% dplyr::arrange(dplyr::desc(ratio)) %>% dplyr::pull(province)
# The tick label carries the ratio AND its denominator. Ordering seven provinces by a
# proportion estimated from a handful of resolved alerts would otherwise read as a
# ranking, which the printed denominator makes impossible to miss.
.prov_tick <- prov_act %>%
  dplyr::transmute(province,
                   tick = sprintf("%s\n%.0f%% (%d/%d)", province, 100 * ratio, confirmed, resolved))
prov_tick_lab <- stats::setNames(c(.prov_tick$tick[match(prov_levels, .prov_tick$province)],
                                   SILENT_LAB), c(prov_levels, SILENT_LAB))
rank_prov <- lfo_f %>%
  dplyr::filter(horizon == H_MAIN) %>%
  dplyr::mutate(prov_grp = factor(dplyr::if_else(province %in% prov_levels, province, SILENT_LAB),
                                  levels = c(prov_levels, SILENT_LAB)),
                outcome = factor(dplyr::if_else(is_new_invasion == 1,
                                                "Invaded within the window", "Not invaded"),
                                 levels = c("Not invaded", "Invaded within the window")))
prov_n <- rank_prov %>% dplyr::group_by(prov_grp) %>%
  dplyr::summarise(n = dplyr::n(), ev = sum(is_new_invasion),
                   ev_zones = dplyr::n_distinct(health_zone[is_new_invasion == 1]),
                   .groups = "drop") %>%
  dplyr::mutate(lab = sprintf("%s zone-forecasts\n%d invaded zone%s", scales::comma(n),
                              ev_zones, ifelse(ev_zones == 1L, "", "s")))

pA5c <- ggplot(rank_prov, aes(prov_grp, rank_avg)) +
  geom_boxplot(data = dplyr::filter(rank_prov, outcome == "Not invaded"),
               width = 0.55, outlier.shape = NA, fill = "grey88", colour = MUTED,
               linewidth = 0.35) +
  geom_jitter(data = dplyr::filter(rank_prov, outcome == "Invaded within the window"),
              width = 0.13, height = 0, colour = FIT_RED, size = 1.5, shape = 21,
              fill = "white", stroke = 0.7) +
  geom_text(data = prov_n, aes(x = prov_grp, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.1, size = 1.95, colour = MUTED, family = base_family, lineheight = 1.05) +
  scale_y_continuous(transform = "log10", breaks = c(1, 3, 10, 30, 100, 300),
                     expand = expansion(mult = c(0.05, 0.22))) +
  scale_x_discrete(labels = prov_tick_lab) +
  labs(x = paste0("Province, ordered by confirmed share of laboratory-resolved alerts (highest first)\n",
                  "grey boxes: at-risk zone-forecasts that were not invaded;  red rings: the realised invasions\n",
                  "the last group can contain NO invasion by construction: a confirmed case is itself a laboratory-resolved alert"),
       y = sprintf("Rank in the %d-week\nwatch-list", H_MAIN)) +
  theme_pub() + theme(axis.text.x = element_text(size = 5.9, lineheight = 1.0))

# D. The live watch-list against laboratory activity: where an operational list points.
live_h <- live %>% dplyr::filter(horizon == H_MAIN) %>%
  dplyr::mutate(rank_live = rank(-p_case_invasion, ties.method = "min"))
TOP_N <- 25L
# slice_min(with_ties = FALSE) would break a tie at the K boundary by row order, which is
# zone name after the joins. Assert there is no tie to break rather than hide one.
.p_sorted <- sort(live_h$p_case_invasion, decreasing = TRUE)
stopifnot(length(.p_sorted) > TOP_N, .p_sorted[TOP_N] > .p_sorted[TOP_N + 1L])
# The bars are coloured by whether the ZONE ITSELF has ever had an alert investigated to a
# laboratory outcome, which is the granularity an operational list is acted on; the province
# is in the tick label. Colouring by the PROVINCE's laboratory activity was uninformative
# here, because every zone in the live top-25 sits in a province that has some.
live_top <- live_h %>% dplyr::slice_min(rank_live, n = TOP_N, with_ties = FALSE) %>%
  dplyr::left_join(dplyr::select(prov_pos, province, prov_ratio = ratio,
                                 prov_resolved = resolved), by = "province") %>%
  dplyr::mutate(lab_status = factor(dplyr::if_else(ever_investigated,
                                                   "Zone has had an alert laboratory-resolved",
                                                   "Zone has had none"),
                                    levels = c("Zone has had an alert laboratory-resolved",
                                               "Zone has had none")),
                zone_f = forcats::fct_reorder(sprintf("%s (%s)", health_zone, province),
                                              p_case_invasion))
message(sprintf("[ascert] live top-%d at h=%d: %d zone(s) have had an alert laboratory-resolved, %d have not; %d of %d sit in a province with any laboratory activity",
                TOP_N, H_MAIN, sum(live_top$ever_investigated), sum(!live_top$ever_investigated),
                sum(live_top$prov_resolved > 0), nrow(live_top)))

pA5d <- ggplot(live_top, aes(p_case_invasion, zone_f, fill = lab_status)) +
  geom_col(width = 0.7, alpha = 0.88) +
  scale_fill_manual(values = unname(c(OKABE[1], OKABE[2])), name = NULL, drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = sprintf("Live %d-week invasion probability (%s scale)", H_MAIN, PROB_SCALE),
       y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank(),
                      axis.text.y = element_text(size = 6.2),
                      legend.position = "top", legend.direction = "horizontal")

FigA5 <- (pA5a | pA5b) / pA5c / pA5d +
  patchwork::plot_layout(heights = c(1, 0.9, 1.15)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = INK))
save_dual(pA5a, "FigureA5a_investigated_by_rank_band", 4.4, 3.4)
save_dual(pA5b, "FigureA5b_province_ratio_vs_rank",    4.6, 3.6)
save_dual(pA5c, "FigureA5c_rank_by_province",          8.4, 3.2)
save_dual(pA5d, "FigureA5d_live_watchlist_lab_status", 5.4, 4.0)
save_dual(FigA5, "FigureA5", 9.6, 11.0, dir = FIG_DIR)

# -----------------------------------------------------------------------------
# 10. WRITE EVERY NUMBER THAT APPEARS IN A PANEL
#     A figure is not a result until the values behind it can be quoted, so each
#     panel's underlying table is written out beside the figures.
# -----------------------------------------------------------------------------
message("\n[ascert] writing panel data")
.write_panel <- function(x, name) {
  readr::write_csv(x, file.path(FIG_DIR, name), na = "")
  message(sprintf("  wrote %-42s %d row(s)", name, nrow(x)))
}

.write_panel(zone_ascert %>%
     dplyr::select(health_zone, province, invaded, epi_travel_h, pop_count,
                   healthsite_count, healthsite_density, log_hs_density, log_pop,
                   surveillance_gap, healthcare_gap, social_vulnerability, ccvi,
                   alerts_total, n_confirmed, n_not_a_case, alerts_resolved,
                   confirm_ratio, confirm_lo, confirm_hi, resolution_frac,
                   neg_per_100k, log_neg_100k, resolved_per_100k, ever_investigated,
                   n_delay, delay_median, delay_p75, delay_over_cap, delay_censored,
                   coverage_ratio, coverage_usable,
                   first_onset, arrival_days, travel_time_h, cases_total),
   "zone_ascertainment.csv")
.write_panel(av_forest %>% dplyr::select(var, label, kind, set, slope_sd, lo, hi, p, dR2, n),
   "A2_arrival_added_variable.csv")
.write_panel(tibble::tibble(panel = "A2a", stat = c("slope_per_sd", "ci_lo", "ci_hi", "p", "delta_R2",
                                          "spearman_rho", "rho_lo", "rho_hi", "rho_p", "n"),
                  value = c(fit_sg$slope_sd, fit_sg$lo, fit_sg$hi, fit_sg$p, fit_sg$dR2,
                            rho_sg$rho, rho_sg$lo, rho_sg$hi, rho_sg$p, fit_sg$n)),
   "A2_arrival_avp_surveillance.csv")
.write_panel(inv_t %>% dplyr::select(health_zone, province, arrival_days, surveillance_gap,
                           delay_median, n_delay, tercile = ter),
   "A2_arrival_by_zone.csv")
.write_panel(cd_forest %>% dplyr::arrange(var) %>%
               dplyr::select(var, label, delta, lo, hi, n_a, n_b),
   "A3_cliffs_delta.csv")
.write_panel(q_stat %>% dplyr::select(q_lab, delta, lo, hi, n_a, n_b),
   "A3_cliffs_delta_within_travel_quartile.csv")
.write_panel(or_forest %>% dplyr::select(var, label, model, or_sd, lo, hi, p, n, n_pos, ci_type),
   "A3_logistic_or.csv")
.write_panel(strat_tab, "A4_stratified_discrimination.csv")
.write_panel(rec_strat %>% dplyr::select(horizon, sg_stratum, k, caught, events, n_zones,
                                        recall, lo, hi),
   "A4_watchlist_recall_by_stratum.csv")
.write_panel(rho_rank, "A4_rank_vs_surveillance_zone_level.csv")
.write_panel(rho_ev %>% dplyr::select(h_lab, rho, lo, hi, p, n_zones = n, n_rows_total),
   "A4_rank_of_truth_vs_surveillance.csv")
.write_panel(dec_med %>% dplyr::rename(surveillance_gap_median = x, rank_median = y),
   "A4_median_rank_by_surveillance_decile.csv")
.write_panel(watch_band, "A5_investigation_by_rank_band.csv")
.write_panel(prov_join %>% dplyr::select(province, n_zones, n_invaded, confirmed, resolved,
                               ratio, lo, hi, mean_rank, n_zonefolds, n_events),
   "A5_province_positivity_and_rank.csv")
.write_panel(prov_n %>% dplyr::rename(province_group = prov_grp), "A5_rank_by_province_counts.csv")
.write_panel(ev %>% dplyr::select(horizon, health_zone, surveillance_gap,
                                  mean_rank_of_truth = rank_avg, n_rows),
   "A4_rank_of_truth_by_zone.csv")
.write_panel(live_top %>% dplyr::select(health_zone, province, rank_live, p_case_invasion,
                              surveillance_gap, ever_investigated, neg_alerts,
                              prov_ratio, prov_resolved),
   "A5_live_top_watchlist.csv")

.input_stamp <- function(paths) {
  paths <- paths[file.exists(paths)]
  stats::setNames(lapply(paths, function(f) list(
    mtime = format(file.info(f)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    bytes = as.numeric(file.info(f)$size),
    md5   = unname(tools::md5sum(f)))), basename(paths))
}
INPUT_FILES <- c(fs_risk_csv(OUT),
                 file.path(OUT, "key_outputs", "arrival_predictors.csv"),
                 file.path(OUT, "key_outputs", "coverage_ratio_zone.csv"),
                 file.path(OUT, "key_outputs", "model_selection.json"),
                 file.path(OUT, "forecasts", "lfo_results.rds"))
.lfo_stamp <- attr(lfo_all, "lfo_stamp")

provenance <- list(
  module = "45_ascertainment_figures.R",
  inputs = .input_stamp(INPUT_FILES),
  lfo_stamp = if (is.null(.lfo_stamp)) NULL else
    .lfo_stamp[vapply(.lfo_stamp, function(v) length(v) == 1L && !is.list(v), logical(1))],
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  analysis_date = format(ANALYSIS_DATE), outbreak_start = format(OUTBREAK_START),
  linelist_folder = jsonlite::fromJSON(LINELIST_JSON)$folder,
  featured_model = FEATURED, risk_table_method = RISK_METHOD,
  forecast_scale = FORECAST_SCALE, lfo_probability_column = PCOL,
  risk_prob_scale = PROB_SCALE,
  conf_level = CONF_LEVEL, n_boot = N_BOOT, random_seed = RANDOM_SEED,
  delay_cap_days = DELAY_CAP_DAYS,
  dropped_degenerate_pillars = DROPPED_PILLARS,
  n_zones = nrow(zone_ascert), n_invaded = sum(zone_ascert$invaded),
  n_zones_no_road_route = as.integer(.n_na_travel),
  n_synthetic_sitrep_rows_excluded = as.integer(n_synth),
  n_alerts_in_window = nrow(ll_win),
  n_resolved_alerts = sum(zone_ascert$alerts_resolved),
  n_zones_ever_investigated = sum(zone_ascert$ever_investigated),
  n_never_invaded_never_investigated =
    sum(!zone_ascert$ever_investigated & !zone_ascert$invaded),
  n_provinces_with_lab_activity = as.integer(nrow(prov_pos) - N_PROV_SILENT),
  n_provinces_silent = as.integer(N_PROV_SILENT),
  surveillance_stratum_split = SG_SPLIT,
  lfo_folds = as.integer(dplyr::n_distinct(lfo_f$fold_id)),
  lfo_zone_folds = nrow(lfo_f),
  lfo_zone_folds_watched_prior = sum(lfo_f$watched_prior),
  live_top_n = as.integer(TOP_N),
  live_top_never_investigated = sum(!live_top$ever_investigated))
writeLines(jsonlite::toJSON(provenance, auto_unbox = TRUE, pretty = TRUE, null = "null"),
           file.path(FIG_DIR, "ascertainment_run_info.json"))
message("  wrote ascertainment_run_info.json")

message("\n[ascert] done — figures and panel data in ", FIG_DIR)

# -----------------------------------------------------------------------------
# 11. REPORT — the headline numbers, composed from the objects above.
#     Every value is interpolated from a computed object; nothing here is typed
#     by hand, so the report cannot drift from the figures.
# -----------------------------------------------------------------------------
# EVERY DIRECTIONAL AND QUALITATIVE WORD BELOW IS DERIVED, not typed. An earlier version
# interpolated the numbers but hard-coded the adjectives ("LOWER", "HIGHER", "falls
# steeply"), and when the pipeline was re-run one of them became the opposite of what the
# data said while the sentence still read as generated output. A report that can state the
# wrong direction with the right number is worse than no report.
.dir <- function(a, b, hi = "higher", lo = "lower", eq = "the same") {
  if (!is.finite(a) || !is.finite(b)) return("not comparable")
  if (isTRUE(all.equal(a, b))) eq else if (a > b) hi else lo
}
.ci_overlap <- function(lo1, hi1, lo2, hi2)
  all(is.finite(c(lo1, hi1, lo2, hi2))) && lo1 <= hi2 && lo2 <= hi1
.rho_word <- function(rho, p, alpha = 0.05) {
  if (!is.finite(rho)) return("not estimated")
  if (!is.finite(p) || p >= alpha) return("not distinguishable from zero")
  if (rho > 0) "positive" else "negative"
}
.excl <- function(lo, hi, null_val) if (!all(is.finite(c(lo, hi)))) "not estimated" else
  if (lo > null_val || hi < null_val) "excludes" else "includes"

.rep_sg  <- rho_rank %>% dplyr::filter(horizon == H_MAIN)
.cd_q1   <- q_stat %>% dplyr::slice(1)
.or_u    <- or_forest %>% dplyr::filter(var == "surveillance_gap", model == "Unadjusted")
.or_a    <- or_forest %>% dplyr::filter(var == "surveillance_gap", model != "Unadjusted")
.st      <- strat_tab
.ter     <- inv_t %>% dplyr::group_by(ter) %>%
  dplyr::summarise(n = dplyr::n(), med = stats::median(arrival_days), .groups = "drop")
.band1   <- watch_band %>% dplyr::slice(1)
.bandlast<- watch_band %>% dplyr::slice(dplyr::n())
.q_last_inv <- max(za_q$qi[za_q$invaded])
.cmp <- .st %>%
  dplyr::mutate(k = dplyr::if_else(grepl("^Sparser", as.character(sg_stratum)), "sp", "dn")) %>%
  dplyr::select(horizon, k, auc_pr, ap_lo, ap_hi, auc_roc, roc_lo, roc_hi) %>%
  tidyr::pivot_wider(names_from = k,
                     values_from = c(auc_pr, ap_lo, ap_hi, auc_roc, roc_lo, roc_hi)) %>%
  dplyr::arrange(horizon)
.cmp_sent <- function(est_sp, est_dn, lo_sp, hi_sp, lo_dn, hi_dn, d = 3) {
  paste(sprintf(paste0("h=%d %s in the sparser stratum (%.", d, "f vs %.", d, "f; the two intervals %s)"),
                .cmp$horizon,
                mapply(.dir, est_sp, est_dn),
                est_sp, est_dn,
                ifelse(mapply(.ci_overlap, lo_sp, hi_sp, lo_dn, hi_dn), "overlap", "do not overlap")),
        collapse = "; ")
}
.any_separation <- any(!mapply(.ci_overlap, .cmp$ap_lo_sp, .cmp$ap_hi_sp, .cmp$ap_lo_dn, .cmp$ap_hi_dn)) ||
                   any(!mapply(.ci_overlap, .cmp$roc_lo_sp, .cmp$roc_hi_sp, .cmp$roc_lo_dn, .cmp$roc_hi_dn))

report <- c(
"# Detection versus onset: spatial ascertainment and the invasion forecast",
"",
sprintf("Generated %s from the %s line-list snapshot, analysis date %s. Featured model `%s`; LFO probabilities on the `%s` column; all intervals %.0f%%. The exact input files, with modification times and MD5 sums, are recorded in `ascertainment_run_info.json` — the pipeline artifacts this report reads are rewritten by `run_all.R` and `43_spread_kinematics.R`, so a number quoted from here must be reconciled against that stamp.",
        format(Sys.Date()), jsonlite::fromJSON(LINELIST_JSON)$folder, format(ANALYSIS_DATE),
        FEATURED, PCOL, 100 * CONF_LEVEL),
"",
"## 1. The estimand",
"",
sprintf("The modelled outcome is the onset week of a zone's first **confirmed** case. Of %s line-list records in the outbreak window, %s reached a laboratory outcome (%s confirmed, %s discarded); %s did not reach one (still suspected or probable, or never classified). Confirmation ratios below are computed on the resolved alerts only, excluding the %d synthetic sitrep top-up rows, which are a count reconciliation against the official cumulative and carry no test.",
        scales::comma(nrow(ll_win)),
        scales::comma(sum(zone_ascert$alerts_resolved)),
        scales::comma(sum(zone_ascert$n_confirmed)),
        scales::comma(sum(zone_ascert$n_not_a_case)),
        scales::comma(nrow(ll_win) - sum(zone_ascert$alerts_resolved)), n_synth),
"",
"## 2. Ascertainment is spatially structured (Figure A1)",
"",
sprintf("- %d of %d health zones never had a single laboratory-resolved alert during the outbreak, including %d of the %d never-invaded zones. In those zones an invasion could not have been detected.",
        sum(!zone_ascert$ever_investigated), nrow(zone_ascert),
        sum(!zone_ascert$ever_investigated & !zone_ascert$invaded), sum(!zone_ascert$invaded)),
sprintf("- %d of %d provinces resolved no alert at all. Among the %d that did, the confirmed share of resolved alerts runs from %.0f%% (%s, %d/%d) to %.0f%% (%s, %d/%d); the denominators matter more than the ordering.",
        N_PROV_SILENT, nrow(prov_pos), nrow(prov_active),
        100 * min(prov_active$ratio), prov_active$province[which.min(prov_active$ratio)],
        prov_active$confirmed[which.min(prov_active$ratio)], prov_active$resolved[which.min(prov_active$ratio)],
        100 * max(prov_active$ratio), prov_active$province[which.max(prov_active$ratio)],
        prov_active$confirmed[which.max(prov_active$ratio)], prov_active$resolved[which.max(prov_active$ratio)]),
"",
"## 3. Surveillance density and the timing of the first confirmed case (Figure A2)",
"",
sprintf("- Added to a baseline of road travel time from the epicentre region, the surveillance gap changes arrival time by %s days per SD (%s), changing R-squared by %+.3f over %d invaded zones. Not plotted: the added-variable panels were dropped from Figure A2; the fits are in `A2_arrival_avp_surveillance.csv` and `A2_arrival_added_variable.csv`.",
        fmt_ci(fit_sg$slope_sd, fit_sg$lo, fit_sg$hi, 1), fmt_p(fit_sg$p), fit_sg$dR2, fit_sg$n),
sprintf("- Median arrival time by surveillance-gap tercile: %s.",
        paste(sprintf("%s = %.0f d (n = %d)", gsub("\n", " ", .ter$ter), .ter$med, .ter$n), collapse = "; ")),
sprintf("- The zone's own median onset-to-confirmation delay, as an ordering of arrival time, is %s (Spearman rho %s, %s, n = %d zones)%s",
        .rho_word(rho_del$rho, rho_del$p),
        fmt_ci(rho_del$rho, rho_del$lo, rho_del$hi), fmt_p(rho_del$p), nrow(del),
        if (identical(.rho_word(rho_del$rho, rho_del$p), "not distinguishable from zero"))
          " — on these data the arrival ordering above is therefore not explained by a zone-level reporting delay."
        else " — this is NOT null, so the arrival ordering above must not be presented as free of reporting delay."),
"- CAVEAT, stated as a limit of the design: facility density also proxies remoteness and urbanicity, so this is consistent with ascertainment delay but is not identified against genuinely later epidemiological arrival.",
"",
"## 4. Invaded and never-invaded zones are not comparable (Figure A3)",
"",
sprintf("- Unadjusted, the surveillance gap of invaded zones is %s than that of never-invaded zones (Cliff's delta %s, which %s zero; Mann-Whitney %s).",
        .dir(cd_sg$delta, 0, "higher", "lower", "no different"),
        fmt_ci(cd_sg$delta, cd_sg$lo, cd_sg$hi), .excl(cd_sg$lo, cd_sg$hi, 0), fmt_p(mw_sg$p.value)),
sprintf("- No invaded zone lies beyond quartile %d of the %d epicentre-travel-time quartiles, and %d of %d are in Q1, so the contrast is identified only in the near field. Inside Q1 the delta is %s, which %s zero (n = %d invaded vs %d never invaded).",
        .q_last_inv, TT_Q, sum(za_q$invaded & za_q$qi == 1L), sum(za_q$invaded),
        fmt_ci(.cd_q1$delta, .cd_q1$lo, .cd_q1$hi), .excl(.cd_q1$lo, .cd_q1$hi, 0),
        .cd_q1$n_a, .cd_q1$n_b),
sprintf("- Odds ratio for being invaded, per SD of surveillance gap: %s unadjusted (interval %s 1), %s adjusted for epicentre travel time (interval %s 1).",
        fmt_ci(.or_u$or_sd, .or_u$lo, .or_u$hi), .excl(.or_u$lo, .or_u$hi, 1),
        fmt_ci(.or_a$or_sd, .or_a$lo, .or_a$hi), .excl(.or_a$lo, .or_a$hi, 1)),
"",
"## 5. The forecast and the surveillance gradient (Figure A4)",
"",
sprintf("- The association between a zone's surveillance gap and its mean rank in the %d-week watch-list is %s (zone-level Spearman rho %s, %s, n = %d zones).%s",
        H_MAIN, .rho_word(.rep_sg$rho, .rep_sg$p),
        fmt_ci(.rep_sg$rho, .rep_sg$lo, .rep_sg$hi), fmt_p(.rep_sg$p), .rep_sg$n,
        if (identical(.rho_word(.rep_sg$rho, .rep_sg$p), "positive"))
          " A larger rank number is a lower place on the list, so zones with sparser facilities are ranked lower."
        else if (identical(.rho_word(.rep_sg$rho, .rep_sg$p), "negative"))
          " A larger rank number is a lower place on the list, so zones with sparser facilities are ranked higher."
        else ""),
sprintf("- Among the zones that were invaded, the association between the surveillance gap and the rank of the truth is: %s.",
        paste(sprintf("%s %s (rho %s, %s, %d zones)", rho_ev$h_lab,
                      mapply(.rho_word, rho_ev$rho, rho_ev$p),
                      fmt_ci(rho_ev$rho, rho_ev$lo, rho_ev$hi), fmt_p(rho_ev$p), rho_ev$n),
              collapse = "; ")),
"- Within-stratum discrimination must be read on a scale whose ceiling does not move with the base rate. AUC-PR skill is average precision divided by the stratum's own base rate, so 1 is the no-skill value in both strata, but the ATTAINABLE MAXIMUM is 1 / base rate and is far higher in the sparser one. The raw skill numbers must not be read as a ratio of discrimination.",
sprintf("  - AUC-PR skill (attainable maximum in brackets): %s.",
        paste(sprintf("h=%d %s %s (max %.0f)", .st$horizon, .st$sg_stratum,
                      fmt_ci(.st$auc_pr_skill, .st$skill_lo, .st$skill_hi, 1),
                      .st$skill_ceiling), collapse = "; ")),
sprintf("  - Average precision, the same quantity with that ceiling divided out: %s.",
        .cmp_sent(.cmp$auc_pr_sp, .cmp$auc_pr_dn, .cmp$ap_lo_sp, .cmp$ap_hi_sp,
                  .cmp$ap_lo_dn, .cmp$ap_hi_dn, 3)),
sprintf("  - AUC-ROC, which shares a [0, 1] scale and a 0.5 null across strata: %s.",
        .cmp_sent(.cmp$auc_roc_sp, .cmp$auc_roc_dn, .cmp$roc_lo_sp, .cmp$roc_hi_sp,
                  .cmp$roc_lo_dn, .cmp$roc_hi_dn, 3)),
sprintf("  - Supported by %s distinct invaded zones. %s",
        paste(unique(sprintf("%s %d", .st$sg_stratum, .st$n_pos_zones)), collapse = " and "),
        if (.any_separation)
          "At least one scale-free comparison separates the strata; read the per-horizon intervals above before drawing a conclusion."
        else
          "No scale-free comparison separates the strata at this event count: every interval overlaps, so the defensible statement is that discrimination is not detectably degraded where surveillance is thin — not that it differs in either direction."),
"",
"## 6. Where the surveillance system was looking (Figure A5)",
"",
sprintf("- Of %s at-risk zone-forecasts, %s (%.1f%%) were for a zone that already had a discarded alert notified at the forecast origin.",
        scales::comma(nrow(lfo_f)), scales::comma(sum(lfo_f$watched_prior)),
        100 * mean(lfo_f$watched_prior)),
sprintf("- Moving down the watch-list, that share is %s at the bottom than at the top: %.1f%% [%.1f, %.1f] in rank band %s against %.1f%% [%.1f, %.1f] in band %s.",
        .dir(.bandlast$share, .band1$share, "higher", "lower", "unchanged"),
        100 * .band1$share, 100 * .band1$lo, 100 * .band1$hi, .band1$band,
        100 * .bandlast$share, 100 * .bandlast$lo, 100 * .bandlast$hi, .bandlast$band),
sprintf("- In the live %d-week watch-list, %d of the top %d zones have never had a single alert laboratory-resolved.",
        H_MAIN, sum(!live_top$ever_investigated), TOP_N),
"",
"## 7. What these figures do and do not establish",
"",
"- They establish that ascertainment varies strongly in space, and that invaded and never-invaded zones differ on capacity in a way that is largely explained by proximity to the epicentre.",
"- They do NOT identify ascertainment delay against genuinely later epidemiological arrival: no measurement of unobserved transmission exists in these data, and every proxy used here is also a proxy for remoteness.",
"- Quantities generated BY the outbreak (confirmation ratio, alert counts, confirmation delay) are never used in an invaded-versus-not comparison, where they would be mechanically tied to the outcome; they appear only among already-invaded zones or re-measured at each forecast origin.",
"- In Figure A5C the final province group can contain no invasion by construction: a confirmed case is itself a laboratory-resolved alert, so a province with an invasion necessarily has laboratory activity.")

writeLines(report, file.path(FIG_DIR, "ASCERTAINMENT_REPORT.md"))
message("  wrote ASCERTAINMENT_REPORT.md")
