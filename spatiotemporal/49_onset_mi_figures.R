# =============================================================================
# 49_onset_mi_figures.R — the onset multiple-imputation DIAGNOSTIC figure suite
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# These are NOT the pooled Figures 2/3/4 — those are drawn by make_manuscript_figures.R and
# 40_cascade_next_dominoes.R with MANU_IN_ROOT pointed at the pooled mirror. This suite draws
# the five figures that answer the question the manuscript has to answer directly: HOW MUCH
# of what we report is onset-imputation uncertainty?
#
#   FigureMI1  skill across completions: MI estimate, total and between-completion spread
#   FigureMI2  per-zone invasion-WEEK movement — the outcome itself moving
#   FigureMI3  live-forecast rank stability across completions
#   FigureMI4  Figure 2A's recall curve: within-arm band vs MI total band
#   FigureMI5  cascade expected-added: within-arm interval vs MI total interval
#
# READ ONLY. Every number comes from the tables 48_onset_mi_pool.R publishes; nothing is
# fitted, pooled or re-derived here. A figure must not fit a model.
#
# CAPTIONS ARE WRITTEN HERE TOO, from the same tables and the pooled run_info.json: every
# number, date and seed in FigureMI<n>_caption.md is read, none is typed. The hand-written
# captions they replace described an imputation mechanism and a fold count that no longer
# existed, and nothing regenerated them. A builder that fails or skips removes its caption,
# so a caption never sits beside a figure it does not describe.
#
# Run:  Rscript spatiotemporal/49_onset_mi_figures.R
# Out:  outputs/key_outputs/manuscript_figures/mi/FigureMI<n>_<stem>.{pdf,png},
#       FigureMI<n>_caption.md
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(scales); library(here)
})

HERE   <- file.path(here::here(), "spatiotemporal")
OUT    <- file.path(HERE, "outputs")
MI_OUT <- Sys.getenv("ONSET_MI_OUT", file.path(OUT, "onset_mi"))
TBL    <- file.path(MI_OUT, "tables")
POOL   <- file.path(MI_OUT, "pooled")

source(file.path(HERE, "00_config.R"))
stop_if_data_definition_sensitivity("49_onset_mi_figures.R")   # hard-coded published paths
source(file.path(HERE, "forecast_scale.R"))
FIG_DIR   <- fs_out_dir(file.path(OUT, "key_outputs", "manuscript_figures", "mi"))
PANEL_DIR <- fs_out_dir(file.path(OUT, "key_outputs", "manuscript_figures", "mi", "panels"))
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
source(file.path(HERE, "figure_style.R"))
# mi_tcrit(): the SAME critical value 48 used for the published bands, so a within-component
# band drawn here cannot come from a different reference distribution than the total.
source(file.path(HERE, "onset_mi_pooling.R"))


# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("49_onset_mi_figures.R")) {
  message("[mi-fig] sourced, not run: nothing is written. Use `Rscript 49_onset_mi_figures.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

.req <- function(p, what) {
  if (!file.exists(p))
    stop("[mi-fig] ", what, " not found (", p, "). Run 48_onset_mi_pool.R first.", call. = FALSE)
  suppressWarnings(readr::read_csv(p, show_col_types = FALSE))
}
if (file.exists(file.path(POOL, ".pool_incomplete")))
  stop("[mi-fig] the last pooling (48_onset_mi_pool.R) did not complete, so ", POOL, " mixes ",
       "this sweep's products with an earlier one's; re-run 48.", call. = FALSE)
RI <- jsonlite::fromJSON(file.path(POOL, "key_outputs", "run_info.json"))
MI <- RI$onset_mi
if (is.null(MI)) stop("[mi-fig] the pooled run_info.json carries no `onset_mi` block.", call. = FALSE)
# Every field the stamp and the captions read. A pooled run_info written before 48 recorded
# them cannot describe the current imputation, so it is a refusal, not a caption about the old.
.mi_need <- c("n_imputations", "seed_min", "seed_max", "random_seed", "stratum", "method",
              "analysis_date", "linelist_snapshot", "onset_impute", "draws_per_arm",
              "prob_scale", "eval_ci", "generated_at")
.mi_miss <- .mi_need[!vapply(.mi_need, function(f)
  length(MI[[f]]) > 0L && !all(is.na(unlist(MI[[f]]))), logical(1))]
if (length(.mi_miss))
  stop("[mi-fig] the pooled run_info.json onset_mi block lacks ", paste(.mi_miss, collapse = ", "),
       "; it predates the current 48_onset_mi_pool.R. Re-run 48 first.", call. = FALSE)
M <- as.integer(MI$n_imputations)
SUB <- sprintf(paste0("M = %d onset completions (onset seeds %d-%d varied, model seed %d held ",
                      "fixed; stratum '%s'); model %s"),
               M, as.integer(MI$seed_min), as.integer(MI$seed_max), as.integer(MI$random_seed),
               MI$stratum, MI$method)
message("[mi-fig] ", SUB)
# THE MI STAMP. save_dual() stamps a figure (figure_style.R::.fig_mi_stamp()) only when
# MI_STAMP_TXT is set, and nothing here set it, so these figures left their folder with no
# statement that they describe an imputation sweep. Set before any builder runs.
MI_STAMP_TXT <- paste0("Onset multiple imputation: ", SUB, ".")

# =============================================================================
# CAPTION HELPERS — formatting only; every value passed in is read from a table
# =============================================================================
.or <- function(a, b) if (is.null(a) || !length(a) || all(is.na(a))) b else a
.cap_path <- function(n) file.path(FIG_DIR, sprintf("FigureMI%d_caption.md", n))
.fdate <- function(d) {
  d <- as.Date(d)
  sprintf("%d %s %s", as.integer(format(d, "%d")), format(d, "%B"), format(d, "%Y"))
}
.sig <- function(x, d = 4) vapply(as.numeric(x), function(v)
  if (is.finite(v)) format(signif(v, d), trim = TRUE, scientific = FALSE) else "\u2014",
  character(1))
.pct <- function(x, d = 1) vapply(as.numeric(x), function(v)
  if (is.finite(v)) sprintf(paste0("%.", d, "f%%"), 100 * v) else "\u2014", character(1))
.int <- function(x) format(as.integer(round(as.numeric(x))), big.mark = ",", trim = TRUE)
.df_txt <- function(x) vapply(as.numeric(x), function(v)
  if (is.infinite(v)) "\u221e" else if (is.finite(v)) format(round(v)) else "\u2014", character(1))
.md_table <- function(df) {
  m <- as.matrix(as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE,
                               check.names = FALSE))
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(m, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
# One sentence-list: "a", "a and b", "a, b and c".
.and <- function(x) if (length(x) <= 1L) paste(x) else
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])

# THE IMPUTATION, as the arms ran it: the mechanism string load_linelist() published
# (01_data_prep.R .ONSET_IMPUTE_DESC, carried by every arm and verified identical by 48), the
# donor floor, the sample-date source, the counts from the data-layer table, and both seeds.
.imputation_text <- function() {
  dl <- .req(file.path(TBL, "onset_mi_datalayer.csv"), "the data-layer table")
  nc <- unique(dl$n_confirmed); ni <- unique(dl$n_onset_imputed)
  if (length(nc) != 1L || length(ni) != 1L)
    stop("[mi-fig] the confirmed or imputed-onset count differs between completions; the ",
         "imputation moves onset DATES, never which records are imputed.", call. = FALSE)
  oi <- MI$onset_impute
  if (!identical(as.character(oi$mode), "impute"))
    stop("[mi-fig] the arms ran onset mode '", oi$mode, "', not an imputation.", call. = FALSE)
  sds <- as.character(oi$sample_date_source)
  sd_txt <- switch(sds,
    earliest = paste("the earliest diagnostic sample date recorded for the case across the",
                     "line-list snapshots up to the analysis date"),
    extract  = "the diagnostic sample date in the line-list extract",
    sprintf("the diagnostic sample date (source '%s')", sds))
  # n_onset_imputed counts confirmed records with a sample date but no usable onset
  # (01_data_prep.R: onset_imputed = !onset_usable & !is.na(date_of_sample_collection)). It
  # includes the sitrep top-up rows, which are a count reconciliation dated by the situation
  # report, not case records with a diagnostic sample, so the two are described apart.
  ns <- unique(dl$n_sitrep_confirmed); nl <- unique(dl$n_onset_imputed_linelist)
  if (length(ns) != 1L || length(nl) != 1L || !is.finite(ns) || !is.finite(nl))
    stop("[mi-fig] the data-layer table does not separate the sitrep top-up rows from the line-list ",
         "records (arms written before those counts were recorded); re-run the sweep (47, then 48).",
         call. = FALSE)
  if (nl + ns != ni)
    stop(sprintf("[mi-fig] line-list (%d) + sitrep (%d) imputed records do not add to the %d imputed.",
                 as.integer(nl), as.integer(ns), as.integer(ni)), call. = FALSE)
  sprintf(paste0("%s of the %s confirmed line-list cases (%.0f%%) carry a diagnostic sample date but ",
                 "no usable symptom-onset date, and in each completion every one of them is imputed ",
                 "as %s; a sample week with fewer than ",
                 "%d such donor records borrows the adjacent weeks. The %s sitrep top-up rows (a ",
                 "reconciliation of the line list to the situation reports' confirmed counts, dated ",
                 "by report date) carry no onset and are imputed the same way. The sample date is %s. Each ",
                 "record's draw is a deterministic function of its identifier and the onset ",
                 "seed (ONSET_IMPUTE_SEED), so the M = %d completions vary that seed alone ",
                 "(%d\u2013%d) and hold the model seed (RANDOM_SEED = %d), which seeds the model ",
                 "fitting, bootstrap and recalibration draws, fixed."),
          .int(nl), .int(nc - ns), 100 * nl / (nc - ns), as.character(oi$desc),
          as.integer(oi$donor_min), .int(ns),
          sd_txt, M, as.integer(MI$seed_min), as.integer(MI$seed_max),
          as.integer(MI$random_seed))
}
.rubin_text <- paste0("Estimates are combined by Rubin's rules: Q\u0304 is the mean of the M ",
                      "per-completion estimates, B the between-completion variance, \u016a the ",
                      "mean within-completion variance and T = \u016a + (1 + 1/M)B the total.")
.lambda_text <- paste0("\u03bb = (1 + 1/M)B/T is the between-completion share of the total ",
                       "variance. It slightly overstates the share attributable to the onset ",
                       "dates: the completions hold the model seed fixed, but each completion's ",
                       "model fits, bootstraps and recalibration draws run on a different ",
                       "completed data set, so their finite-draw Monte Carlo error enters B too.")
.provenance <- function() sprintf(paste0(
  "Provenance: pooled run_info.json written %s; analysis date %s; line-list snapshot `%s`; ",
  "onset seeds %d\u2013%d; model seed %d; stratum `%s`; model `%s`; %s posterior draws retained ",
  "per completion; current-forecast probability scale `%s`."),
  MI$generated_at, .fdate(MI$analysis_date), MI$linelist_snapshot, as.integer(MI$seed_min),
  as.integer(MI$seed_max), as.integer(MI$random_seed), MI$stratum, MI$method,
  .int(MI$draws_per_arm), MI$prob_scale)
.write_caption <- function(n, title, paras, values, sources) {
  sources <- unique(c(sources, "onset_mi_datalayer.csv"))   # .imputation_text() reads it
  txt <- c(sprintf("# Figure MI%d \u2014 caption", n), "",
           sprintf(paste0("Written by 49_onset_mi_figures.R from the pooled run_info.json and ",
                          "%s (M = %d). Every number below is read from them; none is typed ",
                          "by hand."), .and(sprintf("`%s`", sources)), M),
           "", "---", "", "## Caption (manuscript / SI)", "",
           paste(c(sprintf("**%s**", title), paras), collapse = " "),
           "", "---", "", "## Values shown (for checking the typeset figure)", "",
           values, "", .provenance())
  con <- file(.cap_path(n), open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(txt), con)
  message("  saved ", basename(.cap_path(n)))
  invisible(.cap_path(n))
}

ACCENT <- "#0072B2"; WARM <- "#D55E00"; GREEN <- "#009E73"

# =============================================================================
# MI1 — skill across completions
# =============================================================================
# Draws the VARIANCE DECOMPOSITION, not a histogram of arms: the MI estimate with the Rubin
# total SD sqrt(T) as the thin bar and the BETWEEN-COMPLETION COMPONENT OF T,
# sqrt(T - U) = sqrt((1 + 1/M)B), as the thick one, so the thick bar is exactly the part of the
# thin one that is between completions. (It used to be sqrt(B), which is not a component of T:
# T carries (1 + 1/M)B.) Neither
# lambda nor FMI is drawn or printed on the figure (see NO IN-FIGURE ANNOTATION below); both
# are in onset_mi_evaluation_variability.csv. The between-completion spread slightly OVERSTATES
# the onset-imputation share: the completions share one model seed, but each one's fits run on
# different data, so their finite-draw Monte Carlo error enters it too (47's header).
build_mi1 <- function() {
  unlink(.cap_path(1L))
  ev <- .req(file.path(TBL, "onset_mi_evaluation_variability.csv"), "the evaluation variability table")
  # FILTER TO THE FEATURED MODEL. The evaluation table also carries the three structural
  # baselines (they are scored, not merely drawn), and faceting over all four would put a
  # rank-only baseline's AUC-PR beside the model's under one axis label.
  ev <- ev %>% filter(method == MI$method)
  if (!nrow(ev)) stop("[mi-fig] the variability table carries no rows for ", MI$method,
                      call. = FALSE)
  keep <- c("auc_pr_skill", "auc_pr", "auc_roc", "mean_rank_of_truth", "base_rate")
  LBL <- c(auc_pr_skill = "AUC-PR skill (x base rate)", auc_pr = "AUC-PR",
           auc_roc = "AUC-ROC", mean_rank_of_truth = "Mean rank of truth",
           base_rate = "Base rate (invaded share)")
  # NO IN-FIGURE ANNOTATION. House convention (see 40_cascade_next_dominoes.R): the figure
  # carries its identity in the strips and axis titles; everything else is the caption's job.
  # The caption for this figure is FigureMI1_caption.md, written below beside the PDF.
  #
  # WHICH METRICS CAN BE DECOMPOSED IS STILL VISIBLE, without a marker.
  # invasion_evaluation.csv publishes a zone-cluster bootstrap interval (evaluate_invasion()
  # resamples HEALTH ZONES) only for the metrics 48 lists in onset_mi$eval_ci$metrics (the stems
  # that carry a bootstrap interval; read from there, not restated); every other metric arrives with no
  # within-completion variance: U = 0 by construction, T = (1 + 1/M)B, and lambda is EXACTLY 1
  # whatever the data say.
  # Those panels therefore get ONE bar (the between-completion spread, which is also the
  # total); only a decomposable panel gets the thick inner bar. The distinction reads off the
  # geometry and the caption names it — which is better than printing "lambda 1.00" beside a
  # measured "lambda 0.21", where the two look like comparable measurements and are not.
  d0 <- ev %>% filter(metric %in% keep)
  .dec <- d0 %>% group_by(metric) %>% summarise(decomposable = all(u_known), .groups = "drop")
  d <- d0 %>% left_join(.dec, by = "metric") %>%
    mutate(metric = factor(LBL[metric], levels = unname(LBL[keep])),
           horizon = factor(sprintf("%d week ahead", horizon)))
  p <- ggplot(d, aes(x = horizon, y = q, colour = horizon)) +
    geom_linerange(aes(ymin = q - sqrt(pmax(T, 0)), ymax = q + sqrt(pmax(T, 0))),
                   linewidth = 0.8, alpha = 0.45) +
    # The between-completion component OF T: T - U = (1 + 1/M)B.
    geom_linerange(data = ~ dplyr::filter(.x, decomposable),
                   aes(ymin = q - sqrt(pmax(T - U, 0)), ymax = q + sqrt(pmax(T - U, 0))),
                   linewidth = 2.2, alpha = 0.85) +
    geom_point(size = 2.6) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = c(ACCENT, WARM), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    labs(x = NULL, y = "MI estimate") +
    theme_pub(10) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_dual(p, "FigureMI1_skill_across_imputations", 12.4, 4.2, dir = FIG_DIR)

  # ---- caption ----
  shown <- unname(LBL[intersect(keep, unique(d0$metric))])
  dec_m <- intersect(keep, .dec$metric[.dec$decomposable %in% TRUE])
  hz <- sort(unique(d0$horizon))
  nf <- ev %>% filter(metric == "n_folds") %>% pull(q) %>% unique()
  rounds <- if (length(nf) == 1L && is.finite(nf))
    sprintf("%s leave-future-out cross-validation rounds", .int(nf)) else
    "the leave-future-out cross-validation rounds"
  lam_txt <- vapply(dec_m, function(m) {
    r <- d0 %>% filter(metric == m) %>% arrange(horizon)
    sprintf("%s: \u03bb = %s", LBL[[m]],
            .and(sprintf("%.2f at %d week%s", r$lambda, r$horizon, ifelse(r$horizon == 1L, "", "s"))))
  }, character(1))
  ci_lvl <- as.numeric(MI$eval_ci$level)
  dec_txt <- if (length(dec_m))
    sprintf(paste0("Of the metrics shown, only %s carr%s a within-completion interval (each ",
                   "completion's %.0f%% zone-cluster bootstrap interval, resampling health zones), ",
                   "so only %s can be decomposed (%s). For the other metrics \u016a is not ",
                   "available, the two components cannot be separated, and the single bar is the ",
                   "between-completion spread alone (\u03bb = 1 by construction, not by ",
                   "measurement)."),
            .and(unname(LBL[dec_m])), if (length(dec_m) == 1L) "ies" else "y", 100 * ci_lvl,
            if (length(dec_m) == 1L) "that panel" else "those panels", paste(lam_txt, collapse = "; "))
  else paste0("None of the metrics shown carries a within-completion interval, so each single ",
              "bar is the between-completion spread alone (\u03bb = 1 by construction).")
  br_txt <- if (any(c("base_rate", "auc_pr_skill") %in% d0$metric))
    paste0("AUC-PR skill and the base rate are not fixed targets: the base rate is computed from ",
           "the imputed labels, so their spread across completions is uncertainty about what ",
           "happened, not sampling error about a known truth.") else NULL
  vt <- d0 %>% left_join(.dec, by = "metric") %>%
    mutate(o = match(metric, keep)) %>% arrange(o, horizon) %>%
    transmute(metric = unname(LBL[metric]), h = horizon, estimate = .sig(q, 4),
              sd_total = .sig(sqrt(pmax(T, 0)), 3),
              sd_between = .sig(sqrt(pmax(T - U, 0)), 3),
              lam = ifelse(decomposable %in% TRUE, sprintf("%.3f", lambda), "\u2014"),
              df = .df_txt(df))
  names(vt) <- c("metric", "h", "estimate", "total SD (\u221aT)",
                 "between SD (\u221a(T \u2212 \u016a))", "\u03bb", "df")
  .write_caption(1L,
    title = sprintf(paste0("Figure MI1. Cross-validated performance of the featured invasion ",
                           "model (%s) under multiple imputation of missing symptom-onset ",
                           "dates."), MI$method),
    paras = c(.imputation_text(),
              sprintf(paste0("Because a health zone's invasion date is the first confirmed ",
                             "onset among its cases, the imputation moves both the model's ",
                             "inputs and the outcome it is scored against, so the evaluation ",
                             "over %s (analysis date %s) was repeated on every completion."),
                      rounds, .fdate(MI$analysis_date)),
              .rubin_text,
              sprintf(paste0("Panels: %s. Points are Q\u0304 at the %s forecast horizon%s. ",
                             "Thin bars are \u00b1\u221aT, one total standard deviation. The ",
                             "thick inner bar, drawn only where a within-completion interval ",
                             "exists, is \u00b1\u221a(T \u2212 \u016a) = \u00b1\u221a((1 + 1/M)B), ",
                             "the between-completion component of T."),
                      .and(shown), .and(sprintf("%d-week", hz)), if (length(hz) > 1L) "s" else ""),
              dec_txt, .lambda_text, br_txt),
    values = .md_table(vt),
    sources = "onset_mi_evaluation_variability.csv")
  invisible(p)
}

# =============================================================================
# MI2 — the OUTCOME moving: per-zone invasion-week range across completions
# =============================================================================
build_mi2 <- function() {
  unlink(.cap_path(2L))
  w <- .req(file.path(TBL, "onset_mi_invasion_week_shift.csv"), "the invasion-week table")
  moved <- sum(w$week_range > 0); n <- nrow(w)
  pa <- ggplot(w, aes(week_range)) +
    geom_histogram(binwidth = 1, fill = ACCENT, colour = "white", linewidth = 0.3) +
    scale_x_continuous(breaks = scales::pretty_breaks(6)) +
    labs(x = "Invasion-week range across completions (weeks)", y = "Health zones",
         subtitle = sprintf("%d of %d zones (%.0f%%) are dated to more than one week",
                            moved, n, 100 * moved / n)) +
    theme_pub(11) + theme(plot.subtitle = element_text(size = 8, colour = MUTED))
  # The 25 most volatile zones, as an interval per zone. These are the labels the model is
  # scored against, so a zone with a four-week range is a zone whose fold membership — and
  # therefore whether a forecast counted as a hit — changes between completions.
  top <- w %>% arrange(desc(week_range), week_median) %>% head(25) %>%
    mutate(zone = factor(health_zone, levels = rev(health_zone)))
  pb <- ggplot(top, aes(y = zone)) +
    geom_linerange(aes(xmin = week_min, xmax = week_max), linewidth = 1.6,
                   colour = WARM, alpha = 0.55) +
    geom_point(aes(x = week_median), size = 2, colour = INK) +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    labs(x = "Epidemic week of first confirmed onset", y = NULL,
         subtitle = "25 most imputation-sensitive zones (point = median completion)") +
    theme_pub(11) + theme(panel.grid.major.y = element_blank(),
                          axis.text.y = element_text(size = 7.6),
                          plot.subtitle = element_text(size = 8, colour = MUTED))
  fig <- (pa | pb) + plot_layout(widths = c(1, 1.05)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureMI2_invasion_week_shift", 12.0, 6.4, dir = FIG_DIR)

  # ---- caption ----
  rng_max <- max(w$week_range)
  dist <- table(factor(w$week_range, levels = 0:rng_max))
  n3 <- sum(w$week_range >= 3)
  every <- all(w$n_arms == M)
  ob_start <- .or(RI$outbreak_start, format(OUTBREAK_START))
  wide <- w %>% filter(week_range == rng_max) %>% arrange(health_zone)
  every_txt <- if (every)
    sprintf(paste0("Every one of the %d zones is invaded in all %d completions, so the panels ",
                   "show uncertainty in when each zone was invaded, not in whether it was."), n, M)
  else sprintf(paste0("%d of the %d zones are invaded in only some completions; their ranges are ",
                      "over the completions in which they are invaded."), sum(w$n_arms < M), n)
  vt1 <- data.frame(a = "health zones", t(as.integer(dist)), check.names = FALSE)
  names(vt1) <- c("invasion-week range (weeks)", names(dist))
  vt2 <- top %>% transmute(zone = health_zone,
                           weeks = sprintf("%d\u2013%d", as.integer(week_min), as.integer(week_max)),
                           median = .sig(week_median, 3), range = as.integer(week_range),
                           invaded = sprintf("%d of %d", as.integer(n_arms), M))
  names(vt2) <- c("zone", "epidemic weeks (min\u2013max)", "median week", "range (weeks)",
                  "completions in which invaded")
  .write_caption(2L,
    title = sprintf(paste0("Figure MI2. The outcome itself moves: health-zone invasion weeks ",
                           "across M = %d onset completions."), M),
    paras = c(paste0("A health zone's invasion date is the date of the first confirmed symptom ",
                     "onset among its cases, so imputing missing onsets moves not only the ",
                     "model's inputs but the labels it is scored against."),
              .imputation_text(),
              sprintf(paste0("**(A)** Distribution, across the %d invaded health zones, of the ",
                             "range of epidemic weeks of first confirmed onset spanned by the %d ",
                             "completions: %d zones are dated to the same week in every ",
                             "completion and %d (%.0f%%) are not; %d span three weeks or more, ",
                             "and the widest range is %d week%s (%s)."),
                      n, M, as.integer(dist[["0"]]), moved, 100 * moved / n, n3, rng_max,
                      if (rng_max == 1L) "" else "s", .and(wide$health_zone)),
              sprintf(paste0("**(B)** The %d most imputation-sensitive zones: the bar is the ",
                             "range of epidemic weeks of first confirmed onset across ",
                             "completions and the point is the median completion."), nrow(top)),
              every_txt,
              paste0("The invasion week determines which cross-validation round a zone becomes ",
                     "an event in, and therefore whether a forecast counts as a hit, so this is ",
                     "the mechanism behind the between-completion variance in Figures MI1 and MI4."),
              sprintf(paste0("Epidemic weeks are counted on the analysis week grid, week 1 being ",
                             "the week that contains the outbreak start (%s)."), .fdate(ob_start))),
    values = c(.md_table(vt1), "", .md_table(vt2)),
    sources = "onset_mi_invasion_week_shift.csv")
  invisible(fig)
}

# =============================================================================
# MI3 — live-forecast rank stability
# =============================================================================
build_mi3 <- function() {
  unlink(.cap_path(3L))
  r <- .req(file.path(TBL, "onset_mi_rank_stability.csv"), "the rank-stability table")
  rho <- .req(file.path(TBL, "onset_mi_rank_agreement.csv"), "the rank-agreement table")
  # share_top20 is over COMPLETIONS (48: a completion in which the zone is already invaded
  # counts as not top-20, the D5 rule), so the axis below says what it is. A table written
  # before that fix divided by the at-risk completions only and has no n_completions column.
  if (!all(c("n_completions", "share_top20") %in% names(r)) || !"ref_arm" %in% names(rho))
    stop("[mi-fig] the rank tables predate the completion-denominator fix; re-run ",
         "48_onset_mi_pool.R.", call. = FALSE)
  top <- r %>% arrange(rank_median) %>% head(25) %>%
    mutate(zone = factor(health_zone, levels = rev(health_zone)))
  pa <- ggplot(top, aes(y = zone)) +
    geom_linerange(aes(xmin = rank_min, xmax = rank_max), linewidth = 1.6,
                   colour = ACCENT, alpha = 0.5) +
    geom_point(aes(x = rank_median), size = 2, colour = INK) +
    geom_vline(xintercept = 20.5, linetype = "22", colour = MUTED, linewidth = 0.4) +
    scale_x_continuous(breaks = scales::pretty_breaks(6)) +
    labs(x = "2-week invasion-risk rank (1 = highest)", y = NULL,
         subtitle = "top 25 by median rank; dashed = the top-20 cut Figure 3B draws") +
    theme_pub(11) + theme(panel.grid.major.y = element_blank(),
                          axis.text.y = element_text(size = 7.6),
                          plot.subtitle = element_text(size = 8, colour = MUTED))
  # How reliably does a zone make the watch-list at all? A zone in the top 20 in 11 of 20
  # completions is a different operational recommendation from one in all 20, and the
  # single-draw figure cannot distinguish them.
  # The reference completion is the PUBLISHED draw only when its onset seed is the published
  # run's; otherwise the subtitle names the arm (as the caption below does).
  .ref_lab <- if (nrow(rho)) local({
    ref <- unique(rho$ref_arm)
    dl <- .req(file.path(TBL, "onset_mi_datalayer.csv"), "the data-layer table")
    ref_seed <- dl$seed[match(ref, dl$arm)]
    pub_seed <- suppressWarnings(as.integer(RI$onset_imputation$onset_impute_seed))
    if (length(pub_seed) == 1L && !is.na(pub_seed) && identical(as.integer(ref_seed), pub_seed))
      "the published draw" else sprintf("arm %d", as.integer(ref))
  }) else NA_character_
  pb <- ggplot(top, aes(share_top20, zone)) +
    geom_col(fill = GREEN, width = 0.7, alpha = 0.85) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.04))) +
    labs(x = "Share of completions in which the zone is top-20", y = NULL,
         subtitle = if (nrow(rho))
             sprintf("ranking of the other %d completion(s) vs %s: Spearman %s",
                     nrow(rho), .ref_lab,
                     if (nrow(rho) == 1L) sprintf("%.3f", rho$spearman_vs_arm0)
                     else sprintf("%.3f-%.3f", min(rho$spearman_vs_arm0, na.rm = TRUE),
                                  max(rho$spearman_vs_arm0, na.rm = TRUE)))
           else "only one completion carried risk scores") +
    theme_pub(11) + theme(panel.grid.major.y = element_blank(),
                          axis.text.y = element_blank(),
                          plot.subtitle = element_text(size = 8, colour = MUTED))
  fig <- (pa | pb) + plot_layout(widths = c(1.35, 1)) + plot_annotation(tag_levels = "A")
  save_dual(fig, "FigureMI3_rank_stability", 11.0, 6.4, dir = FIG_DIR)

  # ---- caption ----
  nC <- unique(r$n_completions)
  if (length(nC) != 1L) stop("[mi-fig] n_completions is not constant in the rank table.", call. = FALSE)
  n_all <- sum(top$share_top20 >= 1 - 1e-9)
  rho_txt <- if (nrow(rho)) {
    ref <- unique(rho$ref_arm)
    dl <- .req(file.path(TBL, "onset_mi_datalayer.csv"), "the data-layer table")
    ref_seed <- dl$seed[match(ref, dl$arm)]
    pub_seed <- suppressWarnings(as.integer(RI$onset_imputation$onset_impute_seed))
    is_pub <- length(pub_seed) == 1L && !is.na(pub_seed) && identical(as.integer(ref_seed), pub_seed)
    sprintf(paste0("Across the whole at-risk set, Spearman's \u03c1 between the ranking of each ",
                   "of the other %d completion(s) and that of arm %d (onset seed %d%s) ranges ",
                   "from %.4f to %.4f (median %.4f)."),
            nrow(rho), as.integer(ref), as.integer(ref_seed),
            if (is_pub) ", the published run's onset seed" else "",
            min(rho$spearman_vs_arm0, na.rm = TRUE), max(rho$spearman_vs_arm0, na.rm = TRUE),
            stats::median(rho$spearman_vs_arm0, na.rm = TRUE))
  } else "Only one completion carried risk scores, so no rank agreement is reported."
  vt <- top %>% transmute(zone = health_zone, p = .sig(p_mean, 3), sd = .sig(p_sd, 2),
                          med = .sig(rank_median, 3),
                          rr = sprintf("%d\u2013%d", as.integer(rank_min), as.integer(rank_max)),
                          at_risk = sprintf("%d of %d", as.integer(n_arms), as.integer(nC)),
                          share = .pct(share_top20, 0))
  names(vt) <- c("zone", "mean P(invasion, 2 wk) where at risk", "SD", "median rank",
                 "rank min\u2013max", "completions at risk", "top-20 share")
  .write_caption(3L,
    title = sprintf("Figure MI3. Stability of the live watch-list across M = %d onset completions.", M),
    paras = c(paste0("The operational output of the model is a ranking: the zones with the highest ",
                     "2-week invasion probability are those recommended for monitoring, and ",
                     "Figure 3B draws the top 20."),
              .imputation_text(),
              paste0("The current forecast was refitted on every completion and its 2-week ranking ",
                     "recomputed among the zones not yet invaded in that completion."),
              sprintf(paste0("**(A)** The %d zones with the best median rank: the bar spans the ",
                             "lowest and highest rank over the completions in which the zone was ",
                             "still at risk, the point is the median, and the dashed line marks ",
                             "the top-20 cut."), nrow(top)),
              sprintf(paste0("**(B)** The share of the %d completions in which the zone is inside ",
                             "that cut; a completion in which the zone had already been invaded ",
                             "counts as outside it, its invasion probability there being 0. %d of ",
                             "these %d zones are inside the cut in every completion."),
                      as.integer(nC), n_all, nrow(top)),
              rho_txt),
    values = .md_table(vt),
    sources = c("onset_mi_rank_stability.csv", "onset_mi_rank_agreement.csv",
                "onset_mi_datalayer.csv"))
  invisible(fig)
}

# =============================================================================
# MI4 — Figure 2A's curve: within-arm band vs MI total band
# =============================================================================
# The visual statement of what single imputation was hiding. The pooled detection table
# carries both components, so the two ribbons are the SAME estimator with and without the
# between-imputation term — not two different curves.
build_mi4 <- function() {
  unlink(.cap_path(4L))
  d <- .req(file.path(POOL, "key_outputs", "detection_curves.csv"), "the pooled detection curves")
  # The featured method comes from the POOLED run_info's onset_mi block, which 48 wrote from
  # what the arms actually fitted — not from a second read of model_selection.json that could
  # resolve a different field.
  feat <- MI$method
  d2 <- d %>% filter(method == feat, horizon == 2L)
  if (!nrow(d2)) stop("[mi-fig] no pooled h=2 curve for the featured model ", feat,
                      ".", call. = FALSE)
  # BOTH BANDS ARE READ, NOT COMPUTED. 48 publishes recall_within_lo/hi on the SAME t
  # critical value as the total, so this comparison isolates the between-imputation variance
  # instead of also showing the gap between two reference distributions.
  if (!all(c("recall_within_lo", "recall_within_hi", "recall_lambda") %in% names(d2)))
    stop("[mi-fig] the pooled detection table predates the within-band columns; re-run ",
         "48_onset_mi_pool.R.", call. = FALSE)
  d2 <- d2 %>% mutate(within_lo = recall_within_lo, within_hi = recall_within_hi)
  p <- ggplot(d2, aes(k, recall)) +
    geom_ribbon(aes(ymin = recall_lo, ymax = recall_hi, fill = "MI total (Rubin)"),
                alpha = 0.22, colour = NA) +
    geom_ribbon(aes(ymin = within_lo, ymax = within_hi, fill = "Within-completion only"),
                alpha = 0.35, colour = NA) +
    geom_line(colour = ACCENT, linewidth = 0.95) +
    geom_line(aes(k, recall_random), linetype = "22", colour = FAINT, linewidth = 0.7) +
    annotate("text", x = 15, y = 0.075, label = "random watch-list",
             colour = MUTED, size = 3.2, angle = 8) +
    scale_fill_manual(values = c("MI total (Rubin)" = WARM,
                                 "Within-completion only" = ACCENT), name = NULL) +
    scale_y_continuous(labels = percent_format(1), limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Zones actively monitored per round (K)", y = "Share of true invasions caught",
         # TWO LINES. At 6.4 in the one-line version ran off the right edge and the caption
         # lost its closing clause — the part naming the reference distribution.
         # lambda is the BETWEEN-COMPLETION share. It slightly overstates the onset-date share:
         # each completion's fits run on different data, so their Monte Carlo error enters it.
         subtitle = sprintf(paste0("2-week horizon; 90%% bands on Rubin's t.\n",
                                   "Median \u03bb across K: %.2f \u2014 between-completion share ",
                                   "of total variance (onset dates + fits' Monte Carlo error)."),
                            stats::median(d2$recall_lambda, na.rm = TRUE))) +
    theme_pub(11) + theme(legend.position = "top",
                          plot.subtitle = element_text(size = 8, colour = MUTED))
  save_dual(p, "FigureMI4_detection_within_vs_total", 6.4, 5.4, dir = FIG_DIR)

  # ---- caption ----
  lam <- d2$recall_lambda[is.finite(d2$recall_lambda)]
  nf <- unique(d2$n_folds)
  rounds <- if (length(nf) == 1L && is.finite(nf) && abs(nf - round(nf)) < 1e-9)
    sprintf(" over %s leave-future-out rounds", .int(nf)) else ""
  ks <- intersect(as.integer(get0("TOPK_VALUES", ifnotfound = integer(0))), d2$k)
  if (!length(ks)) ks <- range(d2$k)
  atr <- mean(d2$n_atrisk_mean, na.rm = TRUE)
  vt <- d2 %>% filter(k %in% ks) %>% arrange(k) %>%
    transmute(K = as.integer(k), recall = .pct(recall),
              within = sprintf("%s\u2013%s", .pct(within_lo), .pct(within_hi)),
              total = sprintf("%s\u2013%s", .pct(recall_lo), .pct(recall_hi)),
              lam = sprintf("%.2f", recall_lambda), df = .df_txt(recall_df),
              random = .pct(recall_random))
  names(vt) <- c("K", "recall", "within-completion 90%", "MI total 90%", "\u03bb", "df",
                 "random list")
  .write_caption(4L,
    title = paste0("Figure MI4. The operational detection curve with and without the ",
                   "onset-imputation variance."),
    paras = c(sprintf(paste0("The curve is Figure 2A's quantity \u2014 the share of true invasions ",
                             "caught had the K highest-ranked zones been monitored in each round, ",
                             "averaged over rounds \u2014 for the featured model (%s) at the 2-week ",
                             "horizon%s, recomputed on each of the M = %d completions."),
                      feat, rounds, M),
              .imputation_text(), .rubin_text,
              paste0("The line is Q\u0304. The inner band is the within-completion component ",
                     "alone (\u016a, recovered from each completion's 90% fold-bootstrap interval, ",
                     "which resamples cross-validation rounds) \u2014 what a single-imputation ",
                     "analysis reports; the outer band is the 90% interval on the total T. Both ",
                     "use the same critical value from Rubin's t, so they differ only by the ",
                     "between-completion term (1 + 1/M)B."),
              sprintf("Across K = %d\u2013%d the median \u03bb is %.2f (range %.2f\u2013%.2f).",
                      as.integer(min(d2$k)), as.integer(max(d2$k)), stats::median(lam),
                      min(lam), max(lam)),
              .lambda_text,
              sprintf(paste0("The dashed line is the expectation under a randomly ordered ",
                             "watch-list of the same length: K divided by the mean number of ",
                             "at-risk zones per round (%.1f). Values at K = %s are tabulated ",
                             "below."), atr, .and(ks))),
    values = .md_table(vt),
    sources = "pooled/key_outputs/detection_curves.csv")
  invisible(p)
}

# =============================================================================
# MI5 — cascade: within-arm interval vs MI total interval
# =============================================================================
build_mi5 <- function() {
  # A skipped figure has no caption: a stale one describing a figure that was not drawn is
  # removed with it.
  unlink(.cap_path(5L))
  # The pooled run_info records whether THIS sweep pooled a cascade; a file left by an earlier
  # sweep must not be drawn under this one's provenance.
  if (!isTRUE(MI$cascade)) {
    message("[mi-fig] this sweep pooled no cascade (ONSET_MI_CASCADE=0) — FigureMI5 skipped.")
    return(invisible(NULL))
  }
  f <- file.path(POOL, "cascade", "tables", "urban_scenario_summary.csv")
  if (!file.exists(f)) {
    message("[mi-fig] no pooled cascade summary — FigureMI5 skipped (was the sweep run with ",
            "ONSET_MI_CASCADE=0?).")
    return(invisible(NULL))
  }
  # THE HORIZON expected_added is measured at: max(CASCADE_REPORT_HORIZONS)
  # (38_urban_scenarios.R::urban_impact_metrics), which is the largest horizon in the pooled
  # reach table 48 writes beside this summary. Read, not typed.
  rc <- .req(file.path(POOL, "cascade", "tables", "cascade_reach_scores_all_zones.csv"),
             "the pooled cascade reach table")
  H_WK <- suppressWarnings(max(as.integer(rc$horizon), na.rm = TRUE))
  if (!is.finite(H_WK)) stop("[mi-fig] the pooled reach table carries no horizon.", call. = FALSE)
  d <- suppressWarnings(readr::read_csv(f, show_col_types = FALSE))
  # The primary specification: the earliest seed week, province-pooled seeded R (seed_r NA) —
  # the same slice 40_cascade_next_dominoes.R ranks its cities on.
  if ("seed_r" %in% names(d)) d <- d %>% filter(is.na(seed_r))
  d <- d %>% filter(seed_week == min(seed_week, na.rm = TRUE))
  # Same rule as MI4: the within band uses the SAME critical value as the total, so the two
  # differ only by the between-imputation term.
  if (!"expected_added_df" %in% names(d))
    stop("[mi-fig] the pooled urban summary predates the df column; re-run ",
         "48_onset_mi_pool.R.", call. = FALSE)
  d <- d %>%
    mutate(.tc = mi_tcrit(expected_added_df),
           within_lo = expected_added - .tc * expected_added_within_sd,
           within_hi = expected_added + .tc * expected_added_within_sd,
           city = factor(city, levels = city[order(expected_added)]))
  p <- ggplot(d, aes(y = city)) +
    geom_linerange(aes(xmin = expected_added_lo, xmax = expected_added_hi,
                       colour = "MI total (Rubin)"), linewidth = 1.5, alpha = 0.55) +
    geom_linerange(aes(xmin = within_lo, xmax = within_hi,
                       colour = "Within-completion only"), linewidth = 3.4, alpha = 0.8) +
    geom_point(aes(x = expected_added), size = 2.4, colour = INK) +
    scale_colour_manual(values = c("MI total (Rubin)" = WARM,
                                   "Within-completion only" = ACCENT), name = NULL) +
    labs(x = sprintf("Expected additional health zones invaded within %d weeks", H_WK), y = NULL,
         # The within interval is the Monte-Carlo SE of the expected count (sd of the draw-group
         # totals / sqrt(D), 32_cascade_simulator.R), so the MI share here is relative to that
         # precision, not to posterior uncertainty.
         subtitle = paste("Interval on the expected count: Monte-Carlo SE within completions",
                          "vs the Rubin total across completions.")) +
    theme_pub(11) + theme(legend.position = "top", panel.grid.major.y = element_blank(),
                          plot.subtitle = element_text(size = 8, colour = MUTED))
  save_dual(p, "FigureMI5_cascade_mi_band", 7.2, 4.4, dir = FIG_DIR)

  # ---- caption ----
  sw <- unique(d$seed_week)
  spec <- if ("seed_r_label" %in% names(d)) .and(unique(as.character(d$seed_r_label))) else NULL
  lam <- d$expected_added_lambda[is.finite(d$expected_added_lambda)]
  vt <- d %>% arrange(dplyr::desc(expected_added)) %>%
    transmute(city = as.character(city), added = .sig(expected_added, 3),
              within = sprintf("%s\u2013%s", .sig(within_lo, 3), .sig(within_hi, 3)),
              total = sprintf("%s\u2013%s", .sig(expected_added_lo, 3), .sig(expected_added_hi, 3)),
              lam = sprintf("%.2f", expected_added_lambda), df = .df_txt(expected_added_df))
  names(vt) <- c("city", "expected added", "within-completion 90%", "MI total 90%",
                 "\u03bb", "df")
  .write_caption(5L,
    title = paste0("Figure MI5. Three levels of variance in the cascade projection: posterior ",
                   "\u00d7 process \u00d7 onset imputation."),
    paras = c(sprintf(paste0("Expected number of additional health zones invaded within %d weeks ",
                             "if each named city were seeded in projection week %s%s, from the ",
                             "cascade re-run on each of the M = %d completions."),
                      H_WK, .and(sw), if (length(spec)) sprintf(" (seeded-city R: %s)", spec) else "",
                      M),
              .imputation_text(), .rubin_text,
              paste0("The point is Q\u0304. The thick inner bar is the within-completion component ",
                     "alone (\u016a, from each completion's standard error of the expected count), ",
                     "the variance a single-imputation cascade reports; the thin outer bar is the ",
                     "90% interval on the total T. Both use the same critical value from Rubin's ",
                     "t, so they differ only by the between-completion term (1 + 1/M)B."),
              if (length(lam))
                sprintf("Across the %d cities the median \u03bb is %.2f (range %.2f\u2013%.2f).",
                        nrow(d), stats::median(lam), min(lam), max(lam)),
              .lambda_text),
    values = .md_table(vt),
    sources = c("pooled/cascade/tables/urban_scenario_summary.csv",
                "pooled/cascade/tables/cascade_reach_scores_all_zones.csv"))
  invisible(p)
}

# -----------------------------------------------------------------------------
# THREE OUTCOMES, NOT TWO. A builder that deliberately SKIPS (MI5 without a cascade sweep)
# returns NULL without erroring; reporting that as TRUE beside genuinely drawn figures says
# a figure exists when none was written.
run <- function(f, nm) tryCatch({
    message("== ", nm, " ==")
    r <- f()
    if (is.null(r)) NA else TRUE
  }, error = function(e) { message("!! ", nm, " FAILED: ", conditionMessage(e)); FALSE })
ok <- c(MI1 = run(build_mi1, "FigureMI1"), MI2 = run(build_mi2, "FigureMI2"),
        MI3 = run(build_mi3, "FigureMI3"), MI4 = run(build_mi4, "FigureMI4"),
        MI5 = run(build_mi5, "FigureMI5"))
.lbl <- ifelse(is.na(ok), "skipped", ifelse(ok, "drawn", "FAILED"))
message(sprintf("\n[mi-fig] %s -> %s",
                paste(sprintf("%s:%s", names(ok), .lbl), collapse = "  "), FIG_DIR))
# A FAILED PANEL MUST FAIL THE SCRIPT, or the driver reports success while the previous
# run's PDFs stay on disk (the same trap make_manuscript_figures.R closes at its foot).
# A skipped one (NA) is not a failure.
.failed <- names(ok)[!is.na(ok) & !ok]
if (length(.failed))
  stop(length(.failed), " of ", length(ok), " MI figure(s) failed: ",
       paste(.failed, collapse = ", "),
       ". Their files on disk are from an earlier run.", call. = FALSE)

}   # end of the script-run guard

