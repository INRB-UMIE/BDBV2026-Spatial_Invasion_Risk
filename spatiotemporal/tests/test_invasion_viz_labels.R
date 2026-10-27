# =============================================================================
# test_invasion_viz_labels.R — what the invasion figure suite (17_invasion_viz.R) and the
# probability-scale helpers (forecast_scale.R) draw and say
#   * the featured model is always in the skill-over-time panel; sensitivity arms and the other
#     time-varying forms never fill the "top" pool
#   * the delta-stability figure labels its intervals at the level the table carries, says which
#     pool it shows, and is written to OUT_DIAGNOSTICS (it used an undefined OUT_FIGURES)
#   * the LFO tile marks the rounds that carry no fitted recalibration factor
#   * the invasion report states each probability's scale and the uncorrected rounds
#   * fs_apply_lfo_scale() refuses an expression; fs_caption() stamps the scale it resolved, and
#     save_dual()'s append keeps a caption the figure already shows
# Synthetic data only.
# =============================================================================
library(testthat)
.ST <- file.path(here::here(), "spatiotemporal")
suppressMessages(suppressWarnings({
  source(file.path(.ST, "00_config.R"))
  load_pipeline_modules("model", dir = .ST)
  source(file.path(.ST, "17_invasion_viz.R"))
  source(file.path(.ST, "forecast_scale.R"))
}))

test_that("Bayesian ensembles are coloured as ensembles", {
  expect_identical(.invasion_family(c("Bayes-ens-mean", "Bayes-ens-median", "Bayes-M14-fill-geo",
                                      "Ensemble-mean", "Distance-B1")),
                   c("Ensemble", "Ensemble", "Bayesian renewal", "Ensemble", "Baseline (naive)"))
})

.st_fake <- function() {
  ms <- c("Bayes-A-geo", paste0("Bayes-", LETTERS[2:7], "-geo"), "Bayes-A-geo-tvrw1",
          "Bayes-A-gtlong", "Bayes-A-geo-tvfrailty")
  sk <- c(0.10, 0.50, 0.49, 0.48, 0.47, 0.46, 0.45, 0.90, 0.95, 0.85)
  do.call(rbind, lapply(seq_along(ms), function(i) data.frame(
    method = ms[i], horizon = 1L, fold_id = 1:3, cutoff = as.Date("2026-07-04") + 7 * (0:2),
    n_events = c(2L, 3L, 1L), auc_pr_skill = sk[i], hit_at_k = 1L, stringsAsFactors = FALSE)))
}

test_that("skill over time always draws the featured model and keeps arms and other forms out", {
  st <- .st_fake()
  p <- plot_skill_over_time(st, featured = "Bayes-A-geo", save = FALSE)
  drawn <- levels(p$data$method)
  expect_true("Bayes-A-geo" %in% drawn)                       # lowest skill, still drawn
  expect_length(drawn, 6L)
  expect_false(any(c("Bayes-A-geo-tvrw1", "Bayes-A-gtlong", "Bayes-A-geo-tvfrailty") %in% drawn))
  expect_setequal(setdiff(drawn, "Bayes-A-geo"), paste0("Bayes-", LETTERS[2:6], "-geo"))
  # A featured FORM is drawn; the other forms stay out.
  p2 <- plot_skill_over_time(st, featured = "Bayes-A-geo-tvrw1", save = FALSE)
  expect_true("Bayes-A-geo-tvrw1" %in% levels(p2$data$method))
  expect_false("Bayes-A-geo-tvfrailty" %in% levels(p2$data$method))
  # "Prevalence-normalised" describes auc_pr_skill only.
  expect_match(p$labels$subtitle, "Prevalence-normalised")
  p3 <- plot_skill_over_time(st, metric = "hit_at_k", featured = "Bayes-A-geo", save = FALSE)
  expect_false(grepl("Prevalence-normalised", p3$labels$subtitle))
  # A featured model absent from the table is an error, not a silent omission.
  expect_error(plot_skill_over_time(st, featured = "Bayes-Z-geo", save = FALSE), "featured model")
})

test_that("the AUC-PR decomposition needs the featured model named", {
  expect_error(plot_aucpr_decomposition(.st_fake()), "needs `method`")
})

.dstab_fake <- function(z, h = 1L) {
  pf <- data.frame(method = "M", horizon = h, fold_id = 1:5,
                   cutoff = as.Date("2026-07-04") + 7 * (0:4), n_events = 2L,
                   delta = c(0.5, 0.8, 0.6, 1.1, 0.7), se_log = c(0.6, 0.5, 0.7, 0.4, 0.5))
  pf$lo <- pf$delta * exp(-z * pf$se_log); pf$hi <- pf$delta * exp(z * pf$se_log)
  sm <- data.frame(method = "M", horizon = h, n_folds = 5L, pooled_delta = 0.7,
                   pooled_lo = 0.7 * exp(-z * 0.2), pooled_hi = 0.7 * exp(z * 0.2), tau2 = 0.1,
                   I2 = 0.4, Q = 6, Q_df = 4L, Q_p = 0.2, trend_slope = 0.01, trend_p = 0.8,
                   extent_slope = 0, extent_p = 1, median_se_log = 0.5)
  list(pf = pf, sm = sm)
}

test_that("delta stability: interval level read from the table, pool named, written to diagnostics", {
  t95 <- .dstab_fake(stats::qnorm(0.975))
  p <- plot_delta_stability(t95$pf, t95$sm, horizon = 1L, method = "M", save = FALSE)
  expect_match(p[[1]]$labels$subtitle, "95%")
  expect_match(p[[1]]$labels$subtitle, "not the deployed factor")
  expect_match(p$patches$annotation$caption, "trained on")
  t90 <- .dstab_fake(stats::qnorm(0.95), h = 2L)
  p2 <- plot_delta_stability(t90$pf, t90$sm, horizon = 2L, method = "M", save = FALSE)
  expect_match(p2[[1]]$labels$x, "90% interval")
  expect_match(gsub("\n", " ", p2$patches$annotation$caption), "not the deployment pool")
  # Saved under OUT_DIAGNOSTICS, and the retained-figure gate keeps it.
  old <- OUT_DIAGNOSTICS; td <- tempfile("diag"); dir.create(td)
  assign("OUT_DIAGNOSTICS", td, envir = globalenv())
  on.exit(assign("OUT_DIAGNOSTICS", old, envir = globalenv()), add = TRUE)
  suppressMessages(plot_delta_stability(t95$pf, t95$sm, horizon = 1L, method = "M", save = TRUE))
  expect_true(file.exists(file.path(td, "delta_stability_h1.pdf")))
})

test_that("the LFO tile marks rounds without a fitted recalibration factor", {
  cut <- as.Date("2026-07-04") + 7 * (0:3)
  d <- expand.grid(health_zone = c("za", "zb"), k = 1:4, stringsAsFactors = FALSE)
  d$method <- "M"; d$horizon <- 1L; d$cutoff <- cut[d$k]; d$fold_id <- d$k
  d$p_invasion <- 0.1; d$is_new_invasion <- as.integer(d$k == 3 & d$health_zone == "za")
  d$delta_estimable <- d$k >= 3; d$p_recal <- ifelse(d$delta_estimable, 0.05, 0.1)
  p <- plot_lfo_forecast_vs_outcome(d, "M", horizon = 1L, save = FALSE)
  expect_match(p$scales$get_scales("fill")$name, "dashed = raw")
  expect_match(gsub("\n", " ", p$labels$subtitle),
               sprintf("origins before %s", format(lfo_origin(cut[3]), "%d %b")))
  d$delta_estimable <- TRUE; d$p_recal <- 0.05
  p2 <- plot_lfo_forecast_vs_outcome(d, "M", horizon = 1L, save = FALSE)
  expect_false(grepl("dashed", p2$scales$get_scales("fill")$name))
})

test_that("the invasion report names each probability scale and the uncorrected rounds", {
  ev <- data.frame(method = "M", horizon = 1L, partial_cv = FALSE, n_atrisk = 900L, n_invasions = 9L,
                   base_rate = 0.01, auc_pr_skill = 20, auc_pr_skill_lo = 12, auc_pr_skill_hi = 31,
                   auc_pr_skill_date = 19, auc_pr_skill_date_lo = 12, auc_pr_skill_date_hi = 30,
                   sel_auc_pr_skill_date = 18, sel_auc_pr_skill_date_lo = 11,
                   sel_auc_pr_skill_date_hi = 29, mean_rank_of_truth = 15, log_score = 0.03,
                   calibration_in_large = 1.5, calibration_in_large_recal = 1.2, auc_roc = 0.9,
                   n_folds_scored = 10L, sel_n_folds_delta_fitted = 6L, recall_at_15 = 0.7,
                   stringsAsFactors = FALSE)
  rs <- data.frame(health_zone = c("za", "zb"), province = "Ituri", horizon = 1L,
                   was_active_before = FALSE, p_case_invasion = c(0.2, 0.1), rr_nat = c(2, 1),
                   prob_scale = "raw")
  zw <- data.frame(health_zone = c("za", "zb", "zc"), week_start = as.Date("2026-07-04"),
                   confirmed = c(0L, 0L, 3L))
  f <- tempfile(fileext = ".md")
  suppressMessages(write_invasion_report(ev, rs, zone_week = zw, training_cutoff = as.Date("2026-07-04"),
                                         best_bayes_method = "M", out_path = f))
  L <- readLines(f)
  expect_true(any(grepl("**Training data through:** 2026-07-10", L, fixed = TRUE)))
  expect_true(any(grepl("P(case), raw", L, fixed = TRUE)))
  expect_false(any(grepl("only probability", L)))
  expect_true(any(grepl("Log score (raw, all rounds)", L, fixed = TRUE)))
  # all-round per-date skill first (the reported performance), the selection rounds' beside it
  expect_true(any(grepl("| M | 19.0 [12.0-30.0] | 18.0 [11.0-29.0] | 20.0 [12.0-31.0] |", L, fixed = TRUE)))
  expect_true(any(grepl("4 of the 10 rounds carry no fitted factor", L, fixed = TRUE)))
  expect_true(any(grepl("Averaged over the back-test rounds with an invasion", L, fixed = TRUE)))
  expect_false(any(grepl("bayes_ensemble_|model_vs_baseline", L)))
})

test_that("fs_apply_lfo_scale refuses an expression and swaps a named frame", {
  lfo <- data.frame(method = "a", p_invasion = c(0.1, 0.2), p_recal = c(0.05, 0.1),
                    was_active_before = FALSE, delta_preq = 0.5)
  expect_error(suppressMessages(fs_apply_lfo_scale(lfo[lfo$method == "a", ], scale = "recalibrated")),
               "must be a variable name")
  expect_equal(lfo$p_invasion, c(0.1, 0.2))                    # untouched
  a_long_frame_name_beyond_forty_characters_xx <- lfo
  col <- suppressMessages(fs_apply_lfo_scale(a_long_frame_name_beyond_forty_characters_xx,
                                             scale = "recalibrated"))
  expect_identical(col, "p_recal")
  expect_equal(a_long_frame_name_beyond_forty_characters_xx$p_invasion, c(0.05, 0.1))
  expect_equal(a_long_frame_name_beyond_forty_characters_xx$p_invasion_raw, c(0.1, 0.2))
})

test_that("fs_caption stamps the resolved scale; save_dual's append keeps a shown caption once", {
  old <- FS_STAMP_FIGURES
  assign("FS_STAMP_FIGURES", "always", envir = globalenv())
  on.exit(assign("FS_STAMP_FIGURES", old, envir = globalenv()), add = TRUE)
  g <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point() +
    ggplot2::labs(caption = "WINDOW") + ggplot2::theme(plot.caption = ggplot2::element_text())
  expect_match(fs_caption(g, context = "lfo", scale = "recalibrated")$labels$caption,
               "^Probabilities: recalibrated")
  # A frame that fell back to the raw column is stamped RAW, not with the requested scale.
  raw_lfo <- data.frame(method = "a", p_invasion = 0.1, was_active_before = FALSE)
  expect_match(suppressWarnings(fs_caption(g, context = "lfo", scale = "recalibrated",
                                           lfo = raw_lfo))$labels$caption, "RAW")
  # The caller (save_dual(), figure_style.R) appends through .fig_stamp_append(): the figure's
  # own caption must survive, exactly once.
  if (!exists(".fig_stamp_append", mode = "function"))
    suppressMessages(suppressWarnings(try(source(file.path(.ST, "figure_style.R")), silent = TRUE)))
  if (!exists(".fig_stamp_append", mode = "function")) testthat::skip("figure_style.R not loaded")
  cap <- .fig_stamp_append(g, function(x) fs_caption(x, context = "lfo", scale = "recalibrated"))$labels$caption
  expect_match(cap, "^WINDOW\nProbabilities: recalibrated")
  expect_equal(lengths(regmatches(cap, gregexpr("WINDOW", cap))), 1L)
  pw <- patchwork::wrap_plots(g, g) + patchwork::plot_annotation(caption = "PW")
  cap2 <- .fig_stamp_append(pw, function(x) fs_caption(x, context = "lfo", scale = "recalibrated"))$patches$annotation$caption
  expect_match(cap2, "^PW\nProbabilities: recalibrated")
  expect_equal(lengths(regmatches(cap2, gregexpr("PW", cap2))), 1L)
})
