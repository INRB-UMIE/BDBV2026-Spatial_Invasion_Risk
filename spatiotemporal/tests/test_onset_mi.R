# =============================================================================
# tests/test_onset_mi.R — onset multiple-imputation sweep
#
# Exercises the parts that can be checked WITHOUT running the sweep: the pooling
# arithmetic, the equal-weighting invariants, the D5 already-invaded rule, and — most
# importantly — the GUARDS that stop an MI pass from overwriting the published figures or
# from silently repeating one onset draw M times.
#
# The one thing these tests cannot check is acceptance test #1 of PLAN_ONSET_MI.md (arm 0
# reproducing the published evaluation), which needs a real arm. `run_onset_mi.sh` runs arm 0
# first for that reason.
# =============================================================================

suppressPackageStartupMessages({library(testthat); library(dplyr)})

.ST <- file.path(here::here(), "spatiotemporal")
source(file.path(.ST, "onset_mi_pooling.R"))

# ---------------------------------------------------------------------------
# 1. Rubin's rules
# ---------------------------------------------------------------------------
test_that("mi_rubin reproduces the textbook decomposition", {
  q <- c(0.20, 0.24, 0.22, 0.26, 0.18)
  u <- c(0.0009, 0.0010, 0.0011, 0.0009, 0.0012)
  r <- mi_rubin(q, u)
  m <- length(q)
  expect_equal(r$M, m)
  expect_equal(r$q, mean(q))
  expect_equal(r$B, var(q))
  expect_equal(r$U, mean(u))
  expect_equal(r$T, mean(u) + (1 + 1/m) * var(q))
  expect_equal(r$se, sqrt(r$T))
  # lambda is the PROPORTION OF TOTAL VARIANCE from the missing data...
  expect_equal(r$lambda, (1 + 1/m) * var(q) / r$T)
  # ...and riv is its odds form, so lambda = riv/(1+riv) must hold exactly.
  expect_equal(r$lambda, r$riv / (1 + r$riv))
  # Rubin's df, in both algebraic forms
  expect_equal(r$df, (m - 1) / r$lambda^2)
  expect_equal(r$df, (m - 1) * (1 + 1/r$riv)^2)
  # FMI is lambda PLUS the small-sample correction, hence strictly larger here
  expect_equal(r$fmi, r$lambda + (1 - r$lambda) * 2 / (r$df + 3))
  expect_gt(r$fmi, r$lambda)
  expect_true(r$u_known)
  expect_equal(r$n_u, m)
})

test_that("the fmi - lambda gap stated in onset_mi_pooling.R is the true one", {
  # gap = (1 - lambda) * 2/(nu + 3) with nu = (M - 1)/lambda^2: zero at lambda = 0 and 1 (so
  # NOT 2/22 at nu = 19, where lambda = 1), and at most 0.0146 at M = 20.
  gap <- function(l, M = 20) (1 - l) * 2 / ((M - 1) / l^2 + 3)
  o <- optimize(gap, c(0, 1), maximum = TRUE)
  expect_equal(o$objective, 0.0146, tolerance = 5e-4)
  expect_equal(o$maximum, 0.652, tolerance = 2e-3)
  expect_equal(gap(1), 0)
  r <- mi_rubin(rnorm(20))            # lambda = 1, nu = 19
  expect_equal(r$fmi, r$lambda)
  txt <- paste(readLines(file.path(.ST, "onset_mi_pooling.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("2/22 = 0.09", txt, fixed = TRUE))
  expect_true(grepl("maximum is 0.0146", txt, fixed = TRUE))
})

test_that("the interval is Rubin's t, not the normal", {
  # THE REGRESSION THIS PINS: a normal interval is materially too narrow at these M. With an
  # unknown within component lambda = 1, so df = M - 1 and t_{19,.95}/z_{.95} = 1.051.
  set.seed(4)
  q <- rnorm(20, 0.25, 0.02)
  r <- mi_rubin(q)
  expect_equal(r$df, 19)
  expect_equal(r$lambda, 1)
  expect_equal(r$hi - r$lo, 2 * qt(0.95, 19) * r$se)
  expect_gt((r$hi - r$lo) / (2 * qnorm(0.95) * r$se), 1.04)
  # ...and it collapses to the normal when the between component is negligible
  r2 <- mi_rubin(q, rep(1e6, 20))
  expect_equal((r2$hi - r2$lo) / (2 * qnorm(0.95) * r2$se), 1, tolerance = 1e-6)
})

test_that("a single completion carries no between-imputation information", {
  # B = 0 at M = 1 (var() is NA there). T must collapse to the single-draw variance, not NA:
  # a one-arm 'sweep' should reproduce the ordinary answer exactly.
  r <- mi_rubin(0.3, 0.01)
  expect_equal(r$B, 0)
  expect_equal(r$T, 0.01)
  expect_equal(r$lambda, 0)
  expect_equal(r$fmi, 0)
  expect_true(is.infinite(r$df))          # no between-imputation df from one completion
  expect_equal(r$hi - r$lo, 2 * qnorm(0.95) * sqrt(0.01))
})

test_that("an unknown within variance is flagged, not silently treated as measured zero", {
  r <- mi_rubin(c(1, 2, 3))
  expect_false(r$u_known)
  expect_equal(r$U, 0)
  expect_equal(r$lambda, 1)       # all apparent variance is between-imputation, BY CONSTRUCTION
  expect_equal(r$fmi, 1)
  expect_equal(r$df, 2)           # M - 1
  expect_equal(r$n_u, 0L)
  r2 <- mi_rubin(c(1, 2, 3), c(NA, NA, NA))
  expect_false(r2$u_known)
})

test_that("mi_rubin drops non-finite estimates together with their within variances", {
  r <- mi_rubin(c(1, NA, 3, Inf), c(0.1, 0.2, 0.3, 0.4))
  expect_equal(r$M, 2L)
  expect_equal(r$q, 2)
  expect_equal(r$U, mean(c(0.1, 0.3)))   # NOT mean(c(0.1,0.2)) — the pairing must survive
})

test_that("a length mismatch between q and u is fatal, never recycled", {
  # Recycling or subsetting a wrongly-sized u silently pairs one arm's estimate with
  # another arm's within variance, and nothing downstream could detect it.
  expect_error(mi_rubin(c(1, 2, 3), c(0.1, 0.2)), "parallel")
})

test_that("mi_u_from_ci inverts a 90% normal interval", {
  se <- 0.037
  q  <- 0.5
  expect_equal(mi_u_from_ci(q - qnorm(0.95) * se, q + qnorm(0.95) * se), se^2)
  expect_equal(mi_u_from_ci(q - qnorm(0.95) * se, q + qnorm(0.95) * se, level = 0.90), se^2)
  # the level is honoured: a 95% interval of the same se inverts to the same U at level 0.95
  expect_equal(mi_u_from_ci(q - qnorm(0.975) * se, q + qnorm(0.975) * se, level = 0.95), se^2)
  expect_error(mi_u_from_ci(0, 1, level = 90))
  expect_true(is.na(mi_u_from_ci(NA, 1)))
  expect_true(is.na(mi_u_from_ci(1, NA)))
})

test_that("the evaluation intervals 48 inverts are 90% (5%/95%) bootstrap quantiles", {
  # 48 recovers U from every <stem>_lo/_hi pair at EV_CI_LEVEL = 0.90. That is right only while
  # evaluate_invasion() forms them as the 5%/95% quantiles of its zone-cluster bootstrap.
  ev_src <- paste(readLines(file.path(.ST, "16_invasion_eval.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("ci <- function(v) if (length(v) > 10) quantile(v, c(0.05, 0.95)",
                    ev_src, fixed = TRUE))
  expect_true(grepl("Cluster bootstrap by ZONE", ev_src, fixed = TRUE))
  p48 <- paste(readLines(file.path(.ST, "48_onset_mi_pool.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("EV_CI_LEVEL <- 0.90", p48, fixed = TRUE))
  expect_true(grepl("level = EV_CI_LEVEL)", p48, fixed = TRUE))
  expect_true(grepl(".ev_u <- mi_ci_stems(names(ev_long))", p48, fixed = TRUE))
  # the old fixed list is gone, and the report no longer calls the band a fold bootstrap
  expect_false(grepl('.ev_u <- c(auc_pr = "auc_pr", log_score = "log_score")', p48, fixed = TRUE))
  expect_false(grepl("fold-cluster", p48, fixed = TRUE))
})

test_that("mi_ci_stems finds every metric with a published interval, and only those", {
  nms <- c("method", "horizon", "auc_pr", "auc_pr_lo", "auc_pr_hi", "auc_pr_skill",
           "auc_pr_skill_date", "auc_pr_skill_date_chance", "auc_pr_skill_date_lo",
           "auc_pr_skill_date_hi", "log_score", "log_score_lo", "log_score_hi",
           "orphan_lo", "half_lo", "half")
  expect_setequal(mi_ci_stems(nms), c("auc_pr", "auc_pr_skill_date", "log_score"))
  expect_length(mi_ci_stems(c("a", "b")), 0L)
})

test_that("mi_pool_by works without dplyr attached", {
  # The helper is sourced by callers that do not attach dplyr; a bare %>% made it depend on
  # the caller's search path.
  out <- callr_ok <- NULL
  src <- file.path(.ST, "onset_mi_pooling.R")
  res <- suppressWarnings(system2("Rscript", c("-e", shQuote(sprintf(
    "source('%s'); d <- data.frame(k=c(1,1,2,2), a=1:4, v=c(1,2,3,4), uu=rep(0.01,4)); o <- mi_pool_by(d,'k','v','uu'); cat(nrow(o), round(o$q,4))", src))),
    stdout = TRUE, stderr = TRUE))
  expect_true(is.null(attr(res, "status")) || attr(res, "status") == 0L,
              info = paste(res, collapse = " | "))
  expect_true(any(grepl("^2 1.5 3.5", res)), info = paste(res, collapse = " | "))
})

test_that("the MI total variance is never below the within variance", {
  # Mechanical, but it is the property FigureMI4 draws: the MI ribbon must contain the
  # within-completion ribbon at every K, and acceptance test #6 checks it on real output.
  set.seed(1)
  for (i in 1:50) {
    q <- rnorm(8); u <- runif(8, 0, 0.5)
    r <- mi_rubin(q, u)
    expect_gte(r$T, r$U)
  }
})

test_that("mi_pool_by pools per key and counts arms", {
  df <- expand.grid(k = 1:3, arm = 1:4) %>%
    mutate(v = k + arm / 10, uu = 0.01)
  out <- mi_pool_by(df, "k", "v", "uu")
  expect_equal(nrow(out), 3L)
  expect_true(all(out$n_arms == 4L))
  expect_equal(out$q[out$k == 2], mean(2 + (1:4) / 10))
  expect_true(all(c("lambda", "fmi", "df", "riv", "n_u") %in% names(out)))
  expect_true(all(out$fmi >= out$lambda, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# 2. The D5 already-invaded rule
# ---------------------------------------------------------------------------
test_that("mi_zero_affected sets p to exactly zero for that completion's invaded zones", {
  zones <- c("A", "B", "C")
  drl <- list(h1 = list(zones = zones, horizon = 1,
                        cum = matrix(runif(9, 0.1, 0.5), nrow = 3, dimnames = list(NULL, zones)),
                        p   = matrix(runif(9, 0.1, 0.5), nrow = 3, dimnames = list(NULL, zones))))
  out <- mi_zero_affected(drl, "B")
  expect_true(all(out$h1$cum[, 2] == 0))
  expect_true(all(out$h1$p[, 2] == 0))
  # p = 1 - exp(-cum) must still hold at the zeroed column
  expect_equal(unname(out$h1$p[, 2]), 1 - exp(-unname(out$h1$cum[, 2])))
  # untouched elsewhere
  expect_equal(out$h1$cum[, c(1, 3)], drl$h1$cum[, c(1, 3)])
  expect_identical(mi_zero_affected(drl, character(0)), drl)
})

# ---------------------------------------------------------------------------
# 3. Equal weighting — the requirement that makes it multiple imputation
# ---------------------------------------------------------------------------
test_that("mi_weight gives every completion equal weight regardless of rows scored", {
  # Arms score different numbers of cells because the at-risk support moves with the imputed
  # labels. A plain mean over the stacked frame would weight arm 2 four times as heavily here.
  # The weights come from the pipeline's own function (onset_mi_pooling.R), one slice.
  M <- 2L
  n <- c(100L, 400L)
  d <- dplyr::bind_rows(lapply(seq_len(M), function(i)
    tibble::tibble(arm = i, err2 = rep(i * 1.0, n[i]))))
  d$mi_weight <- mi_slice_weights(rep(TRUE, nrow(d)), d$arm, list(rep("F", nrow(d))), M)
  expect_equal(sum(d$mi_weight), 1)
  expect_equal(weighted.mean(d$err2, d$mi_weight), 1.5)   # equal weight per arm
  expect_false(isTRUE(all.equal(mean(d$err2), 1.5)))      # row-weighted would be 1.8
})

test_that("mi_slice_weights equalises the arms WITHIN each (method, horizon) slice", {
  # A row-level consumer averages within (method, horizon). Normalising over an arm's whole
  # frame equalised the arms only across the stack.
  M <- 2L
  d <- dplyr::bind_rows(
    tibble::tibble(arm = 1L, method = "F", horizon = 1L, e = 1, sc = TRUE, .rep = 1:10),
    tibble::tibble(arm = 1L, method = "F", horizon = 2L, e = 1, sc = TRUE, .rep = 1:90),
    tibble::tibble(arm = 2L, method = "F", horizon = 1L, e = 3, sc = TRUE, .rep = 1:50),
    tibble::tibble(arm = 2L, method = "F", horizon = 2L, e = 3, sc = TRUE, .rep = 1:50),
    tibble::tibble(arm = 2L, method = "F", horizon = 2L, e = 9, sc = FALSE, .rep = 1:5))
  d$w <- mi_slice_weights(d$sc, d$arm, list(d$method, d$horizon), M)
  expect_true(all(is.na(d$w[!d$sc])))
  s <- d[d$sc, ] %>% dplyr::group_by(arm, method, horizon) %>%
    dplyr::summarise(w = sum(w), .groups = "drop")
  expect_true(all(abs(s$w - 1 / M) < 1e-12))
  for (h in 1:2) {
    x <- d[d$sc & d$horizon == h, ]
    expect_equal(weighted.mean(x$e, x$w), 2)            # mean of the two arms' slice means
  }
  # The old whole-frame weight: both arms score 100 rows in total, so every row got
  # 1/(2*100) and the h = 1 slice (10 rows vs 50) averaged to 2.67, not 2.
  x <- d[d$sc & d$horizon == 1L, ]
  expect_equal(weighted.mean(x$e, rep(1 / (M * 100), nrow(x))), (10 * 1 + 50 * 3) / 60)
  # 48 builds the weight with it, on the consumer's own filter and slice
  p48 <- paste(readLines(file.path(.ST, "48_onset_mi_pool.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("lfo_pooled$mi_weight <- mi_slice_weights(", p48, fixed = TRUE))
  expect_true(grepl("slice = list(lfo_pooled$method, lfo_pooled$horizon)", p48, fixed = TRUE))
})

test_that("mi_deploy_scale refuses a mixture of recalibrated and raw completions", {
  d <- c(`1` = 0.5, `2` = 0.4)
  expect_identical(mi_deploy_scale(list(a = d, b = d), 1:2), "recalibrated-deployed-mi")
  expect_identical(mi_deploy_scale(list(a = NULL, b = NULL), 1:2), "raw-mi")
  expect_identical(mi_deploy_scale(list(a = c(h1 = .5, h2 = .4), b = 0.5), 1:2),
                   "recalibrated-deployed-mi")       # "h<h>" names; one unnamed scalar
  expect_error(mi_deploy_scale(list(a = d, b = NULL), 1:2), "arm\\(s\\) b")
  expect_error(mi_deploy_scale(list(a = d, b = c(`1` = 0.5)), 1:2), "arm\\(s\\) b")
  expect_error(mi_deploy_scale(list(a = d, b = c(`1` = 0.5, `2` = NA)), 1:2), "arm\\(s\\) b")
  expect_error(mi_deploy_scale(list(a = d, b = c(`1` = 0.5, `2` = 0)), 1:2), "arm\\(s\\) b")
  p48 <- paste(readLines(file.path(.ST, "48_onset_mi_pool.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("rs_pooled$prob_scale <- POOL_PROB_SCALE", p48, fixed = TRUE))
  expect_false(grepl('rs_pooled$prob_scale <- "recalibrated-deployed-mi"', p48, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 4. The guards
# ---------------------------------------------------------------------------
test_that("RANDOM_SEED is honoured from the environment and a typo is fatal", {
  # This is what makes the replicates real: 00_config.R is re-sourced by five modules and by
  # run_cascade.R, so a post-source() override cannot survive a subprocess chain.
  rs <- function(env) {
    out <- suppressWarnings(system2("Rscript",
      c("-e", shQuote("suppressMessages(source(file.path(here::here(),'spatiotemporal','00_config.R'))); cat(RANDOM_SEED)")),
      env = env, stdout = TRUE, stderr = FALSE))
    list(out = paste(out, collapse = ""), status = attr(out, "status") %||% 0L)
  }
  expect_equal(rs("RANDOM_SEED=20260711")$out, "20260711")
  expect_equal(rs("RANDOM_SEED=")$out, "20260704")
  bad <- rs("RANDOM_SEED=not-an-integer")
  expect_true(bad$status != 0L)   # must FAIL, not fall back to the default
})

test_that("ONSET_IMPUTE_SEED is honoured from the environment and a non-integer is fatal", {
  # The MI arms vary THIS seed (47 sets it before sourcing 00_config.R), so it must survive
  # every re-source through the environment, default to RANDOM_SEED when unset, and refuse a
  # value it cannot read rather than fall back to the default (which would repeat arm 0's draw).
  os <- function(env) {
    out <- suppressWarnings(system2("Rscript",
      c("-e", shQuote("suppressMessages(source(file.path(here::here(),'spatiotemporal','00_config.R'))); cat(ONSET_IMPUTE_SEED, RANDOM_SEED, is.integer(ONSET_IMPUTE_SEED))")),
      env = env, stdout = TRUE, stderr = FALSE))
    list(out = paste(out, collapse = ""), status = attr(out, "status") %||% 0L)
  }
  expect_equal(os(c("ONSET_IMPUTE_SEED=20260711", "RANDOM_SEED="))$out, "20260711 20260704 TRUE")
  # unset -> the model seed, so an ordinary run is unchanged
  expect_equal(os(c("ONSET_IMPUTE_SEED=", "RANDOM_SEED=20260801"))$out, "20260801 20260801 TRUE")
  bad <- os(c("ONSET_IMPUTE_SEED=not-an-integer", "RANDOM_SEED="))
  expect_true(bad$status != 0L)   # must FAIL, not fall back to the default
})

test_that("ONSET_IMPUTE_SEED rejects a fractional value (strict integer)", {
  # as.integer("1.5") is 1L, not NA, so a fractional seed would be silently TRUNCATED — two
  # different requested completions could then draw the same onsets. 00_config.R must refuse it.
  out <- suppressWarnings(system2("Rscript",
    c("-e", shQuote("suppressMessages(source(file.path(here::here(),'spatiotemporal','00_config.R'))); cat(ONSET_IMPUTE_SEED)")),
    env = c("ONSET_IMPUTE_SEED=20260704.5", "RANDOM_SEED="), stdout = TRUE, stderr = FALSE))
  expect_true((attr(out, "status") %||% 0L) != 0L,
              info = paste("00_config.R accepted ONSET_IMPUTE_SEED=20260704.5 as", paste(out, collapse = "")))
})

test_that("the LFO cache stamp fingerprints the seed", {
  # Without this, REUSE_CACHED_LFO=1 serves one completion's cross-validation to every arm and
  # the sweep reports a zero-width interval that looks like precision.
  src <- readLines(file.path(.ST, "run_all.R"), warn = FALSE)
  i <- grep("^\\.lfo_stamp <- list\\(", src)
  expect_length(i, 1L)
  blk <- paste(src[i:(i + 90L)], collapse = " ")
  expect_true(grepl("random_seed", blk, fixed = TRUE))
  expect_true(grepl("onset_impute_stratum", blk, fixed = TRUE))
  # PARSE the assignment (grow the line range until it is one complete expression) and check
  # the stamp's own FIELDS, not a substring that a comment could satisfy: the onset seed is
  # what an MI arm varies, so without `onset_impute_seed` every arm would hit one cache entry.
  ex <- NULL
  for (j in i:min(length(src), i + 300L)) {
    ex <- tryCatch(parse(text = src[i:j]), error = function(e) NULL)
    if (!is.null(ex)) break
  }
  expect_false(is.null(ex))
  call <- ex[[1]]
  expect_identical(as.character(call[[2]]), ".lfo_stamp")
  flds <- as.list(call[[3]])[-1]
  expect_true(all(c("random_seed", "onset_impute_seed") %in% names(flds)))
  expect_true(grepl("ONSET_IMPUTE_SEED", paste(deparse(flds$onset_impute_seed), collapse = ""),
                    fixed = TRUE))
  expect_true(grepl("RANDOM_SEED", paste(deparse(flds$random_seed), collapse = ""), fixed = TRUE))
})

test_that("the figure suites refuse pooled inputs without a figure subdirectory", {
  # The single worst thing an MI pass can do is overwrite the PUBLISHED Figures 2/3/4 with
  # pooled twins carrying the same basenames. Both scripts must refuse outright.
  for (f in c("make_manuscript_figures.R", "40_cascade_next_dominoes.R")) {
    r <- suppressWarnings(system2("Rscript", c(file.path(.ST, f)),
                                  env = c("MANU_IN_ROOT=/tmp", "MANU_FIG_SUB="),
                                  stdout = TRUE, stderr = TRUE))
    st <- attr(r, "status")
    expect_true(!is.null(st) && st != 0L,
                info = paste(f, "must exit non-zero when MANU_IN_ROOT is set without MANU_FIG_SUB"))
    expect_true(any(grepl("MANU_FIG_SUB", r)),
                info = paste(f, "must name MANU_FIG_SUB in its refusal"))
  }
})

test_that("the MI diagnostic figures are on the retained-figure allow-list", {
  # save_dual() silently DROPS any figure whose stem is not on FIGURE_KEEP, and the retained-
  # figure sweep deletes strays: an MI figure missing from the list would be built and thrown
  # away with no error anywhere.
  source(file.path(.ST, "00_config.R"), local = TRUE)
  for (stem in c("FigureMI1_skill_across_imputations", "FigureMI2_invasion_week_shift",
                 "FigureMI3_rank_stability", "FigureMI4_detection_within_vs_total",
                 "FigureMI5_cascade_mi_band"))
    expect_true(figure_is_kept(stem), info = stem)
  # The pooled Figures 2/3/4 need no entry of their own — they keep their published stems —
  # but they must survive the gate when written under the mi/ subdirectory.
  for (stem in c("Figure2", "Figure3", "Figure4"))
    expect_true(figure_is_kept(file.path("mi", stem)), info = stem)
})

test_that("CASCADE_OUT_ROOT sandboxes every cascade product", {
  # Without it, each of M cascade arms overwrites the published cascade tables and key_outputs
  # bundle in turn, and the canonical tree ends up holding whichever arm finished last.
  td <- file.path(tempdir(), "mi_cascade_root")
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  # 30_projection_config.R normalizePath()s the root, and on macOS tempdir() lives under a
  # /var -> /private/var symlink, so the returned path is resolved. Compare resolved paths.
  td <- normalizePath(td, mustWork = TRUE)
  out <- suppressWarnings(system2("Rscript",
    c("-e", shQuote(paste0("suppressMessages({source(file.path(here::here(),'spatiotemporal','00_config.R'));",
                           "source(file.path(here::here(),'spatiotemporal','30_projection_config.R'))});",
                           "cat(OUT_CASCADE)"))),
    env = paste0("CASCADE_OUT_ROOT=", td), stdout = TRUE, stderr = FALSE))
  expect_true(any(grepl(td, out, fixed = TRUE)),
              info = paste("OUT_CASCADE did not follow CASCADE_OUT_ROOT:", paste(out, collapse = " ")))
})

test_that("run_cascade.R gates its figure scripts on CASCADE_RUN_FIGURES", {
  # 40_cascade_next_dominoes.R DELETES the published Figure4 when a panel fails, so an arm
  # must never reach it.
  src <- paste(readLines(file.path(.ST, "run_cascade.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("CASCADE_RUN_FIGURES", src, fixed = TRUE))
  i <- regexpr("\\.run_cascade_figs", src)
  expect_true(i > 0)
  # both figure scripts must sit INSIDE the gate
  blk <- substring(src, i, nchar(src))
  gate_end <- regexpr("\\n\\} else \\{", blk)
  expect_true(gate_end > 0)
  inside <- substring(blk, 1, gate_end)
  expect_true(grepl("37_cascade_figure4.R", inside, fixed = TRUE))
  expect_true(grepl("40_cascade_next_dominoes.R", inside, fixed = TRUE))
})

test_that("the arm runner pins the seeds and the stratum before sourcing anything", {
  # TWO SEEDS (00_config.R): the arm varies ONSET_IMPUTE_SEED, the onset draw's only seed, and
  # holds RANDOM_SEED (Stan, EpiNow2, bootstraps) at the published run's value.
  src <- readLines(file.path(.ST, "47_onset_mi_arm.R"), warn = FALSE)
  i_oseed <- grep("Sys.setenv\\(ONSET_IMPUTE_SEED = as.character\\(MI_SEED\\)\\)", src)
  i_rseed <- grep("Sys.setenv\\(RANDOM_SEED|Sys.unsetenv\\(\"RANDOM_SEED\"\\)", src)
  i_str   <- grep("Sys.setenv\\(ONSET_IMPUTE_STRATUM", src)
  i_src   <- grep('source\\(file.path\\(ST_DIR, "00_config.R"\\)\\)', src)
  expect_length(i_src, 1L)
  expect_true(length(i_oseed) == 1L && i_oseed < i_src)
  expect_true(length(i_rseed) >= 1L && max(i_rseed) < i_src)
  expect_true(length(i_str)  >= 1L && min(i_str)  < i_src)
  # the per-arm seed is NOT the model seed any more
  expect_false(any(grepl("Sys.setenv(RANDOM_SEED = as.character(MI_SEED))", src, fixed = TRUE)))
  # ...asserts afterwards that the override actually stuck, and records both seeds
  expect_true(any(grepl("ONSET_IMPUTE_SEED is %s after sourcing", src, fixed = TRUE)))
  txt <- paste(src, collapse = "\n")
  expect_true(grepl("onset_impute_seed = as.integer(ONSET_IMPUTE_SEED), random_seed = as.integer(RANDOM_SEED)",
                    txt, fixed = TRUE))
  # the cascade subprocess gets both, the model seed unvaried
  expect_true(grepl("RANDOM_SEED = as.character(RANDOM_SEED),", txt, fixed = TRUE))
  expect_true(grepl("ONSET_IMPUTE_SEED = as.character(MI_SEED),", txt, fixed = TRUE))
})

test_that("the pool requires one model seed and distinct onset seeds", {
  txt <- paste(readLines(file.path(.ST, "48_onset_mi_pool.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl('POOL_RANDOM_SEED <- .one(function(a) a$random_seed', txt, fixed = TRUE))
  expect_true(grepl("as.integer(a$onset_impute_seed)", txt, fixed = TRUE))
  expect_true(grepl('"random_seed", "onset_impute_seed"', txt, fixed = TRUE))
})

test_that("43's imputation replicates vary the onset seed, not the model seed", {
  txt <- paste(readLines(file.path(.ST, "43_spread_kinematics.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl('assign("ONSET_IMPUTE_SEED", ONSET_SEED_BASE + r, envir = globalenv())',
                    txt, fixed = TRUE))
  expect_false(grepl('assign("RANDOM_SEED", RANDOM_SEED_BASE + r', txt, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 5. The arm deploys what run_all.R deploys; the pool refuses arms of another run
# ---------------------------------------------------------------------------
test_that("the arm recalibrates its live forecast with delta's Gamma posterior, as run_all does", {
  # Selecting c("method", "horizon", "delta") alone dropped gamma_shape/gamma_rate, so every arm
  # FIXED delta while the published forecast draws it per posterior draw — and arm 0 no longer
  # reproduced the published table.
  src <- paste(readLines(file.path(.ST, "47_onset_mi_arm.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl('recal_tbl[, c("method", "horizon", "delta")]', src, fixed = TRUE))
  expect_true(grepl('intersect(c("method", "horizon", "delta", "gamma_shape",', src, fixed = TRUE))
})

test_that("the arm's run_all guard fails when run_all stops deploying the Gamma columns", {
  src <- readLines(file.path(.ST, "47_onset_mi_arm.R"), warn = FALSE)
  a <- grep("^\\.assert_matches_run_all <- function\\(\\) \\{", src)
  b <- grep("^\\.assert_matches_run_all\\(\\)$", src)
  expect_length(a, 1L); expect_length(b, 1L)
  env <- new.env()
  eval(parse(text = src[a:(b - 1L)]), envir = env)
  ra <- readLines(file.path(.ST, "run_all.R"), warn = FALSE)
  env$ST_DIR <- .ST
  expect_silent(env$.assert_matches_run_all())          # the live run_all.R passes
  td <- tempfile("ra"); dir.create(td)
  ra2 <- sub('"gamma_rate"), names(recal_tbl))]', '"zz"), names(recal_tbl))]', ra, fixed = TRUE)
  expect_false(identical(ra, ra2))
  writeLines(ra2, file.path(td, "run_all.R"))
  env$ST_DIR <- td
  expect_error(env$.assert_matches_run_all(), "gamma_rate")
})

test_that("the arm resolves run_all.R's featured spec (a time-varying form included) and fits it alone", {
  src <- paste(readLines(file.path(.ST, "47_onset_mi_arm.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("bayes_default_grid(mobility_matrices, gt_arms = FALSE)", src, fixed = TRUE))
  expect_false(grepl("bayes_default_grid(mobility_matrices))", src, fixed = TRUE))
  # A featured form is its own spec (stage 2), never its fixed-beta parent.
  expect_true(grepl("bayes_spec_for_label(FEATURED_METHOD", src, fixed = TRUE))
  expect_true(grepl("bayes_tv_form_specs(.s1)", src, fixed = TRUE))
  # The back-test model comes from the suite's one factory, and the live fit is on the model's
  # OWN full design, as fit_bayes_suite() fits it (no grid-wide row alignment).
  expect_true(grepl("make_suite_lfo_model(", src, fixed = TRUE))
  expect_false(grepl("common_keys", src, fixed = TRUE))
  # The cascade arm pins the form as it pins the kernel.
  expect_true(grepl("CASCADE_FORM = sp$form", src, fixed = TRUE))
  # and the featured model through the suite's one reader
  expect_true(grepl("read_featured_model(", src, fixed = TRUE))
})

test_that("the pool refuses arms of another analysis date, featured model or data definition", {
  src <- readLines(file.path(.ST, "48_onset_mi_pool.R"), warn = FALSE)
  txt <- paste(src, collapse = "\n")
  for (pat in c("POOL_ADATE", ".ri_pub$analysis_date", "PUB_FEATURED <- read_featured_model(",
                "identical(POOL_METHOD, PUB_FEATURED)", "nondefault_data_definition()",
                "a$linelist_snapshot", "a$onset_impute"))
    expect_true(grepl(pat, txt, fixed = TRUE), info = pat)
  # ...and every refusal precedes the first file the pool writes.
  first_write <- min(grep("write_csv\\(|saveRDS\\(|write_json\\(|writeLines\\(|file.copy\\(", src))
  last_check  <- max(grep("\\.prov_fail\\(", src))
  expect_lt(last_check, first_write)
})

test_that("the pooled cascade carries a LOGICAL `elevated` test and integer ranks", {
  # 38/40 read `elevated %in% TRUE`; a share there admitted only zones elevated in EVERY arm.
  txt <- paste(readLines(file.path(.ST, "48_onset_mi_pool.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("elevated   = mean(elevated %in% TRUE)", txt, fixed = TRUE))
  expect_true(grepl("elevated_share = mean(elevated %in% TRUE)", txt, fixed = TRUE))
  expect_true(grepl("elevated = is.finite(delta_lo) & delta_lo > 0", txt, fixed = TRUE))
  expect_true(grepl("rr_nat_rank = dplyr::min_rank(", txt, fixed = TRUE))
  expect_true(grepl("priority_rank = dplyr::min_rank(", txt, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 6. The pool's provenance, cascade mirror and rank table; the MI figure suite
# ---------------------------------------------------------------------------
test_that("the pool resolves the model seed like 47 and checks seeds and draws files", {
  txt <- paste(readLines(file.path(.ST, "48_onset_mi_pool.R"), warn = FALSE), collapse = "\n")
  # no restated literal: the published seed, else 00_config.R's own default
  expect_false(grepl('Sys.setenv(RANDOM_SEED = "20260704")', txt, fixed = TRUE))
  expect_true(grepl('Sys.unsetenv("RANDOM_SEED")', txt, fixed = TRUE))
  expect_true(grepl('get0("RANDOM_SEED", ifnotfound = NA_integer_)', txt, fixed = TRUE))
  expect_true(grepl("identical(.prs, EXPECTED_MODEL_SEED)", txt, fixed = TRUE))
  # a missing model seed is a refusal, not a skipped comparison
  expect_true(grepl("is.na(.prs)", txt, fixed = TRUE))
  # arm 0 must be the published completion when run_info records its onset seed
  expect_true(grepl(".ri_pub$onset_imputation$onset_impute_seed", txt, fixed = TRUE))
  # each draws file matches its arm file on the fingerprint and both seeds
  expect_true(grepl('c("fingerprint", "arm", "seed", "onset_impute_seed", "random_seed")', txt,
                    fixed = TRUE))
})

test_that("the pooled cascade mirror copies the hub table 40 actually reads, or stops", {
  txt <- paste(readLines(file.path(.ST, "48_onset_mi_pool.R"), warn = FALSE), collapse = "\n")
  # the reader and the writer agree on the name...
  r40 <- paste(readLines(file.path(.ST, "40_cascade_next_dominoes.R"), warn = FALSE), collapse = "\n")
  w38 <- paste(readLines(file.path(.ST, "38_urban_scenarios.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl('file.path(tbl_dir, "urban_hub_selection.csv")', r40, fixed = TRUE))
  expect_true(grepl('file = "urban_hub_selection.csv"', w38, fixed = TRUE))
  # ...and 48 copies that file, refusing (not skipping) when it is missing
  expect_true(grepl('for (rel in "urban_hub_selection.csv")', txt, fixed = TRUE))
  # (the retired names, which no script writes, survive only in the explanatory comment)
  expect_false(grepl('for (rel in c("urban_hub_zones.csv", "urban_seed_zones.csv"))', txt, fixed = TRUE))
  expect_false(grepl("if (!length(fs)) next", txt, fixed = TRUE))
})

test_that("the rank-stability share is over completions, as FigureMI3's axis says", {
  txt <- paste(readLines(file.path(.ST, "48_onset_mi_pool.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("share_top20 = sum(rank <= 20L) / .n_rk", txt, fixed = TRUE))
  expect_false(grepl("share_top20 = mean(rank <= 20L)", txt, fixed = TRUE))
  f49 <- paste(readLines(file.path(.ST, "49_onset_mi_figures.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("Share of completions in which the zone is top-20", f49, fixed = TRUE))
})

test_that("the MI figure suite stamps its figures, draws the between component of T, writes captions", {
  src <- readLines(file.path(.ST, "49_onset_mi_figures.R"), warn = FALSE)
  txt <- paste(src, collapse = "\n")
  i_stamp <- grep("^MI_STAMP_TXT <- ", src)
  i_run   <- grep("^ok <- c\\(MI1 = run\\(build_mi1", src)
  expect_length(i_stamp, 1L); expect_length(i_run, 1L)
  expect_lt(i_stamp, i_run)
  expect_false(grepl("SEED0", txt, fixed = TRUE))
  # the thick bar is sqrt(T - U) = sqrt((1 + 1/M)B), not sqrt(B)
  expect_true(grepl("ymin = q - sqrt(pmax(T - U, 0)), ymax = q + sqrt(pmax(T - U, 0))", txt, fixed = TRUE))
  expect_false(grepl("ymin = q - sqrt(pmax(B, 0))", txt, fixed = TRUE))
  # every builder writes (or, failing/skipping, removes) its own caption
  for (n in 1:5) {
    expect_true(grepl(sprintf("unlink(.cap_path(%dL))", n), txt, fixed = TRUE), info = n)
    expect_true(grepl(sprintf(".write_caption(%dL,", n), txt, fixed = TRUE), info = n)
  }
  expect_false(grepl("fold-cluster", txt, fixed = TRUE))
  # no typed dates or seeds in the caption writer
  expect_false(grepl("2026-", txt, fixed = TRUE))
  expect_false(grepl("20260704", txt, fixed = TRUE))
})
