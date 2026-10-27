# =============================================================================
# tests/test_mobility_value_rank.R — the zone-cluster bootstrap of the MEAN RANK OF TRUTH
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# FigureS25/FigureS26 quote a 90% interval on the mean rank of an invaded zone. That interval
# is not computed by re-ranking every resample: mv_boot_mean_rank() exploits the fact that a
# rank within a forecast origin is LINEAR in the zones' resample multiplicities and gets all
# 10,000 resamples from one matrix product (50_mobility_value.R). The saving is real — a loop
# would re-rank every scored row once per resample, method and horizon — but the identity has
# to hold exactly, including the
# parts that are easy to get wrong and impossible to see in the figure:
#
#   * TIES. The published estimator uses rank(ties.method = "average"), so the fast form must
#     add half the tied mass, and must compare scores with `==` rather than a tolerance.
#   * MULTIPLICITY. A zone drawn twice contributes its rank twice to its origin's mean.
#   * DROPPED TRUTHS. A resample that retains no invasion is not a draw of zero; it is not a
#     draw at all, and an origin whose truths were all dropped must not pull its origin mean.
#   * ORIGIN LOOKUP. mv_rank_setup() indexes origins positionally, which is only safe because
#     it splits on a factor with every origin as a level.
#
# Each is pinned below against a brute-force loop over mv_mean_rank(), the estimator
# 16_invasion_eval.R publishes.
# =============================================================================

suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

# 50_mobility_value.R is not in run_tests.R's live-module list (it is a standalone appendix
# module). Source it here; the run itself is guarded by is_script_run(), so nothing executes.
# suppressWarnings() only around the LOAD: attaching dplyr/ggplot2 here re-emits the
# "package built under R x.y.z" notes the suite has already printed once, and they would be
# reported against this file as failures-in-waiting. Nothing in the tests themselves is
# suppressed.
if (!exists("mv_boot_mean_rank", mode = "function"))
  try(suppressWarnings(source(file.path(here::here(), "spatiotemporal",
                                        "50_mobility_value.R"))), silent = TRUE)

.mvr_fns <- function()
  all(vapply(c("mv_boot_mean_rank", "mv_rank_setup", "mv_mean_rank"), exists,
             logical(1), mode = "function"))

# A small support with everything that makes the estimator awkward: origins of unequal size
# (three cells are missing), scores rounded to one decimal so exact ties are everywhere, one
# origin carrying two invasions and one carrying none.
.mvr_fixture <- function(seed = 7) {
  set.seed(seed)
  zones <- sprintf("Z%02d", 1:9); folds <- sprintf("F%d", 1:5)
  g <- expand.grid(fold_id = folds, health_zone = zones, stringsAsFactors = FALSE)
  g <- g[-c(3, 11, 27), ]
  g$p_invasion <- round(runif(nrow(g)), 1)
  g$is_new_invasion <- 0L
  g$is_new_invasion[c(1, 6, 14, 22, 30)] <- 1L
  cells <- g %>% mutate(cell = paste(fold_id, health_zone, sep = "\r")) %>% arrange(cell)
  zl <- split(seq_len(nrow(cells)), cells$health_zone)
  list(cells = cells, zl = zl, nz = length(zl),
       st = mv_rank_setup(match(cells$fold_id, sort(unique(cells$fold_id))),
                          match(cells$health_zone, names(zl)),
                          cells$is_new_invasion, length(zl)))
}

# Resample multiplicities, the (zones x resamples) form the fast estimator consumes.
.mvr_M <- function(picks, nz) {
  M <- apply(picks, 2L, tabulate, nbins = nz); storage.mode(M) <- "double"; M
}

# The estimator the figure would use if speed were free.
.mvr_brute <- function(cells, zl, picks) {
  vapply(seq_len(ncol(picks)), function(b) {
    g <- cells[unlist(zl[picks[, b]], use.names = FALSE), ]
    if (sum(g$is_new_invasion) == 0L) return(NA_real_)
    mv_mean_rank(g)
  }, numeric(1))
}

test_that("at unit multiplicities the fast form IS the published estimator", {
  if (!.mvr_fns()) testthat::skip("50_mobility_value.R not sourced")
  f <- .mvr_fixture()
  one <- mv_boot_mean_rank(f$cells$p_invasion, f$st, matrix(1, f$nz, 1))[1]
  expect_equal(one, mv_mean_rank(f$cells), tolerance = 1e-12)
})

test_that("every resample matches a brute-force re-ranking, ties and all", {
  if (!.mvr_fns()) testthat::skip("50_mobility_value.R not sourced")
  f <- .mvr_fixture()
  # The fixture rounds scores to one decimal precisely so this is not a tie-free check.
  expect_true(any(duplicated(f$cells$p_invasion)))
  set.seed(101)
  picks <- matrix(sample.int(f$nz, f$nz * 200L, replace = TRUE), nrow = f$nz)
  fast <- mv_boot_mean_rank(f$cells$p_invasion, f$st, .mvr_M(picks, f$nz))
  slow <- .mvr_brute(f$cells, f$zl, picks)
  expect_equal(fast, slow, tolerance = 1e-10)
})

test_that("a resample that retains no invasion is dropped, not scored as zero", {
  if (!.mvr_fns()) testthat::skip("50_mobility_value.R not sourced")
  f <- .mvr_fixture()
  # One truth, in one zone: roughly a third of resamples will not contain it at all.
  cells <- f$cells; cells$is_new_invasion <- 0L
  cells$is_new_invasion[which(cells$health_zone == "Z04")[1]] <- 1L
  st <- mv_rank_setup(match(cells$fold_id, sort(unique(cells$fold_id))),
                      match(cells$health_zone, names(f$zl)), cells$is_new_invasion, f$nz)
  set.seed(102)
  picks <- matrix(sample.int(f$nz, f$nz * 200L, replace = TRUE), nrow = f$nz)
  fast <- mv_boot_mean_rank(cells$p_invasion, st, .mvr_M(picks, f$nz))
  slow <- .mvr_brute(cells, f$zl, picks)
  expect_gt(sum(is.na(slow)), 0)                    # the branch is actually exercised
  expect_identical(is.na(fast), is.na(slow))        # ... and both agree on which draws it is
  expect_equal(fast, slow, tolerance = 1e-10)
})

test_that("a duplicated zone counts twice, and its rank stays tie-averaged", {
  if (!.mvr_fns()) testthat::skip("50_mobility_value.R not sourced")
  f <- .mvr_fixture()
  # Every zone drawn three times: each origin's candidate set triples, so a truth's rank moves
  # from r to 3r - 1 (3(r - 1) zones now strictly above it, plus the middle of its own tied
  # triple). The mean over origins inherits the same map, which a brute-force run confirms.
  M3 <- matrix(3, f$nz, 1)
  trip <- mv_boot_mean_rank(f$cells$p_invasion, f$st, M3)[1]
  picks <- matrix(rep(seq_len(f$nz), each = 3L), ncol = 1L)
  expect_equal(trip, .mvr_brute(f$cells, f$zl, picks)[1], tolerance = 1e-10)
  expect_equal(trip, 3 * mv_mean_rank(f$cells) - 1, tolerance = 1e-10)
})

# mobility_value_panel.R (sourced by 50_mobility_value.R) derives the ladder's like-for-like
# kernel rows from the featured label's covariate set, and the featured model can be a
# time-varying FORM of the stage-1 winner ("<stage1>-tvrw1" etc.). The form suffix is not a
# covariate set: left on the label it would be read as one, and every kernel row would be asked
# for at a covariate set ("tvfrailty") that was never scored.
test_that("a time-varying form is parsed from its stage-1 specification", {
  if (!exists("mv_featured_parts", mode = "function") || !exists("TV_FORM_SUFFIXES"))
    testthat::skip("mobility_value_panel.R / 00_config.R not loaded")
  for (s in TV_FORM_SUFFIXES) {
    p <- mv_featured_parts(paste0("Bayes-M14-fill-geo-", s))
    expect_identical(c(p$kernel, p$covset, p$base), c("M14-fill", "geo", "Bayes-M14-fill-geo"))
  }
  p <- mv_featured_parts("Bayes-M13c-fill-med")
  expect_identical(c(p$kernel, p$covset, p$base), c("M13c-fill", "med", "Bayes-M13c-fill-med"))
  expect_error(mv_featured_parts("Gravity-B4"), "not a Bayes")
  # Without the configured suffix list the parse refuses rather than guess one.
  f <- mv_featured_parts; environment(f) <- new.env(parent = baseenv())
  expect_error(f("Bayes-M14-fill-geo-tvfrailty"), "TV_FORM_SUFFIXES")
})
