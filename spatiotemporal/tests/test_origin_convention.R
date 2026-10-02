# =============================================================================
# tests/test_origin_convention.R — one forecast-origin convention, everywhere
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# A forecast origin is the LAST DAY OF THE LAST TRAINING WEEK: `cut + 6`, never
# `cut + 7`. The weekly grid is anchored (00_config.R WEEK_ANCHOR) so the final
# week ENDS on ANALYSIS_DATE and the deployed as-of date IS that day, so any
# retrospective origin has to be dated the same way to reproduce deployment.
#
# This has drifted twice. The LFO used `cut + 7` until 2026-09-19; the cascade
# calibration still used it afterwards. The consequence is not cosmetic: the
# nowcast multiplier on the most recent week is 2.631x at cut+6 and 2.185x at
# cut+7 — a 20.4% weaker import force. Because `delta` is FITTED on those origins
# and APPLIED to a deployed layer built at 2.631x, the mismatch cannot be absorbed
# by delta; it lands in the published projection.
#
# These tests pin the arithmetic AND the source text, because the arithmetic is
# only reachable through a full fit while the convention is a one-character edit.
# =============================================================================

test_that("the nowcast multiplier at cut+6 is the deployed one, and cut+7 is materially weaker", {
  skip_if_missing <- function(fn) if (!exists(fn, mode = "function")) testthat::skip(paste0(fn, "() unavailable"))
  skip_if_missing("compute_truncation_weights")
  skip_if_missing("effective_onset_sample_delay")
  dly <- tryCatch(effective_onset_sample_delay(), error = function(e) NULL)
  if (is.null(dly) || !is.finite(dly$mean)) testthat::skip("delay params unavailable")

  cut   <- as.Date("2026-08-10")
  weeks <- seq(cut - 7 * 4, cut, by = "week")
  w6 <- suppressMessages(compute_truncation_weights(weeks, analysis_date = cut + 6, delay = dly))
  w7 <- suppressMessages(compute_truncation_weights(weeks, analysis_date = cut + 7, delay = dly))

  last6 <- w6[length(w6)]; last7 <- w7[length(w7)]
  expect_true(is.finite(last6) && last6 > 0)
  # An extra day of reporting makes the most recent week look MORE complete, so its
  # weight rises and the correction multiplier falls.
  expect_gt(last7, last6)
  # The gap is large enough to matter, which is the whole point of the guard.
  expect_gt((1 / last6) / (1 / last7) - 1, 0.10)
})

test_that("no module passes cut + 7 as an as-of / issue / analysis date", {
  st <- file.path(here::here(), "spatiotemporal")
  files <- list.files(st, pattern = "\\.R$", full.names = TRUE)
  files <- files[!grepl("/tests?/", files)]
  bad <- list()
  # Only the ARGUMENTS that define an origin. `cutoff + 7` is legitimate elsewhere — it is
  # the FIRST PREDICTED DAY in a window caption, for instance (.window_caption in
  # 20_forecast_detail.R prints "cutoff + 6 (last trained)" then "cutoff + 7 (first
  # predicted)"), and week_start + 7 is simply the next week. What must never be cut + 7 is
  # the date a reconstruction, a nowcast or an R(t) fit CENSORS ON.
  pat <- paste0("(reaggregate_asof\\s*\\([^)]*|",
                "analysis_date\\s*=\\s*|issue_date\\s*=\\s*|asof\\s*=\\s*)",
                "\\b(cut|cutoff|cutoff_date|origin|origin_date|fold_cut)\\s*\\+\\s*7\\b")
  for (f in files) {
    ln <- readLines(f, warn = FALSE)
    code <- sub("#.*$", "", ln)          # ignore commented-out history
    hit <- grep(pat, code)
    if (length(hit)) bad[[basename(f)]] <- paste0(hit, ": ", trimws(ln[hit]))
  }
  expect_equal(length(bad), 0L,
               info = paste0("an origin must censor at cut + 6:\n",
                             paste(names(bad), unlist(bad), sep = " | ", collapse = "\n")))
})

test_that("the modules that date an origin all use cut + 6", {
  st <- file.path(here::here(), "spatiotemporal")
  # Each of these dates a retrospective origin and must agree with deployment.
  want <- c("16_invasion_eval.R", "33b_cascade_calibration.R")
  for (f in want) {
    p <- file.path(st, f)
    if (!file.exists(p)) next
    code <- sub("#.*$", "", readLines(p, warn = FALSE))
    expect_true(any(grepl("\\b(cut|cutoff_date)\\s*\\+\\s*6\\b", code)),
                info = paste(f, "does not date its origin at cut + 6"))
  }
})

# =============================================================================
# Executable scripts must not run when they are source()d
# =============================================================================
# Several files here are SCRIPTS, not function libraries: sourcing them runs a full
# analysis and overwrites published outputs. run_all.R launches each as its own
# subprocess. A tooling sweep that loads "every module" — or a test reaching for one
# helper — must not trigger the run, and source() gives no warning that it has.
# This happened: a module-load check rewrote six files in key_outputs/ silently.
test_that("scripts that write outputs at top level are guarded against source()", {
  st <- file.path(here::here(), "spatiotemporal")
  # Files whose top level writes published artifacts.
  scripts <- c("43_spread_kinematics.R", "41_kindiv_sweep.R", "40_cascade_next_dominoes.R")
  for (f in scripts) {
    p <- file.path(st, f)
    if (!file.exists(p)) next
    src <- paste(readLines(p, warn = FALSE), collapse = "\n")
    expect_true(grepl("is_script_run\\(", src, fixed = FALSE),
                info = paste(f, "has top-level writes but no is_script_run() guard"))
  }
})

test_that("is_script_run() is FALSE unless R was launched on that exact file", {
  if (!exists("is_script_run", mode = "function")) testthat::skip("is_script_run() unavailable")
  # Under testthat there is no --file= for these names, so every call must be FALSE.
  expect_false(is_script_run("43_spread_kinematics.R"))
  expect_false(is_script_run("41_kindiv_sweep.R"))
  expect_false(is_script_run("40_cascade_next_dominoes.R"))
  # And it must not be fooled by a path: only the basename is compared, but a DIFFERENT
  # basename must never match.
  expect_false(is_script_run("definitely_not_the_running_file.R"))
})
