# =============================================================================
# tests/run_tests.R — Master test runner
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Usage (from repo root):
#   Rscript spatiotemporal/tests/run_tests.R
# Or interactively:
#   source("spatiotemporal/tests/run_tests.R")
#
# All individual test files live in the same directory.  Only functions
# that exist in the environment at test-run time are exercised; test files
# wrap missing-function calls in tryCatch / skip() so that partially-built
# analysis scripts do not cause the whole suite to fail.
# =============================================================================

suppressPackageStartupMessages({
  library(testthat)
  library(here)
})

# ---------------------------------------------------------------------------
# Load configuration (always required)
# ---------------------------------------------------------------------------
source(file.path(here::here(), "spatiotemporal", "00_config.R"))
# Tests that load the real line list must not write its earliest-sample cache into the
# published output tree (00_config.R, .earliest_sample_dates()).
EARLIEST_SAMPLE_CACHE_DIR <- file.path(tempdir(), "earliest_sample_cache")
dir.create(EARLIEST_SAMPLE_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Source the LIVE analysis modules that run_all.R sources (in the same order).
# Each source() is wrapped in tryCatch so that a module that fails to load
# (e.g. an optional dependency such as deSolve absent on a CI runner) does not
# prevent the remaining modules — and their unit tests — from running.
#
# NB: scripts 09-14 were DELETED, not archived (this note used to point at an archive
# directory that does not exist). The metric functions that lived in
# 11_metrics.R / 12_calibration.R are now provided by the live modules — 16_invasion_eval.R
# (.auc_pr / .log_score / .brier), 19_spacetime_eval.R and 20_forecast_detail.R
# (spatiotemporal skill, detection curve, balanced accuracy) — and the suite tests those.
# ---------------------------------------------------------------------------

safe_source <- function(path, label) {
  tryCatch(
    source(path),
    error = function(e) {
      message(sprintf(
        "[run_tests] WARNING: could not source %s (%s). ",
        label, conditionMessage(e)
      ), appendLF = FALSE)
      message("Tests that require functions from this script will be skipped.")
    }
  )
}

# THE MODULE LIST COMES FROM 00_config.R (load_pipeline_modules / PIPELINE_MODULES), so the
# suite cannot test a different module set from the one the pipeline loads. This file used to
# hand-maintain its own list and it had drifted: 22_daily_reissue.R was missing, which
# silently SKIPPED every test that needs linelist_observation_date() — the one function whose
# absence also makes .trunc_fit_asof_cohort() return NULL without a word.
#
# `extra` covers what the TESTS need beyond the pipeline layers:
#   04c  the interval-censored delay fitters (test_dhis2_delay.R, test_delay_strata.R)
#   17   the viz layer (its scoring-adjacent helpers are asserted on)
live_modules <- c(unlist(PIPELINE_MODULES[c("data", "model", "cascade")], use.names = FALSE),
                  "04c_dhis2_delay_windows.R", "17_invasion_viz.R")
for (m in live_modules) {
  safe_source(file.path(here::here(), "spatiotemporal", m), m)
}

# ---------------------------------------------------------------------------
# Run all test_*.R files in this directory
# ---------------------------------------------------------------------------
message("\n=== Running spatiotemporal test suite ===\n")

test_dir(
  path     = file.path(here::here(), "spatiotemporal", "tests"),
  reporter = "progress"
)

message("\n=== Test suite complete ===\n")
