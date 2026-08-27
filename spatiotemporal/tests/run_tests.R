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

# ---------------------------------------------------------------------------
# Source the LIVE analysis modules that run_all.R sources (in the same order).
# Each source() is wrapped in tryCatch so that a module that fails to load
# (e.g. an optional dependency such as deSolve absent on a CI runner) does not
# prevent the remaining modules — and their unit tests — from running.
#
# NB: scripts 09–14 were retired to spatiotemporal/stale/ and are NEVER sourced.
# The metric functions that used to live in 11_metrics.R / 12_calibration.R are
# now provided by the live modules (16_invasion_eval.R: .auc_pr/.log_score/.brier;
# 04b_epinowcast.R: the WIS `wis_one`; 19_spacetime_eval.R / 20_forecast_detail.R:
# spatiotemporal skill, detection curve, balanced accuracy). The test suite has
# been ported to those live equivalents.
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

# Live module list — mirrors run_all.R's source order (09–14 deliberately absent).
live_modules <- c(
  "01_data_prep.R", "02_epi_params.R", "03_mobility_matrices.R",
  "04_nowcasting.R", "04b_epinowcast.R", "04c_dhis2_delay_windows.R",
  "05_baseline_models.R",
  "06_simple_models.R", "07_hhh4_model.R", "08_stochastic_seir.R",
  "15_workhorse.R", "16_invasion_eval.R", "18_ensemble.R",
  "19_spacetime_eval.R", "20_forecast_detail.R", "21_bayesian_renewal.R",
  "17_invasion_viz.R"
)
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
