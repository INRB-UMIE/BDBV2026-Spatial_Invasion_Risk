# =============================================================================
# test_ensemble.R - 18_ensemble.R (ensemble_forecasts / append_ensembles)
# =============================================================================
# The module had NO test coverage, yet it has four pipeline call sites and produces
# bayes_ensemble_risk_scores_all_zones.csv, which is in the CI publish allowlist.
# =============================================================================

# Local helper, matching the convention in every other file in this suite (each defines its
# own; there is no shared helper file).
skip_if_missing <- function(fn_name) {
  if (!existsFunction(fn_name))
    skip(paste("Function", fn_name, "not available - source the relevant analysis script first"))
}

test_that("ensemble_forecasts: mean pooling, and every column uses ONE member set", {
  skip_if_missing("ensemble_forecasts")
  fc <- data.frame(
    health_zone = rep("A", 3), horizon = rep(1L, 3),
    method      = c("m1", "m2", "m3"),
    p_invasion  = c(0.1, 0.2, NA_real_),      # m3 contributes no usable probability
    mu_forecast = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE)
  out <- ensemble_forecasts(fc, members = c("m1", "m2", "m3"),
                            combine = "mean", label = "E", min_members = 2L)
  expect_equal(nrow(out), 1L)
  expect_equal(out$p_invasion[1], 0.15, tolerance = 1e-12)
  # The defect: mu_forecast used to be the 3-member mean (0.2) while p_invasion was the
  # 2-member mean (0.15), so one row described two different ensembles.
  expect_equal(out$mu_forecast[1], 0.15, tolerance = 1e-12,
               info = "every value column must aggregate the SAME member set")
})

test_that("ensemble_forecasts: min_members is a true DISTINCT-member guard", {
  skip_if_missing("ensemble_forecasts")
  dup <- data.frame(health_zone = rep("A", 2), horizon = rep(1L, 2),
                    method = c("m1", "m1"), p_invasion = c(0.1, 0.3),
                    stringsAsFactors = FALSE)
  # The "only 1 of 2 members present" warning is the documented behaviour; pin it rather than
  # letting it leak into the suite's warning count.
  expect_warning(
    out <- ensemble_forecasts(dup, members = c("m1", "m2"), combine = "mean",
                              label = "E", min_members = 2L),
    "only 1 of 2 members present")
  expect_null(out, info = "one member emitting two rows is not a 2-member ensemble")
})

test_that("ensemble_forecasts: median pooling and quantile monotonicity", {
  skip_if_missing("ensemble_forecasts")
  fc <- data.frame(health_zone = rep("A", 3), horizon = rep(1L, 3),
                   method = c("m1", "m2", "m3"),
                   p_invasion = c(0.10, 0.20, 0.60),
                   q05 = c(0.01, 0.02, 0.03), q95 = c(0.5, 0.6, 0.9),
                   stringsAsFactors = FALSE)
  out <- ensemble_forecasts(fc, members = c("m1","m2","m3"), combine = "median",
                             label = "E", min_members = 2L)
  expect_equal(out$p_invasion[1], 0.20, tolerance = 1e-12)   # median, not mean (0.30)
  expect_lte(out$q05[1], out$p_invasion[1])
  expect_gte(out$q95[1], out$p_invasion[1])
})

test_that("ensemble_forecasts: carry columns take the first NON-NA value", {
  skip_if_missing("ensemble_forecasts")
  fc <- data.frame(health_zone = rep("A", 3), horizon = rep(1L, 3),
                   method = c("m1","m2","m3"), p_invasion = c(0.1, 0.2, 0.3),
                   is_new_invasion = c(NA_integer_, 1L, 1L),
                   stringsAsFactors = FALSE)
  out <- ensemble_forecasts(fc, members = c("m1","m2","m3"), combine = "mean",
                             label = "E", min_members = 2L)
  expect_equal(out$is_new_invasion[1], 1L,
               info = "a leading NA outcome must not make the ensemble row unscorable")
})

test_that("append_ensembles is a no-op when fewer than two members are present", {
  skip_if_missing("append_ensembles")
  fc <- data.frame(health_zone = "A", horizon = 1L, method = "m1", p_invasion = 0.1,
                   stringsAsFactors = FALSE)
  expect_equal(nrow(suppressWarnings(suppressMessages(append_ensembles(fc, character(0))))), nrow(fc))
  expect_equal(nrow(suppressWarnings(suppressMessages(append_ensembles(fc, "m1")))), nrow(fc))
})

test_that("ensemble: each column is pooled independently (the Jensen gap is intentional)", {
  skip_if_missing("ensemble_forecasts")
  # Two members whose mu_forecast and p_invasion are BOTH genuine posterior means of their
  # own functional, and therefore do not satisfy p = 1 - exp(-mu) individually. Linear
  # pooling must preserve each column's own mean; deriving one from the other would replace
  # a correct posterior mean with a Jensen-biased transform. This is pinned because the
  # "inconsistency" looks like a bug and has been proposed as one.
  fc <- tibble::tibble(
    method      = c("A", "B"),
    health_zone = c("Z", "Z"),
    horizon     = c(1L, 1L),
    mu_forecast = c(0.50, 20.00),   # right-skewed hazard posterior across members
    p_invasion  = c(0.30, 0.40)     # NOT 1 - exp(-mu): both are posterior means
  )
  out <- ensemble_forecasts(fc, members = c("A", "B"), combine = "mean",
                            label = "Ensemble-mean", min_members = 2L)
  expect_equal(nrow(out), 1L)
  # Each column is the mean of its own member values.
  expect_equal(out$mu_forecast[1], mean(c(0.50, 20.00)))
  expect_equal(out$p_invasion[1],  mean(c(0.30, 0.40)))
  # And the row deliberately does NOT satisfy p = 1 - exp(-mu).
  expect_gt(abs((1 - exp(-out$mu_forecast[1])) - out$p_invasion[1]), 0.5)
})
