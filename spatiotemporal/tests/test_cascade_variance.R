# =============================================================================
# test_cascade_variance.R - cascade_reach_table()'s variance decomposition
# =============================================================================
# The published per-zone reach credible interval is formed by splitting total
# between-iteration variance into a PARAMETER part (between draw groups) and a
# PROCESS part (within a group), and keeping only the former:
#     s2_param = var(p_draw) - mean(p_draw*(1-p_draw))/(n_rep-1)
#
# That makes the LEVEL AT WHICH A NOISE SOURCE IS DRAWN decide whether it reaches
# the interval at all. CASCADE_R_RW_SIGMA (the mean-reverting walk on log R) exists
# specifically to widen the interval with horizon, but it was drawn PER ITERATION,
# which put it entirely inside the process term and subtracted it back out — at large
# sigma it made the published interval NARROWER. It is now drawn per parameter group.
#
# These tests pin the property, not the implementation: a per-group shock must widen
# the interval and a per-iteration shock must not. Without them the placement could
# silently revert and the only symptom would be a too-tight published CrI.
# =============================================================================

skip_if_missing <- function(fn_name) {
  if (!existsFunction(fn_name))
    skip(paste("Function", fn_name, "not available - source the relevant analysis script first"))
}

# Build a synthetic sim whose reach indicator carries a shock applied either once per
# parameter group or independently per iteration.
.mk_sim <- function(level, sigma, nz = 300L, D = 150L, R = 8L, seed = 7L) {
  set.seed(seed)
  M  <- D * R
  dg <- rep(seq_len(D), each = R)
  base  <- runif(nz, 0.05, 0.5)
  shock <- if (identical(level, "group")) rep(stats::rnorm(D, 0, sigma), each = R)
           else stats::rnorm(M, 0, sigma)
  pr  <- pmin(pmax(outer(base, rep(1, M)) * exp(matrix(shock, nz, M, byrow = TRUE)), 0), 1)
  tau <- ifelse(matrix(runif(nz * M), nz, M) < pr, 1L, NA_integer_)
  list(zones_all = sprintf("Z%03d", seq_len(nz)), n_mc = M, n_rep = R, draw_group = dg,
       tau = tau, est_tau = tau, affected0 = rep(FALSE, nz),
       new_by_week = matrix(0, 13L, M))
}

.mean_width <- function(sim) {
  t <- cascade_reach_table(sim, horizons = 1L)
  mean(t$p_hi - t$p_lo, na.rm = TRUE)
}

test_that("cascade_reach_table: a PER-GROUP shock widens the credible interval", {
  skip_if_missing("cascade_reach_table")
  w0  <- .mean_width(.mk_sim("group", 0))
  w30 <- .mean_width(.mk_sim("group", 0.30))
  w60 <- .mean_width(.mk_sim("group", 0.60))
  expect_gt(w30, w0)
  expect_gt(w60, w30)
  # It is the dominant term by sigma = 0.6, not a marginal effect.
  expect_gt(w60, 3 * w0)
})

test_that("cascade_reach_table: a PER-ITERATION shock is absorbed by the process term", {
  skip_if_missing("cascade_reach_table")
  # At sigma = 0 the two placements are the SAME experiment, so they must agree exactly.
  expect_equal(.mean_width(.mk_sim("iteration", 0)),
               .mean_width(.mk_sim("group", 0)), tolerance = 1e-12)
  # The property that matters: the same shock reaches the interval overwhelmingly more when
  # drawn per parameter group than per iteration. Asserting "per-iteration never widens at
  # all" would be wrong — `within` models the process term as p(1-p)/(n_rep-1), i.e. assumes
  # pure Bernoulli within-group variation, so a per-iteration shock partly leaks through.
  # Measured: +21% (per-iteration) against +938% (per-group) at sigma 0.6.
  for (sg in c(0.30, 0.60)) {
    wi <- .mean_width(.mk_sim("iteration", sg))
    wg <- .mean_width(.mk_sim("group", sg))
    expect_gt(wg, 3 * wi)
  }
  # And per-iteration stays near its sigma = 0 width rather than tracking the shock.
  expect_lt(.mean_width(.mk_sim("iteration", 0.60)),
            1.6 * .mean_width(.mk_sim("iteration", 0)))
})

test_that("simulate_cascade draws the R walk per parameter group, not per iteration", {
  skip_if_missing("simulate_cascade")
  src <- deparse(args(simulate_cascade))
  body_txt <- paste(deparse(body(simulate_cascade)), collapse = "\n")
  # The group-level array must exist and be indexed by the group id `g`.
  expect_true(grepl("Z_Rw_grp", body_txt, fixed = TRUE))
  expect_true(grepl("Z_Rw_grp[, , g]", body_txt, fixed = TRUE))
  # And the walk must NOT be refilled from a fresh rnorm inside the per-iteration block.
  expect_false(grepl("Z_Rw[] <- stats::rnorm", body_txt, fixed = TRUE))
})
