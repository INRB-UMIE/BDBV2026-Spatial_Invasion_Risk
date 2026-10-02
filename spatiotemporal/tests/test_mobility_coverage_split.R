# =============================================================================
# tests/test_mobility_coverage_split.R — relocation-OD coverage fill and cohort origin split
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Covers (03_mobility_matrices.R):
#   * cover_relocation_od(): uncovered origins take the base row; covered origins keep their
#     measured shape on covered destinations and take the base on uncovered ones.
#   * split_cohort_rows(): per-origin rows that reproduce the pooled cohort profile exactly under
#     population weights, keep each origin's own geography, keep the pooled share where the base
#     carries no information, and fall back to the pooled profile loudly when a split is impossible.
# Run:  Rscript spatiotemporal/tests/run_tests.R
#   or: testthat::test_file("spatiotemporal/tests/test_mobility_coverage_split.R")
# =============================================================================

suppressPackageStartupMessages(library(testthat))

.st <- if (requireNamespace("here", quietly = TRUE)) file.path(here::here(), "spatiotemporal") else "."
if (!exists("cover_relocation_od", mode = "function") || !exists("split_cohort_rows", mode = "function")) {
  if (!exists("OUT_DIR")) try(suppressMessages(source(file.path(.st, "00_config.R"))), silent = TRUE)
  try(suppressWarnings(suppressMessages(source(file.path(.st, "03_mobility_matrices.R")))), silent = TRUE)
}
.need <- function(fn) if (!exists(fn, mode = "function")) skip(paste(fn, "not available"))

# ---- cover_relocation_od -----------------------------------------------------------------------
.cover_fixture <- function() {
  zones <- c("A", "B", "C", "D", "E")
  raw <- matrix(0, 5, 5, dimnames = list(zones, zones))
  raw["A", "B"] <- 10; raw["A", "C"] <- 30; raw["C", "A"] <- 5   # B, D, E: no observed outflow
  M3 <- make_row_stochastic(raw)
  attr(M3, "raw") <- raw
  attr(M3, "od_origins") <- c("A", "B", "C")
  attr(M3, "od_dests")   <- c("A", "B", "C")                     # D and E are not in the table
  base <- matrix(0, 5, 5, dimnames = list(zones, zones))
  base["A", c("B", "C", "D", "E")] <- c(0.3, 0.3, 0.2, 0.2)
  base["B", c("A", "C", "D")]      <- c(0.2, 0.3, 0.5)
  base["C", c("A", "D")]           <- c(0.5, 0.5)
  base["D", "A"] <- 1
  base["E", "A"] <- 1
  list(zones = zones, M3 = M3, base = base)
}

test_that("cover_relocation_od fills exactly the coverage gaps", {
  .need("cover_relocation_od")
  f <- .cover_fixture()
  W <- suppressMessages(cover_relocation_od(f$M3, f$base, f$zones, "T"))
  # uncovered origins take the base row
  for (o in c("B", "D", "E")) expect_equal(unname(W[o, ]), unname(f$base[o, ]))
  # A: measured shape (B 0.25, C 0.75) scaled by 1 - q, q = base mass on {D, E} = 0.4
  expect_equal(unname(W["A", c("B", "C", "D", "E")]), c(0.6 * 0.25, 0.6 * 0.75, 0.2, 0.2))
  # C: measured shape (A 1) scaled by 1 - q, q = base mass on {D, E} = 0.5
  expect_equal(unname(W["C", c("A", "D", "E")]), c(0.5, 0.5, 0))
  # destinations absent from the table now receive weight; rows are stochastic
  expect_true(sum(W[, "D"]) > 0)
  expect_equal(unname(rowSums(W)), rep(1, 5))
  expect_true(all(diag(W) == 0))
})

test_that("cover_relocation_od refuses an M3 without its coverage attributes", {
  .need("cover_relocation_od")
  f <- .cover_fixture()
  bare <- f$M3; attr(bare, "od_dests") <- NULL
  expect_error(cover_relocation_od(bare, f$base, f$zones, "T"), "od_dests")
})

# ---- split_cohort_rows -------------------------------------------------------------------------
.split_fixture <- function() {
  zones <- c("O1", "O2", "D1", "D2", "D3", "X")
  A <- c("D1", "D2", "D3")
  Mc <- matrix(0, 6, 6, dimnames = list(zones, zones))
  Mc["O1", A] <- c(0.5, 0.3, 0.2)
  Mc["O2", A] <- c(0.5, 0.3, 0.2)                # pooled: identical rows
  attr(Mc, "measured")       <- list(O1 = A, O2 = A)
  attr(Mc, "source_origins") <- list(O1 = c("O1", "O2"), O2 = c("O1", "O2"))
  attr(Mc, "cohort_origins") <- c("O1", "O2")
  base <- matrix(0, 6, 6, dimnames = list(zones, zones))
  base["O1", c("O2", "D1", "D2", "D3")] <- c(0.2, 0.5, 0.2, 0.1)   # O1 leans to D1
  base["O2", c("O1", "D1", "D2", "D3")] <- c(0.2, 0.1, 0.2, 0.5)   # O2 leans to D3
  for (z in c("D1", "D2", "D3", "X")) base[z, "O1"] <- 1
  pop <- c(O1 = 3, O2 = 1, D1 = 1, D2 = 1, D3 = 1, X = 1)
  list(zones = zones, A = A, Mc = Mc, base = base, pop = pop,
       src = list(c = c("O1", "O2")), p = c(D1 = 0.5, D2 = 0.3, D3 = 0.2), w = c(0.75, 0.25))
}

test_that("split_cohort_rows reproduces the pooled profile and keeps each origin's geography", {
  .need("split_cohort_rows")
  f <- .split_fixture()
  S <- suppressMessages(split_cohort_rows(f$Mc, f$base, f$src, f$pop, f$zones, character(0)))
  R <- S[c("O1", "O2"), f$A]
  expect_equal(unname(rowSums(R)), c(1, 1), tolerance = 1e-12)
  expect_equal(unname(colSums(f$w * R)), unname(f$p), tolerance = 1e-9)   # exact reproduction
  expect_gt(R["O1", "D1"], R["O2", "D1"])                                 # O1's base favours D1
  expect_gt(R["O2", "D3"], R["O1", "D3"])                                 # O2's base favours D3
  expect_true(all(S[c("O1", "O2"), c("O1", "O2", "X")] == 0))             # only measured cells
  expect_identical(attr(S, "measured"), attr(f$Mc, "measured"))
})

test_that("a destination the base cannot reach keeps the pooled share in every row", {
  .need("split_cohort_rows")
  f <- .split_fixture()
  f$base[c("O1", "O2"), "D3"] <- 0
  S <- suppressMessages(split_cohort_rows(f$Mc, f$base, f$src, f$pop, f$zones, character(0)))
  R <- S[c("O1", "O2"), f$A]
  expect_equal(unname(R[, "D3"]), c(0.2, 0.2))
  expect_equal(unname(rowSums(R)), c(1, 1), tolerance = 1e-12)
  expect_equal(unname(colSums(f$w * R)), unname(f$p), tolerance = 1e-9)
})

test_that("single-origin cohorts and impossible splits keep the pooled profile", {
  .need("split_cohort_rows")
  f <- .split_fixture()
  one <- f$Mc; attr(one, "cohort_origins") <- "O1"
  S1 <- suppressMessages(split_cohort_rows(one, f$base, list(c = "O1"), f$pop, f$zones, character(0)))
  expect_equal(unname(S1["O1", f$A]), unname(f$p))
  bad <- f$pop; bad["O2"] <- 0
  expect_warning(S2 <- suppressMessages(split_cohort_rows(f$Mc, f$base, f$src, bad, f$zones,
                                                          character(0))), "POOLED")
  expect_equal(unname(S2["O1", f$A]), unname(f$p))
  expect_equal(unname(S2["O2", f$A]), unname(f$p))
  bare <- f$Mc; attr(bare, "measured") <- NULL
  expect_error(split_cohort_rows(bare, f$base, f$src, f$pop, f$zones, character(0)), "measured")
})

test_that("the split source composes with the fill: within-cohort links come from the base", {
  .need("split_cohort_rows"); .need("compose_epicentre")
  f <- .split_fixture()
  S <- suppressMessages(split_cohort_rows(f$Mc, f$base, f$src, f$pop, f$zones, character(0)))
  W <- suppressMessages(compose_epicentre(f$base, S, c("O1", "O2"), f$zones, "T",
                                          fill = "unmeasured"))
  expect_equal(unname(W["O1", "O2"]), unname(f$base["O1", "O2"]))      # 0.2, not 0
  expect_equal(unname(W["O2", "O1"]), unname(f$base["O2", "O1"]))
  expect_equal(unname(rowSums(W[c("O1", "O2"), ])), c(1, 1))
  expect_false(isTRUE(all.equal(W["O1", f$A], W["O2", f$A])))           # per-origin rows
})
