# =============================================================================
# tests/test_mobility.R — Unit tests for 03_mobility_matrices.R functions
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Tests cover:
#   - make_row_stochastic()    row-stochastic normalisation with diagonal zeroing
#   - assert_mobility_matrix() invariant checker (must pass/fail correctly)
#   - build_M5()               radiation model structural properties
#   - build_M6()               travel-time decay matrix
#   - build_M8()               composite epicentre/gravity matrix
#   - .embed_od_to_519()       name-indexed OD embedding (shared by build_M3/build_M15)
#   - build_M15()              symmetrised relocation-OD kernel O + t(O) (NOT inflow-informed:
#                              no independent inflow table exists)
#   - build_M17()              all-kernel consensus ensemble (convex mean + overlaid source rows)
#
# Most tests use small in-memory matrices so no data files are required; the one
# build_M15 real-data test skips when the Flowminder OD file is absent.
# =============================================================================

existsFunction <- function(fn_name) {
  exists(fn_name, mode = "function", envir = .GlobalEnv, inherits = TRUE)
}

skip_if_missing <- function(fn_name) {
  if (!existsFunction(fn_name)) {
    skip(paste("Function", fn_name, "not available — source 03_mobility_matrices.R first"))
  }
}

# =============================================================================
# 1. make_row_stochastic()
# =============================================================================

test_that("make_row_stochastic: rows sum to 1 (or 0) and diagonal is zeroed", {
  skip_if_missing("make_row_stochastic")

  M_raw <- matrix(
    c(0, 3, 1,
      2, 0, 4,
      1, 1, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )
  M_rs <- make_row_stochastic(M_raw)

  row_sums <- rowSums(M_rs)
  expect_true(
    all(abs(row_sums - 1) < 1e-9 | row_sums == 0),
    info = "All rows must sum to 1 (or 0 if all-zero outflow)"
  )
  expect_true(all(M_rs >= 0),          info = "No negative weights after normalisation")
  expect_true(all(diag(M_rs) == 0),    info = "Diagonal must be zeroed before normalisation")
})

test_that("make_row_stochastic: all-zero row remains all-zero after normalisation", {
  skip_if_missing("make_row_stochastic")

  M <- matrix(
    c(0, 5, 3,
      0, 0, 0,   # zero-outflow row
      2, 1, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )
  M_rs <- make_row_stochastic(M)

  expect_equal(unname(rowSums(M_rs["B", , drop = FALSE])), 0,
               info = "Zero-outflow row must remain all-zero (not NaN)")
  expect_equal(unname(rowSums(M_rs["A", , drop = FALSE])), 1, tolerance = 1e-9,
               info = "Active row must sum to 1")
})

test_that("make_row_stochastic: non-zero diagonal is zeroed before normalisation", {
  skip_if_missing("make_row_stochastic")

  # All weight on the diagonal — after zeroing, rows should all be zero
  M_diag <- diag(c(10, 20, 30))
  dimnames(M_diag) <- list(c("A", "B", "C"), c("A", "B", "C"))
  M_rs <- make_row_stochastic(M_diag)

  expect_equal(rowSums(M_rs), c(A = 0, B = 0, C = 0),
               info = "Diagonal-only matrix → all rows zero after zeroing diagonal")
})

test_that("make_row_stochastic: relative weights are preserved", {
  skip_if_missing("make_row_stochastic")

  M <- matrix(
    c(0, 1, 3,
      2, 0, 2,
      6, 3, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )
  M_rs <- make_row_stochastic(M)

  # Row A: off-diagonal = [1, 3] → proportions [0.25, 0.75]
  expect_equal(M_rs["A", "B"], 0.25, tolerance = 1e-9,
               info = "Row A, col B: expected proportion 0.25")
  expect_equal(M_rs["A", "C"], 0.75, tolerance = 1e-9,
               info = "Row A, col C: expected proportion 0.75")
})

# =============================================================================
# 2. assert_mobility_matrix()
# =============================================================================

test_that("assert_mobility_matrix: passes for valid row-stochastic matrix", {
  skip_if_missing("assert_mobility_matrix")

  zones <- c("A", "B", "C")
  W <- matrix(
    c(0,   0.6, 0.4,
      0.3, 0,   0.7,
      0.5, 0.5, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  expect_silent(assert_mobility_matrix(W, zones, "W_valid"))
})

test_that("assert_mobility_matrix: fails for non-zero diagonal", {
  skip_if_missing("assert_mobility_matrix")

  zones <- c("A", "B", "C")
  W_bad <- matrix(
    c(0.1, 0.5, 0.4,
      0.3, 0.0, 0.7,
      0.5, 0.5, 0.0),
    nrow = 3, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  expect_error(
    assert_mobility_matrix(W_bad, zones, "W_bad_diag"),
    info = "Non-zero diagonal must trigger error"
  )
})

test_that("assert_mobility_matrix: fails for row sums > 1", {
  skip_if_missing("assert_mobility_matrix")

  zones <- c("A", "B", "C")
  W_bad <- matrix(
    c(0.0, 0.8, 0.8,   # row sum = 1.6 > 1
      0.3, 0.0, 0.7,
      0.5, 0.5, 0.0),
    nrow = 3, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  expect_error(
    assert_mobility_matrix(W_bad, zones, "W_bad_rs"),
    info = "Row sums > 1 must trigger error"
  )
})

test_that("assert_mobility_matrix: fails for wrong dimension", {
  skip_if_missing("assert_mobility_matrix")

  zones <- c("A", "B", "C", "D")
  W_small <- matrix(
    c(0, 0.6, 0.4,
      0.3, 0, 0.7,
      0.5, 0.5, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A","B","C"), c("A","B","C"))
  )
  expect_error(
    assert_mobility_matrix(W_small, zones, "W_wrong_dim"),
    info = "Wrong dimension must trigger error"
  )
})

# =============================================================================
# 3. build_M6() — OSRM travel-time decay matrix
# =============================================================================

test_that("build_M6 exponential: diagonal is zero and rows sum to 1", {
  skip_if_missing("build_M6")

  zones <- c("A", "B", "C")
  osrm  <- matrix(
    c(0,  60, 120,
      60,  0,  90,
      120, 90,  0),
    nrow = 3, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  M6a <- build_M6(osrm, zones, variant = "exp", kappa = 120)

  expect_true(all(diag(M6a) == 0),  info = "M6a diagonal must be zero")
  expect_true(all(dim(M6a) == c(3, 3)), info = "M6a must be 3×3")
  rs <- rowSums(M6a)
  expect_true(all(rs >= 0 & rs <= 1 + 1e-9),
              info = "M6a row sums must be in [0, 1]")
  # Active rows must sum to 1
  active_rows <- rs > 1e-9
  expect_true(all(abs(rs[active_rows] - 1) < 1e-9),
              info = "M6a active rows must sum exactly to 1")
})

test_that("build_M6 exponential: closer zone receives higher weight", {
  skip_if_missing("build_M6")

  zones <- c("Origin", "Near", "Far")
  osrm  <- matrix(
    c(0,   30, 300,
      30,   0, 270,
      300, 270,  0),
    nrow = 3, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  M6a <- build_M6(osrm, zones, variant = "exp", kappa = 120)

  # From "Origin": Near (30 min) must get higher weight than Far (300 min)
  expect_true(
    M6a["Origin", "Near"] > M6a["Origin", "Far"],
    info = "Exponential decay: closer zone must have higher weight"
  )
})

test_that("build_M6 power-law: diagonal is zero, row sums in [0,1]", {
  skip_if_missing("build_M6")

  zones <- c("A", "B", "C")
  osrm  <- matrix(
    c(0,  45, 180,
      45,  0, 135,
      180,135,  0),
    nrow = 3, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  M6b <- build_M6(osrm, zones, variant = "power", gamma = 1)

  expect_true(all(diag(M6b) == 0),    info = "M6b diagonal must be zero")
  expect_true(all(M6b >= 0),          info = "M6b must be non-negative")
  rs <- rowSums(M6b)
  expect_true(all(rs >= 0 & rs <= 1 + 1e-9), info = "M6b row sums in [0,1]")
})

test_that("build_M6: both variants produce valid zone names in row/colnames", {
  skip_if_missing("build_M6")

  zones <- c("Alpha", "Beta", "Gamma")
  osrm  <- matrix(
    c(0,  60, 90,
      60,  0, 45,
      90, 45,  0),
    nrow = 3, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  for (variant in c("exp", "power")) {
    M6 <- build_M6(osrm, zones, variant = variant, kappa = 60, gamma = 1)
    expect_identical(rownames(M6), zones,
                     info = paste(variant, ": rownames must match zones"))
    expect_identical(colnames(M6), zones,
                     info = paste(variant, ": colnames must match zones"))
  }
})

# =============================================================================
# 4. build_M5() — Radiation model
# =============================================================================

test_that("build_M5: diagonal is zero, row sums in [0, 1]", {
  skip_if_missing("build_M5")

  zones <- c("A", "B", "C", "D")
  pop   <- c(A = 100000, B = 50000, C = 30000, D = 20000)
  osrm  <- matrix(
    c(  0,  30,  60, 120,
       30,   0,  45,  90,
       60,  45,   0,  60,
      120,  90,  60,   0),
    nrow = 4, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  M5 <- build_M5(pop, osrm, zones)

  expect_equal(dim(M5), c(4L, 4L),   info = "M5 must be 4×4")
  expect_true(all(diag(M5) == 0),    info = "M5 diagonal must be zero")
  expect_true(all(M5 >= 0),          info = "M5 must be non-negative")
  rs <- rowSums(M5)
  expect_true(all(rs <= 1 + 1e-9 & rs >= 0),
              info = "M5 row sums must be in [0, 1]")
})

test_that("build_M5: largest population zone sends most outflow to nearest neighbour", {
  skip_if_missing("build_M5")

  # Zone A (largest pop) and Zone B (nearest, medium pop) vs Zone C (far, tiny pop)
  zones <- c("A", "B", "C")
  pop   <- c(A = 500000, B = 200000, C = 5000)
  osrm  <- matrix(
    c( 0, 30, 300,
      30,  0, 270,
      300,270,  0),
    nrow = 3, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  M5 <- build_M5(pop, osrm, zones)

  # From zone A: zone B (near, large pop) should get higher flux than zone C (far, tiny pop)
  expect_true(
    M5["A", "B"] > M5["A", "C"],
    info = "Radiation model: near + large pop zone gets more outflow than far + tiny pop"
  )
})

test_that("build_M5: zone names are preserved in output matrix", {
  skip_if_missing("build_M5")

  zones <- c("Bunia", "Nizi", "Irumu")
  pop   <- c(Bunia = 300000, Nizi = 80000, Irumu = 120000)
  osrm  <- matrix(
    c( 0, 45, 90,
      45,  0, 60,
      90, 60,  0),
    nrow = 3, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  M5 <- build_M5(pop, osrm, zones)

  expect_identical(rownames(M5), zones, info = "M5 rownames must match zones")
  expect_identical(colnames(M5), zones, info = "M5 colnames must match zones")
})

# =============================================================================
# 5. build_M8() — Composite epicentre/gravity matrix
# =============================================================================

test_that("build_M8: epicentre rows come from M1, non-epicentre from M4", {
  skip_if_missing("build_M8")

  zones <- c("Bunia", "Mongbalu", "Rwampara", "Lita", "Nizi", "Bambu")
  n     <- length(zones)

  # Construct two distinguishable row-stochastic matrices
  set.seed(20260704L)
  make_stochastic_mock <- function(zn) {
    M <- matrix(runif(length(zn)^2), length(zn), length(zn),
                dimnames = list(zn, zn))
    diag(M) <- 0
    M / rowSums(M)
  }
  M1_mock <- make_stochastic_mock(zones)
  M4_mock <- make_stochastic_mock(zones)
  # Ensure M1 and M4 differ substantially
  while (max(abs(M1_mock - M4_mock)) < 0.1) {
    M4_mock <- make_stochastic_mock(zones)
  }

  epicentre_zones <- c("Bunia", "Mongbalu", "Rwampara")
  M8 <- build_M8(M1_mock, M4_mock, epicentre_zones, zones)

  # Epicentre rows must match M1 rows exactly (after possible re-normalisation)
  for (z in epicentre_zones) {
    # Re-normalise M1 row to compare against M8 (build_M8 calls make_row_stochastic)
    m1_row_norm <- M1_mock[z, ]
    m1_row_norm[z] <- 0
    m1_row_norm <- m1_row_norm / sum(m1_row_norm)
    expect_equal(
      unname(M8[z, ]), unname(m1_row_norm),
      tolerance = 1e-9,
      info = paste(z, "row in M8 must equal the normalised M1 row")
    )
  }

  # Non-epicentre rows must match M4 rows (after re-normalisation)
  for (z in setdiff(zones, epicentre_zones)) {
    m4_row_norm <- M4_mock[z, ]
    m4_row_norm[z] <- 0
    m4_row_norm <- m4_row_norm / sum(m4_row_norm)
    expect_equal(
      unname(M8[z, ]), unname(m4_row_norm),
      tolerance = 1e-9,
      info = paste(z, "row in M8 must equal the normalised M4 row")
    )
  }
})

test_that("build_M8: output passes assert_mobility_matrix invariants", {
  skip_if_missing("build_M8")
  skip_if_missing("assert_mobility_matrix")

  zones <- c("Bunia", "Mongbalu", "Rwampara", "Lita", "Nizi")
  set.seed(42L)
  make_stochastic_mock <- function(zn) {
    M <- matrix(runif(length(zn)^2) + 0.01, length(zn), length(zn),
                dimnames = list(zn, zn))
    diag(M) <- 0
    M / rowSums(M)
  }
  M1m <- make_stochastic_mock(zones)
  M4m <- make_stochastic_mock(zones)
  M8  <- build_M8(M1m, M4m, c("Bunia", "Mongbalu", "Rwampara"), zones)

  # assert_mobility_matrix should not error (expect_silent takes no `info` arg).
  expect_silent(
    assert_mobility_matrix(M8, zones, "M8_test")
  )
})

test_that("build_M8: all diagonal elements are zero", {
  skip_if_missing("build_M8")

  zones <- c("X", "Y", "Z")
  M_base <- matrix(
    c(0, 0.6, 0.4,
      0.3, 0, 0.7,
      0.5, 0.5, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(zones, zones)
  )
  M8 <- build_M8(M_base, M_base, c("X"), zones)
  expect_true(all(diag(M8) == 0), info = "M8 diagonal must be all-zero")
})

test_that("build_M8: output has correct dimensions and zone names", {
  skip_if_missing("build_M8")

  zones <- c("Alpha", "Beta", "Gamma", "Delta")
  set.seed(99L)
  make_stochastic_mock <- function(zn) {
    M <- matrix(runif(length(zn)^2) + 0.01, length(zn), length(zn),
                dimnames = list(zn, zn))
    diag(M) <- 0
    M / rowSums(M)
  }
  M1m <- make_stochastic_mock(zones)
  M4m <- make_stochastic_mock(zones)
  M8  <- build_M8(M1m, M4m, c("Alpha", "Beta"), zones)

  expect_equal(dim(M8), c(4L, 4L),       info = "M8 must be 4×4")
  expect_identical(rownames(M8), zones,   info = "M8 rownames must match zones")
  expect_identical(colnames(M8), zones,   info = "M8 colnames must match zones")
})

# =============================================================================
# 6. .embed_od_to_519() — OD embedding shared by build_M3 / build_M15
# =============================================================================

test_that(".embed_od_to_519: embeds by NAME (not position), zero diagonal, absent zones zeroed", {
  skip_if_missing(".embed_od_to_519")

  zones <- c("A", "B", "C", "D")
  # 3-zone raw OD whose row order and column order DIFFER (as in the real Flowminder files);
  # assign flows BY NAME so the test checks that embedding is name-indexed, not positional.
  M <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("B", "C", "A")))
  M["A", "B"] <- 5   # A -> B
  M["C", "A"] <- 4   # C -> A
  M["A", "C"] <- 2   # A -> C
  E <- .embed_od_to_519(M, zones, aliases = character(0), label = "test")

  expect_equal(dim(E), c(4L, 4L))
  expect_identical(rownames(E), zones)
  expect_identical(colnames(E), zones)
  expect_true(all(diag(E) == 0), info = "embedded OD must have a zero diagonal")
  # Values placed BY NAME despite the differing raw row/col order.
  expect_equal(E["A", "B"], 5)
  expect_equal(E["C", "A"], 4)
  expect_equal(E["A", "C"], 2)
  # Zone D is absent from the OD -> zero row and zero column
  expect_true(all(E["D", ] == 0) && all(E[, "D"] == 0))
})

# =============================================================================
# 7. M15 symmetrisation property (S = O + t(O)) and the real-data builder
# =============================================================================

test_that("M15 symmetrisation: S = O + t(O) is symmetric and yields a valid kernel", {
  skip_if_missing("make_row_stochastic")
  skip_if_missing("assert_mobility_matrix")

  zones <- c("A", "B", "C", "D")
  O <- matrix(0, 4, 4, dimnames = list(zones, zones))
  O["A", "B"] <- 10; O["B", "C"] <- 4; O["C", "A"] <- 7   # directed, asymmetric flow
  S <- O + t(O); diag(S) <- 0

  expect_true(isSymmetric(S), info = "the symmetrised kernel O + t(O) must be symmetric")
  # symmetrisation fills the reciprocal edge B->A that the directed matrix lacked
  expect_equal(S["A", "B"], S["B", "A"])
  expect_true(S["B", "A"] > 0 && O["B", "A"] == 0)

  W <- make_row_stochastic(S)
  expect_silent(assert_mobility_matrix(W, zones, "M15_test"))
})

test_that("build_M15: real Flowminder data -> valid kernel, distinct from M3, no sparser", {
  skip_if_missing("build_M15")
  skip_if_missing("build_M3")
  skip_if_missing("load_aliases")
  skip_if_missing("load_worldpop")
  f_out <- file.path(FLOWMINDER_DIR, "flowminder__outflow__static.matrix.csv")
  skip_if_not(file.exists(f_out), "Flowminder OD static file not present")

  zones <- names(load_worldpop())
  al    <- load_aliases()
  M15 <- build_M15(zones, al)
  expect_silent(assert_mobility_matrix(M15, zones, "M15"))

  M3 <- build_M3(zones, al)
  expect_false(isTRUE(all.equal(unclass(M15), unclass(M3))),
               info = "symmetrised M15 must differ from the directed M3")
  # M15 fills reciprocal edges, so it has no MORE zero-outflow origin rows than the directed M3
  expect_lte(sum(rowSums(M15) < 1e-9), sum(rowSums(M3) < 1e-9))
})

# =============================================================================
# 8. build_M17() — all-kernel consensus ensemble
# =============================================================================

test_that("build_M17: convex mean is row-stochastic with source rows overlaid", {
  skip_if_missing("build_M17")
  skip_if_missing("compose_epicentre")
  skip_if_missing("assert_mobility_matrix")

  zones <- c("Bunia", "Lita", "Nizi", "Beni", "Katwa")
  set.seed(7L)
  mk <- function() {
    M <- matrix(runif(length(zones)^2) + 0.01, length(zones), length(zones),
                dimnames = list(zones, zones))
    diag(M) <- 0
    M / rowSums(M)
  }
  bases   <- list(a = mk(), b = mk(), c = mk())
  epi     <- mk()
  origins <- c("Bunia", "Beni")
  M17 <- build_M17(bases, epi, origins, zones, "M17_test")

  expect_silent(assert_mobility_matrix(M17, zones, "M17_test"))
  # Source rows overlaid from the empirical kernel (already row-stochastic -> pass through).
  for (z in origins)
    expect_equal(unname(M17[z, ]), unname(epi[z, ]), tolerance = 1e-9,
                 info = paste(z, "source row must equal the overlaid empirical kernel row"))
  # A non-source row is a genuine blend of the members (differs from any single base).
  expect_false(isTRUE(all.equal(unname(M17["Lita", ]), unname(bases$a["Lita", ]))),
               info = "non-source rows must be a blend, not a single member")
  # Guard: fewer than two valid bases must error.
  expect_error(build_M17(list(a = mk()), epi, origins, zones, "M17_bad"))
})

# =============================================================================
# 9. .resolve_source_labels() — raw source labels -> canonical zone names
# =============================================================================

test_that(".resolve_source_labels: exact names, province disambiguation, aliases, case", {
  skip_if_missing(".resolve_source_labels")

  zones <- c("Bunia", "Kilo", "Ngandajika", "Massa",
             "Lubunga (Tshopo)", "Lubunga (Kasaï-Central)",
             "Bili (Bas-Uele)", "Bili (Nord-Ubangi)")

  # exact canonical name wins
  expect_equal(.resolve_source_labels("Kilo", "Ituri", zones), "Kilo")

  # colliding bare name is resolved BY PROVINCE, and punctuation/accents in the
  # province label do not matter ("Kasaï Central" vs the canonical "Kasaï-Central")
  expect_equal(.resolve_source_labels("Lubunga", "Tshopo", zones), "Lubunga (Tshopo)")
  expect_equal(.resolve_source_labels("Lubunga", "Kasaï Central", zones),
               "Lubunga (Kasaï-Central)")
  expect_equal(.resolve_source_labels("Bili", "Nord-Ubangi", zones), "Bili (Nord-Ubangi)")

  # a colliding bare name whose province matches NEITHER candidate is left
  # unresolved — never folded into an arbitrary one (the upstream defect)
  expect_true(is.na(.resolve_source_labels("Lubunga", "Haut-Uele", zones)))
  expect_true(is.na(.resolve_source_labels("Lubunga", NA, zones)))

  # committed label alias (province-verified upstream), then data/aliases.csv
  expect_equal(.resolve_source_labels("Central Massa", "Kongo", zones), "Massa")
  expect_equal(.resolve_source_labels("Foo", "X", zones, c(Foo = "Bunia")), "Bunia")

  # unambiguous case/punctuation-insensitive last resort
  expect_equal(.resolve_source_labels("NGandajika", "Lomami", zones), "Ngandajika")

  # genuinely unknown stays NA
  expect_true(is.na(.resolve_source_labels("Nowhere", "Nowhere", zones)))

  # vectorised, and order preserved
  expect_equal(.resolve_source_labels(c("Kilo", "Nowhere", "Bunia"),
                                      c("Ituri", "X", "Ituri"), zones),
               c("Kilo", NA, "Bunia"))
})

# =============================================================================
# 10. compose_epicentre(fill=) — refilling cells the SOURCE could not observe
# =============================================================================

.fill_fixture <- function() {
  zones <- c("A", "B", "C", "D", "E")          # A = the source origin
  # Base kernel: row A spreads over B..E; other rows are arbitrary but valid.
  base <- matrix(0, 5, 5, dimnames = list(zones, zones))
  base["A", c("B", "C", "D", "E")] <- c(0.4, 0.3, 0.2, 0.1)
  base["B", c("A", "C")] <- c(0.5, 0.5)
  base["C", c("A", "B")] <- c(0.5, 0.5)
  base["D", c("A")] <- 1
  base["E", c("A")] <- 1
  # Source kernel: origin A measured only D and E (B is a co-origin it cannot
  # observe; C was never measured).
  epi <- matrix(0, 5, 5, dimnames = list(zones, zones))
  epi["A", c("D", "E")] <- c(0.75, 0.25)
  attr(epi, "measured")       <- list(A = c("D", "E"))
  attr(epi, "source_origins") <- list(A = c("A", "B"))
  list(zones = zones, base = base, epi = epi)
}

test_that("compose_epicentre(fill='none') is the plain row replacement", {
  skip_if_missing("compose_epicentre")
  f <- .fill_fixture()
  W <- compose_epicentre(f$base, f$epi, "A", f$zones, "T", fill = "none")
  expect_equal(unname(W["A", ]), c(0, 0, 0, 0.75, 0.25))
  expect_equal(unname(W["B", ]), unname(f$base["B", ]))   # non-source rows untouched
})

test_that("compose_epicentre(fill='origins') refills only the source's own origins", {
  skip_if_missing("compose_epicentre")
  f <- .fill_fixture()
  W <- compose_epicentre(f$base, f$epi, "A", f$zones, "T", fill = "origins")
  q <- f$base["A", "B"]                       # U = {B}; q = base mass on U = 0.4
  expect_equal(unname(W["A", "B"]), q)                     # filled from the base
  expect_equal(unname(W["A", "C"]), (1 - q) * 0)           # measured-zero stays zero
  expect_equal(unname(W["A", "D"]), (1 - q) * 0.75)        # data rescaled by (1-q)
  expect_equal(unname(W["A", "E"]), (1 - q) * 0.25)
  expect_equal(sum(W["A", ]), 1)
  expect_equal(W["A", "A"], 0)
})

test_that("compose_epicentre(fill='unmeasured') also refills never-measured zones", {
  skip_if_missing("compose_epicentre")
  f <- .fill_fixture()
  W <- compose_epicentre(f$base, f$epi, "A", f$zones, "T", fill = "unmeasured")
  q <- sum(f$base["A", c("B", "C")])          # U = {B, C}
  expect_equal(unname(W["A", "B"]), unname(f$base["A", "B"]))
  expect_equal(unname(W["A", "C"]), unname(f$base["A", "C"]))
  expect_equal(unname(W["A", "D"]), (1 - q) * 0.75)
  expect_equal(sum(W["A", ]), 1)
})

test_that("compose_epicentre: the q cap RESCALES the filled block so the row still sums to 1", {
  skip_if_missing("compose_epicentre")
  f <- .fill_fixture()
  # Base row A puts 95% of its mass on the refilled set {B, C}.
  f$base["A", ] <- 0
  f$base["A", c("B", "C", "D", "E")] <- c(0.6, 0.35, 0.03, 0.02)
  W <- compose_epicentre(f$base, f$epi, "A", f$zones, "T", fill = "unmeasured",
                         q_max = 0.5)
  expect_equal(sum(W["A", ]), 1)
  expect_equal(sum(W["A", c("B", "C")]), 0.5)              # capped total
  # ...and the cap preserves the base's RELATIVE split inside the filled block
  expect_equal(unname(W["A", "B"] / W["A", "C"]), 0.6 / 0.35)
  expect_equal(unname(W["A", "D"]), 0.5 * 0.75)            # data share is 1 - 0.5
})

test_that("compose_epicentre: degenerate inputs fall back rather than producing a bad row", {
  skip_if_missing("compose_epicentre")
  f <- .fill_fixture()

  # (a) the source has no mass outside the refilled set -> keep the base row
  f2 <- f; f2$epi["A", ] <- 0; f2$epi["A", "B"] <- 1
  attr(f2$epi, "measured")       <- list(A = "B")
  attr(f2$epi, "source_origins") <- list(A = c("A", "B"))
  W <- suppressWarnings(compose_epicentre(f2$base, f2$epi, "A", f2$zones, "T",
                                          fill = "origins"))
  expect_equal(unname(W["A", ]), unname(f2$base["A", ]))

  # (b) no metadata on the source kernel -> no fill, no error
  f3 <- f; attr(f3$epi, "measured") <- NULL; attr(f3$epi, "source_origins") <- NULL
  W3 <- suppressWarnings(compose_epicentre(f3$base, f3$epi, "A", f3$zones, "T",
                                           fill = "origins"))
  expect_equal(unname(W3["A", ]), c(0, 0, 0, 0.75, 0.25))

  # (c) an all-zero base row -> q = 0, so the source profile survives unchanged
  f4 <- f; f4$base["A", ] <- 0
  W4 <- compose_epicentre(f4$base, f4$epi, "A", f4$zones, "T", fill = "origins")
  expect_equal(unname(W4["A", ]), c(0, 0, 0, 0.75, 0.25))
})

test_that("compose_epicentre(fill=) keeps every mobility-matrix invariant", {
  skip_if_missing("compose_epicentre")
  skip_if_missing("assert_mobility_matrix")
  f <- .fill_fixture()
  for (lv in c("none", "origins", "unmeasured")) {
    W <- compose_epicentre(f$base, f$epi, "A", f$zones, "T", fill = lv)
    expect_silent(assert_mobility_matrix(W, f$zones, "T"))
  }
})

# =============================================================================
# 11. Left-censored gravity fit — the three-state form reduces to the two-state
# =============================================================================

test_that(".fit_gravity_censored: an explicit all-zeros mask with cens_lo=0 equals the implicit fit", {
  skip_if_missing(".fit_gravity_censored")
  set.seed(42)
  n  <- 200
  td <- data.frame(pop_i = runif(n, 1e4, 1e6), pop_j = runif(n, 1e4, 1e6),
                   dist_ij = runif(n, 10, 600))
  mu <- exp(-8 + 0.5 * log(td$pop_i) + 0.6 * log(td$pop_j) - 1.2 * log(td$dist_ij + 1))
  td$flow <- rpois(n, mu)
  td$flow[td$flow > 0 & td$flow < 15] <- 0          # mimic the suppression pattern
  fml <- flow ~ log(pop_i) + log(pop_j) + log(dist_ij + 1)
  start <- suppressWarnings(glm(fml, family = poisson(), data = td))
  a <- .fit_gravity_censored(td, fml, 15L, start, "T")                       # implicit
  b <- .fit_gravity_censored(td, fml, 15L, start, "T",
                             cens_vec = td$flow < 15L, cens_lo = 0L)         # explicit
  skip_if(is.null(a) || is.null(b), "censored fit did not converge on the fixture")
  expect_equal(a$coefficients, b$coefficients, tolerance = 1e-8)
  expect_equal(a$loglik, b$loglik, tolerance = 1e-8)
  # cens_lo = 1 removes the P(Y = 0) mass from each censored cell's interval, so the
  # same data can never score HIGHER under the narrower interval.
  cc <- .fit_gravity_censored(td, fml, 15L, start, "T",
                              cens_vec = td$flow < 15L, cens_lo = 1L)
  skip_if(is.null(cc), "three-state fit did not converge on the fixture")
  expect_lt(cc$loglik, a$loglik + 1e-6)
})

# =============================================================================
# 12. summarise_matrix() — per-origin columns are keyed on the ORIGIN names
# =============================================================================

test_that("summarise_matrix: top-5 columns are named from the origins and are populated", {
  skip_if_missing("summarise_matrix")
  zones <- c("Bunia", "Mongbwalu", "Rwampara", "Lita")
  W <- matrix(0, 4, 4, dimnames = list(zones, zones))
  W["Bunia", c("Lita", "Rwampara")]     <- c(0.7, 0.3)
  W["Mongbwalu", c("Lita", "Bunia")]    <- c(0.6, 0.4)
  W["Rwampara", c("Lita")]              <- 1
  W["Lita", c("Bunia")]                 <- 1
  s <- summarise_matrix(W, "T", zones, c("Bunia", "Mongbwalu", "Rwampara"))
  expect_true(all(c("top5_Bunia", "top5_Mongbwalu", "top5_Rwampara") %in% names(s)))
  expect_false(any(vapply(s[grep("^top5_", names(s))], function(x) is.na(x[1]), logical(1))))
  expect_false("top5_Mongbalu" %in% names(s))   # the pre-2026-07 spelling is gone
})

# =============================================================================
# 13. Built kernels — regressions on the two data defects (skip without a build)
# =============================================================================

.mob_built <- function(id) {
  p <- file.path(get0("OUT_MOBILITY", ifnotfound = ""), sprintf("mobility_%s.rds", id))
  if (nzchar(p) && file.exists(p)) readRDS(p) else NULL
}

test_that("built kernels: Kilo is reachable from the epicentre (name-resolution regression)", {
  W1 <- .mob_built("M1")
  skip_if(is.null(W1), "mobility matrices not built in this environment")
  skip_if(!all(c("Bunia", "Kilo") %in% rownames(W1)), "zone spine does not carry these zones")
  expect_gt(W1["Bunia", "Kilo"], 0)
  for (id in c("M8", "M13", "M14")) {
    W <- .mob_built(id)
    if (!is.null(W)) expect_gt(W["Bunia", "Kilo"], 0)
  }
})

test_that("built kernels: the fill variants restore the provider-excluded origin cells", {
  Wf <- .mob_built("M14-fill"); W0 <- .mob_built("M14")
  skip_if(is.null(Wf) || is.null(W0), "fill variants not built in this environment")
  pairs <- list(c("Beni", "Butembo"), c("Makiso Kisangani", "Mangobo"))
  for (p in pairs) {
    if (all(p %in% rownames(W0))) {
      expect_equal(unname(W0[p[1], p[2]]), 0)   # classic kernel asserts zero
      expect_gt(Wf[p[1], p[2]], 0)              # fill variant does not
    }
  }
  expect_equal(unname(rowSums(Wf)[rowSums(Wf) > 0]),
               rep(1, sum(rowSums(Wf) > 0)), tolerance = 1e-9)
})

test_that("built kernels: every manifest entry has a kernel file on disk", {
  man_p <- file.path(get0("OUT_MOBILITY", ifnotfound = ""), "mobility_manifest.csv")
  skip_if(!nzchar(man_p) || !file.exists(man_p), "no manifest in this environment")
  man <- readr::read_csv(man_p, show_col_types = FALSE)
  fs  <- sub("^mobility_(.*)\\.rds$", "\\1",
             list.files(dirname(man_p), pattern = "^mobility_.*\\.rds$"))
  # SUBSET, not set equality: the manifest exists precisely so a consumer can ignore
  # kernel files left behind by an earlier run (a kernel since switched off, or from
  # another branch), so extra .rds files on disk are the expected case, not a failure.
  expect_true(all(man$matrix_id %in% fs),
              info = paste("manifest entries with no file:",
                           paste(setdiff(man$matrix_id, fs), collapse = ", ")))
  expect_gt(nrow(man), 0L)
})
