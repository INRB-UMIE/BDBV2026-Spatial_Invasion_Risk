# =============================================================================
# tests/test_ascertainment_alert_split.R — the per-source / mapped-zone alert split in 45
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# 45_ascertainment_figures.R reports alert totals over the MAPPED zones only (zones in the
# 519-zone spine), and splits them into line-list records and unplaced (notification-only)
# records by calling the shared zone_alert_counts() once per source, with the other source
# passed as an EMPTY frame. Pinned here, on a synthetic line list (no real records):
#   * an empty half is accepted and contributes nothing;
#   * the two single-source counts add up, zone by zone, to the combined count;
#   * mapped + unmapped = all alerts, where "unmapped" is a zone name absent from the spine.
#     The synthetic unmapped name is a bare "Bili", which the spine splits into two
#     disambiguated zones; on real data load_linelist() resolves such names by province, so in
#     45 the unmapped branch is a guard, and zone_alert_counts() itself does no name resolution.
# =============================================================================

test_that("zone_alert_counts splits by source and by mapped zone without losing an alert", {
  if (!exists("zone_alert_counts", mode = "function")) testthat::skip("zone_alert_counts() unavailable")
  from <- as.Date("2026-01-10"); to <- as.Date("2026-01-20"); N <- as.Date("2026-01-15")
  ll <- tibble::tribble(
    ~health_zone, ~date_index,           ~date_of_notification, ~alert_id,          ~final_mve_case_classification,
    "A",          as.Date("2026-01-12"), N,                     "a1",               CONFIRMED_STATUS,
    "A",          as.Date("2026-01-13"), N,                     "a2",               NOT_A_CASE_STATUS,
    "Bili",       as.Date("2026-01-14"), N,                     "b1",               NOT_A_CASE_STATUS,   # unmapped name
    "Bili",       as.Date("2026-01-14"), N,                     "b2",               NA_character_,       # unmapped name
    "C",          as.Date("2026-01-14"), N,                     "SITREP-CONF-C-01", CONFIRMED_STATUS)    # top-up: never an alert
  unplaced <- tibble::tribble(
    ~health_zone, ~date_of_notification, ~date_of_symptom_onset, ~date_of_sample_collection, ~alert_id, ~final_mve_case_classification,
    "A",          N,                     as.Date(NA),            as.Date(NA),                "u1",      NA_character_,
    "Bili",       N,                     as.Date(NA),            as.Date(NA),                "u2",      NA_character_)
  spine_zones <- c("A", "Bili (Bas-Uele)", "Bili (Nord-Ubangi)", "C")

  both <- zone_alert_counts(ll, unplaced, from = from, to = to)
  only_ll <- zone_alert_counts(ll, unplaced[0, , drop = FALSE], from = from, to = to)
  only_un <- zone_alert_counts(ll[0, , drop = FALSE], unplaced, from = from, to = to)

  # an empty half contributes nothing
  expect_identical(attr(only_ll, "n_unplaced"), 0L)
  expect_identical(attr(only_un, "n_line_list"), 0L)
  expect_identical(attr(only_ll, "n_line_list"), attr(both, "n_line_list"))
  expect_identical(attr(only_un, "n_unplaced"), attr(both, "n_unplaced"))

  # zone by zone, the single-source counts add up to the combined count
  tot <- function(z, k) sum(z$alerts_total[z$health_zone == k])
  for (k in both$health_zone) expect_identical(tot(only_ll, k) + tot(only_un, k), tot(both, k))

  # mapped + unmapped = all alerts; the unmapped zone is exactly the ambiguous bare name
  unmapped <- setdiff(both$health_zone, spine_zones)
  expect_identical(unmapped, "Bili")
  mapped_total <- function(z) sum(z$alerts_total[z$health_zone %in% spine_zones])
  n_unmapped <- sum(both$alerts_total[both$health_zone %in% unmapped])
  expect_identical(n_unmapped, 3L)                                   # b1, b2, u2
  expect_identical(mapped_total(only_ll) + mapped_total(only_un), 3L)   # a1, a2, u1
  expect_identical(mapped_total(only_ll) + mapped_total(only_un) + n_unmapped,
                   sum(both$alerts_total))
})

test_that("the forecast origin used for 'already notified' is cutoff + LFO_ORIGIN_OFFSET_DAYS", {
  if (!exists("lfo_origin", mode = "function")) testthat::skip("lfo_origin() unavailable")
  cu <- as.Date("2026-07-04")
  expect_identical(lfo_origin(cu), cu + LFO_ORIGIN_OFFSET_DAYS)
  # an alert notified after the cutoff but on or before the origin counts as already notified
  alert_date <- cu + c(0L, LFO_ORIGIN_OFFSET_DAYS, LFO_ORIGIN_OFFSET_DAYS + 1L)
  expect_identical(alert_date <= lfo_origin(cu), c(TRUE, TRUE, FALSE))
})
