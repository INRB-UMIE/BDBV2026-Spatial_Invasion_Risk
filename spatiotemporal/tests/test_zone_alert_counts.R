# =============================================================================
# tests/test_zone_alert_counts.R — the shared per-zone alert definition
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# zone_alert_counts() (01_data_prep.R) is the ONE definition of "alerts" and "laboratory-
# resolved alerts" behind two published products: the ascertainment appendix
# (45_ascertainment_figures.R -> zone_ascertainment.csv) and panel D of
# Figure3_rank_surveillance (make_manuscript_figures.R). Pinned here, on a synthetic line
# list, are the rules either product would silently change if they drifted:
#   * the sitrep top-up rows are NOT alerts (they carry no test and no laboratory date);
#   * a record is in the window when its case date (date_index) OR its notification date is,
#     inclusive at BOTH ends — so a notification typo does not drop a dated alert, and a
#     pre-outbreak (or impossible) onset does not drop an alert notified in the window;
#   * the records load_linelist() leaves out (no in-outbreak onset or specimen date) are
#     still alerts, placed by notification;
#   * laboratory-resolved = confirmed + not-a-case, and nothing else (suspected, probable and
#     never-classified alerts are unresolved);
#   * a zone with no alert has no row (callers coalesce to 0);
#   * the zone an alert is counted in: a bare name the spine carries only as "Name (Province)"
#     (Bili, Lubunga) is resolved by the record's province (.resolve_zone_by_province), or the
#     alerts of Bili (Bas-Uele) are counted in no zone at all.
# =============================================================================

test_that(".record_in_outbreak is load_linelist()'s keep rule", {
  if (!exists(".record_in_outbreak", mode = "function")) testthat::skip(".record_in_outbreak() unavailable")
  o <- OUTBREAK_START
  on <- c(o, o - 1, NA, NA, o - 1, NA)
  sa <- c(NA, o, o, o - 1, o - 1, NA)
  expect_identical(.record_in_outbreak(on, sa), c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("zone_alert_counts applies the either-date window, the sitrep exclusion and the resolution rule", {
  if (!exists("zone_alert_counts", mode = "function")) testthat::skip("zone_alert_counts() unavailable")
  from <- as.Date("2026-01-10"); to <- as.Date("2026-01-20")
  N <- as.Date("2026-01-15")      # an in-window notification date
  ll <- tibble::tribble(
    ~health_zone, ~date_index,            ~date_of_notification,  ~alert_id,          ~final_mve_case_classification,
    "A",          as.Date("2026-01-10"),  N,                      "a1",               CONFIRMED_STATUS,    # case date on the lower bound
    "A",          as.Date("2026-01-20"),  N,                      "a2",               NOT_A_CASE_STATUS,   # case date on the upper bound
    "A",          as.Date("2026-01-15"),  N,                      "a3",               SUSPECTED_STATUS,    # unresolved
    "A",          as.Date("2026-01-15"),  N,                      "a4",               NA_character_,       # unresolved
    "A",          as.Date("2026-01-15"),  N,                      "a5",               "probable_case",     # unresolved
    "A",          as.Date("2026-01-15"),  as.Date("2027-01-15"),  "a6",               CONFIRMED_STATUS,    # notification typo: kept by case date
    "A",          as.Date("2026-01-15"),  N,                      "SITREP-CONF-A-01", CONFIRMED_STATUS,    # top-up: excluded
    "A",          as.Date("2026-01-09"),  as.Date("2026-01-09"),  "a7",               NOT_A_CASE_STATUS,   # before window on both
    "A",          as.Date("2026-01-21"),  as.Date("2026-01-21"),  "a8",               CONFIRMED_STATUS,    # after window on both
    "A",          as.Date("2026-01-21"),  as.Date("2026-01-19"),  "a9",               NOT_A_CASE_STATUS,   # case date out, notified in: kept
    NA_character_, as.Date("2026-01-15"), N,                      "a10",              CONFIRMED_STATUS,    # no zone
    "B",          as.Date("2026-01-12"),  N,                      NA_character_,      NOT_A_CASE_STATUS,   # NA id is a real alert
    "C",          as.Date("2026-01-12"),  N,                      "SITREP-CONF-C-01", CONFIRMED_STATUS)    # only a top-up

  # Records load_linelist() leaves out: placed by notification only.
  unplaced <- tibble::tribble(
    ~health_zone, ~date_of_notification, ~date_of_symptom_onset,  ~date_of_sample_collection, ~alert_id, ~final_mve_case_classification,
    "A",          as.Date("2026-01-10"), as.Date(NA),             as.Date(NA),                "u1",      NA_character_,   # no dates; on the lower bound
    "A",          as.Date("2026-01-21"), as.Date(NA),             as.Date(NA),                "u2",      NA_character_,   # notified after window
    "A",          as.Date(NA),           as.Date(NA),             as.Date(NA),                "u3",      NA_character_,   # unplaceable
    "A",          as.Date("2026-01-12"), OUTBREAK_START - 400,    as.Date(NA),                "u4",      NA_character_,   # pre-outbreak onset, notified in window
    "D",          as.Date("2026-01-20"), as.Date(NA),             as.Date(NA),                "u5",      NA_character_)   # zone with only an unplaced alert

  z <- zone_alert_counts(ll, unplaced, from = from, to = to)
  expect_identical(z$health_zone, c("A", "B", "D"))   # C has only a top-up row, so no row at all
  expect_identical(attr(z, "n_line_list"), 8L)        # a1-a6, a9, B's alert
  expect_identical(attr(z, "n_unplaced"), 3L)         # u1, u4, u5
  a <- z[z$health_zone == "A", ]
  expect_identical(a$alerts_total, 9L)                # a1-a6, a9 + u1, u4
  expect_identical(a$n_confirmed, 2L)                 # a1, a6
  expect_identical(a$n_not_a_case, 2L)                # a2, a9
  expect_identical(a$alerts_resolved, 4L)
  b <- z[z$health_zone == "B", ]
  expect_identical(c(b$alerts_total, b$n_confirmed, b$n_not_a_case, b$alerts_resolved),
                   c(1L, 0L, 1L, 1L))
  d <- z[z$health_zone == "D", ]
  expect_identical(c(d$alerts_total, d$alerts_resolved), c(1L, 0L))
  expect_true(all(z$alerts_resolved == z$n_confirmed + z$n_not_a_case))
  expect_true(all(z$alerts_resolved <= z$alerts_total))
})

test_that("zone_alert_counts refuses a bad window, a missing column and a kept record passed as unplaced", {
  if (!exists("zone_alert_counts", mode = "function")) testthat::skip("zone_alert_counts() unavailable")
  ll <- tibble::tibble(health_zone = "A", date_index = as.Date("2026-01-15"),
                       date_of_notification = as.Date("2026-01-15"),
                       alert_id = "a1", final_mve_case_classification = CONFIRMED_STATUS)
  und <- tibble::tibble(health_zone = "A", date_of_notification = as.Date("2026-01-15"),
                        date_of_symptom_onset = as.Date(NA), date_of_sample_collection = as.Date(NA),
                        alert_id = "u1", final_mve_case_classification = NA_character_)
  f <- as.Date("2026-01-10"); t <- as.Date("2026-01-20")
  expect_error(zone_alert_counts(ll, und, from = NA, to = t), "from <= to")
  expect_error(zone_alert_counts(ll, und, from = t, to = f), "from <= to")
  expect_error(zone_alert_counts(ll[, -4], und, from = f, to = t))
  # A record load_linelist() KEEPS, passed again as unplaced, would be counted twice.
  und_bad <- und; und_bad$date_of_symptom_onset <- OUTBREAK_START
  expect_error(zone_alert_counts(ll, und_bad, from = f, to = t), "load_unplaced_alerts")
})

test_that("a same-named zone's bare name resolves by the record's province, and only then", {
  if (!exists(".resolve_zone_by_province", mode = "function"))
    testthat::skip(".resolve_zone_by_province() unavailable")
  spine <- c("Bili (Bas-Uele)", "Bili (Nord-Ubangi)", "Lubunga (Kasa\u00ef-Central)", "Lubunga (Tshopo)",
             "Bunia", "Foo", "Foo (Bar)")
  zone <- c("Bili", "Bili", "Bili", "Bili", "Lubunga", "Bunia", "Bili (Bas-Uele)", NA, "Foo", "Elsewhere", "Bili")
  prov <- c("Bas-Uele", "nord ubangi", NA, "Tshopo", "Kasai-Central", "Ituri", "Ituri", "Bas-Uele", "Bar",
            "Bas-Uele", " BAS-UELE ")
  expect_identical(.resolve_zone_by_province(zone, prov, spine),
                   c("Bili (Bas-Uele)", "Bili (Nord-Ubangi)",
                     "Bili",                   # no province: left for the aliases / the spine warning
                     "Bili",                   # a province matching neither entry: never folded into one
                     "Lubunga (Kasa\u00ef-Central)",  # accent/punctuation-insensitive province match
                     "Bunia", "Bili (Bas-Uele)",  # names the spine carries as written are untouched
                     NA,
                     "Foo",                    # a bare name that is itself a spine zone stays as written
                     "Elsewhere", "Bili (Bas-Uele)"))
  expect_identical(.resolve_zone_by_province(zone, NULL, spine), zone)       # no province column
  expect_identical(.resolve_zone_by_province(zone, prov, NULL), zone)        # no spine
})
