# =============================================================================
# tests/test_spread_kinematics.R — the envelope-saturation verdict of 43_spread_kinematics.R
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# spread_kinematics_summary.csv publishes a "saturated / not saturated" verdict on the spread
# envelope. It used to be computed on the RUNNING MAXIMUM R_max(t), which ONE distant zone
# sets: on the 2026-10-09 frame a single 1,606 km seeding at week 18 made week 18 "the week
# the extent reached 90%", made the speed "to 90%" that zone's distance divided by 18, and made
# the slope "after 90%" zero by construction — and the word SATURATED was typed into the CSV.
# compute_kinematics() now measures the envelope on the 90th-percentile reach R_p90(t), fits
# both phases by OLS, needs >= 3 fit-window weeks on each side, and derives the verdict text.
#
# These tests drive compute_kinematics() on SYNTHETIC geography (no line list, no shapefile):
# zones on a line east of the index zone, so every distance is known by construction.
# =============================================================================

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr); library(sf); library(geosphere)
}))

# Only 43's top-level ASSIGNMENTS are evaluated (its functions and its parameter block): the
# script also re-sources 00_config.R, 01 and 03 and its analysis sits inside a guarded `if`,
# none of which is needed or wanted here.
.k43 <- local({
  f <- file.path(here::here(), "spatiotemporal", "43_spread_kinematics.R")
  if (!file.exists(f)) return(NULL)
  e <- new.env(parent = globalenv())
  for (ex in parse(f, keep.source = FALSE))
    # try(): a path assignment such as KEY_DIR needs the config, which these tests do not
    if (is.call(ex) && identical(ex[[1]], as.name("<-"))) try(eval(ex, envir = e), silent = TRUE)
  e
})

# n zones on a parallel at 1N, the index zone Z01 at 30E; road = 1.3 x great-circle, time at
# 50 km/h. `inv_week` gives each zone's first confirmed week (NA = never invaded).
.synthetic <- function(km_east, inv_week, n_weeks) {
  n <- length(km_east)
  zones <- sprintf("Z%02d", seq_len(n))
  lon <- 30 + km_east / 111.32; lat <- rep(1, n)
  crs <- .k43$.laea(30, 1)
  xy <- sf::st_coordinates(sf::st_transform(
    sf::st_sfc(lapply(seq_len(n), function(i) sf::st_point(c(lon[i], lat[i]))), crs = 4326), crs))
  gc <- outer(seq_len(n), seq_len(n), Vectorize(function(i, j)
    geosphere::distGeo(c(lon[i], lat[i]), c(lon[j], lat[j])) / 1000))
  dimnames(gc) <- list(zones, zones)
  geo <- list(ctr = tibble::tibble(health_zone = zones, province = "P", x = xy[, 1], y = xy[, 2],
                                   lon = lon, lat = lat),
              road = 1.3 * gc, time = 1.3 * gc / 50, crs_laea = crs)
  weeks <- as.Date("2026-05-02") + 7 * (0:(n_weeks - 1L))
  zw <- tidyr::expand_grid(health_zone = zones, week_start = weeks) %>%
    dplyr::mutate(wk = as.integer((week_start - weeks[1]) / 7),
                  iw = inv_week[match(health_zone, zones)],
                  # the index zone carries the seed-window cases; every other zone one case a
                  # week from its invasion week on
                  confirmed = dplyr::case_when(health_zone == "Z01" ~ 5L,
                                               !is.na(iw) & wk >= iw ~ 1L, TRUE ~ 0L)) %>%
    dplyr::select(health_zone, week_start, confirmed)
  list(geo = geo, zw = zw, weeks = weeks)
}

test_that("one distant zone cannot set the envelope verdict, and the verdict is the data's", {
  if (is.null(.k43)) testthat::skip("43_spread_kinematics.R not present")
  # 20 near zones spread outward over weeks 1-10, then nearer zones keep arriving to week 20,
  # and ONE zone 1,500 km out arrives at week 18 — the running maximum's jump.
  km  <- c(0, seq(20, 400, length.out = 20), seq(10, 100, length.out = 9), 1500)
  inv <- c(0, 1:20 %/% 2 + 1, 11:19, 18)
  s <- .synthetic(km, inv, n_weeks = 22)
  K <- .k43$compute_kinematics(s$zw, s$geo, fit_cutoff = s$weeks[22], seed_cutoff = s$weeks[3])
  S <- K$stats
  expect_equal(K$epi_zone, "Z01")
  # The running maximum jumps at week 18 (the 1,500 km zone)...
  w <- K$weekly
  expect_gt(w$reach_max_km[w$week == 18], 1.5 * w$reach_max_km[w$week == 17])
  # ...but the p90 envelope reaches 90% of its final value BEFORE that week, so the outlier
  # neither dates the threshold nor makes the post-threshold slope zero by construction.
  expect_lt(S$week_p90_reach_90pct, 18)
  expect_gte(S$n_weeks_env_late, 3L)
  expect_true(is.finite(S$v_env_p90_early) && is.finite(S$v_env_p90_late))
  # The verdict is DERIVED: it names the two slopes it was decided on and follows the stated rule.
  sat <- abs(S$v_env_p90_late) <= 0.1 * abs(S$v_env_p90_early)
  expect_true(startsWith(S$envelope_verdict, if (sat) "SATURATED" else "NOT SATURATED"))
  expect_true(grepl(sprintf("%+.1f km/wk", S$v_env_p90_late), S$envelope_verdict, fixed = TRUE))
  # The post-threshold new-invasion summary sees the distant zone (the envelope hides it).
  expect_equal(S$newmax_max_after, max(1.3 * geosphere::distGeo(c(30, 1),
                                                               c(30 + 1500 / 111.32, 1)) / 1000),
               tolerance = 1e-6)
})

test_that("a steadily expanding envelope is NOT called saturated", {
  if (is.null(.k43)) testthat::skip("43_spread_kinematics.R not present")
  # Two zones a week, each pair further out than the last, for 30 weeks (long enough that
  # the post-threshold phase has the >= 3 weeks a slope needs).
  km  <- c(0, rep(seq(30, 900, length.out = 30), each = 2))
  inv <- c(0, rep(1:30, each = 2))
  s <- .synthetic(km, inv, n_weeks = 31)
  S <- .k43$compute_kinematics(s$zw, s$geo, fit_cutoff = s$weeks[31],
                               seed_cutoff = s$weeks[3])$stats
  expect_gte(S$n_weeks_env_late, 3L)
  expect_true(startsWith(S$envelope_verdict, "NOT SATURATED (still expanding)"))
  expect_gt(S$v_env_p90_late, 0.1 * S$v_env_p90_early)
})

test_that("too few weeks after the threshold gives UNDETERMINED, not a default verdict", {
  if (is.null(.k43)) testthat::skip("43_spread_kinematics.R not present")
  # The envelope only reaches 90% of its final value in the last fit-window week.
  km  <- c(0, seq(30, 900, length.out = 12))
  inv <- c(0, 1:12)
  s <- .synthetic(km, inv, n_weeks = 13)
  S <- .k43$compute_kinematics(s$zw, s$geo, fit_cutoff = s$weeks[13],
                               seed_cutoff = s$weeks[3])$stats
  expect_lt(S$n_weeks_env_late, 3L)
  expect_true(is.na(S$v_env_p90_late))
  expect_true(startsWith(S$envelope_verdict, "UNDETERMINED"))
})

# The completeness ECDF must CAP delays beyond DELAY_CAP_DAYS (count them as not yet arrived
# below the cap, arrived at it), not DROP them: dropping conditions the ECDF on D <= cap and
# over-states completeness at every follow-up shorter than the cap. The docstring always said
# "capped"; the body used to drop.
test_that("week_completeness caps long delays rather than conditioning them away", {
  if (is.null(.k43)) testthat::skip("43_spread_kinematics.R not present")
  asof <- as.Date("2026-09-30"); cap <- .k43$DELAY_CAP_DAYS
  onset <- asof - cap - 10L                    # complete follow-up out to the cap for everyone
  # 90 records known 5 days after onset, 10 known 40 days after (beyond the cap)
  del <- c(rep(5L, 90L), rep(40L, 10L))
  ll <- tibble::tibble(confirmed = TRUE, date_index = rep(onset, 100L),
                       lab_analysis_date = onset + del,
                       reporting_date = as.Date(NA), date_of_notification = as.Date(NA))
  # week A: every onset day has 14-20 days of follow-up (between 5 and the cap);
  # week B: every onset day has >= cap days of follow-up
  wA <- asof - 20L; wB <- asof - cap - 20L
  comp <- .k43$week_completeness(ll, c(wB, wA), asof, cap = cap)
  expect_equal(as.numeric(comp[2]), 0.9)       # 90 of 100 known by 14-20 d: the 10 long ones count
  expect_equal(as.numeric(comp[1]), 1)         # at the cap every delay counts as arrived
  expect_equal(attr(comp, "share_over_cap"), 0.1)
  expect_identical(attr(comp, "n_pairs"), 100L)
})
