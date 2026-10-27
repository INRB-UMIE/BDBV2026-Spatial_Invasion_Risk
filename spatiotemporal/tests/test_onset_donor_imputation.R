# =============================================================================
# test_onset_donor_imputation.R — missing-onset imputation from same-class, same-week donors
# (01_data_prep.R: .impute_onset_delays(), .uniform_from_ids()). Synthetic data only.
# =============================================================================
library(testthat)
suppressMessages(source(file.path(here::here(), "spatiotemporal", "00_config.R")))
suppressMessages(source(file.path(here::here(), "spatiotemporal", "01_data_prep.R")))

wk0 <- as.Date("2026-06-06")                          # a Saturday: a week start under week_start = 6

test_that("the per-record uniform is deterministic, in (0, 1), order-independent and seed-dependent", {
  ids <- sprintf("A%03d", 1:200)
  u1 <- .uniform_from_ids(ids, 1L)
  expect_true(all(u1 > 0 & u1 < 1))
  expect_identical(u1, .uniform_from_ids(ids, 1L))
  o <- sample(200)
  expect_identical(.uniform_from_ids(ids[o], 1L), u1[o])
  expect_false(isTRUE(all.equal(u1, .uniform_from_ids(ids, 2L))))
})

test_that("a recipient draws only from donors of its own class sampled in its own week", {
  # week 1: confirmed donors with intervals 1..30; week 2: confirmed donors all 50;
  # non-confirmed donors in week 1 all 99.
  samp  <- c(rep(wk0 + 3, 30), rep(wk0 + 10, 30), rep(wk0 + 3, 30), wk0 + 4, wk0 + 11, wk0 + 5)
  rec   <- c(samp[1:30] - 1:30, samp[31:60] - 50, samp[61:90] - 99, NA, NA, NA)
  grp   <- c(rep(TRUE, 60), rep(FALSE, 30), TRUE, TRUE, FALSE)
  need  <- c(rep(FALSE, 90), TRUE, TRUE, TRUE)
  usable <- !need
  out <- .impute_onset_delays(samp, rec, usable, grp, need, ids = sprintf("r%d", 1:93), seed = 1L,
                              week_start = 6L, min_donors = 20L)
  expect_true(all(is.na(out[1:90])))
  expect_true(out[91] %in% 1:30)          # confirmed, week 1
  expect_identical(out[92], 50L)          # confirmed, week 2
  expect_identical(out[93], 99L)          # not confirmed, week 1
})

test_that("the draw is the empirical quantile at u, so it is a uniform pick among the donors", {
  samp <- c(rep(wk0 + 2, 25), wk0 + 2)
  rec  <- c(samp[1:25] - (0:24), NA)
  ids  <- "only-recipient"
  u <- .uniform_from_ids(ids, 3L)
  out <- .impute_onset_delays(samp, rec, c(rep(TRUE, 25), FALSE), rep(TRUE, 26), c(rep(FALSE, 25), TRUE),
                              ids = c(sprintf("d%d", 1:25), ids), seed = 3L, week_start = 6L, min_donors = 20L)
  expect_identical(out[26], as.integer((0:24)[floor(u * 25) + 1L]))
})

test_that("a thin week borrows adjacent weeks symmetrically until it has min_donors", {
  # recipient week has 2 donors (interval 7); the week before has 10 (interval 3), the week after
  # 10 (interval 5); two weeks away 50 (interval 40) must NOT be reached with min_donors = 20.
  samp <- c(rep(wk0 + 1, 2), rep(wk0 - 6, 10), rep(wk0 + 8, 10), rep(wk0 + 15, 50), wk0 + 2)
  ivl  <- c(rep(7, 2), rep(3, 10), rep(5, 10), rep(40, 50))
  rec  <- c(samp[1:72] - ivl, NA)
  need <- c(rep(FALSE, 72), TRUE)
  draws <- vapply(1:200, function(k)
    .impute_onset_delays(samp, rec, !need, rep(TRUE, 73), need, ids = c(sprintf("d%d", 1:72), "x"),
                         seed = k, week_start = 6L, min_donors = 20L)[73], integer(1))
  expect_true(all(draws %in% c(3L, 5L, 7L)))
  expect_true(all(c(3L, 5L, 7L) %in% draws))
})

test_that("records sampled after the as-of date are not donors; no donor of a class gives NA", {
  samp <- c(rep(wk0 + 1, 25), rep(wk0 + 2, 25), wk0 + 3, wk0 + 3)
  rec  <- c(samp[1:25] - 4, samp[26:50] - 9, NA, NA)
  grp  <- c(rep(TRUE, 50), TRUE, FALSE)
  need <- c(rep(FALSE, 50), TRUE, TRUE)
  out <- .impute_onset_delays(samp, rec, !need, grp, need, ids = sprintf("r%d", 1:52), seed = 1L,
                              asof = wk0 + 1, week_start = 6L, min_donors = 20L)
  expect_identical(out[51], 4L)           # the wk0 + 2 donors postdate the as-of date
  expect_true(is.na(out[52]))             # no non-confirmed donor exists
})

test_that("the draw of a record does not move when other records are added or rows reordered", {
  samp <- c(rep(wk0 + 2, 40), wk0 + 3, wk0 + 4)
  rec  <- c(samp[1:40] - rep(0:9, 4), NA, NA)
  need <- c(rep(FALSE, 40), TRUE, TRUE)
  ids  <- c(sprintf("d%d", 1:40), "x", "y")
  a <- .impute_onset_delays(samp, rec, !need, rep(TRUE, 42), need, ids, seed = 5L, week_start = 6L)
  o <- rev(seq_along(samp))
  b <- .impute_onset_delays(samp[o], rec[o], (!need)[o], rep(TRUE, 42), need[o], ids[o], seed = 5L, week_start = 6L)
  expect_identical(b[o], a)
  c2 <- .impute_onset_delays(samp[-42], rec[-42], (!need)[-42], rep(TRUE, 41), need[-42], ids[-42],
                             seed = 5L, week_start = 6L)
  expect_identical(c2[41], a[41])
})

test_that("a missing id is refused, not hashed: every NA would share ONE uniform", {
  expect_error(.uniform_from_ids(c("a", NA), 1L), "missing id")
  samp <- c(rep(wk0 + 2, 25), wk0 + 3)
  rec  <- c(samp[1:25] - (0:24), NA)
  need <- c(rep(FALSE, 25), TRUE)
  expect_error(.impute_onset_delays(samp, rec, !need, rep(TRUE, 26), need, ids = NULL, seed = 1L,
                                    week_start = 6L), "one entry per record")
})

test_that(".onset_draw_keys keeps alert_id and gives each id-less record its own order-free key", {
  df <- data.frame(alert_id = c("a1", NA, NA, "a2", NA), health_zone = c("Z", "Z", "Z", "Y", "Z"),
                   date_of_notification = c("2026-06-01", "2026-06-02", "2026-06-03", NA, "2026-06-02"),
                   province = c("P", "P", "P", "Q", NA), stringsAsFactors = FALSE)
  k <- .onset_draw_keys(df)
  expect_identical(k[c(1, 4)], c("a1", "a2"))
  expect_false(anyNA(k))
  expect_true(all(startsWith(k[c(2, 3, 5)], "NOID-")))
  expect_identical(anyDuplicated(k), 0L)                      # distinct records, distinct keys
  o <- c(5, 3, 1, 4, 2)
  expect_identical(.onset_draw_keys(df[o, ]), k[o])           # row order does not move a key
  expect_identical(.onset_draw_keys(df[-2, ]), k[-2])         # nor do other records
  # the key ignores a previously attached key column (it is computed on the extracted fields)
  df2 <- df; df2$.onset_draw_key <- "stale"
  expect_identical(.onset_draw_keys(df2), k)
  # id-less records in one donor week now take their own draws
  u <- .uniform_from_ids(k, 1L)
  expect_identical(anyDuplicated(u), 0L)
})
