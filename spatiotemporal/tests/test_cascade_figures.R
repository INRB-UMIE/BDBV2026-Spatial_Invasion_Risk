# =============================================================================
# test_cascade_figures.R - regression tests for the cascade figure scripts (35/37/38/39/40/42)
# =============================================================================
# Each block pins a property that a 2026-10 review found violated:
#   * a cascade product is drawn only when it is the latest COMPLETED run's (or the live run's)
#     for this analysis date, and an uncalibrated run says so                     [stamp]
#   * the reach CSV carries the run id / analysis date / calibration flag the check reads
#   * Panel B's label survives a fractional (MI-mean) count, and .panel() catches render-time
#     errors rather than letting them abort save_dual() mid-write                  [label]
#   * Panel C's province palette has no shared hue                                 [palette]
#   * the seed-zone provenance is written from the RESOLVED hubs only               [hubs]
#   * the 39 capture curve is the per-origin mean recall, with English month labels [39]
#   * 42 draws the sweep only when its sidecar matches the stamp and the configuration [42]
# None of these tests writes outside tempdir().
# =============================================================================

.ST <- file.path(here::here(), "spatiotemporal")
if (!exists("ANALYSIS_DATE"))
  suppressMessages(suppressWarnings(source(file.path(.ST, "00_config.R"))))
if (!exists("prov_palette", mode = "function"))
  suppressMessages(source(file.path(.ST, "figure_style.R")))
if (!exists("cascade_stamp_check", mode = "function"))
  suppressMessages(source(file.path(.ST, "35_cascade_viz.R")))

# Pull ONE top-level (or guard-nested) `name <- function(...)` definition out of a script without
# running the script: 38 and 40 execute on source (they read data at top level).
.def_from <- function(file, name, envir = new.env(parent = globalenv())) {
  ex <- parse(file.path(.ST, file), keep.source = FALSE)
  hit <- NULL
  walk <- function(e) {
    if (!is.null(hit) || !is.call(e)) return(invisible())
    if (identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(name))) {
      hit <<- e; return(invisible())
    }
    # by index: an empty argument (x[, 1]) is the missing symbol and cannot be bound to a variable
    for (i in seq_along(e)[-1]) if (is.call(e[[i]])) walk(e[[i]])
  }
  for (e in ex) walk(e)
  if (is.null(hit)) stop("definition of ", name, " not found in ", file)
  eval(hit, envir)
  get(name, envir)
}

.write_stamp <- function(dir, id, status = "complete", date = format(as.Date(ANALYSIS_DATE))) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(list(cascade_run_id = id, analysis_date = date, kernel = "M14-fill",
                            form = "fixed", status = status),
                       file.path(dir, "cascade_run_stamp.json"), auto_unbox = TRUE)
}

test_that("cascade_stamp_check accepts only the latest completed (or live) run for this date [stamp]", {
  d <- file.path(tempdir(), "casc_stamp_chk"); unlink(d, recursive = TRUE)
  ad <- as.Date(ANALYSIS_DATE)
  expect_false(cascade_stamp_check("r1", d, ad, live_run_id = "")$ok)          # no stamp
  .write_stamp(d, "r1", "complete")
  expect_true(cascade_stamp_check(c("r1", "r1"), d, ad, live_run_id = "")$ok)
  expect_false(cascade_stamp_check("r0", d, ad, live_run_id = "")$ok)          # older run's table
  expect_false(cascade_stamp_check(NULL, d, ad, live_run_id = "")$ok)          # unstamped table
  expect_false(cascade_stamp_check(NA_character_, d, ad, live_run_id = "")$ok)
  expect_false(cascade_stamp_check(c("r1", "r0"), d, ad, live_run_id = "")$ok) # mixed
  expect_false(cascade_stamp_check("r1", d, ad + 7, live_run_id = "")$ok)      # another date
  .write_stamp(d, "r1", "failed: calibration - no settled origin")
  r <- cascade_stamp_check("r1", d, ad, live_run_id = "r1")
  expect_false(r$ok); expect_match(r$why, "did not complete")
  .write_stamp(d, "r1", "running")
  expect_false(cascade_stamp_check("r1", d, ad, live_run_id = "")$ok)          # abandoned run
  expect_false(cascade_stamp_check("r1", d, ad, live_run_id = "r2")$ok)        # someone else's
  expect_true(cascade_stamp_check("r1", d, ad, live_run_id = "r1")$ok)         # the live run
})

test_that("cascade_live_run_id reads the environment, then run_cascade.R's global [stamp]", {
  old <- Sys.getenv("CASCADE_RUN_ID", NA)
  on.exit(if (is.na(old)) Sys.unsetenv("CASCADE_RUN_ID") else Sys.setenv(CASCADE_RUN_ID = old))
  Sys.setenv(CASCADE_RUN_ID = "env-id")
  expect_identical(cascade_live_run_id(), "env-id")
  Sys.unsetenv("CASCADE_RUN_ID")
  had <- exists("CASCADE_RUN_ID", envir = globalenv(), inherits = FALSE)
  if (!had) {
    assign("CASCADE_RUN_ID", "glob-id", envir = globalenv())
    on.exit(rm("CASCADE_RUN_ID", envir = globalenv()), add = TRUE)
    expect_identical(cascade_live_run_id(), "glob-id")
  }
})

test_that("cascade_calibration_note labels an uncalibrated or unrecorded run, in ASCII [stamp]", {
  expect_null(cascade_calibration_note(c(TRUE, TRUE), "calibrated"))
  n <- cascade_calibration_note(FALSE, "not calibrated (disabled: CASCADE_CALIBRATE=0)")
  expect_match(n, "^UNCALIBRATED cascade run \\(not calibrated \\(disabled: CASCADE_CALIBRATE=0\\)\\)")
  expect_match(cascade_calibration_note(NULL), "not recorded")
  expect_match(cascade_calibration_note(c(TRUE, FALSE)), "not recorded")
  expect_true(all(utf8ToInt(n) < 128L))
})

test_that("cascade_write_reach_csv stamps the run id, date and calibration flag", {
  f <- tempfile(fileext = ".csv")
  rs <- tibble::tibble(health_zone = c("A", "B"), province = "P", horizon = 13L,
                       was_active_before = FALSE, p_case_invasion = c(2.5e-5, 0.123456789),
                       p_lo = c(0, 0.1), p_hi = c(1e-4, 0.15), p_ci_at_mc_floor = c(TRUE, FALSE))
  out <- suppressMessages(cascade_write_reach_csv(rs, f, run_id = "rid", calibrated = FALSE,
                                                  calibration_status = "not calibrated"))
  back <- readr::read_csv(f, show_col_types = FALSE)
  expect_true(all(c("cascade_run_id", "analysis_date", "calibrated", "calibration_status",
                    "p_case_lo", "p_case_hi") %in% names(back)))
  expect_identical(unique(back$cascade_run_id), "rid")
  expect_identical(unique(as.character(back$analysis_date)), format(as.Date(ANALYSIS_DATE)))
  expect_identical(unique(back$calibrated), FALSE)
  expect_equal(back$p_case_invasion[back$health_zone == "A"], 2.5e-5)          # signif, not round(4)
  # a zero-row table still writes, with the columns
  out0 <- suppressMessages(cascade_write_reach_csv(rs[0, ], f, run_id = "rid", calibrated = TRUE))
  expect_true("cascade_run_id" %in% names(out0)); expect_equal(nrow(out0), 0L)
})

test_that("Panel B's label survives a fractional MI mean, and .panel() catches render-time errors [label]", {
  lab <- .def_from("40_cascade_next_dominoes.R", "newly_above_label")
  expect_identical(lab(c(3, 0, 1.4, NA)), c("+3 newly >20%", "+0 newly >20%", "+1.4 newly >20%", ""))
  # the old label is an error only at BUILD time, not when the plot object is constructed
  d <- data.frame(x = 1, y = "a", n = 1.4)
  p_bad <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("+%d newly >20%%", n)))
  expect_s3_class(p_bad, "ggplot")
  # .panel() as 40 defines it, with its failure log redirected into a local environment
  e <- new.env(parent = globalenv()); e$.panel_fail <- character(0)
  pan <- .def_from("40_cascade_next_dominoes.R", ".panel", envir = e)
  res <- suppressMessages(pan(p_bad, "cityrank"))
  expect_length(e$.panel_fail, 1L)
  expect_match(e$.panel_fail, "^cityrank \\(")
  p_ok <- ggplot2::ggplot(d, ggplot2::aes(x, y)) + ggplot2::geom_text(ggplot2::aes(label = lab(n)))
  expect_identical(suppressMessages(pan(p_ok, "cityrank")), p_ok)
  expect_length(e$.panel_fail, 1L)
})

test_that("Panel C's palette is built over the provinces it draws, so no two share a hue [palette]", {
  src <- paste(deparse(.def_from("40_cascade_next_dominoes.R", "build_watchlists")), collapse = "\n")
  expect_true(grepl("PCOL <- prov_palette(d$province)", src, fixed = TRUE))
  expect_true(grepl("anyDuplicated(unname(PCOL))", src, fixed = TRUE))
  # the defect: ~26 provinces through an 8-hue spare pool recycle hues
  many <- c(names(PROV_COL)[1:3], sprintf("Prov%02d", 1:20))
  expect_true(anyDuplicated(unname(prov_palette(many))) > 0)
  few <- c("Kinshasa", "Kwilu", "Haut-Katanga", "Kasaï", "Kasaï-Central", "Kasaï-Oriental")
  expect_equal(anyDuplicated(unname(prov_palette(few))), 0L)
  # 37 uses the same scale (no private OKABE-spare loop)
  s37 <- paste(readLines(file.path(.ST, "37_cascade_figure4.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("prov_cols   <- prov_palette(d$province)", s37, fixed = TRUE))
  expect_false(grepl("setdiff(OKABE, unname(PROV_COL))", s37, fixed = TRUE))
})

test_that("the seed-zone provenance is written from the RESOLVED hubs only [hubs]", {
  e <- new.env(parent = globalenv())
  wsel <- .def_from("38_urban_scenarios.R", "urban_write_hub_selection", envir = e)
  wtab <- .def_from("38_urban_scenarios.R", "urban_write_tables", envir = e)
  res  <- .def_from("38_urban_scenarios.R", "resolve_urban_hubs", envir = e)
  # no default: an unresolved URBAN_HUBS can no longer be written by omission
  expect_true(is.symbol(formals(wsel)$hubs) && !nzchar(as.character(formals(wsel)$hubs)))
  expect_false(grepl("urban_write_hub_selection", paste(deparse(body(wtab)), collapse = " "),
                     fixed = TRUE))
  cand <- tibble::tibble(city = c("C1", "C1", "C2"), Nom = c("z1", "z2", "z3"),
                         PROVINCE = "P", selected = TRUE)
  hubs <- structure(list(C1 = c("z1", "z2"), C2 = "z3"), candidates = cand, derived = TRUE)
  hubs_u <- suppressMessages(res(hubs, zones_all = c("z1", "z3")))   # z2 off the spine
  d <- file.path(tempdir(), "casc_hubsel"); dir.create(d, showWarnings = FALSE)
  suppressMessages(wsel(d, hubs_u))
  s <- readr::read_csv(file.path(d, "urban_hub_selection.csv"), show_col_types = FALSE)
  expect_true("seeded" %in% names(s))
  expect_identical(s$Nom[s$seeded], c("z1", "z3"))
})

test_that("39's capture curve is the per-origin mean recall, with English month labels [39]", {
  e <- new.env(parent = globalenv())
  suppressMessages(sys.source(file.path(.ST, "39_cascade_eval_figure.R"), envir = e))
  # origin 1: 1 invasion, caught at k = 1; origin 2: 3 invasions, 1 caught at k = 1
  det <- tibble::tibble(cutoff = rep(as.Date(c("2026-06-01", "2026-06-08")), each = 4),
                        health_zone = rep(c("a", "b", "c", "d"), 2),
                        p_invasion = c(.9, .1, .05, .01, .9, .5, .4, .3),
                        y = c(1L, 0L, 0L, 0L, 1L, 0L, 1L, 1L))
  cc <- e$.casc_capture_curve(det, ks = 1:2)
  expect_equal(cc$recall[cc$k == 1], mean(c(1, 1/3)))
  expect_equal(cc$recall_pooled[cc$k == 1], 2 / 4)
  old <- Sys.getlocale("LC_TIME"); on.exit(Sys.setlocale("LC_TIME", old), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_TIME", "fr_FR.UTF-8"))
  lab <- e$.casc_fold_labels(c("2026-06-01", "2026-05-25"))
  expect_true(all(grepl("(May|Jun)$", lab)))
  expect_identical(names(lab), c("2026-05-25", "2026-06-01"))
})

test_that("42 draws the sweep only when its sidecar matches the stamp and this configuration [42]", {
  root <- file.path(tempdir(), "casc42"); unlink(root, recursive = TRUE)
  diag <- file.path(root, "cascade", "diagnostics"); dir.create(diag, recursive = TRUE)
  ad <- format(as.Date(ANALYSIS_DATE))
  ks <- c(0.1, 0.3, 1)
  readr::write_csv(tibble::tibble(k_indiv = ks, spearman = c(.9, 1, .95), kendall = c(.8, 1, .9),
                                  top15_jaccard = c(.7, 1, .8), mean_reach_h4 = .01,
                                  mean_reach_h8 = .02, mean_reach_h13 = .03, mean_estab_h13 = .01),
                   file.path(diag, "cascade_kindiv_sweep.csv"))
  readr::write_csv(tidyr::expand_grid(k_indiv = ks, health_zone = c("A", "B")) |>
                     dplyr::mutate(p_invasion = 0.1 * k_indiv, province = "Ituri"),
                   file.path(diag, "cascade_kindiv_sweep_byzone.csv"))
  .meta <- function(...) {
    m <- utils::modifyList(list(cascade_run_id = "r42", analysis_date = ad, calibrated = TRUE,
                                kernel = "M14-fill", form = "fixed", k_base = 0.30,
                                scenario = "S1_status_quo",
                                seed = as.integer(get0("RANDOM_SEED", ifnotfound = 20260704L))),
                           list(...))
    jsonlite::write_json(m, file.path(diag, "cascade_kindiv_sweep_meta.json"), auto_unbox = TRUE)
  }
  .run <- function() suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), shQuote(file.path(.ST, "42_kindiv_sweep_figure.R")),
    env = c(paste0("CASCADE_OUT_ROOT=", root), paste0("ANALYSIS_DATE=", ad),
            "CASCADE_KERNEL=M14-fill", "CASCADE_FORM=fixed", "CASCADE_RUN_ID=",
            paste0("CASCADE_ST_DIR=", .ST)),
    stdout = TRUE, stderr = TRUE))
  fig <- file.path(root, "key_outputs", "Figure_kindiv_sweep.pdf")
  .write_stamp(diag, "r42", "complete"); .meta(form = "rw1")
  out <- .run()
  expect_true(is.null(attr(out, "status")) || attr(out, "status") == 0L)
  expect_true(any(grepl("Skipping", out)) && any(grepl("form \\(sweep rw1", out)))
  expect_false(file.exists(fig))
  .meta(); .write_stamp(diag, "r42", "failed: calibration - x")
  out <- .run()
  expect_true(any(grepl("did not complete", out))); expect_false(file.exists(fig))
  .write_stamp(diag, "r42", "complete")
  out <- .run()
  expect_true(any(grepl("sweep provenance OK", out)), info = paste(tail(out, 5), collapse = "\n"))
  expect_true(file.exists(fig))
})
