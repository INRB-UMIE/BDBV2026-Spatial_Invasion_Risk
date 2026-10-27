# =============================================================================
# 48_onset_mi_pool.R — pool the onset-imputation arms written by 47
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Reads outputs/onset_mi/arm_*.rds + draws_*.rds (+ each arm's cascade sandbox) and writes
#   (a) a MIRROR of outputs/ holding MI-POOLED artifacts in the EXACT schemas the figure
#       scripts read, so make_manuscript_figures.R and 40_cascade_next_dominoes.R draw the
#       pooled Figures 2/3/4 with MANU_IN_ROOT pointed here and no other change;
#   (b) the variability tables behind the MI report and 49's diagnostic figures.
# Design and the reasoning behind every pooling rule: PLAN_ONSET_MI.md.
#
# TWO POOLING ESTIMATORS, AND THE CHOICE IS NOT COSMETIC.
#
#   DRAW MIXTURE (Figure 3). The MI posterior is
#       p(theta | Y_obs) = INT p(theta | Y_obs, Y_mis) p(Y_mis | Y_obs) dY_mis,
#   and with posterior draws per completion that mixture is EXACT: stack the draws. This
#   runs through the pipeline's own .mix_and_summarise_draws(), which resamples JOINT draws
#   so cross-zone correlation — and therefore the rank credible interval — survives pooling.
#   Rubin's rules here would force a normal approximation onto probabilities that live near
#   zero with hard right skew, and would produce intervals crossing zero.
#
#   RUBIN (Figures 2 and 4, and every skill metric). Where the artifact is already a summary
#   carrying a within-arm variance, Q_bar = mean(Q_m), B = var(Q_m), T = U_bar + (1+1/M)B.
#   lambda = (1+1/M)B/T, the share of the total variance that is between-arm, is reported
#   throughout beside fmi (lambda with Rubin's small-sample df correction; onset_mi_pooling.R).
#   It slightly OVERSTATES the share due to not knowing the onset dates. The arms vary only the
#   onset seed (ONSET_IMPUTE_SEED) and hold the model seed (RANDOM_SEED) fixed, and section 1
#   refuses arms that do otherwise; but each arm's Stan / EpiNow2 / epinowcast fits and
#   bootstraps run on a DIFFERENT completed data set, so their finite-draw Monte Carlo error
#   still enters B (see 47's header).
#
# WHAT RUBIN DOES **NOT** LICENSE HERE. auc_pr_skill = auc_pr / base_rate, and base_rate is
# computed from the IMPUTED labels. The spread across arms is therefore uncertainty about
# WHAT THE TRUTH WAS, not sampling error about a fixed quantity. The decomposition is still
# well defined and is reported, but base_rate, n_invasions and n_atrisk are carried beside
# every skill number so a change in skill can never be read as a change in discrimination
# without checking the denominator.
#
# Env:
#   ONSET_MI_OUT   arm/pooled root (default outputs/onset_mi)
#   ONSET_MI_MIN_ARMS  refuse to pool fewer than this many arms (default 2)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(jsonlite)
})

ROOT   <- "/Users/mills/Documents/BDBV2026-Analysis"
ST_DIR <- file.path(ROOT, "spatiotemporal")
OUT    <- file.path(ST_DIR, "outputs")
MI_OUT <- Sys.getenv("ONSET_MI_OUT", file.path(OUT, "onset_mi"))
POOL   <- file.path(MI_OUT, "pooled")
MIN_ARMS <- suppressWarnings(as.integer(Sys.getenv("ONSET_MI_MIN_ARMS", "2")))

# The pooled mirror must NOT be produced with a stray RANDOM_SEED / stratum in the
# environment: 21's .mix_and_summarise_draws() seeds its resample from RANDOM_SEED, so does
# the per-draw Gamma delta in .scale_invasion_draws(), and the risk-score layer re-reads the
# stratum. RANDOM_SEED is resolved EXACTLY as 47 resolves it for the arms: the published run's
# (run_info.json onset_imputation$random_seed) where recorded, otherwise 00_config.R's default
# (the variable is unset so the config's own default applies — no literal is restated here).
# Environment changes only when this file RUNS as the pool (Rscript, or FORCE_RUN=1 -- the rule of
# the script_run_or_forced() guard below, which cannot be called before 00_config.R is sourced):
# a source() must not change the caller's environment.
.run_pool <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  (length(f) == 1L && identical(basename(f), "48_onset_mi_pool.R")) ||
    tolower(trimws(Sys.getenv("FORCE_RUN", ""))) %in% c("1", "true", "t", "yes", "y")
})
PUB_MODEL_SEED <- local({
  .r <- tryCatch(jsonlite::fromJSON(file.path(OUT, "key_outputs", "run_info.json"),
                                    simplifyVector = TRUE)$onset_imputation$random_seed,
                 error = function(e) NULL)
  .v <- if (length(.r) == 1L) suppressWarnings(as.integer(.r)) else NA_integer_
  if (.run_pool) { if (is.na(.v)) Sys.unsetenv("RANDOM_SEED") else Sys.setenv(RANDOM_SEED = as.character(.v)) }
  .v
})
if (.run_pool) Sys.setenv(ONSET_IMPUTE_STRATUM = "confirmed")
source(file.path(ST_DIR, "00_config.R"))
# The model seed every arm must have held: the published one, else the config default this
# session just resolved.
EXPECTED_MODEL_SEED <- if (!is.na(PUB_MODEL_SEED)) PUB_MODEL_SEED else
  suppressWarnings(as.integer(get0("RANDOM_SEED", ifnotfound = NA_integer_)))
if (.run_pool && (length(EXPECTED_MODEL_SEED) != 1L || is.na(EXPECTED_MODEL_SEED) ||
    !identical(as.integer(get0("RANDOM_SEED", ifnotfound = NA_integer_)), EXPECTED_MODEL_SEED)))
  stop("[mi-pool] could not resolve the model seed (RANDOM_SEED) this session must hold.",
       call. = FALSE)
stop_if_data_definition_sensitivity("48_onset_mi_pool.R")   # hard-coded published paths
# ONE module list, from 00_config.R. Nine files used to hand-maintain their own and all nine
# differed — see load_pipeline_modules().
# `extra`: the Rubin estimators live in their own file so tests/test_onset_mi.R can exercise
# them without loading a sweep, and forecast_scale.R supplies fs_lfo_col() — WHICH probability
# column a figure draws — which section 6 must publish Figure 2B on.
load_pipeline_modules("model", extra = c("onset_mi_pooling.R", "forecast_scale.R"),
                      dir = ST_DIR)


# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("48_onset_mi_pool.R")) {
  message("[mi-pool] sourced, not run: nothing is written. Use `Rscript 48_onset_mi_pool.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

for (d in c(POOL, file.path(POOL, "forecasts"), file.path(POOL, "diagnostics"),
            file.path(POOL, "key_outputs"), file.path(POOL, "cascade", "tables"),
            file.path(MI_OUT, "tables")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD + THE VALIDATION GAUNTLET
# =============================================================================
arm_files <- sort(list.files(MI_OUT, "^arm_[0-9]+\\.rds$", full.names = TRUE))
if (length(arm_files) < MIN_ARMS)
  stop(sprintf("[mi-pool] found %d arm file(s) in %s but at least %d are required.",
               length(arm_files), MI_OUT, MIN_ARMS), call. = FALSE)
arms <- lapply(arm_files, readRDS)
M <- length(arms)
message(sprintf("[mi-pool] %d arms: %s", M,
                paste(vapply(arms, function(a) a$tag, character(1)), collapse = ", ")))

.one <- function(f, what, fmt = as.character) {
  v <- unique(vapply(arms, function(a) fmt(f(a))[1], character(1)))
  if (length(v) != 1L)
    stop(sprintf(paste0("[mi-pool] the arms disagree on %s: %s. They are not completions of ",
                        "one study and must not be pooled."),
                 what, paste(v, collapse = " / ")),
         call. = FALSE)
  v
}
POOL_METHOD  <- .one(function(a) a$method,        "the featured model")
POOL_STRATUM <- .one(function(a) a$stratum,       "the onset-imputation stratum")
POOL_ADATE   <- .one(function(a) a$analysis_date, "the analysis date")
.one(function(a) a$training_cutoff, "the training cutoff")
POOL_NCONF   <- .one(function(a) a$n_confirmed,   "the confirmed-case count")
.one(function(a) a$n_draws,         "the retained draw count")

# ---- PROVENANCE: the arms must be completions of THE PUBLISHED RUN ----------------------------
# The pooled mirror is published UNDER THE PUBLISHED RUN'S LABELS: section 9 copies today's
# model_selection.json and run_info.json into it, and Figure 2B's legend dates come from
# run_info's training_window_end. Agreement AMONG the arms (.one() above) says nothing about
# agreement WITH those labels, and nothing checked it: the 20 arms on disk on 2026-10-15 were
# built on 2026-09-25 for analysis date 2026-09-07 and Bayes-M14-fill-geo, under the parametric
# onset draw, and would have been pooled and published as a completion of the 2026-10-09 run.
# Every check below refuses to pool rather than warn, and runs before anything is written.
.prov_fail <- function(...) stop("[mi-pool] ", ..., " These arms are not completions of the ",
                                 "published run: re-run 47 for it, or point ONSET_MI_OUT at the ",
                                 "sweep that matches it. Nothing was written.", call. = FALSE)
.ri_pub <- jsonlite::fromJSON(file.path(OUT, "key_outputs", "run_info.json"), simplifyVector = TRUE)
# The featured model through the suite's ONE reader (figure_style.R, sourced privately so its
# palette and theme objects stay out of this session).
read_featured_model <- local({
  .fs <- new.env()
  sys.source(file.path(ST_DIR, "figure_style.R"), envir = .fs)
  .fs$read_featured_model
})
PUB_FEATURED <- read_featured_model(file.path(OUT, "key_outputs", "model_selection.json"))
stopifnot(length(.ri_pub$analysis_date) == 1L, !is.na(.ri_pub$analysis_date))
if (!identical(as.character(as.Date(POOL_ADATE)), as.character(as.Date(.ri_pub$analysis_date))))
  .prov_fail(sprintf("The arms were built for analysis date %s but run_info.json describes %s.",
                     POOL_ADATE, .ri_pub$analysis_date))
if (!identical(POOL_METHOD, PUB_FEATURED))
  .prov_fail(sprintf("The arms fitted %s but model_selection.json features %s.",
                     POOL_METHOD, PUB_FEATURED))
# THE DATA DEFINITION. Arms written before these fields existed cannot be verified, which is a
# refusal, not a pass.
.missing_fields <- setdiff(c("data_definition", "linelist_snapshot", "onset_impute",
                             "random_seed", "onset_impute_seed"),
                           Reduce(intersect, lapply(arms, names)))
if (length(.missing_fields))
  .prov_fail("At least one arm predates the provenance field(s) ",
             paste(.missing_fields, collapse = ", "), ", so its data definition cannot be verified.")
local({
  .dd <- unique(lapply(arms, function(a) a$data_definition))
  if (length(.dd) != 1L || !identical(.dd[[1]], nondefault_data_definition()))
    .prov_fail("The arms were built under a different data definition (00_config.R switches) ",
               "from this session's published one.")
})
# TWO SEEDS: ONE MODEL SEED, M DISTINCT ONSET SEEDS. A completion is defined by its onset seed
# (ONSET_IMPUTE_SEED) alone; the model seed (RANDOM_SEED: Stan, EpiNow2, bootstraps, the delta
# draws) must be the SAME on every arm, or B would also carry deliberately different Monte-Carlo
# seeds, and it must be the published run's where run_info.json records it, else the
# 00_config.R default (EXPECTED_MODEL_SEED, resolved at the top of this file). A missing seed is
# a refusal: .one() collapses an all-NA field to a single NA, which used to skip the comparison.
POOL_RANDOM_SEED <- .one(function(a) a$random_seed, "the model seed (RANDOM_SEED)")
local({
  .prs <- suppressWarnings(as.integer(POOL_RANDOM_SEED))
  if (length(.prs) != 1L || is.na(.prs))
    .prov_fail("The arms do not record an integer model seed (RANDOM_SEED = ",
               POOL_RANDOM_SEED, ").")
  if (!identical(.prs, EXPECTED_MODEL_SEED))
    .prov_fail(sprintf("The arms held RANDOM_SEED = %d but the published run used %d (%s).",
                       .prs, EXPECTED_MODEL_SEED,
                       if (is.na(PUB_MODEL_SEED)) "00_config.R's default; run_info.json predates the field"
                       else "run_info.json onset_imputation$random_seed"))
  .os <- vapply(arms, function(a) as.integer(a$onset_impute_seed)[1], integer(1))
  if (anyNA(.os) || anyDuplicated(.os))
    .prov_fail("The arms' onset seeds (ONSET_IMPUTE_SEED) are missing or repeated (",
               paste(.os, collapse = ", "), "); each completion needs its own.")
  if (!identical(.os, vapply(arms, function(a) as.integer(a$seed), integer(1))))
    .prov_fail("An arm's recorded onset seed differs from the completion seed it was run as.")
  # ARM 0 IS THE PUBLISHED COMPLETION (run_onset_mi.sh; the report and FigureMI3 compare every
  # other arm against it). Where run_info.json records the published onset seed, arm 0 must
  # have drawn exactly that completion.
  .pub_os <- suppressWarnings(as.integer(.ri_pub$onset_imputation$onset_impute_seed))
  if (length(.pub_os) == 1L && !is.na(.pub_os)) {
    .a0 <- which(vapply(arms, function(a) identical(as.integer(a$arm), 0L), logical(1)))
    if (length(.a0) != 1L)
      .prov_fail(sprintf(paste0("run_info.json records the published onset seed (%d) but the ",
                                "sweep has no arm 0, the completion that must reproduce it."),
                         .pub_os))
    if (!identical(.os[.a0], .pub_os))
      .prov_fail(sprintf(paste0("Arm 0 drew onset seed %d but the published run used ",
                                "ONSET_IMPUTE_SEED = %d, so it is not the published completion."),
                         .os[.a0], .pub_os))
  }
})
POOL_SNAPSHOT <- .one(function(a) a$linelist_snapshot, "the line-list snapshot")
.pub_snapshot <- paste(.ri_pub$linelist_snapshot$folder, .ri_pub$linelist_snapshot$processed_at,
                       sep = "|")
if (!identical(POOL_SNAPSHOT, .pub_snapshot))
  .prov_fail(sprintf("The arms read line-list snapshot %s but run_info.json describes %s.",
                     POOL_SNAPSHOT, .pub_snapshot))
.impute_sig <- function(x) paste(x$mode, x$desc, x$donor_min, x$sample_date_source, sep = " | ")
POOL_IMPUTE <- .one(function(a) .impute_sig(a$onset_impute), "the onset-imputation mechanism")
# The mechanism THIS session's load_linelist() runs, and the published confirmed count, from the
# same data layer the vulnerability index is built from below (built here, once, so the check
# precedes every write). RANDOM_SEED is pinned above; the count does not depend on it.
.dat <- prep_all_data()
.cur_impute <- .impute_sig(list(mode = get0(".ONSET_IMPUTE_MODE", ifnotfound = NA_character_),
                                desc = get0(".ONSET_IMPUTE_DESC", ifnotfound = NA_character_),
                                donor_min = get0("ONSET_DONOR_MIN", ifnotfound = NA),
                                sample_date_source = get0("SAMPLE_DATE_SOURCE",
                                                          ifnotfound = NA_character_)))
if (!identical(POOL_IMPUTE, .cur_impute))
  .prov_fail("The arms imputed onsets as [", POOL_IMPUTE, "] but this session's load_linelist() ",
             "does [", .cur_impute, "].")
if (!identical(as.integer(POOL_NCONF), as.integer(sum(.dat$ll$confirmed %in% TRUE))))
  .prov_fail(sprintf("The arms counted %s confirmed records; this session's line list has %d.",
                     POOL_NCONF, sum(.dat$ll$confirmed %in% TRUE)))
message(sprintf("[mi-pool] provenance: analysis date %s, featured %s, snapshot %s, default data definition",
                POOL_ADATE, POOL_METHOD, POOL_SNAPSHOT))
# n_folds is pooled as a numeric MEAN along with every other numeric column, so a differing
# fold count would publish a fractional "n_folds" in the evaluation table and silently mean
# the arms were cross-validated over different rounds.
local({
  nf <- unique(unlist(lapply(arms, function(a) unique(a$eval$n_folds))))
  if (length(nf) != 1L)
    stop("[mi-pool] the arms were cross-validated over different fold counts (",
         paste(sort(nf), collapse = ", "), "). Their evaluations are not comparable.",
         call. = FALSE)
  message(sprintf("[mi-pool] all arms: %d CV folds", nf))
})

if (!identical(POOL_STRATUM, "confirmed"))
  stop("[mi-pool] the arms were run on stratum ", shQuote(POOL_STRATUM),
       "; this sweep is defined on 'confirmed'.", call. = FALSE)

# THE DELAY FIT IS IMPUTATION-INVARIANT — it is estimated on the COMPLETE onset+sample pairs
# among confirmed records, which are by construction the records that are NOT imputed. Any
# variation across arms therefore means an arm resolved a different stratum or a different
# params file, not a different completion.
.dlf <- do.call(rbind, lapply(arms, function(a)
  data.frame(family = a$delay_imputation$family, mean = round(a$delay_imputation$mean, 9),
             sd = round(a$delay_imputation$sd, 9), est = a$delay_imputation$estimator,
             n = a$delay_imputation$n_fit, stringsAsFactors = FALSE)))
if (nrow(unique(.dlf)) != 1L)
  stop("[mi-pool] the arms did not all read the SAME onset->sample delay (the nowcast's fallback ",
       "delay; the imputation draws donor intervals, not this fit):\n",
       paste(utils::capture.output(print(unique(.dlf))), collapse = "\n"),
       "\nThat fit is imputation-invariant, so this means an arm resolved a different ",
       "stratum or a different params file.", call. = FALSE)

# THE AS-OF FOLD TRUNCATION IS INVARIANT TOO, and for the same structural reason:
# .trunc_fit_asof_cohort() measures onset -> observation from the date_index of records whose
# recorded onset is USABLE (.trunc_asof_onset(), 02_epi_params.R) and leaves imputed records
# out — precisely the records this sweep does not impute. So every completion's folds are corrected by ONE shared distribution, which is what
# makes their evaluations comparable. Asserted rather than assumed: if this ever starts
# varying, the estimator has changed and the arms are no longer scored on a common footing.
local({
  .tr <- do.call(rbind, lapply(arms, function(a)
    data.frame(family = a$trunc_asof$family, mean = signif(a$trunc_asof$mean, 10),
               sd = signif(a$trunc_asof$sd, 10), est = a$trunc_asof$estimator,
               stringsAsFactors = FALSE)))
  if (nrow(unique(.tr)) != 1L)
    stop("[mi-pool] the arms were cross-validated under DIFFERENT as-of fold truncations:\n",
         paste(utils::capture.output(print(unique(.tr))), collapse = "\n"),
         "\nThat estimator is fitted only on records with a usable RECORDED onset, so it must not ",
         "move with the imputation. Their fold nowcasts are not comparable.", call. = FALSE)
  message(sprintf("[mi-pool] shared as-of fold truncation: %s (mean %.4f d) on all %d arms",
                  .tr$family[1], .tr$mean[1], M))
})

# EVERY ROUND'S TRAINING SOURCE MUST AGREE ACROSS ARMS. A snapshot round is nowcast with the
# extract truncation fitted on that date, which is refitted per arm (it replays the vintages under
# each arm's onset seed) and gated by its sanity check; a gate that flips in one arm switches that
# round to the re-censored estimator, and the arms would then not score the same round the same
# way.
local({
  .ts <- lapply(arms, function(a) {
    if (!"training_source" %in% names(a$lfo)) return(NULL)
    u <- unique(a$lfo[, c("fold_id", "training_source")])
    u[order(u$fold_id, u$training_source), , drop = FALSE]
  })
  if (all(vapply(.ts, is.null, logical(1)))) return(invisible(NULL))
  if (any(vapply(.ts, is.null, logical(1))) ||
      !all(vapply(.ts[-1], function(x) isTRUE(all.equal(x, .ts[[1]], check.attributes = FALSE)), logical(1))))
    stop("[mi-pool] the arms do not use the same training source (archived snapshot vs re-censored) ",
         "in every round; their back-tests are not comparable round by round.", call. = FALSE)
})

# THE FINGERPRINTS MUST BE PAIRWISE DISTINCT. Each is a hash of that arm's per-zone FIRST
# CONFIRMED ONSET — the quantity the imputation moves. A repeat means two "imputations" are
# the same completion, which does not merely weaken the result: it DEFLATES B toward zero
# and reports a narrower interval than a single draw would honestly support. This is the one
# failure mode of the whole design, so it is checked here regardless of the seed assertions
# each arm already passed.
.fp <- vapply(arms, function(a) a$fingerprint, character(1))
if (anyDuplicated(.fp)) {
  .d <- vapply(arms, function(a) a$tag, character(1))[duplicated(.fp) | duplicated(.fp, fromLast = TRUE)]
  stop("[mi-pool] arms ", paste(.d, collapse = ", "), " share an onset fingerprint: they are ",
       "the SAME completion. Pooling them would understate the between-imputation variance. ",
       "Check that ONSET_IMPUTE_SEED actually varied (47 asserts this) and that no cache served a ",
       "stale completion.", call. = FALSE)
}
.seeds <- vapply(arms, function(a) as.integer(a$seed), integer(1))
if (anyDuplicated(.seeds))
  stop("[mi-pool] duplicate arm seeds: ", paste(.seeds[duplicated(.seeds)], collapse = ", "),
       call. = FALSE)
message(sprintf("[mi-pool] %d distinct onset completions, onset seeds %d-%d, RANDOM_SEED %s on every arm",
                length(unique(.fp)), min(.seeds), max(.seeds), POOL_RANDOM_SEED))

# partial_cv means an arm was SCORED on fewer cells, not fitted on fewer folds. Pooling one
# into a mean silently mixes two supports.
.pc <- vapply(arms, function(a) any(a$eval$partial_cv %in% TRUE), logical(1))
if (any(.pc))
  stop("[mi-pool] arm(s) ", paste(vapply(arms, function(a) a$tag, character(1))[.pc], collapse = ", "),
       " are flagged partial_cv (scored on fewer cells than the full support). Re-run them ",
       "or drop them deliberately; they must not be averaged in silently.", call. = FALSE)

# THE COMPARATOR SET MUST BE THE SAME ON EVERY ARM. 00_config.R is a live file — it was
# edited mid-sweep once already — and an arm run before such an edit scores a different set
# from one run after it, which would silently give Figure 2A series averaged over different
# numbers of completions.
local({
  .nz <- function(x) if (is.null(x)) character(0) else as.character(x)
  .fb <- lapply(arms, function(a) sort(.nz(a$figure_baselines)))
  if (length(unique(.fb)) != 1L)
    stop("[mi-pool] the arms were run under DIFFERENT FIGURE_BASELINE_METHODS settings:\n  ",
         paste(unique(vapply(.fb, paste, character(1), collapse = ", ")), collapse = "\n  "),
         "\n00_config.R changed during the sweep. Re-run the arms that predate the change.",
         call. = FALSE)
  .dm <- lapply(arms, function(a) sort(.nz(a$dc_methods)))
  if (length(unique(.dm)) != 1L)
    stop("[mi-pool] the arms scored different detection-curve method sets.", call. = FALSE)
  message(sprintf("[mi-pool] comparator set on all %d arms: %s", M,
                  paste(.dm[[1]], collapse = ", ")))
})

HAS_CASCADE <- all(!is.na(vapply(arms, function(a) a$cascade_dir, character(1))))
message(sprintf("[mi-pool] cascade arms: %s", if (HAS_CASCADE) "all present" else "ABSENT — Figure 4 will not be pooled"))
# 40_cascade_next_dominoes.R draws Figure 4 panel A from urban_hub_selection.csv (written by
# 38_urban_scenarios.R) and, without it, silently falls back to a frozen hub list that need not
# match the zones the cascade seeded. Section 8 copies it into the mirror; checked here so a
# missing one refuses before anything is written rather than after half the mirror exists.
if (HAS_CASCADE) local({
  .hf <- file.path(vapply(arms, function(a) a$cascade_dir, character(1)),
                   "cascade", "tables", "urban_hub_selection.csv")
  if (!all(file.exists(.hf)))
    stop("[mi-pool] urban_hub_selection.csv is missing from cascade arm(s) ",
         paste(vapply(arms, function(a) a$tag, character(1))[!file.exists(.hf)], collapse = ", "),
         ". The pooled Figure 4 panel A would fall back to the frozen hub list. Re-run those ",
         "arms' cascades. Nothing was written.", call. = FALSE)
})

# ---- THE DRAWS FILES: one per arm, and each the SAME completion as its arm file ------------
# The Figure-3 mixture reads the draws files and every other product reads the arm files, so
# a draws file left over from another sweep (or a re-run arm whose draws were not rewritten)
# would put a different completion into Figure 3 under this arm's tag. Checked here, before
# the first write, on the fingerprint and both seeds.
draw_files <- sort(list.files(MI_OUT, "^draws_[0-9]+\\.rds$", full.names = TRUE))
if (length(draw_files) != M)
  stop(sprintf(paste0("[mi-pool] %d arm files but %d draws files; every arm must contribute ",
                      "its posterior draws to the mixture."),
               M, length(draw_files)), call. = FALSE)
dws <- lapply(draw_files, readRDS)
.arm_tags <- vapply(arms, function(a) as.character(a$tag)[1], character(1))
.dw_tags  <- vapply(dws, function(d) as.character(d$tag)[1], character(1))
if (anyNA(.dw_tags) || !identical(sort(.dw_tags), sort(.arm_tags)))
  stop("[mi-pool] the arm and draws files do not describe the same arms.", call. = FALSE)
local({
  .f1 <- function(v) if (is.null(v) || !length(v) || is.na(v[1])) "<missing>" else as.character(v[1])
  for (d in dws) {
    a <- arms[[match(d$tag, .arm_tags)]]
    for (f in c("fingerprint", "arm", "seed", "onset_impute_seed", "random_seed")) {
      x <- .f1(d[[f]]); y <- .f1(a[[f]])
      if (identical(x, "<missing>") || !identical(x, y))
        .prov_fail(sprintf(paste0("draws_%s.rds records %s = %s but arm_%s.rds records %s: the ",
                                  "draws are not this arm's completion."), d$tag, f, x, d$tag, y))
    }
  }
})
# THE PROBABILITY SCALE OF THE MIXTURE. Every arm's draws are put on that arm's own deployed
# scale before mixing (section 7). An arm whose deploy delta is NULL — 47 only WARNS when
# invasion_delta_table() fails — would contribute RAW draws to a mixture published as
# recalibrated. All arms recalibrated, or none (then the scale is labelled "raw-mi").
.draw_h <- sort(unique(unlist(lapply(dws, function(d)
  vapply(d$draws, function(el) as.numeric(el$horizon), numeric(1))))))
POOL_PROB_SCALE <- mi_deploy_scale(stats::setNames(lapply(dws, function(d) d$deploy_delta), .dw_tags),
                                   .draw_h)
message(sprintf("[mi-pool] Figure-3 mixture probability scale: %s (horizons %s)",
                POOL_PROB_SCALE, paste(.draw_h, collapse = ", ")))

# =============================================================================
# 2. RUBIN
# =============================================================================
# The estimators themselves are in onset_mi_pooling.R (sourced above):
#   mi_rubin()          Q_bar, B, U_bar, T, FMI for one estimand
#   mi_u_from_ci()      a published 90% bootstrap interval -> a within variance
#   mi_pool_by()        one Rubin result per key combination of a long per-arm frame
#   mi_zero_affected()  the "already invaded in THIS completion -> p = 0" rule (D5)
.rubin      <- mi_rubin
.u_from_ci  <- mi_u_from_ci
.pool_by    <- mi_pool_by

# =============================================================================
# 3. POOLED LFO FRAME  (published for row-level consumers)
# =============================================================================
# The row-wise stack of every arm, carrying `mi_arm` and `mi_weight`. The WEIGHT is the point:
# arms do NOT all score the same number of cells (the at-risk support moves with the imputed
# labels), so a consumer that averages rows would weight each completion by how many rows it
# happened to contribute. (Figure 2A's block does not read it: make_manuscript_figures.R scores
# each arm on that arm's shared support and recalibrated rounds and averages the arms equally.)
# A consumer averaging WITHIN (method, horizon) is served by a weight of 1 / (M * n_scored) with
# n_scored counted within (arm, method, horizon), and every completion sums to 1/M in every
# such slice. (Counting n_scored over an arm's whole frame equalised the arms only across the
# stack, not within the slice the consumer averages.) "Scored" is the consumer's own filter:
# a finite probability and outcome on a zone not already invaded.
lfo_pooled <- dplyr::bind_rows(lapply(arms, function(a) {
  d <- a$lfo
  d$mi_arm  <- a$arm
  d$mi_seed <- a$seed
  d
}))
lfo_pooled$mi_weight <- mi_slice_weights(
  scored = is.finite(lfo_pooled$p_invasion) & is.finite(lfo_pooled$is_new_invasion) &
           !(lfo_pooled$was_active_before %in% TRUE),
  arm = lfo_pooled$mi_arm, slice = list(lfo_pooled$method, lfo_pooled$horizon), M = M)
attr(lfo_pooled, "onset_mi") <- list(n_imputations = M, method = POOL_METHOD,
                                     stratum = POOL_STRATUM, seeds = .seeds)

# THE RECALIBRATED COLUMN MUST SURVIVE THE STACK, on every scored row that can carry one.
# fs_apply_lfo_scale() demotes the ENTIRE frame to RAW probabilities if any such row lacks
# `p_recal`, and fs_stamp() keys on the REQUESTED scale — so a single arm whose 16b attach
# failed would publish both pooled figures on raw probabilities under captions saying
# "recalibrated", with nothing anywhere saying so. That exact mismatch has happened in this
# tree before. Check it here, where the arm responsible can still be named.
local({
  if (!"p_recal" %in% names(lfo_pooled))
    stop("[mi-pool] the pooled LFO frame carries no `p_recal`: no arm attached the 16b ",
         "recalibration. Both pooled figures would silently fall back to RAW probabilities ",
         "while their captions claim otherwise. Was INVASION_RECALIBRATE on for the arms?",
         call. = FALSE)
  scored <- is.finite(lfo_pooled$p_invasion) &
            !(as.logical(lfo_pooled$was_active_before) %in% TRUE)
  # A rank-only baseline has no p_recal BY CONSTRUCTION, not by failure: 16b deliberately
  # never fits one a delta, and they flag themselves prob_calibrated = FALSE.
  if ("prob_calibrated" %in% names(lfo_pooled))
    scored <- scored & !(lfo_pooled$prob_calibrated %in% FALSE)
  bad <- scored & !is.finite(lfo_pooled$p_recal)
  if (any(bad)) {
    by_arm <- table(lfo_pooled$mi_arm[bad])
    stop(sprintf(paste0("[mi-pool] %d of %d scored rows have no recalibrated probability, in ",
                        "arm(s) %s. The pooled figures would be drawn on RAW probabilities ",
                        "under a 'recalibrated' caption. Re-run those arms."),
                 sum(bad), sum(scored), paste(names(by_arm), collapse = ", ")), call. = FALSE)
  }
  message(sprintf("[mi-pool] recalibrated column intact on all %d scored rows", sum(scored)))
})
# THE MIRROR IS WRITTEN ALL-OR-NOTHING. A marker goes down before the first file and is removed
# only after the last, so a pool that refuses part-way (several checks below fire after this
# point) leaves a mirror 49 will not draw from. Without a cascade, the cascade part of an
# EARLIER sweep's mirror is removed: it would otherwise be drawn as this sweep's.
.pool_marker <- file.path(POOL, ".pool_incomplete")
dir.create(POOL, recursive = TRUE, showWarnings = FALSE)
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), .pool_marker)
if (!HAS_CASCADE) unlink(file.path(POOL, "cascade"), recursive = TRUE)
saveRDS(lfo_pooled, file.path(POOL, "forecasts", "lfo_results.rds"))
message(sprintf("[mi-pool] pooled LFO frame: %d rows (%d arms x ~%d)",
                nrow(lfo_pooled), M, round(nrow(lfo_pooled) / M)))

# =============================================================================
# 4. POOLED EVALUATION TABLE
# =============================================================================
ev_long <- dplyr::bind_rows(lapply(arms, function(a)
  a$eval %>% dplyr::mutate(mi_arm = a$arm)))
.ev_keys <- c("method", "horizon")
# Every metric with a published within-arm interval gets a real U; everything else pools
# with U = 0 and FMI = 1, which is flagged in the variability table's `u_known` column.
# DERIVED FROM THE COLUMNS, not listed: any stem carrying both <stem>_lo and <stem>_hi. A
# fixed c(auc_pr, log_score) went stale when evaluate_invasion() added auc_pr_skill_date_lo/hi,
# whose bounds were then averaged as metrics and whose estimand got U = 0.
# THE INTERVAL LEVEL. Every such interval in evaluate_invasion() (16_invasion_eval.R) is the
# 5%/95% quantile pair of its zone-cluster bootstrap — a 90% interval — so U is recovered at
# level 0.90 (tests/test_onset_mi.R pins that coupling against 16's source).
EV_CI_LEVEL <- 0.90
.ev_u <- mi_ci_stems(names(ev_long))
.ev_u <- stats::setNames(.ev_u, .ev_u)
message(sprintf("[mi-pool] evaluation metrics with a within-arm interval: %s",
                if (length(.ev_u)) paste(.ev_u, collapse = ", ") else "none"))
# THEIR _lo/_hi COLUMNS ARE NOT ESTIMANDS. They are each arm's within-arm interval, consumed
# below as U. Pooling them as metrics of their own published the MEAN of the per-arm bounds as
# the pooled interval — a within-arm band that leaves out B entirely, under the label of the MI
# result. The pooled interval for each of these metrics is Rubin's (see ev_pooled below).
.ev_ci_cols <- as.vector(outer(unname(.ev_u), c("_lo", "_hi"), paste0))
.ev_num  <- setdiff(names(ev_long)[vapply(ev_long, is.numeric, logical(1))],
                    c(.ev_keys, "mi_arm", .ev_ci_cols))
ev_var <- dplyr::bind_rows(lapply(.ev_num, function(cn) {
  .has_ci <- cn %in% names(.ev_u) && all(paste0(cn, c("_lo", "_hi")) %in% names(ev_long))
  # The estimate and its interval are carried through the SAME select and renamed together,
  # so the within variance cannot be paired with the wrong row. Assigning `.u` from a second
  # read of ev_long would have relied on select() preserving row order — true today, and
  # exactly the kind of assumption that silently breaks an estimate-to-variance pairing.
  d <- ev_long %>%
    dplyr::select(dplyr::all_of(c(.ev_keys, "mi_arm", cn)),
                  dplyr::any_of(if (.has_ci) paste0(cn, c("_lo", "_hi")) else character(0)))
  names(d)[names(d) == cn] <- ".q"
  if (.has_ci) {
    d$.u <- .u_from_ci(d[[paste0(cn, "_lo")]], d[[paste0(cn, "_hi")]], level = EV_CI_LEVEL)
    .pool_by(d, .ev_keys, ".q", ".u") %>% dplyr::mutate(metric = cn, .before = 1)
  } else {
    .pool_by(d, .ev_keys, ".q") %>% dplyr::mutate(metric = cn, .before = 1)
  }
}))
readr::write_csv(ev_var %>% dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 8))),
                 file.path(MI_OUT, "tables", "onset_mi_evaluation_variability.csv"))

# The pooled invasion_evaluation.csv itself: same schema, every numeric column replaced by
# its MI mean, so Figure 2's annotation block reads pooled numbers through its existing
# `gv()` lookups and no call site in that script has to know MI happened.
ev_pooled <- ev_var %>%
  dplyr::select(dplyr::all_of(c("metric", .ev_keys)), q) %>%
  tidyr::pivot_wider(names_from = metric, values_from = q)
# The pooled INTERVAL columns, in the schema's own names: Rubin's 90% t interval (U from the
# per-arm bootstrap bands, plus (1+1/M)B), so a consumer reading auc_pr_lo/auc_pr_hi reads the
# MI total. AUC-PR is a proportion and is clamped to [0, 1], as the detection-curve bands are.
for (.m in intersect(names(.ev_u), unique(ev_var$metric))) {
  .ci <- ev_var %>% dplyr::filter(metric == .m) %>%
    dplyr::select(dplyr::all_of(.ev_keys), lo, hi)
  if (identical(.m, "auc_pr")) .ci <- .ci %>% dplyr::mutate(lo = pmax(0, lo), hi = pmin(1, hi))
  names(.ci)[names(.ci) == "lo"] <- paste0(.m, "_lo")
  names(.ci)[names(.ci) == "hi"] <- paste0(.m, "_hi")
  ev_pooled <- ev_pooled %>% dplyr::left_join(.ci, by = .ev_keys)
}
.ev_chr <- setdiff(names(arms[[1]]$eval), c(.ev_num, .ev_keys, .ev_ci_cols))
if (length(.ev_chr))
  ev_pooled <- ev_pooled %>%
    dplyr::left_join(arms[[1]]$eval %>% dplyr::select(dplyr::all_of(c(.ev_keys, .ev_chr))),
                     by = .ev_keys)
ev_pooled$n_imputations <- M
readr::write_csv(ev_pooled, file.path(POOL, "diagnostics", "invasion_evaluation.csv"))
message(sprintf("[mi-pool] pooled evaluation: %d method x horizon rows, %d metrics",
                nrow(ev_pooled), length(.ev_num)))

# =============================================================================
# 5. POOLED DETECTION CURVES  (Figure 2A)
# =============================================================================
dc_long <- dplyr::bind_rows(lapply(arms, function(a)
  a$detection %>% dplyr::mutate(mi_arm = a$arm)))
# recall_lo/hi are the 5%/95% quantiles of compute_detection_curve()'s fold bootstrap
# (20_forecast_detail.R::.detc_boot_ci): a 90% interval.
dc_long$.u <- .u_from_ci(dc_long$recall_lo, dc_long$recall_hi, level = 0.90)
dc_keys <- c("method", "horizon", "k")
dc_r <- .pool_by(dc_long, dc_keys, "recall", ".u")
if (any(dc_r$n_arms != M))
  stop("[mi-pool] the detection curves do not cover every arm at every (method, horizon, k); ",
       "Figure 2A would draw series averaged over different numbers of completions.",
       call. = FALSE)

# recall_random is a property of the FOLD DESIGN (k / n_atrisk), and n_atrisk moves with the
# imputed labels, so it is pooled too — as a mean over arms. Figure 2A asserts it is
# identical across METHODS within a row set; that assertion still holds here because every
# arm's curve carries one value per (horizon, k) shared by all its methods, and the mean of
# equal values is equal. The SPREAD across arms is reported in the variability table rather
# than drawn, because a band on the random-targeting reference would suggest the null is
# uncertain in a way the panel does not mean.
# NOTE the column ORDER here is load-bearing. summarise() evaluates its expressions in
# sequence and later ones SEE the earlier results, so computing `n_atrisk_sd` after creating
# `n_atrisk_mean` would take the sd of a single scalar — silently NA for every row, in the
# one column that reports how much the at-risk support moved between completions. The sd is
# therefore taken from the raw column first, under a name that cannot shadow it.
# ONE ROW PER (ARM, horizon, k) FIRST. n_atrisk and recall_random are properties of the FOLD
# DESIGN, so every METHOD within an arm carries the same value; summarising dc_long directly
# would take the spread over M x n_methods rows in which each arm appears n_methods times.
# The mean is unaffected by that duplication but the SD is not: with M = 2 and 4 methods it
# came out a factor sqrt(4/7) = 0.76 too small, in the one column that reports how far the
# at-risk support moved between completions.
.dc_design <- dc_long %>%
  dplyr::distinct(mi_arm, horizon, k, n_atrisk_mean, recall_random)
if (nrow(.dc_design) != M * dplyr::n_distinct(dc_long$horizon) * dplyr::n_distinct(dc_long$k))
  warning("[mi-pool] n_atrisk / recall_random are not constant across methods within an arm; ",
          "the support-spread table below is taken over the distinct values that exist.",
          call. = FALSE)
dc_rand <- .dc_design %>% dplyr::group_by(horizon, k) %>%
  dplyr::summarise(n_atrisk_sd = stats::sd(n_atrisk_mean, na.rm = TRUE),
                   n_atrisk_min = min(n_atrisk_mean, na.rm = TRUE),
                   n_atrisk_max = max(n_atrisk_mean, na.rm = TRUE),
                   n_arms = dplyr::n(),
                   recall_random = mean(recall_random, na.rm = TRUE),
                   n_atrisk_mean = mean(n_atrisk_mean, na.rm = TRUE), .groups = "drop")
.dc_other <- dc_long %>% dplyr::group_by(dplyr::across(dplyr::all_of(dc_keys))) %>%
  dplyr::summarise(dplyr::across(dplyr::any_of(c("recall_pooled", "precision", "n_monitored",
                                                 "caught", "n_events", "n_folds")),
                                 ~ mean(.x, na.rm = TRUE)), .groups = "drop")
detection_pooled <- dc_r %>%
  dplyr::transmute(method, horizon, k,
                   recall = q,
                   # THE BAND IS RUBIN TOTAL, not the within-arm bootstrap. That is the whole
                   # point of the pooled panel: the published single-draw ribbon carries only
                   # U, and the MI ribbon carries U + (1+1/M)B. Clamped to [0,1] because a
                   # normal-theory interval on a proportion can leave the unit interval.
                   recall_lo = pmax(0, lo), recall_hi = pmin(1, hi),
                   recall_within_sd = sqrt(pmax(U, 0)),
                   recall_between_sd = sqrt(pmax(B, 0)),
                   # The within-component band, drawn beside the total on FigureMI4. It is
                   # published HERE, on the SAME t critical value the total uses, so the
                   # comparison shows the between-imputation variance rather than the gap
                   # between two reference distributions — and so the figure computes nothing.
                   recall_within_lo = pmax(0, q - mi_tcrit(df) * sqrt(pmax(U, 0))),
                   recall_within_hi = pmin(1, q + mi_tcrit(df) * sqrt(pmax(U, 0))),
                   recall_lambda = lambda, recall_fmi = fmi, recall_df = df,
                   n_imputations = n_arms) %>%
  dplyr::left_join(.dc_other, by = dc_keys) %>%
  dplyr::left_join(dc_rand %>% dplyr::select(horizon, k, recall_random, n_atrisk_mean),
                   by = c("horizon", "k"))
readr::write_csv(detection_pooled %>%
                   dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 12))),
                 file.path(POOL, "key_outputs", "detection_curves.csv"))
readr::write_csv(dc_rand, file.path(MI_OUT, "tables", "onset_mi_support_spread.csv"))
message(sprintf("[mi-pool] pooled detection curves: %d method(s), median lambda %.3f (FMI %.3f)",
                dplyr::n_distinct(detection_pooled$method),
                stats::median(detection_pooled$recall_lambda, na.rm = TRUE),
                stats::median(detection_pooled$recall_fmi, na.rm = TRUE)))

# =============================================================================
# 6. FIGURE 2B DATA — bar = MI-mean probability, fill = share of imputations per class
# =============================================================================
# THE ESTIMAND IS CONDITIONAL ON BEING IN THE AT-RISK SET, and that is a deliberate
# difference from the Figure-3 rule below. Figure 3 forecasts EVERY zone, so a zone already
# invaded in some completion genuinely has p = 0 there and belongs in the mixture at zero.
# Panel 2B instead RANKS the at-risk zones of a round and colours them by an outcome class
# that is only defined for an at-risk zone — a zone already invaded before the cutoff is not
# in that population at all in that arm, and has no class. Everything here therefore
# conditions on "at risk in this arm at this round", and n_arms_at_risk is published so the
# denominator is visible rather than assumed.
H2 <- 2L
# THE PROBABILITY SCALE. make_manuscript_figures.R calls fs_apply_lfo_scale() at load time,
# which REPLACES lfo$p_invasion with lfo$p_recal; every expression in build_fig2 that reads
# `p_invasion` is therefore reading the RECALIBRATED forecast. This table is consumed in
# place of those expressions, so it must be built on the same column — an earlier revision
# published the RAW p_invasion here and the pooled Panel 2B drew bars a factor of ~2.2 too
# long in the later rounds (raw 0.352 against recalibrated 0.165 for Jiba at the 04 Aug
# round), under a caption and beside a Panel 2A that both describe the recalibrated scale.
# The ORDERING is unaffected — the transform is monotone within a fold — so only the bar
# lengths and the x axis were wrong, which is exactly the kind of error a ranking-based
# sanity check does not catch.
# The pool publishes the DEPLOYED-scale risk table only (no raw twin), so a raw pass would draw
# Figure 2B raw beside a recalibrated Figure 3; refused rather than mixed.
if (fs_is_raw())
  stop("[mi-pool] FORECAST_SCALE=raw is not supported here: the pooled Figure-3 risk table is ",
       "published on the deployed scale only.", call. = FALSE)
.pcol <- unique(vapply(arms, function(a) fs_lfo_col(a$lfo), character(1)))
if (length(.pcol) != 1L)
  stop("[mi-pool] the arms resolve DIFFERENT probability scales (",
       paste(.pcol, collapse = ", "), "); Figure 2B would mix them in one panel.",
       call. = FALSE)
message(sprintf("[mi-pool] Figure 2B probability scale: '%s' (FORECAST_SCALE=%s)",
                .pcol, FORECAST_SCALE))
.f2b_arm <- dplyr::bind_rows(lapply(arms, function(a) {
  a$lfo %>%
    dplyr::filter(method == POOL_METHOD, horizon == H2, is.finite(.data[[.pcol]]),
                  !(was_active_before %in% TRUE)) %>%
    dplyr::transmute(mi_arm = a$arm, cutoff = as.Date(cutoff), health_zone,
                     p_draw = .data[[.pcol]], is_new_invasion,
                     affected_ever = health_zone %in% a$affected_now)
}))
# Rounds: the common fold sequence, then first / middle / last AMONG THE ROUNDS WHOSE
# RECALIBRATION FACTOR WAS ESTIMABLE — the rule make_manuscript_figures.R (build_fig2) applies to
# the single-draw panel, which the MI twin must mirror so the two pick the same rounds.
.folds_by_arm <- split(.f2b_arm$cutoff, .f2b_arm$mi_arm)
all_folds <- sort(Reduce(intersect, lapply(.folds_by_arm, function(x) sort(unique(x)))))
all_folds <- as.Date(all_folds, origin = "1970-01-01")
if (!length(all_folds)) stop("[mi-pool] the arms share no CV fold cutoffs.", call. = FALSE)
.n_per_arm <- vapply(.folds_by_arm, function(x) dplyr::n_distinct(x), integer(1))
if (any(.n_per_arm != length(all_folds)))
  warning(sprintf(paste0("[mi-pool] the arms do not all carry the same fold set (%s); Figure ",
                         "2B uses the %d cutoffs common to every arm."),
                  paste(range(.n_per_arm), collapse = "-"), length(all_folds)), call. = FALSE)
# ESTIMABLE ROUNDS, read the way the figure reads them: per cutoff, whether ANY of the featured
# model's h = 2 rows carries delta_estimable TRUE (the figure takes .est_by_cut over the same
# stacked frame in MI mode). The earliest rounds have no closed training fold to fit delta on,
# are drawn UNCORRECTED, and would put a second probability scale on the panel's shared x axis.
# On a raw-scale pass nothing is recalibrated, so every round is eligible; and when no round
# qualifies (or the frame predates the column) the figure falls back to every round, as here.
# round_index stays the position in the FULL common fold sequence (all_folds), not 1..3.
.elig_folds <- if (identical(.pcol, "p_recal")) {
  .est_cut <- lfo_pooled %>%
    dplyr::filter(method == POOL_METHOD, horizon == H2) %>%
    dplyr::mutate(cutoff = as.Date(cutoff)) %>%
    dplyr::group_by(cutoff) %>%
    dplyr::summarise(est = if ("delta_estimable" %in% names(lfo_pooled))
                             any(delta_estimable %in% TRUE) else TRUE, .groups = "drop")
  all_folds[.est_cut$est[match(all_folds, .est_cut$cutoff)] %in% TRUE]
} else all_folds
if (!length(.elig_folds)) .elig_folds <- all_folds
.sel <- unique(round(seq(1, length(.elig_folds), length.out = 3)))
fold_ord <- .elig_folds[.sel]
message(sprintf("[mi-pool] Figure 2B rounds chosen among %d of %d common rounds with an estimable delta",
                length(.elig_folds), length(all_folds)))

LV <- c("within", "later", "never")     # mapped to the figure's labels below
.f2b <- .f2b_arm %>%
  dplyr::filter(cutoff %in% fold_ord) %>%
  dplyr::mutate(cls = dplyr::case_when(is_new_invasion == 1 ~ LV[1],
                                       affected_ever       ~ LV[2],
                                       TRUE                ~ LV[3]))
.f2b_zone <- .f2b %>%
  dplyr::group_by(cutoff, health_zone) %>%
  dplyr::summarise(p_mi = mean(p_draw), n_arms_at_risk = dplyr::n(), .groups = "drop")
# A cell at risk in only a minority of completions is not a stable member of the panel's
# population; it is dropped from the top-20 candidate set and the count is reported.
.min_arms <- ceiling(M / 2)
.dropped <- sum(.f2b_zone$n_arms_at_risk < .min_arms)
.f2b_zone <- .f2b_zone %>% dplyr::filter(n_arms_at_risk >= .min_arms)
top20 <- .f2b_zone %>% dplyr::group_by(cutoff) %>%
  dplyr::slice_max(p_mi, n = 20, with_ties = FALSE) %>% dplyr::ungroup()
f2b_out <- .f2b %>%
  dplyr::semi_join(top20, by = c("cutoff", "health_zone")) %>%
  dplyr::count(cutoff, health_zone, cls, name = "n_cls") %>%
  dplyr::left_join(top20, by = c("cutoff", "health_zone")) %>%
  dplyr::mutate(share = n_cls / n_arms_at_risk,
                round_index = match(cutoff, all_folds),
                # THE SWEEP SIZE, not a per-cell count. The figure flags a cell whose bar is
                # an average over fewer than M completions, so it needs M itself; inferring
                # it as max(n_arms_at_risk) would understate it whenever EVERY drawn cell
                # happens to be partial, and then flag nothing.
                n_imputations = M) %>%
  tidyr::complete(tidyr::nesting(cutoff, health_zone, p_mi, n_arms_at_risk, n_imputations,
                                 round_index),
                  cls = LV, fill = list(n_cls = 0L, share = 0)) %>%
  dplyr::arrange(cutoff, dplyr::desc(p_mi), cls)
# The class LABELS must be the ones make_manuscript_figures.R defines, and it builds them
# from run_info's training_window_end. Read the SAME field rather than re-deriving the date,
# or the legend and the data would disagree about the cutoff day.
.ri <- jsonlite::fromJSON(file.path(OUT, "key_outputs", "run_info.json"))
DATA_END_LBL <- format(as.Date(.ri$training_window_end), "%d %b")
f2b_out$outcome <- c(within = "Invaded within round",
                     later  = sprintf("Invaded later (by %s)", DATA_END_LBL),
                     never  = sprintf("Not invaded (by %s)", DATA_END_LBL))[f2b_out$cls]
readr::write_csv(f2b_out %>% dplyr::select(cutoff, round_index, health_zone, p_mi, outcome,
                                           share, n_cls, n_arms_at_risk, n_imputations),
                 file.path(POOL, "key_outputs", "Figure2B_mi_data.csv"))
message(sprintf("[mi-pool] Figure2B data: rounds (forecast origins) %s; %d cells dropped (<%d arms at risk)",
                paste(format(lfo_origin(fold_ord), "%d %b"), collapse = ", "), .dropped, .min_arms))

# =============================================================================
# 7. POOLED CURRENT FORECAST — the DRAW MIXTURE  (Figure 3)
# =============================================================================
# `dws` (the draws files) was loaded and matched to the arm files — fingerprint, both seeds —
# in section 1, and POOL_PROB_SCALE resolved there from every arm's deployed delta.
.zref <- dws[[1]]$zones
for (d in dws) if (!identical(d$zones, .zref))
  stop("[mi-pool] arm ", d$tag, " carries a different zone vector. .mix_and_summarise_draws() ",
       "takes the zone names from the FIRST component only, so mixing these would silently ",
       "mislabel every column.", call. = FALSE)
.nd <- unique(vapply(dws, function(d) nrow(d$draws[[1]]$cum), integer(1)))
if (length(.nd) != 1L)
  stop("[mi-pool] the arms carry different draw counts (", paste(.nd, collapse = ", "),
       "). Equal draws per completion is a REQUIREMENT of the mixture: an arm with more ",
       "draws would weight its completion more heavily.", call. = FALSE)

# D5 — A ZONE ALREADY INVADED IN THIS COMPLETION CONTRIBUTES p = 0, NOT NA. The estimand is
# "P(this zone records its FIRST confirmed case in the next two weeks)", which is genuinely
# zero once the zone has been invaded. Zeroing the cumulative hazard is exactly that:
# p = 1 - exp(-0) = 0, and mu_forecast = colMeans(cum) picks up a zero. The rejected
# alternative — conditioning on at-risk-in-every-arm — silently changes the population and is
# selection on the outcome.
# Each arm is put on ITS OWN deployed scale first: delta is fitted per completion, so
# rescaling after mixing would apply one arm's factor to another's draws. Section 1 refused a
# sweep in which some arms carry a deployed delta and others do not (POOL_PROB_SCALE).
comp <- lapply(seq_len(M), function(i) {
  d <- dws[[i]]
  aff <- d$affected_now
  .scale_invasion_draws(mi_zero_affected(d$draws, aff), d$deploy_delta)
})
names(comp) <- vapply(dws, function(d) d$tag, character(1))
# Already invaded in EVERY completion -> genuinely not at risk -> masked to NA and excluded
# from the ranking, exactly as a single-draw run masks its affected set.
aff_all <- Reduce(intersect, lapply(dws, function(d) d$affected_now))
aff_any <- Reduce(union,     lapply(dws, function(d) d$affected_now))
message(sprintf("[mi-pool] already invaded in every arm: %d zones; in at least one: %d",
                length(aff_all), length(aff_any)))

preds_pooled <- .mix_and_summarise_draws(
  comp_draws = comp, weights = rep(1, M), affected_zones = aff_all,
  total_draws = M * .nd)
if (is.null(preds_pooled)) stop("[mi-pool] the draw mixture returned NULL.", call. = FALSE)
.wu <- attr(preds_pooled, "draws_per_component")
if (length(unique(.wu)) != 1L)
  stop("[mi-pool] the mixture drew unequal numbers from the components (",
       paste(range(.wu), collapse = "-"), "), so the completions are not equally weighted.",
       call. = FALSE)
preds_pooled$method <- POOL_METHOD
message(sprintf("[mi-pool] mixture: %d components x %d draws = %d pooled draws",
                M, .wu[1], sum(.wu)))

# n_arms_at_risk per zone — published so Figure 3B can flag a zone that is not at risk in
# every completion, and so a reader can see which rows the mixture averaged over a zero.
.n_atrisk_zone <- tibble::tibble(health_zone = .zref,
                                 n_arms_at_risk = vapply(.zref, function(z)
                                   sum(vapply(dws, function(d) !(z %in% d$affected_now), logical(1))),
                                   integer(1)))

# THE RISK LAYER IS BUILT BY THE PIPELINE'S OWN HELPERS, exactly as run_all.R builds it:
# compute_risk_scores() then add_risk_indices(). The relative-risk index and the priority
# score on the pooled map are therefore produced by the same estimator as the published one,
# and the only thing that differs between the two figures is the probabilities they were
# given. Hand-joining a V column instead would have been a second definition of the index.
#
# The covariate layer (population, CCVI, health sites, travel time) carries no line-list
# input, so the vulnerability index is imputation-invariant and is built ONCE here. That
# invariance is asserted against an arm's own index below rather than assumed.
province_map <- load_province_map()
# `.dat` is the data layer built in section 1 (the provenance checks read it).
osrm_mat_pool <- load_osrm()
vuln_index <- compute_vulnerability_index(.dat$covariates, names(.dat$pop),
                                          osrm_mat = osrm_mat_pool)
if (is.null(vuln_index))
  stop("[mi-pool] the vulnerability index could not be built; the pooled priority and ",
       "rr01 columns Figure 3 draws would be missing.", call. = FALSE)
local({
  a1 <- arms[[1]]$risk_scores
  if (is.null(a1) || !"V" %in% names(a1)) return(invisible())
  v1 <- a1 %>% dplyr::distinct(health_zone, V) %>% dplyr::filter(is.finite(V))
  v2 <- tibble::tibble(health_zone = vuln_index$nom, V2 = vuln_index$V)
  cmp <- dplyr::inner_join(v1, v2, by = "health_zone")
  d <- suppressWarnings(max(abs(cmp$V - cmp$V2), na.rm = TRUE))
  if (is.finite(d) && d > 1e-9)
    stop(sprintf(paste0("[mi-pool] the vulnerability index differs from arm %s's by %.2e. It ",
                        "is built from covariates that carry no line-list input, so it MUST ",
                        "be imputation-invariant; a difference means the covariate layer is ",
                        "not what the arms saw."), arms[[1]]$tag, d), call. = FALSE)
  message(sprintf("[mi-pool] vulnerability index matches arm %s (max |diff| %.2e over %d zones)",
                  arms[[1]]$tag, d, nrow(cmp)))
})
rs_pooled <- compute_risk_scores(preds_pooled, province_map)
rs_pooled <- add_risk_indices(rs_pooled, vuln_index)
rs_pooled <- rs_pooled %>% dplyr::left_join(.n_atrisk_zone, by = "health_zone")

# Same column set / rename / rounding as run_all.R's .write_risk_csv(), so the pooled file is
# the same file on the MI scale and the figure needs no new column knowledge.
.keep <- c("health_zone", "province", "horizon", "was_active_before", "prob_scale",
           "p_case_invasion", "p_lo", "p_hi", "p_median", "rank_med", "rank_lo", "rank_hi",
           "mu_forecast", "rr_nat", "rr_nat_rank", "rr01_nat", "rr_ituri", "rr_ituri_rank",
           "rr_nordkivu", "rr_nordkivu_rank", "rr_hautuele", "rr_hautuele_rank", "V",
           "priority", "priority_rank", "surveillance_gap", "healthcare_gap", "access_gap",
           "social_vulnerability", "healthcare_travel_min", "method", "n_arms_at_risk")
# The scale the mixture is ACTUALLY on (section 1): every arm recalibrated by its own deployed
# delta ("recalibrated-deployed-mi"), or none ("raw-mi"); a mixture of the two was refused.
rs_pooled$prob_scale <- POOL_PROB_SCALE
rs_out <- rs_pooled %>% dplyr::select(dplyr::any_of(.keep)) %>%
  dplyr::rename(dplyr::any_of(c(p_case_lo = "p_lo", p_case_hi = "p_hi"))) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 6))) %>%
  dplyr::arrange(horizon, dplyr::desc(dplyr::coalesce(p_case_invasion, -1)))
readr::write_csv(rs_out, file.path(POOL, "key_outputs", "bayes_risk_scores_all_zones.csv"))
message(sprintf("[mi-pool] pooled risk scores: %d rows, %d zones",
                nrow(rs_out), dplyr::n_distinct(rs_out$health_zone)))

# =============================================================================
# 8. POOLED CASCADE TABLES  (Figure 4)
# =============================================================================
if (HAS_CASCADE) {
  .ctab <- function(rel) lapply(arms, function(a) {
    f <- file.path(a$cascade_dir, "cascade", "tables", rel)
    if (!file.exists(f)) return(NULL)
    suppressWarnings(readr::read_csv(f, show_col_types = FALSE)) %>% dplyr::mutate(mi_arm = a$arm)
  })
  .need_all <- function(lst, what) {
    if (any(vapply(lst, is.null, logical(1))))
      stop("[mi-pool] ", what, " is missing from at least one cascade arm; Figure 4 would ",
           "pool over a different number of completions per panel.", call. = FALSE)
    dplyr::bind_rows(lst)
  }

  # -- reach scores: p_case_invasion, within SD published as p_sd ---------------
  rc <- .need_all(.ctab("cascade_reach_scores_all_zones.csv"), "cascade_reach_scores_all_zones.csv")
  rc$.u <- rc$p_sd^2
  rc_keys <- intersect(c("health_zone", "province", "horizon", "method"), names(rc))
  rc_p <- .pool_by(rc, rc_keys, "p_case_invasion", ".u")
  # `calibrated` (35 stamps each arm's reach table): TRUE only when EVERY arm's cascade was
  # calibrated, FALSE when any was not, NA when an arm does not record it -- so 40 labels an
  # uncalibrated pool rather than drawing it unannounced.
  .rc_has_cal <- "calibrated" %in% names(rc)
  rc_other <- rc %>% dplyr::group_by(dplyr::across(dplyr::all_of(rc_keys))) %>%
    dplyr::summarise(dplyr::across(dplyr::any_of(c("p_establishment", "first_passage_week",
                                                   "mu_forecast", "rr_nat", "rr01_nat", "V",
                                                   "priority", "rr_nat_rank", "priority_rank")),
                                   ~ mean(.x, na.rm = TRUE)),
                     was_active_before = any(was_active_before %in% TRUE),
                     calibrated = if (.rc_has_cal) {
                       .v <- as.logical(calibrated)
                       if (all(.v %in% TRUE)) TRUE else if (any(.v %in% FALSE)) FALSE else NA
                     } else NA,
                     .groups = "drop")
  reach_pooled <- rc_p %>%
    dplyr::transmute(dplyr::across(dplyr::all_of(rc_keys)),
                     p_case_invasion = q, p_case_lo = pmax(0, lo), p_case_hi = pmin(1, hi),
                     p_sd = se, p_within_sd = sqrt(pmax(U, 0)), p_between_sd = sqrt(pmax(B, 0)),
                     p_lambda = lambda, p_fmi = fmi, p_df = df, n_imputations = n_arms) %>%
    dplyr::left_join(rc_other, by = rc_keys)
  # RANKS ARE RECOMPUTED FROM THE POOLED SCORES, NOT AVERAGED. rc_other's mean of the per-arm
  # rr_nat_rank / priority_rank is not a rank (it is fractional, and two zones can share or
  # skip a place), and a mean rank is not the rank of the pooled score. Each is rebuilt by the
  # estimator that produced it per arm — min_rank(desc(.)) of mu_forecast within (method,
  # horizon) for rr_nat_rank (compute_risk_scores(), 15_workhorse.R) and of priority within
  # horizon for priority_rank (add_risk_indices(), 20_forecast_detail.R) — over the zones at
  # risk in the pooled table (was_active_before is TRUE where any completion had the zone
  # already invaded, and those zones are not ranked, as an active zone is not per arm).
  .rr_by <- intersect(c("method", "horizon", "training_cutoff"), names(reach_pooled))
  if (all(c("mu_forecast", "rr_nat_rank") %in% names(reach_pooled)))
    reach_pooled <- reach_pooled %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(.rr_by))) %>%
      dplyr::mutate(rr_nat_rank = dplyr::min_rank(dplyr::desc(
        ifelse(was_active_before %in% TRUE, NA_real_, mu_forecast)))) %>%
      dplyr::ungroup()
  if (all(c("priority", "priority_rank") %in% names(reach_pooled)))
    reach_pooled <- reach_pooled %>%
      dplyr::group_by(dplyr::across(dplyr::any_of("horizon"))) %>%
      dplyr::mutate(priority_rank = dplyr::min_rank(dplyr::desc(
        ifelse(was_active_before %in% TRUE, NA_real_, priority)))) %>%
      dplyr::ungroup()
  readr::write_csv(reach_pooled,
                   file.path(POOL, "cascade", "tables", "cascade_reach_scores_all_zones.csv"))

  # -- urban scenario summary: expected_added carries a real SE -----------------
  us <- .need_all(.ctab("urban_scenario_summary.csv"), "urban_scenario_summary.csv")
  us$.u <- us$expected_added_se^2
  us_keys <- intersect(c("city", "seed_week", "seed_r", "seed_r_label"), names(us))
  us_p <- .pool_by(us, us_keys, "expected_added", ".u")
  us_other <- us %>% dplyr::group_by(dplyr::across(dplyr::all_of(us_keys))) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric) & !dplyr::any_of(c("mi_arm", ".u")),
                                   ~ mean(.x, na.rm = TRUE)),
                     paired = all(paired %in% TRUE), .groups = "drop")
  urban_pooled <- us_p %>%
    dplyr::transmute(dplyr::across(dplyr::all_of(us_keys)),
                     expected_added = q, expected_added_se = se,
                     expected_added_lo = lo, expected_added_hi = hi,
                     expected_added_within_sd = sqrt(pmax(U, 0)),
                     expected_added_between_sd = sqrt(pmax(B, 0)),
                     expected_added_lambda = lambda, expected_added_fmi = fmi,
                     expected_added_df = df, n_imputations = n_arms) %>%
    dplyr::left_join(us_other %>% dplyr::select(-dplyr::any_of(
                       c("expected_added", "expected_added_se",
                         "expected_added_lo", "expected_added_hi"))), by = us_keys)
  readr::write_csv(urban_pooled,
                   file.path(POOL, "cascade", "tables", "urban_scenario_summary.csv"))

  # -- per-city impact tables: p_base / p_cond / delta ---------------------------
  .city_files <- unique(unlist(lapply(arms, function(a)
    basename(list.files(file.path(a$cascade_dir, "cascade", "tables"),
                        pattern = "^urban_impact_.*\\.csv$")))))
  .city_files <- .city_files[!grepl("summary", .city_files)]
  if (!length(.city_files))
    stop("[mi-pool] no urban_impact_*.csv in the cascade arms; Figure 4 panel C needs them.",
         call. = FALSE)
  for (cf in .city_files) {
    im <- .need_all(.ctab(cf), cf)
    im_keys <- intersect(c("city", "health_zone", "province", "seed_week"), names(im))
    .p1 <- .pool_by(im %>% dplyr::mutate(.u = .u_from_ci(p_base_lo, p_base_hi)),
                    im_keys, "p_base", ".u")
    .p2 <- .pool_by(im %>% dplyr::mutate(.u = .u_from_ci(p_cond_lo, p_cond_hi)),
                    im_keys, "p_cond", ".u")
    .p3 <- .pool_by(im %>% dplyr::mutate(.u = delta_se^2), im_keys, "delta", ".u")
    im_other <- im %>% dplyr::group_by(dplyr::across(dplyr::all_of(im_keys))) %>%
      dplyr::summarise(rr = mean(rr, na.rm = TRUE),
                       dplyr::across(dplyr::any_of(c("cases_base", "cases_cond", "cases_delta")),
                                     ~ mean(.x, na.rm = TRUE)),
                       exact_zero = all(exact_zero %in% TRUE),
                       # share of completions in which the contrast was flagged elevated.
                       # NOT named `elevated`: 38/40 read that column as the LOGICAL test
                       # (`elevated %in% TRUE`), so a share there admitted only zones elevated
                       # in every completion. The pooled test is set below from the MI interval.
                       elevated_share = mean(elevated %in% TRUE),
                       .groups = "drop")
    # RANK WITHIN (arm, city, seed_week) — NOT within (arm, city). The impact tables carry
    # one block per seed week, so grouping without it ranks a city's zones across seed weeks
    # together and the "made the top 5" share stops describing any single scenario.
    .rk_grp <- unique(c("mi_arm", setdiff(im_keys, c("health_zone", "province"))))
    .rank_share <- im %>% dplyr::group_by(dplyr::across(dplyr::all_of(.rk_grp))) %>%
      dplyr::arrange(dplyr::desc(delta), .by_group = TRUE) %>%
      dplyr::mutate(.rk = dplyr::row_number()) %>% dplyr::ungroup() %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(im_keys))) %>%
      dplyr::summarise(top5_share = mean(.rk <= 5L), .groups = "drop")
    im_pooled <- .p1 %>%
      dplyr::transmute(dplyr::across(dplyr::all_of(im_keys)), p_base = q,
                       p_base_lo = pmax(0, lo), p_base_hi = pmin(1, hi)) %>%
      dplyr::left_join(.p2 %>% dplyr::transmute(dplyr::across(dplyr::all_of(im_keys)),
                         p_cond = q, p_cond_lo = pmax(0, lo), p_cond_hi = pmin(1, hi)),
                       by = im_keys) %>%
      dplyr::left_join(.p3 %>% dplyr::transmute(dplyr::across(dplyr::all_of(im_keys)),
                         delta = q, delta_se = se, delta_lo = lo, delta_hi = hi,
                         delta_within_sd = sqrt(pmax(U, 0)),
                         delta_between_sd = sqrt(pmax(B, 0)),
                         delta_lambda = lambda, delta_fmi = fmi, delta_df = df,
                         n_imputations = n_arms),
                       by = im_keys) %>%
      dplyr::left_join(im_other, by = im_keys) %>%
      dplyr::left_join(.rank_share, by = im_keys) %>%
      # THE POOLED "MATERIALLY ELEVATED" TEST, the per-arm rule (38_urban_scenarios.R: the 90%
      # paired interval lies above zero) applied to the MI total interval of delta.
      dplyr::mutate(elevated = is.finite(delta_lo) & delta_lo > 0)
    readr::write_csv(im_pooled, file.path(POOL, "cascade", "tables", cf))
  }
  # The table 40 reads but does not pool — the seeded-hub provenance, urban_hub_selection.csv,
  # written by 38_urban_scenarios.R::urban_write_hub_selection() and read by
  # 40_cascade_next_dominoes.R::read_urban_hubs() for Figure 4 panel A. Copied from the first
  # arm (arm 0, whose hub set is the published one). The names this loop used to copy
  # (urban_hub_zones.csv / urban_seed_zones.csv) are written by no script, so it skipped
  # silently and the pooled panel A fell back to 40's frozen URBAN_HUBS_FALLBACK. A missing
  # file is now a refusal (also checked in section 1, before anything was written).
  for (rel in "urban_hub_selection.csv") {
    fs <- file.path(vapply(arms, function(a) a$cascade_dir, character(1)),
                    "cascade", "tables", rel)
    if (!all(file.exists(fs)))
      stop("[mi-pool] ", rel, " is missing from cascade arm(s) ",
           paste(vapply(arms, function(a) a$tag, character(1))[!file.exists(fs)], collapse = ", "),
           "; the pooled Figure 4 panel A would fall back to the frozen hub list.", call. = FALSE)
    .tabs <- lapply(fs, function(f) as.data.frame(readr::read_csv(f, show_col_types = FALSE)))
    # The hub derivation reads population, density and mobility only — no line-list input — so
    # the table is imputation-invariant. as.data.frame() strips readr's `spec`/`problems`
    # attributes, which differ whenever column types are inferred differently and would
    # otherwise trigger a spurious "differs between cascade arms" warning on identical content.
    .h <- unique(vapply(.tabs, rlang::hash, character(1)))
    if (length(.h) != 1L)
      warning("[mi-pool] ", rel, " differs between cascade arms; copying arm ",
              arms[[1]]$tag, "'s.", call. = FALSE)
    # What panel A actually uses is the SEEDED set (40 keeps rows flagged `seeded`, else
    # `selected`). If that differs, the cascades seeded different places per arm and Figure 4
    # would pool unlike scenarios: refuse.
    .seeded <- unique(lapply(.tabs, function(s) {
      flag <- if ("seeded" %in% names(s)) "seeded" else "selected"
      if (!all(c("city", "Nom", flag) %in% names(s))) return(NA_character_)
      k <- which(s[[flag]] %in% TRUE)
      sort(paste(s$city[k], s$Nom[k], sep = " / "))
    }))
    if (length(.seeded) != 1L)
      stop("[mi-pool] the cascade arms seeded DIFFERENT hub zones (", rel, "); Figure 4 ",
           "would pool unlike scenarios.", call. = FALSE)
    if (!file.copy(fs[1], file.path(POOL, "cascade", "tables", rel), overwrite = TRUE))
      stop("[mi-pool] could not copy ", rel, " into the pooled mirror.", call. = FALSE)
  }
  message("[mi-pool] pooled cascade tables written")
} else {
  message("[mi-pool] cascade tables NOT pooled (no cascade arms) — Figure 4 cannot be drawn.")
}

# =============================================================================
# 9. PROVENANCE — run_info.json + model_selection.json for the pooled mirror
# =============================================================================
# model_selection.json is COPIED, not recomputed: the featured model is pinned across arms by
# design (D3), and the figure scripts read this file to learn which method to draw.
file.copy(file.path(OUT, "key_outputs", "model_selection.json"),
          file.path(POOL, "key_outputs", "model_selection.json"), overwrite = TRUE)
ri <- jsonlite::fromJSON(file.path(OUT, "key_outputs", "run_info.json"), simplifyVector = TRUE)
ri$onset_mi <- list(
  n_imputations = M,
  seed_base = as.integer(arms[[1]]$seed_base),
  seed_min = min(.seeds), seed_max = max(.seeds), seeds = .seeds,
  # `seeds` are the ONSET seeds (ONSET_IMPUTE_SEED), one per completion; the model seed is one.
  onset_impute_seeds = .seeds, random_seed = as.integer(POOL_RANDOM_SEED),
  stratum = POOL_STRATUM,
  method = POOL_METHOD,
  # What section 1 verified against the published run before pooling.
  analysis_date = as.character(POOL_ADATE), linelist_snapshot = POOL_SNAPSHOT,
  onset_imputation = POOL_IMPUTE,
  # The same mechanism, structured, so 49's captions describe what the arms actually ran
  # (identical on every arm and in this session; checked in section 1).
  onset_impute = arms[[1]]$onset_impute,
  # The probability scale of the Figure-3 mixture (section 1): "recalibrated-deployed-mi" or
  # "raw-mi".
  prob_scale = POOL_PROB_SCALE,
  # Which evaluation metrics carry a within-arm U, and the level of the interval it came from.
  eval_ci = list(level = EV_CI_LEVEL, metrics = unname(.ev_u),
                 resampling = "health zones (zone-cluster bootstrap, evaluate_invasion())"),
  draws_per_arm = .nd,
  pooled_draws = sum(.wu),
  delay_fit = arms[[1]]$delay_imputation,
  pooling = list(
    figure3 = "exact posterior draw mixture (.mix_and_summarise_draws, equal weight per completion)",
    figure2 = "Rubin: T = Ubar + (1 + 1/M) B; Ubar from the published 90% bootstrap band",
    figure4 = "Rubin: T = Ubar + (1 + 1/M) B; Ubar from p_sd / expected_added_se / delta_se"),
  already_invaded_rule = paste("a zone already invaded in a completion contributes p = 0 to",
                              "the mixture; masked to NA only where it is invaded in ALL"),
  cascade = HAS_CASCADE,
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
jsonlite::write_json(ri, file.path(POOL, "key_outputs", "run_info.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 8)

# =============================================================================
# 10. VARIABILITY TABLES
# =============================================================================
# -- data layer ---------------------------------------------------------------
.dw_by_tag <- stats::setNames(dws, vapply(dws, function(d) d$tag, character(1)))
.arm_delta <- function(a, h) {
  d <- a$deploy_delta
  if (is.null(d)) d <- .dw_by_tag[[a$tag]]$deploy_delta
  if (is.null(d) || !h %in% names(d)) NA_real_ else unname(d[[h]])
}
dl <- dplyr::bind_rows(lapply(arms, function(a) tibble::tibble(
  arm = a$arm, seed = a$seed, fingerprint = substr(a$fingerprint, 1L, 12L),
  n_confirmed = a$n_confirmed, n_onset_imputed = a$n_onset_imputed_confirmed,
  # NA for an arm written before these were recorded; 49 then refuses to describe them.
  n_sitrep_confirmed = a$n_sitrep_confirmed %||% NA_integer_,
  n_onset_imputed_linelist = a$n_onset_imputed_linelist %||% NA_integer_,
  n_zones_invaded = nrow(a$first_onset), n_affected_now = length(a$affected_now),
  trunc_asof_mean = a$trunc_asof$mean,
  # The deployed (extract) truncation is refitted per arm (each replays the vintages under its own
  # onset seed); NA for an arm written before it was recorded.
  trunc_deployed_mean = if (is.null(a$trunc_deployed$mean)) NA_real_ else a$trunc_deployed$mean,
  # Read from the ARM file where present; arms written before deploy_delta was recorded
  # there fall back to the draws file, which has always carried it.
  delta_h1 = .arm_delta(a, "1"), delta_h2 = .arm_delta(a, "2"),
  minutes_lfo = a$minutes_lfo, minutes_cascade = a$minutes_cascade,
  minutes_total = a$minutes_total)))
readr::write_csv(dl, file.path(MI_OUT, "tables", "onset_mi_datalayer.csv"))

# -- per-zone invasion-week movement ------------------------------------------
# The OUTCOME itself, across completions. This is the quantity that makes the whole study
# necessary: a zone whose invasion week moves between completions is a zone whose label the
# model is being scored against differently in each one.
fo <- dplyr::bind_rows(lapply(arms, function(a)
  a$first_onset %>% dplyr::mutate(arm = a$arm))) %>%
  dplyr::mutate(week = as.integer(floor(as.numeric(first_onset -
                  lubridate::floor_date(OUTBREAK_START, "week",
                                        week_start = get0("WEEK_ANCHOR", ifnotfound = 1L))) / 7)) + 1L)
wk_shift <- fo %>% dplyr::group_by(health_zone) %>%
  dplyr::summarise(n_arms = dplyr::n(),
                   week_min = min(week), week_max = max(week),
                   week_median = stats::median(week),
                   week_sd = stats::sd(week),
                   week_range = max(week) - min(week),
                   onset_min = min(first_onset), onset_max = max(first_onset),
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(week_range))
readr::write_csv(wk_shift, file.path(MI_OUT, "tables", "onset_mi_invasion_week_shift.csv"))
.moved <- sum(wk_shift$week_range > 0)
message(sprintf("[mi-pool] invasion week moves in %d/%d zones (median range %.1f wk, max %d)",
                .moved, nrow(wk_shift), stats::median(wk_shift$week_range),
                max(wk_shift$week_range)))

# -- live-forecast rank stability ---------------------------------------------
rk <- dplyr::bind_rows(lapply(arms, function(a) {
  r <- a$risk_scores
  if (is.null(r)) return(NULL)
  r %>% dplyr::filter(horizon == 2L, !(was_active_before %in% TRUE),
                      is.finite(p_case_invasion)) %>%
    dplyr::distinct(health_zone, .keep_all = TRUE) %>%
    dplyr::mutate(rank = rank(-p_case_invasion, ties.method = "min"), arm = a$arm) %>%
    dplyr::select(arm, health_zone, p_case_invasion, rank)
}))
# THE SHARE IS OVER COMPLETIONS, as FigureMI3's axis says ("share of completions"), and it
# follows the D5 rule: a zone already invaded in a completion has p = 0 there and is NOT on
# that completion's watch-list, so it counts as not top-K rather than dropping out of the
# denominator. (mean(rank <= 20) over the arms where the zone was still at risk inflated the
# share of exactly the zones whose invasion status moves between completions.) The
# denominator is the completions that carried risk scores — every arm unless one's risk layer
# failed (47 only warns), which is reported. n_arms = completions in which the zone was at
# risk; the rank columns describe those completions only.
.n_rk <- dplyr::n_distinct(rk$arm)
if (.n_rk < M)
  warning(sprintf(paste0("[mi-pool] only %d of %d arms carried risk scores; the rank-stability ",
                         "shares are over those %d completions."), .n_rk, M, .n_rk),
          call. = FALSE)
rank_stab <- rk %>% dplyr::group_by(health_zone) %>%
  dplyr::summarise(n_arms = dplyr::n(), n_completions = .n_rk,
                   p_mean = mean(p_case_invasion), p_sd = stats::sd(p_case_invasion),
                   rank_median = stats::median(rank), rank_min = min(rank), rank_max = max(rank),
                   share_top10 = sum(rank <= 10L) / .n_rk, share_top20 = sum(rank <= 20L) / .n_rk,
                   .groups = "drop") %>%
  dplyr::arrange(rank_median)
readr::write_csv(rank_stab, file.path(MI_OUT, "tables", "onset_mi_rank_stability.csv"))
# Spearman of each arm's ranking against the lowest-numbered arm with risk scores, a single number for
# "how much does the watch-list move".
# ARM 0 IS EXCLUDED FROM ITS OWN COMPARISON. Its Spearman against itself is 1 by
# construction, and leaving it in made the reported range read "1.000-1.000" whenever the
# other arms happened to agree closely — a summary that cannot distinguish perfect stability
# from a sweep of size one.
.arm0 <- min(rk$arm)
.a0 <- rk %>% dplyr::filter(arm == .arm0) %>% dplyr::select(health_zone, r0 = rank)
rho <- rk %>% dplyr::filter(arm != .arm0) %>%
  dplyr::inner_join(.a0, by = "health_zone") %>%
  dplyr::group_by(arm) %>%
  dplyr::summarise(spearman_vs_arm0 = stats::cor(rank, r0, method = "spearman"),
                   n = dplyr::n(), .groups = "drop") %>%
  # The reference arm, so a reader (and 49's caption) never has to assume it is arm 0.
  dplyr::mutate(ref_arm = .arm0)
if (!nrow(rho))
  rho <- tibble::tibble(arm = integer(0), spearman_vs_arm0 = numeric(0), n = integer(0),
                        ref_arm = integer(0))
readr::write_csv(rho, file.path(MI_OUT, "tables", "onset_mi_rank_agreement.csv"))

# =============================================================================
# 11. SUMMARY
# =============================================================================
# THE FEATURED MODEL ONLY. ev_var also carries the three structural baselines (they are
# scored, not merely drawn), so without this filter every metric appeared FOUR times with no
# column saying which method each row belonged to — four different AUC-PRs stacked under one
# label. The baselines keep their rows in the CSV; the report is about the featured model.
.hl <- ev_var %>%
  dplyr::filter(method == POOL_METHOD,
                metric %in% c("auc_pr", "auc_pr_skill", "auc_roc", "mean_rank_of_truth",
                              "base_rate", "n_invasions", "n_atrisk")) %>%
  dplyr::arrange(horizon, metric)
md <- c(
  "# Onset multiple imputation — variability of the featured model",
  "",
  sprintf("Generated %s. M = %d completions, onset seeds (ONSET_IMPUTE_SEED) %d-%d with the model seed RANDOM_SEED = %s held on every arm, stratum `%s`, model `%s`.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), M, min(.seeds), max(.seeds), POOL_RANDOM_SEED,
          POOL_STRATUM, POOL_METHOD),
  sprintf("Analysis date %s; line-list snapshot `%s`. Onset imputation (identical on every arm and in this session): %s.",
          POOL_ADATE, POOL_SNAPSHOT, POOL_IMPUTE),
  sprintf("Onset->sample delay recorded by every arm (identical across arms; the nowcast's fallback delay, not used by the onset imputation): %s, mean %.3f d, SD %.3f d (%s, n = %s).",
          arms[[1]]$delay_imputation$family, arms[[1]]$delay_imputation$mean,
          arms[[1]]$delay_imputation$sd, arms[[1]]$delay_imputation$estimator,
          format(arms[[1]]$delay_imputation$n_fit)),
  "",
  "## What moves",
  "",
  sprintf("- The invasion WEEK moves in **%d of %d** zones (median range %.1f weeks, max %d).",
          .moved, nrow(wk_shift), stats::median(wk_shift$week_range), max(wk_shift$week_range)),
  sprintf("- The confirmed-case count is constant at %s across every completion, as it must be:",
          format(arms[[1]]$n_confirmed, big.mark = ",")),
  "  the imputation moves onset DATES, never the number of records.",
  if (nrow(rho)) local({
    # "The published draw" only when the reference arm carries the published run's onset seed
    # (as 49's caption decides); arm 0 need not have carried risk scores.
    .ref_seed <- arms[[match(.arm0, vapply(arms, function(a) as.integer(a$arm), integer(1)))]]$seed
    .pub_seed <- suppressWarnings(as.integer(.ri_pub$onset_imputation$onset_impute_seed))
    .is_pub <- length(.pub_seed) == 1L && !is.na(.pub_seed) && identical(as.integer(.ref_seed), .pub_seed)
    sprintf(paste0("- Live 2-week ranking of the other %d completion(s) against arm %d (onset seed ",
                   "%d%s): Spearman %.3f-%.3f."),
            nrow(rho), .arm0, as.integer(.ref_seed),
            if (.is_pub) ", the published draw" else "",
            min(rho$spearman_vs_arm0, na.rm = TRUE), max(rho$spearman_vs_arm0, na.rm = TRUE))
  })
  else "- Live 2-week ranking: only one completion carried risk scores; no comparison.",
  "",
  "## Headline metrics (Rubin)",
  "",
  "`q` = MI estimate; `B` = between-imputation variance; `U` = mean within-imputation",
  "variance; `T` = U + (1 + 1/M)B, the total. `lambda` = (1 + 1/M)B / T is the PROPORTION OF",
  "TOTAL VARIANCE that is BETWEEN completions; `fmi` is that with Rubin's small-sample",
  "degrees-of-freedom correction, `fmi = lambda + (1 - lambda) * 2/(df + 3)`, which is never",
  "the smaller of the two. Intervals are 90% on t with `df` degrees of freedom",
  "(Rubin 1987 §3.1), NOT normal: at df = 19 that is 1.729 against 1.645, so a normal",
  "interval would be 5% too narrow.",
  "",
  "**`lambda` slightly overstates the onset-date share.** The completions differ only in their",
  "onset seed; the model seed is the same on every arm. But each completion's Stan chains,",
  "EpiNow2 R(t) and epinowcast fits, evaluation bootstraps and recalibration-delta draws run on",
  "a different completed data set, so their finite-draw Monte Carlo error enters B along with",
  "the onset imputation.",
  "",
  "Metrics with `u_known = FALSE` have no published within-arm interval, so U is taken as 0,",
  "T = (1 + 1/M)B and lambda = 1 BY CONSTRUCTION rather than by measurement: for those rows",
  "the two components cannot be separated and the interval is the between-imputation spread",
  sprintf(paste0("alone. The metrics that carry a real U are those with a published within-arm ",
                 "interval: %s (each from its %.0f%% zone-cluster bootstrap interval in ",
                 "evaluate_invasion(), which resamples health zones)."),
          if (length(.ev_u)) paste(sprintf("`%s`", .ev_u), collapse = ", ") else "none",
          100 * EV_CI_LEVEL),
  "",
  "```",
  # WIDEN THE CONSOLE FOR THE CAPTURE. print.data.frame wraps at getOption("width") (80 by
  # default), which split the last column into a second block of 56 bare TRUE/FALSE lines
  # with nothing joining them back to their rows. Restored afterwards.
  paste(local({
    .w <- getOption("width"); on.exit(options(width = .w), add = TRUE)
    options(width = 200)
    utils::capture.output(print(as.data.frame(
      .hl %>% dplyr::transmute(horizon, metric,
                               q = signif(q, 5), B = signif(B, 4), U = signif(U, 4),
                               T = signif(T, 4), lambda = round(lambda, 3),
                               fmi = round(fmi, 3), df = round(df, 1), u_known)),
      row.names = FALSE))
  }), collapse = "\n"),
  "```",
  "",
  "## Reading the skill numbers",
  "",
  "`auc_pr_skill = auc_pr / base_rate`, and `base_rate` is computed from the IMPUTED labels,",
  "so its spread across completions is uncertainty about WHAT HAPPENED, not sampling error",
  "about a fixed truth. Read every skill figure beside its `base_rate`, `n_invasions` and",
  "`n_atrisk` rows above.",
  "",
  "## Products",
  "",
  sprintf("- Pooled artifact mirror: `%s`", POOL),
  sprintf("- Variability tables: `%s`", file.path(MI_OUT, "tables")),
  "- Figures: `outputs/key_outputs/manuscript_figures/mi/` (49_onset_mi_figures.R)")
writeLines(md, file.path(MI_OUT, "ONSET_MI_REPORT.md"))
unlink(.pool_marker)   # every pooled product is written

message("\n[mi-pool] DONE")
message("  pooled mirror : ", POOL)
message("  tables        : ", file.path(MI_OUT, "tables"))
message("  report        : ", file.path(MI_OUT, "ONSET_MI_REPORT.md"))

}   # end of the script-run guard

