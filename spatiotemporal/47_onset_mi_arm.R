# =============================================================================
# 47_onset_mi_arm.R — ONE arm of the onset multiple-imputation sweep
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Runs the featured model end to end on ONE onset completion and writes everything
# 48_onset_mi_pool.R needs to pool it. `run_onset_mi.sh` drives the arms;
# 48 pools them; 49 draws the MI diagnostics. Design: PLAN_ONSET_MI.md.
#
# WHAT AN ARM IS. load_linelist() (01_data_prep.R) imputes the confirmed onsets it cannot
# use (the realised share is published as .PCT_ONSET_IMPUTED; ~24% on the 2026-10-09 frame)
# as the sample date minus an interval drawn from DONORS: records of the same class sampled
# in the same week whose onset is recorded. Each record's draw is a deterministic hash of
# (alert_id, ONSET_IMPUTE_SEED), so changing ONSET_IMPUTE_SEED gives an independent completion.
# One arm = one onset seed = one completed data set. The draw moves the PREDICTORS (onset-week counts ->
# Lambda, R(t), the nowcast), the OUTCOME (a zone's invasion week IS its first confirmed onset
# week) and the SUPPORT (which cells are at risk) simultaneously, which is why the arm re-runs
# the model rather than re-scoring one set of predictions. The mechanism this arm actually
# ran under is recorded in its file (onset_impute) and 48 refuses to pool arms that differ.
#
# TWO SEEDS, AND ONLY ONE OF THEM VARIES. ONSET_IMPUTE_SEED (00_config.R) is the only seed the
# missing-onset draw reads; RANDOM_SEED seeds every Monte-Carlo step downstream of the line
# list — the brms/Stan chains (21), the EpiNow2 R(t) and epinowcast fits (02, 04b), the
# evaluation bootstraps (16) and the per-draw recalibration delta (21). Arm m sets
# ONSET_IMPUTE_SEED = base + m and holds RANDOM_SEED at the PUBLISHED run's value for every arm
# (they used to vary RANDOM_SEED, which moved both at once). WHAT B STILL CONTAINS: each arm's
# fits run on a DIFFERENT completed data set, so their finite-draw Monte Carlo error does not
# cancel between arms and enters the between-arm spread too. lambda = (1+1/M)B/T therefore
# slightly OVERSTATES the share of variance due to not knowing the onset dates.
#
# WHY A SEPARATE HARNESS AND NOT run_all.R. The sweep needs ONE model under M completions.
# run_all.R cross-validates the whole grid and then runs the figure, report and
# cascade layers, none of which this question touches, at ~100 min an arm. This harness
# mirrors run_all.R's wiring call-site by call-site (.assert_matches_run_all() below checks the
# back-test, forms back-test, deployed-suite and recalibration call sites and fails
# loudly if the two drift) and stops at the artifacts 48 pools.
#
# WHY THE SEED GOES IN THE ENVIRONMENT. 00_config.R carries no load guard and is re-sourced
# by most pipeline modules and by run_cascade.R. An override applied after source() survives
# only inside a process that controls every source() itself — and this arm shells out to
# run_cascade.R. ONSET_IMPUTE_SEED and RANDOM_SEED are therefore read from the environment by
# 00_config.R, set here BEFORE the first source(), and asserted after it. If the onset-seed
# assertion ever fails, every arm is secretly the same draw and the whole sweep reports a
# zero-width interval that looks like precision.
#
# Env:
#   ONSET_MI_ARM        arm index m (integer >= 0; REQUIRED). Arm 0 = the published draw.
#   ONSET_MI_SEED_BASE  onset-seed base (default 20260704, the project seed). Arm m sets
#                       ONSET_IMPUTE_SEED = base + m; RANDOM_SEED is the published run's.
#   ONSET_MI_OUT        output root (default outputs/onset_mi)
#   ONSET_MI_DRAWS      posterior draws retained per arm (default 1000)
#   ONSET_MI_CASCADE    1 (default) to run the cascade arm, 0 to skip
#   ONSET_MI_MOBILITY   shared mobility cache rds (built once, reused by every arm)
# =============================================================================

suppressPackageStartupMessages({library(dplyr)})

ROOT   <- "/Users/mills/Documents/BDBV2026-Analysis"
ST_DIR <- file.path(ROOT, "spatiotemporal")

# ---- arm identity ------------------------------------------------------------
.arm_env <- trimws(Sys.getenv("ONSET_MI_ARM", ""))
if (!nzchar(.arm_env)) stop("[mi-arm] ONSET_MI_ARM is required.", call. = FALSE)
MI_ARM <- suppressWarnings(as.integer(.arm_env))
if (is.na(MI_ARM) || MI_ARM < 0L)
  stop("[mi-arm] ONSET_MI_ARM must be a non-negative integer, got ", shQuote(.arm_env),
       call. = FALSE)
MI_SEED_BASE <- suppressWarnings(as.integer(Sys.getenv("ONSET_MI_SEED_BASE", "20260704")))
if (is.na(MI_SEED_BASE)) stop("[mi-arm] ONSET_MI_SEED_BASE is not an integer.", call. = FALSE)
MI_SEED   <- MI_SEED_BASE + MI_ARM
MI_OUT    <- Sys.getenv("ONSET_MI_OUT", file.path(ST_DIR, "outputs", "onset_mi"))
MI_DRAWS  <- suppressWarnings(as.integer(Sys.getenv("ONSET_MI_DRAWS", "1000")))
if (is.na(MI_DRAWS) || MI_DRAWS < 50L)
  stop("[mi-arm] ONSET_MI_DRAWS must be an integer >= 50.", call. = FALSE)
MI_CASCADE <- !identical(trimws(Sys.getenv("ONSET_MI_CASCADE", "1")), "0")
MOB_CACHE  <- Sys.getenv("ONSET_MI_MOBILITY", file.path(MI_OUT, "_mobility_cache.rds"))
# SIDE EFFECTS ONLY WHEN THIS FILE RUNS AS AN ARM (Rscript, or FORCE_RUN=1 -- the rule of the
# script_run_or_forced() guard below, which cannot be called before 00_config.R is sourced). A
# source() -- a test, an ad-hoc check -- must not set the caller's environment variables, its
# future plan, or create directories.
.run_arm <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  (length(f) == 1L && identical(basename(f), "47_onset_mi_arm.R")) ||
    tolower(trimws(Sys.getenv("FORCE_RUN", ""))) %in% c("1", "true", "t", "yes", "y")
})
if (.run_arm) dir.create(MI_OUT, recursive = TRUE, showWarnings = FALSE)

# SEED + STRATUM INTO THE ENVIRONMENT, BEFORE ANY source(). The stratum is pinned
# EXPLICITLY rather than left to 00_config.R's default: it is an environment variable whose
# two settings produce two different LABEL SETS ("a switch and not a fix"), so a stray
# export in the calling shell — or inheritance from another sensitivity run in the same
# session — would put some arms on one estimand and some on another, and the pooled result
# would be a mixture of two studies with nothing in the output saying so.
# THE ONSET SEED VARIES PER ARM; THE MODEL SEED DOES NOT. RANDOM_SEED is the published run's,
# read from its run_info.json (write_run_info.R records onset_imputation$random_seed); a
# run_info that predates that field leaves it to 00_config.R's default, which is what an
# ordinary run uses. A stray RANDOM_SEED exported in the calling shell is therefore overridden
# either way, so it can neither vary between arms nor differ from the published model seed.
MODEL_SEED <- local({
  .r <- tryCatch(jsonlite::fromJSON(file.path(ST_DIR, "outputs", "key_outputs", "run_info.json"),
                                    simplifyVector = TRUE)$onset_imputation$random_seed,
                 error = function(e) NULL)
  .v <- if (length(.r) == 1L) suppressWarnings(as.integer(.r)) else NA_integer_
  if (.run_arm) { if (is.na(.v)) Sys.unsetenv("RANDOM_SEED") else Sys.setenv(RANDOM_SEED = as.character(.v)) }
  .v                                   # NA: resolved by 00_config.R's default, asserted below
})
if (.run_arm) {
  Sys.setenv(ONSET_IMPUTE_SEED = as.character(MI_SEED))
  Sys.setenv(ONSET_IMPUTE_STRATUM = "confirmed")
  # The LFO result cache is keyed by .lfo_stamp, which now carries random_seed — but this
  # harness never reads or writes that cache anyway, and saying so is cheaper than relying on
  # a default staying FALSE.
  Sys.setenv(REUSE_CACHED_LFO = "0")
}

source(file.path(ST_DIR, "00_config.R"))
stop_if_data_definition_sensitivity("47_onset_mi_arm.R")   # hard-coded published paths

if (.run_arm && isTRUE(get0("PARALLEL_JOBS", ifnotfound = 1L) > 1L) &&
    requireNamespace("future", quietly = TRUE)) {
  future::plan(future::multicore, workers = PARALLEL_JOBS)
  options(future.globals.maxSize = 12 * 1024^3)
  message(sprintf("[parallel] multicore: %d workers", PARALLEL_JOBS))
}

# 22_daily_reissue.R is REQUIRED, not optional: it defines linelist_observation_date() and
# reaggregate_asof(), without which the as-of truncation silently returns NULL
# (02_epi_params.R guards on exists()) and the leakage-free folds cannot be rebuilt.
# ONE module list, from 00_config.R. Nine files used to hand-maintain their own and all nine
# differed — see load_pipeline_modules().
load_pipeline_modules("model", dir = ST_DIR)


# THIS FILE IS A SCRIPT, NOT A LIBRARY: everything below EXECUTES and writes published
# outputs. It is launched as its own subprocess. Guarded so a source() — a test, an
# ad-hoc check, a sweep that loads "every module" — defines what it needs and stops,
# instead of silently rewriting files. FORCE_RUN=1 forces a run from a source().
if (!script_run_or_forced("47_onset_mi_arm.R")) {
  message("[mi-arm] sourced, not run: nothing is written. Use `Rscript 47_onset_mi_arm.R`",
          " (or FORCE_RUN=1) to run it.")
} else {

# ---- THE ASSERTIONS THAT MAKE THE REPLICATES REAL ---------------------------
stopifnot(is.numeric(RANDOM_SEED), length(RANDOM_SEED) == 1L, !is.na(RANDOM_SEED),
          exists("ONSET_IMPUTE_SEED"), length(ONSET_IMPUTE_SEED) == 1L, !is.na(ONSET_IMPUTE_SEED))
if (!identical(as.integer(ONSET_IMPUTE_SEED), as.integer(MI_SEED)))
  stop(sprintf(paste0("[mi-arm] ONSET_IMPUTE_SEED is %s after sourcing but this arm asked for %d. ",
                      "A later source() reset it, so this arm would repeat another arm's ",
                      "onset draw and the sweep's between-imputation variance would be ",
                      "understated. Refusing to run."), ONSET_IMPUTE_SEED, MI_SEED), call. = FALSE)
# The model seed must be the published run's on EVERY arm (48 asserts they agree).
if (!is.na(MODEL_SEED) && !identical(as.integer(RANDOM_SEED), MODEL_SEED))
  stop(sprintf("[mi-arm] RANDOM_SEED is %s after sourcing but the published run used %d.",
               RANDOM_SEED, MODEL_SEED), call. = FALSE)
.arm_stratum <- get0("ONSET_IMPUTE_STRATUM", ifnotfound = NULL)
if (is.null(.arm_stratum) || !identical(as.character(.arm_stratum)[1], "confirmed"))
  stop("[mi-arm] ONSET_IMPUTE_STRATUM resolved to ",
       if (is.null(.arm_stratum)) "the POOLED fit (NULL)" else shQuote(as.character(.arm_stratum)[1]),
       " but this sweep is defined on the 'confirmed' stratum. The two settings produce ",
       "different label sets and must not be mixed within one pooled result.", call. = FALSE)
TAG <- sprintf("%03d", MI_ARM)
message(sprintf("\n=== onset-MI arm %s | ONSET_IMPUTE_SEED = %d | RANDOM_SEED = %d (held) | stratum = %s | draws = %d ===",
                TAG, ONSET_IMPUTE_SEED, RANDOM_SEED, .arm_stratum, MI_DRAWS))

t_arm <- Sys.time()

# ---- data layer (ARM-SPECIFIC: the imputed onsets live here) ----------------
dat        <- prep_all_data()
zones_all  <- names(dat$pop)
pop_vec    <- dat$pop
covariates <- dat$covariates

# THE ARM FINGERPRINT. A hash of this arm's per-zone FIRST CONFIRMED ONSET — the quantity
# the imputation actually moves and the outcome every fold is scored against. 48 asserts all
# M are pairwise distinct. This is the belt-and-braces behind the seed assertion above: it
# catches a silently repeated draw no matter which mechanism caused it (an unexpected
# re-source, a cache serving a stale completion, a driver bug passing the same arm twice).
.first_onset <- dat$ll %>%
  dplyr::filter(confirmed %in% TRUE, !is.na(date_index), !is.na(health_zone)) %>%
  dplyr::group_by(health_zone) %>%
  dplyr::summarise(first_onset = min(date_index), .groups = "drop") %>%
  dplyr::arrange(health_zone)
ONSET_FINGERPRINT <- rlang::hash(list(zone = .first_onset$health_zone,
                                      first = as.integer(.first_onset$first_onset)))
message(sprintf("[mi-arm] onset fingerprint %s (%d zones with a confirmed onset)",
                substr(ONSET_FINGERPRINT, 1L, 12L), nrow(.first_onset)))

# ---- mobility (ARM-INVARIANT) ------------------------------------------------
# Kernels are built from population, OSRM and the Flowminder files; no line-list input, so
# one build serves every arm. Cached on a key covering the inputs, settings and code that can
# change it (below), so an edit cannot be served a stale kernel set.
osrm_mat      <- load_osrm()
# The road-distance matrix exactly as run_all.R hands it to the build (03's single rule: loaded
# only when a road-distance family is switched on).
osrm_dist_mat <- load_osrm_dist_for_mobility()
# THE KEY COVERS EVERYTHING THAT GOVERNS THE KERNELS: the settings 03_mobility_matrices.R records
# in its own config fingerprint (mobility_manifest.csv), the remaining switches
# build_all_mobility_matrices() reads, and the content of every input it builds from
# (mobility_input_fingerprint(), 03: population, OSRM, the relocation OD table, the short-trip
# and cohort profiles as resolved at the analysis date, the cohort origins, the epicentre
# zones). The values are read from THIS session rather than from the manifest on disk: the
# manifest describes whatever run last wrote outputs/mobility/, and a cache key must describe
# the kernels this arm would build.
.mob_cfg <- lapply(stats::setNames(nm = c(
  "FLOWMINDER_OD_FILE", "MOBILITY_SOURCE_FILL", "MOBILITY_FILL_QMAX", "COHORT_WINDOW",
  "OSRM_GAP_FILL", "M4_CENSORED_FIT", "EPICENTRE_ZONES",
  "COHORT_SOURCES", "FLOWMINDER_ST_TAGS", "MOBILITY_PRIMARY",
  "INCLUDE_OSRM_DIST_MODELS", "INCLUDE_COHORT_MODELS", "INCLUDE_FLOWSTATIC_MODELS",
  "INCLUDE_COHORT_SPLIT_MODELS", "INCLUDE_COHORT_GRAVITY_MODELS", "INCLUDE_SOURCEFILL_MODELS",
  "INCLUDE_UNFILLED_MODELS", "INCLUDE_LEGACY_KERNELS", "INCLUDE_M9_MODELS", "INCLUDE_M15_MODELS")),
  function(nm) get0(nm, ifnotfound = NA))
.mob_key <- rlang::hash(list(zones = zones_all, pop = pop_vec,
                             analysis_date = ANALYSIS_DATE,
                             # the CODE that builds the kernels: an edit to 03 must not be served
                             # kernels built by the previous code.
                             code = unname(tools::md5sum(file.path(ST_DIR, "03_mobility_matrices.R"))),
                             inputs = mobility_input_fingerprint(zones_all, pop_vec, ANALYSIS_DATE,
                                                                 osrm_mat, osrm_dist_mat),
                             cfg = .mob_cfg))
mobility_matrices <- NULL
if (file.exists(MOB_CACHE)) {
  .mc <- readRDS(MOB_CACHE)
  if (identical(.mc$key, .mob_key)) {
    mobility_matrices <- .mc$matrices
    message("[mobility] reused cached matrices: ", paste(names(mobility_matrices), collapse = ", "))
  } else message("[mobility] cache key mismatch — rebuilding.")
}
if (is.null(mobility_matrices)) {
  # build_all_mobility_matrices() SAVES every kernel, its summary and the manifest under the
  # global OUT_MOBILITY -- the PUBLISHED outputs/mobility/ tree, which the cascade arms and 46
  # read. The sweep's build goes to its own directory instead; the published kernels stay those
  # of the run that wrote them.
  .om <- OUT_MOBILITY
  OUT_MOBILITY <- file.path(MI_OUT, "_mobility")
  dir.create(OUT_MOBILITY, recursive = TRUE, showWarnings = FALSE)
  mobility_matrices <- tryCatch(build_all_mobility_matrices(
    zones_all = zones_all, pop_vec = pop_vec, osrm_mat = osrm_mat,
    analysis_date = ANALYSIS_DATE, osrm_dist_mat = osrm_dist_mat),
    finally = OUT_MOBILITY <- .om)
  saveRDS(list(key = .mob_key, matrices = mobility_matrices), MOB_CACHE)
}
gt_pmfs <- compute_all_gt_pmfs()

# ---- the two truncation regimes, exactly as run_all.R resolves them ---------
# extract : the DEPLOYED nowcast + R(t) truncation (onset -> appearance in the extract),
#           deconvolved from the vintage archive. NOT arm-invariant: the panel is the
#           date_index of every vintage replayed through load_linelist_vintage(), which imputes
#           with this arm's ONSET_IMPUTE_SEED, and the fit's cache is keyed on that series, so
#           each arm gets its own Stan fit (recorded in the arm file, trunc_deployed).
# asof    : the FOLD truncation (onset -> observation), measured by onset cohort on the
#           date_index of records whose recorded onset is USABLE (.trunc_asof_onset(),
#           02_epi_params.R) -- exactly the records this sweep does NOT impute -- so it is
#           ARM-INVARIANT: every arm's re-censored folds are corrected by one shared
#           distribution, which is what makes them comparable. It is refitted per arm anyway
#           (sub-second) and 48 ASSERTS the invariance rather than trusting this note.
.trunc_deployed <- tryCatch(
  .trunc_as_delay_spec(epinow2_truncation(issue_date = NULL, ll = dat$ll), "extract"),
  error = function(e) { warning("[mi-arm] deployed truncation unavailable (", conditionMessage(e),
                                "); falling back to the onset->sample delay.", call. = FALSE); NULL })
.trunc_asof <- tryCatch(
  .trunc_as_delay_spec(epinow2_truncation(regime = "asof", ll = dat$ll), "asof"),
  error = function(e) { warning("[mi-arm] as-of truncation unavailable: ",
                                conditionMessage(e), call. = FALSE); NULL })
if (is.null(.trunc_asof))
  stop("[mi-arm] the as-of truncation did not fit. run_all.R would fall back to the ",
       "onset->sample delay here and warn; for a pooled comparison that is not acceptable, ",
       "because the arms could then be corrected by different estimators.", call. = FALSE)
message(sprintf("[mi-arm] fold truncation: %s (mean %.2f d, sd %.2f d, %s)",
                .trunc_asof$family, .trunc_asof$mean, .trunc_asof$sd, .trunc_asof$estimator))

zone_week_nc <- apply_nowcast_correction(dat$zone_week, analysis_date = ANALYSIS_DATE,
                                         delay = .trunc_deployed)
all_weeks       <- sort(unique(zone_week_nc$week_start))
training_cutoff <- max(all_weeks)
affected_now    <- affected_zones(zone_week_nc, training_cutoff)
message(sprintf("[mi-arm] %d/%d zones already affected; %d at risk",
                length(affected_now), length(zones_all), length(zones_all) - length(affected_now)))

# ---- the featured model (PINNED, read not re-selected) ----------------------
# D3 of PLAN_ONSET_MI.md: re-selecting per arm would cost the whole candidate grid per arm,
# and a selection that moves between arms would make the pooled figures a mixture of
# different models. Read from the shipped model_selection.json so "the featured model" is
# the pipeline's own CV pick rather than a second hard-coded copy.
# READ BY THE SUITE'S ONE READER, read_featured_model() (figure_style.R): `bayesian` first,
# `headline` as the fallback, no typed fallback model. This file used to carry its own copy of
# that rule; a second copy reading a different field would let the arm FIT one model while the
# pooled figures LABEL another, and nothing downstream compares the two. figure_style.R is
# sourced into a private environment so its palette and theme objects do not enter this session.
read_featured_model <- local({
  .fs <- new.env()
  sys.source(file.path(ST_DIR, "figure_style.R"), envir = .fs)
  .fs$read_featured_model
})
FEATURED_METHOD <- read_featured_model(file.path(ST_DIR, "outputs", "key_outputs",
                                                 "model_selection.json"))
# THE CROSS-VALIDATED GRID, built exactly as run_all.R builds BAYES_LFO_SPECS: gt_arms = FALSE.
# The default (gt_arms = TRUE) appends the generation-time arms of the PREVIOUS run's pick; run_all.R
# instead refits the arms of THIS run's featured model after selection (below).
BAYES_SPECS <- bayes_default_grid(mobility_matrices, gt_arms = FALSE)
# A FEATURED TIME-VARYING FORM (stage 2 of the selection) is resolved to its own spec -- the
# stage-1 specification refitted by the Stan engine -- never to its fixed-beta parent.
.sel_json <- jsonlite::fromJSON(file.path(ST_DIR, "outputs", "key_outputs", "model_selection.json"),
                                simplifyVector = TRUE)
.s1 <- bayes_spec_for_label(.sel_json$featured$bayesian$stage1_method %||% NA_character_, BAYES_SPECS)
sp <- bayes_spec_for_label(FEATURED_METHOD,
                           c(BAYES_SPECS, bayes_tv_forms_of(.s1)))
if (is.null(sp))
  stop(sprintf("[mi-arm] featured model '%s' is neither in this run's Bayesian grid (%d specs) nor a time-varying form of its stage-1 pick.",
               FEATURED_METHOD, length(BAYES_SPECS)), call. = FALSE)
message(sprintf("[mi-arm] featured model: %s (mob=%s gt=%s cov=[%s] link=%s engine=%s form=%s)",
                sp$label, sp$mob, sp$gt, paste(sp$cov, collapse = ","),
                sp$link %||% "cloglog", sp$engine %||% "brms", sp$form %||% "fixed"))
# TWO run_all.R PATHS THIS ARM DOES NOT MIRROR, refused rather than silently diverged from.
if (isTRUE(get0("GT_MARGINALISE_FEATURED", ifnotfound = FALSE)))
  stop("[mi-arm] GT_MARGINALISE_FEATURED is on: run_all.R then publishes a GT-marginalised ",
       "featured forecast, which this arm does not reproduce (it deploys the fixed-GT one).",
       call. = FALSE)
if (grepl("^M11", sp$mob))
  stop("[mi-arm] the featured kernel ", sp$mob, " is built on demand by run_all.R ",
       "(build_inward_contact_matrix()), not by build_all_mobility_matrices(); this arm cannot ",
       "mirror it.", call. = FALSE)

# ---- GUARD: this harness must call run_invasion_lfo() the way run_all.R does -
# A silent divergence would make the sweep measure the harness rather than the imputation.
.assert_matches_run_all <- function() {
  src <- readLines(file.path(ST_DIR, "run_all.R"), warn = FALSE)
  i <- grep("^  run_invasion_lfo\\(zone_week_outbreak", src)
  if (length(i) != 1L) {
    warning("[mi-arm] could not locate run_all.R's LFO call site; the argument guard did ",
            "not run.", call. = FALSE); return(invisible())
  }
  blk <- paste(src[i:(i + 20L)], collapse = " ")
  need <- c("delay = .trunc_asof", "horizons = LFO_HORIZONS", "analysis_date = ANALYSIS_DATE",
            "min_atrisk_events = 0L", "nowcast_fn = apply_nowcast_correction",
            # THE MODEL SET, not only the named arguments. This guard checked every keyword
            # argument and none of the positional ones, so it passed while the arm scored a
            # strictly smaller model set than run_all does — which is how two comparators
            # Figure 2A requires came to be missing from the pooled curves.
            "c(BASELINE_MODELS, bayes_lfo_models)")
  miss <- need[!vapply(need, function(p) grepl(p, blk, fixed = TRUE), logical(1))]
  if (length(miss))
    stop("[mi-arm] run_all.R's LFO call no longer passes: ", paste(miss, collapse = "; "),
         ". Re-sync this harness before trusting the sweep.", call. = FALSE)
  # THE FORMS' OWN BACK-TEST. A featured time-varying form's published rows come from run_all.R's
  # separate forms call, so it is checked against the same arguments (and the Stan effort, which
  # make_suite_lfo_model() below passes as TV_FORM_ITER).
  fi <- grep("^\\s+run_invasion_lfo\\(zone_week_outbreak, \\.form_models", src)
  if (length(fi) != 1L)
    stop("[mi-arm] could not locate run_all.R's forms back-test (run_invasion_lfo(zone_week_outbreak, ",
         ".form_models, ...)); a featured form's arm could be cross-validated differently.", call. = FALSE)
  blkf <- paste(src[fi:(fi + 6L)], collapse = " ")
  needf <- c(setdiff(need, "c(BASELINE_MODELS, bayes_lfo_models)"),
             "min_train_weeks = LFO_MIN_TRAINING_WEEKS")
  missf <- needf[!vapply(needf, function(p) grepl(p, blkf, fixed = TRUE), logical(1))]
  if (length(missf))
    stop("[mi-arm] run_all.R's forms back-test no longer passes: ", paste(missf, collapse = "; "),
         ". Re-sync this harness.", call. = FALSE)
  if (!any(grepl("iter_stan = TV_FORM_ITER, diag_dir = .form_diag_dir", src, fixed = TRUE)))
    stop("[mi-arm] run_all.R's forms are no longer cross-validated at iter_stan = TV_FORM_ITER; ",
         "re-sync make_suite_lfo_model() below.", call. = FALSE)
  # The deployed suite must also still be fitted with the AS-OF truncation and the rolling
  # predictors; the live forecast below mirrors fit_bayes_suite()'s single-model path, and
  # a change to how run_all.R calls it would leave this arm deploying a different estimator.
  j <- grep("fit_bayes_suite\\(zone_week_nc", src)
  if (length(j) == 1L) {
    blk2 <- paste(src[j:(j + 18L)], collapse = " ")
    # The live fit below mirrors fit_bayes_suite()'s single-model path: every model on its own
    # full design, rolling predictors on the as-of truncation, iter = 2000.
    need2 <- c("trunc_delay = .trunc_asof", "linelist = dat$ll", "analysis_date = ANALYSIS_DATE",
               "iter = 2000L")
    miss2 <- need2[!vapply(need2, function(p) grepl(p, blk2, fixed = TRUE), logical(1))]
    if (length(miss2))
      stop("[mi-arm] run_all.R's fit_bayes_suite() call no longer passes: ",
           paste(miss2, collapse = "; "), ". Re-sync the live-forecast block below.",
           call. = FALSE)
  } else warning("[mi-arm] could not locate run_all.R's fit_bayes_suite() call site.",
                 call. = FALSE)
  # THE GRID AND THE FORMS. BAYES_SPECS and the featured form above copy run_all.R's
  # constructions; if either changes, this arm could resolve a different featured spec.
  for (pat in c("BAYES_LFO_SPECS <- bayes_default_grid(mobility_matrices, gt_arms = FALSE)",
                ".form_specs <- bayes_tv_form_specs(.s1)"))
    if (!any(grepl(pat, src, fixed = TRUE)))
      stop("[mi-arm] run_all.R no longer contains `", pat, "`; the featured spec this arm ",
           "resolves may no longer match the published run's. Re-sync.",
           call. = FALSE)
  # THE DEPLOYED RECALIBRATION CARRIES delta's GAMMA POSTERIOR. run_all.R hands the live suite
  # recal_tbl's delta AND gamma_shape/gamma_rate, so .scale_invasion_draws() draws one delta per
  # posterior draw; the arm below must pass the same columns or its forecast fixes delta while
  # the published one integrates over it — and arm 0 would then not reproduce the published table.
  r <- grep("\\.deploy_delta <- recal_tbl\\[", src)
  if (length(r) != 1L)
    stop("[mi-arm] could not locate run_all.R's deployed-recalibration block (.deploy_delta <- ",
         "recal_tbl[...]); the columns the live forecast is recalibrated with are unverified.",
         call. = FALSE)
  blk3 <- paste(src[r:(r + 2L)], collapse = " ")
  miss3 <- c("\"delta\"", "\"gamma_shape\"", "\"gamma_rate\"")[!vapply(
    c("\"delta\"", "\"gamma_shape\"", "\"gamma_rate\""),
    function(p) grepl(p, blk3, fixed = TRUE), logical(1))]
  if (length(miss3))
    stop("[mi-arm] run_all.R's deployed recalibration no longer passes ",
         paste(miss3, collapse = ", "), " to the live suite; re-sync .deploy_delta below.",
         call. = FALSE)
  invisible()
}
.assert_matches_run_all()

# =============================================================================
# 1. LFO-CV — featured model, then the structural baselines
# =============================================================================
LEAKAGE_FREE_LFO <- TRUE     # run_all.R sets this locally, not in 00_config.R
zone_week_outbreak <- dat$zone_week %>%
  dplyr::filter(week_start >= lubridate::floor_date(OUTBREAK_START, "week",
                                                    week_start = get0("WEEK_ANCHOR", ifnotfound = 1L)))

bayes_model <- list()
# The same factory and sampling effort as run_all.R's back-test (600 for brms, TV_FORM_ITER for a
# time-varying form).
bayes_model[[sp$label]] <- make_suite_lfo_model(
  sp, mobility_matrices, gt_pmfs, covariates, osrm_mat, zones_all, linelist = dat$ll,
  trunc_delay = .trunc_asof, iter_brms = 600L, iter_stan = TV_FORM_ITER)

# THE FREQUENTIST BASELINES ARE SCORED THROUGH THE LFO, exactly as run_all.R scores them
# (`c(BASELINE_MODELS, bayes_lfo_models)`). They are cheap deterministic forecasts, and two
# of them — Distance-B1 and Adjacency-B7 — are NAMED IN FIGURE_BASELINE_METHODS, so Figure 2A
# cannot be drawn without them. An earlier revision passed the featured model alone; the
# panel then failed at render time, 110 minutes after the arms had finished, with "2 series
# are absent from detection_curves.csv". The assertion below now catches that at arm time.
# BASELINE_MODELS IS DEFINED IN run_all.R, NOT IN A SOURCED MODULE, so this harness has to
# carry its own copy — which is a drift risk, and therefore guarded: .assert_baselines_match()
# below checks run_all.R's own block still names the same three labels and calls the same
# three forecast_* functions. The bodies themselves come from 05_baseline_models.R, which IS
# sourced, so only the wiring is duplicated here.
BASELINE_MODELS <- Filter(Negate(is.null), list(
  `Gravity-B4` = function(zw, ti, hz, cut) {
    Yw <- zone_week_to_wide(zw, zones_all)
    forecast_B4(Yw, pop_vec, osrm_mat, ti, hz, zones_all) %>% dplyr::mutate(method = "Gravity-B4")
  },
  `Distance-B1` = function(zw, ti, hz, cut) {
    Yw <- zone_week_to_wide(zw, zones_all)
    forecast_B1(Yw, osrm_mat, ti, hz, zones_all, alpha = 1) %>% dplyr::mutate(method = "Distance-B1")
  },
  `Adjacency-B7` = function(zw, ti, hz, cut) {
    Yw <- zone_week_to_wide(zw, zones_all)
    forecast_B7_adjacency(Yw, osrm_mat, ti, hz, zones_all, method = "distance") %>%
      dplyr::mutate(method = "Adjacency-B7")
  }
))
.assert_baselines_match <- function() {
  src <- readLines(file.path(ST_DIR, "run_all.R"), warn = FALSE)
  i <- grep("^BASELINE_MODELS <- Filter", src)
  if (length(i) != 1L) {
    warning("[mi-arm] could not locate run_all.R's BASELINE_MODELS block; the copy below ",
            "is unguarded.", call. = FALSE); return(invisible())
  }
  blk <- paste(src[i:(i + 20L)], collapse = " ")
  need <- c("`Gravity-B4`", "`Distance-B1`", "`Adjacency-B7`",
            "forecast_B4(", "forecast_B1(", "forecast_B7_adjacency(")
  miss <- need[!vapply(need, function(p) grepl(p, blk, fixed = TRUE), logical(1))]
  if (length(miss))
    stop("[mi-arm] run_all.R's BASELINE_MODELS block no longer matches this harness's copy ",
         "(missing: ", paste(miss, collapse = ", "), "). Re-sync before trusting the sweep.",
         call. = FALSE)
  invisible()
}
.assert_baselines_match()
.lfo_models <- c(BASELINE_MODELS, bayes_model)
message(sprintf("[mi-arm] LFO over %d model(s): %s", length(.lfo_models),
                paste(names(.lfo_models), collapse = ", ")))

t0 <- Sys.time()
lfo <- run_invasion_lfo(zone_week_outbreak, .lfo_models,
                        delay = .trunc_asof,
                        horizons = LFO_HORIZONS, analysis_date = ANALYSIS_DATE,
                        min_atrisk_events = 0L,
                        nowcast_fn = apply_nowcast_correction,
                        linelist = if (isTRUE(LEAKAGE_FREE_LFO)) dat$ll else NULL)
lfo_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
if (is.null(lfo)) stop("[mi-arm] run_invasion_lfo() returned NULL for arm ", TAG, call. = FALSE)
message(sprintf("[mi-arm] LFO done in %.1f min: %d rows, %d folds", lfo_min, nrow(lfo),
                dplyr::n_distinct(lfo$fold_id)))

# ---- ORDER MATTERS, AND IT IS run_all.R's ORDER -----------------------------
#   (1) attach the prequential recalibration to the MODEL rows
#   (2) THEN append the structural baselines (which blank the recal columns, because their
#       rows are a copy of a model's rows with only the score replaced and would otherwise
#       inherit that model's prequential factor)
#   (3) THEN evaluate
# Appending before attaching would fit a delta to a rank-only baseline; evaluating before
# appending would leave the baselines unscored, which is the defect run_all.R fixed on
# 2026-09-23. Neither is a detail the pooled figures could survive.
recal_tbl <- NULL
if (isTRUE(INVASION_RECALIBRATE)) {
  lfo <- tryCatch(attach_invasion_recalibration(lfo),
    error = function(e) { warning("[mi-arm] recal attach failed: ", e$message,
                                  " — continuing with raw probabilities only."); lfo })
  recal_tbl <- tryCatch(invasion_delta_table(lfo),
    error = function(e) { warning("[mi-arm] delta table failed: ", e$message); NULL })
}

.bl_specs <- list(
  list(lbl = "Baseline-gravity",
       sc  = function() naive_epicentre_inflow_scores(mobility_matrices[["M4"]],
                          EPICENTRE_ZONES, pop_vec, zones_all)),
  list(lbl = "Baseline-flowminder-inflow",
       sc  = function() {
         .W <- mobility_matrices[["M_cohort"]]
         if (is.null(.W)) {
           warning("[mi-arm] M_cohort unavailable; the Flowminder-inflow baseline falls back ",
                   "to the pooled short-trip kernel M1.", call. = FALSE)
           .W <- mobility_matrices[["M1"]]
         }
         naive_epicentre_inflow_scores(.W, EPICENTRE_ZONES, pop_vec, zones_all)
       }),
  list(lbl = "Baseline-travel-time",
       sc  = function() epicentre_travel_time_scores(osrm_mat, EPICENTRE_ZONES, zones_all)))
.bl_ok <- character(0)
for (.b in .bl_specs) {
  .sc <- tryCatch(.b$sc(), error = function(e) {
    warning(sprintf("[mi-arm] %s scores failed: %s", .b$lbl, conditionMessage(e)),
            call. = FALSE); NULL })
  if (is.null(.sc) || !any(is.finite(.sc) & .sc > 0)) next
  lfo <- append_naive_detection_curve_model(lfo, .sc, .b$lbl)
  if (.b$lbl %in% lfo$method) .bl_ok <- c(.bl_ok, .b$lbl)
}
# A MISSING BASELINE IS FATAL HERE. Figure 2A's claim is that the model beats three
# structural nulls; an arm short of one would be pooled into a panel drawing a mean over a
# different number of arms per series, with nothing on the figure saying so.
if (length(.bl_ok) != 3L)
  stop("[mi-arm] only ", length(.bl_ok), "/3 structural baselines are available (",
       paste(.bl_ok, collapse = ", "), "). Figure 2A needs all three on every arm.",
       call. = FALSE)

# THE ARM'S SHARED SUPPORT EQUALS THE PUBLISHED ONE ONLY IF NO PUBLISHED METHOD LOST CELLS. The
# arm fits the featured model and the comparators, so its support is anchored on those; run_all.R
# anchors on every method it scored. When every published method covers every cell (coverage 1,
# no partial CV) the two coincide; otherwise arm 0 would be scored, recalibrated and deployed on
# a different support from the published run, so it is refused.
local({
  .pe <- tryCatch(readr::read_csv(file.path(ST_DIR, "outputs", "diagnostics", "invasion_evaluation.csv"),
                                  show_col_types = FALSE), error = function(e) NULL)
  if (is.null(.pe) || !all(c("method", "coverage", "partial_cv") %in% names(.pe)))
    stop("[mi-arm] the published invasion_evaluation.csv (with coverage and partial_cv) is ",
         "needed to check the shared support; run run_all.R first.", call. = FALSE)
  .bad <- unique(.pe$method[(.pe$coverage < 1) %in% TRUE | .pe$partial_cv %in% TRUE])
  if (length(.bad))
    stop("[mi-arm] published method(s) ", paste(.bad, collapse = ", "), " did not cover every ",
         "at-risk cell, so the published shared support is narrower than this arm's; arm 0 ",
         "could not reproduce the published numbers.", call. = FALSE)
})
ev <- evaluate_invasion(lfo)

# ---- detection curves, on THIS ARM'S OWN shared support ---------------------
# The support is resolved per arm and NOT intersected across arms: `invasion_common_cells()`
# is a property of one cross-validation, and each arm's evaluation is internally consistent
# with its own labels. Intersecting would score every arm on the cells that survived the
# unluckiest completion — a different, smaller study. 48 pools at the recall level, on the
# shared k axis, and reports the spread in n_atrisk so the reader can see the denominators.
# run_all.R's OWN .dc_methods construction: the featured model, the appended structural
# baselines, and whichever of the LFO-scored comparators are present.
.dc_methods <- unique(c(sp$label, .bl_ok,
                        intersect(c("Gravity-B4", "Distance-B1", "Adjacency-B7",
                                    "Bayes-ens-mean", "Bayes-ens-median"),
                                  unique(lfo$method))))
.dc_methods <- .dc_methods[!is.na(.dc_methods)]
# FIGURE 2A READS FIGURE_BASELINE_METHODS FROM THE CONFIG, so whatever that names must be
# scored here. Checking it inside the arm means a config change that adds a comparator
# fails in six minutes rather than after the whole sweep, at render time.
.fig_bl <- get0("FIGURE_BASELINE_METHODS", ifnotfound = character(0))
.fig_missing <- setdiff(c(FEATURED_METHOD, .fig_bl), .dc_methods)
if (length(.fig_missing))
  stop("[mi-arm] FIGURE_BASELINE_METHODS names comparator(s) this arm does not score: ",
       paste(.fig_missing, collapse = ", "), ". Figure 2A needs every one of them, and it ",
       "reads that constant from 00_config.R. Add them to the LFO model set or to the ",
       "appended structural baselines before running the sweep.", call. = FALSE)
message(sprintf("[mi-arm] detection curves for %d method(s): %s",
                length(.dc_methods), paste(.dc_methods, collapse = ", ")))
detection <- dplyr::bind_rows(lapply(LFO_HORIZONS, function(hh) {
  .cells <- tryCatch(invasion_common_cells(lfo, hh), error = function(e) NULL)
  if (!length(.cells))
    warning(sprintf("[mi-arm] no shared support at h=%s; curves are UNRESTRICTED.", hh),
            call. = FALSE)
  dplyr::bind_rows(lapply(.dc_methods, function(m)
    tryCatch(compute_detection_curve(lfo, m, hh, ks = 1:25,
                                     common_support = FALSE, support_cells = .cells),
             error = function(e) NULL)))
}))
if (!nrow(detection)) stop("[mi-arm] no detection curves were computed.", call. = FALSE)

# =============================================================================
# 2. LIVE forecast for the featured model — WITH the posterior draws
# =============================================================================
# This mirrors fit_bayes_suite()'s single-model path rather than calling it, because the
# suite fits every grid model by MCMC and only one is needed here. The suite fits every model on
# its OWN full design (only its stacking LOO uses the rows common to the stacked models), so
# the featured model's design is all that is built.
.roll_on <- isTRUE(get0("ROLLING_PREDICTORS", ifnotfound = TRUE))
des <- tryCatch(build_invasion_design(zone_week_nc, mobility_matrices, gt_pmfs,
         covariates, osrm_mat, zones_all, mob = sp$mob, gt = sp$gt,
         candidates = unique(c(sp$cov, "log_pop", "ccvi", "d_min")),
         rolling_ll    = if (.roll_on) dat$ll else NULL,
         rolling_delay = if (.roll_on) .trunc_asof else NULL),
         error = function(e) NULL)
if (is.null(des)) stop("[mi-arm] the featured model's design could not be built.", call. = FALSE)
.tv <- sp$tv %||% "none"
fit <- if (.spec_is_stan(sp)) {
  .ex <- invasion_exposure(zone_week_nc, dat$ll, zones_all, mobility_matrices, gt_pmfs,
                           covariates, osrm_mat, sp$mob, sp$gt, intersect(sp$cov, des$feat),
                           des$center, des$scale, .trunc_asof)
  stan_invasion_fit(des, .ex, zones_all, sp$cov, sp$use_rw, sp$use_frailty, iter = 2000L,
                    label = sp$label)
} else if (.tv_is_structural(.tv)) {
  fit_bayes_renewal_tv(des, cov_spec = sp$cov, iter = 2000L,
                       link = sp$link %||% "cloglog", tv = .tv)
} else {
  fit_bayes_renewal(des, cov_spec = sp$cov, iter = 2000L,
                    link = sp$link %||% "cloglog")
}
if (is.null(fit)) stop("[mi-arm] the featured model's deployed fit returned NULL.", call. = FALSE)

rt_draws <- if (max(LFO_HORIZONS) >= 2L)
  bayes_rt_week_draws(dat$ll, gt_pmfs[[sp$gt]],
                      week_start = max(zone_week_nc$week_start),
                      issue_date = ANALYSIS_DATE) else NULL
off <- bayes_forecast_offsets(zone_week_nc, mobility_matrices, gt_pmfs, covariates, osrm_mat,
                              zones_all, sp$mob, sp$gt, sp$cov, LFO_HORIZONS,
                              rt_draws = rt_draws)
if (is.null(off)) stop("[mi-arm] bayes_forecast_offsets() returned NULL.", call. = FALSE)

# THE DEPLOYED RECALIBRATION FACTOR, this arm's own. run_all.R deploys the POOLED delta
# (fitted on all folds, which for a live forecast genuinely is the past) into the live
# forecast; the prequential per-fold delta stays on the scored rows. Both are per arm.
# WITH ITS GAMMA POSTERIOR, exactly the columns run_all.R passes (asserted in
# .assert_matches_run_all()): .suite_delta_for() attaches gamma_shape/gamma_rate and
# .scale_invasion_draws() then draws delta_s ~ Gamma(shape, rate) per posterior draw, as the
# published forecast does. Selecting c("method", "horizon", "delta") alone used to fix delta
# at its point value in every arm, so arm 0 no longer reproduced the published table.
.deploy_delta <- if (isTRUE(INVASION_RECALIBRATE_DEPLOY) && !is.null(recal_tbl) && nrow(recal_tbl))
  .suite_delta_for(recal_tbl[, intersect(c("method", "horizon", "delta", "gamma_shape",
                                           "gamma_rate"), names(recal_tbl))],
                   sp$label, LFO_HORIZONS) else NULL
message(sprintf("[mi-arm] deployed delta: %s%s",
                if (is.null(.deploy_delta)) "none" else
                  paste(sprintf("h%s=%.4f", names(.deploy_delta), .deploy_delta), collapse = " "),
                if (is.null(.deploy_delta)) ""
                else if (!is.null(attr(.deploy_delta, "gamma_shape")))
                  sprintf(" (drawn per posterior draw from Gamma(%.4g, %.4g))",
                          attr(.deploy_delta, "gamma_shape"), attr(.deploy_delta, "gamma_rate"))
                else " (fixed: no Gamma posterior in the delta table)"))

# POSTERIOR DRAWS, NOT A SUMMARY. The MI posterior is the MIXTURE over completions, and with
# draws in hand that mixture is exact — 48 stacks them through the pipeline's own
# .mix_and_summarise_draws(), which resamples JOINT draws so cross-zone correlation, and
# therefore the rank credible interval, survives pooling. Summarising here and pooling the
# summaries would force a normal approximation onto probabilities that live near zero.
#
# delta = NULL: the draws are stored on the RAW hazard scale and the recalibration is
# applied by 48 with .scale_invasion_draws(). Storing them pre-scaled would bake one arm's
# delta into the pooled object and make the raw twin unrecoverable.
# The tv branch mirrors predict_bayes_invasion_tv_scales() EXACTLY — .tv_offsets() on the
# offsets AND .tv_linpred_fn(fit) as the linear predictor. Either half alone is wrong: the
# offsets carry the forecast week's level and the linpred_fn carries the rw1 walk, and the
# featured model being tv = "none" today is not a reason to leave the other path latent and
# broken for whichever model wins selection next.
drl <- if (.spec_is_stan(sp)) {
  stan_invasion_forecast_draws(fit, off, LFO_HORIZONS)       # RAW, as the suite's stan path
} else {
  .off_used <- if (.tv_is_structural(.tv)) .tv_offsets(fit, off, des) else off
  .lp_fn    <- if (.tv_is_structural(.tv)) .tv_linpred_fn(fit) else NULL
  .bayes_invasion_draws(fit, .off_used, des, LFO_HORIZONS, delta = NULL, linpred_fn = .lp_fn)
}
if (is.null(drl)) stop("[mi-arm] .bayes_invasion_draws() returned NULL.", call. = FALSE)

# THE ARM'S OWN SUMMARY IS TAKEN FROM THE FULL POSTERIOR, BEFORE ANY SUBSAMPLING. The
# subsample exists only to make the POOLED mixture equally weighted across completions; there
# is no reason to throw half the draws away from this arm's own risk scores, which feed the
# rank-stability table and the per-arm variability report. Summarising the full posterior
# also means arm 0's risk scores reproduce the published ones exactly, which is a second
# reproduction check for free.
province_map <- load_province_map()
vuln_index   <- tryCatch(compute_vulnerability_index(covariates, zones_all, osrm_mat = osrm_mat),
                         error = function(e) { warning("[mi-arm] vulnerability index: ", e$message); NULL })
preds_arm <- .finalise_invasion_summary(.scale_invasion_draws(drl, .deploy_delta), affected_now)
preds_arm$method <- sp$label
rs_arm <- tryCatch({
  r <- compute_risk_scores(preds_arm, province_map)
  if (!is.null(vuln_index)) r <- add_risk_indices(r, vuln_index)
  r
}, error = function(e) { warning("[mi-arm] risk scores: ", e$message); NULL })

# ---- NOW subsample, for the PERSISTED draws only -----------------------------
# EQUAL DRAWS PER ARM IS A REQUIREMENT OF MI, not a convenience: the mixture weights each
# completion by how many draws it contributes, so an arm that happened to sample more would
# silently count for more. Subsample to exactly MI_DRAWS with a stream derived from this
# arm's seed, and let 48 assert every arm matched.
.n_avail <- nrow(drl[[1]]$cum)
if (.n_avail < MI_DRAWS)
  stop(sprintf(paste0("[mi-arm] the posterior carries %d draws but ONSET_MI_DRAWS=%d. Lower ",
                      "it, or raise iter; padding by resampling would give this arm ",
                      "duplicated draws and understate its within-imputation variance."),
               .n_avail, MI_DRAWS), call. = FALSE)
# on.exit inside local() so the caller's RNG stream is restored even if sample.int throws.
.ix <- local({
  .old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(.old)) assign(".Random.seed", .old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set_fixed_seed(MI_SEED + 977L)   # per-arm subsample of the posterior draws (the onset draw uses no stream)
  sort(sample.int(.n_avail, MI_DRAWS, replace = FALSE))
})
# ONE index across horizons, exactly as .mix_and_summarise_draws() does: a per-horizon
# resample would leave h=1 and h=2 on different posterior draws of the same fit, breaking
# the within-draw ordering p(<=1) <= p(<=2).
drl <- lapply(drl, function(el) {
  el$cum <- el$cum[.ix, , drop = FALSE]
  el$p   <- el$p[.ix, , drop = FALSE]
  el
})
message(sprintf("[mi-arm] draws: %d of %d retained for the mixture (the arm summary above used all %d)",
                MI_DRAWS, .n_avail, .n_avail))

# =============================================================================
# 3. Cascade arm (optional) — tables only, into this arm's sandbox
# =============================================================================
# A SUBPROCESS, NEVER source(). run_cascade.R executes on load and writes a key_outputs
# bundle; sourcing it here would also inherit this session's objects. CASCADE_OUT_ROOT
# sandboxes every product to this arm; CASCADE_RUN_FIGURES=0 stops it drawing Figure 4 into
# the canonical manuscript_figures (40_cascade_next_dominoes.R DELETES the published
# Figure4 when a panel fails, so this is not only a tidiness question).
# ONE THING THE SANDBOX DOES NOT COVER: the cascade reads 16b's short-horizon factor table from
# the canonical outputs/diagnostics (cascade_fit_delta()), i.e. the PUBLISHED run's, for its
# import-term delta_imp and the week-1 transfer check. Both are diagnostics only -- the seeding
# level is refitted per arm on K-week counts -- but an arm's report quotes the published value.
cascade_dir <- file.path(MI_OUT, sprintf("cascade_%s", TAG))
cascade_ok  <- FALSE
if (MI_CASCADE) {
  dir.create(cascade_dir, recursive = TRUE, showWarnings = FALSE)
  # The cascade subprocess gets the SAME two seeds: this arm's onset completion, and the model
  # seed every arm holds.
  .cenv <- c(RANDOM_SEED = as.character(RANDOM_SEED),
             ONSET_IMPUTE_SEED = as.character(MI_SEED),
             ONSET_IMPUTE_STRATUM = "confirmed",
             CASCADE_OUT_ROOT = cascade_dir,
             CASCADE_RUN_FIGURES = "0",
             # The KERNEL is pinned for the same reason the model is (D3): 30_projection_config.R
             # otherwise re-derives it from model_selection.json per arm, and a selection that
             # moved between arms would silently change what the pooled Figure 4 projects with.
             CASCADE_KERNEL = mobility_kernel_from_method(sp$label),
             # ...and the featured model's FORM of the import coefficient, pinned likewise.
             CASCADE_FORM = sp$form %||% "fixed",
             # Validation arms are OFF: the backtest, SBC and kernel-sensitivity sweeps answer
             # questions about the cascade's construction, not about onset uncertainty, and
             # they are ~40% of its runtime. URBAN is ON — Figure 4 panels B and C are built
             # from urban_scenario_summary.csv and urban_impact_*.csv.
             CASCADE_RUN_BACKTEST = "0", CASCADE_RUN_SBC = "0", CASCADE_RUN_SENS = "0",
             CASCADE_RUN_URBAN = "1",
             CASCADE_LAYER_CACHE = file.path(cascade_dir, "layer.rds"))
  message("[mi-arm] cascade -> ", cascade_dir)
  t1 <- Sys.time()
  # stderr MUST be the same PATH as stdout, not the string "2>&1". system2() treats a
  # character stderr as a FILE NAME unless it is identical to stdout, so "2>&1" created a
  # junk file of that literal name in the working directory and the cascade's error output —
  # the only thing that explains a failed arm — never reached the log the message below
  # points the user to.
  .clog <- file.path(MI_OUT, sprintf("cascade_log_%s.txt", TAG))
  .rc <- system2("Rscript", c(file.path(ST_DIR, "run_cascade.R")),
                 env = paste0(names(.cenv), "=", .cenv),
                 stdout = .clog, stderr = .clog)
  cascade_min <- as.numeric(difftime(Sys.time(), t1, units = "mins"))
  .need <- file.path(cascade_dir, "cascade", "tables",
                     c("cascade_reach_scores_all_zones.csv", "urban_scenario_summary.csv"))
  cascade_ok <- .rc == 0L && all(file.exists(.need))
  if (!cascade_ok)
    stop(sprintf(paste0("[mi-arm] the cascade arm failed (exit %d, %.1f min) or did not write ",
                        "its tables: %s. Pooling M-1 arms silently is not an option, so this ",
                        "arm is not written — see %s."),
                 .rc, cascade_min, paste(basename(.need[!file.exists(.need)]), collapse = ", "),
                 .clog), call. = FALSE)
  message(sprintf("[mi-arm] cascade done in %.1f min", cascade_min))
} else {
  cascade_min <- NA_real_
  message("[mi-arm] ONSET_MI_CASCADE=0 — cascade skipped (Figure 4 will not be poolable).")
}

# =============================================================================
# 4. Persist
# =============================================================================
# Draws go in their own file: they are ~8 MB an arm and 48 reads them only for the Figure-3
# mixture, while the variability tables read the (small) summaries.
saveRDS(list(arm = MI_ARM, seed = MI_SEED, onset_impute_seed = as.integer(ONSET_IMPUTE_SEED),
             random_seed = as.integer(RANDOM_SEED), tag = TAG, draws = drl,
             zones = drl[[1]]$zones, n_draws = MI_DRAWS,
             affected_now = affected_now, deploy_delta = .deploy_delta,
             fingerprint = ONSET_FINGERPRINT),
        file.path(MI_OUT, sprintf("draws_%s.rds", TAG)))

saveRDS(list(
  # BOTH SEEDS. `seed` (= onset_impute_seed) is the completion; random_seed is the model seed,
  # which 48 requires to be ONE value across arms while the onset seeds are all distinct.
  arm = MI_ARM, seed = MI_SEED, seed_base = MI_SEED_BASE, tag = TAG,
  onset_impute_seed = as.integer(ONSET_IMPUTE_SEED), random_seed = as.integer(RANDOM_SEED),
  method = sp$label, stratum = as.character(.arm_stratum),
  # The comparator set Figure 2A will ask for. Recorded per arm because 00_config.R is a
  # live file: it was edited DURING the first production sweep, and an arm run before the
  # edit scores a different comparator set from one run after it. 48 asserts all M agree.
  figure_baselines = .fig_bl, dc_methods = .dc_methods,
  fingerprint = ONSET_FINGERPRINT,
  analysis_date = ANALYSIS_DATE, training_cutoff = training_cutoff,
  # BOTH nuisance estimators below are imputation-INVARIANT, for the same structural reason:
  # each is fitted only on records whose onset is actually RECORDED, which are precisely the
  # records this sweep does not impute. `delay_imputation` is the onset->sample delay read from
  # the params file (the nowcast's FALLBACK delay; the onset imputation itself draws donor
  # intervals, not this fit -- the field keeps its name so earlier arms still pool); the as-of
  # truncation uses the date_index of records with a usable recorded onset. 48 asserts both are
  # identical across arms — any variation means an arm resolved a different stratum, a
  # different params file or a different estimator, NOT a different completion.
  delay_imputation = effective_onset_sample_delay(stratum = .arm_stratum)[
    c("family", "mean", "sd", "estimator", "n_fit")],
  trunc_asof = .trunc_asof[c("family", "mean", "sd", "estimator")],
  # The DEPLOYED (extract) truncation, which is NOT imputation-invariant (each arm replays the
  # vintages under its own onset seed): recorded so the data-layer table can show it.
  trunc_deployed = if (is.null(.trunc_deployed)) NULL else
    .trunc_deployed[c("family", "mean", "sd", "estimator")],
  data_definition = nondefault_data_definition(),
  # PROVENANCE 48 CHECKS BEFORE POOLING (beside analysis_date and method above): the line-list
  # snapshot this arm read, and the onset-imputation MECHANISM load_linelist() actually ran —
  # the donor-interval imputation replaced the parametric delay draw on 2026-10-15, and an arm
  # built under the old mechanism is a completion of a different data definition.
  linelist_snapshot = tryCatch({
    .lj <- jsonlite::fromJSON(LINELIST_JSON, simplifyVector = TRUE)
    paste(.lj$folder, .lj$processed_at, sep = "|")
  }, error = function(e) NA_character_),
  onset_impute = list(mode = get0(".ONSET_IMPUTE_MODE", ifnotfound = NA_character_),
                      desc = get0(".ONSET_IMPUTE_DESC", ifnotfound = NA_character_),
                      donor_min = get0("ONSET_DONOR_MIN", ifnotfound = NA),
                      sample_date_source = get0("SAMPLE_DATE_SOURCE", ifnotfound = NA_character_)),
  lfo = lfo, eval = ev, detection = detection, recal = recal_tbl,
  preds = preds_arm, risk_scores = rs_arm,
  first_onset = .first_onset, affected_now = affected_now,
  n_confirmed = sum(dat$ll$confirmed %in% TRUE),
  n_onset_imputed_confirmed = sum(dat$ll$onset_imputed %in% TRUE & dat$ll$confirmed %in% TRUE),
  # The sitrep top-up rows (alert_id "SITREP-CONF-*", 01_data_prep.R) are confirmed and imputed
  # by construction but are a count reconciliation, not case records with a diagnostic sample:
  # recorded apart so the captions can describe each for what it is.
  n_sitrep_confirmed = sum(dat$ll$confirmed %in% TRUE & grepl("^SITREP-CONF-", as.character(dat$ll$alert_id))),
  n_onset_imputed_linelist = sum(dat$ll$onset_imputed %in% TRUE & dat$ll$confirmed %in% TRUE &
                                   !grepl("^SITREP-CONF-", as.character(dat$ll$alert_id))),
  # The deployed recalibration factor. It was previously stored ONLY in the draws file, so
  # 48's data-layer table — which reads the arm file — reported it as NA for every arm while
  # the mixture (which reads the draws file) applied it correctly. A reporting gap, not a
  # scientific one, but a published table of "the delta each completion deployed" that is
  # entirely NA is worse than no table.
  deploy_delta = .deploy_delta,
  cascade_dir = if (cascade_ok) cascade_dir else NA_character_,
  minutes_lfo = lfo_min, minutes_cascade = cascade_min,
  minutes_total = as.numeric(difftime(Sys.time(), t_arm, units = "mins")),
  run_at = Sys.time()
), file.path(MI_OUT, sprintf("arm_%s.rds", TAG)))

print(as.data.frame(ev[, intersect(c("method", "horizon", "auc_pr", "auc_pr_skill",
                                     "mean_rank_of_truth", "log_score", "base_rate",
                                     "n_invasions", "n_atrisk", "n_folds", "partial_cv"),
                                   names(ev))]))
message(sprintf("[mi-arm] arm %s complete in %.1f min -> %s", TAG,
                as.numeric(difftime(Sys.time(), t_arm, units = "mins")),
                file.path(MI_OUT, sprintf("arm_%s.rds", TAG))))

}   # end of the script-run guard

