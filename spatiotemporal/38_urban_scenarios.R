# =============================================================================
# 38_urban_scenarios.R — CASCADE: urban-hub invasion scenarios (timing-aware)
# BDBV 2026 DRC · implements the "what if a major city is invaded" suite.
#
# Generalises the Kisangani conditional analysis to the country's largest, most
# connected metros. For each hub city and each seeding WEEK, the hub's urban-core
# health zones are FORCE-SEEDED at that week; the cascade then projects the hub's
# onward spread through the mobility network. Impact is measured against the
# baseline (no urban seeding) cascade:
#   - per-zone attributable increase   dP_i = P_cond(i) - P_base(i)
#   - relative risk                    RR_i = P_cond(i) / P_base(i)
#   - threshold crossings              # zones the seeding pushes above P>t (material,
#                                      non-noise upward crossings), and % of baseline
#   - expected additional invasions    E[added] = sum_i dP_i  (i != hub), with a paired
#                                      standard error and interval
#   - provinces newly exposed
# TIMING is a first-class axis: earlier seeding leaves more weeks for onward
# spread, so impact by week 13 is larger for earlier seeds.
#
# TWO THINGS THIS FILE GETS RIGHT THAT AN EARLIER VERSION DID NOT.
#   (1) EVERY CONTRAST IS PAIRED. The conditional run and its baseline share one
#       random-number stream (common random numbers, 32_cascade_simulator.R), so the
#       difference is taken iteration by iteration. Differenced across two independently
#       summarised runs, the flicker of the ~350 zones the seeding cannot reach contributed
#       a summed error the size of the effect itself — which is how the 2026-09-07 run came
#       to report negative expected additions for most cities and to rank zones 1,500 km
#       away as the most affected. Paired, 94% of those zones return exactly zero.
#       "Materially elevated" is correspondingly a TEST — the 90% paired interval excludes
#       zero — not a fixed 0.03 cut.
#   (2) THE SEEDED CITY'S R IS SWEPT, NOT ASSUMED. Left at its province pool, a city with no
#       case history inherits the national value, at which a one-to-two case seed usually
#       fades out — so the headline result was a statement
#       about stochastic fade-out wearing the clothes of a statement about connectivity.
#       The pooled arm remains the primary specification; the sweep makes the premise
#       visible (CASCADE_SEED_R_SWEEP), and applies to the FORCED hub zones only.
#
# Sourced AFTER 30-35 (uses simulate_cascade, cascade_reach_table, and the 35
# shapefile/plot helpers .cascade_shape/.cascade_join/.cascade_ggsave/OUT_CASCADE).
# ASSUMPTIONS documented inline. All quantities are downstream of, and consistent
# with, the confirmed-case-scale cascade (time index = symptom onset).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

# ---- which zones get seeded: an explicit, data-derived rule -----------------
# "[City] invaded" is operationalised as seeding the city's urban-core health zones;
# the cascade then carries spread to the rest of the metro and its catchment. The
# seed set is DERIVED (not hand-listed) in three transparent steps, so a reviewer can
# audit both what was chosen and what was rejected (see urban_hub_selection.csv):
#
#   1. CITY = the union of the health-system `antenne` groups below, from the GRID3/
#      CIESIN health-zone shapefile. This is the only city-like administrative
#      grouping in the pipeline. NOTE: an antenne is a health-system catchment named
#      after its hub city, NOT the city — the Lubumbashi antenne reaches Sakania
#      (20,659 km2) and Pweto; the Mbuji-Mayi antenne reaches Kabeya Kamwanga. Step 2
#      is what removes the hinterland.
#   2. URBAN CORE = zones whose WorldPop population DENSITY clears
#           min(URBAN_CORE_DENS_ABS, URBAN_CORE_DENS_FRAC x the city's max density).
#      The ABSOLUTE arm (1,500 /km2) is the DEGURBA / EU-UN "urban centre" density and
#      does the work for the dense cities, where it lands inside a real cliff
#      (Lubumbashi 1,522 -> 87 /km2; Mbuji-Mayi 3,798 -> 773 /km2). The RELATIVE arm
#      rescues Kananga and Tshikapa, whose health zones bundle the town with its rural
#      hinterland so that no zone reaches 1,500 /km2 and no density cliff exists;
#      taking the LOWER of the two means the relative arm only ever binds there.
#      Relative-only was rejected: the max is set by whichever tiny CBD zone happens to
#      be smallest (Lubumbashi's Kamalondo is 1 km2 at 12,735 /km2), which would push
#      genuinely urban Ruashi/Kenya/Kisanga out of the core.
#      Density (not population count) is essential: DRC health zones span 2-20,659 km2,
#      so ranking by raw population count selects RURAL zones — the most populous zone
#      of the Lubumbashi antenne is Sakania (355k over 20,659 km2, 18/km2, ~300 km out).
#   3. HUB = the top URBAN_HUB_K core zones by Flowminder RELOCATION volume, ranked by
#      URBAN_HUB_RANK_BY. Default "inflow" = "where an importation into this city would
#      most plausibly land", which is the scenario being posed. Inflow is preferred over
#      outflow to avoid ranking seeds by an outbound quantity that parallels the one
#      generating the measured impact; the two rankings agree here in 4 of 5 cities.
#      Population is NOT used to rank: a force-seeded zone receives a zone-independent
#      n_seed_fun() draw and a province-pool R (32_cascade_simulator.R), so the seed's
#      own population never enters the dynamics — only its row/column of W does.
URBAN_CITY_ANTENNES <- list(
  Kinshasa     = c("Kinshasa-Ouest", "Kinshasa-Centre", "Kinshasa-Est"),
  Lubumbashi   = "Lubumbashi",
  `Mbuji-Mayi` = "Mbuji-Mayi",
  Kananga      = "Kananga",
  Tshikapa     = "Tshikapa")
URBAN_CORE_DENS_ABS  <- 1500               # DEGURBA urban-centre density (people/km2)
URBAN_CORE_DENS_FRAC <- 0.25               # relative arm, for cities below the absolute cut
URBAN_HUB_K          <- c(Kinshasa = 3L, Lubumbashi = 3L, `Mbuji-Mayi` = 3L,
                          Kananga = 3L, Tshikapa = 2L)   # seeds per city
URBAN_HUB_RANK_BY    <- Sys.getenv("URBAN_HUB_RANK_BY", unset = "inflow")  # inflow|outflow|pop
# Frozen fallback, used ONLY if the shapefile / WorldPop density / Flowminder inputs are
# unavailable (the derivation is then skipped with a warning). This is the hand-picked
# list the module used before the rule above was introduced; it is kept so the module
# still runs on a partial data checkout, not as the preferred definition.
URBAN_HUBS_FALLBACK <- list(
  Kinshasa     = c("Gombe", "Limete", "Kimbanseke"),
  Lubumbashi   = c("Lubumbashi", "Kampemba", "Ruashi"),
  `Mbuji-Mayi` = c("Diulu", "Bonzola", "Nzaba"),
  Kananga      = c("Kananga", "Katoka", "Lukonga"),
  Tshikapa     = c("Tshikapa", "Kanzala"))

URBAN_SEED_WEEKS <- c(1L, 4L, 8L)          # timing axis (weeks into the projection)
URBAN_THRESHOLDS <- c(0.10, 0.20, 0.50)    # includes the P>50% example metric
# RR is undefined below this baseline probability. The floor is expressed as a MINIMUM NUMBER OF
# MONTE-CARLO ITERATIONS behind the denominator, not as a bare probability: at URBAN_N_MC = 4000
# the old 5e-4 admitted a denominator of FOUR iterations, whose Poisson relative SE is +-50%, and
# `max_rr` then takes the MAXIMUM over 50-90 eligible zones — i.e. it selects the luckiest
# denominator. The shipped table reported Lubumbashi "Max RR 61x" off p_base = 0.00100 = 4/4000.
# 25 iterations gives a relative SE of +-20%, which is defensible for a published ratio.
URBAN_ATTR_MIN_ITER <- 25L
URBAN_ATTR_FLOOR <- 5e-4                   # absolute backstop; the iteration rule dominates
# WHY M IS HIGHER HERE THAN IN THE MAIN CASCADE, and what changed.
#
# This block used to end with "proper fix (not yet implemented): common random numbers".
# It is now implemented (32_cascade_simulator.R), and the numbers below are the before/after
# so a reviewer can see the size of what it fixed rather than take it on trust.
#
# BEFORE. The baseline and the conditional run were separate stochastic runs sharing only a
# master seed, so their streams desynchronised as soon as the seeding paths differed in
# length. The ~350 zones this LOCAL seeding cannot reach therefore flickered: per-zone
# attributable-difference SD ~0.008 at M=2000, ~0.006 at M=4000. Worse, those flickers do
# NOT divide down by sqrt(n) when summed into the headline, because they are positively
# correlated within a run. Measured at M=4000 over 6 replicate master seeds:
#   * baseline vs baseline, INDEPENDENT streams (true effect exactly 0): SE = 1.78 zones.
#   * the then-production structure (shared master seed only):           SE = 0.28-0.42 zones.
# An effect of 0.5 zones sits inside that, which is exactly how the 2026-09-07 grid produced
# negative expected additions (Kananga -0.49, Tshikapa -0.20 at seed week 8) and a Kinshasa
# watch-list topped by zones beside the actual outbreak 1,500 km away. The six Kananga
# week-8 replicates ran +0.04 to +1.22 (mean +0.51): the published -0.49 was a low draw
# around a small POSITIVE effect.
#
# AFTER. Every draw is now an inverse CDF indexed by ZONE, from one stream per iteration, so
# the two runs stay coupled for the whole horizon. Measured on the same far-zone set:
# per-zone SD 0.017 -> 0.0006, and 94% of those zones return EXACTLY zero rather than a small
# non-zero flicker (against 9% before). The headline carries a paired interval from
# cascade_paired_contrast(), and "materially elevated" is that interval excluding zero rather
# than a fixed cut.
#
# M STAYS HIGH ANYWAY. Pairing removes the noise from the DIFFERENCE, not from the levels:
# p_base, p_cond, RR and the threshold counts are still ordinary Monte-Carlo estimates, and
# the credible intervals still need enough parameter draws. Do not lower M on the strength
# of the pairing.
URBAN_N_MC <- as.integer(Sys.getenv("URBAN_N_MC",
                unset = if (identical(Sys.getenv("CASCADE_SMOKE"), "1")) "40" else "4000"))
# FALLBACK ONLY. With a paired contrast, "materially elevated" is a test — the 90% paired
# interval excludes zero (see cascade_paired_contrast) — and this fixed cut is used only
# when no paired contrast is available. It was never a noise floor nor a pre-registered
# effect size, and under pairing the per-zone SE is ~0.0006, so a 0.03 cut would discard
# real effects fifty times larger than the noise.
URBAN_DELTA_MIN <- 0.03

# Which zones a figure should show as materially elevated. Prefers the paired test and
# falls back to the fixed cut, so every panel in the suite draws the same set.
.urban_material <- function(d) {
  if ("elevated" %in% names(d)) return(d$elevated %in% TRUE)
  is.finite(d$delta) & d$delta > URBAN_DELTA_MIN
}

# Validate hub zones against the spine; drop (and warn on) any that are absent.
resolve_urban_hubs <- function(hubs, zones_all) {
  out <- lapply(names(hubs), function(nm) {
    ex <- intersect(hubs[[nm]], zones_all); miss <- setdiff(hubs[[nm]], zones_all)
    if (length(miss)) message("[urban] ", nm, ": not in spine, dropped: ",
                              paste(miss, collapse = ", "))
    ex
  })
  names(out) <- names(hubs)
  out <- out[vapply(out, length, integer(1)) > 0]
  # Carry the derivation provenance through, and record what SURVIVED the spine check
  # as `seeded`, so the audit CSV reflects the zones actually force-seeded (not just
  # the zones the rule picked).
  cand <- attr(hubs, "candidates")
  if (!is.null(cand)) {
    cand$seeded <- mapply(function(city, nom) city %in% names(out) && nom %in% out[[city]],
                          cand$city, cand$Nom)
    attr(out, "candidates") <- cand
    attr(out, "derived") <- attr(hubs, "derived")
  }
  out
}

# Build the per-zone candidate table backing the hub rule: every zone of every city
# antenne, with the three quantities the rule can see (density, population, mobility)
# and the flags that decide selection. Returns NULL if any required input is missing.
urban_hub_candidates <- function(city_antennes = URBAN_CITY_ANTENNES,
                                 dens_abs = URBAN_CORE_DENS_ABS,
                                 dens_frac = URBAN_CORE_DENS_FRAC,
                                 k = URBAN_HUB_K, rank_by = URBAN_HUB_RANK_BY) {
  dens_path <- file.path(WORLDPOP_DIR, "worldpop__pop_density__static.csv")
  pop_path  <- file.path(WORLDPOP_DIR, "worldpop__pop_count__static.csv")
  od_path   <- file.path(FLOWMINDER_DIR,
                         get0("FLOWMINDER_OD_FILE",
                              ifnotfound = "flowminder__outflow__static.matrix.csv"))
  need <- c(SHAPEFILE_PATH, dens_path, pop_path, od_path)
  if (!all(file.exists(need)) || !requireNamespace("sf", quietly = TRUE)) return(NULL)
  shp <- sf::st_drop_geometry(sf::st_read(SHAPEFILE_PATH, quiet = TRUE))
  if (!all(c("Nom", "antenne", "PROVINCE") %in% names(shp))) return(NULL)
  dens <- readr::read_csv(dens_path, show_col_types = FALSE)
  pop  <- readr::read_csv(pop_path,  show_col_types = FALSE)
  names(dens)[2] <- "pop_density"
  # Directed RELOCATION table (Flowminder monthly home-location changes, NOT trips and NOT
  # the model kernel W): rows = origin, cols = destination. It is used ONLY to rank candidate
  # hub zones; the cascade itself propagates with CASCADE_KERNEL. Column sums (arrivals) are
  # the importation proxy used to rank; row sums are the outbound counterpart.
  od <- as.matrix(utils::read.csv(od_path, row.names = 1, check.names = FALSE))
  outflow <- rowSums(od, na.rm = TRUE); inflow <- colSums(od, na.rm = TRUE)

  d <- shp[, c("PROVINCE", "antenne", "Nom")]
  d$city <- NA_character_
  for (nm in names(city_antennes)) d$city[d$antenne %in% city_antennes[[nm]]] <- nm
  d <- d[!is.na(d$city), ]
  if (!nrow(d)) return(NULL)
  d$pop_count   <- pop$pop_count[match(d$Nom, pop$nom)]
  d$pop_density <- dens$pop_density[match(d$Nom, dens$nom)]
  d$inflow      <- unname(inflow[match(d$Nom, names(inflow))])
  d$outflow     <- unname(outflow[match(d$Nom, names(outflow))])
  d$inflow[is.na(d$inflow)]   <- 0        # zone absent from the OD table = no observed flow
  d$outflow[is.na(d$outflow)] <- 0
  # step 2: urban core = density >= min(absolute DEGURBA cut, dens_frac x city max)
  d <- do.call(rbind, lapply(split(d, d$city), function(g) {
    mx <- suppressWarnings(max(g$pop_density, na.rm = TRUE))
    thr <- if (is.finite(mx)) min(dens_abs, dens_frac * mx) else NA_real_
    g$city_max_density <- mx
    g$core_threshold <- thr
    g$urban_core <- is.finite(g$pop_density) & is.finite(thr) & g$pop_density >= thr
    g
  }))
  # step 3: rank the core by the stated mobility criterion; top-k are the seeds
  score <- switch(rank_by, inflow = d$inflow, outflow = d$outflow, pop = d$pop_count,
                  stop("[urban] URBAN_HUB_RANK_BY must be inflow|outflow|pop: ", rank_by))
  d$rank_score <- score; d$rank_by <- rank_by
  d <- do.call(rbind, lapply(split(d, d$city), function(g) {
    g$core_rank <- NA_integer_
    ic <- which(g$urban_core)
    if (length(ic)) g$core_rank[ic[order(-g$rank_score[ic])]] <- seq_along(ic)
    kk <- if (!is.na(k[g$city[1]])) as.integer(k[g$city[1]]) else 3L
    g$selected <- !is.na(g$core_rank) & g$core_rank <= kk
    g[order(!g$urban_core, g$core_rank, -g$rank_score), ]
  }))
  rownames(d) <- NULL
  tibble::as_tibble(d)
}

# Derive URBAN_HUBS from the candidate table. Falls back (with a warning) to the frozen
# hand-picked list if the inputs are unavailable, so a partial checkout still runs.
derive_urban_hubs <- function(...) {
  cand <- tryCatch(urban_hub_candidates(...), error = function(e) {
    message("[urban] hub derivation failed (", conditionMessage(e), ")"); NULL })
  if (is.null(cand) || !any(cand$selected)) {
    warning("[urban] hub inputs unavailable — using frozen URBAN_HUBS_FALLBACK; ",
            "seed zones are NOT data-derived in this run.", call. = FALSE)
    out <- URBAN_HUBS_FALLBACK
    attr(out, "candidates") <- NULL; attr(out, "derived") <- FALSE
    return(out)
  }
  sel <- cand[cand$selected, ]
  out <- lapply(names(URBAN_CITY_ANTENNES), function(nm)
    sel$Nom[sel$city == nm][order(sel$core_rank[sel$city == nm])])
  names(out) <- names(URBAN_CITY_ANTENNES)
  out <- out[vapply(out, length, integer(1)) > 0]
  attr(out, "candidates") <- cand; attr(out, "derived") <- TRUE
  for (nm in names(out))
    message("[urban] ", nm, " seeds (top-", length(out[[nm]]), " of ",
            sum(cand$urban_core & cand$city == nm), " core zones by ",
            URBAN_HUB_RANK_BY, "): ", paste(out[[nm]], collapse = ", "))
  out
}

URBAN_HUBS <- derive_urban_hubs()

# Run one urban scenario: force-seed hub_zones at seed_week; return reach + sim.
run_urban_scenario <- function(prep, scenario, hub_zones, seed_week, delta, psi, n_mc,
                               seed_r = NA_real_, seed = CASCADE_SEED, pop_vec = NULL) {
  # `seed` is passed EXPLICITLY and is the same value the baseline used: that shared stream
  # is what makes the contrast paired. seed_r applies to the forced hub zones only (see
  # seed_r_scope in simulate_cascade), so the sweep asks about this city's transmissibility
  # rather than silently re-parameterising every zone the front later reaches.
  # pop_vec is passed so that these contrasts run the SAME model as the production scenario
  # tables, which have always supplied it. NOTE: within-zone susceptible depletion was REMOVED
  # (2026-09-19), so pop_vec no longer arms anything — it is accepted and unused. It is still
  # passed because a contrast whose two sides, or whose scenario table and urban table, call
  # simulate_cascade() with different arguments is not a contrast.
  # (This comment previously asserted and denied depletion in one breath, a botched edit that
  # left a reader unable to tell which model had been calibrated.)
  sim <- simulate_cascade(prep, scenario, n_mc = n_mc, delta = delta, psi = psi,
                          force_seed = hub_zones, force_seed_week = seed_week,
                          seed_r = seed_r, seed = seed, pop_vec = pop_vec)
  list(sim = sim, reach = cascade_reach_table(sim, horizons = CASCADE_REPORT_HORIZONS))
}

# Impact of a scenario vs baseline, at one horizon. The seeded hub zones and any
# already-affected (masked) zones are EXCLUDED from the impact summaries (the hub
# is invaded by assumption, not prediction), so metrics reflect DOWNSTREAM effect.
urban_impact_metrics <- function(reach_base, reach_cond, hub_zones, province_map,
                                 horizon = max(CASCADE_REPORT_HORIZONS),
                                 thresholds = URBAN_THRESHOLDS,
                                 # n_mc is an explicit formal because attr_floor's default is a
                                 # function of it: the RR denominator must be backed by at least
                                 # URBAN_ATTR_MIN_ITER iterations, not by a bare probability.
                                 n_mc = get0("URBAN_N_MC", ifnotfound = 4000L),
                                 attr_floor = max(URBAN_ATTR_FLOOR,
                                                  URBAN_ATTR_MIN_ITER / max(n_mc, 1L)),
                                 delta_min = URBAN_DELTA_MIN,
                                 sim_base = NULL, sim_cond = NULL) {
  # Carry the reach table's 90% CREDIBLE interval (p_lo/p_hi) through to the per-zone
  # impact, so downstream figures can show uncertainty on the conditional (and baseline)
  # invasion probabilities. Falls back gracefully if a reach table lacks the CI columns.
  gb <- function(rt) {
    keep <- intersect(c("health_zone", "p_case_invasion", "p_lo", "p_hi"), names(rt))
    rt[rt$horizon == horizon, keep]
  }
  b <- gb(reach_base)
  names(b) <- sub("^p_case_invasion$", "p_base", names(b))
  names(b) <- sub("^p_lo$", "p_base_lo", names(b)); names(b) <- sub("^p_hi$", "p_base_hi", names(b))
  c2 <- gb(reach_cond)
  names(c2) <- sub("^p_case_invasion$", "p_cond", names(c2))
  names(c2) <- sub("^p_lo$", "p_cond_lo", names(c2)); names(c2) <- sub("^p_hi$", "p_cond_hi", names(c2))
  m <- merge(b, c2, by = "health_zone")
  m <- m[!(m$health_zone %in% hub_zones) & is.finite(m$p_base) & is.finite(m$p_cond), ]
  m$delta <- m$p_cond - m$p_base
  # rr is a ratio of LEVELS, not a paired quantity, and is reported only where the baseline
  # clears attr_floor — so it never divides two numbers that are both noise.
  m$rr    <- ifelse(m$p_base > attr_floor, m$p_cond / m$p_base, NA_real_)
  if (!is.null(province_map))
    m$province <- province_map$province[match(m$health_zone, province_map$nom)]

  # PAIRED contrast where both simulations are available.
  #
  # WHAT PAIRING DOES AND DOES NOT DO, stated precisely because it is easy to overclaim.
  # With balanced draw groups the paired POINT estimate is arithmetically identical to
  # p_cond - p_base: the mean of the group means is the overall mean. Pairing changes what
  # the two runs have in COMMON. Var(A - B) = Var(A) + Var(B) - 2Cov(A, B), and sharing one
  # random-number stream drives that covariance almost to the variance itself, so the
  # difference's variance collapses — measured, 25-30x on the zones the seeding cannot
  # reach, with 94% of them returning EXACTLY zero instead of flickering. So the gain is not
  # a better number; it is an interpretable one, with a standard error and a test for which
  # zones moved at all. Recomputing delta here rather than reusing p_cond - p_base is
  # therefore not a correction but a guarantee that delta, its SE and its interval are all
  # built from the same group decomposition.
  # CASES, the burden question the invasion metric cannot answer. It is computed over ALL
  # zones and split in-city vs elsewhere: the invasion metric excludes the seeded hub (invaded
  # by assumption) and the already-infected zones (cannot be invaded again), and BOTH
  # exclusions are wrong for cases — most of the burden lands inside the seeded city, and an
  # already-infected zone can still receive extra cases through the import force. Measured on
  # this outbreak, the in-city burden outweighs the downstream burden by one to two orders of
  # magnitude at every arm of the seeded-R sweep — so the zone-invasion headline and the case
  # headline are answering genuinely different questions. Current figures live in
  # urban_scenario_summary.csv and the report's case table rather than in this comment, which
  # would drift; note only that the ratio is large and that the in-city column carries the
  # assumed seed cases as well as the epidemic they start.
  #
  # `horizon` deliberately does not enter: cases accumulate over the whole projection window,
  # so this is the 13-week total by construction, whereas the invasion contrast is evaluated
  # at a named horizon.
  cas <- NULL
  if (!is.null(sim_base) && !is.null(sim_cond) && !is.null(sim_base$cases_zone))
    cas <- cascade_paired_cases(sim_base, sim_cond, in_city = hub_zones)

  ctr <- NULL
  if (!is.null(sim_base) && !is.null(sim_cond)) {
    ctr <- cascade_paired_contrast(sim_base, sim_cond, horizon = horizon, exclude = hub_zones)
    j <- match(m$health_zone, ctr$per_zone$health_zone)
    m$delta       <- ctr$per_zone$delta[j]
    m$delta_se    <- ctr$per_zone$delta_se[j]
    m$delta_lo    <- ctr$per_zone$delta_lo[j]
    m$delta_hi    <- ctr$per_zone$delta_hi[j]
    m$exact_zero  <- ctr$per_zone$exact_zero[j]
    m$elevated    <- ctr$per_zone$elevated[j]
    # Carried so a reader can SEE which zones the saturation term pushed significantly
    # negative, rather than having to re-derive it from the interval columns.
    m$reduced     <- ctr$per_zone$reduced[j]
  }
  if (!is.null(cas)) {
    k <- match(m$health_zone, cas$per_zone$health_zone)
    m$cases_base  <- cas$per_zone$cases_base[k]
    m$cases_cond  <- cas$per_zone$cases_cond[k]
    m$cases_delta <- cas$per_zone$delta[k]
    m$cases_lo    <- cas$per_zone$delta_lo[k]
    m$cases_hi    <- cas$per_zone$delta_hi[k]
  }
  # Materiality is a TEST when the pairing supports one — a zone counts as elevated if its
  # 90% paired interval excludes zero — and only falls back to the fixed delta_min cut when
  # no paired contrast is available. The fixed cut was neither a noise floor nor a
  # pre-registered effect size; under pairing it discards real small effects.
  material <- if (!is.null(ctr)) (m$elevated %in% TRUE)
              else is.finite(m$delta) & m$delta > delta_min
  thr <- dplyr::bind_rows(lapply(thresholds, function(t) {
    # zones pushed from <=t to >t by the seeding, with a material (non-noise) change;
    # far zones unaffected by this local seeding are excluded from the count.
    newly <- sum(m$p_base <= t & m$p_cond > t & material, na.rm = TRUE)
    nb    <- sum(m$p_base > t, na.rm = TRUE)
    tibble::tibble(threshold = t, n_base = nb, newly_above = newly,
                   pct_increase = 100 * newly / max(nb, 1))
  }))
  # provinces newly exposed: a province with a materially-elevated zone above 0.20 that
  # had no baseline zone above 0.20.
  prov_new <- NA_integer_
  if (!is.null(province_map)) {
    pb <- unique(m$province[m$p_base > 0.20])
    pc <- unique(m$province[m$p_cond > 0.20 & material])
    prov_new <- length(setdiff(pc[!is.na(pc)], pb[!is.na(pb)]))
  }
  # THE WITHIN-CITY SIGN INVERSION — a HISTORICAL note, and no longer expected to appear.
  #
  # While psi was fitted (it reached 11.25 on one frame) some of the seeded city's OWN
  # remaining zones came out significantly NEGATIVE: on the 2026-09-01 run, 8 of Kinshasa's
  # zones, led by Kokolo at -0.022 [-0.026, -0.018]. The mechanism was frontier saturation,
  # sat = exp(-psi * max(f_inv - f_inv0, 0)), where f_inv is the share of a zone's INBOUND
  # mobility coming from invaded zones: seeding three Kinshasa zones sharply raised f_inv for
  # their closest mobility neighbours — which are the other Kinshasa zones — and damped their
  # hazard. Mechanistically that is backwards; `sat` was a phenomenological term fitted to stop
  # the cascade over-predicting the multi-week count, not a transmission mechanism, and at a
  # psi that large it outweighed the direct effect locally.
  #
  # psi is now FIXED AT 0 (30_projection_config.R), so sat == 1 and the damping channel does
  # not exist. Seeding a city can therefore only ADD import force, and a within-city negative
  # can no longer be attributed to saturation. If material negatives still appear here, they
  # are Monte-Carlo noise, the frontier-covariate posterior tail, or the re-seeding channel —
  # NOT protection, and not this artefact. Treat their reappearance as a signal to check the
  # pairing rather than as something to explain away.
  #
  # Structural expectation, stated precisely rather than as a blanket impossibility claim.
  # With psi = 0 the coupling is exactly monotone — an extra source can only raise the
  # import force and flip a seeding indicator from 0 to 1 — so no zone's difference can be
  # negative. With psi > 0 frontier saturation is a genuine NEGATIVE feedback: seeding a
  # city raises its neighbours' invaded fraction and can lower a third zone's hazard, so
  # small negatives are a model property, not an error, and must not be clipped away. What
  # is still suspect is a negative TOTAL whose paired interval covers zero: that is noise.
  .exp_added <- if (!is.null(ctr)) ctr$total else sum(m$delta, na.rm = TRUE)
  .exp_se    <- if (!is.null(ctr)) ctr$se else NA_real_
  .exp_lo    <- if (!is.null(ctr)) ctr$lo else NA_real_
  .exp_hi    <- if (!is.null(ctr)) ctr$hi else NA_real_
  if (is.finite(.exp_added) && .exp_added < 0) {
    if (!is.null(ctr) && is.finite(.exp_hi) && .exp_hi > 0)
      message(sprintf(paste0("[urban] expected_added = %.2f [%.2f, %.2f] at horizon %d: negative ",
                             "but its paired interval covers zero, i.e. no detectable downstream ",
                             "effect at M=%d, not a negative one."),
                      .exp_added, .exp_lo, .exp_hi, horizon, URBAN_N_MC))
    else
      warning(sprintf(paste0("[urban] expected_added = %.2f [%.2f, %.2f] at horizon %d is ",
                             "NEGATIVE with an interval excluding zero. At psi > 0 frontier ",
                             "saturation can genuinely produce this; check psi before reading ",
                             "it as an error."),
                      .exp_added, .exp_lo, .exp_hi, horizon), call. = FALSE)
  }

  list(per_zone = m[order(-m$delta), ],
       thresholds = thr,
       expected_added = .exp_added,                          # E[extra zones invaded by H]
       expected_added_se = .exp_se, expected_added_lo = .exp_lo, expected_added_hi = .exp_hi,
       paired = !is.null(ctr) && isTRUE(ctr$paired),
       contrast = ctr, cases = cas,
       # A case count is uninterpretable without its denominator: "+10 cases" means one thing
       # against a 200-case baseline and another against 2000. The baseline projection over
       # the same window is carried alongside — MEAN and MEDIAN both, because the per-
       # iteration distribution is right-skewed (15.7k vs 13.2k on this frame) and quoting
       # either alone would misrepresent it. The percentage in the report divides mean by
       # mean, so the two estimators are never mixed.
       cases_baseline_total = if (is.null(sim_base)) NA_real_
                              else sum(rowMeans(sim_base$cases_zone)),
       cases_baseline_median = if (is.null(sim_base)) NA_real_
                               else unname(stats::median(colSums(sim_base$cases_zone))),
       cases_baseline_pred_lo = if (is.null(sim_base)) NA_real_
                                else unname(stats::quantile(colSums(sim_base$cases_zone), 0.05)),
       cases_baseline_pred_hi = if (is.null(sim_base)) NA_real_
                                else unname(stats::quantile(colSums(sim_base$cases_zone), 0.95)),
       n_elevated = sum(material, na.rm = TRUE),             # # zones materially elevated
       # max over an EMPTY set (no materially-elevated zone) is -Inf, which then leaks
       # into the summary CSV and any downstream max()/plot. The relative risk is
       # undefined here, not negatively infinite — return NA.
       max_rr = local({
         v <- m$rr[material]; v <- v[is.finite(v)]
         if (length(v)) max(v) else NA_real_
       }),
       provinces_newly_exposed = prov_new,
       horizon = horizon)
}

# Run the full hub x seed-week grid. Returns a per-scenario summary table plus,
# for the EARLIEST seed week (max impact), the per-zone impact + sim for detailed
# figures. reach_base = the baseline (no urban seeding) reach table.
urban_scenario_grid <- function(prep, scenario, hubs, seed_weeks, delta, psi, n_mc,
                                reach_base, province_map, sim_base = NULL,
                                seed = CASCADE_SEED, pop_vec = NULL,
                                seed_r_sweep = get0("CASCADE_SEED_R_SWEEP",
                                                    ifnotfound = NA_real_)) {
  summ <- list(); detail <- list()
  # THE SEEDED CITY'S REPRODUCTION NUMBER IS A PREMISE, NOT A FINDING. Left at the province
  # pool, a capital with no case history inherits the national value, at which a one-to-two
  # case seed usually fades out — so the analysis
  # answers "what if a capital is seeded and then behaves like the current national average",
  # and quite correctly finds little happens. Sweeping it makes the premise visible and
  # reportable. NA is the pooled arm and is the primary specification.
  sr_sweep <- unique(seed_r_sweep)
  if (!length(sr_sweep)) sr_sweep <- NA_real_
  if (!any(is.na(sr_sweep))) sr_sweep <- c(NA_real_, sr_sweep)   # keep the pooled arm
  for (city in names(hubs)) {
    message("[urban]   scenario: ", city, " (", length(seed_weeks), " seed-weeks, ",
            length(sr_sweep), " seeded-R arms at week ", min(seed_weeks), " -> ",
            length(seed_weeks) + length(sr_sweep) - 1L, " runs, M=", n_mc, ")")
    # The seeded-R sweep runs at the PRIMARY seed week only. Crossing it with the timing
    # axis would quadruple a 13-week, M=4000 grid to answer a question nobody asks: the
    # sweep exists to expose the premise behind the primary specification, and the timing
    # profile is read at the pooled arm. Every other week therefore runs the pooled arm
    # alone, which is exactly what urban_fig_timing() and the report's timing table use.
    sw_primary <- min(seed_weeks)
    # `==` rather than identical(): seed_weeks may arrive as integer or double depending on
    # the caller, and identical(1L, 1) is FALSE — which would silently run the sweep at no
    # week at all and leave the seeded-R panel empty.
    for (sw in seed_weeks) for (sr in (if (isTRUE(sw == sw_primary)) sr_sweep else NA_real_)) {
      sc <- run_urban_scenario(prep, scenario, hubs[[city]], sw, delta, psi, n_mc,
                               seed_r = sr, seed = seed, pop_vec = pop_vec)
      # Pass this grid's own n_mc so the RR denominator floor tracks the actual simulation size
      # (CASCADE_SMOKE drops it to 40, where a 25-iteration floor is most of the sample).
      im <- urban_impact_metrics(reach_base, sc$reach, hubs[[city]], province_map,
                                 n_mc = n_mc, sim_base = sim_base, sim_cond = sc$sim)
      t50 <- im$thresholds[im$thresholds$threshold == 0.50, ]
      t20 <- im$thresholds[im$thresholds$threshold == 0.20, ]
      summ[[paste(city, sw, sr)]] <- tibble::tibble(
        city = city, seed_week = sw, seed_r = as.numeric(sr),
        seed_r_label = if (is.na(sr)) "province pool" else sprintf("R = %.1f", sr),
        n_hub_zones = length(hubs[[city]]),
        expected_added = im$expected_added, expected_added_se = im$expected_added_se,
        expected_added_lo = im$expected_added_lo, expected_added_hi = im$expected_added_hi,
        paired = im$paired,
        # CASES: the burden split. `in_city` is dominated by the assumed 1-2 seed cases at the
        # pooled R (most introductions fade) and by the onward epidemic at a swept R, which is
        # why the two are reported separately rather than summed.
        cases_in_city    = im$cases$in_city$estimate   %||% NA_real_,
        cases_in_city_lo = im$cases$in_city$lo         %||% NA_real_,
        cases_in_city_hi = im$cases$in_city$hi         %||% NA_real_,
        cases_elsewhere    = im$cases$elsewhere$estimate %||% NA_real_,
        cases_elsewhere_lo = im$cases$elsewhere$lo       %||% NA_real_,
        cases_elsewhere_hi = im$cases$elsewhere$hi       %||% NA_real_,
        cases_total    = im$cases$total$estimate %||% NA_real_,
        cases_total_lo = im$cases$total$lo       %||% NA_real_,
        cases_total_hi = im$cases$total$hi       %||% NA_real_,
        # Median and 90% PREDICTIVE band across iterations — the spread of outcomes, as
        # distinct from the Monte-Carlo precision of the mean above.
        cases_total_median  = im$cases$total$median  %||% NA_real_,
        cases_total_pred_lo = im$cases$total$pred_lo %||% NA_real_,
        cases_total_pred_hi = im$cases$total$pred_hi %||% NA_real_,
        cases_in_city_median   = im$cases$in_city$median   %||% NA_real_,
        cases_elsewhere_median = im$cases$elsewhere$median %||% NA_real_,
        cases_baseline_total   = im$cases_baseline_total   %||% NA_real_,
        cases_baseline_median  = im$cases_baseline_median  %||% NA_real_,
        cases_baseline_pred_lo = im$cases_baseline_pred_lo %||% NA_real_,
        cases_baseline_pred_hi = im$cases_baseline_pred_hi %||% NA_real_,
        n_elevated = im$n_elevated,
        newly_above_20 = t20$newly_above, n_base_gt50 = t50$n_base,
        newly_above_50 = t50$newly_above, pct_increase_gt50 = t50$pct_increase,
        max_rr = im$max_rr, provinces_newly_exposed = im$provinces_newly_exposed)
      # Detail figures use the PRIMARY specification: earliest seed week, pooled R.
      if (sw == min(seed_weeks) && is.na(sr)) {
        im$per_zone$city <- city; im$per_zone$seed_week <- sw
        detail[[city]] <- list(impact = im, reach = sc$reach)
      }
    }
  }
  # Pairing diagnostics measured on THIS run, so the report can state what the coupling
  # actually achieved rather than quote a figure from a past run that will drift out of date.
  # "Far" = a zone the seeding never reached in either arm under any iteration is not a
  # useful definition (that is what we are measuring), so far-ness is defined structurally:
  # zones outside the seeded city's province, which is where the spurious pre-CRN watch-list
  # entries came from.
  .city1 <- if (length(detail)) names(detail)[1] else NA_character_
  pz <- tryCatch(detail[[1]]$impact$per_zone, error = function(e) NULL)
  diag <- NULL
  if (!is.null(pz) && "exact_zero" %in% names(pz) && !is.na(.city1)) {
    # The hub's province comes from its SEEDED ZONES, not from the city name: "Mbuji-Mayi"
    # and "Tshikapa" are cities, and matching them against the health-zone spine would miss
    # and silently label the city's own neighbours as far zones. The hub zones themselves are
    # excluded from `pz`, so the province is looked up in the province map instead.
    hub_prov <- if (!is.null(province_map))
      unique(province_map$province[match(hubs[[.city1]], province_map$nom)]) else NA_character_
    hub_prov <- hub_prov[!is.na(hub_prov)]
    far <- if ("province" %in% names(pz) && length(hub_prov)) !(pz$province %in% hub_prov)
           else rep(TRUE, nrow(pz))
    far[is.na(far)] <- TRUE
    diag <- list(city = .city1, n_zones = nrow(pz),
                 # EXPOSED so consumers can restrict "within-city" counts to the hub's own
                 # province. 36_report.R counted `reduced` over EVERY eligible zone while calling
                 # the result "the seeded city's own remaining health zones" — on the shipped run
                 # that put Mungindu (Kwilu) inside Kinshasa.
                 hub_prov = hub_prov,
                 n_far = sum(far),
                 pct_exact_zero_far = 100 * mean(pz$exact_zero[far], na.rm = TRUE),
                 pct_exact_zero_all = 100 * mean(pz$exact_zero, na.rm = TRUE),
                 se_far = stats::median(pz$delta_se[far], na.rm = TRUE),
                 se_far_max = suppressWarnings(max(pz$delta_se[far], na.rm = TRUE)),
                 n_negative = sum(pz$delta < -1e-12, na.rm = TRUE))
  }
  list(summary = dplyr::bind_rows(summ), detail = detail,
       seed_weeks = seed_weeks, min_week = min(seed_weeks), n_mc = n_mc,
       seed_r_sweep = sr_sweep, pairing_diag = diag)
}

# ---- visualisations --------------------------------------------------------
# 1. Attributable-risk choropleth: where the hub's onward spread concentrates.
urban_map_attributable <- function(impact, hub_label, seed_week,
                                   file_prefix = "urban_attributable") {
  d <- impact$per_zone
  # Show only MATERIALLY-elevated zones (paired interval above zero); grey out the rest,
  # which carry only Monte-Carlo flicker (this local seeding does not reach them). This
  # matches the impact metric and removes visual noise so the real catchment stands out.
  d$delta_show <- ifelse(.urban_material(d), d$delta, NA_real_)
  sp <- .cascade_join(d, "delta_show"); if (is.null(sp)) return(invisible(NULL))
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = "grey80", linewidth = 0.04) +
    scale_fill_viridis_c(option = "rocket", direction = -1, na.value = "grey93",
                         name = "Attributable\nincrease in\nP(invasion)", limits = c(0, NA)) +
    labs(title = sprintf("BDBV 2026 — downstream invasion risk attributable to invasion of %s",
                         hub_label),
         subtitle = sprintf("Materially-elevated zones only: increase in 13-week reach vs baseline if %s is seeded at week %d (projection; upper bound).",
                            hub_label, seed_week)) +
    theme_void(base_size = 11) + theme(plot.subtitle = element_text(size = 7.5))
  .cascade_ggsave(p, sprintf("%s_%s", file_prefix, gsub("[^A-Za-z0-9]+", "_", hub_label)))
}

# 2. Top downstream zones by attributable increase (with relative risk).
urban_fig_top_zones <- function(impact, hub_label, top_n = 15L,
                                file_prefix = "urban_top_zones") {
  d <- impact$per_zone[is.finite(impact$per_zone$delta) &
                       .urban_material(impact$per_zone), ]          # material zones only
  if (!nrow(d)) return(invisible(NULL))
  d <- head(d[order(-d$delta), ], top_n)
  d$health_zone <- factor(d$health_zone, levels = rev(d$health_zone))
  p <- ggplot(d, aes(delta, health_zone)) +
    geom_col(fill = "#b2182b") +
    geom_text(aes(label = ifelse(is.finite(rr), sprintf("RR %.1f", rr), "")),
              hjust = -0.1, size = 3, colour = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.18))) +
    labs(title = sprintf("BDBV 2026 — zones most affected by invasion of %s", hub_label),
         subtitle = "Attributable increase in 13-week reach probability; RR = conditional / baseline.",
         x = "attributable increase in P(invasion)", y = NULL) +
    theme_minimal(base_size = 11)
  .cascade_ggsave(p, sprintf("%s_%s", file_prefix, gsub("[^A-Za-z0-9]+", "_", hub_label)),
                  w = 8, h = 6)
}

# 3. Cross-hub summary bar at the earliest seed week (expected additional invasions).
urban_fig_summary_bar <- function(summary, min_week, file = "urban_impact_summary") {
  # The PRIMARY specification: earliest seed week, seeded city at its province pool. The
  # seeded-R arms are a separate panel, because mixing them here would silently average a
  # premise the reader has not been shown.
  d <- summary[summary$seed_week == min_week, ]
  if ("seed_r" %in% names(d)) d <- d[is.na(d$seed_r), ]
  if (!nrow(d)) return(invisible(NULL))
  d$city <- factor(d$city, levels = d$city[order(d$expected_added)])
  has_ci <- "expected_added_lo" %in% names(d) && any(is.finite(d$expected_added_lo))
  p <- ggplot(d, aes(expected_added, city)) +
    geom_col(fill = "#2166ac") +
    {if (has_ci) geom_linerange(aes(xmin = expected_added_lo, xmax = expected_added_hi),
                                colour = "grey25", linewidth = 0.5)} +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d zones elevated; +%d newly >20%%",
                                  n_elevated, newly_above_20)),
              hjust = -0.02, size = 2.9, colour = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
    labs(title = sprintf("BDBV 2026 — downstream impact if a major city is invaded (seeded week %d)", min_week),
         subtitle = paste("Expected additional health zones invaded within 13 weeks, excluding the",
                          "seeded city itself.\nPaired (common random numbers) contrast against the",
                          "same-seed baseline; bars are 90% intervals across posterior draws.",
                          "\nThe seeded city transmits at its province pool — see the seeded-R panel."),
         x = "expected additional zones invaded (13 weeks)", y = NULL) +
    theme_minimal(base_size = 11)
  .cascade_ggsave(p, file, w = 9, h = 5.5)
}

# 3b. The premise made visible: impact against the SEEDED CITY'S reproduction number.
urban_fig_seed_r <- function(summary, min_week, file = "urban_seed_r_sweep") {
  if (!"seed_r" %in% names(summary)) return(invisible(NULL))
  d <- summary[summary$seed_week == min_week, ]
  if (!length(unique(d$seed_r)) || all(is.na(d$seed_r))) return(invisible(NULL))
  dp <- d[is.na(d$seed_r), ]          # the pooled arm, drawn as a reference line
  ds <- d[!is.na(d$seed_r), ]         # the swept arms
  if (!nrow(ds)) return(invisible(NULL))
  # Ribbon first so the fitted line and points sit on top of their own uncertainty band.
  p <- ggplot(ds, aes(seed_r, expected_added, colour = city, group = city)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    {if ("expected_added_lo" %in% names(ds))
       geom_ribbon(aes(ymin = expected_added_lo, ymax = expected_added_hi, fill = city),
                   alpha = 0.12, colour = NA)} +
    {if (nrow(dp)) geom_hline(data = dp, aes(yintercept = expected_added, colour = city),
                              linetype = "dashed", linewidth = 0.5)} +
    geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
    labs(title = "BDBV 2026 — the seeded city's transmissibility is the premise, not the finding",
         subtitle = paste("Expected additional zones invaded by 13 weeks against the reproduction",
                          "number assumed for the seeded city.\nDashed lines: the province-pool",
                          "arm (the primary specification), in which the city inherits the current",
                          "national value.\nAt the pooled value a one-to-two case seed usually",
                          "fades out, which is why that arm shows little effect."),
         x = "reproduction number of the seeded city", y = "expected additional zones (by 13 wk)",
         colour = NULL, fill = NULL) +
    theme_minimal(base_size = 11)
  .cascade_ggsave(p, file, w = 9, h = 6)
}

# 3c. CASE BURDEN: the question the zone-invasion metric cannot answer.
#
# Two panels because the two quantities differ by an order of magnitude and by kind. Cases
# IN the seeded city are mostly the assumed seed plus its local outbreak — large, and driven
# almost entirely by the assumed R. Cases ELSEWHERE are the genuine downstream burden —
# small, and the quantity a national planner is actually trading off. Plotting them on one
# axis would render the second invisible, which is exactly the mistake this panel exists to
# avoid.
urban_fig_cases <- function(summary, min_week, file = "urban_case_burden") {
  if (!"cases_total" %in% names(summary) || all(is.na(summary$cases_total)))
    return(invisible(NULL))
  d <- summary[summary$seed_week == min_week, ]
  if (!nrow(d)) return(invisible(NULL))
  d$city <- as.character(d$city)
  d$arm <- ifelse(is.na(d$seed_r), "province pool", sprintf("R = %.1f", d$seed_r))
  arm_lv <- c("province pool", sprintf("R = %.1f", sort(unique(d$seed_r[!is.na(d$seed_r)]))))
  d$arm <- factor(d$arm, levels = arm_lv)
  # City order by the PRIMARY specification (pooled arm), so this panel, the seeded-R panel
  # and the invasion figures all rank the cities the same way. Guarded: if the pooled arm is
  # absent, factor(levels = character(0)) would turn EVERY city into NA and silently empty
  # the plot, so fall back to the order the summary arrived in.
  pooled <- is.na(d$seed_r)
  ord <- if (any(pooled)) d$city[pooled][order(d$cases_total[pooled])] else unique(d$city)
  d$city <- factor(d$city, levels = ord)
  if (all(is.na(d$city))) return(invisible(NULL))

  # TWO PLOTS, NOT TWO FACETS, because the panels need DIFFERENT X SCALES — and facets can
  # only share a transform. In-city burden spans more than an order of magnitude between the
  # pooled arm and the highest swept R; on a linear axis that crushes the pooled arm (the PRIMARY
  # specification) against the axis while the assumption-driven R = 2.5 arm fills the panel.
  # A log axis there keeps every arm readable. The downstream panel must stay LINEAR: its
  # values straddle zero, which a log scale cannot represent, and zero is the reference that
  # matters for "is there a detectable effect at all".
  base <- function(dat, xlab) ggplot(dat, aes(est, city, colour = arm)) +
    geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.62),
                   linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.62), size = 1.9) +
    labs(x = xlab, y = NULL, colour = "seeded-city R") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())

  din <- data.frame(city = d$city, arm = d$arm, est = d$cases_in_city,
                    lo = d$cases_in_city_lo, hi = d$cases_in_city_hi)
  dout <- data.frame(city = d$city, arm = d$arm, est = d$cases_elsewhere,
                     lo = d$cases_elsewhere_lo, hi = d$cases_elsewhere_hi)
  # log10 needs strictly positive values. In-city estimates are normally positive by
  # construction (the seed itself is at least one case), but they are NOT guaranteed: if
  # every hub zone is already infected at t0 the intervention is null and the column is all
  # zeros. min() over an empty set returns Inf with a warning and would produce a blank
  # panel, so fall back to a linear axis instead of failing quietly.
  pos <- din$est[is.finite(din$est) & din$est > 0]
  use_log <- length(pos) > 0 && all(din$est[is.finite(din$est)] > 0)
  if (use_log) din$lo <- pmax(din$lo, 0.5 * min(pos))
  p_in <- base(din, if (use_log) "additional confirmed cases (log scale)"
                    else "additional confirmed cases") +
    ggtitle("Inside the seeded city")
  if (use_log)
    p_in <- p_in + scale_x_log10(labels = scales::label_number(accuracy = 1))
  else
    p_in <- p_in + geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35)
  p_out <- base(dout, "additional confirmed cases") +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
    ggtitle("Elsewhere in the country")

  p <- (p_in | p_out) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("BDBV 2026 — additional confirmed cases attributable to seeding a city (seeded week %d)", min_week),
      subtitle = paste0(
        "Paired contrast against the same-seed baseline. Bars are 90% intervals on the MEAN ",
        "across posterior draws, not\npredictive intervals, which are far wider (see the ",
        "report's case table). LEFT is on a LOG axis: it spans more than an\norder of magnitude because it is ",
        "driven by the assumed seeded-city R, and it is the burden the zone-invasion metric ",
        "excludes\nby construction. RIGHT is the downstream national burden, on a linear axis ",
        "so that zero stays visible.\nModelled CONFIRMED cases over 13 weeks, not infections: ",
        "this suite estimates no ascertainment fraction and divides by none."),
      theme = ggplot2::theme(legend.position = "bottom")) &
    ggplot2::theme(legend.position = "bottom")
  .cascade_ggsave(p, file, w = 12, h = 6.4)
}

# 4. Timing sensitivity: impact vs seeding week, per hub.
urban_fig_timing <- function(summary, file = "urban_timing_sensitivity") {
  if ("seed_r" %in% names(summary)) summary <- summary[is.na(summary$seed_r), ]
  if (!nrow(summary)) return(invisible(NULL))
  p <- ggplot(summary, aes(seed_week, expected_added, colour = city, group = city)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    {if ("expected_added_lo" %in% names(summary))
       geom_ribbon(aes(ymin = expected_added_lo, ymax = expected_added_hi, fill = city),
                   alpha = 0.12, colour = NA)} +
    geom_line(linewidth = 1) + geom_point(size = 2.4) +
    scale_x_continuous(breaks = sort(unique(summary$seed_week))) +
    labs(title = "BDBV 2026 — earlier invasion of a city drives more onward spread",
         subtitle = "Expected additional zones invaded by 13 weeks, by the week the city is seeded. Later seeding leaves fewer weeks for spread.",
         x = "week the city is invasion-seeded", y = "expected additional zones (by 13 wk)",
         colour = "City") +
    theme_minimal(base_size = 11)
  .cascade_ggsave(p, file, w = 9, h = 5.5)
}

# Provenance for the seed-zone rule: EVERY candidate zone of every city antenne, with
# the quantities the rule sees and the flags that decided it, so a reviewer can audit
# what was rejected as well as what was chosen. Written next to the scenario outputs.
urban_write_hub_selection <- function(dir, hubs = URBAN_HUBS,
                                      file = "urban_hub_selection.csv") {
  cand <- attr(hubs, "candidates")
  if (is.null(cand)) {
    message("[urban] no derivation provenance to write (frozen fallback in use)")
    return(invisible(NULL))
  }
  path <- file.path(dir, file)
  keep <- intersect(c("city", "PROVINCE", "antenne", "Nom", "pop_count",
                      "pop_density", "city_max_density", "core_threshold", "urban_core",
                      "inflow", "outflow", "rank_by", "rank_score",
                      "core_rank", "selected", "seeded"), names(cand))
  readr::write_csv(cand[, keep], path)
  message("  saved ", basename(path))
  invisible(path)
}

# Write per-hub per-zone impact CSVs + the cross-scenario summary table.
urban_write_tables <- function(grid, dir) {
  readr::write_csv(grid$summary, file.path(dir, "urban_scenario_summary.csv"))
  urban_write_hub_selection(dir)
  for (city in names(grid$detail)) {
    tag <- gsub("[^A-Za-z0-9]+", "_", city)
    d <- grid$detail[[city]]$impact$per_zone
    readr::write_csv(d, file.path(dir, sprintf("urban_impact_%s.csv", tag)))
    # CASES go to their OWN file rather than into urban_impact_*.csv, because that frame
    # excludes the seeded hub zones (invaded by assumption) and the zones already infected at
    # t0 — correct for an invasion metric, but it would silently drop the zones carrying most
    # of the case burden. This frame covers every zone and flags `in_city`, so the
    # in-city / elsewhere split is auditable per zone and not only as a pair of totals.
    cz <- grid$detail[[city]]$impact$cases$per_zone
    if (!is.null(cz) && NROW(cz))
      readr::write_csv(cz, file.path(dir, sprintf("urban_cases_%s.csv", tag)))
  }
  invisible(NULL)
}

# ---- publication figures (Figure4-style) -----------------------------------
# Province colour identity (extends the Figure-4 eastern palette to the western/
# southern provinces that host the urban catchments), Okabe-Ito CVD-safe.
.URBAN_OKABE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
.URBAN_PROV_COL <- c(
  "Ituri"="#0072B2","Nord-Kivu"="#D55E00","Haut-Uele"="#009E73","Sud-Kivu"="#CC79A7",
  "Kinshasa"="#E69F00","Haut-Katanga"="#56B4E9","Kasaï-Oriental"="#CC79A7",
  "Kasaï-Central"="#009E73","Kasaï"="#D55E00","Kwilu"="#0072B2",
  "Kwango"="#F0E442","Lomami"="#000000","Mai-Ndombe"="#7A7A7A")
.urban_prov_palette <- function(provs) {
  provs <- unique(as.character(provs[!is.na(provs)]))
  cols <- setNames(character(length(provs)), provs)
  spare <- setdiff(.URBAN_OKABE, unname(.URBAN_PROV_COL)); si <- 1L
  for (p in provs) {
    if (p %in% names(.URBAN_PROV_COL)) cols[p] <- unname(.URBAN_PROV_COL[p])
    else { cols[p] <- spare[((si - 1L) %% max(length(spare), 1L)) + 1L]; si <- si + 1L }
  }
  cols
}
.URBAN_KEYFIG_DIR <- function() file.path(OUT_DIR, "key_outputs", "figures")

# Figure-4 analogue for a city scenario: Panel A = attributable-risk choropleth
# zoomed to the catchment (seeded core in grey), Panel B = baseline -> conditional
# reach dumbbell for the top affected zones, coloured by province. Built from the
# per-zone impact only (regenerable without re-simulating).
urban_figure4 <- function(impact, hub_label, hub_zones, top_n = 15L,
                          file_prefix = "Figure_urban") {
  if (!requireNamespace("patchwork", quietly = TRUE)) return(invisible(NULL))
  sp0 <- .cascade_shape(); if (is.null(sp0)) return(invisible(NULL))
  d <- impact$per_zone
  mat <- d[.urban_material(d), ]
  if (!nrow(mat)) return(invisible(NULL))
  # --- Panel A: attributable map, zoomed to catchment + seeded core ---
  d$dshow <- ifelse(.urban_material(d), d$delta, NA_real_)
  spm <- .cascade_join(d, "dshow"); if (is.null(spm)) return(invisible(NULL))
  hubkey <- tolower(trimws(hub_zones)); matkey <- tolower(trimws(mat$health_zone))
  focus <- spm[spm$.key %in% c(matkey, hubkey), ]
  bb <- sf::st_bbox(focus)
  mx <- 0.18 * (bb[["xmax"]] - bb[["xmin"]]); my <- 0.18 * (bb[["ymax"]] - bb[["ymin"]])
  pA <- ggplot(spm) +
    geom_sf(aes(fill = .val), colour = "grey78", linewidth = 0.06) +
    geom_sf(data = spm[spm$.key %in% hubkey, ], fill = "grey35",
            colour = "white", linewidth = 0.15) +
    scale_fill_viridis_c(option = "rocket", direction = -1, na.value = "grey92",
                         name = "Attributable increase\nin P(invasion)", limits = c(0, NA)) +
    coord_sf(xlim = c(bb[["xmin"]] - mx, bb[["xmax"]] + mx),
             ylim = c(bb[["ymin"]] - my, bb[["ymax"]] + my), expand = FALSE) +
    theme_void(base_size = 11) +
    theme(legend.position = "right", plot.margin = margin(2, 2, 2, 2))
  # --- Panel B: baseline -> conditional dumbbell, top zones by attributable increase ---
  db <- head(mat[order(-mat$delta), ], top_n)
  db$zone <- factor(db$health_zone, levels = rev(db$health_zone))
  pcols <- .urban_prov_palette(db$province)
  lev <- levels(db$zone); ylab_cols <- unname(pcols[as.character(db$province)[match(lev, db$health_zone)]])
  pB <- ggplot(db) +
    geom_segment(aes(x = p_base, xend = p_cond, y = zone, yend = zone),
                 colour = "grey78", linewidth = 1) +
    geom_point(aes(p_base, zone), colour = "grey55", size = 1.9) +
    geom_point(aes(p_cond, zone, colour = province), size = 2.8) +
    scale_colour_manual(values = pcols, name = "Province") +
    scale_x_continuous(labels = scales::percent_format(1), limits = c(0, NA),
                       expand = expansion(mult = c(0.01, 0.06))) +
    labs(x = "P(invasion by 13 wk): baseline (grey) to conditional (colour)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(colour = ylab_cols, face = "bold", size = 9),
          legend.position = "top")
  # No in-plot title/subtitle (caption lives in the manuscript text); keep A/B tags.
  fig <- (patchwork::wrap_elements(full = pA) | patchwork::wrap_elements(full = pB)) +
    patchwork::plot_layout(widths = c(0.55, 0.45)) +
    patchwork::plot_annotation(tag_levels = "A")
  path <- file.path(.URBAN_KEYFIG_DIR(), sprintf("%s_%s", file_prefix,
                    gsub("[^A-Za-z0-9]+", "_", hub_label)))
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !.fk(path)) return(invisible(path))
  ggsave(paste0(path, ".pdf"), fig, width = 12.5, height = 6.2, bg = "white")
  ggsave(paste0(path, ".png"), fig, width = 12.5, height = 6.2, dpi = 150, bg = "white")
  message("  saved ", basename(path)); invisible(path)
}

# Multi-city small-multiples: attributable-risk maps for all cities in one figure
# (national extent), a comprehensive spatial overview. impacts = named list of
# per_zone data frames (one per city).
urban_map_facet <- function(impacts, file = "Figure_urban_all_cities") {
  sp0 <- .cascade_shape(); if (is.null(sp0)) return(invisible(NULL))
  frames <- lapply(names(impacts), function(city) {
    d <- impacts[[city]]
    d$dshow <- ifelse(.urban_material(d), d$delta, NA_real_)
    s <- .cascade_join(d, "dshow"); if (is.null(s)) return(NULL)
    s$city <- city; s
  })
  frames <- frames[!vapply(frames, is.null, logical(1))]
  if (!length(frames)) return(invisible(NULL))
  sp <- do.call(rbind, frames)
  p <- ggplot(sp) + geom_sf(aes(fill = .val), colour = NA) +
    facet_wrap(~city, nrow = 1) +
    scale_fill_viridis_c(option = "rocket", direction = -1, na.value = "grey93",
                         name = "Attributable increase\nin P(invasion)", limits = c(0, NA)) +
    # No in-plot title/subtitle (caption lives in the manuscript text).
    theme_void(base_size = 11) + theme(strip.text = element_text(face = "bold"))
  path <- file.path(.URBAN_KEYFIG_DIR(), file)
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  .fk <- get0("figure_is_kept", ifnotfound = NULL)
  if (is.function(.fk) && !.fk(path)) return(invisible(path))
  ggsave(paste0(path, ".pdf"), p, width = 15, height = 4.6, bg = "white")
  ggsave(paste0(path, ".png"), p, width = 15, height = 4.6, dpi = 150, bg = "white")
  message("  saved ", basename(path)); invisible(path)
}
