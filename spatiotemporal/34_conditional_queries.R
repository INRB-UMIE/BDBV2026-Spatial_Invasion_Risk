# =============================================================================
# 34_conditional_queries.R — 3-MONTH CASCADE: conditional reach & attribution
# BDBV 2026 DRC · implements PLAN_3MONTH_INVASION.md §3.7 (the task's core ask:
# conditional probability of reaching a zone AND propagating onward).
#
#   cascade_gateway_screen()     cheap flux-based expected outward seedings per zone.
#   cascade_knockout()           EXACT Delta_j (expected downstream invasions removed
#                                by barring zone j as a source) for a shortlist.
#   cascade_conditional_seed()   designed experiment: force-seed a hub at week 1 and
#                                measure downstream reach (generalised Kisangani).
#   cascade_conditional_obs()    observational P(tau_j<=H | tau_i<=H') from the joint
#                                MC histories (for high-reach hubs).
#   cascade_flux_network()       zone- and province-level seeding flux (corridors).
#   cascade_transmission_tree()  modal seeding source per zone (most-likely tree).
#
# Sourced AFTER 30/31/32. Reuses simulate_cascade()'s bar_sources / force_seed /
# flux machinery, so nothing new is hallucinated.
# =============================================================================

# Expected number of downstream zones seeded per source over the horizon
# (rowSums of the expected flux matrix). A cheap gateway SCREEN; knockouts confirm.
cascade_gateway_screen <- function(sim) {
  out_flux <- rowSums(sim$flux)                       # source -> sum over dest
  tibble::tibble(health_zone = sim$zones_all, out_flux = out_flux,
                 was_active_before = sim$affected0) |>
    dplyr::arrange(dplyr::desc(out_flux))
}

# EXACT gateway impact Delta_j = E[#invaded | j a source] - E[#invaded | j knocked out].
# Restricted to a shortlist (467 full re-runs are infeasible, PLAN §3.7): default =
# the top-N by the flux screen among currently-infected / high-out-flux zones.
cascade_knockout <- function(prep, scenario, base_sim, delta, psi = CASCADE_PSI,
                             n_mc = NULL, top_n = CASCADE_KNOCKOUT_TOPN,
                             candidates = NULL) {
  n_mc <- n_mc %||% base_sim$n_mc
  # PAIRED baseline: re-run with NO bar at the SAME n_mc & seed. Delta_j is unbiased;
  # the shared seed gives common random numbers at the start of each iteration (partial
  # variance reduction — the streams diverge once the seeding paths differ in length).
  base_ref <- simulate_cascade(prep, scenario, n_mc = n_mc, delta = delta, psi = psi)
  base_new <- sum(rowMeans(base_ref$new_by_week))     # E[# new invasions over horizon]
  if (is.null(candidates)) {
    scr <- cascade_gateway_screen(base_sim)
    # candidate sources: initially-infected hubs + highest out-flux zones
    cand <- unique(c(sim_top_sources(base_sim, top_n), head(scr$health_zone, top_n)))
    candidates <- head(cand, top_n)
  }
  rows <- lapply(candidates, function(z) {
    sk <- simulate_cascade(prep, scenario, n_mc = n_mc, delta = delta, psi = psi,
                           bar_sources = z)
    ko_new <- sum(rowMeans(sk$new_by_week))
    tibble::tibble(health_zone = z, base_new = base_new, knockout_new = ko_new,
                   delta_j = base_new - ko_new,
                   pct_reduction = 100 * (base_new - ko_new) / max(base_new, 1e-9))
  })
  dplyr::bind_rows(rows) |> dplyr::arrange(dplyr::desc(delta_j))
}

# Fallback knockout-candidate list = the currently-infected (affected) zones in spine
# order. This is only used when the caller passes no `candidates`; the driver always
# passes an incidence-ranked shortlist (run_cascade.R), so this is a last-resort default
# and is deliberately NOT incidence-ranked (the sim object carries no incidence).
sim_top_sources <- function(sim, n) {
  aff <- sim$zones_all[sim$affected0]
  if (!length(aff)) return(character(0))
  head(aff, n)
}

# Designed conditional experiment (PLAN §3.7-1): force-seed `hub` at week 1 and
# measure downstream reach. The generalisation of the Kisangani "seed the three
# zones" scenario to any hub. Returns the reach table under that intervention.
cascade_conditional_seed <- function(prep, scenario, hub, delta, psi = CASCADE_PSI,
                                     n_mc = NULL, horizons = CASCADE_REPORT_HORIZONS,
                                     force_seed_week = 1L) {
  n_mc <- n_mc %||% CASCADE_N_MC
  sim <- simulate_cascade(prep, scenario, n_mc = n_mc, delta = delta, psi = psi,
                          force_seed = hub, force_seed_week = force_seed_week)
  reach <- cascade_reach_table(sim, horizons = horizons)
  reach$conditioned_on <- paste(hub, collapse = "+")
  list(reach = reach, sim = sim)
}

# Observational conditional reach P(tau_j <= H | tau_i <= H') from a base sim's
# joint first-passage matrix. Reliable only for high-reach hubs i (PLAN §3.7-1
# reliability note); returns NA if the conditioning set is too small.
cascade_conditional_obs <- function(sim, hub, H = max(CASCADE_REPORT_HORIZONS),
                                    Hprime = H, min_iters = 30L) {
  zi <- match(hub, sim$zones_all)
  if (is.na(zi)) return(NULL)
  reached_i <- !is.na(sim$tau[zi, ]) & sim$tau[zi, ] <= Hprime
  if (sum(reached_i) < min_iters)
    return(list(hub = hub, n_cond = sum(reached_i), reliable = FALSE, reach = NULL))
  sub <- sim$tau[, reached_i, drop = FALSE]
  p_cond <- rowMeans(!is.na(sub) & sub <= H)
  tibble::tibble(health_zone = sim$zones_all, p_cond = p_cond,
                 was_active_before = sim$affected0) |>
    (\(d) list(hub = hub, n_cond = sum(reached_i), reliable = TRUE,
               reach = dplyr::arrange(d, dplyr::desc(p_cond))))()
}

# Seeding flux network (corridors). Returns top zone edges and a province x
# province flux matrix (expected seedings), from the simulator's flux matrix.
cascade_flux_network <- function(sim, zone_province = NULL, top_edges = 40L) {
  fl <- sim$flux
  idx <- which(fl > 0, arr.ind = TRUE)
  edges <- tibble::tibble(
    source = sim$zones_all[idx[, 1]], dest = sim$zones_all[idx[, 2]],
    expected_seedings = fl[idx]) |>
    dplyr::arrange(dplyr::desc(expected_seedings))
  prov_mat <- NULL
  if (!is.null(zone_province)) {
    ps <- zone_province[match(sim$zones_all, names(zone_province))]
    ps[is.na(ps)] <- "Unknown"
    edges$source_prov <- ps[match(edges$source, sim$zones_all)]
    edges$dest_prov   <- ps[match(edges$dest, sim$zones_all)]
    prov_mat <- tapply(edges$expected_seedings,
                       list(edges$source_prov, edges$dest_prov), sum)
    prov_mat[is.na(prov_mat)] <- 0
  }
  list(edges = head(edges, top_edges), all_edges = edges, province_matrix = prov_mat)
}

# Most-likely invasion tree: for each seeded zone, its modal seeding source
# (argmax over the flux column). A directed forest rooted at the epicentre.
cascade_transmission_tree <- function(sim) {
  fl <- sim$flux
  dest_has <- which(colSums(fl) > 0)
  rows <- lapply(dest_has, function(j) {
    col <- fl[, j]; src <- which.max(col)
    tibble::tibble(dest = sim$zones_all[j], modal_source = sim$zones_all[src],
                   share = col[src] / sum(col), expected_seedings = sum(col))
  })
  dplyr::bind_rows(rows) |> dplyr::arrange(dplyr::desc(expected_seedings))
}
