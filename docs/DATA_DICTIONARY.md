# Data dictionary

Definitions for the CSV outputs in `outputs/<date>/`. Method sections referenced as
"§n" point to [`spatiotemporal/METHODS.md`](../spatiotemporal/METHODS.md).

---

## `bayes_risk_scores_all_zones.csv`

One row per health zone per forecast horizon. Risk quantities are computed **only** for
zones with zero confirmed cases to date; already-affected zones carry `NA` (never `0`),
so they never dilute a ranking or an evaluation denominator.

Throughout, `μ_i` is the featured model's expected first-case count for zone *i* over the
horizon, and `p_i` its calibrated invasion probability.

> **Schema note.** This dictionary describes the current 30-column schema. Dates
> **2026-07-31 through 2026-08-10** carry an earlier 26-column version that lacks
> `p_median`, `rank_med`, `rank_lo` and `rank_hi` — posterior-median and rank-interval
> reporting was added on 2026-08-13. Read the header row rather than assuming a fixed
> column order; every other column is unchanged in name and meaning across all dates.

### Identification

| Column | Definition |
|---|---|
| `health_zone` | DRC health zone name, harmonised to the pipeline's 519-zone spine. |
| `province` | Province containing the zone. |
| `horizon` | Forecast horizon in weeks: `1` = next week, `2` = within the next two weeks. |
| `was_active_before` | `TRUE` if the zone had already recorded a confirmed case before the forecast date. Risk columns are `NA` for these rows. |
| `method` | Identifier of the model configuration that produced the row, e.g. `Bayes-M17-med` (model M17, medium generation time). Cross-reference `model_selection.json`. |

### Invasion probability

| Column | Definition |
|---|---|
| `p_case_invasion` | Probability the zone records its first **confirmed** case within the horizon: `p_case = 1 − exp(−μ_i)`. Ascertainment-agnostic — this is a confirmed-case probability, not an infection probability (§8.1). |
| `p_case_lo`, `p_case_hi` | Lower and upper bounds of the credible interval on `p_case_invasion`. |
| `p_median` | Posterior median invasion probability. Close to `p_case_invasion`; differs where the posterior is skewed. |
| `mu_forecast` | The expected first-case count `μ_i` itself, before transformation to a probability. |

### Rank

| Column | Definition |
|---|---|
| `rank_med` | Zone's rank on the posterior median probability — `1` is highest risk. **The primary operational signal**; rankings are more robust than absolute probabilities. |
| `rank_lo`, `rank_hi` | Rank bounds implied by the credible interval. A wide `rank_lo`–`rank_hi` span means the zone's position is not well determined. |

### Relative risk

Relative risk expresses a zone's risk against a reference group's mean, which lets a zone
that is modest nationally still surface as the top concern within its province (§8.2–8.5).

| Column | Definition |
|---|---|
| `rr_nat` | Nationwide relative risk: `μ_i` divided by the mean `μ` across all at-risk zones. `1.0` = average. |
| `rr_nat_rank` | Rank on `rr_nat`. |
| `rr01_nat` | Bounded 0–1 national index: `p_i / max(p_i)`. The top-risk zone scores `1`. Used as the risk term in `priority`. |
| `rr_ituri`, `rr_ituri_rank` | Relative risk and rank **within Ituri** — mean taken over at-risk Ituri zones only. `NA` outside Ituri. |
| `rr_nordkivu`, `rr_nordkivu_rank` | As above, within Nord-Kivu. |
| `rr_hautuele`, `rr_hautuele_rank` | As above, within Haut-Uele. |

### Vulnerability and preparedness priority

Priority is a **resource-targeting** lens, not a forecast: it asks where an introduction
would do most damage, not only where one is most likely (§8, "Preparedness priority").

| Column | Definition |
|---|---|
| `surveillance_gap` | `1 − rank(healthsite_density)`, percentile-ranked. Higher = less detection reach. |
| `healthcare_gap` | `1 − rank(healthsite_count / population)`, percentile-ranked. Higher = fewer facilities per capita. |
| `access_gap` | `rank(travel time to the nearest zone with a health facility)`, percentile-ranked. Higher = harder to reach care. Zones with their own facility score `0`. |
| `social_vulnerability` | `rank(CCVI)` — percentile-ranked Climate and Conflict Vulnerability Index. Higher = more vulnerable. |
| `healthcare_travel_min` | The underlying road travel time in minutes (OSRM) to the nearest zone with any health facility. `0` for zones with their own. |
| `V` | Vulnerability and capacity index in [0,1]: the equal-weight mean of the four percentile-ranked pillars above. Higher = worse-off. |
| `priority` | `rr01_nat × V`, rescaled to [0,1]. **Multiplicative** — a zone must be *both* at material invasion risk *and* under-resourced to score high. |
| `priority_rank` | Rank on `priority`, `1` = highest priority. |

---

## `harmonised_confirmed_cases.csv`

| Column | Definition |
|---|---|
| `health_zone` | Health zone name on the same 519-zone spine. |
| `cumulative_confirmed_cases` | Confirmed cases in the zone up to that run's **training cutoff** (not the analysis date — they differ by the reporting delay). `0` for at-risk zones. |

Two distinct things are harmonised here, and it is worth keeping them apart:

- **Zone names** are mapped through an alias table onto the canonical 519-zone spine
  shared by every output file, which makes this the correct join key for the risk scores.
- **Counts** reconcile the DHIS2 line list with the INSP situation reports, with the
  sitrep treated as a floor on each zone's confirmed total. Observed counts only — not
  nowcast-corrected. The README section
  [where confirmed case counts come from](../README.md#where-confirmed-case-counts-come-from)
  explains the rule and why it exists.

By construction `cumulative_confirmed_cases > 0` exactly when the zone is flagged
`was_active_before` in `bayes_risk_scores_all_zones.csv` — both derive from the same
reconciled zone-week table, so the two files cannot disagree about which zones are
affected.

---

## `bayes_parameters.csv`

Fitted terms from the Bayesian renewal model.

| Column | Definition |
|---|---|
| `model` | Model configuration the term belongs to, e.g. `Bayes-M8-med`. |
| `term` | Parameter name (`Intercept`, covariate names). |
| `is_intercept` | `TRUE` for the intercept row. |
| `effect_scale` | Scale on which `hr`, `lo` and `hi` are reported — normally `hazard ratio`. |
| `hr` | Posterior median effect on the stated scale. On the hazard-ratio scale, `>1` raises invasion hazard and `<1` lowers it. |
| `lo`, `hi` | Credible interval bounds for `hr`. An interval spanning `1` indicates no clear directional effect. |
| `p_dir` | Probability of direction — posterior mass on the side of no effect the median falls on. `NA` for intercepts. |
| `rhat` | Convergence diagnostic. Values near `1.00` indicate converged chains; **treat anything above `1.01` with caution**. |

---

## `run_info.json` / `run_info.md`

Provenance for the run. Read this before interpreting any forecast — in particular the
distinction between the analysis date and the training cutoff, which differ because of
reporting delay.

| Field | Definition |
|---|---|
| Line-list snapshot | Identifier and processing timestamp of the case-data snapshot used. |
| Analysis (as-of) date | The date the forecast was issued for. |
| Line-list cutoff | Training week anchor — the last week treated as complete. |
| Training window end | Last day of data included. |
| Outbreak start | First date in the outbreak record. |
| Data date range | Span of onset/sample dates present in the data. |
| Health zones (spine) | Total zones in the geographic spine (519). |
| Health zones invaded | Zones with at least one confirmed case; the remainder are the at-risk denominator. |
| Confirmed cases | Cumulative confirmed cases at the cutoff. |
| Forecast horizons | Horizons issued, in weeks. |
| CV folds / fold cutoffs | Leave-future-out cross-validation structure behind the evaluation diagnostics. More folds accumulate over time, so later runs have stronger evaluation. |
| Primary generation time | Generation-time assumption used by the featured model (`short`, `medium`, `long`). |
| Invaded zones | Names of all zones already affected. |
| `code_version` | The analysis-repo commit that produced this run, and whether it was `recorded` at run time or `inferred` for pre-2026-08-28 dates. See the README's provenance section. |

---

## `model_selection.json` / `model_selection.md`

Records which model configuration was featured for the run and the comparison that
selected it, so a change in forecasts across dates can be attributed to a data update or
to a model switch.
