# BDBV2026 — Spatial Invasion Risk

Public code and daily forecast outputs for the spatial invasion-risk model used in the
2026 Bundibugyo ebolavirus (BDBV) outbreak in the Democratic Republic of the Congo,
stewarded by the Institut National de Recherche Biomédicale (INRB).

**The question the model answers.** For every DRC health zone that has recorded *no*
confirmed BDBV case so far, what is the probability it records its **first** confirmed
case in the next **1 week** (horizon `h=1`) or within the next **2 weeks** (`h=2`)?
Zones already affected cannot be newly invaded, so they are masked (`NA`, never `0`) in
every forecast, ranking and evaluation denominator.

The forecast is driven by human mobility out of the affected cluster, combined with a
Bayesian renewal model of transmission intensity at each source zone. Outputs are
re-issued each day against the latest line-list snapshot.

---

## Repository layout

```
spatiotemporal/     Model code — the pipeline that produces the forecasts
  00_config.R       Configuration, zone spine, environment toggles
  01–04c            Data preparation, epidemiological parameters, mobility, nowcasting
  05–08             Baseline and comparison models (frequentist suite, off by default)
  15–22             Invasion workhorse, evaluation, visualisation, daily re-issue
  21_bayesian_renewal.R   The primary Bayesian invasion model
  30–40             3-month cascade projection layer
  run_all.R         Orchestrates the daily run
  tests/            Unit tests for delay, mobility, parameter and pipeline logic
  METHODS.md        Full methods write-up — every model, equation and assumption

outputs/            Daily forecast outputs, one folder per analysis date
  latest.json       Pointer to the most recent analysis date
  YYYY-MM-DD/       See "Daily outputs" below

docs/
  DATA_DICTIONARY.md      Column-by-column definitions for the output CSVs

tools/
  filter_public_outputs.sh  Allowlist that decides what is published here
```

---

## Daily outputs

Each `outputs/<date>/` folder is one complete daily re-issue, roughly 425 KB.

| File | What it is |
|---|---|
| `bayes_risk_scores_all_zones.csv` | **The headline product.** One row per health zone per horizon: calibrated invasion probability with credible interval, national and provincial relative risk, rank, vulnerability index and preparedness priority. |
| `harmonised_confirmed_cases.csv` | Cumulative confirmed cases per health zone, reconciling the DHIS2 line list with the INSP situation reports. See [where confirmed case counts come from](#where-confirmed-case-counts-come-from). |
| `run_info.json` / `run_info.md` | Provenance for the run — line-list snapshot ID, analysis date, training cutoff, number of zones and confirmed cases, cross-validation fold structure, list of invaded zones. Read this before interpreting any forecast. |
| `model_selection.json` / `model_selection.md` | Which model configuration was featured for this run, and the comparison that selected it. |
| `bayes_parameters.csv` | Fitted model terms on the hazard-ratio scale, with credible intervals and `rhat` convergence diagnostics. |
| `bayes_topk_precision_h*.pdf`<br>`bayes_discrimination_summary_h*.pdf`<br>`bayes_lfo_forecast_vs_outcome_h*.pdf`<br>`bayes_model_performance_figure3_h*.pdf`<br>`bayes_predicted_vs_observed_over_folds_h*.pdf`<br>`bayes_priority_scatter_h*.pdf`<br>`bayes_invasion_uncertainty_h*.pdf` | Evaluation diagnostics for both horizons — how well the model has actually been forecasting, assessed by leave-future-out cross-validation. Published alongside the forecasts so the archive can be judged on its record, not just its predictions. |

Column definitions are in [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md).

### Where confirmed case counts come from

**File:** `outputs/<date>/harmonised_confirmed_cases.csv` — one row per health zone (519),
giving cumulative confirmed cases at that date's training cutoff. Present from 2026-08-06
onward; the two earliest published dates predate the pipeline emitting it.

"Harmonised" means the counts reconcile **two** sources rather than reporting one:

- the **DHIS2 line list**, one row per case, held in INRB's access-controlled repository;
- the **INSP situation reports**, specifically the *cumulative* confirmed-case file. The
  daily `new_confirmed_cases` stream is deliberately not used — it undercounts badly
  against the cumulative record.

**The sitrep acts as a floor, never a ceiling.** For each zone the pipeline appends the
shortfall — `max(0, sitrep_cumulative − line-list_confirmed)` — as additional confirmed
records, so every zone reaches at least its officially confirmed count. Where the line
list already meets or exceeds the sitrep, nothing is added. **No case is ever removed.**
Per zone, counts are taken as the maximum across spelling variants on each date, then a
running maximum over dates (which absorbs the sitrep's occasional revisions and dips),
then differenced into dated whole-case increments.

**Why this matters.** In this model a zone counts as "invaded" only once it has a
*confirmed* case in the line list. The line list lags the sitrep for some zones, so a
zone the sitrep has already confirmed can still appear as suspect-only — and would then
be scored as never invaded. That is not just a presentational problem: the retrospective
cross-validation that selects each day's featured model would be scoring every candidate
against a ground truth missing real invasions.

Two things the file is **not**. It is not nowcast-corrected — these are observed counts,
not counts adjusted upward for reporting delay. And it is not a case count in the sense
of a situation report: it is the modelling pipeline's reconciled view at a specific
training cutoff, which is earlier than the analysis date. For official case figures, cite
INSP and WHO, not this file.

Reconciled counts are computed before modelling, with guards against sitrep zones outside
the 519-zone spine and against increments dated after the analysis date (so a back-dated
re-run cannot leak future confirmations). Every appended record is tagged in the pipeline
with a traceable `SITREP-CONF-<zone>-NN` identifier, and the whole step can be switched
off with `APPEND_SITREP_CONFIRMED = FALSE` for a sitrep-free sensitivity run. The
implementation is `.build_sitrep_confirmed_appends()` in
[`spatiotemporal/01_data_prep.R`](spatiotemporal/01_data_prep.R); the full treatment,
including assumptions and effect on the current snapshot, is §1 of
[`spatiotemporal/METHODS.md`](spatiotemporal/METHODS.md).

### Reading METHODS.md

`spatiotemporal/METHODS.md` is mirrored **verbatim** from INRB's analysis repository, so
it describes the methods as of the run it was last written against. Where it names a
*specific* featured model (e.g. "in the current run this is `Bayes-M10-med`") or quotes
fold counts, fitted effect sizes or stacking weights, read those as illustrative of one
run rather than current — the actual featured model has since changed. The authoritative
record for any given date is that date's own files: `model_selection.json` (which model
was featured and why), `run_info.json` (data vintage, fold structure and the code commit
that produced it) and `bayes_parameters.csv` (fitted terms).

The model formulations, equations and assumptions in METHODS.md are current; only the
run-specific numbers drift. References in it to `PLAN_extensions.md`,
`METHODS_AND_RESULTS.md` and similar working documents point to files kept in INRB's
private analysis repository.

### How to read a forecast

- **Rankings are more reliable than absolute probabilities.** The model is calibrated,
  but ascertainment is imperfect; treat the ordering of zones as the operational signal.
- **`p_case` is a confirmed-case probability**, not an infection probability. A zone can
  be seeded without confirming a case inside the horizon.
- **Always check `run_info`** for the training cutoff. The analysis date and the last day
  of usable data are not the same, because of reporting delay.
- **Diagnostics are cumulative.** Later dates have more cross-validation folds and so a
  more informative evaluation than early ones.

### What is deliberately not published here

The full pipeline emits considerably more than this. Excluded, by design:

- **National map PDFs** (~43 MB/day) — they embed full vector basemaps. Anyone can
  redraw them from `bayes_risk_scores_all_zones.csv` plus public DRC health-zone
  shapefiles.
- **`bayes_pairwise_import_force.csv`** (~48 MB/day) — the full zone × zone importation
  matrix.
- **Cascade / 3-month projection outputs and manuscript figures** — the code is here, but
  the projection layer is scenario-based rather than an operational forecast, and its
  outputs are not part of the daily record.

`tools/filter_public_outputs.sh` is the single source of truth for what is published. It
is a strict allowlist: any output not named in it stays out, so new pipeline products are
never published by accident.

---

## Data

**This repository contains no case data.** The model runs on de-identified line-list and
contact-tracing data held in a private, access-controlled INRB repository, together with
public reference data (health-zone population and boundaries, mobility matrices,
vulnerability covariates). The published outputs here are model results aggregated to
health-zone level; they contain no individual-level information.

The code is published for transparency and methodological scrutiny. It will not run
end-to-end without the private inputs, which are governed by INRB data-sharing
agreements. Enquiries about data access should go to INRB.

## Reproducing the model

R 4.x with the packages listed in `spatiotemporal/00_config.R` (`REQUIRED_PACKAGES`),
plus `brms`/`cmdstanr` with CmdStan, `epinowcast`, `EpiNow2`, and `sf` (which needs GDAL,
GEOS and PROJ). The entry point is:

```r
Rscript spatiotemporal/run_all.R
```

Unit tests, which do not require the sensitive inputs:

```r
Rscript spatiotemporal/tests/run_tests.R
```

## Automation and provenance

Outputs are regenerated daily by the INRB analysis pipeline and mirrored here
automatically. Each `outputs/<date>/` folder corresponds to one completed pipeline run
that passed a completeness gate — degraded runs are never published.

**Code and outputs are published together.** `spatiotemporal/` is mirrored from INRB's
private analysis repository at the exact commit that produced that day's outputs, and
both land in a single commit. The code you see here is therefore the code that generated
the forecasts alongside it — not a later version of it.

Every `run_info.json` names that commit:

```json
"code_version": {
  "analysis_repo": "INRB-UMIE/BDBV2026-Analysis",
  "commit": "80d7801d085560e23b9ee71b1b83d5852398443e",
  "source": "recorded"
}
```

`source` distinguishes two cases:

- **`recorded`** — the commit was captured by the pipeline at run time. Exact.
- **`inferred`** — the run predates that mechanism. The commit is reconstructed from the
  analysis ref pinned in CI when those outputs were committed, and carries a `note`
  saying so. This applies to the twelve dates from 2026-07-31 to 2026-08-27; treat those
  as strong evidence rather than certainty.

Because the analysis repository's `main` branch can be ahead of the pinned commit, the
code here may lag the latest INRB development work. That is deliberate: a public forecast
archive is only reproducible if the published code is the code that ran.

Two allowlist scripts govern everything that reaches this repository —
`tools/filter_public_outputs.sh` for outputs and `tools/sync_model_code.sh` for code.
Both rebuild their target from scratch, so upstream deletions and renames propagate and
nothing stale can survive. A scheduled job re-checks the mirror against the pinned commit
and fails if the two have drifted.

## Licence and citation

Code and outputs are released under the licence in [`LICENSE`](LICENSE). If you use
either in academic work, please cite via [`CITATION.cff`](CITATION.cff).
