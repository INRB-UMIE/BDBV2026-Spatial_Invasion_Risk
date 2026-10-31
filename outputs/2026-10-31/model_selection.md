# Model selection — spatiotemporal

_Generated 2026-10-31T10:45:49+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_31102026 (processed_at 2026-10-31T12:00:00) |
| Predictions computed on snapshot | LINELIST_31102026 (processed_at 2026-10-31T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 13, margin 10) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 12 of 17; h2: 12 of 17 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 14) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 12) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 68%; `Bayes-M14-fill-geo` 16%; `Bayes-M10-fill-geo` 6% — kernel M14-fill 85%, M10-fill 6%, M16-fill 6%, M13c-fill 2%, M4c 0% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 52.58 | yes | 30.6 / 22.0 | 40.4 / 50.0 | 0.019 / 0.032 | 40.6 / 29.5 |
| 2 | `Bayes-M10-fill-geo` | 47.82 | yes | 27.5 / 20.3 | 41.7 / 51.4 | 0.019 / 0.033 | 37.7 / 27.9 |
| 3 | `Bayes-M13c-fill-geo` | 39.55 | yes | 22.6 / 17.0 | 34.7 / 41.9 | 0.020 / 0.034 | 38.2 / 27.0 |
| 4 | `Bayes-M14-fill-med` | 49.78 | yes | 27.8 / 22.0 | 47.1 / 58.5 | 0.019 / 0.033 | 40.0 / 29.5 |
| 5 | `Bayes-M16-fill-geo` | 38.98 | yes | 21.7 / 17.3 | 42.3 / 50.9 | 0.020 / 0.034 | 36.4 / 26.8 |
| 6 | `Bayes-M10-fill-med` | 48.14 | yes | 27.4 / 20.7 | 49.6 / 62.6 | 0.020 / 0.033 | 38.9 / 28.1 |
| 7 | `Bayes-M17-fill-geo` | 36.62 | yes | 20.6 / 16.0 | 38.1 / 45.8 | 0.020 / 0.035 | 37.7 / 26.9 |
| 8 | `Bayes-M13c-fill-med` | 38.29 | yes | 22.1 / 16.2 | 42.0 / 51.5 | 0.020 / 0.035 | 37.6 / 26.0 |
| 9 | `Bayes-M16-fill-med` | 44.10 | yes | 25.4 / 18.7 | 47.8 / 59.3 | 0.020 / 0.035 | 37.3 / 27.0 |
| 10 | `Bayes-M4c-geo` | 32.14 | yes | 18.5 / 13.7 | 44.5 / 50.7 | 0.021 / 0.037 | 32.4 / 22.9 |
| 11 | `Bayes-M17-fill-med` | 38.80 | yes | 22.2 / 16.6 | 47.9 / 59.6 | 0.021 / 0.036 | 36.9 / 26.8 |
| 12 | `Bayes-M13-fill-geo` | 27.71 | yes | 15.9 / 11.8 | 43.0 / 51.1 | 0.021 / 0.037 | 35.1 / 24.0 |
| 13 | `Bayes-M4c-med` | 30.19 | yes | 17.5 / 12.7 | 48.1 / 55.1 | 0.022 / 0.037 | 30.7 / 22.1 |
| 14 | `Bayes-M13-fill-med` | 38.31 | yes | 22.4 / 15.9 | 52.4 / 64.3 | 0.022 / 0.038 | 36.0 / 25.7 |
| 15 | `Bayes-M4-geo` | 23.76 | yes | 13.8 / 10.0 | 52.1 / 57.9 | 0.022 / 0.038 | 29.4 / 20.3 |
| 16 | `Bayes-M8-fill-geo` | 22.23 | yes | 13.1 / 9.1 | 50.9 / 57.8 | 0.022 / 0.039 | 32.6 / 22.0 |
| 17 | `Bayes-M8-fill-med` | 31.73 | yes | 19.5 / 12.3 | 59.6 / 71.2 | 0.023 / 0.040 | 34.6 / 23.7 |
| 18 | `Bayes-M4-med` | 28.14 | yes | 16.3 / 11.9 | 61.8 / 71.5 | 0.023 / 0.039 | 30.0 / 20.9 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 50.59 | yes | 28.9 / 21.7 | 38.7 / 47.6 | 0.019 / 0.030 | 39.5 / 29.2 |
| 2 | `Bayes-M14-fill-geo` | 52.58 | yes | 30.6 / 22.0 | 40.4 / 50.0 | 0.019 / 0.032 | 40.6 / 29.5 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
