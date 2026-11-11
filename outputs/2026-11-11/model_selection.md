# Model selection — spatiotemporal

_Generated 2026-11-11T10:45:16+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_11112026 (processed_at 2026-11-11T12:00:00) |
| Predictions computed on snapshot | LINELIST_11112026 (processed_at 2026-11-11T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 12, margin 9) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 13 of 18; h2: 13 of 18 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 9) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 7) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 64%; `Bayes-M10-fill-geo` 16%; `Bayes-M16-fill-geo` 12% — kernel M14-fill 67%, M10-fill 16%, M16-fill 14%, M13c-fill 3% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 44.89 | yes | 24.3 / 20.6 | 42.8 / 44.5 | 0.019 / 0.032 | 36.7 / 25.8 |
| 2 | `Bayes-M10-fill-geo` | 41.71 | yes | 22.2 / 19.5 | 43.8 / 46.4 | 0.019 / 0.032 | 33.1 / 24.9 |
| 3 | `Bayes-M13c-fill-geo` | 35.63 | yes | 18.9 / 16.7 | 37.7 / 41.0 | 0.020 / 0.034 | 34.7 / 24.5 |
| 4 | `Bayes-M16-fill-geo` | 35.59 | yes | 18.3 / 17.3 | 42.8 / 44.6 | 0.020 / 0.034 | 31.3 / 24.5 |
| 5 | `Bayes-M14-fill-med` | 41.54 | yes | 22.0 / 19.5 | 53.4 / 57.9 | 0.020 / 0.033 | 34.7 / 26.3 |
| 6 | `Bayes-M17-fill-geo` | 33.55 | yes | 17.8 / 15.8 | 40.7 / 44.2 | 0.021 / 0.035 | 34.3 / 23.3 |
| 7 | `Bayes-M10-fill-med` | 39.72 | yes | 21.2 / 18.5 | 56.4 / 62.8 | 0.020 / 0.033 | 32.7 / 24.8 |
| 8 | `Bayes-M16-fill-med` | 38.28 | yes | 20.7 / 17.5 | 51.7 / 55.2 | 0.021 / 0.035 | 33.7 / 25.5 |
| 9 | `Bayes-M13c-fill-med` | 33.79 | yes | 17.8 / 16.0 | 49.5 / 55.2 | 0.021 / 0.035 | 33.7 / 24.6 |
| 10 | `Bayes-M4c-geo` | 27.71 | yes | 14.3 / 13.4 | 44.7 / 48.1 | 0.021 / 0.035 | 26.4 / 19.5 |
| 11 | `Bayes-M13-fill-geo` | 26.24 | yes | 13.6 / 12.6 | 44.4 / 48.6 | 0.021 / 0.036 | 31.0 / 21.8 |
| 12 | `Bayes-M17-fill-med` | 33.95 | yes | 17.7 / 16.2 | 57.3 / 65.2 | 0.021 / 0.035 | 31.9 / 24.7 |
| 13 | `Bayes-M4-geo` | 20.84 | yes | 11.0 / 9.9 | 50.9 / 54.0 | 0.022 / 0.036 | 20.8 / 15.9 |
| 14 | `Bayes-M4c-med` | 25.57 | yes | 13.3 / 12.3 | 51.8 / 57.0 | 0.022 / 0.036 | 24.8 / 18.7 |
| 15 | `Bayes-M8-fill-geo` | 20.01 | yes | 10.5 / 9.5 | 49.9 / 53.3 | 0.022 / 0.036 | 23.4 / 17.3 |
| 16 | `Bayes-M13-fill-med` | 33.32 | yes | 17.8 / 15.5 | 68.4 / 80.3 | 0.022 / 0.037 | 31.3 / 23.8 |
| 17 | `Bayes-M8-fill-med` | 25.05 | **spiky** | 12.9 / 12.1 | 73.1 / 85.3 | 0.023 / 0.038 | 26.8 / 20.0 |
| 18 | `Bayes-M4-med` | 24.73 | **spiky** | 13.2 / 11.5 | 74.6 / 86.4 | 0.022 / 0.037 | 23.4 / 17.4 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 46.31 | yes | 24.9 / 21.4 | 41.1 / 42.6 | 0.018 / 0.029 | 37.3 / 27.0 |
| 2 | `Bayes-M14-fill-geo` | 44.89 | yes | 24.3 / 20.6 | 42.8 / 44.5 | 0.019 / 0.032 | 36.7 / 25.8 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
