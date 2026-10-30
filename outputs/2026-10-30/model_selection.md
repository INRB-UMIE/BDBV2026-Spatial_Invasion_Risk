# Model selection — spatiotemporal

_Generated 2026-10-30T13:31:05+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_30102026 (processed_at 2026-10-30T12:00:00) |
| Predictions computed on snapshot | LINELIST_30102026 (processed_at 2026-10-30T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M13c-fill-geo` (kernel M13c-fill, form fixed; stage 1 `Bayes-M13c-fill-geo`: composite 16, margin 2) — **differs from the current pick** |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 12 of 17; h2: 12 of 17 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 15, margin 11) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 15, margin 9) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 54%; `Bayes-M14-fill-geo` 28%; `Bayes-M13c-fill-geo` 7% — kernel M14-fill 82%, M16-fill 8%, M13c-fill 7%, M10-fill 3% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 50.18 | yes | 28.1 / 22.1 | 43.5 / 50.7 | 0.019 / 0.032 | 31.5 / 27.1 |
| 2 | `Bayes-M10-fill-geo` | 47.36 | yes | 27.6 / 19.8 | 44.9 / 53.1 | 0.019 / 0.033 | 32.7 / 26.1 |
| 3 | `Bayes-M13c-fill-geo` | 39.30 | yes | 21.9 / 17.4 | 34.8 / 41.3 | 0.020 / 0.035 | 30.7 / 25.3 |
| 4 | `Bayes-M14-fill-med` | 44.63 | yes | 24.8 / 19.9 | 49.3 / 58.2 | 0.019 / 0.033 | 31.5 / 26.3 |
| 5 | `Bayes-M17-fill-geo` | 36.11 | yes | 20.1 / 16.0 | 38.5 / 45.6 | 0.020 / 0.036 | 29.8 / 24.8 |
| 6 | `Bayes-M13c-fill-med` | 35.91 | yes | 20.1 / 15.8 | 42.5 / 50.8 | 0.020 / 0.036 | 30.3 / 24.5 |
| 7 | `Bayes-M16-fill-geo` | 39.24 | yes | 21.4 / 17.9 | 51.2 / 52.3 | 0.020 / 0.034 | 30.1 / 25.2 |
| 8 | `Bayes-M10-fill-med` | 43.29 | yes | 24.4 / 18.9 | 53.1 / 62.7 | 0.020 / 0.034 | 32.0 / 25.8 |
| 9 | `Bayes-M4c-geo` | 29.83 | yes | 16.6 / 13.2 | 43.0 / 50.2 | 0.021 / 0.037 | 26.3 / 21.3 |
| 10 | `Bayes-M16-fill-med` | 41.22 | yes | 22.7 / 18.5 | 57.1 / 59.8 | 0.020 / 0.036 | 30.3 / 25.5 |
| 11 | `Bayes-M17-fill-med` | 35.92 | yes | 19.8 / 16.1 | 49.2 / 59.0 | 0.021 / 0.037 | 29.5 / 24.6 |
| 12 | `Bayes-M13-fill-geo` | 26.59 | yes | 14.5 / 12.1 | 42.7 / 50.2 | 0.022 / 0.038 | 27.9 / 22.8 |
| 13 | `Bayes-M4c-med` | 28.75 | yes | 16.2 / 12.5 | 46.7 / 54.3 | 0.022 / 0.038 | 26.5 / 20.7 |
| 14 | `Bayes-M13-fill-med` | 31.12 | yes | 16.4 / 14.7 | 53.4 / 64.0 | 0.022 / 0.040 | 27.3 / 23.4 |
| 15 | `Bayes-M4-geo` | 23.24 | yes | 13.1 / 10.1 | 50.7 / 57.3 | 0.023 / 0.040 | 22.8 / 19.4 |
| 16 | `Bayes-M8-fill-geo` | 21.56 | yes | 12.3 / 9.3 | 50.9 / 57.3 | 0.023 / 0.040 | 25.2 / 19.8 |
| 17 | `Bayes-M4-med` | 26.87 | yes | 15.1 / 11.8 | 61.1 / 71.3 | 0.023 / 0.041 | 24.5 / 19.3 |
| 18 | `Bayes-M8-fill-med` | 26.23 | yes | 14.9 / 11.3 | 60.4 / 71.1 | 0.023 / 0.041 | 26.5 / 21.1 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 47.44 | yes | 26.8 / 20.6 | 41.3 / 46.6 | 0.018 / 0.031 | 32.6 / 27.0 |
| 2 | `Bayes-M14-fill-geo` | 50.18 | yes | 28.1 / 22.1 | 43.5 / 50.7 | 0.019 / 0.032 | 31.5 / 27.1 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
