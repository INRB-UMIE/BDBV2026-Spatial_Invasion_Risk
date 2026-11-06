# Model selection — spatiotemporal

_Generated 2026-11-06T10:46:24+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_06112026 (processed_at 2026-11-06T12:00:00) |
| Predictions computed on snapshot | LINELIST_06112026 (processed_at 2026-11-06T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 13, margin 2) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 13 of 18; h2: 13 of 18 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 12) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 12) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 60%; `Bayes-M14-fill-geo` 23%; `Bayes-M16-fill-geo` 7% — kernel M14-fill 84%, M16-fill 7%, M10-fill 5%, M13c-fill 4% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 51.13 | yes | 29.6 / 21.6 | 41.0 / 47.9 | 0.018 / 0.031 | 33.8 / 27.0 |
| 2 | `Bayes-M10-fill-geo` | 48.49 | yes | 28.5 / 20.0 | 42.7 / 50.2 | 0.018 / 0.032 | 34.8 / 26.6 |
| 3 | `Bayes-M13c-fill-geo` | 40.42 | yes | 23.0 / 17.4 | 33.6 / 40.4 | 0.019 / 0.033 | 32.9 / 26.0 |
| 4 | `Bayes-M14-fill-med` | 45.61 | yes | 25.9 / 19.7 | 47.0 / 56.4 | 0.018 / 0.032 | 33.7 / 26.9 |
| 5 | `Bayes-M17-fill-geo` | 37.43 | yes | 21.3 / 16.2 | 37.2 / 44.8 | 0.019 / 0.034 | 32.5 / 25.5 |
| 6 | `Bayes-M16-fill-geo` | 40.22 | yes | 22.2 / 18.0 | 48.5 / 49.5 | 0.019 / 0.033 | 31.8 / 25.9 |
| 7 | `Bayes-M10-fill-med` | 44.44 | yes | 25.7 / 18.7 | 50.8 / 61.0 | 0.019 / 0.033 | 34.4 / 26.3 |
| 8 | `Bayes-M13c-fill-med` | 36.75 | yes | 21.0 / 15.7 | 41.1 / 51.3 | 0.019 / 0.035 | 32.6 / 25.1 |
| 9 | `Bayes-M16-fill-med` | 41.98 | yes | 23.4 / 18.6 | 54.2 / 57.4 | 0.019 / 0.034 | 32.2 / 26.0 |
| 10 | `Bayes-M4c-geo` | 30.95 | yes | 17.6 / 13.3 | 41.6 / 48.9 | 0.020 / 0.036 | 28.6 / 21.8 |
| 11 | `Bayes-M17-fill-med` | 36.83 | yes | 20.8 / 16.1 | 47.5 / 59.6 | 0.020 / 0.035 | 31.8 / 25.2 |
| 12 | `Bayes-M13-fill-geo` | 27.76 | yes | 15.5 / 12.2 | 42.1 / 50.1 | 0.020 / 0.037 | 29.7 / 23.4 |
| 13 | `Bayes-M4c-med` | 29.50 | yes | 17.0 / 12.5 | 45.7 / 54.9 | 0.021 / 0.036 | 28.6 / 21.2 |
| 14 | `Bayes-M13-fill-med` | 32.92 | yes | 18.2 / 14.7 | 52.2 / 70.0 | 0.021 / 0.038 | 29.9 / 24.0 |
| 15 | `Bayes-M4-geo` | 23.42 | yes | 13.4 / 10.0 | 50.8 / 57.3 | 0.021 / 0.038 | 24.2 / 19.6 |
| 16 | `Bayes-M8-fill-geo` | 22.37 | yes | 13.1 / 9.3 | 51.0 / 57.2 | 0.021 / 0.038 | 27.2 / 20.4 |
| 17 | `Bayes-M4-med` | 27.90 | yes | 16.1 / 11.8 | 60.7 / 77.1 | 0.022 / 0.039 | 26.4 / 19.7 |
| 18 | `Bayes-M8-fill-med` | 27.20 | yes | 15.9 / 11.3 | 59.6 / 76.8 | 0.022 / 0.040 | 28.8 / 21.6 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 48.73 | yes | 27.7 / 21.0 | 38.7 / 43.2 | 0.017 / 0.029 | 34.3 / 27.6 |
| 2 | `Bayes-M14-fill-geo` | 51.13 | yes | 29.6 / 21.6 | 41.0 / 47.9 | 0.018 / 0.031 | 33.8 / 27.0 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
