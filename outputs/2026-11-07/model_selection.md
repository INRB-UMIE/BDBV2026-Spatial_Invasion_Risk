# Model selection — spatiotemporal

_Generated 2026-11-07T10:53:55+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_07112026 (processed_at 2026-11-07T12:00:00) |
| Predictions computed on snapshot | LINELIST_07112026 (processed_at 2026-11-07T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 14, margin 12) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 13 of 18; h2: 13 of 18 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 10) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 8) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 64%; `Bayes-M14-fill-geo` 24%; `Bayes-M16-fill-geo` 5% — kernel M14-fill 88%, M16-fill 6%, M10-fill 4%, M13c-fill 1%, M4c 0% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 54.58 | yes | 32.3 / 22.3 | 40.2 / 47.1 | 0.018 / 0.031 | 41.4 / 29.7 |
| 2 | `Bayes-M10-fill-geo` | 50.66 | yes | 29.8 / 20.8 | 41.9 / 49.2 | 0.018 / 0.031 | 39.0 / 28.3 |
| 3 | `Bayes-M13c-fill-geo` | 41.71 | yes | 24.2 / 17.5 | 34.9 / 41.1 | 0.019 / 0.033 | 39.5 / 27.3 |
| 4 | `Bayes-M14-fill-med` | 50.97 | yes | 28.8 / 22.1 | 47.2 / 56.7 | 0.018 / 0.032 | 40.6 / 29.5 |
| 5 | `Bayes-M16-fill-geo` | 41.28 | yes | 23.5 / 17.8 | 42.2 / 48.7 | 0.019 / 0.033 | 37.7 / 27.3 |
| 6 | `Bayes-M17-fill-geo` | 38.97 | yes | 22.6 / 16.4 | 38.1 / 45.2 | 0.019 / 0.034 | 38.5 / 27.1 |
| 7 | `Bayes-M10-fill-med` | 49.44 | yes | 28.6 / 20.8 | 49.9 / 60.9 | 0.019 / 0.032 | 39.8 / 28.2 |
| 8 | `Bayes-M16-fill-med` | 46.08 | yes | 26.3 / 19.8 | 47.7 / 56.8 | 0.019 / 0.034 | 37.8 / 27.8 |
| 9 | `Bayes-M13c-fill-med` | 39.65 | yes | 23.3 / 16.4 | 42.3 / 52.3 | 0.019 / 0.034 | 38.3 / 26.5 |
| 10 | `Bayes-M4c-geo` | 34.13 | yes | 20.1 / 14.0 | 44.5 / 49.8 | 0.020 / 0.035 | 32.7 / 23.0 |
| 11 | `Bayes-M17-fill-med` | 39.83 | yes | 23.0 / 16.8 | 48.1 / 60.2 | 0.020 / 0.035 | 37.4 / 26.8 |
| 12 | `Bayes-M13-fill-geo` | 29.44 | yes | 17.3 / 12.2 | 43.2 / 50.9 | 0.020 / 0.036 | 36.0 / 24.4 |
| 13 | `Bayes-M4c-med` | 31.34 | yes | 18.4 / 12.9 | 48.4 / 55.8 | 0.020 / 0.036 | 31.3 / 22.2 |
| 14 | `Bayes-M13-fill-med` | 39.03 | yes | 23.1 / 15.9 | 52.7 / 70.5 | 0.021 / 0.037 | 36.5 / 25.7 |
| 15 | `Bayes-M4-geo` | 24.76 | yes | 14.6 / 10.2 | 52.4 / 57.8 | 0.021 / 0.037 | 30.0 / 20.2 |
| 16 | `Bayes-M8-fill-geo` | 23.30 | yes | 14.0 / 9.3 | 51.4 / 57.2 | 0.021 / 0.037 | 32.7 / 22.0 |
| 17 | `Bayes-M8-fill-med` | 32.75 | yes | 20.4 / 12.4 | 60.4 / 77.2 | 0.022 / 0.039 | 35.2 / 23.7 |
| 18 | `Bayes-M4-med` | 28.66 | yes | 16.7 / 11.9 | 62.3 / 77.6 | 0.022 / 0.038 | 30.4 / 20.9 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 52.54 | yes | 30.9 / 21.7 | 39.0 / 45.3 | 0.017 / 0.029 | 40.6 / 29.3 |
| 2 | `Bayes-M14-fill-geo` | 54.58 | yes | 32.3 / 22.3 | 40.2 / 47.1 | 0.018 / 0.031 | 41.4 / 29.7 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
