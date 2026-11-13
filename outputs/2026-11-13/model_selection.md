# Model selection — spatiotemporal

_Generated 2026-11-13T11:06:30+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_13112026 (processed_at 2026-11-13T12:00:00) |
| Predictions computed on snapshot | LINELIST_13112026 (processed_at 2026-11-13T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 13, margin 4) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 14 of 19; h2: 14 of 19 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 7) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 7) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 64%; `Bayes-M14-fill-geo` 20%; `Bayes-M16-fill-geo` 7% — kernel M14-fill 85%, M16-fill 7%, M10-fill 5%, M13c-fill 2% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 51.10 | yes | 29.5 / 21.6 | 39.5 / 46.2 | 0.017 / 0.030 | 33.6 / 27.1 |
| 2 | `Bayes-M10-fill-geo` | 48.93 | yes | 28.9 / 20.0 | 41.8 / 48.7 | 0.018 / 0.031 | 34.8 / 26.5 |
| 3 | `Bayes-M13c-fill-geo` | 40.24 | yes | 22.8 / 17.4 | 35.3 / 41.6 | 0.019 / 0.032 | 33.0 / 26.1 |
| 4 | `Bayes-M14-fill-med` | 45.43 | yes | 25.7 / 19.7 | 46.5 / 55.6 | 0.018 / 0.031 | 33.7 / 26.8 |
| 5 | `Bayes-M16-fill-geo` | 40.15 | yes | 22.2 / 17.9 | 46.2 / 47.7 | 0.018 / 0.032 | 31.9 / 25.6 |
| 6 | `Bayes-M17-fill-geo` | 36.77 | yes | 20.6 / 16.2 | 38.2 / 45.3 | 0.019 / 0.034 | 32.2 / 25.6 |
| 7 | `Bayes-M10-fill-med` | 44.21 | yes | 25.5 / 18.7 | 50.9 / 60.6 | 0.019 / 0.032 | 34.3 / 26.2 |
| 8 | `Bayes-M16-fill-med` | 41.89 | yes | 23.3 / 18.6 | 52.6 / 56.0 | 0.019 / 0.033 | 32.1 / 25.9 |
| 9 | `Bayes-M13c-fill-med` | 36.53 | yes | 20.9 / 15.7 | 43.5 / 53.0 | 0.019 / 0.034 | 32.5 / 25.0 |
| 10 | `Bayes-M4c-geo` | 30.96 | yes | 17.7 / 13.3 | 42.8 / 49.9 | 0.020 / 0.035 | 28.3 / 21.7 |
| 11 | `Bayes-M17-fill-med` | 36.60 | yes | 20.6 / 16.0 | 50.6 / 62.2 | 0.020 / 0.035 | 31.7 / 25.1 |
| 12 | `Bayes-M13-fill-geo` | 27.59 | yes | 15.4 / 12.2 | 43.7 / 51.1 | 0.020 / 0.036 | 30.1 / 23.3 |
| 13 | `Bayes-M4c-med` | 29.36 | yes | 16.9 / 12.5 | 47.8 / 56.1 | 0.020 / 0.036 | 28.5 / 21.2 |
| 14 | `Bayes-M13-fill-med` | 32.66 | yes | 18.0 / 14.6 | 61.1 / 77.3 | 0.021 / 0.037 | 29.8 / 23.9 |
| 15 | `Bayes-M4-geo` | 24.43 | yes | 14.2 / 10.3 | 51.4 / 57.5 | 0.021 / 0.037 | 24.3 / 19.9 |
| 16 | `Bayes-M8-fill-geo` | 22.18 | yes | 12.9 / 9.3 | 51.9 / 57.4 | 0.021 / 0.037 | 27.6 / 20.3 |
| 17 | `Bayes-M4-med` | 27.71 | **spiky** | 15.9 / 11.8 | 69.6 / 84.2 | 0.022 / 0.038 | 26.2 / 19.7 |
| 18 | `Bayes-M8-fill-med` | 27.02 | **spiky** | 15.8 / 11.3 | 68.4 / 83.8 | 0.022 / 0.039 | 28.7 / 21.6 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 48.57 | yes | 27.5 / 21.1 | 37.7 / 42.4 | 0.017 / 0.029 | 34.7 / 27.6 |
| 2 | `Bayes-M14-fill-geo` | 51.10 | yes | 29.5 / 21.6 | 39.5 / 46.2 | 0.017 / 0.030 | 33.6 / 27.1 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
