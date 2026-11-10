# Model selection — spatiotemporal

_Generated 2026-11-10T10:42:43+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_10112026 (processed_at 2026-11-10T12:00:00) |
| Predictions computed on snapshot | LINELIST_10112026 (processed_at 2026-11-10T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 12, margin 4) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 13 of 18; h2: 13 of 18 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 12, margin 8) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 12, margin 8) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 75%; `Bayes-M10-fill-geo` 10%; `Bayes-M14-fill-geo` 6% — kernel M14-fill 81%, M10-fill 10%, M16-fill 7%, M13c-fill 2%, M13-fill 0% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 57.90 | yes | 34.6 / 23.3 | 51.0 / 48.8 | 0.017 / 0.030 | 50.6 / 30.7 |
| 2 | `Bayes-M10-fill-geo` | 53.92 | yes | 31.7 / 22.2 | 52.8 / 51.5 | 0.018 / 0.030 | 46.5 / 29.5 |
| 3 | `Bayes-M13c-fill-geo` | 45.62 | yes | 26.8 / 18.8 | 44.3 / 44.2 | 0.018 / 0.032 | 49.8 / 29.1 |
| 4 | `Bayes-M16-fill-geo` | 41.02 | yes | 22.7 / 18.3 | 49.9 / 48.2 | 0.018 / 0.032 | 42.7 / 27.7 |
| 5 | `Bayes-M17-fill-geo` | 40.15 | yes | 23.4 / 16.8 | 48.3 / 48.7 | 0.019 / 0.033 | 47.3 / 27.4 |
| 6 | `Bayes-M14-fill-med` | 53.64 | yes | 31.6 / 22.1 | 63.0 / 63.0 | 0.018 / 0.031 | 49.1 / 30.0 |
| 7 | `Bayes-M10-fill-med` | 51.09 | yes | 29.9 / 21.2 | 67.6 / 68.2 | 0.018 / 0.031 | 44.3 / 28.3 |
| 8 | `Bayes-M16-fill-med` | 44.68 | yes | 26.0 / 18.6 | 59.5 / 59.8 | 0.019 / 0.033 | 45.1 / 28.3 |
| 9 | `Bayes-M13c-fill-med` | 39.36 | yes | 22.7 / 16.6 | 54.1 / 56.9 | 0.019 / 0.033 | 46.1 / 28.0 |
| 10 | `Bayes-M4c-geo` | 34.87 | yes | 20.2 / 14.7 | 52.6 / 51.6 | 0.019 / 0.033 | 34.0 / 22.8 |
| 11 | `Bayes-M13-fill-geo` | 29.93 | yes | 17.2 / 12.8 | 52.5 / 53.3 | 0.020 / 0.034 | 43.7 / 24.4 |
| 12 | `Bayes-M17-fill-med` | 39.82 | yes | 22.6 / 17.2 | 64.1 / 67.4 | 0.019 / 0.033 | 43.5 / 27.9 |
| 13 | `Bayes-M4c-med` | 31.05 | yes | 17.8 / 13.2 | 57.7 / 59.2 | 0.020 / 0.034 | 33.6 / 21.7 |
| 14 | `Bayes-M13-fill-med` | 37.98 | yes | 22.8 / 15.2 | 69.6 / 77.4 | 0.020 / 0.035 | 43.2 / 26.0 |
| 15 | `Bayes-M4-geo` | 23.83 | yes | 13.4 / 10.4 | 58.5 / 58.4 | 0.020 / 0.035 | 26.2 / 19.4 |
| 16 | `Bayes-M8-fill-geo` | 23.31 | yes | 13.5 / 9.8 | 58.8 / 58.7 | 0.020 / 0.035 | 33.9 / 20.7 |
| 17 | `Bayes-M4-med` | 30.44 | yes | 18.1 / 12.4 | 75.9 / 83.2 | 0.021 / 0.036 | 32.8 / 20.9 |
| 18 | `Bayes-M8-fill-med` | 28.35 | yes | 15.9 / 12.5 | 76.4 / 83.2 | 0.021 / 0.037 | 37.3 / 22.9 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 57.34 | yes | 34.0 / 23.4 | 49.3 / 47.7 | 0.017 / 0.028 | 51.1 / 31.6 |
| 2 | `Bayes-M14-fill-geo` | 57.90 | yes | 34.6 / 23.3 | 51.0 / 48.8 | 0.017 / 0.030 | 50.6 / 30.7 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
