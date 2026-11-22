# Model selection — spatiotemporal

_Generated 2026-11-22T09:47:30+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_22112026 (processed_at 2026-11-22T12:00:00) |
| Predictions computed on snapshot | LINELIST_22112026 (processed_at 2026-11-22T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 11) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 15 of 20; h2: 15 of 20 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 15) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 13) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 74%; `Bayes-M16-fill-geo` 12%; `Bayes-M14-fill-geo` 10% — kernel M14-fill 84%, M16-fill 14%, M10-fill 2% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 58.30 | yes | 33.4 / 24.9 | 40.0 / 46.5 | 0.017 / 0.028 | 46.3 / 31.4 |
| 2 | `Bayes-M10-fill-geo` | 53.34 | yes | 30.5 / 22.8 | 41.4 / 49.0 | 0.017 / 0.029 | 41.4 / 29.2 |
| 3 | `Bayes-M16-fill-geo` | 48.95 | yes | 27.9 / 21.1 | 40.4 / 47.6 | 0.018 / 0.030 | 41.8 / 28.8 |
| 4 | `Bayes-M14-fill-med` | 55.78 | yes | 32.3 / 23.5 | 47.5 / 57.3 | 0.017 / 0.029 | 45.1 / 30.5 |
| 5 | `Bayes-M13c-fill-geo` | 42.89 | yes | 23.9 / 18.9 | 35.6 / 40.6 | 0.018 / 0.031 | 41.7 / 28.7 |
| 6 | `Bayes-M16-fill-med` | 54.88 | yes | 31.5 / 23.4 | 46.4 / 56.4 | 0.018 / 0.031 | 43.4 / 29.7 |
| 7 | `Bayes-M17-fill-geo` | 40.38 | yes | 22.6 / 17.8 | 39.1 / 44.9 | 0.019 / 0.032 | 41.2 / 27.8 |
| 8 | `Bayes-M10-fill-med` | 51.42 | yes | 28.9 / 22.5 | 50.5 / 61.4 | 0.018 / 0.030 | 40.7 / 29.0 |
| 9 | `Bayes-M13c-fill-med` | 41.35 | yes | 23.5 / 17.8 | 44.6 / 52.2 | 0.019 / 0.032 | 39.8 / 27.0 |
| 10 | `Bayes-M4c-geo` | 33.51 | yes | 18.9 / 14.6 | 44.5 / 48.2 | 0.019 / 0.033 | 32.6 / 23.3 |
| 11 | `Bayes-M17-fill-med` | 43.40 | yes | 24.9 / 18.5 | 51.5 / 61.2 | 0.019 / 0.032 | 41.0 / 27.8 |
| 12 | `Bayes-M13-fill-geo` | 28.94 | yes | 16.5 / 12.4 | 44.1 / 50.5 | 0.020 / 0.034 | 37.5 / 24.6 |
| 13 | `Bayes-M4c-med` | 30.85 | yes | 17.5 / 13.3 | 49.4 / 54.5 | 0.020 / 0.033 | 31.4 / 22.0 |
| 14 | `Bayes-M13-fill-med` | 45.77 | yes | 26.8 / 19.0 | 61.8 / 73.2 | 0.020 / 0.035 | 41.0 / 27.3 |
| 15 | `Bayes-M4-geo` | 24.01 | yes | 13.8 / 10.2 | 51.9 / 55.9 | 0.020 / 0.034 | 28.7 / 19.6 |
| 16 | `Bayes-M8-fill-geo` | 22.83 | yes | 13.2 / 9.6 | 51.5 / 56.1 | 0.020 / 0.035 | 32.1 / 21.7 |
| 17 | `Bayes-M4-med` | 29.91 | yes | 17.4 / 12.5 | 68.9 / 78.1 | 0.021 / 0.035 | 31.0 / 20.7 |
| 18 | `Bayes-M8-fill-med` | 27.96 | yes | 15.6 / 12.4 | 68.2 / 77.6 | 0.021 / 0.036 | 33.0 / 23.1 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 57.81 | yes | 33.2 / 24.6 | 38.5 / 44.9 | 0.016 / 0.027 | 46.1 / 31.2 |
| 2 | `Bayes-M14-fill-geo` | 58.30 | yes | 33.4 / 24.9 | 40.0 / 46.5 | 0.017 / 0.028 | 46.3 / 31.4 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
