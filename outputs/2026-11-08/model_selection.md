# Model selection — spatiotemporal

_Generated 2026-11-08T10:39:18+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_08112026 (processed_at 2026-11-08T12:00:00) |
| Predictions computed on snapshot | LINELIST_08112026 (processed_at 2026-11-08T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 13) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 13 of 18; h2: 13 of 18 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 13) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 11) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 79%; `Bayes-M16-fill-geo` 7%; `Bayes-M10-fill-geo` 6% — kernel M14-fill 84%, M16-fill 9%, M10-fill 6%, M13c-fill 1% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 54.62 | yes | 30.9 / 23.7 | 45.2 / 52.0 | 0.018 / 0.030 | 45.0 / 30.7 |
| 2 | `Bayes-M10-fill-geo` | 51.85 | yes | 29.8 / 22.0 | 46.3 / 54.1 | 0.018 / 0.031 | 41.3 / 28.7 |
| 3 | `Bayes-M13c-fill-geo` | 40.41 | yes | 22.8 / 17.7 | 38.3 / 44.1 | 0.019 / 0.033 | 41.5 / 27.7 |
| 4 | `Bayes-M16-fill-geo` | 41.17 | yes | 22.4 / 18.7 | 45.8 / 53.2 | 0.019 / 0.032 | 39.0 / 27.9 |
| 5 | `Bayes-M14-fill-med` | 52.44 | yes | 30.2 / 22.2 | 52.6 / 63.2 | 0.018 / 0.031 | 44.1 / 29.8 |
| 6 | `Bayes-M17-fill-geo` | 38.46 | yes | 21.6 / 16.8 | 42.0 / 48.6 | 0.019 / 0.033 | 40.9 / 27.1 |
| 7 | `Bayes-M10-fill-med` | 49.99 | yes | 28.4 / 21.6 | 55.4 / 67.3 | 0.018 / 0.031 | 40.5 / 28.6 |
| 8 | `Bayes-M16-fill-med` | 47.26 | yes | 26.7 / 20.5 | 51.8 / 62.3 | 0.019 / 0.033 | 40.8 / 28.1 |
| 9 | `Bayes-M13c-fill-med` | 37.88 | yes | 21.4 / 16.5 | 46.6 / 56.3 | 0.019 / 0.034 | 38.8 / 26.4 |
| 10 | `Bayes-M4c-geo` | 31.66 | yes | 18.0 / 13.7 | 48.3 / 52.3 | 0.020 / 0.034 | 32.1 / 22.8 |
| 11 | `Bayes-M17-fill-med` | 40.09 | yes | 22.8 / 17.2 | 53.6 / 66.2 | 0.020 / 0.034 | 40.0 / 27.1 |
| 12 | `Bayes-M13-fill-geo` | 28.63 | yes | 16.3 / 12.3 | 46.5 / 54.1 | 0.020 / 0.035 | 37.7 / 24.6 |
| 13 | `Bayes-M4c-med` | 30.12 | yes | 17.3 / 12.8 | 51.5 / 58.6 | 0.020 / 0.035 | 31.4 / 21.8 |
| 14 | `Bayes-M13-fill-med` | 38.31 | yes | 22.2 / 16.1 | 58.1 / 76.7 | 0.021 / 0.036 | 38.5 / 25.7 |
| 15 | `Bayes-M4-geo` | 23.94 | yes | 13.8 / 10.1 | 54.8 / 59.9 | 0.021 / 0.036 | 29.1 / 19.7 |
| 16 | `Bayes-M8-fill-geo` | 22.44 | yes | 12.9 / 9.5 | 54.4 / 59.9 | 0.021 / 0.037 | 32.3 / 21.4 |
| 17 | `Bayes-M4-med` | 28.41 | yes | 16.5 / 11.9 | 65.9 / 81.9 | 0.021 / 0.037 | 30.7 / 20.3 |
| 18 | `Bayes-M8-fill-med` | 26.85 | yes | 14.9 / 12.0 | 65.4 / 81.2 | 0.022 / 0.038 | 32.8 / 23.0 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 55.11 | yes | 31.7 / 23.4 | 43.5 / 49.1 | 0.017 / 0.028 | 45.5 / 30.4 |
| 2 | `Bayes-M14-fill-geo` | 54.62 | yes | 30.9 / 23.7 | 45.2 / 52.0 | 0.018 / 0.030 | 45.0 / 30.7 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
