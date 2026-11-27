# Model selection — spatiotemporal

_Generated 2026-11-27T14:40:40+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_26112026 (processed_at 2026-11-26T12:00:00) |
| Predictions computed on snapshot | LINELIST_26112026 (processed_at 2026-11-26T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo` (kernel M14-fill, form fixed; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 11) — **differs from the current pick** |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 16 of 21; h2: 16 of 21 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 5) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 5) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 59%; `Bayes-M14-fill-geo` 22%; `Bayes-M16-fill-geo` 11% — kernel M14-fill 82%, M16-fill 12%, M10-fill 4%, M13c-fill 3% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 51.51 | yes | 27.9 / 23.6 | 41.3 / 41.8 | 0.016 / 0.028 | 32.8 / 28.2 |
| 2 | `Bayes-M10-fill-geo` | 50.15 | yes | 28.2 / 21.9 | 43.5 / 44.2 | 0.017 / 0.029 | 33.0 / 27.2 |
| 3 | `Bayes-M13c-fill-geo` | 41.95 | yes | 22.7 / 19.3 | 36.4 / 37.8 | 0.017 / 0.031 | 31.0 / 26.1 |
| 4 | `Bayes-M16-fill-geo` | 44.62 | yes | 24.1 / 20.5 | 46.8 / 42.2 | 0.017 / 0.030 | 30.9 / 26.1 |
| 5 | `Bayes-M14-fill-med` | 46.06 | yes | 25.1 / 21.0 | 48.1 / 50.5 | 0.017 / 0.030 | 32.2 / 27.0 |
| 6 | `Bayes-M17-fill-geo` | 38.99 | yes | 20.9 / 18.1 | 39.1 / 41.1 | 0.018 / 0.032 | 31.1 / 25.6 |
| 7 | `Bayes-M16-fill-med` | 47.41 | yes | 26.8 / 20.6 | 52.6 / 49.9 | 0.018 / 0.031 | 32.6 / 27.0 |
| 8 | `Bayes-M10-fill-med` | 44.28 | yes | 24.3 / 20.0 | 52.8 / 55.1 | 0.017 / 0.031 | 32.2 / 26.3 |
| 9 | `Bayes-M13c-fill-med` | 38.41 | yes | 20.9 / 17.5 | 44.6 / 48.5 | 0.018 / 0.032 | 30.4 / 25.4 |
| 10 | `Bayes-M4c-geo` | 31.96 | yes | 17.2 / 14.7 | 45.2 / 45.7 | 0.019 / 0.033 | 26.3 / 21.8 |
| 11 | `Bayes-M13-fill-geo` | 28.02 | yes | 14.6 / 13.4 | 43.8 / 46.4 | 0.019 / 0.034 | 28.1 / 23.2 |
| 12 | `Bayes-M17-fill-med` | 38.51 | yes | 20.8 / 17.7 | 51.6 / 56.8 | 0.018 / 0.033 | 30.5 / 25.7 |
| 13 | `Bayes-M4c-med` | 28.30 | yes | 15.6 / 12.7 | 50.3 / 51.8 | 0.019 / 0.034 | 26.3 / 20.9 |
| 14 | `Bayes-M13-fill-med` | 38.09 | yes | 20.6 / 17.5 | 61.3 / 69.9 | 0.019 / 0.035 | 29.9 / 25.5 |
| 15 | `Bayes-M4-geo` | 22.55 | yes | 12.0 / 10.5 | 52.0 / 52.8 | 0.020 / 0.035 | 21.2 / 18.7 |
| 16 | `Bayes-M8-fill-geo` | 22.01 | yes | 12.1 / 10.0 | 52.8 / 52.9 | 0.020 / 0.036 | 25.5 / 20.2 |
| 17 | `Bayes-M4-med` | 26.90 | yes | 14.9 / 12.0 | 69.8 / 76.5 | 0.020 / 0.037 | 24.6 / 19.5 |
| 18 | `Bayes-M8-fill-med` | 26.92 | yes | 14.6 / 12.3 | 68.8 / 76.5 | 0.020 / 0.037 | 26.8 / 21.9 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 51.63 | yes | 28.3 / 23.3 | 39.1 / 35.8 | 0.016 / 0.028 | 33.1 / 28.4 |
| 2 | `Bayes-M14-fill-geo` | 51.51 | yes | 27.9 / 23.6 | 41.3 / 41.8 | 0.016 / 0.028 | 32.8 / 28.2 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
