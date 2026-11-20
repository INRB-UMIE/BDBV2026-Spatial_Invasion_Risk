# Model selection — spatiotemporal

_Generated 2026-11-20T11:35:09+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_20112026 (processed_at 2026-11-20T12:00:00) |
| Predictions computed on snapshot | LINELIST_20112026 (processed_at 2026-11-20T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 14, margin 2) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 15 of 20; h2: 15 of 20 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 7) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 7) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 63%; `Bayes-M14-fill-geo` 22%; `Bayes-M16-fill-geo` 8% — kernel M14-fill 85%, M16-fill 8%, M10-fill 4%, M13c-fill 2% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 52.60 | yes | 30.5 / 22.1 | 39.3 / 43.7 | 0.017 / 0.029 | 41.7 / 27.6 |
| 2 | `Bayes-M10-fill-geo` | 49.01 | yes | 28.1 / 20.9 | 41.6 / 46.2 | 0.017 / 0.030 | 41.3 / 27.1 |
| 3 | `Bayes-M13c-fill-geo` | 40.54 | yes | 22.7 / 17.8 | 35.1 / 39.6 | 0.018 / 0.031 | 42.6 / 26.7 |
| 4 | `Bayes-M14-fill-med` | 46.15 | yes | 25.7 / 20.4 | 46.5 / 52.7 | 0.017 / 0.030 | 40.9 / 27.5 |
| 5 | `Bayes-M16-fill-geo` | 41.86 | yes | 22.4 / 19.5 | 46.5 / 45.2 | 0.018 / 0.030 | 38.5 / 26.7 |
| 6 | `Bayes-M17-fill-geo` | 37.21 | yes | 20.7 / 16.5 | 37.9 / 42.8 | 0.018 / 0.033 | 39.3 / 25.2 |
| 7 | `Bayes-M10-fill-med` | 44.43 | yes | 25.4 / 19.0 | 50.9 / 57.7 | 0.018 / 0.031 | 39.0 / 25.2 |
| 8 | `Bayes-M16-fill-med` | 43.50 | yes | 23.4 / 20.1 | 52.6 / 52.9 | 0.018 / 0.032 | 38.5 / 27.1 |
| 9 | `Bayes-M13c-fill-med` | 37.24 | yes | 20.8 / 16.4 | 43.6 / 50.3 | 0.019 / 0.033 | 39.7 / 24.7 |
| 10 | `Bayes-M4c-geo` | 31.39 | yes | 17.7 / 13.7 | 42.8 / 47.6 | 0.019 / 0.033 | 32.0 / 20.4 |
| 11 | `Bayes-M17-fill-med` | 37.29 | yes | 20.5 / 16.8 | 50.8 / 59.2 | 0.019 / 0.033 | 38.2 / 25.0 |
| 12 | `Bayes-M13-fill-geo` | 27.78 | yes | 15.4 / 12.4 | 43.6 / 48.6 | 0.019 / 0.034 | 35.4 / 21.9 |
| 13 | `Bayes-M4c-med` | 29.49 | yes | 16.9 / 12.6 | 47.7 / 53.6 | 0.020 / 0.034 | 30.8 / 19.2 |
| 14 | `Bayes-M13-fill-med` | 34.39 | yes | 18.0 / 16.4 | 61.0 / 73.3 | 0.020 / 0.036 | 35.2 / 24.4 |
| 15 | `Bayes-M4-geo` | 24.44 | yes | 14.0 / 10.4 | 51.5 / 55.0 | 0.020 / 0.036 | 27.2 / 17.5 |
| 16 | `Bayes-M8-fill-geo` | 22.46 | yes | 13.1 / 9.4 | 51.4 / 55.0 | 0.020 / 0.036 | 31.4 / 17.5 |
| 17 | `Bayes-M4-med` | 28.07 | **spiky** | 15.9 / 12.1 | 69.5 / 79.9 | 0.021 / 0.037 | 27.9 / 18.2 |
| 18 | `Bayes-M8-fill-med` | 27.23 | **spiky** | 15.8 / 11.5 | 68.4 / 79.8 | 0.021 / 0.037 | 32.0 / 19.0 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 49.10 | yes | 27.6 / 21.5 | 37.5 / 40.2 | 0.016 / 0.028 | 41.9 / 28.1 |
| 2 | `Bayes-M14-fill-geo` | 52.60 | yes | 30.5 / 22.1 | 39.3 / 43.7 | 0.017 / 0.029 | 41.7 / 27.6 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
