# Model selection — spatiotemporal

_Generated 2026-11-23T10:23:15+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_23112026 (processed_at 2026-11-23T12:00:00) |
| Predictions computed on snapshot | LINELIST_23112026 (processed_at 2026-11-23T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 9, margin 9) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 15 of 20; h2: 15 of 20 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 9, margin 14) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 9, margin 12) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 84%; `Bayes-M16-fill-geo` 8%; `Bayes-M10-fill-geo` 3% — kernel M14-fill 86%, M16-fill 11%, M10-fill 3%, M13c-fill 1% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 61.08 | yes | 37.1 / 24.0 | 45.4 / 43.3 | 0.017 / 0.029 | 51.5 / 31.4 |
| 2 | `Bayes-M10-fill-geo` | 53.32 | yes | 31.8 / 21.6 | 47.0 / 45.8 | 0.017 / 0.029 | 44.9 / 29.5 |
| 3 | `Bayes-M16-fill-geo` | 49.09 | yes | 28.7 / 20.4 | 46.1 / 43.9 | 0.018 / 0.030 | 45.0 / 29.2 |
| 4 | `Bayes-M13c-fill-geo` | 47.31 | yes | 27.9 / 19.4 | 39.8 / 39.6 | 0.018 / 0.031 | 43.5 / 29.1 |
| 5 | `Bayes-M14-fill-med` | 55.75 | yes | 33.0 / 22.7 | 57.0 / 56.4 | 0.018 / 0.030 | 47.3 / 30.6 |
| 6 | `Bayes-M16-fill-med` | 54.50 | yes | 31.5 / 23.0 | 54.6 / 54.5 | 0.018 / 0.031 | 45.2 / 30.3 |
| 7 | `Bayes-M17-fill-geo` | 44.07 | yes | 26.3 / 17.8 | 43.9 / 43.8 | 0.019 / 0.032 | 44.9 / 27.8 |
| 8 | `Bayes-M10-fill-med` | 50.22 | yes | 29.2 / 21.1 | 61.4 / 61.5 | 0.018 / 0.030 | 42.9 / 29.0 |
| 9 | `Bayes-M13c-fill-med` | 42.05 | yes | 24.2 / 17.9 | 50.6 / 51.5 | 0.019 / 0.032 | 42.8 / 27.5 |
| 10 | `Bayes-M4c-geo` | 34.25 | yes | 19.8 / 14.4 | 48.2 / 47.2 | 0.019 / 0.033 | 33.9 / 24.0 |
| 11 | `Bayes-M17-fill-med` | 44.41 | yes | 25.3 / 19.1 | 60.5 / 61.3 | 0.019 / 0.032 | 42.4 / 28.4 |
| 12 | `Bayes-M13-fill-geo` | 29.77 | yes | 17.2 / 12.6 | 48.8 / 48.6 | 0.020 / 0.033 | 37.9 / 24.6 |
| 13 | `Bayes-M4c-med` | 30.49 | yes | 17.5 / 13.0 | 54.5 / 54.2 | 0.020 / 0.033 | 33.5 / 21.9 |
| 14 | `Bayes-M13-fill-med` | 45.48 | yes | 26.9 / 18.5 | 70.9 / 72.4 | 0.020 / 0.034 | 40.6 / 27.6 |
| 15 | `Bayes-M4-geo` | 22.93 | yes | 13.0 / 9.9 | 55.0 / 53.9 | 0.020 / 0.034 | 29.1 / 19.7 |
| 16 | `Bayes-M8-fill-geo` | 22.79 | yes | 13.2 / 9.6 | 55.0 / 54.3 | 0.020 / 0.034 | 31.3 / 21.2 |
| 17 | `Bayes-M4-med` | 30.95 | yes | 18.4 / 12.6 | 77.4 / 78.2 | 0.021 / 0.035 | 32.5 / 21.2 |
| 18 | `Bayes-M8-fill-med` | 28.69 | yes | 16.0 / 12.6 | 76.9 / 78.0 | 0.021 / 0.036 | 33.6 / 23.1 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 59.19 | yes | 34.4 / 24.8 | 43.3 / 41.7 | 0.016 / 0.027 | 49.7 / 32.3 |
| 2 | `Bayes-M14-fill-geo` | 61.08 | yes | 37.1 / 24.0 | 45.4 / 43.3 | 0.017 / 0.029 | 51.5 / 31.4 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
