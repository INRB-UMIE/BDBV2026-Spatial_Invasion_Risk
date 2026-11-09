# Model selection — spatiotemporal

_Generated 2026-11-09T10:59:27+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_09112026 (processed_at 2026-11-09T12:00:00) |
| Predictions computed on snapshot | LINELIST_09112026 (processed_at 2026-11-09T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 6) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 13 of 18; h2: 13 of 18 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 13) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 11) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 80%; `Bayes-M16-fill-geo` 6%; `Bayes-M10-fill-geo` 6% — kernel M14-fill 84%, M16-fill 8%, M10-fill 6%, M13c-fill 2%, M4c 0% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 58.58 | yes | 35.8 / 22.8 | 51.1 / 48.0 | 0.018 / 0.030 | 50.9 / 31.1 |
| 2 | `Bayes-M10-fill-geo` | 52.41 | yes | 31.4 / 21.0 | 52.8 / 50.9 | 0.018 / 0.031 | 44.2 / 29.3 |
| 3 | `Bayes-M13c-fill-geo` | 45.18 | yes | 26.9 / 18.3 | 43.0 / 43.5 | 0.019 / 0.033 | 43.0 / 28.4 |
| 4 | `Bayes-M16-fill-geo` | 43.06 | yes | 24.1 / 19.0 | 51.3 / 48.3 | 0.019 / 0.033 | 42.1 / 28.8 |
| 5 | `Bayes-M17-fill-geo` | 42.76 | yes | 25.6 / 17.1 | 47.6 / 47.5 | 0.019 / 0.034 | 44.2 / 27.4 |
| 6 | `Bayes-M14-fill-med` | 52.54 | yes | 31.1 / 21.5 | 62.6 / 62.1 | 0.018 / 0.031 | 46.2 / 29.8 |
| 7 | `Bayes-M10-fill-med` | 48.24 | yes | 28.0 / 20.2 | 67.3 / 67.2 | 0.018 / 0.032 | 42.3 / 28.5 |
| 8 | `Bayes-M16-fill-med` | 46.46 | yes | 26.3 / 20.2 | 60.8 / 60.1 | 0.019 / 0.033 | 41.8 / 28.4 |
| 9 | `Bayes-M13c-fill-med` | 38.68 | yes | 22.1 / 16.6 | 52.9 / 55.5 | 0.019 / 0.034 | 41.7 / 26.9 |
| 10 | `Bayes-M4c-geo` | 33.35 | yes | 19.6 / 13.7 | 51.9 / 51.6 | 0.020 / 0.035 | 34.6 / 23.3 |
| 11 | `Bayes-M17-fill-med` | 41.09 | yes | 23.3 / 17.8 | 63.1 / 66.0 | 0.020 / 0.034 | 41.1 / 27.6 |
| 12 | `Bayes-M13-fill-geo` | 29.40 | yes | 17.0 / 12.4 | 52.1 / 52.4 | 0.020 / 0.035 | 38.2 / 24.7 |
| 13 | `Bayes-M4c-med` | 29.88 | yes | 17.3 / 12.6 | 56.7 / 57.8 | 0.020 / 0.035 | 33.6 / 21.8 |
| 14 | `Bayes-M13-fill-med` | 38.25 | yes | 22.4 / 15.8 | 68.1 / 75.6 | 0.021 / 0.036 | 37.7 / 25.7 |
| 15 | `Bayes-M4-geo` | 22.54 | yes | 12.9 / 9.6 | 58.6 / 57.8 | 0.021 / 0.036 | 29.3 / 19.6 |
| 16 | `Bayes-M8-fill-geo` | 22.30 | yes | 13.0 / 9.3 | 58.6 / 58.0 | 0.021 / 0.036 | 31.3 / 21.0 |
| 17 | `Bayes-M4-med` | 29.50 | yes | 17.5 / 12.0 | 75.0 / 81.5 | 0.021 / 0.037 | 32.1 / 20.9 |
| 18 | `Bayes-M8-fill-med` | 27.83 | yes | 15.5 / 12.3 | 74.4 / 81.1 | 0.022 / 0.037 | 33.5 / 22.9 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 56.31 | yes | 33.0 / 23.3 | 48.3 / 46.0 | 0.017 / 0.029 | 49.1 / 31.6 |
| 2 | `Bayes-M14-fill-geo` | 58.58 | yes | 35.8 / 22.8 | 51.1 / 48.0 | 0.018 / 0.030 | 50.9 / 31.1 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
