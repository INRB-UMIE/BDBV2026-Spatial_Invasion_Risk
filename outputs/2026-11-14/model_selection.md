# Model selection — spatiotemporal

_Generated 2026-11-14T09:34:11+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_14112026 (processed_at 2026-11-14T12:00:00) |
| Predictions computed on snapshot | LINELIST_14112026 (processed_at 2026-11-14T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 19, margin 5) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 14 of 19; h2: 14 of 19 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 10) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 8) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 72%; `Bayes-M14-fill-geo` 14%; `Bayes-M16-fill-geo` 7% — kernel M14-fill 86%, M16-fill 8%, M10-fill 5%, M13c-fill 0% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 54.90 | yes | 32.4 / 22.5 | 38.5 / 45.0 | 0.018 / 0.030 | 41.4 / 29.6 |
| 2 | `Bayes-M10-fill-geo` | 50.57 | yes | 29.5 / 21.0 | 39.9 / 46.9 | 0.018 / 0.030 | 38.6 / 28.5 |
| 3 | `Bayes-M13c-fill-geo` | 41.89 | yes | 24.3 / 17.6 | 33.9 / 39.8 | 0.019 / 0.032 | 39.3 / 27.2 |
| 4 | `Bayes-M14-fill-med` | 50.89 | yes | 28.8 / 22.1 | 45.5 / 54.6 | 0.018 / 0.031 | 40.5 / 29.5 |
| 5 | `Bayes-M16-fill-geo` | 42.79 | yes | 24.8 / 18.0 | 39.9 / 46.5 | 0.018 / 0.032 | 38.1 / 27.2 |
| 6 | `Bayes-M17-fill-geo` | 39.07 | yes | 22.6 / 16.5 | 37.2 / 43.9 | 0.019 / 0.033 | 38.7 / 27.2 |
| 7 | `Bayes-M10-fill-med` | 49.30 | yes | 28.5 / 20.8 | 48.2 / 58.9 | 0.018 / 0.031 | 39.6 / 28.1 |
| 8 | `Bayes-M16-fill-med` | 46.37 | yes | 26.5 / 19.9 | 45.6 / 54.5 | 0.019 / 0.032 | 37.8 / 27.8 |
| 9 | `Bayes-M13c-fill-med` | 39.41 | yes | 23.0 / 16.4 | 43.0 / 52.1 | 0.019 / 0.033 | 38.1 / 26.5 |
| 10 | `Bayes-M4c-geo` | 33.83 | yes | 19.8 / 14.1 | 43.3 / 48.2 | 0.020 / 0.034 | 32.6 / 23.2 |
| 11 | `Bayes-M17-fill-med` | 39.52 | yes | 22.7 / 16.8 | 48.8 / 59.8 | 0.019 / 0.034 | 37.2 / 26.7 |
| 12 | `Bayes-M13-fill-geo` | 29.29 | yes | 17.1 / 12.1 | 43.5 / 50.7 | 0.020 / 0.035 | 36.0 / 24.2 |
| 13 | `Bayes-M4c-med` | 31.20 | yes | 18.3 / 12.9 | 48.7 / 55.5 | 0.020 / 0.035 | 31.2 / 22.2 |
| 14 | `Bayes-M13-fill-med` | 38.76 | yes | 22.8 / 15.9 | 59.0 / 74.9 | 0.021 / 0.036 | 36.2 / 25.7 |
| 15 | `Bayes-M4-geo` | 25.02 | yes | 14.9 / 10.2 | 52.4 / 57.3 | 0.021 / 0.036 | 29.8 / 20.4 |
| 16 | `Bayes-M8-fill-geo` | 23.41 | yes | 14.0 / 9.5 | 51.3 / 56.8 | 0.021 / 0.036 | 32.9 / 21.9 |
| 17 | `Bayes-M8-fill-med` | 29.70 | yes | 17.3 / 12.4 | 66.3 / 81.3 | 0.022 / 0.038 | 33.7 / 23.7 |
| 18 | `Bayes-M4-med` | 28.45 | **spiky** | 16.5 / 11.9 | 68.3 / 81.9 | 0.021 / 0.037 | 30.2 / 20.8 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 53.65 | yes | 31.7 / 21.9 | 36.5 / 43.0 | 0.017 / 0.028 | 40.8 / 29.3 |
| 2 | `Bayes-M14-fill-geo` | 54.90 | yes | 32.4 / 22.5 | 38.5 / 45.0 | 0.018 / 0.030 | 41.4 / 29.6 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
