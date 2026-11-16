# Model selection — spatiotemporal

_Generated 2026-11-16T10:56:33+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_16112026 (processed_at 2026-11-16T12:00:00) |
| Predictions computed on snapshot | LINELIST_16112026 (processed_at 2026-11-16T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 9, margin 8) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 14 of 19; h2: 14 of 19 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 9, margin 10) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 9, margin 10) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 79%; `Bayes-M10-fill-geo` 9%; `Bayes-M16-fill-geo` 5% — kernel M14-fill 82%, M10-fill 9%, M16-fill 8%, M13c-fill 1%, M17-fill 0%, M4c 0% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 56.15 | yes | 33.2 / 23.0 | 48.3 / 45.7 | 0.018 / 0.030 | 49.0 / 30.8 |
| 2 | `Bayes-M10-fill-geo` | 51.99 | yes | 30.9 / 21.1 | 50.0 / 48.5 | 0.018 / 0.030 | 43.8 / 29.2 |
| 3 | `Bayes-M13c-fill-geo` | 45.05 | yes | 26.8 / 18.2 | 42.5 / 43.0 | 0.019 / 0.032 | 42.9 / 27.9 |
| 4 | `Bayes-M16-fill-geo` | 42.22 | yes | 24.0 / 18.2 | 49.1 / 46.8 | 0.019 / 0.032 | 41.9 / 27.9 |
| 5 | `Bayes-M14-fill-med` | 52.30 | yes | 30.9 / 21.4 | 60.7 / 60.3 | 0.018 / 0.030 | 44.1 / 29.8 |
| 6 | `Bayes-M17-fill-geo` | 42.99 | yes | 25.8 / 17.2 | 46.4 / 46.7 | 0.019 / 0.033 | 44.2 / 27.4 |
| 7 | `Bayes-M10-fill-med` | 47.99 | yes | 27.8 / 20.2 | 65.3 / 65.7 | 0.018 / 0.031 | 42.0 / 28.5 |
| 8 | `Bayes-M16-fill-med` | 46.62 | yes | 26.4 / 20.2 | 58.5 / 58.4 | 0.019 / 0.032 | 41.8 / 28.4 |
| 9 | `Bayes-M13c-fill-med` | 38.49 | yes | 21.9 / 16.6 | 54.0 / 56.1 | 0.019 / 0.033 | 41.4 / 26.9 |
| 10 | `Bayes-M4c-geo` | 33.06 | yes | 19.4 / 13.7 | 50.8 / 50.4 | 0.020 / 0.034 | 34.6 / 23.6 |
| 11 | `Bayes-M17-fill-med` | 40.80 | yes | 23.0 / 17.8 | 64.5 / 67.0 | 0.019 / 0.033 | 40.9 / 27.6 |
| 12 | `Bayes-M13-fill-geo` | 29.14 | yes | 16.8 / 12.4 | 51.5 / 51.9 | 0.020 / 0.034 | 37.8 / 24.6 |
| 13 | `Bayes-M4c-med` | 29.73 | yes | 17.2 / 12.6 | 57.5 / 58.5 | 0.020 / 0.034 | 33.4 / 21.7 |
| 14 | `Bayes-M4-geo` | 22.50 | yes | 12.9 / 9.6 | 57.7 / 57.2 | 0.021 / 0.035 | 29.1 / 19.6 |
| 15 | `Bayes-M13-fill-med` | 37.96 | yes | 22.2 / 15.8 | 75.7 / 81.5 | 0.021 / 0.035 | 37.4 / 25.7 |
| 16 | `Bayes-M8-fill-geo` | 22.06 | yes | 12.7 / 9.3 | 57.5 / 57.2 | 0.021 / 0.035 | 30.8 / 21.0 |
| 17 | `Bayes-M4-med` | 29.31 | **spiky** | 17.3 / 12.0 | 82.4 / 87.6 | 0.021 / 0.036 | 31.9 / 20.8 |
| 18 | `Bayes-M8-fill-med` | 27.64 | **spiky** | 15.4 / 12.3 | 81.7 / 87.2 | 0.021 / 0.036 | 33.3 / 22.9 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 56.42 | yes | 33.0 / 23.4 | 45.8 / 44.0 | 0.017 / 0.028 | 48.8 / 31.7 |
| 2 | `Bayes-M14-fill-geo` | 56.15 | yes | 33.2 / 23.0 | 48.3 / 45.7 | 0.018 / 0.030 | 49.0 / 30.8 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
