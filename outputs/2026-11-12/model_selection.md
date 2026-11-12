# Model selection — spatiotemporal

_Generated 2026-11-12T11:27:25+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_12112026 (processed_at 2026-11-12T12:00:00) |
| Predictions computed on snapshot | LINELIST_12112026 (processed_at 2026-11-12T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo` (kernel M14-fill, form fixed; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 7) — **differs from the current pick** |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 14 of 19; h2: 14 of 19 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 5) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 5) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 63%; `Bayes-M14-fill-geo` 18%; `Bayes-M16-fill-geo` 8% — kernel M14-fill 81%, M16-fill 8%, M10-fill 7%, M13c-fill 4% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 48.36 | yes | 25.8 / 22.6 | 44.0 / 45.8 | 0.017 / 0.030 | 32.3 / 28.0 |
| 2 | `Bayes-M10-fill-geo` | 48.23 | yes | 27.5 / 20.8 | 46.0 / 48.3 | 0.018 / 0.031 | 32.2 / 26.2 |
| 3 | `Bayes-M13c-fill-geo` | 38.72 | yes | 20.8 / 17.9 | 39.1 / 41.9 | 0.018 / 0.033 | 30.9 / 25.5 |
| 4 | `Bayes-M14-fill-med` | 42.86 | yes | 23.1 / 19.7 | 51.6 / 55.4 | 0.018 / 0.032 | 31.4 / 26.4 |
| 5 | `Bayes-M16-fill-geo` | 37.72 | yes | 20.1 / 17.6 | 50.1 / 46.5 | 0.018 / 0.032 | 29.4 / 24.7 |
| 6 | `Bayes-M17-fill-geo` | 36.59 | yes | 19.8 / 16.8 | 41.5 / 45.0 | 0.019 / 0.034 | 30.3 / 24.8 |
| 7 | `Bayes-M10-fill-med` | 42.43 | yes | 23.4 / 19.1 | 56.4 / 60.4 | 0.018 / 0.033 | 31.8 / 26.0 |
| 8 | `Bayes-M16-fill-med` | 40.64 | yes | 22.8 / 17.8 | 56.5 / 54.7 | 0.019 / 0.034 | 30.6 / 25.3 |
| 9 | `Bayes-M13c-fill-med` | 35.36 | yes | 19.1 / 16.2 | 47.6 / 53.1 | 0.019 / 0.034 | 29.7 / 24.7 |
| 10 | `Bayes-M4c-geo` | 28.82 | yes | 15.7 / 13.1 | 48.4 / 50.5 | 0.020 / 0.035 | 24.4 / 20.6 |
| 11 | `Bayes-M13-fill-geo` | 26.52 | yes | 14.1 / 12.4 | 46.5 / 51.0 | 0.020 / 0.036 | 28.2 / 22.9 |
| 12 | `Bayes-M17-fill-med` | 35.29 | yes | 18.9 / 16.4 | 55.2 / 62.4 | 0.020 / 0.035 | 29.6 / 25.0 |
| 13 | `Bayes-M4c-med` | 27.68 | yes | 15.3 / 12.4 | 53.3 / 56.4 | 0.020 / 0.036 | 25.7 / 20.7 |
| 14 | `Bayes-M13-fill-med` | 31.34 | yes | 16.6 / 14.8 | 65.6 / 77.1 | 0.021 / 0.037 | 28.2 / 24.1 |
| 15 | `Bayes-M4-geo` | 21.38 | yes | 11.7 / 9.7 | 55.0 / 57.9 | 0.021 / 0.037 | 21.2 / 18.1 |
| 16 | `Bayes-M8-fill-geo` | 21.43 | yes | 12.0 / 9.4 | 55.5 / 57.5 | 0.021 / 0.038 | 25.6 / 19.8 |
| 17 | `Bayes-M8-fill-med` | 25.77 | **spiky** | 14.0 / 11.8 | 73.4 / 84.0 | 0.022 / 0.039 | 26.4 / 21.6 |
| 18 | `Bayes-M4-med` | 25.31 | **spiky** | 14.0 / 11.3 | 74.6 / 84.1 | 0.022 / 0.039 | 24.2 / 19.0 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 48.76 | yes | 26.7 / 22.1 | 41.8 / 38.6 | 0.017 / 0.030 | 32.7 / 27.6 |
| 2 | `Bayes-M14-fill-geo` | 48.36 | yes | 25.8 / 22.6 | 44.0 / 45.8 | 0.017 / 0.030 | 32.3 / 28.0 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
