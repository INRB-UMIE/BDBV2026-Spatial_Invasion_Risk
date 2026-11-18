# Model selection — spatiotemporal

_Generated 2026-11-18T11:11:45+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_18112026 (processed_at 2026-11-18T12:00:00) |
| Predictions computed on snapshot | LINELIST_18112026 (processed_at 2026-11-18T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 13, margin 9) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 14 of 19; h2: 14 of 19 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 9) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 9) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 66%; `Bayes-M10-fill-geo` 14%; `Bayes-M16-fill-geo` 12% — kernel M14-fill 69%, M16-fill 16%, M10-fill 14%, M13c-fill 2%, M17-fill 0% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 46.05 | yes | 24.5 / 21.6 | 43.2 / 42.3 | 0.018 / 0.031 | 35.6 / 26.8 |
| 2 | `Bayes-M10-fill-geo` | 42.11 | yes | 22.0 / 20.1 | 43.7 / 43.8 | 0.018 / 0.031 | 32.7 / 25.4 |
| 3 | `Bayes-M13c-fill-geo` | 36.11 | yes | 18.8 / 17.3 | 37.6 / 39.0 | 0.019 / 0.033 | 34.8 / 25.0 |
| 4 | `Bayes-M16-fill-geo` | 36.71 | yes | 18.3 / 18.5 | 43.2 / 42.5 | 0.019 / 0.032 | 30.1 / 25.3 |
| 5 | `Bayes-M14-fill-med` | 42.35 | yes | 22.1 / 20.3 | 53.4 / 54.8 | 0.019 / 0.032 | 34.7 / 26.8 |
| 6 | `Bayes-M17-fill-geo` | 34.14 | yes | 17.9 / 16.2 | 40.1 / 41.7 | 0.020 / 0.033 | 34.3 / 23.9 |
| 7 | `Bayes-M10-fill-med` | 40.45 | yes | 21.1 / 19.3 | 56.4 / 59.3 | 0.019 / 0.032 | 32.6 / 25.3 |
| 8 | `Bayes-M16-fill-med` | 40.07 | yes | 21.0 / 19.1 | 51.7 / 52.2 | 0.020 / 0.033 | 33.9 / 26.6 |
| 9 | `Bayes-M13c-fill-med` | 34.54 | yes | 17.8 / 16.7 | 49.4 / 52.2 | 0.020 / 0.034 | 33.6 / 25.1 |
| 10 | `Bayes-M4c-geo` | 28.61 | yes | 14.4 / 14.2 | 44.8 / 45.5 | 0.020 / 0.034 | 26.8 / 20.3 |
| 11 | `Bayes-M13-fill-geo` | 26.38 | yes | 13.7 / 12.7 | 44.3 / 46.3 | 0.020 / 0.035 | 31.0 / 21.8 |
| 12 | `Bayes-M17-fill-med` | 34.58 | yes | 17.7 / 16.9 | 57.4 / 61.7 | 0.020 / 0.034 | 31.9 / 25.2 |
| 13 | `Bayes-M4-geo` | 20.74 | yes | 10.8 / 10.0 | 50.9 / 51.5 | 0.021 / 0.035 | 20.5 / 16.1 |
| 14 | `Bayes-M4c-med` | 25.84 | yes | 13.3 / 12.5 | 51.8 / 54.1 | 0.021 / 0.035 | 24.8 / 18.9 |
| 15 | `Bayes-M8-fill-geo` | 20.29 | yes | 10.6 / 9.7 | 50.3 / 51.2 | 0.021 / 0.035 | 23.5 / 17.4 |
| 16 | `Bayes-M13-fill-med` | 34.82 | yes | 17.8 / 17.0 | 68.2 / 75.7 | 0.021 / 0.036 | 31.3 / 24.9 |
| 17 | `Bayes-M8-fill-med` | 25.49 | **spiky** | 13.1 / 12.4 | 73.2 / 80.8 | 0.022 / 0.037 | 27.0 / 20.2 |
| 18 | `Bayes-M4-med` | 25.08 | **spiky** | 13.2 / 11.9 | 74.6 / 81.8 | 0.021 / 0.036 | 23.4 / 17.7 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 46.89 | yes | 24.8 / 22.1 | 40.4 / 39.7 | 0.017 / 0.028 | 37.2 / 27.5 |
| 2 | `Bayes-M14-fill-geo` | 46.05 | yes | 24.5 / 21.6 | 43.2 / 42.3 | 0.018 / 0.031 | 35.6 / 26.8 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
