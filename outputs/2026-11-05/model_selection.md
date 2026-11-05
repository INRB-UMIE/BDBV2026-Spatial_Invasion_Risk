# Model selection — spatiotemporal

_Generated 2026-11-05T21:28:52+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_05112026 (processed_at 2026-11-05T12:00:00) |
| Predictions computed on snapshot | LINELIST_05112026 (processed_at 2026-11-05T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo` (kernel M14-fill, form fixed; stage 1 `Bayes-M14-fill-geo`: composite 13, margin 6) — **differs from the current pick** |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 13 of 18; h2: 13 of 18 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 13, margin 7) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 13, margin 7) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 56%; `Bayes-M14-fill-geo` 21%; `Bayes-M10-fill-geo` 8% — kernel M14-fill 77%, M10-fill 8%, M16-fill 8%, M13c-fill 8% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 48.12 | yes | 25.7 / 22.4 | 46.3 / 47.6 | 0.018 / 0.031 | 31.6 / 27.6 |
| 2 | `Bayes-M10-fill-geo` | 48.20 | yes | 27.6 / 20.6 | 47.5 / 49.8 | 0.018 / 0.032 | 32.6 / 25.8 |
| 3 | `Bayes-M13c-fill-geo` | 38.94 | yes | 21.0 / 17.9 | 37.4 / 40.5 | 0.019 / 0.033 | 31.1 / 25.6 |
| 4 | `Bayes-M14-fill-med` | 43.02 | yes | 23.3 / 19.8 | 52.3 / 56.1 | 0.018 / 0.033 | 31.5 / 26.1 |
| 5 | `Bayes-M17-fill-geo` | 36.24 | yes | 19.6 / 16.6 | 40.8 / 44.6 | 0.019 / 0.035 | 30.4 / 24.8 |
| 6 | `Bayes-M16-fill-geo` | 37.39 | yes | 19.8 / 17.6 | 53.0 / 48.4 | 0.019 / 0.033 | 29.2 / 24.6 |
| 7 | `Bayes-M10-fill-med` | 42.60 | yes | 23.5 / 19.1 | 56.8 / 60.7 | 0.019 / 0.033 | 31.9 / 26.0 |
| 8 | `Bayes-M13c-fill-med` | 35.52 | yes | 19.3 / 16.3 | 45.2 / 51.5 | 0.019 / 0.035 | 29.8 / 24.8 |
| 9 | `Bayes-M16-fill-med` | 40.71 | yes | 22.9 / 17.8 | 58.5 / 56.2 | 0.019 / 0.035 | 30.7 / 25.3 |
| 10 | `Bayes-M4c-geo` | 29.08 | yes | 15.9 / 13.2 | 47.2 / 49.4 | 0.020 / 0.036 | 24.8 / 20.4 |
| 11 | `Bayes-M17-fill-med` | 35.56 | yes | 19.2 / 16.4 | 51.9 / 59.8 | 0.020 / 0.036 | 29.7 / 25.1 |
| 12 | `Bayes-M13-fill-geo` | 26.84 | yes | 14.4 / 12.5 | 45.1 / 50.2 | 0.020 / 0.037 | 28.1 / 22.9 |
| 13 | `Bayes-M4c-med` | 27.78 | yes | 15.4 / 12.4 | 51.4 / 55.0 | 0.021 / 0.037 | 25.7 / 20.7 |
| 14 | `Bayes-M13-fill-med` | 31.50 | yes | 16.7 / 14.8 | 55.8 / 70.2 | 0.021 / 0.038 | 28.3 / 24.1 |
| 15 | `Bayes-M4-geo` | 22.46 | yes | 12.5 / 10.0 | 54.5 / 57.7 | 0.021 / 0.038 | 21.6 / 17.9 |
| 16 | `Bayes-M8-fill-geo` | 21.34 | yes | 11.9 / 9.4 | 55.4 / 57.3 | 0.021 / 0.039 | 25.4 / 19.6 |
| 17 | `Bayes-M8-fill-med` | 25.91 | yes | 14.1 / 11.8 | 63.7 / 76.9 | 0.022 / 0.040 | 26.5 / 21.7 |
| 18 | `Bayes-M4-med` | 25.46 | yes | 14.1 / 11.4 | 64.9 / 77.1 | 0.022 / 0.040 | 24.3 / 19.0 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 48.53 | yes | 26.4 / 22.1 | 42.6 / 40.3 | 0.017 / 0.031 | 32.8 / 27.9 |
| 2 | `Bayes-M14-fill-geo` | 48.12 | yes | 25.7 / 22.4 | 46.3 / 47.6 | 0.018 / 0.031 | 31.6 / 27.6 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
