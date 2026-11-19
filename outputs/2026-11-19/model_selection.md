# Model selection — spatiotemporal

_Generated 2026-11-19T12:03:19+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_19112026 (processed_at 2026-11-19T12:00:00) |
| Predictions computed on snapshot | LINELIST_19112026 (processed_at 2026-11-19T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo` (kernel M14-fill, form fixed; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 7) — **differs from the current pick** |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 15 of 20; h2: 15 of 20 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 5) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 5) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 65%; `Bayes-M14-fill-geo` 13%; `Bayes-M16-fill-geo` 9% — kernel M14-fill 78%, M16-fill 10%, M10-fill 9%, M13c-fill 4% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 48.37 | yes | 25.8 / 22.6 | 44.0 / 43.6 | 0.017 / 0.029 | 32.0 / 27.7 |
| 2 | `Bayes-M10-fill-geo` | 49.21 | yes | 27.6 / 21.6 | 46.2 / 46.0 | 0.017 / 0.030 | 33.0 / 27.1 |
| 3 | `Bayes-M13c-fill-geo` | 39.35 | yes | 20.8 / 18.5 | 38.9 / 39.9 | 0.018 / 0.032 | 30.5 / 25.8 |
| 4 | `Bayes-M16-fill-geo` | 39.35 | yes | 20.2 / 19.1 | 50.3 / 43.9 | 0.017 / 0.031 | 29.5 / 25.6 |
| 5 | `Bayes-M14-fill-med` | 43.59 | yes | 23.2 / 20.4 | 51.5 / 52.6 | 0.017 / 0.031 | 31.5 / 26.5 |
| 6 | `Bayes-M17-fill-geo` | 36.98 | yes | 19.6 / 17.3 | 41.7 / 43.1 | 0.018 / 0.033 | 30.8 / 25.5 |
| 7 | `Bayes-M10-fill-med` | 43.36 | yes | 23.5 / 19.9 | 56.3 / 57.3 | 0.018 / 0.031 | 31.9 / 26.3 |
| 8 | `Bayes-M16-fill-med` | 42.18 | yes | 22.9 / 19.3 | 56.4 / 51.8 | 0.018 / 0.032 | 31.0 / 26.3 |
| 9 | `Bayes-M13c-fill-med` | 35.88 | yes | 19.0 / 16.9 | 47.6 / 50.5 | 0.018 / 0.033 | 29.7 / 25.1 |
| 10 | `Bayes-M4c-geo` | 29.86 | yes | 15.9 / 13.9 | 48.4 / 48.0 | 0.019 / 0.034 | 24.9 / 21.8 |
| 11 | `Bayes-M13-fill-geo` | 26.72 | yes | 14.2 / 12.5 | 46.5 / 48.6 | 0.019 / 0.035 | 28.1 / 22.8 |
| 12 | `Bayes-M17-fill-med` | 35.98 | yes | 18.9 / 17.0 | 55.2 / 59.2 | 0.019 / 0.034 | 29.8 / 25.3 |
| 13 | `Bayes-M4c-med` | 27.94 | yes | 15.3 / 12.7 | 53.4 / 53.6 | 0.019 / 0.035 | 26.2 / 20.9 |
| 14 | `Bayes-M13-fill-med` | 32.74 | yes | 16.6 / 16.2 | 65.6 / 73.1 | 0.020 / 0.036 | 28.2 / 24.8 |
| 15 | `Bayes-M4-geo` | 22.36 | yes | 12.5 / 9.9 | 55.2 / 55.6 | 0.020 / 0.036 | 21.6 / 18.4 |
| 16 | `Bayes-M8-fill-geo` | 21.46 | yes | 12.0 / 9.5 | 55.5 / 55.1 | 0.020 / 0.036 | 25.4 / 20.2 |
| 17 | `Bayes-M8-fill-med` | 25.96 | **spiky** | 13.9 / 12.0 | 73.3 / 79.8 | 0.021 / 0.038 | 26.5 / 21.8 |
| 18 | `Bayes-M4-med` | 25.63 | **spiky** | 14.0 / 11.7 | 74.5 / 79.9 | 0.021 / 0.037 | 24.3 / 19.3 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 49.01 | yes | 26.4 / 22.6 | 41.3 / 38.1 | 0.016 / 0.029 | 32.3 / 27.7 |
| 2 | `Bayes-M14-fill-geo` | 48.37 | yes | 25.8 / 22.6 | 44.0 / 43.6 | 0.017 / 0.029 | 32.0 / 27.7 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
