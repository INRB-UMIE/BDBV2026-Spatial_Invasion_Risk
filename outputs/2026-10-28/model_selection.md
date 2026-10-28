# Model selection — spatiotemporal

_Generated 2026-10-28T10:28:58+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_28102026 (processed_at 2026-10-28T12:00:00) |
| Predictions computed on snapshot | LINELIST_28102026 (processed_at 2026-10-28T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo` (kernel M14-fill, form fixed; stage 1 `Bayes-M14-fill-geo`: composite 12, margin 16) — **differs from the current pick** |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 11 of 16; h2: 11 of 16 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 13, margin 9) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 14, margin 7) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 57%; `Bayes-M10-fill-geo` 20%; `Bayes-M16-fill-geo` 10% — kernel M14-fill 61%, M10-fill 20%, M16-fill 13%, M13c-fill 6%, M13-fill 0% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 47.06 | yes | 25.7 / 21.4 | 45.0 / 48.1 | 0.020 / 0.033 | 38.7 / 26.6 |
| 2 | `Bayes-M10-fill-geo` | 43.94 | yes | 23.9 / 20.0 | 45.6 / 49.6 | 0.020 / 0.033 | 35.8 / 24.4 |
| 3 | `Bayes-M13c-fill-geo` | 36.83 | yes | 19.7 / 17.1 | 36.9 / 41.0 | 0.021 / 0.035 | 36.7 / 24.4 |
| 4 | `Bayes-M17-fill-geo` | 34.88 | yes | 18.6 / 16.3 | 40.6 / 44.6 | 0.021 / 0.036 | 36.2 / 24.4 |
| 5 | `Bayes-M16-fill-geo` | 36.41 | yes | 18.9 / 17.5 | 45.2 / 48.9 | 0.021 / 0.035 | 32.6 / 24.9 |
| 6 | `Bayes-M14-fill-med` | 42.62 | yes | 22.7 / 20.0 | 54.9 / 60.0 | 0.020 / 0.034 | 36.4 / 26.2 |
| 7 | `Bayes-M10-fill-med` | 41.23 | yes | 22.2 / 19.0 | 56.9 / 64.0 | 0.020 / 0.034 | 35.0 / 25.5 |
| 8 | `Bayes-M4c-geo` | 29.50 | yes | 15.5 / 14.0 | 43.9 / 47.7 | 0.021 / 0.036 | 28.2 / 19.6 |
| 9 | `Bayes-M13c-fill-med` | 34.60 | yes | 18.4 / 16.2 | 46.8 / 52.0 | 0.021 / 0.036 | 35.6 / 25.1 |
| 10 | `Bayes-M16-fill-med` | 39.56 | yes | 21.7 / 17.9 | 53.6 / 58.8 | 0.021 / 0.036 | 35.4 / 26.1 |
| 11 | `Bayes-M13-fill-geo` | 27.35 | yes | 14.4 / 12.9 | 44.2 / 48.8 | 0.022 / 0.037 | 32.9 / 21.9 |
| 12 | `Bayes-M17-fill-med` | 34.65 | yes | 18.1 / 16.6 | 54.2 / 60.5 | 0.022 / 0.036 | 33.4 / 25.5 |
| 13 | `Bayes-M4c-med` | 27.16 | yes | 14.3 / 12.8 | 49.3 / 54.1 | 0.022 / 0.036 | 26.4 / 19.1 |
| 14 | `Bayes-M4-geo` | 21.45 | yes | 11.3 / 10.2 | 50.9 / 54.7 | 0.022 / 0.037 | 22.1 / 15.6 |
| 15 | `Bayes-M13-fill-med` | 34.26 | yes | 18.3 / 15.9 | 58.3 / 65.1 | 0.023 / 0.038 | 33.1 / 24.5 |
| 16 | `Bayes-M8-fill-geo` | 20.87 | yes | 11.1 / 9.8 | 50.0 / 54.3 | 0.022 / 0.037 | 25.5 / 17.1 |
| 17 | `Bayes-M8-fill-med` | 26.22 | yes | 13.7 / 12.5 | 63.4 / 70.7 | 0.023 / 0.039 | 28.7 / 20.5 |
| 18 | `Bayes-M4-med` | 26.21 | yes | 14.2 / 12.0 | 64.3 / 71.4 | 0.023 / 0.038 | 24.9 / 17.8 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 47.37 | yes | 25.4 / 22.0 | 43.4 / 46.5 | 0.019 / 0.030 | 38.5 / 27.5 |
| 2 | `Bayes-M14-fill-geo` | 47.06 | yes | 25.7 / 21.4 | 45.0 / 48.1 | 0.020 / 0.033 | 38.7 / 26.6 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
