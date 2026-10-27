# Model selection — spatiotemporal

_Generated 2026-10-26T23:09:09+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_26102026 (processed_at 2026-10-26T12:00:00) |
| Predictions computed on snapshot | LINELIST_26102026 (processed_at 2026-10-26T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 6) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 11 of 16; h2: 11 of 16 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 13) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 11) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 71%; `Bayes-M10-fill-geo` 11%; `Bayes-M16-fill-geo` 8% — kernel M14-fill 76%, M10-fill 11%, M16-fill 10%, M13c-fill 4%, M17-fill 0%, M4c 0% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 58.64 | yes | 35.7 / 23.0 | 51.7 / 50.3 | 0.020 / 0.033 | 50.9 / 30.6 |
| 2 | `Bayes-M10-fill-geo` | 53.22 | yes | 32.0 / 21.2 | 52.8 / 52.6 | 0.020 / 0.033 | 45.2 / 29.3 |
| 3 | `Bayes-M13c-fill-geo` | 45.27 | yes | 27.0 / 18.3 | 42.9 / 43.4 | 0.021 / 0.035 | 44.8 / 28.8 |
| 4 | `Bayes-M16-fill-geo` | 44.29 | yes | 25.3 / 19.0 | 51.4 / 50.9 | 0.021 / 0.035 | 43.1 / 28.8 |
| 5 | `Bayes-M14-fill-med` | 52.01 | yes | 30.7 / 21.4 | 62.7 / 63.3 | 0.020 / 0.034 | 43.9 / 29.7 |
| 6 | `Bayes-M17-fill-geo` | 40.19 | yes | 23.0 / 17.2 | 47.6 / 48.0 | 0.021 / 0.036 | 42.7 / 27.7 |
| 7 | `Bayes-M10-fill-med` | 48.20 | yes | 28.0 / 20.2 | 67.4 / 68.2 | 0.020 / 0.034 | 42.9 / 28.6 |
| 8 | `Bayes-M16-fill-med` | 46.40 | yes | 26.3 / 20.1 | 61.0 / 62.2 | 0.021 / 0.036 | 40.4 / 28.4 |
| 9 | `Bayes-M13c-fill-med` | 38.72 | yes | 22.1 / 16.6 | 53.0 / 53.8 | 0.021 / 0.036 | 41.8 / 27.2 |
| 10 | `Bayes-M4c-geo` | 33.30 | yes | 19.5 / 13.8 | 51.8 / 51.1 | 0.022 / 0.037 | 37.2 / 23.7 |
| 11 | `Bayes-M17-fill-med` | 41.10 | yes | 23.2 / 17.9 | 63.2 / 63.7 | 0.021 / 0.037 | 40.0 / 27.7 |
| 12 | `Bayes-M13-fill-geo` | 29.12 | yes | 16.7 / 12.4 | 52.0 / 52.8 | 0.022 / 0.038 | 38.3 / 24.7 |
| 13 | `Bayes-M4c-med` | 30.21 | yes | 17.3 / 12.9 | 56.6 / 56.5 | 0.022 / 0.038 | 33.6 / 22.0 |
| 14 | `Bayes-M13-fill-med` | 38.28 | yes | 22.4 / 15.9 | 68.5 / 68.6 | 0.023 / 0.039 | 37.7 / 25.8 |
| 15 | `Bayes-M4-geo` | 22.52 | yes | 12.8 / 9.7 | 58.5 / 58.4 | 0.023 / 0.039 | 28.9 / 19.5 |
| 16 | `Bayes-M8-fill-geo` | 22.05 | yes | 12.8 / 9.3 | 58.4 / 58.7 | 0.023 / 0.039 | 31.3 / 20.9 |
| 17 | `Bayes-M4-med` | 29.54 | yes | 17.5 / 12.1 | 74.9 / 74.3 | 0.023 / 0.040 | 31.5 / 20.6 |
| 18 | `Bayes-M8-fill-med` | 27.48 | yes | 15.5 / 12.0 | 74.4 / 74.0 | 0.024 / 0.040 | 33.7 / 23.1 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 56.15 | yes | 32.9 / 23.3 | 48.9 / 48.1 | 0.019 / 0.031 | 49.9 / 31.5 |
| 2 | `Bayes-M14-fill-geo` | 58.64 | yes | 35.7 / 23.0 | 51.7 / 50.3 | 0.020 / 0.033 | 50.9 / 30.6 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
