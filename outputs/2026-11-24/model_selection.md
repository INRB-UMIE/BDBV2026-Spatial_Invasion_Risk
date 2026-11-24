# Model selection — spatiotemporal

_Generated 2026-11-24T09:52:14+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_24112026 (processed_at 2026-11-24T12:00:00) |
| Predictions computed on snapshot | LINELIST_24112026 (processed_at 2026-11-24T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 12, margin 5) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 15 of 20; h2: 15 of 20 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 9) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 9) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 72%; `Bayes-M16-fill-geo` 14%; `Bayes-M10-fill-geo` 7% — kernel M14-fill 77%, M16-fill 16%, M10-fill 7%, M13c-fill 1% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 59.26 | yes | 35.0 / 24.2 | 46.1 / 44.5 | 0.017 / 0.028 | 50.6 / 31.3 |
| 2 | `Bayes-M10-fill-geo` | 55.13 | yes | 32.4 / 22.7 | 47.4 / 46.9 | 0.017 / 0.028 | 46.7 / 29.9 |
| 3 | `Bayes-M16-fill-geo` | 49.06 | yes | 28.3 / 20.7 | 44.4 / 43.5 | 0.017 / 0.030 | 46.3 / 29.4 |
| 4 | `Bayes-M13c-fill-geo` | 47.92 | yes | 28.2 / 19.7 | 41.0 / 40.8 | 0.018 / 0.031 | 50.4 / 29.9 |
| 5 | `Bayes-M14-fill-med` | 57.08 | yes | 33.6 / 23.4 | 57.3 / 57.5 | 0.017 / 0.029 | 50.2 / 31.4 |
| 6 | `Bayes-M16-fill-med` | 52.83 | yes | 31.1 / 21.7 | 53.6 / 54.2 | 0.018 / 0.031 | 48.3 / 30.4 |
| 7 | `Bayes-M17-fill-geo` | 41.89 | yes | 24.3 / 17.5 | 44.4 / 44.6 | 0.018 / 0.031 | 47.6 / 27.8 |
| 8 | `Bayes-M10-fill-med` | 52.80 | yes | 30.7 / 22.1 | 61.8 / 62.4 | 0.017 / 0.029 | 44.6 / 29.0 |
| 9 | `Bayes-M13c-fill-med` | 42.91 | yes | 24.9 / 18.0 | 51.7 / 52.9 | 0.018 / 0.032 | 47.2 / 28.8 |
| 10 | `Bayes-M4c-geo` | 35.89 | yes | 20.5 / 15.4 | 48.9 / 47.7 | 0.019 / 0.032 | 34.5 / 23.6 |
| 11 | `Bayes-M17-fill-med` | 43.35 | yes | 24.8 / 18.5 | 61.4 / 62.7 | 0.019 / 0.032 | 44.6 / 28.7 |
| 12 | `Bayes-M13-fill-geo` | 30.44 | yes | 17.6 / 12.9 | 49.0 / 49.4 | 0.019 / 0.033 | 43.2 / 24.6 |
| 13 | `Bayes-M4c-med` | 31.58 | yes | 18.1 / 13.5 | 55.4 / 55.4 | 0.019 / 0.033 | 33.5 / 21.8 |
| 14 | `Bayes-M13-fill-med` | 45.69 | yes | 27.6 / 18.1 | 72.4 / 74.2 | 0.020 / 0.034 | 46.1 / 27.9 |
| 15 | `Bayes-M4-geo` | 24.00 | yes | 13.6 / 10.4 | 55.3 / 54.5 | 0.019 / 0.033 | 26.9 / 19.8 |
| 16 | `Bayes-M8-fill-geo` | 23.63 | yes | 13.6 / 10.0 | 55.3 / 54.8 | 0.020 / 0.033 | 33.1 / 20.7 |
| 17 | `Bayes-M4-med` | 31.96 | yes | 19.0 / 13.0 | 78.3 / 79.8 | 0.020 / 0.034 | 33.1 / 21.2 |
| 18 | `Bayes-M8-fill-med` | 29.72 | yes | 16.8 / 12.9 | 78.7 / 79.9 | 0.020 / 0.035 | 37.6 / 23.3 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 59.93 | yes | 35.5 / 24.5 | 43.8 / 42.8 | 0.016 / 0.027 | 50.9 / 32.0 |
| 2 | `Bayes-M14-fill-geo` | 59.26 | yes | 35.0 / 24.2 | 46.1 / 44.5 | 0.017 / 0.028 | 50.6 / 31.3 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
