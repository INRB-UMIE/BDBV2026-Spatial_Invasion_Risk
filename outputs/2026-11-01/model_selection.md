# Model selection — spatiotemporal

_Generated 2026-11-01T10:45:40+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_01112026 (processed_at 2026-11-01T12:00:00) |
| Predictions computed on snapshot | LINELIST_01112026 (processed_at 2026-11-01T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 12) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 12 of 17; h2: 12 of 17 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 14) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 13) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 72%; `Bayes-M14-fill-geo` 14%; `Bayes-M16-fill-geo` 6% — kernel M14-fill 86%, M16-fill 7%, M10-fill 4%, M13c-fill 3% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 57.15 | yes | 32.7 / 24.4 | 45.1 / 54.5 | 0.018 / 0.031 | 45.0 / 31.1 |
| 2 | `Bayes-M10-fill-geo` | 51.56 | yes | 29.8 / 21.8 | 46.5 / 57.0 | 0.019 / 0.031 | 40.9 / 28.9 |
| 3 | `Bayes-M13c-fill-geo` | 40.59 | yes | 22.7 / 17.9 | 38.5 / 45.2 | 0.020 / 0.034 | 41.5 / 27.8 |
| 4 | `Bayes-M16-fill-geo` | 41.30 | yes | 22.6 / 18.7 | 45.9 / 55.6 | 0.020 / 0.033 | 39.0 / 27.8 |
| 5 | `Bayes-M14-fill-med` | 52.44 | yes | 30.2 / 22.2 | 52.5 / 65.4 | 0.019 / 0.032 | 44.2 / 29.8 |
| 6 | `Bayes-M17-fill-geo` | 38.61 | yes | 21.7 / 16.9 | 42.3 / 49.8 | 0.020 / 0.034 | 40.7 / 27.5 |
| 7 | `Bayes-M10-fill-med` | 49.76 | yes | 28.3 / 21.5 | 55.3 / 69.5 | 0.019 / 0.032 | 40.2 / 28.5 |
| 8 | `Bayes-M16-fill-med` | 47.24 | yes | 26.7 / 20.5 | 51.9 / 65.0 | 0.020 / 0.034 | 40.8 / 28.1 |
| 9 | `Bayes-M13c-fill-med` | 37.91 | yes | 21.4 / 16.5 | 46.6 / 56.0 | 0.020 / 0.034 | 38.8 / 26.3 |
| 10 | `Bayes-M4c-geo` | 31.59 | yes | 17.9 / 13.7 | 48.5 / 53.5 | 0.021 / 0.035 | 32.2 / 22.8 |
| 11 | `Bayes-M17-fill-med` | 39.68 | yes | 22.4 / 17.2 | 53.7 / 65.8 | 0.020 / 0.035 | 39.8 / 27.1 |
| 12 | `Bayes-M13-fill-geo` | 28.48 | yes | 16.2 / 12.3 | 46.5 / 54.6 | 0.021 / 0.036 | 37.7 / 24.5 |
| 13 | `Bayes-M4c-med` | 29.24 | yes | 16.4 / 12.9 | 51.5 / 58.3 | 0.021 / 0.036 | 30.9 / 21.7 |
| 14 | `Bayes-M13-fill-med` | 38.34 | yes | 22.2 / 16.2 | 58.1 / 70.5 | 0.022 / 0.037 | 38.5 / 25.7 |
| 15 | `Bayes-M4-geo` | 23.64 | yes | 13.6 / 10.0 | 54.6 / 60.4 | 0.022 / 0.037 | 29.1 / 19.7 |
| 16 | `Bayes-M8-fill-geo` | 22.24 | yes | 12.9 / 9.4 | 54.2 / 60.5 | 0.022 / 0.038 | 32.5 / 21.4 |
| 17 | `Bayes-M4-med` | 28.45 | yes | 16.5 / 12.0 | 65.9 / 76.2 | 0.022 / 0.038 | 30.8 / 20.5 |
| 18 | `Bayes-M8-fill-med` | 26.90 | yes | 14.9 / 12.0 | 65.4 / 75.8 | 0.023 / 0.039 | 32.8 / 23.2 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 55.22 | yes | 31.8 / 23.4 | 43.2 / 52.3 | 0.018 / 0.029 | 44.6 / 30.3 |
| 2 | `Bayes-M14-fill-geo` | 57.15 | yes | 32.7 / 24.4 | 45.1 / 54.5 | 0.018 / 0.031 | 45.0 / 31.1 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
