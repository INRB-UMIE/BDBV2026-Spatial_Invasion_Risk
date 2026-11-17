# Model selection — spatiotemporal

_Generated 2026-11-17T10:56:33+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_17112026 (processed_at 2026-11-17T12:00:00) |
| Predictions computed on snapshot | LINELIST_17112026 (processed_at 2026-11-17T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 13, margin 1) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 14 of 19; h2: 14 of 19 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 7) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 11, margin 7) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 69%; `Bayes-M10-fill-geo` 12%; `Bayes-M16-fill-geo` 9% — kernel M14-fill 75%, M10-fill 12%, M16-fill 11%, M13c-fill 2% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 57.68 | yes | 33.7 / 24.0 | 49.2 / 46.6 | 0.017 / 0.029 | 49.9 / 31.5 |
| 2 | `Bayes-M10-fill-geo` | 54.76 | yes | 32.0 / 22.8 | 50.4 / 49.1 | 0.017 / 0.029 | 46.0 / 29.7 |
| 3 | `Bayes-M13c-fill-geo` | 46.11 | yes | 26.8 / 19.3 | 43.3 / 42.7 | 0.018 / 0.031 | 49.6 / 29.4 |
| 4 | `Bayes-M16-fill-geo` | 42.29 | yes | 23.1 / 19.2 | 47.5 / 45.8 | 0.018 / 0.031 | 43.0 / 28.7 |
| 5 | `Bayes-M14-fill-med` | 54.07 | yes | 31.4 / 22.7 | 61.1 / 60.5 | 0.018 / 0.030 | 48.8 / 30.4 |
| 6 | `Bayes-M17-fill-geo` | 40.40 | yes | 23.3 / 17.1 | 47.0 / 46.7 | 0.019 / 0.032 | 47.1 / 27.5 |
| 7 | `Bayes-M10-fill-med` | 51.57 | yes | 29.7 / 21.9 | 65.7 / 65.6 | 0.018 / 0.030 | 44.1 / 28.8 |
| 8 | `Bayes-M16-fill-med` | 46.29 | yes | 26.1 / 20.2 | 57.0 / 57.2 | 0.018 / 0.032 | 45.0 / 29.4 |
| 9 | `Bayes-M13c-fill-med` | 39.77 | yes | 22.5 / 17.3 | 55.0 / 55.8 | 0.019 / 0.033 | 45.5 / 28.3 |
| 10 | `Bayes-M4c-geo` | 35.53 | yes | 20.2 / 15.4 | 51.8 / 49.8 | 0.019 / 0.033 | 33.8 / 23.3 |
| 11 | `Bayes-M13-fill-geo` | 29.76 | yes | 17.0 / 12.8 | 52.1 / 51.6 | 0.019 / 0.034 | 42.7 / 24.3 |
| 12 | `Bayes-M17-fill-med` | 40.22 | yes | 22.4 / 17.8 | 65.3 / 66.1 | 0.019 / 0.033 | 43.1 / 28.2 |
| 13 | `Bayes-M4c-med` | 31.15 | yes | 17.7 / 13.5 | 58.5 / 58.0 | 0.020 / 0.034 | 33.4 / 21.7 |
| 14 | `Bayes-M4-geo` | 23.79 | yes | 13.3 / 10.4 | 57.8 / 56.6 | 0.020 / 0.034 | 26.0 / 19.3 |
| 15 | `Bayes-M13-fill-med` | 39.11 | yes | 22.5 / 16.6 | 77.1 / 78.2 | 0.020 / 0.035 | 42.8 / 26.8 |
| 16 | `Bayes-M8-fill-geo` | 22.91 | yes | 13.1 / 9.9 | 58.0 / 57.1 | 0.020 / 0.034 | 32.8 / 20.5 |
| 17 | `Bayes-M4-med` | 30.62 | yes | 17.9 / 12.7 | 83.5 / 84.0 | 0.020 / 0.035 | 32.5 / 21.0 |
| 18 | `Bayes-M8-fill-med` | 28.79 | yes | 16.1 / 12.7 | 83.6 / 84.0 | 0.021 / 0.036 | 37.3 / 22.9 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 58.06 | yes | 34.0 / 24.0 | 46.5 / 44.8 | 0.016 / 0.028 | 50.1 / 31.7 |
| 2 | `Bayes-M14-fill-geo` | 57.68 | yes | 33.7 / 24.0 | 49.2 / 46.6 | 0.017 / 0.029 | 49.9 / 31.5 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
