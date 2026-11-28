# Model selection — spatiotemporal

_Generated 2026-11-28T14:23:05+0000_

**Featured Bayesian model (two-stage CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo-tvrw1`
**Stage 1 (fixed-beta specification):** `Bayes-M14-fill-geo`; **stage 2** compared it with its time-varying forms `Bayes-M14-fill-geo-tvrw1`; the frailty forms `Bayes-M14-fill-geo-tvfrailty`, `Bayes-M14-fill-geo-tvrw1frailty` are a secondary analysis, scored but never candidates
**Best renewal model:** _none cross-validated in this evaluation_
**Headline model:** `Bayes-M14-fill-geo-tvrw1` (the featured model; ensembles and rank-only comparators are not candidates)

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_27112026 (processed_at 2026-11-27T12:00:00) |
| Predictions computed on snapshot | LINELIST_27112026 (processed_at 2026-11-27T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | All three axes are computed on the rounds AFTER the ascertainment burn-in (cutoff week index >= 7), and AUC-PR skill is scored WITHIN each forecast date (per-date AUC-PR) rather than on all rounds pooled. Excluding the burn-in is a deliberate choice: in the discovery phase where first cases were confirmed partly reflected where surveillance looked first, so those rounds are kept out of the choice of model. Performance is REPORTED over all rounds. Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Selection basis | per_date_post_burnin |
| Pick under the PREVIOUS rule (pooled AUC-PR, all rounds) | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 15, margin 1) — same pick |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Post-burn-in rounds with a fitted delta (featured model; the rest score at delta = 1) | h1: 16 of 21; h2: 16 of 21 |
| Pick under RAW axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 10) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo-tvrw1` (kernel M14-fill, form rw1; stage 1 `Bayes-M14-fill-geo`: composite 10, margin 10) |
| Axis sensitivity | none — the same model wins on both axes |
| Selection stability (zone bootstrap) | `Bayes-M14-fill-geo-tvrw1` 66%; `Bayes-M16-fill-geo` 14%; `Bayes-M14-fill-geo` 13% — kernel M14-fill 80%, M16-fill 14%, M10-fill 3%, M13c-fill 2% |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |
| Cross-check (renewal) | no pipeline pick to compare (family not featured) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo-tvrw1` |

### Stage 1: fixed-beta Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo` | 52.34 | yes | 29.8 / 22.5 | 37.5 / 42.1 | 0.016 / 0.028 | 41.3 / 28.0 |
| 2 | `Bayes-M10-fill-geo` | 49.95 | yes | 29.0 / 20.9 | 39.5 / 44.2 | 0.017 / 0.029 | 41.9 / 27.1 |
| 3 | `Bayes-M16-fill-geo` | 48.32 | yes | 27.3 / 21.1 | 43.6 / 43.3 | 0.017 / 0.030 | 42.4 / 27.9 |
| 4 | `Bayes-M13c-fill-geo` | 43.70 | yes | 25.1 / 18.6 | 33.3 / 38.0 | 0.017 / 0.031 | 44.2 / 27.2 |
| 5 | `Bayes-M14-fill-med` | 49.02 | yes | 27.9 / 21.1 | 43.8 / 50.5 | 0.017 / 0.030 | 42.5 / 27.9 |
| 6 | `Bayes-M17-fill-geo` | 39.48 | yes | 22.1 / 17.3 | 36.3 / 41.2 | 0.018 / 0.032 | 40.3 / 25.3 |
| 7 | `Bayes-M16-fill-med` | 49.65 | yes | 28.1 / 21.5 | 49.5 / 50.9 | 0.018 / 0.031 | 42.1 / 28.1 |
| 8 | `Bayes-M10-fill-med` | 46.11 | yes | 26.4 / 19.7 | 48.1 / 55.5 | 0.017 / 0.030 | 39.7 / 25.8 |
| 9 | `Bayes-M13c-fill-med` | 40.19 | yes | 23.1 / 17.1 | 41.1 / 48.4 | 0.018 / 0.032 | 41.3 / 25.2 |
| 10 | `Bayes-M4c-geo` | 34.53 | yes | 20.0 / 14.5 | 40.3 / 45.5 | 0.019 / 0.033 | 33.4 / 20.9 |
| 11 | `Bayes-M17-fill-med` | 40.26 | yes | 22.8 / 17.4 | 47.9 / 56.7 | 0.018 / 0.033 | 39.9 / 25.5 |
| 12 | `Bayes-M13-fill-geo` | 29.24 | yes | 16.0 / 13.2 | 41.3 / 46.3 | 0.019 / 0.034 | 35.6 / 22.5 |
| 13 | `Bayes-M4c-med` | 30.10 | yes | 17.4 / 12.7 | 45.4 / 51.6 | 0.019 / 0.034 | 31.0 / 19.3 |
| 14 | `Bayes-M13-fill-med` | 40.69 | yes | 22.8 / 17.9 | 57.5 / 70.0 | 0.020 / 0.035 | 38.8 / 25.4 |
| 15 | `Bayes-M4-geo` | 25.58 | yes | 14.4 / 11.2 | 49.0 / 52.9 | 0.020 / 0.035 | 27.1 / 18.0 |
| 16 | `Bayes-M8-fill-geo` | 23.03 | yes | 13.2 / 9.8 | 49.5 / 52.8 | 0.020 / 0.035 | 31.6 / 17.9 |
| 17 | `Bayes-M4-med` | 29.53 | **spiky** | 17.0 / 12.5 | 65.6 / 76.4 | 0.020 / 0.036 | 28.7 / 18.4 |
| 18 | `Bayes-M8-fill-med` | 28.34 | **spiky** | 16.6 / 11.8 | 64.6 / 76.4 | 0.020 / 0.037 | 32.5 / 19.2 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

### Stage 2: the stage-1 specification and its time-varying forms leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log score, recalibrated (h1/2) | Per-date AUC-PR skill, ALL rounds — as reported (h1/2) |
|---|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvrw1` | 51.60 | yes | 29.6 / 22.0 | 35.7 / 38.7 | 0.016 / 0.027 | 44.6 / 28.9 |
| 2 | `Bayes-M14-fill-geo` | 52.34 | yes | 29.8 / 22.5 | 37.5 / 42.1 | 0.016 / 0.028 | 41.3 / 28.0 |

_AUC-PR skill, rank-of-truth and log score are the SELECTION axes: per-date AUC-PR, post-burn-in rounds (the burn-in is kept out of the choice deliberately). The last column is the same per-date AUC-PR skill over ALL rounds, the value the figures report._

_(no Best-renewal models in the evaluation)_
