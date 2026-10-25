# Model selection — spatiotemporal

_Generated 2026-10-25T21:49:57+0000_

**Featured Bayesian model (lowest non-spiky CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M14-fill-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_24102026 (processed_at 2026-10-24T12:00:00) |
| Predictions computed on snapshot | LINELIST_24102026 (processed_at 2026-10-24T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 10, margin 2) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 11, margin 3) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 85.86 | yes | 51.5 / 34.4 | 26.3 / 32.2 | 0.019 / 0.032 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 85.46 | yes | 51.5 / 34.0 | 26.1 / 31.9 | 0.019 / 0.032 |
| 3 | `Bayes-M14-fill-geo-gtshort` | 79.42 | yes | 46.5 / 33.0 | 26.7 / 32.7 | 0.020 / 0.035 |
| 4 | `Bayes-M14-fill-geo` | 73.84 | yes | 43.0 / 30.8 | 26.9 / 33.0 | 0.021 / 0.035 |
| 5 | `Bayes-M14-fill-geo-gtlong` | 71.34 | yes | 41.3 / 30.0 | 27.0 / 32.8 | 0.021 / 0.036 |
| 6 | `Bayes-M13c-fill-geo` | 71.33 | yes | 43.1 / 28.2 | 25.5 / 30.4 | 0.022 / 0.039 |
| 7 | `Bayes-M10-fill-geo` | 74.85 | yes | 44.1 / 30.7 | 28.6 / 35.2 | 0.021 / 0.038 |
| 8 | `Bayes-M16-fill-geo` | 72.15 | yes | 43.5 / 28.7 | 28.4 / 32.7 | 0.022 / 0.037 |
| 9 | `Bayes-M10-fill-med` | 79.81 | yes | 47.1 / 32.7 | 33.7 / 43.0 | 0.022 / 0.038 |
| 10 | `Bayes-M17-fill-geo` | 71.66 | yes | 43.0 / 28.6 | 27.3 / 33.3 | 0.023 / 0.040 |
| 11 | `Bayes-M14-fill-med` | 70.77 | yes | 40.8 / 30.0 | 32.0 / 40.2 | 0.021 / 0.036 |
| 12 | `Bayes-M17-fill-med` | 75.81 | yes | 45.0 / 30.8 | 33.2 / 41.5 | 0.023 / 0.041 |
| 13 | `Bayes-M13c-fill-med` | 69.55 | yes | 40.8 / 28.7 | 29.3 / 35.6 | 0.022 / 0.039 |
| 14 | `Bayes-M16-fill-med` | 70.79 | yes | 41.3 / 29.5 | 32.0 / 38.0 | 0.022 / 0.038 |
| 15 | `Bayes-M13-fill-geo` | 69.59 | yes | 42.8 / 26.7 | 31.3 / 37.8 | 0.024 / 0.043 |
| 16 | `Bayes-M13-fill-med` | 74.87 | yes | 45.2 / 29.7 | 36.3 / 45.3 | 0.025 / 0.046 |
| 17 | `Bayes-M4c-geo` | 58.65 | yes | 35.3 / 23.4 | 33.4 / 37.9 | 0.024 / 0.043 |
| 18 | `Bayes-M8-fill-med` | 73.44 | yes | 45.9 / 27.5 | 43.5 / 51.5 | 0.027 / 0.052 |
| 19 | `Bayes-M4c-med` | 59.29 | yes | 35.0 / 24.3 | 35.8 / 40.8 | 0.024 / 0.043 |
| 20 | `Bayes-M8-fill-geo` | 63.68 | yes | 41.0 / 22.7 | 38.6 / 44.5 | 0.026 / 0.048 |
| 21 | `Bayes-M4-med` | 63.45 | yes | 42.0 / 21.5 | 46.1 / 52.5 | 0.027 / 0.052 |
| 22 | `Bayes-M4-geo` | 54.35 | yes | 33.5 / 20.9 | 40.1 / 44.6 | 0.026 / 0.049 |

_(no Best-renewal models in the evaluation)_
