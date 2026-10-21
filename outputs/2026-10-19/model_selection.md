# Model selection — spatiotemporal

_Generated 2026-10-21T01:25:13+0000_

**Featured Bayesian model (lowest non-spiky CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M14-fill-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_19102026 (processed_at 2026-10-19T12:00:00) |
| Predictions computed on snapshot | LINELIST_19102026 (processed_at 2026-10-19T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 9, margin 2) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 10, margin 3) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 78.00 | yes | 45.9 / 32.1 | 19.4 / 26.6 | 0.018 / 0.032 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 78.88 | yes | 46.8 / 32.1 | 19.5 / 26.5 | 0.018 / 0.031 |
| 3 | `Bayes-M14-fill-geo-gtlong` | 72.14 | yes | 41.5 / 30.6 | 19.9 / 26.9 | 0.019 / 0.033 |
| 4 | `Bayes-M14-fill-geo` | 74.04 | yes | 43.3 / 30.7 | 19.9 / 27.0 | 0.019 / 0.034 |
| 5 | `Bayes-M14-fill-geo-gtshort` | 75.12 | yes | 44.6 / 30.5 | 19.7 / 26.6 | 0.019 / 0.033 |
| 6 | `Bayes-M14-fill-med` | 69.81 | yes | 41.3 / 28.5 | 23.0 / 32.7 | 0.019 / 0.034 |
| 7 | `Bayes-M16-fill-geo` | 62.65 | yes | 35.0 / 27.7 | 22.1 / 27.0 | 0.019 / 0.034 |
| 8 | `Bayes-M10-fill-geo` | 68.28 | yes | 40.0 / 28.3 | 21.7 / 28.6 | 0.020 / 0.036 |
| 9 | `Bayes-M13c-fill-geo` | 59.48 | yes | 33.4 / 26.1 | 20.8 / 26.1 | 0.020 / 0.037 |
| 10 | `Bayes-M16-fill-med` | 64.44 | yes | 36.7 / 27.8 | 25.5 / 33.0 | 0.019 / 0.035 |
| 11 | `Bayes-M17-fill-geo` | 60.94 | yes | 34.1 / 26.9 | 23.2 / 29.1 | 0.021 / 0.038 |
| 12 | `Bayes-M17-fill-med` | 63.55 | yes | 36.9 / 26.7 | 25.6 / 33.6 | 0.021 / 0.039 |
| 13 | `Bayes-M13c-fill-med` | 59.08 | yes | 33.7 / 25.3 | 24.1 / 30.5 | 0.021 / 0.037 |
| 14 | `Bayes-M10-fill-med` | 62.24 | yes | 36.8 / 25.5 | 24.4 / 33.9 | 0.020 / 0.038 |
| 15 | `Bayes-M13-fill-geo` | 54.58 | yes | 30.5 / 24.0 | 27.7 / 34.0 | 0.022 / 0.041 |
| 16 | `Bayes-M13-fill-med` | 58.63 | yes | 34.9 / 23.8 | 29.6 / 37.3 | 0.023 / 0.043 |
| 17 | `Bayes-M4c-geo` | 46.51 | yes | 28.0 / 18.5 | 29.2 / 33.8 | 0.023 / 0.045 |
| 18 | `Bayes-M8-fill-geo` | 43.78 | yes | 24.4 / 19.4 | 36.0 / 41.4 | 0.024 / 0.046 |
| 19 | `Bayes-M4c-med` | 45.02 | yes | 27.5 / 17.5 | 30.5 / 35.8 | 0.023 / 0.045 |
| 20 | `Bayes-M4-geo` | 36.93 | yes | 21.5 / 15.4 | 36.3 / 41.6 | 0.025 / 0.049 |
| 21 | `Bayes-M8-fill-med` | 40.39 | yes | 25.1 / 15.3 | 37.7 / 44.0 | 0.026 / 0.053 |
| 22 | `Bayes-M4-med` | 31.83 | yes | 19.9 / 11.9 | 38.3 / 44.6 | 0.027 / 0.057 |

_(no Best-renewal models in the evaluation)_
