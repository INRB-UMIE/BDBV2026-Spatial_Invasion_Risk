# Model selection — spatiotemporal

_Generated 2026-10-20T03:21:32+0000_

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
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 9, margin 4) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 14, margin 3) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 78.32 | yes | 45.9 / 32.4 | 27.7 / 33.7 | 0.019 / 0.034 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 78.38 | yes | 46.1 / 32.3 | 27.2 / 33.6 | 0.019 / 0.033 |
| 3 | `Bayes-M14-fill-geo-gtshort` | 76.67 | yes | 44.7 / 31.9 | 28.5 / 34.6 | 0.019 / 0.035 |
| 4 | `Bayes-M14-fill-geo` | 75.23 | yes | 43.3 / 31.9 | 28.4 / 34.7 | 0.020 / 0.035 |
| 5 | `Bayes-M14-fill-geo-gtlong` | 72.96 | yes | 41.3 / 31.6 | 28.9 / 35.1 | 0.020 / 0.035 |
| 6 | `Bayes-M13c-fill-geo` | 61.39 | yes | 34.6 / 26.8 | 26.4 / 31.8 | 0.021 / 0.039 |
| 7 | `Bayes-M10-fill-geo` | 68.58 | yes | 39.8 / 28.8 | 31.2 / 37.8 | 0.021 / 0.039 |
| 8 | `Bayes-M16-fill-geo` | 64.37 | yes | 35.7 / 28.7 | 31.3 / 36.8 | 0.020 / 0.037 |
| 9 | `Bayes-M14-fill-med` | 68.45 | yes | 40.2 / 28.2 | 34.5 / 43.6 | 0.020 / 0.037 |
| 10 | `Bayes-M17-fill-geo` | 62.68 | yes | 35.6 / 27.1 | 29.3 / 35.2 | 0.022 / 0.041 |
| 11 | `Bayes-M16-fill-med` | 63.24 | yes | 36.1 / 27.2 | 36.2 / 44.5 | 0.021 / 0.038 |
| 12 | `Bayes-M13c-fill-med` | 58.79 | yes | 33.6 / 25.2 | 31.6 / 38.0 | 0.022 / 0.040 |
| 13 | `Bayes-M10-fill-med` | 63.24 | yes | 37.4 / 25.9 | 38.3 / 48.6 | 0.021 / 0.041 |
| 14 | `Bayes-M17-fill-med` | 62.71 | yes | 37.0 / 25.7 | 38.2 / 47.1 | 0.022 / 0.043 |
| 15 | `Bayes-M13-fill-geo` | 55.83 | yes | 31.6 / 24.2 | 33.4 / 39.5 | 0.023 / 0.044 |
| 16 | `Bayes-M4c-geo` | 47.72 | yes | 28.5 / 19.2 | 34.0 / 39.2 | 0.024 / 0.047 |
| 17 | `Bayes-M13-fill-med` | 59.73 | yes | 35.8 / 24.0 | 42.2 / 51.5 | 0.024 / 0.046 |
| 18 | `Bayes-M8-fill-geo` | 49.05 | yes | 28.8 / 20.3 | 39.8 / 45.4 | 0.025 / 0.050 |
| 19 | `Bayes-M4c-med` | 45.08 | yes | 27.8 / 17.2 | 36.8 / 42.5 | 0.024 / 0.049 |
| 20 | `Bayes-M4-geo` | 40.00 | yes | 23.9 / 16.1 | 40.3 / 45.8 | 0.026 / 0.052 |
| 21 | `Bayes-M8-fill-med` | 39.18 | yes | 24.5 / 14.7 | 48.6 / 57.4 | 0.028 / 0.061 |
| 22 | `Bayes-M4-med` | 32.96 | yes | 21.0 / 11.9 | 49.6 / 58.3 | 0.028 / 0.064 |

_(no Best-renewal models in the evaluation)_
