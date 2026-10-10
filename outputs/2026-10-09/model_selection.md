# Model selection — spatiotemporal

_Generated 2026-10-10T02:08:47+0000_

**Featured Bayesian model (lowest non-spiky CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M14-fill-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_09102026 (processed_at 2026-10-09T12:00:00) |
| Predictions computed on snapshot | LINELIST_09102026 (processed_at 2026-10-09T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 9, margin 3) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 13, margin 4) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 64.35 | yes | 39.4 / 24.9 | 23.2 / 23.6 | 0.020 / 0.037 |
| 2 | `Bayes-M14-fill-geo-gtshort` | 61.07 | yes | 37.0 / 24.0 | 23.8 / 24.2 | 0.020 / 0.038 |
| 3 | `Bayes-M14-fill-geo` | 61.08 | yes | 37.0 / 24.1 | 24.1 / 24.5 | 0.021 / 0.039 |
| 4 | `Bayes-M14-fill-geo-tvrw1` | 63.22 | yes | 39.3 / 23.9 | 23.2 / 23.5 | 0.020 / 0.036 |
| 5 | `Bayes-M14-fill-geo-gtlong` | 59.82 | yes | 35.8 / 24.0 | 24.3 / 24.5 | 0.021 / 0.040 |
| 6 | `Bayes-M16-fill-geo` | 59.28 | yes | 36.9 / 22.3 | 29.0 / 27.1 | 0.022 / 0.041 |
| 7 | `Bayes-M13c-fill-geo` | 56.75 | yes | 35.1 / 21.7 | 24.9 / 25.6 | 0.022 / 0.043 |
| 8 | `Bayes-M14-fill-med` | 55.98 | yes | 33.5 / 22.5 | 29.2 / 29.5 | 0.021 / 0.040 |
| 9 | `Bayes-M10-fill-geo` | 57.90 | yes | 36.1 / 21.8 | 28.5 / 28.8 | 0.022 / 0.043 |
| 10 | `Bayes-M17-fill-geo` | 57.32 | yes | 35.6 / 21.7 | 26.9 / 27.2 | 0.023 / 0.045 |
| 11 | `Bayes-M13c-fill-med` | 53.42 | yes | 33.2 / 20.2 | 30.5 / 30.7 | 0.023 / 0.044 |
| 12 | `Bayes-M17-fill-med` | 59.58 | yes | 39.0 / 20.6 | 33.4 / 32.9 | 0.024 / 0.047 |
| 13 | `Bayes-M16-fill-med` | 57.28 | yes | 36.1 / 21.1 | 36.4 / 32.9 | 0.022 / 0.042 |
| 14 | `Bayes-M13-fill-geo` | 55.54 | yes | 37.0 / 18.6 | 31.5 / 32.4 | 0.025 / 0.049 |
| 15 | `Bayes-M10-fill-med` | 55.37 | yes | 34.4 / 20.9 | 36.0 / 35.1 | 0.022 / 0.044 |
| 16 | `Bayes-M13-fill-med` | 54.54 | yes | 36.5 / 18.1 | 37.9 / 37.3 | 0.026 / 0.052 |
| 17 | `Bayes-M4c-geo` | 46.79 | yes | 29.9 / 16.8 | 32.1 / 34.8 | 0.024 / 0.047 |
| 18 | `Bayes-M4c-med` | 42.20 | yes | 26.9 / 15.3 | 34.8 / 36.9 | 0.025 / 0.049 |
| 19 | `Bayes-M8-fill-geo` | 46.45 | yes | 31.2 / 15.3 | 40.4 / 41.8 | 0.027 / 0.055 |
| 20 | `Bayes-M4-geo` | 36.73 | yes | 24.3 / 12.5 | 40.0 / 41.9 | 0.027 / 0.057 |
| 21 | `Bayes-M8-fill-med` | 41.23 | **spiky** | 29.3 / 11.9 | 53.0 / 50.8 | 0.029 / 0.061 |
| 22 | `Bayes-M4-med` | 29.72 | **spiky** | 19.9 / 9.8 | 51.2 / 50.5 | 0.029 / 0.062 |

_(no Best-renewal models in the evaluation)_
