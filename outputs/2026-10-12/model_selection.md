# Model selection — spatiotemporal

_Generated 2026-10-15T00:17:07+0000_

**Featured Bayesian model (lowest non-spiky CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M14-fill-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_12102026 (processed_at 2026-10-12T12:00:00) |
| Predictions computed on snapshot | LINELIST_12102026 (processed_at 2026-10-12T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 10, margin 3) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 14, margin 4) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 70.92 | yes | 43.1 / 27.9 | 23.4 / 31.0 | 0.019 / 0.036 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 71.07 | yes | 43.6 / 27.5 | 23.4 / 30.9 | 0.019 / 0.035 |
| 3 | `Bayes-M14-fill-geo` | 68.31 | yes | 39.8 / 28.5 | 24.2 / 32.0 | 0.020 / 0.037 |
| 4 | `Bayes-M14-fill-geo-gtlong` | 67.57 | yes | 38.9 / 28.7 | 24.2 / 31.9 | 0.020 / 0.037 |
| 5 | `Bayes-M14-fill-geo-gtshort` | 67.99 | yes | 40.3 / 27.7 | 24.3 / 31.8 | 0.020 / 0.037 |
| 6 | `Bayes-M16-fill-geo` | 60.98 | yes | 35.3 / 25.7 | 24.0 / 32.5 | 0.021 / 0.038 |
| 7 | `Bayes-M14-fill-med` | 61.93 | yes | 36.9 / 25.0 | 28.5 / 38.0 | 0.021 / 0.040 |
| 8 | `Bayes-M10-fill-geo` | 61.63 | yes | 36.6 / 25.0 | 28.7 / 37.1 | 0.021 / 0.041 |
| 9 | `Bayes-M17-fill-geo` | 57.53 | yes | 33.4 / 24.2 | 27.0 / 34.7 | 0.022 / 0.042 |
| 10 | `Bayes-M13c-fill-geo` | 55.67 | yes | 32.3 / 23.4 | 25.8 / 32.3 | 0.022 / 0.041 |
| 11 | `Bayes-M16-fill-med` | 58.68 | yes | 34.8 / 23.9 | 29.0 / 39.3 | 0.021 / 0.039 |
| 12 | `Bayes-M13c-fill-med` | 52.77 | yes | 31.1 / 21.6 | 30.0 / 37.8 | 0.022 / 0.043 |
| 13 | `Bayes-M17-fill-med` | 54.31 | yes | 32.8 / 21.6 | 32.4 / 41.2 | 0.023 / 0.045 |
| 14 | `Bayes-M13-fill-geo` | 52.05 | yes | 30.7 / 21.3 | 31.5 / 39.4 | 0.024 / 0.045 |
| 15 | `Bayes-M10-fill-med` | 53.93 | yes | 32.5 / 21.4 | 35.4 / 45.8 | 0.022 / 0.045 |
| 16 | `Bayes-M13-fill-med` | 49.75 | yes | 30.7 / 19.1 | 35.7 / 44.4 | 0.025 / 0.049 |
| 17 | `Bayes-M4c-geo` | 41.47 | yes | 25.6 / 15.9 | 35.1 / 41.5 | 0.026 / 0.053 |
| 18 | `Bayes-M8-fill-geo` | 40.30 | yes | 24.2 / 16.1 | 41.1 / 48.7 | 0.026 / 0.054 |
| 19 | `Bayes-M4c-med` | 38.93 | yes | 24.2 / 14.8 | 36.6 / 43.3 | 0.026 / 0.054 |
| 20 | `Bayes-M4-geo` | 32.91 | yes | 19.5 / 13.4 | 41.8 / 49.1 | 0.027 / 0.057 |
| 21 | `Bayes-M8-fill-med` | 31.78 | **spiky** | 19.4 / 12.4 | 50.7 / 60.3 | 0.029 / 0.067 |
| 22 | `Bayes-M4-med` | 26.46 | **spiky** | 16.5 / 9.9 | 51.9 / 60.8 | 0.031 / 0.075 |

_(no Best-renewal models in the evaluation)_
