# Model selection — spatiotemporal

_Generated 2026-10-19T15:19:06+0000_

**Featured Bayesian model (lowest non-spiky CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M14-fill-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_16102026 (processed_at 2026-10-16T12:00:00) |
| Predictions computed on snapshot | LINELIST_16102026 (processed_at 2026-10-16T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 9, margin 0) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 10, margin 3) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 55.32 | yes | 32.6 / 22.7 | 18.4 / 24.9 | 0.023 / 0.039 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 56.30 | yes | 33.5 / 22.8 | 18.5 / 25.1 | 0.023 / 0.038 |
| 3 | `Bayes-M14-fill-geo-gtshort` | 50.37 | yes | 28.6 / 21.8 | 18.6 / 26.2 | 0.024 / 0.041 |
| 4 | `Bayes-M14-fill-geo-gtlong` | 49.39 | yes | 28.1 / 21.3 | 18.9 / 26.1 | 0.024 / 0.042 |
| 5 | `Bayes-M14-fill-geo` | 49.29 | yes | 27.9 / 21.4 | 18.6 / 25.9 | 0.024 / 0.042 |
| 6 | `Bayes-M10-fill-geo` | 48.46 | yes | 28.3 / 20.1 | 21.0 / 28.5 | 0.025 / 0.044 |
| 7 | `Bayes-M13c-fill-geo` | 47.58 | yes | 27.7 / 19.9 | 20.5 / 26.5 | 0.026 / 0.046 |
| 8 | `Bayes-M16-fill-geo` | 48.48 | yes | 28.1 / 20.3 | 25.1 / 25.8 | 0.025 / 0.044 |
| 9 | `Bayes-M17-fill-geo` | 47.73 | yes | 28.8 / 18.9 | 21.7 / 28.5 | 0.027 / 0.048 |
| 10 | `Bayes-M14-fill-med` | 47.55 | yes | 26.9 / 20.7 | 20.8 / 31.3 | 0.024 / 0.043 |
| 11 | `Bayes-M13-fill-geo` | 46.58 | yes | 28.7 / 17.8 | 26.0 / 32.7 | 0.028 / 0.051 |
| 12 | `Bayes-M17-fill-med` | 47.57 | yes | 28.8 / 18.8 | 24.1 / 33.0 | 0.027 / 0.050 |
| 13 | `Bayes-M10-fill-med` | 46.59 | yes | 27.0 / 19.6 | 23.5 / 34.7 | 0.025 / 0.047 |
| 14 | `Bayes-M13c-fill-med` | 44.94 | yes | 26.5 / 18.5 | 22.8 / 30.5 | 0.026 / 0.047 |
| 15 | `Bayes-M16-fill-med` | 47.21 | yes | 27.0 / 20.2 | 29.4 / 30.8 | 0.025 / 0.045 |
| 16 | `Bayes-M4c-geo` | 37.06 | yes | 21.6 / 15.5 | 28.0 / 34.3 | 0.027 / 0.050 |
| 17 | `Bayes-M13-fill-med` | 46.26 | yes | 29.1 / 17.2 | 27.7 / 37.0 | 0.029 / 0.055 |
| 18 | `Bayes-M4c-med` | 35.47 | yes | 21.1 / 14.4 | 28.3 / 36.0 | 0.027 / 0.052 |
| 19 | `Bayes-M8-fill-geo` | 37.99 | yes | 23.6 / 14.4 | 34.1 / 40.7 | 0.030 / 0.057 |
| 20 | `Bayes-M4-geo` | 31.50 | yes | 20.0 / 11.5 | 34.3 / 41.1 | 0.030 / 0.059 |
| 21 | `Bayes-M8-fill-med` | 36.31 | yes | 24.9 / 11.4 | 36.4 / 45.7 | 0.032 / 0.064 |
| 22 | `Bayes-M4-med` | 26.76 | yes | 17.5 / 9.3 | 36.3 / 46.3 | 0.032 / 0.066 |

_(no Best-renewal models in the evaluation)_
