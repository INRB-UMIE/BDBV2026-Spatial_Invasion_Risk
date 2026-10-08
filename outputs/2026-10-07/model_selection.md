# Model selection — spatiotemporal

_Generated 2026-10-08T12:52:16+0000_

**Featured Bayesian model (lowest non-spiky CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M14-fill-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_07102026 (processed_at 2026-10-07T12:00:00) |
| Predictions computed on snapshot | LINELIST_07102026 (processed_at 2026-10-07T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 10, margin 5) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 16, margin 1) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 49.50 | yes | 31.9 / 17.5 | 22.1 / 24.0 | 0.020 / 0.040 |
| 2 | `Bayes-M14-fill-geo-gtlong` | 49.59 | yes | 31.7 / 17.9 | 22.6 / 24.7 | 0.020 / 0.040 |
| 3 | `Bayes-M14-fill-geo` | 49.44 | yes | 31.9 / 17.5 | 22.4 / 24.4 | 0.020 / 0.041 |
| 4 | `Bayes-M14-fill-geo-tvrw1` | 49.65 | yes | 32.5 / 17.2 | 22.0 / 23.9 | 0.020 / 0.040 |
| 5 | `Bayes-M14-fill-geo-gtshort` | 49.10 | yes | 32.2 / 16.9 | 22.6 / 24.4 | 0.020 / 0.041 |
| 6 | `Bayes-M14-fill-med` | 48.31 | yes | 30.6 / 17.8 | 26.8 / 29.5 | 0.021 / 0.041 |
| 7 | `Bayes-M16-fill-geo` | 44.99 | yes | 28.9 / 16.1 | 22.9 / 25.8 | 0.021 / 0.043 |
| 8 | `Bayes-M13c-fill-geo` | 42.70 | yes | 27.8 / 14.9 | 22.8 / 25.1 | 0.022 / 0.045 |
| 9 | `Bayes-M10-fill-geo` | 46.17 | yes | 31.1 / 15.0 | 25.7 / 28.1 | 0.021 / 0.046 |
| 10 | `Bayes-M16-fill-med` | 44.69 | yes | 28.2 / 16.5 | 27.6 / 31.1 | 0.021 / 0.043 |
| 11 | `Bayes-M13c-fill-med` | 42.80 | yes | 27.0 / 15.8 | 27.0 / 29.6 | 0.022 / 0.044 |
| 12 | `Bayes-M17-fill-geo` | 40.52 | yes | 26.7 / 13.8 | 25.0 / 27.0 | 0.022 / 0.047 |
| 13 | `Bayes-M10-fill-med` | 44.76 | yes | 29.8 / 14.9 | 31.3 / 34.9 | 0.022 / 0.046 |
| 14 | `Bayes-M17-fill-med` | 42.25 | yes | 27.7 / 14.5 | 29.3 / 32.0 | 0.023 / 0.047 |
| 15 | `Bayes-M13-fill-geo` | 36.67 | yes | 24.8 / 11.9 | 28.8 / 31.2 | 0.024 / 0.051 |
| 16 | `Bayes-M13-fill-med` | 38.44 | yes | 26.4 / 12.1 | 33.1 / 36.4 | 0.024 / 0.052 |
| 17 | `Bayes-M4c-geo` | 31.31 | yes | 20.6 / 10.7 | 31.1 / 34.2 | 0.025 / 0.054 |
| 18 | `Bayes-M4c-med` | 30.20 | yes | 19.5 / 10.7 | 33.1 / 36.4 | 0.025 / 0.055 |
| 19 | `Bayes-M8-fill-geo` | 30.70 | yes | 21.5 / 9.2 | 37.7 / 40.6 | 0.026 / 0.061 |
| 20 | `Bayes-M4-geo` | 20.63 | yes | 13.3 / 7.4 | 38.4 / 41.1 | 0.028 / 0.067 |
| 21 | `Bayes-M8-fill-med` | 29.18 | **spiky** | 20.2 / 9.0 | 45.2 / 49.9 | 0.028 / 0.071 |
| 22 | `Bayes-M4-med` | 19.41 | **spiky** | 13.0 / 6.4 | 44.2 / 49.9 | 0.029 / 0.079 |

_(no Best-renewal models in the evaluation)_
