# Model selection — spatiotemporal

_Generated 2026-10-17T00:14:26+0000_

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
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 9, margin 2) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 9, margin 8) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 72.08 | yes | 45.4 / 26.6 | 32.1 / 31.8 | 0.019 / 0.036 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 71.17 | yes | 45.4 / 25.7 | 31.6 / 31.8 | 0.019 / 0.035 |
| 3 | `Bayes-M14-fill-geo-gtshort` | 65.69 | yes | 40.1 / 25.6 | 32.8 / 33.1 | 0.020 / 0.039 |
| 4 | `Bayes-M14-fill-geo` | 65.68 | yes | 40.2 / 25.5 | 33.3 / 33.2 | 0.021 / 0.040 |
| 5 | `Bayes-M14-fill-geo-gtlong` | 64.78 | yes | 39.5 / 25.3 | 33.5 / 33.4 | 0.021 / 0.040 |
| 6 | `Bayes-M13c-fill-geo` | 63.47 | yes | 39.4 / 24.0 | 31.0 / 32.4 | 0.022 / 0.043 |
| 7 | `Bayes-M17-fill-geo` | 63.58 | yes | 40.6 / 23.0 | 33.3 / 34.7 | 0.023 / 0.045 |
| 8 | `Bayes-M16-fill-geo` | 61.56 | yes | 38.2 / 23.4 | 38.5 / 32.1 | 0.022 / 0.041 |
| 9 | `Bayes-M10-fill-geo` | 61.07 | yes | 37.4 / 23.7 | 37.5 / 38.0 | 0.022 / 0.043 |
| 10 | `Bayes-M14-fill-med` | 60.63 | yes | 36.2 / 24.4 | 39.2 / 40.5 | 0.021 / 0.041 |
| 11 | `Bayes-M13c-fill-med` | 58.25 | yes | 36.2 / 22.0 | 37.5 / 38.5 | 0.023 / 0.045 |
| 12 | `Bayes-M13-fill-geo` | 60.21 | yes | 39.9 / 20.4 | 38.3 / 39.7 | 0.025 / 0.049 |
| 13 | `Bayes-M17-fill-med` | 61.36 | yes | 39.4 / 21.9 | 41.6 / 42.0 | 0.024 / 0.048 |
| 14 | `Bayes-M16-fill-med` | 58.80 | yes | 35.8 / 23.0 | 45.8 / 38.5 | 0.022 / 0.043 |
| 15 | `Bayes-M10-fill-med` | 58.20 | yes | 35.8 / 22.4 | 45.5 / 46.6 | 0.022 / 0.045 |
| 16 | `Bayes-M13-fill-med` | 56.56 | yes | 37.4 / 19.1 | 45.7 / 46.3 | 0.026 / 0.054 |
| 17 | `Bayes-M4c-geo` | 45.81 | yes | 28.4 / 17.4 | 39.7 / 42.3 | 0.024 / 0.048 |
| 18 | `Bayes-M4c-med` | 43.41 | yes | 27.4 / 16.0 | 42.6 / 45.1 | 0.025 / 0.050 |
| 19 | `Bayes-M8-fill-geo` | 49.11 | yes | 32.7 / 16.4 | 46.5 / 49.3 | 0.027 / 0.055 |
| 20 | `Bayes-M4-geo` | 38.77 | yes | 26.2 / 12.6 | 47.7 / 50.1 | 0.027 / 0.057 |
| 21 | `Bayes-M8-fill-med` | 43.46 | **spiky** | 30.7 / 12.7 | 60.4 / 60.1 | 0.029 / 0.063 |
| 22 | `Bayes-M4-med` | 32.18 | **spiky** | 21.6 / 10.6 | 60.8 / 60.9 | 0.029 / 0.064 |

_(no Best-renewal models in the evaluation)_
