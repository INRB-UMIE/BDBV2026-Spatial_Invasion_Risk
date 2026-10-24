# Model selection — spatiotemporal

_Generated 2026-10-23T22:59:26+0000_

**Featured Bayesian model (lowest non-spiky CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M14-fill-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_22102026 (processed_at 2026-10-22T12:00:00) |
| Predictions computed on snapshot | LINELIST_22102026 (processed_at 2026-10-22T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 12, margin 3) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 17, margin 3) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 52.47 | yes | 31.3 / 21.1 | 26.4 / 33.7 | 0.024 / 0.046 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 51.76 | yes | 31.4 / 20.3 | 26.1 / 33.3 | 0.024 / 0.045 |
| 3 | `Bayes-M14-fill-geo-gtlong` | 50.37 | yes | 28.9 / 21.5 | 27.2 / 35.8 | 0.025 / 0.047 |
| 4 | `Bayes-M14-fill-geo-gtshort` | 51.30 | yes | 30.1 / 21.2 | 27.4 / 35.9 | 0.025 / 0.046 |
| 5 | `Bayes-M14-fill-geo` | 50.01 | yes | 29.1 / 20.9 | 27.3 / 35.7 | 0.025 / 0.047 |
| 6 | `Bayes-M13c-fill-geo` | 46.36 | yes | 27.8 / 18.6 | 25.8 / 31.8 | 0.027 / 0.052 |
| 7 | `Bayes-M16-fill-geo` | 45.96 | yes | 27.3 / 18.6 | 28.4 / 35.3 | 0.026 / 0.049 |
| 8 | `Bayes-M10-fill-geo` | 47.16 | yes | 28.2 / 19.0 | 30.0 / 37.9 | 0.026 / 0.052 |
| 9 | `Bayes-M17-fill-geo` | 44.82 | yes | 27.6 / 17.2 | 27.9 / 33.9 | 0.028 / 0.056 |
| 10 | `Bayes-M14-fill-med` | 45.53 | yes | 26.0 / 19.6 | 31.3 / 43.0 | 0.026 / 0.048 |
| 11 | `Bayes-M13c-fill-med` | 41.37 | yes | 24.1 / 17.3 | 29.9 / 37.3 | 0.028 / 0.053 |
| 12 | `Bayes-M16-fill-med` | 42.69 | yes | 24.9 / 17.8 | 32.5 / 41.5 | 0.027 / 0.050 |
| 13 | `Bayes-M13-fill-geo` | 41.05 | yes | 25.9 / 15.1 | 31.4 / 37.5 | 0.030 / 0.061 |
| 14 | `Bayes-M10-fill-med` | 42.54 | yes | 25.9 / 16.6 | 34.8 / 46.8 | 0.027 / 0.056 |
| 15 | `Bayes-M17-fill-med` | 41.65 | yes | 25.1 / 16.6 | 33.6 / 43.7 | 0.028 / 0.057 |
| 16 | `Bayes-M4c-geo` | 32.83 | yes | 19.4 / 13.4 | 33.6 / 39.1 | 0.029 / 0.059 |
| 17 | `Bayes-M13-fill-med` | 37.04 | yes | 23.5 / 13.5 | 36.9 / 47.8 | 0.031 / 0.064 |
| 18 | `Bayes-M4c-med` | 30.43 | yes | 18.2 / 12.3 | 35.2 / 42.0 | 0.030 / 0.062 |
| 19 | `Bayes-M8-fill-geo` | 32.57 | yes | 21.0 / 11.6 | 39.9 / 44.8 | 0.032 / 0.069 |
| 20 | `Bayes-M4-geo` | 21.80 | yes | 12.3 / 9.5 | 39.6 / 44.7 | 0.033 / 0.074 |
| 21 | `Bayes-M8-fill-med` | 24.74 | yes | 16.5 / 8.3 | 43.8 / 54.2 | 0.035 / 0.084 |
| 22 | `Bayes-M4-med` | 17.67 | yes | 10.7 / 7.0 | 44.1 / 54.7 | 0.036 / 0.090 |

_(no Best-renewal models in the evaluation)_
