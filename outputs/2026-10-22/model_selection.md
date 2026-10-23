# Model selection — spatiotemporal

_Generated 2026-10-22T23:33:55+0000_

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
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 9, margin 4) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 14, margin 2) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 50.17 | yes | 30.4 / 19.7 | 28.0 / 33.3 | 0.024 / 0.048 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 49.92 | yes | 30.8 / 19.1 | 28.0 / 32.7 | 0.024 / 0.047 |
| 3 | `Bayes-M14-fill-geo-gtlong` | 46.85 | yes | 26.7 / 20.1 | 29.1 / 36.6 | 0.026 / 0.050 |
| 4 | `Bayes-M14-fill-geo-gtshort` | 46.52 | yes | 27.2 / 19.3 | 30.0 / 36.5 | 0.026 / 0.050 |
| 5 | `Bayes-M14-fill-geo` | 45.62 | yes | 26.7 / 18.9 | 29.4 / 36.8 | 0.026 / 0.050 |
| 6 | `Bayes-M13c-fill-geo` | 42.88 | yes | 25.8 / 17.0 | 28.3 / 32.9 | 0.028 / 0.056 |
| 7 | `Bayes-M10-fill-geo` | 43.19 | yes | 26.6 / 16.6 | 31.3 / 39.1 | 0.027 / 0.056 |
| 8 | `Bayes-M16-fill-geo` | 43.07 | yes | 25.8 / 17.3 | 34.3 / 38.0 | 0.027 / 0.052 |
| 9 | `Bayes-M17-fill-geo` | 41.20 | yes | 25.7 / 15.5 | 30.5 / 35.9 | 0.029 / 0.059 |
| 10 | `Bayes-M14-fill-med` | 42.96 | yes | 24.9 / 18.1 | 35.4 / 44.8 | 0.026 / 0.052 |
| 11 | `Bayes-M13c-fill-med` | 39.05 | yes | 23.0 / 16.0 | 34.6 / 40.4 | 0.029 / 0.057 |
| 12 | `Bayes-M10-fill-med` | 40.66 | yes | 25.3 / 15.4 | 38.3 / 48.5 | 0.028 / 0.060 |
| 13 | `Bayes-M16-fill-med` | 40.88 | yes | 24.3 / 16.6 | 42.1 / 45.6 | 0.028 / 0.054 |
| 14 | `Bayes-M13-fill-geo` | 38.32 | yes | 24.7 / 13.6 | 35.7 / 40.6 | 0.030 / 0.065 |
| 15 | `Bayes-M17-fill-med` | 39.40 | yes | 24.2 / 15.2 | 39.2 / 46.2 | 0.029 / 0.061 |
| 16 | `Bayes-M4c-geo` | 30.71 | yes | 18.3 / 12.4 | 32.8 / 40.1 | 0.030 / 0.063 |
| 17 | `Bayes-M4c-med` | 28.96 | yes | 17.5 / 11.5 | 35.8 / 43.8 | 0.030 / 0.067 |
| 18 | `Bayes-M13-fill-med` | 36.03 | yes | 23.5 / 12.6 | 43.8 / 50.0 | 0.031 / 0.068 |
| 19 | `Bayes-M8-fill-geo` | 30.31 | yes | 19.6 / 10.7 | 40.2 / 47.3 | 0.033 / 0.074 |
| 20 | `Bayes-M4-geo` | 20.73 | yes | 11.8 / 8.9 | 39.8 / 46.9 | 0.034 / 0.079 |
| 21 | `Bayes-M8-fill-med` | 25.07 | yes | 16.8 / 8.3 | 45.1 / 56.0 | 0.035 / 0.092 |
| 22 | `Bayes-M4-med` | 17.70 | yes | 10.8 / 6.9 | 44.7 / 55.4 | 0.036 / 0.099 |

_(no Best-renewal models in the evaluation)_
