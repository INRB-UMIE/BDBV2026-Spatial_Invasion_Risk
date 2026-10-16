# Model selection — spatiotemporal

_Generated 2026-10-16T12:48:41+0000_

**Featured Bayesian model (lowest non-spiky CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M14-fill-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_15102026 (processed_at 2026-10-15T12:00:00) |
| Predictions computed on snapshot | LINELIST_15102026 (processed_at 2026-10-15T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 10, margin 3) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 12, margin 3) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 48.62 | yes | 29.1 / 19.5 | 32.4 / 32.0 | 0.025 / 0.049 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 48.57 | yes | 29.6 / 19.0 | 32.4 / 31.9 | 0.025 / 0.047 |
| 3 | `Bayes-M14-fill-geo-gtlong` | 45.83 | yes | 26.2 / 19.6 | 33.6 / 35.1 | 0.027 / 0.050 |
| 4 | `Bayes-M14-fill-geo` | 44.77 | yes | 25.7 / 19.1 | 33.6 / 34.7 | 0.027 / 0.051 |
| 5 | `Bayes-M14-fill-geo-gtshort` | 44.67 | yes | 26.1 / 18.6 | 33.8 / 34.4 | 0.026 / 0.051 |
| 6 | `Bayes-M13c-fill-geo` | 42.53 | yes | 25.4 / 17.1 | 30.4 / 32.7 | 0.029 / 0.055 |
| 7 | `Bayes-M16-fill-geo` | 43.90 | yes | 26.6 / 17.3 | 35.4 / 36.4 | 0.028 / 0.053 |
| 8 | `Bayes-M10-fill-geo` | 42.97 | yes | 25.6 / 17.4 | 37.4 / 38.9 | 0.028 / 0.056 |
| 9 | `Bayes-M17-fill-geo` | 40.80 | yes | 25.3 / 15.5 | 33.1 / 35.6 | 0.030 / 0.060 |
| 10 | `Bayes-M14-fill-med` | 41.27 | yes | 23.6 / 17.7 | 38.8 / 40.9 | 0.027 / 0.053 |
| 11 | `Bayes-M13c-fill-med` | 38.69 | yes | 22.7 / 16.0 | 35.9 / 38.0 | 0.029 / 0.057 |
| 12 | `Bayes-M13-fill-geo` | 39.76 | yes | 25.3 / 14.5 | 37.2 / 39.9 | 0.031 / 0.064 |
| 13 | `Bayes-M16-fill-med` | 40.99 | yes | 24.4 / 16.6 | 40.1 / 42.0 | 0.028 / 0.054 |
| 14 | `Bayes-M10-fill-med` | 38.59 | yes | 23.4 / 15.1 | 44.9 / 45.9 | 0.028 / 0.059 |
| 15 | `Bayes-M4c-geo` | 31.30 | yes | 18.4 / 12.9 | 37.3 / 41.0 | 0.031 / 0.063 |
| 16 | `Bayes-M17-fill-med` | 38.23 | yes | 23.1 / 15.1 | 39.9 / 42.2 | 0.030 / 0.061 |
| 17 | `Bayes-M4c-med` | 29.89 | yes | 17.6 / 12.2 | 39.3 / 43.0 | 0.031 / 0.066 |
| 18 | `Bayes-M13-fill-med` | 35.32 | yes | 22.6 / 12.7 | 43.7 / 45.9 | 0.032 / 0.069 |
| 19 | `Bayes-M8-fill-geo` | 30.56 | yes | 19.6 / 11.0 | 45.0 / 48.6 | 0.034 / 0.076 |
| 20 | `Bayes-M4-geo` | 20.91 | yes | 11.8 / 9.1 | 44.8 / 48.5 | 0.035 / 0.080 |
| 21 | `Bayes-M8-fill-med` | 23.86 | yes | 16.1 / 7.8 | 55.3 / 59.3 | 0.037 / 0.092 |
| 22 | `Bayes-M4-med` | 17.50 | yes | 10.6 / 6.9 | 54.3 / 58.2 | 0.037 / 0.097 |

_(no Best-renewal models in the evaluation)_
