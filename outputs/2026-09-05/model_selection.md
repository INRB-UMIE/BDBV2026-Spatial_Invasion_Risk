# Model selection — spatiotemporal

_Generated 2026-09-05T16:40:48+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M10-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_05092026 (processed_at 2026-09-05T12:00:00) |
| Predictions computed on snapshot | LINELIST_05092026 (processed_at 2026-09-05T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M10-med` | 66.82 | yes | 38.9 / 27.9 | 14.5 / 16.7 | 0.025 / 0.042 |
| 2 | `Bayes-M17-med` | 66.42 | yes | 39.4 / 27.0 | 15.1 / 17.2 | 0.026 / 0.045 |
| 3 | `Bayes-M10-geo` | 65.79 | yes | 38.9 / 26.9 | 13.6 / 15.6 | 0.024 / 0.041 |
| 4 | `Bayes-M17-geo` | 63.02 | yes | 37.4 / 25.6 | 15.5 / 17.2 | 0.025 / 0.043 |
| 5 | `Bayes-M17-dist-med` | 62.94 | yes | 36.8 / 26.1 | 13.5 / 15.4 | 0.025 / 0.043 |
| 6 | `Bayes-M8-dist` | 62.28 | yes | 36.7 / 25.6 | 15.1 / 16.9 | 0.026 / 0.046 |
| 7 | `Bayes-M8-med` | 61.77 | yes | 36.7 / 25.1 | 19.9 / 21.4 | 0.028 / 0.050 |
| 8 | `Bayes-M10-dist` | 59.07 | yes | 33.4 / 25.6 | 13.1 / 15.4 | 0.025 / 0.043 |
| 9 | `Bayes-M8-geo` | 57.92 | yes | 34.8 / 23.1 | 18.7 / 21.2 | 0.027 / 0.048 |
| 10 | `Bayes-M8-dist-geo` | 57.23 | yes | 33.7 / 23.5 | 15.6 / 17.2 | 0.026 / 0.044 |
| 11 | `Bayes-M17-dist-geo` | 55.92 | yes | 32.1 / 23.9 | 13.6 / 15.3 | 0.025 / 0.042 |
| 12 | `Bayes-M10-dist-geo` | 54.72 | yes | 31.8 / 22.9 | 12.7 / 14.7 | 0.025 / 0.042 |
| 13 | `Bayes-M4-dist-geo` | 37.57 | yes | 21.5 / 16.1 | 23.1 / 25.6 | 0.030 / 0.056 |
| 14 | `Bayes-M4-dist` | 35.87 | yes | 19.6 / 16.3 | 22.6 / 24.5 | 0.030 / 0.059 |
| 15 | `Bayes-M4-geo` | 37.97 | **spiky** | 22.3 / 15.7 | 32.1 / 35.3 | 0.031 / 0.060 |
| 16 | `Bayes-M4-med` | 34.15 | **spiky** | 19.1 / 15.1 | 32.9 / 34.7 | 0.033 / 0.064 |

_(no Best-renewal models in the evaluation)_
