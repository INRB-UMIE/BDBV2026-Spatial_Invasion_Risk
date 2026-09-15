# Model selection — spatiotemporal

_Generated 2026-09-15T16:10:05+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M10-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_15092026 (processed_at 2026-09-15T12:00:00) |
| Predictions computed on snapshot | LINELIST_15092026 (processed_at 2026-09-15T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M10-med` | 76.27 | yes | 45.6 / 30.7 | 14.4 / 13.7 | 0.022 / 0.039 |
| 2 | `Bayes-M17-med` | 75.42 | yes | 45.6 / 29.8 | 15.2 / 14.7 | 0.024 / 0.042 |
| 3 | `Bayes-M10-geo` | 74.23 | yes | 44.5 / 29.7 | 13.0 / 13.1 | 0.022 / 0.039 |
| 4 | `Bayes-M17-geo` | 73.13 | yes | 44.6 / 28.5 | 14.9 / 14.7 | 0.023 / 0.041 |
| 5 | `Bayes-M8-dist` | 69.57 | yes | 42.3 / 27.2 | 14.9 / 14.9 | 0.024 / 0.044 |
| 6 | `Bayes-M8-geo` | 69.20 | yes | 42.7 / 26.5 | 19.0 / 18.4 | 0.025 / 0.046 |
| 7 | `Bayes-M8-med` | 68.77 | yes | 41.3 / 27.5 | 19.8 / 19.6 | 0.026 / 0.048 |
| 8 | `Bayes-M17-dist-med` | 68.40 | yes | 40.4 / 28.0 | 13.2 / 13.2 | 0.023 / 0.041 |
| 9 | `Bayes-M8-dist-geo` | 64.97 | yes | 38.8 / 26.2 | 14.8 / 14.8 | 0.024 / 0.043 |
| 10 | `Bayes-M10-dist` | 62.01 | yes | 35.0 / 27.0 | 13.3 / 12.8 | 0.023 / 0.041 |
| 11 | `Bayes-M17-dist-geo` | 61.42 | yes | 35.9 / 25.5 | 12.6 / 13.2 | 0.023 / 0.041 |
| 12 | `Bayes-M10-dist-geo` | 59.12 | yes | 34.1 / 25.0 | 12.0 / 12.3 | 0.023 / 0.040 |
| 13 | `Bayes-M4-dist-geo` | 39.17 | yes | 22.0 / 17.1 | 21.7 / 22.2 | 0.030 / 0.056 |
| 14 | `Bayes-M4-geo` | 45.24 | **spiky** | 27.2 / 18.0 | 29.6 / 29.5 | 0.031 / 0.058 |
| 15 | `Bayes-M4-med` | 37.52 | **spiky** | 20.9 / 16.6 | 32.7 / 32.8 | 0.032 / 0.063 |
| 16 | `Bayes-M4-dist` | 35.71 | **spiky** | 19.0 / 16.7 | 23.0 / 23.1 | 0.030 / 0.057 |

_(no Best-renewal models in the evaluation)_
