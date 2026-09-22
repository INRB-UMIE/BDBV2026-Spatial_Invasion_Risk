# Model selection — spatiotemporal

_Generated 2026-09-22T18:31:05+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M10-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_20092026 (processed_at 2026-09-20T12:00:00) |
| Predictions computed on snapshot | LINELIST_20092026 (processed_at 2026-09-20T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M10-med` | 75.17 | yes | 43.5 / 31.7 | 14.8 / 14.5 | 0.023 / 0.039 |
| 2 | `Bayes-M17-med` | 74.31 | yes | 43.9 / 30.4 | 15.5 / 15.7 | 0.024 / 0.043 |
| 3 | `Bayes-M17-dist-med` | 71.30 | yes | 42.3 / 29.0 | 14.2 / 14.4 | 0.024 / 0.042 |
| 4 | `Bayes-M10-geo` | 70.96 | yes | 42.3 / 28.6 | 14.3 / 14.2 | 0.022 / 0.039 |
| 5 | `Bayes-M10-dist` | 70.78 | yes | 41.2 / 29.6 | 13.7 / 13.6 | 0.024 / 0.041 |
| 6 | `Bayes-M8-dist` | 69.90 | yes | 40.7 / 29.2 | 15.4 / 15.7 | 0.025 / 0.044 |
| 7 | `Bayes-M8-med` | 69.66 | yes | 40.8 / 28.8 | 19.9 / 20.0 | 0.027 / 0.049 |
| 8 | `Bayes-M17-geo` | 67.32 | yes | 40.0 / 27.3 | 16.5 / 16.5 | 0.023 / 0.041 |
| 9 | `Bayes-M8-geo` | 63.69 | yes | 39.0 / 24.7 | 20.5 / 20.3 | 0.025 / 0.047 |
| 10 | `Bayes-M8-dist-geo` | 63.03 | yes | 37.4 / 25.6 | 16.4 / 16.5 | 0.024 / 0.043 |
| 11 | `Bayes-M17-dist-geo` | 61.51 | yes | 35.5 / 26.0 | 14.7 / 14.9 | 0.023 / 0.041 |
| 12 | `Bayes-M10-dist-geo` | 58.33 | yes | 32.8 / 25.6 | 13.1 / 13.3 | 0.023 / 0.041 |
| 13 | `Bayes-M4-dist` | 42.13 | yes | 23.2 / 18.9 | 21.8 / 22.0 | 0.029 / 0.057 |
| 14 | `Bayes-M4-dist-geo` | 42.04 | yes | 23.8 / 18.2 | 23.8 / 23.8 | 0.029 / 0.054 |
| 15 | `Bayes-M4-geo` | 43.37 | **spiky** | 24.9 / 18.5 | 31.8 / 31.5 | 0.030 / 0.057 |
| 16 | `Bayes-M4-med` | 41.28 | **spiky** | 22.9 / 18.4 | 31.0 / 31.0 | 0.032 / 0.062 |

_(no Best-renewal models in the evaluation)_
