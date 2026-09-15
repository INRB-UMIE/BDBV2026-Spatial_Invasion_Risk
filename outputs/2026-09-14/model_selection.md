# Model selection — spatiotemporal

_Generated 2026-09-15T12:15:50+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_14092026 (processed_at 2026-09-14T12:00:00) |
| Predictions computed on snapshot | LINELIST_14092026 (processed_at 2026-09-14T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-med` | 68.56 | yes | 40.1 / 28.4 | 18.9 / 18.5 | 0.024 / 0.043 |
| 2 | `Bayes-M10-med` | 67.57 | yes | 39.4 / 28.2 | 17.6 / 17.1 | 0.024 / 0.040 |
| 3 | `Bayes-M8-med` | 65.14 | yes | 38.7 / 26.4 | 22.8 / 23.0 | 0.027 / 0.048 |
| 4 | `Bayes-M10-geo` | 65.11 | yes | 37.4 / 27.7 | 16.4 / 15.5 | 0.023 / 0.038 |
| 5 | `Bayes-M17-geo` | 64.44 | yes | 37.0 / 27.4 | 18.6 / 18.0 | 0.024 / 0.040 |
| 6 | `Bayes-M8-dist` | 63.49 | yes | 36.7 / 26.7 | 18.1 / 18.0 | 0.025 / 0.044 |
| 7 | `Bayes-M8-dist-geo` | 61.24 | yes | 35.1 / 26.1 | 18.3 / 17.8 | 0.025 / 0.042 |
| 8 | `Bayes-M8-geo` | 60.99 | yes | 35.4 / 25.6 | 22.4 / 22.0 | 0.026 / 0.045 |
| 9 | `Bayes-M17-dist-geo` | 60.10 | yes | 34.6 / 25.5 | 16.3 / 15.8 | 0.024 / 0.040 |
| 10 | `Bayes-M17-dist-med` | 60.09 | yes | 34.6 / 25.5 | 16.6 / 16.5 | 0.024 / 0.041 |
| 11 | `Bayes-M10-dist` | 58.08 | yes | 33.1 / 25.0 | 16.4 / 16.0 | 0.024 / 0.041 |
| 12 | `Bayes-M10-dist-geo` | 55.48 | yes | 30.3 / 25.2 | 15.0 / 14.3 | 0.024 / 0.040 |
| 13 | `Bayes-M4-dist-geo` | 39.45 | yes | 21.9 / 17.6 | 24.3 / 24.7 | 0.029 / 0.051 |
| 14 | `Bayes-M4-dist` | 36.91 | yes | 20.5 / 16.4 | 23.4 / 23.7 | 0.029 / 0.054 |
| 15 | `Bayes-M4-geo` | 41.57 | **spiky** | 23.6 / 18.0 | 33.4 / 33.4 | 0.030 / 0.054 |
| 16 | `Bayes-M4-med` | 35.97 | **spiky** | 19.5 / 16.5 | 33.0 / 33.5 | 0.031 / 0.060 |

_(no Best-renewal models in the evaluation)_
