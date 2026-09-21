# Model selection — spatiotemporal

_Generated 2026-09-21T07:50:46+0000_

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
| 1 | `Bayes-M10-med` | 77.16 | yes | 45.2 / 31.9 | 14.8 / 14.7 | 0.022 / 0.039 |
| 2 | `Bayes-M17-med` | 76.86 | yes | 46.0 / 30.8 | 15.6 / 16.0 | 0.023 / 0.042 |
| 3 | `Bayes-M17-dist-med` | 72.74 | yes | 43.4 / 29.4 | 14.1 / 14.4 | 0.023 / 0.041 |
| 4 | `Bayes-M10-geo` | 72.62 | yes | 43.8 / 28.8 | 14.4 / 14.3 | 0.022 / 0.039 |
| 5 | `Bayes-M8-dist` | 72.32 | yes | 42.7 / 29.7 | 15.5 / 15.9 | 0.024 / 0.043 |
| 6 | `Bayes-M8-med` | 72.13 | yes | 42.9 / 29.2 | 19.8 / 20.4 | 0.026 / 0.049 |
| 7 | `Bayes-M10-dist` | 71.27 | yes | 41.3 / 30.0 | 13.7 / 13.7 | 0.023 / 0.041 |
| 8 | `Bayes-M17-geo` | 69.50 | yes | 42.0 / 27.5 | 16.6 / 16.5 | 0.023 / 0.041 |
| 9 | `Bayes-M8-geo` | 65.57 | yes | 40.7 / 24.9 | 20.3 / 20.5 | 0.025 / 0.047 |
| 10 | `Bayes-M8-dist-geo` | 65.05 | yes | 39.1 / 26.0 | 16.4 / 16.5 | 0.023 / 0.043 |
| 11 | `Bayes-M17-dist-geo` | 63.18 | yes | 36.8 / 26.4 | 14.7 / 14.9 | 0.023 / 0.041 |
| 12 | `Bayes-M10-dist-geo` | 58.78 | yes | 33.3 / 25.5 | 13.2 / 13.3 | 0.023 / 0.041 |
| 13 | `Bayes-M4-dist` | 43.54 | yes | 24.3 / 19.2 | 21.8 / 22.1 | 0.028 / 0.057 |
| 14 | `Bayes-M4-dist-geo` | 43.53 | yes | 25.0 / 18.5 | 23.8 / 24.0 | 0.028 / 0.054 |
| 15 | `Bayes-M4-geo` | 45.37 | **spiky** | 26.6 / 18.7 | 31.6 / 31.7 | 0.029 / 0.057 |
| 16 | `Bayes-M4-med` | 42.74 | **spiky** | 24.1 / 18.6 | 30.9 / 31.3 | 0.030 / 0.062 |

_(no Best-renewal models in the evaluation)_
