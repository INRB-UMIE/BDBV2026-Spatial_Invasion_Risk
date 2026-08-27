# Model selection — spatiotemporal

_Generated 2026-08-26T06:41:04+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_23082026 (processed_at 2026-08-23T12:00:00) |
| Predictions computed on snapshot | LINELIST_23082026 (processed_at 2026-08-23T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-med` | 79.22 | yes | 48.4 / 30.8 | 18.4 / 18.2 | 0.026 / 0.048 |
| 2 | `Bayes-M10-med` | 78.12 | yes | 47.4 / 30.7 | 18.3 / 18.4 | 0.025 / 0.045 |
| 3 | `Bayes-M8-dist` | 76.85 | yes | 47.2 / 29.6 | 18.5 / 19.0 | 0.027 / 0.049 |
| 4 | `Bayes-M8-med` | 75.73 | yes | 46.4 / 29.3 | 22.4 / 22.3 | 0.028 / 0.054 |
| 5 | `Bayes-M17-dist-med` | 75.47 | yes | 46.1 / 29.4 | 17.0 / 17.0 | 0.026 / 0.046 |
| 6 | `Bayes-M17-geo` | 74.49 | yes | 46.3 / 28.2 | 18.9 / 18.6 | 0.026 / 0.046 |
| 7 | `Bayes-M10-geo` | 74.13 | yes | 46.2 / 28.0 | 17.0 / 17.6 | 0.025 / 0.044 |
| 8 | `Bayes-M10-dist` | 73.31 | yes | 44.5 / 28.8 | 16.7 / 17.1 | 0.026 / 0.046 |
| 9 | `Bayes-M8-geo` | 70.66 | yes | 44.6 / 26.1 | 22.6 / 22.3 | 0.028 / 0.052 |
| 10 | `Bayes-M8-dist-geo` | 69.94 | yes | 43.4 / 26.5 | 19.0 / 19.5 | 0.027 / 0.049 |
| 11 | `Bayes-M17-dist-geo` | 68.06 | yes | 41.3 / 26.8 | 16.5 / 17.4 | 0.026 / 0.046 |
| 12 | `Bayes-M10-dist-geo` | 62.60 | yes | 37.0 / 25.6 | 15.6 / 16.6 | 0.026 / 0.046 |
| 13 | `Bayes-M4-dist-geo` | 49.02 | yes | 29.3 / 19.7 | 26.4 / 27.5 | 0.032 / 0.061 |
| 14 | `Bayes-M4-dist` | 45.75 | yes | 26.1 / 19.6 | 25.9 / 26.1 | 0.032 / 0.065 |
| 15 | `Bayes-M4-geo` | 51.07 | **spiky** | 30.7 / 20.4 | 35.2 / 37.0 | 0.033 / 0.065 |
| 16 | `Bayes-M4-med` | 48.49 | **spiky** | 28.8 / 19.7 | 35.9 / 35.9 | 0.034 / 0.071 |

_(no Best-renewal models in the evaluation)_
