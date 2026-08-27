# Model selection — spatiotemporal

_Generated 2026-08-26T11:54:15+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_26082026 (processed_at 2026-08-26T12:00:00) |
| Predictions computed on snapshot | LINELIST_26082026 (processed_at 2026-08-26T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-geo` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-geo` | 63.80 | yes | 38.0 / 25.8 | 46.0 / 28.0 | 0.025 / 0.045 |
| 2 | `Bayes-M17-med` | 62.93 | yes | 37.3 / 25.6 | 50.0 / 29.7 | 0.025 / 0.046 |
| 3 | `Bayes-M10-med` | 61.42 | yes | 36.4 / 25.0 | 48.9 / 29.7 | 0.024 / 0.045 |
| 4 | `Bayes-M8-med` | 60.51 | yes | 36.4 / 24.2 | 53.9 / 34.0 | 0.027 / 0.052 |
| 5 | `Bayes-M8-geo` | 60.07 | yes | 35.6 / 24.5 | 49.7 / 31.5 | 0.026 / 0.051 |
| 6 | `Bayes-M8-dist` | 58.53 | yes | 34.6 / 23.9 | 52.3 / 31.5 | 0.026 / 0.048 |
| 7 | `Bayes-M10-geo` | 58.29 | yes | 33.6 / 24.6 | 44.5 / 27.0 | 0.024 / 0.044 |
| 8 | `Bayes-M17-dist-med` | 58.14 | yes | 34.0 / 24.1 | 50.3 / 29.5 | 0.025 / 0.046 |
| 9 | `Bayes-M8-dist-geo` | 58.05 | yes | 34.5 / 23.5 | 49.1 / 29.7 | 0.025 / 0.047 |
| 10 | `Bayes-M17-dist-geo` | 55.12 | yes | 31.9 / 23.2 | 46.4 / 27.5 | 0.025 / 0.045 |
| 11 | `Bayes-M10-dist` | 51.33 | yes | 29.1 / 22.2 | 51.0 / 30.1 | 0.025 / 0.046 |
| 12 | `Bayes-M10-dist-geo` | 49.30 | yes | 27.9 / 21.4 | 45.9 / 27.4 | 0.025 / 0.046 |
| 13 | `Bayes-M4-geo` | 36.13 | yes | 20.2 / 16.0 | 58.7 / 46.4 | 0.032 / 0.065 |
| 14 | `Bayes-M4-dist-geo` | 33.12 | yes | 18.3 / 14.8 | 54.0 / 37.9 | 0.031 / 0.062 |
| 15 | `Bayes-M4-dist` | 32.01 | yes | 17.3 / 14.8 | 58.9 / 39.9 | 0.031 / 0.064 |
| 16 | `Bayes-M4-med` | 31.77 | yes | 17.3 / 14.4 | 66.1 / 49.1 | 0.033 / 0.069 |

_(no Best-renewal models in the evaluation)_
