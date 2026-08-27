# Model selection — spatiotemporal

_Generated 2026-08-13T10:52:48+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_13082026 (processed_at 2026-08-13T12:00:00) |
| Predictions computed on snapshot | LINELIST_13082026 (processed_at 2026-08-13T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-med` | 57.61 | yes | 34.6 / 23.0 | 14.6 / 15.3 | 0.031 / 0.056 |
| 2 | `Bayes-M10-med` | 57.40 | yes | 34.0 / 23.4 | 13.6 / 14.6 | 0.030 / 0.053 |
| 3 | `Bayes-M8-dist` | 55.23 | yes | 32.5 / 22.7 | 15.5 / 15.3 | 0.032 / 0.058 |
| 4 | `Bayes-M8-med` | 54.76 | yes | 32.7 / 22.0 | 17.4 / 17.8 | 0.033 / 0.063 |
| 5 | `Bayes-M10-dist` | 54.21 | yes | 30.4 / 23.8 | 13.2 / 14.1 | 0.031 / 0.054 |
| 6 | `Bayes-M17-dist-med` | 53.46 | yes | 30.9 / 22.5 | 14.0 / 14.7 | 0.032 / 0.055 |
| 7 | `Bayes-M10-geo` | 41.44 | yes | 23.0 / 18.5 | 14.9 / 15.1 | 0.033 / 0.058 |
| 8 | `Bayes-M17-geo` | 40.92 | yes | 22.4 / 18.5 | 16.0 / 16.5 | 0.034 / 0.061 |
| 9 | `Bayes-M17-dist-geo` | 38.75 | yes | 20.4 / 18.3 | 15.9 / 15.7 | 0.034 / 0.060 |
| 10 | `Bayes-M8-dist-geo` | 38.53 | yes | 21.0 / 17.5 | 17.7 / 16.9 | 0.035 / 0.063 |
| 11 | `Bayes-M8-geo` | 38.10 | yes | 21.3 / 16.8 | 18.6 / 18.9 | 0.035 / 0.068 |
| 12 | `Bayes-M10-dist-geo` | 37.40 | yes | 19.8 / 17.6 | 14.7 / 15.1 | 0.034 / 0.060 |
| 13 | `Bayes-M4-dist` | 29.33 | yes | 15.8 / 13.5 | 20.9 / 22.2 | 0.039 / 0.085 |
| 14 | `Bayes-M4-dist-geo` | 24.40 | yes | 12.8 / 11.6 | 21.7 / 23.9 | 0.041 / 0.086 |
| 15 | `Bayes-M4-med` | 28.54 | **spiky** | 15.7 / 12.8 | 28.5 / 30.4 | 0.042 / 0.093 |
| 16 | `Bayes-M4-geo` | 24.31 | **spiky** | 12.7 / 11.6 | 28.5 / 31.5 | 0.043 / 0.091 |

_(no Best-renewal models in the evaluation)_
