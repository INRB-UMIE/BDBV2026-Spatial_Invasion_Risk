# Model selection — spatiotemporal

_Generated 2026-08-21T08:34:33+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_21082026 (processed_at 2026-08-21T12:00:00) |
| Predictions computed on snapshot | LINELIST_21082026 (processed_at 2026-08-21T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-med` | 64.56 | yes | 37.6 / 26.9 | 16.1 / 15.5 | 0.030 / 0.052 |
| 2 | `Bayes-M10-med` | 64.47 | yes | 37.0 / 27.4 | 15.8 / 15.4 | 0.029 / 0.049 |
| 3 | `Bayes-M17-dist-med` | 62.19 | yes | 36.3 / 25.8 | 14.6 / 14.7 | 0.030 / 0.052 |
| 4 | `Bayes-M8-dist` | 62.11 | yes | 35.8 / 26.3 | 15.1 / 15.7 | 0.030 / 0.053 |
| 5 | `Bayes-M8-med` | 61.05 | yes | 35.8 / 25.2 | 17.8 / 19.0 | 0.032 / 0.060 |
| 6 | `Bayes-M10-dist` | 59.27 | yes | 33.3 / 25.9 | 15.0 / 14.5 | 0.030 / 0.051 |
| 7 | `Bayes-M17-geo` | 55.19 | yes | 31.8 / 23.4 | 15.3 / 16.1 | 0.031 / 0.054 |
| 8 | `Bayes-M10-geo` | 54.30 | yes | 31.4 / 22.9 | 15.3 / 15.3 | 0.030 / 0.052 |
| 9 | `Bayes-M8-geo` | 50.96 | yes | 30.3 / 20.7 | 17.7 / 19.4 | 0.033 / 0.062 |
| 10 | `Bayes-M8-dist-geo` | 50.17 | yes | 28.3 / 21.9 | 16.0 / 16.5 | 0.032 / 0.057 |
| 11 | `Bayes-M17-dist-geo` | 47.76 | yes | 26.0 / 21.8 | 15.1 / 15.5 | 0.031 / 0.054 |
| 12 | `Bayes-M10-dist-geo` | 46.98 | yes | 26.6 / 20.4 | 14.5 / 14.8 | 0.031 / 0.054 |
| 13 | `Bayes-M4-dist` | 40.71 | yes | 23.3 / 17.4 | 22.0 / 24.0 | 0.038 / 0.083 |
| 14 | `Bayes-M4-med` | 42.43 | **spiky** | 25.0 / 17.4 | 31.1 / 33.5 | 0.042 / 0.091 |
| 15 | `Bayes-M4-geo` | 32.29 | **spiky** | 18.7 / 13.6 | 27.0 / 34.0 | 0.041 / 0.088 |
| 16 | `Bayes-M4-dist-geo` | 31.64 | **spiky** | 18.0 / 13.7 | 21.0 / 25.0 | 0.040 / 0.084 |

_(no Best-renewal models in the evaluation)_
