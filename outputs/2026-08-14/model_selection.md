# Model selection — spatiotemporal

_Generated 2026-08-14T08:58:05+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_14082026 (processed_at 2026-08-14T12:00:00) |
| Predictions computed on snapshot | LINELIST_14082026 (processed_at 2026-08-14T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-med` | 64.16 | yes | 37.4 / 26.8 | 14.9 / 15.3 | 0.030 / 0.055 |
| 2 | `Bayes-M10-med` | 62.41 | yes | 35.3 / 27.1 | 14.8 / 15.1 | 0.030 / 0.052 |
| 3 | `Bayes-M17-dist-med` | 60.08 | yes | 34.3 / 25.8 | 14.3 / 14.7 | 0.030 / 0.054 |
| 4 | `Bayes-M8-dist` | 59.40 | yes | 33.9 / 25.5 | 15.0 / 15.4 | 0.031 / 0.056 |
| 5 | `Bayes-M8-med` | 58.65 | yes | 34.2 / 24.4 | 17.2 / 18.0 | 0.032 / 0.064 |
| 6 | `Bayes-M10-dist` | 58.14 | yes | 32.1 / 26.0 | 14.4 / 14.3 | 0.031 / 0.054 |
| 7 | `Bayes-M10-geo` | 52.22 | yes | 29.4 / 22.9 | 15.6 / 15.2 | 0.031 / 0.056 |
| 8 | `Bayes-M17-geo` | 51.36 | yes | 29.1 / 22.2 | 15.9 / 16.2 | 0.032 / 0.058 |
| 9 | `Bayes-M8-geo` | 47.69 | yes | 27.8 / 19.9 | 18.2 / 18.8 | 0.034 / 0.067 |
| 10 | `Bayes-M8-dist-geo` | 47.41 | yes | 26.3 / 21.1 | 16.4 / 16.7 | 0.033 / 0.060 |
| 11 | `Bayes-M17-dist-geo` | 46.25 | yes | 24.8 / 21.5 | 15.6 / 15.8 | 0.032 / 0.058 |
| 12 | `Bayes-M10-dist-geo` | 45.84 | yes | 24.7 / 21.1 | 15.4 / 14.9 | 0.032 / 0.057 |
| 13 | `Bayes-M4-dist` | 41.02 | yes | 23.5 / 17.6 | 21.4 / 23.0 | 0.038 / 0.085 |
| 14 | `Bayes-M4-med` | 41.55 | **spiky** | 24.7 / 16.8 | 29.8 / 31.7 | 0.041 / 0.094 |
| 15 | `Bayes-M4-dist-geo` | 29.92 | **spiky** | 16.6 / 13.4 | 22.9 / 25.0 | 0.039 / 0.085 |
| 16 | `Bayes-M4-geo` | 29.82 | **spiky** | 17.1 / 12.7 | 30.3 / 33.0 | 0.041 / 0.091 |

_(no Best-renewal models in the evaluation)_
