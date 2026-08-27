# Model selection — spatiotemporal

_Generated 2026-08-18T12:07:01+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M10-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_17082026 (processed_at 2026-08-17T12:00:00) |
| Predictions computed on snapshot | LINELIST_17082026 (processed_at 2026-08-17T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M10-med` | 71.71 | yes | 42.1 / 29.7 | 15.5 / 15.6 | 0.027 / 0.046 |
| 2 | `Bayes-M17-med` | 70.41 | yes | 41.3 / 29.1 | 16.0 / 16.2 | 0.028 / 0.049 |
| 3 | `Bayes-M10-geo` | 68.30 | yes | 40.6 / 27.7 | 14.7 / 14.9 | 0.026 / 0.045 |
| 4 | `Bayes-M8-med` | 67.90 | yes | 40.6 / 27.3 | 18.7 / 19.5 | 0.030 / 0.056 |
| 5 | `Bayes-M8-dist` | 67.73 | yes | 39.8 / 27.9 | 16.2 / 16.8 | 0.029 / 0.051 |
| 6 | `Bayes-M17-geo` | 67.04 | yes | 39.7 / 27.3 | 15.8 / 15.8 | 0.027 / 0.047 |
| 7 | `Bayes-M17-dist-med` | 66.11 | yes | 38.3 / 27.8 | 15.0 / 15.3 | 0.028 / 0.048 |
| 8 | `Bayes-M10-dist` | 65.36 | yes | 38.1 / 27.3 | 14.4 / 14.7 | 0.028 / 0.048 |
| 9 | `Bayes-M17-dist-geo` | 63.94 | yes | 37.7 / 26.2 | 15.2 / 15.3 | 0.027 / 0.047 |
| 10 | `Bayes-M8-geo` | 63.28 | yes | 38.3 / 25.0 | 18.6 / 19.3 | 0.029 / 0.053 |
| 11 | `Bayes-M8-dist-geo` | 63.05 | yes | 37.4 / 25.7 | 16.4 / 16.8 | 0.028 / 0.050 |
| 12 | `Bayes-M10-dist-geo` | 62.60 | yes | 36.9 / 25.7 | 14.3 / 14.4 | 0.028 / 0.047 |
| 13 | `Bayes-M4-dist` | 41.88 | yes | 24.1 / 17.8 | 23.2 / 24.1 | 0.036 / 0.070 |
| 14 | `Bayes-M4-geo` | 48.44 | **spiky** | 28.6 / 19.9 | 33.3 / 34.5 | 0.036 / 0.068 |
| 15 | `Bayes-M4-dist-geo` | 45.51 | **spiky** | 26.6 / 18.9 | 24.7 / 25.4 | 0.034 / 0.064 |
| 16 | `Bayes-M4-med` | 44.73 | **spiky** | 26.6 / 18.1 | 32.1 / 33.3 | 0.038 / 0.077 |

_(no Best-renewal models in the evaluation)_
