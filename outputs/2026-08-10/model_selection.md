# Model selection — spatiotemporal

_Generated 2026-08-10T08:24:22+0000_

**Featured Bayesian model (CV composite):** `Bayes-M10-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_10082026 (processed_at 2026-08-10T12:00:00) |
| Predictions computed on snapshot | LINELIST_10082026 (processed_at 2026-08-10T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Within each horizon, rank by AUC-PR skill (desc), mean rank-of-truth (asc) and log-score (asc); sum the three ranks; pool by summing across horizons; winner = lowest pooled rank-sum among methods covering the most horizons, ties broken by higher total AUC-PR skill. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-med` |

### Bayesian model leaderboard (by pooled CV composite; lower = better)

| Rank | Model | Composite (pooled) | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|
| 1 | `Bayes-M10-med` | 16 | 42.3 / 30.1 | 14.2 / 16.0 | 0.027 / 0.048 |
| 2 | `Bayes-M10-full` | 31 | 40.6 / 27.9 | 13.7 / 15.8 | 0.026 / 0.046 |
| 3 | `Bayes-M17-dist-short` | 38 | 39.8 / 29.7 | 14.2 / 16.0 | 0.027 / 0.049 |
| 4 | `Bayes-M17-med` | 41 | 43.8 / 30.5 | 15.1 / 17.4 | 0.028 / 0.050 |
| 5 | `Bayes-M17-short` | 41 | 42.3 / 30.0 | 15.1 / 17.3 | 0.027 / 0.049 |
| 6 | `Bayes-M17-dist-med` | 47 | 40.6 / 29.7 | 14.3 / 16.2 | 0.028 / 0.050 |
| 7 | `Bayes-M17-short-full` | 50 | 41.6 / 27.9 | 15.3 / 17.4 | 0.027 / 0.048 |
| 8 | `Bayes-M17-dist-long` | 52 | 41.1 / 30.1 | 14.3 / 16.2 | 0.029 / 0.052 |
| 9 | `Bayes-M17-full` | 59 | 41.2 / 28.5 | 15.4 / 17.5 | 0.028 / 0.049 |
| 10 | `Bayes-M17-long` | 60 | 43.0 / 30.5 | 15.2 / 17.4 | 0.029 / 0.053 |
| 11 | `Bayes-M17-dist-short-full` | 64 | 39.2 / 26.8 | 14.9 / 16.7 | 0.027 / 0.049 |
| 12 | `Bayes-M17-dist-full` | 65 | 39.3 / 27.2 | 14.8 / 16.8 | 0.028 / 0.049 |
| 13 | `Bayes-M17-dist-long-full` | 79 | 38.6 / 27.6 | 15.0 / 17.0 | 0.029 / 0.051 |
| 14 | `Bayes-M17-long-full` | 80 | 39.9 / 28.4 | 15.7 / 17.8 | 0.029 / 0.051 |
| 15 | `Bayes-M9-med` | 82 | 42.0 / 29.8 | 17.0 / 19.0 | 0.029 / 0.053 |
| 16 | `Bayes-M8-long` | 105 | 42.1 / 29.4 | 17.4 / 18.9 | 0.031 / 0.059 |
| 17 | `Bayes-M8-med` | 105 | 40.7 / 28.9 | 17.3 / 18.8 | 0.030 / 0.057 |
| 18 | `Bayes-M8-short` | 107 | 40.6 / 28.7 | 17.7 / 19.1 | 0.029 / 0.056 |
| 19 | `Bayes-M9-full` | 109 | 39.3 / 27.2 | 17.3 / 19.6 | 0.029 / 0.052 |
| 20 | `Bayes-M15-short-full` | 110 | 36.9 / 22.5 | 15.9 / 22.9 | 0.028 / 0.052 |
| 21 | `Bayes-M8-full` | 115 | 39.4 / 27.0 | 17.1 / 19.5 | 0.029 / 0.054 |
| 22 | `Bayes-M15-short` | 120 | 40.3 / 25.3 | 24.8 / 27.5 | 0.029 / 0.053 |
| 23 | `Bayes-M15-full` | 122 | 35.3 / 22.2 | 15.7 / 22.8 | 0.029 / 0.052 |
| 24 | `Bayes-M15-med` | 134 | 39.8 / 25.7 | 24.9 / 27.6 | 0.030 / 0.054 |
| 25 | `Bayes-M15-long-full` | 135 | 32.0 / 21.9 | 16.0 / 23.0 | 0.030 / 0.054 |
| 26 | `Bayes-M15-long` | 146 | 38.8 / 26.1 | 25.4 / 27.7 | 0.030 / 0.056 |
| 27 | `Bayes-M8-full-susp` | 163 | 16.5 / 12.4 | 21.8 / 30.6 | 0.031 / 0.058 |
| 28 | `Bayes-M8-susp` | 170 | 19.9 / 11.9 | 23.7 / 39.4 | 0.032 / 0.064 |
| 29 | `Bayes-M4-full` | 171 | 28.9 / 21.3 | 32.2 / 34.6 | 0.036 / 0.069 |
| 30 | `Bayes-M4-med` | 173 | 27.1 / 19.7 | 31.2 / 32.2 | 0.038 / 0.079 |

_(no Best-renewal models in the evaluation)_
