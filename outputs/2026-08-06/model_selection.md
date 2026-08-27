# Model selection — spatiotemporal

_Generated 2026-08-07T01:32:54+0000_

**Featured Bayesian model (CV composite):** `Bayes-M17-dist-short`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-dist-short`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_06082026 (processed_at 2026-08-06T12:00:00) |
| Predictions computed on snapshot | LINELIST_06082026 (processed_at 2026-08-06T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Within each horizon, rank by AUC-PR skill (desc), mean rank-of-truth (asc) and log-score (asc); sum the three ranks; pool by summing across horizons; winner = lowest pooled rank-sum among methods covering the most horizons, ties broken by higher total AUC-PR skill. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-dist-short` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-dist-short` |

### Bayesian model leaderboard (by pooled CV composite; lower = better)

| Rank | Model | Composite (pooled) | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|
| 1 | `Bayes-M17-dist-short` | 18 | 29.4 / 20.8 | 13.6 / 14.6 | 0.032 / 0.057 |
| 2 | `Bayes-M10-med` | 20 | 29.3 / 21.6 | 14.2 / 14.4 | 0.032 / 0.055 |
| 3 | `Bayes-M17-dist-long` | 21 | 29.8 / 21.8 | 13.7 / 14.7 | 0.034 / 0.058 |
| 4 | `Bayes-M17-dist-med` | 25 | 29.3 / 21.1 | 13.7 / 14.7 | 0.033 / 0.057 |
| 5 | `Bayes-M17-short` | 27 | 29.6 / 20.6 | 14.1 / 15.1 | 0.032 / 0.057 |
| 6 | `Bayes-M17-long` | 29 | 30.3 / 21.3 | 14.1 / 14.8 | 0.034 / 0.059 |
| 7 | `Bayes-M9-med` | 39 | 29.7 / 20.9 | 15.9 / 16.0 | 0.034 / 0.059 |
| 8 | `Bayes-M17-med` | 44 | 29.1 / 20.5 | 14.2 / 15.3 | 0.033 / 0.059 |
| 9 | `Bayes-M8-short` | 62 | 29.2 / 19.5 | 15.9 / 16.5 | 0.035 / 0.064 |
| 10 | `Bayes-M8-long` | 70 | 29.5 / 20.8 | 16.3 / 16.8 | 0.036 / 0.066 |
| 11 | `Bayes-M8-med` | 74 | 28.6 / 19.6 | 16.3 / 16.6 | 0.036 / 0.065 |
| 12 | `Bayes-M17-dist-short-full` | 82 | 21.5 / 17.1 | 16.4 / 17.4 | 0.035 / 0.062 |
| 13 | `Bayes-M10-full` | 83 | 21.9 / 17.6 | 18.2 / 17.3 | 0.035 / 0.060 |
| 14 | `Bayes-M17-dist-full` | 87 | 21.1 / 17.9 | 16.6 / 17.4 | 0.036 / 0.062 |
| 15 | `Bayes-M17-dist-long-full` | 92 | 21.1 / 18.1 | 16.9 / 17.6 | 0.036 / 0.063 |
| 16 | `Bayes-M17-full` | 98 | 21.2 / 17.0 | 17.0 / 18.3 | 0.036 / 0.063 |
| 17 | `Bayes-M17-short-full` | 98 | 21.1 / 16.6 | 17.2 / 18.3 | 0.035 / 0.063 |
| 18 | `Bayes-M17-long-full` | 106 | 20.9 / 17.7 | 17.0 / 18.0 | 0.037 / 0.064 |
| 19 | `Bayes-M15-short` | 117 | 25.4 / 17.9 | 35.8 / 33.3 | 0.037 / 0.065 |
| 20 | `Bayes-M9-full` | 117 | 20.9 / 16.5 | 17.8 / 18.9 | 0.037 / 0.065 |
| 21 | `Bayes-M15-med` | 127 | 25.1 / 18.1 | 35.8 / 37.3 | 0.038 / 0.067 |
| 22 | `Bayes-M15-long` | 129 | 24.7 / 18.2 | 35.8 / 37.3 | 0.039 / 0.068 |
| 23 | `Bayes-M8-full` | 129 | 21.2 / 16.0 | 21.1 / 20.5 | 0.038 / 0.070 |
| 24 | `Bayes-M15-short-full` | 147 | 19.6 / 15.3 | 40.7 / 36.2 | 0.038 / 0.066 |
| 25 | `Bayes-M8-susp` | 155 | 10.3 / 9.0 | 19.9 / 22.7 | 0.039 / 0.072 |
| 26 | `Bayes-M15-full` | 156 | 18.5 / 15.2 | 40.8 / 36.5 | 0.039 / 0.068 |
| 27 | `Bayes-M8-full-susp` | 156 | 11.7 / 9.7 | 22.4 / 24.2 | 0.040 / 0.072 |
| 28 | `Bayes-M4-med` | 159 | 14.7 / 12.2 | 27.0 / 29.6 | 0.043 / 0.092 |
| 29 | `Bayes-M15-long-full` | 160 | 18.3 / 15.1 | 39.6 / 36.3 | 0.040 / 0.069 |
| 30 | `Bayes-M4-full` | 163 | 14.0 / 12.0 | 28.8 / 32.3 | 0.043 / 0.090 |

_(no Best-renewal models in the evaluation)_
