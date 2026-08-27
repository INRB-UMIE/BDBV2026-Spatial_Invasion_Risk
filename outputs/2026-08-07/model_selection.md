# Model selection — spatiotemporal

_Generated 2026-08-07T11:58:22+0000_

**Featured Bayesian model (CV composite):** `Bayes-M10-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_07082026 (processed_at 2026-08-07T12:00:00) |
| Predictions computed on snapshot | LINELIST_07082026 (processed_at 2026-08-07T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Within each horizon, rank by AUC-PR skill (desc), mean rank-of-truth (asc) and log-score (asc); sum the three ranks; pool by summing across horizons; winner = lowest pooled rank-sum among methods covering the most horizons, ties broken by higher total AUC-PR skill. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-med` |

### Bayesian model leaderboard (by pooled CV composite; lower = better)

| Rank | Model | Composite (pooled) | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|
| 1 | `Bayes-M10-med` | 17 | 36.7 / 27.7 | 14.7 / 15.8 | 0.031 / 0.053 |
| 2 | `Bayes-M17-dist-short` | 26 | 36.7 / 27.0 | 15.2 / 16.1 | 0.031 / 0.055 |
| 3 | `Bayes-M17-short` | 31 | 37.5 / 27.1 | 15.6 / 17.1 | 0.031 / 0.056 |
| 4 | `Bayes-M17-med` | 33 | 38.1 / 28.3 | 15.7 / 17.3 | 0.032 / 0.057 |
| 5 | `Bayes-M17-dist-long` | 34 | 37.0 / 27.9 | 15.1 / 16.2 | 0.033 / 0.058 |
| 6 | `Bayes-M17-dist-med` | 34 | 36.5 / 26.9 | 15.1 / 16.2 | 0.032 / 0.056 |
| 7 | `Bayes-M17-long` | 48 | 36.9 / 28.3 | 15.7 / 17.3 | 0.033 / 0.059 |
| 8 | `Bayes-M10-full` | 55 | 29.6 / 22.8 | 15.2 / 15.5 | 0.033 / 0.056 |
| 9 | `Bayes-M9-med` | 63 | 36.6 / 26.4 | 15.8 / 18.0 | 0.033 / 0.059 |
| 10 | `Bayes-M17-dist-full` | 77 | 29.9 / 23.0 | 15.9 / 16.5 | 0.034 / 0.059 |
| 11 | `Bayes-M17-dist-short-full` | 81 | 29.6 / 22.4 | 16.1 / 16.5 | 0.033 / 0.058 |
| 12 | `Bayes-M17-short-full` | 83 | 30.6 / 22.4 | 16.6 / 17.4 | 0.033 / 0.058 |
| 13 | `Bayes-M17-full` | 84 | 31.3 / 23.1 | 16.3 / 17.5 | 0.033 / 0.059 |
| 14 | `Bayes-M17-dist-long-full` | 84 | 29.8 / 23.7 | 16.0 / 16.5 | 0.035 / 0.060 |
| 15 | `Bayes-M8-short` | 89 | 36.3 / 25.5 | 16.5 / 18.4 | 0.033 / 0.064 |
| 16 | `Bayes-M17-long-full` | 91 | 29.9 / 23.1 | 16.1 / 17.2 | 0.035 / 0.060 |
| 17 | `Bayes-M8-med` | 92 | 36.2 / 25.9 | 16.6 / 18.1 | 0.034 / 0.065 |
| 18 | `Bayes-M8-long` | 103 | 35.5 / 26.5 | 16.8 / 18.3 | 0.035 / 0.067 |
| 19 | `Bayes-M9-full` | 103 | 30.1 / 22.2 | 16.5 / 18.0 | 0.034 / 0.061 |
| 20 | `Bayes-M15-short` | 111 | 38.4 / 25.2 | 36.5 / 29.3 | 0.038 / 0.065 |
| 21 | `Bayes-M15-med` | 115 | 38.5 / 25.4 | 36.7 / 29.5 | 0.038 / 0.066 |
| 22 | `Bayes-M8-full` | 129 | 29.2 / 20.9 | 16.8 / 18.4 | 0.035 / 0.067 |
| 23 | `Bayes-M15-long` | 134 | 36.7 / 25.1 | 36.9 / 37.6 | 0.039 / 0.070 |
| 24 | `Bayes-M15-full` | 137 | 28.2 / 21.4 | 24.2 / 25.5 | 0.039 / 0.064 |
| 25 | `Bayes-M15-short-full` | 138 | 28.2 / 20.7 | 24.6 / 25.5 | 0.038 / 0.063 |
| 26 | `Bayes-M15-long-full` | 152 | 26.7 / 21.1 | 24.5 / 26.7 | 0.040 / 0.068 |
| 27 | `Bayes-M8-susp` | 152 | 14.8 / 12.9 | 18.4 / 20.8 | 0.038 / 0.071 |
| 28 | `Bayes-M8-full-susp` | 156 | 13.4 / 11.9 | 18.5 / 20.5 | 0.039 / 0.071 |
| 29 | `Bayes-M4-med` | 168 | 22.9 / 16.9 | 26.4 / 30.3 | 0.043 / 0.092 |
| 30 | `Bayes-M4-full` | 170 | 19.2 / 15.1 | 27.7 / 31.6 | 0.042 / 0.085 |

_(no Best-renewal models in the evaluation)_
