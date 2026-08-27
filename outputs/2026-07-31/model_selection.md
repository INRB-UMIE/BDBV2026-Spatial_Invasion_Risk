# Model selection — spatiotemporal

_Generated 2026-07-31T09:40:43+0000_

**Featured Bayesian model (CV composite):** `Bayes-M10-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_31072026 (processed_at 2026-07-31T12:00:00) |
| Predictions computed on snapshot | LINELIST_31072026 (processed_at 2026-07-31T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Within each horizon, rank by AUC-PR skill (desc), mean rank-of-truth (asc) and log-score (asc); sum the three ranks; pool by summing across horizons; winner = lowest pooled rank-sum among methods covering the most horizons, ties broken by higher total AUC-PR skill. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-med` |

### Bayesian model leaderboard (by pooled CV composite; lower = better)

| Rank | Model | Composite (pooled) | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|
| 1 | `Bayes-M10-med` | 8 | 34.6 / 23.4 | 16.8 / 18.2 | 0.033 / 0.060 |
| 2 | `Bayes-M9-med` | 17 | 34.6 / 23.4 | 18.1 / 19.5 | 0.034 / 0.064 |
| 3 | `Bayes-M8-short` | 23 | 34.0 / 22.2 | 18.6 / 18.6 | 0.034 / 0.068 |
| 4 | `Bayes-M10-full` | 25 | 29.1 / 19.5 | 18.1 / 19.5 | 0.035 / 0.063 |
| 5 | `Bayes-M8-long` | 28 | 34.6 / 23.3 | 19.1 / 18.8 | 0.036 / 0.071 |
| 6 | `Bayes-M8-med` | 28 | 33.9 / 22.4 | 18.7 / 18.7 | 0.035 / 0.069 |
| 7 | `Bayes-M9-full` | 41 | 27.7 / 18.5 | 19.7 / 22.1 | 0.037 / 0.068 |
| 8 | `Bayes-M8-full` | 49 | 27.8 / 17.8 | 21.0 / 22.6 | 0.037 / 0.073 |
| 9 | `Bayes-M8-susp` | 55 | 15.0 / 12.1 | 20.9 / 21.1 | 0.040 / 0.077 |
| 10 | `Bayes-M4-med` | 64 | 21.3 / 14.6 | 28.7 / 31.4 | 0.043 / 0.098 |
| 11 | `Bayes-M8-full-susp` | 64 | 13.4 / 10.7 | 22.6 / 24.6 | 0.042 / 0.079 |
| 12 | `Bayes-M4-full` | 66 | 17.5 / 13.2 | 29.4 / 34.4 | 0.043 / 0.091 |

_(no Best-renewal models in the evaluation)_
