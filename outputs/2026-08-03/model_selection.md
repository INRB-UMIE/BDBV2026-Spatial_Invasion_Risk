# Model selection — spatiotemporal

_Generated 2026-08-04T18:45:05+0000_

**Featured Bayesian model (CV composite):** `Bayes-M10-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_03082026 (processed_at 2026-08-03T12:00:00) |
| Predictions computed on snapshot | LINELIST_03082026 (processed_at 2026-08-03T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Within each horizon, rank by AUC-PR skill (desc), mean rank-of-truth (asc) and log-score (asc); sum the three ranks; pool by summing across horizons; winner = lowest pooled rank-sum among methods covering the most horizons, ties broken by higher total AUC-PR skill. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-med` |

### Bayesian model leaderboard (by pooled CV composite; lower = better)

| Rank | Model | Composite (pooled) | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|
| 1 | `Bayes-M10-med` | 8 | 50.5 / 31.9 | 12.9 / 14.2 | 0.024 / 0.046 |
| 2 | `Bayes-M10-full` | 15 | 49.8 / 30.5 | 15.5 / 15.9 | 0.024 / 0.045 |
| 3 | `Bayes-M8-long` | 29 | 50.5 / 30.8 | 16.5 / 17.3 | 0.028 / 0.056 |
| 4 | `Bayes-M9-med` | 29 | 48.9 / 31.4 | 17.3 / 17.4 | 0.027 / 0.051 |
| 5 | `Bayes-M8-med` | 31 | 49.3 / 30.6 | 16.9 / 17.3 | 0.027 / 0.055 |
| 6 | `Bayes-M8-short` | 34 | 48.8 / 30.1 | 17.1 / 17.2 | 0.027 / 0.053 |
| 7 | `Bayes-M9-full` | 34 | 48.8 / 30.0 | 19.9 / 19.8 | 0.026 / 0.049 |
| 8 | `Bayes-M8-full` | 36 | 50.1 / 29.8 | 20.4 / 20.5 | 0.026 / 0.052 |
| 9 | `Bayes-M4-full` | 62 | 34.7 / 22.0 | 37.9 / 35.1 | 0.034 / 0.071 |
| 10 | `Bayes-M8-full-susp` | 62 | 23.7 / 14.0 | 31.1 / 41.7 | 0.030 / 0.059 |
| 11 | `Bayes-M4-med` | 64 | 30.3 / 19.6 | 36.5 / 33.9 | 0.036 / 0.080 |
| 12 | `Bayes-M8-susp` | 64 | 28.7 / 14.9 | 36.0 / 48.6 | 0.032 / 0.065 |

_(no Best-renewal models in the evaluation)_
