# Model selection — spatiotemporal

_Generated 2026-09-04T18:29:07+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_04092026 (processed_at 2026-09-04T12:00:00) |
| Predictions computed on snapshot | LINELIST_04092026 (processed_at 2026-09-04T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-med` | 65.03 | yes | 38.6 / 26.4 | 17.8 / 17.5 | 0.028 / 0.052 |
| 2 | `Bayes-M10-med` | 63.17 | yes | 37.3 / 25.9 | 16.6 / 16.9 | 0.027 / 0.049 |
| 3 | `Bayes-M8-dist` | 62.66 | yes | 37.0 / 25.7 | 16.3 / 17.0 | 0.029 / 0.053 |
| 4 | `Bayes-M17-dist-med` | 61.67 | yes | 36.7 / 25.0 | 15.9 / 16.1 | 0.028 / 0.051 |
| 5 | `Bayes-M8-med` | 61.43 | yes | 37.4 / 24.1 | 20.4 / 21.1 | 0.030 / 0.060 |
| 6 | `Bayes-M10-dist` | 56.97 | yes | 32.4 / 24.5 | 15.6 / 15.8 | 0.028 / 0.050 |
| 7 | `Bayes-M17-geo` | 53.93 | yes | 31.6 / 22.4 | 16.3 / 17.0 | 0.029 / 0.052 |
| 8 | `Bayes-M10-geo` | 53.89 | yes | 31.5 / 22.4 | 14.9 / 15.9 | 0.028 / 0.051 |
| 9 | `Bayes-M8-dist-geo` | 49.01 | yes | 28.3 / 20.7 | 16.0 / 17.2 | 0.030 / 0.055 |
| 10 | `Bayes-M8-geo` | 48.30 | yes | 29.2 / 19.1 | 18.6 / 20.9 | 0.031 / 0.061 |
| 11 | `Bayes-M17-dist-geo` | 47.46 | yes | 26.6 / 20.9 | 15.5 / 15.8 | 0.029 / 0.052 |
| 12 | `Bayes-M10-dist-geo` | 44.24 | yes | 24.4 / 19.8 | 14.2 / 15.1 | 0.029 / 0.053 |
| 13 | `Bayes-M4-dist` | 36.89 | yes | 21.1 / 15.8 | 20.8 / 23.5 | 0.034 / 0.073 |
| 14 | `Bayes-M4-dist-geo` | 31.93 | yes | 18.3 / 13.7 | 19.2 / 24.2 | 0.034 / 0.070 |
| 15 | `Bayes-M4-med` | 37.94 | **spiky** | 22.9 / 15.0 | 29.5 / 32.8 | 0.036 / 0.079 |
| 16 | `Bayes-M4-geo` | 31.33 | **spiky** | 18.0 / 13.3 | 25.7 / 33.0 | 0.036 / 0.073 |

_(no Best-renewal models in the evaluation)_
