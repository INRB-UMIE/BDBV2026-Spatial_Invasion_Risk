# Model selection — spatiotemporal

_Generated 2026-08-28T10:07:34+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M10-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_28082026 (processed_at 2026-08-28T12:00:00) |
| Predictions computed on snapshot | LINELIST_28082026 (processed_at 2026-08-28T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M10-med` | 66.60 | yes | 38.5 / 28.2 | 15.8 / 16.4 | 0.027 / 0.047 |
| 2 | `Bayes-M17-med` | 65.30 | yes | 37.7 / 27.6 | 15.8 / 16.0 | 0.027 / 0.049 |
| 3 | `Bayes-M8-med` | 63.35 | yes | 37.5 / 25.9 | 19.1 / 20.1 | 0.030 / 0.057 |
| 4 | `Bayes-M17-dist-med` | 62.62 | yes | 36.0 / 26.6 | 14.6 / 14.7 | 0.027 / 0.048 |
| 5 | `Bayes-M8-dist` | 61.94 | yes | 35.4 / 26.6 | 15.6 / 16.2 | 0.028 / 0.050 |
| 6 | `Bayes-M10-dist` | 57.56 | yes | 31.9 / 25.7 | 15.3 / 15.0 | 0.027 / 0.048 |
| 7 | `Bayes-M10-geo` | 53.69 | yes | 30.8 / 22.9 | 15.8 / 16.2 | 0.028 / 0.049 |
| 8 | `Bayes-M17-geo` | 53.49 | yes | 30.8 / 22.7 | 17.1 / 16.8 | 0.029 / 0.050 |
| 9 | `Bayes-M8-geo` | 50.33 | yes | 29.5 / 20.8 | 19.9 / 20.5 | 0.030 / 0.057 |
| 10 | `Bayes-M8-dist-geo` | 48.48 | yes | 27.2 / 21.3 | 17.1 / 17.4 | 0.029 / 0.052 |
| 11 | `Bayes-M17-dist-geo` | 47.07 | yes | 25.6 / 21.4 | 15.6 / 15.6 | 0.029 / 0.050 |
| 12 | `Bayes-M10-dist-geo` | 45.58 | yes | 24.7 / 20.9 | 15.0 / 15.4 | 0.029 / 0.051 |
| 13 | `Bayes-M4-dist` | 37.34 | yes | 21.2 / 16.2 | 22.4 / 24.4 | 0.033 / 0.068 |
| 14 | `Bayes-M4-med` | 37.37 | **spiky** | 21.9 / 15.4 | 31.6 / 34.1 | 0.035 / 0.075 |
| 15 | `Bayes-M4-dist-geo` | 30.31 | **spiky** | 16.3 / 14.1 | 23.7 / 25.9 | 0.034 / 0.067 |
| 16 | `Bayes-M4-geo` | 29.47 | **spiky** | 16.0 / 13.4 | 31.9 / 35.4 | 0.035 / 0.071 |

_(no Best-renewal models in the evaluation)_
