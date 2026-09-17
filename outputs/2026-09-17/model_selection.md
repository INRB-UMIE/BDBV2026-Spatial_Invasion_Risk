# Model selection — spatiotemporal

_Generated 2026-09-17T15:31:59+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_17092026 (processed_at 2026-09-17T12:00:00) |
| Predictions computed on snapshot | LINELIST_17092026 (processed_at 2026-09-17T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-med` | 53.29 | yes | 31.0 / 22.3 | 15.1 / 16.4 | 0.029 / 0.051 |
| 2 | `Bayes-M10-med` | 53.14 | yes | 30.4 / 22.7 | 14.2 / 15.2 | 0.028 / 0.047 |
| 3 | `Bayes-M8-med` | 51.40 | yes | 30.6 / 20.8 | 19.4 / 20.4 | 0.032 / 0.058 |
| 4 | `Bayes-M8-dist` | 51.34 | yes | 29.6 / 21.8 | 16.2 / 16.2 | 0.030 / 0.052 |
| 5 | `Bayes-M17-dist-med` | 50.35 | yes | 28.3 / 22.0 | 14.4 / 15.1 | 0.029 / 0.050 |
| 6 | `Bayes-M10-dist` | 49.35 | yes | 26.7 / 22.6 | 13.9 / 14.1 | 0.029 / 0.049 |
| 7 | `Bayes-M10-geo` | 45.60 | yes | 25.0 / 20.6 | 14.3 / 14.0 | 0.029 / 0.049 |
| 8 | `Bayes-M17-geo` | 44.14 | yes | 24.8 / 19.3 | 16.3 / 16.3 | 0.030 / 0.052 |
| 9 | `Bayes-M10-dist-geo` | 43.11 | yes | 22.6 / 20.6 | 14.0 / 13.4 | 0.030 / 0.050 |
| 10 | `Bayes-M8-dist-geo` | 42.56 | yes | 23.5 / 19.1 | 17.4 / 16.1 | 0.031 / 0.054 |
| 11 | `Bayes-M17-dist-geo` | 40.79 | yes | 21.3 / 19.5 | 15.1 / 14.8 | 0.030 / 0.051 |
| 12 | `Bayes-M8-geo` | 40.64 | yes | 22.7 / 18.0 | 19.1 / 19.3 | 0.032 / 0.058 |
| 13 | `Bayes-M4-dist` | 29.38 | yes | 15.7 / 13.6 | 21.3 / 22.3 | 0.033 / 0.067 |
| 14 | `Bayes-M4-dist-geo` | 28.39 | yes | 15.1 / 13.3 | 21.6 / 22.5 | 0.034 / 0.065 |
| 15 | `Bayes-M4-med` | 28.36 | **spiky** | 15.6 / 12.7 | 29.1 / 30.9 | 0.036 / 0.074 |
| 16 | `Bayes-M4-geo` | 26.71 | **spiky** | 14.1 / 12.6 | 28.9 / 30.2 | 0.036 / 0.070 |

_(no Best-renewal models in the evaluation)_
