# Model selection — spatiotemporal

_Generated 2026-08-27T09:16:02+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_27082026 (processed_at 2026-08-27T12:00:00) |
| Predictions computed on snapshot | LINELIST_27082026 (processed_at 2026-08-27T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-med` | 50.73 | yes | 29.7 / 21.1 | 23.1 / 24.9 | 0.033 / 0.058 |
| 2 | `Bayes-M10-med` | 50.45 | yes | 29.3 / 21.2 | 23.1 / 24.4 | 0.032 / 0.055 |
| 3 | `Bayes-M8-med` | 49.01 | yes | 29.3 / 19.7 | 26.0 / 28.4 | 0.035 / 0.065 |
| 4 | `Bayes-M8-dist` | 48.34 | yes | 28.0 / 20.4 | 23.4 / 25.6 | 0.034 / 0.059 |
| 5 | `Bayes-M17-dist-med` | 47.47 | yes | 26.9 / 20.6 | 22.1 / 24.3 | 0.033 / 0.057 |
| 6 | `Bayes-M10-dist` | 46.30 | yes | 25.7 / 20.6 | 22.5 / 24.4 | 0.033 / 0.056 |
| 7 | `Bayes-M10-geo` | 39.14 | yes | 21.4 / 17.7 | 22.2 / 23.1 | 0.034 / 0.057 |
| 8 | `Bayes-M17-geo` | 38.99 | yes | 21.9 / 17.1 | 23.6 / 24.6 | 0.034 / 0.060 |
| 9 | `Bayes-M17-dist-geo` | 36.93 | yes | 19.5 / 17.5 | 22.5 / 23.8 | 0.035 / 0.059 |
| 10 | `Bayes-M10-dist-geo` | 36.66 | yes | 19.4 / 17.2 | 21.9 / 23.6 | 0.035 / 0.059 |
| 11 | `Bayes-M8-geo` | 36.05 | yes | 19.7 / 16.3 | 26.3 / 27.7 | 0.036 / 0.067 |
| 12 | `Bayes-M8-dist-geo` | 35.96 | yes | 19.6 / 16.3 | 24.7 / 25.8 | 0.036 / 0.062 |
| 13 | `Bayes-M4-dist` | 25.12 | yes | 13.3 / 11.8 | 27.0 / 32.1 | 0.038 / 0.079 |
| 14 | `Bayes-M4-dist-geo` | 23.66 | yes | 12.2 / 11.4 | 27.0 / 32.1 | 0.039 / 0.077 |
| 15 | `Bayes-M4-med` | 24.01 | **spiky** | 12.9 / 11.1 | 35.0 / 40.8 | 0.041 / 0.086 |
| 16 | `Bayes-M4-geo` | 22.46 | **spiky** | 11.5 / 11.0 | 34.4 / 39.0 | 0.041 / 0.082 |

_(no Best-renewal models in the evaluation)_
