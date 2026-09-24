# Model selection — spatiotemporal

_Generated 2026-09-24T13:09:42+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M10-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_23092026 (processed_at 2026-09-23T12:00:00) |
| Predictions computed on snapshot | LINELIST_23092026 (processed_at 2026-09-23T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-geo` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M10-geo` | 67.15 | yes | 39.7 / 27.5 | 31.5 / 21.9 | 0.021 / 0.037 |
| 2 | `Bayes-M10-med` | 64.36 | yes | 37.6 / 26.7 | 35.2 / 24.5 | 0.022 / 0.038 |
| 3 | `Bayes-M17-geo` | 64.20 | yes | 37.5 / 26.7 | 29.3 / 21.5 | 0.022 / 0.039 |
| 4 | `Bayes-M17-med` | 63.96 | yes | 36.8 / 27.2 | 35.2 / 24.8 | 0.023 / 0.041 |
| 5 | `Bayes-M8-geo` | 61.71 | yes | 36.1 / 25.6 | 30.5 / 24.1 | 0.024 / 0.045 |
| 6 | `Bayes-M8-med` | 61.46 | yes | 36.0 / 25.4 | 37.9 / 28.3 | 0.025 / 0.047 |
| 7 | `Bayes-M8-dist-geo` | 61.04 | yes | 35.9 / 25.1 | 23.9 / 19.7 | 0.023 / 0.041 |
| 8 | `Bayes-M8-dist` | 60.43 | yes | 35.0 / 25.4 | 26.4 / 21.2 | 0.023 / 0.042 |
| 9 | `Bayes-M17-dist-med` | 59.27 | yes | 34.0 / 25.2 | 26.9 / 20.5 | 0.023 / 0.040 |
| 10 | `Bayes-M10-dist` | 54.98 | yes | 30.8 / 24.2 | 27.5 / 20.2 | 0.022 / 0.040 |
| 11 | `Bayes-M10-dist-geo` | 54.65 | yes | 30.2 / 24.4 | 24.3 / 18.5 | 0.022 / 0.039 |
| 12 | `Bayes-M17-dist-geo` | 54.44 | yes | 29.8 / 24.6 | 24.1 / 18.8 | 0.022 / 0.039 |
| 13 | `Bayes-M4-geo` | 43.42 | yes | 25.9 / 17.6 | 40.1 / 35.1 | 0.028 / 0.055 |
| 14 | `Bayes-M4-dist-geo` | 38.57 | yes | 21.4 / 17.2 | 29.8 / 26.6 | 0.027 / 0.052 |
| 15 | `Bayes-M4-dist` | 35.78 | yes | 19.0 / 16.7 | 30.0 / 26.1 | 0.027 / 0.053 |
| 16 | `Bayes-M4-med` | 36.39 | **spiky** | 20.3 / 16.1 | 44.9 / 37.5 | 0.029 / 0.058 |

_(no Best-renewal models in the evaluation)_
