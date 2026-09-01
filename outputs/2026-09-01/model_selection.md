# Model selection — spatiotemporal

_Generated 2026-09-01T09:45:55+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M10-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M10-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_01092026 (processed_at 2026-09-01T12:00:00) |
| Predictions computed on snapshot | LINELIST_01092026 (processed_at 2026-09-01T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M10-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M10-geo` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M10-geo` | 79.40 | yes | 46.7 / 32.7 | 13.5 / 14.3 | 0.021 / 0.037 |
| 2 | `Bayes-M17-geo` | 78.65 | yes | 46.5 / 32.2 | 15.4 / 15.8 | 0.022 / 0.038 |
| 3 | `Bayes-M10-med` | 75.17 | yes | 44.5 / 30.7 | 15.4 / 15.3 | 0.023 / 0.039 |
| 4 | `Bayes-M17-med` | 74.89 | yes | 44.5 / 30.4 | 15.7 / 15.9 | 0.024 / 0.041 |
| 5 | `Bayes-M8-geo` | 73.50 | yes | 44.0 / 29.5 | 20.2 / 20.1 | 0.024 / 0.044 |
| 6 | `Bayes-M8-med` | 71.12 | yes | 42.4 / 28.8 | 20.9 / 20.8 | 0.026 / 0.047 |
| 7 | `Bayes-M8-dist-geo` | 71.12 | yes | 42.0 / 29.1 | 15.4 / 16.4 | 0.023 / 0.040 |
| 8 | `Bayes-M8-dist` | 69.98 | yes | 41.3 / 28.7 | 16.4 / 16.3 | 0.024 / 0.042 |
| 9 | `Bayes-M17-dist-geo` | 65.58 | yes | 38.1 / 27.5 | 12.8 / 13.9 | 0.022 / 0.038 |
| 10 | `Bayes-M17-dist-med` | 65.26 | yes | 37.4 / 27.9 | 13.8 / 14.4 | 0.023 / 0.040 |
| 11 | `Bayes-M10-dist-geo` | 64.44 | yes | 36.7 / 27.7 | 11.5 / 12.9 | 0.022 / 0.038 |
| 12 | `Bayes-M10-dist` | 60.52 | yes | 33.8 / 26.7 | 13.7 / 14.1 | 0.023 / 0.040 |
| 13 | `Bayes-M4-geo` | 46.16 | **spiky** | 26.3 / 19.9 | 34.1 / 35.4 | 0.030 / 0.057 |
| 14 | `Bayes-M4-dist-geo` | 43.20 | **spiky** | 24.1 / 19.1 | 23.4 / 25.5 | 0.028 / 0.054 |
| 15 | `Bayes-M4-med` | 38.59 | **spiky** | 20.3 / 18.3 | 35.5 / 35.9 | 0.032 / 0.061 |
| 16 | `Bayes-M4-dist` | 35.68 | **spiky** | 18.6 / 17.1 | 25.2 / 25.7 | 0.029 / 0.056 |

_(no Best-renewal models in the evaluation)_
