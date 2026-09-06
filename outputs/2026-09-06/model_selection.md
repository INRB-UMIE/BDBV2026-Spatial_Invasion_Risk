# Model selection — spatiotemporal

_Generated 2026-09-06T19:08:24+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_06092026 (processed_at 2026-09-06T12:00:00) |
| Predictions computed on snapshot | LINELIST_06092026 (processed_at 2026-09-06T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-med` | 77.85 | yes | 46.7 / 31.1 | 14.3 / 16.0 | 0.024 / 0.042 |
| 2 | `Bayes-M10-med` | 77.51 | yes | 46.1 / 31.4 | 13.7 / 15.4 | 0.023 / 0.040 |
| 3 | `Bayes-M10-geo` | 76.87 | yes | 46.4 / 30.5 | 12.3 / 14.5 | 0.022 / 0.038 |
| 4 | `Bayes-M17-geo` | 75.57 | yes | 45.5 / 30.1 | 14.2 / 16.3 | 0.023 / 0.039 |
| 5 | `Bayes-M8-dist` | 74.98 | yes | 44.6 / 30.4 | 14.5 / 16.1 | 0.025 / 0.043 |
| 6 | `Bayes-M8-med` | 74.03 | yes | 44.5 / 29.6 | 18.8 / 20.7 | 0.027 / 0.048 |
| 7 | `Bayes-M17-dist-med` | 74.03 | yes | 44.3 / 29.7 | 12.9 / 14.3 | 0.024 / 0.041 |
| 8 | `Bayes-M8-geo` | 71.71 | yes | 44.0 / 27.7 | 17.3 / 20.1 | 0.025 / 0.045 |
| 9 | `Bayes-M8-dist-geo` | 71.61 | yes | 43.0 / 28.6 | 14.2 / 16.3 | 0.024 / 0.041 |
| 10 | `Bayes-M10-dist` | 71.42 | yes | 42.1 / 29.3 | 12.2 / 13.9 | 0.024 / 0.040 |
| 11 | `Bayes-M17-dist-geo` | 68.28 | yes | 40.0 / 28.2 | 12.5 / 14.6 | 0.023 / 0.039 |
| 12 | `Bayes-M10-dist-geo` | 64.72 | yes | 36.9 / 27.9 | 11.3 / 13.2 | 0.023 / 0.039 |
| 13 | `Bayes-M4-geo` | 50.43 | **spiky** | 29.9 / 20.5 | 29.0 / 33.0 | 0.030 / 0.056 |
| 14 | `Bayes-M4-dist-geo` | 48.22 | **spiky** | 28.2 / 20.0 | 21.5 / 24.5 | 0.028 / 0.052 |
| 15 | `Bayes-M4-med` | 45.73 | **spiky** | 26.4 / 19.4 | 32.0 / 34.5 | 0.032 / 0.062 |
| 16 | `Bayes-M4-dist` | 43.77 | **spiky** | 24.8 / 19.0 | 22.3 / 24.3 | 0.029 / 0.056 |

_(no Best-renewal models in the evaluation)_
