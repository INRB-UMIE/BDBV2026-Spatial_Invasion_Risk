# Model selection — spatiotemporal

_Generated 2026-09-16T16:22:40+0000_

**Featured Bayesian model (CV: best non-spiky AUC-PR skill):** `Bayes-M17-med`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M17-med`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_16092026 (processed_at 2026-09-16T12:00:00) |
| Predictions computed on snapshot | LINELIST_16092026 (processed_at 2026-09-16T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons, drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; the winner is the highest total AUC-PR skill (pooled across horizons), ties broken by lower total rank-of-truth then lower log-score. Gate falls back to the ungated set if it would leave nothing. |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M17-med` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M17-med` |

### Bayesian model leaderboard (non-spiky, by total AUC-PR skill; higher = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M17-med` | 58.27 | yes | 34.2 / 24.1 | 14.7 / 16.9 | 0.025 / 0.045 |
| 2 | `Bayes-M10-med` | 57.43 | yes | 33.3 / 24.1 | 13.3 / 15.4 | 0.024 / 0.043 |
| 3 | `Bayes-M8-dist` | 55.05 | yes | 32.1 / 23.0 | 14.2 / 16.5 | 0.026 / 0.047 |
| 4 | `Bayes-M17-dist-med` | 54.78 | yes | 31.8 / 23.0 | 12.5 / 14.9 | 0.025 / 0.045 |
| 5 | `Bayes-M8-med` | 54.64 | yes | 32.3 / 22.3 | 18.7 / 21.2 | 0.028 / 0.052 |
| 6 | `Bayes-M10-geo` | 52.60 | yes | 29.5 / 23.1 | 12.4 / 13.9 | 0.024 / 0.042 |
| 7 | `Bayes-M8-dist-geo` | 52.26 | yes | 30.4 / 21.8 | 14.6 / 15.9 | 0.026 / 0.047 |
| 8 | `Bayes-M17-geo` | 51.92 | yes | 29.0 / 22.9 | 15.0 / 16.3 | 0.025 / 0.045 |
| 9 | `Bayes-M10-dist` | 51.16 | yes | 29.1 / 22.1 | 12.2 / 14.1 | 0.025 / 0.044 |
| 10 | `Bayes-M17-dist-geo` | 50.31 | yes | 28.6 / 21.7 | 12.8 / 14.3 | 0.025 / 0.044 |
| 11 | `Bayes-M8-geo` | 49.84 | yes | 28.3 / 21.6 | 18.3 / 19.9 | 0.027 / 0.051 |
| 12 | `Bayes-M10-dist-geo` | 47.51 | yes | 26.2 / 21.3 | 11.5 / 13.0 | 0.025 / 0.044 |
| 13 | `Bayes-M4-dist-geo` | 31.01 | yes | 16.3 / 14.8 | 21.6 / 23.6 | 0.031 / 0.060 |
| 14 | `Bayes-M4-dist` | 30.38 | yes | 16.3 / 14.1 | 21.2 / 23.1 | 0.032 / 0.061 |
| 15 | `Bayes-M4-geo` | 33.28 | **spiky** | 18.1 / 15.2 | 28.8 / 32.2 | 0.033 / 0.063 |
| 16 | `Bayes-M4-med` | 29.91 | **spiky** | 16.4 / 13.5 | 30.0 / 32.6 | 0.034 / 0.067 |

_(no Best-renewal models in the evaluation)_
