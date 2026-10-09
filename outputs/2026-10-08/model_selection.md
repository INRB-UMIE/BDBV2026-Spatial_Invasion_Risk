# Model selection — spatiotemporal

_Generated 2026-10-09T09:43:22+0000_

**Featured Bayesian model (lowest non-spiky CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M14-fill-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_08102026 (processed_at 2026-10-08T12:00:00) |
| Predictions computed on snapshot | LINELIST_08102026 (processed_at 2026-10-08T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 12, margin 0) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 13, margin 7) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 48.71 | yes | 30.3 / 18.4 | 22.4 / 26.3 | 0.025 / 0.052 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 47.27 | yes | 29.9 / 17.4 | 22.3 / 25.9 | 0.025 / 0.050 |
| 3 | `Bayes-M14-fill-geo` | 45.62 | yes | 27.2 / 18.5 | 22.8 / 27.1 | 0.026 / 0.054 |
| 4 | `Bayes-M14-fill-geo-gtlong` | 45.57 | yes | 26.8 / 18.8 | 22.8 / 27.3 | 0.027 / 0.053 |
| 5 | `Bayes-M14-fill-geo-gtshort` | 45.57 | yes | 27.4 / 18.2 | 22.9 / 26.7 | 0.026 / 0.054 |
| 6 | `Bayes-M14-fill-med` | 43.60 | yes | 25.3 / 18.3 | 26.0 / 31.7 | 0.027 / 0.053 |
| 7 | `Bayes-M13c-fill-geo` | 41.49 | yes | 25.9 / 15.6 | 23.6 / 26.9 | 0.029 / 0.061 |
| 8 | `Bayes-M16-fill-geo` | 42.07 | yes | 25.5 / 16.5 | 28.9 / 29.7 | 0.028 / 0.056 |
| 9 | `Bayes-M10-fill-geo` | 44.24 | yes | 28.9 / 15.3 | 26.4 / 31.2 | 0.028 / 0.061 |
| 10 | `Bayes-M17-fill-geo` | 39.65 | yes | 25.5 / 14.1 | 25.3 / 29.5 | 0.030 / 0.065 |
| 11 | `Bayes-M13c-fill-med` | 40.50 | yes | 24.0 / 16.5 | 26.8 / 31.9 | 0.029 / 0.058 |
| 12 | `Bayes-M16-fill-med` | 41.12 | yes | 24.3 / 16.8 | 32.7 / 34.4 | 0.028 / 0.056 |
| 13 | `Bayes-M17-fill-med` | 40.17 | yes | 24.8 / 15.3 | 28.9 / 34.4 | 0.030 / 0.063 |
| 14 | `Bayes-M10-fill-med` | 40.71 | yes | 25.2 / 15.5 | 30.8 / 37.1 | 0.028 / 0.062 |
| 15 | `Bayes-M13-fill-geo` | 38.34 | yes | 25.5 / 12.8 | 30.0 / 34.1 | 0.031 / 0.072 |
| 16 | `Bayes-M13-fill-med` | 38.03 | yes | 25.0 / 13.0 | 32.1 / 38.5 | 0.032 / 0.072 |
| 17 | `Bayes-M4c-geo` | 31.07 | yes | 18.8 / 12.2 | 32.4 / 36.2 | 0.031 / 0.067 |
| 18 | `Bayes-M4c-med` | 30.10 | yes | 18.1 / 12.0 | 33.8 / 38.4 | 0.031 / 0.069 |
| 19 | `Bayes-M8-fill-geo` | 32.85 | yes | 22.3 / 10.5 | 40.9 / 43.6 | 0.034 / 0.082 |
| 20 | `Bayes-M4-geo` | 22.22 | yes | 13.2 / 9.0 | 40.8 / 43.7 | 0.036 / 0.086 |
| 21 | `Bayes-M8-fill-med` | 25.86 | **spiky** | 17.7 / 8.2 | 45.6 / 52.3 | 0.036 / 0.097 |
| 22 | `Bayes-M4-med` | 19.07 | **spiky** | 11.9 / 7.2 | 44.3 / 51.9 | 0.038 / 0.104 |

_(no Best-renewal models in the evaluation)_
