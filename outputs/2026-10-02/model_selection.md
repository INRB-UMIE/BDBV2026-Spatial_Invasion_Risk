# Model selection — spatiotemporal

_Generated 2026-10-02T17:13:07+0000_

**Featured Bayesian model (lowest non-spiky CV composite, recalibrated log-score axis):** `Bayes-M14-fill-geo`
**Best renewal model:** _none cross-validated in this evaluation_
**Headline (best over all families):** `Bayes-M14-fill-geo`

## Provenance

| Key | Value |
|---|---|
| Selection evaluation | in-memory (this run's CV) |
| Evaluation path | in-memory |
| Selected on data snapshot | LINELIST_02102026 (processed_at 2026-10-02T12:00:00) |
| Predictions computed on snapshot | LINELIST_02102026 (processed_at 2026-10-02T12:00:00) |
| Horizons pooled | 1, 2 |
| Selection rule | Among methods covering the most horizons (partial-CV methods excluded), drop 'spiky' models whose worst-horizon mean rank-of-truth exceeds 1.5x the field median; among the survivors the winner is the lowest CV COMPOSITE — within each horizon models are ranked by AUC-PR skill (desc), mean rank-of-truth (asc) and log score (asc), the three ranks are summed within horizon and then summed across horizons — with ties broken by higher total AUC-PR skill, then lower total rank-of-truth, then lower log score. Gates fall back to the ungated set if they would leave nothing. Log-score axis: RECALIBRATED (16b prequential delta; INVASION_SELECT_ON_RECAL=TRUE), so the axis measures refinement rather than calibration-in-the-large. |
| Log-score axis | recalibrated (INVASION_SELECT_ON_RECAL = TRUE) |
| Pick under RAW axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 9, margin 2) |
| Pick under RECALIBRATED axis | `Bayes-M14-fill-geo` (kernel M14-fill, composite 11, margin 5) |
| Axis sensitivity | none — the same model wins on both axes |
| Cross-check (Bayesian) | OK — matches pipeline pick `Bayes-M14-fill-geo` |
| Cross-check (renewal) | not run this session (standalone) |
| Cross-check (headline) | OK — matches pipeline pick `Bayes-M14-fill-geo` |

### Bayesian model leaderboard (non-spiky, by CV composite; LOWER composite = better)

| Rank | Model | AUC-PR skill (Σ) | Non-spiky | AUC-PR skill (h1/2) | Rank-of-truth (h1/2) | Log-score (h1/2) |
|---|---|---|---|---|---|---|
| 1 | `Bayes-M14-fill-geo-tvar1` | 53.24 | yes | 31.2 / 22.0 | 22.1 / 24.8 | 0.024 / 0.042 |
| 2 | `Bayes-M14-fill-geo-tvrw1` | 52.75 | yes | 31.3 / 21.4 | 22.0 / 24.9 | 0.024 / 0.041 |
| 3 | `Bayes-M14-fill-geo-gtshort` | 50.30 | yes | 28.5 / 21.8 | 22.9 / 25.7 | 0.025 / 0.043 |
| 4 | `Bayes-M14-fill-geo` | 49.29 | yes | 28.0 / 21.3 | 22.7 / 26.0 | 0.026 / 0.044 |
| 5 | `Bayes-M14-fill-geo-gtlong` | 48.90 | yes | 27.5 / 21.4 | 22.6 / 25.9 | 0.026 / 0.045 |
| 6 | `Bayes-M16-fill-geo` | 48.45 | yes | 29.1 / 19.4 | 26.2 / 28.6 | 0.027 / 0.046 |
| 7 | `Bayes-M13c-fill-geo` | 45.55 | yes | 26.7 / 18.9 | 23.4 / 26.6 | 0.027 / 0.049 |
| 8 | `Bayes-M10-fill-geo` | 47.17 | yes | 27.7 / 19.5 | 25.5 / 29.5 | 0.026 / 0.048 |
| 9 | `Bayes-M17-fill-geo` | 46.99 | yes | 28.2 / 18.7 | 24.7 / 28.1 | 0.028 / 0.051 |
| 10 | `Bayes-M14-fill-med` | 46.90 | yes | 26.5 / 20.4 | 26.0 / 30.1 | 0.026 / 0.045 |
| 11 | `Bayes-M17-fill-med` | 47.65 | yes | 29.3 / 18.3 | 29.5 / 33.6 | 0.028 / 0.052 |
| 12 | `Bayes-M16-fill-med` | 46.09 | yes | 27.5 / 18.6 | 29.2 / 32.6 | 0.027 / 0.047 |
| 13 | `Bayes-M13c-fill-med` | 43.42 | yes | 25.7 / 17.7 | 27.3 / 30.9 | 0.028 / 0.049 |
| 14 | `Bayes-M13-fill-geo` | 44.29 | yes | 28.1 / 16.2 | 28.8 / 32.6 | 0.030 / 0.056 |
| 15 | `Bayes-M10-fill-med` | 44.57 | yes | 26.2 / 18.4 | 29.7 / 34.7 | 0.027 / 0.050 |
| 16 | `Bayes-M4c-geo` | 36.85 | yes | 22.5 / 14.3 | 30.1 / 34.9 | 0.029 / 0.053 |
| 17 | `Bayes-M13-fill-med` | 42.46 | yes | 26.9 / 15.6 | 33.2 / 37.6 | 0.031 / 0.058 |
| 18 | `Bayes-M4c-med` | 33.79 | yes | 20.5 / 13.3 | 32.1 / 37.0 | 0.029 / 0.055 |
| 19 | `Bayes-M8-fill-geo` | 35.32 | yes | 22.8 / 12.5 | 37.9 / 42.0 | 0.032 / 0.063 |
| 20 | `Bayes-M4-geo` | 26.68 | yes | 16.5 / 10.2 | 36.8 / 41.7 | 0.033 / 0.067 |
| 21 | `Bayes-M8-fill-med` | 30.48 | **spiky** | 20.9 / 9.6 | 44.8 / 50.3 | 0.034 / 0.070 |
| 22 | `Bayes-M4-med` | 22.13 | **spiky** | 14.2 / 7.9 | 44.0 / 49.9 | 0.034 / 0.072 |

_(no Best-renewal models in the evaluation)_
