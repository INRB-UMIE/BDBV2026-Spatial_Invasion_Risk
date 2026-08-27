# Submission format for comparable invasion projections

This folder defines a **file format** for spatial invasion projections, so that
forecasts from other models can be compared directly with the ones published in
`outputs/`. It follows the [hubverse](https://hubverse.io) model-output standard used by
the US COVID-19 and FluSight forecast hubs.

**This is a format specification, not a running hub.** Nothing here collects, validates
or scores submissions automatically. It exists so that anyone modelling this outbreak
produces output that lines up with ours — which makes comparison a join rather than a
negotiation, and leaves the door open to standing up a real hub later.

## What you are predicting

For each health zone with **no confirmed BDBV case** as of the projection date, the
probability that it records its **first** confirmed case within 1 week and within
2 weeks.

Three things trip people up, so they are worth stating plainly:

1. **Horizons are cumulative.** Both horizons start the day after `origin_date`.
   Horizon 2 means *within two weeks*, not *during the second week*. For
   `origin_date` 2026-08-27, horizon 1 covers 28 Aug – 3 Sep and horizon 2 covers
   28 Aug – 10 Sep.
2. **Already-affected zones are not a valid target.** A zone that has already recorded a
   confirmed case cannot be newly invaded. **Omit those zones entirely** — do not submit
   a probability of 0 for them. Scoring a zone that could not possibly have been invaded
   inflates apparent skill.
3. **Time is indexed by symptom onset**, not report date or sample date.

## Files

| File | Purpose |
|---|---|
| `tasks.json` | Machine-readable schema in hubverse format: task IDs, output types, target metadata. |
| `locations.csv` | The 519 canonical health-zone names and their provinces. **Use these strings exactly** — zone naming is the single most common source of mismatch. |
| `template.csv` | A blank submission for `origin_date` 2026-08-27: every row you need, `value` left empty. |
| `make_template.py` | Generates a blank template for any date, and converts INRB-format risk scores into a submission. Standard library only. |

For a **complete worked example**, see any date's
[`outputs/<date>/submission_format/`](../outputs/) — the operational INRB model publishes
its own projection in this exact format every day, so there is always a current,
valid file to compare against.

## The format

One row per prediction, long format:

```
origin_date,target,horizon,location,target_end_date,output_type,output_type_id,value
2026-08-27,first confirmed case,1,Aba,2026-09-03,pmf,invasion,0.0301
2026-08-27,first confirmed case,1,Aba,2026-09-03,pmf,no invasion,0.9699
```

| Column | Meaning |
|---|---|
| `origin_date` | Date the projection was made, `YYYY-MM-DD`. Only data available on this date may be used. |
| `target` | Always `first confirmed case`. |
| `horizon` | `1` or `2` — weeks ahead, cumulative (see above). |
| `location` | Health-zone name, exactly as in `locations.csv`. |
| `target_end_date` | Derived: `origin_date + horizon × 7` days. |
| `output_type` | Always `pmf`. |
| `output_type_id` | `invasion` or `no invasion`. |
| `value` | Probability in [0,1]. **The two rows for a zone-horizon must sum to 1.** |

A binary outcome as a two-row `pmf` looks redundant — you are writing `p` and `1−p`. It
is the hubverse convention for binary targets, and it means
[hubEvals](https://hubverse-org.github.io/hubEvals/) and
[scoringutils](https://epiforecasts.io/scoringutils/) will score your file with no
adaptation: Brier score, log score and calibration come for free.

**Naming:** `<origin_date>-<team_abbr>-<model_abbr>.csv`, e.g.
`2026-08-27-INRB-bayes_renewal.csv`. Keep `team_abbr` and `model_abbr` stable across
dates so a time series of your model can be assembled.

## Producing a submission

```bash
# blank template for a given date
python3 make_template.py template ../outputs/2026-08-27/bayes_risk_scores_all_zones.csv \
    2026-08-27 my-template.csv

# convert INRB-format risk scores into a submission
python3 make_template.py convert ../outputs/2026-08-27/bayes_risk_scores_all_zones.csv \
    2026-08-27 2026-08-27-myteam-mymodel.csv
```

The risk-scores file is read in both modes because it defines which zones were still at
risk on that date. If you are working from your own model, the equivalent in R is a
reshape — one row per zone-horizon becomes two:

```r
library(dplyr); library(tidyr)
submission <- my_forecast |>                      # health_zone, horizon, p_invasion
  filter(!already_affected) |>
  transmute(origin_date  = as.Date("2026-08-27"),
            target       = "first confirmed case",
            horizon, location = health_zone,
            target_end_date = origin_date + horizon * 7,
            invasion = p_invasion, `no invasion` = 1 - p_invasion) |>
  pivot_longer(c(invasion, `no invasion`),
               names_to = "output_type_id", values_to = "value") |>
  mutate(output_type = "pmf", .before = output_type_id)
```

Before sharing a file, check: every zone-horizon has exactly two rows summing to 1, all
values are in [0,1], every `location` appears in `locations.csv`, and no already-affected
zone is present.

## Scoring against what actually happened

Ground truth comes from `outputs/<date>/harmonised_confirmed_cases.csv`: a zone was
invaded by a given date if its cumulative confirmed count first became non-zero on or
before that date. Read
[where confirmed case counts come from](../README.md#where-confirmed-case-counts-come-from)
before scoring — those counts reconcile the DHIS2 line list with INSP situation reports,
and the reconciliation matters for what counts as an invasion.

One caveat worth understanding. The reconciliation takes a **running maximum** over sitrep
revisions, so a zone's first-case date can in principle move if the situation reports are
revised. Ground truth is therefore *as-of* the run that produced it, not immutable. For a
fair retrospective comparison, fix one `harmonised_confirmed_cases.csv` vintage, state
which date you used, and score every model against that same file.

## Extending this

The format deliberately covers only what the operational model currently produces.
Natural extensions, in rough order of usefulness:

- **`sample` output type** for joint draws across zones. Invasion events are correlated
  through the mobility network, and a marginal probability per zone cannot express that.
  Sample output is how hubverse represents joint predictive distributions, and it would
  let a model be scored on its spatial dependence structure rather than only its
  per-zone calibration.
- **Longer horizons** to match the 3-month cascade projection layer.
- **A real hub.** `tasks.json` is already in hubverse config form, so the path is:
  add `admin.json` and `model-metadata/`, then
  [`hubValidations`](https://hubverse-org.github.io/hubValidations/) gives PR-based
  submission checking through a GitHub Action.

> **Note:** `tasks.json` was written against hubverse schema
> [v6.0.0](https://github.com/hubverse-org/schemas/blob/main/NEWS.md) but has not been
> machine-validated. Before treating it as authoritative, run
> `hubAdmin::validate_hub_config()` against it.

## Reference

Consortium of Infectious Disease Modeling Hubs. *Coordinating collaborative infectious
disease modeling projects with the hubverse.* Documentation at
[docs.hubverse.io](https://docs.hubverse.io/), source at
[github.com/hubverse-org](https://github.com/hubverse-org).
