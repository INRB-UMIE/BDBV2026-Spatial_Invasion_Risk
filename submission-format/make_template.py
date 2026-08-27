#!/usr/bin/env python3
"""Build hubverse-format submission templates and convert model output to them.

Two jobs:

  template  Write a blank submission template for one origin date — every row a
            submitter needs to fill in, with an empty `value` column.
  convert   Convert an INRB `bayes_risk_scores_all_zones.csv` into a valid
            submission file. Doubles as the reference implementation of the
            mapping from a model's own output to this format.

Usage:
  make_template.py template <risk_scores.csv> <origin_date> <out.csv>
  make_template.py convert  <risk_scores.csv> <origin_date> <out.csv>

The risk-scores file is used in BOTH modes, because it defines which zones were
still at risk on that date — the set of rows a submission should contain. Zones
already affected are omitted entirely rather than given a probability of 0.

No third-party dependencies: standard library only, so anyone can run it.
"""

import csv
import datetime as dt
import sys

TARGET = "first confirmed case"
CATEGORIES = ("invasion", "no invasion")
HORIZONS = (1, 2)
COLUMNS = [
    "origin_date",
    "target",
    "horizon",
    "location",
    "target_end_date",
    "output_type",
    "output_type_id",
    "value",
]


def target_end_date(origin_date: str, horizon: int) -> str:
    """Hubverse derived task ID: origin_date + horizon * 7 days.

    Both horizons share a start (the day after origin_date); horizon 2 is
    CUMULATIVE — 'within two weeks', not 'during week two'. Only the end date
    differs, which is what this returns.
    """
    d = dt.date.fromisoformat(origin_date)
    return (d + dt.timedelta(days=7 * horizon)).isoformat()


def at_risk_zones(risk_scores_path: str):
    """Zones with no confirmed case yet, in canonical order.

    `was_active_before` is TRUE for already-affected zones, which cannot be
    newly invaded and so are not a valid prediction target.
    """
    seen, zones = set(), []
    with open(risk_scores_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            zone = row["health_zone"]
            if zone in seen:
                continue
            if row.get("was_active_before", "").strip().upper() == "TRUE":
                continue
            seen.add(zone)
            zones.append(zone)
    return sorted(zones)


def invasion_probabilities(risk_scores_path: str):
    """{(zone, horizon): p} for at-risk zones with a usable probability."""
    out = {}
    with open(risk_scores_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row.get("was_active_before", "").strip().upper() == "TRUE":
                continue
            raw = row.get("p_case_invasion", "").strip()
            if raw in ("", "NA"):
                continue
            p = float(raw)
            if not 0.0 <= p <= 1.0:
                raise ValueError(
                    f"p_case_invasion out of range for "
                    f"{row['health_zone']} h{row['horizon']}: {p}"
                )
            out[(row["health_zone"], int(row["horizon"]))] = p
    return out


def rows(origin_date: str, zones, probs=None):
    """Emit submission rows. probs=None gives a blank template."""
    for horizon in HORIZONS:
        ted = target_end_date(origin_date, horizon)
        for zone in zones:
            if probs is None:
                values = {c: "" for c in CATEGORIES}
            else:
                p = probs.get((zone, horizon))
                if p is None:
                    continue
                values = {"invasion": round(p, 6), "no invasion": round(1.0 - p, 6)}
            for category in CATEGORIES:
                yield {
                    "origin_date": origin_date,
                    "target": TARGET,
                    "horizon": horizon,
                    "location": zone,
                    "target_end_date": ted,
                    "output_type": "pmf",
                    "output_type_id": category,
                    "value": values[category],
                }


def write(path: str, records):
    records = list(records)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(records)
    return len(records)


def main(argv):
    if len(argv) != 5 or argv[1] not in ("template", "convert"):
        print(__doc__, file=sys.stderr)
        return 2
    mode, risk_scores, origin_date, out_path = argv[1:]
    dt.date.fromisoformat(origin_date)  # fail early on a malformed date

    zones = at_risk_zones(risk_scores)
    probs = invasion_probabilities(risk_scores) if mode == "convert" else None
    n = write(out_path, rows(origin_date, zones, probs))
    print(f"{mode}: {n} rows across {len(zones)} at-risk zones -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
