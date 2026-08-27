#!/usr/bin/env bash
# Strict allowlist filter: sensitive-repo spatiotemporal bundle -> public output set.
#
# Usage: tools/filter_public_outputs.sh <SRC_BUNDLE> <DEST_DIR>
#   SRC_BUNDLE = <sensitive repo>/outputs/<DATE>/spatiotemporal
#   DEST_DIR   = <public repo>/outputs/<DATE>
#
# Design notes
#   * ALLOWLIST ONLY. Never `cp -r` a source directory. A file reaches the public
#     repo only if it is named explicitly below. New pipeline outputs are excluded
#     by default and must be added here deliberately.
#   * Deliberately excluded, and why:
#       bayes_pairwise_import_force.csv      48 MB/day zone x zone matrix
#       bayes_*_national_*.pdf               43 MB/day of vector basemap
#       cascade_*, urban_*, effective_*      3-month projection layer, not the
#                                            operational 1-2 week headline
#       manuscript_figures/, Figure[24]_*    manuscript artefacts, not daily products
#       bayes_ensemble_risk_scores_*.csv     secondary to the primary Bayesian scores
#   * Fails loudly if a REQUIRED file is missing, so a degraded run is never
#     published as if it were complete. Optional files are skipped silently —
#     the pipeline gained outputs over time and early dates legitimately lack some.

set -euo pipefail

SRC="${1:?SRC_BUNDLE required}"
DEST="${2:?DEST_DIR required}"

KEY="$SRC/key_outputs"
REPORTS="$SRC/reports"

[[ -d "$KEY" ]] || { echo "error: no key_outputs under $SRC" >&2; exit 1; }

# --- Required: the run is not publishable without these -----------------------
REQUIRED_KEY=(
  bayes_risk_scores_all_zones.csv     # per-zone invasion probability, CI, rank, priority
  run_info.json                       # provenance: snapshot, cutoff, zone/case counts
  run_info.md
)

# --- Optional: published when present ----------------------------------------
OPTIONAL_KEY=(
  harmonised_confirmed_cases.csv      # cumulative confirmed cases per zone
  model_selection.json                # which model won, and the comparison table
  model_selection.md
  # Evaluation diagnostics (small PDFs, both horizons)
  bayes_topk_precision_h1.pdf
  bayes_topk_precision_h2.pdf
  topk_precision_h1.pdf
  topk_precision_h2.pdf
  bayes_discrimination_summary_h1.pdf
  bayes_discrimination_summary_h2.pdf
  bayes_lfo_forecast_vs_outcome_h1.pdf
  bayes_lfo_forecast_vs_outcome_h2.pdf
  bayes_model_performance_figure3_h1.pdf
  bayes_model_performance_figure3_h2.pdf
  bayes_predicted_vs_observed_over_folds_h1.pdf
  bayes_predicted_vs_observed_over_folds_h2.pdf
  bayes_priority_scatter_h1.pdf
  bayes_priority_scatter_h2.pdf
  bayes_invasion_uncertainty_h1.pdf
  bayes_invasion_uncertainty_h2.pdf
)

OPTIONAL_REPORTS=(
  bayes_parameters.csv                # fitted hazard-ratio terms, CIs, rhat
)

mkdir -p "$DEST"

# copy_one <src> <dest> — refuses symlinks, which could otherwise point outside
# the bundle at a sensitive file while carrying an allowlisted name.
copy_one() {
  local src="$1" dst="$2"
  if [[ -L "$src" ]]; then
    echo "error: $src is a symlink — refusing to follow it" >&2
    exit 1
  fi
  cp "$src" "$dst"
}

missing=0
for f in "${REQUIRED_KEY[@]}"; do
  if [[ -f "$KEY/$f" && ! -L "$KEY/$f" ]]; then
    copy_one "$KEY/$f" "$DEST/$f"
  else
    echo "error: required output missing: key_outputs/$f" >&2
    missing=1
  fi
done
[[ "$missing" -eq 0 ]] || { echo "Refusing to publish an incomplete run." >&2; exit 1; }

for f in "${OPTIONAL_KEY[@]}"; do
  [[ -f "$KEY/$f" ]] && copy_one "$KEY/$f" "$DEST/$f"
done

for f in "${OPTIONAL_REPORTS[@]}"; do
  [[ -f "$REPORTS/$f" ]] && copy_one "$REPORTS/$f" "$DEST/$f"
done

# --- Stamp the analysis commit that produced this run ------------------------
# Written by the pipeline into code_version.json alongside the outputs. Without
# it a published forecast cannot be tied to the code that generated it, so this
# is required for any run from 2026-08-28 onward. The 12 backfilled dates before
# that carry an inferred stamp instead (source="inferred").
CODE_SHA=""
if [[ -f "$SRC/code_version.json" ]]; then
  CODE_SHA="$(jq -r '.commit // empty' "$SRC/code_version.json")"
fi
if [[ -z "$CODE_SHA" && -n "${ANALYSIS_REF:-}" ]]; then
  CODE_SHA="$ANALYSIS_REF"   # fallback: the ref CI checked out
fi
if [[ -z "$CODE_SHA" ]]; then
  echo "error: no code_version.json in the bundle and ANALYSIS_REF unset —" >&2
  echo "       refusing to publish a forecast that cannot be tied to its code." >&2
  exit 1
fi

# --- Redact internal pipeline paths from run_info.json -----------------------
# The upstream file records the container-internal location of the line-list
# pointer. Harmless in itself, but there is no reason to publish the private
# pipeline's directory layout. Snapshot folder + timestamp are kept: they are
# the provenance that makes a published forecast interpretable.
if [[ -f "$DEST/run_info.json" ]]; then
  tmp="$DEST/run_info.json.tmp"
  if jq --arg sha "$CODE_SHA" '
        (if .linelist_snapshot.source_pointer then
           .linelist_snapshot.source_pointer = "<redacted: internal path>"
         else . end)
        | .code_version = {
            analysis_repo: "INRB-UMIE/BDBV2026-Analysis",
            commit: $sha,
            source: "recorded"
          }' "$DEST/run_info.json" > "$tmp" 2>/dev/null; then
    mv "$tmp" "$DEST/run_info.json"
  else
    rm -f "$tmp"
    echo "error: could not rewrite run_info.json (is jq installed?)" >&2
    exit 1
  fi
fi

# --- Guard: nothing oversized should ever land in the public repo -------------
# The agreed footprint is well under 1 MB/day. Anything above 5 MB means the
# allowlist has drifted or an upstream file changed shape.
while IFS= read -r -d '' big; do
  echo "error: $(basename "$big") is $(du -h "$big" | cut -f1) — exceeds the 5 MB public-output limit" >&2
  exit 1
done < <(find "$DEST" -type f -size +5M -print0)

n=$(find "$DEST" -type f | wc -l | tr -d ' ')
echo "Published $n files ($(du -sh "$DEST" | cut -f1)) -> $DEST"
