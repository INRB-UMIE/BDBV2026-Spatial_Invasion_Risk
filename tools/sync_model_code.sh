#!/usr/bin/env bash
# Mirror spatiotemporal/ from the private analysis repo into this public repo.
#
# Usage: tools/sync_model_code.sh <ANALYSIS_CHECKOUT> <PUBLIC_REPO_ROOT> [--check]
#   ANALYSIS_CHECKOUT = a checkout of INRB-UMIE/BDBV2026-Analysis at the ref that
#                       produced the published outputs (the pinned ANALYSIS_REF,
#                       NOT necessarily main)
#   PUBLIC_REPO_ROOT  = root of this repository
#   --check           = report drift and exit non-zero if any; change nothing
#
# Design notes
#   * ALLOWLIST ONLY, by extension and location. The public tree is rebuilt from
#     the source each time, so upstream deletions and renames propagate — a stale
#     file cannot survive in the public repo.
#   * Deliberately excluded, and why:
#       outputs/                      run products; published via the outputs filter
#       *.md except METHODS.md        BAYESIAN_INVASION_REPORT, METHODS_AND_RESULTS,
#                                     PLAN_3MONTH_INVASION, PLAN_extensions and
#                                     SPATIAL_INVASION_3MONTH_REPORT are internal
#                                     working documents pending review
#       everything outside spatiotemporal/   Docker, ci/, old_pipelines/,
#                                     spatiotemporal_conditional/ (not run in CI)
#   * METHODS.md is mirrored VERBATIM. Public-facing caveats live in README.md, so
#     a drift check is a plain diff and can never raise a false alarm.

set -euo pipefail

SRC_REPO="${1:?ANALYSIS_CHECKOUT required}"
PUB_REPO="${2:?PUBLIC_REPO_ROOT required}"
MODE="${3:-sync}"

SRC="$SRC_REPO/spatiotemporal"
DST="$PUB_REPO/spatiotemporal"

[[ -d "$SRC" ]] || { echo "error: no spatiotemporal/ under $SRC_REPO" >&2; exit 1; }

# --- Build the allowlisted file list -----------------------------------------
# Top-level .R, tests/*.R, and METHODS.md. Nothing else, at any depth.
list_public_files() {
  local root="$1"
  ( cd "$root" && {
      find . -maxdepth 1 -type f -name '*.R' -printf '%P\n'
      find ./tests -maxdepth 1 -type f -name '*.R' -printf 'tests/%f\n' 2>/dev/null || true
      [[ -f METHODS.md ]] && echo METHODS.md
    } | sort )
}

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

while IFS= read -r rel; do
  [[ -n "$rel" ]] || continue
  if [[ -L "$SRC/$rel" ]]; then
    echo "error: $SRC/$rel is a symlink — refusing to follow it" >&2
    exit 1
  fi
  mkdir -p "$STAGE/$(dirname "$rel")"
  cp "$SRC/$rel" "$STAGE/$rel"
done < <(list_public_files "$SRC")

# --- Guard: the mirror must contain only code and METHODS.md -----------------
while IFS= read -r -d '' f; do
  case "$f" in
    *.R|*/METHODS.md) ;;
    *) echo "error: unexpected file in mirror: ${f#$STAGE/}" >&2; exit 1 ;;
  esac
done < <(find "$STAGE" -type f -print0)

# No data files should ever ride along in a code sync.
if find "$STAGE" -type f -size +2M | grep -q .; then
  echo "error: oversized file in code mirror — data may have leaked in:" >&2
  find "$STAGE" -type f -size +2M -exec ls -lh {} \; >&2
  exit 1
fi

# --- Check mode: report drift, change nothing --------------------------------
if [[ "$MODE" == "--check" ]]; then
  if diff -rq "$STAGE" "$DST" >/dev/null 2>&1; then
    echo "IN SYNC: spatiotemporal/ matches the analysis repo at this ref."
    exit 0
  fi
  echo "DRIFT DETECTED between the analysis repo and the public repo:"
  diff -rq "$STAGE" "$DST" 2>&1 \
    | sed -e "s#$STAGE#analysis#g" -e "s#$DST#public#g"
  exit 1
fi

# --- Sync mode: rebuild the public tree from the staged mirror ---------------
rm -rf "$DST"
mkdir -p "$DST"
cp -R "$STAGE/." "$DST/"

n=$(find "$DST" -type f | wc -l | tr -d ' ')
echo "Synced $n files into spatiotemporal/ ($(du -sh "$DST" | cut -f1))"
