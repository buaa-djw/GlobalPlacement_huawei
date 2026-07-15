#!/usr/bin/env bash
set -euo pipefail
ROOT="${BENCHMARK_ROOT:-testbench/ispd2005}"
DESIGN_DIR="${1:-$ROOT/adaptec2}"
BIN="${GLOBAL_PLACER_BIN:-./build/global_placer}"
TMP=""
cleanup() { [[ -n "$TMP" ]] && rm -rf "$TMP"; }
trap cleanup EXIT
if [[ -d "$DESIGN_DIR" && -f "$DESIGN_DIR/adaptec2.aux" ]]; then
  AUX="$DESIGN_DIR/adaptec2.aux"
elif [[ -d "$DESIGN_DIR" && -f "$DESIGN_DIR/adaptec2.aux.gz" ]]; then
  TMP="$(mktemp -d)"
  for ext in aux nodes nets pl scl wts; do
    gzip -cd "$DESIGN_DIR/adaptec2.$ext.gz" > "$TMP/adaptec2.$ext"
  done
  AUX="$TMP/adaptec2.aux"
elif [[ -f "$DESIGN_DIR" ]]; then
  AUX="$DESIGN_DIR"
else
  echo "adaptec2 benchmark not found under: $DESIGN_DIR" >&2
  exit 1
fi
"$BIN" --aux "$AUX" --bins 64 64 --target-density 0.9
