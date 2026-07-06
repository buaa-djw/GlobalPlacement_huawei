#!/usr/bin/env bash
set -euo pipefail
LIMBO_ROOT="${LIMBO_ROOT:-/home/djw/Desktop/GlobalPlacement_huawei/thrid_party/Limbo/limbo}"
AUX_PATH="${1:-/home/djw/Desktop/GlobalPlacement_huawei/testbench/ispd2005/adaptec1/adaptec1.aux}"
BUILD_DIR="${BUILD_DIR:-build}"
cmake -S . -B "$BUILD_DIR" -DLIMBO_ROOT="$LIMBO_ROOT"
cmake --build "$BUILD_DIR" -j"$(nproc)"
"$BUILD_DIR/global_placer" --aux "$AUX_PATH"
