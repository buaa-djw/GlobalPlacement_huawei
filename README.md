# GlobalPlacement Huawei Baseline

This repository is currently reduced to a minimal C++17 placement-analysis baseline:

- Bookshelf parser front end through the existing Limbo adapter, with a compatibility parser fallback.
- `PlacementDB` with continuous zero-based cell/net/pin/row IDs.
- Exact HPWL evaluation from the input `.pl` placement.
- Exact cell-bin overlap density evaluation for the first numerical density stage of the TCAD 2015 nonsmooth global-placement flow.

The current stage intentionally does **not** include initial placement, clustering, multilevel placement, HPWL or density subgradients, Polak--Ribiere, Moreau, global optimization, legalization, or detailed placement. Running the program does not modify any cell coordinates.

## Build

```bash
rm -rf build
cmake -S . -B build \
  -DLIMBO_ROOT=/home/djw/Desktop/GlobalPlacement_huawei/thrid_party/Limbo
cmake --build build -j$(nproc)
```

ZLIB is required. If Limbo headers/libraries are unavailable, CMake keeps building with the compatibility Bookshelf reader.

## Run

```bash
./build/global_placer --aux <benchmark.aux> --bins 64 64 --target-density 0.9
```

The executable parses Bookshelf, prints database statistics, computes input HPWL, prints orientation counts, and reports exact density metrics. It performs no packing, projection, optimization, legalization, or coordinate update.

## Tests

```bash
ctest --test-dir build --output-on-failure
```

A benchmark smoke helper is available:

```bash
BENCHMARK_ROOT=testbench/ispd2005 scripts/run_adaptec2_baseline.sh
```

## Geometry and Bookshelf Semantics

Coordinates use Bookshelf database units. `cell.x` and `cell.y` are the lower-left corner. Pin positions are currently evaluated as N-orientation semantics:

```text
pin_x = cell.x + 0.5 * cell.width  + offset_x
pin_y = cell.y + 0.5 * cell.height + offset_y
```

Non-`N` orientations are counted and warned about, but no orientation transform or coordinate modification is performed in this stage.

## Exact Density Definition

The core is divided into uniform `nx * ny` bins. For each cell and bin, exact rectangle overlap is computed with double precision:

```text
overlap_area(i,b) = max(0, min(cell.ux, bin.ux) - max(cell.lx, bin.lx))
                  * max(0, min(cell.uy, bin.uy) - max(cell.ly, bin.ly))
D_b = sum_i overlap_area(i,b)
capacity_b = target_density * bin_area
overflow_b = max(0, D_b - capacity_b)
quadratic_penalty = sum_b overflow_b^2
OFR = sum_b overflow_b / total_movable_cell_area
```

Fixed and movable cells both contribute to `D_b`. The OFR denominator uses only total movable cell area; when no movable cells exist, OFR is reported as zero. Cells outside the core are clipped for density accounting and are never moved.
