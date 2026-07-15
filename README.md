# GlobalPlacement Huawei Baseline

This repository is currently a minimal C++17 placement-analysis baseline for reproducing the paper-defined metrics and exact density model from **“Nonsmooth Optimization Method for VLSI Global Placement”**:

- Bookshelf parser front end through the existing Limbo adapter, with a compatibility parser fallback.
- `PlacementDB` with continuous zero-based cell/net/pin/row IDs.
- Exact HPWL evaluation from the input `.pl` placement.
- Equations (9)--(12) exact density model based on exact cell-bin overlap lengths and quadratic overflow penalty.
- Equation (18) OFR.
- Section III-B exact density subgradient with paper-defined random subderivatives at nondifferentiable boundaries.

The current stage intentionally does **not** include initial placement, clustering, multilevel placement, HPWL subgradients, objective-function combination, Polak--Ribiere optimization, global optimization, legalization, or detailed placement. Running the program does not modify any cell coordinates.

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

The default result output is restricted to the paper placement metrics:

```text
==================[Placement Result]===================
HPWL <value>
OFR <value>
```

## Tests

```bash
ctest --test-dir build --output-on-failure
```

A benchmark smoke helper is available:

```bash
BENCHMARK_ROOT=testbench/ispd2005 scripts/run_adaptec2_baseline.sh
```

## Coordinates

Bookshelf database coordinates in `PlacementDB` store each cell lower-left corner:

```text
cell.x, cell.y = lower-left corner
```

The paper density formulas use cell centers, so the density code converts before applying equations (9)--(12):

```text
x_i = cell.x + 0.5 * cell.width
y_i = cell.y + 0.5 * cell.height
x_b = 0.5 * (bin.lx + bin.ux)
y_b = 0.5 * (bin.ly + bin.uy)
```

Pin positions are currently evaluated as N-orientation semantics:

```text
pin_x = cell.x + 0.5 * cell.width  + offset_x
pin_y = cell.y + 0.5 * cell.height + offset_y
```

Non-`N` orientations are warned about by the parser path when relevant, but no orientation transform or coordinate modification is performed in this stage.

## Exact Density Definition

For each axis, equations (9) and (10) are implemented as the exact overlap length between a cell and a bin using center coordinates. For each bin:

```text
D_b = sum_i overlap_x(i,b) * overlap_y(i,b)
capacity_b = target_density * bin_width * bin_height
overflow_b = max(0, D_b - capacity_b)
P = sum_b overflow_b^2
```

The OFR metric follows equation (18):

```text
OFR = sum_b overflow_b / sum_i (cell.width_i * cell.height_i)
```

The denominator is the raw area of all cells before clipping and includes movable standard cells, fixed standard cells, `terminal`, `terminal_NI`, and cells completely outside the placement core.

Fixed and non-movable objects participate in `D_b`, overflow, penalty, and the OFR denominator. Their density subgradient entries are always zero; only movable standard cells receive equation (12) penalty gradients.

At nondifferentiable overlap boundaries, Section III-B is followed directly: positive-side boundaries use `theta2 ~ Uniform[pi/2, pi]`, negative-side boundaries use `theta3 ~ Uniform[0, pi/2]`, and the equal-size coincident case uses `theta ~ Uniform[0, pi]`. The caller supplies the `std::mt19937_64` seed so repeated runs with the same seed are reproducible.
